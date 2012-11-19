# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import numpy

from Oligotyping.lib import fastalib as u
from Oligotyping.lib.entropy import entropy_analysis

from Oligotyping.utils.utils import append_file
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils.utils import unique_and_store_alignment

class Topology:
    def __init__(self, nodes_output_directory = None):
        self.nodes = {}
        self.alive_nodes = None
        self.final_nodes = None

        self.nodes_output_directory = nodes_output_directory
        
        self.outliers = {}
        self.outlier_reasons = []
        self.next_available_node_id = None
        
        self.average_read_length = None
        self.alignment_length    = None


    def get_new_node_id(self):
        if not self.next_available_node_id:
            new_node_id = 1
            self.next_available_node_id = new_node_id + 1
        else:
            new_node_id = self.next_available_node_id
            self.next_available_node_id += 1

        return '%.9d' % new_node_id

    
    def add_new_node(self, node_id, alignment_path, root = False):
        if not self.nodes_output_directory:
            raise ConfigError, "Nodes output directory has to be declared before adding new nodes"

        node = Node(node_id, self.nodes_output_directory)

        node.alignment = alignment_path
        node.pretty_id = self.get_pretty_id(node_id)
        
        alignment = u.SequenceSource(node.alignment, lazy_init = False)
        node.size = alignment.total_seq

        if root:
            # things to initialize if this is the root node
            node.level = 0

            alignment.next()
            self.alignment_length = len(alignment.seq)
            alignment.reset()
            
            # compute and store the average read length
            # FIXME: this is taking forever for large datasets. there must be a smarter way
            # to do this.
            read_lengths = []
            while alignment.next():
                read_lengths.append(len(alignment.seq.replace('-', '')))
            self.average_read_length = int(round(numpy.mean(read_lengths)))

        alignment.close()

        self.nodes[node_id] = node

        return node

    def gen_alignment_path(self, node_id):
        return os.path.join(self.nodes_output_directory, node_id + '.fa')

    def get_pretty_id(self, node_id):
        pretty_id = node_id
        while 1:
            if pretty_id[0] == '0':
                pretty_id = pretty_id[1:]
            else:
                break
        return pretty_id


    def get_node(self, node_id):
        return self.nodes[node_id]


    def print_node(self, node_id):
        node = self.nodes[node_id]
        print
        print 'Node "%s"' % node
        print '---------------------------------'
        print 'Alive     : %s' % (not node.killed)
        print 'Dirty     : %s' % node.dirty
        print 'Size      : %d' % node.size
        print 'Parent    : %s' % node.parent
        print 'Children  : ', node.children
        print 'Alignment : %s' % node.alignment
        print


    def get_final_count(self):
        return sum([self.nodes[node_id].size for node_id in self.final_nodes])
            

    def get_siblings(self, node_id):
        node = self.nodes[node_id]
        siblings = []
        
        parent_nodes_to_analyze = [self.get_node(node.parent)]
        
        while len(parent_nodes_to_analyze):
            parent = parent_nodes_to_analyze.pop(0)
            
            for candidate in [self.get_node(c) for c in parent.children if c != node.node_id]:
                if candidate.children:
                    parent_nodes_to_analyze.append(candidate)
                else:
                    siblings.append((candidate.size, candidate.node_id))

        # sort siblings least abundant to most
        siblings.sort()
        
        return [s[1] for s in siblings]


    def remove_node(self, node_id, store_content_in_outliers_dict = False, reason = None):
        node = self.nodes[node_id]
        parent = self.nodes[node.parent]
        
        parent.children.remove(node_id)

        if store_content_in_outliers_dict:
            alignment = u.SequenceSource(node.alignment)

            while alignment.next():
                self.store_outlier(alignment.id, alignment.seq, reason)
            alignment.close()

        # get rid of node files.
        self.remove_node_files(node_id)
 
        # it is always sad to pop things
        self.nodes.pop(node_id)

        # and this recursion right here scares the shit out of me:
        if parent.node_id != 'root' and not parent.children:
            self.remove_node(parent.node_id)


    def store_outlier(self, _id, seq, reason = 'unknown_reason'):
        if reason not in self.outlier_reasons:
            self.outlier_reasons.append(reason)
            self.outliers[reason] = []
            
        self.outliers[reason].append((_id, seq),)


    def remove_node_files(self, node_id):
        node = self.nodes[node_id]

        try:
            os.remove(node.alignment)
            os.remove(node.entropy_file)
            os.remove(node.unique_alignment)
        except:
            pass


    def absorb_sibling(self, absorber_node_id, absorbed_node_id):
        # absorbed node gets merged into the absorber node
        absorber = self.get_node(absorber_node_id)
        absorbed = self.get_node(absorbed_node_id)
        
        absorber_parent = self.get_node(absorber.parent)
        absorbed_parent = self.get_node(absorbed.parent)
        
        # append absorbed stuff to the absorber node:
        absorber.read_ids += absorbed.read_ids
        append_file(absorber.alignment, absorbed.alignment)
        absorber.dirty = True
        
        # remove absorbed from the topology
        absorbed_parent.children.remove(absorbed.node_id)
        
        if absorber_parent.node_id != absorbed_parent.node_id and absorbed_parent.node_id != 'root':
            absorbed_parent.size -= absorbed.size
        
        # did the absorbed node's parent just become a final node?
        # if that's the case we gotta remove this from the topology, because reads in this intermediate
        # node was previously split between its child nodes. if all children were absorbed by other nodes
        # this node has no place in the topology anymore.
        if not absorbed_parent.children:
            self.remove_node(absorbed_parent.node_id)

        self.final_nodes.remove(absorbed.node_id)
        self.alive_nodes.remove(absorbed.node_id)
        
        # kill the absorbed
        absorbed.killed = True
        self.remove_node_files(absorbed.node_id)
        
        # refresh the absorbing node
        absorber.refresh()

        
    def get_parent_node(self, node_id):
        return self.nodes[self.nodes[node_id].parent]


    def update_final_nodes(self):
        self.alive_nodes = [n for n in sorted(self.nodes.keys()) if not self.nodes[n].killed]

        # get final nodes sorted by abundance        
        final_nodes_tpls = [(self.nodes[n].size, n) for n in self.alive_nodes if not self.nodes[n].children]
        final_nodes_tpls.sort(reverse = True)
        self.final_nodes = [n[1] for n in final_nodes_tpls]


    def recompute_nodes(self):
        for node_id in self.final_nodes:
            node = self.get_node(node_id)
            if node.dirty:
                node.refresh()


class Node:
    def __init__(self, node_id, output_directory):
        self.node_id            = node_id
        self.pretty_id          = None
        self.representative_seq = None
        self.killed             = False
        self.dirty              = False
        self.entropy            = None
        self.entropy_file       = None
        self.entropy_tpls       = None
        self.parent             = None
        self.children           = []
        self.discriminants      = None
        self.max_entropy        = None
        self.average_entropy    = None
        self.alignment          = None
        self.unique_alignment   = None
        self.read_ids           = []
        self.unique_read_counts = []
        self.size               = 0
        self.level              = None
        self.density            = None
        self.output_directory   = output_directory
        self.file_path_prefix   = os.path.join(self.output_directory, node_id)
        self.freq_curve_img_path = None
        self.competing_unique_sequences_ratio = None
        
    def __str__(self):
        return self.pretty_id

    def do_unique(self):
        self.unique_alignment = self.file_path_prefix + '.unique'
        self.read_ids, self.unique_read_counts, self.representative_seq = \
                    unique_and_store_alignment(self.alignment, output_path = self.unique_alignment)
        self.size = sum(self.unique_read_counts)

    def do_entropy(self):
        self.entropy_file = self.file_path_prefix + '.entropy'
        self.entropy = entropy_analysis(self.unique_alignment, verbose = False, uniqued = True, output_file = self.entropy_file)
        self.entropy_tpls = [(self.entropy[i], i) for i in range(0, len(self.entropy))]
        self.average_entropy = numpy.mean([e for e in self.entropy if e > 0.05] or [0])

    def do_competing_unique_sequences_ratio_and_density(self):
        if len(self.unique_read_counts) == 1:
            self.competing_unique_sequences_ratio = 0
        else:
            self.competing_unique_sequences_ratio = self.unique_read_counts[1] * 1.0 / self.unique_read_counts[0]
                    
        self.density = self.unique_read_counts[0] * 1.0 / sum(self.unique_read_counts)

    def refresh(self):
        self.do_unique()
        self.do_entropy()
        self.do_competing_unique_sequences_ratio_and_density()
        self.dirty = False
