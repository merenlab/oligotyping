#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

__version__ = '0.1'

import os
import sys
import copy
import numpy
import shutil
import cPickle

from Oligotyping.lib import fastalib as u
from Oligotyping.lib.entropy import quick_entropy
from Oligotyping.lib.entropy import entropy_analysis
from Oligotyping.utils.utils import Run
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import get_date
from Oligotyping.utils.utils import append_file
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import same_but_gaps
from Oligotyping.utils.utils import human_readable_number
from Oligotyping.utils.utils import generate_MATRIX_files 
from Oligotyping.utils.utils import homopolymer_indel_exists
from Oligotyping.utils.utils import generate_ENVIRONMENT_file
from Oligotyping.utils.utils import unique_and_store_alignment
from Oligotyping.utils.utils import get_unit_counts_and_percents
from Oligotyping.utils.utils import get_units_across_datasets_dicts
from Oligotyping.utils.cosine_similarity import cosine_distance
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve

class Topology:
    def __init__(self):
        self.nodes = {}
        self.alive_nodes = None
        self.final_nodes = None
        
        self.outliers = {}
        self.outlier_reasons = []

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


class Decomposer:
    def __init__(self, args = None):
        self.alignment = None
        self.min_entropy = 0.2
        self.number_of_discriminants = 3
        self.min_actual_abundance = 0
        self.min_substantive_abundance = 4
        self.output_directory = None
        self.project = None
        self.dataset_name_separator = '_'
        self.generate_sets = False
        self.generate_frequency_curves = False
        self.debug = False
        self.skip_removing_outliers = False
        self.skip_agglomerating_nodes = False
        self.relocate_outliers = False
        self.maximum_variation_allowed = sys.maxint
        self.store_full_topology = False
        self.merge_homopolymer_splits = False
         
        if args:
            self.alignment = args.alignment
            self.min_entropy = args.min_entropy or 0.2
            self.number_of_discriminants = args.number_of_discriminants or 3
            self.min_actual_abundance = args.min_actual_abundance
            self.min_substantive_abundance = args.min_substantive_abundance
            self.output_directory = args.output_directory
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.dataset_name_separator = args.dataset_name_separator
            self.generate_frequency_curves = args.generate_frequency_curves
            self.skip_removing_outliers = args.skip_removing_outliers
            self.skip_agglomerating_nodes = args.skip_agglomerating_nodes
            self.relocate_outliers = args.relocate_outliers
            self.store_full_topology = args.store_full_topology
            self.merge_homopolymer_splits = args.merge_homopolymer_splits
            self.debug = args.debug
        
        self.decomposition_depth = -1

        # there is a difference between 'average read length' and 'alignment length',
        # therefore there are two different variables to keep that information. the first
        # one is the average length of reads when gaps are removed (if there are any):
        # this may vary from organism to organism. 'alignment length', on the other hand,
        # must be the same for each read in the file.
        self.average_read_length = None
        self.alignment_length = None

        self.run = Run()
        self.progress = Progress()

        self.root = None
        self.topology = Topology()
        
        # A recursive method could have solved the puzzle entirely in a dataset,
        # however there are a couple of reasons to not approach this problem with
        # recursion. This list is going to keep all leafs of the topology that needs
        # to be analyzed. Things will be added and removed to the list, while
        # self.topology is being formed, and main loop will check this variable its
        # every cycle.
        self.node_ids_to_analyze = None
        self.next_node_id = None

        self.datasets_dict = {}
        self.datasets = []
        self.unit_counts = None
        self.unit_percents = None
        self.across_datasets_sum_normalized = {}
        self.across_datasets_max_normalized = {}


    def sanity_check(self):
        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise ConfigError, "Alignment file is not accessible: '%s'" % self.alignment

        # check output associated stuff
        if not self.output_directory:
            self.output_directory = os.path.join(os.getcwd(), '-'.join([self.project.replace(' ', '_'), self.get_prefix()]))
        
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except:
                raise ConfigError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                                          (self.output_directory)
        if not os.access(self.output_directory, os.W_OK):
            raise ConfigError, "You do not have write permission for the output directory: '%s'" % self.output_directory

        self.nodes_directory = self.generate_output_destination('NODES', directory = True)


    def init_topology(self):
        self.root = Node('root', self.nodes_directory)
        self.root.pretty_id = 'root'
        self.root.level = 0       
        
        self.root.alignment = self.alignment
        alignment = u.SequenceSource(self.alignment, lazy_init = False)
        self.root.size = alignment.total_seq

        alignment.next()
        self.alignment_length = len(alignment.seq)
        alignment.reset()
        
        # compute and store the average read length
        read_lengths = []
        while alignment.next():
            read_lengths.append(len(alignment.seq.replace('-', '')))
        self.average_read_length = numpy.mean(read_lengths)
        alignment.close()
 
        if self.root.size < self.min_actual_abundance:
            raise ConfigError, "The number of reads in alignment file (%d) is smaller than --min-actual-abundance (%d)" % \
                                                                (self.root.size, self.min_actual_abundance)

        self.node_ids_to_analyze = ['root']
        self.next_node_id = 1 
        self.topology.nodes['root'] = self.root
            
            
    def get_prefix(self):
        prefix = 'm%.2f-A%d-M%d-d%d' % (self.min_entropy,
                                        self.min_actual_abundance,
                                        self.min_substantive_abundance,
                                        self.number_of_discriminants)

        return prefix


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def decompose(self):
        self.sanity_check()
        self.init_topology()

        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(self.info_file_path)

        self.run.info('project', self.project)
        self.run.info('run_date', get_date())
        self.run.info('version', __version__)
        self.run.info('root_alignment', self.alignment)
        self.run.info('skip_agglomerating_nodes', self.skip_agglomerating_nodes)
        self.run.info('merge_homopolymer_splits', self.merge_homopolymer_splits)
        self.run.info('skip_removing_outliers', self.skip_removing_outliers)
        self.run.info('store_full_topology', self.store_full_topology)
        self.run.info('total_seq', pretty_print(self.topology.nodes['root'].size))
        self.run.info('A', self.min_actual_abundance)
        self.run.info('M', self.min_substantive_abundance)
        self.run.info('output_directory', self.output_directory)
        self.run.info('nodes_directory', self.nodes_directory)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))

        # business time.
        self._generate_raw_topology()
        
        if not self.skip_agglomerating_nodes:
            self._agglomerate_nodes()
            
        if self.merge_homopolymer_splits:
            self._merge_homopolymer_splits()
        
        if not self.skip_removing_outliers:
            self._remove_outliers()
            
        if self.relocate_outliers:
            self._relocate_outliers()
        
        if self.store_full_topology:
            self._store_topology_dict()

        if self.generate_frequency_curves:
            self._generate_frequency_curves()


        # ready for final numbers.
        for reason in self.topology.outlier_reasons:
            self.run.info(reason, pretty_print(len(self.topology.outliers[reason])))

        self.run.info('num_sequences_after_qc', pretty_print(self.topology.get_final_count()))
        self.run.info('num_final_nodes', pretty_print(len(self.topology.final_nodes)))

        self._store_light_topology_dict()
        self._store_topology_text()
        self._generate_datasets_dict()
        self._get_unit_counts_and_percents()
        
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()
 
        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        self.run.store_info_dict(info_dict_file_path)
        self.run.quit()
        

    def _generate_datasets_dict(self):
        self.progress.new('Computing Samples Dict')
        
        for node_id in self.topology.final_nodes:
            node = self.topology.nodes[node_id]
            self.progress.update('Analyzing Node ID: "%s" (size: %d)'\
                                                        % (node_id, node.size))
        
            for read_id in node.read_ids:
                dataset = self.dataset_name_from_defline(read_id)
        
                if not self.datasets_dict.has_key(dataset):
                    self.datasets_dict[dataset] = {}
                    self.datasets.append(dataset)

                if self.datasets_dict[dataset].has_key(node_id):
                    self.datasets_dict[dataset][node_id] += 1
                else:
                    self.datasets_dict[dataset][node_id] = 1

        self.datasets.sort()
        self.progress.end()


    def _generate_ENVIRONMENT_file(self):
        self.progress.new('ENVIRONMENT File')
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        self.progress.update('Being generated')
        
        generate_ENVIRONMENT_file(self.datasets,
                                  self.datasets_dict,
                                  environment_file_path)

        self.progress.end()
        self.run.info('environment_file_path', environment_file_path)        


    def _get_unit_counts_and_percents(self):
        self.progress.new('Unit counts and percents')
        self.progress.update('Data is being generated')
            
        self.unit_counts, self.unit_percents = get_unit_counts_and_percents(self.topology.final_nodes, self.datasets_dict)
            
        self.progress.end()


    def _generate_MATRIX_files(self):
        self.progress.new('Matrix Files')
        self.progress.update('Being generated')
            
        matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")    
            
        generate_MATRIX_files(self.topology.final_nodes,
                              self.datasets,
                              self.unit_counts,
                              self.unit_percents,
                              matrix_count_file_path,
                              matrix_percent_file_path)
            
        self.progress.end()
        self.run.info('matrix_count_file_path', matrix_count_file_path)
        self.run.info('matrix_percent_file_path', matrix_percent_file_path)


    def dataset_name_from_defline(self, defline):
        return self.dataset_name_separator.join(defline.split('|')[0].split(self.dataset_name_separator)[0:-1])


    def _generate_raw_topology(self):
        self.progress.new('Raw Topology')
        # main loop
        while 1:
            self.decomposition_depth += 1
            if not len(self.node_ids_to_analyze):
                self.progress.end()
                break

            # following for loop will go through all nodes that are stored in
            # self.node_ids_to_analyze list. while those nodes are being decomposed,
            # new nodes will appear and need to be analyzed next round. following
            # variable will keep track of the new nodes that emerge, and replace
            # self.node_ids_to_analyze for the next cycle of the main loop.
            new_node_ids_to_analyze = []

            for node_id in self.node_ids_to_analyze:
  
                node = self.topology.nodes[node_id]
                
                
                p = 'LEVEL: %d: Number of nodes to analyze: %d, Analyzing node id: %s (#%d, size: %d)'\
                                                         % (self.decomposition_depth,
                                                            len(self.node_ids_to_analyze),
                                                            node_id,
                                                            self.node_ids_to_analyze.index(node_id),
                                                            node.size)
                self.progress.update(p)

                # FIXME: This part is extremely inefficient, take care of it. Maybe it shouldn't be working with
                #        files but everything should be stored in memory..
                #
                # 1. unique all reads in the node and store them for entropy analysis.
                # 2. store unique read counts.
                # 3. store the most abundant unique read in the node as representative.
                node.do_unique()
                                
                # if the most abundant unique read in a node is smaller than self.min_actual_abundance kill the node
                # and store read information into self.topology.outliers
                if node.unique_read_counts[0] < self.min_substantive_abundance:
                    if node.node_id == 'root':
                        self.progress.end()
                        raise ConfigError, "Number of unique reads in the root node (%d) is less than the declared minimum (%d)." \
                                                % (node.unique_read_counts[0],
                                                   self.min_substantive_abundance)

                    else:
                        # remove the node and store its content.
                        self.topology.remove_node(node.node_id, True, 'min_substantive_abundance_reason')
                        continue

                if node.size < self.min_actual_abundance:
                    # remove the node and store its content.
                    self.topology.remove_node(node.node_id, True, 'min_actual_abundance_reason')                    
                    continue
                
                # competing_unique_sequences_ratio refers to the ratio between the most abundant unique
                # read count and the second most abundant unique read count in a node. smaller the number,
                # better the level of decomposition. however it is important to consider that one organism
                # might be overprinting, increasing the ratio over a closely related organism that trapped
                # in the same node.
                #
                # 'node density' refers to the ratio of most abundant unique read count to all reads
                # that are accumulated in the node. higher the number, lower the variation within the
                # node.
                node.do_competing_unique_sequences_ratio_and_density()

                p += ' CUSR: %.2f / D: %.2f' % (node.competing_unique_sequences_ratio, node.density)
                self.progress.update(p)

                if node.competing_unique_sequences_ratio < 0.025 or node.density > 0.85:
                    # Finalize this node.
                    continue

                # find out about the entropy distribution in the given node:
                node.do_entropy()

                p += ' / ME: %.2f / AE: %.2f' % (max(node.entropy), node.average_entropy)
                self.progress.update(p)
                
                # IF all the unique reads in the node are smaller than the self.min_substantive_abundance,
                # there is no need to further compose this node.
                if sorted(node.unique_read_counts, reverse = True)[1] < self.min_substantive_abundance:
                    # we are done with this node.
                    continue

                # discriminants for this node are being selected from the list of entropy tuples:
                # entropy_tpls look like this:
                #
                #   [(0.01018386415077845, 0), (0.045599806834330125, 1), (0.0093838641507784544, 2), ... ]
                #
                # Probably a function should be called here to make sure discriminants are not high entropy
                # locations driven by homopolymer region associated indels, or dynamicaly set the number of 
                # discriminants for a given node. for instance, if there is one base left in a node that is
                # to define two different organisms, this process should be able to *overwrite* the parameter
                # self.number_of_discriminants.
                node.discriminants = [d[1] for d in sorted(node.entropy_tpls, reverse = True)[0:self.number_of_discriminants] if d[0] > self.min_entropy]

                if not len(node.discriminants):
                    # FIXME: Finalize this node.
                    continue
                
                alignment = u.SequenceSource(node.alignment, lazy_init = False)
                
                # before we go through the parent alignment to find new set of nodes, we need to keep
                # track of new nodes.
                node_oligo_mapping_dict = {}
                
                new_nodes = {}
                new_node_ids = set([])

                # go through the parent alignment
                while alignment.next():
                    if alignment.pos % 1000 == 0 or alignment.pos == alignment.total_seq - 1:
                        self.progress.update(p + ' / %.1f%%' % (alignment.pos * 100.0 / alignment.total_seq))
                    oligo = ''.join([alignment.seq[d] for d in node.discriminants])
                   
                    # new_node_id has to be unique, so things wouldn't overwrite each other.
                    # it was like this at the beginning, and it caused problems:
                    #
                    #   new_node_id = '%s_%s' % ('-'.join(map(lambda d: str(d), node.discriminants)), oligo)
                    if node_oligo_mapping_dict.has_key(oligo):
                        new_node_id = node_oligo_mapping_dict[oligo]
                        new_nodes[new_node_id]['file_obj'].store(alignment, split = False)
                        new_nodes[new_node_id]['node_obj'].read_ids.append(alignment.id)
                    else:
                        new_node_id = self.get_new_node_id()
                        node_oligo_mapping_dict[oligo] = new_node_id

                        new_nodes[new_node_id] = {}
                        new_node_ids = list(new_node_ids)
                        new_node_ids.append(new_node_id)
                        new_node_ids = set(new_node_ids)

                        new_node = Node(new_node_id, self.nodes_directory)
                        pretty_id = new_node_id
                        while 1:
                            if pretty_id[0] == '0':
                                pretty_id = pretty_id[1:]
                            else:
                                break
                        new_node.pretty_id = pretty_id
                        new_node.alignment = os.path.join(self.nodes_directory, new_node_id + '.fa')
                        new_node.read_ids.append(alignment.id)
                        new_node.level = node.level + 1
                       
                        new_node.parent = node.node_id
                        node.children.append(new_node_id)

                        new_nodes[new_node_id]['file_obj'] = u.FastaOutput(new_node.alignment)
                        new_nodes[new_node_id]['file_obj'].store(alignment, split = False)

                        new_nodes[new_node_id]['node_obj'] = new_node

                # all reads in the parent alignment are analyzed.
                new_node_ids = new_nodes.keys()
                for i in range(0, len(new_node_ids)):
                    new_node_id = new_node_ids[i]
                    self.progress.update(p + ' / new nodes are being stored in the topology: %d of %d ' % (i + 1, len(new_node_ids)))
                    # finalize new nodes
                    new_nodes[new_node_id]['file_obj'].close()
                    new_nodes[new_node_id]['node_obj'].size = len(new_nodes[new_node_id]['node_obj'].read_ids)
                    
                    # and add them into the topology:
                    self.topology.nodes[new_node_id] = copy.deepcopy(new_nodes[new_node_id]['node_obj'])
                    new_node_ids_to_analyze.append(new_node_id)

            # this is time to set new nodes for the analysis.
            self.node_ids_to_analyze = [n for n in new_node_ids_to_analyze]

        
        #finally:
        self.topology.update_final_nodes()
        self.progress.end()

        if self.debug:
            for node_id in self.topology.final_nodes:
                vis_freq_curve(self.node.unique_alignment, output_file = self.node.unique_alignment + '.png') 

        self.run.info('num_raw_nodes', pretty_print(len(self.topology.final_nodes)))

        # fin.


    def _refresh_topology(self):
        self.progress.new('Refreshing the topology')
        self.progress.update('Updating final nodes...')

        self.topology.update_final_nodes()

        dirty_nodes = []
        for node_id in self.topology.final_nodes:
            node = self.topology.get_node(node_id)
            if node.dirty:
                dirty_nodes.append(node)
                
        for node in dirty_nodes:
            self.progress.update('Synchronizing dirty node (%d of %d)' % (dirty_nodes.index(node) + 1, len(dirty_nodes)))
            node.refresh()

        self.progress.end()


    def _merge_homopolymer_splits(self):
        self.progress.new('Merge homopolymer splits')
        
        final_nodes = copy.deepcopy(self.topology.final_nodes)
        
        while final_nodes:
            node = self.topology.nodes[final_nodes.pop(0)]
            self.progress.update('Processing node ID: "%s" (remaining: %d)' % (node.pretty_id, len(final_nodes)))

            siblings = self.topology.get_siblings(node.node_id)
            
            while len(siblings):
                sibling = self.topology.nodes[siblings.pop(0)]
                
                if same_but_gaps(node.representative_seq, sibling.representative_seq):
                    if homopolymer_indel_exists(node.representative_seq, sibling.representative_seq):
                        self.topology.absorb_sibling(node.node_id, sibling.node_id)
                    
                        if sibling.node_id in final_nodes:
                            final_nodes.remove(sibling.node_id)
            
        self.progress.end()
        
        # reset temporary stuff
        self.datasets = []
        self.datasets_dict = {}
        self.unit_counts = {}
        self.unit_percents = {}
        
        self._refresh_topology()
 

    def _agglomerate_nodes(self):
        # since we generate nodes based on selected components immediately, some of the nodes will obviously
        # be error driven, since systematic errors can inflate entropy to a point where a column can create
        # its own entropy peak to be selected among all other things. most of the time sequencing errors are
        # random (except some platform-dependent systematic errors, such as homopolymer region indels). so,
        # the errors should not change beta diversity, and frequency distributin patterns of error driven nodes
        # should follow their parent node very tightly (because if a node is very abundant in a dataset, the erronous
        # node stemmed from the same parent will also be relatively abundant in the same dataset). a metric that
        # has the ability to describe similar patterns of frequency distribution, such as cosine similarity, can 
        # be used to agglomerate nodes that diverge from each other by one base and has almost the exact frequency
        # distribution pattern across datasets. this type of process would agglomerate operons, and sometimes very
        # closely related taxa that consistently co-occur, but since minimum entropy decomposition is not interested
        # in diversity much, I am not sure whether this is a bad thing.
        
        self.progress.new('Agglomerating nodes')
        
        # generating a temporary datasets dict.
        self._generate_datasets_dict()
        self._get_unit_counts_and_percents()
        
        self.across_datasets_sum_normalized, self.across_datasets_max_normalized =\
                get_units_across_datasets_dicts(self.topology.final_nodes, self.datasets, self.unit_percents) 

        # compare final nodes with their parent nodes.
        final_nodes = copy.deepcopy(self.topology.final_nodes)
        
        while final_nodes:
            node = self.topology.nodes[final_nodes.pop(0)]
            
            self.progress.update('Processing node ID: "%s" (remaining: %d)' % (node.pretty_id, len(final_nodes)))

            siblings = self.topology.get_siblings(node.node_id)

            while len(siblings):
                sibling = self.topology.nodes[siblings.pop(0)]

                e = quick_entropy([node.representative_seq, sibling.representative_seq])
                
                if len(e) == 1:
                    d = cosine_distance(self.across_datasets_max_normalized[node.node_id], self.across_datasets_max_normalized[sibling.node_id])
                    if d < 0.1:
                        self.topology.absorb_sibling(node.node_id, sibling.node_id)

                        if sibling.node_id in final_nodes:
                            final_nodes.remove(sibling.node_id)
            
        self.progress.end()
        
        # reset temporary stuff
        self.datasets = []
        self.datasets_dict = {}
        self.unit_counts = {}
        self.unit_percents = {}

        self._refresh_topology()


    def _remove_outliers(self):
        # there are potential issues with the raw topology generated. 
        #
        # when one organism dominates a given node (when there are a lot of reads in the node from one template),
        # the entropy peaks indicating variation from organisms with lower abundances may be buried deep, and not
        # meet the minimum entropy for further decomposition criteria. therefore the node would not be fully
        # decomposed. In some cases reads from very distant organisms may end up together in one node.
        #
        # one solution might be (1) going through every read and compare these reads with the representative read of
        # the node, (2) binning reads that are distant from the representative read, and (3) re-assigning them by
        # comparing each read to previously found final nodes.
        
        # to decide at what level should algorithm be concerned about divergent reads in a node, there has to be a
        # threshold that defines what is the maximum variation from the most abundant unique sequence. following
        # variable serves for this purpose. FIXME: this is a very crude way to do it. 
        self.maximum_variation_allowed = int(round(self.average_read_length * 1.0 / 100)) or 1

        self.progress.new('Refined Topology: Removing Outliers')
        for i in range(0, len(self.topology.final_nodes)):
            node = self.topology.nodes[self.topology.final_nodes[i]]
            outlier_seqs = []

            self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.topology.final_nodes)))

            unique_alignment = u.SequenceSource(node.unique_alignment)
            
            # skip the most abundant one
            unique_alignment.next()

            while unique_alignment.next():
                e = quick_entropy([node.representative_seq, unique_alignment.seq])

                if len(e) > self.maximum_variation_allowed:
                    # this read does not belong in this node.
                    outlier_seqs.append(unique_alignment.seq)

            unique_alignment.close()

            if not len(outlier_seqs):
                # no outlier whatsoever. move on to the next.
                continue
            else:
                node.dirty = True
            
            # we have all the ids for this node to be removed. these reads should be remove from the actual alignment.
            alignment = u.SequenceSource(node.alignment)
            alignment_temp = u.FastaOutput(node.alignment + '.temp')
            
            self.progress.append(' / screening node to remove %d outliers' % len(outlier_seqs))
            outlier_seqs = set(outlier_seqs)
            
            while alignment.next():
                if alignment.seq in outlier_seqs:
                    self.topology.store_outlier(alignment.id, alignment.seq, 'maximum_variation_allowed_reason')
                    outlier_seqs.remove(alignment.seq)
                else:
                    alignment_temp.store(alignment, split=False)

            alignment.close()
            alignment_temp.close()

            # now all the reads that we wanted to keep stored in the '.temp' fasta, and all the reads we wanted
            # to remove are stored in outliers dict, we can overwrite original fasta with '.temp'.
            shutil.move(node.alignment + '.temp', node.alignment)
            
        self.progress.end()
        
        self._refresh_topology()
        
    
    def _relocate_outliers(self):
        self.progress.new('Refined Topology: Processing Outliers')
        
        # FIXME: this needs to be re-implemented.
        
        counter = 0
        relocated = 0
        for seq in self.outliers:
            counter += 1
            self.progress.update('Processing %d of %d / relocated: %d' % (counter, len(self.outliers), relocated))
            scores = []
            for node_id in self.topology.final_nodes:
                node = self.topology.nodes[node_id]
                scores.append((len(quick_entropy([seq, node.representative_seq])),
                               node.size, node_id))
            
            best_score = sorted(scores)[0]
            
            entropy_score = best_score[0]
            suggested_node = best_score[2]
            
            if suggested_node != self.outliers[seq]['from'] and entropy_score <= self.maximum_variation_allowed:
                # FIXME: put this guy in suggested_node..
                relocated += 1
                pass
            
        self.progress.end()

        self._refresh_topology()


    def _refresh_final_nodes(self):
        """Refresh nodes using the alignment files"""
        self.progress.new('Refreshing Nodes')
        
        for i in range(0, len(self.topology.final_nodes)):
            node = self.topology.nodes[self.topology.final_nodes[i]]
            self.progress.update('Processing %d of %d (current size: %d)' % (i + 1, len(self.topology.final_nodes), node.size))
            node.refresh()
            
        self.progress.end()
        
        num_sequences_after_qc = sum([self.topology.nodes[node_id].size for node_id in self.topology.final_nodes])
        self.run.info('num_sequences_after_qc', pretty_print(num_sequences_after_qc))
        
    def _generate_frequency_curves(self):
        self.progress.new('Generating mini entropy figures')
        for i in range(0, len(self.topology.alive_nodes)):
            node = self.topology.nodes[self.topology.alive_nodes[i]]
            
            if node.killed:
                continue
            
            self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.topology.alive_nodes)))
                                 
            node.freq_curve_img_path = node.unique_alignment + '.png'
            vis_freq_curve(node.unique_alignment, output_file = node.freq_curve_img_path, mini = True,\
                           title = '%s\n(%s)' % (node.pretty_id, human_readable_number(node.size))) 
        
        self.progress.end()


    def get_new_node_id(self):
        new_node_id = self.next_node_id
        self.next_node_id += 1
        return '%.12d' % new_node_id


    def _store_topology_dict(self):
        self.progress.new('Generating topology dict (full)')
        self.progress.update('Processing %d nodes' % len(self.topology.nodes))
        
        topology_dict_file_path = self.generate_output_destination('TOPOLOGY.cPickle')
        cPickle.dump(self.topology.nodes, open(topology_dict_file_path, 'w'))
        
        self.progress.end()
        
        self.run.info('topology_dict', topology_dict_file_path)


    def _store_light_topology_dict(self):
        self.progress.new('Generating topology dict (lightweight)')
        lightweight_topology_dict = {}
        
        self.progress.update('Processing %d nodes' % len(self.topology.nodes))
        for node_id in self.topology.nodes:
            node = self.topology.nodes[node_id]
            if node.killed:
                continue

            new_node = copy.deepcopy(node)
            new_node.read_ids = None
            new_node.entropy_tpls = None
            
            lightweight_topology_dict[node_id] = new_node

        self.progress.end()
        
        topology_dict_file_path = self.generate_output_destination('TOPOLOGY-LIGHT.cPickle')
        cPickle.dump(lightweight_topology_dict, open(topology_dict_file_path, 'w'))
        self.run.info('topology_light_dict', topology_dict_file_path)


    def _store_topology_text(self):
        topology_text_file_path = self.generate_output_destination('TOPOLOGY.txt')
        topology_text_file_obj = open(topology_text_file_path, 'w')
        for node_id in self.topology.alive_nodes:
            node = self.topology.nodes[node_id]
            topology_text_file_obj.write('%s\t%d\t%s\t%d\t%s\n' \
                                               % (node.node_id,
                                                  node.size,
                                                  node.parent or '',
                                                  node.level,
                                                  ','.join(node.children) or ''))
        topology_text_file_obj.close()
        self.run.info('topology_text', topology_text_file_path)


if __name__ == '__main__':
    pass

