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
from Oligotyping.lib.entropy import entropy_analysis
from Oligotyping.utils.utils import Run
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import get_date
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import human_readable_number
from Oligotyping.utils.utils import generate_MATRIX_files 
from Oligotyping.utils.utils import generate_ENVIRONMENT_file 
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve


class Node:
    def __init__(self, node_id):
        self.node_id            = node_id
        self.pretty_id          = None
        self.representative_seq = None
        self.killed             = False
        self.entropy            = None
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
        self.freq_curve_img_path = None
        self.competing_unique_sequences_ratio = None


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

        self.root = Node('root')
        self.root.pretty_id = 'root'
        self.root.size = sys.maxint
        self.root.level = 0

        self.topology = {'root': self.root}

        # A recursive method could have solved the puzzle entirely in a dataset,
        # however there are a couple of reasons to not approach this problem with
        # recursion. This list is going to keep all leafs of the topology that needs
        # to be analyzed. Things will be added and removed to the list, while
        # self.topology is being formed, and main loop will check this variable its
        # every cycle.
        self.node_ids_to_analyze = ['root']
        self.next_node_id = 1

        self.datasets_dict = {}
        self.datasets = []
        self.alive_nodes = None
        self.final_nodes = None


    def sanity_check(self):
        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise ConfigError, "Alignment file is not accessible: '%s'" % self.alignment

        self.topology['root'].alignment = self.alignment
        alignment = u.SequenceSource(self.alignment, lazy_init = False)
        self.topology['root'].size = alignment.total_seq

        alignment.next()
        self.alignment_length = len(alignment.seq)
        alignment.reset()
        
        # compute and store the average read length
        read_lengths = []
        while alignment.next():
            read_lengths.append(len(alignment.seq.replace('-', '')))
        self.average_read_length = numpy.mean(read_lengths)
        alignment.close()

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
 
        if self.topology['root'].size < self.min_actual_abundance:
            raise ConfigError, "The number of reads in alignment file (%d) is smaller than --min-actual-abundance (%d)" % \
                                                                (self.topology['root'].size, self.min_actual_abundance)
            
    def get_prefix(self):
        prefix = 'm%.2f-A%d-d%d' % (self.min_entropy,
                                    self.min_actual_abundance,
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

        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(self.info_file_path)

        self.run.info('project', self.project)
        self.run.info('run_date', get_date())
        self.run.info('version', __version__)
        self.run.info('root_alignment', self.alignment)
        self.run.info('total_seq', pretty_print(self.topology['root'].size))
        self.run.info('A', self.min_actual_abundance)
        self.run.info('M', self.min_substantive_abundance)
        self.run.info('output_directory', self.output_directory)
        self.run.info('nodes_directory', self.nodes_directory)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))

        # business time.
        self.generate_raw_topology()
        
        if self.generate_frequency_curves:
            self._generate_frequency_curves()

        self.store_topology_dict()
        self.store_topology_text()
        self._generate_datasets_dict()
        generate_ENVIRONMENT_file(self)
        generate_MATRIX_files(self.final_nodes, self)

        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        self.run.store_info_dict(info_dict_file_path)
        self.run.quit()


    def store_unique_alignment(self, alignment_path, output_path):
        output = u.FastaOutput(output_path)
        alignment = u.SequenceSource(alignment_path, unique = True)
        unique_read_counts = []
        while alignment.next():
            unique_read_counts.append(len(alignment.ids))
            output.store(alignment, split = False)
        output.close()
        alignment.close()
        return unique_read_counts


    def _generate_datasets_dict(self):
        self.progress.new('Computing Samples Dict')
        for node_id in self.final_nodes:
            node = self.topology[node_id]
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
        # generate environment file
        self.progress.new('ENVIRONMENT File')
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        f = open(environment_file_path, 'w')
        self.progress.update('Being generated')
        for dataset in self.datasets:
            for node in self.datasets_dict[dataset]:
                f.write("%s\t%s\t%d\n" % (node, dataset, self.datasets_dict[dataset][node]))
        f.close()
        self.progress.end()
        self.run.info('environment_file_path', environment_file_path)


    def dataset_name_from_defline(self, defline):
        return self.dataset_name_separator.join(defline.split('|')[0].split(self.dataset_name_separator)[0:-1])


    def generate_raw_topology(self):
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
  
                node = self.topology[node_id]
                
                self.progress.update('LEVEL: %d: Number of nodes to analyze: %d, Analyzing node id: %s (#%d, size: %d)'\
                                                         % (self.decomposition_depth,
                                                            len(self.node_ids_to_analyze),
                                                            node_id,
                                                            self.node_ids_to_analyze.index(node_id),
                                                            node.size))



                if node.size <= self.min_actual_abundance:
                    # FIXME: Finalize this node.
                    continue

                # FIXME: This part is extremely inefficient, take care of it. Maybe it shouldn't be working with
                #        files but everything should be stored in memory..
                node_file_path_prefix = os.path.join(self.nodes_directory, node_id)
                node.unique_alignment = node_file_path_prefix + '.unique'
                node.unique_read_counts = self.store_unique_alignment(node.alignment, output_path = node.unique_alignment)
                                
                # if the most abundant unique read in a node is smaller than self.min_actual_abundance kill the node.
                if node.unique_read_counts[0] < self.min_substantive_abundance:
                    node.killed = True
                    continue

                # assign the most abundant unique read in the node as representative. it will be useful
                # during the refinement step.
                node_unique_alignment = u.SequenceSource(node.unique_alignment)
                node_unique_alignment.next()
                node.representative_seq = node_unique_alignment.seq
                node_unique_alignment.close()
                
                # competing_unique_sequences_ratio refers to the ratio between the most abundant unique
                # read count and the second most abundant unique read count in a node. smaller the number,
                # better the level of decomposition. however it is important to consider that one organism
                # might be overprinting, increasing the ratio over a closely related organism that trapped
                # in the same node.
                if len(node.unique_read_counts) == 1:
                    node.competing_unique_sequences_ratio = 0
                else:
                    node.competing_unique_sequences_ratio = node.unique_read_counts[1] * 1.0 / node.unique_read_counts[0]

                # 'node density' refers to the ratio of most abundant unique read count to all reads
                # that are accumulated in the node. higher the number, lower the variation within the
                # node.
                node.density = node.unique_read_counts[0] * 1.0 / sum(node.unique_read_counts)

                self.progress.append(' CUSR: %.2f / D: %.2f' % (node.competing_unique_sequences_ratio, node.density))

                if node.competing_unique_sequences_ratio < 0.025 or node.density > 0.85:
                    # Finalize this node.
                    continue

                # find out about the entropy distribution in the given node:
                node_entropy_output_path = node_file_path_prefix + '.entropy'
                node.entropy = entropy_analysis(node.unique_alignment, verbose = False, uniqued = True, output_file = node_entropy_output_path)
                node.entropy_tpls = [(node.entropy[i], i) for i in range(0, self.alignment_length)]
                node.average_entropy = numpy.mean([e for e in node.entropy if e > 0.05])

                self.progress.append(' / ME: %.2f / AE: %.2f' % (max(node.entropy), node.average_entropy))
                
                # IF all the unique reads in the node are smaller than the self.min_substantive_abundance,
                # there is no need to further compose this node.
                if sorted(node.unique_read_counts, reverse = True)[1] < self.min_substantive_abundance:
                    # we are done with this node.
                    continue

                if self.debug:
                    print
                    print node.unique_read_counts[0:10], node.competing_unique_sequences_ratio, node.density
                
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
                
                alignment = u.SequenceSource(node.alignment)
                
                # before we go through the parent alignment to find new set of nodes, we need to keep
                # track of new nodes.
                node_oligo_mapping_dict = {}
                
                new_nodes = {}
                new_node_ids = set([])

                # go through the parent alignment
                while alignment.next():
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

                        new_node = Node(new_node_id)
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
                for new_node_id in new_nodes:
                    # finalize new nodes
                    new_nodes[new_node_id]['file_obj'].close()
                    new_nodes[new_node_id]['node_obj'].size = len(new_nodes[new_node_id]['node_obj'].read_ids)
                    
                    # and add them into the topology:
                    self.topology[new_node_id] = copy.deepcopy(new_nodes[new_node_id]['node_obj'])
                    new_node_ids_to_analyze.append(new_node_id)

            # this is time to set new nodes for the analysis.
            self.node_ids_to_analyze = [n for n in new_node_ids_to_analyze]

        
        #finally:
        self.alive_nodes = [n for n in sorted(self.topology.keys()) if not self.topology[n].killed]
        self.final_nodes = [n for n in self.alive_nodes if not self.topology[n].children]
        self.progress.end()

        if self.debug:
            for node_id in self.final_nodes:
                vis_freq_curve(self.node.unique_alignment, output_file = self.node.unique_alignment + '.png') 

        num_sequences_after_qc = sum([self.topology[node_id].size for node_id in self.final_nodes])
        self.run.info('num_sequences_after_qc', pretty_print(num_sequences_after_qc))
        self.run.info('num_final_nodes', pretty_print(len(self.final_nodes)))

        # fin.


    def _generate_frequency_curves(self):
        self.progress.new('Generating mini entropy figures')
        for i in range(0, len(self.alive_nodes)):
            node = self.topology[self.alive_nodes[i]]
            
            if node.killed:
                continue
            
            self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.alive_nodes)))
                                 
            node.freq_curve_img_path = node.unique_alignment + '.png'
            vis_freq_curve(node.unique_alignment, output_file = node.freq_curve_img_path, mini = True,\
                           title = '%s\n(%s)' % (node.pretty_id, human_readable_number(node.size))) 
        
        self.progress.end()


    def get_new_node_id(self):
        new_node_id = self.next_node_id
        self.next_node_id += 1
        return '%.12d' % new_node_id


    def store_topology_dict(self):
        topology_dict_file_path = self.generate_output_destination('TOPOLOGY.cPickle')
        cPickle.dump(self.topology, open(topology_dict_file_path, 'w'))
        self.run.info('topology_dict', topology_dict_file_path)


    def store_topology_text(self):
        topology_text_file_path = self.generate_output_destination('TOPOLOGY.txt')
        topology_text_file_obj = open(topology_text_file_path, 'w')
        for node_id in self.alive_nodes:
            node = self.topology[node_id]
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

