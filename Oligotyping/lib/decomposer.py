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
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import generate_MATRIX_files 
from Oligotyping.utils.utils import generate_ENVIRONMENT_file 
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve


class Node:
    def __init__(self, node_id):
        self.node_id          = node_id
        self.entropy          = None
        self.entropy_tpls     = None
        self.parent           = None
        self.children         = []
        self.discriminants    = None
        self.max_entropy      = None
        self.alignment        = None
        self.unique_alignment = None
        self.read_ids         = []
        self.size             = 0
        self.level            = None


class Decomposer:
    def __init__(self, args = None):
        self.alignment = None
        self.min_entropy = 0.2
        self.number_of_discriminants = 3
        self.min_actual_abundance = 10
        self.output_directory = None
        self.project = None
        self.dataset_name_separator = '_'
        self.generate_sets = False
 
        if args:
            self.alignment = args.alignment
            self.min_entropy = args.min_entropy or 0.2
            self.number_of_discriminants = args.number_of_discriminants or 3
            self.min_actual_abundance = args.min_actual_abundance
            self.output_directory = args.output_directory
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.dataset_name_separator = args.dataset_name_separator
        
        self.decomposition_depth = -1

        # there is a difference between 'average read length' and 'alignment length',
        # therefore there are two different variables to keep that information. the first
        # one is the average length of reads when gaps are removed (if there are any):
        # this may vary from organism to organism. 'alignment length', on the other hand,
        # must be the same for each read in the file.
        self.average_read_length = None
        self.alignment_length = None

        # this variable is going to be used together with 'average read length' and
        # dataset size to calculate the maximum expected number of error driven unique
        # reads:
        self.expected_error = 1 / 250.0

        self.run = Run()
        self.progress = Progress()

        self.root = Node('root')
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
        alignment.reset()

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
        prefix = 'm%.2f-A%d' % (self.min_entropy, self.min_actual_abundance)

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
        self.run.info('root alignment', self.alignment)
        self.run.info('total_seq', pretty_print(self.topology['root'].size))
        self.run.info('output_directory', self.output_directory)
        self.run.info('nodes_directory', self.nodes_directory)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))

        # business time.
        self.generate_raw_topology()
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
        total_reads = 0
        while alignment.next():
            total_reads += len(alignment.ids)
            output.store(alignment, split = False)
        output.close()
        alignment.close()
        return total_reads


    def estimate_expected_max_frequency_of_an_erronous_unique_sequence(self, number_of_reads, average_read_length, expected_error = 1/250.0):
        # maximum number of occurence of an error driven unique sequence among N reads.
        # of course this maximum assumes that all reads are coming from one template,
        # substitution probabilities are homogeneous and there are no systemmatical errors,
        # so it is a mere approximation, but for our purpose, it is going to be enough:

        return round((expected_error * (1 / 3.0)) * ((1 - expected_error) ** (average_read_length - 1)) * number_of_reads) 


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
            # self.node_ids_to_analyze list. while those nodes are being decompoesed,
            # new nodes will appear and need to be analyzed next round. following
            # variable will keep track of the new nodes that emerge, and replace
            # self.node_ids_to_analyze for the next cycle of the main loop.
            new_node_ids_to_analyze = []

            for node_id in self.node_ids_to_analyze:
                
                node = self.topology[node_id]
                
                self.progress.update('Number of nodes to analyze: %d, Analyzing node id: %s (#%d, size: %d)'\
                                                         % (len(self.node_ids_to_analyze),
                                                            node_id,
                                                            self.node_ids_to_analyze.index(node_id),
                                                            node.size))

                if node.size <= self.min_actual_abundance:
                    # FIXME: Finalize this node.
                    continue

                # FIXME: This part is extremely inefficient, take care of it.
                node_file_path_prefix = os.path.join(self.nodes_directory, node_id)
                
                node.unique_alignment = node_file_path_prefix + '.unique'
                self.store_unique_alignment(node.alignment, output_path = node.unique_alignment)

                expected_max_frequency_of_an_erronous_unique_sequence = self.estimate_expected_max_frequency_of_an_erronous_unique_sequence(node.size,
                                                                                                                                            self.average_read_length,
                                                                                                                                            self.expected_error)

                node_entropy_output_path = node_file_path_prefix + '.entropy'
                node.entropy = entropy_analysis(node.unique_alignment, verbose = False, uniqued = True, output_file = node_entropy_output_path)
                node.entropy_tpls = [(node.entropy[i], i) for i in range(0, self.alignment_length)]
                #vis_freq_curve(node.unique_alignment, output_file = node_file_path_prefix + '.png') 

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
                #
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
        self.final_nodes = [n for n in sorted(self.topology.keys()) if not self.topology[n].children]

        # fin.

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
        for node_id in self.topology:
            node = self.topology[node_id]
            topology_text_file_obj.write('%s\t%d\t%s\t%d\t%s\n' % (node.node_id,
                                                               node.size,
                                                               node.parent or '',
                                                               node.level,
                                                               ','.join(node.children) or ''))
        topology_text_file_obj.close()
        self.run.info('topology_text', topology_text_file_path)





if __name__ == '__main__':
    parser = parsers.decomposer()
    decomposer = Decomposer(parser.parse_args())

