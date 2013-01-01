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

__version__ = '0.1-alpha'

import os
import sys
import copy
import time
import shutil
import cPickle
import logging

from Oligotyping.utils import blast
from Oligotyping.lib import fastalib as u
from Oligotyping.lib.topology import Topology
from Oligotyping.utils.utils import Multiprocessing
from Oligotyping.utils.utils import Run
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import get_date
from Oligotyping.utils.utils import ConfigError
from Oligotyping.utils.utils import run_command 
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import get_pretty_name
from Oligotyping.utils.utils import check_input_alignment
from Oligotyping.utils.utils import human_readable_number
from Oligotyping.utils.utils import generate_MATRIX_files 
from Oligotyping.utils.utils import store_filtered_matrix
from Oligotyping.utils.utils import get_temporary_file_name
from Oligotyping.utils.utils import get_sample_mapping_dict
from Oligotyping.utils.utils import homopolymer_indel_exists
from Oligotyping.utils.utils import mapping_file_simple_check
from Oligotyping.utils.utils import generate_ENVIRONMENT_file
from Oligotyping.utils.utils import get_read_objects_from_file
from Oligotyping.utils.utils import get_unit_counts_and_percents
from Oligotyping.utils.utils import get_units_across_datasets_dicts
from Oligotyping.utils.utils import trim_uninformative_columns_from_alignment
from Oligotyping.utils.utils import get_temporary_file_names_for_BLAST_search
from Oligotyping.utils.utils import get_percent_identity_for_N_base_difference
from Oligotyping.utils.cosine_similarity import cosine_distance
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve


class Decomposer:
    def __init__(self, args = None):
        self.alignment = None
        self.min_entropy = 0.3
        self.number_of_discriminants = 4
        self.min_actual_abundance = 0
        self.min_substantive_abundance = 4
        self.output_directory = None
        self.project = None
        self.dataset_name_separator = '_'
        self.generate_sets = False
        self.generate_frequency_curves = False
        self.skip_refining_topology = False # FIXME: ADD THIS IN PARSERS!
        self.skip_removing_outliers = False
        self.skip_agglomerating_nodes = False
        self.relocate_outliers = False
        self.maximum_variation_allowed = None
        self.store_full_topology = False
        self.merge_homopolymer_splits = False
        self.no_threading = False
        self.number_of_threads = None
        self.log_file_path = None
        self.keep_tmp = False
        self.gen_html = True
        self.skip_figures = False
        self.skip_check_input_file = False
        self.sample_mapping = None
         
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
            self.maximum_variation_allowed = args.maximum_variation_allowed
            self.no_threading = args.no_threading
            self.number_of_threads = args.number_of_threads
            self.keep_tmp = args.keep_tmp
            self.skip_figures = args.skip_figures
            self.skip_check_input_file = args.skip_check_input_file
            self.sample_mapping = args.sample_mapping
            self.gen_html = args.gen_html

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
        self.logger = None

        self.root = None
        self.topology = Topology()
        
        # A recursive method could have solved the puzzle entirely in a dataset,
        # however there are a couple of reasons to not approach this problem with
        # recursion. This list is going to keep all leafs of the topology that needs
        # to be analyzed. Things will be added and removed to the list, while
        # self.topology is being formed, and main loop will check this variable its
        # every cycle.
        self.node_ids_to_analyze = None

        self.datasets_dict = {}
        self.datasets = []
        self.unit_counts = None
        self.unit_percents = None
        self.across_datasets_sum_normalized = {}
        self.across_datasets_max_normalized = {}

        # be smart, turn the threading on if necessary.
        if self.number_of_threads:
            self.no_threading = False

    def check_apps(self):
        try:
            blast.LocalBLAST(None, None, None)
        except blast.ModuleVersionError:
            raise ConfigError, blast.version_error_text
        except blast.ModuleBinaryError:
            raise ConfigError, blast.missing_binary_error_text

        # FIXME: check R modules here.


    def check_dirs(self):
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

        self.tmp_directory = self.generate_output_destination('TMP', directory = True)
        self.nodes_directory = self.generate_output_destination('NODES', directory = True)
        self.figures_directory = self.generate_output_destination('FIGURES', directory = True)
        self.outliers_directory = self.generate_output_destination('OUTLIERS', directory = True)


    def check_input_files(self):
        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise ConfigError, "Alignment file is not accessible: '%s'" % self.alignment

        if self.sample_mapping:
            if (not os.path.exists(self.sample_mapping)) or (not os.access(self.sample_mapping, os.R_OK)):
                raise ConfigError, "Sample mapping file is not accessible: '%s'" % self.sample_mapping

        if self.sample_mapping:
            mapping_file_simple_check(self.sample_mapping)

        if self.skip_check_input_file:
            return 

        self.progress.new('Checking the input FASTA')
        samples = check_input_alignment(self.alignment, self.dataset_name_from_defline, self.progress)
        if not samples:
            raise ConfigError, 'Exiting.'

    def _init_logger(self):
        self.logger = logging.getLogger('decomposer')
        self.log_file_path = self.generate_output_destination('RUNINFO.log')
        
        if os.path.exists(self.log_file_path):
            os.remove(self.log_file_path)
        
        hdlr = logging.FileHandler(self.log_file_path)
        formatter = logging.Formatter('%(asctime)s\t%(levelname)s\t%(message)s')
        hdlr.setFormatter(formatter)
        self.logger.addHandler(hdlr) 
        self.logger.setLevel(logging.DEBUG)


    def _init_topology(self):
        self.progress.new('Initializing topology')
        self.progress.update('May take a while depending on the number of reads...')

        self.topology.nodes_output_directory = self.nodes_directory
        
        reads = get_read_objects_from_file(self.alignment)
        
        self.root = self.topology.add_new_node('root', reads, root = True)
        
        if self.root.size < self.min_actual_abundance:
            raise ConfigError, "The number of reads in alignment file (%d) is smaller than --min-actual-abundance (%d)" % \
                                                                (self.root.size, self.min_actual_abundance)

        self.node_ids_to_analyze = ['root']

        self.progress.end()
            
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
        self.check_apps()
        self.check_dirs()

        # we have just enough to start logging.        
        self._init_logger()
        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(self.info_file_path)

        self.check_input_files()
        
        # we're in business.
        self.run.info('project', self.project)
        self.run.info('run_date', get_date())
        self.run.info('version', __version__)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))
        self.run.info('multi_threaded', not self.no_threading)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('log_file_path', self.log_file_path)
        self.run.info('root_alignment', self.alignment)
        self.run.info('sample_mapping', self.sample_mapping)
        self.run.info('skip_agglomerating_nodes', self.skip_agglomerating_nodes)
        self.run.info('merge_homopolymer_splits', self.merge_homopolymer_splits)
        self.run.info('skip_removing_outliers', self.skip_removing_outliers)
        self.run.info('relocate_outliers', self.relocate_outliers)
        self.run.info('store_full_topology', self.store_full_topology)
        self.run.info('skip_figures', self.skip_figures)
        self.run.info('m', self.min_entropy)
        self.run.info('d', self.number_of_discriminants)
        self.run.info('A', self.min_actual_abundance)
        self.run.info('M', self.min_substantive_abundance)

        # set number of threads to be used
        if not self.number_of_threads:
            self.number_of_threads = Multiprocessing(None).num_thread

        self._init_topology()

        # to decide at what level should algorithm be concerned about divergent reads in a node, there has to be a
        # threshold that defines what is the maximum variation from the most abundant unique sequence. if user did
        # not define this value, it is being set here as follows (FIXME: this is a very crude way to do it): 
        if not self.maximum_variation_allowed:
            self.maximum_variation_allowed = int(round(self.topology.average_read_length * 1.0 / 100)) or 1

        self.run.info('maximum_variation_allowed', self.maximum_variation_allowed)
        self.run.info('total_seq', pretty_print(self.topology.nodes['root'].size))
        self.run.info('average_read_length', self.topology.average_read_length)
        self.run.info('alignment_length', self.topology.alignment_length)
        self.run.info('output_directory', self.output_directory)
        self.run.info('nodes_directory', self.nodes_directory)
        self.run.info('tmp_directory', self.tmp_directory)
        if not self.skip_figures:
            self.run.info('figures_directory', self.nodes_directory)

        # business time.
        self._generate_raw_topology()

        if not self.skip_refining_topology:
            self._refine_topology()
           
        if self.relocate_outliers:
            self._relocate_all_outliers()

        self._generate_datasets_dict()
        self._get_unit_counts_and_percents()

        # all done.        
        self._report_final_numbers()
         
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()

        if self.generate_frequency_curves:
            self._generate_frequency_curves()
 
        self._store_light_topology_dict()
        self._store_topology_text()
        self._store_final_nodes()
        self._store_all_outliers()
        self._store_node_representatives()
        
        if self.store_full_topology:
            self._store_topology_dict()

        for node_id in self.topology.final_nodes:
            self.logger.info('final node: %s (%d)' % (node_id,
                                                      self.topology.nodes[node_id].size))

        self.run.info('end_of_run', get_date())
                
        if (not self.skip_figures):
            self._generate_default_figures()
        
        if (not self.skip_figures) and self.sample_mapping:
            self._generate_exclusive_figures()

        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        self.run.store_info_dict(info_dict_file_path)

        if (not self.keep_tmp):
            shutil.rmtree(self.tmp_directory)

        self.logger.info('fin.')
        self.run.quit()

        # finally:
        if self.gen_html:
            self._generate_html_output()

    def _store_final_nodes(self):
        self.progress.new('Storing final nodes')

        total_final_nodes = len(self.topology.final_nodes)
        
        if self.no_threading:
            for i in range(0, total_final_nodes):
                self.progress.update('%s of %s' % (pretty_print(i + 1),
                                                   pretty_print(total_final_nodes)))
                node_id = self.topology.final_nodes[i]
                node = self.topology.get_node(node_id)
                node.store()

        else:
            
            def worker(data_chunk, shared_counter):
                for node_id in data_chunk:
                    node = self.topology.get_node(node_id)
                    node.store()
                    shared_counter.set(shared_counter.value + 1)

            mp = Multiprocessing(worker, self.number_of_threads)
            data_chunks = mp.get_data_chunks(self.topology.final_nodes, spiral = True)
            shared_counter = mp.get_shared_integer()
            
            for chunk in data_chunks:
                args = (chunk, shared_counter)
                mp.run(args)
        
            while 1:
                num_processes = len([p for p in mp.processes if p.is_alive()])
                
                if not num_processes:
                    # all threads are done.
                    break
        
                self.progress.update('Storing final nodes in %d threads: %d of %d' % (num_processes,
                                                                                      shared_counter.value,
                                                                                      total_final_nodes))
                time.sleep(1)


        self.progress.end()
        

    def _store_outliers(self, reason, output_file_path):
        self.progress.new('Storing outliers')
        output = u.FastaOutput(output_file_path)
            
        self.progress.update('Storing reads removed due to "%s" (size: %d)'\
                                            % (reason, len(self.topology.outliers[reason])))
 
        for read_object in self.topology.outliers[reason]:
            for read_id in read_object.ids:
                output.write_id(read_id)
                output.write_seq(read_object.seq, split = False)
            
        output.close()
        self.progress.end()


    def _store_all_outliers(self):
        for reason in self.topology.outlier_reasons:
            output_file_path = os.path.join(self.outliers_directory, reason + '.fa')
            self._store_outliers(reason, output_file_path)


    def _generate_datasets_dict(self):
        self.progress.new('Computing Samples Dict')
        
        for node_id in self.topology.final_nodes:
            node = self.topology.nodes[node_id]
            self.progress.update('Analyzing Node ID: "%s" (size: %d)'\
                                                        % (node_id, node.size))
        
            for read in node.reads:
                for read_id in read.ids:
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
            
        self.matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        self.matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")    
            
        generate_MATRIX_files(self.topology.final_nodes,
                              self.datasets,
                              self.unit_counts,
                              self.unit_percents,
                              self.matrix_count_file_path,
                              self.matrix_percent_file_path)
            
        self.progress.end()
        self.run.info('matrix_count_file_path', self.matrix_count_file_path)
        self.run.info('matrix_percent_file_path', self.matrix_percent_file_path)


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
                
                
                p = '[LVL %d] Analyzing %d of %d / ID: %s / SIZE: %d'\
                                                         % (self.decomposition_depth,
                                                            self.node_ids_to_analyze.index(node_id) + 1,
                                                            len(self.node_ids_to_analyze),
                                                            node.pretty_id,
                                                            node.size)
                                                         
                self.logger.info('analyzing node id: %s (%d)' % (node_id, node.size))
                self.progress.update(p)

                # if the most abundant unique read in a node is smaller than self.min_actual_abundance kill the node
                # and store read information into self.topology.outliers
                if node.reads[0].frequency < self.min_substantive_abundance:
                    if node.node_id == 'root':
                        self.progress.end()
                        raise ConfigError, "Number of unique reads in the root node (%d) is less than the declared minimum (%d)." \
                                                % (node.reads[0].frequency,
                                                   self.min_substantive_abundance)

                    else:
                        # remove the node and store its content.
                        self.logger.info('remove node (MSA): %s' % node_id)
                        self.topology.remove_node(node.node_id, True, 'min_substantive_abundance_reason')
                        continue

                if node.size < self.min_actual_abundance:
                    # remove the node and store its content.
                    self.topology.remove_node(node.node_id, True, 'min_actual_abundance_reason')                    
                    self.logger.info('remove node (MAA): %s' % node_id)
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

                p += ' / CUSR: %.2f / D: %.2f' % (node.competing_unique_sequences_ratio, node.density)
                self.progress.update(p)

                if node.competing_unique_sequences_ratio < 0.025 or node.density > 0.85:
                    # Finalize this node.
                    self.logger.info('finalize node (CUSR/ND): %s' % node_id)
                    continue

                # find out about the entropy distribution in the given node:
                node.do_entropy()

                p += ' / ME: %.2f / AE: %.2f' % (max(node.entropy), node.average_entropy)
                self.progress.update(p)
                
                # IF the abundance of the second most abundant unique read in the node is smaller than 
                # the self.min_substantive_abundance criteria, there is no need to further decompose
                # this node. because anything spawns from here, will end up in the outlier bin except
                # the most abundant unique read. of course by not decomposing any further we are losing
                # the opportunity to 'purify' this node further, but we are not worried about it,
                # because 'max_allowed_variation' outliers will be removed from this node later on.  
                if node.reads[1].frequency < self.min_substantive_abundance:
                    # we are done with this node.
                    self.logger.info('finalize node (SMA < MSA): %s' % node_id)
                    continue

                # discriminants for this node are being selected from the list of entropy tuples:
                # entropy_tpls look like this:
                #
                #   [(125, 2.0464393446710156), (131, 1.895461844238322), (118, 1.8954618442383218), ... ]
                #
                # Probably a function should be called here to make sure discriminants are not high entropy
                # locations driven by homopolymer region associated indels, or dynamicaly set the number of 
                # discriminants for a given node. for instance, if there is one base left in a node that is
                # to define two different organisms, this process should be able to *overwrite* the parameter
                # self.number_of_discriminants.                    
                node.discriminants = [d[0] for d in node.entropy_tpls[0:self.number_of_discriminants] if d[1] > self.min_entropy]

                if not len(node.discriminants):
                    # FIXME: Finalize this node.
                    self.logger.info('finalize node (ND): %s' % node_id)
                    continue
                else:
                    self.logger.info('using %d D (%s) to decompose: %s'\
                                     % (len(node.discriminants),
                                        ','.join([str(d) for d in node.discriminants]),
                                        node_id))
                
                # before we go through the parent reads to find new set of nodes, we need a variable to keep
                # track of them:
                new_nodes_dict = {}

                # go through the parent reads
                counter = 0
                while 1:
                    counter += 1
                    if counter % 1000 == 0 or counter == node.size - 1:
                        self.progress.update(p + ' / %.1f%%' % (counter * 100.0 / node.size))
                    
                    if not node.reads:
                        # all reads were processed, break while loop.
                        break
                    
                    read = node.reads.pop()

                    oligo = ''.join([read.seq[d] for d in node.discriminants])
                   
                    if new_nodes_dict.has_key(oligo):
                        new_nodes_dict[oligo]['reads'].append(read)
                    else:
                        new_node_id = self.topology.get_new_node_id()
                        new_nodes_dict[oligo] = {}
                        new_nodes_dict[oligo]['node_id'] = new_node_id
                        new_nodes_dict[oligo]['reads'] = [read]
                
                
                # all reads in the parent reads are analyzed. time to add spawned nodes into the topology.
                oligos = new_nodes_dict.keys()
                len_oligos = len(oligos)
                for i in range(0, len_oligos):
                    self.progress.update(p + ' / new nodes %d of %d ' % (i + 1, len_oligos))

                    new_node = self.topology.add_new_node(new_nodes_dict[oligos[i]]['node_id'],
                                                          new_nodes_dict[oligos[i]]['reads'],
                                                          parent_id = node.node_id)

                    new_node_ids_to_analyze.append(new_node.node_id)
                    self.logger.info('new node: %s' % new_node.node_id)

            # this is time to set new nodes for the analysis.
            self.node_ids_to_analyze = [n for n in new_node_ids_to_analyze]

        
        #finally:
        self.progress.end()
        self.topology.update_final_nodes()

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

        if self.no_threading:
            for node in dirty_nodes:
                self.progress.update('Synchronizing dirty nodes (%d of %d)' % (dirty_nodes.index(node) + 1, len(dirty_nodes)))
                node.refresh()
        else:
            # worker function..
            def worker(data_chunk, shared_counter, results_array):
                for node in data_chunk:
                    node.refresh()
                    results_array.append(node)
                    shared_counter.set(shared_counter.value + 1)

            mp = Multiprocessing(worker, self.number_of_threads)
            data_chunks = mp.get_data_chunks(dirty_nodes, spiral = True)
            shared_counter = mp.get_shared_integer()
            results_array = mp.get_empty_shared_array()
            
            for chunk in data_chunks:
                args = (chunk, shared_counter, results_array)
                mp.run(args)
        
            while 1:
                num_processes = len([p for p in mp.processes if p.is_alive()])
                
                if not num_processes:
                    # all threads are done. replace nodes with resulting nodes.
                    for node in results_array:
                        self.topology.nodes[node.node_id] = node
                    break
        
                self.progress.update('Processing in %d threads: %d of %d' % (num_processes,
                                                                             shared_counter.value,
                                                                             len(dirty_nodes)))
                time.sleep(1)

        self.progress.end()



    def _refine_topology(self):
        # FIXME: this is the most sophisticated and second most important part of the algorithm.
        # explain it nicely, meren, kthxbye (OP will surely deliver).

        iteration = 0
        it_is_OK_to_pass_this = lambda: iteration > 0 and (not len(self.topology.zombie_nodes))
        
        while 1:
            self.logger.info('refine topology iteration: %d' % iteration)

            if it_is_OK_to_pass_this():
                self.logger.info('end of refine topology')

                # report the number of outliers removed during the refinement step
                removed_outliers_total = 0
                for reason in self.topology.outlier_reasons:
                    count = sum([read_obj.frequency for read_obj in self.topology.outliers[reason]])
                    removed_outliers_total += count
                    self.run.info('removed_%s' % reason, pretty_print(count))
                self.run.info('removed_outliers_total', pretty_print(removed_outliers_total))

                break

            if not self.skip_agglomerating_nodes:
                if it_is_OK_to_pass_this():
                    pass
                else:
                    self._agglomerate_nodes(iteration)

            if self.merge_homopolymer_splits:
                if it_is_OK_to_pass_this():
                    pass
                else:
                    self._merge_homopolymer_splits(iteration)

            # if iteration != 0: all zombie nodes that were resulted from the iteration-- considered by functions that refines
            # the topology. whatever is left here as zombie, must actually not be a zombie, but a respectable member of the
            # society. so. it is time to reset the topology.zombie.
            self.topology.zombie_nodes = []
            
            if not self.skip_removing_outliers:
                self._remove_outliers(iteration, standby_bin_only = iteration > 0)

                # if iteration == 0, read this first:
                # we removed all the outliers that violate the maximum variation and made sure our final nodes are
                # 'pure'. nice. but now there is another problem to fix: the frequency of some of the outlier
                # read_objects, especially the ones that are coming from very large nodes may actually be larger
                # than the minimum_subtantive_abundance criteria. so, they are the reads that were supposed to be
                # a node, but were betrayed by the nature of the dataset or the sequencing platform. either due to
                # noise, or the identity to the representative of the node they were trapped in based on the
                # discriminant that happened to be chosen to decompose that branch of the topology. now we are going
                # to identify them, and move them out of the outliers bin just to attach them to the topology as a new node.
                # but since they haven't gone through previous steps of refinements, we can't treat them as a regular
                # node. so they will be marked as 'zombie' nodes. we are actually in a 'while True' loop and the only
                # condition to break is to make sure the zombie_nodes bin is empty. in an ideal world it shouldn't take
                # more than two cycles (zombie bins are introduced, loop goes back to the beginning, they are taken care
                # of, and when we are here the second time no more zombie bins are found, we're golden). but for the sake
                # of robustness, I didn't want to rely on this and implement this part of the algorithm as a complete
                # state machine.
                
                abundant_reads_in_outlier_bin = []
                
                if self.topology.outliers.has_key('maximum_variation_allowed_reason'):    
                    abundant_reads_in_outlier_bin = [read_object for read_object in \
                                                        self.topology.outliers['maximum_variation_allowed_reason'] \
                                                            if read_object.frequency > self.min_substantive_abundance]
                
                self.progress.new('Abundant Outliers Bin; ITER %d' % (iteration))
                number_of_abundant_reads_in_outlier_bin = len(abundant_reads_in_outlier_bin)
                for i in range(0, number_of_abundant_reads_in_outlier_bin):
                    self.progress.update('%d of %d' % (i + 1, number_of_abundant_reads_in_outlier_bin))

                    read_object = abundant_reads_in_outlier_bin[i]

                    new_node_id = self.topology.get_new_node_id()
                    self.topology.add_new_node(new_node_id, [read_object], parent_id = 'root')
                    self.topology.zombie_nodes.append(new_node_id)
                    self.topology.outliers['maximum_variation_allowed_reason'].remove(read_object)

                    self.topology.final_nodes.append(new_node_id)
                    self.topology.alive_nodes.append(new_node_id)

                    self.logger.info('new zombie: %s' % new_node_id)
                self.progress.end()
                
            iteration += 1


    def _merge_homopolymer_splits(self, iteration):
        # FIXME: this actually traverses the sibling nodes in the topology to merge homopoymer splits, but
        # it doesn't make any sense? who says homopolymer splits are going to be in the same branch of the
        # topology? they might have been split way before and ended up extremely distant places in the topology
        # and this has to be fixed at some point (which unfortunately will fuck up the performance drastically
        # once it is fixed, but whatever. science > my pride).
       
        nz = pretty_print(len(self.topology.zombie_nodes))
        self.progress.new('Merging HP splits :: ITER %d%s' % (iteration,
                                                              ' #Z: %s' % (nz if nz else '')))

        dealing_with_zombie_nodes = False

        if self.topology.zombie_nodes:
            nodes = copy.deepcopy(self.topology.zombie_nodes)
            dealing_with_zombie_nodes = True
        else:
            nodes = copy.deepcopy(self.topology.final_nodes)

        # ---
        # use blastn to get the similarity dict for gaps.
        query, target, output = get_temporary_file_names_for_BLAST_search(prefix = "HPS_%d_" % iteration,\
                                                                          directory = self.tmp_directory)

        self.topology.store_node_representatives(nodes, query)
        self.topology.store_node_representatives(self.topology.final_nodes, target)

        
        min_percent_identity = get_percent_identity_for_N_base_difference(self.topology.average_read_length, N = 1)
        params = "-perc_identity %.2f" % (min_percent_identity)
        b = self._perform_blast(query, target, output, params, job = 'HPS')
        
        self.progress.update('Generating similarity dict from blastn results')
        similarity_dict = b.get_results_dict(mismatches = 0, gaps = 1)

        node_ids = set(similarity_dict.keys())
        nodes_to_skip = set()

        while node_ids:
            node = self.topology.nodes[node_ids.pop()]

            if node.node_id in nodes_to_skip:
                continue

            self.progress.update('Processing node ID: "%s" (remaining: %d)' % (node.pretty_id, len(node_ids)))
            
            for sibling_id in similarity_dict[node.node_id]:
                if sibling_id in nodes_to_skip:
                    continue
                
                sibling = self.topology.nodes[sibling_id]


                if homopolymer_indel_exists(node.representative_seq, sibling.representative_seq):
                    if dealing_with_zombie_nodes:
                        self.topology.merge_nodes(sibling.node_id, node.node_id)
                        self.topology.standby_bin.append(sibling.node_id)
                        nodes_to_skip.add(node.node_id)
                        self.logger.info('zombie node merged (HPS): %s -> %s' % (node.node_id, sibling.node_id))
                        break
                    else:
                        self.topology.merge_nodes(node.node_id, sibling.node_id)
                        nodes_to_skip.add(sibling.node_id)
                        self.logger.info('nodes merged (HPS): %s -> %s' % (node.node_id, sibling.node_id))

                        # this shouldn't be necessary, but just in case:
                        if sibling.node_id in node_ids:
                            node_ids.remove(sibling.node_id)

        # reset temporary stuff
        self.datasets = []
        self.datasets_dict = {}
        self.unit_counts = {}
        self.unit_percents = {}
        
        self.progress.end()
        self._refresh_topology()


    def _agglomerate_nodes(self, iteration):
        # since we generate nodes based on selected components immediately, some of the nodes will obviously
        # be error driven, since systematic errors can inflate entropy to a point where a column can create
        # its own entropy peak to be selected among all other things. most of the time sequencing errors are
        # random (except some platform-dependent systematic errors, such as homopolymer region indels). so,
        # the errors should not change beta diversity, and frequency distribution patterns of error driven nodes
        # should follow their parent node very tightly (because if a node is very abundant in a dataset, the erroneous
        # node stemmed from the same parent will also be relatively abundant in the same dataset). a metric that
        # has the ability to describe similar patterns of frequency distribution, such as cosine similarity, can 
        # be used to agglomerate nodes that diverge from each other by one base and has almost the exact frequency
        # distribution pattern across datasets. this type of process would agglomerate operons, and sometimes very
        # closely related taxa that consistently co-occur, but since minimum entropy decomposition is not interested
        # in diversity much, I am not sure whether this is a bad thing.


        # generating a temporary datasets dict.
        self._generate_datasets_dict()
        self.logger.info('temp datasets dict is ready; %d samples found' % (len(self.datasets)))
        self._get_unit_counts_and_percents()
        self.logger.info('temp unit counts and percents dicts are ready')
        
        self.across_datasets_sum_normalized, self.across_datasets_max_normalized =\
                get_units_across_datasets_dicts(self.topology.final_nodes, self.datasets, self.unit_percents) 
        
        nz = pretty_print(len(self.topology.zombie_nodes))
        self.progress.new('Agglomerating nodes :: ITER %d%s' % (iteration,
                                                                ' #Z: %s' % (nz if nz else '')))
 
        dealing_with_zombie_nodes = False

        if self.topology.zombie_nodes:
            nodes = copy.deepcopy(self.topology.zombie_nodes)
            dealing_with_zombie_nodes = True
        else:
            nodes = copy.deepcopy(self.topology.final_nodes)

        # ---
        # use blastn to get the similarity dict. I don't like it here, but I'll take care of it later.        
        query, target, output = get_temporary_file_names_for_BLAST_search(prefix = "AN_%d_" % iteration,\
                                                                          directory = self.tmp_directory)
        self.topology.store_node_representatives(nodes, query)
        self.topology.store_node_representatives(self.topology.final_nodes, target)

        min_percent_identity = get_percent_identity_for_N_base_difference(self.topology.average_read_length)
        params = "-perc_identity %.2f" % (min_percent_identity)

        b = self._perform_blast(query, target, output, params, job = 'AN')
        
        self.progress.update('Generating similarity dict from blastn results')
        similarity_dict = b.get_results_dict(mismatches = 1, gaps = 0)

        node_ids = set(similarity_dict.keys())
        nodes_to_skip = set()

        while node_ids:
            node = self.topology.nodes[node_ids.pop()]

            if node.node_id in nodes_to_skip:
                continue

            self.progress.update('Processing node ID: "%s" (remaining: %d)' % (node.pretty_id, len(node_ids)))
            
            for sibling_id in similarity_dict[node.node_id]:
                if sibling_id in nodes_to_skip:
                    continue
                
                sibling = self.topology.nodes[sibling_id]

                d = cosine_distance(self.across_datasets_max_normalized[node.node_id], self.across_datasets_max_normalized[sibling.node_id])

                if d < 0.1:
                    if dealing_with_zombie_nodes:
                        self.topology.merge_nodes(sibling.node_id, node.node_id)
                        self.topology.standby_bin.append(sibling.node_id)
                        nodes_to_skip.add(node.node_id)
                        self.logger.info('zombie node merged (AN, d: %.2f): %s -> %s' % (d, node.node_id, sibling.node_id))
                        break
                    else:
                        self.topology.merge_nodes(node.node_id, sibling.node_id)
                        nodes_to_skip.add(sibling.node_id)
                        self.logger.info('nodes merged (AN, d: %.2f): %s -> %s' % (d, node.node_id, sibling.node_id))

                        # this shouldn't be necessary, but just in case:
                        if sibling.node_id in node_ids:
                            node_ids.remove(sibling.node_id)


        # reset temporary stuff
        self.datasets = []
        self.datasets_dict = {}
        self.unit_counts = {}
        self.unit_percents = {}

        self.progress.end()
        self._refresh_topology()


    def _remove_outliers(self, iteration, standby_bin_only = False):
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


        sb = pretty_print(len(self.topology.standby_bin))
        self.progress.new('Removing Outliers :: ITER %d%s' % (iteration,
                                                              ' #SB: %s' % (sb if sb else '')))

        if standby_bin_only:
            node_list = copy.deepcopy(self.topology.standby_bin)
            self.topology.standby_bin = []
        else:
            node_list = copy.deepcopy(self.topology.final_nodes)


        max_percent_identity = get_percent_identity_for_N_base_difference(self.topology.average_read_length,
                                                                          self.maximum_variation_allowed)

        if self.no_threading:
            # no threading
            for i in range(0, len(node_list)):
                node_id = node_list[i]
                node = self.topology.nodes[node_id]
                
                self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.topology.final_nodes)))

                job = 'XO_%s_' % node_id
                query, target, output = get_temporary_file_names_for_BLAST_search(prefix = job,\
                                                                                  directory = self.tmp_directory)
                id_to_read_object_dict = {}
                for read_obj in node.reads[1:]:
                    id_to_read_object_dict[read_obj.md5id] = read_obj

                query_obj = u.FastaOutput(query)
                for _id in id_to_read_object_dict:
                    query_obj.write_id(_id)
                    query_obj.write_seq(id_to_read_object_dict[_id].seq.replace('-', ''), split = False)
                query_obj.close()

                target_obj = u.FastaOutput(target)
                target_obj.write_id(node.reads[0].md5id)
                target_obj.write_seq(node.reads[0].seq.replace('-', ''), split = False)            
                target_obj.close()

                b = self._perform_blast(query, target, output, params = '', job = job)
                    
                similarity_dict = b.get_results_dict(max_identity = max_percent_identity)
                               
                if len(similarity_dict):
                    node.dirty = True
                else:
                    continue
    
                self.progress.append(' / screening node to remove %d outliers' % len(similarity_dict))

                for _id in similarity_dict:
                    outlier_read_object = id_to_read_object_dict[_id]
                    node.reads.remove(outlier_read_object)
                    self.topology.store_outlier(outlier_read_object, 'maximum_variation_allowed_reason')
 
                self.logger.info('%d outliers removed from node: %s'\
                            % (sum([id_to_read_object_dict[_id].frequency for _id in similarity_dict]),
                               node_id))
 
        else:
            
            total_number_of_nodes_to_analyze = len(node_list)

            def worker(data_chunk, shared_outlier_seqs_list, shared_dirty_nodes_list, shared_counter):
                for node_id in data_chunk:
                    node = self.topology.nodes[node_id]
                    shared_counter.set(shared_counter.value + 1)
                    
                    job = 'XO_%s_' % node_id
                    query, target, output = get_temporary_file_names_for_BLAST_search(prefix = job,\
                                                                                      directory = self.tmp_directory)
                    id_to_read_object_dict = {}
                    for read_obj in node.reads[1:]:
                        id_to_read_object_dict[read_obj.md5id] = read_obj

                    query_obj = u.FastaOutput(query)
                    for _id in id_to_read_object_dict:
                        query_obj.write_id(_id)
                        query_obj.write_seq(id_to_read_object_dict[_id].seq.replace('-', ''), split = False)
                    query_obj.close()

                    target_obj = u.FastaOutput(target)
                    target_obj.write_id(node.reads[0].md5id)
                    target_obj.write_seq(node.reads[0].seq.replace('-', ''), split = False)            
                    target_obj.close()

                    b = self._perform_blast(query, target, output, params = '', job = job, no_threading = True)
                    
                    similarity_dict = b.get_results_dict(max_identity = max_percent_identity)
            
                    for _id in similarity_dict:
                        outlier_read_object = id_to_read_object_dict[_id]
                        node.reads.remove(outlier_read_object)
                        shared_outlier_seqs_list.append(outlier_read_object)
                    
                    if len(similarity_dict):
                        node.dirty = True
                        shared_dirty_nodes_list.append(node)
                
                        self.logger.info('%d outliers removed from node: %s'\
                            % (sum([id_to_read_object_dict[_id].frequency for _id in similarity_dict]),
                               node_id))

            mp = Multiprocessing(worker, self.number_of_threads)
            data_chunks = mp.get_data_chunks(node_list, spiral = True)
            shared_dirty_nodes_list = mp.get_empty_shared_array()
            shared_outlier_seqs_list = mp.get_empty_shared_array()
            shared_counter = mp.get_shared_integer()

            for chunk in data_chunks:
                args = (chunk, shared_outlier_seqs_list, shared_dirty_nodes_list, shared_counter)
                mp.run(args)
        
            while 1:
                num_processes = len([p for p in mp.processes if p.is_alive()])
                
                if not num_processes:
                    # all threads are done
                    break
        
                self.progress.update('%d of %d done in %d threads' % (shared_counter.value,
                                                                      total_number_of_nodes_to_analyze,
                                                                      num_processes,))
                time.sleep(1)


            for node in shared_dirty_nodes_list:
                self.topology.nodes[node.node_id] = node
            
            for outlier_read_object in shared_outlier_seqs_list:
                self.topology.store_outlier(outlier_read_object, 'maximum_variation_allowed_reason')

        self.progress.end()

        self._refresh_topology()


    def _relocate_all_outliers(self):    
        total_relocated_outliers = 0
        
        if not self.topology.outliers:
            self.progress.end()
            self.run.info('relocated_outliers', total_relocated_outliers)
            return
        
        for reason in self.topology.outlier_reasons:
            total_relocated_outliers += self._relocate_outliers(reason, refresh_final_nodes = False)
        
        self.run.info('relocated_outliers_total', pretty_print(total_relocated_outliers))
        self._refresh_final_nodes()


    def _relocate_outliers(self, reason, refresh_final_nodes = True):
        self.progress.new('Processing %s' % get_pretty_name(reason))
        query, target, output = get_temporary_file_names_for_BLAST_search(prefix = "RO_%s_" % reason,
                                                                          directory = self.tmp_directory)

        outliers = self.topology.outliers[reason]

        id_to_read_object_dict = {}
        for read_obj in outliers:
            id_to_read_object_dict[read_obj.md5id] = read_obj

        query_obj = u.FastaOutput(query)
        for _id in id_to_read_object_dict:
            query_obj.write_id(_id)
            query_obj.write_seq(id_to_read_object_dict[_id].seq.replace('-', ''), split = False)
        query_obj.close()

        self.topology.store_node_representatives(self.topology.final_nodes, target)

        min_percent_identity = get_percent_identity_for_N_base_difference(self.topology.average_read_length,
                                                                          self.maximum_variation_allowed)
        params = "-perc_identity %.2f -max_target_seqs 1" % (min_percent_identity)
        b = self._perform_blast(query, target, output, params)

        self.progress.update('Generating similarity dict from blastn results')
        similarity_dict = b.get_results_dict(min_identity = min_percent_identity)

        outliers_relocated = len(similarity_dict)
        counter = 0
        for _id in similarity_dict:
            counter += 1
            self.progress.update('Relocating outliers: %d of %d' % (counter,
                                                                    outliers_relocated))
            self.topology.relocate_outlier(id_to_read_object_dict[_id],
                                           similarity_dict[_id].pop(),
                                           reason)

        self.progress.end()
        self.run.info('relocated_%s' % reason, pretty_print(outliers_relocated))

        if refresh_final_nodes:
            self._refresh_final_nodes()
        
        return outliers_relocated


    def _refresh_final_nodes(self):
        self.progress.new('Refreshing Nodes')
        
        for i in range(0, len(self.topology.final_nodes)):
            node = self.topology.nodes[self.topology.final_nodes[i]]
            self.progress.update('Processing %d of %d (current size: %d)' % (i + 1, len(self.topology.final_nodes), node.size))
            node.refresh()
            
        self.progress.end()
        

    def _report_final_numbers(self):
        self.run.info('num_datasets_in_fasta', pretty_print(len(self.datasets)))
        self.run.info('num_final_nodes', pretty_print(len(self.topology.final_nodes)))
        self.run.info('num_sequences_after_qc', pretty_print(self.topology.get_final_count()))

        final_outliers_total = 0
        for reason in self.topology.outlier_reasons:
            count = sum([read_obj.frequency for read_obj in self.topology.outliers[reason]])
            final_outliers_total += count
            self.run.info('final_%s' % reason, pretty_print(count))

        self.run.info('final_outliers_total', pretty_print(final_outliers_total))

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


    def _store_node_representatives(self): 
        # store representative sequences per oligotype if they are computed
        self.progress.new('Representative Sequences FASTA File')
        node_representatives_file_path = self.generate_output_destination("NODE-REPRESENTATIVES.fasta")
        f = open(node_representatives_file_path, 'w')
        self.progress.update('Being generated')
        for node_id in self.topology.final_nodes:
            node = self.topology.get_node(node_id)
            f.write('>%s|size:%d\n' % (node.node_id, node.size))
            f.write(node.representative_seq + '\n')
        f.close()
        self.progress.end()
        trim_uninformative_columns_from_alignment(node_representatives_file_path)
        self.run.info('node_representatives_file_path', node_representatives_file_path)



    def _perform_blast(self, query, target, output, params, no_threading = False, job = "NONE"):
        s = blast.LocalBLAST(query, target, output)

        s.make_blast_db()
        self.logger.info('makeblastdb for %s: %s' % (job, s.makeblastdb_cmd))
    
        if self.no_threading or no_threading:
            s.params = params 
            s.search()
            self.logger.info('blastn for %s: %s' % (job, s.search_cmd))
        else:
            s.params = params + " -num_threads %d" % (self.number_of_threads)
            s.search_parallel(self.number_of_threads, 2000)
            self.logger.info('parallel blastn for %s: %s' % (job, s.search_cmd))

        return s

    def _generate_html_output(self):
        from Oligotyping.utils.html.error import HTMLError
        try:
            from Oligotyping.utils.html.for_decomposition import generate_html_output
        except HTMLError, e:
            sys.stdout.write('\n\n\t%s\n\n' % e)
            sys.exit()

        self.progress.new('HTML Output')
        output_directory_for_html = self.generate_output_destination("HTML-OUTPUT", directory = True)
        self.progress.update('Generating')
        index_page = generate_html_output(self.run.info_dict, html_output_directory = output_directory_for_html)
        self.progress.end()
        sys.stdout.write('\n\n\tView results in your browser: "%s"\n\n' % index_page)


    def _generate_default_figures(self):
        if len(self.datasets) < 3:
            return None

        self.progress.new('Figures')

        import Oligotyping
        scripts_dir_path = os.path.dirname(Oligotyping.__file__)

        figures_dict = {}
        figures_dict['00_default'] = {}
        
        for (analysis, script, output_dir) in [('Cluster Analysis', '../Scripts/R/cluster-analysis.R', 'cluster_analysis'),
                                               ('NMDS Analysis', '../Scripts/R/metaMDS-analysis.R', 'nmds_analysis')]:
            figures_dict['00_default'][output_dir] = {}
                    
            target_dir = self.generate_output_destination('%s/__default__/%s' \
                                                                % (os.path.basename(self.figures_directory), output_dir),
                                                          directory = True)
            
            for (distance_metric, matrix_file) in [("canberra", self.matrix_percent_file_path),
                                                   ("kulczynski", self.matrix_percent_file_path),
                                                   ("jaccard", self.matrix_percent_file_path),
                                                   ("horn", self.matrix_percent_file_path),
                                                   ("chao", self.matrix_count_file_path)]:
                output_prefix = os.path.join(target_dir, distance_metric)
                cmd_line = ("%s %s %s %s %s >> %s 2>&1" % 
                                        (os.path.join(scripts_dir_path, script),
                                         matrix_file,
                                         distance_metric,
                                         self.project,
                                         output_prefix,
                                         self.log_file_path))
                self.progress.update('%s "%s" ...' % (analysis, distance_metric))
                self.logger.info('figure 00_default %s %s: %s' % (output_prefix,
                                                                  distance_metric,
                                                                  cmd_line))
                run_command(cmd_line)
                figures_dict['00_default'][output_dir][distance_metric] = output_prefix
    
        figures_dict_file_path = self.generate_output_destination("FIGURES.cPickle")
        cPickle.dump(figures_dict, open(figures_dict_file_path, 'w'))
        self.progress.end()
        self.run.info('figures_dict_file_path', figures_dict_file_path)


    def _generate_exclusive_figures(self):
        self.progress.new('Exclusive Figures')

        import Oligotyping
        scripts_dir_path = os.path.dirname(Oligotyping.__file__)
        exclusive_figures_dict = {}

        sample_mapping_dict = get_sample_mapping_dict(self.sample_mapping)
        
        for category in sample_mapping_dict:
            exclusive_figures_dict[category] = {}
            samples = sample_mapping_dict[category].keys()
            
            # double filter: first makes sure sample was not removed from the analysis due to losing all its reads during the
            # refinement, second makes sure that sample was actually mapped to something in the sample mapping file.
            samples = filter(lambda s: sample_mapping_dict[category][s], filter(lambda s: s in self.datasets, samples))
            samples.sort()

            mapping_file_path = get_temporary_file_name('%s-' % category, '-mapping.txt', self.tmp_directory)
            mapping_file = open(mapping_file_path, 'w')
            mapping_file.write('samples\t%s\n' % (category))
            
            for sample in samples:
                mapping_file.write('%s\t%s\n' % (sample, sample_mapping_dict[category][sample]))
            mapping_file.close()

            if samples == self.datasets:
                matrix_percent_path = self.matrix_percent_file_path
                matrix_count_path = self.matrix_count_file_path
            else:
                matrix_percent_path = get_temporary_file_name('%s-' % category, '-matrix-percent.txt', self.tmp_directory)
                matrix_count_path = get_temporary_file_name('%s-' % category, '-matrix-count.txt', self.tmp_directory)

                if store_filtered_matrix(self.matrix_percent_file_path, matrix_percent_path, samples) < 3:
                    self.logger.info("skipping exclusive figs for '%s'; less than 3 samples were left in MP"\
                                             % (category))
                    continue
                if store_filtered_matrix(self.matrix_count_file_path, matrix_count_path, samples) < 3:
                    self.logger.info("skipping exclusive figs for '%s'; less than 3 samples were left in MC"\
                                             % (category))
                    continue

            # ready to roll.
            self.logger.info("exclusive figs for '%s' with %d samples; mapping: '%s', MP: '%s', MC: '%s'"\
                                 % (category, len(samples), mapping_file_path, matrix_percent_path, matrix_count_path))


            for (analysis, script, output_dir) in [('NMDS Analysis', '../Scripts/R/metaMDS-analysis-with-metadata.R', 'nmds_analysis')]:
                exclusive_figures_dict[category][output_dir] = {}
                        
                target_dir = self.generate_output_destination('%s/%s/%s' % (os.path.basename(self.figures_directory),
                                                                            category,
                                                                            output_dir),
                                                              directory = True)
                
                for (distance_metric, matrix_file) in [("canberra", matrix_percent_path),
                                                       ("kulczynski", matrix_percent_path),
                                                       ("jaccard", matrix_percent_path),
                                                       ("horn", matrix_percent_path),
                                                       ("chao", matrix_count_path)]:
                    output_prefix = os.path.join(target_dir, distance_metric)
                    cmd_line = ("%s %s %s %s %s %s %s >> %s 2>&1" % 
                                            (os.path.join(scripts_dir_path, script),
                                             matrix_file,
                                             mapping_file_path,
                                             distance_metric,
                                             category,
                                             self.project,
                                             output_prefix,
                                             self.log_file_path))
                    self.progress.update('%s "%s" for "%s" ...' % (analysis, distance_metric, category))
                    self.logger.info('exclusive figure %s %s %s: %s' % (category,
                                                                        output_prefix,
                                                                        distance_metric,
                                                                        cmd_line))
                    run_command(cmd_line)
                    exclusive_figures_dict[category][output_dir][distance_metric] = output_prefix
    
        exclusive_figures_dict_file_path = self.generate_output_destination("EXCLUSIVE-FIGURES.cPickle")
        cPickle.dump(exclusive_figures_dict, open(exclusive_figures_dict_file_path, 'w'))
        self.progress.end()
        self.run.info('exclusive_figures_dict_file_path', exclusive_figures_dict_file_path)


            

if __name__ == '__main__':
    pass

