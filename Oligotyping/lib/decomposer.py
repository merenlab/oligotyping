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
import numpy
import shutil
import cPickle
import logging

from Oligotyping.lib import fastalib as u
from Oligotyping.lib.topology import Topology
from Oligotyping.lib.shared import generate_default_figures
from Oligotyping.lib.shared import generate_exclusive_figures

from Oligotyping.utils import blast
from Oligotyping.utils import utils 
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve


class Decomposer:
    def __init__(self, args = None):
        self.analysis = 'decomposition'
        self.alignment = None
        self.min_entropy = 0.0965
        self.normalize_m = True
        self.number_of_discriminants = 4
        self.min_actual_abundance = 0
        self.min_substantive_abundance = None
        self.output_directory = None
        self.project = None
        self.sample_name_separator = '_'
        self.generate_sets = False
        self.generate_frequency_curves = False
        self.skip_refining_topology = False # FIXME: ADD THIS IN PARSERS!
        self.skip_removing_outliers = False
        self.relocate_outliers = False
        self.maximum_variation_allowed = None
        self.store_topology_dict = False
        self.merge_homopolymer_splits = False
        self.no_threading = False
        self.number_of_threads = None
        self.log_file_path = None
        self.keep_tmp = False
        self.gen_html = True
        self.skip_gen_figures = False
        self.skip_check_input_file = False
        self.sample_mapping = None
        self.skip_gexf_files = False
        self.skip_basic_analyses = False
        self.quick = False
         
        if args:
            self.alignment = args.alignment
            self.min_entropy = args.min_entropy or 0.3
            self.normalize_m = not args.skip_m_normalization
            self.number_of_discriminants = args.number_of_discriminants or 3
            self.min_actual_abundance = args.min_actual_abundance
            self.min_substantive_abundance = args.min_substantive_abundance
            self.output_directory = args.output_directory
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.sample_name_separator = args.sample_name_separator
            self.generate_frequency_curves = args.generate_frequency_curves
            self.skip_removing_outliers = args.skip_removing_outliers
            self.relocate_outliers = args.relocate_outliers
            self.store_topology_dict = args.store_topology_dict
            self.merge_homopolymer_splits = args.merge_homopolymer_splits
            self.maximum_variation_allowed = args.maximum_variation_allowed
            self.no_threading = args.no_threading
            self.number_of_threads = args.number_of_threads
            self.keep_tmp = args.keep_tmp
            self.skip_gen_figures = args.skip_gen_figures
            self.skip_basic_analyses = args.skip_gen_figures
            self.skip_check_input_file = args.skip_check_input_file
            self.sample_mapping = args.sample_mapping
            self.gen_html = args.gen_html
            self.skip_gexf_files = args.skip_gexf_files
            self.quick = args.quick

        self.decomposition_depth = -1

        if self.quick:
            self.skip_gexf_files = True
            self.skip_basic_analyses = True
            self.skip_check_input_file = True
            self.skip_gen_figures = True
            self.skip_refining_topology = True
            self.skip_removing_outliers = True
            self.gen_html = False

        # there is a difference between 'average read length' and 'alignment length',
        # therefore there are two different variables to keep that information. the first
        # one is the average length of reads when gaps are removed (if there are any):
        # this may vary from organism to organism. 'alignment length', on the other hand,
        # must be the same for each read in the file.
        self.average_read_length = None
        self.alignment_length = None

        self.run = utils.Run()
        self.progress = utils.Progress()
        self.logger = None

        self.root = None
        self.topology = Topology()
        
        # A recursive method could have solved the puzzle entirely in a sample,
        # however there are a couple of reasons to not approach this problem with
        # recursion. This list is going to keep all leafs of the topology that needs
        # to be analyzed. Things will be added and removed to the list, while
        # self.topology is being formed, and main loop will check this variable its
        # every cycle.
        self.node_ids_to_analyze = None

        self.samples_dict = {}
        self.samples = []
        self.unit_counts = None
        self.unit_percents = None
        self.across_samples_sum_normalized = {}
        self.across_samples_max_normalized = {}

        # be smart, turn the threading on if necessary.
        if self.number_of_threads:
            self.no_threading = False

    def check_apps(self):
        try:
            blast.LocalBLAST(None, None, None)
        except blast.ModuleVersionError:
            raise utils.ConfigError, blast.version_error_text
        except blast.ModuleBinaryError:
            raise utils.ConfigError, blast.missing_binary_error_text

        # FIXME: check R modules here.


    def check_dirs(self):
        # check output associated stuff
        if not self.output_directory:
            self.output_directory = os.path.join(os.getcwd(), '-'.join([self.project.replace(' ', '_'), self.get_prefix()]))
        
        if not os.path.exists(self.output_directory):
            try:
                os.makedirs(self.output_directory)
            except:
                raise utils.ConfigError, "Output directory does not exist (attempt to create one failed as well): '%s'" % \
                                                                          (self.output_directory)
        if not os.access(self.output_directory, os.W_OK):
            raise utils.ConfigError, "You do not have write permission for the output directory: '%s'" % self.output_directory

        self.tmp_directory = self.generate_output_destination('TMP', directory = True)
        self.nodes_directory = self.generate_output_destination('NODES', directory = True)
        self.figures_directory = self.generate_output_destination('FIGURES', directory = True)
        self.outliers_directory = self.generate_output_destination('OUTLIERS', directory = True)


    def check_input_files(self):
        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise utils.ConfigError, "Alignment file is not accessible: '%s'" % self.alignment

        if self.sample_mapping:
            if (not os.path.exists(self.sample_mapping)) or (not os.access(self.sample_mapping, os.R_OK)):
                raise utils.ConfigError, "Sample mapping file is not accessible: '%s'" % self.sample_mapping

        samples = None
        if not self.skip_check_input_file:
            self.progress.new('Checking the input FASTA')
            samples = utils.check_input_alignment(self.alignment, self.sample_name_separator, self.progress)
            if not samples:
                raise utils.ConfigError, 'Exiting.'
            self.progress.end()

        if self.sample_mapping:
            utils.mapping_file_simple_check(self.sample_mapping, samples)
            sample_mapping_new_destination = self.generate_output_destination("SAMPLE-MAPPING.txt")
            shutil.copy(self.sample_mapping, sample_mapping_new_destination)
            self.sample_mapping = sample_mapping_new_destination


    def _init_logger(self, path = None):
        self.logger = logging.getLogger('decomposer')
        self.topology.logger = self.logger
        
        if path:
            self.log_file_path = path 
        else:
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
        
        reads = utils.get_read_objects_from_file(self.alignment)
        
        self.root = self.topology.add_new_node('root', reads, root = True)
        
        if self.root.size < self.min_actual_abundance:
            raise utils.ConfigError, "The number of reads in alignment file (%d) is smaller than --min-actual-abundance (%d)" % \
                                                                (self.root.size, self.min_actual_abundance)

        self.node_ids_to_analyze = ['root']

        self.progress.end()

            
    def get_prefix(self):
        prefix = 'm%.2f-A%d-M%d-d%d' % (self.min_entropy,
                                        self.min_actual_abundance,
                                        self.min_substantive_abundance,
                                        self.number_of_discriminants)

        return prefix


    def set_min_substantive_abundance(self):
        suggested_M = round(self.topology.nodes['root'].size / 5000.0)
        self.min_substantive_abundance = suggested_M if suggested_M > 4 else 4


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
        
        if self.sample_mapping:
            self.sample_mapping_dict = utils.get_sample_mapping_dict(self.sample_mapping)
        else:
            self.sample_mapping_dict = None

        # we're in business.
        self.run.info('project', self.project)
        self.run.info('run_date', utils.get_date())
        self.run.info('version', __version__)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))
        self.run.info('multi_threaded', not self.no_threading)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('log_file_path', self.log_file_path)
        self.run.info('root_alignment', self.alignment)
        self.run.info('sample_mapping', self.sample_mapping)
        self.run.info('quick', self.quick)
        self.run.info('merge_homopolymer_splits', self.merge_homopolymer_splits)
        self.run.info('skip_removing_outliers', self.skip_removing_outliers)
        self.run.info('relocate_outliers', self.relocate_outliers)
        self.run.info('store_topology_dict', self.store_topology_dict)
        self.run.info('skip_gen_figures', self.skip_gen_figures)
        self.run.info('m', self.min_entropy)
        self.run.info('normalize_m', self.normalize_m)
        self.run.info('d', self.number_of_discriminants)
        self.run.info('A', self.min_actual_abundance)

        # set number of threads to be used
        if not self.number_of_threads:
            self.number_of_threads = utils.Multiprocessing(None).num_thread

        # If --min-substantive abundance is not zero, use it. otherwise a default value is
        # going to be computed for it after the topology is initiated and the number of reads
        # in the dataset is known.
        if self.min_substantive_abundance:
            self.run.info('M', self.min_substantive_abundance)

        # initializing the topology.
        self._init_topology()

        if not self.min_substantive_abundance:
            self.set_min_substantive_abundance()
            self.run.info('M', self.min_substantive_abundance)


        # to decide at what level should algorithm be concerned about divergent reads in a node, there has to be a
        # threshold that defines what is the maximum variation from the most abundant unique sequence. if user did
        # not define this value, it is being set here as follows (FIXME: this is a very crude way to do it): 
        if not self.maximum_variation_allowed:
            self.maximum_variation_allowed = int(round(self.topology.average_read_length * 1.0 / 100)) or 1

        self.run.info('maximum_variation_allowed', self.maximum_variation_allowed)
        self.run.info('total_seq', utils.pretty_print(self.topology.nodes['root'].size))
        self.run.info('average_read_length', self.topology.average_read_length)
        self.run.info('alignment_length', self.topology.alignment_length)
        self.run.info('output_directory', self.output_directory)
        self.run.info('nodes_directory', self.nodes_directory)
        self.run.info('tmp_directory', self.tmp_directory)
        self.run.info('figures_directory', self.figures_directory)

        # business time.
        self._generate_raw_topology()

        if not self.skip_refining_topology:
            self._refine_topology()
           
        if self.relocate_outliers:
            self._relocate_all_outliers()

        self._generate_samples_dict()
        self._get_unit_counts_and_percents()

        # all done.        
        self._report_final_numbers()
         
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()

        if self.generate_frequency_curves:
            self._generate_frequency_curves()
 
        self._store_topology()
        self._store_final_nodes()
        self._store_all_outliers()
        self._store_node_representatives()
        self._store_read_distribution_table()
        
        if self.store_topology_dict:
            self._store_topology_dict()

        for node_id in self.topology.final_nodes:
            self.logger.info('final node: %s (%d)' % (node_id,
                                                      self.topology.nodes[node_id].size))

        self.run.info('end_of_run', utils.get_date())

        if not self.skip_gexf_files:
            self._generate_gexf_network_file()

        if not self.skip_gen_figures:
            self._generate_default_figures()
        
        if (not self.skip_gen_figures) and self.sample_mapping:
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
                        raise utils.ConfigError, "Number of unique reads in the root node (%d) is less than the declared minimum (%d)." \
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

                if node.competing_unique_sequences_ratio < 0.0005 or node.density > 0.85:
                    # Finalize this node.
                    self.logger.info('finalize node (CUSR/ND): %s' % node_id)
                    continue

                # find out about the entropy distribution in the given node:
                node.do_entropy()
                
                # normalize m if the user hasn't opted out.
                if self.normalize_m:
                    node.set_normalized_m(self.min_entropy, self.topology.frequency_of_the_most_abundant_read)
                    self.logger.info('normalized m (NM) for %s: %.3f ' % (node_id, node.normalized_m))

                p += ' / ME: %.2f / AE: %.2f / NM: %s' % (max(node.entropy),
                                                          node.average_entropy,
                                                          ('%.3f' % node.normalized_m) if self.normalize_m else None)
                self.progress.update(p)
                
                # IF the abundance of the second most abundant unique read in the node is smaller than 
                # the self.min_substantive_abundance criteria, there is no need to further decompose
                # this node. because anything spawns from here, will end up in the outlier bin except
                # the most abundant unique read. of course by not decomposing any further we are losing
                # the opportunity to 'purify' this node further, but we are not worried about it,
                # because 'max_allowed_variation' outliers will be removed from this node later on.  
                # UPDATE: Well, this causes some serious purity issues. For instance a node with,
                # 
                # >Read_1|frequency:957
                # >Read_2|frequency:120
                # >Read_3|frequency:57
                # >Read_4|frequency:7
                #
                # is finalized due to SMA < MSA although the entropy looked like this:
                #
                #    http://i.imgur.com/ctFnJE2.png
                #
                # when M = 300. This begs for a FIXME.
                #
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
                if self.normalize_m:
                    node.discriminants = [d[0] for d in node.entropy_tpls[0:self.number_of_discriminants] if d[1] > node.normalized_m]
                else:
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
                
                
                # all reads in the parent node are analyzed. time to add spawned nodes into the topology.
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

        self.run.info('num_raw_nodes', utils.pretty_print(len(self.topology.final_nodes)))

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

            mp =  utils.Multiprocessing(worker, self.number_of_threads)
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
                    self.run.info('removed_%s' % reason, utils.pretty_print(count))
                self.run.info('removed_outliers_total', utils.pretty_print(removed_outliers_total))

                break

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
                # a node, but were betrayed by the nature of the sample or the sequencing platform. either due to
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
       
        nz = utils.pretty_print(len(self.topology.zombie_nodes))
        self.progress.new('Merging HP splits :: ITER %d%s' % (iteration,
                                                              ' #Z: %s' % (nz if nz else '')))

        self.progress.update('Running blastn')
        dealing_with_zombie_nodes = False

        if self.topology.zombie_nodes:
            nodes = copy.deepcopy(self.topology.zombie_nodes)
            dealing_with_zombie_nodes = True
        else:
            nodes = copy.deepcopy(self.topology.final_nodes)

        # ---
        # use blastn to get the similarity dict for gaps.
        query, target, output = utils.get_temporary_file_names_for_BLAST_search(prefix = "HPS_%d_" % iteration,\
                                                                          directory = self.tmp_directory)

        self.topology.store_node_representatives(nodes, query)
        self.topology.store_node_representatives(self.topology.final_nodes, target)

        min_percent_identity = utils.get_percent_identity_for_N_base_difference(self.topology.average_read_length) - 1
        params = "-perc_identity %.2f" % (min_percent_identity)
        b = self._perform_blast(query, target, output, params, job = 'HPS', no_threading = True)
        
        self.progress.update('Generating similarity dict from blastn results')
        similarity_dict = b.get_results_dict(mismatches = 0, gaps = 1)
        
        node_ids = set(similarity_dict.keys())

        # we first generate list of lists that, each of which contains node ids that are
        # split by a homopolymer region-associated in/del (HPS).
        merge_clusters = []
        for source_node_id in node_ids:
            related_ids = [source_node_id]
            node = self.topology.nodes[source_node_id]
            
            for target_node_id in similarity_dict[source_node_id]:
                sibling = self.topology.nodes[target_node_id]
                if utils.homopolymer_indel_exists(node.representative_seq, sibling.representative_seq):
                    related_ids.append(target_node_id)

            if len(related_ids) > 1:
                merge_clusters.append(set(related_ids))

        # previous step generates many clusters that are connected by shared nodes, such as this one:
        #
        # [set([000000007', '000000006']), set(['000000005', '000000006'])]
        #
        # so we finalize it here by absorbing connected clusters into one cluster, which
        # generates this from the one above:
        #
        # set(['000000005', '000000007', '000000006'])]
        #
        hps_merge_final = []
        while len(merge_clusters):
            cluster_absorbing = merge_clusters.pop()

            redo = True
            while redo:
                redo = False
                for cluster in merge_clusters:
                    if cluster_absorbing.intersection(cluster):
                        cluster_absorbing.update(cluster)
                        merge_clusters.remove(cluster)
                        redo = True

            hps_merge_final.append(cluster_absorbing)

        # code above populates merge_clusters, which looks like this at the end:
        #
        #    [set(['000000009', '000000002']), set(['000000005', '000000007', '000000006'])]
        #
        # the one liner below sorts these ladies appropriately. so it is always the lower abundance
        # is merged with the higher abundance:
        #
        #    [[(4, '000000009'), (2, '000000002')], [(5, '000000007'), (4, '000000005'), (3, '000000006')]]
        hps_merge_final = [sorted([(self.topology.nodes[t].size, t) for t in x], reverse = True) for x in hps_merge_final]
        self.logger.info('merge clusters: %s' % hps_merge_final.__str__())

        # go through every merge cluster, perform merging the most left node
        for i in range(0, len(hps_merge_final)):
            merge_cluster = hps_merge_final[i]
            self.progress.update('Processing merge cluster %d of %d' % (i, len(merge_clusters)))
            
            node_id = merge_cluster[0][1]
            
            for sibling_id in [m[1] for m in merge_cluster[1:]]:

                if dealing_with_zombie_nodes:
                    self.logger.info('zombie node merged (HPS): %s <<< %s' % (node_id, sibling_id))
                    self.topology.merge_nodes(node_id, sibling_id)
                    self.topology.standby_bin.append(node_id)
                else:
                    self.logger.info('nodes merged (HPS): %s <<< %s' % (node_id, sibling_id))
                    self.topology.merge_nodes(node_id, sibling_id)


        # reset temporary stuff
        self.samples = []
        self.samples_dict = {}
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


        sb = utils.pretty_print(len(self.topology.standby_bin))
        self.progress.new('Removing Outliers :: ITER %d%s' % (iteration,
                                                              ' #SB: %s' % (sb if sb else '')))

        if standby_bin_only:
            node_list = copy.deepcopy(self.topology.standby_bin)
            self.topology.standby_bin = []
        else:
            node_list = copy.deepcopy(self.topology.final_nodes)

        
        min_percent_identity = utils.get_percent_identity_for_N_base_difference(self.topology.average_read_length,
                                                                          self.maximum_variation_allowed)
        param = "-perc_identity %.2f" % (min_percent_identity)

        if self.no_threading:
            # no threading
            for i in range(0, len(node_list)):
                node_id = node_list[i]
                node = self.topology.nodes[node_id]
                
                self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.topology.final_nodes)))

                job = 'XO_%s_' % node_id
                query, target, output = utils.get_temporary_file_names_for_BLAST_search(prefix = job,\
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

                b = self._perform_blast(query, target, output, params = param, job = job)

                # while I feel ashamed for this redundancy, you go ahead and read the description
                # written in the worker function below on this voodoo stuff.
                similarity_dict = b.get_results_dict(min_identity = min_percent_identity)
                read_ids_to_keep = set(similarity_dict.keys())
                all_read_ids = set(id_to_read_object_dict.keys())
                outliers = all_read_ids.difference(read_ids_to_keep)

                if len(outliers):
                    node.dirty = True
                else:
                    continue
    
                self.progress.append(' / screening node to remove %d outliers' % len(outliers))

                for _id in outliers:
                    outlier_read_object = id_to_read_object_dict[_id]
                    node.reads.remove(outlier_read_object)
                    self.topology.store_outlier(outlier_read_object, 'maximum_variation_allowed_reason')
 
                self.logger.info('%d outliers removed from node: %s'\
                            % (sum([id_to_read_object_dict[_id].frequency for _id in outliers]),
                               node_id))
 
        else:
            
            def worker(node_id, shared_outlier_seqs_list, shared_dirty_nodes_list):
                node = self.topology.nodes[node_id]
                    
                job = 'XO_%s_' % node_id
                query, target, output = utils.get_temporary_file_names_for_BLAST_search(prefix = job,\
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

                b = self._perform_blast(query, target, output, params = param, job = job, no_threading = True)

                # something semi-smart: get all the read ids that are more similar to the rep_seq
                # than allowed max_variation; keep them, remove anything that doesn't show up here.
                # the other option would be to search for low similarity guys, but it would have
                # required much more computational investment. 
                similarity_dict = b.get_results_dict(min_identity = min_percent_identity)
                read_ids_to_keep = set(similarity_dict.keys())
                all_read_ids = set(id_to_read_object_dict.keys())
                outliers = all_read_ids.difference(read_ids_to_keep)

                for _id in outliers:
                    outlier_read_object = id_to_read_object_dict[_id]
                    node.reads.remove(outlier_read_object)
                    shared_outlier_seqs_list.append(outlier_read_object)
                
                if len(outliers):
                    node.dirty = True
                    shared_dirty_nodes_list.append(node)
                
                    self.logger.info('%d outliers removed from node: %s (max frequency: %d; mean frequency: %.2f)'\
                        % (sum([id_to_read_object_dict[_id].frequency for _id in outliers]),
                           node_id,
                           max([id_to_read_object_dict[_id].frequency for _id in outliers]),
                           numpy.mean([id_to_read_object_dict[_id].frequency for _id in outliers]),))

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            shared_dirty_nodes_list = mp.get_empty_shared_array()
            shared_outlier_seqs_list = mp.get_empty_shared_array()

            # arrange processes
            processes_to_run = []
            for node in node_list:
                processes_to_run.append((node, shared_outlier_seqs_list, shared_dirty_nodes_list),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)

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
        
        self.run.info('relocated_outliers_total', utils.pretty_print(total_relocated_outliers))
        self._refresh_final_nodes()


    def _relocate_outliers(self, reason, refresh_final_nodes = True):
        self.progress.new('Processing %s' % utils.get_pretty_name(reason))
        query, target, output = utils.get_temporary_file_names_for_BLAST_search(prefix = "RO_%s_" % reason,
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

        min_percent_identity = utils.get_percent_identity_for_N_base_difference(self.topology.average_read_length,
                                                                          self.maximum_variation_allowed)

        self.progress.update('Running blastn (query: %s, target: %s)' % (utils.pretty_print(len(outliers)),
                                                                         utils.pretty_print(len(self.topology.final_nodes))))
        params = "-perc_identity %.2f -max_target_seqs 1" % (min_percent_identity)
        b = self._perform_blast(query, target, output, params, job = 'RO_%s_' % reason)

        self.progress.update('Generating similarity dict from blastn results')
        similarity_dict = b.get_results_dict(min_identity = min_percent_identity)

        num_outlier_objects = len(similarity_dict)
        num_outliers_relocated = sum([id_to_read_object_dict[_id].frequency for _id in similarity_dict])

        counter = 0
        for _id in similarity_dict:
            counter += 1
            self.progress.update('Relocating outliers: %d of %d' % (counter,
                                                                    num_outlier_objects))
            self.topology.relocate_outlier(id_to_read_object_dict[_id],
                                           similarity_dict[_id].pop(),
                                           reason)

        self.progress.end()
        self.run.info('relocated_%s' % reason, utils.pretty_print(num_outliers_relocated))

        if refresh_final_nodes:
            self._refresh_final_nodes()
        
        return num_outliers_relocated


    def _refresh_final_nodes(self):
        self.progress.new('Refreshing Nodes')
        
        for i in range(0, len(self.topology.final_nodes)):
            node = self.topology.nodes[self.topology.final_nodes[i]]
            self.progress.update('Processing %d of %d (current size: %d)' % (i + 1, len(self.topology.final_nodes), node.size))
            node.refresh()
            
        self.progress.end()
        

    def _generate_frequency_curves(self):
        self.progress.new('Generating mini entropy figures')
        for i in range(0, len(self.topology.alive_nodes)):
            node = self.topology.nodes[self.topology.alive_nodes[i]]
            
            if node.killed:
                continue
            
            self.progress.update('Node ID: "%s" (%d of %d)' % (node.pretty_id, i + 1, len(self.topology.alive_nodes)))
                                 
            node.freq_curve_img_path = node.unique_alignment + '.png'
            vis_freq_curve(node.unique_alignment, output_file = node.freq_curve_img_path, mini = True,\
                           title = '%s\n(%s)' % (node.pretty_id, utils.human_readable_number(node.size))) 
        
        self.progress.end()


    def _store_read_distribution_table(self):
        self.progress.new('Read distribution table')
        self.read_distribution_table_path = self.generate_output_destination("READ-DISTRIBUTION.txt")

        def get_dict_entry_tmpl():
            d = {'represented_reads': 0}
            for reason in self.topology.outlier_reasons:
                d[reason] = 0
            return d

        read_distribution_dict = {}
        
        self.progress.update('Processing reads that were represented in results')
        for sample in self.samples_dict:
            if not read_distribution_dict.has_key(sample):
                read_distribution_dict[sample] = get_dict_entry_tmpl()

            read_distribution_dict[sample]['represented_reads'] = sum(self.samples_dict[sample].values())
            
        for reason in self.topology.outlier_reasons:
            self.progress.update('Processing outliers (%s)' % (reason))
            for read_object in self.topology.outliers[reason]:
                for read_id in read_object.ids:
                    sample = utils.get_sample_name_from_defline(read_id, self.sample_name_separator)
                    
                    if not read_distribution_dict.has_key(sample):
                        read_distribution_dict[sample] = get_dict_entry_tmpl()
                
                    read_distribution_dict[sample][reason] += 1
        
        self.progress.update('Storing...')
        utils.generate_TAB_delim_file_from_dict(read_distribution_dict,
                                          self.read_distribution_table_path,
                                          order = ['represented_reads'] + self.topology.outlier_reasons)

        self.progress.end()
        self.run.info('read_distribution_table_path', self.read_distribution_table_path)


    def _store_final_nodes(self):
        self.progress.new('Storing final nodes')

        total_final_nodes = len(self.topology.final_nodes)
        
        if self.no_threading:
            for i in range(0, total_final_nodes):
                self.progress.update('%s of %s' % (utils.pretty_print(i + 1),
                                                   utils.pretty_print(total_final_nodes)))
                node_id = self.topology.final_nodes[i]
                node = self.topology.get_node(node_id)
                node.store()

        else:
            
            def worker(node_id):
                node = self.topology.get_node(node_id)
                node.store()

            mp = utils.Multiprocessing(worker, self.number_of_threads)
            
            # arrange processes
            processes_to_run = []
            for node_id in self.topology.final_nodes:
                processes_to_run.append((node_id,),)

            # start the main loop to run all processes
            mp.run_processes(processes_to_run, self.progress)

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


    def _generate_samples_dict(self):
        self.progress.new('Computing Samples Dict')
        
        for node_id in self.topology.final_nodes:
            node = self.topology.nodes[node_id]
            self.progress.update('Analyzing Node ID: "%s" (size: %d)'\
                                                        % (node_id, node.size))
        
            for read in node.reads:
                for read_id in read.ids:
                    sample = utils.get_sample_name_from_defline(read_id, self.sample_name_separator)
            
                    if not self.samples_dict.has_key(sample):
                        self.samples_dict[sample] = {}
                        self.samples.append(sample)
    
                    if self.samples_dict[sample].has_key(node_id):
                        self.samples_dict[sample][node_id] += 1
                    else:
                        self.samples_dict[sample][node_id] = 1

        self.samples.sort()
        self.progress.end()


    def _generate_ENVIRONMENT_file(self):
        self.progress.new('ENVIRONMENT File')
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        self.progress.update('Being generated')
        
        utils.generate_ENVIRONMENT_file(self.samples,
                                        self.samples_dict,
                                        environment_file_path)

        self.progress.end()
        self.run.info('environment_file_path', environment_file_path)        


    def _get_unit_counts_and_percents(self):
        self.progress.new('Unit counts and percents')
        self.progress.update('Data is being generated')
            
        self.unit_counts, self.unit_percents = utils.get_unit_counts_and_percents(self.topology.final_nodes, self.samples_dict)
            
        self.progress.end()


    def _generate_MATRIX_files(self):
        self.progress.new('Matrix Files')
        self.progress.update('Being generated')
            
        self.matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        self.matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")    
            
        utils.generate_MATRIX_files(self.topology.final_nodes,
                                    self.samples,
                                    self.unit_counts,
                                    self.unit_percents,
                                    self.matrix_count_file_path,
                                    self.matrix_percent_file_path)
            
        self.progress.end()
        self.run.info('matrix_count_file_path', self.matrix_count_file_path)
        self.run.info('matrix_percent_file_path', self.matrix_percent_file_path)


    def _store_topology_dict(self):
        self.progress.new('Generating topology dict (lightweight)')
        topology_dict = {}
        
        self.progress.update('Processing %d nodes' % len(self.topology.nodes))
        for node_id in self.topology.nodes:
            node = self.topology.nodes[node_id]
            if node.killed:
                continue

            new_node = copy.deepcopy(node)
            new_node.entropy_tpls = None
            
            topology_dict[node_id] = new_node

        self.progress.end()
        
        topology_dict_file_path = self.generate_output_destination('TOPOLOGY-LIGHT.cPickle')
        cPickle.dump(topology_dict, open(topology_dict_file_path, 'w'))
        self.run.info('topology_light_dict', topology_dict_file_path)


    def _store_topology(self):
        self.progress.new('Generating output files for topology')
        topology_text_file_path = self.generate_output_destination('TOPOLOGY.txt')
        topology_text_file_obj = open(topology_text_file_path, 'w')
        
        nodes_dict = {}
        for node_id in self.topology.alive_nodes:
            node = self.topology.nodes[node_id]
            topology_text_file_obj.write('%s\t%d\t%s\t%d\t%s\n' \
                                               % (node.node_id,
                                                  node.size,
                                                  node.parent or '',
                                                  node.level,
                                                  ','.join(node.children) or ''))

            nodes_dict[node_id] = {'size': node.size,
                                   'level': node.level,
                                   'parent': node.parent,
                                   'children': ','.join(node.children) or None,
                                   'final_node': 'Yes' if node.children else 'No',
                                   'max_entropy': node.max_entropy,
                                   'average_entropy': node.average_entropy,
                                   'density': node.density,
                                   'num_comps_larger_than_m': len([True for tpl in node.entropy_tpls if tpl[1] > self.min_entropy]),
                                   'normalized_m': node.normalized_m}


        topology_text_file_obj.close()

        if not self.skip_gexf_files:
            topology_gexf_file_path = self.generate_output_destination('TOPOLOGY.gexf')
            attribute_types_dict = {'size': "int",
                                    'level': "int",
                                    'max_entropy': "float",
                                    'average_entropy': "float",
                                    'density': "float",
                                    'num_comps_larger_than_m': "int",
                                    'normalized_m': "float"}

            utils.generate_gexf_network_file_for_nodes_topology(nodes_dict,
                                                                topology_gexf_file_path,
                                                                attribute_types_dict = attribute_types_dict)

        self.progress.end()

        self.run.info('topology_text', topology_text_file_path)

        if not self.skip_gexf_files:
            self.run.info('topology_gexf', topology_gexf_file_path)


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
        utils.trim_uninformative_columns_from_alignment(node_representatives_file_path)
        self.run.info('node_representatives_file_path', node_representatives_file_path)



    def _perform_blast(self, query, target, output, params, no_threading = False, job = "NONE"):
        s = blast.LocalBLAST(query, target, output, log = self.generate_output_destination('BLAST.log'))
        self.logger.info('local blast request for job "%s": (q) %s (t) %s (o) %s (p) %s (th) %s'\
                                               % (job, query, target, output, params, not no_threading))

        s.make_blast_db()
        self.logger.info('makeblastdb for %s: %s' % (job, s.makeblastdb_cmd))

        # OK. this is what it was like for params within the else clause below:
        #
        #    s.params = params + " -num_threads %d" % (self.number_of_threads) 
        #
        # I think there is a problem with asking blastn to run multi-threaded in the context
        # of this pipeline, since it handles multiprocessing by itself (even though I know it sounds
        # funny, my home-made poor man's parallel blast is faster than blast with -num_threads option).
        # if search_parallel is being called, blast process shouldn't be bothered with -num_threads.
        # because I believe it makes everything extremely slower (since N blast search start each with
        # N threads, we exceed the number of cores quickly and the overhead from the scheduler kills it).
        # no_threading is only coming from remove_outliers, and the process that calls with no_threading True
        # is already multi-threaded there as well. so either of these cases require an extra -num_threads
        # directive to blastn. therefore I removed it from both (I'll do more tests, I might put it back). 

        if self.no_threading or no_threading:
            s.params = params
            s.search()
            self.logger.info('blastn for %s: %s' % (job, s.search_cmd))
        else:
            s.params = params
            s.search_parallel(self.number_of_threads, 2000, keep_parts = self.keep_tmp)
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


    def _generate_gexf_network_file(self):
        self.gexf_network_file_path = self.generate_output_destination("NETWORK.gexf")

        self.progress.new('GEXF Network File')
       
        utils.generate_gexf_network_file(self.topology.final_nodes,
                                         self.samples_dict,
                                         self.unit_percents,
                                         self.gexf_network_file_path,
                                         sample_mapping_dict = self.sample_mapping_dict,
                                         project = self.project)

        self.progress.end()
        self.run.info('gexf_network_file_path', self.gexf_network_file_path)


    def _generate_default_figures(self):
        if len(self.samples) < 3:
            return None

        self.progress.new('Figures')

        figures_dict = generate_default_figures(self)
        figures_dict_file_path = self.generate_output_destination("FIGURES.cPickle")
        cPickle.dump(figures_dict, open(figures_dict_file_path, 'w'))

        self.progress.end()
        self.run.info('figures_dict_file_path', figures_dict_file_path)


    def _generate_exclusive_figures(self):
        if len(self.samples) < 3:
            return None

        self.progress.new('Exclusive Figures')

        exclusive_figures_dict = generate_exclusive_figures(self)
        exclusive_figures_dict_file_path = self.generate_output_destination("EXCLUSIVE-FIGURES.cPickle")
        cPickle.dump(exclusive_figures_dict, open(exclusive_figures_dict_file_path, 'w'))

        self.progress.end()
        self.run.info('exclusive_figures_dict_file_path', exclusive_figures_dict_file_path)

            
    def _report_final_numbers(self):
        self.run.info('num_samples_in_fasta', utils.pretty_print(len(self.samples)))
        self.run.info('num_final_nodes', utils.pretty_print(len(self.topology.final_nodes)))
        self.run.info('num_sequences_after_qc', utils.pretty_print(self.topology.get_final_count()))

        final_outliers_total = 0
        for reason in self.topology.outlier_reasons:
            count = sum([read_obj.frequency for read_obj in self.topology.outliers[reason]])
            final_outliers_total += count
            self.run.info('final_%s' % reason, utils.pretty_print(count))

        self.run.info('final_outliers_total', utils.pretty_print(final_outliers_total))

if __name__ == '__main__':
    pass

