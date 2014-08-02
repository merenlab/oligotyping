#!/usr/bin/python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

__version__ = '1.0'

import os
import sys
import copy
import shutil
import cPickle
import logging
import itertools
import math

from Oligotyping.utils import utils
from Oligotyping.utils import blast
from Oligotyping.utils.random_colors import random_colors
from Oligotyping.utils.random_colors import get_color_shade_dict_for_list_of_values
from Oligotyping.lib import fastalib as u
from Oligotyping.lib.shared import generate_default_figures
from Oligotyping.lib.shared import generate_exclusive_figures
from Oligotyping.visualization.frequency_curve_and_entropy import vis_freq_curve
from Oligotyping.visualization.oligotype_sets_distribution import vis_oligotype_sets_distribution
from Oligotyping.visualization.oligotype_network_structure import oligotype_network_structure
from Oligotyping.visualization.oligotype_distribution_stack_bar import oligotype_distribution_stack_bar
from Oligotyping.visualization.oligotype_distribution_across_samples import oligotype_distribution_across_samples


class Oligotyping:
    def __init__(self, args = None):
        self.analysis = 'oligotyping'
        self.entropy   = None
        self.alignment = None
        self.quals_dict = None
        self.min_base_quality = None
        self.number_of_auto_components = 5
        self.selected_components = None
        self.limit_oligotypes_to = None
        self.exclude_oligotypes = None
        self.min_number_of_samples = 1
        self.min_percent_abundance = 0.0
        self.min_actual_abundance = 0
        self.min_substantive_abundance = 4
        self.project = None
        self.output_directory = None
        self.sample_name_separator = '_'
        self.limit_representative_sequences = sys.maxint
        self.quick = False
        self.no_figures = False
        self.no_display = False
        self.keep_tmp = False
        self.blast_ref_db = None
        self.skip_blast_search = False
        self.gen_html = False
        self.gen_sample_oligo_networks = False
        self.colors_list_file = None
        self.generate_sets = False
        self.cosine_similarity_threshold = 0.1
        self.sample_mapping = None
        self.log_file_path = None
        self.skip_check_input_file = False
        self.skip_basic_analyses = False
        self.skip_gexf_network_file = False
        self.no_threading = False
        self.number_of_threads = None

        Absolute = lambda x: os.path.join(os.getcwd(), x) if not x.startswith('/') else x 

        if args:
            self.entropy = Absolute(args.entropy)
            self.alignment = Absolute(args.alignment)
            self.quals_dict = utils.process_command_line_args_for_quality_files(args, _return = 'quals_dict')
            self.min_base_quality = args.min_base_quality
            self.number_of_auto_components = args.number_of_auto_components
            self.selected_components = args.selected_components
            self.limit_oligotypes_to = args.limit_oligotypes_to
            self.exclude_oligotypes = args.exclude_oligotypes
            self.min_number_of_samples = args.min_number_of_samples
            self.min_percent_abundance = args.min_percent_abundance
            self.min_actual_abundance = args.min_actual_abundance
            self.min_substantive_abundance = args.min_substantive_abundance
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.output_directory = args.output_directory
            self.sample_name_separator = args.sample_name_separator
            self.limit_representative_sequences = args.limit_representative_sequences or sys.maxint
            self.quick = args.quick
            self.no_figures = args.no_figures
            self.no_display = args.no_display
            self.keep_tmp = args.keep_tmp
            self.blast_ref_db = Absolute(args.blast_ref_db) if args.blast_ref_db else None
            self.skip_blast_search = args.skip_blast_search
            self.gen_html = args.gen_html
            self.gen_sample_oligo_networks = args.gen_sample_oligo_networks
            self.colors_list_file = args.colors_list_file
            self.cosine_similarity_threshold = args.cosine_similarity_threshold
            self.generate_sets = args.generate_sets
            self.sample_mapping = args.sample_mapping
            self.skip_check_input_file = args.skip_check_input_file
            self.skip_basic_analyses = args.skip_basic_analyses
            self.skip_gexf_network_file = args.skip_gexf_network_file
            self.no_threading = args.no_threading
            self.number_of_threads = args.number_of_threads
        
        self.run = utils.Run()
        self.progress = utils.Progress()

        self.samples_dict = {}
        self.sample_mapping_dict = {}
        self.excluded_read_ids_tracker = {}
        self.representative_sequences_per_oligotype = {}
        self.across_samples_sum_normalized = {}
        self.across_samples_max_normalized = {}
        self.unit_counts = None
        self.unit_percents = None
        self.oligotype_sets = None
        self.samples = []
        self.abundant_oligos = []

        self.final_oligo_counts_dict = {}
        self.final_oligo_entropy_distribution_dict = {}
        self.final_oligo_unique_distribution_dict = {}
        self.final_purity_score_dict = {}
        self.total_purity_score_dict = {}
        self.colors_dict = None
        self.score_color_dict = None
        self.figures_directory = None

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
        if self.number_of_auto_components != None and self.selected_components != None:
            raise utils.ConfigError, "You either have to declare 'auto components' (-c) or 'selected components' (-C)."

        if self.number_of_auto_components == None and self.selected_components == None:
            raise utils.ConfigError, "Both 'auto components' (-c), and 'selected components' (-C) were declared."

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
        self.figures_directory = self.generate_output_destination('FIGURES', directory = True)


    def check_input(self):
        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise utils.ConfigError, "Alignment file is not accessible: '%s'" % self.alignment
        
        if (not os.path.exists(self.entropy)) or (not os.access(self.entropy, os.R_OK)):
            raise utils.ConfigError, "Entropy file is not accessible: '%s'" % self.entropy

        if self.sample_mapping:
            if (not os.path.exists(self.sample_mapping)) or (not os.access(self.sample_mapping, os.R_OK)):
                raise utils.ConfigError, "Sample mapping file is not accessible: '%s'" % self.sample_mapping

        if self.colors_list_file:
            if not os.path.exists(self.colors_list_file):
                raise utils.ConfigError, "Colors list file does not exist: '%s'" % self.colors_list_file
            first_characters = list(set([c.strip()[0] for c in open(self.colors_list_file)]))
            if len(first_characters) != 1 or first_characters[0] != '#':
                raise utils.ConfigError, "Colors list file does not seem to be correctly formatted"

        # set the alignment lentgh (it will be necessary to check certain params)
        alignment = u.SequenceSource(self.alignment)
        alignment.next()
        self.alignment_length = len(alignment.seq)
        alignment.close()

        # now we know that input files are OK, lets check input params before we go any further.
        self.check_params()

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


    def check_params(self):
        if self.selected_components:
            try:
                self.selected_components = [int(c) for c in self.selected_components.split(',')]
            except:
                raise utils.ConfigError, "Selected components should be comma separated integer values (such as '4,8,15,25,47')."
        
            if max(self.selected_components) >= self.alignment_length:
                raise utils.ConfigError, "There is at least one component ('%d') that is bigger than the alignment length."\
                                                                        % max(self.selected_components) 
        
            if min(self.selected_components) < 0:
                raise utils.ConfigError, "Selected components can't be smaller than 0"

            components_declared_more_than_once = [c[0] for c in itertools.groupby(sorted(self.selected_components))\
                                                                        if len(list(c[1])) > 1]
            N = len(components_declared_more_than_once)
            if N:
                raise utils.ConfigError, "You declared %s component%s (%s) more than once."\
                                             % ('a' if N == 1 else '%s' % str(N), 
                                                's' if N > 1 else '',
                                                ', '.join([str(c) for c in components_declared_more_than_once]))

        if self.min_base_quality:
            try:
                self.min_base_quality = int(self.min_base_quality)
                assert(self.min_base_quality >= 0 and self.min_base_quality <= 40)
            except:
                raise utils.ConfigError, "Minimum base quality must be an integer between 0 and 40."

        if self.limit_oligotypes_to:
            self.limit_oligotypes_to = [o.strip().upper() for o in self.limit_oligotypes_to.split(',')]
            if len(self.limit_oligotypes_to) == 1:
                raise utils.ConfigError, "There must be more than one oligotype for --limit-oligotypes parameter."

            if len([n for n in ''.join(self.limit_oligotypes_to) if n not in ['A', 'T', 'C', 'G', '-']]):
                raise utils.ConfigError, "Oligotypes defined by --limit-oligotypes parameter seems to have ambiguous characters."

        if self.exclude_oligotypes:
            self.exclude_oligotypes = [o.strip().upper() for o in self.exclude_oligotypes.split(',')]
            
            if len([n for n in ''.join(self.exclude_oligotypes) if n not in ['A', 'T', 'C', 'G', '-']]):
                raise utils.ConfigError, "Oligotypes defined by --exclude-oligotypes parameter seems to have ambiguous characters."
            
        return True


    def _init_logger(self, path = None):
        self.logger = logging.getLogger('oligotyping')
        
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


    def get_prefix(self):
        prefix = 's%d-a%.1f-A%d-M%d' % (self.min_number_of_samples,
                                        self.min_percent_abundance,
                                        self.min_actual_abundance,
                                        self.min_substantive_abundance)

        if self.selected_components:
           
            # I don't have any desire to solve dependencies of the initialization steps properly, so
            # please have a cup of ugly hack:
            if type(self.selected_components) == type(''):
                num_sc = len(self.selected_components.split(','))
            else:
                num_sc = len(self.selected_components)

            prefix = 'sc%d-%s' % (num_sc, prefix)
        else:
            prefix = 'c%d-%s' % (self.number_of_auto_components, prefix)
        
        if self.quals_dict:
            prefix = '%s-q%d' % (prefix, self.min_base_quality)

        return prefix


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def run_all(self):
        self.check_apps()
        self.check_dirs()

        # ready to init logging        
        self._init_logger()
        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.run.init_info_file_obj(self.info_file_path)

        self.check_input()

        self.progress.new('Initializing')
        self.progress.update('Reading the input FASTA')
        self.fasta = u.SequenceSource(self.alignment, lazy_init = False)
        self.progress.end()

        self.column_entropy = [int(x.strip().split()[0]) for x in open(self.entropy).readlines()]
        
        if self.sample_mapping:
            self.sample_mapping_dict = utils.get_sample_mapping_dict(self.sample_mapping)

        self.run.info('project', self.project)
        self.run.info('run_date', utils.get_date())
        self.run.info('version', __version__)
        self.run.info('multi_threaded', not self.no_threading)
        self.run.info('alignment', self.alignment)
        self.run.info('entropy', self.entropy)
        self.run.info('sample_mapping', self.sample_mapping)
        self.run.info('output_directory', self.output_directory)
        self.run.info('tmp_directory', self.tmp_directory)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('quals_provided', True if self.quals_dict else False)
        self.run.info('cmd_line', utils.get_cmd_line(sys.argv))
        self.run.info('total_seq', self.fasta.total_seq)
        self.run.info('alignment_length', self.alignment_length)
        self.run.info('number_of_auto_components', self.number_of_auto_components or 0)
        self.run.info('number_of_selected_components', len(self.selected_components) if self.selected_components else 0)
        self.run.info('generate_sets', self.generate_sets)
        self.run.info('skip_basic_analyses', self.skip_basic_analyses)
        if self.generate_sets:
            self.run.info('T', self.cosine_similarity_threshold)
        self.run.info('s', self.min_number_of_samples)
        self.run.info('a', self.min_percent_abundance)
        self.run.info('A', self.min_actual_abundance)
        self.run.info('M', self.min_substantive_abundance)
        if self.quals_dict:
            self.run.info('q', self.min_base_quality)
        if self.limit_oligotypes_to:
            self.run.info('limit_oligotypes_to', self.limit_oligotypes_to)
        if self.exclude_oligotypes:
            self.run.info('exclude_oligotypes', self.exclude_oligotypes)
        
        if self.number_of_auto_components:
            # locations of interest based on the entropy scores
            self.bases_of_interest_locs = sorted([self.column_entropy[i] for i in range(0, self.number_of_auto_components)])
            self.run.info('bases_of_interest_locs',', '.join([str(x) for x in self.bases_of_interest_locs]))
        elif self.selected_components:
            self.bases_of_interest_locs = sorted(self.selected_components)
            self.run.info('bases_of_interest_locs',', '.join([str(x) for x in self.bases_of_interest_locs]))

        if self.blast_ref_db:
            self.run.info('blast_ref_db', self.blast_ref_db)

        # set number of threads to be used
        if not self.number_of_threads:
            self.number_of_threads = utils.Multiprocessing(None).num_thread
        
        self._construct_samples_dict()
        self._contrive_abundant_oligos()
        self._refine_samples_dict()
        self._get_unit_counts_and_percents()
        self._get_units_across_samples_dicts()
        self._generate_random_colors()
        
        self._generate_FASTA_file()
        self._generate_NEXUS_file()        
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()
        self._store_read_distribution_table()

        if self.generate_sets:
            self._generate_MATRIX_files_for_units_across_samples()
            self._agglomerate_oligos_based_on_cosine_similarity()
            self._generate_MATRIX_files_for_oligotype_sets()       
             
        if ((not self.no_figures) and (not self.quick)) and self.gen_sample_oligo_networks:
            self._generate_sample_oligotype_network_figures()
        if (not self.no_figures) and self.generate_sets:
            self._generate_stack_bar_figure_with_agglomerated_oligos()
            self._generate_oligos_across_samples_figure()
            self._generate_sets_across_samples_figure()
            

        if not self.quick:
            self._generate_representative_sequences()

        if self.representative_sequences_per_oligotype:
            self._generate_representative_sequences_FASTA_file()

        if ((not self.no_figures) and (not self.quick)):
            self._generate_default_figures()

        if ((not self.no_figures) and (not self.quick)) and self.sample_mapping:
            self._generate_exclusive_figures()
            
        if (not self.skip_gexf_network_file) and (not self.quick):
            self._generate_gexf_network_file()

        # store the final information about oligos
        self.run.info('final_oligos', self.abundant_oligos, quiet = True)
        self.run.info('final_oligo_counts_dict', self.final_oligo_counts_dict, quiet = True)
        self.run.info('final_oligo_entropy_distribution_dict', self.final_oligo_entropy_distribution_dict, quiet = True)
        self.run.info('final_oligo_unique_distribution_dict', self.final_oligo_unique_distribution_dict, quiet = True)
        self.run.info('final_purity_score_dict', self.final_purity_score_dict, quiet = True)
        self.run.info('total_purity_score_dict', self.total_purity_score_dict)
        self.run.info('end_of_run', utils.get_date())

        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        self.run.store_info_dict(info_dict_file_path)

        if (not self.keep_tmp):
            shutil.rmtree(self.tmp_directory)

        self.run.quit()

        if self.gen_html:
            self._generate_html_output()

    def _construct_samples_dict(self):
        """This is where oligotypes are being genearted based on bases of each
           alignment at the location of interest"""

        self.progress.new('Sample Dict Construction')

        if self.quals_dict:
            num_reads_eliminated_due_to_min_base_quality = 0

        self.fasta.reset()
        while self.fasta.next():
            if self.fasta.pos % 1000 == 0:
                self.progress.update('Analyzing: %s' \
                                    % (utils.pretty_print(self.fasta.pos)))

            sample = utils.get_sample_name_from_defline(self.fasta.id, self.sample_name_separator)
            
            if not self.samples_dict.has_key(sample):
                self.samples_dict[sample] = {}
                self.samples.append(sample)

            if self.quals_dict:
                # if qual_dicts is available, each base of interest will be tested
                # against --min-base-quality parameter to make sure that it is above
                # the expected quality score. 
                quality_scores = self.quals_dict[self.fasta.id]
                quality_scores_of_bases_of_interest = [quality_scores[o] for o in self.bases_of_interest_locs if not quality_scores[o] == None]
               
                min_base_quality = min([base_quality for base_quality in quality_scores_of_bases_of_interest if base_quality] or [0])

                if min_base_quality < self.min_base_quality:
                    # if True, discard the read
                    # FIXME: Discarded reads should be stored somewhere else for further analysis
                    num_reads_eliminated_due_to_min_base_quality += 1
                    continue
                else:
                    oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
                
            else:
                # if quals_dict is not available, oligotypes will be generated without
                # checking the base qualities
                oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
        
            if self.samples_dict[sample].has_key(oligo):
                self.samples_dict[sample][oligo] += 1
            else:
                self.samples_dict[sample][oligo] = 1
       
        self.samples.sort()
        self.progress.end()
        self.run.info('num_samples_in_fasta', len(self.samples_dict))

        if self.quals_dict:
            self.run.info('num_reads_eliminated_due_to_min_base_quality', num_reads_eliminated_due_to_min_base_quality)
            if self.fasta.total_seq == num_reads_eliminated_due_to_min_base_quality:
                raise utils.ConfigError, "All reads were eliminated due to --min-base-quality (%d) rule" % self.min_base_quality
        

    def _register_removal(self, oligo, reason = 'unknown'):
        if not self.excluded_read_ids_tracker.has_key(reason):
            self.excluded_read_ids_tracker[reason] = {}
            
        for sample in self.samples:
            if self.samples_dict[sample].has_key(oligo):
                if not self.excluded_read_ids_tracker[reason].has_key(sample):
                    self.excluded_read_ids_tracker[reason][sample] = self.samples_dict[sample][oligo]
                else:
                    self.excluded_read_ids_tracker[reason][sample] += self.samples_dict[sample][oligo]

        
    def _contrive_abundant_oligos(self):
        # cat oligos | uniq
        self.progress.new('Contriving Abundant Oligos')

        # a performance optimization workaround in order to lessen the 
        # number of expensive 'keys()' calls on samples_dict to be made
        oligos_in_samples_dict = {}
        for sample in self.samples:
            oligos_in_samples_dict[sample] = set(self.samples_dict[sample].keys())
        
        oligos_set = []
        for sample in self.samples:
            self.progress.update('Unique Oligos: ' + utils.P(self.samples.index(sample), len(self.samples)))
            for oligo in oligos_in_samples_dict[sample]:
                if oligo not in oligos_set:
                    oligos_set.append(oligo)
        self.progress.end()
        self.run.info('num_unique_oligos', len(oligos_set))
       

        self.progress.new('Computing Oligo Abundances')
        # count oligo abundance
        oligo_sample_abundance = []
        for i in range(0, len(oligos_set)):
            oligo = oligos_set[i]
            
            if i % 100 == 0 or i == len(oligos_set) - 1:
                self.progress.update(utils.P(i, len(oligos_set)))
            
            count = 0
            for sample in self.samples:
                if oligo in oligos_in_samples_dict[sample]:
                    count += 1
            oligo_sample_abundance.append((count, oligo),)
        oligo_sample_abundance.sort()
        self.progress.end()

        # eliminate oligos based on the number of samples they appear
        # (any oligo required to appear in at least 'self.min_number_of_samples'
        # samples)
        self.progress.new('Applying -s parameter')
        non_singleton_oligos = []
        for i in range(0, len(oligo_sample_abundance)):
            if i % 100 == 0 or i == len(oligo_sample_abundance) - 1:
                self.progress.update(utils.P(i, len(oligo_sample_abundance)))
            tpl = oligo_sample_abundance[i]
            if tpl[0] >= self.min_number_of_samples:
                non_singleton_oligos.append(tpl[1])
            else:
                self._register_removal(tpl[1], 'failed_s')

        self.progress.end()
        self.run.info('num_oligos_after_s_elim', len(non_singleton_oligos))


        # sample_sums keeps the actual number of oligos that are present in non_singleton_oligos list,
        # for each sample. computing it here once is more optimized.
        sample_sums = {}
        SUM = lambda sample: sum([self.samples_dict[sample][o] for o in non_singleton_oligos \
                                                                if self.samples_dict[sample].has_key(o)])
        for sample in self.samples:
            sample_sums[sample] = SUM(sample)

        # eliminate very rare oligos (the percent abundance of every oligo should be
        # more than 'self.min_percent_abundance' percent in at least one sample)
        self.progress.new('Applying -a parameter')
        for i in range(0, len(non_singleton_oligos)):
            oligo = non_singleton_oligos[i]
            if i % 100 == 0 or i == len(non_singleton_oligos) - 1:
                self.progress.update(utils.P(i, len(non_singleton_oligos)))
            
            percent_abundances = []
            for sample in self.samples:
                if self.samples_dict[sample].has_key(oligo):
                    percent_abundances.append((self.samples_dict[sample][oligo] * 100.0 / sample_sums[sample],
                                               self.samples_dict[sample][oligo],
                                               sample_sums[sample],
                                               sample))

            percent_abundances.sort(reverse = True)

            # NOTE: if a sample has less than 100 sequences, percent abundance doesn't mean much.
            #       if user wants to eliminate oligotypes that doesn't appear in at least one sample
            #       more than 1% abundance, a singleton of that oligotype that appears in a sample
            #       which has 50 sequences would make that oligotype pass the filter. I think if an
            #       oligotype passes the percent filter, sample size and actual count of the oligotype
            #       should also be considered before considering it as an abundant oligotype:
            for abundance_percent, abundance_count, sample_size, sample in percent_abundances:
                PercentAbundance_OK = abundance_percent >= self.min_percent_abundance
                DatesetSize_OK      = sample_size > 100 or abundance_count > self.min_percent_abundance

                if PercentAbundance_OK and DatesetSize_OK:
                    self.abundant_oligos.append((sum([x[1] for x in percent_abundances]), oligo))
                    break
                else:
                    self._register_removal(oligo, 'failed_a')

        self.progress.end()
        self.run.info('num_oligos_after_a_elim', len(self.abundant_oligos))
        
        self.abundant_oligos = [x[1] for x in sorted(self.abundant_oligos, reverse = True)]


        # eliminate very rare oligos (the ACTUAL ABUNDANCE, which is the sum of oligotype in all samples
        # should should be more than 'self.min_actual_abundance'.
        self.progress.new('Applying -A parameter')
        if self.min_actual_abundance > 0:
            oligos_for_removal = []
            for i in range(0, len(self.abundant_oligos)):
                oligo = self.abundant_oligos[i]

                if i % 100 == 0 or i == len(self.abundant_oligos) - 1:
                    self.progress.update(utils.P(i, len(non_singleton_oligos)))

                oligo_actual_abundance = sum([self.samples_dict[sample][oligo] for sample in self.samples_dict\
                                                        if self.samples_dict[sample].has_key(oligo)])
                if self.min_actual_abundance > oligo_actual_abundance:
                    oligos_for_removal.append(oligo)

            for oligo in oligos_for_removal:
                self.abundant_oligos.remove(oligo)
                self._register_removal(oligo, 'failed_A')

        self.progress.end()
        self.run.info('num_oligos_after_A_elim', len(self.abundant_oligos))


        # eliminate oligos based on -M / --min-substantive-abundance parameter.
        #
        # Here is a pesky problem. -A parameter eliminates oligotypes based on the number of sequences
        # represented by them. But this is not a very reliable way to eliminate noise, and sometimes it
        # eliminates more signal than noise. Here is an example: Say Oligotype #1 and Oligotype #2 both
        # represent 20 reads. But O#1 has only one unique sequence, so all reads that are being
        # represented by O#1 are actually the same. Meanwhile O#2 has 20 unique reads in it. So each
        # read differs from each other at bases that are not being used by oligotyping. Simply one could
        # argue that O#2 is full of noise, while O#1 is a robust oligotype that probably represents one
        # and only one organism. If you set -A to 25, both will be eliminated. But if there would be a
        # parameter that eliminates oligotypes based on the number of most abundant unique sequence
        # they entail, it could be set to, say '5', and O#1 would have survived that filter while O#2
        # the crappy oligotype would be filtered out. 
        #
        # Following function, _get_unique_sequence_distributions_within_abundant_oligos, returns the
        # dictionary that can be used to do that.
        #
        # And here is the ugly part about implementing this: This has to be done before the generation
        # of representative sequences. Upto the section where we generate representative sequences,
        # we only work with 'abundances' and we don't actually know what is the distribution of unique
        # sequences an oligotype conceals. This information is being computed when the representative
        # sequences are being computed. But in order to compute representative sequences we need to
        # know 'abundant' oligotypes first, and in order to finalize 'abundant' oligotypes
        # we need to run this cool filter. Chicken/egg. It is extremely inefficient, and I hate
        # to do this but this somewhat redundant step is mandatory and I can't think of any better
        # solution... And if you read this comment all the way here you either must be very bored or
        # very interested in using this codebase properly. Thanks.

        self.progress.new('Applying -M parameter')
        if self.min_substantive_abundance:
            oligos_for_removal = []
            unique_sequence_distributions = self._get_unique_sequence_distributions_within_abundant_oligos()

            num_abundant_oligos = len(self.abundant_oligos)

            for i in range(0, num_abundant_oligos):
                self.progress.update(utils.P(i, num_abundant_oligos))
                oligo = self.abundant_oligos[i]
                if max(unique_sequence_distributions[oligo]) < self.min_substantive_abundance:
                    oligos_for_removal.append(oligo)

            for oligo in oligos_for_removal:
                self._register_removal(oligo, 'failed_M')
                self.abundant_oligos.remove(oligo)

        self.progress.end()
        self.run.info('num_oligos_after_M_elim', len(self.abundant_oligos))


        # if 'limit_oligotypes_to' is defined, eliminate all other oligotypes
        if self.limit_oligotypes_to:
            self.abundant_oligos = [oligo for oligo in self.abundant_oligos if oligo in self.limit_oligotypes_to]
            
            for oligo in [oligo for oligo in self.abundant_oligos if not oligo in self.limit_oligotypes_to]:
                self._register_removal(oligo, 'failed_limit')
            
            self.run.info('num_oligos_after_l_elim', len(self.abundant_oligos))
            if len(self.abundant_oligos) == 0:
                raise utils.ConfigError, "\n\n\t--limit-oligotypes parameter eliminated all oligotypes.\
                                          \n\tPlease make sure --limit-oligotypes matches with actual oligotypes.\n\n\tQuiting.\n"

        # if 'exclude_oligotypes' is defined, remove them from analysis if they are present
        if self.exclude_oligotypes:
            self.abundant_oligos = [oligo for oligo in self.abundant_oligos if not oligo in self.exclude_oligotypes]
            
            for oligo in self.exclude_oligotypes:
                self._register_removal(oligo, 'excluded')
            
            self.run.info('num_oligos_after_e_elim', len(self.abundant_oligos))


        # storing final counts
        for oligo in self.abundant_oligos:
            self.final_oligo_counts_dict[oligo] = sum([self.samples_dict[sample][oligo] for sample in self.samples_dict\
                                                        if self.samples_dict[sample].has_key(oligo)])

        # in case no oligos left
        if not len(self.abundant_oligos):
            raise utils.ConfigError, "\n\n\tAll oligotypes were discarded during the noise removal step.\
                                      \n\tPlease check your parameters.\n\n\tQuiting.\n"

        # if there is only one oligotype left, skip basic analyses
        if len(self.abundant_oligos) == 1:
            self.skip_basic_analyses = True
            self.run.info('skip_basic_analyses', self.skip_basic_analyses)


    def _refine_samples_dict(self):
        # removing oligos from samples dictionary that didn't pass
        # MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE and
        # MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS filters.
        self.progress.new('Refining Samples Dict')

        self.progress.update('Deep-copying the dictionary .. ')
        samples_dict_copy = copy.deepcopy(self.samples_dict)
        self.progress.append('done')

        samples_to_remove = []
        for i in range(0, len(self.samples)):
            sample = self.samples[i]

            self.progress.update('Analyzing samples: ' + utils.P(i + 1, len(self.samples)))
            
            for oligo in samples_dict_copy[sample]:
                if oligo not in self.abundant_oligos:
                    self.samples_dict[sample].pop(oligo)
            if not self.samples_dict[sample]:
                samples_to_remove.append(sample)
        for sample in samples_to_remove:
            self.samples.remove(sample)
            self.samples_dict.pop(sample)

        self.num_sequences_after_qc = sum([sum(self.samples_dict[sample].values()) for sample in self.samples_dict]) 

        self.progress.end()
        self.run.info('num_sequences_after_qc', self.num_sequences_after_qc)

        if len(samples_to_remove):
            self.run.info('samples_removed_after_qc', samples_to_remove)               
            
        if len(self.samples) < 3:
            self.skip_basic_analyses = True
            self.run.info('skip_basic_analyses', self.skip_basic_analyses)
        

    def _generate_FASTA_file(self): 
        # store abundant oligos
        self.progress.new('FASTA File')
        oligos_fasta_file_path = self.generate_output_destination("OLIGOS.fasta")
        f = open(oligos_fasta_file_path, 'w')
        self.progress.update('Being generated')
        for oligo in self.abundant_oligos:
            f.write('>' + oligo + '\n')
            f.write(oligo + '\n')
        f.close()
        self.progress.end()
        self.run.info('oligos_fasta_file_path', oligos_fasta_file_path)
 

    def _generate_representative_sequences_FASTA_file(self): 
        # store representative sequences per oligotype if they are computed
        self.progress.new('Representative Sequences FASTA File')
        representative_seqs_fasta_file_path = self.generate_output_destination("OLIGO-REPRESENTATIVES.fasta")
        f = open(representative_seqs_fasta_file_path, 'w')
        self.progress.update('Being generated')
        for oligo in self.abundant_oligos:
            f.write('>' + oligo + '\n')
            f.write(self.representative_sequences_per_oligotype[oligo] + '\n')
        f.close()
        self.progress.end()
        self.run.info('representative_seqs_fasta_file_path', representative_seqs_fasta_file_path)
        
        
    def _generate_NEXUS_file(self):
        # generate NEXUS file of oligos
        self.progress.new('NEXUS File')
        oligos_nexus_file_path = self.generate_output_destination("OLIGOS.nexus")
        f = open(oligos_nexus_file_path, 'w')
        f.write("""begin data;
            dimensions ntax=%d nchar=%d;
            format datatype=dna interleave=no gap=-;
            matrix\n""" % (len(self.abundant_oligos), len(self.abundant_oligos[0])))
        self.progress.update('Being generated')
        for oligo in self.abundant_oligos:
            f.write('    %.40s %s\n' % (oligo, oligo))
        f.write('    ;\n')
        f.write('end;\n')
        f.close()
        self.progress.end()
        self.run.info('oligos_nexus_file_path', oligos_nexus_file_path)


    def _get_unit_counts_and_percents(self):
        self.progress.new('Unit counts and percents')
        self.progress.update('Data is being generated')
            
        self.unit_counts, self.unit_percents = utils.get_unit_counts_and_percents(self.abundant_oligos, self.samples_dict)
            
        self.progress.end()


    def _generate_MATRIX_files_for_units_across_samples(self):
        self.progress.new('Oligos across samples')
        self.progress.update('Matrix files are being generated')

        across_samples_MN_file_path = self.generate_output_destination("OLIGOS-ACROSS-DATASETS-MAX-NORM.txt")
        across_samples_SN_file_path = self.generate_output_destination("OLIGOS-ACROSS-DATASETS-SUM-NORM.txt")
             
        utils.generate_MATRIX_files_for_units_across_samples(self.abundant_oligos,
                                                             self.samples,
                                                             across_samples_MN_file_path,
                                                             across_samples_SN_file_path,
                                                             self.across_samples_max_normalized,
                                                             self.across_samples_sum_normalized)

        self.progress.end()
        self.run.info('across_samples_MN_file_path', across_samples_MN_file_path)
        self.run.info('across_samples_SN_file_path', across_samples_SN_file_path)


    def _get_units_across_samples_dicts(self):
        self.progress.new('Oligos across samples')
        self.progress.update('Data is being generated')

        self.across_samples_sum_normalized, self.across_samples_max_normalized =\
                utils.get_units_across_samples_dicts(self.abundant_oligos, self.samples, self.unit_percents) 
            
        self.progress.end()

 
    def _generate_ENVIRONMENT_file(self):
        self.progress.new('ENVIRONMENT File')
        self.environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        self.progress.update('Being generated')
        
        utils.generate_ENVIRONMENT_file(self.samples,
                                        self.samples_dict,
                                        self.environment_file_path)

        self.progress.end()
        self.run.info('environment_file_path', self.environment_file_path)

    def _generate_MATRIX_files(self):
        self.progress.new('Matrix Files')
        self.progress.update('Being generated')
            
        self.matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        self.matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")    
            
        utils.generate_MATRIX_files(self.abundant_oligos,
                                   self.samples,
                                   self.unit_counts,
                                   self.unit_percents,
                                   self.matrix_count_file_path,
                                   self.matrix_percent_file_path)
            
        self.progress.end()
        self.run.info('matrix_count_file_path', self.matrix_count_file_path)
        self.run.info('matrix_percent_file_path', self.matrix_percent_file_path)


    def _store_read_distribution_table(self):
        self.progress.new('Read distribution table')
        self.read_distribution_table_path = self.generate_output_destination("READ-DISTRIBUTION.txt")

        def get_dict_entry_tmpl():
            d = {'represented_reads': 0}
            for reason in self.excluded_read_ids_tracker:
                d[reason] = 0
            return d

        read_distribution_dict = {}
        
        self.progress.update('Processing reads that were represented in results')
        for sample in self.samples_dict:
            if not read_distribution_dict.has_key(sample):
                read_distribution_dict[sample] = get_dict_entry_tmpl()

            read_distribution_dict[sample]['represented_reads'] = sum(self.samples_dict[sample].values())

        for reason in self.excluded_read_ids_tracker:
            self.progress.update('Processing excluded oligos (%s)' % (reason))
            for sample in self.excluded_read_ids_tracker[reason]:
                if not read_distribution_dict.has_key(sample):
                    read_distribution_dict[sample] = get_dict_entry_tmpl()
                        
                read_distribution_dict[sample][reason] = self.excluded_read_ids_tracker[reason][sample]

    
        self.progress.update('Storing...')
        utils.generate_TAB_delim_file_from_dict(read_distribution_dict,
                                                self.read_distribution_table_path,
                                                order = ['represented_reads'] + sorted(self.excluded_read_ids_tracker.keys()))

        self.progress.end()
        self.run.info('read_distribution_table_path', self.read_distribution_table_path)


    def _generate_random_colors(self):
        self.colors_file_path = self.generate_output_destination('COLORS')
        if self.colors_list_file:
            # it means user provided a list of colors to be used for oligotypes
            colors = [c.strip() for c in open(self.colors_list_file).readlines()]
            if len(colors) < len(self.abundant_oligos):
                raise utils.ConfigError, "Number of colors defined in colors file (%d),\
                                          is smaller than the number of abundant oligotypes (%d)" % \
                                                        (len(colors), len(self.abundant_oligos))
            colors_dict = {}
            for i in range(0, len(self.abundant_oligos)):
                colors_dict[self.abundant_oligos[i]] = colors[i]

            self.colors_dict = colors_dict
            
            # generate COLORS file derived from --colors-list-file
            colors_file = open(self.colors_file_path, 'w')
            for oligotype in self.abundant_oligos:
                colors_file.write('%s\t%s\n' % (oligotype, self.colors_dict[oligotype]))
            colors_file.close()

        else:
            self.colors_dict = random_colors(self.abundant_oligos, self.colors_file_path)
        self.run.info('colors_file_path', self.colors_file_path)

    
    def _agglomerate_oligos_based_on_cosine_similarity(self):
        from Oligotyping.utils.cosine_similarity import get_oligotype_sets
        self.progress.new('Agglomerating Oligotypes into Sets')
        oligotype_sets_file_path = self.generate_output_destination("OLIGOTYPE-SETS.txt")
        self.progress.update('Computing')
        self.oligotype_sets = get_oligotype_sets(self.abundant_oligos,
                                                 self.across_samples_sum_normalized,
                                                 self.cosine_similarity_threshold,
                                                 oligotype_sets_file_path)
        
        self.progress.end()
        self.run.info('oligotype_sets_file_path', oligotype_sets_file_path)
        self.run.info('oligotype_sets_info', '%d oligotypes agglomerated into %d sets'\
                                            % (len(self.abundant_oligos), len(self.oligotype_sets)))


        self.progress.new('Generating data objects for newly generated oligotype sets')
        self.progress.update('New Colors')
        self.oligotype_set_ids = range(0, len(self.oligotype_sets))
        
        self.colors_dict_for_oligotype_sets = {}
        for set_id in self.oligotype_set_ids:
            self.colors_dict_for_oligotype_sets[set_id] = self.colors_dict[self.oligotype_sets[set_id][0]]

        self.progress.update('New Samples Dict')
        self.samples_dict_with_agglomerated_oligos = {}
        for sample in self.samples:
            self.samples_dict_with_agglomerated_oligos[sample] = {}

        for set_id in self.oligotype_set_ids:
            oligotype_set = self.oligotype_sets[set_id]
            for sample in self.samples:
                self.samples_dict_with_agglomerated_oligos[sample][set_id] = 0
                for oligo in self.samples_dict[sample]:
                    if oligo in oligotype_set:
                        self.samples_dict_with_agglomerated_oligos[sample][set_id] += self.samples_dict[sample][oligo]

        self.progress.end()


    def _generate_MATRIX_files_for_oligotype_sets(self):
        self.progress.new('Matrix Files for Oligotype Sets')
        counts_file_path = self.generate_output_destination("MATRIX-COUNT-OLIGO-SETS.txt")
        percents_file_path = self.generate_output_destination("MATRIX-PERCENT-OLIGO-SETS.txt")
        
        d = self.samples_dict_with_agglomerated_oligos
        oligotype_set_percents = {}
        oligotype_set_counts = {}


        self.progress.update('Generating the data')
        for oligotype_set_id in self.oligotype_set_ids:
            counts = []
            percents = []
            for sample in self.samples:
                if d[sample].has_key(oligotype_set_id):
                    counts.append(d[sample][oligotype_set_id])
                    percents.append(d[sample][oligotype_set_id] * 100.0 / sum(d[sample].values()))
                else:
                    counts.append(0)
                    percents.append(0.0)

            oligotype_set_percents[oligotype_set_id] = percents
            oligotype_set_counts[oligotype_set_id] = counts
        
        self.progress.update('Generating files')
        counts_file = open(counts_file_path, 'w')
        percents_file = open(percents_file_path, 'w')       
        
        counts_file.write('\t'.join([''] + self.samples) + '\n')
        percents_file.write('\t'.join([''] + self.samples) + '\n')

        for oligotype_set_id in self.oligotype_set_ids:
            counts_file.write('\t'.join(['Set_' + str(oligotype_set_id)] + [str(c) for c in oligotype_set_counts[oligotype_set_id]]) + '\n')
            percents_file.write('\t'.join(['Set_' + str(oligotype_set_id)] + [str(p) for p in oligotype_set_percents[oligotype_set_id]]) + '\n')
        
        counts_file.close()
        percents_file.close()

        self.progress.end()
        self.run.info('matrix_count_oligo_sets_file_path', counts_file_path)
        self.run.info('matrix_percent_oligo_sets_file_path', percents_file_path)


    def _get_unique_sequence_distributions_within_abundant_oligos(self):
        # compute and return the unique sequence distribution within per oligo
        # dictionary. see the explanation where the function is called. oligos
        # listed in this dictionary MAY NOT be the final oligos once the noise
        # filtering step has ended.

        temp_unique_distributions = dict(zip(self.abundant_oligos, [{} for x in range(0, len(self.abundant_oligos))]))

        self.fasta.reset()
        while self.fasta.next():
            if self.progress and self.fasta.pos % 1000 == 0:
                self.progress.update('Computing sequence distributions: %.2f%%' \
                                                % (self.fasta.pos * 100.0 / self.fasta.total_seq))
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
            if oligo in self.abundant_oligos:
                try:
                    temp_unique_distributions[oligo][self.fasta.seq] += 1
                except KeyError:
                    temp_unique_distributions[oligo][self.fasta.seq] = 1

        for oligo in self.abundant_oligos:
            temp_unique_distributions[oligo] = sorted(temp_unique_distributions[oligo].values(), reverse = True)


        return temp_unique_distributions


    def _generate_representative_sequences(self):
        # create a fasta file with a representative full length consensus sequence for every oligotype

        # this is what is going on here: we go through all oligotypes, gather sequences that are being
        # represented by a particular oligotype, unique them and report the top ten unique sequences
        # ordered by the frequency.
        self.progress.new('Representative Sequences')

        output_directory_for_reps = self.generate_output_destination("OLIGO-REPRESENTATIVES", directory = True)

        fasta_files_dict = {}
        unique_files_dict = {}
        for oligo in self.abundant_oligos:
            if oligo not in fasta_files_dict:
                try:
                    fasta_file_path = os.path.join(output_directory_for_reps, '%.5d_' % self.abundant_oligos.index(oligo) + oligo)
                    fasta_files_dict[oligo] = {'file': open(fasta_file_path, 'w'),
                                               'path': fasta_file_path}
                    unique_files_dict[oligo] = {'file': open(fasta_file_path + '_unique', 'w'),
                                                'path': fasta_file_path + '_unique'}
                except IOError:
                    print '\n\t'.join(['',
                                       'WARNING: Oligotyping process has reached the maximum number of open files',
                                       'limit defined by the operating system. There are "%d" oligotypes to be'\
                                                                 % len(self.abundant_oligos),
                                       'stored. You can learn the actual limit by typing "ulimit -n" in the.run.',
                                       '',
                                       'You can increase this limit temporarily by typing "ulimit -n NUMBER", and',
                                       'restart the process. It seems using %d as NUMBER might be a good start.'\
                                                                % (len(self.abundant_oligos) * 1.1),
                                       '',
                                       'Until this issue is solved, representative sequences are not going to be',
                                       'computed.',
                                       ''])

                    # clean after yourself. close every file, delete directory, exit.
                    [map(lambda x: x.close(), [g[o]['file'] for o in g]) for g in [fasta_files_dict, unique_files_dict]]
                    shutil.rmtree(output_directory_for_reps)
                    sys.exit()

        self.fasta.reset()
        while self.fasta.next():
            if self.fasta.pos % 1000 == 0:
                self.progress.update('Generating Individual FASTA Files: %.2f%%' \
                                                % (self.fasta.pos * 100.0 / self.fasta.total_seq))
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
            if oligo in self.abundant_oligos:
                fasta_files_dict[oligo]['file'].write('>%s\n' % (self.fasta.id))
                fasta_files_dict[oligo]['file'].write('%s\n' % self.fasta.seq)
        
        self.progress.end()

        self.progress.new('Representative Sequences')
        for oligo in self.abundant_oligos:
            fasta_files_dict[oligo]['file'].close()


            self.progress.update('Unique reads for %s (%d of %d)' \
                                        % (oligo,
                                           self.abundant_oligos.index(oligo) + 1,
                                           len(self.abundant_oligos)))
            fasta_file_path = fasta_files_dict[oligo]['path']
            fasta = u.SequenceSource(fasta_file_path, lazy_init = False, unique = True)
          
            # this dict is going to hold the information of how unique sequences within an oligotype
            # is distributed among samples:
            distribution_among_samples = {}

            fasta.next()
            # this is the first read in the unique reads list, which is the most abundant unique sequence
            # for the oligotype. so we are going to store it in a dict to generate
            # representative sequences FASTA file:
            self.representative_sequences_per_oligotype[oligo] = fasta.seq
            fasta.reset()


            # FIXME: I am going to come back to this and fix it at some point. Storing 'distribution_among_samples'
            # information in separate cPickle files per oligo is not the smartest thing to do.
            self.final_oligo_unique_distribution_dict[oligo] = []
            while fasta.next() and fasta.pos <= self.limit_representative_sequences:
                unique_files_dict[oligo]['file'].write('>%s_%d|freq:%d\n'\
                                                                     % (oligo,
                                                                        fasta.pos,
                                                                        len(fasta.ids)))
                unique_files_dict[oligo]['file'].write('%s\n' % fasta.seq)
                
                # store only the first 20
                if not fasta.pos > 20:
                    self.final_oligo_unique_distribution_dict[oligo].append(len(fasta.ids))

                for sample_id in fasta.ids:
                    sample_name = utils.get_sample_name_from_defline(sample_id, self.sample_name_separator)
                    if not distribution_among_samples.has_key(sample_name):
                        distribution_among_samples[sample_name] = {}
                    d = distribution_among_samples[sample_name]
                    if not d.has_key(fasta.pos):
                        d[fasta.pos] = 1
                    else:
                        d[fasta.pos] += 1
                
            fasta.close()
            unique_files_dict[oligo]['file'].close()

            unique_fasta_path = unique_files_dict[oligo]['path']
            distribution_among_samples_dict_path = unique_fasta_path + '_distribution.cPickle'
            cPickle.dump(distribution_among_samples, open(distribution_among_samples_dict_path, 'w'))

        self.progress.end()
        

        self._get_purity_score()
        self._get_total_purity_score()

        
        self.progress.new('Generating Entropy Figures')
        if (not self.quick) and (not self.no_figures):
            if self.no_threading:
                for oligo in self.abundant_oligos:
                    self.progress.update('%s (%d of %d)' % (oligo,
                                                            self.abundant_oligos.index(oligo) + 1,
                                                            len(self.abundant_oligos)))
                    unique_fasta_path = unique_files_dict[oligo]['path']
                    self._generate_entropy_figure_for_abundant_oligotype(oligo, unique_fasta_path, self.final_oligo_entropy_distribution_dict)
            else:
                mp = utils.Multiprocessing(self._generate_entropy_figure_for_abundant_oligotype, self.number_of_threads)
                entropy_per_oligo_shared_dict = mp.get_empty_shared_dict()

                # arrange processes
                processes_to_run = []
                for oligo in self.abundant_oligos:
                    unique_fasta_path = unique_files_dict[oligo]['path']
                    processes_to_run.append((oligo, unique_fasta_path, entropy_per_oligo_shared_dict),)
    
                # start the main loop to run all processes
                mp.run_processes(processes_to_run, self.progress)
                
                self.final_oligo_entropy_distribution_dict = copy.deepcopy(entropy_per_oligo_shared_dict)
        self.progress.end()


        if (not self.quick) and (not self.skip_blast_search):
            self.progress.new('Performing %s BLAST search for representative sequences'\
                                            % ('LOCAL' if self.blast_ref_db else 'REMOTE'))

            if self.blast_ref_db:
                # if there is a local db to search representative sequences against,
                # just perform the blast search in one thread
                self._perform_local_BLAST_search_for_oligo_representative(unique_files_dict)
            else:
                # if the search is going to be on NCBI, parallelize it:
                if self.no_threading:
                    for oligo in self.abundant_oligos:
                        self.progress.update('%s (%d of %d)' % (oligo,
                                                                self.abundant_oligos.index(oligo) + 1,
                                                                len(self.abundant_oligos)))
                        self._perform_remote_BLAST_search_for_oligo_representative(oligo, unique_files_dict)
                else:
                    mp = utils.Multiprocessing(self._perform_remote_BLAST_search_for_oligo_representative, self.number_of_threads)
        
                    # arrange processes
                    processes_to_run = []
                    for oligo in self.abundant_oligos:
                        unique_fasta_path = unique_files_dict[oligo]['path']
                        processes_to_run.append((oligo, unique_files_dict,),)
        
                    # start the main loop to run all processes
                    mp.run_processes(processes_to_run, self.progress)
            
            self.progress.end()

        self.run.info('output_directory_for_reps', output_directory_for_reps) 

    def _get_purity_score(self):
        for oligo in self.final_oligo_unique_distribution_dict:
            freq_dict = self.final_oligo_unique_distribution_dict[oligo]
            if len(self.final_oligo_unique_distribution_dict[oligo]) > 1:
                bp = (freq_dict[1] / (freq_dict[0] * 1.0))
                self.final_purity_score_dict[oligo] = 1 - bp
            else:
                self.final_purity_score_dict[oligo] = 1.00


    def _get_total_purity_score(self):
        sorted_scores = sorted(self.final_purity_score_dict.values())
        last_quarter = sorted_scores[:int(math.ceil(len(sorted_scores)/4.0))] # take the last quarter of the unique sequences   
        final_total = reduce(lambda x, y: x + y, last_quarter) / len(last_quarter) # take the average of these sequences 
        
        self.total_purity_score_dict =  "%.2f" %final_total

    def _perform_local_BLAST_search_for_oligo_representative(self, unique_files_dict):            
        query, target, output = utils.get_temporary_file_names_for_BLAST_search(prefix = "REPS_", directory = self.tmp_directory)
                
        representative_fasta_entries = []
        for oligo in self.abundant_oligos:
            self.progress.update('Storing representative sequences for "%s" ...' % oligo)
            unique_fasta_path = unique_files_dict[oligo]['path']
            unique_fasta = u.SequenceSource(unique_fasta_path)
            unique_fasta.next()
            representative_fasta_entries.append((oligo, unique_fasta.seq),)
            unique_fasta.close()
        utils.append_reads_to_FASTA(representative_fasta_entries, query)


        self.progress.update('Generating a copy of target BLAST db ...')
        self.logger.info('copying blastdb from "%s" to %s' % (self.blast_ref_db, target))
        shutil.copy(self.blast_ref_db, target)
        utils.mask_defline_whitespaces_in_FASTA(target, '<$!$>')
        utils.mask_defline_whitespaces_in_FASTA(query, '<$!$>')
    
        params = "-perc_identity 90"
        job = 'reps'
        s = blast.LocalBLAST(query, target, output, log = self.generate_output_destination('BLAST.log'))
        self.logger.info('local blast request for job "%s": (q) %s (t) %s (o) %s (p) %s'\
                                               % (job, query, target, output, params))
        
        
        self.progress.update('Running makeblastdb ...')
        s.make_blast_db()
        self.logger.info('makeblastdb for %s: %s' % (job, s.makeblastdb_cmd))
        
        s.params = params
        self.progress.update('Performing blastn ...')
        s.search()
        self.logger.info('blastn for %s: %s' % (job, s.search_cmd))
        
        self.progress.update('Processing BLAST results ...')
        fancy_results_dict = s.get_fancy_results_dict(defline_white_space_mask = '<$!$>')

        self.progress.update('Storing BLAST results ...')
        for oligo in self.abundant_oligos:
            unique_fasta_path = unique_files_dict[oligo]['path']
            fancy_blast_result_output_path = unique_fasta_path + '_BLAST.cPickle'
            if fancy_results_dict.has_key(oligo):
                cPickle.dump(fancy_results_dict[oligo], open(fancy_blast_result_output_path, 'w'))
            else:
                cPickle.dump([], open(fancy_blast_result_output_path, 'w'))


    def _perform_remote_BLAST_search_for_oligo_representative(self, oligo, unique_files_dict):
        # will perform remote BLAST
        r = blast.RemoteBLAST()
        
        unique_fasta_path = unique_files_dict[oligo]['path']
        unique_fasta = u.SequenceSource(unique_fasta_path)
        unique_fasta.next()
        blast_output_xml = unique_fasta_path + '_BLAST.xml'
        blast_output_dict = unique_fasta_path + '_BLAST.cPickle'

        # FIXME: this value should be paramaterized
        max_blast_attempt = 3

        def blast_search_wrapper(seq, xml_path, pickle_path):
            try:
                results = r.search(seq, xml_path)
                results_list = r.get_fancy_results_list(results)
                cPickle.dump(results_list, open(pickle_path, 'w'))
                return True
            except:
                return False

        for blast_attempt in range(0, max_blast_attempt):
            self.progress.update('searching for "%s" (%d of %d) (attempt #%d)' % (oligo,
                                                                                  self.abundant_oligos.index(oligo) + 1,
                                                                                  len(self.abundant_oligos), blast_attempt + 1))
                            
            if blast_search_wrapper(unique_fasta.seq, blast_output_xml, blast_output_dict):
                break
            else:
                continue

        unique_fasta.close()

        return True

    def _generate_entropy_figure_for_abundant_oligotype(self, oligo, unique_fasta_path, final_oligo_entropy_distribution_dict):
        entropy_file_path = unique_fasta_path + '_entropy'
        color_per_column_path  = unique_fasta_path + '_color_per_column.cPickle'

        # generate entropy output at 'entropy_file_path' along with the image
        vis_freq_curve(unique_fasta_path, output_file = unique_fasta_path + '.png', entropy_output_file = entropy_file_path)

        # use entropy output to generate a color shade for every columns in alignment
        # for visualization purposes
        entropy_values_per_column = [0] * self.alignment_length
        for column, entropy in [x.strip().split('\t') for x in open(entropy_file_path)]:
            entropy_values_per_column[int(column)] = float(entropy)
        
        final_oligo_entropy_distribution_dict[oligo] = entropy_values_per_column
            
        color_shade_dict = get_color_shade_dict_for_list_of_values(entropy_values_per_column)

        color_per_column = [0] * self.alignment_length
        for i in range(0, self.alignment_length):
            color_per_column[i] = color_shade_dict[entropy_values_per_column[i]]        

        cPickle.dump(color_per_column, open(color_per_column_path, 'w'))
    
    def _generate_sample_oligotype_network_figures(self):
        output_directory_for_samples = self.generate_output_destination("DATASETS", directory = True)
        oligotype_network_structure(self.run.info_dict['environment_file_path'], output_dir = output_directory_for_samples)
        self.run.info('output_directory_for_samples', output_directory_for_samples) 
 

    def _generate_oligos_across_samples_figure(self):
        self.progress.new('Oligotypes Across Samples Figure')
        oligos_across_samples_file_path = self.generate_output_destination('OLIGOS-ACROSS-DATASETS.png')
        self.progress.update('Generating')
        oligos = copy.deepcopy(self.abundant_oligos)
        oligotype_distribution_across_samples(self.samples_dict, self.colors_dict, oligos_across_samples_file_path,\
                                               oligos = oligos, project_title = self.project, display = False)
        self.progress.end()
        self.run.info('oligos_across_samples_file_path', oligos_across_samples_file_path)


    def _generate_sets_across_samples_figure(self):
        self.progress.new('Oligotype Sets Across Samples Figure')
        figure_path = self.generate_output_destination('OLIGO-SETS-ACROSS-DATASETS.png')
        self.progress.update('Generating')
        vis_oligotype_sets_distribution(self.oligotype_sets, self.across_samples_sum_normalized, self.samples,\
                               display = False, colors_dict = self.colors_dict, output_file = figure_path,\
                               project_title = 'Oligotype Sets Across Samples for "%s", at Cosine Similarity Threshold of %.4f'\
                                        % (self.project, self.cosine_similarity_threshold), legend = False)
        self.progress.end()
        self.run.info('oligotype_sets_across_samples_figure_path', figure_path)


    def _generate_stack_bar_figure_with_agglomerated_oligos(self):
        self.progress.new('Stackbar Figure with Agglomerated Oligos')
        stack_bar_file_path = self.generate_output_destination('STACKBAR-AGGLOMERATED-OLIGOS.png')
        self.progress.update('Generating')

        oligotype_distribution_stack_bar(self.samples_dict_with_agglomerated_oligos, self.colors_dict_for_oligotype_sets,\
                                         stack_bar_file_path, oligos = self.oligotype_set_ids, project_title = self.project,\
                                         display = not self.no_display)
        self.progress.end()
        self.run.info('stack_bar_with_agglomerated_oligos_file_path', stack_bar_file_path)


    def _generate_default_figures(self):

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


    def _generate_gexf_network_file(self):
        self.gexf_network_file_path = self.generate_output_destination("NETWORK.gexf")

        self.progress.new('GEXF Network File')
       
        utils.generate_gexf_network_file(self.abundant_oligos,
                                         self.samples_dict,
                                         self.unit_percents,
                                         self.gexf_network_file_path,
                                         sample_mapping_dict = self.sample_mapping_dict,
                                         project = self.project)
        
        self.progress.end()
        self.run.info('gexf_network_file_path', self.gexf_network_file_path)


    def _generate_html_output(self):
        if self.no_figures:
            sys.stdout.write('\n\n\t"--no-figures" parameter is given, skipping HTML output...\n\n')
            return
        
        if self.quick:
            sys.stdout.write('\n\n\t"--quick" parameter is given, skipping HTML output...\n\n')
            return

        from Oligotyping.utils.html.error import HTMLError
        try:
            from Oligotyping.utils.html.for_oligotyping import generate_html_output
        except HTMLError, e:
            sys.stdout.write('\n\n\t%s\n\n' % e)
            sys.exit()

        self.progress.new('HTML Output')
        output_directory_for_html = self.generate_output_destination("HTML-OUTPUT", directory = True)
        self.progress.update('Generating')
        index_page = generate_html_output(self.run.info_dict, html_output_directory = output_directory_for_html)
        self.progress.end()
        sys.stdout.write('\n\n\tView results in your browser: "%s"\n\n' % index_page)


if __name__ == '__main__':
    pass
