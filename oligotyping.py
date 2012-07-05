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

import os
import sys
import copy
import shutil
import cPickle
import tempfile
import operator


from lib import fastalib as u
from visualization.frequency_curve_and_entropy import vis_freq_curve
from visualization.oligotype_distribution_stack_bar import oligotype_distribution_stack_bar
from visualization.oligotype_network_structure import oligotype_network_structure
from utils.random_colors import random_colors
from utils.random_colors import get_color_shade_dict_for_list_of_values
from utils.utils import P
from utils.utils import Run
from utils.utils import Progress
from utils.utils import pretty_print
from utils.utils import get_terminal_size
from utils.utils import process_command_line_args_for_quality_files

# FIXME: test whether Biopython is installed or not here.
from utils.blast_interface import remote_blast_search, local_blast_search


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Config Error: %s' % self.e


class Oligotyping:
    def __init__(self, args = None):
        self.entropy   = None
        self.alignment = None
        self.quals_dict = None
        self.min_base_quality = None
        self.number_of_auto_components = 5
        self.selected_components = None
        self.limit_oligotypes_to = None
        self.min_number_of_datasets = 5
        self.min_percent_abundance = 0.0
        self.project = None
        self.output_directory = None
        self.dataset_name_separator = '_'
        self.limit_representative_sequences = sys.maxint
        self.quick = False
        self.no_figures = False
        self.no_display = False
        self.blast_ref_db = None
        self.skip_blast_search = False
        self.gen_html = False
        self.gen_dataset_oligo_networks = False

        Absolute = lambda x: os.path.join(os.getcwd(), x) if not x.startswith('/') else x 

        if args:
            self.entropy = Absolute(args.entropy)
            self.alignment = Absolute(args.alignment)
            self.quals_dict = process_command_line_args_for_quality_files(args, _return = 'quals_dict')
            self.min_base_quality = args.min_base_quality
            self.number_of_auto_components = args.number_of_auto_components
            self.selected_components = args.selected_components
            self.limit_oligotypes_to = args.limit_oligotypes_to
            self.min_number_of_datasets = args.min_number_of_datasets
            self.min_percent_abundance = args.min_percent_abundance
            self.min_actual_abundance = args.min_actual_abundance
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.output_directory = args.output_directory
            self.dataset_name_separator = args.dataset_name_separator
            self.limit_representative_sequences = args.limit_representative_sequences or sys.maxint
            self.quick = args.quick
            self.no_figures = args.no_figures
            self.no_display = args.no_display
            self.blast_ref_db = Absolute(args.blast_ref_db) if args.blast_ref_db else None
            self.skip_blast_search = args.skip_blast_search
            self.gen_html = args.gen_html
            self.gen_dataset_oligo_networks = args.gen_dataset_oligo_networks
        
        self.run = None
        self.progress = Progress()

        self.datasets_dict = {}
        self.datasets = []
        self.abundant_oligos = []
        self.colors_dict = None

    def sanity_check(self):
        if self.number_of_auto_components != None and self.selected_components != None:
            raise ConfigError, "Both 'auto components' (-c) and 'selected components' (-C) has been declared."
        
        if self.number_of_auto_components == None and self.selected_components == None:
            raise ConfigError, "Either only 'auto components' (-c), or only 'selected components' (-C) can be declared."

        if self.selected_components:
            try:
                self.selected_components = [int(c) for c in self.selected_components.split(',')]
            except:
                raise ConfigError, "Selected components should be comma separated integer values (such as '4,8,15,25,47')."

        if self.min_base_quality:
            try:
                self.min_base_quality = int(self.min_base_quality)
                assert(self.min_base_quality >= 0 and self.min_base_quality <= 40)
            except:
                raise ConfigError, "Minimum base quality must be an integer between 0 and 40."

        if self.limit_oligotypes_to:
            self.limit_oligotypes_to = [o.strip().upper() for o in self.limit_oligotypes_to.split(',')]
            if len(self.limit_oligotypes_to) == 1:
                raise ConfigError, "There must be more than one oligotype for --limit-oligotypes parameter."

            if len([n for n in ''.join(self.limit_oligotypes_to) if n not in ['A', 'T', 'C', 'G', '-']]):
                raise ConfigError, "Oligotypes defined by --limit-oligotypes parameter seems to have ambiguours characters."
        
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

        if (not os.path.exists(self.alignment)) or (not os.access(self.alignment, os.R_OK)):
            raise ConfigError, "Alignment file is not accessible: '%s'" % self.alignment
        
        if (not os.path.exists(self.entropy)) or (not os.access(self.entropy, os.R_OK)):
            raise ConfigError, "Entropy file is not accessible: '%s'" % self.entropy

        return True


    def dataset_name_from_defline(self, defline):
        return self.dataset_name_separator.join(defline.split('|')[0].split(self.dataset_name_separator)[0:-1])


    def get_prefix(self):
        prefix = 's%d-a%.1f-A%d' % (self.min_number_of_datasets, self.min_percent_abundance, self.min_actual_abundance)

        if self.selected_components:
            prefix = 'sc%d-%s' % (len(self.selected_components), prefix)
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
        self.sanity_check()
        
        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.run = Run(self.info_file_path)

        self.fasta = u.SequenceSource(self.alignment, lazy_init = False)
        self.column_entropy = [int(x.strip().split()[0]) for x in open(self.entropy).readlines()]

        self.fasta.next()
        self.alignment_length = len(self.fasta.seq)
        self.fasta.reset()

        self.run.info('project', self.project)
        self.run.info('alignment', self.alignment)
        self.run.info('entropy', self.entropy)
        self.run.info('output_directory', self.output_directory)
        self.run.info('info_file_path', self.info_file_path)
        self.run.info('quals_provided', True if self.quals_dict else False)
        self.run.info('cmd_line', ' '.join(sys.argv).replace(', ', ','))
        self.run.info('total_seq', pretty_print(self.fasta.total_seq))
        self.run.info('alignment_length', pretty_print(self.alignment_length))
        self.run.info('number_of_auto_components', self.number_of_auto_components or 0)
        self.run.info('number_of_selected_components', len(self.selected_components) if self.selected_components else 0)
        self.run.info('s', self.min_number_of_datasets)
        self.run.info('a', self.min_percent_abundance)
        self.run.info('A', self.min_actual_abundance)
        if self.quals_dict:
            self.run.info('q', self.min_base_quality)
        if self.limit_oligotypes_to:
            self.run.info('limit_oligotypes_to', self.limit_oligotypes_to)
        
        if self.number_of_auto_components:
            # locations of interest based on the entropy scores
            self.bases_of_interest_locs = sorted([self.column_entropy[i] for i in range(0, self.number_of_auto_components)])
            self.run.info('bases_of_interest_locs',', '.join([str(x) for x in self.bases_of_interest_locs]))
        elif self.selected_components:
            self.bases_of_interest_locs = sorted(self.selected_components)
            self.run.info('bases_of_interest_locs',', '.join([str(x) for x in self.bases_of_interest_locs]))

        if self.blast_ref_db:
            self.run.info('blast_ref_db', self.blast_ref_db)

        self._construct_datasets_dict()
        self._contrive_abundant_oligos()
        self._refine_datasets_dict()
        self._generate_FASTA_file()
        self._generate_NEXUS_file()
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()
        self._generate_viamics_datasets_dict()
        if not self.quick:
            self._generate_representative_sequences()
        self._generate_random_colors()
        if ((not self.no_figures) and (not self.quick)) and self.gen_dataset_oligo_networks:
            self._generate_dataset_oligotype_network_figures()
        if not self.no_figures:
            self._generate_stack_bar_figure()

        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        self.run.store_info_dict(info_dict_file_path)
        self.run.quit()

        if self.gen_html:
            self._generate_html_output()


    def _construct_datasets_dict(self):
        """This is where oligotypes are being genearted based on bases of each
           alignment at the location of interest"""

        self.progress.new('Dataset Dict Construction')

        if self.quals_dict:
            num_reads_eliminated_due_to_min_base_quality = 0

        self.fasta.reset()
        while self.fasta.next():
            if self.fasta.pos % 1000 == 0:
                self.progress.update('Analyzing: %s' \
                                    % (pretty_print(self.fasta.pos)))

            dataset = self.dataset_name_from_defline(self.fasta.id)
            
            if not self.datasets_dict.has_key(dataset):
                self.datasets_dict[dataset] = {}
                self.datasets.append(dataset)

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
        
            if self.datasets_dict[dataset].has_key(oligo):
                self.datasets_dict[dataset][oligo] += 1
            else:
                self.datasets_dict[dataset][oligo] = 1
       
        self.progress.reset()
        self.run.info('num_datasets_in_fasta', pretty_print(len(self.datasets_dict)))

        if self.quals_dict:
            self.run.info('num_reads_eliminated_due_to_min_base_quality', pretty_print(num_reads_eliminated_due_to_min_base_quality))
            if self.fasta.total_seq == num_reads_eliminated_due_to_min_base_quality:
                raise ConfigError, "All reads were eliminated due to --min-base-quality (%d) rule" % self.min_base_quality
        

    def _contrive_abundant_oligos(self):
        # cat oligos | uniq
        self.progress.new('Contriving Abundant Oligos')

        # a performance optimization workaround in order to lessen the 
        # number of expensive 'keys()' calls on datasets_dict to be made
        oligos_in_datasets_dict = {}
        for dataset in self.datasets:
            oligos_in_datasets_dict[dataset] = set(self.datasets_dict[dataset].keys())
        
        oligos_set = []
        for dataset in self.datasets:
            self.progress.update('Unique Oligos: ' + P(self.datasets.index(dataset), len(self.datasets)))
            for oligo in oligos_in_datasets_dict[dataset]:
                if oligo not in oligos_set:
                    oligos_set.append(oligo)
        self.progress.reset()
        self.run.info('num_unique_oligos', pretty_print(len(oligos_set)))
       

        # count oligo abundance
        oligo_dataset_abundance = []
        for i in range(0, len(oligos_set)):
            oligo = oligos_set[i]
            
            if i % 100 == 0 or i == len(oligos_set) - 1:
                self.progress.update('Counting oligo abundance: ' + P(i, len(oligos_set)))
            
            count = 0
            for dataset in self.datasets:
                if oligo in oligos_in_datasets_dict[dataset]:
                    count += 1
            oligo_dataset_abundance.append((count, oligo),)
        oligo_dataset_abundance.sort()
        self.progress.reset()

        # eliminate oligos based on the number of datasets they appear
        # (any oligo required to appear in at least 'self.min_number_of_datasets'
        # datasets)
        non_singleton_oligos = []
        for i in range(0, len(oligo_dataset_abundance)):
            if i % 100 == 0 or i == len(oligo_dataset_abundance) - 1:
                 self.progress.update('Eliminating singletons: ' + P(i, len(oligo_dataset_abundance)))
            tpl = oligo_dataset_abundance[i]
            if tpl[0] >= self.min_number_of_datasets:
                non_singleton_oligos.append(tpl[1])
        self.progress.reset()
        self.run.info('num_oligos_after_s_elim', pretty_print(len(non_singleton_oligos)))
        

        # dataset_sums keeps the actual number of oligos that are present in non_singleton_oligos list,
        # for each dataset. computing it here once is more optimized.
        dataset_sums = {}
        SUM = lambda dataset: sum([self.datasets_dict[dataset][o] for o in non_singleton_oligos \
                                                                if self.datasets_dict[dataset].has_key(o)])
        for dataset in self.datasets:
            dataset_sums[dataset] = SUM(dataset)

        # eliminate very rare oligos (the percent abundance of every oligo should be
        # more than 'self.min_percent_abundance' percent in at least one dataset)
        for i in range(0, len(non_singleton_oligos)):
            oligo = non_singleton_oligos[i]
            if i % 100 == 0 or i == len(non_singleton_oligos) - 1:
                 self.progress.update('Applying -a parameter: ' + P(i, len(non_singleton_oligos)))
            
            percent_abundances = []
            for dataset in self.datasets:
                if self.datasets_dict[dataset].has_key(oligo):
                    percent_abundances.append((self.datasets_dict[dataset][oligo] * 100.0 / dataset_sums[dataset],
                                               self.datasets_dict[dataset][oligo],
                                               dataset_sums[dataset],
                                               dataset))

            percent_abundances.sort(reverse = True)

            # NOTE: if a dataset has less than 100 sequences, percent abundance doesn't mean much.
            #       if user wants to eliminate oligotypes that doesn't appear in at least one dataset
            #       more than 1% abundance, a singleton of that oligotype that appears in a dataset
            #       which has 50 sequences would make that oligotype pass the filter. I think if an
            #       oligotype passes the percent filter, dataset size and actual count of the oligotype
            #       should also be considered before considering it as an abundant oligotype:
            for abundance_percent, abundance_count, dataset_size, dataset in percent_abundances:
                PercentAbundance_OK = abundance_percent >= self.min_percent_abundance
                DatesetSize_OK      = dataset_size > 100 or abundance_count > self.min_percent_abundance

                if PercentAbundance_OK and DatesetSize_OK:
                    self.abundant_oligos.append((sum([x[1] for x in percent_abundances]), oligo))
                    break

        self.progress.reset()
        self.run.info('num_oligos_after_a_elim', pretty_print(len(self.abundant_oligos)))
        
        self.abundant_oligos = [x[1] for x in sorted(self.abundant_oligos, reverse = True)]

        # eliminate very rare oligos (the ACTUAL ABUNDANCE, which is the sum of oligotype in all datasets
        # should should be more than 'self.min_actual_abundance'.
        oligos_for_removal = []
        for i in range(0, len(self.abundant_oligos)):
            oligo = self.abundant_oligos[i]

            if i % 100 == 0 or i == len(self.abundant_oligos) - 1:
                 self.progress.update('Applying -A parameter: ' + P(i, len(non_singleton_oligos)))

            oligo_actual_abundance = sum([self.datasets_dict[dataset][oligo] for dataset in self.datasets_dict if self.datasets_dict[dataset].has_key(oligo)])
            if self.min_actual_abundance > oligo_actual_abundance:
                oligos_for_removal.append(oligo)

        for oligo in oligos_for_removal:
            self.abundant_oligos.remove(oligo)
        self.progress.reset()
        self.run.info('num_oligos_after_A_elim', pretty_print(len(self.abundant_oligos)))


        # if 'limit_oligotypes_to' is defined, eliminate all other oligotypes
        if self.limit_oligotypes_to:
            self.abundant_oligos = [oligo for oligo in self.abundant_oligos if oligo in self.limit_oligotypes_to]
            self.run.info('num_oligos_after_l_elim', pretty_print(len(self.abundant_oligos)))
            if len(self.abundant_oligos) == 0:
                raise ConfigError, "Something is wrong; all oligotypes were eliminated with --limit-oligotypes. Quiting."

        self.progress.end()


    def _refine_datasets_dict(self):
        # removing oligos from datasets dictionary that didn't pass
        # MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE and
        # MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS filters.
        self.progress.new('Refining Datasets Dict')

        self.progress.update('Deepcopying datasets dict .. ')
        datasets_dict_copy = copy.deepcopy(self.datasets_dict)
        self.progress.append('done')

        datasets_to_remove = []
        for i in range(0, len(self.datasets)):
            dataset = self.datasets[i]

            self.progress.update('Analyzing datasets: ' + P(i, len(self.datasets) - 1))
            
            for oligo in datasets_dict_copy[dataset]:
                if oligo not in self.abundant_oligos:
                    self.datasets_dict[dataset].pop(oligo)
            if not self.datasets_dict[dataset]:
                datasets_to_remove.append(dataset)
        for dataset in datasets_to_remove:
            self.datasets.remove(dataset)
            self.datasets_dict.pop(dataset)

        self.progress.end()
        
        number_of_reads_in_datasets_dict = sum([sum(self.datasets_dict[dataset].values()) for dataset in self.datasets_dict]) 

        self.run.info('num_sequences_after_qc', '%s of %s (%.2f%%)'\
                            % (pretty_print(number_of_reads_in_datasets_dict),
                               pretty_print(self.fasta.total_seq),
                               number_of_reads_in_datasets_dict * 100.0 / self.fasta.total_seq))

        if len(datasets_to_remove):
            self.run.info('datasets_removed_after_qc', datasets_to_remove)


    def _generate_FASTA_file(self): 
        # store abundant oligos
        self.progress.new('FASTA File')
        oligos_fasta_file_path = self.generate_output_destination("OLIGOS.fasta")
        f = open(oligos_fasta_file_path, 'w')
        self.progress.update('Generating')
        for oligo in self.abundant_oligos:
            f.write('>' + oligo + '\n')
            f.write(oligo + '\n')
        f.close()
        self.progress.end()
        self.run.info('oligos_fasta_file_path', oligos_fasta_file_path)
        
        
    def _generate_NEXUS_file(self):
        # generate NEXUS file of oligos
        self.progress.new('NEXUS File')
        oligos_nexus_file_path = self.generate_output_destination("OLIGOS.nexus")
        f = open(oligos_nexus_file_path, 'w')
        f.write("""begin data;
            dimensions ntax=%d nchar=%d;
            format datatype=dna interleave=no gap=-;
            matrix\n""" % (len(self.abundant_oligos), len(self.abundant_oligos[0])))
        self.progress.update('Generating')
        for oligo in self.abundant_oligos:
            f.write('    %.40s %s\n' % (oligo, oligo))
        f.write('    ;\n')
        f.write('end;\n')
        f.close()
        self.progress.end()
        self.run.info('oligos_nexus_file_path', oligos_nexus_file_path)
    

    def _generate_ENVIRONMENT_file(self):
        # generate environment file
        self.progress.new('ENVIRONMENT File')
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        f = open(environment_file_path, 'w')
        self.progress.update('Generating')
        for dataset in self.datasets:
            for oligo in self.datasets_dict[dataset]:
                f.write("%s\t%s\t%d\n" % (oligo, dataset, self.datasets_dict[dataset][oligo]))
        f.close()
        self.progress.end()
        self.run.info('environment_file_path', environment_file_path)


    def _generate_MATRIX_files(self):
        # generate matrices..
        self.progress.new('Matrix Files')
        matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")
        count_file = open(matrix_count_file_path, 'w')
        percent_file = open(matrix_percent_file_path, 'w')
        
        count_file.write('\t'.join([''] + self.datasets) + '\n')
        percent_file.write('\t'.join([''] + self.datasets) + '\n')
        self.progress.update('Generating')
        for oligo in self.abundant_oligos:
            counts = []
            percents = []
            for dataset in self.datasets:
                if self.datasets_dict[dataset].has_key(oligo):
                    counts.append(str(self.datasets_dict[dataset][oligo]))
                    percents.append(str(self.datasets_dict[dataset][oligo] * 100.0 / sum(self.datasets_dict[dataset].values())))
                else:
                    counts.append('0')
                    percents.append('0.0')
            count_file.write('\t'.join([oligo] + counts) + '\n')
            percent_file.write('\t'.join([oligo] + percents) + '\n')
        count_file.close()
        percent_file.close()
        self.progress.end()
        self.run.info('matrix_count_file_path', matrix_count_file_path)
        self.run.info('matrix_percent_file_path', matrix_percent_file_path)
        
   
    def _generate_viamics_datasets_dict(self):
        self.progress.new('VIAMICS Datasets Dict')
        # generate viamics datasets dict 
        viamics_datasets_dict = {}
        viamics_datasets_dict_file_path = self.generate_output_destination("ENVIRONMENT.cPickle")
        self.progress.update('Generating')
        for dataset in self.datasets_dict:
            viamics_datasets_dict[dataset] = {}
            viamics_datasets_dict[dataset]['species'] = {}
            for oligo in self.datasets_dict[dataset]:
                if oligo in self.abundant_oligos:
                    viamics_datasets_dict[dataset]['species'][oligo] = self.datasets_dict[dataset][oligo]
        
        for dataset in viamics_datasets_dict:
            viamics_datasets_dict[dataset]['tr'] = sum(viamics_datasets_dict[dataset]['species'].values())
            viamics_datasets_dict[dataset]['bases_of_interest_locs'] = self.bases_of_interest_locs
        self.progress.update('Storing')
        cPickle.dump(viamics_datasets_dict, open(viamics_datasets_dict_file_path, 'w'))
        self.progress.end()

        self.run.info('viamics_datasets_dict_file_path', viamics_datasets_dict_file_path)


    def _generate_representative_sequences(self):
        # create a fasta file with a representative full length consensus sequence for every oligotype

        # this is what is going on here: we go through all oligotypes, gather sequences that are being
        # represented by a particular oligotype, unique them and report the top ten unique sequences
        # ordered by the frequency.
        self.progress.new('Represenative Sequences')

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
                                                % (self.fasta.pos * 100 / self.fasta.total_seq))
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
            if oligo in self.abundant_oligos:
                fasta_files_dict[oligo]['file'].write('>%s\n' % (self.fasta.id))
                fasta_files_dict[oligo]['file'].write('%s\n' % self.fasta.seq)
        
        self.progress.end()

        for oligo in self.abundant_oligos:
            self.progress.new('Representative Sequences | %s (%d of %d)'\
                % (oligo, self.abundant_oligos.index(oligo) + 1, len(self.abundant_oligos)))
            fasta_files_dict[oligo]['file'].close()

            fasta_file_path = fasta_files_dict[oligo]['path']
            fasta = u.SequenceSource(fasta_file_path, lazy_init = False, unique = True)
          
            # this dict is going to hold the information of how unique sequences within an oligotype
            # is distributed among datasets:
            distribution_among_datasets = {}

            self.progress.update('Unique reads in FASTA ..') 

            while fasta.next() and fasta.pos <= self.limit_representative_sequences:
                unique_files_dict[oligo]['file'].write('>%s_%d|freq:%d\n'\
                                                                     % (oligo,
                                                                        fasta.pos,
                                                                        len(fasta.ids)))
                unique_files_dict[oligo]['file'].write('%s\n' % fasta.seq)

                for dataset_id in fasta.ids:
                    dataset_name = self.dataset_name_from_defline(dataset_id)
                    if not distribution_among_datasets.has_key(dataset_name):
                        distribution_among_datasets[dataset_name] = {}
                    d = distribution_among_datasets[dataset_name]
                    if not d.has_key(fasta.pos):
                        d[fasta.pos] = 1
                    else:
                        d[fasta.pos] += 1
                
            fasta.close()
            unique_files_dict[oligo]['file'].close()

            unique_fasta_path = unique_files_dict[oligo]['path']
            distribution_among_datasets_dict_path = unique_fasta_path + '_distribution.cPickle'
            cPickle.dump(distribution_among_datasets, open(distribution_among_datasets_dict_path, 'w'))

            if (not self.quick) and (not self.skip_blast_search):
                # perform BLAST search and store results
                unique_fasta = u.SequenceSource(unique_fasta_path)
                unique_fasta.next()
               
                if self.blast_ref_db:
                    # if self.blast_ref_db is set, then perform a local BLAST search 
                    # against self.blast_ref_db
                    oligo_representative_blast_output = unique_fasta_path + '_BLAST.txt'

                    self.progress.update('Local BLAST Search..')
                        
                    local_blast_search(unique_fasta.seq, self.blast_ref_db, oligo_representative_blast_output)

                else:
                    # if self.blast_ref_db is not set, perform a BLAST search on NCBI
                    oligo_representative_blast_output = unique_fasta_path + '_BLAST.xml'

                    # FIXME: this value should be paramaterized
                    max_blast_attempt = 3

                    def blast_search_wrapper(seq, blast_output):
                        try:
                            remote_blast_search(seq, blast_output)
                            return True
                        except:
                            return False

                    for blast_attempt in range(0, max_blast_attempt):
                        self.progress.update('NCBI BLAST search (attempt #%d)' % (blast_attempt + 1))
                            
                        if blast_search_wrapper(unique_fasta.seq, oligo_representative_blast_output):
                            break
                        else:
                            continue

                unique_fasta.close()

            if (not self.quick) and (not self.no_figures):
                entropy_file_path = unique_fasta_path + '_entropy'
                color_per_column_path  = unique_fasta_path + '_color_per_column.cPickle'

                # generate entropy output at 'entropy_file_path' along with the image
                self.progress.update('Generating entropy figure')
                vis_freq_curve(unique_fasta_path, output_file = unique_fasta_path + '.png', entropy_output_file = entropy_file_path)

                # use entropy output to generate a color shade for every columns in alignment
                # for visualization purposes
                entropy_values_per_column = [0] * self.alignment_length
                for column, entropy in [x.strip().split('\t') for x in open(entropy_file_path)]:
                    entropy_values_per_column[int(column)] = float(entropy)
                color_shade_dict = get_color_shade_dict_for_list_of_values(entropy_values_per_column)

                color_per_column = [0] * self.alignment_length
                for i in range(0, self.alignment_length):
                    color_per_column[i] = color_shade_dict[entropy_values_per_column[i]]        

                cPickle.dump(color_per_column, open(color_per_column_path, 'w'))
       
        self.progress.end()
        self.run.info('output_directory_for_reps', output_directory_for_reps) 


    def _generate_random_colors(self):
        random_color_file_path = self.generate_output_destination('COLORS')
        self.colors_dict = random_colors(copy.deepcopy(self.abundant_oligos), random_color_file_path)
        self.run.info('random_color_file_path', random_color_file_path)


    def _generate_dataset_oligotype_network_figures(self):
        output_directory_for_datasets = self.generate_output_destination("DATASETS", directory = True)
        oligotype_network_structure(self.run.info_dict['environment_file_path'], output_dir = output_directory_for_datasets)
        self.run.info('output_directory_for_datasets', output_directory_for_datasets) 
 

    def _generate_stack_bar_figure(self):
        self.progress.new('Stackbar Figure')
        stack_bar_file_path = self.generate_output_destination('STACKBAR.png')
        self.progress.update('Generating')
        oligotype_distribution_stack_bar(self.datasets_dict, self.colors_dict, stack_bar_file_path, oligos = self.abundant_oligos, project_title = self.project, display = not self.no_display)
        self.progress.end()
        self.run.info('stack_bar_file_path', stack_bar_file_path)


    def _generate_html_output(self):
        from utils.html.error import HTMLError
        try:
            from utils.html.generate import generate_html_output
        except HTMLError, e:
            sys.stdout.write('\n\n\t%s\n\n' % e)
            sys.exit()

        self.progress.new('HTML Output')
        output_directory_for_html = self.generate_output_destination("HTML-OUTPUT", directory = True)
        self.progress.update('Generating')
        index_page = generate_html_output(self.run.info_dict, html_output_directory = output_directory_for_html)
        self.progress.end()
        if not self.no_display:
            sys.stdout.write('\n\n\tView results in your browser: "%s"\n\n' % index_page)


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Start an Oligotyping Process')
    parser.add_argument('alignment', metavar = 'INPUT ALIGNMENT',
                        help = 'Alignment file that contains all datasets and sequences in FASTA format')
    parser.add_argument('entropy', metavar = 'ENTROPY',
                        help = 'File that contains the columns and the entropy values computer previously')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-c', '--number-of-auto-components', type=int, default=None,
                        help = 'Number of components to use from alignment to generate oligotypes. Default\
                                is "5", which is a completely arbitrary value. Number of components should\
                                be determined after a careful examination of entropy figure.')
    parser.add_argument('--qual-scores-file', metavar = 'QUAL SCORES FILE',
                        help = 'FASTA formatted file that contains PHRED base call values\
                                for each read in the alignment file')
    parser.add_argument('--qual-scores-dict', metavar = 'QUAL SCORES DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                                PHRED base call values for each read in the alignment file. If you\
                                provide --qual-scores-file, that file will be used to recompute this\
                                dictionary and the file you refer with this parameter will\
                                not be ignored')
    parser.add_argument('--qual-stats-dict', metavar = 'QUAL STATS DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                                PHRED base call quality score statistics for the alignment file. If\
                                you provide --qual-scores-dict, it will be used to recompute this\
                                dictionary and the file you refer to with this parameter will\
                                actually not be used')
    parser.add_argument('-q', '--min-base-quality', type=int, default=15,
                        help = 'Minimum quality score for each base in locations of interest of a read to be\
                                considered in an oligotype. When base quality score files are provided, this\
                                value makes sure that low quality bases that are more likely to be the result\
                                of random sequencing errors do not create artificial oligotypes. Any read that has\
                                less quality score than the given value, will simply be discarded. This parameter\
                                only in effect when --qual-scores-file or --qual-scores-dict parameters are used. \
                                Defeault --min-base-quality is 15.')
    parser.add_argument('-C', '--selected-components', type=str, default=None,
                        help = 'Comma separated entropy components to be used during the oligotyping process.')
    parser.add_argument('-s', '--min-number-of-datasets', type=int, required=True,
                        help = 'Minimum number of datasets oligotype expected to appear. The deafult is "5", which\
                                is another completely arbitrary value. This parameter should be defined based\
                                on the number of datasets included in the analysis. If there are 10 datasets,\
                                3 might be a good choice, if there are 5 datasets, 1 would be a better one\
                                depending on the study.')
    parser.add_argument('-a', '--min-percent-abundance', type=float, default=0.0,
                        help = 'Minimum percent abundance of an oligotype in at least one dataset. The default\
                                is "0.0". Just like --min-number-of-datasets parameter, this parameter too is\
                                to eliminate oligotypes that are formed by sequencing errors occured at the\
                                component of interest. The value should be decided based on the average number\
                                of sequences every dataset has.')
    parser.add_argument('-A', '--min-actual-abundance', type=int, default=0,
                        help = 'Minimum total abundance of an oligotype in all datastes. The default\
                                is "0". If the total abundance of an oligotype is smaller than the number given\
                                with this parameter, oligotype would be eliminated and not included in downstream\
                                analyses.')
    parser.add_argument('-t', '--dataset-name-separator', type=str, default='_',
                        help = 'Character that separates dataset name from unique info in the defline. For insatnce\
                                if the defline says >dataset-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')
    parser.add_argument('-l', '--limit-representative-sequences', type=int, default=None,
                        help = 'At the end of the oligotyping sequences that are being represented by the same\
                                oligotype are being uniqued and stored in separate files. The number of sequences\
                                to keep from the frequency ordered list can be defined with this parameter (e.g.\
                                -l 10 would make it possible that only first 10 sequence would be stored). Default\
                                is 0, which stores everything, but when the dataset size is too big, this could\
                                take up disk space.')
    parser.add_argument('--limit-oligotypes-to', type = str, default = None,
                        help = 'Comma separated list of oligotypes to be taken into account during the analysis.\
                                All other oligotypes will be discarded if a list of oligotypes is being speficied\
                                with this parameter.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'Some relatively insignificant parts of the analysis may take a lot of time, such as\
                                generating figures for representative sequences. When this parameter is set, all\
                                trivial steps would be skipped to give results as soon as possible.')
    parser.add_argument('--no-figures', action = 'store_true', default = False,
                        help = 'When set, no figures will be generated or displayed.')
    parser.add_argument('--blast-ref-db', default = None, type=str,
                        help = 'When set, BLAST search will be done locally against the ref db (local BLAST search\
                                requires USEARCH)')
    parser.add_argument('--skip-blast-search', action = 'store_true', default = False,
                        help = 'When set, BLAST search step will not be performed.')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                        help = 'When set, no figures will be shown.')
    parser.add_argument('--gen-html', action = 'store_true', default = False,
                        help = 'Generate static HTML output to browse analysis results.')
    parser.add_argument('--gen-dataset-oligo-networks', action = 'store_true', default = False,
                        help = 'Generate oligotype network structure figures for each dataset.')
    parser.add_argument('--project', default = None, type=str,
                        help = 'When a project name is set, given name will be used in figures whenever possible.')


    oligotyping = Oligotyping(parser.parse_args())

    try:
        oligotyping.run_all()
    except ConfigError, e:
        print e
        sys.exit()

 
