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

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'lib'))
import fastalib as u

from visualization.frequency_curve_and_entropy import vis_freq_curve
from visualization.oligotype_distribution_stack_bar import oligotype_distribution_stack_bar
from utils.random_colors import random_colors
from utils.constants import pretty_names

def pp(n):
    """Pretty print function for very big numbers.."""
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

def trim_uninformative_columns_from_alignment(input_file_path):
    input_fasta = u.SequenceSource(input_file_path, lazy_init = False)
    input_fasta.next()
    invalid_columns = range(0, len(input_fasta.seq))
    input_fasta.reset()
    
    while input_fasta.next():
        for i in invalid_columns:
            if input_fasta.seq[i] != '-':
                invalid_columns.remove(i)
    
    columns_to_keep = [x for x in range(0, invalid_columns[-1]) if x not in invalid_columns]
    
    input_fasta.reset()

    temp_file = tempfile.NamedTemporaryFile(delete = False)
    temp_file_path = temp_file.name
    temp_file.close()

    temp_file = u.FastaOutput(temp_file_path)

    while input_fasta.next():
        new_seq = ''
        for i in columns_to_keep:
            new_seq += input_fasta.seq[i]
        temp_file.write_id(input_fasta.id)
        temp_file.write_seq(new_seq, split = False)
    
    temp_file.close()

    # overwrite the original file with trimmed content
    shutil.move(temp_file_path, input_file_path)


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Config Error: %s' % self.e

class Oligotyping:
    def __init__(self, args = None):
        self.alignment = None
        self.entropy   = None
        self.project = None
        self.output_directory = None
        self.number_of_components = 5
        self.min_number_of_datasets = 5
        self.min_percent_abundance = 1.0
        self.dataset_name_separator = '_'
        self.limit_representative_sequences = sys.maxint

        Absolute = lambda x: os.path.join(os.getcwd(), x) if not x.startswith('/') else x 

        if args:
            self.entropy = Absolute(args.entropy)
            self.alignment = Absolute(args.alignment)
            self.number_of_components = args.number_of_components
            self.min_number_of_datasets = args.min_number_of_datasets
            self.min_percent_abundance = args.min_percent_abundance
            self.project = args.project or os.path.basename(args.alignment).split('.')[0]
            self.output_directory = args.output_directory or os.path.join(os.getcwd(), '-'.join([self.project, self.get_prefix()]))
            self.dataset_name_separator = args.dataset_name_separator
            self.limit_representative_sequences = args.limit_representative_sequences or sys.maxint
            self.quick = args.quick
            self.no_figures = args.no_figures
            self.no_display = args.no_display
            self.gen_html = args.gen_html

        self.datasets_dict = {}
        self.datasets = []
        self.abundant_oligos = []
        self.colors_dict = None
        self.run_info_dict = {}

    def sanity_check(self):
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


    def dataset_name_from_defline(self, defline):
        return self.dataset_name_separator.join(defline.split('|')[0].split(self.dataset_name_separator)[0:-1])


    def get_prefix(self):
        return 'C%d-S%d-A%.1f' % (self.number_of_components,
                                  self.min_number_of_datasets,
                                  self.min_percent_abundance)


    def generate_output_destination(self, postfix, directory = False):
        return_path = os.path.join(self.output_directory, postfix)

        if directory == True:
            if os.path.exists(return_path):
                shutil.rmtree(return_path)
            os.makedirs(return_path)

        return return_path


    def info(self, key, value):
        if pretty_names.has_key(key):
            label = pretty_names[key]
        else:
            label = key

        self.run_info_dict[key] = value

        info_line = "%s %s: %s\n" % (label, '.' * (65 - len(label)), str(value))
        self.info_file_obj.write(info_line + '\n')

        if not self.no_display:
            sys.stdout.write(info_line)


    def run_all(self):
        self.sanity_check()
        
        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.info_file_obj = open(self.info_file_path, 'w')

        self.fasta = u.SequenceSource(self.alignment, lazy_init = False)
        self.column_entropy = [int(x.strip().split()[0]) for x in open(self.entropy).readlines()]

        self.info('project', self.project)
        self.info('alignment', self.alignment)
        self.info('entropy', self.entropy)
        self.info('output_directory', self.output_directory)
        self.info('info_file_path', self.info_file_path)
        self.info('cmd_line', ' '.join(sys.argv))
        self.info('total_seq', pp(self.fasta.total_seq))
        self.info('number_of_components', self.number_of_components)
        self.info('s', self.min_number_of_datasets)
        self.info('a', self.min_percent_abundance)
        
        # locations of interest based on the entropy scores
        self.bases_of_interest_locs = sorted([self.column_entropy[i] for i in range(0, self.number_of_components)])
        self.info('bases_of_interest_locs', ', '.join([str(x) for x in self.bases_of_interest_locs]))

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
        if not self.no_figures:
            self._generate_stack_bar_figure()

        self.info_file_obj.close()
        self._store_run_info_dict()

        if self.gen_html:
            self._generate_html_output()


    def _store_run_info_dict(self):
        info_dict_file_path = self.generate_output_destination("RUNINFO.cPickle")
        cPickle.dump(self.run_info_dict, open(info_dict_file_path, 'w'))


    def _construct_datasets_dict(self):
        self.fasta.reset()
        while self.fasta.next():
            dataset = self.dataset_name_from_defline(self.fasta.id)
            
            if not self.datasets_dict.has_key(dataset):
                self.datasets_dict[dataset] = {}
                self.datasets.append(dataset)
        
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
        
            if self.datasets_dict[dataset].has_key(oligo):
                self.datasets_dict[dataset][oligo] += 1
            else:
                self.datasets_dict[dataset][oligo] = 1
        self.info('num_datasets_in_fasta', pp(len(self.datasets_dict)))

    
    def _contrive_abundant_oligos(self):
        # cat oligos | uniq
        oligos_set = []
        for dataset in self.datasets:
            for oligo in self.datasets_dict[dataset].keys():
                if oligo not in oligos_set:
                    oligos_set.append(oligo)
        self.info('num_unique_oligos', pp(len(oligos_set)))
        
        # count oligo abundance
        oligo_abundance = []
        for oligo in oligos_set:
            count = 0
            for dataset in self.datasets:
                if oligo in self.datasets_dict[dataset].keys():
                    count += 1
            oligo_abundance.append((count, oligo),)
        oligo_abundance.sort()

        # eliminate singleton/doubleton oligos (any oligo required to appear in at least
        # 'self.min_number_of_datasets' datasets)
        non_singleton_oligos = []
        for tpl in oligo_abundance:
            if tpl[0] >= self.min_number_of_datasets:
                non_singleton_oligos.append(tpl[1])
        self.info('num_oligos_after_s_elim', pp(len(non_singleton_oligos)))
        
        # eliminate very rare oligos (the percent abundance of every oligo should be
        # more than 'self.min_percent_abundance' percent in at least one dataset)
        SUM = lambda dataset: sum([self.datasets_dict[dataset][o] for o in non_singleton_oligos \
                                                                if self.datasets_dict[dataset].has_key(o)])
        for oligo in non_singleton_oligos:
            percent_abundances = []
            for dataset in self.datasets:
                if self.datasets_dict[dataset].has_key(oligo):
                    percent_abundances.append((self.datasets_dict[dataset][oligo] * 100.0 / SUM(dataset), self.datasets_dict[dataset][oligo], SUM(dataset), dataset))
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

        self.abundant_oligos = [x[1] for x in sorted(self.abundant_oligos, reverse = True)]

        self.info('num_oligos_after_a_elim', pp(len(self.abundant_oligos)))
      

    def _refine_datasets_dict(self):
        # removing oligos from datasets dictionary that didn't pass
        # MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE and
        # MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS filters.
        datasets_dict_copy = copy.deepcopy(self.datasets_dict)
        datasets_to_remove = []
        for dataset in self.datasets:
            for oligo in datasets_dict_copy[dataset]:
                if oligo not in self.abundant_oligos:
                    self.datasets_dict[dataset].pop(oligo)
            if not self.datasets_dict[dataset]:
                datasets_to_remove.append(dataset)
        for dataset in datasets_to_remove:
            self.datasets.remove(dataset)
            self.datasets_dict.pop(dataset)

        number_of_reads_in_datasets_dict = sum([sum(self.datasets_dict[dataset].values()) for dataset in self.datasets_dict]) 

        self.info('num_sequences_after_qc', '%s of %s (%.2f%%)'\
                            % (pp(number_of_reads_in_datasets_dict),
                               pp(self.fasta.total_seq),
                               number_of_reads_in_datasets_dict * 100.0 / self.fasta.total_seq))

        if len(datasets_to_remove):
            self.info('datasets_removed_after_qc', datasets_to_remove)
        else:
            self.info('datasets_removed_after_qc', datasets_to_remove)


    def _generate_FASTA_file(self): 
        # store abundant oligos
        oligos_fasta_file_path = self.generate_output_destination("OLIGOS.fasta")
        f = open(oligos_fasta_file_path, 'w')
        for oligo in self.abundant_oligos:
            f.write('>' + oligo + '\n')
            f.write(oligo + '\n')
        f.close()
        self.info('oligos_fasta_file_path', oligos_fasta_file_path)
        
        
    def _generate_NEXUS_file(self):
        # generate NEXUS file of oligos
        oligos_nexus_file_path = self.generate_output_destination("OLIGOS.nexus")
        f = open(oligos_nexus_file_path, 'w')
        f.write("""begin data;
            dimensions ntax=%d nchar=%d;
            format datatype=dna interleave=no;
            matrix\n""" % (len(self.abundant_oligos), len(self.abundant_oligos[0])))
        for oligo in self.abundant_oligos:
            f.write('    %.20s %s\n' % (oligo, oligo))
        f.write('    ;\n')
        f.write('end;\n')
        f.close()
        self.info('oligos_nexus_file_path', oligos_nexus_file_path)
    

    def _generate_ENVIRONMENT_file(self):
        # generate environment file
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        f = open(environment_file_path, 'w')
        for dataset in self.datasets:
            for oligo in self.datasets_dict[dataset]:
                f.write("%s\t%s\t%d\n" % (oligo, dataset, self.datasets_dict[dataset][oligo]))
        f.close()
        self.info('environment_file_path', environment_file_path)


    def _generate_MATRIX_files(self):
        # generate matrices..
        matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")
        count_file = open(matrix_count_file_path, 'w')
        percent_file = open(matrix_percent_file_path, 'w')
        
        count_file.write('\t'.join([''] + self.datasets) + '\n')
        percent_file.write('\t'.join([''] + self.datasets) + '\n')
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
        self.info('matrix_count_file_path', matrix_count_file_path)
        self.info('matrix_percent_file_path', matrix_percent_file_path)
        
   
    def _generate_viamics_datasets_dict(self):
        # generate viamics datasets dict 
        viamics_datasets_dict = {}
        viamics_datasets_dict_file_path = self.generate_output_destination("ENVIRONMENT.cPickle")
        for dataset in self.datasets_dict:
            viamics_datasets_dict[dataset] = {}
            viamics_datasets_dict[dataset]['species'] = {}
            for oligo in self.datasets_dict[dataset]:
                if oligo in self.abundant_oligos:
                    viamics_datasets_dict[dataset]['species'][oligo] = self.datasets_dict[dataset][oligo]
        
        for dataset in viamics_datasets_dict:
            viamics_datasets_dict[dataset]['tr'] = sum(viamics_datasets_dict[dataset]['species'].values())
            viamics_datasets_dict[dataset]['bases_of_interest_locs'] = self.bases_of_interest_locs
        cPickle.dump(viamics_datasets_dict, open(viamics_datasets_dict_file_path, 'w'))
        self.info('viamics_datasets_dict_file_path', viamics_datasets_dict_file_path)

    def _generate_representative_sequences(self):
        # create a fasta file with a representative full length consensus sequence for every oligotype

        # this is what is going on here: we go through all oligotypes, gather sequences that are being
        # represented by a particular oligotype, unique them and report the top ten unique sequences
        # ordered by the frequency.

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
                                       'stored. You can learn the actual limit by typing "ulimit -n" in the console.',
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
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
            if oligo in self.abundant_oligos:
                fasta_files_dict[oligo]['file'].write('>%s\n' % (self.fasta.id))
                fasta_files_dict[oligo]['file'].write('%s\n' % self.fasta.seq)
            
        for oligo in self.abundant_oligos:
            fasta_files_dict[oligo]['file'].close()

            fasta_file_path = fasta_files_dict[oligo]['path']
            fasta = u.SequenceSource(fasta_file_path, lazy_init = False, unique = True)

            while fasta.next() and fasta.pos <= self.limit_representative_sequences:
                unique_files_dict[oligo]['file'].write('>%s_%d|freq:%d\n'\
                                                                     % (oligo,
                                                                        fasta.pos,
                                                                        len(fasta.ids)))
                unique_files_dict[oligo]['file'].write('%s\n' % fasta.seq)
        
            fasta.close()
            unique_files_dict[oligo]['file'].close()

            if (not self.quick) and (not self.no_figures):
                unique_fasta_path = unique_files_dict[oligo]['path']
                vis_freq_curve(unique_fasta_path, output_file = unique_fasta_path + '.png')
        
        self.info('output_directory_for_reps', output_directory_for_reps) 


    def _generate_random_colors(self):
        random_color_file_path = self.generate_output_destination('COLORS')
        self.colors_dict = random_colors(copy.deepcopy(self.abundant_oligos), random_color_file_path)
        self.info('random_color_file_path', random_color_file_path)


    def _generate_stack_bar_figure(self):
        stack_bar_file_path = self.generate_output_destination('STACKBAR.png')
        oligotype_distribution_stack_bar(self.datasets_dict, self.colors_dict, stack_bar_file_path, oligos = self.abundant_oligos, project_title = self.project, display = not self.no_display)
        self.info('stack_bar_file_path', stack_bar_file_path)

    def _generate_html_output(self):
        from utils.html.error import HTMLError
        try:
            from utils.html.generate import generate_html_output
        except HTMLError, e:
            sys.stdout.write('\n\n\t%s\n\n' % e)
            sys.exit()

        output_directory_for_html = self.generate_output_destination("HTML-OUTPUT", directory = True)
        index_page = generate_html_output(self.run_info_dict, html_output_directory = output_directory_for_html)
        if not self.no_display:
            sys.stdout.write('\n\n\tView results in your browser: "%s"\n\n' % index_page)

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Start an Oligotyping Process')
    parser.add_argument('-i', '--alignment', required=True, metavar = 'INPUT ALIGNMENT',
                        help = 'Alignment file that contains all datasets and sequences in FASTA format')
    parser.add_argument('-e', '--entropy', required=True, metavar = 'ENTROPY',
                        help = 'File that contains the columns and the entropy values computer previously')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-c', '--number-of-components', type=int, required=True,
                        help = 'Number of components to use from alignment to generate oligotypes. Default\
                                is "5", which is a completely arbitrary value. Number of components should\
                                be determined after the careful examination of entropy figure.')
    parser.add_argument('-s', '--min-number-of-datasets', type=int, required=True,
                        help = 'Minimum number of datasets oligotype expected to appear. The deafult is "5", which\
                                is another completely arbitrary value. This parameter should be defined based\
                                on the number of datasets included in the analysis. If there are 10 datasets,\
                                3 might be a good choice, if there are 5 datasets, 1 would be a better one\
                                depending on the study.')
    parser.add_argument('-a', '--min-percent-abundance', type=float, required=True,
                        help = 'Minimum percent abundance of an oligotype in at least one dataset. The default\
                                is "1.0". Just like --min-number-of-datasets parameter, this parameter too is\
                                to eliminate oligotypes that are formed by sequencing errors occured at the\
                                component of interest. The value should be decided based on the average number\
                                of sequences every dataset has.')
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
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'Some relatively insignificant parts of the analysis may take a lot of time, such as\
                                generating figures for representative sequences. When this parameter is set, all\
                                trivial steps would be skipped to give results as soon as possible.')
    parser.add_argument('--no-figures', action = 'store_true', default = False,
                        help = 'When set, no figures will be generated or displayed.')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                        help = 'When set, no figures will be shown.')
    parser.add_argument('--gen-html', action = 'store_true', default = False,
                        help = 'Generate static HTML output to browse analysis results.')
    parser.add_argument('--project', default = None, type=str,
                        help = 'When a project name is set, given name will be used in figures whenever possible.')


    oligotyping = Oligotyping(parser.parse_args())

    try:
        oligotyping.run_all()
    except ConfigError, e:
        print e
        sys.exit()

 
