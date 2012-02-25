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

import fastalib as u

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
        self.output_directory    = None
        self.number_of_components = 5
        self.min_number_of_samples = 5
        self.min_percent_abundance = 1.0
        self.dataset_name_separator = '_'

        if args:
            self.alignment = args.alignment
            self.entropy = args.entropy
            self.output_directory = args.output_directory or os.path.dirname(args.alignment)
            self.number_of_components = args.number_of_components
            self.min_number_of_samples = args.min_number_of_samples
            self.min_percent_abundance = args.min_percent_abundance
            self.dataset_name_separator = args.dataset_name_separator

        self.samples_dict = {}
        self.samples = []
        self.abundant_oligos = []

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


    def generate_output_destination(self, postfix):
        output_file_prefix = '%s-C%d-S%d-A%.1f' % (os.path.basename(self.alignment).split('.')[0],
                                                   self.number_of_components,
                                                   self.min_number_of_samples,
                                                   self.min_percent_abundance)

        return os.path.join(self.output_directory, output_file_prefix + '-' + postfix)


    def info(self, label, value):
        info_line = "%s %s: %s" % (label, '.' * (60 - len(label)), str(value))
        self.info_file_obj.write(info_line + '\n')
        print info_line


    def run_all(self):
        self.sanity_check()
        
        self.info_file_path = self.generate_output_destination('RUNINFO')
        self.info_file_obj = open(self.info_file_path, 'w')

        self.fasta = u.SequenceSource(self.alignment, lazy_init = False)
        self.column_entropy = [int(x.strip().split()[0]) for x in open(self.entropy).readlines()]
       
        self.info('Output directory', self.output_directory)
        self.info('Extraction info output file', self.info_file_path)
        self.info('Input FASTA file', self.alignment)
        self.info('Input entropy file', self.entropy)
        self.info('Number of sequences in FASTA', pp(self.fasta.total_seq))
        self.info('Number of entropy components to use', self.number_of_components)
        self.info('Min number of samples oligotype appears', self.min_number_of_samples)
        self.info('Min % abundance of oligotype in at least one sample', self.min_percent_abundance)
        
        # locations of interest based on the entropy scores
        self.bases_of_interest_locs = sorted([self.column_entropy[i] for i in range(0, self.number_of_components)])
        self.info('Bases of interest', ', '.join([str(x) for x in self.bases_of_interest_locs]))
       
        self._construct_samples_dict()
        self._contrive_abundant_oligos()
        self._refine_samples_dict()
        self._generate_NEXUS_file()
        self._generate_ENVIRONMENT_file()
        self._generate_MATRIX_files()
        self._generate_viamics_samples_dict()
        self._generate_representative_concensus_sequences()

        self.info_file_obj.close()


    def _construct_samples_dict(self):
        self.fasta.reset()
        while self.fasta.next():
            sample = self.dataset_name_from_defline(self.fasta.id)
            
            if not self.samples_dict.has_key(sample):
                self.samples_dict[sample] = {}
                self.samples.append(sample)
        
            oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
        
            if self.samples_dict[sample].has_key(oligo):
                self.samples_dict[sample][oligo] += 1
            else:
                self.samples_dict[sample][oligo] = 1
        self.info('Number of samples in FASTA', pp(len(self.samples_dict)))

    
    def _contrive_abundant_oligos(self):
        # cat oligos | uniq
        oligos_set = []
        for sample in self.samples:
            for oligo in self.samples_dict[sample].keys():
                if oligo not in oligos_set:
                    oligos_set.append(oligo)
        self.info('Number of unique oligotypes', pp(len(oligos_set)))
        
        # count oligo abundance
        oligo_abundance = []
        for oligo in oligos_set:
            count = 0
            for sample in self.samples:
                if oligo in self.samples_dict[sample].keys():
                    count += 1
            oligo_abundance.append((count, oligo),)
        oligo_abundance.sort()
        
        # eliminate singleton/doubleton oligos (any oligo required to appear in at least
        # 'self.min_number_of_samples' samples)
        non_singleton_oligos = []
        for tpl in oligo_abundance:
            if tpl[0] >= self.min_number_of_samples:
                non_singleton_oligos.append(tpl[1])
        self.info('Oligotypes after "min number of samples" elimination', pp(len(non_singleton_oligos)))
        
        # eliminate very rare oligos (the percent abundance of every oligo should be
        # more than 'self.min_percent_abundance' percent in at least one sample)
        for oligo in non_singleton_oligos:
            percent_abundances = []
            for sample in self.samples:
                if self.samples_dict[sample].has_key(oligo):
                    percent_abundances.append(self.samples_dict[sample][oligo] * 100.0 / sum([self.samples_dict[sample][o] for o in non_singleton_oligos if self.samples_dict[sample].has_key(o)]))
            percent_abundances.sort()
            if percent_abundances[-1] >= self.min_percent_abundance:
                self.abundant_oligos.append(oligo)
        self.info('Oligotypes after "min % abundance in a sample" elimination', pp(len(self.abundant_oligos)))
 
        # store abundant oligos
        abundant_oligos_file_path = self.generate_output_destination("OLIGOS.fasta")
        f = open(abundant_oligos_file_path, 'w')
        for oligo in sorted(self.abundant_oligos):
            f.write('>' + oligo + '\n')
            f.write(oligo + '\n')
        f.close()
        self.info('Abundant oligotypes file path', abundant_oligos_file_path)
       

    def _refine_samples_dict(self):
        # removing oligos from samples dictionary that didn't pass
        # MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE and
        # MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS filters.
        samples_dict_copy = copy.deepcopy(self.samples_dict)
        samples_to_remove = []
        for sample in self.samples:
            for oligo in samples_dict_copy[sample]:
                if oligo not in self.abundant_oligos:
                    self.samples_dict[sample].pop(oligo)
            if not self.samples_dict[sample]:
                samples_to_remove.append(sample)
        for sample in samples_to_remove:
            self.samples.remove(sample)
            self.samples_dict.pop(sample)
        if len(samples_to_remove):
            self.info('Samples removed for having 0 oligotypes left after filtering', ', '.join(samples_to_remove))
        
        
    def _generate_NEXUS_file(self):
        # generate NEXUS file of oligos
        oligos_nexus_file_path = self.generate_output_destination("OLIGOS.nexus")
        f = open(oligos_nexus_file_path, 'w')
        f.write("""begin data;
            dimensions ntax=%d nchar=%d;
            format datatype=dna interleave=no;
            matrix\n""" % (len(self.abundant_oligos), len(self.abundant_oligos[0])))
        for oligo in sorted(self.abundant_oligos):
            f.write('    %.20s %s\n' % (oligo, oligo))
        f.write('    ;\n')
        f.write('end;\n')
        f.close()
        self.info('NEXUS file for oligotypes', oligos_nexus_file_path)
    

    def _generate_ENVIRONMENT_file(self):
        # generate environment file
        environment_file_path = self.generate_output_destination("ENVIRONMENT.txt")
        f = open(environment_file_path, 'w')
        for sample in self.samples:
            for oligo in self.samples_dict[sample]:
                f.write("%s\t%s\t%d\n" % (oligo, sample, self.samples_dict[sample][oligo]))
        f.close()
        self.info('Environment file for Viamics/UniFrac analysis', environment_file_path)


    def _generate_MATRIX_files(self):
        # generate matrices..
        matrix_count_file_path = self.generate_output_destination("MATRIX-COUNT.txt")
        matrix_percent_file_path = self.generate_output_destination("MATRIX-PERCENT.txt")
        count_file = open(matrix_count_file_path, 'w')
        percent_file = open(matrix_percent_file_path, 'w')
        
        count_file.write('\t'.join([''] + self.samples) + '\n')
        percent_file.write('\t'.join([''] + self.samples) + '\n')
        for oligo in self.abundant_oligos:
            counts = []
            percents = []
            for sample in self.samples:
                if self.samples_dict[sample].has_key(oligo):
                    counts.append(str(self.samples_dict[sample][oligo]))
                    percents.append(str(self.samples_dict[sample][oligo] * 100.0 / sum(self.samples_dict[sample].values())))
                else:
                    counts.append('0')
                    percents.append('0.0')
            count_file.write('\t'.join([oligo] + counts) + '\n')
            percent_file.write('\t'.join([oligo] + percents) + '\n')
        count_file.close()
        percent_file.close()
        self.info('Data matrix (counts)', matrix_count_file_path)
        self.info('Data matrix (percents)', matrix_percent_file_path)
        
   
    def _generate_viamics_samples_dict(self):
        # generate viamics samples dict 
        viamics_samples_dict = {}
        viamics_samples_dict_file_path = self.generate_output_destination("ENVIRONMENT.cPickle")
        for sample in self.samples_dict:
            viamics_samples_dict[sample] = {}
            viamics_samples_dict[sample]['species'] = {}
            for oligo in self.samples_dict[sample]:
                if oligo in self.abundant_oligos:
                    viamics_samples_dict[sample]['species'][oligo] = self.samples_dict[sample][oligo]
        
        for sample in viamics_samples_dict:
            viamics_samples_dict[sample]['tr'] = sum(viamics_samples_dict[sample]['species'].values())
            viamics_samples_dict[sample]['bases_of_interest_locs'] = self.bases_of_interest_locs
        cPickle.dump(viamics_samples_dict, open(viamics_samples_dict_file_path, 'w'))
        self.info('Serialized Viamics samples dictionary', viamics_samples_dict_file_path)


    def _generate_representative_concensus_sequences(self):
        # create a fasta file with a representative full length consensus sequence for every oligotype

        #
        # FIXME: some of these sequences, especially the ones that are being represented by rare oligotypes
        #        can be chimeric. it is very normal given the nature of the process. but there should be a way
        #        to test, and communicate results back to the user.
        #
        #
        representative_oligotypes_file_path = self.generate_output_destination("REPRESENTATIVE-OLIGO-SEQS.fasta")
        f = open(representative_oligotypes_file_path, 'w')
        for abundant_oligo in self.abundant_oligos:
            self.fasta.reset()
            counter = 0
            rep_sequences = []
            while self.fasta.next() and counter < 5000:
                oligo = ''.join(self.fasta.seq[o] for o in self.bases_of_interest_locs)
                if oligo == abundant_oligo:
                    rep_sequences.append(self.fasta.seq)
                    counter += 1
            
            consensus_dict = {}
            consensus_sequence = ''
            alignment_length = len(rep_sequences[0])

            for i in range(0, alignment_length):
                consensus_dict[i] = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '-': 0}
            for sequence in rep_sequences:
                for pos in range(0, alignment_length):
                    consensus_dict[pos][sequence[pos]] += 1
            for pos in range(0, alignment_length):
                consensus_sequence += sorted(consensus_dict[pos].iteritems(), key=operator.itemgetter(1), reverse=True)[0][0]
            f.write('>' + abundant_oligo + '\n')
            f.write(consensus_sequence + '\n')
        f.close()
        
        # remove uninformative columns from representative full length consensus sequences
        trim_uninformative_columns_from_alignment(representative_oligotypes_file_path)
        self.info('Representative sequences for oligotypes', representative_oligotypes_file_path) 
        
        
 
if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Convert FastQ to FASTA')
    parser.add_argument('-i', '--alignment', required=True, metavar = 'INPUT ALIGNMENT',
                        help = 'Alignment file that contains all samples and sequences in FASTA format')
    parser.add_argument('-e', '--entropy', required=True, metavar = 'ENTROPY',
                        help = 'File that contains the columns and the entropy values computer previously')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = '.')
    parser.add_argument('-c', '--number-of-components', type=int, required=True,
                        help = 'Number of components to use from alignment to generate oligotypes. Default\
                                is "5", which is a completely arbitrary value. Number of components should\
                                be determined after the careful examination of entropy figure.')
    parser.add_argument('-s', '--min-number-of-samples', type=int, required=True,
                        help = 'Minimum number of samples oligotype expected to appear. The deafult is "5", which\
                                is another completely arbitrary value. This parameter should be defined based\
                                on the number of datasets included in the analysis. If there are 10 datasets,\
                                3 might be a good choice, if there are 5 datasets, 1 would be a better one\
                                depending on the study.')
    parser.add_argument('-a', '--min-percent-abundance', type=float, required=True,
                        help = 'Minimum percent abundance of an oligotype in at least one dataset. The default\
                                is "1.0". Just like --min-number-of-samples parameter, this parameter too is\
                                to eliminate oligotypes that are formed by sequencing errors occured at the\
                                component of interest. The value should be decided based on the average number\
                                of sequences every sample has.')
    parser.add_argument('-t', '--dataset-name-separator', type=str, default='_',
                        help = 'Character that separates dataset name from unique info in the defline. For insatnce\
                                if the defline says >dataset-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')

    oligotyping = Oligotyping(parser.parse_args())

    try:
        oligotyping.run_all()
    except ConfigError, e:
        print e
        sys.exit()

 
