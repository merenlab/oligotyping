# -*- coding: utf-8 -*-

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
import time
import math
import fcntl
import shutil
import struct
import termios 
import cPickle
import tempfile
import numpy as np

from Oligotyping.lib import fastalib as u
from Oligotyping.utils.constants import pretty_names

P = lambda x, y: '%.2f%%' % (x * 100.0 / y)


class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Config Error: %s' % self.e


def generate_MATRIX_files(units, obj):
    # generate matrices..
    obj.progress.new('Matrix Files')
    matrix_count_file_path = obj.generate_output_destination("MATRIX-COUNT.txt")
    matrix_percent_file_path = obj.generate_output_destination("MATRIX-PERCENT.txt")
   
    unit_percents = {}
    unit_counts = {}

    obj.progress.update('Generating the data')
    for dataset in obj.datasets:
        counts = []
        percents = []
        for unit in units:
            if obj.datasets_dict[dataset].has_key(unit):
                counts.append(obj.datasets_dict[dataset][unit])
                percents.append(obj.datasets_dict[dataset][unit] * 100.0 / sum(obj.datasets_dict[dataset].values()))
            else:
                counts.append(0)
                percents.append(0.0)
                
        unit_counts[dataset] = counts
        unit_percents[dataset] = percents
    

    obj.progress.update('Generating Matrix Counts and Percents')
    count_file = open(matrix_count_file_path, 'w')
    percent_file = open(matrix_percent_file_path, 'w')       
    
    count_file.write('\t'.join(['samples'] + units) + '\n')
    percent_file.write('\t'.join(['samples'] + units) + '\n')

    for dataset in obj.datasets:
        count_file.write('\t'.join([dataset] + [str(c) for c in unit_counts[dataset]]) + '\n')
        percent_file.write('\t'.join([dataset] + [str(p) for p in unit_percents[dataset]]) + '\n')
    
    count_file.close()
    percent_file.close()

    obj.progress.end()
    obj.run.info('matrix_count_file_path', matrix_count_file_path)
    obj.run.info('matrix_percent_file_path', matrix_percent_file_path)

    #
    # yes, I know git has branches. I am just being lazy. let me be.
    #
    #obj.progress.new('Cytoscape Network Files')
    #obj.progress.new('Being generated')
    #cytoscape_edges_file_path = obj.generate_output_destination("CYTOSCAPE-EDGES.txt")
    #cytoscape_nodes_file_path = obj.generate_output_destination("CYTOSCAPE-NODES.txt")

    #cytoscape_edges_file = open(cytoscape_edges_file_path, 'w')
    #cytoscape_nodes_file = open(cytoscape_nodes_file_path, 'w')

    #cytoscape_edges_file.write('from\tto\tweight\tconsensus\n')
    #for dataset in obj.datasets:
    #    for i in range(0, len(units)):
    #        if unit_percents[dataset][i]:
    #            cytoscape_edges_file.write('%s\t%s\t%.4f\t%s\n' % (dataset,
    #                                                               obj.abundant_units[i],
    #                                                               unit_percents[dataset][i],
    #                                                               units[i]))

    #cytoscape_nodes_file.write('source\tinteraction\ttarget\n')
    #for source in obj.datasets:
    #    for unit in obj.datasets_dict[source]:
    #        for target in obj.datasets:
    #            if target == source:
    #                continue
    #            if obj.datasets_dict[target].has_key(unit):
    #                cytoscape_nodes_file.write('%s\t%s\t%s\n' % (source,
    #                                                             unit,
    #                                                             target))
    #for unit in units: 
    #    cytoscape_nodes_file.write('%s\t%s\t%s\n' % (unit,
    #                                                   'unit',
    #                                                   obj.final_unit_counts_dict[unit]))
    #                                                   
    #cytoscape_edges_file.close()
    #cytoscape_nodes_file.close()
    #obj.progress.end()
    #obj.run.info('cytoscape_edges_file_path', cytoscape_edges_file_path)
    #obj.run.info('cytoscape_nodes_file_path', cytoscape_nodes_file_path)

    if obj.generate_sets:
        unit_type = 'Oligos' if obj.analysis == 'oligotyping' else 'Nodes'
        obj.progress.new('Matrix Files For %s Across Datasets' % unit_type)
        across_datasets_MN_file_path = obj.generate_output_destination("%s-ACROSS-DATASETS-MAX-NORM.txt" % unit_type.upper())
        across_datasets_SN_file_path = obj.generate_output_destination("%s-ACROSS-DATASETS-SUM-NORM.txt" % unit_type.upper())

        across_datasets_MN_file = open(across_datasets_MN_file_path, 'w')
        across_datasets_SN_file = open(across_datasets_SN_file_path, 'w')

        across_datasets_MN_file.write('\t'.join(['sample'] + units) + '\n')
        across_datasets_SN_file.write('\t'.join(['sample'] + units) + '\n')

        obj.progress.update('Generating data for %s across datasets' % unit_type)
  
        for unit in units:
            obj.across_datasets_sum_normalized[unit] = []
            obj.across_datasets_max_normalized[unit] = []

        for i in range(0, len(units)):
            unit = units[i]
            sum_across_datasets = sum([unit_percents[dataset][i] for dataset in obj.datasets])
            max_across_datasets = max([unit_percents[dataset][i] for dataset in obj.datasets])
            for dataset in obj.datasets:
                obj.across_datasets_sum_normalized[unit].append(unit_percents[dataset][i]  * 100.0 / sum_across_datasets)
                obj.across_datasets_max_normalized[unit].append(unit_percents[dataset][i]  * 100.0 / max_across_datasets)

        obj.progress.update('Generating files')
        for i in range(0, len(obj.datasets)):
            dataset = obj.datasets[i]
            across_datasets_MN_file.write('\t'.join([dataset] + [str(obj.across_datasets_max_normalized[unit][i]) for unit in units]) + '\n')
            across_datasets_SN_file.write('\t'.join([dataset] + [str(obj.across_datasets_sum_normalized[unit][i]) for unit in units]) + '\n')
        
        across_datasets_MN_file.close()
        across_datasets_SN_file.close()

        obj.progress.end()
        obj.run.info('across_datasets_MN_file_path', across_datasets_MN_file_path)
        obj.run.info('across_datasets_SN_file_path', across_datasets_SN_file_path)

 
def unique_and_store_alignment(alignment_path, output_path):
    output = u.FastaOutput(output_path)
    alignment = u.SequenceSource(alignment_path, unique = True)
        
    alignment.next()
    most_abundant_unique_read = alignment.seq
    alignment.reset()
 
    unique_read_counts = []
    while alignment.next():
        unique_read_counts.append(len(alignment.ids))
        output.store(alignment, split = False)
            
    output.close()
    alignment.close()
        
    return (unique_read_counts, most_abundant_unique_read)


def generate_ENVIRONMENT_file(obj):
    # generate environment file
    obj.progress.new('ENVIRONMENT File')
    environment_file_path = obj.generate_output_destination("ENVIRONMENT.txt")
    f = open(environment_file_path, 'w')
    obj.progress.update('Being generated')
    for dataset in obj.datasets:
        for unit in obj.datasets_dict[dataset]:
            f.write("%s\t%s\t%d\n" % (unit, dataset, obj.datasets_dict[dataset][unit]))
    f.close()
    obj.progress.end()
    obj.run.info('environment_file_path', environment_file_path)


def get_unique_sequences_from_FASTA(alignment, limit = 10):
    unique_sequences = []

    fasta = u.SequenceSource(alignment, unique = True, lazy_init = False)

    while fasta.next() and fasta.pos < limit:
        unique_sequences.append((fasta.seq, len(fasta.ids), len(fasta.ids) / float(fasta.total_seq)))

    return unique_sequences


def get_oligos_sorted_by_abundance(datasets_dict, oligos = None):
    datasets = datasets_dict.keys()
    datasets.sort()

    if oligos == None:
        oligos = []
        map(lambda o: oligos.extend(o), [v.keys() for v in datasets_dict.values()])
        oligos = list(set(oligos))

    abundant_oligos = []
    
    SUM = lambda dataset: sum([datasets_dict[dataset][o] for o in oligos \
                                                if datasets_dict[dataset].has_key(o)])
    for oligo in oligos:
        percent_abundances = []

        for dataset in datasets:
            if datasets_dict[dataset].has_key(oligo):
                percent_abundances.append((datasets_dict[dataset][oligo] * 100.0 / SUM(dataset),\
                                           datasets_dict[dataset][oligo], SUM(dataset), dataset))

        percent_abundances.sort(reverse = True)

        # FIXME: excuse me, WTF is going on here?:
        for abundance_percent, abundance_count, dataset_size, dataset in percent_abundances:
            abundant_oligos.append((sum([x[1] for x in percent_abundances]), oligo))
            break

    return [x[1] for x in sorted(abundant_oligos)]


def get_vectors_from_oligotypes_across_datasets_matrix(file_path):
    oligotypes_across_datasets_file_obj = open(file_path)
   
    oligos = []
    vectors = {}

    for line in oligotypes_across_datasets_file_obj.readlines()[1:]:
        fields = line.strip().split('\t')
        
        oligo = fields[0]
        oligos.append(oligo)
        vectors[oligo] = [float(c) for c in fields[1:]]

    return (oligos, vectors)


def get_qual_stats_dict(quals_dict, output_file_path = None, verbose = True):
    """This function takes quals dict (which can be obtained by calling the
       utils.utils.get_quals_dict function) and returns a dictionary that
       simply contains the summary of quality scores per location in the
       alignment"""

    # FIXME: get_quals_dict and get_qual_stats_dict functions are only for
    #        454 technology at this moment.

    progress = Progress()
    progress.verbose = verbose
    progress.new('Summary of quality scores per column is being computed')
    
    qual_stats_dict = {}
    alignment_length = len(quals_dict[quals_dict.keys()[0]])
    for pos in range(0, alignment_length):
        progress.update('Position: %d of %d' % (pos + 1, alignment_length))

        qual_stats_dict[pos] = {}
        quals_for_pos = [q[pos] for q in quals_dict.values() if q[pos]]
        if not quals_for_pos:
            qual_stats_dict[pos] = None
            continue
        qual_stats_dict[pos]['mean']  = np.mean(quals_for_pos)
        qual_stats_dict[pos]['std']   = np.std(quals_for_pos)
        qual_stats_dict[pos]['max']   = np.max(quals_for_pos)
        qual_stats_dict[pos]['min']   = np.min(quals_for_pos)
        qual_stats_dict[pos]['count'] = len(quals_for_pos)
    
    if output_file_path:
        cPickle.dump(quals_dict, open(output_file_path, 'w'))

    progress.end()
    return qual_stats_dict
  
def get_quals_dict(quals_file, alignment_file, output_file_path = None, verbose = True):
    """This function takes qual scores file in FASTA format, expands each
       entry to match base calls in the corresponding aligned read in the
       FASTA file (which requires deflines to be identical), and finally
       returns a dictionary that contains qual scores as a list of integer
       values that are bound to deflines as key/value pairs"""

    quals_dict = {}
    quals_aligned_dict = {}

    progress = Progress()
    progress.verbose = verbose
    progress.new('Quality scores dictionary is being generated')
 
    alignment = u.SequenceSource(alignment_file)
    qual = u.QualSource(quals_file)

    while qual.next():
        if qual.pos % 1000 == 0:
            progress.update('Step 1 of 2 :: Quality scores read: %s' % (pretty_print(qual.pos)))
        quals_dict[qual.id] = qual.quals_int

    while alignment.next():
        if alignment.pos % 1000 == 0:
            progress.update('Step 2 of 2 :: Alignments matched: %s' % (pretty_print(alignment.pos)))
            sys.stderr.flush()

        matching_qual = quals_dict[alignment.id] 

        qual_aligned = []
        for i in range(0, len(alignment.seq)):
            if alignment.seq[i] != '-':
                qual_aligned.append(matching_qual.pop(0))
            else:
                qual_aligned.append(None)

        quals_aligned_dict[alignment.id] = qual_aligned
    progress.end()

    if output_file_path:
        cPickle.dump(quals_aligned_dict, open(output_file_path, 'w'))

    return quals_aligned_dict

def process_command_line_args_for_quality_files(args, _return = 'qual_stats_dict', verbose = True):
    """this function computes and returns the dictionary of interest (indicated
       with '_return' parameter if qual score files were provided via the command
       line interface.
       
       _return value expected to be either 'qual_stats_dict', or 'quals_dict'.

       """

    progress = Progress()
    progress.verbose = verbose

    if _return not in ['qual_stats_dict', 'quals_dict']:
        return None

    if args.qual_scores_file:
        quals_dict = get_quals_dict(args.qual_scores_file,\
                                    args.alignment,\
                                    output_file_path = args.qual_scores_file + '.cPickle')
        if _return == 'quals_dict':
            return quals_dict

        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                                              output_file_path = args.qual_scores_file + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_scores_dict:
        quals_dict = cPickle.load(open(args.qual_scores_dict))

        if _return == 'quals_dict':
            return quals_dict

        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                            output_file_path = args.qual_scores_dict.split('.cPickle')[0] + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_stats_dict:
        qual_stats_dict = cPickle.load(open(args.qual_stats_dict))
        
        if _return == 'qual_stats_dict':
            return qual_stats_dict
    
    else:
        return None
 

def get_datasets_dict_from_environment_file(environment_file_path):
    datasets_dict = {}
    for oligo, dataset, count in [l.strip().split('\t') for l in open(environment_file_path).readlines()]:
        if datasets_dict.has_key(dataset):
            datasets_dict[dataset][oligo] = int(count)
        else:
            datasets_dict[dataset] = {oligo: int(count)}
    return datasets_dict


def human_readable_number(n):
    postfix = ['','K','M']
    level = max(0, min(len(postfix) - 1, int(math.floor(math.log10(abs(n))/3.0))))
    return '%.0f%s' % (n/10 ** (3 * level), postfix[level])


def pretty_print(n):
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

def get_date():
    return time.strftime("%d %b %Y, %H:%M:%S", time.localtime())

def get_terminal_size():
    """function was taken from http://stackoverflow.com/a/566752"""
    def ioctl_GWINSZ(fd):
        try:
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])


def estimate_expected_max_frequency_of_an_erronous_unique_sequence(number_of_reads, average_read_length, expected_error = 1/250.0):
    # maximum number of occurence of an error driven unique sequence among N reads.
    # of course this maximum assumes that all reads are coming from one template,
    # substitution probabilities are homogeneous and there are no systemmatical errors,
    # so it is a mere approximation, but for our purpose, it is going to be enough:

    return round((expected_error * (1 / 3.0)) * ((1 - expected_error) ** (average_read_length - 1)) * number_of_reads) 


def HTMLColorToRGB(colorstring):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    return (r / 255.0, g / 255.0, b / 255.0)


def colorize(txt):
    return '\033[0;30m\033[46m%s\033[0m' % txt


class Progress:
    def __init__(self):
        self.pid = None
        self.verbose = True

    def new(self, pid):
        if self.pid:
            self.end()
        self.pid = pid

    def write(self, c):
        if self.verbose:
            sys.stderr.write(c)
            sys.stderr.flush()

    def reset(self):
        self.write('\r' + ' ' * get_terminal_size()[0])
        self.write('\r')
    
    def append(self, msg):
        self.write(colorize('%s' % (msg)))

    def update(self, msg):
        self.write('\r' + colorize(' ' * get_terminal_size()[0]))
        self.write(colorize('\r[%s] %s' % (self.pid, msg)))
    
    def end(self):
        self.reset()
        self.pid = None


class Run:
    """a class that keeps info about an oligotyping run, and deal with the console output"""
    def __init__(self, info_file_path = None, verbose = True):
        if info_file_path:
            self.init_info_file_obj(info_file_path)
        else:
            self.info_file_obj = None

        self.info_dict = {}
        self.verbose = verbose

    def init_info_file_obj(self, info_file_path):
            self.info_file_obj = open(info_file_path, 'w')

    def info(self, key, value):
        if pretty_names.has_key(key):
            label = pretty_names[key]
        else:
            label = key

        self.info_dict[key] = value

        info_line = "%s %s: %s\n" % (label, '.' * (65 - len(label)), str(value))
        if self.info_file_obj:
            self.info_file_obj.write(info_line)

        if self.verbose:
            sys.stderr.write(info_line)

    def store_info_dict(self, destination):
        cPickle.dump(self.info_dict, open(destination, 'w'))

    def quit(self):
        if self.info_file_obj:
            self.info_file_obj.close()


