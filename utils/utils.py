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
import fcntl
import time
import shutil
import termios 
import cPickle
import struct
import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from lib import fastalib as u
from constants import pretty_names

P = lambda x, y: '%.2f%%' % (x * 100.0 / y)


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

        for abundance_percent, abundance_count, dataset_size, dataset in percent_abundances:
            abundant_oligos.append((sum([x[1] for x in percent_abundances]), oligo))
            break

    return [x[1] for x in sorted(abundant_oligos)]


def get_qual_stats_dict(quals_dict, output_file_path = None):
    """This function takes quals dict (which can be obtained by calling the
       utils.utils.get_quals_dict function) and returns a dictionary that
       simply contains the summary of quality scores per location in the
       alignment"""

    # FIXME: get_quals_dict and get_qual_stats_dict functions are only for
    #        454 technology at this moment.

    qual_stats_dict = {}
    alignment_length = len(quals_dict[quals_dict.keys()[0]])
    for pos in range(0, alignment_length):

        sys.stderr.write('\r    Qual stats are being computed: %d of %d' % (pos + 1, alignment_length))
        sys.stderr.flush()

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

    sys.stderr.write('\n')
    return qual_stats_dict
  
def get_quals_dict(quals_file, alignment_file, output_file_path = None):
    """This function takes qual scores file in FASTA format, expands each
       entry to match base calls in the corresponding aligned read in the
       FASTA file (which requires deflines to be identical), and finally
       returns a dictionary that contains qual scores as a list of integer
       values that are bound to deflines as key/value pairs"""

    quals_dict = {}
    quals_aligned_dict = {}
 
    alignment = u.SequenceSource(alignment_file)
    qual = u.QualSource(quals_file)

    while qual.next():
        if qual.pos % 1000 == 0:
            sys.stderr.write('\r    Quals dict is being generated. Step 1 of 2; pos: %s' % (pretty_print(qual.pos)))
            sys.stderr.flush()
        quals_dict[qual.id] = qual.quals_int
    sys.stderr.write('\n')

    while alignment.next():
        if alignment.pos % 1000 == 0:
            sys.stderr.write('\r    Quals dict is being generated. Step 2 of 2; pos: %s' % (pretty_print(alignment.pos)))
            sys.stderr.flush()

        matching_qual = quals_dict[alignment.id] 

        qual_aligned = []
        for i in range(0, len(alignment.seq)):
            if alignment.seq[i] != '-':
                qual_aligned.append(matching_qual.pop(0))
            else:
                qual_aligned.append(None)

        quals_aligned_dict[alignment.id] = qual_aligned
    sys.stderr.write('\n')

    if output_file_path:
        cPickle.dump(quals_aligned_dict, open(output_file_path, 'w'))

    return quals_aligned_dict

def process_command_line_args_for_quality_files(args, _return = 'qual_stats_dict'):
    """this function computes and returns the dictionary of interest (indicated
       with '_return' parameter if qual score files were provided via the command
       line interface.
       
       _return value expected to be either 'qual_stats_dict', or 'quals_dict'.

       """

    if _return not in ['qual_stats_dict', 'quals_dict']:
        return None

    if args.qual_scores_file:
        sys.stderr.write('* Generating quality scores dictionary..\n')
        quals_dict = get_quals_dict(args.qual_scores_file,\
                                    args.alignment,\
                                    output_file_path = args.qual_scores_file + '.cPickle')
        if _return == 'quals_dict':
            return quals_dict

        sys.stderr.write('* Computing quality stats dictionary file from quality scores dictionary.\n')
        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                                              output_file_path = args.qual_scores_file + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_scores_dict:
        sys.stderr.write('* Reading quality scores dictionary..\n')
        quals_dict = cPickle.load(open(args.qual_scores_dict))

        if _return == 'quals_dict':
            return quals_dict

        sys.stderr.write('* Computing qual stats dict file from quals dict..\n')
        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                            output_file_path = args.qual_scores_dict.split('.cPickle')[0] + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_stats_dict:
        sys.stderr.write('* Reading qual stats dictionary..\n')
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
    return time.strftime("%d %b %Y, %H:%M:%S", time.gmtime())

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
            cr = (env['LINES'], env['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])


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


