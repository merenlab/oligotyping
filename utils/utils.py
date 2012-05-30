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
import shutil
import fcntl
import termios 
import cPickle
import struct
import numpy as np

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from lib import fastalib as u

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
        quals_dict[qual.id] = qual.quals_int

    while alignment.next():
        matching_qual = quals_dict[alignment.id] 

        qual_aligned = []
        for i in range(0, len(alignment.seq)):
            if alignment.seq[i] != '-':
                qual_aligned.append(matching_qual.pop(0))
            else:
                qual_aligned.append(None)

        quals_aligned_dict[alignment.id] = qual_aligned

    if output_file_path:
        cPickle.dump(quals_aligned_dict, open(output_file_path, 'w'))

    return quals_aligned_dict


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
