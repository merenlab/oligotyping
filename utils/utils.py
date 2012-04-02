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

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
from lib import fastalib as u


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
