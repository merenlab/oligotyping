#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

__version__ = '0.2' # Nov 08, 2012

import numpy
import operator
from scipy import log2 as log
from numpy import sqrt

import Oligotyping.lib.fastalib as u
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import Run
from Oligotyping.utils.utils import P


class EntropyError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Error: %s' % self.e


VALID_CHARS = {'nucleotide': set(['A', 'T', 'C', 'G', '-']),
               'amino_acid': set(['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V', '-'])}

run = Run()


def entropy(l, l_qual = None, expected_qual_score = 40, amino_acid_sequences = False, sqrt_norm = False):
    l = l.upper() 
    
    valid_chars = VALID_CHARS['amino_acid'] if amino_acid_sequences else VALID_CHARS['nucleotide']

    if sqrt_norm:
        l_normalized = ''
        for char in valid_chars:
            l_normalized += char * int(round(sqrt(l.count(char))))
        l = l_normalized


    E_Cs = []
    for char in valid_chars:
        P_C = (l.count(char) * 1.0 / len(l)) + 0.0000000000000000001
        E_Cs.append(P_C * log(P_C))
   
    if l_qual:
        # return weighted entropy
        return -(sum(E_Cs) * (l_qual['mean'] / expected_qual_score))
    else:
        # return un-weighted entropy
        return -(sum(E_Cs))


def entropy_analysis(alignment_path, output_file = None, verbose = True, uniqued = False, freq_from_defline = None, weighted = False, qual_stats_dict = None, amino_acid_sequences = False):
    if freq_from_defline == None:
        freq_from_defline = lambda x: int([t.split(':')[1] for t in x.split('|') if t.startswith('freq')][0])

    lines = []
    previous_alignment_length = None

    progress = Progress()
    progress.verbose = verbose
   
    alignment = u.SequenceSource(alignment_path)

    progress.new('Processing the Alignment')

    # processing the alignment file..
    while alignment.next():
        # check the alignment lengths along the way:
        if previous_alignment_length:
            if previous_alignment_length != len(alignment.seq):
                raise EntropyError, "Not all reads have the same length."

        # print out process info
        if alignment.pos % 10000 == 0:
            progress.update('Reads processed: %s' % (pretty_print(alignment.pos)))
        
        # fill 'lines' variable
        if not uniqued:
            lines.append(alignment.seq)
        else:
            try:
                frequency = freq_from_defline(alignment.id)
            except IndexError:
                raise EntropyError, "Reads declared as unique, but they do not have proper deflines. See help for --uniqued."
                
            for i in range(0, frequency):
                lines.append(alignment.seq)

        previous_alignment_length = len(alignment.seq)

    progress.end()
    if verbose:
        run.info('Number of reads', pretty_print(alignment.pos))

    alignment.close()


    # entropy analysis
    progress.new('Entropy Analysis')
    entropy_tpls = []

    for position in range(0, len(lines[0])):
        progress.update(P(int(position + 1), len(lines[0])))
   
        if len(set([x[position] for x in lines])) == 1:
            entropy_tpls.append((position, 0.0),)
        else:
            column = "".join([x[position] for x in lines])

            if weighted:
                if not qual_stats_dict: 
                    raise EntropyError, "Weighted entropy is selected, but no qual stats are provided"
                e = entropy(column, l_qual = qual_stats_dict[position], amino_acid_sequences = amino_acid_sequences)
            else:
                e = entropy(column, amino_acid_sequences = amino_acid_sequences)

            if e < 0.00001:
                entropy_tpls.append((position, 0.0),)
            else:
                entropy_tpls.append((position, e),)

    sorted_entropy_tpls = sorted(entropy_tpls, key=operator.itemgetter(1), reverse=True)
    entropy_components_larger_than_0 = [e[1] for e in entropy_tpls if e[1] > 0]

    progress.end()
    if verbose:
        run.info('Entropy analysis', 'Done (total of %d components greater than 0, mean: %.2f, max: %.2f, min: %.2f).' \
                                                        % (len(entropy_components_larger_than_0),
                                                           numpy.mean(entropy_components_larger_than_0),
                                                           numpy.max(entropy_components_larger_than_0),
                                                           numpy.min(entropy_components_larger_than_0)))
   
    if output_file:
        entropy_output = open(output_file, 'w')
        for _component, _entropy in sorted_entropy_tpls:
            entropy_output.write('%d\t%.4f\n' % (_component, _entropy))
        if verbose:
            run.info('Entropy analysis output file path', output_file)
        entropy_output.close()
    
    return [x[1] for x in entropy_tpls]


def quick_entropy(l, amino_acid_sequences = False):
    if len(set([len(x) for x in l])) != 1:
        raise EntropyError, "Not all vectors have the same length."
    
    entropy_tpls = []
    for position in range(0, len(l[0])):
        if len(set([x[position] for x in l])) == 1:
            entropy_tpls.append((position, 0.0),)
        else:
            column = "".join([x[position] for x in l])

            e = entropy(column, amino_acid_sequences = amino_acid_sequences)

            if e < 0.00001:
                entropy_tpls.append((position, 0.0),)
            else:
                entropy_tpls.append((position, e),)

    return [x[1] for x in entropy_tpls if x[1] > 0]
