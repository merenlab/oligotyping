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

import os
import sys
import operator
import cPickle
from scipy import log2 as log

import lib.fastalib as u
from utils.utils import process_command_line_args_for_quality_files
from utils.random_colors import get_list_of_colors
from utils.utils import pretty_print

class EntropyError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Error: %s' % self.e

COLORS = {'A': 'red',
          'T': 'blue', 
          'C': 'green', 
          'G': 'purple', 
          'N': 'white', 
          '-': '#CACACA'}


def entropy(l):
    P = lambda n: (len([x for x in l if x.upper() == n.upper()]) * 1.0 / len(l)) + 0.0000000000000000001
    return -(sum([P(N) * log(P(N)) for N in ['A', 'T', 'C', 'G', '-']]))

def weighted_entropy(column, column_qual, expected_qual_score = 40):
    P = lambda n: (len([x for x in column if x.upper() == n.upper()]) * 1.0 / len(column)) + 0.0000000000000000001
    if column_qual:
        return -(sum([P(N) * log(P(N)) for N in ['A', 'T', 'C', 'G', '-']]) * (column_qual['mean'] / expected_qual_score))
    else:
        return entropy(column)


def entropy_analysis(alignment_path, output_file = None, verbose = True, uniqued = False, freq_from_defline = None, weighted = False, qual_stats_dict = None):
    if freq_from_defline == None:
        freq_from_defline = lambda x: int([t.split(':')[1] for t in x.split('|') if t.startswith('freq')][0])

    lines = []
    previous_alignment_length = None
   
    alignment = u.SequenceSource(alignment_path)

    # processing the alignment file..
    while alignment.next():
        # check the alignment lengths along the way:
        if previous_alignment_length:
            if previous_alignment_length != len(alignment.seq):
                raise EntropyError, "Not all reads have the same length."

        # print out process info
        if verbose:
            if alignment.pos % 10000 == 0:
                sys.stderr.write('\rReading FASTA into memory; reads processed: %s' \
                                % (pretty_print(alignment.pos)))
                sys.stderr.flush()
        
        # fill 'lines' variable
        if not uniqued:
            lines.append(alignment.seq)
        else:
            frequency = freq_from_defline(alignment.id)
            for i in range(0, frequency):
                lines.append(alignment.seq)

        previous_alignment_length = len(alignment.seq)

    alignment.close()
    if verbose:
        sys.stderr.write('\n')

    entropy_tpls = [] 
   
    for position in range(0, len(lines[0])):
        if verbose:
            sys.stderr.write('\rPerforming entropy analysis: %d%%' \
                                % (int((position + 1) * 100.0 / len(lines[0]))))
            sys.stderr.flush()
   
        if len(set([x[position] for x in lines])) == 1:
            entropy_tpls.append((position, 0.0),)
        else:
            column = "".join([x[position] for x in lines])

            if weighted:
                if not qual_stats_dict: 
                    raise EntropyError, "Weighted entropy is selected, but no qual stats are provided"
                e = weighted_entropy(column, qual_stats_dict[position])
            else:
                e = entropy(column)

            if e < 0.00001:
                entropy_tpls.append((position, 0.0),)
            else:
                entropy_tpls.append((position, e),)
    
    if verbose:
        print
   
    if output_file:
        entropy_output = open(output_file, 'w')
        for _component, _entropy in sorted(entropy_tpls, key=operator.itemgetter(1), reverse=True):
            entropy_output.write('%d\t%.4f\n' % (_component, _entropy))
        entropy_output.close()
    
    return [x[1] for x in entropy_tpls]


def get_unique_sequences(alignment, limit = 10):
    unique_sequences = []

    fasta = u.SequenceSource(alignment, unique = True, lazy_init = False)

    while fasta.next() and fasta.pos < limit:
        unique_sequences.append((fasta.seq, len(fasta.ids), len(fasta.ids) / float(fasta.total_seq)))

    return unique_sequences


def visualize_distribution(alignment, entropy_values, output_file, quick = False, no_display = False, qual_stats_dict = None, weighted = False):
    import matplotlib.pyplot as plt
    import numpy as np

    y_maximum = max(entropy_values) + (max(entropy_values) / 10.0)
    number_of_uniques_to_show = int(y_maximum * 100)
    unique_sequences = get_unique_sequences(alignment, limit = number_of_uniques_to_show)

    fig = plt.figure(figsize = (len(unique_sequences[0][0]) / 20, 10))

    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.050, top = 0.92, bottom = 0.10)

    ax = fig.add_subplot(111)

    if not quick:
        current = 0
        for y in range(number_of_uniques_to_show - 1, 0, -3):
            unique_sequence = unique_sequences[current][0]
            count = unique_sequences[current][1]
            frequency = unique_sequences[current][2]
            for i in range(0, len(unique_sequence)):
                plt.text(i, y / 100.0, unique_sequence[i],\
                                    fontsize = 5, color = COLORS[unique_sequence[i]])

            percent = int(round(frequency * len(unique_sequence))) or 1
            plt.fill_between(range(0, percent), (y + 1.15) / 100.0, (y - 0.85) / 100.0, color="green", alpha = 0.2)
            plt.text(percent + 0.8, (y - 1.2) / 100.0, count, fontsize = 5, color = 'gray')

            current += 1
            if current + 1 > len(unique_sequences):
                break

    if not quick and qual_stats_dict:
        # add mean quality values in the background of the figure.
        colors = get_list_of_colors(21, colormap="RdYlGn")
        colors = [colors[0] for _ in range(0, 20)] + colors

        max_count = max([qual_stats_dict[q]['count'] for q in qual_stats_dict if qual_stats_dict[q]])

        for pos in range(0, len(entropy_values)):
            if not qual_stats_dict[pos]:
                continue

            mean = int(round(qual_stats_dict[pos]['mean']))
            count = qual_stats_dict[pos]['count']
            plt.fill_between([pos, pos + 1], y1 = 0, y2 = y_maximum, color = colors[mean], alpha = (log(count) / log(max_count)) / 5)

    ind = np.arange(len(entropy_values))
    ax.bar(ind, entropy_values, color = 'black', lw = 0.5)
    ax.set_xlim([0, len(unique_sequences[0][0])])
    ax.set_ylim([0, y_maximum])
    plt.xlabel('Position in the Alignment')
    if weighted:
        plt.ylabel('Weighted Shannon Entropy')
    else:
        plt.ylabel('Shannon Entropy')
    plt.savefig(output_file)

    if not no_display:
        plt.show()


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Entropy Analysis')
    parser.add_argument('alignment', metavar = 'ALIGNMENT', help = 'Alignment file\
                         that contains all datasets and sequences in FASTA format')
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
    parser.add_argument('--weighted', action = 'store_true', default = False,
                        help = 'When set, entropy computation per column will use\
                        mean quality score for each column.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'When set, entropy values will be shown as fast as\
                                possible (some visualization steps will be skipped).')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                                help = 'When set, no figures will be shown.')

    args = parser.parse_args()


    # process qual scores if provided
    qual_stats_dict = process_command_line_args_for_quality_files(args, _return = 'qual_stats_dict')       

    output_text_path = args.alignment + '%s-ENTROPY' % ('-WEIGHTED' if args.weighted else '')
    output_img_path = args.alignment + '%s-ENTROPY.png' % ('-WEIGHTED' if args.weighted else '')

    entropy_values = entropy_analysis(args.alignment,
                                      output_file = output_text_path,
                                      weighted = args.weighted,
                                      qual_stats_dict = qual_stats_dict)

    visualize_distribution(args.alignment,
                           entropy_values,
                           output_file = output_img_path,
                           quick = args.quick,
                           no_display = args.no_display,
                           qual_stats_dict = qual_stats_dict,
                           weighted = args.weighted)

