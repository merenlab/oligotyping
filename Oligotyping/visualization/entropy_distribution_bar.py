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

import numpy as np
from scipy import log2 as log
import matplotlib.pyplot as plt

from Oligotyping.utils.random_colors import get_list_of_colors
from Oligotyping.utils.utils import get_unique_sequences_from_FASTA
from Oligotyping.utils.utils import Progress
from Oligotyping.utils.utils import Run
from Oligotyping.utils.utils import NUCL_COLORS


run = Run()
progress = Progress()


def entropy_distribution_bar(alignment, entropy_values, output_file, quick = False, no_display = False, qual_stats_dict = None, weighted = False, verbose = False):
    progress.verbose = verbose
    progress.new('Entropy Distribution Figure')
    progress.update('Computing ')

    y_maximum = max(entropy_values) + (max(entropy_values) / 10.0)
    y_maximum = 1 if y_maximum < 1 else y_maximum

    number_of_uniques_to_show = int(y_maximum * 100)

    if alignment == None:
        quick = True

    colors_dict = {}
    if not quick:
        unique_sequences = get_unique_sequences_from_FASTA(alignment, limit = number_of_uniques_to_show)
        
        chars = []
        for seq in unique_sequences:
            chars += seq[0]
        chars = set(chars)
      
        colors_dict = NUCL_COLORS
        
        missing_chars = [char for char in chars if char not in NUCL_COLORS.keys()]
            
        if missing_chars:
            colors_for_missing_chars = get_list_of_colors(len(missing_chars), colormap="RdYlGn")
            for i in range(0, len(missing_chars)):
                char = missing_chars[i]
                colors_dict[char] = colors_for_missing_chars[i]
    else:
        unique_sequences = None

    fig = plt.figure(figsize = (len(entropy_values) / 20, 10))

    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.050, top = 0.92, bottom = 0.10)

    ax = fig.add_subplot(111)

    if not quick:
        current = 0
        for y in range(number_of_uniques_to_show - 1, 0, -3):
            progress.append('.')
            unique_sequence = unique_sequences[current][0].upper()
            count = unique_sequences[current][1]
            frequency = unique_sequences[current][2]
            for i in range(0, len(unique_sequence)):
                plt.text(i, y / 100.0, unique_sequence[i],\
                                    fontsize = 5, color = colors_dict[unique_sequence[i]])

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
    ax.set_xlim([0, len(entropy_values)])
    ax.set_ylim([0, y_maximum])
    plt.xlabel('Position in the Alignment')
    if weighted:
        plt.ylabel('Weighted Shannon Entropy')
    else:
        plt.ylabel('Shannon Entropy')

    progress.update('Saving into "%s"' % output_file)
    plt.savefig(output_file + '.png')
    plt.savefig(output_file + '.pdf')

    if verbose:
        progress.clear()
        run.info('Entropy figure output path', output_file + '.{png, pdf}')

    if not no_display:
        try:
            progress.update('Entropy figure is being shown (you do not have display? you can avoid this step by using --no-display))')
            plt.show()
        except:
            pass

    progress.end()

if __name__ == '__main__':
    pass
