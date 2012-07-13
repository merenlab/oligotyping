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
import cPickle
import numpy as np
import matplotlib.pyplot as plt

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))

import lib.fastalib as u
from lib.entropy import entropy_analysis

def vis_freq_curve(fasta_file_path, output_file = None, x_limit = 20, display = False, freq_from_defline = None, entropy_output_file = None, verbose = False):
    if freq_from_defline == None:
        freq_from_defline = lambda x: int([t.split(':')[1] for t in x.split('|') if t.startswith('freq')][0])

    fasta = u.SequenceSource(fasta_file_path)

    frequency_list = []
    while fasta.next():
        try:
            frequency_list.append(freq_from_defline(fasta.id)) 
        except:
            print 'frequency info can not be read from defline.'
            sys.exit()

    frequency_list_to_plot = frequency_list[0:x_limit] + [0] * (x_limit - len(frequency_list) \
                                            if len(frequency_list) < x_limit else 0)

    
    fig = plt.figure(figsize=(24, 10))

    plt.subplots_adjust(left=0.05, bottom = 0.15, top = 0.95, right = 0.99)
  
    plt.subplot(2, 1, 1)
    plt.grid(True) 
    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)
  
    plt.plot(frequency_list_to_plot, lw = 3, c = 'black')
 
    plt.xlabel('Order in the File', size = 'x-large')
    plt.ylabel('Frequency of the Unique Sequence', size = 'x-large')
    plt.title('Frequency Distribution of Unique Sequences in %s' % os.path.basename(fasta_file_path))
    plt.ylim(ymin = -max(frequency_list_to_plot) * 0.05, ymax = max(frequency_list_to_plot) * 1.05)
    plt.xlim(xmin = -0.05, xmax = x_limit - 1)
    plt.xticks(range(0, x_limit), [str(i) for i in range(1, x_limit + 1)], rotation=90, size='small')
    

    plt.subplot(2, 1, 2)
    plt.subplots_adjust(left=0.05, bottom = 0.1, top = 0.95, right = 0.99)

    try:
        plt.grid(axis='y') 
    except:
        plt.grid(True)
    plt.rcParams.update({'axes.linewidth' : 0.9})
    plt.rc('grid', color='0.40', linestyle='-', linewidth=0.1)

    entropy_values = entropy_analysis(fasta_file_path, output_file = entropy_output_file, verbose = verbose, uniqued = True)

    y_maximum = max(entropy_values) * 1.1
    y_maximum = 1.1 if y_maximum < 1 else y_maximum
    ind = np.arange(len(entropy_values))
    plt.bar(ind, entropy_values, color = 'black', lw = 0.5)
    plt.xlim([0, len(entropy_values)])
    plt.ylim([0, y_maximum])
    plt.xticks( range(0, len(entropy_values), 5), rotation=90, size = 'x-small')
    
    plt.xlabel('Position in the Alignment', size = 'x-large')
    plt.ylabel('Shannon Entropy', size = 'x-large')

    if output_file:
        plt.savefig(output_file)
    if display:
        plt.show()

    plt.clf()
    plt.close('all')

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate Distribution of Unique Sequences Figure')
    parser.add_argument('fasta', metavar = 'FASTA_FILE', help = 'Sequences file in FASTA format')
    parser.add_argument('-x', '--x-limit', default = 20, type = int, metavar = 'X_LIMIT',\
                        help = 'Number of items to show from frequency list')

    fasta_file_path = parser.parse_args().fasta
    x_limit = parser.parse_args().x_limit
    vis_freq_curve(fasta_file_path, x_limit = x_limit, output_file = fasta_file_path + '.png', display = True, verbose = True)

