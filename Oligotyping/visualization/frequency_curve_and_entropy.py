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
import numpy as np
import matplotlib.pyplot as plt

import Oligotyping.lib.fastalib as u
from Oligotyping.lib.entropy import entropy_analysis

def vis_freq_curve(fasta_file_path, output_file = None, x_limit = 20, display = False, freq_from_defline = None, entropy_output_file = None, verbose = False, mini = False, title = None):
    if freq_from_defline == None:
        freq_from_defline = lambda x: int([t.split(':')[1] for t in x.split('|') if t.startswith('freq')][0])

    fasta = u.SequenceSource(fasta_file_path)

    frequency_list = []
    while next(fasta):
        try:
            frequency_list.append(freq_from_defline(fasta.id)) 
        except:
            print('frequency info can not be read from defline.')
            sys.exit()

    frequency_list_to_plot = frequency_list[0:x_limit] + [0] * (x_limit - len(frequency_list) \
                                            if len(frequency_list) < x_limit else 0)


    entropy_values = entropy_analysis(fasta_file_path, output_file = entropy_output_file, verbose = verbose, uniqued = True)
    
    if mini:
        plt.figure(figsize=(2, 2))
        plt.subplots_adjust(left=0.01, bottom = 0, top = 1, right = 1)
        plt.subplot(1, 1, 1)
        plt.grid(False)
        plt.xticks([])
        plt.yticks([])
        
        ax=plt.gca()
        plt.setp(ax, frame_on=False)
        
        y_maximum = 1.1
        x_maximum = len(entropy_values)
        ind = np.arange(len(entropy_values))

        text_x, text_y = x_maximum / 2, y_maximum / 2
        
        plt.text(text_x, text_y, title if title else 'title',
                        horizontalalignment='center',
                        verticalalignment='center',
                        backgroundcolor='white',
                        fontsize=40, color='red')
        
        plt.ylim(ymax = y_maximum)
        plt.xlim(xmax = x_maximum)
        
        plt.bar(ind, entropy_values, color = 'black', lw = 0.5)
        
    else:
        plt.figure(figsize=(24, 10))
        plt.subplots_adjust(left=0.05, bottom = 0.15, top = 0.95, right = 0.99)
        plt.subplot(2, 1, 1)
        plt.grid(True) 
        plt.rcParams.update({'axes.linewidth' : 0.9})
        plt.rc('grid', color='0.50', linestyle='-', linewidth=0.1)
        plt.xticks( list(range(0, len(entropy_values), 5)), rotation=90, size = 'x-small')
      
        plt.plot(frequency_list_to_plot, lw = 3, c = 'black')
     
        plt.xlabel('Order in the File', size = 'x-large')
        plt.ylabel('Frequency of the Unique Sequence', size = 'x-large')
        if title:
            plt.title(title)
        else:
            plt.title('Frequency Distribution of Unique Sequences in %s' % os.path.basename(fasta_file_path))
        plt.ylim(ymin = -max(frequency_list_to_plot) * 0.05, ymax = max(frequency_list_to_plot) * 1.05)
        plt.xlim(xmin = -0.05, xmax = x_limit - 1)
        plt.xticks(list(range(0, x_limit)), [str(i) for i in range(1, x_limit + 1)], rotation=90, size='small')
        
    
        plt.subplot(2, 1, 2)
        plt.subplots_adjust(left=0.05, bottom = 0.1, top = 0.95, right = 0.99)
    
        try:
            plt.grid(axis='y') 
        except:
            plt.grid(True)
        plt.rcParams.update({'axes.linewidth' : 0.9})
        plt.rc('grid', color='0.40', linestyle='-', linewidth=0.1)
        
        y_maximum = max(entropy_values) * 1.1
        y_maximum = 1.1 if y_maximum < 1 else y_maximum
        ind = np.arange(len(entropy_values))
        plt.bar(ind, entropy_values, color = 'black', lw = 0.5)
        plt.xlim([0, len(entropy_values)])
        plt.ylim([0, y_maximum])
        plt.xticks( list(range(0, len(entropy_values), 5)), rotation=90, size = 'x-small')
        
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
    parser.add_argument('--mini', action = 'store_true', default = False,
                        help = 'Generate mini images')
    parser.add_argument('--title', metavar = 'TITLE_TEXT', help = 'Title to appear on top of the\
                         figure')

    args = parser.parse_args()
    
    vis_freq_curve(args.fasta, x_limit = args.x_limit, display = True, verbose = True, mini = args.mini, title = args.title)

