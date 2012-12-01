# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2011, A. Murat Eren
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
import numpy as np
import matplotlib.pyplot as plt
import cPickle

from Oligotyping.utils.random_colors import random_colors
from Oligotyping.utils.utils import HTMLColorToRGB
from Oligotyping.utils.utils import get_oligos_sorted_by_abundance
from Oligotyping.utils.utils import get_datasets_dict_from_environment_file


def oligotype_distribution_stack_bar(datasets_dict, colors_dict, output_file = None, legend = False,\
                                     colors_export = None, project_title = None, display = True, oligos = None):
    datasets = datasets_dict.keys()
    datasets.sort()
   
    if oligos == None:
        oligos = get_oligos_sorted_by_abundance(datasets_dict, oligos)
    else:
        oligos.reverse()
 
    if colors_dict == None:
        colors_dict = random_colors(copy.deepcopy(oligos))

    datasets_oligo_vectors = {}
    for dataset in datasets:
        vector = []
        for oligo in oligos:
            if datasets_dict[dataset].has_key(oligo):
                vector.append(datasets_dict[dataset][oligo])
            else:
                vector.append(0)
        datasets_oligo_vectors[dataset] = vector
    
    datasets_oligo_vectors_percent_normalized = {}
    for dataset in datasets:
        total_oligos_in_dataset = sum(datasets_oligo_vectors[dataset])
        vector = []
        for oligo_abundance in datasets_oligo_vectors[dataset]:
            vector.append(oligo_abundance * 100.0 / total_oligos_in_dataset)
        datasets_oligo_vectors_percent_normalized[dataset] = vector
   
    # figure.. 
    fig = plt.figure(figsize=(20, 10))
    
    if legend:
        plt.subplots_adjust(left=0.03, bottom = 0.15, top = 0.97, right = 0.90)
    else:
        plt.subplots_adjust(left=0.03, bottom = 0.15, top = 0.97, right = 0.99)
    
    
    N = len(datasets)
    ind = np.arange(N)
    width = 0.75
    
    bars = []
    colors_list = []

    for i in range(0, len(oligos)):
        values = [datasets_oligo_vectors_percent_normalized[dataset][i] for dataset in datasets]
        bottom = [sum(datasets_oligo_vectors_percent_normalized[dataset][0:i]) for dataset in datasets]
        try:
            color = HTMLColorToRGB(colors_dict[oligos[i]])
            colors_list.append(colors_dict[oligos[i]])
        except:
            color = 'black'
            colors_list.append('#000000')
   

        p = plt.bar(ind, values, width, bottom=bottom, color=color)
        bars.append(p)

    if colors_export:
        colors_list = reversed(colors_list)
        colors_file = open(colors_export, 'w')
        for c in colors_list:
            colors_file.write('%s\n' % c)
        colors_file.close()

    plt.ylabel('Oligotype Distribution', size='large')
    plt.title('Stacked Bar Charts of Oligotype Distribution %s' \
                 % (('for "%s"' % project_title) if project_title else ''))

    plt.xticks(ind+width/2., datasets, rotation=90, size='small')
    plt.yticks([])
    plt.ylim(ymax = 100)
    plt.xlim(xmin = -(width) / 2, xmax = len(datasets))
    
    if legend:
        plt.legend([b[0] for b in bars][::-1], oligos[::-1], bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0, shadow=True, fancybox=True)
        
        leg = plt.gca().get_legend()
        ltext  = leg.get_texts()
        llines = leg.get_lines()
        frame  = leg.get_frame()
        
        frame.set_facecolor('0.80')
        plt.setp(ltext, fontsize='small', fontname='arial', family='monospace')
        plt.setp(llines, linewidth=1.5)
    
    if output_file:
        plt.savefig(output_file)
    if display:
        try:
            plt.show()
        except:
            pass

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Stack Bar Representation of Oligotype Distribution')
    parser.add_argument('environment_file', metavar = 'ENVIRONMENT_FILE',\
                        help = 'Oligotype distribution in datasets')
    parser.add_argument('--colors-file', metavar = 'COLORS_FILE', default = None,\
                        help = 'Two column file that contains colors for oligotypes')
    parser.add_argument('--color-list-file', metavar = 'COLORS_FILE', default = None,\
                        help = 'Single column file that contains a list of colors')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name for the figure to be stored. File name\
                                must end with "png", "jpg", or "tiff".')
    parser.add_argument('--legend', action = 'store_true', default = False,
                        help = 'Turn on legend')
    parser.add_argument('--colors-export', metavar = 'COLORS_LIST_FILE',
                        help = 'Store the color list into a file')
    parser.add_argument('--project-title', default = None, metavar = 'PROJECT_TITLE',\
                        help = 'Project name for the datasets.')


    args = parser.parse_args()

    datasets_dict = get_datasets_dict_from_environment_file(args.environment_file)

    if args.colors_file:
        colors_dict = {}
        for oligotype, color in [line.strip().split('\t') for line in open(args.colors_file).readlines()]:
            colors_dict[oligotype] = color
    elif args.color_list_file:
        colors_dict = {}
        colors = [line.strip() for line in open(args.color_list_file).readlines()]
        oligos = get_oligos_sorted_by_abundance(datasets_dict, None)
        oligos.reverse()
        if len(oligos) > len(colors):
            sys.stderr.write('Error: number of colors in file is less than number of oligos. Quiting.\n')
            sys.exit()
        for oligo in oligos:
            colors_dict[oligo] = colors[oligos.index(oligo)]
    else:
        colors_dict = None

    oligotype_distribution_stack_bar(datasets_dict,
                                     colors_dict,
                                     output_file = args.output_file,
                                     legend = args.legend,
                                     colors_export = args.colors_export,
                                     project_title = args.project_title)










