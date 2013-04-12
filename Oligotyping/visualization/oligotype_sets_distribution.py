# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import numpy as np
import matplotlib.pyplot as plt
import cPickle

from Oligotyping.utils.random_colors import get_list_of_colors
from Oligotyping.utils.utils import HTMLColorToRGB
from Oligotyping.utils.utils import get_vectors_from_oligotypes_across_samples_matrix


def vis_oligotype_sets_distribution(partitions, vectors, samples, colors_dict = None, output_file = None, legend = False, project_title = None, display = True):
    if colors_dict == None:
        colors_dict = {}
        list_of_colors = get_list_of_colors(len(partitions), colormap = 'Accent')
        for i in range(0, len(partitions)):
            colors_dict[partitions[i][0]] = list_of_colors[i]
            #for oligo in partitions[i]:
            #    print '%s,%s' % (oligo, list_of_colors[i])

    # figure.. 
    plt.figure(figsize=(20, 7))
    
    if legend:
        plt.subplots_adjust(left=0.03, bottom = 0.25, top = 0.95, right = 0.87)
    else:
        plt.subplots_adjust(left=0.03, bottom = 0.25, top = 0.95, right = 0.99)

    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True) 

    N = len(samples)
    ind = np.arange(N)
    width = 0.75
    
    number_of_dimensions = len(vectors.values()[0])

    for i in range(0, len(partitions)):
        group = partitions[i]
        vector = [] 
        mins = []
        maxs = []

        for d in range(0, number_of_dimensions):
            vector.append(np.mean([vectors[oligo][d] for oligo in group]))
            mins.append(np.min([vectors[oligo][d] for oligo in group]))
            maxs.append(np.max([vectors[oligo][d] for oligo in group]))
            
        try:
            color = HTMLColorToRGB(colors_dict[group[0]])
        except:
            color = 'black'

        plt.fill_between(range(0, len(vector)), maxs, mins, color=color, alpha = 0.1)
        plt.plot(vector, color=color, linewidth = 1, alpha = 0.95, label = 'Set #%d' % i)
        if len(vector) < 50:
            plt.plot(vector, color=color, linewidth = 3, alpha = 0.7, label = '_nolegend_')
            plt.plot(vector, color=color, linewidth = 7, alpha = 0.6, zorder = i, label = '_nolegend_')
    
    plt.ylabel('Oligotype Set Abundance', size='large')
    plt.title(project_title if project_title else 'Oligotype Sets Across Samples')

    plt.xticks(ind, samples, rotation=90, size='small')
    plt.yticks([])
    plt.ylim(ymax = 100)
    plt.xlim(xmin = -(width) / 2, xmax = len(samples) - 0.5)
    
    if legend:
        plt.legend(bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0, shadow=True, fancybox=True)
        
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

    parser = argparse.ArgumentParser(description='Visualize Oligotype Sets Across Samples')
    parser.add_argument('environment_file', metavar = 'ENVIRONMENT_FILE',\
                        help = 'Environment file')
    parser.add_argument('sets_file', metavar = 'PARTITIONS_FILE',\
                        help = 'Serialized list of lists for oligotype sets')
    parser.add_argument('oligotypes_across_samples', metavar = 'OLIGOTYPES_ACROSS_DATASETS',\
                        help = 'A TAB-delimited matrix file that contains normalized\
                                oligotype frequencies across samples')
    parser.add_argument('--colors-file', metavar = 'COLORS_FILE', default = None,\
                        help = 'File that contains random colors for oligotypes')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name for the figure to be stored. File name\
                                must end with "png", "jpg", or "tiff".')
    parser.add_argument('--legend', action = 'store_true', default = False,
                        help = 'Turn on legend')
    parser.add_argument('--project-title', default = None, metavar = 'PROJECT_TITLE',\
                        help = 'Project name for the samples.')


    args = parser.parse_args()
    
    samples = []
    for oligotype, sample, count in [line.strip().split('\t') for line in open(args.environment_file).readlines()]:
        if sample in samples:
            continue
        else:
            samples.append(sample)
    samples.sort()

    if args.colors_file:
        colors_dict = {}
        for oligotype, color in [line.strip().split('\t') for line in open(args.colors_file).readlines()]:
            colors_dict[oligotype] = color
    else:
        colors_dict = None

    partitions = cPickle.load(open(args.partitions_file))

    oligos, vectors = get_vectors_from_oligotypes_across_samples_matrix(args.oligotypes_across_samples)

    vis_oligotype_sets_distribution(partitions, vectors, samples, colors_dict, output_file = args.output_file,\
                           legend = args.legend, project_title = args.project_title, display = True)
