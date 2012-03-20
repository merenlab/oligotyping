# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

from matplotlib import mpl
import numpy as np
import matplotlib.pyplot as plt
import pylab
import cPickle
import sys
from scipy import stats

def HTMLColorToRGB(colorstring):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    return (r / 255.0, g / 255.0, b / 255.0)

def oligotype_distribution_stack_bar(datasets_dict, colors_dict, output_file = None, legend = False, project_title = None, display = True, oligos = None):
    datasets = datasets_dict.keys()
    datasets.sort()
   
    if oligos == None:
        oligos = []
        map(lambda o: oligos.extend(o), [v.keys() for v in datasets_dict.values()])
        oligos = sorted(list(set(oligos)))
    else:
        oligos.reverse()
  
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
    
    f = []
    for i in range(0, len(oligos)):
        values = [datasets_oligo_vectors_percent_normalized[dataset][i] for dataset in datasets]
        bottom = [sum(datasets_oligo_vectors_percent_normalized[dataset][0:i]) for dataset in datasets]
        try:
            color = HTMLColorToRGB(colors_dict[oligos[i]])
        except:
            color = 'black'

        p = plt.bar(ind, values, width, bottom=bottom, color=color)
        bars.append(p)
    
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
    parser.add_argument('environment_file', metavar = 'ENVIRONMENT_FILE', help = 'Oligotype distribution in datasets')
    parser.add_argument('colors_file', metavar = 'COLORS_FILE', help = 'File that contains random colors for oligotypes')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File name for the figure to be stored. File name must end with "png", "jpg", or "tiff".')
    parser.add_argument('--legend', action = 'store_true', default = False,
                        help = 'Turn on legend')
    parser.add_argument('--project-title', default = None, metavar = 'PROJECT_TITLE',\
                        help = 'Project name for the datasets.')


    args = parser.parse_args()

    datasets_dict = {}
    for oligotype, dataset, count in [line.strip().split('\t') for line in open(args.environment_file).readlines()]:
        if datasets_dict.has_key(dataset):
            if datasets_dict[dataset].has_key(oligotype):
                datasets_dict[dataset][oligotype] += int(count)
            else:
                datasets_dict[dataset][oligotype] = int(count)
        else:
            datasets_dict[dataset] = {}
            datasets_dict[dataset][oligotype] = int(count)

    colors_dict = {}
    for oligotype, color in [line.strip().split('\t') for line in open(args.colors_file).readlines()]:
        colors_dict[oligotype] = color

    oligotype_distribution_stack_bar(datasets_dict, colors_dict, legend = args.legend, project_title = args.project_title)










