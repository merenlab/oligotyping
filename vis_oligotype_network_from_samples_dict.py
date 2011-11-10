#!/usr/bin/python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2011, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import pylab
import matplotlib.pyplot as plt
import numpy as np
import sys
import copy


env = sys.argv[1]

tuples = [line.strip().split("\t") for line in open(env)]

original_samples_dict = {}

base_pos = {'-': 5, 'A': 4, 'T': 3, 'C': 2, 'G': 1}
base_colors = {'-': 'white', 'A': 'red', 'T': 'green', 'C': 'blue', 'G': 'yellow'}

for tpl in tuples:
    oligo, sample, count = tpl

    m = []
    for base in oligo:
        m.append(base_pos[base])

    if original_samples_dict.has_key(sample):
        original_samples_dict[sample][oligo] = (m, int(count))
    else:
        original_samples_dict[sample] = {oligo: (m, int(count))}

def get_REAL_samples_dict(samples_dict, title, condition):
    filtered_samples_dict = {}
    filtered_samples_dict[title] = {}

    temp = {}
    for sample in samples_dict:
        if condition(sample):
            temp[sample] = samples_dict[sample]

    oligos = []
    for l in [x.keys() for x in temp.values()]:
        for o in l:
            oligos.append(o) 

    for oligo in set(oligos):
        total_oligo = sum([t[oligo][1] for t in [x for x in temp.values()] if t.has_key(oligo)])
        m = []
        for base in oligo:
            m.append(base_pos[base])

        filtered_samples_dict[title][oligo] = (m, total_oligo)

    return filtered_samples_dict

def get_PERCENT_samples_dict(samples_dict):
    filtered_samples_dict = {}

    for sample in samples_dict:
        temp = {}
        total_reads = sum([t[1] for t in samples_dict[sample].values()])
        for oligo in samples_dict[sample]:
            temp[oligo] = (samples_dict[sample][oligo][0], (samples_dict[sample][oligo][1] * 100 / total_reads) + 1)
        filtered_samples_dict[sample] = temp

    return filtered_samples_dict


samples_dict = original_samples_dict

print samples_dict

for sample in samples_dict:
    total_reads = sum([x[1] for x in samples_dict[sample].values()])
    print sample, total_reads

    N = len(samples_dict[sample].keys()[0])
    ind = np.arange(N) + 1

    fig = plt.figure(figsize = (N + 2, 6))

    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.025, top = 0.92, bottom = 0.05)

      #left  = 0.125   the left side of the subplots of the figure
      #right = 0.9     the right side of the subplots of the figure
      #bottom = 0.1    the bottom of the subplots of the figure
      #top = 0.9       the top of the subplots of the figure
      #wspace = 0.2    the amount of width reserved for blank space between subplots
      #hspace = 0.2    the amount of height reserved for white space between subplots

    ax = fig.add_subplot(111)
    ax.plot([0], [0], visible = False)
    ax.plot([N], [0], visible = False)
    ax.plot([0], [6], visible = False)


    for pos in range(0, len(oligo)):
        bases = {}
        for oligo in samples_dict[sample]:
            base = oligo[pos]
            if bases.has_key(base):
                bases[base] += samples_dict[sample][oligo][1]
            else:
                bases[base] = samples_dict[sample][oligo][1]
        for base in bases:
            ratio = bases[base] * 1.0 / total_reads
            ax.plot([pos + 1], [base_pos[base]], 'o', c = 'white', lw = 1, ls="--", alpha = 0.75, ms = ratio * 100)
            ax.plot([pos + 1], [base_pos[base]], 'o', c = base_colors[base], lw = 1, alpha = ratio / 5, ms = ratio * 100)


    for oligo in samples_dict[sample]:
        ratio = samples_dict[sample][oligo][1] * 1.0 / total_reads
        ax.plot(ind, samples_dict[sample][oligo][0], c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 8, alpha = ratio / 5)
        ax.plot(ind, samples_dict[sample][oligo][0], c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 6, alpha = ratio / 3)
        ax.plot(ind, samples_dict[sample][oligo][0], c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 4, alpha = ratio)
        ax.plot(ind, samples_dict[sample][oligo][0], c = 'black', solid_capstyle = "round", solid_joinstyle = "round", lw = 1, alpha = 0.1)
        ax.plot(ind, samples_dict[sample][oligo][0], c = 'white', solid_capstyle = "round", solid_joinstyle = "round", lw = 1, alpha = ratio)

    pylab.yticks(np.arange(6), ('', 'G', 'C', 'T', 'A', '--'), size = 'x-large')
    pylab.title(sample + " (total reads: %s)" % total_reads)

    locs = range(0, N + 2)
    pylab.xticks(locs, [''] + ["VL " + str(x) for x in range(0, len(locs))[1:-1]] + [''])

    plt.savefig(sys.argv[1] + '-' + sample + ".png")
