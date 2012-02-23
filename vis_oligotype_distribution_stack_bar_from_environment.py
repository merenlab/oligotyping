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

# samples_dict_serialized as a parameter..
samples_dict = cPickle.load(open(sys.argv[1])) 

# taxon_color_dict from viamics
oligos_color_dict = cPickle.load(open(sys.argv[2]))
# FIXME: This viamics dependency should not be here.

samples = [s for s in samples_dict.keys()]
samples.sort()

oligos = []
for oligo_list in [d['species'].keys() for d in samples_dict.values()]:
    oligos += oligo_list

oligos = list(set(oligos))
oligos.sort()

samples_oligo_vectors = {}
for sample in samples:
    o = []
    for oligo in oligos:
        if samples_dict[sample]['species'].has_key(oligo):
            o.append(samples_dict[sample]['species'][oligo])
        else:
            o.append(0)
    samples_oligo_vectors[sample] = o


samples_oligo_vectors_percent_normalized = {}
for sample in samples:
    total_oligos_in_sample = sum(samples_oligo_vectors[sample])
    o = []
    for oligo_abundance in samples_oligo_vectors[sample]:
        o.append(oligo_abundance * 100.0 / total_oligos_in_sample)
    samples_oligo_vectors_percent_normalized[sample] = o

fig = plt.figure(figsize=(20, 10))

plt.subplots_adjust(left=0.03, bottom = 0.15, top = 0.97, right = 0.90)

fig.canvas.set_window_title('Sjiz')

N = len(samples)
ind = np.arange(N)
width = 0.75

bars = []

f = []
for i in range(0, len(oligos)):
    values = [samples_oligo_vectors_percent_normalized[sample][i] for sample in samples]
    bottom = [sum(samples_oligo_vectors_percent_normalized[sample][0:i]) for sample in samples]
    try:
        c = oligos_color_dict['species'][oligos[i]]
    except:
        c = 'black'
    p = plt.bar(ind, values, width, bottom=bottom, color=c)
    bars.append(p)

plt.ylabel('Oligotype Distribution')
plt.title('Stacked Bar Charts of Oligotype Distribution Among Samples')
plt.xticks(ind+width/2., samples, rotation=90, size='xx-small')
plt.yticks([])
plt.ylim(ymax = 100)
plt.xlim(xmin = -(width) / 2, xmax = len(samples))

plt.legend([b[0] for b in bars][::-1], oligos[::-1], bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0, shadow=True, fancybox=True)

leg = plt.gca().get_legend()
ltext  = leg.get_texts()
llines = leg.get_lines()
frame  = leg.get_frame()

frame.set_facecolor('0.80')
plt.setp(ltext, fontsize='small', fontname='arial', family='monospace')
plt.setp(llines, linewidth=1.5)

plt.savefig(sys.argv[1] + '-SBC.png')
plt.show()
