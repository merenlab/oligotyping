# -*- coding: utf-8 -*-

from matplotlib import mpl
import numpy as np
import matplotlib.pyplot as plt
import pylab
import cPickle
import sys
from scipy import stats

samples_dict = cPickle.load(open(sys.argv[1])) #samples_dict_serialized as a parameter..
oligos_color_dict = cPickle.load(open(sys.argv[2])) # taxon_color_dict from viamics

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


# lets seee..
fig = plt.figure(figsize=(20, 10))

plt.subplots_adjust(left=0.03, bottom = 0.15, top = 0.97, right = 0.80)

fig.canvas.set_window_title('Sjiz')

N = len(samples)
ind = np.arange(N)    # the x locations for the groups
width = 0.75  # the width of the bars: can also be len(x) sequence

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
plt.xticks(ind+width/2., samples, rotation=90)
plt.yticks([])
plt.ylim(ymax = 100)
plt.xlim(xmin = -(width) / 2, xmax = len(samples))

plt.legend([b[0] for b in bars][::-1], oligos[::-1], bbox_to_anchor=(1.01, 1), loc=2, borderaxespad=0.0, shadow=True, fancybox=True)

leg = plt.gca().get_legend()
ltext  = leg.get_texts()  # all the text.Text instance in the legend
llines = leg.get_lines()  # all the lines.Line2D instance in the legend
frame  = leg.get_frame()  # the patch.Rectangle instance surrounding the legend

# see text.Text, lines.Line2D, and patches.Rectangle for more info on
# the settable properties of lines, text, and rectangles
frame.set_facecolor('0.80')      # set the frame face color to light gray
plt.setp(ltext, fontsize='small', fontname='arial', family='monospace')    # the legend text fontsize
plt.setp(llines, linewidth=1.5)      # the legend linewidth


plt.show()
