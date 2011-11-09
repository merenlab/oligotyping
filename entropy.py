#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import operator
from scipy import log2 as log
import numpy as np

def entropy(l):
    p_a = (len([x for x in l if x.upper() == 'A']) * 1.0 / len(l)) + 0.0000000000000000001
    p_t = (len([x for x in l if x.upper() == 'T']) * 1.0 / len(l)) + 0.0000000000000000001
    p_c = (len([x for x in l if x.upper() == 'C']) * 1.0 / len(l)) + 0.0000000000000000001
    p_g = (len([x for x in l if x.upper() == 'G']) * 1.0 / len(l)) + 0.0000000000000000001
    p_n = (len([x for x in l if x.upper() == '-']) * 1.0 / len(l)) + 0.0000000000000000001

    return -(p_a*log(p_a) + p_t*log(p_t) + p_c*log(p_c) + p_g*log(p_g) + p_n*log(p_n))

lines = [l for l in open(sys.argv[1]) if not l.startswith('>')]

entropy_tpls = []

start = 0
end = len(lines[0])

for i in range(start, end):
    sys.stderr.write('\rENTROPY: %d/%d' % (i, end))
    sys.stderr.flush()

    if set([x[i] for x in lines]) == set(['.']) or set([x[i] for x in lines]) == set(['-']):
        entropy_tpls.append((i, 0.0),)
    else:
        column = "".join([x[i] for x in lines])
        e = entropy(column)
        if e < 0.00001:
            entropy_tpls.append((i, 0.0),)
        else:
            entropy_tpls.append((i, e),)
            #print i, e
print

stuff = [x[1] for x in entropy_tpls]

def draw(stuff):
    import matplotlib.pyplot as plt
    import numpy as np

    fig = plt.figure(figsize = (10, 2))

    plt.rcParams.update({'axes.linewidth' : 0.1})
    plt.rc('grid', color='0.70', linestyle='-', linewidth=0.1)
    plt.grid(True)

    plt.subplots_adjust(hspace = 0, wspace = 0, right = 0.995, left = 0.050, top = 0.92, bottom = 0.10)

    ax = fig.add_subplot(111)
    ind = np.arange(len(stuff))
    ax.bar(ind, stuff, color = 'black', lw = 0.5)
    ax.set_xlim([0, end])
    plt.savefig(sys.argv[1] + '-ENTROPY.png')
    plt.show()

e = sorted(entropy_tpls, key=operator.itemgetter(1), reverse=True)

entropy_output = open(sys.argv[1] + '-ENTROPY.txt', 'w')

for _component, _entropy in sorted(entropy_tpls, key=operator.itemgetter(1), reverse=True):
    entropy_output.write('%d\t%.4f\n' % (_component, _entropy))

entropy_output.close()

draw(stuff)
