#!/usr/bin/env python

# get me some alignments and qual scores:
#
# python $me.py fasta.fa fasta.qual
#

import sys
import cPickle
from scipy import log2 as log

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import Oligotyping.lib.fastalib as u

from Oligotyping.utils.utils import get_quals_dict
from Oligotyping.utils.random_colors import get_list_of_colors

alignment_file = sys.argv[1]
quals_file = sys.argv[2]

quals_dict = get_quals_dict(quals_file, alignment_file)

cPickle.dump(quals_dict, open(alignment_file + '-QUALS-DICT', 'w'))

print 'output file:'
print '  - "%s"' % (alignment_file + '-QUALS-DICT')
