#!/usr/bin/env python

# get me some alignments and qual scores:
#
# python $me.py fasta.fa fasta.qual
#

import sys
import cPickle

from Oligotyping.utils.utils import get_quals_dict

alignment_file = sys.argv[1]
quals_file = sys.argv[2]

quals_dict = get_quals_dict(quals_file, alignment_file)

cPickle.dump(quals_dict, open(alignment_file + '-QUALS-DICT', 'w'))

print 'output file:'
print '  - "%s"' % (alignment_file + '-QUALS-DICT')
