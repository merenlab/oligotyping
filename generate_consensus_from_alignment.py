# -*- coding: utf-8 -*-

import sys
import os
import operator

sys.path.append('..')
import fastalib as u

consensus_sequence = ''

fasta = u.SequenceSource(sys.argv[1])
fasta.next()
alignment_length = len(fasta.seq)
consensus_dict = {}
for i in range(0, alignment_length):
    consensus_dict[i] = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '-': 0}

fasta.reset()

while fasta.next():
    for pos in range(0, alignment_length):
        consensus_dict[pos][fasta.seq[pos]] += 1

for pos in range(0, alignment_length):
    consensus_sequence += sorted(consensus_dict[pos].iteritems(), key=operator.itemgetter(1), reverse=True)[0][0]

print '>' + os.path.basename(sys.argv[1])
print consensus_sequence
