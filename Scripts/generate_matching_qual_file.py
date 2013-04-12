# -*- coding: utf-8 -*-

import sys

import Oligotyping.lib.fastalib as u

alignment = u.SequenceSource(sys.argv[1])
quals = u.SequenceSource(sys.argv[2])

alignment.next()
quals.next()

qual = [int(q) for q in quals.seq.split()]
qual_aligned = []
for i in range(0, len(alignment.seq)):
    if alignment.seq[i] != '-':
        qual_aligned.append(qual.pop(0))
    else:
        qual_aligned.append(None)
print alignment.seq
print qual_aligned
