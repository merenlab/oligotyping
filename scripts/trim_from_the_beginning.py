# -*- coding: utf-8 -*-

import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../lib'))
import fastalib as u

fasta = u.SequenceSource(sys.argv[1], lazy_init = False)
offset = int(sys.argv[2])
output = u.FastaOutput(sys.argv[1] + '-TRIMMED_BEGINNING')

while fasta.next():
    output.write_id(fasta.id)
    output.write_seq(fasta.seq[offset:], split = False)

fasta.close()
output.close()
