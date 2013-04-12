# -*- coding: utf-8 -*-

import sys

import Oligotyping.lib.fastalib as u

fasta = u.SequenceSource(sys.argv[1], lazy_init = False)
output = u.FastaOutput(sys.argv[1] + '-TRIMMED')

trim_from = int(sys.argv[2])
trim_to = int(sys.argv[3]) if len(sys.argv) == 4 else None

while fasta.next():
    output.write_id(fasta.id)
    if trim_to:
        output.write_seq(fasta.seq[trim_from:trim_to], split = False)
    else:
        output.write_seq(fasta.seq[trim_from:], split = False)

fasta.close()
output.close()
