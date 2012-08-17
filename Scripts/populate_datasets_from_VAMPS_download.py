# -*- coding: utf-8 -*-

import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../'))
import lib.fastalib as u

fasta = u.SequenceSource(sys.argv[1])
taxon = sys.argv[2]

output = u.FastaOutput(sys.argv[2].replace(';', ''))

while fasta.next():
    if fasta.id.find(taxon) > -1:
        acc = fasta.id.split('|')[0]
        project = fasta.id.split('|')[1].split('=')[1]
        dataset = fasta.id.split('|')[2].split('=')[1]
        new_id = project + '_' + dataset + '_' + acc

        abundance = int(fasta.id.split('|')[7].split('=')[1])

        for i in range(0, abundance):
            output.write_id('%s-%s|%s' % (new_id, str(i), fasta.id))
            output.write_seq(fasta.seq, split = False)

fasta.close()
output.close()
