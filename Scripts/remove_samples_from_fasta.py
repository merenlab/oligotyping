# removes samples from FASTA file:
#
# ./me FASTA_FILE sample_1,sample_2,[...],sample_N
#

import sys
import operator

import Oligotyping.lib.fastalib as u
from Oligotyping.utils.utils import pretty_print as pp

fasta = u.SequenceSource(sys.argv[1])
output = u.FastaOutput(sys.argv[1] + '-SAMPLES-REMOVED.fa')
samples_to_be_removed = sys.argv[2].split(',')

while fasta.next():
    if fasta.pos % 1000 == 0:
        sys.stderr.write('\rreads processed so far: %d' % (fasta.pos))
        sys.stderr.flush()
    sample_name = '_'.join(fasta.id.split('_')[:-1])

    if sample_name in samples_to_be_removed:
        continue

    output.store(fasta, split=False)

sys.stderr.write('\rNew FASTA file .............: %s\n' % (sys.argv[1] + '-SAMPLES-REMOVED.fa'))
fasta.close()
output.close()
