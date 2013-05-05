import sys
import operator

import Oligotyping.lib.fastalib as u
from Oligotyping.utils.utils import pretty_print as pp

fasta = u.SequenceSource(sys.argv[1])

samples = {}
while fasta.next():
    if fasta.pos % 1000 == 0:
        sys.stderr.write('\rreads processed so far: %d' % (fasta.pos))
        sys.stderr.flush()
    sample_name = '_'.join(fasta.id.split('_')[:-1])

    if samples.has_key(sample_name):
        samples[sample_name] += 1
    else:
        samples[sample_name] = 1

sys.stderr.write('\rSamples and read counts found in the FASTA file:\n')
for sample, read_count in sorted(samples.iteritems(), key=operator.itemgetter(1), reverse = True):
    print '%-30s %s' % (sample, pp(read_count)) 

print
print
print 'Total number of samples: ', pp(len(samples))
print 'Total number of reads: ', pp(fasta.pos)
print
fasta.close()
