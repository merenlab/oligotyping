import sys

import Oligotyping.lib.fastalib as u

fasta = u.SequenceSource(sys.argv[1])
output = u.FastaOutput(sys.argv[1] + '-PADDED-WITH-GAPS')

longest_read = 0
while fasta.next():
    if len(fasta.seq) > longest_read:
        longest_read = len(fasta.seq)

fasta.reset()

while fasta.next():
    if fasta.pos % 10000 == 0:
        sys.stdout.write('\rreads processed so far: %d' % (fasta.pos))
        sys.stdout.flush()

    gaps = longest_read - len(fasta.seq)
    
    output.write_id(fasta.id)
    output.write_seq(fasta.seq + '-' * gaps, split = False)


fasta.close()
print
