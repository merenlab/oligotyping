import os
import sys

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../lib'))
import fastalib as u

taxon = sys.argv[1]

PATH_TO_OTU_ID_TO_GREENGENES_FILE = '/Users/meren/Desktop/MBL/Oligotyping/GreenGenes/otu_id_to_greengenes.txt'
PATH_TO_GREENGENES_ALIGNMENT_FILE = '/Users/meren/Desktop/MBL/Oligotyping/GreenGenes/gg_97_otus_6oct2010_aligned.fasta'

ids = []

for id, tax in [line.strip().split('\t') for line in open(PATH_TO_OTU_ID_TO_GREENGENES_FILE).readlines()]:
    if tax.find(taxon) > 0:
        ids.append(id)

ids = list(set(ids))
print '%d ids found for %s.' % (len(ids), taxon)

template = u.FastaOutput('%s.tmpl' % taxon)
fasta = u.SequenceSource(PATH_TO_GREENGENES_ALIGNMENT_FILE)
while fasta.next():
    if fasta.id in ids:
        template.store(fasta, split = False)
        ids.remove(fasta.id)

fasta.close()
template.close()
