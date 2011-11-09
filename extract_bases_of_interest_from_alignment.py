#!/usr/bin/python
# -*- coding: utf-8
# ~-

import os
import sys
sys.path.append('..')
import copy
import cPickle
import operator

import fastalib as u

# <EDIT> ##########################################################################################################
NUMBER_OF_COMPONENTS_TO_USE                               = 5
MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS                   = 5
MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE = 3
SAMPLE_FROM_DEFLINE                                       = lambda x: "_".join(x.split('|')[0].split("_")[0:-1])
# </EDIT> #########################################################################################################


OUTPUT_DIR         = os.path.dirname(sys.argv[1]) or '.'
OUTPUT_FILE_PREFIX = '%s-C%d-S%d-A%d' % (os.path.basename(sys.argv[1]).split('.')[0],
                                         NUMBER_OF_COMPONENTS_TO_USE,
                                         MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS,
                                         MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE)

GEN_OUTPUT_DEST = lambda postfix: os.path.join(OUTPUT_DIR, OUTPUT_FILE_PREFIX + '-' + postfix)

info_file_path = GEN_OUTPUT_DEST('RUNINFO')
info_file_obj = open(info_file_path, 'w')

def info(label, value, file_obj = None):
    info_line = "%s %s: %s" % (label, '.' * (60 - len(label)), str(value))
    if file_obj:
        info_file_obj.write(info_line + '\n')
    print info_line

def pp(n):
    """Pretty print function for very big numbers.."""
    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

# ~

fasta = u.SequenceSource(sys.argv[1])
column_entropy = [int(x.strip().split()[0]) for x in open(sys.argv[2]).readlines()]

info('Output directory', OUTPUT_DIR)
info('Extraction info output file', info_file_path, info_file_obj)
info('Input FASTA file', sys.argv[1], info_file_obj)
info('Input entropy file', sys.argv[2], info_file_obj)
info('Number of sequences in FASTA', pp(fasta.total_seq), info_file_obj)
info('Number of entropy components to use', NUMBER_OF_COMPONENTS_TO_USE, info_file_obj)
info('Min number of samples oligotype appears', MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS, info_file_obj)
info('Min % abundance of oligotype in at least one sample', MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE, info_file_obj)

# locations of interest based on the entropy scores
bases_of_interest_locs = sorted([column_entropy[i] for i in range(0, NUMBER_OF_COMPONENTS_TO_USE)])
info('Bases of interest', ', '.join([str(x) for x in bases_of_interest_locs]), info_file_obj)

# constructing samples_dict
samples_dict = {}
samples = []
while fasta.next():
    sample = SAMPLE_FROM_DEFLINE(fasta.id)
    
    if not samples_dict.has_key(sample):
        samples_dict[sample] = {}
        samples.append(sample)

    oligo = ''.join(fasta.seq[o] for o in bases_of_interest_locs)

    if samples_dict[sample].has_key(oligo):
        samples_dict[sample][oligo] += 1
    else:
        samples_dict[sample][oligo] = 1
info('Number of samples in FASTA', pp(len(samples_dict)), info_file_obj)

# cat oligos | uniq
oligos_set = []
for sample in samples:
    for oligo in samples_dict[sample].keys():
        if oligo not in oligos_set:
            oligos_set.append(oligo)
info('Number of unique oligotypes', pp(len(oligos_set)), info_file_obj)

# count oligo abundance
oligo_abundance = []
for oligo in oligos_set:
    count = 0
    for sample in samples:
        if oligo in samples_dict[sample].keys():
            count += 1
    oligo_abundance.append((count, oligo),)
oligo_abundance.sort()

# eliminate singleton/doubleton oligos (any oligo required to appear in at least
# MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS samples)
non_singleton_oligos = []
for tpl in oligo_abundance:
    if tpl[0] >= MIN_NUMBER_OF_SAMPLES_OLIGOTYPE_APPEARS:
        non_singleton_oligos.append(tpl[1])
info('Oligotypes after "min number of samples" elimination', pp(len(non_singleton_oligos)), info_file_obj)

# eliminate very rare oligos (the percent abundance of every oligo should be
# more than MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE percent
# in at least one sample)
abundant_oligos = []
for oligo in non_singleton_oligos:
    percent_abundances = []
    for sample in samples:
        if samples_dict[sample].has_key(oligo):
            percent_abundances.append(samples_dict[sample][oligo] * 100.0 / sum([samples_dict[sample][o] for o in non_singleton_oligos if samples_dict[sample].has_key(o)]))
    percent_abundances.sort()
    if percent_abundances[-1] >= MIN_PERCENT_ABUNDANCE_OF_OLIGOTYPE_IN_AT_LEAST_ONE_SAMPLE:
        abundant_oligos.append(oligo)
info('Oligotypes after "min % abundance in a sample" elimination', pp(len(abundant_oligos)), info_file_obj)

# store abundant oligos
abundant_oligos_file_path = GEN_OUTPUT_DEST("OLIGOS.fasta")
f = open(abundant_oligos_file_path, 'w')
for oligo in sorted(abundant_oligos):
    f.write('>' + oligo + '\n')
    f.write(oligo + '\n')
f.close()
info('Abundant oligotypes file path', abundant_oligos_file_path, info_file_obj)

# generate NEXUS file of oligos
oligos_nexus_file_path = GEN_OUTPUT_DEST("OLIGOS.nexus")
f = open(oligos_nexus_file_path, 'w')
f.write("""begin data;
    dimensions ntax=%d nchar=%d;
    format datatype=dna interleave=no;
    matrix\n""" % (len(abundant_oligos), len(abundant_oligos[0])))
for oligo in sorted(abundant_oligos):
    f.write('    %.20s %s\n' % (oligo, oligo))
f.write('    ;\n')
f.write('end;\n')
f.close()
info('NEXUS file for oligotypes', oligos_nexus_file_path, info_file_obj)

# generate environment file
environment_file_path = GEN_OUTPUT_DEST("ENVIRONMENT.txt")
f = open(environment_file_path, 'w')
for sample in samples:
    for oligo in samples_dict[sample]:
        if oligo in abundant_oligos and samples_dict[sample][oligo] > 0:
            f.write("%s\t%s\t%d\n" % (oligo, sample, samples_dict[sample][oligo]))
f.close()
info('Environment file for Viamics/UniFrac analysis', environment_file_path, info_file_obj)

# generate viamics samples dict 
viamics_samples_dict = {}
viamics_samples_dict_file_path = GEN_OUTPUT_DEST("ENVIRONMENT.cPickle")
for sample in samples_dict:
    viamics_samples_dict[sample] = {}
    viamics_samples_dict[sample]['species'] = {}
    for oligo in samples_dict[sample]:
        if oligo in abundant_oligos:
            viamics_samples_dict[sample]['species'][oligo] = samples_dict[sample][oligo]

for sample in viamics_samples_dict:
    viamics_samples_dict[sample]['tr'] = sum(viamics_samples_dict[sample]['species'].values())
    viamics_samples_dict[sample]['bases_of_interest_locs'] = bases_of_interest_locs
cPickle.dump(viamics_samples_dict, open(viamics_samples_dict_file_path, 'w'))
info('Serialized Viamics samples dictionary', viamics_samples_dict_file_path, info_file_obj)

# create a fasta file with a representative full length consensus sequence for every oligotype
representative_oligotypes_file_path = GEN_OUTPUT_DEST("REPRESENTATIVE-OLIGO-SEQS.fasta")
f = open(representative_oligotypes_file_path, 'w')
for abundant_oligo in abundant_oligos:
    fasta.reset()
    counter = 0
    rep_sequences = []
    while fasta.next() and counter < 5000:
        oligo = ''.join(fasta.seq[o] for o in bases_of_interest_locs)
        if oligo == abundant_oligo:
            rep_sequences.append(fasta.seq)
            counter += 1
    
    consensus_dict = {}
    consensus_sequence = ''
    alignment_length = len(rep_sequences[0])
    for i in range(0, alignment_length):
        consensus_dict[i] = {'A': 0, 'T': 0, 'C': 0, 'G': 0, '-': 0}
    for sequence in rep_sequences:
        for pos in range(0, alignment_length):
            consensus_dict[pos][sequence[pos]] += 1
    for pos in range(0, alignment_length):
        consensus_sequence += sorted(consensus_dict[pos].iteritems(), key=operator.itemgetter(1), reverse=True)[0][0]
    f.write('>' + abundant_oligo + '\n')
    f.write(consensus_sequence + '\n')
f.close()
os.system('python trim_uninformative_columns_from_alignment.py %s' % representative_oligotypes_file_path)
os.system('mv %s %s' % (representative_oligotypes_file_path + '-TRIMMED', representative_oligotypes_file_path))
info('Representative sequences for oligotypes', representative_oligotypes_file_path, info_file_obj) 

info_file_obj.close()
# done.
