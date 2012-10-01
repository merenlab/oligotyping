#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import sys
import argparse

import Oligotyping.lib.fastalib as u


def main(fasta_file_path, min_percent = 95.0, output_file_path = None):
    fasta = u.SequenceSource(fasta_file_path)
    
    fasta.next()
    alignment_length = len(fasta.seq)
    fasta.reset()
    
    positions = {}
    
    while fasta.next():
        if fasta.pos % 1000 == 0:
            sys.stderr.write('\rAnalyzing all reads; pos: %d' % fasta.pos)
            sys.stderr.flush()
        for i in range(0, alignment_length):
            if fasta.seq[i] != '-':
                for j in range(i, alignment_length):
                    try:
                        positions[j] += 1
                    except:
                        positions[j] = 1
                break
    
    fasta.reset()
    sys.stderr.write('\n')
    
    num_reads = positions[alignment_length - 1]
    trim_location = 0
    
    for i in range(0, alignment_length):
        pct_reads_will_survive = positions[i] * 100.0 / num_reads
        if pct_reads_will_survive >= min_percent and not trim_location:
            trim_location = i
            trim_location_pct_reads_survive = pct_reads_will_survive
        if pct_reads_will_survive == 100:
            print
            print 'All reads are going to be trimmed from the %dth position.' % (trim_location_pct_reads_survive)
            
            if 100 - trim_location_pct_reads_survive:
                print
                print '%d reads that do not reach to this locaition will be eliminated.' % ((100 - trim_location_pct_reads_survive) / 100.0 * num_reads)
            
            if min_percent < 100:
                print
                print 'If all reads were to be retained, alignments should have been trimmed from'
                print 'the %dth location, however, this would have required all reads to lose %d' % (i, i - trim_location)
                print 'bases'
            print
            break


    output = u.FastaOutput(output_file_path if output_file_path else sys.argv[1] + '-TRIMMED')

    while fasta.next():
        if fasta.pos % 1000 == 0:
            sys.stderr.write('\rStoring trimmed reads; pos: %d' % fasta.pos)
            sys.stderr.flush()

        if fasta.seq[trim_location:].startswith('-'):
            continue
        else:
            output.write_id(fasta.id)
            output.write_seq(fasta.seq[trim_location:], split = False)

    sys.stderr.write('\n')
    sys.stderr.write('\n')
    print 'Trimmed reads stored: "%s"\n' % (output_file_path if output_file_path else sys.argv[1] + '-TRIMMED')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Smart trim ragged ends from the beginning of an alignment')
    parser.add_argument('fasta_file', metavar = 'FASTA FILE',
                        help = 'Alignment to be trimmed')
    parser.add_argument('--min-percent', metavar = 'PERCENT', default=95.0, type=float,
                        help = 'Even if there is only one read that is too short and therefore full of gap characters,\
                                the first location in the alignment file that *every* read has a base would have to match\
                                the length of that short read. With this percentage you can specify what is the percentage\
                                of reads you expect to pass while this trimming script tries to maximize the remaining\
                                read length after trimming. Default is %(default).2f')
    parser.add_argument('-o', '--output', help = 'Output file name', default = None)


    args = parser.parse_args()
    main(args.fasta_file, args.min_percent, args.output)
