# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import os
import sys
import cStringIO

try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except:
    print "This script requires BioPython (http://biopython.org/wiki/Biopython) to be installed."
    sys.exit()

def blast(sequence, output_file = None):
    result_handle = NCBIWWW.qblast("blastn", "nt", sequence)

    result = result_handle.read()

    if output_file:
       open(output_file, "w").write(result)
 
    return cStringIO.StringIO(result)

def get_blast_results_dict(blast_results, num_results = 5):
    blast_results_dict = {}

    blast_record = list(NCBIXML.parse(blast_results))[0]
    num_results = len(blast_record.alignments) if len(blast_record.alignments) < num_results else num_results

    for i in range(0, num_results):
        b = {}
        b['query_length'] = int(blast_record.query_length)
        
        alignment = blast_record.alignments[i]
        hsp = alignment.hsps[0]
        
        b['hit_def'] = alignment.hit_def
        b['accession'] = alignment.accession
        b['ncbi_link'] = 'http://www.ncbi.nlm.nih.gov/nuccore/%s' % b['accession']
        b['hsp_query'] = hsp.query
        b['hsp_match'] = hsp.match
        b['hsp_subject'] = hsp.sbjct

        b['identity'] = len([x for x in hsp.match if x == '|']) * 100.0 / len(hsp.query)
        b['coverage'] = len(hsp.query) * 100.0 / b['query_length']

        blast_results_dict[i] = b

    return blast_results_dict

def blast_pretty_print(blast_results_dict, num_results = 1):
    num_results = len(blast_results_dict) if len(blast_results_dict) < num_results else num_results

    for i in range(0, num_results):
        b = blast_results_dict[i]
        print
        print '#######################################################################################'
        print b['hit_def']
        print '#######################################################################################'
        print
        print 'Accession: %s (%s)' % (b['accession'], b['ncbi_link'])
        print 'Coverage : %.2f' % b['coverage']
        print 'Identity : %.2f' % b['identity']
        print
        for j in range(0, len(b['hsp_query']) - 80, 80):
            print b['hsp_query'][j:j+80]
            print b['hsp_match'][j:j+80]
            print b['hsp_subject'][j:j+80]
            print
        print b['hsp_query'][j+80:]
        print b['hsp_match'][j+80:]
        print b['hsp_subject'][j+80:]
            
        

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Query a sequence in NCBI')
    parser.add_argument('sequence', metavar = 'SEQUENCE', help = 'Sequences to be queried in NCBI')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File to store XML output')

    args = parser.parse_args()

    blast_results = blast(args.sequence, args.output_file)
    blast_results_dict = get_blast_results_dict(blast_results)
    blast_pretty_print(blast_results_dict)
