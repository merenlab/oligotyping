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
import tempfile
import cStringIO
import subprocess

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '..'))
import lib.fastalib as u 

# FIXME: carry this class into utils (and don't be quiet about the exception):
class LocalSearchResult:
    def __init__(self):
        self.query_label = None
        self.target_label = None
        self.percent_identity = None
        self.alignment_length = None
        self.number_of_mismatches = None
        self.number_of_gaps = None
        self.start_in_query = None
        self.end_in_query = None
        self.start_in_target = None 
        self.end_in_target = None
        self.query_length = None
        self.coverage = None
        self.query_aligned = None
        self.target_aligned = None
        self.hsp_match = None
         
    def init_with_b6output_line(self, line):
        self.b6output_line = line
        self.query_label, self.target_label, self.percent_identity,\
        self.alignment_length, self.number_of_mismatches, self.number_of_gaps,\
        self.start_in_query, self.end_in_query, self.start_in_target, self.end_in_target,\
        self.e_value, self.bit_score = line.strip().split('\t')

    def init_with_enhanced_line(self, line):
        self.enhanced_line = line
        self.query_label, self.target_label, self.percent_identity,\
        self.alignment_length, self.number_of_mismatches, self.number_of_gaps,\
        self.start_in_query, self.end_in_query, self.start_in_target, self.end_in_target,\
        self.e_value, self.bit_score, self.query_length, self.coverage,\
        self.query_aligned, self.target_aligned, self.hsp_match = line.strip().split('\t')
    
    def get_enhanced_results_line(self):
        enhanced_results_line = '%s\t%d\t%.2f\t%s\t%s\t%s\n' % (self.b6output_line, self.query_length, self.coverage,\
                                               self.query_aligned, self.target_aligned, self.hsp_match)
        return enhanced_results_line


# FIXME: this is not a good place to do this:
try:
    from cogent.align.algorithm import nw_align
except:
    print '''
    You need 'cogent' module in your Python path to run this software.

    You can get more information about the installation here:

        http://pycogent.sourceforge.net/install.html

'''
    sys.exit(-1)

try:
    from Bio.Blast import NCBIWWW
    from Bio.Blast import NCBIXML
except:
    print "This script requires BioPython (http://biopython.org/wiki/Biopython) to be installed."
    sys.exit()


def local_blast_search(sequence, reference_db, output_file = None):

    # FIXME: carry this function into utils..
    def get_read_by_read_id(ref_id):
        ref = u.SequenceSource(reference_db)
        while ref.next():
            if ref.id == ref_id:
                return ref.seq

    sequence = sequence.replace('-', '')
    query_tmp_fasta = tempfile.NamedTemporaryFile(delete=False)
    query_tmp_fasta_path = query_tmp_fasta.name
    results_tmp_path = tempfile.NamedTemporaryFile(delete=False).name
    
    query_tmp_fasta.write('>sequence\n%s\n' % sequence)
    query_tmp_fasta.close()

    usearch_process = ['usearch', '-query', query_tmp_fasta_path, '-db', reference_db, '-evalue', '0.0001', '-blast6out', results_tmp_path, '-maxaccepts', '20']
    if subprocess.call(usearch_process, stderr=open('/dev/null', 'w')) == 0:
        results = open(results_tmp_path).readlines()

        # enhance USEARCH results with more information
        updated_results = []
        for result in results:
            result = result.strip()

            # read all fields from USEARCH blast6output
            r = LocalSearchResult()
            r.init_with_b6output_line(result)

            query = sequence
            target = get_read_by_read_id(r.target_label)

            # parts that were aligned during the search are being aligned to each other to generate
            # hsp_match data to include into results
            r.query_aligned, r.target_aligned = nw_align(query[int(r.start_in_query) - 1:int(r.end_in_query)],\
                                                         target[int(r.start_in_target) - 1:int(r.end_in_target)]) 

            r.query_length = len(sequence)
            r.coverage = len(r.query_aligned.replace('-', '')) * 100.0 / r.query_length
            r.hsp_match = ''.join(['|' if r.query_aligned[i] == r.target_aligned[i] else ' ' for i in range(0, len(r.query_aligned))])
           
            updated_results.append(r.get_enhanced_results_line()) 
             
        if output_file:
            open(output_file, 'w').write(''.join(updated_results))

        os.remove(query_tmp_fasta_path)
        os.remove(results_tmp_path)

        return updated_results
    else:
        print 'Something went wrong while performing the local search :/'
        sys.exit(-2)


def get_local_blast_results_dict(blast_results, num_results = 20):

    blast_results_dict = {}

    for i in range(0, num_results if num_results < len(blast_results) else len(blast_results)):
        b = {}
        r = LocalSearchResult()
        r.init_with_enhanced_line(blast_results[i])

        b['query_length'] = int(r.query_length)
        b['hit_def'] = r.target_label
        b['accession'] = None
        b['ncbi_link'] = None
        b['hsp_query'] = r.query_aligned
        b['hsp_match'] = r.hsp_match
        b['hsp_subject'] = r.target_aligned
        b['identity'] = float(r.percent_identity)
        b['coverage'] = float(r.coverage)

        blast_results_dict[i] = b

    return blast_results_dict


def remote_blast_search(sequence, output_file = None):
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

    try:
        blast_results.close()
    except:
        pass

    return blast_results_dict


def blast_pretty_print(blast_results_dict, num_results = 1):
    num_results = len(blast_results_dict) if len(blast_results_dict) < num_results else num_results

    for i in range(0, num_results if num_results < len(blast_results_dict) else len(blast_results_dict)):
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

    parser = argparse.ArgumentParser(description='Query a sequence either in NCBI or search it against a local reference DB')
    parser.add_argument('sequence', metavar = 'SEQUENCE', help = 'Sequence to be queried')
    parser.add_argument('--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'File to store search results')
    parser.add_argument('--refdb', default = None, metavar = 'LOCAL_SEARCH',\
                        help = 'Perform a local search')

    args = parser.parse_args()

    if args.refdb == None:
        print 'performing online search'
        blast_results = remote_blast_search(args.sequence, args.output_file)
        blast_results_dict = get_blast_results_dict(blast_results)
        blast_pretty_print(blast_results_dict)
    else:
        print 'performing local search (db: %s)' % args.refdb
        blast_results = local_blast_search(args.sequence, args.refdb, args.output_file)
        blast_results_dict = get_local_blast_results_dict(blast_results)
        blast_pretty_print(blast_results_dict, 10)
