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
import time
import math
import fcntl
import shutil
import struct
import hashlib
import random
import string
import termios 
import cPickle
import tempfile
import subprocess
import numpy as np
import multiprocessing

from Oligotyping.lib import fastalib as u
from Oligotyping.utils.constants import pretty_names
from Oligotyping.utils.aligner import nw_align

P = lambda x, y: '%.2f%%' % (x * 100.0 / y)

NUCL_COLORS = {'A': 'red', 'T': 'blue',  'C': 'green', 'G': 'purple', 
               'N': 'white', '-': '#CACACA', 'K': '#CACACA', 'R': '#CACACA',
               'Y': '#CACACA', 'W': '#CACACA', 'S': '#CACACA', 'M': '#CACACA'}

class ConfigError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Config Error: %s' % self.e


class LibError(Exception):
    def __init__(self, e = None):
        Exception.__init__(self)
        self.e = e
        return
    def __str__(self):
        return 'Lib Error: %s' % self.e


def get_unit_counts_and_percents(units, samples_dict):
    # this function returns two dictionaries that contain unit counts and percents
    # across samples. that can be used for agglomeration, as well as generation of
    # environment and matrix files.

    unit_percents = {}
    unit_counts = {}

    sample_totals = {}
    for sample in samples_dict:
        sample_totals[sample] = sum(samples_dict[sample].values())

    for sample in samples_dict:
        counts = []
        percents = []
        for unit in units:
            if samples_dict[sample].has_key(unit):
                counts.append(samples_dict[sample][unit])
                percents.append(samples_dict[sample][unit] * 100.0 / sample_totals[sample])
            else:
                counts.append(0)
                percents.append(0.0)
                
        unit_counts[sample] = counts
        unit_percents[sample] = percents

    return (unit_counts, unit_percents)


def import_error(e):
    print '''
    Sorry. It seems you are missing a module that is required by
    the oligotyping pipeline. Here is the original import error:

        "%s"

    Please make sure you have all the requirements installed. This
    post may be useful to set up the right environment for the
    oligotyping pipeline:

    http://oligotyping.org/2012/05/11/oligotyping-pipeline-explained/

    \n''' % e
    sys.exit()


def generate_MATRIX_files(units, samples, unit_counts, unit_percents, matrix_count_file_path, matrix_percent_file_path):
    count_file = open(matrix_count_file_path, 'w')
    percent_file = open(matrix_percent_file_path, 'w')       
    
    count_file.write('\t'.join(['samples'] + units) + '\n')
    percent_file.write('\t'.join(['samples'] + units) + '\n')

    for sample in samples:
        count_file.write('\t'.join([sample] + [str(c) for c in unit_counts[sample]]) + '\n')
        percent_file.write('\t'.join([sample] + [str(p) for p in unit_percents[sample]]) + '\n')
    
    count_file.close()
    percent_file.close()


def get_units_across_samples_dicts(units, samples, unit_percents):
    across_samples_sum_normalized = {}
    across_samples_max_normalized = {}
  
    for unit in units:
        across_samples_sum_normalized[unit] = []
        across_samples_max_normalized[unit] = []

    for i in range(0, len(units)):
        unit = units[i]
        sum_across_samples = sum([unit_percents[sample][i] for sample in samples])
        max_across_samples = max([unit_percents[sample][i] for sample in samples])
            
        for sample in samples:
            across_samples_sum_normalized[unit].append(unit_percents[sample][i]  * 100.0 / sum_across_samples)
            across_samples_max_normalized[unit].append(unit_percents[sample][i]  * 100.0 / max_across_samples)

    return (across_samples_sum_normalized, across_samples_max_normalized)


def generate_MATRIX_files_for_units_across_samples(units, samples, MN_fp, SN_fp, MN_dict, SN_dict):
        across_samples_MN_file = open(MN_fp, 'w')
        across_samples_SN_file = open(SN_fp, 'w')

        across_samples_MN_file.write('\t'.join(['sample'] + units) + '\n')
        across_samples_SN_file.write('\t'.join(['sample'] + units) + '\n')

        for i in range(0, len(samples)):
            sample = samples[i]
            across_samples_MN_file.write('\t'.join([sample] + [str(MN_dict[unit][i]) for unit in units]) + '\n')
            across_samples_SN_file.write('\t'.join([sample] + [str(SN_dict[unit][i]) for unit in units]) + '\n')
        
        across_samples_MN_file.close()
        across_samples_SN_file.close()


def get_num_nt_diff_between_two_aligned_sequences(seq1, seq2):
    if len(seq1) != len(seq2):
        raise LibError, "Two sequences are not equal in length:\n\t%s\n\t%s" % (seq1, seq2)

    return len(["diff" for i in range(0, len(seq1)) if seq1[i] != seq2[i]])


def homopolymer_indel_exists(seq1, seq2):
    seq1, seq2 = trim_uninformative_gaps_from_sequences(seq1, seq2)
    
    # sometimes alignments look like this:
    #
    #    CCCGAAAAAA--TAT
    #    CCCGAAA---AATAT
    #
    # where the correct alignment should look like this:
    #
    #    CCCGAAAAAATAT
    #    CCCGAAAAA-TAT
    # 
    # causes this function to return false. in order to fix that problem
    # we perform needleman-wunch alignment here:
    if sum([seq1.count('-'), seq2.count('-')]) > 1:
        seq1, seq2 = nw_align(seq1.replace('-', ''), seq2.replace('-', ''))

    gap_index = seq1.find('-')
    if gap_index == -1:
        gap_index = seq2.find('-')
        
        # so the gap is in seq2. replace seq1 and 2 so it would be certain
        # that the sequence with gap is seq1:
        seq1, seq2 = seq2, seq1
        
    if gap_index == -1:
        return False

    isHP = lambda x: len(set(x)) == 1
    isHPindel = lambda (s, e): seq1[s:e] == seq2[s:e] and isHP(seq1[s:e]) == 1 and seq2[gap_index] == seq2[s]
    
    def DownStream(sequence):
        i = 3
        while isHP(sequence[gap_index - i - 1:gap_index]):
            i += 1
        return (gap_index - i, gap_index)

    def UpStream(sequence):
        i = 4
        while isHP(sequence[gap_index + 1:gap_index + i + 1]):
            i += 1
        return (gap_index + 1, gap_index + i)

    # check downstream of the gap
    if gap_index >= 3:
        if isHPindel(DownStream(seq1)):
            return True
        
    # check upstream of the gap
    if len(seq1) - gap_index > 3:
        if isHPindel(UpStream(seq1)):
            return True
        
    return None


def append_file(target_path, source_path, remove_source = True):
    target = open(target_path, 'a')
    source = open(source_path, 'r')
    
    for line in source.readlines():
        target.write(line)
    
    target.close()
    source.close()
    
    if remove_source:
        os.remove(source_path)


def append_reads_to_FASTA(read_id_sequence_tuples_list, fasta_file_path):
    fasta_file = open(fasta_file_path, 'a')
    
    for read_id, sequence in read_id_sequence_tuples_list:
        fasta_file.write('>%s\n' % read_id)
        fasta_file.write('%s\n' % sequence)
    
    fasta_file.close()


def remove_white_space_mask_from_B6_entry(entry, defline_white_space_mask = '<$!$>'):
    # a stupid workaround due to the stupid behavior of blastn. if there are white spaces in the
    # defline of a FASTA file, it removes anything after the first white space. so if any input file
    # underwent of any cleaning process, this is an attempt to put those white spaces back in their
    # place.
    entry.subject_id = entry.subject_id.replace(defline_white_space_mask, ' ')
    entry.hit_def = entry.hit_def.replace(defline_white_space_mask, ' ')
    entry.query_id = entry.query_id.replace(defline_white_space_mask, ' ')
    entry.raw_line = entry.raw_line.replace(defline_white_space_mask, ' ')

    return entry

    

def mask_defline_whitespaces_in_FASTA(fasta_file_path, defline_white_space_mask = '<$!$>'):
    temp_file_path = fasta_file_path + '.tmp'
    fasta = u.SequenceSource(fasta_file_path)
    output = u.FastaOutput(fasta_file_path + '.tmp')
    
    while fasta.next():
        output.write_id(fasta.id.replace(' ', defline_white_space_mask))
        output.write_seq(fasta.seq, split = False)

    shutil.move(temp_file_path, fasta_file_path)

def unique_and_store_alignment(alignment_path, output_path):
    output = u.FastaOutput(output_path)
    alignment = u.SequenceSource(alignment_path, unique = True)
        
    alignment.next()
    most_abundant_unique_read = alignment.seq
    alignment.reset()
 
    read_ids = []
    unique_read_counts = []
    while alignment.next():
        read_ids += alignment.ids
        unique_read_counts.append(len(alignment.ids))
        output.store(alignment, split = False)
            
    output.close()
    alignment.close()
        
    return (read_ids, unique_read_counts, most_abundant_unique_read)


def generate_TAB_delim_file_from_dict(data_dict, output_file_path, order, first_column = 'samples'):
    f = open(output_file_path, 'w')
    
    f.write('%s\n' % ('\t'.join([first_column] + order)))
    for item in data_dict:
        line = [item]
        for column in order:
            if not data_dict[item].has_key(column):
                line.append('')
            else:
                line.append(str(data_dict[item][column]))
        f.write('%s\n' % '\t'.join(line))

    f.close()


def generate_ENVIRONMENT_file(samples, samples_dict, environment_file_path):
    # generate environment file
    f = open(environment_file_path, 'w')
    for sample in samples:
        for unit in samples_dict[sample]:
            f.write("%s\t%s\t%d\n" % (unit, sample, samples_dict[sample][unit]))
    f.close()


def get_unique_sequences_from_FASTA(alignment, limit = 10):
    unique_sequences = []

    fasta = u.SequenceSource(alignment, unique = True, lazy_init = False)

    while fasta.next() and fasta.pos < limit:
        unique_sequences.append((fasta.seq, len(fasta.ids), len(fasta.ids) / float(fasta.total_seq)))

    return unique_sequences



def get_purity_score(final_oligo_unique_distribution_dict):
    for oligo in final_oligo_unique_distribution_dict:
        freq_dict = final_oligo_unique_distribution_dict[oligo]
        if len(final_oligo_unique_distribution_dict[oligo]) > 1:
            bp = (freq_dict[1] / (freq_dict[0] * 1.0))
            final_purity_score_dict[oligo] = 1 - bp
        else:
            final_purity_score_dict[oligo] = 1.00

def get_total_purity_score(final_purity_score_dict):
    sorted_scores = sorted(final_purity_score_dict.values())
    last_quarter = sorted_scores[:int(math.ceil(len(sorted_scores)/4.0))] # take the last quarter of the unique sequences   
    final_total = reduce(lambda x, y: x + y, last_quarter) / len(last_quarter) # take the average of these sequences 
        
    total_purity_score_dict =  "%.2f" %final_total

def get_oligos_sorted_by_abundance(samples_dict, oligos = None, min_abundance = 0):
    samples = samples_dict.keys()
    samples.sort()

    if oligos == None:
        oligos = []
        map(lambda o: oligos.extend(o), [v.keys() for v in samples_dict.values()])
        oligos = list(set(oligos))

    abundant_oligos = []
    
    for oligo in oligos:
        percent_abundances = []

        for sample in samples:
            sum_sample = sum(samples_dict[sample].values())
            if samples_dict[sample].has_key(oligo):
                percent_abundances.append((samples_dict[sample][oligo] * 100.0 / sum_sample,\
                                           samples_dict[sample][oligo], sum_sample, sample))

        percent_abundances.sort(reverse = True)

        # FIXME: excuse me, WTF is going on here?:
        for abundance_percent, abundance_count, sample_size, sample in percent_abundances:
            abundant_oligos.append((sum([x[1] for x in percent_abundances]), oligo))
            break

    return [x[1] for x in sorted(abundant_oligos) if x[0] > min_abundance]


def get_vectors_from_oligotypes_across_samples_matrix(file_path):
    oligotypes_across_samples_file_obj = open(file_path)
   
    oligos = []
    vectors = {}

    for line in oligotypes_across_samples_file_obj.readlines()[1:]:
        fields = line.strip().split('\t')
        
        oligo = fields[0]
        oligos.append(oligo)
        vectors[oligo] = [float(c) for c in fields[1:]]

    return (oligos, vectors)


def generate_gexf_network_file(units, samples_dict, unit_percents, output_file, sample_mapping_dict = None,
                               unit_mapping_dict = None, project = None, sample_size=8, unit_size=2,
                               skip_sample_labels = False, skip_unit_labels = False):
    output = open(output_file, 'w')

    samples = sorted(samples_dict.keys())
    sample_mapping_categories = sorted([k for k in sample_mapping_dict.keys() if k != 'colors']) if sample_mapping_dict else None
    unit_mapping_categories = sorted([k for k in unit_mapping_dict.keys() if k not in ['colors', 'labels']]) if unit_mapping_dict else None
    
    output.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
    output.write('''<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" xmlns="http://www.gexf.net/1.2draft" version="1.2">\n''')
    output.write('''<meta lastmodifieddate="2010-01-01+23:42">\n''')
    output.write('''    <creator>Oligotyping pipeline</creator>\n''')
    if project:
        output.write('''    <creator>Network description for %s</creator>\n''' % (project))
    output.write('''</meta>\n''')
    output.write('''<graph type="static" defaultedgetype="undirected">\n\n''')

    if sample_mapping_dict:
        output.write('''<attributes class="node" type="static">\n''')
        for i in range(0, len(sample_mapping_categories)):
            category = sample_mapping_categories[i]
            output.write('''    <attribute id="%d" title="%s" type="string" />\n''' % (i, category))
        output.write('''</attributes>\n\n''')

    # FIXME: IDK what the hell is this one about:
    if unit_mapping_dict:
        output.write('''<attributes class="edge">\n''')
        for i in range(0, len(unit_mapping_categories)):
            category = unit_mapping_categories[i]
            output.write('''    <attribute id="%d" title="%s" type="string" />\n''' % (i, category))
        output.write('''</attributes>\n\n''')

    output.write('''<nodes>\n''')
    for sample in samples:
        if skip_sample_labels:
            output.write('''    <node id="%s">\n''' % (sample))
        else:
            output.write('''    <node id="%s" label="%s">\n''' % (sample, sample))
        output.write('''        <viz:size value="%d"/>\n''' % sample_size)
        if sample_mapping_dict and sample_mapping_dict.has_key('colors'):
            output.write('''        <viz:color r="%d" g="%d" b="%d" a="1"/>\n''' %\
                                             HTMLColorToRGB(sample_mapping_dict['colors'][sample], scaled = False))

        if sample_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(sample_mapping_categories)):
                category = sample_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, sample_mapping_dict[category][sample]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    for unit in units:
        if skip_unit_labels:
            output.write('''    <node id="%s">\n''' % (unit))
        else:
            if unit_mapping_dict and unit_mapping_dict.has_key('labels'):
                output.write('''    <node id="%s" label="%s">\n''' % (unit, unit_mapping_dict['labels'][unit]))
            else:
                output.write('''    <node id="%s">\n''' % (unit))
        output.write('''        <viz:size value="%d" />\n''' % unit_size)

        if unit_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(unit_mapping_categories)):
                category = unit_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    output.write('''</nodes>\n''')
    
    edge_id = 0
    output.write('''<edges>\n''')
    for sample in samples:
        for i in range(0, len(units)):
            unit = units[i]
            if unit_percents[sample][i] > 0.0:
                if unit_mapping_dict:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f">\n''' % (edge_id, unit, sample, unit_percents[sample][i]))
                    output.write('''        <attvalues>\n''')
                    for i in range(0, len(unit_mapping_categories)):
                        category = unit_mapping_categories[i]
                        output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, unit_mapping_dict[category][unit]))
                    output.write('''        </attvalues>\n''')
                    output.write('''    </edge>\n''')
                else:
                    output.write('''    <edge id="%d" source="%s" target="%s" weight="%f" />\n''' % (edge_id, unit, sample, unit_percents[sample][i]))


                edge_id += 1
    output.write('''</edges>\n''')
    output.write('''</graph>\n''')
    output.write('''</gexf>\n''')
    
    output.close()


def generate_gexf_network_file_for_nodes_topology(nodes_dict, output_file, attribute_types_dict = {}, project = None, skip_node_labels = False):
    output = open(output_file, 'w')

    nodes = sorted(nodes_dict.keys())
    nodes_mapping_categories = sorted([k for k in nodes_dict[nodes[0]].keys() if k not in ['children', 'parent']]) if nodes_dict else None

    output.write('''<?xml version="1.0" encoding="UTF-8"?>\n''')
    output.write('''<gexf xmlns:viz="http:///www.gexf.net/1.1draft/viz" xmlns="http://www.gexf.net/1.2draft" version="1.2">\n''')
    output.write('''<meta lastmodifieddate="2010-01-01+23:42">\n''')
    output.write('''    <creator>Oligotyping pipeline</creator>\n''')
    if project:
        output.write('''    <creator>Topology description for %s</creator>\n''' % (project))
    output.write('''</meta>\n''')
    output.write('''<graph type="static" defaultedgetype="undirected">\n\n''')

    if nodes_mapping_categories:
        output.write('''<attributes class="node" type="static">\n''')
        for i in range(0, len(nodes_mapping_categories)):
            category = nodes_mapping_categories[i]
            attr_type = attribute_types_dict[category] if attribute_types_dict.has_key(category) else "string"
            output.write('''    <attribute id="%d" title="%s" type="%s" />\n''' % (i, category, attr_type))
        output.write('''</attributes>\n\n''')

    output.write('''<nodes>\n''')
    for node in nodes:
        if skip_node_labels:
            output.write('''    <node id="%s">\n''' % (node))
        else:
            output.write('''    <node id="%s" label="%s">\n''' % (node, node))
        output.write('''        <viz:size value="%.2f"/>\n''' % (math.log10(nodes_dict[node]['size'] + 1) + 1))

        #if nodes_dict and nodes_dict.has_key('colors'):
        #    output.write('''        <viz:color r="%d" g="%d" b="%d" a="1"/>\n''' %\
        #                                     HTMLColorToRGB(nodes_dict['colors'][node], scaled = False))

        if nodes_mapping_categories:
            output.write('''        <attvalues>\n''')
            for i in range(0, len(nodes_mapping_categories)):
                category = nodes_mapping_categories[i]
                output.write('''            <attvalue id="%d" value="%s"/>\n''' % (i, nodes_dict[node][category]))
            output.write('''        </attvalues>\n''')

        output.write('''    </node>\n''')

    output.write('''</nodes>\n''')

    edge_id = 0
    output.write('''<edges>\n''')
    for node in nodes:
        if node == 'root':
            continue
        output.write('''    <edge id="%d" source="%s" target="%s" weight="1" />\n''' % (edge_id, node, nodes_dict[node]['parent']))
        edge_id += 1

    output.write('''</edges>\n''')
    output.write('''</graph>\n''')
    output.write('''</gexf>\n''')

    output.close()



def get_qual_stats_dict(quals_dict, output_file_path = None, verbose = True):
    """This function takes quals dict (which can be obtained by calling the
       utils.utils.get_quals_dict function) and returns a dictionary that
       simply contains the summary of quality scores per location in the
       alignment"""

    # FIXME: get_quals_dict and get_qual_stats_dict functions are only for
    #        454 technology at this moment.

    progress = Progress()
    progress.verbose = verbose
    progress.new('Summary of quality scores per column is being computed')
    
    qual_stats_dict = {}
    alignment_length = len(quals_dict[quals_dict.keys()[0]])
    for pos in range(0, alignment_length):
        progress.update('Position: %d of %d' % (pos + 1, alignment_length))

        qual_stats_dict[pos] = {}
        quals_for_pos = [q[pos] for q in quals_dict.values() if q[pos]]
        if not quals_for_pos:
            qual_stats_dict[pos] = None
            continue
        qual_stats_dict[pos]['mean']  = np.mean(quals_for_pos)
        qual_stats_dict[pos]['std']   = np.std(quals_for_pos)
        qual_stats_dict[pos]['max']   = np.max(quals_for_pos)
        qual_stats_dict[pos]['min']   = np.min(quals_for_pos)
        qual_stats_dict[pos]['count'] = len(quals_for_pos)
    
    if output_file_path:
        cPickle.dump(quals_dict, open(output_file_path, 'w'))

    progress.end()
    return qual_stats_dict
  
def get_quals_dict(quals_file, alignment_file, output_file_path = None, verbose = True):
    """This function takes qual scores file in FASTA format, expands each
       entry to match base calls in the corresponding aligned read in the
       FASTA file (which requires deflines to be identical), and finally
       returns a dictionary that contains qual scores as a list of integer
       values that are bound to deflines as key/value pairs"""

    quals_dict = {}
    quals_aligned_dict = {}

    progress = Progress()
    progress.verbose = verbose
    progress.new('Quality scores dictionary is being generated')
 
    alignment = u.SequenceSource(alignment_file)
    qual = u.QualSource(quals_file)

    while qual.next():
        if qual.pos % 1000 == 0:
            progress.update('Step 1 of 2 :: Quality scores read: %s' % (pretty_print(qual.pos)))
        quals_dict[qual.id] = qual.quals_int

    while alignment.next():
        if alignment.pos % 1000 == 0:
            progress.update('Step 2 of 2 :: Alignments matched: %s' % (pretty_print(alignment.pos)))
            sys.stderr.flush()

        matching_qual = quals_dict[alignment.id] 

        qual_aligned = []
        for i in range(0, len(alignment.seq)):
            if alignment.seq[i] != '-':
                qual_aligned.append(matching_qual.pop(0))
            else:
                qual_aligned.append(None)

        quals_aligned_dict[alignment.id] = qual_aligned
    progress.end()

    if output_file_path:
        cPickle.dump(quals_aligned_dict, open(output_file_path, 'w'))

    return quals_aligned_dict

def process_command_line_args_for_quality_files(args, _return = 'qual_stats_dict', verbose = True):
    """this function computes and returns the dictionary of interest (indicated
       with '_return' parameter if qual score files were provided via the command
       line interface.
       
       _return value expected to be either 'qual_stats_dict', or 'quals_dict'.

       """

    progress = Progress()
    progress.verbose = verbose

    if _return not in ['qual_stats_dict', 'quals_dict']:
        return None

    if args.qual_scores_file:
        quals_dict = get_quals_dict(args.qual_scores_file,\
                                    args.alignment,\
                                    output_file_path = args.qual_scores_file + '.cPickle')
        if _return == 'quals_dict':
            return quals_dict

        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                                              output_file_path = args.qual_scores_file + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_scores_dict:
        quals_dict = cPickle.load(open(args.qual_scores_dict))

        if _return == 'quals_dict':
            return quals_dict

        qual_stats_dict = get_qual_stats_dict(quals_dict,\
                            output_file_path = args.qual_scores_dict.split('.cPickle')[0] + '.STATS.cPickle')

        if _return == 'qual_stats_dict':
            return qual_stats_dict

    elif args.qual_stats_dict:
        qual_stats_dict = cPickle.load(open(args.qual_stats_dict))
        
        if _return == 'qual_stats_dict':
            return qual_stats_dict
    
    else:
        return None
 

def get_filtered_samples_dict(units, samples, samples_dict):
    filtered_samples_dict = {}

    for sample in samples:
        filtered_samples_dict[sample] = {}
        for unit in units:
            if samples_dict[sample].has_key(unit):
                filtered_samples_dict[sample][unit] = samples_dict[sample][unit]

    return filtered_samples_dict


def get_samples_dict_from_environment_file(environment_file_path):
    samples_dict = {}
    for oligo, sample, count in [l.strip().split('\t') for l in open(environment_file_path).readlines()]:
        if samples_dict.has_key(sample):
            if samples_dict[sample].has_key(oligo):
                samples_dict[sample][oligo] += int(count)
            else:
                samples_dict[sample][oligo] = int(count)
        else:
            samples_dict[sample] = {oligo: int(count)}
    return samples_dict


def human_readable_number(n):
    postfix = ['','K','M']
    level = max(0, min(len(postfix) - 1, int(math.floor(math.log10(abs(n))/3.0))))
    return '%.0f%s' % (n/10 ** (3 * level), postfix[level])


def pretty_print(n):
    """Pretty print function for very big integers"""
    if type(n) != int:
        return n

    ret = []
    n = str(n)
    for i in range(len(n) - 1, -1, -1):
        ret.append(n[i])
        if (len(n) - i) % 3 == 0:
            ret.append(',')
    ret.reverse()
    return ''.join(ret[1:]) if ret[0] == ',' else ''.join(ret)

def same_but_gaps(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        raise ValueError, "Alignments have different lengths"
    
    for i in range(0, len(sequence1)):
        if sequence1[i] == '-' or sequence2[i] == '-':
            continue
        if sequence1[i] != sequence2[i]:
            return False
    
    return True

def trim_uninformative_gaps_from_sequences(sequence1, sequence2):
    if len(sequence1) != len(sequence2):
        raise ValueError, "Alignments have different lengths"
    
    columns_to_discard = []
    
    for i in range(0, len(sequence1)):
        if set([sequence1[i], sequence2[i]]) == set(['-']):
            columns_to_discard.append(i)
        else:
            continue
            
    s1 = ''.join([sequence1[i] for i in range(0, len(sequence1)) if i not in columns_to_discard])
    s2 = ''.join([sequence2[i] for i in range(0, len(sequence1)) if i not in columns_to_discard])
    
    return (s1, s2)


def get_temporary_file_names_for_BLAST_search(prefix, directory):
    query  = get_temporary_file_name(prefix='%s' % prefix, suffix='.fa', directory=directory)
    target = get_temporary_file_name(prefix='%s' % prefix, suffix='.db', directory=directory)
    output = get_temporary_file_name(prefix='%s' % prefix, suffix='.b6', directory=directory)
        
    return (query, target, output)


def get_percent_identity_for_N_base_difference(average_read_length, N = 1):
    percent_identity = (average_read_length - N) * 100.0 / average_read_length
    if float('%.2f' % percent_identity) == 1.00:
        percent_identity = 0.99
        
    return percent_identity


def is_program_exist(program):
    IsExe = lambda p: os.path.isfile(p) and os.access(p, os.X_OK)

    fpath, fname = os.path.split(program)

    if fpath:
        if IsExe(program):
            return True
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if IsExe(exe_file):
                return True

    return None

def trim_uninformative_columns_from_alignment(input_file_path):
    input_fasta = u.SequenceSource(input_file_path, lazy_init = False)
    input_fasta.next()
    fasta_read_len = len(input_fasta.seq)
    invalid_columns = range(0, fasta_read_len)
    input_fasta.reset()
    
    while input_fasta.next():
        cols_not_invalid = []
        for i in invalid_columns:
            if input_fasta.seq[i] != '-':
                cols_not_invalid.append(i)
        invalid_columns = [i for i in invalid_columns if not i in cols_not_invalid]

    columns_to_keep = [x for x in range(0, fasta_read_len) if x not in invalid_columns]

    input_fasta.reset()

    temp_file = tempfile.NamedTemporaryFile(delete = False)
    temp_file_path = temp_file.name
    temp_file.close()

    temp_file = u.FastaOutput(temp_file_path)

    while input_fasta.next():
        new_seq = ''
        for i in columns_to_keep:
            new_seq += input_fasta.seq[i]
        temp_file.write_id(input_fasta.id)
        temp_file.write_seq(new_seq, split = False)
    
    temp_file.close()

    # overwrite the original file with trimmed content
    shutil.move(temp_file_path, input_file_path)

def get_temporary_file_name(prefix = '', suffix = '', directory = '/tmp'):
    temp_file_obj = tempfile.NamedTemporaryFile(prefix=prefix, suffix=suffix, dir=directory, delete=False)
    temp_file_path = temp_file_obj.name
    temp_file_obj.close()
    return temp_file_path

def get_date():
    return time.strftime("%d %b %y %H:%M:%S", time.localtime())

def get_terminal_size():
    """function was taken from http://stackoverflow.com/a/566752"""
    def ioctl_GWINSZ(fd):
        try:
            cr = struct.unpack('hh', fcntl.ioctl(fd, termios.TIOCGWINSZ,
        '1234'))
        except:
            return None
        return cr
    cr = ioctl_GWINSZ(0) or ioctl_GWINSZ(1) or ioctl_GWINSZ(2)
    if not cr:
        try:
            fd = os.open(os.ctermid(), os.O_RDONLY)
            cr = ioctl_GWINSZ(fd)
            os.close(fd)
        except:
            pass
    if not cr:
        try:
            cr = (os.environ['LINES'], os.environ['COLUMNS'])
        except:
            cr = (25, 80)
    return int(cr[1]), int(cr[0])


def estimate_expected_max_frequency_of_an_erronous_unique_sequence(number_of_reads, average_read_length, expected_error = 1/250.0):
    # maximum number of occurence of an error driven unique sequence among N reads.
    # of course this maximum assumes that all reads are coming from one template,
    # substitution probabilities are homogeneous and there are no systemmatical errors,
    # so it is a mere approximation, but for our purpose, it is going to be enough:

    return round((expected_error * (1 / 3.0)) * ((1 - expected_error) ** (average_read_length - 1)) * number_of_reads) 


def HTMLColorToRGB(colorstring, scaled = True):
    """ convert #RRGGBB to an (R, G, B) tuple """
    colorstring = colorstring.strip()
    if colorstring[0] == '#': colorstring = colorstring[1:]
    if len(colorstring) != 6:
        raise ValueError, "input #%s is not in #RRGGBB format" % colorstring
    r, g, b = colorstring[:2], colorstring[2:4], colorstring[4:]
    r, g, b = [int(n, 16) for n in (r, g, b)]

    if scaled:
        return (r / 255.0, g / 255.0, b / 255.0)
    else:
        return (r, g, b)


def run_command(cmdline):
    try:
        if subprocess.call(cmdline, shell = True) < 0:
            raise ConfigError, "command was terminated: '%s'" % (cmdline)
    except OSError, e:
        raise ConfigError, "command was failed for the following reason: '%s' ('%s')" % (e, cmdline) 


def check_command_output(cmdline):
    try:
        return subprocess.check_output(cmdline.split())
    except AttributeError:
        # python 2.6 workaround:
        p = subprocess.Popen(cmdline.split(), stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        output=''
        while(True):
            retcode = p.poll()
            output += p.stdout.readline() + '\n'
            if(retcode is not None):
                break
        return output


def get_sample_name_from_defline(defline, sample_name_separator = '_'):
    if sample_name_separator != '|':
        return sample_name_separator.join(defline.split('|')[0].split(sample_name_separator)[0:-1])
    else:
        return sample_name_separator.join(defline.split(sample_name_separator)[0:-1])


def check_input_alignment(alignment_path, sample_name_separator, progress_func = None):
    alignment = u.SequenceSource(alignment_path)
    samples = set([])
    previous_alignment_length = None

    while alignment.next():
        if progress_func and alignment.pos % 5000 == 0:
            progress_func.update('Reading input; %s, %s samples found'\
                                        % (pretty_print(alignment.pos),
                                           pretty_print(len(samples))))

        sample = get_sample_name_from_defline(alignment.id, sample_name_separator)
        if sample not in samples:
            samples.add(sample)
    
        # check the alignment lengths along the way:
        if previous_alignment_length:
            if previous_alignment_length != len(alignment.seq):
                raise ConfigError, "Not all reads have the same length."

        previous_alignment_length = len(alignment.seq)

    # if the number of samples we find in the alignment is more than half of the number of
    # reads in the alignment, we might be in trouble.
    if len(samples) * 2 > alignment.pos:
        sys.stderr.write("\n\n")
        sys.stderr.write("Number of samples in the alignment is more than half of the number of reads.\n")
        sys.stderr.write("This usually indicates that the sample name recovery from the defline is not\n")
        sys.stderr.write("working properly. If you believe this is normal, and your sample names\n")
        sys.stderr.write("expected to look like these, you can bypass this check with --skip-check-input\n")
        sys.stderr.write("parameter:\n\n")
                
        counter = 0
        for sample in samples:
            if counter == 10:
                break
            sys.stderr.write('\t- %s\n' % sample)
            counter += 1
        if len(samples) > 10:
            sys.stderr.write('\t- (%s more)\n' % pretty_print(len(samples) - 10))
        sys.stderr.write("\n\n")
        sys.stderr.write("If there is a problem with the recovery of sample names, please refer\n")
        sys.stderr.write("to the tutorial for the proper formatting of FASTA deflines.")
        sys.stderr.write("\n\n")
            
        alignment.close()
        return None
    if len(samples) == 1:
        sys.stderr.write("\n\n")
        sys.stderr.write("There is only one sample found in the alignment file during the initial check.\n")
        sys.stderr.write("If this is expected, and the following sample is the only sample in the file,\n")
        sys.stderr.write("please bypass this check by declaring --skip-check-input parameter:\n\n")
                
        for sample in samples:
            sys.stderr.write('\t- %s\n' % sample)

        sys.stderr.write("\n\n")
        sys.stderr.write("If there is a problem with the recovery of sample names, please refer\n")
        sys.stderr.write("to the tutorial for the proper formatting of FASTA deflines.")
        sys.stderr.write("\n\n")
            
        alignment.close()
        return None
    else:
        alignment.close()
        return samples


def mapping_file_simple_check(mapping_file_path, samples_expected = None):
    mapping_file = open(mapping_file_path)
    header_line = mapping_file.readline()
    
    if header_line.find('\t') < 0:
        raise ConfigError, "Mapping file doesn't seem to be a TAB delimited file"

    header_fields = header_line.strip('\n').split('\t')

    if len(header_fields) < 2:
        raise ConfigError, "No categories were found in the mapping file"
    if header_fields[0] != 'samples':
        raise ConfigError, "First column of the first row of mapping file must be 'samples'"
    if len(header_fields) != len(set(header_fields)):
        raise ConfigError, "In the mapping file, every category must be unique"

    samples_found = []
    num_entries = 0
    for line in mapping_file.readlines():
        num_entries += 1
        
        fields = line.strip('\n').split('\t')
        
        samples_found.append(fields[0])
 
        if len(fields) != len(header_fields):
            raise ConfigError, "Not every line in the mapping file has the same number of fields " +\
                                "(line %d has %d columns)" % (num_entries + 1, len(fields))
        for field in fields[1:]:
            if field == "":
                continue
            if field[0] in '0123456789':
                raise ConfigError, "Categories in the mapping file cannot start with digits: '%s'" % field

    if samples_expected:
        samples_missing = [] 
        for sample in samples_expected:
            if not sample in samples_found:
                samples_missing.append(sample)
        
        if samples_missing:
            raise ConfigError, "Mapping file seems to be missing %d sample(s) that appear in the FASTA file:\n\n- %s\n\n"\
                                      % (len(samples_missing), ', '.join(samples_missing))

    if num_entries < 3:
        raise ConfigError, "Mapping file seems to have less than three samples"

    mapping_file.close()
    return True


def get_sample_mapping_dict(mapping_file_path):
    mapping_dict = {}
    mapping_file = open(mapping_file_path)
    
    header_line = mapping_file.readline()
    categories = header_line.strip('\n').split('\t')[1:]
    for category in categories:
        mapping_dict[category] = {}
    
    for fields in [line.strip('\n').split('\t') for line in mapping_file.readlines()]:
        sample = fields[0]
        mappings = fields[1:]
        
        for i in range(0, len(categories)):
            category = categories[i]
            mapping = mappings[i]
            
            if mapping == '':
                mapping_dict[categories[i]][sample] = None
                continue
            else:
                mapping_dict[categories[i]][sample] = mapping        
            
    mapping_file.close()
    return mapping_dict

def store_filtered_matrix(old_matrix_path, new_matrix_path, samples):
    new_matrix = open(new_matrix_path, 'w')
    old_matrix = open(old_matrix_path, 'r')

    num_lines_written = 0                
    new_matrix.write(old_matrix.readline())
    for line in old_matrix.readlines():
        if line.split('\t')[0] in samples:
            num_lines_written += 1
            new_matrix.write(line)
    
    new_matrix.close()
    old_matrix.close()

    return num_lines_written


class Progress:
    def __init__(self):
        self.pid = None
        self.verbose = True
        self.terminal_width = None

        self.get_terminal_width()
        self.color_prefix = '\033[0;30m\033[46m'
        self.color_postfix = '\033[0m'
        
        self.currently_shown = None


    def get_terminal_width(self):
        try:
            self.terminal_width = get_terminal_size()[0]
        except:
            self.terminal_width = 80 


    def new(self, pid):
        if self.pid:
            raise LibError, "Progress.new() can't be called before ending the previous one (Existing: '%s', Competing: '%s')." % (self.pid, pid)

        if not self.verbose:
            return

        self.pid = '%s %s' % (get_date(), pid)
        self.get_terminal_width()
        self.currently_shown = None


    def write(self, c):
        surpass = self.terminal_width - len(c)
        
        if surpass < 0:
            c = c[0:-(-surpass + 4)] + ' (...)'
        else:
            self.currently_shown = c
            c = c + ' ' * surpass

        if self.verbose:
            sys.stderr.write(self.color_prefix + c + self.color_postfix)
            sys.stderr.flush()
            

    def reset(self):
        self.clear()

    def clear(self):
        if not self.verbose:
            return
        null = '\r' + ' ' * (self.terminal_width) 
        sys.stderr.write(null)
        sys.stderr.write('\r')
        sys.stderr.flush()
        self.currently_shown = None


    def append(self, msg):
        if not self.verbose:
            return
        self.write('%s%s' % (self.currently_shown, msg))


    def update(self, msg):
        if not self.verbose:
            return
        self.clear()
        self.write('\r[%s] %s' % (self.pid, msg))

    
    def end(self):
        self.pid = None
        if not self.verbose:
            return
        self.clear()


def get_pretty_name(key):
    if pretty_names.has_key(key):
        return pretty_names[key]
    else:
        return key

    
def get_cmd_line(argv):
    c_argv = []
    for i in argv:
        if ' ' in i:
            c_argv.append('"%s"' % i)
        else:
            c_argv.append(i)
    return ' '.join(c_argv)


class Run:
    """a class that keeps info about an oligotyping run, and deal with the console output"""
    def __init__(self, info_file_path = None, verbose = True):
        if info_file_path:
            self.init_info_file_obj(info_file_path)
        else:
            self.info_file_obj = None

        self.info_dict = {}
        self.verbose = verbose


    def init_info_file_obj(self, info_file_path):
            self.info_file_obj = open(info_file_path, 'w')


    def info(self, key, value, quiet = False):
        self.info_dict[key] = value
        
        if quiet:
            return True
        
        if type(value) == int:
            value = pretty_print(value)

        label = get_pretty_name(key)

        info_line = "%s %s: %s\n" % (label, '.' * (65 - len(label)), str(value))
        if self.info_file_obj:
            self.info_file_obj.write(info_line)

        if self.verbose:
            sys.stderr.write(info_line)


    def store_info_dict(self, destination):
        cPickle.dump(self.info_dict, open(destination, 'w'))


    def quit(self):
        if self.info_file_obj:
            self.info_file_obj.close()

def get_read_objects_from_file(input_file_path):
    input_fasta = u.SequenceSource(input_file_path, unique = True)
    read_objects = []
    
    while input_fasta.next():
        read_objects.append(UniqueFASTAEntry(input_fasta.seq, input_fasta.ids))

    input_fasta.close()
    return read_objects


def split_fasta_file(input_file_path, dest_dir, prefix = 'part', num_reads_per_file = 5000):
    input_fasta = u.SequenceSource(input_file_path)
    
    parts = []
    next_part = 1
    part_obj = None

    while input_fasta.next():
        if (input_fasta.pos - 1) % num_reads_per_file == 0:
            if part_obj:
                part_obj.close()

            rand_bit = ''.join([random.choice(string.ascii_letters + string.digits) for n in xrange(8)])
            file_path = os.path.join(dest_dir, '%s-%d-%s.fa' % (prefix, next_part, rand_bit))
            parts.append(file_path)
            next_part += 1
            part_obj = u.FastaOutput(file_path)

        part_obj.store(input_fasta, split = False)
  
    if part_obj:
        part_obj.close()

    return parts


class Multiprocessing:
    def __init__(self, target_function, num_thread = None):
        self.cpu_count = multiprocessing.cpu_count()
        self.num_thread = num_thread or (self.cpu_count - (int(round(self.cpu_count / 10.0)) or 1))
        self.target_function = target_function
        self.processes = []
        self.manager = multiprocessing.Manager()


    def get_data_chunks(self, data_array, spiral = False):
        data_chunk_size = (len(data_array) / self.num_thread) or 1
        data_chunks = []
        
        if len(data_array) <= self.num_thread:
            return [[chunk] for chunk in data_array]

        if spiral:
            for i in range(0, self.num_thread):
                data_chunks.append([data_array[j] for j in range(i, len(data_array), self.num_thread)])
            
            return data_chunks
        else:
            for i in range(0, self.num_thread):
                if i == self.num_thread - 1:
                    data_chunks.append(data_array[i * data_chunk_size:])
                else:
                    data_chunks.append(data_array[i * data_chunk_size:i * data_chunk_size + data_chunk_size])

        return data_chunks

                
    def run(self, args, name = None):
        t = multiprocessing.Process(name = name,
                                    target = self.target_function,
                                    args = args)
        self.processes.append(t)
        t.start()


    def get_empty_shared_array(self):
        return self.manager.list()


    def get_empty_shared_dict(self):
        return self.manager.dict()

    
    def get_shared_integer(self):
        return self.manager.Value('i', 0)


    def run_processes(self, processes_to_run, progress_obj = None):
        tot_num_processes = len(processes_to_run)
        sent_to_run = 0
        while 1:
            NumRunningProceses = lambda: len([p for p in self.processes if p.is_alive()])
            
            if NumRunningProceses() < self.num_thread and processes_to_run:
                for i in range(0, self.num_thread - NumRunningProceses()):
                    if len(processes_to_run):
                        sent_to_run += 1
                        self.run(processes_to_run.pop())

            if not NumRunningProceses() and not processes_to_run:
                # let the blastn program finish writing all output files.
                # FIXME: this is ridiculous. find a better solution.
                time.sleep(1)
                break

            if progress_obj:
                progress_obj.update('%d of %d done in %d threads (currently running processes: %d)'\
                                                         % (sent_to_run - NumRunningProceses(),
                                                            tot_num_processes,
                                                            self.num_thread,
                                                            NumRunningProceses()))
            time.sleep(1)


class UniqueFASTAEntry:
    def __init__(self, seq, ids):
        self.seq = seq
        self.ids = ids
        self.md5id = hashlib.md5(self.seq).hexdigest()
        self.frequency = len(ids)
