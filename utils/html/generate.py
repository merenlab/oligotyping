#!/usr/bin/python
# -*- coding: utf-8 -*-

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
import shutil
import copy

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../..'))
import lib.fastalib as u
from utils.constants import pretty_names

from error import HTMLError

try:
    from django.template.loader import render_to_string
    from django.conf import settings
except ImportError:
    raise HTMLError, 'You need to have Django module (http://djangoproject.com) installed on your system to generate HTML output.'


# a snippet from StackOverflow for dict lookups in
# templates
from django.template.defaultfilters import register
@register.filter(name='lookup')
def lookup(dict, index):
    if index in dict:
        return dict[index]
    return ''

@register.filter(name='multiply') 
def multiply(value, arg):
    return int(value) * int(arg) 

absolute = os.path.join(os.path.dirname(os.path.realpath(__file__)))
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, DEFAULT_CHARSET='utf-8', TEMPLATE_DIRS = (os.path.join(absolute, 'templates'),))

from django.template.loader import get_template
t = get_template('index.tmpl')

def generate_html_output(run_info_dict, html_output_directory = None):
    if not html_output_directory:    
        html_output_directory = os.path.join(run_info_dict['output_directory'], 'HTML-OUTPUT')
        
    if not os.path.exists(html_output_directory):
        os.makedirs(html_output_directory)
    
    html_dict = copy.deepcopy(run_info_dict)

    shutil.copy2(os.path.join(absolute, 'static/style.css'), os.path.join(html_output_directory, 'style.css'))
    shutil.copy2(os.path.join(absolute, 'static/header.png'), os.path.join(html_output_directory, 'header.png'))

    def copy_as(source, dest_name):
        dest = os.path.join(html_output_directory, dest_name)
        shutil.copy2(source, dest)
        return dest

    # embarrassingly ad-hoc:
    html_dict['entropy_figure'] = copy_as(os.path.join(run_info_dict['entropy'].split('.')[0] + '.png'), 'entropy.png')
    html_dict['stackbar_figure'] = copy_as(run_info_dict['stack_bar_file_path'], 'stackbar.png')
    html_dict['matrix_count_file_path'] = copy_as(run_info_dict['matrix_count_file_path'], 'matrix_counts.txt')
    html_dict['matrix_percent_file_path'] = copy_as(run_info_dict['matrix_percent_file_path'], 'matrix_percents.txt')
    html_dict['environment_file_path'] = copy_as(run_info_dict['environment_file_path'], 'environment.txt')
    html_dict['oligos_fasta_file_path'] = copy_as(run_info_dict['oligos_fasta_file_path'], 'oligos.fa.txt')
    html_dict['oligos_nexus_file_path'] = copy_as(run_info_dict['oligos_nexus_file_path'], 'oligos.nex.txt')

    # include pretty names
    html_dict['pretty_names'] = pretty_names
    # get colors dict
    html_dict['color_dict'] = get_colors_dict(run_info_dict['random_color_file_path'])
    # get abundant oligos list
    html_dict['oligos'] = get_oligos_list(run_info_dict['oligos_fasta_file_path'])
    # get unique sequence dict (which will contain the most frequent unique sequence for given oligotype)
    if html_dict.has_key('output_directory_for_reps'):
        html_dict['sequence_dict'] = get_unique_sequences_dict(html_dict)
        html_dict['alignment_length'] = len(html_dict['sequence_dict'].values()[0])
    else:
        html_dict['sequence_dict'] = dict(zip(html_dict['oligos'], ['(representative sequences were not computed during the analysis)'] * len(html_dict['oligos'])))


    index_page = os.path.join(html_output_directory, 'index.html')
    rendered = render_to_string('index.tmpl', html_dict)

    open(index_page, 'w').write(rendered.encode("utf-8"))

    return index_page

def get_colors_dict(random_color_file_path):
    colors_dict = {}
    for oligo, color in [line.strip().split('\t') for line in open(random_color_file_path).readlines()]:
        colors_dict[oligo] = color
    return colors_dict

def get_oligos_list(oligos_file_path):
    oligos_list = []
    fasta = u.SequenceSource(oligos_file_path)
    while fasta.next():
        oligos_list.append(fasta.seq)
    return oligos_list

def get_unique_sequences_dict(html_dict):
    oligos, rep_dir = html_dict['oligos'], html_dict['output_directory_for_reps']
    entropy_components = [int(x) for x in html_dict['bases_of_interest_locs'].split(',')]

    sequence_dict = {}
    
    for i in range(0, len(oligos)):
        unique_file_path = os.path.join(rep_dir, '%.5d_' % i + oligos[i] + '_unique')
        f = u.SequenceSource(unique_file_path)
        f.next()
        sequence_dict[oligos[i]] = ''.join(map(lambda j: '<span class="c">%s</span>' % f.seq[j] if j in entropy_components else f.seq[j], [j for j in range(len(f.seq))]))
        f.close()
    return sequence_dict


if __name__ == '__main__':
    import cPickle
    import argparse

    parser = argparse.ArgumentParser(description='Generate Static HTML Output from Oligotyping Run')
    parser.add_argument('run_info_dict_path', metavar = 'DICT', help = 'Serialized run info dictionary (RUNINFO.cPickle)')
    parser.add_argument('-o', '--output-directory', default = None, metavar = 'OUTPUT_DIR',\
                        help = 'Output directory for HTML output to be stored')

    args = parser.parse_args()
   
    run_info_dict = cPickle.load(open(args.run_info_dict_path))

    index_page = generate_html_output(run_info_dict, args.output_directory) 

    print '\n\n\tHTML output is ready: "%s"\n' % index_page
