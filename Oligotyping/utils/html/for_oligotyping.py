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
import copy
import shutil
import cPickle

from Oligotyping.lib import fastalib as u
from Oligotyping.utils.constants import pretty_names
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import get_samples_dict_from_environment_file
from Oligotyping.utils.random_colors import get_list_of_colors
from error import HTMLError


try:
    from django.template.loader import render_to_string
    from django.conf import settings
    from django.template.defaultfilters import register
except ImportError:
    raise HTMLError, 'You need to have Django module (http://djangoproject.com) installed on your system to generate HTML output.'

@register.filter(name='diffs')
def diffs(l, index):
    """this filter is being used in oligo.tmpl, takes first 20 unique sequence
       as an argument (l) and compares l[index] to l[0] to return differences"""
    _diffs = []
    for i in range(0, len(l[0])):
        if l[0][i] != l[index][i]:
            _diffs.append(i)

    return ', '.join([str(i) for i in _diffs])

@register.filter(name='lookup')
def lookup(d, index):
    if index in d:
        return d[index]
    return ''

@register.filter(name='get_list_item')
def get_list_item(l, index):
    if index < len(l):
        return l[index]
    return ''

@register.filter(name='has_perfect_hit')
def has_perfect_hit(b6_entry_list):
    return b6_entry_list[0].identity == 100.0 if len(b6_entry_list) else False

@register.filter(name='get_blast_hits')
def get_blast_hits(d, max_num = 8):
    '''gets a dictionary of BLAST results, returns
       the target_labels where 100% identity is
       achieved'''

    num_show = len(d) if len(d) < max_num else max_num

    if num_show == 0:
        return ''
        
    ret_line = '<p><b>BLAST search results at a glance</b> (%d of %d total hits are shown):' %\
                                            (num_show, len(d))
    for i in range(0, num_show):
        entry = d[i]
        if entry.identity == 100.0:
            ret_line += '<p>* %s (<b><i>identity: %.2f%%, query coverage: %.2f%%</i></b>)' \
                                    % (entry.hit_def.replace("'", '"'),
                                       entry.identity,
                                       entry.coverage)
        else:
            ret_line += '<p>* %s (<i>identity: %.2f%%, query coverage: %.2f%%</i>)' \
                                    % (entry.hit_def.replace("'", '"'),
                                       entry.identity,
                                       entry.coverage)
    return ret_line

@register.filter(name='as_percentage_of') 
def as_percentage_of(part, whole):
    try:
        return "%.2f%%" % (float(part) / whole * 100.0)
    except (ValueError, ZeroDivisionError):
        return ""

@register.filter(name='percentify') 
def percentify(l):
    total = sum(l)
    if total:
        return [p * 100.0 / total for p in l]
    else:
        return [0] * len(l)

@register.filter(name='presicion') 
def presicion(value, arg):
    if value == 0:
        return '0.' + '0' * arg
    else:
        t = '%' + '.%d' % arg + 'f'
        return t % value

@register.filter(name='sorted_by_value') 
def sorted_by_value(d):
    return sorted(d, key=d.get, reverse=True)

@register.filter(name='get_colors') 
def get_colors(number_of_colors):
    return get_list_of_colors(number_of_colors, colormap="Dark2")

@register.filter(name='values') 
def values(d):
    return d.values()

@register.filter(name='mod') 
def mod(value, arg):
    return value % arg 

@register.filter(name='pretify') 
def pretify(arg):
    return pretty_print(arg) 

@register.filter(name='multiply') 
def multiply(value, arg):
    return int(value) * int(arg) 

@register.filter(name='var') 
def var(arg):
    return 'x_' + arg.replace(' ', '_').replace('-', '_').replace('+', '_').replace('.', '_') 

@register.filter(name='cleangaps') 
def celangaps(arg):
    return arg.replace('-', '')

@register.filter(name='sumvals') 
def sumvals(arg, clean = None):
    if clean:
        return sum(arg.values())
    return pretty_print(sum(arg.values()))

@register.filter(name='mklist') 
def mklist(arg):
    return range(0, int(arg))

absolute = os.path.join(os.path.dirname(os.path.realpath(__file__)))
settings.configure(DEBUG=True, TEMPLATE_DEBUG=True, DEFAULT_CHARSET='utf-8', TEMPLATE_DIRS = (os.path.join(absolute, 'templates'),))

from django.template.loader import get_template
t = get_template('index_for_oligo.tmpl')

def generate_html_output(run_info_dict, html_output_directory = None, entropy_figure = None):
    if not html_output_directory:    
        html_output_directory = os.path.join(run_info_dict['output_directory'], 'HTML-OUTPUT')
        
    if not os.path.exists(html_output_directory):
        os.makedirs(html_output_directory)
    
    html_dict = copy.deepcopy(run_info_dict)

    shutil.copy2(os.path.join(absolute, 'static/style.css'), os.path.join(html_output_directory, 'style.css'))
    shutil.copy2(os.path.join(absolute, 'static/header_1.png'), os.path.join(html_output_directory, 'header.png'))
    shutil.copy2(os.path.join(absolute, 'static/missing_image.png'), os.path.join(html_output_directory, 'missing.png'))
    shutil.copy2(os.path.join(absolute, 'static/colorbar.png'), os.path.join(html_output_directory, 'colorbar.png'))
    shutil.copy2(os.path.join(absolute, 'scripts/jquery-1.7.1.js'), os.path.join(html_output_directory, 'jquery-1.7.1.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/popup.js'), os.path.join(html_output_directory, 'popup.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/g.pie.js'), os.path.join(html_output_directory, 'g.pie.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/g.raphael.js'), os.path.join(html_output_directory, 'g.raphael.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/raphael.js'), os.path.join(html_output_directory, 'raphael.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/morris.js'), os.path.join(html_output_directory, 'morris.js'))

    def copy_as(source, dest_name, essential = True):
        dest = os.path.join(html_output_directory, dest_name)

        if essential:
            shutil.copy2(source, dest)
        else:
            # it is ok if you fail to copy files that are not
            # essential.. 
            try:
                shutil.copy2(source, dest)
            except:
                sys.stderr.write('\n\n[HTML] Warning: Source file not found\n\tSource: "%s"\n\tDest: "%s\n\n"' % (source, dest))

        return os.path.basename(dest)

    # embarrassingly ad-hoc:
    if entropy_figure:
        if entropy_figure.endswith('.pdf') or entropy_figure.endswith('.png'):
            entropy_figure = entropy_figure[:-4]
            
    CP = lambda e, o:  copy_as(os.path.join(e + ('.%s' % ext)), o, essential = True if ext == 'png' else False)
    for ext in ['png', 'pdf']:
        output_file = 'entropy.%s' % ext
        if entropy_figure:
            html_dict['entropy_figure_%s' % ext] = CP(entropy_figure, output_file)
        else:
            try:
                html_dict['entropy_figure_%s' % ext] = CP(run_info_dict['entropy'], output_file)
            except:
                html_dict['entropy_figure_%s' % ext] = CP(run_info_dict['entropy'][:-4], output_file)

 
    if run_info_dict['gexf_network_file_path']:
        html_dict['gexf_network_file_path'] = copy_as(run_info_dict['gexf_network_file_path'], 'network.gexf')

    if run_info_dict['sample_mapping']:
        html_dict['sample_mapping'] = copy_as(run_info_dict['sample_mapping'], 'sample_mapping.txt')
    else:
        html_dict['sample_mapping'] = None

    html_dict['matrix_count_file_path'] = copy_as(run_info_dict['matrix_count_file_path'], 'matrix_counts.txt')
    html_dict['matrix_percent_file_path'] = copy_as(run_info_dict['matrix_percent_file_path'], 'matrix_percents.txt')
    html_dict['read_distribution_table_path'] = copy_as(run_info_dict['read_distribution_table_path'], 'read_distribution.txt')
    html_dict['environment_file_path'] = copy_as(run_info_dict['environment_file_path'], 'environment.txt')
    html_dict['oligos_fasta_file_path'] = copy_as(run_info_dict['oligos_fasta_file_path'], 'oligos.fa.txt')
    html_dict['oligos_nexus_file_path'] = copy_as(run_info_dict['oligos_nexus_file_path'], 'oligos.nex.txt')


    def get_figures_dict(html_dict_prefix):
        html_dict_key = '%s_file_path' % html_dict_prefix
        if html_dict.has_key(html_dict_key):
            figures_dict = cPickle.load(open(html_dict[html_dict_key]))
            for _map in figures_dict:
                for _func in figures_dict[_map]:
                    for _op in figures_dict[_map][_func]:
                        if os.path.exists(figures_dict[_map][_func][_op] + '.pdf') and os.path.exists(figures_dict[_map][_func][_op] + '.png'):
                            prefix = copy_as(figures_dict[_map][_func][_op] + '.pdf', '%s.pdf' % '-'.join([_map, _func, _op]))
                            prefix = copy_as(figures_dict[_map][_func][_op] + '.png', '%s.png' % '-'.join([_map, _func, _op]))
                            figures_dict[_map][_func][_op] = '.'.join(prefix.split('.')[:-1])
                        else:
                            figures_dict[_map][_func][_op] = None
            return figures_dict
        else:
            return None
        
    
    html_dict['figures_dict'] = get_figures_dict('figures_dict')
    html_dict['exclusive_figures_dict'] = get_figures_dict('exclusive_figures_dict')


    if html_dict['generate_sets']:
        html_dict['across_samples_MN_file_path'] = copy_as(run_info_dict['across_samples_MN_file_path'], 'across_samples_max_normalized.txt')
        html_dict['across_samples_SN_file_path'] = copy_as(run_info_dict['across_samples_SN_file_path'], 'across_samples_sum_normalized.txt')
        html_dict['oligo_sets_stackbar_figure'] = copy_as(run_info_dict['stack_bar_with_agglomerated_oligos_file_path'], 'stackbar_with_oligo_sets.png')
        html_dict['oligos_across_samples_figure'] = copy_as(run_info_dict['oligos_across_samples_file_path'], 'oligos_across_samples.png')
        html_dict['oligotype_sets_figure'] = copy_as(run_info_dict['oligotype_sets_across_samples_figure_path'], 'oligotype_sets.png')
        html_dict['matrix_count_oligo_sets_file_path'] = copy_as(run_info_dict['matrix_count_oligo_sets_file_path'], 'matrix_counts_oligo_sets.txt')
        html_dict['matrix_percent_oligo_sets_file_path'] = copy_as(run_info_dict['matrix_percent_oligo_sets_file_path'], 'matrix_percents_oligo_sets.txt')
        html_dict['oligotype_sets_file'] = copy_as(run_info_dict['oligotype_sets_file_path'], 'oligotype_sets.txt')
        html_dict['oligotype_sets'] = [l.strip().split('\t')[1].split(',') for l in open(run_info_dict['oligotype_sets_file_path'])]
 
    if html_dict.has_key('representative_seqs_fasta_file_path'):
        html_dict['representative_seqs_fasta_file_path'] = copy_as(run_info_dict['representative_seqs_fasta_file_path'], 'oligo-representatives.fa.txt')
    else:
        html_dict['representative_seqs_fasta_file_path'] = None
    if run_info_dict.has_key('blast_ref_db') and os.path.exists(run_info_dict['blast_ref_db']):
        html_dict['blast_ref_db_path'] = copy_as(run_info_dict['blast_ref_db'], 'reference_db.fa')
    html_dict['entropy_components'] = [int(x) for x in html_dict['bases_of_interest_locs'].split(',')]
    html_dict['samples_dict'] = get_samples_dict_from_environment_file(run_info_dict['environment_file_path'])
    html_dict['samples'] = sorted(html_dict['samples_dict'].keys())
    html_dict['blast_results_found'] = False

    # get alignment length
    html_dict['alignment_length'] = get_alignment_length(run_info_dict['alignment'])
    # include pretty names
    html_dict['pretty_names'] = pretty_names
    # get purity score colors dict
    html_dict['score_color_dict'] = {}
    gradient = get_list_of_colors(26, colormap = 'RdYlGn')
    for oligo in run_info_dict['final_purity_score_dict']:
        html_dict['score_color_dict'][oligo] = gradient[int(run_info_dict['final_purity_score_dict'][oligo] * 25)]
    # get total purity score color dict
    html_dict['total_score_color'] = gradient[int(float(run_info_dict['total_purity_score_dict']) * 25)]
    # get colors dict
    html_dict['color_dict'] = get_colors_dict(run_info_dict['colors_file_path'])
    # get abundant oligos list
    html_dict['oligos'] = get_oligos_list(run_info_dict['oligos_fasta_file_path'])
    # get oligo frequencies
    html_dict['frequency'] = {}
    for oligo in html_dict['oligos']:
        html_dict['frequency'][oligo] = pretty_print(sum([d[oligo] for d in html_dict['samples_dict'].values() if d.has_key(oligo)]))
    # get purity score
    html_dict['purity_score'] = run_info_dict['final_purity_score_dict']
    # get total purity score
    html_dict['total_purity_score'] = run_info_dict['total_purity_score_dict']
    # get unique sequence dict (which will contain the most frequent unique sequence for given oligotype)
    if html_dict.has_key('output_directory_for_reps'):
        html_dict['rep_oligo_seqs_clean_dict'], html_dict['rep_oligo_seqs_fancy_dict'] = get_unique_sequences_dict(html_dict)
        html_dict['oligo_reps_dict'] = get_oligo_reps_dict(html_dict, html_output_directory)
        html_dict['component_reference'] = ''.join(['<a onmouseover="popup(\'\#%d\', 50)" href="">|</a>' % i for i in range(0, html_dict['alignment_length'])])

    # get javascript code for sample pie-charts
    html_dict['pie_charts_js'] = render_to_string('pie_charts_js.tmpl', html_dict)

    # FIXME: code below is very inefficient and causes a huge
    # memory issue. fix it by not using deepcopy.
    # generate individual oligotype pages
    if html_dict.has_key('output_directory_for_reps'):
        for i in range(0, len(html_dict['oligos'])):
            oligo = html_dict['oligos'][i]
            tmp_dict = copy.deepcopy(html_dict)
            tmp_dict['oligo'] = oligo
            tmp_dict['distribution'] = get_oligo_distribution_dict(oligo, html_dict)
            oligo_page = os.path.join(html_output_directory, 'oligo_%s.html' % oligo)
            
            tmp_dict['index'] = i + 1
            tmp_dict['total'] = len(html_dict['oligos'])
            tmp_dict['prev'] = None
            tmp_dict['next'] = None
            if i > 0:
                tmp_dict['prev'] = 'oligo_%s.html' % html_dict['oligos'][i - 1]
            if i < (len(html_dict['oligos']) - 1):
                tmp_dict['next'] = 'oligo_%s.html' % html_dict['oligos'][i + 1]
            
            rendered = render_to_string('single_oligo.tmpl', tmp_dict)
    
            open(oligo_page, 'w').write(rendered.encode("utf-8"))


    # generate index
    index_page = os.path.join(html_output_directory, 'index.html')
    rendered = render_to_string('index_for_oligo.tmpl', html_dict)

    open(index_page, 'w').write(rendered.encode("utf-8"))

    return index_page

def get_colors_dict(colors_file_path):
    colors_dict = {}
    for oligo, color in [line.strip().split('\t') for line in open(colors_file_path).readlines()]:
        colors_dict[oligo] = color
    return colors_dict

def get_oligos_list(oligos_file_path):
    oligos_list = []
    fasta = u.SequenceSource(oligos_file_path)
    while fasta.next():
        oligos_list.append(fasta.seq)
    return oligos_list

def get_oligo_distribution_dict(oligo, html_dict):
    rep_dir = html_dict['output_directory_for_reps']
    oligo_distribution_dict = cPickle.load(open(os.path.join(rep_dir, '%.5d_'\
        % html_dict['oligos'].index(oligo) + oligo + '_unique_distribution.cPickle')))
    
    ret_dict = {}

    for sample in oligo_distribution_dict:
        ret_dict[sample] = [0] * 20
        for i in range(0, 20):
            if oligo_distribution_dict[sample].has_key(i + 1):
                ret_dict[sample][i] = oligo_distribution_dict[sample][i + 1]

    return ret_dict


def get_oligo_reps_dict(html_dict, html_output_directory):
    oligos, rep_dir = html_dict['oligos'], html_dict['output_directory_for_reps']

    oligo_reps_dict = {}
    oligo_reps_dict['imgs'] = {}
    oligo_reps_dict['fancy_seqs'] = {}
    oligo_reps_dict['clear_seqs'] = {}
    oligo_reps_dict['frequency'] = {}
    oligo_reps_dict['component_references'] = {}
    oligo_reps_dict['blast_results'] = {}

    for i in range(0, len(oligos)):
        oligo = oligos[i]

        alignment_base_path = os.path.join(rep_dir, '%.5d_' % i + oligo)

        diversity_image_path =  alignment_base_path + '_unique.png'
        diversity_image_dest = os.path.join(html_output_directory, os.path.basename(diversity_image_path))
        shutil.copy2(diversity_image_path, diversity_image_dest)
        oligo_reps_dict['imgs'][oligo] = os.path.basename(diversity_image_dest)

        unique_sequences_path = alignment_base_path + '_unique'
        uniques = u.SequenceSource(unique_sequences_path)
        oligo_reps_dict['fancy_seqs'][oligo] = []
        oligo_reps_dict['clear_seqs'][oligo] = []
        oligo_reps_dict['frequency'][oligo] = []
        while uniques.next() and uniques.pos <= 20:
            oligo_reps_dict['clear_seqs'][oligo].append(uniques.seq)
            oligo_reps_dict['fancy_seqs'][oligo].append(get_decorated_sequence(uniques.seq, html_dict['entropy_components']))
            oligo_reps_dict['frequency'][oligo].append(pretty_print(uniques.id.split('|')[1].split(':')[1]))

        entropy_file_path = alignment_base_path + '_unique_entropy'
        entropy_values_per_column = [0] * html_dict['alignment_length']
        for column, entropy in [x.strip().split('\t') for x in open(entropy_file_path)]:
            entropy_values_per_column[int(column)] = float(entropy)

        color_per_column = cPickle.load(open(alignment_base_path + '_unique_color_per_column.cPickle'))
        oligo_reps_dict['component_references'][oligo] = ''.join(['<span style="background-color: %s;"><a onmouseover="popup(\'\column: %d<br />entropy: %.4f\', 100)" href="">|</a></span>' % (color_per_column[i], i, entropy_values_per_column[i]) for i in range(0, html_dict['alignment_length'])])

        blast_results_dict = alignment_base_path + '_unique_BLAST.cPickle'
        if os.path.exists(blast_results_dict):
            html_dict['blast_results_found'] = True
            oligo_reps_dict['blast_results'][oligo] = cPickle.load(open(blast_results_dict))
        else:
            oligo_reps_dict['blast_results'][oligo] = None

    return oligo_reps_dict

def get_alignment_length(alignment_path):
    alignment = u.SequenceSource(alignment_path)
    alignment.next()
    return len(alignment.seq)

def get_unique_sequences_dict(html_dict):
    oligos, rep_dir = html_dict['oligos'], html_dict['output_directory_for_reps']

    rep_oligo_seqs_clean_dict = {}
    rep_oligo_seqs_fancy_dict = {}
    
    for i in range(0, len(oligos)):
        unique_file_path = os.path.join(rep_dir, '%.5d_' % i + oligos[i] + '_unique')
        f = u.SequenceSource(unique_file_path)
        f.next()
        rep_oligo_seqs_clean_dict[oligos[i]] = f.seq
        rep_oligo_seqs_fancy_dict[oligos[i]] = get_decorated_sequence(f.seq, html_dict['entropy_components'])
        f.close()
    return (rep_oligo_seqs_clean_dict, rep_oligo_seqs_fancy_dict)

def get_decorated_sequence(seq, components):
    """returns sequence with html decorations"""
    return ''.join(map(lambda j: '<span class="c">%s</span>' % seq[j] if j in components else seq[j], [j for j in range(len(seq))]))

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate Static HTML Output from Oligotyping Run')
    parser.add_argument('run_info_dict_path', metavar = 'DICT', help = 'Serialized run info dictionary (RUNINFO.cPickle)')
    parser.add_argument('-o', '--output-directory', default = None, metavar = 'OUTPUT_DIR',\
                        help = 'Output directory for HTML output to be stored')
    parser.add_argument('--entropy-figure', default = None, metavar = 'ENTROPY_FIGURE',\
                        help = 'Path for entropy figure *without* the file extension (e.g. only "/path/to/entropy" \
                               for "/path/to/entropy.png")')

    args = parser.parse_args()
   
    run_info_dict = cPickle.load(open(args.run_info_dict_path))

    index_page = generate_html_output(run_info_dict, args.output_directory, args.entropy_figure) 

    print '\n\tHTML output is ready: "%s"\n' % index_page
