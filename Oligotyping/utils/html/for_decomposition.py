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
import copy
import shutil
import cPickle

from Oligotyping.utils.constants import pretty_names
from Oligotyping.utils.utils import pretty_print
from Oligotyping.utils.utils import get_datasets_dict_from_environment_file
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
    for i in d.keys()[0:num_show]:
        if d[i]['identity'] == 100.0:
            ret_line += '<p>* %s (<b><i>identity: %.2f%%, query coverage: %.2f%%</i></b>)' \
                                    % (d[i]['hit_def'].replace("'", '"'),
                                       d[i]['identity'],
                                       d[i]['coverage'])
        else:
            ret_line += '<p>* %s (<i>identity: %.2f%%, query coverage: %.2f%%</i>)' \
                                    % (d[i]['hit_def'].replace("'", '"'),
                                       d[i]['identity'],
                                       d[i]['coverage'])
    return ret_line

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
t = get_template('index_for_decomposition.tmpl')

def generate_html_output(run_info_dict, html_output_directory = None):
    if not html_output_directory:
        html_output_directory = os.path.join(run_info_dict['output_directory'], 'HTML-OUTPUT')
        
    if not os.path.exists(html_output_directory):
        os.makedirs(html_output_directory)
    
    html_dict = copy.deepcopy(run_info_dict)

    shutil.copy2(os.path.join(absolute, 'static/style.css'), os.path.join(html_output_directory, 'style.css'))
    shutil.copy2(os.path.join(absolute, 'static/header_2.png'), os.path.join(html_output_directory, 'header.png'))
    shutil.copy2(os.path.join(absolute, 'scripts/jquery-1.7.1.js'), os.path.join(html_output_directory, 'jquery-1.7.1.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/popup.js'), os.path.join(html_output_directory, 'popup.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/g.pie.js'), os.path.join(html_output_directory, 'g.pie.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/g.raphael.js'), os.path.join(html_output_directory, 'g.raphael.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/raphael.js'), os.path.join(html_output_directory, 'raphael.js'))
    shutil.copy2(os.path.join(absolute, 'scripts/morris.js'), os.path.join(html_output_directory, 'morris.js'))

    def copy_as(source, dest_name):
        dest = os.path.join(html_output_directory, dest_name)
        shutil.copy2(source, dest)
        return os.path.basename(dest)

    html_dict['matrix_count_file_path'] = copy_as(run_info_dict['matrix_count_file_path'], 'matrix_counts.txt')
    html_dict['matrix_percent_file_path'] = copy_as(run_info_dict['matrix_percent_file_path'], 'matrix_percents.txt')
    html_dict['environment_file_path'] = copy_as(run_info_dict['environment_file_path'], 'environment.txt')

    if html_dict.has_key('node_representatives_file_path'):
        html_dict['node_representatives_file_path'] = copy_as(run_info_dict['node_representatives_file_path'], 'node-representatives.fa.txt')
    else:
        html_dict['node_representatives_file_path'] = None

    if run_info_dict.has_key('blast_ref_db') and os.path.exists(run_info_dict['blast_ref_db']):
        html_dict['blast_ref_db_path'] = copy_as(run_info_dict['blast_ref_db'], 'reference_db.fa')

    html_dict['datasets_dict'] = get_datasets_dict_from_environment_file(run_info_dict['environment_file_path'])
    html_dict['datasets'] = sorted(html_dict['datasets_dict'].keys())
    html_dict['blast_results_found'] = False

    # include pretty names
    html_dict['pretty_names'] = pretty_names

    # get javascript code for dataset pie-charts
    html_dict['pie_charts_js'] = render_to_string('pie_charts_js.tmpl', html_dict)


    # generate index
    index_page = os.path.join(html_output_directory, 'index.html')
    rendered = render_to_string('index_for_decomposition.tmpl', html_dict)

    open(index_page, 'w').write(rendered.encode("utf-8"))

    return index_page

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Generate Static HTML Output from Oligotyping Run')
    parser.add_argument('run_info_dict_path', metavar = 'DICT', help = 'Serialized run info dictionary (RUNINFO.cPickle)')
    parser.add_argument('-o', '--output-directory', default = None, metavar = 'OUTPUT_DIR',\
                        help = 'Output directory for HTML output to be stored')

    args = parser.parse_args()
   
    run_info_dict = cPickle.load(open(args.run_info_dict_path))

    index_page = generate_html_output(run_info_dict, args.output_directory) 

    print '\n\tHTML output is ready: "%s"\n' % index_page
