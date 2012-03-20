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

from error import HTMLError

try:
    from django.template.loader import render_to_string
    from django.conf import settings
except ImportError:
    raise HTMLError, 'You need to have Django module (http://djangoproject.com) installed on your system to generate HTML output.'

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

    # embarrassingly ad-hoc:
    entropy_figure = os.path.join(html_output_directory, 'entropy.png')
    shutil.copy2(os.path.join(run_info_dict['entropy'].split('.')[0] + '.png'), entropy_figure)
    html_dict['entropy_figure'] = entropy_figure

    stackbar_figure = os.path.join(html_output_directory, 'stackbar.png')
    shutil.copy2(os.path.join(run_info_dict['stack_bar_file_path']), stackbar_figure)
    html_dict['stackbar_figure'] = stackbar_figure 

    index_page = os.path.join(html_output_directory, 'index.html')
    rendered = render_to_string('index.tmpl', html_dict)
    open(index_page, 'w').write(rendered.encode("utf-8"))

    return index_page

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
