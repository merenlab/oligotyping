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

absolute = os.path.join(os.path.dirname(os.path.realpath(__file__)))


def generate_html_output(run_info_dict, output_directory = None):
    if not output_directory:    
        output_directory = os.getcwd()
    
    html_dict = copy.deepcopy(run_info_dict)

    shutil.copy2(os.path.join(absolute, 'static/style.css'), os.path.join(output_directory, 'style.css'))
    shutil.copy2(os.path.join(absolute, 'static/header.png'), os.path.join(output_directory, 'header.png'))

    # embarrassingly ad-hoc:
    entropy_figure = os.path.join(output_directory, 'entropy.png')
    shutil.copy2(os.path.join(run_info_dict['entropy'].split('.')[0] + '.png'), entropy_figure)
    html_dict['entropy_figure'] = entropy_figure

    stackbar_figure = os.path.join(output_directory, 'stackbar.png')
    shutil.copy2(os.path.join(run_info_dict['stack_bar_file_path']), stackbar_figure)
    html_dict['stackbar_figure'] = stackbar_figure 

    info_section = open(os.path.join(absolute, 'templates/info.tmpl')).read() % html_dict
    html_dict['info_section'] = info_section

    base = open(os.path.join(absolute, 'templates/base.tmpl')).read() % html_dict

    index_page = os.path.join(output_directory, 'index.html')

    open(index_page, 'w').write(base)

    return index_page

if __name__ == '__main__':
    pass
