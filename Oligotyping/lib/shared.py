#!/usr/bin/python
# -*- coding: utf-8

# Copyright (C) 2013, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

#
# some shared functions between decomposition and oligotyping, that require
# 'self' to be passed.
#

import os

from Oligotyping.utils.utils import run_command
from Oligotyping.utils.utils import store_filtered_matrix
from Oligotyping.utils.utils import get_sample_mapping_dict
from Oligotyping.utils.utils import get_temporary_file_name


def generate_default_figures(_object):
    figures_dict = {}
    figures_dict['basic_analyses'] = {}
    figures_dict['basic_reports'] = {}

    #
    # basic reports
    #

    for (analysis, script, output_dir) in [('Stackbar', 'o-stackbar.R', 'stackbar')]:
        #¬†generating a stackbar is only feasible for oligotyping:
        if _object.analysis != 'oligotyping':
            continue

        figures_dict['basic_reports'][output_dir] = {}

        target_dir = _object.generate_output_destination('%s/__default__/%s' \
                                                            % (os.path.basename(_object.figures_directory), output_dir),
                                                      directory = True)
            
        output_prefix = os.path.join(target_dir, output_dir)
        cmd_line = ('%s "%s" --title "%s" -o "%s" --colors_file "%s" >> "%s" 2>&1' \
                                          % (script,
                                             _object.environment_file_path,
                                             _object.project,
                                             output_prefix,
                                             _object.colors_file_path,
                                             _object.log_file_path))
        _object.progress.update('%s ...' % (analysis))
        _object.logger.info('figure basic_reports: %s' % (cmd_line))
        run_command(cmd_line)
        figures_dict['basic_reports'][output_dir][output_dir] = output_prefix


    for (analysis, script, output_dir) in [('Read Distribution Lines', 'o-lines-for-each-column.R', 'lines'),
                                           ('Read Distribution Bars', 'o-bars-for-each-column.R', 'bars')]:
        figures_dict['basic_reports'][output_dir] = {}

        target_dir = _object.generate_output_destination('%s/__default__/%s' \
                                                            % (os.path.basename(_object.figures_directory), output_dir),
                                                      directory = True)
            
        output_prefix = os.path.join(target_dir, output_dir)
        cmd_line = ('%s "%s" "%s" >> "%s" 2>&1' % (script,
                                             _object.read_distribution_table_path,
                                             output_prefix,
                                             _object.log_file_path))
        _object.progress.update('%s ...' % (analysis))
        _object.logger.info('figure basic_reports: %s' % (cmd_line))
        run_command(cmd_line)
        figures_dict['basic_reports'][output_dir][output_dir] = output_prefix

        
    #
    # basic analyses
    #
    if not _object.skip_basic_analyses:
        for (analysis, script, output_dir) in [('Cluster Analysis', 'o-cluster-analysis.R', 'cluster_analysis'),
                                               ('NMDS Analysis', 'o-metaMDS-analysis.R', 'nmds_analysis')]:
            figures_dict['basic_analyses'][output_dir] = {}
                        
            target_dir = _object.generate_output_destination('%s/__default__/%s' \
                                                                % (os.path.basename(_object.figures_directory), output_dir),
                                                          directory = True)
                
            for (distance_metric, matrix_file) in [("canberra", _object.matrix_percent_file_path),
                                                   ("kulczynski", _object.matrix_percent_file_path),
                                                   ("jaccard", _object.matrix_percent_file_path),
                                                   ("horn", _object.matrix_percent_file_path),
                                                   ("bray", _object.matrix_percent_file_path)]:
                output_prefix = os.path.join(target_dir, distance_metric)
                cmd_line = ('%s "%s" %s "%s" "%s" >> "%s" 2>&1' % 
                                        (script,
                                         matrix_file,
                                         distance_metric,
                                         _object.project,
                                         output_prefix,
                                         _object.log_file_path))
                _object.progress.update('%s "%s" ...' % (analysis, distance_metric))
                _object.logger.info('figure basic_analyses: %s' % (cmd_line))
                run_command(cmd_line)
                figures_dict['basic_analyses'][output_dir][distance_metric] = output_prefix
            
    return figures_dict


def generate_exclusive_figures(_object):
    exclusive_figures_dict = {}

    sample_mapping_dict = get_sample_mapping_dict(_object.sample_mapping)
        
    for category in sample_mapping_dict:
        exclusive_figures_dict[category] = {}
        samples = list(sample_mapping_dict[category].keys())
            
        # double filter: first makes sure sample was not removed from the analysis due to losing all its reads during the
        #¬†refinement, second makes sure that sample was actually mapped to something in the sample mapping file.
        samples = [s for s in [s for s in samples if s in _object.samples] if sample_mapping_dict[category][s]]
        samples.sort()

        mapping_file_path = get_temporary_file_name('%s-' % category, '-mapping.txt', _object.tmp_directory)
        mapping_file = open(mapping_file_path, 'w')
        mapping_file.write('samples\t%s\n' % (category))
            
        for sample in samples:
            mapping_file.write('%s\t%s\n' % (sample, sample_mapping_dict[category][sample]))
        mapping_file.close()

        if samples == _object.samples:
            matrix_percent_path = _object.matrix_percent_file_path
            matrix_count_path = _object.matrix_count_file_path
        else:
            matrix_percent_path = get_temporary_file_name('%s-' % category, '-matrix-percent.txt', _object.tmp_directory)
            matrix_count_path = get_temporary_file_name('%s-' % category, '-matrix-count.txt', _object.tmp_directory)

            if store_filtered_matrix(_object.matrix_percent_file_path, matrix_percent_path, samples) < 3:
                _object.logger.info("skipping exclusive figs for '%s'; less than 3 samples were left in MP"\
                                         % (category))
                continue
            if store_filtered_matrix(_object.matrix_count_file_path, matrix_count_path, samples) < 3:
                _object.logger.info("skipping exclusive figs for '%s'; less than 3 samples were left in MC"\
                                         % (category))
                continue

        # ready to roll.
        _object.logger.info("exclusive figs for '%s' with %d samples; mapping: '%s', MP: '%s', MC: '%s'"\
                             % (category, len(samples), mapping_file_path, matrix_percent_path, matrix_count_path))


        for (analysis, script, output_dir) in [('NMDS Analysis', 'o-metaMDS-analysis-with-metadata.R', 'nmds_analysis')]:
            exclusive_figures_dict[category][output_dir] = {}
                        
            target_dir = _object.generate_output_destination('%s/%s/%s' % (os.path.basename(_object.figures_directory),
                                                                        category,
                                                                        output_dir),
                                                                        directory = True)
                
            for (distance_metric, matrix_file) in [("canberra", matrix_percent_path),
                                                   ("kulczynski", matrix_percent_path),
                                                   ("jaccard", matrix_percent_path),
                                                   ("horn", matrix_percent_path),
                                                   ("bray", matrix_percent_path)]:
                output_prefix = os.path.join(target_dir, distance_metric)
                cmd_line = ('%s -o "%s" -d "%s" -m "%s" --title "%s" "%s" "%s" >> "%s" 2>&1' % 
                                        (script,
                                         output_prefix,
                                         distance_metric,
                                         category,
                                         _object.project,
                                         matrix_file,
                                         mapping_file_path,
                                         _object.log_file_path))
                _object.progress.update('%s "%s" for "%s" ...' % (analysis, distance_metric, category))
                _object.logger.info('exclusive figure: %s' % (cmd_line))
                run_command(cmd_line)
                exclusive_figures_dict[category][output_dir][distance_metric] = output_prefix


        # heatmap
        for (analysis, script, output_dir) in [('Heatmap Analysis', 'o-heatmap.R', 'heatmap_analysis')]:
            exclusive_figures_dict[category][output_dir] = {}
                        
            target_dir = _object.generate_output_destination('%s/%s/%s' % (os.path.basename(_object.figures_directory),
                                                                        category,
                                                                        output_dir),
                                                                        directory = True)
                
            for (distance_metric, matrix_file) in [("canberra", matrix_percent_path),
                                                   ("kulczynski", matrix_percent_path),
                                                   ("jaccard", matrix_percent_path),
                                                   ("horn", matrix_percent_path),
                                                   ("bray", matrix_percent_path)]:
                output_prefix = os.path.join(target_dir, distance_metric)
                cmd_line = ('%s "%s" -m "%s" -d %s --title "%s" -o "%s" >> "%s" 2>&1' % 
                                        (script,
                                         matrix_file,
                                         mapping_file_path,
                                         distance_metric,
                                         _object.project,
                                         output_prefix,
                                         _object.log_file_path))
                _object.progress.update('%s "%s" for "%s" ...' % (analysis, distance_metric, category))
                _object.logger.info('exclusive figure: %s' % (cmd_line))
                run_command(cmd_line)
                exclusive_figures_dict[category][output_dir][distance_metric] = output_prefix

    return exclusive_figures_dict


if __name__ == '__main__':
    pass
