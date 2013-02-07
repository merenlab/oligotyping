# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

pretty_names = {'output_directory': 'Output directory',
                'project': 'Project',
                'run_date': 'Run date',
                'end_of_run': 'End of run',
                'version': 'Library version',
                'cmd_line': 'Command line',
                'info_file_path': 'Extraction info output file',
                'alignment': 'Input alignment',
                'entropy': 'Input entropy file',
                'total_seq': 'Number of sequences analyzed',
                'alignment_length': 'Number of characters in each alignment',
                'number_of_auto_components': 'Number of entropy components to be selected automatically',
                'number_of_selected_components': 'Number of entropy components chosen by the user',
                'num_reads_eliminated_due_to_min_base_quality': 'Number of reads eliminated due to --min-base-quality',
                'num_oligos_after_l_elim': 'Number of oligotypes left after --limit-oligotypes-to elimination',
                'num_oligos_after_e_elim': 'Number of oligotypes left after --exclude-oligotypes elimination',
                's': 'Min number of datasets oligotype expected to appear',
                'a': 'Min % abundance of oligotype in at least one dataset',
                'A': 'Min total abundance of oligotype in all datasets',
                'M': 'Min substantive abundance of an oligotype (-M)',
                'T': 'Cosine similarty threshold to generate oligotype sets',
                'q': 'Min PHRED score for each base of interes in every read',
                'm': 'Min entropy for a component to be picked for decomposition',
                'd': 'Max number of discriminants to use for decomposition',
                'maximum_variation_allowed': 'Maximum variation allowed in each node (-V)',
                'average_read_length': 'Average read length (without gaps)',
                'quals_provided': 'Quality scores were provided',
                'blast_ref_db_provided': 'Reference DB were provided for local BLAST search',
                'blast_ref_db': 'Reference DB for local BLAST search',
                'limit_oligotypes_to': 'Discarded all other oligotypes except',
                'exclude_oligotypes': 'Oligotypes excluded from the analysis',
                'bases_of_interest_locs': 'Base locations of interest in the alignment',
                'num_datasets_in_fasta': 'Number of datasets found',
                'num_unique_oligos': 'Number of unique oligotypes (raw)',
                'num_oligos_after_s_elim': 'Oligotypes after "min number of datasets" elimination',
                'num_oligos_after_a_elim': 'Oligotypes after "min % abundance in a dataset" elimination',
                'num_oligos_after_A_elim': 'Oligotypes after "min total abundance (-A)" elimination',
                'num_oligos_after_M_elim': 'Oligotypes after "min substantive abundance (-M)" elimination',
                'num_sequences_after_qc': 'Number of sequences represented after quality filtering',
                'datasets_removed_after_qc': 'Datasets removed for having 0 oligotypes left after filtering',
                'oligos_fasta_file_path': 'FASTA file for abundant oligotypes',
                'representative_seqs_fasta_file_path': 'Representative sequences per oligotype',
                'oligos_nexus_file_path': 'NEXUS file for abundant oligotypes',
                'environment_file_path': 'Environment file',
                'matrix_count_file_path': 'Sample/oligotype abundance data matrix (counts)',
                'matrix_percent_file_path': 'Sample/oligotype abundance data matrix (percents)',
                'generate_sets': 'Oligotype sets were requested to be generated',
                'matrix_count_oligo_sets_file_path': 'Abundance data matrix for oligotype sets (counts)',
                'matrix_percent_oligo_sets_file_path': 'Abundance data matrix for oligotype sets (percents)',
                'across_datasets_MN_file_path': 'Oligotypes across datasets matrix (MAX normalized)',
                'across_datasets_SN_file_path': 'Oligotypes across datasets matrix (SUM normalized)',
                'oligotype_sets_file_path': 'Groups of oligotypes',
                'oligotype_sets_info': 'Sets of oligotypes based on frequency patterns',
                'oligotype_sets_figure_path': 'Oligotype sets figure',
                'oligos_across_datasets_file_path': 'Oligotypes across datasets figure',
                'output_directory_for_reps': 'Representative sequences for oligotypes directory',
                'colors_file_path': 'Random colors for oligotypes',
                'stack_bar_file_path': 'Oligotype distribution stack-bar figure',
                'stack_bar_with_agglomerated_oligos_file_path': 'Stack-bar figure with oligotype sets',
                'output_directory_for_datasets': 'Directory to store datasets associted information',
                'removed_min_substantive_abundance_reason': 'Outliers removed due to -M',
                'removed_maximum_variation_allowed_reason': 'Outliers removed due to -V',
                'removed_outliers_total': 'Total number of outliers removed during the refinement',
                'relocated_min_substantive_abundance_reason': 'Relocated outliers originally removed due to -M',
                'relocated_maximum_variation_allowed_reason': 'Relocated outliers originally removed due to -V',
                'relocated_outliers_total': 'Total number of relocated outliers',
                'final_min_substantive_abundance_reason': 'Final number of outliers due to -M',
                'final_maximum_variation_allowed_reason': 'Final number of outliers due to -V',
                'final_outliers_total': 'Final total number of outliers',
                'num_raw_nodes': 'Number of raw nodes (before the refinement)',
                'num_final_nodes': 'Number of final nodes (after the refinement)',
                'skip_agglomerating_nodes': 'Skip agglomerating nodes',
                'merge_homopolymer_splits': 'Merge homopolymer splits',
                'skip_removing_outliers': 'Skip removing outliers',
                'relocate_outliers': 'Try to relocate outliers',
                'read_distribution_table_path': 'Read distribution among datasets table',
                'node_representatives_file_path': 'Representative sequences per node',
                'sample_mapping': 'Mapping file',
                'skip_basic_analyses': 'Skip performing basic analyses'
                }
