#!/usr/bin/python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import argparse

def decomposer():
    parser = argparse.ArgumentParser(description='Minimum Entropy Decomposition')
    parser.add_argument('alignment', metavar = 'FILEPATH',
                        help = 'Alignment file that contains all samples and sequences in FASTA format')
    parser.add_argument('-m', '--min-entropy', type=float, default=0.0965, metavar="FLOAT",
                        help = 'Minimum entropy for a component to have in order to be picked as a\
                                discriminant. Defeault: %(default)f')
    parser.add_argument('-X', '--skip-m-normalization', action = 'store_true', default = False,
                        help = 'Skip the m normalization heuristics. Normalization heuristics normalize --min-entropy\
                                parameter per node before decomposition depending on the node size, and the most\
                                abundant sequence in the dataset. Please see the documentation for more\
                                comprehensive description and best practices. Default: %(default)s')
    parser.add_argument('-d', '--number-of-discriminants', type=int, default=4, metavar="INTEGER",
                        help = 'Number of discriminant locations to be used for entropy decomposition\
                                discriminant. Defeault: %(default)d')
    parser.add_argument('-A', '--min-actual-abundance', type=int, default=0, metavar="INTEGER",
                        help = 'Minimum number of reads in a node for decomposition to continue. Decomposition\
                                will continue for any node that has more reads than this number as far as they\
                                present an entropy that is larger than --min-entropy. This number should be\
                                chosen carefully depending on the size of the sample. Althought this parameter\
                                is available to you for historical reasons, For noise filtering, you are\
                                encouraged to use --min-substantive-abundance parameter instead.')
    parser.add_argument('-M', '--min-substantive-abundance', type=int, default=0, metavar = "INTEGER",
                        help = 'Unlike "actual" abundance, "substantive" abundance is interested in the abundance\
                                of the most abundant read in a node. If the abundance of the most abundant\
                                unique sequence in a node is smaller than the number given with this parameter\
                                the node will be eliminated and not included in downstream analyses. This is the\
                                most appropriate, and most cases, the only noise filtering parameter. If the user\
                                does not set a value for minimum substantive abudnance, MED algorithm will set\
                                one by default by dividing the number of reads in the input dataset by 5,000.')
    parser.add_argument('-V', '--maximum-variation-allowed', type=int, default=None, metavar = 'INTEGER',
                        help = 'This parameter is being used to remove "outliers" from nodes. The similarity of a\
                                read in a node is less than --maximum-variation-allowed than the representative\
                                sequence of the node, it is identified as an outlier. If not set, this value is \
                                being computed depending on the average read length.')
    parser.add_argument('-t', '--sample-name-separator', type=str, default='_', metavar = "CHAR",
                        help = 'Character that separates sample name from unique info in the defline. For insatnce\
                                if the defline says >sample-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-p', '--project', default = None, type=str, metavar = "STR",
                        help = 'When a project name is set, given name will be used in figures whenever possible.')
    parser.add_argument('-g', '--generate-frequency-curves', action = 'store_true', default = False,
                        help = 'When set, figure with frequency curve for unique reads and entropy distribution\
                                will be generated for each node. Depending on the number of nodes, this might\
                                be a time consuming step.')
    parser.add_argument('-S', '--skip-removing-outliers', action = 'store_true', default = False,
                        help = 'When set, outliers will not be removed from nodes.')
    parser.add_argument('-H', '--merge-homopolymer-splits', action = 'store_true', default = False,
                        help = 'When set, nodes that differ from each other by only one nucleotide that happens\
                                to be observed as an insertion at the upstream or downstream of a homopolymer\
                                region will be merged.')
    parser.add_argument('-R', '--relocate-outliers', action = 'store_true', default = False,
                        help = 'Outliers are identified in two places: (1) during the raw topology computation\
                                and (2) during the refinement step where distant reads are removed from nodes.\
                                This parameter, when set, makes the pipeline go through each read identified as\
                                an outlier and try to find the best nodes for them. Please read the documentation\
                                for details. This step might take a long time. Default: %(default)s')
    parser.add_argument('-F', '--store-topology-dict', action = 'store_true', default = False,
                        help = 'When set, topology dict with read ids will be generated. This may take a very large\
                                disk space and computation time for large data sets')
    parser.add_argument('-K', '--keep-tmp', action = 'store_true', default = False,
                        help = 'When set, directory with temporary BLAST results will not be deleted at the end of the\
                                run. It may be necessary to debug the results')
    parser.add_argument('-T', '--no-threading', action = 'store_true', default = False,
                        help = 'When set, decomposer does not spawn multiple threads. Default behavior is\
                                multi-threaded.')
    parser.add_argument('-N', '--number-of-threads', type=int, default = None, metavar = "INTEGER",
                        help = 'Number of threads to use. It is a good idea to keep this number smaller than the number\
                                of CPU cores available. If not set, this number will be set to 90%% of available cores,\
                                or (available cores - 1) if 10%% of the cores is a number smaller than 1')    
    parser.add_argument('-E', '--sample-mapping', metavar = 'FILEPATH', default = None,
                        help = 'TAB delimited categorical mapping of samples to be used for post-analysis\
                                visualizations. Refer to the tutorial for the file format')
    parser.add_argument('--gen-html', action = 'store_true', default = False,
                        help = 'When set, decomposer will generate a static HTML output to browse analysis results')
    parser.add_argument('--skip-gen-figures', action = 'store_true', default = False,
                        help = 'When set, decomposer will not attempt to generate figures post analysis')
    parser.add_argument('--skip-check-input-file', action = 'store_true', default = False,
                        help = 'When set, input FASTA will not be checked for potential errors')
    parser.add_argument('--skip-gexf-files', action = 'store_true', default = False,
                        help = 'When set, GEXF files for network and topology will not be generated')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'When set, the pipeline will do only the essential steps, skipping anything\
                                auxiliary, even if other parameters require otherwise. Please do not use it other than\
                                benchmarking or testing purposes')


    return parser

def oligotyping():
    parser = argparse.ArgumentParser(description='Start an Oligotyping Process')
    parser.add_argument('alignment', metavar = 'INPUT ALIGNMENT',
                        help = 'Alignment file that contains all samples and sequences in FASTA format')
    parser.add_argument('entropy', metavar = 'ENTROPY',
                        help = 'File that contains the columns and the entropy values computer previously')
    parser.add_argument('-o', '--output-directory', help = 'Output directory', default = None)
    parser.add_argument('-c', '--number-of-auto-components', type=int, default=None,
                        help = 'Number of components to use from alignment to generate oligotypes. Default\
                                is "5", which is a completely arbitrary value. Number of components should\
                                be determined after a careful examination of entropy figure.')
    parser.add_argument('--qual-scores-file', metavar = 'QUAL SCORES FILE',
                        help = 'FASTA formatted file that contains PHRED base call values\
                                for each read in the alignment file')
    parser.add_argument('--qual-scores-dict', metavar = 'QUAL SCORES DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                                PHRED base call values for each read in the alignment file. If you\
                                provide --qual-scores-file, that file will be used to recompute this\
                                dictionary and the file you refer with this parameter will\
                                not be ignored')
    parser.add_argument('--qual-stats-dict', metavar = 'QUAL STATS DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                                PHRED base call quality score statistics for the alignment file. If\
                                you provide --qual-scores-dict, it will be used to recompute this\
                                dictionary and the file you refer to with this parameter will\
                                actually not be used')
    parser.add_argument('-q', '--min-base-quality', type=int, default=15,
                        help = 'Minimum quality score for each base in locations of interest of a read to be\
                                considered in an oligotype. When base quality score files are provided, this\
                                value makes sure that low quality bases that are more likely to be the result\
                                of random sequencing errors do not create artificial oligotypes. Any read that has\
                                less quality score than the given value, will simply be discarded. This parameter\
                                only in effect when --qual-scores-file or --qual-scores-dict parameters are used. \
                                Defeault is %(default)d.')
    parser.add_argument('-C', '--selected-components', type=str, default=None,
                        help = 'Comma separated entropy components to be used during the oligotyping process.')
    parser.add_argument('-s', '--min-number-of-samples', type=int, default=1,
                        help = 'Minimum number of samples oligotype expected to appear. The deafult is "5", which\
                                is another completely arbitrary value. This parameter should be defined based\
                                on the number of samples included in the analysis. If there are 10 samples,\
                                3 might be a good choice, if there are 5 samples, 1 would be a better one\
                                depending on the study. Default is %(default)d.')
    parser.add_argument('-a', '--min-percent-abundance', type=float, default=0.0,
                        help = 'Minimum percent abundance of an oligotype in at least one sample. The default\
                                is "0.0". Just like --min-number-of-samples parameter, this parameter too is\
                                to eliminate oligotypes that are formed by sequencing errors occured at the\
                                component of interest. The value should be decided based on the average number\
                                of sequences every sample has.')
    parser.add_argument('-A', '--min-actual-abundance', type=int, default=0,
                        help = 'Minimum total abundance of an oligotype in all datastes. The default\
                                is "0". If the total abundance of an oligotype is smaller than the number given\
                                with this parameter, oligotype would be eliminated and not included in downstream\
                                analyses. Default is %(default)d.')
    parser.add_argument('-M', '--min-substantive-abundance', type=int, default=0,
                        help = 'Unlike "actual" abundance, "substantive" abundance is interested in the abundance\
                                of the most abundant read in an oligotype. If the abundance of the most abundant\
                                unique sequence in an oligotype smaller than the number given with this parameter\
                                the oligotype will be eliminated and not included in downstream analyses. Default\
                                is %(default)d.')
    parser.add_argument('-t', '--sample-name-separator', type=str, default='_',
                        help = 'Character that separates sample name from unique info in the defline. For insatnce\
                                if the defline says >sample-1_GD7BRW402IVMZE, the separator should be set to "_"\
                                (which is the default character).')
    parser.add_argument('-l', '--limit-representative-sequences', type=int, default=None,
                        help = 'At the end of the oligotyping sequences that are being represented by the same\
                                oligotype are being uniqued and stored in separate files. The number of sequences\
                                to keep from the frequency ordered list can be defined with this parameter (e.g.\
                                -l 10 would make it possible that only first 10 sequence would be stored). Default\
                                is 0, which stores everything, but when the sample size is too big, this could\
                                take up disk space.')
    parser.add_argument('--limit-oligotypes-to', type = str, default = None,
                        help = 'Comma separated list of oligotypes to be taken into account during the analysis.\
                                All other oligotypes will be discarded if a list of oligotypes is being speficied\
                                with this parameter.')
    parser.add_argument('-e', '--exclude-oligotypes', type = str, default = None,
                        help = 'Comma separated list of oligotypes to be excluded from the the analysis.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'Some relatively insignificant parts of the analysis may take a lot of time, such as\
                                generating figures for representative sequences. When this parameter is set, all\
                                trivial steps would be skipped to give results as soon as possible.')
    parser.add_argument('--no-figures', action = 'store_true', default = False,
                        help = 'When set, no figures will be generated or displayed.')
    parser.add_argument('--blast-ref-db', default = None, type=str,
                        help = 'When set, BLAST search will be done locally against the ref db (local BLAST search\
                                requires NCBI+ tools)')
    parser.add_argument('--colors-list-file', default = None, type=str,
                        help = 'Optional file that contains HTML color codes in each line to color oligotypes. Number\
                                of colors in the file has to be equal or greater than the number of abundant\
                                oligotypes, for which colors are going to be used for.')
    parser.add_argument('--skip-blast-search', action = 'store_true', default = False,
                        help = 'When set, BLAST search step will not be performed.')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                        help = 'When set, no figures will be shown.')
    parser.add_argument('--gen-html', action = 'store_true', default = False,
                        help = 'Generate static HTML output to browse analysis results.')
    parser.add_argument('--generate-sets', action = 'store_true', default = False,
                        help = 'Agglomerate oligotypes into oligotype sets when their frequency patterns across\
                        samples are similar. Oligotype sets simply put oligotypes into the same set if they \
                        co-occur in samples consistenly.')
    parser.add_argument('-K', '--keep-tmp', action = 'store_true', default = False,
                        help = 'When set, directory with temporary results will not be deleted at the end of the\
                                run. It may be necessary to debug the results')
    parser.add_argument('-S', '--cosine-similarity-threshold', default = 0.1, type=float, metavar='COS_SIM_TR',\
                        help = 'This value is used to agglomerate oligotypes into higher order groups. The higher\
                                the threshold is, the more oligotypes will be pulled together. Cosine similarity\
                                would return 0 for perfectly similar two vectors. Default is %(default)f.')
    parser.add_argument('--gen-sample-oligo-networks', action = 'store_true', default = False,
                        help = 'Generate oligotype network structure figures for each sample.')
    parser.add_argument('-E', '--sample-mapping', metavar = 'FILEPATH', default = None,
                        help = 'TAB delimited categorical mapping of samples to be used for post-analysis\
                                visualizations. Refer to the tutorial for the file format')
    parser.add_argument('--project', default = None, type=str,
                        help = 'When a project name is set, given name will be used in figures whenever possible.')
    parser.add_argument('--skip-check-input-file', action = 'store_true', default = False,
                        help = 'When set, input FASTA will not be checked for potential errors')
    parser.add_argument('--skip-basic-analyses', action = 'store_true', default = False,
                        help = 'When set, basic analyses, such as basic NMDS plots and clustering, will be\
                                skipped')
    parser.add_argument('--skip-gexf-network-file', action = 'store_true', default = False,
                        help = 'When set, GEXF network file will not be generated')
    parser.add_argument('-T', '--no-threading', action = 'store_true', default = False,
                        help = 'When set, oligotyping will not spawn multiple threads. Default behavior is\
                                multi-threaded whenever possible.')
    parser.add_argument('-N', '--number-of-threads', type=int, default = None, metavar = "INTEGER",
                        help = 'Number of threads to use. It is a good idea to keep this number smaller than the number\
                                of CPU cores available. If not set, this number will be set to 90%% of available cores,\
                                or (available cores - 1) if 10%% of the cores is a number smaller than 1')    

    return parser


def entropy():
    parser = argparse.ArgumentParser(description='Entropy Analysis')
    parser.add_argument('alignment', metavar = 'ALIGNMENT', help = 'Alignment file\
                         that contains all samples and sequences in FASTA format')
    parser.add_argument('--qual-scores-file', metavar = 'QUAL SCORES FILE',
                        help = 'FASTA formatted file that contains PHRED base call values\
                         for each read in the alignment file')
    parser.add_argument('--qual-scores-dict', metavar = 'QUAL SCORES DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                        PHRED base call values for each read in the alignment file. If you\
                        provide --qual-scores-file, that file will be used to recompute this\
                        dictionary and the file you refer with this parameter will\
                        not be ignored')
    parser.add_argument('--qual-stats-dict', metavar = 'QUAL STATS DICT',
                        help = 'Previously computed and serialized dictionary that contains\
                        PHRED base call quality score statistics for the alignment file. If\
                        you provide --qual-scores-dict, it will be used to recompute this\
                        dictionary and the file you refer to with this parameter will\
                        actually not be used')
    parser.add_argument('--uniqued', action = 'store_true', default = False,
                        help = 'When set, entropy computation will assume that the reads\
                        in FASTA file are unique. Frequency information of unique reads\
                        must be stored in the deflines. Every defline in the FASTA file\
                        must present the frequency information in this format:\
                        "freq:NUMBER", e.g. ">Read_ID|X|Y|freq:42", or ">Read_ID|freq:42|X|Y"')
    parser.add_argument('--weighted', action = 'store_true', default = False,
                        help = 'When set, entropy computation per column will use\
                        mean quality score for each column.')
    parser.add_argument('--amino-acid-sequences', action = 'store_true', default = False,
                        help = 'If sequences are composed of amino acids, instead of\
                                nucleotides.')
    parser.add_argument('--quick', action = 'store_true', default = False,
                        help = 'When set, entropy values will be shown as fast as\
                                possible (some visualization steps will be skipped).')
    parser.add_argument('--no-display', action = 'store_true', default = False,
                                help = 'When set, no figures will be shown.')

    return parser
