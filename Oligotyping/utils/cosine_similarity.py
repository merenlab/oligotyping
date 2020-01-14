#-*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

# This program helps to agglomorate oligotypes based on their frequency
# patterns across samples, assuming that most of them are co-occurring
# and don't functionally somewhat anologous.
#
# Here is an example input file:
#
#     OLIGOTYPES        	DATASET_NAME1   DATASET_NAME2   DATASET_NAME3   DATASET_NAME4   DATASET_NAME5   DATASET_NAME6   DATASET_NAME7   DATASET_NAME8   DATASET_NAME1   
#     TAGATAAAAA-C-GTTAC	20.4509573117	2.31599654326	10.6854411688	18.141183236	2.75361300492	13.4206516811	12.5506081754	9.18451529079	10.497033588
#     TAGATAGAAA-C-GTTAC	12.7128172113	0.36018939857	15.0573756759	18.9353393378	0.458036304291	19.9503623429	13.3957224948	5.47629339621	13.6538638382
#     TAGAGAGAAA-C-GTTAC	18.0246980365	0.615297579484	14.7570747071	12.8939331142	0.36652632687	17.9981258947	14.3154911844	5.35317215042	15.6756810063
#     AAGATAGAAG-C-GCTAC	13.0686419929	0.878617519887	11.6350565111	19.9689321105	0.473692712916	17.4658459333	16.597647458	6.57207683968	13.3395099821
#     TGGAGAGGAC-T-TGCGT	0.475162605924	25.571525922	9.43118736779	1.15559558043	41.5513135103	2.97259932294	6.80887034513	1.44177620772	10.5919691377
#     TAGAGAAAATCT-AATAC	13.7423501569	0.640768984652	14.6270833878	10.7464043591	0.566629189838	19.4586949365	24.4534185641	5.40165386795	10.3629965532
#     CGGTTAGAAG-C-CCCTC	0.489691111705	29.1574841354	11.0672069019	1.51011971333	39.6255850137	1.87336134179	4.81021673323	0.518058382573	10.9482766663
#     GAGTTAGGAG-C-GCTAC	0.198166962328	25.8465099248	28.0636126996	1.03592770249	28.6071336884	7.38720775725	3.41467006542	0.900050098907	4.54672110072
#     TAGAGAGAAA-T-GTTAC	14.4873726757	0.264949578898	15.5784524907	13.5772103533	0.386489951693	17.1867602335	17.2095305909	2.24129536992	19.0679387553
#     TACTTAGAACTT-GGTAT	18.9421855944	1.3822243641	10.3208847966	14.7093902734	0.650727124538	8.22547255033	19.0420299574	12.9040938859	13.8229914534
#     TATAGTGGACCT-AGTAT	12.1090224002	0.26665440018	19.3617011303	15.1000904089	0.459966563452	20.2927274957	12.8589434862	4.85846741089	14.6924267042

import sys
from scipy import spatial

from Oligotyping.utils.utils import get_vectors_from_oligotypes_across_samples_matrix
from Oligotyping.visualization.oligotype_sets_distribution import vis_oligotype_sets_distribution

def cosine_distance(v1, v2):
    v1_abs, v2_abs = [], []
    for i in range(0, len(v1)):
        v1_abs.append(v1[i] * 100.0 / ((v1[i] + v2[i]) or 1))
        v2_abs.append(v2[i] * 100.0 / ((v1[i] + v2[i]) or 1))

    return spatial.distance.cosine(v1_abs, v2_abs)

def get_oligotype_sets_greedy(oligos, vectors, cosine_similarity_threshold, output_file = None):
    next_set_id = 0
    set_representatives = {}
    set_ids = {}
    len_oligos = len(oligos)
    
    for i in range(0, len_oligos):
        if i % 10 == 0 or i == len_oligos - 1:
            sys.stderr.write('\rOligo %d of %d :: num sets: %d' % (i, len_oligos, len(set_ids)))
            sys.stderr.flush()
        oligo = oligos[i]
        vector = vectors[oligo]
        
        shortest_distance_set_ID = None
        shortest_distance = sys.maxsize
        for set_representative in set_representatives:
            distance = cosine_distance(set_representatives[set_representative], vector)
            
            # this is a terrible thing to do, but when there are 1 million units
            # you really need to speed things up.
            if distance < 0.01:
                set_ids[set_representative].add(oligo)
                break
    
            if distance < shortest_distance:
                shortest_distance = distance
                shortest_distance_set_ID = set_representative
        
        if not shortest_distance_set_ID or not (shortest_distance < cosine_similarity_threshold):
            set_representatives['Set_%d' % next_set_id] = vector
            set_ids['Set_%d' % next_set_id] = set([oligo])
            next_set_id += 1
        else:
            set_ids[shortest_distance_set_ID].add(oligo)
        
    sys.stderr.write('\n')

    if output_file:
        f = open(output_file, 'w')
        for set_id in set_ids:
            f.write('%s\t%s\n' % (set_id, ','.join(set_ids[set_id])))
    f.close()

    return set_ids
 
 
def get_oligotype_sets(oligos, vectors, cosine_similarity_threshold, output_file = None):
    oligotype_sets = []
    distances = {}
    
    for i in range(0, len(oligos)):
        if oligos[i] not in distances:
            distances[oligos[i]] = {}
        for j in range(i, len(oligos)):
            if oligos[j] not in distances:
                distances[oligos[j]] = {}
            
            distances[oligos[i]][oligos[j]] = cosine_distance(vectors[oligos[i]], vectors[oligos[j]])
            distances[oligos[j]][oligos[i]] = cosine_distance(vectors[oligos[i]], vectors[oligos[j]])
    
    ids = list(range(0, len(oligos)))
    while 1:
        if not len(ids):
            break
        
        new_oligotype_set = [ids[0]]
        seed = oligos[ids[0]]
    
        for _id in ids[1:]:
            candidate = oligos[_id]
            if distances[seed][candidate] <= cosine_similarity_threshold:
                new_oligotype_set.append(_id)
    
        for _id in new_oligotype_set:
            ids.remove(_id)
    
        oligotype_sets.append(new_oligotype_set)


    oligotype_sets_final = [[oligos[i] for i in oligotype_set] for oligotype_set in oligotype_sets]

    if output_file:
        f = open(output_file, 'w')
        for i in range(0, len(oligotype_sets_final)):
            oligotype_set = oligotype_sets_final[i]
            f.write('Set_%d\t%s\n' % (i, ','.join(oligotype_set)))

    return oligotype_sets_final


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Group Oligotypes Based on Cosine Similarity')
    parser.add_argument('oligotypes_across_samples', metavar = 'OLIGOTYPES_ACROSS_DATASETS',\
                        help = 'A TAB-delimited matrix file that that contains normalized\
                                oligotype frequencies across samples. See the source code\
                                for an example format')
    parser.add_argument('-c', '--cosine-similarity-threshold', default = 0.1, type=float,\
                        metavar = 'COS_SIM_THRESHOLD', help = 'The higher the threshold is,\
                                the more oligotypes will be pulled together. Cosine similarity\
                                would return 0 for perfectly similar two vectors. Default is %(default)f.')
    parser.add_argument('-o', '--output-file', default = None, metavar = 'OUTPUT_FILE',\
                        help = 'Serialized list of lists object that contains oligotype sets')


    args = parser.parse_args()
 

    def get_samples():
        return [d.strip() for d in open(args.oligotypes_across_samples).readline().split('\t')[1:]]

    oligos, vectors = get_vectors_from_oligotypes_across_samples_matrix(args.oligotypes_across_samples)
    
    partitions = get_oligotype_sets(oligos, vectors, args.cosine_similarity_threshold, args.output_file)

    samples = get_samples()


    print('\n\t%d oligotypes split into %d partitions based on cosine similarity of %f. Here how they were distributed:\n'\
                        % (len(oligos), len(partitions), args.cosine_similarity_threshold))
    
    for partition in partitions:
        print('    - %s\n' % (', '.join(partition)))

    vis_oligotype_sets_distribution(partitions, vectors, samples, legend = True,\
        project_title = 'Cosine Similarity Threshold %.4f' % args.cosine_similarity_threshold)
