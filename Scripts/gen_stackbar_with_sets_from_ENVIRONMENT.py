# -*- coding: utf-8 -*-
# takes an environment file and a cosine similarity threshold as a parameter,
# generates an environment file with sets of units defined by the similarity
# they possess with respect to the frequency distribution patterns among
# samples.

import sys
from Oligotyping.utils.utils import get_samples_dict_from_environment_file
from Oligotyping.utils.utils import get_oligos_sorted_by_abundance
from Oligotyping.utils.utils import get_units_across_samples_dicts
from Oligotyping.utils.utils import get_unit_counts_and_percents
from Oligotyping.utils.cosine_similarity import get_oligotype_sets
from Oligotyping.visualization.oligotype_distribution_stack_bar import oligotype_distribution_stack_bar
from Oligotyping.utils.utils import generate_ENVIRONMENT_file 

samples_dict = get_samples_dict_from_environment_file(sys.argv[1])
oligos = get_oligos_sorted_by_abundance(samples_dict)
unit_counts, unit_percents = get_unit_counts_and_percents(oligos, samples_dict)
samples = samples_dict.keys()

across_samples_sum_normalized, across_samples_max_normalized = get_units_across_samples_dicts(oligos, samples_dict.keys(), unit_percents) 
oligotype_sets = get_oligotype_sets(oligos,
                                    across_samples_sum_normalized,
                                    float(sys.argv[2]))

print '%d sets from %d units' % (len(oligotype_sets), len(oligos))

oligotype_set_ids = range(0, len(oligotype_sets))

samples_dict_with_agglomerated_oligos = {}

for sample in samples:
    samples_dict_with_agglomerated_oligos[sample] = {}

for set_id in oligotype_set_ids:
    oligotype_set = oligotype_sets[set_id]
    for sample in samples:
        samples_dict_with_agglomerated_oligos[sample][set_id] = 0
        for oligo in samples_dict[sample]:
            if oligo in oligotype_set:
                samples_dict_with_agglomerated_oligos[sample][set_id] += samples_dict[sample][oligo]
    
    print 'set %s: %s' % (set_id, ', '.join(oligotype_set))

oligotype_distribution_stack_bar(samples_dict_with_agglomerated_oligos, None)
generate_ENVIRONMENT_file(samples,
                          samples_dict_with_agglomerated_oligos,
                          sys.argv[1] + '-cos-%s-SETS-ENVIRON' % sys.argv[2])
