# -*- coding: utf-8 -*-
# takes an environment file and a generates matching percent and count matrices.

import sys
from Oligotyping.utils.utils import get_samples_dict_from_environment_file
from Oligotyping.utils.utils import get_oligos_sorted_by_abundance
from Oligotyping.utils.utils import get_units_across_samples_dicts
from Oligotyping.utils.utils import get_unit_counts_and_percents
from Oligotyping.utils.utils import generate_MATRIX_files

samples_dict = get_samples_dict_from_environment_file(sys.argv[1])
oligos = get_oligos_sorted_by_abundance(samples_dict)
unit_counts, unit_percents = get_unit_counts_and_percents(oligos, samples_dict)
samples = samples_dict.keys()

generate_MATRIX_files(oligos, samples, unit_counts, unit_percents, sys.argv[1] + '-MATRIX-COUNT',  sys.argv[1] + '-MATRIX-PERCENT')
