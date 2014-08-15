#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

# takes an environment file and a generates matching percent and count matrices.

import sys

from Oligotyping.utils.utils import get_samples_dict_from_environment_file
from Oligotyping.utils.utils import get_oligos_sorted_by_abundance
from Oligotyping.utils.utils import get_units_across_samples_dicts
from Oligotyping.utils.utils import get_unit_counts_and_percents
from Oligotyping.utils.utils import generate_MATRIX_files

samples_dict = get_samples_dict_from_environment_file(sys.argv[1])
oligos = get_oligos_sorted_by_abundance(samples_dict)
oligos.reverse()
unit_counts, unit_percents = get_unit_counts_and_percents(oligos, samples_dict)
samples = sorted(samples_dict.keys())

generate_MATRIX_files(oligos, samples, unit_counts, unit_percents, sys.argv[1] + '-MATRIX-COUNT',  sys.argv[1] + '-MATRIX-PERCENT')
