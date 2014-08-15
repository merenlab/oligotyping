#!/usr/bin/env python
# -*- coding: utf-8

# Copyright (C) 2013, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

import sys

matrix_counts = open(sys.argv[1])
environment = open(sys.argv[1] + '-ENV', 'w')

units = matrix_counts.readline().strip().split('\t')[1:]

for line in matrix_counts.readlines():
    fields = line.strip().split('\t')
    sample = fields[0]
    for i in range(0, len(units)):
        unit = units[i]
        value = int(fields[i + 1])
        if value:
            environment.write('%s\t%s\t%d\n' % (unit, sample, value))

matrix_counts.close()
environment.close()
