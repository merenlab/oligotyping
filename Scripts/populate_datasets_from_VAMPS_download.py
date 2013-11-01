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

import sys

import Oligotyping.lib.fastalib as u


fasta = u.SequenceSource(sys.argv[1])
taxon = sys.argv[2]

output = u.FastaOutput(sys.argv[2].replace(';', ''))

while fasta.next():
    if fasta.id.find(taxon) > -1:
        acc = fasta.id.split('|')[0]
        project = fasta.id.split('|')[1].split('=')[1]
        sample = fasta.id.split('|')[2].split('=')[1]
        new_id = project + '_' + sample + '_' + acc

        abundance = int(fasta.id.split('|')[7].split('=')[1])

        for i in range(0, abundance):
            output.write_id('%s-%s|%s' % (new_id, str(i), fasta.id))
            output.write_seq(fasta.seq, split = False)

fasta.close()
output.close()
