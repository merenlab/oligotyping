#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import inspect
import unittest

from Oligotyping.utils import blast as blast

my_path = os.path.realpath(os.path.abspath(os.path.split(inspect.getfile(inspect.currentframe()))[0]))

def files_are_the_same(file1, file2):
    lines1 = open(file1).readlines()
    lines2 = open(file2).readlines()

    if len(lines1) != len(lines2):
        return False

    for i in range(0, len(lines1)):
        if lines1[i] != lines2[i]:
            return False

    return True

class Tests(unittest.TestCase):
    def setUp(self):
        self.query  = os.path.join(my_path, 'files/unaligned-25K-illumina-test.fa')
        self.target = os.path.join(my_path, 'files/unaligned-25K-illumina-target.db')
        self.expected_output = os.path.join(my_path, 'files/unaligned-25K-illumina-result.b6')
        self.output = os.path.join(my_path, 'files/unaligned-25K-illumina-test.b6')
        self.min_percent_identity = 97.0

    def tearDown(self):
        pass
       
    def test_01_ParallelBlast(self):
        s = blast.LocalBLAST(self.query, 
                             self.target, 
                             self.output)
        s.make_blast_db()
        s.params = '-perc_identity %f' % self.min_percent_identity 
        s.search_parallel(num_processes = 2, num_reads_per_process = 1000)
        similarity_dict = s.get_results_dict(min_identity = self.min_percent_identity)

        self.assertTrue(len(similarity_dict) == 21076)

        self.assertTrue(files_are_the_same(self.expected_output, self.output))


    def test_99_CleanUp(self):
        for output in ['unaligned-25K-illumina-target.db.nhr',
                       'unaligned-25K-illumina-target.db.nin',
                       'unaligned-25K-illumina-target.db.nsq',
                       'unaligned-25K-illumina-test.b6']:
            os.remove(os.path.join(my_path, 'files/%s' % output))
        pass
