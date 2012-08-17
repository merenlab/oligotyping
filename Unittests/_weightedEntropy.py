#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import inspect
import unittest

from Oligotyping.lib.entropy import entropy_analysis
from Oligotyping.utils.utils import get_quals_dict
from Oligotyping.utils.utils import get_qual_stats_dict

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
        self.output_directory_path = os.path.join(my_path, 'test-weighted-entropy')
        if os.path.exists(self.output_directory_path):
            shutil.rmtree(self.output_directory_path)
        os.makedirs(self.output_directory_path)
        self.alignment = os.path.join(my_path, 'files/500-V6V4-Pelagibacter.fasta')
        self.qual_scores_file = os.path.join(my_path, 'files/500-V6V4-Pelagibacter.qual')
        self.expected_result = os.path.join(my_path, 'files/500-V6V4-Pelagibacter-WEIGHTED-ENTROPY')

    def tearDown(self):
        pass

    def test_01_RunWeightedEntropy(self):
        output_file = os.path.join(self.output_directory_path, 'entropy.txt')
        QD = get_quals_dict(self.qual_scores_file, self.alignment, output_file_path = os.path.join(self.output_directory_path, 'QUALS_DICT'), verbose = False)
        QSD = get_qual_stats_dict(QD, output_file_path = os.path.join(self.output_directory_path, 'QUAL_STATS_DICT'), verbose = False)
        entropy_analysis(self.alignment, output_file = output_file, verbose = False, weighted = True, qual_stats_dict = QSD)
        self.assertTrue(files_are_the_same(self.expected_result, output_file))

    def test_99_CleanUp(self):
        shutil.rmtree(self.output_directory_path)
        pass
