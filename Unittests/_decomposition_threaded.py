#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import unittest

from Oligotyping.lib.decomposer import Decomposer

my_path = os.path.dirname(os.path.realpath(__file__))

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
        self.output_directory_path = os.path.join(my_path, 'test-decomposition-threaded')
        self.decomposer = Decomposer()
        self.decomposer.alignment = os.path.join(my_path, 'files/reads-noisy.fa')
        self.decomposer.min_entropy = 0.3
        self.decomposer.min_actual_abundance = 0
        self.decomposer.min_substantive_abundance = 2
        self.decomposer.skip_check_input_file = True
        self.decomposer.number_of_discriminants = 1
        self.decomposer.progress.verbose = False
        self.decomposer.run.verbose = False
        self.decomposer.skip_removing_outliers = False
        self.decomposer.relocate_outliers = True
        self.decomposer.skip_agglomerating_nodes = True
        self.decomposer.threading = True
        self.decomposer.output_directory = self.output_directory_path 


    def tearDown(self):
        pass

    def test_02_Decompose(self):
        self.decomposer.decompose()

    def test_02_Environment(self):
        self.assertTrue(files_are_the_same(os.path.join(my_path, 'files/reads-noisy-environment.txt'),
                                           os.path.join(self.output_directory_path, 'ENVIRONMENT.txt')))

    def test_03_MatrixPercent(self):
        self.assertTrue(files_are_the_same(os.path.join(my_path, 'files/reads-noisy-matrix-percent.txt'),
                                           os.path.join(self.output_directory_path, 'MATRIX-PERCENT.txt')))

    def test_99_CleanUp(self):
        shutil.rmtree(self.output_directory_path)
