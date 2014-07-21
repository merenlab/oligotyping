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
        self.output_directory_path = os.path.join(my_path, 'test-decomposition')
        self.decomposer = Decomposer()
        self.decomposer.alignment = os.path.join(my_path, 'files/clone43-v6v4.fa')
        self.decomposer.min_entropy = 0.2
        self.decomposer.min_actual_abundance = 1
        self.decomposer.min_substantive_abundance = 1
        self.decomposer.skip_check_input_file = True
        self.decomposer.number_of_discriminants = 3
        self.decomposer.progress.verbose = False
        self.decomposer.run.verbose = False
        self.decomposer.skip_removing_outliers = True
        self.decomposer.skip_agglomerating_nodes = True
        self.decomposer.output_directory = self.output_directory_path 


    def tearDown(self):
        pass

    def test_01_PrefixGenerator(self):
        prefix = self.decomposer.get_prefix()
        self.assertTrue(prefix == 'm0.20-A1-M1-d3')

    def test_02_Decompose(self):
        self.decomposer.decompose()
        self.assertTrue(self.decomposer.topology.nodes['root'].discriminants == [292, 296, 293])
        self.assertTrue(self.decomposer.topology.nodes['root'].entropy_tpls[0:5] == [(292, 2.2031591775819033), (296, 2.1539606786564707), (293, 2.1279277605456857), (298, 2.1170394670626229), (300, 2.107517097722345)])
        self.assertTrue(len(self.decomposer.topology.nodes) == 54)
        self.assertTrue(len(self.decomposer.topology.final_nodes) == 43)

    def test_03_Environment(self):
        self.assertTrue(files_are_the_same(os.path.join(my_path, 'files/clone43-v6v4-environment.txt'),
                                           os.path.join(self.output_directory_path, 'ENVIRONMENT.txt')))

    def test_04_MatrixPercent(self):
        self.assertTrue(files_are_the_same(os.path.join(my_path, 'files/clone43-v6v4-matrix-percent.txt'),
                                           os.path.join(self.output_directory_path, 'MATRIX-PERCENT.txt')))

    def test_99_CleanUp(self):
        shutil.rmtree(self.output_directory_path)
