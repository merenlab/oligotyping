#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import inspect
import unittest

from Oligotyping.lib.oligotyping import Oligotyping

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
        self.output_directory_path = os.path.join(my_path, 'test-oligotyping')
        self.oligotyping = Oligotyping()
        self.oligotyping.alignment = os.path.join(my_path, 'files/unaligned-25K-illumina-test.fa')
        self.oligotyping.entropy = os.path.join(my_path, 'files/unaligned-25K-illumina-test-entropy.txt')
        self.oligotyping.number_of_auto_components = 20
        self.oligotyping.min_percent_abundance = 1
        self.oligotyping.min_actual_abundance = 100
        self.oligotyping.min_number_of_datasets = 1
        self.oligotyping.project = 'Unaligned 25K Illumina Test'
        self.oligotyping.quick = True
        self.oligotyping.output_directory = self.output_directory_path 
        self.oligotyping.no_display = True
        self.oligotyping.progress.verbose = False
        self.oligotyping.run.verbose = False

    def tearDown(self):
        pass

    def test_01_PrefixGenerator(self):
        prefix = self.oligotyping.get_prefix()
        self.assertTrue(prefix == 'c20-s1-a1.0-A100-M0')

    def test_02_Oligotyping(self):
        self.oligotyping.run_all()

    def test_03_EnvironmentFile(self):
        self.assertTrue(files_are_the_same(os.path.join(my_path, 'files/unaligned-25K-illumina-ENVIRONMENT.txt'),
                                           os.path.join(self.output_directory_path, 'ENVIRONMENT.txt')))
       
    def test_99_CleanUp(self):
        shutil.rmtree(self.output_directory_path)
