#!/usr/bin/python
# -*- coding: utf-8 -*-

import os
import shutil
import inspect
import unittest

from Oligotyping.lib.entropy import entropy_analysis
from Oligotyping.lib.entropy import quick_entropy

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
        self.output_directory_path = os.path.join(my_path, 'test-entropy')
        if os.path.exists(self.output_directory_path):
            shutil.rmtree(self.output_directory_path)
        os.makedirs(self.output_directory_path)
        self.alignment = os.path.join(my_path, 'files/unaligned-25K-illumina-test.fa')
        self.unique_alignment = os.path.join(my_path, 'files/unaligned-unique-25K-illumina-test.fa')
        self.expected_result = os.path.join(my_path, 'files/unaligned-25K-illumina-test-entropy.txt')

    def tearDown(self):
        pass

    def test_01_RunEntropy(self):
        output_file = os.path.join(self.output_directory_path, 'entropy.txt')
        entropy_analysis(self.alignment, output_file = output_file, verbose = False)
        self.assertTrue(files_are_the_same(self.expected_result, output_file))

    def test_02_RunEntropyOnUnique(self):
        output_file = os.path.join(self.output_directory_path, 'entropy.txt')
        entropy_analysis(self.unique_alignment, output_file = output_file, uniqued = True, verbose = False)
        self.assertTrue(files_are_the_same(self.expected_result, output_file))

    def test_03_RunQuickEntropy(self):
        n = len(quick_entropy(['ATCGATCGATCG', 'AACGATCGATGG']))
        self.assertTrue(n == 2)
        
    def test_99_CleanUp(self):
        shutil.rmtree(self.output_directory_path)
