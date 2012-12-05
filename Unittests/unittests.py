#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

# unittest declerations
import _entropy
import _weightedEntropy
import _oligotyping
import _decomposition
import _decomposition_threaded

def __suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(_entropy.Tests))
    suite.addTest(unittest.makeSuite(_weightedEntropy.Tests))
    suite.addTest(unittest.makeSuite(_oligotyping.Tests))
    suite.addTest(unittest.makeSuite(_decomposition.Tests))
    suite.addTest(unittest.makeSuite(_decomposition_threaded.Tests))

    return suite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=3).run(__suite())
