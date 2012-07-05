#!/usr/bin/python
# -*- coding: utf-8 -*-

import unittest

#unittest declerations
import _oligotyping

def __suite():
    suite = unittest.TestSuite()
    suite.addTest(unittest.makeSuite(_oligotyping.Tests))

    return suite

if __name__ == "__main__":
    unittest.TextTestRunner(verbosity=3).run(__suite())
