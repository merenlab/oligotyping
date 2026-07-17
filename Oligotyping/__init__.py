# -*- coding: utf-8

import os
import sys

# Make sure the Python environment hasn't changed since the installation (happens more often than you'd think
# on systems working with multiple Python installations that are managed through modules):
try:
    if sys.version_info < (3, 7, 0):
        v =  '.'.join([str(x) for x in sys.version_info[0:3]])
        sys.stderr.write("Your active Python version is '%s'. Anything less than '3.7.0' will not do it for the oligotyping pipeline :/\n" % v)
        sys.exit(-1)
except Exception:
    sys.stderr.write("(oligotyping pipeline failed to learn about your Python version, but it will pretend as if nothing happened)\n\n")

from Oligotyping.utils.utils import Run

run = Run()

def set_version():
    try:
        # `importlib.metadata` is a part of the standard library as of Python 3.8, and is the
        # modern replacement for the deprecated `pkg_resources` (which is no longer installed
        # by default with recent Python / setuptools versions).
        from importlib.metadata import version
        __version__ = version("oligotyping")
    except Exception:
        # maybe it is not installed but being run from the codebase dir?
        try:
            __version__ = open(os.path.normpath(os.path.dirname(os.path.abspath(__file__))) + '/../VERSION').read().strip()
        except Exception:
            __version__ = 'unknown'

    return __version__


def print_version():
    run.info("Oligotyping Pipeline Version", __version__, mc = 'green')

__version__ = set_version()

if '-v' in sys.argv or '--version' in sys.argv:
    print_version()
    sys.exit()
