#!/usr/bin/python
# -*- coding: utf-8

# Copyright (C) 2010 - 2012, A. Murat Eren
#
# This program is free software; you can redistribute it and/or modify it under
# the terms of the GNU General Public License as published by the Free
# Software Foundation; either version 2 of the License, or (at your option)
# any later version.
#
# Please read the COPYING file.

from lib.oligotyping import Oligotyping
from lib.oligotyping import ConfigError
from utils import parsers


if __name__ == '__main__':
    parser = parsers.oligotyping()
    oligotyping = Oligotyping(parser.parse_args())

    try:
        oligotyping.run_all()
    except ConfigError, e:
        print e
        sys.exit()

