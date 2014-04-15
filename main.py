#!/usr/bin/env python

"""
Script to execute the visdata module using a configuration file
:Author:
    Toni Sagrista Selles
:Organization:
    Astronomisches Rechen-Institut - Zentrum fur Astronomie Heidelberg - UNIVERSITAT HEIDELBERG
:Version:
    0.1
"""

from os import path
import sys

import visdata as vd


# Check args
if(len(sys.argv) != 2):
    print("usage: %s <config_file>" % path.basename(sys.argv[0]))
    exit()

vd.run(sys.argv[1])
