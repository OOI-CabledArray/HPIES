#!/usr/bin/env python2
# ftdb.py -- make file time data base

import os
import sys
from optparse import OptionParser
import time
from time import strptime
from calendar import timegm
from datetime import datetime, timedelta
import glob
import scipy.io
import numpy as np

def main():
  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose", dest="verbose", 
    action="store_true", default=False,
    help="print debug info to stdout")

  parser.add_option("-r", "--rsn", dest="rsn", 
    action="store_true", default=False,
    help="for rsn defaults")

  parser.add_option("-o", "--okmc", dest="okmc", 
    action="store_true", default=False,
    help="for okmc defaults")

  (options, args) = parser.parse_args()

  if options.rsn and options.okmc:
    print('-r and -o are mutually exclusive -- exiting')
    parser.print_help()
    sys.exit(1)

  


