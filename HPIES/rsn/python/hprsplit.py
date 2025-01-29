#!/usr/bin/env python2
# hprsplit.py -- split old python files into daily files
# John Dunlap, APL/UW

from __future__ import print_function

import sys
import os
import gzip
import time
from calendar import timegm
from datetime import datetime, timedelta
from optparse import OptionParser


def mymkdir(mydir):
  try: 
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      cleanup()

def gz_open(ifile):
  if ifile[-3:] == '.gz':
    try:
      ifd = gzip.open(ifile,'r')
    except:
      print('1 cannot open ifile=',ifile)
      sys.exit(1)
  else:
    try:
      ifd = gzip.open(ifile+'.gz','r')
    except:
      try:
        ifd = open(ifile,'rt')
      except:
        print('2 cannot open ifile=',ifile)
        sys.exit(1)
  return ifd


def main():

  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 1.0")

  parser.add_option("-b", "--basedir", dest="basedir", 
    metavar='D', default='./xxx',
    help="base directory [default: %default]")

  parser.add_option("-n", "--name", dest="hpiesname", 
    metavar='NAM', default=None,
    help="HPIES name, e.g., hpr001")

  (options, args) = parser.parse_args()

  if options.hpiesname == None:
    print('need hpiesname')
    parser.print_help()
    sys.exit(1)

  if len(args) == 0:
    print('no input files specified')
    parser.print_help()
    sys.exit(1)

  ymd = '19700101'
  ofd = None

  for ifile in args:
    print('ifile=',ifile)

    ifd = gz_open(ifile)

    for line in ifd:
      if len(line) >= 16 and line[8] == 'T' and line[15] == ' ':
        if line[0:8] != ymd:
          if ofd:
            ofd.close()
          ymd = line[0:8]
          odir = options.basedir + '/' + options.hpiesname
          mymkdir(odir)
          ofile =  odir + '/' + options.hpiesname + '_' + ymd + 'UTC.txt'
          print('ofile=',ofile)
          ofd = open(ofile,'a')
      print(line,file=ofd,end='')

  ofd.close()
        


if __name__ == '__main__':
  main()


