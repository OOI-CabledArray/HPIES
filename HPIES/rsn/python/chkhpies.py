#!/usr/bin/env python3
# chkhpies.py -- check for missing hpies files

from datetime import datetime, timedelta
# import scipy.io
# import os

from optparse import OptionParser
from glob import glob


def main():

  global printed
  global pattern
  global idx1
  global idx2
  global mis1
  global mis2

  parser = OptionParser(usage='%prog [Options] stations')

  parser.add_option('-v', '--verbose',
    action='count', dest='verbose', default=0,
    help='print status messages to stdout')

  parser.add_option('--datebeg', metavar='YMD',
    dest='datebeg', default=None,
    help='yyyy-mm-dd of first day')

  parser.add_option('--dateend', metavar='YMD',
    dest='dateend', default=None,
    help='yyyy-mm-dd of last day')

  parser.add_option('--idir', metavar='D',
    dest='idir', default='/data/rsn/hpies',
    help='input directory')

  (options, args) = parser.parse_args()

  if options.datebeg is None:
    datebeg = datetime.utcnow() - timedelta(days=1)
  else:
    datebeg = datetime.strptime(options.datebeg, '%Y-%m-%d')
  print('datebeg=',datebeg.strftime('%Y-%m-%d'))

  if options.dateend is None:
    dateend = datetime.utcnow() - timedelta(days=1)
  else:
    dateend = datetime.strptime(options.dateend, '%Y-%m-%d')
  print('dateend=',dateend.strftime('%Y-%m-%d'))

  prefs = ['HPIESA101_10.33.5.5_2101_','HPIESA301_10.31.5.5_2101_']
  suf = '_UTC.dat'

  print('idir=',options.idir)
  print('prefs=',prefs)
  print('suf=',suf)
  print('------------------------------------------------------')

  for pref in prefs:
    mis1 = None
    mis2 = None
    idx = 0
    printed = False
    filedate = datebeg
    while filedate <= dateend:
      idx += 1
      pattern = options.idir + '/' + pref + filedate.strftime('%Y%m%dT????' + suf)
      files = glob(pattern)
      if options.verbose > 1:
        print(pattern)
        print(files)
      missdate = filedate.strftime('%Y-%m-%d')
      if len(files) == 0:
        if options.verbose:
          print('                                   ',end=' ')
          print('idx0=',idx,'misdate= ',missdate)
        idx2 = idx
        mis2 = missdate
        if mis1 is None:
          idx1 = idx
          mis1 = missdate
      else:
        if mis1 is not None:
          doprint()
          mis1 = None

      filedate += timedelta(days=1)
    if mis1 is not None:
      doprint()
  print('------------------------------------------------------')

def doprint():
  global printed
  if not printed:
    print(pattern)
    printed = True
  print('nmis={0:3d}'.format(idx2-idx1+1),end=' ')
  print('{0:3d} {1:3d}'.format(idx1,idx2),mis1,mis2)

if __name__ == '__main__':
  main()
