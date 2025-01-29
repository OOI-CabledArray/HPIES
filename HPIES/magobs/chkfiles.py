#! /usr/bin/env python
# chkfiles.py -- see which files are missing

from __future__ import print_function
from datetime import datetime, timedelta
import scipy.io
from optparse import OptionParser
import os
import gzip

def main():

  parser = OptionParser(usage='%prog [Options] stations')

  parser.add_option('--datebeg', metavar='YMD',
    dest='datebeg', default=None,
    help='yyyy-mm-dd of first day')

  parser.add_option('--dateend', metavar='YMD',
    dest='dateend', default=None,
    help='yyyy-mm-dd of last day')

  parser.add_option('-v', '--verbose',
    action='count', dest='verbose', default=0,
    help='print status messages to stdout')

  parser.add_option('--idir', metavar='D',
    dest='idir', default='/data/magobs/usgs',
    help='input directory')

  parser.add_option('-q', '--quick', dest='quick', default=False,
    action='store_true',help='quick check of .gz file size')

  parser.add_option('-e', '--endsame', dest='endsame', default=False,
    action='store_true',help='make dateend same as datebeg')

  (options, args) = parser.parse_args()

  if len(args) == 0:
    parser.print_help()
    exit(1)

  if options.datebeg == None:
    datebeg = datetime.utcnow() - timedelta(days=1)
  else:
    datebeg = datetime.strptime(options.datebeg, '%Y-%m-%d')
  print('datebeg=',datebeg.strftime('%Y-%m-%d'))

  if options.dateend == None:
    if options.endsame:
      dateend = datebeg
    else:
      dateend = datetime.utcnow() - timedelta(days=1)
  else:
    dateend = datetime.strptime(options.dateend, '%Y-%m-%d')
  print('dateend=',dateend.strftime('%Y-%m-%d'))

  # filedate = datetime.strptime('2012-06-01 17:05', '%Y-%m-%d %H:%M')
  nfound = 0
  nfoundgz = 0
  ntried = 0
  missing = 0

  for stn in args:
    print('testing',stn)

    filedate = datebeg
    while filedate <= dateend:
      filename = str.format('{0}/{1}/{2}{3}vsec.sec', \
                 options.idir,stn.upper(),stn.lower(), \
                 filedate.strftime('%Y%m%d'))
      
#     ifd = None
#     try:
#       ifd = open(filename,'r')
#     except IOError:
#       try:
#         ifd = open(filename+'.gz','r')
#       except IOError:
#         print(filename + ' is missing')

#     if ifd != None:
#       ifd.close()
#       nfound += 1

      found = False
      if os.path.exists(filename):
        nfound += 1
        found = True
        if os.stat(filename).st_size != 6135749:
          print(filename,'   is not correct, len=',os.stat(filename).st_size)

      filename += '.gz'
      if os.path.exists(filename):
        nfoundgz += 1
        found = True
        if options.quick:
          if os.stat(filename).st_size < 2000:
            print(filename,'is too short, len=',os.stat(filename).st_size)
        else:
          with gzip.open(filename) as fp:
            lins = fp.readlines()
          nbytes = 0
          for lin in lins:
            nbytes += len(lin)
          if nbytes != 6135749 and nbytes != 6135891 and nbytes != 6135607:
            print(filename,'nbytes=',nbytes,'should be 6135749 or 6135891 or 6135607')

      if not found:
        print(filename,'is missing')
        missing += 1

      ntried += 1
      filedate += timedelta(days=1)

  print(str.format('#    files found: {0}',nfound))
  print(str.format('# gz files found: {0}',nfoundgz))
  print(str.format('#    files tried: {0}',ntried))
  print(str.format('#  files missing: {0}',missing))

if __name__ == '__main__':
  main()
