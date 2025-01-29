#!/usr/bin/env python2
# getmagobs.py

# ./getmagobs.py -v --datebeg 2021-05-28 --dateend 2021-11-11

# 2021-11-12: ws/edge replaced by ws/data/
# wget http://geomag.usgs.gov/ws/data/\
#\?id=BOU\
#\&starttime=2016-11-01T00:00:00Z\
#\&endtime=2016-11-01T23:59:59\
#\&elements=X,Y,Z,F\
#\&sampling_period=1\
# -O boutst

# new way
# wget http://geomag.usgs.gov/ws/edge/\
#\?id=BOU\
#\&starttime=2016-11-01T00:00:00Z\
#\&endtime=2016-11-01T23:59:59\
#\&elements=X,Y,Z,F\
#\&sampling_period=1\
# -O boutst

# old way
# from ohm:~dunlap/magobs/usgs/wget1s.sh
# for S in GUA HON FRD FRN NEW SHU TUC SIT BSL SJG
# do
#   s=`echo $S | tr 'A-Z' 'a-z'`
#   d=`date +%Y%m%d`
#   H="http://geomag.usgs.gov/data/magnetometer/"$S"/OneSecond/"$s$d"vsec.sec"
#   wget -nv -P ~dunlap/magobs/usgs/$S $H
# done

from __future__ import print_function
import sys, os
from optparse import OptionParser
from datetime import datetime, timedelta
import urllib
from time import sleep

def mymkdir(mydir):
  try:
    os.makedirs(mydir)
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:', mydir)
      sys.exit(1)


def main():

  print('getmagobs.py starting at',datetime.utcnow(),'UTC')

  parser = OptionParser(usage='%prog [Options] stations')

  parser.add_option('-v', '--verbose', 
    action='count', dest='verbose', default=0,
    help='print status messages to stdout')

  parser.add_option('--odir', metavar='D',
    dest='odir', default='/data/magobs/usgs',
    help='output directory')

  parser.add_option('--datebeg', metavar='YMD',
    dest='datebeg', default=None,
    help='yyyy-mm-dd of first day')

  parser.add_option('--dateend', metavar='YMD',
    dest='dateend', default=None,
    help='yyyy-mm-dd of last day')

  parser.add_option('--ndays', metavar='N',
    dest='ndays', default=1,
    help='number of days')

  parser.add_option('--delay', metavar='N',
    type='float', dest='delay', default=60.0,
    help='seconds to delay between requests')

  (options, args) = parser.parse_args()

  if len(args) == 0:
    parser.print_help()
    sys.exit(1)

  if options.datebeg is not None and options.dateend is None:
    options.dateend = options.datebeg

  if options.datebeg == None:
    datebeg = datetime.utcnow() - timedelta(days=1)
  else:
    try:
      datebeg = datetime.strptime(options.datebeg, '%Y-%m-%d')
    except:
      datebeg = datetime.strptime(options.datebeg, '%Y%m%d')
  print('datebeg=',datebeg.strftime('%Y-%m-%d'))

  if options.dateend == None:
    dateend = datetime.utcnow() - timedelta(days=1)
  else:
    try:
      dateend = datetime.strptime(options.dateend, '%Y-%m-%d')
    except:
      dateend = datetime.strptime(options.dateend, '%Y%m%d')
  print('dateend=',dateend.strftime('%Y-%m-%d'))

  dateget = datebeg
  while dateget <= dateend:
    print('dateget=',dateget.strftime('%Y-%m-%d'))

    for stn in args:
      print('  stn=',stn)
      url = 'http://geomag.usgs.gov/ws/data/' + \
            '?id=' + stn.upper() + \
            '&starttime=' + dateget.strftime('%Y-%m-%d') + 'T00:00:00Z' + \
            '&endtime='  + dateget.strftime('%Y-%m-%d') + 'T23:59:59Z' + \
            '&elements=H,D,Z,F' + \
            '&sampling_period=1'
      if options.verbose:
        print('  url=',url)

      stndir = '{0:s}/{1:s}'.format(options.odir,stn.upper())
      opath = '{0:s}/{1:s}{2:s}vsec.sec'.\
              format(stndir, stn.lower(), dateget.strftime('%Y%m%d'))
      if options.verbose:
        print('  opath=',opath)

      mymkdir(stndir)
      for itry in range(3):
        try:
          urllib.urlretrieve(url,filename=opath)
        except:
          print('  urlretrieve failed:')
          print('    itry=',itry)
          print('    url=',url)
          print('    opath=',opath)
          print('    will try again in 600 s')
          sleep(600)
          print('  stn=',stn)

      cmd = 'gzip --force {0}'.format(opath)
      print('  cmd=',cmd)
      os.system(cmd)

      if len(args) > 1 or datebeg != dateend:
        sleep(options.delay) # wait between requests so as to not hog the server

    # end for stn in args:

    dateget += timedelta(days=1)

  # end while dateget <= dateend:

  for stn in args:
    stndir = '{0:s}/{1:s}'.format(options.odir,stn.upper())
    cmd = 'rsync -aL -e "ssh -i /home/dunlap/.ssh/o2k_id_rsa -o BatchMode=yes" {0} kirin:/var/www/html/magobs/usgs'.format(stndir)
    print('cmd=',cmd)
    os.system(cmd)

  print('getmagobs.py finished at',datetime.utcnow(),'UTC')

if __name__ == '__main__':
  main()
