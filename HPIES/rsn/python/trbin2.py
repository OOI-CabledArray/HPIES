#!/usr/bin/env python2
# trbin.py -- test read HPIES with binary time stamps
# https://confluence.oceanobservatories.org/display/syseng/CIAD+MI+Port+Agent+Design#CIADMIPortAgentDesign-Binary

"""
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150721T0001_UTC.dat'
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150723T0001_UTC.dat'
./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150721T0000_UTC.dat'


./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150724T0000_UTC.dat'

./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150722T2302_UTC.dat # turn on
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150724T0001_UTC.dat
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150609T0000_UTC.dat
./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150622T0000_UTC.dat # A -> B
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150703T0000_UTC.dat # stopping
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150703T2021_UTC.dat # stopping
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150709T1627_UTC.dat

./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150723T0000_UTC.dat > x
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150724T0001_UTC.dat > x # ies_opto_on
./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150724T0000_UTC.dat > x
./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150714T2225_UTC.dat > x
./trbin.py /data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150711T0000_UTC.dat > x
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20141229T0000_UTC.dat > x
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150725T0000_UTC.dat > x

./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20160715T0000_UTC.dat > x0
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20160716T0245_UTC.dat > x1 # IES 177
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20160802T1830_UTC.dat > x2 # IES 177
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20160808T2020_UTC.dat > x3 # IES 308 installed
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20170117T1949_UTC.dat > x4 # restart after outage
./trbin.py /data/rsn/hpies/HPIESA101_10.33.5.5_2101_20170409T0000_UTC.dat > x5 # April 9, 2017

"""


from __future__ import print_function

import sys
import numpy as np
from string import find
from datetime import datetime, timedelta
from collections import namedtuple
from optparse import OptionParser

def uxt2str(uxt):
  pyt = datetime(1970,1,1) + timedelta(0,uxt)
  return pyt.strftime('%Y-%m-%d %H:%M:%S')

def nextlineinit(fn):
  NL = namedtuple('NEXTLINE',[])

  ifp = open(fn,'rb')
  NL.all = np.fromfile(ifp,dtype=np.uint8)
  ifp.close()

  NL.fn = fn
  NL.dat = ""
  NL.sv = ""
  NL.timestamp = ""
  NL.uxt = 0
  NL.ndat = -1
  NL.count_j2xs = 0
  NL.count_nolf = 0

# NL.ib = find(NL.all.tostring(),'\xA3\x9D\x7A')
# print('init ib=',NL.ib)
# print('init NL.all[NL.ib:NL.ib+3].tostring()=',repr(NL.all[NL.ib:NL.ib+3].tostring()))

  NL.ib = 0


  return NL

def nextpacket(NL):
  resync = 0

  if NL.ib == len(NL.all):
    return -1  # end of file
  if NL.ib + 16 >= len(NL.all):
    return -2  # end of file

  # check for 0xA39D7A magic of the binary time stamp
  jb = find(NL.all[NL.ib:NL.ib+3].tostring(),'\xA3\x9D\x7A')
  if jb != 0:
    jb = find(NL.all[NL.ib:].tostring(),'\xA3\x9D\x7A')
    NL.dat = NL.all[NL.ib:NL.ib+jb].tostring()
#   print('skipped: ',repr(NL.dat))
    NL.ib += jb
    NL.timestamp = 'N ....-..-.. ..:..:..'
    if find(NL.dat,'\n#### at ') == 0:
      j = find(NL.dat,',')
      uxt = float(NL.dat[9:j])
      NL.timestamp = 'N ' + uxt2str(uxt)
    return jb

  hdr = NL.all[NL.ib:NL.ib+16]
  npak = hdr[4]*256 + hdr[5]
  pak = NL.all[NL.ib:NL.ib+npak]

  chkcalc = 0
  for i in range(npak):
    if i< 6 or i> 7:
      if pak[i] < 0:
        print('pak[i]<0, ib=',NL.ib,'i=',i)
        sys.exit(1)
      chkcalc = (chkcalc ^ pak[i]) & 0xFFFF
  chkbuf = pak[6]*256 + pak[7]

  if chkcalc != chkbuf:
#   print('bad check sum: calc=',chkcalc,'buf=',chkbuf,'ib=',NL.ib,'npak=',npak)
    jb = find(NL.all[NL.ib:].tostring(),'\xA3\x9D\x7A')
    NL.dat = NL.all[NL.ib:NL.ib+jb].tostring()
#   print('bad chksum: nl=',loop,'synced jb=',jb)
#   print('skipped=',repr(NL.dat))
    NL.ib += jb
    NL.timestamp = 'C ....-..-.. ..:..:..'
    return jb

  secs = ((long(hdr[ 8]) * 256 + long(hdr[ 9])) * 256 + long(hdr[10])) * 256 + long(hdr[11])
  pico = ((long(hdr[12]) * 256 + long(hdr[13])) * 256 + long(hdr[14])) * 256 + long(hdr[15])
  NL.uxt = (secs - 25567 * 86400) + pico * 1e-9
  pyt = datetime(1970,1,1) + timedelta(0,NL.uxt)
  NL.timestamp = pyt.strftime('B %Y-%m-%d %H:%M:%S')

  NL.ndat = npak - 16
# print('hdr=',repr(hdr),'NL.ndat=',NL.ndat)

  NL.dat = NL.all[NL.ib+16:NL.ib+16+NL.ndat].tostring()
  NL.ib += NL.ndat + 16

  return NL.ndat

def main():

  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 0.1")

  parser.add_option("-v", "--verbose", dest="verbose",
    action="store_true", default=False,
    help="print debug info to stdout")

  (options, args) = parser.parse_args()

  if len(args) < 1:
    print('no input files specified')
    parser.print_help()
    sys.exit()

  for fn in args:
    print('fn=',fn)

    NL = nextlineinit(fn)

    while True:
      nget = nextpacket(NL)
      if nget == -1:
        break
      elif nget < -1:
        print('ERROR: ib=',NL.ib,'file bytes=',len(NL.all))
        sys.exit(1)
      print('dat:',NL.timestamp, repr(NL.dat))

    if options.verbose:
      print('fn=',fn)
      print('count_j2xs=',NL.count_j2xs)
      print('count_nolf=',NL.count_nolf)
      print('ndat=      ',NL.ndat)

if __name__ == '__main__':
  main()
