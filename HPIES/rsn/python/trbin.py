#!/usr/bin/env python2
# trbin.py -- test read binary

"""
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150721T0001_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150723T0001_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150721T0000_UTC.dat'


  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150724T0000_UTC.dat'

  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150722T2302_UTC.dat' # turn on
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150724T0001_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150609T0000_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150622T0000_UTC.dat' # A -> B
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150703T0000_UTC.dat' # stopping
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150703T2021_UTC.dat' # stopping
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150709T1627_UTC.dat'

  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150723T0000_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150724T0001_UTC.dat' # ies_opto_on
  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150724T0000_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150714T2225_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA301_10.31.5.5_2101_20150711T0000_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20141229T0000_UTC.dat'
  fn = '/data/rsn/hpies/HPIESA101_10.33.5.5_2101_20150725T0000_UTC.dat'
"""


from __future__ import print_function

import sys
import numpy as np
from string import find
from datetime import datetime, timedelta
from collections import namedtuple
from optparse import OptionParser
import ctypes

def nextlineinit(fn):
  NL = namedtuple('NEXTLINE',[])

  ifp = open(fn,'rb')
  NL.all = np.fromfile(ifp,dtype=np.uint8)
  ifp.close()

  NL.fn = fn
  NL.ib = 0
  NL.dat = ""
  NL.sv = ""
  NL.timestamp = ""
  NL.uxt = 0
  NL.maxjlf = 0
  NL.ndat = -1
  NL.count_j2xs = 0
  NL.count_nolf = 0

  return NL

def nextline(NL):
  resync = 0

  if NL.ib == len(NL.all):
    return -1  # end of file
  if NL.ib + 16 >= len(NL.all):
    return -2  # end of file

  # check for 0xA39D7A magic of the binary time stamp
  jb = find(NL.all[NL.ib:NL.ib+16].tostring(),'\xA3\x9D\x7A')

  while True:
    if jb >= 0:
      if jb > 0:
        print('changing NL.ib by jb=',jb,end='')
        print(', skipped=',repr(NL.all[NL.ib:NL.ib+jb].tostring()))
        NL.ib += jb;
      NL.jlf = find(NL.dat,'\n')
#     print('NL.jlf=',NL.jlf)
      while NL.jlf < 0 and NL.ib < len(NL.all):
        hdr = NL.all[NL.ib:NL.ib+16]
        npak = long(hdr[4])*256 + long(hdr[5])
        pak = NL.all[NL.ib:NL.ib+npak]
        chkcalc = 0
        for i in range(npak):
          if i< 6 or i> 7:
            if pak[i] < 0:
              print('pak[i]<0')
              sys.exit(1)
            chkcalc = (chkcalc ^ pak[i]) & 0xFFFF
        chkbuf = long(pak[6]*256+pak[7])
        if chkcalc != chkbuf:
#         print('bad check sum: calc=',chkcalc,'buf=',chkbuf,'ib=',NL.ib,'npak=',npak)
          jb = find(NL.all[NL.ib:].tostring(),'\xA3\x9D\x7A')
          print('bad check sum: found next sync in',jb,'bytes.  skipped:',repr(NL.all[NL.ib:NL.ib+jb].tostring()))
          NL.ib += jb
          jb = 0
          resync = 1
          break # from while NL.jlf < 0 and NL.ib < len(NL.all):
        else:
          resync = 0
          NL.ndat = long(hdr[4])*256 + long(hdr[5]) - 16
#         print('hdr=',repr(hdr),'NL.ndat=',NL.ndat)
          NL.ib += 16
          secs = ((long(hdr[ 8]) * 256 + long(hdr[ 9])) * 256 + long(hdr[10])) * 256 + long(hdr[11])
          pico = ((long(hdr[12]) * 256 + long(hdr[13])) * 256 + long(hdr[14])) * 256 + long(hdr[15])
          NL.uxt = (secs - 25567 * 86400) + pico * 1e-9
          pyt = datetime(1970,1,1) + timedelta(0,NL.uxt)
          NL.timestamp = pyt.strftime('B %Y-%m-%d %H:%M:%S')
          NL.dat = NL.dat + NL.all[NL.ib:NL.ib+NL.ndat].tostring()
          NL.ib += NL.ndat
          # find first line feed
          NL.jlf = find(NL.dat, '\n')
#         print(NL.timestamp,'L jb=',jb, 'ib=',NL.ib, 'ndat=',NL.ndat)
      if NL.jlf < 0 and NL.ib >= len(NL.all):
        return -1  # end of file
#     print(NL.timestamp,'len(NL.dat)=',len(NL.dat),'NL.jlf=',NL.jlf)
      linein = NL.dat[:NL.jlf]
      NL.dat = NL.dat[NL.jlf+1:] # pushdown

    else:
#     print('checking for ascii time stamp')
      NL.jlf = find(NL.dat,'\n')
      while NL.jlf < 0:
        if NL.ib == len(NL.all):
          return -1  # end of file
        ie = NL.ib + 256
        if ie > len(NL.all):
          ie = len(NL.all)
        chunk = NL.all[NL.ib:ie].tostring()

        j1 = find(chunk,'<OOI-TS ')
        j2ts = find(chunk,' TS>\r\n')
        j2xs = find(chunk,' XS>\r\n')   # why XS> ????
        j3 = find(chunk, '<\\OOI-TS>\r\n')

        if j1>=0 and (j2ts>=0 or j2xs>=0) and j3>=0:
#         print('found ascii time stamp')
          if j2ts != -1 and j2xs == -1:
            j2 = j2ts
          elif j2ts != -1 and j2xs != -1:
            if j2xs < j2ts:
              j2 = j2xs
              NL.count_j2xs += 1
            else:
              j2 = j2ts
          elif j2ts == -1 and j2xs != -1:
            j2 = j2xs
            NL.count_j2xs += 1
          else:
            print('nextline(): ERROR: j2ts & j2xs both missing')
            print('  j1=',j1,'j2ts=',j2ts,'j2xs=',j2xs,'j3=',j3)
            print('  chunk=',repr(chunk))
            return -4

          if j1 == 0 and j2 > j1 and j3 > j2:
            NL.timestamp = chunk[j1:j2+4]
            NL.timestamp = 'A '+chunk[j1+8:j2-8].replace('T',' ')
            NL.dat = NL.dat + chunk[:j1] + chunk[j2+6:j3]
            NL.ib += j3+11
            NL.jlf = find(NL.dat, '\n')
          else:
            print('nextline(): ERROR: j1, j2, j3 not in sequence')
            print('  j1=',j1,'j2ts=',j2ts,'j2xs=',j2xs,'j3=',j3)
            print('  chunk=',repr(chunk))
            return -3 # ERROR

        else:
#         print('no time stamp found')
          NL.timestamp = 'N 0000-00-00 00:00:00'
          j = find(chunk, '\n')
          if j>=0:
            NL.dat = NL.dat + chunk[:j+1];
            NL.ib += j+1
          else:
            NL.count_nolf += 1
            NL.dat = NL.dat + chunk;
            NL.ib += len(chunk)
          NL.jlf = find(NL.dat, '\n')

#     if NL.jlf > NL.maxjlf:
#       NL.maxjlf = NL.jlf
      linein = NL.dat[:NL.jlf];
      NL.dat = NL.dat[NL.jlf+1:] # pushdown

    if not resync:
      break

  if NL.jlf > NL.maxjlf:
    NL.maxjlf = NL.jlf

  if len(linein) and (linein[-1] == '\r' or linein[-1] == '\n'):
    linein = linein[:-1]
  NL.oline = NL.timestamp + '  ' + linein
  return len(NL.oline)

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
      nget = nextline(NL)
      if nget == -1:
        break
      elif nget < -1:
        print('ERROR: ib=',NL.ib,'file bytes=',len(NL.all))
        sys.exit(1)
      print(repr(NL.oline))

    if options.verbose:
      print('fn=',fn)
      print('maxjlf=    ',NL.maxjlf)
      print('count_j2xs=',NL.count_j2xs)
      print('count_nolf=',NL.count_nolf)
      print('ndat=      ',NL.ndat)

if __name__ == '__main__':
  main()
