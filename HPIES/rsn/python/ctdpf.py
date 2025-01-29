#!/usr/bin/env python2
# ctdpf.py -- decode, plot and process hex data

# ./ctdpf.py /data/rsn/ctd/CTDPFB301_10.31.5.12_2101_20140831T0000_UTC.dat

from __future__ import print_function

import os
import sys
from optparse import OptionParser
import time
import math
import numpy as np
import matplotlib
matplotlib.use('Agg') # needed when used with cron or no DISPLAY or slow remote connection
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
# import matplotlib.ticker as mticker
import matplotlib.dates as mdates
import collections
import gzip
from datetime import datetime, timedelta
import glob
from time import strptime
from calendar import timegm
import copy

def mydispfig(pdfdir,pltnam):
  mymkdir(pdfdir)
  pdffile = pdfdir + '/' + pltnam + '.pdf'
  print(pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pdfdir)

def mymkdir(mydir):
  try: 
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      sys.exit(1)

def fix_xdates(axlist,nss):
  fmt = mdates.DateFormatter('%b %d %H%MZ')
  for ax in axlist:
    xax = ax.get_xaxis()
    xax.set_major_formatter( fmt )

    dx = 1.0/24.0
    xlim = ax.get_xlim()
    if xlim[-1] - xlim[0] > 4.0/24.0:
      ax.set_xlim([np.floor(xlim[0]+dx)-dx,np.ceil(xlim[-1]-dx)+dx])

    if nss>1:
      xticks = ax.get_xticks()
      newxticks = []
      for i in range(0,len(xticks),nss):
        newxticks.append(xticks[i])
      ax.set_xticks(newxticks)


def mymain():

  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose", dest="verbose", 
    action="store_true", default=False,
    help="print debug info to stdout")

  (options, args) = parser.parse_args()

  ifiles=[]
  for arg in args:
    for ifile in glob.glob(arg):
      ifiles.append(ifile)
  for ifile in sorted(ifiles):
    print('ifile=',ifile)

    leafname = os.path.basename(ifile)

    ifd = open(ifile,'r')

    linepie=''
    lineno=0
    uxt  = []
    var1 = []
    var2 = []
    var3 = []
    var4 = []
    var5 = []
    var6 = []
    var7 = []

    tgm = timegm(strptime('20140901T010203','%Y%m%dT%H%M%S'))
    ymdhms = '19700101T000000'

    for lineraw in ifd:
      while len(lineraw)>1 and (lineraw[-1] == '\r' or lineraw[-1] == '\n'):
        lineraw = lineraw[0:-1]
      if len(lineraw)<1:
        continue

      # replace OOI-TS time stamps with hprtst.py-style time stamps
      i1a = lineraw.find('<OOI-TS ')
      if i1a>0:
        print('error: <OOI-TS not at front of line')
        print('i1a=',i1a,'lineraw=',lineraw)
        sys.exit()
      if i1a==0:
        if len(lineraw) != 39:
          print('len(lineraw)=',len(lineraw),'not right')
          print('lineraw=',lineraw)
          sys.exit()
        i1b = lineraw.find('Z TN>')
        i1c = lineraw.find('Z TS>')
        i1d = lineraw.find('Z XS>')
        if i1b != 34 and i1c != 34 and i1d != 34:
          print('i1b=',i1b,'or i1c=',i1c,'not right')
          print('lineraw=',lineraw)
          sys.exit()
        ymdhms = lineraw[8:12] + lineraw[13:15] + lineraw[16:21] + lineraw[22:24] + lineraw[25:27]
  #       print('ymdhms:',ymdhms)
        continue
      
      i2 = lineraw.find('<\\OOI-TS>')
      if i2 > 0:
        if i2 != len(lineraw) - 9:
          print('error: <\\OOI-TS> not at end of line')
          print('i2=',i2,'len(lineraw)=',len(lineraw),'lineraw=',lineraw)
          sys.exit()
        linepie += lineraw[0:i2]
        continue
      if i2 == 0:
        continue

      if i1a < 0 and i2 < 0:
        linein = ymdhms + ' ' + linepie + lineraw
        linepie = ''
      else:
        print('unexpected: i1a=',i1a,'i2=',i2)
        sys.exit()

      lineno += 1

      if len(linein) != 52:
        print('len(linein)=',len(linein))
        continue

      try:
        secs = timegm(strptime(linein[0:15],'%Y%m%dT%H%M%S'))
      except:
        print('cannot decode time from linein=',linein)
        continue

      if secs < 10:
        continue

      uxt.append (secs)
      var1.append(int(linein[16:22],16))
      var2.append(int(linein[22:28],16))
      var3.append(int(linein[28:34],16))
      var4.append(int(linein[34:38],16))
      var5.append(int(linein[38:42],16))
      var6.append(int(linein[42:46],16))
      var7.append(int(linein[46:52],16))

    ifd.close()

    mlt = np.array(uxt,dtype='double') / 86400 + 719529 - 366
    var1 = np.array(var1,dtype='double')
    var2 = np.array(var2,dtype='double')
    var3 = np.array(var3,dtype='double')
    var4 = np.array(var4,dtype='double')
    var5 = np.array(var5,dtype='double')
    var6 = np.array(var6,dtype='double')
    var7 = np.array(var7,dtype='double')

    print('linein=',linein)
    print('mlt=',mlt[-1])
    print('var1=',var1[-1])
    print('var2=',var2[-1])

    # Temperature coefficients for CTDPFB301
    TA0 =  1.305824e-03
    TA1 =  2.632108e-04
    TA2 = -2.236076e-07
    TA3 =  1.438215e-07
    TOFFSET = 0.000000e+00

  # v = var1 / 13107.0
  # T_L1 = var1 / 10000.0 - 10.0   # from DPS version 1-00, Jan 2012
  # from cal sheet:
  # T90 = 1/{a0 + a1[ln(n)] + a2[ln^2(n)] + a3[ln^3(n)]} - 273.15
  # lnn = np.log(var1 / 100.0)
  # T_L1 = 1.0 / (TA0 + lnn * (TA1 + lnn * (TA2 + lnn * TA3))) - 273.15 + TOFFSET

  # from version 1-04 of DPS
  # for SBE 16plus V2, Output Format 0
    MV = (var1 - 524288) / 1.6e7
    R = (MV * 2.900e9 + 1.024e8) / (2.048e4 - MV * 2.0e5)
    lnR = np.log(R)
    T_L1 = 1.0/(TA0 + lnR*(TA1 + lnR*(TA2 + lnR*TA3))) - 273.15

    # Pressure coefficients for CTDPFB301
    PA0     =  5.466934e-01
    PA1     =  1.566520e-02
    PA2     = -6.428222e-10
    PTCA0   =  5.244772e+05
    PTCA1   =  2.172837e+00
    PTCA2   = -2.416201e-02
    PTCB0   =  2.504512e+01
    PTCB1   =  4.250000e-04
    PTCB2   =  0.000000e+00
    PTEMPA0 = -6.560674e+01
    PTEMPA1 =  5.257837e+01
    PTEMPA2 = -2.586355e-01
    POFFSET =  0.000000e+00
    PRANGE  =  5.076000e+03

    tv = var4 / 13107.0
    t = PTEMPA0 + tv * (PTEMPA1 + tv * PTEMPA2)
    x = var3 - PTCA0 - t * (PTCA1 + t * PTCA2)
    n = x * PTCB0 / (PTCB0 + t * (PTCB1 + t * PTCB2))
    P_psia = PA0 + n * (PA1 + n * PA2)
    P_L1 = P_psia * 0.689475729 - 10.1325
    TofP = t
    
    # Conductivity coefficients for CTDPFB301
    G     = -9.741097e-01
    H     =  1.359033e-01
    I     = -1.468517e-04
    J     =  2.867985e-05
    CPCOR = -9.570000e-08
    CTCOR =  3.250000e-06
    CSLOPE = 1.000000e+00

    f = var2 / 256.0 / 1000.0 # kHz
    C_L1 = G + f*f*(H + f*(I + f*J)) / (1.0 + CTCOR*T_L1 + CPCOR*P_L1)


    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()

    fig.suptitle(leafname)

    axlist = []

    ax = fig.add_subplot(4,2,1)
    axlist.append(ax)
    ax.plot_date(mlt,T_L1,'b-')
    plt.ylabel('T, deg C')

    ax = fig.add_subplot(4,2,2)
    axlist.append(ax)
    ax.plot_date(mlt,C_L1,'b-')
    plt.ylabel('C, S/m')

    ax = fig.add_subplot(4,2,3)
    axlist.append(ax)
    ax.plot_date(mlt,P_L1,'b-')
    plt.ylabel('P, dbar')

    ax = fig.add_subplot(4,2,4)
    axlist.append(ax)
    ax.plot_date(mlt, TofP ,'b-')
    plt.ylabel('TofP, deg C')

    ax = fig.add_subplot(4,2,5)
    axlist.append(ax)
    ax.plot_date(mlt,var5 / 13107.0 ,'b-')
    plt.ylabel('Volts0')

    ax = fig.add_subplot(4,2,6)
    axlist.append(ax)
    ax.plot_date(mlt,var6 / 13107.0 ,'b-')
    plt.ylabel('Volts1')

    ax = fig.add_subplot(4,2,7)
    axlist.append(ax)
    ax.plot_date(mlt,var7 / 10000.0 - 10.0,'b-')
    plt.ylabel('Optode, um/l')

    fig.autofmt_xdate()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, \
      wspace=0.3, hspace=None)
    fix_xdates(axlist,2)
    mydispfig('./pdf-ctdpf',leafname)

'''
getcc
<CalibrationCoefficients DeviceType = 'SBE19plus' SerialNumber = '01950031'>
   <Calibration format = 'TEMP1' id = 'Main Temperature'>
      <SerialNum>01950031</SerialNum>
      <CalDate>07-May-14</CalDate>
      <TA0>1.305824e-03</TA0>
      <TA1>2.632108e-04</TA1>
      <TA2>-2.236076e-07</TA2>
      <TA3>1.438215e-07</TA3>
      <TOFFSET>0.000000e+00</TOFFSET>
   </Calibration>
   <Calibration format = 'WBCOND0' id = 'Main Conductivity'>
      <SerialNum>01950031</SerialNum>
      <CalDate>07-May-14</CalDate>
      <G>-9.741097e-01</G>
      <H>1.359033e-01</H>
      <I>-1.468517e-04</I>
      <J>2.867985e-05</J>
      <CPCOR>-9.570000e-08</CPCOR>
      <CTCOR>3.250000e-06</CTCOR>
      <CSLOPE>1.000000e+00</CSLOPE>
   </Calibration>
   <Calibration format = 'VOLT0' id = 'Volt 0'>
      <OFFSET>-4.828842e-02</OFFSET>
      <SLOPE>1.251476e+00</SLOPE>
   </Calibration>
   <Calibration format = 'VOLT0' id = 'Volt 1'>
      <OFFSET>-4.879263e-02</OFFSET>
      <SLOPE>1.251588e+00</SLOPE>
   </Calibration>
   <Calibration format = 'STRAIN0' id = 'Main Pressure'>
      <SerialNum>3898298</SerialNum>
      <CalDate>05-May-14</CalDate>
      <PA0>5.466934e-01</PA0>
      <PA1>1.566520e-02</PA1>
      <PA2>-6.428222e-10</PA2>
      <PTCA0>5.244772e+05</PTCA0>
      <PTCA1>2.172837e+00</PTCA1>
      <PTCA2>-2.416201e-02</PTCA2>
      <PTCB0>2.504512e+01</PTCB0>
      <PTCB1>4.250000e-04</PTCB1>
      <PTCB2>0.000000e+00</PTCB2>
      <PTEMPA0>-6.560674e+01</PTEMPA0>
      <PTEMPA1>5.257837e+01</PTEMPA1>
      <PTEMPA2>-2.586355e-01</PTEMPA2>
      <POFFSET>0.000000e+00</POFFSET>
      <PRANGE>5.076000e+03</PRANGE>
'''

if __name__ == '__main__':
  mymain()
