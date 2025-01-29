#!/usr/bin/env python2
# vel3d.py -- decode, plot and process VEL3DB301...

# ./vel3d.py /data/rsn/vel3d/VEL3DB301_10.31.2.12_2101_20140903T0000_UTC.dat

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

    try:
      ifd = open(ifile,'r')
    except:
      print('ifile does not exist')
      continue

    linepie=''
    lineno=0
    uxt  = []
    VA = []
    VB = []
    VC = []
    VD = []
    E = []
    N = []
    W = []
    T = []
    MX = []
    MY = []
    Pitch = []
    Roll = []

    tgm = timegm(strptime('20140901T010203','%Y%m%dT%H%M%S'))

    for linein in ifd:
      while len(linein)>1 and (linein[-1] == '\r' or linein[-1] == '\n'):
        linein = linein[0:-1]
      if len(linein)<1:
        continue

      lineno += 1

      x = linein.split()
      if len(x) != 18:
        continue
      if x[2] != '2014':
        continue

      ymdhms = x[2]+x[0]+x[1]+'T'+x[3]+x[4]+x[5][0:2]
      try:
        secs = timegm(strptime(ymdhms,'%Y%m%dT%H%M%S'))
      except:
        print('cannot decode time from linein=',linein)
        print('ymdhms=',ymdhms)
        sys.exit(1)

      uxt.append (secs)
      VA.append(float(x[6]))
      VB.append(float(x[7]))
      VC.append(float(x[8]))
      VD.append(float(x[9]))
      E.append(float(x[10]))
      N.append(float(x[11]))
      W.append(float(x[12]))
      T.append(float(x[13]))
      MX.append(float(x[14]))
      MY.append(float(x[15]))
      Pitch.append(float(x[16]))
      Roll.append(float(x[17]))

    ifd.close()

    mlt = np.array(uxt,dtype='double') / 86400 + 719529 - 366
    VA = np.array(VA,dtype='double')
    VB = np.array(VB,dtype='double')
    VC = np.array(VC,dtype='double')
    VD = np.array(VD,dtype='double')
    E = np.array(E,dtype='double')
    N = np.array(N,dtype='double')
    W = np.array(W,dtype='double')
    T = np.array(T,dtype='double')
    MX = np.array(MX,dtype='double')
    MY = np.array(MY,dtype='double')
    Pitch = np.array(Pitch,dtype='double')
    Roll = np.array(Roll,dtype='double')

    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()

    fig.suptitle(leafname)

    axlist = []

    ax = fig.add_subplot(4,3,1)
    axlist.append(ax)
    ax.plot_date(mlt,VA,'b-')
    plt.ylabel('VA')

    ax = fig.add_subplot(4,3,2)
    axlist.append(ax)
    ax.plot_date(mlt,VB,'b-')
    plt.ylabel('VB')

    ax = fig.add_subplot(4,3,3)
    axlist.append(ax)
    ax.plot_date(mlt,VC,'b-')
    plt.ylabel('VC')

    ax = fig.add_subplot(4,3,4)
    axlist.append(ax)
    ax.plot_date(mlt, VD ,'b-')
    plt.ylabel('VD')

    ax = fig.add_subplot(4,3,5)
    axlist.append(ax)
    ax.plot_date(mlt,E ,'b-')
    plt.ylabel('E')

    ax = fig.add_subplot(4,3,6)
    axlist.append(ax)
    ax.plot_date(mlt,N ,'b-')
    plt.ylabel('N')

    ax = fig.add_subplot(4,3,7)
    axlist.append(ax)
    ax.plot_date(mlt,W,'b-')
    plt.ylabel('W')

    ax = fig.add_subplot(4,3,8)
    axlist.append(ax)
    ax.plot_date(mlt,T,'b-')
    plt.ylabel('T')

    ax = fig.add_subplot(4,3,9)
    axlist.append(ax)
    ax.plot_date(mlt,W,'b-')
    plt.ylabel('W')

    ax = fig.add_subplot(4,3,10)
    axlist.append(ax)
    ax.plot_date(mlt,MX,'b-')
    plt.ylabel('MX')

    ax = fig.add_subplot(4,3,11)
    axlist.append(ax)
    ax.plot_date(mlt,MY,'b-')
    plt.ylabel('MY')

    ax = fig.add_subplot(4,3,12)
    axlist.append(ax)
    ax.hold(True)
    ax.plot_date(mlt,Pitch,'b-')
    ax.plot_date(mlt,Roll,'r-')
    ax.hold(False)
    plt.ylabel('P&R')

    fig.autofmt_xdate()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=None, \
      wspace=0.3, hspace=None)
    fix_xdates(axlist,2)
    mydispfig('./pdf-vel3d',leafname)

if __name__ == '__main__':
  mymain()
