#!/usr/bin/env python2
# gpgga.py -- decode & plot R/V Thompson ftp indian cd /ship-data/scs/data/NAV

# ./gpgga.py /data/rsn/nav/CNAV3050-GGA-RAW_*.Raw -p rsn -m rsn

from __future__ import print_function

import os
import sys
from optparse import OptionParser
import time
import math
import numpy as np
import matplotlib
matplotlib.use('Agg') # needed when used with cron 
                      # or no DISPLAY 
                      # or slow remote connection
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
import scipy.io

def writepdf(pdfdir, namewoext):
  mymkdir(pdfdir)
  pdffile = os.path.join(pdfdir,namewoext+'.pdf')
  print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pdfdir)

def writemat(matdir,namewoext,matdict):
  mymkdir(matdir)
  matfile = os.path.join(matdir,namewoext+'.mat')
  print('matfile=',matfile)
  scipy.io.savemat(matfile, matdict, format='4')

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

  parser.add_option("-p", "--pltdir", dest="pltdir", 
    default=None, help="plot directory")

  parser.add_option("-m", "--matdir", dest="matdir", 
    default=None, help="matlab output directory")

  (options, args) = parser.parse_args()

  if options.matdir == 'rsn':
    options.matdir = '/data/rsn/mat/gpgga'

  if options.pltdir == 'rsn':
    options.pltdir = '/data/rsn/plots/gpgga'

  ifiles=[]
  for arg in args:
    for ifile in glob.glob(arg):
      ifiles.append(ifile)

  for ifile in sorted(ifiles):
    print('ifile=',ifile)
    (idir,iname) = os.path.split(ifile)
    iname = os.path.splitext(iname)[0]

    try:
      ifd = open(ifile,'r')
    except:
      print('ifile does not exist')
      continue

    lineno=0
    uxt  = []
    lat = []
    lon = []

    for linein in ifd:
      while len(linein)>1 and (linein[-1] == '\r' or linein[-1] == '\n'):
        linein = linein[0:-1]
      if len(linein)<1:
        continue

      lineno += 1

      x = linein.split(',')
      if len(x) != 17:
        print('wrong len: linein=',linein)
        sys.exit(1)
      if x[2] != "$GPGGA":
        print('not $GPGGA: linein=',linein)
        sys.exit(1)

      ymdhms = x[0]+'_'+x[1][0:8]
      try:
        secs = timegm(strptime(ymdhms,'%m/%d/%Y_%H:%M:%S'))
      except:
        print('cannot decode time from linein=',linein)
        print('ymdhms=',ymdhms)
        sys.exit(1)

      latd = int(x[4][0:2])
      latm = float(x[4][2:])
      lath = x[5]
      lond = int(x[6][0:3])
      lonm = float(x[6][3:])
      lonh = x[7]

      latx = latd + latm / 60.0
      if lath == 'S':
        latx *= -1.0
      lonx = lond + lonm / 60.0
      if lonh == 'W':
        lonx *= -1.0

      uxt.append (secs)
      lat.append(latx)
      lon.append(lonx)

    ifd.close()

    mlt = np.array(uxt,dtype='double') / 86400 + 719529 - 366
    lat = np.array(lat,dtype='double')
    lon = np.array(lon,dtype='double')

    # HPIES 1AB -- 2014-08-07 1958Z 45 49.4371 129 45.5859
    # HPIES 2SB -- 2014-08-27 1550Z 44 31.3346 125 22.8260
    lat0 =   44.0 + 31.3346 / 60.0
    lon0 = -125.0 - 22.8260 / 60.0

    xsf = np.cos(lat0 * np.pi / 180.0)

    y = (lat - lat0) * 1852.0 * 60.0
    x = (lon - lon0) * 1852.0 * 60.0 * xsf
    
    dist = np.sqrt(x**2 + y**2) / 1000.0

    if options.matdir:
      matdict = {'ifile':ifile,'uxt':uxt,'lat':lat,'lon':lon}
      writemat(options.matdir,iname,matdict)

    if options.pltdir:

      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()

      fig.suptitle(iname + '\nLatitude & Longitude')

      axlist = []

      ax = fig.add_subplot(3,1,1)
      axlist.append(ax)
      ax.plot_date(mlt,lat,'b-')
      ax.grid(True)
      plt.ylabel('Lat, deg')

      ax = fig.add_subplot(3,1,2)
      axlist.append(ax)
      ax.plot_date(mlt,lon,'b-')
      ax.grid(True)
      plt.ylabel('Lon, deg')

      ax = fig.add_subplot(3,1,3)
      axlist.append(ax)
      ax.plot_date(mlt,dist,'b-')
      ax.grid(True)
      ax.set_ylim(-1.0,50.0)
      plt.ylabel('Dist, km')

      fig.autofmt_xdate()
      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, \
        wspace=0.4, hspace=0.2)
      fix_xdates(axlist,1)

      writepdf(options.pltdir,iname)

if __name__ == '__main__':
  mymain()
