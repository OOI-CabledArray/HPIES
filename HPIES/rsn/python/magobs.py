#!/usr/bin/env python2
# magobs.py -- decode, plot and process 1-s ~/magobs/usgs/...

# rsync -a ohm:~dunlap/magobs/usgs/{NEW,FRN}  ~dunlap/magobs/usgs
# ./magobs.py ~/magobs/usgs/NEW/new20140831vsec.sec -p rsn -m rsn
# ./magobs.py ~/magobs/usgs/NEW/new2014????vsec.sec -p rsn -m rsn
# ./magobs.py /data/okmc/mag/one-second/txt/GUA/gua20120703vsec.sec -p okmc -m okmc
# ./magobs.py /data/okmc/mag/one-second/txt/GUA/gua201*vsec.sec.gz -p okmc -m okmc

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

  parser.add_option("-p", "--pltdir",
    dest="pltdir", default=None,
    help="plotting directory")

  parser.add_option('-m', '--matdir',
    dest='matdir', default=None,
    help='matlab output directory')

  (options, args) = parser.parse_args()

  if options.matdir == 'rsn':
    options.matdir = '/data/rsn/mat/magobs'

  if options.pltdir == 'rsn':
    options.pltdir = '/data/rsn/plots/magobs'

  if options.matdir == 'okmc':
    options.matdir = '/data/okmc/mat/magobs'

  if options.pltdir == 'okmc':
    options.pltdir = '/data/okmc/plots/magobs'

  ifiles=[]
  for arg in args:
    for ifile in glob.glob(arg):
      ifiles.append(ifile)

  for ifile in sorted(ifiles):
    print('ifile=',ifile)
    (idir,iname) = os.path.split(ifile)

    oname = iname
    if oname[-3:] == '.gz':
      oname = oname[:-3]
    if oname[-4:] == '.sec':
      oname = oname[:-4]
    print('oname=',oname)

    obsstn = oname[:3]

    try:
      if ifile[-3:] == '.gz':
        ifd = gzip.open(ifile,'r')
      else:
        ifd = open(ifile,'r')
    except:
      print('cannot open ifile')
      continue

    lineno=0
    uxt  = []
    H = []
    D = []
    Z = []
    F = []

    for linein in ifd:
      while len(linein)>1 and (linein[-1] == '\r' or linein[-1] == '\n'):
        linein = linein[0:-1]
      if len(linein)<1:
        continue

      lineno += 1

      x = linein.split()
      if x[-1] == '|':  # header
        continue
      if len(x) != 7:
        print('wrong len: linein=',linein)
        sys.exit()
      if x[0][0:3] != '201':
        print('not a year: linein=',linein)
        sys.exit()

      ymdhms = x[0]+'T'+x[1][0:8]
      try:
        secs = timegm(strptime(ymdhms,'%Y-%m-%dT%H:%M:%S'))
      except:
        print('cannot decode time from linein=',linein)
        print('ymdhms=',ymdhms)
        sys.exit(1)

      uxt.append (secs)
      H.append(float(x[3]))
      D.append(float(x[4]))
      Z.append(float(x[5]))
      F.append(float(x[6]))

    ifd.close()

    mlt = np.array(uxt,dtype='double') / 86400 + 719529 - 366
    H = np.array(H,dtype='double')
    D = np.array(D,dtype='double')
    Z = np.array(Z,dtype='double')
    F = np.array(F,dtype='double')

    i = np.nonzero(H == 99999.0)[0]
    j = np.nonzero(H != 99999.0)[0]
    H[i] = np.interp(i,j,H[j])
    print('intepolated',len(i),'H values')
    i = np.nonzero(Z == 99999.0)[0]
    j = np.nonzero(Z != 99999.0)[0]
    Z[i] = np.interp(i,j,Z[j])
    print('intepolated',len(i),'Z values')

    # filter by factor of m
    m = 10
    n = int(len(H)/m)  # new length
    Hf = np.mean(np.reshape(H[0:m*n],(m,n),order='F'),0)
    Zf = np.mean(np.reshape(Z[0:m*n],(m,n),order='F'),0)
    uxtf = np.mean(np.reshape(uxt[0:m*n],(m,n),order='F'),0)

    if options.matdir:
      matdict = {'ifile':ifile,'oname':oname,'obsstn':obsstn,'uxt':uxtf,'H':Hf,'Z':Zf}
      writemat(options.matdir,oname,matdict)

    if options.pltdir:

      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()

      fig.suptitle(oname + '\nHorizontal & Vertical Magnetic Field')

      axlist = []
      ylim = [-100.0,100.0]
      ym = 10.0
      Hm = np.round(np.median(H)/ym)*ym
      Dm = np.round(np.median(D)/ym)*ym
      Zm = np.round(np.median(Z)/ym)*ym
      Fm = np.round(np.median(F)/ym)*ym

      ax = fig.add_subplot(2,1,1)
      axlist.append(ax)
      ax.plot_date(mlt,H-Hm,'b-')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('H - {0:.0f} nT'.format(Hm))

      ax = fig.add_subplot(2,1,2)
      axlist.append(ax)
      ax.plot_date(mlt,Z-Zm,'b-')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('Z - {0:.0f} nT'.format(Zm))

      fig.autofmt_xdate()
      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, \
        wspace=0.4, hspace=0.2)
      fix_xdates(axlist,1)
      writepdf(options.pltdir, oname)


if __name__ == '__main__':
  mymain()
