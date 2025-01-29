#!/usr/bin/env python2
# hprdif.py -- look at differences between A & B channels

"""
./hprdif.py /data/rsn/mat/hprpro/HPIES-20141222-2SB-hprpro.mat
"""

from __future__ import print_function

import os
from sys import exit
from optparse import OptionParser
import math
import numpy as np
import matplotlib
# matplotlib.use('Agg') # for non-gui plotting
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import collections
from scipy.io import loadmat
from calendar import timegm
#import time
from time import strptime, strftime
from datetime import datetime,timedelta
import scipy.io

def date2uxt(yr,mo,da,hr,mi,se):
  pyt = datetime(yr,mo,da,hr,mi,se)
  dif = pyt - datetime(1970,1,1,0,0,0)
  uxt = dif.days * 86400 + dif.seconds
  return uxt

def writemat(matdir,namewoext,matdict):
  mymkdir(matdir)
  matfile = os.path.join(matdir,namewoext+'.mat')
  print('matfile=',matfile)
  scipy.io.savemat(matfile, matdict, format='4')

def writepdf(pdfdir, namewoext):
  mymkdir(pdfdir)
  pdffile = os.path.join(pdfdir,namewoext+'.pdf')
  print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pdfdir)

def writepng(pngdir, ifile):
  basename = os.path.basename(ifile)
  rootname = os.path.splitext(basename)[0]
  pngfile = os.path.join(options.pngdir,rootname+'.png')
  print('pngfile=',pngfile)
  mymkdir(pngdir)
  plt.savefig(pngfile)
  os.system('/home/dunlap/bin/updateframe.run ' + pngdir)

def mymkdir(mydir):
  try: 
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      sys.exit(1)

def adj_ylim(ax):
  ylim = np.array(ax.get_ylim(),'double')
  adj = (ylim[1] - ylim[0]) * 0.02
  ylim[0] -= adj
  ylim[1] += adj
  ax.set_ylim(ylim)

def blkavg(a,m):
  n = int(len(a) / m)
  return np.mean(np.reshape(a[0:n*m],(m,n),order='F'),0)

def mylim(x,m):
  j = np.nonzero(x >  m)[0]; x[j] =  m
  j = np.nonzero(x < -m)[0]; x[j] = -m
  return x

def main():

  parser = OptionParser(
    usage="%prog [Options] file[s]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose",
    action="count", dest="verbose", default=0,
    help="print status messages to stdout")

  parser.add_option("--pngdir",
    dest="pngdir", default=None,
    help="PNG directory")

  global options
  (options, args) = parser.parse_args()

  if len(args) == 0:
    parser.print_help()
    print('no files given')
    exit(1)

  for ifile in args:
    print('ifile=',ifile)

    try:
      M = loadmat(ifile)
    except:
      print('warning: cannot loadmat('+ifile+') -- skipped')
      continue

    uxt = M['dem_uxt'][0]
    abm = M['dem_abm']
    e1  = M['dem_e1_amp'][0]
    e2  = M['dem_e2_amp'][0]

    ja = np.nonzero(abm == 'a')[0]
    jb = np.nonzero(abm == 'b')[0]

    ta = uxt[ja] - uxt[0]
    tb = uxt[jb] - uxt[0]
    e1a = e1[ja]
    e1b = e1[jb]
    e2a = e2[ja]
    e2b = e2[jb]

    e1ai = np.interp(tb, ta, e1a)
    e1bi = np.interp(ta, tb, e1b)
    e2ai = np.interp(tb, ta, e2a)
    e2bi = np.interp(ta, tb, e2b)

    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    fig.suptitle(ifile)

    ylim = [-1.5,+1.5]

    ax = fig.add_subplot(2,1,1)
    plt.plot(ta, e1a-e1bi, label='e1a-e1bi, std={0:.2f}, mean={1:.2f}'.\
      format(np.std(e1a-e1bi),np.mean(e1a-e1bi)))
    plt.plot(tb, e1ai-e1b, label='e1ai-e1b, std={0:.2f}, mean={1:.2f}'.\
      format(np.std(e1ai-e1b),np.mean(e1ai-e1b)))
    plt.legend()
    plt.gca().set_ylim(ylim)
    plt.grid(True)
    plt.ylabel('uV')
    plt.title('e1a-e1b diffs')

    ax = fig.add_subplot(2,1,2)
#   plt.plot(ta, e2a - e2bi, label='std(e2a-e2bi)={0:.2f}'.format(np.std(e2a-e2bi)))
#   plt.plot(tb, e2ai - e2b, label='std(e2ai-e2b)={0:.2f}'.format(np.std(e2ai-e2b)))
    plt.plot(ta, e2a-e2bi, label='e2a-e2bi, std={0:.2f}, mean={1:.2f}'.\
      format(np.std(e2a-e2bi),np.mean(e2a-e2bi)))
    plt.plot(tb, e2ai-e2b, label='e2ai-e2b, std={0:.2f}, mean={1:.2f}'.\
      format(np.std(e2ai-e2b),np.mean(e2ai-e2b)))
    plt.legend()
    plt.gca().set_ylim(ylim)
    plt.grid(True)
    plt.ylabel('uV')
    plt.xlabel('time, s')
    plt.title('e2a-e2b diffs')

    fig.show()

    if options.pngdir:
      writepng(options.pngdir,ifile)
    else:
      raw_input("Press Enter to continue...")


if __name__ == '__main__':
  main()
