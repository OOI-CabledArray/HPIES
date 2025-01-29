#! /usr/bin/env python2
# rotcoh.py -- coherence between pairs of hefs

"""
D=20140831
F=/data/rsn/mat/hprpro/HPIES-"$D"T0000-1AB-hprpro.mat
G=/data/rsn/mat/hprpro/HPIES-"$D"T0000-2SB-hprpro.mat
O=HPIES-"$D"-1AB-2SB
./rotcoh.py --ts --auto --rot $F $G --opref $O
"""

from __future__ import print_function

import sys
import os
from optparse import OptionParser
import numpy as np
import matplotlib
matplotlib.use('Agg') # needed when used with cron 
                      # or no DISPLAY 
                      # or slow remote connection
from scipy.io import loadmat
import matplotlib.pyplot as plt
import collections
from calendar import timegm
import glob
from datetime import datetime, timedelta
import matplotlib.dates as mdates

from spectrajhd import nfft235, autospec, twospec, rotcoh

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

def writepdf(pdfdir, namewoext):
  mymkdir(pdfdir)
  pdffile = os.path.join(pdfdir,namewoext+'.pdf')
  print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pdfdir)

def mymkdir(mydir):
  try: 
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      sys.exit(1)

def uxt2pyt(uxt):
  return datetime(1970,1,1) + timedelta(0,uxt)
  
def uxt2str(uxt):
  pyt = datetime(1970,1,1) + timedelta(0,uxt)
  return pyt.strftime('%Y-%m-%d %H:%M:%S')

def mylog10(x):
  j = np.nonzero(x > 0.0)
  y = np.tile(np.nan,len(x))
  y[j] = np.log10(x[j])
  return y
  

def main():

  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose",
    action="store_true", dest="verbose", default=False,
    help="print status messages to stdout")

  parser.add_option("--ts",
    action="store_true", dest="plt_ts", default=False,
    help="time series")

  parser.add_option("--auto",
    action="store_true", dest="plt_auto", default=False,
    help="auto-spectra")

  parser.add_option("--rot",
    action="store_true", dest="plt_rot", default=False,
    help="rotary coherence")

  parser.add_option("-p", "--pltdir",
    dest="pltdir", default='/data/rsn/plots/rotcoh',
    help="plotting directory")

  parser.add_option("-o", "--opref",
    dest="opref", default=None,
    help="ouput name prefix")

  parser.add_option("--tfft", dest="tfft", default=10.0, type="float",
    help="time of each fft, hours")

  (options, args) = parser.parse_args()

  if len(args) != 2:
    print('must have two files')
    sys.exit(1)

  file1 = args[0]
  file2 = args[1]

  print('file1=',file1)
  print('file2=',file2)

  (idir1,bname1) = os.path.split(file1)
  name1 = os.path.splitext(bname1)[0]

  (idir2,bname2) = os.path.split(file2)
  name2 = os.path.splitext(bname2)[0]

  # default output filename prefix
  if options.opref:
    opref = options.opref
  else:
    # shorter opref for RSN file name patterns
    s1 = name1.split('-')
    s2 = name2.split('-')
    if len(s1)==4 and len(s2)==4 and \
     s1[0]==s2[0] and s1[1]==s2[1] and s1[3]==s2[3]:
      opref = '{0}-{1}-{2}-{3}-{4}'.format(s1[0],s1[1],s1[2],s2[2],s1[3])
    else:
      opref = name1 + '_' + name2

  try:
    F1 = loadmat(file1)
  except:
    print('cannot load file1=',file1)
    sys.exit(1)

  try:
    F2 = loadmat(file2)
  except:
    print('cannot load file2=',file2)
    sys.exit(1)

  # x1 = file1.split('-'); D1 = x1[1]; L1 = x1[2]
  # x2 = file2.split('-'); D2 = x2[1]; L2 = x2[2]

  L1 = F1['runid'][0]
  L2 = F2['runid'][0]

  uxtmin = F1['uxtmin'][0]
  uxtmax = F1['uxtmax'][0]

  uxt1 = F1['uxt'][0]
  uxt2 = F2['uxt'][0]

  # common time interval
  dt = 240.0
  uxtmin = np.max([uxt1[0],uxt2[0]])
  uxtmax = np.min([uxt1[-1],uxt2[-1]])

  uxt = np.arange(np.ceil(uxtmin/dt)*dt, np.floor(uxtmax/dt)*dt, dt)
  mlt = uxt / 86400.0 + 719529 - 366

  tsamp = np.mean(np.diff(mlt))*24 # hours
  fsamp = 1.0 / tsamp              # cycles/hour

  nfft = np.round(options.tfft / tsamp)
  nfft = nfft235(int(nfft/2))*2
  if nfft < 4:
    print('len(mlt)=',len(mlt),'nfft=',nfft,'too small -- skipped')
    sys.exit(1)

  tfft = tsamp * nfft

  print('tsamp=',tsamp,'hours,','fsamp=',fsamp,'cph')
  print('nfft=',nfft,'tfft=',tfft)

  ja1 = np.nonzero(F1['abm'] == 'a')[0]
  jb1 = np.nonzero(F1['abm'] == 'b')[0]
  ja2 = np.nonzero(F2['abm'] == 'a')[0]
  jb2 = np.nonzero(F2['abm'] == 'b')[0]

  # first HEF
  e1a1 = np.interp(uxt, uxt1[ja1], F1['e1o'][0][ja1])
  e1b1 = np.interp(uxt, uxt1[jb1], F1['e1o'][0][jb1])
  e2a1 = np.interp(uxt, uxt1[ja1], F1['e2o'][0][ja1])
  e2b1 = np.interp(uxt, uxt1[jb1], F1['e2o'][0][jb1])

  # second HEF
  e1a2 = np.interp(uxt, uxt2[ja2], F2['e1o'][0][ja2])
  e1b2 = np.interp(uxt, uxt2[jb2], F2['e1o'][0][jb2])
  e2a2 = np.interp(uxt, uxt2[ja2], F2['e2o'][0][ja2])
  e2b2 = np.interp(uxt, uxt2[jb2], F2['e2o'][0][jb2])

  if options.plt_ts:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltname = opref + '-ts'

    fig.suptitle(pltname + '\nred=A, blu=B')

    ylim = [-15.0,15.0]

    axlist = []

    ax = fig.add_subplot(4,1,1)
    axlist.append(ax)
    ax.plot_date(mlt,e1a1,'r-')
    ax.plot_date(mlt,e1b1,'b-')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('e1a1,e1b1')

    ax = fig.add_subplot(4,1,2)
    axlist.append(ax)
    ax.plot_date(mlt,e2a1,'r-')
    ax.plot_date(mlt,e2b1,'b-')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('e2a1,e2b1')

    ax = fig.add_subplot(4,1,3)
    axlist.append(ax)
    ax.plot_date(mlt,e1a2,'r-')
    ax.plot_date(mlt,e1b2,'b-')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('e1a2,e1b2')

    ax = fig.add_subplot(4,1,4)
    axlist.append(ax)
    ax.plot_date(mlt,e2a2,'r-')
    ax.plot_date(mlt,e2b2,'b-')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('e2a2,e2b2')

    fig.autofmt_xdate()
    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, \
      wspace=0.4, hspace=0.2)
    fix_xdates(axlist,1)

    writepdf(options.pltdir, pltname)

  if options.plt_auto:

    # file1
    Se1a1 = autospec(e1a1,nfft,fsamp,'mean')
    Se1b1 = autospec(e1b1,nfft,fsamp,'mean')
    Se2a1 = autospec(e2a1,nfft,fsamp,'mean')
    Se2b1 = autospec(e2b1,nfft,fsamp,'mean')

    # file2
    Se1a2 = autospec(e1a2,nfft,fsamp,'mean')
    Se1b2 = autospec(e1b2,nfft,fsamp,'mean')
    Se2a2 = autospec(e2a2,nfft,fsamp,'mean')
    Se2b2 = autospec(e2b2,nfft,fsamp,'mean')

    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltname = opref + '-auto'

    fig.suptitle(pltname + \
      '\nnfft={0}'.format(nfft) + \
      ', npie={0}'.format(Se1a1.npie))

    axlist = []
    ylim = [-60,30]

    ax = fig.add_subplot(4,1,1)
    axlist.append(ax)
    ax.plot(Se1a1.freq,10.0*mylog10(Se1a1.spec),'r-')
    ax.plot(Se1b1.freq,10.0*mylog10(Se1b1.spec),'b-')
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('Se11, dB')

    ax = fig.add_subplot(4,1,2)
    axlist.append(ax)
    ax.plot(Se2a1.freq,10.0*mylog10(Se2a1.spec),'r-')
    ax.plot(Se2b1.freq,10.0*mylog10(Se2b1.spec),'b-')
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('Se21, dB')

    ax = fig.add_subplot(4,1,3)
    axlist.append(ax)
    ax.plot(Se1a2.freq,10.0*mylog10(Se1a2.spec),'r-')
    ax.plot(Se1b2.freq,10.0*mylog10(Se1b2.spec),'b-')
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('Se12, dB')

    ax = fig.add_subplot(4,1,4)
    axlist.append(ax)
    ax.plot(Se2a2.freq,10.0*mylog10(Se2a2.spec),'r-')
    ax.plot(Se2b2.freq,10.0*mylog10(Se2b2.spec),'b-')
    ax.set_xscale('log')
    ax.set_ylim(ylim)
    ax.grid(True)
    plt.ylabel('Se22, dB')

    plt.xlabel('Frequency, cycles/hour')

    writepdf(options.pltdir, pltname)

  if options.plt_rot:

    for witch in [\
        'Raaaa','Raaab','Raaba','Raabb',\
        'Rabaa','Rabab','Rabba','Rabbb',\
        'Rbaaa','Rbaab','Rbaba','Rbabb',\
        'Rbbaa','Rbbab','Rbbba','Rbbbb' ]:

      if witch == 'Raaaa':
        R = rotcoh(e1a1,e2a1,e1a2,e2a2,nfft,fsamp)
        ttl = L1+'(e1a,e2a) '+L2+'(e1a,e2a)'
      elif witch == 'Raaab':
        R = rotcoh(e1a1,e2a1,e1a2,e2b2,nfft,fsamp)
        ttl = L1+'(e1a,e2a) '+L2+'(e1a,e2b)'
      elif witch == 'Raaba':
        R = rotcoh(e1a1,e2a1,e1b2,e2a2,nfft,fsamp)
        ttl = L1+'(e1a,e2a) '+L2+'(e1b,e2a)'
      elif witch == 'Raabb':
        R = rotcoh(e1a1,e2a1,e1b2,e2b2,nfft,fsamp)
        ttl = L1+'(e1a,e2a) '+L2+'(e1b,e2b)'

      elif witch == 'Rabaa':
        R = rotcoh(e1a1,e2b1,e1a2,e2a2,nfft,fsamp)
        ttl = L1+'(e1a,e2b) '+L2+'(e1a,e2a)'
      elif witch == 'Rabab':
        R = rotcoh(e1a1,e2b1,e1a2,e2b2,nfft,fsamp)
        ttl = L1+'(e1a,e2b) '+L2+'(e1a,e2b)'
      elif witch == 'Rabba':
        R = rotcoh(e1a1,e2b1,e1b2,e2a2,nfft,fsamp)
        ttl = L1+'(e1a,e2b) '+L2+'(e1b,e2a)'
      elif witch == 'Rabbb':
        R = rotcoh(e1a1,e2b1,e1b2,e2b2,nfft,fsamp)
        ttl = L1+'(e1a,e2b) '+L2+'(e1b,e2b)'

      elif witch == 'Rbaaa':
        R = rotcoh(e1b1,e2a1,e1a2,e2a2,nfft,fsamp)
        ttl = L1+'(e1b,e2a) '+L2+'(e1a,e2a)'
      elif witch == 'Rbaab':
        R = rotcoh(e1b1,e2a1,e1a2,e2b2,nfft,fsamp)
        ttl = L1+'(e1b,e2a) '+L2+'(e1a,e2b)'
      elif witch == 'Rbaba':
        R = rotcoh(e1b1,e2a1,e1b2,e2a2,nfft,fsamp)
        ttl = L1+'(e1b,e2a) '+L2+'(e1b,e2a)'
      elif witch == 'Rbabb':
        R = rotcoh(e1b1,e2a1,e1b2,e2b2,nfft,fsamp)
        ttl = L1+'(e1b,e2a) '+L2+'(e1b,e2b)'

      elif witch == 'Rbbaa':
        R = rotcoh(e1b1,e2b1,e1a2,e2a2,nfft,fsamp)
        ttl = L1+'(e1b,e2b) '+L2+'(e1a,e2a)'
      elif witch == 'Rbbab':
        R = rotcoh(e1b1,e2b1,e1a2,e2b2,nfft,fsamp)
        ttl = L1+'(e1b,e2b) '+L2+'(e1a,e2b)'
      elif witch == 'Rbbba':
        R = rotcoh(e1b1,e2b1,e1b2,e2a2,nfft,fsamp)
        ttl = L1+'(e1b,e2b) '+L2+'(e1b,e2a)'
      elif witch == 'Rbbbb':
        R = rotcoh(e1b1,e2b1,e1b2,e2b2,nfft,fsamp)
        ttl = L1+'(e1b,e2b) '+L2+'(e1b,e2b)'

      else:
        print('unknown witch=',witch)
        sys.exit(1)

      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltname = opref + '-' + witch

      fig.suptitle(pltname + \
        '\nnfft={0}'.format(nfft) + \
        ', npie={0}'.format(R.npie) + \
        ', red:+freq, blu:-freq' + \
        '\n{0}'.format(ttl))

      axlist = []

      ax = fig.add_subplot(4,1,1)
      axlist.append(ax)
      ax.plot(R.F,R.gam2_p12,'r-')
      ax.plot(R.F,R.gam2_n12,'b-')
      ax.set_xscale('log')
      ax.set_ylim([0,1])
      ax.grid(True)
      plt.ylabel('Inner R^2')

      ax = fig.add_subplot(4,1,2)
      axlist.append(ax)
      ax.plot(R.F,R.phi_p12 * 180 / np.pi,'r-')
      ax.plot(R.F,R.phi_n12 * 180 / np.pi,'b-')
      ax.set_xscale('log')
      ax.set_ylim([-180,180])
      ax.set_yticks([-180,-90,0,90,180])
      ax.grid(True)
      plt.ylabel('Inner Phase')

      ax = fig.add_subplot(4,1,3)
      axlist.append(ax)
      ax.plot(R.F,R.lam2_p12,'r-')
      ax.plot(R.F,R.lam2_n12,'b-')
      ax.set_xscale('log')
      ax.set_ylim([0,1])
      ax.grid(True)
      plt.ylabel('Outer R^2')

      ax = fig.add_subplot(4,1,4)
      axlist.append(ax)
      ax.plot(R.F,R.psi_p12 * 180 / np.pi,'r-')
      ax.plot(R.F,R.psi_n12 * 180 / np.pi,'b-')
      ax.set_xscale('log')
      ax.set_ylim([-180,180])
      ax.set_yticks([-180,-90,0,90,180])
      ax.grid(True)
      plt.ylabel('Outer Phase')

      plt.xlabel('Frequency, cycles/hour')

      writepdf(options.pltdir, pltname)
  print()

if __name__ == '__main__':
  main()
