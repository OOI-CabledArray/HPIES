#! /usr/bin/env python2
# coher.py -- coherence between hef and magobs

"""
rsync -a ohm:~dunlap/magobs/usgs/{NEW,FRN}  ~dunlap/magobs/usgs
./magobs.py ~/magobs/usgs/NEW/new2014091[2345]vsec.sec -p rsn -m rsn
./coher.py --tfft 8 --rsn --coh --ts /data/rsn/mat/hprpro/HPIES-20140912-1AB-hprpro.mat
==>
fileno= 1 ifile= /data/rsn/mat/hprpro/HPIES-20140912-1AB-hprpro.mat
runid= 1AB
oname= HPIES-20140912-1AB
tsamp= 0.0769726247999 hours, fsamp= 12.991631799 cph
nfft= 100 tfft= 7.69726247999 hrs
magfile= /data/rsn/mat/magobs/new20140912vsec.mat
len(H)= 8640
pdffile= /data/rsn/plots/coher/HPIES-20140912-1AB-ts.pdf
pdffile= /data/rsn/plots/coher/HPIES-20140912-1AB-coh-fh-e1.pdf
pdffile= /data/rsn/plots/coher/HPIES-20140912-1AB-coh-fh-e2.pdf
pdffile= /data/rsn/plots/coher/HPIES-20140912-1AB-coh-fz-e1.pdf
pdffile= /data/rsn/plots/coher/HPIES-20140912-1AB-coh-fz-e2.pdf
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
from spectrajhd import nfft235, autospec, twospec
from time import strptime

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# rc('text', usetex=True)

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
    help="time series EF & MF")

  parser.add_option("--spec",
    action="store_true", dest="plt_spec", default=False,
    help="autospectra of EF & MF")

  parser.add_option("--coh",
    action="store_true", dest="plt_coh", default=False,
    help="coherence between EF & MF")

  parser.add_option("--rsn",
    action="store_true", dest="rsn", default=False,
    help="sets RSN defaults")

  parser.add_option("--okmc",
    action="store_true", dest="okmc", default=False,
    help="sets OKMC defaults")

  parser.add_option("--cat",
    action="store_true", dest="cat", default=False,
    help="concatonate files")

  parser.add_option("--pltdir",
    dest="pltdir", default=None,
    help="plotting directory")

  parser.add_option("--tfft", dest="tfft", default=8.0, type="float",
    help="time of each fft, hours")

  (options, args) = parser.parse_args()

  if options.rsn and options.okmc:
    print('both --rsn and --okmc cannot both be set')
    parser.print_help()
    sys.exit(1)

  if (not options.rsn) and (not options.okmc):
    print('either --rsn or --okmc must be set')
    parser.print_help()
    sys.exit(1)

  if options.pltdir == None:
    if options.rsn:
      options.pltdir = '/data/rsn/plots/coher'
    if options.okmc:
      options.pltdir = '/data/okmc/plots/coher'

# print(args)

  ifiles=[]
  for arg in args:
    for ifile in glob.glob(arg):
      ifiles.append(ifile)

  fileno=0
  for ifile in sorted(ifiles):
    fileno += 1
    print('fileno=',fileno,'ifile=',ifile)
    (idir,bname) = os.path.split(ifile)
    iname = os.path.splitext(bname)[0]

    # assume hprpro name style
    x = iname.split('-')
#   print('x=',x)
    oname = "{0}-{1}-{2}".format(x[0],x[1],x[2])
#   print('oname=',oname)

#   print("strptime(x[1],'%Y%m%d')=",strptime(x[1],'%Y%m%d'))
#   print("timegm(strptime(x[1],'%Y%m%d'))=",timegm(strptime(x[1],'%Y%m%d')))
    if options.cat:
      if fileno == 1:
        uxtfn1 = timegm(strptime(x[1],'%Y%m%d'))
    else:
      uxtfn1 = timegm(strptime(x[1],'%Y%m%d'))

    if options.cat:
      HEF = loadmat(ifile)
      if fileno == 1:
        uxthef = np.array([])
        uxtmin = np.array([])
        uxtmax = np.array([])
        abm    = np.array([])
        e1o    = np.array([])
        e2o    = np.array([])
      uxthef = np.append(uxthef, HEF['dem_uxt'][0])
      uxtmin = np.append(uxtmin, HEF['dem_uxtmin'][0])
      uxtmax = np.append(uxtmax, HEF['dem_uxtmax'][0])
      abm    = np.append(abm   , HEF['dem_abm'])
      e1o    = np.append(e1o   , HEF['dem_e1_amp'][0])
      e2o    = np.append(e2o   , HEF['dem_e2_amp'][0])
      runid  = HEF['runid'][0]
      if fileno < len(ifiles):
        continue
    else:
      HEF = loadmat(ifile)
      uxthef = HEF['dem_uxt'][0]
      uxtmin = HEF['dem_uxtmin'][0]
      uxtmax = HEF['dem_uxtmax'][0]
      abm    = HEF['dem_abm']
      e1o    = HEF['dem_e1_amp'][0]
      e2o    = HEF['dem_e2_amp'][0]
      runid  = HEF['runid'][0]

    print('runid=',runid)
    print('oname=',oname)

    toff = uxtfn1 - uxthef[0]
    print('toff=',toff)
    if uxthef[0] < uxtfn1:
      uxthef += toff

    mlt = uxthef / 86400.0 + 719529 - 366

    tsamp = (mlt[-1] - mlt[0]) / (len(mlt)-1) * 24
    tsamp = tsamp * 2 # because of ja,jb splitting e1o into alternate parts
    fsamp = 1.0 / tsamp              # cycles/hour

    nfft = np.round(options.tfft / tsamp)
    nfft = nfft235(int(nfft/2))*2
    tfft = tsamp * nfft

    print('tsamp=',tsamp,'hours,','fsamp=',fsamp,'cph')
    print('nfft=',nfft,'tfft=',tfft,'hrs')

    if nfft < 4:
      print('len(mlt)=',len(mlt),'nfft=',nfft,'too small -- skipped')
      sys.exit(1)

    ja = np.nonzero(abm == 'a')[0]
    jb = np.nonzero(abm == 'b')[0]

    e1a    = e1o[ja]
    e1b    = e1o[jb]
    e2a    = e2o[ja]
    e2b    = e2o[jb]

    uxta = uxthef[ja]
    uxtb = uxthef[jb]

    mlta = uxta / 86400.0 + 719529 - 366
    mltb = uxtb / 86400.0 + 719529 - 366

    # pick data from mag files with matching times
    uxtf1 = np.floor(np.min(uxthef) / 86400) * 86400
    uxtf2 = np.ceil(np.max(uxthef) / 86400) * 86400
    print('np.min(uxthef)=',np.min(uxthef))
    print('np.max(uxthef)=',np.max(uxthef))
    print('uxtf1=',uxtf1)
    print('uxtf2=',uxtf2)
    ndays = 0
    H = np.array([])
    Z = np.array([])
    uxtmag = np.array([])
    for uxtf in np.arange(uxtf1,uxtf2,86400):
      ndays += 1
      pytf = datetime(1970,1,1) + timedelta(0,uxtf)
      ymdf = pytf.strftime('%Y%m%d')
      print('ymdf=',ymdf)
      if options.okmc:
        magfile = '/data/okmc/mat/magobs/gua' + ymdf + 'vsec.mat'
      if options.rsn:
        magfile = '/data/rsn/mat/magobs/new' + ymdf + 'vsec.mat'
      print('magfile=',magfile)

      # extract magnetic data
      try:
        MAG = loadmat(magfile)
      except:
        print('magfile=',magfile,'does not exist -- exit')
        sys.exit()
      uxtmag = np.append(uxtmag,MAG['uxt'][0])
      H = np.append(H,MAG['H'][0])
      Z = np.append(Z,MAG['Z'][0])

    print('len(H)=',len(H))

    # select fh and fz which go with HEF times
    fha = np.tile(np.nan,len(ja))
    fza = np.tile(np.nan,len(ja))
    for i in range(len(ja)):
      ji = ja[i]
      j = np.nonzero((uxtmin[ji] <= uxtmag) & (uxtmag <= uxtmax[ji]))[0]
      if len(j) > 0:
        fha[i] = np.mean(H[j])
        fza[i] = np.mean(Z[j])
    if not np.isfinite(fha.all()):
      print('some fha,fza missing... exiting')
      sys.exit(1)

#   print('len(fha)=',len(fha))

    fhb = np.tile(np.nan,len(jb))
    fzb = np.tile(np.nan,len(jb))
    for i in range(len(jb)):
      ji = jb[i]
      j = np.nonzero((uxtmin[ji] <= uxtmag) & (uxtmag <= uxtmax[ji]))[0]
      if len(j) > 0:
        fhb[i] = np.mean(H[j])
        fzb[i] = np.mean(Z[j])
    if not np.isfinite(fhb.all()):
      print('some fhb,fzb missing... exiting')
      sys.exit(1)

#   print('len(fhb)=',len(fhb))

    if options.plt_ts:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()

      pltname = oname + '-ts'
      if ndays > 1:
        pltname += '-{0:02d}d'.format(ndays)

    # fig.suptitle(pltname + '\n' + 'Electric Field (' + r'$\micro$' + 'V) & Magnetic Field (nT)')
      fig.suptitle(pltname + '\n' + 'Electric Field (' + 'u' + 'V) & Magnetic Field (nT), red=A, blue=B')

      ym = 10.0
      fhm = np.round(np.median(fha)/ym)*ym
      fzm = np.round(np.median(fza)/ym)*ym

      ylim_ef = [-25.0,25.0]
      ylim_mf = [-150.0,150.0]

      axlist = []

      ax = fig.add_subplot(4,1,1)
      axlist.append(ax)
      ax.plot_date(mlta,e1a,'r-')
      ax.plot_date(mltb,e1b,'b-')
      ax.set_ylim(ylim_ef)
      ax.grid(True)
#     plt.ylabel(r'$e1o \micro V$')
      plt.ylabel('e1o, ' + 'u' + 'V')

      ax = fig.add_subplot(4,1,2)
      axlist.append(ax)
      ax.plot_date(mlta,e2a,'r-')
      ax.plot_date(mltb,e2b,'b-')
      ax.set_ylim(ylim_ef)
      ax.grid(True)
#     plt.ylabel('$e2o, ' + r'$\micro$' + 'V')
      plt.ylabel('e2o, ' + 'u' + 'V')

      ax = fig.add_subplot(4,1,3)
      axlist.append(ax)
      ax.plot_date(mlta,fha-fhm,'r-')
      ax.plot_date(mltb,fhb-fhm,'b-')
      ax.set_ylim(ylim_mf)
      ax.grid(True)
      plt.ylabel('fh - {0:.0f}'.format(fhm))

      ax = fig.add_subplot(4,1,4)
      axlist.append(ax)
      ax.plot_date(mlta,fza-fzm,'r-')
      ax.plot_date(mltb,fzb-fzm,'b-')
      ax.set_ylim(ylim_mf)
      ax.grid(True)
      plt.ylabel('fz - {0:.0f}'.format(fzm))

      fig.autofmt_xdate()
      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, \
        wspace=0.4, hspace=0.2)
      fix_xdates(axlist,1)

      writepdf(options.pltdir, pltname)

    if options.plt_spec:

      Se1a = autospec(e1a,nfft,fsamp,'mean')
      Se1b = autospec(e1b,nfft,fsamp,'mean')
      Se2a = autospec(e2a,nfft,fsamp,'mean')
      Se2b = autospec(e2b,nfft,fsamp,'mean')
      Sfha = autospec(fha,nfft,fsamp,'mean')
      Sfza = autospec(fza,nfft,fsamp,'mean')
      Sfhb = autospec(fhb,nfft,fsamp,'mean')
      Sfzb = autospec(fzb,nfft,fsamp,'mean')
      print('Se1a.npie=',Se1a.npie)

      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()

      pltname = oname + '-mag-spec'
      if ndays > 1:
        pltname += '-{0:02d}d'.format(ndays)

      fig.suptitle(pltname + '\ntfft={0:.1f} hrs'.format(tfft) + ', npie={0}'.format(Se1o.npie))

      axlist = []
      ylim = [-60,30]

      ax = fig.add_subplot(4,1,1)
      axlist.append(ax)
      ax.plot(Se1a.freq,10.0*mylog10(Se1a.spec),'r-')
      ax.plot(Se1b.freq,10.0*mylog10(Se1b.spec),'b-')
      ax.set_xscale('log')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('Se1o, dB')

      ax = fig.add_subplot(4,1,2)
      axlist.append(ax)
      ax.plot(Se2b.freq,10.0*mylog10(Se2b.spec),'r-')
      ax.plot(Se2b.freq,10.0*mylog10(Se2b.spec),'b-')
      ax.set_xscale('log')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('Se2o, dB')

      ax = fig.add_subplot(4,1,3)
      axlist.append(ax)
      ax.plot(Sfha.freq,10.0*mylog10(Sfha.spec),'r-')
      ax.plot(Sfhb.freq,10.0*mylog10(Sfhb.spec),'b-')
      ax.set_xscale('log')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('Sfh, dB')

      ax = fig.add_subplot(4,1,4)
      axlist.append(ax)
      ax.plot(Sfza.freq,10.0*mylog10(Sfza.spec),'r-')
      ax.plot(Sfzb.freq,10.0*mylog10(Sfzb.spec),'b-')
      ax.set_ylim(ylim)
      ax.grid(True)
      plt.ylabel('Sfz, dB')
      ax.set_xscale('log')

      plt.xlabel('Frequency, cycles/hour')

      writepdf(options.pltdir, pltname)

    if options.plt_coh:
      
      witch_list = ['fh-e1','fh-e2','fz-e1','fz-e2']
      
      for witch in witch_list:
        if   witch == 'fh-e1':
          xa = fha
          ya = e1a
          xb = fhb
          yb = e1b
        elif witch == 'fh-e2':
          xa = fha
          ya = e2a
          xb = fhb
          yb = e2b
        elif witch == 'fz-e1':
          xa = fza
          ya = e1a
          xb = fzb
          yb = e1b
        elif witch == 'fz-e2':
          xa = fza
          ya = e2a
          xb = fzb
          yb = e2b
        else:
          print('unknown witch=',witch,'exiting')
          sys.exit(1)

        Sa = twospec(xa,ya,nfft,fsamp)
        Sb = twospec(xb,yb,nfft,fsamp)

        fig = plt.figure(num=1,figsize=(10, 7))
        fig.clf()

        pltname = oname + '-coh-' + witch
        if ndays > 1:
          pltname += '-{0:02d}d'.format(ndays)

        fig.suptitle(pltname + \
            '\ntfft={0:.1f} hrs'.format(tfft) + ', npie={0}'.format(Sa.npie) + \
            '\nred:A, blu:B')

        axlist = []
        ylim = [-60,30]

        ax = fig.add_subplot(2,1,1)
        axlist.append(ax)
        ax.plot(Sa.freq,Sa.rxy,'r-')
        ax.plot(Sb.freq,Sb.rxy,'b-')
        ax.set_xscale('log')
        ax.set_ylim([0,1])
        ax.grid(True)
        plt.ylabel('Rxy^2')

        ax = fig.add_subplot(2,1,2)
        axlist.append(ax)
        ax.plot(Sa.freq,Sa.pxy * 180 / np.pi,'r-')
        ax.plot(Sb.freq,Sb.pxy * 180 / np.pi,'b-')
        ax.set_xscale('log')
        ax.set_ylim([-180,180])
        ax.set_yticks([-180,-90,0,90,180])
        ax.grid(True)
        plt.ylabel('Pxy^2')

        plt.xlabel('Frequency, cycles/hour')

        writepdf(options.pltdir, pltname)


if __name__ == '__main__':
  main()
