#!/usr/bin/env python2
# hprwsr.py -- estimate water switch resistances to 
#              make ocean data from demod data from hprpro.py

"""
./hprwsr.py -s 1AB 
./hprwsr.py -s 2SB 
./hprwsr.py -s H1  
./hprwsr.py -s H2  
./hprwsr.py -s H3  
./hprwsr.py -s H4  
./hprwsr.py -s H5 

for i in H1 H2 H3 H4 H5; do ./hprwsr.py -s $i; done 
for i in 1AB 2SB;        do ./hprwsr.py -s $i; done

"""

from __future__ import print_function

import os
from sys import exit
from optparse import OptionParser
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import collections
from scipy.io import loadmat
from calendar import timegm
#import time
from time import strptime, strftime
from datetime import datetime,timedelta
import scipy.io
import glob

def writemat(matodir,namewoext,matdict):
  mymkdir(matodir)
  matofile = os.path.join(matodir,namewoext+'.mat')
  if options.verbose:
    print('matofile=',matofile)
  scipy.io.savemat(matofile, matdict, format='4')

def writepdf(pltdir, namewoext):
  mymkdir(pltdir)
  pdffile = os.path.join(pltdir,namewoext+'.pdf')
  if options.verbose:
    print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pltdir)

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

def mylim(x,m):
  j = np.nonzero(x >  m)[0]; x[j] =  m
  j = np.nonzero(x < -m)[0]; x[j] = -m
  return x

def main():

  parser = OptionParser(
    usage="%prog [Options]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose",
    action="count", dest="verbose", default=0,
    help="print status messages to stdout")

  parser.add_option("-d","--date",
    dest="date", default=None,
    help="YYYYMMDD, e.g., 20121001")

  parser.add_option("-s","--site",
    dest="site", default=None,
    help="deployment location, e.g., H1,...,H5 or 1AB or 2SB")

  parser.add_option("--do_cal_ylim", dest="do_cal_ylim", 
    action="store_true", default=False,
    help="apply manual ylimit to cal plots")

  parser.add_option("-p", "--make_plots", dest="make_plots", action="store_true", default=False,
    help="make plots")

  global options
  (options, args) = parser.parse_args()

  if options.make_plots:
    do_plt_R1      = True
    do_plt_R2      = True
    do_plt_cal     = True
    do_plt_ocean   = True
  else:
    do_plt_R1      = False
    do_plt_R2      = False
    do_plt_cal     = False
    do_plt_ocean   = False

  do_plt_demod   = False
  do_plt_mot     = False


  tlbl = 'hours'
  tdiv = 3600
  xlim = [0,24]
  redmrk = 'r.-'
  blumrk = 'b.-'
  grnmrk = 'g.-'
  redmrkmot = 'r.-'
  blumrkmot = 'b.-'
  grnmrkmot = 'g.-'
  mrktt = 'b.-'
  mrksz = 1

  redmrk_eftt = redmrk
  blumrk_eftt = blumrk
  grnmrk_eftt = grnmrk
  mrktt_eftt = mrktt
  mrksz_eftt = mrksz

  o = options.site;
  if o=='H1' or o=='H2' or o=='H3' or o=='H4' or o=='H5' :
    deployment = 'OKMC'
    pyidir = '/data/okmc/mat2/hprpro' # hprpro.py tskip=40, twait=10
    pyidir = '/data/okmc/mat/hprpro'  # hprpro.py tskip=15, twait=7
    pyidir = '/data/okmc/mat3/hprpro' # hprpro.py tskip=40, twait=10
    pyidir = '/data/okmc/mat4/hprpro'
    matodir = '/data/okmc/mat4/hprwsr'
    pltdir = '/data/okmc/plots/hprwsr'
    pltdirq = './plots'
  elif o=='1AB' or o=='2SB':
    deployment = 'RSN'
    pyidir = '/data/rsn/mat/hprpro'  # hprpro.py tskip=15, twait=7
    matodir = '/data/rsn/mat/hprwsr'
    pltdir = '/data/rsn/plots/hprwsr'
    pltdirq = './plots'
  else:
    print('unknown site=',options.site)
    exit(1)

  ifiles=[]
  for ifile in glob.glob(pyidir + '/HPIES-*-' + options.site + '-hprpro.mat'):
    ifiles.append(ifile)
  if options.verbose:
      print('len(ifiles)=',len(ifiles))
  if len(ifiles) == 0:
    exit(1)

  for ifile in sorted(ifiles):

    (idir,bname) = os.path.split(ifile)
    iname = os.path.splitext(bname)[0]

    # assume hprpro name style
    sname = iname.split('-')
    try:
      uxtfn = timegm(strptime(sname[1],'%Y%m%d'))
    except:
      print('cannot decode date from ifile=',ifile)
      exit(1)

    pyt = datetime(1970,1,1,0,0,0) + timedelta(0,uxtfn)
    pytstr = pyt.strftime('%Y%m%d')
    pyfile = pyidir+'/HPIES-'+pytstr+'-'+options.site+'-hprpro.mat'
    if options.verbose:
      print('pyfile=',pyfile)

    try:
      PY = loadmat(pyfile)
    except:
      print('warning: cannot loadmat(',pyfile,')')
      continue

    dem_uxt     = PY['dem_uxt'][0]
    dem_abm     = PY['dem_abm']
    dem_e1_amp   = PY['dem_e1_amp'][0]
    dem_e2_amp   = PY['dem_e2_amp'][0]

    if len(PY['cal_uxt']) == 0:
      print('not enough cal data.  Skipped ifile=',ifile)
      continue

    cal_uxt = PY['cal_uxt'][0]
    cal_abu = PY['cal_abu']
    cal_e1a_amp  = PY['cal_e1a_amp'][0]
    cal_e1b_amp  = PY['cal_e1b_amp'][0]
    cal_e1c_amp  = PY['cal_e1c_amp'][0]
    cal_e2a_amp  = PY['cal_e2a_amp'][0]
    cal_e2b_amp  = PY['cal_e2b_amp'][0]
    cal_e2c_amp  = PY['cal_e2c_amp'][0]

    mot_abu  = PY['mot_abu']
    mot_uxt1 = PY['mot_uxt1'][0]
    mot_uxt2 = PY['mot_uxt2'][0]
    mot_cur1 = PY['mot_cur1'][0]
    mot_cur2 = PY['mot_cur2'][0]
    mot_dur1 = PY['mot_dur1'][0]
    mot_dur2 = PY['mot_dur2'][0]

    if len(dem_uxt) == 0:
      print('no data found -- exit')
      exit(1)

    uxt_ref = uxtfn

    if options.verbose >= 2:
      print('len(dem_uxt)=',len(dem_uxt))
      print('len(cal_uxt)=',len(cal_uxt))
      print('len(mot_uxt1)=',len(mot_uxt1))
      print('len(mot_uxt2)=',len(mot_uxt2))

    jm1 = np.searchsorted(mot_uxt1, dem_uxt)
    jm2 = np.searchsorted(mot_uxt2, dem_uxt)

    for i in range(1,len(dem_uxt)):
      j1=jm1[i]
      j2=jm2[i]
      if j1-2 < 0 or j1 >= len(mot_uxt1):
        print('mot_uxt1 probs: i=',i,'j1=',j1,'skipped')
        continue
      if j2-2 < 0 or j2 >= len(mot_uxt1):
        print('mot_uxt1 probs: i=',i,'j2=',j2,'skipped')
        continue

    ja = np.nonzero(dem_abm == 'a')[0]
    jb = np.nonzero(dem_abm == 'b')[0]

    if options.verbose >= 2:
      print('len(dem_abm)=',len(dem_abm))
      print('len(ja)=',len(ja))
      print('len(jb)=',len(jb))


    # apply factor of six to e1b for H3 -- JHD
    if options.site == 'H3':
      sf = 6.0
      print('H3 e1b,e1c sf=',sf)
      dem_e1_amp[jb]  *= sf
      cal_e1b_amp      *= sf
      cal_e1c_amp      *= sf

    # OKMC salt bridge resistances assuming conductivity = 3.5 S/m
    Rec = 95 # from "end+" to "C+" electrode
    Rca = 55 # from "C+" electrode to "A+" electrode
    Rab = 50 # from "A-" electrode to "B+" electrode in back side of water switch
    Rbc = 55 # from "B-" electrode to "C-" electrode
    Rce = 95 # from "C-" electrode to "end-"
    Rx = Rca + Rab + Rbc
    Ry = Rx + Rec + Rce

    # compute water switch resistances from cal data
    # note: e1b, e2b signs were already inverted in hprpro.py
    ediffmin = 1  # microvolts

    # A-pinched
    jak = np.nonzero(cal_abu == 'a')[0]
    tak = cal_uxt[jak] - uxt_ref
    e1a1k = cal_e1a_amp[jak]
    e1b1k = cal_e1b_amp[jak]
    e1c1k = cal_e1c_amp[jak]
    ediff = e1c1k - e1a1k - e1b1k
    j = np.nonzero(ediff < ediffmin)
    ediff[j] = ediffmin
    R1ap = Rx * e1a1k / ediff # R1ap
    R1bu = Rx * e1b1k / ediff # R1bu
    e2a1k = cal_e2a_amp[jak]
    e2b1k = cal_e2b_amp[jak]
    e2c1k = cal_e2c_amp[jak]
    ediff = e2c1k - e2a1k - e2b1k
    j = np.nonzero(ediff < ediffmin)
    ediff[j] = ediffmin
    R2ap = Rx * e2a1k / ediff # R2ap
    R2bu = Rx * e2b1k / ediff # R2bu

    # B-pinched
    jbk = np.nonzero(cal_abu == 'b')[0]
    tbk = cal_uxt[jbk] - uxt_ref
    e1a2k = cal_e1a_amp[jbk]
    e1b2k = cal_e1b_amp[jbk]
    e1c2k = cal_e1c_amp[jbk]
    ediff = e1c2k - e1a2k - e1b2k
    j = np.nonzero(ediff < ediffmin)
    ediff[j] = ediffmin
    R1au = Rx * e1a2k / ediff # R1au
    R1bp = Rx * e1b2k / ediff # R1bp
    e2a2k = cal_e2a_amp[jbk]
    e2b2k = cal_e2b_amp[jbk]
    e2c2k = cal_e2c_amp[jbk]
    ediff = e2c2k - e2a2k - e2b2k
    j = np.nonzero(ediff < ediffmin)
    ediff[j] = ediffmin
    R2au = Rx * e2a2k / ediff # R1au
    R2bp = Rx * e2b2k / ediff # R1bp

    Rpmax = 1e5
    Rumax = 2e3
    f = 1.02
    Rplim = [0,Rpmax*f]
    Rulim = [0,Rumax*f]

    R1ap = mylim(R1ap,Rpmax)
    R1au = mylim(R1au,Rumax)
    R1bu = mylim(R1bu,Rumax)
    R1bp = mylim(R1bp,Rpmax)

    R2ap = mylim(R2ap,Rpmax)
    R2au = mylim(R2au,Rumax)
    R2bu = mylim(R2bu,Rumax)
    R2bp = mylim(R2bp,Rpmax)

    # seconds from uxt_ref
    ta = dem_uxt[ja] - uxt_ref
    tb = dem_uxt[jb] - uxt_ref
    if options.verbose >= 2:
      print('uxt_ref=',uxt_ref)

    dta = np.diff(ta)
    dtb = np.diff(tb)
    dt = np.concatenate((dta,dtb))
    j = np.nonzero(np.isfinite(dt))[0]
    dtmed = np.median(dt[j])
    if options.verbose >= 2:
      print('dtmed=',dtmed)

    dem_e1a_amp = dem_e1_amp[ja]
    dem_e1b_amp = dem_e1_amp[jb]
    dem_e2a_amp = dem_e2_amp[ja]
    dem_e2b_amp = dem_e2_amp[jb]

    # cause pen to lift for gaps larger than for cal
    options.ngap = 5
    if options.ngap > 0:
      jza = np.nonzero(dta > dtmed * options.ngap)[0]
      jzb = np.nonzero(dtb > dtmed * options.ngap)[0]
      ta[jza] = np.nan
      tb[jzb] = np.nan

    e1a_ocean = np.tile(np.nan,len(ta))
    for i in range(len(ta)):
      Ra1 = np.interp(ta[i], tak, R1ap)
      Ra2 = np.interp(ta[i], tbk, R1au)
      Rb1 = np.interp(ta[i], tak, R1bu)
      Rb2 = np.interp(ta[i], tbk, R1bp)
      ead = dem_e1a_amp[i]
      ebd = np.interp(ta[i], tb, dem_e1b_amp)
      try:
        c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
      except:
        c = [np.nan, np.nan]
      e1a_ocean[i] = c[0] * (Ry + Ra1 + Rb1)

    e1b_ocean = np.tile(np.nan,len(tb))
    for i in range(len(tb)):
      Ra1 = np.interp(tb[i], tak, R1ap)
      Ra2 = np.interp(tb[i], tbk, R1au)
      Rb1 = np.interp(tb[i], tak, R1bu)
      Rb2 = np.interp(tb[i], tbk, R1bp)
      ead = np.interp(tb[i], ta, dem_e1a_amp)
      ebd = dem_e1b_amp[i]
      try:
        c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
      except:
        c = [np.nan, np.nan]
      e1b_ocean[i] = c[1] * (Ry + Ra2 + Rb2)

    e2a_ocean = np.tile(np.nan,len(ta))
    for i in range(len(ta)):
      Ra1 = np.interp(ta[i], tak, R2ap)
      Ra2 = np.interp(ta[i], tbk, R2au)
      Rb1 = np.interp(ta[i], tak, R2bu)
      Rb2 = np.interp(ta[i], tbk, R2bp)
      ead = dem_e2a_amp[i]
      ebd = np.interp(ta[i], tb, dem_e2b_amp)
      try:
        c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
      except:
        c = [np.nan, np.nan]
      e2a_ocean[i] = c[0] * (Ry + Ra1 + Rb1)

    e2b_ocean = np.tile(np.nan,len(tb))
    for i in range(len(tb)):
      Ra1 = np.interp(tb[i], tak, R2ap)
      Ra2 = np.interp(tb[i], tbk, R2au)
      Rb1 = np.interp(tb[i], tak, R2bu)
      Rb2 = np.interp(tb[i], tbk, R2bp)
      ead = np.interp(tb[i], ta, dem_e2a_amp)
      ebd = dem_e2b_amp[i]
      try:
        c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
      except:
        c = [np.nan, np.nan]
      e2b_ocean[i] = c[1] * (Ry + Ra2 + Rb2)

    pytfn = datetime(1970,1,1,0,0,0) + timedelta(0,uxtfn)
    datestr = pytfn.strftime('%Y%m%d')
    tlbl = tlbl + ' from ' + pytfn.strftime('%Y-%m-%d 0000 UTC')

    matdict = {                \
      'site':options.site,     \
      'uxt_ref':uxt_ref,       \
      'ta':ta,           \
      'tb':tb,           \
      'e1a_ocean':e1a_ocean,     \
      'e1b_ocean':e1b_ocean,     \
      'e2a_ocean':e2a_ocean,     \
      'e2b_ocean':e2b_ocean,     \
      'dem_e1a_amp':dem_e1a_amp, \
      'dem_e1b_amp':dem_e1b_amp, \
      'dem_e2a_amp':dem_e2a_amp, \
      'dem_e2b_amp':dem_e2b_amp  }
    matname = 'hprwsr-'+datestr+'-'+options.site
    writemat(matodir, matname, matdict)

    if do_plt_ocean:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'hprwsr-'+datestr+'-'+options.site+'-ocean'
      pltnamq = 'hprwsr-ocean'
      if options.verbose:
        print('pltnam=',pltnam)
      fig.suptitle(pltnam+'\nblu:demod, red=ocean, microvolts')

      ylim = [-20,20]

      ax = fig.add_subplot(4,1,1)
      ax.hold(True)
      ax.plot(ta/tdiv,dem_e1a_amp,blumrk)
      ax.plot(ta/tdiv,e1a_ocean,redmrk)
      ax.hold(False)
      ax.set_xlim(xlim)
      ax.set_ylim(ylim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e1a')

      ax = fig.add_subplot(4,1,2)
      ax.hold(True)
      ax.plot(tb/tdiv,dem_e1b_amp,blumrk)
      ax.plot(tb/tdiv,e1b_ocean,redmrk)
      ax.hold(False)
      ax.set_xlim(xlim)
      ax.set_ylim(ylim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e1b')

      ax = fig.add_subplot(4,1,3)
      ax.hold(True)
      ax.plot(ta/tdiv,dem_e2a_amp,blumrk)
      ax.plot(ta/tdiv,e2a_ocean,redmrk)
      ax.hold(False)
      ax.set_xlim(xlim)
      ax.set_ylim(ylim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e2a')

      ax = fig.add_subplot(4,1,4)
      ax.hold(True)
      ax.plot(tb/tdiv,dem_e2b_amp,blumrk)
      ax.plot(tb/tdiv,e2b_ocean,redmrk)
      ax.hold(False)
      ax.set_xlim(xlim)
      ax.set_ylim(ylim)
      ax.grid(True)
  #   ax.xaxis.set_ticklabels([])
      plt.ylabel('e2b')

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')

    if do_plt_R1:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'hprwsr-'+datestr+'-'+options.site+'-R1'
      pltnamq = 'hprwsr-R1'
      if options.verbose:
        print('pltnam=',pltnam)
      fig.suptitle(pltnam+'\nR1, ohms')

      ax = fig.add_subplot(2,2,1)
      ax.plot(tak/tdiv,R1ap,redmrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rplim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('R1ap')
      plt.title('A-pinched')

      ax = fig.add_subplot(2,2,2)
      ax.plot(tbk/tdiv,R1au,blumrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rulim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('R1au')
      plt.title('B-pinched')

      ax = fig.add_subplot(2,2,3)
      ax.plot(tak/tdiv,R1bu,blumrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rulim)
      ax.grid(True)
      plt.ylabel('R1bu')
      plt.xlabel(tlbl)

      ax = fig.add_subplot(2,2,4)
      ax.plot(tbk/tdiv,R1bp,redmrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rplim)
      ax.grid(True)
      plt.ylabel('R1bp')
      plt.xlabel(tlbl)

      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')

    if do_plt_R2:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'hprwsr-'+datestr+'-'+options.site+'-R2'
      pltnamq = 'hprwsr-R2'
      if options.verbose:
        print('pltnam=',pltnam)
      fig.suptitle(pltnam+'\nR2, ohms')

      ax = fig.add_subplot(2,2,1)
      ax.plot(tak/tdiv,R2ap,redmrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rplim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('R2ap')
      plt.title('A-pinched')

      ax = fig.add_subplot(2,2,2)
      ax.plot(tbk/tdiv,R2au,blumrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rulim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('R2au')
      plt.title('B-pinched')

      ax = fig.add_subplot(2,2,3)
      ax.plot(tak/tdiv,R2bu,blumrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rulim)
      ax.grid(True)
      plt.ylabel('R2bu')
      plt.xlabel(tlbl)

      ax = fig.add_subplot(2,2,4)
      ax.plot(tbk/tdiv,R2bp,redmrk)
      ax.set_xlim(xlim)
  #   ax.set_ylim(Rplim)
      ax.grid(True)
      plt.ylabel('R2bp')
      plt.xlabel(tlbl)

      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')


    if do_plt_demod:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'hprwsr-'+datestr+'-'+options.site+'-demod'
      pltnamq = 'hprwsr-demod'
      if options.verbose:
        print('pltnam=',pltnam)
      fig.suptitle(pltnam+'\nDemod, microvolts')

      ylim = [-20,20]

      ax = fig.add_subplot(4,1,1)
      plt.sca(ax)
      ax.plot(ta/tdiv,dem_e1a_amp,blumrk,ms=mrksz_eftt)
      ax.set_ylim(ylim)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e1a')

      ax = fig.add_subplot(4,1,2)
      plt.sca(ax)
      ax.plot(tb/tdiv,dem_e1b_amp,blumrk,ms=mrksz_eftt)
      ax.set_ylim(ylim)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e1b')

      ax = fig.add_subplot(4,1,3)
      plt.sca(ax)
      ax.plot(ta/tdiv,dem_e2a_amp,blumrk,ms=mrksz_eftt)
      ax.set_ylim(ylim)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e2a')

      ax = fig.add_subplot(4,1,4)
      plt.sca(ax)
      ax.plot(tb/tdiv,dem_e2b_amp,blumrk,ms=mrksz_eftt)
      ax.set_ylim(ylim)
      ax.set_xlim(xlim)
      ax.grid(True)
      plt.ylabel('e2b')
      plt.xlabel(tlbl)

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')

    if do_plt_cal:

      for Nint in range(2):
        Nstr = str(Nint+1)

        fig = plt.figure(num=1,figsize=(10, 7))
        fig.clf()
        pltnam = 'hprwsr-'+datestr+'-'+options.site+'-cal'+Nstr
        pltnamq = 'hprwsr-cal'+Nstr
        if options.verbose:
          print('pltnam=',pltnam)
        fig.suptitle(pltnam+'\ne'+Nstr+' cal, microvolts')

        if options.site[0] == 'H':
          ylim_big = [0,600]
          ylim_small = [-1,5]
        else:
          ylim_big = [0,1500]
          ylim_small = [-1,15]
        ylim_ab_ratio = [-0.01,0.10]
        ylim_abc_ratio = [-0.20,0.20]

        ja = np.nonzero(cal_abu == 'a')[0]
        jb = np.nonzero(cal_abu == 'b')[0]

        ta = cal_uxt[ja] - uxt_ref
        tb = cal_uxt[jb] - uxt_ref

        dta = np.diff(ta)
        dtb = np.diff(tb)
        dta_med = np.median(dta)
        dtb_med = np.median(dtb)
        jza = np.nonzero(dta > dta_med * 1.5)[0]
        jzb = np.nonzero(dtb > dtb_med * 1.5)[0]
        ta[jza] = np.nan
        tb[jzb] = np.nan

        if Nstr == '1':
          eNa_amp = cal_e1a_amp
          eNb_amp = cal_e1b_amp
          eNc_amp = cal_e1c_amp
        elif Nstr == '2':
          eNa_amp = cal_e2a_amp
          eNb_amp = cal_e2b_amp
          eNc_amp = cal_e2c_amp
        else:
          print('bad Nstr=',Nstr)
          exit(1)

        ax = fig.add_subplot(5,2,1)
        ax.plot(ta/tdiv,eNa_amp[ja],redmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_big)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'ap')
        plt.title('A pinched')

        ax = fig.add_subplot(5,2,2)
        ax.plot(tb/tdiv,eNa_amp[jb],blumrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_small)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'au')
        plt.title('B pinched')

        ax = fig.add_subplot(5,2,3)
        ax.plot(ta/tdiv,eNb_amp[ja],blumrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_small)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'bu')


        ax = fig.add_subplot(5,2,4)
        ax.plot(tb/tdiv,eNb_amp[jb],redmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_big)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'bp')

        ax = fig.add_subplot(5,2,5)
        ax.plot(ta/tdiv,eNc_amp[ja],grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_big)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'c')

        ax = fig.add_subplot(5,2,6)
        ax.plot(tb/tdiv,eNc_amp[jb],grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_big)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('e'+Nstr+'c')

        ax = fig.add_subplot(5,2,7)
        ax.plot(ta/tdiv,(eNa_amp[ja]+eNb_amp[ja])/eNc_amp[ja]-1,grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_abc_ratio)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('(a+b)/c-1')

        ax = fig.add_subplot(5,2,8)
        ax.plot(tb/tdiv,(eNa_amp[jb]+eNb_amp[jb])/eNc_amp[jb]-1,grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_abc_ratio)
        ax.set_xlim(xlim)
        ax.grid(True)
        ax.xaxis.set_ticklabels([])
        plt.ylabel('(a+b)/c-1')

        ax = fig.add_subplot(5,2,9)
        ax.plot(ta/tdiv,eNb_amp[ja]/eNa_amp[ja],grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_ab_ratio)
        ax.set_xlim(xlim)
        ax.grid(True)
        plt.ylabel('b/a')
        plt.xlabel(tlbl)

        ax = fig.add_subplot(5,2,10)
        ax.plot(tb/tdiv,eNa_amp[jb]/eNb_amp[jb],grnmrk,ms=mrksz)
        if options.do_cal_ylim:
          ax.set_ylim(ylim_ab_ratio)
        ax.set_xlim(xlim)
        ax.grid(True)
        plt.ylabel('a/b')
        plt.xlabel(tlbl)

        fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

        writepdf(pltdir, pltnam)
        os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')

    if do_plt_mot:
      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'hprwsr-'+datestr+'-'+options.site+'-mot'
      pltnamq = 'hprwsr-mot'
      if options.verbose:
        print('pltnam=',pltnam)
      fig.suptitle(pltnam+'\nMotor current & duration, red=A, blu=B')

      ja = np.nonzero(mot_abu == 'a')[0]
      jb = np.nonzero(mot_abu == 'b')[0]

      ta = mot_uxt1[ja] - uxt_ref
      tb = mot_uxt1[jb] - uxt_ref

      dta = np.diff(ta)
      dtb = np.diff(tb)
      dta_med = np.median(dta)
      dtb_med = np.median(dtb)
      jza = np.nonzero(dta > dta_med * 1.5)[0]
      jzb = np.nonzero(dtb > dtb_med * 1.5)[0]
      ta[jza] = np.nan
      tb[jzb] = np.nan

      ax = fig.add_subplot(4,2,1)
      plt.sca(ax)
      plt.title('Motor 1')
      ax.plot(ta/tdiv,mot_cur1[ja],redmrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('cur1a')

      ax = fig.add_subplot(4,2,3)
      plt.sca(ax)
      ax.plot(tb/tdiv,mot_cur1[jb],blumrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('cur1b')

      ax = fig.add_subplot(4,2,5)
      plt.sca(ax)
      ax.plot(ta/tdiv,mot_dur1[ja],redmrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('dur1a')

      ax = fig.add_subplot(4,2,7)
      plt.sca(ax)
      ax.plot(tb/tdiv,mot_dur1[jb],blumrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      plt.ylabel('dur1b')
      plt.xlabel(tlbl)

      ax = fig.add_subplot(4,2,2)
      plt.sca(ax)
      plt.title('Motor 2')
      ax.plot(ta/tdiv,mot_cur2[ja],redmrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('cur2a')

      ax = fig.add_subplot(4,2,4)
      plt.sca(ax)
      ax.plot(tb/tdiv,mot_cur2[jb],blumrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('cur2b')

      ax = fig.add_subplot(4,2,6)
      plt.sca(ax)
      ax.plot(ta/tdiv,mot_dur2[ja],redmrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('dur2a')

      ax = fig.add_subplot(4,2,8)
      plt.sca(ax)
      ax.plot(tb/tdiv,mot_dur2[jb],blumrkmot)
      adj_ylim(ax)
      ax.set_xlim(xlim)
      ax.grid(True)
      plt.ylabel('dur2b')
      plt.xlabel(tlbl)

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf' + ' ' + pltdirq + '/' + pltnamq + '.pdf')

if __name__ == '__main__':
  main()
