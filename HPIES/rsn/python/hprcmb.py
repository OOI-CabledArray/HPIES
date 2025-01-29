#!/usr/bin/env python2
# hprcmb.py -- combine hpies data
#   collect ndays of data
#   estimate switch resistances
#   make ocean data
#   filter results

"""
./hprcmb.py -d 20140822 -s 1AB --ndays=160 --filtdays=1 -m -p
./hprcmb.py -d 20140822 -s 2SB --ndays=160 --filtdays=1 -m -p
./hprcmb.py -d 20120615 -s H1  --ndays=350 --filtdays=1 -m -p
./hprcmb.py -d 20120615 -s H2  --ndays=350 --filtdays=1 -m -p
./hprcmb.py -d 20120615 -s H3  --ndays=350 --filtdays=1 -m -p
./hprcmb.py -d 20120615 -s H4  --ndays=350 --filtdays=1 -m -p
./hprcmb.py -d 20120615 -s H5  --ndays=350 --filtdays=1 -m -p

for i in H1 H2 H3 H4 H5; do ./hprcmb.py -d 20120615 -s $i --ndays=350 --filtdays=1; done 
for i in 1AB 2SB;        do ./hprcmb.py -d 20140822 -s $i --ndays=140 --filtdays=1; done

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

def uxt2str(uxt):
  pyt = datetime(1970,1,1) + timedelta(0,uxt)
  return pyt.strftime('%Y-%m-%d %H:%M:%S')

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

def writepdf(pltdir, namewoext):
  mymkdir(pltdir)
  pdffile = os.path.join(pltdir,namewoext+'.pdf')
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

def blkavg(a,m):
  n = int(len(a) / m)
  return np.mean(np.reshape(a[0:n*m],(m,n),order='F'),0)

def mylim(x,m):
  j = np.nonzero(x >  m)[0]; x[j] =  m
  j = np.nonzero(x < -m)[0]; x[j] = -m
  return x

def main():

  if False:
    # test blkavg()
    a = [1,2,1,2,1,2,1,2,1,2,1,2]
    n = 6
    b = blkavg(a,n)
    print('a=',a)
    print('n=',n)
    print('b=',b)
    exit(0)


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

  parser.add_option("-n","--ndays", dest="ndays",
    type="int", metavar="N", default=1,
    help="number of days to plot")

  parser.add_option("-f","--filtdays", dest="filtdays",
    type="float", metavar="D", default=0.0,
    help="number of days to filter")

  parser.add_option("-g","--ngap", dest="ngap",
    type="int", metavar="G", default=0.0,
    help="insert NaN if gap is G times median of time diffs")

  parser.add_option("-m","--do_mat", dest="do_mat", 
    action="store_true", default=False,
    help="make matlab file")

  parser.add_option("--do_cal_ylim", dest="do_cal_ylim", 
    action="store_true", default=False,
    help="apply manual ylimit to cal plots")

  parser.add_option("-p", "--mkplt", dest="mkplt", 
    action="store_true", default=False,
    help="make plots")

  parser.add_option("--ttpr", dest="plt_tt_pr", 
    action="store_true", default=False,
    help="plot travel times and pressure")

  parser.add_option("--plot_date", dest="do_plot_date", 
    action="store_true", default=False,
    help="plot dates")

  global options
  (options, args) = parser.parse_args()

  if options.mkplt:
    do_plt_R1      = True
    do_plt_R2      = True
    do_plt_cal     = True
    do_plt_ocean   = True
    do_plt_uv      = True
    do_plt_tt_pr   = True
  else:
    do_plt_R1      = False
    do_plt_R2      = False
    do_plt_cal     = False
    do_plt_ocean   = False
    do_plt_uv      = False
    do_plt_tt_pr   = False

  do_plt_demod   = False
  do_plt_std     = False
  do_plt_self    = False
  do_plt_mot     = False

  if options.plt_tt_pr:
    do_plt_tt_pr = True

  filtsecs = options.filtdays * 86400

  add_ndays      = str.format('-{0}d',options.ndays)
  add_ndays_filt = str.format('-{0}d-f{1}d',options.ndays,options.filtdays)
  add_filtdays   = ', filter: {0} days'.format(options.filtdays)


  if options.ndays < 5:
    tlbl = 'hours'
    tdiv = 3600
    xlim = [0,options.ndays*24]
    redmrk = 'r.-'
    blumrk = 'b.-'
    grnmrk = 'g.-'
    redmrkmot = 'r.-'
    blumrkmot = 'b.-'
    grnmrkmot = 'g.-'
    mrktt = 'b.-'
    mrksz = 1
  else:
    tlbl = 'days'
    tdiv = 86400
    xlim = [0,options.ndays]
    redmrk = 'r-'
    blumrk = 'b-'
    grnmrk = 'g-'
    redmrkmot = 'r-'
    blumrkmot = 'b-'
    grnmrkmot = 'g-'
    mrktt = 'b.-'
    mrksz = 1

  if filtsecs == 0:
    redmrk_eftt = redmrk
    blumrk_eftt = blumrk
    grnmrk_eftt = grnmrk
    mrktt_eftt = mrktt
    mrksz_eftt = mrksz
  else:
    redmrk_eftt = 'ro-'
    blumrk_eftt = 'bo-'
    grnmrk_eftt = 'go-'
    mrktt_eftt = 'b.-'
    mrksz_eftt = 7

  if options.date == None or options.site == None:
    parser.print_help()
    exit(1)

  if options.site[0] == 'H':
    deployment = 'OKMC'
    pydir = '/data/okmc/mat2/hprpro' # hprpro.py tskip=40, twait=10
    pydir = '/data/okmc/mat/hprpro'  # hprpro.py tskip=15, twait=7
    pydir = '/data/okmc/mat3/hprpro' # hprpro.py tskip=40, twait=10
    pydir = '/data/okmc/mat4/hprpro'
    pltdir = '/data/okmc/plots/hprcmb'
    matdir = '/data/okmc/mat4/hprcmb'
  else:
    deployment = 'RSN'
    pydir = '/data/rsn/mat/hprpro'  # hprpro.py tskip=15, twait=7
    pltdir = '/data/rsn/plots/hprcmb'
    matdir = '/data/rsn/mat/hprcmb'
    do_plt_demodMP = False
    do_plt_demodAB = False

  pltdirq = './plots'

  try:
    uxtbeg = timegm(strptime(options.date,'%Y%m%d'))
  except:
    print('cannot decode --date=',options.date)
    exit(1)

  dem_uxt     = []
  dem_abm     = []
  e1dem = []
  e2dem = []
  e1std   = []
  e2std   = []
  e1self  = []
  e2self  = []

  cal_uxt   = []
  cal_abu   = []
  cal_e1a_amp  = []
  cal_e1b_amp  = []
  cal_e1c_amp  = []
  cal_e2a_amp  = []
  cal_e2b_amp  = []
  cal_e2c_amp  = []

  mot_abu  = []
  mot_uxt1 = []
  mot_uxt2 = []
  mot_cur1 = []
  mot_cur2 = []
  mot_dur1 = []
  mot_dur2 = []

  aux_uxt = []
  tt1  = []
  tt2  = []
  tt3  = []
  tt4  = []
  pres = []

  n_dem = 0
  n_cal = 0
  n_mot = 0
  n_aux = 0

  for iday in range(options.ndays):
    pyt = datetime(1970,1,1,0,0,0) + timedelta(0,uxtbeg + 86400 * iday)
    pytstr = pyt.strftime('%Y%m%d')
    pyfile = pydir+'/HPIES-'+pytstr+'-'+options.site+'-hprpro.mat'
    if options.verbose:
      print('pyfile=',pyfile)

    try:
      PY = loadmat(pyfile)
    except:
      if options.verbose:
        print('warning: cannot loadmat(',pyfile,')')
      continue

    if len(PY['dem_uxt']) > 0:
      dem_uxt = np.append(dem_uxt, PY['dem_uxt'][0])
      dem_abm = np.append(dem_abm, PY['dem_abm'])
      e1dem   = np.append(e1dem,   PY['dem_e1_amp'][0])
      e2dem   = np.append(e2dem,   PY['dem_e2_amp'][0])
      e1std   = np.append(e1std,   PY['dem_e1_std'][0])
      e2std   = np.append(e2std,   PY['dem_e2_std'][0])
      e1self  = np.append(e1self,  PY['dem_e1_self'][0])
      e2self  = np.append(e2self,  PY['dem_e2_self'][0])
      n_dem += 1

    if len(PY['cal_uxt']) > 0:
      cal_uxt = np.append(cal_uxt,    PY['cal_uxt'][0])
      cal_abu = np.append(cal_abu,    PY['cal_abu'])
      cal_e1a_amp  = np.append(cal_e1a_amp, PY['cal_e1a_amp'][0])
      cal_e1b_amp  = np.append(cal_e1b_amp, PY['cal_e1b_amp'][0])
      cal_e1c_amp  = np.append(cal_e1c_amp, PY['cal_e1c_amp'][0])
      cal_e2a_amp  = np.append(cal_e2a_amp, PY['cal_e2a_amp'][0])
      cal_e2b_amp  = np.append(cal_e2b_amp, PY['cal_e2b_amp'][0])
      cal_e2c_amp  = np.append(cal_e2c_amp, PY['cal_e2c_amp'][0])
      n_cal += 1

    if len(PY['mot_abu']) > 0:
      mot_abu  = np.append(mot_abu,  PY['mot_abu'])
      mot_uxt1 = np.append(mot_uxt1, PY['mot_uxt1'][0])
      mot_uxt2 = np.append(mot_uxt2, PY['mot_uxt2'][0])
      mot_cur1 = np.append(mot_cur1, PY['mot_cur1'][0])
      mot_cur2 = np.append(mot_cur2, PY['mot_cur2'][0])
      mot_dur1 = np.append(mot_dur1, PY['mot_dur1'][0])
      mot_dur2 = np.append(mot_dur2, PY['mot_dur2'][0])
      n_mot += 1

    if len(PY['aux_uxt']) > 0:
      aux_uxt = np.append(aux_uxt, PY['aux_uxt'][0])
      tt1  = np.append(tt1 , PY['aux_tt1'][0])
      tt2  = np.append(tt2 , PY['aux_tt2'][0])
      tt3  = np.append(tt3 , PY['aux_tt3'][0])
      tt4  = np.append(tt4 , PY['aux_tt4'][0])
      pres = np.append(pres, PY['aux_pres'][0])
      n_aux += 1

  # finished accumulating input files
  print('n_dem=',n_dem,'n_cal=',n_cal,'n_mot=',n_mot,'n_aux=',n_aux)

  if len(dem_uxt) == 0:
    print('no dem data found -- exit')
    exit(1)

  uxt_ref = uxtbeg

  if options.verbose > 1:
    print('len(dem_uxt)=',len(dem_uxt))
    print('len(cal_uxt)=',len(cal_uxt))
    print('len(mot_uxt1)=',len(mot_uxt1))
    print('len(mot_uxt2)=',len(mot_uxt2))
    print('len(aux_uxt)=',len(aux_uxt))

  taux = aux_uxt - uxt_ref

  jm1 = np.searchsorted(mot_uxt1, dem_uxt)
  jm2 = np.searchsorted(mot_uxt2, dem_uxt)

  for i in range(1,len(dem_uxt)):
    j1=jm1[i]
    j2=jm2[i]
    if j1 == len(mot_uxt1):
      j1 -= 1
    if j2 == len(mot_uxt2):
      j2 -= 1
    if j1-2 < 0 or j1 >= len(mot_uxt1):
      print('mot_uxt1 j1 probs: i=',i,'j1=',j1,'len(mot_uxt1)=',len(mot_uxt1),'skipped')
      continue
    if j2-2 < 0 or j2 >= len(mot_uxt2):
      print('mot_uxt2 j2 probs: i=',i,'j2=',j2,'len(mot_uxt2)=',len(mot_uxt2),'skipped')
      continue
    te  = dem_uxt[i] - uxt_ref
    tm11 = mot_uxt1[j1-2] - uxt_ref - te
    tm12 = mot_uxt1[j1-1] - uxt_ref - te
    tm13 = mot_uxt1[j1]   - uxt_ref - te
    tm21 = mot_uxt2[j2-2] - uxt_ref - te
    tm22 = mot_uxt2[j2-1] - uxt_ref - te
    tm23 = mot_uxt2[j2]   - uxt_ref - te

    do_zap_tmo = False
    if do_zap_tmo:
      if options.site == 'H5':
        tmo = 1.5
      else:
        tmo = 4.5
      if mot_dur1[j1-2] > tmo or mot_dur1[j1-1] > tmo or mot_dur1[j1] > tmo:
        if False:
          print('i=',i,'j1=',j1,'zap e1','ab=',dem_abm[i],end=' ')
          print('dur1=',mot_dur1[j1-2], mot_dur1[j1-1],mot_dur1[j1])
        e1dem[i] = np.nan
      if mot_dur2[j2-2] > tmo or mot_dur2[j2-1] > tmo or mot_dur2[j2] > tmo:
        if False:
          print('i=',i,'j2=',j2,'zap e2','ab=',dem_abm[i],end=' ')
          print('dur2=',mot_dur2[j2-2], mot_dur2[j2-1],mot_dur2[j2])
        e2dem[i] = np.nan

    ok = True
    if tm11 > 0 or tm12 > 0 or tm13 < 0:
      print('mot_uxt1 tm probs, i=',i,'dem_uxt[i]=',uxt2str(dem_uxt[i]))
      if options.verbose:
        print('  i={0:3d}'.format(i))
        print('  tm11={0:7.0f}'.format(tm11))
        print('  tm12={0:7.0f}'.format(tm12))
        print('  tm13={0:7.0f}'.format(tm13))
        print('  tm12-tm11={0:7.0f}'.format(tm12-tm11))
        print('  tm13-tm12={0:7.0f}'.format(tm13-tm12))
        print('  te=',te,'dem_uxt[i]=',dem_uxt[i],'uxt_ref=',uxt_ref)
    if tm21 > 0 or tm22 > 0 or tm23 < 0:
      print('mot_uxt2 tm probs, i=',i,'dem_uxt[i]=',uxt2str(dem_uxt[i]))
      if options.verbose:
        print('  i={0:3d}'.format(i))
        print('  tm21={0:7.0f}'.format(tm21))
        print('  tm22={0:7.0f}'.format(tm22))
        print('  tm24={0:7.0f}'.format(tm23))
        print('  tm22-tm21={0:7.0f}'.format(tm22-tm21))
        print('  tm23-tm22={0:7.0f}'.format(tm23-tm22))
        print('  te=',te,'dem_uxt[i]=',dem_uxt[i],'uxt_ref=',uxt_ref)

    if options.verbose > 1:
      print('len(jm)=',len(jm))
      print('len(dem_uxt)=',len(dem_uxt))
      print('len(mot_uxt1)=',len(mot_uxt1))
      print('len(mot_uxt2)=',len(mot_uxt2))

  ja = np.nonzero(dem_abm == 'a')[0]
  jb = np.nonzero(dem_abm == 'b')[0]
  if options.verbose > 1:
    print('len(dem_abm)=',len(dem_abm))
    print('len(ja)=',len(ja))
    print('len(jb)=',len(jb))
    if len(ja) != len(jb):
      print('len(ja) != len(jb)')
    print('ja[0] =',ja[0], 'ja[-1]=',ja[-1])
    print('jb[0] =',jb[0], 'jb[-1]=',jb[-1])
    print('dem_uxt[ja[0]] =',dem_uxt[ja[0]], 'dem_uxt[ja[-1]]=',dem_uxt[ja[-1]])
    print('dem_uxt[jb[0]] =',dem_uxt[jb[0]], 'dem_uxt[jb[-1]]=',dem_uxt[jb[-1]])

  jab = sorted(set(ja).intersection(set(jb+1)))
  jba = sorted(set(jb).intersection(set(ja+1)))
  if options.verbose:
    print('len(jab)=',len(jab))
    print('len(jba)=',len(jba))
  if len(jab) >= len(jba):
    jan = np.array(jab)
    jbn = np.array(jab)-1
  if len(jba) > len(jab):
    jan = np.array(jba)-1
    jbn = np.array(jba)
  jat = set(ja).difference(jan)
  jbt = set(jb).difference(jbn)
  if len(jat)>0 or len(jbt)>0:
    print('keeping only adjacent a&b: tossed: a:',len(jat),'b:',len(jbt))
  ja = jan
  jb = jbn

  # apply factor of six to e1b and e1c for H3 -- JHD
  if options.site == 'H3':
    e1bsf = 6.0
    print('H3 e1b,e1c scale factor=',e1bsf)
    e1dem[jb]  *= e1bsf
    e1std[jb]  *= e1bsf
    e1self[jb] *= e1bsf
    cal_e1b_amp *= e1bsf
    cal_e1c_amp *= e1bsf

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
  ta_cal = cal_uxt[jak] - uxt_ref
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
  tb_cal = cal_uxt[jbk] - uxt_ref
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
  dem_ta = dem_uxt[ja] - uxt_ref
  dem_tb = dem_uxt[jb] - uxt_ref
  if options.verbose > 1:
    print('uxt_ref=',uxt_ref)

  dta = np.diff(dem_ta)
  dtb = np.diff(dem_tb)
  dt = np.concatenate((dta,dtb))
  j = np.nonzero(np.isfinite(dt))[0]
  dt_med = np.median(dt[j])
  if options.verbose > 1:
    print('dt_med=',dt_med)

  e1a_dem = e1dem[ja]
  e1b_dem = e1dem[jb]
  e2a_dem = e2dem[ja]
  e2b_dem = e2dem[jb]

  e1a_std = e1std[ja]
  e1b_std = e1std[jb]
  e2a_std = e2std[ja]
  e2b_std = e2std[jb]

  e1a_self = e1self[ja]
  e1b_self = e1self[jb]
  e2a_self = e2self[ja]
  e2b_self = e2self[jb]

  # cause pen to lift for gaps larger than for cal
  if options.ngap > 0:
    jza = np.nonzero(dta > dt_med * options.ngap)[0]
    jzb = np.nonzero(dtb > dt_med * options.ngap)[0]
    dem_ta[jza] = np.nan
    dem_tb[jzb] = np.nan

  # filtering
  if filtsecs > 0:
    nfilt_ef = int(filtsecs / dt_med)
    if options.verbose > 1:
      print('nfilt_ef=',nfilt_ef)

    dem_ta = blkavg(dem_ta,nfilt_ef)
    dem_tb = blkavg(dem_tb,nfilt_ef)

    e1a_dem = blkavg(e1a_dem,nfilt_ef)
    e1b_dem = blkavg(e1b_dem,nfilt_ef)
    e2a_dem = blkavg(e2a_dem,nfilt_ef)
    e2b_dem = blkavg(e2b_dem,nfilt_ef)

    e1a_std = blkavg(e1a_std,nfilt_ef)
    e1b_std = blkavg(e1b_std,nfilt_ef)
    e2a_std = blkavg(e2a_std,nfilt_ef)
    e2b_std = blkavg(e2b_std,nfilt_ef)

    e1a_self = blkavg(e1a_self,nfilt_ef)
    e1b_self = blkavg(e1b_self,nfilt_ef)
    e2a_self = blkavg(e2a_self,nfilt_ef)
    e2b_self = blkavg(e2b_self,nfilt_ef)

    if False:
      dtaux_med = np.median(np.diff(taux))
      nfilt_aux = int(filtsecs / dtaux_med)
      if options.verbose > 1:
        print('nfilt_aux=',nfilt_aux)
      taux = blkavg(taux,nfilt_aux)
      tt1  = blkavg(tt1, nfilt_aux)
      tt2  = blkavg(tt2, nfilt_aux)
      tt3  = blkavg(tt3, nfilt_aux)
      tt4  = blkavg(tt4, nfilt_aux)
      pres = blkavg(pres,nfilt_aux)

  e1a_ocean = np.tile(np.nan,len(dem_ta))
  for i in range(len(dem_ta)):
    Ra1 = np.interp(dem_ta[i], ta_cal, R1ap)
    Ra2 = np.interp(dem_ta[i], tb_cal, R1au)
    Rb1 = np.interp(dem_ta[i], ta_cal, R1bu)
    Rb2 = np.interp(dem_ta[i], tb_cal, R1bp)
    ead = e1a_dem[i]
    ebd = np.interp(dem_ta[i], dem_tb, e1b_dem)
    try:
      c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
    except:
      print('note: set e1a_ocean arm currents to NaN, i=',i)
      c = [np.nan, np.nan]
    e1a_ocean[i] = c[0] * (Ry + Ra1 + Rb1)

  e1b_ocean = np.tile(np.nan,len(dem_tb))
  for i in range(len(dem_tb)):
    Ra1 = np.interp(dem_tb[i], ta_cal, R1ap)
    Ra2 = np.interp(dem_tb[i], tb_cal, R1au)
    Rb1 = np.interp(dem_tb[i], ta_cal, R1bu)
    Rb2 = np.interp(dem_tb[i], tb_cal, R1bp)
    ead = np.interp(dem_tb[i], dem_ta, e1a_dem)
    ebd = e1b_dem[i]
    try:
      c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
    except:
      print('note: set e1b_ocean arm currents to NaN, i=',i)
      c = [np.nan, np.nan]
    e1b_ocean[i] = c[1] * (Ry + Ra2 + Rb2)

  e2a_ocean = np.tile(np.nan,len(dem_ta))
  for i in range(len(dem_ta)):
    Ra1 = np.interp(dem_ta[i], ta_cal, R2ap)
    Ra2 = np.interp(dem_ta[i], tb_cal, R2au)
    Rb1 = np.interp(dem_ta[i], ta_cal, R2bu)
    Rb2 = np.interp(dem_ta[i], tb_cal, R2bp)
    ead = e2a_dem[i]
    ebd = np.interp(dem_ta[i], dem_tb, e2b_dem)
    try:
      c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
    except:
      print('note: set e2a_ocean arm currents to NaN, i=',i)
      c = [np.nan, np.nan]
    e2a_ocean[i] = c[0] * (Ry + Ra1 + Rb1)

  e2b_ocean = np.tile(np.nan,len(dem_tb))
  for i in range(len(dem_tb)):
    Ra1 = np.interp(dem_tb[i], ta_cal, R2ap)
    Ra2 = np.interp(dem_tb[i], tb_cal, R2au)
    Rb1 = np.interp(dem_tb[i], ta_cal, R2bu)
    Rb2 = np.interp(dem_tb[i], tb_cal, R2bp)
    ead = np.interp(dem_tb[i], dem_ta, e2a_dem)
    ebd = e2b_dem[i]
    try:
      c = np.linalg.solve([[Ra1, -Ra2],[Rb1, -Rb2]], [ead, -ebd])
    except:
      print('note: set e2b_ocean arm currents to NaN, i=',i)
      c = [np.nan, np.nan]
    e2b_ocean[i] = c[1] * (Ry + Ra2 + Rb2)

  if deployment == 'OKMC':
    esep = 2.25   # m
    Fz = -17700   # nT
    magvar = -2.6 # degrees
  elif deployment == 'RSN':
    esep = 2.55
  else:
    print('unknown deployment=',deployment)
    exit(1)

  # convert electric potential difference to electric field
  ef1a = -e1a_ocean / esep  
  ef1b = -e1b_ocean / esep 
  ef2a = -e2a_ocean / esep 
  ef2b = -e2b_ocean / esep 

  # combine A and B preamps
  if options.site == 'H1':
    ef1 = (ef1a + ef1b) * 0.5  # old ws25, sign ok
    ef2 = (ef2a + ef2b) * 0.5  # old ws26, sign ok
    dem_tab = (dem_ta + dem_tb) * 0.5
    hdg = 74.2
  elif options.site == 'H2':
    ef1 = -(ef1a + ef1b) * 0.5 # new ws11, flip sign 
    ef2 = -(ef2a + ef2b) * 0.5 # new ws13, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    hdg = 86.5
  elif options.site == 'H3':
    ef1 = ef1a            # old ws23, sign ok
    ef2 = ef2b            # old ws27, sign ok
    dem_tab = (dem_ta + dem_tb) * 0.5
    hdg = 182.8
    j = np.nonzero(dem_tab < date2uxt(2012,12,03,13,0,0) - uxt_ref)[0]
    ef2[j] = np.nan
  elif options.site == 'H4':
    ef1 = -(ef1a + ef1b) * 0.5  # new ws10, flip sign
    ef2 =  (ef2a + ef2b) * 0.5  # old ws28, sign ok
    dem_tab = (dem_ta + dem_tb) * 0.5
    hdg = 58.4
    j = np.nonzero(dem_tab < date2uxt(2012,12,01,0,0,0) - uxt_ref)[0]
    ef2[j] = np.nan
  elif options.site == 'H5':
    ef1 =  (ef1a + ef1b) * 0.5   # old ws24, sign ok
    ef2 = -(ef2a + ef2b) * 0.5   # new ws12, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    hdg = 49.7
  elif options.site == '1AB':
    ef1 = -(ef1a + ef1b) * 0.5   # new ws, flip sign
    ef2 = -(ef2a + ef2b) * 0.5   # new ws, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    magvar = 16.5
    Fz = -47263
    hdg = 78.6
  elif options.site == 'AB':
    ef1 = -(ef1a + ef1b) * 0.5   # new ws, flip sign
    ef2 = -(ef2a + ef2b) * 0.5   # new ws, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    magvar = 16.5
    Fz = -47263
    hdg = np.NaN
  elif options.site == '2SB':
    ef1 = -(ef1a + ef1b) * 0.5   # new ws, flip sign
    ef2 = -(ef2a + ef2b) * 0.5   # new ws, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    magvar = 15.8
    Fz = -47334
    hdg = 166
  elif options.site == 'SB':
    ef1 = -(ef1a + ef1b) * 0.5   # new ws, flip sign
    ef2 = -(ef2a + ef2b) * 0.5   # new ws, flip sign
    dem_tab = (dem_ta + dem_tb) * 0.5
    magvar = 15.8
    Fz = -47334
    hdg = np.NaN
  else:
    print('unknown options.site')
    exit(1)

  # rotate electric field to east and north components
  ang = (hdg + magvar - 180) * math.pi / 180.0
  c = math.cos(ang)
  s = math.sin(ang)
  efx = c * ef1 + s * ef2
  efy = c * ef2 - s * ef1

  u = efy * 1e3 / Fz
  v = efx * 1e3 / Fz

  if options.site == 'H1':
    ttref = 1.850
    prref = 1401.5
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'H2':
    ttref = 1.795
    prref = 1360.5
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'H3':
    prref = 2110
    ttref = 2.798
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'H4':
    prref = 4658
    ttref = 6.090
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'H5':
    prref = 2186
    ttref = 2.890
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == '1AB':
    ttref = 3.518
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prref = 2658.2
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'AB':
    ttref = 3.518
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prref = 2658.2
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == '2SB':
    ttref = 3.905
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    prref = 2956.0
    ylim_pres = [-0.5, +0.5]
    prlo = prref - 5.0
    prhi = prref + 5.0
  elif options.site == 'SB':
    ttref = 3.904
    ttlo = ttref - 0.03
    tthi = ttref + 0.03
    ylim_tt = [-0.003,0.003]
    prref = 2954.8
    ylim_pres = [-0.2, +0.2]
    prlo = prref - 5.0
    prhi = prref + 5.0
  else:
    print('unknown site=',options.site,'when filtering tt and pr')
    exit(1)

  # make pressure and travel time averages align with EF averages
  if options.filtdays > 0:
    if options.verbose:
      print('filtering by',options.filtdays,'days')
    ftaux = np.tile(np.nan,len(dem_tab))
    ftt1  = np.tile(np.nan,len(dem_tab))
    ftt2  = np.tile(np.nan,len(dem_tab))
    ftt3  = np.tile(np.nan,len(dem_tab))
    ftt4  = np.tile(np.nan,len(dem_tab))
    fpres = np.tile(np.nan,len(dem_tab))
    for i in range(len(dem_tab)):
      tbeg = dem_tab[i] - options.filtdays * 86400 / 2.0
      tend = dem_tab[i] + options.filtdays * 86400 / 2.0
      jt = np.nonzero((taux >= tbeg) & (taux < tend))[0]
      j1 = np.nonzero((tt1[jt] >= ttlo) & (tt1[jt] <= tthi))[0]
      j2 = np.nonzero((tt2[jt] >= ttlo) & (tt2[jt] <= tthi))[0]
      j3 = np.nonzero((tt3[jt] >= ttlo) & (tt3[jt] <= tthi))[0]
      j4 = np.nonzero((tt4[jt] >= ttlo) & (tt4[jt] <= tthi))[0]
      jp = np.nonzero((pres[jt] >= prlo) & (pres[jt] <= prhi))[0]

#     print('len(jt)=',len(jt))
#     print('len(pres[jt])=',len(pres[jt]))
#     print('len(jp)=',len(jp))
#     print('len(jt[jp])=',len(jt[jp]))
#     print('len(pres[jt[jp]])=',len(pres[jt[jp]]))
      
#     ok = True
#     if len(jt) == 0:
#       print('len(jt)==0')
#       ok = False
#     if len(taux)==0:
#       print('len(taux)==0')
#       ok = False
#     if len(j1)==0:
#       print('len(j1)==0')
#       ok = False
#     if ok == False:
#       print('len(taux)=',len(taux))
#       print('taux[0]=',taux[0])
#       print('taux[-1]=',taux[-1])
#       print('tbeg=',tbeg)
#       print('tend=',tend)
      if len(jt) > 0:
        ftaux[i] = np.mean(taux[jt])
        if len(j1) > 0:
          ftt1[i]  = np.mean(tt1[jt[j1]])
        if len(j2) > 0:
          ftt2[i]  = np.mean(tt2[jt[j2]])
        if len(j3) > 0:
          ftt3[i]  = np.mean(tt3[jt[j3]])
        if len(j4) > 0:
          ftt4[i]  = np.mean(tt4[jt[j4]])
        if len(jp) > 0:
          fpres[i] = np.mean(pres[jt[jp]])
    taux = ftaux
    tt1  = ftt1
    tt2  = ftt2
    tt3  = ftt3
    tt4  = ftt4
    pres = fpres
    if options.verbose:
      j = np.nonzero(np.isfinite(pres))[0]
      print('np.mean(pres)=',np.mean(pres[j]))
      print('np.median(pres)=',np.median(pres[j]))
      print('filtering finished')

  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,uxtbeg)
  tlbl = tlbl + ' from ' + pytbeg.strftime('%Y-%m-%d 0000 UTC')

  if options.do_mat:
    matdict = {
      'site':options.site,     
      'nfilt_ef':nfilt_ef,
#     'nfilt_aux':nfilt_aux,
      'ngap':options.ngap,
      'filtdays':options.filtdays,
      'dt_med':dt_med,
#     'dtaux_med':dtaux_med,
      'uxt_ref':uxt_ref,       
      'dem_ta':dem_ta, 'dem_tb':dem_tb,           
      'e1a_dem':e1a_dem, 'e1b_dem':e1b_dem, 'e2a_dem':e2a_dem, 'e2b_dem':e2b_dem, 
      'e1a_ocean':e1a_ocean,  'e1b_ocean':e1b_ocean,  'e2a_ocean':e2a_ocean,  'e2b_ocean':e2b_ocean,     
      'dem_tab':dem_tab,'ef1':ef1, 'ef2':ef2, 'efx':efx, 'efy':efy, 'u':u, 'v':v,
      't_units':'seconds after uxt_ref', 'e_units':'uV', 'ef_units':'uV/m', 'uv_units':'m/s', 
      'esep':esep, 'Fz':Fz, 'hdg':hdg,
      'ta_cal':ta_cal, 'tb_cal':tb_cal, 
      'R1ap':R1ap, 'R1au':R1au, 'R1bp':R1bp, 'R1bu':R1bu, 
      'R2ap':R2ap, 'R2au':R2au, 'R2bp':R2bp, 'R2bu':R2bu, 
      'taux':taux,'tt1':tt1, 'tt2':tt2, 'tt3':tt3, 'tt4':tt4,  'pres':pres,
    }
    matname = 'hprcmb-'+options.date+'-'+options.site+add_ndays_filt
    writemat(matdir, matname, matdict)

  if do_plt_R1:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-R1'+add_ndays_filt
    pltnamq = 'HPIES-R1'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nR1, ohms' + add_filtdays)

    ax = fig.add_subplot(2,2,1)
    ax.plot(ta_cal/tdiv,R1ap,redmrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rplim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('R1ap')
    plt.title('A-pinched')

    ax = fig.add_subplot(2,2,2)
    ax.plot(tb_cal/tdiv,R1au,blumrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rulim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('R1au')
    plt.title('B-pinched')

    ax = fig.add_subplot(2,2,3)
    ax.plot(ta_cal/tdiv,R1bu,blumrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rulim)
    ax.grid(True)
    plt.ylabel('R1bu')
    plt.xlabel(tlbl)

    ax = fig.add_subplot(2,2,4)
    ax.plot(tb_cal/tdiv,R1bp,redmrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rplim)
    ax.grid(True)
    plt.ylabel('R1bp')
    plt.xlabel(tlbl)

    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_R2:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-R2'+add_ndays_filt
    pltnamq = 'HPIES-R2'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nR2, ohms' + add_filtdays)

    ax = fig.add_subplot(2,2,1)
    ax.plot(ta_cal/tdiv,R2ap,redmrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rplim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('R2ap')
    plt.title('A-pinched')

    ax = fig.add_subplot(2,2,2)
    ax.plot(tb_cal/tdiv,R2au,blumrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rulim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('R2au')
    plt.title('B-pinched')

    ax = fig.add_subplot(2,2,3)
    ax.plot(ta_cal/tdiv,R2bu,blumrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rulim)
    ax.grid(True)
    plt.ylabel('R2bu')
    plt.xlabel(tlbl)

    ax = fig.add_subplot(2,2,4)
    ax.plot(tb_cal/tdiv,R2bp,redmrk)
    ax.set_xlim(xlim)
#   ax.set_ylim(Rplim)
    ax.grid(True)
    plt.ylabel('R2bp')
    plt.xlabel(tlbl)

    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')


  if do_plt_demod:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-demod'+add_ndays_filt
    pltnamq = 'HPIES-demod'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nDemod, microvolts' + add_filtdays)

    ylim = [-20,20]

    ax = fig.add_subplot(4,1,1)
    plt.sca(ax)
    ax.plot(ta/tdiv,e1a_dem,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a')

    ax = fig.add_subplot(4,1,2)
    plt.sca(ax)
    ax.plot(tb/tdiv,e1b_dem,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b')

    ax = fig.add_subplot(4,1,3)
    plt.sca(ax)
    ax.plot(ta/tdiv,e2a_dem,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a')

    ax = fig.add_subplot(4,1,4)
    plt.sca(ax)
    ax.plot(tb/tdiv,e2b_dem,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    plt.ylabel('e2b')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_std:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-std'+add_ndays_filt
    pltnamq = 'HPIES-std'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nStandard Deviation, microvolts' + add_filtdays)


    ylim = [0,10]

    ax = fig.add_subplot(4,1,1)
    plt.sca(ax)
    ax.plot(ta/tdiv,e1a_std,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a')

    ax = fig.add_subplot(4,1,2)
    plt.sca(ax)
    ax.plot(tb/tdiv,e1b_std,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b')

    ax = fig.add_subplot(4,1,3)
    plt.sca(ax)
    ax.plot(ta/tdiv,e2a_std,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a')

    ax = fig.add_subplot(4,1,4)
    plt.sca(ax)
    ax.plot(tb/tdiv,e2b_std,blumrk,ms=mrksz_eftt)
    ax.set_ylim(ylim)
    ax.set_xlim(xlim)
    ax.grid(True)
    plt.ylabel('e2b')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_self:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-self'+add_ndays_filt
    pltnamq = 'HPIES-self'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nSelf Potential, microvolts' + add_filtdays)

    ax = fig.add_subplot(4,1,1)
    plt.sca(ax)
    ax.plot(ta/tdiv,e1a_self,blumrk,ms=mrksz_eftt)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a')

    ax = fig.add_subplot(4,1,2)
    plt.sca(ax)
    ax.plot(tb/tdiv,e1b_self,blumrk,ms=mrksz_eftt)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b')

    ax = fig.add_subplot(4,1,3)
    plt.sca(ax)
    ax.plot(ta/tdiv,e2a_self,blumrk,ms=mrksz_eftt)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a')

    ax = fig.add_subplot(4,1,4)
    plt.sca(ax)
    ax.plot(tb/tdiv,e2b_self,blumrk,ms=mrksz_eftt)
    ax.set_xlim(xlim)
    ax.grid(True)
    plt.ylabel('e2b')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_cal:

    for Nint in range(2):
      Nstr = str(Nint+1)

      fig = plt.figure(num=1,figsize=(10, 7))
      fig.clf()
      pltnam = 'HPIES-'+options.date+'-'+options.site+'-cal'+Nstr+add_ndays
      pltnamq = 'HPIES-cal'+Nstr
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

      cal_ta = cal_uxt[ja] - uxt_ref
      cal_tb = cal_uxt[jb] - uxt_ref

      dta = np.diff(cal_ta)
      dtb = np.diff(cal_tb)
      dta_med = np.median(dta)
      dtb_med = np.median(dtb)
      jza = np.nonzero(dta > dta_med * 1.5)[0]
      jzb = np.nonzero(dtb > dtb_med * 1.5)[0]
      cal_ta[jza] = np.nan
      cal_tb[jzb] = np.nan

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
      ax.plot(cal_ta/tdiv,eNa_amp[ja],redmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_big)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'ap')
      plt.title('A pinched')

      ax = fig.add_subplot(5,2,2)
      ax.plot(cal_tb/tdiv,eNa_amp[jb],blumrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_small)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'au')
      plt.title('B pinched')

      ax = fig.add_subplot(5,2,3)
      ax.plot(cal_ta/tdiv,eNb_amp[ja],blumrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_small)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'bu')


      ax = fig.add_subplot(5,2,4)
      ax.plot(cal_tb/tdiv,eNb_amp[jb],redmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_big)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'bp')

      ax = fig.add_subplot(5,2,5)
      ax.plot(cal_ta/tdiv,eNc_amp[ja],grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_big)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'c')

      ax = fig.add_subplot(5,2,6)
      ax.plot(cal_tb/tdiv,eNc_amp[jb],grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_big)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('e'+Nstr+'c')

      ax = fig.add_subplot(5,2,7)
      ax.plot(cal_ta/tdiv,(eNa_amp[ja]+eNb_amp[ja])/eNc_amp[ja]-1,grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_abc_ratio)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('(a+b)/c-1')

      ax = fig.add_subplot(5,2,8)
      ax.plot(cal_tb/tdiv,(eNa_amp[jb]+eNb_amp[jb])/eNc_amp[jb]-1,grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_abc_ratio)
      ax.set_xlim(xlim)
      ax.grid(True)
      ax.xaxis.set_ticklabels([])
      plt.ylabel('(a+b)/c-1')

      ax = fig.add_subplot(5,2,9)
      ax.plot(cal_ta/tdiv,eNb_amp[ja]/eNa_amp[ja],grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_ab_ratio)
      ax.set_xlim(xlim)
      ax.grid(True)
      plt.ylabel('b/a')
      plt.xlabel(tlbl)

      ax = fig.add_subplot(5,2,10)
      ax.plot(cal_tb/tdiv,eNa_amp[jb]/eNb_amp[jb],grnmrk,ms=mrksz)
      if options.do_cal_ylim:
        ax.set_ylim(ylim_ab_ratio)
      ax.set_xlim(xlim)
      ax.grid(True)
      plt.ylabel('a/b')
      plt.xlabel(tlbl)

      fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

      writepdf(pltdir, pltnam)
      os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_mot:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-mot'+add_ndays
    pltnamq = 'HPIES-mot'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nMotor current & duration, red=A, blu=B')

    ja = np.nonzero(mot_abu == 'a')[0]
    jb = np.nonzero(mot_abu == 'b')[0]

    mot_ta = mot_uxt1[ja] - uxt_ref
    mot_tb = mot_uxt1[jb] - uxt_ref

    dta = np.diff(mot_ta)
    dtb = np.diff(mot_tb)
    dta_med = np.median(dta)
    dtb_med = np.median(dtb)
    jza = np.nonzero(dta > dta_med * 1.5)[0]
    jzb = np.nonzero(dtb > dtb_med * 1.5)[0]
    mot_ta[jza] = np.nan
    mot_tb[jzb] = np.nan

    ax = fig.add_subplot(4,2,1)
    plt.sca(ax)
    plt.title('Motor 1')
    ax.plot(mot_ta/tdiv,mot_cur1[ja],redmrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('cur1a')

    ax = fig.add_subplot(4,2,3)
    plt.sca(ax)
    ax.plot(mot_tb/tdiv,mot_cur1[jb],blumrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('cur1b')

    ax = fig.add_subplot(4,2,5)
    plt.sca(ax)
    ax.plot(mot_ta/tdiv,mot_dur1[ja],redmrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('dur1a')

    ax = fig.add_subplot(4,2,7)
    plt.sca(ax)
    ax.plot(mot_tb/tdiv,mot_dur1[jb],blumrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    plt.ylabel('dur1b')
    plt.xlabel(tlbl)

    ax = fig.add_subplot(4,2,2)
    plt.sca(ax)
    plt.title('Motor 2')
    ax.plot(mot_ta/tdiv,mot_cur2[ja],redmrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('cur2a')

    ax = fig.add_subplot(4,2,4)
    plt.sca(ax)
    ax.plot(mot_tb/tdiv,mot_cur2[jb],blumrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('cur2b')

    ax = fig.add_subplot(4,2,6)
    plt.sca(ax)
    ax.plot(mot_ta/tdiv,mot_dur2[ja],redmrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('dur2a')

    ax = fig.add_subplot(4,2,8)
    plt.sca(ax)
    ax.plot(mot_tb/tdiv,mot_dur2[jb],blumrkmot)
    adj_ylim(ax)
    ax.set_xlim(xlim)
    ax.grid(True)
    plt.ylabel('dur2b')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_tt_pr:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-tt-pr'+add_ndays_filt
    pltnamq = 'HPIES-tt-pr'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam + add_filtdays)
#     ', ttref={0}'.format(ttref) + ', prref={0}'.format(prref))

    # ttmed = np.round(np.median(tt1)/0.01) * 0.01
    # ylim_tt = [ttref-0.005,ttref+0.005]
    # ylim_tt = [-0.005,0.005]
    
    # prmed = np.round(np.median(pres)/1.0) * 1.0
    # ylim_pres = [prref-0.5,prref+0.5]
    # ylim_pres = [-0.5,0.5]

    mltplt = (taux + uxt_ref) / 86400 + 719529 - 366

    ax = fig.add_subplot(2,1,1)
    ax1 = ax
    plt.sca(ax)
    if options.do_plot_date:
      ax.plot_date(mltplt,(tt1-ttref)*1000.0,'b.-',ms=mrksz_eftt)
      ax.plot_date(mltplt,(tt2-ttref)*1000.0,'r.-',ms=mrksz_eftt)
      ax.plot_date(mltplt,(tt3-ttref)*1000.0,'g.-',ms=mrksz_eftt)
      ax.plot_date(mltplt,(tt4-ttref)*1000.0,'k.-',ms=mrksz_eftt)
    else:
      ax.plot(taux/tdiv,(tt1-ttref)*1000.0,'b.-',ms=mrksz_eftt)
      ax.plot(taux/tdiv,(tt2-ttref)*1000.0,'r.-',ms=mrksz_eftt)
      ax.plot(taux/tdiv,(tt3-ttref)*1000.0,'g.-',ms=mrksz_eftt)
      ax.plot(taux/tdiv,(tt4-ttref)*1000.0,'k.-',ms=mrksz_eftt)
      ax.set_xlim(xlim)
#   ax.set_ylim(ylim_tt)
    ax.grid(True)
#   ax.xaxis.set_ticklabels([])
    ax.set_title('Travel Times 1,2,3,4'+' minus {0} ms'.format(ttref*1000.0))
    ax.set_ylabel('ms')

#   ax = fig.add_subplot(5,1,2)
#   plt.sca(ax)
#   ax.plot(taux/tdiv,tt2-ttref,mrktt_eftt,ms=mrksz_eftt)
#   ax.set_ylim(ylim_tt)
#   ax.set_xlim(xlim)
#   ax.grid(True)
#   ax.xaxis.set_ticklabels([])
#   plt.ylabel('tt2'+'-{0}'.format(ttref))

#   ax = fig.add_subplot(5,1,3)
#   plt.sca(ax)
#   ax.plot(taux/tdiv,tt3-ttref,mrktt_eftt,ms=mrksz_eftt)
#   ax.set_ylim(ylim_tt)
#   ax.set_xlim(xlim)
#   ax.grid(True)
#   ax.xaxis.set_ticklabels([])
#   plt.ylabel('tt3'+'-{0}'.format(ttref))

#   ax = fig.add_subplot(5,1,4)
#   plt.sca(ax)
#   ax.plot(taux/tdiv,tt4-ttref,mrktt_eftt,ms=mrksz_eftt)
#   ax.set_ylim(ylim_tt)
#   ax.set_xlim(xlim)
#   ax.grid(True)
#   ax.xaxis.set_ticklabels([])
#   plt.ylabel('tt4'+'-{0}'.format(ttref))

    ax = fig.add_subplot(2,1,2, sharex=ax1)
    plt.sca(ax)
    if options.do_plot_date:
      ax.plot_date(mltplt,pres-prref,mrktt_eftt,ms=mrksz_eftt)
    else:
      ax.plot(taux/tdiv,pres-prref,mrktt_eftt,ms=mrksz_eftt)
  #   ax.set_ylim(ylim_pres)
      ax.set_xlim(xlim)
    ax.grid(True)
    ax.set_title('Pressure'+' minus {0} dbar'.format(prref))
    ax.set_ylabel('dbar')

    if options.do_plot_date:
      fig.autofmt_xdate()
    else:
      plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')
   
  if do_plt_ocean:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-ocean'+add_ndays_filt
    pltnamq = 'HPIES-ocean'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam+'\nblu:demod, red=ocean, microvolts' + add_filtdays)

    ylim = [-20,20]

    ax = fig.add_subplot(4,1,1)
    ax.hold(True)
    print('len(dem_ta)=',len(dem_ta),'len(e1a_dem)=',len(e1a_dem))
    ax.plot(dem_ta/tdiv,e1a_dem,blumrk)
    ax.plot(dem_ta/tdiv,e1a_ocean,redmrk)
    ax.hold(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a')

    ax = fig.add_subplot(4,1,2)
    ax.hold(True)
    ax.plot(dem_tb/tdiv,e1b_dem,blumrk)
    ax.plot(dem_tb/tdiv,e1b_ocean,redmrk)
    ax.hold(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b')

    ax = fig.add_subplot(4,1,3)
    ax.hold(True)
    ax.plot(dem_ta/tdiv,e2a_dem,blumrk)
    ax.plot(dem_ta/tdiv,e2a_ocean,redmrk)
    ax.hold(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a')

    ax = fig.add_subplot(4,1,4)
    ax.hold(True)
    ax.plot(dem_tb/tdiv,e2b_dem,blumrk)
    ax.plot(dem_tb/tdiv,e2b_ocean,redmrk)
    ax.hold(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
#   ax.xaxis.set_ticklabels([])
    plt.ylabel('e2b')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

  if do_plt_uv:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.date+'-'+options.site+'-uv'+add_ndays_filt
    pltnamq = 'HPIES-uv'
    if options.verbose:
      print('pltnam=',pltnam)
    fig.suptitle(pltnam + add_filtdays)

    if deployment == 'OKMC':
      ylim = [-0.5,0.5]
    else:
      ylim = [-0.1,0.1]

    ax = fig.add_subplot(2,1,1)
    ax.plot(dem_tab/tdiv,u,blumrk)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
    ax.xaxis.set_ticklabels([])
    plt.ylabel('u, m/s')

    ax = fig.add_subplot(2,1,2)
    ax.plot(dem_tab/tdiv,v,blumrk)
    ax.hold(False)
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    ax.grid(True)
#   ax.xaxis.set_ticklabels([])
    plt.ylabel('v, m/s')
    plt.xlabel(tlbl)

    writepdf(pltdir, pltnam)
    os.system('cp ' + pltdir + '/' + pltnam + '.pdf ' + pltdirq + '/' + pltnamq + '.pdf')

if __name__ == '__main__':
  main()
