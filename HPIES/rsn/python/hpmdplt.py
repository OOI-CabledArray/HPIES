#! /usr/bin/env python2
# hpmdplt.py -- hpies multi-day plots of hprpro.py .mat output

from __future__ import print_function

"""
./hpmdplt.py --plt_motor --plt_clk_diffs --plt_ef_pr_tt --png --runid 2018 -s AB
./hpmdplt.py --plt_motor --plt_clk_diffs --plt_ef_pr_tt --png --runid 2018 -s SB
"""

import os
import sys
import time
import math
import numpy as np
from collections import namedtuple
from datetime import datetime,timedelta
# import glob
from time import strptime
from calendar import timegm
# from string import find
from optparse import OptionParser
from scipy.io import loadmat

def mymkdir(mydir):
  try:
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      sys.exit(1)

def writepng(pngdir, namewoext):
  mymkdir(pngdir)
  pngfile = os.path.join(pngdir,namewoext+'.png')
  print('pngfile=',pngfile)
  plt.savefig(pngfile)
  os.system('/home/dunlap/bin/updateframe.run ' + pngdir)

def writepdf(pltdir, namewoext):
  mymkdir(pltdir)
  pdffile = os.path.join(pltdir,namewoext+'.pdf')
  print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pltdir)

def main():
  global plt
  press_enter = 'Press Enter key for next plot, Ctrl-\\ to quit: '

  parser = OptionParser(usage='%prog [Options]')

  parser.add_option('-v', '--verbose',
    action='count', dest='verbose', default=0,
    help='print status messages to stdout')

  parser.add_option('--matdir', dest='matdir', metavar='D',
    default='/data/rsn/mat/hprpro/',
    help='directory with input .mat files')

# parser.add_option('--pngdir', dest='pngdir', metavar='D',
#   default='png', help='directory for .png files')

  parser.add_option('--pltdir', dest='pltdir', metavar='D',
    default='/data/hpies/rsn/plots/hpmdplt', help='directory for .pdf files')

  parser.add_option('--pdf',
    action='store_true', dest='pdf', default=False,
    help='make pdf files')

  parser.add_option('--png',
    action='store_true', dest='png', default=False,
    help='make png files')

  parser.add_option('--date_beg', dest='date_beg', metavar='D',
    default=None,
    help='beginning datetime, yyyy-mm-dd[_hh:mm:ss]')

  parser.add_option('--date_end', dest='date_end', metavar='D',
    default=None,
    help='ending datetime, yyyy-mm-dd[_hh:mm:ss]')

  parser.add_option('-d','--runid', dest='runid', metavar='D',
    default=None, help='run name')

  parser.add_option('-s','--site', dest='site', metavar='S',
    default=None, help='site name')

  parser.add_option('--plt_clk_diffs',
    action='store_true', dest='plt_clk_diffs', default=False,
    help='plot clock differences')

  parser.add_option('--plt_clk_diff_pres',
    action='store_true', dest='plt_clk_diff_pres', default=False,
    help='plot clock difference and pressure')

  parser.add_option('--plt_ef_pr_tt',
    action='store_true', dest='plt_ef_pr_tt', default=False,
    help='plot electric field, pressure and travel times')

  parser.add_option('--plt_motor',
    action='store_true', dest='plt_motor', default=False,
    help='plot motor and temperature')

  (options, args) = parser.parse_args()

  if options.pdf or options.png:
    import matplotlib
    matplotlib.use('Agg') # needed when used with cron or no DISPLAY or slow remote connection
    import matplotlib.pyplot as plt
    plot_interactive = False
  else:
    import matplotlib
    import matplotlib.pyplot as plt
    plot_interactive = True

  if len(args) != 0:
    parser.print_help()
    sys.exit(1)

  if options.runid is None and \
    (options.date_beg is None or options.date_end is None):
    parser.print_help()
    print('please provide --runid or both of --date_beg and --date_end')
    exit(1)

  ies_uxt_offset = []
  uobeg = []
  uoend = []

  if options.runid == 'test':
    options.date_beg = '2016-08-08'
    options.date_end = '2016-09-30'

  if options.runid == '2014':
    if options.site == 'AB':
      options.date_beg = '2014-07-31'
      options.date_end = '2015-07-22'
    if options.site == 'SB':
      options.date_beg = '2014-08-24'
      options.date_end = '2015-07-21'

  if options.runid == '2015':
    if options.site == 'AB':
      options.date_beg = '2015-07-28'
#     options.date_end = '2016-08-08'
      options.date_end = '2016-07-24'
    if options.site == 'SB':
      options.date_beg = '2015-07-23'
      options.date_end = '2016-08-08'

  if options.runid == '2016':
    if options.site == 'AB':
      options.date_beg = '2016-07-24'
      options.date_end = '2017-08-12'
    if options.site == 'SB':
      options.date_beg = '2016-08-08'
      options.date_end = '2017-08-12'

  if options.runid == '2017':
    if options.site == 'AB':
      options.date_beg = '2017-08-12'
      options.date_end = '2018-07-09'  # replaced
    if options.site == 'SB':
      options.date_beg = '2017-08-12'
      options.date_end = '2018-07-09'  # left in place

  if options.runid == '2018':
    if options.site == 'AB':
      options.date_beg = '2018-07-09'  # replaced
      options.date_end = '2019-06-17 21:00'
    if options.site == 'SB':
      options.date_beg = '2018-07-09'  # second year on same hardware
      options.date_end = '2019-06-17 21:00'

  if options.runid == '2019':
    if options.site == 'AB':
      options.date_beg = '2019-06-17 21:00'  # replaced
      options.date_end = '2020-08-05'
    if options.site == 'SB':
      options.date_beg = '2019-06-17 21:00'
      options.date_end = '2020-08-05' # not replaced this year

  if options.runid == '2020':
    if options.site == 'AB':
      options.date_beg = '2020-08-05'
      options.date_end = '2021-08-24' # not replaced this year
      uobeg = [timegm(strptime('2020-08-05', '%Y-%m-%d')),\
               timegm(strptime('2021-01-15', '%Y-%m-%d'))]
      uoend = [timegm(strptime('2020-08-07', '%Y-%m-%d')),\
               timegm(strptime('2021-08-24', '%Y-%m-%d'))]
      ies_uxt_offset = [76000000 - 680000 + 1600,\
                        90000000 - 590000 + 5300]
    if options.site == 'SB':
      options.date_beg = '2020-08-05' # not replaced last year
      options.date_end = '2021-08-04'

  if options.runid == '2021':
    if options.site == 'AB':
      options.date_beg = '2021-08-24' # continuosly running since 2020
      options.date_end = '2022-08-29'
      uobeg = [timegm(strptime(options.date_beg,'%Y-%m-%d'))]
      uoend = [timegm(strptime('2021-11-02',    '%Y-%m-%d'))]
      ies_uxt_offset = [90000000 - 600000 + 15000 + 150 + 12]
    if options.site == 'SB':
      options.date_beg = '2021-08-24'
      options.date_end = '2022-08-29'

  if options.runid == '2022':
    if options.site == 'AB':
      options.date_beg = '2022-08-29' # new instrument
      options.date_end = '2023-09-01'
    if options.site == 'SB':
      options.date_beg = '2022-08-29'
      options.date_end = '2023-09-01'


  print('hpmdplt.py runid=',options.runid,'site=',options.site,
        'beg=',options.date_beg,'end=',options.date_end)

  if len(options.date_beg) == 10:
    pyt_beg = strptime(options.date_beg,'%Y-%m-%d')
  elif len(options.date_beg) == 16:
    pyt_beg = strptime(options.date_beg,'%Y-%m-%d %H:%M')
  else:
    print("unknown options.date_beg=",options.date_beg)
    exit(1)

  if len(options.date_end) == 10:
    pyt_end = strptime(options.date_end,'%Y-%m-%d')
  elif len(options.date_end) == 16:
    pyt_end = strptime(options.date_end,'%Y-%m-%d %H:%M')
  else:
    print("unknown options.date_end=",options.date_end)
    exit(1)

  uxt_beg = timegm(pyt_beg)
  uxt_end = timegm(pyt_end)

  SOD = 86400
  day_first = int(uxt_beg / SOD)
  day_last  = int(uxt_end / SOD)

# /data/rsn/mat/hprpro/HPIES-20160924-SB-hprpro.mat

  dem_uxt     = np.array([], dtype='double')
  dem_e1_amp  = np.array([], dtype='double')
  dem_e2_amp  = np.array([], dtype='double')
  dem_abm     = np.array([], dtype='int')
  clk_diff_stm_ies = np.array([], dtype='double')
  clk_diff_stm_rsn = np.array([], dtype='double')

  aux_uxt  = np.array([], dtype='double')
  aux_pres = np.array([], dtype='double')
  aux_tt1  = np.array([], dtype='double')
  aux_tt2  = np.array([], dtype='double')
  aux_tt3  = np.array([], dtype='double')
  aux_tt4  = np.array([], dtype='double')

  count_mot12 = 0
  mot_abu   = np.array([], dtype='double')
  mot1_uxt  = np.array([], dtype='double')
  mot2_uxt  = np.array([], dtype='double')
  mot1_cur_max  = np.array([], dtype='double')
  mot2_cur_max  = np.array([], dtype='double')
  mot1_cur_med  = np.array([], dtype='double')
  mot2_cur_med  = np.array([], dtype='double')
  mot1_dur  = np.array([], dtype='double')
  mot2_dur  = np.array([], dtype='double')

  dem_nskip = 0
  dem_nuse  = 0
  aux_nskip = 0
  aux_nuse  = 0
  mot1_nuse  = 0
  mot2_nuse  = 0
  mot1_nskip  = 0
  mot2_nskip  = 0

  for day in range(day_first,day_last+1):
    uxt = day * SOD
    pyt = datetime(1970,1,1) + timedelta(0,uxt)
    ymd = pyt.strftime('%Y%m%d')
    file = 'HPIES-' + ymd + '-' + options.site + '-hprpro.mat'
#   print('file=',file)

    pathname = os.path.join(options.matdir,file)
    try:
      M = loadmat(pathname)
    except:
      M = None
#     print('cannot load',pathname)

    if M is not None:
      try:
        x = M['dem_uxt'][0]
        dem_ok = True
      except:
        dem_ok = False
        ok = False
        if options.verbose >= 2:
          print(pathname,'cannot get dem_uxt')

      try:
        x = M['aux_uxt_xfr'][0]
        aux_ok = True
      except:
        aux_ok = False
        if options.verbose >= 2:
          print(pathname,'cannot get aux_uxt_xfr')

      try:
        x = M['mot_uxt1'][0]
        mot1_ok = True
      except:
        mot1_ok = False
        if options.verbose >= 2:
          print(pathname,'cannot get mot_uxt1')

      try:
        x = M['mot_uxt2'][0]
        mot2_ok = True
      except:
        mot2_ok = False
        if options.verbose >= 2:
          print(pathname,'cannot get mot_uxt1')
    else:
      dem_ok = False
      aux_ok = False
      mot1_ok = False
      mot2_ok = False

    if dem_ok:
      dem_nuse += 1
    else:
      dem_nskip += 1

    if aux_ok:
      aux_nuse += 1
    else:
      aux_nskip += 1

    if mot1_ok:
      mot1_nuse += 1
    else:
      mot1_nskip += 1

    if mot2_ok:
      mot2_nuse += 1
    else:
      mot2_nskip += 1

    if dem_ok:
      dem_uxt =          np.append(dem_uxt,          M['dem_uxt'][0])
      dem_e1_amp =       np.append(dem_e1_amp,       M['dem_e1_amp'][0])
      dem_e2_amp =       np.append(dem_e2_amp,       M['dem_e2_amp'][0])
      dem_abm =          np.append(dem_abm,          M['dem_abm'][0])
      if len(M['dem_uxt'][0]) != len(M['dem_abm'][0]):
        print('warning: len(dem_uxt)=',len(dem_uxt),\
                   'and len(dem_abm)=',len(dem_abm),\
                   'not same length, file=',file)

    if aux_ok:
#     aux_uxt =          np.append(aux_uxt,          M['aux_uxt'][0])
      aux_uxt =          np.append(aux_uxt,          M['aux_uxt_xfr'][0])
      aux_pres =         np.append(aux_pres,         M['aux_pres'][0])
      aux_tt1 =          np.append(aux_tt1,          M['aux_tt1'][0])
      aux_tt2 =          np.append(aux_tt2,          M['aux_tt2'][0])
      aux_tt3 =          np.append(aux_tt3,          M['aux_tt3'][0])
      aux_tt4 =          np.append(aux_tt4,          M['aux_tt4'][0])
      clk_diff_stm_ies = np.append(clk_diff_stm_ies, M['clk_diff_stm_ies'][0])
      clk_diff_stm_rsn = np.append(clk_diff_stm_rsn, M['clk_diff_stm_rsn'][0])

    if mot1_ok and mot2_ok:
      count_mot12 += 1
      mot_abu  =         np.append(mot_abu , M['mot_abu'][0])
      mot1_uxt =         np.append(mot1_uxt, M['mot_uxt1'][0])
      mot2_uxt =         np.append(mot2_uxt, M['mot_uxt2'][0])
      try:
        mot1_cur_max =     np.append(mot1_cur_max, M['mot_cur1'][0])
        mot2_cur_max =     np.append(mot2_cur_max, M['mot_cur2'][0])
      except:
        mot1_cur_max =     np.append(mot1_cur_max, M['mot_cur1_max'][0])
        mot2_cur_max =     np.append(mot2_cur_max, M['mot_cur2_max'][0])
      try:
        mot1_cur_med =     np.append(mot1_cur_med, M['mot_cur1_med'][0])
        mot2_cur_med =     np.append(mot2_cur_med, M['mot_cur2_med'][0])
      except:
        mot1_cur_med =     np.append(mot1_cur_med, np.tile(np.NaN,np.shape(M['mot_uxt1'][0])))
        mot2_cur_med =     np.append(mot2_cur_med, np.tile(np.NaN,np.shape(M['mot_uxt2'][0])))
      mot1_dur =         np.append(mot1_dur, M['mot_dur1'][0])
      mot2_dur =         np.append(mot2_dur, M['mot_dur2'][0])

  if options.verbose:
    print('dem_nuse=',dem_nuse,'dem_nskip=',dem_nskip)
    print('aux_nuse=',aux_nuse,'aux_nskip=',aux_nskip)
    print('mot1_nuse=',mot1_nuse,'mot1_nskip=',mot1_nskip)
    print('mot2_nuse=',mot2_nuse,'mot2_nskip=',mot2_nskip)
    print('count_mot12=',count_mot12)



# print('shape(dem_abm)=',np.shape(dem_abm))
# print('shape(dem_uxt)=',np.shape(dem_uxt))
  ja = np.all([(dem_abm == ord('a')),(dem_uxt >= uxt_beg),(dem_uxt <= uxt_end)],axis=0)
  jb = np.all([(dem_abm == ord('b')),(dem_uxt >= uxt_beg),(dem_uxt <= uxt_end)],axis=0)
  if len(ja)<2 or len(jb)<2:
    print('error: len(ja)=',len(ja),'or len(jb)=',len(jb),'dem too short')
    return

  e1a = dem_e1_amp[ja]
  e1b = dem_e1_amp[jb]
  e2a = dem_e2_amp[ja]
  e2b = dem_e2_amp[jb]

  j = np.nonzero((aux_uxt <= uxt_beg) | (aux_uxt >= uxt_end))[0]
  if len(j) > 0:
    print('number aux_uxt out of bounds=',len(j))
    if options.verbose:
      auj = aux_uxt[j]
      print('  min(auj) =', datetime.fromtimestamp(np.min(auj)))
      print('  max(auj) =', datetime.fromtimestamp(np.max(auj)))
      dauj = np.diff(auj)
      print('  med(dauj)=',np.median(dauj))
      print('  min(dauj)=',np.min(dauj))
      print('  max(dauj)=',np.max(dauj))
      j = np.nonzero(dauj > np.median(dauj)*1.1)[0]
      for i in j:
        print('dauj[{0}] = {1}'.format(i,dauj[i]))
        print('  auj[{0}] = {1}'.format(i,datetime.fromtimestamp(auj[i])))
        print('  auj[{0}] = {1}'.format(i+1,datetime.fromtimestamp(auj[i+1])))

  j = np.nonzero((aux_uxt >= uxt_beg) & (aux_uxt <= uxt_end))[0]
  aux_uxt  = aux_uxt[j]
  aux_pres = aux_pres[j]
  aux_tt1  = aux_tt1[j]
  aux_tt2  = aux_tt2[j]
  aux_tt3  = aux_tt3[j]
  aux_tt4  = aux_tt4[j]
  clk_diff_stm_ies = clk_diff_stm_ies[j]
  clk_diff_stm_rsn = clk_diff_stm_rsn[j]

  if len(uobeg)>0:
    for i in range(len(uobeg)):
      print('i=',i,'ies_uxt_offset=',ies_uxt_offset[i])
      j = np.nonzero((aux_uxt >= uobeg[i]) & (aux_uxt <= uoend[i]))[0]
      clk_diff_stm_ies[j] -= ies_uxt_offset[i]
  

  clk_offset_ies = clk_diff_stm_ies - clk_diff_stm_rsn

  ttarr = [aux_tt1,aux_tt2,aux_tt3,aux_tt4]
  ttmed = np.median(ttarr,axis=0)

  mlta     = dem_uxt[ja] / 86400 + 719529 - 366
  mltb     = dem_uxt[jb] / 86400 + 719529 - 366
  aux_mlt  = aux_uxt     / 86400 + 719529 - 366
  mot1_mlt = mot1_uxt    / 86400 + 719529 - 366
  mot2_mlt = mot2_uxt    / 86400 + 719529 - 366

  
  try:
    mltbeg = np.min([np.min(mlta),np.min(mltb),np.min(aux_mlt),np.min(mot1_mlt),np.min(mot1_mlt)])
  except:
    mltbeg = np.min(mlta)

  try:
    mltend = np.max([np.max(mlta),np.max(mltb),np.max(aux_mlt),np.max(mot1_mlt),np.max(mot1_mlt)])
  except:
    mltend = np.max(mlta)

  r = mltend - mltbeg
  xlim = [mltbeg-r*0.02,mltend+r*0.02]

  if options.verbose >= 2:
    print('np.shape(dem_uxt)=',np.shape(dem_uxt))
    print('np.min(dem_uxt)=',np.min(dem_uxt))
    print('np.max(dem_uxt)=',np.max(dem_uxt))

    print('np.shape(ja)=',np.shape(ja))
    print('np.shape(jb)=',np.shape(jb))

    print('np.shape(mlta)=',np.shape(mlta))
    print('np.shape(mltb)=',np.shape(mltb))

  if options.verbose >= 1:
    dem_pyt_mina = datetime(1970,1,1) + timedelta(0,np.min(dem_uxt[ja]))
    dem_pyt_maxa = datetime(1970,1,1) + timedelta(0,np.max(dem_uxt[ja]))
    dem_pyt_minb = datetime(1970,1,1) + timedelta(0,np.min(dem_uxt[jb]))
    dem_pyt_maxb = datetime(1970,1,1) + timedelta(0,np.max(dem_uxt[jb]))
    
    print('dem_mina=',dem_pyt_mina.strftime('%Y-%m-%d_%H:%M:%S'))
    print('dem_maxa=',dem_pyt_maxa.strftime('%Y-%m-%d_%H:%M:%S'))
    print('dem_minb=',dem_pyt_minb.strftime('%Y-%m-%d_%H:%M:%S'))
    print('dem_maxb=',dem_pyt_maxb.strftime('%Y-%m-%d_%H:%M:%S'))

    if len(aux_uxt):
      aux_pyt_min = datetime(1970,1,1) + timedelta(0,np.min(aux_uxt))
      aux_pyt_max = datetime(1970,1,1) + timedelta(0,np.max(aux_uxt))

      print('aux_min =',aux_pyt_min.strftime('%Y-%m-%d_%H:%M:%S'))
      print('aux_max =',aux_pyt_max.strftime('%Y-%m-%d_%H:%M:%S'))

  mot1_ja = np.all([(mot_abu == ord('a')),(mot1_uxt >= uxt_beg),(mot1_uxt <= uxt_end)],axis=0)
  mot1_jb = np.all([(mot_abu == ord('b')),(mot1_uxt >= uxt_beg),(mot1_uxt <= uxt_end)],axis=0)
  if len(mot1_ja)<2 or len(mot1_jb)<2:
    print('error: len(mot1_ja)=',len(mot1_ja),'or len(mot1_jb)=',len(mot1_jb),\
          'mot1 ja,jb too short')
    return

  mot1a_mlt     = mot1_uxt    [mot1_ja] / 86400 + 719529 - 366
  mot1b_mlt     = mot1_uxt    [mot1_jb] / 86400 + 719529 - 366
  mot1a_cur_max = mot1_cur_max[mot1_ja]
  mot1b_cur_max = mot1_cur_max[mot1_jb]
  mot1a_cur_med = mot1_cur_med[mot1_ja]
  mot1b_cur_med = mot1_cur_med[mot1_jb]
  mot1a_dur     = mot1_dur    [mot1_ja]
  mot1b_dur     = mot1_dur    [mot1_jb]

  if options.verbose >= 1:
    mot1_pyt_min = datetime(1970,1,1) + timedelta(0,np.min(mot1_uxt))
    mot1_pyt_max = datetime(1970,1,1) + timedelta(0,np.max(mot1_uxt))
    mot2_pyt_min = datetime(1970,1,1) + timedelta(0,np.min(mot2_uxt))
    mot2_pyt_max = datetime(1970,1,1) + timedelta(0,np.max(mot2_uxt))

    print('mot1_min=',mot1_pyt_min.strftime('%Y-%m-%d_%H:%M:%S'))
    print('mot1_max=',mot1_pyt_max.strftime('%Y-%m-%d_%H:%M:%S'))
    print('mot2_min=',mot2_pyt_min.strftime('%Y-%m-%d_%H:%M:%S'))
    print('mot2_max=',mot2_pyt_max.strftime('%Y-%m-%d_%H:%M:%S'))

  mot2_ja = np.all([(mot_abu == ord('a')),(mot2_uxt >= uxt_beg),(mot2_uxt <= uxt_end)],axis=0)
  mot2_jb = np.all([(mot_abu == ord('b')),(mot2_uxt >= uxt_beg),(mot2_uxt <= uxt_end)],axis=0)
  if len(mot2_ja)<2 or len(mot2_jb)<2:
    print('error: len(mot2_ja)=',len(mot2_ja),'or len(mot2_jb)=',len(mot2_jb),\
          'mot2 ja,jb too short')
    return

  mot2a_mlt = mot2_uxt[mot2_ja] / 86400 + 719529 - 366
  mot2b_mlt = mot2_uxt[mot2_jb] / 86400 + 719529 - 366
  mot2a_cur_max = mot2_cur_max[mot2_ja]
  mot2b_cur_max = mot2_cur_max[mot2_jb]
  mot2a_cur_med = mot2_cur_med[mot2_ja]
  mot2b_cur_med = mot2_cur_med[mot2_jb]
  mot2a_dur = mot2_dur[mot2_ja]
  mot2b_dur = mot2_dur[mot2_jb]

  if options.verbose >= 2:
    print('len(mot1_ja)=',len(mot1_ja))
    print('len(mot1_jb)=',len(mot1_jb))
    print('len(mot2_ja)=',len(mot2_ja))
    print('len(mot2_jb)=',len(mot2_jb))

  if options.plt_clk_diffs:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.runid+'-'+options.site+'-clk-diffs'

    fig.suptitle(pltnam+'\n Clock diffs')

    ax1 = fig.add_subplot(3,1,1)
    ax2 = fig.add_subplot(3,1,2, sharex=ax1)
    ax3 = fig.add_subplot(3,1,3, sharex=ax1)

    ax1.plot_date(aux_mlt,clk_diff_stm_rsn,'b.',markersize=1)
    ax1.set_ylabel('clk diff stm-rsn, s')
    ax1.grid(True)
    ax1.set_xlim(xlim)

    ax2.plot_date(aux_mlt,clk_diff_stm_ies,'b.',markersize=1)
    ax2.set_ylabel('clk diff stm-ies, s')
    ax2.grid(True)

    ax3.plot_date(aux_mlt,clk_offset_ies,'b.',markersize=1)
    ax3.set_ylabel('clk offset rsn-ies, s')
    ax3.grid(True)

    fig.autofmt_xdate()
#   axlist = {ax1,ax2,ax3}
#   fix_xdates(axlist,1)
#   adj_xdates(axlist,aux_mlt)

    if options.pdf:
      writepdf(options.pltdir,pltnam)
    if options.png:
      writepng(options.pltdir,pltnam)
    if not options.pdf and not options.png:
      plt.show()
  
  if options.plt_clk_diff_pres:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.runid+'-'+options.site+'-clk-diff-pres'

    fig.suptitle(pltnam)
    fig.suptitle(pltnam+'\nClock diff and pressure')

    ax1 = fig.add_subplot(2,1,1)
    ax2 = fig.add_subplot(2,1,2, sharex=ax1)

    ax1.plot_date(aux_mlt, clk_offset_ies, 'b.', markersize=1)
    ax1.set_ylabel('clk diff rsn-ies, s')
    ax1.grid(True)
    ax1.set_xlim(xlim)

    ax2.plot_date(aux_mlt,aux_pres,'b.',markersize=1)
    ax2.set_ylabel('pressure, dbar')
    if options.site == 'SB':
      ax2.set_ylim([2940,2960])
    if options.site == 'AB':
      ax2.set_ylim([2655,2665])
    ax2.grid(True)

    fig.autofmt_xdate()

    if options.pdf:
      writepdf(options.pltdir,pltnam)
    if options.png:
      writepng(options.pltdir,pltnam)
    if not options.pdf and not options.png:
      plt.show()

  if options.plt_ef_pr_tt:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.runid+'-'+options.site+'-ef-pr-tt'

    fig.suptitle(pltnam+'\nHEF & PR unfiltered, TT median')

    ax1 = fig.add_subplot(4,1,1)
    ax2 = fig.add_subplot(4,1,2, sharex=ax1)
    ax3 = fig.add_subplot(4,1,3, sharex=ax1)
    ax4 = fig.add_subplot(4,1,4, sharex=ax1)

    ax1.plot_date(mlta, e1a, 'b.', markersize=1)
    ax1.plot_date(mltb, e1b, 'r.', markersize=1)
    ax1.set_ylabel('e1, uV')
    ax1.grid(True)
    ax1.set_xlim(xlim)
    ax1.set_ylim([-30,30])

    ax2.plot_date(mlta, e2a, 'b.', markersize=1)
    ax2.plot_date(mltb, e2b, 'r.', markersize=1)
    ax2.set_ylabel('e2, uV')
    ax2.grid(True)
    ax2.set_ylim([-30,30])

    ax3.plot_date(aux_mlt,aux_pres,'b.',markersize=1)
    ax3.set_ylabel('pressure, dbar')
    ax3.grid(True)
    if options.site == 'SB':
      ax3.set_ylim([2950,2960])
    if options.site == 'AB':
      ax3.set_ylim([2653,2663])

    ax4.plot_date(aux_mlt, ttmed, 'b.', markersize=1)
    ax4.set_ylabel('TT, s')
    ax4.grid(True)
    if options.site == 'SB':
      ax4.set_ylim([3.86,3.92])
    if options.site == 'AB':
      ax4.set_ylim([3.45,3.55])

    fig.autofmt_xdate()
    
    if options.pdf:
      writepdf(options.pltdir,pltnam)
    if options.png:
      writepng(options.pltdir,pltnam)
    if not options.pdf and not options.png:
      plt.show()

  if options.plt_motor:
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    pltnam = 'HPIES-'+options.runid+'-'+options.site+'-motor'

    fig.suptitle(pltnam+'\nMotor Current & Duration\nblu=A max, red=B max, cyan=A median, magenta=B median')

    ax1 = fig.add_subplot(6,1,1)
    ax2 = fig.add_subplot(6,1,2, sharex=ax1)
    ax3 = fig.add_subplot(6,1,3, sharex=ax1)
    ax4 = fig.add_subplot(6,1,4, sharex=ax1)
    ax5 = fig.add_subplot(6,1,5, sharex=ax1)
    ax6 = fig.add_subplot(6,1,6, sharex=ax1)

    if options.verbose >= 2:
      print('nok_mot1a_mlt=',len(np.nonzero(np.isfinite(mot1a_mlt))[0]))
      print('nok_mot1a_cur_max=',len(np.nonzero(np.isfinite(mot1a_cur_max))[0]))
      print('nok_mot1a_cur_med=',len(np.nonzero(np.isfinite(mot1a_cur_med))[0]))

    ax1.plot_date(mot1a_mlt, mot1a_cur_max, 'b.', markersize=1)
    ax1.plot_date(mot1b_mlt, mot1b_cur_max, 'r.', markersize=1)
    ax1.plot_date(mot1a_mlt, mot1a_cur_med, 'c.', markersize=1)
    ax1.plot_date(mot1b_mlt, mot1b_cur_med, 'm.', markersize=1)
    ax1.set_ylabel('mot1\ncur mA')
    ax1.grid(True)
    ax1.set_xlim(xlim)
    ax1.set_ylim([-100,1000])

    ax2.plot_date(mot2a_mlt, mot2a_cur_max, 'b.', markersize=1)
    ax2.plot_date(mot2b_mlt, mot2b_cur_max, 'r.', markersize=1)
    ax2.plot_date(mot2a_mlt, mot2a_cur_med, 'c.', markersize=1)
    ax2.plot_date(mot2b_mlt, mot2b_cur_med, 'm.', markersize=1)
    ax2.set_ylabel('mot2\ncur mA')
    ax2.grid(True)
    ax2.set_ylim([-100,1000])

    ax3.plot_date(mot1a_mlt, mot1a_dur, 'b.', markersize=1)
    ax3.plot_date(mot1b_mlt, mot1b_dur, 'r.', markersize=1)
    ax3.set_ylabel('mot1\ndur s')
    ax3.grid(True)
    ax3.set_ylim([-0.5,5.5])

    ax4.plot_date(mot2a_mlt, mot2a_dur, 'b.', markersize=1)
    ax4.plot_date(mot2b_mlt, mot2b_dur, 'r.', markersize=1)
    ax4.set_ylabel('mot2\ndur s')
    ax4.grid(True)
    ax4.set_ylim([-0.5,5.5])

    ax5.plot_date(mot1a_mlt[1:], np.diff(mot1a_mlt)*86400, 'b.', markersize=1)
    ax5.plot_date(mot1b_mlt[1:], np.diff(mot1b_mlt)*86400, 'r.', markersize=1)
    ax5.set_ylabel('mot1\ndt s')
    ax5.grid(True)
    ax5.set_ylim([0,1000])

    ax6.plot_date(mot2a_mlt[1:], np.diff(mot2a_mlt)*86400, 'b.', markersize=1)
    ax6.plot_date(mot2b_mlt[1:], np.diff(mot2b_mlt)*86400, 'r.', markersize=1)
    ax6.set_ylabel('mot2\ndt s')
    ax6.grid(True)
    ax6.set_ylim([0,1000])

    fig.autofmt_xdate()
    
    if options.pdf:
      writepdf(options.pltdir,pltnam)
    if options.png:
      writepng(options.pltdir,pltnam)
    if plot_interactive:
      fig.show()
      raw_input(press_enter)

if __name__ == '__main__':
  main()
