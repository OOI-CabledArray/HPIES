#!/usr/bin/env python2
# hpproc.py -- decode, plot and process ascii files of raw data from hprtst.py and hefcf2.py
# hprtst.py is a real-time program for hpies-rsn
# hefcf2.py converts binary compact flash files from hefrun.c to ascii files like hprtst.py makes
# see ../man/users-manual.txt for a description of the files
# see ../man/processing.txt for a description of the processing
# John Dunlap, dunlap@apl.uw.edu
# updated May 20, 2014
# started life as hprdec.py then added processing

"""
48-hr tests:
./hpproc.py -d -k -a -s 15 -w 15 -y 7  \
  --swap_e2bc --swap_e2ab \
  --pltnambase ws30-31-48hr \
  hprtst-48hr

./hpproc.py -d -k -a -s 15 -w 15 -y 7 \
  --pltnambase ws32-33-48hr \
  hprtst.out

./hpproc.py -e -s 15 -l 60 -n 15 -p 3 \
  hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
  Plots 15 half-cycles (each 240 s) of raw ef data in each frame.
  Also plots above minus 3rd-order polynomial fit to unpinched data.
  First 15 s not shown.
  Note that 15 s is 5 time-constants of preamp which has 1 pole RC at 3 s.
  pdf/raw-ef-cat

./hpproc.py -c -s 15 hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
  plots resistance check (aka cal), one half-cycle per pdf
  pdf/raw-cal/

./hpproc.py -e -s 15 -l 10 hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
  plots raw ef, one half-cycle per pdf
  pdf/raw-ef/

./hpproc.py -m -l 10 hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
  plots raw motor current, one half-cycle per pdf
  pdf/raw-mot/

./hpproc.py -n 10 -m -l 50 hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
  plots concatenated raw motor current
  pdf/raw-mot-cat/

./hpproc.py -a -s 15 hefcf2-odir/H[12345]-{010,050,100,150,200,250,300,350}.gz
./hpproc.py -a -s 15 hefcf2-odir/H2-010.gz
  plots averages, skipping first 15 s after each pinch
  pdf/hef-avg/

./hpproc.py run1
./hpproc.py -s 15 -o hefcf2-odir/H2-010

./hpproc.py -m -l 10 tests/apr15a       # raw motor current
./hpproc.py -c -s 3 tests/apr15a       # raw cal
./hpproc.py -e -s 3 tests/apr15a       # raw ef
./hpproc.py -e -s 3 -n 5 tests/apr15a  # raw ef concatenated for several HCY
  Tests Apr 15, 2014 mods to stm and hef firmware on HEF-006, IES-309

./hpproc.py -g hppdir hefcf2-odir/H2-010 # output gzip files

./hpprocx.py hprtst.outz -k -a -w 5 -y 7 # testing July 11, 2014

"""

from __future__ import print_function

import os
import sys
from optparse import OptionParser
import time
import math
import numpy as np
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

from crclib import crc3kerm
from crclib import chksumnmea
from get_hpies_info import get_hpies_info

def find_all(a_str, sub):
  start = 0
  while True:
    start = a_str.find(sub, start)
    if start == -1:
      return
    yield start
    start += len(sub)

def check_crc_aux(line):
  return True
  i = line.find('AUX,')
  if i<0:
    return
  j = list(find_all(line, ','))
  print('check_crc_aux: i=',i,'j=',j)
  print('line=',line)
  sys.exit(0)
  return True

def check_crc(line):
  global count_bad_crc_check, count_crc_ok, count_crc_missing
  global count_bad_crc_decode

  i = line.find('#')
  j = line.rfind('*')
  if i<0 or j<0:
    count_crc_missing += 1
    return False
  sect = line[i:j]
  # print('i=',i,'j=',j,'tst=['+tst+']')
  crc_computed = crc3kerm(sect)

  try:
    crc_rcv = int(line[j+1:j+5],16)
  except:
    crc_rcv = None
    print('bad crc decode, lineno=',lineno,'linein=['+linein+']')
    count_bad_crc_decode += 1
    return False

  if crc_rcv != None:
    if crc_computed != crc_rcv:
      print('crc error, lineno=',lineno,'linein=[' + linein + ']')
      count_bad_crc_check += 1
      return False
    else:
      count_crc_ok += 1
      return True

def split_aux(line):
  global aux_split
  global ies_secs_diff

  i = line.find('AUX,')
  j = line.rfind('*')
  if i<0 or j<0:
    aux_split = None
#   print('bad AUX line')
#   print_info_and_exit()
    return
  aux_split = line[i:j].split(',')

  if len(aux_split) != 13:
    print('AUX wrong split len:',len(aux_split))
    print('lineno=',lineno)
    print('line=',line)
    print('aux_split=',aux_split)
    aux_split = None
    return False

  try:
    ies_secs_diff = int(aux_split[1]) - int(aux_split[12])
  except:
    print('cannot decode aux split')
    ies_secs_diff = None
    return False

  return True
  
def split_tod(line):
  global tod_split
  global stm_secs_diff
  i = line.find('#')
  j = line.rfind('*')
  if i<0 or j<0:
    print('bad TOD line')
    print_info_and_exit()
#   tod_split = None
#   return
  tod_split = line[i:j].split(',')
  if len(tod_split) == 3:
    try:
      stm_secs_diff = int(tod_split[2]) - int(tod_split[1])
    except:
      print('cannot decode TOD line')
      print_info_and_exit()

def split_hef(line):
  global hef_split
  global count_removed

  i = line.find('#')
  j = line.rfind('*')
  if i<0 or j<0:
    print('bad HEF line')
    print_info_and_exit()
#   hef_split = None
#   return
  while j>2 and line[j-2] == '\\' and (line[j-1] == 'r' or line[j-1] == 'n'):
    j -= 2
    count_removed += 1
  if i<j:
    hef_split = line[i:j].split()
  else:
    print('hpproc: i<j')
    print_info_and_exit()

def decode_hef_hdr():
  global HEF_HDR
  global count_EC
  global count_M, count_E
  global count_C, count_hef_unknown

  if len(hef_split) < 1:
    print('len(hef_split)=',len(hef_split))
    print_info_and_exit()

  # check if this is a header for EF, cal or motor
  if hef_split[0].find('#3__HE') == 0:
    if HEF_HDR.hef_iscan != -1:
      print('decode_hef_hdr(): hef_iscan=',HEF_HDR.hef_iscan)
      return
      # print_info_and_exit()

    global hef_ntoss
    if hef_ntoss > 0:
      print('previous hef_ntoss=',hef_ntoss)
    hef_ntoss = 0

    if hef_split[0] == '#3__HE05':
      if len(hef_split) != 12:
        print('HE05 wrong split len:',len(hef_split))
        print('hef_split:',hef_split)
        print_info_and_exit()

      if options.verbose:
        print('hpproc: HE:',linein)
        sys.stdout.flush()

      try:
        HEF_HDR.typ      = hef_split[1]
        HEF_HDR.abu       = hef_split[2]
        HEF_HDR.ibeg     = int(hef_split[3])
        HEF_HDR.iend     = int(hef_split[4])
        HEF_HDR.hcno     = int(hef_split[5])
        HEF_HDR.uxt      = int(hef_split[6])
        HEF_HDR.ticks    = int(hef_split[7])
        HEF_HDR.navg_mot = int(hef_split[8])
        HEF_HDR.navg_ef  = int(hef_split[9])
        HEF_HDR.navg_cal = int(hef_split[10])
      except:
        print('error decoding HE05')
        print_info_and_exit()

      # global nskip_ef, nskip_cal
      HEF_HDR.nskip_ef  = int(0.5 + META.tskip / (HEF_HDR.navg_ef  * META.tsamp_ef))
      HEF_HDR.nskip_cal = int(0.5 + META.tskip / (HEF_HDR.navg_cal * META.tsamp_ef))
      HEF_HDR.nwait_cal = int(0.5 + META.twait / (HEF_HDR.navg_cal * META.tsamp_ef))

      # global tskip_ef, tskip_cal
      HEF_HDR.tskip_ef  = HEF_HDR.nskip_ef  * HEF_HDR.navg_ef   * META.tsamp_ef
      HEF_HDR.tskip_cal = HEF_HDR.nskip_cal * HEF_HDR.navg_cal  * META.tsamp_ef

      # print('nwait_cal=',HEF_HDR.nwait_cal)
      # print('navg_cal=',HEF_HDR.navg_cal)
      # print('tsamp_ef=',META.tsamp_ef)
      HEF_HDR.twait_cal = HEF_HDR.nwait_cal * HEF_HDR.navg_cal  * META.tsamp_ef

      vref = 2.5       # ADC Vref
      preamp = 800     # preamp gain
      divider = 0.1773 # resistor network at ADC input
      HEF_HDR.off_ef  = pow(2,15) * HEF_HDR.navg_ef
      HEF_HDR.off_cal = pow(2,15) * HEF_HDR.navg_cal
      HEF_HDR.sf_ef  = vref / pow(2,16) / preamp / divider * 1e6 / HEF_HDR.navg_ef  # uV
      HEF_HDR.sf_cal = vref / pow(2,16) / preamp / divider * 1e6 / HEF_HDR.navg_cal # uV
      HEF_HDR.sf_mot = vref / pow(2,16) * 1e3 / HEF_HDR.navg_mot                    # mA

      if HEF_HDR.uxt < 1372636800:                    # July 1, 2013
        HEF_HDR.uxt_ibeg = HEF_HDR.uxt - (HEF_HDR.iend-HEF_HDR.ibeg) * META.tsamp_ef # OKMC & before
      elif HEF_HDR.uxt > 1396310400:                  # April 1, 2014
        HEF_HDR.uxt_ibeg = HEF_HDR.uxt                        # RSN & later
      else:
        print('not sure how to compute uxt_ibeg between OKMC and RSN')
        print_info_and_exit()

      if INFO.uxt_poly_ref:
        # align HEF time to IES time using polynomial fits
        HEF_HDR.uxt_ibeg = HEF_HDR.uxt_ibeg - np.polyval(INFO.uxt_poly_coef, HEF_HDR.uxt_ibeg - INFO.uxt_poly_ref)

      if HEF_HDR.typ == 'f' or HEF_HDR.typ == 'r':
        count_M += 1
        navg = HEF_HDR.navg_mot
      elif HEF_HDR.typ == 'E':
        count_E += 1
        count_EC += 1
        navg = HEF_HDR.navg_ef
      elif HEF_HDR.typ == 'C':
        count_C += 1
        count_EC += 1
        navg = HEF_HDR.navg_cal
      else:
        count_hef_unknown += 1
        print('unknown HEF hdr type')
        print_info_and_exit()

      if HEF_HDR.abu != 'a' and HEF_HDR.abu != 'b' and HEF_HDR.abu != 'u':
        print('HEF_HDR.abu should be a, b or u')
        print_info_and_exit()

      if HEF_HDR.navg_mot < 1 or HEF_HDR.navg_ef < 1 or HEF_HDR.navg_cal < 1:
        print('HEF_HDR: navg_mot, navg_ef, navg_cal must all be > 0')
        print_info_and_exit()

      HEF_HDR.hef_nscan = int((HEF_HDR.iend - HEF_HDR.ibeg) / navg)
      HEF_HDR.hef_iscan = 0

    else:
      print('unknown hef_hdr')
      print_info_and_exit()

def decode_cal_status():
  global CAL_STATUS, plot_raw_cal_flag

  if hef_split[0].find('#3__SC') == 0:
    if options.verbose:
      print('found #3_SC')
      sys.stdout.flush()
    CAL_STATUS = collections.namedtuple('CAL_STATUS',[])
    if len(hef_split) != 4:
      print('cal status header record should have 4 tokens')
      print('cal status hdr:',hef_split)
      print_info_and_exit()
    CAL_STATUS.i0 = int(hef_split[1])
    CAL_STATUS.ns = int(hef_split[2])
    CAL_STATUS.nl = int(hef_split[3])
    CAL_STATUS.lno = 0
    CAL_STATUS.j = []
    CAL_STATUS.s = []
    
  if hef_split[0].find('#3__Sc') == 0:
    if len(hef_split) != 3:
      print('cal status data record should have 3 tokens')
      print('cal status data:',hef_split)
      print_info_and_exit()
    CAL_STATUS.lno += 1
    CAL_STATUS.j.append(int(hef_split[1]))
    CAL_STATUS.s.append(int(hef_split[2]))

    if CAL_STATUS.lno == CAL_STATUS.nl:
      if options.verbose:
        print('before expand_cal_status()')
        sys.stdout.flush()
      expand_cal_status()
      if options.verbose:
        print('finished expand_cal_status()')
        sys.stdout.flush()
      compute_cal_demod()
      if options.verbose:
        print('finished compute_cal_demod()')
        sys.stdout.flush()
      if plot_raw_cal_flag:
        plot_raw_cal_flag = False
        if options.limEC == 0 or count_EC <= options.limEC:
          if options.ncat == 0:
            plot_raw_cal()
          else:
            plot_raw_cal_cat()

def expand_cal_status():
  global CAL_STATUS

  CAL_STATUS.stat_fast = np.tile(np.nan,CAL_STATUS.ns)
  jn = 0
  sn = CAL_STATUS.s[0]
  for i in range(1,CAL_STATUS.nl):
    jp = jn
    jn = CAL_STATUS.j[i]
    sp = sn
    sn = CAL_STATUS.s[i]
    CAL_STATUS.stat_fast[range(jp,jn)] = np.tile(sp,jn-jp)
  CAL_STATUS.stat_fast[range(jn,CAL_STATUS.ns)] = np.tile(sn,CAL_STATUS.ns-jn)
  if len(CAL_STATUS.stat_fast) != CAL_STATUS.ns:
    print('CAL_STATUS.stat_fast has wrong length')
    print_info_and_exit()
  jnz = np.nonzero(np.isnan(CAL_STATUS.stat_fast))[0]
  if len(jnz) > 0:
    print('not all CAL_STATUS.stat_fast filled in')
    print('  len(jnz)=',len(jnz))
    print('  jnz=',jnz)
    print_info_and_exit()
  try:
    CAL_STATUS.hcyc_status = CAL_STATUS.stat_fast[HCYC.ind];
  except:
    print('cannot compute CAL_STATUS.hcyc_status')
    print_info_and_exit()
  CAL_STATUS.neg = np.nonzero(CAL_STATUS.hcyc_status == 16)[0]
  CAL_STATUS.pos = np.nonzero(CAL_STATUS.hcyc_status == 17)[0]
  if len(CAL_STATUS.neg) + len(CAL_STATUS.pos) != len(CAL_STATUS.hcyc_status):
    print('some unknown values in CAL_STATUS.hcyc_status')
    print_info_and_exit()
    
def decode_mot_status():
  global MOT_STATUS, plot_raw_mot_flag

  if hef_split[0].find('#3__SM') == 0:
    MOT_STATUS = collections.namedtuple('MOT_STATUS',[])
    if len(hef_split) != 4:
      print('mot status header record should have 4 tokens')
      print('mot status hdr:',hef_split)
      print_info_and_exit()
    MOT_STATUS.i0 = int(hef_split[1])
    MOT_STATUS.ns = int(hef_split[2])
    MOT_STATUS.nl = int(hef_split[3])
    MOT_STATUS.lno = 0
    MOT_STATUS.j = []
    MOT_STATUS.s = []
    
  if hef_split[0].find('#3__Sm') == 0:
    if len(hef_split) != 3:
      print('mot status data record should have 3 tokens')
      print('mot status data:',hef_split)
      print_info_and_exit()
    MOT_STATUS.lno += 1
    # j is index where s first appears:
    MOT_STATUS.j.append(int(hef_split[1]))
    MOT_STATUS.s.append(int(hef_split[2]))

    if MOT_STATUS.lno == MOT_STATUS.nl:
      expand_mot_status()
      accum_mot_cur_dur()
      if plot_raw_mot_flag:
        plot_raw_mot_flag = False
        if options.limEC == 0 or count_EC <= options.limEC:
          if options.ncat == 0:
            plot_raw_mot()
          else:
            plot_raw_mot_cat()

def expand_mot_status():
  global MOT_STATUS

  MOT_STATUS.stat_fast = np.tile(np.nan,MOT_STATUS.ns)
  jn = 0
  sn = MOT_STATUS.s[0]
  for i in range(1,MOT_STATUS.nl):
    jp = jn
    jn = MOT_STATUS.j[i]
    sp = sn
    sn = MOT_STATUS.s[i]
    # now jn is index where sn first appears
    MOT_STATUS.stat_fast[range(jp,jn)] = np.tile(sp,jn-jp)
  MOT_STATUS.stat_fast[range(jn,MOT_STATUS.ns)] = np.tile(sn,MOT_STATUS.ns-jn)
  if len(MOT_STATUS.stat_fast) != MOT_STATUS.ns:
    print('MOT_STATUS.stat_fast has wrong length')
    print_info_and_exit()
  jnz = np.nonzero(np.isnan(MOT_STATUS.stat_fast))[0]
  if len(jnz) > 0:
    print('not all MOT_STATUS.stat_fast filled in')
    print('  len(jnz)=',len(jnz))
    print('  jnz=',jnz)
    print_info_and_exit()
  try:
    MOT_STATUS.stat_pick = MOT_STATUS.stat_fast[HCYM.ind];
  except:
    print('cannot compute MOT_STATUS.stat_pick')
    print_info_and_exit()

def append_hef_data():
  global mot_ind_a, mot_cur_a
  global mot_ind_b, mot_cur_b
  global mot_ind_u, mot_cur_u
  global hef_ind_a, hef_e1a_a, hef_e1b_a, hef_e1c_a, hef_e2a_a, hef_e2b_a, hef_e2c_a
  global hef_ind_b, hef_e1a_b, hef_e1b_b, hef_e1c_b, hef_e2a_b, hef_e2b_b, hef_e2c_b
  global cal_ind_a, cal_e1a_a, cal_e1b_a, cal_e1c_a, cal_e2a_a, cal_e2b_a, cal_e2c_a
  global cal_ind_b, cal_e1a_b, cal_e1b_b, cal_e1c_b, cal_e2a_b, cal_e2b_b, cal_e2c_b
# global hef_iscan
  global HCYE  # half-cycle EF
  global HCYC  # half-cycle CAL
  global HCYM  # half-cycle MOT
  global plot_raw_cal_flag
  global plot_raw_mot_flag

  if hef_split[0].find('#3__D')==0 and HEF_HDR.hef_iscan < 0:
    global hef_ntoss
    hef_ntoss += 1
    if hef_ntoss == 1:
      print('append_hef_data(): tossed data record without header record')
      print('  lineno=',lineno)
      print('  linein=',linein)
    return

  # header for EF, cal or motor
  if hef_split[0].find('#3__HE')==0:
    # add gaps at each HEF header to signify a break in the sampling
    mot_ind_a.append(None)
    mot_cur_a.append(None)

    hef_ind_a.append(None)
    hef_e1a_a.append(None)
    hef_e1b_a.append(None)
    hef_e1c_a.append(None)
    hef_e2a_a.append(None)
    hef_e2b_a.append(None)
    hef_e2c_a.append(None)

    cal_ind_a.append(None)
    cal_e1a_a.append(None)
    cal_e1b_a.append(None)
    cal_e1c_a.append(None)
    cal_e2a_a.append(None)
    cal_e2b_a.append(None)
    cal_e2c_a.append(None)

    mot_ind_b.append(None)
    mot_cur_b.append(None)
    mot_cur_u.append(None)

    hef_ind_b.append(None)
    hef_e1a_b.append(None)
    hef_e1b_b.append(None)
    hef_e1c_b.append(None)
    hef_e2a_b.append(None)
    hef_e2b_b.append(None)
    hef_e2c_b.append(None)

    cal_ind_b.append(None)
    cal_e1a_b.append(None)
    cal_e1b_b.append(None)
    cal_e1c_b.append(None)
    cal_e2a_b.append(None)
    cal_e2b_b.append(None)
    cal_e2c_b.append(None)

  # motor current data
  if hef_split[0] == '#3__DM':

    try:
      ind = int(hef_split[1])
      cur = int(hef_split[2])
    except:
      print('cannot decode motor current data')
      print_info_and_exit()

    cur *= HEF_HDR.sf_mot

    if HEF_HDR.abu == 'a':
      mot_ind_a.append(ind)
      mot_cur_a.append(cur)
    elif HEF_HDR.abu == 'b':
      mot_ind_b.append(ind)
      mot_cur_b.append(cur)
    elif HEF_HDR.abu == 'u':
      mot_ind_u.append(ind)
      mot_cur_u.append(cur)
    else:
      print('unknown abu=',HEF_HDR.abu)
      print_info_and_exit()

    # half-cycle arrays
    if HEF_HDR.hef_iscan == 0:
      HCYM = collections.namedtuple('HCYM',[])
      HCYM.uxt = []
      HCYM.ind = []
      HCYM.cur = []
    HCYM.ind.append(ind)
    HCYM.cur.append(cur)

    HEF_HDR.hef_iscan += 1

    if HEF_HDR.hef_iscan == HEF_HDR.hef_nscan:
      HEF_HDR.hef_iscan = -1

    if HEF_HDR.hef_iscan == -1:
      HCYM.abu = HEF_HDR.abu

      # HEF UXT of start of motor move:
      HCYM.uxt0 = HEF_HDR.uxt_ibeg

      # elapsed seconds of each data point since start of motor move:
      HCYM.secs =  np.array(HCYM.ind,dtype='double') * META.tsamp_mot

      # UXT of of each sample:
      HCYM.uxt =  HCYM.secs + HCYM.uxt0
      HCYM.cur = np.array(HCYM.cur,dtype='double')

    if HEF_HDR.hef_iscan == -1 and options.do_plot_raw_mot:
      if options.limEC == 0 or count_EC <= options.limEC:
        plot_raw_mot_flag = True # defer plot_raw_mot() until get status

  # EF data
  elif hef_split[0] == '#3__DE':
    try:
      # convert ADC counts to microvolts at preamp input
      # flip sign of "b" because "b" preamp input opposite of "a"
      ind =   int(hef_split[1])
      e1a =  (int(hef_split[3]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
      e1b = -(int(hef_split[4]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
      e1c =  (int(hef_split[2]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
      e2a =  (int(hef_split[6]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
      e2b = -(int(hef_split[7]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
      e2c =  (int(hef_split[5]) - HEF_HDR.off_ef) * HEF_HDR.sf_ef
    except:
      print('EF data decode failed, lineno=',lineno,'linein=',linein)
      print('hef_split=',hef_split)
      print_info_and_exit()

    if options.swap_e1bc:
      x = -e1b
      e1b = -e1c
      e1c = x
      
    if options.swap_e1ab:
      x = e1b
      e1b = e1a
      e1a = x
      
    if options.swap_e2bc:
      x = -e2b
      e2b = -e2c
      e2c = x
      
    if options.swap_e2ab:
      x = e2b
      e2b = e2a
      e2a = x
      

    if ind != (HEF_HDR.hef_iscan + 1) * HEF_HDR.navg_ef - 1:
      print('hpproc.py: ind != (HEF_HDR.hef_iscan + 1) * HEF_HDR.navg_ef - 1')
      print('  ind=',ind)
      print('  hef_iscan=',HEF_HDR.hef_iscan)
      print('  navg_ef=',HEF_HDR.navg_ef)
      HEF_HDR.hef_iscan = -1
      return
      # print_info_and_exit()

    # file-length arrays
    if HEF_HDR.abu == 'a':
      hef_ind_a.append(ind)
      hef_e1a_a.append(e1a)
      hef_e1b_a.append(e1b)
      hef_e1c_a.append(e1c)
      hef_e2a_a.append(e2a)
      hef_e2b_a.append(e2b)
      hef_e2c_a.append(e2c)
    elif HEF_HDR.abu == 'b':
      hef_ind_b.append(ind)
      hef_e1a_b.append(e1a)
      hef_e1b_b.append(e1b)
      hef_e1c_b.append(e1c)
      hef_e2a_b.append(e2a)
      hef_e2b_b.append(e2b)
      hef_e2c_b.append(e2c)
    elif HEF_HDR.abu == 'u':
      hef_ind_u.append(ind)
      hef_e1a_u.append(e1a)
      hef_e1b_u.append(e1b)
      hef_e1c_u.append(e1c)
      hef_e2a_u.append(e2a)
      hef_e2b_u.append(e2b)
      hef_e2c_u.append(e2c)
    else:
      print('hpproc.py: unknown abu=',HEF_HDR.abu)
      print_info_and_exit()

    # half-cycle arrays
    if HEF_HDR.hef_iscan == 0:
      HCYE = collections.namedtuple('HCYE',[])
      HCYE.ind = []
      HCYE.e1a = []
      HCYE.e1b = []
      HCYE.e1c = []
      HCYE.e2a = []
      HCYE.e2b = []
      HCYE.e2c = []
    try:
      HCYE.ind.append(ind)
      HCYE.e1a.append(e1a)
      HCYE.e1b.append(e1b)
      HCYE.e1c.append(e1c)
      HCYE.e2a.append(e2a)
      HCYE.e2b.append(e2b)
      HCYE.e2c.append(e2c)
    except:
      print('hpproc.py: cannot append HCYE:')
      print('  hef_iscan=',HEF_HDR.hef_iscan)
      print('  lineno=',lineno)
      print('  linein=',linein)
      sys.exit(1)

    HEF_HDR.hef_iscan += 1

    if HEF_HDR.hef_iscan == HEF_HDR.hef_nscan:
      HEF_HDR.hef_iscan = -1

    if HEF_HDR.hef_iscan == -1:
      HCYE.abu = HEF_HDR.abu
      HCYE.hcno = HEF_HDR.hcno

      # HEF UXT of end of motor move:
      HCYE.uxt0 = HEF_HDR.uxt_ibeg

      # elapsed seconds of each data point since motor move:
      HCYE.secs =  np.array(HCYE.ind,dtype='double') * META.tsamp_ef

      # HEF UXT of of each sample:
      HCYE.uxt =  HCYE.secs + HCYE.uxt0

      HCYE.e1a = np.array(HCYE.e1a,dtype='double')
      HCYE.e1b = np.array(HCYE.e1b,dtype='double')
      HCYE.e1c = np.array(HCYE.e1c,dtype='double')
      HCYE.e2a = np.array(HCYE.e2a,dtype='double')
      HCYE.e2b = np.array(HCYE.e2b,dtype='double')
      HCYE.e2c = np.array(HCYE.e2c,dtype='double')

      HCYE.cou = count_EC

      compute_hef_demod(HCYE)

    if HEF_HDR.hef_iscan == -1 and options.do_extra:
      print_uxt_ies_chk()

    if HEF_HDR.hef_iscan == -1 and options.do_plot_raw_ef:
      if options.limEC == 0 or count_EC <= options.limEC:
        if options.ncat == 0:
          compute_HCYE()
          plot_raw_ef()
        else:
          plot_raw_ef_cat()

#   if HEF_HDR.hef_iscan == -1 and options.do_plot_hist_res:
#     compute_HCYE()
#     accum_histogram()

    # compute averages of HEF data
    if HEF_HDR.hef_iscan == -1 and HEF_HDR.nskip_ef < HEF_HDR.hef_nscan:
#     print('hpproc.py: hef_iscan=',HEF_HDR.hef_iscan,'hef_nscan=',HEF_HDR.hef_nscan)

      j = np.array(range(HEF_HDR.nskip_ef,HEF_HDR.hef_nscan))

      global HEF_AVG
      HEF_AVG.hcno.append(HEF_HDR.hcno)         # half-cycle number
      HEF_AVG.abu.append  (HEF_HDR.abu)           # A or B pinched
      HEF_AVG.navg.append(HEF_HDR.navg_ef)      # number of hardware samples averaged
      HEF_AVG.nuse.append(len(j))               # number of points used in averages
      HEF_AVG.uxt.append (np.mean(HCYE.uxt[j])) # unix seconds
      HEF_AVG.e1am.append(np.mean(HCYE.e1a[j])) # mean
      HEF_AVG.e1bm.append(np.mean(HCYE.e1b[j]))
      HEF_AVG.e1cm.append(np.mean(HCYE.e1c[j]))
      HEF_AVG.e2am.append(np.mean(HCYE.e2a[j]))
      HEF_AVG.e2bm.append(np.mean(HCYE.e2b[j]))
      HEF_AVG.e2cm.append(np.mean(HCYE.e2c[j]))
      HEF_AVG.e1as.append(np.std(HCYE.e1a[j]))  # std dev
      HEF_AVG.e1bs.append(np.std(HCYE.e1b[j]))
      HEF_AVG.e1cs.append(np.std(HCYE.e1c[j]))
      HEF_AVG.e2as.append(np.std(HCYE.e2a[j]))
      HEF_AVG.e2bs.append(np.std(HCYE.e2b[j]))
      HEF_AVG.e2cs.append(np.std(HCYE.e2c[j]))

      t = HCYE.uxt[j] - np.mean(HCYE.uxt[j])
      e1a_poly = np.polyfit(t,HCYE.e1a[j],1)
      e1b_poly = np.polyfit(t,HCYE.e1b[j],1)
      e1c_poly = np.polyfit(t,HCYE.e1c[j],1)
      e2a_poly = np.polyfit(t,HCYE.e2a[j],1)
      e2b_poly = np.polyfit(t,HCYE.e2b[j],1)
      e2c_poly = np.polyfit(t,HCYE.e2c[j],1)
      e1a_fit = np.polyval(e1a_poly,t)
      e1b_fit = np.polyval(e1b_poly,t)
      e1c_fit = np.polyval(e1c_poly,t)
      e2a_fit = np.polyval(e2a_poly,t)
      e2b_fit = np.polyval(e2b_poly,t)
      e2c_fit = np.polyval(e2c_poly,t)
      HCYE.e1a_res = HCYE.e1a[j] - e1a_fit
      HCYE.e1b_res = HCYE.e1b[j] - e1b_fit
      HCYE.e1c_res = HCYE.e1c[j] - e1c_fit
      HCYE.e2a_res = HCYE.e2a[j] - e2a_fit
      HCYE.e2b_res = HCYE.e2b[j] - e2b_fit
      HCYE.e2c_res = HCYE.e2c[j] - e2c_fit

      # fs suffix: fancy standard deviations
      HEF_AVG.e1afs.append(np.std(HCYE.e1a_res))
      HEF_AVG.e1bfs.append(np.std(HCYE.e1b_res))
      HEF_AVG.e1cfs.append(np.std(HCYE.e1c_res))
      HEF_AVG.e2afs.append(np.std(HCYE.e2a_res))
      HEF_AVG.e2bfs.append(np.std(HCYE.e2b_res))
      HEF_AVG.e2cfs.append(np.std(HCYE.e2c_res))

    if HEF_HDR.hef_iscan > HEF_HDR.hef_nscan:
      print('hpproc.py: hef_iscan > hef_nscan')
      print_info_and_exit()

  # calibration data (resistance check)
  elif hef_split[0] == '#3__DC':
#   if options.verbose:
#     print('#3__DC: hef_iscan=',HEF_HDR.hef_iscan)
    try:
      # convert ADC counts to microvolts at preamp input
      # flip sign of "b" because "b" preamp input opposite of "a"
      ind =  int(hef_split[1])
      e1a =  (int(hef_split[3]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
      e1b = -(int(hef_split[4]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
      e1c =  (int(hef_split[2]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
      e2a =  (int(hef_split[6]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
      e2b = -(int(hef_split[7]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
      e2c =  (int(hef_split[5]) - HEF_HDR.off_cal) * HEF_HDR.sf_cal
    except:
      print('hpproc.py: cannot decode cal data')
      print_info_and_exit()
      
    if options.swap_e1bc:
      x = -e1b
      e1b = -e1c
      e1c = x
      
    if options.swap_e1ab:
      x = e1b
      e1b = e1a
      e1a = x
      
    if options.swap_e2bc:
      x = -e2b
      e2b = -e2c
      e2c = x
      
    if options.swap_e2ab:
      x = e2b
      e2b = e2a
      e2a = x
      
    # file-length arrays
    if HEF_HDR.abu == 'a':
      cal_ind_a.append(ind)
      cal_e1a_a.append(e1a)
      cal_e1b_a.append(e1b)
      cal_e1c_a.append(e1c)
      cal_e2a_a.append(e2a)
      cal_e2b_a.append(e2b)
      cal_e2c_a.append(e2c)
    elif HEF_HDR.abu == 'b':
      cal_ind_b.append(ind)
      cal_e1a_b.append(e1a)
      cal_e1b_b.append(e1b)
      cal_e1c_b.append(e1c)
      cal_e2a_b.append(e2a)
      cal_e2b_b.append(e2b)
      cal_e2c_b.append(e2c)
    elif HEF_HDR.abu == 'u':
      cal_ind_u.append(ind)
      cal_e1a_u.append(e1a)
      cal_e1b_u.append(e1b)
      cal_e1c_u.append(e1c)
      cal_e2a_u.append(e2a)
      cal_e2b_u.append(e2b)
      cal_e2c_u.append(e2c)
    else:
      print('hpproc.py: unknown abu=',HEF_HDR.abu)
      print_info_and_exit()

    # half-cycle arrays of cal data
    if HEF_HDR.hef_iscan == 0:
      HCYC = collections.namedtuple('HCYC',[])
      HCYC.ind = []
      HCYC.e1a = []
      HCYC.e1b = []
      HCYC.e1c = []
      HCYC.e2a = []
      HCYC.e2b = []
      HCYC.e2c = []
    try:
      HCYC.ind.append(ind)
      HCYC.e1a.append(e1a)
      HCYC.e1b.append(e1b)
      HCYC.e1c.append(e1c)
      HCYC.e2a.append(e2a)
      HCYC.e2b.append(e2b)
      HCYC.e2c.append(e2c)
    except:
      print('hpproc.py: cannot append HCYC:')
      print('  hef_iscan=',HEF_HDR.hef_iscan)
      print('  lineno=',lineno)
      print('  linein=',linein)
      sys.exit(1)

    HEF_HDR.hef_iscan += 1

    if HEF_HDR.hef_iscan == HEF_HDR.hef_nscan:
      HEF_HDR.hef_iscan = -1

    if HEF_HDR.hef_iscan == -1:

#     print('hpproc.py: len(HCYC.ind)=',len(HCYC.ind))

      HCYC.abu = HEF_HDR.abu
      HCYC.uxt0 = HEF_HDR.uxt_ibeg
      HCYC.secs =  np.array(HCYC.ind,dtype='double') * META.tsamp_ef
      HCYC.uxt =  HCYC.secs + HCYC.uxt0
      HCYC.e1a = np.array(HCYC.e1a,dtype='double')
      HCYC.e1b = np.array(HCYC.e1b,dtype='double')
      HCYC.e1c = np.array(HCYC.e1c,dtype='double')
      HCYC.e2a = np.array(HCYC.e2a,dtype='double')
      HCYC.e2b = np.array(HCYC.e2b,dtype='double')
      HCYC.e2c = np.array(HCYC.e2c,dtype='double')


    if HEF_HDR.hef_iscan == -1 and options.do_plot_raw_cal:
      if options.limEC == 0 or count_EC <= options.limEC:
        plot_raw_cal_flag = True # defer plot_raw_cal() until get status

# defer compute averages and demod CAL data until get cal status
#   if HEF_HDR.hef_iscan == -1 and HEF_HDR.nskip_cal < HEF_HDR.hef_nscan:
#     print('hpproc.py: hef_iscan=',HEF_HDR.hef_iscan,'hef_nscan=',HEF_HDR.hef_nscan)
#     demod_cal_flag = True; # defer demod until receive cal status

    if HEF_HDR.hef_iscan > HEF_HDR.hef_nscan:
      print('hpproc.py: HEF_HDR.hef_iscan > hef_nscan')
      print_info_and_exit()

# else:
#   print('hpproc.py: unknown data, lineno=',lineno,'linein',linein)
#   print_info_and_exit()

def compute_cal_demod():
  global HCYC, CAL_DEMOD

  n = len(HCYC.uxt)
  if n <= 0:
    print('hpproc.py: len demod cal arrays == 0')
    return

  if len(CAL_STATUS.hcyc_status) != n:
    print('hpproc.py: demod_cal() len error')
    print_info_and_exit()

  # three basis functions
  squarewave = np.tile(0,(1,n))[0]
  squarewave[CAL_STATUS.neg] = -1.0
  squarewave[CAL_STATUS.pos] =  1.0
  trend = np.linspace(-1.0,1.0,num=n)
  constant = np.ones(n)

  if len(np.nonzero(squarewave == 0)[0]) > 0:
    print('hpproc.py: squarewave should have no zero value elements')
    print_info_and_exit()

  # find ju, indices to use
  # start with all indices available
  wt = np.ones(n)

  # toss indices before nskip_cal
  # to allow pinch tubes and preamps to stabilize
  j = range(HEF_HDR.nskip_cal)
  wt[j] = 0

  # toss nwait_cal indices after squarewave polarity change
  # to allow preamp to stabilize
  jflip = np.nonzero(np.diff(squarewave))[0]
  for iflip in jflip:
    j = np.arange(iflip,iflip+HEF_HDR.nwait_cal)
    i = np.nonzero(j < n)[0]
    if len(i) > 0:
      j = j[i]
    wt[j] = 0

  # use indices not tossed above
  ju = np.nonzero(wt != 0)[0]

  if len(ju) <= 3:
    print('hpproc: warning: ju=',ju,'len(wt)=',len(wt),\
      'HEF_HDR.nskip_cal=',HEF_HDR.nskip_cal,\
      'HEF_HDR.nwait_cal=',HEF_HDR.nwait_cal,
      'jflip=',jflip)
    sys.stdout.flush()
    
  HCYC.ju = ju
  HCYC.tu = HCYC.secs[ju]

  # find indices which occur during IES comms with HEF
  # HEF data xfr'd to IES 1:41 (101 s) after IES 10 minute mark
  uxtmj = np.mod(HCYC.uxt[j],600)
  HCYC.jjm = np.nonzero(np.logical_and(uxtmj > 98, uxtmj < 108))[0]

  HCYC.tujm = HCYC.tu[HCYC.jjm]
  HCYC.jjjm = ju[HCYC.jjm]

  # basis matrix for least squares fitting using np.linalg.lstsq()

  if len(ju)> 3:
    Bju = np.vstack([squarewave[ju], trend[ju], constant[ju]]).T
    # averages of synchronous demodulation
    HCYC.e1a_amp, e1a_trnd, e1a_con = np.linalg.lstsq(Bju, HCYC.e1a[ju])[0]
    HCYC.e1b_amp, e1b_trnd, e1b_con = np.linalg.lstsq(Bju, HCYC.e1b[ju])[0]
    HCYC.e1c_amp, e1c_trnd, e1c_con = np.linalg.lstsq(Bju, HCYC.e1c[ju])[0]
    HCYC.e2a_amp, e2a_trnd, e2a_con = np.linalg.lstsq(Bju, HCYC.e2a[ju])[0]
    HCYC.e2b_amp, e2b_trnd, e2b_con = np.linalg.lstsq(Bju, HCYC.e2b[ju])[0]
    HCYC.e2c_amp, e2c_trnd, e2c_con = np.linalg.lstsq(Bju, HCYC.e2c[ju])[0]

    # residuals to fits
    HCYC.e1a_res = HCYC.e1a - (HCYC.e1a_amp * squarewave + e1a_trnd * trend + e1a_con)
    HCYC.e1b_res = HCYC.e1b - (HCYC.e1b_amp * squarewave + e1b_trnd * trend + e1b_con)
    HCYC.e1c_res = HCYC.e1c - (HCYC.e1c_amp * squarewave + e1c_trnd * trend + e1c_con)
    HCYC.e2a_res = HCYC.e2a - (HCYC.e2a_amp * squarewave + e2a_trnd * trend + e2a_con)
    HCYC.e2b_res = HCYC.e2b - (HCYC.e2b_amp * squarewave + e2b_trnd * trend + e2b_con)
    HCYC.e2c_res = HCYC.e2c - (HCYC.e2c_amp * squarewave + e2c_trnd * trend + e2c_con)

    # compute standard deviations of residuals
    HCYC.e1a_std = np.std(HCYC.e1a_res[ju])
    HCYC.e1b_std = np.std(HCYC.e1b_res[ju])
    HCYC.e1c_std = np.std(HCYC.e1c_res[ju])
    HCYC.e2a_std = np.std(HCYC.e2a_res[ju])
    HCYC.e2b_std = np.std(HCYC.e2b_res[ju])
    HCYC.e2c_std = np.std(HCYC.e2c_res[ju])
  else:
    HCYC.e1a_amp = e1a_trnd = e1a_con = np.nan
    HCYC.e1b_amp = e1b_trnd = e1b_con = np.nan
    HCYC.e1c_amp = e1c_trnd = e1c_con = np.nan
    HCYC.e2a_amp = e2a_trnd = e2a_con = np.nan
    HCYC.e2b_amp = e2b_trnd = e2b_con = np.nan
    HCYC.e2c_amp = e2c_trnd = e2c_con = np.nan

    HCYC.e1a_res = np.nan
    HCYC.e1b_res = np.nan
    HCYC.e1c_res = np.nan
    HCYC.e2a_res = np.nan
    HCYC.e2b_res = np.nan
    HCYC.e2c_res = np.nan

    HCYC.e1a_std = np.nan
    HCYC.e1b_std = np.nan
    HCYC.e1c_std = np.nan
    HCYC.e2a_std = np.nan
    HCYC.e2b_std = np.nan
    HCYC.e2c_std = np.nan

  # append amplitudes and standard deviations for duration of file
  CAL_DEMOD.hcno.append(HEF_HDR.hcno)
  CAL_DEMOD.uxt.append(np.mean(HCYC.uxt))
  CAL_DEMOD.abu.append(HCYC.abu)
  CAL_DEMOD.navg.append(HEF_HDR.navg_cal)
  CAL_DEMOD.nuse.append(len(ju))
  CAL_DEMOD.e1a_amp.append(HCYC.e1a_amp)
  CAL_DEMOD.e1b_amp.append(HCYC.e1b_amp)
  CAL_DEMOD.e1c_amp.append(HCYC.e1c_amp)
  CAL_DEMOD.e2a_amp.append(HCYC.e2a_amp)
  CAL_DEMOD.e2b_amp.append(HCYC.e2b_amp)
  CAL_DEMOD.e2c_amp.append(HCYC.e2c_amp)
  CAL_DEMOD.e1a_std.append(HCYC.e1a_std)
  CAL_DEMOD.e1b_std.append(HCYC.e1b_std)
  CAL_DEMOD.e1c_std.append(HCYC.e1c_std)
  CAL_DEMOD.e2a_std.append(HCYC.e2a_std)
  CAL_DEMOD.e2b_std.append(HCYC.e2b_std)
  CAL_DEMOD.e2c_std.append(HCYC.e2c_std)

def plot_cal_demod():

  pdfdir = options.pdfroot + '/hef-avg/'
  mymkdir(pdfdir)

  ylim_big = [-800,800]
  ylim_small = [-5,5]

  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  abu = np.array(CAL_DEMOD.abu)
  ja = np.nonzero(abu == 'a')[0]
  jb = np.nonzero(abu == 'b')[0]

  uxt = np.array(CAL_DEMOD.uxt)
  hrsja = (uxt[ja] - uxt[0]) / 3600.0
  hrsjb = (uxt[jb] - uxt[0]) / 3600.0
  global hrsref
  hrsref = uxt[0] / 3600.0
  
  mltja = uxt[ja] / 86400 + 719529 - 366 
  mltjb = uxt[jb] / 86400 + 719529 - 366 


  e1a_amp_ja = np.array(CAL_DEMOD.e1a_amp)[ja]
  e1b_amp_ja = np.array(CAL_DEMOD.e1b_amp)[ja]
  e1c_amp_ja = np.array(CAL_DEMOD.e1c_amp)[ja]

  e1a_amp_jb = np.array(CAL_DEMOD.e1a_amp)[jb]
  e1b_amp_jb = np.array(CAL_DEMOD.e1b_amp)[jb]
  e1c_amp_jb = np.array(CAL_DEMOD.e1c_amp)[jb]

  e2a_amp_ja = np.array(CAL_DEMOD.e2a_amp)[ja]
  e2b_amp_ja = np.array(CAL_DEMOD.e2b_amp)[ja]
  e2c_amp_ja = np.array(CAL_DEMOD.e2c_amp)[ja]

  e2a_amp_jb = np.array(CAL_DEMOD.e2a_amp)[jb]
  e2b_amp_jb = np.array(CAL_DEMOD.e2b_amp)[jb]
  e2c_amp_jb = np.array(CAL_DEMOD.e2c_amp)[jb]

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-e1-cal-demod'.format(pltnambase)

  if options.interactive:
    fig.canvas.set_window_title('E1 Cal Demod')

  fig.suptitle('{0}\ntskip={1:.1f}, twait={2:.1f}, red:pinched, blu:unpinched\n{3:s}\n{4}'.
    format(pltnam,HEF_HDR.tskip_cal,HEF_HDR.twait_cal,pytstr,wsnams))

  axlist = []

  ax = fig.add_subplot(3,2,1)
  axlist.append(ax)
  ax.plot_date(mltja, e1a_amp_ja,'r.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e1ap, uV')
  plt.title('A pinched')

  ax = fig.add_subplot(3,2,2)
  axlist.append(ax)
  ax.plot_date(mltjb, e1a_amp_jb,'b.-')
  ax.set_ylim(ylim_small)
  ax.grid(True)
  plt.ylabel('e1au, uV')
  plt.title('B pinched')

  ax = fig.add_subplot(3,2,3)
  axlist.append(ax)
  ax.plot_date(mltja, e1b_amp_ja,'b.-')
  ax.set_ylim(ylim_small)
  ax.grid(True)
  plt.ylabel('e1bu, uV')

  ax = fig.add_subplot(3,2,4)
  axlist.append(ax)
  ax.plot_date(mltjb, e1b_amp_jb,'r.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e1bp, uV')

  ax = fig.add_subplot(3,2,5)
  axlist.append(ax)
  ax.plot_date(mltja, e1c_amp_ja,'g.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e1cap, uV')

  ax = fig.add_subplot(3,2,6)
  axlist.append(ax)
  ax.plot_date(mltjb, e1c_amp_jb,'g.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e1cbp, uV')

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  fix_xdates(axlist,2)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-e2-cal-demod'.format(pltnambase)

  if options.interactive:
    fig.canvas.set_window_title('E2 Cal Demod')

  fig.suptitle('{0}\ntskip={1:.1f}, twait={2:.1f}, red:pinched, blu:unpinched\n{3:s}\n{4}'.
    format(pltnam,HEF_HDR.tskip_cal,HEF_HDR.twait_cal,pytstr,wsnams))

  axlist = []

  ax = fig.add_subplot(3,2,1)
  axlist.append(ax)
  ax.plot_date(mltja, e2a_amp_ja,'r.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e2ap, uV')
  plt.title('A pinched')

  ax = fig.add_subplot(3,2,2)
  axlist.append(ax)
  ax.plot_date(mltjb, e2a_amp_jb,'b.-')
  ax.set_ylim(ylim_small)
  ax.grid(True)
  plt.ylabel('e2au, uV')
  plt.title('B pinched')

  ax = fig.add_subplot(3,2,3)
  axlist.append(ax)
  ax.plot_date(mltja, e2b_amp_ja,'b.-')
  ax.set_ylim(ylim_small)
  ax.grid(True)
  plt.ylabel('e2bu, uV')

  ax = fig.add_subplot(3,2,4)
  axlist.append(ax)
  ax.plot_date(mltjb, e2b_amp_jb,'r.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e2bp, uV')

  ax = fig.add_subplot(3,2,5)
  axlist.append(ax)
  ax.plot_date(mltja, e2c_amp_ja,'g.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e2cap, uV')
# plt.xlabel('UTC')

  ax = fig.add_subplot(3,2,6)
  axlist.append(ax)
  ax.plot_date(mltjb, e2c_amp_jb,'g.-')
  ax.set_ylim(ylim_big)
  ax.grid(True)
  plt.ylabel('e2cbp, uV')
# plt.xlabel('UTC')

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  fix_xdates(axlist,2)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

def fix_xdates(axlist,nss):
  fmt = mdates.DateFormatter('%b %d %H%MZ')
  for ax in axlist:
    xax = ax.get_xaxis()
    xax.set_major_formatter( fmt )

    if nss>1:
      xticks = ax.get_xticks()
      newxticks = []
      for i in range(0,len(xticks),nss):
        newxticks.append(xticks[i])
      ax.set_xticks(newxticks)

def accum_mot_cur_dur():
  global MOT_CUR

  if len(HCYM.uxt) == 0:
    return

  if len(MOT_STATUS.stat_pick) != len(HCYM.uxt):
    print('hpproc.py: accum_mot_cur_dur: lengths not right')
    print_info_and_exit(1)

  jf1 = np.nonzero(MOT_STATUS.stat_fast == 48)[0] # motor 1 running
  jf2 = np.nonzero(MOT_STATUS.stat_fast == 3 )[0] # motor 2 running

  jp1 = np.nonzero(MOT_STATUS.stat_pick == 48)[0] # motor 1 running
  jp2 = np.nonzero(MOT_STATUS.stat_pick == 3 )[0] # motor 2 running

  if len(jp1)>0 and len(jp2)>0:
    MOT_CUR.hcno.append(HEF_HDR.hcno)
    MOT_CUR.uxt1.append(HCYM.uxt[jp1[0]])
    MOT_CUR.uxt2.append(HCYM.uxt[jp2[0]])
    MOT_CUR.abu.append(HEF_HDR.abu)
    MOT_CUR.dur1.append(len(jf1) * META.tsamp_mot)
    MOT_CUR.dur2.append(len(jf2) * META.tsamp_mot)
    MOT_CUR.cur1.append(np.max(HCYM.cur[jp1]))
    MOT_CUR.cur2.append(np.max(HCYM.cur[jp2]))

def plot_mot_cur_dur():
  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  abu   = np.array(MOT_CUR.abu)
  uxt1  = np.array(MOT_CUR.uxt1,dtype='double')
  uxt2  = np.array(MOT_CUR.uxt2,dtype='double')
  cur1  = np.array(MOT_CUR.cur1,dtype='double')
  cur2  = np.array(MOT_CUR.cur2,dtype='double')
  dur1  = np.array(MOT_CUR.dur1,dtype='double')
  dur2  = np.array(MOT_CUR.dur2,dtype='double')

  ja = np.nonzero(abu == 'a')
  jb = np.nonzero(abu == 'b')

  uxtref = np.min([uxt1[0],uxt2[0]])
  hrs1  = (uxt1 - uxtref)/3600.0
  hrs2  = (uxt2 - uxtref)/3600.0

  mltoff = 719529 - 366 + uxtref / 86400
  mlt1 = hrs1 / 24.0 + mltoff
  mlt2 = hrs2 / 24.0 + mltoff

  pdfdir = options.pdfroot + '/hef-avg/'
  mymkdir(pdfdir)

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-mot-cur-dur'.format(pltnambase)

  fig.suptitle('{0}, red=A pinched, blu=B pinched\n{1}\n{2}'.
    format(pltnam,pytstr,wsnams))

  ylim_cur = [0,700]
  ylim_dur = [0,5.5]

  axlist = []

  ax = fig.add_subplot(2,2,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mlt1[ja], cur1[ja],'r.-')
  ax.plot_date(mlt1[jb], cur1[jb],'b.-')
  ax.hold(False)
  ax.set_ylim(ylim_cur)
  ax.grid(True)
  plt.ylabel('cur1')
  
  ax = fig.add_subplot(2,2,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mlt2[ja], cur2[ja],'r.-')
  ax.plot_date(mlt2[jb], cur2[jb],'b.-')
  ax.hold(False)
  ax.set_ylim(ylim_cur)
  ax.grid(True)
  plt.ylabel('cur2')
  
  ax = fig.add_subplot(2,2,3)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mlt1[ja], dur1[ja],'r.-')
  ax.plot_date(mlt1[jb], dur1[jb],'b.-')
  ax.hold(False)
  ax.set_ylim(ylim_dur)
  ax.grid(True)
  plt.ylabel('dur1')
# plt.xlabel('hours')
  
  ax = fig.add_subplot(2,2,4)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mlt2[ja], dur2[ja],'r.-')
  ax.plot_date(mlt2[jb], dur2[jb],'b.-')
  ax.hold(False)
  ax.set_ylim(ylim_dur)
  ax.grid(True)
  plt.ylabel('dur2')
# plt.xlabel('hours')

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  fix_xdates(axlist,2)

  if options.interactive:
    fig.show()
    print('hpproc.py: Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)


def check_hef_time():
  global  hef_secs_diff

  # header for EF, cal or motor
  if hef_split[0] == '#3__HE05':
    try:
      hef_secs_diff = int(hef_split[6]) - int(hef_split[11])
    except:
      print('hpproc.py: bad HE05 hef_secs_diff')
      print_info_and_exit()

def append_compass():
  global comp_uxt, comp_hdg, comp_pitch, comp_roll, comp_temp
  global count_compass

  if len(hef_split) < 1:
    return

  # compass data
  if hef_split[0].find('#3__HC') == 0:
    if hef_split[0] == '#3__HC03':
      if len(hef_split) != 8:
        print('hpproc.py: HC03 wrong split len:',len(hef_split))
        print('  hef_split:',hef_split)
        print_info_and_exit()
      count_compass += 1
    else:
      print('hpproc.py: unknown compass header')
      print_info_and_exit()

    try:
      uxt_HC  = int(hef_split[1])
      hdg   = int(hef_split[3])
      pitch = int(hef_split[4])
      roll  = int(hef_split[5])
      temp  = int(hef_split[6])
    except:
      print('hpproc.py: error in decoding compass data')
      print_info_and_exit()

    comp_uxt.append(uxt_HC)
    comp_hdg.append(hdg)
    comp_pitch.append(pitch)
    comp_roll.append(roll)
    comp_temp.append(temp)

def append_aux():
  global aux_uxt, aux_tt1, aux_tt2, aux_tt4, aux_tt4
  global aux_pres, aux_temp, aux_btemp, aux_bfreq
  global aux_uxt_xfr
  global uxt_aux, aux_flag

  if aux_split == None:
    return

  if len(aux_split) != 13:
    print('hpproc.py: len(aux_split) wrong: ',len(aux_split))
    print('aux_split:',aux_split)

  try:
    aux_flag = True
    uxt_aux = int(aux_split[1]) # IES time of beginning of 10-min sequence
    ntt = int(aux_split[2])
    tt1 = int(aux_split[3])
    tt2 = int(aux_split[4])
    tt3 = int(aux_split[5])
    tt4 = int(aux_split[6])
    pres = int(aux_split[7])
    temp = int(aux_split[8])
    btemp = int(aux_split[9])
    bfreq = float(aux_split[10])
    uxt_xfr = int(aux_split[12]) # STM time of data reception from IES's AUX
  except:
    print('hpproc.py: cannot decode aux port data')
    print_info_and_exit()

  if ntt != 4:
    print('hpproc.py: ntt=',ntt,'is wrong')
    print_info_and_exit()

  aux_uxt.append(uxt_aux)
  aux_tt1.append(tt1)
  aux_tt2.append(tt2)
  aux_tt3.append(tt3)
  aux_tt4.append(tt4)
  aux_pres.append(pres)
  aux_temp.append(temp)
  aux_btemp.append(btemp)
  aux_bfreq.append(bfreq)
  aux_uxt_xfr.append(uxt_xfr)

def plot_test():
  x = [1,2,3,4,5]
  y = [1,2,3,4,5]
  x  = np.array(x,dtype='double')
  y  = np.array(y,dtype='double')
  y = y * 5

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('Test')
  fig.suptitle('testing some matplotlib calls')

  axlist = []

  ax = fig.add_subplot(1,1,1)
  axlist.append(ax)
  ax.plot(x,y,'ro-')
  ax.grid(True)
  plt.xlabel('x')
  plt.title('plot test')
  plt.draw()

def compute_hef_avg():
  global HEF_AVG
  if len(HEF_AVG.uxt) == 0:
    init_HEF_AVG()
    return

  HEF_AVG.uxt = np.array(HEF_AVG.uxt)   # average unix time
  HEF_AVG.abu  = np.array(HEF_AVG.abu)    # A or B pinched

  HEF_AVG.e1am = np.array(HEF_AVG.e1am) # e1am in ../man/processing.txt
  HEF_AVG.e1bm = np.array(HEF_AVG.e1bm)
  HEF_AVG.e1cm = np.array(HEF_AVG.e1cm)
  HEF_AVG.e2am = np.array(HEF_AVG.e2am)
  HEF_AVG.e2bm = np.array(HEF_AVG.e2bm)
  HEF_AVG.e2cm = np.array(HEF_AVG.e2cm)

  HEF_AVG.e1as = np.array(HEF_AVG.e1as) # e1as in ../man/processing.txt
  HEF_AVG.e1bs = np.array(HEF_AVG.e1bs)
  HEF_AVG.e1cs = np.array(HEF_AVG.e1cs)
  HEF_AVG.e2as = np.array(HEF_AVG.e2as)
  HEF_AVG.e2bs = np.array(HEF_AVG.e2bs)
  HEF_AVG.e2cs = np.array(HEF_AVG.e2cs)

  HEF_AVG.e1afs = np.array(HEF_AVG.e1afs)  # fancy std dev
  HEF_AVG.e1bfs = np.array(HEF_AVG.e1bfs)
  HEF_AVG.e1cfs = np.array(HEF_AVG.e1cfs)
  HEF_AVG.e2afs = np.array(HEF_AVG.e2afs)
  HEF_AVG.e2bfs = np.array(HEF_AVG.e2bfs)
  HEF_AVG.e2cfs = np.array(HEF_AVG.e2cfs)

  HEF_AVG.t = HEF_AVG.uxt - HEF_AVG.uxt[0]
  HEF_AVG.ja = np.nonzero(HEF_AVG.abu == 'a')[0] # indices when A tube is pinched
  HEF_AVG.jb = np.nonzero(HEF_AVG.abu == 'b')[0] # indices when B tube is pinched

  t  = HEF_AVG.t
  ja = HEF_AVG.ja
  jb = HEF_AVG.jb

  # connect to ocean when pinch tube across preamp input
  HEF_AVG.e1apm = HEF_AVG.e1am[ja] # e1apm in ../man/processing.txt
  HEF_AVG.e1bpm = HEF_AVG.e1bm[jb]
  HEF_AVG.e2apm = HEF_AVG.e2am[ja]
  HEF_AVG.e2bpm = HEF_AVG.e2bm[jb]
  HEF_AVG.e1aps = HEF_AVG.e1as[ja] # e1aps in ../man/processing.txt
  HEF_AVG.e1bps = HEF_AVG.e1bs[jb]
  HEF_AVG.e2aps = HEF_AVG.e2as[ja]
  HEF_AVG.e2bps = HEF_AVG.e2bs[jb]
  HEF_AVG.e1apfs = HEF_AVG.e1afs[ja]
  HEF_AVG.e1bpfs = HEF_AVG.e1bfs[jb]
  HEF_AVG.e2apfs = HEF_AVG.e2afs[ja]
  HEF_AVG.e2bpfs = HEF_AVG.e2bfs[jb]

  # self potential when unpinch on same side as preamp
  HEF_AVG.e1aum = HEF_AVG.e1am[jb] # e1aum in ../man/processing.txt
  HEF_AVG.e1bum = HEF_AVG.e1bm[ja]
  HEF_AVG.e2aum = HEF_AVG.e2am[jb]
  HEF_AVG.e2bum = HEF_AVG.e2bm[ja]
  HEF_AVG.e1aus = HEF_AVG.e1as[jb] # e1aus in ../man/processing.txt
  HEF_AVG.e1bus = HEF_AVG.e1bs[ja]
  HEF_AVG.e2aus = HEF_AVG.e2as[jb]
  HEF_AVG.e2bus = HEF_AVG.e2bs[ja]
  HEF_AVG.e1aufs = HEF_AVG.e1afs[jb]
  HEF_AVG.e1bufs = HEF_AVG.e1bfs[ja]
  HEF_AVG.e2aufs = HEF_AVG.e2afs[jb]
  HEF_AVG.e2bufs = HEF_AVG.e2bfs[ja]

  # interpolate self potentials to ocean times
  HEF_AVG.e1aui = np.interp(t[ja],t[jb], HEF_AVG.e1aum)
  HEF_AVG.e1bui = np.interp(t[jb],t[ja], HEF_AVG.e1bum)
  HEF_AVG.e2aui = np.interp(t[ja],t[jb], HEF_AVG.e2aum)
  HEF_AVG.e2bui = np.interp(t[jb],t[ja], HEF_AVG.e2bum)

  # ocean estimate = pinch minus interpolated unpinch
  HEF_AVG.e1ao = HEF_AVG.e1apm - HEF_AVG.e1aui # t[ja] # e1ao in ../man/processing.txt
  HEF_AVG.e1bo = HEF_AVG.e1bpm - HEF_AVG.e1bui # t[jb]
  HEF_AVG.e2ao = HEF_AVG.e2apm - HEF_AVG.e2aui # t[ja]
  HEF_AVG.e2bo = HEF_AVG.e2bpm - HEF_AVG.e2bui # t[jb]

  HEF_AVG.e1o = np.tile(np.nan,(1,len(t)))[0]
  HEF_AVG.e2o = np.tile(np.nan,(1,len(t)))[0]
  HEF_AVG.e1o[ja] = HEF_AVG.e1ao
  HEF_AVG.e1o[jb] = HEF_AVG.e1bo
  HEF_AVG.e2o[ja] = HEF_AVG.e2ao
  HEF_AVG.e2o[jb] = HEF_AVG.e2bo

  # don't use last point -- interpolation may extrapolate poorly
  HEF_AVG.e1o[-1] = np.nan
  HEF_AVG.e1o[-1] = np.nan
  HEF_AVG.e2o[-1] = np.nan
  HEF_AVG.e2o[-1] = np.nan

  # new water switches have reversed sign
  # which is coded as the sign of esep
  HEF_AVG.E1 = -HEF_AVG.e1o / INFO.esep1      # uV/m
  HEF_AVG.E2 = -HEF_AVG.e2o / INFO.esep2      # uV/m
    
  # rotate to east (Ex) and north (Ey) components
  ang = COMP.hdg_avg + INFO.magvar - INFO.compass_angle
  c = np.cos(ang * np.pi / 180.0)
  s = np.sin(ang * np.pi / 180.0)
  HEF_AVG.Ex = c * HEF_AVG.E1 + s * HEF_AVG.E2      # uV/m
  HEF_AVG.Ey = c * HEF_AVG.E2 - s * HEF_AVG.E1      # uV/m

def plot_scatter():
  pdfdir = options.pdfroot + '/scatter/'
  mymkdir(pdfdir)

  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  m = 6
  n1a = int(len(HEF_AVG.e1ao)/m)
  n1b = int(len(HEF_AVG.e1bo)/m)
  n2a = int(len(HEF_AVG.e2ao)/m)
  n2b = int(len(HEF_AVG.e2bo)/m)

  if n1a != n2a:
    print('n1a & n2a should be the same')
    sys.exit(1)
  if n1b != n2b:
    print('n1b & n2b should be the same')
    sys.exit(1)

  na = n1a
  nb = n1b

  # make na == nb so can scatter-plot A against B
  if na > nb:
    na = nb
  if nb > na:
    nb = na

  # filter the ocean signals
  e1aof = np.mean(np.reshape(HEF_AVG.e1ao[0:m*na],(m,na),order='F'),0)
  e1bof = np.mean(np.reshape(HEF_AVG.e1bo[0:m*nb],(m,nb),order='F'),0)
  e2aof = np.mean(np.reshape(HEF_AVG.e2ao[0:m*na],(m,na),order='F'),0)
  e2bof = np.mean(np.reshape(HEF_AVG.e2bo[0:m*nb],(m,nb),order='F'),0)

  # filter the unpinched self potentials
  e1aumf = np.mean(np.reshape(HEF_AVG.e1aum[0:m*nb],(m,nb),order='F'),0)
  e1bumf = np.mean(np.reshape(HEF_AVG.e1bum[0:m*na],(m,na),order='F'),0)
  e2aumf = np.mean(np.reshape(HEF_AVG.e2aum[0:m*nb],(m,nb),order='F'),0)
  e2bumf = np.mean(np.reshape(HEF_AVG.e2bum[0:m*na],(m,na),order='F'),0)

  # scatter-plot filtered ocean data
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-ocean-scatter'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  ylim = [-options.ylim_ocean,options.ylim_ocean]
  xlim = [-options.ylim_ocean,options.ylim_ocean]

  axlist = []

  ax = fig.add_subplot(2,2,1,aspect='equal')
  axlist.append(ax)
  ax.plot(e1aof,e1bof,'b.')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.xlabel('e1a, uV')
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(2,2,2,aspect='equal')
  axlist.append(ax)
  ax.plot(e2aof,e2bof,'b.')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.xlabel('e2a, uV')
  plt.ylabel('e2b, uV')

  ax = fig.add_subplot(2,2,3,aspect='equal')
  axlist.append(ax)
  ax.plot(e1aof,e2aof,'b.')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.xlabel('e1a, uV')
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(2,2,4,aspect='equal')
  axlist.append(ax)
  ax.plot(e1bof,e2bof,'b.')
  ax.set_ylim(ylim)
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.xlabel('e1b, uV')
  plt.ylabel('e2b, uV')

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  if options.interactive:
    fig.show()
    print('hpproc.py: Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # scatter-plot filtered self-potential
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-self-pot-scatter'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  axlist = []

  ax = fig.add_subplot(2,2,1,aspect='equal')
  axlist.append(ax)
  ax.plot(e1aumf,e1bumf,'b.')
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e1a, uV')
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(2,2,2,aspect='equal')
  axlist.append(ax)
  ax.plot(e2aumf,e2bumf,'b.')
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e2a, uV')
  plt.ylabel('e2b, uV')

  ax = fig.add_subplot(2,2,3,aspect='equal')
  axlist.append(ax)
  ax.plot(e1aumf,e2aumf,'b.')
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e1a, uV')
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(2,2,4,aspect='equal')
  axlist.append(ax)
  ax.plot(e1bumf,e2bumf,'b.')
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e1b, uV')
  plt.ylabel('e2b, uV')

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  if options.interactive:
    fig.show()
    print('hpproc.py: Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # scatter-plot ocean data against self-potential
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-ocean-self-pot-scatter'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s\n{2}'.format(pltnam,HEF_HDR.tskip_ef,wsnams))

  axlist = []

  ax = fig.add_subplot(2,2,1)
  axlist.append(ax)
  ax.plot(e1aof,e1aumf,'b.')
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e1a ocean, uV')
  plt.ylabel('e1a self, uV')

  ax = fig.add_subplot(2,2,2)
  axlist.append(ax)
  ax.plot(e2aof,e2aumf,'b.')
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e2a ocean, uV')
  plt.ylabel('e2a self, uV')

  ax = fig.add_subplot(2,2,3)
  axlist.append(ax)
  ax.plot(e1bof,e1bumf,'b.')
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e1b ocean, uV')
  plt.ylabel('e1b self, uV')

  ax = fig.add_subplot(2,2,4)
  axlist.append(ax)
  ax.plot(e2bof,e2bumf,'b.')
  ax.set_xlim(xlim)
  ax.grid(True)
  plt.locator_params(nbins=4)
  plt.xlabel('e2b ocean, uV')
  plt.ylabel('e2b self, uV')

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  if options.interactive:
    fig.show()
    print('hpproc.py: Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

def plot_ef_avg_std():
  pdfdir = options.pdfroot + '/hef-avg/'
  mymkdir(pdfdir)
  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  ylim = [-options.ylim_ocean,options.ylim_ocean]
  mrkrsz1 = 10
  mrkrsz1 = 3
  mrkrsz2 = 5

  # plot ocean data
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-ocean'.format(pltnambase)

  if options.interactive:
    fig.canvas.set_window_title('HEF Ocean Averages')

  fig.suptitle('{0}\ntskip={1:.1f} s, red:A, blu:B\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  ja = HEF_AVG.ja
  jb = HEF_AVG.jb
  tc = HEF_AVG.t

# print('len(ja)=',len(ja))
# print('len(jb)=',len(jb))
# print('len(tc)=',len(tc))

# print('len(HEF_AVG.e1ao)=',len(HEF_AVG.e1ao))
# print('len(HEF_AVG.e1bo)=',len(HEF_AVG.e1bo))
# print('len(HEF_AVG.e2ao)=',len(HEF_AVG.e2ao))
# print('len(HEF_AVG.e2bo)=',len(HEF_AVG.e2bo))

  hrsja = tc[ja] / 3600
  hrsjb = tc[jb] / 3600
  hrsjc = tc    / 3600

# print('len(np.isfinite(hrsja))=',len(np.isfinite(hrsja)))
# print('len(np.isfinite(hrsjb))=',len(np.isfinite(hrsjb)))
# print('len(np.isfinite(hrsjc))=',len(np.isfinite(hrsjc)))

  m = 6
  na = int(len(HEF_AVG.e1ao)/m)
  nb = int(len(HEF_AVG.e1bo)/m)
  nc = int(len(tc)/m)

# print('na=',na,'nb=',nb,'nc=',nc)

  hrsjaf = np.mean(np.reshape(hrsja[0:m*na],(m,na),order='F'),0)
  hrsjbf = np.mean(np.reshape(hrsjb[0:m*nb],(m,nb),order='F'),0)
  hrsjcf = np.mean(np.reshape(hrsjc[0:m*nc],(m,nc),order='F'),0)

  mltoff = HEF_AVG.uxt[0] / 86400 + 719529 - 366
  mltja = hrsja / 24.0 + mltoff
  mltjb = hrsjb / 24.0 + mltoff
  mltjc = hrsjc / 24.0 + mltoff
  mltjaf = hrsjaf / 24.0 + mltoff
  mltjbf = hrsjbf / 24.0 + mltoff
  mltjcf = hrsjcf / 24.0 + mltoff

  HEF_AVG.e1aof = np.mean(np.reshape(HEF_AVG.e1ao[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e1bof = np.mean(np.reshape(HEF_AVG.e1bo[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2aof = np.mean(np.reshape(HEF_AVG.e2ao[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2bof = np.mean(np.reshape(HEF_AVG.e2bo[0:m*nb],(m,nb),order='F'),0)

  HEF_AVG.e1ausf = np.mean(np.reshape(HEF_AVG.e1aus[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e1busf = np.mean(np.reshape(HEF_AVG.e1bus[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2ausf = np.mean(np.reshape(HEF_AVG.e2aus[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2busf = np.mean(np.reshape(HEF_AVG.e2bus[0:m*na],(m,na),order='F'),0)

  HEF_AVG.e1aufsf = np.mean(np.reshape(HEF_AVG.e1aufs[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e1bufsf = np.mean(np.reshape(HEF_AVG.e1bufs[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2aufsf = np.mean(np.reshape(HEF_AVG.e2aufs[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2bufsf = np.mean(np.reshape(HEF_AVG.e2bufs[0:m*na],(m,na),order='F'),0)

  HEF_AVG.e1apsf = np.mean(np.reshape(HEF_AVG.e1aps[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e1bpsf = np.mean(np.reshape(HEF_AVG.e1bps[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2apsf = np.mean(np.reshape(HEF_AVG.e2aps[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2bpsf = np.mean(np.reshape(HEF_AVG.e2bps[0:m*nb],(m,nb),order='F'),0)

  HEF_AVG.e1apfsf = np.mean(np.reshape(HEF_AVG.e1apfs[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e1bpfsf = np.mean(np.reshape(HEF_AVG.e1bpfs[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2apfsf = np.mean(np.reshape(HEF_AVG.e2apfs[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2bpfsf = np.mean(np.reshape(HEF_AVG.e2bpfs[0:m*nb],(m,nb),order='F'),0)

  HEF_AVG.e1aumf = np.mean(np.reshape(HEF_AVG.e1aum[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e1bumf = np.mean(np.reshape(HEF_AVG.e1bum[0:m*na],(m,na),order='F'),0)
  HEF_AVG.e2aumf = np.mean(np.reshape(HEF_AVG.e2aum[0:m*nb],(m,nb),order='F'),0)
  HEF_AVG.e2bumf = np.mean(np.reshape(HEF_AVG.e2bum[0:m*na],(m,na),order='F'),0)

  HEF_AVG.e1cmf = np.mean(np.reshape(HEF_AVG.e1cm[0:m*nc],(m,nc),order='F'),0)
  HEF_AVG.e2cmf = np.mean(np.reshape(HEF_AVG.e2cm[0:m*nc],(m,nc),order='F'),0)

# print('len(np.isfinite(hrsjaf))=',len(np.isfinite(hrsjaf)))
# print('len(HEF_AVG.e1ao)=',len(HEF_AVG.e1ao))
# print('len(np.isfinite(HEF_AVG.e1ao))=',len(np.isfinite(HEF_AVG.e1ao)))
# print('len(np.isfinite(HEF_AVG.e1aof))=',len(np.isfinite(HEF_AVG.e1aof)))

# print('hrsja=',hrsja)
# print('hrsjaf=',hrsjaf)

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e1ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e1bo,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e1aof,'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e1bof,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  ax.grid(True)
  plt.ylabel(options.e1ws + ' e1o, uV')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e2ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e2bo,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e2aof,'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e2bof,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  ax.grid(True)
  plt.ylabel(options.e2ws + ' e2o, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  fix_xdates(axlist,1)
    
  if options.interactive:
    fig.show()
    print('hpproc.py: Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # plot standard deviations
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-std-dev'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s, red:pinched, blu:unpinched\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  ylim = [0,10]

  # print('tc[jb]=',tc[jb])

  axlist = []
  
  ax = fig.add_subplot(4,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e1aus,'b.',markersize=mrkrsz1)
  ax.plot_date(mltja,HEF_AVG.e1aps,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e1ausf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjaf,HEF_AVG.e1apsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e1a, uV')

  ax = fig.add_subplot(4,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e1bus,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e1bps,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e1busf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e1bpsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(4,1,3)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e2aus,'b.',markersize=mrkrsz1)
  ax.plot_date(mltja,HEF_AVG.e2aps,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e2ausf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjaf,HEF_AVG.e2apsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(4,1,4)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e2bus,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e2bps,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e2busf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e2bpsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e2b, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  fix_xdates(axlist,1)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # plot standard deviations of data minus polynomial fit
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-std-dev-poly'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s, red:pinched, blu:unpinched\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  ylim = [0,10]
  
  axlist = []

  ax = fig.add_subplot(4,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e1aufs,'b.',markersize=mrkrsz1)
  ax.plot_date(mltja,HEF_AVG.e1apfs,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e1aufsf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjaf,HEF_AVG.e1apfsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e1a, uV')

  ax = fig.add_subplot(4,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e1bufs,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e1bpfs,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e1bufsf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e1bpfsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(4,1,3)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e2aufs,'b.',markersize=mrkrsz1)
  ax.plot_date(mltja,HEF_AVG.e2apfs,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e2aufsf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjaf,HEF_AVG.e2apfsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(4,1,4)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e2bufs,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e2bpfs,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e2bufsf,'b.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e2bpfsf,'r.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim)
  fixlims(ax)
  ax.grid(True)
  plt.ylabel('e2b, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  fix_xdates(axlist,1)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # plot self potentials
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-self-pot'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))

  axlist = []
  
  ax = fig.add_subplot(6,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e1aum,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e1aumf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e1a')

  ax = fig.add_subplot(6,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e1bum,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e1bumf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e1b')

  ax = fig.add_subplot(6,1,3)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjc,HEF_AVG.e1cm,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjcf,HEF_AVG.e1cmf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e1c')
  plt.xlabel('hours, ' + pytstr)

  ax = fig.add_subplot(6,1,4)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e2aum,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e2aumf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e2a')

  ax = fig.add_subplot(6,1,5)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e2bum,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e2bumf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e2b')

  ax = fig.add_subplot(6,1,6)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjc,HEF_AVG.e2cm,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjcf,HEF_AVG.e2cmf,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e2c')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  fix_xdates(axlist,1)

  # plt.draw()

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # plot self potential differences 
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-self-pot-diffs'.format(pltnambase)

  fig.suptitle('{0}\ntskip={1:.1f} s\n{2}'.
    format(pltnam,HEF_HDR.tskip_ef,wsnams))
  
  HEF_AVG.e1aumi = np.interp(hrsja,hrsjb, HEF_AVG.e1aum)
  HEF_AVG.e1bumi = np.interp(hrsjb,hrsja, HEF_AVG.e1bum)
  HEF_AVG.e2aumi = np.interp(hrsja,hrsjb, HEF_AVG.e2aum)
  HEF_AVG.e2bumi = np.interp(hrsjb,hrsja, HEF_AVG.e2bum)

  HEF_AVG.e1aumfi = np.interp(hrsjaf,hrsjbf, HEF_AVG.e1aumf)
  HEF_AVG.e1bumfi = np.interp(hrsjbf,hrsjaf, HEF_AVG.e1bumf)
  HEF_AVG.e2aumfi = np.interp(hrsjaf,hrsjbf, HEF_AVG.e2aumf)
  HEF_AVG.e2bumfi = np.interp(hrsjbf,hrsjaf, HEF_AVG.e2bumf)

  ylimo = [-10,10]

  axlist = []

  ax = fig.add_subplot(4,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e1aum - HEF_AVG.e1bumi,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e1aumf - HEF_AVG.e1bumfi,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e1a-e1b self')

  ax = fig.add_subplot(4,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e1ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e1bo,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e1aof,'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e1bof,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylimo)
  ax.grid(True)
  plt.ylabel(options.e2ws + ' e1o, uV')

  ax = fig.add_subplot(4,1,3)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltjb,HEF_AVG.e2aum - HEF_AVG.e2bumi,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjbf,HEF_AVG.e2aumf - HEF_AVG.e2bumfi,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.yaxis.set_major_locator(MaxNLocator(4))
  ax.grid(True)
  plt.ylabel('e2a-e2b self')

  ax = fig.add_subplot(4,1,4)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,HEF_AVG.e2ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,HEF_AVG.e2bo,'b.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,HEF_AVG.e2aof,'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjbf,HEF_AVG.e2bof,'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylimo)
  ax.grid(True)
  plt.ylabel(options.e2ws + ' e2o, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  fix_xdates(axlist,1)

  # plt.draw()

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

def plot_tt():

  if len(aux_uxt) == 0:
    return

  uxt = np.array(aux_uxt,dtype='double')
  age = uxt - uxt[-1]
  tt1  = np.array(aux_tt1,dtype='double') * 1e-5
  tt2  = np.array(aux_tt2,dtype='double') * 1e-5
  tt3  = np.array(aux_tt3,dtype='double') * 1e-5
  tt4  = np.array(aux_tt4,dtype='double') * 1e-5

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('Travel Times')
  fig.suptitle(pltnambase + ' plot_tt')

  axlist = []

  ax = fig.add_subplot(4,1,1)
  axlist.append(ax)
  ax.plot(age,tt1,'bo-')
  ax.grid(True)
  plt.ylabel('TT1, s')
  plt.title('IES Travel Times')

  ax = fig.add_subplot(4,1,2)
  axlist.append(ax)
  ax.plot(age,tt2,'bo-')
  ax.grid(True)
  plt.ylabel('TT2, s')

  ax = fig.add_subplot(4,1,3)
  axlist.append(ax)
  ax.plot(age,tt3,'bo-')
  ax.grid(True)
  plt.ylabel('TT3, s')

  ax = fig.add_subplot(4,1,4)
  axlist.append(ax)
  ax.plot(age,tt4,'bo-')
  ax.grid(True)
  plt.ylabel('TT4, s')
  plt.xlabel('Time, s')

  plt.draw()

def plot_pres_temp():

  if len(aux_uxt) == 0:
    return

  uxt = np.array(aux_uxt,dtype='double')
  age = uxt - uxt[-1]
  pres  = np.array(aux_pres,dtype='double') * 1e-3 # fixed May 23, 2014
  temp  = np.array(aux_temp,dtype='double') * 1e-3

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('IES pressure & temperature')
  fig.suptitle(pltnambase + ' plot_pres_temp')

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.plot(age,pres,'bo-')
  ax.grid('on')
  plt.ylabel('P, dbar')
  plt.title('IES Pressure & Temperature')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.plot(age,temp,'bo-')
  ax.grid('on')
  plt.xlabel('Time, s')
  plt.ylabel('T, C')
  plt.draw()

def plot_bliley():

  if len(aux_uxt) == 0:
    return

  uxt = np.array(aux_uxt,dtype='double')
  age = uxt - uxt[-1]
  btemp  = np.array(aux_btemp,dtype='double') * 0.001
  bfreq  = np.array(aux_bfreq,dtype='double')

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('IES Bliley Oscillator')
  fig.suptitle(pltnambase + ' plot_bliley')

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.plot(age,btemp,'bo-')
  ax.grid('on')
  plt.ylabel('T, C')
  plt.title('Bliley Oscillator')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.plot(age,bfreq-4e6,'bo-')
  ax.grid('on')
  plt.ylabel('Freq - 4e6, Hz')
  plt.xlabel('Time, s')
  plt.draw()

def plot_mot_overlay():

  if len(mot_ind_a) == 0 and len(mot_ind_b) == 0:
    return

  secs_a = np.array(mot_ind_a,dtype='double') * META.tsamp_mot
  secs_b = np.array(mot_ind_b,dtype='double') * META.tsamp_mot
  secs_u = np.array(mot_ind_u,dtype='double') * META.tsamp_mot
  cur_a = np.array(mot_cur_a,dtype='double')
  cur_b = np.array(mot_cur_b,dtype='double')
  cur_u = np.array(mot_cur_u,dtype='double')


  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('HEF Motor')
  fig.suptitle(pltnambase + ' plot_mot_overlay')

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.plot(secs_a,cur_a,'b.-')
  ax.grid('on')
  plt.ylabel('A, mA')
  plt.title('Motor Current')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.plot(secs_b,cur_b,'b.-')
  ax.grid('on')
  plt.ylabel('B, mA')
  plt.xlabel('time, s ('+ifile+')')

  plt.draw()


# For OKMC, check that computed time of IES cycle start
# is same as the time from the AUX record
# The last column printed should be integer multiple of 600
def print_uxt_ies_chk():
  uxt_ies_fake = np.arange(0,240,1.024) + HCYE.uxt0
  tmod_fake = np.mod(uxt_ies_fake,600)
  j0_fake = np.nonzero(np.diff(tmod_fake) < 0)[0] + 1
  if len(j0_fake) == 1 and aux_flag:
    uxt_ies_j0 = uxt_ies_fake[j0_fake][0]
    print('hpproc.py: uxt_ies_j0={0:.1f}'.format(uxt_ies_j0),\
          'uxt_aux={0:.1f}'.format(uxt_aux),\
          'uxt_ies_j0-uxt_aux={0:5.1f}'.format(uxt_ies_j0-uxt_aux))

def plot_raw_mot():
# pyt = datetime(1970,1,1,0,0,0) + timedelta(0,HCYM.uxt0)
# pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'
  pltnam = str.format('{0}-{1:04d}-raw-mot',pltnambase,count_EC)

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('Raw Motor')
  fig.suptitle(pltnam + ' ' + pytstr + '\n' + wsnams)

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.plot(HCYM.secs,HCYM.cur,'b.-')
  fixlims(ax)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('MotCur, mA')
  plt.title('Raw Motor Current, tavg={0:.3f} s, red=A pinched, blu=B pinched'.format(HEF_HDR.navg_mot * META.tsamp_mot))

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  t = np.arange(0,len(MOT_STATUS.stat_fast)) * META.tsamp_mot
  ax.plot(t, MOT_STATUS.stat_fast,'b.-')
  fixlims(ax)
  ax.grid('on')
  plt.ylabel('stat_fast')
  plt.xlabel('time, s')

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = options.pdfroot + '/raw-mot/'+pltnam+'.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run pdf/raw-mot')

# plot concatinated motor current to see trends
def plot_raw_mot_cat():
  global NCYM
  global count_raw_mot_cat

  HCYM.ABU = HCYM.abu.upper()

  uxt = HCYM.uxt
  cur = HCYM.cur

  abi = np.tile(ord(HCYM.ABU),len(HCYM.uxt))
  sta = MOT_STATUS.stat_pick
  tpl = uxt - HCYM.uxt0

  count_raw_mot_cat += 1

  if len(uxt) < 2:
    print('hpproc.py: plot_raw_mot_cat() len(uxt) too short')

  if count_raw_mot_cat == 1:
    NCYM = collections.namedtuple('NCYM',[])
    NCYM.uxt0 = HCYM.uxt0
    NCYM.uxt_cat = uxt
    NCYM.cur_cat = cur

    NCYM.abi_cat = abi
    NCYM.sta_cat = sta
    NCYM.tpl_cat = tpl
  else:
    if len(NCYM.tpl_cat) > 0 and np.isfinite(NCYM.tpl_cat[-1]):
      tpl_off = NCYM.tpl_cat[-1] + 5
    else:
      tpl_off = 0

    if len(uxt) > 0:
      # insert nan to lift pen during motor runs
      NCYM.uxt_cat = np.append(NCYM.uxt_cat, np.nan)
      NCYM.cur_cat = np.append(NCYM.cur_cat, np.nan)

      NCYM.tpl_cat = np.append(NCYM.tpl_cat, np.nan)
      NCYM.abi_cat = np.append(NCYM.abi_cat, np.nan)
      NCYM.sta_cat = np.append(NCYM.sta_cat, np.nan)

      NCYM.uxt_cat = np.append(NCYM.uxt_cat, uxt)
      NCYM.cur_cat = np.append(NCYM.cur_cat, cur)

      NCYM.abi_cat = np.append(NCYM.abi_cat, abi)
      NCYM.sta_cat = np.append(NCYM.sta_cat, sta)
      NCYM.tpl_cat = np.append(NCYM.tpl_cat, tpl + tpl_off)

  if count_raw_mot_cat == options.ncat:
    count_raw_mot_cat = 0

#   pyt = datetime(1970,1,1,0,0,0) + timedelta(0,uxt0)
#   pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
    pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
    pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
    pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
             pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'
    pltnam = str.format('{0}-{1:04d}-raw-mot-cat{2}',pltnambase,count_EC,options.ncat)

    # zap cur_cat unless running
    jz = np.nonzero((sta_cat != 3) & (sta_cat != 48))
    NCYM.cur_cat[jz] = np.nan

    # include the nan values to lift pen between motor runs
    ja = np.nonzero((NCYM.abi_cat == ord('A')) | np.isnan(NCYM.abi_cat))[0]
    jb = np.nonzero((NCYM.abi_cat == ord('B')) | np.isnan(NCYM.abi_cat))[0]

    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    if options.interactive:
      fig.canvas.set_window_title('Motor Current')
    fig.suptitle(pltnam + ' ' + pytstr + '\n' + wsnams)

    axlist = []

    ax = fig.add_subplot(1,1,1)
    ax.hold(True)
    ax.plot(NCYM.tpl_cat[ja],NCYM.cur_cat[ja],'r.-')
    ax.plot(NCYM.tpl_cat[jb],NCYM.cur_cat[jb],'b.-')
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    # ax.xaxis.set_ticklabels([])
    plt.ylabel('mA')
    plt.title('Raw Motor Current, tavg={0:.3f} s, red=A pinched, blu=B pinched'.format(HEF_HDR.navg_mot * META.tsamp_mot))
    plt.xlabel('seconds')

    # put bigger gaps between plots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.4)

    if options.interactive:
      fig.show()
      print('Click figure window to continue')
      plt.waitforbuttonpress()
    else:
      pdfdir  = options.pdfroot + '/raw-mot-cat{0}'.format(options.ncat)
      mymkdir(pdfdir)
      pdffile = pdfdir + '/' + pltnam + '.pdf'
      print(pdffile)
      plt.savefig(pdffile)
      os.system('updateframe.run ' + pdfdir)

def compute_HCYE():
  global HCYE

  if HCYE.abu == 'a':
    HCYE.ma = 'r'
    HCYE.mb = 'b'
  elif HCYE.abu == 'b':
    HCYE.ma = 'b'
    HCYE.mb = 'r'
  else:
    HCYE.ma = 'g'
    HCYE.mb = 'g'

  # indices which lop off first nskip points
  j = np.arange(HEF_HDR.nskip_ef,len(HCYE.secs),dtype=int)

  # find indices which occur during IES comms with HEF
  # HEF data xfr'd to IES 1:41 (101 s) after IES 10 minute mark
  uxtmj = np.mod(HCYE.uxt[j],600)
  HCYE.jjm = np.nonzero(np.logical_and(uxtmj > 98, uxtmj < 108))[0]

  HCYE.tj = HCYE.secs[j]
  HCYE.tjjm = HCYE.tj[HCYE.jjm]

  HCYE.jjjm = j[HCYE.jjm]

  if len(HCYE.tj) < META.fit_order + 1:
    print('hpproc.py: len(HCYE.tj) too small, pltnam=',pltnam)
    return

  e1a_poly = np.polyfit(HCYE.tj,HCYE.e1a[j],META.fit_order)
  e1b_poly = np.polyfit(HCYE.tj,HCYE.e1b[j],META.fit_order)
  e1c_poly = np.polyfit(HCYE.tj,HCYE.e1c[j],META.fit_order)
  e2a_poly = np.polyfit(HCYE.tj,HCYE.e2a[j],META.fit_order)
  e2b_poly = np.polyfit(HCYE.tj,HCYE.e2b[j],META.fit_order)
  e2c_poly = np.polyfit(HCYE.tj,HCYE.e2c[j],META.fit_order)

  e1a_fit = np.polyval(e1a_poly,HCYE.tj)
  e1b_fit = np.polyval(e1b_poly,HCYE.tj)
  e1c_fit = np.polyval(e1c_poly,HCYE.tj)
  e2a_fit = np.polyval(e2a_poly,HCYE.tj)
  e2b_fit = np.polyval(e2b_poly,HCYE.tj)
  e2c_fit = np.polyval(e2c_poly,HCYE.tj)

  HCYE.e1a_res = HCYE.e1a[j] - e1a_fit
  HCYE.e1b_res = HCYE.e1b[j] - e1b_fit
  HCYE.e1c_res = HCYE.e1c[j] - e1c_fit
  HCYE.e2a_res = HCYE.e2a[j] - e2a_fit
  HCYE.e2b_res = HCYE.e2b[j] - e2b_fit
  HCYE.e2c_res = HCYE.e2c[j] - e2c_fit

  nbins = 50
  HCYE.e1a_hist = np.histogram(HCYE.e1a_res, bins=nbins)
  HCYE.e1b_hist = np.histogram(HCYE.e1b_res, bins=nbins)
  HCYE.e1c_hist = np.histogram(HCYE.e1c_res, bins=nbins)
  HCYE.e2a_hist = np.histogram(HCYE.e2a_res, bins=nbins)
  HCYE.e2b_hist = np.histogram(HCYE.e2b_res, bins=nbins)
  HCYE.e2c_hist = np.histogram(HCYE.e2c_res, bins=nbins)

def plot_raw_ef():
  pyt = datetime(1970,1,1,0,0,0) + timedelta(0,HCYE.uxt0)
  pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')

# pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
# pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
# pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
#          pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  pltnam = str.format('{0}-{1:04d}-raw-ef',pltnambase,count_EC)

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('HEF with IES pings')

  fig.suptitle(pltnam + '\nred:pinched, blu:unpinched\n' + pytstr + '\n' + wsnams)

  # indices which lop off first nskip points
  j = np.arange(HEF_HDR.nskip_ef,len(HCYE.secs),dtype=int)
  HCYE.tj = HCYE.secs[j]

  axlist = []

  ax = fig.add_subplot(6,2,1)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1a[j],HCYE.ma)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1a[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1a, uV')
  plt.title('raw EF, ' + HCYE.abu.upper() + ' pinched')

  ax = fig.add_subplot(6,2,3)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1b[j],HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1b[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(6,2,5)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1c[j],HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1c[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1c, uV')

  ax = fig.add_subplot(6,2,7)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2a[j],HCYE.ma)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2a[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(6,2,9)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2b[j],HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2b[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  plt.ylabel('e2b, uV')

  ax = fig.add_subplot(6,2,9)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2c[j],HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2c[HCYE.jjjm],'go')
  ax.hold(False)
  ax.grid('on')
  plt.ylabel('e2c, uV')
  plt.xlabel('time, s')

  ax = fig.add_subplot(6,2,2)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1a_res,HCYE.ma)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1a_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1a-fit, uV')
  plt.title('raw EF minus fit(order={0})'.format(META.fit_order))

  ax = fig.add_subplot(6,2,4)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1b_res,HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1b_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1b-fit, uV')

  ax = fig.add_subplot(6,2,6)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e1c_res,HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e1c_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1c-fit, uV')

  ax = fig.add_subplot(6,2,8)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2a_res,HCYE.ma)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2a_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2a-fit, uV')

  ax = fig.add_subplot(6,2,10)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2b_res,HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2b_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  plt.ylabel('e2b-fit, uV')

  ax = fig.add_subplot(6,2,12)
  axlist.append(ax)
  ax.plot(HCYE.tj,HCYE.e2c_res,HCYE.mb)
  ax.hold(True)
  ax.plot(HCYE.tjjm,HCYE.e2c_res[HCYE.jjm],'go')
  ax.hold(False)
  ax.grid('on')
  plt.ylabel('e2c-fit, uV')
  plt.xlabel('time, s')

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  # plt.draw()

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdffile = options.pdfroot + '/raw-ef/'+pltnam+'.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run pdf/raw-ef')

# plot concatinated EF data to see trends
def plot_raw_ef_cat():
  global NCYE, count_raw_ef_cat


  # indices which lop off first nskip_ef points
  j = np.arange(HEF_HDR.nskip_ef,len(HCYE.secs),dtype=int)
  uxt = HCYE.uxt[j]
  e1a = HCYE.e1a[j]
  e1b = HCYE.e1b[j]
  e1c = HCYE.e1c[j]
  e2a = HCYE.e2a[j]
  e2b = HCYE.e2b[j]
  e2c = HCYE.e2c[j]

  abi = np.tile(ord(HCYE.abu),len(uxt))

  # print('count_EC=',count_EC)

  count_raw_ef_cat += 1

# if len(uxt) < 2:
#   print('uxt too short: count_raw_ef_cat=',count_raw_ef_cat)

  if count_raw_ef_cat == 1:
    NCYE = collections.namedtuple('NCYE',[])
    NCYE.uxt0 = HCYE.uxt0

    NCYE.uxt_cat = uxt
    NCYE.e1a_cat = e1a
    NCYE.e1b_cat = e1b
    NCYE.e1c_cat = e1c
    NCYE.e2a_cat = e2a
    NCYE.e2b_cat = e2b
    NCYE.e2c_cat = e2c

    NCYE.abi_cat = abi
  else:
    if len(uxt) > 0:
      # insert nan to lift pen during motor runs
      NCYE.uxt_cat = np.append(NCYE.uxt_cat, np.nan)
      NCYE.e1a_cat = np.append(NCYE.e1a_cat, np.nan)
      NCYE.e1b_cat = np.append(NCYE.e1b_cat, np.nan)
      NCYE.e1c_cat = np.append(NCYE.e1c_cat, np.nan)
      NCYE.e2a_cat = np.append(NCYE.e2a_cat, np.nan)
      NCYE.e2b_cat = np.append(NCYE.e2b_cat, np.nan)
      NCYE.e2c_cat = np.append(NCYE.e2c_cat, np.nan)

      NCYE.abi_cat = np.append(NCYE.abi_cat, np.nan)

      NCYE.uxt_cat = np.append(NCYE.uxt_cat, uxt)
      NCYE.e1a_cat = np.append(NCYE.e1a_cat, e1a)
      NCYE.e1b_cat = np.append(NCYE.e1b_cat, e1b)
      NCYE.e1c_cat = np.append(NCYE.e1c_cat, e1c)
      NCYE.e2a_cat = np.append(NCYE.e2a_cat, e2a)
      NCYE.e2b_cat = np.append(NCYE.e2b_cat, e2b)
      NCYE.e2c_cat = np.append(NCYE.e2c_cat, e2c)

      NCYE.abi_cat = np.append(NCYE.abi_cat, abi)

  if count_raw_ef_cat == options.ncat:
    count_raw_ef_cat = 0

    # find indices which occur during IES comms with HEF
    # HEF data xfr'd to IES 1:41 (101 s) after IES 10 minute mark
    tmod = np.mod(NCYE.uxt_cat,600)
    jies = np.nonzero(np.logical_and(tmod > 98, tmod < 108))[0]

    j = np.nonzero(np.isfinite(NCYE.uxt_cat))[0]
    ja = np.nonzero(NCYE.abi_cat == ord('a'))[0]
    jb = np.nonzero(NCYE.abi_cat == ord('b'))[0]
    if len(j) != len(ja) + len(jb):
      print('hpproc.py: plot_raw_cal_cat(): len(j) should equal len(ja) + len(jb)')
      print_info_and_exit()

    tj_cat = NCYE.uxt_cat - NCYE.uxt0

#   pyt = datetime(1970,1,1,0,0,0) + timedelta(0,NCYE.uxt0)
#   pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
    pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,NCYE.uxt_cat[0])
    pytend = datetime(1970,1,1,0,0,0) + timedelta(0,NCYE.uxt_cat[-1])
    pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
             pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'
    pltnam = str.format('{0}-{1:04d}-raw-ef-cat{2}',pltnambase,count_EC,options.ncat)

    if (len(jb) < META.fit_order + 1) | (len(ja) < META.fit_order + 1):
      print('hpproc.py: too few ja or jb, pltnam=',pltnam)
      return

    e1a_poly = np.polyfit(tj_cat[jb],NCYE.e1a_cat[jb],META.fit_order)
    e1b_poly = np.polyfit(tj_cat[ja],NCYE.e1b_cat[ja],META.fit_order)
    e1c_poly = np.polyfit(tj_cat[ja],NCYE.e1c_cat[ja],META.fit_order)
    e2a_poly = np.polyfit(tj_cat[jb],NCYE.e2a_cat[jb],META.fit_order)
    e2b_poly = np.polyfit(tj_cat[ja],NCYE.e2b_cat[ja],META.fit_order)
    e2c_poly = np.polyfit(tj_cat[ja],NCYE.e2c_cat[ja],META.fit_order)
    e1a_fit = np.tile(np.nan,len(tj_cat))
    e1b_fit = np.tile(np.nan,len(tj_cat))
    e1c_fit = np.tile(np.nan,len(tj_cat))
    e2a_fit = np.tile(np.nan,len(tj_cat))
    e2b_fit = np.tile(np.nan,len(tj_cat))
    e2c_fit = np.tile(np.nan,len(tj_cat))
    e1a_fit[j] = np.polyval(e1a_poly,tj_cat[j])
    e1b_fit[j] = np.polyval(e1b_poly,tj_cat[j])
    e1c_fit[j] = np.polyval(e1c_poly,tj_cat[j])
    e2a_fit[j] = np.polyval(e2a_poly,tj_cat[j])
    e2b_fit[j] = np.polyval(e2b_poly,tj_cat[j])
    e2c_fit[j] = np.polyval(e2c_poly,tj_cat[j])
    NCYE.e1a_res_cat = NCYE.e1a_cat - e1a_fit
    NCYE.e1b_res_cat = NCYE.e1b_cat - e1b_fit
    NCYE.e1c_res_cat = NCYE.e1c_cat - e1c_fit
    NCYE.e2a_res_cat = NCYE.e2a_cat - e2a_fit
    NCYE.e2b_res_cat = NCYE.e2b_cat - e2b_fit
    NCYE.e2c_res_cat = NCYE.e2c_cat - e2c_fit

    ja = np.nonzero((NCYE.abi_cat == ord('a')) | np.isnan(NCYE.abi_cat))[0]
    jb = np.nonzero((NCYE.abi_cat == ord('b')) | np.isnan(NCYE.abi_cat))[0]

    mrkp = 'r'
    mrku = 'c'
    mrki = 'k.'
    mrkc = 'g'
    
    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    if options.interactive:
      fig.canvas.set_window_title('HEF with IES pings')
    fig.suptitle(pltnam + '\nred:pinched, cyan:unpinched, blk:IES comms\n' + pytstr + '\n' + wsnams)

    axlist = []

    ax = fig.add_subplot(6,2,1)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1a_cat[jies],mrki)
    ax.plot(tj_cat[ja],  NCYE.e1a_cat[ja],mrkp)
    ax.plot(tj_cat[jb],  NCYE.e1a_cat[jb],mrku)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a, uV')
    plt.title('raw EF cat, tskip={0:.1f} s'.format(options.tskip))

    ax = fig.add_subplot(6,2,3)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1b_cat[jies],mrki)
    ax.plot(tj_cat[ja],  NCYE.e1b_cat[ja],mrku)
    ax.plot(tj_cat[jb],  NCYE.e1b_cat[jb],mrkp)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b, uV')

    ax = fig.add_subplot(6,2,5)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1c_cat[jies],mrki)
    ax.plot(tj_cat[ja],  NCYE.e1c_cat[ja],mrkc)
    ax.plot(tj_cat[jb],  NCYE.e1c_cat[jb],mrkc)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1c, uV')

    ax = fig.add_subplot(6,2,7)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2a_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2a_cat[ja],mrkp)
    ax.plot(tj_cat[jb],NCYE.e2a_cat[jb],mrku)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a, uV')

    ax = fig.add_subplot(6,2,9)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2b_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2b_cat[ja],mrku)
    ax.plot(tj_cat[jb],NCYE.e2b_cat[jb],mrkp)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    plt.ylabel('e2b, uV')
    plt.xlabel('time, s')

    ax = fig.add_subplot(6,2,11)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2c_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2c_cat[ja],mrkc)
    ax.plot(tj_cat[jb],NCYE.e2c_cat[jb],mrkc)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    plt.ylabel('e2c, uV')
    plt.xlabel('time, s')

    ax = fig.add_subplot(6,2,2)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1a_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e1a_res_cat[ja],mrkp)
    ax.plot(tj_cat[jb],NCYE.e1a_res_cat[jb],mrku)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a-fit, uV')
    plt.title('raw EF cat minus fit(order={0})'.format(META.fit_order))
    ylim2 = ax.get_ylim()
    ax2 = ax

    ax = fig.add_subplot(6,2,4)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1b_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e1b_res_cat[ja],mrku)
    ax.plot(tj_cat[jb],NCYE.e1b_res_cat[jb],mrkp)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b-fit, uV')
    ylim4 = ax.get_ylim()
    ax4 = ax

    ax = fig.add_subplot(6,2,6)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e1c_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e1c_res_cat[ja],mrkc)
    ax.plot(tj_cat[jb],NCYE.e1c_res_cat[jb],mrkc)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1c-fit, uV')
    ylim6 = ax.get_ylim()
    ax6 = ax

    ylim = [min(ylim2[0],ylim4[0]),max(ylim2[1],ylim4[1])]
    ax2.set_ylim(ylim)
    ax4.set_ylim(ylim)

    ax = fig.add_subplot(6,2,8)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2a_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2a_res_cat[ja],mrkp)
    ax.plot(tj_cat[jb],NCYE.e2a_res_cat[jb],mrku)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a-fit, uV')
    ylim8 = ax.get_ylim()
    ax8 = ax

    ax = fig.add_subplot(6,2,10)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2b_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2b_res_cat[ja],mrku)
    ax.plot(tj_cat[jb],NCYE.e2b_res_cat[jb],mrkp)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    plt.ylabel('e2b-fit, uV')
    plt.xlabel('time, s')
    ylim10 = ax.get_ylim()
    ax10 = ax

    ax = fig.add_subplot(6,2,12)
    axlist.append(ax)
    ax.hold(True)
    ax.plot(tj_cat[jies],NCYE.e2c_res_cat[jies],mrki)
    ax.plot(tj_cat[ja],NCYE.e2c_res_cat[ja],mrkc)
    ax.plot(tj_cat[jb],NCYE.e2c_res_cat[jb],mrkc)
    ax.hold(False)
    fixlims(ax)
    ax.grid('on')
    plt.ylabel('e2c-fit, uV')
    plt.xlabel('time, s')
    ylim12 = ax.get_ylim()
    ax12 = ax

    ylim = [min(ylim8[0],ylim12[0]),max(ylim8[1],ylim12[1])]
    ax8.set_ylim(ylim)
    ax12.set_ylim(ylim)

#   for ax in [ax2,ax4,ax6,ax8]:
    for ax in axlist:
      xlim = ax.get_xlim()
      ylim = ax.get_ylim()
      i0 = int((xlim[0]+NCYE.uxt0-101)/600)+1
      i1 = int((xlim[1]+NCYE.uxt0-101)/600)+1
      ax.hold(True)
      for i in range(i0,i1):
        tplt = i * 600 - NCYE.uxt0 + 101
        ax.plot(tplt,ylim[0],'kd')
        ax.plot(tplt,ylim[1],'kd')
      ax.hold(False)
      fixlims(ax)

    # put bigger gaps between plots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

    if options.interactive:
      fig.show()
      print('Click figure window to continue')
      plt.waitforbuttonpress()
    else:
      pdfdir  = options.pdfroot + '/raw-ef-cat{0}'.format(options.ncat)
      mymkdir(pdfdir)
      pdffile = pdfdir + '/' + pltnam + '.pdf'
      print(pdffile)
      plt.savefig(pdffile)
      os.system('updateframe.run ' + pdfdir)

def plot_raw_cal():

  if HCYC.abu == 'a':
    HCYC.ma = 'r'
    HCYC.mb = 'b'
  elif HCYC.abu == 'b':
    HCYC.ma = 'b'
    HCYC.mb = 'r'
  else:
    HCYC.ma = 'g'
    HCYC.mb = 'g'

  pyt = datetime(1970,1,1,0,0,0) + timedelta(0,HCYC.uxt0)
  pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
# pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
# pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
# pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
#          pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  pltnam = str.format('{0}-{1:04d}-raw-cal',pltnambase,count_EC)

  # j = np.arange(HEF_HDR.nskip_cal,len(HCYC.secs),dtype=int)
  ju = HCYC.ju

  if HCYC.tu == None:
    print('hpproc.py: HCYC.tu == None, pltnam=',pltnam)
    return

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('CAL with IES pings')
  fig.suptitle(pltnam + '\nred:pinched, blu:unpinched\n' + pytstr + '\n' + wsnams)

  print('hpproc.py: len(ju)=',len(ju),'len(HCYC.tu)=',len(HCYC.tu),'len(HCYC.e1a)=',len(HCYC.e1a))

  axlist = []

  ax = fig.add_subplot(7,2,1)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e1a[ju],HCYC.ma)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1a[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1a, uV')
  plt.title('raw CAL, ' + HCYC.abu.upper() + ' pinched')

  ax = fig.add_subplot(7,2,3)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e1b[ju],HCYC.mb)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1b[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(7,2,5)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e1c[ju],'k')
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1c[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1c, uV')

  ax = fig.add_subplot(7,2,7)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2a[ju],HCYC.ma)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2a[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(7,2,9)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2b[ju],HCYC.mb)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2b[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2b, uV')

  ax = fig.add_subplot(7,2,11)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2c[ju],'k')
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2c[HCYC.jjjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  plt.ylabel('e2c, uV')
  plt.xlabel('time, s')

  ax = fig.add_subplot(7,2,2)
  axlist.append(ax)
  print('len(HCYC.tu)=',len(HCYC.tu),'len(HCYC.e1a_res)=',len(HCYC.e1a_res))
  ax.plot(HCYC.tu,HCYC.e1a_res[ju],HCYC.ma)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1a_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1a-fit, uV')
  plt.title('raw CAL minus fit(order={0})'.format(META.fit_order))

  ax = fig.add_subplot(7,2,4)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e1b_res[ju],HCYC.mb)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1b_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1b-fit, uV')

  ax = fig.add_subplot(7,2,6)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e1c_res[ju],'k')
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e1c_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e1c-fit, uV')

  ax = fig.add_subplot(7,2,8)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2a_res[ju],HCYC.ma)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2a_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2a-fit, uV')

  ax = fig.add_subplot(7,2,10)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2b_res[ju],HCYC.mb)
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2b_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  ax.xaxis.set_ticklabels([])
  plt.ylabel('e2b-fit, uV')

  ax = fig.add_subplot(7,2,12)
  axlist.append(ax)
  ax.plot(HCYC.tu,HCYC.e2c_res[ju],'k')
  ax.hold(True)
  ax.plot(HCYC.tujm,HCYC.e2c_res[HCYC.jjm],'go')
  ax.hold(False)
  fixlims(ax)
  ax.yaxis.set_major_locator(MaxNLocator(5))
  ax.grid('on')
  plt.ylabel('e2c-fit, uV')
  plt.xlabel('time, s')

  t = np.arange(0,len(CAL_STATUS.stat_fast)) * META.tsamp_ef
  ax.plot(t, CAL_STATUS.stat_fast,'b')
  ax.grid('on')
  plt.ylabel('stat_fast')
  plt.xlabel('time, s, ' + pytstr)
  fixlims(ax)

  ax = fig.add_subplot(7,2,14)
  axlist.append(ax)
  ax.plot(HCYC.secs, CAL_STATUS.hcyc_status,'b')
  ax.grid('on')
  plt.ylabel('hcyc_status')
  plt.xlabel('time, s, ' + pytstr)
  fixlims(ax)

  # put bigger gaps between plots
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    pdfdir = options.pdfroot + '/raw-cal'
    mymkdir(pdfdir)
    pdffile = pdfdir + '/' + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

def plot_raw_cal_cat():
  global NCYC
  global count_raw_cal_cat

  # indices which lop off first nskip_cal points
  j = np.arange(HEF_HDR.nskip_cal,len(HCYC.secs),dtype=int)
  uxt = HCYC.uxt[j]
  e1a = HCYC.e1a[j]
  e1b = HCYC.e1b[j]
  e1c = HCYC.e1c[j]
  e2a = HCYC.e2a[j]
  e2b = HCYC.e2b[j]
  e2c = HCYC.e2c[j]


  abi = np.tile(ord(HCYC.abu),len(j))
  sta = CAL_STATUS.hcyc_status[j]
  tpl = uxt - HCYC.uxt0
  stf = CAL_STATUS.stat_fast
  tpf = np.arange(0,len(stf)) * META.tsamp_ef


  count_raw_cal_cat += 1

# if len(uxt) < 2:
#   print('uxt too short: count_raw_cal_cat=',count_raw_cal_cat)
#   print('  nskip_cal=',HEF_HDR.nskip_cal)
#   print('  len(HCYC.secs)=',len(HCYC.secs))
#   print('  HCYC.secs=',HCYC.secs)
#   print('  len(j)=',len(j))
#   print('  tsamp_ef=',META.tsamp_ef)
#   print('  tskip_cal=',tskip_cal)
#   print('  navg_cal=',HEF_HDR.navg_cal)

  if count_raw_cal_cat == 1:
    NCYC = collections.namedtuple('NCYC',[])
    NCYC.uxt0 = HCYC.uxt0

    NCYC.uxt_cat = uxt
    NCYC.e1a_cat = e1a
    NCYC.e1b_cat = e1b
    NCYC.e1c_cat = e1c
    NCYC.e2a_cat = e2a
    NCYC.e2b_cat = e2b
    NCYC.e2c_cat = e2c

    NCYC.abi_cat = abi
    NCYC.sta_cat = sta
    NCYC.stf_cat = stf
    NCYC.tpl_cat = tpl
    NCYC.tpf_cat = tpf
  else:
    if len(NCYC.tpl_cat) > 0 and np.isfinite(NCYC.tpl_cat[-1]):
      tpl_off = NCYC.tpl_cat[-1] + 5
      tpf_off = NCYC.tpf_cat[-1] + 5
    else:
      tpl_off = 0
      tpf_off = 0

    if len(uxt) > 0:
      # insert nan to lift pen during motor runs
      NCYC.uxt_cat = np.append(NCYC.uxt_cat, np.nan)
      NCYC.e1a_cat = np.append(NCYC.e1a_cat, np.nan)
      NCYC.e1b_cat = np.append(NCYC.e1b_cat, np.nan)
      NCYC.e1c_cat = np.append(NCYC.e1c_cat, np.nan)
      NCYC.e2a_cat = np.append(NCYC.e2a_cat, np.nan)
      NCYC.e2b_cat = np.append(NCYC.e2b_cat, np.nan)
      NCYC.e2c_cat = np.append(NCYC.e2c_cat, np.nan)

      NCYC.abi_cat = np.append(NCYC.abi_cat, np.nan)
      NCYC.sta_cat = np.append(NCYC.sta_cat, np.nan)
      NCYC.stf_cat = np.append(NCYC.stf_cat, np.nan)
      NCYC.tpl_cat = np.append(NCYC.tpl_cat, np.nan)
      NCYC.tpf_cat = np.append(NCYC.tpf_cat, np.nan)

      NCYC.uxt_cat = np.append(NCYC.uxt_cat, uxt)
      NCYC.e1a_cat = np.append(NCYC.e1a_cat, e1a)
      NCYC.e1b_cat = np.append(NCYC.e1b_cat, e1b)
      NCYC.e1c_cat = np.append(NCYC.e1c_cat, e1c)
      NCYC.e2a_cat = np.append(NCYC.e2a_cat, e2a)
      NCYC.e2b_cat = np.append(NCYC.e2b_cat, e2b)
      NCYC.e2c_cat = np.append(NCYC.e2c_cat, e2c)

      NCYC.abi_cat = np.append(NCYC.abi_cat, abi)
      NCYC.sta_cat = np.append(NCYC.sta_cat, sta)
      NCYC.stf_cat = np.append(NCYC.stf_cat, stf)
      NCYC.tpl_cat = np.append(NCYC.tpl_cat, tpl + tpl_off)
      NCYC.tpf_cat = np.append(NCYC.tpf_cat, tpf + tpf_off)

  if count_raw_cal_cat == options.ncat:
    count_raw_cal_cat = 0

    # find indices which occur during IES comms with HEF
    # HEF data xfr'd to IES 1:41 (101 s) after IES 10 minute mark
    tmod = np.mod(NCYC.uxt_cat,600)
    jies = np.nonzero(np.logical_and(tmod > 98, tmod < 108))[0]

    j = np.nonzero(np.isfinite(NCYC.uxt_cat))[0]
    ja = np.nonzero(NCYC.abi_cat == ord('a'))[0]
    jb = np.nonzero(NCYC.abi_cat == ord('b'))[0]
    if len(j) != len(ja) + len(jb):
      print('hpproc.py: plot_raw_cal_cat(): len(j) should equal len(ja) + len(jb)')
      print_info_and_exit()

    tj_cat = NCYC.uxt_cat - NCYC.uxt0

#   pyt = datetime(1970,1,1,0,0,0) + timedelta(0,NCYC.uxt0)
#   pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
    pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,NCYC.uxt_cat[0])
    pytend = datetime(1970,1,1,0,0,0) + timedelta(0,NCYC.uxt_cat[-1])
    pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
             pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'
    pltnam = str.format('{0}-{1:04d}-raw-cal-cat{2}',pltnambase,count_EC,options.ncat)

    if (len(jb) < META.fit_order + 1) | (len(ja) < META.fit_order + 1):
      print('hpproc.py: too few ja or jb, pltnam=',pltnam)
      return

    e1a_poly = np.polyfit(tj_cat[jb],NCYC.e1a_cat[jb],META.fit_order)
    e1b_poly = np.polyfit(tj_cat[ja],NCYC.e1b_cat[ja],META.fit_order)
    e1c_poly = np.polyfit(tj_cat[ja],NCYC.e1c_cat[ja],META.fit_order)
    e2a_poly = np.polyfit(tj_cat[jb],NCYC.e2a_cat[jb],META.fit_order)
    e2b_poly = np.polyfit(tj_cat[ja],NCYC.e2b_cat[ja],META.fit_order)
    e2c_poly = np.polyfit(tj_cat[ja],NCYC.e2c_cat[ja],META.fit_order)
    e1a_fit = np.tile(np.nan,len(tj_cat))
    e1b_fit = np.tile(np.nan,len(tj_cat))
    e1c_fit = np.tile(np.nan,len(tj_cat))
    e2a_fit = np.tile(np.nan,len(tj_cat))
    e2b_fit = np.tile(np.nan,len(tj_cat))
    e2c_fit = np.tile(np.nan,len(tj_cat))
    e1a_fit[j] = np.polyval(e1a_poly,tj_cat[j])
    e1b_fit[j] = np.polyval(e1b_poly,tj_cat[j])
    e1c_fit[j] = np.polyval(e1c_poly,tj_cat[j])
    e2a_fit[j] = np.polyval(e2a_poly,tj_cat[j])
    e2b_fit[j] = np.polyval(e2b_poly,tj_cat[j])
    e2c_fit[j] = np.polyval(e2c_poly,tj_cat[j])
    NCYC.e1a_res_cat = NCYC.e1a_cat - e1a_fit
    NCYC.e1b_res_cat = NCYC.e1b_cat - e1b_fit
    NCYC.e1c_res_cat = NCYC.e1c_cat - e1c_fit
    NCYC.e2a_res_cat = NCYC.e2a_cat - e2a_fit
    NCYC.e2b_res_cat = NCYC.e2b_cat - e2b_fit
    NCYC.e2c_res_cat = NCYC.e2c_cat - e2c_fit

    # include the nan values to lift pen between motor runs
    ja = np.nonzero((NCYC.abi_cat == ord('a')) | np.isnan(NCYC.abi_cat))[0]
    jb = np.nonzero((NCYC.abi_cat == ord('b')) | np.isnan(NCYC.abi_cat))[0]


    mrkp = 'r'  # marker for pinched EF
    mrku = 'c'  # marker for unpinched EF
    mrki = 'k.' # marker for IES comms

    fig = plt.figure(num=1,figsize=(10, 7))
    fig.clf()
    if options.interactive:
      fig.canvas.set_window_title('CAL with IES pings')
    fig.suptitle(pltnam + '\nred:pinched, cyan:unpinched, blk:IES comms\n' + pytstr + '\n' + wsnams)

    axlist = []

    ax = fig.add_subplot(7,2,1)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1a_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e1a_cat[ja],mrkp)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e1a_cat[jb],mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a, uV')
    plt.title('raw CAL cat, tskip={0:.1f} s'.format(options.tskip))
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    xlim1 = ax.get_xlim()
    axlist.append(ax)

    ax = fig.add_subplot(7,2,3)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1b_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e1b_cat[ja],mrku)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e1b_cat[jb],mrkp)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,5)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1c_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat,      NCYC.e1c_cat,mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1c, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,7)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2a_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e2a_cat[ja],mrkp)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e2a_cat[jb],mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,9)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2b_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e2b_cat[ja],mrku)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e2b_cat[jb],mrkp)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2b, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,11)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2c_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat,      NCYC.e2c_cat,mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2c, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    xlim11 = ax.get_xlim()
    axlist.append(ax)

    ax = fig.add_subplot(7,2,2)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1a_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e1a_res_cat[ja],mrkp)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e1a_res_cat[jb],mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1a-fit, uV')
    plt.title('raw CAL cat minus fit(order={0})'.format(META.fit_order))
    fixlims(ax)
    ax.yaxis.set_major_locator(MaxNLocator(5))
    axlist.append(ax)

    ax = fig.add_subplot(7,2,4)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1b_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e1b_res_cat[ja],mrku)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e1b_res_cat[jb],mrkp)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1b-fit, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,6)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e1c_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat,      NCYC.e1c_res_cat,mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e1c-fit, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,8)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2a_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e2a_res_cat[ja],mrkp)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e2a_res_cat[jb],mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2a-fit, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,10)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2b_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat[ja],  NCYC.e2b_res_cat[ja],mrku)
    ax.plot(NCYC.tpl_cat[jb],  NCYC.e2b_res_cat[jb],mrkp)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2b-fit, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,12)
    ax.hold(True)
    ax.plot(NCYC.tpl_cat[jies],NCYC.e2c_res_cat[jies],mrki)
    ax.plot(NCYC.tpl_cat,      NCYC.e2c_res_cat,mrku)
    ax.hold(False)
    ax.grid('on')
    ax.xaxis.set_ticklabels([])
    plt.ylabel('e2c-fit, uV')
    ax.yaxis.set_major_locator(MaxNLocator(5))
    fixlims(ax)
    axlist.append(ax)

    ax = fig.add_subplot(7,2,13)
#   print('len(NCYC.tpf_cat)=',len(NCYC.tpf_cat),'len(NCYC.stf_cat)=',len(NCYC.stf_cat))
    ax.plot(NCYC.tpf_cat, NCYC.stf_cat,'b')
    ax.grid('on')
    plt.ylabel('stat_fast')
    plt.xlabel('time, s, ' + pytstr)
    ax.set_xlim(xlim11)
    fixlims(ax)
    axlist.append(ax)
    

    ax = fig.add_subplot(7,2,14)
    ax.plot(NCYC.tpl_cat, NCYC.sta_cat,'b')
    ax.grid('on')
    plt.ylabel('stat_pick')
    plt.xlabel('time, s, ' + pytstr)
    fixlims(ax)
    axlist.append(ax)

    for ax in axlist:
      ax.set_xlim(xlim1)

    # put bigger gaps between plots
    fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.3, hspace=None)

    if options.interactive:
      fig.show()
      print('Click figure window to continue')
      plt.waitforbuttonpress()
    else:
      pdffile = options.pdfroot + '/raw-cal-cat/'+pltnam+'.pdf'
      print(pdffile)
      mymkdir(options.pdfroot + '/raw-cal-cat')
      plt.savefig(pdffile)
      os.system('updateframe.run pdf/raw-cal-cat')

def plot_hef_overlay():

  if len(hef_ind_a) == 0 and len(hef_ind_b) == 0:
    return

  # convert ADC counts to microvolts at preamp input
  # flip sign of "b" because "b" preamp input opposite of "a"

  secs_a = np.array(hef_ind_a,dtype='double') * META.tsamp_ef
  e1a_a = np.array(hef_e1a_a,dtype='double')
  e1b_a = np.array(hef_e1b_a,dtype='double')
  e1c_a = np.array(hef_e1c_a,dtype='double')
  e2a_a = np.array(hef_e2a_a,dtype='double')
  e2b_a = np.array(hef_e2b_a,dtype='double')
  e2c_a = np.array(hef_e2c_a,dtype='double')

  secs_b = np.array(hef_ind_b,dtype='double') * META.tsamp_ef
  e1a_b = np.array(hef_e1a_b,dtype='double')
  e1b_b = np.array(hef_e1b_b,dtype='double')
  e1c_b = np.array(hef_e1c_b,dtype='double')
  e2a_b = np.array(hef_e2a_b,dtype='double')
  e2b_b = np.array(hef_e2b_b,dtype='double')
  e2c_b = np.array(hef_e2c_b,dtype='double')

  secs_u = np.array(hef_ind_u,dtype='double') * META.tsamp_ef
  e1a_u = np.array(hef_e1a_u,dtype='double')
  e1b_u = np.array(hef_e1b_u,dtype='double')
  e1c_u = np.array(hef_e1c_u,dtype='double')
  e2a_u = np.array(hef_e2a_u,dtype='double')
  e2b_u = np.array(hef_e2b_u,dtype='double')
  e2c_u = np.array(hef_e2c_u,dtype='double')

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('HEF Ocean')
  fig.suptitle(pltnambase + '\nred:A, blu:B' + '\n' + wsnams)

  axlist = []

  ax = fig.add_subplot(6,2,1)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e1a_a,'r.-')
  ax.grid('on')
  plt.ylabel('e1a, uV')
  plt.title('A pinched')

  ax = fig.add_subplot(6,2,3)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e1b_a,'b.-')
  ax.grid('on')
  plt.ylabel('e1b, uV')

  ax = fig.add_subplot(6,2,5)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e1c_a,'b.-')
  ax.grid('on')
  plt.ylabel('e1c, uV')

  ax = fig.add_subplot(6,2,7)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e2a_a,'r.-')
  ax.grid('on')
  plt.ylabel('e2a, uV')

  ax = fig.add_subplot(6,2,9)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e2b_a,'b.-')
  ax.grid('on')
  plt.ylabel('e2b, uV')
  plt.xlabel('time, s, ' + pytstr)

  ax = fig.add_subplot(6,2,11)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e2c_a,'b.-')
  ax.grid('on')
  plt.ylabel('e2c, uV')
  plt.xlabel('time, s, ' + pytstr)

  ax = fig.add_subplot(6,2,2)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e1a_b,'r.-')
  ax.grid('on')
  plt.title('B pinched')

  ax = fig.add_subplot(6,2,4)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e1b_b,'b.-')
  ax.grid('on')

  ax = fig.add_subplot(6,2,6)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e1c_b,'b.-')
  ax.grid('on')

  ax = fig.add_subplot(6,2,8)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e2a_b,'r.-')
  ax.grid('on')

  ax = fig.add_subplot(6,2,10)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e2b_b,'b.-')
  ax.grid('on')

  ax = fig.add_subplot(6,2,12)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e2c_b,'b.-')
  ax.grid('on')
  plt.xlabel('time, s, ' + pytstr)

  plt.draw()


def plot_cal():

  if len(cal_ind_a) == 0 and len(cal_ind_b) == 0:
    return

  # convert ADC counts to microvolts at preamp input
  # flip sign of "b" because "b" preamp input opposite of "a"

  secs_a = np.array(cal_ind_a,dtype='double') * META.tsamp_ef
  e1c_a = np.array(cal_e1c_a,dtype='double')
  e1a_a = np.array(cal_e1a_a,dtype='double')
  e1b_a = np.array(cal_e1b_a,dtype='double')
  e2c_a = np.array(cal_e2c_a,dtype='double')
  e2a_a = np.array(cal_e2a_a,dtype='double')
  e2b_a = np.array(cal_e2b_a,dtype='double')

  secs_b = np.array(cal_ind_b,dtype='double') * META.tsamp_ef
  e1c_b = np.array(cal_e1c_b,dtype='double')
  e1a_b = np.array(cal_e1a_b,dtype='double')
  e1b_b = np.array(cal_e1b_b,dtype='double')
  e2c_b = np.array(cal_e2c_b,dtype='double')
  e2a_b = np.array(cal_e2a_b,dtype='double')
  e2b_b = np.array(cal_e2b_b,dtype='double')

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('HEF Cal Data')
  fig.suptitle(pltnambase + '\nred:A, blu:B, grn:C' + '\n' + wsnams)

  axlist = []

  ax = fig.add_subplot(2,2,1)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e1a_a,'r.-')
  ax.hold('on')
  ax.plot(secs_a,HEF_AVG.e1b_a,'b.-')
  ax.plot(secs_a,HEF_AVG.e1c_a,'g.-')
  ax.hold('off')
  ax.grid('on')
  plt.ylabel('e1, uV')
  plt.title('A pinched')

  ax = fig.add_subplot(2,2,2)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e1a_b,'r.-')
  ax.hold('on')
  ax.plot(secs_b,HEF_AVG.e1b_b,'b.-')
  ax.plot(secs_b,HEF_AVG.e1c_b,'g.-')
  ax.hold('off')
  ax.grid('on')
  plt.title('B pinched')

  ax = fig.add_subplot(2,2,3)
  axlist.append(ax)
  ax.plot(secs_a,HEF_AVG.e2a_a,'r.-')
  ax.hold('on')
  ax.plot(secs_a,HEF_AVG.e2b_a,'b.-')
  ax.plot(secs_a,HEF_AVG.e2c_a,'g.-')
  ax.hold('off')
  ax.grid('on')
  plt.ylabel('e2, uV')

  ax = fig.add_subplot(2,2,4)
  axlist.append(ax)
  ax.plot(secs_b,HEF_AVG.e2a_b,'r.-')
  ax.hold('on')
  ax.plot(secs_b,HEF_AVG.e2b_b,'b.-')
  ax.plot(secs_b,HEF_AVG.e2c_b,'g.-')
  ax.hold('off')
  ax.grid('on')

  plt.draw()

def plot_compass():

  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  if options.interactive:
    fig.canvas.set_window_title('Compass')
  fig.suptitle(pltnambase + ' plot_compass')

  secs = COMP.uxt
  secs = secs - secs[-1]

  axlist = []

  ax = fig.add_subplot(4,1,1)
  axlist.append(ax)
  ax.plot(secs,COMP.hdg,'ro-')
  ax.grid('on')
  plt.ylabel('Hdg')

  ax = fig.add_subplot(4,1,2)
  axlist.append(ax)
  ax.plot(secs,COMP.pitch,'bo-')
  ax.grid('on')
  plt.ylabel('Pitch')

  ax = fig.add_subplot(4,1,3)
  axlist.append(ax)
  ax.plot(secs,COMP.roll,'bo-')
  ax.grid('on')
  plt.ylabel('Roll')

  ax = fig.add_subplot(4,1,4)
  axlist.append(ax)
  ax.plot(secs,COMP.temp,'bo-')
  ax.grid('on')
  plt.ylabel('Temp')

  plt.xlabel('Time, s, ' + pytstr)
  plt.draw()

def initbufs():
  global aux_uxt, aux_tt1, aux_tt2, aux_tt3, aux_tt4
  global aux_pres, aux_temp, aux_btemp, aux_bfreq
  global aux_uxt_xfr
  aux_uxt   = []
  aux_tt1    = []
  aux_tt2    = []
  aux_tt3    = []
  aux_tt4    = []
  aux_pres   = []
  aux_temp   = []
  aux_btemp  = []
  aux_bfreq  = []
  aux_uxt_xfr = []

  global comp_uxt, comp_hdg, comp_pitch, comp_roll, comp_temp
  comp_uxt  = []
  comp_hdg   = []
  comp_pitch = []
  comp_roll  = []
  comp_temp  = []

  global mot_ind_a, mot_cur_a, mot_ind_b, mot_cur_b, mot_ind_u, mot_cur_u
  mot_ind_a = []
  mot_ind_b = []
  mot_ind_u = []
  mot_cur_a = []
  mot_cur_b = []
  mot_cur_u = []

  global hef_ind_a, hef_ind_b
  hef_ind_a = []
  hef_ind_b = []
  hef_ind_u = []

  global hef_e1a_a, hef_e1b_a, hef_e1c_a, hef_e2a_a, hef_e2b_a, hef_e2c_a
  hef_e1a_a = []
  hef_e1b_a = []
  hef_e1c_a = []
  hef_e2a_a = []
  hef_e2b_a = []
  hef_e2c_a = []

  global hef_e1a_b, hef_e1b_b, hef_e1c_b, hef_e2a_b, hef_e2b_b, hef_e2c_b
  hef_e1a_b = []
  hef_e1b_b = []
  hef_e1c_b = []
  hef_e2a_b = []
  hef_e2b_b = []
  hef_e2c_b = []

  global cal_ind_a, cal_ind_b, cal_ind_u
  cal_ind_a = []
  cal_ind_b = []
  cal_ind_u = []

  global cal_e1a_a, cal_e1b_a, cal_e1c_a
  cal_e1a_a = []
  cal_e1b_a = []
  cal_e1c_a = []

  global cal_e2a_a, cal_e2b_a, cal_e2c_a
  cal_e2a_a = []
  cal_e2b_a = []
  cal_e2c_a = []

  global cal_e1a_b, cal_e1b_b, cal_e1c_b
  cal_e1a_b = []
  cal_e1b_b = []
  cal_e1c_b = []

  global cal_e2a_b, cal_e2b_b, cal_e2c_b
  cal_e2a_b = []
  cal_e2b_b = []
  cal_e2c_b = []

  global cal_e1a_u, cal_e1b_u, cal_e1c_u
  cal_e1a_u = []
  cal_e1b_u = []
  cal_e1c_u = []

  global cal_e2a_u, cal_e2b_u, cal_e2c_u
  cal_e2a_u = []
  cal_e2b_u = []
  cal_e2c_u = []

def initcounts():
  global count_tod, count_hef, count_aux
  count_tod = 0
  count_hef = 0
  count_aux = 0

  global count_crc_ok, count_bad_crc_check, count_bad_crc_decode
  count_crc_ok = 0
  count_bad_crc_check = 0
  count_bad_crc_decode = 0

  global count_crc_missing; count_crc_missing  = 0

  global count_M
  count_M = 0

  global count_EC
  count_EC = 0

  global count_E
  count_E = 0

  global count_compass
  count_compass = 0

  global count_C
  count_C = 0

  global count_hef_unknown
  count_hef_unknown = 0

  global lineno, count_lineok, count_skipped
  lineno = 0
  count_lineok = 0
  count_skipped = 0

  global count_removed
  count_removed = 0

def printcounts():
  print('hpproc.py: lineno              =',lineno)
  print('  count_skipped      =',count_skipped)
  print('  count_lineok        =',count_lineok)
  print('  count_crc_ok        =',count_crc_ok)
  print('  count_crc_missing   =',count_crc_missing)
  print('  count_bad_crc_check =',count_bad_crc_check)
  print('  count_bad_crc_decode=',count_bad_crc_decode)
  print('  count_tod=',count_tod)
  print('  count_hef=',count_hef)
  print('  count_aux=',count_aux)
  print('  count_EC =',count_EC)
  print('  count_E  =',count_E)
  print('  count_C  =',count_C)
  print('  count_M  =',count_M)
  print('  count_compass=',count_compass)
  print('  count_hef_unknown=',count_hef_unknown)
  print('  count_removed=',count_removed)

def printsizes():
  print('hpproc.py: len(aux_uxt  )=',len( aux_uxt   ))
  if options.verbose:
    print('len(aux_tt1   )=',len( aux_tt1    ))
    print('len(aux_tt2   )=',len( aux_tt2    ))
    print('len(aux_tt3   )=',len( aux_tt3    ))
    print('len(aux_tt4   )=',len( aux_tt4    ))
    print('len(aux_pres  )=',len( aux_pres   ))
    print('len(aux_temp  )=',len( aux_temp   ))
    print('len(aux_btemp )=',len( aux_btemp  ))
    print('len(aux_bfreq )=',len( aux_bfreq  ))

  print('len(comp_uxt )=',len( comp_uxt  ))
  if options.verbose:
    print('len(comp_hdg  )=',len( comp_hdg   ))
    print('len(comp_pitch)=',len( comp_pitch ))
    print('len(comp_roll )=',len( comp_roll  ))
    print('len(comp_temp )=',len( comp_temp  ))

  print('len(mot_ind_a )=',len(mot_ind_a ))
  if options.verbose:
    print('len(mot_cur_a )=',len(mot_cur_a ))

  print('len(mot_ind_b )=',len(mot_ind_b ))
  if options.verbose:
    print('len(mot_cur_b )=',len(mot_cur_b ))

  print('len(hef_ind_a )=',len(hef_ind_a ))
  if options.verbose:
    print('len(hef_e1a_a )=',len(hef_e1a_a ))
    print('len(hef_e1b_a )=',len(hef_e1b_a ))
    print('len(hef_e1c_a )=',len(hef_e1c_a ))
    print('len(hef_e2a_a )=',len(hef_e2a_a ))
    print('len(hef_e2b_a )=',len(hef_e2b_a ))
    print('len(hef_e2c_a )=',len(hef_e2c_a ))

  print('len(hef_ind_b )=',len(hef_ind_b ))
  if options.verbose:
    print('len(hef_e1a_b )=',len(hef_e1a_b ))
    print('len(hef_e1b_b )=',len(hef_e1b_b ))
    print('len(hef_e1c_b )=',len(hef_e1c_b ))
    print('len(hef_e2a_b )=',len(hef_e2a_b ))
    print('len(hef_e2b_b )=',len(hef_e2b_b ))
    print('len(hef_e2c_b )=',len(hef_e2c_b ))

  print('len(cal_ind_a )=',len(cal_ind_a ))
  if options.verbose:
    print('len(cal_e1a_a )=',len(cal_e1a_a ))
    print('len(cal_e1b_a )=',len(cal_e1b_a ))
    print('len(cal_e1c_a )=',len(cal_e1c_a ))
    print('len(cal_e2a_a )=',len(cal_e2a_a ))
    print('len(cal_e2b_a )=',len(cal_e2b_a ))
    print('len(cal_e2c_a )=',len(cal_e2c_a ))

  print('len(cal_ind_b )=',len(cal_ind_b ))
  if options.verbose:
    print('len(cal_e1a_b )=',len(cal_e1a_b ))
    print('len(cal_e1b_b )=',len(cal_e1b_b ))
    print('len(cal_e1c_b )=',len(cal_e1c_b ))
    print('len(cal_e2a_b )=',len(cal_e2a_b ))
    print('len(cal_e2b_b )=',len(cal_e2b_b ))
    print('len(cal_e2c_b )=',len(cal_e2c_b ))

def print_info_and_exit():
  print('ERROR:')
  print('  ifile=',ifile)
  print('  lineno=',lineno)
  print('  linein=',linein)
  sys.exit(1)

def gz_open(ifile):
  if ifile[-3:] == '.gz':
    try:
      ifd = gzip.open(ifile,'r')
    except:
      print('cannot open ifile=',ifile)
      sys.exit(1)
  else:
    try:
      ifd = gzip.open(ifile+'.gz','r')
    except:
      try:
        ifd = open(ifile,'rt')
      except:
        print('cannot open ifile=',ifile)
        sys.exit(1)
  return ifd

def larger_axlim( axlim ):
    """ argument axlim expects 2-tuple 
        returns slightly larger 2-tuple """
    axmin,axmax = axlim
    axrng = axmax - axmin
    new_min = axmin - 0.03 * axrng
    new_max = axmax + 0.03 * axrng
    return new_min,new_max

def fixlims(ax):
  ax.set_xlim( larger_axlim( ax.get_xlim() ) )
  ax.set_ylim( larger_axlim( ax.get_ylim() ) )

def init_HEF_FIFO():
  global HEF_FIFO
  HEF_FIFO = collections.deque(maxlen=options.fifolen)

def init_HEF_DEMOD():
  global HEF_DEMOD
  HEF_DEMOD = collections.namedtuple('HEF_DEMOD',[])
  HEF_DEMOD.nout = 0
  HEF_DEMOD.hcno = []
  HEF_DEMOD.uxt = []
  HEF_DEMOD.abm = []
  HEF_DEMOD.e1o = []
  HEF_DEMOD.e2o = []
  HEF_DEMOD.e1std = []
  HEF_DEMOD.e2std = []

  HEF_DEMOD.e1co = []
  HEF_DEMOD.e2co = []
  HEF_DEMOD.e1c_std = []
  HEF_DEMOD.e2c_std = []

def compute_hef_demod(x):
  global HEF_FIFO, HEF_DEMOD
  ok = True;

  if x.abu == 'a' or x.abu == 'b':
    HEF_FIFO.append(x)
  else:
    return

  if len(HEF_FIFO) == options.fifolen:
    npp = []
    hcno = []
    uxt = []
    ind = []
    abu = []
    e1a = []
    e1b = []
    e1c = []
    e2a = []
    e2b = []
    e2c = []
    sqw = []
    wts = []
    cou = []
    for i in range(options.fifolen):
      n = len(HEF_FIFO[i].uxt)
      npp = np.append(npp,n)
      uxt = np.append(uxt,HEF_FIFO[i].uxt)
      hcno = np.append(hcno,HEF_FIFO[i].hcno)
      ind = np.append(ind,HEF_FIFO[i].ind)
      abu = np.append(abu,HEF_FIFO[i].abu)
      e1a = np.append(e1a,HEF_FIFO[i].e1a)
      e1b = np.append(e1b,HEF_FIFO[i].e1b)
      e1c = np.append(e1c,HEF_FIFO[i].e1c)
      e2a = np.append(e2a,HEF_FIFO[i].e2a)
      e2b = np.append(e2b,HEF_FIFO[i].e2b)
      e2c = np.append(e2c,HEF_FIFO[i].e2c)
      cou = np.append(cou,HEF_FIFO[i].cou)

      if HEF_FIFO[i].abu == 'a':
        sqw = np.append(sqw,+0.5*np.ones(n))
      if HEF_FIFO[i].abu == 'b':
        sqw = np.append(sqw,-0.5*np.ones(n))

      t = HEF_FIFO[i].uxt - HEF_FIFO[i].uxt[0]
      k = np.nonzero(t < options.tskip)[0]
      w = np.ones(n, dtype=np.int)
      w[k] = 0
      wts = np.append(wts,w)

    if not all(np.diff(hcno) == 1):
#     print('hpproc: skipped: hcno must increment by 1, diff(hcno)=',np.diff(hcno))
      ok = False
      # return

    if ok:
      tnd = uxt - np.mean(uxt)
      tnd = tnd / (tnd[-1] - tnd[0])
      con = np.ones(len(tnd))

  #   if options.fifolen == 3:
  #     ma = ['b','a','b']
  #     mb = ['a','b','a']
  #   elif options.fifolen == 3:
  #     ma = ['a','b','a','b','a']
  #     mb = ['b','a','b','a','b']
  #   else:
  #     print('error: hpproc: not coded yet for fifolen=',options.fifolen)
  #     sys.exit(1)

      ma = []
      mb = []
      for i in range(options.fifolen):
        if i % 2 == 0:
          ma.append('b')
          mb.append('a')
        else:
          ma.append('a')
          mb.append('b')

      ju = np.nonzero(wts > 0)[0]

      if len(ju) > 3:
        uxtmean = np.mean(uxt[ju])

        # basis matrix for least squares fitting using np.linalg.lstsq()
        Bju = np.vstack([sqw[ju], tnd[ju], con[ju]]).T

        if all(abu == ma):
          abm = 'a'
          e1sqw, e1tnd, e1con = np.linalg.lstsq(Bju, e1a[ju])[0]
          e2sqw, e2tnd, e2con = np.linalg.lstsq(Bju, e2a[ju])[0]
          e1res = e1a[ju] - (e1sqw * sqw[ju] + e1tnd * tnd[ju] + e1con)
          e2res = e2a[ju] - (e2sqw * sqw[ju] + e2tnd * tnd[ju] + e2con)
        elif all(abu == mb):
          abm = 'b'
          # flipped signs here because e1b and e2b signs flipped much earlier
          e1sqw, e1tnd, e1con = np.linalg.lstsq(Bju, -e1b[ju])[0]
          e2sqw, e2tnd, e2con = np.linalg.lstsq(Bju, -e2b[ju])[0]
          e1res = -e1b[ju] - (e1sqw * sqw[ju] + e1tnd * tnd[ju] + e1con)
          e2res = -e2b[ju] - (e2sqw * sqw[ju] + e2tnd * tnd[ju] + e2con)
        else:
          print('hpproc: skipped abu=',abu,'must be ma=',ma,'or mb=',mb)
          ok = False
          #return

        mf = int(options.fifolen / 2)
        if abm != abu[mf]:
          print('error: hpproc: abm != abu[mf]')
          sys.exit(1)

    if ok:

      e1std = np.std(e1res)
      e2std = np.std(e2res)

      e1c_sqw, e1c_tnd, e1c_con = np.linalg.lstsq(Bju, e1c[ju])[0]
      e2c_sqw, e2c_tnd, e2c_con = np.linalg.lstsq(Bju, e2c[ju])[0]
      e1c_res = e1c[ju] - (e1c_sqw * sqw[ju] + e1c_tnd * tnd[ju] + e1c_con)
      e2c_res = e2c[ju] - (e2c_sqw * sqw[ju] + e2c_tnd * tnd[ju] + e2c_con)
      e1c_std = np.std(e1c_res)
      e2c_std = np.std(e2c_res)

    else:
      uxtmean = np.nan
      abm = 'u'
      e1std = e1sqw = e1tnd = e1con = np.nan
      e2std = e2sqw = e2tnd = e2con = np.nan
      e1c_std = e1c_sqw = e1c_tnd = e1c_con = np.nan
      e2c_std = e2c_sqw = e2c_tnd = e2c_con = np.nan

    HEF_DEMOD.nout += 1
    HEF_DEMOD.uxt.append(uxtmean)
    HEF_DEMOD.abm.append(abm)
    HEF_DEMOD.e1o.append(e1sqw)
    HEF_DEMOD.e2o.append(e2sqw)
    HEF_DEMOD.e1std.append(e1std)
    HEF_DEMOD.e2std.append(e2std)

    HEF_DEMOD.e1co.append(e1c_sqw)
    HEF_DEMOD.e2co.append(e2c_sqw)
    HEF_DEMOD.e1c_std.append(e1c_std)
    HEF_DEMOD.e2c_std.append(e2c_std)

def plot_hef_demod():
  pdfdir = options.pdfroot + '/hef-avg/'
  mymkdir(pdfdir)

# pyt = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_DEMOD.uxt[0])
# pytstr = pyt.strftime('%Y-%m-%d %H:%M:%S')
  pytbeg = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[0])
  pytend = datetime(1970,1,1,0,0,0) + timedelta(0,HEF_AVG.uxt[-1])
  pytstr = pytbeg.strftime('%Y-%m-%d %H:%M:%S') + ' to ' + \
           pytend.strftime('%Y-%m-%d %H:%M:%S') + ' UTC'

  ylim_o = [-options.ylim_ocean,options.ylim_ocean]
  ylim_s = [0.0,2 * options.ylim_ocean]
  mrkrsz1 = 10
  mrkrsz1 = 3
  mrkrsz2 = 5
  lw2 = 1

  abm = np.array(HEF_DEMOD.abm)
  e1o = np.array(HEF_DEMOD.e1o)
  e2o = np.array(HEF_DEMOD.e2o)

  e1std = np.array(HEF_DEMOD.e1std)
  e2std = np.array(HEF_DEMOD.e2std)

  e1co = np.array(HEF_DEMOD.e1co)
  e2co = np.array(HEF_DEMOD.e2co)

  e1c_std = np.array(HEF_DEMOD.e1c_std)
  e2c_std = np.array(HEF_DEMOD.e2c_std)

  ja = np.nonzero(abm == 'a')[0]
  jb = np.nonzero(abm == 'b')[0]
  tc  = HEF_DEMOD.uxt - HEF_DEMOD.uxt[0]

  hrsja = tc[ja] / 3600
  hrsjb = tc[jb] / 3600

  e1ao = e1o[ja]
  e1bo = e1o[jb]
  e2ao = e2o[ja]
  e2bo = e2o[jb]

  e1coa = e1co[ja]
  e1cob = e1co[jb]
  e2coa = e2co[ja]
  e2cob = e2co[jb]

  m = 6
  na = int(len(ja)/m)
  nb = int(len(jb)/m)

  hrsjaf = np.mean(np.reshape(hrsja[0:m*na],(m,na),order='F'),0)
  hrsjbf = np.mean(np.reshape(hrsjb[0:m*nb],(m,nb),order='F'),0)

  e1aof = np.mean(np.reshape(e1ao[0:m*na],(m,na),order='F'),0)
  e1bof = np.mean(np.reshape(e1bo[0:m*nb],(m,nb),order='F'),0)
  e2aof = np.mean(np.reshape(e2ao[0:m*na],(m,na),order='F'),0)
  e2bof = np.mean(np.reshape(e2bo[0:m*nb],(m,nb),order='F'),0)

  e1coaf = np.mean(np.reshape(e1coa[0:m*na],(m,na),order='F'),0)
  e1cobf = np.mean(np.reshape(e1cob[0:m*nb],(m,nb),order='F'),0)
  e2coaf = np.mean(np.reshape(e2coa[0:m*na],(m,na),order='F'),0)
  e2cobf = np.mean(np.reshape(e2cob[0:m*nb],(m,nb),order='F'),0)


  # plot ocean data
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-demod-ocean'.format(pltnambase)

  if options.interactive:
    fig.canvas.set_window_title('HEF Ocean Averages')

  fig.suptitle('{0}\ntskip={1:.1f} s, fifolen={2}, red:A, blu:B'.
    format(pltnam,HEF_HDR.tskip_ef,options.fifolen) + '\n' + wsnams)

  mltoff = 719529 - 366 + HEF_DEMOD.uxt[0] / 86400
  mltja = hrsja / 24.0 + mltoff
  mltjb = hrsjb / 24.0 + mltoff
  mltjaf = hrsjaf / 24.0 + mltoff
  mltjbf = hrsjbf / 24.0 + mltoff

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.hold(True)

  ax.plot_date(mltja,e1ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,e1bo,'b.',markersize=mrkrsz1)
  if options.do_plot_ec_demod:
    ax.plot_date(mltja,e1coa,'m.',markersize=mrkrsz1)
    ax.plot_date(mltjb,e1cob,'c.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,e1aof,'r.-',markersize=mrkrsz2,linewidth=lw2)
  ax.plot_date(mltjbf,e1bof,'b.-',markersize=mrkrsz2,linewidth=lw2)
  if options.do_plot_ec_demod:
    ax.plot_date(mltjaf,e1coaf,'m.-',markersize=mrkrsz2,linewidth=lw2)
    ax.plot_date(mltjbf,e1cobf,'c.-',markersize=mrkrsz2,linewidth=lw2)
  ax.hold(False)
  ax.set_ylim(ylim_o)
  ax.xaxis.set_ticklabels([])
  ax.grid(True)
  plt.ylabel(options.e1ws + ' e1o, uV')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,e2ao,'r.',markersize=mrkrsz1)
  ax.plot_date(mltjb,e2bo,'b.',markersize=mrkrsz1)
  if options.do_plot_ec_demod:
    ax.plot_date(mltja,e2coa,'m.',markersize=mrkrsz1)
    ax.plot_date(mltjb,e2cob,'c.',markersize=mrkrsz1)
  ax.plot_date(mltjaf,e2aof,'r.-',markersize=mrkrsz2,linewidth=lw2)
  ax.plot_date(mltjbf,e2bof,'b.-',markersize=mrkrsz2,linewidth=lw2)
  if options.do_plot_ec_demod:
    ax.plot_date(mltjaf,e2coaf,'m.-',markersize=mrkrsz2,linewidth=lw2)
    ax.plot_date(mltjbf,e2cobf,'c.-',markersize=mrkrsz2,linewidth=lw2)
  ax.hold(False)
  ax.set_ylim(ylim_o)
  ax.grid(True)
  plt.ylabel(options.e2ws + ' e2o, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.2)

  fix_xdates(axlist,1)

  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

  # plot std dev data
  fig = plt.figure(num=1,figsize=(10, 7))
  fig.clf()
  pltnam = '{0}-hef-demod-std-dev'.format(pltnambase)

  if options.interactive:
    fig.canvas.set_window_title('HEF Standard Deviations')

  fig.suptitle('{0}\ntskip={1:.1f} s, fifolen={2}, red:A, blu:B\n{3}'.
    format(pltnam,HEF_HDR.tskip_ef,options.fifolen,wsnams))

  axlist = []

  ax = fig.add_subplot(2,1,1)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,e1std[ja],'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjb,e1std[jb],'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim_s)
  ax.grid(True)
  plt.ylabel('e1std, uV')

  ax = fig.add_subplot(2,1,2)
  axlist.append(ax)
  ax.hold(True)
  ax.plot_date(mltja,e2std[ja],'r.-',markersize=mrkrsz2)
  ax.plot_date(mltjb,e2std[jb],'b.-',markersize=mrkrsz2)
  ax.hold(False)
  ax.set_ylim(ylim_s)
  ax.grid(True)
  plt.ylabel('e2std, uV')
  plt.xlabel(pytstr)

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=0.83, wspace=0.4, hspace=0.2)

  fix_xdates(axlist,1)
    
  if options.interactive:
    fig.show()
    print('Click figure window to continue')
    plt.waitforbuttonpress()
  else:
    mymkdir(pdfdir)
    pdffile = pdfdir + pltnam + '.pdf'
    print(pdffile)
    plt.savefig(pdffile)
    os.system('updateframe.run ' + pdfdir)

def init_HEF_AVG():
  global HEF_AVG
  HEF_AVG = collections.namedtuple('HEF_AVG',[])
  HEF_AVG.hcno = []
  HEF_AVG.uxt = []
  HEF_AVG.abu = []
  HEF_AVG.navg = []
  HEF_AVG.nuse = []
  HEF_AVG.e1am = []
  HEF_AVG.e1bm = []
  HEF_AVG.e1cm = []
  HEF_AVG.e2am = []
  HEF_AVG.e2bm = []
  HEF_AVG.e2cm = []
  HEF_AVG.e1as = []
  HEF_AVG.e1bs = []
  HEF_AVG.e1cs = []
  HEF_AVG.e2as = []
  HEF_AVG.e2bs = []
  HEF_AVG.e2cs = []
  HEF_AVG.e1afs = []
  HEF_AVG.e1bfs = []
  HEF_AVG.e1cfs = []
  HEF_AVG.e2afs = []
  HEF_AVG.e2bfs = []
  HEF_AVG.e2cfs = []

def init_CAL_DEMOD():
  global CAL_DEMOD
  CAL_DEMOD = collections.namedtuple('CAL_DEMOD',[])
  CAL_DEMOD.hcno = []
  CAL_DEMOD.uxt = []
  CAL_DEMOD.abu = []
  CAL_DEMOD.navg = []
  CAL_DEMOD.nuse = []
  CAL_DEMOD.e1a_amp = []
  CAL_DEMOD.e1b_amp = []
  CAL_DEMOD.e1c_amp = []
  CAL_DEMOD.e2a_amp = []
  CAL_DEMOD.e2b_amp = []
  CAL_DEMOD.e2c_amp = []
  CAL_DEMOD.e1a_std = []
  CAL_DEMOD.e1b_std = []
  CAL_DEMOD.e1c_std = []
  CAL_DEMOD.e2a_std = []
  CAL_DEMOD.e2b_std = []
  CAL_DEMOD.e2c_std = []

def init_MOT_CUR():
  global MOT_CUR
  MOT_CUR = collections.namedtuple('MOT_CUR',[])
  MOT_CUR.hcno = []
  MOT_CUR.uxt1 = []
  MOT_CUR.uxt2 = []
  MOT_CUR.abu = []
  MOT_CUR.cur1 = []
  MOT_CUR.cur2 = []
  MOT_CUR.dur1 = []
  MOT_CUR.dur2 = []

def writegz(pltnambase):
  ofile = options.odir + '/' + pltnambase + '-hpp.gz'
  try:
    ofp = gzip.open(ofile,'w')
  except:
    print('cannot open output file=',ofile)
    sys.exit(1)

  print('ofile=',ofile)
  if options.verbose:
    print('len(HEF_AVG.uxt)=',len(HEF_AVG.uxt))
    print('len(CAL_DEMOD.uxt)=',len(CAL_DEMOD.uxt))

  print('#META#',file=ofp)
  print('# prog          = hpproc.py',file=ofp)
  print('# runid         =',INFO.runid,file=ofp)
  print('# pltnambase      =',pltnambase,file=ofp)
  print('# ofile         =',ofile,file=ofp)
  print('# tskip         =',META.tskip,'s',file=ofp)
  print('# twait         =',META.twait,'s',file=ofp)
  print('# tsamp_ef      =',META.tsamp_ef,'s',file=ofp)
  print('# tsamp_mot     =',META.tsamp_mot,'s',file=ofp)
  print('# fit_order     =',META.fit_order,file=ofp)
  print('# esep1         =',INFO.esep1,'m',file=ofp)
  print('# esep2         =',INFO.esep2,'m',file=ofp)
  print('# magvar        =',INFO.magvar,'deg',file=ofp)
  print('# compass_angle =',INFO.compass_angle,'deg',file=ofp)
  print('# info_heading  =',INFO.heading,'deg',file=ofp)
  print('# COMP.hdg_avg  = {0:.2f} deg'.format(COMP.hdg_avg),file=ofp)
  print('# uxt_poly_ref  =',INFO.uxt_poly_ref,'s',file=ofp)
  print('# uxt_poly_coef =',INFO.uxt_poly_coef,file=ofp)
  print('# note: units "uV" are microvolts',file=ofp)
  print('# note: vars "uxt" are Unix/Linux time -- seconds since Jan 1, 1970 00:00:00 UTC',file=ofp)

  print('#HDR#',file=ofp)
  print('# name  = HEF',file=ofp)
  print('# vars  = index,hcno,uxt,abu,navg,nuse,Ex,Ey,E1,E2,e1am,e1bm,e2am,e2bm,e1as,e1bs,e2as,e2bs,e1afs,e1bfs,e2afs,e2bfs',file=ofp)
  print('# units = none,none,s,none,none,none,uV/m,uV/m,uV/m,uV/m,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV',file=ofp)
  print('# ndat  =',len(HEF_AVG.uxt),file=ofp)
  print('#DATA#',file=ofp)
  for i in range(len(HEF_AVG.uxt)):
    print('HEF,{0}'.format(i),                 file=ofp,end='')
    print(',{0}'.format(HEF_AVG.hcno[i]),      file=ofp,end='')
    print(',{0:.0f}'.format(HEF_AVG.uxt[i]),   file=ofp,end='')
    print(',{0}'.format(HEF_AVG.abu[i]),        file=ofp,end='')
    print(',{0}'.format(HEF_AVG.navg[i]),      file=ofp,end='')
    print(',{0}'.format(HEF_AVG.nuse[i]),      file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.Ex[i]),    file=ofp,end='') # ocean electric field, earth coords
    print(',{0:.3f}'.format(HEF_AVG.Ey[i]),    file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.E1[i]),    file=ofp,end='') # ocean electric field, arm coords
    print(',{0:.3f}'.format(HEF_AVG.E2[i]),    file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e1am[i]),  file=ofp,end='') # means
    print(',{0:.3f}'.format(HEF_AVG.e1bm[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2am[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2bm[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e1as[i]),  file=ofp,end='') # std dev
    print(',{0:.3f}'.format(HEF_AVG.e1bs[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2as[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2bs[i]),  file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e1afs[i]), file=ofp,end='') # fancy std dev
    print(',{0:.3f}'.format(HEF_AVG.e1bfs[i]), file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2afs[i]), file=ofp,end='')
    print(',{0:.3f}'.format(HEF_AVG.e2bfs[i]), file=ofp,end='')
    print(file=ofp) # end of line

  print('#HDR#',file=ofp)
  print('# name  = CAL',file=ofp)
  print('# vars  = index,hcno,uxt,abu,navg,nuse,e1am,e1bm,e1cm,e2am,e2bm,e2cm,e1as,e1bs,e1cs,e2as,e2bs,e2cs',file=ofp)
  print('# units = none,none,s,none,none,none,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV,uV',file=ofp)
  print('# ndat  =',len(CAL_DEMOD.uxt),file=ofp)
  print('#DATA#',file=ofp)
  for i in range(len(CAL_DEMOD.uxt)):
    print('CAL,{0}'.format(i),                     file=ofp,end='')
    print(',{0}'.format(CAL_DEMOD.hcno[i]),        file=ofp,end='')
    print(',{0:.0f}'.format(CAL_DEMOD.uxt[i]),     file=ofp,end='')
    print(',{0}'.format(CAL_DEMOD.abu[i]),          file=ofp,end='')
    print(',{0}'.format(CAL_DEMOD.navg[i]),        file=ofp,end='')
    print(',{0}'.format(CAL_DEMOD.nuse[i]),        file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1am[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1bm[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1cm[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2am[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2bm[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2cm[i]),     file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1as[i]), file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1bs[i]), file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e1cs[i]), file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2as[i]), file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2bs[i]), file=ofp,end='')
    print(',{0:.3f}'.format(CAL_DEMOD.e2cs[i]), file=ofp,end='')
    print(file=ofp) # end of line

  print('#HDR#',file=ofp)
  print('# name  = MOT',file=ofp)
  print('# vars  = index,hcno,uxt1,uxt2,abu,dur1,dur2,cur1,cur2',file=ofp)
  print('# units = none,none,s,s,none,s,s,mA,mA',file=ofp)
  print('# ndat  =',len(MOT_CUR.uxt1),file=ofp)
  print('#DATA#',file=ofp)
  for i in range(len(MOT_CUR.uxt1)):
    print('MOT,{0}'.format(i),               file=ofp,end='')
    print(',{0}'.format(MOT_CUR.hcno[i]),    file=ofp,end='')
    print(',{0:.0f}'.format(MOT_CUR.uxt1[i]),file=ofp,end='')
    print(',{0:.0f}'.format(MOT_CUR.uxt2[i]),file=ofp,end='')
    print(',{0}'.format(MOT_CUR.abu[i]),      file=ofp,end='')
    print(',{0:.3f}'.format(MOT_CUR.dur1[i]),file=ofp,end='')
    print(',{0:.3f}'.format(MOT_CUR.dur2[i]),file=ofp,end='')
    print(',{0:.0f}'.format(MOT_CUR.cur1[i]),file=ofp,end='')
    print(',{0:.0f}'.format(MOT_CUR.cur2[i]),file=ofp,end='')
    print(file=ofp) # end of line

  print('#HDR#',file=ofp)
  print('# name  = IES',file=ofp)
  print('# vars  = index,uxt_ies,tt1,tt2,tt3,tt4,pres,temp,btemp,bfreq,uxt_xfr',file=ofp)
  print('# units = none,s,s,s,s,s,dbar,degC,degC,Hz,s',file=ofp)
  print('# ndat  =',len(aux_uxt),file=ofp)
  print('#DATA#',file=ofp)
  for i in range(len(aux_uxt)):
    print('IES,{0}'.format(i)                    ,file=ofp,end='')
    print(',{0:.0f}'.format(aux_uxt[i]   * 0.1)  ,file=ofp,end='')
    print(',{0:.5f}'.format(aux_tt1[i]   * 1e-5) ,file=ofp,end='')
    print(',{0:.5f}'.format(aux_tt2[i]   * 1e-5) ,file=ofp,end='')
    print(',{0:.5f}'.format(aux_tt3[i]   * 1e-5) ,file=ofp,end='')
    print(',{0:.5f}'.format(aux_tt4[i]   * 1e-5) ,file=ofp,end='')
    print(',{0:.3f}'.format(aux_pres[i]  * 1e-3) ,file=ofp,end='')
    print(',{0:.3f}'.format(aux_temp[i]  * 1e-3) ,file=ofp,end='')
    print(',{0:.3f}'.format(aux_btemp[i] * 1e-3) ,file=ofp,end='')
    print(',{0}'.format(aux_bfreq[i])            ,file=ofp,end='')
    print(',{0:.0f}'.format(aux_uxt_xfr[i])      ,file=ofp,end='')
    print(file=ofp) # end of line

  print('#HDR#',file=ofp)
  print('# name  = COMPASS',file=ofp)
  print('# vars  = index,uxt,heading,pitch,roll,temperature',file=ofp)
  print('# units = none,s,deg,deg,deg,degC',file=ofp)
  print('# ndat  =',len(comp_uxt),file=ofp)
  print('# info_heading  = ',INFO.heading,file=ofp)
  print('#DATA#',file=ofp)
  for i in range(len(comp_uxt)):
    print('COMPASS,{0}'.format(i),          file=ofp,end='')
    print(',{0:.0f}'.format(COMP.uxt[i]),   file=ofp,end='')
    print(',{0:.1f}'.format(COMP.hdg[i]),   file=ofp,end='')
    print(',{0:.1f}'.format(COMP.pitch[i]), file=ofp,end='')
    print(',{0:.1f}'.format(COMP.roll[i]),  file=ofp,end='')
    print(',{0:.1f}'.format(COMP.temp[i]),  file=ofp,end='')
    print(file=ofp) # end of line

  print('#EOD#',file=ofp)
  ofp.close()

def compute_compass():
  global COMP
  COMP = collections.namedtuple('COMP',[])

  COMP.uxt   = np.array(comp_uxt)
  COMP.hdg   = np.array(comp_hdg)   * 360.0 / 4096.0
  COMP.pitch = np.array(comp_pitch) *  90.0 / 4096.0
  COMP.roll  = np.array(comp_roll ) * 180.0 / 4096.0
  COMP.temp  = np.array(comp_temp ) * 0.1

  jz = np.nonzero((COMP.hdg   <    0.0) | (COMP.hdg   > 360.0) | \
                  (COMP.pitch <  -90.0) | (COMP.pitch >  90.0) | \
                  (COMP.roll  < -180.0) | (COMP.roll  > 180.0))[0]
  if len(jz) > 0:
    COMP.hdg  [jz] = np.nan
    COMP.pitch[jz] = np.nan
    COMP.roll [jz] = np.nan
    COMP.temp [jz] = np.nan

  j = np.nonzero(np.isfinite(COMP.hdg))[0]
  COMP.hdg_avg   = np.mean(COMP.hdg[j])
  COMP.pitch_avg = np.mean(COMP.pitch[j])
  COMP.roll_avg  = np.mean(COMP.roll[j])
  COMP.temp_avg  = np.mean(COMP.temp[j])

def decode_time(txt):
  try:
    return timegm(strptime(txt,'%Y-%m-%d_%H:%M:%S'))
  except:
    pass
  try:
    return timegm(strptime(txt,'%Y%m%d_%H%M%S'))
  except:
    pass
  try:
    return timegm(strptime(txt,'%Y%m%dT%H%M%S'))
  except:
    pass
  return None

def mymkdir(mydir):
  try: 
    os.makedirs(mydir) # like mkdir -p
  except OSError:
    if not os.path.isdir(mydir):
      print('cannot make directory:',mydir)
      sys.exit(1)

# main
if __name__ == '__main__':

  HEF_HDR = collections.namedtuple('HEF_HDR',[])


  aux_flag = False
  plot_raw_cal_flag = False
  plot_raw_mot_flag = False

  parser = OptionParser(
    usage="%prog [Options] ifile[s]",
    version="%prog 1.0")

  parser.add_option("-v", "--verbose", dest="verbose", 
    action="store_true", default=False,
    help="print debug info to stdout")

  parser.add_option("-x", "--extra", dest="do_extra",
    action="store_true", default=False,
    help="print counts and sizes")

  parser.add_option("-i", "--interactive", dest="interactive", 
    action="store_true", default=False,
    help="interactive plotting")

  parser.add_option("-o", "--plot_overlay", dest="do_plot_overlay", 
    action="store_true", default=False,
    help="plot data after each pinch on same axis")

  parser.add_option("-e", "--plot_raw_ef", dest="do_plot_raw_ef", 
    action="store_true", default=False,
    help="plot raw hef, green = IES data request times")

  parser.add_option("-c", "--plot_raw_cal", dest="do_plot_raw_cal", 
    action="store_true", default=False,
    help="plot raw cal, green = IES data request times")

  parser.add_option("-d", "--plot_cal_demod", dest="do_plot_cal_demod", 
    action="store_true", default=False,
    help="plot demodulated calibrator output")

  parser.add_option("--plot_mot_cur_dur", dest="do_plot_mot_cur_dur", 
    action="store_true", default=False,
    help="plot motor current and duration")

  parser.add_option("-m", "--plot_raw_mot", dest="do_plot_raw_mot", 
    action="store_true", default=False,
    help="plot raw motor current")

  parser.add_option("-k", "--plot_hef_demod", dest="do_plot_hef_demod", 
    action="store_true", default=False,
    help="plot HEF demod")

  parser.add_option("--plot_ec_demod", dest="do_plot_ec_demod", 
    action="store_true", default=False,
    help="plot HEF C-channel demod")

  parser.add_option("-a", "--plot_avgs", dest="do_plot_avgs", 
    action="store_true", default=False,
    help="plot avg data [default: %default]")

  parser.add_option("--plot_scatter", dest="do_plot_scatter", 
    action="store_true", default=False,
    help="scatter plot A & B self potentials and ocean signal")

  parser.add_option("-n", "--ncat", dest="ncat",
    type="int", metavar="N", default=5, 
    help="number of half cycles for raw plots [default: %default]")

  parser.add_option("-l", "--limEC", dest="limEC",
    type="int", metavar="N", default=0, 
    help="Limit to count_EC per file for raw plots [default: %default]")

  parser.add_option("-s", "--tskip", dest="tskip", 
    type="float", metavar="T", default=15.0, 
    help="skip T first seconds of EF and CAL data [default: %default]")

  parser.add_option("-w", "--twait", dest="twait", 
    type="float", metavar="T", default=15.0, 
    help="wait T seconds after cal state change [default: %default]")

  parser.add_option("-p","--order", dest="fit_order", 
    type="int",  metavar="N", default=1, 
    help="polynomial fitting order for ef_raw residual plots [default: %default]")

  parser.add_option("-g","--odir", dest="odir", metavar='DIR',
    default=None, help="output directory for gz files")

  parser.add_option("-y", "--ylim_ocean", dest="ylim_ocean", 
    type="float", metavar="Y", default=20.0, 
    help="y-limit on hef-ocean plots [default: %default]")

  parser.add_option('-f',"--fifolen", dest="fifolen", metavar='L', 
    type="int", default=3, 
    help="HEF demod FIFO length [default: %default]")

  parser.add_option("--swap_e1bc", dest="swap_e1bc", 
    action="store_true", default=False,
    help="swap e1b and e1c")

  parser.add_option("--swap_e1ab", dest="swap_e1ab", 
    action="store_true", default=False,
    help="swap e1a and e1b")

  parser.add_option("--swap_e2bc", dest="swap_e2bc", 
    action="store_true", default=False,
    help="swap e2b and e2c")

  parser.add_option("--swap_e2ab", dest="swap_e2ab", 
    action="store_true", default=False,
    help="swap e2a and e2b")

  parser.add_option("--e1ws", dest="e1ws", metavar='WS',
    default="", help="name for water switch attached to e1")

  parser.add_option("--e2ws", dest="e2ws", metavar='WS',
    default="", help="name for water switch attached to e2")

  parser.add_option("--begdate", dest="begdate", metavar='DATE',
    default='1970-01-01_00:00:00',
    help="beginning datetime, yyyy-mm-dd_hh:mm:ss")

  parser.add_option("--enddate", dest="enddate", metavar='DATE',
    default='1970-01-01_00:00:00',
    help="ending datetime, yyyy-mm-dd_hh:mm:ss")

  parser.add_option("--pltnambase", dest="pltnambase",  metavar='NB',
    default=None, help="plot name base")

  parser.add_option("--pdfroot", dest="pdfroot", metavar='DN',
    default='./pdf', help="root name of directory for pdf files")

  (options, args) = parser.parse_args()

  wsnams = ''
  if len(options.e1ws) > 0:
    wsnams = wsnams + ' e1=' + options.e1ws
  if len(options.e2ws) > 0: 
    wsnams = wsnams + ' e2=' + options.e2ws

  mf = int(options.fifolen / 2)
  if mf * 2 + 1 != options.fifolen:
    print('error: hpproc: fifolen=',options.fifolen,'must be odd')
    sys.exit(1)


  META = collections.namedtuple('META',[])
  META.tsamp_mot = 0.025
  META.tsamp_ef  = 0.1024
  META.tskip     = options.tskip
  META.twait     = options.twait
  META.fit_order = options.fit_order

  # HEF_HDR.nwait_cal = None
  # HEF_HDR.twait_cal = None

  # HEF_HDR.nskip_cal = None
  # HEF_HDR.tskip_cal = None
  # HEF_HDR.nskip_ef  = None
  # HEF_HDR.tskip_ef  = None
  # HEF_HDR.tavg_ef   = None
  # HEF_HDR.tavg_cal  = None
  # HEF_HDR.tavg_mot  = None

  if options.interactive:
    try:
      __IPYTHON__
    except:
      print('cannot use "--interactive" unless running in ipython')
      sys.exit()

  if len(args) < 1:
    print('not enough args')
    parser.print_help()
    sys.exit()

  begsecs = decode_time(options.begdate)
  if begsecs == None:
    print('cannot decode time: begdate=',options.begdate)
    sys.exit(1)

  endsecs = decode_time(options.enddate)
  if endsecs == None:
    print('cannot decode time: enddate=',options.enddate)
    sys.exit(1)

  if options.verbose:
    print('begdate=',options.begdate, 'begsecs=',begsecs)
    print('enddate=',options.enddate, 'endsecs=',endsecs)
    print('tskip=',META.tskip,'seconds')
    print('fit_order=',META.fit_order)
    print('twait=',META.twait,'seconds')

  HEF_HDR.hef_iscan = -1
  HEF_HDR.hef_nscan = 0
  hef_ntoss = 0

  if options.odir != None:
    mymkdir(options.odir)

  ifiles = []
  for arg in args:
    for ifile in glob.glob(arg):
      ifiles.append(ifile)
    for ifile in glob.glob(arg + '.gz'):
      ifiles.append(ifile)

  if options.verbose:
    print('len(ifiles)=',len(ifiles))

  for ifile in sorted(ifiles):
    print('ifile=',ifile)

    ifd = gz_open(ifile)

    leafname = os.path.basename(ifile)

    i = leafname.find('.gz')
    if i > 0:
      leafname = leafname[:i]

    INFO = get_hpies_info(leafname)

    if options.pltnambase == None:
      pltnambase = leafname
    else:
      pltnambase = options.pltnambase
      

    if options.verbose:
      print('runid=',INFO.runid)

    count_raw_mot_cat = 0
    count_raw_cal_cat = 0
    count_raw_ef_cat = 0

    initbufs()
    initcounts()
    init_HEF_AVG()
    init_HEF_FIFO()
    init_HEF_DEMOD()
    init_CAL_DEMOD()
    init_MOT_CUR()

    for linein in ifd:
      lineno += 1

      if len(linein) < 19:
        print('hpproc: short line: lineno=',lineno,linein)
        count_skipped += 1
        continue

      if linein.find('hprsn:') >= 0:
        print('hpproc: lineno=',lineno,linein)
        count_skipped += 1
        continue

      # use only those lines with iso8601 time followed by #N_
      # if linein.find('T') != 8 or linein.find(' #') != 15 or linein.find('_') != 18:
      if linein[8] != 'T' or linein[15:17] != ' #' or linein[18] != '_':
        count_skipped += 1
        continue

      while linein[-1] == '\r' or linein[-1] == '\n':
        linein = linein[0:-1]

      iso = linein[:15]
      line = linein[16:]

      secs = timegm(strptime(iso,'%Y%m%dT%H%M%S'))
      if secs < begsecs:
        continue
      if endsecs > 0 and secs > endsecs:
        continue

      count_lineok += 1

      if options.verbose and count_lineok == 1:
        print('first linein=',linein)

      # time of day from RSN
      if line[0:3] == '#2_':
        if check_crc(line):
          split_tod(line)
          count_tod += 1

      # HEF data
      if line[0:3] == '#3_':
        if check_crc(line):
          split_hef(line)
          if hef_split != None:
            decode_hef_hdr()
            check_hef_time()
            append_hef_data()   # also plots raw ef
            decode_cal_status() # also plots raw cal
            decode_mot_status() # also plots raw mot
            append_compass()
            count_hef += 1

      # IES data from its AUX2 port
      if line[0:3] == '#5_':
        if check_crc(line):
          if check_crc_aux(line):
            if split_aux(line) == True:
              append_aux()
              count_aux += 1
            if options.verbose:
              print('hpproc: AUX2:',linein)

    ifd.close()
    if options.verbose:
      print('closed ifd')

    compute_compass()
    compute_hef_avg()

    if options.odir != None:
      writegz(pltnambase)

    if options.do_extra:
      printcounts()
      printsizes()

    if count_lineok == 0:
      print('no data found')
      sys.exit(1)

    if options.do_plot_overlay:
      plot_mot_overlay()
      plot_hef_overlay()
      plot_cal()
      plot_compass()
      plot_tt()
      plot_pres_temp()
      plot_bliley()
      plt.show()

    if options.do_plot_avgs:
      plot_ef_avg_std()

    if options.do_plot_scatter:
      plot_scatter()

    if options.do_plot_hef_demod:
      plot_hef_demod()

    if options.do_plot_cal_demod:
      plot_cal_demod()

    if options.do_plot_mot_cur_dur:
      plot_mot_cur_dur()

    if count_removed > 0:
      print('count_removed=',count_removed)

    print('ifile=',ifile,'finished')
