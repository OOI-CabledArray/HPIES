#!/usr/bin/env python2
# iescf1.py -- read IES CF1 files made by iesjhd.c

from __future__ import print_function

import sys
import os
from optparse import OptionParser
import collections
import time
from calendar import timegm
from datetime import datetime, timedelta
# from pytz import timezone
import matplotlib.dates as mdates
import math
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import gzip
import struct

from crclib import crc3kerm

def fixlims(ax):
  ax.set_xlim( larger_axlim( ax.get_xlim() ) )
  ax.set_ylim( larger_axlim( ax.get_ylim() ) )

def fixdateticks(ax):
  ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
  lo,hi = ax.get_xlim()
  if hi-lo < 0.4:
    inc = 1
  else:
    inc = 4
  ax.xaxis.set_major_locator(mdates.HourLocator(interval=inc))

def fixthings(ax):
  fixlims(ax)
  fixdateticks(ax)
  
def larger_axlim( axlim ):
    """ argument axlim expects 2-tuple 
        returns slightly larger 2-tuple """
    axmin,axmax = axlim
    axrng = axmax - axmin
    new_min = axmin - 0.02 * axrng
    new_max = axmax + 0.02 * axrng
    return new_min,new_max

def average(arr, n):
    end =  n * int(len(arr)/n)
    return numpy.mean(arr[:end].reshape(-1, n), 1)

def accum_motor_stats(pyt,AB,adc,navg,state):
  global maxcur1A, maxcur1B, maxcur2A, maxcur2B
  global secs1A, secs1B, secs2A, secs2B
  global pytA, pytB
  dt = 0.025
  adc   = np.array(adc,dtype=np.float)
  state = np.array(state,dtype=np.int)
  state1 = (state >> 4 ) & 15
  state2 = state & 15
  js1 = np.nonzero((state1 >= 1) & (state1 <= 3))[0]
  js2 = np.nonzero((state2 >= 1) & (state2 <= 3))[0]
  ja1 = np.unique(js1 / navg)
  ja2 = np.unique(js2 / navg)

  vref = 2.5
  cur1 = adc[ja1] * vref / 65536 * 1000 # mA
  cur2 = adc[ja2] * vref / 65536 * 1000 # mA

  if AB == 'a':
    pytA.append(pyt)
    secs1A.append(len(js1) * dt)
    secs2A.append(len(js2) * dt)
    maxcur1A.append(max(cur1))
    maxcur2A.append(max(cur2))
  elif AB == 'b':
    pytB.append(pyt)
    secs1B.append(len(js1) * dt)
    secs2B.append(len(js2) * dt)
    maxcur1B.append(max(cur1))
    maxcur2B.append(max(cur2))
  else:
    print('error: AB not a or b')
    sys.exit(1)
  
def plot_mottim():
  fig = plt.figure(1)
  fig.clf()
  pltnam = 'hefcf2-'+runid+'-'+filenostr+'-mottim'

  ax = fig.add_subplot(221)
  ax.plot(pytA,secs1A,'.')
  fixthings(ax)
  plt.ylabel('m1A, s')

  ttl = str.format('Motor Times {0} {1}',runid,prefix)
  plt.title(ttl)
  
  ax = fig.add_subplot(222)
  ax.plot(pytA,secs2A,'.')
  fixthings(ax)
  plt.ylabel('m2A, s')
  plt.title(pltnam)
  
  ax = fig.add_subplot(223)
  ax.plot(pytB,secs1B,'.')
  fixthings(ax)
  plt.ylabel('m1B, s')

  ax = fig.add_subplot(224)
  ax.plot(pytB,secs2B,'.')
  fixthings(ax)
  plt.ylabel('m2B, s')

  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

  if options.interactive:
    fig.show()
    print('Click figure window for next plot')
    plt.waitforbuttonpress()
  else:
    pdffile = './pdf/daily/'+pltnam+'.pdf'
    print(pdffile)
    plt.savefig(pdffile)

def plot_maxcur():
# fig = plt.figure(1, figsize=(10,8), dpi=300)
  fig = plt.figure(1)
  fig.clf()
  pltnam = 'hefcf2-'+runid+'-'+filenostr+'-maxcur'

  ax = fig.add_subplot(221)
  ax.plot(pytA,maxcur1A,'.')
  fixthings(ax)
  plt.ylabel('m1A, mA')

  ttl = str.format('Max Current {0} {1}',runid,prefix)
  plt.title(ttl)
  
  ax = fig.add_subplot(223)
  ax.plot(pytB,maxcur1B,'.')
  fixthings(ax)
  plt.ylabel('m1B, mA')
  
  ax = fig.add_subplot(222)
  ax.plot(pytA,maxcur2A,'.')
  fixthings(ax)
  plt.ylabel('m2A, mA')
  plt.title(pltnam)
  
  ax = fig.add_subplot(224)
  ax.plot(pytB,maxcur2B,'.')
  fixthings(ax)
  plt.ylabel('m2B, mA')
  
  fig.autofmt_xdate()
  fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

  if options.interactive:
    fig.show()
    print('Click figure window for next plot')
    plt.waitforbuttonpress()
  else:
    pdffile = './pdf/daily/'+pltnam+'.pdf'
    print(pdffile)
#   plt.savefig(pdffile,dpi=fig.dpi)
    plt.savefig(pdffile)

def plot_motor_each(HE, adc, state):
  global nhicur

  pytstr = HE.pyt.strftime('%Y-%m-%d %H:%M:%S')
  ttl = str.format('{0} {1} {2} {3}',runid,ifile,pytstr,HE.AB)
  ttl = ifile + ' ' + pytstr + ' ' + HE.AB.upper()

  dt = 0.025
  adc   = np.array(adc,dtype=np.float)
  state = np.array(state,dtype=np.int)

  state1 = (state >> 4 ) & 15
  state2 = state & 15
  js1 = np.nonzero((state1 >= 1) & (state1 <= 3))[0]
  js2 = np.nonzero((state2 >= 1) & (state2 <= 3))[0]
  ja1 = np.unique(js1 / HE.navg_mot)
  ja2 = np.unique(js2 / HE.navg_mot)

  vref = 2.5
  cur1 = adc[ja1] * vref / 65536 * 1000 # mA
  cur2 = adc[ja2] * vref / 65536 * 1000 # mA

  if hicur > 0 and max(cur1) > hicur:
    nhicur = 5
    ishicur = True
  else:
    ishicur = False

  if nhicur == 0:
    return

  nhicur -= 1


  state1 = state1[js1]
  state2 = state2[js2]

  tcur1 = np.array(range(0,len(cur1)),dtype=np.float) * dt * HE.navg_mot
  tcur2 = np.array(range(0,len(cur2)),dtype=np.float) * dt * HE.navg_mot
  tstate1 = np.array(range(0,len(state1)),dtype=np.float) * dt
  tstate2 = np.array(range(0,len(state2)),dtype=np.float) * dt


# fig = plt.figure(1, figsize=(10,8), dpi=300)
  fig = plt.figure(1)
  fig.clf()
  pltnam = str.format('hefcf2-{0}-{1}-{2:04d}',\
           runid,filenostr,count_mot_per_file)
  
  tmax = 6     # s
  curmax = 700 # mA

  ax = fig.add_subplot(221)
  ax.plot(tcur1,cur1,'.-')
  ax.set_xlim(0,tmax)
  ax.set_ylim(0,curmax)
  ax.grid(True)
  plt.ylabel('Motor1, mA')
  plt.title(ttl)

  ax = fig.add_subplot(222)
  ax.plot(tcur2,cur2,'.-')
  ax.set_xlim(0,tmax)
  ax.set_ylim(0,curmax)
  ax.grid(True)
  plt.ylabel('Motor2, mA')
  plt.title(pltnam)

  ax = fig.add_subplot(223)
  ax.plot(tstate1,state1)
  ax.set_xlim(0,tmax)
  ax.set_ylim(-1,5)
  ax.grid(True)
  plt.ylabel('Motor1 State')
  plt.xlabel('time, s')

  ax = fig.add_subplot(224)
  ax.plot(tstate2,state2)
  ax.set_xlim(0,tmax)
  ax.set_ylim(-1,5)
  ax.grid(True)
  plt.ylabel('Motor2 State')
  plt.xlabel('time, s')

  fig.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0.3, hspace=None)

  if options.interactive:
    fig.show()
    print('Click figure window for next plot')
    plt.waitforbuttonpress()
  else:
    if ishicur:
      pltnam = pltnam + '-hicur'
    pdffile = './pdf/each/'+pltnam+'.pdf'
    print(pdffile)
    plt.savefig(pdffile)

def init_motor_stats():
  global maxcur1A, maxcur1B, maxcur2A, maxcur2B
  global secs1A, secs1B, secs2A, secs2B
  global pytA, pytB
  maxcur1A = []; maxcur1B = []; maxcur2A = []; maxcur2B = []
  secs1A = []; secs1B = []; secs2A = []; secs2B = []
  pytA = []; pytB = []

def print_rsn_crc(pyt,channo,s):
  # add channel number and CRC like hprsn.c on STM would do
  s = str.format('#{0:d}_{1:s}',channo,s)
  s = str.format('{0:s}*{1:04x}',s,crc3kerm(s))
  iso = pyt.strftime('%Y%m%dT%H%M%S ')
  if ofp != None:
    ofp.write(iso + s + '\r\n')
  else:
    print(iso +  s)

def print_rsn_compass(COMP):
  s = str.format('_{0} {1} {2} {3} {4} {5} {6} {1}', \
        'HC03',COMP.secs,COMP.ticks,COMP.heading, \
        COMP.pitch,COMP.roll,COMP.temperature)
  print_rsn_crc(COMP.pyt,3,s)

def print_rsn_hef_hdr(HE):
  # simulate HEF
  s = str.format('_{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {6}', \
                 HE.typ,HE.MEC,HE.AB,HE.ibeg,HE.iend,HE.hcno,\
                 HE.secs,HE.ticks,HE.navg_mot,HE.navg_ef,HE.navg_cal)
  print_rsn_crc(HE.pyt,3,s)

def print_rsn_state(pyt,styp, state):
  # state hdr
  s = str.format('_S{0} {1} {2}',styp.upper(),0,len(state))
  print_rsn_crc(pyt,3,s)

  # state changes
  ltyp = styp.lower()
  si = state[0]
  for i in range(len(state)):
    sp = si
    si = state[i]
    if i==0 or si != sp:
      s = str.format('_S{0} {1} {2}',ltyp,i,si)
      print_rsn_crc(pyt,3,s)

def print_rsn_mot(pyt,MOT):
  # motor current
  na = MOT.navg
  for ia in range(len(MOT.adc)):
    s = str.format('_DM {0} {1:.0f}',(ia+1)*na-1,MOT.adc[ia]*na)
    print_rsn_crc(pyt,3,s)
  print_rsn_state(pyt,'M',MOT.state)

def print_rsn_avgs(pyt,typ, ia,avgs,navg):
  s = str.format('_D{0} {1}',typ,(ia+1)*navg-1)
  for i in range(len(avgs)):
    s = s + str.format(' {0:.0f}',avgs[i]*navg)
  print_rsn_crc(pyt,3,s)

def print_rsn_cal():
  s = str.format('_DC {0}',(ia+1)*navg-1)
  for i in range(len(avgs)):
    s = s + str.format(' {0:.0f}',avgs[i]*navg)
  print_rsn_crc(pyt,3,s)

def print_info_and_exit():
  print('ERROR:')
  print('  ifile=',ifile)
  print('  lineno=',lineno)
  print('  line=',line)
  sys.exit(1)

# main

parser = OptionParser(
  usage="%prog [Options] runid[s]", 
  version="%prog 1.0")

parser.add_option("-v", "--verbose",
  action="store_true", dest="verbose", default=False,
  help="print debug info to stdout")

(options, args) = parser.parse_args()

if len(args) < 1:
  parser.print_help()
  sys.exit()

for runid in args:

  if   runid == 'H1':
    ifile  = '/data/okmc/recovery/H1/IES-176/DATA/C176_2.DAT'
  elif runid == 'H2':
    ifile  = '/data/okmc/recovery/H2/IES-179/DATA/C179_1.DAT'
  elif runid == 'H3':
    ifile  = '/data/okmc/recovery/H3/IES-175/DATA/C175_1.DAT'
  elif runid == 'H4':
    ifile  = '/data/okmc/recovery/H4/IES-177/DATA/C177_1.DAT'
  elif runid == 'H5':
    ifile  = '/data/okmc/recovery/H5/IES-178/DATA/C178_1.DAT'
  else:
    print('unknown runid=',runid)
    sys.exit(1)

  print('ifile=',ifile)

  try:
    ifp = open(ifile,'rb')
  except:
    print('cannot open ifile=',ifile)
    sys.exit(1)

  pyt_ies = []
  uxt_ies = []
  uxt_hef = []

  lineno = 0
  for line in ifp:
    lineno += 1
#   if lineno == 1:
#     continue
    if len(line) <= 1:
      continue
    if line.find('NO HEF DATA AVAILABLE') > -1:
      continue
    toks1 = line.split()
    if len(toks1) != 7:
      print('len(toks1) should be 7, is=',len(toks1))
      print_info_and_exit()
    toks2 = toks1[6].split(',')
    if len(toks2) != 30:
      print('len(toks2) should be 30, is=',len(toks2))
      print_info_and_exit()

    pyt = datetime(int(toks1[0]),int(toks1[1]),int(toks1[2]),\
                   int(toks1[3]),int(toks1[4]),int(toks1[5]))
    pyt_ies.append(pyt)
    uxt_ies.append(timegm(pyt.timetuple()))
    uxt_hef.append(int(toks2[29]))

  uxt_ies = np.array(uxt_ies,dtype=np.double)
  uxt_hef = np.array(uxt_hef,dtype=np.double)

  terr = uxt_hef - uxt_ies
  t = uxt_hef - uxt_hef[0]
  coef = np.polyfit(t,terr,2)
  fit = np.polyval(coef,t)
  res = terr - fit

  print('runid=',runid,', uxt_hef0=',uxt_hef[0],', coef=',coef)

  ifp.close()

  fig = plt.figure()
  fig.clf()
  ax = fig.add_subplot(111)
  ax.plot(pyt_ies, res,'b.')
  fig.autofmt_xdate()
  plt.title(runid+' HEF-IES time diff minus quadratic fit')
  plt.ylabel('s')
  plt.draw()

  # H1 is the worst with residual range of +-3 s
  # H2, H3, H4, H5 are better at +-1.5 s

plt.show()
