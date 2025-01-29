#!/usr/bin/env python2
# hefcf2.py -- read HEF raw CF2 files

"""
./hefcf2.py run4 -r -t -c     -o /data/okmc/hefcf2
./hefcf2.py run3 -f 5 -p -H 550
./hefcf2.py run4 -r -p -H 550 -o /data/okmc/hefcf2
./hefcf2.py H1 -f 100 -r      -o /data/okmc/hefcf2

./hefcf2.py H2 -f 100 -r      -o /data/okmc/hefcf2
./hefcf2.py H3 -f 400 -r      -o /data/okmc/hefcf2

for h in 1 2 3 4 5; do
for f in 10 50 100 150 200 250 300 350; do 
./hefcf2.py H$h -f$f -r -o /data/okmc/hefcf2
done
done

Aug 1, 2014:
./hefcf2.py H1 100 -r -o /data/okmc/hefcf2

"""

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
from matplotlib.ticker import MaxNLocator

from crclib import crc3kerm

def fixlims(ax):
  ax.set_xlim( larger_axlim( ax.get_xlim() ) )
  ax.set_ylim( larger_axlim( ax.get_ylim() ) )


def fixdateticks(ax):
  ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m-%d %H:%M'))
# ax.xaxis.set_major_locator(MaxNLocator(6))
# above is crude, below is nicer
  ax.xaxis.set_major_locator(mdates.AutoDateLocator())
  ticks = ax.get_xticks()
  n = len(ticks)//4
  ax.set_xticks(ticks[::n])

# lo,hi = ax.get_xlim()
# if hi-lo < 0.1:
#   ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=10),MaxNLocator(5))
# elif hi-lo < 0.4:
#   ax.xaxis.set_major_locator(mdates.HourLocator(interval=1))
# else:
#   ax.xaxis.set_major_locator(mdates.HourLocator(interval=4))
# ax.xaxis.set_major_locator(mdates.AutoDateLocator(minticks=3,maxticks=5))
# ax.xaxis.set_major_locator(mdates.AutoDateLocator())

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
# cur1 = 4300.0 / 5000.0 *adc[ja1] * vref / 65536 * 1000 # mA
# cur2 = 4300.0 / 5000.0 *adc[ja2] * vref / 65536 * 1000 # mA
  cur1 = 5000.0 / 5000.0 *adc[ja1] * vref / 65536 * 1000 # mA
  cur2 = 5000.0 / 5000.0 *adc[ja2] * vref / 65536 * 1000 # mA

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

  fig.suptitle(comment)

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

  fig.suptitle(comment)

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
# cur1 = 4300.0 / 5000.0 *adc[ja1] * vref / 65536 * 1000 # mA
# cur2 = 4300.0 / 5000.0 *adc[ja2] * vref / 65536 * 1000 # mA
  cur1 = 5000.0 / 5000.0 *adc[ja1] * vref / 65536 * 1000 # mA
  cur2 = 5000.0 / 5000.0 *adc[ja2] * vref / 65536 * 1000 # mA

  if hicur > 0 and max(cur1) > hicur:
    nhicur = 5
    ishicur = True
  else:
    ishicur = False

  global isubsamp

  if subsamp > 0 and isubsamp >= subsamp:
    isubsamp = 0
    doss = True
  else:
    isubsamp += 1
    doss = False

  if nhicur == 0 and doss == False:
    return

  if nhicur > 0:
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

def dailyopen(opref, pyt, s):
  global uxtfn, ofp
  donew = False
  if ofp == None:
    print(pyt)
    uxt = timegm(pyt.timetuple())
    uxt0 = int(uxt/86400)*86400
    donew = True
  elif not (s.find('#3__DM ')==0 or s.find('#3__DE ')==0 or \
            s.find('#3__DC ')==0 or s.find('#3__Sm ')==0 or \
            s.find('#3__s')==0):
    uxt = timegm(pyt.timetuple())
    uxt0 = int(uxt/86400)*86400
    if uxt0 != uxtfn:
      ofp.close()
      donew = True
  if donew:
    uxtfn = uxt0
    pytfn = datetime(1970,1,1) + timedelta(0,uxtfn)
    ofile = opref + '_' + pytfn.strftime('%Y%m%d')
    print('ofile=',ofile)
    try:
      ofp = open(ofile,'ab')
    except:
      print('cannot open ofile=',ofile)
      sys.exit(1)
    return

def print_rsn_crc(pyt,channo,s0):
  # add channel number and CRC like hprsn.c on STM would do
  s1 = str.format('#{0:d}_{1:s}',channo,s0)
  s2 = str.format('{0:s}*{1:04x}',s1,crc3kerm(s1))
  iso = pyt.strftime('%Y%m%dT%H%M%S')
  s3 = iso + ' ' + s2
  if options.opref:
    dailyopen(options.opref, pyt, s1)
    ofp.write(s3 + '\r\n')
  else:
    print(s3)

def print_rsn_compass(COMP):
  s = str.format('_{0} {1} {2} {3} {4} {5} {6} {1}', \
        'HC03',COMP.secs,COMP.ticks,COMP.heading, \
        COMP.pitch,COMP.roll,COMP.temperature)
  print_rsn_crc(COMP.pyt,3,s)

def load_ies_files():
# global ies_data, ies_hef_uxt, ies_hef_ind
  global IES
  
  count_no_hef_data_available = 0

  # read whole file of HEF data
  try:
    with open(ies_Cfile) as ifp:
      ies_C = ifp.readlines()
    ifp.close()
  except:
    print('load_ies_files: cannot read Cfile=',ies_Cfile)
    sys.exit(1)

  # read whole file of TT data
  try:
    with open(ies_Tfile) as ifp:
      ies_TT = ifp.readlines()
    ifp.close()
  except:
    print('load_ies_files: cannot read Tfile=',ies_Tfile)
    sys.exit(1)

  # read whole file of Pressure & Temperature data
  try:
    with open(ies_Pfile) as ifp:
      ies_PT = ifp.readlines()
    ifp.close()
  except:
    print('load_ies_files: cannot read Pfile=',ies_Pfile)
    sys.exit(1)

  uxtCi = []
  indCh = []
  uxtCh = []
  for i in range(len(ies_C)):
    sp = ies_C[i].split(',')
    appended = False
    if len(sp) == 30:
      sp0 = sp[0].split()
      if len(sp0) == 7:
#       s = '{0:04d}{1:02d}{2:02d}T{3:02d}{4:02d}{5:02d}'.\
#         format(int(sp0[0]),int(sp0[1]),int(sp0[2]),\
#                int(sp0[3]),int(sp0[4]),int(sp0[5]))
        try:
          pyt = time.strptime(sp[0],'%Y %m %d %H %M %S HEF2')
          uxt = timegm(pyt)
          ok = True
        except:
          print('load_ies_files: ies_C: cannot decode sp[0]=',sp[0])
          print('i=',i,'ies_Cfile=',ies_Cfile)
          ok = False
        if ok:
          uxtCi.append(int(uxt / 600) * 600)
          uxtCh.append(int(sp[29]))
          indCh.append(i)
          appended = True
    if not appended:
      if ies_C[i].find('NO HEF DATA AVAILABLE'):
        count_no_hef_data_available += 1
      elif i > 0:
        print('i=',i,'ies_Cfile=',ies_Cfile)
        print('load_ies_files: ies_C: discarded line=',ies_C[i])

  uxtTT = []
  TT1= []
  TT2= []
  TT3= []
  TT4= []
  for i in range(len(ies_TT)):
    sp = ies_TT[i].split()
    if len(sp) == 25:
      hr0 = int(sp[0])
      for j in range(6):
        j4 = j * 4
        uxtTT.append(hr0 * 3600 + j * 600)
        TT1.append(int(sp[j4+1]))
        TT2.append(int(sp[j4+2]))
        TT3.append(int(sp[j4+3]))
        TT4.append(int(sp[j4+4]))
    elif i>1:
      print('i=',i,'ies_Tfile=',ies_Tfile)
      print('load_ies_files: ies_TT: discarded line=',ies_TT[i])


  uxtPT = []
  PR = []
  TE = []
  for i in range(len(ies_PT)):
    sp = ies_PT[i].split()
    if len(sp) == 13:
      hr0 = int(sp[0])
      for j in range(6):
        j2 = j * 2
        uxtPT.append(hr0 * 3600 + j * 600)
        PR.append(int(sp[j2+1]))
        TE.append(int(sp[j2+2]))
    elif i > 1:
      print('i=',i,'ies_Tfile=',ies_Tfile)
      print('load_ies_files: ies_PT: discarded line=',ies_PT[i])

  if count_no_hef_data_available>0:
    print('load_ies_files(): count_no_hef_data_available=',\
            count_no_hef_data_available)

  uxtCi = np.array(uxtCi,dtype='int')
  uxtCh = np.array(uxtCh,dtype='int')
  indCh = np.array(indCh,dtype='int')

  uxtCoefRef = np.mean(uxtCh)
  uxtCoef = np.polyfit(uxtCh-uxtCoefRef,uxtCh-uxtCi,2)
# print('uxtCoef=',uxtCoef)
  uxtChAdj = uxtCh - np.polyval(uxtCoef,uxtCh-uxtCoefRef)

# print('uxtTT:    len=',len(uxtTT),'[0]=',uxtTT[0],'[-1]=',uxtTT[-1])
# print('uxtPT:    len=',len(uxtPT),'[0]=',uxtPT[0],'[-1]=',uxtPT[-1])
# print('uxtCi:    len=',len(uxtCi),'[0]=',uxtCi[0],'[-1]=',uxtCi[-1])
# print('uxtCh:    len=',len(uxtCh),'[0]=',uxtCh[0],'[-1]=',uxtCh[-1])
# print('uxtChAdj: len=',len(uxtChAdj),'[0]=',uxtChAdj[0],'[-1]=',uxtChAdj[-1])

# print('uxtChAdj[0]  - uxtCi[0] =',uxtChAdj[0]  - uxtCi[0])
# print('uxtChAdj[-1] - uxtCi[-1]=',uxtChAdj[-1] - uxtCi[-1])

  IES = collections.namedtuple('IES',[])

  IES.uxtCi = np.array(uxtCi)
  IES.uxtCh = np.array(uxtCh)
  IES.uxtChAdj = np.array(uxtChAdj)
  IES.indCh = np.array(indCh)

  IES.uxtTT = np.array(uxtTT)
  IES.TT1 = np.array(TT1)
  IES.TT2 = np.array(TT2)
  IES.TT3 = np.array(TT3)
  IES.TT4 = np.array(TT4)

  IES.uxtPT = np.array(uxtPT)
  IES.PR = np.array(PR)
  IES.TE = np.array(TE)


def pick_ies_data(uxt):
  x = np.abs(IES.uxtCh - uxt)
  jC = x.argmin() # index of IES.uxtC* variables
  if x[jC] < 1:
    uxtCi = IES.uxtCi[jC]
  else:
    print('pick_ies_data() bad: jC=',jC,'x[jC]=',x[jC],
      'uxt=',uxt,'IES.uxtCh[jC]=',IES.uxtCh[jC])
    uxtCi = 0
    TT1 = 0
    TT2 = 0
    TT3 = 0
    TT4 = 0
    PR  = 0
    TE  = 0

  if uxtCi > 0:
    x = np.abs(IES.uxtTT - uxtCi)
    jT = x.argmin()
    if x[jT] < 1:
      TT1 = IES.TT1[jT]
      TT2 = IES.TT2[jT]
      TT3 = IES.TT3[jT]
      TT4 = IES.TT4[jT]
    else:
      print('pick_ies_data() bad: jT=',jT,
        'uxt=',uxt,'uxt-uxtCi=',uxt-uxtCi,
        'uxtCi-IES.uxtTT[jT]=',uxtCi-IES.uxtTT[jT])
      TT1 = 0
      TT2 = 0
      TT3 = 0
      TT4 = 0

    x = np.abs(IES.uxtPT - uxtCi)
    jP = x.argmin()
    if x[jP] < 1:
      PR  = IES.PR[jP]
      TE  = IES.TE[jP]
    else:
      print('pick_ies_data() bad: jP=',jP,
        'uxt=',uxt,'uxt-uxtCi=',uxt-uxtCi,
        'uxtCi-IES.uxtPT[jP]=',uxtCi-IES.uxtPT[jP])
      PR  = 0
      TE  = 0

  # fake values for Bliley not in standard IES files
  BlileyT = TE
  BlileyF = 4000000.000

  s = 'AUX,{0:d},04,{1:06d},{2:06d},{3:06d},{4:06d},{5:07d},{6:06d},{7:06d},{8:012.3f}'.\
    format(uxtCi,TT1,TT2,TT3,TT4,PR,TE,BlileyT,BlileyF)
  crc = crc3kerm(s)
  s = '\\r\\r' + s + ',{0:04X}'.format(crc)
  return s

def print_rsn_aux(uxt):
  # To show time of data transfer from HPIES-OKMC HEF to IES
  uxt0 = uxt - uxt % 600
  iesstr = pick_ies_data(uxt)
  s = '{0},{1}'.format(iesstr,uxt)
  pyt = datetime(1970,1,1,0,0,0) + timedelta(0,uxt)
  print_rsn_crc(pyt, 5, s)

def print_rsn_hef_hdr(HE):
  # simulate HEF
  s = str.format('_{0} {1} {2} {3} {4} {5} {6} {7} {8} {9} {10} {6}', \
                 HE.typ,HE.MEC,HE.AB,HE.ibeg,HE.iend,HE.hcno,\
                 HE.secs,HE.ticks,HE.navg_mot,HE.navg_ef,HE.navg_cal)
  print_rsn_crc(HE.pyt,3,s)

def print_rsn_state(pyt,styp, state):

  # count state changes
  count = 0
  ltyp = styp.lower()
  si = state[0]
  for i in range(len(state)):
    sp = si
    si = state[i]
    if i==0 or si != sp:
      count += 1

  # write state hdr
  s = str.format('_S{0} {1} {2} {3}',styp.upper(),0,len(state),count)
  print_rsn_crc(pyt,3,s)

  # write state changes
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

##########################################################
# main

parser = OptionParser(
  usage="%prog [Options] runid", 
  version="%prog 1.0")

parser.add_option("-v", "--verbose",
  action="store_true", dest="verbose", default=False,
  help="print debug info to stdout")

parser.add_option("-p", "--do_plt_each",
  action="store_true", dest="do_plt_each", default=False,
  help="plot each motor run separately")

parser.add_option("-b", "--nbeg",
  dest="nbeg", default="0", metavar='NBEG',
  help="motor run count to begin motor current plots")

parser.add_option("-f", "--thefile",
  dest="thefile", default="-1", metavar='THEFILE',
  help="single file number to process")

parser.add_option("-H", "--hicur",
  dest="hicur", default="0", metavar='HICUR',
  help="label high current plots with '-hicur'")

parser.add_option("-S", "--subsamp",
  dest="subsamp", default="0", metavar='SUBSAMP',
  help="subsample plots (each)")

parser.add_option("-e", "--nend",
  dest="nend", default="0", metavar='NEND',
  help="motor run count to end motor current plots, zero indicates to end of files")

parser.add_option("-t","--plt_mottim",
  action="store_true", dest="do_plt_mottim", default=False,
  help="plot motor times")

parser.add_option('-c',"--plt_motcur",
  action="store_true", dest="do_plt_maxcur", default=False,
  help="plot motor current")

parser.add_option('-r',"--rsn",
  action="store_true", dest="do_rsn", default=False,
  help="make files like hprtst.py using compact flash files from hefrun.c")

parser.add_option("-o", "--opref",
  dest="opref", default=None, metavar='ODIR',
  help="output director for -r")

parser.add_option('-i',"--interactive",
  action="store_true", dest="interactive", default=False,
  help="interactive plotting -- run using ipython")

(options, args) = parser.parse_args()

nbeg = int(options.nbeg)
nend = int(options.nend)
thefile = int(options.thefile)
hicur = float(options.hicur)
subsamp = int(options.subsamp)
ofp = None


if len(args) < 1:
  print('no runid(s) given')
  parser.print_help()
  sys.exit(1)

for runid in args:

  if runid == 'H1':
    comment = 'OKMC-H1'
    prefix = 'JUN8A'
    idir   = '/data/okmc/recovery/H1/HEF002/DATA'
    ies_Cfile  = '/data/okmc/recovery/H1/IES-176/DATA/C176_2.DAT'
    ies_Tfile  = '/data/okmc/recovery/H1/IES-176/DATA/T176_2.DAT'
    ies_Pfile  = '/data/okmc/recovery/H1/IES-176/DATA/P176_2.DAT'
  elif runid == 'H2':
    comment = 'OKMC-H2'
    prefix = 'JUN08B'
    idir   = '/data/okmc/recovery/H2/HEF005/DATA'
    ies_Cfile  = '/data/okmc/recovery/H2/IES-179/DATA/C179_1.DAT'
    ies_Tfile  = '/data/okmc/recovery/H2/IES-179/DATA/T179_1.DAT'
    ies_Pfile  = '/data/okmc/recovery/H2/IES-179/DATA/P179_1.DAT'
  elif runid == 'H3':
    comment = 'OKMC-H3'
    prefix = 'JUN09A'
    idir   = '/data/okmc/recovery/H3/HEF001/DATA'
    ies_Cfile  = '/data/okmc/recovery/H3/IES-175/DATA/C175_1.DAT'
    ies_Tfile  = '/data/okmc/recovery/H3/IES-175/DATA/T175_1.DAT'
    ies_Pfile  = '/data/okmc/recovery/H3/IES-175/DATA/P175_1.DAT'
  elif runid == 'H4':
    comment = 'OKMC-H4'
    prefix = 'JUN09B'
    idir   = '/data/okmc/recovery/H4/HEF003/DATA'
    ies_Cfile  = '/data/okmc/recovery/H4/IES-177/DATA/C177_1.DAT'
    ies_Tfile  = '/data/okmc/recovery/H4/IES-177/DATA/T177_1.DAT'
    ies_Pfile  = '/data/okmc/recovery/H4/IES-177/DATA/P177_1.DAT'
  elif runid == 'H5':
    comment = 'OKMC-H5'
    prefix = 'JUN11A'
    idir   = '/data/okmc/recovery/H5/HEF004/DATA'
    ies_Cfile  = '/data/okmc/recovery/H5/IES-178/DATA/C178_1.DAT'
    ies_Tfile  = '/data/okmc/recovery/H5/IES-178/DATA/T178_1.DAT'
    ies_Pfile  = '/data/okmc/recovery/H5/IES-178/DATA/P178_1.DAT'
  elif runid == 'run1':
    comment = 'm1: old on bench, m2: not connected'
    prefix = 'FEB12A'
    idir   = '/data/hpies/5kpsi/run1'
  # myzone = timezone('US/Pacific')
  elif runid == 'run2':
    comment = 'm1: new 5000 psig, m2=old on bench'
    prefix = '12FEBB'
    idir   = '/data/hpies/5kpsi/run2'
    if nbeg == 0:
      nbeg = 50
  elif runid == 'run3':
    comment = 'm1: new 5000 psig, m2=old on bench'
    prefix = 'FEB18A'
    idir   = '/data/hpies/5kpsi/run3'
  elif runid == 'run4':
    comment = 'm1: new 5000 psig, m2: old on bench'
    prefix = 'FEB24A'
    idir   = '/data/hpies/5kpsi/run4'
  elif runid == 'run5':
    comment = 'm1: new 5000 -> 0 psig, m2: old on bench'
    prefix = 'MAR03A'
    idir   = '/data/hpies/5kpsi/run5'
    if nbeg == 0:
      nbeg = 2
  # 2014-03-03 -- run5: short with several pressures
  # Time     psig
  # 11:43 -- 5000
  # 11:45 -- 4000 
  # 11:50 -- 3000
  # 11:55 -- 2000
  # 12:00 -- 1000
  # 12:05 --  500
  # 12:10 --    0
  # 12:30 -- stop
  else:
    print('unknown runid=',runid)
    sys.exit(1)

  init_motor_stats()

  count_HE = 0
  count_comp = 0
  count_ef = 0
  count_cal = 0
  count_mot = 0
  count_mot_per_file = 0
  count_dcs = 0
  count_None = 0

  load_ies_files()

  for ifile in sorted(os.listdir(idir)):
    tok = ifile.split('.')
    if len(tok) < 2 or len(tok) > 3:
  #   print('skipped ifile=',ifile)
      continue
    if tok[0] != prefix:
  #   print('skipped ifile=',ifile)
      continue
    if len(tok[1]) != 3:
  #   print('skipped ifile=',ifile)
      continue
    try:
      filenostr = tok[1]
      fileno = int(filenostr)
    except:
  #   print('skipped ifile=',ifile)
      continue
    if thefile >= 0 and fileno != thefile:
      continue
    ipath = idir + '/' + ifile
    print('ipath=',ipath)

    if len(tok) == 3: 
      if tok[2] != 'gz':
        print('unknown file=',file)
        continue
      try:
        ifp = gzip.open(ipath,'rb')
      except:
        ifp = None

    if len(tok) == 2:
      try:
        ifp = open(ipath,'rb')
      except:
        ifp = None

    if ifp == None:
      print('cannot open ipath=',ipath)
      count_None += 1
      continue

    """
    if options.odir != None:
      dailyopen(None,None,None)
      ofile = str.format('{0:s}/{1:s}-{2:s}',options.odir,runid,filenostr)
      print('ofile=',ofile)
      try:
        ofp = open(ofile,'wt')
      except:
        print('cannot open ofile=',ofile)
        sys.exit(1)
    else:
      ofp = None
    """


    ibuf = ifp.read()
    ifp.close()

    isubsamp = 0
    nhicur = 0

    ii = 0
    while ii < len(ibuf) - 4:
      typ = struct.unpack('>4s',ibuf[ii+0:ii+4])[0]
      if typ == 'HC02':
        count_comp += 1
        COMP = collections.namedtuple('COMP',[])
        COMP.typ = typ
        COMP.secs        = struct.unpack('>I',ibuf[ii+4 :ii+8])[0]
        COMP.ticks       = struct.unpack('>H',ibuf[ii+8 :ii+10])[0]
        COMP.heading     = struct.unpack('>h',ibuf[ii+10:ii+12])[0]
        COMP.pitch       = struct.unpack('>h',ibuf[ii+12:ii+14])[0]
        COMP.roll        = struct.unpack('>h',ibuf[ii+14:ii+16])[0]
        COMP.temperature = struct.unpack('>h',ibuf[ii+16:ii+18])[0]
        COMP.pyt = datetime(1970,1,1,0,0,0) + timedelta(0,COMP.secs)
        ii += 18
  #     print(typ, COMP.secs)
        if options.do_rsn:
          print_rsn_compass(COMP)
      elif typ == 'HE04':
        count_HE += 1
        hdr = ibuf[ii:ii+100]
        ii += 100
        i = hdr.find('\n')
        if i > 0 and i < 96:
          hdr = hdr[0:i]
  #       print('hdr=',hdr)
        else:
          print('no LF in hdr, ii=',ii)
          sys.exit(1)
        a = hdr.split()
        if len(a) != 11:
          print('should have 11 tokens in HE04 hdr, ii=',ii)
          print(a)
          sys.exit(1)
        HE = collections.namedtuple('HE',[])
        HE.typ = a[0]
        HE.MEC = a[1]
        HE.AB  = a[2]
        HE.ibeg = int(a[3])
        HE.iend = int(a[4])
        HE.hcno = int(a[5])
  #         count_prev = count_master
  #         count_master = HE.hcno + count_wrap
        HE.secs  = int(a[6])
        HE.ticks = int(a[7])
        HE.navg_mot = int(a[8])
        HE.navg_ef  = int(a[9])
        HE.navg_cal = int(a[10])
        nscan = HE.iend - HE.ibeg

  #     fmt = '%Y-%m-%d %H:%M:%S %Z%z'
  #     pyt = datetime(1970,1,1,0,0,0,tzinfo=myzone) + timedelta(0,HE.secs)
        HE.pyt = datetime(1970,1,1,0,0,0) + timedelta(0,HE.secs)

        print_rsn_hef_hdr(HE)

        if HE.MEC == 'f' or HE.MEC == 'r':
          count_mot += 1
          count_mot_per_file += 1
          if nend != 0 and count_mot > nend:
            sys.exit(1)
          nss = int(nscan / HE.navg_mot)
          fmt = str.format('>{0:d}f',nss)
          MOT = collections.namedtuple('MOT',[])
          MOT.navg = HE.navg_mot
          MOT.adc = struct.unpack(fmt,ibuf[ii:ii+nss*4])
          ii += nss * 4      # one channel motor current (float32)
          fmt = str.format('>{0:d}B',nscan)
          MOT.state = struct.unpack(fmt,ibuf[ii:ii+nscan])
          ii += nscan        # motor state (uint8)
          if count_mot >= nbeg:
            if options.do_plt_each:
              plot_motor_each(HE, MOT.adc, MOT.state)
            accum_motor_stats(HE.pyt,HE.AB,MOT.adc,HE.navg_mot, MOT.state)
          if options.do_rsn:
            print_rsn_mot(HE.pyt,MOT)
        elif HE.MEC == 'E':
          count_ef += 1
          nss = int(nscan / HE.navg_ef)
          # EF data
          fmt = str.format('>{0:d}f',4)
          for ia in range(nss):
            avgs = struct.unpack(fmt,ibuf[ii:ii+16])
            ii += 16 # 4 channels of float32
            if options.do_rsn:
              print_rsn_avgs(HE.pyt,'E',ia,avgs,HE.navg_ef)
        elif HE.MEC == 'C':
          count_cal += 1
          nss = int(nscan / HE.navg_cal)
          # cal data
          fmt = str.format('>{0:d}f',6)
          for ia in range(nss):
            avgs = struct.unpack(fmt,ibuf[ii:ii+24])
            ii += 24  # 6 channels ADC (float32)
            if options.do_rsn:
              print_rsn_avgs(HE.pyt,'C',ia,avgs,HE.navg_cal)
          # cal state
          fmt = str.format('>{0:d}B',nscan)
          state = struct.unpack(fmt,ibuf[ii:ii+nscan])
          ii += nscan        # cal_run_state_fast (uint8)
          if options.do_rsn:
            print_rsn_state(HE.pyt,'C',state)
        else:
          print('unknown MEC=',HE.MEC,'ii=',ii)
          sys.exit(1)

  #         print(typ,HE.secs, HE.MEC, nscan, nss)
        
      # HEF data sent to IES in stand-alone HPIES
      # should not see this in HPIES-RSN
      elif typ == 'DCS1':
        count_dcs += 1
        dcs1 = ibuf[ii+4:ii+4+512]
        ii += 4 + 512
        i = dcs1.find(chr(0))
        dcs1 = dcs1[:i]
        dcs1_toks = dcs1.split(',')
        # print('len(dcs1_toks)=',len(dcs1_toks))
        if len(dcs1_toks) != 30:
          print('bad dcs1, ntoks={0}'.format(len(dcs1_toks)))
          dcs1_uxt = None
          sys.exit(1)
        else:
          dcs1_uxt = int(dcs1_toks[29])
        print_rsn_aux(dcs1_uxt)
        
        

      else:
        print('unknown typ=',typ)
        print('  ii=',ii)
        sys.exit(1)

    if options.do_plt_mottim:
      plot_mottim()
    if options.do_plt_maxcur:
      plot_maxcur()

    init_motor_stats()
    count_mot_per_file = 0

  if options.verbose:
    print('runid=',runid)
    print('  count_HE=',count_HE)
    print('  count_mot=',count_mot)
    print('  count_ef=',count_ef)
    print('  count_cal=',count_cal)
    print('  count_comp=',count_comp)
    print('  count_dcs=',count_dcs)
  if count_None>0:
    print('  count_None=',count_None)

  if options.do_plt_each:
    os.system('updateframe.run pdf/each')
  if options.do_plt_maxcur or options.do_plt_mottim:
    os.system('updateframe.run pdf/daily')

  print('')
