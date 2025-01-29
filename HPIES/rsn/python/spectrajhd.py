#!/usr/bin/env python2
# spectrajhd.py

from __future__ import print_function

import os
import sys
import numpy as np
from collections import namedtuple

def error(s):
  print(s)
  sys.exit(1)

def nfft235(nfft_first):
# NFFT = NFFT235(NFFT_MAX) returns an integer which is
# a product of only 2, 3 and 5.  NFFT is equal to or less
# than NFFT_MAX.  NFFT235 is useful with the Fast Fourrier
# Transform

  nfft = 0
  nfft_first = int(nfft_first)
  for ntry in range(nfft_first/2+1):
    ntst = nfft_first - ntry
    n = ntst
    while n > 1:
      if n % 5 == 0:
        n /= 5
      elif n % 3 == 0:
        n /= 3
      elif n % 2 == 0:
        n /= 2
      else:
        break
    if n == 1:
      nfft = ntst
      break
  return nfft

def test_nfft235():
  print('nfft235(12)=',nfft235(12))
  print('nfft235(11)=',nfft235(11))
  print('nfft235(111)=',nfft235(111))
  print('nfft235(1024)=',nfft235(1024))
  print('nfft235(1023)=',nfft235(1023))

def autospec(x,nfft,fsamp,dflag):
# [s,f,npieces] = autospec(x,nfft,fsamp,dflag)
# x = input data must be real and longer than nfft
# nfft = number of points in each FFT'd piece.
# fsamp = sample rate in Hz
# dflag = 'none','mean','linear', detrending applied per piece
# s = autospectrum normalized so that integral of s(f) == std(x)^2
# f = frequency Hz corresponding to s
# npieces = number of pieces
# input data 50% overlapped
# hanning window applied to each piece after dflag operation

  if nfft <= 0:
    error('nfft must be positive')
  if nfft % 2 != 0:
    error('nfft must be even')
  if not np.isreal(x.all()):
    error('x must be real')
  if x.ndim != 1:
    error('x must be one-dimensional')
  if fsamp <= 0:
    error('fsamp must be positive')

  # Note difference from hanning.m.
  # This only spreads line spectra to the two adjacent frequency bins.
  # and normalize the window
  wind = (1 - np.cos(2.0*np.pi*np.arange(nfft)/float(nfft)))/2.0
  wind /= np.sqrt(sum(np.power(wind,2)))

  noverlap = nfft/2
  nslide = nfft - noverlap
  npieces = int((len(x)-noverlap)/nslide)

  spec = np.zeros(nfft/2)
  kpiece = 0
  for ipiece in range(npieces):
    xx = x[np.arange(nfft)+ipiece*nslide]

    if not np.isfinite(xx.all()):
      print('autospec(): not all x finite: skipped ipiece={0}'.format(ipiece))
      continue

    kpiece += 1

    if dflag == 'mean':
      xx = xx - np.mean(xx)
    elif dflag == 'linear':
      i = np.arange(nfft) - nfft / 2 #  # zero mean
      coefs = np.polyfit(i,xx,1)
      fit = np.polyval(coefs,i)
      xx -= fit
      del i
      del coefs
      del fit
    elif dflag == 'none':
      pass
    else:
      error('dflag={} is unknown'.format(dflag))

    xx *= wind

    y = np.fft.fft(xx,nfft)
    spec += np.power(np.absolute(y[np.arange(nfft/2)]),2)

  if kpiece == 0:
    print('autospec: no good pieces -- output all zero')
  else:
    spec *= 2.0 / fsamp / kpiece
  
  SPEC = namedtuple('SPEC',[])
  SPEC.spec  = spec
  SPEC.df    = fsamp / float(nfft)
  SPEC.freq  = SPEC.df * np.arange(nfft/2)
  SPEC.intgl = sum(SPEC.spec) * SPEC.df
  SPEC.var   = np.power(np.std(x),2)
  SPEC.kpie  = kpiece
  SPEC.npie  = npieces
  SPEC.dflag = dflag
  SPEC.nfft  = nfft 
  SPEC.fsamp = fsamp 

  return SPEC

def test_autospec():
  nfft  = 150
  fsamp = 1.0
  tsamp = 1.0 / fsamp
  t = np.arange(nfft*3) * tsamp
  freq  = fsamp / 8.0
  freq  = fsamp / nfft * 3
  x = np.sin(np.pi * 2.0 * freq * t)
  S = autospec(x, nfft, fsamp, 'none')

# print('x=',x)
# print('spec=',S.spec)
  print('   len(x)=',len(x))
  print('len(spec)=',len(S.spec))
  print('len(freq)=',len(S.freq))
  print('      var=',S.var)
  print('    intgl=',S.intgl)
  print('       df=',S.df)
  print('     kpie=',S.kpie)
  print('     npie=',S.npie)


def twospec(x,y,nfft,fsamp):

#
# see spec_rxy_pxy.m
#  rxy = (cxy.^2 + qxy.^2) ./ (axx.*ayy) # coherence squared
#  pxy = -atan2(qxy,cxy)                 # phase
#
# see spec_rotary.m
#  acw = (axx+ayy)/2 - qxy               # anticlockwise
#  cw  = (axx+ayy)/2 + qxy               # clockwise
#
# see spec_rot_coh.m
#  gccwcw = complex(axx-ayy, 2*cxy) ./ sqrt((axx+ayy).^2 - (2*qxy).^2)
#  rotary_coherence_amplitude = abs(gccwcw)
#  rotary_coherence_phase = atan2(2 * cxy, axx - ayy)
#
# complex transfer function:
#  txy = qxy ./ axx

  if fsamp <= 0.0:
    error('fsamp must be positive')
  if nfft <= 0:
    error('nfft must be positive')
  if nfft % 2 != 0:
    error('nfft must be even')
  if not np.isreal(x.all()):
    error('x must be real')
  if not np.isreal(y.all()):
    error('y must be real')
  if x.ndim != 1:
    error('x must be one-dimensional')
  if y.ndim != 1:
    error('y must be one-dimensional')
  if len(x) != len(y):
    error('x & y must be same length')
  if len(x) < nfft:
    error('length of x & y must not be less than nfft')

  S2 = namedtuple('S2',[])

  # normalized Hanning window
  wind = (1.0 - np.cos(2.0*np.pi*np.arange(nfft)/float(nfft)))/2.0
  wind /= np.sqrt(sum(np.power(wind,2)))

  # 50% overlapped
  noverlap = nfft/2
  nslide = nfft - noverlap
  S2.npie = int((len(x)-noverlap)/nslide)

  # skips zero and nyquist
  kk = np.arange(1,nfft/2)
  jj = nfft - kk

  S2.axx = np.zeros(len(kk))
  S2.ayy = np.zeros(len(kk))
  S2.cxy = np.zeros(len(kk))
  S2.qxy = np.zeros(len(kk))

  S2.npie_ng = 0
  S2.npie_ok = 0

  for ipie in range(S2.npie):

    i = np.arange(nfft) + (ipie-1)*nslide
    if (not np.isfinite(x[i].all())) | (not np.isfinite(y[i].all())):
      S2.npie_ng += 1
      continue
    else:
      S2.npie_ok += 1

#   xi = x[i]
#   yi = y[i]
    xi = (x[i] - np.mean(x[i])) * wind
    yi = (y[i] - np.mean(y[i])) * wind

    z = xi + 1j * yi

    Z = np.fft.fft(z,nfft)

    # cos and sin coefficients
    xc =  np.real(Z[kk]) + np.real(Z[jj])
    xs = -np.imag(Z[kk]) + np.imag(Z[jj])
    yc =  np.imag(Z[kk]) + np.imag(Z[jj])
    ys =  np.real(Z[kk]) - np.real(Z[jj])

    S2.axx += (xc*xc + xs*xs)
    S2.ayy += (yc*yc + ys*ys)
    S2.cxy += (xc*yc + xs*ys)
    S2.qxy += (yc*xs - xc*ys)

  if S2.npie_ok > 0:
    S2.axx *= 0.5 / fsamp / S2.npie_ok
    S2.ayy *= 0.5 / fsamp / S2.npie_ok
    S2.cxy *= 0.5 / fsamp / S2.npie_ok # in-phase
    S2.qxy *= 0.5 / fsamp / S2.npie_ok # quadrature

  S2.df = fsamp / float(nfft)
  S2.freq = S2.df * kk

  S2.intgl_axx = sum(S2.axx) * S2.df
  S2.intgl_ayy = sum(S2.ayy) * S2.df
  S2.varx = np.power(np.std(x),2)
  S2.vary = np.power(np.std(y),2)

  # coherence squared and phase
  # i = np.nonzero(np.logical_and(S2.axx>=0,S2.ayy>=0))[0]
  # S2.rxy = (S2.cxy[i]*S2.cxy[i] + S2.qxy[i]*S2.qxy[i]) / (S2.axx[i]*S2.ayy[i])
  # S2.pxy = -np.arctan2(S2.qxy[i],S2.cxy[i])
  S2.rxy = (S2.cxy*S2.cxy + S2.qxy*S2.qxy) / (S2.axx*S2.ayy)
  S2.pxy = -np.arctan2(S2.qxy,S2.cxy)

  # rotatary spectra
  S2.acw = (S2.axx+S2.ayy)/2.0 - S2.qxy  # anticlockwise
  S2.cw  = (S2.axx+S2.ayy)/2.0 + S2.qxy  # clockwise

  # rotary coherence
  gccwcw = (S2.axx-S2.ayy + 2j * S2.cxy) / \
    np.sqrt((S2.axx+S2.ayy)*(S2.axx+S2.ayy) - (2.0*S2.qxy)*(2.0*S2.qxy))
  S2.rot_coh_amp = np.absolute(gccwcw)
  S2.rot_coh_ph = np.arctan2(2.0 * S2.cxy, S2.axx - S2.ayy)

  if(S2.npie_ok == 0):
    print('warning: no good pieces.  ng='.S2.npie_ng)

  return S2


def rotcoh(u1, v1, u2, v2, nfft, fsamp):



# Input:
#   u1, v1 -- first  velocity series, u=east & v=north components
#   u2, v2 -- second velocity series
#   nfft   -- number of points in each FFT
#   fsamp  -- sample frequency, Hz
#
# Output:
#   F                 -- frequency, Hz
#   S11_cw, S11_ccw   -- rotary spectrum of u1, v1
#   S22_cw, S22_ccw   -- rotary spectrum of u2, v2
#   S12_cw, S12_ccw   -- inner rotary cross spectrum of u1, v1, u2, v2
#   gam2_cw, gam2_ccw -- inner rotary coherence squared
#   phi_cw, phi_ccw   -- inner rotary coherence phase
#   lam2_1ccw_2cw, lam2_1cw_2ccw -- outer rotary coherence squared
#   psi_1ccw_2cw, psi_1cw_2ccw   -- outer rotary coherence phase
#
# Where:
#   inner -- rotating in same direction
#   outer -- rotating in opposite directions
#   ccw -- counter clockwise
#   cw  -- clockwise
#
# Equations from section 5.8.8.1 of
# "Data Analysis Methods in Physical Oceanography"
# by William J. Emery and Richard E. Thomsom,
# 1997, 1998 & 2001 Elsevier B.V., Amsterdam, The Netherlands

  if nfft % 2 != 0:
    error('nfft must be even')

  noverlap = nfft/2
  nslide   = nfft - noverlap
  npie     = int((len(u1)-noverlap)/nslide)
  npie_ok = 0
  npie_ng = 0

  kk = np.arange(1,nfft/2)
  jj = nfft - kk
  df = fsamp / nfft
  F = df * kk

  Sp12 = np.tile(0.0+0.0j,len(F))
  Sn12 = np.tile(0.0+0.0j,len(F))
  Sp11 = np.tile(0.0+0.0j,len(F))
  Sn11 = np.tile(0.0+0.0j,len(F))
  Sp22 = np.tile(0.0+0.0j,len(F))
  Sn22 = np.tile(0.0+0.0j,len(F))

  gam2_p12_1 = np.zeros(len(F))
  gam2_p12_2 = np.zeros(len(F))
  gam2_p12_3 = np.zeros(len(F))
  gam2_p12_4 = np.zeros(len(F))
  gam2_p12 = np.tile(np.nan,(1,len(F)))
  phi_p12 = np.tile(np.nan,(1,len(F)))

  gam2_n12_1 = np.zeros(len(F))
  gam2_n12_2 = np.zeros(len(F))
  gam2_n12_3 = np.zeros(len(F))
  gam2_n12_4 = np.zeros(len(F))
  gam2_n12 = np.tile(np.nan,(1,len(F)))
  phi_n12 = np.tile(np.nan,(1,len(F)))

  Yp12 = np.zeros(len(F))
  Yn12 = np.zeros(len(F))

  lam2_p12_1 = np.zeros(len(F))
  lam2_p12_2 = np.zeros(len(F))
  lam2_p12_3 = np.zeros(len(F))
  lam2_p12_4 = np.zeros(len(F))
  lam2_p12_5 = np.zeros(len(F))
  lam2_p12 = np.tile(np.nan,(1,len(F)))

  lam2_n12_1 = np.zeros(len(F))
  lam2_n12_2 = np.zeros(len(F))
  lam2_n12_3 = np.zeros(len(F))
  lam2_n12_4 = np.zeros(len(F))
  lam2_n12_5 = np.zeros(len(F))
  lam2_n12 = np.tile(np.nan,(1,len(F)))

  psi_p12_1 = np.zeros(len(F))
  psi_p12_2 = np.zeros(len(F))
  psi_p12 = np.tile(np.nan,(1,len(F)))

  psi_n12_1 = np.zeros(len(F))
  psi_n12_2 = np.zeros(len(F))
  psi_n12 = np.tile(np.nan,(1,len(F)))

  # normalized hanning window
  window = (1.0 - np.cos(2.0*np.pi*np.arange(nfft)/nfft))/2.0
  window /= np.sqrt(np.mean(window*window))


  for ipie in range(npie):

    j = np.arange(nfft) + ipie*nslide

    u1j = u1[j]
    v1j = v1[j]
    u2j = u2[j]
    v2j = v2[j]

    all_ok = np.isfinite(u1j.all()) & \
             np.isfinite(v1j.all()) & \
             np.isfinite(u2j.all()) & \
             np.isfinite(v2j.all())

    if not all_ok:
      npie_ng += 1
  

    if all_ok:

      npie_ok += 1

      u1j = (u1j - np.mean(u1j)) * window
      v1j = (v1j - np.mean(v1j)) * window
      u2j = (u2j - np.mean(u2j)) * window
      v2j = (v2j - np.mean(v2j)) * window

      w1j = u1j + 1j*v1j
      w2j = u2j + 1j*v2j

      Z1 = np.fft.fft(w1j,nfft) / nfft
      Z2 = np.fft.fft(w2j,nfft) / nfft

      # cos and sin coefficients
      Bu1 =  np.real(Z1[kk]) + np.real(Z1[jj])
      Au1 = -np.imag(Z1[kk]) + np.imag(Z1[jj])
      Bv1 =  np.imag(Z1[kk]) + np.imag(Z1[jj])
      Av1 =  np.real(Z1[kk]) - np.real(Z1[jj])

      Bu2 =  np.real(Z2[kk]) + np.real(Z2[jj])
      Au2 = -np.imag(Z2[kk]) + np.imag(Z2[jj])
      Bv2 =  np.imag(Z2[kk]) + np.imag(Z2[jj])
      Av2 =  np.real(Z2[kk]) - np.real(Z2[jj])

      # amplitudes for positive and negative frequencies 
      Ap1 = np.sqrt(np.power(Bv1+Au1,2) + np.power(Av1-Bu1,2))/2.0
      An1 = np.sqrt(np.power(Bv1-Au1,2) + np.power(Av1+Bu1,2))/2.0
      Ap2 = np.sqrt(np.power(Bv2+Au2,2) + np.power(Av2-Bu2,2))/2.0
      An2 = np.sqrt(np.power(Bv2-Au2,2) + np.power(Av2+Bu2,2))/2.0

      # phases
   #  Pp1 = np.arctan2(Au1-Bu1,Au1+Bv1) # as per 5.8.47a
      Pp1 = np.arctan2(Av1-Bu1,Au1+Bv1) # corrected
      Pn1 = np.arctan2(Bu1+Av1,Bv1-Au1)
   #  Pp2 = np.arctan2(Au2-Bu2,Au2+Bv2) # as per 5.8.47a
      Pp2 = np.arctan2(Av2-Bu2,Au2+Bv2) # corrected
      Pn2 = np.arctan2(Bu2+Av2,Bv2-Au2)

      # complex rotary spectra
      use_W=1
      if use_W:
        Wp1 = Ap1 * np.exp(-Pp1 * 1j)
        Wn1 = An1 * np.exp(-Pn1 * 1j)
        Wp2 = Ap2 * np.exp(-Pp2 * 1j)
        Wn2 = An2 * np.exp(-Pn2 * 1j)
    

      # inner cross spectrum, S
      if use_W:
        Sp12 += Wp1 * np.conj(Wp2) # this Sp12 may have wrong sign
        Sn12 += np.conj(Wn1) * Wn2
      else:
        Sp12 += Ap1 * Ap2 * np.exp(-(Pp1-Pp2) * 1j)  # omega > 0
        Sn12 += An1 * An2 * np.exp( (Pn1-Pn2) * 1j)  # omega < 0
    

      Sp11 += Ap1 * Ap1
      Sn11 += An1 * An1
      Sn22 += An2 * An2
      Sp22 += Ap2 * Ap2

      # for inner rotary coherence squared, gamma squared, gam2
      gam2_p12_1 = gam2_p12_1 + Ap1 * Ap2 * np.cos(Pp1-Pp2)
      gam2_p12_2 = gam2_p12_2 + Ap1 * Ap2 * np.sin(Pp1-Pp2)
      gam2_p12_3 = gam2_p12_3 + Ap1*Ap1
      gam2_p12_4 = gam2_p12_4 + Ap2*Ap2

      gam2_n12_1 = gam2_n12_1 + An1 * An2 * np.cos(Pn1-Pn2)
      gam2_n12_2 = gam2_n12_2 + An1 * An2 * np.sin(Pn1-Pn2)
      gam2_n12_3 = gam2_n12_3 + An1*An1
      gam2_n12_4 = gam2_n12_4 + An2*An2

      # for outer rotary cross spectrum, Y
      if use_W:
        Yp12 = Yp12 + Wn1 * Wp2
        Yn12 = Yn12 + Wp1 * Wn2
      else:
        Yp12 = Yp12 + An1 * Ap2 * np.exp((Pp2-Pn1) * 1j)
    #   Yn12 = Yn12 + Ap1 * An2 * np.exp((Pp1-Pp2) * 1j) # from book
        Yn12 = Yn12 + Ap1 * An2 * np.exp((Pp1-Pn2) * 1j) # corrected

      # outer rotary coherence squared, lambda squared, lam2
      lam2_p12_1 = lam2_p12_1 + An1 * Ap2
      lam2_p12_2 = lam2_p12_2 + np.cos(Pp2-Pn1)
    # lam2_p12_3 = lam2_p12_3 + np.sin(Pp2-Pp1) # as per 5.8.63
      lam2_p12_3 = lam2_p12_3 + np.sin(Pp2-Pn1) # corrected
      lam2_p12_4 = lam2_p12_4 + Ap2*Ap2
      lam2_p12_5 = lam2_p12_5 + An1*An1

      lam2_n12_1 = lam2_n12_1 + Ap1 * An2
      lam2_n12_2 = lam2_n12_2 + np.cos(Pp1-Pn2)
    # lam2_n12_3 = lam2_n12_3 + np.sin(Pp1-Pp2) # as per 5.8.63
      lam2_n12_3 = lam2_n12_3 + np.sin(Pp1-Pn2) # corrected
      lam2_n12_4 = lam2_n12_4 + Ap1*Ap1
      lam2_n12_5 = lam2_n12_5 + An2*An2

      # outer rotary coherence phase
      psi_p12_1 = psi_p12_1 + An1 * Ap2 * np.sin(Pn1-Pp2)
      psi_p12_2 = psi_p12_2 + An1 * Ap2 * np.cos(Pn1-Pp2)

      psi_n12_1 = psi_n12_1 + Ap1 * An2 * np.sin(Pn2-Pp1)
      psi_n12_2 = psi_n12_2 + Ap1 * An2 * np.cos(Pn2-Pp1)


  if npie_ok>0:
    Sp12 = Sp12 / npie_ok * 2 / fsamp / 2 * nfft
    Sn12 = Sn12 / npie_ok * 2 / fsamp / 2 * nfft

    Sp11 = Sp11 / npie_ok * 2 / fsamp / 2 * nfft
    Sn11 = Sn11 / npie_ok * 2 / fsamp / 2 * nfft
    Sp22 = Sp22 / npie_ok * 2 / fsamp / 2 * nfft
    Sn22 = Sn22 / npie_ok * 2 / fsamp / 2 * nfft

    gam2_p12_1 = gam2_p12_1 / npie_ok
    gam2_p12_2 = gam2_p12_2 / npie_ok
    gam2_p12_3 = gam2_p12_3 / npie_ok
    gam2_p12_4 = gam2_p12_4 / npie_ok

    gam2_n12_1 = gam2_n12_1 / npie_ok
    gam2_n12_2 = gam2_n12_2 / npie_ok
    gam2_n12_3 = gam2_n12_3 / npie_ok
    gam2_n12_4 = gam2_n12_4 / npie_ok

    lam2_p12_1 = lam2_p12_1 / npie_ok
    lam2_p12_2 = lam2_p12_2 / npie_ok
    lam2_p12_3 = lam2_p12_3 / npie_ok
    lam2_p12_4 = lam2_p12_4 / npie_ok
    lam2_p12_5 = lam2_p12_5 / npie_ok

    lam2_n12_1 = lam2_n12_1 / npie_ok
    lam2_n12_2 = lam2_n12_2 / npie_ok
    lam2_n12_3 = lam2_n12_3 / npie_ok
    lam2_n12_4 = lam2_n12_4 / npie_ok
    lam2_n12_5 = lam2_n12_5 / npie_ok

    psi_p12_1 = psi_p12_1 / npie_ok
    psi_p12_2 = psi_p12_2 / npie_ok

    psi_n12_1 = psi_n12_1 / npie_ok
    psi_n12_2 = psi_n12_2 / npie_ok

    # inner rotary coherence squared, gamma squared
    gam2_p12 = (np.power(gam2_p12_1,2) + np.power(gam2_p12_2,2)) / (gam2_p12_3 * gam2_p12_4)
    gam2_n12 = (np.power(gam2_n12_1,2) + np.power(gam2_n12_2,2)) / (gam2_n12_3 * gam2_n12_4)

    # inner rotary coherence phase lag, phi
    phi_p12 = np.arctan2(-np.imag(Sp12),np.real(Sp12))
    phi_n12 = np.arctan2(-np.imag(Sn12),np.real(Sn12))

    # should be identical
    phi_p12_chk = np.arctan2( gam2_p12_2,gam2_p12_1)
    phi_n12_chk = np.arctan2(-gam2_n12_2,gam2_n12_1)

    # check if identical
    if max(np.abs(phi_p12-phi_p12_chk) % 2.0*np.pi)>1e-9:
      print('chk=',max(np.abs(phi_p12-phi_p12_chk) % 2.0*np.pi))
      error("phi_p12 doesn't check")
  
    if max(np.abs(phi_n12-phi_n12_chk) % 2.0*np.pi)>1e-9:
      print('chk=',max(np.abs(phi_n12-phi_n12_chk) % 2.0*np.pi))
      error("phi_n12 doesn't check")
  

    # outer rotary coherence squared, lambda squared, lam2
    lam2_p12 = np.power(lam2_p12_1,2) * (np.power(lam2_p12_2,2) + np.power(lam2_p12_3,2)) / \
      (lam2_p12_4 * lam2_p12_5)
    lam2_n12 = np.power(lam2_n12_1,2) * (np.power(lam2_n12_2,2) + np.power(lam2_n12_3,2)) / \
      (lam2_n12_4 * lam2_n12_5)

    # outer rotary coherence phase
    psi_p12 = np.arctan2(psi_p12_1,psi_p12_2)
    psi_n12 = np.arctan2(psi_n12_1,psi_n12_2)

  R = namedtuple('RotCoh',[])
  R.F    = F
  R.Sp11 = Sp11
  R.Sn11 = Sn11
  R.Sp22 = Sp22
  R.Sn22 = Sn22
  R.Sp12 = Sp12
  R.Sn12 = Sn12
  R.gam2_p12 = gam2_p12
  R.gam2_n12 = gam2_n12
  R.phi_p12  = phi_p12
  R.phi_n12  = phi_n12
  R.lam2_p12 = lam2_p12
  R.lam2_n12 = lam2_n12
  R.psi_p12  = psi_p12
  R.psi_n12  = psi_n12
  R.npie  = npie_ok
  R.nbad  = npie_ng

  return R

def test_twospec():
  nfft  = 16
  fsamp = 1.0
  tsamp = 1.0 / fsamp
  t = np.arange(nfft*3) * tsamp
  freq  = fsamp / nfft * 3
  x = np.cos(np.pi * 2.0 * freq * t)
  y = np.sin(np.pi * 2.0 * freq * t)

  nfft = 1024
  nfft = 16

  fsamp = 10000.0
  tsamp = 1.0/fsamp
  T = nfft * tsamp

  freq = 0.0            # zero
  freq = (nfft/2) / T   # nyquist
  freq = 4.0/T 
  freq = 5.0/T 
  freq = 1.0/T            # lowest non-zero frequency

  npoints = 1 * nfft

  t = tsamp * np.arange(npoints) 

  omega = 2.0 * np.pi * freq
  x = 1.0 * np.cos(omega*t) + 2.0 * np.sin(omega*t)
  y = 3.0 * np.cos(omega*t) - 4.0 * np.sin(omega*t)

  S2 = twospec(x, y, nfft, fsamp)

# print('x=',x)
# print('y=',y)
# print('z=',S2.z)
  print('intgl_axx=',S2.intgl_axx)
  print('intgl_ayy=',S2.intgl_ayy)
  print('varx=',S2.varx)
  print('vary=',S2.vary)
  print('axx=',S2.axx)
  print('ayy=',S2.ayy)
  print('cxy=',S2.cxy)
  print('qxy=',S2.qxy)
  print('freq=',S2.freq)

def test_rotcoh():
  nfft  = 64
  fn    = 3.0
  fsamp = 1000
  nrep  = 100

  npie = nrep*2-1
  npts = nrep * nfft
  nslide = nfft/2

  npts = nfft * nrep
  tsamp = 1.0/fsamp
  T = nfft * tsamp
  df = 1.0/T
  freq = df * fn
  omega = 2.0 * np.pi * freq
  t = tsamp * np.arange(npts)

  c = np.cos(omega*t)
  s = np.sin(omega*t)
  u1 = 0.5*c + 0.0*s
  v1 = 0.0*c + 0.5*s # ccw

  # u2,v2 rotated from u1, v1
  ang = 60.0 * np.pi/180.0
  u2 = u1 * np.cos(ang) - v1 * np.sin(ang)
  v2 = v1 * np.cos(ang) + u1 * np.sin(ang)

  runno = 4
  if runno == 1:
    pass
  if runno == 2:
    v1 = -v1
    v2 = -v2
  if runno == 3:
    v2 = -v2
  if runno == 4:
    v1 = -v1

  noise = 1.0 / np.sqrt(12) / np.power(2.0,16) # ADC converter noise
  u1 += noise * np.random.randn(len(u1))
  v1 += noise * np.random.randn(len(v1))
  u2 += noise * np.random.randn(len(u2))
  v2 += noise * np.random.randn(len(v2))

  R = rotcoh(u1,v1,u2,v2,nfft,fsamp)

  print('R.gam2_p12=\n',R.gam2_p12)
  print('R.gam2_n12=\n',R.gam2_n12)
  print('R.lam2_p12=\n',R.lam2_p12)
  print('R.lam2_n12=\n',R.lam2_n12)


if __name__ == '__main__':
# test_autospec()
# test_twospec()
# test_rotcoh()
  test_nfft235()

# end of file
