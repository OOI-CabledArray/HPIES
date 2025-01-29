#!/bin/usr/env python2

from __future__ import print_function

import numpy as np
import spectrajhd

print('xyz')

nfft  = 150
fsamp = 1.0
tsamp = 1.0 / fsamp
t = np.arange(nfft*3) * tsamp
freq  = fsamp / 8.0
freq  = fsamp / nfft * 3
x = np.sin(np.pi * 2.0 * freq * t)
S = spectrajhd.autospec(x, nfft, fsamp, 'none')

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
