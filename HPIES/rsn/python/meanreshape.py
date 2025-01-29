#!/usr/bin/env python2
# meanreshap.py -- test mean(reshape())

from __future__ import print_function

import os
import sys
import numpy as np

# main
if __name__ == '__main__':

  a = [1,2,3,4,5,6,7,8,9,10,11,12,13]
  m = 3
  n = int(len(a) / m)
  af = np.mean(np.reshape(a[0:n*m],(m,n),order='F'),0)
  print('n=',n,'len(af)=',len(af))
  print('af=',af)
