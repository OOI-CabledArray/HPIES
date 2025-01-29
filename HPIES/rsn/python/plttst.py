#!/usr/bin/env python2
# plttst.py -- fiddle with plot_date limits

from __future__ import print_function

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import numpy as np

def fix_xdates(axlist,nss):
  fmt = mdates.DateFormatter('%b %d %H%MZ')
  for ax in axlist:
    xax = ax.get_xaxis()
    xax.set_major_formatter( fmt )

    dx = 0.03
    xlim = ax.get_xlim()
    if xlim[-1] - xlim[0] > 0.5:
      ax.set_xlim([np.floor(xlim[0]+dx)-dx,np.ceil(xlim[-1]-dx)+dx])

    if nss>1:
      xticks = ax.get_xticks()
      newxticks = []
      for i in range(0,len(xticks),nss):
        newxticks.append(xticks[i])
      ax.set_xticks(newxticks)

pi = np.pi
twopi = 2.0 * pi


fig = plt.figure(num=1,figsize=(10, 7))
fig.clf()

mlt = np.array([0.1,0.4,0.7,0.9]) + 100
x = np.cos(mlt * twopi)
y = np.sin(mlt * twopi)

print('mlt=',mlt)

axlist=[]
ax = fig.add_subplot(2,1,1)
axlist.append(ax)
ax.plot_date(mlt,x,'r.-')

ax = fig.add_subplot(2,1,2)
axlist.append(ax)
ax.plot_date(mlt,y,'b.-')

fig.autofmt_xdate()
fix_xdates(axlist,1)

plt.savefig('plttst.pdf')

