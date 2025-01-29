#! /usr/bin/env python2

from __future__ import print_function
import collections
import numpy

fifo = collections.deque(maxlen=3)

a = numpy.array([11,12,13,14,15])
b = numpy.array([21,22,23,24,25,26])
c = numpy.array([31,32,33,34,35,36,37])
d = numpy.array([41,42,43,44,45,46,47,48])

E = collections.namedtuple('Enam',[])
E.e = [0,1,2]
E.f = [3,4,5,6]


fifo.append(a); print('len(fifo)=',len(fifo))
fifo.append(b); print('len(fifo)=',len(fifo))
fifo.append(c); print('len(fifo)=',len(fifo))
fifo.append(d); print('len(fifo)=',len(fifo))
fifo.append(E); print('len(fifo)=',len(fifo))

print('len(fifo)',len(fifo))
print('fifo[0]=',fifo[0])
print('fifo[1]=',fifo[1])
print('fifo[2].e=',fifo[2].e)
print('fifo[2].f=',fifo[2].f)

