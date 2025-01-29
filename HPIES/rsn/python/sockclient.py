#!/usr/bin/env python2.6
# sockclient.py -- socket client

from __future__ import print_function

import os
import sys
import socket
import subprocess
from collections import deque

# main

port = 29346
host = socket.gethostname()
print('sockserv.py: host=',host,'port=',port)

s = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
s.connect((host,port))
sys.stdout.write('%')

while True:
  line = sys.stdin.readline()
  if line == ' ':
    break
  s.send(line)
  ibuf = s.recv(1024)
  sys.stdout.write(ibuf)
  sys.stdout.write('%')
s.close()

  
s.send("hello")
print('received:',ibuf)
s.close
