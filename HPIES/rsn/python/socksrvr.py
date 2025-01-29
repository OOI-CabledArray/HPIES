#!/usr/bin/env python2.6
# socksrvr.py -- socket server

from __future__ import print_function

import os
import sys
import socket
import subprocess
from collections import deque

# main

s = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
host = socket.gethostname()
print('sockserv.py: host=',host)
port = 29346
s.bind((host,port))
s.listen(5)

while True:
  c, addr = s.accept()
  data = c.recv(1024)
  if data:
    print('data=',data)
    c.send(data)
  c.close()
