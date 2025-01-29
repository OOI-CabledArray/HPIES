#!/usr/bin/env python2.6
# digiterm.py -- terminal emulator with three ports using select
#                console, digi-serial and digi-tcp
#                console port code is linux specific

from __future__ import print_function

import sys
import socket
import select
import Queue
import serial
import termios
import os
import tty

# main

ttyname  = '/dev/ttyRP4'
baud = 9600

host = '10.180.80.174'
port = 2101

tcp = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
tcp.connect((host,port))
tcp.setblocking(0)

ser = serial.Serial(ttyname, baud , timeout=0)

coni = sys.stdin.fileno()
cono = sys.stdout.fileno()

old = termios.tcgetattr(coni)

new = termios.tcgetattr(coni)
new[3] = new[3] & ~termios.ICANON & ~termios.ECHO & ~termios.ISIG
new[6][termios.VMIN] = 1
new[6][termios.VTIME] = 0
termios.tcsetattr(coni, termios.TCSANOW, new)
tty.setraw(coni)

def getkey():
  return os.read(coni,1)

def finish():
  termios.tcsetattr(coni, termios.TCSAFLUSH, old)
  sys.exit(0)

inputs = [ tcp, ser, coni ]
outputs = [ ]

# outgoing message queues
message_queues = {}
message_queues[tcp] = Queue.Queue()
message_queues[ser] = Queue.Queue()
message_queues[cono] = Queue.Queue()


while inputs:
  readable, writable, exceptional = select.select(inputs, outputs, inputs)

  for s in readable:
    if s is tcp:
      tcpdat = tcp.recv(1024)
      if tcpdat:
#       print('tcpdat=',tcpdat)
        message_queues[ser].put(tcpdat)
        if ser not in outputs:
          outputs.append(ser)
        message_queues[cono].put(tcpdat)
        if cono not in outputs:
          outputs.append(cono)
      else:
        print('s.recv indicates closed.  exiting.')
        finish()
    if s is ser:
      serdat = ser.read(9999)
      if serdat:
#       print('serdat=',serdat)
        message_queues[tcp].put(serdat)
        if tcp not in outputs:
          outputs.append(tcp)
    if s is coni:
      condat = os.read(coni,1)
      if condat:
        if condat == chr(3):
          print('console got Ctrl-C.  exiting.\r\n')
          finish()
#       print('ord(condat)=',ord(condat))
        message_queues[tcp].put(condat)
        if tcp not in outputs:
          outputs.append(tcp)

  for s in writable:
    try:
      next_msg = message_queues[s].get_nowait()
    except Queue.Empty:
      outputs.remove(s)
    else:
      if s is tcp:
        tcp.send(next_msg)
      if s is ser:
        ser.write(next_msg)
      if s is cono:
        os.write(cono,next_msg)

  for s in exceptional:
    if s is tcp:
      print('exceptional = tcp')
    if s is ser:
      print('exceptional = ser')
    if s is cono:
      print('exceptional = cono')
    if s is coni:
      print('exceptional = coni')

finish()
