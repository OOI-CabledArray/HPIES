#!/usr/bin/env python2.6
# digiterm.py -- terminal emulator to digi-tcp using select
# console port code is linux specific

from __future__ import print_function

import sys
import socket
import select
import Queue
import termios
import os
import tty
from optparse import OptionParser
import serial

from locktty import lock_tty, unlock_tty


def finish():
  try:
    termios.tcsetattr(coni, termios.TCSAFLUSH, saved)
  except:
    pass
  sys.exit('finished')

if __name__ == '__main__':

  ctrlv = False

  parser = OptionParser(usage="%prog [Options]",version="%prog 1.0")

  parser.add_option("-v", "--verbose",
  action="store_true", dest="verbose", default=False,
  help="print status messages to stdout")

  parser.add_option("-t", "--tty",
  action="store_true", dest="do_tty", default=False,
  help="use /dev/ttyR7")

  # for DigiConnect in office use: -j
  parser.add_option("-j", "--jbox", 
  action="store_true", dest="do_jbox", default=False,
  help="use 10.180.80.174:2101")

  (options, args) = parser.parse_args()

  if options.do_jbox and options.do_tty:
    print('only -j or -t can be used')
    parser.print_help()
    finish()

  if not options.do_jbox and not options.do_tty:
    print('either -j or -t must be used')
    parser.print_help()
    finish()

  if options.do_jbox:
    host = '10.180.80.174'
    port = 2101
    try:
      tcp = socket.socket (socket.AF_INET, socket.SOCK_STREAM)
      tcp.connect((host,port))
      tcp.setblocking(0)
    except:
      print('cannot open host=',host,'port=',port)
      finish()

  if options.do_tty:
    dev = '/dev/ttyRP7'
    baud = 38400
    lock_tty(dev)
    try:
      ser = serial.Serial(dev, baud , timeout=0)
    except:
      print('cannot open dev=',dev,'baud=',baud)
      finish()

  # setup console
  coni = sys.stdin.fileno()
  cono = sys.stdout.fileno()
  saved = termios.tcgetattr(coni)
  new = termios.tcgetattr(coni)
  new[3] = new[3] & ~termios.ICANON & ~termios.ECHO & ~termios.ISIG
  new[6][termios.VMIN] = 1
  new[6][termios.VTIME] = 0
  termios.tcsetattr(coni, termios.TCSANOW, new)
  tty.setraw(coni)

  # setup for select()
  message_queues = {}
  message_queues[cono] = Queue.Queue()
  outputs = []
  if options.do_tty:
    inputs  = [ser, coni]
    message_queues[ser] = Queue.Queue()
  if options.do_jbox:
    inputs  = [tcp, coni]
    message_queues[tcp] = Queue.Queue()

  count_select = 0

  while inputs:
    count_select += 1

    readable, writable, exceptional = select.select(inputs, outputs, inputs)

    for rdbl in readable:
      if options.do_tty and rdbl is ser:
        serdat = ser.read(9999)
        if len(serdat):
          message_queues[cono].put(serdat)
          if cono not in outputs:
            outputs.append(cono)

      if options.do_jbox and rdbl is tcp:
        tcpdat = tcp.recv(1024)
        if tcpdat:
          message_queues[cono].put(tcpdat)
          if cono not in outputs:
            outputs.append(cono)
        else:
          print('s.recv indicates closed.  exiting.')
          finish()

      if rdbl is coni:
        condat = os.read(coni,1)
        if condat:
          oc = ord(condat)
          if oc == 3:
            print('console got Ctrl-C.  exiting.\r\n')
            print('count_select=',count_select)
            finish()

          if ctrlv == False:
            if oc != 22: # Ctrl-V
              constr = condat
            else:
              ctrlv = True
              contmp = ''
              constr = ''
          else:
            if oc != 13: # CR
              contmp = contmp + condat
            else:
              ctrlv = False
              num = int(contmp)
              if num >= 0 and num <= 255:
                print('\r\nsending num=',num,'decimal\r\n')
                constr = chr(num)
              else:
                print('\r\num=',num,'decimal: not in range [0,255]\r\n')
                constr = ''

          if len(constr):
            if options.do_tty:
              message_queues[ser].put(constr)
              if ser not in outputs:
                outputs.append(ser)
            if options.do_jbox:
              message_queues[tcp].put(constr)
              if tcp not in outputs:
                outputs.append(tcp)

    for wrtbl in writable:
      try:
        next_msg = message_queues[wrtbl].get_nowait()
      except Queue.Empty:
        outputs.remove(wrtbl)
      else:
        if options.do_jbox and wrtbl is tcp:
          tcp.send(next_msg)
        if options.do_tty and wrtbl is ser:
          ser.write(next_msg)
        if wrtbl is cono:
          for i in range(len(next_msg)):
            c = next_msg[i]
            co = ord(c)
            if co > 32 and co < 127:
              cs = c
            else:
              cs = '\r\nreceived num= {0} decimal\r\n'.format(co)
            os.write(cono,cs)

    for xcptnl in exceptional:
      if xcptnl is ser:
        print('exceptional = ser')
      if xcptnl is tcp:
        print('exceptional = tcp')
      if xcptnl is cono:
        print('exceptional = cono')
      if xcptnl is coni:
        print('exceptional = coni')
      finish()

  finish()
