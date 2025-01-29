from __future__ import print_function

import sys
from calendar import timegm
from time import strptime

def get_hpies_info(leafname):
  import collections
  INFO = collections.namedtuple('HPIES_INFO',[])

  tok = leafname.split('_')
# tok= ['HPIESA301', '10.31.5.5', '2101', '20150601T0000', 'UTC.dat']
  uxt_fn = timegm(strptime(tok[3],'%Y%m%dT%H%M'))


  INFO.flip_e1a = False
  INFO.flip_e1b = False
  INFO.flip_e1c = False
  INFO.flip_e2a = False
  INFO.flip_e2b = False
  INFO.flip_e2c = False

  INFO.ws1 = None
  INFO.ws2 = None
  INFO.e1bsf = None

  if leafname.find('H1') == 0:
    INFO.deployment = 'OKMC'
    INFO.runid = 'H1'
    INFO.uxt_poly_ref = 1339106500.0
    INFO.uxt_poly_coef = [  1.58670248e-14,  3.69147464e-07,  9.01547311e-01]
    INFO.esep1 = 2.25 # m
    INFO.esep2 = 2.25 # m
    INFO.magvar = -2.6
    INFO.heading = 74.2
    INFO.compass_angle = 180.0 # e2 points south when heading = 0
    INFO.ws1 = 'ws25'
    INFO.ws2 = 'ws26'
    INFO.flip_e1a = True
    INFO.flip_e1b = True
    INFO.flip_e1c = True
    INFO.flip_e2a = True
    INFO.flip_e2b = True
    INFO.flip_e2c = True
  elif leafname.find('H2') == 0:
    INFO.deployment = 'OKMC'
    INFO.runid = 'H2'
    INFO.uxt_poly_ref = 1339128702.0
    INFO.uxt_poly_coef = [ -3.19513404e-15,  1.54736995e-05,  1.29811396e+00]
    INFO.esep1 = -2.25 # m -- flipped sign for new WS
    INFO.esep2 = -2.25 # m -- flipped sign for new WS
    INFO.magvar = -2.6
    INFO.heading = 86.5
    INFO.compass_angle = 180.0 # e2 points south when heading = 0
    INFO.ws1 = 'ws11'
    INFO.ws2 = 'ws13'
#   INFO.flip_e1a = True
#   INFO.flip_e1b = True
#   INFO.flip_e1c = True
#   INFO.flip_e2a = True
#   INFO.flip_e2b = True
    INFO.flip_e2c = True
  elif leafname.find('H3') == 0:
    INFO.deployment = 'OKMC'
    INFO.runid = 'H3'
    INFO.uxt_poly_ref = 1339190501.0
    INFO.uxt_poly_coef = [  6.54613352e-15,  9.76169175e-07,  7.30209617e-01]
    INFO.esep1 = 2.25 # m
    INFO.esep2 = 2.25 # m
    INFO.magvar = -2.6
    INFO.heading = 182.8
    INFO.compass_angle = 180.0 # e2 points south when heading = 0
    INFO.ws1 = 'ws23'
    INFO.ws2 = 'ws27'
    INFO.flip_e1a = True
    INFO.flip_e1b = True
    INFO.flip_e1c = True
  elif leafname.find('H4') == 0:
    INFO.deployment = 'OKMC'
    INFO.runid = 'H4'
    INFO.uxt_poly_ref = 1339217500.0
    INFO.uxt_poly_coef = [ -2.38816069e-15, -3.17793271e-05,  9.89743586e-02]
    INFO.esep1 = -2.25 # m -- flipped sign for new WS
    INFO.esep2 =  2.25 # m
    INFO.magvar = -2.6
    INFO.heading = 58.4
    INFO.compass_angle = 180.0 # e2 points south when heading = 0
    INFO.ws1 = 'ws10'
    INFO.ws2 = 'ws28'
#   INFO.flip_e1a = True
#   INFO.flip_e1b = True
    INFO.flip_e1c = True
  elif leafname.find('H5') == 0:
    INFO.deployment = 'OKMC'
    INFO.runid = 'H5'
    INFO.uxt_poly_ref = 1339380101.0
    INFO.uxt_poly_coef = [  3.28728101e-15, -1.75745563e-05,  7.82347793e-01]
    INFO.esep1 =  2.25 # m
    INFO.esep2 = -2.25 # m -- flipped sign for new WS
    INFO.magvar = -2.6
    INFO.heading = 49.7 
    INFO.compass_angle = 180.0 # e2 points south when heading = 0
    INFO.ws1 = 'ws24'
    INFO.ws2 = 'ws12'
    INFO.flip_e1a = True
    INFO.flip_e1b = True
    INFO.flip_e1c = True
  elif leafname.find('HPIESA301') == 0:
    INFO.deployment = 'RSN'
    INFO.uxt_poly_ref = 0.0
    INFO.uxt_poly_coef = [0.0, 0.0, 0.0]
    INFO.esep1 = 2.25
    INFO.esep2 = 2.25
    INFO.magvar = 0.0
    INFO.heading = 0.0
    INFO.compass_angle = 180.0
    INFO.flip_e1a = True
    INFO.flip_e1b = True
    INFO.flip_e1c = False
    INFO.flip_e2a = True
    INFO.flip_e2b = True
    INFO.flip_e2c = False
    if uxt_fn < timegm(strptime('20150725T0000','%Y%m%dT%H%M')):
      INFO.runid = 'AB1'
      INFO.ws1 = 'ws32'
      INFO.ws2 = 'ws34'
    else:
      INFO.runid = 'AB4'
      INFO.ws1 = ''
      INFO.ws2 = ''

  elif leafname.find('HPIESA101') == 0:
    INFO.deployment = 'RSN'
    INFO.uxt_poly_ref = 0.0
    INFO.uxt_poly_coef = [0.0, 0.0, 0.0]
    INFO.esep1 = 2.25
    INFO.esep2 = 2.25
    INFO.magvar = 0.0
    INFO.heading = 0.0
    INFO.compass_angle = 180.0
    if uxt_fn < timegm(strptime('20150722T0000','%Y%m%dT%H%M')):
      INFO.runid = 'SB2'
      INFO.e1bsf = 6.0
      INFO.ws1 = 'ws31'
      INFO.ws2 = 'ws33'
    else:
      INFO.runid = 'SB3'
      INFO.e1bsf = None
      INFO.ws1 = ''
      INFO.ws2 = ''
      INFO.flip_e1c = True
      INFO.flip_e2c = True

  elif leafname.find('hef009a') == 0:
    INFO.deployment = 'LAB'
    INFO.runid = 'hef009a'
    INFO.uxt_poly_ref = 0.0
    INFO.uxt_poly_coef = [0.0, 0.0, 0.0]
    INFO.esep1 = 2.25
    INFO.esep2 = 2.25
    INFO.magvar = 0.0
    INFO.heading = 0.0
    INFO.compass_angle = 180.0
    INFO.ws1 = 'ws30'
    INFO.ws2 = 'ws40'
  elif leafname.find('hef009b') == 0:
    INFO.deployment = 'LAB'
    INFO.runid = 'hef009b'
    INFO.uxt_poly_ref = 0.0
    INFO.uxt_poly_coef = [0.0, 0.0, 0.0]
    INFO.esep1 = 2.25
    INFO.esep2 = 2.25
    INFO.magvar = 0.0
    INFO.heading = 0.0
    INFO.compass_angle = 180.0
    INFO.ws1 = 'ws41'
    INFO.ws2 = 'ws35'
  else:
    print('unknown runid=',runid,'leafname=',leafname)
    sys.exit(1)
    INFO.deployment = None
    INFO.runid = leafname
    INFO.uxt_poly_ref = 0.0
    INFO.uxt_poly_coef = [0.0, 0.0, 0.0]
    INFO.esep1 = 2.25
    INFO.esep2 = 2.25
    INFO.magvar = 0.0
    INFO.heading = 0.0
    INFO.compass_angle = 180.0

  INFO.ymd = None
  if INFO.deployment == 'OKMC' or INFO.deployment == 'LAB':
    # assume filenames from hefcf2.py
    x = leafname.split('_') 
    if len(x) == 2:
      INFO.ymd = x[1]
  if INFO.deployment == 'RSN':
    # assume filenames from APL engineering data
    x = leafname.split('_')
    if len(x) == 5:
      y = x[3].split('T')
      if len(y)==2 and y[1]=='0000':
        INFO.ymd = y[0]
      else:
        INFO.ymd = x[3]
  if INFO.ymd == None:
    print('get_hpies_info.py: warning: no ymd time found in leafname=',leafname)
      

  return INFO
