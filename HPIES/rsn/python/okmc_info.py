def okmc_info(leafname):

  if leafname.find('H1') >= 0:
    runid = 'H1'
    uxt_hef0 = 1339106500.0
    uxt_coef = [  1.58670248e-14,  3.69147464e-07,  9.01547311e-01]
  elif leafname.find('H2') >= 0:
    runid = 'H2'
    uxt_hef0 = 1339128702.0
    uxt_coef = [ -3.19513404e-15,  1.54736995e-05,  1.29811396e+00]
  elif leafname.find('H3') >= 0:
    runid = 'H3'
    uxt_hef0 = 1339190501.0
    uxt_coef = [  6.54613352e-15,  9.76169175e-07,  7.30209617e-01]
  elif leafname.find('H4') >= 0:
    runid = 'H4'
    uxt_hef0 = 1339217500.0
    uxt_coef = [ -2.38816069e-15, -3.17793271e-05,  9.89743586e-02]
  elif leafname.find('H5') >= 0:
    runid = 'H5'
    uxt_hef0 = 1339380101.0
    uxt_coef = [  3.28728101e-15, -1.75745563e-05,  7.82347793e-01]
  else:
    runid = None
    uxt_hef0 = None
    uxt_coef = None

  return runid, uxt_hef0, uxt_coef
