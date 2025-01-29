def writepdf(pdfdir, namewoext):
  mymkdir(pdfdir)
  pdffile = os.path.join(pdfdir,namewoext+'.pdf')
  print('pdffile=',pdffile)
  plt.savefig(pdffile)
  os.system('/home/dunlap/bin/updateframe.run ' + pdfdir)

def writemat(matdir,namewoext,matdict):
  mymkdir(matdir)
  matfile = os.path.join(matdir,namewoext+'.mat')
  print('matfile=',matfile)
  scipy.io.savemat(matfile, matdict, format='4')

