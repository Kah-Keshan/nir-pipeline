from ds9disp import dodisp, panto, wcs_mark
import os
from datetime import datetime as d
wd=os.getcwd()
os.chdir('/Volumes/BigDalek/khaled/')
from pyraf import iraf
os.chdir(wd)
iraf.stsdas()
iraf.analysis()
iraf.isophote()
iraf.fitting()
import ds9 as pyds9
import pyfits
import pyfits as pf
import subprocess as sub
import asciidata as ad
import matplotlib
matplotlib.use("Agg")
import pylab as pl
import time
import pywcs
import aplpy as ap
import scipy.special as spec
import scipy.optimize as opt
from aplpy import make_rgb_image as RGB
matplotlib.rc('verbose',level='silent')
import matplotlib.cm as cmap
import astLib.astCoords as ac
from iraf import daophot
iraf.unlearn(iraf.isophote)
import sys


DISP = True
DISP = False
if DISP:
  ds9 = pyds9.ds9()
  #  ds9 = pysao.ds9()
  ds9.set('scale zscale')
  ds9.set('tile yes')
  ds9.set('tile grid')

HICAT=ad.open('/Volumes/BigDalek/khaled/nir_pipeline_NHK/lists/hicatlist')
hidat = ad.open('/Volumes/BigDalek/khaled/nir_pipeline_NHK/lists/masterHIv',delimiter='|')


Csystem = "galactic"
# SEx paramters
SEconfig = './SEconfig'
SEparam = './SE.param'
SEnnw = './default.nnw'

#other variables
#tmppath= "/dev/shm/"
tmppath= './tmp/'

# paths
imagepath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/'
catpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/hicats/'
psfpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/psf/'
stamppath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/stamps/'
ellpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/reduction/phot/SBPs/'
outcatpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/reduction/phot/cats/'

# global image parameters
#expos = 24.	# exposure time
ps = 0.45	# pixel scale
sig_mzp = 0.02	# error on magzp determinatiion ### NB this is a guess for now ###
area = ps*ps


bands=['j','h','k']
pfmts=['b,','g,','r,']
lfmts=['b','g','r']
cols=['blue','green','red']


# routine for finding index of search term in AsciiColumn
def afind(find, incolumn):
    index=[]
    for i in range(len(incolumn)):
	if find==incolumn[i]:
	    index=index+[i]
    return index

###################################################################
def set_ellipse_def():

  iraf.unlearn(iraf.isophote)
  ## IRAF ELLIPSE DEFAULTS ##
  iraf.ellipse.interactive='no'
  iraf.ellipse.verbose='no'

  #iraf.geompar.minsma = 0.0
  iraf.geompar.maxsma = 'INDEF'
  iraf.geompar.step=0.02
  #iraf.geompar.linear = 'no'
  #iraf.geompar.maxrit = 'INDEF'
  iraf.geompar.recenter = 'no'
  iraf.geompar.xylearn = 'no'
  #iraf.geompar.physical = 'yes'

  iraf.controlpar.conver = 0.2
  iraf.controlpar.minit = 10
  iraf.controlpar.maxit = 100
  iraf.controlpar.hcenter = 'yes' 
  iraf.controlpar.hellip = 'no'
  iraf.controlpar.hpa = 'no'
  #iraf.controlpar.wander = 2.
  iraf.controlpar.maxgerr = 0.75
  iraf.controlpar.olthresh = 0.
  iraf.controlpar.soft = 'yes'

  #iraf.samplepar.integrmode = 'bi-linear'
  #iraf.samplepar.usclip = 3.0
  #iraf.samplepar.lsclip = 3.0 
  iraf.samplepar.nclip = 10
  #iraf.samplepar.fflag = 0.5
  #iraf.samplepar.sdevice = 'none'
  #iraf.samplepar.tsample = 'none'
  #iraf.samplepar.absangle = 'yes'
  #iraf.samplepar.harmonics = 'none'

  # mag pars: m = mag0 - 2.5log[(intens-xerolevle)/refer]
  iraf.magpar.mag0 = 0.0
  #iraf.magpar.refer = 1.0
  #iraf.magpar.zerolevel = 0.0
  return
def smooth(x,window_len=11,window='hanning'):
    import numpy
    """smooth the data using a window with requested size.
    
    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal 
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.
    
    input:
        x: the input signal 
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal
        
    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)
    
    see also: 
    
    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter
 
    TODO: the window parameter could be the window itself if an array instead of a string   
    """

    if x.ndim != 1:
        raise ValueError, "smooth only accepts 1 dimension arrays."

    if x.size < window_len:
        raise ValueError, "Input vector needs to be bigger than window size."


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError, "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'"


    s=numpy.r_[2*x[0]-x[window_len:1:-1],x,2*x[-1]-x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=ones(window_len,'d')
    else:
        w=eval('numpy.'+window+'(window_len)')

    y=numpy.convolve(w/w.sum(),s,mode='same')
    return y[window_len-1:-window_len+1]

#############################
def findgalcoor(im,x,y,r,dr=10.):
  '''
# find intensity weighted mean position (X,Y) in image in a box 2r around x,y
  '''
  data = pyfits.getdata(im)
  head = pyfits.getheader(im)
  ny=head.get('NAXIS1')
  nx=head.get('NAXIS2')
  print 'findgalcoor size %.1f'%(r)
  #print nx,ny
  print data.shape
  print x,y,r
  r = int(r)

  box = pl.array(data[max(int(x-r),1):min(int(x+r)+1,nx),max(int(y-r),1):min(int(y+r)+1,ny)])

  xind = 1.0*pl.arange(max(int(x-r),1),min(int(x+r)+1,nx))
  yind = 1.0*pl.arange(max(int(y-r),1),min(int(y+r)+1,ny))
  Y = pl.sum((box)*yind)/pl.sum(box)+0.5
  X = pl.sum(pl.transpose(box)*xind)/pl.sum(box)+0.5

  if pl.sqrt((x-X)**2. + (y-Y)**2.) > dr:
    print 'error: new coords >%.1f pix from initial guess' %(dr)
    print 'x0,y0 = %.1f,%.1f / X,Y = %.1f,%.1f' %(x,y,X,Y)
    X = x
    Y = y
  Peak = box.max()

  return X,Y

def ellopen(ellfile,null=['INDEF'],comments='#',delimiter=''):
  '''
  ellopen: open ellipse file and return array of data, nulls=NaN
  '''
  f = open(ellfile,'r')
  X = []
  for line in f:
    if line[0]!=comments:
      if len(delimiter) < 1: row = [val for val in line.split()]
      else: row = [val for val in line.split(delimiter)]
      r=[]
      for i in range(len(row)):
	try: r.append(float(row[i]))
	except ValueError:
	  if (row[i] in null): r.append(pl.NaN)
	  else: r.append(-9999)
      X.append(r)
  X = pl.array(X)
  f.close()
  return X.transpose()

####################################################

def fit_ell(name,fignum,band,image,s,ss,x,y,ell0,pa0,startrad='INDEF',minrad=10.0,maxrad='INDEF',col='k'):
  '''
  fit_ell(image,s,ss,x,y,ell0,pa0,startrad='INDEF',minrad=10.0,maxrad='INDEF')
  run IRAF ellipse on image (sky = s, sigma = ss)
  - x,y = central coords (fixed)
  - ell0,pa0 = inital guess for ellipse
  - startrad = start at sma startrad
  - minrad = go inwards to minrad
  - maxrad = go out to maxrad ('INDEF' stop when gradient/SNR conditions met)
  '''
  outtbl = tmppath+'tmp_ellip'
  outtbldat = tmppath+'tmp_ellip.dat'
  if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
  if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))
#run ellipse fitting
  set_ellipse_def()  # set ellipse defaults
  if ell0 < 0.05: ell0 = 0.05
  if startrad < 5.0: startrad = 5.0
  minrad = startrad
  ellsuccess = True
  print"test -----------------------------------"
  try:
    print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(image,outtbl,x,y,ell0,pa0,startrad,minrad,maxrad)
    iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell0,pa0=pa0,sma0=startrad,minsma=minrad,maxsma=maxrad)
  except:
    iraf.samplepar.nclip = 1
    print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(image,outtbl,x,y,ell0,pa0,startrad,minrad,maxrad)
    iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell0,pa0=pa0,sma0=startrad,minsma=minrad,maxsma=maxrad)
  if ellsuccess:#try:
    iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,ELLIP,ELLIP_ERR,PA,PA_ERR,X0,X0_ERR,Y0,Y0_ERR')
    # plot output
    ell = ellopen(outtbldat)
    rad = ell[0]*ps
    pa=ell[5]
    ellip=ell[3]

    ell[1] -= s
    I1 = pl.sum(ell[1]>1.5*ss)
    I2 = pl.sum(ell[1]>2.5*ss)
    if I1 == len(ell[1]): I1 -= 1
    if I2 == len(ell[1]): I2 -= 1
    print ' I2, I1 = %.1f, %.1f' %(I2,I1)
    print ' sky = %.3f %.3f' %(s,ss)
    print ' cnts = %.2f %.2f' %(ell[1][I1],ell[1][I2])
    
    R1s = rad[I1]
    R2s = rad[I2]
    print ' R2, R1 = %.1f, %.1f arcsec' %(R2s,R1s)
    
    # average pa and ell
    pa1s = pa[I1]+90.
    pa2s = pa[I2]+90.
    PA = pl.average(pa[I2:I1+1])
    PAs = pl.sqrt(pl.sum(ell[6][I2:I1+1]**2.))/len(pa[I2:I1+1]) 
    #normalise PA to (-90,90]
    if PA<-90:
      th = PA+90
      PA = 90 - th
    if PA > 90:
      th = PA - 90
      PA = -90 + th
	
    el1s = ell[3][I1]
    el2s = ell[3][I2]
    ELL = pl.average(ell[3][I2:I1+1])
    ELLs = pl.sqrt(pl.sum(ell[4][I2:I1+1]**2.))/len(ell[3][I2:I1+1])#pl.std(ell[3][I2:I1+1])

    if DISP:
      if col=='b': color='blue'
      elif col=='g': color='green'
      elif col=='r': color='red'
      wcs_mark(ds9,x,y,R2s,(el2s-1)*R2s,pa2s,colour=color)
      wcs_mark(ds9,x,y,R1s,(el1s-1)*R1s,pa1s,colour=color)

    # if the ellipse is drastically different take the original 
    # for pos angle ignore circular galaxies
    '''
    if ELL < 0.8:
      print 'checking pa'
      print PA,pa0,PA-pa0
      if abs(PA-pa0) > 20:
	if PA < 0 : tPA = PA + 180.
	else: tPA = PA
	if pa0 < 0: tpa0 = pa0 + 180.
	else: tpa0 = pa0
	print tPA,tpa0
	if abs(tPA-tpa0) > 20:
	  print 'pos. angle differs too greatly from initial input'
	  print ' pa_in = %.1f, pa_out = %.1f, choosing %.1f +/- 10' %(pa0,PA,pa0)
	  PA = pa0
	  PAs = 10.
    if abs(ELL-ell0) > 0.2:
      print 'ellipticity differs too greatly from initial input'
      print ' ell_in = %.3f, ell_out = %.3f, choosing %.3f +/- 0.2' %(ell0,ELL,ell0)
      ELL = ell0
      ELLs = 0.2   # guess the guestimate error to be 0.2
    '''

    print ' PA(avg) = %.2f +/- %.2f deg ccw from y' %(PA,PAs)
    print ' ELL(avg) = %.2f +/- %.2f' %(ELL,ELLs)

    # add to plot
    fig = pl.figure(fignum)
    if band == 'j': n=0
    if band == 'h': n=1
    if band == 'k': n=2
    ax1 = pl.subplot(3,3,1+n)
    ax2 = pl.subplot(3,3,4+n,sharex=ax1)
    ax3 = pl.subplot(3,3,7+n,sharex=ax1)
    if band == 'j':
      ax1.set_ylabel('$\phi$ [deg]')
      ax2.set_ylabel(r'$\epsilon$')
      ax3.set_ylabel(r'$I(a)$')
    ax3.set_xlabel('$a$ [arcsec]')

    ax1.errorbar(rad,pa,ell[6],fmt=col+',',capsize=0)
    ax2.errorbar(rad,ell[3],ell[4],fmt=col+',',capsize=0)
    ax3.errorbar(rad,ell[1],ell[2],fmt=col+',',capsize=0)
    x1l,x1u = ax1.get_xlim()
    x2l,x2u = ax2.get_xlim()
    x3l,x3u = ax3.get_xlim()
    y3l,y3u = ax3.get_ylim()
    # vertiacl R2, R3 lines on all plots
    #y1l,y1u = ax1.get_ylim()
    ax1.set_ylim(-90,90)
    ax1.hlines(PA,x1l+0.1,x1u-0.1,'k','solid')
    ax1.vlines(R1s,-89,89,col,'solid')
    ax1.vlines(R2s,-89,89,col,'dashed')
    #y2l,y2u = ax2.get_ylim()
    ax2.hlines(ELL,x1l+0.1,x1u-0.1,'k','solid')
    ax2.vlines(R1s,0,1,col,'dashed')
    ax2.vlines(R2s,0,1,col,'solid')

    ax3.hlines(0,x3l+0.1,x3u-0.1,'k','solid')
    ax3.hlines(ss,x3l+0.1,x3u-0.1,'k','dashed')
    ax3.vlines(R1s,0.,0.99*y3u,col,'dashed')
    ax3.vlines(R2s,0.,0.99*y3u,col,'solid')

#print '########### SS = ',ss

    #pl.text(,,'$\phi=%.2f\pm%.2f$'%(),ax1.transAxes)

    ax1.set_ylim(-90,90)
    ax2.set_ylim(0,1)
    #if band == 'k':
    #  pl.subplots_adjust(left=0.05, right = 0.95)
    ellplotname = '%s%s_ell'%(ellpath,name)
    pl.savefig(ellplotname+'.png',dpi=200)
    pl.savefig(ellplotname+'.eps',dpi=200)
    if DISP and col=='k':
      os.system('xv %s &'%(ellplotname))
      pl.clf()
  
  else:#except:
    print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(image,outtbl,x,y,ell0,pa0,startrad,minrad,maxrad)
    
    print 'ERROR in FIT_ELL'

    PA = pa0
    PAs = 10.
    ELL = ell0
    ELLs = 0.2 
    ell = pl.zeros((10,11))

    fig = pl.figure(fignum)
    if band == 'j': n=0
    if band == 'h': n=1
    if band == 'k': n=2
    ax1 = pl.subplot(3,3,1+n)
    ax2 = pl.subplot(3,3,4+n,sharex=ax1)
    ax3 = pl.subplot(3,3,7+n,sharex=ax1)
    if band == 'j':
      ax1.set_ylabel('$\phi$ [deg]')
      ax2.set_ylabel(r'$\epsilon$')
      ax3.set_ylabel(r'$I(a)$')
    ax3.set_xlabel('$a$ [arcsec]')

    x1l,x1u = ax1.get_xlim()
    x2l,x2u = ax2.get_xlim()
    x3l,x3u = ax3.get_xlim()
    y3l,y3u = ax3.get_ylim()
    # vertiacl R2, R3 lines on all plots
    #y1l,y1u = ax1.get_ylim()
    ax1.set_ylim(-90,90)
    ax1.hlines(PA,x1l+0.1,x1u-0.1,'k','solid')
    #y2l,y2u = ax2.get_ylim()
    ax2.hlines(ELL,x1l+0.1,x1u-0.1,'k','solid')

    ax3.hlines(0,x3l+0.1,x3u-0.1,'k','solid')
    ax3.hlines(ss,x3l+0.1,x3u-0.1,'k','dashed')

    ax1.set_ylim(-90,90)
    ax2.set_ylim(0,1)
    #if band == 'k':
    #  pl.subplots_adjust(left=0.05, right = 0.95)
    ellplotname = '%s%s_ell'%(ellpath,name)
    pl.savefig(ellplotname+'.png',dpi=160)
    pl.savefig(ellplotname+'.eps',dpi=160)
    if DISP and col=='k':
      os.system('xv %s &'%(ellplotname))
      pl.clf()
 


  return ell,[PA,PAs],[ELL,ELLs],R1s,R2s

def ellip_major(x,y,x0,y0,pa,ba):
  '''
  return semi-major axis length of position (x,y) on an ellipse defined by (x0,y0,eps,pa)
  pa is ccw from +x
  '''
  phi = pa * pl.pi/180.  # convert pa to radians
  dx = x-x0
  dy = y-y0
  cp = pl.cos(phi)
  sp = pl.sin(phi)
  a = pl.sqrt(( dx*cp - dy*sp )**2. + (( dx*sp + dy*cp )/ba)**2.)
  return a

def modelcor(modelname,(nx,ny),x0,y0,pa,ba,C,rout):
  x0 -= 1
  y0 -= 1
  pa = pa*-1
  field = pl.zeros((ny,nx))
  print 'generating model galaxy, this may take a while'
  print 'pa = ',pa
  print 'ba = ',ba
  r,c = field.shape
  x1 = int(max(0,x0-rout))
  x2 = int(min(r,x0+rout))
  y1 = int(max(0,y0-rout))
  y2 = int(min(c,y0+rout))
  if nx != ny:
    x1 = 0
    x2 = r
    y1 = 0
    y2 = c
    
  print 'x0,y0 = ',x0,y0
  print 'nx,ny = ',nx,ny
  print 'r,c = ',r,c
  print 'x1,x2,y1,y2 = ',x1,x2,y1,y2
  for xi in range(x1,x2):
    for yi in range(y1,y2):
      ai = ellip_major(c-1-yi,r-1-xi,x0,y0,pa,ba)
      #if ai < rout:
      val = fitfun(C,ai)
      if val <= 0: field[r-1-xi,c-1-yi] += val
  #modelgal((nx,ny),Gx,Gy,Gpa,Gba,GI0,Grh)
  hdu=pf.PrimaryHDU(field)
  hdulist = pf.HDUList([hdu])
  if os.path.isfile(modelname): os.system('rm '+modelname)
  hdulist.writeto(modelname)
  return 



def fitfun(c,x):
  #y = pl.zeros(len(x))
  if len(c)==3:
    y = c[0]*pl.exp(-x/c[1])+c[2]
  elif len(c)==2:
    #simple linear fit
    y = c[1]*x+c[0]
  else:
    try: y = pl.zeros(len(x))
    except: y = 0.
  return y
def fitfunres(c,x,y):
  yres = y-fitfun(c,x)
  return yres

def dithercor(r,I,fun='lin'):
  flag = 'o'
  #measue outer value
  Isky = pl.average(I[-11:-10])
  #make sure outer value goes to zero
  I = I-Isky
  Is = smooth(I,20,'hanning')

  dI = (I[1:]-I[:-1])/(r[1:]-r[:-1])
  dIs = (Is[1:]-Is[:-1])/(r[1:]-r[:-1])

  #i1 = 0
  #i2 = len(I)-1
  #fl1 = False
  #fl2 = False
  #for i in range(len(dIs)-1):
    #if dIs[i]*dIs[i+1] < 0:
      #if not fl1:
	#i1 = i
      #else:
	#i2 = i
      #break
  #rin = r[i1]
  #rout = r[i2]
  #rin_ind = pl.sum(r < rin)
  #rout_ind = pl.sum(r < rout)

  #find minumum and rminumum
  Imin = I.min()
  Imin_ind = pl.find(I==Imin)
  Ismin = Is.min()
  Ismin_ind = pl.find(Is==Ismin)
  rin = r[Ismin_ind]
  rin_ind = pl.sum(r < rin)
  rout = r[-1]
  rout_ind = pl.sum(r < rout)

  #rin = float(raw_input('Enter inner fitting radius: '))
  #print pl.sum(r < rin)
  #rin_ind = pl.sum(r < rin)
  #rout = float(raw_input('Enter outer fitting radius: '))
  #print pl.sum(r < rout)
  #rout_ind = pl.sum(r < rout)

  rmin = rin

  print 'range: ',Ismin_ind, rin_ind, rout_ind

  rf = r[rin_ind:rout_ind]
  If = Is[rin_ind:rout_ind]

  Itest = Is[(rin_ind+rout_ind)/2:-1]

  dIsig = pl.average(Itest)-3*pl.std(Itest)
  print 'Iav, Isig: ',pl.average(Itest),pl.std(Itest)
  print 'dI, Imin: ',dIsig, Ismin
  #pl.plot(rf,If)
  #pl.show()
  #rmin = r[Imin_ind]
  #fit range: rmin to end
  #rf = r[Imin_ind:rout_ind]
  #If = I[Imin_ind:rout_ind]
  if fun == 'exp':
    c0 = [-2,rmin,0]
  elif fun == 'lin':
    c0 = [-2,0.01]

  if Ismin < dIsig:
    #fit
    try:
      sol,cov,info,msg,success = opt.leastsq(fitfunres,c0,args=(rf,If),full_output=1)
      C = 0.8*sol
      #C = C*0.75
    except:
      print 'ERROR in fit'
      C = []
      flag = 'e'
  else:
    print 'too flat'
    C = []
    flag = 'n'

  return C,rin,rout,flag

def dithercorfun(image,x,y,ell,pa,s,ss,magzp,minrad=5,maxrad='INDEF',col='k',name=''):
  outtbl = tmppath+'tmp_ellip'
  outtbldat = tmppath+'tmp_ellip.dat'
  if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
  if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))

  #run ellipse fitting
  set_ellipse_def()  # set ellipse defaults
  startrad = 0.5*maxrad
  #print '########## MAXRAD = ',maxrad
  #print '########## STARTRAD = ',startrad
  if ell < 0.05: ell = 0.05
  iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')#
  #print '### iraf.ellipse(',im,',',outtbl,',x0=',x,'y0=',y,',ellip0=',ell,',pa0=',pa,',minsma=',minrad,',maxsma=',maxrad,',hellip=yes,hpa=yes)'#,sma0=startrad
  iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,TFLUX_E')
  # plot output
  ell = ellopen(outtbldat)
  #ell[1] -= s
  rad = ell[0]*ps

  intpersqarc = ell[1]/area
  interr = abs(ell[2]/area)

  rad = pl.ma.masked_where(pl.isnan(intpersqarc),rad).compressed()
  intpersqarc = pl.ma.masked_where(pl.isnan(intpersqarc),intpersqarc).compressed()
  interr = pl.ma.masked_where(pl.isnan(intpersqarc),interr).compressed()

  rad = pl.ma.masked_where(pl.isnan(interr),rad).compressed()
  intpersqarc = pl.ma.masked_where(pl.isnan(interr),intpersqarc).compressed()
  interr = pl.ma.masked_where(pl.isnan(interr),interr).compressed()

  #C=dithercor(rad,intpersqarc)

  dr = rad[1:]
  dI = (intpersqarc[1:]-intpersqarc[:-1])/(rad[1:]-rad[:-1])
  

  #pl.figure()
  #pl.errorbar(dr,dI,fmt='.',color='r')

  #pl.xlabel('radius [arcsec]')
  #pl.ylabel('delta Intensity per sq arcsec [cnts]')
  ##pl.ylim([s/area-50,s/area+50])
  #pl.minorticks_on()
  #ddcorplotname = '%s%s_deriv.png'%(ellpath,name)
  #pl.savefig(ddcorplotname,dpi=100)
  #if DISP:
    #os.system('xv '+ddcorplotname+' &')

  #pl.figure()
  #pl.errorbar(rad,intpersqarc,interr,fmt='.',color='r')

  #pl.xlabel('radius [arcsec]')
  #pl.ylabel('Intensity per sq arcsec [cnts]')
  #pl.ylim([s/area-50,s/area+50])
  #pl.minorticks_on()
  #dcorplotname = '%s%s.png'%(ellpath,name)
  #pl.savefig(dcorplotname,dpi=100)
  #if DISP:
    #os.system('xv '+dcorplotname+' &')

  #pl.show()

  #confirm = raw_input('Apply correction? (y) ')
  confirm = 'y'
  if confirm == 'n':
    C =[]
    ru = rad[-1]
  else: 
    fitcont = True
    func = 'lin'
    while fitcont:
      C,rl,ru,flg=dithercor(rad,intpersqarc,func)
      print 'rl, ru = ',rl,ru
      pl.figure(figsize=(5.5,3.5))
      ax1 = pl.subplot(111)
      pl.errorbar(rad,intpersqarc,interr,fmt='.',color='r')
      pl.text(0.95,0.95,r'%s' %(flg),transform=ax1.transAxes)
      pl.xlabel('radius [arcsec]')
      pl.ylabel('Intensity per sq arcsec [cnts]')
      pl.ylim([s/area-50,s/area+50])
      pl.minorticks_on()
      pl.errorbar(rad,intpersqarc- fitfun(C,rad),interr,fmt='.',color='b')
      xl,xu=pl.xlim()
      pl.hlines(s/area,xl,xu)
      pl.plot(rad,fitfun(C,rad)+s/area)
      pl.vlines(rl,s/area-50,s/area+50)
      pl.vlines(ru,s/area-50,s/area+50)
      pl.ylim([s/area-50,s/area+50])

      pl.subplots_adjust(left = 0.15, bottom = 0.12)
      dcorplotname = '%s%s'%(ellpath,name)
      pl.savefig(dcorplotname+'.png',dpi=200)
      pl.savefig(dcorplotname+'.eps',dpi=200)
      if DISP:
	os.system('xv '+dcorplotname+' &')

      pl.clf()
      #confirmcor = raw_input('Is the correction ok? (y/lin/q) ')
      confirmcor = 'y'
      if confirmcor == 'q':
	print 'DIE'
	C = []
	fitcont = False
      elif confirmcor =='y':
	fitcont = False
      else:
	func = confirmcor
      if len(C)>0 :
	if C[1] < 0:
	  C = []
	  fitcont = False

  return C,ru

def radprof(image,x,y,ell,pa,s,ss,magzp,minrad='INDEF',maxrad='INDEF',col='k'):
  '''
radprof(image,x,y,ell,pa,s,ss,magzp,startrad='INDEF',minrad='INDEF',maxra
d='INDEF')
  '''
  outtbl = tmppath+'tmp_ellip'
  outtbldat = tmppath+'tmp_ellip.dat'
  if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
  if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))

  #run ellipse fitting
  set_ellipse_def()  # set ellipse defaults
  startrad = 0.5*maxrad
  #print '########## MAXRAD = ',maxrad
  #print '########## STARTRAD = ',startrad
  if ell < 0.05: ell = 0.05

  try:
    print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(image,outtbl,x,y,ell,pa,startrad,minrad,maxrad)
    iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')

  except:
    iraf.samplepar.nclip = 0

    print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(image,outtbl,x,y,ell,pa,startrad,minrad,maxrad) 
    iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')
#
  #print '### iraf.ellipse(',im,',',outtbl,',x0=',x,'y0=',y,',ellip0=',ell,',pa0=',pa,',minsma=',minrad,',maxsma=',maxrad,',hellip=yes,hpa=yes)'#,sma0=startrad
  iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,TFLUX_E')
  # plot output
  ell = ellopen(outtbldat)
  ell[1] -= s
  rad = ell[0]*ps

  intpersqarc = ell[1]/area

  interr = abs(ell[2]/area)
  sb = -2.5*pl.log10(intpersqarc)+magzp

  sig = -2.5*pl.log10(ss) + magzp
  i = 0
  while sb[i] < sig:
    if i >= len(sb)-1:
      break
    i += 1
  I1 = i
  #I1 = pl.sum(sb<sig)
  if I1 == len(sb): I1 -= 1
  r1s = rad[I1]
  #print 'r1s = ',r1s
  #print 'I1 = ',I1
  Rout = 1.25*r1s
  #print ' Rout = %.1f arcsec' %(Rout)
  I = pl.sum(rad<Rout)
  if I == len(sb): I -= 1

  sberru = -2.5*pl.log10((intpersqarc+interr))+magzp
  sberrl = -2.5*pl.log10((intpersqarc-interr))+magzp
  #restrict output range to R1s for now we want the full range
  rad = rad[:I]
  sb = sb[:I]
  sberru = sberru[:I]
  sberrl = sberrl[:I]
  sberr = [sb-sberrl,sberru-sb]


  t = open('tempSBP','w')
  for it in range(len(rad)):
    t.write(str(rad[it])+'\t'+str(ell[1][it])+'\t'+str(ell[2][it]) +'\n' )
  t.close()

  return r1s,rad,[sb,sberr]

#def radprof(image,ACOR,x,y,ell,pa,s,ss,magzp,minrad='INDEF',maxrad='INDEF',col='k'):
  #'''
#radprof(image,x,y,ell,pa,s,ss,magzp,startrad='INDEF',minrad='INDEF',maxra
#d='INDEF')
  #'''
  #outtbl = tmppath+'tmp_ellip'
  #outtbldat = tmppath+'tmp_ellip.dat'
  #if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
  #if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))

  ##run ellipse fitting
  #set_ellipse_def()  # set ellipse defaults
  #startrad = 0.5*maxrad
  ##print '########## MAXRAD = ',maxrad
  ##print '########## STARTRAD = ',startrad
  #if ell < 0.05: ell = 0.05
  #iraf.ellipse(image,outtbl,x0=x,y0=y,ellip0=ell,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')#
  ##print '### iraf.ellipse(',im,',',outtbl,',x0=',x,'y0=',y,',ellip0=',ell,',pa0=',pa,',minsma=',minrad,',maxsma=',maxrad,',hellip=yes,hpa=yes)'#,sma0=startrad
  #iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,TFLUX_E')
  ## plot output
  #ell = ellopen(outtbldat)
  #ell[1] -= s
  #rad = ell[0]*ps

  #intpersqarc = ell[1]/area

  #intpersqarc = intpersqarc - fitfun(ACOR,rad)

  #interr = abs(ell[2]/area)
  #sb = -2.5*pl.log10(intpersqarc)+magzp

  #sig = -2.5*pl.log10(ss) + magzp
  #I1 = pl.sum(sb<sig)
  #if I1 == len(sb): I1 -= 1
  #r1s = rad[I1]
  ##print 'r1s = ',r1s
  ##print 'I1 = ',I1
  #Rout = 1.25*r1s
  ##print ' Rout = %.1f arcsec' %(Rout)
  #I = pl.sum(rad<Rout)
  #if I == len(sb): I -= 1

  #sberru = -2.5*pl.log10((intpersqarc+interr))+magzp
  #sberrl = -2.5*pl.log10((intpersqarc-interr))+magzp
  ## restrict output range to R1s for now we want the full range
  ##rad = rad[:I]
  ##sb = sb[:I]
  ##sberru = sberru[:I]
  ##sberrl = sberrl[:I]
  #sberr = [sb-sberrl,sberru-sb]


  #t = open('tempSBP','w')
  #for it in range(len(rad)):
    #t.write(str(rad[it])+'\t'+str(ell[1][it])+'\t'+str(ell[2][it]) +'\n' )
  #t.close()

  #return r1s,rad,[sb,sberr]


def sbfitform(x,c):
  if len(c)==2:
    #print 'e'
    c1,c2 = c
    return c1+1.0857*((x/c2))
  elif len(c)==3:
    #print 's'
    c1,c2,c3 = c
    return c1+1.0857*((x/c2)**c3)
  elif len(c)==5:
    ##print '2e'
    c1,c2,c3,c4,c5 = c
    return -2.5*pl.log10(10**(-0.4*(c1+1.0857*((x/c2)**c3))) + 10**(-0.4*(c4+1.0857*(x/c5))))
    #print 't'
  elif len(c)==6:
    #print '2e'
    c1,c2,c3,c4,c5,c6 = c
    return -2.5*pl.log10(10**(-0.4*(c1+1.0857*((x/c2)**c3))) + 10**(-0.4*(c4+1.0857*((x/c5)**c6))))
  else: return x

def wresidual(c,x,y,yerr):
  return (y-sbfitform(x,c))/(yerr)

def fitsbindiv(sma,sb,sberr,c0):
  '''
fitsb - fits sersic function to sb profile
  '''
  flag = 0
  ncoeff = len(c0)
  try:
    sol,cov,info,msg,success = opt.leastsq(wresidual,c0,args=(sma,sb,sberr),full_output=1)
    C = sol
    C_err = pl.zeros(ncoeff)
    #chi2 = pl.sum(wresidual(C,sma,sb,sberr))
    chisq = pl.sum(info["fvec"]*info["fvec"])
    dof=len(sma)-len(sol)
    rms = pl.sqrt(chisq/dof)
    if cov is not None:
      for i in range(ncoeff):
	C_err[i] = pl.sqrt(cov[i][i]*chisq/dof)
  except:
    #print sb,sberr
    C = c0
    C_err = pl.zeros(len(c0))
    rms = 999
    flag = 1
  #for i,pmin in enumerate(sol):
  #  print "%2i %-10s %12f +/- %10f" %(i,'c',pmin,pl.sqrt(cov[i,i])*pl.sqrt(chisq/dof)) 

  return C, C_err, rms, flag

def fitsb(sma,sb,rlim,col='k'):
  '''
fitsb - fits sersic function to sb profile
  '''
  ce = [20.0,20.4]
  cs = [20.0,20.4,0.8]
  cse = [20.0,15.0,1.0,22.0,15.0]
  c2e = [20.0,15.0,1.0,22.0,15.0,1.0]

  i1 = pl.sum(sma<rlim[0])
  i2 = pl.sum(sma<rlim[1])
  i1 = (i1+i2)/2
  sberr = pl.sqrt(sb[1][0]**2.+sb[1][1]**2)
  Cs = []
  Cs_err = []
  chis = []
  fitflg = []
  for c0 in [ce,cs,cse,c2e]:
    c,c_err,chi,flg = fitsbindiv(sma[i1:i2],sb[0][i1:i2],sberr[i1:i2],c0)
    Cs += [c]
    Cs_err += [c_err]
    chis += [chi]
    fitflg += [flg]
  chis=pl.array(chis)
  fitflg=pl.array(fitflg)
  print 'chisq: ',chis
  temp = chis + fitflg*100
  print 'chisq flagged: ',temp
  it = pl.find(temp==temp.min())[0]
  print 'SB fits: ',chis,it

  ind = 0
  while sb[0][ind] < 21.:
    ind += 1
    if ind >= len(sb[0])-1:
      ind -= 1
      break
  rout = sma[ind]

  coeff = Cs[it]
  ## check value of dm!!
  if len(coeff)==2:
    coeff = pl.append(coeff,1.)
  #print len(coeff)
  if len(coeff)==3:
    n = coeff[2]
    rh = coeff[1]
    Imis = spec.gammainc(2./n,(rout/rh)**n)
    Itot = spec.gamma(2./n)
  if len(coeff)==5:
    coeff = pl.append(coeff,1.)
  if len(coeff)==6:
    I0 = pl.log10(coeff[0])
    n0 = coeff[2]
    rh0 = coeff[1]
    I1 = pl.log10(coeff[3])
    n1 = coeff[5]
    rh1 = coeff[4]
    Imis = (rh0**2.*I0/n0)*spec.gammainc(2./n0,(rout/rh0)**n0) + (rh1**2.*I1/n1)*spec.gammainc(2./n1,(rout/rh1)**n1)
    Itot = (rh0**2.*I0/n0)*spec.gamma(2./n0)+ (rh1**2.*I1/n1)*spec.gamma(2./n1)
  dm = -2.5*pl.log10(Imis/Itot)
  # if too big fit exponential to outer disk
  if dm > 1:
    i1 = pl.sum(sma<rlim[0])
    i2 = pl.sum(sma<rlim[1])
    i1 = (i1+2*i2)/3
    Cs = []
    Cs_err = []
    chis = []
    fitflg = []
    for c0 in [ce]:
      c,c_err,chi,flg = fitsbindiv(sma[i1:i2],sb[0][i1:i2],sberr[i1:i2],c0)
      Cs += [c]
      Cs_err += [c_err]
      chis += [chi]
      fitflg += [flg]
    chis=pl.array(chis)
    fitflg=pl.array(fitflg)
    print 'chisq: ',chis
    temp = chis + fitflg*100
    print 'chisq flagged: ',temp
    it = pl.find(temp==temp.min())[0]
    print 'SB fits: ',chis,it

  return Cs[it], Cs_err[it], chis[it], fitflg[it]
  print 'sersic index=%s'% Cs[it]

def plotsb(fitsimage,sky,ss,sig,galx,galy,ell,pa,sma,sb,coeff,coeff_err,goodfit='',name='',band='',Col='k'):

  ellplotname = '%s%s_%s'%(ellpath,name,band)

  # add to plot
  fig = pl.figure(figsize=(4,6))
  fig.suptitle(band+' Surface Brightness\n'+name)
  ax1 = fig.add_axes([0.12,0.4,0.8,0.5]) #left,bottom,width,height

  ax2 = fig.add_axes([0.12,0.12,0.8,0.2],sharex=ax1)
  ax1.set_ylabel('$\mu$ [mag arcsec$^{-2}$]')
  ax2.set_ylabel(r'residual')
  ax2.set_xlabel('$a$ [arcsec]')

  v1 = sky-3*ss
  data = pf.getdata(fitsimage)
  boxsize=10
  col,row = data.shape
  x1=(row/2-1)-boxsize
  x2=(row/2-1)+boxsize
  y1=(col/2-1)-boxsize
  y2=(col/2-1)+boxsize
  #sub=data[x1:x2,y1:y2]
  v2=data[x1:x2,y1:y2].max()  #  peak pixel
  #image
  #ax3.imshow(pl.log10(data),cmap=cmap.gray)

  ax3 = ap.FITSFigure(fitsimage,figure=fig,subplot=[0.675,0.7,0.25,0.2])
  ax3.hide_axis_labels()
  ax3.hide_tick_labels()
  ax3.show_grayscale(stretch='log',invert=True,vmin=v1,vmax=v2)
  galra,galdec=ax3.pixel2world(galx,galy)


  xmin,xmax=ax2.get_xlim()
  ax1.hlines(sig,0,xmax,color='k',linestyles='dashed')
  # data and fit
  print len(sma)
  print len(sb[0])
  print len(sb[1][0])
  print len(sb[1][1])

  if len(sma)>0:
    ax1.errorbar(sma,sb[0],sb[1],fmt=Col+',',capsize=0)

    SBfit = sbfitform(sma,coeff)
    ax1.plot(sma,SBfit,Col)
    ax1.set_ylim(22,12)

    #residuals
    ax2.errorbar(sma,SBfit-sb[0],sb[1],fmt=Col+',',capsize=0)
    xmin,xmax=ax2.get_xlim()
    ymin,ymax=ax2.get_ylim()
    ax1.hlines(sig,0,xmax,color='k',linestyles='dashed')
    ax2.set_ylim(-1.,1.)
    ax2.hlines(0.,0.,xmax,'k')

    it = 0
    while sb[0][it] < sig:
      it += 1
      if it >= len(sb[0])-1:
	it -= 1
	break
    I1 = it

    #I1 = pl.sum(sb[0]<sig)
    if I1 == len(sb[0]): I1 -= 1
    r1s = sma[I1]

    ax3.recenter(galra,galdec,radius=1.5*r1s/3600.)
    #I removed ec='w' from the coming line
    ax3.show_ellipses([galra], [galdec], r1s/3600., (ell[0]-1)*r1s/3600., angle=pa[0]+90.)


    ax1.vlines(r1s,26,12,color='k',linestyles='dashed')
    ax2.vlines(r1s,-0.99,0.99,color='k',linestyles='dashed')

    # text showing fit params
    ctext = ['$\mu_0$','$r_h$','$n$','$\mu_{02}$','$r_{h2}$','$n_2$']
    ofile = open('%s_sersicn.txt'%ellplotname, 'w')
    for i in range(len(coeff)):
      pl.text(0.0125,0.92-i*0.07,r'%s = $%.1f \pm %.1f $' %(ctext[i],coeff[i],coeff_err[i]),transform=ax1.transAxes)
      print 'surface brightness index=%.1f' % coeff[i]
    ##### output file
    
      print >>ofile, ctext[i], coeff[i],coeff_err[i], '\n'
    ofile.close()

    pl.text(0.0125,0.83,r'$rms$ = $%.1f$' %(goodfit),transform=ax2.transAxes)

  else:
    print 'ERROR in plotsb'
  pl.text(0.625,0.55,r'$\phi$ = $%.1f \pm %.1f$'%(pa[0],pa[1]),transform=ax1.transAxes)
  pl.text(0.625,0.5,r'$b/a$ = $%.1f \pm %.1f$'%(ell[0],ell[1]),transform=ax1.transAxes)

  ax1.set_ylim(26,12)
  ax2.set_ylim(-1.,1.)
  ax1.set_xlim(xmin=0)
  fig.canvas.draw()

  pl.savefig(ellplotname+'.png',dpi=250)
  pl.savefig(ellplotname+'.eps',dpi=250)
  if DISP: os.system('xv %s &'%(ellplotname))
  pl.clf()
  return 

def elphot(im,Riso,pa,ellip,x,y,s,ss,mzp):
  outtbl = tmppath+'tmp_ellip'
  outtbldat = tmppath+'tmp_ellip.dat'
  if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
  if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))

  set_ellipse_def()  # set ellipse defaults

  # mag pars: m = mag0 - 2.5log[(intens-xerolevle)/refer]
  iraf.magpar.mag0 = 0.
  iraf.magpar.refer = 10.**(mzp/2.5)
  iraf.magpar.zerolevel = s

  startrad = Riso/ps
  #print 'sma0=',startrad
  if startrad < 5.: startrad = 5.
  minrad = startrad
  maxrad = startrad
  if ellip < 0.05: ellip = 0.05
  if ellip > 1.: ellip = 1.

  try:
    try:
      print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(im,outtbl,x,y,ellip,pa,startrad,minrad,maxrad)
      iraf.ellipse(im,outtbl,x0=x,y0=y,ellip0=ellip,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')#
    except:
      iraf.samplepar.nclip = 0
      print 'iraf.ellipse(%s,%s,x0=%f,y0=%f,ellip0=%f,pa0=%f,sma0=%f,minsma=%f,maxsma=%f)' %(im,outtbl,x,y,ellip,pa,startrad,minrad,maxrad) 
      iraf.ellipse(im,outtbl,x0=x,y0=y,ellip0=ellip,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')#

    iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,TFLUX_E,TMAG_E,NPIX_E')
    #os.system('cat '+outtbldat)
    ell = ellopen(outtbldat)
    Npix = ell[5][0]
    Tflux = ell[3][0]
    #isomag = -2.5*pl.log10(ell[3])+mzp
    isosky = Npix*s
    isoflux = (Tflux-isosky)
    isoflux_err = pl.sqrt(isoflux+Npix*ss**2.)
    #print 'isoflux = %.0f +/- %.0f dn' %(isoflux,isoflux_err)
    isomag = -2.5*pl.log10(isoflux) +mzp
    isomag_err = pl.sqrt((2.5*isoflux_err/(isoflux*2.30259))**2.+(sig_mzp)**2.)
  #print 'isomag = %.3f +/- %.3f mag' %(isomag,isomag_err)
  except:
    print 'ERROR in elphot'
    isomag = pl.nan
    isomag_err = pl.nan
  return isomag,isomag_err

def skyann(im,x,y,pa,ellip,radius2):
  F = []
  N = []
  for rad in radius2:
    outtbl = tmppath+'tmp_ellip'
    outtbldat = tmppath+'tmp_ellip.dat'
    if os.path.isfile(outtbl+'.tab'): os.system('rm %s.tab' %(outtbl))
    if os.path.isfile(outtbldat): os.system('rm %s' %(outtbldat))
    set_ellipse_def()  # set ellipse defaults
    startrad = rad/ps
    if startrad < 5.: startrad = 5.
    minrad = startrad
    maxrad = startrad
    if ellip < 0.05: ellip = 0.05
    if ellip > 1.: ellip = 1.
    iraf.ellipse(im,outtbl,x0=x,y0=y,ellip0=ellip,pa0=pa,minsma=minrad,maxsma=maxrad,sma0=startrad,hellip='yes',hpa='yes')#
    iraf.tdump(outtbl,datafile=outtbldat,cdfile='',pfile='',columns='SMA,INTENS,INT_ERR,TFLUX_E,TMAG_E,NPIX_E')
    #os.system('cat '+outtbldat)
    ell = ellopen(outtbldat)
    #print 'npix=',ell[5][0]
    #print 'f=',ell[3][0]
  
    N += [ell[5][0]]
    F += [ell[3][0]]
  Fs = F[1]-F[0]
  Ns = N[1]-N[0]
  s = Fs/Ns
  
  return s

#def isophote(im,x,y,ell,pa,isoval,s,ss,mzp,r,sb,mu0,rh,n):
def isophote(im,x,y,ell,pa,isoval,s,ss,mzp,sma,sb):#,mu0,rh,n):
  # get isophotal index
  #SBfit = sersic(r,mu0[0],rh[0],n[0])
  I1 = pl.sum(sb[0]<isoval)
  N = 5
  yval = sb[0][I1-N:I1+N]
  xval = sma[I1-N:I1+N]
  #print xval
  #print yval
  change = False
  for i in range(len(yval)-1,0,-1):
    if pl.isnan(yval[i]):
      lim = i
      change = True
  if change:
    xval = xval[:lim]
    yval = yval[:lim]
  print xval
  print yval
  try:
    a,b = pl.polyfit(xval,yval,1)
    Riso = (isoval-b)/a
  except:
    if I1 == len(sma): I1 -= 1
    if len(sma)>0: Riso = sma[I1]
    else: Riso = pl.nan

  if I1 == len(sma): I1 -= 1
  if len(sma)> 0: Riso2 = sma[I1]
  else: Riso2 =pl.nan
  print 'R_%.2f = %.2f arcsec' %(isoval,Riso)
  print 'R_%.2f(2) = %.2f arcsec' %(isoval,Riso2)
  if pl.isnan(Riso):
    Riso = Riso2
  if pl.isnan(Riso):
    Riso = Riso2
    Mag,Mag_err = pl.nan,pl.nan
  else:
    Mag,Mag_err = elphot(im,Riso,pa,ell,x,y,s,ss,mzp)

  #pl.figure()
  #pl.plot(ell[0],ell[4])
  return Riso,[Mag,Mag_err]

def mtotal(im,x,y,ell,pa,s,ss,mzp,sma,sb,coeff):

  iso_out=21.
  rout,M = isophote(im,x,y,ell,pa,iso_out,s,ss,mzp,sma,sb)
  miso = M[0]
  miso_err = M[1]
  #print len(coeff)
  if len(coeff)==2:
    coeff = pl.append(coeff,1.)
  #print len(coeff)
  if len(coeff)==3:
    n = coeff[2]
    rh = coeff[1]
    Imis = spec.gammainc(2./n,(rout/rh)**n)
    Itot = spec.gamma(2./n)
  if len(coeff)==5:
    coeff = pl.append(coeff,1.)
  if len(coeff)==6:
    I0 = pl.log10(coeff[0])
    n0 = coeff[2]
    rh0 = coeff[1]
    I1 = pl.log10(coeff[3])
    n1 = coeff[5]
    rh1 = coeff[4]
    Imis = (rh0**2.*I0/n0)*spec.gammainc(2./n0,(rout/rh0)**n0) + (rh1**2.*I1/n1)*spec.gammainc(2./n1,(rout/rh1)**n1)
    Itot = (rh0**2.*I0/n0)*spec.gamma(2./n0)+ (rh1**2.*I1/n1)*spec.gamma(2./n1)
  dm = -2.5*pl.log10(Imis/Itot)
  mt = miso-dm

  mt_err = pl.sqrt(miso_err**2. + (0.1*dm**2.))

  print 'Total mag: %.3f \pm %.3f' %(mt,mt_err)
      #print 'sersic parameters :%s %s' %(n,n0)
  #mterr = miso_err

  return mt, dm, mt_err


def maperture(im,x,y,s,ss,mzp):
  aps = [3.,5.,7.,10.,15.]
  m_ap = pl.zeros(len(aps))
  msig_ap = pl.zeros(len(aps))
  sb_core = 0.

  ellip = 0.05
  pa = 0.0
  for i in range(len(aps)):
    apa = aps[i]
    ap = apa/ps
    head = pyfits.getheader(im)
    data = pyfits.getdata(im)
    r,c = data.shape
    mask = pl.zeros((r,c))
    X = x-1
    Y = y-1
    Npix = 0
    success = True
    print X-ap-2,X+ap+3
    for xi in range(int(X-ap-2),int(X+ap+3)):
      for yi in range(int(Y-ap-2),int(Y+ap+3)):
	if (X-xi)**2.+(Y-yi)**2 < (ap+0.5)**2:
	  try:
	    mask[yi,xi]=1
	  except:
	    success = False
	  Npix += 1

    if success:
      apdat = data*mask 
      Tflux = pl.sum(apdat)
      isosky = Npix*s
      isoflux = (Tflux-isosky)
      isoflux_err = pl.sqrt(isoflux+Npix*ss**2.)
      isomag = -2.5*pl.log10(isoflux) +mzp
      isomag_err = pl.sqrt((2.5*isoflux_err/(isoflux*2.30259))**2.+(sig_mzp)**2.)
      m_ap[i] = isomag
      msig_ap[i] = isomag_err
    else:
      m_ap[i] = pl.nan
      msig_ap[i] = pl.nan
  print 'test_this_place'

    #print '*Npix = ',Npix
    #print '*TF_%i    = %.2f' %(apa,Tflux)
    #print '*F_%i    = %.2f' %(apa,isoflux)
    #print '*m_%i    = %.3f' %(apa,m_ap[i])
    #print '*msig_%i = %.3f' %(apa,msig_ap[i])
  sb_core = m_ap[1]+2.5*pl.log10(25.)
    
  return m_ap,msig_ap,sb_core

# calculate stellar density in field
def calc_field_dens(im):
  os.system('rm tmp*')
  iraf.unlearn(iraf.daophot)
  imgname = 'tmp.fits'
  head = pf.getheader(im)
  nx = head.get('NAXIS1')
  ny = head.get('NAXIS2')
  #exclude outer 200 pixels
  d = 200
  x1 = d
  y1 = d
  x2 = nx-d
  y2 = ny-d
  imcut = '%s[%i:%i,%i:%i]' %(im.replace('.fits',''),x1,x2,y1,y2)
  iraf.imcopy(imcut,imgname)

  ps = 0.45
  psfrad=12

  iraf.datapar.readnoise=147
  iraf.datapar.epadu=132
  iraf.datapar.datamax=30000
  iraf.fitskyp.salgori="mode"
  iraf.fitskyp.annulus=12
  iraf.fitskyp.dannulu=3
  iraf.fitskyp.sloclip=0 
  iraf.fitskyp.shiclip=0
  iraf.daopars.function="auto"
  iraf.daopars.varorder=2  # change from 1
  iraf.daofind.interac='no'
  iraf.daofind.verify='no' 
  iraf.daofind.verbose= 'no'
  iraf.phot.interac='no'
  iraf.phot.verify='no'
  iraf.phot.verbose='no'

  # measure sky and rms from image
  skylevel = float(iraf.imstat(images=imgname, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
  stdev = float(iraf.imstat(images=imgname, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
  iraf.datapars.sigma=stdev
  iraf.datapars.datamin = skylevel - 5*iraf.datapars.sigma

  head = pf.getheader(imgname)
  seeing = head.get('SEEING')
  if seeing is None:
    raise PsfError('seeing')
    #return
    
  fwhm   = seeing / ps
  iraf.datapars.fwhmpsf = fwhm

  # read magnitude zero point (MAGZP)
  magzp = head.get('MAGZP')
  iraf.photpars.zmag = magzp
  if magzp is None:
    raise PsfError('magzp')
    return
  expos = head.get('EXPOS')
  iraf.datapars.itime = expos
  if expos is None:
    raise PsfError('expos')
    return
  nx = head.get('NAXIS1')
  ny = head.get('NAXIS2')
  fov=(nx*ny*ps**2)/(60.**4)

  iraf.photpars.aperture=iraf.datapars.fwhmpsf
  iraf.daopars.fitrad=iraf.datapars.fwhmpsf

  iraf.daofind(image=imgname, output="tmp1.coo.1", threshold=3., verify='no')
  #iraf.pselect(infiles="tmp1.coo.1", outfiles="tmp1.coo.2", expr="SHARPNESS != INDEF")
  iraf.phot(image=imgname, coords="tmp1.coo.1", output="tmp1.mag.1", verify='no')
  iraf.pselect("tmp1.mag.1", "tmp1.mag.2", expr="MAG != INDEF")
  iraf.pdump("tmp1.mag.2", "MAG", 'yes', Stdout="tmp1.mag.1.pdump")

  mag=pl.loadtxt('tmp1.mag.1.pdump')
  
  N14 = pl.log10(pl.sum(mag<14)/fov)
  N16 = pl.log10(pl.sum(mag<16)/fov)
  print 'N14 = %.3f ' %(N14)
  print 'N16 = %.3f ' %(N16)
  return N14,N16

# need hicat info ra0, dec0, a, b, Gpa
def photloop(Jimage_nir,ra0,dec0,a,b,Gpa):
  Ga = 0.3*a    # initial ellipse fit
  Gell = 1-b/a
  print
  print ' ### %s ###' %(Jimage_nir)
  os.system('rm %s*' %(tmppath))
  print 'INITIAL ELLIPSE'
  print ' # ra0, dec0 = (%f,%f)' %(ra0,dec0)
  print ' # a0, pa0, ell0 = %.2f, %.1f, %.2f ' %(a,Gpa,Gell)
  #iraf.geompar.maxsma = 2*a

  Jim = imagepath+'gal_j'+Jimage_nir+'_starsub.fits'
  Him = imagepath+'gal_h'+Jimage_nir+'_starsub.fits'
  Kim = imagepath+'gal_k'+Jimage_nir+'_starsub.fits'

  SD14,SD16 = calc_field_dens(imagepath+'k'+Jimage_nir.split('_')[0]+'.fits')
  indebv = afind(Jimage_nir.split('_')[0],hidat['NAME'])[0]
  Ebv = hidat['E(B-V)'][indebv]

  print 'SD = %.3f' %(SD14)
  print 'Ebv = %.3f' %(Ebv)

  # CHECK IF NEAR EDGES
  Jbig = imagepath+'j'+Jimage_nir.split('_')[0]+'.fits'
  head = pf.getheader(Jbig)
  wcs = pywcs.WCS(head)
  blub=wcs.wcs_sky2pix(pl.array([ra0]),pl.array([dec0]),1)
  x0=float(blub[0])
  y0=float(blub[1])
  print x0,y0
  nx=head.get('NAXIS1')
  ny=head.get('NAXIS2')
  Photflag = 0
  dE = 150
  if (x0 < dE) or (y0 < dE) or (x0 > nx-dE) or (y0 > ny-dE):
    Photflag += 1
    print 'EDGE!!!'
  if a > 500:		# larg gal
    Photflag += 2

  if not os.path.isfile(Jim):
    print 'missing image %s' %(Jim)
    return
  if not os.path.isfile(Him):
    print 'missing image %s' %(Him)
    return
  if not os.path.isfile(Kim):
    print 'missing image %s' %(Kim)
    return

  # add J+H+K = SUP image - super co-add
  SUPim = tmppath+'sup.fits'
  iraf.imarith(Jim, "+", Him, SUPim)
  iraf.imarith(SUPim, "+", Kim, SUPim)

  #get image coords in all bands - SHOULD be the same as images are pix registered
  head = pf.getheader(Jim)
  wcs = pywcs.WCS(head)
  blub=wcs.wcs_sky2pix(pl.array([ra0]),pl.array([dec0]),1)
  x0=float(blub[0])
  y0=float(blub[1])
  nx=head.get('NAXIS1')
  ny=head.get('NAXIS2')
  print ' # x0, y0 = (%.2f,%.2f)' %(x0,y0)

  #find the galaxy center (centroid)
  frad = a
  if frad < 25: frad = 25
  x,y=findgalcoor(SUPim,x0,y0,frad)
  print ' # X, Y = (%.2f,%.2f)' %(x,y)
  #get ra/dec, glat/glon
  blub = wcs.wcs_pix2sky(pl.array([x]),pl.array([y]),1)
  ra = float(blub[0])
  dec = float(blub[1])
  print ' # RA, DEC = (%.5f,%.5f)' %(ra,dec)
  glat,glon = ac.convertCoords('J2000','GALACTIC',ra,dec,2000)
  print ' # GLAT, GLON = (%.5f,%.5f)' %(glat,glon)
  ellfig = pl.figure(figsize=(10.0,10.0))  #pa/ellipticity fig
  print "test_1"  
  #initialise arrays
  sky_vals = pl.zeros([3,2])
  sky_marc = pl.zeros(3)
  magzps = pl.zeros(3)
  sees = pl.zeros(3)
  exposs = pl.zeros(3)
  posang = pl.zeros([3,2])
  ellipticity = pl.zeros([3,2])
  r20 = pl.zeros(3)
  m20 = pl.zeros([3,2])
  r21 = pl.zeros(3)
  m21 = pl.zeros([3,2])
  r22 = pl.zeros(3)
  m22 = pl.zeros([3,2])
  sb_core = pl.zeros(3)
  peak = pl.zeros(3)
  mT = pl.zeros(3)
  mT_err = pl.zeros(3)
  dmT = pl.zeros(3)
  SBflag = pl.zeros(3)
  m_ap = pl.zeros([3,5])
  msig_ap = pl.zeros([3,5])
  print "test_2"
  images = [tmppath+'j.fits',tmppath+'h.fits',tmppath+'k.fits']
  images2 = images
  print "test_3"
  for i in range(3):		#TEMP k only
    band = bands[i]
    lfmt = lfmts[i]  #line formats
    F = str(i+1)
    print 
    print '## running %s%s photometry ##' %(band,Jimage_nir)


    # image
    im_in = imagepath+'gal_'+band+Jimage_nir+'_starsub.fits'

    xband,yband=findgalcoor(im_in,x0,y0,frad)

    #check images exist
    if not os.path.isfile(im_in):
      print 'no %s image: %s' %(band,im_in)
      break
    im = images[i]  # make temp image

    os.system('cp %s %s' %(im_in,im))
    if DISP:
      ds9.set('frame '+F)
      ds9.set('file '+im)
      wcs_mark(ds9,x,y,colour='black')
      wcs_mark(ds9,x,y,a,b,Gpa+90.)


    head = pf.getheader(im)
    magzpp = head.get('MAGZP')  # for counts/sec <-> mag
    see = head.get('SEEING')
    sees[i] = see
    expos = head.get('EXPOS')
    exposs[i] = expos

    magzp = magzpp + 2.5*pl.log10(expos)    #for total counts <-> mag

    # measure sky in postage
    sky=float(iraf.imstatistics(im,fields='mean',nclip=30,format=0,Stdout=1)[0])
    skysig=float(iraf.imstatistics(im,fields='stddev',nclip=30,format=0,Stdout=1)[0])
    sky_vals[i][0] = sky
    sky_vals[i][1] = skysig

    #iraf.imarith(im, "-", sky, im0)

    sig = -2.5*pl.log10(skysig) + magzp
    print ' %s 1sig_SB = %.3f mag/arcsec^2' %(band,sig)
    sky_marc[i] = sig
    magzps[i] = magzp

    #peak pixel
    data = pf.getdata(im)
    box = pl.array(data[int(x-10):int(x+10)+1,int(y-10):int(y+10)+1])
    peak[i] = -2.5*pl.log10(box.max()) + magzp

    ell1,PA,ELL,R1s,R2s = fit_ell(Jimage_nir,ellfig.number,band,im,sky,skysig,xband,yband,Gell,Gpa,startrad=Ga,minrad=see/ps,maxrad=3*a,col=lfmt)
    rad=ell1[0]*ps

    posang[i] = PA
    ellipticity[i] = ELL

    # remeasure sky in annulus
    #print 'IMAGE: sky, skysig = ',sky,skysig
    #if DISP:
      #wcs_mark(ds9,x,y,1.5*R1s,(ELL[0]-1)*1.5*R1s,PA[0]+90.,colour='magenta')
      #wcs_mark(ds9,x,y,2.5*R1s,(ELL[0]-1)*2.5*R1s,PA[0]+90.,colour='magenta')
    #sky2 = skyann(im,x,y,PA[0],ELL[0],[1.5*R1s,2.5*R1s])
    #print 'ANNULUS: sky = ',sky2,skysig
    #dsky=sky-sky2
    #print 'DIFFERENCE: dsky = ',dsky
    #skyf=sky-3*dsky

    # big galaxies: do dither correction
    if (R1s > 10) and (R1s*(1-ELL[0]) > 10):
	#try:
	Dcorrect, Rout = dithercorfun(im,x,y,ELL[0],PA[0],sky,skysig,magzp,minrad=0.5*see/ps,maxrad=6*R1s/ps,col=lfmt,name=Jimage_nir+'_dcor_'+band)
    
	if len(Dcorrect) >0:
	  if Dcorrect[0] > -25:
	    modelcor(tmppath+band+'dcormod.fits',(nx,ny),xband,yband,PA[0]+90,1-ELL[0],Dcorrect,Rout)
	    if os.path.isfile(tmppath+band+'dcor.fits'): os.system('rm '+tmppath+band+'dcor.fits')
	    iraf.imarith(im,'-',tmppath+band+'dcormod.fits',tmppath+band+'dcor.fits')
	    im = tmppath+band+'dcor.fits'
	    images2[i] = im
	    if DISP:
	      ds9.set('frame %i' %(3+int(F)))
	      ds9.set('file '+tmppath+band+'dcor.fits')
	#except:
	#print 'ERROR doing dither correction'
    else:
      Dcorrect=[]

    #print 'R1s = ',R1s
    startrad = 3*see/ps

    R1s,rad,SB = radprof(im,x,y,ELL[0],PA[0],sky,skysig,magzp,minrad=0.5*see/ps,maxrad=3*R1s/ps,col=lfmt)
    # mask nans
    rad = pl.ma.masked_where(pl.isnan(SB[0]),rad).compressed()
    SBdat = pl.ma.masked_where(pl.isnan(SB[0]),SB[0]).compressed()
    SBerr1 = pl.ma.masked_where(pl.isnan(SB[0]),SB[1][0]).compressed()
    SBerr2 = pl.ma.masked_where(pl.isnan(SB[0]),SB[1][1]).compressed()
    #for itr in range(len(rad)):
    #print '%f\t%f\t%f\t%f' %(rad[itr],SBdat[itr],SBerr1[itr],SBerr2[itr])

    SB = [SBdat,[SBerr1,SBerr2]]

    #R1s,rad,SB = radprof(im,Dcorrect,x,y,ELL[0],PA[0],sky,skysig,magzp,minrad=0.5*see/ps,maxrad=3*R1s/ps,col=lfmt)
    sbc,sbc_err,chi,SBflag[i]=fitsb(rad,SB,[see,R1s],col=lfmt)
    print 'sersic index=%s' % sbc

    plotsb(im,sky,skysig,sig,x,y,ELL,PA,rad,SB,sbc,sbc_err,goodfit=chi,name=Jimage_nir,band=band,Col=lfmt)

    #iraf.imarith(im,'-',skysig,im)
    r20[i],m20[i] = isophote(im,x,y,ELL[0],PA[0],20.,sky,skysig,magzp,rad,SB)
    print 'R20 = %.2f arcsec' %(r20[i])
    print 'm20 = %.2f mag' %(m20[i][0])
    r21[i],m21[i] = isophote(im,x,y,ELL[0],PA[0],21.,sky,skysig,magzp,rad,SB)
    r22[i],m22[i] = isophote(im,x,y,ELL[0],PA[0],21.,sky,skysig,magzp,rad,SB)

    mT[i],dmT[i],mT_err[i] = mtotal(im,x,y,ELL[0],PA[0],sky,skysig,magzp,rad,SB,sbc)

    m_ap[i],msig_ap[i],sb_core[i] = maperture(im,x,y,sky,skysig,magzp)

    print 'm_[3,5,7,10,15,20] = ',m_ap[i]
  #end for i in range(3) - bands


  #fiducial isophotal magnitudes
  r20_k20f = r20[2]
  m20_k20f = pl.zeros([3,2])
  r21_j21f = r21[0]
  m21_j21f = pl.zeros([3,2])
  for i in range(3):
    print 'k fiducial: band %i' %(i)
    #print images2[i],r20_k20f,posang[2][0],ellipticity[2][0],x,y,sky_vals[i][0],sky_vals[i][1],magzps[i]
    m20_k20f[i] = elphot(images2[i],r20_k20f,posang[2][0],ellipticity[2][0],x,y,sky_vals[i][0],sky_vals[i][1],magzps[i])
    print 'm_20 = ',m20_k20f[i]
    m21_j21f[i] = elphot(images2[i],r21_j21f,posang[2][0],ellipticity[2][0],x,y,sky_vals[i][0],sky_vals[i][1],magzps[i])
    print 'm_21 = ',m20_k20f[i]


  #print 'SRC/FIELD: '
  rahms = ac.decimal2hms(ra,'')
  decdms = ac.decimal2dms(dec,'')
  #catfile = open(outcatpath+Jimage_nir+'.nircat1','w')
  catfile = open(outcatpath+Jimage_nir+'.nircat1','w')
    
#write to cat file
  Sline = '%s\t%s\t%s\t' %('ZOA'+rahms+decdms,Jimage_nir.split('_')[0],head.get('DATE_UTC'))
  Sline += '%.5f\t%.5f\t%.5f\t%.5f\t'%(ra,dec,glat,glon)
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(ellipticity[0][0],ellipticity[0][1],posang[0][0],posang[0][1],ellipticity[1][0],ellipticity[1][1],posang[1][0],posang[1][1],ellipticity[2][0],ellipticity[2][1],posang[2][0],posang[2][1])
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(r20[0],m20[0][0],m20[0][1],r20[1],m20[1][0],m20[1][1],r20[2],m20[2][0],m20[2][1])
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(r21[0],m21[0][0],m21[0][1],r21[1],m21[1][0],m21[1][1],r21[2],m21[2][0],m21[2][1])
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(r22[0],m22[0][0],m22[0][1],r22[1],m22[1][0],m22[1][1],r22[2],m22[2][0],m22[2][1])
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(r20_k20f,m20_k20f[0][0],m20_k20f[0][1],m20_k20f[1][0],m20_k20f[1][1],m20_k20f[2][0],m20_k20f[2][1])
  Sline += '%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t' %(r21_j21f,m21_j21f[0][0],m21_j21f[0][1],m21_j21f[1][0],m21_j21f[1][1],m21_j21f[2][0],m21_j21f[2][1])
  for i in range(3):
    Sline += '%.3f\t' %(mT[i])
  for i in range(3):
    Sline += '%.3f\t' %(mT_err[i])
  for i in range(3):
    Sline += '%.3f\t' %(dmT[i])
  for i in range(3):
    for j in range(5):
      Sline += '%.3f\t%.3f\t' %(m_ap[i][j],msig_ap[i][j])
  for i in range(3):
    Sline += '%.3f\t' %(sb_core[i])
  for i in range(3):
    Sline += '%.3f\t' %(peak[i])
  for i in range(3):
    Sline += '%.3f\t' %(sky_marc[i])
  Sline += '%i\t%i\t' %(nx,ny)
  Sline += '%i\t' %(Photflag)
  for i in range(3):
    Sline += '%i\t' %(SBflag[i])
  for i in range(3):
    Sline += '%.3f\t' %(magzps[i])
  for i in range(3):
    Sline += '%.3f\t' %(sees[i])
  for i in range(3):
    Sline += '%.1f\t' %(exposs[i])
  Sline += '%.3f\t' %(Ebv)
  Sline += '%.3f\t' %(SD14)

  catfile.write(Sline+'\n')
  catfile.close()

  return Sline

def hicat(field):

  global nircat1file,NIRLOGfile
  if len(field.split('_'))==1:
    ind = afind(field,HICAT[0])
  else:
    C = field.split('_')
    inds = afind(C[0],HICAT[0])
    ind = [inds[int(C[1])-1]]
  for i in range(len(ind)):
    # wcs coords
    Jfield_ind = HICAT[2][ind[i]]
    rah = HICAT[3][ind[i]]
    dech = HICAT[4][ind[i]]
    # initial a,b,phi - in image coords
    ah = HICAT[5][ind[i]]
    bh = HICAT[6][ind[i]]
    pah = HICAT[7][ind[i]]+90.  # hicat stores ccw from image x/ellipse requires ccw from image y

    #convert a,b to pix
    if bh > ah:
      tab = ah
      ah = bh
      bh = tab
      pah = pah + 90.
    ah = ah*3600/ps
    bh = bh*3600/ps
    # normalise pa to (-90,90]
    if pah > 270: pah = -(360-pah)
    elif pah > 180: pah = pah - 180
    elif pah > 90: pah = -(180-pah)

    nircat1 = open(nircat1file,'a')
    NIRLOG = open(NIRLOGfile,'a')
    try:
    #if True:
      photline = photloop(Jfield_ind,rah,dech,ah,bh,pah)
      print photline
      #dummy = raw_input('press enter')
      nircat1.write(photline+'\n')
      NIRLOG.write(Jfield_ind+'\tdone'+'\n')
    except Exception as e:
    #else:
      print e
      print 'Photometry failed for %s' %(Jfield_ind)
      NIRLOG.write(Jfield_ind+'\tFAIL'+'\n')
    nircat1.close()
    NIRLOG.close()
  pl.close('all')
  return

########################################################################
##########################   MAIN   ####################################
########################################################################


#ds9 = pysao.ds9()
Csystem = "galactic"

#flist="galaxyfields.list.1"
flist = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated_list'
#flist = '2masscomplist'
#flist = 'calibrated_list_fail'
#flist = '2masscal_list'
#flist = 'nswcal_list'
if os.path.isfile(flist):
  infile = open(flist,'r')
  fieldlist = []
  for line in infile:
    if line[0] == '#': continue
    fieldlist += [line.strip()]
  infile.close()
else:
  print 'no fieldlist: %s' %(flist)

lastfield = 'lastfield'
if os.path.isfile(lastfield):
  lastN = open(lastfield,'r')
  n = int(lastN.readline().strip())
  lastN.close()
else:
  n = 0


#fieldlist=['J0712-09','J0716-18C','J0727-23','J0730-25','J0730-28','J0739-24','J0740-22','J0740-30A','J0741-22','J0744-13','J0746-18A_1']'J0740-32_1'['J0740-32']
#fieldlist = ['J0712-09']

#fieldlist = ['J0740-22']\
#fieldlist=['J0748-26B']

nircat1file = 'nircat1'
NIRLOGfile = 'nirphot.log'
if os.path.isfile(nircat1file): nircat1 = open(nircat1file,'a')
else: nircat1 = open(nircat1file,'w')
nircat1.close()
if os.path.isfile(NIRLOGfile): NIRLOG = open(NIRLOGfile,'a')
else: NIRLOG = open(NIRLOGfile,'w')
NIRLOG.close()

for field in fieldlist:
  print field
  hicat(field)
  #cont = raw_input('q to quit')
  #if cont == 'q': break
