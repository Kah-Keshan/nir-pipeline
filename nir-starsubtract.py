import os
wd=os.getcwd()
os.chdir('/Volumes/BigDalek/khaled/')
from pyraf import iraf
os.chdir(wd)
#iraf.sirius()
iraf.stsdas()
iraf.analysis()
iraf.isophote()
iraf.unlearn(iraf.daophot)
from datetime import datetime as d
import ds9 as pyds9
#from ds9 import*
import pyfits as pf
import subprocess as sub
import asciidata as ad
import pylab as pl
import math
import time
import pywcs
import aplpy
from aplpy import make_rgb_image as RGB
from ds9disp import dodisp, panto, wcs_mark, wcs_box
from skyvgrad import remove_vgrad
from iraf import daophot
from astropy.wcs import WCS
import numpy as np
#run ds9 display
Display = True
AUTO = False
#AUTO = True
#print
#raw_input('running AUTO ',AUTO,': )
if AUTO: Display = False

# SEx paramters
SEconfig = './SEconfig'
SEparam = './SE.param'
SEnnw = './default.nnw'

#other variables
#tmppath= "/dev/shm/"
tmppath= './tmp/'

class PsfError(Exception):
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return `self.value`



##################################################################

def mkpsf(corename,outname, verbose):
    '''
        mkpsf - make psf for image corename.fits
        
        usage: mkpsf(corename,verbose)
        corename (string): name of fits file without the fits extension (can include full path)
        verbose (string): 'yes' or 'no' controls output of iraf DAOPHOT functions
        
        '''
    #data_max = 30000
    # Parameters for DAOPHOT
    psfrad=12
    #psf_data_max=10000
    
    loop = 0
    
    os.system('rm %stmp*' %(tmppath))
    
    imgname= corename+".fits"
    psfname= corename+".psf.fits"
    tmplog = tmppath+"tmp.log"
    print "### making psf for %s" %(imgname)
    
    iraf.datapar.readnoise=147
    iraf.datapar.epadu=132
    iraf.datapar.datamax=30000
    
    iraf.fitskyp.salgori="mode"
    iraf.fitskyp.annulus=12
    iraf.fitskyp.dannulu=3
    iraf.fitskyp.sloclip=0
    iraf.fitskyp.shiclip=0
    
    iraf.daopars.function="auto"
    iraf.daopars.varorder=1  # change from 1
    
    iraf.daofind.interac='no'
    iraf.daofind.verify='no'
    iraf.daofind.verbose= verbose
    iraf.phot.interac='no'
    iraf.phot.verify='no'
    iraf.phot.verbose=verbose
    iraf.allstar.verify='no'
    iraf.allstar.verbose=verbose
    iraf.substar.verify='no'
    iraf.substar.verbose=verbose
    
    # measure sky and rms from image
    skylevel = float(iraf.imstat(images=imgname, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    stdev = float(iraf.imstat(images=imgname, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    iraf.datapars.sigma=stdev
    iraf.datapars.datamin = skylevel - 5*iraf.datapars.sigma
    
    head = pf.getheader(imgname)
    seeing = head.get('SEEING')
    if seeing is None:
        raise PsfError('seeing')
        return
    
    fwhm   = seeing / 0.45
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
    
    print "skylevel: ",skylevel
    print "stddev: ",iraf.datapars.sigma
    print "data_min: ", iraf.datapars.datamin
    print "seeing: ", seeing
    print "fwhm: ", fwhm
    print "magxp: ", magzp
    print "expos: ", expos
    
    iraf.daofind(image=imgname, output=tmppath+"tmp1.coo.1", threshold = 3, verify='no')
    iraf.pselect(infiles=tmppath+"tmp1.coo.1", outfiles=tmppath+"tmp1.coo.2", expr="SHARPNESS != INDEF")
    
    iraf.photpars.aperture=iraf.datapars.fwhmpsf
    iraf.daopars.fitrad=iraf.datapars.fwhmpsf
    if(iraf.datapars.fwhmpsf < 3):
        iraf.photpars.aperture=3
        iraf.daopars.fitrad=3
    
    iraf.phot(image=imgname, coords=tmppath+"tmp1.coo.2", output=tmppath+"tmp1.mag.1", verify='no')
    iraf.pselect(tmppath+"tmp1.mag.1", tmppath+"tmp1.sel.1", expr="MAG>11 && MAG<14 && MERR>=0 && MERR<0.05 && STDEV<"+str(iraf.datapar.sigma))
    
    iraf.pstselect(image=imgname, photfile=tmppath+"tmp1.sel.1", pstfile=tmppath+"tmp1.pst.1", maxnpsf=50, datamax=15000, interac='no', verify='no', verbose=verbose)
    
    os.system("awk 'NF==5 {print $2,$3}' %stmp1.pst.1 > %stmp1.pst.1.dump" %(tmppath,tmppath))
    
    print "# psf fit 1:"
    iraf.psf(image=imgname, photfile=tmppath+"tmp1.sel.1", pstfile=tmppath+"tmp1.pst.1", psfimage=tmppath+"tmp1.psf.1", opstfile=tmppath+"tmp1.pst.2", groupfil=tmppath+"tmp1.psg.1",datamax=15000, psfrad=psfrad, interac='no',verify='no', Stdout=tmppath+"tmp1.psf.log.1")
    os.system('tail -20 %stmp1.psf.log.1' %(tmppath))
    
    print "# allstar 1:"
    iraf.allstar(image=imgname, photfile=tmppath+"tmp1.mag.1", psfimage=tmppath+"tmp1.psf.1",allstarfile=tmppath+"tmp1.als.1", rejfile=tmppath+"tmp1.arj.1", subimage=tmppath+"tmp1.sub.1",psfrad=psfrad )
    
    iraf.daofind(image=tmppath+"tmp1.sub.1", output=tmppath+"tmp1.coo.5", threshold = 2.25)
    iraf.pselect(infiles=tmppath+"tmp1.coo.5", outfiles=tmppath+"tmp1.coo.6", expr="SHARPNESS != INDEF")
    iraf.phot(image=tmppath+"tmp1.sub.1", coords=tmppath+"tmp1.coo.6", output=tmppath+"tmp1.mag.2")
    
    print "# allstar 2:"
    iraf.allstar(image=tmppath+"tmp1.sub.1", photfile=tmppath+"tmp1.mag.2", psfimage=tmppath+"tmp1.psf.1",allstarfile=tmppath+"tmp1.als.2", rejfile=tmppath+"tmp1.arj.2", subimage=tmppath+"tmp1.sub.2",psfrad=psfrad )
    iraf.pconcat(infiles=tmppath+"tmp1.als.1,"+tmppath+"tmp1.als.2", outfile=tmppath+"tmp1.als.3")
    
    iraf.pdump(tmppath+"tmp1.als.1", "XCENTER,YCENTER,MAG,MERR", 'yes', Stdout=tmppath+"tmp1.als.1.pdump")
    os.system("awk '$3>11 && $3<14 && $4>0 && $4<0.03' %stmp1.als.1.pdump > %stmp1.als.2.pdump" %(tmppath,tmppath))
    iraf.pselect(tmppath+"tmp1.als.1", tmppath+"tmp1.sel.2", expr="MAG>11 && MAG<14 && MERR>=0 && MERR<0.02")
    os.system("iso_star_select %stmp1.sel.2 30 > %stmp1.sel.3" %(tmppath,tmppath))
    
    print "# substar 1:"
    iraf.substar(image=imgname, photfile=tmppath+"tmp1.als.3", exfile=tmppath+"tmp1.sel.3",psfimage=tmppath+"tmp1.psf.1", subimage=tmppath+"tmp2.sub.1", psfrad=psfrad)
    
    iraf.pdump(tmppath+"tmp1.sel.3", "XCENTER,YCENTER,MAG,MERR", 'yes',Stdout=tmppath+"tmp1.sel.3.pdump")
    iraf.phot(image=tmppath+"tmp2.sub.1", coords=tmppath+"tmp1.sel.3.pdump", output=tmppath+"tmp2.mag.1")
    iraf.pselect(tmppath+"tmp2.mag.1", tmppath+"tmp2.sel.1",expr="MAG>11 && MAG<14 && MERR>=0 && MERR<0.05 && STDEV<"+str(stdev))
    
    iraf.pdump(tmppath+"tmp2.sel.1", "XCENTER,YCENTER,MAG,MERR", 'yes',Stdout=tmppath+"tmp2.sel.1.pdump")
    iraf.pstselect(image=tmppath+"tmp2.sub.1", photfile=tmppath+"tmp2.sel.1", pstfile=tmppath+"tmp2.pst.1", maxnpsf=50, datamax=15000, interac='no', verify='no', verbose=verbose)
    
    print "# psf fit 2:"
    iraf.psf(image=tmppath+"tmp2.sub.1", photfile=tmppath+"tmp2.sel.1", pstfile=tmppath+"tmp2.pst.1",psfimage=tmppath+"tmp2.psf.1", opstfile=tmppath+"tmp2.pst.2", groupfil=tmppath+"tmp2.psg.1",datamax=15000, psfrad=psfrad, interac='no', verify='no',Stdout=tmppath+"tmp2.psf.log.2")
    os.system("tail -20 %stmp2.psf.log.2" %(tmppath))
    
    os.system("awk 'NF==5 {print $2,$3}' %stmp2.pst.2 > %stmp2.pst.1.dump" %(tmppath,tmppath))
    
    print "# allstar 3:"
    iraf.allstar(image=imgname, photfile=tmppath+"tmp1.als.3", psfimage=tmppath+"tmp2.psf.1", allstarfile=tmppath+"tmp2.als.1", rejfile=tmppath+"tmp2.arj.1", subimage=tmppath+"tmp2.sub.2", psfrad=psfrad )
    
    print "# substar 2:"
    iraf.substar(image=imgname, photfile=tmppath+"tmp1.als.3", exfile=tmppath+"tmp1.sel.3", psfimage=tmppath+"tmp2.psf.1", subimage=tmppath+"tmp3.sub.1", psfrad=psfrad)
    
    print "psf fit 3:"
    iraf.psf(image=tmppath+"tmp3.sub.1", photfile=tmppath+"tmp2.sel.1", pstfile=tmppath+"tmp2.pst.1", psfimage=tmppath+"tmp3.psf.1", opstfile=tmppath+"tmp3.pst.1", groupfil=tmppath+"tmp3.psg.1",datamax=15000, psfrad=psfrad, interac='no',verify='no',Stdout=tmppath+"tmp3.psf.log.1")
    os.system("tail -20 %stmp3.psf.log.1" %(tmppath))
    
    print "# allstar 4:"
    iraf.allstar(image=imgname, photfile=tmppath+"tmp2.als.1", psfimage=tmppath+"tmp3.psf.1",allstarfile=tmppath+"tmp3.als.1", rejfile=tmppath+"tmp3.arj.1", subimage=tmppath+"tmp3.sub.2",psfrad=psfrad )
    
    print "copying files"
    os.system("cp %s %s.psf.fits" %(tmppath+'tmp3.psf.1.fits',outname))
    os.system("cp %s %s.pst" %(tmppath+'tmp3.pst.1',outname))
    os.system("cp %s %s.psg" %(tmppath+'tmp3.psg.1',outname))
    
    print 'psf file %s.psf.fits created' %(outname)
    
    return

# turn coo file to reg file for easy display
def coo2reg(coo,reg,col='white',r = 3.,wcs='image'):
    '''
        coo2reg - converts DAOFIND coo file to ds9 region file
        (internal function of psfpart)
        
        usage: coo2reg(coo,reg,col='white',r = 3.,wcs='image')
        coo (str) - .coo file
        reg (str) - ds9 region file
        col (str) - ds9 region colors
        r (float) - ds9 region size
        wcs (str) - coords in wcs 'image'/'fk5'
        '''
    coof = open(coo,'r')
    regf = open(reg,'w')
    regf.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    regf.write('physical \n')
    for line in coof:
        if line[0]!='#':
            C = line.split()
            x = float(C[0])
            y = float(C[1])
            if wcs == 'image':
                Sline = 'image; circle(%f,%f,%f) # color=%s\n' %(x,y,r,col)
            Sline = 'physical; circle(%f,%f,%f) # color=%s\n' %(x,y,r,col)
            regf.write(Sline)
    coof.close()
    regf.close()
##################################################################


### use Source Extractor to find stars in images ###
def findstars(im,x0,y0,thresh,coofile,lim,sky):
    '''
        findstars - use Source Extractor to find stars in images
        (internal function of psfpart)
        
        usage: findstars(im,x0,y0,thresh,coofile,lim,sky)
        '''
    #print im,x0,y0,thresh,coofile,lim
    im1 = im +'.fits'
    fits=pf.open(im1)
    magzp=fits[0].header['MAGZP']
    seeing=fits[0].header['SEEING']
    xm = fits[0].header['NAXIS1']/2
    ym = fits[0].header['NAXIS1']/2
    imcut = '%s[%i:%i,%i:%i]' %(im,lim[0],lim[1],lim[2],lim[3])
    im2 = im+'_ps.fits'
    iraf.imutil.imcopy(imcut,im2)
    # check for required files for SE to run and run SE on Hband
    if not os.path.isfile(SEparam) and os.path.isfile(SEnnw) and os.path.isfile(SEconfig):
        print "missing required SE config files"
        return

    sub.Popen(["sex",im2,"-c",SEconfig,"-CATALOG_NAME","image.cat","-PARAMETERS_NAME",SEparam,"-MAG_ZEROPOINT",str(magzp),'-SEEING_FWHM',str(seeing),"-STARNNW_NAME",SEnnw,"-DETECT_THRESH",str(thresh),'-BACK_VALUE',str(sky)]).wait()
    #os.system('rm %s' %(im2))
    
    # open ds9 region file
    reg = open('stars.reg','w')
    coo = open(coofile,'w')
    reg.write('global color=green font="helvetica 10 normal" select=1 highlite=1 edit=1 move=1 delete=1 include=1 fixed=0 source\n')
    reg.write('physical\n')
    print 'FINDSTARS threshold %.0f: (red - nonstar, blue - star, green - faint star?)' %(thresh)



    if os.path.isfile('image.cat'):
    

        SEcat = ad.open('image.cat')



        #calculate mean mu



	cen_sb = SEcat['MU_MAX'].tonumpy()
	mean_sb = SEcat['MAG_ISOCOR'].tonumpy()  + 2.5*pl.log10(0.45*pl.pi*SEcat['A_IMAGE'].tonumpy()*SEcat['B_IMAGE'].tonumpy())


        # select 'stars' and write to starlist and mark SE detections



        for i in range(SEcat.nrows):
        
            x = SEcat['X_IMAGE'][i]
            y = SEcat['Y_IMAGE'][i]
            a = SEcat['A_IMAGE'][i]
            b = SEcat['B_IMAGE'][i]
            if SEcat['MAG_ISOCOR'][i] < 90.:
                star = True
                col = 'blue'
                if (b/a < 0.5) and (a > 5.):
                    col = 'red'
                    star = False
                #dont find the core
                if pl.sqrt(((x+lim[0]-1-x0))**2 + (y+lim[2]-1-y0)**2) < 5.:
                    col = 'white'
                    star = False
                Sline = 'image; ellipse(%f,%f,%f,%f,%f) # color=%s\n' %(x+lim[0]-1,y+lim[2]-1,a,b,-SEcat['THETA_IMAGE'][i],col)
                reg.write(Sline)
                if star:
                    Sline = "%f\t%f\n" %(x+lim[0]-1,y+lim[2]-1)
                    coo.write(Sline)
        reg.close()
        coo.close()
        #load region file in all frames



        if Display: ds9.set('regions load stars.reg')
        ds9.set('regions load stars.reg')
 
    return

def cleanstarlist(coofile,xr=pl.array([]),yr=pl.array([])):
    '''
        cleanstarlist - interactively remove misidentified stars on galaxy
        (internal function of psfpart)
        '''
    menu = ''' # clean stars #
        (r) remove non-stars
        (a) add non-detected stars
        (q) done
        '''
    cleanstill = True
    #xr = []
    #yr = []
    if len(xr) > 0:
        for i in range(len(xr)):
            wcs_mark(ds9,xr[i],yr[i],colour='red')
    xa = pl.array([])
    ya = pl.array([])
    while cleanstill:
        clean = raw_input(menu)
        if clean == 'r':
            # deselect stars
            x0 = 9999
            y0 = 9999
            print 'click on source that is not a star: (click on same source twice to quit) '
            while True:
                crd = ds9.get('imexam coordinate ' +'image')
                C = crd.split()
                x = float(C[0])
                y = float(C[1])
                if pl.sqrt((x-x0)**2+(y-y0)**2) < 2:
                    print 'done'
                    break
                xr=pl.append(xr,x)
                yr=pl.append(yr,y)
                wcs_mark(ds9,x,y,colour='red')
                x0 = x
                y0 = y
            print "%i non-stars selected for removal from coo file" %(len(xr))
        
        # add stars
        elif clean == 'a':
            x0 = 9999
            y0 = 9999
            
            print 'click on source that is a star: (click on same source twice to quit) '
            while True:
                crd = ds9.get('imexam coordinate ' +'image')
                C = crd.split()
                x = float(C[0])
                y = float(C[1])
                if pl.sqrt((x-x0)**2+(y-y0)**2) < 2:
                    print 'done'
                    break
                xa=pl.append(xa,x)
                ya=pl.append(ya,y)
                wcs_mark(ds9,x,y,colour='blue')
                x0 = x
                y0 = y
            print "%i stars selected for addition to coo file" %(len(xa))
        
        # done cleaning
        else: #clean =='q':
            
            # read coo file and write new temp file
            if len(xr)>0:
                coo = open(coofile,'r')
                cooout = open(tmppath+'newcoo.tmp','w')
                i = 0
                print '%i stars to be removed' %(len(xr))
                for line in coo:
                    if line[0] == '#': cooout.write(line)
                    else:
                        C = line.split()
                        x = float(C[0])
                        y = float(C[1])
                        sep = pl.sqrt((x-xr)**2. + (y-yr)**2)
                        if not any(sep < 4): cooout.write(line)
                        else: i += 1
                print "%i non-stars removed from coo file" %(i)
                coo.close()
            else:
                os.system('cp -v %s %s' %(coofile,tmppath+'newcoo.tmp'))
                cooout = open(tmppath+'newcoo.tmp','a')
            
            # add coords to coofile
            for i in range(len(xa)):
                line = '%.5f\t%.5f\n' %(xa[i],ya[i])
                cooout.write(line)
            cooout.close()
            print "%i stars added to coo file" %(len(xa))
            
            #copy temp coo to coo
            if os.path.isfile(tmppath+'newcoo.tmp'):
                os.system('cp -v %snewcoo.tmp %s' %(tmppath,coofile))
            print 'done cleaning'
            cleanstill = False
    
    return xr,yr

# routine to match daofind and se lists
def matchdaose(daoin,sein,out,first,ellipse):
    '''
        matchdaose - routine to match daofind and se lists
        (internal function of psfpart)
        '''
    #print ellipse
    rad = 2. ##pix
    daof = open(daoin,'r')
    xd = []
    yd = []
    for line in daof:
        if line[0] != '#':
            C = line.split()
            xd += [float(C[0])]
            yd += [float(C[1])]
    xd = pl.array(xd)
    yd = pl.array(yd)
    daof.close()
    firstf = open(first,'r')
    xf = []
    yf = []
    for line in firstf:
        if line[0] != '#':
            C = line.split()
            xf += [float(C[0])]
            yf += [float(C[1])]
    xf = pl.array(xf)
    yf = pl.array(yf)
    firstf.close()
    sef = open(sein,'r')
    xs = []
    ys = []
    for line in sef:
        if line[0] != '#':
            C = line.split()
            xs += [float(C[0])]
            ys += [float(C[1])]
    xs = pl.array(xs)
    ys = pl.array(ys)
    sef.close()
    outf = open(out,'w')
    x0_e, y0_e, a_e, eps_e, phi_e = ellipse
    #ellipse_phi is deg cntclock from +y - ellipse eqn needs rad cntclock from +x
    phi_e = (90-phi_e)*pl.pi/180.
    #phi_e = (phi_e)*pl.pi/180.
    
    for i in range(len(xs)):
        r2 = pl.sqrt((xs[i]-xf)**2.+(ys[i]-yf)**2.) #first - already removed
        if not any(r2 < rad):
            #print i,xs[i],ys[i],lim, (xs[i] > lim[0]) ,(xs[i] < lim[1]) ,(ys[i] > lim[2]), (ys[i] < lim[3])
            # if in the initial ellipse, must match daof
            ai=pl.sqrt(((xs[i]-x0_e)*pl.cos(phi_e)-(ys[i]-y0_e)*pl.sin(phi_e))**2.+(((xs[i]-x0_e)*pl.sin(phi_e)+(ys[i]-y0_e)*pl.cos(phi_e))/eps_e)**2.)
            if ai < a_e:
                r = pl.sqrt((xs[i]-xd)**2.+(ys[i]-yd)**2.) #daofind
                if any(r < rad):
                    outf.write('%f\t%f\n'%(xs[i],ys[i]))
            else:
                #print out
                outf.write('%f\t%f\n'%(xs[i],ys[i]))
    outf.close()

def psfpart(corename,outname,psfname,galx,galy,size,verbose,def_psf,maxsma,nloop,a,b,pa,semiunit='deg'):
    '''
        psfpart - this routine does starsubtraction on single image
        # requires psf, galaxy coordinates (pix)
        # galaxy initial a,b , pa (world)
        
        usage: psfpart(corename,outname,psfname,galx,galy,size,verbose,def_psf,maxsma,nloop,a,b,pa)
        corename - input image name (exclude .fits) (may include full path if not in cwd)
        outname (str) - output image name (include .fits)
        psfname (str) - psf fits name
        galx,galy (float) - initial x,y coords of galaxy (pixels)
        size (int) - size in pixels of box within which to remove stars
        verbose ('yes','no')- controls amount of output for DAOPHOT
        def_psf ('yes','no')- use default psf (yes for corename.psf.fits, no for specified psf)
        maxsma (float) - maximum radius of galaxy (for ellipse fitting)
        nloop (int) - number of iterations
        a,b (float) - rough semi-minor -major axes in semiunit
        pa (float) - rough position angle of galaxy
        semiunit (str='deg'/'pix') - set to pix for a/b in pixels
        
        '''
    psfrad = 12.
    loop = 1
    os.system('rm %s*' %(tmppath))
    os.system('rm %s_*' %(tmppath))
    
    imgname = corename+".fits"
    print
    print '*** Running star-subtraction on image %s ***' %(imgname)
    print " - x,y = %f %f" %(galx,galy)
    
    if(def_psf=='yes'):
        psfname = corename+".psf.fits"
    tmplog = tmppath+"tmp.log"
    
    #convert a,b to pix
    if b > a:
        tab = a
        a = b
        b = tab
        pa = pa + 90.
    if semiunit=='deg':
        a = a*3600/0.45
        b = b*3600/0.45
    ell = 1-b/a
    eps = b/a
    ell = max(0.05,ell)
    # normalise pa to (-90,90]
    pa = pa+90.
    # normalise pa to (-90,90]
    if pa > 270: pa = -(360-pa)
    elif pa > 180: pa = pa - 180
    elif pa > 90: pa = -(180-pa)
    
    print a,b,ell,eps,pa
    # image size
    
    head = pf.getheader(imgname)
    
    ximgsize = head.get('NAXIS1')
    yimgsize = head.get('NAXIS2')
    
    # subimage limits: check if runs over edges
    xlim1 =  galx - (size/2)
    if(xlim1<1):
        xlim1=1
    xlim2 =  galx + (size/2)
    if(xlim2>ximgsize):
        xlim2=ximgsize
    ylim1 =  galy - (size/2)
    if(ylim1<1):
        ylim1=1
    ylim2 =  galy + (size/2)
    if(ylim2>yimgsize):
        ylim2=yimgsize
    
    xlim1 = int(xlim1)
    ylim1 = int(ylim1)
    xlim2 = int(xlim2)
    ylim2 = int(ylim2)
    
    # 4x4 pixels around galaxy core
    xlim3 = int(galx-2)
    xlim4 = int(galx+2)
    ylim3 = int(galy-2)
    ylim4 = int(galy+2)
    
    # get skylevel and rms
    skylevel = float(iraf.imstat(images=imgname, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    stdev = float(iraf.imstat(images=imgname, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    iraf.datapars.sigma = stdev
    iraf.datapars.datamin = skylevel - 5*iraf.datapars.sigma
    print ' - skylevel is %f' %(skylevel)
    print ' - sky rms is %f' %(stdev)
    
    seeing = head.get('SEEING')
    iraf.datapars.fwhmpsf = seeing / 0.45
    print ' - seeing is %f arcsec' %(seeing)
    
    iraf.centerpars.calgorithm = 'centroid'
    
    g_imgname = imgname
    
    #subtract sky
    iraf.imarith(imgname, "-", skylevel, tmppath+"_orgskysub")
    # iterate
    subloop = True
    disp = False
    frm = 0
    while subloop:
        
        

        if Display: ds9.set('frame delete all')
        if loop >= 3 and Display:
            disp = True
        print '######################'
        print 'ITERATION %i' %(loop)
        loop_name =tmppath+"_tmp.gsub"+str(loop)
        
        iraf.imarith(g_imgname, "-", skylevel, tmppath+"tmp_skysub")
        if disp:
            ds9.set('frame 1')
            frm = int(ds9.get('frame'))
            ds9.set('file %s.fits' %(tmppath+"tmp_skysub"))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale zscale')
            ds9.set('scale linear')
        print '# FRAME %i: loop image for modelling %i' %(frm,loop)
        
        # model galaxy
        
        iraf.ellipse(tmppath+"tmp_skysub", output=tmppath+"tmp_tab",x0=galx, y0=galy, ellip0=ell, pa0=pa, maxsma=maxsma, hcenter='yes', hellip='no', hpa='no', recenter='yes', olthres=0.1, nclip=10, step=0.02,verbose='no',interactive='no')
        print "test_1"
	if disp: wcs_mark(ds9,galx,galy,a,b,pa+90.)
        print "test_2"
	iraf.bmodel(table=tmppath+"tmp_tab",output=tmppath+"tmp_model1", parent="")
	print "test_3"
        iraf.boxcar(tmppath+"tmp_model1",output=tmppath+"tmp_model",xwindow=3,ywindow=3)
       	print "test_4"
	if disp:
            ds9.set('frame 2')
            frm = int(ds9.get('frame'))
            ds9.set('file %stmp_model.fits' %(tmppath))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale minmax')
            ds9.set('scale linear')
        print '# FRAME %i: model galaxy: %i' %(frm,loop)
        
        #subtract model galaxy
        iraf.imarith(tmppath+"_orgskysub","-",tmppath+"tmp_model",tmppath+"tmp_gsub")
        iraf.imarith(tmppath+"tmp_gsub", "+", skylevel, tmppath+"tmp_gsub2")
        
        iraf.imreplace(tmppath+"tmp_gsub2["+str(xlim3)+":"+str(xlim4)+","+str(ylim3)+":"+str(ylim4)+"]", skylevel, lower=iraf.INDEF, upper=iraf.INDEF)
        
        if disp:
            ds9.set('frame 3')
            frm = int(ds9.get('frame'))
            ds9.set('file %stmp_gsub2.fits' %(tmppath))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale zscale')
            ds9.set('scale linear')
        print '# FRAME %i: model subtracted, pixels replaced: %i' %(frm,loop)
        
        #run SE



        findstars(tmppath+"tmp_gsub2",galx,galy,3.5,tmppath+"tmp1.coo.2",[xlim1,xlim2,ylim1,ylim2],skylevel)


        # remove non-stars on gal
        if not AUTO:
            if loop > nloop - 2:
                xrej,yrej=cleanstarlist(tmppath+"tmp1.coo.2",[],[])
        
        iraf.photpars.aperture=iraf.datapars.fwhmpsf
        iraf.daopars.fitrad=iraf.datapars.fwhmpsf * 2
        if(iraf.datapars.fwhmpsf < 3):
            iraf.photpars.aperture=3
            iraf.daopars.fitrad=3
        
        iraf.phot(image=tmppath+"tmp_gsub2", coords=tmppath+"tmp1.coo.2", output=tmppath+"tmp1.mag.1",verbose='no',interactive='no',verify='no')
        
        iraf.allstar(image=tmppath+"tmp_gsub2", photfile=tmppath+"tmp1.mag.1", psfimage=psfname,allstarfile=tmppath+"tmp1.als.1", subimage=tmppath+"tmp1.sub.1", rejfile=tmppath+"tmp1.arj.1",critsnratio=0.1,psfrad=psfrad,verbose='no',verify='no')
        if os.path.isfile(tmppath+'allstar1'):
            os.system('rm '+tmppath+'allstar1')
        if os.path.isfile(tmppath+'allstar_rej1'):
            os.system('rm '+tmppath+'allstar_rej1')
        iraf.pdump(tmppath+"tmp1.als.1",'XCEN,YCEN','yes',Stdout=tmppath+'allstar1')
        iraf.pdump(tmppath+"tmp1.arj.1",'XCEN,YCEN','yes',Stdout=tmppath+'allstar_rej1')
        
        if disp:
            coo2reg(tmppath+'allstar1',tmppath+'allstar1.reg','black',r=1)
            coo2reg(tmppath+'allstar_rej1',tmppath+'allstar_rej1.reg','cyan',r=1)
            ds9.set('frame 4')
            frm = int(ds9.get('frame'))
            ds9.set('file %stmp1.sub.1.fits' %(tmppath))
            ds9.set('regions load %s' %(tmppath+'allstar1.reg'))
            ds9.set('regions load %s' %(tmppath+'allstar_rej1.reg'))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale zscale')
            ds9.set('scale linear')
        print '# FRAME %i: allstar (black - good, cyan - reject): %i' %(frm,loop)
        
        iraf.imreplace(tmppath+"tmp1.sub.1", skylevel, lower='INDEF', upper=iraf.datapars.datamin)
        
        iraf.daofind(image=tmppath+"tmp1.sub.1", output=tmppath+"tmp1.coo.3", threshold = 1.8, verify='no',verbose='no')
        
        iraf.pselect(infiles=tmppath+"tmp1.coo.3", outfiles=tmppath+"tmp1.coo.4",expr="XCENTER >"+str(galx-a)+" && XCENTER <"+str(galx+a)+"&& YCENTER >"+str(galy-a)+" && YCENTER <"+str(galy+a)  )
        
        print '# FRAME %i: daofind threshold 2 (magenta): %i' %(frm,loop)
        # SE to find stars


        findstars(tmppath+"tmp1.sub.1",galx,galy,1.5,tmppath+"tmp1.coo.5",[xlim1,xlim2,ylim1,ylim2],skylevel)


        
        # match daof and se sources and clean any detectins of previous residuals
        matchdaose(tmppath+"tmp1.coo.4",tmppath+"tmp1.coo.5",tmppath+"tmp1.coo.6",tmppath+"tmp1.coo.2",[galx,galy,1.25*a,eps,pa])
	

#        tmpLIST6 = tmppath + "tmp1.coo.4"
        
        if disp:
            coo2reg(tmppath+"tmp1.coo.6",tmppath+"tmp_reg",'green',r=4)
            # remove stars on gal
            #wcs_box(ds9,galx,galy,2*a,2*a)
            wcs_mark(ds9,galx,galy,a,b,pa+90.)
            ds9.set('regions load %stmp_reg\n' %(tmppath))
        
        # remove non-stars on gal
        if not AUTO:
            if loop > nloop -2:
                xrej2,yrej2=cleanstarlist(tmppath+"tmp1.coo.6",xrej,yrej)
        
        iraf.phot(image=tmppath+"tmp1.sub.1", coords=tmppath+"tmp1.coo.6", output=tmppath+"tmp1.mag.2",verbose='no',interactive='no',verify='no')
# 	iraf.phot(image=tmppath+"tmp1.sub.1", coords=tmpLIST6, output=tmppath+"tmp1.mag.2",verbose='no',interactive='no',verify='no')

        
        iraf.allstar(image=tmppath+"tmp1.sub.1", photfile=tmppath+"tmp1.mag.2", psfimage=psfname,allstarfile=tmppath+"tmp1.als.2", subimage=tmppath+"tmp1.sub.2", rejfile=tmppath+"tmp1.arj.2",psfrad=psfrad ,verbose='no')
        if os.path.isfile(tmppath+'allstar2'):
            os.system('rm '+tmppath+'allstar2')
        if os.path.isfile(tmppath+'allstar_rej2'):
            os.system('rm '+tmppath+'allstar_rej2')
        iraf.pdump(tmppath+"tmp1.als.2",'XCEN,YCEN','yes',Stdout=tmppath+'allstar2')
        iraf.pdump(tmppath+"tmp1.arj.2",'XCEN,YCEN','yes',Stdout=tmppath+'allstar_rej2')
        if disp:
            coo2reg(tmppath+'allstar2',tmppath+'allstar2.reg','black',r=1)
            coo2reg(tmppath+'allstar_rej2',tmppath+'allstar_rej2.reg','yellow',r=1)
            ds9.set('frame 5')
            frm = int(ds9.get('frame'))
            ds9.set('file %stmp1.sub.2.fits' %(tmppath))
            ds9.set('regions load %s' %(tmppath+'allstar2.reg'))
            ds9.set('regions load %s' %(tmppath+'allstar_rej2.reg'))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale zscale')
            ds9.set('scale linear')
        print '# FRAME %i: allstar (black - good, yellow - reject): %i' %(frm,loop)
        
        iraf.pconcat(infiles=tmppath+"tmp1.als.1,"+tmppath+"tmp1.als.2", outfile=tmppath+"tmp1.als.3")
        
        #print "substar"
        iraf.substar(image=imgname, photfile=tmppath+"tmp1.als.3", psfimage=psfname,subimage=loop_name, psfrad=psfrad,exfile="",verbose='no',verify='no')
        
        #automatically clean residuals outside of masked galaxy
        # only clean the bright stars removed
        if loop > 2:
            if os.path.isfile(tmppath+'starlist'): os.system('rm %sstarlist' %(tmppath))
            dat = open(tmppath+'allstar1','r')
            xrem = []
            yrem = []
            for line in dat:
                C= line.split()
                xrem += [float(C[0])]
                yrem += [float(C[1])]
            dat.close()
            autocleanresiduals(loop_name,xrem,yrem,[galx,galy,1.5*a,eps,pa])#[galx-1.*a,galx+1.*a,galy-1.*a,galy+1.*a])
        
        if disp:
            ds9.set('frame 6')
            frm = int(ds9.get('frame'))
            ds9.set('file %s.fits' %(loop_name))
            ds9.set('pan to %i %i image' %(galx,galy))
            ds9.set('zoom to 2')
            ds9.set('scale zscale')
            ds9.set('scale linear')
        print '# FRAME %i: substar done (new g_img): %i' %(frm,loop)
        
        if loop < nloop -1:
            print 'LOOP %i complete' %(loop)
            dummy = 'c'
        else:
            if not AUTO: dummy = raw_input('LOOP %i complete, will finish now unless you loop again (l)' %(loop))
            else: dummy = 'q'
        if dummy == 'q': subloop = False
        loop = loop+1
        if loop > nloop -1 and dummy != 'l': subloop = False
        g_imgname = loop_name
        os.system("rm %stmp*" %(tmppath))
        
        
        skylevel = float(iraf.imstat(images=loop_name+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xlim1,xlim2,ylim1,ylim2), fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        stdev = float(iraf.imstat(images=imgname, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        iraf.datapars.sigma = stdev
        iraf.datapars.datamin = skylevel - 5*iraf.datapars.sigma
        print 'skylevel = ',skylevel
    
    print "loop: "+loop_name
    print "out: "+outname
    
    inp = loop_name+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xlim1,xlim2,ylim1,ylim2)
    if os.path.isfile(outname):
        os.system('rm %s' %(outname))
    iraf.imutil.imcopy(inp,outname,verbose='no')
    
    return

### clean residuals ###
def cleanresidual(inp,out,x1,y1):
    '''
        cleanresidual - interactively clean residuals
        '''
    print "cleaning %s" %(inp)
    tempinp = tmppath+'inp.fits'
    tmpout = tmppath+'tmpedit.fits'
    os.system('cp %s %s' %(inp,tempinp))
    
    ds9.set('frame delete all')
    ds9.set('frame 1')
    ds9.set('scale zscale')
    ds9.set('file %s' %(inp))
    ds9.set('zoom to 4')
    #ds9.set('frame single')
    panto(ds9,x1,y1)
    
    #open temporary imedit cursor file#
    imcursor = open(tmppath+"tmp.cursor",'w')
    sline=''':aperture circular
        :search -2.
        :radius 2.
        :buffer 0.
        :width 2.
        :value 0.
        :sigma INDEF
        :xorder 2
        :yorder 2
        '''
    imcursor.write(sline)
    # Input image test_daof_sub
    sline = '# Input image %s' %(tempinp)
    imcursor.write(sline)
    
    menu = ''' # clean residuals #
        (n) select negative blemishes
        (p) select positive blemishes
        (v) constant value replacement
        (c) clean
        (u) undo all cleaning
        (q) done
        '''
    cleanstill = True
    cleaned = False
    while cleanstill:
        clean = raw_input(menu)
        # remove negative blemishes
        if clean == 'n':
            imcursor.write(':search 2.\n')
            imcursor.write(':value 0.\n')
            x0 = 9999.
            y0 = 9999.
            print 'click on residual/blemish position then click on residual radius'
            while True:
                crd = ds9.get('imexam coordinate ' +'image')
                C = crd.split()
                x = float(C[0])
                y = float(C[1])
                dx = x-x0
                dy = y-y0
                dr = (dx**2.+dy**2.)**0.5
                # stop if we click on the same star
                if (dr < 2.0):
                    break
                crd = ds9.get('imexam coordinate ' + 'image')
                C = crd.split()
                xr = float(C[0])
                yr = float(C[1])
                rad = ((x-xr)**2.+(y-yr)**2.)**0.5
                wcs_mark(ds9,x,y,major=str(rad),minor=str(rad),colour='red')
                sline = ':radius %i\n%i %i 1 b\n' %(rad,x,y)
                imcursor.write(sline)
                x0 = x
                y0 = y
        # remove positive blemishes
        elif clean == 'p':
            imcursor.write(':search -2.\n')
            imcursor.write(':value 0.\n')
            x0 = 9999.
            y0 = 9999.
            print 'click on positive blemish position then click on residual radius'
            while True:
                crd = ds9.get('imexam coordinate ' +'image')
                C = crd.split()
                x = float(C[0])
                y = float(C[1])
                dx = x-x0
                dy = y-y0
                dr = (dx**2.+dy**2.)**0.5
                # stop if we click on the same star
                if (dr < 2.0):
                    break
                crd = ds9.get('imexam coordinate ' + 'image')
                C = crd.split()
                xr = float(C[0])
                yr = float(C[1])
                rad = ((x-xr)**2.+(y-yr)**2.)**0.5
                wcs_mark(ds9,x,y,major=str(rad),minor=str(rad),colour='blue')
                
                sline = ':radius %i\n%i %i 1 b\n' %(rad,x,y)
                imcursor.write(sline)
                x0 = x
                y0 = y
        # constant value repacement
        elif clean == 'v':
            skyval= float(iraf.imstat(images=tempinp, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
            imcursor.write(':search 2.\n')
            imcursor.write(':value %f\n' %(skyval))
            #imcursor.write(':sigma %f\n' %(skyval))
            x0 = 9999.
            y0 = 9999.
            print 'click on residual/blemish position then click on residual radius'
            while True:
                crd = ds9.get('imexam coordinate ' +'image')
                C = crd.split()
                x = float(C[0])
                y = float(C[1])
                dx = x-x0
                dy = y-y0
                dr = (dx**2.+dy**2.)**0.5
                # stop if we click on the same star
                if (dr < 2.0):
                    break
                crd = ds9.get('imexam coordinate ' + 'image')
                C = crd.split()
                xr = float(C[0])
                yr = float(C[1])
                rad = ((x-xr)**2.+(y-yr)**2.)**0.5
                wcs_mark(ds9,x,y,major=str(rad),minor=str(rad),colour='red')
                sline = ':radius %i\n%i %i 1 e\n' %(rad,x,y)
                imcursor.write(sline)
                x0 = x
                y0 = y
        
        # replace inp image with original
        elif clean == 'u':
            if os.path.isfile(tempinp):
                os.system('cp %s %s' %(inp,tempinp))
                ds9.set('frame delete all')
                ds9.set('frame 1')
                ds9.set('file %s' %(tempinp))
            if os.path.isfile(tmppath+"tmp.cursor"):
                os.system('rm '+tmppath+"tmp.cursor")
        
        # run imedit
        elif clean == 'c':
            skyval= float(iraf.imstat(images=tempinp, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
            cleaned = True
            # command to write file
            imcursor.write(':write\n')
            imcursor.close()
            # have we run it before
            if os.path.isfile(tmpout):
                os.system('rm %s' %(tmpout))
            print 'runnning imedit now...'
            os.system('cat '+tmppath+"tmp.cursor")
            try:
                iraf.imedit(tempinp,tmpout,cursor=tmppath+"tmp.cursor",display='no',autodisplay='no',aperture='circular',search=-2.,buffer=0.,width=2.,xorder=2,yorder=2,sigma=skyval)
                ds9.set('frame new')
                ds9.set('file %s' %(tmpout))
                panto(ds9,x1,y1)
                os.system('cp -v %s %s' %(tmpout,out))
            except:
                'imedit error!!!'
            os.system('cp %s %s' %(tmpout,tempinp))
            # open new imcursor for further editing
            imcursor = open(tmppath+"tmp.cursor",'w')
            sline=''':aperture circular
                :search -2.
                :radius 2.
                :buffer 0.
                :width 2.
                :value 0.
                :sigma INDEF
                :xorder 2
                :yorder 2
                '''
            imcursor.write(sline)
            # Input image test_daof_sub
            sline = '# Input image %s' %(tempinp)
            imcursor.write(sline)
        
        
        else: #clean =='q':
            if not cleaned:
                conf = raw_input('You have not cleaned, continue? (y)')
                if conf == 'n':
                    continue
            imcursor.close()
            print 'done cleaning blemishes'
            # remove the temporary file
            if os.path.isfile(tmpout):
                os.system('rm %s' %(tmpout))
            cleanstill = False
    
    return

def autocleanresiduals(inp,xlist,ylist,ellipse):
    '''
        ### automatically clean residuals outside of masked galaxy###
        # input : inp - input image
        #	  x1,y1,maskrad - galaxy coords and square mask size
        # 	  xlist,ylist - list of coords where star has been removed
        '''
    #print 'inp=',inp
    #print 'xlist=',xlist
    #print 'ylist=',ylist
    #print 'ellipse=',ellipse
    tempinp = tmppath+'inp.fits'
    tmpout = tmppath+'tmpedit.fits'
    os.system('cp %s.fits %s' %(inp,tempinp))
    
    psfrad = 12.
    boxsize = 2
    skyval= float(iraf.imstat(images=inp, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    skysig= float(iraf.imstat(images=inp, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
    
    head = pf.getheader(tempinp)
    nx = head.get('NAXIS1')
    ny = head.get('NAXIS2')
    
    x0_e, y0_e, a_e, eps_e, phi_e = ellipse
    #ellipse_phi is deg cntclock from +y!! - ellipse eqn needs rad cntclock from +x
    phi_e = (90-phi_e)*pl.pi/180.
    #phi_e = (phi_e)*pl.pi/180.
    
    xres = []
    yres = []
    data = pf.getdata(tempinp)
    for i in range(len(xlist)):
        xd = int(xlist[i])
        yd = int(ylist[i])
        ai=pl.sqrt(((xd-x0_e)*pl.cos(phi_e)-(yd-y0_e)*pl.sin(phi_e))**2.+(((xd-x0_e)*pl.sin(phi_e)+(yd-y0_e)*pl.cos(phi_e))/eps_e)**2.)
        if ai < a_e:
            #print 'in ellipse'
            continue
        elif xd < 20 or xd > nx-20 or yd<20 or yd > ny-20:
            #print 'edge'
            continue
        else:
            sumval = 0.
            nval = 0.
            smin = skyval
            smax = skyval
            x1 = max(xd-1-boxsize,0.)
            y1 = max(yd-1-boxsize,0.)
            x2 = min(xd+boxsize,nx)
            y2 = min(yd+boxsize,ny)
            for x in range(x1,x2):
                for y in range(y1,y2):
                    if data[y,x] > smax: smax = data[y,x]
                    if data[y,x] < smin: smin = data[y,x]
            if smin < skyval-2.5*skysig: snfl = True
            else: snfl = False
            if smax > skyval+2.5*skysig: smfl = True
            else: smfl = False
            if smfl or snfl:
                xres += [xd]
                yres += [yd]
    print "autocleaning %s: removing %i residuals" %(inp,len(xres))
    
    #dummy = raw_input('wait')
    
    #open temporary imedit cursor file#
    imcursor = open(tmppath+"tmp.cursor",'w')
    # Input image test_daof_sub
    S=''':aperture circular
        :search 2.
        :radius 15.
        :buffer 0.
        :width 2.
        :xorder 2
        :yorder 2
        # Input image %s
        ''' %(tempinp)
    
    imcursor.write(':search 2.\n')
    for i in range(len(xres)):
        #negative blemishes
        sline = ':radius %i\n%i %i 1 b\n' %(psfrad,xres[i],yres[i])
        #print sline
        imcursor.write(sline)
    
    # run imedit
    imcursor.write(':write\n')
    imcursor.close()
    try:
        iraf.imedit(tempinp,tmpout,cursor=tmppath+"tmp.cursor",display='no',autodisplay='no',aperture='circular',radius=psfrad,search=-3.,buffer=3.,width=3.,value=skyval,sigma='INDEF')#skysig)
    except:
        'imedit error!!! did not remove all residuals'
    os.system('cp -v %s %s.fits' %(tmpout,inp))
    
    return

##################################################################
## make a postage stamp from j,h,k images at postiion ra,dec (pix)
def postage(image,x,y,N,rgb=True,rgbout='',inp='y',boxsize=10.):
    '''
        postage(field,xcenter,ycenter,backsize,rgb=True,rgbout=stamppath+field+'_'+str(linenum+1)+'.png')
        '''
    print 'making postage stamp for '+image
    C = image.split('_')
    field = C[0]
    num = C[1]
    print 'x=%.5f, y=%.5f, N=%i' %(x,y,N)
    headj=pf.getheader(imagepath+'gal_j'+image+'.fits')
    headh=pf.getheader(imagepath+'gal_h'+image+'.fits')
    headk=pf.getheader(imagepath+'gal_k'+image+'.fits')
    
    ximgsize = headj.get('NAXIS1')
    yimgsize = headj.get('NAXIS2')
    
    magzpj = headj.get('MAGZP')
    magzph = headh.get('MAGZP')
    magzpk = headk.get('MAGZP')
    
    offcentre = False
    # subimage limits: check if runs over edges
    xlim1 =  x - (N/2)
    if(xlim1<1):
        xlim1=1
        offcentre=True
    xlim2 =  x + (N/2)
    if(xlim2>ximgsize):
        xlim2=ximgsize
        offcentre=True
    ylim1 =  y - (N/2)
    if(ylim1<1):
        ylim1=1
        offcentre=True
    ylim2 =  y + (N/2)
    if(ylim2>yimgsize):
        offcentre=True
        ylim2=yimgsize
    
    xl = int(xlim1)
    yl = int(ylim1)
    xu = int(xlim2)
    yu = int(ylim2)
    print xl,yl,xu,yu
    print 'postage stamp is %i x %i pixels' %(xu-xl,yu-yl)
    # get xlimits, ylimits
    # cutout
    inpj = imagepath+'gal_j'+image+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    inph = imagepath+'gal_h'+image+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    inpk = imagepath+'gal_k'+image+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    outj = stamppath+'gal_%s%s_%s_starsub.fits' %('j',field,num)
    outh = stamppath+'gal_%s%s_%s_starsub.fits' %('h',field,num)
    outk = stamppath+'gal_%s%s_%s_starsub.fits' %('k',field,num)
    outjsub = stamppath+'galsub_%s%s_%s_starsub.fits' %('j',field,num)
    outhsub = stamppath+'galsub_%s%s_%s_starsub.fits' %('h',field,num)
    outksub = stamppath+'galsub_%s%s_%s_starsub.fits' %('k',field,num)
    
    # starry cutouts
    inpsj = imagepath+'gal_j'+field+'_'+num+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    inpsh = imagepath+'gal_h'+field+'_'+num+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    inpsk = imagepath+'gal_k'+field+'_'+num+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
    outsj = stamppath+'gal_%s%s_%s.fits' %('j',field,num)
    outsh = stamppath+'gal_%s%s_%s.fits' %('h',field,num)
    outsk = stamppath+'gal_%s%s_%s.fits' %('k',field,num)
    outsjsub = stamppath+'galsub_%s%s_%s.fits' %('j',field,num)
    outshsub = stamppath+'galsub_%s%s_%s.fits' %('h',field,num)
    outsksub = stamppath+'galsub_%s%s_%s.fits' %('k',field,num)
    
    newstamp = True
    if os.path.isfile(outj):
        if inp == 'y' and not AUTO: conf = raw_input('stamp exists, redo (y)')
        else: conf = 'y'
        if conf == 'n': newstamp = False
        else: os.system('rm %s*%s_%s*' %(stamppath,field,num))
    
    if newstamp:
        iraf.imutil.imcopy(inpj,outj)
        iraf.imutil.imcopy(inph,outh)
        iraf.imutil.imcopy(inpk,outk)
        iraf.imutil.imcopy(inpsj,outsj)
        iraf.imutil.imcopy(inpsh,outsh)
        iraf.imutil.imcopy(inpsk,outsk)
        skyj=float(iraf.imstatistics(outj, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        skyh=float(iraf.imstatistics(outh, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        skyk=float(iraf.imstatistics(outk, fields="mode", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        sigj=float(iraf.imstatistics(outj, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        sigh=float(iraf.imstatistics(outh, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        sigk=float(iraf.imstatistics(outk, fields="stddev", nclip=30, lsigma=3, usigma=3,  format=0,Stdout=1)[0])
        iraf.imarith(outj,'-',skyj,outjsub)
        iraf.imarith(outh,'-',skyh,outhsub)
        iraf.imarith(outk,'-',skyk,outksub)
        iraf.imarith(outsj,'-',skyj,outsjsub)
        iraf.imarith(outsh,'-',skyh,outshsub)
        iraf.imarith(outsk,'-',skyk,outsksub)

    if rgb:
        if len(rgbout)==0:
            outpic = stamppath+'rgb_'+field+'_'+num+'_starsub.jpg'
            outpics = stamppath+'rgb_'+field+'_'+num+'.jpg'
        else:
            outpic = rgbout
            outpics = rgbout
        makergb = True
        if os.path.isfile(outpic):
            if inp == 'y':
                if not AUTO: conf = raw_input('rgb stamp exists, redo (y)')
                else: conf = 'y'
            else: conf = 'n'
            if conf == 'n':
                makergb = False
        if os.path.isfile(outpics):
            if inp == 'y':
                if not AUTO: conf = raw_input('starry rgb stamp exists, redo (y)')
                else: conf = 'y'
            else: conf = 'n'
            if conf == 'n':
                makergb = False
        if makergb:
            #iraf.imgeom.magnify(outk,outkm,mag,mag)
            # get colour scaling based on 'normal' galaxy colours
            dat=pf.getdata(outksub)
            print outksub
            col,row = dat.shape
            print col,row
            x1=int((row/2-1)-boxsize)
            x2=int((row/2-1)+boxsize)
            y1=int((col/2-1)-boxsize)
            y2=int((col/2-1)+boxsize)
            sub=dat[x1:x2,y1:y2]
            print dat[x1:x2,y1:y2],boxsize,x1,x2,y1,y2
            maxsub=dat[x1:x2,y1:y2].max()  # K peak pixel
            
            xmax=0
            ymax=0
            for row in range(x1,x2):
                for col in range(y1,y2):
                    if dat[row,col]==maxsub:
                        xmax=row
                        ymax=col
            galcore=dat[xmax-2:xmax+3,ymax-2:ymax+3]
            #print galcore
            Kpeak = pl.average(galcore)
            Kpeakmag = -2.5*pl.log10(Kpeak)+magzpk
            JK=1.0  # J-H
            HK=0.28 # H-K
            Hpeakmag = HK+Kpeakmag
            Jpeakmag = JK+Kpeakmag
            Hpeak=10**((Hpeakmag-magzph)/-2.5)
            Jpeak=10**((Jpeakmag-magzpj)/-2.5)
            print Kpeak,Hpeak,Jpeak
            
            S='log'		#strech for rgb
            vlr=-0.075*Kpeak
            vur2=sigk
            vur=Kpeak
            vur=max(vur,vur2)
            vlg=-0.075*Hpeak
            vug2=sigh
            vug=Hpeak
            vug=max(vug,vug2)
            vlb=-0.075*Jpeak
            vub2=sigj
            vub=Jpeak
            vub=max(vub,vub2)
            print 'limits: ',vlr,vur,vlg,vug,vlb,vub,Kpeak,Hpeak,Jpeak
            print 'making RGB image: ',outksub
            RGB((outksub,outhsub,outjsub), outpic, stretch_r=S, vmin_r=vlr, vmax_r=vur, stretch_g=S, vmin_g=vlg, vmax_g=vug, stretch_b=S, vmin_b=vlb, vmax_b=vub)
            print 'test22-08'
            RGB((outsksub,outshsub,outsjsub), outpics, stretch_r=S, vmin_r=vlr, vmax_r=vur, stretch_g=S, vmin_g=vlg, vmax_g=vug, stretch_b=S, vmin_b=vlb, vmax_b=vub) #vmid_r
            if not AUTO:
                os.system('xv %s &' %(outpic))
                os.system('xv %s &' %(outpics))
    return


def subloop(field):
    '''
        subloop run psf making, starsubtraction and cleaning on each band ###
        input:  field name
        implicitly requires: hicat file (ws separated ra,dec,a,b,pa; comments |) listing sources in field
        j,h,k fits images in imagepath must be pixel ref'ed
        '''
    hicat = catpath+field+".hicat"
    print hicat
    if os.path.isfile(hicat):
        Hfile = open(hicat,'r')
        linenum = 0
        # read start info from file hicat and loop thru sources
        for gline in Hfile:
            if gline[0]=='|': continue  # ignore comments
            label = str(linenum+1)
            C = gline.split()
            ra = float(C[0])
            dec = float(C[1])
            a = float(C[2])
            b = float(C[3])
            pa = float(C[4])
            
            print '#### SOURCE %s ####' %(label)
            
            # pixel coords - ASSUMING J,H,K bands are pixel referenced
            bandname=['j','h','k']
            xP=np.zeros(3); yP=np.zeros(3)
            for ind in range(3):
                head = pf.getheader(imagepath+bandname[ind]+field+'.fits')
                #    headh = pf.getheader(imagepath+'h'+field+'.fits')
                #headk = pf.getheader(imagepath+'k'+field+'.fits')

                wcs = WCS(head)
                D = wcs.wcs_world2pix(pl.array([ra]),pl.array([dec]),1)
                xP[ind] = float(D[0])
                yP[ind] = float(D[1])

                print "DEBUG: ",ra,dec,a,b," ",imagepath+bandname[ind]+field+'.fits'    #  changed j to k
                print xP[ind], yP[ind]


                #print '(x,y) = (%.1f,%.1f)'%(xP,yP)
                radius1 = a*3600./0.45  # pixels
                radius2 = b*3600./0.45  # pixels
                #radius11 = (a+b)/2.  NH:used in postage stamps, not needed
                gradius=radius1*1.5	# maxsma for ellipse fitting
                backsize = 8*radius1	# size of the cutout
                if backsize < 200.:
                  backsize = 200.
                print "backsize: ",backsize
            
                # IMAGE CUTOUTS in each band #
                print 'making cutouts...'
                ximgsize = head.get('NAXIS1')
                yimgsize = head.get('NAXIS2')
                print 'ximgsize, yimgsize: ',ximgsize,yimgsize
                # subimage limits: check if runs over edges
                N = backsize
                xlim1 =  xP[ind] - (N/2)
                if(xlim1<1):
                  xlim1=1
                xlim2 =  xP[ind] + (N/2)
                if(xlim2>ximgsize):
                  xlim2=ximgsize
                ylim1 =  yP[ind] - (N/2)
                if(ylim1<1):
                  ylim1=1
                ylim2 =  yP[ind] + (N/2)
                if(ylim2>yimgsize):
                  ylim2=yimgsize
            
                xl = int(xlim1)
                yl = int(ylim1)
                xu = int(xlim2)
                yu = int(ylim2)
                print "limits (%s): "%(bandname[ind]),xl,xu,yl,yu

            # get xlimits, ylimits
            # cutout
            #for band in ['j','h','k']:
                inp = imagepath+bandname[ind]+field+'_smooth'
                if not os.path.isfile(inp+'.fits'): inp = imagepath+bandname[ind]+field+'[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
                else: inp = imagepath+bandname[ind]+field+'_smooth[%0.0f:%0.0f,%0.0f:%0.0f]' %(xl,xu,yl,yu)
                #inp = imagepath+band+field
                out = imagepath+'gal_%s%s_%s.fits' %(bandname[ind],field,label)

                print bandname[ind]," ",inp,out

                if os.path.isfile(out):
                    os.system('rm %s' %(out))
                iraf.imutil.imcopy(inp,out,verbose='no')

            '''
            print "DEBUG -- headers "
            headj = pf.getheader(imagepath+'gal_j'+field+'_'+label+'.fits')
                #      headh = pf.getheader(imagepath+'gal_h'+field+'_'+label+'.fits')
                #headk = pf.getheader(imagepath+'gal_k'+field+'_'+label+'.fits')

            wcs = pywcs.WCS(headj)
            D = wcs.wcs_sky2pix(pl.array([ra]),pl.array([dec]),1)
            x0 = float(D[0])
            y0 = float(D[1])
            print 'DDBUG (x,y) = (%.1f,%.1f)'%(x0,y0)
            '''


#            Display = True

            #dx = x0-xP
            #dy = y0-yP
            # display in ds9 and pan
            if Display:
                print " DO DISPLAY "
                dodisp(ds10,imagepath+'gal_j'+field+'_'+label+'.fits',0.,imagepath+'gal_h'+field+'_'+label+'.fits',0.,imagepath+'gal_k'+field+'_'+label+'.fits',0.,zoom='0.5')
                panto(ds10,ra,dec,csys='fk5')
                ds10.set('zoom to 2')
                ds10.set('single')

            print "DEBUG pstamp_sub = "
            pstamp_sub = stamppath+'rgb_'+field+'_'+label+'_starsub.jpg'
            pstamp = stamppath+'rgb_'+field+'_'+label+'.jpg'
            print pstamp
            print "HERE DEBUG: ",pstamp
            #if os.path.isfile(pstamp_sub): os.system('xv %s &' %(pstamp_sub))
            #if os.path.isfile(pstamp): os.system('xv %s &' %(pstamp))
            
            # STAR SUBTRACTION in each band #
            print 'starting starsubtraction'
            frame = 1
            for ind in range(3):
                print '\n ### ', bandname[ind]+field,' ###'
                n=0
                os.system('rm -r %s*' %(tmppath))
                
                ### automatic starsubtraction ###
                starsubin = imagepath+bandname[ind]+field+'_smooth'
                if not os.path.isfile(starsubin+'.fits'): starsubin = imagepath+bandname[ind]+field
                starsubout = imagepath+'gal_'+bandname[ind]+field+'_'+label+'_starsub.fits'
                inp = 'n'
                #if AUTO: inp = 'y'
                #else:
                if os.path.isfile(starsubout):
                    if AUTO: inp = 'y'
                    else: inp = str(raw_input('starsubtraction done, redo? (n) '))
                else:
                    if AUTO: inp = 'y'
                    else: inp = str(raw_input('no starsub file, run starsubtraction? (n) '))
                if inp == 'q': break
                if inp == 'y':
                    print "DEBUG running"

                    psfpart(corename=starsubin,outname=starsubout,psfname=psfpath+bandname[ind]+field+".psf.fits",galx=xP[ind],galy=yP[ind],size=backsize,verbose='no',def_psf='no',maxsma=gradius,nloop=5,a=a,b=b,pa=pa)

                    print "DEBUG psfpart done"

                    if Display:
                        ds10.set('frame '+str(frame))
                        ds10.set('file '+starsubout)
                        ds10.set('scale zscale')
                        panto(ds10,ra,dec,csys='fk5')
                        ds10.set('zoom to 2')
                        frame += 1
                        print "after display 1"
            if Display:
                dodisp(ds10,imagepath+'gal_j'+field+'_'+label+'_starsub.fits',0.,imagepath+'gal_h'+field+'_'+label+'_starsub.fits',0.,imagepath+'gal_k'+field+'_'+label+'_starsub.fits',0.)
                panto(ds10,ra,dec,csys='fk5')
                ds10.set('zoom to 4')
                ds10.set('tile no')
            
                print "after display 2"
            
            # RESIDUAL CLEANING in each band #
            if not AUTO: inp = raw_input('clean images (y) ')
            else: inp = 'n'
            if inp == 'q': break
            if inp != 'n':
                for ind in range(3):
                    starsubout = imagepath+'gal_'+bandname[ind]+field+'_'+label+'_starsub.fits'
                    if os.path.isfile(starsubout):
                        #### interactive cleaning of images for bad psf residuals ###
                        cleanout = imagepath+'gal_'+bandname[ind]+field+'_'+label+'_starsub.fits'
                        inp2 = 'y'
                        inp2 = str(raw_input('run starsubtraction cleaning? (y,n) '))
                        if inp2 != 'n':
                            if os.path.isfile(starsubout): cleanresidual(starsubout,cleanout,xP[ind],yP[ind])
                    else:
                        print 'no star subtracted image'
            
            if Display:
                dodisp(ds10,imagepath+'gal_j'+field+'_'+label+'_starsub.fits',0.,imagepath+'gal_h'+field+'_'+label+'_starsub.fits',0.,imagepath+'gal_k'+field+'_'+label+'_starsub.fits',0.)
                panto(ds10,ra,dec,csys='fk5')
            if not AUTO: savepost = raw_input('save postage stamp? (y) ')
            else: savepost = 'y'
            if savepost == 'y' or len(savepost) < 1:
                if radius1 < 30: postrad = 90.
                else: postrad = radius1*3
                postage(field+'_'+label+'_starsub',xP[ind],yP[ind],postrad,rgb=True)
            
            # increment catalogue line number
            linenum += 1      
    
    else:
        print 'no galaxy catalogue %s' %(hicat)
    return

def selectfield(field):
    #display
    if Display:
        dodisp(ds10,imagepath+'j'+field+'.fits',0.,imagepath+'h'+field+'.fits',0.,imagepath+'k'+field+'.fits',0.)
    
    # check for psfs, make them if nec
    print "checking for psf's"
    for band in ['j','h','k']:
        #make psf for field if not already made
        psf = psfpath+band+field+'.psf.fits'
        inp = 'n'
        if os.path.isfile(psf):
            if not AUTO: inp = str(raw_input('makepsf done, redo? (n) '))
            else: inp = 'n'
        #inp = 'n'
        else:
            print 'no psf file %s' %(psf)
            if not AUTO: inp = str(raw_input('run makepsf? (n) '))
            else: inp = 'y'
        if inp == 'y':
            if os.path.isfile(psf):
                os.system('rm ' + psf)
            try:
                mkpsf(corename=imagepath+band+field,outname=psfpath+band+field,verbose='no')
            except: print 'psf error'
    
    # SKY CLEANING in each band #
    if not AUTO: inp = raw_input('enter bands for sky gradient cleaning (jhk/none) ')
    else: inp = ''
    if len(inp) > 0:
        for band in ['j','h','k']:
            gradin = imagepath+band+field+'.fits'
            if os.path.isfile(gradin):
                ## removal of vertical feature in sky ##
                skyout = imagepath+band+field+'_smooth.fits'
                if band not in inp:
                    gradin = imagepath+band+field+'.fits'
                    skyout = imagepath+band+field+'_smooth.fits'
                    os.system('cp -v %s %s' %(gradin,skyout))
                else:
                    inp2 = 'n'
                    if os.path.isfile(skyout):
                        inp2 = str(raw_input('%s sky gradients done, redo? (y,n,N to save as done) ' %(band)))
                    else:
                        inp2 = 'y'
                    if inp2 == 'y' or len(inp2) <1:
                        if os.path.isfile(skyout): os.system('rm %s' %(skyout))
                        sky = float(iraf.imstatistics(gradin,fields='mode',nclip=30,format=0,Stdout=1)[0])
                        skysig = float(iraf.imstatistics(gradin,fields='stddev',nclip=30,format=0,Stdout=1)[0])
                        remove_vgrad(gradin,skyout,sky,skysig)
                    if inp2 == 'N':
                        if os.path.isfile(skyout):
                            os.system('cp -v %s %s' %(gradin,skyout))
            else:
                print 'no image:',gradin
    if Display:
        dodisp(ds10,imagepath+'j'+field+'_smooth.fits',0.,imagepath+'h'+field+'_smooth.fits',0.,imagepath+'k'+field+'_smooth.fits',0.,zoom='0.5')
    if not AUTO: dummy = raw_input('any key to continue')
    return


############################################################################
#############		  MAIN PROGRAM				############
############################################################################

if Display:
    ds9 = pyds9.ds9()
    ds10 = pyds9.ds9()
    Csystem = "galactic"
    ds9.set('tile yes')
    ds9.set('tile grid')
    ds10.set('tile yes')
    ds10.set('tile grid')

imagepath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/'
catpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/hicats/'
psfpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/psf/'
stamppath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/stamps/'

#inlist = 'objects.list'

#flist="galaxyfields.list"
#flist = 'calibrated_list'
#fieldlist1 = []
#infile = open(flist,'r')
#for line in infile:
#fieldlist1 += [line.strip()]
#infile.close()

flist='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated_list'
#flist = 'autocleanfail'
if os.path.isfile(flist):
    infile = open(flist,'r')
    fieldlist = []
    for line in infile:
        if line.strip() not in fieldlist:
            fieldlist += [line.strip()]
    infile.close()
    fieldlist = pl.array(fieldlist)
else:
    print 'no fieldlist: %s' %(flist)

lastfield = 'lastfield'
if os.path.isfile(lastfield):
    lastN = open(lastfield,'r')
    n = int(lastN.readline().strip())
    lastN.close()
else:
    n = 0

errorlog = 'subpsf.error'
if os.path.isfile(errorlog): errorfile = open(errorlog,'a')
else: errorfile = open(errorlog,'w')
errorfile.write('#started at %s\n' %(d.now()))
errorfile.close()

passlog = 'subpsf.log'
if os.path.isfile(passlog): passfile = open(passlog,'a')
else: passfile = open(passlog,'w')
passfile.write('#started at %s\n'%(d.now()))
passfile.close()

menu = '''Field Viewer
    n - select field by name
    b - select field by number
    a - automatically run makepsf
    p - automatically run postage
    q - quit
    option: 
    '''
while True:
    # loop by number
    if not AUTO: menu_inp = raw_input(menu)
    else: menu_inp = 'b'
    
    ### b - loop by number ###
    if menu_inp == 'b':
        print 'last field: %i ' %(n)
        if not AUTO: inp2 = raw_input('Enter field number: ')
        else: inp2 = ''
        if len(inp2) < 1: n += 1
        else:  n = int(inp2)
        field = fieldlist[n]
        print "field %i - %s selected" %(n,field)
        selectfield(field)
        # starsub loop for this field
        try:
            subloop(field)
            passfile = open(passlog,'a')
            passfile.write('%s : %s\n' %(d.now(),field))
            passfile.close()
        except:
            print 'field %s failed' %(field)
            errorfile = open(errorlog,'a')
            errorfile.write('%s : %s\n' %(d.now(),field))
            errorfile.close()
        # save last
        lastN = open(lastfield,'w')
        lastN.write('%i' %(n))
        lastN.close()
    
    ### n - loop by name ###
    if menu_inp == 'n':
        inp2 = raw_input('Enter field name: ')
        if inp2 not in fieldlist:
            print 'not in list'
            continue
        else:
            field = inp2
        n =pl.find(fieldlist==field)
        #field = fieldlist[n]
        print "field %i - %s selected" %(n,field)
        selectfield(field)
        # starsub loop for this field
        subloop(field)
    
    ### p - postage for all ###
    if menu_inp=='p':
        print '### running postage for all fields in galaxyfieldlist ###'    
        # loop thru fieldlist
        for field in fieldlist:
            #print '### %s ###' %(field)
            hicat = catpath+field+".hicat"
            print hicat
            if os.path.isfile(hicat):
                print hicat
                Hfile = open(hicat,'r')
                linenum = 0
                for line in Hfile:
                  if line[0]!='|':
                    #C2 = Hfile.readlines()[linenum].split()
                    C2 = line.split()
                    xcenter = float(C2[0])
                    ycenter = float(C2[1])
                    radius11 = float(C2[2])*3600/0.45
                    label = linenum+1
                    gradius=radius11*1.5
                    backsize = 4*radius11
                    #print line
                    jim=imagepath+'gal_j'+field+'_'+str(linenum+1)+'_starsub.fits'
                    him=imagepath+'gal_h'+field+'_'+str(linenum+1)+'_starsub.fits'
                    kim=imagepath+'gal_k'+field+'_'+str(linenum+1)+'_starsub.fits'
                    print jim
                    if os.path.isfile(jim) and os.path.isfile(him) and os.path.isfile(kim):
                        print 'making stamp'
                        postage(field+'_'+str(linenum+1)+'_starsub',xcenter,ycenter,backsize,rgb=True,inp='n')
                    linenum += 1
    
    
    ### a - psf for all ###
    if menu_inp=='a':
        print '### running mkpsf for all fields in galaxyfieldlist ###'
        successlogF = 'psfsuccess.log'
        faillogF = 'psffail.log'
        psflogF = 'psf.log'
        successlog = open(successlogF,'w')
        faillog = open(faillogF,'w')
        psflog = open(psflogF,'w')
        faillog.close()
        successlog.close()
        psflog.close()
        start = time.time()
        print start
        
        # loop thru fieldlist
        for field in fieldlist:
            successlog = open(successlogF,'a')
            faillog = open(faillogF,'a')
            psflog = open(psflogF,'a')
            print
            print '### %s ###' %(field)
            psfexist=[0,0,0] # no j,h,k psfs
            ### loop thru bands
            i=0
            for band in ['j','h','k']:
                print ' # ', band+field,' #'
                
                #make psf for field if not already made
                psf = psfpath+band+field+'.psf.fits'
                inp = 'n'
                if os.path.isfile(psf):
                    print 'psf made, continuing'
                    psfexist[i]=1
                else:
                    try:
                        mkpsf(corename=imagepath+band+field,outname=psfpath+band+field,verbose='no')
                        if os.path.isfile(psf):
                            print 'psf made, continuing'
                            psfexist[i]=1
                    except PsfError, pe:
                        print 'psf failed'
                        faillog.write(band+field+'\t'+pe.value+'\n')
                    except:
                        print 'psf failed'
                        faillog.write(band+field+'\t'+'unknown'+'\n')
                i+=1
            step = time.time()
            print "Time elapsed = %.2f minutes" %((step-start)/60.)
            if sum(psfexist) == 3:
                successlog.write(field+'\n')
            psflog.write("%s\t%i\t%i\t%i\n"%(field,psfexist[0],psfexist[1],psfexist[2]))
            faillog.close()
            successlog.close()
            psflog.close()
    
    
    ### q - quit ###
    elif menu_inp == 'q':
        os.system("killall xv")
        print 'Quit'
        break


