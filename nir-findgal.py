import os
wd=os.getcwd()
os.chdir('/Volumes/BigDalek/khaled/')
from pyraf import iraf
os.chdir(wd)

#import pysao
from ds9 import*
import pylab as pl
import math
import asciidata as ad
#from useful import afind


Csystem='fk5'

def afind(find, incolumn):
    index=[]
    for i in range(len(incolumn)):
     if find==incolumn[i]:
        index=index+[i]
    return index

def showfield(name):
    '''
    showfield - displays j,h,k images for field 'name'
    input: - name
    '''
    print 'displaying field ',name

    kname=imagepath+'k'+name+'.fits'
    hname=imagepath+'h'+name+'.fits'
    jname=imagepath+'j'+name+'.fits'
    if os.path.isfile(kname) and os.path.isfile(hname) and os.path.isfile(jname):

	ds9.set('frame delete all')
	# get sky vals
	ksky=float(iraf.imstatistics(kname,fields='mode',nclip=30,format=0,Stdout=1)[0]) # compute and print image pixel statistics
	hsky=float(iraf.imstatistics(hname,fields='mode',nclip=30,format=0,Stdout=1)[0])
	jsky=float(iraf.imstatistics(jname,fields='mode',nclip=30,format=0,Stdout=1)[0])

	#set limits
	Rlim=str(ksky-7)+' '+str(ksky+100)
	Glim=str(hsky-8)+' '+str(hsky+150)
	Blim=str(jsky-8)+' '+str(jsky+100)

	scale='sqrt'

	ds9.set('single')
	ds9.set('zoom 1')
	ds9.set('rgb')
	ds9.set('rgb red')
	ds9.set('scale '+scale)
	ds9.set('scale limits '+Rlim)
	ds9.set('file '+kname)
	ds9.set('rgb green')
	ds9.set('scale '+scale)
	ds9.set('scale limits '+Glim)
	ds9.set('file '+hname)
	ds9.set('rgb blue')
	ds9.set('scale '+scale)
	ds9.set('scale limits '+Blim)
	ds9.set('file '+jname)
    else:
	print "missing files"
    return


def wcs_mark(x,y,major='0.0025',minor='0.0025',angle='0',frame='1',colour='white',csys=Csystem):
    '''
    wcs_mark - mark ellipse on ds9
    '''
    S=csys+"; ellipse "+x+' '+y+' '+major+' '+minor+' '+angle+" # color = "+ colour
    ds9.set('frame '+frame)
    ds9.set('regions ',S)
    return

def showmarkers(name):
    '''
    showmarkers - load ds9 regions for field name
    '''
    catname=catpath+name+'.cat'
    if os.path.isfile(catname):
	print 'load catalogue: ',catname
	frm=ds9.get('frame').strip()
    radius=0.0025 #in degrees 20pixels=9arcsecond/3600
	major=str(radius)
	minor=str(radius)
	angle="0"
	file1=open(catname,'r')
	dat=file1.readlines()
	for line in dat:
	    C=line.split()
	    if len(C)>0 and C[0]!='|' and C[0]!='#':
		wcs_mark(C[0],C[1],major,minor,angle,frm,'cyan')

	file1.close()
    else:
	print 'no catalogue file: ',catname
    
    hicatname=catpath+name+'.hicat'
    if os.path.isfile(hicatname):
	print 'load catalogue: ',hicatname
	frm=ds9.get('frame').strip()
	file1=open(hicatname,'r')
	dat=file1.readlines()
	for line in dat:
	    C=line.split()
	    if len(C)>0 and C[0]!='|' and C[0]!='#':
		wcs_mark(C[0],C[1],C[2],C[3],C[4],frm,'blue')

	file1.close()
    else:
	print 'no HI catalogue file: ',hicatname
    return

def showHIdat(field):
    '''
    showHIdat - print HI data for field fieldname
    '''
    ind=afind(field,masterdat['NAME'])
    if len(ind)==0:
	print 'error, not in masterlist'
	return
    ind=ind[0]
    pos=str(masterdat['L'][ind])+', '+str(masterdat['B'][ind])
    print '************************************************'
    print '%9s: %-15s  %9s: %-15s' %('Name',masterdat['NAME'][ind],'flux',masterdat['FLUX'][ind])
    print '%9s: %-15s  %9s: %-15s' %('l,b',pos,'HI mass',masterdat['H1MASS'][ind])
    print '%9s: %-15s  %9s: %-15s' %('E(B-V)',masterdat['E(B-V)'][ind],'velocity',masterdat['VHEL'][ind])
    print '%9s: %-15s  %9s: %-15s' %('A(V)',masterdat['AV'][ind],'w50',masterdat['W50'][ind])
    print
    #if masterdat['SPECTRUM'][ind]=='yes':
    spfile='/Volumes/BigDalek/khaled/nir_pipeline_NHK/lists/HI/spectra/'+field+'.ps'
    if os.path.isfile(spfile):
	os.system('xv -rotate -90 -expand 0.5 '+spfile+' &') 
    else:
	print 'no HI spectrum (%s)' %(masterdat['CUBE'][ind])
    if masterdat['2MASS'][ind]!='\N':
	list2mass=masterdat['2MASS'][ind].split(';')
	for i in range(len(list2mass)):
	    print list2mass[i]
	    tmfile='/Volumes/BigDalek/khaled/nir_pipeline_NHK/lists/2MASS/images/'+list2mass[i]+'.jpg'
	    os.system('xv '+tmfile+' &') 
	    xra=list2mass[i][6:14]
	    xdec=list2mass[i][15:]
	    sra=str(15*(int(xra[:2])+int(xra[2:4])/60.+int(xra[4:6])/3600.+float(xra[6:])/360000.))
	    sdec=str((int(xdec[:2])+int(xdec[2:4])/60.+int(xdec[4:6])/3600.+float(xdec[6:])/360000.))
	    if list2mass[i][14]=='-':
		sdec='-'+sdec
	    frm=ds9.get('frame').strip()
	    wcs_mark(sra,sdec,major='0.003',minor='0.003',frame=frm,colour='magenta',csys='fk5')
    else:
	print 'no 2mass counterparts'
    print '************************************************'
    return

def plotHIZOA(field):
    '''
    plotHIZOA - make image showing position of target in galactic coords on HIZOA grid
    '''
    ind=afind(field,masterdat['NAME'])
    if len(ind)==0:
	print 'error, not in masterlist'
	return
    ind=ind[0]
    pl.figure(1,figsize=(7,3))
    pl.clf()
    pl.subplots_adjust(right=0.95,bottom=0.2)
    pl.xlabel('Galactic Longitude $l$ [deg]')
    pl.ylabel('Galactic Latitude $b$ [deg]')
    l=masterdat['L'].tonumpy()
    l=l+360.*(l<=150)
    b=masterdat['B'].tonumpy()
    pl.plot(l,b,'b.')
    LI=masterdat['L'][ind]+360.*(masterdat['L'][ind]<=150)
    pl.plot(LI,masterdat['B'][ind],'ro')
    pl.xlim(420,190)
    pl.savefig('tempfield.png',dpi=60)
    os.system('xv tempfield.png &')
    return


def loadfieldlist(fieldlistname):
    '''
    loadfieldlist -  read list of fields from file fieldlist if it exists
    '''
    fieldlist=[]
    if os.path.isfile(fieldlistname):
	fieldlistflag=True
	fieldlistF=open(fieldlistname,'r')
	for line in fieldlistF:
	    fieldlist+=[line.split()[0]]
	print  '%i fieldnames from file "%s" loaded' %(len(fieldlist),fieldlistname)
    else:
	print 'error: no such file ',fieldlistname
    return fieldlist

######################################################

# open ds9
ds9 = ds9()
ds9.set('wcs sky '+Csystem)
ds9.set('wcs skyformat degrees')

if len(os.sys.argv)==1:
    imagepath='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/'
    catpath='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated/hicats/'
elif len(os.sys.argv)==2:
    imagepath="./"
    field=os.sys.argv[1]
    showfield(field)
    catpath="./"
else:
    print "too many arguments"
    os.sys.exit()

pspath='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/ps/'

#masterlist='/Volumes/BigDalek/khaled/nir_pipeline_NHK/lists/masterHI'
#masterdat=ad.open(masterlist,delimiter='|')

Lastfield=-1
fieldlistname='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/calibrated_list'
lastfieldlist = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/lastfield'
if os.path.isfile(lastfieldlist):
    filelast = open(lastfieldlist,'r')
    sLast = filelast.readline()
    try:
	Lastfield = int(sLast)
    except ValueError:
	print 'error reading last field'
	Lastfield = -1
    filelast.close()
    savefieldlast=True
else:
    Lastfield = -1

fieldlist=loadfieldlist(fieldlistname)
savefieldlast=True
fieldlistflag=(len(fieldlist)>0)


#################################################
# menu
text='''### Galaxy Locator ###
n - new field
t - change target list
c - add comment to field
l - change limits
g - mark galaxy
h - mark likely HI counterpart
p - print source list
q - quit
'''



while True:
    option=str(raw_input(text))
    
    ### new field ###
    if option=='n':
	os.system('killall xv')
	if fieldlistflag:
	    print 'Last field ',Lastfield
	    T=str(raw_input('Enter field number: '))
	    i = len(T)
	    if i <= 0:
		fieldnum=Lastfield+1
	    else:
		fieldnum=int(T)
	    if fieldnum>len(fieldlist)-1:
		print 'no more fields in list'
		continue
	    field=fieldlist[fieldnum]
	    Lastfield=fieldnum
	else:
	    field=str(raw_input('Enter field name: '))
	
	showfield(field)
	showmarkers(field)
#	showHIdat(field)
#	plotHIZOA(field)
        print "debug 2"

    ### change targetlist ###
    elif option == 't':
	fieldlistname = str(raw_input('New fieldlist: '))
	if fieldlistname == 'all':
	    fieldlistname = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/images1007/list1007'
	    lastfieldlist = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/images1007/list1007_last'
	    catpath = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/images1007/'
	if fieldlistname in ['norma','Norma']:
	    fieldlistname='/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRAC/IRSF/normafieldlist'
	    savelastfield=False
	    Lastfield=-1
	if os.pat.isfile(lastfieldlist):
	    filelast = open(lastfieldlist,'r')
	    sLast = filelast.readline()
	    Lastfield = int(sLast)
	    filelast.close()
	    savefieldlast=True
	else:
	    Lastfield = -1
	    savefieldlast=True
	fieldlist = loadfieldlist(fieldlistname)
	fieldlistflag = ( len(fieldlist) > 0 )

    ### change limits ###
    elif option=='l':
	print 'Current limits'
	ds9.set('rgb red')
	print 'red ',ds9.get('scale limits').strip()
	ds9.set('rgb green')
	print 'green ',ds9.get('scale limits').strip()
	ds9.set('rgb blue')
	print 'blue ',ds9.get('scale limits').strip()
	color=str(raw_input('Enter color (rgb): '))
	lim=str(raw_input('Enter new limits (z1 z2): '))
	if color in ['r','red']:
	    ds9.set('rgb red')
	    Rlim=lim
	    ds9.set('scale limits '+Rlim)
	elif color in ['g','green']:
	    ds9.set('rgb green')
	    Glim=lim
	    ds9.set('scale limits '+Glim)
	elif color in ['b','blue']:
	    ds9.set('rgb blue')
	    Blim=lim
	    ds9.set('scale limits '+Blim)

    ### add field comment ###
    elif option=='c':
	if len(ds9.get('file').strip())==0:
	    print 'no field loaded'
	    continue
	cmmt=str(raw_input('Field comment (focus, sky, seeing, combine, other): '))
	ind=afind(field,masterdat['NAME'])
	if len(ind)==0:
	    print 'error, not in masterlist'
	    continue
	ind=ind[0]
	masterdat['IRSF_COMMENT'][ind]=cmmt
	masterdat.flush()

    ### mark galaxies ###
    elif option=='g':
	if len(ds9.get('file').strip())==0:
	    print 'no field loaded'
	    continue

	frm=ds9.get('frame').strip()
	print 'choose galaxy'
    
	catname=catpath+field+'.cat'
	if os.path.isfile(catname):
	    file1=open(catname,'a')
	else:
	    file1=open(catname,'w')
	    file1.write('|   ra    |     dec   |\n')
	
	ra0 = 0.
	dec0 = 0.
	while True:
	    crd=' '
	    print 'click on galaxy center'
	    crd=ds9.get('imexam coordinate '+Csystem)
	    L=len(crd)

	    C=crd.split()
	    if len(C)==0:
		print 'imexam break'
		break
	    ra=float(C[0])
	    dec=float(C[1])
	    dra=(ra-ra0)*3600.*math.cos(dec/57.2957795)
	    ddec=(dec-dec0)*3600.
	    dr=(dra**2 + ddec**2)**0.5
	    if (dr<2.0):
		break

	    radius=0.0025
	    major=str(radius)
	    minor=str(radius)
	    angle="0"
	    wcs_mark(C[0],C[1],major,minor,angle,frm,"cyan")
	    print 'ra,dec = ',ra,dec

	    ra0 = ra
	    dec0 = dec
	    Sline = ' %.5f  %.5f\n'%(ra,dec)
	    file1.write(Sline)

	file1.close()


    ### mark HI counterpart ###
    elif option=='h':
	if len(ds9.get('file').strip())==0:
	    print 'no field loaded'
	    continue

	frm=ds9.get('frame').strip()
	print 'choose HI counterpart'
    
	catname=catpath+field+'.hicat'
	if os.path.isfile(catname):
	    file1=open(catname,'a')
	else:
	    file1=open(catname,'w')
	    file1.write('|   ra    |     dec   |     a     |     b     |    phi    |\n')
	
	ra0 = 0.
	dec0 = 0.
	while True:
	    crd=' '
	    print 'click on galaxy center'
	    crd=ds9.get('imexam coordinate '+Csystem)
	    L=len(crd)

	    C=crd.split()
	    ra=float(C[0])
	    dec=float(C[1])
	    dra=(ra-ra0)*3600.*math.cos(dec/57.2957795)
	    ddec=(dec-dec0)*3600.
	    dr=(dra**2 + ddec**2)**0.5
	    if (dr<2.0):
		break
	    
	    print 'click on semi-major axis'
	    crd=ds9.get('imexam coordinate '+Csystem)
	    L=len(crd)
	    C=crd.split()
	    ara=float(C[0])
	    adec=float(C[1])
	    adra=(ara-ra)*math.cos(dec/57.2957795)
	    addec=(adec-dec)
	    phi=math.atan(abs(addec/adra))*57.2957795
	    phi+=90.
	    if addec/adra<0:
		phi=90.+phi
	    else:
		phi=90.-phi
	    print adra,addec,phi
	    a=(adra**2 + addec**2)**0.5

	    print 'click on semi-minor axis'
	    crd=ds9.get('imexam coordinate '+Csystem)
	    L=len(crd)
	    C=crd.split()
	    bra=float(C[0])
	    bdec=float(C[1])
	    bdra=(ra-bra)*math.cos(dec/57.2957795)
	    bddec=(dec-bdec)
	    b=(bdra**2 + bddec**2)**0.5

	    radius=0.0025
	    major=str(a)
	    minor=str(b)
	    angle=str(phi)
	    wcs_mark(str(ra),str(dec),major,minor,angle,frm,"blue")
	    print 'ra,dec = ',ra,dec

	    ra0 = ra
	    dec0 = dec
	    Sline = ' %.5f  %.5f  %.5f  %.5f  %.5f\n'%(ra,dec,a,b,phi)
	    file1.write(Sline)

	file1.close()

    ### print the source catalogue ###
    elif option=='p':
	if len(ds9.get('file').strip())==0:
	    print 'no field loaded'
	    continue
	catname=catpath+field+'.cat'
	if os.path.isfile(catname):
	    os.system('cat '+catname)

    ### quit ###
    elif option=='q':
	print 'bye...'
	os.system('killall xv')
	if savefieldlast:
	    print 'saving last field'
	    lastfieldlist = '/Volumes/BigDalek/khaled/nir_pipeline_NHK/IRSFdata/images1007/list1007_last'
	    filelast = open(lastfieldlist,'w')
	    filelast.write(str(Lastfield))
	    filelast.close()
	break
    
    else:
	continue
