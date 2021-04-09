################SETTINGS######################################
import sys
root         = sys.argv[1]
pathname     = '/home/tbarend/FAP/FAP_20200622_Fo-Fo_with_pdb_6ZH7/Fo-Fo/'
mtzinname    = pathname+root               # mtz file name
amplitude    = 'FOFCWT'
phase        = 'PHFOFCWT'
pdbinname    = pathname+'6ZH7.pdb'               # input coordinate file
mapoutname   = pathname+root+'.av.xplor'              				  # output map file
pdboutname   = pathname+root+'.av.pdb'
a            = 61.39                                     # unit cell parameters
b            = 60.01
c            = 182.90
alpha        = 90.0
beta         = 90.6
gamma        = 90.0
radius       = 30                                          # half the (box length - 1) in grid points
step         = 0.5                                        # grid spacing in Angstrom
chains       = ['A','B']                                       # the chain identifier
RESs         = [801]                                      # a python list of residue numbers

####################import modules############################
#!/usr/bin/env python
import numpy as np
import string
#from mpl_toolkits.mplot3d import Axes3D
#import matplotlib.pyplot as plt
#from matplotlib import cm
import scipy.io as io
import os

##################ReadPdb##########################
# read a pdb file (duh...)
#
#
#ATOM    293  CG  PHE A  38      19.222  -0.685  10.342  1.00  8.64           C
#0123456789012345678901234567890
#0         1         2         3
def readpdb(infile):
    pdb=open(infile,'r')
    x=[]
    y=[]
    z=[]
    name=[]
    element=[]
    resn=[]
    resid=[]
    chain=[]
    atid=[]
    for line in pdb:
        if line[0:4]=='ATOM' or line[0:6]=='HETATM':
            #print line
            x.append(float(line[28:38]))
            y.append(float(line[38:46]))
            z.append(float(line[46:54]))
            name.append(line[13:16])
            element.append(line[76:78])
            resn.append(line[17:20])
            chain.append(line[21])
            resid.append(float(line[22:26]))
	    atid.append(line[:26])
    x=np.asarray(x)
    y=np.asarray(y)
    z=np.asarray(z)
    C=np.zeros((len(x),3))
    C[:,0]=x
    C[:,1]=y
    C[:,2]=z
    pdb.close()
    return C,name,resn,chain,resid,element,atid 



##################ElectronDensity##################
# Calculate the (unnormalized) electron density
# at a point cfrac (a one dimensional vector of
# shape (3) containing x,y,z in fractional
# coordinates) from a list of h,k,l,F and phi
# No correction for unit cell volume is made at
# this point!
def ElectronDensity(cfrac,h,k,l,F,phi):
    hprod=h*cfrac[0]
    kprod=k*cfrac[1]
    lprod=l*cfrac[2]

    prod=2.0*np.pi*(hprod+kprod+lprod)
    shift=np.cos(phi-prod)

    rho=np.sum(np.multiply(F,shift))
    return rho

###################################################

#######################skew########################
# Calculate the skew or scale matrix from unit cell
# parameters
#
def skew(a,b,c,alpha,beta,gamma):
    pi=3.14159265359
    alpha=pi*alpha/180.0
    beta=pi*beta/180.0
    gamma=pi*gamma/180.0

    vsq=1.0 - (np.cos(alpha)*np.cos(alpha)) - (np.cos(beta)*np.cos(beta)) - (np.cos(gamma)*np.cos(gamma)) - (2*np.cos(alpha)*np.cos(beta)*np.cos(gamma))
    v=np.sqrt(vsq) # volume of unit parallelepiped

    
    s=np.zeros((3,3))
    s[0,0]=1.0/a
    s[0,1]=-1.0*np.cos(gamma)/(a*np.sin(gamma))
    s[0,2]=(np.cos(alpha)*np.cos(gamma)-np.cos(beta)) / (a*v*np.sin(gamma))
    s[1,1]=1.0/(b*np.sin(gamma))
    s[1,2]=(np.cos(beta)*np.cos(gamma)-np.cos(alpha)) / (b*v*np.sin(gamma))
    s[2,2]=np.sin(gamma)/(c*v)
    return s

###################################################

######################CellVolume###################
# Calculate the unit cell volume from unit cell
# parameters
#
def CellVolume(a,b,c,alpha,beta,gamma):
    pi=3.14159265359
    alpha=pi*alpha/180.0
    beta=pi*beta/180.0
    gamma=pi*gamma/180.0

    left=(1.0 - np.square(np.cos(alpha)) - np.square(np.cos(beta))- np.square(np.cos(gamma))  )
    right=2.0*np.cos(alpha)*np.cos(beta)*np.cos(gamma)


    V=a*b*c*np.sqrt(left+right)

    return V
###################################################

####################ReadHkl########################
# Read a hkl file with format h k l F phi
#
def readhkl(hklinname):
    hklin=open(hklinname,'r')            # open map coefficients

    h=[]                                 # declare empty lists for storage
    k=[]
    l=[]
    F=[]
    phi=[]


    for line in hklin:                  # loop over lines in file and read
        elements=line.split()           # map coefficients
        h.append(int(elements[0]))
        k.append(int(elements[1]))
        l.append(int(elements[2]))
        F.append(float(elements[3]))
        phi.append(float(elements[4]))

    hklin.close()                       # close file again

    h=np.asarray(h)                     # turn python lists into numpy arrays
    k=np.asarray(k)
    l=np.asarray(l)
    F=np.asarray(F)
    phi=np.asarray(phi)
    return h,k,l,F,phi

#######
#
#######
def readmtzcolumns_P1(mtzin, columnnames):
    columnnamesstring=" "
    for n in range(len(columnnames)):
        columnnamesstring=columnnamesstring+columnnames[n]+" "
    os.system("rm -f temporary.hkl")
    infile=open('sftoolstmpin.com', 'w')
    infile.write('#!/bin/bash \n')
    infile.write('sftools <<eof >tmp.out\n')
    infile.write('read '+mtzin+' \n')
    infile.write('Y \n')
    infile.write('expand \n')
    infile.write('calc col mult = M\n')
    infile.write('calc col stol = STOL\n')
    infile.write('write temporary.hkl col mult stol'+columnnamesstring+'\n')
    
    infile.write('format(3i4,'+str(int(2+len(columnnames)))+'f16.7)\n')
    infile.write('DEGREES'+'\n')
    infile.write('STOP'+'\n')
    infile.write('eof'+'\n')
    infile.close()
    os.chmod('sftoolstmpin.com',int(0b111111000)) # rwxrwx---
    os.system("./sftoolstmpin.com")
    #os.system("rm -f sftoolstmpin.com tmp.out")
    infile=open("temporary.hkl","r")
    linesin=infile.readlines()
    infile.close()
    h=[]
    k=[]
    l=[]
    m=[]
    stol=[]
    values=np.empty([len(columnnames),len(linesin)])
    v=0
    for line in linesin:
        elements=line.split()
        h.append(int(elements[0]))
        k.append(int(elements[1]))
        l.append(int(elements[2]))
        m.append(float(elements[3]))
        stol.append(float(elements[4]))
        for c in range(len(columnnames)):
            values[c,v]=float(elements[c+5])    
        v=v+1
    #os.system("rm -f temporary.hkl sftoolstmpin.com tmp.out")
    h=np.asarray(h)
    k=np.asarray(k)
    l=np.asarray(l)
    m=np.asarray(m)
    stol=np.asarray(stol)
    return h,k,l,m,stol,values







#######################writexplormap#############
#
# Write an x-plor format map to file
def writexplormap(mapoutname,totalmap,a,b,c,alpha,beta,gamma,na,mina,maxa,nb,minb,maxb,nc,minc,maxc):
    size=np.shape(totalmap)
    mapout=open(mapoutname,'w')
    mapout.write('      2 !NTITLE\n')
    mapout.write('REMARKS FILENAME=  \n')
    mapout.write('REMARKS DATE:xx-Xxx-xx  xx:xx:xx       created by user: \n')
    NA="%8i" % na
    AMIN="%8i" % mina
    AMAX="%8i" % maxa
    NB="%8i" % nb
    BMIN="%8i" % minb
    BMAX="%8i" % maxb
    NC="%8i" % nc
    CMIN="%8i" % minc
    CMAX="%8i" % maxc
    mapout.write(NA+AMIN+AMAX+NB+BMIN+BMAX+NC+CMIN+CMAX+"\n")
    A="%12.5e" % a
    B="%12.5e" % b
    C="%12.5e" % c
    ALPHA="%12.5e" % alpha
    BETA="%12.5e" % beta
    GAMMA="%12.5e" % gamma
    mapout.write(A+B+C+ALPHA+BETA+GAMMA+"\n")
    mapout.write("ZYX\n")



    sectionsize=size[0]*size[1]
    sectionparts=range(0,sectionsize,6)
    last=sectionparts.pop()

    for z in range(size[2]):
        SECT="%8i" % z
        mapout.write(SECT+'\n')
        section=totalmap[:,:,z]
        section=np.reshape(section,(sectionsize),order='F')

        for n in sectionparts:
            line=""
            for m in range(n,n+6):
                pixel="%12.5e" % section[m]
                line=line+pixel
            mapout.write(line+"\n")
        line=""
        for m in range(last,sectionsize):
            pixel="%12.5e" % section[m]
            line=line+pixel
        mapout.write(line+"\n")

    mapout.write('-9999\n')
    SIGMA="%12.4e" % np.std(totalmap)
    AVER="%12.4e" % np.average(totalmap)
    mapout.write(SIGMA+AVER+'\n')       
    mapout.close()
##################################################

########################MAIN PROGRAM##############
print '###################AVERAGE A DENSITY MAP AROUND PEPTIDES####################'
print ''

scale=skew(a,b,c,alpha,beta,gamma)   # calculate scale matrix
print "The scale matrix is"
print scale
print ''
V=CellVolume(a,b,c,alpha,beta,gamma) # calculate unit cell volume
print "The unit cell volume is",V,"cubic Angstrom"
print ''

#h,k,l,F,phi=readhkl(hklinname)
print 'I will work with',mtzinname
print 'and use columns',amplitude,phase
h,k,l,m,stol,values=readmtzcolumns_P1(mtzinname,[amplitude,phase])

F=values[-2,:]
phi=values[-1,:]

print 'Read',len(F),'structure factor amplitudes from',mtzinname
print 'The first record on the hkl file looks like this:',h[0],k[0],l[0],F[0],phi[0]

print ""
avphi=np.average(phi)
if (avphi>170) and (avphi<190):
    print 'The numerical value of the average phase angle is',avphi
    print '--> these are probably degrees, I will convert them...'
    print ""
    phi=3.14159265359*(phi/180.0) #if it is in degrees


if (avphi>3.0) and (avphi<3.3):
    print 'The numerical value of the average phase angle is',avphi
    print '--> these are probably radians, no conversion needed.'
    print ""

if (avphi>0.9) and (avphi<1.1):
    print 'The numerical value of the average phase angle is',avphi
    print '--> these are probably radians/pi, I will convert them.'
    print ""
    phi=np.pi*phi

avphi=np.average(phi)

print 'The average angle of the phases in memory is',np.average(phi)/np.pi,'pi radians.'
print 'This value should be around 1.0'
print ""







C,name,resn,chain,resid,element,atid=readpdb(pdbinname)

mapsize=(2.0*radius)+1.0
mapsize=mapsize*step
print 'The map will have sides of',mapsize,'Angstrom'
print 'I will use coefficients',amplitude,phase
totalmap=np.zeros((2*radius+1,2*radius+1,2*radius+1))
for chain_search in chains:
    for resid_search in RESs:
        Cphe=np.zeros((4,3))
        for n in range(len(name)):
            if resid[n]==resid_search and chain_search==chain[n]:
                if name[n]=='C1 ':
                    Cphe[1,0]=C[n,0]
                    Cphe[1,1]=C[n,1]
                    Cphe[1,2]=C[n,2]
                if name[n]=='C2 ':
                    Cphe[2,0]=C[n,0]
                    Cphe[2,1]=C[n,1]
                    Cphe[2,2]=C[n,2]
                if name[n]=='O1 ':
                    Cphe[0,0]=C[n,0]
                    Cphe[0,1]=C[n,1]
                    Cphe[0,2]=C[n,2]
                if name[n]=='O2 ':
                    Cphe[3,0]=C[n,0]
                    Cphe[3,1]=C[n,1]
                    Cphe[3,2]=C[n,2]
                   
        print 'I found the atoms of ',chain_search,resid_search,'.'
        start,end=1,2
        vector1=np.asarray([ Cphe[end,0]-Cphe[start,0] , Cphe[end,1]-Cphe[start,1] , Cphe[end,2]-Cphe[start,2]])
        xvec=vector1/np.linalg.norm(vector1)
        start,end=1,0
        vector2=np.asarray([ Cphe[end,0]-Cphe[start,0] , Cphe[end,1]-Cphe[start,1] , Cphe[end,2]-Cphe[start,2]])
        
        approx_yvec=vector2/np.linalg.norm(vector2)
        zvec=np.cross(xvec,approx_yvec)
        zvec=zvec/np.linalg.norm(zvec)
        yvec=-1.0*np.cross(xvec,zvec)
        yvec=yvec/np.linalg.norm(yvec)
        print 'My x vector (from C1 to C2) will be',xvec
        print 'My y-prime vector (from O1 to O2 ) will be',approx_yvec
        print 'My z vector (cross product of x and y-prime) will be',zvec
        print 'My final y vector (cross product of x and z) will be',yvec
        xmid=np.average(Cphe[:,0])
        ymid=np.average(Cphe[:,1])
        zmid=np.average(Cphe[:,2])
        mid=np.asarray([xmid,ymid,zmid])
        print 'The center of the CH2COOH is at',mid

        emap=np.zeros((2*radius+1,2*radius+1,2*radius+1))
        for xg in range(2*radius+1):
            for yg in range(2*radius+1):
                for zg in range(2*radius+1):
                    xgm=xg-radius
                    ygm=yg-radius
                    zgm=zg-radius
                    pos=mid+(xgm*step*xvec)+(ygm*step*yvec)+(zgm*step*zvec)
                    cfrac=np.mod(np.dot(scale,pos),1.0)
                    emap[xg,yg,zg]=ElectronDensity(cfrac,h,k,l,F,phi)                            
        totalmap=totalmap+emap
        
totalmap=totalmap/float(len(RESs))
totalmap=totalmap/np.std(totalmap)
#totalmap=totalmap-np.average(totalmap)


size=np.shape(totalmap)
na=int(size[0])-1
nb=int(size[1])-1
nc=int(size[2])-1
mina=0
minb=0
minc=0
maxa=na
maxb=nb
maxc=nc

writexplormap(mapoutname,totalmap,mapsize,mapsize,mapsize,90.0,90.0,90.0,na,mina,maxa,nb,minb,maxb,nc,minc,maxc)

##########NEW FUNCTION###############
contourlevel=3.0
negcopy=np.array(totalmap, copy=True)
poscopy=np.array(totalmap, copy=True)
negcopy[negcopy>(-1.0*contourlevel)]=0.0
poscopy[poscopy<contourlevel]=0.0
negint=np.sum(negcopy, axis=None)
posint=np.sum(poscopy, axis=None)
print 'Using a contour level of',contourlevel
print 'The integrated positive value is',posint
print 'And the integrated negative value is',negint
#####################################




print 'Wrote XPLOR format map',mapoutname

#names=['ATOM      1  C6  FAD A   1','ATOM      2  C9  FAD A   1','ATOM      3  C9  FAD A   1','ATOM      4  C5X FAD A   1']

pdbout=open(pdboutname,'w')
origin=np.asarray([xmid,ymid,zmid])
for n in range(len(atid)):
    atom=np.asarray([C[n,0],C[n,1],C[n,2]])
    atom=atom-origin
    atomx=(radius+.5)*step+np.dot(atom,xvec)
    atomy=(radius+.5)*step+np.dot(atom,yvec)
    atomz=(radius+.5)*step+np.dot(atom,zvec)
    na=atid[n]
    x="%12.3f" % atomx
    y="%8.3f"  % atomy
    z="%8.3f"  % atomz
    o="  1.00"
    bval=20.0
    b="%6.2f" % bval
    s="          "
    e=" "
    outstring=na+x+y+z+o+b+s+e+'\n'
    outlist=list(outstring)
    #outlist[16]='A'
    outstring=''.join(outlist)
    pdbout.write(outstring)
pdbout.close()
print 'Write pdb file',pdboutname

print ''
print 'To load into pymol enter the following commands:'
print '------------------------------------------------'
print 'reinitialize'
print 'load',pdboutname,',pdb'
print 'load',mapoutname,',avmap'   



