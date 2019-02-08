####################import modules############################
#!/usr/bin/env python
#
import numpy as np
import string
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt


##################ReadPdb##########################
# read a pdb file (duh...)
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
    x=np.asarray(x)
    y=np.asarray(y)
    z=np.asarray(z)
    C=np.zeros((len(x),3))
    C[:,0]=x
    C[:,1]=y
    C[:,2]=z
    pdb.close()
    return C,name,resn,chain,resid,element 



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

    prod=2*np.pi*(hprod+kprod+lprod)
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

    left=(1 - np.square(np.cos(alpha)) - np.square(np.cos(beta))- np.square(np.cos(gamma))  )
    right=2*np.cos(alpha)*np.cos(beta)*np.cos(gamma)


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


for root in ['20','40','60','80','100']:

    hklinname    = 'D:/PUMPPROBE/FOFO_SAMENUMBERIMAGES/THAU/HKL/thau_'+root+'fs.hkl'                 # text file with h k l F phi in P1
    pdbinname    = 'D:/PUMPPROBE/PUMPROBE_13.1-coot-0_refmac1.pdb'   # input coordinate file
    mapoutname   = 'D:/PUMPPROBE/FOFO_SAMENUMBERIMAGES/THAU/thau_phe_'+root+'fs.xplor'               # output map file
    pdboutname   = 'D:/PUMPPROBE/FOFO_SAMENUMBERIMAGES/THAU/thau_phe.pdb'
    pseoutname   = 'D:/PUMPPROBE/FOFO_SAMENUMBERIMAGES/THAU/thau_phe_'+root+'fs.pse'
    a            = 58.0                                      # unit cell parameters
    b            = 58.0 
    c            = 150.10
    alpha        = 90.0
    beta         = 90.0
    gamma        = 90.0
    radius       = 4                                          # half the (box length - 1) in grid points
    step         = 0.75                                        # grid spacing in Angstrom
    chain_search = 'A'                                       # the chain identifier

    PHEs=[3,58,80,90,108,114,152,172,173,181,203]

    ########################MAIN PROGRAM##############

    scale=skew(a,b,c,alpha,beta,gamma)   # calculate scale matrix
    print "The scale matrix is"
    print scale
    V=CellVolume(a,b,c,alpha,beta,gamma) # calculate unit cell volume
    print "The unit cell volume is",V

    h,k,l,F,phi=readhkl(hklinname)
    print 'Read',len(F),'structure factor amplitudes from',hklinname
    print 'The first record looks like this:'
    print h[0],k[0],l[0],F[0],phi[0]
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

    C,name,resn,chain,resid,element=readpdb(pdbinname)

    mapsize=(2*radius)+1
    mapsize=mapsize*step
    print 'The map will have sides of',mapsize,'Angstrom'



    Cphe=np.zeros((6,3))
    totalmap=np.zeros((2*radius+1,2*radius+1,2*radius+1))

    for resid_search in PHEs:

        for n in range(len(name)):
            if resid[n]==resid_search and chain_search==chain[n]:
                if name[n]=='CG ':
                    Cphe[0,0]=C[n,0]
                    Cphe[0,1]=C[n,1]
                    Cphe[0,2]=C[n,2]
                if name[n]=='CD1':
                    Cphe[1,0]=C[n,0]
                    Cphe[1,1]=C[n,1]
                    Cphe[1,2]=C[n,2]
                if name[n]=='CD2':
                    Cphe[2,0]=C[n,0]
                    Cphe[2,1]=C[n,1]
                    Cphe[2,2]=C[n,2]
                if name[n]=='CE1':
                    Cphe[3,0]=C[n,0]
                    Cphe[3,1]=C[n,1]
                    Cphe[3,2]=C[n,2]
                if name[n]=='CE2':
                    Cphe[4,0]=C[n,0]
                    Cphe[4,1]=C[n,1]
                    Cphe[4,2]=C[n,2]
                if name[n]=='CZ ':
                    Cphe[5,0]=C[n,0]
                    Cphe[5,1]=C[n,1]
                    Cphe[5,2]=C[n,2]
        print 'I found the atoms of ',chain_search,resid_search,'.'
        vector1=np.asarray([ Cphe[5,0]-Cphe[0,0] , Cphe[5,1]-Cphe[0,1] , Cphe[5,2]-Cphe[0,2]])
        xvec=vector1/np.linalg.norm(vector1)
        #print 'The vector between CG and CZ is determined as',vector1
        vector2=np.asarray([ Cphe[4,0]-Cphe[3,0] , Cphe[4,1]-Cphe[3,1] , Cphe[4,2]-Cphe[3,2]])
        #print 'The vector between CE2 and CE1 is determined as',vector2
        yvec=vector2/np.linalg.norm(vector2)
        zvec=np.cross(xvec,yvec)
        print 'My x vector will be',xvec
        print 'My y vector will be',yvec
        print 'My z vector will be',zvec
        xmid=np.average(Cphe[:,0])
        ymid=np.average(Cphe[:,1])
        zmid=np.average(Cphe[:,2])
        mid=np.asarray([xmid,ymid,zmid])
        print 'The center of the ring is at',mid

        x=[]
        y=[]
        z=[]
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
                    
                    x.append(pos[0])
                    y.append(pos[1])
                    z.append(pos[2])
                                
        #fig = plt.figure()
        #ax = fig.add_subplot(111, projection='3d')
        #ax.scatter(x,y,z)
        #ax.scatter(Cphe[:,0],Cphe[:,1],Cphe[:,2],c='r')
        #plt.show()



       
        #emap=emap-np.average(emap)

        totalmap=totalmap+emap


               
    totalmap=totalmap/float(len(PHEs))
    totalmap=totalmap/np.std(totalmap)

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

    #plt.imshow(totalmap[:,:,int(nc/2.0)])
    #plt.show()

    x=[]
    y=[]
    outpdb=open(pdboutname,'w')
    for a in range(0,360,60):
        angle=np.pi*a/180.0
        xn=(radius+.5)*step+((2.77/2.0)*np.cos(angle))
        yn=(radius+.5)*step+(-1.0*(2.77/2.0)*np.sin(angle))
        zn=(radius+.5)*step+0.0


        xs="%12.3f" % xn
        ys="%8.3f"  % yn
        zs="%8.3f"  % zn
        os="  1.00"
        bval=20.0
        bs="%6.2f" % bval
        ss="          "
        es=" "
        outstring='ATOM    438  C   PHE A   1'+xs+ys+zs+os+bs+ss+es+'\n'
        outpdb.write(outstring)
    #write Cbeta coordinates
    angle=np.pi*180.0/180.0
    xn=(radius+.5)*step+((2.77/2.0)*np.cos(angle))
    yn=(radius+.5)*step+(-1.0*(2.77/2.0)*np.sin(angle))
    zn=(radius+.5)*step+0.0
    xn=xn-1.51
    xs="%12.3f" % xn
    ys="%8.3f"  % yn
    zs="%8.3f"  % zn
    os="  1.00"
    bval=20.0
    bs="%6.2f" % bval
    ss="          "
    es=" "
    outstring='ATOM    438  C   PHE A   1'+xs+ys+zs+os+bs+ss+es+'\n'
    outpdb.write(outstring)

    #done with pdb output file
    outpdb.close()
    print 'Wrote pdb file',pdboutname
    print ''
    print 'To load into pymol enter the following commands:'
    print '------------------------------------------------'
    print 'reinitialize'
    print 'load '+pdboutname+',pdbin'
    print 'load '+mapoutname+',mapin'
    print 'isomesh pos,mapin,3.0'
    print 'color green, pos'
    print 'isomesh neg, mapin,-3.0'
    print 'color red, neg'
    print 'save '+pseoutname
