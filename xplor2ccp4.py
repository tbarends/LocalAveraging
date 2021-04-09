import struct
import numpy as np
import glob

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
#######################WRITECCP4MAP_FULL####################
#
def writeccp4map_full(mapoutname,SIZE,START,INTERVALS,CELL,ORDER,SKEW,SKEWTRN,SYMOPS,maparray):
    HEADER=bytes()
    
    print 'Now in writeccp4map function'
    print 'SIZE looks like:',SIZE
    for value in SIZE:
        HEADER += struct.pack('i',value)
    print 'Setting the mode to 2...'    
    HEADER += struct.pack('i',int(2)) #SET THE MODE TO 2
    print 'I think START looks like:',START
    for value in START:
        HEADER += struct.pack('i',value)
    print 'I think INTERVALS looks like:',INTERVALS
    for value in INTERVALS:
        HEADER += struct.pack('i',value)
    print 'I think CELL looks like:',CELL

    for value in CELL:
        HEADER += struct.pack('f',value)
    print 'I think ORDER looks like:',ORDER
    for value in ORDER:
        HEADER += struct.pack('i',value)
        
    AMIN=np.min(maparray)
    AMAX=np.max(maparray)
    AMEAN=np.average(maparray)
    print 'Min,max,av:',AMIN,AMAX,AMEAN
    HEADER += struct.pack('3f',AMIN,AMAX,AMEAN)
    #print 'I think INTERVALS looks like:',INTERVALS
    #for value in INTERVALS:
    #    HEADER += struct.pack('i',value)
    
    ISPG=1
    NSYMBT=0
    LSKFLG=1 # IF SKEW IS TO BE WRITTEN OUT
    #LSKFLG=0
    
    HEADER += struct.pack('3i',ISPG,NSYMBT,LSKFLG)

    SKEW=np.reshape(SKEW,(9))
    #SKEW=np.zeros((9))
    SKEWTRN=np.zeros((3))
    print 'I think the SKEW looks like:',SKEW
    for value in SKEW:
        #print value
        HEADER += struct.pack('f',value)

    print 'I think the SKEWTRN looks like:',SKEWTRN
    for value in SKEWTRN:
        HEADER += struct.pack('f',value)
    print 'adding 15 zeros...'
    for n in range(15):
        print n,
        HEADER += struct.pack('i',int(0))
    print ''
    print 'Adding the word MAP_'
    HEADER += struct.pack('4c','M','A','P',' ')
    HEADER += struct.pack('4c',chr(0x44), chr(0x41), chr(0x00), chr(0x00)) 
    # little-endian machine stamp 'D','A',0,0
    

    ARMS=np.std(maparray)
    HEADER += struct.pack('f',ARMS)
    HEADER += struct.pack('i',3)
    for n in range(800):
        HEADER += struct.pack('c',' ')

    mapoutfile=open(mapoutname, "wb")
    mapoutfile.write(HEADER)
    mapsize=np.shape(maparray)
    print 'The map size is',mapsize
    maplength=(mapsize[0]*mapsize[1]*mapsize[2])
    print 'The map length is',maplength
    maparray=np.reshape(maparray,(maplength))
    for n in range(maplength):
        voxelout=struct.pack('f',maparray[n])
        mapoutfile.write(voxelout)
    mapoutfile.close()
##############################################################################   
    
##################################################
#
# Read an XPLOR format map file
def readxplormap(mapinname):
    mapin=open(mapinname,'r')
    title=mapin.readline()
    #print 'Map Title:',title
    
    line1=mapin.readline() # works for XPLOR maps from PHENIX, some programs write a bigger header...
    line2=mapin.readline() # need to make general...
    
    #print 'Encountered end of remarks...'
    #read number of slices and their start and end numbers
    line=mapin.readline()
    elements=line.split()
    #print elements
    NA=int(elements[0])
    AMIN=int(elements[1])
    AMAX=int(elements[2])
    NB=int(elements[3])
    BMIN=int(elements[4])
    BMAX=int(elements[5])
    NC=int(elements[6])
    CMIN=int(elements[7])
    CMAX=int(elements[8])
    
    #read unit cell parameters:
    line=mapin.readline()
    elements=line.split()
    A=float(elements[0])
    B=float(elements[1])
    C=float(elements[2])
    ALPHA=float(elements[3])
    BETA=float(elements[4])
    GAMMA=float(elements[5])

    #read dimension order:
    line=mapin.readline()

    first=line[0]
    if first=='X':
        NSECT=1+AMAX-AMIN
    if first=='Y':
        NSECT=1+BMAX-BMIN
    if first=='Z':
        NSECT=1+CMAX-CMIN
 
    second=line[1]
    if second=='X':
        MEDIUM=1+AMAX-AMIN
    if second=='Y':
        MEDIUM=1+BMAX-BMIN
    if second=='Z':
        MEDIUM=1+CMAX-CMIN

    third=line[2]
    if third=='X':
        FAST=1+AMAX-AMIN
    if third=='Y':
        FAST=1+BMAX-BMIN
    if third=='Z':
        FAST=1+CMAX-CMIN

    sectionsize=(MEDIUM)*(FAST)
    densmap=np.zeros((FAST,MEDIUM,NSECT))

    for n in range(NSECT):
        section=[]
        line=mapin.readline()
        #print 'Now reading section #',line
        incomplete=True
        while (incomplete):
            line=mapin.readline()
            nfields=int(len(line)/12)
            for m in range(nfields):
                start=m*12
                end=start+12
                pixel=float(str(line[start:end]))
                section.append(pixel)
                if len(section)==sectionsize:
                    incomplete=False
                    
        #print 'Read',len(section),'voxels into this section'
        section=np.reshape(section,(FAST,MEDIUM),order='F')
        densmap[:,:,n]=section
    mapin.close()
    return densmap,A,B,C,ALPHA,BETA,GAMMA,NA,AMIN,AMAX,NB,BMIN,BMAX,NC,CMIN,CMAX  

################################################################################

################################MAIN###################################
path='./'
filenames=glob.glob((path+'*.*'))
print filenames
for filename in filenames:
    if filename[-6:]=='.xplor':
        xplorin=filename
        ccp4out=filename[:-6]+'.ccp4'


        totalmap,a,b,c,alpha,beta,gamma,na,mina,maxa,nb,minb,maxb,nc,minc,maxc = readxplormap(xplorin)
        SKEW=skew(a,b,c,alpha,beta,gamma)
        SIZE=[1+maxa-mina,1+maxb-minb,1+maxc-minc]
        START=[mina,minb,minc]
        INTERVALS=[na,nb,nc]
        CELL=[a,b,c,alpha,beta,gamma]
        ORDER=[3,2,1] # XPLOR ZYX
        SKEWTRN=[0.,0.,0.]
        SYMOPS=''

        print 'Writing out the ccp4 map',ccp4out
        writeccp4map_full(ccp4out,SIZE,START,INTERVALS,CELL,ORDER,SKEW,SKEWTRN,SYMOPS,totalmap)
        print 'Appears to be done :)'
