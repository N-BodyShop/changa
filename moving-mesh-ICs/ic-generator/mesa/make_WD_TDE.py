#import matplotlib
#matplotlib.use("Agg")
import math
import numpy as np
import subprocess
import pickle
import io
from scipy.interpolate import interp1d

import sys
sys.path.append('..')

from readWriteTipsy import *

import argparse
G = 6.67e-8
mSun = 1.99e33
rSun = 7e10
rho_max = 4e6
mStar = 0.
DEFAULT_NUM_PARTICLES = 100000
DEFAULT_BETA = 2

def rescaleExt( rStar, rParticle, rMax=3000.) : 
    scaleR = 2**np.minimum(np.maximum((rParticle/rStar - 1.)*6,0),np.log2(rMax))
    scaleR = np.maximum( scaleR, 1.)
    scaleR = np.minimum( scaleR, rMax)
    return scaleR

def filterMonotonicity( rho) :
    boolArray = np.ones(rho.size, dtype=bool)
    return boolArray
    mono = True
    for i in range(rho.size-1) :
        if( mono and rho[i] > rho[i+1]) :
            boolArray[i] = True
        else :
            mono = False
            boolArray[i] = False

    boolArray[-1] = boolArray[-2]
    return boolArray

def makeStarMesh(profileFile, numParticles, posArray, minRho = 0) :
    mcore, rcore, mModel, rModel, rhoModel, eModel = buildSimpleStar()
    if( profileFile != "simple") : 
        print("Loading " + profileFile)
        mcore, rcore, mModel, rModel, rhoModel, eModel = buildStar( profileFile, True)
        mcore = mcore*mSun
        rcore = rcore*rSun
        mModel = mModel*mSun
        rModel = rModel*rSun

    mass=mModel.max()
    rho_max = rhoModel.max()

    rArray = np.linalg.norm( posArray, axis=1)
    rCircle = rArray.max()/math.sqrt(3.)

    # find the number of particles in the circle
    boolArray = np.argsort(rArray)
    rArray = rArray[boolArray]
    posArray = posArray[boolArray]
    boolArray = rArray < rCircle
    rArray = rArray[boolArray]
    posArray = posArray[boolArray]

    eArray = np.zeros(rArray.size)
    vArray = np.zeros([rArray.size, 3])
    XArray = np.zeros(rArray.size)

    # now rescale
    rMax = rModel.max()
    rArray = rArray*rMax/rCircle*(rArray.size/numParticles)**0.33333
    mArray = np.ones(rArray.size)
    interpolate = interp1d(rModel,rhoModel,kind = 'linear')
    interpolateE = interp1d(rModel,eModel,kind = 'linear')
    #interpolateX = interp1d(rModel,XModel,kind = 'linear')

    lastModelP = 0
    for i in range(numParticles) :
        if( rArray[i] > rMax) :
            break
        vol = 4.*math.pi/3.*rMax**3/numParticles
        mArray[i] = max( interpolate( rArray[i]), minRho)*vol
        eArray[i] = interpolateE( rArray[i])
        #XArray[i] = interpolateX( rArray[i])
        posArray[i,:] = posArray[i,:]/np.linalg.norm(posArray[i,:])*rArray[i]
        lastModelP = lastModelP+1

    return posArray[:lastModelP], vArray[:lastModelP], eArray[:lastModelP], XArray[:lastModelP], mArray[:lastModelP], rho_max, mcore, rcore, mass

def buildStar( profileFile="profile15.data", makePlot=False) :
    import mesa_reader as mr
    import matplotlib.pyplot as pl
    mcore = None
    rcore = None
    m = None
    r = None
    rho = None
    T = None
    p = None
    Xout = None
    data = mr.MesaData(profileFile)
    logT   = data.logT
    mass   = data.mass
    logR   = data.logR
    logRho = data.logRho
    logP   = data.logP

    R = 1e1**logR
    T = np.maximum(1e1**logT, 1e6)
    mtot = mass.max()
    maxR = R.max()

    rho = 1e1**logRho
    if( makePlot ) :
        pl.clf()
        pl.semilogy(R, 1e1**logRho, lw=2, color="blue")

    meanRho = mtot*mSun/(4.*math.pi/3.*(maxR*rSun)**3)

    mcore, rcore = 0, 0 

    return mcore, rcore, mass, R, rho, T

def buildSimpleStar() : 
    N=int(1e5)
    mue = 2
    gamma=5/3
    k=1e13/(mue**gamma)
    rhoc=1.5e6
    Pc=k*rhoc**gamma
    dr=1e6
    kB = 1.38e16
    mp = 1.6e-24
    
    mcore, rcore = 0, 0 
    #initialize radius, pressure, density, and mass arrays
    r=np.zeros(N) 
    P=np.zeros(N)  
    M=np.zeros(N)
    rho=np.zeros(N)

    #starting points
    P[0]=Pc
    P[1]=Pc
    r[0]=0
    r[1]=dr
    rho[1]=(P[1]/k)**(1/gamma)
    rho[0]=(P[0]/k)**(1/gamma)
    M[1]=4*np.pi/3*rho[1]*dr**3

    #fill in the rest
    i=2

    while P[i-1]>(1e-6*Pc):
        P[i]=P[i-1]-rho[i-1]*G*M[i-1]/(r[i-1]**2)*dr
        M[i]=M[i-1]+4*np.pi*r[i-1]**2*rho[i-1]*dr
        rho[i]=(P[i]/k)**(1/gamma)
        r[i]=r[i-1]+dr
        i=i+1

    #resize arrays, find mass and radius
    P=P[0:i]
    M=M[0:i]
    rho=rho[0:i]
    r=r[0:i]
    Mtot=M[-1]
    R=r[-1]
    T = mue*mp*P/rho/kB
    print(rho)
    return mcore, rcore, M, r, rho, T


parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('profilefile', default="profile70.data", help="<profilename>")
parser.add_argument('--num_particles', default=DEFAULT_NUM_PARTICLES, type=int, help='number of particles for star')
parser.add_argument('--beta', default=DEFAULT_BETA, type=float, help='beta of tde')
parser.add_argument('--no_atm', action='store_true', help='no atmosphere surrounding object')
parser.add_argument('--no_TDE', action='store_true', help="use the pickle", default=False)
parser.add_argument('--no_BH', action='store_true', help="use the pickle", default=False)
args = parser.parse_args()

N = 24 
posArray, velArray = makePosVel( N, useGlass=True)

# now center
numParticles = args.num_particles
massParticleInGram = mSun/numParticles

xmax = 2e11
rhoExt = 1e-2
minRhoExt = 1e-5
tempExt = 5e6

ymax = xmax
zmax = xmax

p1, v1, e1, X1, m1, rho_max, mcore, rcore, mStar = makeStarMesh( args.profilefile, numParticles, posArray, minRho = rhoExt)

pos  = p1
vel  = v1
ener = e1
mass = m1
Hfrac = X1

r1Max = np.linalg.norm(p1[-1])
print( "Number of star particles: {0} with size {1}".format(e1.size, r1Max))
r1 = np.array([0.,0.,0.])

# check xmax
if( xmax < r1Max) :
    print( 'warning')
# do the external gas

if( not args.no_atm) :
    # find the mean separation
    dr1 = np.linalg.norm( p1[0:-2,:] - p1[-1,:], axis=1)
    dr1.sort()
    #h = 0.1*r1Max #
    h = dr1[0:8].mean()

    massParticleExtInGram = 4.*3.1415/3.*rhoExt*h**3
    massParticleExtInCode = massParticleExtInGram
    minMassParticleExtInGram = 4.*3.1415/3.*minRhoExt*h**3
    minMassParticleExtInCode = massParticleExtInGram

    # set the fiducial external mass -- get the spacing right near the star
    totalExtMass = rhoExt * 8. * r1Max**3
    numberExtParticles = totalExtMass/massParticleExtInCode

    #find the scaling
    print( "ratio={0}".format((posArray.size/3/numberExtParticles)**0.3333))
    scaling = (posArray.size/3/numberExtParticles)**0.3333
    posArray = posArray*r1Max/(0.5*N)*scaling
    rPosArray = np.linalg.norm(posArray - r1, axis=1)

    # now rescale to get the rest correct 
    scaling = rescaleExt( r1Max, rPosArray) #, xmax/r1Max)

    ExtPosArray = posArray*scaling[:,np.newaxis]
    # cut the cube first
    boolArray = ExtPosArray[:,0] < xmax
    boolArray = np.logical_and( boolArray, ExtPosArray[:,0] > -xmax)

    boolArray = np.logical_and( boolArray, ExtPosArray[:,1] < ymax)
    boolArray = np.logical_and( boolArray, ExtPosArray[:,1] > -ymax)

    boolArray = np.logical_and( boolArray, ExtPosArray[:,2] < zmax)
    boolArray = np.logical_and( boolArray, ExtPosArray[:,2] > -zmax)

    ExtPosArray = ExtPosArray[boolArray]
    scaling = scaling[boolArray]

    # now exclude the regions around star

    boolArray = np.linalg.norm(ExtPosArray - r1, axis=1) > r1Max

    ExtPosArray = ExtPosArray[boolArray]
    print( "h = {0}".format(h/7e10))
    print( scaling)
    print( r1Max/7e10)
    print( xmax/7e10)
    # create the other arrays
    NumExtParticles = int(ExtPosArray.size/3)
    ExtVelArray = np.zeros([NumExtParticles,3])
    ExtMArray = np.maximum(massParticleExtInCode*np.ones(NumExtParticles), minMassParticleExtInCode*np.ones(NumExtParticles)*scaling[boolArray]**3)
    ExtEArray = tempExt*np.ones(NumExtParticles)
    ExtXArray = 0.7*np.ones( NumExtParticles)

    print( ExtPosArray)
    print( ExtMArray.min(), ExtMArray.max())
    # now append
    pos = np.append(pos, ExtPosArray, axis=0)
    vel = np.append(vel, ExtVelArray, axis=0)
    ener = np.append(ener, ExtEArray)
    mass = np.append(mass, ExtMArray)
    Hfrac = np.append( Hfrac, ExtXArray)

x =  pos[:,0]
y =  pos[:,1]
z =  pos[:,2]
vx = vel[:,0]
vy = vel[:,1]
vz = vel[:,2]

print( "Number of total particles: {0}".format(mass.size))

xboolArray = np.logical_and( x < xmax, x > -xmax)
yboolArray = np.logical_and( y < ymax, y > -ymax)
zboolArray = np.logical_and( z < zmax, z > -zmax)

boolArray = np.logical_and( np.logical_and(xboolArray, yboolArray), zboolArray)
x = x[boolArray]
y = y[boolArray]
z = z[boolArray]

vx = vx[boolArray]
vy = vy[boolArray]
vz = vz[boolArray]

mass = mass[boolArray]
#rhos = rhoArray[boolArray]
temp = ener[boolArray]
metals = Hfrac[boolArray]
#pArray = pArray[boolArray]
print( "rho_max = ", rho_max)
print( "mcore = ", mcore)

xcore = None
ycore = None
zcore = None
vxcore = None
vycore = None
vzcore = None


alpha = 25.5
beta = args.beta
eta = 2.5
bhsplit = 8

mbh = alpha**3*mStar
rtidal = alpha * r1Max
rbh = 0.1*rtidal  # 0.1 of rtidal
tdyn = G**-0.5
vstar = math.sqrt(2.*G*mStar/r1Max)*tdyn
print( "mbh = {0} {1} {2}".format( mbh, rbh, vstar/1e5/tdyn))
r0 = eta * rtidal
zbh = 0. #-r0
ybh = 0.
xbh = 0.
vxbh =0. #alpha*vstar/eta/math.sqrt(beta)
vybh = 0.
vzbh = 0. #alpha*vstar*math.sqrt( (eta*beta - 1.)/(eta*eta*beta))
print( "z = {0} vx = {1}  vz ={2}".format(-r0, vxbh, vzbh))
if( mcore <= 0.) :
    mcore = np.ones([bhsplit])*mbh/bhsplit
    rcore = np.ones([bhsplit])*rbh
    xcore = np.ones([bhsplit])*xbh
    ycore = np.ones([bhsplit])*ybh
    zcore = np.ones([bhsplit])*zbh
    vxcore = np.ones([bhsplit])*vxbh
    vycore = np.ones([bhsplit])*vybh
    vzcore = np.ones([bhsplit])*vzbh
    
else :
    mcore = np.array([mcore, mbh])
    rcore = np.array([rcore, rbh])
    xcore = np.array([0.,xbh])
    ycore = np.array([0.,ybh])
    zcore = np.array([0.,zbh])
    vxcore = np.array([0.,vxbh])
    vycore = np.array([0.,vybh])
    vzcore = np.array([0.,vzbh])


#mcore = np.array([mcore])
#mcore = np.array([])
if not args.no_TDE : 
    z[:] += r0
    z[z>=zmax] -= 2.*zmax
    vx[:] = alpha*vstar/eta/math.sqrt(beta)
    vz[:] = -alpha*vstar*math.sqrt( (eta*beta - 1.)/(eta*eta*beta))
    r = np.sqrt(x*x + y*y + z*z)
    vz[numParticles:] = -np.sqrt( 2*G*mbh/np.maximum(r[numParticles:],rbh))*tdyn
    print(vz)

print("Number of gas particles: ", z.size)
if args.no_BH or args.no_TDE:
    # reset the black hole
    mcore = np.array([])

print("mcore ", mcore)
#zero out the gas
#x = np.array([])
#y = x
#z = x
#vx = x
#vy = x
#vz = x
#mass = x
#temp = x
#metals = x
if( mcore.size > 0.) :
    write_tipsy( "tde", x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals,
        xdm=xcore, ydm=ycore, zdm=zcore, vxdm=vxcore, vydm=vycore, vzdm=vzcore, mdm=mcore, dmeps=rcore)
else : 
    write_tipsy( "tde", x, y, z, vx, vy, vz, mass, temp, rho_max, gasmetals=metals)
