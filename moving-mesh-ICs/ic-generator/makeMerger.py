import numpy as np
from scipy.interpolate import interp1d
import myhelmholtz as helm
import scipy.integrate as si
import math
import argparse
import h5py

from readWriteTipsy import *

G = 6.67384e-8
kB = 1.38e-16
mp = 1.6e-24

def adbDerivs( r, y, K, gamma) :
    P = y[0]
    m = y[1]
    rho = (P/K)**(1./gamma)

    dPdr = -rho*G*m/(r*r)
    dMdr = 4.*math.pi*rho*r*r

    return [dPdr, dMdr]

def makeAdStar( mass, gamma=1.667, temp=1e7, tol=1e-3, n=1) :
    rhoc = 30.
    # compute central pressure

    m = 0.
    r = 0.
    mArray = None
    rArray = None
    eArray = None
    rhoArray = None
    pArray = None
    Tc = temp
    gamma = (1. + n)/(1.*n)
    while abs((m-mass)/mass) > tol :
        integral = si.ode( adbDerivs).set_integrator("dop853")
        Pc = kB*Tc/mp*rhoc
        K = Pc/rhoc**gamma
        r = 1e8
        m = 4.*math.pi/3. * rhoc*r**3
        P = Pc - 2.*math.pi/3.*G*r**2*rhoc**2
        y0 = np.array([P, m])
        dr = r*0.1
        integral.set_initial_value(y0, r).set_f_params( K, gamma)
        mArray = [0.]
        rArray = [0.]
        eArray = [Tc]
        rhoArray = [rhoc]
        pArray = [Pc]
        while integral.successful() and P > 1e-4*Pc :
            P, m = integral.integrate( integral.t+dr)
            r = integral.t
            rho = (P/K)**(1./gamma)
            T = (P/rho*mp/kB)
            rArray.append(r)
            mArray.append(m)
            eArray.append(T)
            rhoArray.append(rho)
            pArray.append(rho)

            #compute new dr
            dPdr, dMdr = adbDerivs( r, [P,m], K, gamma)
            dr = 0.01*max(min(abs(P/dPdr), abs(m/dMdr)), r)

        rhoc = rhoc * (1.0 - 1.0*(mass-m)/mass)
        print "new rhoc = {0} {1} {4} {2} {3:5.3e}".format(rhoc, rho, m, r, T)

    return mArray, rArray, eArray, rhoArray, pArray

def makeStarMesh(mass, massParticleInGram, posArray, n=3) :
    mModel, rModel, eModel, rhoModel, pModel = makeAdStar( mass, n=n)

    massParticleInCode = massParticleInGram
    numParticles = int(mass/massParticleInGram)
    rArray = np.linalg.norm( posArray, axis=1)
    boolArray = np.argsort(rArray)
    rArray = rArray[boolArray]
    posArray = posArray[boolArray]
    eArray = np.ones(rArray.size)
    vArray = np.ones([rArray.size, 3])

    mArray = np.arange(1,rArray.size)*massParticleInGram
    interpolate = interp1d(mModel,rModel,kind = 'linear')
    interpolateE = interp1d(mModel,eModel,kind = 'linear')

    lastModelP = 0
    
    for i in range(numParticles) :
        if( mArray[i] > mModel[-1] ) :
            break
        lastModelP = i
        r = interpolate( mArray[i])
        eArray[i] = interpolateE( mArray[i])
        posArray[i,:] = posArray[i,:]/rArray[i] * r

    rMax = np.linalg.norm(posArray[lastModelP,:])    
    lastModelP = lastModelP+1
    return posArray[:lastModelP], vArray[:lastModelP], eArray[:lastModelP], massParticleInCode*np.ones(lastModelP)

def findRandV( m1, r1, m2, r2) :
    # compute the hills radius
    r = max(r1+r2, max(r1*(2.*m2/m1)**0.333, r2*(2.*m1/m2)**0.333))
    mtot = m1 + m2
    x1 = m2/mtot * r
    x2 = -m1/mtot * r
    
    v = math.sqrt(G*mtot/r)
    
    pos1 = np.array([x1, 0, 0])
    pos2 = np.array([x2, 0, 0])
    vel1 = np.array([0, v*m2/mtot, 0])
    vel2 = np.array([0, -v*m1/mtot, 0])
   
    tdyn = G**-0.5 
    return pos1, pos2, vel1*tdyn, vel2*tdyn

parser = argparse.ArgumentParser(prog='PROG')
#parser.add_argument('type', default="sedov", help="ic to make (sedov, evrard, star, merger)")
parser.add_argument('--no_atm', action='store_true', help='no atmosphere surrounding object')
args = parser.parse_args()

#pos, vel = makePosVel()
#pos, vel = readSquare()
N = 10
posArray, velArray = makePosVel( N)

# now center
mSun = 1.99e33
mass1 = 1.0*mSun
mass2 = 0.5*mSun
massParticleInGram = mSun/1e6

xmax = 0.5e12
rhoExt = 1e-5
tempExt = 1e5
rho_max = 10

ymax = xmax
zmax = xmax

p1, v1, e1, m1 = makeStarMesh( mass1, massParticleInGram, posArray)
p2, v2, e2, m2 = makeStarMesh( mass2, massParticleInGram, posArray)

print "Number of star particles: {0} + {1} = ".format(e1.size, e2.size, e1.size+e2.size)
 
r1Max = np.linalg.norm(p1[-1])
r2Max = np.linalg.norm(p2[-1])

# find the separation
r1, r2, rdot1, rdot2 = findRandV( mass1, r1Max, mass2, r2Max)
pos  = np.append(p1+r1, p2+r2, axis=0)
vel  = np.append(v1+rdot1, v2+rdot2, axis=0)
ener = np.append(e1, e2)
mass = np.append(m1, m2)

# check xmax
if( 2.*xmax < 2*(r1Max+np.linalg.norm(r1+r2)+r2Max)) : 
    print 'warning'
# do the external gas

if( not args.no_atm) :
    # find the mean separation
    dr1 = np.linalg.norm( p1[0:-2,:] - p1[-1,:], axis=1)
    dr1.sort()
    h1 = dr1[0:16].mean()
 
    dr2 = np.linalg.norm( p2[0:-2,:] - p2[-1,:], axis=1)
    dr2.sort()
    h2 = dr2[0:16].mean()

    h = 0.5*(h1 + h2)

    massParticleExtInGram = 4.*3.1415/3.*rhoExt*h**3
    massParticleExtInCode = massParticleExtInGram

    # find the total mass
    totalExtMass = rhoExt * 8. * xmax**3
    numberExtParticles = totalExtMass/massParticleExtInCode

    #find the scaling
    scaling = (posArray.size/3/numberExtParticles)**0.3333/(0.5*N)
    
    ExtPosArray = posArray*scaling*xmax

    # cut the cube first
    boolArray = ExtPosArray[:,0] < xmax
    boolArray = np.logical_and( boolArray, ExtPosArray[:,0] > -xmax)

    boolArray = np.logical_and( boolArray, ExtPosArray[:,1] < ymax)
    boolArray = np.logical_and( boolArray, ExtPosArray[:,1] > -ymax)

    boolArray = np.logical_and( boolArray, ExtPosArray[:,2] < zmax)
    boolArray = np.logical_and( boolArray, ExtPosArray[:,2] > -zmax)

    ExtPosArray = ExtPosArray[boolArray]
    
    # now exclude the regions around star

    boolArray = np.linalg.norm(ExtPosArray - r1, axis=1) > r1Max
    boolArray = np.logical_and(boolArray, np.linalg.norm(ExtPosArray - r2, axis=1) > r2Max)

    ExtPosArray = ExtPosArray[boolArray]

    # create the other arrays
    NumExtParticles = ExtPosArray.size/3 
    ExtVelArray = np.zeros([NumExtParticles,3])
    ExtMArray = massParticleExtInCode*np.ones(NumExtParticles)
    ExtEArray = tempExt*np.ones(NumExtParticles)

    # now append
    pos = np.append(pos, ExtPosArray, axis=0)
    vel = np.append(vel, ExtVelArray, axis=0)
    ener = np.append(ener, ExtEArray)
    mass = np.append(mass, ExtMArray)
    
x =  pos[:,0]
y =  pos[:,1]
z =  pos[:,2]
vx = vel[:,0]
vy = vel[:,1]
vz = vel[:,2]

print "Number of total particles: {0}".format(mass.size)

xboolArray = np.logical_and( x <= xmax, x > -xmax)
yboolArray = np.logical_and( y <= ymax, y > -ymax)
zboolArray = np.logical_and( z <= zmax, z > -zmax)

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
#pArray = pArray[boolArray]
rho_max = 10
write_tipsy( "merger", x, y, z, vx, vy, vz, mass, temp, rho_max)
#write_gadget( "test.hdf5", x, y, z, vx, vy, vz, ms, pArray/rhos*1.5, rhos, hsoft, xmax, ymax, zmax)
