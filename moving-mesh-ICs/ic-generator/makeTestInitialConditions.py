import numpy as np
from scipy.interpolate import interp1d
try :
    import myhelmholtz as helm
except: 
    print( "cannot import myhelmholtz")

import scipy.integrate as si
import math
import linear_wave as lw
import argparse

from readWriteTipsy import *

kB = 1.38e-16
mp = 1.6e-24
G = 6.67384e-8

SPHERICAL_RESCALING=2
LINEAR_RESCALING=1
rescale = SPHERICAL_RESCALING

def getpressure(dens, temp, abar, zbar):
	pressure,energy,soundspeed,gammaout,entropy = helm.gethelmholtzeos(temp,dens,abar,zbar)
	return pressure

#def getdensity(press, temp, abar, zbar) :
#    dens, energy, soundspeed, gammaout,entropy = helm.getptinversion( temp, abar, zbar, press)
#    return dens

def getderiv(dens, temp, abar, zbar):
	step = dens*1e-4
	pval1 = getpressure(dens + 2*step, temp, abar, zbar)
	pval2 = getpressure(dens + step, temp, abar, zbar)
	pval3 = getpressure(dens - step, temp, abar, zbar)
	pval4 = getpressure(dens - 2*step, temp, abar, zbar)
	output = (-pval1 + 8.0*pval2 - 8.0*pval3 + pval4)/(12.0*step)
	return output

def getdensity(press, temp, abar, zbar):
	muelectron = abar/zbar	#Nucleons per electron
	Kval = 1.0e13		#K for polytrope estimate

	guess = muelectron*(press/Kval)**0.6
	tol = 1e-6
	currenttol = 1.0

	thisguess = guess

	i = 0
	while currenttol > tol:
		nextguess = thisguess - (getpressure(thisguess, temp, abar, zbar)-press)/getderiv(thisguess, temp, abar, zbar)
		currenttol = abs(nextguess - thisguess)/thisguess
		thisguess = nextguess
		i += 1
		if i > 100:
			print( "Error: Newton's Mtd. did not converge!")
			break

	return nextguess

def derivs( r, y, temp, abar, zbar) :
    P = y[0]
    m = y[1]
    rho = getdensity(P, temp, abar, zbar)

    dPdr = -rho*G*m/(r*r)
    dMdr = 4.*math.pi*rho*r*r

    return [dPdr, dMdr]

def makeWD( mass, temp=1e7, abar=12., zbar=6., tol=1e-3) :
    helm.initializehelmholtz()
    # guess a density
    rhoc = 1e6
    #compute the central pressure

    m = 0.
    r = 0.
    Pc = getpressure( rhoc, temp, abar, zbar)
    #print Pc, getdensity(Pc, temp, abar, zbar)
    #return
    mArray = None
    rArray = None
    while abs((m-mass)/mass) > tol :
        integral = si.ode( derivs).set_integrator("dop853")

        Pc = getpressure( rhoc, temp, abar, zbar)
        r = 1e7
        m = 4.*math.pi/3. * rhoc*r**3
        P = Pc - 2.*math.pi/3.*G*r**2*rhoc**2
        y0 = np.array([P, m])
        dr = r*0.1
        integral.set_initial_value(y0, r).set_f_params( temp, abar, zbar)
        mArray = [0.]
        rArray = [0.]
        while integral.successful() and P > 1e-7*Pc :
            #print y
            P, m = integral.integrate( integral.t+dr)
            #y = si.odeint( derivs, y, [r, r+dr], args=(temp, abar, zbar))
            r = integral.t
            rArray.append(r)
            mArray.append(m)

            #compute new dr
            dPdr, dMdr = derivs( r, [P,m], temp, abar, zbar)
            dr = 0.25*min(abs(P/dPdr), abs(m/dMdr))
            #print y, r
        rhoc = rhoc * (1.0 + 0.5*(mass-m)/mass)
        print( "new rhoc = {0:5.3e} {1:5.3e} {2:5.3e}".format(rhoc, m, np.array(rArray).max()))

    #print m, mass, r, rhoc
    return np.array(mArray), np.array(rArray), temp*np.ones(len(mArray))

def adbDerivs( r, y, K, gamma) :
    P = y[0]
    m = y[1]
    rho = (P/K)**(1./gamma)

    dPdr = -rho*G*m/(r*r)
    dMdr = 4.*math.pi*rho*r*r

    return [dPdr, dMdr]

def makeAdStar( mass, gamma=1.667, temp=1e7, tol=1e-3, n=1) :
    rhoc = 1e2
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
            #print y
            P, m = integral.integrate( integral.t+dr)
            #y = si.odeint( derivs, y, [r, r+dr], args=(temp, abar, zbar))
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
            #print y, r
        rhoc = rhoc * (1.0 - 1.0*(mass-m)/mass)
        print( "new rhoc = {0} {1} {4} {2} {3:5.3e}".format(rhoc, rho, m, r, T))

    #print m, mass, r, rhoc
    return mArray, rArray, eArray, rhoArray, pArray
def makeConstantDensity( mass, rmax=None, temp=1e1,sedov=True) :
    if( rmax == None) : 
    	rmax = 8e11
    rArray = np.arange(0.,rmax,rmax*0.001)
    mArray = mass*(rArray/rmax)**3.
    eArray = np.ones(rArray.size)*temp

    if(sedov) :
        rSedov = rmax*0.02

        eArray[ rArray < rSedov] = 1e10
    return mArray, rArray, eArray
    
def makeEvrard( mass, rmax) :
    rArray = np.arange(0.,rmax,rmax*0.001)
    mArray = mass*(rArray/rmax)**2
    eArray = np.ones(rArray.size)*1e3
    
    return mArray, rArray, eArray



parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('type', default="sedov", help="ic to make (sedov, evrard, star, merger)")
parser.add_argument('--no_atm', action='store_true', help='no atmosphere surrounding object')
parser.add_argument('--net_v', action='store_true', help='include netv')
parser.add_argument('--nx', type=int, default=64, help='include netv')
args = parser.parse_args()
#pos, vel = makePosVel()
#pos, vel = readSquare()
N = 8
NX = N
NY = N
NZ = NY
posArray, velArray = makePosVel( NX, NY, NZ, useGlass = False)

# now center
mSun = 1.99e33
massWD = 1.0*mSun
massParticleInGram = mSun/1e5
tempWD = 1e7

xmax = 1e12

rhoExt = 1e-8
tempExt = 1e5
tempGas = None
metalsDefault = 0.
no_atm = args.no_atm

mModel = None
rModel = None
eModel = None
vModel = None
metalsModel = None

r2Shadow = None
a2Shadow = 1.5*1.0
b2Shadow = 1.0**2

scaleHeight = None
gravScaleHeight = None
massModel = massWD

if( args.type == "sedov") : 
    xmax = 1e12
    rhoExt = 1e-8
    tempExt = 1e5
    rho_max = 10

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    massWD = ExtMass
    massModel = ExtMass
    massParticleInGram = ExtMass/1e5
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax)
elif( args.type == "shadow") : 
    xmax = 2.4e19
    rhoExt = 3e-22
    tempExt = 1e2
    rho_max = 1e-21

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    r2Shadow = (rmax*0.2)**2
    massWD = ExtMass
    massModel = ExtMass
    massParticleInGram = ExtMass/1e7
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, sedov=False)
elif ( args.type == "turbBox") :
    xmax = 2.4e19
    rhoExt = 3e-22
    tempExt = 1e2
    rho_max = 1e-21

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    massParticleInGram = ExtMass/2e7
    massWD = ExtMass
    massModel = ExtMass
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, sedov=False)
elif ( args.type == "linear_wave") :
    xmax = 1.5e10
    rhoExt = 1e-10
    tempGas = 100.
    tempExt = tempGas
    rho_max = 3e-18

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    Nparticles = 10000
    massParticleInGram = ExtMass/Nparticles
    ymax = 0.5*xmax#/(1*max(1.,(Nparticles/1e4)**0.333))
    zmax = ymax
    massWD = ExtMass
    massModel = ExtMass
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, temp=tempGas, sedov=False)    
    rescale=LINEAR_RESCALING
elif ( args.type == "thermal") :
    xmax = 1.5e18
    rhoExt = 1e-15
    tempExt = 1e2
    rho_max = 1e-14

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    massParticleInGram = ExtMass/1e5
    massWD = ExtMass
    massModel = ExtMass
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, sedov=False)
elif ( args.type == "davis") : 
    Tstar = 82.
    tempGas = Tstar
    tempExt = Tstar
    mu = 2.33
    grav = 1.46e-6
    cs = math.sqrt(kB*Tstar/(mu*mp))
    hstar = (cs*cs/grav)
    kappaRstar = 0.0316*(Tstar/10)**2
    taustar = 10.
    Sigmastar = taustar/kappaRstar
    rhostar = Sigmastar/hstar
    ymax = 100*hstar
    xmax = ymax
    zmax = ymax

    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    xmax = ymax/2
    zmax = xmax

    rhoExt = 1e-10*rhostar
    rho_max = rhostar
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    massParticleInGram = ExtMass/1e7
    massWD = ExtMass
    massModel = ExtMass
    scaleHeight = hstar
    gravScaleHeight = 2*hstar

    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, temp=tempGas, sedov=False)

    print( "xmax, ymax, zmax= ", xmax, ymax, zmax)
    print( "dxmax, dymax, dzmax= ", 2*xmax, 2*ymax, 2*zmax)
    print( "hstar/cs", hstar/cs/G**-0.5, hstar/cs, cs)
    print( "ExtMass = ", ExtMass)
    print( "hstar = ", hstar)
    print( "rhostar, rhomin = ", rhostar, rhoExt)

elif ( args.type == "radShock") :
    xmax = 1.5e18
    rhoExt = 1e-15
    tempExt = 1e2
    rho_max = 1e-14

    ymax = 1.5e18
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    massParticleInGram = ExtMass/1e6
    massWD = ExtMass
    massModel = ExtMass
    mModel, rModel, eModel = makeConstantDensity( ExtMass, rmax, sedov=False)
elif ( args.type == "wd") :
    xmax = 2.5e9
    tempWD = 1e7
    massWD = 0.6*mSun
    massModel = massWD
    rhoExt = 1e-2
    tempExt = 1e5
    rho_max = 3e6

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    mModel, rModel, eModel = makeWD( massWD, temp=tempWD)

elif ( args.type == "star") :
    xmax = 0.5e12
    rhoExt = 1e-5
    tempExt = 1e5
    rho_max = 10

    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    mModel, rModel, eModel, rhoModel, pModel = makeAdStar( massWD, n=3)

elif ( args.type == "evrard") :
    xmax = 1e12
    rhoExt = 1e-8
    tempExt = 1e5
    rho_max = 10
 
    ymax = xmax
    zmax = xmax
    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3

    massModel = 1e3*mSun
    massParticleInGram = massModel/1e6
    mModel, rModel, eModel = makeEvrard( massModel, rmax)

massParticleInCode = massParticleInGram
rArray = np.linalg.norm( posArray, axis=1)
boolArray = np.argsort(rArray)
rArray = rArray[boolArray]
posArray = posArray[boolArray]
eArray = np.ones(rArray.size)
metalsArray = np.ones(rArray.size) *metalsDefault
#rhoArray = np.ones(rArray.size)
#pArray = np.ones(rArray.size)
mArray = np.arange(1,rArray.size)*massParticleInGram
interpolate = interp1d(mModel,rModel,kind = 'linear')
#interpolateRho = interp1d(mModel,rhoModel,kind = 'linear')
#interpolateP = interp1d(mModel,pModel,kind = 'linear')
interpolateE = interp1d(mModel,eModel,kind = 'linear')
interpolateV = None
if( vModel != None) :
    interpolateV = interp1d(mModel,vModel,kind = 'linear')

interpolateMetals = None
if( metalsModel != None) :
    interpolateMetals = interp1d(mModel,metalsModel,kind = 'linear')

lastModelP = 0
# for i in range(numParticlesWD) :

#     if( mArray[i] > mModel[-1] ) :
#         break
#     lastModelP = i
#     #print mModel[0]
#     #print mArray[i]/mModel[-1]
#     r = interpolate( mArray[i])
#     #rhoArray[i] = interpolateRho( mArray[i])
#     #pArray[i] = interpolateP( mArray[i])
#     eArray[i] = interpolateE( mArray[i])
#     posArray[i,:] = posArray[i,:]/rArray[i] * r
    
#     if( vModel != None) :
#         v = interpolateV( mArray[i])
#         velArray[i,:] = posArray[i,:]/r*v

#     if( metalsModel != None) :
#         metalsArray[i] = interpolateMetals( mArray[i])

#print mModel[0]
#print mArray[i]/mModel[-1]
numParticlesWD = None
rModelMax = None
if( rescale == SPHERICAL_RESCALING) :   

    numParticlesWD = int(mModel[-1]/massParticleInGram)
    print( "Building stellar particles with {0} particles".format( numParticlesWD))
    print( mArray.size)

    r = interpolate( mArray[0:numParticlesWD])
    eArray[0:numParticlesWD] = interpolateE( mArray[0:numParticlesWD])
    posArray[0:numParticlesWD,0] = posArray[0:numParticlesWD,0]/rArray[0:numParticlesWD] * r[0:numParticlesWD]
    posArray[0:numParticlesWD,1] = posArray[0:numParticlesWD,1]/rArray[0:numParticlesWD] * r[0:numParticlesWD]
    posArray[0:numParticlesWD,2] = posArray[0:numParticlesWD,2]/rArray[0:numParticlesWD] * r[0:numParticlesWD]
   
    velArray[:,:] = 0 
    if( vModel != None) :
        v = interpolateV( mArray[0:numParticlesWD])
        velArray[0:numParticlesWD,:] = posArray[0:numParticlesWD,:]/r[:,np.newaxis]*v[:,np.newaxis]

    if( metalsModel != None) :
        metalsArray[0:numParticlesWD] = interpolateMetals( mArray[0:numParticlesWD])

    rModelMax = np.linalg.norm(posArray[numParticlesWD-1,:])
else :
    numberParticles = massModel/massParticleInCode
    # find the aspect ratio
    sizes = np.array([xmax,ymax,zmax])
    sorted = sizes.argsort()
    aspectRatio = sizes[0]*sizes[1]*sizes[2]/sizes.max()**3
    currentNumber = posArray.size/3*aspectRatio

    linearScaling = (numberParticles/currentNumber)**0.3333
    print( "aspect ratio {0} linear rescale {1} {2}".format(aspectRatio, linearScaling, 1./linearScaling*sizes.max()/(0.5*N)))
    posArray *= 1./linearScaling*sizes.max()/(0.5*N)
    lastP = posArray.size/3
    numParticlesWD = int(lastP)
    print( "number particle {0}".format(numberParticles))
       
lastModelP = numParticlesWD-1
# do the external gas
lastP = lastModelP

if( not no_atm and rModelMax < rmax) :
    # find the mean separation

    dr = np.linalg.norm(posArray[0:lastModelP,:] - posArray[lastModelP,:], axis=1)
    dr.sort()
    h = dr[0:8].mean()
    massParticleExtInGram = 4.*3.1415/3.*rhoExt*h**3
    massParticleExtInCode = massParticleExtInGram
    print( h, massParticleExtInGram)

    #subtract the filled in mass
    rModel = np.linalg.norm( posArray[lastModelP,:])
    mExtArray = np.arange(0,rArray.size-lastModelP)*massParticleExtInGram
    rExtModel = np.arange(rModel,rmax,rModel*0.01)
    mExtModel = rhoExt*4.*3.1415/3.*(rExtModel**3-rModel**3)
    #print mExtModel
    interpolateExt = interp1d(mExtModel,rExtModel,kind = 'linear')
    totNumParticles = posArray.shape[0]
    for i in range(0,totNumParticles-lastModelP) :
        if( mExtArray[i] > mExtModel[-1]) :
            break
        #print "building {0}".format(mExtArray[i])
        lastP = i+lastModelP+1
        r = interpolateExt( mExtArray[i])
        posArray[lastP,:] = posArray[lastP,:]/rArray[lastP] * r
    print( "lastP = {0} lastModelP={1}".format(lastP,lastModelP), mExtModel[-1]/massParticleExtInGram)

x =  posArray[0:lastP,0]
y =  posArray[0:lastP,1]
z =  posArray[0:lastP,2]
vx = velArray[0:lastP,0]
vy = velArray[0:lastP,1]
vz = velArray[0:lastP,2]

ms = massParticleInCode*np.ones(lastP)
if( args.type == "shadow") :
    ms[x*x/a2Shadow + y*y/b2Shadow < r2Shadow]= 1000.*massParticleInCode

if( args.type == "radShock") : 
    boolArray = y < 0 
    vshock = 1e6
    vy[:] = -vshock
    vy[y < 0.] = vshock


print( ms.max(), ms.min())
temp = eArray[0:lastP]
metals = metalsArray[0:lastP]

if( lastModelP < lastP) :
    ms[lastModelP+1:lastP] = massParticleExtInCode
    temp[lastModelP+1:lastP] = tempExt
    #pArray[lastModelP+1:lastP] = kB*tempExt*rhoExt/mp
    #rhoArray[lastModelP+1:lastP] = rhoExt
    print( "Number of WD particles: {0}".format(ms[ms > massParticleExtInCode].size))
else :
    print( "Number of WD 2-particles: {0}".format(ms.size))

#print temp[0:lastP]
    # clean up for fit in box

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
if(args.net_v) :
    vx = vx + 3e5/G**0.5

ms = ms[boolArray]
#rhos = rhoArray[boolArray]
temp = temp[boolArray]
metals = metals[boolArray]
#pArray = pArray[boolArray]

if( args.type == "linear_wave") : #setup linear wave in x-direction
    Tgas = 1e7 #tempGas
    rho = 1e-2
    mu = 1
    mp = 1.6e-24
    kB = 1.38e-16
    Pgas = kB*Tgas/(mu*mp)*rho
    gamma = 1.666667
    cs = math.sqrt(Pgas/rho)
    arad = 7.56e-15
    Trad = Tgas
    Prad = arad*Trad**4
    Prat = Prad/Pgas
    cspeed = 3e10
    cred = 1e-2*cspeed 
    crat = cspeed/cs
    kappa = 3.33333e-7
    r = cred/cspeed
    #r = 1.
    L = 2*xmax
    wavelength=L
    ka = 2*math.pi/(wavelength/L)
    #sigma = cspeed*kappa*rho*(L/cs)
    sigma = cred*kappa*rho*(L/cs)
    sigma = kappa*rho*wavelength
    
    print( "ka={0:5.3e} sigma={1:5.3e} Prat = {2:5.3e} crat = {3:5.3e} r={4:5.3e} {5:5.3e}".format(ka, sigma, Prat, crat, r, crat*r/sigma))
    print( "rho={0:5.3e} cs={1:5.3e} cr={2:5.3e} kappa={3:5.3e} Trad={4:5.3e}".format( rho, cs, cred, kappa, Trad))
    omega, drho, dv, dT, dp, dEr, dFr = lw.eigenvec(ka=ka,crat=crat,prat=Prat,sigma=sigma,gamma=gamma,r=r)
    delta = 1e-1
    kr = ka*x/L
    coskx = np.cos(kr)
    sinkx = np.sin(kr)
    tdyn = G**-0.5
    t0 = L/cs/tdyn
    print( "time to integrate {0:7.4e} {1:7.4e}".format(2*math.pi/omega.real*t0, 2*math.pi/omega.real*t0/tdyn))
    print( "dxperiod {0:12.7e} dyperiod {1:12.7e} dzperiod {2:12.7e}".format(2*xmax, 2*ymax, 2*zmax))
    #ms += (coskx*drho.real + sinkx*drho.imag)*delta*ms
    #temp += (coskx*dT.real + sinkx*dT.imag)*delta*temp
    #vx += (coskx*dv.real + sinkx*dv.imag)*delta*cs*tdyn

    # redo everything
    nx = args.nx
    dx = 2*xmax/nx
    xf, yf, zf = None, None, None
    frac = 0.5*1e-6
    for zpos in np.arange(-dx, dx*1.5, dx) :
        for ypos in np.arange(-dx, dx*1.5, dx) :
            x = np.arange(-xmax+dx*0.5, xmax, dx)
            x += 2*(np.random.random(x.size)-0.5)*dx*frac
            y = ypos*np.ones(x.size) + 2*(np.random.random(x.size)-0.5)*dx*frac
            z = zpos*np.ones(x.size) + 2*(np.random.random(x.size)-0.5)*dx*frac
            if xf is None : 
                xf = x.copy()
                yf = y.copy()
                zf = z.copy()
            else : 
                xf = np.append( xf, x)
                yf = np.append( yf, y)
                zf = np.append( zf, z)

    x = xf
    y = yf
    z = zf
    vx = np.zeros(xf.size)
    vy = np.zeros(xf.size)
    vz = np.zeros(xf.size)
    ms = np.ones(xf.size)
    temp = np.zeros(xf.size)
    metals = None
    print( "Period = {0:5.3e}".format(2*dx*1.5))
    rho_max = 1e-10
    
if( args.type == "davis") : 
    rho = rhoExt
    meanHeight2 = scaleHeight*gravScaleHeight
    norm = 1. + 0.25*np.sin(2.*math.pi*x/xmax)*np.sin(2.*math.pi*z/zmax)
    ms *= norm*np.maximum(np.ones(ms.size), rho_max/rho*np.maximum(math.exp(-gravScaleHeight/scaleHeight)*np.exp(-np.abs(y)/scaleHeight), np.exp(-y*y/meanHeight2)))
    #y = np.arange(-ymax, ymax, 0.01*ymax)
    #print np.maximum(np.ones(y.size), rho_max/rho*np.maximum(math.exp(-gravScaleHeight/scaleHeight)*np.exp(-np.abs(y)/scaleHeight), np.exp(-y*y/meanHeight2)))
    #quit()

write_tipsy( args.type, x, y, z, vx, vy, vz, ms, temp, rho_max, gasmetals = metals)
