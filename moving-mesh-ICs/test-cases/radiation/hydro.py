import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import math
from numba import jit,njit, prange
import linear_wave as lw
import scipy.optimize as opt
import argparse
import sys

noLorentz = False
#Geometry
CARTESIAN = 0
SPHERICAL = 2
CYLINDRICAL = 1
RIEMANN_TEST = 1

useRadiation = True
geometry = CARTESIAN

delta = 1e-3

SOUND_WAVE = 0
RADIATION_PROP = 1
THERMAL_EQ = 2
LINEAR_WAVE = 3
PROBLEM = LINEAR_WAVE

if( PROBLEM == SOUND_WAVE) : 
    useRadiation = False


#constants
pi = math.pi
kB = 1.3806e-16
mp = 1.6726e-24
kms = 1e5
FOUR_THIRD_PI_G = 4./3.*math.pi*6.6743e-8
xmax = 1.5e10
cspeed = 2.9979e10
fourPi = 4*pi
aRad = 7.5657e-15
nradx = 1./np.sqrt(3.)


NVARS=3
if( useRadiation) :
    NVARS=5

IRHO = 0
IMOM = 1
IENER = 2
IRADL = 3
IRADR = 4

#run parameters
courant = 0.8
gamma = 5./3.
NGRID= 130
rho0 = 10000.*mp
T0 = 1e2
Trad0 = 1e1
Tmin = 1e0
rhoMin = 0.1*mp
Tmax = 1e6
cs0 = math.sqrt(gamma*kB*T0/mp)
cred = 2.*cs0
kappa = 1e11

#T0 = .485e7/2.15 #tempGas
T0 = .485e7 #tempGas
T0 = 10.85e10 #tempGas
#T0 = 10*1.46**2*0.485e10
#T0 = 0.2252e7 #tempGas
#T0 = 2.15*1.045e7 #tempGas

#rho0 = 0.9685e11
rho0 = 1e-2
rho0 = 20*5.88e9
cs0 = math.sqrt(kB*T0/mp)
Trad0 = T0
Prat = aRad*Trad0**4/(rho0*kB*T0/mp)
cred = 1e0*cspeed 
kappa = 3.3333e-9
L = 2*xmax
wavelength=L
tcross = L/cs0

writePlots = True

def findEigen( quiet=False) : 
    ka = 2*math.pi/(wavelength/L)
    crat = cspeed/cs0
    r = cred/cspeed
    sigma = kappa*rho0*wavelength
    if not quiet :
        print("crat = {0:5.3e} ka = {1:5.3e} prat = {2:5.3e} sigma = {3:5.3e} r = {4:5.3e}".format(crat, ka, Prat, sigma, r))
        print("kappa={0:5.3e} rho={1:5.3e} cs={2:5.3e}".format(kappa, rho0, cspeed/crat))
    quiet=False
    omega, drho, dv, dT, dp, dEr, dFr = lw.eigenvec(ka=ka,crat=crat,prat=Prat,sigma=sigma,gamma=gamma,r=r,quiet=quiet)

    if not quiet :
        print( "dWaveRhoReal = {0:9.6e}\ndWaveRhoImag = {1:9.6e}".format( drho.real, drho.imag))
        print( "dWaveVxReal = {0:9.6e}\ndWaveVxImag = {1:9.6e}".format( dv.real, dv.imag))
        print( "dWaveIeReal = {0:9.6e}\ndWaveIeImag = {1:9.6e}".format( dT.real, dT.imag))
        print( "dWaveErReal = {0:9.6e}\ndWaveErImag = {1:9.6e}".format( dEr.real, dEr.imag))
        print( "dWaveFrReal = {0:9.6e}\ndWaveFrImag = {1:9.6e}".format( dFr.real, dFr.imag))
    return omega, drho, dv, dT, dp, dEr, dFr


@njit
def geoFactor(x):
#
# Geometry factors 
#
    return 1 # only cartesian
    if(geometry == CARTESIAN) :
        return 1.
    elif(geometry == SPHERICAL) : 
        return x*x
    elif(geometry == CYLINDRICAL) : 
        return x

def generateGrid( rmin, rmax, rhoInit=rho0, vInit = 0., Tinit=T0, Trad=Trad0, quiet=False) :
# Generates the grid and sets initial conditions on primitives
    dx = np.ones(NGRID)*(rmax-rmin)/(NGRID-2)
    xCoord = dx*(np.arange(NGRID)+0.5 - 1)+rmin

    rho = rhoInit*np.ones(NGRID)
    v   = vInit*np.ones(NGRID) 
    ie  = kB*Tinit/mp/(gamma-1.)*np.ones(NGRID)
    radL = aRad/fourPi*Trad**4*np.ones(NGRID)
    radR = aRad/fourPi*Trad**4*np.ones(NGRID)

    if( PROBLEM == SOUND_WAVE) :
        rho += delta*rhoInit*np.cos( math.pi*xCoord/rmax)
        v    = delta*math.sqrt(gamma)*cs0*np.cos( math.pi*xCoord/rmax)
        ie  += delta*kB*Tinit/mp*np.cos( math.pi*xCoord/rmax)
    elif( PROBLEM == RADIATION_PROP) :
        radL = np.zeros(NGRID)#aRad*Tinit**4*np.ones(NGRID)
        radR = np.zeros(NGRID)
        radR = aRad*Tinit**4*np.exp(-xCoord**2/(0.09*xmax**2))
    elif( PROBLEM == THERMAL_EQ) :
        radR = aRad/fourPi*Trad**4*np.ones(NGRID)
        radL = aRad/fourPi*Trad**4*np.ones(NGRID)
    elif( PROBLEM == LINEAR_WAVE) : 
        omega, drho, dv, dT, dp, dEr, dFr = findEigen(quiet=quiet)
        x = xCoord
        rho += rhoInit*delta*(drho.real*np.cos(pi*x/rmax) + drho.imag*np.sin(pi*x/rmax))
        v = cs0*delta*(dv.real*np.cos(pi*x/rmax) + dv.imag*np.sin(pi*x/rmax))
        ie += ie*delta*(dT.real*np.cos(pi*x/rmax) + dT.imag*np.sin(pi*x/rmax))
        deltaE = aRad*Trad**4/fourPi*delta*(dEr.real*np.cos(pi*x/rmax) + dEr.imag*np.sin(pi*x/rmax))
        deltaF = aRad*Trad**4/fourPi*delta*(dFr.real*np.cos(pi*x/rmax) + dFr.imag*np.sin(pi*x/rmax))
        radR = aRad/fourPi*Trad**4*np.ones(NGRID) + deltaE + deltaF/(2*nradx)
        radL = aRad/fourPi*Trad**4*np.ones(NGRID) + deltaE - deltaF/(2*nradx)
    return dx, xCoord, Prim2ConGrid( xCoord, rho, v, ie, radL, radR)


def Prim2ConGrid(xcoord, rho, v, ie, radL = None, radR = None) :
#
#
    grid = np.zeros([NGRID,NVARS])
    for i in range(xcoord.size) :
        if( radL is None) :
            grid[i,:] =  Prim2Con(xcoord[i],rho[i],v[i],ie[i])
        else :
            grid[i,:] =  Prim2Con(xcoord[i],rho[i],v[i],ie[i], radL[i], radR[i])

    return grid

@njit
def Prim2Con( r, rho, v, ie, radL = 0., radR = 0.) : 
#
# Primitive to conservative solver
#
    u = np.zeros(NVARS)
    gFactor = geoFactor(r)
    u[IRHO] = rho*gFactor
    u[IMOM] = rho*v*gFactor
    u[IENER] = rho*(ie+0.5*v*v)*gFactor
    if( useRadiation) :
        u[IRADL] = radL
        u[IRADR] = radR
    return u

@njit
def Con2PrimGrid(xcoord, grid) :
#
# Conservative to Primitive Solver on the grid
#
    rho = []
    v = []
    ie = []
    radL = []
    radR = []
    for i in range(grid.shape[0]) : 
        dens, vel, iener, IradL, IradR = Con2Prim(xcoord[i], grid[i,:])
        rho.append(dens)
        v.append(vel)
        ie.append(iener)
        radL.append( IradL)
        radR.append( IradR)
    return np.array(rho), np.array(v), np.array(ie), np.array(radL), np.array(radR)

@njit
def Con2Prim(r, u) :
#
# Conservative to primite solver
#
    gFactor = geoFactor(r)
    rho = u[IRHO]/gFactor
    v = u[IMOM]/(rho*gFactor)
    ie = (u[IENER]/rho - 0.5*v*v)/gFactor
    radL = 0
    radR = 0
    if( useRadiation) :    
        radL = u[IRADL]
        radR = u[IRADR]

    return rho, v, ie, radL, radR

@njit
def HLLE( r, csL, csR, FL, FR, UL, UR) :
#
# HLLE Riemann solver
#
    vL = FL[IRHO]/UL[IRHO]
    vR = FR[IRHO]/UR[IRHO]
    sL = min(vL, vR) - max(csL, csR)
    sR = max(vL, vR) + max(csL, csR)

    if( sL >= 0.) :
        return FL
    elif ( sR <= 0.):
        return FR
    elif ((sL <= 0.) and (sR >= 0.)) :
        FHLL = (FL * sR - FR * sL + (UR - UL) * sL * sR) / (sR - sL)
        return FHLL

    print( sL, sR, csL, csR)

    return FL

@njit
def alpha( rho, T) : 
    kappaPlanck, kappaMean, kappaSca = computeOpacity( rho, T)
    tau = (L/NGRID*10.*(kappaMean + kappaSca)*rho)**2
    factor = np.sqrt((1.-math.exp(-tau))/tau)
    return 1
    #print(tau, factor)
    if tau < 1e-2 : 
        return 1

    return factor


@njit 
def upwindRad( uL, uR) : 
#calculate the upwind flux for rad
    rho = uL[IRHO]
    T = T0

    alphaC = cred*alpha(rho, T)
    fluxRad = np.zeros(uL.shape)
    fluxRad[IRADL] = -uR[IRADL] * nradx * alphaC
    fluxRad[IRADR] = uL[IRADR] * nradx * alphaC
    return fluxRad

@njit
def SoundSpeed( ie) :
#
# Get Sound speed for gamma equation of state
#  
    cs = np.sqrt((gamma-1.)*ie*gamma)
    return cs

@njit
def TimeStep( dr, rho, v, ie) :
#
# return the minimum timestep
#
    cs = SoundSpeed( ie)
    alphaC = alpha( rho0, T0)*cred*nradx
    return courant*np.min( dr/(np.abs(v) + cs + alphaC))

@njit
def calFlux( r, uL, uR) : 
#
# Calculate the flux in based in reconstructed left and right states
#
    rhoL, vL, iEL, radL, radR = Con2Prim(r, uL)
    rhoR, vR, iER, radL, radR = Con2Prim(r, uR)
    FL = np.zeros(NVARS)
    FR = np.zeros(NVARS)
    gFactor = geoFactor(r)
    FL[IRHO] = uL[IMOM]
    FR[IRHO] = uR[IMOM]
    csL = SoundSpeed( iEL)
    PL = rhoL*csL*csL/gamma
    FL[IMOM] = vL*uL[IMOM] + PL*gFactor
    FL[IENER] = (uL[IENER] + PL*gFactor)*vL
    csR = SoundSpeed( iER)
    PR = rhoR*csR*csR/gamma
    FR[IMOM] = vR*uR[IMOM] + PR*gFactor
    FR[IENER] = (uR[IENER] + PR*gFactor)*vR
    return FL, csL, FR, csR



@njit
def reconstruct( dx, xcoord, grid, firstOrder=False) :
#
# reconstruction -- piecewise constant for now
#
    TINY = 1e-20
    rho, v, ie, radL, radR = Con2PrimGrid(xcoord, grid)
    uL = np.zeros(grid.shape)
    uR = np.zeros(grid.shape)
    uL[1:-1] = grid[1:-1]
    uR[1:-1] = grid[1:-1]
    if( not firstOrder) : 
        for i in range(1, xcoord.size-1) :
            dudx = (grid[i+1,:] - grid[i-1,:])/(2.*dx[i])
            phi = limiter( grid[i-1,:], grid[i,:], grid[i+1,:])
            uL[i,:] -= dudx[:]*dx[i]*0.5*phi
            uR[i,:] += dudx[:]*dx[i]*0.5*phi
    uL = bc(xcoord, uL)
    uR = bc(xcoord, uR)

    return uL, uR

@njit
def limiter( uminus, u, uplus):
#minmod limiter
    #print(r)
    TINY = 1e-30
    r = (u - uminus)/(uplus - u)
    r[ np.abs(uplus-u)<TINY] = 1.
    #return 0.
    #return np.maximum( 0., np.minimum(1., r)).min() #minmod
    return (r+np.abs(r))/(1.+np.abs(r)) #van Leer
    #beta = 1.5
    #return np.maximum(0., np.minimum(r, beta))

@njit
def bc(xcoord,grid) :
#
# impose periodic boundary conditions
#
    grid[-1,:] = grid[1,:]
    grid[0,:] = grid[-2,:]

    return grid

@njit
def validateGrid( r, grid) : 
    rhoGrid, vGrid, ieGrid, radLGrid, radRGrid = Con2PrimGrid( r, grid) 
    for i in range(NGRID) :
        resetGrid = False
        if(rhoGrid[i] < rhoMin) :
            resetGrid = True
            rhoGrid[i] = rhoMin
        if( ieGrid[i] < kB*Tmin/mp/(gamma-1)) :
            resetGrid = True
            ieGrid[i] = kB*Tmin/mp/(gamma-1)
        if( radLGrid[i] < 0.) :
            resetGrid = True
            radLGrid[i] = 0.
        if( radRGrid[i] < 0.) : 
            resetGrid = True
            radRGrid[i] = 0.
        if( resetGrid) : 
            grid[i,:] = Prim2Con(r[i], rhoGrid[i],vGrid[i],ieGrid[i],radLGrid[i],radRGrid[i])

@njit
def AdvanceTimeStep( t, dt, dr, r, grid) :
#
# make a timestep -- second order
#
    grid1 = grid.copy()
    
    for iOrder in range(2) : 
        grid0 = grid1
        grid1 = grid.copy()
        dt1 = dt
        firstOrder = False

        if( iOrder == 0) : # first half step
            dt1 = dt*0.5
            firstOrder = True        
        uL, uR = reconstruct( dr, r, grid0, firstOrder=firstOrder)
        rhoGrid, vGrid, ieGrid, radLGrid, radRGrid = Con2PrimGrid( r, grid0) 

        for i in range(1,grid0.shape[0]-1) :
            x = r[i]
            dx = dr[i]
            ul = uR[i]
            ur = uL[i+1]
            fl, csl, fr, csr = calFlux(x+dx, ul, ur)
            fluxR = HLLE( x, csl, csr, fl, fr, ul, ur) 
            if( useRadiation) :
                temp = upwindRad( ul, ur)
                fluxR[IRADL] = temp[IRADL]
                fluxR[IRADR] = temp[IRADR]

            ul = uR[i-1]
            ur = uL[i]
            fl, csl, fr, csr = calFlux(x-dx, ul, ur)
            fluxL = HLLE( x, csl, csr, fl, fr, ul, ur)
            if( useRadiation) :
                temp = upwindRad( ul, ur)

                fluxL[IRADL] = temp[IRADL]
                fluxL[IRADR] = temp[IRADR]
            gFactor = geoFactor(x)

            rho, v, ie, radL, radR = Con2Prim(x, grid0[i,:])
        
            grid1[i] -= (fluxR-fluxL)/dr[i]*dt1 
            #grid1[i,IMOM] += geometry*gFactor/x*P*dt1
            #grid1[i] += sources(t+dt1-dt*0.5, x, grid0[i,:])*dt1
            if useRadiation:
                grid1[i,:] = RadiationSource( dt1, x, grid1[i,:])
            bc(r, grid1)
        validateGrid(r, grid1)


    
    gridNew = grid1
    return gridNew


@njit
def sources( t, r, u) :
#
# Impose sources
#
    source = np.zeros(NVARS)
    return source

@njit
def FourthPolyRoot(coef4, tconst) :
# First, get the double root of
# z^3-4*tconst/coef4 * z - 1/coef4^2==0
    asquar = coef4 * coef4
    acubic = coef4 * asquar
    ccubic = tconst * tconst * tconst
    delta1 = 0.25 - 64.0 * ccubic * coef4/27.0
    if(delta1 < 0.0) :
        return None
    else :
        delta1 = math.sqrt(delta1)

    zroot = 0.0
    if(delta1 > 1.e11) : 
    #to avoid small number cancellation
        zroot = delta1**(-2.0/3.0)/3.0
    else:
        zroot = (0.5 + delta1)**(1.0/3.0) - (-0.5 + delta1)**(1.0/3.0)
    zroot *= coef4**(-2.0/3.0)
    rcoef = math.sqrt(zroot)
    delta2 = -zroot + 2.0/(coef4*rcoef)

    if(delta2 < 0.0) : 
        return None
    else :
        delta2 = math.sqrt(delta2)

    root = 0.5 * (delta2 - rcoef)

    if(root < 0.0) :
        return None
    return root

@njit
def computeOpacity( rho, Tgas) : 
    global kappa
    kappaPlanck = kappa
    kappaMean = kappa
    kappaSca = 0.
    return kappaPlanck, kappaMean, kappaSca

@njit
def JRad( radL, radR) : 
    return 0.5*(radL + radR)
    
@njit
def computeNewTgas( dt, rho, v, ie, radL, radR) :
# compute the old J from intensities
    Told = ie*(gamma-1.)*mp/kB
    Tgas = Told
    RTerm = rho*kB/mp/(gamma-1.)
    beta = 0.
    lorentz2 = 1.

    J = 0. 
    Jprime = 0.
    invJfactorT4 = 0.

    kappaPlanck = 0.
    kappaMean = 0.
    kappaSca = 0.

    kappaPlanck, kappaMean, kappaSca = computeOpacity( rho, Tgas)

    sigmaAbsPl =   rho*kappaPlanck
    sigmaAbsMean = rho*kappaMean
    sigmaSca =     rho*kappaSca

    rads = np.array([radL, radR])
    nxs = np.array([-nradx, nradx])
    J = JRad( radL, radR)
    #  double vnc = 1. - dotProduct(beta, HydroUtils::instance->radNormals[i]);
    #  double Gamma = vnc*sqrt(lorentz2);
    #get angular dependent gamma 
    beta = v/cspeed
    if( noLorentz) :
        beta = 0.

    lorentz2 = 1./(1.0 - beta*beta)
    Jprime = 0
    invJfactorT4 = 0.
    SigmaI = 0.
    SigmaGamma = 0.
    for i in range(2) : 
        nx = nxs[i]
        rad = rads[i]
        vnc = 1 - beta*nx
        Gamma = vnc/math.sqrt(lorentz2)
        sigmaAbsMeanGamma = sigmaAbsMean*Gamma
        csigmaGammadt = cred*sigmaAbsMeanGamma*dt
        Jfactor = 1. + csigmaGammadt
        SigmaI += 0.5*rad/Jfactor
        SigmaGamma += 0.5*Gamma/Jfactor
#        Jprime += rad/Jfactor*0.5
#        invJfactorT4 += Gamma/Jfactor*0.5
    
    # fourth order root solve
    credsigmadt = cred*sigmaAbsMean*dt
    csigmadt = cspeed*sigmaAbsMean*dt#/math.sqrt(3.)
    csigmaPldt = cspeed*sigmaAbsPl*dt
    eta = csigmadt/(1.+credsigmadt)/RTerm
#    coef4 = csigmaPldt/RTerm*aRad*(1. - credsigmadt*invJfactorT4)
#    tconst = -csigmadt/RTerm*Jprime*fourPi - Told
    coef4 = csigmaPldt/RTerm*aRad*(1. - credsigmadt*SigmaGamma)
    tconst = -csigmadt/RTerm*SigmaI*fourPi - Told

    #coef4 = eta*aRad
    #tconst = -eta*J*fourPi - Told
    Tnew = FourthPolyRoot( coef4, tconst)
#    print( (Tnew-Told)/Told/(dt/tcross))
#    print(Tnew, coef4, tconst, J, Told)
    if( Tnew is None) :
        Tnew = Told
        deltaE = 0.

    T4 = Tnew*Tnew*Tnew*Tnew
    Jold = J
    J = credsigmadt*aRad/fourPi*T4*SigmaGamma + SigmaI
    deltaE = RTerm*(Tnew - Told)

    #print(deltaE, (J-Jold)*fourPi )

    return J, Tnew, deltaE

@njit
def RadiationSource( dt, r, uNew) :
    uOldLab = uNew.copy()
    uOld = uNew.copy()
    #return uNew
    uComOld = labToCom( uOld)
    uComNew = uComOld.copy()
    rho, v, ie, radL, radR = Con2Prim( r, uComOld) 
  
    rCSpeed = cred
    cSpeed = cspeed

    # implicit update
    T = ie*(gamma-1)*mp/kB
    J, Tnew, DeltaE = computeNewTgas( dt, rho, v, ie, radL, radR)
    kappaPlanck, kappaMean, kappaSca = computeOpacity( rho, T)
    sigmaAbsPl =   rho*kappaPlanck
    sigmaAbsMean = rho*kappaMean
    sigmaSca =     rho*kappaSca
       
    rads = np.array([radL, radR])
    nxs = np.array([-nradx, nradx])
    newRads = rads
    beta = v/cspeed
    if( noLorentz) :
        beta = 0.

    lorentz2 = 1./(1.0 - beta*beta)
    
    for i in range(2) : 
        nx = nxs[i]
        rad = rads[i]
        vnc = 1. - beta*nx
        Gamma = vnc*math.sqrt(lorentz2)

        cSigmaAbsMeandt = dt*rCSpeed*sigmaAbsMean*Gamma
        cSigmaAbsPldt = dt*rCSpeed*sigmaAbsPl*Gamma
        cSigmaScatdt = dt*rCSpeed*sigmaSca*Gamma
        invJfactor = 1./(1. + cSigmaAbsMeandt + cSigmaScatdt)
      
        deltaRad = cSigmaAbsPldt*aRad*Tnew**4/fourPi# + cSigmaScatdt*J
        newRads[i] = (rad + deltaRad)*invJfactor

    
    uComNew[IRADL] = newRads[0]
    uComNew[IRADR] = newRads[1]

    uNewLab = comToLab(uComNew)
    flux = (uNewLab[IRADR] - uNewLab[IRADL])*nradx*0.5
    oldFlux = (uOldLab[IRADR] - uOldLab[IRADL])*0.5*nradx
    #include changes in matter source terms.
    uNewLab[IENER] += DeltaE

    DeltaP = (oldFlux-flux)/cred*fourPi
    uNewLab[IMOM] += DeltaP
    DeltaKE = 0.5*(uNewLab[IMOM]**2 - uOldLab[IMOM]**2)/rho
    uNewLab[IENER] += DeltaKE

    dUdt = (uNewLab - uOldLab)/(dt/tcross)/uOldLab
    #print(r/L, dUdt)
    return uNewLab

@njit
def transCoeff( v) : 
    beta = v/cspeed
    lorentz2 = 1./(1.0 - beta*beta)
    vncL = 1. + beta*nradx
    vncR = 1. - beta*nradx
    angCoeffL = vncL*vncL*lorentz2
    angCoeffR = vncR*vncR*lorentz2
    return angCoeffL*angCoeffL, angCoeffR*angCoeffR
    
@njit
def labToCom( uLab) :
    uCom = uLab
    if( noLorentz) : 
        return uCom 
    tcL, tcR = transCoeff( uLab[IMOM]/uLab[IRHO])
    uCom[IRADL] = uLab[IRADL]*tcL
    uCom[IRADR] = uLab[IRADR]*tcR

    return uCom
    
@njit
def comToLab( uCom) :
    uLab = uCom
    if( noLorentz) : 
        return uLab
    tcL, tcR = transCoeff( uCom[IMOM]/uCom[IRHO])
    uLab[IRADL] = uCom[IRADL]/tcL
    uLab[IRADR] = uCom[IRADR]/tcR

    return uLab
  


def runProblem( sigma=1., quiet=False, showPlots=True) : 
    global writePlots
    writePlots = showPlots
    global kappa, rho0
    kappa = sigma/(rho0*wavelength)

    dr, r, grid = generateGrid(rmin=-xmax, rmax=xmax, quiet=quiet)
    t = 0
    tEnd = tcross
    dtFrame = 0.05*tEnd
    if PROBLEM == THERMAL_EQ :
        dtFrame = L/NGRID/(cs0 + cred)*0.3
        tEnd = 20*dtFrame
    step = 0
    frame = 0
    summaryData = []

    while( t < tEnd) :
        if( t>=dtFrame*frame) : 
            plot(t/tcross,r,grid, frame, quiet=quiet, summaryData=summaryData)
            frame+=1

        rho, v, ie, radL, radR = Con2PrimGrid( r, grid)
        dt = TimeStep( dr, rho, v, ie)
        if( t+dt > frame*dtFrame) : 
            dt = frame*dtFrame - t
            #print('step {0} t = {1:5.2E} dt = {2:5.2E}'.format(step, t/tEnd, dt/tEnd))

        if(not quiet) :
            print('step {0} t = {1:5.2E} dt = {2:5.2E}'.format(step, t/tEnd, dt/tEnd))
        step+=1
        grid = AdvanceTimeStep( t, dt, dr, r, grid)
        bc( r, grid)

        t = t+dt
        
    plot(t/tcross,r,grid, frame, quiet=quiet, summaryData=summaryData)
    summaryData = np.array(summaryData)
    if( PROBLEM == THERMAL_EQ and writePlots) :
        pl.clf()
        pl.scatter( summaryData[:,0], summaryData[:,1])
        pl.scatter( summaryData[:,0], summaryData[:,2])
        pl.xlim(summaryData[0,0], summaryData[-1,0])
        pl.savefig("summary.pdf")
    return summaryData

def main() : 
    if False and PROBLEM == LINEAR_WAVE : 
        data = []



    else : 
        parser = argparse.ArgumentParser(description='Process some integers.')
        parser.add_argument('sigma', metavar='N', type=float,
                    help='an integer for the accumulator')
        parser.add_argument('--quiet', action='store_true', default=False)
        parser.add_argument('--no_plots', action='store_true', default=False)
        parser.add_argument('--grid', metavar='N', type=int, default=256 )

        args = parser.parse_args()
        global NGRID
        NGRID=int(args.grid+2)

        summaryData = runProblem( sigma=args.sigma, quiet=args.quiet, showPlots=not args.no_plots)
        np.savetxt(sys.stdout, summaryData)


def waveFunction( time=0.,isCos=True) : 
    def wavef ( x, amp, vph) :
        kL = 1.
        t = time
        if( isCos) :
            return amp*np.cos(math.pi*(kL*x - 2.*vph*t))#*math.cos(2.*math.pi*t)
        else :
            return amp*np.cos(math.pi*kL*x - phase)#*math.sin(2.*math.pi*t)
    
    return wavef

def errorFunction( x, dRho, t = 0., isCos=True) : 
    omega, drho, dv, dT, dp, dEr, dFr = findEigen(quiet=True)
    gamma = omega.imag
    kL = 1.
    vph = omega.real/(2*pi*kL)
    #print(gamma)
    analytic = np.exp(-gamma*t) * np.cos(pi*(kL*x - 2.*vph*t))
    l1Norm = np.abs( dRho/delta - analytic).sum()/x.size
    return l1Norm 

def plot( t, x, grid, i=0, summaryData = None, quiet=False):
    rho, v, ie, radL, radR = Con2PrimGrid( x, grid)

    filename = "movie_frame{0:04d}.png".format(i)
    if( not quiet) : 
        print("Writing file {0}".format(filename))
    ratio = rho/rho0
    M2 = 2.*ratio/((gamma+1.) -(gamma - 1.)*ratio)
    T = (gamma-1.)*ie*mp/kB
    x = x/xmax
    jrad = JRad(radL, radR)
    jrad0 = jrad.mean()
    ie0 = ie.mean()
    if( writePlots and (PROBLEM == SOUND_WAVE or PROBLEM == LINEAR_WAVE)) : 
        pl.clf()
        pl.plot( x, rho/rho0-1, lw=2,label=r"$\rho$")
        pl.plot( x, v/cs0, ls="-.", lw=2, label=r"$v/c_{s,0}$")
        pl.plot( x, jrad/jrad0 - 1., ls="-", lw=2, label=r"$J_r$")
        pl.plot( x, ie/ie0 - 1., ls=":", lw=2, label=r"$\epsilon$")
        pl.text(0.25, -0.75*delta, r"$t = {0:6.4f}$".format(t),fontsize=18)
        pl.ylim(-1.1*delta,delta*1.1)
    elif( writePlots and PROBLEM == RADIATION_PROP) :
        pl.clf()
        pl.plot( x, radL, lw=2, label="$I_1$")
        pl.plot( x, radR, lw=2, label="$I_2$")
        pl.ylabel(r"$I$",fontsize=20)
    elif( PROBLEM == THERMAL_EQ) :
        Trad = (fourPi*JRad(radL, radR)/aRad)**0.25
        if( not summaryData is None) :
            summaryData.append( [t, T.mean(), Trad.mean()])
        print("T={0:5.3e} {1:5.3e}".format(T.mean(), Trad.mean()))
    if( PROBLEM == LINEAR_WAVE or PROBLEM == SOUND_WAVE) :
        popt, pcov = opt.curve_fit(waveFunction(t),x,rho/rho0-1,p0=[delta,1.])
        pl.plot( x, waveFunction(t)(x, popt[0], popt[1]), lw=2,label=r"$\rho$")
        l1Norm = errorFunction( x, rho/rho0-1, t = t, isCos=True) 
        if( not quiet) : 
            print( t, popt)
        summaryData.append([t, popt[0], popt[1], l1Norm])

    if(writePlots) : 
        pl.xlabel("$r$ [kpc]",fontsize=20)
        pl.legend(loc=1,fontsize=18)
        pl.savefig( filename)        

if __name__== "__main__":
    main()
