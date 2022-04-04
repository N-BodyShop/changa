import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl
import math
from numba import jit,njit, prange

#Geometry
CARTESIAN = 0
SPHERICAL = 2
CYLINDRICAL = 1
RIEMANN_TEST = 1
geometry = SPHERICAL
#geometry = CARTESIAN
mode=0
#mode=RIEMANN_TEST

#constants
pi = math.pi
Rsun = 7e10
kB = 1.38e-16
mp = 1.6e-24
yrs = 24*3600*365.24
kms = 1e5

NVARS = 3
IRHO = 0
IMOM = 1
IENER = 2

#run parameters
courant = 0.2
gamma = 5./3.
NGRID=1000
mu = 0.588*mp
T_atm = 1e4
Tmax =1e7
rho_atm = 1e-14

writeFrames = True

@njit
def geoFactor(x):
#
# Geometry factors 
#
    if(geometry == CARTESIAN) :
        return 1.
    elif(geometry == SPHERICAL) : 
        return x*x
    elif(geometry == CYLINDRICAL) : 
        return x

def initial_conditions( r) : 
    rho, vr, temp = 0,0,0

    rho1 = 1e1**1.259*r**-0.62
    rho2 = 1e1**30.42*r**-2.89

    rho = 1./np.sqrt(1/rho1**2 + 1/rho2**2)

    vr = np.maximum(-3.63e5 + r*9.588e-8, 0)

    T1 = 1e1**20.35*r**-1.239
    T2 = 1e1**8*r**-0.287
    temp = 1/np.sqrt(1/T1**2 + 1/T2**2)

    rho[r > 1000*Rsun] = rho_atm
    vr[r > 1000*Rsun] =  0
    temp[r > 1000*Rsun] = T_atm

    return rho, vr, temp

def generateGrid( rmin=0.001*Rsun, rmax=1000*Rsun) :
# Generates the grid and sets initial conditions on primitives
    dx = np.ones(NGRID,dtype=np.double)*rmax/NGRID
    if( geometry == CARTESIAN) :
        rmin = 0.
    xCoord = dx*(np.arange(NGRID)+0.5)+rmin

    rho, v, T = initial_conditions( xCoord)
    ie  = kB*T/mu/(gamma-1.)

    return dx, xCoord, Prim2ConGrid( xCoord, rho, v, ie)


def Prim2ConGrid(xcoord, rho, v, ie) :
#
# Primitive to conservative solver (for the grid)
#
    grid = np.zeros([NGRID,NVARS],dtype=np.double)
    for i in range(xcoord.size) :
        grid[i,:] =  Prim2Con(xcoord[i],rho[i],v[i],ie[i])
    return grid

@njit
def Prim2Con( r, rho, v, ie) : 
#
# Primitive to conservative solver
#
    u = np.zeros(NVARS,dtype=np.double)
    gFactor = geoFactor(r)
    u[IRHO] = rho*gFactor
    u[IMOM] = rho*v*gFactor
    u[IENER] = rho*(ie+0.5*v*v)*gFactor
    return u

@njit
def Con2PrimGrid(xcoord, grid) :
#
# Conservative to Primitive Solver on the grid
#
    rho = np.zeros(grid.shape[0],dtype=np.double)
    v = np.zeros(grid.shape[0],dtype=np.double)
    ie = np.zeros(grid.shape[0],dtype=np.double)
    for i in range(grid.shape[0]) : 
        dens, vel, iener = Con2Prim(xcoord[i], grid[i,:])
        rho[i] = dens
        v[i] = vel
        ie[i] = iener
    return rho, v, ie

@njit
def Con2Prim(r, u) :
#
# Conservative to primite solver
#
    gFactor = geoFactor(r)
    rho = u[IRHO]/gFactor
    v = u[IMOM]/(rho*gFactor)
    ie = (u[IENER]/gFactor/rho - 0.5*v*v)
    return rho, v, ie

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
    	FHLL = FL
    elif ( sR <= 0.):
        FHLL = FR
    elif ((sL < 0.) and (sR >= 0.)) :
    	FHLL = (FL * sR - FR * sL + (UR - UL) * sL * sR) / (sR - sL)
    return FHLL

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
    return courant*np.min( dr/(np.abs(v) +cs))

@njit
def calFlux( r, uL, uR) : 
#
# Calculate the flux in based in reconstructed left and right states
#
    rhoL, vL, iEL = Con2Prim(r, uL)
    rhoR, vR, iER = Con2Prim(r, uR)
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
def reconstruct( dx, xcoord, grid) :
#
# reconstruction -- piecewise constant for now
#
    rho, v, ie = Con2PrimGrid(xcoord, grid)
    uL = grid*0.
    uR = grid*0.
    for i in range(xcoord.size) :
        uL[i,:] = Prim2Con( xcoord[i] - dx[i]*0.5, rho[i], v[i], ie[i])
        uR[i,:] = Prim2Con( xcoord[i] + dx[i]*0.5, rho[i], v[i], ie[i])
    return uL, uR

@njit
def bc(xcoord,grid) :
#
# impose boundary conditions
#
    rho0, v0, ie0 = Con2Prim(xcoord[1], grid[1,:])
    v0 = 0.  

    grid[0,:] = Prim2Con(xcoord[0],rho0,v0,ie0)
    #rho, v, ie = Con2Prim(xcoord[-2], grid[-2,:])
    ie  = kB*T_atm/mu/(gamma-1.)
    grid[-1,:] = Prim2Con(xcoord[-1],rho_atm,0,ie)

    return grid

@njit
def AdvanceTimeStep( t, dt, dr, r, grid) :
#
# make a timestep -- first order for now
#
    uL, uR = reconstruct( dr, r, grid) 

    gridNew = grid.copy()
    for i in range(1,grid.shape[0]-1) :
        x = r[i]
        dx = dr[i]
        ul = uR[i]
        ur = uL[i+1]
        fl, csl, fr, csr = calFlux(x+dx, ul, ur)
        fluxR = HLLE( x, csl, csr, fl, fr, ul, ur) 
        ul = uR[i-1]
        ur = uL[i]
        fl, csl, fr, csr = calFlux(x-dx, ul, ur)
        fluxL = HLLE( x, csl, csr, fl, fr, ul, ur)
        gFactor = geoFactor(x)
        rho, v, ie = Con2Prim(x, grid[i,:])
        P = rho*ie*(gamma-1.)
        if( fluxR.shape[0] < 3) :
            print("Here", i, dt, fluxR, fluxL)
        gridNew[i] -= (fluxR-fluxL)/dr[i]*dt 
        gridNew[i,IMOM] += geometry*gFactor/x*P*dt
        gridNew[i] += sources(t, x, grid[i,:])*dt
    
    return gridNew

@njit
def sources( t, r, u) :
#
# Impose sources
#
    source = np.zeros(u.shape)
    #if(mode != RIEMANN_TEST and geometry==CYLINDRICAL and r < 10*rjet and t < tjet) :
    #    gFactor = geoFactor(r)
    #    source[IENER] = Ljet/Dpp/(pi*rjet*rjet)*(2.*pi*gFactor)*math.exp(-r*r/(rjet*rjet))
    return source

def main() : 
    dr, r, grid = generateGrid(rmin = 10*Rsun, rmax = 1e4*Rsun)
    t = 0
    tEnd = 10*yrs
    frame = 0
    dtFrame = 1e-1*yrs
    step = 0
    while( t < tEnd) :
        if( t>=dtFrame*frame) : 
            plot(t,r,grid, frame)
            frame+=1

        rho, v, ie = Con2PrimGrid( r, grid)
        dt = TimeStep( dr, rho, v, ie)
        if( t+dt > frame*dtFrame) : 
            dt = frame*dtFrame - t
        if( step % 100 == 0) : 
            print('step {0} t = {1:5.2E} dt = {2:5.2E} '.format(step, t/yrs, dt/yrs))
        step+=1
        grid = AdvanceTimeStep( t, dt, dr, r, grid)
        t = t+dt
        grid = bc( r, grid)
        
        # clean data
        #rho, v, ie = Con2PrimGrid( r, grid)
        #ie_floor = kB*T_atm/mu/(gamma-1.)
        #ie = np.maximum(ie, ie_floor)
        #grid1 = Prim2ConGrid(r, rho, v, ie) 
        #print((grid[:,2]-grid1[:,2])/grid[:,2])
    plot(t,r,grid, frame)

def plot( t, x, grid, i=0):
    if not writeFrames:
        return
    rho, v, ie = Con2PrimGrid( x, grid)

    filename = "movie_frame{0:04d}.png".format(i)
    print( "Writing file {0}".format(filename))
    pl.clf()
    T = (gamma-1.)*ie*mu/kB
    x = x/Rsun
    pl.loglog( x/Rsun, rho/rho_atm, lw=2,label=r"$\rho$")
    pl.loglog( x/Rsun, T/1e3, ls="dotted", lw=2,label=r"$\log_{10}(T_3)$")
    pl.semilogx( x/Rsun, v/kms, ls="-.", lw=2, label="$v$")
    pl.ylabel(r"$\rho$, $M$, $\log_{10}(T_3)$, $v/c_{s,0}$",fontsize=20)
    pl.xlabel("$r$ [kpc]",fontsize=20)
    #pl.text(2500, 0.25, r"$r_{\rm jet}"+r" = {0:4.2f}$ Mpc".format(rjet/Mpc),fontsize=18)
    #pl.text(2500, 0.5, r"$t={0:3.2f}$ Gyrs".format(t/Gyrs),fontsize=18)
    pl.legend(loc=1,fontsize=18)
    #pl.ylim(0,3)
    #pl.xlim(0,1e4)
    pl.savefig( filename)        

if __name__== "__main__":
    main()
