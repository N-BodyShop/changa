import sys
sys.path.append("..")
import starmaker as sm
import scipy.integrate as si
import math
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as pl

import numpy as np
try :
	import nuc_eos
except : 
	print( "cannot import nuc_eos")
import argparse



kB = 1.38e-16
mp = 1.6e-24
G = 6.67384e-8

def NSderivs( r, y, temp, ye) :
    P = y[0]
    m = y[1]
    #print(temp/1e8)
    rho = nuc_eos.rho_find(P, temp, ye)

    dPdr = -rho*G*m/(r*r)
    dMdr = 4.*math.pi*rho*r*r

    return [dPdr, dMdr]


def makeSN() : 
    sm.rscaling = 2
    sm.maxRatio = 1000
    sm.N = 20
    rModel, rhoModel, eModel, yeModel, vModel = np.loadtxt( "s20WH07_NW_15.19ms_LS220.flash", skiprows=78, usecols=[0,1,2,3,4], unpack=True)
    vModel = vModel*G**-0.5
    mModel = np.zeros(rModel.size)
    mModel[1:] = 4.*3.1415/3.*(rModel[1:]**3 - rModel[:-1]**3) * np.sqrt(rhoModel[1:]*rhoModel[:-1])
    mModel[0] = 4.*3.1415/3.*rModel[0]**3*rhoModel[0]
    mModel = np.cumsum(mModel)
    # convert to code units
    massModel = mModel.max()

    pModel = []
    for rho, T, ye in zip( rhoModel, eModel, yeModel) :
        pModel.append( nuc_eos.press( rho, T, ye, p0=None)[0])
    pModel = np.array(pModel)
    dPdr = -(pModel[1:] - pModel[:-1])/(rModel[1:] - rModel[:-1])
    gModel = G*mModel/rModel**2
    rhog = (rhoModel*gModel)[1:]
    pl.loglog( rModel[1:], dPdr)
    pl.loglog( rModel[1:], rhog)
    pl.loglog( rModel, pModel*1e-5)
    #for r, x, y, rho, m  in zip( rModel[1::10], dPdr[::10], rhog[::10], rhoModel[::10], mModel[::10]) :
    #    print(r,x, y, rho, m)
    #print(dPdr, rhog)
    #pl.loglog( rModel, rhoModel/rhoModel[0],label=r"$\rho$")
    #pl.loglog( rModel, rhoModel,label=r"$P$")
    #pl.loglog( rModel, mModel/massModel,label=r"$M$")
    #pl.loglog( rModel, -vModel*G**0.5/3e10,label=r"$-v$")
    pl.ylim(1e23, 1e29)
    pl.xlim(1e5,1e7)
    pl.legend(loc="best")
    pl.savefig("sn.pdf")
    np.savetxt("sn.out", np.array(list(zip(rModel, pModel, rhoModel*gModel))))
    print( massModel)
    #return mModel, rModel, eModel, vModel
    #mModel, rModel, eModel = makeNS( massNS, temp=tempNS)

    return 0, 0, mModel, rModel, rhoModel, eModel, yeModel, vModel

def makeNS( mass, temp=1e7, ye=0.1, tol=1e-3) :
    # guess a density
    rhoc = 1e14
    #compute the central pressure

    m = 0.
    r = 0.
    Pc,epsc,entc,cs2c = nuc_eos.press( rhoc, temp, ye)
    print( Pc)
    #return
    mArray = None
    rArray = None
    rhocPrev = -1
    massPrev = -1

    while abs((m-mass)/mass) > tol :
        integral = si.ode( NSderivs).set_integrator("dop853")

        Pc,epsc,entc,cs2c = nuc_eos.press( rhoc, temp, ye)
        #print("central pressure", Pc, rhoc/1e15)
        r = 1e1
        m = 4.*math.pi/3. * rhoc*r**3
        P = Pc - 2.*math.pi/3.*G*r**2*rhoc**2
        y0 = np.array([P, m])
        dr = r*0.1
        integral.set_initial_value(y0, r).set_f_params( temp, ye)
        mArray = [0.]
        rArray = [0.]
        rhoArray = [rhoc]
        rho = rhoc
        while integral.successful() and P > 1e-5*Pc and rho > 1e9:
            #print P, m
            P, m = integral.integrate( integral.t+dr)
            
            #y = si.odeint( derivs, y, [r, r+dr], args=(temp, abar, zbar))
            r = integral.t
            rArray.append(r)
            mArray.append(m)

            #compute new dr
            dPdr, dMdr = NSderivs( r, [P,m], temp, ye)
            dr = min(0.25*min(abs(P/dPdr), abs(m/dMdr)),0.1*r)
            rho = nuc_eos.rho_find(P, temp, ye)
            #print P, abs(P/dPdr)/r, dr/r, rho, m
            rhoArray.append(rho)
            #print y, r
            
        if( massPrev < 0) :
            rhocPrev = rhoc
            massPrev = m
            rhoc = rhoc * 1.1 #(1.0 + 0.1*(mass-m)/mass)
        else :
            rhocTemp = rhoc
            DeltaRho = -(rhoc - rhocPrev)/(m - massPrev) * (m-mass)
            if( abs(DeltaRho)/rhoc < 0.1) :  
                rhoc = rhoc + DeltaRho
            else :
                if( DeltaRho > 0 ) :
                    rhoc = 1.1 * rhoc
                else :
                    rhoc = 0.9 * rhoc
            
            rhocPrev = rhocTemp
            massPrev = m
            
        print( "new rhoc = {0:5.3e} {1:5.3e} {2:5.3e}".format(rhoc, m, np.array(rArray).max()))
        
    return 0, 0, np.array(mArray), np.array(rArray), np.array(rhoArray), temp*np.ones(len(mArray)), ye*np.ones(len(mArray))


mSun = 1.99e33

nuc_eos.initialize()
xmax = 25e6
tempNS = 1e9
massNS = 1.4*mSun
massParticleInGram = mSun/2e4

massModel = massNS
rhoExt = 1e9
tempExt = 2e9
rho_max = 1e14
ymax = xmax
zmax = xmax

parser = argparse.ArgumentParser(prog='PROG')
parser.add_argument('type', default="ns", help="ic to make (sedov, evrard, star, merger)")

args = parser.parse_args()

mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel = None,None,None,None,None,None,None, None
rMax = None

if( args.type == "sn") :
 
    xmax = 1e9
    ymax = xmax
    zmax = xmax

    rhoExt = 150
    tempExt = 1e9
    rho_max = 1e14

    rmax = 1.05*math.sqrt(xmax**2 + ymax**2 + zmax**2)
    ExtMass = 4.*3.1415/3.*rhoExt*rmax**3
    mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel = makeSN()
    rMax=1e7

else :
    mcore, rcore, mModel, rModel, rhoModel, eModel, XModel = makeNS(massNS, temp=tempNS, ye=0.1)

sm.maxRatio = 2000
sm.makeStar(mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, vModel, xmax, ymax, zmax, massParticleInGram, rhoExt, tempExt, rMax=rMax,outfile=args.type) 
