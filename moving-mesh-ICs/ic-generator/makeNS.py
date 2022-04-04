import starmaker as sm
import scipy.integrate as si
import math

import numpy as np
try :
	import nuc_eos
except : 
	print( "cannot import nuc_eos")

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
massParticleInGram = mSun/2e5

massModel = massNS
rhoExt = 1e9
tempExt = 1e9
rho_max = 1e14

ymax = xmax
zmax = xmax

mcore, rcore, mModel, rModel, rhoModel, eModel, XModel = makeNS(massNS, temp=tempNS, ye=0.1)
sm.makeStar(mcore, rcore, mModel, rModel, rhoModel, eModel, XModel, xmax, ymax, zmax, massParticleInGram, rhoExt, tempExt, outfile="ns") 
