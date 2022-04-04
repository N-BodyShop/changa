from ctypes import *
import scipy.optimize as so
libnuc = None

LENGTHGF=6.77269222552442e-06
TIMEGF=2.03040204956746e05
RHOGF=1.61887093132742e-18
PRESSGF=1.80123683248503e-39
EPSGF=1.11265005605362e-21
INVRHOGF=6.17714470405638e17
INVEPSGF=8.98755178736818e20
INVPRESSGF=5.55174079257738e38
K_TO_MEV=8.625e-11
import math 


def initialize(file="LS220_234r_136t_50y_analmu_20091212_SVNr26.h5") : 
    global libnuc
    cdll.LoadLibrary("./libnuc_eos.so")
    libnuc = CDLL("libnuc_eos.so")
    libnuc.nuc_eos_energy_shift_for_c.restype = c_double
    libnuc.nuc_eos_C_ReadTable_for_c( file.encode('ascii'))

def rho_find( pres, temp, ye, rho0 = None) :
    rho = 0.
    lgP = math.log10(pres)
    if( rho0 == None) : 
        if( lgP < 32.5) :
            lgRho = (lgP - 24.5)/1.3 + 8
        else :
            lgRho = (lgP - 32.5)/3 + 14
        rho0 = 1e1**lgRho
    try : 
        #print("Here: ", pres, rho0/1e15, temp*K_TO_MEV)
        rho = so.newton( press, rho0, rtol=1e-4, args=(temp,ye,pres),maxiter=50)
    except( RuntimeError ):
        print( pres, temp*K_TO_MEV, ye)
        #print( e.message)
        raise
    return rho

def pressRhoE( rho, E, ye) : 
    global libnuc
    n=c_int(1)
    energy_shift = libnuc.nuc_eos_energy_shift_for_c()

    xrho=c_double(rho*RHOGF)
    xtemp=c_double(1e8*K_TO_MEV)
    xye=c_double(ye)
    xeps=c_double(E*EPSGF-energy_shift)
    xprs=c_double()
    xent=c_double()
    xcs2=c_double()
    xdedt=c_double()
    xdpderho=c_double()
    xdpdrhoe=c_double()
    xmunu=c_double()
    xprec=c_double(1e-5)

    keyerr=c_int()
    anyerr=c_int()

    libnuc.nuc_eos_m_kt0_short(byref(n),byref(xrho),byref(xtemp),byref(xye), 
	byref(xeps),byref(xprs),byref(xent),byref(xcs2),byref(xdedt),
	byref(xdpderho),byref(xdpdrhoe),byref(xmunu),byref(xprec),
		      byref(keyerr),byref(anyerr))

    if(keyerr.value != 0) :
        print( keyerr.value, rho, E, ye)
    #energy_shift = c_double.in_dll(libnuc, "nuc_eos_energy_shift_for_c").value
    return xprs.value*INVPRESSGF,(xeps.value+energy_shift)*INVEPSGF,xent.value,xcs2.value,xtemp.value/K_TO_MEV

def press( rho, temp, ye, p0=None) : 
    global libnuc
    n=c_int(1)

    xrho=c_double(rho*RHOGF)
    xtemp=c_double(temp*K_TO_MEV)
    xye=c_double(ye)
    xeps=c_double()
    xprs=c_double()
    xent=c_double()
    xcs2=c_double()
    xdedt=c_double()
    xdpderho=c_double()
    xdpdrhoe=c_double()
    xmunu=c_double()

    keyerr=c_int()
    anyerr=c_int()
    libnuc.nuc_eos_m_kt1_short(byref(n),byref(xrho),byref(xtemp),byref(xye), 
	byref(xeps),byref(xprs),byref(xent),byref(xcs2),byref(xdedt),
	byref(xdpderho),byref(xdpdrhoe),byref(xmunu),
		      byref(keyerr),byref(anyerr))
    if(keyerr.value != 0) :
        print( keyerr.value, rho, temp, ye)
    if( p0) : 
        return xprs.value*INVPRESSGF/p0 - 1.0
    energy_shift = libnuc.nuc_eos_energy_shift_for_c()
    return xprs.value*INVPRESSGF,(xeps.value + energy_shift)*INVEPSGF,xent.value,xcs2.value

if __name__ == "__main__":
    initialize()
    import numpy as np
    import math
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as pl
    pArray = []
    c2 = 9e20
    lgRhoArray = np.arange( 10.0,15.0,0.5) 
    #print( press( 1.154e+13, 1e8, 0.1))
    #print( pressRhoE( 5.04e14, 1.949e21, 0.1))
    for lgrho in lgRhoArray :
        p,eps, ent,cs2 = press(1e1**lgrho,1e9, 0.1)
        p2, eps2, ent2, _, T2 = pressRhoE( 1e1**lgrho, eps, 0.1)
        gamma = cs2*c2/(p/1e1**lgrho)
        print( "{0:5.3e} {1:5.3e} {2:5.3e} {3:5.3e} {4:5.3e} {5:5.3e}".format( lgrho, p, p2, eps2/eps,  1e1**lgrho*eps*(gamma-1), gamma))
        pArray.append(math.log10(p))
    pl.plot( lgRhoArray, pArray)
    pAnalytic = 10.**(24.5 + 1.3*(lgRhoArray - 8)) + 10**(32.5 + 3*(lgRhoArray - 14))
    pl.plot( lgRhoArray, np.log10(pAnalytic))
    pl.savefig( "test.pdf")
