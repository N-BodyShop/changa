from ctypes import *
import scipy.optimize as so
libnuc = None

def initialize(file="LS220_234r_136t_50y_analmu_20091212_SVNr26.h5") : 
    global libnuc
    cdll.LoadLibrary("libnuc_eos.so")
    libnuc = CDLL("libnuc_eos.so")
    libnuc.nuc_eos_C_ReadTable_for_c(file.encode('ascii'))

def rho_find( pres, temp, ye) :
    rho0 = 1e14
    PRESSGF=1.80123683248503e-39
    EPSGF=1.11265005605362e-21

    return so.newton( press, rho0, tol=1e-4, args=(temp,ye,pres*PRESSGF))

def press( rho, temp, ye, p0=None) : 
    n=c_int(1)
    RHOGF=1.61887093132742e-18
    xrho=c_double(rho*RHOGF)
    xtemp=c_double(temp)
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
    global libnuc
    libnuc.nuc_eos_m_kt1_short_for_c(byref(n),byref(xrho),byref(xtemp),byref(xye), 
	byref(xeps),byref(xprs),byref(xent),byref(xcs2),byref(xdedt),
	byref(xdpderho),byref(xdpdrhoe),byref(xmunu),
		      byref(keyerr),byref(anyerr))
    if( p0) : 
        return xprs.value - p0
    energy_shift = c_double.in_dll(libnuc, "nuc_eos_energy_shift_for_c").value * 0.
    print(energy_shift)
    INVRHOGF = 6.17714470405638e17
    INVEPSGF = 8.98755178736818e20
    INVPRESSGF = 5.55174079257738e38
 
    return xprs.value*INVPRESSGF,(xeps.value+energy_shift)*INVPRESSGF,xent.value,xcs2.value

if __name__ == "__main__":
    initialize() 
    p,eps, ent,cs2 = press(2.00003e10,0.1,0.1)
    print(p, eps)
    print(rho_find(p, 0.1, 0.1)/1e10)
