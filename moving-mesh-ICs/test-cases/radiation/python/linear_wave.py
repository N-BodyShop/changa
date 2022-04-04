
#import athena_read
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import scipy.optimize as opt
import math

class container:
    """Simple class to store input"""
    pass


def plot_conv(infile='error_py.out'):
    """
    Plots convergence from error file generated with get_error_tab()
    """
    data = np.loadtxt(infile)
    print( data.shape)
    #nx=[32,64,128,256,512]
    #err=[1.18e-8,3.91e-9,1.57e-9,7.15e-10,3.44e-10]
    plt.plot(data[:,0],data[:,1],'s')

    plt.yscale('log')
    plt.xscale('log')
    plt.savefig("convergence.pdf",bbox_inches='tight')
    plt.close()

def read_error_file(fname):

    infile = open(fname,"r")
    data = container()
    line = infile.readline().split()
    data.prat  = float(line[0])
    data.crat  = float(line[1])
    data.sigma = float(line[2])
    line = infile.readline().split()
    data.omegar  = float(line[0])
    data.omegai  = float(line[1])
    line = infile.readline().split()
    data.delvr  = float(line[0])
    data.delvi  = float(line[1])
    line = infile.readline().split()
    data.delpr  = float(line[0])
    data.delpi  = float(line[1])
    line = infile.readline().split()
    data.delerr  = float(line[0])
    data.deleri  = float(line[1])
    line = infile.readline().split()
    data.delfr  = float(line[0])
    data.delfi  = float(line[1])
    line = infile.readline().split()
    data.time  = float(line[0])
    data.derr  = float(line[1])
    data.verr = float(line[2])
    data.perr = float(line[3])

    return data

def get_error_tab(infile='errorfile',athfile='radwave.out5.00010.tab',outfile='error.pdf',
                  d0=1,p0=1,amp=1.0e-6,dimensions=2):
    """
    Reads in error output and tabulated ouptput of Athena.  Error contains eigenvalues
    and eigenvectors for radiative linear wave which are used to construct solution at
    end of run. L1 error is computed using this solution and the athena formatted data
    output.
    """

    # read in Athena error output
    er = read_error_file(infile)

    # read in Athena data
    data = None #athena_read.tab(athfile,dimensions=dimensions)

    nx3=data.shape[0]
    nx2=data.shape[1]
    nx1=data.shape[2]

    rho = data[:,:,:,2]
    x = data[:,:,:,0]

    time = er.time
    kw = 2.*np.pi
    vp = er.omegar/kw
    ampt = amp/np.exp(er.omegai*time)
    count = nx1 * nx2 *nx3
    rhot = np.empty([nx3,nx2,nx1])
    x1v = np.empty(nx1)
    error = 0.0
    for k in xrange(nx3):
        for j in xrange(nx2):
            for i in xrange(nx1):
                x0=x[k,j,i]
                x1v[i] = x0
                cosx = np.cos(kw*(vp*er.time-x0))
                sinx = np.sin(kw*(vp*er.time-x0))
                rhot[k,j,i] = d0+ampt*cosx
                #if (j == 0):
                    #print x0,rho[k,j,i]-d0,rhot[k,j,i]-d0
                error = error + abs(rhot[k,j,i]-rho[k,j,i])

    # compute error and write to file
    error = error/(nx1*nx2*nx3)
    f = open('error_py.out', 'a')
    f.write("{:d} {:e}\n".format(nx1,error))
    f.close()
    print( error,er.derr)

    # generate plots
    plt.plot(x1v,rho[0,0,:]-d0)
    plt.plot(x1v,rhot[0,0,:]-d0,'--')

    plt.savefig(outfile,bbox_inches='tight')
    plt.close()


def sinusoid(x, wrt, wit):
    """
    Returns a sinusoidal function for fitting routine
    """
    amp = 1.0e-6
    time = 1.0
    ai = 0.
    return amp*np.exp(-wit)*(np.cos(wrt*time-2.*np.pi*x)+
                                 ai*np.sin(wrt*time-2.*np.pi*x))

def fit_wave(infile='errorfile',athfile='radwave.block0.out5.00010.tab',
                  d0=1,p0=1,amp=1.0e-6,dimensions=2,plot=False):
    """
    Function to read in Athena tabulated output and fit it to a sinusoidal
    function.  Uses scipy.optimize
    """

    # read in error output from Athena
    er = read_error_file(infile)

    # read in data from athena
    data = None #athena_read.tab(athfile,dimensions=dimensions)

    nx3=data.shape[0]
    nx2=data.shape[1]
    nx1=data.shape[2]

    # set rho and x from athena output
    drho = data[:,:,:,2].flatten()-d0
    x = data[:,:,:,0].flatten()

    #insert fitting stuff here


    popt, pcov = opt.curve_fit(sinusoid,x,drho,p0=[2.*np.pi,0.0])
    print( popt)
    #print pcov

    if (plot):
        x1v = np.empty(nx1)
        rhof = np.empty(nx1)
        for i in xrange(nx1):
            x1v[i] = data[0,0,i,0]
            rhof[i] = sinusoid(x1v[i],popt[0],popt[1])

        plt.rc('text',usetex=True)
        plt.plot(x1v,data[0,0,:,2]-d0)
        plt.plot(x1v,rhof,'--')
        plt.xlabel(r"$x$")
        plt.ylabel(r"$\delta \rho$")
        plt.savefig('fit.pdf',bbox_inches='tight')
        plt.clf()

    return popt

def dispersion_solution(ka=2.*np.pi,sigma=1.0,crat=10.0,prat=1.0,gamma=1.666666667,r=1.,all=False):
    """
    Solve the dispersion relation for a radiative linear wave
    """

    gam1 = gamma-1.
    #ka = 2.*np.pi
    k2 = ka*ka
    k4 = k2*k2
    c2 = crat*crat
    c3 = c2*crat
    sig2 = sigma*sigma
    r2 = r * r

    coeff = np.empty(6,dtype=complex)
    # Define coefficients in terms of physical variables
    coeff[5] = 4*c3*k4/3.*r2
    coeff[4] = 1j*(c2*k4*gamma/(3*prat*gam1*sigma)*r2+4.*(3.+2.*r)*r/3.*c2*k2*sigma
                   +16./9.*c2*k2*r2*prat*sigma+c2*k2*r2*gamma*sigma/(prat*gam1))
    coeff[3] = -(4*crat*k2+4.*c3*k2/3.*r2+4./9.*crat*k2*r2/gam1+2*crat*r*k2*gamma/(prat*gam1))
    coeff[2] = -1j*(c2*k2*r2+3.*k2*gamma+3.*c2*r2*sig2+4.*prat*r*sig2+12.*c2*r*prat*sig2*gam1
                    +16.*prat*prat*sig2*gam1)/(3.*prat*gam1*sigma)
    coeff[1] = (4.*crat+4./(3.*crat*gam1)+2.*crat*r/(prat*gam1))
    coeff[0] = 1j/(prat*gam1*sigma)

    roots = np.roots(coeff)
    #print roots

    if (all):
        return roots

    first = True
    for root in roots:
        if ((root.real > 0.25*ka) and (root.imag > 0.0)):
            if (first):
                rootc = root
                first = False
            else:
                if (root.real < rootc.real):
                    rootc = root

    if (first):
        print( roots)
    return rootc

   # for root in roots:
   #     if ((root.real > 0.25*ka) and (root.real < 2.2*ka) and (root.imag > 0.0)):
   #         return root

def eigenvec(ka=2.*np.pi,sigma=1.0,crat=10.0,prat=1.0,gamma=1.666666667,r=1., quiet=False):

   omega = dispersion_solution(ka=ka,sigma=sigma,crat=crat,prat=prat,gamma=gamma,r=r)
   gam1 = gamma-1.
   #ka = 2.*np.pi
   k2 = ka*ka
   k4 = k2*k2
   c2 = crat*crat
   c3 = c2*crat
   o2 = omega*omega
   o3 = omega*o2
   sig2 = sigma*sigma
   r2 = r * r
   drho = 1.
   dv = omega/ka*drho
   dp = (12*1j*prat*sigma*o3-2*1j*c2*r*sigma*omega*(2*k2*prat*r-9*o2)+3*crat*o2*(4.*prat*r*sig2
         -3*o2)+3*c3*r2*(3*sig2*o2+k2*(4*prat*sig2+o2)))/(3*crat*k2*(c2*r2*(k2+(3+4*prat)*sig2)
         +6*1j*crat*r*sigma*omega-3*o2))
   dT = dp - drho
#   dT = 1j*gam1*omega*(-3.*omega*omega+6.*1j*crat*omega*r*sigma+c2*r2*
#         (k2+(3.+4.*prat)*sig2))/(-3.*1j*omega*omega*omega+4.*c3*
#         gam1*k2*prat*r2*sigma-6.*crat*omega*omega*(r+2.*gam1*prat)*sigma
#         +1j*c2*omega*r*(k2*r+3.*(r+4.*gam1*prat)*sig2))
#   dp = dT + drho
   dEr = 4*dT+(1j*omega*dp-1j*gamma*ka*dv)/(gam1*sigma*crat*prat)
   dFr = ((1j*omega/r+crat*sigma)/crat*dEr-4*sigma*dT)/(1j*ka)
   if not quiet : 
      print( "omega: ", omega.real,omega.imag, omega.imag/omega.real)
      print( "dv:  ", dv.real,dv.imag, dv.imag/dv.real)
      print( "dT:  ", dT.real,dT.imag, dT.imag/dT.real)
      print( "dp:  ",dp.real,dp.imag, dp.imag/dp.real)
      print( "dEr: ", dEr.real,dEr.imag, dEr.imag/dEr.real)
      print( "dFr: ",dFr.real,dFr.imag, dFr.imag/dFr.real)
   return omega, drho, dv, dT, dp, dEr, dFr

def read_dispersion(fname='disp.tab'):

    data = np.loadtxt(fname)
    print( data)


def plot_dispersion(crat=[1.e3],prat=1.,r=[1.],trange=[1.e-3,1.e3],ka=2.*np.pi,gamma=1.666666667,
                    outfile='disp.pdf'):
    """
    Plot dispersion as function of tau
    """
    nx=100
    tmax = np.log10(trange[1]/trange[0])
    sigma = np.empty(nx)
    sigma = 10**(np.arange(0,nx)/(nx-1.)*tmax)*trange[0]
    ncr = len(crat)
    nr = len(r)
    wr = np.empty((nx,ncr,nr))
    wi = np.empty((nx,ncr,nr))


    for k in xrange(ncr):
        for j in xrange(nr):
            for i in xrange(nx):
                root = dispersion_solution(ka=ka,crat=crat[k],prat=prat,r=r[j],gamma=gamma,sigma=sigma[i])
                wr[i,k,j] = root.real
                wi[i,k,j] = root.imag

    sb1 = plt.subplot(211)
    for k in xrange(ncr):
        for j in xrange(nr):
            plt.plot(sigma,wr[:,k,j]/ka)

    plt.xscale('log')
    sb2 = plt.subplot(212)
    for k in xrange(ncr):
        for j in xrange(nr):
            plt.plot(sigma,wi[:,k,j])
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(outfile)
    plt.close()


def plot_dispersion_reduce(crat=1.e3,prat=1.,r=[1.],trange=[1.e-3,1.e3],ka=2.*np.pi,gamma=1.666666667,
                           outfile='disp.pdf',fname='disp.tab',orlim=None,oilim=None):
    """
    Plot dispersion as function of tau
    """
    nx=100
    tmax = np.log10(trange[1]/trange[0])
    sigma = np.empty(nx)
    sigma = 10**(np.arange(0,nx)/(nx-1.)*tmax)*trange[0]
    nr = len(r)
    wr = np.empty((nx,nr))
    wi = np.empty((nx,nr))

    ath1 = np.loadtxt('disp_1.tab')
    ath2 = np.loadtxt('disp_0.1.tab')
    ath3 = np.loadtxt('disp_0.01.tab')

    for j in xrange(nr):
        for i in xrange(nx):
            root = dispersion_solution(ka=ka,crat=crat,prat=prat,r=r[j],gamma=gamma,sigma=sigma[i])
            wr[i,j] = root.real
            wi[i,j] = root.imag

    plt.rc('text',usetex=True)
    plt.rc('font',family='serif')
    pcl = ['k','r','b','g']
    label = ['1','0.1','0.01']
    sb1 = plt.subplot(211)
    for j in xrange(nr):
        plt.plot(sigma,wr[:,j]/ka,pcl[j],linewidth=2,label=label[j])
    plt.xscale('log')
    plt.xlim(trange)
    plt.text(100.,1.4,'C=1000, P=1',fontsize=16)
    if (orlim is not None):
        plt.ylim(orlim)
    plt.plot(ath1[:,2],ath1[:,4]/ka,'ok')
    plt.plot(ath2[:,2],ath2[:,4]/ka,'or')
    plt.plot(ath3[:,2],ath3[:,4]/ka,'ob')
    plt.ylabel(r"$\omega_{\rm R}/k$",fontsize=16)
    sb1.axes.get_xaxis().set_ticklabels({})
    plt.legend(loc=0)
    sb2 = plt.subplot(212)
    for j in xrange(nr):
        plt.plot(sigma,wi[:,j],pcl[j],linewidth=2)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(trange)
    plt.ylabel(r"$\omega_{\rm I}$",fontsize=16)
    plt.xlabel(r"$\tau$",fontsize=16)
    if (oilim is not None):
        plt.ylim(oilim)
    plt.plot(ath1[:,2],ath1[:,5],'ok')
    plt.plot(ath2[:,2],ath2[:,5],'or')
    plt.plot(ath3[:,2],ath3[:,5],'ob')
    #plt.plot(ath[:,2],ath[:,5],'o')
    plt.tight_layout(h_pad=0.1,w_pad=0.5)
    plt.savefig(outfile)
    plt.close()

def plot_dispersion_all(trange=[1.e-3,1.e3],ka=2.*np.pi,crat=10.0,prat=1.0,gamma=1.666666667,
                    outfile='disp.pdf'):
    """
    Plot dispersion as function of tau
    """
    nx=100
    tmax = np.log10(trange[1]/trange[0])
    sigma = np.empty(nx)
    sigma = 10**(np.arange(0,nx)/(nx-1.)*tmax)*trange[0]
    wr = np.empty((nx,5))
    wi = np.empty((nx,5))
    realr = np.empty(5)
    imagr = np.empty(5)
    for i in xrange(nx):
        roots = dispersion_solution(ka=ka,crat=crat,prat=prat,gamma=gamma,sigma=sigma[i],all=True)
        rootss = sorted(roots, key=lambda roots: roots.real)

        #for j in xrange(5):
        #    realr[j] = roots[j].real
        #realr = np.sort(realr)
        #for j in xrange(5):
        #    for root in roots:
        #        if (realr[j] == root.real):
        #            imagr[j] = root.imag
        #print realr,imagr
        #print roots
        for j in xrange(5):
            wr[i,j] = rootss[j].real
            wi[i,j] = rootss[j].imag
            #wr[i,j] = realr[j]
            #wi[i,j] = imagr[j]

    #print wr

    sb1 = plt.subplot(211)
    for i in xrange(5):
        plt.plot(sigma,wr[:,i]/ka)

    plt.xscale('log')
    sb2 = plt.subplot(212)
    for i in xrange(5):
        plt.plot(sigma,wi[:,i])
    plt.xscale('log')
    plt.yscale('log')
    plt.savefig(outfile)
    plt.close()

def main() :
    Tgas = 1e7
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
    cred = 3e8
    cspeed = 3e8
    crat = cspeed/cs
    kappa = 3.3333e-9
    r = cred/cspeed
    #r = 1.
    L = 3e10
    wavelength=L
    ka = 2*math.pi/(wavelength/L)#/(r/10)
    sigma = kappa*rho*wavelength
    print( "cs = {0}, Prat = {1}, crat = {2} sigma={3} r={4:5.3e}".format(cs,Prat, crat,sigma,r))
    print( ka, sigma)
    omega, drho, dv, dT, dp, dEr, dFr = eigenvec(ka=ka,crat=crat,prat=Prat,sigma=sigma,gamma=gamma,r=r)
    print( drho.imag)

if( __name__ == "__main__") :
    main()
