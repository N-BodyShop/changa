import hydro
import numpy as np
import matplotlib
import math
import matplotlib.pyplot as pl
import subprocess
import multiprocessing
import subprocess
import io

delta=hydro.delta

def work(sigma):
    cmd = "python hydro.py {0} --quiet --no_plots".format(sigma)
    completedProcess = subprocess.run(cmd, shell=True, capture_output=True,text=True)
    data = np.loadtxt( io.StringIO(completedProcess.stdout))
    return sigma, data

def plotDisp( lgSigmaLow=-2., lgSigmaHigh=3, data = None) : 

    sigmaArray = []
    omegaArray = []
    for lgSigma in np.arange( lgSigmaLow, lgSigmaHigh, 0.1) : 
        sigma = 1e1**lgSigma
        hydro.kappa = sigma/(hydro.rho0*hydro.wavelength)
        omega, drho, dv, dT, dp, dEr, dFr = hydro.findEigen() 
        sigmaArray.append(sigma)
        omegaArray.append(omega)
        #print(sigma, omega.real, omega.imag)
            
    matplotlib.rc('text', usetex=True)
    f, (ax1, ax2) = pl.subplots(2, sharex=True)
    #ax2.semilogx( sigmaArray. np.real(np.array( omegaArray)))
    ax2.loglog( sigmaArray, np.imag(np.array( omegaArray)))
    ax1.semilogx( sigmaArray, np.real(np.array( omegaArray))/(2.*math.pi))
    if not data is None : 
        sigmaData = data[:,0]
        vph = data[:,1]
        omegai = np.log(np.abs(data[:,2])/delta)
        print(sigmaData, vph, omegai)
        ax2.scatter( sigmaData, -omegai)
        ax1.scatter( sigmaData, vph)

    pl.xlim(0.01,1000.)
    pl.ylim(1e-3,1.0)
    pl.setp([a.get_xticklabels() for a in f.axes[:-1]], visible=False)
    f.subplots_adjust(hspace=0)


    pl.savefig("disp.pdf")

    return 

#plotDisp(-2., 3)
#exit()
count = 4#multiprocessing.cpu_count()
pool = multiprocessing.Pool(processes=count)
lgSigmas = np.arange(-3, 3, 0.5)
results = pool.map(work, 1e1**lgSigmas)

data = []
for result in results :
    sigma = result[0]
    dataArray = result[1]
    amp = dataArray[-1,1]
    vph = dataArray[-1,2]
    data.append([sigma, vph, amp])

plotDisp(-2., 3, np.array(data))

