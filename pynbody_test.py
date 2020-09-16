import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import cython
import pytest
import scipy
import astropy
import astropy.units as u
import gala.dynamics as gd
from gala.units import UnitSystem
from gala.units import galactic
import gala.integrate as gi
import superfreq
from superfreq.core import SuperFreqResult
from collections import Counter
from timeit import default_timer as timer

import pynbody

#qy = 0.9
#qz = 0.7
#core = np.sqrt(0.1)

# Define a grid of initial conditions for logarithmic halo potential
#num_xyz_values = 120
#num_range = np.linspace(0.000001,np.sqrt(np.e-core**2)*max(1,qy)-0.000001,num_xyz_values) # Start at zero to avoid double-counting
#x0 = np.repeat(num_range,num_xyz_values)
#y0 = np.tile(num_range,num_xyz_values)
#x0y0 = np.transpose([x0,y0])
#x0y0 = np.array([p for p in x0y0 if p[0]**2 + (p[1]/qy)**2 < np.e-core**2]) # Only consider points inside the ellipse in the xy-plane
#x0 = np.array([p[0] for p in x0y0])
#y0 = np.array([p[1] for p in x0y0])
#z0 = np.abs(np.sqrt(np.round(np.e-core**2 - x0**2 - (y0/qy)**2,14))*qz)

#f = pynbody.snapshot.new(dm=10122)
#f['mass'] += 0.0000000001
#f['pos'] += np.transpose([x0,y0,z0])
#print(f['mass'])
#print(f['pos'])
#print(f['vel'])
#f.properties['a'] = 0.0
#f.write(filename="zero_vel_conditions.tbin",fmt=pynbody.tipsy.TipsySnap)

nSteps = 100000 # Number of time-steps to look at (= number of simsnap files produced)
nParticles = 10122 # Number of particles to follow

x = np.empty((int(nSteps/10),nParticles))
y = np.empty((int(nSteps/10),nParticles))
z = np.empty((int(nSteps/10),nParticles))
vx = np.empty((int(nSteps/10),nParticles))
vy = np.empty((int(nSteps/10),nParticles))
vz = np.empty((int(nSteps/10),nParticles))
start = timer()
for i in range(int(nSteps/10)):
    n = str(10*(i+1)) # Corresponds to iOutInterval = 10
    # List of simsnaps (one item per unit of time):
    snapshot = pynbody.load('/phys/users/jesse1a/Desktop/Python/integrated_orbits.'+n.zfill(6))
    x[i] = snapshot['x']
    y[i] = snapshot['y']
    z[i] = snapshot['z']
    vx[i] = snapshot['vx']
    vy[i] = snapshot['vy']
    vz[i] = snapshot['vz']
    print('Time to run time step',i,'(in sec):',timer()-start)
# Make each row of position/velocity matrix correspond to time series for single particle
x = np.transpose(x)
y = np.transpose(y)
z = np.transpose(z)
vx = np.transpose(vx)
vy = np.transpose(vy)
vz = np.transpose(vz)

# Construct a grid of frequency ratios (general code applies for any potential used)
dDelta = 0.01 # Amount of time in a single timestep
T = np.linspace(0,dDelta*nSteps,int(nSteps/10))
freq_class = superfreq.SuperFreq(T,keep_calm=True)
fy_ratio0 = np.zeros(nParticles)
fz_ratio0 = np.zeros(nParticles)
start = timer()
for k in range(nParticles):
    #fig,(ax1,ax2,ax3) = plt.subplots(1,3)
    #ax1.plot(x[k],y[k])
    #ax1.set_title('y vs. x')
    #ax2.plot(x[k],z[k])
    #ax2.set_title('z vs. x')
    #ax3.plot(y[k],z[k])
    #ax3.set_title('z vs. y')
    #plt.show()
    
    xfreqs,xamps,xphis = freq_class.frecoder(x[k] + 1j*vx[k])
    yfreqs,yamps,yphis = freq_class.frecoder(y[k] + 1j*vy[k])
    zfreqs,zamps,zphis = freq_class.frecoder(z[k] + 1j*vz[k])
    if k == 0:
        print('Time elapsed (in sec) to obtain frequencies of first particle`s orbit =',timer()-start)
    fy_ratio0[k] = yfreqs[yamps.argmax()]/xfreqs[xamps.argmax()]
    fz_ratio0[k] = zfreqs[zamps.argmax()]/xfreqs[xamps.argmax()]
    if np.remainder(k,500) == 0:
        print('k =',k,'Time elapsed (in min) =',(timer()-start)/60)
print('Time to run frequency analysis with p0 = 0 (in min):',(timer()-start)/60)

# Plot frequency map
print('Number of unique fy_ratio0 values:', len(Counter(fy_ratio0).keys()))
print('Number of unique fz_ratio0 values:', len(Counter(fz_ratio0).keys()))
plt.scatter(fy_ratio0,fz_ratio0,s=1)
plt.tight_layout()
plt.show()
