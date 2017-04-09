#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Analyzes the results from the glass damping force test.  The test is just
the spherical collapse test from ChaNGa but with every particle given a large
velocity in the z direction.  The mean z velocity should damp exponentially
with the inverse timescale given by the runtime parameter dGlassDamper

Created on Wed Jan 25 15:06:15 2017

@author: ibackus
"""
import numpy as np
import matplotlib.pyplot as plt
import pynbody
SimArray = pynbody.array.SimArray
import diskpy

paramname = 'sphere.param'

# Parse param file
param = diskpy.utils.configparser(paramname, 'param')
dampingTime = 1.0/param['dGlassDamper']
fprefix = param['achOutName']
dDelta = param['dDelta']
# Load simulation outputs
sim = [pynbody.load(f, paramfile=paramname) \
       for f in diskpy.pychanga.get_fnames(fprefix)]
# Analyze
vel = np.array([f['vel'].mean(0) for f in sim])
v0 = vel[0, 2]
t = dDelta * np.arange(len(sim))
vpred = v0 * np.exp(-t/dampingTime)
# Plot the results
plt.clf()
for i in range(3):
    
    plt.plot(t, vel[:,i], 'o')
    
plt.plot(t, vpred, 'k')
plt.legend(['vx', 'vy', 'vz', 'vPred'], loc='best')
plt.title('Glass damping force test')
plt.ylabel('Mean particle velocity')
plt.xlabel('time')
plt.savefig('results.png')
plt.show()
