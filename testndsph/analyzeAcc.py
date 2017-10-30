#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  2 12:44:44 2017

@author: ibackus
"""
import matplotlib.pyplot as plt
import numpy as np
import sys
import pynbody

def load_acc(fname):
    """
    Loads the 3D acceleration auxiliary array (for a tipsy snap)
    
    fname is the snapshot filename (i.e. leaving out the .acc2 extension)
    """
    acc_name = fname + '.acc2'
    acc0 = np.genfromtxt(acc_name)
    nParticles = int(acc0[0])
    acc = acc0[1:].reshape([nParticles, 3], order='F')
    
    if np.any(np.isnan(acc)):
        
        raise ValueError, "NaN encountered in acceleration file " + acc_name
        
    return acc

if __name__ == '__main__':
    
    # Settings
#    fname = 'acceleration_1d/cube_1d_nx1000.000000'
#    if len(sys.argv) > 1:
#        fname = sys.argv[1] + '.000000'
    fname = sys.argv[1] + '.000000'
    sigma = float(sys.argv[2])
    figname = 'acceleration_test_results.png'
#    sigma = 0.05
    cs = 1.
    
    # Load
    f = pynbody.load(fname)
    r = f['r']
    rho = f['rho']
    acc = load_acc(fname)
    a = np.sqrt((acc**2).sum(1))
    # Strip units for simplicity
    r.units = '1'
    rho.units = '1'
    
    # Calculate
    # We estimate the background density to focus on density gradient estimates.
    # We don't want to be sensitive to estimating the density wrong by some
    # fixed factor
    mask = (r > 5 * sigma)
    rho0 = f['rho'][mask].mean()
    rho0.units = '1'
    apred = cs**2 * r * (rho - rho0)/rho/sigma**2
    
    # Plot results
    plt.clf()
    plt.plot(r, a, 'o', r, apred, '.')
    plt.legend(['actual', 'predicted'])
    plt.xlabel('r')
    plt.ylabel('acceleration')
    plt.title('Acceleration (gradient estimate) test')
    plt.draw()
    plt.savefig(figname)
    print 'Test results saved to', figname
    plt.show()
