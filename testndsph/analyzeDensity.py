#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 17 16:03:40 2017

@author: ibackus
"""

import numpy as np
import sys

if __name__ == "__main__":
    
    fprefix = sys.argv[1]
    fname = fprefix + ".000000.gasden"
    rho_array = np.genfromtxt(fname, skip_header=1)
    
    print "Target density: 1"
    print "Average density:", rho_array.mean()