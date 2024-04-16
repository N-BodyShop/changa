#run with the catalogue as an input file
#   python3 mf2.py "path to catalogue"



import sys
import numpy as np

#previously determined cumulative mass function against which the input file will be compared
test_logm = np.array([14.1744, 14.4744, 14.7744, 15.0744, 15.3744])
test_logn = np.array([-5.1217, -5.5176, -6.0164, -6.4314, -7.4314])
cube_side = 300


#define a function that takes in an array of halos and returns an array that only containts those with a mass greater than M
def greater_than_mass(M, mn_array):
    return mn_array[:,1][(mn_array[:,1]>M)]

#load in the data from the input file
catalogue=np.loadtxt(sys.argv[1])


#pull out just the mass and the number of each halo
mvsn = catalogue[:,0:3]


#find the cumulative mass function at values of the mass that match the test function
mfc=np.zeros(test_logm.size)
for i in range(test_logm.size):
    n_greater=greater_than_mass(10**test_logm[i], mvsn).size
    mfc[i]=n_greater

#convert number into log of number density
mf_units = mfc/cube_side**3
logn=np.log10(mf_units)

#find percent error from the test function
logn_pe = np.max(np.abs((logn-test_logn)/test_logn)*100)

if logn_pe<5:
    print("Congrats! This simulation produces a mass function as predicted!")
    
else:
    print("PROBLEM!  This mass function differs from predicted by more than 5%!")

print("The largest error is %.2f%%.\n" % (logn_pe))
