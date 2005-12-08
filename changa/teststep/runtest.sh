# simple test of integrator
# check for lines with "Energy" in the output
# run Parallel gravity for 100 timesteps of .01.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

charmrun +p 2 ParallelGravity -v -T 0.01 -p 2 -n 100 king_soft.bin

# pkdgrav test_pkd.param
