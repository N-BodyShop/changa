# simple test of integrator
# check for lines with "Energy" in the output
# run ChaNGa for 100 timesteps of .01.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

../charmrun +p 2 ../ChaNGa -v 1 -p 2 test_pg.param

# pkdgrav test_pkd.param
