# simple test of integrator
# run ChaNGa for 10 timesteps of .1.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

../charmrun +p 2 ../ChaNGa -v 1 test_pg.param ++local

# pkdgrav test_pkd.param
