# simple test of integrator
# run ChaNGa for 10 timesteps of .1.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

../charmrun +p 8 ../changa ++local -D 3 -iPush 15000 test_pg.param

# pkdgrav test_pkd.param
