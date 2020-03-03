# simple test of integrator
# run ChaNGa for 10 timesteps of .1.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

jsrun -n2 -a1 -c1 -g1 -K1 -r2 ../ChaNGa -v 1 test_pg.param ++ppn 1

# pkdgrav test_pkd.param
