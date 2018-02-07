#!/bin/sh
#SBATCH -J Cambridge           # job name
#SBATCH -o CAM_result.o%j       # output and error file name (%j expands to jobID)
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH -p gpu     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:1:00        # run time (hh:mm:ss) - 1.5 hours

# simple test of integrator
# run ChaNGa for 10 timesteps of .1.  This is roughly
# one dynamical time in the center of the cluster.
# compare this with numbers in "pkdtest.log"

../charmrun +p 8 ++mpiexec ++remote-shell "ibrun -o 0" ../ChaNGa -v 1 test_pg.param

# pkdgrav test_pkd.param
