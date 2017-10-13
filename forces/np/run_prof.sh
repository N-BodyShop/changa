#!/bin/bash
#SBATCH -J Cambridge           # job name
#SBATCH -o CAM_result.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -N 8
#SBATCH -p gpu     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:1:00        # run time (hh:mm:ss) - 1.5 hours

../../charmrun +p2 ++mpiexec ++remote-shell "ibrun -o 0" ../../ChaNGa -v 1 lambs_30K.param >& DIAG
