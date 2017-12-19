#!/bin/sh
#SBATCH -J Cambridge           # job name
#SBATCH -o CAM_result.o%j       # output and error file name (%j expands to jobID)
#SBATCH -n 32              # total number of mpi tasks requested
#SBATCH -N 8
#SBATCH -p gpu     # queue (partition) -- normal, development, etc.
#SBATCH -t 00:1:00        # run time (hh:mm:ss) - 1.5 hours

test: testenergy testrestart

testenergy: pgtest_ms.log
	@echo -n 'energy should be: '
	@tail -1 pkdtest.log | awk '{print $$3;}'
	@echo -n 'energy is: '
	@tail -1 pgtest_ms.log | awk '{print $$4;}'
	@tail -1 pgtest_ms.log | awk '{if(sqrt(($$4 + 32.19)^2) < 0.005) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'
	
	
pgtest_ms.log:
	sh runtest.sh #../charmrun +p 2 ../ChaNGa test_pg.param
	
testrestart:
	../charmrun +p 8 ++mpiexec ++remote-shell "ibrun -o 0" ../ChaNGa -dt 0.001 -binout 6 -v 1 test_pg.param
	../charmrun +p 8 ++mpiexec ++remote-shell "ibrun -o 0" ../ChaNGa -n 15 +restart pgtest_ms.chk0
clean:
	rm -rf *.chk? pgtest_ms.log pgtest_ms.000* lastcheckpoint
