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
	jsrun -n2 -a1 -c1 -g1 -K1 -r2 ../ChaNGa -dt 0.001 -binout 6 -v 1 test_pg.param ++ppn 1
	jsrun -n2 -a1 -c1 -g1 -K1 -r2 ../ChaNGa -n 15 +restart pgtest_ms.chk0 ++ppn 1
clean:
	rm -rf *.chk? pgtest_ms.log pgtest_ms.000* lastcheckpoint
