test: pgtest_ms.log
	@echo -n 'energy should be: '
	@tail -1 pkdtest.log | awk '{print $$3;}'
	@echo -n 'energy is: '
	@tail -1 pgtest_ms.log | awk '{print $$4;}'
	@tail -1 pgtest_ms.log | awk '{if(sqrt(($$4 + 32.19)^2) < 0.005) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'
	
	
pgtest_ms.log:
	sh runtest.sh #../charmrun +p 2 ../ChaNGa test_pg.param
	
clean:
	rm -rf *.chk? pgtest_ms.log pgtest_ms.000* lastcheckpoint
