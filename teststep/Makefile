test: testenergy testrestart

testenergy: pgtest_ms.log
	@echo -n 'energy should be: '
	@tail -1 pkdtest.log | awk '{print $$3;}'
	@echo -n 'energy is: '
	@tail -1 pgtest_ms.log | awk '{print $$4;}'
	
pgtest_ms.log:
	sh runtest.sh #../charmrun +p 2 ../ChaNGa test_pg.param
	
testrestart:
	../charmrun +p 2 ../ChaNGa -dt 0.001 -binout 6 -v 1 test_pg.param
	../charmrun +p 2 ../ChaNGa -n 15 +restart pgtest_ms.chk0
clean:
	rm -rf *.chk? pgtest_ms.log pgtest_ms.000* lastcheckpoint
