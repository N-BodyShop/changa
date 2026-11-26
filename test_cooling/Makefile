test: testh2 testmetal testcosmo testboley

testh2: h2.log
	@echo -n 'energy should be: '
	@tail -1 h2_master.log | awk '{print $$9;}'
	@echo -n 'energy is: '
	@tail -1 h2.log | awk '{print $$9;}'
	@tail -1 h2.log | awk '{if(sqrt(($$9 - 1.3476e-10)^2) < 1e-24) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'

testmetal: metal.log
	@echo -n 'energy should be: '
	@tail -1 metal_master.log | awk '{print $$9;}'
	@echo -n 'energy is: '
	@tail -1 metal.log | awk '{print $$9;}'
	@tail -1 metal.log | awk '{if(sqrt(($$9 - 1.34749e-10)^2) < 1e-24) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'

testcosmo: cosmo.log
	@echo -n 'energy should be: '
	@tail -1 cosmo_master.log | awk '{print $$9;}'
	@echo -n 'energy is: '
	@tail -1 cosmo.log | awk '{print $$9;}'
	@tail -1 cosmo.log | awk '{if(sqrt(($$9 - 1.34225e-10)^2) < 1e-24) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'

testboley: boley.log
	@echo -n 'energy should be: '
	@tail -1 boley_master.log | awk '{print $$5;}'
	@echo -n 'energy is: '
	@tail -1 boley.log | awk '{print $$5;}'
	@tail -1 boley.log | awk '{if(sqrt(($$5 - 4.88461e-23)^2) < 1e-24) { print "TEST PASSED"; exit 0;} else {print "TEST FAILED"; exit -1;}}'
