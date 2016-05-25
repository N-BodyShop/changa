test: cha.prof
	@echo comparing the contents of cha.prof with gas.prof
	python compare.py	

cha.prof: collapse.000030
	tipsy -nodisp < tipsy.mac

collapse.000030:
	../charmrun +p 2 ../ChaNGa collapse.param

clean:
	rm -rf cha.prof collapse.0000* *.chk?  .tipsy_history
