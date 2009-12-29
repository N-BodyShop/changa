test: cha.prof
	@echo comparing the contents of cha.prof with gas.prof
	python compare.py	

cha.prof: adiabtophat_glass_28768_Std.000030
	tipsy -nodisp < tipsy.mac

adiabtophat_glass_28768_Std.000030:
	../charmrun +p 2 ../ChaNGa adiabtophat_glass_28768_Std.param

clean:
	rm -rf cha.prof adiabtophat_glass_28768_Std.000030 *.chk?  .tipsy_history
