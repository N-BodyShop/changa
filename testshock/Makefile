test: pres.out
	@echo compare the contents of pres.out with pres_gasoline.out

pressplot: pressplot.o
	$(CC) -o pressplot pressplot.o -lm

pres.out: pressplot shocktube.000300
	./pressplot -1 1 100 < shocktube.000300 > pres.out

shocktube.000300:
	../charmrun ../ChaNGa shocktube.param

clean:
	rm -rf pressplot pressplot.o shocktube.000300 pres.out *.chk?  shocktube.log lastcheckpoint
