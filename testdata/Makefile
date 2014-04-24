# Change the N here to make a new size dataset.
N=1000

testdata.tbin: ppartt Makefile
	./ppartt $N -.5 .5 -1 > testdata.tbin

testdata: ppartnc Makefile
	./ppartnc $N -.5 .5 -1 testdata

ppartt: ppartt.c tipsydefs.h
	$(CC) -o ppartt ppartt.c -lm

pgaspart: pgaspart.c tipsydefs.h
	$(CC) -o pgaspart pgaspart.c -lm

ppartnc: ppartnc.C
	$(CXX) -I../../utility/structures -o ppartnc ppartnc.C -lm

clean:
	-rm -rf pgaspart ppartt ppartnc testdata.tbin testdata

