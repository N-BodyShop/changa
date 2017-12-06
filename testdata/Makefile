# Change the N here to make a new size dataset.
N=1000

testdata.tbin: ppartt Makefile
	./ppartt $N -.5 .5 -1 > testdata.tbin

testdata: ppartnc Makefile
	./ppartnc $N -.5 .5 -1 testdata

ppartt: ppartt.cpp tipsydefs.h
	$(CC) -o ppartt ppartt.cpp -lm

pgaspart: pgaspart.cpp tipsydefs.h
	$(CC) -o pgaspart pgaspart.cpp -lm

ppartnc: ppartnc.cpp
	$(CXX) -O3 -I../../utility/structures -o ppartnc ppartnc.cpp -lm

pgaspartnc: pgaspartnc.cpp
	$(CXX) -O3 -I../../utility/structures -o pgaspartnc pgaspartnc.cpp -lm

clean:
	-rm -rf pgaspart pgaspartnc ppartt ppartnc testdata.tbin testdata

