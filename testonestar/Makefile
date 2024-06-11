CXXFLAGS = -g -O0 -DSNIA_TST -DWENDLAND -I.. -I../../utility/structures -DCOOLING_NONE -Wall -I/usr/include/tirpc
CXX=../../charm/bin/charmc
CHARMC=$(CXX)
LDFLAGS = -g -language charm++ ../romberg.o

default:
	run.sh

imf: imf.o
	$(CXX) $(LDFLAGS) imf.o -o imf

imf.o: imf.cpp ../imf.h testimf.decl.h

supernovaia: supernovaia.o imf.o starlifetime.o supernova.o
	$(CXX) $(LDFLAGS) supernovaia.o imf.o starlifetime.o supernova.o -lck -lconverse -o supernovaia

supernovaia.o: supernovaia.cpp ../imf.h testsnia.decl.h

clean:
	rm -f *.o supernovaia *.decl.h *.def.h onestar.000130* onestar.starlog
	rm -f onestar.log onestar.timings onestar.out lastcheckpoint
	rm -f SNII.out SNIa.out winds.out eff.out
	rm -rf onestar.chk?

%.decl.h %.def.h : %.ci
	$(CHARMC) $<
