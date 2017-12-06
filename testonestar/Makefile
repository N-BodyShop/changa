CXXFLAGS = -g -O0 -DSNIA_TST -I.. -I../../utility/structures -DCOOLING_NONE -Wall
CXX=../../charm/bin/charmc
CHARMC=$(CXX)
LDFLAGS = -g -language charm++ ../romberg.o

default:
	run.sh

imf: imf.o
	$(CXX) $(LDFLAGS) imf.o -o imf

imf.o: imf.cpp ../imf.h testimf.decl.h

supernovaia: supernovaia.o imf.o starlifetime.o supernova.o
	$(CXX) $(LDFLAGS) supernovaia.o imf.o starlifetime.o supernova.o -o supernovaia

supernovaia.o: supernovaia.cpp ../imf.h testsnia.decl.h

%.decl.h %.def.h : %.ci
	$(CHARMC) $<
