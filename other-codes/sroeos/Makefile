
include make.inc

SOURCES=driver.cc \
        nuc_eos_press.cc \
        readtable.cc \
        nuc_eos_press_cs2.cc \
        nuc_eos_dpdrhoe_dpderho.cc \
        nuc_eos_short.cc \
        nuc_eos_full.cc \
        nuc_eos_for_c.cc
CSOURCES= 
INCLUDES=nuc_eos.hh helpers.hh 
OBJECTS=$(SOURCES:.cc=.o )
COBJECTS=$(CSOURCES:.c=.o )

all: driver libnuc_eos.a 

driver: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o driver $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) 


$(F90OBJ): %.o: %.f90
	$(F90) -c $< -o $@

$(OBJECTS): %.o: %.cc $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(HDF5INCS) -c $< -o $@

$(COBJECTS): %.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) $(HDF5INCS) -c $< -o $@

libnuc_eos.a: $(OBJECTS)
	ar r libnuc_eos.a nuc_eos*.o readtable.o

libnuc_eos.so: $(OBJECTS)
	$(CXX) -shared -o libnuc_eos.so nuc_eos*.o readtable.o $(HDF5LIBS) $(EXTRALIBS)
clean:
	rm -f *.o
	rm -rf *.cc.*
	rm -f *.a *.so
	rm -f driver test
