SOURCES=
CSOURCES=athena_util.c 

ifdef MHD
	CSOURCES+=hlld.c hlle.c
else
  	CSOURCES+=hllc.c hlle.c
endif

INCLUDES=athena_defs.h copyright.h
OBJECTS=$(SOURCES:.cc=.o )
COBJECTS=$(CSOURCES:.c=.o )
all: libathena.a 

driver: $(OBJECTS) $(COBJECTS) $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -o driver $(OBJECTS) $(COBJECTS) $(HDF5LIBS) $(EXTRALIBS) 


$(OBJECTS): %.o: %.cc $(INCLUDES)
	$(CXX) $(CXXFLAGS) $(HDF5INCS) -c $< -o $@

$(COBJECTS): %.o: %.c $(INCLUDES)
	$(CC) $(CFLAGS) $(HDF5INCS) -c $< -o $@

libathena.a: $(OBJECTS) $(COBJECTS)
	ar r libathena.a $(COBJECTS) $(OBJECTS)

clean:
	rm -f *.o *.a 
