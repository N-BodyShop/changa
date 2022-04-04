SROEOS_DIR = ../../../other-codes/sroeos
HDF5LIBS = -L$(HOME)/Code/local/lib -lhdf5
HDF5INCS = -I$(HOME)/Code/local/include
ifdef TACC_HDF5_DIR
	HDF5LIBS=-L$(TACC_HDF5_LIB) -lhdf5
	HDF5INCS=-I$(TACC_HDF5_INC)
endif

CFLAGS=-g -I/usr/include/malloc

all: libnuc_eos.so

libnuc_eos.so:
	cd $(SROEOS_DIR); make clean; make HDF5INCS="$(HDF5INCS)" HDF5LIBS="$(HDF5LIBS)" libnuc_eos.so; cp libnuc_eos.so $(PWD)

clean:
	rm -f $(TOOLSDIR)/*.o
	rm -f $(TOOLS)
	rm -f $(TOOLS)
