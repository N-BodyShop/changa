ifndef MESA_DIR
$(error MESA_DIR must be set)
endif

LDOPTS = -L$(MESA_DIR)/lib -leos -lchem -linterp_2d -linterp_1d -lnum -lf2crlibm -lcrlibm -lmtx -lconst -lutils -lmesalapack -lmesablas -lm

CFLAGS=-g
FC=gfortran

all: hydro

hydro : hydro.f90
	$(FC) -fopenmp -O2 -o hydro hydro.f90 -I$(MESA_DIR)/include $(LDOPTS)

clean:
	rm -f *.o
	rm -f hydro
