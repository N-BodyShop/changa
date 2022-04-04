# C compiler
CC=gcc
# C++ compiler
CXX=g++
F90=gfortran -fPIC -g

# CXX compiler flags (for optimization)
CXXFLAGS=-g -fPIC -ansi -Wall -Winline -mtune=corei7 -fopenmp -msse4.2 -std=c++11 -Drestrict=__restrict__ -Ofast $(EXTRAINCS) 

#LDFLAGS=-static

# C compiler flags (for optimization)
CFLAGS=-g -fPIC -O3 -std=gnu99 -ansi -Wall -ftree-vectorizer-verbose=6  -Winline -mtune=corei7 -fopenmp -msse2

EXTRALIBS=-lm 

# path to HDF5 directory and includes
#HDF5DIR=$(HOME)/Code/local
#HDF5DIR=$(HDF5HOME)
HDF5INCS=-I$(HDF5DIR)/include
HDF5LIBS=-L$(HDF5DIR)/lib -lhdf5 -lz



