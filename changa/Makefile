# Flags that can be passed during compilation
#FLAG_PRINT = -DCOSMO_PRINT=1       #enable really verbose prints
FLAG_STATISTICS = -DCOSMO_STATS=1  # enable statistics collection
#FLAG_DEBUG = -DCOSMO_DEBUG=2       # enable debugging information
DEFINE_FLAGS = $(FLAG_PRINT) $(FLAG_STATISTICS) $(FLAG_DEBUG)

STRUCTURES_PATH = ../structures
CXXFLAGS +=  -g -I$(STRUCTURES_PATH) -I. -I../GenericTree $(DEFINE_FLAGS)
LDFLAGS += -language charm++ -module commlib -module RefineLB -module RefineCommLB -module GreedyLB -module GreedyCommLB -module OrbLB -tracemode projections #-module EveryLB #-memory paranoid 
LDLIBS +=  $(STRUCTURES_PATH)/SFC.o ../GenericTree/GenericTreeNode.o ../structures/libTipsy.a

CHARMC = $(HOME)/peak/src/charm/bin/charmc $(OPTS)
CXX = $(CHARMC)
CC = $(CXX)

TARGET = ParallelGravity
VERSION = 0.2
all: $(TARGET)

$(TARGET).o: $(TARGET).decl.h CacheManager.decl.h Makefile

$(TARGET): TreePiece.o CacheManager.o Sorter.o DataManager.o Reductions.o

CacheManager.o: CacheManager.C CacheManager.h CacheManager.decl.h\
	../GenericTree/MultipoleMoments.h ParallelGravity.h Makefile

TreePiece.o ParallelGravity.o: ParallelGravity.h TreeNode.h CacheManager.h\
	Reductions.decl.h\
	../structures/tree_xdr.h\
	../GenericTree/MultipoleMoments.h ../GenericTree/GenericTreeNode.h Makefile
Reductions.o: Reductions.h Reductions.decl.h

%.decl.h %.def.h : %.ci
	$(CHARMC) $<

docs:
	doxygen $(TARGET).doxygen

dist:
	mkdir $(TARGET)-$(VERSION)
	cp Makefile $(TARGET).doxygen *.h *.cpp *.ci $(TARGET)-$(VERSION)/
	tar zcf $(TARGET)-$(VERSION).tar.gz $(TARGET)-$(VERSION)
	rm -Rf $(TARGET)-$(VERSION)

clean:
	rm -f core* *.o *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 

.PHONY: all docs dist clean
