# Flags that can be passed during compilation
#FLAG_PRINT = -DCOSMO_PRINT=1       #enable really verbose prints
FLAG_STATISTICS = -DCOSMO_STATS=1  # enable statistics collection
#FLAG_DEBUG = -DCOSMO_DEBUG=1       # enable debugging information
DEFINE_FLAGS = $(FLAG_PRINT) $(FLAG_STATISTICS) $(FLAG_DEBUG)

STRUCTURES_PATH = ../structures
CXXFLAGS +=  -I$(STRUCTURES_PATH) -I. -I../GenericTree -O3 -Wall $(DEFINE_FLAGS)
LDFLAGS += -language charm++ -module commlib -module RefineLB -module RefineCommLB -module MetisLB -module GreedyLB -module GreedyCommLB -module OrbLB -tracemode projections #-module EveryLB #-memory paranoid 
LDLIBS +=  $(STRUCTURES_PATH)/SFC.o ../GenericTree/GenericTreeNode.o

CHARMC = ../../charm/bin/charmc $(OPTS)
CXX = $(CHARMC)
CC = $(CXX)

TARGET = ParallelGravity
VERSION = 0.2
all: $(TARGET)

$(TARGET).o: $(TARGET).decl.h CacheManager.decl.h

$(TARGET): TreePiece.o CacheManager.o

CacheManager.o: CacheManager.C CacheManager.h CacheManager.decl.h\
	../GenericTree/MultipoleMoments.h ParallelGravity.h

TreePiece.o ParallelGravity.o: ParallelGravity.h TreeNode.h CacheManager.h\
	../GenericTree/MultipoleMoments.h ../GenericTree/GenericTreeNode.h

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
