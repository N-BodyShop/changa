# Flags that can be passed during compilation
#FLAG_PRINT = -DCOSMO_PRINT=1       #enable really verbose prints
FLAG_STATISTICS = -DCOSMO_STATS=1  # enable statistics collection
FLAG_DEBUG = -DCOSMO_DEBUG=2       # enable debugging information
CACHE_TREE = -DCACHE_TREE          # enable processor tree inside the cache
#INTERLIST = -DINTERLIST_VER=2      # enable interaction lists
DEFINE_FLAGS = $(FLAG_PRINT) $(FLAG_STATISTICS) $(FLAG_DEBUG) $(CACHE_TREE) $(INTERLIST)

STRUCTURES_PATH = ../structures
CXXFLAGS += -O3 -I$(STRUCTURES_PATH) -I../ParallelGravity -I../GenericTree  -Wall $(DEFINE_FLAGS)
LDFLAGS += -language charm++ -module commlib -module RefineLB -module RefineCommLB -module GreedyLB -module GreedyCommLB -module OrbLB -module RotateLB #-tracemode projections #-module EveryLB #-memory paranoid
LDLIBS += ../structures/libTipsy.a

CHARMC = ../../charm/bin/charmc $(OPTS)
CXX = $(CHARMC)
CC = $(CXX)
CXX_DEPEND = $(CXX) -M -MM -MG $(CXXFLAGS)

OBJECTS = DataManager.o Reductions.o TreePiece.o CacheManager.o Sorter.o \
	  param.o ../GenericTree/GenericTreeNode.o ParallelGravity.o
SRSC = Reductions.cpp DataManager.cpp Sorter.cpp TreePiece.cpp CacheManager.C \
	param.c ../GenericTree/GenericTreeNode.C ParallelGravity.cpp

TARGET = ParallelGravity
VERSION = 0.2
all: $(TARGET)

$(TARGET): $(OBJECTS)

# $(TARGET).o: $(TARGET).decl.h CacheManager.decl.h Makefile

# $(TARGET): TreePiece.o CacheManager.o Sorter.o DataManager.o Reductions.o \
# 	param.o ../GenericTree/GenericTreeNode.o

# ParallelGravity.o: ParallelGravity.h param.h parameters.h

# CacheManager.o: CacheManager.C CacheManager.h CacheManager.decl.h\
# 	../GenericTree/MultipoleMoments.h ParallelGravity.h Makefile

# TreePiece.o ParallelGravity.o: ParallelGravity.h TreeNode.h CacheManager.h\
# 	Reductions.decl.h\
# 	../structures/tree_xdr.h\
# 	../GenericTree/MultipoleMoments.h ../GenericTree/GenericTreeNode.h Makefile
# Reductions.o: Reductions.h Reductions.decl.h

# ../GenericTree/GenericTreeNode.o: ../GenericTree/GenericTreeNode.C ../GenericTree/GenericTreeNode.h Makefile

%.decl.h %.def.h : %.ci
	$(CHARMC) $<

%.o: Makefile

docs:
	doxygen $(TARGET).doxygen

dist:
	mkdir $(TARGET)-$(VERSION)
	cp Makefile $(TARGET).doxygen *.h *.cpp *.ci $(TARGET)-$(VERSION)/
	tar zcf $(TARGET)-$(VERSION).tar.gz $(TARGET)-$(VERSION)
	rm -Rf $(TARGET)-$(VERSION)

clean:
	rm -f core* $(OBJECTS) *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 

# depend:
# 	$(CXX_DEPEND) $(SRSC) > Makefile.dep

# The following line is a script usable to regenerate the dependace file,
# without the inclusion of charm headers.
#../../charm/bin/charmc  -M -MM -MG -O3 -I../structures -I../ParallelGravity -I../GenericTree  -Wall  -DCOSMO_STATS=1   -DCOSMO_DEBUG=2        -DCACHE_TREE            Reductions.cpp DataManager.cpp Sorter.cpp TreePiece.cpp CacheManager.C param.c ../GenericTree/GenericTreeNode.C ParallelGravity.cpp|while read i;do echo $i| awk -F' ' '{for (i=1;i<NF;++i) print $i" \\"}';echo;done|grep -v "charm/bin" > Makefile.dep


.PHONY: all docs dist clean depend

include Makefile.dep
