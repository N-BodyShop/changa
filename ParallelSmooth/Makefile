STRUCTURES_PATH = $(HOME)/projects/structures
CXXFLAGS += -g -O3 -Wall -I$(STRUCTURES_PATH)
LDFLAGS += -L$(STRUCTURES_PATH) -language charm++
LDLIBS += -lpopt -lTipsy

CHARMC = $(HOME)/projects/charm/bin/charmc
CXX = $(CHARMC)
CC = $(CXX)

TARGET = ParallelSmooth
VERSION = 0.2
all:: $(TARGET)

Reductions.o: Reductions.decl.h

NChilada.o: TreePiece.h Sorter.h DataManager.h Tree.h
NChilada.o DataManager.o TreePiece.o Sorter.o: NChilada.h NChilada.decl.h Particle.h
TreePiece.o: DataManager.h Tree.h
NCHILADA_OBJS = NChilada.o DataManager.o TreePiece.o Sorter.o

ParallelSmooth.decl.h: NChilada.decl.h
Smooth_TreePiece.o ParallelSmooth.o: ParallelSmooth.decl.h Particle.h
Smooth_TreePiece.o: TreePiece.h DataManager.h RequestResponse.h SPH_Requests.h
ParallelSmooth.o: Smooth_TreePiece.h Reductions.decl.h
SPH_Requests.o: SPH_Requests.h RequestResponse.h

$(TARGET): $(NCHILADA_OBJS) Smooth_TreePiece.o $(TARGET).o RequestResponse.o SPH_Requests.o Reductions.o SFC.o

%.decl.h %.def.h : %.ci
	$(CHARMC) $<

docs::
	doxygen $(TARGET).doxygen

dist::
	mkdir $(TARGET)-$(VERSION)
	cp Makefile $(TARGET).doxygen *.h *.cpp *.ci $(TARGET)-$(VERSION)/
	tar zcf $(TARGET)-$(VERSION).tar.gz $(TARGET)-$(VERSION)
	rm -Rf $(TARGET)-$(VERSION)

clean::
	rm -f core* *.o *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 
