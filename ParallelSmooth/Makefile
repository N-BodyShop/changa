OPT = -Wall

HOST = /host/hermes/home/gwl
CHARMC = $(HOST)/projects/charm/bin/charmc
STRUCTURES_PATH = $(HOST)/projects/structures

INCLUDES = -I$(STRUCTURES_PATH)
CXXFLAGS = $(OPT) $(INCLUDES)
CXX = $(CHARMC)

TARGET = ParallelSmooth
VERSION = 0.2
all:: $(TARGET)

$(TARGET): Makefile

NChilada.o: TreePiece.h Sorter.h DataManager.h Tree.h
NChilada.o DataManager.o TreePiece.o Sorter.o: NChilada.h NChilada.decl.h Particle.h
TreePiece.o: DataManager.h Tree.h
NCHILADA_OBJS = NChilada.o DataManager.o TreePiece.o Sorter.o $(STRUCTURES_PATH)/TipsyFile.o

NChilada.decl.h: NChilada.ci
	$(CHARMC) $<
	

ParallelSmooth.decl.h: NChilada.decl.h
Smooth_TreePiece.o ParallelSmooth.o: ParallelSmooth.decl.h Particle.h
Smooth_TreePiece.o: TreePiece.h DataManager.h RequestResponse.h SPH_Requests.h
ParallelSmooth.o: Smooth_TreePiece.h
SPH_Requests.o: SPH_Requests.h RequestResponse.h

OBJS = $(NCHILADA_OBJS) Smooth_TreePiece.o $(TARGET).o RequestResponse.o SPH_Requests.o


$(TARGET): $(OBJS)
	$(CHARMC) $(CXXFLAGS) -o $@ $(OBJS) -language charm++ -lpopt

$(TARGET).decl.h: $(TARGET).ci
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
