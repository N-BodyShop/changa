STRUCTURES_PATH = ../structures
CXXFLAGS +=   -Wall -I$(STRUCTURES_PATH) -I../GenericTree  -O3 #-I../CacheManager
LDFLAGS += -language charm++ -module commlib -tracemode projections #-memory paranoid 
LDLIBS += -lpopt $(STRUCTURES_PATH)/SFC.o

CHARMC = $(HOME)/installations/checkedout/charm/bin/charmc
CXX = $(CHARMC)
CC = $(CXX)

TARGET = ParallelGravity
VERSION = 0.1
all: $(TARGET)

$(TARGET).o: $(TARGET).decl.h CacheManager.decl.h

$(TARGET): TreePiece.o CacheManager.o

CacheManager.o: CacheManager.C CacheManager.h CacheManager.decl.h

TreePiece.o ParallelGravity.o: ParallelGravity.h TreeNode.h CacheManager.h

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
