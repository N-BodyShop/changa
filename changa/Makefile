STRUCTURES_PATH = ../structures
CXXFLAGS += -g -O3 -Wall -I$(STRUCTURES_PATH) -I../GenericTree -I../CacheManager
LDFLAGS += -language charm++ -module commlib #-memory paranoid #-tracemode projections
LDLIBS += -lpopt $(STRUCTURES_PATH)/SFC.o

CHARMC = $(HOME)/Projects/charm/bin/charmc
CXX = $(CHARMC)
CC = $(CXX)

TARGET = ParallelGravity
VERSION = 0.1
all:: $(TARGET)

$(TARGET).o: $(TARGET).decl.h

$(TARGET): TreePiece.o

TreePiece.o ParallelGravity.o: ParallelGravity.h TreeNode.h

%.decl.h %.def.h : %.ci
	$(CHARMC) $<

docs::
	doxygen $(TARGET).doxygen

dist::cd
	mkdir $(TARGET)-$(VERSION)
	cp Makefile $(TARGET).doxygen *.h *.cpp *.ci $(TARGET)-$(VERSION)/
	tar zcf $(TARGET)-$(VERSION).tar.gz $(TARGET)-$(VERSION)
	rm -Rf $(TARGET)-$(VERSION)

clean::
	rm -f core* *.o *~ $(TARGET) *.decl.h *.def.h charmrun conv-host 
