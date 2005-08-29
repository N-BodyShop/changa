CXXFLAGS = -O3 -g -I../structures -I../ParallelGravity -I.
CHARMC = ../../charm/bin/charmc $(OPTS)
CXX = $(CHARMC)
CC = $(CXX)

TARGETS = GenericTreeNode.o

.PHONY: all clean

all: $(TARGETS)

clean:
	rm -f $(TARGETS)
