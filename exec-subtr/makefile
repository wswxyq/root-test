include ../Makefile.common
VPATH = $(DBIN):$(DOBJ)/subtr:$(DEXE)
CC=g++ #compiler:g++
CFLAGS=-c -g -Wall `root-config --cflags` #give the path to c header
LDFLAGS=`root-config --glibs` #optimizing parameters

LNKOPTS=$(LDFLAGS) $(DOBJ)/subtr/*.o -o $(DBIN)/$@ 
ECOPTS=$(CFLAGS) $< -o $(DOBJ)/subtr/$@


EXECUTABLE=	subtr 

all: $(EXECUTABLE)

subtr: subtr.o
	$(CC) $(LNKOPTS) 
subtr.o: subtr.cc subtr_exec.h subtr_common.h
	$(CC) $(ECOPTS) 


