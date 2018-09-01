include ../Makefile.common
VPATH = $(DBIN):$(DOBJ)/subtrx:$(DEXE)
CC=g++ #compiler:g++
CFLAGS=-c -g -Wall `root-config --cflags` #give the path to c header
LDFLAGS=`root-config --glibs` #optimizing parameters

LNKOPTS=$(LDFLAGS) $(DOBJ)/subtrx/*.o -o $(DBIN)/$@ 
ECOPTS=$(CFLAGS) $< -o $(DOBJ)/subtrx/$@


EXECUTABLE=	subtrx

all: $(EXECUTABLE)

subtr: subtrx.o
	$(CC) $(LNKOPTS) 
subtr.o: subtrx.cc subtr_execx.h subtr_commonx.h
	$(CC) $(ECOPTS) 


