include ../Makefile.common
VPATH = $(DBIN):$(DOBJ)/fit:$(DEXE)
CC=g++ #compiler:g++
CFLAGS=-c -g -Wall `root-config --cflags` #give the path to c header
LDFLAGS=`root-config --ldflags --glibs` -lRooFit -lRooFitCore -lMinuit -lRooStats

LNKOPTS=$(LDFLAGS) $(DOBJ)/fit/*.o -o $(DBIN)/$@ 
ECOPTS=$(CFLAGS) $< -o $(DOBJ)/fit/$@


EXECUTABLE=	fit2D_M 

all: $(EXECUTABLE)

fit2D_M: fit2D_M.o
	$(CC) $(LNKOPTS) 
fit2D_M.o: fit2D_M.cc 
	$(CC) $(ECOPTS) 




