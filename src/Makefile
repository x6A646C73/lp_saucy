# Makefile
# Project maintenance for saucy
#
# by Paul T. Darga <pdarga@umich.edu>
# and Mark Liffiton <liffiton@umich.edu>
# and Hadi Katebi <hadik@umich.edu>
#
# Copyright (C) 2004, The Regents of the University of Michigan
# See the LICENSE file for details.

COINDIR=/Users/jschrock/CoinUtils/build
IDIR =../inc
CC=g++
#CC=gcc
CFLAGS=-pedantic -Wall -O3 -I$(IDIR)
#LDFLAGS=-shared
#CFLAGS=-pedantic -Wall -O0 -g -I$(IDIR)
LDLIBS=-lz
LPLIB=-L$(COINDIR)/lib -lCoinUtils -lbz2 -lz -lm
LPINC=-I$(COINDIR)/include/coin
#TARGET_LIB=libsaucy.so

ODIR=obj
BIN=../bin

_DEPS = saucy.h amorph.h util.h platform.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

#_OBJ0 = util.o saucy0.o saucyio.o main.o
#_OBJ2-5 = util.o saucy2-5.o saucyio.o main.o
_OBJ2-5 = util.o saucy2-5.o amorph.o LP2Graph.o
_OBJPBIL = util.o saucy2-5.o amorph.o findAlmostPBIL.o
_OBJDYN = util.o saucy2-5.o amorph.o findAlmostDYN.o
#OBJ0 = $(patsubst %,$(ODIR)/%,$(_OBJ0))
OBJ2-5 = $(patsubst %,$(ODIR)/%,$(_OBJ2-5))
OBJPBIL = $(patsubst %,$(ODIR)/%,$(_OBJPBIL))
OBJDYN = $(patsubst %,$(ODIR)/%,$(_OBJDYN))

saucy2-5: $(OBJ2-5)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LDLIBS) $(LPLIB) $(LPINC)

findAlmostPBIL: $(OBJPBIL)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LDLIBS) $(LPLIB) $(LPINC)

findAlmostDYN: $(OBJDYN)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LDLIBS) $(LPLIB) $(LPINC)

$(ODIR)/util.o: util.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/saucy2-5.o: saucy2-5.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

$(ODIR)/amorph.o: amorph.cpp $(DEPS)
	$(CC) -c -o $@ $(LPLIB) $(LPINC) $< $(CFLAGS)

$(ODIR)/LP2Graph.o: LP2Graph.cpp $(DEPS)
	$(CC) -c -o $@ $(LPLIB) $(LPINC) $< $(CFLAGS)

$(ODIR)/findAlmostPBIL.o: findAlmostPBIL_saucy.cpp $(DEPS)
	$(CC) -c -o $@ $(LPLIB) $(LPINC) $< $(CFLAGS)

$(ODIR)/findAlmostDYN.o: findAlmostDYN_saucy.cpp $(DEPS)
	$(CC) -c -o $@ $(LPLIB) $(LPINC) $< $(CFLAGS)

all: saucy2-5 findAlmostPBIL findAlmostDYN

.PHONY: all clean

clean:
	rm -f $(ODIR)/*.o $(BIN)/* 
