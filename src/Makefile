# Makefile
# Project maintenance for saucy
#
# by Paul T. Darga <pdarga@umich.edu>
# and Mark Liffiton <liffiton@umich.edu>
# and Hadi Katebi <hadik@umich.edu>
#
# Changes by Jonathan Schrock UTK <jschroc1@utk.edu>
# Made to work with edge-weighted graphs and CoinUtils
#
# Copyright (C) 2004, The Regents of the University of Michigan
# See the LICENSE file for details.

IDIR=../inc
CoinInc=/Users/jschrock/CoinUtils/build/include/coin
CoinLib=/Users/jschrock/CoinUtils/build/lib
CC=g++
LDLIBS=-lz -lbz2 -lm -lCoinUtils
CFLAGS=-ansi -pedantic -Wall -Wno-write-strings -O3 -I$(IDIR) -I$(CoinInc) -L$(CoinLib)

ODIR=obj
BIN=../bin

_DEPS = saucy.h lp_amorph.h util.h platform.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))

_OBJ = util.o lp_amorph.o saucy.o P2Graph.o
OBJ = $(patsubst %,$(ODIR)/%,$(_OBJ))

_OBJ2-5 = util.o lp_amorph.o saucy2-5.o P2Graph.o
OBJ2-5 = $(patsubst %,$(ODIR)/%,$(_OBJ2-5))


$(ODIR)/%.o: %.cpp $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

all: saucy saucy2-5 

#Build saucy with weighted edges and refinement technique 2:
#   count number of different weights between vertex and refining cell
saucy: $(OBJ)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LDLIBS)

saucy2-5: $(OBJ2-5)
	$(CC) -o $(BIN)/$@ $^ $(CFLAGS) $(LDLIBS)

.PHONY: all clean

clean:
	rm -f $(ODIR)/*.o $(BIN)/* 
