CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2fasta cif2chain

cif2fasta: cif2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2chain: cif2chain.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
