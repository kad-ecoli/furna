CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2fasta cif2chain fasta2miss subset_GOA obo2csv

cif2fasta: cif2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2chain: cif2chain.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2miss: fasta2miss.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

subset_GOA: subset_GOA.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

obo2csv: obo2csv.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
