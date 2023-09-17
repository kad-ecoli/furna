CC=g++
CFLAGS=-O3 -std=c++11
LDFLAGS=-static

all: cif2fasta cif2chain fasta2miss fasta2taxon subset_GOA obo2csv cif2pdb fasta2freq NWalign

cif2fasta: cif2fasta.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2chain: cif2chain.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2miss: fasta2miss.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2taxon: fasta2taxon.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

fasta2freq: fasta2freq.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}


subset_GOA: subset_GOA.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

obo2csv: obo2csv.cpp StringTools.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cif2pdb: cif2pdb.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

NWalign: NWalign.cpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}
