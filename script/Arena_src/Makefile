CC=clang++
CFLAGS=-O3
LDFLAGS=-static

Arena: Arena.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp PDBParser.hpp pstream.h GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp AtomicClashes.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

install: Arena
	cp Arena ../Arena
