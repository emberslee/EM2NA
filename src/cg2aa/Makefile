CC=g++
CFLAGS=-O3
LDFLAGS=-static

cg2aa: Arena.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp PDBParser.hpp pstream.h GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp AtomicClashes.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

cg2aa_dna: Arena_dna.cpp MissingDNAatom.hpp Superpose.hpp IdealRNA.hpp PDBParser.hpp pstream.h GeometryTools.hpp cssr.hpp BondLengths.hpp BaseConformation.hpp AtomicClashes.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

Arena_counter: Arena_counter.cpp MissingRNAatom.hpp Superpose.hpp IdealRNA.hpp PDBParser.hpp pstream.h GeometryTools.hpp cssr.hpp AtomicClashes_counter.hpp
	${CC} ${CFLAGS} $@.cpp -o $@ ${LDFLAGS}

all: cg2aa cg2aa_dna

clean:
	rm cg2aa cg2aa_dna -f
