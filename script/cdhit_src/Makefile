CC = g++ -Wall -ggdb
CC = g++ -pg
CC = g++

# default with OpenMP
# with OpenMP
# in command line: 
# make openmp=yes
ifeq ($(openmp),no)
  CCFLAGS = -DNO_OPENMP
else
  CCFLAGS = -fopenmp
endif

LDFLAGS = -static
#LDFLAGS = /usr/lib/x86_64-linux-gnu/libz.a -o

# default with zlib
# without zlib
# in command line:
# make zlib=no
ifeq ($(zlib),no)
  CCFLAGS += 
  LDFLAGS += -o
else
  CCFLAGS += -DWITH_ZLIB
  LDFLAGS += -lz -o
endif

# support debugging
# in command line:
# make debug=yes
# make openmp=yes debug=yes
ifeq ($(debug),yes)
CCFLAGS += -ggdb
else
CCFLAGS += -O2
endif

ifdef MAX_SEQ
CCFLAGS += -DMAX_SEQ=$(MAX_SEQ)
endif

PROGS = cd-hit-est

# Propagate hardening flags
CCFLAGS := $(CPPFLAGS) $(CCFLAGS) $(CXXFLAGS)

.c++.o:
	$(CC) $(CCFLAGS) -c $<

all: $(PROGS)

clean:
	rm -f *.o $(PROGS)

# programs

cd-hit-est: cdhit-common.o cdhit-utility.o cdhit-est.o
	$(CC) $(CCFLAGS) cdhit-est.o cdhit-common.o cdhit-utility.o $(LDFLAGS) cd-hit-est

# objects
cdhit-common.o: cdhit-common.c++ cdhit-common.h
	$(CC) $(CCFLAGS) cdhit-common.c++ -c

cdhit-utility.o: cdhit-utility.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-utility.c++ -c

cdhit.o: cdhit.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit.c++ -c

cdhit-est.o: cdhit-est.c++ cdhit-utility.h
	$(CC) $(CCFLAGS) cdhit-est.c++ -c

PREFIX ?= ../script

install:
	for prog in $(PROGS); do \
		install -m 0755 $$prog $(PREFIX); \
	done
