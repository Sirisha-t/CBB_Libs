C++      = g++
CFLAGS   = -std=c++11 -Wall -Wno-deprecated -fopenmp
OPT      = -O2

LIBPATH = -L/usr/local/lib -L/usr/local/lib64
LIBS     = -lrt -lm -lz -lgomp -lpthread -lboost_iostreams -lboost_program_options -lboost_filesystem -lboost_system



OBJS    = main_index.o seq_index.o loader.o bio_alphabet.o
EXE     = ../bin/mm_index


.cpp.o:
	$(C++) $(CFLAGS) $(OPT) -c $<

install:  mm_index

mm_index: $(OBJS)
	$(C++) $(OPT) -o $(EXE) $(OBJS) $(LIBPATH) $(LIBS)
	@echo

all: clean install

clean:
	rm -f ../bin/* *.o *core
