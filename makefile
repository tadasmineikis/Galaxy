CC    = g++
CURENT = $(OPTIM)
CFLAGS= -c -Wall $(INCLUDE) $(CURENT) -fopenmp
LFLAGS= -Wall $(CURENT) 
DEBUG = -g3 -ggdb3
OPTIM = -O2  -march=native
OBJS  = main.o galaxy.o util.o neighbours.o
EXE = galaxy_2.0
RUN_DIR =
INCLUDE = 
LIB     = -larmadillo -lgsl -lgslcblas -lgomp

all: galemo

galemo: $(OBJS)
	$(CC) $(LFLAGS) $(OBJS) -o $(EXE) $(LIB)

main.o: galaxy.h main.cpp
	$(CC) $(CFLAGS) main.cpp

galaxy.o: galaxy.h galaxy.cpp util.h
	$(CC) $(CFLAGS) galaxy.cpp

neighbours.o: neighbours.h neighbours.cpp galaxy.h
	$(CC) $(CFLAGS) neighbours.cpp

util.o: util.h util.cpp
	$(CC) $(CFLAGS) util.cpp

clean:
	rm -f *.o *~ $(EXE)

tar:
	tar cfv $(EXE).tar *.h *.cpp makefile

zip:
	gzip -f $(EXE).tar

arch: tar zip
