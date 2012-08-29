CFLAGS=-march=native -ggdb3 -O2 -I./include -std=c99 -fopenmp
LDFLAGS=-lm -lglfw
MPIFLAGS= -DMPI

all: nbody

mpi: nbody_mpi

util.o: util.c
	gcc -c util.c -o util.o $(CFLAGS)
fileio.o: fileio.c
	gcc -c fileio.c -o fileio.o $(CFLAGS)
tree.o: tree.c
	gcc -c tree.c -o tree.o $(CFLAGS)
draw.o: draw.c
	gcc -c draw.c -o draw.o $(CFLAGS)
nbody.o: nbody.c
	gcc -c nbody.c -o nbody.o $(CFLAGS)
nbody_mpi.o: nbody.c
	mpicc -c nbody.c -o nbody_mpi.o $(CFLAGS) $(MPIFLAGS)
nbody: nbody.o tree.o draw.o fileio.o util.o
	gcc nbody.o tree.o draw.o fileio.o util.o -o nbody $(CFLAGS) $(LDFLAGS)
nbody_mpi: nbody_mpi.o tree.o draw.o fileio.o util.o
	mpicc nbody_mpi.o tree.o draw.o fileio.o util.o -o nbody $(CFLAGS) $(LDFLAGS) $(MPIFLAGS)

clean:
	rm -rf *.o nbody


	
