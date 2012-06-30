CFLAGS=-march=native -ggdb3 -O2 -I./include -std=c99 -fopenmp 
LDFLAGS=-lm -lglfw

all: nbody

fileio.o: fileio.c
	gcc -c fileio.c -o fileio.o $(CFLAGS)
tree.o: tree.c
	gcc -c tree.c -o tree.o $(CFLAGS)
draw.o: draw.c
	gcc -c draw.c -o draw.o $(CFLAGS)
nbody.o: nbody.c
	gcc -c nbody.c -o nbody.o $(CFLAGS)
nbody: nbody.o tree.o draw.o fileio.o
	gcc nbody.o tree.o draw.o fileio.o -o nbody $(CFLAGS) $(LDFLAGS)  

clean:
	rm -rf *.o nbody