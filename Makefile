CFLAGS=-march=native -ggdb3 -O2 -I./include -std=c99 -fopenmp 
LDFLAGS=-lm -lglfw

all: nbody

tree.o: tree.c
	gcc -c tree.c -o tree.o $(CFLAGS)
draw.o: draw.c
	gcc -c draw.c -o draw.o $(CFLAGS)
nbody.o: nbody.c
	gcc -c nbody.c -o nbody.o $(CFLAGS)
nbody: nbody.o tree.o draw.o 
	gcc nbody.o tree.o draw.o -o nbody $(CFLAGS) $(LDFLAGS)  

clean:
	rm -rf *.o nbody
