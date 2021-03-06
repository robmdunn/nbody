/*
Gravitational N-body Simulation
Copyright 2012 Robert Dunn

This file is part of Nbody.

Nbody is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

Nbody is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with Nbody.  If not, see <http://www.gnu.org/licenses/>.
*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "nbody.h"
#include "tree.h"
#include "draw.h"
#include "fileio.h"
#include "util.h"

#ifdef MPI
#include <mpi.h>
int  numtasks, rank, rc; 

int broadcastaccel(struct body * bodies, int srcrank)
{
	MPI_Bcast(&(bodies->ax),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->ay),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
}

int broadcastbody(struct body * bodies, int srcrank)
{
	MPI_Bcast(&(bodies->m),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->x),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->y),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->vx),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
	MPI_Bcast(&(bodies->vy),1,MPI_DOUBLE,srcrank,MPI_COMM_WORLD);
}

void broadcastbodies(int * nbodies, struct body ** bodies)
{
	MPI_Bcast(nbodies,1,MPI_INT,0,MPI_COMM_WORLD);
	if(rank!=0)
	{
		if(!(*bodies = malloc( (*nbodies)*sizeof(struct body))))
		{
			printf("Error: rank %d failed to allocate memory for %d bodies, Abort\n", rank, *nbodies);
			MPI_Abort(MPI_COMM_WORLD,0);
		}
	}
	for(int i = 0; i < *nbodies; i++)
	{
		broadcastbody((*bodies)+i,0);
	}
	return;
}

#endif

struct body * randinitbodies(const int nbodies, const double mass, const double spin, const double mzero)
//randomly initialize bodies in r and theta, and give them a spin around mzero 
{
	struct body * bodies;
	double r, theta;
	
	if(!(bodies = malloc(nbodies*sizeof(struct body))))  //allocate mem
	{
		printf("Error: failed to allocate memory for %d bodies, exit\n",nbodies);
		return NULL;
	}

	for(int i = 0; i < nbodies; i++)  //initialize with random position
	{
		bodies[i].m = mass;
		bodies[i].vx = 0.0;
		bodies[i].vy = 0.0;
		bodies[i].ax = 0.0;
		bodies[i].ay = 0.0;
		
		r = randf();
		theta = 2*PI*randf();
		
		bodies[i].x = r*cos(theta);
		bodies[i].y = r*sin(theta);
		
		if(spin!=0.0) //apply spin if specified
		{
			bodies[i].vx = sin(theta)*(spin+spin*0.1*randf())/(1+r);
			bodies[i].vy = -cos(theta)*(spin+spin*0.1*randf())/(1+r);
		}
		
		//printf("body %d: x=%f y=%f\n",i,bodies[i].x,bodies[i].y);
	}
	
	bodies[0].m = mzero;  //put bodies[0] at the center and give it a different mass
	bodies[0].x = 0.0;
	bodies[0].y = 0.0;
	bodies[0].vx = 0.0;
	bodies[0].vy = 0.0;
	
	return bodies;
	
}

void freebodies(struct body * bodies)
{
	free(bodies);
	return;
}

int runtimestep(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio)
//recalculate forces on bodies, then integrate force and velocity 
{
	struct node * rootnode;
	double xmin, xmax, ymin, ymax;
	
	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;
		
	for(int i = 0; i < nbodies; i++)  //reset accel
	{
		bodies[i].ax = 0.0;
		bodies[i].ay = 0.0;	
		xmin=min(xmin,bodies[i].x);
		xmax=max(xmax,bodies[i].x);
		ymin=min(ymin,bodies[i].y);
		ymax=max(ymax,bodies[i].y);
	}
	
	rootnode = createnode(bodies+0,xmin,xmax,ymin,ymax);
	
	for(int i = 1; i < nbodies; i++)
	{
		insertbody(bodies+i, rootnode);
	}
	
	#pragma omp parallel 
	{	
		#pragma omp for
		#ifdef MPI
		for(int i = rank; i < nbodies; i+=numtasks)  //sum accel
		#else
		for(int i = 0; i < nbodies; i++)  //sum accel
		#endif
		{
			treesum(rootnode, bodies+i, G, fudge, treeratio);
		}
	
		#ifdef MPI
		for(int i = 0; i < nbodies; i++)
		{
			broadcastaccel(bodies+i,i%numtasks);
		}
		#endif
	
		#pragma omp for
		for(int i = 0; i < nbodies; i++)
		{			
			bodies[i].vx += bodies[i].ax * timestep; // integrate accel into velcoity
			bodies[i].vy += bodies[i].ay * timestep;
			
			bodies[i].x += timestep*bodies[i].vx; //integrate velocity into position
			bodies[i].y += timestep*bodies[i].vy;
		} 		
	}
	
	#ifdef MPI
	if(rank==0)
	#endif
	{
		draw(bodies, nbodies, rootnode);
	}
	
	destroytree(rootnode);
	return 0;
}

int simulateloop(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio, const int write_interval, const char * outfile, const int graphics)
//simulation loop,  increments the timer and runs timesteps while the window remains open, writes bodies to file at specified intervals.
{
	double simtime = 0.0;
	int stepnum = 0;
	int run = 1;
	
	while (run)
	{

		simtime += timestep;
			
		runtimestep(bodies, nbodies, timestep, G, fudge, treeratio);
		stepnum++;
		if(write_interval!=0 && outfile && stepnum % write_interval == 0)
		{
			writebodies(outfile, bodies, nbodies, timestep, G, fudge, treeratio);
		}
		
		#ifdef MPI
		if(rank==0)
		#endif
		{
			printf("time: %f\r",simtime);
			if( windowopen() == 0 && graphics == 1) 
			{
				run = 0;
			}
			/*if(stepnum > 500)
			{
				run = 0;
			}*/
		}
		#ifdef MPI
		MPI_Bcast(&run,1,MPI_INT,0,MPI_COMM_WORLD); 
		#endif
		
	} 
	
	return 1;
	
}

void printhelp()
{
	printf("This program solves 2d F=G*m1*m2/(r^2+s) for n bodies\nAlgorthm based on Barnes-Hut tree\n");
	printf("Usage: nbody [options]\n");
	printf("List of options:\n");
	printf("\t-h\t\tshow this help output\n");
	printf("\t-dt <num>\tspecify timestep\n");
	printf("\t-tr <num>\tspecify tree ratio\n\t\t\t(distance to center of mass/quadrant diagonal size)\n");
	printf("\t-g <num>\tspecify G in force equation [6.67e-11]\n");
	printf("\t-sf <num>\tspecify softening factor s in force equation [0.005]\n");
	printf("\t-x <x> \tspecify window x size [600]\n");
	printf("\t-y <y> \tspecify window y size [600]\n");
	printf("\t-p <num>\tspecify body point size [1.0]\n");
	printf("\t-r <filename>\tresume from specified filename\n\t\t\t(if this is not set a random distribution is used)\n");
	printf("\t-o <filename>\toutput to file <filename>\n");
	printf("\t-nsteps <num>\tspecify output frequency in steps [100]\n");
	printf("\t-n <num>\tspecify number of bodies for random distribution [100]\n");
	printf("\t-m <num>\tspecify mass for random distribution [2000]\n");
	printf("\t-s <num>\tspecify spin for random distribution [0.05]\n");
	printf("\t-mz <num>\tspecify mass of body starting at 0,0 [1e7]\n");
	printf("\t-ng \tno graphics\n");

	return;
}

int main(int argc, char * argv[])
//main collects command line arguments, allocs/initializes bodies, then enters main simulation loop
{
	int nbodies = 100;
	struct body * bodies;
	double mass = 2000; 
	double G = 6.67384e-11;
	double timestep = 0.1;
	double fudge = 0.005;  //"Softening Factor" to prevent singularity in force calculation
	double spin = 0.05;
	double mzero = 10000000;
	char * infile = NULL;
	char * outfile = NULL;
	int write_interval = 100;
	double treeratio = 3; // threshold for ratio of distance to quadrant size for treesum
	int help = 0;
	int x = 600, y = 600;
	float pointsize = 1;
	int graphics = 1;
	
	#ifdef MPI	
	rc = MPI_Init(&argc,&argv);
	if (rc != MPI_SUCCESS) 
	{
		printf("Error in MPI_Init, Abort.");
		MPI_Abort(MPI_COMM_WORLD, rc);
	}
	MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	printf("nbody: MPI init on rank %d of %d\n",rank, numtasks); 
	#endif
	
	
	for(int i = 1; i < argc; i++)
	{
		if(!strcmp("-h",argv[i])) { help = 1; }
		if(!strcmp("-n",argv[i])) { nbodies = atoi(argv[i+1]); i++; }
		if(!strcmp("-dt",argv[i])) { timestep = atof(argv[i+1]); i++; }
		if(!strcmp("-tr",argv[i])) { treeratio = atof(argv[i+1]); i++; }
		if(!strcmp("-g",argv[i])) { G = atof(argv[i+1]); i++; }
		if(!strcmp("-sf",argv[i])) { fudge = atof(argv[i+1]); i++; }
		if(!strcmp("-m",argv[i])) { mass = atof(argv[i+1]); i++; }
		if(!strcmp("-s",argv[i])) { spin = atof(argv[i+1]); i++; }
		if(!strcmp("-mz",argv[i])) { mzero = atof(argv[i+1]); i++; }
		if(!strcmp("-r",argv[i])) { infile = argv[i+1]; i++; }
		if(!strcmp("-o",argv[i])) { outfile = argv[i+1]; i++; }
		if(!strcmp("-nsteps",argv[i])) { write_interval = atoi(argv[i+1]); i++; }
		if(!strcmp("-x", argv[i])) { x = atoi(argv[i+1]); i+=1; }
		if(!strcmp("-y", argv[i])) { y = atoi(argv[i+1]); i+=1; }
		if(!strcmp("-p", argv[i])) { pointsize = (float) atof(argv[i+1]); i+=1; }
		if(!strcmp("-ng", argv[i])) { graphics = 0; }
	}
	
	srand(1);  //lets seed RNG with a constant to make this easier to debug
	
	#ifdef MPI
	if(rank==0) 
	#endif
	{
		if(help)
		{
			printhelp();
			return 0;
		}
		
		if(infile) //if infile was specified...
		{
			if(!readbodies(infile, &bodies, &nbodies, &timestep, &G, &fudge, &treeratio)) //read file
			{
				printf("failed to read file %s\n",infile);
				return 1;
			}
		} else if(!(bodies = randinitbodies(nbodies, mass, spin, mzero)))  //else allocate + initialize bodies
		{
			printf("Failed to allocate memory for %d bodies\n",nbodies);
			return 1;
		}
		
		if(graphics != 0)
		{
			if(!initwindow(x, y, pointsize)) //initialize window
			{
				printf("Failed to initialize GL, exit\n");
				return 1;
			}
		}
		
		printf("nbodies = %d\n", nbodies);
		printf("timestep = %e\n", timestep);
		printf("tree threshold ratio = %f\n", treeratio);
		printf("G = %e\n",G);
		printf("softening factor = %f\n", fudge);
	}
	
	#ifdef MPI
	broadcastbodies(&nbodies,&bodies);
	#endif
	
	simulateloop(bodies, nbodies, timestep, G, fudge, treeratio, write_interval, outfile, graphics);
		
	closewindow();
	
	freebodies(bodies);
	
	#ifdef MPI
	MPI_Finalize();
	#endif
	return 0;
}
