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
#include <omp.h>
#include "nbody.h"
#include "tree.h"
#include "draw.h"
#include "fileio.h"

double randf() 
{
	return (double)rand() / (double)RAND_MAX;
}

struct body * randinitbodies(const int nbodies, const double mass, const double spin, const double mzero)
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
		bodies[i].fx = 0.0;
		bodies[i].fy = 0.0;
		
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
	printf("Free mem\n");
	free(bodies);
	return;
}

int runtimestep(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio)
{
	double dx, dy, r, rsqr; //x distance, y distance, distance, distance^2
	double force;
	double G_m1;
	double f_over_r;
	double dx_f_over_r, dy_f_over_r;
	struct node * rootnode;
	double xmin, xmax, ymin, ymax;
	
	xmin = 0.0;
	xmax = 0.0;
	ymin = 0.0;
	ymax = 0.0;
		
	for(int i = 0; i < nbodies; i++)  //reset forces
	{
		bodies[i].fx = 0.0;
		bodies[i].fy = 0.0;	
		xmin=min(xmin,bodies[i].x);
		xmax=max(xmax,bodies[i].x);
		ymin=min(ymin,bodies[i].y);
		ymax=max(ymax,bodies[i].y);
	}
	
	printf("xmin=%f xmax=%f ymin=%f ymax=%f ",xmin,xmax,ymin,ymax);
	
	rootnode = createnode(bodies+0,xmin,xmax,ymin,ymax);
	
	for(int i = 1; i < nbodies; i++)
	{
		insertbody(bodies+i, rootnode);
	}
	
	#pragma omp parallel private(dx,dy,r,rsqr,force,G_m1,f_over_r,dx_f_over_r,dy_f_over_r)
	{		
		#pragma omp for schedule(static,1)
		for(int i = 0; i < nbodies; i++)  //sum forces
		{			
			treesum(rootnode, bodies+i, G, fudge, treeratio);
			//G_m1 = G*bodies[i].m;
			
			//for(int j = i+1; j < nbodies; j ++)
			//{ 

				//dx = bodies[j].x - bodies[i].x;
				//dy = bodies[j].y - bodies[i].y;
				//rsqr = pow(dx,2) + pow(dy,2);
			
				//force = (G_m1*bodies[j].m)/(fudge+rsqr);
								
				//r = sqrt(rsqr);
				//f_over_r = force/r;
				
				//dx_f_over_r = f_over_r*dx;
				//dy_f_over_r = f_over_r*dy;
	
				//bodies[i].fx += dx_f_over_r;
				//bodies[i].fy += dy_f_over_r;
				//bodies[j].fx -= dx_f_over_r;
				//bodies[j].fy -= dy_f_over_r;
				
				////printf("p1=%d p2=%d dist=%e f=%e\n",i,j,dist,force);
			//}
		}

		#pragma omp for
		for(int i = 0; i < nbodies; i++) 
		{

			double timestep_over_m = timestep / bodies[i].m;
			bodies[i].vx += bodies[i].fx * timestep_over_m; // integrate force into velcoity
			bodies[i].vy += bodies[i].fy * timestep_over_m;
			
			bodies[i].x += timestep*bodies[i].vx; //integrate velocity into position
			bodies[i].y += timestep*bodies[i].vy;
			//printf("p=%d vx=%e vy=%e x=%e y=%e\n", i, bodies[i].vx, bodies[i].vy, bodies[i].x, bodies[i].y);
		} 
		
	}	
	draw(bodies, nbodies, rootnode);
	
	destroytree(rootnode);
	return 0;
}

int simulateloop(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio)
{
	double simtime = 0.0;
	
	while (windowopen())
	{
		
		simtime += timestep;
		printf("time: %f\r",simtime);
				
		runtimestep(bodies, nbodies, timestep, G, fudge, treeratio);
					
	} 
	
	return 1;
	
}

int main(int argc, char * argv[])
{
	int nbodies;
	struct body * bodies;
	double mass = 2000; 
	double G = 6.67384e-11;
	double timestep = 0.1;
	double fudge = 0.005;  //"Softening Factor" to prevent singularity in force calculation
	double spin = 0.05;
	double mzero = 10000000;

	double treeratio = 3; // threshold for ratio of distance to quadrant size for treesum
	
	if(argc!=2 && argc!= 5 && argc!=7 && argc!=10)
	{
		printf("usage: %s <nbodies> [<mass> <timestep> <tree ratio>[<G> <fudge> [<spin> <mzero>]]]\n",argv[0]);
		return 1;
	}
	
	nbodies = atoi(argv[1]);  // get a pile of command line arguments
	if(argc==5 || argc == 7 || argc==9)
	{
		mass = atof(argv[2]);
		timestep = atof(argv[3]);
		treeratio = atof(argv[4]);
		if(argc==7 || argc==9)
		{
			G = atof(argv[5]);
			fudge = atof(argv[6]);
			if(argc==9)
			{
				spin = atof(argv[7]);
				mzero = atof(argv[8]);
			}
		}
	}
		
	printf("Using nbodies = %d\n", nbodies);
	printf("mass = %e kg\n", mass);
	printf("timestep = %e sec\n", timestep);
	printf("tree threshold ratio = %f\n", treeratio);
	printf("G = %e N (m/kg)^2\n",G);
	printf("fudge = %f\n", fudge);
	printf("spin = %f\n", spin);
	printf("mzero mass = %e\n",mzero);
	
	srand(1);  //lets seed RNG with a constant to make this easier to debug
	
	
	if(!(bodies = randinitbodies(nbodies, mass, spin, mzero)))  //allocate + initialize bodies
	{
		printf("Failed to allocate memory for %d bodies\n",nbodies);
		return 1;
	}
	
	if(!initwindow(600, 600)) //initialize window
	{
		printf("Failed to initialize GL, exit\n");
		return 1;
	}
	writebodies("test.bdy", bodies, nbodies, timestep, G, fudge, treeratio);
	printf("Entering simulation loop\n");
	simulateloop(bodies, nbodies, timestep, G, fudge, treeratio);
		
	closewindow();
	
	freebodies(bodies);
	
	return 0;
}
