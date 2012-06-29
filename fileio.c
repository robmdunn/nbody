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

#include <stdio.h>
#include "nbody.h"


int readbodies(const char * infilename, struct body * bodies, int *nbodies, double *timestep, double *G, double *fudge, double *treeratio)
{
	FILE * infile;
	if(!(infile = fopen(infilename, "r")));
	{
		return 0;
	}
	fscanf(infile, "%le %le %le %le %d",timestep, G, fudge, treeratio, nbodies);
	for(int i = 0; i < *nbodies; i++)
	{
		fscanf(infile,"%le %le %le %le %le", &bodies[i].m, &bodies[i].x, &bodies[i].y, &bodies[i].vx, &bodies[i].vy);
	}
	return 1;
}

int writebodies(const char * outfilename, const struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio)
{
	FILE * outfile; 
	if(!(outfile = fopen(outfilename, "w")))
	{
		return 0;
	}
	fprintf(outfile,"%1.16e\n%1.16e\n%1.16e\n%1.16e\n%d\n",timestep, G, fudge, treeratio, nbodies);
	for(int i = 0; i < nbodies; i++)
	{
		fprintf(outfile,"%1.16e %1.16e %1.16e %1.16e %1.16e\n", bodies[i].m, bodies[i].x, bodies[i].y, bodies[i].vx, bodies[i].vy);
	}
	return 1;
	
}
