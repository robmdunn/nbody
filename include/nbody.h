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

#ifndef PI
#define PI 3.141592653589793238462
#endif

#ifndef NBODY
#define NBODY

struct body {
	double m; //mass
	double x, y; //x,y position
	double vx, vy; //x,y velocity
	double fx, fy; //x,y force
};

double randf();
struct body * randinitbodies(const int nbodies, const double mass, const double spin, const double mzero);
int runtimestep(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio);
int simulateloop(struct body * bodies, const int nbodies, const double timestep, const double G, const double fudge, const double treeratio);
void freebodies(struct body * bodies);


#endif
