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

#ifndef TREE
#define TREE

enum quadrant { q1, q2, q3, q4 };
enum method { method_allpoints, method_tree };

struct node {
	double totalmass;	
	double centerx, centery;
	double xmin, xmax;
	double ymin, ymax;
	struct body * bodyp;
	struct node * q1;
	struct node * q2;
	struct node * q3;
	struct node * q4;
};

double max (double a, double b);
double min (double a, double b);
enum quadrant getquadrant(double x, double y, double xmin, double xmax, double ymin, double ymax);
struct node * createnode(struct body * bodyp, double xmin, double xmax, double ymin, double ymax);
void updatecenterofmass(struct node * nodep, struct body * bodyp);
void insertbody(struct body * insbody, struct node * nodep);
void treesum(struct node * nodep, struct body * bodyp, double G, double fudge, double ratiothreshold );
void destroytree(struct node * nodep);
	
#endif
	
