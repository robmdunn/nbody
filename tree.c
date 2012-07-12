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
#include <stdlib.h>
#include <stdio.h>
#include "tree.h"
#include "nbody.h"
#include "util.h"

enum quadrant { q1, q2, q3, q4 };

enum quadrant getquadrant(double x, double y, double xmin, double xmax, double ymin, double ymax)
//makes a rectangle with bounds of xmin,xmax,ymin,ymax, and returns the quadrant that (x,y) is in
{
	double midx, midy;
	
	midx = xmin + 0.5*(xmax - xmin);
	midy = ymin + 0.5*(ymax - ymin);
	
	if(y > midy)
	{
		if(x > midx)
		{
			return q1;
		} else
		{
			return q2;
		}
	} else {
		if(x > midx)
		{
			return q4;
		} else
		{
			return q3;
		}		
	}
	
}



struct node * createnode(struct body * bodyp, double xmin, double xmax, double ymin, double ymax)
//creates a leaf node to insert into the tree
{
	struct node * rootnode;
	if(!(rootnode=malloc(sizeof(struct node))))
	{
		printf("Unable to allocate node, exit");
	}
	
	rootnode->totalmass = bodyp->m;
	rootnode->centerx = bodyp->x;
	rootnode->centery = bodyp->y;
	rootnode->xmin = xmin;
	rootnode->xmax = xmax;
	rootnode->ymin = ymin;
	rootnode->ymax = ymax;
		
	rootnode->bodyp = bodyp;
	rootnode->q1 = NULL;
	rootnode->q2 = NULL;
	rootnode->q3 = NULL;
	rootnode->q4 = NULL;
	
	return rootnode;
}

void updatecenterofmass(struct node * nodep, struct body * bodyp)
//updates the center of mass after inserting a point into a branch
{
	nodep->centerx = (nodep->totalmass*nodep->centerx + bodyp->m*bodyp->x)/(nodep->totalmass + bodyp->m);
	nodep->centery = (nodep->totalmass*nodep->centery + bodyp->m*bodyp->y)/(nodep->totalmass + bodyp->m);
	nodep->totalmass += bodyp->m;
	return;
}

void insertbody(struct body * insbody, struct node * nodep)
//inserts a body into the tree, converting leaf nodes into branches if necessary
{
	enum quadrant existingquad, newquad;
	double xmid, ymid;
	
	xmid = nodep->xmin + 0.5*(nodep->xmax - nodep->xmin);
	ymid = nodep->ymin + 0.5*(nodep->ymax - nodep->ymin);
		
	if(nodep->bodyp != NULL) //if the node is a leaf convert to a branch by inserting the leaf point into one of its subquadrants
	{
		existingquad = getquadrant(nodep->bodyp->x, nodep->bodyp->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
			
		switch (existingquad)
		{
			case q1:
				nodep->q1 = createnode(nodep->bodyp, xmid, nodep->xmax, ymid, nodep->ymax);
				break;
			case q2:
				nodep->q2 = createnode(nodep->bodyp, nodep->xmin, xmid, ymid, nodep->ymax);
				break;
			case q3:
				nodep->q3 = createnode(nodep->bodyp, nodep->xmin, xmid, nodep->ymin, ymid);
				break;
			case q4:
				nodep->q4 = createnode(nodep->bodyp, xmid, nodep->xmax, nodep->ymin, ymid);
				break;
		}
		nodep->bodyp = NULL;
	}
	
	newquad = getquadrant(insbody->x, insbody->y, nodep->xmin, nodep->xmax, nodep->ymin, nodep->ymax);
	
	updatecenterofmass(nodep,insbody);
	
	switch (newquad) //insert the new point into one of the quadrants if empty, otherwise recurse deeper into tree
	{
	case q1:
		if(nodep->q1 == NULL)
		{		
			nodep->q1 = createnode(insbody, xmid, nodep->xmax, ymid, nodep->ymax);		
		} else {
			insertbody(insbody,nodep->q1);
		}
		break;
	case q2:
		if(nodep->q2 == NULL)
		{			
			nodep->q2 = createnode(insbody, nodep->xmin, xmid, ymid, nodep->ymax);
		} else {			
			insertbody(insbody,nodep->q2);
		}
		break;
	case q3:
		if(nodep->q3 == NULL)
		{			
			nodep->q3 = createnode(insbody, nodep->xmin, xmid, nodep->ymin, ymid);
		} else {			
			insertbody(insbody,nodep->q3);
		}
		break;
	case q4:
		if(nodep->q4 == NULL)
		{			
			nodep->q4 = createnode(insbody, xmid, nodep->xmax, nodep->ymin, ymid);
		} else {			
			insertbody(insbody,nodep->q4);
		}
		break;
	}
		
}

void treesum(struct node * nodep, struct body * bodyp, double G, double fudge, double ratiothreshold )
//sum the forces on body bodyp from points in tree with root node nodep
{
	double dx, dy, r, rsqr; //x distance, y distance, distance, distance^2
	double force;
	double f_over_r;
	double quad_diag_sqr, quad_diag;
		
	dx = nodep->centerx - bodyp->x;
	dy = nodep->centery - bodyp->y;
	
	rsqr = pow(dx,2) + pow(dy,2);
	r = sqrt(rsqr);
	
	quad_diag_sqr = ( pow(nodep->xmax - nodep->xmin, 2) + pow(nodep->ymax - nodep->ymin, 2) );
	quad_diag = sqrt(quad_diag_sqr);
	
	if( (((r/quad_diag) > ratiothreshold) || (nodep->bodyp))&&(nodep->bodyp!=bodyp) )
	{
		force = (G * nodep->totalmass * bodyp->m)/(fudge + rsqr);

		f_over_r = force/r;
		
		bodyp->fx += f_over_r*dx;
		bodyp->fy += f_over_r*dy;		
		
	} else {
		if(nodep->q1) { treesum(nodep->q1, bodyp, G, fudge, ratiothreshold); }
		if(nodep->q2) { treesum(nodep->q2, bodyp, G, fudge, ratiothreshold); }
		if(nodep->q3) { treesum(nodep->q3, bodyp, G, fudge, ratiothreshold); }
		if(nodep->q4) { treesum(nodep->q4, bodyp, G, fudge, ratiothreshold); }
	}
	return;
}

void destroytree(struct node * nodep)
//recursively delete subnodes, then delete self
{
	if(nodep != NULL)
	{
		if(nodep->q1 != NULL) 
		{ 
			destroytree(nodep->q1); 
		}
		if(nodep->q2 != NULL) 
		{ 
			destroytree(nodep->q2); 
		}
		if(nodep->q3 != NULL) 
		{ 
			destroytree(nodep->q3); 
		}
		if(nodep->q4 != NULL) 
		{ 
			destroytree(nodep->q4); 
		}
		free(nodep);
	}
}
	
