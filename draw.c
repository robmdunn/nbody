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
#include <GL/glfw.h>
#include "nbody.h"
#include "tree.h"
#include "draw.h"

int initwindow(int x, int y, float pointsize)
{
	
	if(!glfwInit())	//initialize GLFW
	{
		printf("Failed to initialize GLFW\n");
		return 0;
	}
	
	if( !glfwOpenWindow( x, y, 0,0,0,0, 32,0, GLFW_WINDOW ) )
	{
		printf("Failed to open GLFW window\n");
		glfwTerminate();
		return 0;
	}

	glfwSetWindowTitle("Nbody");
	
	glClearColor(0.0f, 0.0f, 0.1f, 0.0f);
	glColor4f(1.0f, 1.0f, 1.0f,0.0f);
		
	glfwEnable( GLFW_STICKY_KEYS );
	
	glDisable(GL_DEPTH_TEST);

	glOrtho (-1.5f, 1.5f, -1.5f, 1.5f, 0, 1);
	
	glPointSize(pointsize);
	glClear(GL_COLOR_BUFFER_BIT);
	
	//printf("GL Init\n");
	
	return 1;
}

void closewindow()
{
	//printf("\nTerminate GL\n");
	glfwTerminate();
}

void drawbodies(struct body * bodies, const int nbodies)
{
	glColor4f(1.0f, 1.0f, 1.0f,0.0f);

	glBegin(GL_POINTS);
	for(int i=0; i < nbodies; i++)
	{
		glVertex2f(bodies[i].x,bodies[i].y);
	}
	glEnd();
}

void drawtreelines(struct node * nodep)
{
	glVertex2f(nodep->xmin,nodep->ymin);
	glVertex2f(nodep->xmax,nodep->ymin);
	
	glVertex2f(nodep->xmax,nodep->ymin);
	glVertex2f(nodep->xmax,nodep->ymax);
	
	glVertex2f(nodep->xmax,nodep->ymax);
	glVertex2f(nodep->xmin,nodep->ymax);
	
	glVertex2f(nodep->xmin,nodep->ymax);
	glVertex2f(nodep->xmin,nodep->ymin);

	
	if(nodep->q1){ drawtreelines(nodep->q1); }
	if(nodep->q2){ drawtreelines(nodep->q2); }
	if(nodep->q3){ drawtreelines(nodep->q3); }
	if(nodep->q4){ drawtreelines(nodep->q4); }	
}

void drawtree(struct node * nodep)
{
	glColor4f(0.3f, 0.3f, 0.3f, 0.0f);
	
	glBegin(GL_LINES);

	drawtreelines(nodep);

	glEnd();
}

void draw(struct body * bodies, const int nbodies, struct node * nodep)
{
	glClear(GL_COLOR_BUFFER_BIT);

	drawtree(nodep);	
	drawbodies(bodies, nbodies);
	
	glfwSwapBuffers();
}

int windowopen()
{
	return glfwGetKey( GLFW_KEY_ESC ) != GLFW_PRESS && glfwGetWindowParam( GLFW_OPENED ) ;
}
