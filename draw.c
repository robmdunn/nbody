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

int initwindow(int x, int y)
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
	
	glPointSize(1);
	glClear(GL_COLOR_BUFFER_BIT);
	
	printf("GL Init\n");
	
	return 1;
}

void closewindow()
{
	printf("\nTerminate GL\n");
	glfwTerminate();
}

void drawbodies(struct body * bodies, const int nbodies)
{
	glColor4f(1.0f, 1.0f, 1.0f,0.0f);
	for(int i=0; i < nbodies; i++)
	{
		glBegin(GL_POINTS);
		glVertex2f(bodies[i].x,bodies[i].y);
		glEnd();
	}
}

void drawtree(struct node * nodep)
{
	glColor4f(0.3f, 0.3f, 0.3f, 0.0f);
	
	glBegin(GL_LINES);
	glVertex2f(nodep->xmin,nodep->ymin);
	glVertex2f(nodep->xmax,nodep->ymin);
	glEnd();
	
	glBegin(GL_LINES);
	glVertex2f(nodep->xmax,nodep->ymin);
	glVertex2f(nodep->xmax,nodep->ymax);
	glEnd();
	
	glBegin(GL_LINES);
	glVertex2f(nodep->xmax,nodep->ymax);
	glVertex2f(nodep->xmin,nodep->ymax);
	glEnd();
	
	glBegin(GL_LINES);
	glVertex2f(nodep->xmin,nodep->ymax);
	glVertex2f(nodep->xmin,nodep->ymin);
	glEnd();
	
	if(nodep->q1){ drawtree(nodep->q1); }
	if(nodep->q2){ drawtree(nodep->q2); }
	if(nodep->q3){ drawtree(nodep->q3); }
	if(nodep->q4){ drawtree(nodep->q4); }	
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
