#ifndef DEFINES_H
#define DEFINES_H

#include <fstream>
#include <string>
#include <iostream>
#include <cmath>
#include <vector>
#include <ctime>
#include <boost/tuple/tuple.hpp>
#define GNUPLOT_ENABLE_PTY
#include "gnuplot-iostream.h"
#include <omp.h>
#include "Obstacle.h" // Obstacle class definition

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%     MACROS           %%%%%%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


#define PI 3.14159265359
#define IX(i,j) ((i)+(N+2)*(j))
#define SWAP(x0,x) {float *tmp=x0;x0=x;x=tmp;}


//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%   SIMULATION PARAMETERS  %%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


const int N = 200; // simulation grid size
const int maxDraw = 50; // define maximal number of drawn velocity field vectors in each dimension
const int drawRatio = N/maxDraw; // scaling factor for drawn velocity field vectors
const int size = (N+2)*(N+2); // grid size incl. boundaries

const float dt = 0.0001; // incremental time step length
float t = 0.; // current simulation time
float tmax = 1; // maximal simulation time
int counter = 0; // time iteration counter

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%   PREALLOCATION          %%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Gnuplot densPlot; // Gnuplot Object for density plot
Gnuplot velPlot; // Gnuplot Object for velocity field plot
Gnuplot absVelPlot; // Gnuplot Object for maginutude of velocity field plot

float u[size]; // fluid field variables
float v[size];
float u_prev[size];
float v_prev[size];
float dens[size];
float dens_prev[size];

bool occupiedGrid[size]; // define flow obstacles
float Nx[size]; // coordinate variable x
float Ny[size]; // coordinate variable y

int nObstacles = 0; // number of obstacles in fluid



//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%     FLUID PROPERTIES     %%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


float visc = 0.001; // viscosity
float diff = 0.01; // diffusion rate

//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
//%%%%%%%%%   FUNCTION PROTOTYPES    %%%
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

// Initialization functions
void initializeGrid(); // initialize grid variables
Obstacle defineObstacles(); // define fluid obstacle
void initializeFluid(float *u,float *v,float *dens); // initialize fluid velocities and density distribution at t=0


// Fluid dynamics functions

void get_from_UI(float *dens, float *u, float *v, float *dens_prev, float *u_prev, float *v_prev , float t,Gnuplot &); // static density input
void vel_step(int N, float *u, float *v, float *u_prev, float *v_prev, float visc, float dt); // determine velocity vectors in next time step
void dens_step(int N, float *dens, float *dens_prev, float *u, float *v, float diff, float dt); // determine fluid field in next time step

void project(int N, float *u, float *v, float *u0, float *v0 );
void add_source(int N, float *x, float *x0, float dt );
void set_bnd(int N, int b, float *x);
void diffuse(int N, int,float *x,float *x0, float diff, float dt);
void advect( int N, int b, float * d, float * d0, float * u, float * v, float dt );

Obstacle obstacle_step(Obstacle &particle); // update obstacle position and velocity according to fluid
void updateGrid(Obstacle &particle); // update occupied cell-vector

// Visualization
void drawVel(); // Gnuplot of velocity and density
void drawDens(); //...
void drawAbsVel();
//void drawObstacle();

#endif
