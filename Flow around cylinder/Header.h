//
//  Header.h
//  Assignment 5
//
//  Created by Sanketh on 11/10/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#ifndef Header_h
#define Header_h


// Declare constants

const int nx = 400;            // Grid points along X
const int ny = 300;            // Grid points along Y
const int q = 9;               // Discretised velocity directions
const double ro = 1.0;         // Density initial guess
const double rho_out = 1.0;    // Density at outlet (Pressure boundary condition)
const double h = ny;           // Length of the domain in X and Y
const double dy = 1.0/(ny-1);  // Grid size in Y
const double dx = 1.0/(nx-1);  // Grid size in X
const int Re = 46;             // Reynolds number
const double radius = 10.0;    // Radius of the cylinder
const double u_in = 0.05;      // Velocity at inlet (Velocity boundary condition)
const double nu   = u_in*2.0*radius/Re; // Kinematic viscosity
const double tau  = 3.0*nu+0.5;         // Relaxation time scale
const int x_0    = nx/3;       // x coordinate of centre of cylinder
const int y_0    = ny/2;       // y coordinate of centre of cylinder
const int coll_dir[] = {0,3,4,1,2,7,8,5,6};// Collision direction
const int n_iter = 200000;      // Number of time steps for simulation
const int n_skip = 10;          // Skip time steps for convergence
double err_vel = 1;             // Convergence error in velocity
double err_f = 1;               // Convergence error in pdf



// D2Q9
const double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0}; // X velocities
const double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0}; // Y velocities
const double w[]  = {16.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};  // Weights
int n = 0;                      // Number of iterations


// Declare macroscopic variables

double rho[nx][ny] ;        // Density array
double f[q][nx][ny];        // Probability density functions
double f_old[q][nx][ny];    // Probability density functions for last time step
double f_eq[q][nx][ny];     // Equilibrium probability density functions
double f_temp[q][nx][ny];   // Probability density functions (at previous time step)
double ux[nx][ny];          // Macroscopic X velocity
double uy[nx][ny];          // Macroscopic Y velocity
double u_abs[nx][ny];       // Macroscopic absolute velocity
double u_abs_temp[nx][ny];  // Macroscopic absolute velocity (at previous time step)
double y[ny];               // Y grid
double x[nx];               // X grid
int    counter  = 0;        // Counter for time steps for force calculation
int    n_force_save = 10;   // Save force after these many time steps
int    arr_size  = n_iter/n_force_save; // Size of tiem dependent force array
double Fx[20000];        // Time varying force in X direction
double Fy[20000];        // Time varying force in Y direction
double Cd[20000];        // Drag coefficient
double Cl[20000];        // Lift coefficient
int int_nodes[nx][ny];      // Location of interior nodes of the cylinder
int bd_nodes[nx][ny];       // Location of boundary nodes of the cylinder

// Declare functions

void initialise();
void detect_nodes();
void update();
void collide();
void stream();
void force();
void write_force();
void write_nodes();
void write_grid();
void write_velocity();
void convergence_error();



#endif /* Header_h */

