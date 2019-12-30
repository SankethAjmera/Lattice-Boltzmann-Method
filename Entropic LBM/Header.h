//
//  Header.h
//  Entropic LBM
//
//  Created by Sanketh on 18/11/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#ifndef Header_h
#define Header_h


// Declare constants

const int nx = 100;            // Grid points along X
const int ny = 100;            // Grid points along Y
const int q = 9;              // Discretised velocity directions
const double ro = 1.0;        // Density initial guess
const double Re = 1000.0;     // Reynolds number
const double u_wall = 0.01;  // Wall velocity
const double nu   = u_wall*ny/Re;    // Kinematic viscosity
const double beta = 1.0/(6.0*nu+1.0);// Beta (0,1)
const double h = ny;          // Length of the domain in X and Y
const double dy = 1.0/ny;     // Grid size in X and Y

const double error_max = 1e-5; // Maximum error for convergence
// D2Q9
const double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0}; // X velocities
const double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0}; // Y velocities
const double w[]  =  {16.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};  // Weights
double err_con   = 1.0;      // Convergence error
double err_alpha = 1e-6;     // Convergence error for alpha
double alpha_old[nx][ny] ;   // Alpha at each location
double H;                    // H-function at a given location
double H_der;                // H-function derivative at a given location
double H_alpha;              // H-function at a given location, as a function of alpha
double H_der_alpha;          // H-function derivative at a given location, as a function of alpha
int n_iter = 250000;         // Number of iterations
int n = 0;                   // Iterations counter

// Declare macroscopic variables

double rho[nx][ny] ;        // Density array
double f[q][nx][ny];        // Probability density functions
double f_old[q][nx][ny];    // Probability density functions for last time step
double f_eq[q][nx][ny];     // Equilibrium probability density functions
double f_temp[q][nx][ny];   // Temporary probability density functions
double ux[nx][ny];          // Macroscopic X velocity
double uy[nx][ny];          // Macroscopic Y velocity
double u_abs[nx][ny];       // Macroscopic absolute velocity
double u_abs_temp[nx][ny];  // Temporary variable for macroscopic absolute velocity
double y[ny];               // Y grid
double x[nx];               // X grid

// Declare functions

void initialise();
void update();
void collide();
void stream();
void convergence_error();
void absolute_error();
void calculate_alpha();
void entropy_function(int x, int y);
void write_grid();
void write_velocity();



#endif /* Header_h */

