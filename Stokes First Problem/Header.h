//
//  Header.h
//  Assignment 4
//
//  Created by Sanketh on 20/09/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#ifndef Header_h
#define Header_h


// Declare constants

const int nx = 10;              // Grid points along X
const int ny = 200;            // Grid points along Y [CHOOSE ODD NUMBER OF POINTS]
const int q = 9;               // Discretised velocity directions
const double ro = 1.0;        // Density initial guess
const double tau = 1;          // Relaxation time scale
const double h = ny;           // Length of the domain in X and Y
const double dy = 1.0/ny;      // Grid size in X and Y
const double nu   = (tau-0.5)/3.0;  // Dynamic viscosity
const double u_wall = 0.01;         // Wall velocity
const double len_scale = nu/u_wall; // Diffusion length scale
const double time_scale = pow(h, 2.0)/nu; // Diffusion time scale
const int n_iter = 5000;                  // Number of time steps


// D2Q9
const double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0}; // X velocities
const double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0}; // Y velocities
const double w[]  =  {16.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};  // Weights
double err_con = 100.0; // Convergence error
int n = 0;              // Number of iterations


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
double x[ny];               // X grid
int    counter  = 0;        // Counter for time steps for force calculation
double sf[n_iter];   // Shear force transient

// Declare functions

void initialise();
void update();
void collide();
void stream();
void force();
void write_force();
void write_velocity(int n_iter);
void write_grid(int nx,int ny);


#endif /* Header_h */
