//
//  Header.h
//  Assignment 6
//
//  Created by Sanketh on 06/11/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#ifndef Header_h
#define Header_h


// Declare constants

const int nx = 100;            // Grid points along X
const int ny = 100;            // Grid points along Y
const int q = 9;               // Discretised velocity directions
const double rho0 = 1.0;       // Liquid density
double rho_c = rho0*log(2.0);  // critical density
double tau = 1.0;              // relaxation time
double G_c = -2.0/(9.0*rho0);  // critical G
double r = 1.375;              // Choose ratio of G and G_c
double g   = -r*4.0/(9.0*rho0);// Get g
const double h = ny;           // Length of the domain in X and Y
const double dy = 1.0/(ny-1);  // Grid size in Y
const double dx = 1.0/(nx-1);  // Grid size in X
const int n_iter = 20000;      // Number of time steps for simulation
const int n_skip = 10;         // Skip time steps for convergence
const int n_save = 100;        // Save time snaps
double err_vel = 1.0;          // Convergence error in velocity
double err_rho = 1.0;          // Convergence error in pdf
double psi[nx][ny];            // Value of psi at a given node

// D2Q9
const double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0}; // X velocities
const double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0}; // Y velocities
const double w[]  = {16.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};  // Weights
int n = 0;                      // Number of iterations


// Declare macroscopic variables

double rho[nx][ny] ;        // Density array
double rho_temp[nx][ny] ;   // Density array from last time step
double f[q][nx][ny];        // Probability density functions
double f_eq[q][nx][ny];     // Equilibrium probability density functions
double f_temp[q][nx][ny];   // Probability density functions (at previous time step)
double ux[nx][ny];          // Macroscopic X velocity
double uy[nx][ny];          // Macroscopic Y velocity
double u_abs[nx][ny];       // Macroscopic absolute velocity
double u_abs_temp[nx][ny];  // Macroscopic absolute velocity (at previous time step)
double y[ny];               // Y grid
double x[nx];               // X grid
double Fx[nx][ny];          // Time varying force in X direction
double Fy[nx][ny];          // Time varying force in Y direction

// Declare functions

void initialise();
void update();
void collide();
void stream();
void force();
void convergence_error();
void write_data();



#endif /* Header_h */

