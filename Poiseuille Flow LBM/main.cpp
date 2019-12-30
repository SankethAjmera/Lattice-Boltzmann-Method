//
//  main.cpp
//  Assignment2_LBM
//
//  Created by Sanketh on 02/09/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#include <iostream>
#include <fstream> // for file stream
#include <sstream> // for string stream
#include <cstdlib> // standard library
#include <cmath>   // math library

// Define constants

const int nx = 4;           // Grid points along X
const int ny = 51;           // Grid points along Y [CHOOSE ODD NUMBER OF POINTS]
const int q = 9;             // Discretised velocity directions
const double ro = 36;         // Density guess
const double tau = 0.53;      // Relaxation time scale
const double u_max = 0.1;    // Maximum velocity
const double error_max = 1e-8; // Maximum error for convergence
const double ex[] = {0,1.0,0,-1.0,0,1.0,-1.0,-1.0,1.0}; // X velocities
const double ey[] = {0,0,1.0,0,-1.0,1.0,1.0,-1.0,-1.0}; // Y velocities
const double w[]  =  {16.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,4.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};  // Weights
const double L =  2.0; // Length of the domain
const double dy =  L/(ny-1.0); // Grid size
double err_con; // Convergence error
double err_abs; // Absolute error


// Dependent variables
const double mu   = ro*(tau-0.5)/3.0;           // Dynamic viscosity
const double dpdx = u_max*8.0*mu/pow(ny,2.0);   // Body force

double rho[nx][ny] ;        // Density array
double f[q][nx][ny];        // Probability density functions
double f_eq[q][nx][ny];     // Equilibrium probability density functions
double f_temp[q][nx][ny];   // Temporary probability density functions
double ux[nx][ny];          // Macroscopic X velocity
double uy[nx][ny];          // Macroscopic Y velocity
double u_abs[nx][ny];       // Macroscopic velocity
double u_abs_temp[nx][ny];  // Temporary macroscopic velocity
double y[ny];               // Y grid
double u_actual[ny];        // Analytical X velocity
double bc = 0;              // Switch for handling BCs : 1 for mid-grid; 0 for on-grid


void initialise();
void update();
void collide();
void stream();
void convergence_error();
void absolute_error();
void write();


int main() {

    initialise();
    int n = 0; // Number of iterations
    
    do {
        n += 1;
        update();
        collide();
        stream();
        convergence_error();
        absolute_error();
        if (n%100 ==0) {
            std::cout << "[# Iter "<< n <<"] Convergence error is " << err_con << " Absolute error is " << err_abs << " Max velocity is "<< ux[2][(ny-1)/2]<< "\n";
        }
        
    } while (err_con>error_max);
    
//    std::cout<< err_abs<<" "<< n ;
//    write();
    }

// Initialise variables
void initialise(){
    
    // Define initial values of paramteres
    err_con = 100.0;
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            // Macro variables
            rho[i][j] = ro;
            ux[i][j] = 0;
            uy[i][j] = 0;
            // Distribution functions
            for (int a =0 ; a<q; a++) {
                f[a][i][j]      = w[a]*rho[i][j];
                f_eq[a][i][j]   = w[a]*rho[i][j];
                f_temp[a][i][j] = 0;
            }
        }}
    
    // Define analytical solution for on-grid 
    for (int i=0; i<ny; i++) {
        if (i==0) {
            y[i] = -1;}
        else{ y[i] = y[i-1] + dy;
    }
        u_actual[i] = u_max*(1-pow(y[i],2.0));
}
    return;
}


// Update macroscopic variables
void update(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            rho[i][j] = 0;
            ux[i][j]  = 0;
            uy[i][j]  = 0;
            for (int a=0; a<q; a++) {
                rho[i][j] =  rho[i][j] + f[a][i][j];
                ux[i][j]  =  ux[i][j]  + ex[a]*f[a][i][j];
                uy[i][j]  =  uy[i][j]  + ey[a]*f[a][i][j];
            }
            u_abs_temp[i][j] = u_abs[i][j];
            ux[i][j] = (ux[i][j]+dpdx/2.0)/rho[i][j];
            uy[i][j] = (uy[i][j])/rho[i][j];
            u_abs[i][j] = pow(pow(ux[i][j],2.0)+pow(uy[i][j],2.0), 0.5);
        }
    }
    return;
}

// Perform collision
void collide(){
    
    for (int i=0; i<nx; i++) {
        if (bc == 1) { // Mid grid
            for (int j=0; j<ny; j++) {
                for (int a=0; a<q; a++) {
                    f_eq[a][i][j] = w[a]*rho[i][j]*(1.0+3.0*(ux[i][j]*ex[a]+uy[i][j]*ey[a])+4.5*pow(ux[i][j]*ex[a]+uy[i][j]*ey[a],2.0)-1.5*(pow(ux[i][j],2.0)+pow(uy[i][j],2.0)));
                    double omega = -(f[a][i][j]- f_eq[a][i][j])/tau;
                    double t1 = ex[a]-ux[i][j];
                    double t2 = ex[a]*ux[i][j]+ey[a]*uy[i][j];
                    double s  = (1.0-0.5/tau)*w[a]*(3.0*t1+9.0*t2*ex[a])*dpdx;
                    f_temp[a][i][j] = f[a][i][j]+omega+s;
                }
            }
        }else if (bc == 0){ // On grid
            for (int j=0; j<ny; j++) {
                for (int a=0; a<q; a++) {
                    if (j==0 || j == ny-1) {
                        f_temp[a][i][j] = f[a][i][j];
                    }else{
                    f_eq[a][i][j] = w[a]*rho[i][j]*(1.0+3.0*(ux[i][j]*ex[a]+uy[i][j]*ey[a])+4.5*pow(ux[i][j]*ex[a]+uy[i][j]*ey[a],2.0)-1.5*(pow(ux[i][j],2.0)+pow(uy[i][j],2.0)));
                    double omega = -(f[a][i][j]- f_eq[a][i][j])/tau;
                    double t1 = ex[a]-ux[i][j];
                    double t2 = ex[a]*ux[i][j]+ey[a]*uy[i][j];
                    double s  = (1.0-0.5/tau)*w[a]*(3.0*t1+9.0*t2*ex[a])*dpdx;
                    f_temp[a][i][j] = f[a][i][j]+omega+s;
                }
                }
            }
            
        }

        }
    return;
}


// Perform streaming
void stream(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int a=0; a<q; a++) {
                int ix = i+ex[a];
                int iy = j+ey[a];
                // Midway Bounceback boundary condition
                if (iy<0) {
                    f[2][i][0] = f_temp[4][i][0];
                    f[5][i][0] = f_temp[7][i][0];
                    f[6][i][0] = f_temp[8][i][0];
                }else if (iy>ny-1){
                    f[4][i][ny-1] = f_temp[2][i][ny-1];
                    f[7][i][ny-1] = f_temp[5][i][ny-1];
                    f[8][i][ny-1] = f_temp[6][i][ny-1];
                }else{
                    if (ix<0) {
                        ix = nx-1;
                    }if (ix>nx-1){
                        ix = 0;
                    }
                f[a][ix][iy] = f_temp[a][i][j];
            }
        }
    }
    }
    return;
}
        

// Calculate convergence error using L2 norm

void convergence_error(){
// Intialise temporary vaiables
    double temp1 = 0;
    double temp2 = 0;

    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            temp1 = temp1+pow(u_abs[i][j]-u_abs_temp[i][j], 2.0);
            temp2 = temp2+pow(u_abs[i][j], 2.0);
        }
    }
    err_con = pow(temp1/temp2, 0.5);
    return;
}

// Calculate absolute error using L2 norm
void absolute_error(){
    // Intialise temporary vaiables
    double temp1 = 0;
    double temp2 = 0;
    
//    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            temp1 = temp1+pow(ux[2][j]-u_actual[j], 2.0);
            temp2 = temp2+pow(u_actual[j], 2.0);
        }
//    }
    err_abs = pow(temp1/temp2, 0.5);
    return;
}

// Write files

void write(){
    std::stringstream file_name;
//    if (bc == 0) {
//        file_name  <<"Ongrid_"<<ny<<".dat";
//    }else if (bc ==1){
//        file_name  <<"Midgrid_"<<ny<<".dat";
//    }
    
    file_name << "abc.dat";
    
    std::ofstream file_out;
    file_out.open(file_name.str());
    
    for (int i = 0; i<ny; i++) {
        file_out << ux[2][i] << "\t"<< y[i];
    }
    file_out.close();
    return;
}


