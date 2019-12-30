//
//  main.cpp
//  Assignment 6 - Phase change via Shen-Chang model
//
//  Created by Sanketh on 06/11/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//

#include <iostream>
#include <fstream> // for writing files
#include <sstream> // for string handling
#include <cstdlib> // standard library
#include <cmath>   // math library
#include <iomanip> // for setting precision
#include <string>
#include "Header.h"

using namespace std;

int main() {
    initialise();
    cout << "Starting simulation for G_ratio = "<< r <<endl;
    do {
        n += 1;
        force();
        update();
        collide();
        stream();
        convergence_error();
        write_data();
    } while (n<n_iter);
    cout << "Simulation for G_ratio = "<< r << " completed at " << n_iter << " iterations"<<endl;
    return 0;
}


// Initialise variables
void initialise(){
    
    
    // On grid
    y[0] = 0;  // Y
    for (int i=1; i<ny; i++) {
        y[i] = y[i-1]+dy;
    }
    x[0] = 0;  // X
    for (int i=1; i<nx; i++) {
        x[i] = x[i-1]+dx;
    }
    
    // Define initial values of paramteres
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            // Temporary variables for the last time step data
            double R =(double)rand()/RAND_MAX*2.0-1.0; // Random number between -1 and 1
            u_abs_temp[i][j] = 0.0;
//            rho[i][j] = rho_c;
//            if ((i == nx/2 || i == nx/3 || i == 2*nx/3) &&( j == ny/2  || j == nx/3 || j == 2*nx/3 )) {
                rho[i][j] = rho_c*(1.0+pow(10.0,-3.0)*R);
//            }
            
            Fx[i][j] = 0.0;
            Fy[i][j] = 0.0;
            // Distribution functions
            for (int a =0 ; a<q; a++) {
                // Initial condition corresponds to zero velocity throughout the domain
                f[a][i][j]      = w[a]*rho[i][j];
                f_eq[a][i][j]   = w[a]*rho[i][j];
                f_temp[a][i][j] = f_eq[a][i][j];
            }
        }
    }
}



// Calcuate force

void force(){
    
    // Get rho from f
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            rho_temp[i][j]   = rho[i][j];
            rho[i][j] = 0.0;
            for (int a=0; a<q; a++) {
                rho[i][j] += f[a][i][j];
            }
            psi[i][j]  =  rho0*(1.0-exp(-rho[i][j]/rho0));
        }
    }
    
    
    
    
    double dpsi_x;
    double dpsi_y;
    double c1 = 1.0/3.0;
    double c2 = 1.0/12.0;
    
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
    
            
            dpsi_x = 0.0;
            dpsi_y = 0.0;
            
            
            // Neighbouring nodes
            int i_p = i+1;
            int i_n = i-1;
            int j_p = j+1;
            int j_n = j-1;
            
            // Periodic condition
            if (i_p > nx-1) {
                i_p = 0;
            }
            if (i_n < 0) {
                 i_n = nx-1;
             }
            if (j_p > ny-1) {
                 j_p = 0;
             }
            if (j_n < 0) {
                 j_n = ny-1;
             }
            
            // Calculate the gradient of psi
            
            dpsi_x = c1*(psi[i_p][j]-psi[i_n][j])+c2*(psi[i_p][j_p]-psi[i_n][j_p]+psi[i_p][j_n]-psi[i_n][j_n]);
            dpsi_y =
            c1*(psi[i][j_p]-psi[i][j_n])+c2*(psi[i_p][j_p]-psi[i_p][j_n]+psi[i_n][j_p]-psi[i_n][j_n]);

            
            Fx[i][j] = -3.0*g*psi[i][j]*dpsi_x;
            Fy[i][j] = -3.0*g*psi[i][j]*dpsi_y;
            
//            for (int a = 0; a<9; a++) {
//                if (a<5) {
//                    fx_temp += -3.0*g*ex[a]*psi(i,j)*dpsi_x[i][j];
//                    fy_temp += -3.0*g*ey[a]*psi(i,j)*dpsi_y[i][j];
//                }else{
//                    fx_temp += -g*ex[a]*psi(i,j)*dpsi_x[i][j]/4.0;
//                    fy_temp += -g*ey[a]*psi(i,j)*dpsi_y[i][j]/4.0;
//                }
//            }
            

        }
    }
}


// Update macroscopic variables

void update(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
//            rho_temp[i][j]   = rho[i][j];
//            rho[i][j] = 0;
            ux[i][j]  = 0;
            uy[i][j]  = 0;
            for (int a=0; a<q; a++) {
//                rho[i][j] += f[a][i][j];
                ux[i][j]  += ex[a]*f[a][i][j];
                uy[i][j]  += ey[a]*f[a][i][j];
            }
            // Update variables at new time step
            u_abs_temp[i][j] = u_abs[i][j];
            ux[i][j] = (ux[i][j]+Fx[i][j]/2.0)/rho[i][j];
            uy[i][j] = (uy[i][j]+Fy[i][j]/2.0)/rho[i][j];
            u_abs[i][j] = pow(pow(ux[i][j],2.0)+pow(uy[i][j],2.0), 0.5);
        }
    }
    return;
}


// Perform collision
void collide(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int a=0; a<q; a++) {
                f_eq[a][i][j] = w[a]*rho[i][j]*(1.0+3.0*(ux[i][j]*ex[a]+uy[i][j]*ey[a])+4.5*pow(ux[i][j]*ex[a]+uy[i][j]*ey[a],2.0)-1.5*(pow(ux[i][j],2.0)+pow(uy[i][j],2.0)));
                double omega = -(f[a][i][j]- f_eq[a][i][j])/tau;
                double t1 = ex[a]-ux[i][j];
                double t2 = ex[a]*ux[i][j]+ey[a]*uy[i][j];
                double t3 = ey[a]-uy[i][j];
                double s1  = (1.0-0.5/tau)*w[a]*(3.0*t1+9.0*t2*ex[a])*Fx[i][j];
                double s2  = (1.0-0.5/tau)*w[a]*(3.0*t3+9.0*t2*ey[a])*Fy[i][j];
                f_temp[a][i][j] = f[a][i][j]+omega+s1+s2;
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
                
                int ix = i+(int)ex[a];
                int iy = j+(int)ey[a];
                
                // Periodic condition
                if (ix < 0) {
                    ix = nx-1;
                }
                if (iy < 0) {
                     iy = ny-1;
                 }
                if (ix > nx-1) {
                     ix = 0;
                 }
                if (iy > ny-1) {
                     iy = 0;
                 }
                
                // Streaming at all locations
                
                f[a][ix][iy] = f_temp[a][i][j];
            }
        }
    }

    return;
}


// Check convergence

void convergence_error(){
    // Intialise temporary vaiables
    if (n%n_skip == 0){
        double temp1 = 0;
        double temp2 = 0;
        double temp3 = 0;
        double temp4 = 0;
        
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                
                temp1 += pow(u_abs[i][j]-u_abs_temp[i][j], 2.0);
                temp2 += pow(u_abs[i][j], 2.0);
                
                temp3 += pow(rho[i][j]-rho_temp[i][j], 2.0);
                temp4 += pow(rho[i][j], 2.0);
            }
        }
        err_vel = pow(temp1/temp2, 0.5);
        err_rho = pow(temp3/temp4, 0.5);
        // Calculate convergence error once in N iterations
        cout<<"[# Iter "<< n <<"] Con. error in velocity is " << err_vel << "\t";
        cout<<" Con. error in density is " << err_rho << endl;
    }
    return;
}


void write_data(){
    if (n%n_save == 0) {
            ofstream testfile_1;
        testfile_1.open("/Users/sanketh/github/LBM/Assignment_6/Assignment_6/Assignment_6/data/R_"+to_string(r)+"_rho_"+to_string(n)+".dat");
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                testfile_1 << rho[i][j]<< "\t";
            }
            testfile_1 << endl;
        }
        testfile_1.close();
        
        ofstream testfile_2;
        testfile_2.open("/Users/sanketh/github/LBM/Assignment_6/Assignment_6/Assignment_6/data/R_"+to_string(r)+"_u_"+to_string(n)+".dat");
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                testfile_2 << ux[i][j]<< "\t";
            }
            testfile_2 << endl;
        }
        testfile_2.close();
        
        ofstream testfile_3;
        testfile_3.open("/Users/sanketh/github/LBM/Assignment_6/Assignment_6/Assignment_6/data/R_"+to_string(r)+"_v_"+to_string(n)+".dat");
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                testfile_3 << uy[i][j]<< "\t";
            }
            testfile_3 << endl;
        }
        testfile_3.close();

    }
    }









