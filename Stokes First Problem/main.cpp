//
//  Solve Stoke's first porblem
//  Assignment 4
//
//  Created by Sanketh on 20/09/19.
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
    
    do {
        n+=1;
        update();
        collide();
        stream();
        force();
        if (n%50 == 0) {
            cout << "[Iter # "<< n<< "]"<<endl;
        }
    } while (n< n_iter);
    
    write_grid(nx, ny);
    write_velocity(n_iter);
    write_force();
    
    return 0;
}



// Initialise variables
void initialise(){
    
    
    // Halfway grid
    y[0] = dy/2; x[0] = dy/2;
    for (int i=1; i<ny; i++) {
        y[i] = y[i-1]+dy;
        x[i] = y[i];
    }
    
    // Define initial values of paramteres
    
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
        }
    }
}


// Update macroscopic variables
void update(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            rho[i][j] = 0;
            ux[i][j]  = 0;
            uy[i][j]  = 0;
            for (int a=0; a<q; a++) {
                rho[i][j] += f[a][i][j];
                ux[i][j]  += ex[a]*f[a][i][j];
                uy[i][j]  += ey[a]*f[a][i][j];
            }
            u_abs_temp[i][j] = u_abs[i][j];
            ux[i][j] = ux[i][j]/rho[i][j];
            uy[i][j] = uy[i][j]/rho[i][j];
            u_abs[i][j] = pow(pow(ux[i][j],2.0)+pow(uy[i][j],2.0), 0.5);
            
            for (int a=0; a<q; a++) {
                f_old[a][i][j] = f[a][i][j];
            }
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
                f_temp[a][i][j] = f[a][i][j]+omega;
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
                
                if (ix<0){     // Left wall
                    ix = nx-1;
                }
                else if (ix>nx-1){  // Right wall
                    ix = 0;
                }
                
                
                // Bottom wall (Midway bounceback)
                if (iy< 0) {
                    f[2][i][j] = f_temp[4][i][j];
                    f[5][i][0] = f_temp[7][i][0]+6*w[5]*rho[i][j]*u_wall;
                    f[6][i][0] = f_temp[8][i][0]-6*w[6]*rho[i][j]*u_wall;
                    
                // Top boundary (Specular bounceback)
                }else if (iy>ny-1){
                    f[4][i][ny-1] = f_temp[2][i][ny-1];
                    if (i== nx-1){
                        f[7][i][ny-1] = f_temp[6][0][ny-1];
                    }else{
                        f[7][i][ny-1] = f_temp[6][i+1][ny-1];
                    }
                    if (i == 0){
                        f[8][i][ny-1] = f_temp[5][nx-1][ny-1];
                    }else{
                        f[8][i][ny-1] = f_temp[5][i-1][ny-1];
                    }
                }else{
                    f[a][ix][iy] = f_temp[a][i][j];
                }
                
            }
        }
    }
    return;
}

// Calcuate force on wall

void force(){
    // Calculating f equiliruim along direction 6
    double t1 = -u_wall;
    double t2 = pow(t1,2.0);
    double t3 = pow(u_wall, 2.0);
    double f_eq_6 = w[6]*(1+3.0*t1+4.5*t2-1.5*t3);
    // Calculating f equiliruim along direction 5
    t1 = u_wall;
    t2 = pow(t1,2.0);
    t3 = pow(u_wall, 2.0);
    double f_eq_5 = w[5]*(1+3.0*t1+4.5*t2-1.5*t3);
    
    
    for (int i=0; i<nx; i++) {
        double sf_7 = 2.0*(f_temp[7][i][0]-f_eq_5+u_wall/6.0);
        double sf_8 = 2.0*(f_temp[8][i][0]-f_eq_6-u_wall/6.0);
        sf[counter] += (sf_8-sf_7);
    }
    
    counter +=1;
}


void write_force(){
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_4/Assigment_4/Assigment_4/shear_force.dat");
    for (int j=0; j<n_iter; j++) {
        testfile_1 << sf[j]<< "\t";
    }
    testfile_1.close();
}




void write_velocity(int n_iter){
    // X velocity
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_4/Assigment_4/Assigment_4/test"+to_string(n_iter)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << ux[i][j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
}

void write_grid(int nx,int ny){
    // Y grid
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_4/Assigment_4/Assigment_4/y_grid.dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << y[j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
}
