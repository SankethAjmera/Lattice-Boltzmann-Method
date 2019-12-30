//
//  Solve Lid driven cavity flow using LBM
//  Assignment 3
//
//  Created by Sanketh on 06/09/19.
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
        update(n);
        collide();
        stream();
        convergence_error(n);

    } while (err_con>error_max);
    write_grid(nx, ny);
    write(nx,ny,Re);

//    cout<<u_wall;
    
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
void update(int n){
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
            if (n%50 ==0) {
                for (int a=0; a<q; a++) {
                    f_old[a][i][j] = f[a][i][j];
                }
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
                // Left wall
                if (ix < 0) {
                    f[1][0][j] = f_temp[3][0][j];
                    f[5][0][j] = f_temp[7][0][j];
                    f[8][0][j] = f_temp[6][0][j];
                    // Right wall
                }else if (ix>nx-1){
                    f[3][nx-1][j] = f_temp[1][nx-1][j];
                    f[7][nx-1][j] = f_temp[5][nx-1][j];
                    f[6][nx-1][j] = f_temp[8][nx-1][j];
                    // Bottom wall
                }else if (iy<0 && i!= 0 && i!=nx-1){
                    f[2][i][0] = f_temp[4][i][0];
                    f[5][i][0] = f_temp[7][i][0];
                    f[6][i][0] = f_temp[8][i][0];
                    // Top wall
                    // Top left and right corners are made to remain fixed
                }else if (iy>ny-1 && i!= 0 && i!=nx-1){
                    f[4][i][ny-1] = f_temp[2][i][ny-1];
                    f[7][i][ny-1] = f_temp[5][i][ny-1]-6*w[7]*rho[i][j]*u_wall;
                    f[8][i][ny-1] = f_temp[6][i][ny-1]+6*w[8]*rho[i][j]*u_wall;
                    // Interior nodes
                }else{
                    f[a][ix][iy] = f_temp[a][i][j];
                }
            }
        }
    }
    return;
}

// Calculate convergence error using L2 norm

void convergence_error(int n){
    // Intialise temporary vaiables
    if (n%50 == 0){
        double temp1 = 0;
        double temp2 = 0;
        
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                for (int a=0; a<q; a++) {
                    temp1 = temp1+pow(f[a][i][j]-f_old[a][i][j], 2.0);
                    temp2 = temp2+pow(f[a][i][j], 2.0);
                }
            }
        }
        err_con = pow(temp1/temp2, 0.5);
        // Calculate convergence error once in N iterations
        cout<<"[# Iter "<< n <<"] Convergence error is " << err_con << endl;
    }
    return;
}

void write(int nx,int ny,int Re){
    // X velocity
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_3/Assn_3/Assn_3/ux_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string(Re)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << ux[i][j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
    // Y velocity
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_3/Assn_3/Assn_3/uy_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string(Re)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_2 << uy[i][j]<< "\t";
        }
        testfile_2 <<"\n";
    }
    testfile_2.close();
}

void write_grid(int nx,int ny){
    // X grid
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_3/Assn_3/Assn_3/x_"+to_string(nx)+"x"+to_string(ny)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << x[j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
    // Y grid
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_3/Assn_3/Assn_3/y_"+to_string(nx)+"x"+to_string(ny)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_2 << y[i]<< "\t";
        }
        testfile_2 <<"\n";
    }
    testfile_2.close();
}
