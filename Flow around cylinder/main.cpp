//
//  main.cpp
//  Assignment 5 - Flow over a cylinder
//
//  Created by Sanketh on 11/10/19.
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
    
    cout<< "Solving flow over a cylinder for Re = "<< Re<<endl;
    
    initialise();
    detect_nodes();
    //    write_grid();
    //    write_nodes();
    
    do {
        n += 1;
        update();
        collide();
        stream();
        convergence_error();
        force();
    } while (n<n_iter);
    
    write_velocity();
    write_force();
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
            u_abs_temp[i][j] = 0;
            rho[i][j] = 1;
            Fx[0] = 0.0;
            Fy[0] = 0.0;
            Cd[0] = 0.0;
            Cd[0] = 0.0;
            // Distribution functions
            for (int a =0 ; a<q; a++) {
                // Initial condition corresponds to uniform velocity throughout the domain
                f[a][i][j]      = w[a]*rho[i][j]*(1.0+3.0*ex[a]*u_in+4.5*pow(ex[a]*u_in,2.0)-1.5*u_in*u_in);
                f_eq[a][i][j]   = w[a]*rho[i][j]*(1.0+3.0*ex[a]*u_in+4.5*pow(ex[a]*u_in,2.0)-1.5*u_in*u_in);
                f_temp[a][i][j] = f_eq[a][i][j];
            }
        }
    }
}


// Detect locations of interior nodes of the cylinder

void detect_nodes(){
    
    // Interior nodes
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            double dist = pow(i-x_0,2)+pow(j-y_0,2);
            int_nodes[i][j] = 0;
            if (dist <= pow(radius,2.0)) {
                int_nodes[i][j] = 1;
            }
        }
    }
    
    // Boundary nodes
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int a=0; a<q; a++) {
                int ix = i+ex[a];
                int iy = j+ey[a];
                bd_nodes[i][j] = 0;
                if (ix>=0 && ix<nx && iy>=0 && iy<ny && int_nodes[i][j] == 0 && int_nodes[ix][iy] == 1) {
                    bd_nodes[i][j] = 1; // Boundary nodes
                    break;
                }
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
            
            // Update variables at previous time step
            if (n%n_skip == 0){
                u_abs_temp[i][j] = u_abs[i][j];
                for (int a=0; a<q; a++) {
                    f_old[a][i][j] = f[a][i][j];
                }
            }
            
            
            // Update variables at new time step
            ux[i][j] = ux[i][j]/rho[i][j];
            uy[i][j] = uy[i][j]/rho[i][j];
            u_abs[i][j] = pow(pow(ux[i][j],2.0)+pow(uy[i][j],2.0), 0.5);
        }
    }
    return;
}


// Perform collision
void collide(){
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            //            if (int_nodes[i][j] == 0) {
            for (int a=0; a<q; a++) {
                f_eq[a][i][j] = w[a]*rho[i][j]*(1.0+3.0*(ux[i][j]*ex[a]+uy[i][j]*ey[a])+4.5*pow(ux[i][j]*ex[a]+uy[i][j]*ey[a],2.0)-1.5*(pow(ux[i][j],2.0)+pow(uy[i][j],2.0)));
                double omega = -(f[a][i][j]-f_eq[a][i][j])/tau;
                f_temp[a][i][j] = f[a][i][j]+omega;
            }
            //            }
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
                
                // Streaming at all loations initially
                
                if (ix>=0 && ix<nx && iy>=0 && iy<ny) {
                    f[a][ix][iy] = f_temp[a][i][j];
                }
            }
        }
    }
    
    // Overwrite Boundary condition on the surface of the cylinder
    
    for (int i=0; i<nx; i++) {
        for (int j=0; j<ny; j++) {
            for (int a=0; a<q; a++) {
                
                int ix = i+ex[a];
                int iy = j+ey[a];
                
                if (int_nodes[i][j]==0 && int_nodes[ix][iy]==1) {
                    // If the stream directions reach the boundary from outside, apply halfway bounceback
                    f[coll_dir[a]][i][j] = f_temp[a][i][j];
                }else if (int_nodes[i][j]==1 && int_nodes[ix][iy]==0) {
                    // If the stream directions reach the boundary from inside, apply halfway bounceback
                    f[coll_dir[a]][i][j] = f_temp[a][i][j];
                }
            }
        }
    }
    
    
    // Overwrite boundary conditions on four sides
    
    for (int i = 1; i<nx-1; i++) {
        // Top boundary condition (Specular bounceback)
        f[4][i][ny-1] = f_temp[2][i][ny-1];
        f[7][i][ny-1] = f_temp[6][i+1][ny-1];
        f[8][i][ny-1] = f_temp[5][i-1][ny-1];
        
        // Bottom boundary condition (Specular bounceback)
        f[2][i][0] = f_temp[4][i][0];
        f[6][i][0] = f_temp[7][i+1][0];
        f[5][i][0] = f_temp[8][i-1][0];
    }
    
    for (int j = 0; j<ny; j++) {
        // Inlet boundary condition
        double rho_in = 1.0/(1.0-u_in)*(f[0][0][j]+f[2][0][j]+f[4][0][j]+2.0*(f[3][0][j]+f[6][0][j]+f[7][0][j]));
        f[1][0][j] = f[3][0][j] + 2.0/3.0*rho_in*u_in;
        f[5][0][j] = f[7][0][j] + 0.5*(f[4][0][j]-f[2][0][j]) + 1.0/6.0*rho_in*u_in;
        f[8][0][j] = f[6][0][j] - 0.5*(f[4][0][j]-f[2][0][j]) + 1.0/6.0*rho_in*u_in;
        // Outlet boundary condition
        double u_out = -1.0 + (f[0][nx-1][j]+f[2][nx-1][j]+f[4][nx-1][j]+2.0*(f[1][nx-1][j]+f[5][nx-1][j]+f[8][nx-1][j]))/rho_out;
        f[3][nx-1][j] = f[1][nx-1][j] - 2.0/3.0*rho_out*u_out;
        f[6][nx-1][j] = f[8][nx-1][j] + 0.5*(f[4][nx-1][j]-f[2][nx-1][j]) - 1.0/6.0*rho_out*u_out;
        f[7][nx-1][j] = f[5][nx-1][j] - 0.5*(f[4][nx-1][j]-f[2][nx-1][j]) - 1.0/6.0*rho_out*u_out;
    }
    return;
}

// Calcuate force on wall

void force(){
    
    if (n%n_skip == 0) {
        
        counter += 1;
        double Fx_temp;
        double Fy_temp;
        Fx[counter]=0.0;
        Fy[counter]=0.0;
        for (int i=0; i<nx; i++) {
            for (int j=0; j<ny; j++) {
                if (bd_nodes[i][j]==1) {
                    for (int a=0; a<q; a++) {
                        int ix = i+ex[a];
                        int iy = j+ey[a];
                        if (int_nodes[ix][iy]==1) {
                            Fx[counter] += 2.0*ex[a]*(f_temp[a][i][j]-f_temp[coll_dir[a]][ix][iy]);
                            Fy[counter] += 2.0*ey[a]*(f_temp[a][i][j]-f_temp[coll_dir[a]][ix][iy]);
                        }
                    }
                }
            }
        }
        // Caluclate force in integer times
        Fx_temp = (Fx[counter-1]+Fx[counter])/2.0;
        Fy_temp = (Fy[counter-1]+Fy[counter])/2.0;
    
        
        // Caluclate force coefficeints
        Cd[counter]=Fx_temp/(radius*u_in*u_in*rho_out);
        Cl[counter]=Fy_temp/(radius*u_in*u_in*rho_out);
        
        double err_force = abs((Cl[counter]-Cl[counter-1])/Cl[counter]);
        
        cout<< "Con. error in lift force = "<< err_force<<"\t";
        cout<< "Cd = "<< Cd[counter]<<"\t";
        cout<< "Cl = "<< Cl[counter]<<endl;
    }
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
                for (int a=0; a<q; a++) {
                    temp3 += pow(f_old[a][i][j]-f[a][i][j], 2.0);
                    temp4 += pow(f[a][i][j], 2.0);
                    
                }
                temp1 += pow(u_abs[i][j]-u_abs_temp[i][j], 2.0);
                temp2 += pow(u_abs[i][j], 2.0);
            }
        }
        err_vel = pow(temp1/temp2, 0.5);
        err_vel = pow(temp3/temp4, 0.5);
        // Calculate convergence error once in N iterations
        cout<<"[# Iter "<< n <<"] Con. error in velocity is " << err_vel << "\t";
    }
    return;
}


void write_force(){
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/Cl_"+to_string(Re)+".dat");
    for (int j=1; j<arr_size; j++) {
        testfile_1 << Cl[j]<< "\t";
    }
    testfile_1.close();
    
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/Cd_"+to_string(Re)+".dat");
    for (int j=1; j<arr_size; j++) {
        testfile_2 << Cd[j]<< "\t";
    }
    testfile_2.close();
}



void write_nodes(){
    
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/nodes_interior.dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << int_nodes[i][j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/nodes_boundary.dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_2 << bd_nodes[i][j]<< "\t";
        }
        testfile_2 <<"\n";
    }
    testfile_2.close();
}


void write_grid(){
    
    // Y grid
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/y_grid.dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << y[j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
    // X grid
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/x_grid.dat");
    for (int i=0; i<ny; i++) {
        for (int j= 0; j<nx; j++) {
            testfile_2 << x[j]<< "\t";
        }
        testfile_2 <<"\n";
    }
    testfile_2.close();
    
}

void write_velocity(){
    // Y velocity
    ofstream testfile_1;
    testfile_1.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/ux_"+to_string(Re)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_1 << ux[i][j]<< "\t";
        }
        testfile_1 <<"\n";
    }
    testfile_1.close();
    
    
    // X velocity
    ofstream testfile_2;
    testfile_2.open("/Users/sanketh/github/LBM/Assignment_5/Assignment 5/Assignment 5/uy_"+to_string(Re)+".dat");
    for (int i=0; i<nx; i++) {
        for (int j= 0; j<ny; j++) {
            testfile_2 << uy[i][j]<< "\t";
        }
        testfile_2 <<"\n";
    }
    testfile_2.close();
}
