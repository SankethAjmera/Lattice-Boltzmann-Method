//
//  main.cpp
//  Entropic LBM
//
//  Created by Sanketh on 18/11/19.
//  Copyright © 2019 Sanketh. All rights reserved.
//


#include <iostream>
#include <fstream> // for writing files
#include <sstream> // for string handling
#include <cstdlib> // standard library
#include <cmath>   // math library
#include <iomanip> // for setting precision
#include <string>
#include <chrono>
#include "Header.h"

using namespace std::chrono; 
using namespace std;

int main() {
         
         cout << " Entropic Lattice Boltzmann Method with single relaxation time on D2Q9." << endl;
         cout << " Solving Lid driven cavity flow for Re = "<< Re << endl;
         cout << " Velocity of wall is " << u_wall<< endl;
         cout << " Grid size is "  << nx << " x " << ny << endl;
         cout << " Beta value is " << beta << endl;
         cout << " Kinematic viscosity is " << nu << endl;
//         cout << " Running simulation until convergence error of " << error_max <<  "."<<endl;
         cout << " Running simulation for " << n_iter <<  " iterations."<<endl;
         
         // Time counter
           
         // Get starting timepoint
           auto start = high_resolution_clock::now();

         
         // Intialise parameters
         
         initialise();
         
         // Start time loop
         do {
                  n+=1;
                  update();
                  calculate_alpha();
                  collide();
                  stream();
                  convergence_error();
         } while (n < n_iter);
         // End time loop
         
         cout << " Entropic Lattice Boltzmann Method with single relaxation time on D2Q9." << endl;
         cout << " Solved Lid driven cavity flow for Re = "<< Re << endl;
         cout << " Velocity of wall is " << u_wall<< endl;
         cout << " Grid size is "  << nx << " x " << ny << endl;
         cout << " Beta value is " << beta << endl;
         cout << " Kinematic viscosity is " << nu << endl;
         cout << " Simulation complete after " << n << " iterations for Re = " << Re <<"."<<endl;
         
         // Save data
//         write_grid();
//         write_velocity();
         
         
         
           // Get ending timepoint
           auto stop = high_resolution_clock::now();
         
           // Get duration. Substart timepoints to
           // get durarion. To cast it to proper unit
           // use duration cast method
           auto duration = duration_cast<microseconds>(stop - start);
         
           cout << "Time taken by simulation: "
                << duration.count() << " microseconds" << endl;
         
         
         return 0;
}


// Initialise variables
void initialise(){
         
         // Define halfway grid
         
         y[0] = dy/2; x[0] = dy/2;
         for (int i=1; i<ny; i++) {
                  y[i] = y[i-1]+dy;
                  x[i] = y[i];
         }
         
         // Define initial values of paramteres
         
         for (int i=0; i<nx; i++) {
                  for (int j=0; j<ny; j++) {
                           // Macro variables
                           
                           rho[i][j] = ro;    // Initialise density
                           ux[i][j] = 0;      // Initialise velocity
                           uy[i][j] = 0;
                           alpha_old[i][j] = 2.0; // Initial guess for alpha
                           
                           // Initialise distribution functions
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
                           rho[i][j] = 0.0;
                           ux[i][j]  = 0.0;
                           uy[i][j]  = 0.0;
                           // Compute density and velocity
                           for (int a=0; a<q; a++) {
                                    rho[i][j] += f[a][i][j];
                                    ux[i][j]  += ex[a]*f[a][i][j];
                                    uy[i][j]  += ey[a]*f[a][i][j];
                           }
                           u_abs_temp[i][j] = u_abs[i][j];
                           ux[i][j] = ux[i][j]/rho[i][j];
                           uy[i][j] = uy[i][j]/rho[i][j];
                           u_abs[i][j] = pow(pow(ux[i][j],2.0)+pow(uy[i][j],2.0), 0.5);
                  }
         }
         
         
         // Compute feq
         
         for (int i=0; i<nx; i++) {
                  for (int j=0; j<ny; j++) {
                           for (int a =0 ; a<q; a++) {
                                    double f_x = (2-pow(1.0+3.0*ux[i][j]*ux[i][j],0.5))*pow(((2.0*ux[i][j]+pow(1.0+3.0*ux[i][j]*ux[i][j],0.5))/(1-ux[i][j])), ex[a]);
                                    
                                    double f_y = (2-pow(1.0+3.0*uy[i][j]*uy[i][j],0.5))*pow(((2.0*uy[i][j]+pow(1.0+3.0*uy[i][j]*uy[i][j],0.5))/(1-uy[i][j])), ey[a]);
                                    
                                    f_eq[a][i][j]   = w[a]*rho[i][j]*f_x*f_y; // feq
                           }
                  }
         }
         
         return;
}

// Calculate the discrete H-values at each location

void entropy_function(int x, int y){
         // Returns the value of H(f),H(f+alpha(feq-f)) and their derivatives for a given alpha
         
         // Calculate H(f) and its derivative
         
         H = 0.0;
         H_der = 0.0;
         
         for (int a = 0; a<q; a++) {
                  H     += f[a][x][y]*log(f[a][x][y]/w[a]);
                  H_der += 1+log(f[a][x][y]/w[a]);
         }
         
         
         // Calculate H(f+alpha(feq-f)) and its derivative at a given location
         
         H_alpha = 0.0;
         H_der_alpha = 0.0;
         
         for (int a = 0; a<q; a++) {
                  double f_new = f[a][x][y]*(1-alpha_old[x][y]) + alpha_old[x][y]*f_eq[a][x][y];
                  H_alpha     += f_new*log(f_new/w[a]);
                  H_der_alpha += (1-alpha_old[x][y])*(1+log(f_new/w[a]));
         }
         
}


// Calculate the value of alpha at each location for a given time iteration

void calculate_alpha(){
         
         double alpha_new = 0.0;
         double error = 1.0;
         
         do {
                  
                  for (int i=0; i<nx; i++) {
                           for (int j=0; j<ny; j++) {
                                    
                                    // Calculate entropy and its derivative
                                    entropy_function(i,j);
                                    // Get the equation
                                    double H_function     = H - H_alpha;
                                    double H_function_der = H_der - H_der_alpha;
                                    alpha_new = alpha_old[i][j] - H_function/H_function_der; // Update via Newton-Raphson
                                    
                                    // Trivial root case
                                    if (alpha_new == 0.0) {
                                             double temp_array[q];
                                             for (int a=0; a<q; a++) {
                                                      temp_array[a] = abs(f[a][i][j]/(f_eq[a][i][j]-f[a][i][j]));
                                             }
                                             double max_alpha = temp_array[0];
                                             // Get maximum guess
                                             for (int a=1; a<q; a++) {
                                                      if (temp_array[a]>temp_array[a-1]) {
                                                               max_alpha = temp_array[a];
                                                      }
                                             }
                                             alpha_new = max_alpha;
                                    }
                                    
                                    
                                    error = abs(alpha_new-alpha_old[i][j]); // Calculate error
                                    alpha_old[i][j] = alpha_new; // Update
                           }
                  }
                  
         } while (error>err_alpha);
         
}


// Perform collision
void collide(){
         for (int i=0; i<nx; i++) {
                  for (int j=0; j<ny; j++) {
                           for (int a=0; a<q; a++) {
                                    // Collision operator
                                    double omega = -(f[a][i][j]- f_eq[a][i][j])*alpha_old[i][j]*beta;
                                    // Collide
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

void convergence_error(){
         // Calculate convergence in velocity
         if (n%50 == 0){
                  double temp1 = 0;
                  double temp2 = 0;
                  
                  for (int i=0; i<nx; i++) {
                           for (int j=0; j<ny; j++) {
                                    temp1 += pow(u_abs[i][j]-u_abs_temp[i][j], 2.0);
                                    temp2 += pow(u_abs[i][j], 2.0);
                           }
                  }
                  err_con = pow(temp1/temp2, 0.5);
                  // Calculate convergence error once in N iterations
                  cout<<" [# Iter "<< n <<", Re = "<< Re<<"] Convergence error is " << err_con << endl;
         }
         return;
}

void write_velocity(){
         // X velocity
         ofstream testfile_1;
         testfile_1.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/2_u_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string((int)Re)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_1 << ux[i][j]<< "\t";
                  }
                  testfile_1 <<"\n";
         }
         testfile_1.close();
         
         // Y velocity
         ofstream testfile_2;
         testfile_2.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/2_v_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string((int)Re)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_2 << uy[i][j]<< "\t";
                  }
                  testfile_2 <<"\n";
         }
         testfile_2.close();
         
         // Density
         ofstream testfile_3;
         testfile_3.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/2_rho_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string((int)Re)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_3 << rho[i][j]<< "\t";
                  }
                  testfile_3 <<"\n";
         }
         testfile_3.close();
         
         // Alpha
         ofstream testfile_4;
         testfile_4.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/2_alpha_"+to_string(nx)+"x"+to_string(ny)+"_Re_"+to_string((int)Re)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_4 << alpha_old[i][j]<< "\t";
                  }
                  testfile_4 <<"\n";
         }
         testfile_4.close();
}


void write_grid(){
         // X grid
         ofstream testfile_1;
         testfile_1.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/x_"+to_string(nx)+"x"+to_string(ny)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_1 << x[j]<< "\t";
                  }
                  testfile_1 <<"\n";
         }
         testfile_1.close();
         
         // Y grid
         ofstream testfile_2;
         testfile_2.open("/Users/sanketh/github/LBM/Project/Entropic LBM/Entropic LBM/y_"+to_string(nx)+"x"+to_string(ny)+".dat");
         for (int i=0; i<nx; i++) {
                  for (int j= 0; j<ny; j++) {
                           testfile_2 << y[i]<< "\t";
                  }
                  testfile_2 <<"\n";
         }
         testfile_2.close();
}
