import numpy as np
import matplotlib.pyplot as plt
from tdma import tdma

# Determine error in FDM and FVM approaches for POISUELLE flow (Poisson equation)
# Create a tridiagonal matrix
def tridiag(a, b, c, k1, k2, k3): # ki is the ith diagonal relative to principal diagonal
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

L = 1 # length of the domain
error_fdm = np.zeros(256-4);
error_fvm = np.zeros(256-4);
for N in range(4,256):

    # FINITE DIFFERENCE FORMULATION

    h = L/(N-1); # spacing
    x = np.linspace(0,L,N) # grid

    # Coefficent matrix
    a = np.ones(N-3);    # upper and lower diagonal
    b = -2*np.ones(N-2); # principal diagonal
    A = tridiag(a, b, a, -1, 0, 1)

    # Constant matrix
    X = -8*np.ones((N-2,1))*h*h; # Poisson equation RHS

    # Solve for variables at nodes
    T_num = tdma(A,X); # TDMA
    T_fdm = np.append(np.append(0,T_num),0) # Add Dirichlet BCs
    T_actual = -4*x*(x-1)

    err_fdm = np.sqrt(np.sum(np.power(T_fdm-T_actual,2)))/N;
    error_fdm[N-4] = err_fdm;

    # FINITE VOLUME FORMULATION

    dy = L/(N-1);   # spacing (dyn = dys = dy)
    dys = dy;
    dyn = dy;

    # Coefficent matrix
    a = np.ones(N-3)/dyn;    # upper diagonal
    b = -np.ones(N-2)*(1/dys+1/dyn); # principal diagonal
    a = np.ones(N-3)/dys;    # lower diagonal
    A = tridiag(a, b, a, -1, 0, 1)

    # Constant matrix
    X = -8*np.ones((N-2,1))*dy; # Poisson equation RHS

    # Solution
    T_num = tdma(A,X); # TDMA
    T_fvm = np.append(np.append(0,T_num),0) # Append BCs

    # Error
    err_fvm = np.sqrt(np.sum(np.power(T_fdm-T_actual,2)))/N;
    error_fvm[N-4] = err_fvm;

# Plot error
plt.plot((error_fdm),'r*')
plt.plot((error_fvm),'g.',markersize = 3)
plt.legend(('finite difference','finite volume'));
plt.xlabel('number of nodes');plt.ylabel('$L_2$  error');
plt.show()
