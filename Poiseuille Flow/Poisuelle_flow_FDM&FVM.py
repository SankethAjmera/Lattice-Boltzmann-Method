import numpy as np
import matplotlib.pyplot as plt
from tdma import tdma

# Solve for POISUELLE flow (Poisson equation)
################### FINITE DIFFERENCE FORMULATION ##################
L = 1 # length of the domain
N = 64; #  number of nodes
h = L/(N-1); # spacing
x = np.linspace(0,L,N) # grid

# Create a tridiagonal matrix
def tridiag(a, b, c, k1, k2, k3): # ki is the ith diagonal relative to principal diagonal
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

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


error = np.sum(np.power(T_fdm-T_actual,2));


###################FINITE VOLUME FORMULATION ##################

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



# Plot
plt.plot(x,T_fdm,'r-')
plt.plot(x,T_fvm,'g^')
plt.plot(x,T_actual,'bo',markersize = 3)
plt.xlabel('position');plt.ylabel('velocity'); plt.axis([0,1,0,1])
plt.legend(('finite difference solution','finite volume solution','exact solution','')); #plt.title('Solution to Couette flow')
m = 10; counter = np.linspace(0,m,m+1);counter = counter.astype(int)# number of labels in Y
plt.grid()
plt.show()


print(error)
