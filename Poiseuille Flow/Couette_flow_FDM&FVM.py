import numpy as np
import matplotlib.pyplot as plt
from tdma import tdma

# Solve for COUETTE FLOW (Laplace equation)

################## FINITE DIFFERENCE FORMULATION ##################
L = 1      # length of the domain
N =  64   #  number of nodes
Tl = -0.5; # Left Dir BC
Tr = +0.5; # Right Dir BC
h = L/(N-1);   # spacing
x = np.linspace(0,L,N) # grid

# Create a tridiagonal matrix
def tridiag(a, b, c, k1, k2, k3): # ki is the ith diagonal relative to principal diagonal
    return np.diag(a, k1) + np.diag(b, k2) + np.diag(c, k3)

# Coefficent matrix
a = np.ones(N-3);    # upper and lower diagonal
b = -2*np.ones(N-2); # principal diagonal
A = tridiag(a, b, a, -1, 0, 1)

# Constant matrix
X = np.zeros((N-2,1));X[0] = -Tl; X[-1] = -Tr;

# Solution
T_actual = x-0.5;
T_num = tdma(A,X); # TDMA
T_fdm = np.append(np.append(Tl,T_num),Tr) # Append BCs

################### FINITE VOLUME FORMULATION ##################
dy = L/(N-1);   # spacing (dyn = dys = dy)
dys = dy;
dyn = dy;

# Coefficent matrix
a = np.ones(N-3)/dyn;    # upper diagonal
b = -np.ones(N-2)*(1/dys+1/dyn); # principal diagonal
a = np.ones(N-3)/dys;    # lower diagonal
A = tridiag(a, b, a, -1, 0, 1)

# Constant matrix
X = np.zeros((N-2,1));X[0] = -Tl; X[-1] = -Tr; X = X/dy;

# Solution
T_actual = x-0.5;
T_num = tdma(A,X); # TDMA
T_fvm = np.append(np.append(Tl,T_num),Tr) # Append BCs


# Plot
plt.plot(x,T_fdm,'r-')
plt.plot(x,T_fvm,'g^')
plt.plot(x,T_actual,'bo',markersize = 3)
plt.xlabel('position');plt.ylabel('velocity'); plt.axis([0,1,-0.5,0.5])
plt.legend(('finite difference solution','finite volume solution','exact solution','')); #plt.title('Solution to Couette flow')
m = 10; counter = np.linspace(0,m,m+1);counter = counter.astype(int)# number of labels in Y
y_values = np.linspace(Tl,Tr,m+1); y_values = [round(i,2) for i in y_values]
new_ticks = np.append(np.append("$V_{lower}$", list(y_values[1:m])),"$V_{upper}$")
plt.yticks(y_values,new_ticks)
plt.grid()
plt.show()
