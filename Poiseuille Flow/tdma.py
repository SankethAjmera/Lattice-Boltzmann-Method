import numpy as np
# Solve the system of equaitons CX = B, given A and B
def tdma(C,B):
    N = np.size(C,1)
    A = np.zeros((N,N+1));
    A[:,:-1] = C; A[:,-1] = np.transpose(B) # Combine LHS and RHS A = [C|B]
# Matrix row operations
    A[0,:] = A[0,:]/A[0,0] # first row
    for i in range(1,N):
        A[i,:] = A[i,:]- A[i-1,:]*(A[i,i-1]/A[i-1,i-1]) # row subtraction
        A[i,:] = A[i,:]/A[i,i]                          # normalisation
# Back substitution
    X = np.zeros((N,1)) # solution
    X[-1] = A[-1][-1]   # last element
    for i in range(N-2,-1,-1):
        X[i] = A[i][-1]-A[i][i+1]*X[i+1]
    return X
# T_actual = np.linalg.solve(C,np.transpose(B))
