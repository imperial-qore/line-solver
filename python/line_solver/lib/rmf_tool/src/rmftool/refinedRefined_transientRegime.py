import numpy as np
from numpy import tensordot as dot
from itertools import permutations
from scipy.special import factorial

def symetric_tensor(T):
    n = len(T.shape)
    newT = np.zeros(shape=(T.shape),dtype=np.float32)
    for sigma in permutations(range(n)):
        newT += np.transpose(T,axes=list(sigma))
    return(newT/factorial(n))

def drift_r(Fp,Fpp,Q,V,W):
    dV = dot(Fp,V,1) + dot(Fpp,W,2)/2
    dW = 2*symetric_tensor(dot(Fp,W,1))+Q
    return(dV,dW)

def drift_r_vector(X,n,computeF,computeFp,computeFpp,computeQ):
    x = X[0:n]
    V = X[n:2*n]
    W = X[2*n:].reshape((n,n))
    
    if n==1:
        parameters=x[0]
    else:
        parameters = x
    
    F   = computeF(parameters)
    Fp  = computeFp(parameters)
    Fpp = computeFpp(parameters)
    
    Q   = computeQ(parameters)
    
    dV,dW = drift_r(Fp,Fpp,Q,V,W)
    dX = np.zeros(2*n+n*n)
    dX[0:n] = F
    dX[n:2*n] = dV
    dX[2*n:] = dW.reshape(n**2)
    return(dX)
   
def drift_rr(Fp,Fpp,Fppp,Fpppp,V,W,A,B,C,D,Q,Qp,Qpp,R):
    dA = dot(Fp,A,1)+dot(Fpp,B,2)/2+dot(Fppp,C,3)/6+dot(Fpppp,D,4)/24
    dB = symetric_tensor( 2*dot(Fp,B,1) + dot(Fpp,C,2) + dot(Fppp,D,3)/3 + dot(Qp,V,1)+dot(Qpp,W,2)/2)
    dC = symetric_tensor( 3*dot(Fp,C,1) + dot(Fpp,D,2)*3/2 + 3*dot(Q,V,0)+3*dot(Qp,W,1) + R)
    dD = symetric_tensor( 4*dot(Fp,D,1) + 6*dot(Q,W,0) )
    return(dA,dB,dC,dD)

def drift_rr_vector(X,n,computeF,computeFp,computeFpp,computeQ,
                    computeFppp,computeFpppp,computeQp,computeQpp,computeR):
    x = X[0:n]
    V = X[n:2*n]
    W = X[2*n:2*n+n**2].reshape((n,n))
    A = X[2*n+n**2:3*n+n**2]
    B = X[3*n+n**2:3*n+2*n**2].reshape((n,n))
    C = X[3*n+2*n**2:3*n+2*n**2+n**3].reshape((n,n,n))
    D = X[3*n+2*n**2+n**3:].reshape((n,n,n,n))

    if n==1:
        parameters=x[0]
    else:
        parameters = x

    F   = computeF(parameters)
    Fp  = computeFp(parameters)
    Fpp = computeFpp(parameters)
    
    Q   = computeQ(parameters)
    
    
    Fp = computeFp(parameters)
    Fpp = computeFpp(parameters)
    Fppp = computeFppp(parameters)
    Fpppp = computeFpppp(parameters)
    
    Q = computeQ(parameters)
    Qp = computeQp(parameters)
    Qpp = computeQpp(parameters)
    
    R = computeR(parameters)
    
    dA,dB,dC,dD = drift_rr(Fp,Fpp,Fppp,Fpppp,V,W,A,B,C,D,Q,Qp,Qpp,R)
    
    dX = np.zeros(3*n+2*n**2+n**3+n**4)
    dX[0:2*n+n**2] = drift_r_vector(X[0:2*n+n**2],n,computeF,computeFp,computeFpp,computeQ)
    dX[2*n+n**2:3*n+n**2] = dA
    dX[3*n+n**2:3*n+2*n**2] = dB.reshape(n**2)
    dX[3*n+2*n**2:3*n+2*n**2+n**3]=dC.reshape(n**3)
    dX[3*n+2*n**2+n**3:] = dD.reshape(n**4)
    
    return(dX)
