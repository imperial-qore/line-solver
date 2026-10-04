import numpy as np

from scipy.linalg import solve_continuous_lyapunov, inv
from scipy.special import factorial
from scipy import sparse
from scipy.sparse import kron 

from itertools import permutations


def symetric_tensor(T):
    n = len(T.shape)
    newT = np.zeros(shape=(T.shape),dtype=np.float32)
    for sigma in permutations(range(n)):
        newT += np.transpose(T,axes=list(sigma))
    return(newT/factorial(n))

def computeW(Fp,Q):
    return( solve_continuous_lyapunov ( Fp,-Q ) )

def computeV(Fp,Fpp,W):
    return(-np.tensordot(inv(Fp),np.tensordot(Fpp,W/2,2),1))

def computeD(Fp,Q,W):
    n = Q.shape[0]
    QW = symetric_tensor(np.tensordot(W,Q,0))
    QWvec = QW.reshape((n**4))
    Id = sparse.eye(n)
    M = sparse.csr_matrix(kron(kron(kron(Fp,Id),Id),Id)
         +kron(kron(kron(Id,Fp),Id),Id)
         +kron(kron(kron(Id,Id),Fp),Id)
         +kron(kron(kron(Id,Id),Id),Fp))
    Dvec = sparse.linalg.lgmres(M,-6*QWvec, atol=1e-7)
    return(Dvec[0].reshape((n,n,n,n)))

def computeC(Fp,Fpp,Q,Qp,R,D,W,V):
    n = Fp.shape[0]    

    QV = symetric_tensor(np.tensordot(Q,V,0))
    QpW = symetric_tensor(np.tensordot(Qp,W,1))
    FppD = symetric_tensor(np.tensordot(Fpp,D,2))
    
    Id = sparse.eye(n)
    M = sparse.csr_matrix( kron(kron(Fp,Id),Id)
                          + kron(kron(Id,Fp),Id)
                          + kron(kron(Id,Id),Fp) )
    b = (- (3*FppD/2+3*QV+3*QpW+R)).reshape( (n**3))
    Cvec = sparse.linalg.lgmres(M, b, atol=1e-7 )
    return( Cvec[0].reshape((n,n,n)) )

def computeB(Fp,Fpp,Fppp,Qp,Qpp,V,W,D,C):
    FppC = symetric_tensor(np.tensordot(Fpp,C,2))
    FpppD = symetric_tensor(np.tensordot(Fppp,D,3))
    QpV = symetric_tensor(np.tensordot(Qp,V,1))
    QppW = symetric_tensor(np.tensordot(Qpp,W,2))

    b=-( FppC + FpppD/3 + QpV + QppW/2)
    return(solve_continuous_lyapunov ( Fp, b) )

def computeA(Fp,Fpp,Fppp,Fpppp,B,C,D):
    FppB = np.tensordot(Fpp,B,2)
    FpppC = np.tensordot(Fppp,C,3)
    FppppD = np.tensordot(Fpppp,D,4)
    return(-np.tensordot(inv(Fp),FppB/2+FpppC/6+FppppD/24,1) )

def computePiV(pi, Fp,Fpp, Q):
    """ Returns the constants V and W (1/N-term for the steady-state)
    
    This function assumes that Fp is invertible.
    """
    W = computeW(Fp,Q)
    V = computeV(Fp,Fpp,W)
    return(pi,V, (V,W))

def computePiVA(pi, Fp,Fpp,Fppp,Fpppp, Q,Qp,Qpp, R):
    W = computeW(Fp,Q)
    V = computeV(Fp,Fpp,W)

    D = computeD(Fp,Q,W)
    C = computeC(Fp,Fpp,Q,Qp,R,D,W,V)
    B = computeB(Fp,Fpp,Fppp,Qp,Qpp,V,W,D,C)
    A = computeA(Fp,Fpp,Fppp,Fpppp,B,C,D)

    return(pi,V,A, (V,W,A,B,C,D))
