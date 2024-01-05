from numpy import *

def pfqn_nrp(L = None,N = None,Z = None,alpha = None,options = None):
        if len(L) != 0 and len(L[0]) != 0:        M, R = L.shape
    else
        M = 0
        R = N.shape[1]
    Nt = sum(N)
    if sum(Z) > 0:
        L = array([[L],[Z]])
        alpha[end() + 1,arange[1,Nt+1]] = arange(0,Nt+1)

    if M == 1 and sum(Z) == 0:
        __,lG = pfqn_gld(L,N,alpha)
        return lG
    else:
        Lmax = amax(L)

    L = L / matlib.repmat(Lmax,L.shape[0],1)

    x0 = zeros((1,R))
    __,__,lG = laplaceapprox(lambda x = None: infradius_hnorm(x,L,N,alpha),x0)
    lG = lG + N * log(transpose(Lmax))
    return lG