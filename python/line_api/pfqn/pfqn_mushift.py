from numpy import *
def pfqn_mushift(mu = None,iset = None):
    # shifts the service rate vector
    M,N = mu.shape
    for i in transpose(iset):
        for m in arange(0,M+1):
            if m == i:
                mushifted[m,arange[1,[N - 1]+1]] = mu(m,arange(2,N+1))
            else:
                mushifted[m,arange[1,[N - 1]+1]] = mu(m,arange(0,(N - 1)+1))

    return mushifted