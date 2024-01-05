from numpy import *

def pfqn_xzgsbup(L = None,N = None,Z = None):
    M = L.shape[0]
    R = Z + sum(L) + amax(L) * (N - 1)
    for i in arange(0,M):
        if L[i] < amax(L):
            R = R + (L[i] - amax(L)) * pfqn_qzgbup(L,N - 1,Z,i)

    X = 2 * N / (R + sqrt(R ** 2 - 4 * Z * amax(L) * N))
    return X


def pfqn_qzgbup(L = None,N = None,Z = None,i = None):
    sigma = sum(L ** 2) / sum(L)
    from . import pfqn_xzabaup
    Yi = L[i] * amin(array([1 / amax(L),N / (Z + sum(L) + sigma * (N - 1 - Z * pfqn_xzabaup(L,N - 1,Z)))]))
    if Yi < 1:
        Qgb = Yi / (1 - Yi) - (Yi ** (N + 1)) / (1 - Yi)
    else:
        Qgb = N

    return Qgb