from numpy import *

def pfqn_xzgsblow(L = None,N = None,Z = None):
    M = L.shape[0]
    R = Z + sum(L) + amax(L) * (N - 1)
    for i in arange(0,M):
        if L[i] < amax(L):
            R = R + (L[i] - amax(L)) * pfqn_qzgblow(L,N - 1,Z,i)

    X = 2 * N * (1 / (R + sqrt(R ** 2 - 4 * Z * amax(L) * (N - 1))))
    return X


def pfqn_qzgblow(L = None,N = None,Z = None,i = None):
    yi = N * L[i] / (Z + sum(L) + amax(L) * N)
    Qgb = yi / (1 - yi) - (yi ** (N + 1)) / (1 - yi)

    return Qgb