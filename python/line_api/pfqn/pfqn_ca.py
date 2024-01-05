from numpy import *
from ..util import *


def pfqn_ca(L=None, N=None, Z=None):

    if len(L) != 0 and len(L[0]) != 0:
        M, R = L.shape
    else:
        M = 0
        R = N.shape[0]

    if M == 0:
        lGn = - sum(factln(N)) + sum(multiply(N, log(sum(Z, 0))))
        Gn = exp(lGn)
        return Gn, lGn

    if sum(N) == 0:
        Gn = 1
        lGn = 0
        return Gn, lGn

    if amin(N) < 0:
        Gn = 0
        lGn = -inf
        return Gn, lGn

    if Z is None or len(Z) == 0:
        Z = zeros((R))
    elif ndim(Z) == 2:
        Z = sum(Z, 0)

    G = ones((M + 1, prod(N + 1)))
    n = pprod(N)
    while sum(n) != - 1:
        idxn = hashpop(n, N)
        G[0, idxn] = pfqn_pffz(Z, n)
        for m in arange(1, M + 1):
            G[m, idxn] = G[m - 1, idxn]
            for r in arange(0, R):
                if n[r] >= 1:
                    n[r] = n[r] - 1
                    idxn_1r = hashpop(n, N)
                    n[r] = n[r] + 1
                    G[m, idxn] += L[m - 1, r] * G[m, idxn_1r]
        n = pprod(n, N)

    Gn = G[-1, -1]
    lGn = log(Gn)
    return Gn, lGn


def pfqn_pffz(Z=None, n=None):
    # F=FZ(Z,N)

    R = n.size
    if sum(n) == 0:
        f = 1
        return f

    f = 0
    for r in arange(0, R):
        if Z[r] == 0:
            f = 0
            return f
        f = f + log(Z[r]) * n[r]
        f = f - factln(n[r])

    f = exp(f)
    return f
