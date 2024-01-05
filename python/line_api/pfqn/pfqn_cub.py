import numpy
from numpy import *
import numpy.matlib as matlib
from scipy.special import factorial

from ..util import *


def pfqn_cub(L=None, N=None, Z=None, order=None):
    # [GN,LGN]=PFQN_CUB(L,N,Z,ORDER)

    # PFQN_GM Exact and approximate solution of closed product-form queueing
    # networks by Grundmann-Moeller cubature rules

    # [Gn,lGn]=pfqn_gm(L,N,Z,S)
    # Input:
    # L : MxR demand matrix. L[i,r] is the demand of class-r at queue i
    # N : 1xR population vector. N[r] is the number of jobs in class r
    # Z : 1xR think time vector. Z[r] is the total think time of class r
    # S : degree of the cubature rule. Exact if S=ceil((sum(N)-1)/2).

    # Output:
    # Gn : estimated normalizing constat
    # lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
    #      will be correctly estimated.

    # Reference:
    # G. Casale. Accelerating performance inference over closed systems by
    # asymptotic methods. ACM SIGMETRICS 2017.
    # Available at: http://dl.acm.org/citation.cfm?id=3084445

    if len(L) != 0 and len(L[0]) != 0:
        M, R = L.shape
    else:
        M = 0
        R = N.shape[1]

    if len(L) == 0 or len(N) == 0 or sum(N) == 0:
        Gn = 1
        return Gn, lGn

    if order is None:
        order = numpy.int8(ceil((sum(N) - 1) / 2))

    if Z is None or sum(Z) == 0:
        Nt = sum(N)
        beta = N / Nt
        f = lambda x=None: matmul(transpose(array([[x], [1 - sum(x)]], dtype=object)), L)
        I, Q, ns = simplexquad(lambda x=None: prod(power(f(x), N)), M - 1, order, 1e-08)
        Gn = Q * exp(factln(sum(N) + M - 1) - sum(factln(N)))
        Gn = Gn[-1]
    else:
        steps = 1e4
        Nt = sum(N)
        beta = N / Nt
        Gn = 0
        vmax = Nt * 10
        tol = 1e-10
        dv = vmax / steps
        for v in arange(0, vmax + dv, dv):
            Lv = L * v + matlib.repmat(Z, M, 1)
            h = lambda x=None: (matmul(transpose(array([[x], [1 - sum(x)]], dtype=object)), Lv))
            I, __, __ = simplexquad(lambda x=None: prod(power(h(x), N)), M - 1, order, tol)
            dG = exp(- v) * (v ** (M - 1)) * I * dv
            Gn = Gn + dG
            if v > 0 and dG < tol * Gn:
                break
        Gn = Gn * exp(- sum(factln(N)))

    lGn = log(Gn)
    return Gn, lGn


def simplexquad(f=None, n=None, order=None, tol=None):
    # [I,Q,NV]=SIMPLEXQUAD(F,N,ORDER,TOL)

    Q, nv = grnmol(lambda x=None: f(x), eye(n, n + 1), order, tol)
    I = Q[-1]
    return I, Q, nv


def grnmol(f=None, V=None, s=None, tol=None):
    # [Q,NV] = GRNMOL( F, V, S , TOL)

    #   Q = grnmol( f, V )
    #     computes approximation to the integral of f over an s-simplex
    #     with vertices as columns of V, an n x (n+1) matrix, using
    #     order 1, 2, ..., s (degree 2s+1) Grundmann-Moler rules.
    #   Output Q is a vector approximations of degree 1, 3, ... 2s+1.
    #   Example: # final two results should be 2/(11x10x9x8x7x6)
    #    n = 4;  grnmol(@(x)x(1)^2*x(n)^5,eye(n,n+1),4)
    #     Reference:
    #       "Invariant Integration Formulas for the N-Simplex by
    #        Combinatorial Methods", A. Grundmann and H. M. Moller,
    #        SIAM J Numer. Anal. 15(1978), pp. 282-290

    n, __ = V.shape  # n is the dimension
    d = 0  # order of the polynomial

    Q = zeros((s + 1, 1))
    Qv = Q.copy()
    Vol = 1 / factorial(n)

    nv = 0
    while 1:

        m = n + 2 * d + 1
        al = ones((n, 1))
        alz = 2 * d + 1
        Qs = 0
        while 1:
            Qs = Qs + f(matmul(V, array([[alz], al], dtype=object)) / m)
            nv = nv + 1
            for j in arange(0, n):
                alz = alz - 2
                if alz > 0:
                    al[j] = al[j] + 2
                    break
                alz = alz + al[j] + 1
                al[j] = 1
            if alz == 2 * d + 1:
                break

        d = d + 1
        Qv[d - 1] = Vol * Qs
        Q[d - 1] = 0
        p = 2 / (prod(2 * array([arange(n + 1, m + 1)])))
        for i in arange(1, d + 1):
            Q[d - 1] = Q[d - 1] + ((m + 2 - 2 * i) ** (2 * d - 1)) * p * Qv[d + 1 - i - 1]
            p = - p * (m + 1 - i) / i
        if d > s or (d > 1 and abs(Q[d - 1] - Q[d - 2]) < tol * Q[d - 2]):
            delete(Q, arange(d, Q.shape[0]), 0)
            break

    return Q, nv
# %%
