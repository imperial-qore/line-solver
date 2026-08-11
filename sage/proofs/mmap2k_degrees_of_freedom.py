"""How many invariants does a general MMAP(2,K) have, and how many does the
marked ACYCLIC canonical form (m3a's MAMAP(2,K)) reach?

Method: numerical rank of the Jacobian of a large invariant vector with respect
to the free parameters, at a random interior point. The rank is the dimension
of the image, i.e. the number of independent characteristics the family can
realize, independent of any hand count of similarity transformations.
"""
import itertools
from sage.all import RealField, matrix, vector, identity_matrix, RDF

RF = RealField(300)


def invariants(D0, D1c, K, nmom=6, nlag=4):
    """moments, lag-k acf, class probabilities, forward/backward moments,
    one- and two-step class transition probabilities, joint class moments"""
    n = 2
    D1 = sum(D1c[1:], D1c[0])
    A = (-D0).inverse()
    P = A * D1
    one = vector(D0.base_ring(), [1] * n)
    M = (P.transpose() - identity_matrix(D0.base_ring(), n))
    M[n - 1, :] = matrix(D0.base_ring(), 1, n, [1] * n)
    pie = M.solve_right(vector(D0.base_ring(), [0] * (n - 1) + [1]))
    out = []
    fact = 1
    for k in range(1, nmom + 1):
        fact *= k
        out.append(fact * ((pie * (A ** k)) * one))
    Pk = identity_matrix(D0.base_ring(), n)
    for k in range(1, nlag + 1):
        Pk = Pk * P
        out.append((pie * A * Pk * A) * one)
    for c in range(K):
        pc = (pie * (A * D1c[c])) * one
        out.append(pc)
        out.append(((pie * (A * A * D1c[c])) * one))          # p_c B_c
        out.append(((pie * (A * D1c[c] * A)) * one))          # p_c F_c
        for d in range(K):
            out.append((pie * (A * D1c[c]) * (A * D1c[d])) * one)      # sigma
            out.append((pie * (A * D1c[c]) * A * (A * D1c[d])) * one)  # joint moment
    return vector(D0.base_ring(), out)


def jac_rank(param_to_mmap, x0, K):
    h = RF(1e-20)
    base = invariants(*param_to_mmap(x0), K=K)
    cols = []
    for i in range(len(x0)):
        xp = list(x0); xp[i] = xp[i] + h
        cols.append((invariants(*param_to_mmap(xp), K=K) - base) / h)
    J = matrix(RDF, [[float(c[r]) for c in cols] for r in range(len(base))])
    sv = J.singular_values()
    tol = max(sv) * 1e-9 if sv else 0
    return sum(1 for s in sv if s > tol)


for K in (1, 2, 3, 4, 5, 6):
    # (a) general MMAP(2,K): D0 off-diagonals + all entries of every D1c
    def general(x, K=K):
        p = list(x)
        d01, d10 = p[0], p[1]
        D1c = []
        idx = 2
        for c in range(K):
            D1c.append(matrix(RF, [[p[idx], p[idx + 1]], [p[idx + 2], p[idx + 3]]]))
            idx += 4
        D1 = sum(D1c[1:], D1c[0])
        D0 = matrix(RF, [[0, d01], [d10, 0]])
        for i in range(2):
            D0[i, i] = -(sum(D0[i, j] for j in range(2) if j != i) + sum(D1[i, j] for j in range(2)))
        return D0, D1c
    x0 = [RF(0.3), RF(0.4)] + [RF(0.5 + 0.13 * i) for i in range(4 * K)]
    rg = jac_rank(general, x0, K)

    # (b) marked acyclic canonical form: h1, h2, r1, r2 and 3 fractions per class
    def canon(x, K=K):
        # marking fractions must sum to one in each phase, so only K-1 classes
        # are free; the last one closes the generator
        h1, h2, r1, r2 = x[:4]
        q = [list(x[4 + 3 * c: 7 + 3 * c]) for c in range(K - 1)]
        q.append([1 - sum(q[c][j] for c in range(K - 1)) for j in range(3)])
        D0 = matrix(RF, [[-1 / h1, r1 / h1], [0, -1 / h2]])
        D1 = matrix(RF, [[(1 - r1) / h1, 0], [(1 - r2) / h2, r2 / h2]])
        D1c = [matrix(RF, [[D1[0, 0] * q[c][0], 0],
                           [D1[1, 0] * q[c][1], D1[1, 1] * q[c][2]]]) for c in range(K)]
        return D0, D1c
    y0 = [RF(0.7), RF(1.9), RF(0.35), RF(0.45)] + [RF(0.17 + 0.09 * i) for i in range(3 * (K - 1))]
    rc = jac_rank(canon, y0, K)
    print('K = %d:  general MMAP(2,K) parameters %d -> independent invariants %d'
          % (K, 2 + 4 * K, rg))
    print('        marked acyclic canonical parameters %d -> independent invariants %d'
          % (4 + 3 * (K - 1), rc))
    print('        dimension gap: %d' % (rg - rc))
