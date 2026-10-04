"""At order 4 the square system came out rank deficient. Is that the chosen
characteristics, or is one marking direction invisible to every characteristic
that is linear in the class matrix?

For a single class, the linear functionals of D1^(c) available are
    pie A^a D1^(c) A^b 1,   a >= 1, b >= 0
(a is the backward order, b the forward order). Their span is what any
single-class linear characteristic can see. Rank of that span, versus the z
marking fractions, decides identifiability.
"""
import time
from sage.all import PolynomialRing, QQ, Matrix, vector

for n in (2, 3, 4, 5):
    z = n + 1
    names = ['h%d' % i for i in range(n)] + ['r%d' % i for i in range(n - 1)] + ['s']
    names += ['q%d' % j for j in range(z)]
    R = PolynomialRing(QQ, names, order='lex')
    F = R.fraction_field()
    g = [F(x) for x in R.gens()]
    h = g[:n]; r = g[n:2 * n - 1]; s = g[2 * n - 1]; q = g[2 * n:]

    D0 = Matrix(F, n, n, 0)
    for i in range(n):
        D0[i, i] = -1 / h[i]
        if i + 1 < n:
            D0[i, i + 1] = r[i] / h[i]
    D1 = Matrix(F, n, n, 0)
    for i in range(n):
        D1[i, 0] = (1 - (r[i] if i + 1 < n else s)) / h[i]
    D1[n - 1, n - 1] = s / h[n - 1]
    Dc = Matrix(F, n, n, 0)
    for i in range(n):
        Dc[i, 0] = D1[i, 0] * q[i]
    Dc[n - 1, n - 1] = D1[n - 1, n - 1] * q[n]

    one = vector(F, [1] * n)
    A = (-D0).inverse()
    P = A * D1
    T = (P.transpose() - Matrix.identity(F, n)); T[n - 1, :] = Matrix(F, 1, n, [1] * n)
    pie = T.solve_right(vector(F, [0] * (n - 1) + [1]))

    rows = []
    labels = []
    for a in range(1, n + 3):
        for b in range(0, n + 3):
            y = (pie * (A ** a) * Dc * (A ** b)) * one
            rows.append([F(F(y).derivative(R.gen(2 * n + j))) for j in range(z)])
            labels.append((a, b))
    J = Matrix(F, rows)
    rk = J.rank()
    print('n = %d: marking fractions z = %d, span of ALL linear characteristics has rank %d  %s'
          % (n, z, rk, 'identifiable' if rk == z else 'DEFICIENT by %d' % (z - rk)), flush=True)
