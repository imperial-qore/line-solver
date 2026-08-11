"""Minimal independent characteristic set per order, and the size of the
pre-computed marking inverse that would have to be hard-coded."""
import time
from sage.all import PolynomialRing, QQ, Matrix, vector

for n in (2, 3, 4):
    t0 = time.time()
    z = n + 1
    names = ['h%d' % i for i in range(n)] + ['r%d' % i for i in range(n - 1)] + ['s']
    names += ['q%d' % j for j in range(z)]
    R = PolynomialRing(QQ, names, order='lex')
    F = R.fraction_field()
    g = [F(x) for x in R.gens()]
    h = g[:n]; r = g[n:2 * n - 1]; s = g[2 * n - 1]
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
        Dc[i, 0] = D1[i, 0] * g[2 * n + i]
    Dc[n - 1, n - 1] = D1[n - 1, n - 1] * g[2 * n + n]

    one = vector(F, [1] * n)
    A = (-D0).inverse()
    P = A * D1
    T = (P.transpose() - Matrix.identity(F, n)); T[n - 1, :] = Matrix(F, 1, n, [1] * n)
    pie = T.solve_right(vector(F, [0] * (n - 1) + [1]))

    # greedy: smallest (a+b) first, keep a characteristic when it raises the rank
    cand = sorted([(a, b) for a in range(1, n + 3) for b in range(0, n + 3)],
                  key=lambda ab: (ab[0] + ab[1], ab[0]))
    rows, chosen = [], []
    for (a, b) in cand:
        y = (pie * (A ** a) * Dc * (A ** b)) * one
        row = [F(F(y).derivative(R.gen(2 * n + j))) for j in range(z)]
        M = Matrix(F, rows + [row])
        if M.rank() > len(rows):
            rows.append(row)
            chosen.append((a, b))
        if len(rows) == z:
            break
    M = Matrix(F, rows)
    t_build = time.time() - t0
    t1 = time.time()
    Mi = M.inverse()
    sizes = [len(str(F(Mi[i, j]))) for i in range(z) for j in range(z)]
    print('n = %d  z = %d  characteristics (backward a, forward b): %s' % (n, z, chosen), flush=True)
    print('        build %.1fs  inverse %.1fs  entries: max %d chars, total %d chars'
          % (t_build, time.time() - t1, max(sizes), sum(sizes)), flush=True)
