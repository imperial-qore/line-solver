"""SageMath proof of map2_fit (Heindl-Horvath-Gross explicit inverse of the MAP(2)).

The fit is written in the normalized moments

    r1 = e1, r2 = e2/2, r3 = e3/6,
    h2 = (r2 - r1^2)/r1^2, h3 = (r3 r1 - r2^2)/r1^4,
    b  = h3 + h2^2 - h2,  c = sqrt(b^2 + 4 h2^3),

and every branch is rational in (e1, h2, h3, g2, c). c is therefore kept as an
indeterminate and each identity is reduced modulo c^2 - (b^2 + 4 h2^3): writing
the numerator as A + c*B, the identity holds on BOTH branches of the square
root iff A = B = 0.

Three branches are checked: the two that share the diagonal-D0 form (hyper,
b >= 0 and b < 0 with g2 >= 0), the correlated hyper branch (b < 0, g2 < 0) and
the hypo branch, which flips the sign of c.
"""
from sage.all import PolynomialRing, QQ, Matrix, vector, factorial

R = PolynomialRing(QQ, ['c', 'e1', 'h2', 'h3', 'g2', 'sq'], order='lex')
K = R.fraction_field()
c, e1, h2, h3, g2, sq = [K(g) for g in R.gens()]

r1 = e1
b = h3 + h2 ** 2 - h2
CREL = R((c ** 2 - (b ** 2 + 4 * h2 ** 3)).numerator())   # c^2 - (b^2 + 4 h2^3)

# target moments implied by (e1, h2, h3)
E1 = e1
E2 = 2 * (h2 * r1 ** 2 + r1 ** 2)
E3 = 6 * ((h3 * r1 ** 4 + (E2 / 2) ** 2) / r1)


def characteristics(D0, D1):
    one = vector(K, [1, 1])
    A = (-D0).inverse()
    P = A * D1
    T = (P.transpose() - Matrix.identity(K, 2))
    T[1, :] = Matrix(K, 1, 2, [1, 1])
    pie = T.solve_right(vector(K, [0, 1]))
    mom = [factorial(k) * ((pie * (A ** k)) * one) for k in (1, 2, 3)]
    return mom, P.trace() - 1


def check(tag, expr, extra_rel=None):
    """expr must vanish modulo c^2 - (b^2 + 4h2^3) (and optionally sq^2 + h3)."""
    num = K(expr).numerator()
    rels = [CREL] if extra_rel is None else [CREL, extra_rel]
    red = num.reduce(rels)
    poly = red.polynomial(R.gen(0)) if red.degree(R.gen(0)) > 0 else None
    ok = (red == 0)
    if not ok and poly is not None:
        ok = all(co == 0 for co in poly.list())
    print('  %-52s %s' % (tag, 'IDENTITY' if ok else 'NOT ZERO'))
    if not ok:
        print('       residual:', red)
    return ok


res = []

print('== map2_fit hyper branch (diagonal D0) ==')
D0 = (1 / (2 * r1 * h3)) * Matrix(K, [[-(2 * h2 + b - c), 0], [0, -(2 * h2 + b + c)]])
D1 = (1 / (4 * r1 * h3)) * Matrix(K, [
    [(2 * h2 + b - c) * (1 - b / c + g2 * (1 + b / c)), (2 * h2 + b - c) * (1 + b / c) * (1 - g2)],
    [(2 * h2 + b + c) * (1 - b / c) * (1 - g2), (2 * h2 + b + c) * (1 + b / c + g2 * (1 - b / c))]])
mom, gam = characteristics(D0, D1)
res.append(check('E1 - e1', mom[0] - E1))
res.append(check('E2 - e2', mom[1] - E2))
res.append(check('E3 - e3', mom[2] - E3))
res.append(check('gamma2 - g2', gam - g2))

print('== map2_fit correlated branch (b < 0, g2 < 0) ==')
a = (h3 + h2 ** 2) / h2
d1 = ((1 - a) * (2 * h2 * g2 + b - c) + g2 * (b + c) - (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
d2 = ((g2 - 1) * (b - c)) / ((1 - a) * (2 * h2 + b - c) + 2 * c)
D0 = (1 / (2 * r1 * h3)) * Matrix(K, [[-(2 * h2 + b - c), (2 * h2 + b - c) * (1 - a)],
                                      [0, -(2 * h2 + b + c)]])
D1 = (1 / (2 * r1 * h3)) * Matrix(K, [[(2 * h2 + b - c) * d1, (2 * h2 + b - c) * (a - d1)],
                                      [(2 * h2 + b + c) * d2, (2 * h2 + b + c) * (1 - d2)]])
mom, gam = characteristics(D0, D1)
res.append(check('E1 - e1', mom[0] - E1))
res.append(check('E2 - e2', mom[1] - E2))
res.append(check('E3 - e3', mom[2] - E3))
res.append(check('gamma2 - g2', gam - g2))

print('== map2_fit hypo branch (c -> -c, a from sqrt(-h3)) ==')
# sq stands for sqrt(-h3), so sq^2 + h3 = 0 is the extra relation
SQREL = R((sq ** 2 + h3).numerator())
a = (2 * h2 + b - c) * (h2 + sq) / (2 * h2 * sq)
cc = -c
d1 = ((1 - a) * (2 * h2 * g2 + b - cc) + g2 * (b + cc) - (b - cc)) / ((1 - a) * (2 * h2 + b - cc) + 2 * cc)
d2 = ((g2 - 1) * (b - cc)) / ((1 - a) * (2 * h2 + b - cc) + 2 * cc)
D0 = (1 / (2 * r1 * h3)) * Matrix(K, [[-(2 * h2 + b - cc), (2 * h2 + b - cc) * (1 - a)],
                                      [0, -(2 * h2 + b + cc)]])
D1 = (1 / (2 * r1 * h3)) * Matrix(K, [[(2 * h2 + b - cc) * d1, (2 * h2 + b - cc) * (a - d1)],
                                      [(2 * h2 + b + cc) * d2, (2 * h2 + b + cc) * (1 - d2)]])
mom, gam = characteristics(D0, D1)
res.append(check('E1 - e1', mom[0] - E1, SQREL))
res.append(check('E2 - e2', mom[1] - E2, SQREL))
res.append(check('E3 - e3', mom[2] - E3, SQREL))
res.append(check('gamma2 - g2', gam - g2, SQREL))

print('\nmap2_fit reproduces (e1, e2, e3, g2) on every branch:', all(res))
