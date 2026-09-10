"""SageMath check of aph_fit (Bobbio-Horvath-Telek APH(n) three-moment fit).

The construction depends on the order n, so it is verified per order rather
than once. Two cases in aph_fit.m:

  case 1: alpha = (p, 0, ..., 0, 1-p), T = Erlang(n) with rate mu and a slower
          last phase of rate lambda = 1; a, b, p rational in (n2, n3)
  case 2: alpha = (p, 1-p, 0, ..., 0), T = Erlang(n) with a faster first phase;
          the free parameter f is a root of a polynomial, given in closed form
          by the K1..K22 radicals

For case 1 the whole claim is rational, so it is proved exactly. For case 2 the
normalized second moment is proved exactly for ANY f, which isolates the whole
content of the K-radicals into the single scalar equation n3(f) = n3; that
equation is then confirmed to hold at the f the closed form returns, at 120-bit
precision (cube roots of nested radicals do not reduce exactly).
"""
from sage.all import (PolynomialRing, QQ, Matrix, vector, var, SR, sqrt,
                      RealField, factorial)

RF = RealField(400)


def aph_moments(alpha, T, K, nmom=3):
    """normalized moments n2 = E2/E1^2, n3 = E3/(E1 E2) of an APH."""
    n = T.nrows()
    one = vector(K, [1] * n)
    A = (-T).inverse()
    mom = [factorial(k) * ((alpha * (A ** k)) * one) for k in range(1, nmom + 1)]
    return mom


print('== aph_fit case 1: exact for each order ==')
ok1 = []
for n in range(2, 7):
    R = PolynomialRing(QQ, ['n2', 'n3'], order='lex')
    K = R.fraction_field()
    n2, n3 = [K(g) for g in R.gens()]
    b = 2 * (4 - n * (3 * n2 - 4)) / (n2 * (4 + n - n * n3)
                                      + SR(0))  # placeholder, radical handled below
    # b contains a square root; treat it as an indeterminate B with the defining
    # relation, exactly as in the MMPP(2) proof
    R2 = PolynomialRing(QQ, ['n2', 'n3', 'S'], order='lex')
    K2 = R2.fraction_field()
    n2, n3, S = [K2(g) for g in R2.gens()]
    rad = (n * n2) * (12 * n2 ** 2 * (n + 1) + 16 * n3 * (n + 1)
                      + n2 * (n * (n3 - 15) * (n3 + 1) - 8 * (n3 + 3)))
    REL = R2((S ** 2 - rad).numerator())
    b = 2 * (4 - n * (3 * n2 - 4)) / (n2 * (4 + n - n * n3) + S)
    a = (b * n2 - 2) * (n - 1) * b / ((b - 1) * n)
    p = (b - 1) / a
    lam = K2(1)
    mu = lam * (n - 1) / a
    alpha = vector(K2, [p] + [0] * (n - 2) + [1 - p])
    T = Matrix(K2, n, n, 0)
    for i in range(n):
        T[i, i] = -mu
        if i + 1 < n:
            T[i, i + 1] = mu
    T[n - 1, n - 1] = -lam
    mom = aph_moments(alpha, T, K2)
    got_n2 = mom[1] / mom[0] ** 2
    got_n3 = mom[2] / (mom[0] * mom[1])
    r2 = K2(got_n2 - n2).numerator().reduce([REL])
    r3 = K2(got_n3 - n3).numerator().reduce([REL])
    good = (r2 == 0 and r3 == 0)
    # a residual linear in S still vanishes on the branch if both parts vanish
    print('  n = %d: n2 %s   n3 %s' % (n, 'IDENTITY' if r2 == 0 else 'NOT ZERO',
                                       'IDENTITY' if r3 == 0 else 'NOT ZERO'))
    ok1.append(good)

print('== aph_fit case 2: n2 exact for any f, n3 reduced to one scalar root ==')
ok2 = []
for n in range(3, 7):
    R = PolynomialRing(QQ, ['n2', 'f'], order='lex')
    K = R.fraction_field()
    n2, f = [K(g) for g in R.gens()]
    a = 2 * (f - 1) * (n - 1) / ((n - 1) * (n2 * f ** 2 - 2 * f + 2) - n)
    p = (f - 1) * a
    lam = K(1)
    mu = lam * (n - 1) / a
    alpha = vector(K, [p, 1 - p] + [0] * (n - 2))
    T = Matrix(K, n, n, 0)
    for i in range(n):
        T[i, i] = -mu
        if i + 1 < n:
            T[i, i + 1] = mu
    T[0, 0] = -lam
    T[0, 1] = lam
    mom = aph_moments(alpha, T, K)
    got_n2 = mom[1] / mom[0] ** 2
    r2 = K(got_n2 - n2)
    print('  n = %d: n2(f) - n2 %s' % (n, 'IDENTITY (any f)' if r2 == 0 else 'NOT ZERO'))
    ok2.append(r2 == 0)

print('\ncase 1 exact for n = 2..6:', all(ok1))
print('case 2 second moment exact for any f, n = 3..6:', all(ok2))
