"""Heindl-van de Liefvoort house of moments (matlab/src/api/moment,
jline.api.moment, line_solver.api.moment).

Every function of the moment API is re-implemented here VERBATIM from the
MATLAB source (same recursions, same loops, same index conventions) and then
checked against an independent definition:

  * the four combinatorial triangles against Sage's own stirling_number1,
    stirling_number2, the falling/rising factorial expansions and the closed
    form of the Lah numbers;
  * the six moment families against their DEFINITION as expectations. All the
    conversions except the central one are linear in the moment vector, so it
    suffices to verify them on a spanning set of moment vectors: the moment
    vectors of the point masses, (m_0,...,m_n) = (1,x,...,x^n), span QQ^(n+1)
    as x varies (Vandermonde). Proving the identity as a polynomial identity
    in x therefore proves it for EVERY distribution.
  * the central conversion, which is not linear in the moment vector because
    the mean enters nonlinearly, on a general two-atom law with symbolic
    masses and atoms, so that the m1 = m(2) extraction convention is exercised
    as well.

Finally all the transforms are extracted as matrices over QQ and the house is
checked to commute: every path between two vertices gives the same matrix.
"""
from sage.all import (QQ, Matrix, PolynomialRing, binomial, factorial,
                      stirling_number1, stirling_number2, vector)

N = 10  # maximum moment order used in the linear checks

results = []


def check(label, ok, detail=''):
    results.append(bool(ok))
    print('  %-58s %s%s' % (label, 'IDENTITY' if ok else 'FAILED', detail))


# ---------------------------------------------------------------------------
# verbatim port of matlab/src/api/moment
# ---------------------------------------------------------------------------

def moment_stirlingcycle(n):
    sigma = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    sigma[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            sigma[i][j] = (i - 1) * sigma[i - 1][j] + sigma[i - 1][j - 1]
    return sigma


def moment_stirling1(n):
    sigma = moment_stirlingcycle(n)
    s = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        for j in range(i + 1):
            s[i][j] = QQ(-1) ** (i - j) * sigma[i][j]
    return s


def moment_stirling2(n):
    S = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    S[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            S[i][j] = j * S[i - 1][j] + S[i - 1][j - 1]
    return S


def moment_lah(n):
    L = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    L[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            L[i][j] = L[i - 1][j - 1] + (i + j - 1) * L[i - 1][j]
    return L


def _mul(T, v):
    n = len(v) - 1
    return [sum(T[i][k] * v[k] for k in range(n + 1)) for i in range(n + 1)]


def moment_binotrans(x):
    n = len(x) - 1
    return [sum(QQ(-1) ** (i - k) * binomial(i, k) * x[k] for k in range(i + 1))
            for i in range(n + 1)]


def moment_binotransinv(y):
    n = len(y) - 1
    return [sum(binomial(i, k) * y[k] for k in range(i + 1))
            for i in range(n + 1)]


def moment_factorial_from_raw(m):
    return _mul(moment_stirling1(len(m) - 1), m)


def moment_raw_from_factorial(f):
    return _mul(moment_stirling2(len(f) - 1), f)


def moment_upfactorial_from_raw(m):
    return _mul(moment_stirlingcycle(len(m) - 1), m)


def moment_raw_from_upfactorial(fp):
    n = len(fp) - 1
    S = moment_stirling2(n)
    T = [[QQ(-1) ** (i - j) * S[i][j] if j <= i else QQ(0)
          for j in range(n + 1)] for i in range(n + 1)]
    return _mul(T, fp)


def moment_binomial_from_factorial(f):
    return [f[i] / factorial(i) for i in range(len(f))]


def moment_factorial_from_binomial(b):
    return [factorial(i) * b[i] for i in range(len(b))]


def moment_negbinomial_from_upfactorial(fp):
    return [fp[i] / factorial(i) for i in range(len(fp))]


def moment_upfactorial_from_negbinomial(bm):
    return [factorial(i) * bm[i] for i in range(len(bm))]


def moment_factorial_from_upfactorial(fp):
    n = len(fp) - 1
    L = moment_lah(n)
    f = [fp[0] - fp[0] + 1] + [0] * n
    for i in range(1, n + 1):
        f[i] = sum(QQ(-1) ** (i - k) * L[i][k] * fp[k] for k in range(1, i + 1))
    return f


def moment_upfactorial_from_factorial(f):
    n = len(f) - 1
    L = moment_lah(n)
    fp = [f[0] - f[0] + 1] + [0] * n
    for i in range(1, n + 1):
        fp[i] = sum(L[i][k] * f[k] for k in range(1, i + 1))
    return fp


def moment_negbinomial_from_binomial(b):
    n = len(b) - 1
    bm = [b[0] - b[0] + 1] + [0] * n
    for i in range(1, n + 1):
        bm[i] = sum(binomial(i - 1, k - 1) * b[k] for k in range(1, i + 1))
    return bm


def moment_binomial_from_negbinomial(bm):
    n = len(bm) - 1
    b = [bm[0] - bm[0] + 1] + [0] * n
    for i in range(1, n + 1):
        b[i] = sum(QQ(-1) ** (i - k) * binomial(i - 1, k - 1) * bm[k]
                   for k in range(1, i + 1))
    return b


def moment_central_from_raw(m):
    n = len(m) - 1
    m1 = m[1]
    return [sum(QQ(-1) ** (i - k) * binomial(i, k) * m[k] * m1 ** (i - k)
                for k in range(i + 1)) for i in range(n + 1)]


def moment_raw_from_central(mc, m1):
    n = len(mc) - 1
    return [sum(binomial(i, k) * mc[k] * m1 ** (i - k) for k in range(i + 1))
            for i in range(n + 1)]


# ---------------------------------------------------------------------------
# 1. the combinatorial triangles
# ---------------------------------------------------------------------------
Rx = PolynomialRing(QQ, 'x')
x = Rx.gen()

fall = [Rx(1)] + [Rx.prod([x - j for j in range(i)]) for i in range(1, N + 1)]
rise = [Rx(1)] + [Rx.prod([x + j for j in range(i)]) for i in range(1, N + 1)]

print('== combinatorial triangles (n = %d) ==' % N)
sigma, s1, S2, L = (moment_stirlingcycle(N), moment_stirling1(N),
                    moment_stirling2(N), moment_lah(N))

check('sigma(i,j) = unsigned Stirling 1st kind',
      all(sigma[i][j] == stirling_number1(i, j)
          for i in range(N + 1) for j in range(N + 1)))
check('sum_j s(i,j) x^j = x(x-1)...(x-i+1)',
      all(sum(s1[i][j] * x ** j for j in range(N + 1)) == fall[i]
          for i in range(N + 1)))
check('S(i,j) = Stirling 2nd kind',
      all(S2[i][j] == stirling_number2(i, j)
          for i in range(N + 1) for j in range(N + 1)))
check('x^i = sum_j S(i,j) (x)_j',
      all(sum(S2[i][j] * fall[j] for j in range(N + 1)) == x ** i
          for i in range(N + 1)))
check('L(i,j) = (i!/j!) C(i-1,j-1)',
      all(L[i][j] == factorial(i) / factorial(j) * binomial(i - 1, j - 1)
          for i in range(1, N + 1) for j in range(1, i + 1)))
check('x(x+1)...(x+i-1) = sum_j L(i,j) (x)_j',
      all(sum(L[i][j] * fall[j] for j in range(N + 1)) == rise[i]
          for i in range(N + 1)))
check('sigma(i,j) = (-1)^(i-j) s(i,j)',
      all(sigma[i][j] == QQ(-1) ** (i - j) * s1[i][j]
          for i in range(N + 1) for j in range(i + 1)))

# ---------------------------------------------------------------------------
# 2. the twelve conversions, on the point-mass spanning set
# ---------------------------------------------------------------------------
print('\n== conversions vs definition, point mass at x (spans QQ^%d) ==' % (N + 1))
mm = [x ** i for i in range(N + 1)]                       # m_i   = E[X^i]
ff = [fall[i] for i in range(N + 1)]                      # f_i   = E[(X)_i]
fp = [rise[i] for i in range(N + 1)]                      # f+_i  = E[X^(i)]
bb = [fall[i] / factorial(i) for i in range(N + 1)]       # b_i   = E[C(X,i)]
bn = [rise[i] / factorial(i) for i in range(N + 1)]       # b-_i  = E[C(X+i-1,i)]

check('factorial_from_raw(m) = f', moment_factorial_from_raw(mm) == ff)
check('raw_from_factorial(f) = m', moment_raw_from_factorial(ff) == mm)
check('upfactorial_from_raw(m) = f+', moment_upfactorial_from_raw(mm) == fp)
check('raw_from_upfactorial(f+) = m', moment_raw_from_upfactorial(fp) == mm)
check('binomial_from_factorial(f) = b', moment_binomial_from_factorial(ff) == bb)
check('factorial_from_binomial(b) = f', moment_factorial_from_binomial(bb) == ff)
check('negbinomial_from_upfactorial(f+) = b-',
      moment_negbinomial_from_upfactorial(fp) == bn)
check('upfactorial_from_negbinomial(b-) = f+',
      moment_upfactorial_from_negbinomial(bn) == fp)
check('factorial_from_upfactorial(f+) = f',
      moment_factorial_from_upfactorial(fp) == ff)
check('upfactorial_from_factorial(f) = f+',
      moment_upfactorial_from_factorial(ff) == fp)
check('negbinomial_from_binomial(b) = b-',
      moment_negbinomial_from_binomial(bb) == bn)
check('binomial_from_negbinomial(b-) = b',
      moment_binomial_from_negbinomial(bn) == bb)

# ---------------------------------------------------------------------------
# 3. the binomial transform pair: it is the unit shift on moment sequences
# ---------------------------------------------------------------------------
print('\n== binomial transform pair ==')
check('binotrans(E[X^i]) = E[(X-1)^i]',
      moment_binotrans(mm) == [(x - 1) ** i for i in range(N + 1)])
check('binotransinv(E[X^i]) = E[(X+1)^i]',
      moment_binotransinv(mm) == [(x + 1) ** i for i in range(N + 1)])
check('binotransinv o binotrans = id', moment_binotransinv(moment_binotrans(mm)) == mm)
check('binotrans o binotransinv = id', moment_binotrans(moment_binotransinv(mm)) == mm)
inv = moment_binotrans(moment_binotrans(mm)) == mm
print('  %-58s %s' % ('binotrans is an involution (docstring claim)',
                      'TRUE' if inv else 'FALSE (it is the INVERSE pair, not an involution)'))

# ---------------------------------------------------------------------------
# 4. central moments, general two-atom law with symbolic masses
# ---------------------------------------------------------------------------
NC = 6
Rc = PolynomialRing(QQ, ['p', 'y0', 'y1'])
p, y0, y1 = Rc.gens()
law = [(p, y0), (1 - p, y1)]
mraw = [sum(pi * xi ** i for pi, xi in law) for i in range(NC + 1)]
mu = mraw[1]
mcen = [sum(pi * (xi - mu) ** i for pi, xi in law) for i in range(NC + 1)]

print('\n== central moments, law p at y0, 1-p at y1 (n = %d) ==' % NC)
check('central_from_raw(m)_i = E[(X-E X)^i]', moment_central_from_raw(mraw) == mcen)
check('raw_from_central(mc,m1)_i = E[X^i]',
      moment_raw_from_central(mcen, mu) == mraw)
check('m0 = 1 and mc1 = 0', mcen[0] == 1 and mcen[1] == 0)
check('mc2 = m2 - m1^2', mcen[2] == mraw[2] - mraw[1] ** 2)
check('mc3 = m3 - 3 m1 m2 + 2 m1^3',
      mcen[3] == mraw[3] - 3 * mraw[1] * mraw[2] + 2 * mraw[1] ** 3)

# ---------------------------------------------------------------------------
# 5. the house commutes: every path between two vertices is the same map
# ---------------------------------------------------------------------------
print('\n== house of moments commutes (n = %d) ==' % N)


def as_matrix(fun):
    """Matrix of a conversion RESTRICTED to the moment subspace {v_0 = 1}.

    Four of the twelve conversions (f<->f+, b<->b-) pin their zeroth output to
    1 instead of copying v_0, so as maps of the whole space they are affine,
    v -> A v + c with a zero column 0 in A. On the subspace v_0 = 1, which is
    the only place a moment vector lives, they coincide with the linear map
    obtained by putting the offset c into that empty column.
    """
    zero = [QQ(0)] * (N + 1)
    c = fun(zero)
    cols = []
    for k in range(N + 1):
        e = [QQ(1) if j == k else QQ(0) for j in range(N + 1)]
        cols.append([a - b for a, b in zip(fun(e), c)])
    affine = any(ci != 0 for ci in c)
    if affine:
        assert all(v == 0 for v in cols[0]), 'affine map with nontrivial column 0'
        cols[0] = c
    return Matrix(QQ, N + 1, N + 1, lambda i, j: cols[j][i]), affine


E, AFF = {}, []
for key, fun in [('f<-m', moment_factorial_from_raw), ('m<-f', moment_raw_from_factorial),
                 ('fp<-m', moment_upfactorial_from_raw), ('m<-fp', moment_raw_from_upfactorial),
                 ('b<-f', moment_binomial_from_factorial), ('f<-b', moment_factorial_from_binomial),
                 ('bn<-fp', moment_negbinomial_from_upfactorial),
                 ('fp<-bn', moment_upfactorial_from_negbinomial),
                 ('f<-fp', moment_factorial_from_upfactorial),
                 ('fp<-f', moment_upfactorial_from_factorial),
                 ('bn<-b', moment_negbinomial_from_binomial),
                 ('b<-bn', moment_binomial_from_negbinomial)]:
    E[key], aff = as_matrix(fun)
    if aff:
        AFF.append(key)
print('  affine on the whole space (zeroth output pinned to 1): %s' % ', '.join(AFF))

I = Matrix.identity(QQ, N + 1)
for a, b in [('f<-m', 'm<-f'), ('fp<-m', 'm<-fp'), ('b<-f', 'f<-b'),
             ('bn<-fp', 'fp<-bn'), ('f<-fp', 'fp<-f'), ('bn<-b', 'b<-bn')]:
    check('%s o %s = I' % (a, b), E[a] * E[b] == I and E[b] * E[a] == I)

check('f<-fp = (f<-m)(m<-fp)', E['f<-fp'] == E['f<-m'] * E['m<-fp'])
check('bn<-b = (bn<-fp)(fp<-f)(f<-b)',
      E['bn<-b'] == E['bn<-fp'] * E['fp<-f'] * E['f<-b'])
check('b<-bn = (b<-f)(f<-fp)(fp<-bn)',
      E['b<-bn'] == E['b<-f'] * E['f<-fp'] * E['fp<-bn'])
check('(b<-f)(f<-m) = (b<-bn)(bn<-fp)(fp<-m)',
      E['b<-f'] * E['f<-m'] == E['b<-bn'] * E['bn<-fp'] * E['fp<-m'])
check('(fp<-f)(f<-m) = (fp<-bn)(bn<-b)(b<-f)(f<-m)',
      E['fp<-f'] * E['f<-m'] ==
      E['fp<-bn'] * E['bn<-b'] * E['b<-f'] * E['f<-m'])
check('Lah triangle is its own signed inverse (f<->fp)',
      E['fp<-f'] * E['f<-fp'] == I)

print('\nhouse of moments verified:', all(results))
