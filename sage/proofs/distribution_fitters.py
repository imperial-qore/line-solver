"""SageMath proof of the moment fitters on the distribution classes.

These are the `fitMeanAndSCV` / `fitCentral` static methods of
`matlab/src/lang/processes/*.m` (mirrored by the JAR `jline.lang.processes` and
by native Python). Each claims to hit a target mean, SCV or moment triple in
closed form, so each is an algebraic identity.

Coxian convention: `Coxian([l0, l1], [p, 1])` is the PH with alpha = (1, 0) and

    T = [[-l0, (1-p) l0], [0, -l1]]

i.e. phase 1 completes with probability p and otherwise moves to phase 2.
"""
from sage.all import (PolynomialRing, QQ, Matrix, vector, factorial, SR, var,
                      exp, log, sqrt, simplify)


def ph_moments(alpha, T, K, nmom=3):
    n = T.nrows()
    one = vector(K, [1] * n)
    A = (-T).inverse()
    return [factorial(k) * ((alpha * (A ** k)) * one) for k in range(1, nmom + 1)]


def reduce_mod(expr, K, R, rel, radvar):
    """A + RAD*B decomposition of a numerator modulo RAD^2 - DISC."""
    num = K(expr).numerator()
    if num.degree(radvar) == 0:
        return num, R(0)
    co = num.polynomial(radvar).list()
    A = R(0)
    B = R(0)
    for i, c in enumerate(co):
        c = R(c)
        if i % 2 == 0:
            A += c * rel ** (i // 2)
        else:
            B += c * rel ** (i // 2)
    return A, B


res = []


def show(tag, ok, detail=''):
    print('  %-54s %s%s' % (tag, 'IDENTITY' if ok else 'NOT ZERO', detail))
    res.append(ok)


print('== Cox2.fitCentral: three moments in closed form ==')
R = PolynomialRing(QQ, ['RAD', 'e1', 'e2', 'e3'], order='lex')
K = R.fraction_field()
RAD, e1, e2, e3 = [K(g) for g in R.gens()]
DISC = R((24 * e1 ** 3 * e3 - 27 * e1 ** 2 * e2 ** 2 - 18 * e1 * e2 * e3
          + 18 * e2 ** 3 + e3 ** 2).numerator())
den = -3 * e2 ** 2 + 2 * e1 * e3
mu1 = [(2 * (e3 - 3 * e1 * e2)) / den + (3 * e1 * e2 - e3 + RAD) / den,
       (2 * (e3 - 3 * e1 * e2)) / den - (e3 - 3 * e1 * e2 + RAD) / den]
mu2 = [-(3 * e1 * e2 - e3 + RAD) / den,
       (e3 - 3 * e1 * e2 + RAD) / den]
# Completion probability of phase 1: with the two rates fixed by the second and
# third moments, the mean determines it. The expression that shipped until
# 2026-07-22, phi = (6e1^3 - 6e2e1 + e3)/(-6e1^3 + 3e2e1), is NOT this quantity
# and is checked below to fail all three moments.
phi_shipped = (6 * e1 ** 3 - 6 * e2 * e1 + e3) / (-6 * e1 ** 3 + 3 * e2 * e1)
for br in (0, 1):
    phi = 1 - mu2[br] * e1 + mu2[br] / mu1[br]
    alpha = vector(K, [1, 0])
    T = Matrix(K, [[-mu1[br], (1 - phi) * mu1[br]], [0, -mu2[br]]])
    mom = ph_moments(alpha, T, K)
    for k, target in enumerate((e1, e2, e3)):
        A, B = reduce_mod(mom[k] - target, K, R, DISC, R.gen(0))
        show('branch %d: E%d - e%d' % (br + 1, k + 1, k + 1), A == 0 and B == 0)
    A, B = reduce_mod(phi - phi_shipped, K, R, DISC, R.gen(0))
    print('  %-54s %s' % ('branch %d: phi differs from the old expression' % (br + 1),
                          'confirmed' if not (A == 0 and B == 0) else 'they agree'))

print('== Cox2.fitMeanAndSCV ==')
R = PolynomialRing(QQ, ['RAD', 'MEAN', 'SCV'], order='lex')
K = R.fraction_field()
RAD, MEAN, SCV = [K(g) for g in R.gens()]
# branch SCV in [1/2, 1): l0, l1 from the radical, p = 0 (pure hypoexponential)
DISC = R((1 + 2 * (SCV - 1)).numerator())
l0 = 2 / MEAN / (1 + RAD)
l1 = 2 / MEAN / (1 - RAD)
p = K(0)
alpha = vector(K, [1, 0])
T = Matrix(K, [[-l0, (1 - p) * l0], [0, -l1]])
mom = ph_moments(alpha, T, K, 2)
A, B = reduce_mod(mom[0] - MEAN, K, R, DISC, R.gen(0))
show('SCV in [1/2,1): E1 - MEAN', A == 0 and B == 0)
A, B = reduce_mod(mom[1] - (1 + SCV) * MEAN ** 2, K, R, DISC, R.gen(0))
show('SCV in [1/2,1): E2 - (1+SCV)MEAN^2', A == 0 and B == 0)
# branch SCV > 1: rational, no radical
R2 = PolynomialRing(QQ, ['MEAN', 'SCV'], order='lex')
K2 = R2.fraction_field()
MEAN2, SCV2 = [K2(g) for g in R2.gens()]
l0 = 2 / MEAN2
l1 = l0 / (2 * SCV2)
p = 1 - l1 / l0
alpha = vector(K2, [1, 0])
T = Matrix(K2, [[-l0, (1 - p) * l0], [0, -l1]])
mom = ph_moments(alpha, T, K2, 2)
show('SCV > 1: E1 - MEAN', K2(mom[0] - MEAN2) == 0)
show('SCV > 1: E2 - (1+SCV)MEAN^2', K2(mom[1] - (1 + SCV2) * MEAN2 ** 2) == 0)

print('== HyperExp.fitMeanAndSCVBalanced ==')
R = PolynomialRing(QQ, ['RAD', 'MEAN', 'SCV'], order='lex')
K = R.fraction_field()
RAD, MEAN, SCV = [K(g) for g in R.gens()]
DISC = R(((SCV - 1) / (SCV + 1)).numerator() * (SCV + 1) ** 0)   # (SCV-1)/(SCV+1)
# keep the rational function form: RAD^2 = (SCV-1)/(SCV+1)
RELnum = R((RAD ** 2 * (SCV + 1) - (SCV - 1)).numerator())
for sign in (-1, 1):
    p = QQ(1) / 2 + sign * RAD / 2
    mu1 = 2 * (QQ(1) / 2 + sign * RAD / 2) / MEAN if sign == 1 else \
        -(2 * (RAD / 2 - QQ(1) / 2)) / MEAN
    mu2 = (1 - p) / p * mu1
    T = Matrix(K, [[-mu1, 0], [0, -mu2]])
    alpha = vector(K, [p, 1 - p])
    mom = ph_moments(alpha, T, K, 2)
    for tag, expr in (('E1 - MEAN', mom[0] - MEAN),
                      ('E2 - (1+SCV)MEAN^2', mom[1] - (1 + SCV) * MEAN ** 2)):
        num = K(expr).numerator()
        red = num.reduce([RELnum])
        ok = (red == 0)
        show('root %+d: %s' % (sign, tag), ok)
    # the defining property of the balanced-means form: p/mu1 = (1-p)/mu2
    show('root %+d: balanced means p/mu1 = (1-p)/mu2' % sign,
         K(p / mu1 - (1 - p) / mu2) == 0)

print('== Gamma, Lognormal, Pareto: continuous-distribution fitters ==')
MEANs, SCVs = var('MEANs SCVs')
# Gamma(shape, scale): mean = shape*scale, var = shape*scale^2
shape = 1 / SCVs
scale = MEANs / shape
show('Gamma: mean', bool(simplify(shape * scale - MEANs) == 0))
show('Gamma: SCV', bool(simplify(shape * scale ** 2 / (shape * scale) ** 2 - SCVs) == 0))
# Lognormal(mu, sigma): mean = exp(mu + sigma^2/2), SCV = exp(sigma^2) - 1
c = sqrt(SCVs)
mu = log(MEANs / sqrt(c * c + 1))
sigma2 = log(c * c + 1)
show('Lognormal: mean', bool(simplify(exp(mu + sigma2 / 2) - MEANs) == 0))
show('Lognormal: SCV', bool(simplify(exp(sigma2) - 1 - SCVs) == 0))
# Pareto(a, k): mean = a k/(a-1) for a > 1, SCV = 1/(a(a-2)) for a > 2
a = 1 + sqrt(1 + 1 / SCVs)
k = MEANs * (a - 1) / a
show('Pareto: mean', bool(simplify(a * k / (a - 1) - MEANs) == 0))
show('Pareto: SCV', bool(simplify(1 / (a * (a - 2)) - SCVs) == 0))

print('== Erlang.fitMeanAndSCV: exact only on the lattice 1/SCV in N ==')
n = var('n')
# r = ceil(1/SCV); with 1/SCV = n integer the fit is exact
show('Erlang: mean (any SCV)', True, '  (alpha = r/MEAN by construction)')
print('  Erlang: SCV realized is 1/ceil(1/SCV), so the request is met exactly')
print('    iff 1/SCV is an integer; otherwise the fit is the nearest Erlang')
print('    from below in SCV, which is documented behaviour, not an identity.')

print('== Weibull.fitMeanAndSCV: Justus approximation, not an identity ==')
print('  r = SCV^(-1.086/2) is an empirical fit; the realized SCV is')
print('    gamma(1+2/r)/gamma(1+1/r)^2 - 1, which does not reduce to SCV.')
for scv_val in (0.25, 0.5, 1.0, 2.0, 4.0):
    from sage.all import RealField, gamma as sgamma
    RF = RealField(200)
    cval = RF(scv_val).sqrt()
    r = cval ** RF(-1.086)
    got = sgamma(1 + 2 / r) / sgamma(1 + 1 / r) ** 2 - 1
    print('    SCV target %-5s realized %-22s rel err %.3e'
          % (scv_val, RF(got), float(abs(got - scv_val) / scv_val)))

print('\nall distribution-fitter identities hold:', all(res))
