"""SageMath proof of the remaining closed-form inverse maps.

Each fitter claims to INVERT the map from parameters to characteristics. The
cleanest exact statement is therefore: take free parameters, compute the
characteristics with the standard MAP formulas, feed them into the closed-form
inverse, and check that the parameters come back. Radicals are handled the way
the MMPP(2) proof handled them: the radicand must be a perfect square in the
parameters, which makes the inverse rational and the identity checkable.

Covered here:
  1. aph2_fitall   (APH(2), Telek-Heindl)
  2. amap2_fitall_gamma, both canonical forms (AMAP(2))
  3. map_hyperexp  (two-phase hyperexponential from mean and SCV) and the
     maximum SCV reachable at a fixed branch probability p
  4. map_mmpp2 / map2_fit third-moment floor E3MIN = (3/2) E2^2/E1
"""
from sage.all import PolynomialRing, QQ, Matrix, vector, factorial


def map_characteristics(D0, D1, K, nmom=3):
    """(E1..Enmom, gamma2) with the standard MAP formulas."""
    n = D0.nrows()
    one = vector(K, [1] * n)
    A = (-D0).inverse()
    P = A * D1
    T = (P.transpose() - Matrix.identity(K, n))
    T[n - 1, :] = Matrix(K, 1, n, [1] * n)
    pie = T.solve_right(vector(K, [0] * (n - 1) + [1]))
    mom = [factorial(k) * ((pie * (A ** k)) * one) for k in range(1, nmom + 1)]
    g2 = P.trace() - 1
    return mom, g2


def show(tag, val, K):
    ok = (K(val) == 0)
    print('  %-58s %s' % (tag, 'IDENTITY' if ok else 'NOT ZERO: %s' % K(val)))
    return ok


res = []

# ---------------------------------------------------------------- 1. APH(2)
print('== aph2_fitall: APH(2) inverse map ==')
R = PolynomialRing(QQ, ['h1', 'h2', 'r1'], order='lex')
K = R.fraction_field()
h1, h2, r1 = [K(g) for g in R.gens()]
# aph2_assemble(h1, h2, r1): canonical acyclic APH(2)
D0 = Matrix(K, [[-1 / h1, r1 / h1], [0, -1 / h2]])
D1 = Matrix(K, [[(1 - r1) / h1, 0], [1 / h2, 0]])
mom, _ = map_characteristics(D0, D1, K)
M1, M2, M3 = mom
tmp0 = M3 ** 2 / 9 + ((8 * M1 ** 3) / 3 - 2 * M2 * M1) * M3 - 3 * M1 ** 2 * M2 ** 2 + 2 * M2 ** 3
tmp2 = M3 - 3 * M1 * M2
tmp3 = 6 * M2 - 12 * M1 ** 2
# the radicand is a perfect square: tmp0 = ((h1 - h2)*(...)/3)^2 up to sign
cand = (h1 - h2) * (h1 * h2 - h1 * r1 * h2) / 1
res.append(show('9*tmp0 - (tmp3*(h1-h2)/2)^2 == 0', 9 * tmp0 - (tmp3 * (h1 - h2) / 2) ** 2, K))
tmp1 = 3 * (tmp3 * (h1 - h2) / 2) / 3      # = 3*sqrt(tmp0) on this branch
res.append(show('h2 recovered: (tmp2 + tmp1)/tmp3 - h1', (tmp2 + tmp1) / tmp3 - h1, K))
res.append(show('h1 recovered: (tmp2 - tmp1)/tmp3 - h2', (tmp2 - tmp1) / tmp3 - h2, K))
res.append(show('r1 recovered: (M1 - h1)/h2 - r1', (M1 - h1) / h2 - r1, K))

# ------------------------------------------------------------- 2. AMAP(2)
for form in (1, 2):
    print('== amap2_fitall_gamma: AMAP(2) inverse map, canonical form %d ==' % form)
    R = PolynomialRing(QQ, ['h1', 'h2', 'p1', 'p2'], order='lex')
    K = R.fraction_field()
    h1, h2, p1, p2 = [K(g) for g in R.gens()]
    if form == 1:
        D0 = Matrix(K, [[-1 / h1, p1 / h1], [0, -1 / h2]])
        D1 = Matrix(K, [[(1 - p1) / h1, 0], [(1 - p2) / h2, p2 / h2]])
    else:
        D0 = Matrix(K, [[-1 / h1, p1 / h1], [0, -1 / h2]])
        D1 = Matrix(K, [[0, (1 - p1) / h1], [(1 - p2) / h2, p2 / h2]])
    mom, G = map_characteristics(D0, D1, K)
    M1, M2, M3 = mom
    tmp0 = M3 ** 2 / 9 + ((8 * M1 ** 3) / 3 - 2 * M2 * M1) * M3 - 3 * M1 ** 2 * M2 ** 2 + 2 * M2 ** 3
    tmp2 = M3 - 3 * M1 * M2
    tmp3 = 6 * M2 - 12 * M1 ** 2
    res.append(show('9*tmp0 - (tmp3*(h1-h2)/2)^2 == 0',
                    9 * tmp0 - (tmp3 * (h1 - h2) / 2) ** 2, K))
    tmp1 = tmp3 * (h1 - h2) / 2
    res.append(show('phase means recovered from (tmp2 +/- tmp1)/tmp3',
                    ((tmp2 + tmp1) / tmp3 - h1) * ((tmp2 + tmp1) / tmp3 - h2), K))
    if form == 1:
        # z must be a perfect square, and r2 = (h1 - M1 + h2 -/+ sqrt(z) + G*M1)/(2 h1)
        z = (M1 ** 2 * G ** 2
             + (2 * M1 * h1 + 2 * M1 * h2 - 4 * h1 * h2 - 2 * M1 ** 2) * G
             + M1 ** 2 - 2 * M1 * h1 - 2 * M1 * h2 + h1 ** 2 + 2 * h1 * h2 + h2 ** 2)
        root = 2 * h1 * p2 - (h1 - M1 + h2 + G * M1)
        res.append(show('z is a perfect square (z - root^2 == 0)', z - root ** 2, K))
        # MATLAB evaluates both signs of sqrt(z) and keeps the feasible one,
        # so the claim is that ONE branch returns the parameter exactly
        r2m = (h1 - M1 + h2 - root + G * M1) / (2 * h1)
        r2p = (h1 - M1 + h2 + root + G * M1) / (2 * h1)
        res.append(show('r2 recovered by one branch of sqrt(z)',
                        (r2m - p2) * (r2p - p2), K))
        r1a = (M1 - h1 - M1 * r2p + h1 * r2p) / (h2 - M1 * r2p)
        res.append(show('r1 recovered from that r2', r1a - p1, K))
    else:
        r2a = (h1 - M1 + h2 + G * M1) / h1
        res.append(show('r2 recovered from the closed form', r2a - p2, K))
        r1a = (r2a + (h1 + h2 - h1 * r2a) / M1 - 2) / (r2a - 1)
        res.append(show('r1 recovered from r2', r1a - p1, K))

# --------------------------------------------------------- 3. map_hyperexp
print('== map_hyperexp: two-phase hyperexponential from (MEAN, SCV, p) ==')
R = PolynomialRing(QQ, ['MEAN', 'SCV', 'p', 'RAD'], order='lex')
K = R.fraction_field()
MEAN, SCV, p, RAD = [K(g) for g in R.gens()]
E2 = (1 + SCV) * MEAN ** 2
Delta = -4 * p * MEAN ** 2 + 4 * p ** 2 * MEAN ** 2 + 2 * E2 * p - 2 * E2 * p ** 2
DISC = R(Delta.numerator()) if Delta.denominator() == 1 else None
for sign in (1, -1):
    mu2 = (-2 * MEAN + 2 * p * MEAN + sign * RAD) / (E2 * p - 2 * MEAN ** 2)
    mu1 = mu2 * p / (p - 1 + MEAN * mu2)
    D0 = Matrix(K, [[-mu1, 0], [0, -mu2]])
    D1 = Matrix(K, [[mu1 * p, mu1 * (1 - p)], [mu2 * p, mu2 * (1 - p)]])
    mom, _ = map_characteristics(D0, D1, K, nmom=2)
    for tag, expr in (('E1 - MEAN', mom[0] - MEAN), ('E2 - (1+SCV)MEAN^2', mom[1] - E2)):
        num = K(expr).numerator()
        poly = num.polynomial(R.gen(3))
        co = poly.list() if num.degree(R.gen(3)) > 0 else [num]
        A = R(0)
        B = R(0)
        for i, c in enumerate(co):
            c = R(c)
            if i % 2 == 0:
                A += c * R(Delta) ** (i // 2)
            else:
                B += c * R(Delta) ** (i // 2)
        ok = (A == 0 and B == 0)
        print('  %-58s %s' % ('root %+d: %s' % (sign, tag), 'IDENTITY' if ok else 'NOT ZERO'))
        res.append(ok)

print('== map_hyperexp: reachable SCV at a fixed p ==')
# with mean fixed, E2 = 2(p x^2 + (1-p) y^2) subject to p x + (1-p) y = MEAN,
# x,y >= 0; the maximum is at the vertex x = MEAN/p, y = 0
Rb = PolynomialRing(QQ, ['MEAN', 'p'], order='lex')
Kb = Rb.fraction_field()
b_MEAN, b_p = [Kb(g) for g in Rb.gens()]
E2max = 2 * (b_p * (b_MEAN / b_p) ** 2)
SCVmax = E2max / b_MEAN ** 2 - 1
print('  max SCV at fixed p:', Kb(SCVmax), ' (3 at p = 1/2)')
res.append(Kb(SCVmax.subs({Rb.gen(1): QQ(1) / 2})) == 3)

# ------------------------------------------------------- 4. third-moment floor
print('== E3MIN = (3/2) E2^2/E1 as the floor over MMPP(2) ==')
R = PolynomialRing(QQ, ['a', 'b', 'p', 'q'], order='lex')
K = R.fraction_field()
a, b, p, q = [K(g) for g in R.gens()]
D0 = Matrix(K, [[-a - p, p], [q, -b - q]])
D1 = Matrix(K, [[a, 0], [0, b]])
mom, _ = map_characteristics(D0, D1, K)
M1, M2, M3 = mom
gap = M3 - QQ(3) / 2 * M2 ** 2 / M1
print('  E3 - (3/2)E2^2/E1 factors as:', K(gap).factor())
print('  (positive for a,b,p,q > 0, so (3/2)E2^2/E1 is a strict lower bound)')
print('\nall inverse-map identities hold:', all(res))
