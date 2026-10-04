"""Third factorial moment of the MMPP(2) counting process, i.e. the xg3t
expression that mmpp2_fitc uses to solve for h and that mmpp2_fitc_approx puts
in its objective."""
from sage.all import SR, var, matrix, vector, exp, integrate, limit, oo, assume as sage_assume

l1, l2, r1, r2, t, s, u, v = var('l1 l2 r1 r2 t s u v')
for w in (l1, l2, r1, r2, t, s, u, v):
    sage_assume(w > 0)

d = r1 + r2
Q = matrix(SR, [[-r1, r1], [r2, -r2]])
L = matrix(SR, [[l1, 0], [0, l2]])
one = vector(SR, [1, 1])
pi = vector(SR, [r2 / d, r1 / d])
lam = (pi * L) * one


def expQ(w):
    e = exp(-d * w)
    return matrix(SR, [[(r2 + r1 * e) / d, (r1 - r1 * e) / d],
                       [(r2 - r2 * e) / d, (r1 + r2 * e) / d]])


integ = ((pi * expQ(v)) * L * expQ(u - v) * L * expQ(s - u) * L) * one
i1 = integrate(integ.simplify_full(), v, 0, u)
i2 = integrate(i1.simplify_full(), u, 0, s)
m3fact = 6 * integrate(i2.simplify_full(), s, 0, t)

i1b = integrate((((pi * expQ(u)) * L * expQ(s - u) * L) * one).simplify_full(), u, 0, s)
m2fact = 2 * integrate(i1b.simplify_full(), s, 0, t)
EN = lam * t
VarN = (m2fact + EN - EN ** 2).simplify_full()
IDC = (VarN / EN).simplify_full()
binf = limit(IDC, t=oo).simplify_full()

# E[N^3] from the factorial moments, then the third CENTRAL moment
EN2 = m2fact + EN
EN3 = m3fact + 3 * m2fact + EN
m3c = (EN3 - 3 * EN2 * EN + 2 * EN ** 3).simplify_full()

# mmpp2_fitc / _approx claim (xg3t and the m3t2 relation)
xa = lam
p = (l1 - l2) * (r1 - r2)
xg3t = (xa ** 3 * t ** 3 + 3 * xa ** 2 * (binf - 1) * t ** 2
        + 3 * xa * (binf - 1) / d * (p / d - xa) * t
        + 3 * xa / d ** 2 * (binf - 1) * (p + xa * d) * t * exp(-t * d)
        - 6 * xa / d ** 3 * (binf - 1) * p * (1 - exp(-t * d)))
xm3 = xg3t - 3 * xa * t * (xa * t - 1) * IDC - xa * t * (xa * t - 1) * (xa * t - 2)

diff = (m3c - xm3).simplify_full()
try:
    ok = bool(diff == 0)
except Exception:
    ok = False
print('  third central moment of counts - mmpp2_fitc_approx xm3t2 ->', ok)
if not ok:
    print('  residual (should be 0):', diff)
