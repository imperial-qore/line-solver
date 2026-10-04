"""SageMath derivation of the MMPP(2) counting-process formulas used by
mmpp2_fitc and mmpp2_fitc_approx.

Both fitters work from the index of dispersion for counts

    IDC(t) = Var[N(t)] / E[N(t)]

and mmpp2_fitc inverts it through the Lambert W function on the claim that

    IDC(t) = binf - (binf - 1) * (1 - exp(-x t)) / (x t),      x = r1 + r2

with binf the limit as t -> infinity. That shape is what makes c = (binf-1) /
(binf-bt1) invertible; if it is wrong the whole fit is wrong even though its
round trip is self-consistent.

Var[N(t)] is derived here from first principles for the MMPP(2):

    E[N] = lambda t
    E[N(N-1)] = 2 int_0^t int_0^s pi exp(Q u) L exp(Q (s-u)) L 1 du ds

with Q the modulating generator, L = diag(l1, l2) and pi its stationary vector,
integrated in closed form by Maxima.
"""
from sage.all import SR, var, matrix, vector, exp, integrate, limit, oo


def sbool(x):
    """truth of a symbolic relation"""
    try:
        return bool(x)
    except Exception:
        return False

from sage.all import assume as sage_assume
l1, l2, r1, r2, t, s, u, x = var('l1 l2 r1 r2 t s u x')
for v in (l1, l2, r1, r2, t):
    sage_assume(v > 0)

Q = matrix(SR, [[-r1, r1], [r2, -r2]])
L = matrix(SR, [[l1, 0], [0, l2]])
one = vector(SR, [1, 1])
pi = vector(SR, [r2 / (r1 + r2), r1 / (r1 + r2)])       # stationary of Q
lam = (pi * L) * one                                     # arrival rate

d = r1 + r2
# exp(Q v) in closed form for the 2-state generator
def expQ(v):
    e = exp(-d * v)
    return matrix(SR, [[(r2 + r1 * e) / d, (r1 - r1 * e) / d],
                       [(r2 - r2 * e) / d, (r1 + r2 * e) / d]])

print('  exp(Q v) check: rows sum to one ->',
      sbool((expQ(u) * one - one).simplify_full() == vector(SR, [0, 0])))
print('  exp(Q v) check: pi is stationary ->',
      sbool(((pi * expQ(u)) - pi).simplify_full() == vector(SR, [0, 0])))

# E[N(t)(N(t)-1)] = 2 int_0^t int_0^s pi exp(Qu) L exp(Q(s-u)) L 1 du ds
integrand = ((pi * expQ(u)) * L * expQ(s - u) * L) * one
inner = integrate(integrand.simplify_full(), u, 0, s)
m2fact = 2 * integrate(inner.simplify_full(), s, 0, t)

EN = lam * t
VarN = (m2fact + EN - EN ** 2).simplify_full()
IDC = (VarN / EN).simplify_full()
print('\n  IDC(t) =', IDC.factor())

binf = limit(IDC, t=oo)
print('  IDC(inf) =', binf.simplify_full().factor())

# the shape mmpp2_fitc inverts
claim = binf - (binf - 1) * (1 - exp(-d * t)) / (d * t)
print('\n  IDC(t) - [binf - (binf-1)(1-exp(-x t))/(x t)] with x = r1+r2 ->',
      sbool((IDC - claim).simplify_full() == 0))

# closed form of binf used by the fitters
binf_claim = 1 + 2 * (l1 - l2) ** 2 * r1 * r2 / ((r1 + r2) ** 2 * (l1 * r2 + l2 * r1))
print('  binf - [1 + 2(l1-l2)^2 r1 r2/((r1+r2)^2 (l1 r2 + l2 r1))] ->',
      sbool((binf - binf_claim).simplify_full() == 0))

# the xbt1 expression optimized by mmpp2_fitc_approx, at factor = 1
xb = ((r1 * (2 * l1 ** 2 * r2 ** 2 * t - 2 * l2 ** 2 * r2 - 2 * l1 ** 2 * r2
             + 2 * l2 ** 2 * r2 ** 2 * t + 4 * l1 * l2 * r2
             + 2 * l1 ** 2 * r2 * exp(-r1 * t - r2 * t)
             + 2 * l2 ** 2 * r2 * exp(-r1 * t - r2 * t)
             - 4 * l1 * l2 * r2 ** 2 * t
             - 4 * l1 * l2 * r2 * exp(-r1 * t - r2 * t))
       + r1 ** 2 * (2 * r2 * t * l1 ** 2 - 4 * r2 * t * l1 * l2 + 2 * r2 * t * l2 ** 2))
      / (t * (r1 + r2) ** 3 * (l1 * r2 + l2 * r1)) + 1)
print('  IDC(t) - xbt1 expression of mmpp2_fitc_approx ->',
      sbool((IDC - xb).simplify_full() == 0))
