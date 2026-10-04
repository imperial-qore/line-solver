"""Cumulants and factorial cumulants: the vertex the house of moments was
missing (matlab/src/api/moment, jline.api.moment, line_solver.api.moment).

`moment_cumulant_from_raw` / `moment_raw_from_cumulant` and their factorial
counterparts are ported here VERBATIM from the MATLAB source and checked
against two INDEPENDENT definitions:

  * the generating function. The moments are taken as FREE indeterminates
    m_1,...,m_n, so nothing about a particular law is assumed, and the exact
    exponential generating series is formed over that polynomial ring. The
    cumulants are then read off log(sum m_k s^k/k!) and the moments off
    exp(sum kappa_k s^k/k!), which is the definition;
  * the Leonov-Shiryaev set-partition formulas, summed over SetPartitions(n),
    which is the closed form the recursion is the fast evaluation of.

Because the moments are indeterminates, an identity that holds here holds for
EVERY moment sequence, not just for the sequences of some family of laws.

The factorial cumulants are the same transform applied to the factorial
moments: with z = 1+u the probability generating function is E[(1+u)^N] =
sum_k f_k u^k/k!, so log E[z^N] = sum_k kappa^[]_k (z-1)^k/k! is the same
exponential relation. That is verified rather than assumed, and the Poisson
case (all factorial cumulants beyond the first vanish) is checked symbolically
in lambda.
"""
from sage.all import (QQ, PolynomialRing, PowerSeriesRing, SetPartitions,
                      binomial, factorial)

N = 8  # maximum order

results = []


def check(label, ok):
    results.append(bool(ok))
    print('  %-58s %s' % (label, 'IDENTITY' if ok else 'FAILED'))


# ---------------------------------------------------------------------------
# verbatim port of the implemented recursions
# ---------------------------------------------------------------------------

def moment_cumulant_from_raw(m):
    n = len(m) - 1
    kappa = [0] * (n + 1)
    for i in range(1, n + 1):
        acc = 0
        for k in range(1, i):
            acc += binomial(i - 1, k - 1) * kappa[k] * m[i - k]
        kappa[i] = m[i] - acc
    return kappa


def moment_raw_from_cumulant(kappa):
    n = len(kappa) - 1
    m = [0] * (n + 1)
    m[0] = 1
    for i in range(1, n + 1):
        acc = 0
        for k in range(1, i + 1):
            acc += binomial(i - 1, k - 1) * kappa[k] * m[i - k]
        m[i] = acc
    return m


# moment_factcumulant_from_factorial and moment_factorial_from_factcumulant
# delegate to the two routines above in all three codebases.

# ---------------------------------------------------------------------------
# 1. cumulants from the generating function, with FREE moments
# ---------------------------------------------------------------------------
Rm = PolynomialRing(QQ, ['m%d' % i for i in range(1, N + 1)])
mgen = [Rm(1)] + list(Rm.gens())
Ss = PowerSeriesRing(Rm, 's', default_prec=N + 1)
s = Ss.gen()

egf = Ss(sum(mgen[k] * s ** k / factorial(k) for k in range(N + 1))).add_bigoh(N + 1)
cgf = egf.log()
kappa_gf = [Rm(0)] + [Rm(cgf[k] * factorial(k)) for k in range(1, N + 1)]

print('== cumulants vs log of the moment generating function (n = %d) ==' % N)
kappa_rec = moment_cumulant_from_raw(mgen)
check('kappa_n = n! [s^n] log(sum m_k s^k/k!)',
      all(kappa_rec[i] == kappa_gf[i] for i in range(N + 1)))
check('kappa_1 = m1, kappa_2 = m2-m1^2',
      kappa_rec[1] == mgen[1] and kappa_rec[2] == mgen[2] - mgen[1] ** 2)
check('kappa_3 = m3 - 3 m1 m2 + 2 m1^3',
      kappa_rec[3] == mgen[3] - 3 * mgen[1] * mgen[2] + 2 * mgen[1] ** 3)
check('kappa_4 = m4 - 4 m1 m3 - 3 m2^2 + 12 m1^2 m2 - 6 m1^4',
      kappa_rec[4] == (mgen[4] - 4 * mgen[1] * mgen[3] - 3 * mgen[2] ** 2
                       + 12 * mgen[1] ** 2 * mgen[2] - 6 * mgen[1] ** 4))

# ---------------------------------------------------------------------------
# 2. moments from the cumulants, with FREE cumulants
# ---------------------------------------------------------------------------
Rk = PolynomialRing(QQ, ['k%d' % i for i in range(1, N + 1)])
kgen = [Rk(0)] + list(Rk.gens())
Sk = PowerSeriesRing(Rk, 's', default_prec=N + 1)
sk = Sk.gen()
mgf = Sk(sum(kgen[k] * sk ** k / factorial(k) for k in range(1, N + 1))).add_bigoh(N + 1).exp()
m_gf = [Rk(1)] + [Rk(mgf[k] * factorial(k)) for k in range(1, N + 1)]

print('\n== moments vs exp of the cumulant generating function (n = %d) ==' % N)
m_rec = moment_raw_from_cumulant(kgen)
check('m_n = n! [s^n] exp(sum kappa_k s^k/k!)',
      all(m_rec[i] == m_gf[i] for i in range(N + 1)))
check('raw_from_cumulant o cumulant_from_raw = id',
      moment_cumulant_from_raw(moment_raw_from_cumulant(kgen)) == kgen)
check('cumulant_from_raw o raw_from_cumulant = id',
      moment_raw_from_cumulant(moment_cumulant_from_raw(mgen)) == mgen)

# ---------------------------------------------------------------------------
# 3. the Leonov-Shiryaev partition formulas
# ---------------------------------------------------------------------------
print('\n== Leonov-Shiryaev set-partition formulas (n = %d) ==' % N)
ok_km, ok_mk = True, True
for n in range(1, N + 1):
    acc_k, acc_m = Rm(0), Rm(0)
    for pi in SetPartitions(n):
        r = len(pi)
        prod_m = Rm(1)
        for B in pi:
            prod_m *= mgen[len(B)]
        acc_k += (-1) ** (r - 1) * factorial(r - 1) * prod_m
        prod_k = Rm(1)
        for B in pi:
            prod_k *= kappa_rec[len(B)]
        acc_m += prod_k
    ok_km = ok_km and acc_k == kappa_rec[n]
    ok_mk = ok_mk and acc_m == mgen[n]
check('kappa_n = sum_pi (-1)^(r-1) (r-1)! prod m_|B|', ok_km)
check('m_n = sum_pi prod kappa_|B|', ok_mk)

# ---------------------------------------------------------------------------
# 4. factorial cumulants: same transform on the factorial moments
# ---------------------------------------------------------------------------
print('\n== factorial cumulants vs log of the probability generating function ==')
Rf = PolynomialRing(QQ, ['f%d' % i for i in range(1, N + 1)])
fgen = [Rf(1)] + list(Rf.gens())
Su = PowerSeriesRing(Rf, 'u', default_prec=N + 1)
u = Su.gen()
# E[z^N] with z = 1+u is sum_k E[C(N,k)] u^k = sum_k f_k u^k/k!
pgf = Su(sum(fgen[k] * u ** k / factorial(k) for k in range(N + 1))).add_bigoh(N + 1)
lpgf = pgf.log()
kf_gf = [Rf(0)] + [Rf(lpgf[k] * factorial(k)) for k in range(1, N + 1)]
kf_rec = moment_cumulant_from_raw(fgen)   # what the three codebases call
check('kappa^[]_n = n! [u^n] log(sum f_k u^k/k!)',
      all(kf_rec[i] == kf_gf[i] for i in range(N + 1)))
check('factorial_from_factcumulant o factcumulant_from_factorial = id',
      moment_raw_from_cumulant(kf_rec) == fgen)

Rl = PolynomialRing(QQ, 'lam')
lam = Rl.gen()
pois_f = [lam ** k for k in range(N + 1)]
pois_kf = moment_cumulant_from_raw(pois_f)
check('Poisson(lam): kappa^[]_1 = lam, kappa^[]_n = 0 for n >= 2',
      pois_kf[1] == lam and all(pois_kf[i] == 0 for i in range(2, N + 1)))
check('Poisson(lam): every cumulant equals lam',
      moment_cumulant_from_raw(
          moment_raw_from_cumulant([Rl(0)] + [lam] * N))[1:] == [lam] * N)
# and the two vertices agree: raw moments from the all-lambda cumulants are the
# Touchard polynomials, which is what the Stirling edge gives from f_k = lam^k
stirling2 = [[QQ(0)] * (N + 1) for _ in range(N + 1)]
stirling2[0][0] = QQ(1)
for i in range(1, N + 1):
    for j in range(1, i + 1):
        stirling2[i][j] = j * stirling2[i - 1][j] + stirling2[i - 1][j - 1]
check('Poisson(lam): cumulant vertex agrees with the Stirling edge',
      moment_raw_from_cumulant([Rl(0)] + [lam] * N) ==
      [sum(stirling2[i][j] * pois_f[j] for j in range(N + 1)) for i in range(N + 1)])

print('\ncumulant conversions verified:', all(results))
