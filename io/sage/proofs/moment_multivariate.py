"""The multivariate house of moments (matlab/src/api/moment, jline.api.moment,
line_solver.api.moment): the joint conversions of a random vector, the joint
cumulants, and the two counting-process identities behind class marking.

Every routine is ported here VERBATIM from the MATLAB source and checked
against an INDEPENDENT definition:

  * the twelve SEPARABLE conversions, against the definition of the six joint
    moment families. They are linear in the joint moment array, and the arrays
    of the product point masses, m_a = prod_j x_j^(a_j), span the whole tensor
    space as the atoms vary (a Vandermonde per mode), so a polynomial identity
    in (x_1,...,x_d) proves the conversion for EVERY joint law;
  * the joint CENTRAL conversion, on a point mass with the means left FREE.
    Both sides are linear in the law at fixed means, so this again covers every
    joint law once the free means are set to the true ones;
  * the joint CUMULANTS, against the joint cumulant generating function. The
    joint moments are taken as FREE indeterminates and log(sum m_a x^a/a!) is
    formed by exact truncated series arithmetic, so nothing about a particular
    law is assumed. The multivariate Leonov-Shiryaev partition formula is
    checked against the same recursion;
  * the MARKING formula, from the composition rule of the probability
    generating function under i.i.d. multinomial marking, with the aggregate
    factorial moments left free;
  * the AGGREGATION formula, from the Vandermonde convolution of falling
    factorials, again on point masses, so it holds for any joint law of the
    parts, marked or not.
"""
from itertools import product

from sage.all import (QQ, PolynomialRing, SetPartitions, binomial, factorial,
                      prod)

results = []


def check(label, ok):
    results.append(bool(ok))
    print('  %-58s %s' % (label, 'IDENTITY' if ok else 'FAILED'))


def boxes(dims):
    """All multi-indices of the array, in lexicographic order."""
    return list(product(*[range(n) for n in dims]))


# ---------------------------------------------------------------------------
# verbatim port of the implemented routines
# ---------------------------------------------------------------------------

def moment_lah(n):
    L = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    L[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            L[i][j] = L[i - 1][j - 1] + (i + j - 1) * L[i - 1][j]
    return L


def moment_stirling2(n):
    S = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    S[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            S[i][j] = j * S[i - 1][j] + S[i - 1][j - 1]
    return S


def moment_stirlingcycle(n):
    sig = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    sig[0][0] = QQ(1)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            sig[i][j] = (i - 1) * sig[i - 1][j] + sig[i - 1][j - 1]
    return sig


def moment_housematrix(edge, n):
    if edge == 'upfactorial_from_raw':
        return moment_stirlingcycle(n)
    if edge == 'factorial_from_raw':
        sig = moment_stirlingcycle(n)
        return [[QQ(-1) ** (i - j) * sig[i][j] for j in range(n + 1)] for i in range(n + 1)]
    if edge == 'raw_from_factorial':
        return moment_stirling2(n)
    if edge == 'raw_from_upfactorial':
        S = moment_stirling2(n)
        return [[QQ(-1) ** (i - j) * S[i][j] for j in range(n + 1)] for i in range(n + 1)]
    T = [[QQ(0)] * (n + 1) for _ in range(n + 1)]
    if edge in ('binomial_from_factorial', 'negbinomial_from_upfactorial'):
        for i in range(n + 1):
            T[i][i] = QQ(1) / factorial(i)
        return T
    if edge in ('factorial_from_binomial', 'upfactorial_from_negbinomial'):
        for i in range(n + 1):
            T[i][i] = QQ(factorial(i))
        return T
    T[0][0] = QQ(1)
    if edge in ('factorial_from_upfactorial', 'upfactorial_from_factorial'):
        L = moment_lah(n)
        for i in range(1, n + 1):
            for k in range(1, i + 1):
                T[i][k] = (L[i][k] if edge == 'upfactorial_from_factorial'
                           else QQ(-1) ** (i - k) * L[i][k])
        return T
    if edge in ('negbinomial_from_binomial', 'binomial_from_negbinomial'):
        for i in range(1, n + 1):
            for k in range(1, i + 1):
                c = binomial(i - 1, k - 1)
                T[i][k] = (c if edge == 'negbinomial_from_binomial'
                           else QQ(-1) ** (i - k) * c)
        return T
    raise ValueError('unknown edge %s' % edge)


def moment_tensortrans(A, dims, T, mode):
    out = {}
    for a in boxes(dims):
        acc = 0
        for k in range(dims[mode]):
            b = list(a)
            b[mode] = k
            acc += T[a[mode]][k] * A[tuple(b)]
        out[a] = acc
    return out


def moment_jointtrans(A, dims, edge):
    out = A
    for mode in range(len(dims)):
        out = moment_tensortrans(out, dims, moment_housematrix(edge, dims[mode] - 1), mode)
    return out


def moment_joint_central_from_raw_mean(A, dims, mu):
    out = A
    for mode in range(len(dims)):
        n = dims[mode] - 1
        T = [[binomial(i, k) * (-mu[mode]) ** (i - k) if k <= i else 0
              for k in range(n + 1)] for i in range(n + 1)]
        out = moment_tensortrans(out, dims, T, mode)
    return out


def moment_joint_raw_from_central(A, dims, mu):
    out = A
    for mode in range(len(dims)):
        n = dims[mode] - 1
        T = [[binomial(i, k) * mu[mode] ** (i - k) if k <= i else 0
              for k in range(n + 1)] for i in range(n + 1)]
        out = moment_tensortrans(out, dims, T, mode)
    return out


def moment_joint_cumulant_from_raw(m, dims):
    d = len(dims)
    kappa = {}
    for a in boxes(dims):
        if sum(a) == 0:
            kappa[a] = 0
            continue
        j = next(l for l in range(d) if a[l] > 0)
        acc = 0
        for b in product(*[range(ai + 1) for ai in a]):
            if sum(b) == 0 or b == a or b[j] < 1:
                continue
            c = prod([binomial(a[l] - (1 if l == j else 0), b[l] - (1 if l == j else 0))
                      for l in range(d)])
            if c == 0:
                continue
            acc += c * kappa[b] * m[tuple(a[l] - b[l] for l in range(d))]
        kappa[a] = m[a] - acc
    return kappa


def moment_joint_raw_from_cumulant(kappa, dims):
    d = len(dims)
    m = {}
    for a in boxes(dims):
        if sum(a) == 0:
            m[a] = 1
            continue
        j = next(l for l in range(d) if a[l] > 0)
        acc = 0
        for b in product(*[range(ai + 1) for ai in a]):
            if sum(b) == 0 or b[j] < 1:
                continue
            c = prod([binomial(a[l] - (1 if l == j else 0), b[l] - (1 if l == j else 0))
                      for l in range(d)])
            if c == 0:
                continue
            acc += c * kappa[b] * m[tuple(a[l] - b[l] for l in range(d))]
        m[a] = acc
    return m


def moment_joint_marking(f, p, dims):
    return {a: prod([p[j] ** a[j] for j in range(len(p))]) * f[sum(a)]
            for a in boxes([n + 1 for n in dims])}


def moment_joint_aggregate(F, dims):
    nmax = min(dims) - 1
    out = [0] * (nmax + 1)
    for a in boxes(dims):
        n = sum(a)
        if n <= nmax:
            out[n] += factorial(n) / prod([factorial(aj) for aj in a]) * F[a]
    return out


# ---------------------------------------------------------------------------
# 1. the twelve separable conversions on product point masses
# ---------------------------------------------------------------------------
for DIMS in [(5, 5), (4, 3, 3)]:
    d = len(DIMS)
    R = PolynomialRing(QQ, ['x%d' % j for j in range(d)])
    xs = R.gens()
    fall = [[R.prod([xs[j] - t for t in range(i)]) for i in range(DIMS[j])] for j in range(d)]
    rise = [[R.prod([xs[j] + t for t in range(i)]) for i in range(DIMS[j])] for j in range(d)]
    pw = [[xs[j] ** i for i in range(DIMS[j])] for j in range(d)]

    def tens(tab):
        return {a: R.prod([tab[j][a[j]] for j in range(d)]) for a in boxes(DIMS)}

    m = tens(pw)
    f = tens(fall)
    fp = tens(rise)
    b = {a: f[a] / prod([factorial(aj) for aj in a]) for a in boxes(DIMS)}
    bn = {a: fp[a] / prod([factorial(aj) for aj in a]) for a in boxes(DIMS)}

    print('== separable joint conversions, product point mass, dims %s ==' % (DIMS,))
    for edge, src, dst in [('factorial_from_raw', m, f), ('raw_from_factorial', f, m),
                           ('upfactorial_from_raw', m, fp), ('raw_from_upfactorial', fp, m),
                           ('binomial_from_factorial', f, b), ('factorial_from_binomial', b, f),
                           ('negbinomial_from_upfactorial', fp, bn),
                           ('upfactorial_from_negbinomial', bn, fp),
                           ('factorial_from_upfactorial', fp, f),
                           ('upfactorial_from_factorial', f, fp),
                           ('negbinomial_from_binomial', b, bn),
                           ('binomial_from_negbinomial', bn, b)]:
        got = moment_jointtrans(src, DIMS, edge)
        check('joint_%s' % edge, all(got[a] == dst[a] for a in boxes(DIMS)))

    # 2. the central conversion, with FREE means
    Rc = PolynomialRing(QQ, ['x%d' % j for j in range(d)] + ['u%d' % j for j in range(d)])
    xc = Rc.gens()[:d]
    uc = Rc.gens()[d:]
    mc_def = {a: Rc.prod([(xc[j] - uc[j]) ** a[j] for j in range(d)]) for a in boxes(DIMS)}
    m_c = {a: Rc.prod([xc[j] ** a[j] for j in range(d)]) for a in boxes(DIMS)}
    print('== joint central moments, free means, dims %s ==' % (DIMS,))
    got = moment_joint_central_from_raw_mean(m_c, DIMS, uc)
    check('joint_central_from_raw_mean = E[prod (X_j-mu_j)^a_j]',
          all(got[a] == mc_def[a] for a in boxes(DIMS)))
    back = moment_joint_raw_from_central(mc_def, DIMS, uc)
    check('joint_raw_from_central inverts it',
          all(back[a] == m_c[a] for a in boxes(DIMS)))
    if d == 2:
        # the mean-reading entry point, on a genuine two-atom law: only there
        # does the (1,1) entry become the covariance, since a point mass has
        # none. Both sides are linear in the law, so this closes the argument.
        Rl = PolynomialRing(QQ, ['q', 'a0', 'a1', 'b0', 'b1'])
        q, a0, a1, b0, b1 = Rl.gens()
        law = [(q, a0, b0), (1 - q, a1, b1)]
        mlaw = {a: sum(w * u ** a[0] * v ** a[1] for w, u, v in law) for a in boxes(DIMS)}
        mu0, mu1 = mlaw[(1, 0)], mlaw[(0, 1)]
        mclaw = {a: sum(w * (u - mu0) ** a[0] * (v - mu1) ** a[1] for w, u, v in law)
                 for a in boxes(DIMS)}
        gotl = moment_joint_central_from_raw_mean(mlaw, DIMS, [mu0, mu1])
        check('joint_central_from_raw on a two-atom law',
              all(gotl[a] == mclaw[a] for a in boxes(DIMS)))
        check('entry (1,1) of the central array is the covariance',
              gotl[(1, 1)] == mlaw[(1, 1)] - mu0 * mu1)

# ---------------------------------------------------------------------------
# 3. joint cumulants vs the joint cumulant generating function, FREE moments
# ---------------------------------------------------------------------------
for DIMS in [(3, 3), (3, 2, 2)]:
    d = len(DIMS)
    idx = [a for a in boxes(DIMS) if sum(a) > 0]
    Rm = PolynomialRing(QQ, ['m' + '_'.join(str(t) for t in a) for a in idx])
    mfree = dict(zip(idx, Rm.gens()))
    mfree[tuple([0] * d)] = Rm(1)
    P = PolynomialRing(Rm, ['x%d' % j for j in range(d)])
    xs = P.gens()

    def trunc(poly):
        out = P(0)
        for mon, co in poly.dict().items():
            if all(mon[j] < DIMS[j] for j in range(d)):
                out += co * P.prod([xs[j] ** mon[j] for j in range(d)])
        return out

    ser = trunc(sum(mfree[a] / prod([factorial(aj) for aj in a])
                    * P.prod([xs[j] ** a[j] for j in range(d)]) for a in idx))
    logser, pwr = P(0), P(1)
    for r in range(1, sum(DIMS) + 1):
        pwr = trunc(pwr * ser)
        if pwr == 0:
            break
        logser += QQ(-1) ** (r + 1) / r * pwr
    kappa_gf = {a: prod([factorial(aj) for aj in a])
                * logser.monomial_coefficient(P.prod([xs[j] ** a[j] for j in range(d)]))
                for a in idx}

    print('== joint cumulants vs log of the joint mgf, free moments, dims %s ==' % (DIMS,))
    kappa = moment_joint_cumulant_from_raw(mfree, DIMS)
    check('kappa_a = a! [x^a] log(sum m_b x^b/b!)',
          all(kappa[a] == kappa_gf[a] for a in idx))
    check('joint_raw_from_cumulant inverts joint_cumulant_from_raw',
          all(moment_joint_raw_from_cumulant(kappa, DIMS)[a] == mfree[a] for a in boxes(DIMS)))
    e1 = tuple([1] + [0] * (d - 1))
    e12 = tuple([1, 1] + [0] * (d - 2))
    check('kappa_(1,0,..) = m_(1,0,..) and kappa_(1,1,0,..) is the covariance',
          kappa[e1] == mfree[e1] and
          kappa[e12] == mfree[e12] - mfree[e1] * mfree[tuple([0, 1] + [0] * (d - 2))])

    # multivariate Leonov-Shiryaev over the set partitions of the index multiset
    ok_ls, ok_inv = True, True
    for a in idx:
        if sum(a) > 4:
            continue
        labels = []
        for j in range(d):
            labels += [j] * a[j]
        n = len(labels)
        acc_k, acc_m = Rm(0), Rm(0)
        for pi in SetPartitions(n):
            r = len(pi)
            pk, pm = Rm(1), Rm(1)
            for B in pi:
                cnt = [0] * d
                for pos in B:
                    cnt[labels[pos - 1]] += 1
                pm *= mfree[tuple(cnt)]
                pk *= kappa[tuple(cnt)]
            acc_k += QQ(-1) ** (r - 1) * factorial(r - 1) * pm
            acc_m += pk
        ok_ls = ok_ls and acc_k == kappa[a]
        ok_inv = ok_inv and acc_m == mfree[a]
    check('kappa_a = sum_pi (-1)^(r-1)(r-1)! prod m_(block)', ok_ls)
    check('m_a = sum_pi prod kappa_(block)', ok_inv)

# ---------------------------------------------------------------------------
# 4. marking: the per-class counts of a multinomially marked count
# ---------------------------------------------------------------------------
print('\n== multinomial marking, from the composition rule of the pgf ==')
NF, DM = 5, (2, 2)
d = len(DM)
Rp = PolynomialRing(QQ, ['f%d' % k for k in range(1, NF + 1)] + ['p%d' % j for j in range(1, d)])
fs = [Rp(1)] + list(Rp.gens()[:NF])
ps = list(Rp.gens()[NF:])
ps = ps + [1 - sum(ps)]                      # the marking probabilities sum to 1
Pu = PolynomialRing(Rp, ['u%d' % j for j in range(d)])
us = Pu.gens()


def trunc_u(poly):
    out = Pu(0)
    for mon, co in poly.dict().items():
        if all(mon[j] <= DM[j] for j in range(d)):
            out += co * Pu.prod([us[j] ** mon[j] for j in range(d)])
    return out


# with z_c = 1+u_c and sum_c p_c = 1, the aggregate argument is w = 1 + sum p_c u_c,
# so the joint pgf is G(w) = sum_k f_k (sum_c p_c u_c)^k / k!
w = sum(ps[j] * us[j] for j in range(d))
H, pw = Pu(0), Pu(1)
for k in range(NF + 1):
    H = trunc_u(H + fs[k] / factorial(k) * pw)
    pw = trunc_u(pw * w)
F_pgf = {a: prod([factorial(aj) for aj in a])
         * H.monomial_coefficient(Pu.prod([us[j] ** a[j] for j in range(d)]))
         for a in boxes([n + 1 for n in DM])}
F_impl = moment_joint_marking(fs, ps, DM)
check('E[prod (N_c)_(a_c)] = (prod p_c^a_c) f_|a|',
      all(Rp(F_pgf[a] - F_impl[a]) == 0 for a in F_pgf))

# ---------------------------------------------------------------------------
# 5. aggregation: Vandermonde convolution of falling factorials
# ---------------------------------------------------------------------------
print('\n== aggregation of the parts, on product point masses ==')
DA = (5, 5, 5)
d = len(DA)
Ra = PolynomialRing(QQ, ['y%d' % j for j in range(d)])
ys = Ra.gens()
Fpt = {a: Ra.prod([Ra.prod([ys[j] - t for t in range(a[j])]) for j in range(d)])
       for a in boxes(DA)}
tot = sum(ys)
fall_tot = [Ra.prod([tot - t for t in range(n)]) for n in range(min(DA))]
check('f_n = sum_(|a|=n) (n!/prod a_j!) F_a = (X_1+...+X_d)_n',
      moment_joint_aggregate(Fpt, DA) == fall_tot)
check('aggregate o marking = the aggregate factorial moments',
      moment_joint_aggregate(moment_joint_marking(fs, ps, DM),
                             (DM[0] + 1, DM[1] + 1)) == [Rp(v) for v in fs[:DM[0] + 1]])

print('\nmultivariate house of moments verified:', all(results))
