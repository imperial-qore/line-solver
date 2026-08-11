"""Regression tests for the api/moment moment-conversion package.

The package implements the conversions of the "house of moments" of

  A. Heindl and A. van de Liefvoort. Moment conversions for discrete
  distributions. PMCCS, 2003.

between six families of moments of a discrete random variable N: the power
(raw) moments m_n = E[N^n], the central moments m_n^c = E[(N-m_1)^n], the
factorial moments f_n = E[N(N-1)...(N-n+1)], the binomial moments
b_n = E[C(N,n)], and the two families introduced by the reference, the
upward-factorial moments f_n^+ = E[N(N+1)...(N+n-1)] and the negative-binomial
moments b_n^- = E[C(N+n-1,n)].

All vectors use the harmonized 0-based order convention of the reference: a
vector of length n+1 holds the moments of order 0..n, so element 0 is the
order-0 moment, which equals 1 for every family.

The conversions are refereed by two independent oracles, neither of which knows
anything about Stirling or Lah numbers: a brute-force oracle evaluating every
family by direct summation over the pmf, and closed forms that hold for
specific distributions (f_n = lambda^n for Poisson, f_n = (K!/(K-n)!)*p^n for
Binomial(K,p)).

Values agree with the MATLAB and JAR codebases; see test_moment.m in
line-test.git and jline.api.moment.MomentTest.
"""

from math import comb, factorial, exp

import numpy as np
import pytest

from line_solver.api.moment import (
    moment_stirling1,
    moment_stirling2,
    moment_stirlingcycle,
    moment_lah,
    moment_binotrans,
    moment_binotransinv,
    moment_factorial_from_raw,
    moment_raw_from_factorial,
    moment_upfactorial_from_raw,
    moment_raw_from_upfactorial,
    moment_binomial_from_factorial,
    moment_factorial_from_binomial,
    moment_negbinomial_from_upfactorial,
    moment_upfactorial_from_negbinomial,
    moment_binomial_from_negbinomial,
    moment_negbinomial_from_binomial,
    moment_factorial_from_upfactorial,
    moment_upfactorial_from_factorial,
    moment_central_from_raw,
    moment_raw_from_central,
    moment_cumulant_from_raw,
    moment_raw_from_cumulant,
    moment_factcumulant_from_factorial,
    moment_factorial_from_factcumulant,
    moment_joint_factorial_from_raw,
    moment_joint_raw_from_factorial,
    moment_joint_upfactorial_from_raw,
    moment_joint_raw_from_upfactorial,
    moment_joint_binomial_from_factorial,
    moment_joint_factorial_from_binomial,
    moment_joint_negbinomial_from_upfactorial,
    moment_joint_upfactorial_from_negbinomial,
    moment_joint_factorial_from_upfactorial,
    moment_joint_upfactorial_from_factorial,
    moment_joint_negbinomial_from_binomial,
    moment_joint_binomial_from_negbinomial,
    moment_joint_central_from_raw,
    moment_joint_central_from_raw_mean,
    moment_joint_raw_from_central,
    moment_joint_cumulant_from_raw,
    moment_joint_raw_from_cumulant,
    moment_joint_factcumulant_from_factorial,
    moment_joint_factorial_from_factcumulant,
    moment_joint_marking,
    moment_joint_aggregate,
    moment_housematrix,
    moment_tensortrans,
)

TOL = 1e-10


def relerr(a, b):
    """Normwise relative error max|a-b| / max(1,max|b|).

    A componentwise relative error is not a usable criterion for these
    conversions. The triangles carry alternating signs and large binomial
    coefficients, so an entry that is exactly zero in exact arithmetic is
    reached by cancellation between terms of the magnitude of the input: for
    Uniform(0..5) the falling factorial annihilates every atom at order 6, so
    f_6 = 0 exactly, yet it is computed by cancelling terms of size m_6 = 3419
    and lands on ~5e-13 of roundoff. Dividing that by |b_i| = 0 is meaningless.
    Measuring the residual against the norm of the vector is the standard
    conditioning-aware criterion for a linear transform and still resolves a
    genuine defect to ~1e-14 here, four orders inside the TOL gate.
    """
    a = np.asarray(a, dtype=float).ravel()
    b = np.asarray(b, dtype=float).ravel()
    return np.max(np.abs(a - b)) / max(1.0, np.max(np.abs(b)))


def pmf_oracle(kk, pk, N):
    """Every moment family by direct summation over the pmf, for orders 0..N.

    Deliberately naive: it mirrors the definitions of Section 2 of the
    reference and shares no code with api/moment.
    """
    kk = np.asarray(kk, dtype=float).ravel()
    pk = np.asarray(pk, dtype=float).ravel()
    m = np.zeros(N + 1)
    f = np.zeros(N + 1)
    fp = np.zeros(N + 1)
    b = np.zeros(N + 1)
    bm = np.zeros(N + 1)
    mc = np.zeros(N + 1)
    for n in range(N + 1):
        for k, p in zip(kk, pk):
            m[n] += k ** n * p
            pf = 1.0
            pfp = 1.0
            for j in range(n):
                pf *= k - j        # falling factorial k(k-1)...(k-n+1)
                pfp *= k + j       # rising  factorial k(k+1)...(k+n-1)
            f[n] += pf * p
            fp[n] += pfp * p
            ki = int(round(k))
            if ki >= n:
                b[n] += comb(ki, n) * p
            if n == 0:
                bm[n] += p
            else:
                # the k=0 term contributes C(n-1,n) = 0, cf. eq. (7)
                bm[n] += comb(ki + n - 1, n) * p
    m1 = m[1]
    for n in range(N + 1):
        for k, p in zip(kk, pk):
            mc[n] += (k - m1) ** n * p
    return dict(m=m, f=f, fp=fp, b=b, bm=bm, mc=mc, m1=m1)


def pmf_cases():
    """Four pmfs with distinct structure.

    All have mass at k=0, which exercises the k=0 term of the negative-binomial
    moments, and all have finite support so the oracle is exact (the Poisson
    case is truncated far beyond its mean).
    """
    cases = []
    kk = np.arange(11)
    cases.append(("Binomial(10,0.3)", kk,
                  np.array([comb(10, int(k)) * 0.3 ** k * 0.7 ** (10 - k) for k in kk])))
    kk = np.arange(6)
    cases.append(("DiscreteUniform(0..5)", kk, np.ones(6) / 6))
    cases.append(("Irregular", np.arange(6),
                  np.array([0.10, 0.20, 0.05, 0.30, 0.15, 0.20])))
    kk = np.arange(61)
    cases.append(("Poisson(2) truncated", kk,
                  np.array([exp(-2) * 2.0 ** int(k) / factorial(int(k)) for k in kk])))
    return cases


CASES = pmf_cases()
CASE_IDS = [c[0] for c in CASES]


# ---------- triangles ----------------------------------------------------

def test_stirling1_triangle():
    """Signed Stirling numbers of the first kind, textbook rows 0..4."""
    s = moment_stirling1(4)
    np.testing.assert_allclose(s[0], [1, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(s[1], [0, 1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(s[2], [0, -1, 1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(s[3], [0, 2, -3, 1, 0], atol=1e-12)
    np.testing.assert_allclose(s[4], [0, -6, 11, -6, 1], atol=1e-12)
    assert np.count_nonzero(np.triu(s, 1)) == 0


def test_stirling2_triangle():
    """Stirling numbers of the second kind, textbook rows 0..4."""
    S = moment_stirling2(4)
    np.testing.assert_allclose(S[0], [1, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(S[1], [0, 1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(S[2], [0, 1, 1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(S[3], [0, 1, 3, 1, 0], atol=1e-12)
    np.testing.assert_allclose(S[4], [0, 1, 7, 6, 1], atol=1e-12)
    assert np.count_nonzero(np.triu(S, 1)) == 0


def test_stirlingcycle_triangle():
    """Stirling cycle numbers (unsigned Stirling numbers of the first kind)."""
    sigma = moment_stirlingcycle(4)
    np.testing.assert_allclose(sigma[0], [1, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(sigma[1], [0, 1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(sigma[2], [0, 1, 1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(sigma[3], [0, 2, 3, 1, 0], atol=1e-12)
    np.testing.assert_allclose(sigma[4], [0, 6, 11, 6, 1], atol=1e-12)
    assert np.count_nonzero(np.triu(sigma, 1)) == 0
    # row sums equal n!, since the permutations of n elements are partitioned
    # by their number of cycles
    for n in range(5):
        assert abs(sigma[n].sum() - factorial(n)) < 1e-12


def test_lah_triangle():
    """Lah numbers, textbook rows 0..4."""
    L = moment_lah(4)
    np.testing.assert_allclose(L[0], [1, 0, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(L[1], [0, 1, 0, 0, 0], atol=1e-12)
    np.testing.assert_allclose(L[2], [0, 2, 1, 0, 0], atol=1e-12)
    np.testing.assert_allclose(L[3], [0, 6, 6, 1, 0], atol=1e-12)
    np.testing.assert_allclose(L[4], [0, 24, 36, 12, 1], atol=1e-12)
    assert np.count_nonzero(np.triu(L, 1)) == 0


def test_lah_matches_explicit_factorial_form():
    """The recursion must agree with L(n,k) = (n!/k!)*C(n-1,k-1)."""
    n = 8
    L = moment_lah(n)
    for i in range(1, n + 1):
        for j in range(1, i + 1):
            expected = (factorial(i) / factorial(j)) * comb(i - 1, j - 1)
            assert abs(L[i, j] - expected) <= 1e-12 * max(1.0, abs(expected))


def test_stirling_sign_relation():
    """s(n,k) = (-1)^(n-k) * sigma(n,k), eq. (12)."""
    n = 7
    s = moment_stirling1(n)
    sigma = moment_stirlingcycle(n)
    for i in range(n + 1):
        for j in range(i + 1):
            assert abs(s[i, j] - (-1) ** (i - j) * sigma[i, j]) <= 1e-12 * max(1.0, abs(s[i, j]))


def test_stirling_triangles_are_inverse():
    """Eqs. (10) and (11) make the two Stirling triangles mutually inverse."""
    n = 7
    s = moment_stirling1(n)
    S = moment_stirling2(n)
    assert relerr(S @ s, np.eye(n + 1)) < TOL
    assert relerr(s @ S, np.eye(n + 1)) < TOL


@pytest.mark.parametrize("x", [-3, -0.5, 0, 1, 2.5, 4, 9])
def test_stirling1_generates_falling_factorial(x):
    """Defining identity (10): sum_k s(n,k) x^k = x(x-1)...(x-n+1)."""
    n = 6
    s = moment_stirling1(n)
    for i in range(n + 1):
        lhs = sum(s[i, j] * x ** j for j in range(i + 1))
        rhs = 1.0
        for j in range(i):
            rhs *= x - j
        assert abs(lhs - rhs) <= 1e-8 + 1e-10 * abs(rhs)


@pytest.mark.parametrize("x", [-3, -0.5, 0, 1, 2.5, 4, 9])
def test_stirling2_expands_power_into_falling_factorials(x):
    """Defining identity (11): x^n = sum_k S(n,k) x(x-1)...(x-k+1)."""
    n = 6
    S = moment_stirling2(n)
    for i in range(n + 1):
        rhs = 0.0
        for j in range(i + 1):
            ff = 1.0
            for l in range(j):
                ff *= x - l
            rhs += S[i, j] * ff
        assert abs(x ** i - rhs) <= 1e-8 + 1e-10 * abs(x ** i)


def test_triangle_input_validation():
    for fn in (moment_stirling1, moment_stirling2, moment_stirlingcycle, moment_lah):
        with pytest.raises(ValueError):
            fn(-1)
        with pytest.raises(ValueError):
            fn(2.5)


# ---------- binomial transform -------------------------------------------

def test_binotrans_known_values():
    """Eq. (8) evaluated by hand for x = [1,2,5,15].

    y_0 =  1
    y_1 = -1 + 2            =  1
    y_2 =  1 - 4 + 5        =  2
    y_3 = -1 + 6 - 15 + 15  =  5
    """
    y = moment_binotrans([1, 2, 5, 15])
    np.testing.assert_allclose(y, [1, 1, 2, 5], rtol=1e-12, atol=1e-12)


def test_binotrans_and_its_inverse_undo_one_another():
    """Eq. (8) and its inverse (9) undo one another in both orders."""
    x = np.array([1, 2.5, 7, -3, 11, 0.5])
    assert relerr(moment_binotransinv(moment_binotrans(x)), x) < TOL
    assert relerr(moment_binotrans(moment_binotransinv(x)), x) < TOL


# ---------- conversions against the pmf oracle ---------------------------

@pytest.mark.parametrize("name,kk,pk", CASES, ids=CASE_IDS)
def test_all_conversions_vs_pmf_oracle(name, kk, pk):
    """The acceptance gate: all 14 conversions refereed by direct summation."""
    o = pmf_oracle(kk, pk, 6)
    # power <-> factorial, eq. (13)
    assert relerr(moment_factorial_from_raw(o['m']), o['f']) < TOL
    assert relerr(moment_raw_from_factorial(o['f']), o['m']) < TOL
    # power <-> upward-factorial, via the Stirling cycle numbers
    assert relerr(moment_upfactorial_from_raw(o['m']), o['fp']) < TOL
    assert relerr(moment_raw_from_upfactorial(o['fp']), o['m']) < TOL
    # factorial <-> binomial
    assert relerr(moment_binomial_from_factorial(o['f']), o['b']) < TOL
    assert relerr(moment_factorial_from_binomial(o['b']), o['f']) < TOL
    # upward-factorial <-> negative-binomial, eq. (7)
    assert relerr(moment_negbinomial_from_upfactorial(o['fp']), o['bm']) < TOL
    assert relerr(moment_upfactorial_from_negbinomial(o['bm']), o['fp']) < TOL
    # binomial <-> negative-binomial, the shifted binomial transform of eq. (14)
    assert relerr(moment_binomial_from_negbinomial(o['bm']), o['b']) < TOL
    assert relerr(moment_negbinomial_from_binomial(o['b']), o['bm']) < TOL
    # factorial <-> upward-factorial, via the Lah numbers
    assert relerr(moment_factorial_from_upfactorial(o['fp']), o['f']) < TOL
    assert relerr(moment_upfactorial_from_factorial(o['f']), o['fp']) < TOL
    # power <-> central
    assert relerr(moment_central_from_raw(o['m']), o['mc']) < TOL
    assert relerr(moment_raw_from_central(o['mc'], o['m1']), o['m']) < TOL


@pytest.mark.parametrize("name,kk,pk", CASES, ids=CASE_IDS)
def test_order_zero_and_one_invariants(name, kk, pk):
    """Every family agrees at orders 0 and 1; m_1^c = 0 (cf. Section 2)."""
    o = pmf_oracle(kk, pk, 5)
    f = moment_factorial_from_raw(o['m'])
    fp = moment_upfactorial_from_raw(o['m'])
    b = moment_binomial_from_factorial(f)
    bm = moment_negbinomial_from_upfactorial(fp)
    mc = moment_central_from_raw(o['m'])
    for v in (o['m'], f, fp, b, bm, mc):
        assert abs(v[0] - 1.0) < 1e-12
    m1 = o['m1']
    for v in (f, fp, b, bm):
        assert abs(v[1] - m1) <= 1e-10 * max(1.0, abs(m1))
    assert abs(mc[1]) < 1e-10


# ---------- conversions against closed forms -----------------------------

@pytest.mark.parametrize("lam", [0.5, 2.0, 7.0])
def test_poisson_factorial_moments_closed_form(lam):
    """Poisson(lambda) has f_n = lambda^n exactly."""
    N = 6
    kk = np.arange(121)
    pk = np.array([exp(-lam) * lam ** int(k) / factorial(int(k)) for k in kk])
    m = pmf_oracle(kk, pk, N)['m']
    f = moment_factorial_from_raw(m)
    assert relerr(f, lam ** np.arange(N + 1)) < 1e-8


def test_binomial_factorial_moments_closed_form():
    """Binomial(K,p) has f_n = (K!/(K-n)!)*p^n, and f_n = 0 for n > K."""
    N, K, p = 6, 5, 0.4
    kk = np.arange(K + 1)
    pk = np.array([comb(K, int(k)) * p ** k * (1 - p) ** (K - k) for k in kk])
    m = pmf_oracle(kk, pk, N)['m']
    f = moment_factorial_from_raw(m)
    expected = np.zeros(N + 1)
    for n in range(N + 1):
        # the falling factorial annihilates orders above K
        expected[n] = (factorial(K) / factorial(K - n)) * p ** n if n <= K else 0.0
    assert relerr(f, expected) < 1e-8
    assert abs(f[K + 1]) < 1e-8


def test_bernoulli_all_families_by_hand():
    """Bernoulli(p): m_n = p, f_n = 0 (n>=2), f_n^+ = n!*p, b_n^- = p."""
    N, p = 5, 0.3
    o = pmf_oracle([0, 1], [1 - p, p], N)
    f = moment_factorial_from_raw(o['m'])
    fp = moment_upfactorial_from_raw(o['m'])
    b = moment_binomial_from_factorial(f)
    bm = moment_negbinomial_from_upfactorial(fp)
    np.testing.assert_allclose(o['m'], [1] + [p] * N, rtol=1e-12)
    np.testing.assert_allclose(f, [1, p] + [0] * (N - 1), atol=1e-10)
    np.testing.assert_allclose(b, [1, p] + [0] * (N - 1), atol=1e-10)
    np.testing.assert_allclose(fp, [1] + [factorial(n) * p for n in range(1, N + 1)], rtol=1e-10)
    np.testing.assert_allclose(bm, [1] + [p] * N, rtol=1e-10)


def test_deterministic_point_mass():
    """A point mass at c reduces the conversions to the defining identities."""
    N, c = 6, 4
    o = pmf_oracle([c], [1.0], N)
    f = moment_factorial_from_raw(o['m'])
    fp = moment_upfactorial_from_raw(o['m'])
    for n in range(N + 1):
        ff, rf = 1.0, 1.0
        for j in range(n):
            ff *= c - j
            rf *= c + j
        assert abs(f[n] - ff) <= 1e-8 + 1e-10 * abs(ff)
        assert abs(fp[n] - rf) <= 1e-8 + 1e-10 * abs(rf)


# ---------- structural properties ----------------------------------------

@pytest.mark.parametrize("name,kk,pk", CASES, ids=CASE_IDS)
def test_house_of_moments_commutes(name, kk, pk):
    """Figure 1 is a commuting diagram: two routes to a family must agree.

    This catches a wrong triangle on one edge that a pure round-trip test
    (which traverses the same edge forwards and backwards) cannot see.
    """
    o = pmf_oracle(kk, pk, 6)
    f = moment_factorial_from_raw(o['m'])
    fp = moment_upfactorial_from_raw(o['m'])
    # raw -> upfactorial directly, vs raw -> factorial -> upfactorial
    assert relerr(moment_upfactorial_from_factorial(f), fp) < TOL
    # raw -> factorial directly, vs raw -> upfactorial -> factorial
    assert relerr(moment_factorial_from_upfactorial(fp), f) < TOL
    # binomial by the n! route vs by the shifted-binomial-transform route
    bm = moment_negbinomial_from_upfactorial(fp)
    assert relerr(moment_binomial_from_negbinomial(bm), moment_binomial_from_factorial(f)) < TOL
    # negative-binomial by the n! route vs by the shifted-transform route
    b = moment_binomial_from_factorial(f)
    assert relerr(moment_negbinomial_from_binomial(b), bm) < TOL


@pytest.mark.parametrize("name,kk,pk", CASES, ids=CASE_IDS)
def test_roundtrips(name, kk, pk):
    """Each conversion pair composes to the identity."""
    o = pmf_oracle(kk, pk, 6)
    m = o['m']
    assert relerr(moment_raw_from_factorial(moment_factorial_from_raw(m)), m) < TOL
    assert relerr(moment_raw_from_upfactorial(moment_upfactorial_from_raw(m)), m) < TOL
    assert relerr(moment_raw_from_central(moment_central_from_raw(m), m[1]), m) < TOL
    f = moment_factorial_from_raw(m)
    assert relerr(moment_factorial_from_binomial(moment_binomial_from_factorial(f)), f) < TOL
    assert relerr(moment_factorial_from_upfactorial(moment_upfactorial_from_factorial(f)), f) < TOL
    fp = moment_upfactorial_from_raw(m)
    assert relerr(moment_upfactorial_from_negbinomial(moment_negbinomial_from_upfactorial(fp)), fp) < TOL
    b = moment_binomial_from_factorial(f)
    assert relerr(moment_binomial_from_negbinomial(moment_negbinomial_from_binomial(b)), b) < TOL


def test_central_conversion_holds_for_continuous_rv():
    """Section 5: the power <-> central rules also hold for continuous r.v.s.

    Refereed on Gamma(k,theta), whose raw moments are
    m_n = theta^n * gamma(k+n)/gamma(k) and whose variance is k*theta^2.
    """
    N, k, theta = 5, 3, 2
    m = np.array([theta ** n * factorial(k + n - 1) / factorial(k - 1) for n in range(N + 1)])
    mc = moment_central_from_raw(m)
    assert abs(mc[0] - 1.0) < 1e-12
    assert abs(mc[1]) < 1e-10
    assert abs(mc[2] - k * theta ** 2) <= 1e-10 * k * theta ** 2      # variance
    assert abs(mc[3] - 2 * k * theta ** 3) <= 1e-10 * 2 * k * theta ** 3
    assert relerr(moment_raw_from_central(mc, m[1]), m) < TOL


# ---------- interface behaviour ------------------------------------------

def test_scalar_input_is_order_zero_only():
    """A length-1 input carries only the order-0 moment and passes through."""
    for fn in (moment_factorial_from_raw, moment_upfactorial_from_raw,
               moment_binomial_from_factorial, moment_binomial_from_negbinomial,
               moment_factorial_from_upfactorial):
        out = fn([1.0])
        assert abs(out[0] - 1.0) < 1e-12


def test_central_from_raw_requires_the_mean():
    """m_n^c is defined relative to m_1, so a length-1 input is rejected."""
    with pytest.raises(ValueError):
        moment_central_from_raw([1.0])


def test_list_and_array_inputs_agree():
    """A plain list and a numpy array must give the same result."""
    ml = [1, 2, 6, 22]
    ma = np.array(ml, dtype=float)
    for fn in (moment_factorial_from_raw, moment_upfactorial_from_raw,
               moment_central_from_raw, moment_binotrans, moment_binotransinv):
        np.testing.assert_allclose(fn(ml), fn(ma), rtol=1e-12, atol=1e-12)


# ---------------------------------------------------------------------------
# Cumulants and factorial cumulants
# ---------------------------------------------------------------------------

def _joint_oracle(law, dims):
    """Brute-force joint moment families of a discrete law given as a list of
    (probability, value tuple) pairs, by direct summation."""
    def ex(g):
        return sum(p * g(v) for p, v in law)

    def fall(x, i):
        out = 1.0
        for t in range(i):
            out *= x - t
        return out

    def rise(x, i):
        out = 1.0
        for t in range(i):
            out *= x + t
        return out

    m = np.zeros(dims)
    f = np.zeros(dims)
    fp = np.zeros(dims)
    for a in np.ndindex(*dims):
        m[a] = ex(lambda v, a=a: np.prod([v[j] ** a[j] for j in range(len(a))]))
        f[a] = ex(lambda v, a=a: np.prod([fall(v[j], a[j]) for j in range(len(a))]))
        fp[a] = ex(lambda v, a=a: np.prod([rise(v[j], a[j]) for j in range(len(a))]))
    return m, f, fp


LAW3 = [(0.4, (1, 2)), (0.3, (3, 0)), (0.2, (2, 5)), (0.1, (0, 1))]


def test_cumulants_of_the_poisson_are_all_lambda():
    """Every cumulant of a Poisson variable equals its rate."""
    lam = 2.0
    m = moment_raw_from_cumulant([0.0] + [lam] * 5)
    np.testing.assert_allclose(m, [1, 2, 6, 22, 94, 454], rtol=1e-12)
    np.testing.assert_allclose(moment_cumulant_from_raw(m), [0.0] + [lam] * 5,
                               rtol=1e-12, atol=1e-12)


def test_low_order_cumulants_are_the_textbook_ones():
    """kappa_1, kappa_2 and kappa_3 are the mean, the variance and the third
    central moment."""
    m = [1.0, 2.0, 6.0, 22.0, 94.0]
    k = moment_cumulant_from_raw(m)
    mc = moment_central_from_raw(m)
    assert abs(k[0]) < TOL
    assert abs(k[1] - m[1]) < TOL
    assert abs(k[2] - (m[2] - m[1] ** 2)) < TOL
    assert abs(k[2] - mc[2]) < TOL
    assert abs(k[3] - mc[3]) < TOL
    # the fourth cumulant is the excess, not the fourth central moment
    assert abs(k[4] - (mc[4] - 3 * mc[2] ** 2)) < TOL


def test_factorial_cumulants_of_the_poisson_vanish_beyond_the_first():
    """f_n = lambda^n for a Poisson, so all its factorial cumulants but the
    first are zero. This is what makes them a measure of non-Poissonness."""
    lam = 1.7
    f = [lam ** n for n in range(6)]
    kf = moment_factcumulant_from_factorial(f)
    assert abs(kf[1] - lam) < TOL
    assert max(abs(kf[n]) for n in range(2, 6)) < TOL
    np.testing.assert_allclose(moment_factorial_from_factcumulant(kf), f,
                               rtol=1e-10, atol=1e-12)


def test_cumulant_roundtrips():
    """Each cumulant pair composes to the identity."""
    m = [1.0, 1.5, 4.0, 13.0, 55.0, 260.0]
    np.testing.assert_allclose(moment_raw_from_cumulant(moment_cumulant_from_raw(m)), m,
                               rtol=1e-10)
    f = moment_factorial_from_raw(m)
    np.testing.assert_allclose(
        moment_factorial_from_factcumulant(moment_factcumulant_from_factorial(f)), f,
        rtol=1e-10)


def test_cumulants_preserve_orientation():
    """A row input gives a row output, as everywhere else in the package."""
    m = np.array([[1.0, 2.0, 6.0, 22.0]])
    assert moment_cumulant_from_raw(m).shape == (4,)


# ---------------------------------------------------------------------------
# Joint (multivariate) conversions
# ---------------------------------------------------------------------------

def test_joint_conversions_against_the_brute_force_oracle():
    """Every separable joint conversion reproduces the definition of the target
    family, evaluated by direct summation over the pmf."""
    dims = (5, 5)
    m, f, fp = _joint_oracle(LAW3, dims)
    fact = np.array([factorial(i) for i in range(dims[0])])
    b = f / np.outer(fact, fact)
    bn = fp / np.outer(fact, fact)
    np.testing.assert_allclose(moment_joint_factorial_from_raw(m), f, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_raw_from_factorial(f), m, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_upfactorial_from_raw(m), fp, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_raw_from_upfactorial(fp), m, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_binomial_from_factorial(f), b, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_factorial_from_binomial(b), f, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_negbinomial_from_upfactorial(fp), bn,
                               rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_upfactorial_from_negbinomial(bn), fp,
                               rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_factorial_from_upfactorial(fp), f,
                               rtol=1e-9, atol=1e-8)
    np.testing.assert_allclose(moment_joint_upfactorial_from_factorial(f), fp,
                               rtol=1e-9, atol=1e-8)
    np.testing.assert_allclose(moment_joint_negbinomial_from_binomial(b), bn,
                               rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_binomial_from_negbinomial(bn), b,
                               rtol=1e-9, atol=1e-9)


def test_joint_conversions_reduce_to_the_univariate_ones():
    """A one-dimensional joint array must give exactly what the univariate
    conversion gives, otherwise the two APIs would disagree at d = 1."""
    m = np.array([1.0, 2.0, 6.0, 22.0, 94.0])
    np.testing.assert_allclose(moment_joint_factorial_from_raw(m),
                               moment_factorial_from_raw(m), rtol=1e-12)
    np.testing.assert_allclose(moment_joint_upfactorial_from_raw(m),
                               moment_upfactorial_from_raw(m), rtol=1e-12)
    f = moment_factorial_from_raw(m)
    np.testing.assert_allclose(moment_joint_upfactorial_from_factorial(f),
                               moment_upfactorial_from_factorial(f), rtol=1e-12)
    np.testing.assert_allclose(moment_joint_central_from_raw(m),
                               moment_central_from_raw(m), rtol=1e-12, atol=1e-12)


def test_joint_central_moments_and_the_covariance():
    """The (1,1) entry of the joint central array is the covariance, and the
    conversion inverts."""
    dims = (4, 4)
    m, _, _ = _joint_oracle(LAW3, dims)
    mc = moment_joint_central_from_raw(m)
    assert abs(mc[1, 1] - (m[1, 1] - m[1, 0] * m[0, 1])) < TOL
    assert abs(mc[1, 0]) < TOL and abs(mc[0, 1]) < TOL
    assert abs(mc[2, 0] - (m[2, 0] - m[1, 0] ** 2)) < TOL
    np.testing.assert_allclose(moment_joint_raw_from_central(mc, [m[1, 0], m[0, 1]]), m,
                               rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_central_from_raw_mean(m, [m[1, 0], m[0, 1]]), mc,
                               rtol=1e-12, atol=1e-12)


def test_joint_cumulants_vanish_off_the_axes_under_independence():
    """For independent components the joint moment array is an outer product
    and every mixed cumulant is zero. This is the defining property of the
    joint cumulants and no separable transform has it."""
    dims = (4, 4)
    m1 = moment_raw_from_cumulant([0.0, 2.0, 2.0, 2.0])       # Poisson(2)
    m2 = moment_raw_from_cumulant([0.0, 1.0, 3.0, 5.0])
    m = np.outer(m1, m2)
    k = moment_joint_cumulant_from_raw(m)
    for i in range(1, dims[0]):
        for j in range(1, dims[1]):
            assert abs(k[i, j]) < 1e-8, 'mixed cumulant (%d,%d)' % (i, j)
    np.testing.assert_allclose(k[:, 0], moment_cumulant_from_raw(m1), rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(k[0, :], moment_cumulant_from_raw(m2), rtol=1e-9, atol=1e-9)


def test_joint_cumulant_roundtrip_and_covariance():
    dims = (4, 4)
    m, f, _ = _joint_oracle(LAW3, dims)
    k = moment_joint_cumulant_from_raw(m)
    assert abs(k[1, 1] - (m[1, 1] - m[1, 0] * m[0, 1])) < TOL
    np.testing.assert_allclose(moment_joint_raw_from_cumulant(k), m, rtol=1e-9, atol=1e-9)
    kf = moment_joint_factcumulant_from_factorial(f)
    np.testing.assert_allclose(moment_joint_factorial_from_factcumulant(kf), f,
                               rtol=1e-9, atol=1e-9)


def test_joint_conversions_in_three_dimensions():
    """Nothing in the joint machinery is limited to d = 2."""
    dims = (3, 3, 3)
    law = [(0.5, (1, 2, 0)), (0.3, (3, 0, 1)), (0.2, (2, 1, 4))]
    m, f, fp = _joint_oracle(law, dims)
    np.testing.assert_allclose(moment_joint_factorial_from_raw(m), f, rtol=1e-9, atol=1e-9)
    np.testing.assert_allclose(moment_joint_upfactorial_from_raw(m), fp, rtol=1e-9, atol=1e-9)
    k = moment_joint_cumulant_from_raw(m)
    np.testing.assert_allclose(moment_joint_raw_from_cumulant(k), m, rtol=1e-9, atol=1e-9)


def test_marking_of_a_poisson_gives_independent_poissons():
    """A Poisson count marked independently splits into independent Poisson
    streams, so every joint factorial cumulant of order two or more vanishes."""
    lam, p = 3.0, [0.25, 0.75]
    f = [lam ** n for n in range(7)]
    F = moment_joint_marking(f, p, [3, 3])
    kf = moment_joint_factcumulant_from_factorial(F)
    assert abs(kf[1, 0] - lam * p[0]) < TOL
    assert abs(kf[0, 1] - lam * p[1]) < TOL
    for a in np.ndindex(*F.shape):
        if sum(a) >= 2:
            assert abs(kf[a]) < 1e-8, 'factorial cumulant %s' % (a,)


def test_marking_and_aggregation_are_inverse():
    """Summing the marked classes recovers the aggregate factorial moments."""
    f = [1.0, 2.0, 4.5, 11.0, 30.0]
    p = [0.2, 0.3, 0.5]
    F = moment_joint_marking(f, p, [1, 1, 1])
    np.testing.assert_allclose(moment_joint_aggregate(F), f[:2], rtol=1e-12)
    F2 = moment_joint_marking(f, [0.4, 0.6], [2, 2])
    np.testing.assert_allclose(moment_joint_aggregate(F2), f[:3], rtol=1e-12)


def test_aggregation_holds_for_a_dependent_joint_law():
    """The aggregation identity is Vandermonde, so it needs no independence:
    the factorial moments of X+Y follow from the joint ones of (X,Y)."""
    dims = (5, 5)
    _, f, _ = _joint_oracle(LAW3, dims)
    tot = [(0.4, (3,)), (0.3, (3,)), (0.2, (7,)), (0.1, (1,))]
    _, fs, _ = _joint_oracle(tot, (5,))
    np.testing.assert_allclose(moment_joint_aggregate(f), fs, rtol=1e-9, atol=1e-9)


def test_house_matrix_and_tensor_transform():
    """moment_housematrix must return the matrix of the univariate edge, and
    moment_tensortrans must apply it along the requested mode only."""
    T = moment_housematrix('factorial_from_raw', 4)
    np.testing.assert_allclose(T, moment_stirling1(4), rtol=1e-12)
    m = np.arange(20, dtype=float).reshape(4, 5) + 1.0
    got = moment_tensortrans(m, moment_housematrix('raw_from_factorial', 4), 1)
    for i in range(4):
        np.testing.assert_allclose(got[i, :], moment_raw_from_factorial(m[i, :]), rtol=1e-12)
    with pytest.raises(ValueError):
        moment_tensortrans(m, T, 5)
    with pytest.raises(ValueError):
        moment_housematrix('no_such_edge', 3)


def test_joint_error_paths():
    with pytest.raises(ValueError):
        moment_joint_central_from_raw(np.ones((3, 1)))
    with pytest.raises(ValueError):
        moment_joint_raw_from_central(np.ones((3, 3)), [1.0])
    with pytest.raises(ValueError):
        moment_joint_marking([1.0, 2.0], [0.5, 0.5], [2, 2])
    with pytest.raises(ValueError):
        moment_joint_marking([1.0, 2.0], [0.5, 0.5], [1])
