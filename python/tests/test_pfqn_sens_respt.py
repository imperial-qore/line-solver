"""Regression tests for pfqn_sens_respt.

pfqn_sens_respt returns the exact raw moments E[W_(i,l)^t] of the sojourn time of
a class-l job at an FCFS b-server center of a closed product-form network
(Strelen 1990, Theorem 4.1, equations (4.1)-(4.5) and Remarks 4.2-4.3).

These tests mirror the MATLAB validator pfqn_sens_respt_validate.m. Five
references:

  A. brute-force enumeration. This is the strongest check because it shares none
     of Theorem 4.1's algebra. The equilibrium product form is enumerated to get
     the exact arrival-theorem marginals p_i(j, N-1_l); the sojourn time
     conditioned on finding j jobs is known in closed form (Exp(mu) if j < b,
     otherwise an Erlang(j-b+1, b*mu) queueing delay plus an Exp(mu) service), so
     its moments are formed directly and mixed over j. Neither the coefficients
     a_(t,tau)(0) nor the recursion (4.2) enter, so the agreement tests both;
  B. the published table of Example 3.4 (continued) of the reference, which
     prints E(W_i) and sigma^2_(W_i) for the Kobayashi model;
  C. the internal identity W(i,l) = w_i(l)/V(i,l): the t = 1 case of (4.5) must
     reproduce the MVA residence time divided by the visit ratio, which is a
     completely different expression;
  D. pfqn_mva for the base measures in the single-server case;
  E. the single-job network, where an arriving job always finds an empty station,
     so W is exactly Exp(mu) and every moment is known in closed form.
"""

import itertools
from math import comb, factorial

import numpy as np
import pytest

from line_solver.api.pfqn import (
    pfqn_mva,
    pfqn_sens_respt,
)


# =========================================================================
# helpers
# =========================================================================

def _relerr(a, b):
    """Mirrors the MATLAB validators' relerr: absolute error scaled by
    max(1, max(|a|,|b|)), so small entries are held to an absolute bound."""
    a = np.asarray(a, float).flatten()
    b = np.asarray(b, float).flatten()
    if a.size == 0:
        return 0.0
    scale = np.maximum(1.0, np.maximum(np.abs(a), np.abs(b)))
    return float(np.max(np.abs(a - b) / scale))


def _compositions_leq(n, M):
    """All nonnegative integer vectors of length M summing to at most n."""
    if M == 1:
        return [(k,) for k in range(n + 1)]
    out = []
    for first in range(n + 1):
        for sub in _compositions_leq(n - first, M - 1):
            out.append((first,) + sub)
    return out


def _brute_marginals(S, V, N, Z, b):
    """P[Q_i = j] for every station, by enumerating the closed product form of a
    network of FCFS b-server stations:

      f_i(q_i) = q_i! prod_l (a(i,l)^q_il / q_il!) prod_{j=1}^{q_i} 1/min(j,b_i)

    with a(i,l) = S(i)*V(i,l), plus the delay term for the think times.
    """
    from scipy.special import gammaln

    V = np.asarray(V, float)
    S = np.asarray(S, float).flatten()
    N = np.asarray(N, int).flatten()
    Z = np.asarray(Z, float).flatten()
    b = np.asarray(b, int).flatten()
    M, R = V.shape

    a = np.zeros((M, R))
    for i in range(M):
        a[i, :] = S[i] * V[i, :]

    per = [_compositions_leq(int(N[r]), M) for r in range(R)]
    weights = []
    tots = []
    for combo in itertools.product(*per):
        nir = np.array(combo, dtype=int).T          # (M, R)
        lw = 0.0
        ok = True
        for i in range(M):
            ni = int(nir[i, :].sum())
            lw += gammaln(ni + 1)
            for j in range(1, ni + 1):
                lw -= np.log(min(j, b[i]))
            for r in range(R):
                if nir[i, r] > 0:
                    if a[i, r] <= 0:
                        ok = False
                        break
                    lw += nir[i, r] * np.log(a[i, r]) - gammaln(nir[i, r] + 1)
            if not ok:
                break
        if ok:
            for r in range(R):
                n0r = int(N[r] - nir[:, r].sum())
                if n0r > 0:
                    if Z[r] <= 0:
                        ok = False
                        break
                    lw += n0r * np.log(Z[r]) - gammaln(n0r + 1)
        weights.append(np.exp(lw) if ok else 0.0)
        tots.append(nir.sum(axis=1))

    w = np.array(weights, float)
    w /= w.sum()
    tot = np.stack(tots, axis=0)                    # (K, M)

    maxj = max(int(np.sum(N)), 1)
    pj = np.zeros((M, 1 + maxj))
    for k in range(len(w)):
        for i in range(M):
            pj[i, tot[k, i]] += w[k]
    return pj


def _cond_moment(j, b, mu, t):
    """E[(W|j)^t] where a job arriving to find j jobs at an FCFS b-server station
    waits an Erlang(max(0, j-b+1), b*mu) and is then served for an Exp(mu)."""
    k = max(0, j - b + 1)
    v = 0.0
    for s in range(t + 1):
        # E[X^s] with X ~ Exp(mu)
        EX = factorial(s) / mu ** s
        p = t - s
        # E[Y^p] with Y ~ Erlang(k, b*mu), and Y = 0 when k = 0
        if k == 0:
            EY = 1.0 if p == 0 else 0.0
        else:
            theta = b * mu
            EY = 1.0
            for aa in range(p):
                EY *= (k + aa)
            EY /= theta ** p
        v += comb(t, s) * EX * EY
    return v


def _brute_respt(S, V, N, Z, b, tmax):
    """Sojourn-time moments from first principles: enumerate the product form to
    get the exact arrival-theorem marginals p_i(j, N-e_l), then mix the
    conditional sojourn-time moments over j. Uses none of Theorem 4.1."""
    V = np.asarray(V, float)
    S = np.asarray(S, float).flatten()
    N = np.asarray(N, int).flatten()
    b = np.asarray(b, int).flatten()
    M, R = V.shape
    bmax = int(np.max(b))
    WM = np.zeros((M, R, tmax))
    for l in range(R):
        if N[l] == 0:
            continue
        Nl = N.copy()
        Nl[l] -= 1
        pj = _brute_marginals(S, V, Nl, Z, b)      # pj[i, j] at population N - e_l
        for i in range(M):
            if V[i, l] <= 0:
                continue
            mu = 1.0 / S[i]
            for t in range(1, tmax + 1):
                acc = 0.0
                for j in range(pj.shape[1]):
                    acc += pj[i, j] * _cond_moment(j, b[i], mu, t)
                WM[i, l, t - 1] = acc
    pAll = _brute_marginals(S, V, N, Z, b)
    pN = np.zeros((M, bmax))
    for i in range(M):
        for j in range(bmax):
            pN[i, j] = pAll[i, j]
    return WM, pN


def _random_models(seed=5, trials=36):
    """Mirrors the model sweep of pfqn_sens_respt_validate.m."""
    rng = np.random.default_rng(seed)
    out = []
    for trial in range(1, trials + 1):
        M = int(rng.integers(1, 4))
        R = int(rng.integers(1, 3))
        S = 0.2 + rng.random(M)
        V = 0.3 + rng.random((M, R))
        if trial % 4 == 0 and M > 1:
            V[0, 0] = 0.0          # a class that skips a station
        N = rng.integers(1, 4, size=R).astype(float)
        Z = 0.3 + rng.random(R) if trial % 2 == 0 else np.zeros(R)
        b = (rng.integers(1, 4, size=M) if trial % 3 == 0
             else np.ones(M, dtype=int))
        out.append((S, V, N, Z, b))
    return out


MODELS = _random_models()


# =========================================================================
# pfqn_sens_respt
# =========================================================================

def test_sens_respt_matches_brute_force_arrival_theorem():
    """A. Ground truth: the exact arrival-theorem marginals from enumeration,
    mixed with the closed-form conditional sojourn-time moments. Shares none of
    Theorem 4.1's algebra."""
    err = 0.0
    n = 0
    for S, V, N, Z, b in MODELS:
        if np.prod(N + 1) > 24 or V.shape[0] > 3:
            continue
        res = pfqn_sens_respt(S, V, N, Z, b, 3)
        Wb, pb = _brute_respt(S, V, N, Z, b, 3)
        err = max(err, _relerr(res.WM, Wb))
        # .p is ragged: station i only defines j = 0..b_i-1, the range the
        # b-server recursion needs, and the rest of the row is zero padding out
        # to max(b). Comparing the padding against the true marginal would be
        # comparing against something the algorithm never claims to compute.
        for i in range(V.shape[0]):
            err = max(err, _relerr(res.p[i, :b[i]], pb[i, :b[i]]))
        n += 1
    assert n > 0
    assert err <= 1e-9, f"brute-force disagreement {err:.3e} over {n} models"


def test_sens_respt_identity_w_equals_residence_over_visits():
    """C. W(i,l) = w_i(l)/V(i,l): the t = 1 case of (4.5) must reproduce the MVA
    residence time divided by the visit ratio, a completely different
    expression."""
    err = 0.0
    for S, V, N, Z, b in MODELS:
        res = pfqn_sens_respt(S, V, N, Z, b, 3)
        M, R = V.shape
        for i in range(M):
            for l in range(R):
                if V[i, l] > 0 and N[l] > 0:
                    err = max(err, _relerr(res.W[i, l],
                                           res.Wresid[i, l] / V[i, l]))
    assert err <= 1e-10, f"identity W = w/V violated by {err:.3e}"


def test_sens_respt_base_measures_match_pfqn_mva():
    """D. The primal must reproduce pfqn_mva entry by entry in the single-server
    case, otherwise the sojourn moments are of the wrong network."""
    err = 0.0
    n = 0
    for S, V, N, Z, b in MODELS:
        if not np.all(b == 1):
            continue
        M, R = V.shape
        L = np.zeros((M, R))
        for i in range(M):
            L[i, :] = S[i] * V[i, :]
        res = pfqn_sens_respt(S, V, N, Z, b, 3)
        r = pfqn_mva(L, N, Z)
        XN, QN, UN = r[0], r[2], r[3]
        err = max(err, _relerr(res.X, XN), _relerr(res.Q, QN),
                  _relerr(res.U, UN))
        n += 1
    assert n > 0
    assert err <= 1e-10, f"base measures differ from pfqn_mva by {err:.3e}"


def test_sens_respt_single_job_is_exactly_exponential():
    """E. With one job the arriving job always finds the station empty, so
    W ~ Exp(mu) exactly and E[W^t] = t!/mu^t. An analytic check independent of
    every other reference."""
    S1 = np.array([0.4, 0.25])
    V1 = np.array([[1.0], [2.0]])
    r1 = pfqn_sens_respt(S1, V1, np.array([1.0]), np.array([0.7]),
                         np.array([1, 1]), 3)
    err = 0.0
    for i in range(2):
        mu = 1.0 / S1[i]
        err = max(err, _relerr(r1.WM[i, 0, 0], 1.0 / mu),
                  _relerr(r1.WM[i, 0, 1], 2.0 / mu ** 2),
                  _relerr(r1.WM[i, 0, 2], 6.0 / mu ** 3),
                  _relerr(r1.WVar[i, 0], 1.0 / mu ** 2))
    assert err <= 1e-12, f"single-job Exp(mu) disagreement {err:.3e}"


def test_sens_respt_strelen_example_34_published_sojourn():
    """B. Strelen Example 3.4 (continued): the Kobayashi central-server model at
    n = 3. The paper prints E(W_i) and sigma^2_(W_i) to five decimals."""
    xs = np.concatenate([np.full(9, 0.0215), [0.104, 0.104, 0.019]])
    es = np.concatenate([np.full(9, 9.333), [10.5, 10.5, 105.0]])
    rk = pfqn_sens_respt(xs, es.reshape(-1, 1), np.array([3.0]),
                         np.array([0.0]), np.ones(12, dtype=int), 3)
    paperW = np.array([0.02275, 0.14178, 0.03322])
    paperWV = np.array([0.00052, 0.01846, 0.00083])
    gotW = np.array([rk.W[0, 0], rk.W[9, 0], rk.W[11, 0]])
    gotWV = np.array([rk.WVar[0, 0], rk.WVar[9, 0], rk.WVar[11, 0]])
    err = max(_relerr(gotW, paperW), _relerr(gotWV, paperWV))
    assert err <= 5e-4, (f"Strelen Example 3.4 sojourn disagreement {err:.3e} "
                         f"(paper E(W_12)={paperW[2]:.5f} got {gotW[2]:.5f}; "
                         f"sigma2={paperWV[2]:.5f} got {gotWV[2]:.5f})")


def test_sens_respt_multiserver_marginals_sum_below_one():
    """The ragged .p row must be a valid partial marginal: nonnegative and
    summing to at most one over j = 0..b_i-1."""
    for S, V, N, Z, b in MODELS:
        res = pfqn_sens_respt(S, V, N, Z, b, 3)
        for i in range(V.shape[0]):
            row = res.p[i, :b[i]]
            assert np.all(row >= -1e-12), f"negative marginal {row}"
            assert np.sum(row) <= 1.0 + 1e-12, f"marginal mass {np.sum(row)} > 1"


def test_sens_respt_tmax_one_and_two_agree_with_tmax_three():
    """tmax only truncates the reported moments; the ones that are reported must
    not depend on how many were asked for."""
    S = np.array([0.4, 0.7, 0.3])
    V = np.array([[1.0, 0.5], [2.0, 1.5], [0.5, 1.0]])
    N = np.array([2.0, 1.0])
    Z = np.array([0.3, 0.0])
    b = np.array([1, 2, 3])
    r3 = pfqn_sens_respt(S, V, N, Z, b, 3)
    r2 = pfqn_sens_respt(S, V, N, Z, b, 2)
    r1 = pfqn_sens_respt(S, V, N, Z, b, 1)
    np.testing.assert_allclose(r1.WM[:, :, 0], r3.WM[:, :, 0], rtol=1e-14)
    np.testing.assert_allclose(r2.WM[:, :, :2], r3.WM[:, :, :2], rtol=1e-14)
    np.testing.assert_allclose(r2.WVar, r3.WVar, rtol=1e-14)


def test_sens_respt_empty_population_is_zero():
    S = np.array([0.4, 0.7])
    V = np.array([[1.0], [2.0]])
    res = pfqn_sens_respt(S, V, np.array([0.0]), np.array([1.0]))
    np.testing.assert_allclose(res.W, 0.0, atol=0)
    np.testing.assert_allclose(res.WM, 0.0, atol=0)
    np.testing.assert_allclose(res.Q, 0.0, atol=0)
    np.testing.assert_allclose(res.m, 0.0, atol=0)


def test_sens_respt_rejects_open_population():
    with pytest.raises(ValueError):
        pfqn_sens_respt(np.array([1.0]), np.array([[1.0]]), np.array([np.inf]))


def test_sens_respt_rejects_bad_tmax():
    with pytest.raises(ValueError):
        pfqn_sens_respt(np.array([1.0]), np.array([[1.0]]), np.array([2.0]),
                        None, None, 4)


def test_sens_respt_rejects_zero_servers():
    with pytest.raises(ValueError):
        pfqn_sens_respt(np.array([1.0]), np.array([[1.0]]), np.array([2.0]),
                        None, np.array([0]))


def test_sens_respt_rejects_nonpositive_service_time():
    with pytest.raises(ValueError):
        pfqn_sens_respt(np.array([0.0]), np.array([[1.0]]), np.array([2.0]))
