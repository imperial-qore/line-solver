"""Regression tests for pfqn_sens_mva, pfqn_sens_ldmx_ec and pfqn_sens_mvaldmx.

pfqn_sens_mva returns the exact per-station queue-length variances and
covariances of a closed product-form network by the MVA-like moment recursion of
de Souza e Silva and Muntz (1988), Corollary 1. pfqn_sens_mvaldmx is its mixed
load-dependent counterpart, the moment analysis of Akyildiz and Strelen (1991),
which additionally reaches the cross-station covariances.

These tests mirror the MATLAB validators pfqn_sens_mva_validate.m and
pfqn_sens_mvaldmx_validate.m. The references are chosen so that every channel of
the derivation is exercised by something that does not share its code:

  * brute-force enumeration of the product-form equilibrium distribution
    (ground truth);
  * the differentiated-MVA Jacobian of pfqn_sens, via the identity
    Cov[n(i,r),n(i,s)] = L(i,s) dQ(i,r)/dL(i,s). The reference is read off the
    raw Jacobian dQ, NOT off sens.QCov: pfqn_sens now sources its same-station
    blocks from pfqn_sens_mva, so comparing against sens.QCov would compare the
    recursion with itself;
  * pfqn_mva / pfqn_mvaldmx for the base measures. The primal must be
    reproduced entry by entry, otherwise the derivative is of the wrong
    function;
  * central finite differences of pfqn_mvaldmx with respect to the
    demand-scaling parameter y(j,s);
  * the raw asymmetry of the covariance before symmetrization. The recursion
    computes W(k,j;t,j) and W(t,j;k,j) by numerically distinct expressions, so
    their agreement is a nontrivial check of the formula.
"""

import itertools

import numpy as np
import pytest

from line_solver.api.pfqn import (
    pfqn_mva,
    pfqn_mvaldmx,
    pfqn_sens_mva,
    pfqn_sens_mvaldmx,
    pfqn_ldmx_ec,
    pfqn_sens_ldmx_ec,
    pfqn_sens,
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


def _brute_moments(L, N, Z):
    """Exact moments by enumerating the closed product-form distribution.

    Stations 0..M-1 are single-server fixed-rate centers; the think time Z is an
    infinite-server station that carries no moment:
      p(n) ~ prod_i [ n_i! prod_r L(i,r)^n(i,r)/n(i,r)! ]
             * prod_r Z(r)^n(0,r)/n(0,r)!
    """
    from scipy.special import gammaln

    L = np.asarray(L, float)
    N = np.asarray(N, int)
    Z = np.asarray(Z, float)
    M, R = L.shape

    per = [_compositions_leq(int(N[r]), M) for r in range(R)]
    states = []
    weights = []
    for combo in itertools.product(*per):
        # nir[i, r] = jobs of class r at station i
        nir = np.array(combo, dtype=int).T          # (M, R)
        lw = 0.0
        ok = True
        for i in range(M):
            ni = int(nir[i, :].sum())
            lw += gammaln(ni + 1)
            for r in range(R):
                if nir[i, r] > 0:
                    if L[i, r] <= 0:
                        ok = False
                        break
                    lw += nir[i, r] * np.log(L[i, r]) - gammaln(nir[i, r] + 1)
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
        states.append(nir)
        weights.append(np.exp(lw) if ok else 0.0)

    w = np.array(weights, float)
    w /= w.sum()
    S = np.stack(states, axis=0)                    # (K, M, R)

    Q = np.tensordot(w, S, axes=(0, 0))             # (M, R)
    QCov = np.zeros((M, R, R))
    for i in range(M):
        for r in range(R):
            for s in range(R):
                m2 = float(np.sum(w * S[:, i, r] * S[:, i, s]))
                QCov[i, r, s] = m2 - Q[i, r] * Q[i, s]
    return Q, QCov


def _mu_at(mu, i, j):
    """Limited load dependence: the rate saturates at its last tabulated value."""
    return mu[i, j - 1] if j <= mu.shape[1] else mu[i, -1]


def _brute_ldmx(lam, D, N, Z, mu, Kopen):
    """Exact moments by enumerating the mixed load-dependent product form.

      p(n) ~ prod_i [ n_i! prod_r a(i,r)^n(i,r)/n(i,r)!
                      prod_{j=1}^{n_i} 1/mu(i,j) ]
             * prod_{closed c} Z(c)^n(0,c)/n(0,c)!

    with a(i,r) = D(i,r) for a closed class and a(i,r) = lambda(r)*D(i,r) for an
    open class, n_i the total population at station i, and the closed classes
    constrained to sum to N. Open classes are truncated at Kopen jobs.
    """
    from scipy.special import gammaln

    D = np.atleast_2d(np.asarray(D, float))
    N = np.asarray(N, float).flatten()
    Z = np.asarray(Z, float).flatten()
    lam = np.asarray(lam, float).flatten()
    mu = np.atleast_2d(np.asarray(mu, float))
    M, R = D.shape

    a = np.zeros((M, R))
    for r in range(R):
        if np.isinf(N[r]):
            a[:, r] = lam[r] * D[:, r]
        else:
            a[:, r] = D[:, r]

    alloc = []
    for r in range(R):
        if np.isinf(N[r]):
            alloc.append(_compositions_leq(Kopen, M))
        else:
            alloc.append(_compositions_leq(int(N[r]), M))

    states = []
    weights = []
    for combo in itertools.product(*alloc):
        nir = np.array(combo, dtype=int).T          # (M, R)
        lw = 0.0
        ok = True
        for i in range(M):
            ni = int(nir[i, :].sum())
            lw += gammaln(ni + 1)
            for j in range(1, ni + 1):
                lw -= np.log(_mu_at(mu, i, j))
            for r in range(R):
                if nir[i, r] > 0:
                    if a[i, r] <= 0:
                        ok = False
                        break
                    lw += nir[i, r] * np.log(a[i, r]) - gammaln(nir[i, r] + 1)
            if not ok:
                break
        if ok:
            for c in range(R):
                if np.isinf(N[c]):
                    continue
                n0c = int(N[c] - nir[:, c].sum())
                if n0c > 0:
                    if Z[c] <= 0:
                        ok = False
                        break
                    lw += n0c * np.log(Z[c]) - gammaln(n0c + 1)
        states.append(nir)
        weights.append(np.exp(lw) if ok else 0.0)

    w = np.array(weights, float)
    w /= w.sum()
    S = np.stack(states, axis=0)

    Q = np.tensordot(w, S, axes=(0, 0))
    QCov = np.zeros((M, R, R))
    for i in range(M):
        for r in range(R):
            for s in range(R):
                m2 = float(np.sum(w * S[:, i, r] * S[:, i, s]))
                QCov[i, r, s] = m2 - Q[i, r] * Q[i, s]
    return Q, QCov


def _sens_cov_ref(L, N, Z, mi):
    """Same-station covariance read off the raw pfqn_sens Jacobian:
    Cov[n(i,r),n(i,s)] = L(i,s) * dQ(i,r)/dL(i,s)."""
    L = np.asarray(L, float)
    M, R = L.shape
    sens = pfqn_sens(L, N, Z, mi)
    pL = -np.ones((M, R), dtype=int)
    for p, pr in enumerate(sens.params):
        if pr['type'] == 'L':
            pL[pr['station'], pr['jobclass']] = p
    ref = np.zeros((M, R, R))
    for i in range(M):
        for r in range(R):
            for s in range(R):
                if pL[i, s] >= 0:
                    ref[i, r, s] = L[i, s] * np.asarray(sens.dQ)[i, r, pL[i, s]]
    return ref


def _random_closed_models(seed=0, trials=40):
    """Mirrors the model sweep of pfqn_sens_mva_validate.m."""
    rng = np.random.default_rng(seed)
    out = []
    for trial in range(1, trials + 1):
        M = int(rng.integers(1, 4))
        R = int(rng.integers(1, 4))
        L = 0.2 + rng.random((M, R))
        N = rng.integers(0, 4, size=R)
        if not np.any(N > 0):
            N[0] = 2
        Z = 0.3 + rng.random(R) if trial % 2 == 0 else np.zeros(R)
        # exercise a zero-demand column now and then: class 0 never visits
        # station 0
        if trial % 5 == 0 and M > 1:
            L[0, 0] = 0.0
        out.append((L, N.astype(float), Z))
    return out


CLOSED_MODELS = _random_closed_models()


# =========================================================================
# pfqn_sens_mva
# =========================================================================

def test_sens_mva_base_measures_match_pfqn_mva():
    """C. The primal must reproduce pfqn_mva entry by entry."""
    err = 0.0
    for L, N, Z in CLOSED_MODELS:
        mom = pfqn_sens_mva(L, N, Z)
        # pfqn_mva returns (XN, CN, QN, UN, RN, TN, AN); RN is the (M x R)
        # residence time, i.e. the MATLAB pfqn_mva 4th output that mom.R mirrors
        r = pfqn_mva(L, N, Z)
        XN, QN, UN, RN = r[0], r[2], r[3], r[4]
        err = max(err, _relerr(mom.X, XN), _relerr(mom.Q, QN),
                  _relerr(mom.U, UN), _relerr(mom.R, RN))
    assert err <= 1e-10, f"base measures differ from pfqn_mva by {err:.3e}"


def test_sens_mva_matches_brute_force_product_form():
    """A. Ground truth: enumeration of the closed product-form distribution."""
    err = 0.0
    n = 0
    for L, N, Z in CLOSED_MODELS:
        if np.prod(N + 1) > 64 or L.shape[0] > 3:
            continue
        mom = pfqn_sens_mva(L, N, Z)
        Qb, QCovb = _brute_moments(L, N, Z)
        err = max(err, _relerr(mom.Q, Qb), _relerr(mom.QCov, QCovb))
        n += 1
    assert n > 0
    assert err <= 1e-9, f"brute-force disagreement {err:.3e} over {n} models"


def test_sens_mva_matches_pfqn_sens_jacobian():
    """B. The identity Cov[n(i,r),n(i,s)] = L(i,s) dQ(i,r)/dL(i,s), including
    mi > 1, which brute force does not cover."""
    rng = np.random.default_rng(7)
    err = 0.0
    n = 0
    for L, N, Z in CLOSED_MODELS:
        M = L.shape[0]
        for micase in range(2):
            mi = np.ones(M) if micase == 0 else rng.integers(1, 4, size=M).astype(float)
            mom = pfqn_sens_mva(L, N, Z, mi)
            ref = _sens_cov_ref(L, N, Z, mi)
            err = max(err, _relerr(mom.QCov, ref))
            # Theorem 3: variance of the total queue length at a station
            tot_ref = np.array([ref[i, :, :].sum() for i in range(M)])
            err = max(err, _relerr(mom.QTotVar, tot_ref))
            n += 1
    assert n > 0
    assert err <= 1e-9, f"Jacobian disagreement {err:.3e} over {n} models"


def test_sens_mva_raw_covariance_is_symmetric():
    """D. The two triangles are numerically distinct expressions; their
    agreement before symmetrization is an independent residual."""
    rng = np.random.default_rng(7)
    err = 0.0
    for L, N, Z in CLOSED_MODELS:
        M = L.shape[0]
        for micase in range(2):
            mi = np.ones(M) if micase == 0 else rng.integers(1, 4, size=M).astype(float)
            err = max(err, pfqn_sens_mva(L, N, Z, mi).QCovAsym)
    assert err <= 1e-9, f"raw QCov asymmetry {err:.3e}"


def test_sens_mva_qvar_is_bernoulli_at_unit_population():
    """With N=1 the per-station queue length is Bernoulli, so Var = q(1-q).
    An analytic check independent of every other reference."""
    L = np.array([[1.0], [0.5]])
    mom = pfqn_sens_mva(L, np.array([1.0]), np.array([0.0]))
    q = np.asarray(mom.Q).flatten()
    np.testing.assert_allclose(np.asarray(mom.QVar).flatten(), q * (1.0 - q),
                               rtol=1e-12, atol=1e-14)


def test_sens_mva_empty_population_is_zero():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    mom = pfqn_sens_mva(L, np.array([0.0, 0.0]), np.array([1.0, 1.0]))
    np.testing.assert_allclose(mom.Q, 0.0, atol=0)
    np.testing.assert_allclose(mom.QCov, 0.0, atol=0)
    np.testing.assert_allclose(mom.QTotVar, 0.0, atol=0)
    assert mom.QCovAsym == 0.0


def test_sens_mva_rejects_open_population():
    with pytest.raises(ValueError):
        pfqn_sens_mva(np.array([[1.0]]), np.array([np.inf]), np.array([0.0]))


# =========================================================================
# pfqn_sens wiring: same-station blocks must come from the recursion
# =========================================================================

def test_pfqn_sens_sources_same_station_blocks_from_sens_mva():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    sens = pfqn_sens(L, N, Z)
    mom = pfqn_sens_mva(L, N, Z)
    M, R = L.shape
    for i in range(M):
        np.testing.assert_allclose(np.asarray(sens.QCov)[i, :, i, :],
                                   mom.QCov[i, :, :], rtol=0, atol=0)
    np.testing.assert_allclose(np.asarray(sens.QVar), mom.QVar, rtol=0, atol=0)
    np.testing.assert_allclose(np.asarray(sens.QTotVar), mom.QTotVar,
                               rtol=0, atol=0)
    assert sens.QCovAsym == mom.QCovAsym


def test_pfqn_sens_cross_station_blocks_still_from_jacobian():
    """Cross-station blocks are unreachable by the same-station recursion and
    must remain L(j,s) dQ(i,r)/dL(j,s)."""
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([3.0, 4.0]), np.array([1.0, 3.0])
    sens = pfqn_sens(L, N, Z)
    M, R = L.shape
    pL = np.zeros((M, R), dtype=int)
    for p, pr in enumerate(sens.params):
        if pr['type'] == 'L':
            pL[pr['station'], pr['jobclass']] = p
    for i in range(M):
        for j in range(M):
            if i == j:
                continue
            for r in range(R):
                for s in range(R):
                    np.testing.assert_allclose(
                        np.asarray(sens.QCov)[i, r, j, s],
                        L[j, s] * np.asarray(sens.dQ)[i, r, pL[j, s]],
                        rtol=0, atol=0)


# =========================================================================
# pfqn_sens_ldmx_ec
# =========================================================================

def _random_ld_models(seed=1, trials=12):
    """Mirrors the mixed load-dependent sweep of pfqn_sens_mvaldmx_validate.m."""
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(trials):
        M = int(rng.integers(1, 3))
        Ropen = int(rng.integers(0, 2))
        R = 1 + Ropen
        D = 0.2 + 0.6 * rng.random((M, R))
        N = np.zeros(R)
        N[0] = rng.integers(1, 4)                 # closed class
        lam = np.zeros(R)
        if Ropen == 1:
            N[1] = np.inf                         # open class
            lam[1] = 0.05 + 0.15 * rng.random()   # keep station loads modest
        Z = np.zeros(R)
        Z[0] = 0.5 * rng.random()
        NCtot = int(np.sum(N[np.isfinite(N)]))
        # limited load dependence: rates grow up to level b then saturate
        b = int(rng.integers(1, 4))
        mu = np.zeros((M, max(NCtot, 1)))
        for i in range(M):
            for n in range(mu.shape[1]):
                mu[i, n] = min(n + 1, b) * (0.8 + 0.4 * rng.random())
        Lo = np.array([float(np.dot(lam, D[i, :])) for i in range(M)])
        # keep the geometric tail of the limited load dependence stable
        if np.any(Lo / mu[:, -1] > 0.6):
            continue
        out.append((lam, D, N, Z, mu))
    return out


LD_MODELS = _random_ld_models()


def test_sens_ldmx_ec_primal_matches_pfqn_ldmx_ec():
    """The differentiated effective capacities must carry the same primal as
    pfqn_ldmx_ec, otherwise the derivative is of the wrong function."""
    err = 0.0
    for lam, D, N, Z, mu in LD_MODELS:
        mux = np.hstack([mu, mu[:, -1:]])
        EC, E, Eprime, Lo = pfqn_ldmx_ec(lam, D, mux)
        ECs, Es, Eps, Los = pfqn_sens_ldmx_ec(lam, D, mux)[:4]
        err = max(err, _relerr(EC, ECs), _relerr(E, Es),
                  _relerr(Eprime, Eps), _relerr(Lo, Los))
    assert err <= 1e-12, f"effective capacity primal differs by {err:.3e}"


def test_sens_ldmx_ec_matches_finite_differences_in_lo():
    """dEC/dLo, dE/dLo and dEprime/dLo against central differences taken by
    scaling the open-class arrival rate, the only channel into Lo."""
    err = 0.0
    n = 0
    for lam, D, N, Z, mu in LD_MODELS:
        if not np.any(lam > 0):
            continue
        mux = np.hstack([mu, mu[:, -1:]])
        _, _, _, Lo, dEC, dE, dEprime = pfqn_sens_ldmx_ec(lam, D, mux)
        M = D.shape[0]
        h = 1e-7
        for i in range(M):
            if Lo[i] <= 0:
                continue
            # perturb Lo(i) alone by scaling station i's open demands
            Dp, Dm = D.copy(), D.copy()
            openIdx = np.where(np.isinf(N))[0]
            Dp[i, openIdx] *= (1 + h)
            Dm[i, openIdx] *= (1 - h)
            ECp, Ep, Epp = pfqn_sens_ldmx_ec(lam, Dp, mux)[:3]
            ECm, Em, Epm = pfqn_sens_ldmx_ec(lam, Dm, mux)[:3]
            dLo = 2 * h * Lo[i]
            err = max(err, _relerr((ECp[i] - ECm[i]) / dLo, dEC[i]))
            err = max(err, _relerr((Ep[i] - Em[i]) / dLo, dE[i]))
            err = max(err, _relerr((Epp[i] - Epm[i]) / dLo, dEprime[i]))
            n += 1
    assert n > 0
    assert err <= 1e-6, f"dEC/dLo disagreement {err:.3e} over {n} stations"


# =========================================================================
# pfqn_sens_mvaldmx
# =========================================================================

def test_sens_mvaldmx_base_measures_match_pfqn_mvaldmx():
    """A. The primal must reproduce pfqn_mvaldmx entry by entry."""
    err = 0.0
    for lam, D, N, Z, mu in LD_MODELS:
        M = D.shape[0]
        mom = pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M))
        XN, QN, UN, CN = pfqn_mvaldmx(lam, D, N, Z, mu, np.ones(M))[:4]
        err = max(err, _relerr(mom.X, XN), _relerr(mom.Q, QN),
                  _relerr(mom.U, UN), _relerr(mom.R, CN))
    assert err <= 1e-12, f"base measures differ from pfqn_mvaldmx by {err:.3e}"


def test_sens_mvaldmx_matches_finite_differences_of_mvaldmx():
    """B. The differentiated recursion itself, including the load-dependent
    channel dEC/dLo and the open-class channel dLo/dy of eq. (21)."""
    err = 0.0
    n = 0
    h = 1e-6
    for lam, D, N, Z, mu in LD_MODELS:
        M, R = D.shape
        mom = pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M))
        for j in range(M):
            for s in range(R):
                if D[j, s] <= 0:
                    continue
                Dp, Dm = D.copy(), D.copy()
                Dp[j, s] = D[j, s] * (1 + h)
                Dm[j, s] = D[j, s] * (1 - h)
                QNp = pfqn_mvaldmx(lam, Dp, N, Z, mu, np.ones(M))[1]
                QNm = pfqn_mvaldmx(lam, Dm, N, Z, mu, np.ones(M))[1]
                fd = (QNp - QNm) / (2 * h)          # d nbar / dy at y=1
                an = mom.QCovFull[:, :, j, s]
                err = max(err, _relerr(an, fd))
                n += 1
    assert n > 0
    assert err <= 1e-6, f"finite-difference disagreement {err:.3e} over {n} params"


def test_sens_mvaldmx_covariance_is_symmetric():
    """E. Cov[n(i,r),n(j,s)] and Cov[n(j,s),n(i,r)] are computed by
    differentiating two different classes' equations."""
    err = 0.0
    for lam, D, N, Z, mu in LD_MODELS:
        M = D.shape[0]
        err = max(err, pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M)).QCovAsym)
    assert err <= 1e-8, f"raw QCovFull asymmetry {err:.3e}"


def test_sens_mvaldmx_closed_load_independent_limit_matches_sens_mva():
    """C. With lambda=0 and mu=1 the recursion must collapse onto the
    independently validated de Souza e Silva and Muntz recursion."""
    rng = np.random.default_rng(2)
    err = 0.0
    n = 0
    for _ in range(12):
        M = int(rng.integers(1, 4))
        R = int(rng.integers(1, 3))
        D = 0.2 + rng.random((M, R))
        N = rng.integers(1, 4, size=R).astype(float)
        Z = 0.4 * rng.random(R)
        lam = np.zeros(R)
        mu = np.ones((M, int(N.sum())))
        mom = pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M))
        ref = pfqn_sens_mva(D, N, Z)
        err = max(err, _relerr(mom.Q, ref.Q), _relerr(mom.QCov, ref.QCov),
                  _relerr(mom.QVar, ref.QVar), _relerr(mom.QTotVar, ref.QTotVar))
        n += 1
    assert n > 0
    assert err <= 1e-9, f"closed load-independent limit differs by {err:.3e}"


def test_sens_mvaldmx_matches_brute_force_closed_load_dependent():
    """D1. Exact enumeration of the closed load-dependent product form."""
    rng = np.random.default_rng(3)
    err = 0.0
    n = 0
    for _ in range(10):
        M, R = 2, 1
        D = 0.3 + 0.5 * rng.random((M, R))
        N = rng.integers(2, 5, size=R).astype(float)
        Z = 0.3 * rng.random(R)
        lam = np.zeros(R)
        b = int(rng.integers(2, 4))
        mu = np.zeros((M, int(N.sum())))
        for i in range(M):
            for k in range(mu.shape[1]):
                mu[i, k] = min(k + 1, b) * (0.8 + 0.4 * rng.random())
        mom = pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M))
        Qb, QCovb = _brute_ldmx(lam, D, N, Z, mu, 0)
        err = max(err, _relerr(mom.Q, Qb), _relerr(mom.QCov, QCovb))
        n += 1
    assert n > 0
    assert err <= 1e-8, f"brute-force closed LD disagreement {err:.3e} over {n} models"


def test_sens_mvaldmx_matches_brute_force_mixed_load_dependent():
    """D2. Truncated enumeration of the mixed load-dependent product form. This
    is the only check that closes the loop on the identity
    Cov = d nbar / dy in the mixed load-dependent case. The open populations are
    truncated, which converges geometrically, hence the looser tolerance."""
    rng = np.random.default_rng(4)
    err = 0.0
    n = 0
    for _ in range(6):
        M, R = 2, 2
        D = 0.3 + 0.4 * rng.random((M, R))
        N = np.array([float(rng.integers(1, 3)), np.inf])
        Z = np.array([0.3 * rng.random(), 0.0])
        lam = np.array([0.0, 0.05 + 0.1 * rng.random()])
        b = int(rng.integers(1, 3))
        mu = np.zeros((M, int(np.sum(N[np.isfinite(N)]))))
        for i in range(M):
            for k in range(mu.shape[1]):
                mu[i, k] = min(k + 1, b) * (1.0 + 0.3 * rng.random())
        Lo = np.array([float(np.dot(lam, D[i, :])) for i in range(M)])
        if np.any(Lo / mu[:, -1] > 0.4):
            continue
        mom = pfqn_sens_mvaldmx(lam, D, N, Z, mu, np.ones(M))
        Qb, QCovb = _brute_ldmx(lam, D, N, Z, mu, 60)
        err = max(err, _relerr(mom.Q, Qb), _relerr(mom.QCov, QCovb))
        n += 1
    assert n > 0
    assert err <= 5e-5, f"brute-force mixed LD disagreement {err:.3e} over {n} models"


def test_sens_mvaldmx_rejects_all_open_model():
    with pytest.raises(ValueError):
        pfqn_sens_mvaldmx(np.array([0.5]), np.array([[1.0]]), np.array([np.inf]),
                         np.array([0.0]), np.ones((1, 1)), np.ones(1))


def test_sens_mvaldmx_rejects_arrival_rate_on_closed_class():
    with pytest.raises(ValueError):
        pfqn_sens_mvaldmx(np.array([0.5]), np.array([[1.0]]), np.array([2.0]),
                         np.array([0.0]), np.ones((1, 2)), np.ones(1))
