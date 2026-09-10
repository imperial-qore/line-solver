"""Regression tests for pfqn_sens_mom.

pfqn_sens_mom returns the exact moments E[Q_i], E[Q_i^2], E[Q_i^3] and the
covariances Cov[Q_i,Q_j] of the per-station TOTAL queue lengths of a closed
product-form network, by second-order forward-mode differentiation of the MVA
recursion (Strelen 1990, Theorems 2.1/3.1/3.2/3.5 and equation (3.2)).

These tests mirror the MATLAB validator pfqn_sens_mom_validate.m. Five
references, chosen so that no channel of the derivation is checked against
itself:

  A. brute-force enumeration of the closed product-form distribution, which is
     ground truth for m, Var, Cov, E[Q^2] and E[Q^3];
  B. pfqn_sens_mva. Summing its per-class covariance matrix at station i over all
     class pairs must give Var[Q_i], since Var[sum_r n(i,r)] = sum_{r,s}
     Cov[n(i,r),n(i,s)]. This ties the per-station-total moments of Strelen to
     the finer per-class moments of de Souza e Silva and Muntz;
  C. pfqn_mva for the base measures;
  D. Cov symmetry: x_j dm_i/dx_j and x_i dm_j/dx_i are computed by different
     derivative tracks and must agree;
  E. the published table of Example 3.4 of the reference (the Kobayashi
     central-server model), which pins the second derivative against numbers the
     author printed rather than against our own code.
"""

import itertools

import numpy as np
import pytest

from line_solver.api.pfqn import (
    pfqn_mva,
    pfqn_sens_mva,
    pfqn_sens_mom,
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


def _brute_totals(L, N, Z):
    """Moments of the per-station TOTAL queue lengths by enumerating the closed
    product-form equilibrium distribution:

      p(n) ~ prod_i [ n_i! prod_r L(i,r)^n(i,r)/n(i,r)! ]
             * prod_r Z(r)^n(0,r)/n(0,r)!

    This shares none of the differentiated-MVA algebra, so it is ground truth for
    every moment pfqn_sens_mom reports.
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
    tot = S.sum(axis=2)                             # (K, M)

    m = np.array([float(np.sum(w * tot[:, i])) for i in range(M)])
    M2 = np.array([float(np.sum(w * tot[:, i] ** 2)) for i in range(M)])
    M3 = np.array([float(np.sum(w * tot[:, i] ** 3)) for i in range(M)])
    Var = M2 - m ** 2
    Cov = np.zeros((M, M))
    for i in range(M):
        for j in range(M):
            Cov[i, j] = float(np.sum(w * tot[:, i] * tot[:, j])) - m[i] * m[j]
    return m, Var, Cov, M2, M3


def _brute_perclass(L, N, Z):
    """Per-class moments of n(i,r) by enumerating the closed product form. This
    is the ground truth for the ``groups = 1..R`` setting, i.e. Akyildiz and
    Strelen's Theorem 1 with T = {r}."""
    from scipy.special import gammaln

    L = np.asarray(L, float)
    N = np.asarray(N, int)
    Z = np.asarray(Z, float)
    M, R = L.shape

    per = [_compositions_leq(int(N[r]), M) for r in range(R)]
    states = []
    weights = []
    for combo in itertools.product(*per):
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
    m = np.tensordot(w, S, axes=(0, 0))
    M2 = np.tensordot(w, S.astype(float) ** 2, axes=(0, 0))
    M3 = np.tensordot(w, S.astype(float) ** 3, axes=(0, 0))
    Var = M2 - m ** 2
    return m, Var, M2, M3


def _random_closed_models(seed=3, trials=40):
    """Mirrors the model sweep of pfqn_sens_mom_validate.m."""
    rng = np.random.default_rng(seed)
    out = []
    for trial in range(1, trials + 1):
        M = int(rng.integers(1, 4))
        R = int(rng.integers(1, 3))
        L = 0.2 + rng.random((M, R))
        N = rng.integers(0, 4, size=R)
        if not np.any(N > 0):
            N[0] = 2
        Z = 0.3 + rng.random(R) if trial % 2 == 0 else np.zeros(R)
        mi = (rng.integers(1, 4, size=M).astype(float) if trial % 3 == 0
              else np.ones(M))
        out.append((L, N.astype(float), Z, mi))
    return out


CLOSED_MODELS = _random_closed_models()


# =========================================================================
# pfqn_sens_mom
# =========================================================================

def test_sens_mom_base_measures_match_pfqn_mva():
    """C. The primal must reproduce pfqn_mva entry by entry, otherwise the
    derivative is of the wrong function."""
    err = 0.0
    for L, N, Z, mi in CLOSED_MODELS:
        mom = pfqn_sens_mom(L, N, Z, mi)
        # pfqn_mva returns (XN, CN, QN, UN, RN, TN, AN); RN is the (M x R)
        # residence time, i.e. the MATLAB pfqn_mva 4th output that mom.R mirrors
        r = pfqn_mva(L, N, Z, mi)
        XN, QN, UN, RN = r[0], r[2], r[3], r[4]
        err = max(err, _relerr(mom.X, XN), _relerr(mom.Q, QN),
                  _relerr(mom.U, UN), _relerr(mom.R, RN))
    assert err <= 1e-10, f"base measures differ from pfqn_mva by {err:.3e}"


def test_sens_mom_matches_brute_force_product_form():
    """A. Ground truth: enumeration of the closed product-form distribution,
    covering m, Var, Cov, E[Q^2] and E[Q^3]."""
    err = 0.0
    n = 0
    for L, N, Z, mi in CLOSED_MODELS:
        if np.prod(N + 1) > 32 or L.shape[0] > 3 or not np.all(mi == 1):
            continue
        mom = pfqn_sens_mom(L, N, Z, mi)
        mb, Varb, Covb, M2b, M3b = _brute_totals(L, N, Z)
        err = max(err, _relerr(mom.m, mb), _relerr(mom.Var, Varb),
                  _relerr(mom.Cov, Covb), _relerr(mom.M2, M2b),
                  _relerr(mom.M3, M3b))
        n += 1
    assert n > 0
    assert err <= 1e-9, f"brute-force disagreement {err:.3e} over {n} models"


def test_sens_mom_total_variance_matches_sens_mva():
    """B. Var[sum_r n(i,r)] = sum_{r,s} Cov[n(i,r),n(i,s)], which ties the
    per-station-total moments of Strelen to the per-class moments of de Souza e
    Silva and Muntz. Includes mi > 1, which brute force does not cover."""
    err = 0.0
    for L, N, Z, mi in CLOSED_MODELS:
        mom = pfqn_sens_mom(L, N, Z, mi)
        ref = pfqn_sens_mva(L, N, Z, mi)
        err = max(err, _relerr(mom.Var, ref.QTotVar))
    assert err <= 1e-9, f"Var vs pfqn_sens_mva QTotVar differs by {err:.3e}"


def test_sens_mom_perclass_grouping_matches_brute_force_and_sens_mva():
    """F. groups = 1..R scales one class at a time, which is Akyildiz and
    Strelen's Theorem 1 with T = {r}. It must reproduce the per-class moments of
    the brute-force distribution, INCLUDING the third, and its second moments
    must equal pfqn_sens_mva's exactly."""
    err = 0.0
    n = 0
    for L, N, Z, mi in CLOSED_MODELS:
        if np.prod(N + 1) > 32 or L.shape[0] > 3 or not np.all(mi == 1):
            continue
        R = L.shape[1]
        momc = pfqn_sens_mom(L, N, Z, mi, np.arange(1, R + 1))
        mc, Varc, M2c, M3c = _brute_perclass(L, N, Z)
        err = max(err, _relerr(momc.m, mc), _relerr(momc.Var, Varc),
                  _relerr(momc.M2, M2c), _relerr(momc.M3, M3c))
        # the second moments must agree with the independently validated
        # de Souza e Silva and Muntz recursion
        err = max(err, _relerr(momc.Var, pfqn_sens_mva(L, N, Z, mi).QVar))
        n += 1
    assert n > 0
    assert err <= 1e-9, f"per-class grouping disagreement {err:.3e} over {n} models"


def test_sens_mom_grouping_shapes_and_degenerate_partitions():
    """G. A grouping that puts every class in ONE group must reproduce the
    default per-station totals exactly, and the group axis must appear only when
    G > 1."""
    err = 0.0
    n = 0
    for L, N, Z, mi in CLOSED_MODELS:
        M, R = L.shape
        mom = pfqn_sens_mom(L, N, Z, mi)
        # the single-group case is the documented collapsed shape
        assert mom.m.shape == (M,)
        assert mom.Cov.shape == (M, M)
        assert mom.dm.shape == (M, M)
        if R == 2:
            momg = pfqn_sens_mom(L, N, Z, mi, np.array([1, 1]))
            err = max(err, _relerr(momg.m, mom.m), _relerr(momg.Var, mom.Var),
                      _relerr(momg.M3, mom.M3), _relerr(momg.Cov, mom.Cov))
            n += 1
        momc = pfqn_sens_mom(L, N, Z, mi, np.arange(1, R + 1))
        if R > 1:
            assert momc.m.shape == (M, R)
            assert momc.Cov.shape == (M, R, M, R)
            assert momc.dm.shape == (M, R, M, R)
            # the per-class means must sum back to the station totals
            err = max(err, _relerr(np.sum(momc.m, axis=1), mom.m))
    assert n > 0
    assert err <= 1e-9, f"degenerate grouping disagreement {err:.3e}"


def test_sens_mom_rejects_malformed_groups():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    N, Z = np.array([2.0, 1.0]), np.array([0.0, 0.0])
    with pytest.raises(ValueError):
        pfqn_sens_mom(L, N, Z, None, np.array([1]))            # wrong length
    with pytest.raises(ValueError):
        pfqn_sens_mom(L, N, Z, None, np.array([0, 1]))         # not 1-based
    with pytest.raises(ValueError):
        pfqn_sens_mom(L, N, Z, None, np.array([1, 3]))         # empty group 2


def test_sens_mom_raw_covariance_is_symmetric():
    """D. x_j dm_i/dx_j and x_i dm_j/dx_i are distinct expressions carried along
    different derivative tracks; their agreement before symmetrization is an
    independent residual of the recursion."""
    err = 0.0
    for L, N, Z, mi in CLOSED_MODELS:
        err = max(err, pfqn_sens_mom(L, N, Z, mi).CovAsym)
    assert err <= 1e-9, f"raw Cov asymmetry {err:.3e}"


def test_sens_mom_strelen_example_34_published_table():
    """E. Strelen Example 3.4, the Kobayashi central-server model: 12 type-1
    queues, one class. Queues 1-9: x=0.0215, e=9.333; queues 10,11: x=0.104,
    e=10.5; queue 12: x=0.019, e=105. The paper prints E(Q_i) and sigma^2_{Q_i}
    for n = 3, 2, 1 to five significant digits.

    This pins the second derivative against numbers the author printed rather
    than against our own code.
    """
    xs = np.concatenate([np.full(9, 0.0215), [0.104, 0.104, 0.019]])
    es = np.concatenate([np.full(9, 9.333), [10.5, 10.5, 105.0]])
    Lk = (xs * es).reshape(-1, 1)          # demands, 12 x 1
    # rows: queues 1-9, 10-11, 12; columns: n = 3, 2, 1
    paperM = np.array([[0.07606, 0.05835, 0.03353],
                       [0.53316, 0.36327, 0.18246],
                       [1.24917, 0.74835, 0.33334]])
    paperVar = np.array([[0.07893, 0.05873, 0.03240],
                         [0.57689, 0.34341, 0.14917],
                         [1.02546, 0.56250, 0.22222]])
    err = 0.0
    for col in range(3):
        n_jobs = 3 - col                   # col 0 -> n=3, col 1 -> n=2, col 2 -> n=1
        mk = pfqn_sens_mom(Lk, np.array([float(n_jobs)]), np.array([0.0]))
        got = np.array([mk.m[0], mk.m[9], mk.m[11]])
        gotV = np.array([mk.Var[0], mk.Var[9], mk.Var[11]])
        err = max(err, _relerr(got, paperM[:, col]),
                  _relerr(gotV, paperVar[:, col]))
        # the nine identical queues must be identical, and so must 10 and 11
        err = max(err, _relerr(mk.m[0:9], np.full(9, mk.m[0])),
                  _relerr(mk.m[9], mk.m[10]))
    assert err <= 5e-5, f"Strelen Example 3.4 disagreement {err:.3e}"


def test_sens_mom_qvar_is_bernoulli_at_unit_population():
    """With N=1 the per-station total queue length is Bernoulli, so
    Var = m(1-m) and E[Q^k] = m for every k. An analytic check independent of
    every other reference."""
    L = np.array([[1.0], [0.5]])
    mom = pfqn_sens_mom(L, np.array([1.0]), np.array([0.0]))
    m = np.asarray(mom.m).flatten()
    np.testing.assert_allclose(mom.Var, m * (1.0 - m), rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(mom.M2, m, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(mom.M3, m, rtol=1e-12, atol=1e-14)


def test_sens_mom_empty_population_is_zero():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    mom = pfqn_sens_mom(L, np.array([0.0, 0.0]), np.array([1.0, 1.0]))
    np.testing.assert_allclose(mom.Q, 0.0, atol=0)
    np.testing.assert_allclose(mom.m, 0.0, atol=0)
    np.testing.assert_allclose(mom.Var, 0.0, atol=0)
    np.testing.assert_allclose(mom.Cov, 0.0, atol=0)
    np.testing.assert_allclose(mom.M2, 0.0, atol=0)
    np.testing.assert_allclose(mom.M3, 0.0, atol=0)
    assert mom.CovAsym == 0.0
    assert np.all(np.isnan(mom.Skew))


def test_sens_mom_rejects_open_population():
    with pytest.raises(ValueError):
        pfqn_sens_mom(np.array([[1.0]]), np.array([np.inf]), np.array([0.0]))
