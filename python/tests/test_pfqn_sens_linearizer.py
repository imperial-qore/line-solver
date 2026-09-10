"""Regression tests for pfqn_sens_linearizer.

pfqn_sens_linearizer approximates the moments E[Q_i], Var[Q_i], Cov[Q_i,Q_j],
E[Q_i^2] and E[Q_i^3] of the per-station total queue lengths of a closed
product-form network by the LINEARIZER-2 / LINEARIZER-3 algorithms of Strelen
(1990), Section 5, equations (5.1)-(5.8).

These tests mirror the MATLAB validator pfqn_sens_linearizer_validate.m. The
routine is an APPROXIMATION, so it must not be held to machine precision. The
checks are therefore of three kinds:

  A. accuracy bands against the exact pfqn_sens_mom, on models small enough for
     the exact lattice. The reference reports relative errors below 2.1% on E[Q],
     4.1% on E[Q^2] and 6.2% on E[Q^3] over its own 51 networks; the bands
     asserted here are of that order. This is the only meaningful statement of
     correctness for an approximation: it must track the exact answer, not equal
     it;
  B. exactness where the approximation degenerates. At a population of one job
     the CORE estimate of the queue lengths one job down is identically zero
     whatever the delta terms are, so the Linearizer equations coincide with the
     exact MVA and every moment must match pfqn_sens_mom to roundoff. This pins
     the derivative algebra independently of the heuristic;
  C. structural invariants that hold for any population: the mean queue lengths
     conserve the population, and they agree with LINE's own pfqn_linearizer,
     which runs the same heuristic without derivatives.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import (
    pfqn_linearizer,
    pfqn_sens_mom,
    pfqn_sens_linearizer,
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


def _random_closed_models(seed=7, trials=60):
    """Mirrors the model sweep of pfqn_sens_linearizer_validate.m."""
    rng = np.random.default_rng(seed)
    out = []
    for trial in range(1, trials + 1):
        M = int(rng.integers(2, 5))
        R = int(rng.integers(1, 3))
        L = 0.2 + rng.random((M, R))
        N = rng.integers(1, 5, size=R).astype(float)
        Z = 0.5 + rng.random(R) if trial % 2 == 0 else np.zeros(R)
        out.append((L, N, Z))
    return out


CLOSED_MODELS = _random_closed_models()


def _one_job_models(seed=8, trials=12):
    """Mirrors the degenerate one-job sweep of the MATLAB validator."""
    rng = np.random.default_rng(seed)
    out = []
    for trial in range(1, trials + 1):
        M = int(rng.integers(2, 5))
        L = 0.2 + rng.random((M, 1))
        Z = np.array([0.4 + rng.random()]) if trial % 2 == 0 else np.array([0.0])
        out.append((L, Z))
    return out


ONE_JOB_MODELS = _one_job_models()


# The bands are the accuracy the reference itself claims in Section 5 over its 51
# networks, so this asserts that the port reproduces the paper's own accuracy
# statement rather than some slacker figure. The model seed is fixed, so the
# model set is deterministic and these are hard regression guards.
BAND_M = 0.021    # r(Q)   < 2.1% in the reference
BAND_M2 = 0.041   # r(Q^2) < 4.1% in the reference
BAND_M3 = 0.062   # r(Q^3) < 6.2% in the reference
# The reference does not report an error on the variance. It is naturally larger
# than the one on E[Q^2] because Var = E[Q^2] - E[Q]^2 is a difference of larger
# numbers, so relative error is amplified; banded here only to catch regressions.
BAND_VAR = 0.08


# =========================================================================
# pfqn_sens_linearizer
# =========================================================================

def test_sens_linearizer_tracks_exact_moments_within_paper_bands():
    """A. The approximation must track the exact pfqn_sens_mom to the accuracy
    the reference claims: 2.1% on E[Q], 4.1% on E[Q^2], 6.2% on E[Q^3]."""
    worst_m = worst_m2 = worst_m3 = worst_var = 0.0
    n = 0
    for L, N, Z in CLOSED_MODELS:
        app = pfqn_sens_linearizer(L, N, Z)
        ex = pfqn_sens_mom(L, N, Z)
        worst_m = max(worst_m, _relerr(app.m, ex.m))
        worst_m2 = max(worst_m2, _relerr(app.M2, ex.M2))
        worst_m3 = max(worst_m3, _relerr(app.M3, ex.M3))
        worst_var = max(worst_var, _relerr(app.Var, ex.Var))
        n += 1
    assert n > 0
    assert worst_m <= BAND_M, (f"E[Q] error {100*worst_m:.3f}% exceeds band "
                               f"{100*BAND_M:.1f}% over {n} models")
    assert worst_var <= BAND_VAR, (f"Var[Q] error {100*worst_var:.3f}% exceeds "
                                   f"band {100*BAND_VAR:.1f}% over {n} models")
    assert worst_m2 <= BAND_M2, (f"E[Q^2] error {100*worst_m2:.3f}% exceeds band "
                                 f"{100*BAND_M2:.1f}% over {n} models")
    assert worst_m3 <= BAND_M3, (f"E[Q^3] error {100*worst_m3:.3f}% exceeds band "
                                 f"{100*BAND_M3:.1f}% over {n} models")


def test_sens_linearizer_is_exact_at_unit_population():
    """B. At one job the CORE estimate of the queue lengths one job down is
    identically zero whatever the delta terms are, so the Linearizer equations
    degenerate to the exact MVA and every moment must match pfqn_sens_mom to
    roundoff. This pins the derivative algebra independently of the heuristic."""
    err = 0.0
    for L, Z in ONE_JOB_MODELS:
        app = pfqn_sens_linearizer(L, np.array([1.0]), Z)
        ex = pfqn_sens_mom(L, np.array([1.0]), Z)
        err = max(err, _relerr(app.m, ex.m), _relerr(app.Var, ex.Var),
                  _relerr(app.M2, ex.M2), _relerr(app.M3, ex.M3),
                  _relerr(app.Cov, ex.Cov))
    assert err <= 1e-9, f"one-job case differs from exact by {err:.3e}"


def test_sens_linearizer_conserves_population():
    """C. Every job is either queued at a station or thinking in the delay."""
    err = 0.0
    for L, N, Z in CLOSED_MODELS:
        app = pfqn_sens_linearizer(L, N, Z)
        R = L.shape[1]
        for r in range(R):
            in_net = float(np.sum(app.Q[:, r]))
            in_delay = float(app.X[0, r] * Z[r])
            err = max(err, abs(in_net + in_delay - N[r]) / max(1.0, N[r]))
    assert err <= 1e-8, f"population conservation violated by {err:.3e}"


def test_sens_linearizer_means_match_pfqn_linearizer():
    """C. The mean queue lengths must agree with LINE's own pfqn_linearizer,
    which runs the same heuristic without derivatives."""
    err = 0.0
    for L, N, Z in CLOSED_MODELS:
        M = L.shape[0]
        app = pfqn_sens_linearizer(L, N, Z)
        QL = pfqn_linearizer(L, N, Z, ['PS'] * M, 1e-10, 500)[0]
        err = max(err, _relerr(np.sum(app.Q, axis=1), np.sum(QL, axis=1)))
    assert err <= 5e-2, f"means differ from pfqn_linearizer by {err:.3e}"


def test_sens_linearizer_covariance_is_symmetric_after_packing():
    """The reported Cov is symmetrized, as a covariance matrix must be; the raw
    asymmetry is reported separately in CovAsym and is a genuine error indicator
    here rather than a roundoff residual, so it is only sanity-bounded."""
    for L, N, Z in CLOSED_MODELS:
        app = pfqn_sens_linearizer(L, N, Z)
        np.testing.assert_allclose(app.Cov, app.Cov.T, rtol=0, atol=1e-15)
        assert app.CovAsym >= 0.0
        assert np.isfinite(app.CovAsym)


def test_sens_linearizer_variance_is_the_cov_diagonal():
    """Internal consistency of the packing: Var[Q_i] = Cov[Q_i,Q_i]."""
    for L, N, Z in CLOSED_MODELS:
        app = pfqn_sens_linearizer(L, N, Z)
        np.testing.assert_allclose(app.Var, np.diag(app.Cov), rtol=1e-12,
                                   atol=1e-14)


def test_sens_linearizer_empty_population_is_zero():
    L = np.array([[1.0, 2.0], [0.5, 0.3]])
    app = pfqn_sens_linearizer(L, np.array([0.0, 0.0]), np.array([1.0, 1.0]))
    np.testing.assert_allclose(app.Q, 0.0, atol=0)
    np.testing.assert_allclose(app.m, 0.0, atol=0)
    np.testing.assert_allclose(app.Var, 0.0, atol=0)
    np.testing.assert_allclose(app.M2, 0.0, atol=0)
    np.testing.assert_allclose(app.M3, 0.0, atol=0)
    assert app.CovAsym == 0.0
    assert app.iter == 0


def test_sens_linearizer_rejects_open_population():
    with pytest.raises(ValueError):
        pfqn_sens_linearizer(np.array([[1.0]]), np.array([np.inf]),
                             np.array([0.0]))
