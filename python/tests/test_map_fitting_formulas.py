"""Drift guard for the closed-form MAP/MMPP/APH fitting formulas.

Two kinds of assertion, both required:

- ROUND TRIP: the fitted process must reproduce the statistics it was asked to
  match. This catches a mis-transcribed expression without needing a golden.
- PINNED VALUES: a few rates from the MATLAB reference, whose MMPP(2) closed
  form was verified symbolically in SageMath. This catches a change that is
  self-consistent but is no longer the same formula.

See _kb/03-api-layer.md, "MMPP(2) moment fitting".
"""
import math

import numpy as np
import pytest

from line_solver.api.mam.map_analysis import (
    map2_fit, map_erlang, map_gamma2, map_hyperexp, map_mean, map_moment,
    map_mmpp2, map_scv,
)
from line_solver.lib.kpctoolbox.aph import aph_fit
from line_solver.api.mam.aph2_fitting import aph2_fit
from line_solver.lib.kpctoolbox.mmpp import mmpp2_fit3


def _pair(M):
    if isinstance(M, dict):
        return np.asarray(M['D0'], float), np.asarray(M['D1'], float)
    return np.asarray(M[0], float), np.asarray(M[1], float)


def _moments(D0, D1, k=3):
    return [map_moment(D0, D1, i) for i in range(1, k + 1)]


MMPP2_CASES = [(1, 5, 45, 0.3), (0.5, 1, 3.6, 0.8), (2, 12, 216, 0.05), (1, 21, 1000, 1e-4)]


@pytest.mark.parametrize("e1,e2,e3,g2", MMPP2_CASES)
def test_mmpp2_fit3_matches_moments_and_decay_rate(e1, e2, e3, g2):
    D0, D1 = mmpp2_fit3(e1, e2, e3, g2)
    m = _moments(D0, D1)
    assert m[0] == pytest.approx(e1, rel=1e-9)
    assert m[1] == pytest.approx(e2, rel=1e-9)
    assert m[2] == pytest.approx(e3, rel=1e-9)
    assert map_gamma2(D0, D1) == pytest.approx(g2, rel=1e-8, abs=1e-12)
    # an MMPP(2) has non-negative rates by construction
    assert D1[0, 0] >= 0 and D1[1, 1] >= 0 and D0[0, 1] >= 0 and D0[1, 0] >= 0


def test_mmpp2_fit3_matches_the_matlab_closed_form():
    # Pinned against MATLAB mmpp2_fit3.m; the JAR and the C++ header agree.
    D0, D1 = mmpp2_fit3(1, 5, 45, 0.3)
    assert D1[0, 0] == pytest.approx(0.11835708875252297, abs=1e-12)
    assert D1[1, 1] == pytest.approx(3.0416429112474770, abs=1e-12)
    assert D0[0, 1] == pytest.approx(0.25333822637151965, abs=1e-12)
    assert D0[1, 0] == pytest.approx(0.58666177362848038, abs=1e-12)


def test_map_mmpp2_realizes_the_requested_lag1_autocorrelation():
    # rho1 = gamma2*(1-1/SCV)/2 holds for every MMPP(2) (proved in SageMath)
    mean, scv, skew, acf1 = 1.0, 4.0, 3.5, 0.2
    D0, D1 = _pair(map_mmpp2(mean, scv, skew, acf1))
    assert map_mean(D0, D1) == pytest.approx(mean, rel=1e-9)
    assert map_scv(D0, D1) == pytest.approx(scv, rel=1e-8)
    assert map_gamma2(D0, D1) == pytest.approx(acf1 / ((1 - 1 / scv) / 2), rel=1e-8)


@pytest.mark.parametrize("e1,e2,e3,g2", [(1, 5, 45, .3), (1, 3, 20, .1), (.5, 1, 3.6, .5), (2, 12, 216, .05)])
def test_map2_fit_matches_moments_and_decay_rate(e1, e2, e3, g2):
    M, _ = map2_fit(e1, e2, e3, g2)
    D0, D1 = _pair(M)
    m = _moments(D0, D1)
    assert m[0] == pytest.approx(e1, rel=1e-9)
    assert m[1] == pytest.approx(e2, rel=1e-9)
    assert m[2] == pytest.approx(e3, rel=1e-9)
    assert map_gamma2(D0, D1) == pytest.approx(g2, rel=1e-8, abs=1e-12)


@pytest.mark.parametrize("e1,e2,e3", [(1, 3, 15), (1, 2, 6.5), (1, 5, 60), (2, 9, 60), (1, 1.2, 1.7)])
def test_aph_fit_matches_three_moments(e1, e2, e3):
    M, _ = aph_fit(e1, e2, e3)
    D0, D1 = _pair(M)
    m = _moments(D0, D1)
    assert m[0] == pytest.approx(e1, rel=1e-9)
    assert m[1] == pytest.approx(e2, rel=1e-9)
    assert m[2] == pytest.approx(e3, rel=1e-9)


@pytest.mark.parametrize("e1,e2,e3", [(1, 3, 15), (1, 2.5, 10), (0.5, 0.75, 1.8)])
def test_aph2_fit_matches_feasible_moments(e1, e2, e3):
    D0, D1 = _pair(aph2_fit(e1, e2, e3))
    m = _moments(D0, D1)
    assert m[0] == pytest.approx(e1, rel=1e-9)
    assert m[1] == pytest.approx(e2, rel=1e-9)
    assert m[2] == pytest.approx(e3, rel=1e-9)


@pytest.mark.parametrize("k", [2, 3, 5])
def test_map_erlang_matches_mean_and_scv(k):
    D0, D1 = _pair(map_erlang(1.5, k))
    assert map_mean(D0, D1) == pytest.approx(1.5, rel=1e-9)
    assert map_scv(D0, D1) == pytest.approx(1.0 / k, rel=1e-8)


def test_map_hyperexp_moment_form_returns_the_feasible_second_root():
    # p=0.5 with SCV=4 is only reachable after the fallback to a smaller p, and
    # only through the second root of the quadratic. Returning None here means
    # the feasible root is being discarded.
    res = map_hyperexp(1.0, 4.0, 0.5)
    assert res is not None
    D0, D1 = _pair(res)
    assert map_mean(D0, D1) == pytest.approx(1.0, rel=1e-8)
    assert map_scv(D0, D1) == pytest.approx(4.0, rel=1e-8)
    # the array form of the same name must keep working
    D0b, D1b = _pair(map_hyperexp([0.5, 0.5], [1.0, 2.0]))
    assert D0b.shape == (2, 2)


def test_generic_mmpp2_identities_hold_numerically():
    # Identities proved in SageMath and relied on by mmpp2_fit/fit1/fit4 and
    # map_mmpp2: g2 = tr(P)-1, SCV >= 1, rho1 = (g2/2)(1-1/SCV), rho_k geometric.
    rng = np.random.default_rng(12345)
    for _ in range(20):
        mu00, mu11, q01, q10 = rng.uniform(0.1, 5.0, size=4)
        D0 = np.array([[-mu00 - q01, q01], [q10, -mu11 - q10]])
        D1 = np.array([[mu00, 0.0], [0.0, mu11]])
        A = np.linalg.inv(-D0)
        P = A @ D1
        w, V = np.linalg.eig(P.T)
        pi = np.real(V[:, np.argmin(abs(w - 1))])
        pi = pi / pi.sum()
        one = np.ones(2)
        e1 = pi @ A @ one
        e2 = 2 * pi @ A @ A @ one
        scv = (e2 - e1 ** 2) / e1 ** 2
        g2 = np.trace(P) - 1
        rho1 = (pi @ A @ P @ A @ one - e1 ** 2) / (e2 - e1 ** 2)
        rho2 = (pi @ A @ P @ P @ A @ one - e1 ** 2) / (e2 - e1 ** 2)
        assert g2 == pytest.approx(mu00 * mu11 / (mu00 * mu11 + mu00 * q10 + mu11 * q01), rel=1e-12)
        assert scv >= 1 - 1e-12
        assert rho1 == pytest.approx(g2 * (1 - 1 / scv) / 2, rel=1e-9)
        assert rho2 == pytest.approx(g2 * rho1, rel=1e-9)


def test_cox2_fit_central_matches_three_moments():
    # The exit probability is determined by the mean once the two rates are
    # fixed: phi = 1 - mu2*e1 + mu2/mu1. The expression that shipped before
    # 2026-07-22 was a different quantity and matched no moment at all.
    from line_solver.distributions.markovian import Cox2
    for e1, e2, e3 in ((1.9, 7.5, 44.85),
                       (0.7333333333333333, 1.2888888888888888, 3.688888888888889),
                       (0.45, 0.43, 0.633)):
        var = e2 - e1 ** 2
        skew = (e3 - 3 * e1 * e2 + 2 * e1 ** 3) / var ** 1.5
        cx = Cox2.fit_central(e1, var, skew)
        D0, D1 = _pair(cx.getRepresentation())
        alpha = np.zeros(D0.shape[0])
        alpha[0] = 1.0
        A = np.linalg.inv(-D0)
        one = np.ones(D0.shape[0])
        m = [math.factorial(k) * float(alpha @ np.linalg.matrix_power(A, k) @ one)
             for k in (1, 2, 3)]
        assert m[0] == pytest.approx(e1, rel=1e-9)
        assert m[1] == pytest.approx(e2, rel=1e-9)
        assert m[2] == pytest.approx(e3, rel=1e-9)


def test_phase_type_skewness_is_not_reported_as_zero():
    # Distribution.getSkew defaults to 0.0; the Markovian classes must override
    # it, otherwise every phase-type silently reports zero skewness.
    from line_solver.distributions.markovian import Coxian, MMPP2
    cx = Coxian([2.0, 0.5], [0.3, 1.0])
    e1, e2, e3 = 1.9, 7.5, 44.85
    var = e2 - e1 ** 2
    assert cx.getSkew() == pytest.approx((e3 - 3 * e1 * e2 + 2 * e1 ** 3) / var ** 1.5, rel=1e-9)
    assert MMPP2(0.5, 2.0, 0.3, 0.7).getSkew() != 0.0
