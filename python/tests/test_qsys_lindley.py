"""
Tests for the native Python conditional Lindley recursion algorithms.

The expected values are the ones the MATLAB and JAR twins produce, and were
themselves validated against 4e6 Monte Carlo replications of the corresponding
one-step experiment.
"""

import numpy as np
import pytest

from line_solver.api import qsys


class TestMm1Lindley:
    def test_moments_match_theorem_1(self):
        r = qsys.qsys_mm1_lindley(0.8, 1.0, 2.0, 3)
        assert r['moments'][0, 0] == pytest.approx(1.890206, abs=1e-6)
        assert r['moments'][0, 1] == pytest.approx(5.274485, abs=1e-6)
        assert r['moments'][0, 2] == pytest.approx(18.220680, abs=1e-6)
        assert r['mean'][0] == pytest.approx(1.890206, abs=1e-6)
        assert r['var'][0] == pytest.approx(
            r['moments'][0, 1] - r['moments'][0, 0] ** 2, abs=1e-12)

    def test_moments_from_an_empty_queue(self):
        r = qsys.qsys_mm1_lindley(0.8, 1.0, 0.0, 3)
        assert r['moments'][0, 0] == pytest.approx(0.444444, abs=1e-6)
        assert r['moments'][0, 1] == pytest.approx(0.888889, abs=1e-6)
        assert r['moments'][0, 2] == pytest.approx(2.666667, abs=1e-6)

    def test_corollary_2_agrees_with_theorem_1_at_order_one(self):
        wn = [0.0, 0.1, 1.0, 3.7, 25.0]
        for lam in (0.3, 0.8, 1.0, 1.9):
            for mu in (0.5, 1.0, 2.0):
                r = qsys.qsys_mm1_lindley(lam, mu, wn)
                assert np.allclose(r['mean'], r['moments'][:, 0], atol=1e-9)
                assert np.all(r['var'] >= -1e-12)

    def test_finite_above_saturation(self):
        # the moments are a one-step expectation, so they exist for any load
        r = qsys.qsys_mm1_lindley(2.0, 1.0, 1.0, 3)
        assert np.all(np.isfinite(r['moments']))
        assert r['var'][0] > 0.0

    def test_against_direct_monte_carlo(self):
        rng = np.random.default_rng(12345)
        n = 2_000_000
        lam, mu, wn = 0.8, 1.0, 2.0
        a = rng.exponential(1.0 / lam, n)
        s = rng.exponential(1.0 / mu, n)
        w = np.maximum(wn + s - a, 0.0)
        r = qsys.qsys_mm1_lindley(lam, mu, wn, 3)
        assert r['moments'][0, 0] == pytest.approx(w.mean(), rel=5e-3)
        assert r['moments'][0, 1] == pytest.approx((w ** 2).mean(), rel=1e-2)
        assert r['moments'][0, 2] == pytest.approx((w ** 3).mean(), rel=2e-2)

    def test_vector_input_shape(self):
        wn = [0.0, 1.5, 10.0]
        r = qsys.qsys_mm1_lindley(0.5, 2.0, wn, 4)
        assert r['moments'].shape == (3, 4)
        assert r['mean'].shape == (3,)
        assert r['var'].shape == (3,)

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            qsys.qsys_mm1_lindley(0.8, 1.0, -1.0)
        with pytest.raises(ValueError):
            qsys.qsys_mm1_lindley(-0.8, 1.0, 1.0)
        with pytest.raises(ValueError):
            qsys.qsys_mm1_lindley(0.8, 1.0, 1.0, 0)


class TestHh1Lindley:
    def test_moments_match_theorem_3(self):
        r = qsys.qsys_hh1_lindley([0.5, 2.0], [0.4, 0.6], [1.0, 4.0], [0.7, 0.3],
                                  1.0, 2)
        assert r['mean'][0] == pytest.approx(1.048425, abs=1e-6)
        assert r['moments'][0, 1] == pytest.approx(2.141581, abs=1e-6)

    def test_degenerate_mixture_reduces_to_mm1(self):
        wn = [0.0, 1.0, 2.0, 7.5]
        hh = qsys.qsys_hh1_lindley([0.8], [1.0], [1.0], [1.0], wn, 3)
        mm = qsys.qsys_mm1_lindley(0.8, 1.0, wn, 3)
        assert np.allclose(hh['moments'], mm['moments'], atol=1e-12)

    def test_variance_exceeds_the_mixture_of_phase_variances(self):
        # the phase is itself random, so the between-phase spread of the means adds
        lam, pa = [0.5, 2.0], [0.4, 0.6]
        mu, ps = [1.0, 4.0], [0.7, 0.3]
        mixed = qsys.qsys_hh1_lindley(lam, pa, mu, ps, 1.0)
        avg = 0.0
        for i in range(2):
            for j in range(2):
                avg += pa[i] * ps[j] * qsys.qsys_mm1_lindley(lam[i], mu[j], 1.0)['var'][0]
        assert mixed['var'][0] > avg

    def test_against_direct_monte_carlo(self):
        rng = np.random.default_rng(21)
        n = 2_000_000
        lam, pa = np.array([0.5, 2.0]), [0.4, 0.6]
        mu, ps = np.array([1.0, 4.0]), [0.7, 0.3]
        wn = 1.0
        ia = rng.choice(2, n, p=pa)
        js = rng.choice(2, n, p=ps)
        a = rng.exponential(1.0, n) / lam[ia]
        s = rng.exponential(1.0, n) / mu[js]
        w = np.maximum(wn + s - a, 0.0)
        r = qsys.qsys_hh1_lindley(lam, pa, mu, ps, wn, 2)
        assert r['moments'][0, 0] == pytest.approx(w.mean(), rel=5e-3)
        assert r['moments'][0, 1] == pytest.approx((w ** 2).mean(), rel=1e-2)

    def test_rejects_unnormalized_probabilities(self):
        with pytest.raises(ValueError):
            qsys.qsys_hh1_lindley([1.0, 2.0], [0.5, 0.9], [1.0], [1.0], 0.0)


class TestTandemConditionalMean:
    def test_reference_values(self):
        assert qsys.qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 0.0, 0.0)['mean'][0] \
            == pytest.approx(0.345679, abs=1e-6)
        assert qsys.qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 2.0, 3.0)['mean'][0] \
            == pytest.approx(2.906058, abs=1e-6)
        assert qsys.qsys_mm1_tandem_lindley(0.8, 1.0, 1.5, 5.0, 0.5)['mean'][0] \
            == pytest.approx(0.527153, abs=1e-6)
        assert qsys.qsys_mm1_tandem_lindley(0.5, 2.0, 1.2, 1.0, 4.0)['mean'][0] \
            == pytest.approx(3.486517, abs=1e-6)

    def test_equal_rate_branch_agrees_with_the_limit(self):
        target = qsys.qsys_mm1_tandem_lindley(1.0, 1.0, 1.0, 1.0, 1.0)['mean'][0]
        assert target == pytest.approx(1.084585, abs=1e-6)
        near = qsys.qsys_mm1_tandem_lindley(1.0 + 1e-6, 1.0, 1.0, 1.0, 1.0)['mean'][0]
        assert target == pytest.approx(near, abs=1e-6)

    def test_interdeparture_mean(self):
        # an always-busy upstream server passes on its own service time
        far = qsys.qsys_mm1_tandem_lindley(0.9, 1.0, 1.0, 200.0, 0.0)
        assert far['interdepMean'][0] == pytest.approx(1.0, abs=1e-9)
        assert far['idleProb'][0] == pytest.approx(0.0, abs=1e-9)
        # from an empty upstream queue it is E[max(A,S1)]
        zero = qsys.qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, 0.0, 0.0)
        assert zero['interdepMean'][0] == pytest.approx(1 / 0.8 + 1.0 - 1 / 1.8, abs=1e-9)

    def test_against_direct_monte_carlo_of_the_isolated_step(self):
        rng = np.random.default_rng(31)
        n = 2_000_000
        for lam, mu1, mu2, x, y in [(0.8, 1.0, 1.0, 0.0, 0.0),
                                    (0.8, 1.0, 1.0, 2.0, 3.0),
                                    (1.0, 1.0, 1.0, 1.0, 1.0)]:
            a = rng.exponential(1.0 / lam, n)
            s1 = rng.exponential(1.0 / mu1, n)
            s1p = rng.exponential(1.0 / mu1, n)
            s2 = rng.exponential(1.0 / mu2, n)
            d = np.maximum(a - x - s1, 0.0) + s1p
            w = np.maximum(y + s2 - d, 0.0)
            got = qsys.qsys_mm1_tandem_lindley(lam, mu1, mu2, x, y)['mean'][0]
            assert got == pytest.approx(w.mean(), rel=5e-3)

    def test_rejects_mismatched_shapes(self):
        with pytest.raises(ValueError):
            qsys.qsys_mm1_tandem_lindley(0.8, 1.0, 1.0, [1.0, 2.0], [1.0, 2.0, 3.0])


class TestTandemPathRecursion:
    @staticmethod
    def direct_tandem(a, s):
        """Event-driven reference: arrival at station k is departure from k-1."""
        n, k = s.shape
        epoch = np.concatenate([[0.0], np.cumsum(a[:n - 1])])
        w = np.zeros((n, k))
        dep = np.zeros((n, k))
        for i in range(n):
            arrival = epoch[i]
            for j in range(k):
                prev = dep[i - 1, j] if i > 0 else -np.inf
                w[i, j] = max(prev - arrival, 0.0)
                dep[i, j] = max(arrival, prev) + s[i, j]
                arrival = dep[i, j]
        return w, dep

    def test_matches_direct_event_simulation(self):
        rng = np.random.default_rng(11)
        n, k = 4000, 4
        a = rng.exponential(1 / 0.8, n)
        s = rng.exponential(1.0, (n, k))
        r = qsys.qsys_tandem_lindley(a, s)
        wd, dep = self.direct_tandem(a, s)
        assert np.abs(r['W'] - wd).max() < 1e-9
        assert np.abs(r['departure'] - dep).max() < 1e-9

    def test_first_station_is_plain_lindley(self):
        rng = np.random.default_rng(3)
        n = 500
        a = rng.exponential(1 / 0.7, n)
        s = rng.exponential(1.0, (n, 2))
        r = qsys.qsys_tandem_lindley(a, s)
        w = 0.0
        for i in range(n - 1):
            w = max(w + s[i, 0] - a[i], 0.0)
            assert r['W'][i + 1, 0] == pytest.approx(w, abs=1e-12)

    def test_the_two_interdeparture_forms_agree(self):
        rng = np.random.default_rng(19)
        n, k = 800, 3
        a = rng.exponential(1 / 0.75, n)
        s = rng.exponential(1.0, (n, k))
        r = qsys.qsys_tandem_lindley(a, s)
        w, g = r['W'], r['G']
        for j in range(k - 1):
            idle = np.maximum(g[:n - 1, j] - w[:n - 1, j] - s[:n - 1, j], 0.0) \
                + s[1:, j]
            diff = g[:n - 1, j] + w[1:, j] - w[:n - 1, j] + s[1:, j] - s[:n - 1, j]
            assert np.allclose(g[:n - 1, j + 1], idle, atol=1e-12)
            assert np.allclose(g[:n - 1, j + 1], diff, atol=1e-12)

    def test_published_proposition_1_is_wrong(self):
        # the reference drops S[n+1,k] - S[n,k]; station 1 survives, downstream
        # stations do not. This pins the defect so a future edit cannot reintroduce it
        rng = np.random.default_rng(11)
        n, k = 2000, 2
        a = rng.exponential(1 / 0.8, n)
        s = rng.exponential(1.0, (n, k))
        wd, _ = self.direct_tandem(a, s)

        w = np.zeros((n, k))
        for i in range(n - 1):
            gap = a[i]
            for j in range(k):
                w[i + 1, j] = max(w[i, j] + s[i, j] - gap, 0.0)
                gap = w[i + 1, j] - w[i, j] + gap      # the published form
        assert np.allclose(w[:, 0], wd[:, 0], atol=1e-9), "station 1 is unaffected"
        assert np.abs(w[:, 1] - wd[:, 1]).max() > 1.0, "downstream stations are wrong"

    def test_honours_an_initial_state(self):
        rng = np.random.default_rng(5)
        n = 100
        a = rng.exponential(1.0, n)
        s = rng.exponential(1.0, (n, 3))
        r = qsys.qsys_tandem_lindley(a, s, [2.0, 1.0, 0.5])
        assert np.allclose(r['W'][0, :], [2.0, 1.0, 0.5])

    def test_rejects_mismatched_shapes(self):
        with pytest.raises(ValueError):
            qsys.qsys_tandem_lindley([1.0, 2.0], np.ones((3, 2)))
        with pytest.raises(ValueError):
            qsys.qsys_tandem_lindley([1.0, -2.0], np.ones((2, 2)))
