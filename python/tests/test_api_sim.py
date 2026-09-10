"""
Tests for the native Python output-analysis package.

The expected values are the ones the MATLAB and JAR twins produce; where an
external reference exists the test compares against it directly, which is what
validates the Shapiro-Wilk port.
"""

import numpy as np
import pytest
from scipy.stats import shapiro as scipy_shapiro
from scipy.stats import t as scipy_t

from line_solver.api import sim


def exponential_sample(n, seed):
    """A reproducible standard exponential stream."""
    return np.random.default_rng(seed).exponential(1.0, n)


class TestDistributions:
    def test_normal_and_t_quantiles(self):
        assert sim.normcdf(1.96) == pytest.approx(0.9750021048517795, abs=1e-12)
        assert sim.norminv(0.975) == pytest.approx(1.959963984540054, abs=1e-10)
        assert sim.tinv(0.975, 10) == pytest.approx(2.2281388519649385, abs=1e-9)
        assert sim.tinv(0.975, 29) == pytest.approx(2.045229642132703, abs=1e-9)
        assert sim.tinv(0.025, 19) == pytest.approx(-2.093024054408309, abs=1e-9)

    def test_norminv_endpoints(self):
        assert sim.norminv(0.0) == float('-inf')
        assert sim.norminv(1.0) == float('inf')
        with pytest.raises(ValueError):
            sim.norminv(1.5)


class TestShapiroWilk:
    def test_against_matlab_reference_values(self):
        x12 = [2.1, 3.4, 1.9, 5.6, 3.3, 4.8, 2.2, 6.1, 3.9, 4.4, 2.8, 5.1]
        r = sim.shapirowilk(x12)
        assert r['W'] == pytest.approx(0.9500369248, abs=1e-9)
        assert r['pvalue'] == pytest.approx(0.6375246422, abs=1e-8)
        assert not r['reject']

        # n = 3 uses the exact null law rather than Royston's transform
        r3 = sim.shapirowilk([1, 2, 5])
        assert r3['W'] == pytest.approx(0.9230769231, abs=1e-9)
        assert r3['pvalue'] == pytest.approx(0.4632628749, abs=1e-8)
        assert np.isnan(r3['zscore'])

        # 4 <= n <= 11 uses the middle branch
        r7 = sim.shapirowilk([1, 2, 3, 4, 5, 6, 20])
        assert r7['W'] == pytest.approx(0.7086833563, abs=1e-9)
        assert r7['pvalue'] == pytest.approx(0.0046079816, abs=1e-9)
        assert r7['reject']

    def test_against_scipy_across_all_three_branches(self):
        rng = np.random.default_rng(7)
        worst_w = worst_p = 0.0
        for n in [3, 4, 5, 6, 7, 8, 11, 12, 16, 32, 100, 500, 2000]:
            for trial in range(4):
                x = (rng.standard_normal(n) if trial % 2 == 0
                     else rng.exponential(1.0, n))
                mine = sim.shapirowilk(x)
                sw, sp = scipy_shapiro(x)
                worst_w = max(worst_w, abs(mine['W'] - sw))
                if n > 3:
                    worst_p = max(worst_p, abs(mine['pvalue'] - sp))
        assert worst_w < 1e-8, "W disagrees with scipy by %g" % worst_w
        assert worst_p < 1e-6, "the p-value disagrees with scipy by %g" % worst_p

    def test_weights_are_antisymmetric(self):
        for n in [3, 4, 5, 6, 12, 32, 101]:
            a = sim.shapirowilk_weights(n)
            assert np.allclose(a, -a[::-1], atol=1e-12)
            if n % 2 == 1:
                assert abs(a[n // 2]) < 1e-12

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            sim.shapirowilk([1, 2])
        with pytest.raises(ValueError):
            sim.shapirowilk(np.ones(10))
        with pytest.raises(ValueError):
            sim.shapirowilk(np.zeros(6000))


class TestVonNeumann:
    def test_against_matlab_reference_values(self):
        r = sim.vonneumann([1, 5, 2, 8, 3, 9, 4, 7, 2, 6, 3, 8])
        assert r['ratio'] == pytest.approx(2.8285714286, abs=1e-9)
        assert r['zscore'] == pytest.approx(1.5666355475, abs=1e-8)
        assert r['pvalue'] == pytest.approx(0.1171999039, abs=1e-8)
        assert not r['reject']

    def test_null_moments_hold_under_iid_normal(self):
        # E[ratio] = 2 and Var = 4(b-2)/((b-1)(b+1)) were confirmed by Monte Carlo
        rng = np.random.default_rng(3)
        for b in (16, 32):
            ratios = np.array([sim.vonneumann(rng.standard_normal(b))['ratio']
                               for _ in range(4000)])
            claim_var = 4.0 * (b - 2) / ((b - 1) * (b + 1))
            assert ratios.mean() == pytest.approx(2.0, abs=0.03)
            assert ratios.var() == pytest.approx(claim_var, rel=0.10)

    def test_rejects_a_random_walk(self):
        rng = np.random.default_rng(11)
        walk = np.cumsum(rng.standard_normal(60))
        r = sim.vonneumann(walk)
        assert r['reject']
        # positive correlation shrinks the successive differences
        assert r['ratio'] < 2.0

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            sim.vonneumann([1, 2])
        with pytest.raises(ValueError):
            sim.vonneumann(np.ones(8))


class TestStsQuantileAreas:
    def test_recovers_the_variance_parameter_on_iid_data(self):
        # for i.i.d. Exp(1), sigma_p^2 = p(1-p)/f(y_p)^2 = p/(1-p) exactly
        b, m, p = 32, 4000, 0.9
        truth = p / (1.0 - p)
        aps, nps = [], []
        for r in range(20):
            y = exponential_sample(b * m, 1000 + r)
            s = sim.sts_quantile_areas(y, b, m, p)
            aps.append(s['Ap'])
            nps.append(s['Np'])
            assert s['areas'].size == b
            assert s['bqe'].size == b
            assert s['n'] == b * m
        assert np.mean(aps) == pytest.approx(truth, rel=0.15)
        assert np.mean(nps) == pytest.approx(truth, rel=0.15)

    def test_quantile_is_an_order_statistic(self):
        y = exponential_sample(6000, 31)
        p = 0.75
        s = sim.sts_quantile_areas(y, 6, 1000, p)
        assert s['quantile'] == np.sort(y)[int(np.ceil(6000 * p)) - 1]

    def test_combined_estimator_definition(self):
        y = exponential_sample(8000, 17)
        s = sim.sts_quantile_areas(y, 8, 1000, 0.5)
        b = s['b']
        expected = (b * s['Ap'] + (b - 1) * s['Np']) / (2 * b - 1)
        assert s['Vp'] == pytest.approx(expected, abs=1e-12)

    def test_single_batch_has_no_between_batch_degrees_of_freedom(self):
        y = exponential_sample(2000, 5)
        s = sim.sts_quantile_areas(y, 1, 2000, 0.5)
        assert s['areas'].size == 1
        assert np.isfinite(s['Ap'])
        assert np.isnan(s['Np'])
        assert np.isnan(s['Vp'])

    def test_prefix_quantiles_are_exact(self):
        # the Fenwick descent must reproduce a brute-force prefix quantile
        rng = np.random.default_rng(23)
        m, p = 200, 0.7
        batch = rng.exponential(1.0, m)
        s = sim.sts_quantile_areas(batch, 1, m, p)
        full = np.sort(batch)[int(np.ceil(m * p)) - 1]
        acc = 0.0
        for k in range(1, m + 1):
            prefix = np.sort(batch[:k])[int(np.ceil(p * k)) - 1]
            acc += k * (full - prefix)
        expected = sim.DEFAULT_WEIGHT * acc / (m * np.sqrt(m))
        assert s['areas'][0] == pytest.approx(expected, abs=1e-12)

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            sim.sts_quantile_areas(np.zeros(100), 3, 10, 0.5)
        with pytest.raises(ValueError):
            sim.sts_quantile_areas(exponential_sample(100, 1), 10, 10, 1.5)


class TestFquest:
    def test_brackets_the_exact_quantile_on_iid_data(self):
        for p in (0.5, 0.9, 0.99):
            truth = -np.log(1.0 - p)
            covered = 0
            reps = 8
            for r in range(reps):
                y = exponential_sample(200000, 9000 + r)
                res = sim.fquest(y, p, 0.05)
                assert res['estimate'] == pytest.approx(truth, rel=0.05)
                assert res['upper'] > res['lower']
                assert res['truncated'] > 0
                assert res['R'] == 1
                if res['lower'] <= truth <= res['upper']:
                    covered += 1
            assert covered >= reps - 2, "p=%s covered %d/%d" % (p, covered, reps)

    def test_delivered_half_width_is_the_combined_estimator(self):
        y = exponential_sample(60000, 424242)
        res = sim.fquest(y, 0.9, 0.05)
        assert res['halfwidth'] == pytest.approx((res['upper'] - res['lower']) / 2.0)
        if not res['heuristic']:
            expected = scipy_t.ppf(0.975, 2 * res['b'] - 1) \
                * np.sqrt(res['Vp'] / res['n'])
            assert res['halfwidth'] == pytest.approx(expected, abs=1e-12)
            assert res['warnings'] == []

    def test_narrower_interval_from_a_longer_path(self):
        short = sim.fquest(exponential_sample(50000, 1), 0.9, 0.05)
        long_ = sim.fquest(exponential_sample(400000, 1), 0.9, 0.05)
        assert long_['halfwidth'] < short['halfwidth']

    def test_rejects_bad_input(self):
        with pytest.raises(ValueError):
            sim.fquest(exponential_sample(1000, 1), 1.0)
        with pytest.raises(ValueError):
            sim.fquest(exponential_sample(1000, 1), 0.5, 0.0)
        with pytest.raises(ValueError):
            sim.fquest(exponential_sample(1000, 1), 0.5, 0.05, {'nosuchoption': 1})


class TestFirquest:
    def test_brackets_the_exact_quantile_on_iid_data(self):
        p = 0.9
        truth = -np.log(1.0 - p)
        R = 5
        covered = 0
        reps = 6
        for rep in range(reps):
            y = np.column_stack([exponential_sample(40000, 77000 + rep * R + r)
                                 for r in range(R)])
            res = sim.firquest(y, p, 0.05)
            assert res['R'] == R
            assert res['n'] == R * res['b'] * res['m']
            assert res['upper'] > res['lower']
            if res['lower'] <= truth <= res['upper']:
                covered += 1
        assert covered >= reps - 2, "covered %d/%d" % (covered, reps)

    def test_truncates_every_replication(self):
        y = np.column_stack([exponential_sample(30000, 555 + r) for r in range(4)])
        res = sim.firquest(y, 0.5, 0.05)
        assert res['truncated'] > 0
        assert res['b'] * res['m'] <= 30000 - res['truncated']

    def test_default_batch_counts(self):
        assert sim.firquest_batchcounts(2) == [14, 11, 8, 5]
        assert sim.firquest_batchcounts(3) == [10, 8, 6, 4]
        assert sim.firquest_batchcounts(4) == [6, 5, 4, 3]
        assert sim.firquest_batchcounts(5) == [5, 4, 3, 2]
        assert sim.firquest_batchcounts(9) == [5, 4, 3, 2]
        assert sim.firquest_batchcounts(10) == [4, 3, 2, 1]
        assert sim.firquest_batchcounts(16) == [4, 3, 2, 1]
        assert sim.firquest_batchcounts(17) == [3, 2, 1]
        assert sim.firquest_batchcounts(23) == [2, 1]
        assert sim.firquest_batchcounts(33) == [1]
        assert sim.firquest_batchcounts(100) == [1]
        with pytest.raises(ValueError):
            sim.firquest_batchcounts(1)

    def test_rejects_a_single_path(self):
        with pytest.raises(ValueError):
            sim.firquest(exponential_sample(1000, 1), 0.5)
