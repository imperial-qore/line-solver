"""
Tests for the native Python Ciucu-Mehri tandem tail bounds.

The load-bearing check is the exact reduction: for M/M/1 -> ./M/1 the five
inequalities of Lemma 4 hold as equalities, so the bound must return the exact
tails, (1 + theta x) e^{-theta x} for the sojourn time and the Kraemer form for
the waiting time. The remaining expected values are the digits the MATLAB and
JAR twins produce, themselves checked against an exact CTMC reference for the
Erlang(2)/M/1 tandem.
"""

from math import exp

import numpy as np
import pytest

from line_solver.api import qsys


def exp_lst(rate):
    return lambda s: rate / (rate + s)


def exp_dlst(rate):
    return lambda s: rate / (rate + s) ** 2


def det_lst(d):
    return lambda s: exp(-s * d)


def det_dlst(d):
    return lambda s: d * exp(-s * d)


def erlang2_lst(rate):
    return lambda s: (rate / (rate + s)) ** 2


def erlang2_dlst(rate):
    return lambda s: 2.0 * rate ** 2 / (rate + s) ** 3


class TestTandemUbCiucu:
    def test_exact_for_mm1_tandem(self):
        mu, lam, theta = 1.0, 0.5, 0.5
        x = np.array([1.0, 5.0, 10.0])
        r = qsys.qsys_tandem_ub_ciucu(x, exp_lst(lam), 1.0, mu, exp_dlst(lam))
        assert r['theta'] == pytest.approx(theta, abs=1e-12)
        assert r['A'] == pytest.approx(1.0)
        assert r['B'] == pytest.approx(0.0, abs=1e-12)
        assert r['D'] == 0.0
        S = (1.0 + theta * x) * np.exp(-theta * x)
        W = (1.0 - 2.0 * theta ** 2 / (mu * (mu + theta))
             + x * (mu - theta) * theta / (mu + theta)) * np.exp(-theta * x)
        assert r['S'] == pytest.approx(S, abs=1e-9)
        assert r['W'] == pytest.approx(W, abs=1e-9)

    def test_dm1_tandem_matches_matlab(self):
        d = 4.0 / 3.0                       # utilization 3/4, unit service rate
        r = qsys.qsys_tandem_ub_ciucu([5.0, 10.0, 20.0], det_lst(d), 1.0, 1.0, det_dlst(d))
        assert r['theta'] == pytest.approx(0.454394983439251, abs=1e-10)
        # Arrivals less variable than Poisson put the bound on the D > 0 branch.
        assert r['D'] > 0.0
        assert r['S'] == pytest.approx(
            [0.493699069414895, 0.09117765981667, 0.00182565464670795], abs=1e-9)
        assert r['W'][0] == pytest.approx(0.249922122549679, abs=1e-9)

    def test_erlang_m1_tandem_matches_matlab(self):
        r = qsys.qsys_tandem_ub_ciucu([10.0, 20.0], erlang2_lst(1.0), 1.0, 1.0,
                                      erlang2_dlst(1.0))
        assert r['theta'] == pytest.approx(0.618033988749894, abs=1e-10)
        assert r['S'] == pytest.approx([0.0170463897138194, 6.62788942098923e-05], abs=1e-12)

    def test_hyperexponential_service_matches_matlab(self):
        p2, p1, m2 = 0.9, 0.1, 1.69
        m1 = p1 * m2 / (m2 - p2)            # CV(Y) = 2 with E[Y] = 1
        mean_y = p1 / m1 + p2 / m2
        rate = 2.0 / (mean_y / 0.5)         # Erlang(2) arrivals at rho = 1/2
        r = qsys.qsys_tandem_ub_ciucu([10.0, 25.0, 50.0], erlang2_lst(rate), [p1, p2],
                                      [m1, m2], erlang2_dlst(rate))
        assert r['theta'] == pytest.approx(0.1500514416369, abs=1e-10)
        assert r['S'] == pytest.approx(
            [0.524484692983612, 0.0877479694886629, 0.0033336271735292], abs=1e-9)
        assert np.all(np.isnan(r['W']))     # the W form is Exp-service only

    def test_numerical_derivative_matches_analytic(self):
        d = 4.0 / 3.0
        ref = qsys.qsys_tandem_ub_ciucu([5.0, 10.0], det_lst(d), 1.0, 1.0, det_dlst(d))
        num = qsys.qsys_tandem_ub_ciucu([5.0, 10.0], det_lst(d), 1.0, 1.0)
        assert num['alpha'] == pytest.approx(ref['alpha'], rel=1e-10)
        assert num['S'] == pytest.approx(ref['S'], rel=1e-10)

    def test_bound_dominates_a_simulated_tandem(self):
        # D/M/1 -> ./M/1 at rho = 3/4: a Lindley replay must stay under the bound.
        # Sojourn times along one path are strongly correlated, so a single run can
        # overshoot a tail of 1e-2 by a tenth; the check is on the mean of ten
        # independent replications, with a three-standard-error margin.
        n, d, mu, reps = 200000, 4.0 / 3.0, 1.0, 10
        x = np.array([1.0, 2.0, 5.0, 10.0])
        arrivals = np.cumsum(np.full(n, d))
        keep = slice(n // 10, n)
        emp = np.zeros((reps, x.size))
        for rep in range(reps):
            rng = np.random.default_rng(23000 + rep)
            cy = np.cumsum(rng.exponential(1.0 / mu, n))
            cz = np.cumsum(rng.exponential(1.0 / mu, n))
            d1 = np.maximum.accumulate(arrivals - np.concatenate(([0.0], cy[:-1]))) + cy
            d2 = np.maximum.accumulate(d1 - np.concatenate(([0.0], cz[:-1]))) + cz
            sojourn = d2[keep] - arrivals[keep]
            emp[rep] = [np.mean(sojourn > xi) for xi in x]
        r = qsys.qsys_tandem_ub_ciucu(x, det_lst(d), 1.0, mu, det_dlst(d))
        mean = emp.mean(axis=0)
        se = emp.std(axis=0, ddof=1) / np.sqrt(reps)
        assert np.all(r['S'] >= mean - 3.0 * se)

    def test_rejects_unstable_and_malformed_input(self):
        with pytest.raises(ValueError):
            qsys.qsys_tandem_ub_ciucu([1.0], exp_lst(2.0), 1.0, 1.0, exp_dlst(2.0))
        with pytest.raises(ValueError):
            qsys.qsys_tandem_ub_ciucu([-1.0], exp_lst(0.5), 1.0, 1.0, exp_dlst(0.5))
        with pytest.raises(ValueError):
            qsys.qsys_tandem_ub_ciucu([1.0], exp_lst(0.5), [0.3, 0.3], [1.0, 2.0])
