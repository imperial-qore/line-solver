"""Normalizing-constant route for mixed limited load-dependent networks (pfqn_ncldmx).

pfqn_ncldmx states the Bruell-Balbo-Afshari mixed load-dependent network as a
normalizing constant: the closed subnetwork carries the effective-capacity rates
mu_i^eff(n) = 1/EC_i(n), the open classes contribute the separable prefactor
sum_i log E_i(0), and the mean measures come out of normalizing-constant ratios
rather than the population-lattice recursion of pfqn_mvaldmx.

The reference is pfqn_mvaldmx itself, which is exact for these models: the two must
agree on lG, on the throughputs and on the queue lengths of both the closed and the
open classes, on multiserver rates, genuinely nonlinear rates, think time, several
closed chains, and the degenerate no-closed-jobs case.
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_mvaldmx, pfqn_ncldmx

TOL = 1e-9


def _ms(c, n):
    return np.minimum(np.arange(1, n + 1), c).astype(float)


CASES = {
    'multiserver rates': (
        np.array([0.0, 0.2]),
        np.array([[1.0, 0.8], [0.6, 0.4]]),
        np.array([5.0, np.inf]),
        np.zeros(2),
        np.vstack([_ms(3, 6), _ms(2, 6)]),
    ),
    'single server': (
        np.array([0.0, 0.2]),
        np.array([[1.0, 0.8], [0.6, 0.4]]),
        np.array([5.0, np.inf]),
        np.zeros(2),
        np.ones((2, 6)),
    ),
    'two closed chains': (
        np.array([0.0, 0.0, 0.15]),
        np.array([[1.0, 0.5, 0.7], [0.4, 0.9, 0.3], [0.6, 0.2, 0.5]]),
        np.array([3.0, 4.0, np.inf]),
        np.zeros(3),
        np.vstack([_ms(2, 8), _ms(4, 8), np.ones(8)]),
    ),
    'nonlinear rates': (
        np.array([0.0, 0.1]),
        np.array([[1.0, 0.9], [0.7, 0.3]]),
        np.array([5.0, np.inf]),
        np.zeros(2),
        np.array([[1, 1.7, 2.2, 2.5, 2.6, 2.6], [1, 1.4, 1.4, 1.4, 1.4, 1.4]], dtype=float),
    ),
    'think time': (
        np.array([0.0, 0.2]),
        np.array([[1.0, 0.8], [0.6, 0.4]]),
        np.array([6.0, np.inf]),
        np.array([2.0, 0.0]),
        np.vstack([_ms(2, 7), np.ones(7)]),
    ),
    'no closed jobs': (
        np.array([0.0, 0.2]),
        np.array([[1.0, 0.8], [0.6, 0.4]]),
        np.array([0.0, np.inf]),
        np.zeros(2),
        np.ones((2, 2)),
    ),
}


@pytest.mark.parametrize('name', sorted(CASES))
def test_matches_mvaldmx(name):
    lam, D, N, Z, mu = CASES[name]
    M = D.shape[0]
    nc = pfqn_ncldmx(lam, D, N, Z, mu, np.ones(M))
    XN, QN, _, _, lGN, _ = pfqn_mvaldmx(lam, D, N, Z, mu, np.ones(M))
    assert nc.lG == pytest.approx(float(lGN), abs=TOL)
    np.testing.assert_allclose(nc.XN, np.asarray(XN).flatten(), atol=TOL)
    np.testing.assert_allclose(nc.QN, np.asarray(QN), atol=TOL)


def test_open_prefactor_is_the_load_independent_shrink():
    # With unit rates E_i(0) = 1/(1-rho_i), so lGopen = -sum_i log(1-rho_i).
    lam = np.array([0.0, 0.25])
    D = np.array([[1.0, 0.8], [0.5, 0.4]])
    N = np.array([6.0, np.inf])
    mu = np.ones((2, 7))
    nc = pfqn_ncldmx(lam, D, N, np.zeros(2), mu, np.ones(2))
    rho = lam @ D.T
    assert nc.lGopen == pytest.approx(-np.sum(np.log(1.0 - rho)), abs=1e-12)


def test_single_server_open_queue_length_is_the_classical_formula():
    # b_i = 1 leaves no marginal correction, so the open queue length collapses to
    # lambda_r D_ir (1 + Q_i^closed) / (1 - rho_i).
    lam = np.array([0.0, 0.2])
    D = np.array([[1.0, 0.8], [0.6, 0.4]])
    N = np.array([5.0, np.inf])
    mu = np.ones((2, 6))
    nc = pfqn_ncldmx(lam, D, N, np.zeros(2), mu, np.ones(2))
    rho = lam @ D.T
    for i in range(2):
        expected = lam[1] * D[i, 1] * (1.0 + nc.QN[i, 0]) / (1.0 - rho[i])
        assert nc.QN[i, 1] == pytest.approx(expected, rel=1e-10)


def test_arrival_rate_on_a_closed_class_is_rejected():
    with pytest.raises(ValueError):
        pfqn_ncldmx(np.array([0.3, 0.2]),
                    np.array([[1.0, 0.8], [0.6, 0.4]]),
                    np.array([5.0, np.inf]),
                    np.zeros(2), np.ones((2, 6)), np.ones(2))
