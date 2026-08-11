"""
Validates the DAC (Distribution Analysis by Chain) algorithm against the
numerical results of E. de Souza e Silva, "Distribution Analysis of Product Form
Queueing Networks", UCLA Computer Science Department, CSD-870023, April 1987.

Reference model (Section 3, Figure 1): an availability model of a computer
system with one CPU and three memory modules (two active, one spare), served by
a single repairman. Chain 1 is the CPU, chain 2 the memory modules. Center 1
models CPU failures, center 2 memory failures and center 3 the repair facility.
Center 2 is queue dependent, since at most two modules can fail at a time.
Appendix parameters: L = [5,0; 0,10; 2,1], N = [1,3].
"""

import numpy as np
import pytest

from line_solver.api.pfqn import pfqn_dac, pfqn_mvald

TOL = 1e-4


def _model():
    L = np.array([[5.0, 0.0], [0.0, 10.0], [2.0, 1.0]])
    N = np.array([1, 3])
    return L, N


def _cold_spare_mu():
    # The spare cannot fail while it is not in use, so mu_2(n>=2) = 2.
    return np.array([[1.0, 1.0, 1.0, 1.0],
                     [1.0, 2.0, 2.0, 2.0],
                     [1.0, 1.0, 1.0, 1.0]])


def _hot_standby_mu():
    # The spare is powered and can also fail, so mu_2(3) grows accordingly.
    mu = _cold_spare_mu()
    mu[1, 2] = 2.5
    mu[1, 3] = 2.5
    return mu


def _at(P, states, a, b, c):
    row = np.where((states[:, 0] == a) & (states[:, 1] == b) & (states[:, 2] == c))[0]
    assert len(row) == 1, "state (%d,%d,%d) not enumerated" % (a, b, c)
    return P[row[0]]


def test_cold_spare_joint_distribution():
    L, N = _model()
    P, states, _, _, _, _, pi = pfqn_dac(L, N, None, _cold_spare_mu())

    # Joint distribution, Appendix step k=4.
    assert _at(P, states, 1, 3, 0) == pytest.approx(0.5381, abs=TOL)
    assert _at(P, states, 1, 2, 1) == pytest.approx(0.1076, abs=TOL)
    assert _at(P, states, 1, 1, 2) == pytest.approx(0.02152, abs=TOL)
    assert _at(P, states, 1, 0, 3) == pytest.approx(0.002152, abs=TOL)
    assert _at(P, states, 0, 3, 1) == pytest.approx(0.2153, abs=TOL)
    assert _at(P, states, 0, 2, 2) == pytest.approx(0.08609, abs=TOL)
    assert _at(P, states, 0, 1, 3) == pytest.approx(0.02582, abs=TOL)
    assert _at(P, states, 0, 0, 4) == pytest.approx(0.003442, abs=TOL)

    # States unreachable because a chain does not visit the center.
    assert _at(P, states, 4, 0, 0) == pytest.approx(0.0, abs=TOL)
    assert _at(P, states, 2, 2, 0) == pytest.approx(0.0, abs=TOL)
    assert _at(P, states, 0, 4, 0) == pytest.approx(0.0, abs=TOL)

    assert P.sum() == pytest.approx(1.0, abs=1e-10)

    # Probability that the CPU is working, Section 5 example.
    assert pi[0, 1] == pytest.approx(0.6694, abs=TOL)

    # Availability: at least one CPU and one memory module operational.
    av = (_at(P, states, 1, 3, 0) + _at(P, states, 1, 2, 1) + _at(P, states, 1, 1, 2))
    assert av == pytest.approx(0.6672, abs=TOL)


def test_cold_spare_mean_measures():
    L, N = _model()
    _, _, X, Q, _, _, _ = pfqn_dac(L, N, None, _cold_spare_mu())

    # Appendix step k=4 reports the per-customer queue lengths of the memory
    # chain as L(2)=0.8983 and L(3)=0.1017. The chain holds 3 exchangeable
    # customers, so compare per customer, at the precision the paper prints.
    assert Q[1, 1] / 3.0 == pytest.approx(0.8983, abs=TOL)
    assert Q[2, 1] / 3.0 == pytest.approx(0.1017, abs=TOL)
    assert Q[0, 1] == pytest.approx(0.0, abs=TOL)

    # The memory chain never visits the CPU failure center, and the CPU chain
    # holds a single customer split between its own center and the repairman.
    assert Q[0, 0] + Q[2, 0] == pytest.approx(1.0, abs=1e-10)


def test_hot_standby_joint_distribution():
    L, N = _model()
    P, states, _, _, _, _, _ = pfqn_dac(L, N, None, _hot_standby_mu())

    assert _at(P, states, 1, 3, 0) == pytest.approx(0.5068, abs=TOL)
    assert _at(P, states, 1, 2, 1) == pytest.approx(0.1267, abs=TOL)
    assert _at(P, states, 1, 1, 2) == pytest.approx(0.02535, abs=TOL)
    assert _at(P, states, 1, 0, 3) == pytest.approx(0.002536, abs=TOL)
    assert _at(P, states, 0, 3, 1) == pytest.approx(0.2027, abs=TOL)
    assert _at(P, states, 0, 2, 2) == pytest.approx(0.1014, abs=TOL)
    assert _at(P, states, 0, 1, 3) == pytest.approx(0.03041, abs=TOL)
    assert _at(P, states, 0, 0, 4) == pytest.approx(0.004054, abs=TOL)
    assert P.sum() == pytest.approx(1.0, abs=1e-10)


def test_agrees_with_mvald_on_load_dependent_networks():
    """DAC must reproduce the mean measures of load-dependent MVA exactly."""
    rng = np.random.default_rng(7)
    for _ in range(5):
        L = rng.uniform(0.5, 4.0, (3, 2))
        N = np.array([2, 3])
        K = int(N.sum())
        mu = np.ones((3, K))
        mu[0, :] = np.minimum(np.arange(1, K + 1), 2)  # two-server station
        mu[2, :] = np.arange(1, K + 1)                 # infinite server
        P, _, X, Q, U, _, pi = pfqn_dac(L, N, None, mu)
        Xm, Qm, _, _, _, _, pim = pfqn_mvald(L, N, np.zeros(2), mu)

        assert P.sum() == pytest.approx(1.0, abs=1e-10)
        assert np.allclose(X, np.asarray(Xm).flatten(), atol=1e-8)
        assert np.allclose(Q, np.asarray(Qm), atol=1e-8)
        assert np.allclose(pi, np.asarray(pim), atol=1e-8)


def test_think_time_appends_delay_station():
    """A non-zero think time is modelled as an appended infinite-server center."""
    L = np.array([[1.0, 2.0], [3.0, 1.0]])
    N = np.array([2, 2])
    Z = np.array([4.0, 5.0])
    P, states, X, Q, _, _, _ = pfqn_dac(L, N, Z)
    assert states.shape[1] == 3
    assert P.sum() == pytest.approx(1.0, abs=1e-10)

    # Equivalent model with the delay written out as an explicit IS station.
    Lx = np.array([[1.0, 2.0], [3.0, 1.0], [4.0, 5.0]])
    mux = np.vstack((np.ones((2, 4)), np.arange(1.0, 5.0)))
    Xm, Qm, _, _, _, _, _ = pfqn_mvald(Lx, N, np.zeros(2), mux)
    assert np.allclose(X, np.asarray(Xm).flatten(), atol=1e-8)
    assert np.allclose(Q, np.asarray(Qm)[:2, :], atol=1e-8)


def test_probabilities_are_non_negative():
    """Every state carries non-negative mass: the recursion is numerically stable."""
    L = np.array([[1.0], [1e-3]])
    N = np.array([8])
    P, _, _, _, _, _, _ = pfqn_dac(L, N)
    assert (P >= 0).all()
    assert P.sum() == pytest.approx(1.0, abs=1e-10)
