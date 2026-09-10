"""
Analytic validation of ctmc_fau.

The oracle is the matrix exponential pi0 expm(Q t), computed independently by
scipy, rather than a recorded baseline: a recorded value cannot catch an error
that the reference implementation shares. Two structural properties are
asserted alongside the numbers, because they are what the method promises and
what a plausible-looking wrong answer would break: the result is a
componentwise lower bound on the exact distribution, and its missing mass IS
its L1 error.
"""

import numpy as np
import pytest
import scipy.sparse as sp
from scipy.linalg import expm

from line_solver.api.mc import ctmc_fau, ctmc_uniformization


def mm1k_generator(lam, mu, K):
    """Generator of an M/M/1/K queue."""
    n = K + 1
    Q = sp.lil_matrix((n, n))
    for i in range(n):
        if i < n - 1:
            Q[i, i + 1] = lam
        if i > 0:
            Q[i, i - 1] = mu
    Q.setdiag(-np.asarray(Q.sum(axis=1)).ravel())
    return Q.tocsr()


def exact_transient(pi0, Q, t):
    """pi0 expm(Q t) through a dense matrix exponential."""
    dense = Q.toarray() if sp.issparse(Q) else np.asarray(Q, dtype=np.float64)
    return pi0 @ expm(dense * t)


def test_two_state_chain_matches_matrix_exponential():
    a, b, t = 3.0, 1.5, 0.7
    Q = np.array([[-a, a], [b, -b]])
    pi0 = np.array([1.0, 0.0])
    pit, info = ctmc_fau(pi0, Q, t)
    exact = exact_transient(pi0, Q, t)
    assert np.allclose(pit, exact, atol=1e-6)
    assert info.steps > 1


def test_error_is_the_missing_mass():
    Q = mm1k_generator(1.0, 2.0, 30)
    pi0 = np.zeros(31)
    pi0[0] = 1.0
    pit, info = ctmc_fau(pi0, Q, 5.0)
    exact = exact_transient(pi0, Q, 5.0)
    l1 = float(np.sum(np.abs(pit - exact)))
    # The three approximations only remove mass, so the defect IS the error.
    assert info.error_bound == pytest.approx(l1, rel=1e-6, abs=1e-15)
    assert info.error_bound <= 1e-5


def test_result_is_a_componentwise_lower_bound():
    Q = mm1k_generator(1.0, 2.0, 30)
    pi0 = np.zeros(31)
    pi0[0] = 1.0
    pit, _ = ctmc_fau(pi0, Q, 5.0)
    exact = exact_transient(pi0, Q, 5.0)
    assert np.all(pit <= exact + 1e-14)


def test_tolerance_tightens_the_error():
    Q = mm1k_generator(1.0, 2.0, 20)
    pi0 = np.zeros(21)
    pi0[0] = 1.0
    loose, info_loose = ctmc_fau(pi0, Q, 4.0, epsilon=1e-4)
    tight, info_tight = ctmc_fau(pi0, Q, 4.0, epsilon=1e-10)
    exact = exact_transient(pi0, Q, 4.0)
    assert info_tight.error_bound < info_loose.error_bound
    assert info_tight.steps > info_loose.steps
    assert np.sum(np.abs(tight - exact)) < np.sum(np.abs(loose - exact))


def test_fast_states_carrying_no_mass_do_not_set_the_cost():
    # Six slow states, then four states of rate 1e6 that the initial
    # distribution cannot reach within the horizon. Ordinary uniformization
    # pays for the fast ones, adaptive uniformization does not.
    ns, nf, t = 6, 4, 1.0
    n = ns + nf
    Q = np.zeros((n, n))
    for i in range(ns - 1):
        Q[i, i + 1] = 0.5
        Q[i + 1, i] = 0.4
    for i in range(ns, n - 1):
        Q[i, i + 1] = 1e6
        Q[i + 1, i] = 1e6
    Q[n - 1, ns] = 1e6
    Q -= np.diag(Q.sum(axis=1))
    pi0 = np.zeros(n)
    pi0[0] = 1.0

    pit, info = ctmc_fau(pi0, Q, t)
    exact = exact_transient(pi0, Q, t)
    assert np.allclose(pit, exact, atol=1e-6)
    assert info.lambda_max <= 1.0
    assert info.uniform_rate >= 1e6
    # The step count follows the visited rate, not the global one.
    assert info.steps < 50
    reference = ctmc_uniformization(pi0, Q, t)
    assert np.allclose(pit, reference, atol=1e-6)


def test_occupancy_threshold_bounds_the_support():
    # A wide chain in which the far states hold negligible mass at t.
    K, t = 400, 4.0
    Q = mm1k_generator(1.0, 3.0, K)
    pi0 = np.zeros(K + 1)
    pi0[0] = 1.0
    pit, info = ctmc_fau(pi0, Q, t, epsilon=1e-8, delta=1e-10)
    exact = exact_transient(pi0, Q, t)
    assert info.support_max < K + 1
    assert np.sum(np.abs(pit - exact)) <= info.error_bound + 1e-14
    assert np.sum(np.abs(pit - exact)) < 1e-7


def test_absorbing_chain_terminates_and_is_exact():
    Q = np.array([[-2.0, 2.0, 0.0], [0.0, -1.0, 1.0], [0.0, 0.0, 0.0]])
    pi0 = np.array([1.0, 0.0, 0.0])
    pit, info = ctmc_fau(pi0, Q, 3.0)
    exact = exact_transient(pi0, Q, 3.0)
    assert np.allclose(pit, exact, atol=1e-6)
    assert info.absorbed
    assert info.steps == 3


def test_zero_horizon_returns_the_initial_distribution():
    Q = mm1k_generator(1.0, 2.0, 5)
    pi0 = np.zeros(6)
    pi0[2] = 1.0
    pit, info = ctmc_fau(pi0, Q, 0.0)
    assert np.array_equal(pit, pi0)
    assert info.error_bound == 0.0


def test_sparse_and_dense_generators_agree():
    Q = mm1k_generator(1.0, 2.0, 15)
    pi0 = np.zeros(16)
    pi0[0] = 1.0
    sparse_pit, _ = ctmc_fau(pi0, Q, 3.0)
    dense_pit, _ = ctmc_fau(pi0, Q.toarray(), 3.0)
    assert np.allclose(sparse_pit, dense_pit, rtol=0, atol=1e-15)


def _mm1_model(cutoff, timespan, config):
    """M/M/1 at rho = 0.5 behind SolverCTMC, truncated at cutoff."""
    from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                             SchedStrategy, SolverCTMC)
    from line_solver.solvers.solver_ctmc.solver_ctmc import SolverCTMCOptions
    model = Network('mm1_transient')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(0.5))
    queue.setService(oclass, Exp(1.0))
    model.link(Network.serialRouting(source, queue, sink))
    options = SolverCTMCOptions(cutoff=cutoff, timespan=timespan, verbose=False,
                                config=config)
    return SolverCTMC(model, options)


def _queue_trajectory(solver):
    QNt, _, _ = solver.getTranAvg()
    tran = QNt[1][0]
    return np.asarray(tran.t, dtype=float), np.asarray(tran.metric, dtype=float)


def test_solver_transient_method_fau_matches_the_matrix_exponential():
    # The solver path: options.config['transient_method'] = 'fau' advances the
    # same forward equation by fast adaptive uniformization. The oracle is the
    # analyzer's own expm trajectory on the same generator and the same grid.
    t_ref, q_ref = _queue_trajectory(_mm1_model(6, [0, 10], {}))
    t_fau, q_fau = _queue_trajectory(
        _mm1_model(6, [0, 10], {'transient_method': 'fau'}))
    assert np.allclose(t_ref, t_fau)
    assert np.max(np.abs(q_fau - q_ref)) < 1e-6


def test_solver_transient_method_fau_tightens_with_epsilon():
    _, q_ref = _queue_trajectory(_mm1_model(6, [0, 10], {}))
    _, q_loose = _queue_trajectory(
        _mm1_model(6, [0, 10], {'transient_method': 'fau'}))
    _, q_tight = _queue_trajectory(
        _mm1_model(6, [0, 10], {'transient_method': 'fau', 'fau_epsilon': 1e-12}))
    assert np.max(np.abs(q_tight - q_ref)) < np.max(np.abs(q_loose - q_ref))
    assert np.max(np.abs(q_tight - q_ref)) < 1e-12


def test_solver_rejects_an_unknown_transient_method():
    with pytest.raises(ValueError):
        _queue_trajectory(_mm1_model(6, [0, 10], {'transient_method': 'uniformization'}))
