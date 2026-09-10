"""``pfqn_qdlin``: the array-level twin of SolverMVA's ``qdlin`` method.

The kernel exists to answer, on a plain demand matrix, exactly what the solver
answers on the equivalent ``Network``, so the sharpest test is the solver
itself: every case below builds the model whose chain demands are ``L`` and
whose think times are ``Z`` and asserts that the two agree to machine
precision, iteration path included.

Two of the checks pin PROPERTIES rather than numbers, and those are the ones a
bad transcription trips:

- With ONE class the correction the reference applies coincides with the
  queue-dependent AMVA form, because its class-aggregate gamma lands in the one
  slice that is read. With two or more classes it does not, and ``qdlin`` then
  parts company with ``lin``. Both halves are asserted.
- ``mu`` and ``nservers`` are DIFFERENT mechanisms here, unlike in
  ``pfqn_qdamva``, which folds the multiserver curve into ``mu``.

``pfqn_qdlin`` used to be the Wang-Sevcik QDLIN, which scaled the own-class
contribution by (N_r - 1)/N_r and was therefore Bard-Schweitzer under another
name. That it is no longer ``pfqn_bs`` is asserted too.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, GlobalConstants, MVA, Network,
                         Queue, SchedStrategy, VerboseLevel)
from line_solver.api.pfqn import pfqn_bs, pfqn_qdlin

GlobalConstants.set_verbose(VerboseLevel.SILENT)


def _build(L, N, Z, nservers=None, mu=None):
    """The Network whose chain demands are L and whose think times are Z.

    The solver recovers a demand as 1/rate, so the caller must hand the kernel
    the same doubles it hands the model: see ``_demands``.
    """
    M, K = L.shape
    model = Network('qdlin')
    nodes = []
    delay = None
    if np.any(np.asarray(Z) > 0):
        delay = Delay(model, 'Think')
        nodes.append(delay)
    queues = []
    for k in range(M):
        q = Queue(model, 'Q%d' % (k + 1), SchedStrategy.PS)
        if nservers is not None and nservers[k] > 1:
            q.setNumberOfServers(int(nservers[k]))
        if mu is not None:
            q.set_load_dependence(list(mu[k, :]))
        queues.append(q)
        nodes.append(q)
    for r in range(K):
        cls = ClosedClass(model, 'C%d' % (r + 1), int(N[r]), nodes[0], 0)
        if delay is not None:
            delay.set_service(cls, Exp(1.0 / Z[r]))
        for k in range(M):
            queues[k].set_service(cls, Exp(1.0 / L[k, r]))
    model.link(Network.serialRouting(*nodes))
    return model, (1 if delay is not None else 0)


def _demands(L, Z):
    """The demands as the solver will see them, 1/(1/x)."""
    Ze = 1.0 / (1.0 / Z) if np.any(np.asarray(Z) > 0) else np.asarray(Z, dtype=float)
    return 1.0 / (1.0 / L), Ze


L1 = np.array([[0.5], [0.3]])
N1 = np.array([3.0])
Z1 = np.array([1.0])

L2 = np.array([[0.4, 0.9], [0.6, 0.1], [0.2, 0.5]])
N2 = np.array([4.0, 3.0])
Z2 = np.array([2.0, 1.0])


@pytest.mark.parametrize('L,N,Z,nservers,mu', [
    (L1, N1, Z1, None, None),
    (L2, N2, Z2, None, None),
    (L1, N1, Z1, np.array([2.0, 1.0]), None),
    (L2, N2, Z2, np.array([3.0, 1.0, 2.0]), None),
    (L2, N2, np.zeros(2), None, None),
    (L1, N1, Z1, None, np.array([[1.0, 2.0, 2.0, 2.0], [1.0, 1.0, 1.0, 1.0]])),
])
def test_matches_solvermva_qdlin(L, N, Z, nservers, mu):
    model, off = _build(L, N, Z, nservers, mu)
    Le, Ze = _demands(L, Z)
    Q, U, R, X, C, it = pfqn_qdlin(Le, N, Ze, mu, nservers)

    M, K = L.shape
    solver = MVA(model, 'qdlin')
    Qs = np.asarray(solver.getAvgQLen(), dtype=float).reshape(M + off, K)[off:, :]
    Us = np.asarray(solver.getAvgUtil(), dtype=float).reshape(M + off, K)[off:, :]
    Rs = np.asarray(solver.getAvgRespT(), dtype=float).reshape(M + off, K)[off:, :]
    Xs = np.asarray(solver.getAvgTput(), dtype=float).reshape(M + off, K)[off:, :]

    np.testing.assert_allclose(Q, Qs, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(U, Us, rtol=1e-12, atol=1e-14)
    np.testing.assert_allclose(R, Rs, rtol=1e-12, atol=1e-14)
    for r in range(K):
        np.testing.assert_allclose(np.full(M, X[0, r]), Xs[:, r], rtol=1e-12, atol=1e-14)
    assert it > 0


def test_littles_law_holds():
    Le, Ze = _demands(L2, Z2)
    Q, U, R, X, C, it = pfqn_qdlin(Le, N2, Ze)
    for r in range(N2.size):
        assert Q[:, r].sum() + X[0, r] * Z2[r] == pytest.approx(N2[r], rel=1e-6)
        assert C[0, r] == pytest.approx(R[:, r].sum() + Z2[r], rel=1e-6)


def test_one_class_coincides_with_the_aggregate_form_and_two_do_not():
    # The class-aggregate gamma lands in slice 0, the only slice a single-chain
    # model reads, so qdlin and lin are the SAME computation there. With two
    # chains the untouched slices are read as zero and they part company.
    model1, off1 = _build(L1, N1, Z1, np.array([2.0, 1.0]), None)
    q1 = np.asarray(MVA(model1, 'qdlin').getAvgQLen(), dtype=float)
    l1 = np.asarray(MVA(model1, 'lin').getAvgQLen(), dtype=float)
    np.testing.assert_allclose(q1, l1, rtol=0, atol=0)

    model2, off2 = _build(L2, N2, Z2, np.array([3.0, 1.0, 2.0]), None)
    q2 = np.asarray(MVA(model2, 'qdlin').getAvgQLen(), dtype=float)
    l2 = np.asarray(MVA(model2, 'lin').getAvgQLen(), dtype=float)
    assert np.max(np.abs(q2 - l2)) > 1e-6


def test_mu_and_nservers_are_different_mechanisms():
    Le, Ze = _demands(L1, Z1)
    mu = np.array([[1.0, 2.0, 2.0, 2.0], [1.0, 1.0, 1.0, 1.0]])
    via_mu = pfqn_qdlin(Le, N1, Ze, mu, None)[0]
    via_servers = pfqn_qdlin(Le, N1, Ze, None, np.array([2.0, 1.0]))[0]
    # Close, because both describe the same two-server station, but not equal:
    # the server count goes through the softmin, the lattice through interp.
    assert not np.array_equal(via_mu, via_servers)
    assert np.max(np.abs(via_mu - via_servers)) < 1e-4


def test_is_no_longer_bard_schweitzer():
    # The Wang-Sevcik QDLIN this name used to carry reproduced pfqn_bs to
    # iteration tolerance. The Linearizer arm does not.
    Le, Ze = _demands(L2, Z2)
    Q = pfqn_qdlin(Le, N2, Ze)[0]
    Qbs = np.asarray(pfqn_bs(Le, N2, Ze)[1], dtype=float).reshape(Q.shape)
    assert np.max(np.abs(Q - Qbs)) > 1e-3


def test_a_class_with_no_jobs_reports_zero_throughput():
    # X carries a 1/sum(ST) seed for every class; the sweeps only write the
    # classes with jobs, so an uncleared seed would be reported as the empty
    # class's throughput. Caught by the JAR twin's test, present in three ports.
    L = np.array([[0.5, 0.2], [0.3, 0.7]])
    Q, U, R, X, C, it = pfqn_qdlin(L, np.array([2.0, 0.0]), np.zeros(2))
    assert X[0, 1] == 0.0
    assert np.all(Q[:, 1] == 0.0)
    assert np.all(U[:, 1] == 0.0)
    assert C[0, 1] == 0.0
    assert X[0, 0] > 0.0


def test_empty_population_and_refusals():
    Q, U, R, X, C, it = pfqn_qdlin(L1, np.zeros(1), Z1)
    assert it == 0
    assert np.all(Q == 0.0)

    with pytest.raises(ValueError):
        pfqn_qdlin(L1, np.array([1.0, 2.0]), Z1)
    with pytest.raises(ValueError):
        pfqn_qdlin(L1, N1, np.array([1.0, 2.0]))
    with pytest.raises(ValueError):
        pfqn_qdlin(L1, N1, Z1, None, np.array([1.0, 1.0, 1.0]))
    with pytest.raises(ValueError):
        pfqn_qdlin(L1, np.array([np.inf]), Z1)
