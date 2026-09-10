"""Neuse-Chandy SCAT (pfqn_scat) and its SolverMVA dispatch.

SCAT is Linearizer with ONE Delta refresh instead of three, so what pins it is
not a closed form but its POSITION: it must sit strictly between Bard-Schweitzer
(no refresh) and Linearizer (three), and equal neither. The absolute figures
below are MATLAB's own, printed by pfqn_scat.m and by SolverMVA with
options.method = 'scat' on the same two models the C++ method matrix uses, so a
drift in the refresh count moves them at once.
"""

import numpy as np

from line_solver import (ClosedClass, Delay, Exp, Network, Queue,
                         SchedStrategy, SolverMVA)
from line_solver.api.pfqn import pfqn_bs, pfqn_linearizer, pfqn_mva, pfqn_scat

SCAT_L = np.array([[1.0, 0.5], [0.4, 1.2], [0.8, 0.3]])
SCAT_N = np.array([4.0, 3.0])
SCAT_Z = np.array([1.0, 2.0])
# MATLAB pfqn_scat(L,N,Z,[PS;PS;PS],1e-8,1000)
SCAT_X = np.array([0.617697394, 0.431056566])


def test_pfqn_scat_matches_matlab_reference():
    Q, U, W, T, C, X, iters = pfqn_scat(SCAT_L, SCAT_N, SCAT_Z, ['PS'] * 3)
    assert np.allclose(np.asarray(X).flatten(), SCAT_X, atol=1e-9)


def test_pfqn_scat_sits_between_bard_schweitzer_and_linearizer():
    _, _, _, _, _, Xs, _ = pfqn_scat(SCAT_L, SCAT_N, SCAT_Z, ['PS'] * 3)
    _, _, _, _, _, Xl, _ = pfqn_linearizer(SCAT_L, SCAT_N, SCAT_Z, ['PS'] * 3)
    Xb = np.asarray(pfqn_bs(SCAT_L, SCAT_N, SCAT_Z)[0]).flatten()
    Xe = np.asarray(pfqn_mva(SCAT_L, SCAT_N, SCAT_Z)[0]).flatten()
    Xs = np.asarray(Xs).flatten()
    Xl = np.asarray(Xl).flatten()
    for r in range(2):
        e_scat = abs(Xs[r] - Xe[r])
        e_lin = abs(Xl[r] - Xe[r])
        e_bs = abs(Xb[r] - Xe[r])
        assert e_lin < e_scat, f"linearizer must beat scat on class {r}"
        assert e_scat < e_bs, f"scat must beat bard-schweitzer on class {r}"


def test_pfqn_scat_is_cheaper_than_linearizer():
    _, _, _, _, _, _, its = pfqn_scat(SCAT_L, SCAT_N, SCAT_Z, ['PS'] * 3)
    _, _, _, _, _, _, itl = pfqn_linearizer(SCAT_L, SCAT_N, SCAT_Z, ['PS'] * 3)
    # one refresh pass against three; the ratio is not exactly 1/3 because the
    # inner Core loop converges in a different number of steps at each pass
    assert its < itl


def _model_a():
    """Delay(Z=1) + FCFS Queue(D=0.5), one closed class of 3."""
    model = Network('A')
    d = Delay(model, 'Delay')
    q = Queue(model, 'Queue', SchedStrategy.FCFS)
    c = ClosedClass(model, 'C1', 3, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def _model_b():
    """The same with two servers at the queue and a population of 5."""
    model = Network('B')
    d = Delay(model, 'Delay')
    q = Queue(model, 'Queue', SchedStrategy.FCFS)
    q.setNumberOfServers(2)
    c = ClosedClass(model, 'C1', 5, d)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    model.link(Network.serialRouting(d, q))
    return model


def _metrics(model, method):
    solver = SolverMVA(model, method)
    res = solver.getAvg()
    Q, U, R = res[0], res[1], res[2]
    return np.array([Q[0, 0], Q[1, 0], U[1, 0], R[1, 0]])


def test_solver_mva_dispatches_scat():
    model = _model_a()
    expected = np.array([1.593017294, 1.406982706, 0.796508647, 0.883218727])
    assert np.allclose(_metrics(model, 'scat'), expected, atol=1e-8)
    # the advertised amva.* spelling must resolve to the same algorithm
    assert np.allclose(_metrics(model, 'amva.scat'), expected, atol=1e-8)
    methods = SolverMVA(model).listValidMethods()
    assert 'scat' in methods and 'amva.scat' in methods


def test_scat_takes_a_multiserver_model_through_seidmann():
    # unlike aql/qsa/tay, SCAT is not refused on a multiserver model: it reaches
    # the algorithm Seidmann-scaled, exactly as bs does
    model = _model_b()
    expected = np.array([2.883835883, 2.116164117, 0.720958971, 0.733801854])
    assert np.allclose(_metrics(model, 'scat'), expected, atol=1e-8)
    assert 'scat' in SolverMVA(model).listValidMethods()
