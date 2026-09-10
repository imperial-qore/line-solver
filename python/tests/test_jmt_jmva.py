"""SolverJMT analytical (JMVA) methods.

The eleven jmva* aliases had no test in any codebase, which is how the native
Python path came to advertise them, treat them as analytical in
isStochasticMethod, and then run the JSIM SIMULATION for every one of them.
These tests pin the three properties that separate the two engines: an
analytical solve is exact on a product-form model, it does not move with the
seed, and it reports the method that was asked for.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClassSwitch, Source, Sink,
                         ClosedClass, OpenClass, Exp, SchedStrategy, SolverJMT,
                         SolverMVA)


def _two_chain_model():
    """Closed product-form model, two chains of one class each."""
    m = Network('jmva_two_chain')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', 3, d)
    c2 = ClosedClass(m, 'C2', 2, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(0.5))
    q.setService(c1, Exp(2.0))
    q.setService(c2, Exp(3.0))
    P = m.initRoutingMatrix()
    for c in (c1, c2):
        P[c] = Network.serialRouting(d, q)
    m.link(P)
    return m


def _one_chain_two_class_model():
    """Closed model whose two classes share ONE chain (class switch).

    This is the model that exercises the chain-to-class disaggregation: JMVA
    solves a single aggregated chain and the per-class split has to be
    reconstructed from alpha and the chain service times.
    """
    m = Network('jmva_one_chain')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.PS)
    cs = ClassSwitch(m, 'CS', [[0, 1], [1, 0]])
    a = ClosedClass(m, 'A', 4, d)
    b = ClosedClass(m, 'B', 0, d)
    d.setService(a, Exp(1.0))
    d.setService(b, Exp(2.0))
    q.setService(a, Exp(3.0))
    q.setService(b, Exp(1.5))
    P = m.initRoutingMatrix()
    for c in (a, b):
        P[c] = Network.serialRouting(d, q, cs)
    m.link(P)
    return m


def _multiserver_model():
    m = Network('jmva_multiserver')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.PS)
    q.setNumberOfServers(3)
    c = ClosedClass(m, 'C1', 6, d)
    d.setService(c, Exp(0.5))
    q.setService(c, Exp(1.5))
    m.link(Network.serialRouting(d, q))
    return m


def _open_model():
    m = Network('jmva_open')
    s = Source(m, 'Src')
    q = Queue(m, 'Q', SchedStrategy.PS)
    k = Sink(m, 'Snk')
    c = OpenClass(m, 'C1')
    s.setArrival(c, Exp(0.8))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def _solve(model, method, seed=23000):
    solver = SolverJMT(model, method, seed=seed, samples=20000, verbose=0)
    solver._table_silent = True
    return solver


@pytest.mark.parametrize('method', ['jmva', 'jmva.mva'])
def test_jmva_is_exact_on_product_form(method):
    """Exact MVA under either name agrees with LINE's own exact MVA."""
    exact = SolverMVA(_two_chain_model(), 'exact')
    Q_exact = np.asarray(exact.getAvgQLen())
    T_exact = np.asarray(exact.getAvgTput())

    solver = _solve(_two_chain_model(), method)
    assert np.allclose(np.asarray(solver.getAvgQLen()), Q_exact, atol=1e-9)
    assert np.allclose(np.asarray(solver.getAvgTput()), T_exact, atol=1e-9)


def test_jmva_disaggregates_one_chain_into_two_classes():
    """A chain holding two classes is split back onto the classes.

    Both classes must carry a share of the chain, and the closed population
    must be conserved: a disaggregation that dropped the split would leave one
    class holding everything.
    """
    exact = SolverMVA(_one_chain_two_class_model(), 'exact')
    Q_exact = np.asarray(exact.getAvgQLen())

    solver = _solve(_one_chain_two_class_model(), 'jmva')
    Q = np.asarray(solver.getAvgQLen())

    assert np.allclose(Q, Q_exact, atol=1e-6)
    assert Q.sum() == pytest.approx(4.0, abs=1e-6)
    assert (Q.sum(axis=0) > 0).all()


def test_jmva_multiserver_matches_exact_mva():
    exact = SolverMVA(_multiserver_model(), 'exact')
    Q_exact = np.asarray(exact.getAvgQLen())
    solver = _solve(_multiserver_model(), 'jmva')
    assert np.allclose(np.asarray(solver.getAvgQLen()), Q_exact, atol=1e-6)


def test_jmva_open_model_matches_exact_mva():
    """An open model reaches JMVA without its Source, which carries no queue."""
    exact = SolverMVA(_open_model(), 'exact')
    Q_exact = np.asarray(exact.getAvgQLen())
    solver = _solve(_open_model(), 'jmva')
    assert np.allclose(np.asarray(solver.getAvgQLen()), Q_exact, atol=1e-6)


def test_jmva_is_seed_independent():
    """An analytical solve must not move with the seed; a simulation does."""
    a = np.asarray(_solve(_two_chain_model(), 'jmva', seed=1).getAvgQLen())
    b = np.asarray(_solve(_two_chain_model(), 'jmva', seed=999999).getAvgQLen())
    assert np.array_equal(a, b)


def test_jsim_is_seed_dependent():
    """The counterpart: the simulation DOES move with the seed.

    Without this, test_jmva_is_seed_independent would also pass if jmva were
    served by a simulation that ignored the seed.
    """
    a = np.asarray(_solve(_two_chain_model(), 'jsim', seed=1).getAvgQLen())
    b = np.asarray(_solve(_two_chain_model(), 'jsim', seed=999999).getAvgQLen())
    assert not np.array_equal(a, b)


@pytest.mark.parametrize('method', ['jmva', 'jmva.mva', 'jmva.amva'])
def test_jmva_reports_the_method_it_ran(method):
    solver = _solve(_two_chain_model(), method)
    solver.getAvgQLen()
    assert solver.getMethod() == method


def test_amva_is_approximate_and_differs_from_exact():
    """Bard-Schweitzer is an approximation, so it must NOT reproduce exact MVA.

    This is the guard that the algorithm method name reaches JMVA at all: if every
    jmva* alias resolved to the same engine, this would fail.
    """
    exact = np.asarray(SolverMVA(_two_chain_model(), 'exact').getAvgQLen())
    amva = np.asarray(_solve(_two_chain_model(), 'jmva.amva').getAvgQLen())
    assert not np.allclose(amva, exact, atol=1e-6)
    assert np.allclose(amva, exact, atol=0.5)


def test_approximate_algorithms_refuse_multiserver():
    """MATLAB writeJMVA.m rejects these on a multi-server model; so must this."""
    with pytest.raises(Exception, match='multi-server'):
        _solve(_multiserver_model(), 'jmva.amva').getAvgQLen()


def test_norm_const_is_refused_for_simulation_and_returned_for_jmva():
    with pytest.raises(NotImplementedError):
        _solve(_two_chain_model(), 'jsim').getProbNormConstAggr()
    solver = _solve(_two_chain_model(), 'jmva')
    solver.getAvgQLen()
    # JMVA reports <normconst logValue>, which is NaN for algorithms that do
    # not compute one; the call must resolve either way rather than raise.
    float(solver.getProbNormConstAggr())
