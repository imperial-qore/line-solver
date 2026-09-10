"""
options.config['fes_stations'] on SolverCTMC.

Flow-equivalent aggregation collapses a station subset into one load-dependent
station whose rate is the isolated subnetwork's throughput. The transform has
existed in all four codebases with NO solver consumer at all: it was exercised
by examples and tests only, so nothing in the solver stack depended on it. This
is that consumer, and the collapsed stations' own metrics come back through the
Chandy-Herzog-Woo conditional sum E[Q_i] = sum_n P(N_fes = n) * Q_i(n).

The oracle is an identity, not a golden. On a product-form model the
decomposition is EXACT, so the reduced solve plus the conditioning must
reproduce the full chain's table station by station, INCLUDING the collapsed
stations. A wrong FES rate moves the surviving stations, a wrong conditional sum
moves only the collapsed ones, and a wrong visit ratio moves only their
throughput.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue,
                         SchedStrategy, SolverCTMC)


def cycle(n=4):
    """Think -> Q1 -> Q2 -> Q3 -> Think, a closed product-form cycle."""
    m = Network('fes')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    q3 = Queue(m, 'Q3', SchedStrategy.PS)
    c = ClosedClass(m, 'C1', n, d, 0)
    d.setService(c, Exp(1.0))
    q1.setService(c, Exp(2.0))
    q2.setService(c, Exp(3.0))
    q3.setService(c, Exp(4.0))
    m.link(Network.serialRouting(d, q1, q2, q3))
    return m


def solve(fes=None):
    cfg = {'fes_stations': fes} if fes else {}
    return SolverCTMC(cycle(), config=cfg)


def test_reduced_solve_reproduces_the_exact_table():
    exact = solve()
    fes = solve([2, 3])
    for name in ('getAvgQLen', 'getAvgUtil', 'getAvgTput'):
        a = np.asarray(getattr(exact, name)())
        b = np.asarray(getattr(fes, name)())
        assert a.shape == b.shape, name
        assert np.allclose(a, b, atol=1e-9), (name, a, b)


def test_population_is_conserved_by_the_conditional_split():
    fes = solve([2, 3])
    QN = np.asarray(fes.getAvgQLen())
    assert QN.sum() == pytest.approx(4.0, abs=1e-9)


def test_a_different_subset_gives_the_same_answer():
    # The choice of subset is the caller's and changes only which stations are
    # enumerated, so an exact decomposition must be invariant to it.
    exact = np.asarray(solve().getAvgQLen())
    for subset in ([1, 2], [2, 3]):
        got = np.asarray(solve(subset).getAvgQLen())
        assert np.allclose(got, exact, atol=1e-9), subset


def test_the_reduction_names_itself():
    fes = solve([2, 3])
    fes.getAvgQLen()
    assert '/fes' in str(fes._result.method)


def test_a_subset_that_saves_nothing_is_refused_by_name():
    with pytest.raises(ValueError, match='at least two stations'):
        solve([2]).getAvgQLen()
    with pytest.raises(ValueError, match='every station'):
        solve([0, 1, 2, 3]).getAvgQLen()
    with pytest.raises(ValueError, match='0-based station'):
        solve([2, 99]).getAvgQLen()
