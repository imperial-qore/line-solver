"""
options.config['chain_aggregation'] on SolverCTMC.

ModelAdapter.aggregate_chains collapses every chain onto a single class and
sn_deaggregate_chain_results maps chain-level metrics back through alpha. Both
existed in all four codebases with NO solver consumer at all: the transform was
exercised by examples and tests only, so nothing in the solver stack depended on
it and a defect in it could not surface as a wrong answer.

The oracle is an identity rather than a golden. On a PRODUCT-FORM model the
chain is the unit MVA and convolution already solve in, so the aggregation is
exact and the aggregated solve must reproduce the exact multiclass CTMC table,
station by station and class by class. That is a statement about the transform
and the deaggregation together, and it is what fails if either mis-derives
alpha.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Exp, Network, Queue,
                         SchedStrategy, SolverCTMC)


def two_class_one_chain(n=3):
    """Delay -> Q1 -> Q2 -> Delay, class switch on the Q1 -> Q2 link."""
    m = Network('agg')
    d = Delay(m, 'Think')
    q1 = Queue(m, 'Q1', SchedStrategy.PS)
    q2 = Queue(m, 'Q2', SchedStrategy.PS)
    c1 = ClosedClass(m, 'C1', n, d, 0)
    c2 = ClosedClass(m, 'C2', 0, d, 0)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q1.setService(c1, Exp(2.0))
    q1.setService(c2, Exp(2.0))
    q2.setService(c1, Exp(3.0))
    q2.setService(c2, Exp(3.0))
    P = m.initRoutingMatrix()
    P.set(c1, c1, d, q1, 1.0)
    P.set(c1, c2, q1, q2, 1.0)
    P.set(c2, c1, q2, d, 1.0)
    m.link(P)
    return m


def one_class():
    m = Network('plain')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.PS)
    c = ClosedClass(m, 'C1', 3, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    return m


def test_the_model_really_has_two_classes_in_one_chain():
    sn = two_class_one_chain().getStruct()
    assert int(sn.nclasses) == 2
    assert int(sn.nchains) == 1


def test_aggregated_solve_reproduces_the_exact_table():
    exact = SolverCTMC(two_class_one_chain())
    agg = SolverCTMC(two_class_one_chain(), config={'chain_aggregation': True})
    for name in ('getAvgQLen', 'getAvgUtil', 'getAvgTput'):
        a = np.asarray(getattr(exact, name)())
        b = np.asarray(getattr(agg, name)())
        assert a.shape == b.shape, name
        assert np.allclose(a, b, atol=1e-9), name


def test_flow_is_conserved_through_the_deaggregation():
    agg = SolverCTMC(two_class_one_chain(), config={'chain_aggregation': True})
    TN = np.asarray(agg.getAvgTput())
    per_station = TN.sum(axis=1)
    assert per_station[0] > 0
    assert np.allclose(per_station, per_station[0], atol=1e-9)


def test_population_is_conserved_through_the_deaggregation():
    agg = SolverCTMC(two_class_one_chain(5), config={'chain_aggregation': True})
    QN = np.asarray(agg.getAvgQLen())
    assert QN.sum() == pytest.approx(5.0, abs=1e-9)


def test_the_option_is_declined_where_the_transform_is_the_identity():
    # nchains == nclasses, so the guard sends the model down the ordinary path.
    exact = np.asarray(SolverCTMC(one_class()).getAvgQLen())
    agg = SolverCTMC(one_class(), config={'chain_aggregation': True})
    assert np.allclose(np.asarray(agg.getAvgQLen()), exact, atol=1e-12)
    assert 'chainaggr' not in str(getattr(agg._result, 'method', ''))


def test_the_aggregated_solve_names_itself():
    agg = SolverCTMC(two_class_one_chain(), config={'chain_aggregation': True})
    agg.getAvgQLen()
    assert 'chainaggr' in str(agg._result.method)
