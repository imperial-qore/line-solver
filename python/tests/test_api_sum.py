"""
Tests for the summation method (SUM/ESUM) and the closing method,
validated against the worked examples in Bolch, Greiner, de Meer,
Trivedi, "Queueing Networks and Markov Chains", 2nd ed., Wiley, 2006.
"""

import numpy as np

from line_solver.api.sum import sum_closed, sum_closing


def test_example_9_5_product_form_sum():
    # Bolch Example 9.5: closed PF network, N=4 nodes, K=3, node 1 is -/M/2
    XN, QN, _, _, _ = sum_closed([0.5, 0.3, 0.4, 1.0], 3, 0,
                                 [2, 1, 1, np.inf], None, 1e-3)
    assert abs(XN[0] - 1.193) < 1e-3
    assert np.allclose(QN.ravel(), [0.637, 0.470, 0.700, 1.193], atol=1e-3)


def test_example_10_10_esum():
    # Bolch Example 10.10: closed NPF network, N=5 nodes, K=17 (Table 10.15)
    XN, _, _, _, _ = sum_closed(
        [1 / 13.5, 0.2 / 1.15, 0.4 / 1.2, 0.3 / 1.2, 0.1 / 1.7], 17, 0,
        [1, 3, np.inf, 4, 1], [1.0, 0.8, 1.0, 3.0, 1.6])
    assert abs(XN[0] - 11.34) < 5e-3
    assert abs(17 / XN[0] - 1.50) < 5e-3


def test_example_10_13_closing_mixed():
    # Bolch Example 10.13: mixed NPF network, class 1 closed (K1=9),
    # class 2 open (lambda=5, ca2=0.7), closed with K2=500 (Table 10.21)
    e1 = np.array([1, 1.428571, 1.428571, 0.428571, 0.428571])
    e2 = np.array([1, 1.25, 1.428571, 0.25, 0.428571])
    mu = np.array([4, 4, 6, 4, 5.0])
    c2 = np.array([0.4, 0.3, 0.3, 0.4, 0.5])
    L = np.column_stack([e1 / mu, e2 / mu])
    XN, QN, UN, _, _, _ = sum_closing(
        [0, 5], [1, 0.7], L, [3, 4, 3, 2, 2], np.column_stack([c2, c2]),
        [9, np.inf], [0, 0], 500)
    assert abs(XN[1] - 5.0) < 1e-2
    assert np.allclose(UN.sum(axis=1), [0.84, 0.85, 0.80, 0.43, 0.43], atol=5e-3)
    assert np.allclose(QN.sum(axis=1), [5.2, 5.8, 4.1, 1.0, 1.0], atol=5e-2)


def test_closing_converges_to_arrival_rate():
    # closing method: open-class throughput approaches lambda0 from below
    e = np.array([1, 3 / 7, 7 / 9, 1]) / 0.9
    L = e / np.array([9, 10, 12, 4])
    XN, _, _, _, _, _ = sum_closing(3, 1.5, L, [1, 1, 1, 1],
                                    [0.5, 0.8, 2.4, 4.0], None, None, 5000)
    assert abs(XN[0] - 3.0) < 1e-3


def test_solver_mva_method_sum():
    # SolverMVA(method='sum') wiring; reference values from MATLAB solver_mva_sum
    from line_solver import (Network, Delay, Queue, Source, Sink, ClosedClass,
                             OpenClass, SchedStrategy, Exp, Erlang, HyperExp,
                             SolverMVA)
    model = Network('sum_closed')
    think = Delay(model, 'Think')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    q2.setNumberOfServers(2)
    c = ClosedClass(model, 'C', 10, think, 0)
    think.setService(c, Exp(1))
    q1.setService(c, Erlang.fitMeanAndSCV(0.3, 0.5))
    q2.setService(c, HyperExp.fitMeanAndSCV(0.8, 3))
    model.link(Network.serialRouting(think, q1, q2))
    t = SolverMVA(model, 'sum').getAvgTable()
    assert abs(t.QLen[0] - 2.2234) < 1e-3
    assert abs(t.QLen[1] - 1.5257) < 1e-3
    assert abs(t.QLen[2] - 6.2509) < 1e-3
    assert abs(t.Util[2] - 0.88934) < 1e-3

    model = Network('sum_mixed')
    src = Source(model, 'Src')
    th = Delay(model, 'Th')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    snk = Sink(model, 'Snk')
    cc = ClosedClass(model, 'C', 5, th, 0)
    oc = OpenClass(model, 'O')
    src.setArrival(oc, Exp(0.4))
    th.setService(cc, Exp(1))
    q1.setService(cc, Erlang.fitMeanAndSCV(0.4, 0.5))
    q1.setService(oc, Erlang.fitMeanAndSCV(0.5, 0.5))
    P = model.initRoutingMatrix()
    P.set(cc, Network.serialRouting(th, q1))
    P.set(oc, Network.serialRouting(src, q1, snk))
    model.link(P)
    t = SolverMVA(model, 'sum').getAvgTable()
    rows = {(t.Station[i], t.JobClass[i]): i for i in range(len(t.Station))}
    assert abs(t.QLen[rows[('Th', 'C')]] - 1.622) < 1e-3
    assert abs(t.QLen[rows[('Q1', 'C')]] - 3.378) < 1e-3
    assert abs(t.QLen[rows[('Q1', 'O')]] - 1.0413) < 1e-3
    assert abs(t.Util[rows[('Q1', 'O')]] - 0.2) < 1e-3
