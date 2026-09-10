"""
Round-robin dispatching in the two decomposition methods, MVA 'qna' and MAM
'mna'.

A round-robin node splits the departure stream one-in-k, so each destination
sees the k-fold convolution of the interarrival time: with a Poisson source and
k=2 the arrivals at each queue are Erlang-2 (SCV 1/2) and the queue is shorter
than under Bernoulli routing at the same rate. Reference values from the MATLAB
twin (solver_qna.m, solver_mna_open.m); RAND must return the exact M/M/1 value,
which pins the k=1 branch to the original Markovian formula.
"""

import numpy as np
import pytest

from line_solver import (Network, Source, Router, Queue, Sink, OpenClass, Exp,
                         SchedStrategy, RoutingStrategy, SolverMVA, SolverMAM)
from line_solver.api.npfqn import npfqn_traffic_split_rr
from line_solver.api.sn import sn_rt_stations

TOL = 1e-3


def rr_model(strategy, lam=1.4, mu=1.0):
    model = Network('rr_split')
    source = Source(model, 'Source')
    router = Router(model, 'Router')
    q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    oclass = OpenClass(model, 'Class1')
    source.setArrival(oclass, Exp(lam))
    q1.setService(oclass, Exp(mu))
    q2.setService(oclass, Exp(mu))
    model.addLink(source, router)
    model.addLink(router, q1)
    model.addLink(router, q2)
    model.addLink(q1, sink)
    model.addLink(q2, sink)
    router.setRouting(oclass, strategy)
    return model


def test_split_degree_and_station_projection():
    sn = rr_model(RoutingStrategy.RROBIN).getStruct()
    # the source feeds a round-robin router with two outlinks
    assert npfqn_traffic_split_rr(sn).ravel().tolist() == [2.0, 1.0, 1.0]
    # RAND leaves every split Markovian
    sn_rand = rr_model(RoutingStrategy.RAND).getStruct()
    assert npfqn_traffic_split_rr(sn_rand).ravel().tolist() == [1.0, 1.0, 1.0]
    # sn.rt is stateful-indexed (4x4 here): the station projection absorbs the router
    rt, V = sn_rt_stations(sn)
    assert rt.shape == (3, 3)
    assert rt[0, 1] == pytest.approx(0.5)
    assert rt[0, 2] == pytest.approx(0.5)
    assert V.ravel().tolist() == [1.0, 0.5, 0.5]


def test_qna_separates_round_robin_from_bernoulli():
    qrr = np.asarray(SolverMVA(rr_model(RoutingStrategy.RROBIN), 'qna').getAvgQLen()).ravel()
    assert qrr[1] == pytest.approx(1.9250, abs=TOL)
    assert qrr[2] == pytest.approx(1.9250, abs=TOL)

    qrand = np.asarray(SolverMVA(rr_model(RoutingStrategy.RAND), 'qna').getAvgQLen()).ravel()
    # Bernoulli split of a Poisson stream is Poisson: exact M/M/1
    assert qrand[1] == pytest.approx(2.3333, abs=TOL)
    assert qrand[2] == pytest.approx(2.3333, abs=TOL)
    assert qrr[1] < qrand[1]


def test_mna_separates_round_robin_from_bernoulli():
    qrr = np.asarray(SolverMAM(rr_model(RoutingStrategy.RROBIN), 'mna').getAvgQLen()).ravel()
    # E_2/M/1 solved by MMAPPH1FCFS on the two-moment fit of the split flow
    assert qrr[1] == pytest.approx(1.8204, abs=TOL)
    assert qrr[2] == pytest.approx(1.8204, abs=TOL)

    qrand = np.asarray(SolverMAM(rr_model(RoutingStrategy.RAND), 'mna').getAvgQLen()).ravel()
    assert qrand[1] == pytest.approx(2.3333, abs=TOL)
    assert qrand[2] == pytest.approx(2.3333, abs=TOL)
    assert qrr[1] < qrand[1]


def test_round_robin_rejected_outside_qna_and_mna():
    with pytest.raises(RuntimeError):
        SolverMVA(rr_model(RoutingStrategy.RROBIN), 'amva').getAvgQLen()
