"""SolverMAM on a discrete (slotted) time scale.

The model is recognized as discrete-time from its distributions alone and solved
by the Q-MAM discrete-time queues under the late arrival system with delayed
access. Targets are the Geo/Geo/1 closed form and the MATLAB numbers recorded in
_kb/06-solver-catalog.md, which LDES slotted independently corroborated.
"""

import numpy as np
import pytest

from line_solver import (DMAP, DiscreteUniform, Det, Exp, Geometric, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source, SolverMAM)


def _single(arrival, service, name='DT1'):
    model = Network(name)
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'C1')
    source.setArrival(job_class, arrival)
    queue.setService(job_class, service)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def test_geo_geo_1_matches_the_closed_form():
    a, s = 0.2, 0.5
    table = SolverMAM(_single(Geometric(a), Geometric(s))).getAvgTable()
    # LAS-DA: E[N] = A(1-A)/(S-A), E[T] = (1-A)/(S-A), U = A/S
    assert table.QLen[1] == pytest.approx(a * (1 - a) / (s - a), abs=1e-6)
    assert table.RespT[1] == pytest.approx((1 - a) / (s - a), abs=1e-6)
    assert table.Util[1] == pytest.approx(a / s, abs=1e-6)
    assert table.Tput[1] == pytest.approx(a, abs=1e-6)


def test_det_service_stays_on_the_lattice():
    # MATLAB reports 0.466667, LDES slotted 0.466615
    table = SolverMAM(_single(Geometric(0.2), Det(2))).getAvgTable()
    assert table.QLen[1] == pytest.approx(0.466667, abs=1e-5)
    assert table.Util[1] == pytest.approx(0.4, abs=1e-6)


def test_discrete_uniform_service():
    # MATLAB reports 0.465000, a direct LAS-DA slot recursion 0.464552
    table = SolverMAM(_single(Geometric(0.15), DiscreteUniform(1, 4))).getAvgTable()
    assert table.QLen[1] == pytest.approx(0.465000, abs=1e-5)
    assert table.Util[1] == pytest.approx(0.375, abs=1e-6)


def test_dmap_arrivals_reach_the_struct_and_are_solved():
    # Before the discrete-time path existed a DMAP reached sn.proc as None.
    # MATLAB reports 0.700000, LDES slotted 0.698332.
    D0 = np.array([[0.5, 0.2], [0.1, 0.6]])
    D1 = np.array([[0.25, 0.05], [0.1, 0.2]])
    table = SolverMAM(_single(DMAP(D0, D1), Geometric(0.6))).getAvgTable()
    assert table.QLen[1] == pytest.approx(0.700000, abs=1e-5)
    assert table.Util[1] == pytest.approx(0.5, abs=1e-6)


def test_tandem_reproduces_the_discrete_burke_result():
    # The stationary departure stream of a Geo/Geo/1 queue is Bernoulli, so the
    # second queue is EXACT: 0.8 = A(1-A)/(S-A) with A = 0.2, S = 0.4.
    model = Network('DTtandem')
    source = Source(model, 'Source')
    q1 = Queue(model, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Q2', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    job_class = OpenClass(model, 'C1')
    source.setArrival(job_class, Geometric(0.2))
    q1.setService(job_class, Geometric(0.5))
    q2.setService(job_class, Geometric(0.4))
    model.link(Network.serialRouting(source, q1, q2, sink))

    table = SolverMAM(model).getAvgTable()
    assert table.QLen[1] == pytest.approx(0.533333, abs=1e-4)
    assert table.QLen[2] == pytest.approx(0.800000, abs=1e-3)
    assert table.Util[1] == pytest.approx(0.4, abs=1e-6)
    assert table.Util[2] == pytest.approx(0.5, abs=1e-6)


def test_continuous_models_still_take_the_continuous_path():
    # A guard against the detection swallowing ordinary models
    table = SolverMAM(_single(Exp(0.2), Exp(0.5), 'MM1')).getAvgTable()
    assert table.QLen[1] == pytest.approx(0.4 / 0.6, abs=1e-4)


def test_mixing_a_dmap_with_a_continuous_law_is_refused():
    # A DMAP has no continuous-time reading: its (D0,D1) are probability
    # matrices, so the continuous machinery would return a wrong number silently.
    D0 = np.array([[0.5, 0.2], [0.1, 0.6]])
    D1 = np.array([[0.25, 0.05], [0.1, 0.2]])
    with pytest.raises(RuntimeError, match="slotted time scale"):
        SolverMAM(_single(DMAP(D0, D1), Exp(0.6))).getAvgTable()
