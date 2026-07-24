"""Regression test for the native Robust Queueing Network Analyzer (RQNA,
solver_rqna), MVA method='rqna'. RQNA characterizes each flow by its index of
dispersion for counts (IDC) and bounds the mean workload by robust queueing,
capturing bursty non-renewal (MMPP/MAP) arrivals.

Values are cross-validated to 10 digits against the MATLAB and JAR
implementations. Reference: W. Whitt and W. You (2018), "A Robust Queueing
Network Analyzer Based on Indices of Dispersion"."""

import numpy as np

from line_solver import (
    Network, Source, Queue, Sink, OpenClass, MMPP2, Exp, SchedStrategy, SolverMVA,
)


def _mmpp2_single_queue():
    m = Network('RqnaMMPP2')
    s = Source(m, 'Source')
    q = Queue(m, 'Queue', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'Class1')
    s.setArrival(c, MMPP2(0.5, 2.0, 0.5, 0.1))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def _mmpp2_tandem():
    m = Network('RqnaTandem')
    s = Source(m, 'Source')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'Class1')
    s.setArrival(c, MMPP2(0.5, 2.0, 0.5, 0.1))
    q1.setService(c, Exp(3.0))
    q2.setService(c, Exp(2.5))
    m.link(Network.serialRouting(s, q1, q2, k))
    return m


def _mm1():
    m = Network('RqnaMM1')
    s = Source(m, 'Source')
    q = Queue(m, 'Queue', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'Class1')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(s, q, k))
    return m


def test_rqna_mmpp2_single_queue():
    Q, U, R, T, _, _ = SolverMVA(_mmpp2_single_queue(), method='rqna').getAvg()
    assert abs(Q[1] - 8.9877675841) < 1e-6
    assert abs(U[1] - 0.8750000000) < 1e-9
    assert abs(R[1] - 5.1358671909) < 1e-6
    assert abs(T[1] - 1.7500000000) < 1e-9


def test_rqna_mmpp2_tandem():
    Q, U, R, T, _, _ = SolverMVA(_mmpp2_tandem(), method='rqna').getAvg()
    assert abs(Q[1] - 1.5277217866) < 1e-6
    assert abs(Q[2] - 2.6353748404) < 1e-6
    assert abs(R[1] - 0.8729838781) < 1e-6
    assert abs(R[2] - 1.5059284802) < 1e-6


def test_rqna_mm1_exact():
    # RQNA is exact on the M/M/1 special case (Poisson arrivals, IDC == 1).
    Q, U, R, T, _, _ = SolverMVA(_mm1(), method='rqna').getAvg()
    assert abs(Q[1] - 1.0) < 1e-6
    assert abs(R[1] - 1.0) < 1e-6


def test_rqna_default_autodispatch_multiqueue():
    # The 'default' method auto-selects RQNA for a multi-queue bursty open net.
    Q, U, R, T, _, _ = SolverMVA(_mmpp2_tandem()).getAvg()
    assert abs(Q[1] - 1.5277217866) < 1e-6
    assert abs(Q[2] - 2.6353748404) < 1e-6
