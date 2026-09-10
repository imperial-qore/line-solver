"""Utilization convention on G-networks solved by SolverCTMC.

Utilization at a station exposed to G-network signals is the CARRIED load, i.e.
the mean busy-server fraction. With exponential service the completion rate is
mu times the mean number of busy servers, so E[busy]/c = T*E[S]/c holds exactly
even though signals remove jobs without a completion; an arrival-based estimator
lambda*E[S]/c would instead measure the OFFERED load (that is the defect fixed
in the MATLAB and JAR analyzers by ctmc_signal_lossy / CtmcSignalLossy; native
Python computes only T*E[S]/c and was already correct).

Cases:
  1. single-server FCFS, untargeted negative customers  (exact Gelenbe)
  2. single-server PS, untargeted negative customers    (exact, same rho)
  3. single-server FCFS, catastrophes                   (exact, quadratic root)
  4. PS, two positive classes, signal targeting one     (victim class only)
  5. multiserver c=2, negative customers                (MATLAB CTMC golden)
  6. Geometric batch removal                            (MATLAB CTMC golden)
  7. Erlang-2 service, negative customers               (busy fraction > T*E[S])

Cases 5 to 7 have no closed form; their goldens are the MATLAB CTMC values,
themselves validated against the LDES sample path in
line-test.git/test/testsAdvFeatures/des/test_gnetwork_ctmc_util.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import math

import pytest

from line_solver import (Network, Source, Queue, Sink, OpenClass, Signal,
                         SignalType, SchedStrategy, Exp, Geometric,
                         RemovalPolicy, Erlang, SolverCTMC)

MU = 1.0
TOL = 1e-6

# EVERY CUTOFF HERE IS THE SMALLEST ONE THAT IS CONVERGED, not a round number.
# The queue length decays geometrically at the case's own rate, so the truncated
# tail is what bounds the error, and each cutoff below is set where that tail
# sits at least two orders of magnitude inside the assertion it has to satisfy:
# 20 leaves 8.5e-9 on the rho = 0.357 cases against TOL, 30 leaves 5.0e-7 on the
# batch case against its 1e-4. A LARGER CUTOFF IS NOT FREE ON A SIGNAL MODEL --
# the space is linear in it only because a signal class is capped out of the
# station buffer, and the C++ engine, which did not apply that cap, enumerated
# the buffer orderings instead and was OOM-killed at 43.6 GB on cutoff 30.


def _gnetwork1(sched, nservers, lambda_pos, mu, lambda_neg, signal_type,
               rem_dist=None, rem_policy=None):
    """Single queue with one positive class and one signal class."""
    model = Network('GNetworkUtil')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', sched)
    queue.setNumberOfServers(nservers)
    sink = Sink(model, 'Sink')

    pos = OpenClass(model, 'Positive')
    source.setArrival(pos, Exp(lambda_pos))
    queue.setService(pos, Exp(mu))

    if rem_dist is None:
        neg = Signal(model, 'Negative', signal_type)
    else:
        neg = Signal(model, 'Negative', signal_type, 0, rem_dist, rem_policy)
    source.setArrival(neg, Exp(lambda_neg))
    queue.setService(neg, Exp(mu))

    P = model.initRoutingMatrix()
    P.set(pos, pos, source, queue, 1.0)
    P.set(pos, pos, queue, sink, 1.0)
    P.set(neg, neg, source, queue, 1.0)
    P.set(neg, neg, queue, sink, 1.0)
    model.link(P)
    return model


def test_fcfs_exact():
    """rho = lambda+/(mu + lambda-) by flow balance: every arrival either
    completes or is removed while the server is busy."""
    lambda_pos, lambda_neg = 0.5, 0.4
    model = _gnetwork1(SchedStrategy.FCFS, 1, lambda_pos, MU, lambda_neg,
                       SignalType.NEGATIVE)
    solver = SolverCTMC(model, cutoff=20)
    rho = lambda_pos / (MU + lambda_neg)
    assert solver.getAvgUtil()[1][0] == pytest.approx(rho, abs=TOL)
    assert solver.getAvgQLen()[1][0] == pytest.approx(rho / (1 - rho), abs=TOL)
    assert solver.getAvgTput()[1][0] == pytest.approx(MU * rho, abs=TOL)
    # The offered load lambda+/mu = 0.5 is what an arrival-based estimator gives.
    assert abs(solver.getAvgUtil()[1][0] - lambda_pos / MU) > 0.1


def test_ps_exact():
    """The Gelenbe product form is discipline-invariant here, so PS must return
    the same rho and mean queue length as FCFS."""
    lambda_pos, lambda_neg = 0.5, 0.4
    model = _gnetwork1(SchedStrategy.PS, 1, lambda_pos, MU, lambda_neg,
                       SignalType.NEGATIVE)
    solver = SolverCTMC(model, cutoff=20)
    rho = lambda_pos / (MU + lambda_neg)
    assert solver.getAvgUtil()[1][0] == pytest.approx(rho, abs=TOL)
    assert solver.getAvgQLen()[1][0] == pytest.approx(rho / (1 - rho), abs=TOL)


def test_catastrophe_exact():
    """Catastrophes flush the station: p_n = (1-r) r^n with
    mu r^2 - (lambda + mu + delta) r + lambda = 0, so Util = P(busy) = r."""
    lambda_pos, delta = 0.5, 0.4
    model = _gnetwork1(SchedStrategy.FCFS, 1, lambda_pos, MU, delta,
                       SignalType.CATASTROPHE)
    solver = SolverCTMC(model, cutoff=20)
    b = lambda_pos + MU + delta
    r = (b - math.sqrt(b * b - 4 * lambda_pos * MU)) / (2 * MU)
    assert solver.getAvgUtil()[1][0] == pytest.approx(r, abs=TOL)
    assert solver.getAvgQLen()[1][0] == pytest.approx(r / (1 - r), abs=TOL)
    assert solver.getAvgTput()[1][0] == pytest.approx(MU * r, abs=TOL)


def test_targeted_class():
    """A targeted signal downgrades ONLY its victim class: C1 is never removed,
    so all its arrivals complete and its utilization stays at lambda1/mu."""
    lambda1, lambda2, lambda_neg = 0.2, 0.2, 0.5
    model = Network('GNetworkTargeted')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.PS)
    sink = Sink(model, 'Sink')
    c1 = OpenClass(model, 'C1')
    c2 = OpenClass(model, 'C2')
    source.setArrival(c1, Exp(lambda1))
    queue.setService(c1, Exp(MU))
    source.setArrival(c2, Exp(lambda2))
    queue.setService(c2, Exp(MU))
    neg = Signal(model, 'Negative', SignalType.NEGATIVE).forJobClass(c2)
    source.setArrival(neg, Exp(lambda_neg))
    queue.setService(neg, Exp(MU))
    P = model.initRoutingMatrix()
    P.set(c1, c1, source, queue, 1.0)
    P.set(c1, c1, queue, sink, 1.0)
    P.set(c2, c2, source, queue, 1.0)
    P.set(c2, c2, queue, sink, 1.0)
    P.set(neg, neg, source, queue, 1.0)
    P.set(neg, neg, queue, sink, 1.0)
    model.link(P)

    solver = SolverCTMC(model, cutoff=12)
    UN = solver.getAvgUtil()
    TN = solver.getAvgTput()
    assert UN[1][0] == pytest.approx(lambda1 / MU, abs=1e-4)
    assert UN[1][1] < lambda2 / MU - 0.01
    assert UN[1][1] == pytest.approx(TN[1][1] / MU, abs=1e-4)
    assert UN[1][1] == pytest.approx(0.12482827, abs=1e-4)


def test_multiserver():
    """No closed form; golden is the MATLAB CTMC value (LDES-validated)."""
    lambda_pos, lambda_neg, c = 1.2, 0.3, 2
    model = _gnetwork1(SchedStrategy.FCFS, c, lambda_pos, MU, lambda_neg,
                       SignalType.NEGATIVE)
    solver = SolverCTMC(model, cutoff=25)
    UN = solver.getAvgUtil()
    TN = solver.getAvgTput()
    assert UN[1][0] == pytest.approx(TN[1][0] / (c * MU), abs=1e-8)
    assert UN[1][0] < lambda_pos / (c * MU) - 0.05
    assert UN[1][0] == pytest.approx(0.50119329, abs=1e-4)


def test_phase_type_service_busy_fraction():
    """With Erlang-2 service the destroyed jobs' partial service is real busy
    time with no completion, so T*E[S]/c (0.34941) under-counts the busy
    fraction; the state-space occupancy (0.37648) is the exact value, and it is
    what the LDES sample path measures (0.37696 at 5e5 samples)."""
    lambda_pos, lambda_neg = 0.5, 0.4
    model = Network('GNetworkPH')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    pos = OpenClass(model, 'Positive')
    source.setArrival(pos, Exp(lambda_pos))
    queue.setService(pos, Erlang.fitMeanAndOrder(1.0, 2))
    neg = Signal(model, 'Negative', SignalType.NEGATIVE)
    source.setArrival(neg, Exp(lambda_neg))
    queue.setService(neg, Exp(MU))
    P = model.initRoutingMatrix()
    P.set(pos, pos, source, queue, 1.0)
    P.set(pos, pos, queue, sink, 1.0)
    P.set(neg, neg, source, queue, 1.0)
    P.set(neg, neg, queue, sink, 1.0)
    model.link(P)

    solver = SolverCTMC(model, cutoff=20)
    UN = solver.getAvgUtil()
    TN = solver.getAvgTput()
    QN = solver.getAvgQLen()
    assert TN[1][0] == pytest.approx(0.34940753, abs=1e-6)
    assert QN[1][0] == pytest.approx(0.55459104, abs=1e-6)
    assert UN[1][0] == pytest.approx(0.37648116, abs=1e-6)
    # Strictly above the carried-load value: partial service is not a completion.
    assert UN[1][0] > TN[1][0] + 0.02
    # Flow balance: arrivals = completions + removals, removals = lambda- * Util.
    assert TN[1][0] + lambda_neg * UN[1][0] == pytest.approx(lambda_pos, abs=1e-6)


def test_batch_removal():
    """Geometric batch size clipped at the population; golden is the MATLAB
    CTMC value (LDES-validated)."""
    lambda_pos, lambda_neg = 0.8, 0.3
    model = _gnetwork1(SchedStrategy.FCFS, 1, lambda_pos, MU, lambda_neg,
                       SignalType.NEGATIVE, Geometric(0.5), RemovalPolicy.RANDOM)
    solver = SolverCTMC(model, cutoff=30)
    UN = solver.getAvgUtil()
    TN = solver.getAvgTput()
    assert UN[1][0] == pytest.approx(TN[1][0] / MU, abs=1e-8)
    assert UN[1][0] < lambda_pos / MU - 0.05
    assert UN[1][0] == pytest.approx(0.56421833, abs=1e-4)
