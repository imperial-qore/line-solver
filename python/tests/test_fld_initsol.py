"""The fluid initial condition: decoding sn.state, whatever layout it is in.

WHY THIS FILE EXISTS. `sn.state[isf]` IS NOT A PER-CLASS COUNT VECTOR. Its
layout depends on the station's scheduling: PS and INF carry per-class counts,
but FCFS, HOL, LCFS and SIRO carry the BUFFER ORDERING -- one entry per job,
holding that job's class -- so entry r there is the class of the r-th queued
job, not the number of class-r jobs.

Reading it as a count vector is a SILENT wrong answer, not an error, and it is
the worst kind: the fluid ODE conserves whatever population it is handed, so a
closed class whose REFERENCE STATION is FCFS started with one job instead of N
and every reported mean came back scaled by roughly 1/N with no warning. A
Queue(FCFS) <-> Delay cycle at N=6 returned [0.333 0.667] against the exact
[4.024 1.976]. The MATLAB reference has always decoded with State.toMarginal
(`solver_fluid_initsol.m`), as do the JAR and C++ ports; only native python
did not.

The assertion that catches this is CONSERVATION, not agreement with CTMC: the
fluid mean is an approximation and may legitimately differ, but the total
population it starts from is not up for approximation.
"""

import numpy as np
import pytest

from line_solver import (Network, Delay, Queue, ClosedClass, Erlang, Exp,
                         SchedStrategy, SolverFLD, SolverCTMC)


def _cycle(N, sched, service=None, nservers=1):
    """A closed two-station cycle whose REFERENCE STATION is the queue.

    The reference station is what makes this a regression test: sn.state places
    the whole population there, so a station whose state is buffer-ordered is
    exactly where the decode has to be right.
    """
    m = Network('cycle')
    q = Queue(m, 'Q1', sched)
    q.setNumberOfServers(nservers)
    d = Delay(m, 'Think')
    c = ClosedClass(m, 'C', N, q, 0)
    q.setService(c, service if service is not None else Exp(2.0))
    d.setService(c, Exp(1.0))
    m.link(Network.serialRouting(q, d))
    return m


# Every scheduling whose state is buffer-ordered, plus the count-layout ones as
# controls: the decode must be right for BOTH layouts, not swapped.
# HOL is buffer-ordered too, but SolverFLD's featset refuses it outright, so it
# has no initial condition to get wrong here.
@pytest.mark.parametrize('sched', [
    SchedStrategy.FCFS, SchedStrategy.LCFS, SchedStrategy.SIRO, SchedStrategy.PS,
])
@pytest.mark.parametrize('N', [1, 6, 20])
def test_a_queue_reference_station_starts_with_the_whole_population(sched, N):
    QN = np.asarray(SolverFLD(_cycle(N, sched), 'matrix').getAvgQLen()).ravel()
    assert QN.sum() == pytest.approx(N, rel=1e-6), (
        f'{sched} N={N}: the fluid ODE conserves what it is handed, so a total '
        f'of {QN.sum()} means the initial condition, not the dynamics, is wrong')


def test_a_multiserver_queue_reference_station_is_no_different():
    for nservers in (2, 4):
        QN = np.asarray(SolverFLD(_cycle(10, SchedStrategy.FCFS, nservers=nservers),
                                  'matrix').getAvgQLen()).ravel()
        assert QN.sum() == pytest.approx(10, rel=1e-6)


def test_multi_phase_service_at_a_buffered_reference_station():
    """The phase branch of the decode: jobs in the WAITING BUFFER carry no phase
    and are restarted in phase 1, only those in service carry kir. With one
    server and N=6 at an Erlang(3), five of the six are in the buffer."""
    m = _cycle(6, SchedStrategy.FCFS, service=Erlang.fitMeanAndOrder(0.5, 3))
    QN = np.asarray(SolverFLD(m, 'matrix').getAvgQLen()).ravel()
    assert QN.sum() == pytest.approx(6, rel=1e-6)


def test_a_multiclass_buffered_reference_station_keeps_each_class_apart():
    """Buffer ordering is where the layouts differ MOST: the state is a list of
    class ids, so a decode that reads position r as 'class r' mixes the classes
    as well as losing the count."""
    m = Network('mc')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    d = Delay(m, 'Think')
    a = ClosedClass(m, 'A', 3, q, 0)
    b = ClosedClass(m, 'B', 2, q, 0)
    q.setService(a, Exp(2.0)); q.setService(b, Exp(3.0))
    d.setService(a, Exp(1.0)); d.setService(b, Exp(1.0))
    P = m.initRoutingMatrix()
    P[a] = Network.serialRouting(q, d)
    P[b] = Network.serialRouting(q, d)
    m.link(P)
    QN = np.asarray(SolverFLD(m, 'matrix').getAvgQLen())
    assert QN.sum() == pytest.approx(5, rel=1e-6)


def test_the_answer_itself_and_not_only_its_total():
    """Conservation alone would pass on an initial condition that put the whole
    population at the WRONG station, so pin the mean against exact CTMC. The
    tolerance is the fluid approximation error, which is why it is not tight."""
    m = _cycle(6, SchedStrategy.FCFS)
    QN = np.asarray(SolverFLD(m, 'matrix').getAvgQLen()).ravel()
    exact = np.asarray(SolverCTMC(m, 'exact').getAvgQLen()).ravel()
    assert np.max(np.abs(QN - exact)) < 0.05, f'fluid {QN} vs exact {exact}'


def test_the_reference_station_is_what_decides_it():
    """The same two stations with the population referenced at the DELAY were
    correct all along, because a delay's state is already a count vector. Both
    orientations must now give the same physics."""
    m = Network('delayref')
    d = Delay(m, 'Think')
    q = Queue(m, 'Q1', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 6, d, 0)
    d.setService(c, Exp(1.0)); q.setService(c, Exp(2.0))
    m.link(Network.serialRouting(d, q))
    delay_ref = np.asarray(SolverFLD(m, 'matrix').getAvgQLen()).ravel()
    queue_ref = np.asarray(SolverFLD(_cycle(6, SchedStrategy.FCFS), 'matrix').getAvgQLen()).ravel()
    assert delay_ref.sum() == pytest.approx(6, rel=1e-6)
    # the cycle is the same model with the two stations swapped in the index order
    assert delay_ref[::-1] == pytest.approx(queue_ref, rel=1e-4)
