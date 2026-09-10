"""Regression tests for the derived START and PREEMPT event tags.

EventType.START marks a job BEGINNING or RESUMING its hold on a server;
EventType.PREEMPT marks a job holding a server being pushed back into the
buffer. Neither carries a clock: they are instantaneous tags on the arc of the
ARV or DEP that causes them, so they add no state, no synchronization and no
numerical change to anything the solvers already reported.

A SOURCE HAS NO SERVER TO SEIZE, so its row of both derived rates is zero: a job
is CREATED there rather than admitted to service, and the identity below is a
statement about stations that hold jobs. Every assertion here therefore reads the
queue's row, never the Source's.

Oracles, in increasing strength:
  - a non-preemptive station starts one service per departure, so
    getStartRate == getAvgTput exactly and getPreemptRate == 0;
  - at a preemptive station startRate == TN + preemptRate, because every job
    starts service once per entry into a server and every preemption is followed
    by exactly one later resume or restart;
  - preempt-resume and preempt-independent report the SAME preemption rate:
    which phase the displaced job resumes in is not a property of how often it
    is displaced;
  - the tags are NOT synchronizations, so the event filtration getGenerator
    returns keeps its length and the generator keeps its entries;
  - SolverSSA estimates the same identity over its sampled path.

The MATLAB twin is line-test.git/test/testsCTMC/test_ctmc_event_filtration.m, the JAR twin is
SolverCTMCEventFiltrationTest.java and the C++ twin is
cpp/tests/test_ctmc_event_filtration.cpp.
"""
import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, EventType, Exp, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source, SolverCTMC, SolverSSA)


def open_single_class(sched, nservers=1):
    """Source -> Queue(sched) -> Sink, one open class."""
    model = Network('mm1')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', sched)
    sink = Sink(model, 'Sink')
    cls = OpenClass(model, 'Class1')
    source.setArrival(cls, Exp(0.5))
    queue.setService(cls, Exp(1.0))
    if nservers != 1:
        queue.setNumberOfServers(nservers)
    model.link(Network.serialRouting(source, queue, sink))
    return model


def open_two_class_prio(sched):
    """Source -> Queue(sched) -> Sink, an urgent and a normal open class."""
    model = Network('mm1prio')
    source = Source(model, 'Source')
    queue = Queue(model, 'Queue', sched)
    sink = Sink(model, 'Sink')
    urgent = OpenClass(model, 'Urgent', 0)
    normal = OpenClass(model, 'Normal', 1)
    source.setArrival(urgent, Exp(0.4))
    source.setArrival(normal, Exp(0.4))
    queue.setService(urgent, Exp(1.0))
    queue.setService(normal, Exp(1.0))
    P = model.initRoutingMatrix()
    P[urgent] = Network.serialRouting(source, queue, sink)
    P[normal] = Network.serialRouting(source, queue, sink)
    model.link(P)
    return model


def ctmc_at(model, cutoff):
    return SolverCTMC(model, cutoff=cutoff)


@pytest.mark.parametrize('sched,nservers', [
    (SchedStrategy.FCFS, 1),
    (SchedStrategy.PS, 1),
    (SchedStrategy.FCFS, 2),
])
def test_nonpreemptive_starts_one_service_per_departure(sched, nservers):
    """Promotion without displacement: one start per completion, no preemption."""
    solver = ctmc_at(open_single_class(sched, nservers), 4)
    TN = solver.getAvgTput()
    startN = solver.getStartRate()
    preemptN = solver.getPreemptRate()
    np.testing.assert_allclose(startN[1, :], TN[1, :], atol=1e-12)
    np.testing.assert_allclose(preemptN[1, :], np.zeros(TN.shape[1]), atol=1e-12)
    # the Source creates jobs rather than admitting them to service
    np.testing.assert_allclose(startN[0, :], np.zeros(TN.shape[1]), atol=1e-12)


def test_preemptive_satisfies_the_start_identity():
    solver = ctmc_at(open_two_class_prio(SchedStrategy.FCFSPRPRIO), 4)
    TN = solver.getAvgTput()
    startN = solver.getStartRate()
    preemptN = solver.getPreemptRate()
    np.testing.assert_allclose(startN[1, :], TN[1, :] + preemptN[1, :], atol=1e-9)
    assert preemptN[1, 0] == pytest.approx(0.0, abs=1e-12)  # urgent is never displaced
    assert preemptN[1, 1] > 0.0                             # normal must be


def test_resume_and_independent_report_the_same_preemption_rate():
    """PR and PI differ in the resumed phase, not in how often a job is displaced."""
    pr = ctmc_at(open_two_class_prio(SchedStrategy.FCFSPRPRIO), 4).getPreemptRate()
    pi = ctmc_at(open_two_class_prio(SchedStrategy.FCFSPIPRIO), 4).getPreemptRate()
    np.testing.assert_allclose(pi, pr, atol=1e-12)


def test_closed_cyclic_satisfies_the_start_identity():
    model = Network('cyclic')
    think = Delay(model, 'Think')
    queue = Queue(model, 'Queue', SchedStrategy.FCFS)
    jobs = ClosedClass(model, 'Jobs', 3, think)
    think.setService(jobs, Exp(1.0))
    queue.setService(jobs, Exp(2.0))
    model.link(Network.serialRouting(think, queue))
    solver = SolverCTMC(model)
    TN = solver.getAvgTput()
    preemptN = solver.getPreemptRate()
    # both stations hold jobs here, so the identity holds on the full matrix
    np.testing.assert_allclose(solver.getStartRate(), TN + preemptN, atol=1e-9)
    np.testing.assert_allclose(preemptN, np.zeros_like(TN), atol=1e-12)


def test_the_tags_stay_out_of_the_event_filtration():
    """The derived filtration rides on arcs the generator already carries."""
    model = open_two_class_prio(SchedStrategy.FCFSPRPRIO)
    solver = ctmc_at(model, 3)
    infGen, eventFilt = solver.getGenerator()
    assert infGen.shape[0] == 38
    assert abs(np.abs(np.asarray(infGen.todense() if hasattr(infGen, 'todense')
                                 else infGen)).sum() - 106.0) < 1e-9
    filt = solver.getEventFiltration(EventType.PREEMPT)
    assert len(filt) == model.getNumberOfStations()
    assert len(filt[0]) == model.getNumberOfClasses()
    with pytest.raises(ValueError):
        solver.getEventFiltration(EventType.DEP)


@pytest.mark.parametrize('method', ['serial', 'nrm'])
def test_both_ssa_engines_agree_at_a_non_preemptive_station(method):
    """The counters must not depend on options.method.

    The serial engine accumulates the tag RATE of the enabled transitions, the
    NRM engine COUNTS the events it fires and divides by the simulated time. The
    two estimators agree in the limit and differ by simulation error at any
    finite budget, so at an FCFS station both must land on startRate == TN with
    no preemption.
    """
    solver = SolverSSA(open_single_class(SchedStrategy.FCFS),
                       samples=100000, seed=23000, method=method)
    solver.getAvgTable()
    TN = solver.getAvgTput()
    startN = solver.getStartRate()
    preemptN = solver.getPreemptRate()
    assert startN[1, 0] > 0.0
    np.testing.assert_allclose(startN[1, 0], TN[1, 0], atol=5e-3)
    assert preemptN[1, 0] == pytest.approx(0.0, abs=1e-12)
    assert startN[0, 0] == pytest.approx(0.0, abs=1e-12)  # the Source seizes nothing


def test_ssa_estimates_the_same_identity():
    """The sampled path carries the same identity, within simulation error.

    Native Python satisfies this at a PREEMPTIVE station. The MATLAB twin
    asserts it at a non-preemptive one instead, because MATLAB's SolverSSA
    mis-splits the per-class throughput of an FCFSPRPRIO station -- a
    pre-existing defect of that engine, not of the tags, recorded in
    _kb/06-solver-catalog.md.
    """
    model = open_two_class_prio(SchedStrategy.FCFSPRPRIO)
    solver = SolverSSA(model, samples=100000, seed=23000, method='serial', cutoff=6)
    solver.getAvgTable()
    TN = solver.getAvgTput()
    startN = solver.getStartRate()
    preemptN = solver.getPreemptRate()
    np.testing.assert_allclose(startN[1, :], TN[1, :] + preemptN[1, :], atol=1e-2)
    # the high class is never displaced and starts service once per arrival
    assert preemptN[1, 0] == pytest.approx(0.0, abs=1e-12)
    assert startN[1, 0] == pytest.approx(0.4, abs=2e-2)
    assert preemptN[1, 1] > 0.0
