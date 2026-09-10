"""Drift guard for the defaults and checks Network._sanitize applies at link().

Port of the second half of MATLAB `sanitize.m`, which native Python previously
carried only in part (the "Queue with no service" check). Every expectation here
was read off the MATLAB reference on the same model:

    SEPT schedparam row = 3 1 2
    LEPT schedparam row = 1 3 2
    SEPT duplicate means: SEPT does not support identical service time means.
    missing class service: Job class 'C2' has no service configured at any
                           station. ...

See _kb/11-conventions-and-gotchas.md.
"""
import numpy as np
import pytest

from line_solver import (ClosedClass, Delay, Disabled, Exp, Network, OpenClass,
                         Queue, SchedStrategy, Sink, Source)
from line_solver.lang.base import RoutingStrategy


def _cycle(model, classes, a, b):
    P = model.initRoutingMatrix()
    for c in classes:
        P.set(c, c, a, b, 1.0)
        P.set(c, c, b, a, 1.0)
    model.link(P)
    return model


def _three_class_queue(sched):
    m = Network('x')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', sched)
    cs = [ClosedClass(m, 'C%d' % i, 1, d) for i in (1, 2, 3)]
    for c in cs:
        d.setService(c, Exp(1.0))
    for c, rate in zip(cs, (1 / 3.0, 1.0, 1 / 2.0)):   # means 3, 1, 2
        q.setService(c, Exp(rate))
    return _cycle(m, cs, d, q), q


@pytest.mark.parametrize('sched,expected', [(SchedStrategy.SEPT, [3.0, 1.0, 2.0]),
                                            (SchedStrategy.LEPT, [1.0, 3.0, 2.0])])
def test_sept_lept_rank_the_class_means(sched, expected):
    # sn.schedparam carries the RANK of each class mean, ascending for SEPT and
    # descending for LEPT; the values are MATLAB's on the same model.
    model, _ = _three_class_queue(sched)
    sn = model.getStruct()
    assert np.asarray(sn.schedparam)[1, :3].tolist() == expected


@pytest.mark.parametrize('sched', [SchedStrategy.SEPT, SchedStrategy.LEPT])
def test_sept_lept_refuse_identical_means(sched):
    # Two classes with the same mean leave the order undefined; the reference
    # refuses rather than break the tie arbitrarily.
    m = Network('dup')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', sched)
    cs = [ClosedClass(m, 'C%d' % i, 1, d) for i in (1, 2)]
    for c in cs:
        d.setService(c, Exp(1.0))
        q.setService(c, Exp(0.5))
    with pytest.raises(ValueError, match='identical service time means'):
        _cycle(m, cs, d, q)


def test_a_closed_class_needs_service_somewhere():
    m = Network('ns')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', 1, d)
    c2 = ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    q.setService(c1, Exp(2.0))
    with pytest.raises(ValueError, match="Job class 'C2' has no service"):
        _cycle(m, [c1, c2], d, q)


def test_an_open_class_may_route_straight_to_the_sink():
    # An open class is exempt from the "service somewhere" check by design.
    m = Network('open')
    s = Source(m, 'S')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    k = Sink(m, 'K')
    c1 = OpenClass(m, 'C1')
    c2 = OpenClass(m, 'C2')
    s.setArrival(c1, Exp(1.0))
    s.setArrival(c2, Exp(0.5))
    q.setService(c1, Exp(2.0))
    P = m.initRoutingMatrix()
    P.set(c1, c1, s, q, 1.0)
    P.set(c1, c1, q, k, 1.0)
    P.set(c2, c2, s, k, 1.0)
    m.link(P)          # must not raise
    assert m.getStruct() is not None


def test_the_checks_are_skipped_when_checks_are_off():
    m = Network('ns')
    m.setChecks(False)
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', 1, d)
    c2 = ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    q.setService(c1, Exp(2.0))
    _cycle(m, [c1, c2], d, q)   # must not raise


def test_defaults_fill_the_unconfigured_slots():
    m = Network('defaults')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    c1 = ClosedClass(m, 'C1', 1, d)
    c2 = ClosedClass(m, 'C2', 1, d)
    d.setService(c1, Exp(1.0))
    d.setService(c2, Exp(1.0))
    q.setService(c1, Exp(2.0))          # C2 is NOT served at Q
    # C2 IS NOT ROUTED TO Q EITHER. Routing a class into a station that cannot
    # serve it is a flow sink -- it arrives and never leaves -- and link() now
    # refuses it: on the D <-> Q cycle this test used to build, MVA reported
    # Q/C2 ArvR 1 against Tput 0, CTMC dropped C2 altogether and SSA returned
    # D/C2 QLen 2e-06 with the Q rows missing. What this test is about is the
    # DEFAULTS sanitize fills in for a class a station does not serve, and
    # those are filled whether or not the class is routed there, so C2 keeps
    # its own self-loop at its reference station and every assertion below
    # still sees the Disabled service, the zero class cap and the DISABLED
    # routing strategy at Q.
    P = m.initRoutingMatrix()
    P.set(c1, c1, d, q, 1.0)
    P.set(c1, c1, q, d, 1.0)
    P.set(c2, c2, d, d, 1.0)
    m.link(P)

    from line_solver.distributions import Disabled
    assert isinstance(q._service_process[c2], Disabled)
    assert q._class_capacity[c2] == 0
    # a class a station cannot serve must not be routed to it
    strat = q._routing_strategies.get(c2)
    assert getattr(strat, 'name', None) == RoutingStrategy.DISABLED.name
    # and the Sink routes nothing at all
    sn = m.getStruct()
    assert sn is not None


def test_the_sink_routes_nothing():
    m = Network('sink')
    s = Source(m, 'S')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    k = Sink(m, 'K')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(Network.serial_routing([s, q, k]))
    strat = k._routing_strategies.get(c)
    assert getattr(strat, 'name', None) == RoutingStrategy.DISABLED.name


def test_a_class_routed_to_a_station_that_cannot_serve_it_is_refused():
    """An unservable class with incoming flow is a flow sink, not a modelling choice.

    Source -> Queue -> Sink with two open classes and a service for one of them.
    The Source releases both, so class B arrives at a Queue that cannot serve it
    and never leaves: SolverMVA used to report Q/B with ArvR 1 against Tput 0
    and sn.rates(Q,B) = NaN, with no warning of any kind.
    """
    m = Network('unservable')
    source = Source(m, 'Source')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    sink = Sink(m, 'Sink')
    a = OpenClass(m, 'A')
    b = OpenClass(m, 'B')
    source.setArrival(a, Exp(1.0))
    source.setArrival(b, Exp(1.0))
    q.setService(a, Exp(2.0))           # B is left unset
    m.addLink(source, q)
    m.addLink(q, sink)
    with pytest.raises(ValueError, match="routed to it"):
        m.getStruct()


def test_an_explicitly_disabled_class_routed_in_is_refused_too():
    """`Disabled()` declares "this class does not come here"; routing contradicts it.

    An absent slot and an explicit Disabled() reach the solvers as the same
    struct and produce the same wrong numbers, so the guard cannot let the
    spelling decide. What decides is whether anything routes the class in.
    """
    m = Network('unservable-explicit')
    source = Source(m, 'Source')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    sink = Sink(m, 'Sink')
    a = OpenClass(m, 'A')
    b = OpenClass(m, 'B')
    source.setArrival(a, Exp(1.0))
    source.setArrival(b, Exp(1.0))
    q.setService(a, Exp(2.0))
    q.setService(b, Disabled())
    m.addLink(source, q)
    m.addLink(q, sink)
    with pytest.raises(ValueError, match="routed to it"):
        m.getStruct()


def test_a_disabled_class_that_never_visits_is_accepted():
    """The class-switching idiom must keep working.

    Each queue serves exactly one class and marks the others Disabled; none of
    them is routed to a queue that cannot serve it, so nothing is a flow sink
    and the model must still build. This is the shape SolverCTMC's
    class-switching chain is written in.
    """
    m = Network('chain')
    source = Source(m, 'Source')
    q1 = Queue(m, 'Q1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Q2', SchedStrategy.FCFS)
    sink = Sink(m, 'Sink')
    a = OpenClass(m, 'A')
    b = OpenClass(m, 'B')
    source.setArrival(a, Exp(0.3))
    source.setArrival(b, Disabled())
    q1.setService(a, Exp(1.0))
    q1.setService(b, Disabled())
    q2.setService(a, Disabled())
    q2.setService(b, Exp(0.9))
    P = m.initRoutingMatrix()
    P.set(a, a, source, q1, 1.0)
    P.set(a, b, q1, q2, 1.0)
    P.set(b, b, q2, sink, 1.0)
    m.link(P)
    sn = m.getStruct()
    assert sn is not None


def test_the_guard_is_off_under_set_checks_false():
    """model.set_checks(False) must still be the escape hatch, as for every other check."""
    m = Network('unservable-unchecked')
    source = Source(m, 'Source')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    sink = Sink(m, 'Sink')
    a = OpenClass(m, 'A')
    b = OpenClass(m, 'B')
    source.setArrival(a, Exp(1.0))
    source.setArrival(b, Exp(1.0))
    q.setService(a, Exp(2.0))
    m.set_checks(False)
    m.addLink(source, q)
    m.addLink(q, sink)
    assert m.getStruct() is not None
