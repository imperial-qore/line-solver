"""sn_has_product_form must read the finite buffers.

Until 2026-08-16 no conjunct of sn_has_product_form looked at sn.cap,
sn.classcap or sn.droprule: its gates were scheduling, heterogeneous FCFS,
priorities, fork-join, state-dependent routing and the FCFS SCV. So
cqn_bas_blocking -- Queue1 blocking after service into a Queue2 that holds one
of the two circulating jobs -- reported a product form on a network whose
truncation couples the two station occupancies. The conjunct is now
`not sn_has_blocking(sn)`.

The predicate must NOT fire on a declared-but-unreachable buffer, nor on the
single-station M/M/1/K loss system whose truncated geometric distribution is a
product form over its one station.

Mirrored by line-test.git/test_sn_has_blocking.m, the SnApiTest case
bindingBufferExcludesProductForm and cpp/tests/test_sn_api.cpp.
"""

from line_solver import (Network, Queue, Delay, Source, Sink, ClosedClass,
                         OpenClass, Exp, SchedStrategy, DropStrategy, SolverMVA)
from line_solver.api.sn import sn_has_blocking, sn_has_product_form, sn_is_mm1k_loss


def build_bas():
    model = Network('bas')
    q1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(model, 'Class1', 2, q1, 0)
    q1.setService(c, Exp(1.0))
    q2.setService(c, Exp(0.8))
    q2.setCap(1)
    q1.setDropRule(c, DropStrategy.BAS)
    model.link(Network.serialRouting(q1, q2))
    return model


def build_closed(cap=None):
    model = Network('pf')
    d = Delay(model, 'Think')
    q = Queue(model, 'Q', SchedStrategy.PS)
    c = ClosedClass(model, 'C1', 3, d, 0)
    d.setService(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    if cap is not None:
        q.setCap(cap)
    model.link(Network.serialRouting(d, q))
    return model


def build_mm1k(k=3):
    model = Network('mm1k')
    s = Source(model, 'Source')
    q = Queue(model, 'Queue', SchedStrategy.FCFS)
    k_sink = Sink(model, 'Sink')
    q.setNumberOfServers(1)
    q.setCapacity(k)
    oc = OpenClass(model, 'Class1', 0)
    s.setArrival(oc, Exp(0.8))
    q.setService(oc, Exp(1.0))
    model.link(Network.serialRouting(s, q, k_sink))
    return model


def test_binding_buffer_excludes_product_form():
    sn = build_bas().getStruct()
    assert sn_has_blocking(sn)
    assert not sn_has_product_form(sn)


def test_bas_model_still_solves_through_sqd():
    # the product-form guard must not reach the BAS dispatch: the model routes
    # to Smith queue decomposition before it is consulted
    qn = SolverMVA(build_bas(), 'sqd').getAvgQLen()
    assert abs(float(qn.sum()) - 2.0) < 1e-6           # population conserved
    assert abs(float(qn.flatten()[1]) - 1.0822) < 1e-3


def test_non_binding_capacity_keeps_product_form():
    sn = build_closed(3).getStruct()
    assert not sn_has_blocking(sn)
    assert sn_has_product_form(sn)


def test_binding_capacity_without_blocking_rule_also_fires():
    sn = build_closed(1).getStruct()
    assert sn_has_blocking(sn)
    assert not sn_has_product_form(sn)


def test_single_station_loss_system_is_exempt():
    sn = build_mm1k(3).getStruct()
    assert sn_is_mm1k_loss(sn)
    assert not sn_has_blocking(sn)
    assert sn_has_product_form(sn)
