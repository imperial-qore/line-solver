"""Derivation of the QRF BAS blocking tables from an sn.

The tables used to be hand-built by the caller and the solver refused a model
without them, which made ``qrf.bas`` unreachable from a plain
``SolverBA(model, 'qrf.bas')``. Everything the enumeration needs is implied by
the model, so the expected values below are not a record of what this
implementation happens to produce: they are what the blocking structure of each
model forces.
"""

import numpy as np
import pytest

from line_solver import (ClosedClass, DropStrategy, Exp, Network, Queue,
                         SchedStrategy, SolverBA)
from line_solver.api.sn.qrf_blocking import sn_to_qrf_blocking, sn_to_qrf_capacity


class _Opt:
    def __init__(self, **config):
        self.config = config


def _cqn_bas_blocking():
    """Queue1 -BAS-> Queue2(cap 1), N = 2. The in-tree example."""
    m = Network('cqn_bas_blocking')
    q1 = Queue(m, 'Queue1', SchedStrategy.FCFS)
    q2 = Queue(m, 'Queue2', SchedStrategy.FCFS)
    c = ClosedClass(m, 'Class1', 2, q1, 0)
    q1.set_service(c, Exp(1.0))
    q2.set_service(c, Exp(0.8))
    q2.setCap(1)
    q1.setDropRule(c, DropStrategy.BAS)
    m.link(Network.serialRouting(q1, q2))
    return m


def _two_feeders(N=4, cap3=1):
    """Q1 -> {Q2, Q3(cap 1)}, Q2 -> Q3, Q3 -> Q1: two queues can block behind Q3."""
    m = Network('twoFeeders')
    qs = [Queue(m, 'Q%d' % (i + 1), SchedStrategy.FCFS) for i in range(3)]
    c = ClosedClass(m, 'C', N, qs[0], 0)
    for q, r in zip(qs, (1.0, 0.9, 1.2)):
        q.set_service(c, Exp(r))
    qs[2].setCap(cap3)
    for q in qs[:2]:
        q.setDropRule(c, DropStrategy.BAS)
    P = m.initRoutingMatrix()
    P.set(c, c, qs[0], qs[1], 0.5)
    P.set(c, c, qs[0], qs[2], 0.5)
    P.set(c, c, qs[1], qs[2], 1.0)
    P.set(c, c, qs[2], qs[0], 1.0)
    m.link(P)
    return m


def test_cqn_bas_blocking_tables():
    blk, msg = sn_to_qrf_blocking(_cqn_bas_blocking().getStruct())
    assert msg == ''
    assert blk['f'] == 2          # Queue2 is the one binding buffer
    assert blk['MR'] == 2         # empty configuration, plus Queue1 blocked
    assert blk['ZM'] == 1
    assert blk['blockers'] == [1]
    # F is an occupancy bound: unbounded stations are capped by the population
    assert list(blk['F']) == [2, 1]
    # Invariant 1: configuration 1 is the empty one, which ZERO4/5/7/8 assume
    assert list(blk['BB'][0]) == [0, 0] and blk['ZZ'][0] == 0
    assert blk['MM'][0][0] == 0
    # Configuration 2: Queue1 blocked, and it heads the FIFO blocking order
    assert list(blk['BB'][1]) == [1, 0] and blk['ZZ'][1] == 1
    assert blk['MM'][1][0] == 1
    # MM1 is indexed by the queue that BECOMES blocked, not by f
    assert list(blk['MM1'][0]) == [2, 0]
    assert list(blk['MM1'][1]) == [0, 0]


def test_capacity_reports_only_binding_buffers():
    F, binding, msg = sn_to_qrf_capacity(_cqn_bas_blocking().getStruct())
    assert msg == ''
    assert list(F) == [2, 1]
    # Queue1's capacity is the chain population refreshCapacity derives, which
    # can never refuse a job, so it must not read as a buffer.
    assert list(binding) == [False, True]


def test_every_blocking_order_is_enumerated():
    blk, msg = sn_to_qrf_blocking(_two_feeders().getStruct())
    assert msg == ''
    assert blk['f'] == 3
    assert blk['ZM'] == 2
    assert blk['MR'] == 5  # empty + 2 singletons + 2 orders of {Q1, Q2}
    assert blk['blockers'] == [1, 2]
    for m in range(blk['MR']):
        assert blk['ZZ'][m] == int(np.sum(blk['BB'][m]))
        assert blk['BB'][m][blk['f'] - 1] == 0  # f is never blocked behind itself


def test_successor_map_preserves_the_head():
    """Invariant 3: blocking appends at the tail, so MM1 keeps MM(m,0)."""
    blk, msg = sn_to_qrf_blocking(_two_feeders().getStruct())
    assert msg == ''
    for m in range(blk['MR']):
        if blk['ZZ'][m] >= blk['ZM']:
            continue
        for j in blk['blockers']:
            if blk['BB'][m][j - 1]:
                continue
            mp = blk['MM1'][m][j - 1]
            assert mp >= 1, 'every feeder not yet blocked must have a successor'
            assert blk['ZZ'][mp - 1] == blk['ZZ'][m] + 1
            assert blk['BB'][mp - 1][j - 1] == 1
            if blk['ZZ'][m] > 0:
                assert blk['MM'][mp - 1][0] == blk['MM'][m][0]


def test_refuses_more_than_one_binding_buffer():
    sn = _two_feeders().getStruct()
    sn.cap = np.array([4.0, 1.0, 1.0])
    sn.classcap = np.array([[4.0], [1.0], [1.0]])
    blk, msg = sn_to_qrf_blocking(sn)
    assert blk is None
    # qrf.bas carries a scalar f; qrf.rsrd sums PBB over every full queue
    assert 'single finite-capacity queue' in msg
    assert 'qrf.rsrd' in msg


def test_oversized_enumeration_is_refused_not_truncated():
    blk, msg = sn_to_qrf_blocking(_two_feeders().getStruct(), _Opt(qrf_maxvars=10))
    assert blk is None
    assert 'cannot be truncated' in msg


def test_no_reachable_blocking_state():
    # f can hold all 3 jobs, so nothing can ever be held behind it.
    blk, msg = sn_to_qrf_blocking(_two_feeders(N=3, cap3=3).getStruct())
    assert msg == ''
    assert blk['ZM'] == 0 and blk['MR'] == 1 and blk['ZZ'][0] == 0


@pytest.mark.parametrize('method', ['qrf.bas', 'qrf.rsrd'])
def test_solver_runs_without_hand_built_tables(method):
    """The regression this whole change exists for: a plain call must work."""
    table = SolverBA(_cqn_bas_blocking(), method).getAvgTable()
    assert table is not None and len(table) == 2


def test_qrf_bas_matches_the_exact_chain_on_cqn_bas_blocking():
    """The derived LP pins both utilizations, and they are the exact ones.

    States (n1, n2, blocked): (2,0) -1-> (1,1) -1-> (1,1)+Q1 blocked, with Q2
    completing at 0.8 in both directions, so pi = (0.262295, 0.327869, 0.409836).
    Queue1 is SERVING in the first two states and BLOCKED in the third, holding a
    job it has already finished, so its utilization is pi0 + pi1 = 0.590164.
    Queue2 is busy in the last two, pi1 + pi2 = 0.737705. Both are what
    SolverCTMC reports.

    These used to be [1, 0.737705]: the objective summed p2 over every blocking
    configuration, which is P(n >= 1) -- occupancy, counting a blocked server as
    busy. The objective is the e variables now, which UEFF restricts to the
    configurations where the station is NOT blocked, so Queue1 reads 0.590164
    instead of a vacuous 1.
    """
    util = np.asarray(SolverBA(_cqn_bas_blocking(), 'qrf.bas').getAvgUtil()).ravel()
    assert util[0] == pytest.approx(0.590164, abs=1e-5)
    assert util[1] == pytest.approx(0.737705, abs=1e-5)


def test_qrf_bas_throughput_stops_inheriting_the_occupancy_inflation():
    """U = X*V*s holds for a single server only when U is the BUSY fraction.

    With occupancy in U the derived throughput came out [1, 1] against an exact
    0.590164 -- a valid upper bound, but 69% high. With utilization it is exact.
    """
    tput = np.asarray(SolverBA(_cqn_bas_blocking(), 'qrf.bas').getAvgTput()).ravel()
    assert tput[0] == pytest.approx(0.590164, abs=1e-5)
    assert tput[1] == pytest.approx(0.590164, abs=1e-5)


# --------------------------------------------------------------------------
# 'default' on a blocked model. This is the call in the original report:
# SolverBA(model) used to refuse outright, because it means the geometric upper
# bound and that bound is blind to the buffer. It now means 'qrf.bas', which
# models the buffer and derives its own tables.
# --------------------------------------------------------------------------

def _delay_blocked():
    """A blocked model outside the QRF shape: the delay station rules it out."""
    from line_solver import Delay
    m = Network('delayBlocked')
    d = Delay(m, 'D')
    q = Queue(m, 'Q', SchedStrategy.FCFS)
    c = ClosedClass(m, 'C', 3, d, 0)
    d.set_service(c, Exp(1.0))
    q.set_service(c, Exp(0.8))
    q.setCap(1)
    m.link(Network.serialRouting(d, q))
    return m


def test_default_routes_to_qrf_bas_on_a_blocked_model():
    from line_solver.solvers.solver_ba.solver_ba_analyzer import ba_blocking_default
    sn = _cqn_bas_blocking().getStruct()
    method, why = ba_blocking_default(sn)
    assert method == 'qrf.bas'
    assert why == ''

    routed = np.asarray(SolverBA(_cqn_bas_blocking()).getAvgUtil()).ravel()
    named = np.asarray(SolverBA(_cqn_bas_blocking(), 'qrf.bas').getAvgUtil()).ravel()
    assert np.allclose(routed, named)


def test_default_is_offered_back_by_the_method_list():
    valid = SolverBA(_cqn_bas_blocking()).list_valid_methods()
    assert 'default' in valid
    # and nothing blocking-blind survives beside it
    from line_solver.solvers.solver_ba.solver_ba_analyzer import (ba_ignores_blocking,
                                                                  ba_resolve_method)
    for m in valid:
        if m == 'default':
            continue
        assert not ba_ignores_blocking(ba_resolve_method(m))


def test_unroutable_blocked_model_still_refuses_and_says_why():
    from line_solver.solvers.solver_ba.solver_ba_analyzer import ba_blocking_default
    method, why = ba_blocking_default(_delay_blocked().getStruct())
    assert method == ''
    assert 'delay station' in why
    with pytest.raises(ValueError) as e:
        SolverBA(_delay_blocked()).getAvgTable()
    assert 'does not support finite-buffer blocking' in str(e.value)
    assert 'do not apply here either' in str(e.value)


def test_get_bounds_refuses_the_routed_default_as_one_sided():
    # qrf.bas is solved in the 'max' direction alone, so there is no bracket --
    # and saying "does not support blocking" here would contradict the run.
    with pytest.raises(ValueError) as e:
        SolverBA(_cqn_bas_blocking()).get_bounds()
    assert 'UPPER-only' in str(e.value)
