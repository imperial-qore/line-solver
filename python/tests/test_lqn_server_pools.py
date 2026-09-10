"""
Heterogeneous server pools with a class-compatibility graph on a LayeredNetwork
server, and their lowering onto the layer station by SolverLN.

The pools are declared on a Processor (operands = its tasks) or on a Task
(operands = its entries) with ``addServerType``, and SolverLN lowers them to the
activated-server rate of ``sn_compat_rate``, carried as a joint dependence. See
_kb/06-solver-catalog.md (LN section).
"""

import numpy as np
import pytest

from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                         Exp, SchedStrategy, LN, MVA)
from line_solver.api.sn import sn_compat_rate, sn_compat_peak, sn_compat_scaling
from line_solver.lang.constant.server_type import ServerType

# s1 -> op0 only, s2 -> {op0, op1}, s3 -> op1 only
COMPAT = [[1, 0], [1, 1], [0, 1]]
COUNTS = [1, 1, 1]
RATES = [1, 1, 1]


def test_compat_rate_integer_states():
    """The activated-server rate at the integer states of the 3x2 pool."""
    assert sn_compat_rate(COMPAT, COUNTS, RATES, [2, 0]) == pytest.approx(2.0)
    assert sn_compat_rate(COMPAT, COUNTS, RATES, [2, 1]) == pytest.approx(3.0)
    assert sn_compat_rate(COMPAT, COUNTS, RATES, [0, 1]) == pytest.approx(2.0)
    assert sn_compat_rate(COMPAT, COUNTS, RATES, [0, 0]) == pytest.approx(0.0)
    assert sn_compat_peak(COUNTS, RATES) == pytest.approx(3.0)


def test_compat_rate_matches_hard_indicator_on_the_lattice():
    """
    min(1, .) and the hard indicator agree at EVERY integer state.

    This is the invariant that keeps the order-independent law intact: the
    fractional relaxation exists only so a mean-value solver can see the
    structure, and CTMC and simulation, which evaluate at integer states alone,
    must not be able to tell the two apart.
    """
    compat = np.asarray(COMPAT, dtype=float)
    for n0 in range(4):
        for n1 in range(4):
            n = np.array([n0, n1], dtype=float)
            hard = float(np.sum(np.asarray(COUNTS)[np.any((compat != 0) & (n > 0), axis=1)]))
            assert sn_compat_rate(COMPAT, COUNTS, RATES, n) == pytest.approx(hard)


def test_compat_rate_scales_below_one_job():
    """Below one compatible job the pool contributes proportionally, not fully."""
    # only op0 present, at half a job: pools s1 and s2 reach it, each at weight 0.5
    assert sn_compat_rate(COMPAT, COUNTS, RATES, [0.5, 0.0]) == pytest.approx(1.0)
    # a pool sees the COMBINED load of the operands it can reach
    assert sn_compat_rate([[1, 1]], [1], [1], [0.4, 0.7]) == pytest.approx(1.0)
    assert sn_compat_rate([[1, 1]], [1], [1], [0.4, 0.2]) == pytest.approx(0.6)


def _full_pool_eta(n, nservers=3):
    """
    The scaling a fully-compatible pool of `nservers` imposes, written out.

    Deliberately NOT a call to sn_compat_scaling: a test that compared the
    implementation with itself would have stayed green through the 2026-08-28
    change of the denominator, which is the one thing here worth pinning.
    """
    total = float(np.sum(np.maximum(np.asarray(n, dtype=float), 0.0)))
    if total <= 0:
        return 1.0
    return min(1.0, total) / min(1.0, total / nservers)


def _build(kind, mult=3):
    model = LayeredNetwork('LQNpools')
    P0 = Processor(model, 'P0', 1, SchedStrategy.PS)
    P1 = Processor(model, 'P1', mult, SchedStrategy.PS)
    T1 = Task(model, 'T1', 6, SchedStrategy.REF).on(P0).setThinkTime(Exp(1.0))
    T2 = Task(model, 'T2', 6, SchedStrategy.FCFS).on(P1)
    T3 = Task(model, 'T3', 6, SchedStrategy.FCFS).on(P1)
    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T3)
    Activity(model, 'A1', Exp(2.0)).on(T1).boundTo(E1).synchCall(E2, 1).synchCall(E3, 1)
    Activity(model, 'A2', Exp(3.0)).on(T2).boundTo(E2).repliesTo(E2)
    Activity(model, 'A3', Exp(2.0)).on(T3).boundTo(E3).repliesTo(E3)
    if kind == 'compat':
        P1.addServerType(ServerType('S1', 1, [T2]))
        P1.addServerType(ServerType('S2', 1, [T2, T3]))
        P1.addServerType(ServerType('S3', 1, [T3]))
    elif kind == 'full':
        P1.addServerType(ServerType('All', 3, [T2, T3]))
    elif kind == 'eta1':
        P1.setJointDependence(lambda n: 1.0, np.array([1.0, 1.0]))
    elif kind == 'etafull':
        P1.setJointDependence(_full_pool_eta, np.array([1.0, 1.0]))
    elif kind == 'unserved':
        P1.addServerType(ServerType('OnlyT2', 3, [T2]))
    elif kind == 'badmult':
        P1.addServerType(ServerType('TooFew', 2, [T2, T3]))
    return model


def _qlen(kind):
    t = LN(_build(kind), lambda mm: MVA(mm), method='srvn.cs').getAvgTable()
    r = t[t['Node'].isin(['T2', 'T3'])]
    return np.array([r['QLen'].iloc[0], r['QLen'].iloc[1]])


def test_full_pool_scaling_is_the_redundancy_speed_up():
    """
    The law itself, at the states that distinguish it from the neutral one.

    eta = min(1, N) / min(1, N/S): S below one job, S/N from one job to S, and 1
    from S jobs up. It is ABOVE ONE at low occupancy -- the speed-up of servers
    that would otherwise be idle -- which is what tells the activated-server law
    apart from a plain multiserver station.
    """
    for n, want in [([0.2, 0.0], 3.0), ([1.0, 0.0], 3.0), ([1.0, 1.0], 1.5),
                    ([2.0, 1.0], 1.0), ([3.0, 0.0], 1.0), ([5.0, 5.0], 1.0)]:
        assert sn_compat_scaling([[1, 1]], [3], [1], n) == pytest.approx(want), n


def test_full_pool_equals_its_own_declared_scaling():
    """
    A fully-compatible pool of m servers lowers to the scaling its law declares.

    The pool machinery must add nothing beyond that scaling, which is what
    comparing it against the SAME eta declared by hand pins.

    IT IS NOT THE NEUTRAL eta == 1, and until 2026-08-28 it was: the denominator
    damped by min(1, N) rather than min(1, N/S), which cancelled the redundancy
    speed-up and collapsed a full pool onto the plain multiserver station. This
    test asserted that collapse. The second half below is what would have caught
    the change had it been written this way first.
    """
    np.testing.assert_allclose(_qlen('full'), _qlen('etafull'), rtol=1e-12)
    assert not np.allclose(_qlen('full'), _qlen('eta1'), atol=1e-6)


def test_compatibility_structure_changes_the_answer():
    """
    The compatibility graph must not be invisible.

    Under the graph each task reaches only 2 of the 3 servers, so the station
    clears less than the fully-compatible pool and throughput falls.
    """
    full, compat = _qlen('full'), _qlen('compat')
    assert not np.allclose(full, compat, atol=1e-6)


def test_srvn_ph_refuses_a_pool():
    """The composed phase-type law cannot carry a station rate law, so it says so."""
    with pytest.raises(Exception) as ei:
        LN(_build('compat'), lambda mm: MVA(mm), method='srvn.ph').getAvgTable()
    assert 'srvn.cs' in str(ei.value)


def test_pool_and_joint_dependence_are_exclusive():
    """Two rate laws on one server is a modelling error, not a composition."""
    model = LayeredNetwork('X')
    P1 = Processor(model, 'P1', 2, SchedStrategy.PS)
    T2 = Task(model, 'T2', 2, SchedStrategy.FCFS).on(P1)
    Entry(model, 'E2').on(T2)
    P1.addServerType(ServerType('All', 2, [T2]))
    with pytest.raises(ValueError, match='server pools'):
        P1.setJointDependence(lambda n: 1.0, np.array([1.0]))


def test_pool_sizes_must_match_the_multiplicity():
    """The pools partition the declared servers; a mismatch is refused."""
    with pytest.raises(ValueError, match='multiplicity'):
        _build('badmult').getStruct()


def test_every_operand_must_be_servable():
    """An operand no pool can serve would never complete, so it is refused."""
    with pytest.raises(ValueError, match='compatible with no server pool'):
        _build('unserved').getStruct()
