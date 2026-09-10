"""Every LQN solver must agree about WHICH CELLS OF THE TABLE EXIST.

A NaN in an LQN average table is not a failed computation: it says the quantity
is not defined for that element kind. A processor has no queue length, response
time, residence or throughput of its own; an entry has no residence; nothing
reports an arrival rate. Two solvers that disagree about the mask disagree about
the MODEL rather than about arithmetic, and a tolerance-based comparison cannot
see it -- `parity-static/_compare_vendor.py` passes any cell where either side
is NaN, which is why the goldens carry the string "NaN" beside the numbers.

Twin of `jar/src/test/java/jline/solvers/ln/LqnNaNMaskParityTest.java`; the C++
and MATLAB twins are `cpp/tests/test_lqn_nan_mask_parity.cpp` and
`line-test.git/test/testsMisc/test_lqn_nan_mask_parity.m`. LQNS is deliberately
NOT in the panel: it reports no residence time and no arrival rate at all, and
leaves both undefined by decision rather than by omission (see
`_kb/06-solver-catalog.md`, LQNS wrapper) -- and its binary is not always
present.
"""

import math
import os

import pytest

import line_solver
# The worktree copy of line_solver must be the one under test; a global install
# would silently validate the wrong code.
assert os.path.realpath(__file__).rsplit('/python/', 1)[0] in os.path.realpath(
    line_solver.__file__), (
    "test must import the worktree line_solver, got %s" % line_solver.__file__)

from line_solver import (Activity, Entry, Erlang, Exp, GlobalConstants, LayeredNetwork,
                         Processor, SchedStrategy, SolverLN, SolverMVA, SolverNC, Task,
                         VerboseLevel)

METRICS = ('QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput')


def _multi_solvers():
    """The model of `matlab/examples/basic/layeredModel/lqn_multi_solvers.m`."""
    m = LayeredNetwork('LQN1')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    t1 = Task(m, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(m, 'E1').on(t1)
    p2 = Processor(m, 'P2', 1, SchedStrategy.INF)
    t2 = Task(m, 'T2', 1, SchedStrategy.INF).on(p2)
    e2 = Entry(m, 'E2').on(t2)
    t1.set_think_time(Erlang.fit_mean_and_order(0.0001, 2))
    Activity(m, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(e2, 3)
    Activity(m, 'A2', Exp(1.0)).on(t2).bound_to(e2).replies_to(e2)
    return m


def _ignored_component():
    """The same model plus a component NOTHING REACHES: no reference task, no
    caller.

    SolverLN marks such a component `ignore` and reports zero for it -- which is
    right for the measures its elements HAVE, and was being written over the ones
    they do not, flattening the mask. The `lqn_ofbiz` golden carries exactly this
    shape (its USAGE_DELAY component), which is why the defect survived the
    2026-08-21 audit: no golden pairs LQNS with an LN row on a model that has one.
    """
    m = LayeredNetwork('LQNignore')
    p1 = Processor(m, 'P1', 1, SchedStrategy.INF)
    t1 = Task(m, 'T1', 1, SchedStrategy.REF).on(p1)
    e1 = Entry(m, 'E1').on(t1)
    p2 = Processor(m, 'P2', 1, SchedStrategy.INF)
    t2 = Task(m, 'T2', 1, SchedStrategy.INF).on(p2)
    e2 = Entry(m, 'E2').on(t2)
    t1.set_think_time(Exp(1.0))
    Activity(m, 'A1', Exp(1.0)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(m, 'A2', Exp(1.0)).on(t2).bound_to(e2).replies_to(e2)
    # unreachable: T3 is not a reference task and nobody calls E3
    p3 = Processor(m, 'P3', 1, SchedStrategy.INF)
    t3 = Task(m, 'T3', 1, SchedStrategy.FCFS).on(p3)
    e3 = Entry(m, 'E3').on(t3)
    Activity(m, 'A3', Exp(1.0)).on(t3).bound_to(e3).replies_to(e3)
    return m


def _mask_of(table):
    """`{"P1.QLen": True, "P1.Util": False, ...}` -- the mask, spelled so a
    failure reads."""
    mask = {}
    names = [str(v) for v in table['Node']]
    for metric in METRICS:
        col = list(table[metric])
        for i, name in enumerate(names):
            if i >= len(col):
                break
            try:
                mask['%s.%s' % (name, metric)] = math.isnan(float(col[i]))
            except (TypeError, ValueError):
                mask['%s.%s' % (name, metric)] = False
    return mask


def _ln_table(model, layer=None, method=None):
    GlobalConstants.set_verbose(VerboseLevel.SILENT)
    opts = SolverLN.default_options()
    opts.verbose = 0
    if method is not None:
        opts.method = method
    if layer is None:
        return SolverLN(model, opts).avg_table()
    return SolverLN(model, layer, opts).avg_table()


def _mva_layer(net):
    o = SolverMVA.default_options()
    o.verbose = 0
    return SolverMVA(net, o)


def _nc_layer(net):
    o = SolverNC.default_options()
    o.verbose = 0
    return SolverNC(net, o)


def _assert_same_mask(ref_name, ref, other_name, other):
    disagree = []
    for key, want in ref.items():
        got = other.get(key)
        if got is None:
            continue  # a solver that omits a column is not a mask defect
        if got != want:
            disagree.append('%s (%s=%s, %s=%s)' % (
                key, ref_name, 'NaN' if want else 'value',
                other_name, 'NaN' if got else 'value'))
    assert not disagree, ('LQN table NaN mask differs between %s and %s: %s'
                          % (ref_name, other_name, '; '.join(disagree)))


def test_every_solver_agrees_on_which_cells_exist():
    mva = _mask_of(_ln_table(_multi_solvers(), _mva_layer))
    nc = _mask_of(_ln_table(_multi_solvers(), _nc_layer))
    _assert_same_mask('LN(MVA)', mva, 'LN(NC)', nc)


def test_the_mask_is_the_documented_one():
    """The mask itself, pinned.

    Without this the test above would still pass if every solver regressed the
    same way.
    """
    mask = _mask_of(_ln_table(_multi_solvers(), _mva_layer))

    # A processor has no queue, no response time, no residence, no throughput of
    # its own -- only a utilization. TN is set NaN explicitly in getEnsembleAvg
    # "for consistency with LQNS".
    for p in ('P1', 'P2'):
        assert mask['%s.QLen' % p], '%s.QLen must be NaN' % p
        assert mask['%s.RespT' % p], '%s.RespT must be NaN' % p
        assert mask['%s.ResidT' % p], '%s.ResidT must be NaN' % p
        assert mask['%s.Tput' % p], '%s.Tput must be NaN' % p
        assert not mask['%s.Util' % p], '%s.Util must be a value' % p
    # An entry has a response time but no residence.
    for e in ('E1', 'E2'):
        assert mask['%s.ResidT' % e], '%s.ResidT must be NaN' % e
        assert not mask['%s.RespT' % e], '%s.RespT must be a value' % e
    # Tasks and activities carry both.
    for x in ('T1', 'T2', 'A1', 'A2'):
        assert not mask['%s.ResidT' % x], '%s.ResidT must be a value' % x
    # Nobody reports an arrival rate on an LQN.
    for x in ('P1', 'T1', 'E1', 'A1'):
        assert mask['%s.ArvR' % x], '%s.ArvR must be NaN' % x


def test_an_ignored_component_keeps_the_mask():
    """A component no reference task reaches is IDLE, not undefined.

    Its elements report zero for every measure their kind has, and NaN for the
    ones it does not -- the same mask a reachable element carries. Until
    2026-08-25 the ignore branch wrote a flat zero across all six columns, so an
    unreachable processor claimed a queue length of 0 and an arrival rate of 0
    where every solver, LQNS included, reports neither.
    """
    mask = _mask_of(_ln_table(_ignored_component(), _mva_layer))
    reachable = _mask_of(_ln_table(_multi_solvers(), _mva_layer))

    # the unreachable rows carry the SAME mask as the reachable ones
    for unreached, reached, metric in (
            ('P3', 'P1', 'QLen'), ('P3', 'P1', 'RespT'), ('P3', 'P1', 'ResidT'),
            ('P3', 'P1', 'Tput'), ('P3', 'P1', 'Util'),
            ('T3', 'T2', 'RespT'), ('T3', 'T2', 'ResidT'),
            ('E3', 'E2', 'RespT'), ('E3', 'E2', 'ResidT'),
            ('A3', 'A2', 'ResidT')):
        assert mask['%s.%s' % (unreached, metric)] == reachable['%s.%s' % (reached, metric)], (
            '%s.%s must carry the mask of %s.%s' % (unreached, metric, reached, metric))

    # spelled out, so the loop above cannot pass by both sides regressing
    assert mask['P3.QLen'], 'P3.QLen must be NaN'
    assert mask['P3.RespT'], 'P3.RespT must be NaN'
    assert mask['P3.Tput'], 'P3.Tput must be NaN'
    assert not mask['P3.Util'], 'P3.Util must be a value'
    assert mask['T3.RespT'], 'T3.RespT must be NaN'
    assert not mask['T3.ResidT'], 'T3.ResidT must be a value'
    assert mask['E3.ResidT'], 'E3.ResidT must be NaN'
    for x in ('P3', 'T3', 'E3', 'A3'):
        assert mask['%s.ArvR' % x], '%s.ArvR must be NaN' % x

    # and every layer engine still agrees about it
    nc = _mask_of(_ln_table(_ignored_component(), _nc_layer))
    _assert_same_mask('LN(MVA)', mask, 'LN(NC)', nc)


def test_the_ph_encoding_masks_an_ignored_component_the_same_way():
    """`srvn.ph` assembles the table in its own routine, and had its own copy of
    the flat-zero branch.

    The two encodings rebuild every figure differently -- one reads class rows
    off the ensemble, the other composes a phase-type law -- so agreeing on the
    mask is a claim about the table, not about shared code.
    """
    default = _mask_of(_ln_table(_ignored_component(), _mva_layer))
    ph = _mask_of(_ln_table(_ignored_component(), _mva_layer, method='srvn.ph'))
    _assert_same_mask('LN(MVA)', default, "LN(MVA, srvn.ph)", ph)
    assert ph['P3.QLen'], 'P3.QLen must be NaN under srvn.ph'
    assert ph['P3.ArvR'], 'P3.ArvR must be NaN under srvn.ph'
    assert not ph['P3.Util'], 'P3.Util must be a value under srvn.ph'
