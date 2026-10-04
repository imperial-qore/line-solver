"""
Validates the LQN conservation-law enumerator (line_solver.api.lqn.balance_equations).

Two things are checked, and they are different in kind.

STRUCTURE, on a two-tier client-server LQN with no solution supplied: the right
relations are emitted with the right term sets, branches and constants. That part
depends on the model only.

CONSISTENCY, on a converged SolverLN: every relation the enumerator emits is a law
the solution must satisfy, so each residual must vanish to the solver's own
tolerance. This is the real content of the test -- it turns the enumerator into a
conservation check on SolverLN itself. A residual that does NOT vanish is either a
documented convention difference (mult against the njobs SolverLN iterates on) or a
defect in the solution, and the models used here have neither.

Twin of the MATLAB test in line-test.git test/testsLN/test_lqn_balance_equations.m.
"""

import numpy as np
import pytest

from line_solver import (Activity, Entry, Exp, Immediate, LayeredNetwork, Processor,
                         SchedStrategy, SolverLN, Task)
from line_solver.api.lqn import balance_equations

TOL = 1e-6


def _client_server():
    """T1 (ref, 50 threads, think 2) -> T2 (50 threads) -> T3 (25 threads)."""
    m = LayeredNetwork('cs')
    P1 = Processor(m, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 3, SchedStrategy.PS)
    T1 = Task(m, 'T1', 50, SchedStrategy.REF).on(P1).setThinkTime(Exp(1.0 / 2))
    T2 = Task(m, 'T2', 50, SchedStrategy.FCFS).on(P1).setThinkTime(Exp(1.0 / 3))
    T3 = Task(m, 'T3', 25, SchedStrategy.FCFS).on(P2).setThinkTime(Exp(1.0 / 4))
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    E3 = Entry(m, 'E3').on(T3)
    Activity(m, 'AS1', Exp(10)).on(T1).boundTo(E1).synchCall(E2, 1)
    Activity(m, 'AS2', Exp(20)).on(T2).boundTo(E2).synchCall(E3, 5).repliesTo(E2)
    Activity(m, 'AS3', Exp(50)).on(T3).boundTo(E3).repliesTo(E3)
    return m


def _serial():
    """Two tasks with a serial activity pair each, so visit counts are exercised."""
    from line_solver import ActivityPrecedence
    m = LayeredNetwork('serial')
    P1 = Processor(m, 'P1', 1, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 1, SchedStrategy.PS)
    T1 = Task(m, 'T1', 10, SchedStrategy.REF).on(P1).setThinkTime(Exp.fitMean(100))
    T2 = Task(m, 'T2', 1, SchedStrategy.FCFS).on(P2).setThinkTime(Immediate())
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    A1 = Activity(m, 'AS1', Exp.fitMean(1.6)).on(T1).boundTo(E1)
    A2 = Activity(m, 'AS2', Immediate()).on(T1).synchCall(E2, 1)
    A3 = Activity(m, 'AS3', Exp.fitMean(5)).on(T2).boundTo(E2)
    A4 = Activity(m, 'AS4', Exp.fitMean(1)).on(T2).repliesTo(E2)
    T1.addPrecedence(ActivityPrecedence.Serial([A1, A2]))
    T2.addPrecedence(ActivityPrecedence.Serial([A3, A4]))
    return m


def _two_callers():
    """
    Two reference tasks with DIFFERENT call counts into one server entry, kept well
    below saturation. This is the case the per-class decomposition exists for: the
    server's busy threads split over two caller classes that must be told apart.
    """
    m = LayeredNetwork('twocallers')
    P1 = Processor(m, 'P1', 4, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 2, SchedStrategy.PS)
    C1 = Task(m, 'C1', 3, SchedStrategy.REF).on(P1).setThinkTime(Exp.fitMean(1.0))
    C2 = Task(m, 'C2', 4, SchedStrategy.REF).on(P1).setThinkTime(Exp.fitMean(2.0))
    S = Task(m, 'S', 6, SchedStrategy.FCFS).on(P2)
    E1 = Entry(m, 'E1').on(C1)
    E2 = Entry(m, 'E2').on(C2)
    ES = Entry(m, 'ES').on(S)
    Activity(m, 'A1', Exp.fitMean(0.3)).on(C1).boundTo(E1).synchCall(ES, 1)
    Activity(m, 'A2', Exp.fitMean(0.5)).on(C2).boundTo(E2).synchCall(ES, 2)
    Activity(m, 'AS', Exp.fitMean(0.05)).on(S).boundTo(ES).repliesTo(ES)
    return m


def _by(out, kind, name):
    for r in out.eqs:
        if r.kind == kind and r.targetname == name:
            return r
    raise AssertionError('no %s relation on %s' % (kind, name))


# ----------------------------------------------------------------------
# structure
# ----------------------------------------------------------------------

def test_families_and_counts():
    out = balance_equations(_client_server().getStruct())
    kinds = [r.kind for r in out.eqs]
    assert kinds.count('little') == 3        # one per task
    assert kinds.count('callflow') == 2      # one per call
    assert kinds.count('entryflow') == 2     # E2, E3; E1 is driven by the ref cycle
    assert kinds.count('actflow') == 3       # one per activity
    assert kinds.count('hostutil') == 2      # one per processor


def test_little_branches_and_constants():
    lqn = _client_server().getStruct()
    out = balance_equations(lqn)
    r1 = _by(out, 'little', 'T1')
    assert r1.branch == 'ref'
    assert r1.const == 50
    assert r1.termisentry == [True]          # its own entry drives its cycle
    r2 = _by(out, 'little', 'T2')
    assert r2.branch == 'queueing'
    assert r2.const == 50
    assert r2.termisentry == [False]         # one call class
    assert r2.scaled                         # U is normalized to [0,1] here
    r3 = _by(out, 'little', 'T3')
    assert r3.branch == 'queueing'
    assert r3.const == 25


def test_hostutil_coefficients_are_the_host_demands():
    out = balance_equations(_client_server().getStruct())
    r = _by(out, 'hostutil', 'P1')
    assert r.const == 2                      # the declared multiplicity of P1
    assert sorted(np.round(r.coeff, 10)) == [0.05, 0.1]
    r2 = _by(out, 'hostutil', 'P2')
    assert r2.const == 3
    assert np.allclose(r2.coeff, [0.02])


def test_callflow_carries_the_mean_call_count():
    out = balance_equations(_client_server().getStruct())
    ys = dict((r.targetname, r.coeff[0]) for r in out.eqs if r.kind == 'callflow')
    assert ys['AS1=>E2'] == pytest.approx(1.0)
    assert ys['AS2=>E3'] == pytest.approx(5.0)


def test_visits_are_one_on_a_serial_chain():
    lqn = _serial().getStruct()
    out = balance_equations(lqn)
    for r in out.eqs:
        if r.kind == 'actflow':
            assert r.coeff[0] == pytest.approx(1.0)


def test_incidence_matrices_have_the_call_classes_only():
    lqn = _client_server().getStruct()
    out = balance_equations(lqn)
    # rows are tasks and columns calls, both 0-based over their own band, so
    # T2 is row 1 and the call into it column 0
    assert out.A_little[1, 0] == 1.0
    assert out.A_little[2, 1] == 1.0
    assert out.A_little[0].sum() == 0.0      # T1 is entry-driven, no call class
    assert out.A_flow.shape == (lqn.ncalls, lqn.nidx)


def test_symbolic_report_is_printable_without_a_solution():
    out = balance_equations(_client_server().getStruct())
    txt = str(out)
    assert 'thread-pool Little' in txt
    assert 'host utilization law' in txt
    assert np.isnan(out.maxresidual)


# ----------------------------------------------------------------------
# consistency against a converged SolverLN
# ----------------------------------------------------------------------

def test_two_callers_split_the_server_thread_pool_per_class():
    """
    A server called by two tasks gets one class per CALL, not one per task or one in
    total, and the classes carry different utilizations because the call counts
    differ. The thread-pool law must still close over their sum.
    """
    model = _two_callers()
    lqn = model.getStruct()
    out = balance_equations(lqn)
    r = _by(out, 'little', 'S')
    assert len(r.terms) == 2
    assert r.termisentry == [False, False]        # both are call classes
    assert r.const == 6

    solver = SolverLN(model, verbose=False)
    QN, UN, RN, TN, AN, WN = solver.get_ensemble_avg()
    out = balance_equations(lqn, solver, UN)
    r = _by(out, 'little', 'S')
    assert len(r.perclassutil) == 2
    assert r.perclassutil[0] != pytest.approx(r.perclassutil[1], rel=1e-3)
    assert sum(r.perclassutil) < 1.0             # unsaturated, so the law is an equality
    assert not r.clamped
    assert abs(r.residual) <= TOL


@pytest.mark.parametrize('build', [_client_server, _serial, _two_callers],
                         ids=['clientserver', 'serial', 'twocallers'])
def test_every_relation_vanishes_on_a_converged_solve(build):
    model = build()
    lqn = model.getStruct()
    solver = SolverLN(model, verbose=False)
    QN, UN, RN, TN, AN, WN = solver.get_ensemble_avg()
    out = balance_equations(lqn, solver, UN)
    bad = ['%s %s: residual %.3e (lhs %.6g, rhs %.6g)'
           % (r.kind, r.targetname, r.residual, r.lhs, r.rhs)
           for r in out.eqs
           if not r.degenerate and not np.isnan(r.residual) and abs(r.residual) > TOL]
    assert not bad, 'conservation violated:\n' + '\n'.join(bad)
    assert out.maxresidual <= TOL


def test_per_class_utilization_reproduces_the_solver_iterate():
    """
    The per-call decomposition of the busy threads must add up to the utilization
    SolverLN itself carries for that task, which is what makes the emitted
    per-class U a refinement of the solver's aggregate rather than a new quantity.

    The two reach the same number by different routes through the same iterate --
    this one as sum_c X(c)*S(e(c))/mult, the solver's as the layer station's own UN
    -- so they agree to the fixed point's tolerance and not to machine precision.
    """
    model = _client_server()
    lqn = model.getStruct()
    solver = SolverLN(model, verbose=False)
    solver.get_ensemble_avg()
    out = balance_equations(lqn, solver)
    util = np.asarray(solver.util).ravel()
    for name in ('T2', 'T3'):
        r = _by(out, 'little', name)
        assert sum(r.perclassutil) == pytest.approx(util[r.target], rel=1e-6)


def test_hostutil_is_not_instantiated_without_the_reported_utilization():
    """
    The host law needs the REPORTED utilization, which the `util` iterate does not
    carry. Without it the record must stay symbolic rather than report a residual
    against a zero it never meant.
    """
    model = _client_server()
    lqn = model.getStruct()
    solver = SolverLN(model, verbose=False)
    solver.get_ensemble_avg()
    out = balance_equations(lqn, solver)
    for r in out.eqs:
        if r.kind == 'hostutil':
            assert np.isnan(r.residual)
        elif not r.degenerate:
            assert not np.isnan(r.residual)
