"""
Regression tests for the convergence test of the SolverNC non-exponential
correction loop.

The loop residual mirrors MATLAB solver_nc.m line 66,

    while max(abs(1-eta./eta_1)) > options.iter_tol & it < options.iter_max

whose NaN semantics are load-bearing: a station whose eta is identically zero
in both iterates yields 0/0 = NaN, and MATLAB's max omits NaN, so that station
is excluded from the residual rather than holding the loop open. Porting this
with a regularized denominator, abs(1 - eta/(eta_1 + tol)), turns the NaN into
0 and hence a permanent residual of exactly 1, which never falls below
iter_tol: the loop then runs to iter_max on every call.

That failure mode is invisible to a value-based golden, because the extra
iterations are a fixed point and leave every reported metric unchanged. It is
only observable as work, so the integration test below counts calls to
npfqn_nonexp_approx instead of comparing numbers. Before the fix, lqn_ofbiz
took 405s instead of 3.3s with identical output.
"""
import os

import numpy as np
import pytest

from line_solver import (Activity, ClosedClass, Delay, Entry, Erlang, Exp,
                         GlobalConstants, LayeredNetwork, LN, NC, Network,
                         Processor, Queue, SchedStrategy, SolverNC, Task,
                         VerboseLevel)
from line_solver.api.solvers.nc.handler import _eta_residual


class TestEtaResidual:
    """Unit tests for the residual itself, against MATLAB's max/NaN semantics."""

    def test_all_zero_pair_is_excluded_not_scored_as_one(self):
        # The defect in one line: station 1 is zero in both iterates. MATLAB
        # scores it NaN and drops it, so the residual is station 0's alone.
        assert _eta_residual(np.array([0.5, 0.0]), np.array([0.5, 0.0])) == 0.0

    def test_converged_component_survives_a_zero_component(self):
        eta = np.array([0.8, 0.0, 0.4])
        eta_1 = np.array([0.8, 0.0, 0.2])
        # Only station 2 is unconverged: |1 - 0.4/0.2| = 1.0
        assert _eta_residual(eta, eta_1) == pytest.approx(1.0)

    def test_matches_matlab_formula_when_no_zeros_present(self):
        eta = np.array([0.5, 0.25, 1.0])
        eta_1 = np.array([0.4, 0.25, 0.8])
        expected = np.max(np.abs(1.0 - eta / eta_1))
        assert _eta_residual(eta, eta_1) == pytest.approx(expected)

    def test_all_nan_reports_zero_not_nan(self):
        # Nothing is left to converge; a NaN residual would compare False
        # against iter_tol and silently exit, so report 0 explicitly.
        assert _eta_residual(np.zeros(3), np.zeros(3)) == 0.0

    def test_first_iterate_against_zero_initializer_does_not_converge(self):
        # solver_nc seeds eta_1 = zeros(1,M), eta = ones(1,M); the loop must
        # enter, i.e. the residual must exceed any sane iter_tol.
        assert _eta_residual(np.ones(3), np.zeros(3)) > 1e-6

    def test_zero_denominator_with_nonzero_numerator_stays_infinite(self):
        # 1/0 = Inf in MATLAB too, and Inf is not NaN: this station is
        # genuinely unconverged and must keep the loop open.
        assert _eta_residual(np.array([0.3, 1.0]), np.array([0.3, 0.0])) == np.inf


def _count_nonexp_calls(monkeypatch, solve):
    """Run solve() and return how many times npfqn_nonexp_approx was called."""
    import line_solver.api.npfqn as npfqn_pkg
    import line_solver.api.npfqn.nonexp as nonexp_mod

    calls = {'n': 0, 'saw_zero_eta': False}
    original = nonexp_mod.npfqn_nonexp_approx

    def counting(*args, **kwargs):
        calls['n'] += 1
        result = original(*args, **kwargs)
        if np.any(np.asarray(result.eta) == 0.0):
            calls['saw_zero_eta'] = True
        return result

    # The handler imports the symbol inside the loop body, so both the package
    # re-export and the defining module have to be patched.
    monkeypatch.setattr(nonexp_mod, 'npfqn_nonexp_approx', counting)
    monkeypatch.setattr(npfqn_pkg, 'npfqn_nonexp_approx', counting)
    solve()
    return calls


def _mini_lqn():
    """Smallest LQN whose NC layers exhibit a permanently zero eta entry.

    A processor layer contains stations for tasks that do not call it in that
    layer; such a station has zero utilization, so eta = rho = 0 there on every
    iteration.
    """
    model = LayeredNetwork('mini')
    p1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    p2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    t1 = Task(model, 'T1', 5, SchedStrategy.REF).on(p1).set_think_time(Exp(1.0))
    t2 = Task(model, 'T2', 3, SchedStrategy.FCFS).on(p2)
    e1 = Entry(model, 'E1').on(t1)
    e2 = Entry(model, 'E2').on(t2)
    Activity(model, 'A1', Exp(2.0)).on(t1).bound_to(e1).synch_call(e2, 1)
    Activity(model, 'A2', Exp(3.0)).on(t2).bound_to(e2).replies_to(e2)
    return model


# The two integration tests below count calls to the native npfqn_nonexp_approx.
# Under lang='java' the solve is delegated to jline.jar over JSON, so the native
# correction loop is never entered and the count is identically zero: the loop
# they guard is the one in the JAR, covered on that side by its own tests.
_skip_java_dispatch = pytest.mark.skipif(
    os.environ.get('LINE_SOLVER_LANG') == 'java',
    reason="counts calls to the native npfqn_nonexp_approx; the java-dispatch "
           "path (solve delegated to jline.jar) never enters that loop")


class TestNonexpLoopIterationCount:
    @_skip_java_dispatch
    def test_lqn_layers_do_not_burn_iter_max(self, monkeypatch):
        GlobalConstants.set_verbose(VerboseLevel.SILENT)
        model = _mini_lqn()
        calls = _count_nonexp_calls(
            monkeypatch, lambda: LN(model, lambda x: NC(x)).get_avg_table())

        # Guards the premise: without a zero eta entry this model would not
        # exercise the defect at all and the count below would prove nothing.
        assert calls['saw_zero_eta'], (
            "no layer produced a zero eta entry; the model no longer exercises "
            "the NaN path and this test has stopped being a regression test")

        # Measured 160 with the correct residual and 1158 with the regularized
        # one. The threshold sits between the two with wide margin, so ordinary
        # drift in layer count or LN iterations cannot trip it.
        assert calls['n'] < 400, (
            "npfqn_nonexp_approx called %d times; the non-exponential "
            "correction loop is running to iter_max instead of converging"
            % calls['n'])

    @_skip_java_dispatch
    def test_closed_fcfs_network_converges_in_few_iterations(self, monkeypatch):
        GlobalConstants.set_verbose(VerboseLevel.SILENT)
        model = Network('nonexp')
        delay = Delay(model, 'D')
        queue = Queue(model, 'Q', SchedStrategy.FCFS)
        cls = ClosedClass(model, 'C', 4, delay)
        delay.set_service(cls, Exp(1.0))
        queue.set_service(cls, Erlang.fit_mean_and_scv(1.0, 0.5))
        model.link(Network.serial_routing(delay, queue))

        calls = _count_nonexp_calls(
            monkeypatch, lambda: SolverNC(model).avg_table())

        # Erlang service at an FCFS station does enter the correction loop, so
        # the count is nonzero, but it converges in a handful of iterations.
        assert calls['n'] > 0
        assert calls['n'] < 50
