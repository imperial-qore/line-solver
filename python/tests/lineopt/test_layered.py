"""
Tests for LayeredNetwork (LQN) support in line_solver.opt.

Covers model-type detection, LQN decision variables (apply / getLayer /
currentValue), problem validation, the SolverLN evaluation path, the three
gradient modes, and layer freezing (explicit + adaptive). The solve-based tests
are marked slow because each LQN evaluation runs a full SolverLN fixed point.
"""

import numpy as np
import pytest

from line_solver.opt.layered import (
    is_layered, elem_name, resolve_activity, resolve_task, resolve_processor,
    activity_processor_name, dist_mean,
)
from line_solver.opt.variables import (
    HostDemand, TaskThinkTime, ActivityThinkTime, TaskMultiplicity,
    TaskReplication, ProcessorMultiplicity, ServiceRate,
)
from line_solver.opt.problem import OptimizationProblem


def _build_lqn():
    """3-task client/server LQN (gallery_lqn_basic)."""
    from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                             Exp, SchedStrategy)
    m = LayeredNetwork('LQN-Basic')
    P1 = Processor(m, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(m, 'P2', 3, SchedStrategy.PS)
    T1 = Task(m, 'T1', 50, SchedStrategy.REF).on(P1).set_think_time(Exp(1 / 2))
    T2 = Task(m, 'T2', 50, SchedStrategy.FCFS).on(P1).set_think_time(Exp(1 / 3))
    T3 = Task(m, 'T3', 25, SchedStrategy.FCFS).on(P2).set_think_time(Exp(1 / 4))
    E1 = Entry(m, 'E1').on(T1)
    E2 = Entry(m, 'E2').on(T2)
    E3 = Entry(m, 'E3').on(T3)
    Activity(m, 'AS1', Exp(10)).on(T1).bound_to(E1).synch_call(E2, 1)
    Activity(m, 'AS2', Exp(20)).on(T2).bound_to(E2).synch_call(E3, 5).replies_to(E2)
    Activity(m, 'AS3', Exp(50)).on(T3).bound_to(E3).replies_to(E3)
    return m


def _build_flat():
    from line_solver import (Network, Source, Queue, Sink, OpenClass, Exp,
                             SchedStrategy)
    m = Network('OQN')
    s = Source(m, 'Source')
    q = Queue(m, 'Queue', SchedStrategy.FCFS)
    k = Sink(m, 'Sink')
    c = OpenClass(m, 'C')
    s.setArrival(c, Exp(1.0))
    q.setService(c, Exp(2.0))
    m.link(m.serialRouting(s, q, k))
    return m


# ---- detection and element resolution ---------------------------------------

def test_is_layered_true_for_lqn():
    assert is_layered(_build_lqn()) is True


def test_is_layered_false_for_flat():
    assert is_layered(_build_flat()) is False


def test_element_resolution_and_linkage():
    m = _build_lqn()
    assert elem_name(resolve_activity(m, 'AS1')) == 'AS1'
    assert elem_name(resolve_task(m, 'T2')) == 'T2'
    assert elem_name(resolve_processor(m, 'P2')) == 'P2'
    # AS1 -> T1 -> P1; AS3 -> T3 -> P2
    assert activity_processor_name(m, 'AS1') == 'P1'
    assert activity_processor_name(m, 'AS3') == 'P2'


def test_dist_mean_handles_scalar_dist_and_none():
    from line_solver import Exp
    assert dist_mean(None) is None
    assert dist_mean(2.0) == 2.0
    assert abs(dist_mean(Exp(4.0)) - 0.25) < 1e-9


# ---- decision variable apply / getLayer / currentValue ----------------------

def test_hostdemand_apply_and_current_value():
    m = _build_lqn()
    var = HostDemand('AS1', bounds=(0.02, 0.2))
    # decode midpoint, apply, read back
    d = var.decode(np.array([0.5]))
    assert abs(d - 0.11) < 1e-9
    var.apply(m, d)
    assert abs(resolve_activity(m, 'AS1').getHostDemandMean() - 0.11) < 1e-6
    assert abs(var.currentValue(m) - 0.11) < 1e-6
    # host layer is the processor the activity runs on
    assert var.getLayer(m) == ['P1']
    # sensitivity row key: (processor, activity)
    assert var.sensKey(m) == ('P1', 'AS1')


def test_multiplicity_and_replication_apply():
    m = _build_lqn()
    tm = TaskMultiplicity('T2', bounds=(1, 100))
    tm.apply(m, 40)
    assert resolve_task(m, 'T2').multiplicity == 40
    assert tm.currentValue(m) == 40
    pm = ProcessorMultiplicity('P1', bounds=(1, 8))
    pm.apply(m, 5)
    assert resolve_processor(m, 'P1').multiplicity == 5
    assert pm.getLayer(m) == ['P1']
    tr = TaskReplication('T3', bounds=(1, 4))
    tr.apply(m, 3)
    assert tr.currentValue(m) == 3


def test_thinktime_apply():
    m = _build_lqn()
    tt = TaskThinkTime('T1', bounds=(0.1, 2.0))
    tt.apply(m, 1.5)
    assert abs(dist_mean(resolve_task(m, 'T1').think_time) - 1.5) < 1e-6
    at = ActivityThinkTime('AS2', bounds=(0.0, 1.0))
    at.apply(m, 0.3)
    assert abs(dist_mean(resolve_activity(m, 'AS2').think_time) - 0.3) < 1e-6


# ---- problem validation -----------------------------------------------------

def test_validation_rejects_flat_var_on_lqn():
    from line_solver import Queue, SchedStrategy
    m = _build_lqn()
    p = OptimizationProblem(m)
    assert p.isLayered() is True
    # a flat ServiceRate variable on an LQN model must be flagged
    fake_station = Queue.__new__(Queue)
    fake_station._name = 'P1'
    p._variables.append(ServiceRate.__new__(ServiceRate))
    p._variables[-1]._station = fake_station
    p._variables[-1]._jobclass = fake_station
    p._variables[-1]._name = 'bad'
    from line_solver.opt.objectives import MinimizeSystemResponseTime
    p.setObjective(MinimizeSystemResponseTime('T1'))
    errors = p.validate()
    assert any('flat-network variable' in e for e in errors)


def test_validation_rejects_lqn_var_on_flat():
    m = _build_flat()
    p = OptimizationProblem(m)
    assert p.isLayered() is False
    p.addVariable(HostDemand('AS1', bounds=(0.1, 0.2)))
    from line_solver.opt.objectives import MinimizeSystemResponseTime
    p.setObjective(MinimizeSystemResponseTime())
    errors = p.validate()
    assert any('LayeredNetwork variable' in e for e in errors)


# ---- evaluation (SolverLN path) ---------------------------------------------

@pytest.mark.slow
def test_lqn_evaluation_populates_node_metrics():
    from line_solver.opt.evaluator import LineEvaluator
    m = _build_lqn()
    ev = LineEvaluator(m, [HostDemand('AS1', bounds=(0.02, 0.2))])
    assert ev.is_layered is True
    res = ev.evaluateValues({'AS1_hostdemand': 0.1})
    assert res.feasible
    # per-node metrics keyed by LQN node name
    assert res.getUtilization('P1') > 0
    assert res.getResponseTime('E1') > 0
    # reference-task system metrics derived
    assert res.getSystemThroughput('T1') > 0


# ---- gradient modes ---------------------------------------------------------

@pytest.mark.slow
@pytest.mark.parametrize('mode', ['fd', 'partial_sens', 'partial_plus_fd'])
def test_gradient_modes_reduce_objective(mode):
    m = _build_lqn()
    p = OptimizationProblem(m)
    p.addVariable(HostDemand('AS1', bounds=(0.02, 0.2)))
    p.addVariable(HostDemand('AS2', bounds=(0.01, 0.1)))
    from line_solver.opt.objectives import (MinimizeSystemResponseTime,
                                            UtilizationConstraint)
    p.setObjective(MinimizeSystemResponseTime(
        'T1', subject_to=[UtilizationConstraint('P1', max_value=0.95)]))
    # The stopping rule that defines this test is max_iterations/convergence, NOT
    # the clock: the assertions below are about where the optimizer LANDS. At the
    # former time_limit=180 the clock was the binding rule for 'fd', which needs
    # two extra solves per variable per iteration and converges in 183.7s on an
    # idle host -- 2% OVER its own budget. So the test passed or failed on host
    # speed alone, and it duly failed on a loaded cluster node while passing
    # locally. The limit is now a safety valve with real headroom (measured 184s,
    # allowed 900s), which does not slow the passing case down: 'fd' stops itself
    # at convergence and the other two modes at ~114s.
    res = p.solve(optimizer='gradient', lqn_gradient=mode,
                  max_iterations=8, gradient_restarts=1, seed=1, time_limit=900)
    # If the clock ever becomes binding again, say so here rather than let it
    # surface as an unexplained variable value below.
    assert res.terminated_by != 'time_limit', (
        'the %s gradient run was cut off by its time limit, so the assertions '
        'below would be testing host speed, not the optimizer' % mode)
    assert res.feasible
    # minimizing latency drives both demands to their lower bounds
    assert res.variable_values['AS1_hostdemand'] < 0.1
    assert res.variable_values['AS2_hostdemand'] < 0.05


# ---- layer freezing ---------------------------------------------------------

@pytest.mark.slow
def test_explicit_layer_freeze_holds_variable():
    m = _build_lqn()
    p = OptimizationProblem(m)
    p.addVariable(HostDemand('AS1', bounds=(0.02, 0.2)))
    p.addVariable(HostDemand('AS3', bounds=(0.005, 0.05)))  # lives on P2
    from line_solver.opt.objectives import MinimizeSystemResponseTime
    p.setObjective(MinimizeSystemResponseTime('T1'))
    res = p.solve(optimizer='gradient', lqn_gradient='fd',
                  frozen_layers=['P2'], max_iterations=6,
                  gradient_restarts=1, seed=1, time_limit=180)
    # AS3 is on frozen layer P2 -> not optimized (absent from free results)
    assert 'AS3_hostdemand' not in res.variable_values
    assert 'AS1_hostdemand' in res.variable_values


@pytest.mark.slow
def test_adaptive_auto_freeze_converges():
    m = _build_lqn()
    p = OptimizationProblem(m)
    p.addVariable(HostDemand('AS1', bounds=(0.02, 0.2)))
    p.addVariable(HostDemand('AS2', bounds=(0.01, 0.1)))
    p.addVariable(HostDemand('AS3', bounds=(0.005, 0.05)))
    from line_solver.opt.objectives import (MinimizeSystemResponseTime,
                                            UtilizationConstraint)
    p.setObjective(MinimizeSystemResponseTime(
        'T1', subject_to=[UtilizationConstraint('P1', max_value=0.95)]))
    wf = p.decompose()
    # Bound the per-layer subproblem solvers so the test stays tractable.
    wf.setSolverOptions(optimizer='gradient', lqn_gradient='fd',
                        max_iterations=3, gradient_restarts=1,
                        time_limit=60, seed=1)
    r = wf.solveLayered(max_cycles=4, tolerance=1e-3, auto_freeze=True,
                        freeze_tol=1e-2)
    assert np.isfinite(r.final_objective)
    assert hasattr(r, 'frozen_layers')
    assert hasattr(r, 'model_evaluations')
    # auto-freeze must actually freeze converged layers on this model
    assert len(r.frozen_layers) >= 1
