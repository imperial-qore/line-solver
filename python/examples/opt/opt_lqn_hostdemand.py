"""
Tutorial 9: LayeredNetwork (LQN) host-demand tuning with layer freezing.

Optimize the mean host demands of two activities in a 3-task client/server
LQN to minimize the reference task's end-to-end response time, subject to a
utilization cap on the shared processor P1. Demonstrates:

  * passing a LayeredNetwork to OptimizationProblem (solved with SolverLN),
  * LQN decision variables (HostDemand) and node-keyed objectives/constraints,
  * the three LQN gradient sources (lqn_gradient),
  * explicit layer freezing (frozen_layers) and adaptive auto_freeze.

Minimizing latency drives both tunable demands to their lower bounds.
"""

from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                         Exp, SchedStrategy)
from line_solver import (OptimizationProblem, HostDemand,
                         MinimizeSystemResponseTime, UtilizationConstraint)


def build():
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


def make_problem(model):
    p = OptimizationProblem(model)
    p.add_variable(HostDemand('AS1', bounds=(0.02, 0.2)))
    p.add_variable(HostDemand('AS2', bounds=(0.01, 0.1)))
    p.set_objective(MinimizeSystemResponseTime(
        'T1', subject_to=[UtilizationConstraint('P1', max_value=0.95)]))
    return p


if __name__ == '__main__':
    for mode in ('fd', 'partial_sens', 'partial_plus_fd'):
        res = make_problem(build()).solve(
            optimizer='gradient', lqn_gradient=mode, seed=1,
            max_iterations=8, gradient_restarts=1, time_limit=180)
        vals = {k: round(v, 4) for k, v in res.variable_values.items()}
        print(f"[{mode:16s}] obj={res.objective_value:.5f} "
              f"feasible={res.feasible} evals={res.model_evaluations} {vals}")

    # Explicit freeze: hold the P2 layer (AS3 lives there) fixed.
    res = make_problem(build()).solve(
        optimizer='gradient', lqn_gradient='fd', frozen_layers=['P2'],
        seed=1, max_iterations=8, gradient_restarts=1, time_limit=180)
    print(f"[frozen_layers=P2 ] obj={res.objective_value:.5f} "
          f"free={list(res.variable_values)}")

    # Adaptive auto-freeze via the layer-wise decomposition.
    wf = make_problem(build()).decompose()
    wf.setSolverOptions(optimizer='gradient', lqn_gradient='fd',
                        max_iterations=4, gradient_restarts=1,
                        time_limit=60, seed=1)
    r = wf.solveLayered(max_cycles=6, tolerance=1e-3, auto_freeze=True,
                        freeze_tol=1e-2)
    print(f"[auto_freeze      ] obj={r.final_objective:.5f} "
          f"cycles={r.cycles_completed} frozen={r.frozen_layers} "
          f"evals={r.model_evaluations}")
