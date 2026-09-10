#!/usr/bin/env python3
"""Method 'srvn.ph': the activity graph of an entry as a phase-type server law.

The default LN method turns every activity graph into routing: one class per
task, entry, activity and call, plus Fork, Join, Router and ClassSwitch nodes.
Method 'srvn.ph' composes each entry graph into a single phase-type law by the
exact series-parallel reduction of Workflow.toPH, so a layer becomes a
two-station cycle, Delay('Clients') + Queue(server), with one closed class per
caller task. The sequencing survives as a distribution rather than as routing.

The model below exercises the two constructs the reduction handles exactly and
the default routing encoding only approximates: an AND fork-join, whose branch
times are a maximum and not a sum, and a geometric loop.
"""

import time

from line_solver import (Activity, ActivityPrecedence, Entry, Exp,
                         GlobalConstants, LayeredNetwork, LN, MVA, Processor,
                         SchedStrategy, Task, VerboseLevel)
from line_solver.solvers import SolverLNOptions


def lqn_srvnph():
    model = LayeredNetwork('srvnphExample')

    P1 = Processor(model, 'P1', 1, SchedStrategy.INF)
    T1 = Task(model, 'T1', 20, SchedStrategy.REF).on(P1)
    T1.setThinkTime(Exp.fitMean(1.0))
    E1 = Entry(model, 'E1').on(T1)

    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T2 = Task(model, 'T2', 5, SchedStrategy.FCFS).on(P2)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T2)

    P3 = Processor(model, 'P3', 1, SchedStrategy.PS)
    T3 = Task(model, 'T3', 3, SchedStrategy.FCFS).on(P3)
    E4 = Entry(model, 'E4').on(T3)

    # client
    Activity(model, 'A1', Exp.fitMean(0.1)).on(T1).boundTo(E1).synchCall(E2, 1).synchCall(E3, 1)

    # E2: AND fork-join over two branches, one of which calls a third tier
    A20 = Activity(model, 'A20', Exp.fitMean(0.2)).on(T2).boundTo(E2)
    A21 = Activity(model, 'A21', Exp.fitMean(0.3)).on(T2)
    A22 = Activity(model, 'A22', Exp.fitMean(0.2)).on(T2).synchCall(E4, 1)
    A23 = Activity(model, 'A23', Exp.fitMean(0.1)).on(T2).repliesTo(E2)
    T2.addPrecedence(ActivityPrecedence.AndFork(A20, [A21, A22]))
    T2.addPrecedence(ActivityPrecedence.AndJoin([A21, A22], A23))

    # E3: a geometric loop of mean 3 around a body of two activities
    A30 = Activity(model, 'A30', Exp.fitMean(0.1)).on(T2).boundTo(E3)
    A31 = Activity(model, 'A31', Exp.fitMean(0.2)).on(T2)
    A32 = Activity(model, 'A32', Exp.fitMean(0.1)).on(T2)
    A33 = Activity(model, 'A33', Exp.fitMean(0.1)).on(T2).repliesTo(E3)
    T2.addPrecedence(ActivityPrecedence.Loop(A30, [A31, A32, A33], 3))

    # third tier
    Activity(model, 'A4', Exp.fitMean(0.15)).on(T3).boundTo(E4).repliesTo(E4)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.SILENT)
    model = lqn_srvnph()

    lnoptions = SolverLNOptions()
    lnoptions.verbose = False

    t0 = time.time()
    solver_default = LN(model, lambda l: MVA(l, verbose=False), lnoptions)
    table_default = solver_default.getAvgTable()
    t_default = time.time() - t0
    print('\nLN(default) Results [%.3f s]:' % t_default)
    print(table_default)

    lnoptions2 = SolverLNOptions()
    lnoptions2.verbose = False
    # the method name; bare 'srvnph' is unrecognised and falls back to srvn.cs
    lnoptions2.method = 'srvn.ph'
    t0 = time.time()
    solver_srvnph = LN(model, lambda l: MVA(l, verbose=False), lnoptions2)
    table_srvnph = solver_srvnph.getAvgTable()
    t_srvnph = time.time() - t0
    print('\nLN(srvn.ph) Results [%.3f s]:' % t_srvnph)
    print(table_srvnph)

    nclasses_default = sum(m.getNumberOfClasses() for m in solver_default.ensemble)
    nclasses_srvnph = sum(m.getNumberOfClasses() for m in solver_srvnph.ensemble)
    print('\nLayer classes: default %d, srvn.ph %d. Runtime: %.3fs vs %.3fs (%.2fx).'
          % (nclasses_default, nclasses_srvnph, t_default, t_srvnph,
             t_default / t_srvnph if t_srvnph > 0 else float('inf')))
