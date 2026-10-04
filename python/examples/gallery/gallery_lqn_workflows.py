#!/usr/bin/env python3
"""Gallery Example: gallery_lqn_workflows"""

from line_solver import *

def gallery_lqn_workflows():
    """Layered network with loop, and-fork/join and or-fork/join precedence."""
    model = LayeredNetwork('LQN-Workflows')
    P1 = Processor(model, 'P1', float('inf'), SchedStrategy.INF)
    T1 = Task(model, 'T1', 1, SchedStrategy.REF).on(P1)
    T1.set_think_time(Immediate())
    E1 = Entry(model, 'Entry').on(T1)

    P2 = Processor(model, 'P2', float('inf'), SchedStrategy.INF)
    T2 = Task(model, 'T2', float('inf'), SchedStrategy.INF).on(P2).set_think_time(Immediate())
    E2 = Entry(model, 'E2').on(T2)

    P3 = Processor(model, 'P3', 5, SchedStrategy.PS)
    T3 = Task(model, 'T3', float('inf'), SchedStrategy.INF).on(P3)
    T3.set_think_time(Exp.fit_mean(10))
    E3 = Entry(model, 'E3').on(T3)

    A1 = Activity(model, 'A1', Exp.fit_mean(1)).on(T1).bound_to(E1)
    A2 = Activity(model, 'A2', Exp.fit_mean(2)).on(T1)
    A3 = Activity(model, 'A3', Exp.fit_mean(3)).on(T1).synch_call(E2)

    B1 = Activity(model, 'B1', Exp.fit_mean(0.1)).on(T2).bound_to(E2)
    B2 = Activity(model, 'B2', Exp.fit_mean(0.2)).on(T2)
    B3 = Activity(model, 'B3', Exp.fit_mean(0.3)).on(T2)
    B4 = Activity(model, 'B4', Exp.fit_mean(0.4)).on(T2)
    B5 = Activity(model, 'B5', Exp.fit_mean(0.5)).on(T2)
    B6 = Activity(model, 'B6', Exp.fit_mean(0.6)).on(T2).synch_call(E3).replies_to(E2)

    C1 = Activity(model, 'C1', Exp.fit_mean(0.1)).on(T3).bound_to(E3)
    C2 = Activity(model, 'C2', Exp.fit_mean(0.2)).on(T3)
    C3 = Activity(model, 'C3', Exp.fit_mean(0.3)).on(T3)
    C4 = Activity(model, 'C4', Exp.fit_mean(0.4)).on(T3)
    C5 = Activity(model, 'C5', Exp.fit_mean(0.5)).on(T3).replies_to(E3)

    T1.add_precedence(ActivityPrecedence.Loop(A1, [A2, A3], 3))
    T2.add_precedence(ActivityPrecedence.Serial(B4, B5))
    T2.add_precedence(ActivityPrecedence.AndFork(B1, [B2, B3, B4]))
    T2.add_precedence(ActivityPrecedence.AndJoin([B2, B3, B5], B6))
    T3.add_precedence(ActivityPrecedence.OrFork(C1, [C2, C3, C4], [0.3, 0.3, 0.4]))
    T3.add_precedence(ActivityPrecedence.OrJoin([C2, C3, C4], C5))
    return model



if __name__ == '__main__':
    model = gallery_lqn_workflows()
    print('Model built:', type(model).__name__)
