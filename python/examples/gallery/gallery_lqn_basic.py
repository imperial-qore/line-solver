#!/usr/bin/env python3
"""Gallery Example: gallery_lqn_basic"""

from line_solver import *

def gallery_lqn_basic():
    """Basic 3-task layered queueing network (client/server chain)."""
    model = LayeredNetwork('LQN-Basic')
    P1 = Processor(model, 'P1', 2, SchedStrategy.PS)
    P2 = Processor(model, 'P2', 3, SchedStrategy.PS)
    T1 = Task(model, 'T1', 50, SchedStrategy.REF).on(P1).set_think_time(Exp(1 / 2))
    T2 = Task(model, 'T2', 50, SchedStrategy.FCFS).on(P1).set_think_time(Exp(1 / 3))
    T3 = Task(model, 'T3', 25, SchedStrategy.FCFS).on(P2).set_think_time(Exp(1 / 4))
    E1 = Entry(model, 'E1').on(T1)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T3)
    A1 = Activity(model, 'AS1', Exp(10)).on(T1).bound_to(E1).synch_call(E2, 1)
    A2 = Activity(model, 'AS2', Exp(20)).on(T2).bound_to(E2).synch_call(E3, 5).replies_to(E2)
    A3 = Activity(model, 'AS3', Exp(50)).on(T3).bound_to(E3).replies_to(E3)
    return model



if __name__ == '__main__':
    model = gallery_lqn_basic()
    print('Model built:', type(model).__name__)
