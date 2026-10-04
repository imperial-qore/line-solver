#!/usr/bin/env python3
"""Gallery Example: gallery_multitier (Layered Queueing Network)"""

from line_solver import *

def gallery_multitier():
    """
    Create a 4-tier J2EE Layered Queueing Network (client/app/database).

    Reference LQN with multiple entries per task and or-fork/or-join
    activity precedence. Ported from MATLAB gallery_multitier.

    Returns:
        LayeredNetwork: 3-layer (client, application, database) LQN model.
    """
    model = LayeredNetwork('testLQN3')

    # Layer 1: client
    P0 = Processor(model, 'P0', 1, SchedStrategy.PS)
    T0 = Task(model, 'T0', 1, SchedStrategy.REF).on(P0)
    E0 = Entry(model, 'E0').on(T0)

    # Layer 2: application server
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1)
    E10 = Entry(model, 'E10').on(T1)
    E11 = Entry(model, 'E11').on(T1)
    E12 = Entry(model, 'E12').on(T1)
    E13 = Entry(model, 'E13').on(T1)

    # Layer 3: database
    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T2 = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P2)
    E20 = Entry(model, 'E20').on(T2)
    E21 = Entry(model, 'E21').on(T2)
    E22 = Entry(model, 'E22').on(T2)
    E23 = Entry(model, 'E23').on(T2)

    # Client activities
    A0 = Activity(model, 'A0', Exp(1.0)).on(T0).bound_to(E0).synch_call(E12, 1.0)
    A1 = Activity(model, 'A1', Exp(1.0)).on(T0).synch_call(E10, 1.0)
    A2 = Activity(model, 'A2', Exp(1.0)).on(T0).synch_call(E11, 1.0)
    A3 = Activity(model, 'A3', Exp(1.0)).on(T0).synch_call(E13, 1.0)

    # Application activities
    B0 = Activity(model, 'B0', Exp(1.0)).on(T1).bound_to(E10)
    B1 = Activity(model, 'B1', Exp(1.0)).on(T1).replies_to(E10)
    B2 = Activity(model, 'B2', Exp(1.0)).on(T1).bound_to(E11)
    B3 = Activity(model, 'B3', Exp(1.0)).on(T1).synch_call(E21, 1.0).replies_to(E11)
    B4 = Activity(model, 'B4', Exp(1.0)).on(T1).bound_to(E12).synch_call(E20, 1.0).replies_to(E12)
    B5 = Activity(model, 'B5', Exp(1.0)).on(T1).bound_to(E13)
    B6 = Activity(model, 'B6', Exp(1.0)).on(T1)
    B7 = Activity(model, 'B7', Exp(1.0)).on(T1).synch_call(E22, 1.0)
    B7a = Activity(model, 'B7a', Exp(1.0)).on(T1)
    B7b = Activity(model, 'B7b', Exp(1.0)).on(T1).synch_call(E23, 1.0)
    B8 = Activity(model, 'B8', Exp(1.0)).on(T1).replies_to(E13)

    # Database activities
    C0 = Activity(model, 'C0', Exp(1.0)).on(T2).bound_to(E20)
    C1 = Activity(model, 'C1', Exp(1.0)).on(T2).replies_to(E20)
    C2 = Activity(model, 'C2', Exp(1.0)).on(T2).bound_to(E21).replies_to(E21)
    C3 = Activity(model, 'C3', Exp(1.0)).on(T2).bound_to(E22)
    C4 = Activity(model, 'C4', Exp(1.0)).on(T2)
    C5 = Activity(model, 'C5', Exp(1.0)).on(T2).replies_to(E22)
    C6 = Activity(model, 'C6', Exp(1.0)).on(T2).bound_to(E23).replies_to(E23)

    # Precedences
    T0.add_precedence(ActivityPrecedence.Serial(A0, A1, A2, A3))
    T1.add_precedence(ActivityPrecedence.Serial(B0, B1))
    T1.add_precedence(ActivityPrecedence.Serial(B2, B3))
    T1.add_precedence(ActivityPrecedence.Serial(B5, B6, B7))
    T1.add_precedence(ActivityPrecedence.OrFork(B7, [B7a, B7b], [0.7, 0.3]))
    T1.add_precedence(ActivityPrecedence.OrJoin([B7a, B7b], B8))
    T2.add_precedence(ActivityPrecedence.Serial(C0, C1))
    T2.add_precedence(ActivityPrecedence.Serial(C3, C4, C5))

    return model



if __name__ == '__main__':
    model = gallery_multitier()
    print('Model built:', type(model).__name__)
