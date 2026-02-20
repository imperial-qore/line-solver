"""
Layered Queueing Network with Serial Activities

This example demonstrates:
- LQN with 2 processors and 2 tasks
- Serial activity precedence within tasks
- Task 1 makes synchronous calls to Task 2
- Activities are executed in serial order
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    print('This example illustrates the execution on a layered queueing network model.')
    print('Performance indexes now refer to processors, tasks, entries, and activities.')

    model = LayeredNetwork('myLayeredModel')

    P = np.empty(2, dtype=object)
    P[0] = Processor(model, 'P1', 1, SchedStrategy.PS)
    P[1] = Processor(model, 'P2', 1, SchedStrategy.PS)

    T = np.empty(2, dtype=object)
    T[0] = Task(model, 'T1', 10, SchedStrategy.REF).on(P[0]).set_think_time(Exp.fit_mean(100))
    T[1] = Task(model, 'T2', 1, SchedStrategy.FCFS).on(P[1]).set_think_time(Immediate())

    E = np.empty(2, dtype=object)
    E[0] = Entry(model, 'E1').on(T[0])
    E[1] = Entry(model, 'E2').on(T[1])

    A = np.empty(4, dtype=object)
    A[0] = Activity(model, 'AS1', Exp.fit_mean(1.6)).on(T[0]).bound_to(E[0])
    A[1] = Activity(model, 'AS2', Immediate()).on(T[0]).synch_call(E[1], 1)
    A[2] = Activity(model, 'AS3', Exp.fit_mean(5)).on(T[1]).bound_to(E[1])
    A[3] = Activity(model, 'AS4', Exp.fit_mean(1)).on(T[1]).replies_to(E[1])

    # Serial precedence for activities
    T[0].add_precedence(ActivityPrecedence.serial([A[0], A[1]]))
    T[1].add_precedence(ActivityPrecedence.serial([A[2], A[3]]))

    # Solve using LQNS
    solver = SolverLQNS(model, keep=True)
    avg_table = solver.get_avg_table()
    print('\nLQNS Results:')
    print(avg_table)

    # Note: get_raw_avg_tables() not yet implemented in Python
    # avg_table_raw, call_avg_table = solver.get_raw_avg_tables()
    # print('\nRaw Average Table:')
    # print(avg_table_raw)
    # print('\nCall Average Table:')
    # print(call_avg_table)

    # Solve using LN
    avg_table_ln = SolverLN(model).get_avg_table()
    print('\nLN Results:')
    print(avg_table_ln)
