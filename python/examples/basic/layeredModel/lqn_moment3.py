"""
Layered Queueing Network - LN solver with moment3 method

This example demonstrates:
- LQN with 2 processors and 2 tasks
- Task 1 (reference task) makes multiple synchronous calls
- Task 2 has multiple entries (E2 and E3)
- Serial activity precedence within tasks
- Comparing LN solver default vs moment3 methods
"""

from line_solver import *


def lqn_moment3():

    model = LayeredNetwork('myLayeredModel')

    # Definition of processors, tasks and entries
    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 100, SchedStrategy.REF).on(P1)
    E1 = Entry(model, 'E1').on(T1)

    P2 = Processor(model, 'P2', 1, SchedStrategy.PS)
    T2 = Task(model, 'T2', 1, SchedStrategy.INF).on(P2)
    E2 = Entry(model, 'E2').on(T2)
    E3 = Entry(model, 'E3').on(T2)

    # Definition of activities
    T1.set_think_time(Erlang.fit_mean_and_order(10, 1))

    A1 = Activity(model, 'A1', Exp(1)).on(T1).bound_to(E1).synch_call(E2).synch_call(E3, 1)

    A20 = Activity(model, 'A20', Exp(1)).on(T2).bound_to(E2)
    A21 = Activity(model, 'A21', Exp(1)).on(T2)
    A22 = Activity(model, 'A22', Exp(1)).on(T2).replies_to(E2)
    T2.add_precedence(ActivityPrecedence.serial([A20, A21, A22]))

    A3 = Activity(model, 'A3', Exp(1)).on(T2).bound_to(E3).replies_to(E3)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = lqn_moment3()

    # LN solver with default method (using NC layer solver)
    solver_ln = LN(model, lambda m: NC(m, verbose=False), verbose=False)
    avg_table_ln = solver_ln.get_avg_table()
    print('\nLN Results:')
    print(avg_table_ln)

    # LN solver with moment3 method (using NC layer solver)
    try:
        ln_opts = SolverLNOptions()
        ln_opts.method = 'moment3'
        ln_opts.verbose = False
        solver_ln2 = LN(model, lambda m: NC(m, verbose=False), ln_opts)
        avg_table_ln2 = solver_ln2.get_avg_table()
        print('\nLN(moment3) Results:')
        print(avg_table_ln2)
    except Exception as e:
        import sys
        print(f'\nLN(moment3) failed: {e}', file=sys.stderr)
        print('\nLN(moment3) Results:')
        print(avg_table_ln)
