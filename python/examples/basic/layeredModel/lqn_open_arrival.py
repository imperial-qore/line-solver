"""
Layered Queueing Network with an entry-level open arrival.

This example demonstrates:
- LQN with 1 processor and 1 task
- A single entry receiving an external Poisson arrival stream
- A bound activity with Exp service that processes each request

Expected: Source delivers arrivals at rate 0.2, the server processes them
with mean 1.6, so E1_Open throughput ~0.2 and processor utilization ~0.32.

Note: pure-open LQN layers exercise the Source/Sink/OpenClass plumbing
added for entry-level open arrivals (mirrors MATLAB buildLayersRecursive.m
lines 213-255 and JAR LN.java lines 718-750). The LN loop
re-applies the static arrival distribution each iteration via
arvproc_classes_updmap.
"""

from line_solver import *


def lqn_open_arrival():
    model = LayeredNetwork('openArrivalLQN')

    P1 = Processor(model, 'P1', 1, SchedStrategy.PS)
    T1 = Task(model, 'T1', 1, SchedStrategy.FCFS).on(P1)
    T1.set_think_time(Immediate())

    E1 = Entry(model, 'E1').on(T1)
    E1.setArrival(Exp(0.2))

    Activity(model, 'A1', Exp.fit_mean(1.6)).on(T1).bound_to(E1).replies_to(E1)

    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = lqn_open_arrival()

    avg_table_ln = LN(model).get_avg_table()
    print('\nLN Results:')
    print(avg_table_ln)
