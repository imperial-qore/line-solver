"""
Layered Queueing Network with an entry-level open arrival.

This example demonstrates:
- LQN with 1 processor and 1 task
- A single entry receiving an external Poisson arrival stream
- A bound activity with Exp service that processes each request

Expected: arrivals at rate 0.2 against a mean service of 1.6 take 0.32 of the
host. Nothing else reaches T1, so T1 has no task layer and SolverLN represents
the stream by the thread pool it drives -- a closed chain of mult(T1) jobs whose
surrogate delay is closed on the known rate, the construction a forwarding
target gets (_open_arrival_rate_of). Reported: entry throughput 0.2, host
utilization 0.32, entry response time 1.6, which lqns gives exactly and lqsim
(0.192-0.200) and LDES (0.19986 / 0.31957 / 1.599) confirm. MATLAB, the JAR and
the C++ port agree. Placing an open class on the host layer instead would load
it twice, since that chain has no other delay to cycle against: that was the
earlier reading of 0.425 / 0.68 / 2.3529.

Note: the Source/Sink/OpenClass plumbing for entry-level open arrivals is still
exercised whenever the entry ALSO has a sync/async caller or is a forwarding
target; the LN loop re-applies the static arrival distribution each iteration
via arvproc_classes_updmap.
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
