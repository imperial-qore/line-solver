"""
Minimal Closed Fork-Join Network (native CTMC/SSA)

This example demonstrates:
- Closed network with 1 job
- Delay -> Fork -> Parallel FCFS Queues -> Join -> Delay
- Exact native fork-join analysis via the tag-augmented state space
- Closed form for N=1: X = 1/(Z + E[max(S1,S2)]) = 0.612244898
"""

from line_solver import *


def fj_tiny_closed():

    model = Network('model')

    delay = Delay(model, 'Delay1')
    queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
    fork = Fork(model, 'Fork1')
    join = Join(model, 'Join1', fork)

    jobclass1 = ClosedClass(model, 'Class1', 1, delay)

    delay.set_service(jobclass1, Exp(1.0))
    queue1.set_service(jobclass1, Exp(2.0))
    queue2.set_service(jobclass1, Exp(3.0))

    P = model.init_routing_matrix()
    P.set(jobclass1, jobclass1, delay, fork, 1.0)
    P.set(jobclass1, jobclass1, fork, queue1, 1.0)
    P.set(jobclass1, jobclass1, fork, queue2, 1.0)
    P.set(jobclass1, jobclass1, queue1, join, 1.0)
    P.set(jobclass1, jobclass1, queue2, join, 1.0)
    P.set(jobclass1, jobclass1, join, delay, 1.0)

    model.link(P)
    return model


if __name__ == '__main__':
    GlobalConstants.set_verbose(VerboseLevel.STD)
    model = fj_tiny_closed()

    solver = [CTMC(model),
              SSA(model, seed=23000, samples=50000),
              JMT(model, seed=23000)]

    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name().replace("Solver", "")}')
        print(solver[s].avg_table())
