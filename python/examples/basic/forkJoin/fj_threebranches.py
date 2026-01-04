"""
Fork-Join with Three Branches and Serial Stages

This example demonstrates:
- Fork-join with 3 parallel branches
- Branch 1: Queue1
- Branch 2: Queue2 → Queue3 (serial stages)
- Two classes with different service times
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    delay = Delay(model, 'Delay1')
    queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
    queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
    queue3 = Queue(model, 'Queue3', SchedStrategy.PS)
    fork = Fork(model, 'Fork')
    join = Join(model, 'Join', fork)

    jobclass1 = ClosedClass(model, 'class1', 10, delay, 0)
    jobclass2 = ClosedClass(model, 'class2', 10, delay, 0)

    # Class 1 service times
    queue1.set_service(jobclass1, Exp(1.5))
    queue2.set_service(jobclass1, Exp(1.1))
    queue3.set_service(jobclass1, Exp(2.5))
    delay.set_service(jobclass1, Exp(0.5))

    # Class 2 service times
    queue1.set_service(jobclass2, Exp(2.8))
    queue2.set_service(jobclass2, Exp(3.0))
    queue3.set_service(jobclass2, Exp(1.0))
    delay.set_service(jobclass2, Exp(0.8))

    P = model.init_routing_matrix()

    # Class 1 routing
    P.set(jobclass1, jobclass1, delay, fork, 1.0)
    P.set(jobclass1, jobclass1, fork, queue1, 1.0)
    P.set(jobclass1, jobclass1, fork, queue2, 1.0)
    P.set(jobclass1, jobclass1, queue2, queue3, 1.0)
    P.set(jobclass1, jobclass1, queue3, join, 1.0)
    P.set(jobclass1, jobclass1, queue1, join, 1.0)
    P.set(jobclass1, jobclass1, join, delay, 1.0)

    # Class 2 routing
    P.set(jobclass2, jobclass2, delay, fork, 1.0)
    P.set(jobclass2, jobclass2, fork, queue1, 1.0)
    P.set(jobclass2, jobclass2, fork, queue2, 1.0)
    P.set(jobclass2, jobclass2, queue2, queue3, 1.0)
    P.set(jobclass2, jobclass2, queue3, join, 1.0)
    P.set(jobclass2, jobclass2, queue1, join, 1.0)
    P.set(jobclass2, jobclass2, join, delay, 1.0)

    model.link(P)

    solver = MVA(model)
    print(f'SOLVER: {solver.get_name()}')
    avg_table = solver.avg_table()
