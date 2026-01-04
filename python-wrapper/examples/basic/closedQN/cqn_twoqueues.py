"""
Two Queues with Reducible Routing Matrix

This example demonstrates:
- Closed network with 1 job
- Delay station and 2 queues
- Reducible routing matrix (some nodes can reach only themselves)
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(3, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
    node[2] = Queue(model, 'Queue2', SchedStrategy.FCFS)

    jobclass = ClosedClass(model, 'Class1', 1, node[0], 0)

    node[0].set_service(jobclass, Exp.fit_mean(1.0))  # mean = 1
    node[1].set_service(jobclass, Exp.fit_mean(1.5))  # mean = 1.5
    node[2].set_service(jobclass, Exp.fit_mean(3.0))  # mean = 3.0

    P = model.init_routing_matrix()
    # Reducible routing: Queue1 and Queue2 only route to themselves
    P.set(jobclass, jobclass, node[0], node[0], 0.2)
    P.set(jobclass, jobclass, node[0], node[1], 0.3)
    P.set(jobclass, jobclass, node[0], node[2], 0.5)
    P.set(jobclass, jobclass, node[1], node[1], 1.0)
    P.set(jobclass, jobclass, node[2], node[2], 1.0)
    model.link(P)

    # Solve with MVA
    solver = MVA(model)
    print(f'\nSOLVER: {solver.get_name()}')
    avg_table = solver.avg_table()
