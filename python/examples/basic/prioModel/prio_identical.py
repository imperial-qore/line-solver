"""
GPS Priority with Identical Priority Levels

This example demonstrates:
- GPS (Generalized Processor Sharing) priority scheduling
- 4 classes with different priority levels (0, 1, 1, 2)
- Classes 2 and 3 have identical priority (1)
- Service weights for GPS scheduling
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(2, dtype=object)
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.GPSPRIO)

    jobclass = np.empty(4, dtype=object)
    jobclass[0] = ClosedClass(model, 'Class1', 6, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 4, node[0], 1)
    jobclass[2] = ClosedClass(model, 'Class3', 4, node[0], 1)
    jobclass[3] = ClosedClass(model, 'Class4', 1, node[0], 2)

    node[0].set_service(jobclass[0], Erlang(3.0, 2))
    node[0].set_service(jobclass[1], Exp(1.0))
    node[0].set_service(jobclass[2], Exp(1.0))
    node[0].set_service(jobclass[3], Exp(2.0))

    # Set service with weights for GPS
    w1 = 12
    node[1].set_service(jobclass[0], Exp(30), w1)
    w2 = 3
    node[1].set_service(jobclass[1], Exp(2), w2)
    w3 = 5
    node[1].set_service(jobclass[2], Exp(12), w3)
    w4 = 1
    node[1].set_service(jobclass[3], Exp(1), w4)

    # Serial routing for all classes
    P = model.init_routing_matrix()
    routing = Network.serial_routing(node)
    for r in range(4):
        P.set(jobclass[r], jobclass[r], node[0], node[1], routing[0][1])
        P.set(jobclass[r], jobclass[r], node[1], node[0], routing[1][0])
    model.link(P)

    # Run solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model))
    solver = np.append(solver, JMT(model, seed=23000, samples=30000, keep=True))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
