"""
Processor Sharing with Priority (PS-PRIO)

This example demonstrates:
- PS-PRIO scheduling (Processor Sharing with Priority)
- 2 classes with different priorities (0 and 1)
- HyperExp and Erlang distributions
- Priority affects service order but maintains PS within priority levels
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('MyNetwork')

    node = np.empty(2, dtype=object)
    node[0] = Delay(model, 'SlowDelay')
    node[1] = Queue(model, 'PSPRIOQueue', SchedStrategy.PSPRIO)

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 2, node[0], 1)

    node[0].set_service(jobclass[0], Erlang(3, 2))
    node[0].set_service(jobclass[1], HyperExp(0.5, 3.0, 10.0))

    node[1].set_service(jobclass[0], HyperExp(0.1, 1.0, 10.0))
    node[1].set_service(jobclass[1], Exp(1))

    # Serial routing for both classes
    P = model.init_routing_matrix()
    routing = Network.serial_routing(node)
    P.set(jobclass[0], jobclass[0], node[0], node[1], routing[0][1])
    P.set(jobclass[0], jobclass[0], node[1], node[0], routing[1][0])
    P.set(jobclass[1], jobclass[1], node[0], node[1], routing[0][1])
    P.set(jobclass[1], jobclass[1], node[1], node[0], routing[1][0])
    model.link(P)

    # Run solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model))
    solver = np.append(solver, JMT(model, seed=23000, verbose=True, samples=5000))
    solver = np.append(solver, SSA(model, seed=23000, verbose=True, samples=5000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
