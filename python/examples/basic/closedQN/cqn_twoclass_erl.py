"""
Two-Class Closed Network with Erlang Distributions and Class Switching

This example demonstrates:
- 2 closed classes with different populations (15 and 5)
- Class switching between classes
- Mix of Exp and Erlang service distributions
- Different routing strategies (RAND, RROBIN, WRROBIN)
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('model')

    node = np.empty(4, dtype=object)
    node[0] = ClassSwitch(model, 'CS', [[0, 1], [1, 0]])
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    node[2] = Queue(model, 'Queue2', SchedStrategy.PS)
    node[3] = Queue(model, 'Delay', SchedStrategy.INF)

    jobclass = np.empty(2, dtype=object)
    jobclass[0] = ClosedClass(model, 'Class1', 15, node[3], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 5, node[3], 0)

    node[1].set_service(jobclass[0], Exp.fit_mean(1.5))  # mean = 1.5
    node[1].set_service(jobclass[1], Erlang.fit_mean_and_order(1.5, 2))  # mean = 1.5

    node[2].set_service(jobclass[0], Erlang.fit_mean_and_order(1.5, 2))  # mean = 1.5
    node[2].set_service(jobclass[1], Exp.fit_mean(1.5))  # mean = 1.5

    node[3].set_service(jobclass[0], Exp.fit_mean(1.0))  # mean = 1
    node[3].set_service(jobclass[1], Exp.fit_mean(1.0))  # mean = 1

    # Set routing
    model.add_link(node[1], node[0])
    model.add_link(node[2], node[0])
    model.add_link(node[0], node[3])
    model.add_link(node[3], node[1])
    model.add_link(node[3], node[2])

    # Set routing strategies
    node[0].set_routing(jobclass[0], RoutingStrategy.RAND)
    node[1].set_routing(jobclass[0], RoutingStrategy.RAND)
    node[2].set_routing(jobclass[0], RoutingStrategy.RAND)
    node[3].set_routing(jobclass[0], RoutingStrategy.RROBIN)

    node[0].set_routing(jobclass[1], RoutingStrategy.RAND)
    node[1].set_routing(jobclass[1], RoutingStrategy.RAND)
    node[2].set_routing(jobclass[1], RoutingStrategy.RAND)
    node[3].set_routing(jobclass[1], RoutingStrategy.WRROBIN, node[1], 1)
    node[3].set_routing(jobclass[1], RoutingStrategy.WRROBIN, node[2], 2)

    # Solve
    solver = np.array([], dtype=object)
    solver = np.append(solver, JMT(model, seed=23000, verbose=True))
    solver = np.append(solver, DES(model, seed=23000, verbose=True))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
