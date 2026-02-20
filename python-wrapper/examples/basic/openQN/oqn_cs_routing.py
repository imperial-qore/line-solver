"""
Open Network with Class Switching

This example demonstrates:
- Class switching in open networks
- 3 classes (A, B, C) with class transitions
- Classes A and B arrive from source, switch to C
- Complex routing with class changes
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('myModel')

    # Nodes
    node = np.empty(5, dtype=object)
    node[0] = Source(model, 'Source 1')
    node[1] = Queue(model, 'Queue 1', SchedStrategy.PS)
    node[2] = ClassSwitch(model, 'ClassSwitch 1')
    node[3] = Sink(model, 'Sink 1')
    node[4] = Queue(model, 'Queue 2', SchedStrategy.PS)

    # Classes
    jobclass = np.empty(3, dtype=object)
    jobclass[0] = OpenClass(model, 'Class A', 0)
    jobclass[1] = OpenClass(model, 'Class B', 0)
    jobclass[2] = OpenClass(model, 'Class C', 0)

    # Arrivals
    node[0].set_arrival(jobclass[0], Exp.fit_mean(0.5))
    node[0].set_arrival(jobclass[1], Exp.fit_mean(1.0))
    node[0].set_arrival(jobclass[2], Disabled.get_instance())

    # Service times
    node[1].set_service(jobclass[0], Exp.fit_mean(0.2))
    node[1].set_service(jobclass[1], Exp.fit_mean(0.3))
    node[1].set_service(jobclass[2], Exp.fit_mean(0.333333))

    node[4].set_service(jobclass[0], Exp.fit_mean(1.0))
    node[4].set_service(jobclass[1], Exp.fit_mean(1.0))
    node[4].set_service(jobclass[2], Exp.fit_mean(0.15))

    # Class switching matrix (identity)
    C = node[2].init_class_switch_matrix()
    for i in range(len(C)):
        for j in range(len(C[i])):
            C[i][j] = 1.0 if i == j else 0.0
    node[2].set_class_switching_matrix(C)

    # Routing
    P = model.init_routing_matrix()

    # Class A routing
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[2], 1.0)
    P.set(jobclass[0], jobclass[2], node[2], node[4], 1.0)
    P.set(jobclass[0], jobclass[0], node[4], node[3], 1.0)

    # Class B routing
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[1], node[2], 1.0)
    P.set(jobclass[1], jobclass[2], node[2], node[4], 1.0)
    P.set(jobclass[1], jobclass[1], node[4], node[3], 1.0)

    # Class C routing
    P.set(jobclass[2], jobclass[2], node[0], node[1], 1.0)
    P.set(jobclass[2], jobclass[2], node[1], node[2], 1.0)
    P.set(jobclass[2], jobclass[2], node[2], node[4], 1.0)
    P.set(jobclass[2], jobclass[2], node[4], node[3], 1.0)

    model.link(P)

    # Run multiple solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, keep=True, verbose=True, cutoff=[[1, 1, 0], [3, 3, 0], [0, 0, 3]], seed=23000))
    solver = np.append(solver, FLD(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, MVA(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, MAM(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, NC(model, keep=True, verbose=True, seed=23000))
    solver = np.append(solver, JMT(model, keep=True, verbose=True, seed=23000, samples=100000))
    solver = np.append(solver, SSA(model, keep=True, verbose=True, seed=23000, samples=100000))
    solver = np.append(solver, DES(model, keep=True, verbose=True, seed=23000, samples=100000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name().replace("Solver", "")}')
        avg_table[s] = solver[s].avg_table()
        print(avg_table[s])
