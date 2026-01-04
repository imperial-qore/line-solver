"""
Open Network with HOL (Head-of-Line) Priority

This example demonstrates:
- Open network with HOL priority scheduling
- 3 classes with different priorities (Class2 has priority 1)
- Multiple scheduling strategies: FCFS, SIRO, PS, HOL
- HOL gives strict priority to higher-priority classes
"""

from line_solver import *

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('MyNetwork')

    # Nodes
    node = np.empty(6, dtype=object)
    node[0] = Source(model, 'Source')
    node[1] = Queue(model, 'WebServer', SchedStrategy.FCFS)
    node[2] = Queue(model, 'Storage1', SchedStrategy.SIRO)
    node[3] = Queue(model, 'Storage2', SchedStrategy.PS)
    node[4] = Queue(model, 'Storage3', SchedStrategy.HOL)
    node[5] = Sink(model, 'Sink')

    # Classes (Class2 has priority 1, others have priority 0)
    jobclass = np.empty(3, dtype=object)
    jobclass[0] = OpenClass(model, 'Class1', 0)
    jobclass[1] = OpenClass(model, 'Class2', 1)
    jobclass[2] = OpenClass(model, 'Class3', 0)

    # Arrivals
    node[0].set_arrival(jobclass[0], Exp.fit_mean(10.0))
    node[0].set_arrival(jobclass[1], Exp.fit_mean(10.0))
    node[0].set_arrival(jobclass[2], Exp.fit_mean(10.0))

    # Service times
    node[1].set_service(jobclass[0], Exp.fit_mean(0.3))
    node[1].set_service(jobclass[1], Exp.fit_mean(0.5))
    node[1].set_service(jobclass[2], Exp.fit_mean(0.6))

    node[2].set_service(jobclass[0], Exp.fit_mean(1.1))
    node[2].set_service(jobclass[1], Exp.fit_mean(1.3))
    node[2].set_service(jobclass[2], Exp.fit_mean(1.5))

    node[3].set_service(jobclass[0], Exp.fit_mean(2.0))
    node[3].set_service(jobclass[1], Exp.fit_mean(2.1))
    node[3].set_service(jobclass[2], Exp.fit_mean(1.9))

    node[4].set_service(jobclass[0], Exp.fit_mean(2.5))
    node[4].set_service(jobclass[1], Exp.fit_mean(1.9))
    node[4].set_service(jobclass[2], Exp.fit_mean(4.3))

    # Routing
    P = model.init_routing_matrix()

    # Class 1 routing
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[2], 0.25)
    P.set(jobclass[0], jobclass[0], node[1], node[3], 0.25)
    P.set(jobclass[0], jobclass[0], node[1], node[4], 0.25)
    P.set(jobclass[0], jobclass[0], node[1], node[5], 0.25)
    P.set(jobclass[0], jobclass[0], node[2], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[3], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[4], node[1], 1.0)

    # Class 2 routing
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[1], node[2], 0.25)
    P.set(jobclass[1], jobclass[1], node[1], node[3], 0.25)
    P.set(jobclass[1], jobclass[1], node[1], node[4], 0.25)
    P.set(jobclass[1], jobclass[1], node[1], node[5], 0.25)
    P.set(jobclass[1], jobclass[1], node[2], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[3], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[4], node[1], 1.0)

    # Class 3 routing
    P.set(jobclass[2], jobclass[2], node[0], node[1], 1.0)
    P.set(jobclass[2], jobclass[2], node[1], node[2], 0.25)
    P.set(jobclass[2], jobclass[2], node[1], node[3], 0.25)
    P.set(jobclass[2], jobclass[2], node[1], node[4], 0.25)
    P.set(jobclass[2], jobclass[2], node[1], node[5], 0.25)
    P.set(jobclass[2], jobclass[2], node[2], node[1], 1.0)
    P.set(jobclass[2], jobclass[2], node[3], node[1], 1.0)
    P.set(jobclass[2], jobclass[2], node[4], node[1], 1.0)

    model.link(P)

    # Run solvers
    solver = np.array([], dtype=object)
    solver = np.append(solver, CTMC(model, seed=23000, cutoff=1, samples=10000))
    solver = np.append(solver, MVA(model, seed=23000))
    solver = np.append(solver, JMT(model, seed=23000, samples=10000))
    solver = np.append(solver, SSA(model, seed=23000, samples=10000))

    avg_table = np.empty(len(solver), dtype=object)
    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table[s] = solver[s].avg_table()
