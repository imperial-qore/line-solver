"""
Mixed Queueing Network with Multiple Chains and Class Switching

Demonstrates class switching within both open and closed job classes.
- Chain 1: Closed class with 2 classes that can switch (interactive users)
- Chain 2: Open classes with 2 classes that can switch (batch requests)
- Chain 3: Open class (external load)
"""

from line_solver import *
import numpy as np

if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = Network('mqn_multichain')

    # Block 1: nodes
    node = np.empty(6, dtype=object)
    node[0] = Delay(model, 'ThinkingTime')
    node[1] = Queue(model, 'WebServer', SchedStrategy.FCFS)
    node[2] = Queue(model, 'AppServer', SchedStrategy.PS)
    node[3] = Queue(model, 'DataServer', SchedStrategy.FCFS)
    node[4] = Source(model, 'Source')
    node[5] = Sink(model, 'Sink')

    # Block 2: job classes
    jobclass = np.empty(5, dtype=object)
    # Chain 1: Closed classes (interactive users that can switch types)
    jobclass[0] = ClosedClass(model, 'InteractiveA', 3, node[0])
    jobclass[1] = ClosedClass(model, 'InteractiveB', 2, node[0])
    # Chain 2: Open classes (batch jobs that can switch types)
    jobclass[2] = OpenClass(model, 'BatchA', 0)
    jobclass[3] = OpenClass(model, 'BatchB', 0)
    # Chain 3: Open class (external load)
    jobclass[4] = OpenClass(model, 'ExternalLoad', 0)

    # Block 3: service times
    # ThinkingTime
    node[0].set_service(jobclass[0], Exp.fit_mean(1.5))
    node[0].set_service(jobclass[1], Exp.fit_mean(2.0))

    # WebServer
    node[1].set_service(jobclass[0], Exp.fit_mean(0.4))
    node[1].set_service(jobclass[1], Exp.fit_mean(0.5))
    node[1].set_service(jobclass[2], Exp.fit_mean(0.3))
    node[1].set_service(jobclass[3], Exp.fit_mean(0.35))
    node[1].set_service(jobclass[4], Exp.fit_mean(0.2))

    # AppServer
    node[2].set_service(jobclass[0], Exp.fit_mean(0.8))
    node[2].set_service(jobclass[1], Exp.fit_mean(1.0))
    node[2].set_service(jobclass[2], Exp.fit_mean(0.6))
    node[2].set_service(jobclass[3], Exp.fit_mean(0.7))
    node[2].set_service(jobclass[4], Exp.fit_mean(0.5))

    # DataServer
    node[3].set_service(jobclass[0], Exp.fit_mean(0.5))
    node[3].set_service(jobclass[1], Exp.fit_mean(0.6))
    node[3].set_service(jobclass[2], Exp.fit_mean(1.2))
    node[3].set_service(jobclass[3], Exp.fit_mean(1.5))
    node[3].set_service(jobclass[4], Exp.fit_mean(0.8))

    # Block 4: arrival rates (reduced for stability)
    node[4].set_arrival(jobclass[2], Exp(0.1))
    node[4].set_arrival(jobclass[3], Exp(0.05))
    node[4].set_arrival(jobclass[4], Exp(0.05))

    # Block 5: routing with class switching
    P = model.init_routing_matrix()

    # ===== Chain 1: Closed classes with switching =====
    # InteractiveA routing:
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[0], jobclass[0], node[1], node[2], 0.7)      # WebServer -> AppServer (stay in InteractiveA)
    P.set(jobclass[0], jobclass[1], node[1], node[2], 0.3)      # Class switch to InteractiveB at WebServer
    P.set(jobclass[0], jobclass[0], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[0], jobclass[0], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # InteractiveB routing:
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[1], jobclass[1], node[1], node[2], 0.8)      # WebServer -> AppServer (stay in InteractiveB)
    P.set(jobclass[1], jobclass[0], node[1], node[2], 0.2)      # Class switch to InteractiveA at WebServer
    P.set(jobclass[1], jobclass[1], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[1], jobclass[1], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # ===== Chain 2: Open classes with switching =====
    # BatchA routing:
    P.set(jobclass[2], jobclass[2], node[4], node[1], 1.0)      # Source -> WebServer
    P.set(jobclass[2], jobclass[2], node[1], node[2], 0.6)      # WebServer -> AppServer (stay in BatchA)
    P.set(jobclass[2], jobclass[3], node[1], node[2], 0.4)      # Class switch to BatchB at WebServer
    P.set(jobclass[2], jobclass[2], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[2], jobclass[2], node[3], node[5], 1.0)      # DataServer -> Sink

    # BatchB routing:
    P.set(jobclass[3], jobclass[3], node[4], node[1], 1.0)      # Source -> WebServer
    P.set(jobclass[3], jobclass[3], node[1], node[2], 0.7)      # WebServer -> AppServer (stay in BatchB)
    P.set(jobclass[3], jobclass[2], node[1], node[2], 0.3)      # Class switch to BatchA at WebServer
    P.set(jobclass[3], jobclass[3], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[3], jobclass[3], node[3], node[5], 1.0)      # DataServer -> Sink

    # ===== Chain 3: External load =====
    # ExternalLoad routing:
    P.set(jobclass[4], jobclass[4], node[4], node[1], 1.0)      # Source -> WebServer
    P.set(jobclass[4], jobclass[4], node[1], node[2], 1.0)      # WebServer -> AppServer
    P.set(jobclass[4], jobclass[4], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[4], jobclass[4], node[3], node[5], 1.0)      # DataServer -> Sink

    model.link(P)

    # Block 6: solve with multiple solvers
    print('Mixed QN with 3 chains, 5 classes, and class switching:')
    print('- Chain 1 (Closed): InteractiveA and InteractiveB (switch at WebServer)')
    print('- Chain 2 (Open): BatchA and BatchB (switch at WebServer)')
    print('- Chain 3 (Open): ExternalLoad')
    print()

    solver = np.array([], dtype=object)
    solver = np.append(solver, MVA(model, verbose=True, seed=23000))
    solver = np.append(solver, SSA(model, verbose=True, cutoff=5, seed=23000))

    for s in range(len(solver)):
        print(f'SOLVER: {solver[s].get_name().replace("Solver", "")}')
        avg_table = solver[s].avg_chain_table()
        print(avg_table)
        print()
