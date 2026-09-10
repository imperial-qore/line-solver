"""
Closed Queueing Network with Multiple Chains and Class Switching

This example demonstrates:
- 3 chains, 5 classes with class switching
- Chain 1: High-priority users (HighPriorityA and HighPriorityB)
- Chain 2: Regular users (RegularA and RegularB)
- Chain 3: Background tasks
- Mix of FCFS and PS scheduling
- Product-form (BCMP) network
"""

from line_solver import *
import numpy as np


def cqn_multichain_cs():
    model = Network('cqn_multichain')

    # Block 1: nodes
    node = np.empty(4, dtype=object)
    node[0] = Delay(model, 'ThinkingTime')
    node[1] = Queue(model, 'WebServer', SchedStrategy.FCFS)
    node[2] = Queue(model, 'AppServer', SchedStrategy.PS)
    node[3] = Queue(model, 'DataServer', SchedStrategy.FCFS)

    # Block 2: job classes
    # Chain 1: High-priority users (can switch between HighPriorityA and HighPriorityB)
    jobclass = np.empty(5, dtype=object)
    jobclass[0] = ClosedClass(model, 'HighPriorityA', 2, node[0])
    jobclass[1] = ClosedClass(model, 'HighPriorityB', 1, node[0])
    # Chain 2: Regular users (can switch between RegularA and RegularB)
    jobclass[2] = ClosedClass(model, 'RegularA', 3, node[0])
    jobclass[3] = ClosedClass(model, 'RegularB', 2, node[0])
    # Chain 3: Background tasks
    jobclass[4] = ClosedClass(model, 'Background', 2, node[0])

    # Block 3: service times (all exponential for product-form)
    # ThinkingTime (infinite server node)
    node[0].set_service(jobclass[0], Exp.fit_mean(1.5))
    node[0].set_service(jobclass[1], Exp.fit_mean(1.6))
    node[0].set_service(jobclass[2], Exp.fit_mean(2.0))
    node[0].set_service(jobclass[3], Exp.fit_mean(2.1))
    node[0].set_service(jobclass[4], Exp.fit_mean(3.0))

    # WebServer (FCFS, exponential)
    node[1].set_service(jobclass[0], Exp.fit_mean(0.4))
    node[1].set_service(jobclass[1], Exp.fit_mean(0.42))
    node[1].set_service(jobclass[2], Exp.fit_mean(0.5))
    node[1].set_service(jobclass[3], Exp.fit_mean(0.52))
    node[1].set_service(jobclass[4], Exp.fit_mean(0.8))

    # AppServer (PS, exponential)
    node[2].set_service(jobclass[0], Exp.fit_mean(0.7))
    node[2].set_service(jobclass[1], Exp.fit_mean(0.75))
    node[2].set_service(jobclass[2], Exp.fit_mean(0.9))
    node[2].set_service(jobclass[3], Exp.fit_mean(0.95))
    node[2].set_service(jobclass[4], Exp.fit_mean(1.2))

    # DataServer (FCFS, exponential)
    node[3].set_service(jobclass[0], Exp.fit_mean(0.5))
    node[3].set_service(jobclass[1], Exp.fit_mean(0.52))
    node[3].set_service(jobclass[2], Exp.fit_mean(0.6))
    node[3].set_service(jobclass[3], Exp.fit_mean(0.62))
    node[3].set_service(jobclass[4], Exp.fit_mean(1.0))

    # Block 4: routing with class switching
    M = 4  # number of nodes
    K = 5  # number of classes

    P = model.init_routing_matrix()

    # ===== Chain 1: High-priority users with switching =====
    # HighPriorityA routing:
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[0], jobclass[0], node[1], node[2], 0.7)      # WebServer -> AppServer (stay in HighPriorityA)
    P.set(jobclass[0], jobclass[1], node[1], node[2], 0.3)      # Class switch to HighPriorityB at WebServer
    P.set(jobclass[0], jobclass[0], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[0], jobclass[0], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # HighPriorityB routing:
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[1], jobclass[1], node[1], node[2], 0.75)     # WebServer -> AppServer (stay in HighPriorityB)
    P.set(jobclass[1], jobclass[0], node[1], node[2], 0.25)     # Class switch to HighPriorityA at WebServer
    P.set(jobclass[1], jobclass[1], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[1], jobclass[1], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # ===== Chain 2: Regular users with switching =====
    # RegularA routing:
    P.set(jobclass[2], jobclass[2], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[2], jobclass[2], node[1], node[2], 0.65)     # WebServer -> AppServer (stay in RegularA)
    P.set(jobclass[2], jobclass[3], node[1], node[2], 0.35)     # Class switch to RegularB at WebServer
    P.set(jobclass[2], jobclass[2], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[2], jobclass[2], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # RegularB routing:
    P.set(jobclass[3], jobclass[3], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[3], jobclass[3], node[1], node[2], 0.7)      # WebServer -> AppServer (stay in RegularB)
    P.set(jobclass[3], jobclass[2], node[1], node[2], 0.3)      # Class switch to RegularA at WebServer
    P.set(jobclass[3], jobclass[3], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[3], jobclass[3], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    # ===== Chain 3: Background tasks =====
    # Background routing:
    P.set(jobclass[4], jobclass[4], node[0], node[1], 1.0)      # ThinkingTime -> WebServer
    P.set(jobclass[4], jobclass[4], node[1], node[2], 1.0)      # WebServer -> AppServer
    P.set(jobclass[4], jobclass[4], node[2], node[3], 1.0)      # AppServer -> DataServer
    P.set(jobclass[4], jobclass[4], node[3], node[0], 1.0)      # DataServer -> ThinkingTime

    model.link(P)

    return model


if __name__ == "__main__":
    GlobalConstants.set_verbose(VerboseLevel.STD)

    model = cqn_multichain_cs()

    # Block 5: solve with multiple solvers
    # Note: CTMC is not included because the state space is too large for this model
    # (10 jobs across 4 stations with 5 classes and class switching).
    # Use MVA for analytical solution or SSA for simulation.
    solver = np.array([], dtype=object)
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, SSA(model, seed=23000, samples=200000))

    for s in range(len(solver)):
        print(f'\nSOLVER: {solver[s].get_name()}')
        avg_table = solver[s].avg_chain_table()
        print(avg_table)
