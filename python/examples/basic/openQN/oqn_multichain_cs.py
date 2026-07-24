"""
Open Queueing Network with Multiple Chains and Class Switching

This example demonstrates job class switching within multiple open chains:
- Chain 1: Interactive requests (2 classes that switch)
- Chain 2: Batch jobs (2 classes that switch)
- Chain 3: Real-time tasks (1 class)
"""
from line_solver import *
import numpy as np

if __name__ == "__main__":
    model = Network('oqn_multichain')

    # Block 1: nodes
    node = np.empty(6, dtype=object)
    node[0] = Source(model, 'RequestSource')
    node[1] = Queue(model, 'Router', SchedStrategy.FCFS)
    node[2] = Queue(model, 'WebCache', SchedStrategy.PS)
    node[3] = Queue(model, 'AppServer', SchedStrategy.FCFS)
    node[4] = Queue(model, 'DataServer', SchedStrategy.PS)
    node[5] = Sink(model, 'ResponseSink')

    # Block 2: job classes
    # Chain 1: Interactive requests (can switch between InteractiveA and InteractiveB)
    jobclass = np.empty(5, dtype=object)
    jobclass[0] = OpenClass(model, 'InteractiveA', 0)
    jobclass[1] = OpenClass(model, 'InteractiveB', 0)
    # Chain 2: Batch jobs (can switch between BatchA and BatchB)
    jobclass[2] = OpenClass(model, 'BatchA', 0)
    jobclass[3] = OpenClass(model, 'BatchB', 0)
    # Chain 3: Real-time tasks
    jobclass[4] = OpenClass(model, 'RealTime', 0)

    # Block 3: service times
    # Router
    node[1].set_service(jobclass[0], Exp.fit_mean(0.15))
    node[1].set_service(jobclass[1], Exp.fit_mean(0.18))
    node[1].set_service(jobclass[2], Exp.fit_mean(0.2))
    node[1].set_service(jobclass[3], Exp.fit_mean(0.22))
    node[1].set_service(jobclass[4], Exp.fit_mean(0.1))

    # WebCache
    node[2].set_service(jobclass[0], Exp.fit_mean(0.3))
    node[2].set_service(jobclass[1], Exp.fit_mean(0.35))
    node[2].set_service(jobclass[2], Exp.fit_mean(1.0))
    node[2].set_service(jobclass[3], Exp.fit_mean(1.2))
    node[2].set_service(jobclass[4], Exp.fit_mean(0.2))

    # AppServer
    node[3].set_service(jobclass[0], Exp.fit_mean(0.5))
    node[3].set_service(jobclass[1], Exp.fit_mean(0.6))
    node[3].set_service(jobclass[2], Exp.fit_mean(1.5))
    node[3].set_service(jobclass[3], Exp.fit_mean(1.8))
    node[3].set_service(jobclass[4], Exp.fit_mean(0.3))

    # DataServer
    node[4].set_service(jobclass[0], Exp.fit_mean(0.4))
    node[4].set_service(jobclass[1], Exp.fit_mean(0.45))
    node[4].set_service(jobclass[2], Exp.fit_mean(0.8))
    node[4].set_service(jobclass[3], Exp.fit_mean(1.0))
    node[4].set_service(jobclass[4], Exp.fit_mean(0.25))

    # Block 4: arrival rates (reduced to ensure stability at all stations)
    # Interactive requests
    node[0].set_arrival(jobclass[0], Exp(0.3))
    node[0].set_arrival(jobclass[1], Exp(0.1))
    # Batch jobs
    node[0].set_arrival(jobclass[2], Exp(0.15))
    node[0].set_arrival(jobclass[3], Exp(0.12))
    # Real-time tasks
    node[0].set_arrival(jobclass[4], Exp(0.2))

    # Block 5: routing with class switching
    P = model.init_routing_matrix()

    # ===== Chain 1: Interactive requests with switching =====
    # InteractiveA routing
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1)      # Source -> Router
    P.set(jobclass[0], jobclass[0], node[1], node[2], 0.6)    # Router -> WebCache (stay)
    P.set(jobclass[0], jobclass[1], node[1], node[2], 0.4)    # Class switch to InteractiveB
    P.set(jobclass[0], jobclass[0], node[2], node[3], 1)      # WebCache -> AppServer
    P.set(jobclass[0], jobclass[0], node[3], node[4], 1)      # AppServer -> DataServer
    P.set(jobclass[0], jobclass[0], node[4], node[5], 1)      # DataServer -> Sink

    # InteractiveB routing
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1)      # Source -> Router
    P.set(jobclass[1], jobclass[1], node[1], node[2], 0.7)    # Router -> WebCache (stay)
    P.set(jobclass[1], jobclass[0], node[1], node[2], 0.3)    # Class switch to InteractiveA
    P.set(jobclass[1], jobclass[1], node[2], node[3], 1)      # WebCache -> AppServer
    P.set(jobclass[1], jobclass[1], node[3], node[4], 1)      # AppServer -> DataServer
    P.set(jobclass[1], jobclass[1], node[4], node[5], 1)      # DataServer -> Sink

    # ===== Chain 2: Batch jobs with switching =====
    # BatchA routing
    P.set(jobclass[2], jobclass[2], node[0], node[1], 1)      # Source -> Router
    P.set(jobclass[2], jobclass[2], node[1], node[2], 0.5)    # Router -> WebCache (stay)
    P.set(jobclass[2], jobclass[3], node[1], node[2], 0.5)    # Class switch to BatchB
    P.set(jobclass[2], jobclass[2], node[2], node[3], 1)      # WebCache -> AppServer
    P.set(jobclass[2], jobclass[2], node[3], node[4], 1)      # AppServer -> DataServer
    P.set(jobclass[2], jobclass[2], node[4], node[5], 1)      # DataServer -> Sink

    # BatchB routing
    P.set(jobclass[3], jobclass[3], node[0], node[1], 1)      # Source -> Router
    P.set(jobclass[3], jobclass[3], node[1], node[2], 0.6)    # Router -> WebCache (stay)
    P.set(jobclass[3], jobclass[2], node[1], node[2], 0.4)    # Class switch to BatchA
    P.set(jobclass[3], jobclass[3], node[2], node[3], 1)      # WebCache -> AppServer
    P.set(jobclass[3], jobclass[3], node[3], node[4], 1)      # AppServer -> DataServer
    P.set(jobclass[3], jobclass[3], node[4], node[5], 1)      # DataServer -> Sink

    # ===== Chain 3: Real-time tasks =====
    P.set(jobclass[4], jobclass[4], node[0], node[1], 1)      # Source -> Router
    P.set(jobclass[4], jobclass[4], node[1], node[2], 1)      # Router -> WebCache
    P.set(jobclass[4], jobclass[4], node[2], node[3], 1)      # WebCache -> AppServer
    P.set(jobclass[4], jobclass[4], node[3], node[4], 1)      # AppServer -> DataServer
    P.set(jobclass[4], jobclass[4], node[4], node[5], 1)      # DataServer -> Sink

    model.link(P)

    # Block 6: solve with multiple solvers
    print('Open QN with 3 chains, 5 classes, and class switching:')
    print('- Chain 1: InteractiveA and InteractiveB (switch at Router)')
    print('- Chain 2: BatchA and BatchB (switch at Router)')
    print('- Chain 3: RealTime')
    print()

    solver = np.array([], dtype=object)
    solver = np.append(solver, MVA(model))
    solver = np.append(solver, FLD(model))
    solver = np.append(solver, SSA(model, seed=23000))

    for s in range(len(solver)):
        print(f'SOLVER: {solver[s].get_name()}')
        avg_table = solver[s].avg_table()
        print(avg_table)
        print()
