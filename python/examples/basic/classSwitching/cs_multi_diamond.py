# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Class switching example 2
# Example of class switching controlled by a reducible Markov chain
model = Network('mm1cs')
# %%
# Block 1: nodes
node = np.empty(5, dtype=object)
node[0] = Source(model, 'Source 1')
node[1] = Queue(model, 'Queue 0', SchedStrategy.FCFS)
node[2] = Queue(model, 'Queue 1', SchedStrategy.FCFS)
node[3] = Queue(model, 'Queue 2', SchedStrategy.FCFS)
node[4] = Sink(model, 'Sink 1')

# Block 2: classes
jobclass = np.empty(3, dtype=object)
jobclass[0] = OpenClass(model, 'Class1', 0)
jobclass[1] = OpenClass(model, 'Class2', 0)
jobclass[2] = OpenClass(model, 'Class3', 0)
# %%
# Service configurations
# Note: Exp(rate) creates exponential with given rate (matching MATLAB's Exp(rate))
node[0].set_arrival(jobclass[0], Exp(1.0))   # (Source 1,Class1) rate=1
node[1].set_service(jobclass[0], Exp(10.0))  # (Queue 0,Class1) rate=10
node[2].set_service(jobclass[1], Exp(20.0))  # (Queue 1,Class2) rate=20
node[3].set_service(jobclass[2], Exp(30.0))  # (Queue 2,Class3) rate=30
# %%
# Routing matrix with class switching
P = model.init_routing_matrix()  # initialize routing matrix
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[1], 0.2)
P.set(jobclass[0], jobclass[1], node[1], node[2], 0.3)
P.set(jobclass[0], jobclass[2], node[1], node[3], 0.5)
P.set(jobclass[1], jobclass[1], node[2], node[4], 1.0)
P.set(jobclass[2], jobclass[2], node[3], node[4], 1.0)
model.link(P)
model.print_routing_matrix()
# %%
# Solve
solver = MVA(model)
AvgTable = solver.avg_chain_table()
print(AvgTable)