# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Class switching example 3
# Example of class switching controlled by a reducible Markov chain
model = Network('mm1cs')
# %%
# Block 1: nodes
node = np.empty(3, dtype=object)
node[0] = Delay(model, 'Queue 0')
node[1] = Delay(model, 'Queue 1')
node[2] = Delay(model, 'Queue 2')

# Block 2: classes
jobclass = np.empty(3, dtype=object)
jobclass[0] = ClosedClass(model, 'Class1', 1, node[0])
jobclass[1] = ClosedClass(model, 'Class2', 0, node[0])
jobclass[2] = ClosedClass(model, 'Class3', 0, node[0])
# %%
# Service configurations
node[0].set_service(jobclass[0], Exp.fit_mean(1.000000))  # (Queue 0,Class1)
node[1].set_service(jobclass[1], Exp.fit_mean(2.000000))  # (Queue 1,Class2)
node[2].set_service(jobclass[2], Exp.fit_mean(3.000000))  # (Queue 2,Class3)
# %%
# Routing matrix with class switching
P = model.init_routing_matrix()  # initialize routing matrix
P.set(jobclass[0], jobclass[0], node[0], node[0], 0.2)
P.set(jobclass[0], jobclass[1], node[0], node[1], 0.3)
P.set(jobclass[0], jobclass[2], node[0], node[2], 0.5)
P.set(jobclass[1], jobclass[0], node[1], node[0], 1.0)
P.set(jobclass[2], jobclass[0], node[2], node[0], 1.0)
model.link(P)
model.print_routing_matrix()
# %%
# Solve
solver = MVA(model)
AvgTable = solver.avg_chain_table()
print(AvgTable)