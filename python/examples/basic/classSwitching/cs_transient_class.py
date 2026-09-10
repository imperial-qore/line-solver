# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Class switching example 4
# Example of class switching controlled by a reducible Markov chain
# In this variant the job remains either in class 2 or class 3 forever
model = Network('reducible_cs')
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
node[0].set_service(jobclass[1], Exp.fit_mean(1.000000))  # (Queue 0,Class2)
node[0].set_service(jobclass[2], Exp.fit_mean(1.000000))  # (Queue 0,Class3)
node[1].set_service(jobclass[1], Exp.fit_mean(1.000000))  # (Queue 1,Class2)
node[2].set_service(jobclass[2], Exp.fit_mean(1.000000))  # (Queue 2,Class3)
# %%
# Routing matrix with class switching - transient class example
# Matches MATLAB: Class1 switches to Class2/Class3 when MOVING to Queue1/Queue2
P = model.init_routing_matrix()  # initialize routing matrix

# Class 1 routing: stays at Queue0 or switches class when moving to Queue1/Queue2
P.set(jobclass[0], jobclass[0], node[0], node[0], 0.2)  # Stay in class 1 at Queue 0
P.set(jobclass[0], jobclass[1], node[0], node[1], 0.3)  # Switch to class 2, move to Queue 1
P.set(jobclass[0], jobclass[2], node[0], node[2], 0.5)  # Switch to class 3, move to Queue 2

# Class 2 routing: cycles between Queue 0 and Queue 1
P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)  # Class 2: Queue 0 to Queue 1
P.set(jobclass[1], jobclass[1], node[1], node[0], 1.0)  # Class 2: Queue 1 to Queue 0

# Class 3 routing: cycles between Queue 0 and Queue 2
P.set(jobclass[2], jobclass[2], node[0], node[2], 1.0)  # Class 3: Queue 0 to Queue 2
P.set(jobclass[2], jobclass[2], node[2], node[0], 1.0)  # Class 3: Queue 2 to Queue 0

model.link(P)
# %%
model.print_routing_matrix()
# %%
# Solve
solver = MVA(model)
AvgTable = solver.avg_chain_table()
print(AvgTable)