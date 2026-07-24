# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Class switching example 1
# A basic M/M/1 with explicit definition of a ClassSwitch node
# Recommended ClassSwitch declaration style
model = Network('mm1cs')
# %%
# Block 1: nodes
node = np.empty(4, dtype=object)
node[0] = Source(model, 'Source 1')
node[1] = Queue(model, 'Queue 1', SchedStrategy.FCFS)
node[2] = Sink(model, 'Sink 1')
node[3] = ClassSwitch(model, 'ClassSwitch 1')

# Block 2: classes
jobclass = np.empty(2, dtype=object)
jobclass[0] = OpenClass(model, 'Class1', 0)
jobclass[1] = OpenClass(model, 'Class2', 0)
# %%
# Service configurations
node[0].set_arrival(jobclass[0], Exp.fit_mean(10.000000))  # (Source 1,Class1)
node[0].set_arrival(jobclass[1], Exp.fit_mean(2.000000))   # (Source 1,Class2)
node[1].set_service(jobclass[0], Exp.fit_mean(1.000000))   # (Queue 1,Class1)
node[1].set_service(jobclass[1], Exp.fit_mean(1.000000))   # (Queue 1,Class2)
# %%
# Block 3: topology
# The class switching matrix can now be declared after the classes, so the
# ClassSwitch node can be declared outside Block 1.
csmatrix = node[3].init_class_switch_matrix()  # element (i,j) = probability that class i switches to j

# Get class indices for the matrix (0-based for numpy array indexing)
class1_idx = jobclass[0].get_index0()
class2_idx = jobclass[1].get_index0()

# Set class switching probabilities using array indexing
csmatrix[class1_idx, class1_idx] = 0.3
csmatrix[class1_idx, class2_idx] = 0.7
csmatrix[class2_idx, class1_idx] = 1.0
node[3].set_class_switching_matrix(csmatrix)

P = model.init_routing_matrix()  # initialize routing matrix
P.set(jobclass[0], jobclass[0], node[0], node[3], 1.0)  # (Source 1,Class1) -> (ClassSwitch 1,Class1)
P.set(jobclass[0], jobclass[0], node[1], node[2], 1.0)  # (Queue 1,Class1) -> (Sink 1,Class1)
P.set(jobclass[0], jobclass[0], node[3], node[1], 1.0)  # (ClassSwitch 1,Class1) -> (Queue 1,Class1)
P.set(jobclass[1], jobclass[1], node[0], node[3], 1.0)  # (Source 1,Class2) -> (ClassSwitch 1,Class2)
P.set(jobclass[1], jobclass[1], node[1], node[2], 1.0)  # (Queue 1,Class2) -> (Sink 1,Class2)
P.set(jobclass[1], jobclass[1], node[3], node[1], 1.0)  # (ClassSwitch 1,Class2) -> (Queue 1,Class2)
model.link(P)
model.print_routing_matrix()
# %%
# Solve
solver = MVA(model)
AvgTable = solver.avg_chain_table()
print(AvgTable)