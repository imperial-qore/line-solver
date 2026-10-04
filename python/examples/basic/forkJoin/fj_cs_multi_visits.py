# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join network example 8 - fork-join with multiple visits within the same chain
model = Network('model')
# %%
# Create network
source = Source(model,'Source')
queue1 = Queue(model,'Queue1',SchedStrategy.PS)
queue2 = Queue(model,'Queue2',SchedStrategy.PS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)
sink = Sink(model,'Sink')

jobclass1 = OpenClass(model, 'class1')
jobclass2 = OpenClass(model, 'class2')
# %%
# Service configurations
source.set_arrival(jobclass1, Exp(0.1))
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(1.0))
queue1.set_service(jobclass2, Exp(1.0))
queue2.set_service(jobclass2, Exp(1.0))
# %%
# Set routing matrix with class switching - multiple visits
P = model.init_routing_matrix()

# Class 1 initial routing
P.set(jobclass1, jobclass1, source, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join, 1.0)
P.set(jobclass1, jobclass1, queue2, join, 1.0)

# Class switching: class1 -> class2 after join
P.set(jobclass1, jobclass2, join, fork, 1.0)

# Class 2 routing
P.set(jobclass2, jobclass2, fork, queue1, 1.0)
P.set(jobclass2, jobclass2, fork, queue2, 1.0)
P.set(jobclass2, jobclass2, queue1, join, 1.0)
P.set(jobclass2, jobclass2, queue2, join, 1.0)
P.set(jobclass2, jobclass2, join, sink, 1.0)

model.link(P)
# %%
# Solve with MVA and JMT (matches the MATLAB example)
solver = MVA(model)
print(f'SOLVER: {solver.get_name()}')
avgTable = solver.avg_table()
print(avgTable)

solver = JMT(model, seed=23000)
print(f'SOLVER: {solver.get_name()}')
avgTable = solver.avg_table()
print(avgTable)