# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join queueing network
model = Network('model')
# %%
# Create network 
source = Source(model, 'Source')

queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

fork = Fork(model, 'Fork')
fork.set_tasks_per_link(2)
join = Join(model, 'Join', fork)

sink = Sink(model, 'Sink')
# %%
# Create job classes
jobclass1 = OpenClass(model, 'class1')
jobclass2 = OpenClass(model, 'class2')
# %%
#Set arrival and service processes
# Class 1
source.set_arrival(jobclass1, Exp(0.25))
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(0.75))
# Class 2
source.set_arrival(jobclass2, Exp(0.25))
queue1.set_service(jobclass2, Immediate())
queue2.set_service(jobclass2, Exp(2.0))
# %%
# Set routing matrix
P = model.init_routing_matrix() # Class 1 routing
P.set(jobclass1, jobclass1, source, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join, 1.0)
P.set(jobclass1, jobclass1, queue2, join, 1.0)
P.set(jobclass1, jobclass1, join, sink, 1.0) # Class 2 routing
P.set(jobclass2, jobclass2, source, fork, 1.0)
P.set(jobclass2, jobclass2, fork, queue1, 1.0)
P.set(jobclass2, jobclass2, fork, queue2, 1.0)
P.set(jobclass2, jobclass2, queue1, join, 1.0)
P.set(jobclass2, jobclass2, queue2, join, 1.0)
P.set(jobclass2, jobclass2, join, sink, 1.0)
model.link(P)
# %%
# Solve with multiple methods
solvers = [ JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)