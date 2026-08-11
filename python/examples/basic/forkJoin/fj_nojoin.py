# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join network example 4
model = Network('model')
# %%
# Create network
source = Source(model,'Source')
queue1 = Queue(model,'Queue1',SchedStrategy.PS)
queue2 = Queue(model,'Queue2',SchedStrategy.PS)
queue3 = Queue(model,'Queue3',SchedStrategy.PS)
fork = Fork(model,'Fork')
sink = Sink(model,'Sink')

jobclass1 = OpenClass(model, 'class1')
# %%
# Service configurations
source.set_arrival(jobclass1, Exp(0.5))
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(2.0))
queue3.set_service(jobclass1, Exp(3.0))
# %%
# Pure fork routing (no join)
P = model.init_routing_matrix()
P.set(jobclass1, jobclass1, source, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, fork, queue3, 1.0)
P.set(jobclass1, jobclass1, queue1, sink, 1.0)
P.set(jobclass1, jobclass1, queue2, sink, 1.0)
P.set(jobclass1, jobclass1, queue3, sink, 1.0)
model.link(P)
# %%
# Solve with multiple methods
solvers = [MVA(model), JMT(model, seed=23000)]

for i, solver in enumerate(solvers):
    print(f'\nSOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)