# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join network example 6
model = Network('model')
# %%
# Create network
source = Source(model,'Source')
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS)
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS)
fork1 = Fork(model,'Fork1')
join1 = Join(model,'Join1',fork1)
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS)
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS)
fork2 = Fork(model,'Fork2')
join2 = Join(model,'Join2',fork2)
sink = Sink(model,'Sink')

jobclass1 = OpenClass(model, 'class1')
# %%
# Service configurations
source.set_arrival(jobclass1, Exp(0.4))
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(1.0))
queue3.set_service(jobclass1, Exp(1.0))
queue4.set_service(jobclass1, Exp(1.0))
# %%
# Two sequential fork-join stages routing
P = model.init_routing_matrix()
P.set(jobclass1, jobclass1, source, fork1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join1, 1.0)
P.set(jobclass1, jobclass1, queue2, join1, 1.0)
P.set(jobclass1, jobclass1, join1, fork2, 1.0)
P.set(jobclass1, jobclass1, fork2, queue3, 1.0)
P.set(jobclass1, jobclass1, fork2, queue4, 1.0)
P.set(jobclass1, jobclass1, queue3, join2, 1.0)
P.set(jobclass1, jobclass1, queue4, join2, 1.0)
P.set(jobclass1, jobclass1, join2, sink, 1.0)
model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)