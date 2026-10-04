# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join network example 9
model = Network('model')
# %%
# Create network
delay = Delay(model, 'Delay1')
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS)
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0)
jobclass2 = ClosedClass(model, 'class2', 10, delay, 0)
# %%
# Service configurations
# Class 1
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(2.0))
delay.set_service(jobclass1, Exp(0.5))

# Class 2
queue1.set_service(jobclass2, Exp(1.0))
queue2.set_service(jobclass2, Exp(2.0))
delay.set_service(jobclass2, Exp(0.2))
# %%
# Set routing matrix
P = model.init_routing_matrix()

# Class 1 routing
P.set(jobclass1, jobclass1, delay, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join, 1.0)
P.set(jobclass1, jobclass1, queue2, join, 1.0)
P.set(jobclass1, jobclass1, join, delay, 1.0)

# Class 2 routing
P.set(jobclass2, jobclass2, delay, fork, 1.0)
P.set(jobclass2, jobclass2, fork, queue1, 1.0)
P.set(jobclass2, jobclass2, fork, queue2, 1.0)
P.set(jobclass2, jobclass2, queue1, join, 1.0)
P.set(jobclass2, jobclass2, queue2, join, 1.0)
P.set(jobclass2, jobclass2, join, delay, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)