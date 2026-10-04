# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join with Class Switching
model = Network('model')
# %%
# Create network
delay = Delay(model, 'Delay1')
delay2 = Delay(model, 'Delay2')
queue1 = Queue(model,'Queue1',SchedStrategy.PS)
queue2 = Queue(model,'Queue2',SchedStrategy.PS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0)
jobclass2 = ClosedClass(model, 'class2', 10, delay, 0)
# %%
# Service configurations
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(1.0))
delay.set_service(jobclass1, Exp(0.5))
delay.set_service(jobclass2, Exp(0.5))
delay2.set_service(jobclass2, Exp(2.0))
# %%
# Set routing matrix with class switching
P = model.init_routing_matrix()

# Class switching routing: jobclass1 -> jobclass2 -> jobclass1
P.set(jobclass1, jobclass2, delay, delay2, 1.0)    # Class switch: class1 -> class2 at delay
P.set(jobclass2, jobclass1, delay2, fork, 1.0)     # Class switch: class2 -> class1 at delay2

# Fork-join routing for class1
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join, 1.0)
P.set(jobclass1, jobclass1, queue2, join, 1.0)
P.set(jobclass1, jobclass1, join, delay, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)