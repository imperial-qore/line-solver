# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Fork-Join network example 7
model = Network('model')
# %%
# Create network
delay = Delay(model,'Delay')
fork1 = Fork(model,'Fork1')
join1 = Join(model,'Join1',fork1)
queue1 = Queue(model,'Queue1',SchedStrategy.PS)
queue2 = Queue(model,'Queue2',SchedStrategy.PS)

jobclass1 = ClosedClass(model, 'class1', 1, delay, 0)
jobclass2 = ClosedClass(model, 'class2', 1, delay, 0)
# %%
# Service configurations
delay.set_service(jobclass1, Exp(0.25))
queue1.set_service(jobclass1, Exp(2.0))
queue2.set_service(jobclass1, Exp(2.0))

delay.set_service(jobclass2, Exp(0.25))
queue1.set_service(jobclass2, Exp(2.0))
queue2.set_service(jobclass2, Exp(2.0))
# %%
# Set routing matrix with class switching
P = model.init_routing_matrix()

# Class 1 routing
P.set(jobclass1, jobclass1, delay, fork1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join1, 1.0)
P.set(jobclass1, jobclass1, queue2, join1, 1.0)
P.set(jobclass1, jobclass1, join1, delay, 1.0)

# Class 2 routing (note the class switching from jobclass2 to jobclass1)
P.set(jobclass2, jobclass2, delay, fork1, 1.0)
P.set(jobclass2, jobclass1, fork1, queue1, 1.0)
P.set(jobclass2, jobclass1, fork1, queue2, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [MVA(model), JMT(model, seed=23000)] # JMT has a bug on this one
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)