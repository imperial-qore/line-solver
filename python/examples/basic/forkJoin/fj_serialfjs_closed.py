# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Series of Fork-Join Networks
model = Network('model')
# %%
# Create series of fork-join stages
delay = Delay(model, 'Delay1')
queue1 = Queue(model,'Queue1',SchedStrategy.PS)
queue2 = Queue(model,'Queue2',SchedStrategy.PS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)
queue3 = Queue(model,'Queue3',SchedStrategy.PS)
queue4 = Queue(model,'Queue4',SchedStrategy.PS)
fork2 = Fork(model,'Fork2')
join2 = Join(model,'Join2', fork2)

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0)
# %%
# Service configurations
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(1.0))
delay.set_service(jobclass1, Exp(0.5))
queue3.set_service(jobclass1, Exp(1.0))
queue4.set_service(jobclass1, Exp(1.0))
# %%
# Series fork-join routing
P = model.init_routing_matrix()

# First fork-join stage
P.set(jobclass1, jobclass1, delay, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join, 1.0)
P.set(jobclass1, jobclass1, queue2, join, 1.0)
# Between stages: join -> fork2
P.set(jobclass1, jobclass1, join, fork2, 1.0)
# Second fork-join stage
P.set(jobclass1, jobclass1, fork2, queue3, 1.0)
P.set(jobclass1, jobclass1, fork2, queue4, 1.0)
P.set(jobclass1, jobclass1, queue3, join2, 1.0)
P.set(jobclass1, jobclass1, queue4, join2, 1.0)
P.set(jobclass1, jobclass1, join2, delay, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)