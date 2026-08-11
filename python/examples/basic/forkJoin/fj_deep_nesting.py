# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
# Nested Fork-Join Network
model = Network('model')
# %%
# Create nested fork-join structure
delay = Delay(model, 'Delay1')
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS)
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS)
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS)
fork2 = Fork(model,'Fork2')
join2 = Join(model,'Join2', fork2)

jobclass1 = ClosedClass(model, 'class1', 1, delay, 0)
# %%
# Service configurations
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(1.0))
delay.set_service(jobclass1, Exp(0.5))
queue3.set_service(jobclass1, Exp(2.0))
queue4.set_service(jobclass1, Exp(2.0))
# %%
# Nested fork-join routing
P = model.init_routing_matrix()

P.set(jobclass1, jobclass1, delay, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
# Nested structure: queue1 goes to fork2
P.set(jobclass1, jobclass1, queue1, fork2, 1.0)
P.set(jobclass1, jobclass1, fork2, queue3, 1.0)
P.set(jobclass1, jobclass1, fork2, queue4, 1.0)
P.set(jobclass1, jobclass1, queue3, join2, 1.0)
P.set(jobclass1, jobclass1, queue4, join2, 1.0)
P.set(jobclass1, jobclass1, join2, join, 1.0)
# Queue2 goes directly to join
P.set(jobclass1, jobclass1, queue2, join, 1.0)
P.set(jobclass1, jobclass1, join, delay, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000)]
# Add MVA solver
solvers.append(MVA(model))

for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)