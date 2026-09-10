# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
#Advanced Fork-Join Example 12#Complex fork-join with 2 parallel branches
model = Network('model')
# %%
# Create network
delay = Delay(model, 'Delay1')
queue1 = Queue(model,'Queue1',SchedStrategy.FCFS)
queue2 = Queue(model,'Queue2',SchedStrategy.FCFS)
queue3 = Queue(model,'Queue3',SchedStrategy.FCFS)
queue4 = Queue(model,'Queue4',SchedStrategy.FCFS)
queue5 = Queue(model,'Queue5',SchedStrategy.FCFS)
fork = Fork(model,'Fork')
join = Join(model,'Join', fork)

jobclass1 = ClosedClass(model, 'class1', 10, delay, 0)
# %%
# Service configurations
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(2.0))
queue3.set_service(jobclass1, Exp(1.0))
queue4.set_service(jobclass1, Exp(3.0))
queue5.set_service(jobclass1, Exp(0.8))
delay.set_service(jobclass1, Exp(0.5))
# %%
# Set routing matrix with complex serial routing
P = model.init_routing_matrix()

P.set(jobclass1, jobclass1, delay, fork, 1.0)
P.set(jobclass1, jobclass1, fork, queue1, 1.0)
P.set(jobclass1, jobclass1, fork, queue2, 1.0)
# Complex serial routing:
P.set(jobclass1, jobclass1, queue1, queue4, 1.0)  # queue1 -> queue4
P.set(jobclass1, jobclass1, queue4, queue5, 1.0)  # queue4 -> queue5  
P.set(jobclass1, jobclass1, queue5, join, 1.0)    # queue5 -> join
P.set(jobclass1, jobclass1, queue2, queue3, 1.0)  # queue2 -> queue3
P.set(jobclass1, jobclass1, queue3, join, 1.0)    # queue3 -> join
P.set(jobclass1, jobclass1, join, delay, 1.0)

model.link(P)
# %%
# Solve with multiple methods
solvers = [JMT(model, seed=23000), MVA(model)]
for i, solver in enumerate(solvers):
    print(f'SOLVER {i+1}: {solver.get_name()}')
    avgTable = solver.avg_table()
    print(avgTable)