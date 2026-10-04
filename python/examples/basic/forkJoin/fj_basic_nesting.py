# %%
from line_solver import *
import numpy as np
GlobalConstants.set_verbose(VerboseLevel.STD)
# %%
model = Network('model')

# Nodes
delay = Delay(model, 'Delay')
fork1 = Fork(model, 'Fork1')
fork1.set_tasks_per_link(1)
fork11 = Fork(model, 'Fork1_1')
fork11.set_tasks_per_link(2)
join1 = Join(model, 'Join1', fork1)
join11 = Join(model, 'Join1_1', fork11)
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

# Classes
jobclass1 = ClosedClass(model, 'class1', 5, delay, 0)
jobclass2 = ClosedClass(model, 'class2', 2, delay, 0)

# Service times
delay.set_service(jobclass1, Exp(0.25))
queue1.set_service(jobclass1, Exp(1.0))
queue2.set_service(jobclass1, Exp(0.75))

delay.set_service(jobclass2, Exp(0.25))
queue1.set_service(jobclass2, Exp(2.0))
queue2.set_service(jobclass2, Exp(2.0))
# %%
# Routing
P = model.init_routing_matrix()

# Class 1 routing
P.set(jobclass1, jobclass1, delay, fork1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue1, 1.0)
P.set(jobclass1, jobclass1, fork1, queue2, 1.0)
P.set(jobclass1, jobclass1, queue1, join1, 1.0)
P.set(jobclass1, jobclass1, queue2, join1, 1.0)
P.set(jobclass1, jobclass1, join1, delay, 1.0)

# Class 2 routing
P.set(jobclass2, jobclass2, delay, fork11, 1.0)
P.set(jobclass2, jobclass2, fork11, fork1, 1.0)
P.set(jobclass2, jobclass2, fork1, queue1, 1.0)
P.set(jobclass2, jobclass2, fork1, queue2, 1.0)
P.set(jobclass2, jobclass2, queue1, join1, 1.0)
P.set(jobclass2, jobclass2, queue2, join1, 1.0)
P.set(jobclass2, jobclass2, join1, join11, 1.0)
P.set(jobclass2, jobclass2, join11, delay, 1.0)

model.link(P)
# %%
# Solve with different solvers
solver = []

# JMT solver
jmt_options = JMT.default_options()
jmt_options.seed = 23000
solver.append(JMT(model, jmt_options))

# MVA solver
solver.append(MVA(model))

# Get average tables
AvgTable = []
for s in solver:
    print(f'\nSOLVER: {s.get_name()}')
    avg_table = s.avg_table()
    print(avg_table)
    AvgTable.append(avg_table)