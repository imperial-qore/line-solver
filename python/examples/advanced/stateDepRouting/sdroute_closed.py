"""
State-Dependent Routing - Closed Network

This example demonstrates state-dependent routing in a closed queueing network.
A delay node uses round-robin routing to distribute jobs between two queues.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network with state-dependent routing
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

job_class = ClosedClass(model, 'Class1', 1, delay, 0)

delay.setService(job_class, HyperExp.fit_mean_and_scv(1, 25))
queue1.setService(job_class, Exp(1))
queue2.setService(job_class, Exp(2))

# Create links
model.addLink(delay, delay)
model.addLink(delay, queue1)
model.addLink(delay, queue2)
model.addLink(queue1, delay)
model.addLink(queue2, delay)

# Set state-dependent routing at delay node
delay.setRouting(job_class, RoutingStrategy.RROBIN)
queue1.setProbRouting(job_class, delay, 1.0)
queue2.setProbRouting(job_class, delay, 1.0)

print('=== State-Dependent Routing - Closed Network ===\n')

# Solve with CTMC
print('CTMC (exact):')
solver_ctmc = CTMC(model, keep=True)
avg_table_ctmc = solver_ctmc.getAvgTable()
print(avg_table_ctmc)

# Solve with JMT
print('\nJMT Simulation:')
solver_jmt = JMT(model, samples=int(1e5), seed=23000)
avg_table_jmt = solver_jmt.getAvgTable()
print(avg_table_jmt)

# Solve with SSA
print('\nSSA (Stochastic State-space Analysis):')
solver_ssa = SSA(model, verbose=True, samples=int(1e4), seed=23000)
avg_table_ssa = solver_ssa.getAvgTable()
print(avg_table_ssa)

print('\nNote: The delay node uses round-robin routing to alternate between Queue1 and Queue2.')
print('      Jobs return from queues back to the delay node with probability 1.0.')
