"""
State-Dependent Routing - Closed Network with Two Classes

This example demonstrates state-dependent routing in a closed queueing network
with two job classes. Both classes use round-robin routing at the delay node.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network with two classes and state-dependent routing
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)

job_class1 = ClosedClass(model, 'Class1', 1, delay, 0)
job_class2 = ClosedClass(model, 'Class2', 2, delay, 0)

# Service distributions for Class 1
# Renewal processes
map21 = APH([1, 0], [[-2, 2], [0, -0.5]])
map22 = PH(1, -1)
# Non-renewal process
map31 = MAP([[-20, 0], [0, -1]], [[0, 20], [0.8, 0.2]])
map32 = MAP([[-4, 3], [4, -6]], [[1, 0], [0, 2]])

delay.setService(job_class1, HyperExp.fit_mean_and_scv(1, 25))
queue1.setService(job_class1, map21)
queue2.setService(job_class1, map31)

delay.setService(job_class2, HyperExp.fit_mean_and_scv(1, 25))
queue1.setService(job_class2, map22)
queue2.setService(job_class2, map32)

# Create links
model.addLink(delay, delay)
model.addLink(delay, queue1)
model.addLink(delay, queue2)
model.addLink(queue1, delay)
model.addLink(queue2, delay)

# Set state-dependent routing for both classes
delay.setRouting(job_class1, RoutingStrategy.RROBIN)
queue1.setProbRouting(job_class1, delay, 1.0)
queue2.setProbRouting(job_class1, delay, 1.0)

delay.setRouting(job_class2, RoutingStrategy.RROBIN)
queue1.setProbRouting(job_class2, delay, 1.0)
queue2.setProbRouting(job_class2, delay, 1.0)

print('=== State-Dependent Routing - Closed Network with Two Classes ===\n')

# Solve with CTMC
print('CTMC (exact):')
solver_ctmc = CTMC(model, keep=False)
avg_table_ctmc = solver_ctmc.getAvgTable()
print(avg_table_ctmc)

# Solve with JMT
print('\nJMT Simulation:')
solver_jmt = JMT(model, samples=int(1e5), seed=23000)
avg_table_jmt = solver_jmt.getAvgTable()
print(avg_table_jmt)

# Solve with SSA
print('\nSSA (Stochastic State-space Analysis):')
solver_ssa = SSA(model, verbose=True, samples=int(5e3), seed=23000)
avg_table_ssa = solver_ssa.getAvgTable()
print(avg_table_ssa)

print('\nNote: This example uses various Markovian Arrival Processes (MAP/APH/PH)')
print('      for service distributions, including both renewal and non-renewal processes.')
print('      Both classes use round-robin routing at the delay node.')
