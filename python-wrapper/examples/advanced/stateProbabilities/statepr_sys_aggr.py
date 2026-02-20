"""
State Probabilities - System-Wide Aggregated Probabilities

This example demonstrates computation of system-wide (joint) state probabilities
for the aggregated state space in a multi-class queueing network with class switching.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network model with 4 classes
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
queue2.setNumberOfServers(2)

N = [1, 0, 3, 0]
job_class1 = ClosedClass(model, 'Class1', N[0], delay, 0)
job_class2 = ClosedClass(model, 'Class2', N[1], delay, 0)
job_class3 = ClosedClass(model, 'Class3', N[2], delay, 0)
job_class4 = ClosedClass(model, 'Class4', N[3], delay, 0)

delay.setService(job_class1, Exp(1))
delay.setService(job_class2, Exp(2))
delay.setService(job_class3, Exp(1))
delay.setService(job_class4, Exp(1))

queue1.setService(job_class1, Exp(3))
queue1.setService(job_class2, Exp(4))
queue1.setService(job_class3, Exp(5))
queue1.setService(job_class4, Exp(1))

queue2.setService(job_class1, Exp(1))
queue2.setService(job_class2, Exp(3))
queue2.setService(job_class3, Exp(5))
queue2.setService(job_class4, Exp(2))

K = model.getNumberOfClasses()
P = [[None] * K for _ in range(K)]

P[0][0] = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
P[0][1] = [[0, 0, 0], [0, 0, 0], [1, 0, 0]]
P[0][2] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[0][3] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

P[1][0] = [[0, 0, 0], [0, 0, 0], [1, 0, 0]]
P[1][1] = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
P[1][2] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[1][3] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]

P[2][0] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[2][1] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[2][2] = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]
P[2][3] = [[0, 0, 0], [0, 0, 0], [1, 0, 0]]

P[3][0] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[3][1] = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
P[3][2] = [[0, 0, 0], [0, 0, 0], [1, 0, 0]]
P[3][3] = [[0, 0, 1], [0, 0, 0], [0, 0, 0]]

model.link(P)

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

print('=== State Probabilities - System-Wide Aggregated Probabilities ===\n')
print('This example illustrates the calculation of probabilities via normalizing constants.\n')

# Set a custom initial state (all jobs at station 3)
n = [[0, 0, 0, 0],
     [0, 0, 0, 0],
     [N[0], N[1], N[2], N[3]]]

stations = model.getStations()
for i in range(M):
    stations[i].setState(n[i])

print(f'Query state: All {sum(N)} jobs at station 3 (Queue2)')
print(f'State configuration: {n}\n')

# Solver options
options = {
    'verbose': 1,
    'seed': 23000
}

# Compute system-wide aggregated probabilities
print('Computing getProbSysAggr() with multiple solvers:\n')

solver_ctmc = CTMC(model, options)
pr_ctmc = solver_ctmc.getProbSysAggr()
print(f'CTMC: Pr_ctmc = {pr_ctmc}')

solver_nc = NC(model, dict(options, method='exact'))
pr_nc = solver_nc.getProbSysAggr()
print(f'NC (exact): Pr_nc = {pr_nc}')

solver_jmt = JMT(model, samples=int(1e4), seed=532733)
pr_jmt = solver_jmt.getProbSysAggr()
print(f'JMT: Pr_jmt = {pr_jmt}')

print('\nNote: getProbSysAggr() returns the joint probability of the entire')
print('      system being in the specified aggregated state.')
print('      This is useful for analyzing the likelihood of specific system configurations.')
