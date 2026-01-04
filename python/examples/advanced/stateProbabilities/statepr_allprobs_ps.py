"""
State Probabilities - All Probability Types (PS)

This example demonstrates computation of various types of state probabilities
in a processor-sharing queueing network with class switching.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network model with PS queues
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
queue2.setNumberOfServers(2)

N = [2, 0]
job_class1 = ClosedClass(model, 'Class1', N[0], delay, 0)
job_class2 = ClosedClass(model, 'Class2', N[1], delay, 0)

delay.setService(job_class1, Exp(1))
delay.setService(job_class2, Exp(1))

queue1.setService(job_class1, Exp(3))
queue1.setService(job_class2, Exp(4))

queue2.setService(job_class1, Exp(1))
queue2.setService(job_class2, Exp(3))

K = model.getNumberOfClasses()
P = model.initRoutingMatrix()

P[0][0] = [[0, 1, 0], [0, 0, 0], [1, 0, 0]]
P[0][1] = [[0, 0, 0], [0, 0, 1], [0, 0, 0]]
P[1][0] = [[0, 1, 0], [0, 0, 0], [1, 0, 0]]
P[1][1] = [[0, 0, 0], [0, 0, 1], [0, 0, 0]]

model.link(P)

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

print('=== State Probabilities - All Probability Types (PS) ===\n')
print('This example illustrates the calculation of different types of probabilities.\n')

# Set a custom initial state
n = [[0, 0],
     [1, 0],
     [0, 1]]

stations = model.getStations()
for i in range(M):
    stations[i].setState(n[i])

state = model.getState()
print(f'Query state: {n}\n')

# Solver options
options = {
    'verbose': 1,
    'samples': int(2e4),
    'seed': 23000
}

# Test all probability methods with multiple solvers
print('=== get_prob_aggr (Marginal Aggregated Probabilities) ===')
print('Probability that station M is in the specified aggregated state.\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_aggr(stations[M-1])
print(f'CTMC: Pmarga = {pr}')

solver_nc = NC(model, options)
pr = solver_nc.get_prob_aggr(stations[M-1])
print(f'NC:   Pmarga = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_aggr(stations[M-1])
print(f'SSA:  Pmarga = {pr}')

solver_jmt = JMT(model, options)
pr = solver_jmt.get_prob_aggr(stations[M-1])
print(f'JMT:  Pmarga = {pr}')

print('\n=== get_prob (Marginal Detailed Probabilities) ===')
print('Marginal probabilities for the detailed state space (tracks service phases).\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob(stations[M-1])
print(f'CTMC: Pmarg = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob(stations[M-1])
print(f'SSA:  Pmarg = {pr}')

print('\n=== get_prob_sys_aggr (Joint Aggregated Probabilities) ===')
print('Joint state probabilities for the aggregated state space.\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_sys_aggr()
print(f'CTMC: Pjointa = {pr}')

solver_nc = NC(model, options)
pr = solver_nc.get_prob_sys_aggr()
print(f'NC:   Pjointa = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_sys_aggr()
print(f'SSA:  Pjointa = {pr}')

solver_jmt = JMT(model, options)
pr = solver_jmt.get_prob_sys_aggr()
print(f'JMT:  Pjointa = {pr}')

print('\n=== get_prob_sys (Joint Detailed Probabilities) ===')
print('Joint state probabilities for the detailed state space (with service phases).\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_sys()
print(f'CTMC: Pjoint = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_sys()
print(f'SSA:  Pjoint = {pr}')

print('\nNote: Probability types:')
print('  - Pmarga:  Marginal probability (aggregated state)')
print('  - Pmarg:   Marginal probability (detailed state with phases)')
print('  - Pjointa: Joint probability (aggregated state)')
print('  - Pjoint:  Joint probability (detailed state with phases)')
