"""
State Probabilities - All Probability Types (FCFS)

This example demonstrates computation of various types of state probabilities:
- Marginal aggregated probabilities (get_prob_aggr)
- Marginal detailed probabilities (get_prob)
- Joint aggregated probabilities (get_prob_sys_aggr)
- Joint detailed probabilities (get_prob_sys)

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network model with class switching
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
queue2.setNumberOfServers(2)

N = [2, 0]
job_class1 = ClosedClass(model, 'Class1', N[0], delay, 0)
job_class2 = ClosedClass(model, 'Class2', N[1], delay, 0)

delay.setService(job_class1, Exp(1))
delay.setService(job_class2, Exp(1))

queue1.setService(job_class1, Exp(3))
queue1.setService(job_class2, Exp(4))

queue2.setService(job_class1, Exp(3))
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

print('=== State Probabilities - All Probability Types (FCFS) ===\n')
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
    'samples': int(1e5),
    'seed': 23000
}

# Test all probability methods with multiple solvers
print('=== get_prob_aggr (Marginal Aggregated Probabilities) ===')
print('Probability that station M is in the specified aggregated state.\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_aggr(stations[M-1])
print(f'CTMC: Station {M} is in state {state[M-1]} with probability {pr}')

solver_nc = NC(model, options)
pr = solver_nc.get_prob_aggr(stations[M-1])
print(f'NC:   Station {M} is in state {state[M-1]} with probability {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_aggr(stations[M-1])
print(f'SSA:  Station {M} is in state {state[M-1]} with probability {pr}')

solver_jmt = JMT(model, options)
pr = solver_jmt.get_prob_aggr(stations[M-1])
print(f'JMT:  Station {M} is in state {state[M-1]} with probability {pr}')

print('\n=== get_prob (Marginal Detailed Probabilities) ===')
print('Detailed probability including service phases.\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob(stations[M-1])
print(f'CTMC: Detailed probability = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob(stations[M-1])
print(f'SSA:  Detailed probability = {pr}')

print('\n=== get_prob_sys_aggr (Joint Aggregated Probabilities) ===')
print('Joint probability for the aggregated system state.\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_sys_aggr()
print(f'CTMC: Joint aggregated probability = {pr}')

solver_nc = NC(model, options)
pr = solver_nc.get_prob_sys_aggr()
print(f'NC:   Joint aggregated probability = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_sys_aggr()
print(f'SSA:  Joint aggregated probability = {pr}')

solver_jmt = JMT(model, options)
pr = solver_jmt.get_prob_sys_aggr()
print(f'JMT:  Joint aggregated probability = {pr}')

print('\n=== get_prob_sys (Joint Detailed Probabilities) ===')
print('Joint probability for the detailed system state (with service phases).\n')

solver_ctmc = CTMC(model, options)
pr = solver_ctmc.get_prob_sys()
print(f'CTMC: Joint detailed probability = {pr}')

solver_ssa = SSA(model, options)
pr = solver_ssa.get_prob_sys()
print(f'SSA:  Joint detailed probability = {pr}')

print('\nNote: Aggregated vs. Detailed state spaces:')
print('  - Aggregated: (n_i1, ..., n_iR) - job counts per class at each station')
print('  - Detailed: Also tracks service phases for non-exponential distributions')
print('  - Marginal: Probability for a single station')
print('  - Joint (Sys): Probability for the entire system state')
