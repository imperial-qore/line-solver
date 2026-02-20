"""
State Probabilities - Aggregated State Space

This example demonstrates how to compute marginal state probabilities for
the aggregated state space using normalizing constants (NC) and CTMC solvers.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed network model
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
P = [[None, None], [None, None]]

P[0][0] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]
P[1][1] = [[0, 1, 0], [0, 0, 1], [1, 0, 0]]

model.link(P)

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

print('=== State Probabilities - Aggregated State Space ===\n')
print('This example illustrates the calculation of probabilities via normalizing constants.\n')

# Set a custom state to query
# -1 means ignore that station in probability calculation
n = [[-1, -1],
     [-1, -1],
     [0, 0]]

stations = model.getStations()
for i in range(M):
    if n[i][0] != -1 or n[i][1] != -1:
        stations[i].setState(n[i])

# Solver options
options = {
    'verbose': 1,
    'seed': 23000
}

# Query probability of state at station M (queue2)
i = M - 1  # Python 0-indexing
target_station = stations[i]
target_state = n[i]

print(f'Querying probability for station {i+1} ({target_station.getName()}) in state {target_state}\n')

# Solve with CTMC
solver_ctmc = CTMC(model, options)
pr_ctmc = solver_ctmc.getProbAggr(target_station)
print(f'CTMC: Station {i+1} is in state {target_state} with probability {pr_ctmc}')

# Solve with NC
solver_nc = NC(model, options)
pr_nc = solver_nc.getProbAggr(target_station)
print(f'NC:   Station {i+1} is in state {target_state} with probability {pr_nc}')

print('\nNote: getProbAggr() computes marginal probabilities for the aggregated')
print('      state space, where each station is specified by a tuple (n_i1, ..., n_iR),')
print('      with n_ir being the number of class r jobs at station i.')
print('      Rows set to [-1, -1] are ignored in the calculation.')
