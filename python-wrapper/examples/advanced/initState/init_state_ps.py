"""
Initial State - PS Queue with Class Switching

This example demonstrates initial state specification for a processor-sharing
queue with two classes and class-switching behavior. Uses marginal and started
job specification.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with PS queue and class switching
model = Network('model')

delay = Delay(model, 'InfiniteServer')
queue = Queue(model, 'Queue1', SchedStrategy.PS)
queue.setNumberOfServers(2)

job_class1 = ClosedClass(model, 'Class1', 3, delay, 0)
job_class2 = ClosedClass(model, 'Class2', 1, delay, 0)

delay.setService(job_class1, Exp(3))
delay.setService(job_class2, Exp(0.5))

queue.setService(job_class1, Exp(0.1))
queue.setService(job_class2, Exp(1))

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

# Class-switching routing matrix
P = [[None, None], [None, None]]
P[0][0] = [[0.3, 0.1], [0.2, 0]]
P[0][1] = [[0.6, 0], [0.8, 0]]
P[1][1] = [[0, 1], [0, 0]]
P[1][0] = [[0, 0], [1, 0]]

model.link(P)

# Initialize from marginal and started jobs
# Marginal: [2,1; 1,0] means station 1 has 2 class1 and 1 class2 jobs,
#                      station 2 has 1 class1 and 0 class2 jobs
# Started: [0,0; 1,0] means at station 2, 1 class1 job has started service
model.initFromMarginalAndStarted([[2, 1], [1, 0]], [[0, 0], [1, 0]])

# Solver options
options = {
    'verbose': 1,
    'seed': 23000,
    'samples': int(1e5)
}

print('=== Initial State - PS Queue with Class Switching ===\n')
print('This example shows solver execution on a 2-class 2-node class-switching model')
print('with specified initial state.\n')

print('Initial state specification:')
print('  Marginal: [[2,1], [1,0]] - Job distribution across stations and classes')
print('  Started:  [[0,0], [1,0]] - Number of jobs in service at each station\n')

# Solve with multiple solvers
print('CTMC Solver:')
solver_ctmc = CTMC(model, options)
avg_table_ctmc = solver_ctmc.get_avg_table()
print(avg_table_ctmc)

print('\nJMT Simulation:')
solver_jmt = JMT(model, options)
avg_table_jmt = solver_jmt.get_avg_table()
print(avg_table_jmt)

print('\nSSA (Stochastic State-space Analysis):')
solver_ssa = SSA(model, options)
avg_table_ssa = solver_ssa.get_avg_table()
print(avg_table_ssa)

print('\nFLD (Fluid) Solver:')
solver_fluid = FLD(model, options)
avg_table_fluid = solver_fluid.get_avg_table()
print(avg_table_fluid)

print('\nMVA Solver:')
solver_mva = MVA(model, options)
avg_table_mva = solver_mva.get_avg_table()
print(avg_table_mva)

print('\nNC (Normalizing Constant) Solver:')
solver_nc = NC(model, options)
avg_table_nc = solver_nc.get_avg_table()
print(avg_table_nc)

print('\nNote: Initial state affects transient behavior but not steady-state metrics.')
print('      MVA and NC compute steady-state only, so initial state has no effect.')
print('      CTMC, JMT, SSA, and FLD can compute transient behavior from the initial state.')
