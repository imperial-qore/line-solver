"""
Initial State - FCFS Queue with Non-Exponential Service

This example demonstrates initial state specification with non-exponential
service distributions (Erlang). Shows three different prior distributions
on the initial state.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with FCFS queue and non-exponential service
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue.setNumberOfServers(3)

job_class1 = ClosedClass(model, 'Class1', 3, queue, 0)
job_class2 = ClosedClass(model, 'Class2', 2, queue, 0)

delay.setService(job_class1, Exp(1))
delay.setService(job_class2, Exp(1))
queue.setService(job_class1, Exp(1.2))
queue.setService(job_class2, Erlang.fit_mean_and_scv(1.0, 0.5))

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

# Class-switching routing matrix
P = [[None, None], [None, None]]
P[0][0] = [[0.3, 0.1], [0.2, 0]]
P[0][1] = [[0.6, 0], [0.8, 0]]
P[1][1] = [[0, 1], [0, 0]]
P[1][0] = [[0, 0], [1, 0]]

model.link(P)

# Get transient handles
Qt, Ut, Tt = model.getTranHandles()

# Solver options
options = {
    'verbose': 1,
    'samples': int(1e4),
    'stiff': True,
    'timespan': [0, 5]
}

print('=== Initial State - FCFS Queue with Non-Exponential Service ===\n')

# Initialize solvers
solver_ctmc = CTMC(model, options)
solver_fluid = FLD(model, options)

# Prior 1: Default initialization
print('--- Prior 1: Default initialization ---')
model.initDefault()
sn = model.getStruct()
print('Initial state is:')
print(f'Class 1 state: {sn.space[0][0]}')
print(f'Class 2 state: {sn.space[1][0]}')

QNt_ctmc_1, _, _ = solver_ctmc.getTranAvg(Qt, Ut, Tt)
solver_ctmc.reset()

QNt_fluid_1, _, _ = solver_fluid.getTranAvg(Qt, Ut, Tt)
solver_fluid.reset()

print(f'CTMC queue length at t=0: {QNt_ctmc_1[1][0].metric[0] if len(QNt_ctmc_1[1][0].metric) > 0 else "N/A"}')
print(f'FLD queue length at t=0: {QNt_fluid_1[1][0].metric[0] if len(QNt_fluid_1[1][0].metric) > 0 else "N/A"}')

# Prior 2: Prior on first state with given marginal
print('\n--- Prior 2: First state with marginal [0,0; 4,1] ---')
model.initFromMarginal([[0, 0], [4, 1]])
sn = model.getStruct()
print('Initial state is:')
print(f'Class 1 state: {sn.space[0][0]}')
print(f'Class 2 state: {sn.space[1][0]}')

solver_ctmc.reset()
QNt_ctmc_2, _, _ = solver_ctmc.getTranAvg(Qt, Ut, Tt)
solver_ctmc.reset()

solver_fluid.reset()
QNt_fluid_2, _, _ = solver_fluid.getTranAvg(Qt, Ut, Tt)
solver_fluid.reset()

print(f'CTMC queue length at t=0: {QNt_ctmc_2[1][0].metric[0] if len(QNt_ctmc_2[1][0].metric) > 0 else "N/A"}')
print(f'FLD queue length at t=0: {QNt_fluid_2[1][0].metric[0] if len(QNt_fluid_2[1][0].metric) > 0 else "N/A"}')

# Prior 3: Uniform prior over all states with the same marginal
print('\n--- Prior 3: Uniform prior over states with marginal [0,0; 4,1] ---')
model.initFromMarginal([[0, 0], [4, 1]])
sn = model.getStruct()
print('Initial states include:')
print(f'Number of possible states: {len(sn.space[1])}')

# Set uniform prior
prior = queue.getStatePrior()
if prior is not None:
    uniform_prior = [1.0 / len(prior)] * len(prior)
    queue.setStatePrior(uniform_prior)
    print(f'Set uniform prior over {len(uniform_prior)} states')

solver_ctmc.reset()
QNt_ctmc_3, _, _ = solver_ctmc.getTranAvg(Qt, Ut, Tt)

solver_fluid.reset()
QNt_fluid_3, _, _ = solver_fluid.getTranAvg(Qt, Ut, Tt)

print(f'CTMC queue length at t=0: {QNt_ctmc_3[1][0].metric[0] if len(QNt_ctmc_3[1][0].metric) > 0 else "N/A"}')
print(f'FLD queue length at t=0: {QNt_fluid_3[1][0].metric[0] if len(QNt_fluid_3[1][0].metric) > 0 else "N/A"}')

print('\nNote: This example shows three types of initial state specification:')
print('  1. Default: All jobs at reference stations')
print('  2. Marginal: First state matching the marginal distribution')
print('  3. Uniform: Uniform distribution over all states with given marginal')
