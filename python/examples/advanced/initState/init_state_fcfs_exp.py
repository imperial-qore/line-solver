"""
Initial State - FCFS Queue with Exponential Service

This example demonstrates how to specify initial states in a closed queueing
network with FCFS scheduling and exponential service. Different priors on
the initial state distribution affect transient behavior.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create closed model with FCFS queue
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue1', SchedStrategy.FCFS)

job_class = ClosedClass(model, 'Class1', 5, queue, 0)

delay.setService(job_class, Exp(1))
queue.setService(job_class, Exp(0.7))

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

# Circular routing
P = Network.serial_routing([delay, queue])
model.link(P)

# Get transient handles
Qt, Ut, Tt = model.getTranHandles()

# Solver options
options = {
    'verbose': 0,
    'samples': int(1e4),
    'stiff': True,
    'timespan': [0, 40]
}

print('=== Initial State - FCFS Queue with Exponential Service ===\n')

# Prior 1: Default initialization
print('--- Prior 1: Default initialization ---')
model.initDefault()
state = model.getState()
print('Initial state is:')
print(f'Station 1 (Delay): {state[0][0]}')
print(f'Station 2 (Queue): {state[1][0]}')

solver_ctmc = CTMC(model, options)
QNt_ctmc, UNt_ctmc, TNt_ctmc = solver_ctmc.getTranAvg(Qt, Ut, Tt)

solver_fluid = FLD(model, options)
QNt_fluid, UNt_fluid, TNt_fluid = solver_fluid.getTranAvg(Qt, Ut, Tt)

print('\nTransient queue length at station 2, class 1 computed.')
print(f'CTMC: {len(QNt_ctmc[1][0].t)} time points')
print(f'FLD: {len(QNt_fluid[1][0].t)} time points')

# Prior 2: Prior on state with given marginal (3 jobs in station 2)
print('\n--- Prior 2: Prior on state with 3 jobs in station 2 ---')
model.initFromMarginal([2, 3])
state = model.getState()
print('Initial state is:')
print(f'Station 1 (Delay): {state[0][0]}')
print(f'Station 2 (Queue): {state[1][0]}')

solver_ctmc.reset()
QNt_marg_ctmc, UNt_marg_ctmc, TNt_marg_ctmc = solver_ctmc.getTranAvg(Qt, Ut, Tt)

solver_fluid.reset()
QNt_marg_fluid, UNt_marg_fluid, TNt_marg_fluid = solver_fluid.getTranAvg(Qt, Ut, Tt)

print('\nTransient queue length at station 2, class 1 computed with marginal prior.')
print(f'CTMC: {len(QNt_marg_ctmc[1][0].t)} time points')
print(f'FLD: {len(QNt_marg_fluid[1][0].t)} time points')

print('\nNote: Different initial state priors lead to different transient behaviors.')
print('      Default initialization starts with all jobs at their reference station.')
print('      Marginal initialization sets a specific number of jobs at each station.')
