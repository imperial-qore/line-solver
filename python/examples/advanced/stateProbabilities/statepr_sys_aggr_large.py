"""
State Probabilities - System-Wide Aggregated Probabilities (Large Model)

This example demonstrates computation of system-wide (joint) state probabilities
for a larger multi-class queueing network with class switching.

Expected result from MATLAB: Pr ≈ 0.000348

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *
import numpy as np
import time

# Create closed network model with 4 classes
model = Network('model')

# Block 1: nodes (3 queues, no delay)
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)
queue3 = Queue(model, 'Queue3', SchedStrategy.PS)
queue3.setNumberOfServers(3)

# Block 2: classes
N = [1, 1, 1, 1]
job_class1 = ClosedClass(model, 'Class1', N[0], queue1, 0)
job_class2 = ClosedClass(model, 'Class2', N[1], queue1, 0)
job_class3 = ClosedClass(model, 'Class3', N[2], queue1, 0)
job_class4 = ClosedClass(model, 'Class4', N[3], queue1, 0)

# Set service times for Queue1
queue1.setService(job_class1, Exp(1))
queue1.setService(job_class2, Exp(2))
queue1.setService(job_class3, Exp(1))
queue1.setService(job_class4, Exp(1))

# Set service times for Queue2
queue2.setService(job_class1, Exp(3))
queue2.setService(job_class2, Exp(4))
queue2.setService(job_class3, Exp(5))
queue2.setService(job_class4, Exp(1))

# Set service times for Queue3
queue3.setService(job_class1, Exp(1))
queue3.setService(job_class2, Exp(3))
queue3.setService(job_class3, Exp(5))
queue3.setService(job_class4, Exp(2))

# Block 3: routing with class switching
P = {}

# Class1 routing
P[(job_class1, job_class1)] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
P[(job_class1, job_class2)] = np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]])
P[(job_class1, job_class3)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class1, job_class4)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

# Class2 routing
P[(job_class2, job_class1)] = np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]])
P[(job_class2, job_class2)] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
P[(job_class2, job_class3)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class2, job_class4)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])

# Class3 routing
P[(job_class3, job_class1)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class3, job_class2)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class3, job_class3)] = np.array([[0, 1, 0], [0, 0, 1], [0, 0, 0]])
P[(job_class3, job_class4)] = np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]])

# Class4 routing
P[(job_class4, job_class1)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class4, job_class2)] = np.array([[0, 0, 0], [0, 0, 0], [0, 0, 0]])
P[(job_class4, job_class3)] = np.array([[0, 0, 0], [0, 0, 0], [1, 0, 0]])
P[(job_class4, job_class4)] = np.array([[0, 0, 1], [0, 0, 0], [0, 0, 0]])

model.link(P)

M = model.getNumberOfStations()
K = model.getNumberOfClasses()

print('=== State Probabilities - System-Wide Aggregated (Large Model) ===\n')
print('This example illustrates the calculation of probabilities via normalizing constants.\n')

# Set a custom initial state (all jobs at station 3)
n = [[0, 0, 0, 0],
     [0, 0, 0, 0],
     [N[0], N[1], N[2], N[3]]]

stations = model.getStations()
for i in range(M):
    stations[i].setState(n[i])

print(f'Query state: All {sum(N)} jobs at station 3 (Queue3)')
print(f'State configuration: {n}\n')

# Solver options
options = {'verbose': 1, 'seed': 23000}

print('Computing getProbSysAggr() with CTMC solver:\n')

# CTMC solver
start = time.time()
solver_ctmc = CTMC(model, options)
pr_ctmc = solver_ctmc.getProbSysAggr()
print(f'CTMC: Pr_ctmc = {pr_ctmc}')
print('Pr_ctmc =')
print(pr_ctmc)
print(f'CTMC time: {time.time() - start:.3f}s')
print(f'(Expected MATLAB result: ~0.000348)\n')

print('Note: getProbSysAggr() returns the joint probability of the entire')
print('      system being in the specified aggregated state.')
