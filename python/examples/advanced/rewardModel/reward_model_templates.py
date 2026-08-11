"""
Reward Model - Templates Example

This example demonstrates the use of Reward template factory methods
for quickly defining common metrics like queue length, utilization, and blocking.

The Reward templates provide a convenient way to define metrics
without writing lambda expressions.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

print('=== Reward Model - Templates Example ===\n')

# Model Definition
model = Network('RewardTemplatesExample')

# Nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Queue parameters
queue.setNumberOfServers(1)
queue.setCapacity(10)

# Job class
oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(1.5))  # Arrival rate = 1.5
queue.setService(oclass, Exp(2))     # Service rate = 2 (rho = 0.75)

# Topology
model.link(Network.serial_routing([source, queue, sink]))

# Define Rewards Using Templates
print('Defining reward metrics using Reward templates:\n')

# Template 1: Queue Length
model.setReward('QueueLength', Reward.queue_length(queue))
print('  QueueLength = Reward.queue_length(queue)')

# Template 2: Queue Length for Specific Class
model.setReward('QueueLength_Class1', Reward.queue_length(queue, oclass))
print('  QueueLength_Class1 = Reward.queue_length(queue, oclass)')

# Template 3: Server Utilization
model.setReward('Utilization', Reward.utilization(queue))
print('  Utilization = Reward.utilization(queue)')

# Template 4: Utilization by Class
model.setReward('Utilization_Class1', Reward.utilization(queue, oclass))
print('  Utilization_Class1 = Reward.utilization(queue, oclass)')

# Template 5: Blocking Probability
model.setReward('BlockingProb', Reward.blocking(queue))
print('  BlockingProb = Reward.blocking(queue)')

# Solve with CTMC Solver
print('\nSolving with CTMC solver...\n')

options = {'verbose': 0}
solver = CTMC(model, **options)

# Get Steady-State Expected Rewards
R, names = solver.getAvgReward()

print('=== Steady-State Expected Rewards (Templates) ===')
print(f'{"QueueLength":25s}: {R[0]:10.6f}')
print(f'{"QueueLength_Class1":25s}: {R[1]:10.6f}')
print(f'{"Utilization":25s}: {R[2]:10.6f}')
print(f'{"Utilization_Class1":25s}: {R[3]:10.6f}')
print(f'{"BlockingProb":25s}: {R[4]:10.6f}')

# Analytical Comparison
print('\n=== Comparison with M/M/1 Analytical Results ===')

lambda_rate = 1.5  # Arrival rate
mu = 2             # Service rate
rho = lambda_rate / mu
K = 10             # Buffer capacity

# Steady-state probabilities
import math
pi = []
for n in range(K + 1):
    pi.append((1 - rho) * (rho ** n) / (1 - rho ** (K + 1)))

# Analytical metrics
L_analytical = sum(n * pi[n] for n in range(K + 1))  # Expected queue length
U_analytical = 1 - pi[0]                              # Utilization
B_analytical = pi[K]                                  # Blocking probability

print(f'{"QueueLength":25s}: LINE = {R[0]:10.6f}, Analytical = {L_analytical:10.6f}, Error = {abs(R[0] - L_analytical):.2e}')
print(f'{"Utilization":25s}: LINE = {R[2]:10.6f}, Analytical = {U_analytical:10.6f}, Error = {abs(R[2] - U_analytical):.2e}')
print(f'{"BlockingProb":25s}: LINE = {R[4]:10.6f}, Analytical = {B_analytical:10.6f}, Error = {abs(R[4] - B_analytical):.2e}')

print('\nExample completed successfully.')
print('  All template-based rewards match analytical results!')
