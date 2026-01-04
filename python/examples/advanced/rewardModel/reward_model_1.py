"""
Reward Model - Basic Example

This example demonstrates the setReward feature for defining custom
reward functions on a queueing network model and computing steady-state
expected rewards using the CTMC solver.

Reward functions allow modeling various performance metrics such as:
- Queue lengths
- Server utilization
- Blocking probabilities
- Custom cost functions

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *
import math

# Model Definition
# Create a simple M/M/1/K queue with finite buffer
model = Network('RewardExample')

# Block 1: nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Set finite buffer capacity
queue.setNumberOfServers(1)
queue.setCapacity(3)  # Maximum jobs in the system

# Block 2: job classes
oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(2))  # Arrival rate = 2
queue.setService(oclass, Exp(3))   # Service rate = 3 (utilization ~ 0.67)

# Block 3: topology
model.link(Network.serial_routing([source, queue, sink]))

print('=== Reward Model - Basic Example ===\n')

# Define Reward Functions
# setReward(name, rewardFn) where rewardFn uses state accessor methods

# Reward 1: Queue length (number of jobs in the queue)
model.setReward('QueueLength', lambda state: state.at(queue, oclass))

# Reward 2: Utilization (1 if server busy, 0 if idle)
model.setReward('Utilization', Reward.utilization(queue, oclass))

# Reward 3: Blocking indicator (1 if buffer full, 0 otherwise)
model.setReward('BlockingProb', Reward.blocking(queue))

# Reward 4: Weighted queue cost (quadratic penalty for long queues)
model.setReward('QueueCost', lambda state: state.at(queue, oclass) ** 2)

# Solve with CTMC Solver
print('Solving with CTMC solver...\n')

options = {
    'verbose': 1
}

solver = CTMC(model, options)

# Get Steady-State Expected Rewards
R, names = solver.get_avg_reward()

print('\n=== Steady-State Expected Rewards ===')
for i in range(len(names)):
    print(f'{names[i]:>15s}: {R[i]:.6f}')

# Get Transient Reward Analysis (if available)
try:
    t, V, names, state_space = solver.get_reward()
    print(f'\n=== State Space ===')
    print(f'Number of states: {len(state_space)}')
    if len(state_space) > 0:
        print(f'State dimensions: {len(state_space[0])}')
except:
    print('\n=== Transient Reward Analysis ===')
    print('Transient analysis not available or not yet implemented.')

# Compare with Analytical Results for M/M/1/K
print('\n=== Comparison with M/M/1/K Analytical Results ===')

lambda_rate = 2.0  # Arrival rate
mu = 3.0          # Service rate
rho = lambda_rate / mu
K = 10            # Buffer capacity

# For M/M/1/K queue, steady-state probabilities:
# pi(n) = (1-rho) * rho^n / (1 - rho^(K+1))
if rho != 1:
    pi = []
    for n in range(K + 1):
        pi.append((1 - rho) * (rho ** n) / (1 - rho ** (K + 1)))
else:
    pi = [1.0 / (K + 1)] * (K + 1)

# Analytical expected queue length
L_analytical = sum(n * pi[n] for n in range(K + 1))

# Analytical utilization (P(server busy) = 1 - pi(0))
U_analytical = 1 - pi[0]

# Analytical blocking probability = pi(K)
B_analytical = pi[K]

# Analytical queue cost (E[N^2])
Cost_analytical = sum((n ** 2) * pi[n] for n in range(K + 1))

print(f'{"QueueLength":>15s}: LINE = {R[0]:.6f}, Analytical = {L_analytical:.6f}, Error = {abs(R[0] - L_analytical):.2e}')
print(f'{"Utilization":>15s}: LINE = {R[1]:.6f}, Analytical = {U_analytical:.6f}, Error = {abs(R[1] - U_analytical):.2e}')
print(f'{"BlockingProb":>15s}: LINE = {R[2]:.6f}, Analytical = {B_analytical:.6f}, Error = {abs(R[2] - B_analytical):.2e}')
print(f'{"QueueCost":>15s}: LINE = {R[3]:.6f}, Analytical = {Cost_analytical:.6f}, Error = {abs(R[3] - Cost_analytical):.2e}')

print('\nNote: Reward functions provide a powerful way to define custom performance')
print('      metrics beyond standard measures like utilization and queue length.')
print('      The CTMC solver can compute steady-state expected values and transient')
print('      behavior of these reward functions.')
