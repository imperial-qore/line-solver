"""
Reward Model - Multi-Class Example

This example demonstrates reward-based CTMC analysis in a multi-class queueing
network, showing how to define per-class metrics and analyze class-specific behavior.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

print('=== Reward Model - Multi-Class Example ===\n')

# Model Definition
model = Network('MultiClassRewardExample')

# Nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.PS)  # Processor Sharing
sink = Sink(model, 'Sink')

# Server capacity and buffer limit
queue.setNumberOfServers(2)
queue.setCapacity(6)  # Limit state space

# Job classes with different characteristics
class_interactive = OpenClass(model, 'Interactive')
class_batch = OpenClass(model, 'Batch')

# Different arrival rates
source.setArrival(class_interactive, Exp(2.0))  # Interactive: lambda = 2.0
source.setArrival(class_batch, Exp(1.5))        # Batch: lambda = 1.5

# Different service time requirements
queue.setService(class_interactive, Exp(0.5))   # Interactive: mu = 2.0 (fast)
queue.setService(class_batch, Exp(1.0))         # Batch: mu = 1.0 (slow)

# Topology (same for both classes)
P = model.initRoutingMatrix()
P[class_interactive] = Network.serial_routing([source, queue, sink])
P[class_batch] = Network.serial_routing([source, queue, sink])
model.link(P)

# Define Per-Class Rewards
print('Defining per-class reward metrics:\n')

# === Interactive Class Metrics ===

# Jobs of Interactive class
model.setReward('Interactive_QLen', lambda state: state.at(queue, class_interactive))
print('  Interactive_QLen = state.at(queue, Interactive)')

# Utilization contributed by Interactive class
model.setReward('Interactive_Util', Reward.utilization(queue, class_interactive))
print('  Interactive_Util = Reward.utilization(queue, Interactive)')

# === Batch Class Metrics ===

# Jobs of Batch class
model.setReward('Batch_QLen', lambda state: state.at(queue, class_batch))
print('  Batch_QLen = state.at(queue, Batch)')

# Utilization contributed by Batch class
model.setReward('Batch_Util', Reward.utilization(queue, class_batch))
print('  Batch_Util = Reward.utilization(queue, Batch)')

# === Comparative Metrics ===

# Total queue length (both classes)
model.setReward('Total_QLen', lambda state: state.at(queue).total())
print('  Total_QLen = state.at(queue).total()')

# Total system utilization
model.setReward('Total_Util', Reward.utilization(queue))
print('  Total_Util = Reward.utilization(queue)')

# Ratio of Interactive to Batch jobs
model.setReward('Interactive_Ratio', lambda state:
    state.at(queue, class_interactive) / max(state.at(queue, class_batch), 0.001))
print('  Interactive_Ratio = Interactive / Batch')

# === Priority-Aware Metrics ===

# Weighted response time cost
model.setReward('Weighted_Cost', lambda state:
    3.0 * state.at(queue, class_interactive) +
    1.0 * state.at(queue, class_batch))
print('  Weighted_Cost = 3.0*Interactive + 1.0*Batch')

# Service fairness indicator
model.setReward('Fairness', lambda state:
    float(state.at(queue, class_interactive) > 0 and state.at(queue, class_batch) > 0))
print('  Fairness = (Interactive > 0) AND (Batch > 0)')

# Solve with CTMC Solver
print('\nSolving with CTMC solver...\n')

options = {'verbose': 0, 'cutoff': 6}
solver = CTMC(model, **options)

# Get Steady-State Expected Rewards
R, names = solver.getAvgReward()

print('=== Steady-State Expected Rewards (Multi-Class) ===\n')

print('Interactive Class:')
print(f'  {"Interactive_QLen":25s}: {R[0]:10.6f} (jobs)')
print(f'  {"Interactive_Util":25s}: {R[1]:10.6f} (server capacity)')

print('\nBatch Class:')
print(f'  {"Batch_QLen":25s}: {R[2]:10.6f} (jobs)')
print(f'  {"Batch_Util":25s}: {R[3]:10.6f} (server capacity)')

print('\nComparative Metrics:')
print(f'  {"Total_QLen":25s}: {R[4]:10.6f} (both classes)')
print(f'  {"Total_Util":25s}: {R[5]:10.6f} (2 servers)')
print(f'  {"Interactive_Ratio":25s}: {R[6]:10.6f} (ratio)')

print('\nPriority-Aware Metrics:')
print(f'  {"Weighted_Cost":25s}: {R[7]:10.6f} (weighted cost)')
print(f'  {"Fairness":25s}: {R[8]:10.6f} (fraction of time both present)')

# Analysis
print('\n=== Analysis ===')

lambda_int = 2.0
mu_int = 2.0
lambda_batch = 1.5
mu_batch = 1.0
c = 2  # Number of servers

rho_int = lambda_int / (c * mu_int)     # 2 / (2*2) = 0.5
rho_batch = lambda_batch / (c * mu_batch)  # 1.5 / (2*1) = 0.75
rho_total = rho_int + rho_batch          # 1.25

print('System Characteristics:')
print(f'  Interactive: lambda={lambda_int:.1f}, mu={mu_int:.1f} per server')
print(f'  Batch:       lambda={lambda_batch:.1f}, mu={mu_batch:.1f} per server')
print(f'  Servers: {c}')
print(f'  Total utilization: {rho_total:.3f} ({100*rho_total:.1f}% capacity)')

# Class composition
total_qlen = R[4]
if total_qlen > 0:
    pct_int = 100 * R[0] / total_qlen
    pct_batch = 100 * R[2] / total_qlen
    print('\nQueue Composition:')
    print(f'  Interactive jobs: {pct_int:.2f}% of queue')
    print(f'  Batch jobs: {pct_batch:.2f}% of queue')

# Response time estimation (Little's Law)
total_arrival = lambda_int + lambda_batch
est_resp_time_int = R[0] / lambda_int
est_resp_time_batch = R[2] / lambda_batch

print("\nEstimated Response Times (Little's Law):")
print(f'  Interactive: {est_resp_time_int:.6f} time units')
print(f'  Batch: {est_resp_time_batch:.6f} time units')

print('\nExample completed successfully.')
