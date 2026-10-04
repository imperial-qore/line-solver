"""
Reward Model - Aggregation Example

This example demonstrates the aggregation capabilities of the RewardState API.

The RewardStateView class supports aggregation operations:
  - total()  : Sum all values
  - max()    : Maximum value
  - min()    : Minimum value
  - count()  : Count non-zero entries

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

print('=== Reward Model - Aggregation Example ===\n')

# Model Definition
model = Network('RewardAggregationExample')

# Nodes
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Queue parameters
queue.setNumberOfServers(1)
queue.setCapacity(4)  # Limit state space

# Job classes
class1 = OpenClass(model, 'HighPriority')
class2 = OpenClass(model, 'LowPriority')

# Arrivals
source.setArrival(class1, Exp(1.0))    # High priority arrival rate = 1.0
source.setArrival(class2, Exp(0.8))    # Low priority arrival rate = 0.8

# Service times
queue.setService(class1, Exp(3.0))
queue.setService(class2, Exp(3.0))

# Routing (same for both classes)
P = model.initRoutingMatrix()
P[class1] = Network.serial_routing([source, queue, sink])
P[class2] = Network.serial_routing([source, queue, sink])
model.link(P)

# Define Rewards Using Aggregation
print('Defining reward metrics using aggregation operations:\n')

# Total jobs at Queue (all classes)
model.setReward('TotalJobs', lambda state: state.at(queue).total())
print('  TotalJobs = state.at(queue).total()')

# Maximum population of any class at Queue
model.setReward('MaxClass', lambda state: state.at(queue).max())
print('  MaxClass = state.at(queue).max()')

# Count of classes with jobs at Queue
model.setReward('ClassCount', lambda state: state.at(queue).count())
print('  ClassCount = state.at(queue).count()')

# HighPriority jobs
model.setReward('HP_Jobs', lambda state: state.at(queue, class1))
print('  HP_Jobs = state.at(queue, class1)')

# LowPriority jobs
model.setReward('LP_Jobs', lambda state: state.at(queue, class2))
print('  LP_Jobs = state.at(queue, class2)')

# Weighted load (high priority weighted 2x)
model.setReward('WeightedLoad', lambda state:
    2.0 * state.at(queue, class1) + 1.0 * state.at(queue, class2))
print('  WeightedLoad = 2.0 * HP + 1.0 * LP')

# Solve with CTMC Solver
print('\nSolving with CTMC solver...\n')

options = {'verbose': 0}
solver = CTMC(model, **options)

# Get Steady-State Expected Rewards
R, names = solver.getAvgReward()

print('=== Steady-State Expected Rewards (Aggregation) ===')
print(f'  {"TotalJobs":20s}: {R[0]:10.6f}')
print(f'  {"MaxClass":20s}: {R[1]:10.6f}')
print(f'  {"ClassCount":20s}: {R[2]:10.6f}')
print(f'  {"HP_Jobs":20s}: {R[3]:10.6f}')
print(f'  {"LP_Jobs":20s}: {R[4]:10.6f}')
print(f'  {"WeightedLoad":20s}: {R[5]:10.6f}')

print('\nExample completed successfully.')
