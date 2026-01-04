"""
Cyclic Polling - Exhaustive Service with Deterministic Distributions

This example demonstrates a polling system with exhaustive service policy
using deterministic (constant) arrival and service times.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create M[2]/M[2]/1-Gated polling model
model = Network('M[2]/M[2]/1-Gated')

# Block 1: nodes
source = Source(model, 'mySource')
queue = Queue(model, 'myQueue', SchedStrategy.POLLING)
sink = Sink(model, 'mySink')

# Block 2: classes with deterministic arrivals and services
oclass1 = OpenClass(model, 'myClass1')
source.setArrival(oclass1, Det.fit_mean(1.0))
queue.setService(oclass1, Det.fit_mean(0.001))

oclass2 = OpenClass(model, 'myClass2')
source.setArrival(oclass2, Det.fit_mean(1.0))
queue.setService(oclass2, Det.fit_mean(0.001))

# Set polling policy to exhaustive with immediate switchover
queue.setPollingType(PollingType.EXHAUSTIVE)
queue.setSwitchover(oclass1, Immediate())
queue.setSwitchover(oclass2, Immediate())

# Block 3: topology
P = model.initRoutingMatrix()
P[0] = Network.serial_routing([source, queue, sink])
P[1] = Network.serial_routing([source, queue, sink])
model.link(P)

print('=== Cyclic Polling - Exhaustive Service (Deterministic) ===\n')

# Solve with MVA (approximate solution in general)
print('MVA Solution (approximate in general):')
avg_table_mva = MVA(model).get_avg_table()
print(avg_table_mva)

# Solve with JMT simulation
print('\nJMT Simulation:')
avg_table_jmt = JMT(model, seed=23000, samples=int(1e5)).get_avg_table()
print(avg_table_jmt)

print('\nNote: This example uses deterministic (constant) arrivals and service times')
print('      with immediate (zero-time) switchovers between queues.')
