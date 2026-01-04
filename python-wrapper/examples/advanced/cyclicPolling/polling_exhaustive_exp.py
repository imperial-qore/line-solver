"""
Cyclic Polling - Exhaustive Service with Exponential Distributions

This example demonstrates a polling system with exhaustive service policy,
where the server continues serving a queue until it becomes empty.

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

# Block 2: classes
oclass1 = OpenClass(model, 'myClass1')
source.setArrival(oclass1, Exp(0.1))
queue.setService(oclass1, Exp(1.0))

oclass2 = OpenClass(model, 'myClass2')
source.setArrival(oclass2, Exp(0.1))
queue.setService(oclass2, Exp(1.5))

# Set polling policy to exhaustive
queue.setPollingType(PollingType.EXHAUSTIVE)

# Block 3: topology
P = model.initRoutingMatrix()
P[0] = Network.serial_routing([source, queue, sink])
P[1] = Network.serial_routing([source, queue, sink])
model.link(P)

print('=== Cyclic Polling - Exhaustive Service (Exponential) ===\n')

# Solve with MVA (approximate solution in general)
print('MVA Solution (approximate in general):')
avg_table_mva = MVA(model).get_avg_table()
print(avg_table_mva)

# Solve with JMT simulation
print('\nJMT Simulation:')
avg_table_jmt = JMT(model, seed=23000).get_avg_table()
print(avg_table_jmt)

print('\nNote: In exhaustive polling, the server continues serving jobs from a queue')
print('      until it becomes empty before switching to the next queue.')
