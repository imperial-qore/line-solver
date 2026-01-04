"""
Cyclic Polling - Gated Service Policy

This example demonstrates a polling system with gated service policy,
where the server serves only those jobs present when polling begins.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create M[2]/M[2]/1-Exhaustive polling model with gated policy
model = Network('M[2]/M[2]/1-Exhaustive')

# Block 1: nodes
source = Source(model, 'mySource')
queue = Queue(model, 'myQueue', SchedStrategy.POLLING)
sink = Sink(model, 'mySink')

# Block 2: classes
oclass1 = OpenClass(model, 'myClass1')
source.setArrival(oclass1, Exp(1.0))
queue.setService(oclass1, Exp(4.0))

oclass2 = OpenClass(model, 'myClass2')
source.setArrival(oclass2, Exp(0.8))
queue.setService(oclass2, Exp(1.5))

# Set polling policy to gated with exponential switchover times
queue.setPollingType(PollingType.GATED)
queue.setSwitchover(oclass1, Exp(1.0))
queue.setSwitchover(oclass2, Exp(0.5))

# Block 3: topology
P = model.initRoutingMatrix()
P[0] = Network.serial_routing([source, queue, sink])
P[1] = Network.serial_routing([source, queue, sink])
model.link(P)

print('=== Cyclic Polling - Gated Service Policy ===\n')

# Solve with MVA (exact solution)
print('MVA Solution (exact):')
avg_table_mva = MVA(model).get_avg_table()
print(avg_table_mva)

# Solve with JMT simulation
print('\nJMT Simulation:')
avg_table_jmt = JMT(model, seed=23000, samples=int(1e6)).get_avg_table()
print(avg_table_jmt)

print('\nNote: In gated polling, the server serves only jobs present at the queue')
print('      when it begins polling. Jobs arriving during service must wait for')
print('      the next polling cycle.')
