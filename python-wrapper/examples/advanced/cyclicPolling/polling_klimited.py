"""
Cyclic Polling - K-Limited Service Policy

This example demonstrates a polling system with k-limited service policy,
where the server serves at most k jobs from a queue before switching.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create M[2]/M[2]/1-Gated polling model with k-limited policy
model = Network('M[2]/M[2]/1-Gated')

# Block 1: nodes
source = Source(model, 'mySource')
queue = Queue(model, 'myQueue', SchedStrategy.POLLING)
sink = Sink(model, 'mySink')

# Block 2: classes
oclass1 = OpenClass(model, 'myClass1')
source.setArrival(oclass1, Exp(0.2))
queue.setService(oclass1, Exp(1.0))

oclass2 = OpenClass(model, 'myClass2')
source.setArrival(oclass2, Exp(0.3))
queue.setService(oclass2, Exp(1.5))

# Set polling policy to k-limited with k=1 (serve at most 1 job per visit)
queue.setPollingType(PollingType.KLIMITED, 1)
queue.setSwitchover(oclass1, Exp(1))
queue.setSwitchover(oclass2, Immediate())

# Block 3: topology
P = model.initRoutingMatrix()
P[0] = Network.serial_routing([source, queue, sink])
P[1] = Network.serial_routing([source, queue, sink])
model.link(P)

print('=== Cyclic Polling - K-Limited Service Policy ===\n')

# Solve with MVA (approximate solution)
print('MVA Solution (approximate):')
avg_table_mva = MVA(model).get_avg_table()
print(avg_table_mva)

# Solve with JMT simulation
print('\nJMT Simulation:')
avg_table_jmt = JMT(model, samples=int(1e5), seed=23000).get_avg_table()
print(avg_table_jmt)

print('\nNote: In k-limited polling with k=1, the server serves at most 1 job')
print('      from each queue before switching to the next queue.')
print('      This prevents a busy queue from monopolizing the server.')
