"""
Switchover Times - Basic Example

This example demonstrates how to specify switchover times between job classes
in a queueing system. Switchover times represent the time required for a server
to switch from serving one class to another.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Create M[2]/M[2]/1-Gated model with switchover times
model = Network('M[2]/M[2]/1-Gated')

# Block 1: nodes
source = Source(model, 'mySource')
queue = Queue(model, 'myQueue', SchedStrategy.FCFS)
sink = Sink(model, 'mySink')

# Block 2: classes
oclass1 = OpenClass(model, 'myClass1')
source.setArrival(oclass1, Exp(0.2))
queue.setService(oclass1, Exp(0.1))

oclass2 = OpenClass(model, 'myClass2')
source.setArrival(oclass2, Exp(0.8))
queue.setService(oclass2, Exp(1.5))

# Set switchover times between classes
# Switchover from Class1 to Class2: Exponential with rate 1
queue.setSwitchover(oclass1, oclass2, Exp(1))
# Switchover from Class2 to Class1: Erlang with mean 1 and order 2
queue.setSwitchover(oclass2, oclass1, Erlang(1, 2))

# Block 3: topology
P = model.initRoutingMatrix()
P[0] = Network.serial_routing([source, queue, sink])
P[1] = Network.serial_routing([source, queue, sink])
model.link(P)

print('=== Switchover Times - Basic Example ===\n')

# Solve with JMT simulation
print('JMT Simulation:')
avg_table_jmt = JMT(model, seed=23000, keep=True).get_avg_table()
print(avg_table_jmt)

print('\nNote: Switchover times model the overhead of switching between job classes.')
print('      This is particularly important in systems where server reconfiguration')
print('      or setup times are non-negligible.')
print('\nSwitchover configuration:')
print('  - Class1 -> Class2: Exp(1) - Exponential with mean 1')
print('  - Class2 -> Class1: Erlang(1,2) - Erlang with mean 1 and order 2')
