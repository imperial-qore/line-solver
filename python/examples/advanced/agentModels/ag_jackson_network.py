"""
Jackson Network with MAM

This example demonstrates MAM with RCAT methods on an open Jackson network with
probabilistic routing. The RCAT algorithm models job transfers between
queues as synchronization actions between interacting processes.

Reference: Jackson, J.R. (1957). "Networks of waiting lines"

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
lam = 1.0   # External arrival rate
mu = [2.0, 3.0, 2.5]  # Service rates

# Routing probabilities
p12 = 0.4  # Queue1 -> Queue2
p13 = 0.3  # Queue1 -> Queue3
p1s = 0.3  # Queue1 -> Sink
p21 = 0.2  # Queue2 -> Queue1
p23 = 0.3  # Queue2 -> Queue3
p2s = 0.5  # Queue2 -> Sink
p3s = 1.0  # Queue3 -> Sink

# Create model
model = Network('Jackson-3Q')

source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(lam))
queue1.setService(oclass, Exp(mu[0]))
queue2.setService(oclass, Exp(mu[1]))
queue3.setService(oclass, Exp(mu[2]))

# Set routing matrix
P = model.initRoutingMatrix()
P.set(oclass, oclass, source, queue1, 1.0)
P.set(oclass, oclass, queue1, queue2, p12)
P.set(oclass, oclass, queue1, queue3, p13)
P.set(oclass, oclass, queue1, sink, p1s)
P.set(oclass, oclass, queue2, queue1, p21)
P.set(oclass, oclass, queue2, queue3, p23)
P.set(oclass, oclass, queue2, sink, p2s)
P.set(oclass, oclass, queue3, sink, p3s)
model.link(P)

# Solve with MAM using INAP method
print('=== Jackson Network ===\n')
solver_inap = MAM(model, 'inap')
avg_table_inap = solver_inap.get_avg_table()

print('MAM (method=inap):')
print(avg_table_inap)

# Solve with MVA for comparison
solver_mva = MVA(model)
avg_table_mva = solver_mva.get_avg_table()

print('\nMVA:')
print(avg_table_mva)
