"""
G-Network (Gelenbe Network) with Negative Customers

This example demonstrates MAM with RCAT methods on a G-network with negative customers.
Negative customers (signals) remove jobs from queues when they arrive,
modeling job cancellations or service interrupts.

Network topology:
  - Source generates positive customers (Class1) and negative signals (Class2)
  - Positive customers flow: Source -> Queue1 -> Queue2 -> Sink
  - Negative signals target Queue2, removing jobs from it

Reference: Gelenbe, E. (1991). "Product-form queueing networks with
           negative and positive customers", Journal of Applied Probability

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
lambda_pos = 1.0   # Positive customer arrival rate
lambda_neg = 0.3   # Negative signal arrival rate
mu1 = 2.0          # Service rate at Queue1
mu2 = 3.0          # Service rate at Queue2

# Create model
model = Network('GNetwork-Example')

source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

# Positive customer class (normal jobs)
posClass = OpenClass(model, 'Positive')
source.setArrival(posClass, Exp(lambda_pos))
queue1.setService(posClass, Exp(mu1))
queue2.setService(posClass, Exp(mu2))

# Negative signal class (removes jobs from target queue)
negClass = Signal(model, 'Negative', SignalType.NEGATIVE)
source.setArrival(negClass, Exp(lambda_neg))
queue1.setService(negClass, Exp(mu1))  # Signals also get "served" (trigger)
queue2.setService(negClass, Exp(mu2))

# Set routing using P.set(class_src, class_dst, node_src, node_dst, prob) API
P = model.initRoutingMatrix()
# Positive customers: Source -> Queue1 -> Queue2 -> Sink
P.set(posClass, posClass, source, queue1, 1.0)
P.set(posClass, posClass, queue1, queue2, 1.0)
P.set(posClass, posClass, queue2, sink, 1.0)
# Negative signals: Source -> Queue1 -> Queue2 -> Sink
P.set(negClass, negClass, source, queue1, 1.0)
P.set(negClass, negClass, queue1, queue2, 1.0)
P.set(negClass, negClass, queue2, sink, 1.0)
model.link(P)

# Solve with MAM using INAP method
print('=== G-Network with Negative Customers ===\n')
solver_mam = MAM(model, 'inap')
avg_table = solver_mam.get_avg_table()

print('MAM (method=inap):')
print(avg_table)
