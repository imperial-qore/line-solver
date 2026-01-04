"""
G-Network (Gelenbe Network) with Negative Customers

This example demonstrates MAM with INAP method on a G-network with negative customers.
Negative customers (signals) remove jobs from queues when they arrive,
modeling job cancellations or service interrupts.

Network topology:
  - Source generates positive customers (Positive) and negative signals (Negative)
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
pos_class = OpenClass(model, 'Positive')
source.setArrival(pos_class, Exp(lambda_pos))
queue1.setService(pos_class, Exp(mu1))
queue2.setService(pos_class, Exp(mu2))

# Negative signal class (removes jobs from target queue)
# Using Signal class with SignalType.NEGATIVE for automatic G-network handling
neg_class = Signal(model, 'Negative', SignalType.NEGATIVE)
source.setArrival(neg_class, Exp(lambda_neg))
queue1.setService(neg_class, Exp(mu1))  # Signals also get "served" (trigger)
queue2.setService(neg_class, Exp(mu2))

# Set routing using RoutingMatrix
P = model.init_routing_matrix()
# Positive customers: Source -> Queue1 -> Queue2 -> Sink
P.set(pos_class, pos_class, source, queue1, 1.0)
P.set(pos_class, pos_class, queue1, queue2, 1.0)
P.set(pos_class, pos_class, queue2, sink, 1.0)
# Negative signals: Source -> Queue1 -> Queue2 (removes job) -> Sink
P.set(neg_class, neg_class, source, queue1, 1.0)
P.set(neg_class, neg_class, queue1, queue2, 1.0)
P.set(neg_class, neg_class, queue2, sink, 1.0)
model.link(P)

# Solve with MAM using INAP method
solver_mam = MAM(model, 'inap')
avg_table_mam = solver_mam.get_avg_table()
print(avg_table_mam)
