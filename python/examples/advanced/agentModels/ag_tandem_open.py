"""
Open Tandem Queue with MAM (INAP method)

This example demonstrates MAM with INAP method on an open tandem queueing network
(M/M/1 -> M/M/1) using the RCAT/INAP algorithm.

Reference: Marin and Rota-Bulo', "A Mean-Field Analysis of a Class of
           Interactive Distributed Systems", MASCOTS 2009

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
arrival_rate = 0.5   # Arrival rate
service_rate1 = 1.0  # Service rate at queue 1
service_rate2 = 1.5  # Service rate at queue 2

# Create model: Source -> Queue1 -> Queue2 -> Sink
model = Network('Tandem-MM1')

source = Source(model, 'Source')
queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')

oclass = OpenClass(model, 'Class1')
source.setArrival(oclass, Exp(arrival_rate))
queue1.setService(oclass, Exp(service_rate1))
queue2.setService(oclass, Exp(service_rate2))

model.link(Network.serial_routing([source, queue1, queue2, sink]))

# Analytical M/M/1 solution
U1_exact = arrival_rate / service_rate1
U2_exact = arrival_rate / service_rate2
Q1_exact = U1_exact / (1 - U1_exact)
Q2_exact = U2_exact / (1 - U2_exact)
R1_exact = 1 / (service_rate1 - arrival_rate)
R2_exact = 1 / (service_rate2 - arrival_rate)

print('=== Open Tandem Queue (M/M/1 -> M/M/1) ===\n')
print('Analytical (M/M/1):')
print(f'  Queue1: U={U1_exact:.4f}, Q={Q1_exact:.4f}, R={R1_exact:.4f}')
print(f'  Queue2: U={U2_exact:.4f}, Q={Q2_exact:.4f}, R={R2_exact:.4f}\n')

# Solve with MAM using INAP method
solver_inap = MAM(model, 'inap')
avg_table_inap = solver_inap.get_avg_table()

print('MAM (method=inap):')
print(avg_table_inap)

# Note: 'exact' method (AutoCAT) is not yet implemented in native Python
# It falls back to INAP with a warning in the JAR version

# Solve with MVA for comparison
solver_mva = MVA(model)
avg_table_mva = solver_mva.get_avg_table()

print('\nMVA:')
print(avg_table_mva)
