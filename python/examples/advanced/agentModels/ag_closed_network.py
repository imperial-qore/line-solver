"""
Closed Network with MAM (INAP method)

This example demonstrates MAM with INAP method on a closed queueing network
with processor-sharing (PS) queues using the RCAT/INAP algorithm.

The RCAT (Reversed Compound Agent Theorem) decomposes the network
into interacting stochastic processes and uses fixed-point iteration
to compute equilibrium measures.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
N = 10         # Number of jobs
mu1 = 2.0      # Service rate at queue 1
mu2 = 1.0      # Service rate at queue 2

# Create model: Queue1 <-> Queue2 (closed loop)
model = Network('Closed-2Q')

queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

cclass = ClosedClass(model, 'Class1', N, queue1)
queue1.setService(cclass, Exp(mu1))
queue2.setService(cclass, Exp(mu2))

model.link(Network.serial_routing([queue1, queue2]))

print(f'=== Closed Network (2 PS Queues, {N} jobs) ===\n')

# Solve with MAM using INAP method
solver_inap = MAM(model, 'inap')
avg_table_inap = solver_inap.get_avg_table()

print('MAM (method=inap):')
print(avg_table_inap)

# Solve with MAM using exact method
solver_exact = MAM(model, 'exact')
avg_table_exact = solver_exact.get_avg_table()

print('\nMAM (method=exact):')
print(avg_table_exact)

# Solve with MVA for comparison
solver_mva = MVA(model)
avg_table_mva = solver_mva.get_avg_table()

print('\nMVA:')
print(avg_table_mva)

# Solve with CTMC for exact results (if state space is small)
if N <= 20:
    solver_ctmc = CTMC(model)
    avg_table_ctmc = solver_ctmc.get_avg_table()
    print('\nCTMC (exact):')
    print(avg_table_ctmc)
