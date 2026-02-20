"""
Multiclass Closed Network with MAM

This example demonstrates MAM with RCAT methods on a multiclass closed queueing
network. The RCAT algorithm handles multiple job classes by creating
separate processes for each (station, class) pair.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

# Parameters
N1 = 5       # Number of class 1 jobs
N2 = 3       # Number of class 2 jobs
mu1_q1 = 2.0  # Service rate for class 1 at queue 1
mu2_q1 = 1.5  # Service rate for class 2 at queue 1
mu1_q2 = 1.0  # Service rate for class 1 at queue 2
mu2_q2 = 0.8  # Service rate for class 2 at queue 2

# Create model
model = Network('Multiclass-Closed')

queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

class1 = ClosedClass(model, 'Class1', N1, queue1)
class2 = ClosedClass(model, 'Class2', N2, queue1)

queue1.setService(class1, Exp(mu1_q1))
queue1.setService(class2, Exp(mu2_q1))
queue2.setService(class1, Exp(mu1_q2))
queue2.setService(class2, Exp(mu2_q2))

# Simple routing: each class cycles through both queues
P = model.initRoutingMatrix()
P.set(class1, class1, queue1, queue2, 1.0)
P.set(class1, class1, queue2, queue1, 1.0)
P.set(class2, class2, queue1, queue2, 1.0)
P.set(class2, class2, queue2, queue1, 1.0)
model.link(P)

# Solve with MAM using INAP method
print('=== Multiclass Closed Network ===\n')
solver_inap = MAM(model, 'inap')
avg_table_inap = solver_inap.get_avg_table()

print('MAM (method=inap):')
print(avg_table_inap)

# Solve with MVA for comparison
solver_mva = MVA(model)
avg_table_mva = solver_mva.get_avg_table()

print('\nMVA:')
print(avg_table_mva)
