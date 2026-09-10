"""
Multiclass FES Aggregation with Norton's Theorem Verification

This example demonstrates Norton's theorem for closed multiclass
queueing networks: a subset of stations is replaced by a single
Flow-Equivalent Server (FES) whose class-dependent service rates equal the
throughputs of the isolated subnetwork.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import *
from line_solver.io.model_adapter import ModelAdapter

N1 = 3  # class-1 jobs
N2 = 2  # class-2 jobs

# Original 4-station tandem network with 2 classes
model = Network('OriginalModel')

delay = Delay(model, 'Delay')
q1 = Queue(model, 'Q1', SchedStrategy.PS)
q2 = Queue(model, 'Q2', SchedStrategy.PS)
q3 = Queue(model, 'Q3', SchedStrategy.PS)

jobclass1 = ClosedClass(model, 'Class1', N1, delay, 0)
jobclass2 = ClosedClass(model, 'Class2', N2, delay, 0)

delay.set_service(jobclass1, Exp.fit_mean(1.0))
delay.set_service(jobclass2, Exp.fit_mean(1.5))
q1.set_service(jobclass1, Exp.fit_mean(0.5))
q1.set_service(jobclass2, Exp.fit_mean(0.8))
q2.set_service(jobclass1, Exp.fit_mean(0.3))
q2.set_service(jobclass2, Exp.fit_mean(0.6))
q3.set_service(jobclass1, Exp.fit_mean(0.4))
q3.set_service(jobclass2, Exp.fit_mean(0.7))

P = model.init_routing_matrix()
P[0][0] = model.serial_routing([delay, q1, q2, q3])
P[1][1] = model.serial_routing([delay, q1, q2, q3])
model.link(P)

# Solve original model
print('MVA (original):')
avg_orig = MVA(model, method='exact').getAvgTable()
print(avg_orig)

# Aggregate Q1, Q2, Q3 into FES and solve with NC (convolution)
result = ModelAdapter.aggregate_fes(model, [q1, q2, q3])
fes_model = result['fes_model']
print('FES model created with %d stations' % fes_model.get_number_of_stations())

print('NC (FES model):')
avg_nc = NC(fes_model, method='exact').getAvgTable()
print(avg_nc)
