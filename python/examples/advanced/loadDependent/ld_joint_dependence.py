"""
Load-Dependent Service - Joint Dependence

This example demonstrates joint-dependent (NON-product-form) service rates: the
service rate reads a foreign class marginal (class-1 population), so it is
not expressible as a product-form beta_{i,r}(n_{i,r}). In this example,
multi-server behavior is applied only for class 1 jobs.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

N = 16  # number of jobs in class 1
c = 2   # number of servers

print('=== Load-Dependent Service - Joint Dependence ===\n')

# Joint-dependent (non-product-form) model
jdmodel = Network('model')

delay = Delay(jdmodel, 'Delay')
queue = Queue(jdmodel, 'Queue1', SchedStrategy.PS)

job_class1 = ClosedClass(jdmodel, 'Class1', N, delay, 0)
job_class2 = ClosedClass(jdmodel, 'Class2', N // 2, delay, 0)

delay.set_service(job_class1, Exp.fit_mean(1.0))
delay.set_service(job_class2, Exp.fit_mean(2.0))

queue.set_service(job_class1, Exp.fit_mean(1.5))
queue.set_service(job_class2, Exp.fit_mean(2.5))

# Joint dependence (non-product-form): rate reads the class-1 marginal only,
# shared across classes. ni is a vector where ni[r] is the class-r population.
queue.set_joint_dependence(lambda ni: min(ni[0], c), c)  # peak rate scaling = c (Util = T*S/c)

P = jdmodel.init_routing_matrix()
P[0][0] = jdmodel.serial_routing([delay, queue])
P[1][1] = jdmodel.serial_routing([delay, queue])
jdmodel.link(P)

print('CTMC (exact):')
jd_avg_table_ctmc = CTMC(jdmodel).getAvgTable()
print(jd_avg_table_ctmc)

print('\nMVA with QD method:')
jd_avg_table_jd = MVA(jdmodel, method='qd').getAvgTable()
print(jd_avg_table_jd)

# JMT is not solved here: the JSIM writer has no representation for the
# joint-dependence handle, so SolverJMT rejects the model rather than
# silently solving it unscaled (see SolverJMT.getFeatureSet).

print('\nNote: Joint-dependent service is non-product-form (see set_joint_dependence).')
print('      In this example, service rate scales with Class 1 population only,')
print('      modeling c servers available exclusively for Class 1 jobs.')
print('      Class 2 jobs do not benefit from multi-server parallelism.')
