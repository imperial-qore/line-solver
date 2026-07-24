"""
Load-Dependent Service - Class Dependence (product-form)

This example demonstrates class-dependent (product-form) service rates. Each
class sees a scaling beta_{i,r} that depends ONLY on its own per-class
population n_{i,r} at the station. This is the BCMP product-form case of
QD-AMVA (Casale, Perez, Wang, IFIP PERFORMANCE 2015):
D_{i,r}(n) = theta_{i,r} * beta_{i,r}(n_{i,r}). The handle returns a length-R
vector [beta_{i,1}(n), beta_{i,2}(n)], of which the solver picks class r.
Contrast ld_joint_dependence.py, where the scalar min(ni[0], c) reads a foreign
class marginal and is therefore non-product-form (set_joint_dependence).

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

N = 16  # number of jobs in class 1
c = 2   # number of servers

print('=== Load-Dependent Service - Class Dependence (product-form) ===\n')

cdmodel = Network('model')

delay = Delay(cdmodel, 'Delay')
queue = Queue(cdmodel, 'Queue1', SchedStrategy.PS)

job_class1 = ClosedClass(cdmodel, 'Class1', N, delay, 0)
job_class2 = ClosedClass(cdmodel, 'Class2', N // 2, delay, 0)

delay.set_service(job_class1, Exp.fit_mean(1.0))
delay.set_service(job_class2, Exp.fit_mean(2.0))

queue.set_service(job_class1, Exp.fit_mean(1.5))
queue.set_service(job_class2, Exp.fit_mean(2.5))

# beta_{i,r}(n_{i,r}): class 1 scales up to c servers with its OWN count, class
# 2 is always single-server. Both entries read only their own marginal, so the
# demands satisfy the product-form recurrence. Peak rate scaling [c, 1] per
# class normalizes Util = T*S/peak.
import numpy as _np
queue.set_class_dependence(lambda ni: [min(_np.ravel(ni)[0], c), 1], [c, 1])

P = cdmodel.init_routing_matrix()
P[0][0] = cdmodel.serial_routing([delay, queue])
P[1][1] = cdmodel.serial_routing([delay, queue])
cdmodel.link(P)

print('CTMC (exact):')
cd_avg_table_ctmc = CTMC(cdmodel).getAvgTable()
print(cd_avg_table_ctmc)

print('\nMVA with QD method:')
cd_avg_table_cd = MVA(cdmodel, method='qd').getAvgTable()
print(cd_avg_table_cd)

# JMT is not solved here: the JSIM writer has no representation for the
# class-dependence handle, so SolverJMT rejects the model rather than
# silently solving it unscaled (see SolverJMT.getFeatureSet).

print('\nNote: Class-dependent (product-form) service scales each class by its')
print('      own per-class population, preserving the BCMP product form.')
