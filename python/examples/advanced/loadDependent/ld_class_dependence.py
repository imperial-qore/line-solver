"""
Load-Dependent Service - Class Dependence

This example demonstrates class-dependent service rates, where the service
rate depends on the per-class population at a station. In this example,
multi-server behavior is applied only for class 1 jobs.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

N = 16  # number of jobs in class 1
c = 2   # number of servers

print('=== Load-Dependent Service - Class Dependence ===\n')

# Class-dependent model
cdmodel = Network('model')

delay = Delay(cdmodel, 'Delay')
queue = Queue(cdmodel, 'Queue1', SchedStrategy.PS)

job_class1 = ClosedClass(cdmodel, 'Class1', N, delay, 0)
job_class2 = ClosedClass(cdmodel, 'Class2', N // 2, delay, 0)

delay.set_service(job_class1, Exp.fit_mean(1.0))
delay.set_service(job_class2, Exp.fit_mean(2.0))

queue.set_service(job_class1, Exp.fit_mean(1.5))
queue.set_service(job_class2, Exp.fit_mean(2.5))

# Class dependence: multi-server only for class 1 jobs
# ni is a vector where ni[r] is the number of jobs in class r at the station
# Service rate scales with class 1 population (ni[0]), up to c servers
queue.set_class_dependence(lambda ni: min(ni[0], c))

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

print('\nJMT Simulation:')
cd_avg_table_jmt = JMT(cdmodel, seed=23000).getAvgTable()
print(cd_avg_table_jmt)

print('\nNote: Class-dependent service allows different scaling for different classes.')
print('      In this example, service rate scales with Class 1 population only,')
print('      modeling c servers available exclusively for Class 1 jobs.')
print('      Class 2 jobs do not benefit from multi-server parallelism.')
