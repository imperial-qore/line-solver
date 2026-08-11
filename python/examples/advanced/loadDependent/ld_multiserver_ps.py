"""
Load-Dependent Service - Multi-Server PS Queue

This example demonstrates how to model a multi-server processor-sharing queue
using load-dependent service rates. Compares standard multi-server model with
load-dependent approximation using multiple solvers.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

N = 4  # number of jobs
c = 3  # number of servers

print('=== Load-Dependent Service - Multi-Server PS Queue ===\n')

# Standard multi-server model
print('--- Standard Multi-Server Model ---')
model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

job_class1 = ClosedClass(model, 'Class1', N, delay, 0)
job_class2 = ClosedClass(model, 'Class2', N // 2, delay, 0)

delay.set_service(job_class1, Exp.fit_mean(1.0))
delay.set_service(job_class2, Exp.fit_mean(2.0))

queue1.set_service(job_class1, Exp.fit_mean(1.5))
queue1.set_service(job_class2, Exp.fit_mean(2.5))
queue1.set_number_of_servers(c)

queue2.set_service(job_class1, Exp.fit_mean(3.5))
queue2.set_service(job_class2, Exp.fit_mean(4.5))
queue2.set_number_of_servers(c)

P = model.init_routing_matrix()
P[0][0] = model.serial_routing([delay, queue1, queue2])
P[1][1] = model.serial_routing([delay, queue1, queue2])
model.link(P)

ms_t = MVA(model, 'exact').getAvgTable()
print('MVA (exact):')
print(ms_t)

# Load-dependent model
print('\n--- Load-Dependent Model ---')
ldmodel = Network('ldmodel')

delay = Delay(ldmodel, 'Delay')
queue1 = Queue(ldmodel, 'Queue1', SchedStrategy.PS)
queue2 = Queue(ldmodel, 'Queue2', SchedStrategy.PS)

job_class1 = ClosedClass(ldmodel, 'Class1', N, delay, 0)
job_class2 = ClosedClass(ldmodel, 'Class2', N // 2, delay, 0)

delay.set_service(job_class1, Exp.fit_mean(1.0))
delay.set_service(job_class2, Exp.fit_mean(2.0))

queue1.set_service(job_class1, Exp.fit_mean(1.5))
queue1.set_service(job_class2, Exp.fit_mean(2.5))
# Multi-server with c servers using load dependence
queue1.set_load_dependence([min(i, c) for i in range(1, N + N // 2 + 1)])

queue2.set_service(job_class1, Exp.fit_mean(3.5))
queue2.set_service(job_class2, Exp.fit_mean(4.5))
# Multi-server with c servers using load dependence
queue2.set_load_dependence([min(i, c) for i in range(1, N + N // 2 + 1)])

P = ldmodel.init_routing_matrix()
P[0][0] = ldmodel.serial_routing([delay, queue1, queue2])
P[1][1] = ldmodel.serial_routing([delay, queue1, queue2])
ldmodel.link(P)

print('\nCTMC (exact):')
lld_avg_table_ctmc = CTMC(ldmodel).getAvgTable()
print(lld_avg_table_ctmc)

print('\nNC (exact):')
lld_avg_table_nc = NC(ldmodel, method='exact').getAvgTable()
print(lld_avg_table_nc)

print('\nNC with RD method:')
lld_avg_table_rd = NC(ldmodel, method='rd').getAvgTable()
print(lld_avg_table_rd)

print('\nNC with NRP method:')
lld_avg_table_nrp = NC(ldmodel, method='nrp').getAvgTable()
print(lld_avg_table_nrp)

print('\nNC with NRL method:')
lld_avg_table_nrl = NC(ldmodel, method='nrl').getAvgTable()
print(lld_avg_table_nrl)

print('\nMVA (exact):')
lld_avg_table_mva_ld = MVA(ldmodel, method='exact').getAvgTable()
print(lld_avg_table_mva_ld)

print('\nMVA with QD method:')
lld_avg_table_qd = MVA(ldmodel, method='qd').getAvgTable()
print(lld_avg_table_qd)

print('\nJMT Simulation:')
lld_avg_table_jmt = JMT(ldmodel, seed=23000).getAvgTable()
print(lld_avg_table_jmt)

print('\nNote: Load dependence allows modeling multi-server queues by scaling')
print('      service rates based on the number of jobs present.')
