"""
Load-Dependent Service - Multi-Server FCFS Queue

This example demonstrates load-dependent service for modeling multi-server
FCFS queues, comparing three approaches: standard multi-server, load-dependent,
and class-dependent service rates.

Copyright (c) 2012-2025, Imperial College London
All rights reserved.
"""

from line_solver import *

N = 16  # number of jobs
c = 2   # number of servers

print('=== Load-Dependent Service - Multi-Server FCFS Queue ===\n')

# Standard multi-server model
print('--- Standard Multi-Server Model ---')
model = Network('model')

delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue1', SchedStrategy.FCFS)

job_class = ClosedClass(model, 'Class1', N, delay, 0)

delay.set_service(job_class, Exp.fit_mean(1.0))
queue.set_service(job_class, Exp.fit_mean(1.5))
queue.set_number_of_servers(c)

model.link(model.serial_routing([delay, queue]))

ms_t = NC(model).getAvgTable()
print('NC Solver:')
print(ms_t)

# Load-dependent model
print('\n--- Load-Dependent Model ---')
ldmodel = Network('model')

delay = Delay(ldmodel, 'Delay')
queue = Queue(ldmodel, 'Queue1', SchedStrategy.FCFS)

job_class = ClosedClass(ldmodel, 'Class1', N, delay, 0)

delay.set_service(job_class, Exp.fit_mean(1.0))
queue.set_service(job_class, Exp.fit_mean(1.5))
# Multi-server with c servers using load dependence
queue.set_load_dependence([min(i, c) for i in range(1, N + 1)])

ldmodel.link(ldmodel.serial_routing([delay, queue]))

print('\nCTMC (exact):')
lld_avg_table_ctmc = CTMC(ldmodel).getAvgTable()
print(lld_avg_table_ctmc)

print('\nNC (exact):')
lld_avg_table_nc = NC(ldmodel).getAvgTable()
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

# Class-dependent model
print('\n--- Class-Dependent Model ---')
cdmodel = Network('model')

delay = Delay(cdmodel, 'Delay')
queue = Queue(cdmodel, 'Queue1', SchedStrategy.FCFS)

job_class = ClosedClass(cdmodel, 'Class1', N, delay, 0)

delay.set_service(job_class, Exp.fit_mean(1.0))
queue.set_service(job_class, Exp.fit_mean(1.5))
# Class dependence: service rate depends on per-class populations
# ni is a vector where ni[r] is the number of jobs in class r at station
queue.set_class_dependence(lambda ni: min(sum(ni), c))

cdmodel.link(cdmodel.serial_routing([delay, queue]))

print('\nCTMC (exact):')
cd_avg_table_ctmc = CTMC(cdmodel).getAvgTable()
print(cd_avg_table_ctmc)

print('\nMVA with QD method:')
cd_avg_table_cd = MVA(cdmodel, method='qd').getAvgTable()
print(cd_avg_table_cd)

print('\nJMT Simulation:')
cd_avg_table_jmt = JMT(cdmodel, seed=23000).getAvgTable()
print(cd_avg_table_jmt)
