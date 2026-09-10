"""
Busy Period of a Subnetwork

The busy period of order n for a set of stations is the time from the instant a
job entering the set finds n-1 jobs in it up to the next instant when fewer than
n remain. SolverNC evaluates it exactly from the normalizing constants of the
subnetwork and of its complement (H. Daduna, "Busy Periods for Subnetworks in
Stochastic Networks: Mean Value Analysis", J. ACM 35(3), 1988); SolverLDES
measures the same quantity along a simulated sample path.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from line_solver import *

print('=== Busy Period of a Subnetwork ===\n')

N = 5
model = Network('busyPeriodModel')

queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS)
queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS)
queue3 = Queue(model, 'Queue3', SchedStrategy.FCFS)

jobclass = ClosedClass(model, 'Class1', N, queue1, 0)

queue1.setService(jobclass, Exp(1.5))
queue2.setService(jobclass, Exp(0.9))
queue3.setService(jobclass, Exp(2.0))

P = model.initRoutingMatrix()
P.set(jobclass, jobclass, queue1, queue2, 0.6)
P.set(jobclass, jobclass, queue1, queue3, 0.4)
P.set(jobclass, jobclass, queue2, queue1, 0.7)
P.set(jobclass, jobclass, queue2, queue3, 0.3)
P.set(jobclass, jobclass, queue3, queue1, 0.5)
P.set(jobclass, jobclass, queue3, queue2, 0.5)
model.link(P)

subnets = [[0], [1], [2], [0, 1]]
orders = [1, 3, 5]

# Exact mean value analysis of the busy period
nc = SolverNC(model)
print('Mean busy period of order n (exact, SolverNC)')
print(f"{'subnetwork':<18}{'n=1':>10}{'n=3':>10}{'n=5':>10}")
for subnet in subnets:
    b = nc.getAvgBusyPeriod(subnet, orders)[0]
    print(f'{str(subnet):<18}{b[0]:>10.4f}{b[1]:>10.4f}{b[2]:>10.4f}')

# The same quantity measured by simulation
opts = LDESOptions()
opts.samples = 2000000
opts.seed = 23000
ldes = SolverLDES(model, opts)
print('\nMean busy period of order n (measured, SolverLDES)')
print(f"{'subnetwork':<18}{'n=1':>10}{'n=3':>10}{'n=5':>10}")
for subnet in subnets:
    b, _ = ldes.getAvgBusyPeriod(subnet, -1, orders)
    print(f'{str(subnet):<18}{b[0]:>10.4f}{b[1]:>10.4f}{b[2]:>10.4f}')
