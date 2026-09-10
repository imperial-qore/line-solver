"""
State Probabilities - Joint per-station TOTAL queue lengths

Joint probability of the per-station TOTAL queue lengths, all classes summed
out, from SolverNC.getProbSysMarg.

Compare with statepr_sys_aggr.py, which fixes the PER-CLASS population of every
station: each probability here is the sum of that one over every per-class
table with these row sums. The fibre grows combinatorially, so the quantity is
evaluated as a permanent of the demand matrix replicated once per job rather
than by enumerating it.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
from line_solver import *
from line_solver.api.pfqn import multichoose

model = Network('model')

delay = Delay(model, 'Delay')
queue1 = Queue(model, 'Queue1', SchedStrategy.PS)
queue2 = Queue(model, 'Queue2', SchedStrategy.PS)

N = [2, 1]
job_class1 = ClosedClass(model, 'Class1', N[0], delay, 0)
job_class2 = ClosedClass(model, 'Class2', N[1], delay, 0)

delay.setService(job_class1, Exp(1 / 1.5))
delay.setService(job_class2, Exp(1 / 2.0))
queue1.setService(job_class1, Exp(1 / 0.7))
queue1.setService(job_class2, Exp(1 / 0.4))
queue2.setService(job_class1, Exp(1 / 0.3))
queue2.setService(job_class2, Exp(1 / 0.9))

P = model.initRoutingMatrix()
P.set(job_class1, job_class1, Network.serialRouting(delay, queue1, queue2))
P.set(job_class2, job_class2, Network.serialRouting(delay, queue1, queue2))
model.link(P)

solver = SolverNC(model)

# Every way of splitting the closed population across the stations
states = np.atleast_2d(multichoose(model.getNumberOfStations(), sum(N)))

print('\n  n(Delay) n(Queue1) n(Queue2)        P(n)')
p = np.zeros(states.shape[0])
for j in range(states.shape[0]):
    p[j], _ = solver.getProbSysMarg(states[j, :])
    print('  %8d %9d %9d  %10.6f' % (states[j, 0], states[j, 1], states[j, 2], p[j]))
print('  ' + '-' * 42)
print('  sum %39.6f' % np.sum(p))

# The law is exact, so its first moments are the queue lengths
print('\n  E[n] from the joint law : %s' % np.array2string(p @ states, precision=8))
qlen = np.atleast_2d(np.asarray(SolverCTMC(model).getAvgQLen(), dtype=float))
print('  QLen from SolverCTMC    : %s' % np.array2string(np.sum(qlen, axis=1), precision=8))

# The approximate permanent engines trade accuracy for cost on models whose
# class count makes the exact expansion dear. They need a demand matrix with
# full support and refuse a structural zero rather than flooring it.
pb = np.zeros(states.shape[0])
for j in range(states.shape[0]):
    pb[j], _ = solver.getProbSysMarg(states[j, :], 'bethe')
pb = pb / np.sum(pb)
print('  Bethe engine, mean relative error : %.2f%%' % (100 * np.mean(np.abs(pb - p) / p)))

print(solver.citations())
