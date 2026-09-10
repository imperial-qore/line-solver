"""
Stochastic network calculus: a delay quantile with a certified violation
probability.

Every other solver in LINE answers with a MEAN. The 'snc' family of SolverBA
answers with a TAIL: given a violation probability eps, it returns a delay d for
which P{D > d} <= eps holds, and the guarantee is valid for any work-conserving
scheduling policy at the station. That is the quantity a service-level objective
is written against.

The model is a single M/M/1 station, whose exact tail is known in closed form, so
every number below can be checked. Twin of
matlab/examples/advanced/networkCalculus/snc_delay_quantile.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import math

from line_solver import *
from line_solver.api.snc import (snc_bound_delay, snc_env_poisson, snc_perc_delay,
                                 snc_srv_exp)

print('=== Stochastic network calculus: delay quantile ===\n')

lam = 0.6
mu = 1.0
model = Network('SncQuantile')
source = Source(model, 'Source')
queue = Queue(model, 'Queue', SchedStrategy.FCFS)
sink = Sink(model, 'Sink')
jobclass = OpenClass(model, 'Class1')
source.setArrival(jobclass, Exp(lam))
queue.setService(jobclass, Exp(mu))
model.link(Network.serialRouting(source, queue, sink))

solver = SolverBA(model, 'snc.upper')

# getPercTable is to the snc family what getAvgTable is to a mean solver: one row
# per station and class, reporting the response-time and queue-length quantiles at
# the requested violation probability.
print(solver.getPercTable(1e-3))

# The exact M/M/1 sojourn tail is exp(-(mu-lam)*d) and the exact queue-length tail
# is rho^(n+1). The bound reproduces both DECAY RATES exactly and pays a constant
# prefactor, so the ratio of the bounded quantile to the exact one falls towards 1
# as eps is tightened: the family is at its best exactly where simulation is at
# its worst, deep in the tail.
print('\n%-8s %10s %10s %8s   %10s %10s %8s'
      % ('eps', 'd bound', 'd exact', 'ratio', 'n bound', 'n exact', 'ratio'))
for eps in (1e-2, 1e-3, 1e-6, 1e-9, 1e-12):
    d = solver.getDelayPerc(eps)
    b = solver.getBacklogPerc(eps)
    dexact = -math.log(eps) / (mu - lam)
    nexact = math.log(eps) / math.log(lam / mu) - 1.0
    print('%-8.0e %10.4f %10.4f %8.3f   %10.4f %10.4f %8.3f'
          % (eps, d[1, 0], dexact, d[1, 0] / dexact, b[1, 0], nexact, b[1, 0] / nexact))

# getAvgTable still works: the response time reported by 'snc.upper' is the
# integral of the tail bound, hence an upper bound on the mean. It is loose, and
# deliberately so -- integrating over the whole axis is dominated by the prefactor
# rather than by the decay rate that the family gets right. Use the mean columns to
# bracket, the quantiles to plan.
avg = solver.getAvgTable()
print()
print(avg)
exact_r = 1.0 / (mu - lam)
exact_q = (lam / mu) / (1.0 - lam / mu)
r = float(avg['RespT'].iloc[-1])
q = float(avg['QLen'].iloc[-1])
print('\nexact M/M/1: R = %.4f, Q = %.4f' % (exact_r, exact_q))
print('snc.upper  : R = %.4f (%.1fx), Q = %.4f (%.1fx)'
      % (r, r / exact_r, q, q / exact_q))

# The solver is a thin wrapper over line_solver.api.snc. An arrival envelope and a
# service envelope are callables of the Chernoff parameter theta, and every bound
# is an infimum over theta of a closed-form expression. Note snc_srv_exp, not
# snc_srv_rate: the work unit here is the JOB, so the server is the counting
# process of an Exp(mu) service, and a constant-rate element would model an M/D/1
# and understate the delay.
arv = lambda theta: snc_env_poisson(lam, theta)
srv = lambda theta: snc_srv_exp(mu, theta)
d, theta = snc_perc_delay(arv, srv, 1e-3)
print('\napi: d(1e-3) = %.4f at theta = %.4f' % (d, theta))
print('     the optimal theta approaches log(mu/lam) = %.4f, which is'
      % math.log(mu / lam))
print('     what makes the backlog decay rate exact')
eps_at, theta_at = snc_bound_delay(arv, srv, d)
print('api: P{D > %.4f} <= %.3e  (theta = %.4f), exact tail %.3e'
      % (d, eps_at, theta_at, math.exp(-(mu - lam) * d)))
