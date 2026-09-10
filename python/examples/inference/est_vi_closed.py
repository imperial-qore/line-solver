"""Variational inference for Markovian queueing networks, on a closed loop.

Method 'vi' (I. Perez, G. Casale, Adv. Appl. Prob. 53(3), 2021) infers service
rates from NOISY QUEUE-LENGTH READINGS taken over time: each reading is exact
with probability 1-epsilon and uniform over the remaining feasible values
otherwise. Unlike the other estimators it returns a conjugate Gamma POSTERIOR
per rate, not only a point estimate.

It reads the queue lengths of EVERY station, not only the estimated one: the
transition counts the method is written in are pinned by the whole picture.
"""
import numpy as np

from line_solver import (Network, Delay, Queue, Exp, ClosedClass, SchedStrategy,
                         MetricType, SolverMVA)
from line_solver.inference import ParamEstimator, SampledMetric

N = 10

# define model, with the queue rate to be estimated
model = Network('model')
delay = Delay(model, 'Delay')
queue = Queue(model, 'Queue1', SchedStrategy.FCFS)
cl = ClosedClass(model, 'Class1', N, delay, 0)
delay.setService(cl, Exp(0.5))
queue.setService(cl, Exp(2.0))     # starting point of the estimate
model.link(Network.serialRouting(delay, queue))

# queue-length readings, one per unit time, 10% of them faulty
ts = np.arange(1, 11, dtype=float)
qlen = np.array([1, 2, 3, 2, 4, 3, 5, 4, 3, 4], dtype=float)

# estimate the service rate of the queue
se = ParamEstimator(model)
se.options['method'] = 'vi'
se.options['epsilon'] = 0.1       # probability that a reading is faulty
se.options['prior_shape'] = 2.0   # Gamma prior shape; the rate is set from the model
se.options['ngrid'] = 51          # time grid of the backward and forward passes
se.options['nsamples'] = 32       # lattice points per marginal
se.options['ymax'] = 60           # transition-count truncation
se.options['iter_max'] = 5
se.addSamples(SampledMetric(MetricType.QLen, ts, N - qlen, delay, cl))
se.addSamples(SampledMetric(MetricType.QLen, ts, qlen, queue, cl))
est = np.ravel(se.estimate_at([queue]))
post = np.atleast_2d(se.options['posterior'])

print('Estimated demand: %.8f' % est[0])
print('posterior service rate ~ Gamma(%.4f, %.4f), mean %.4f'
      % (post[0, 0], post[0, 1], post[0, 0] / post[0, 1]))
print('evidence lower bound over the iterations:',
      np.round(np.asarray(se.options['bound']), 3))

# solve the model the estimate has been written into
print(SolverMVA(model).getAvgTable())
