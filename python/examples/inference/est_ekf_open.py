"""Example: EKF estimation on an open network."""
import numpy as np
from line_solver import Network, Delay, Queue, Source, Sink, Exp, OpenClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

model = Network('model')
node = [None] * 4
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
node[2] = Source(model, 'Source')
node[3] = Sink(model, 'Sink')
jobclass = [None]
jobclass[0] = OpenClass(model, 'Class1', 0)
node[0].setService(jobclass[0], Exp.fitMean(1.0))
node[1].setService(jobclass[0], Exp(float('nan')))
node[2].setArrival(jobclass[0], Exp(1.0))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[3], 1.0)
P.set(jobclass[0], jobclass[0], node[2], node[0], 1.0)
model.link(P)

n = 100
ts = np.arange(1, n + 1, dtype=float)
arvr_samples = 1.0 * np.ones(n) + np.random.randn(n) * 0.05
util_samples = 0.4 * np.ones(n) + np.random.randn(n) * 0.02
util_samples = np.clip(util_samples, 0.01, 0.95)
respt_samples = 0.4 / (1 - util_samples)

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node[1], jobclass[0])
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node[1], jobclass[0])
util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

print('\n=== EKF Estimator (Open Network) ===')
options = ParamEstimator.defaultOptions()
options['method'] = 'ekf'
se = ParamEstimator(model, options)
se.addSamples(lambda1)
se.addSamples(respT1)
se.addSamples(util)
se.interpolate()
estVal = se.estimateAt(node[1])
print(f'EKF demand: Class1={estVal[0]:.4f} (true=0.4000)')

solver = MVA(model)
print(f'\nSOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
