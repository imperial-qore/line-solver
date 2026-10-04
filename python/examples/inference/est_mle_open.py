"""Example: MLE estimation on an open network."""
import numpy as np
from line_solver import Network, Delay, Queue, Source, Sink, Exp, OpenClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

model = Network('model')
node = [None] * 4
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
node[2] = Source(model, 'Source')
node[3] = Sink(model, 'Sink')
jobclass = [None] * 2
jobclass[0] = OpenClass(model, 'Class1', 0)
jobclass[1] = OpenClass(model, 'Class2', 0)
node[0].setService(jobclass[0], Exp.fitMean(1.0))
node[0].setService(jobclass[1], Exp.fitMean(1.0))
node[1].setService(jobclass[0], Exp(float('nan')))
node[1].setService(jobclass[1], Exp(float('nan')))
node[2].setArrival(jobclass[0], Exp(1.0))
node[2].setArrival(jobclass[1], Exp(0.5))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[3], 1.0)
P.set(jobclass[0], jobclass[0], node[2], node[0], 1.0)
P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
P.set(jobclass[1], jobclass[1], node[1], node[3], 1.0)
P.set(jobclass[1], jobclass[1], node[2], node[0], 1.0)
model.link(P)

n = 100
ts = np.arange(1, n + 1, dtype=float)
arvr1_samples = 1.0 * np.ones(n) + np.random.randn(n) * 0.05
arvr2_samples = 0.5 * np.ones(n) + np.random.randn(n) * 0.03
util_samples = 0.25 * np.ones(n) + np.random.randn(n) * 0.02
util_samples = np.clip(util_samples, 0.01, 0.95)
respt1_samples = 0.1 / (1 - util_samples)
respt2_samples = 0.3 / (1 - util_samples)

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node[1], jobclass[0])
lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node[1], jobclass[1])
respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node[1], jobclass[0])
respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node[1], jobclass[1])
util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

print('\n=== MLE Estimator (Open Network) ===')
options = ParamEstimator.defaultOptions()
options['method'] = 'mle'
se = ParamEstimator(model, options)
se.addSamples(lambda1)
se.addSamples(lambda2)
se.addSamples(respT1)
se.addSamples(respT2)
se.addSamples(util)
se.interpolate()
estVal = se.estimateAt(node[1])
print(f'MLE demands: Class1={estVal[0]:.4f} (true=0.1000), Class2={estVal[1]:.4f} (true=0.3000)')

solver = MVA(model)
print(f'\nSOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
