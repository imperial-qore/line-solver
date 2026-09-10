"""Example: EKF estimation with autoMethod selection."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

np.random.seed(1)
model = Network('model')
node = [None] * 2
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
jobclass = [None]
jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
node[0].setService(jobclass[0], Exp.fitMean(1.0))
node[1].setService(jobclass[0], Exp(float('nan')))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
model.link(P)

n = 100
ts = np.arange(1, n + 1, dtype=float)
arvr_samples = 1.5 * np.ones(n) - np.random.rand(n) * 0.1
util_samples = 0.4 * arvr_samples
respt_samples = 0.4 / (1 - util_samples)

lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node[1], jobclass[0])
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node[1], jobclass[0])
util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

# EKF
print('\n=== EKF Estimator ===')
options = ParamEstimator.defaultOptions()
options['method'] = 'ekf'
se = ParamEstimator(model, options)
se.addSamples(lambda1)
se.addSamples(respT1)
se.addSamples(util)
se.interpolate()
estVal_ekf = se.estimateAt(node[1])
print(f'EKF demand: Class1={estVal_ekf[0]:.4f} (true=0.4000)')

# autoMethod
print('\n=== autoMethod selection ===')
node[1].setService(jobclass[0], Exp(float('nan')))
model.reset()
options2 = ParamEstimator.defaultOptions()
options2['method'] = 'auto'
se2 = ParamEstimator(model, options2)
se2.addSamples(lambda1)
se2.addSamples(respT1)
se2.addSamples(util)
se2.interpolate()
method = se2.autoMethod()
print(f'autoMethod selected: {method}')
print(f'Required metrics for {method}: {ParamEstimator.getRequiredMetrics(method)}')
estVal_auto = se2.estimateAt(node[1])
print(f'Auto demand: Class1={estVal_auto[0]:.4f} (true=0.4000)')

solver = MVA(model)
print(f'\nSOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
