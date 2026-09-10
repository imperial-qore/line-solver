"""Example: QMLE estimation on a closed network."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

np.random.seed(1)
model = Network('model')
node = [None] * 2
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
jobclass = [None] * 2
jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
jobclass[1] = ClosedClass(model, 'Class2', 3, node[0], 0)
node[0].setService(jobclass[0], Exp.fitMean(1.0))
node[0].setService(jobclass[1], Exp.fitMean(1.0))
node[1].setService(jobclass[0], Exp(float('nan')))
node[1].setService(jobclass[1], Exp(float('nan')))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
P.set(jobclass[1], jobclass[1], node[1], node[0], 1.0)
model.link(P)

n = 1000
ts = np.arange(1, n + 1, dtype=float)
qlen1_samples = 0.3 * np.ones(n) + np.random.rand(n) * 0.1
qlen2_samples = 0.8 * np.ones(n) + np.random.rand(n) * 0.1

options = ParamEstimator.defaultOptions()
options['method'] = 'qmle'
se = ParamEstimator(model, options)
ql1 = SampledMetric(MetricType.QLen, ts, qlen1_samples, node[1], jobclass[0])
ql2 = SampledMetric(MetricType.QLen, ts, qlen2_samples, node[1], jobclass[1])
se.addSamples(ql1)
se.addSamples(ql2)
se.interpolate()
estVal = se.estimateAt(node[1])
print('Estimated demands:', estVal)

solver = MVA(model)
print(f'SOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
