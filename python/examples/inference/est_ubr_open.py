"""Example: UBR estimation on an open network."""
import numpy as np
from line_solver import Network, Delay, Queue, Source, Sink, Exp, HyperExp, OpenClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

model = Network('model')
node = [None] * 4
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.FCFS)
node[2] = Source(model, 'Source')
node[3] = Sink(model, 'Sink')

jobclass = [None]
jobclass[0] = OpenClass(model, 'Class1', 0)

node[0].setService(jobclass[0], HyperExp(0.5, 3.0, 10.0))
node[1].setService(jobclass[0], Exp(float('nan')))
node[2].setArrival(jobclass[0], Exp(0.1))

P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[3], 1.0)
P.set(jobclass[0], jobclass[0], node[2], node[0], 1.0)
model.link(P)

n = 1000
ts = np.arange(1, n + 1, dtype=float)
arvrate_samples = 2 * np.ones(n) - np.random.rand(n) * 0.15
util_samples = np.ones(n) - np.random.rand(n) * 0.05
respt_samples = 0.7 / (1 - util_samples)

estoptions = ParamEstimator.defaultOptions()
estoptions['method'] = 'ubr'
se = ParamEstimator(model, estoptions)

lambda1 = SampledMetric(MetricType.ArvR, ts, arvrate_samples, node[1], jobclass[0])
respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node[1], jobclass[0])
util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

se.addSamples(lambda1)
se.addSamples(respT1)
se.addSamples(util)
se.interpolate()
estVal = se.estimateAt(node[1])
print('Estimated demands:', estVal)

solver = MVA(model)
print(f'SOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
