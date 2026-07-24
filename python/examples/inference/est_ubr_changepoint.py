"""Example: UBR sliding-window change-point detection (open network)."""
import numpy as np
from line_solver import Network, Delay, Queue, Source, Sink, Exp, OpenClass, SchedStrategy, MetricType
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

n = 200
changeT = 100
D1 = 0.3
D2 = 0.6
lam = 1.0
noise = 0.02

ts = np.arange(1, n + 1, dtype=float)
arvr_samples = lam * np.ones(n) + np.random.randn(n) * noise

util_samples = np.zeros(n)
respt_samples = np.zeros(n)
for t in range(n):
    D = D1 if t < changeT else D2
    U = lam * D + np.random.randn() * noise
    U = max(0.01, min(0.95, U))
    util_samples[t] = U
    respt_samples[t] = D / (1 - U) + np.random.randn() * noise * 0.1

W = 30
estimates = np.zeros(n - W + 1)

for t in range(W - 1, n):
    idx = slice(t - W + 1, t + 1)
    node[1].setService(jobclass[0], Exp(float('nan')))
    model.reset()

    estoptions = ParamEstimator.defaultOptions()
    estoptions['method'] = 'ubr'
    se = ParamEstimator(model, estoptions)

    se.addSamples(SampledMetric(MetricType.ArvR, ts[idx], arvr_samples[idx], node[1], jobclass[0]))
    se.addSamples(SampledMetric(MetricType.Util, ts[idx], util_samples[idx], node[1]))
    se.interpolate()
    estVal = se.estimateAt(node[1])
    estimates[t - W + 1] = estVal[0] if hasattr(estVal, '__len__') else estVal

print('\n=== UBR Change-Point Detection ===')
print(f'True demand: D={D1} (t<=100), D={D2} (t>100)')
print(f'Window size: {W}\n')

checkpoints = [W, 50, 80, 100, 110, 120, 130, 150, 180, 200]
print(f'{"t":>6}  {"Estimated":>10}  {"True":>10}')
print(f'{"------":>6}  {"----------":>10}  {"----------":>10}')
for cp in checkpoints:
    if W <= cp <= n:
        est = estimates[cp - W]
        trueD = D1 if cp <= changeT else D2
        print(f'{cp:6d}  {est:10.4f}  {trueD:10.4f}')
