"""Example: UBO sliding-window change-point detection (open network, 2 classes)."""
import numpy as np
from line_solver import Network, Delay, Queue, Source, Sink, Exp, OpenClass, SchedStrategy, MetricType
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

n = 200; changeT = 100
D1_before = 0.1; D1_after = 0.3; D2_val = 0.2
lambda1 = 1.0; lambda2 = 0.5; noise = 0.01
ts = np.arange(1, n + 1, dtype=float)
arvr1_samples = lambda1 * np.ones(n) + np.random.randn(n) * noise
arvr2_samples = lambda2 * np.ones(n) + np.random.randn(n) * noise
util_samples = np.zeros(n)
respt1_samples = np.zeros(n)
respt2_samples = np.zeros(n)
for t in range(n):
    d1 = D1_before if t < changeT else D1_after
    U = lambda1 * d1 + lambda2 * D2_val + np.random.randn() * noise
    U = max(0.01, min(0.95, U))
    util_samples[t] = U
    respt1_samples[t] = d1 / (1 - U) + np.random.randn() * noise * 0.1
    respt2_samples[t] = D2_val / (1 - U) + np.random.randn() * noise * 0.1

W = 30
est1 = np.zeros(n - W + 1)
est2 = np.zeros(n - W + 1)
for t in range(W - 1, n):
    idx = slice(t - W + 1, t + 1)
    node[1].setService(jobclass[0], Exp(float('nan')))
    node[1].setService(jobclass[1], Exp(float('nan')))
    model.reset()
    estoptions = ParamEstimator.defaultOptions()
    estoptions['method'] = 'ubo'
    se = ParamEstimator(model, estoptions)
    se.addSamples(SampledMetric(MetricType.ArvR, ts[idx], arvr1_samples[idx], node[1], jobclass[0]))
    se.addSamples(SampledMetric(MetricType.ArvR, ts[idx], arvr2_samples[idx], node[1], jobclass[1]))
    se.addSamples(SampledMetric(MetricType.RespT, ts[idx], respt1_samples[idx], node[1], jobclass[0]))
    se.addSamples(SampledMetric(MetricType.RespT, ts[idx], respt2_samples[idx], node[1], jobclass[1]))
    se.addSamples(SampledMetric(MetricType.Util, ts[idx], util_samples[idx], node[1]))
    se.interpolate()
    estVal = se.estimateAt(node[1])
    est1[t - W + 1] = estVal[0]
    est2[t - W + 1] = estVal[1] if len(estVal) > 1 else 0

print('\n=== UBO Change-Point Detection (2-class) ===')
print(f'Class 1: D={D1_before} (t<=100), D={D1_after} (t>100)')
print(f'Class 2: D={D2_val} (constant)')
print(f'Window size: {W}\n')
checkpoints = [W, 50, 80, 100, 110, 120, 130, 150, 180, 200]
print(f'{"t":>6}  {"Est D1":>8}  {"True D1":>8}  {"Est D2":>8}  {"True D2":>8}')
for cp in checkpoints:
    if W <= cp <= n:
        trueD1 = D1_before if cp <= changeT else D1_after
        print(f'{cp:6d}  {est1[cp-W]:8.4f}  {trueD1:8.4f}  {est2[cp-W]:8.4f}  {D2_val:8.4f}')
