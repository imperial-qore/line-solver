"""Example: Compare trace-based estimators (MLPS, FMLPS, Gibbs) on a PS station."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, MVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric

np.random.seed(1)
N = 5
model = Network('model')
node = [None] * 2
node[0] = Delay(model, 'Delay')
node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
jobclass = [None]
jobclass[0] = ClosedClass(model, 'Class1', N, node[0], 0)
node[0].setService(jobclass[0], Exp.fitMean(1.0))
node[1].setService(jobclass[0], Exp(float('nan')))
P = model.initRoutingMatrix()
P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
model.link(P)

# Generate synthetic trace data
n = 200
arrival_times = np.sort(np.random.rand(n) * 100)
response_times = 0.5 + np.random.rand(n) * 0.3
tput_samples = np.ones(n) * (N / (1.0 + 0.5))

arvData = SampledMetric(MetricType.ArvR, arrival_times, arrival_times, node[1], jobclass[0])
arvData.setTrace()
rtData = SampledMetric(MetricType.RespT, arrival_times, response_times, node[1], jobclass[0])
rtData.setTrace()
tputData = SampledMetric(MetricType.Tput, arrival_times, tput_samples, node[1], jobclass[0])

# MLPS
print('\n=== MLPS Estimator ===')
options = ParamEstimator.defaultOptions()
options['method'] = 'mlps'
se = ParamEstimator(model, options)
se.addSamples(arvData)
se.addSamples(rtData)
se.interpolate()
try:
    estVal_mlps = se.estimateAt(node[1])
    print(f'MLPS demand: Class1={estVal_mlps[0]:.4f}')
except Exception as e:
    print(f'MLPS skipped: {e}')
    estVal_mlps = float('nan')

# FMLPS
print('\n=== FMLPS Estimator ===')
node[1].setService(jobclass[0], Exp(float('nan')))
model.reset()
options['method'] = 'fmlps'
se = ParamEstimator(model, options)
se.addSamples(arvData)
se.addSamples(rtData)
se.interpolate()
try:
    estVal_fmlps = se.estimateAt(node[1])
    print(f'FMLPS demand: Class1={estVal_fmlps[0]:.4f}')
except Exception as e:
    print(f'FMLPS skipped: {e}')
    estVal_fmlps = float('nan')

# Gibbs
print('\n=== Gibbs Estimator ===')
node[1].setService(jobclass[0], Exp(float('nan')))
model.reset()
options['method'] = 'gibbs'
se = ParamEstimator(model, options)
se.addSamples(arvData)
se.addSamples(rtData)
se.addSamples(tputData)
se.interpolate()
try:
    estVal_gibbs = se.estimateAt(node[1])
    print(f'Gibbs demand: Class1={estVal_gibbs[0]:.4f}')
except Exception as e:
    print(f'Gibbs skipped: {e}')
    estVal_gibbs = float('nan')

# Compare
mlps_val = estVal_mlps[0] if hasattr(estVal_mlps, '__len__') else estVal_mlps
fmlps_val = estVal_fmlps[0] if hasattr(estVal_fmlps, '__len__') else estVal_fmlps
gibbs_val = estVal_gibbs[0] if hasattr(estVal_gibbs, '__len__') else estVal_gibbs
print(f'\n=== Comparison (true demand ~ 0.5) ===')
print(f'MLPS:  {mlps_val:.4f}')
print(f'FMLPS: {fmlps_val:.4f}')
print(f'Gibbs: {gibbs_val:.4f}')

solver = MVA(model)
print(f'\nSOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
