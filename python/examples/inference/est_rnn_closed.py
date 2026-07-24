"""Example: RNN-based estimation from queue-length traces."""
import numpy as np
from line_solver import SchedStrategy, MVA, Exp, MetricType
from line_solver.inference import ParamEstimator, SampledMetric
from line_solver.inference.api import infer_quick_model_rnn
from testutil import infer_generate_qlen_traces

np.random.seed(1)
model = infer_quick_model_rnn(
    False,
    [SchedStrategy.FCFS] * 3,
    [[150, 30, 90]],
    [20, 40, 60],
    [100]
)
node = model.getNodes()
jobclass = model.classes

print('\n=== Generating queue-length traces ===')
times_list = []
traces_list = []
for trial in range(2):
    model_tmp = infer_quick_model_rnn(
        False,
        [SchedStrategy.FCFS] * 3,
        [[150, 30, 90]],
        [20, 40, 60],
        [100]
    )
    timeI, QN = infer_generate_qlen_traces(model_tmp, 20)
    times_list.append(timeI)
    traces_list.append(QN)

# Reset service rates
for i in range(len(node)):
    try:
        node[i].setService(jobclass[0], Exp(float('nan')))
    except Exception:
        pass  # skip Source/Sink

print('\n=== RNN Estimator ===')
options = ParamEstimator.defaultOptions()
options['method'] = 'rnn'
se = ParamEstimator(model, options)

for trial in range(len(times_list)):
    QN = traces_list[trial]
    timeI = times_list[trial]
    for i in range(min(3, len(QN))):
        ql = SampledMetric(MetricType.QLen, timeI, QN[i][0], node[i], jobclass[0])
        se.addSamples(ql)

estVal = se.estimateAt(node)
print('RNN estimated demands:')
for i in range(estVal.shape[0]):
    print(f'  Node {i}: {estVal[i, 0]:.4f}')

solver = MVA(model)
print(f'\nSOLVER: {solver.getName()}')
avgTable = solver.getAvgTable()
print(avgTable)
