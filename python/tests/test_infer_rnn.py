"""Test RNN estimator."""
import pytest
import numpy as np

torch = pytest.importorskip("torch")

from line_solver import SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric
from line_solver.inference.api import infer_quick_model_rnn
from testutil import infer_generate_qlen_traces


def test_infer_rnn():
    np.random.seed(1)

    times_list = []
    traces_list = []
    for n_trial in range(2):
        model = infer_quick_model_rnn(
            False,
            [SchedStrategy.PS] * 5,
            [[150, 30, 90, 45, 60]],
            [20, 40, 60, 20, 35],
            [100]
        )
        timeI, QN = infer_generate_qlen_traces(model, 20)
        times_list.append(timeI)
        traces_list.append(QN)

    node = model.getNodes()
    jobclass = model.classes

    options = ParamEstimator.defaultOptions()
    options['method'] = 'rnn'
    se = ParamEstimator(model, options)

    for n_trial in range(len(times_list)):
        QN = traces_list[n_trial]
        timeI = times_list[n_trial]
        for i in range(5):
            aql = SampledMetric(MetricType.QLen, timeI, QN[i][0], node[i], jobclass[0])
            se.addSamples(aql)

    estVal = se.estimateAt(node)
    print('RNN estimated demands:', estVal)

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_rnn()
