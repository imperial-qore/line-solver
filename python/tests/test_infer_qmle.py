"""Test QMLE estimator."""
import numpy as np
from line_solver import (
    Network, Delay, Queue, Exp, ClosedClass,
    SchedStrategy, SolverMVA, MetricType,
)
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_qmle():
    np.random.seed(1)

    # define model with true demands (2 classes, 2 and 3 jobs)
    trueDemands = np.array([0.2, 0.4])
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None, None]
    jobclass[0] = ClosedClass(model, 'Class1', 2, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 3, node[0], 0)

    node[0].setService(jobclass[0], Exp.fitMean(1.0))
    node[0].setService(jobclass[1], Exp.fitMean(1.0))
    node[1].setService(jobclass[0], Exp.fitMean(trueDemands[0]))
    node[1].setService(jobclass[1], Exp.fitMean(trueDemands[1]))

    P = model.initRoutingMatrix()
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
    P.set(jobclass[1], jobclass[1], node[0], node[1], 1.0)
    P.set(jobclass[1], jobclass[1], node[1], node[0], 1.0)
    model.link(P)

    # get true steady-state queue lengths from MVA
    solver_mva = SolverMVA(model)
    trueQLen = solver_mva.getAvgQLen()
    trueQLen1 = trueQLen[1, 0]  # Queue station (index 1), Class 1 (index 0)
    trueQLen2 = trueQLen[1, 1]  # Queue station (index 1), Class 2 (index 1)

    # generate model-consistent queue length samples
    n = 5000
    ts = np.arange(1, n + 1, dtype=float)
    qlen1_samples = trueQLen1 * np.ones(n) + np.random.rand(n) * 0.02 - 0.01
    qlen2_samples = trueQLen2 * np.ones(n) + np.random.rand(n) * 0.02 - 0.01

    # reset service for estimation
    node[1].setService(jobclass[0], Exp(float('nan')))
    node[1].setService(jobclass[1], Exp(float('nan')))

    # estimate demands using QMLE
    options = ParamEstimator.defaultOptions()
    options['method'] = 'qmle'
    se = ParamEstimator(model, options)

    ql1 = SampledMetric(MetricType.QLen, ts, qlen1_samples, node[1], jobclass[0])
    ql2 = SampledMetric(MetricType.QLen, ts, qlen2_samples, node[1], jobclass[1])

    se.addSamples(ql1)
    se.addSamples(ql2)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    est = np.atleast_1d(estVal).flatten()
    assert np.all(np.abs(est[:2] - trueDemands) / trueDemands < 0.10), \
        f'QMLE: estimated [{est[0]:.4f}, {est[1]:.4f}] too far from true [{trueDemands[0]:.4f}, {trueDemands[1]:.4f}]'
    print(f'Estimated demands: Class1={est[0]:.4f}, Class2={est[1]:.4f}')

    # solve model
    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_qmle()
