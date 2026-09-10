"""Test UBR estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_ubr():
    np.random.seed(1)
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None, None]
    jobclass[0] = ClosedClass(model, 'Class1', 1, node[0], 0)
    jobclass[1] = ClosedClass(model, 'Class2', 2, node[0], 0)

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
    arvr1_samples = 1.5 + np.random.rand(n) * 0.1
    arvr2_samples = 2.0 + np.random.rand(n) * 0.1
    util_samples = 0.2 * arvr1_samples + 0.4 * arvr2_samples

    options = ParamEstimator.defaultOptions()
    options['method'] = 'ubr'
    se = ParamEstimator(model, options)

    lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node[1], jobclass[0])
    lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node[1], jobclass[1])
    util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

    se.addSamples(lambda1)
    se.addSamples(lambda2)
    se.addSamples(util)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    trueDemands = np.array([0.2, 0.4])
    assert np.all(np.abs(estVal - trueDemands) / trueDemands < 0.10), \
        f'UBR: estimated {estVal} too far from true {trueDemands}'

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_ubr()
