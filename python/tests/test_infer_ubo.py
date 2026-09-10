"""Test UBO estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_ubo():
    np.random.seed(1)
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None, None]
    jobclass[0] = ClosedClass(model, 'Class1', 1, node[0], 0)
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
    arvr1_samples = 2.0 * np.ones(n) - np.random.rand(n) * 0.15
    arvr2_samples = 1.0 * np.ones(n) - np.random.rand(n) * 0.15
    util_samples = 0.1 * arvr1_samples + 0.3 * arvr2_samples
    respt1_samples = 0.1 / (1.0 - util_samples)
    respt2_samples = 0.3 / (1.0 - util_samples)

    options = ParamEstimator.defaultOptions()
    options['method'] = 'ubo'
    se = ParamEstimator(model, options)

    lambda1 = SampledMetric(MetricType.ArvR, ts, arvr1_samples, node[1], jobclass[0])
    lambda2 = SampledMetric(MetricType.ArvR, ts, arvr2_samples, node[1], jobclass[1])
    respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node[1], jobclass[0])
    respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node[1], jobclass[1])
    util = SampledMetric(MetricType.Util, ts, util_samples, node[1])

    se.addSamples(lambda1)
    se.addSamples(lambda2)
    se.addSamples(respT1)
    se.addSamples(respT2)
    se.addSamples(util)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    trueDemands = np.array([0.1, 0.3])
    # estVal may have multiple rows (iterations); use last row
    finalEst = estVal[-1, :] if estVal.ndim > 1 else estVal
    assert np.all(np.abs(finalEst - trueDemands) / trueDemands < 0.10), \
        f'UBO: estimated {finalEst} too far from true {trueDemands}'
    print(f'Estimated demands: Class1={finalEst[0]:.4f}, Class2={finalEst[1]:.4f}')

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_ubo()
