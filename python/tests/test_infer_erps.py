"""Test ERPS estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType, EventType
from line_solver.inference import Event
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_erps():
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
    arvr1_samples = 1.0 * np.ones(n) - np.random.rand(n) * 0.15
    arvr2_samples = 2.0 * np.ones(n) - np.random.rand(n) * 0.15
    util_samples = 0.1 * arvr1_samples + 0.3 * arvr2_samples
    respt1_samples = 0.1 / (1.0 - util_samples)
    respt2_samples = 0.3 / (1.0 - util_samples)
    aqlen1_samples = 1.0 + util_samples / (1.0 - util_samples)
    aqlen2_samples = 1.0 + util_samples / (1.0 - util_samples)

    options = ParamEstimator.defaultOptions()
    options['method'] = 'erps'
    se = ParamEstimator(model, options)

    aql1 = SampledMetric(MetricType.QLen, ts, aqlen1_samples, node[1])  # aggregate queue-length
    aql1.setConditional(Event(EventType.ARV, node[1], jobclass[0]))  # conditional on class-1 arrivals

    aql2 = SampledMetric(MetricType.QLen, ts, aqlen2_samples, node[1])  # aggregate queue-length
    aql2.setConditional(Event(EventType.ARV, node[1], jobclass[1]))  # conditional on class-2 arrivals

    respT1 = SampledMetric(MetricType.RespT, ts, respt1_samples, node[1], jobclass[0])
    respT2 = SampledMetric(MetricType.RespT, ts, respt2_samples, node[1], jobclass[1])

    se.addSamples(aql1)
    se.addSamples(aql2)
    se.addSamples(respT1)
    se.addSamples(respT2)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    trueDemands = np.array([0.1, 0.3])
    assert np.all(np.abs(estVal - trueDemands) / trueDemands < 0.10), \
        f'ERPS: estimated {estVal} too far from true {trueDemands}'
    print(f'Estimated demands: Class1={estVal[0]:.4f}, Class2={estVal[1]:.4f}')

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_erps()
