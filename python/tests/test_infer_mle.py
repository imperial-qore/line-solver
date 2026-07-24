"""Test MLE estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_mle():
    np.random.seed(1)
    trueDemands = np.array([0.1, 0.3])
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None, None]
    jobclass[0] = ClosedClass(model, 'Class1', 1, node[0], 0)
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

    # Get true steady-state metrics from MVA
    solver_mva = SolverMVA(model)
    trueRespT = solver_mva.getAvgRespT()
    trueUtil = solver_mva.getAvgUtil()
    trueTput = solver_mva.getAvgTput()

    stIdx = node[1].getStationIndex()
    trueR1 = trueRespT[stIdx - 1, 0]
    trueR2 = trueRespT[stIdx - 1, 1]
    trueU = np.sum(trueUtil[stIdx - 1, :])
    trueX1 = trueTput[stIdx - 1, 0]
    trueX2 = trueTput[stIdx - 1, 1]

    # Generate noisy dataset consistent with model steady state
    n = 1000
    ts = np.arange(1, n + 1, dtype=float)
    noise_scale = 0.05
    arvr1_samples = trueX1 * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueX1
    arvr2_samples = trueX2 * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueX2
    respt1_samples = trueR1 * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueR1
    respt2_samples = trueR2 * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueR2
    util_samples = trueU * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueU

    # Reset service for estimation
    node[1].setService(jobclass[0], Exp(float('nan')))
    node[1].setService(jobclass[1], Exp(float('nan')))

    # Estimate demands
    options = ParamEstimator.defaultOptions()
    options['method'] = 'mle'
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

    assert np.all(np.abs(estVal - trueDemands) / trueDemands < 0.10), \
        f'MLE: estimated {estVal} too far from true {trueDemands}'
    print(f'Estimated demands: Class1={estVal[0]:.4f}, Class2={estVal[1]:.4f}')

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_mle()
