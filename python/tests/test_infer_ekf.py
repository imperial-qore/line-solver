"""Test EKF estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_ekf():
    np.random.seed(1)
    trueDemand = 0.3
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None]
    jobclass[0] = ClosedClass(model, 'Class1', 5, node[0], 0)

    node[0].setService(jobclass[0], Exp.fitMean(1.0))
    node[1].setService(jobclass[0], Exp.fitMean(trueDemand))

    P = model.initRoutingMatrix()
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
    model.link(P)

    # Get true steady-state metrics from MVA
    solver_mva = SolverMVA(model)
    trueRespT = solver_mva.getAvgRespT()
    trueUtil = solver_mva.getAvgUtil()
    trueTput = solver_mva.getAvgTput()

    stIdx = node[1].getStationIndex()
    trueR = trueRespT[stIdx - 1, 0]
    trueU = trueUtil[stIdx - 1, 0]
    trueX = trueTput[stIdx - 1, 0]

    # Generate noisy dataset consistent with model steady state
    n = 1000
    ts = np.arange(1, n + 1, dtype=float)
    noise_scale = 0.05
    arvr_samples = trueX * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueX
    respt_samples = trueR * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueR
    util_samples = trueU * np.ones(n) + (np.random.rand(n) - 0.5) * noise_scale * trueU

    # Reset service for estimation
    node[1].setService(jobclass[0], Exp(float('nan')))

    # Estimate demands
    options = ParamEstimator.defaultOptions()
    options['method'] = 'ekf'
    se = ParamEstimator(model, options)

    lambda1 = SampledMetric(MetricType.ArvR, ts, arvr_samples, node[1], jobclass[0])
    respT1 = SampledMetric(MetricType.RespT, ts, respt_samples, node[1], jobclass[0])
    util = SampledMetric(MetricType.Util, ts, util_samples, node[1])
    se.addSamples(lambda1)
    se.addSamples(respT1)
    se.addSamples(util)

    estVal = se.estimateAt(node[1])

    assert abs(estVal - trueDemand) / trueDemand < 0.10, \
        f'EKF: estimated {estVal:.4f}, relative error {100 * abs(estVal - trueDemand) / trueDemand:.2f}% exceeds 10%'

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_ekf()
