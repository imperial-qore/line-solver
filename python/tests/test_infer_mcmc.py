"""Test MCMC estimator."""
import numpy as np
from line_solver import Network, Delay, Queue, Exp, ClosedClass, SchedStrategy, SolverMVA, MetricType
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_mcmc():
    np.random.seed(1)
    trueDemand = 0.1
    model = Network('model')
    node = [None, None]
    node[0] = Delay(model, 'Delay')
    node[1] = Queue(model, 'Queue1', SchedStrategy.PS)
    jobclass = [None]
    jobclass[0] = ClosedClass(model, 'Class1', 1, node[0], 0)

    node[0].setService(jobclass[0], Exp.fitMean(1.0))
    node[1].setService(jobclass[0], Exp.fitMean(trueDemand))

    P = model.initRoutingMatrix()
    P.set(jobclass[0], jobclass[0], node[0], node[1], 1.0)
    P.set(jobclass[0], jobclass[0], node[1], node[0], 1.0)
    model.link(P)

    # Get true steady-state queue length from MVA
    solver_mva = SolverMVA(model)
    trueQLen = solver_mva.getAvgQLen()
    trueQLen_queue = trueQLen[1]  # Queue station (index 1, 0-based)

    # Generate model-consistent queue length samples
    n = 5000
    ts = np.arange(1, n + 1, dtype=float)
    qlen_samples = trueQLen_queue * np.ones(n) + np.random.rand(n) * 0.005 - 0.0025

    # Reset service for estimation
    node[1].setService(jobclass[0], Exp(float('nan')))

    # Estimate demands using MCMC
    options = ParamEstimator.defaultOptions()
    options['method'] = 'mcmc'
    se = ParamEstimator(model, options)

    ql = SampledMetric(MetricType.QLen, ts, qlen_samples, node[1])  # aggregate queue-length

    se.addSamples(ql)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    est = float(estVal) if np.ndim(estVal) == 0 else float(estVal.flat[0])
    assert abs(est - trueDemand) / trueDemand < 0.10, \
        f'MCMC: estimated {est:.4f}, relative error {100 * abs(est - trueDemand) / trueDemand:.2f}% exceeds 10%'
    print(f'Estimated demand: Class1={est:.4f}')

    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_mcmc()
