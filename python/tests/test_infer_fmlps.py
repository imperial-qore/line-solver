"""Test FMLPS estimator."""
import numpy as np
from line_solver import (
    Network, Delay, Queue, Exp, ClosedClass,
    SchedStrategy, SolverMVA, SolverSSA, MetricType,
)
from line_solver.inference import ParamEstimator, SampledMetric


def test_infer_fmlps():
    np.random.seed(1)

    # define model with true demand (1 job, single class)
    trueDemand = 0.5
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

    # generate trace data from SSA simulation
    # Use fewer events than MATLAB (10000) since Python ODE solver is slower
    solver_ssa = SolverSSA(model, seed=1, samples=500)
    samplePath = solver_ssa.sample(node[1], 500)

    # extract per-class arrival and departure times at the queue
    queue_node_idx = node[1].get_index()
    arv_times = []
    dep_times = []
    for ev in samplePath.event:
        if ev.node == queue_node_idx and ev.class_idx == 1:
            if ev.event == 'ARV':
                arv_times.append(ev.t)
            elif ev.event == 'DEP':
                dep_times.append(ev.t)
    n = min(len(arv_times), len(dep_times))
    arv_times = np.array(arv_times[:n])
    dep_times = np.array(dep_times[:n])
    response_times = dep_times - arv_times
    arrival_times = arv_times

    # reset service for estimation
    node[1].setService(jobclass[0], Exp(float('nan')))

    # create trace-format SampledMetric objects
    arvData = SampledMetric(MetricType.ArvR, arrival_times, arrival_times, node[1], jobclass[0])
    arvData.setTrace()
    rtData = SampledMetric(MetricType.RespT, arrival_times, response_times, node[1], jobclass[0])
    rtData.setTrace()

    # estimate demands
    options = ParamEstimator.defaultOptions()
    options['method'] = 'fmlps'
    se = ParamEstimator(model, options)
    se.addSamples(arvData)
    se.addSamples(rtData)
    se.interpolate()
    estVal = se.estimateAt(node[1])

    est = float(estVal) if np.ndim(estVal) == 0 else float(estVal.flat[0])
    assert abs(est - trueDemand) / trueDemand < 0.10, \
        f'FMLPS: estimated {est:.4f}, relative error {100 * abs(est - trueDemand) / trueDemand:.2f}% exceeds 10%'

    # solve model
    solver = SolverMVA(model)
    avgTable = solver.getAvgTable()
    print(avgTable)


if __name__ == '__main__':
    test_infer_fmlps()
