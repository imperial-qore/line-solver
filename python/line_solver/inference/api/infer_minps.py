import numpy as np

from line_solver.inference.api.infer_mlps import infer_mlps
from line_solver.inference.api.infer_rps import infer_rps


def infer_minps(model, node, rt, cls, ql):
    """MINPS demand estimation method.

    Runs both MLPS and RPS estimators and selects the one with the
    smaller mean demand estimate.

    Args:
        model: LINE Network model
        node: PS queue node
        rt: response time samples (n,)
        cls: class of each sample (n,), 1-based
        ql: queue lengths at arrival (n x R)

    Returns:
        demand_est: 1-D array of estimated demands (R,)

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    V = node.getNumberOfServers()

    demand_mlps = infer_mlps(model, node, rt, cls, ql)
    demand_rps = infer_rps(rt, cls, ql, V)

    if np.mean(demand_mlps) < np.mean(demand_rps):
        return demand_mlps
    else:
        return demand_rps
