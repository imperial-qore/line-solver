import numpy as np
from scipy.optimize import nnls


def infer_rps(rt, cls, ql, V):
    """Regression for Processor Sharing (RPS) demand estimation.

    Based on mean-value analysis for PS stations:
        E[R_r] = E[D_r] * E[Q_bar_A] / V

    Args:
        rt: response time samples (n,)
        cls: class of each sample (n,), 1-based
        ql: queue lengths at arrival (n x R)
        V: number of servers

    Returns:
        demand_est: 1-D array of estimated demands (R,)

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    rt = np.asarray(rt, dtype=float).flatten()
    cls = np.asarray(cls, dtype=int).flatten()
    ql = np.asarray(ql, dtype=float)
    R = int(np.max(cls))

    demand_est = np.zeros(R)
    for r in range(1, R + 1):
        idx = cls == r
        resp_times = rt[idx]
        q_bar_a = (np.sum(ql[idx, :], axis=1) + 1) / V
        demand_est[r - 1], _ = nnls(q_bar_a.reshape(-1, 1), resp_times)
    return demand_est
