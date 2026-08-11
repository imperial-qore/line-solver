import numpy as np
from line_solver.inference.api.infer_compute_ql_at_arrival import infer_compute_ql_at_arrival


def infer_get_qlen_arrival(data):
    """Compute queue lengths at arrival from legacy cell-based data format.

    Wrapper around infer_compute_ql_at_arrival for the legacy format
    where data[3][k] contains arrival times (in ms) and data[4][k]
    contains response times for class k.

    Args:
        data: dict/list of lists in standard format (6 x K+1)
              data[2][k] = arrival times (ms), data[3][k] = response times

    Returns:
        ql: list of K arrays, each num_samples(k) x K

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    K = len(data[2])  # number of classes (excluding last column)

    at_all = []
    rt_all = []
    cls_all = []
    num_obs = []
    for k in range(K):
        nk = len(data[2][k])
        num_obs.append(nk)
        at_all.append(np.asarray(data[2][k], dtype=float) / 1000.0)  # ms -> s
        rt_all.append(np.asarray(data[3][k], dtype=float))
        cls_all.append(np.full(nk, k + 1, dtype=int))

    at = np.concatenate(at_all)
    rt = np.concatenate(rt_all)
    cls = np.concatenate(cls_all)

    n = len(at)
    jobid = np.arange(1, n + 1)
    ql_unsorted = infer_compute_ql_at_arrival(at, jobid, rt, jobid, cls, K)

    # Split into per-class arrays
    ql = []
    counter = 0
    for k in range(K):
        ql.append(ql_unsorted[counter:counter + num_obs[k], :])
        counter += num_obs[k]

    return ql
