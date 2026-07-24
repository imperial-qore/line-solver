import numpy as np
from line_solver.inference.api.infer_get_qlen_arrival import infer_get_qlen_arrival
from line_solver.inference.api.infer_minps import infer_minps


def infer_minps_setup(data, init_sample, sample_size, V, model, node):
    """Setup input data for the MINPS estimation method and call it.

    Args:
        data: list of lists in standard format
        init_sample: first sample index (0-based)
        sample_size: number of samples to use (0 = all)
        V: number of servers
        model: LINE Network model
        node: PS queue node

    Returns:
        demand_est: 1-D array of estimated demands

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    R = len(data[2])

    sample_number = np.array([len(data[2][k]) for k in range(R)])

    # Remove classes without samples
    new_data = [[] for _ in range(6)]
    for k in range(R):
        if sample_number[k] > 0:
            for j in range(6):
                new_data[j].append(data[j][k])

    R = len(new_data[2])
    qls = infer_get_qlen_arrival(new_data)

    rt_all = []
    cls_all = []
    ql_all = []
    at_all = []
    for k in range(R):
        rt_all.append(np.asarray(new_data[3][k], dtype=float))
        cls_all.append(np.full(len(new_data[3][k]), k + 1, dtype=int))
        ql_all.append(qls[k])
        at_all.append(np.asarray(new_data[2][k], dtype=float) / 1000.0)

    rt = np.concatenate(rt_all)
    cls = np.concatenate(cls_all)
    ql = np.vstack(ql_all)
    at = np.concatenate(at_all)

    # Sort by arrival time
    all_times = np.column_stack([at, rt, cls, ql])
    order = np.argsort(all_times[:, 0])
    all_times = all_times[order]

    at = all_times[:, 0]
    rt = all_times[:, 1]
    cls = all_times[:, 2].astype(int)
    ql = all_times[:, 3:]

    if sample_size == 0:
        sample_size = ql.shape[0]

    first_sample = init_sample
    final_sample = init_sample + sample_size
    sample_set = slice(first_sample, final_sample)

    ql_exp = ql[sample_set, :]
    rt_exp = rt[sample_set]
    cls_exp = cls[sample_set]

    # Remove zero response times
    valid = rt_exp > 0
    rt_exp = rt_exp[valid]
    cls_exp = cls_exp[valid]
    ql_exp = ql_exp[valid, :]

    return infer_minps(model, node, rt_exp, cls_exp, ql_exp)
