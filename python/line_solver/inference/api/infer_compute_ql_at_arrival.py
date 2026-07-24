import numpy as np


def infer_compute_ql_at_arrival(at, at_jobid, rt, rt_jobid, cls, R):
    """Compute per-class queue lengths at arrival.

    Reconstructs the queue state seen by each arriving job using arrival
    and departure times. At ties, departures are processed before arrivals.

    Args:
        at: arrival times (n,)
        at_jobid: job IDs for arrival times (n,)
        rt: response times (m,)
        rt_jobid: job IDs for response times (m,)
        cls: class of each arrival sample (n,), 1-based
        R: number of classes

    Returns:
        ql: n x R matrix of per-class queue lengths at each arrival

    Copyright (c) 2012-2026, Imperial College London
    All rights reserved.
    """
    at = np.asarray(at, dtype=float).flatten()
    at_jobid = np.asarray(at_jobid).flatten()
    rt = np.asarray(rt, dtype=float).flatten()
    rt_jobid = np.asarray(rt_jobid).flatten()
    cls = np.asarray(cls, dtype=int).flatten()
    n = len(at)

    # Match response times to arrivals by job ID
    rt_matched = np.zeros(n)
    jobid_to_rt = {}
    for idx, jid in enumerate(rt_jobid):
        jobid_to_rt[jid] = rt[idx]
    for idx, jid in enumerate(at_jobid):
        if jid not in jobid_to_rt:
            raise ValueError('Not all arrival job IDs found in response time job IDs.')
        rt_matched[idx] = jobid_to_rt[jid]

    # Sort arrivals by time
    sort_idx = np.argsort(at, kind='stable')
    at_sorted = at[sort_idx]
    cls_sorted = cls[sort_idx]
    rt_sorted = rt_matched[sort_idx]

    exit_times = at_sorted + rt_sorted

    # Event list: [time, type (-1=dep/+1=arv), sorted_idx, class]
    events = np.zeros((2 * n, 4))
    events[:n, 0] = at_sorted
    events[:n, 1] = 1  # arrival
    events[:n, 2] = np.arange(n)
    events[:n, 3] = cls_sorted

    events[n:, 0] = exit_times
    events[n:, 1] = -1  # departure
    events[n:, 2] = np.arange(n)
    events[n:, 3] = cls_sorted

    # Sort by time; departures (-1) before arrivals (+1) at same time
    order = np.lexsort((events[:, 1], events[:, 0]))
    events = events[order]

    state = np.zeros(R, dtype=int)
    ql_sorted = np.zeros((n, R), dtype=int)
    for i in range(2 * n):
        c = int(events[i, 3]) - 1  # 0-based class index
        if events[i, 1] == 1:  # arrival
            state[c] += 1
            ql_sorted[int(events[i, 2]), :] = state.copy()
        else:  # departure
            state[c] -= 1

    # Unsort back to original input order
    ql = np.zeros((n, R), dtype=int)
    ql[sort_idx, :] = ql_sorted

    return ql
