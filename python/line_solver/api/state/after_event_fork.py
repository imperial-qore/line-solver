"""
afterEvent handler for stateful Fork nodes (FJ-augmented structs only).

Port of matlab/src/lang/+State/afterEventFork.m.
"""

import numpy as np
from ...constants import EventType


def after_event_fork(sn, ind, event, job_class, space_buf, space_srv, space_var,
                     is_simulation=False):
    """
    Handle afterEvent logic for stateful Fork nodes (FJ-augmented structs
    only, see ModelAdapter.fjtag). The fork state is a per-class count of
    parent jobs momentarily held before the fork firing. Arrivals are
    buffered here; the atomic multi-branch emission is not a DEP event but
    a fork firing synchronization (sn.fjsync) handled by after_fj_event.

    Returns (outspace, outrate, outprob).
    """
    space_srv = np.atleast_2d(np.array(space_srv, dtype=float))
    space_var = np.atleast_2d(np.array(space_var, dtype=float)) if np.size(space_var) else \
        np.zeros((space_srv.shape[0], 0))

    outspace = np.zeros((0, 0))
    outrate = np.zeros((0, 0))
    outprob = np.ones((1, 1))

    if event == EventType.ARV:
        srv = space_srv.copy()
        srv[:, job_class] = srv[:, job_class] + 1
        outspace = np.hstack([srv, space_var]) if space_var.shape[1] > 0 else srv
        # passive action, rate is unspecified
        outrate = -1.0 * np.ones((outspace.shape[0], 1))
        outprob = np.ones((outspace.shape[0], 1))
    elif event == EventType.DEP:
        # departures from a Fork are never generated as regular syncs
        # (see refreshSync); they occur only through sn.fjsync firings
        pass

    if is_simulation and outspace.shape[0] > 1:
        tot_rate = np.sum(outrate)
        cum_rate = np.cumsum(outrate) / tot_rate
        firing_ctr = int(np.sum(np.random.rand() > cum_rate.ravel()))
        firing_ctr = min(firing_ctr, outspace.shape[0] - 1)
        outspace = outspace[firing_ctr:firing_ctr + 1, :]
        outrate = np.array([[tot_rate]])
        outprob = outprob[firing_ctr:firing_ctr + 1, :]

    return outspace, outrate, outprob
