"""
Fire a fork synchronization on a global state (FJ-augmented structs only).

Port of matlab/src/lang/+State/afterFJEvent.m.
"""

import numpy as np
from ...constants import EventType, GlobalConstants


def after_fj_event(sn, fjentry, glspace, is_simulation=False):
    """
    Fire a fork synchronization (one entry of sn.fjsync, built by
    ModelAdapter.fjtag) on the global state GLSPACE, a list with the local
    state of every stateful node. The firing atomically consumes one parent
    job of class r held at the (stateful) Fork node and emits one sibling
    per branch, in the auxiliary classes of the entry's tag, at the branch
    head nodes.

    Enabling condition:
     1. the Fork holds at least one class-r parent job;
     2. the entry's tag is the LOWEST free tag for this (fork, class): a tag
        is free iff its auxiliary classes have zero occupancy network-wide.

    Returns (out_global_states, outrate, outprob):
    - out_global_states: list of successor global states (each a list of
      per-stateful-node 1D arrays)
    - outrate: np.ndarray of firing rates (immediate)
    - outprob: np.ndarray of outcome probabilities
    """
    from .after_event import after_event
    from .marginal import toMarginalAggr

    R = int(sn.nclasses)
    f = int(fjentry['fork'])
    r = int(fjentry['class'])
    isf_f = int(sn.nodeToStateful[f])

    empty = ([], np.zeros(0), np.zeros(0))

    # 1. parent job held at the fork
    forkstate = np.asarray(glspace[isf_f], dtype=float).ravel()
    if forkstate[-R + r] < 1:
        return empty

    # 2. lowest-free-tag test: occupancy of each tag's auxiliary classes,
    # scanned across all stateful nodes
    auxall = np.asarray(fjentry['auxall'], dtype=int)  # B x T aux class indices
    B, T = auxall.shape
    t = int(fjentry['tag'])
    nglobal = np.zeros(R)
    for isf in range(int(sn.nstateful)):
        _, nir = toMarginalAggr(sn, int(sn.statefulToNode[isf]), np.atleast_2d(glspace[isf]))
        nglobal = nglobal + np.asarray(nir).ravel()[:R]
    occ = np.sum(nglobal[auxall], axis=0)  # length T
    if occ[t] > 0:
        return empty  # tag in use
    if np.any(occ[:t] == 0):
        return empty  # a lower tag is free: that entry fires instead

    # consume the parent job at the fork
    newgl = [np.array(g, dtype=float).ravel().copy() for g in glspace]
    newgl[isf_f][-R + r] = newgl[isf_f][-R + r] - 1

    # emit fjentry.weight (tasksPerLink) siblings per branch: sequential
    # application over the partial outcome list handles branches sharing the
    # same head node, repeated emissions on the same branch, and expands
    # phase-entry mixtures of non-exponential sibling services
    branchheads = np.asarray(fjentry['branchheads'], dtype=int).ravel()
    auxclasses = np.asarray(fjentry['auxclasses'], dtype=int).ravel()
    weight = int(fjentry.get('weight', 1))
    partials = [newgl]
    partprob = [1.0]
    emissions = list(range(B)) * weight
    for b in emissions:
        bh = int(branchheads[b])
        isf_b = int(sn.nodeToStateful[bh])
        a = int(auxclasses[b])
        newpartials = []
        newpartprob = []
        for pp in range(len(partials)):
            curgl = partials[pp]
            arvspace, _, arvprob, _, _ = after_event(
                sn, bh, np.atleast_2d(curgl[isf_b]), EventType.ARV, a, is_simulation)
            arvspace = np.atleast_2d(arvspace)
            if arvspace.size == 0:
                # sibling arrival blocked (cannot occur under the per-tag
                # auxiliary class capacity invariant); disable the firing
                return empty
            arvprob = np.asarray(arvprob).ravel()
            for io in range(arvspace.shape[0]):
                nextgl = [g.copy() for g in curgl]
                nextgl[isf_b] = arvspace[io, :].copy()
                newpartials.append(nextgl)
                if io < len(arvprob):
                    newpartprob.append(partprob[pp] * arvprob[io])
                else:
                    newpartprob.append(partprob[pp])
        partials = newpartials
        partprob = newpartprob

    out_global_states = partials
    outprob = np.asarray(partprob, dtype=float) * float(fjentry.get('prob', 1.0))
    outrate = GlobalConstants.Immediate * np.ones(len(partials))

    if is_simulation:
        if len(partials) > 1:
            cum_prob = np.cumsum(outprob) / np.sum(outprob)
            firing_ctr = int(np.sum(np.random.rand() > cum_prob))
            firing_ctr = min(firing_ctr, len(partials) - 1)
            out_global_states = [out_global_states[firing_ctr]]
            outrate = np.array([outrate[firing_ctr]])
        # the phase-entry choice has already been sampled inside after_event,
        # so the firing competes at the full immediate rate
        outprob = np.array([1.0])

    return out_global_states, outrate, outprob
