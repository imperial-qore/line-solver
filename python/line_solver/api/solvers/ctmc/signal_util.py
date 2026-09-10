"""Utilization helpers for G-network signals under CTMC.

A job annihilated by a signal leaves the station without a service completion.
Two consequences for the utilization estimators:

* an arrival-based estimator lambda*E[S]/c measures the OFFERED load, so it is
  invalid for a class a signal can remove (native Python never uses it, but the
  MATLAB and JAR analyzers do -- see ctmc_signal_lossy.m / CtmcSignalLossy.java);
* the departure-based estimator T*E[S]/c is exact only under exponential
  service, since a job destroyed mid-service leaves busy time behind with no
  completion. With phase-type service it under-counts (M/Er2/1 with
  lambda+ = 0.5, lambda- = 0.4 gives 0.34941 against a true 0.37696).

`busy_fraction` therefore reads the busy-server occupancy off the enumerated
state space, which is exact for any service process and coincides with
T*E[S]/c when service is exponential.

Mirrors matlab/src/solvers/CTMC/ctmc_signal_lossy.m and ctmc_signal_busy.m.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np

from ...sn import SchedStrategy


def signal_lossy_classes(sn, arv_by_class):
    """Classes a G-network signal can annihilate at the station of interest.

    Args:
        sn: network structure
        arv_by_class: stationary arrival rate of each class at that station

    Returns:
        boolean array of length nclasses
    """
    K = len(arv_by_class)
    lossy = np.zeros(K, dtype=bool)
    issignal = getattr(sn, 'issignal', None)
    if issignal is None:
        return lossy
    issignal = np.asarray(issignal, dtype=bool).flatten()
    if issignal.size < K:
        issignal = np.pad(issignal, (0, K - issignal.size))
    if not issignal.any():
        return lossy
    target = getattr(sn, 'signaltarget', None)
    target = np.asarray(target).flatten() if target is not None else None
    for r in range(K):
        if not issignal[r] or arv_by_class[r] <= 0:
            continue
        # sn.signaltarget is 0-based in native Python; -1 means untargeted.
        tgt = int(target[r]) if (target is not None and r < target.size) else -1
        if 0 <= tgt < K:
            lossy[tgt] = True
        else:
            # Untargeted signals are class-agnostic: they remove any non-signal
            # job, matching signal_removal.py, MAM and LDES.
            lossy[~issignal[:K]] = True
    return lossy


def busy_fraction(sn, ind, ist, sched, nservers, space_isf, state_space_hashed,
                  isf, pi, K):
    """Exact per-class busy-server fraction at a station.

    PS-like disciplines share the servers among all resident jobs, so class k
    gets the weighted share n_k w_k / sum_j n_j w_j of the busy servers; the
    remaining disciplines expose the in-service indicator through toMarginal.
    """
    from ...state.marginal import toMarginal

    UNb = np.zeros(K)
    if space_isf is None:
        return UNb
    is_ps = sched in (SchedStrategy.PS, SchedStrategy.DPS, SchedStrategy.GPS,
                      SchedStrategy.LPS)
    weights = np.ones(K)
    if is_ps and getattr(sn, 'schedparam', None) is not None:
        for r in range(K):
            sp = sn.schedparam[ist, r]
            if sp is not None and float(sp) > 0:
                weights[r] = float(sp)
    for s_idx in range(len(pi)):
        if pi[s_idx] == 0:
            continue
        h = int(state_space_hashed[s_idx, isf])
        ni, nir, sir, _ = toMarginal(sn, ind, space_isf[h:h + 1])
        tot = float(np.sum(np.atleast_1d(ni)))
        if tot <= 0:
            continue
        if is_ps:
            nir = np.atleast_1d(nir).flatten()
            wtot = float(np.sum(nir[:K] * weights))
            if wtot > 0:
                busy = min(tot, nservers) / nservers
                for k in range(K):
                    UNb[k] += pi[s_idx] * (nir[k] * weights[k] / wtot) * busy
        else:
            sir = np.atleast_1d(sir).flatten()
            for k in range(min(K, len(sir))):
                UNb[k] += pi[s_idx] * sir[k] / nservers
    return UNb
