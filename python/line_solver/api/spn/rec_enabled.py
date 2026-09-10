"""
Enabling-degree distribution of one mode of a product-form stochastic Petri net,
by the masked MDD-rec recursion.

S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 5.3.

The enabling degree of a mode in marking m is

    e(m) = min_{l : I_l > 0} floor(m_l / I_l),

zero when any inhibitor threshold is met. P(e >= k) is therefore the mass of the
marking subset in which EVERY input level holds at least k*I_l tokens and no
inhibitor fires, which is a per-level restriction and so exactly what
mdd_rec_masked computes: the paper's second modified recurrence is the same walk
under a different mask, not a second algorithm.

The masses returned are UNNORMALISED, as in the paper; divide by G from mdd_rec
for probabilities. spn_metrics does that and turns them into the transition
measures.

See also: mdd_rec, spn_metrics, spn_mdd.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, Sequence

import numpy as np

from ..io.logging import line_error
from ..mdd.rec import mdd_rec_masked


def spn_rec_enabled(mdds, g: Sequence[Sequence[float]], mde: Dict,
                    nplacelevels: int) -> Dict[str, object]:
    """Unnormalised enabling-degree masses of one mode.

    Parameters
    ----------
    mdds : the reachable set built by spn_mdd
    g : per-level product-form factors, one sequence per level
    mde : the mode, as returned in info['modes'] by spn_mdd
    nplacelevels : how many leading levels are place levels

    Returns
    -------
    dict with 'ge' (mass of {e >= k}), 'eq' (mass of {e == k}) and 'maxDegree',
    the largest enabling degree the place bounds permit (E_j in the paper).
    """
    K = int(mdds.K)
    L = int(nplacelevels)
    if L > K:
        line_error('spn_rec_enabled', 'more place levels than diagram levels')
    enab = np.asarray(mde['enab'], dtype=float).ravel()
    inhib = np.asarray(mde['inhib'], dtype=float).ravel()

    # E_j: the enabling degree cannot exceed what the tightest input place bound
    # allows. A mode with no input place has no bound and is refused rather than
    # silently truncated, matching spn_mdd's own refusal.
    emax = 0
    has_input = False
    for l in range(L):
        if not enab[l] > 0:
            continue
        cap = int(np.floor((int(mdds.domain[l]) - 1) / enab[l]))
        emax = cap if not has_input else min(emax, cap)
        has_input = True
    if not has_input:
        line_error('spn_rec_enabled',
                   'the mode consumes from no place, so its enabling degree is unbounded')

    ge = np.zeros(emax + 2)
    for k in range(emax + 1):
        mask = [np.ones(int(mdds.domain[j]), dtype=bool) for j in range(K)]
        for l in range(L):
            need = enab[l] * k
            v = np.arange(int(mdds.domain[l]), dtype=float)
            # k = 0 asks only that the marking exist, so the inhibitor test
            # belongs to k >= 1: e = 0 covers the inhibited markings too.
            drop = v < need
            if k > 0:
                drop = np.logical_or(drop, v >= inhib[l])
            mask[l][drop] = False
        ge[k] = mdd_rec_masked(mdds, g, mask)
    eq = np.zeros(emax + 2)
    eq[:emax + 1] = ge[:emax + 1] - ge[1:emax + 2]
    return {'ge': ge, 'eq': eq, 'maxDegree': emax}
