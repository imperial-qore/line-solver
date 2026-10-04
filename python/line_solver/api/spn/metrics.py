"""
Stationary measures of a product-form stochastic Petri net from the MDD-rec
masses.

S. Balsamo, A. Marin, I. Stojic, FGCS 111 (2020) 475-490, Sec. 3.1 for the
definitions and Sec. 5.3 for the recursions they are read off.

  n(P_j) = sum_k k P(m_j = k)                    mean tokens
  u(P_j) = 1 - P(m_j = 0)                        place utilization
  u(T_j) = P(e_j >= 1)                           transition utilization
  x(T_j) = sum_k min(k, c_j) W(T_j) P(e_j = k)   throughput
  x(P_j) = sum_T I_j(T) x(T)                     tokens removed per unit time

ONE DEVIATION FROM THE PAPER'S x(T_j), AND IT IS A GENERALISATION. The paper
writes x(T_j) = sum_k k W(T_j) P(e_j = k), which is INFINITE-SERVER firing
semantics -- every enabling set fires in parallel. LINE's own rate law is
min(enabling degree, nmodeservers) * W(T), so c_j above is the mode's server
count: c_j = 1 recovers single-server semantics, x = W(T) P(e >= 1), and
c_j = infinity recovers the paper's formula exactly. Using the paper's form for a
single-server mode would report a throughput that grows with the token
population of a net whose transition can only fire one set at a time.

The measures come out of ONE reachable set and ONE set of g_l, so they are
mutually consistent by construction: no per-measure fixed point, no iteration.

See also: mdd_rec, spn_rec_enabled, spn_mdd.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, Sequence

import numpy as np

from ..io.logging import line_error
from ..mdd.rec import mdd_rec, mdd_rec_marginal
from .rec_enabled import spn_rec_enabled


def spn_metrics(mdds, g: Sequence[Sequence[float]], info: Dict) -> Dict[str, object]:
    """Every measure of Sec. 3.1 from one diagram and one product form.

    Parameters
    ----------
    mdds : the reachable set built by spn_mdd
    g : per-level product-form factors g_l(v)
    info : the metadata spn_mdd returned alongside the diagram

    Returns
    -------
    dict with 'G', 'tokens', 'placeUtil', 'placeTput', 'modeUtil', 'modeTput'
    and 'marginal', the per-level P(m_l = k).
    """
    L = int(info['nplacelevels'])
    modes = info['modes']

    G = mdd_rec(mdds, g)
    if not G > 0:
        line_error('spn_metrics',
                   'the normalising constant is not positive; the g_l passed do not describe a '
                   'product form over this reachable set')

    marginal = []
    tokens = np.zeros(L)
    place_util = np.zeros(L)
    place_tput = np.zeros(L)
    for l in range(L):
        pk = mdd_rec_marginal(mdds, g, l) / G
        marginal.append(pk)
        tokens[l] = float(np.dot(np.arange(pk.size), pk))
        place_util[l] = 1.0 - pk[0]

    E = len(modes)
    mode_util = np.zeros(E)
    mode_tput = np.zeros(E)
    for e in range(E):
        mde = modes[e]
        en = spn_rec_enabled(mdds, g, mde, L)
        mode_util[e] = en['ge'][1] / G
        # W(T) is the scalar firing rate of the mode; a phase-type firing time
        # has no single rate, so its throughput is left to the phase-level
        # marginal rather than reported through this formula.
        if int(mde['nph']) > 1:
            line_error('spn_metrics',
                       'mode %d of node %d has a phase-type firing time, whose throughput is not '
                       'W(T) times an enabling probability; read it from the phase-level marginal '
                       'instead' % (int(mde['mode']) + 1, int(mde['trans']) + 1))
        # The formula below is W(T)*E[min(enabling degree, servers)], which is the
        # rate law only when no marking-dependent multiplier is in play. With one,
        # the firing rate is not a function of the enabling degree at all, so the
        # enabling-degree law is the wrong summary to take it from. The MARGINALS
        # above are unaffected -- they come from the product form, not the rates.
        if mde.get('dep') is not None:
            line_error('spn_metrics',
                       'mode %d of node %d has a marking-dependent firing rate, so its '
                       'throughput is not W(T) times a function of the enabling degree and '
                       'cannot be read from the enabling-degree law. The token marginals are '
                       'still exact' % (int(mde['mode']) + 1, int(mde['trans']) + 1))
        rate = float(np.asarray(mde['D1'], dtype=float).ravel()[0])
        srv = float(mde['srv'])
        k = np.arange(en['eq'].size, dtype=float)
        served = k if np.isinf(srv) else np.minimum(k, srv)
        x = float(np.dot(served * rate, en['eq']) / G)
        mode_tput[e] = x
        enab = np.asarray(mde['enab'], dtype=float).ravel()
        for l in range(L):
            if enab[l] > 0:
                place_tput[l] += enab[l] * x
    return {
        'G': G,
        'tokens': tokens,
        'placeUtil': place_util,
        'placeTput': place_tput,
        'modeUtil': mode_util,
        'modeTput': mode_tput,
        'marginal': marginal,
    }
