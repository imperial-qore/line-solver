"""Moment-relaxation LP bounds for a stochastic Petri net (Liu 1998).

The polytope and the LP are api.spn.spn_lpbnd; this module maps the LINE model
onto them and reads one side of the bracket back per place and class. Native
port of matlab/src/solvers/BA/solver_ba_spnlp_analyzer.m.

METHOD NAMES. Four, in two families:

  spnlp.upper     Markovian LP, upper side          exponential firing
  spnlp.lower     Markovian LP, lower side          exponential firing
  spnlp.op.upper  operational LP, upper side        any phase-type firing
  spnlp.op.lower  operational LP, lower side        any phase-type firing

The operational variant drops the second-moment, covariance and Little's-law
families and the whole E[X_p e_t] block with them, which is what removes the
exponential requirement. It is much looser, and is the reference's own "without
Markovian assumption" column.

BOUND CONVENTION. Q[i,r] is the reported side of the bracket on the mean number
of class-r tokens in place i. T[i,r] is the same side of the bracket on the
token throughput of that place, and R follows by Little's law from the two.
U[i,r] = Q[i,r] DELIBERATELY: a Place is an INF station and LINE reports U = Q
at an infinite server, which is what SolverCTMC and solver_nc_spn_analyzer both
do on the same net. The reference's place utilization 1 - P(m = 0) is a
different quantity and is not this column.

A TRANSITION GETS NO ROW. It is a StatefulNode and not a Station, so it has no
station index; the mode throughputs and enabling probabilities the LP also
brackets stay inside spn_lpbnd's return value, the same way spn_metrics keeps
modeTput and modeUtil off the table.

HOW TIGHT. The reference's own Table 2 measures it on a four-server production
line: the upper side lands 2% to 11% above simulation and the lower side 30% to
40% below it, both comfortably inside the operational bounds it also reports.
Expect a usable upper bound and a weak lower one.

Reference:
    Z. Liu (1998). Performance analysis of stochastic timed Petri nets using
    linear programming approach. IEEE Transactions on Software Engineering
    24(11), 1014-1030.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

import numpy as np

BA_SPNLP = {'spnlp.upper', 'spnlp.lower', 'spnlp.op.upper', 'spnlp.op.lower'}

_SPNLP_SPEC = {
    'spnlp.upper': (True, 1),
    'spnlp.lower': (True, 0),
    'spnlp.op.upper': (False, 1),
    'spnlp.op.lower': (False, 0),
}


def spnlp_bound(sn, method, options=None):
    """Bracket one side of the Petri-net moment relaxation, per place and class."""
    from ...api.spn import spn_lpbnd
    from ...api.sn.network_struct import NodeType, SchedStrategy

    if method not in _SPNLP_SPEC:
        raise ValueError("Unknown SPN bound method '%s'. Valid: spnlp.upper, spnlp.lower, "
                         "spnlp.op.upper, spnlp.op.lower." % method)
    markovian, side = _SPNLP_SPEC[method]

    M = int(sn.nstations)
    K = int(sn.nclasses)
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    C = np.zeros((1, K))
    X = np.zeros((1, K))

    # ----- model gates -----
    # Every other check belongs to spn_lpbnd, which refuses by name on the mode
    # it cannot represent. What must be decided here is only whether this is a
    # Petri net at all, and whether the places carry an embedded queue the
    # relaxation has no variable for.
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    if not np.any(nodetype == int(NodeType.TRANSITION)):
        raise ValueError("Method '%s' bounds a stochastic Petri net; this model has no "
                         "Transition node. Use the queueing-network bound families, or "
                         "SolverMVA/SolverNC." % method)
    nodetostation = np.ravel(np.asarray(sn.nodeToStation, dtype=int))
    for nd in np.nonzero(nodetype == int(NodeType.PLACE))[0]:
        ist = int(nodetostation[nd])
        if ist >= 0 and sn.sched and sn.sched[ist] != SchedStrategy.INF:
            raise ValueError(
                "Method '%s' does not support queueing places: place %s serves under a "
                "non-INF discipline, and the relaxation carries one variable per "
                "(place, class) marking with no notion of an embedded queue."
                % (method, sn.nodenames[nd]))

    lpopt = {'markovian': markovian}
    config = getattr(options, 'config', None) if options is not None else None
    if config is not None:
        # The reference's liveness rows hold only on a live net, which
        # spn_lpbnd cannot certify, so they are opt-in. See its docstring.
        live = config.get('spnlp_assumelive') if isinstance(config, dict) \
            else getattr(config, 'spnlp_assumelive', None)
        if live is not None:
            lpopt['assumelive'] = bool(live)
        init = config.get('spnlp_init') if isinstance(config, dict) \
            else getattr(config, 'spnlp_init', None)
        if init is not None:
            lpopt['init'] = init

    bnd = spn_lpbnd(sn, lpopt)

    Rn = int(bnd['nclasses'])
    for pp, nd in enumerate(bnd['places']):
        ist = int(nodetostation[nd])
        if ist < 0:
            continue
        for k in range(Rn):
            l = pp * Rn + k
            Q[ist, k] = bnd['tokens'][side, l]
            U[ist, k] = Q[ist, k]
            T[ist, k] = bnd['placeTput'][side, l]
            if T[ist, k] > 0:
                R[ist, k] = Q[ist, k] / T[ist, k]

    refstat = np.ravel(np.asarray(sn.refstat, dtype=int))
    for k in range(K):
        ref = int(refstat[k])
        if 0 <= ref < M:
            X[0, k] = T[ref, k]
        Nk = float(Q[:, k].sum())
        if X[0, k] > 0 and Nk > 0:
            C[0, k] = Nk / X[0, k]

    return Q, U, R, T, C, X


__all__ = ['BA_SPNLP', 'spnlp_bound']
