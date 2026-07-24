"""
Maximum Entropy Method (MEM) for Open Queueing Networks.

Parameter-extraction wrapper around :func:`me_oqn` (Kouvatsos 1994). Mirrors
MATLAB ``matlab/src/solvers/NC/solver_nc_mem.m``: it pulls service/arrival
rates, squared coefficients of variation, server counts and per-class routing
probabilities out of the network structure, then runs the
entropy-maximisation iteration.

Reference:
    D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
    Annals of Operations Research, 48:63-126, 1994.
"""
from typing import Any, Optional, Tuple

import numpy as np

from .me_oqn import me_oqn
from .me_oqn_blk import me_oqn_blk, RULE_BAS, RULE_LOSS
from ..sn.predicates import sn_is_open_model, sn_is_closed_model
from ..sn.network_struct import NodeType
from ...lang.base import SchedStrategy, DropStrategy


def _config_get(options: Any, name: str, default: Any) -> Any:
    """Read ``options.config.<name>`` tolerating dict- or object-style options."""
    if options is None:
        return default
    cfg = options.get('config') if isinstance(options, dict) else getattr(options, 'config', None)
    if cfg is None:
        return default
    if isinstance(cfg, dict):
        return cfg.get(name, default)
    return getattr(cfg, name, default)


_SUPPORTED_SCHED = (SchedStrategy.EXT, SchedStrategy.INF, SchedStrategy.FCFS,
                    SchedStrategy.PS, SchedStrategy.SIRO, SchedStrategy.LCFS,
                    SchedStrategy.LCFSPR)

_SUPPORTED_NODES = (int(NodeType.SOURCE), int(NodeType.SINK),
                    int(NodeType.QUEUE), int(NodeType.DELAY))


def solver_nc_mem_supports(sn) -> Tuple[bool, str]:
    """Check whether the Maximum Entropy Method supports the model in ``sn``.

    Returns ``(True, '')`` when the model is either a plain open queueing
    network (Section 3.2: Source, Queue, Delay and Sink nodes; GE/GE/1,
    GE/GE/c and GE/GE/inf building blocks) or a plain closed queueing
    network (Section 3.3: Queue and Delay nodes; G/G/1 and G/G/inf
    building blocks only, so finite multiserver stations are rejected),
    in both cases without class switching and with non-priority scheduling
    disciplines only; otherwise ``(False, reason)`` with the first
    unsupported feature found.
    """
    isopen = sn_is_open_model(sn)
    isclosed = sn_is_closed_model(sn)
    ismixed = not isopen and not isclosed

    src_t, sink_t = int(NodeType.SOURCE), int(NodeType.SINK)
    for ind in range(int(sn.nnodes)):
        nt = int(sn.nodetype[ind])
        if nt not in _SUPPORTED_NODES:
            return False, "MEM supports only Source, Queue, Delay and Sink nodes."
        if isclosed and nt in (src_t, sink_t):
            return False, "MEM supports only Queue and Delay nodes in closed models."

    # Class switching is not part of the Kouvatsos (1994) network model
    R = int(sn.nclasses)
    csmask = np.asarray(sn.csmask, dtype=bool)
    if np.any(csmask & ~np.eye(R, dtype=bool)):
        return False, "MEM does not support class switching."

    # Absorbing self-loops (p_ii=1) make the routing reducible and the
    # geometric feedback transform 1/(1-p_ii) degenerate
    rtnodes = np.asarray(sn.rtnodes, dtype=float)
    st2node = np.asarray(sn.get_station_indices()).astype(int)
    for ist in range(int(sn.nstations)):
        ind = st2node[ist]
        for r in range(R):
            if rtnodes[ind * R + r, ind * R + r] >= 1 - 1e-9:
                return False, ("MEM does not support absorbing self-loop "
                               "routing (reducible network).")

    # Only non-priority disciplines are supported; the PR/HOL constraint
    # formulae are not given in Kouvatsos (1994)
    for ist in range(int(sn.nstations)):
        sched_ist = sn.sched[ist]
        if sched_ist not in _SUPPORTED_SCHED:
            return False, (f"MEM does not support the "
                           f"{SchedStrategy(sched_ist).name} scheduling strategy.")

    if isopen or ismixed:
        # A Source node must be present for the external arrival extraction
        if not any(int(sn.nodetype[ind]) == src_t for ind in range(int(sn.nnodes))):
            return False, "MEM requires a Source node when open classes are present."
    if isclosed or ismixed:
        # Closed classes build on G/G/1 and G/G/inf queues only (Section 3.3)
        nservers = np.asarray(sn.nservers, dtype=float).ravel()
        for ist in range(int(sn.nstations)):
            if np.isfinite(nservers[ist]) and nservers[ist] > 1:
                return False, "MEM does not support multiserver stations in closed or mixed models."

    # see _kb/03-api-layer.md for rationale
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    scv = np.asarray(sn.scv, dtype=float)
    capped = [ist for ist in range(int(sn.nstations))
              if sn.sched[ist] != SchedStrategy.EXT and np.isfinite(sn_get_buffer_size(sn, ist))]
    if capped:
        if not isopen:
            return False, "MEM supports finite station buffers only in open models."
        if R > 1:
            return False, ("MEM supports finite station buffers only in single-class models: "
                           "the censored GE/GE/c/0;N building block is single class.")
        for ist in capped:
            if not np.isfinite(nservers[ist]) or nservers[ist] < 1:
                return False, ("MEM cannot apply a finite buffer to the infinite-server "
                               f"station {ist + 1}.")
            if sn.sched[ist] != SchedStrategy.FCFS:
                return False, ("MEM supports finite station buffers only under FCFS scheduling; "
                               f"station {ist + 1} uses {SchedStrategy(sn.sched[ist]).name}.")
            dr = sn_get_drop_rule(sn, ist)
            if dr not in (DropStrategy.DROP, DropStrategy.BAS):
                return False, ("MEM supports the DROP and BAS drop rules at a finite buffer; "
                               f"station {ist + 1} uses {DropStrategy(dr).name}.")
            if np.isfinite(scv[ist, 0]) and scv[ist, 0] < 1 - 1e-12:
                return False, ("MEM with finite buffers needs a service scv of at least 1 at "
                               f"station {ist + 1}: the GE distribution is not defined below 1.")
        for ist in range(int(sn.nstations)):
            if sn.sched[ist] == SchedStrategy.EXT and np.isfinite(scv[ist, 0]) and scv[ist, 0] < 1 - 1e-12:
                return False, ("MEM with finite buffers needs an external interarrival scv of at "
                               "least 1: the GE distribution is not defined below 1.")
        return True, ""

    return True, ""


def sn_get_buffer_size(sn, ist) -> float:
    """Physical buffer size of a station in jobs, in service included.

    Kendall's K: the tighter of the station capacity sn.cap and the sum of
    the per-class capacities sn.classcap, inf when the station is unbounded.
    Both fields already fold setCapacity, setClassCapacity, a finite orbit
    and the closed-chain population, so this is the single place that decides
    whether a buffer BINDS.

    Only a buffer that can actually BIND is reported. refreshCapacity derives
    a FINITE classcap (the chain population) for EVERY closed model, so a
    plain finiteness test would report a buffer at every station of every
    closed model; a capacity at least as large as the total population can
    never refuse a job and is returned as inf. The population sum is inf as
    soon as one class is open, so any finite capacity reachable by an open
    class binds.
    """
    n = np.inf
    cap = getattr(sn, 'cap', None)
    if cap is not None:
        cap = np.asarray(cap, dtype=float).ravel()
        if ist < cap.size and cap[ist] >= 0:
            n = min(n, cap[ist])
    classcap = getattr(sn, 'classcap', None)
    if classcap is not None:
        classcap = np.asarray(classcap, dtype=float)
        if ist < classcap.shape[0]:
            row = classcap[ist, :]
            row = row[row > 0]
            if row.size:
                n = min(n, float(np.sum(row)))
    njobs = getattr(sn, 'njobs', None)
    if njobs is not None:
        total_jobs = float(np.sum(np.asarray(njobs, dtype=float)))
        if n >= total_jobs:
            n = np.inf  # declared but unreachable: it can never refuse a job
    return n


def sn_get_drop_rule(sn, ist) -> int:
    """Drop rule declared at a station for the first class.

    Defaults to DROP, which is what refreshCapacity assigns to a finite
    buffer reachable by an open class.
    """
    droprule = getattr(sn, 'droprule', None)
    if droprule is None:
        return int(DropStrategy.DROP)
    if isinstance(droprule, np.ndarray):
        if ist < droprule.shape[0] and droprule.shape[1] > 0:
            return int(droprule[ist, 0])
        return int(DropStrategy.DROP)
    try:
        per_class = droprule.get(ist)
        if per_class is None:
            return int(DropStrategy.DROP)
        if isinstance(per_class, dict):
            for key in (0, list(per_class.keys())[0] if per_class else None):
                if key in per_class:
                    return int(per_class[key])
            return int(DropStrategy.DROP)
        return int(per_class)
    except (AttributeError, IndexError, KeyError, TypeError):
        return int(DropStrategy.DROP)


def solver_nc_mem(sn, options: Optional[Any] = None
                  ) -> Tuple[np.ndarray, np.ndarray, np.ndarray,
                             np.ndarray, np.ndarray, np.ndarray, int]:
    """Maximum Entropy Method for open queueing networks.

    Args:
        sn: NetworkStruct from ``model.get_struct()``.
        options: Solver options; MEM-specific fields read from ``options.config``:
            ``mem_tol`` (default 1e-6), ``mem_maxiter`` (default 1000),
            ``mem_verbose`` (default False).

    Returns:
        ``(QN, UN, RN, TN, CN, XN, totiter)`` where QN/UN/RN/TN are ``M x R``
        arrays (queue lengths, utilizations, response times, throughputs) and
        CN/XN are ``1 x R`` system response times and throughputs per class.
    """
    memok, memreason = solver_nc_mem_supports(sn)
    if not memok:
        raise ValueError(memreason)

    M = int(sn.nstations)
    R = int(sn.nclasses)

    if sn_is_closed_model(sn):
        return _solver_nc_mem_closed(sn, options, M, R)
    if not sn_is_open_model(sn):
        return _solver_nc_mem_mixed(sn, options, M, R)

    tol = _config_get(options, 'mem_tol', 1e-6)
    maxiter = _config_get(options, 'mem_maxiter', 1000)
    verbose = _config_get(options, 'mem_verbose', False)
    mem_options = {'tol': tol, 'maxiter': maxiter, 'verbose': verbose}

    rates = np.asarray(sn.rates, dtype=float)
    scv = (np.asarray(sn.scv, dtype=float)
           if getattr(sn, 'scv', None) is not None else np.ones((M, R)))
    nservers_all = np.asarray(sn.nservers, dtype=float).ravel()

    # station index -> node index
    st2node = np.asarray(sn.get_station_indices()).astype(int)
    nodetype = sn.nodetype
    rtnodes = np.asarray(sn.rtnodes, dtype=float)
    src_t = int(NodeType.SOURCE)

    # Locate the source station (external arrivals)
    source_station = None
    for i in range(M):
        if int(nodetype[st2node[i]]) == src_t:
            source_station = i
            break
    if source_station is None:
        raise ValueError("MEM requires a Source node.")

    qs = [i for i in range(M) if i != source_station]  # queueing/delay stations
    Mq = len(qs)

    # Service rates, scvs and server counts (rates/scv are NaN for classes
    # not served at a station)
    mu = np.zeros((Mq, R))
    Cs = np.ones((Mq, R))
    nservers = np.ones(Mq)
    for k, ist in enumerate(qs):
        nservers[k] = nservers_all[ist]
        for r in range(R):
            if np.isfinite(rates[ist, r]) and rates[ist, r] > 0:
                mu[k, r] = rates[ist, r]
                if np.isfinite(scv[ist, r]) and scv[ist, r] > 0:
                    Cs[k, r] = scv[ist, r]

    # External arrivals: distribute the source output along its routing
    lambda0 = np.zeros((Mq, R))
    Ca0 = np.zeros((Mq, R))
    src_node = st2node[source_station]
    for r in range(R):
        if np.isfinite(rates[source_station, r]) and rates[source_station, r] > 0:
            ext_rate = rates[source_station, r]
            ca_ext = 1.0
            if np.isfinite(scv[source_station, r]) and scv[source_station, r] > 0:
                ca_ext = scv[source_station, r]
            for k, ist in enumerate(qs):
                dst_node = st2node[ist]
                route_prob = rtnodes[src_node * R + r, dst_node * R + r]
                if route_prob > 0:
                    lambda0[k, r] = ext_rate * route_prob
                    Ca0[k, r] = ca_ext

    # Routing probabilities between queueing stations
    P = np.zeros((Mq, Mq, R))
    for r in range(R):
        for j, jst in enumerate(qs):
            jn = st2node[jst]
            for k, ist in enumerate(qs):
                inode = st2node[ist]
                P[j, k, r] = rtnodes[jn * R + r, inode * R + r]

    # see _kb/03-api-layer.md for rationale
    nbuf = np.array([sn_get_buffer_size(sn, ist) for ist in qs], dtype=float)
    if np.any(np.isfinite(nbuf)):
        blockrule = np.array([RULE_BAS if sn_get_drop_rule(sn, ist) == int(DropStrategy.BAS)
                              else RULE_LOSS for ist in qs], dtype=int)
        Qb, Wb, Tb, Ub, _, _, pba, _, totiter = me_oqn_blk(
            Mq, lambda0[:, 0], np.where(Ca0[:, 0] > 0, Ca0[:, 0], 1.0), mu[:, 0], Cs[:, 0],
            P[:, :, 0], np.asarray(nservers, dtype=float).ravel(), nbuf, blockrule, mem_options)
        QN = np.zeros((M, R))
        UN = np.zeros((M, R))
        RN = np.zeros((M, R))
        TN = np.zeros((M, R))
        for k, ist in enumerate(qs):
            QN[ist, 0] = Qb[k]
            UN[ist, 0] = Ub[k]
            RN[ist, 0] = Wb[k]
            TN[ist, 0] = Tb[k]
        XN = np.zeros((1, R))
        CN = np.zeros((1, R))
        # see _kb/03-api-layer.md for rationale
        if np.isfinite(rates[source_station, 0]) and rates[source_station, 0] > 0:
            TN[source_station, 0] = rates[source_station, 0]
            XN[0, 0] = rates[source_station, 0]
            accepted = XN[0, 0]
            for k, ist in enumerate(qs):
                if lambda0[k, 0] > 0 and blockrule[k] == RULE_LOSS:
                    accepted -= lambda0[k, 0] * pba[k]
            if accepted > 0:
                CN[0, 0] = float(np.sum(QN[[ist for ist in qs], 0])) / accepted
        # Reported under its own name so that solver.citations() reaches the
        # transfer-blocking reference on top of the base MEM one.
        solver_nc_mem.last_method = 'mem.blocking'
        return QN, UN, RN, TN, CN, XN, totiter

    # Run the Maximum Entropy fixed-point algorithm
    insens = np.array([sn.sched[ist] in (SchedStrategy.PS, SchedStrategy.LCFSPR)
                       for ist in qs])
    L, W, Ca, Cd, lam, rho, totiter = me_oqn(
        Mq, R, lambda0, Ca0, mu, Cs, P, nservers, insens, mem_options)

    # Map results back to station-indexed LINE outputs
    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    RN = np.zeros((M, R))
    TN = np.zeros((M, R))
    for k, ist in enumerate(qs):
        QN[ist, :] = L[k, :]
        UN[ist, :] = rho[k, :]
        RN[ist, :] = W[k, :]
        TN[ist, :] = lam[k, :]

    # Cap utilization of unstable stations at 1 (LINE convention)
    for k, ist in enumerate(qs):
        if np.isfinite(nservers_all[ist]):
            utot = np.sum(UN[ist, :])
            if utot > 1:
                UN[ist, :] = UN[ist, :] / utot

    # Source station: report the external arrival rates as throughputs and
    # derive system metrics by Little's law
    XN = np.zeros((1, R))
    CN = np.zeros((1, R))
    for r in range(R):
        if np.isfinite(rates[source_station, r]) and rates[source_station, r] > 0:
            TN[source_station, r] = rates[source_station, r]
            XN[0, r] = rates[source_station, r]
            CN[0, r] = np.sum(QN[[ist for ist in qs], r]) / XN[0, r]

    return QN, UN, RN, TN, CN, XN, totiter


def _solver_nc_mem_closed(sn, options, M, R):
    """Closed-network MEM (Kouvatsos 1994, Section 3.3): two-stage
    pseudo-open decomposition plus convolution over the population lattice."""
    from .me_cqn import me_cqn

    tol = _config_get(options, 'mem_tol', 1e-6)
    maxiter = _config_get(options, 'mem_maxiter', 1000)
    verbose = _config_get(options, 'mem_verbose', False)
    mem_options = {'tol': tol, 'maxiter': maxiter, 'verbose': verbose}

    rates = np.asarray(sn.rates, dtype=float)
    scv = (np.asarray(sn.scv, dtype=float)
           if getattr(sn, 'scv', None) is not None else np.ones((M, R)))
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    refstat = np.asarray(sn.refstat, dtype=float).ravel().astype(int)

    st2node = np.asarray(sn.get_station_indices()).astype(int)
    rtnodes = np.asarray(sn.rtnodes, dtype=float)

    mu = np.zeros((M, R))
    Cs = np.ones((M, R))
    for ist in range(M):
        for r in range(R):
            if np.isfinite(rates[ist, r]) and rates[ist, r] > 0:
                mu[ist, r] = rates[ist, r]
                if np.isfinite(scv[ist, r]) and scv[ist, r] > 0:
                    Cs[ist, r] = scv[ist, r]

    # Routing probabilities between stations
    P = np.zeros((M, M, R))
    for r in range(R):
        for j in range(M):
            jn = st2node[j]
            for k in range(M):
                inode = st2node[k]
                P[j, k, r] = rtnodes[jn * R + r, inode * R + r]

    insens = np.array([sn.sched[ist] in (SchedStrategy.PS, SchedStrategy.LCFSPR)
                       for ist in range(M)])
    L, W, Ca, Cd, lam, rho, X, totiter = me_cqn(
        M, R, njobs, mu, Cs, P, nservers, refstat, insens, mem_options)

    QN = L
    UN = rho
    RN = W
    TN = lam
    XN = np.zeros((1, R))
    CN = np.zeros((1, R))
    for r in range(R):
        XN[0, r] = X[r]
        if X[r] > 0:
            CN[0, r] = njobs[r] / X[r]  # class cycle time by Little's law

    return QN, UN, RN, TN, CN, XN, totiter


def _solver_nc_mem_mixed(sn, options, M, R):
    """Mixed-network MEM: composition of the open (Section 3.2) and closed
    (Section 3.3) algorithms with product-form-style conditioning."""
    from .me_mqn import me_mqn

    tol = _config_get(options, 'mem_tol', 1e-6)
    maxiter = _config_get(options, 'mem_maxiter', 1000)
    verbose = _config_get(options, 'mem_verbose', False)
    mem_options = {'tol': tol, 'maxiter': maxiter, 'verbose': verbose}

    rates = np.asarray(sn.rates, dtype=float)
    scv = (np.asarray(sn.scv, dtype=float)
           if getattr(sn, 'scv', None) is not None else np.ones((M, R)))
    nservers_all = np.asarray(sn.nservers, dtype=float).ravel()
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    refstat_all = np.asarray(sn.refstat, dtype=float).ravel().astype(int)

    st2node = np.asarray(sn.get_station_indices()).astype(int)
    nodetype = sn.nodetype
    rtnodes = np.asarray(sn.rtnodes, dtype=float)
    src_t = int(NodeType.SOURCE)

    open_cls = np.array([np.isinf(njobs[r]) for r in range(R)])
    N = njobs.copy()

    source_station = None
    for i in range(M):
        if int(nodetype[st2node[i]]) == src_t:
            source_station = i
            break
    if source_station is None:
        raise ValueError("MEM requires a Source node when open classes are present.")

    qs = [i for i in range(M) if i != source_station]
    Mq = len(qs)

    mu = np.zeros((Mq, R))
    Cs = np.ones((Mq, R))
    nservers = np.ones(Mq)
    refstat = np.zeros(R, dtype=int)
    for k, ist in enumerate(qs):
        nservers[k] = nservers_all[ist]
        for r in range(R):
            if np.isfinite(rates[ist, r]) and rates[ist, r] > 0:
                mu[k, r] = rates[ist, r]
                if np.isfinite(scv[ist, r]) and scv[ist, r] > 0:
                    Cs[k, r] = scv[ist, r]
    for r in range(R):
        if not open_cls[r]:
            refstat[r] = qs.index(refstat_all[r])

    # External arrivals of the open classes along the source routing
    lambda0 = np.zeros((Mq, R))
    Ca0 = np.zeros((Mq, R))
    src_node = st2node[source_station]
    for r in range(R):
        if open_cls[r] and np.isfinite(rates[source_station, r]) and rates[source_station, r] > 0:
            ext_rate = rates[source_station, r]
            ca_ext = 1.0
            if np.isfinite(scv[source_station, r]) and scv[source_station, r] > 0:
                ca_ext = scv[source_station, r]
            for k, ist in enumerate(qs):
                dst_node = st2node[ist]
                route_prob = rtnodes[src_node * R + r, dst_node * R + r]
                if route_prob > 0:
                    lambda0[k, r] = ext_rate * route_prob
                    Ca0[k, r] = ca_ext

    # Routing probabilities between queueing stations
    P = np.zeros((Mq, Mq, R))
    for r in range(R):
        for j, jst in enumerate(qs):
            jn = st2node[jst]
            for k, ist in enumerate(qs):
                inode = st2node[ist]
                P[j, k, r] = rtnodes[jn * R + r, inode * R + r]

    insens = np.array([sn.sched[ist] in (SchedStrategy.PS, SchedStrategy.LCFSPR)
                       for ist in qs])
    L, W, Ca, Cd, lam, rho, X, totiter = me_mqn(
        Mq, R, open_cls, lambda0, Ca0, N, mu, Cs, P, nservers, refstat, insens, mem_options)

    QN = np.zeros((M, R))
    UN = np.zeros((M, R))
    RN = np.zeros((M, R))
    TN = np.zeros((M, R))
    for k, ist in enumerate(qs):
        QN[ist, :] = L[k, :]
        UN[ist, :] = rho[k, :]
        RN[ist, :] = W[k, :]
        TN[ist, :] = lam[k, :]

    # Cap utilization of unstable stations at 1 (LINE convention)
    for k, ist in enumerate(qs):
        if np.isfinite(nservers_all[ist]):
            utot = np.sum(UN[ist, :])
            if utot > 1:
                UN[ist, :] = UN[ist, :] / utot

    XN = np.zeros((1, R))
    CN = np.zeros((1, R))
    for r in range(R):
        XN[0, r] = X[r]
        if open_cls[r]:
            TN[source_station, r] = X[r]
            if X[r] > 0:
                CN[0, r] = np.sum(QN[[ist for ist in qs], r]) / X[r]
        else:
            if X[r] > 0:
                CN[0, r] = N[r] / X[r]  # class cycle time

    return QN, UN, RN, TN, CN, XN, totiter
