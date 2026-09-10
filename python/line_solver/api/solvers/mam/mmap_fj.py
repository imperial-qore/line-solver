"""
MAM/MMAP fork-join decomposition ('dec.source.mmap').

Native Python port of the MATLAB family

    matlab/src/solvers/MAM/solver_mam_basic_mmap.m
    matlab/src/solvers/MAM/solver_mam_basic_mmap_inner.m
    matlab/src/solvers/MAM/solver_mam_basic_mmap_closed.m
    matlab/src/solvers/MAM/solver_mam_traffic_mmap.m
    matlab/src/solvers/MAM/mam_truncate_renorm.m
    matlab/src/solvers/MAM/mam_detect_mmck.m

mirrored by the JAR handlers jline.solvers.mam.handlers.Solver_mam_basic_mmap*.

The method refines the per-node departure processes of the network by a
parametric decomposition, superposing and splitting them as marked MAPs and
synchronizing the branches of a fork-join construct with mmap_max. Open
networks call the inner algorithm directly with the arrival rates implied by
the source; closed networks wrap it in a per-class bisection on a surrogate
arrival rate so that the queue lengths reproduce the closed population.

Representation conventions used throughout this module:
  - a MAP is a list [D0, D1] of square ndarrays,
  - an MMAP is a list [D0, D1, Dc1, ..., DcK], i.e. the aggregate D1 at index 1
    and the per-class marking matrices from index 2 on (the MATLAB cell layout),
  - DEP[ind][r] is the class-r departure MAP of node ind,
  - ARV[ind] is the arrival MMAP of node ind, or None when node ind takes no
    arrivals (a Source, or a ClassSwitch hidden by the stochastic complement).
"""

import numpy as np
import numpy.matlib as ml
from dataclasses import replace
from typing import Any, Dict, List, Optional, Tuple

from ...sn import NetworkStruct, NodeType, SchedStrategy, sn_is_open_model
from ...da import da_fpi
from ...mc.dtmc import dtmc_stochcomp
from ...mam.map_analysis import map_acf, map_mean, map_pie
from ...mam.mmap_ops import (
    mmap_compress,
    mmap_hide,
    mmap_lambda,
    mmap_normalize,
    mmap_max,
    mmap_super,
)
from ...mam.qbd_depproc import qbd_depproc_etaqa, qbd_depproc_etaqa_ps
from ...npfqn.traffic import npfqn_traffic_merge, npfqn_traffic_split_cs
from ...qsys import qsys_mmck, qsys_mmapg1k
from ....constants import GlobalConstants
from ....lib.thirdparty.butools.queues import MMAPPH1FCFS
from ....lib.thirdparty.fj import sn_build_fj_sync_map
from .handler import (
    SolverMAMOptions,
    SolverMAMReturn,
    _get_scheduling,
    _get_service_times,
    _get_visits,
    _mam_detect_mmck,
)

__all__ = [
    'solver_mam_basic_mmap',
    'solver_mam_basic_mmap_inner',
    'solver_mam_basic_mmap_closed',
    'solver_mam_traffic_mmap',
    'mam_truncate_renorm',
    'mam_svc_mixture',
    'mam_detect_mmck',
]

# Scheduling disciplines the inner algorithm treats as FCFS-like, i.e. served in
# arrival order by a single queue whose departure process is built from the
# MAP/MAP/1 QBD (MATLAB: the {FCFS, HOL, FCFSPRPRIO} case labels).
_FCFS_LIKE = (SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.FCFSPRPRIO)

# Node types that emit a departure flow into the traffic equations.
_FLOW_SOURCES = (NodeType.SOURCE, NodeType.DELAY, NodeType.QUEUE,
                 NodeType.FORK, NodeType.JOIN)


def _zero_mmap(K: int) -> List[np.ndarray]:
    """MATLAB {[0],[0],[0]}: the null flow carrying no arrivals of any class."""
    return [np.zeros((1, 1)) for _ in range(2 + K)]


def _map_normalize_feasible(D0: np.ndarray, D1: np.ndarray) -> List[np.ndarray]:
    """Make a (D0,D1) pair a feasible MAP, as MATLAB map_normalize does.

    Negative entries are clamped to zero and the diagonal of D0 is rebuilt so
    that D0+D1 closes as an infinitesimal generator. This is NOT the native
    map_analysis.map_normalize, which rescales a MAP to unit mean; the MATLAB
    function of the same name repairs feasibility and leaves the mean alone.
    """
    D0 = np.asarray(D0, dtype=np.float64).copy()
    D1 = np.asarray(D1, dtype=np.float64).copy()
    D0 = np.real(D0)
    D1 = np.real(D1)
    D0[D0 < 0] = 0.0
    D1[D1 < 0] = 0.0
    n = D0.shape[0]
    for i in range(n):
        D0[i, i] = 0.0
        D0[i, i] = -(np.sum(D0[i, :]) + np.sum(D1[i, :]))
    return [D0, D1]


def _map_scale(mapproc: List[np.ndarray], newmean: float) -> List[np.ndarray]:
    """Rescale a MAP to have mean inter-arrival time newmean (MATLAB map_scale).

    The native map_analysis.map_scale multiplies the generator by a factor
    instead of targeting a mean, so it cannot be used directly here.
    """
    D0 = np.asarray(mapproc[0], dtype=np.float64)
    D1 = np.asarray(mapproc[1], dtype=np.float64)
    ratio = map_mean(D0, D1) / newmean
    return _map_normalize_feasible(D0 * ratio, D1 * ratio)


def _map_exponential(mean: float) -> List[np.ndarray]:
    """Poisson process of mean inter-arrival time mean, as a MAP."""
    mu = 1.0 / mean
    return [np.array([[-mu]]), np.array([[mu]])]


def _proc_to_map(entry: Any) -> Optional[List[np.ndarray]]:
    """Normalize an sn.proc entry into a MAP [D0, D1].

    Native models store a general MAP/PH as its (D0,D1) matrices but keep
    Exp/Erlang/HyperExp in compact dict form, whereas MATLAB's sn.proc always
    holds (D0,D1). Dict forms expand to the renewal MAP D0 = T, D1 = (-T e) a.
    Returns None for an absent process (a class disabled at the station).
    """
    if entry is None:
        return None
    # sn.proc stores (D0, D1); proc_to_map also accepts the legacy descriptors.
    from ...sn.proc_form import proc_to_map
    D0, D1 = proc_to_map(entry)
    if D0 is not None:
        return [np.atleast_2d(D0), np.atleast_2d(D1)]
    raise RuntimeError("solver_mam_basic_mmap: unsupported process representation %r"
                       % type(entry))


def _mmap_normalize_list(mmap: List[np.ndarray]) -> List[np.ndarray]:
    """mmap_normalize on the [D0, D1, Dc1..DcK] list layout."""
    D0, D_list = mmap_normalize(mmap[0], list(mmap[1:]))
    return [D0] + D_list


def _mmap_compress_list(mmap: List[np.ndarray]) -> List[np.ndarray]:
    """Compress an MMAP down to order 1.

    All mmap_compress calls in the MATLAB sources reach the
    'default'/'mixture'/'mixture.order1' branch: solver_mam_traffic_mmap and
    solver_mam_basic_mmap_inner pass options.config, which carries no 'method'
    field and so defaults to 'default', while the post-mmap_max call passes
    method=config.compress='mixture.order1'. The native mmap_compress spells the
    same order-1 mixture fit 'mixture'.
    """
    D0, D_list = mmap_compress(mmap[0], list(mmap[1:]), method="mixture")
    return [D0] + D_list


def _mmap_order(mmap: List[np.ndarray]) -> int:
    """Number of phases of an MMAP (MATLAB length(MMAP{1}))."""
    return np.asarray(mmap[0]).shape[0]


def mam_detect_mmck(sn: NetworkStruct, ist: int, K: int,
                    mmap_node: Optional[List[np.ndarray]]) -> Tuple[bool, float]:
    """Decide whether station ist matches the M/M/c/K assumptions.

    Port of matlab/src/solvers/MAM/mam_detect_mmck.m. Returns (True, muRate)
    only when the aggregated arrival MMAP at the node is single-phase (a Poisson
    superposition of the class arrivals) and every active class has exponential
    service at the station at one shared rate.

    Args:
        sn: network struct
        ist: station index
        K: number of classes
        mmap_node: aggregated arrival MMAP [D0, D1, Dc1..DcK] at the node

    Returns:
        (is_mmck, mu_rate); mu_rate is NaN when is_mmck is False.
    """
    if mmap_node is None or len(mmap_node) == 0:
        return False, float('nan')
    if np.asarray(mmap_node[0]).shape[0] != 1:
        # arrivals are not a single-phase Poisson superposition
        return False, float('nan')
    return _mam_detect_mmck(sn, ist, K)


def mam_svc_mixture(D_arr: List[np.ndarray],
                    pie_list: List[np.ndarray],
                    D0_list: List[np.ndarray]) -> Dict[str, Any]:
    """Arrival-weighted phase-type mixture of the per-class service laws.

    Port of matlab/src/solvers/MAM/mam_svc_mixture.m. Returns the service
    descriptor accepted by qsys_mapg1k/qsys_mmapg1k.

    The mixture is PH(alpha_mix, T_mix) with alpha_mix = [w_1*pie_1, ...],
    T_mix = blkdiag(D0_1, ...), and w_k = lambda_k/sum_j lambda_j the fraction
    of arrivals belonging to class k. It is therefore the service law of an
    arbitrary packet, and it reduces to the common law exactly (as a
    distribution) when every class shares one.

    Args:
        D_arr: [D0, D_class1, ..., D_classR], the arrival MMAP
        pie_list: per-class PH initial distributions
        D0_list: per-class PH sub-generators

    Returns:
        Dict service descriptor with keys 'type' ('ph'), 'alpha', 'T'.
    """
    from ...mc.ctmc import ctmc_solve
    n_classes = len(pie_list)
    D0 = np.asarray(D_arr[0], dtype=np.float64)
    # 1-D on BOTH sides: with theta a (1,M) row and e a (M,1) column the product
    # below is a (1,1) array, and numpy 2 no longer coerces a size-1 array to a
    # scalar, so every float() here raised TypeError instead of computing.
    e_arr = np.ones(D0.shape[0])
    Dsum = np.zeros_like(D0)
    for k in range(n_classes):
        Dsum = Dsum + np.asarray(D_arr[k + 1], dtype=np.float64)
    theta = np.asarray(ctmc_solve(D0 + Dsum), dtype=np.float64).reshape(-1)
    lambda_k = np.zeros(n_classes)
    for k in range(n_classes):
        lambda_k[k] = float(theta @ np.asarray(D_arr[k + 1], dtype=np.float64) @ e_arr)
    sumL = float(np.sum(lambda_k))
    if sumL > 0:
        w = lambda_k / sumL
    else:
        w = np.ones(n_classes) / n_classes

    n_k = [np.asarray(p).size for p in pie_list]
    n_total = int(np.sum(n_k))
    alpha_mix = np.zeros(n_total)
    T_mix = np.zeros((n_total, n_total))
    offset = 0
    for k in range(n_classes):
        nk = n_k[k]
        alpha_mix[offset:offset + nk] = w[k] * np.asarray(pie_list[k]).flatten()
        T_mix[offset:offset + nk, offset:offset + nk] = np.asarray(D0_list[k])
        offset += nk
    return {'type': 'ph', 'alpha': alpha_mix, 'T': T_mix}


def mam_truncate_renorm(D_arr: List[np.ndarray],
                        pie_list: List[np.ndarray],
                        D0_list: List[np.ndarray],
                        capK: int) -> Tuple[float, float, np.ndarray]:
    """Finite-buffer marginal for MMAP[K]/PH[K]/1/FCFS.

    Port of matlab/src/solvers/MAM/mam_truncate_renorm.m. Solves the
    infinite-buffer queue via BuTools MMAPPH1FCFS, truncates the marginal queue
    length distribution at the buffer capacity capK and renormalizes. For an
    M/M/1 input the renormalized distribution coincides exactly with the M/M/1/K
    marginal; for a general MMAP/PH it is an ASTA-style approximation.

    Args:
        D_arr: [D0, D_class1, ..., D_classK], the MMAP argument of MMAPPH1FCFS
        pie_list: per-class PH initial distributions
        D0_list: per-class PH sub-generators
        capK: buffer capacity (max number of jobs in system, >= 1)

    Returns:
        (meanQ, lossProb, p_norm) with meanQ the mean number in system clipped
        to [0, capK], lossProb the renormalized boundary mass p(N=capK) and
        p_norm the renormalized truncated marginal of length capK+1.
    """
    n_classes = len(pie_list)

    if n_classes > 1:
        # MMAPPH1FCFS('ncDistr') returns the per-class marginal P(N_k=n), not the
        # joint P(N_total=n) the truncation requires. Aggregate to a single-class
        # MMAP/PH/1 by summing the arrival matrices and building a
        # workload-weighted PH mixture for service.
        from ...mc.ctmc import ctmc_solve
        D0 = np.asarray(D_arr[0], dtype=np.float64)
        Dsum = np.zeros_like(D0)
        for k in range(n_classes):
            Dsum = Dsum + np.asarray(D_arr[k + 1], dtype=np.float64)
        # 1-D on both sides; see mam_svc_mixture for why the (1,1) product broke.
        e_arr = np.ones(D0.shape[0])
        theta = np.asarray(ctmc_solve(D0 + Dsum), dtype=np.float64).reshape(-1)
        lambda_k = np.zeros(n_classes)
        for k in range(n_classes):
            lambda_k[k] = float(theta @ np.asarray(D_arr[k + 1], dtype=np.float64) @ e_arr)
        sumL = float(np.sum(lambda_k))
        if sumL > 0:
            w = lambda_k / sumL
        else:
            w = np.ones(n_classes) / n_classes
        n_k = [np.asarray(p).size for p in pie_list]
        n_total = int(np.sum(n_k))
        alpha_mix = np.zeros((1, n_total))
        T_mix = np.zeros((n_total, n_total))
        offset = 0
        for k in range(n_classes):
            nk = n_k[k]
            alpha_mix[0, offset:offset + nk] = w[k] * np.asarray(pie_list[k]).flatten()
            T_mix[offset:offset + nk, offset:offset + nk] = np.asarray(D0_list[k])
            offset += nk
        D_call = [D0, Dsum]
        pie_call = [alpha_mix]
        D0_call = [T_mix]
    else:
        D_call = list(D_arr)
        pie_call = list(pie_list)
        D0_call = list(D0_list)

    n_levels = capK + 1
    pdistr = MMAPPH1FCFS([ml.matrix(np.asarray(D, dtype=np.float64)) for D in D_call],
                         [ml.matrix(np.asarray(p, dtype=np.float64).reshape(1, -1)) for p in pie_call],
                         [ml.matrix(np.asarray(T, dtype=np.float64)) for T in D0_call],
                         'ncDistr', n_levels)
    pdistr = np.abs(np.asarray(pdistr, dtype=np.float64).flatten())
    if pdistr.size < n_levels:
        pdistr = np.concatenate([pdistr, np.zeros(n_levels - pdistr.size)])
    p_in = pdistr[:n_levels]
    mass_in = float(np.sum(p_in))
    if mass_in <= 0:
        p_norm = np.zeros(n_levels)
        p_norm[0] = 1.0
    else:
        p_norm = p_in / mass_in
    mean_q = max(0.0, min(float(capK), float(np.arange(n_levels) @ p_norm)))
    loss_prob = float(p_norm[-1])
    return mean_q, loss_prob, p_norm


def solver_mam_traffic_mmap(sn: NetworkStruct,
                            DEP: Dict[int, Dict[int, List[np.ndarray]]],
                            config: SolverMAMOptions,
                            fj_sync_map) -> Dict[int, Optional[List[np.ndarray]]]:
    """FJ-aware traffic solver: superpose, split and synchronize the flows.

    Port of matlab/src/solvers/MAM/solver_mam_traffic_mmap.m. Extends the plain
    traffic solver with mmap_max synchronization at the join points.

    Args:
        sn: network struct
        DEP: DEP[ind][r], the class-r departure MAP [D0, D1] of node ind
        config: solver options carrying merge/compress/space_max/fj_sync_q_len
        fj_sync_map: map built by sn_build_fj_sync_map

    Returns:
        ARV[ind], the arrival MMAP at node ind, or None when node ind takes no
        arrival flow.
    """
    I = sn.nnodes
    R = sn.nclasses
    fj_sync_q_len = int(getattr(config, 'fj_sync_q_len', 2) or 2)
    space_max = config.space_max

    # Index over all non-ClassSwitch nodes
    non_cs_classes: List[int] = []
    is_ncs = [False] * I
    node_to_ncs = [-1] * I
    ncs_to_node: List[int] = []
    for ind in range(I):
        if sn.nodetype[ind] != NodeType.CLASSSWITCH:
            non_cs_classes.extend(range(ind * R, (ind + 1) * R))
            is_ncs[ind] = True
            node_to_ncs[ind] = len(ncs_to_node)
            ncs_to_node.append(ind)

    # Hide the nodes that are class switches
    rtncs = np.asarray(dtmc_stochcomp(np.asarray(sn.rtnodes, dtype=np.float64),
                                      np.array(non_cs_classes, dtype=int)),
                       dtype=np.float64)
    Inc = len(ncs_to_node)

    # DEP is indexed by node (ind, r): convert to the MMAP layout, marking the
    # single class of each departure flow (MATLAB MMAP{ind,r}{3} = MMAP{ind,r}{2})
    MMAP: Dict[int, Dict[int, List[np.ndarray]]] = {}
    for ind in range(I):
        row: Dict[int, List[np.ndarray]] = {}
        for r in range(R):
            d = DEP.get(ind, {}).get(r)
            if d is None or len(d) == 0 or np.any(np.isnan(np.asarray(d[0], dtype=np.float64))):
                row[r] = [np.zeros((1, 1)), np.zeros((1, 1)), np.zeros((1, 1))]
            else:
                D0 = np.asarray(d[0], dtype=np.float64)
                D1 = np.asarray(d[1], dtype=np.float64)
                row[r] = [D0, D1, D1.copy()]
        MMAP[ind] = row

    ARV: Dict[int, Optional[List[np.ndarray]]] = {}
    DEP_NCS: Dict[int, List[np.ndarray]] = {}
    LINKS: Dict[int, Dict[int, Optional[List[np.ndarray]]]] = {}

    # Build the nodeSync matrix in NCS indexing
    node_sync_ncs = np.zeros((Inc, Inc), dtype=int)
    for ind in range(I):
        if not is_ncs[ind]:
            continue
        inc = node_to_ncs[ind]
        for jnd in range(I):
            if not is_ncs[jnd]:
                continue
            jnc = node_to_ncs[jnd]
            g = int(fj_sync_map.nodeSync[ind, jnd])
            if g > 0:
                node_sync_ncs[inc, jnc] = g

    # First determine all outgoing flows from all nodes
    for ind in range(I):
        if not is_ncs[ind] or sn.nodetype[ind] not in _FLOW_SOURCES:
            continue
        inc = node_to_ncs[ind]

        if R > 1:
            # Order-preserving bounded superposition of the per-class departure
            # flows. Superposing all R at once builds the full Kronecker product
            # (order = product of the per-class orders), which compounds
            # geometrically across the refinement iterations and exhausts memory.
            # Superposing class by class in class order keeps the marked-class
            # order 1..R that the downstream class-switch split relies on, and is
            # identical to the one-shot product whenever it stays within
            # space_max.
            dep = MMAP[ind][0]
            for rr in range(1, R):
                dep = mmap_super(dep, MMAP[ind][rr])
                if _mmap_order(dep) > space_max:
                    dep = _mmap_compress_list(dep)
        else:
            dep = MMAP[ind][0]
        DEP_NCS[inc] = dep

        Psplit = np.zeros((R, Inc * R))
        for r in range(R):
            for jnd in range(I):
                if not is_ncs[jnd]:
                    continue
                jnc = node_to_ncs[jnd]
                for s in range(R):
                    Psplit[r, jnc * R + s] = rtncs[inc * R + r, jnc * R + s]

        Fsplit = npfqn_traffic_split_cs(DEP_NCS[inc], Psplit)
        LINKS[inc] = {}
        for jnc in range(Inc):
            flow = Fsplit.get(jnc)
            LINKS[inc][jnc] = _mmap_normalize_list(flow) if flow is not None else None

    # Then determine all incoming flows, with FJ synchronization
    for ind in range(I):
        if not is_ncs[ind] or sn.nodetype[ind] == NodeType.SOURCE:
            ARV[ind] = None
            continue
        inc = node_to_ncs[ind]

        # Partition incoming links into sync groups and independent flows
        sync_groups_at_node = sorted({int(g) for g in node_sync_ncs[inc, :] if g > 0})

        independent_flows: List[List[np.ndarray]] = []
        sync_flows: Dict[int, List[List[np.ndarray]]] = {}

        for jnc in range(Inc):
            flow = LINKS.get(jnc, {}).get(inc)
            if flow is None or len(flow) == 0:
                continue
            if float(np.sum(mmap_lambda(flow))) <= GlobalConstants.FineTol:
                continue
            gid = int(node_sync_ncs[inc, jnc])
            if gid == 0:
                independent_flows.append(flow)
            else:
                sync_flows.setdefault(gid, []).append(flow)

        # Process synchronized flows: apply mmap_max iteratively within each group
        sync_results: List[List[np.ndarray]] = []
        for gid in sync_groups_at_node:
            group_flows = sync_flows.get(gid)
            if not group_flows:
                continue
            synced_flow = group_flows[0]
            for f in range(1, len(group_flows)):
                synced_flow = mmap_max(synced_flow, group_flows[f], fj_sync_q_len)
                # mmap_max builds the synchronization state space but does not
                # enforce MMAP feasibility. Normalize before using the result in
                # compression or lambda calculations.
                synced_flow = _mmap_normalize_list(synced_flow)
                if _mmap_order(synced_flow) > space_max:
                    synced_flow = _mmap_compress_list(synced_flow)
            sync_results.append(synced_flow)

        # Merge synced flows with independent flows
        all_flows = sync_results + independent_flows

        if len(all_flows) > 1:
            ARV[ind] = npfqn_traffic_merge({i: f for i, f in enumerate(all_flows)},
                                           config.merge, config.compress)
        elif len(all_flows) == 1:
            ARV[ind] = all_flows[0]
        else:
            # No flows: take the first non-empty link
            fallback = None
            for jnc in range(Inc):
                cand = LINKS.get(jnc, {}).get(inc)
                if cand is not None and len(cand) > 0:
                    fallback = cand
                    break
            ARV[ind] = fallback if fallback is not None else _zero_mmap(R)

    return ARV


def solver_mam_basic_mmap_inner(sn: NetworkStruct,
                                options: SolverMAMOptions,
                                lambda_r: np.ndarray) -> SolverMAMReturn:
    """MAM/MMAP fork-join decomposition at fixed per-class arrival rates.

    Port of matlab/src/solvers/MAM/solver_mam_basic_mmap_inner.m. Performs the
    departure-process refinement loop only; enforcing the closed population is
    the wrapper's responsibility (see solver_mam_basic_mmap_closed).

    Args:
        sn: network struct
        options: solver options
        lambda_r: per-class arrival rate vector of length nclasses

    Returns:
        SolverMAMReturn with the station metrics and the iteration count.
    """
    config = options
    etaqa_n = int(getattr(config, 'etaqa_trunc', 8) or 8)
    space_max = config.space_max

    I = sn.nnodes
    M = sn.nstations
    K = sn.nclasses
    V = _get_visits(sn)
    S = _get_service_times(sn)
    lambda_r = np.asarray(lambda_r, dtype=np.float64).flatten()

    nservers = np.asarray(sn.nservers, dtype=np.float64).flatten()
    rates = np.asarray(sn.rates, dtype=np.float64)
    njobs = np.asarray(sn.njobs, dtype=np.float64).flatten()
    cap = np.asarray(sn.cap, dtype=np.float64).flatten()
    isslc = np.asarray(sn.isslc, dtype=bool).flatten()
    node_to_station = np.asarray(sn.nodeToStation).flatten().astype(int)
    station_to_node = np.asarray(sn.stationToNode).flatten().astype(int)
    isstation = np.asarray(sn.isstation).flatten().astype(bool)

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    XN = np.zeros((1, K))

    fj_sync_map = sn_build_fj_sync_map(sn)

    # Local PH copy: MATLAB rescales PH{ist}{k} in place, so sn.proc must not be
    # mutated here.
    PH: Dict[int, Dict[int, Optional[List[np.ndarray]]]] = {}
    for ist in range(M):
        row: Dict[int, Optional[List[np.ndarray]]] = {}
        proc_ist = sn.proc[ist] if (getattr(sn, 'proc', None) is not None
                                    and ist < len(sn.proc)) else None
        for k in range(K):
            entry = proc_ist[k] if (proc_ist is not None and k < len(proc_ist)) else None
            row[k] = _proc_to_map(entry)
        PH[ist] = row
    pie: Dict[int, Dict[int, np.ndarray]] = {ist: {} for ist in range(M)}
    D0c: Dict[int, Dict[int, np.ndarray]] = {ist: {} for ist in range(M)}

    def prepare_pie(ist: int, k: int) -> None:
        p = PH[ist][k]
        if p is None:
            PH[ist][k] = _map_exponential(GlobalConstants.Immediate)
            pie[ist][k] = np.array([[1.0]])
            D0c[ist][k] = np.array([[-GlobalConstants.Immediate]])
            return
        pie[ist][k] = np.asarray(map_pie(p[0], p[1]), dtype=np.float64).reshape(1, -1)
        d0 = np.asarray(p[0], dtype=np.float64)
        if np.any(np.isnan(d0)):
            PH[ist][k] = _map_exponential(GlobalConstants.Immediate)
            pie[ist][k] = np.array([[1.0]])
            D0c[ist][k] = np.array([[-GlobalConstants.Immediate]])
        else:
            D0c[ist][k] = d0

    # Prepare PH service distributions
    for ist in range(M):
        sched = _get_scheduling(sn, ist)
        if sched == SchedStrategy.EXT:
            for k in range(K):
                r = rates[ist, k]
                TN[ist, k] = 0.0 if np.isnan(r) else r
        elif sched in _FCFS_LIKE or sched == SchedStrategy.PS:
            for k in range(K):
                p = PH[ist][k]
                if p is not None:
                    PH[ist][k] = _map_scale(p, map_mean(p[0], p[1]) / nservers[ist])
                prepare_pie(ist, k)
        elif sched == SchedStrategy.INF:
            for k in range(K):
                prepare_pie(ist, k)

    def ph_mean(ist: int, k: int) -> float:
        p = PH[ist][k]
        if p is None:
            return 0.0
        return float(map_mean(p[0], p[1]))

    def arrival_marks(arv: List[np.ndarray]) -> List[np.ndarray]:
        """MATLAB {ARV{ind}{[1,3:end]}}: D0 plus the per-class marks, dropping
        the aggregate D1 at position 2."""
        return [arv[0]] + [arv[2 + k] for k in range(K)]

    DEP: Dict[int, Dict[int, List[np.ndarray]]] = {}

    def init_dep() -> None:
        DEP.clear()
        for ind in range(I):
            row: Dict[int, List[np.ndarray]] = {}
            is_fork_join = sn.nodetype[ind] in (NodeType.FORK, NodeType.JOIN)
            if isstation[ind] and not is_fork_join:
                ist = node_to_station[ind]
                for r in range(K):
                    p = PH[ist][r]
                    if V[ist, r] > 0 and lambda_r[r] > 0:
                        row[r] = _map_scale(p, 1.0 / (lambda_r[r] * V[ist, r]))
                    elif isslc[r]:
                        # Self-looping classes are handled separately by the
                        # closed-network wrapper (their metrics are pinned at the
                        # reference station). They must not inject a saturating
                        # arrival stream into the shared-queue decomposition, so
                        # give them a vanishing-rate departure process.
                        row[r] = _map_exponential(1.0 / GlobalConstants.Zero)
                    else:
                        row[r] = p
            else:
                for r in range(K):
                    if lambda_r[r] > 0:
                        row[r] = _map_exponential(1.0 / lambda_r[r])
                    else:
                        row[r] = _map_exponential(1.0 / GlobalConstants.Immediate)
            DEP[ind] = row

    def update_dep(ARV: Dict[int, Optional[List[np.ndarray]]]) -> None:
        for ist in range(M):
            ind = station_to_node[ist]
            nt = sn.nodetype[ind]
            sched = _get_scheduling(sn, ist)
            if nt == NodeType.QUEUE:
                arv = ARV.get(ind)
                if arv is None or len(arv) == 0:
                    continue
                is_fcfs_like = sched in _FCFS_LIKE
                is_ps = (sched == SchedStrategy.PS)
                if not is_fcfs_like and not is_ps:
                    continue
                rho = float(np.sum(UN[ist, :]))
                for r in range(K):
                    types = [q for q in range(K) if q != r]
                    A0, A_list = mmap_hide(arv[0], list(arv[2:]), types)
                    A = [A0, sum(A_list)]
                    srv = PH[ist][r]
                    na = np.asarray(A[0]).shape[0]
                    ns = np.asarray(srv[0]).shape[0]
                    etaqa_sz = (etaqa_n + 1) * na * ns
                    if is_fcfs_like:
                        if etaqa_sz <= space_max and rho < 1 - GlobalConstants.FineTol:
                            try:
                                dep = qbd_depproc_etaqa(A, srv, etaqa_n)
                                dep = _map_normalize_feasible(dep[0], dep[1])
                            except Exception:
                                dep = srv
                        else:
                            dep = srv
                        if V[ist, r] > 0 and lambda_r[r] > 0:
                            dep = _map_scale(dep, 1.0 / (lambda_r[r] * V[ist, r]))
                        DEP[ind][r] = dep
                    else:
                        if V[ist, r] > 0 and lambda_r[r] > 0:
                            if etaqa_sz <= space_max and rho < 1 - GlobalConstants.FineTol:
                                try:
                                    dep = qbd_depproc_etaqa_ps(A, srv, etaqa_n)
                                    dep = _map_normalize_feasible(dep[0], dep[1])
                                except Exception:
                                    dep = srv
                            else:
                                dep = srv
                            dep = _map_scale(dep, 1.0 / (lambda_r[r] * V[ist, r]))
                            DEP[ind][r] = dep
            elif nt == NodeType.JOIN:
                for r in range(K):
                    if TN[ist, r] > 0:
                        DEP[ind][r] = _map_exponential(1.0 / TN[ist, r])

    def mmap_dec_sweep(_x, itnum):
        # Initialize departure processes (node-indexed: DEP[ind][r])
        if itnum == 1:
            init_dep()

        # Compute arrival processes with FJ synchronization
        ARV = solver_mam_traffic_mmap(sn, DEP, config, fj_sync_map)

        xref = QN.copy()
        for ist in range(M):
            ind = station_to_node[ist]
            nt = sn.nodetype[ind]
            sched = _get_scheduling(sn, ist)

            if nt == NodeType.JOIN:
                for k in range(K):
                    TN[ist, k] = lambda_r[k]
                    UN[ist, k] = 0.0
                    QN[ist, k] = 0.0
                    RN[ist, k] = 0.0
            elif nt == NodeType.QUEUE:
                arv = ARV.get(ind)
                if arv is None or len(arv) == 0:
                    continue
                if _mmap_order(arv) > space_max:
                    arv = _mmap_compress_list(arv)
                    ARV[ind] = arv

                finite_cap_used = False
                if sched in _FCFS_LIKE:
                    is_finite_cap = bool(np.isfinite(cap[ist]))
                    if is_finite_cap:
                        capK = int(cap[ist])
                        is_mmck, mu_mmck = mam_detect_mmck(sn, ist, K, arv)
                        loss_per_class_fc = None
                        if is_mmck:
                            lam = mmap_lambda(arv)
                            aggr_lambda = float(np.nansum(lam))
                            exact_res = qsys_mmck(aggr_lambda, mu_mmck,
                                                  int(nservers[ist]), capK)
                            mean_q_fc = exact_res['meanQueueLength']
                            loss_prob_fc = exact_res['lossProbability']
                        elif int(nservers[ist]) == 1:
                            # Exact MMAP[K]/G/1/K with per-class loss ratio; the
                            # phase resolution of pKvec tells apart classes of
                            # equal rate but different burstiness, which the
                            # truncate-and-renormalize fallback cannot express.
                            marks = arrival_marks(arv)
                            svc_mix = mam_svc_mixture(
                                marks,
                                [pie[ist][k] for k in range(K)],
                                [D0c[ist][k] for k in range(K)])
                            ex_res = qsys_mmapg1k(marks[0], marks[1:], svc_mix, capK)
                            mean_q_fc = ex_res['meanQueueLength']
                            loss_prob_fc = ex_res['lossAggregate']
                            loss_per_class_fc = np.asarray(ex_res['lossRatio'],
                                                           dtype=np.float64)
                        else:
                            mean_q_fc, loss_prob_fc, _ = mam_truncate_renorm(
                                arrival_marks(arv),
                                [pie[ist][k] for k in range(K)],
                                [D0c[ist][k] for k in range(K)],
                                capK)
                        lambda_inflow = np.asarray(mmap_lambda(arv), dtype=np.float64)
                        lambda_inflow = np.nan_to_num(lambda_inflow, nan=0.0)
                        if loss_per_class_fc is not None:
                            # exact branch: the loss ratio differs by class
                            TN_eff = lambda_inflow * (1 - loss_per_class_fc[:K])
                        else:
                            TN_eff = lambda_inflow * (1 - loss_prob_fc)
                        sum_tn = float(np.sum(TN_eff))
                        S_actual = np.array([ph_mean(ist, k) * nservers[ist]
                                             for k in range(K)])
                        if sum_tn > 0:
                            savg_eff = float(np.nansum(TN_eff * S_actual)) / sum_tn
                            Wq = max(0.0, mean_q_fc / sum_tn - savg_eff)
                        else:
                            Wq = 0.0
                        for k in range(K):
                            TN[ist, k] = TN_eff[k]
                            UN[ist, k] = TN[ist, k] * ph_mean(ist, k)
                            if TN[ist, k] > 0:
                                RN[ist, k] = Wq + S_actual[k]
                                QN[ist, k] = TN[ist, k] * RN[ist, k]
                            else:
                                RN[ist, k] = 0.0
                                QN[ist, k] = 0.0
                        finite_cap_used = True
                    else:
                        lam = np.asarray(mmap_lambda(arv), dtype=np.float64)
                        rho_classes = np.nan_to_num(
                            lam * np.array([ph_mean(ist, k) for k in range(K)]), nan=0.0)
                        # Exclude self-looping classes: they are pinned by the
                        # closed wrapper and must not drive the FCFS station into
                        # the saturation branch.
                        rho_ist = float(np.sum(rho_classes[~isslc]))
                        if rho_ist < 1 - GlobalConstants.FineTol:
                            # MMAPPH1FCFS models service as a renewal phase-type
                            # (marginal), discarding service-time autocorrelation.
                            # For the single-class single-server M/MAP/1 case with
                            # a genuinely correlated service process, use the exact
                            # MAP/MAP/1 queue, which carries the service phase
                            # across departures.
                            use_mapmap1 = (K == 1 and nservers[ist] == 1
                                           and abs(float(np.asarray(
                                               map_acf(PH[ist][0][0], PH[ist][0][1], 1)
                                           ).flatten()[0])) > GlobalConstants.CoarseTol)
                            if use_mapmap1:
                                from ....lib.thirdparty.qmam import (
                                    q_ct_map_map_1, MAPMAP1Options)
                                res_mm = q_ct_map_map_1(
                                    np.asarray(arv[0], dtype=np.float64),
                                    np.asarray(arv[2], dtype=np.float64),
                                    np.asarray(PH[ist][0][0], dtype=np.float64),
                                    np.asarray(PH[ist][0][1], dtype=np.float64),
                                    MAPMAP1Options(max_num_comp=100000))
                                ql = np.asarray(res_mm.queue_length, dtype=np.float64).flatten()
                                QN[ist, 0] = float(np.sum(np.arange(ql.size) * ql))
                            else:
                                qret = MMAPPH1FCFS(
                                    [ml.matrix(np.asarray(D, dtype=np.float64))
                                     for D in arrival_marks(arv)],
                                    [ml.matrix(np.asarray(pie[ist][k], dtype=np.float64).reshape(1, -1))
                                     for k in range(K)],
                                    [ml.matrix(np.asarray(D0c[ist][k], dtype=np.float64))
                                     for k in range(K)],
                                    'ncMoms', 1, 'ncDistr', 2)
                                # BuTools returns the measures class-major:
                                # [ncMoms_1, ncDistr_1, ncMoms_2, ncDistr_2, ...]
                                for k in range(K):
                                    QN[ist, k] = float(np.sum(np.asarray(qret[2 * k],
                                                                         dtype=np.float64)))
                        else:
                            # Saturation: bound queue lengths so the wrapper's
                            # bisection can recognise overload without NaNs
                            # propagating from MMAPPH1FCFS.
                            for k in range(K):
                                if np.isfinite(njobs[k]):
                                    QN[ist, k] = njobs[k]
                                else:
                                    QN[ist, k] = 1.0 / GlobalConstants.FineTol
                        for k in range(K):
                            TN[ist, k] = lam[k]
                elif sched == SchedStrategy.PS:
                    lam = np.asarray(mmap_lambda(arv), dtype=np.float64)
                    for k in range(K):
                        TN[ist, k] = lam[k]
                        UN[ist, k] = TN[ist, k] * S[ist, k]
                    # Self-looping classes are pinned by the closed wrapper and
                    # must not count toward the PS sharing denominator (a single
                    # permanent SLC job would otherwise drive the queue to
                    # saturation).
                    Uden = min(1 - GlobalConstants.FineTol, float(np.sum(UN[ist, ~isslc])))
                    for k in range(K):
                        QN[ist, k] = UN[ist, k] / (1 - Uden)

                if not finite_cap_used:
                    c = nservers[ist]
                    for k in range(K):
                        UN[ist, k] = TN[ist, k] * ph_mean(ist, k)
                        QN[ist, k] = QN[ist, k] + TN[ist, k] * (ph_mean(ist, k) * c) * (c - 1) / c
                        with np.errstate(divide='ignore', invalid='ignore'):
                            RN[ist, k] = QN[ist, k] / TN[ist, k]
            else:
                if sched == SchedStrategy.INF:
                    arv = ARV.get(ind)
                    if arv is not None and len(arv) > 0:
                        lam = np.asarray(mmap_lambda(arv), dtype=np.float64)
                        for k in range(K):
                            TN[ist, k] = lam[k]
                    for k in range(K):
                        if TN[ist, k] > 0:
                            UN[ist, k] = S[ist, k] * TN[ist, k]
                            QN[ist, k] = TN[ist, k] * S[ist, k]
                            RN[ist, k] = S[ist, k]
                # SchedStrategy.EXT: Source, TN already set above

        # Update departure processes
        update_dep(ARV)
        return QN.copy(), xref

    # departure-process fixed point (FJ parametric decomposition), driven on the
    # station queue lengths by the generic DA driver
    def _rel_norm(xn, xr):
        xn = np.asarray(xn, dtype=np.float64).flatten()
        xr = np.asarray(xr, dtype=np.float64).flatten()
        return float(np.max(np.abs(xn - xr) / (xr + GlobalConstants.FineTol)))

    # the legacy loop tested convergence only from the third sweep
    _, totiter, _ = da_fpi(mmap_dec_sweep, QN.copy(), options.iter_max,
                           options.iter_tol, norm=_rel_norm, miniter=3)
    if options.verbose:
        print("MAM FJ parametric decomposition completed in %d iterations." % totiter)

    # Join: derive QN/RN from parallel branch means
    for join_idx in range(I):
        if sn.nodetype[join_idx] != NodeType.JOIN:
            continue
        join_stat = node_to_station[join_idx]
        if join_stat < 0:
            continue
        sync_groups = sorted({int(g) for g in fj_sync_map.nodeSync[join_idx, :] if g > 0})
        for r in range(K):
            if TN[join_stat, r] <= 0:
                continue
            sync_delay = 0.0
            join_arrival_rate = 0.0
            for gid in sync_groups:
                branch_rt: List[float] = []
                branch_tput = 0.0
                for b in range(I):
                    if int(fj_sync_map.nodeSync[join_idx, b]) != gid:
                        continue
                    branch_stat = node_to_station[b]
                    if branch_stat < 0 or RN[branch_stat, r] <= 0:
                        continue
                    branch_rt.append(float(RN[branch_stat, r]))
                    branch_tput += float(TN[branch_stat, r])
                if len(branch_rt) < 2:
                    continue
                lambdai = np.array([1.0 / rt for rt in branch_rt])
                max_branch_rt = _inclusion_exclusion_max_mean(lambdai)
                sync_delay += max(max_branch_rt - float(np.mean(branch_rt)), 0.0)
                join_arrival_rate += branch_tput
            RN[join_stat, r] = sync_delay
            QN[join_stat, r] = join_arrival_rate * sync_delay
            UN[join_stat, r] = 0.0

    CN = np.nansum(RN, axis=0).reshape(1, -1)
    QN = np.nan_to_num(QN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)

    result = SolverMAMReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.method = 'dec.source.mmap'
    result.it = totiter
    return result


def _inclusion_exclusion_max_mean(lambdai: np.ndarray) -> float:
    """Mean of the maximum of independent exponentials of rates lambdai.

    MATLAB: sum_pow (-1)^pow * sum(1./sum(nchoosek(lambdai, pow+1), 2)), i.e.
    the inclusion-exclusion expansion E[max] = sum_{S nonempty} (-1)^{|S|+1} /
    sum_{i in S} lambda_i, enumerated here over the subset bitmasks.
    """
    n = len(lambdai)
    total = 0.0
    for mask in range(1, 1 << n):
        subset_sum = 0.0
        bits = 0
        for i in range(n):
            if mask & (1 << i):
                subset_sum += lambdai[i]
                bits += 1
        if subset_sum <= 0:
            continue
        term = 1.0 / subset_sum
        total += term if (bits % 2 == 1) else -term
    return total


def solver_mam_basic_mmap_closed(sn: NetworkStruct,
                                 options: SolverMAMOptions) -> SolverMAMReturn:
    """Closed-network wrapper around solver_mam_basic_mmap_inner.

    Port of matlab/src/solvers/MAM/solver_mam_basic_mmap_closed.m. Drives a
    per-class bisection on the surrogate arrival rate lambda so that the inner
    solver's queue lengths match the closed population sn.njobs. Mirrors the
    outer-loop structure of solver_mna_closed.
    """
    K = sn.nclasses
    M = sn.nstations
    rates = np.asarray(sn.rates, dtype=np.float64)
    nservers = np.asarray(sn.nservers, dtype=np.float64).flatten()
    njobs = np.asarray(sn.njobs, dtype=np.float64).flatten()
    isslc = np.asarray(sn.isslc, dtype=bool).flatten()
    refstat = np.asarray(sn.refstat).flatten().astype(int)
    S = _get_service_times(sn)

    # Per-class bisection bounds: upper = slowest non-INF station rate for that class
    non_inf_stations = np.where(np.isfinite(nservers))[0]
    inf_stations = np.where(~np.isfinite(nservers))[0]
    lambda_lb = np.zeros(K)
    lambda_ub = np.zeros(K)
    for k in range(K):
        rates_k = np.array([])
        if non_inf_stations.size > 0:
            rates_k = rates[non_inf_stations, k]
            rates_k = rates_k[np.isfinite(rates_k) & (rates_k > 0)]
        if rates_k.size == 0:
            rates_inf = rates[inf_stations, k] if inf_stations.size > 0 else np.array([])
            rates_inf = rates_inf[np.isfinite(rates_inf) & (rates_inf > 0)]
            lambda_ub[k] = 1.0 if rates_inf.size == 0 else float(np.max(rates_inf))
        else:
            lambda_ub[k] = float(np.min(rates_k))

    # open classes contribute 0; only closed populations gate convergence
    QNc = np.where(np.isfinite(njobs), njobs, 0.0)
    QN_chain = np.zeros(K)

    it_out = 0
    lambda_r = lambda_ub.copy()
    # Self-looping classes are pinned by the SLC clamp below; they must not
    # contribute a (saturating) surrogate arrival stream to the inner algorithm.
    lambda_r[isslc] = 0.0

    # MATLAB assigns a struct copy; mutating iter_max/verbose on the caller's
    # options object would leak the inner budget back out to the dispatcher.
    # Cap the MMAP phase-dimension truncation at the same value used by the
    # solver_mna_closed wrapper this routine mirrors. The analyzer default
    # (space_max=128) lets the FJ synchronization/superposition inflate the
    # arrival MMAP to ~128 phases, whose ETAQA/QBD solve is O(dim^3), while the
    # bisection re-enters the inner solve ~29 times: ~154s on a 3-class model
    # that CTMC solves in 0.7s. Measured, the compressed result is bit-identical
    # from space_max=4 up to 128, so the inflated dimension is pure wasted work.
    # Capping at 16 is lossless and matches solver_mna_closed.
    inner_options = replace(options,
                            iter_max=max(20, int(np.ceil(options.iter_max / 10.0))),
                            verbose=False,
                            space_max=min(options.space_max, 16))

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    # Last successful inner-algorithm outputs (fallback if the final trial diverges)
    QN_last, UN_last, RN_last = QN.copy(), UN.copy(), RN.copy()
    TN_last, CN_last, XN_last = TN.copy(), CN.copy(), XN.copy()
    have_good = False
    algorithm_ok = True

    bisect_tol = max(options.iter_tol, 1e-3)

    while np.max(np.abs(QN_chain - QNc)) > bisect_tol and it_out < options.iter_max:
        it_out += 1
        if it_out > 1:
            bracket_collapsed = True
            for k in range(K):
                if not np.isfinite(QNc[k]) or QNc[k] == 0 or isslc[k]:
                    continue
                if QN_chain[k] < QNc[k]:
                    lambda_lb[k] = lambda_r[k]
                else:
                    lambda_ub[k] = lambda_r[k]
                lambda_r[k] = 0.5 * (lambda_lb[k] + lambda_ub[k])
                # Bisection can still refine class k only while its bracket is
                # wider than the precision floor below which lambda cannot move
                # any reported metric.
                if (lambda_ub[k] - lambda_lb[k]) > GlobalConstants.FineTol * max(1.0, abs(lambda_ub[k])):
                    bracket_collapsed = False
            # Bracket-width stagnation break. The loop condition above tests only
            # the residual population gap, which never closes when no surrogate
            # lambda reproduces the closed population (the per-class targets are
            # not simultaneously attainable). The bisection then keeps halving
            # brackets that have already narrowed past FineTol, re-entering the
            # inner algorithm for changes it cannot resolve. Stop once every
            # active bracket has collapsed and keep the outputs in hand.
            if bracket_collapsed:
                it_out -= 1
                break

        try:
            r = solver_mam_basic_mmap_inner(sn, inner_options, lambda_r)
            QN, UN, RN, TN, CN, XN = r.Q, r.U, r.R, r.T, r.C, r.X
            algorithm_ok = True
        except Exception:
            # Inner algorithm diverged (typically MMAPPH1FCFS / lyap NaN under
            # saturation). Treat all chains as overloaded so the bisection drops
            # lambda on its next step.
            algorithm_ok = False

        if algorithm_ok:
            # SLC clamp: all jobs at refstat for self-looping classes
            for k in range(K):
                if isslc[k]:
                    QN[:, k] = 0.0
                    QN[refstat[k], k] = njobs[k]
            QN_chain = np.sum(QN, axis=0)
            QN_chain[~np.isfinite(QN_chain)] = 1.0 / GlobalConstants.FineTol
            QN_last, UN_last, RN_last = QN.copy(), UN.copy(), RN.copy()
            TN_last, CN_last, XN_last = TN.copy(), CN.copy(), XN.copy()
            have_good = True
        else:
            QN_chain = np.ones(K) * (1.0 / GlobalConstants.FineTol)

    # If the last trial diverged, fall back to the most recent successful one
    if not algorithm_ok and have_good:
        QN, UN, RN = QN_last, UN_last, RN_last
        TN, CN, XN = TN_last, CN_last, XN_last

    # Final SLC pass: pin throughput/utilisation at refstat (mirrors solver_mna_closed)
    for k in range(K):
        if isslc[k]:
            QN[:, k] = 0.0
            ist = refstat[k]
            QN[ist, k] = njobs[k]
            TN[ist, k] = njobs[k] * rates[ist, k]
            RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0.0
            UN[ist, k] = S[ist, k] * TN[ist, k]

    # Population redistribution within chain (matches solver_mna_closed)
    inchain_map = getattr(sn, 'inchain', None) or {}
    for c in range(sn.nchains):
        if c not in inchain_map:
            continue
        inchain = np.asarray(inchain_map[c]).flatten().astype(int)
        if inchain.size == 0:
            continue
        if c < njobs.size and np.isfinite(njobs[c]):
            sumQ = float(np.sum(QN[:, inchain]))
            if sumQ > 0:
                QN[:, inchain] = njobs[c] * QN[:, inchain] / sumQ

    # Delay/INF utilisation = mean number of jobs (matches solver_mna_closed)
    for ist in range(M):
        if _get_scheduling(sn, ist) == SchedStrategy.INF:
            UN[ist, :] = QN[ist, :]

    CN = np.nansum(RN, axis=0).reshape(1, -1)
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)

    result = SolverMAMReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.method = 'dec.source.mmap'
    result.it = it_out
    return result


def solver_mam_basic_mmap(sn: NetworkStruct,
                          options: Optional[SolverMAMOptions] = None) -> SolverMAMReturn:
    """Top-level dispatcher for the MAM/MMAP fork-join decomposition.

    Port of matlab/src/solvers/MAM/solver_mam_basic_mmap.m. Open networks call
    solver_mam_basic_mmap_inner directly with the arrival rates derived from the
    source/refstat; closed networks go through solver_mam_basic_mmap_closed,
    which wraps the inner algorithm in an MNA-style bisection on per-class
    throughput.
    """
    import time
    start_time = time.time()

    if options is None:
        options = SolverMAMOptions()

    if sn_is_open_model(sn):
        K = sn.nclasses
        C = sn.nchains
        rates = np.asarray(sn.rates, dtype=np.float64)
        refstat = np.asarray(sn.refstat).flatten().astype(int)
        inchain_map = getattr(sn, 'inchain', None) or {}
        lambda_r = np.zeros(K)
        for c in range(C):
            if c not in inchain_map:
                continue
            inchain = np.asarray(inchain_map[c]).flatten().astype(int)
            if inchain.size == 0:
                continue
            lambdas_inchain = rates[refstat[inchain[0]], inchain]
            lambdas_inchain = lambdas_inchain[np.isfinite(lambdas_inchain)]
            lambda_r[inchain] = float(np.sum(lambdas_inchain))
        result = solver_mam_basic_mmap_inner(sn, options, lambda_r)
    else:
        result = solver_mam_basic_mmap_closed(sn, options)

    result.A = result.T.copy() if result.T is not None else None
    result.W = result.R.copy() if result.R is not None else None
    result.runtime = time.time() - start_time
    return result
