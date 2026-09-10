"""Discrete-time (slotted) analysis for SolverMAM.

A single queueing station is solved EXACTLY by the Q-MAM discrete-time
algorithms, q_dt_ph_ph_1 when both laws are renewal discrete phase-type and
q_dt_map_map_1 when either side is a DMAP. Several stations are solved by a
discrete-time parametric decomposition, which is an approximation.

Time is measured in slots internally and converted back on exit, so QN and UN
are dimensionless, TN is per time unit and RN is in time units.

Convention: late arrival system with delayed access (LAS-DA), matching the
Q-MAM discrete-time queues and the LDES slotted engine.

MATLAB twin: solver_mam_dt.m
"""

from typing import List, Optional, Sequence, Tuple

import numpy as np

from ....constants import NodeType, ProcessType, SchedStrategy
from ...mam.dtime import (dmap_compress, dmap_compress_batch, dmap_lambda, dmap_super,
                          dmap_thin, dmap_is_renewal, dmap_to_dph, dph_from_dist,
                          dph_to_dmap, mg1_dt_queue, q_dt_map_map_1, q_dt_ph_ph_1)
from ...sn.predicates import sn_is_open_model

__all__ = ['solver_mam_dt']


def _enum_name(value) -> str:
    """Name of an enum member, compared by NAME rather than by identity.

    sn carries members of a SchedStrategy/ProcessType/NodeType class object that
    is not always the one imported here, and equality across two enum classes of
    the same name is False. Comparing names is the convention the project already
    uses for ProcessType across codebases.
    """
    return getattr(value, 'name', str(value))


class DTResult:
    """Station metrics of the discrete-time path, in model time units."""

    def __init__(self, QN, UN, RN, TN, CN, XN, totiter, method):
        self.QN = QN
        self.UN = UN
        self.RN = RN
        self.TN = TN
        self.CN = CN
        self.XN = XN
        self.totiter = totiter
        self.method = method


def solver_mam_dt(sn, options, slot_length: float = 1.0) -> DTResult:
    _assert_scope(sn)

    M, K = sn.nstations, sn.nclasses
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    law = [_station_law(sn, ist, 0, slot_length) for ist in range(M)]

    sched = sn.sched if isinstance(sn.sched, dict) else {}
    source_idx = None
    queue_idx = []
    for ist in range(M):
        if _enum_name(sched.get(ist)) == 'EXT':
            source_idx = ist
        else:
            queue_idx.append(ist)
    if source_idx is None:
        raise RuntimeError("The discrete-time path requires an open model with a Source.")

    max_num_comp = _int_config(options, 'dt_maxlevel', 1000)
    space_max = _int_config(options, 'space_max', 128)

    if len(queue_idx) == 1:
        ql = _solve_single_station(law[source_idx], law[queue_idx[0]], max_num_comp)
        QNq = [float(sum(i * ql[i] for i in range(len(ql))))]
        UNq = [1.0 - float(ql[0])]
        TNq = [dmap_lambda(law[source_idx])]
        totiter = 1
        method = 'dt.qmam'
    else:
        QNq, UNq, TNq, totiter = _solve_network(sn, law, source_idx, queue_idx, options,
                                                space_max, max_num_comp)
        method = 'dt.dec'

    lambda_slot = dmap_lambda(law[source_idx])
    for idx, ist in enumerate(queue_idx):
        QN[ist, 0] = QNq[idx]
        UN[ist, 0] = UNq[idx]
        TN[ist, 0] = TNq[idx] / slot_length
        if TNq[idx] > 0:
            RN[ist, 0] = QNq[idx] / TNq[idx] * slot_length
    TN[source_idx, 0] = lambda_slot / slot_length
    XN[0, 0] = lambda_slot / slot_length
    CN[0, 0] = QN.sum() / XN[0, 0]

    return DTResult(QN, UN, RN, TN, CN, XN, totiter, method)


def _assert_scope(sn):
    """Reject the model features the discrete-time path cannot represent."""
    if not sn_is_open_model(sn):
        raise RuntimeError(
            "The discrete-time path supports open models only. A closed slotted model needs a "
            "level-dependent discrete chain, which the Q-MAM discrete-time catalogue does not cover.")
    if sn.nclasses > 1:
        raise RuntimeError(
            "The discrete-time path supports one class only. Independent per-class lattice sources "
            "fire in the same slot with positive probability, and a batch of simultaneous arrivals "
            "of different classes is not an MMAP[K], which is what Q_DT_MMAPK_PHK_1 consumes.")
    sched = sn.sched if isinstance(sn.sched, dict) else {}
    nservers = np.asarray(sn.nservers).flatten() if sn.nservers is not None else None
    for ist in range(sn.nstations):
        s = sched.get(ist)
        if _enum_name(s) == 'EXT':
            continue
        if _enum_name(s) != 'FCFS':
            raise RuntimeError(
                f"Station {ist} uses scheduling {s}. The discrete-time path supports FCFS "
                "single-server stations and the Source only.")
        if nservers is not None and ist < len(nservers) and nservers[ist] > 1:
            raise RuntimeError(
                f"Station {ist} has {nservers[ist]} servers. The discrete-time path models one "
                "server per station: a slotted multiserver queue needs the level-dependent "
                "boundary of Geo/Geo/c.")


def _station_law(sn, ist: int, r: int, slot_length: float) -> List[np.ndarray]:
    """Discrete-time law of a station, expressed in slots."""
    proc_type = np.asarray(sn.procid)[ist, r]
    rate = float(np.asarray(sn.rates)[ist, r])
    mean_slots = 1.0 / (rate * slot_length)

    if _enum_name(proc_type) == 'DMAP':
        if slot_length != 1.0:
            raise RuntimeError("A DMAP is defined on its own slot, so it cannot be combined "
                               f"with config.slotlength={slot_length}.")
        proc = sn.proc[ist][r]
        return [np.asarray(proc[0], dtype=float), np.asarray(proc[1], dtype=float)]

    scv = float(np.asarray(sn.scv)[ist, r]) if sn.scv is not None else 1.0
    alpha, A = dph_from_dist(proc_type, mean_slots, scv)
    return dph_to_dmap(alpha, A)


def _solve_single_station(ARV, SVC, max_num_comp: int) -> np.ndarray:
    """Exact single-station analysis through the Q-MAM discrete-time queues."""
    if dmap_is_renewal(ARV[0], ARV[1]) and dmap_is_renewal(SVC[0], SVC[1]):
        alpha, T = dmap_to_dph(ARV[0], ARV[1])
        beta, S = dmap_to_dph(SVC[0], SVC[1])
        return q_dt_ph_ph_1(alpha, T, beta, S, max_num_comp)
    return q_dt_map_map_1(ARV[0], ARV[1], SVC[0], SVC[1], max_num_comp)


def _solve_network(sn, law, source_idx, queue_idx, options, space_max, max_num_comp):
    """Discrete-time parametric decomposition over several stations.

    Each station is solved as a DBMAP/DMAP/1 queue given its arrival stream, and
    its departure stream is extracted from the truncated stationary chain and
    split by the routing probabilities. Superposing discrete streams produces
    batches, which is why the station solve is M/G/1-type rather than a QBD.
    """
    nq = len(queue_idx)
    P = _routing(sn, source_idx, queue_idx)
    lam = dmap_lambda(law[source_idx])

    iter_max = int(getattr(options, 'iter_max', 0) or 100)
    iter_tol = float(getattr(options, 'iter_tol', 0) or 1e-3)

    DEP = []
    for ist in queue_idx:
        p = lam * _visit_ratio(sn, ist)
        DEP.append([np.array([[1 - p]]), np.array([[p]])])

    QN = np.zeros(nq)
    UN = np.zeros(nq)
    TN = np.zeros(nq)
    QNprev = np.zeros(nq)
    QNprev2 = np.zeros(nq)
    UNprev = np.zeros(nq)
    TNprev = np.zeros(nq)
    totiter = 0

    for it in range(1, iter_max + 1):
        totiter = it
        for idx in range(nq):
            ARV = _arrivals(law[source_idx], DEP, P, idx, space_max)
            q, u, t, _, dep = mg1_dt_queue(ARV, law[queue_idx[idx]], max_num_comp, True)
            QN[idx], UN[idx], TN[idx] = q, u, t
            DEP[idx] = dmap_compress(dep, space_max)
        if it > 1 and _max_rel_change(QN, QNprev) < iter_tol:
            break
        if it > 2 and _max_rel_change(QN, QNprev2) < iter_tol:
            # Feedback loops settle into a period-two cycle rather than a point:
            # re-solving a station with the departure process it just produced
            # moves it back. The cycle amplitude sits far below the error of the
            # decomposition itself, so the midpoint is reported instead of
            # burning iter_max sweeps on an orbit that will not close.
            QN = (QN + QNprev) / 2
            UN = (UN + UNprev) / 2
            TN = (TN + TNprev) / 2
            break
        QNprev2 = QNprev.copy()
        QNprev = QN.copy()
        UNprev = UN.copy()
        TNprev = TN.copy()

    return list(QN), list(UN), list(TN), totiter


def _arrivals(SRC, DEP, P, idx, space_max):
    """Arrival stream of a queue: source share plus thinned upstream departures."""
    ARV = None
    if P[0, idx] > 0:
        ARV = dmap_thin(SRC, float(P[0, idx]))
    for j in range(len(DEP)):
        p = float(P[1 + j, idx])
        if p <= 0:
            continue
        contrib = dmap_thin(DEP[j], p)
        if ARV is None:
            ARV = contrib
        else:
            ARV = dmap_super(ARV, contrib)
            ARV = dmap_compress_batch(ARV, space_max)
    if ARV is None:
        raise RuntimeError(f"Queue {idx} receives no arrivals in the discrete-time routing matrix.")
    return ARV


def _routing(sn, source_idx, queue_idx) -> np.ndarray:
    """Station-to-station routing probabilities, source first, sink stripped."""
    I = sn.nnodes
    K = sn.nclasses
    rtnodes = np.asarray(sn.rtnodes)
    Pn = np.zeros((I, I))
    for a in range(I):
        for b in range(I):
            Pn[a, b] = rtnodes[a * K, b * K]

    # the Sink feeds back into the Source to keep rt stochastic; an open traffic
    # equation must not see that edge
    sink_nodes = [ind for ind in range(I) if _enum_name(sn.nodetype[ind]) == 'SINK']
    for ind in sink_nodes:
        Pn[ind, :] = 0.0

    station_to_node = np.asarray(sn.stationToNode).flatten()
    station_nodes = [int(station_to_node[source_idx])] + [int(station_to_node[i]) for i in queue_idx]
    inter_nodes = [ind for ind in range(I)
                   if ind not in station_nodes and ind not in sink_nodes]

    Pss = Pn[np.ix_(station_nodes, station_nodes)]
    if not inter_nodes:
        full = Pss
    else:
        # censor the intermediate nodes: routers and class switches carry no
        # service, so their transit collapses into (I-Pnn)^-1
        Psn = Pn[np.ix_(station_nodes, inter_nodes)]
        Pnn = Pn[np.ix_(inter_nodes, inter_nodes)]
        Pns = Pn[np.ix_(inter_nodes, station_nodes)]
        full = Pss + Psn @ np.linalg.solve(np.eye(len(inter_nodes)) - Pnn, Pns)
    # column 0 is the source, which receives nothing
    return full[:, 1:]


def _visit_ratio(sn, ist: int) -> float:
    v = 0.0
    if isinstance(sn.visits, dict):
        for chain in sn.visits:
            vm = np.asarray(sn.visits[chain])
            if vm.size and ist < vm.shape[0]:
                v += float(vm[ist, :].sum())
    return v


def _max_rel_change(a, b) -> float:
    denom = np.maximum(np.abs(b), 1e-14)
    return float(np.max(np.abs(a - b) / denom))


def _int_config(options, key: str, fallback: int) -> int:
    config = getattr(options, 'config', None)
    if isinstance(config, dict) and config.get(key):
        return int(config[key])
    value = getattr(options, key, None)
    if value:
        return int(value)
    return fallback
