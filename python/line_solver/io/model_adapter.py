"""
Model Adapter for transforming and adapting queueing network models.

Single home for every model-transformation primitive, mirroring MATLAB's
@ModelAdapter class directory:
- Fork-join transformations: fjtag, mmt, Heidelberger-Trivedi ht, sort_forks,
  the path enumerators (paths, paths_cs, find_paths, find_paths_cs)
- Class-level transformations: tag_chain, aggregate_chains, remove_class
- Flow-equivalent server aggregation: aggregate_fes

Port from:
    - matlab/src/io/@ModelAdapter/ModelAdapter.m
    - matlab/src/io/@ModelAdapter/{fjtag,mmt,ht,sortForks}.m
    - matlab/src/io/@ModelAdapter/{paths,pathsCS,findPaths,findPathsCS}.m
    - matlab/src/io/@ModelAdapter/{tagChain,aggregateChains,removeClass}.m
    - matlab/src/io/@ModelAdapter/aggregateFES.m

line_solver.api.io.model_adapter re-exports this module; it must not carry a
second copy of any of these primitives.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import copy
import numpy as np
from typing import Dict, List, Optional, Tuple, Any
from itertools import combinations
from dataclasses import dataclass

from ..api.io.logging import line_warning
from ..api.sn.network_struct import NodeType
from ..api.sn.transforms import get_chain_for_class


@dataclass
class MMTResult:
    """Result of MMT transformation."""
    nonfjmodel: Any  # Transformed model without fork-joins
    fjclassmap: np.ndarray  # fjclassmap[r] = original class index for auxiliary class r
    fjforkmap: np.ndarray  # fjforkmap[r] = fork node index for auxiliary class r
    fanout: np.ndarray  # fanout[r] = number of output jobs for auxiliary class r
    routing_matrix: Any = None  # P matrix used for relink (avoids expensive refresh_struct)
    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
    service_src: Optional[Dict[Tuple[int, int], int]] = None
    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
    immediate_slots: Optional[set] = None
    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
    arrival_src: Optional[Dict[int, int]] = None
    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
    base_model: Any = None
    fork_lambda_init: Optional[np.ndarray] = None


@dataclass
class HTResult:
    """Result of Heidelberger-Trivedi transformation."""
    nonfjmodel: Any  # Transformed model without fork-joins
    fjclassmap: np.ndarray  # fjclassmap[r] = original class index for auxiliary class r
    fjforkmap: np.ndarray  # fjforkmap[r] = fork node index for auxiliary class r
    fj_auxiliary_delays: Dict[int, int]  # fj_auxiliary_delays[j] = auxiliary delay node index for join j


@dataclass
class DeaggregationInfo:
    """Information needed to deaggregate chain-level results back to class-level."""
    alpha: np.ndarray  # Aggregation factors matrix (M x K)
    inchain: List[List[int]]  # Classes in each chain
    original_sn: Any  # Original NetworkStruct
    is_aggregated: bool
    V_chain: Optional[np.ndarray] = None  # Visit ratios per chain
    ST_chain: Optional[np.ndarray] = None  # Service times per chain
    L_chain: Optional[np.ndarray] = None  # Demands per chain
    SCV_chain: Optional[np.ndarray] = None  # SCVs per chain
    N_chain: Optional[np.ndarray] = None  # Population per chain
    lambda_chain: Optional[np.ndarray] = None  # Arrival rates per chain
    is_open_chain: Optional[np.ndarray] = None  # Open chain flags
    refstat: Optional[np.ndarray] = None  # Reference stations
    refstat_chain: Optional[np.ndarray] = None  # Reference stations per chain
    nclasses: int = 0
    nchains: int = 0


@dataclass
class TaggedModelResult:
    """Result of tagging a chain for response time analysis."""
    model: Any  # Tagged model
    tagged_job: Any  # Tagged job class


@dataclass
class FESAggregationInfo:
    """Information for Flow-Equivalent Server aggregation."""
    original_model: Any
    station_subset: List[Any]
    subset_indices: np.ndarray
    complement_indices: np.ndarray
    throughput_table: List[np.ndarray]
    cutoffs: np.ndarray
    stoch_comp_subset: np.ndarray
    stoch_comp_complement: np.ndarray
    isolated_model: Any
    fes_node_idx: int



class ModelAdapter:
    """
    Static class to transform and adapt queueing network models.

    Provides methods for fork-join network transformations using
    MMT (for open networks) and Heidelberger-Trivedi (for closed networks).
    """

    @staticmethod
    def fjtag(model):
        """
        Build a tag-augmented copy of a closed fork-join model for exact
        native analysis by SolverCTMC/SolverSSA.

        Port of matlab/src/io/@ModelAdapter/fjtag.m. For each (fork f,
        class r) pair with matched join j, branch b and tag t, an auxiliary
        transient closed class A(f,r,b,t) with population 0 is created. The
        tag identifies the origin job: the fork firing (after_fj_event) emits
        one sibling per branch in classes A(f,r,1..B,t) using the lowest free
        tag, and the join fires only when all siblings OF THE SAME TAG are
        buffered, releasing one class-r job (after_event_join).

        Returns (fjmodel, fjsn, fjclassmap) where fjclassmap[a] is the 0-based
        original class of auxiliary class a, or -1 for original classes.
        """
        from ..api.fjnative import sn_fj_validate
        from ..lang.classes import ClosedClass
        from ..lang.base import JoinStrategy

        sn = model.get_struct()
        sn_fj_validate(sn)

        K = int(sn.nclasses)
        I = int(sn.nnodes)
        Vnodes = ModelAdapter._compute_node_visits(sn)
        fork_val = int(NodeType.FORK.value) if hasattr(NodeType.FORK, 'value') else int(NodeType.FORK)
        join_val = int(NodeType.JOIN.value) if hasattr(NodeType.JOIN, 'value') else int(NodeType.JOIN)
        nodetype_vals = np.array([int(nt.value) if hasattr(nt, 'value') else int(nt) for nt in sn.nodetype])

        fjmodel = model.copy()
        fjmodel.allow_replace = True
        # routing over the original (pre-augmentation) class set, as list P[r][s]
        P = model.get_linked_routing_matrix()
        if P is None:
            raise RuntimeError('The native CTMC/SSA fork-join implementation requires the routing '
                               'topology to be generated using Network.link.')
        # deep-copy P into a mutable list-of-lists of node matrices
        P = [[np.array(P[r][s], dtype=float, copy=True) if (r < len(P) and s < len(P[r]) and P[r][s] is not None)
              else np.zeros((I, I)) for s in range(K)] for r in range(K)]
        Ip = P[0][0].shape[0]

        fjmodel.reset_network(True)
        fjmodel.reset_struct()
        fjmodel.fork_stateful = True
        fjmodel.is_fj_augmented = True

        fjclassmap = [-1] * K
        fjforkmap = [-1] * K
        fjjoinmap = [-1] * K
        fjbranchmap = [-1] * K
        fjtagmap = [-1] * K

        forkIndexes = np.where(nodetype_vals == fork_val)[0]
        chains = np.asarray(sn.chains) if sn.chains is not None else np.eye(K, dtype=bool)
        njobs = np.asarray(sn.njobs).ravel()
        classprio = np.asarray(sn.classprio).ravel() if getattr(sn, 'classprio', None) is not None else np.zeros(K)
        refstat = np.asarray(sn.refstat).ravel()
        rtnodes = np.asarray(sn.rtnodes)

        # forkinfo rows: (f, j, r, branchheads, branchsets, auxmatrix (BxT), w)
        forkinfo = []

        for f in forkIndexes:
            f = int(f)
            j = int(np.where(sn.fj[f, :])[0][0])
            # tasks emitted per output link at each fork firing
            w = 1
            if sn.nodeparam is not None and f in sn.nodeparam and isinstance(sn.nodeparam[f], dict):
                fo = sn.nodeparam[f].get('fanOut', None)
                if fo is not None:
                    fo_arr = np.atleast_1d(fo)
                    if fo_arr.size > 0:
                        w = int(round(float(fo_arr.flat[0])))
            for r in np.where(Vnodes[f, :] > 0)[0]:
                r = int(r)
                # branch heads: nodes receiving class r directly from the fork
                branchheads = []
                for jnd in range(I):
                    if rtnodes[f * K + r, jnd * K + r] > 0:
                        branchheads.append(jnd)
                B = len(branchheads)
                if B < 2:
                    raise RuntimeError('Degenerate forks with a single output link are not supported '
                                       'by the native CTMC/SSA fork-join implementation.')
                # branch discovery: class-r BFS closure from each head up to the join
                branchsets = [None] * B
                for b in range(B):
                    visitset = [branchheads[b]]
                    frontier = [branchheads[b]]
                    while frontier:
                        cn = frontier.pop(0)
                        if nodetype_vals[cn] == fork_val:
                            raise RuntimeError('Nested fork-join is not supported by the native '
                                               'CTMC/SSA fork-join implementation.')
                        if nodetype_vals[cn] == join_val and cn != j:
                            raise RuntimeError('Overlapping fork-join pairs are not supported by the '
                                               'native CTMC/SSA fork-join implementation.')
                        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                        if cn >= Ip:
                            raise RuntimeError('Class switching between fork and join is not '
                                               'supported by the native CTMC/SSA fork-join implementation.')
                        for jnd in range(I):
                            for s in range(K):
                                if rtnodes[cn * K + r, jnd * K + s] > 0:
                                    if s != r:
                                        raise RuntimeError('Class switching between fork and join is not '
                                                           'supported by the native CTMC/SSA fork-join implementation.')
                                    if jnd != j and jnd not in visitset:
                                        visitset.append(jnd)
                                        frontier.append(jnd)
                    # trap check: every branch node must reach the join within the branch
                    canreach = [j]
                    changed = True
                    while changed:
                        changed = False
                        for cn in visitset:
                            if cn not in canreach:
                                for jnd in list(canreach):
                                    if rtnodes[cn * K + r, jnd * K + r] > 0:
                                        canreach.append(cn)
                                        changed = True
                                        break
                    if not all(cn in canreach for cn in visitset):
                        raise RuntimeError('Fork branches from which the Join is unreachable are not '
                                           'supported by the native CTMC/SSA fork-join implementation.')
                    branchsets[b] = visitset
                # tag pool size = chain population of class r
                c = ModelAdapter._fj_find_chain(chains, r)
                if c is not None:
                    T = int(round(np.sum(njobs[np.where(chains[c, :])[0]])))
                else:
                    T = int(round(njobs[r]))
                if not np.isfinite(T):
                    raise RuntimeError('Chains with infinite population routed through a Fork are not '
                                       'supported by the native CTMC/SSA fork-join implementation.')
                auxmatrix = np.zeros((B, T), dtype=int)
                for t in range(T):
                    for b in range(B):
                        auxname = '%s_f%d_b%d_t%d' % (sn.classnames[r], f, b, t)
                        auxclass = ClosedClass(fjmodel, auxname, 0,
                                               fjmodel._stations[int(refstat[r])], int(classprio[r]))
                        a = int(auxclass._index)
                        auxmatrix[b, t] = a
                        # extend the fj maps to length a+1
                        while len(fjclassmap) <= a:
                            fjclassmap.append(-1); fjforkmap.append(-1); fjjoinmap.append(-1)
                            fjbranchmap.append(-1); fjtagmap.append(-1)
                        fjclassmap[a] = r
                        fjforkmap[a] = f
                        fjjoinmap[a] = j
                        fjbranchmap[a] = b
                        fjtagmap[a] = t
                        # sibling service on the branch copies the original class
                        for cn in branchsets[b]:
                            if sn.isstation[cn] and nodetype_vals[cn] != join_val:
                                svc = model._nodes[cn].get_service(model._classes[r])
                                fjmodel._nodes[cn].set_service(auxclass, copy.deepcopy(svc))
                        # register the auxiliary class at the join input section
                        fjmodel._nodes[j].set_strategy(auxclass, JoinStrategy.STD)
                        fjmodel._nodes[j].set_required(auxclass, -1)
                        # sibling routing: copy the class-r branch routing; the
                        # auxiliary class terminates at the join (no outgoing row)
                        Paa = np.zeros((Ip, Ip))
                        for cn in branchsets[b]:
                            Paa[cn, :] = P[r][r][cn, :]
                        # grow P to include the new aux class index a
                        ModelAdapter._fj_grow_P(P, a, Ip)
                        P[a][a] = Paa
                forkinfo.append((f, j, r, list(branchheads), branchsets, auxmatrix, w))

        # complete the routing list over all class pairs
        Kaug = len(fjmodel._classes)
        for a in range(Kaug):
            ModelAdapter._fj_grow_P(P, a, Ip)
        for a in range(Kaug):
            for b in range(Kaug):
                if P[a][b] is None or np.size(P[a][b]) == 0:
                    P[a][b] = np.zeros((Ip, Ip))

        fjmodel.relink(P)
        # re-initialize the default state over the enlarged class set
        fjmodel.init_default()
        fjsn = fjmodel.get_struct()

        # ---- sn post-edits ----
        fjsn.isfjaugmented = True

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        classcap = np.asarray(fjsn.classcap, dtype=float)
        for (f, j, r, branchheads, branchsets, auxmatrix, w) in forkinfo:
            B, T = auxmatrix.shape
            for b in range(B):
                support = list(branchsets[b]) + [j]
                for t in range(T):
                    a = int(auxmatrix[b, t])
                    classcap[:, a] = 0
                    for cn in support:
                        if sn.isstation[cn]:
                            ist = int(sn.nodeToStation[cn])
                            classcap[ist, a] = w
            # nodeparam fj blocks read by after_event_join / after_fj_event
            if fjsn.nodeparam is None:
                fjsn.nodeparam = {}
            if f not in fjsn.nodeparam or not isinstance(fjsn.nodeparam[f], dict):
                fjsn.nodeparam[f] = {}
            fjf = fjsn.nodeparam[f].get('fj')
            if fjf is None:
                fjf = {'classes': [], 'joins': [], 'auxmatrix': [None] * Kaug, 'branchheads': [None] * Kaug}
                fjsn.nodeparam[f]['fj'] = fjf
            fjf['classes'].append(r)
            fjf['joins'].append(j)
            fjf['auxmatrix'][r] = auxmatrix
            fjf['branchheads'][r] = list(branchheads)
            if j not in fjsn.nodeparam or not isinstance(fjsn.nodeparam[j], dict):
                fjsn.nodeparam[j] = {}
            fjj = fjsn.nodeparam[j].get('fj')
            if fjj is None:
                fjj = {'fork': f, 'origclasses': [], 'auxmatrix': [None] * Kaug, 'required': [None] * Kaug}
                fjsn.nodeparam[j]['fj'] = fjj
            fjj['origclasses'].append(r)
            fjj['auxmatrix'][r] = auxmatrix
            fjj['required'][r] = w * np.ones(B)  # STD: tasksPerLink siblings per branch
        fjsn.classcap = classcap

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        cap = np.asarray(fjsn.cap, dtype=float).ravel().copy()
        for ist in range(fjsn.nstations):
            rowsum = float(np.sum(classcap[ist, :][np.isfinite(classcap[ist, :])]))
            if np.isfinite(cap[ist]) and cap[ist] < rowsum:
                cap[ist] = rowsum
        fjsn.cap = cap

        # fork firing synchronizations: one entry per (fork, class, tag)
        fjsync = []
        for (f, j, r, branchheads, branchsets, auxmatrix, w) in forkinfo:
            T = auxmatrix.shape[1]
            for t in range(T):
                entry = {
                    'fork': f, 'join': j, 'class': r, 'tag': t,
                    'branchheads': list(branchheads),
                    'auxclasses': auxmatrix[:, t].tolist(),
                    'auxall': auxmatrix,        # B x T, for the tag-occupancy scan
                    'weight': w,                # tasksPerLink: siblings per branch
                    'prob': 1.0,
                }
                fjsync.append(entry)
        fjsn.fjsync = fjsync
        fjsn.fjclassmap = np.array(fjclassmap, dtype=int)

        return fjmodel, fjsn, np.array(fjclassmap, dtype=int)

    @staticmethod
    def _fj_find_chain(chains, r):
        idx = np.where(chains[:, r])[0]
        return int(idx[0]) if idx.size > 0 else None

    @staticmethod
    def _fj_grow_P(P, a, Ip):
        """Grow the routing list-of-lists P so index a is valid (Kaug x Kaug)."""
        target = a + 1
        while len(P) < target:
            P.append([])
        for r in range(len(P)):
            while len(P[r]) < target:
                P[r].append(np.zeros((Ip, Ip)))

    @staticmethod
    def mmt(model, fork_lambda: Optional[np.ndarray] = None) -> MMTResult:
        """
        MMT transformation for open fork-join networks.

        Transforms a queueing network containing fork-join subsystems into an
        equivalent network without fork-join nodes by:
        1. Replacing fork nodes with routers
        2. Replacing join nodes with delays (infinite servers)
        3. Creating auxiliary open classes for each forked class
        4. Setting up arrival rates from source to model fork behavior

        Args:
            model: Network model to transform
            fork_lambda: Optional arrival rates for auxiliary classes

        Returns:
            MMTResult with transformed model and mapping information
        """
        from ..lang.network import Network
        from ..lang.nodes import Source, Sink, Delay, Router
        from ..lang.classes import OpenClass
        from ..distributions.continuous import Exp, Immediate, Disabled

        sn = model._sn
        if sn is None:
            model.refresh_struct()
            sn = model._sn

        nclasses = sn.nclasses

        if fork_lambda is None:
            # Use small default arrival rate
            fork_lambda = 1e-6 * np.ones(nclasses)

        fjclassmap = np.zeros(0, dtype=int)
        fjforkmap = np.zeros(0, dtype=int)
        fanout = np.zeros(0, dtype=float)

        # Find fork and join indices
        # Note: Use helper method for proper enum comparison
        fork_indices = ModelAdapter._find_node_indices(sn.nodetype, NodeType.FORK)
        join_indices = ModelAdapter._find_node_indices(sn.nodetype, NodeType.JOIN)

        if len(fork_indices) == 0:
            # No forks, return original model
            return MMTResult(
                nonfjmodel=model,
                fjclassmap=fjclassmap,
                fjforkmap=fjforkmap,
                fanout=fanout
            )

        # Create a copy of the model
        nonfjmodel = model.copy()
        nonfjmodel.allow_replace = True  # Allow replacing Fork/Join with Router/Delay

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        service_src: Dict[Tuple[int, int], int] = {}
        for i in range(sn.nnodes):
            if not sn.isstation[i]:
                continue
            if sn.nodetype[i] in (NodeType.SOURCE, NodeType.FORK, NodeType.JOIN):
                continue
            for r in range(nclasses):
                if model._nodes[i].get_service(model._classes[r]) is not None:
                    service_src[(i, r)] = int(r)
        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        arrival_src: Dict[int, int] = {}
        _src_node = model.get_source()
        if _src_node is not None:
            for r in range(nclasses):
                if _src_node.get_arrival(model._classes[r]) is not None:
                    arrival_src[r] = int(r)
        immediate_slots = set()

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        P_orig = model.get_linked_routing_matrix()
        if P_orig is None:
            raise ValueError(
                "SolverMVA can process fork-join networks only if their routing "
                "topology has been generated using Network.link()."
            )
        # Deep copy P so modifications don't affect the original model
        P = [[P_orig[r][s].copy() for s in range(len(P_orig[r]))] for r in range(len(P_orig))]

        # Reset network structure
        nonfjmodel.reset_network(True)

        # Compute node visits
        Vnodes = ModelAdapter._compute_node_visits(sn)

        # Store forked classes for each fork
        forked_classes = {}
        orig_fanout = np.zeros((sn.nnodes, nclasses), dtype=int)

        # Replace forks with routers and adjust routing probabilities
        for f in fork_indices:
            fork_node = model._nodes[f]
            for r in range(nclasses):
                # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                destinations = set()
                for s in range(nclasses):
                    nonzero = np.where(P[r][s][f, :] > 1e-12)[0]
                    destinations.update(nonzero.tolist())
                fanout_count = len(destinations)
                if fanout_count > 0:
                    orig_fanout[f, r] = int(fanout_count)
                    # Divide routing probabilities by fanout across all class pairs
                    for s in range(nclasses):
                        P[r][s][f, :] = P[r][s][f, :] / orig_fanout[f, r]
                else:
                    orig_fanout[f, r] = 0

            # Replace fork with router
            router = Router(nonfjmodel, fork_node.name)
            nonfjmodel._nodes[f] = router
            forked_classes[f] = np.where(Vnodes[f, :] > 0)[0]

        # Replace joins with delays
        for j in join_indices:
            join_node = model._nodes[j]
            delay = Delay(nonfjmodel, join_node.name)
            nonfjmodel._nodes[j] = delay
            # Update stations list
            station_idx = model._nodes[j].get_station_index0()
            if station_idx is not None:
                nonfjmodel._stations[station_idx] = delay
            # Set immediate service for all classes
            for c in range(len(nonfjmodel._classes)):
                delay.set_service(nonfjmodel._classes[c], Immediate())
                immediate_slots.add((j, c))

        # Ensure source and sink exist
        source = nonfjmodel.get_source()
        sink = nonfjmodel.get_sink()
        if source is None:
            source = Source(nonfjmodel, 'Source')
        if sink is None:
            sink = Sink(nonfjmodel, 'Sink')

        # Expand P to include new nodes
        nnodes = len(nonfjmodel._nodes)
        for r in range(nclasses):
            for s in range(nclasses):
                old_P = P[r][s]
                new_P = np.zeros((nnodes, nnodes))
                new_P[:old_P.shape[0], :old_P.shape[1]] = old_P
                P[r][s] = new_P

        nonfjmodel.connections = np.zeros((nnodes, nnodes))

        oclass_list = []

        # Process each fork
        for f in fork_indices:
            # Find join associated to fork f
            join_idx = np.where(sn.fj[f, :] > 0)[0]
            if len(join_idx) > 1:
                raise ValueError(
                    "SolverMVA supports at present only a single join station per fork node."
                )
            # MATLAB: joinIdx = find(sn.fj(f,:)); if empty, join-related ops are no-ops
            # We must NOT skip — auxiliary classes still need to be created
            has_join = len(join_idx) > 0
            join_idx = join_idx[0] if has_join else None

            # Find chains associated to forked classes
            # Handles both 1D (chains[r] = chain_id) and 2D (chains[c,r] > 0) formats
            forked_chain_indices = []
            forked = forked_classes.get(f, [])
            for class_idx in forked:
                chain_id = get_chain_for_class(sn.chains, class_idx)
                if chain_id >= 0 and chain_id not in forked_chain_indices:
                    forked_chain_indices.append(chain_id)

            for fc in forked_chain_indices:
                oclass_map = {}
                # Get classes in this chain
                if hasattr(sn, 'inchain') and sn.inchain is not None and fc in sn.inchain:
                    inchain = sn.inchain[fc]
                else:
                    # Fallback: find classes from chains matrix
                    chains_arr = np.asarray(sn.chains)
                    if chains_arr.ndim == 2:
                        inchain = np.where(chains_arr[fc, :] > 0)[0]
                    else:
                        inchain = np.where(chains_arr == fc)[0]

                for r in inchain:
                    # Create auxiliary open class
                    class_name = f"{nonfjmodel._classes[r].name}.{nonfjmodel._nodes[f].name}"
                    new_class = OpenClass(nonfjmodel, class_name)
                    oclass_map[r] = new_class

                    # Update mappings
                    new_idx = new_class.get_index0()
                    fjclassmap = np.append(fjclassmap, r)
                    fjforkmap = np.append(fjforkmap, f)

                    # Compute fanout
                    fork_node = model._nodes[f]
                    tasks_per_link = getattr(fork_node, '_tasks_per_link', None)
                    if tasks_per_link is None:
                        tasks_per_link = 1
                    new_fanout = orig_fanout[f, r] * tasks_per_link
                    fanout = np.append(fanout, new_fanout)

                    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                    if orig_fanout[f, r] == 0:
                        source.set_arrival(new_class, Disabled())
                    else:
                        source.set_arrival(new_class, Exp(fork_lambda[r]))

                    # Set service times at stations
                    for i in range(sn.nnodes):
                        if sn.isstation[i]:
                            node_type = sn.nodetype[i]
                            if node_type == NodeType.JOIN:
                                # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                                nonfjmodel._nodes[i].set_service(new_class, Immediate())
                                immediate_slots.add((i, new_idx))
                            elif node_type not in (NodeType.SOURCE, NodeType.FORK):
                                # Copy service from original class
                                orig_service = model._nodes[i].get_service(model._classes[r])
                                if orig_service is not None:
                                    nonfjmodel._nodes[i].set_service(new_class, copy.deepcopy(orig_service))
                                    service_src[(i, new_idx)] = int(r)

                    oclass_list.append(new_class)

                # Set routing for auxiliary classes
                for r in inchain:
                    for s in inchain:
                        oclass_r = oclass_map.get(r)
                        oclass_s = oclass_map.get(s)
                        if oclass_r is not None and oclass_s is not None:
                            # Extend P to new classes
                            P = ModelAdapter._ensure_p_size(P, max(oclass_r.get_index0(), oclass_s.get_index0()) + 1, nnodes)
                            # Copy routing from original class
                            P[oclass_r.get_index0()][oclass_s.get_index0()] = P[r][s].copy()

                # Set source and sink routing for auxiliary classes
                # MATLAB mmt.m lines 122-128: first zero out source/join rows for ALL pairs
                for r in inchain:
                    for s in inchain:
                        oclass_r = oclass_map.get(r)
                        oclass_s = oclass_map.get(s)
                        if oclass_r is not None and oclass_s is not None:
                            P = ModelAdapter._ensure_p_size(P, max(oclass_r.get_index0(), oclass_s.get_index0()) + 1, nnodes)
                            # Zero out source row for all class pairs
                            P[oclass_r.get_index0()][oclass_s.get_index0()][source.get_index0(), :] = 0.0
                            # Zero out join row (no-op if no join, matching MATLAB empty indexing)
                            if has_join:
                                P[oclass_r.get_index0()][oclass_s.get_index0()][join_idx, :] = 0.0
                # Then set source->fork and join->sink for diagonal only
                for r in inchain:
                    oclass_r = oclass_map.get(r)
                    if oclass_r is not None:
                        P[oclass_r.get_index0()][oclass_r.get_index0()][source.get_index0(), f] = 1.0
                        # Set join->sink routing (no-op if no join)
                        if has_join:
                            P[oclass_r.get_index0()][oclass_r.get_index0()][join_idx, sink.get_index0()] = 1.0

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        from ..lang.nodes import ClassSwitch
        new_nclasses = len(nonfjmodel._classes)
        for node in nonfjmodel._nodes:
            if isinstance(node, ClassSwitch):
                old_matrix = node._switch_matrix
                if old_matrix is not None:
                    old_size = old_matrix.shape[0]
                    if old_size < new_nclasses:
                        new_matrix = np.zeros((new_nclasses, new_nclasses))
                        new_matrix[:old_size, :old_size] = old_matrix
                        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                        for aux_local in range(len(fjclassmap)):
                            orig_r = int(fjclassmap[aux_local])
                            if orig_r < 0:
                                continue
                            abs_aux_r = nclasses + aux_local
                            if abs_aux_r >= new_nclasses:
                                continue
                            # Find what class orig_r switches to
                            for orig_s in range(old_size):
                                if old_matrix[orig_r, orig_s] > 0:
                                    # Find the auxiliary class corresponding
                                    # to orig_s from the same fork as aux_local.
                                    fork_of_aux = int(fjforkmap[aux_local])
                                    target_aux = -1
                                    for other_local in range(len(fjclassmap)):
                                        if (int(fjclassmap[other_local]) == orig_s
                                                and int(fjforkmap[other_local]) == fork_of_aux):
                                            target_aux = nclasses + other_local
                                            break
                                    if target_aux >= 0 and target_aux < new_nclasses:
                                        new_matrix[abs_aux_r, target_aux] += (
                                            old_matrix[orig_r, orig_s]
                                        )
                                    else:
                                        # No matching auxiliary; keep identity
                                        new_matrix[abs_aux_r, abs_aux_r] += (
                                            old_matrix[orig_r, orig_s]
                                        )
                        node._switch_matrix = new_matrix

        # Relink the model with updated routing
        nonfjmodel.relink(P)

        # Mark nonfjmodel to skip recursive _refresh_fork_join_nodevisits
        nonfjmodel._skip_fj_nodevisits = True

        # Note: Routing for auxiliary classes is already set up correctly via P matrix
        # The routing probabilities have been divided by fanout when replacing forks

        return MMTResult(
            nonfjmodel=nonfjmodel,
            fjclassmap=fjclassmap,
            fjforkmap=fjforkmap,
            fanout=fanout,
            routing_matrix=P,
            service_src=service_src,
            immediate_slots=immediate_slots,
            arrival_src=arrival_src,
            base_model=model,
            fork_lambda_init=(None if fork_lambda is None else np.array(fork_lambda, dtype=float)),
        )

    @staticmethod
    def refresh_services_from_base(mmt_result) -> bool:
        """Re-feed a previously built ``mmt_result.nonfjmodel`` from the current
        service and arrival parameters of the base model it was derived from.

        This is what makes an MMT transformation reusable across the iterations
        of an outer fixed point (SolverLN re-solves each layer many times, and
        only the rates change between iterations -- the fork topology, on which
        the transformation depends, does not). Rebuilding instead costs a full
        deep copy of the layer per iteration.

        Slots recorded in ``service_src`` / ``arrival_src`` are rewritten from
        the base model. Slots the transformation owns are never read from the
        base -- doing so is precisely the corruption this indirection exists to
        prevent -- but they are not simply left alone either, because the fork
        loop mutates some of them as it iterates: it overwrites every join with
        its current synchronisation delay. Those (``immediate_slots``, and the
        auxiliary arrivals driven by ``fork_lambda``) are reset to the initial
        value a cold transformation would have set, so that reuse enters the
        fork loop from the same state rather than from the previous outer
        iteration's converged one.

        Returns False when the recorded provenance cannot be applied, in which
        case the caller must fall back to a cold transformation.
        """
        if mmt_result is None or mmt_result.service_src is None:
            return False
        base = mmt_result.base_model
        nonfj = mmt_result.nonfjmodel
        if base is None or nonfj is None:
            return False

        nclasses = len(base._classes)
        for (i, c), r in mmt_result.service_src.items():
            if i >= len(nonfj._nodes) or c >= len(nonfj._classes) or r >= nclasses:
                return False
            svc = base._nodes[i].get_service(base._classes[r])
            if svc is None:
                return False
            # deepcopy mirrors the cold path, which hands every slot its own
            # distribution object rather than aliasing one across classes.
            nonfj._nodes[i].set_service(nonfj._classes[c], copy.deepcopy(svc))

        from ..distributions.continuous import Immediate as _Immediate
        for (i, c) in (mmt_result.immediate_slots or ()):
            if i >= len(nonfj._nodes) or c >= len(nonfj._classes):
                return False
            nonfj._nodes[i].set_service(nonfj._classes[c], _Immediate())

        if mmt_result.arrival_src:
            src_base = base.get_source()
            src_nonfj = nonfj.get_source()
            if src_base is None or src_nonfj is None:
                return False
            for c, r in mmt_result.arrival_src.items():
                arv = src_base.get_arrival(base._classes[r])
                if arv is None:
                    return False
                src_nonfj.set_arrival(nonfj._classes[c], copy.deepcopy(arv))

        # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
        fl = mmt_result.fork_lambda_init
        if fl is not None:
            from ..distributions.continuous import Exp, Disabled
            src_nonfj = nonfj.get_source()
            if src_nonfj is None:
                return False
            for r_aux in range(len(mmt_result.fjclassmap)):
                aux_idx = nclasses + r_aux
                if aux_idx >= len(nonfj._classes):
                    return False
                r = int(mmt_result.fjclassmap[r_aux])
                if mmt_result.fanout[r_aux] == 0:
                    src_nonfj.set_arrival(nonfj._classes[aux_idx], Disabled())
                elif r < len(fl):
                    src_nonfj.set_arrival(nonfj._classes[aux_idx], Exp(fl[r]))
        return True

    @staticmethod
    def ht(model) -> HTResult:
        """
        Heidelberger-Trivedi transformation for closed fork-join networks.

        Transforms a closed queueing network containing fork-join subsystems
        into an equivalent network without fork-join nodes by:
        1. Replacing fork nodes with routers
        2. Replacing join nodes with delays
        3. Creating auxiliary closed classes for each parallel branch
        4. Adding auxiliary delay nodes to model sojourn times

        Based on: PHILIP HEIDELBERGER and KISHOR S. TRIVEDI,
        "Analytic Queueing Models for Programs with Internal Concurrency"

        Args:
            model: Network model to transform

        Returns:
            HTResult with transformed model and mapping information
        """
        from ..lang.network import Network
        from ..lang.nodes import Delay, Router
        from ..lang.classes import ClosedClass
        from ..distributions.continuous import Immediate

        sn = model._sn
        if sn is None:
            model.refresh_struct()
            sn = model._sn

        # fjclassmap[i] = -1 for original classes, = original_class_index for auxiliary classes
        # fjforkmap[i] = fork_index for auxiliary classes
        nclasses_orig = len(model._classes)
        fjclassmap = np.full(nclasses_orig, -1, dtype=int)
        fjforkmap = np.full(nclasses_orig, -1, dtype=int)
        fj_auxiliary_delays = {}

        # Find fork and join indices
        # Note: Use list iteration instead of np.where for enum comparison
        fork_indices = ModelAdapter._find_node_indices(sn.nodetype, NodeType.FORK)

        if len(fork_indices) == 0:
            return HTResult(
                nonfjmodel=model,
                fjclassmap=fjclassmap,
                fjforkmap=fjforkmap,
                fj_auxiliary_delays=fj_auxiliary_delays
            )

        # Create copy and allow replacements
        nonfjmodel = model.copy()
        nonfjmodel.allow_replace = True

        P = nonfjmodel.get_linked_routing_matrix()
        nonfjmodel.reset_network(True)

        # Compute node visits
        Vnodes = ModelAdapter._compute_node_visits(sn)

        forked_classes = {}

        # Replace forks with routers (add_node handles replacement via allow_replace)
        for f in fork_indices:
            Router(nonfjmodel, nonfjmodel._nodes[f].name)
            forked_classes[f] = np.where(Vnodes[f, :] > 0)[0]

        # Replace joins with delays and add auxiliary delays (add_node handles replacement)
        for j in ModelAdapter._find_node_indices(sn.nodetype, NodeType.JOIN):
            delay = Delay(nonfjmodel, nonfjmodel._nodes[j].name)

            # Set immediate service for existing classes
            for c in range(len(nonfjmodel._classes)):
                delay.set_service(nonfjmodel._classes[c], Immediate())

            # Add auxiliary delay node
            aux_delay = Delay(nonfjmodel, f'Auxiliary Delay - {nonfjmodel._nodes[j].name}')
            fj_auxiliary_delays[j] = aux_delay.get_index0()

            # Set immediate service for existing classes at auxiliary delay
            for r in range(len(nonfjmodel._classes)):
                aux_delay.set_service(nonfjmodel._classes[r], Immediate())

            # Update routing: join -> aux_delay, aux_delay -> original successors
            nclasses = len(P)
            for r in range(nclasses):
                for s in range(nclasses):
                    # aux_delay inherits successors of join
                    P = ModelAdapter._ensure_p_size(P, nclasses, aux_delay.get_index0() + 1)
                    P[r][s] = ModelAdapter._ensure_matrix_size(P[r][s], aux_delay.get_index0() + 1)
                    P[r][s][aux_delay.get_index0(), :] = P[r][s][j, :]
                    P[r][s][:, aux_delay.get_index0()] = 0.0
                # join -> aux_delay
                P[r][r][j, :] = 0.0
                P[r][r][j, aux_delay.get_index0()] = 1.0

        nonfjmodel.connections = np.zeros((len(nonfjmodel._nodes), len(nonfjmodel._nodes)))

        # Create auxiliary classes for each fork
        for f in fork_indices:
            join_idx = np.where(sn.fj[f, :] > 0)[0]
            if len(join_idx) > 1:
                raise ValueError(
                    "SolverMVA supports at present only a single join station per fork node."
                )
            if len(join_idx) == 0:
                continue
            join_idx = join_idx[0]

            # Find forked chains - use sn.inchain for Python
            forked_chain_indices = []
            for chain_idx in range(sn.nchains):
                # Get classes in this chain from sn.inchain
                chain_classes = sn.inchain.get(chain_idx, np.array([]))
                # Check if any forked classes are in this chain
                forked = forked_classes.get(f, [])
                if len(chain_classes) > 0 and len(forked) > 0:
                    if np.any(np.isin(forked, chain_classes)):
                        forked_chain_indices.append(chain_idx)

            for fc in forked_chain_indices:
                aux_classes = {}
                # Use sn.inchain to get classes in chain
                inchain = sn.inchain.get(fc, np.array([]))

                for r in inchain:
                    if sn.nodevisits is not None and fc < len(sn.nodevisits):
                        if sn.nodevisits[fc][f, r] == 0:
                            continue

                    # Get parallel branches from routing matrix
                    branches = ModelAdapter._get_fork_branches(model, f, r)
                    if len(branches) == 0:
                        continue

                    parallel_branches = len(branches)

                    for par in range(parallel_branches):
                        fork_node = model._nodes[f]
                        tasks_per_link = fork_node.get_tasks_per_link()
                        if tasks_per_link is None:
                            tasks_per_link = 1
                        elif np.any(tasks_per_link > 1):
                            raise ValueError(
                                "Multiple tasks per link are not supported in H-T."
                            )
                        else:
                            tasks_per_link = 1

                        # Get population for original class
                        orig_class = model._classes[r]
                        if not hasattr(orig_class, 'population'):
                            raise ValueError(
                                "H-T method can be used only on closed models."
                            )

                        aux_population = tasks_per_link * orig_class.population

                        # Create auxiliary closed class
                        class_name = f"{nonfjmodel._classes[r].name}.{nonfjmodel._nodes[f].name}.B{par + 1}"
                        ref_station = nonfjmodel._nodes[fj_auxiliary_delays[join_idx]]
                        aux_class = ClosedClass(nonfjmodel, class_name, aux_population, ref_station, 0)
                        aux_classes[(r, par)] = aux_class

                        # Extend fjclassmap/fjforkmap to accommodate new class index
                        aux_idx = aux_class.get_index0()
                        if aux_idx >= len(fjclassmap):
                            new_size = aux_idx + 1
                            fjclassmap = np.concatenate([fjclassmap, np.full(new_size - len(fjclassmap), -1, dtype=int)])
                            fjforkmap = np.concatenate([fjforkmap, np.full(new_size - len(fjforkmap), -1, dtype=int)])
                        fjclassmap[aux_idx] = r
                        fjforkmap[aux_idx] = f

                        # Set service rates
                        for i in range(sn.nnodes):
                            if sn.isstation[i]:
                                node_type = sn.nodetype[i]
                                if node_type == NodeType.JOIN:
                                    nonfjmodel._nodes[i].set_service(aux_class, Immediate())
                                elif node_type not in (NodeType.SOURCE, NodeType.FORK):
                                    orig_service = model._nodes[i].get_service(model._classes[r])
                                    if orig_service is not None:
                                        nonfjmodel._nodes[i].set_service(aux_class, copy.deepcopy(orig_service))

                        # Set immediate service at auxiliary delay
                        nonfjmodel._nodes[fj_auxiliary_delays[join_idx]].set_service(aux_class, Immediate())

                # Set routing for auxiliary classes
                for r in inchain:
                    if sn.nodevisits is not None and fc < len(sn.nodevisits):
                        if sn.nodevisits[fc][f, r] == 0:
                            continue

                    # Get parallel branches from routing matrix
                    branches = ModelAdapter._get_fork_branches(model, f, r)
                    if len(branches) == 0:
                        continue

                    parallel_branches = len(branches)

                    for s in inchain:
                        if sn.nodevisits is not None and fc < len(sn.nodevisits):
                            if sn.nodevisits[fc][f, s] == 0:
                                continue

                        for par in range(parallel_branches):
                            aux_r = aux_classes.get((r, par))
                            aux_s = aux_classes.get((s, par))
                            if aux_r is None or aux_s is None:
                                continue

                            # Extend P for new classes
                            P = ModelAdapter._ensure_p_size(
                                P,
                                max(aux_r.get_index0(), aux_s.get_index0()) + 1,
                                len(nonfjmodel._nodes)
                            )

                            # Copy routing
                            P[aux_r.get_index0()][aux_s.get_index0()] = P[r][s].copy()

                            # Get the specific branch target for this parallel branch
                            branch_target = branches[par]
                            branch_target_idx = branch_target.get_index0() if hasattr(branch_target, 'get_index0') else branch_target._node_index

                            # Fork -> specific branch target
                            P[aux_r.get_index0()][aux_s.get_index0()][f, :] = 0.0
                            P[aux_r.get_index0()][aux_s.get_index0()][f, branch_target_idx] = 1.0

                            # Join -> aux_delay
                            P[aux_r.get_index0()][aux_s.get_index0()][join_idx, fj_auxiliary_delays[join_idx]] = 1.0

                            # aux_delay -> fork
                            P[aux_r.get_index0()][aux_s.get_index0()][fj_auxiliary_delays[join_idx], :] = 0.0
                            P[aux_r.get_index0()][aux_s.get_index0()][fj_auxiliary_delays[join_idx], f] = 1.0

                        # Route original classes straight to join (MATLAB lines 116-118)
                        P[r][s][f, :] = 0.0
                        P[r][s][f, join_idx] = 1.0

        nonfjmodel.relink(P)

        return HTResult(
            nonfjmodel=nonfjmodel,
            fjclassmap=fjclassmap,
            fjforkmap=fjforkmap,
            fj_auxiliary_delays=fj_auxiliary_delays
        )

    @staticmethod
    def sort_forks(sn, fjforkmap: np.ndarray, fjclassmap: np.ndarray,
                   nonfjmodel, routing_matrix=None) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sort fork nodes and identify outer/nested forks.

        Args:
            sn: Network structure
            fjforkmap: Fork node mapping for auxiliary classes
            fjclassmap: Class mapping for auxiliary classes
            nonfjmodel: Transformed model
            routing_matrix: Optional pre-computed P matrix (avoids expensive refresh_struct)

        Returns:
            Tuple of (forks, parents) arrays where:
            - forks[f, r] = 1 if fork f is an outer fork for class r
            - parents[f] = parent fork index for fork f
        """
        max_class = int(np.max(fjclassmap)) + 1 if len(fjclassmap) > 0 else 1
        forks = np.zeros((len(sn.nodetype), max_class), dtype=int)

        fork_indices = ModelAdapter._find_node_indices(sn.nodetype, NodeType.FORK)
        forks[fork_indices, :] = 1

        parents = np.zeros(len(sn.nodetype), dtype=int)
        parents[fork_indices] = fork_indices

        P = routing_matrix if routing_matrix is not None else nonfjmodel.get_linked_routing_matrix()

        for f in fork_indices:
            join_idx_arr = np.where(sn.fj[f, :] > 0)[0]
            if len(join_idx_arr) == 0:
                continue
            join_idx = join_idx_arr[0]

            fork_aux_classes = np.where(fjforkmap == f)[0]

            for s in fork_aux_classes:
                r = fjclassmap[s]
                discovered_forks = np.zeros(len(sn.nodetype), dtype=int)
                discovered_forks[fork_indices] = 1

                # Get routing matrix for class r
                if r < len(P) and r < len(P[r]):
                    conn = P[r][r]
                else:
                    continue

                nested = ModelAdapter._nested_forks(f, join_idx, conn, discovered_forks.copy(), sn)
                forks[:, r] = forks[:, r] & nested
                parents[nested == 0] = parents[f]

        return forks, parents

    @staticmethod
    def _nested_forks(start_node: int, end_node: int, conn: np.ndarray,
                      forks: np.ndarray, sn) -> np.ndarray:
        """Find nested forks between start and end nodes."""
        if start_node == end_node:
            return forks

        for i in np.where(conn[start_node, :] > 0)[0]:
            if sn.nodetype[i] == NodeType.FORK:
                forks[i] = 0
            forks = forks & ModelAdapter._nested_forks(i, end_node, conn, forks, sn)

        return forks

    @staticmethod
    def find_paths(sn, P: np.ndarray, start: int, end_node: int, r: int,
                   to_merge: List[int], QN: np.ndarray, TN: np.ndarray,
                   current_time: float, fjclassmap: np.ndarray,
                   fjforkmap: np.ndarray, nonfjmodel) -> np.ndarray:
        """
        Find response times along each path from start to end_node.

        Args:
            sn: Network structure
            P: Routing matrix for class r
            start: Starting node index
            end_node: Ending node index
            r: Class index
            to_merge: Classes to merge for queueing/throughput calculation
            QN: Queue length matrix
            TN: Throughput matrix
            current_time: Accumulated response time
            fjclassmap: Class mapping for auxiliary classes
            fjforkmap: Fork mapping for auxiliary classes
            nonfjmodel: Transformed model

        Returns:
            Array of response times for each parallel path
        """
        if start == end_node:
            # Subtract synchronization time at join
            station = sn.nodeToStation[start]
            if not np.isnan(station):
                q_len = np.sum(QN[int(station), to_merge])
                tput = np.sum(TN[int(station), to_merge])
                if tput > 0:
                    return np.array([current_time - q_len / tput])
            return np.array([current_time])

        ri = []
        for i in np.where(P[start, :] > 0)[0]:
            q_len = 0
            tput = 1

            station = sn.nodeToStation[i] if i < len(sn.nodeToStation) else np.nan
            if not np.isnan(station):
                q_len = np.sum(QN[int(station), to_merge])
                tput = np.sum(TN[int(station), to_merge])
                if tput <= 0:
                    tput = 1

            if sn.nodetype[i] == NodeType.FORK:
                # Handle nested fork
                join_idx_arr = np.where(sn.fj[i, :] > 0)[0]
                if len(join_idx_arr) > 0:
                    join_idx = join_idx_arr[0]
                    s_arr = np.where((fjforkmap == i) & (fjclassmap == r))[0]
                    if len(s_arr) > 0:
                        orig_nclasses = len(nonfjmodel._classes) - len(fjclassmap)
                        s = orig_nclasses + s_arr[0]
                    else:
                        s = r

                    # Recursively find paths within nested fork
                    paths = ModelAdapter.find_paths(
                        sn, P, i, join_idx, r,
                        to_merge + [s], QN, TN, 0,
                        fjclassmap, fjforkmap, nonfjmodel
                    )

                    # The inner join fires on the k-th branch completion, k = the
                    # branch count on a standard join and the declared quorum on a
                    # PARTIAL one.
                    from ..api.fjnative import sn_join_quorum
                    d0 = ModelAdapter._compute_sync_delay(
                        paths, sn_join_quorum(sn, int(join_idx), r, len(paths)))

                    # Set service at join for merged classes
                    from ..distributions.continuous import Exp
                    sync_time = max(0, d0 - np.mean(paths))
                    for cls in to_merge + [s]:
                        if cls < len(nonfjmodel._classes):
                            nonfjmodel._nodes[join_idx].set_service(
                                nonfjmodel._classes[cls],
                                Exp.fit_mean(sync_time) if sync_time > 0 else Exp(1e10)
                            )

                    # Continue from join
                    nested_ri = ModelAdapter.find_paths(
                        sn, P, join_idx, end_node, r,
                        to_merge, QN, TN, current_time + d0,
                        fjclassmap, fjforkmap, nonfjmodel
                    )
                    ri.extend(nested_ri)
            else:
                # Regular node, add response time and continue
                nested_ri = ModelAdapter.find_paths(
                    sn, P, i, end_node, r,
                    to_merge, QN, TN, current_time + q_len / tput,
                    fjclassmap, fjforkmap, nonfjmodel
                )
                ri.extend(nested_ri)

        return np.array(ri) if ri else np.array([current_time])

    @staticmethod
    def find_paths_cs(sn, P: np.ndarray, cur_node: int, end_node: int,
                      cur_class: int, to_merge: List[int], QN: np.ndarray,
                      TN: np.ndarray, current_time: float,
                      fjclassmap: np.ndarray, fjforkmap: np.ndarray,
                      nonfjmodel, on_path=None) -> np.ndarray:
        """
        Find response times along paths with class switches.

        Traverses the routing matrix from cur_node to end_node, computing
        response times along each parallel path. This is the core algorithm
        for MMT fork-join analysis.

        Args:
            sn: Network structure (original model)
            P: Combined routing matrix (nclasses*nnodes x nclasses*nnodes)
            cur_node: Current node index (0-based)
            end_node: Ending node index (0-based)
            cur_class: Current class index (0-based)
            to_merge: List of classes to merge for Q/T computation
            QN: Queue length matrix [stations x classes]
            TN: Throughput matrix [stations x classes]
            current_time: Accumulated response time
            fjclassmap: fjclassmap[aux_idx] = original class index
            fjforkmap: fjforkmap[aux_idx] = fork node index
            nonfjmodel: Transformed model
            on_path: (class, node) states already held by the current path

        Returns:
            Array of response times, one per parallel path
        """
        if cur_node == end_node:
            # At destination: subtract synchronization time at join
            station = sn.nodeToStation[cur_node]
            if not np.isnan(station):
                q_len = np.sum(QN[int(station), to_merge])
                tput = np.sum(TN[int(station), to_merge])
                if tput > 0:
                    return np.array([current_time - q_len / tput])
            return np.array([current_time])

        ri = []
        nfjsn = nonfjmodel._sn
        if nfjsn is None:
            nonfjmodel.refresh_struct()
            nfjsn = nonfjmodel._sn

        # Get original number of nodes for indexing into combined routing matrix
        # rtorig contains the original routing matrix dimensions
        if hasattr(nfjsn, 'rtorig') and nfjsn.rtorig is not None and len(nfjsn.rtorig) > 0:
            if len(nfjsn.rtorig[0]) > 0 and nfjsn.rtorig[0][0] is not None:
                orig_nodes = nfjsn.rtorig[0][0].shape[0]
            else:
                orig_nodes = P.shape[0] // max(len(nfjsn.classnames), 1)
        else:
            orig_nodes = P.shape[0] // max(len(nfjsn.classnames), 1)

        # Compute row index in combined routing matrix
        # For 0-based indexing: row = cur_class * orig_nodes + cur_node
        row_idx = cur_class * orig_nodes + cur_node
        if row_idx < 0 or row_idx >= P.shape[0]:
            return np.array([current_time])

        # enumerate SIMPLE paths only: a call mean above one adds a geometric
        # loop through the .Aux class, so the fork-to-join graph is cyclic
        if on_path is None:
            on_path = frozenset()
        if row_idx in on_path:
            return np.array([])
        on_path = on_path | {row_idx}

        # Find all transitions from current (class, node)
        for transition in np.where(P[row_idx, :] > 0)[0]:
            cur_merge = list(to_merge)

            # Decode transition into (next_class, next_node) using 0-based indexing
            next_class = transition // orig_nodes
            next_node = transition % orig_nodes

            # Update merged class list
            cur_merge[0] = next_class

            q_len = 0
            tput = 1

            if next_node < len(sn.nodeToStation):
                station = sn.nodeToStation[next_node]
                if not np.isnan(station):
                    q_len = np.sum(QN[int(station), cur_merge])
                    tput = np.sum(TN[int(station), cur_merge])
                    if tput <= 0:
                        tput = 1

            # Check if next node is a (nested) Fork
            if sn.nodetype[next_node] == NodeType.FORK:
                join_idx_arr = np.where(sn.fj[next_node, :] > 0)[0]
                if len(join_idx_arr) > 0:
                    join_idx = join_idx_arr[0]
                    # see _kb/12-interfaces-and-docs.md (ModelAdapter mmt/fjtag port notes) for rationale
                    s_arr = np.where((fjforkmap == next_node) & (fjclassmap == next_class))[0]
                    if len(s_arr) > 0:
                        orig_nclasses = len(nonfjmodel._classes) - len(fjclassmap)
                        s = orig_nclasses + s_arr[0]
                    else:
                        s = next_class

                    # Recursively find paths within nested fork-join
                    paths = ModelAdapter.find_paths_cs(
                        sn, P, next_node, join_idx, next_class,
                        cur_merge + [s], QN, TN, 0,
                        fjclassmap, fjforkmap, nonfjmodel, on_path
                    )

                    # Compute E[max] for synchronization
                    d0 = ModelAdapter._compute_sync_delay(paths)

                    # Set sync delay at join for nested fork
                    from ..distributions.continuous import Exp
                    sync_time = max(0, d0 - np.mean(paths))
                    for cls in cur_merge + [s]:
                        if cls < len(nonfjmodel._classes):
                            nonfjmodel._nodes[join_idx].set_service(
                                nonfjmodel._classes[cls],
                                Exp.fit_mean(sync_time) if sync_time > 0 else Exp(1e10)
                            )

                    # Continue from join to end_node
                    nested_ri = ModelAdapter.find_paths_cs(
                        sn, P, join_idx, end_node, next_class,
                        cur_merge, QN, TN, current_time + d0,
                        fjclassmap, fjforkmap, nonfjmodel, on_path
                    )
                    ri.extend(nested_ri)
            else:
                # Regular node: add response time and continue
                nested_ri = ModelAdapter.find_paths_cs(
                    sn, P, next_node, end_node, next_class,
                    cur_merge, QN, TN, current_time + q_len / tput,
                    fjclassmap, fjforkmap, nonfjmodel, on_path
                )
                ri.extend(nested_ri)

        return np.array(ri) if ri else np.array([current_time])

    @staticmethod
    def paths(sn, P: np.ndarray, cur_node: int, end_node: int, r: int,
              RN: np.ndarray, current_time: float,
              stats: List[int]) -> Tuple[np.ndarray, List[int], np.ndarray]:
        """
        Find all paths from cur_node to end_node and compute response times.

        Args:
            sn: Network structure
            P: Routing matrix
            cur_node: Current node
            end_node: End node
            r: Class index
            RN: Response time matrix (may be modified for nested forks)
            current_time: Accumulated time
            stats: List of visited stations

        Returns:
            Tuple of (response_times, stations, updated_RN)
        """
        if cur_node == end_node:
            return np.array([current_time]), stats, RN

        ri = []
        stat = []
        current_rn = 0

        if cur_node < len(sn.nodeToStation):
            station = sn.nodeToStation[cur_node]
            if not np.isnan(station):
                current_rn = RN[int(station), r]
                stats = stats + [int(station)]

        for next_node in np.where(P[cur_node, :] > 0)[0]:
            if sn.nodetype[next_node] == NodeType.FORK:
                join_idx_arr = np.where(sn.fj[next_node, :] > 0)[0]
                if len(join_idx_arr) > 0:
                    join_idx = join_idx_arr[0]
                    join_station = int(sn.nodeToStation[join_idx])

                    if RN[join_station, r] == 0:
                        ri1, stat1, RN = ModelAdapter.paths(
                            sn, P, next_node, join_idx, r, RN, 0, []
                        )
                        from ..api.fjnative import sn_join_quorum
                        d0 = ModelAdapter._compute_sync_delay(
                            ri1, sn_join_quorum(sn, int(join_idx), r, len(ri1)))
                        RN[join_station, r] = d0
                        for st in stat1:
                            RN[st, r] = 0

                    ri1, stat1, RN = ModelAdapter.paths(
                        sn, P, join_idx, end_node, r, RN,
                        current_time + current_rn, stats
                    )
            else:
                ri1, stat1, RN = ModelAdapter.paths(
                    sn, P, next_node, end_node, r, RN,
                    current_time + current_rn, stats
                )

            ri.extend(ri1)
            stat.extend(stat1)

        return np.array(ri), stat, RN

    @staticmethod
    def paths_cs(sn, orig_nodes: int, P: np.ndarray, cur_node: int,
                 end_node: int, cur_class: int, RN: np.ndarray,
                 current_time: float,
                 stats: List[int]) -> Tuple[np.ndarray, List[int], np.ndarray]:
        """
        Find paths with class switching.

        Similar to paths() but handles class switches.
        Uses 0-based indexing for Python.

        Args:
            sn: Network structure
            orig_nodes: Number of original nodes
            P: Combined routing matrix
            cur_node: Current node (0-based)
            end_node: End node (0-based)
            cur_class: Current class (0-based)
            RN: Response time matrix
            current_time: Accumulated time
            stats: List of visited stations

        Returns:
            Tuple of (response_times, stations, updated_RN)
        """
        if cur_node == end_node:
            return np.array([current_time]), stats, RN

        ri = []
        stat = []
        current_rn = 0

        if cur_node < len(sn.nodeToStation):
            station = sn.nodeToStation[cur_node]
            if not np.isnan(station):
                current_rn = RN[int(station), cur_class]
                stats = stats + [int(station)]

        # Use 0-based indexing: row = cur_class * orig_nodes + cur_node
        row_idx = cur_class * orig_nodes + cur_node
        if row_idx < 0 or row_idx >= P.shape[0]:
            return np.array([current_time]), stats, RN

        for transition in np.where(P[row_idx, :] > 0)[0]:
            # Decode using 0-based indexing
            next_class = transition // orig_nodes
            next_node = transition % orig_nodes

            if sn.nodetype[next_node] == NodeType.FORK:
                join_idx_arr = np.where(sn.fj[next_node, :] > 0)[0]
                if len(join_idx_arr) > 0:
                    join_idx = join_idx_arr[0]
                    join_station = int(sn.nodeToStation[join_idx])

                    if RN[join_station, next_class] == 0:
                        ri1, stat1, RN = ModelAdapter.paths_cs(
                            sn, orig_nodes, P, next_node, join_idx,
                            next_class, RN, 0, []
                        )
                        d0 = ModelAdapter._compute_sync_delay(ri1)
                        RN[join_station, next_class] = d0
                        for st in stat1:
                            RN[st, next_class] = 0

                    ri1, stat1, RN = ModelAdapter.paths_cs(
                        sn, orig_nodes, P, join_idx, end_node,
                        next_class, RN, current_time + current_rn, stats
                    )
            else:
                ri1, stat1, RN = ModelAdapter.paths_cs(
                    sn, orig_nodes, P, next_node, end_node,
                    next_class, RN, current_time + current_rn, stats
                )

            ri.extend(ri1)
            stat.extend(stat1)

        return np.array(ri), stat, RN

    @staticmethod
    def _compute_sync_delay(path_times: np.ndarray, k: Optional[int] = None) -> float:
        """
        Compute the instant the join fires, by the Heidelberger-Trivedi formula.

        For K parallel branches with response times r_1, ..., r_K, the join fires at the
        k-th of them: the expected maximum E[max(r_1, ..., r_K)] on a standard join, and
        the k-th order statistic under a quorum. Both come from fj_ordstat_exp, whose
        inclusion-exclusion sum reduces to the classical one at k = K.

        Args:
            path_times: Array of response times for parallel paths
            k: quorum; None or K is the standard join

        Returns:
            Expected synchronization point time
        """
        if len(path_times) == 0:
            return 0.0
        if len(path_times) == 1:
            return path_times[0]

        path_times = np.asarray(path_times)
        # Avoid division by zero
        path_times = np.maximum(path_times, 1e-10)
        from ..api.fjnative import fj_ordstat_exp
        return fj_ordstat_exp(path_times, len(path_times) if k is None else k)

    @staticmethod
    def _compute_node_visits(sn) -> np.ndarray:
        """Compute total node visits across all chains."""
        Vnodes = np.zeros((sn.nnodes, sn.nclasses))

        # Try nodevisits first
        valid_nodevisits = False
        if sn.nodevisits is not None and len(sn.nodevisits) > 0:
            for chain_idx in range(len(sn.nodevisits)):
                nv = sn.nodevisits[chain_idx]
                # Check if nodevisits entry is a proper matrix
                if isinstance(nv, np.ndarray) and nv.shape == (sn.nnodes, sn.nclasses):
                    Vnodes += nv
                    valid_nodevisits = True
            if valid_nodevisits:
                return Vnodes

        # Fall back: assume all non-source/sink nodes are visited by all classes
        # This is safe for closed networks where all nodes are on the cycle
        if hasattr(sn, 'nodetype'):
            from ..api.sn.network_struct import NodeType
            for i in range(sn.nnodes):
                nt = sn.nodetype[i]
                nt_val = nt.value if hasattr(nt, 'value') else nt
                if nt_val not in (NodeType.SOURCE.value, NodeType.SINK.value):
                    Vnodes[i, :] = 1.0

        return Vnodes

    @staticmethod
    def _ensure_p_size(P: List[List[np.ndarray]], nclasses: int,
                       nnodes: int) -> List[List[np.ndarray]]:
        """Ensure P has enough dimensions for all classes and nodes."""
        while len(P) < nclasses:
            P.append([])
        for r in range(nclasses):
            while len(P[r]) < nclasses:
                P[r].append(np.zeros((nnodes, nnodes)))
            for s in range(nclasses):
                P[r][s] = ModelAdapter._ensure_matrix_size(P[r][s], nnodes)
        return P

    @staticmethod
    def _ensure_matrix_size(M: np.ndarray, size: int) -> np.ndarray:
        """Ensure matrix is at least size x size."""
        if M.shape[0] >= size and M.shape[1] >= size:
            return M
        new_M = np.zeros((size, size))
        new_M[:M.shape[0], :M.shape[1]] = M
        return new_M

    @staticmethod
    def _find_node_indices(nodetype_list, target_type) -> np.ndarray:
        """
        Find indices of nodes with a specific type.

        This method properly handles enum comparisons with lists of NodeType values,
        avoiding issues with numpy's broadcasting on enum comparisons.

        Args:
            nodetype_list: List or array of NodeType values
            target_type: The NodeType to search for

        Returns:
            NumPy array of indices where nodetype matches target_type
        """
        indices = []
        target_val = target_type.value if hasattr(target_type, 'value') else target_type
        for i, nt in enumerate(nodetype_list):
            nt_val = nt.value if hasattr(nt, 'value') else nt
            if nt_val == target_val:
                indices.append(i)
        return np.array(indices, dtype=int)

    @staticmethod
    def _get_fork_branches(model, fork_idx: int, class_idx: int) -> list:
        """
        Get the branch destinations from a Fork node for a given class.

        Args:
            model: The Network model
            fork_idx: Index of the Fork node
            class_idx: Index of the job class

        Returns:
            List of destination node indices
        """
        rm = model.get_routing_matrix()
        if rm is None:
            return []

        fork_node = model._nodes[fork_idx]
        job_class = model._classes[class_idx]

        # Get destinations from routing matrix
        branches = []
        routes = rm._routes.get((job_class, job_class), {})
        for (from_node, to_node), prob in routes.items():
            if from_node == fork_node and prob > 0:
                branches.append(to_node)

        return branches

    @staticmethod
    def aggregate_fes(model, station_subset, options=None):
        """
        Replace a station subset with a Flow-Equivalent Server (FES).

        Replaces a subset of stations in a closed product-form queueing network
        with a single FES station. The FES has a class-dependence handle (cdscaling)
        service rates where the rate for class-c in state (n1,...,nK) equals the
        throughput of class-c in an isolated subnetwork consisting only of the
        subset stations.

        Args:
            model: Closed product-form Network model.
            station_subset: List of Station objects to aggregate.
            options (dict, optional): Options dict with keys:
                - solver (str): Solver for throughput computation (default: 'mva')
                - cutoffs (numpy.ndarray): Per-class population cutoffs
                - verbose (bool): Enable verbose output

        Returns:
            dict with keys:
                - fes_model: New Network with the FES station replacing the subset.
                - fes_station: Reference to the FES Queue station in the new model.
                - deagg_info (dict): Deaggregation information.
        """
        from ..api.fes import fes_aggregate, fes_aggregate_with_options

        if options is not None:
            solver = options.get('solver', 'mva')
            cutoffs = options.get('cutoffs', None)
            verbose = options.get('verbose', False)
            return fes_aggregate_with_options(model, station_subset, solver=solver,
                                              cutoffs=cutoffs, verbose=verbose)
        else:
            return fes_aggregate(model, station_subset)


    @staticmethod
    def tag_chain(model: Any, chain: Any, jobclass: Any = None,
                  suffix: str = '.tagged') -> TaggedModelResult:
        """
        Create a tagged job model for response time analysis.

        One job of `jobclass` is MOVED out of its own class into a new class of
        population 1, so a solver can follow that single job through the network.

        THIS REPLACED A STUB. The previous body edited a NetworkStruct in place
        (nclasses, njobs, rt) on a deep-copied model and returned `tagged_job` as
        a plain dict. It created no JobClass object, so the native state-space
        generator -- which walks the model's class and node OBJECTS -- had
        nothing to walk; the struct said nclasses=2 while the model still held
        one class; and the population was never moved (njobs stayed [3,0] where
        it should read [2,1]). The result could not be solved at all, which is
        why SolverCTMC.getCdfRespT fell back to fitting an exponential.

        The three things that make a tagged class real: a service process at
        EVERY station cloned from the source class, the source class's routing
        replicated for the new class over the linked routing matrix, and one job
        actually MOVED rather than added.

        Args:
            model: Network model
            chain: Chain containing the class to tag
            jobclass: Specific class to tag (default: first class in chain)
            suffix: Suffix for tagged class names (default: '.tagged')

        Returns:
            TaggedModelResult with the tagged model and the tagged JobClass

        References:
            MATLAB: matlab/src/io/@ModelAdapter/tagChain.m
        """
        from ..lang.classes import ClosedClass

        chain_classes = list(getattr(chain, 'classes', []) or [])
        if jobclass is None and chain_classes:
            jobclass = chain_classes[0]
        if jobclass is None:
            raise ValueError("tag_chain: the chain carries no class to tag")
        if not chain_classes:
            chain_classes = [jobclass]

        tagged_model = model.copy()

        # Resolve the chain's classes inside the COPY: the chain object holds
        # references into the original model, which the copy does not share.
        by_name = {c.getName(): c for c in tagged_model.get_classes()}
        chain_in_tagged = [by_name[c.getName()] for c in chain_classes
                           if c.getName() in by_name]
        source = by_name.get(jobclass.getName())
        if source is None or not chain_in_tagged:
            raise ValueError("tag_chain: the class to tag is not present in the model copy")

        # The linked routing matrix must be read BEFORE the new classes exist,
        # so that it is the original class-pair map and not a half-extended one.
        P = tagged_model.get_linked_routing_matrix()
        if P is None:
            raise ValueError("tag_chain: the model has no linked routing matrix; "
                             "link() it before asking for a tagged copy")
        P = [[np.array(P[r][s], dtype=float, copy=True) for s in range(len(P[r]))]
             for r in range(len(P))]

        nnodes = len(tagged_model.get_nodes())
        stations = tagged_model.get_stations()
        old_index = {c.getName(): i for i, c in enumerate(tagged_model.get_classes())}

        tagged_classes = []
        for src in chain_in_tagged:
            pop = 1 if src is source else 0
            refstat = getattr(src, '_refstat', None)
            if refstat is None and stations:
                refstat = stations[0]
            new_cls = ClosedClass(tagged_model, src.getName() + suffix, pop, refstat,
                                  getattr(src, '_prio', 0))
            # A class with no service process is served nowhere, and the
            # state-space generator then has no transition to build.
            for st in stations:
                if not hasattr(st, 'get_service'):
                    continue
                dist = st.get_service(src)
                if dist is not None:
                    st.set_service(new_cls, copy.deepcopy(dist))
            tagged_classes.append(new_cls)

        # MOVE the job: the tagged class gained one, so the source loses one.
        if hasattr(source, 'setNumberOfJobs'):
            source.setNumberOfJobs(source.getNumberOfJobs() - 1)

        # Replicate the chain's routing for the tagged classes, pair by pair.
        R = len(tagged_model.get_classes())
        Pnew = [[np.zeros((nnodes, nnodes)) for _ in range(R)] for _ in range(R)]
        for r in range(min(len(P), R)):
            for s_ in range(min(len(P[r]), R)):
                Pnew[r][s_] = np.array(P[r][s_], dtype=float, copy=True)
        for ir, src_r in enumerate(chain_in_tagged):
            nr = tagged_classes[ir]._index if hasattr(tagged_classes[ir], '_index') \
                else R - len(tagged_classes) + ir
            orr = old_index[src_r.getName()]
            for is_, src_s in enumerate(chain_in_tagged):
                ns = tagged_classes[is_]._index if hasattr(tagged_classes[is_], '_index') \
                    else R - len(tagged_classes) + is_
                oss = old_index[src_s.getName()]
                if orr < len(P) and oss < len(P[orr]):
                    Pnew[nr][ns] = np.array(P[orr][oss], dtype=float, copy=True)

        tagged_model.reset_network(False)
        tagged_model.link(Pnew)
        tagged_model.refresh_struct()

        return TaggedModelResult(model=tagged_model,
                                 tagged_job=tagged_classes[-1] if tagged_classes else jobclass)

    @staticmethod
    def aggregate_chains(model: Any, suffix: str = '') -> Tuple[Any, np.ndarray, DeaggregationInfo]:
        """
        Transform a multi-class model into an equivalent chain-aggregated model.

        Classes belonging to the same chain (i.e. classes that can switch into
        each other) are merged into one aggregate class, so the returned model
        has one class per chain and no class switching. The aggregate preserves
        the total chain population (closed chains), the total arrival rate (open
        chains), the chain service demands and the chain-level routing.

        Args:
            model: Source Network model with potentially multiple classes per chain
            suffix: Optional suffix for chain class names (default: '')

        Returns:
            Tuple of (chain_model, alpha, deagg_info), where alpha is the
            (M, K) matrix of aggregation factors and deagg_info carries the
            data needed by sn_deaggregate_chain_results to map chain-level
            metrics back to class-level metrics.

        References:
            MATLAB: matlab/src/io/@ModelAdapter/aggregateChains.m
            JAR: jline.lang.ModelAdapter.aggregateChains
        """
        from ..lang.network import Network
        from ..lang.base import NodeType
        from ..lang.nodes import Source, Sink, Queue, Delay, Router, ClassSwitch
        from ..lang.classes import OpenClass, ClosedClass
        from ..distributions import Exp, Det, Erlang, HyperExp, Disabled
        from ..api.sn.demands import sn_get_demands_chain
        from ..constants import GlobalConstants

        if not hasattr(model, 'get_struct') and not hasattr(model, 'getStruct'):
            raise TypeError(
                "ModelAdapter.aggregate_chains requires a Network model, not a bare NetworkStruct.")
        sn = model.getStruct() if hasattr(model, 'getStruct') else model.get_struct()

        M = int(sn.nstations)
        K = int(sn.nclasses)
        C = int(sn.nchains)

        inchain = [np.asarray(sn.inchain[c], dtype=int).ravel() for c in range(C)]

        # If each class is its own chain there is nothing to merge
        if C == K:
            chain_model = model.copy()
            alpha = np.eye(M, K)
            deagg_info = DeaggregationInfo(
                alpha=alpha,
                inchain=[list(inchain[c]) for c in range(C)],
                original_sn=sn,
                is_aggregated=False,
                nclasses=K,
                nchains=C,
            )
            return chain_model, alpha, deagg_info

        demands = sn_get_demands_chain(sn)
        Lchain = np.asarray(demands.Lchain, dtype=float)
        STchain = np.asarray(demands.STchain, dtype=float)
        Vchain = np.asarray(demands.Vchain, dtype=float)
        alpha = np.asarray(demands.alpha, dtype=float)
        Nchain = np.asarray(demands.Nchain, dtype=float).ravel()
        SCVchain = np.asarray(demands.SCVchain, dtype=float)
        refstatchain = np.asarray(demands.refstatchain, dtype=int).ravel()

        # Open versus closed chains, and the aggregate arrival rate of each
        # open chain
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        is_open_chain = np.zeros(C, dtype=bool)
        lambda_chain = np.zeros(C)
        source_station = None
        for i in range(int(sn.nnodes)):
            if int(sn.nodetype[i]) == int(NodeType.SOURCE):
                source_station = int(sn.nodeToStation[i])
                break
        for c in range(C):
            is_open_chain[c] = bool(np.any(np.isinf(njobs[inchain[c]])))
            if is_open_chain[c] and source_station is not None and source_station >= 0:
                rates = np.asarray(sn.rates, dtype=float)
                lambda_chain[c] = float(np.nansum(rates[source_station, inchain[c]]))

        # Rebuild the model with one class per chain
        chain_model = Network(f"{model.getName() if hasattr(model, 'getName') else model.name}_aggregated")

        node_map = {}
        station_map = {}
        for i, node in enumerate(model.getNodes()):
            # auto-added class switch nodes disappear with class switching, and
            # user-defined ones have nothing left to switch
            if isinstance(node, ClassSwitch):
                continue
            if isinstance(node, Source):
                node_map[i] = Source(chain_model, node.name)
            elif isinstance(node, Sink):
                node_map[i] = Sink(chain_model, node.name)
            elif isinstance(node, Delay):
                node_map[i] = Delay(chain_model, node.name)
            elif isinstance(node, Queue):
                new_node = Queue(chain_model, node.name, node.get_sched_strategy())
                nservers = node.get_number_of_servers()
                if np.isfinite(nservers):
                    new_node.setNumberOfServers(int(nservers))
                cap = node.get_capacity()
                if cap is not None and np.isfinite(cap):
                    new_node.setCapacity(int(cap))
                node_map[i] = new_node
            elif isinstance(node, Router):
                node_map[i] = Router(chain_model, node.name)
            else:
                line_warning('aggregate_chains',
                             f"Node type {type(node).__name__} not fully supported in chain aggregation.")
                continue
            ist = int(sn.nodeToStation[i]) if sn.nodeToStation is not None else -1
            if ist >= 0:
                station_map[ist] = node_map[i]

        # One aggregate class per chain
        chain_class = []
        for c in range(C):
            names = [str(sn.classnames[k]) for k in inchain[c]]
            name = names[0] if len(names) == 1 else f"Chain{c + 1}"
            if suffix:
                name = name + suffix
            if is_open_chain[c]:
                chain_class.append(OpenClass(chain_model, name))
            else:
                refstat = station_map.get(int(refstatchain[c]))
                if refstat is None:
                    raise RuntimeError(
                        f"Reference station {int(refstatchain[c])} for chain {c} not found in aggregated model.")
                chain_class.append(ClosedClass(chain_model, name, float(Nchain[c]), refstat))

        if np.any(is_open_chain):
            chain_source = chain_model.getSource()
            for c in range(C):
                if is_open_chain[c] and lambda_chain[c] > 0:
                    chain_source.setArrival(chain_class[c], Exp(lambda_chain[c]))

        # Aggregate service processes, matched on the chain mean and SCV
        for ist in range(M):
            station = station_map.get(ist)
            if station is None or isinstance(station, (Source, Sink)):
                continue
            for c in range(C):
                st = STchain[ist, c]
                if st > 0 and np.isfinite(st):
                    scv = SCVchain[ist, c]
                    if not np.isfinite(scv) or scv <= 0:
                        scv = 1.0
                    if abs(scv - 1.0) < GlobalConstants.FineTol:
                        dist = Exp(1.0 / st)
                    elif scv < 1.0:
                        if scv < GlobalConstants.FineTol:
                            dist = Det(st)
                        else:
                            dist = Erlang.fitMeanAndOrder(st, max(1, int(round(1.0 / scv))))
                    else:
                        dist = HyperExp.fitMeanAndSCV(st, scv)
                    station.setService(chain_class[c], dist)
                else:
                    station.setService(chain_class[c], Disabled())

        # Chain-level routing: aggregate the class-level routing probabilities
        # with the alpha weights, then renormalise
        P = chain_model.initRoutingMatrix()
        rt = np.asarray(sn.rt, dtype=float)
        station_to_stateful = sn.stationToStateful if getattr(sn, 'stationToStateful', None) is not None else None
        for c in range(C):
            rows = {}
            for ist in range(M):
                src = station_map.get(ist)
                if src is None:
                    continue
                if not is_open_chain[c] and isinstance(src, (Source, Sink)):
                    continue
                isf_i = int(station_to_stateful[ist]) if station_to_stateful is not None else ist
                for jst in range(M):
                    dst = station_map.get(jst)
                    if dst is None:
                        continue
                    isf_j = int(station_to_stateful[jst]) if station_to_stateful is not None else jst
                    pij = 0.0
                    for k in inchain[c]:
                        for s in inchain[c]:
                            from_idx = isf_i * K + int(k)
                            to_idx = isf_j * K + int(s)
                            if from_idx < rt.shape[0] and to_idx < rt.shape[1]:
                                p_ks = rt[from_idx, to_idx]
                                if alpha[ist, int(k)] > 0 and p_ks > 0:
                                    pij += alpha[ist, int(k)] * p_ks
                    if pij > GlobalConstants.FineTol:
                        rows.setdefault(src, []).append((dst, pij))
            for src, entries in rows.items():
                total = sum(p for _, p in entries)
                if total > GlobalConstants.FineTol:
                    if abs(total - 1.0) > 0.01:
                        line_warning('aggregate_chains',
                                     f"Large normalization correction at node {src.name} for chain {c}: rowSum={total:.4f}")
                    for dst, p in entries:
                        P.set(chain_class[c], chain_class[c], src, dst, p / total)

        chain_model.link(P)

        deagg_info = DeaggregationInfo(
            alpha=alpha,
            inchain=[list(inchain[c]) for c in range(C)],
            original_sn=sn,
            is_aggregated=True,
            V_chain=Vchain,
            ST_chain=STchain,
            L_chain=Lchain,
            SCV_chain=SCVchain,
            N_chain=Nchain,
            lambda_chain=lambda_chain,
            is_open_chain=is_open_chain,
            refstat=np.asarray(sn.refstat).ravel() if getattr(sn, 'refstat', None) is not None else None,
            refstat_chain=refstatchain,
            nclasses=K,
            nchains=C,
        )

        return chain_model, alpha, deagg_info

    @staticmethod
    def remove_class(model: Any, jobclass: Any) -> Any:
        """
        Return a copy of the model with a job class removed.

        The original model is left untouched, so it stays solvable; use this for
        ablation studies or a per-class decomposition. The removal logic itself
        lives in Network.remove_class and is delegated to here rather than
        duplicated, so that the two entry points cannot drift. Mirrors
        ModelAdapter.removeClass in MATLAB and the JAR.

        Args:
            model: Network model
            jobclass: Job class to remove (class object or 0-based index)

        Returns:
            A copy of the model without the specified class

        References:
            MATLAB: matlab/src/io/@ModelAdapter/removeClass.m
            JAR: jline.lang.ModelAdapter.removeClass
        """
        if not hasattr(model, 'remove_class'):
            raise TypeError(
                "ModelAdapter.remove_class requires a Network model, not a bare NetworkStruct.")
        new_model = model.copy() if hasattr(model, 'copy') else copy.deepcopy(model)
        new_model.remove_class(jobclass)
        return new_model

    # Alias for Java camelCase compatibility
    aggregateFES = aggregate_fes


# Convenience function aliases
def tag_chain(model: Any, chain: Any, jobclass: Any = None,
              suffix: str = '.tagged') -> TaggedModelResult:
    """Create a tagged job model. See ModelAdapter.tag_chain."""
    return ModelAdapter.tag_chain(model, chain, jobclass, suffix)


def aggregate_chains(model: Any, suffix: str = '') -> Tuple[Any, np.ndarray, DeaggregationInfo]:
    """Aggregate chains in a model. See ModelAdapter.aggregate_chains."""
    return ModelAdapter.aggregate_chains(model, suffix)


def remove_class(model: Any, jobclass: Any) -> Any:
    """Remove a class from model. See ModelAdapter.remove_class."""
    return ModelAdapter.remove_class(model, jobclass)


def sort_forks(sn: Any, fjforkmap: np.ndarray, fjclassmap: np.ndarray,
               nonfjmodel: Any) -> Tuple[np.ndarray, np.ndarray]:
    """Sort forks in fork-join network. See ModelAdapter.sort_forks."""
    return ModelAdapter.sort_forks(sn, fjforkmap, fjclassmap, nonfjmodel)


def mmt(model: Any, fork_lambda: Optional[np.ndarray] = None) -> MMTResult:
    """MMT transform for fork-join networks. See ModelAdapter.mmt."""
    return ModelAdapter.mmt(model, fork_lambda)


def ht(model: Any) -> HTResult:
    """HT transform for fork-join networks. See ModelAdapter.ht."""
    return ModelAdapter.ht(model)


def paths(sn: Any, P: np.ndarray, cur_node: int, end_node: int,
          r: int, RN: np.ndarray, current_time: float,
          stats: List[int]) -> Tuple[List[float], List[int], np.ndarray]:
    """Find paths in network. See ModelAdapter.paths."""
    return ModelAdapter.paths(sn, P, cur_node, end_node, r, RN, current_time, stats)


def paths_cs(sn: Any, orignodes: int, P: np.ndarray, cur_node: int,
             end_node: int, cur_class: int, RN: np.ndarray,
             current_time: float, stats: List[int]) -> Tuple[List[float], List[int], np.ndarray]:
    """Find paths with class switching. See ModelAdapter.paths_cs."""
    return ModelAdapter.paths_cs(sn, orignodes, P, cur_node, end_node, cur_class, RN, current_time, stats)


def find_paths(sn: Any, P: np.ndarray, start: int, end_node: int,
               r: int, to_merge: List[int], QN: np.ndarray, TN: np.ndarray,
               current_time: float, fjclassmap: np.ndarray,
               fjforkmap: np.ndarray, nonfjmodel: Any) -> List[float]:
    """Find paths for fork-join analysis. See ModelAdapter.find_paths."""
    return ModelAdapter.find_paths(sn, P, start, end_node, r, to_merge, QN, TN,
                                   current_time, fjclassmap, fjforkmap, nonfjmodel)


def find_paths_cs(sn: Any, P: np.ndarray, cur_node: int, end_node: int,
                  cur_class: int, to_merge: List[int], QN: np.ndarray,
                  TN: np.ndarray, current_time: float, fjclassmap: np.ndarray,
                  fjforkmap: np.ndarray, nonfjmodel: Any) -> List[float]:
    """Find paths with class switching for fork-join. See ModelAdapter.find_paths_cs."""
    return ModelAdapter.find_paths_cs(sn, P, cur_node, end_node, cur_class, to_merge,
                                       QN, TN, current_time, fjclassmap, fjforkmap, nonfjmodel)


def fjtag(model):
    """Fork-join state-space tagging. See ModelAdapter.fjtag."""
    return ModelAdapter.fjtag(model)


def aggregate_fes(model, station_subset, options=None):
    """Flow-equivalent server aggregation. See ModelAdapter.aggregate_fes."""
    return ModelAdapter.aggregate_fes(model, station_subset, options)


__all__ = [
    'ModelAdapter',
    'DeaggregationInfo',
    'TaggedModelResult',
    'FESAggregationInfo',
    'MMTResult',
    'HTResult',
    'tag_chain',
    'aggregate_chains',
    'remove_class',
    'sort_forks',
    'mmt',
    'ht',
    'paths',
    'paths_cs',
    'find_paths',
    'find_paths_cs',
    'fjtag',
    'aggregate_fes',
]
