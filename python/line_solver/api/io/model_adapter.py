"""
Model adaptation and transformation utilities.

This module provides the ModelAdapter class for transforming and adapting
queueing network models, including:
- Creating tagged job models for response time analysis
- Chain aggregation for multi-class models
- Fork-join network transformations (MMT, HT)
- Path finding in networks
- Model preprocessing operations

Port from:
    - matlab/src/io/@ModelAdapter/ModelAdapter.m
    - matlab/src/io/@ModelAdapter/tagChain.m
    - matlab/src/io/@ModelAdapter/aggregateChains.m
    - matlab/src/io/@ModelAdapter/removeClass.m
    - matlab/src/io/@ModelAdapter/mmt.m
    - matlab/src/io/@ModelAdapter/ht.m
    - matlab/src/io/@ModelAdapter/sortForks.m
    - matlab/src/io/@ModelAdapter/paths.m
    - matlab/src/io/@ModelAdapter/pathsCS.m
    - matlab/src/io/@ModelAdapter/findPaths.m
    - matlab/src/io/@ModelAdapter/findPathsCS.m
"""

import copy
import numpy as np
from typing import Any, Optional, Dict, List, Tuple
from dataclasses import dataclass, field

from .logging import line_warning, line_error


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


@dataclass
class MMTResult:
    """Result from MMT (Markov Modulated Transform) for Fork-Join networks."""
    nonfjmodel: Any  # Model with fork-joins replaced
    fjclassmap: np.ndarray  # Maps auxiliary class to original class
    fjforkmap: np.ndarray  # Maps auxiliary class to fork node
    fanout: np.ndarray  # Fanout per auxiliary class


@dataclass
class HTResult:
    """Result from HT (Heidelberger-Trivedi) transform for Fork-Join networks."""
    nonfjmodel: Any  # Model with fork-joins replaced
    fjclassmap: np.ndarray  # Maps auxiliary class to original class
    fjforkmap: np.ndarray  # Maps auxiliary class to fork node
    fj_auxiliary_delays: Dict[int, int]  # Maps join to auxiliary delay


class ModelAdapter:
    """
    Static class for model transformations and adaptations.

    Provides functionality for:
    - Creating tagged job models for response time analysis
    - Chain aggregation (merging classes in the same chain)
    - Fork-join network transformations (MMT, HT)
    - Flow-Equivalent Server (FES) aggregation
    - Path finding in networks
    - Model preprocessing operations

    All methods are static and do not require instantiation.

    References:
        MATLAB: matlab/src/io/@ModelAdapter/ModelAdapter.m
    """

    @staticmethod
    def tag_chain(model: Any, chain: Any, jobclass: Any = None,
                  suffix: str = '.tagged') -> TaggedModelResult:
        """
        Create a tagged job model for response time analysis.

        Args:
            model: Network model
            chain: Chain containing the class to tag
            jobclass: Specific class to tag (default: first class in chain)
            suffix: Suffix for tagged class names (default: '.tagged')

        Returns:
            TaggedModelResult containing tagged model and tagged job class

        References:
            MATLAB: matlab/src/io/@ModelAdapter/tagChain.m
        """
        if jobclass is None and hasattr(chain, 'classes') and len(chain.classes) > 0:
            jobclass = chain.classes[0]

        tagged_model = copy.deepcopy(model)

        if hasattr(tagged_model, 'getStruct'):
            sn = tagged_model.getStruct()
        else:
            sn = tagged_model

        if hasattr(chain, 'index'):
            chain_indexes = chain.index if isinstance(chain.index, list) else [chain.index]
        else:
            chain_indexes = [0]

        K = sn.nclasses
        R = K

        tagged_classes = []
        for r in chain_indexes:
            new_class_idx = K + len(tagged_classes)

            if hasattr(sn, 'classnames'):
                old_name = sn.classnames[r]
                new_name = f"{old_name}{suffix}"
                sn.classnames = list(sn.classnames) + [new_name]

            if hasattr(sn, 'njobs'):
                sn.njobs = np.append(sn.njobs, 0)
                if hasattr(jobclass, 'index') and r == jobclass.index:
                    sn.njobs[new_class_idx] = 1
                    sn.njobs[r] -= 1

            tagged_classes.append(new_class_idx)

        sn.nclasses = K + len(tagged_classes)

        if hasattr(sn, 'rt') and sn.rt is not None:
            old_rt = sn.rt
            new_size = sn.nstations * sn.nclasses
            new_rt = np.zeros((new_size, new_size))
            old_size = sn.nstations * K
            new_rt[:old_size, :old_size] = old_rt

            for ir, r in enumerate(chain_indexes):
                for js, s in enumerate(chain_indexes):
                    for ist in range(sn.nstations):
                        for jst in range(sn.nstations):
                            old_from = ist * K + r
                            old_to = jst * K + s
                            new_from = ist * sn.nclasses + (R + ir)
                            new_to = jst * sn.nclasses + (R + js)
                            if old_from < old_rt.shape[0] and old_to < old_rt.shape[1]:
                                new_rt[new_from, new_to] = old_rt[old_from, old_to]

            sn.rt = new_rt

        tagged_job = {
            'class_indices': tagged_classes,
            'original_class': jobclass.index if hasattr(jobclass, 'index') else 0,
            'suffix': suffix,
        }

        return TaggedModelResult(model=tagged_model, tagged_job=tagged_job)

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
        from ...lang.network import Network
        from ...lang.base import NodeType
        from ...lang.nodes import Source, Sink, Queue, Delay, Router, ClassSwitch
        from ...lang.classes import OpenClass, ClosedClass
        from ...distributions import Exp, Det, Erlang, HyperExp, Disabled
        from ...api.sn.demands import sn_get_demands_chain
        from ...constants import GlobalConstants

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

    @staticmethod
    def sort_forks(sn: Any, fjforkmap: np.ndarray, fjclassmap: np.ndarray,
                   nonfjmodel: Any) -> Tuple[np.ndarray, np.ndarray]:
        """
        Sort fork nodes based on nesting structure.

        Args:
            sn: NetworkStruct
            fjforkmap: Map from auxiliary class to fork node index
            fjclassmap: Map from auxiliary class to original class index
            nonfjmodel: Model without fork-join (from mmt or ht)

        Returns:
            Tuple of (forks, parents)

        References:
            MATLAB: matlab/src/io/@ModelAdapter/sortForks.m
        """
        N = len(sn.nodetype)
        max_class = int(np.max(fjclassmap)) + 1 if len(fjclassmap) > 0 else 1

        forks = np.zeros((N, max_class))
        parents = np.zeros(N, dtype=int)

        fork_indices = []
        for i, nt in enumerate(sn.nodetype):
            if hasattr(nt, 'name') and nt.name == 'FORK':
                fork_indices.append(i)
                forks[i, :] = 1
                parents[i] = i

        if hasattr(nonfjmodel, 'getLinkedRoutingMatrix'):
            P = nonfjmodel.getLinkedRoutingMatrix()
        else:
            P = None

        def nested_forks(start_node: int, end_node: int, conn: np.ndarray,
                         fork_flags: np.ndarray) -> np.ndarray:
            if start_node == end_node:
                return fork_flags

            for i in np.where(conn[start_node, :] > 0)[0]:
                nt = sn.nodetype[i]
                if hasattr(nt, 'name') and nt.name == 'FORK':
                    fork_flags[i] = 0
                fork_flags = fork_flags & nested_forks(i, end_node, conn, fork_flags)

            return fork_flags

        for f in fork_indices:
            if hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[f, :])[0]
                if len(join_indices) == 0:
                    continue
                join_idx = join_indices[0]

                fork_aux_classes = np.where(fjforkmap == f)[0]

                for s in fork_aux_classes:
                    if s >= len(fjclassmap):
                        continue
                    r = int(fjclassmap[s])

                    if P is not None and r < len(P) and r < P.shape[0]:
                        conn = P[r, r] if isinstance(P, dict) else P
                        if isinstance(conn, np.ndarray):
                            discovered_forks = np.zeros(N)
                            for i in fork_indices:
                                discovered_forks[i] = 1

                            nested = nested_forks(f, join_idx, conn, discovered_forks)
                            forks[:, r] = forks[:, r] & nested

                            for i in np.where(nested == 0)[0]:
                                parents[i] = parents[f]

        return forks, parents

    @staticmethod
    def mmt(model: Any, fork_lambda: Optional[np.ndarray] = None) -> MMTResult:
        """
        MMT (Markov Modulated Transform) for Fork-Join networks.

        Args:
            model: Network model with fork-join structures
            fork_lambda: Optional arrival rates for artificial classes

        Returns:
            MMTResult containing transformed model and class mappings

        References:
            MATLAB: matlab/src/io/@ModelAdapter/mmt.m
        """
        sn = model.getStruct() if hasattr(model, 'getStruct') else model

        if fork_lambda is None:
            fork_lambda = 1e-6 * np.ones(sn.nclasses)

        fjclassmap = []
        fjforkmap = []
        fanout = []

        nonfjmodel = copy.deepcopy(model)

        if hasattr(nonfjmodel, 'getLinkedRoutingMatrix'):
            P = nonfjmodel.getLinkedRoutingMatrix()
        else:
            line_error('mmt', 'Model must have linked routing matrix')
            return MMTResult(nonfjmodel, np.array([]), np.array([]), np.array([]))

        if hasattr(sn, 'nodevisits') and sn.nodevisits is not None:
            Vnodes = sum(sn.nodevisits)
        else:
            Vnodes = np.ones((sn.nnodes, sn.nclasses))

        forked_classes = {}
        fork_indices = []

        for i, nt in enumerate(sn.nodetype):
            if hasattr(nt, 'name') and nt.name == 'FORK':
                fork_indices.append(i)

        origfanout = np.zeros((sn.nnodes, sn.nclasses))

        for f in fork_indices:
            for r in range(sn.nclasses):
                if hasattr(model, 'nodes') and f < len(model.nodes):
                    node = model.nodes[f]
                    if hasattr(node, 'output') and hasattr(node.output, 'outputStrategy'):
                        if r < len(node.output.outputStrategy):
                            strategy = node.output.outputStrategy[r]
                            if len(strategy) > 2 and strategy[2] is not None:
                                origfanout[f, r] = len(strategy[2])
                                if P is not None:
                                    for s in range(sn.nclasses):
                                        if isinstance(P, dict) and (r, s) in P:
                                            P[(r, s)][f, :] = P[(r, s)][f, :] / origfanout[f, r]

            forked_classes[f] = np.where(Vnodes[f, :] > 0)[0].tolist()

        for f in fork_indices:
            if hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[f, :])[0]
                if len(join_indices) == 0:
                    continue

                if hasattr(sn, 'chains') and sn.chains is not None:
                    forked_chains = np.where(np.sum(sn.chains[:, forked_classes.get(f, [])], axis=1))[0]
                else:
                    forked_chains = [0]

                for fc in forked_chains:
                    if hasattr(sn, 'chains') and sn.chains is not None:
                        inchain = np.where(sn.chains[fc, :])[0].tolist()
                    else:
                        inchain = list(range(sn.nclasses))

                    for r in inchain:
                        fjclassmap.append(r)
                        fjforkmap.append(f)

                        if hasattr(model, 'nodes') and f < len(model.nodes):
                            node = model.nodes[f]
                            tasks_per_link = 1
                            if hasattr(node, 'output') and hasattr(node.output, 'tasksPerLink'):
                                tasks_per_link = node.output.tasksPerLink
                                if tasks_per_link > 1:
                                    line_warning('mmt', 'No synchronisation delays implemented in MMT for multiple tasks per link')
                            fanout.append(int(origfanout[f, r] * tasks_per_link))
                        else:
                            fanout.append(1)

        return MMTResult(
            nonfjmodel=nonfjmodel,
            fjclassmap=np.array(fjclassmap, dtype=int),
            fjforkmap=np.array(fjforkmap, dtype=int),
            fanout=np.array(fanout, dtype=int)
        )

    @staticmethod
    def ht(model: Any) -> HTResult:
        """
        HT (Heidelberger-Trivedi) transform for Fork-Join networks.

        This method is for closed networks only.

        Args:
            model: Network model with fork-join structures

        Returns:
            HTResult containing transformed model and class mappings

        References:
            MATLAB: matlab/src/io/@ModelAdapter/ht.m
        """
        sn = model.getStruct() if hasattr(model, 'getStruct') else model

        fjclassmap = []
        fjforkmap = []
        fj_auxiliary_delays = {}

        nonfjmodel = copy.deepcopy(model)

        if hasattr(sn, 'nodevisits') and sn.nodevisits is not None:
            Vnodes = sum(sn.nodevisits)
        else:
            Vnodes = np.ones((sn.nnodes, sn.nclasses))

        forked_classes = {}
        fork_indices = []

        for i, nt in enumerate(sn.nodetype):
            if hasattr(nt, 'name') and nt.name == 'FORK':
                fork_indices.append(i)
                forked_classes[i] = np.where(Vnodes[i, :] > 0)[0].tolist()

        for i, nt in enumerate(sn.nodetype):
            if hasattr(nt, 'name') and nt.name == 'JOIN':
                aux_delay_idx = sn.nnodes + len(fj_auxiliary_delays)
                fj_auxiliary_delays[i] = aux_delay_idx

        for f in fork_indices:
            if hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[f, :])[0]
                if len(join_indices) == 0:
                    continue

                if hasattr(sn, 'chains') and sn.chains is not None:
                    forked_chains = np.where(np.sum(sn.chains[:, forked_classes.get(f, [])], axis=1))[0]
                else:
                    forked_chains = [0]

                for fc in forked_chains:
                    if hasattr(sn, 'chains') and sn.chains is not None:
                        inchain = np.where(sn.chains[fc, :])[0].tolist()
                    else:
                        inchain = list(range(sn.nclasses))

                    for r in inchain:
                        if hasattr(sn, 'nodevisits') and fc < len(sn.nodevisits):
                            if sn.nodevisits[fc][f, r] == 0:
                                continue

                        parallel_branches = 1
                        if hasattr(model, 'nodes') and f < len(model.nodes):
                            node = model.nodes[f]
                            if hasattr(node, 'output') and hasattr(node.output, 'outputStrategy'):
                                if r < len(node.output.outputStrategy):
                                    strategy = node.output.outputStrategy[r]
                                    if len(strategy) > 2 and strategy[2] is not None:
                                        parallel_branches = len(strategy[2])

                            if hasattr(node, 'output') and hasattr(node.output, 'tasksPerLink'):
                                if node.output.tasksPerLink > 1:
                                    line_error('ht', 'Multiple tasks per link are not supported in H-T')

                            if hasattr(model, 'classes') and r < len(model.classes):
                                cls = model.classes[r]
                                if hasattr(cls, '__class__') and cls.__class__.__name__ == 'OpenClass':
                                    line_error('ht', 'H-T method can be used only on closed models')

                        for par in range(parallel_branches):
                            fjclassmap.append(r)
                            fjforkmap.append(f)

        return HTResult(
            nonfjmodel=nonfjmodel,
            fjclassmap=np.array(fjclassmap, dtype=int),
            fjforkmap=np.array(fjforkmap, dtype=int),
            fj_auxiliary_delays=fj_auxiliary_delays
        )

    @staticmethod
    def paths(sn: Any, P: np.ndarray, cur_node: int, end_node: int,
              r: int, RN: np.ndarray, current_time: float,
              stats: List[int]) -> Tuple[List[float], List[int], np.ndarray]:
        """
        Find all paths from current node to end node and compute response times.

        Args:
            sn: NetworkStruct
            P: Routing probability matrix
            cur_node: Current node index
            end_node: Target node index
            r: Class index
            RN: Response time matrix (M x K)
            current_time: Accumulated response time so far
            stats: List of stations visited

        Returns:
            Tuple of (response_times, stations, updated_RN)

        References:
            MATLAB: matlab/src/io/@ModelAdapter/paths.m
        """
        from itertools import combinations

        if cur_node == end_node:
            return [current_time], stats, RN

        ri = []
        stat = []
        current_rn = 0

        if hasattr(sn, 'nodeToStation') and not np.isnan(sn.nodeToStation[cur_node]):
            ist = int(sn.nodeToStation[cur_node])
            current_rn = RN[ist, r]
            stats = stats + [ist]

        next_nodes = np.where(P[cur_node, :] > 0)[0]

        for next_node in next_nodes:
            nt = sn.nodetype[next_node]
            is_fork = hasattr(nt, 'name') and nt.name == 'FORK'

            if is_fork and hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[next_node, :])[0]
                if len(join_indices) > 0:
                    join_idx = join_indices[0]

                    ist_join = int(sn.nodeToStation[join_idx]) if hasattr(sn, 'nodeToStation') else join_idx
                    if RN[ist_join, r] == 0:
                        ri1, stat1, RN = ModelAdapter.paths(sn, P, next_node, join_idx, r, RN, 0, [])

                        if len(ri1) > 0 and all(x > 0 for x in ri1):
                            lambdai = [1.0 / x for x in ri1]
                            d0 = 0
                            parallel_branches = len(ri1)

                            for pow_val in range(parallel_branches):
                                current_sum = 0
                                for combo in combinations(lambdai, pow_val + 1):
                                    current_sum += 1.0 / sum(combo)
                                d0 += ((-1) ** pow_val) * current_sum

                            RN[ist_join, r] = d0
                            for s in stat1:
                                RN[s, r] = 0

                    ri1, stat1, RN = ModelAdapter.paths(sn, P, join_idx, end_node, r, RN,
                                                        current_time + current_rn, stats)
                    ri.extend(ri1)
                    stat.extend(stat1)
                else:
                    ri1, stat1, RN = ModelAdapter.paths(sn, P, next_node, end_node, r, RN,
                                                        current_time + current_rn, stats)
                    ri.extend(ri1)
                    stat.extend(stat1)
            else:
                ri1, stat1, RN = ModelAdapter.paths(sn, P, next_node, end_node, r, RN,
                                                    current_time + current_rn, stats)
                ri.extend(ri1)
                stat.extend(stat1)

        return ri, stat, RN

    @staticmethod
    def paths_cs(sn: Any, orignodes: int, P: np.ndarray, cur_node: int,
                 end_node: int, cur_class: int, RN: np.ndarray,
                 current_time: float, stats: List[int]) -> Tuple[List[float], List[int], np.ndarray]:
        """
        Find paths with class switching and compute response times.

        Args:
            sn: NetworkStruct
            orignodes: Number of original nodes
            P: Class-switching routing matrix (flattened)
            cur_node: Current node index
            end_node: Target node index
            cur_class: Current class index
            RN: Response time matrix
            current_time: Accumulated response time
            stats: Stations visited

        Returns:
            Tuple of (response_times, stations, updated_RN)

        References:
            MATLAB: matlab/src/io/@ModelAdapter/pathsCS.m
        """
        from itertools import combinations

        if cur_node == end_node:
            return [current_time], stats, RN

        ri = []
        stat = []
        current_rn = 0

        if hasattr(sn, 'nodeToStation') and not np.isnan(sn.nodeToStation[cur_node]):
            ist = int(sn.nodeToStation[cur_node])
            current_rn = RN[ist, cur_class]
            stats = stats + [ist]

        row_idx = cur_class * orignodes + cur_node
        transitions = np.where(P[row_idx, :] > 0)[0] if row_idx < P.shape[0] else []

        for transition in transitions:
            next_class = transition // orignodes
            next_node = transition % orignodes

            nt = sn.nodetype[next_node]
            is_fork = hasattr(nt, 'name') and nt.name == 'FORK'

            if is_fork and hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[next_node, :])[0]
                if len(join_indices) > 0:
                    join_idx = join_indices[0]

                    ist_join = int(sn.nodeToStation[join_idx]) if hasattr(sn, 'nodeToStation') else join_idx
                    if RN[ist_join, next_class] == 0:
                        ri1, stat1, RN = ModelAdapter.paths_cs(sn, orignodes, P, next_node,
                                                                join_idx, next_class, RN, 0, [])

                        if len(ri1) > 0 and all(x > 0 for x in ri1):
                            lambdai = [1.0 / x for x in ri1]
                            d0 = 0
                            parallel_branches = len(ri1)

                            for pow_val in range(parallel_branches):
                                current_sum = 0
                                for combo in combinations(lambdai, pow_val + 1):
                                    current_sum += 1.0 / sum(combo)
                                d0 += ((-1) ** pow_val) * current_sum

                            RN[ist_join, next_class] = d0
                            for s in stat1:
                                RN[s, next_class] = 0

                    ri1, stat1, RN = ModelAdapter.paths_cs(sn, orignodes, P, join_idx,
                                                            end_node, next_class, RN,
                                                            current_time + current_rn, stats)
                    ri.extend(ri1)
                    stat.extend(stat1)
                else:
                    ri1, stat1, RN = ModelAdapter.paths_cs(sn, orignodes, P, next_node,
                                                            end_node, next_class, RN,
                                                            current_time + current_rn, stats)
                    ri.extend(ri1)
                    stat.extend(stat1)
            else:
                ri1, stat1, RN = ModelAdapter.paths_cs(sn, orignodes, P, next_node,
                                                        end_node, next_class, RN,
                                                        current_time + current_rn, stats)
                ri.extend(ri1)
                stat.extend(stat1)

        return ri, stat, RN

    @staticmethod
    def find_paths(sn: Any, P: np.ndarray, start: int, end_node: int,
                   r: int, to_merge: List[int], QN: np.ndarray, TN: np.ndarray,
                   current_time: float, fjclassmap: np.ndarray,
                   fjforkmap: np.ndarray, nonfjmodel: Any) -> List[float]:
        """
        Find response times along paths for fork-join analysis.

        Args:
            sn: NetworkStruct
            P: Routing matrix
            start: Starting node
            end_node: Ending node
            r: Original class index
            to_merge: Classes to merge for queue length computation
            QN: Queue length matrix
            TN: Throughput matrix
            current_time: Accumulated time
            fjclassmap: Fork-join class map
            fjforkmap: Fork-join fork map
            nonfjmodel: Transformed model

        Returns:
            List of response times for each path

        References:
            MATLAB: matlab/src/io/@ModelAdapter/findPaths.m
        """
        from itertools import combinations

        if start == end_node:
            ist = int(sn.nodeToStation[start]) if hasattr(sn, 'nodeToStation') else start
            q_len = sum(QN[ist, k] for k in to_merge)
            tput = sum(TN[ist, k] for k in to_merge)
            if tput > 0:
                return [current_time - q_len / tput]
            return [current_time]

        ri = []
        next_nodes = np.where(P[start, :] > 0)[0]

        for i in next_nodes:
            q_len = 0
            tput = 1

            if hasattr(sn, 'nodeToStation') and not np.isnan(sn.nodeToStation[i]):
                ist = int(sn.nodeToStation[i])
                q_len = sum(QN[ist, k] for k in to_merge)
                tput = sum(TN[ist, k] for k in to_merge) if sum(TN[ist, k] for k in to_merge) > 0 else 1

            nt = sn.nodetype[i]
            is_fork = hasattr(nt, 'name') and nt.name == 'FORK'

            if is_fork:
                if hasattr(sn, 'fj') and sn.fj is not None:
                    join_indices = np.where(sn.fj[i, :])[0]
                    if len(join_indices) > 0:
                        join_idx = join_indices[0]

                        s_classes = np.where((fjforkmap == i) & (fjclassmap == r))[0].tolist()
                        new_to_merge = to_merge + s_classes

                        paths = ModelAdapter.find_paths(sn, P, i, join_idx, r, new_to_merge,
                                                        QN, TN, 0, fjclassmap, fjforkmap, nonfjmodel)

                        if len(paths) > 0 and all(x > 0 for x in paths):
                            lambdai = [1.0 / x for x in paths]
                            d0 = 0
                            parallel_branches = len(paths)

                            for pow_val in range(parallel_branches):
                                current_sum = 0
                                for combo in combinations(lambdai, pow_val + 1):
                                    current_sum += 1.0 / sum(combo)
                                d0 += ((-1) ** pow_val) * current_sum

                            ri.extend(ModelAdapter.find_paths(sn, P, join_idx, end_node, r, to_merge,
                                                              QN, TN, current_time + d0,
                                                              fjclassmap, fjforkmap, nonfjmodel))
                        else:
                            ri.extend(ModelAdapter.find_paths(sn, P, join_idx, end_node, r, to_merge,
                                                              QN, TN, current_time,
                                                              fjclassmap, fjforkmap, nonfjmodel))
                    else:
                        ri.extend(ModelAdapter.find_paths(sn, P, i, end_node, r, to_merge,
                                                          QN, TN, current_time + q_len / tput,
                                                          fjclassmap, fjforkmap, nonfjmodel))
                else:
                    ri.extend(ModelAdapter.find_paths(sn, P, i, end_node, r, to_merge,
                                                      QN, TN, current_time + q_len / tput,
                                                      fjclassmap, fjforkmap, nonfjmodel))
            else:
                ri.extend(ModelAdapter.find_paths(sn, P, i, end_node, r, to_merge,
                                                  QN, TN, current_time + q_len / tput,
                                                  fjclassmap, fjforkmap, nonfjmodel))

        return ri

    @staticmethod
    def find_paths_cs(sn: Any, P: np.ndarray, cur_node: int, end_node: int,
                      cur_class: int, to_merge: List[int], QN: np.ndarray,
                      TN: np.ndarray, current_time: float, fjclassmap: np.ndarray,
                      fjforkmap: np.ndarray, nonfjmodel: Any) -> List[float]:
        """
        Find response times along paths with class switching for fork-join analysis.

        Args:
            sn: NetworkStruct
            P: Class-switching routing matrix
            cur_node: Current node
            end_node: Target node
            cur_class: Current class
            to_merge: Classes to merge
            QN: Queue lengths
            TN: Throughputs
            current_time: Accumulated time
            fjclassmap: Fork-join class map
            fjforkmap: Fork-join fork map
            nonfjmodel: Transformed model

        Returns:
            List of response times

        References:
            MATLAB: matlab/src/io/@ModelAdapter/findPathsCS.m
        """
        from itertools import combinations

        if cur_node == end_node:
            ist = int(sn.nodeToStation[cur_node]) if hasattr(sn, 'nodeToStation') else cur_node
            q_len = sum(QN[ist, k] for k in to_merge)
            tput = sum(TN[ist, k] for k in to_merge)
            if tput > 0:
                return [current_time - q_len / tput]
            return [current_time]

        ri = []

        nfjsn = nonfjmodel.getStruct(False) if hasattr(nonfjmodel, 'getStruct') else nonfjmodel
        orignodes = len(nfjsn.rtorig[0][0]) if hasattr(nfjsn, 'rtorig') else sn.nnodes

        row_idx = cur_class * orignodes + cur_node
        transitions = np.where(P[row_idx, :] > 0)[0] if row_idx < P.shape[0] else []

        for transition in transitions:
            cur_merge = to_merge.copy()
            next_class = transition // orignodes
            next_node = transition % orignodes
            cur_merge[0] = next_class

            q_len = 0
            tput = 1

            if hasattr(sn, 'nodeToStation') and not np.isnan(sn.nodeToStation[next_node]):
                ist = int(sn.nodeToStation[next_node])
                q_len = sum(QN[ist, k] for k in cur_merge)
                tput = sum(TN[ist, k] for k in cur_merge) if sum(TN[ist, k] for k in cur_merge) > 0 else 1

            nt = sn.nodetype[next_node]
            is_fork = hasattr(nt, 'name') and nt.name == 'FORK'

            if is_fork and hasattr(sn, 'fj') and sn.fj is not None:
                join_indices = np.where(sn.fj[next_node, :])[0]
                if len(join_indices) > 0:
                    join_idx = join_indices[0]

                    s_classes = np.where((fjforkmap == next_node) & (fjclassmap == next_class))[0].tolist()
                    new_merge = cur_merge + s_classes

                    paths = ModelAdapter.find_paths_cs(sn, P, next_node, join_idx, next_class,
                                                        new_merge, QN, TN, 0,
                                                        fjclassmap, fjforkmap, nonfjmodel)

                    if len(paths) > 0 and all(x > 0 for x in paths):
                        lambdai = [1.0 / x for x in paths]
                        d0 = 0
                        parallel_branches = len(paths)

                        for pow_val in range(parallel_branches):
                            current_sum = 0
                            for combo in combinations(lambdai, pow_val + 1):
                                current_sum += 1.0 / sum(combo)
                            d0 += ((-1) ** pow_val) * current_sum

                        ri.extend(ModelAdapter.find_paths_cs(sn, P, join_idx, end_node, next_class,
                                                              cur_merge, QN, TN, current_time + d0,
                                                              fjclassmap, fjforkmap, nonfjmodel))
                    else:
                        ri.extend(ModelAdapter.find_paths_cs(sn, P, join_idx, end_node, next_class,
                                                              cur_merge, QN, TN, current_time,
                                                              fjclassmap, fjforkmap, nonfjmodel))
                else:
                    ri.extend(ModelAdapter.find_paths_cs(sn, P, next_node, end_node, next_class,
                                                          cur_merge, QN, TN, current_time + q_len / tput,
                                                          fjclassmap, fjforkmap, nonfjmodel))
            else:
                ri.extend(ModelAdapter.find_paths_cs(sn, P, next_node, end_node, next_class,
                                                      cur_merge, QN, TN, current_time + q_len / tput,
                                                      fjclassmap, fjforkmap, nonfjmodel))

        return ri


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
]
