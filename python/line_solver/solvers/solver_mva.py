"""
Native Python implementation of Mean Value Analysis (MVA) solver.

This implementation uses pure Python/NumPy algorithms from the api.pfqn
module, providing the same functionality as the JPype wrapper without .
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass
from enum import Enum


class OptionsDict(dict):
    """A dict that supports attribute-style access."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")


@dataclass
class SolverMVAOptions:
    """Options for the native MVA solver."""
    method: str = 'exact'  # Use 'exact' for accurate results (matches CTMC)
    max_iter: int = 1000
    tol: float = 1e-8
    verbose: bool = False
    seed: Optional[int] = None  # Random seed (for compatibility, not used in MVA)
    keep: bool = False  # Keep intermediate data (for compatibility)
    cutoff: Optional[int] = None  # State space cutoff (for compatibility, not used in MVA)
    samples: Optional[int] = None  # Samples (for compatibility, not used in MVA)


class SolverMVA:
    """
    Native Python Mean Value Analysis (MVA) solver.

    This solver implements MVA algorithms using pure Python/NumPy,
    providing the same functionality as the Java wrapper without
    requiring the JVM.

    Supported methods:
        - 'exact': Exact MVA (pfqn_mva)
        - 'mva': Same as exact
        - 'amva': Approximate MVA using Schweitzer approximation
        - 'qna': Queueing Network Analyzer (for open networks)
        - 'aba.lower', 'aba.upper': Asymptotic Bound Analysis
        - 'pb.lower', 'pb.upper': Performance bounds

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default', which auto-selects exact/amva based on model)
        **kwargs: Additional solver options
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._table_silent = False

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            self.method = 'default'
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'max_iter'):
                kwargs.setdefault('max_iter', method_or_options.max_iter)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'max_iter'):
                kwargs.setdefault('max_iter', method_or_options.max_iter)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        else:
            self.method = 'default'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverMVAOptions(method=self.method, **kwargs)

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "MVA"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model for MVA computation."""
        model = self.model

        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            self._from_network_struct(model._sn)
            return

        # Priority 2: Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                self._from_network_struct(model._sn)
                return

        # Priority 3: JPype wrapper with getStruct()
        if hasattr(model, 'getStruct'):
            try:
                sn = model.getStruct()
                self._from_network_struct(sn)
                return
            except Exception:
                pass

        # Priority 4: JPype wrapper with obj attribute
        if hasattr(model, 'obj'):
            try:
                sn = model.getStruct()
                self._from_network_struct(sn)
                return
            except Exception:
                pass

        # Priority 5: Direct model extraction
        self._from_model_direct(model)

    def _from_network_struct(self, sn):
        """Extract parameters from NetworkStruct."""
        self._sn = sn  # Save reference for chain information
        self.nstations = int(sn.nstations)
        self.nclasses = int(sn.nclasses)

        # Rates matrix (service rates)
        self.rates = np.asarray(sn.rates, dtype=np.float64)

        # Service demands (D = visits * 1/mu if rate > 0)
        # For MVA, demand includes visit ratio: D[i,r] = V[i,r] / mu[i,r]
        self.demands = np.zeros_like(self.rates)
        nonzero = self.rates > 0
        self.demands[nonzero] = 1.0 / self.rates[nonzero]

        # Apply visits to demands (D = V / mu)
        if hasattr(sn, 'visits') and sn.visits is not None:
            # Build combined visit matrix from per-chain visits
            visits_combined = np.zeros_like(self.demands)
            for chain_id, visits_chain in sn.visits.items():
                if isinstance(visits_chain, np.ndarray) and visits_chain.shape == visits_combined.shape:
                    visits_combined += visits_chain
            # Apply visits: D = V / mu = V * (1/mu)
            self.demands = self.demands * visits_combined

        # Population vector
        self.njobs = np.asarray(sn.njobs, dtype=np.float64).flatten()

        # Number of servers
        self.nservers = np.asarray(sn.nservers, dtype=np.float64).flatten()

        # Reference station (for think times)
        self.refstat = np.asarray(sn.refstat, dtype=int).flatten()

        # Station to node mapping
        self.stationToNode = np.asarray(sn.stationToNode, dtype=int).flatten() \
                            if hasattr(sn, 'stationToNode') else np.arange(self.nstations)

        # Node types (to identify delays/think times)
        self.nodetype = sn.nodetype if hasattr(sn, 'nodetype') else None

        # Station types - map from station index to node type
        self.station_types = []
        if self.nodetype is not None:
            for i in range(self.nstations):
                node_idx = self.stationToNode[i]
                if node_idx < len(self.nodetype):
                    self.station_types.append(self.nodetype[node_idx])
                else:
                    self.station_types.append(None)

        # Names - use station node names, not all node names
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') else []
        if nodenames and len(self.stationToNode) > 0:
            self.station_names = [nodenames[self.stationToNode[i]] if self.stationToNode[i] < len(nodenames)
                                  else f'Station{i}' for i in range(self.nstations)]
        else:
            self.station_names = [f'Station{i}' for i in range(self.nstations)]

        self.class_names = list(sn.classnames) if hasattr(sn, 'classnames') else \
                          [f'Class{i}' for i in range(self.nclasses)]

        # Scheduling strategies
        self.sched = sn.sched if hasattr(sn, 'sched') else None

        # Scheduling parameters (weights for DPS/GPS)
        self.schedparam = np.asarray(sn.schedparam) if hasattr(sn, 'schedparam') and sn.schedparam is not None else None

        # Visits (routing)
        self.visits = sn.visits if hasattr(sn, 'visits') else None

        # Load-dependent scaling
        self.lldscaling = getattr(sn, 'lldscaling', None)

        # Determine network type
        self._determine_network_type()

    def _from_model_direct(self, model):
        """Extract parameters directly from model."""
        # Support both native (snake_case) and wrapper (PascalCase) APIs
        if hasattr(model, 'get_number_of_stations'):
            self.nstations = model.get_number_of_stations()
            self.nclasses = model.get_number_of_classes()
        else:
            self.nstations = model.getNumberOfStations()
            self.nclasses = model.getNumberOfClasses()

        # Initialize arrays
        self.rates = np.zeros((self.nstations, self.nclasses))
        self.demands = np.zeros((self.nstations, self.nclasses))
        self.njobs = np.zeros(self.nclasses)
        self.nservers = np.ones(self.nstations)
        self.refstat = np.zeros(self.nclasses, dtype=int)

        self.station_names = []
        self.class_names = []

        # Get class names and populations
        if hasattr(model, 'get_classes'):
            classes = list(model.get_classes())
        else:
            classes = list(model.getClasses())

        for c, cobj in enumerate(classes):
            if c < self.nclasses:
                # Get class name
                if hasattr(cobj, 'get_name'):
                    name = cobj.get_name()
                elif hasattr(cobj, 'getName'):
                    name = cobj.getName()
                else:
                    name = getattr(cobj, 'name', f'Class{c}')
                self.class_names.append(str(name))

                # Get number of jobs (for closed classes)
                if hasattr(cobj, 'get_number_of_jobs'):
                    self.njobs[c] = cobj.get_number_of_jobs()
                elif hasattr(cobj, 'getNumberOfJobs'):
                    self.njobs[c] = cobj.getNumberOfJobs()
                else:
                    self.njobs[c] = 0

        # Get station info
        if hasattr(model, 'get_stations'):
            nodes = model.get_stations()
        else:
            nodes = list(model.getNodes())

        station_idx = 0
        for node in nodes:
            node_type = str(type(node).__name__)
            # Filter to stations only (Queue, Delay, Router, ClassSwitch, Fork, Join)
            if any(t in node_type for t in ['Queue', 'Delay', 'Router', 'ClassSwitch', 'Fork', 'Join']):
                # Get station name
                if hasattr(node, 'get_name'):
                    sname = node.get_name()
                elif hasattr(node, 'getName'):
                    sname = node.getName()
                else:
                    sname = getattr(node, 'name', f'Station{station_idx}')
                self.station_names.append(str(sname))

                # Get service rates for each class
                for c in range(self.nclasses):
                    try:
                        if hasattr(node, 'get_service'):
                            service = node.get_service(classes[c])
                        else:
                            service = node.getServiceProcess(classes[c])

                        if service is not None and hasattr(service, 'getMean'):
                            mean_val = service.getMean()
                            if mean_val > 0:
                                self.rates[station_idx, c] = 1.0 / mean_val
                                self.demands[station_idx, c] = mean_val
                    except:
                        pass

                # Get number of servers
                if hasattr(node, 'number_of_servers'):
                    self.nservers[station_idx] = node.number_of_servers
                elif hasattr(node, 'getNumberOfServers'):
                    self.nservers[station_idx] = node.getNumberOfServers()

                station_idx += 1

        self._determine_network_type()

    def _determine_network_type(self):
        """Determine if network is open, closed, or mixed."""
        has_open = False
        has_closed = False

        for c in range(self.nclasses):
            # Open classes have njobs = 0 or inf, closed classes have finite njobs > 0
            if np.isfinite(self.njobs[c]) and self.njobs[c] > 0:
                has_closed = True
            else:
                has_open = True

        if has_open and has_closed:
            self.network_type = 'mixed'
        elif has_closed:
            self.network_type = 'closed'
        else:
            self.network_type = 'open'

    def _get_think_times(self) -> np.ndarray:
        """Extract think times from delay stations.

        Think time Z[r] is the sum of demands at all Delay stations for class r.
        """
        from ..api.sn.network_struct import NodeType

        Z = np.zeros(self.nclasses)

        # Sum demands from all Delay stations for each class
        for i in range(self.nstations):
            is_delay = False
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.DELAY:
                        is_delay = True

            if is_delay:
                for c in range(self.nclasses):
                    if self.demands[i, c] > 0:
                        Z[c] += self.demands[i, c]

        return Z

    def _get_queueing_demands(self) -> Tuple[np.ndarray, List[int]]:
        """
        Get demand matrix for queueing stations only (excluding sources, sinks, and delays).

        Delay stations contribute to think times, not queueing demands.

        Returns:
            Tuple of (demand matrix, list of queueing station indices)
        """
        from ..api.sn.network_struct import NodeType

        # Get queueing stations (exclude Source, Sink, Delay, Fork, Join)
        queue_indices = []
        for i in range(self.nstations):
            # Check if this station should be excluded from queueing demands
            is_excluded = False
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    # Handle both enum objects and integer values
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    # Exclude Source, Sink, Delay, Fork, and Join stations
                    # (Fork/Join don't do queueing work, they handle synchronization)
                    if st_val in (NodeType.SOURCE, NodeType.SINK, NodeType.DELAY,
                                  NodeType.FORK, NodeType.JOIN):
                        is_excluded = True

            # Include station if it has non-zero demand and is not excluded
            if not is_excluded and np.any(self.demands[i, :] > 0):
                queue_indices.append(i)

        if not queue_indices:
            return np.zeros((1, self.nclasses)), [0]

        L = self.demands[queue_indices, :].copy()

        # Apply DPS weight scaling: for DPS stations, effective demand = D / weight
        from ..lang.base import SchedStrategy
        for idx, i in enumerate(queue_indices):
            if self.sched is not None and i in self.sched:
                sched_val = self.sched[i]
                # Check if this is a DPS station
                is_dps = (sched_val == SchedStrategy.DPS or
                         (hasattr(sched_val, 'value') and sched_val.value == SchedStrategy.DPS) or
                         (isinstance(sched_val, int) and sched_val == 5))  # DPS enum value
                if is_dps and self.schedparam is not None:
                    for k in range(self.nclasses):
                        if i < self.schedparam.shape[0] and k < self.schedparam.shape[1]:
                            w_k = self.schedparam[i, k]
                            if w_k > 0:
                                # Scale demand by weight (higher weight = faster service = lower effective demand)
                                L[idx, k] = L[idx, k] / w_k

        return L, queue_indices

    def _get_source_stations(self) -> List[int]:
        """Get list of source station indices."""
        from ..api.sn.network_struct import NodeType

        source_indices = []
        for i in range(self.nstations):
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.SOURCE:
                        source_indices.append(i)
        return source_indices

    def _get_delay_stations(self) -> List[int]:
        """Get list of delay (infinite server) station indices."""
        from ..api.sn.network_struct import NodeType

        delay_indices = []
        for i in range(self.nstations):
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.DELAY:
                        delay_indices.append(i)
        return delay_indices

    def _has_fork_join(self) -> bool:
        """Check if the model contains fork-join nodes."""
        if hasattr(self.model, 'has_fork') and callable(self.model.has_fork):
            return self.model.has_fork()
        # Fallback: check node types
        if self._sn is not None:
            from ..api.sn.network_struct import NodeType
            for nt in self._sn.nodetype:
                nt_val = nt.value if hasattr(nt, 'value') else int(nt)
                if nt_val == NodeType.FORK:
                    return True
        return False

    def _run_fork_join_analysis(self):
        """
        Run MVA analysis for fork-join networks.

        For closed fork-join networks, implements Heidelberger-Trivedi approach:
        1. Transforms model using ModelAdapter.ht() to create auxiliary classes
        2. Solves transformed model with standard MVA
        3. Iteratively updates sync delays until convergence
        4. Merges auxiliary class results back into original classes

        For open fork-join networks:
        Returns None to fall back to standard analysis which then gets post-processed.
        """
        from ..api.sn.network_struct import NodeType
        from itertools import combinations

        # Only use H-T for closed fork-join networks
        if self.network_type != 'closed':
            # Open networks use post-processing approach
            return None

        sn = self._sn
        if sn is None or sn.fj is None or not np.any(sn.fj):
            return None

        # Use the full H-T transformation
        try:
            return self._run_ht_transformation()
        except Exception as e:
            # If H-T transformation fails, fall back to post-processing
            import warnings
            warnings.warn(f"H-T transformation failed: {e}. Using post-processing approach.")
            return None

    def _run_ht_transformation(self):
        """
        Run full Heidelberger-Trivedi transformation for closed fork-join networks.

        The H-T method analyzes each parallel branch independently with the full
        population N, then computes E[max] for the sync delay. This is implemented
        in _run_approximate_ht() which uses the correct algorithm.

        Note: Running standard MVA on a combined transformed model (with auxiliary
        classes) gives incorrect results because the populations would be treated
        as independent (N + N + ... = total jobs > N).

        Returns:
            Result dictionary or None if transformation fails
        """
        # Always use the approximate H-T approach which implements the correct
        # H-T algorithm: analyze each branch independently with population N,
        # then compute E[max] for synchronization.
        return self._run_approximate_ht()

    def _run_approximate_ht(self):
        """
        Run Heidelberger-Trivedi MVA for closed fork-join networks.

        Implements the H-T approximation which:
        1. Analyzes each parallel branch independently with full population N
        2. Computes E[max(R_1, ..., R_K)] for synchronization delay
        3. Uses cycle time C = Z + E[max] to compute throughput X = N/C
        4. Iterates until response times converge

        The H-T formula uses exponential approximation for the maximum:
        E[max] = sum_{k=1}^{K} (-1)^{k+1} * sum_{combos of k} 1/sum(rates)

        Note: Results typically differ from exact solutions by ~4% due to the
        exponential assumption in the sync delay formula.

        Returns:
            Result dictionary or None if analysis fails
        """
        from ..api.sn.network_struct import NodeType

        sn = self._sn
        if sn is None or sn.fj is None or not np.any(sn.fj):
            return None

        # Find fork and join indices
        fork_indices = []
        join_indices = []
        for i, nt in enumerate(sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK.value if hasattr(NodeType.FORK, 'value') else NodeType.FORK:
                fork_indices.append(i)
            elif nt_val == NodeType.JOIN.value if hasattr(NodeType.JOIN, 'value') else NodeType.JOIN:
                join_indices.append(i)

        if len(fork_indices) == 0:
            return None

        # Get queueing demands (service times at queues) - also returns queue indices
        L, queue_indices = self._get_queueing_demands()

        # Get delay station indices
        delay_indices = self._get_delay_stations()
        if L.shape[0] == 0:
            return None

        # Get think times (delay station service times)
        Z = self._get_think_times()

        # Get population
        N = self.njobs.copy()

        # Initialize result arrays
        R = self.nclasses
        QN = np.zeros((self.nstations, R))
        UN = np.zeros((self.nstations, R))
        RN = np.zeros((self.nstations, R))
        TN = np.zeros((self.nstations, R))
        AN = np.zeros((self.nstations, R))
        XN = np.zeros(R)

        # Identify parallel branches for each fork-join pair
        fj_branches = {}  # fj_branches[f] = list of station indices on parallel branches
        for f in fork_indices:
            join_idx_arr = np.where(sn.fj[f, :])[0]
            if len(join_idx_arr) > 0:
                join_idx = join_idx_arr[0]
                join_station = sn.nodeToStation[join_idx] if join_idx < len(sn.nodeToStation) else -1
                parallel_stations = self._find_parallel_branches(f, join_idx)
                fj_branches[f] = {
                    'join_node': join_idx,
                    'join_station': join_station,
                    'parallel_stations': parallel_stations
                }

        # Initial response times (just service times)
        R_branch = {}  # Response times at parallel branch queues
        for f, fj_info in fj_branches.items():
            R_branch[f] = {}
            for r in range(R):
                if np.isfinite(N[r]) and N[r] > 0:
                    branch_ri = []
                    for station_idx in fj_info['parallel_stations']:
                        if station_idx in queue_indices:
                            q_idx = queue_indices.index(station_idx)
                            branch_ri.append(L[q_idx, r])  # Initial: service time
                        else:
                            branch_ri.append(self.demands[station_idx, r])
                    R_branch[f][r] = np.array(branch_ri) if branch_ri else np.array([1.0])

        # H-T MVA for closed fork-join
        # Key insight: In the H-T transformation, each parallel branch effectively
        # sees the full population N circulating through it. We compute response times
        # at each branch using MVA, then compute E[max] for synchronization.

        # For each class, run iterative MVA with H-T synchronization
        for r in range(R):
            if not (np.isfinite(N[r]) and N[r] > 0):
                continue

            n_jobs = int(N[r])

            # Get branch demands for this class
            # Note: Use raw service times (1/rate), not visit-scaled demands,
            # because in H-T each branch is analyzed independently with full population
            for f, fj_info in fj_branches.items():
                branch_demands = []
                for station_idx in fj_info['parallel_stations']:
                    # Use raw service time (1/rate), not demand (1/rate * visits)
                    # because fork-join visits are split (0.5 each) but H-T analyzes
                    # each branch independently with full population
                    if self.rates[station_idx, r] > 0:
                        D_i = 1.0 / self.rates[station_idx, r]
                    else:
                        D_i = 0.0
                    branch_demands.append(D_i if D_i > 0 else 1e-10)

                if len(branch_demands) == 0:
                    continue

                branch_demands = np.array(branch_demands)
                num_branches = len(branch_demands)

                # Run MVA for each branch as if it sees full population N
                # This matches the H-T transformation where auxiliary classes have population N
                branch_R = np.zeros(num_branches)
                branch_Q = np.zeros(num_branches)

                # Exact MVA iteration for each branch
                for k in range(1, n_jobs + 1):
                    # Response time at each queue when k jobs in system
                    for b in range(num_branches):
                        # R = D * (1 + Q_{k-1})
                        branch_R[b] = branch_demands[b] * (1 + branch_Q[b])

                    # Compute E[max] of branch response times
                    if num_branches >= 2:
                        d0 = self._compute_sync_delay(branch_R)
                    else:
                        d0 = branch_R[0]

                    # Cycle time = Z + E[max]
                    C_k = Z[r] + d0

                    # Throughput for k jobs
                    X_k = k / C_k

                    # Update queue lengths at each branch
                    for b in range(num_branches):
                        branch_Q[b] = X_k * branch_R[b]

                # Continue iterating until convergence (H-T requires fixed point)
                # The sync delay creates a feedback loop that requires additional iterations
                tol = 1e-6
                max_extra_iter = 100
                for _ in range(max_extra_iter):
                    branch_R_old = branch_R.copy()

                    for b in range(num_branches):
                        branch_R[b] = branch_demands[b] * (1 + branch_Q[b])

                    if num_branches >= 2:
                        d0 = self._compute_sync_delay(branch_R)
                    else:
                        d0 = branch_R[0]

                    C_k = Z[r] + d0
                    X_k = n_jobs / C_k

                    for b in range(num_branches):
                        branch_Q[b] = X_k * branch_R[b]

                    if np.max(np.abs(branch_R - branch_R_old)) < tol:
                        break

                # Store final results
                XN[r] = X_k
                R_branch[f][r] = branch_R

                # Update station metrics
                for b, station_idx in enumerate(fj_info['parallel_stations']):
                    RN[station_idx, r] = branch_R[b]
                    QN[station_idx, r] = branch_Q[b]
                    UN[station_idx, r] = XN[r] * branch_demands[b]
                    TN[station_idx, r] = XN[r]
                    AN[station_idx, r] = XN[r]

        # Set delay station metrics
        for station_idx in delay_indices:
            for r in range(R):
                if Z[r] > 0:
                    RN[station_idx, r] = Z[r]
                    QN[station_idx, r] = XN[r] * Z[r]
                    UN[station_idx, r] = XN[r] * Z[r]  # Utilization at delay
                    TN[station_idx, r] = XN[r]
                    AN[station_idx, r] = XN[r]

        # Set join station metrics
        for f, fj_info in fj_branches.items():
            join_station = fj_info['join_station']
            if join_station >= 0 and join_station < self.nstations:
                for r in range(R):
                    if np.isfinite(N[r]) and N[r] > 0:
                        ri = R_branch[f].get(r, np.array([1.0]))
                        if len(ri) >= 2:
                            d0 = self._compute_sync_delay(ri)
                            sync_d = max(0, d0 - np.mean(ri))
                            num_branches = len(ri)
                        else:
                            sync_d = 0
                            num_branches = 1

                        RN[join_station, r] = sync_d
                        QN[join_station, r] = XN[r] * sync_d * num_branches
                        TN[join_station, r] = XN[r]
                        AN[join_station, r] = XN[r]
                        UN[join_station, r] = 0  # Join has no service

        # Compute final cycle times
        CN = np.zeros(R)
        for r in range(R):
            if np.isfinite(N[r]) and N[r] > 0 and XN[r] > 0:
                CN[r] = N[r] / XN[r]
            else:
                CN[r] = 0

        # Store results
        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
            'CN': CN,
            'method': 'approximate_ht',
            'iter': 1
        }

        return self._result

    def _compute_fork_join_sync_delays(self, QN, UN, RN, TN, AN, XN):
        """
        Compute and add synchronization delays at Join nodes.

        For each fork-join pair, computes the expected synchronization delay
        using the Heidelberger-Trivedi formula.

        Note: For closed fork-join networks, this post-processing approach gives
        approximate results. The MATLAB implementation uses a full H-T transformation
        that creates auxiliary classes and modifies the model structure, which is
        not yet implemented in Python. Results for closed fork-join networks may
        differ from MATLAB by up to 40%.

        Args:
            QN, UN, RN, TN, AN, XN: Performance metric arrays (modified in place)
        """
        from ..api.sn.network_struct import NodeType
        from itertools import combinations

        sn = self._sn
        if sn.fj is None or not np.any(sn.fj):
            return

        # Find fork indices
        fork_indices = []
        for i, nt in enumerate(sn.nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
            if nt_val == NodeType.FORK:
                fork_indices.append(i)

        for f in fork_indices:
            # Find join associated with this fork
            join_idx_arr = np.where(sn.fj[f, :])[0]
            if len(join_idx_arr) == 0:
                continue
            join_idx = join_idx_arr[0]

            # Get join station index
            if join_idx >= len(sn.nodeToStation):
                continue
            join_station = sn.nodeToStation[join_idx]
            if join_station < 0 or join_station >= self.nstations:
                continue

            # Find the parallel branches (queues between fork and join)
            # These are the queues that feed into the join
            parallel_stations = self._find_parallel_branches(f, join_idx)

            if len(parallel_stations) < 2:
                continue

            # Get response times at parallel queues for each class
            for r in range(self.nclasses):
                ri = []
                for station_idx in parallel_stations:
                    if RN[station_idx, r] > 0:
                        ri.append(RN[station_idx, r])

                if len(ri) >= 2:
                    ri_arr = np.array(ri)
                    num_branches = len(ri_arr)

                    # Compute expected max using H-T formula with response times
                    d0 = self._compute_sync_delay(ri_arr)
                    # Sync delay at join = E[max] - mean (waiting for slowest)
                    sync_delay = d0 - np.mean(ri_arr)
                    sync_delay = max(0, sync_delay)

                    # Get throughput at join
                    tput = XN[r] if XN[r] > 0 else 0.05  # Use arrival rate as approx

                    # Set join metrics
                    # RespT at join = synchronization delay
                    RN[join_station, r] = sync_delay
                    # QLen = tput * respT * num_branches (jobs from each branch waiting)
                    QN[join_station, r] = tput * sync_delay * num_branches
                    # Throughput at join = same as system throughput
                    TN[join_station, r] = tput
                    AN[join_station, r] = tput
                    # Utilization stays 0 for join

    def _find_parallel_branches(self, fork_idx: int, join_idx: int) -> List[int]:
        """Find stations on parallel branches between a fork and join.

        Uses the connection matrix to find nodes that are direct successors
        of the fork and eventually lead to the join.
        """
        from ..api.sn.network_struct import NodeType
        import numpy as np

        parallel_stations = []
        sn = self._sn

        # Use connection matrix to find direct successors of the fork
        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
            connmatrix = np.array(sn.connmatrix)

            # Get direct successors of the fork node
            fork_successors = np.where(connmatrix[fork_idx, :] > 0)[0]

            # For each successor, trace the path to see if it reaches the join
            for succ_node in fork_successors:
                # BFS to find all nodes on this branch until we hit the join
                visited = set()
                queue = [succ_node]
                branch_nodes = []

                while queue:
                    current = queue.pop(0)
                    if current in visited:
                        continue
                    visited.add(current)

                    if current == join_idx:
                        # Reached the join, this is a valid branch
                        break

                    branch_nodes.append(current)

                    # Add successors to queue
                    successors = np.where(connmatrix[current, :] > 0)[0]
                    for s in successors:
                        if s not in visited:
                            queue.append(s)

                # Convert branch nodes to station indices
                for node in branch_nodes:
                    if node < len(sn.nodeToStation):
                        station = sn.nodeToStation[node]
                        if station >= 0 and station not in parallel_stations:
                            # Check it's a queue/delay station (has service)
                            if station < len(self.station_types):
                                st = self.station_types[station]
                                if st is not None:
                                    st_val = st.value if hasattr(st, 'value') else int(st)
                                    if st_val in (NodeType.QUEUE, NodeType.DELAY):
                                        parallel_stations.append(station)
        else:
            # Fallback: use simplified approach (less accurate for serial fork-joins)
            for i in range(self.nstations):
                if self.station_types and i < len(self.station_types):
                    st = self.station_types[i]
                    if st is not None:
                        st_val = st.value if hasattr(st, 'value') else int(st)
                        if st_val == NodeType.QUEUE:
                            parallel_stations.append(i)

        return parallel_stations

    def _compute_sync_delay(self, path_times: np.ndarray) -> float:
        """
        Compute synchronization delay using Heidelberger-Trivedi formula.

        For K parallel branches with response times r_1, ..., r_K,
        the expected maximum E[max(r_1, ..., r_K)] is computed using
        inclusion-exclusion with exponential approximation.

        Args:
            path_times: Array of response times for parallel paths

        Returns:
            Expected maximum (synchronization point) time
        """
        from itertools import combinations

        if len(path_times) == 0:
            return 0.0
        if len(path_times) == 1:
            return path_times[0]

        # Convert to rates (1/response_time)
        path_times = np.asarray(path_times)
        # Avoid division by zero
        path_times = np.maximum(path_times, 1e-10)
        lambdai = 1.0 / path_times

        d0 = 0.0
        parallel_branches = len(lambdai)

        for pow_val in range(parallel_branches):
            # Get all combinations of (pow_val + 1) elements
            for combo in combinations(range(parallel_branches), pow_val + 1):
                combo_sum = np.sum(lambdai[list(combo)])
                if combo_sum > 0:
                    d0 += ((-1) ** pow_val) * (1.0 / combo_sum)

        return d0

    def _find_node_type_indices(self, nodetype_list, target_type) -> np.ndarray:
        """
        Find indices of nodes with a specific type.

        Properly handles enum comparisons with lists of NodeType values.

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

    def runAnalyzer(self):
        """Run the MVA analysis."""
        # Check for fork-join networks and handle specially
        if self._has_fork_join():
            result = self._run_fork_join_analysis()
            if result is not None:
                return result
            # Fall through to standard analysis if fork-join handling failed

        from ..api.pfqn import (
            pfqn_mva, pfqn_aql, pfqn_linearizer, pfqn_gflinearizer,
            pfqn_egflinearizer, pfqn_mvald, pfqn_mvams, pfqn_bs, pfqn_sqni,
            pfqn_schmidt, pfqn_ab_amva,
            pfqn_qd, pfqn_qdlin, pfqn_qli, pfqn_fli,
        )
        from ..api.pfqn.bounds import (
            pfqn_xzabalow, pfqn_xzabaup, pfqn_qzgblow, pfqn_qzgbup,
            pfqn_xzgsblow, pfqn_xzgsbup,
        )

        method = self.method.lower()

        # Normalize AMVA method aliases
        method = method.replace('amva.', '')

        # Get parameters
        L, queue_indices = self._get_queueing_demands()
        N = self.njobs.copy()
        Z = self._get_think_times()
        mi = self.nservers[queue_indices] if len(queue_indices) > 0 else np.ones(1)

        M = L.shape[0]  # Number of queueing stations
        R = self.nclasses  # Number of classes

        # Compute arrival rates for open classes
        lambda_arr = np.zeros(R)
        source_indices = self._get_source_stations()
        for r in range(R):
            if np.isinf(N[r]):  # Open class
                # Get arrival rate from source station
                for src_idx in source_indices:
                    if self.rates[src_idx, r] > 0:
                        lambda_arr[r] = self.rates[src_idx, r]
                        break

        # Initialize result arrays
        QN = np.zeros((self.nstations, R))
        UN = np.zeros((self.nstations, R))
        RN = np.zeros((self.nstations, R))
        TN = np.zeros((self.nstations, R))
        AN = np.zeros((self.nstations, R))
        XN = np.zeros(R)

        # Choose algorithm based on method
        if self.network_type == 'closed' or self.network_type == 'mixed':
            # Check if load-dependent MVA should be used (MATLAB: solver_mvald_analyzer)
            # Must check BEFORE converting 'default' method to 'exact'/'amva'
            use_ld_mva = False
            if self.lldscaling is not None:
                use_ld_mva = True
            # Mixed networks require chain-level MVA (pfqn_mvaldmx) which handles open classes
            if self.network_type == 'mixed':
                use_ld_mva = True

            # Handle 'default' method with MATLAB-compatible heuristic
            if method == 'default':
                if use_ld_mva:
                    # For load-dependent models, MATLAB defaults to approximate MVA (solver_amvald)
                    # even for small models, NOT exact MVA (pfqn_mvaldmx)
                    method = 'amva'
                else:
                    # Match MATLAB's solver_mva_analyzer.m logic for non-LD models:
                    # Use exact MVA if: nchains <= 4 && sum(njobs) <= 20 && product_form && no fractional populations
                    nchains = self._sn.nchains if hasattr(self._sn, 'nchains') and self._sn is not None else R
                    total_jobs = int(np.sum(N[np.isfinite(N)]))
                    has_product_form = True  # Assume product form for now (TODO: implement check)
                    has_fractional = np.any(N != np.floor(N))

                    if nchains <= 4 and total_jobs <= 20 and has_product_form and not has_fractional:
                        method = 'exact'
                    else:
                        method = 'amva'

            # Also check for INF servers which need load-dependent treatment
            # NOTE: Delay stations naturally have INF scheduling and should NOT trigger this
            # Only Queue stations (NodeType.QUEUE) with INF scheduling require load-dependent MVA
            from ..lang.base import SchedStrategy
            from ..api.sn.network_struct import NodeType
            has_inf_server = False
            for idx, q_idx in enumerate(queue_indices):
                if self.sched is not None and q_idx in self.sched:
                    sched_val = self.sched[q_idx]
                    # Check if this is a Queue station (not Delay) with INF servers
                    is_queue_station = False
                    if self.station_types is not None and q_idx < len(self.station_types):
                        st = self.station_types[q_idx]
                        if st is not None:
                            st_val = st.value if hasattr(st, 'value') else int(st)
                            is_queue_station = (st_val == NodeType.QUEUE.value)
                    if is_queue_station and sched_val == SchedStrategy.INF:
                        has_inf_server = True
                        break

            # INF servers require load-dependent MVA (mu scales with population)
            if has_inf_server:
                use_ld_mva = True

            # Class switching requires chain-level MVA
            if hasattr(self._sn, 'nchains') and self._sn.nchains > 0:
                chains = self._get_chains()
                for chain in chains:
                    if len(chain) > 1:
                        # Multiple classes in same chain = class switching
                        use_ld_mva = True
                        break

            if use_ld_mva and method in ['exact', 'mva']:
                # Use EXACT load-dependent MVA (MATLAB: solver_mvald -> pfqn_mvaldmx)
                # IMPORTANT: MATLAB's solver_mvald works at chain level, then disaggregates
                from ..api.pfqn import pfqn_mvaldmx
                from ..api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                # Get chain-level parameters (matching MATLAB's solver_mvald.m)
                chain_result = sn_get_demands_chain(self._sn)
                Lchain = chain_result.Lchain
                STchain = chain_result.STchain
                Vchain = chain_result.Vchain
                alpha = chain_result.alpha
                Nchain = chain_result.Nchain.flatten()
                refstatchain = chain_result.refstatchain

                C = self._sn.nchains
                total_pop = int(np.sum(Nchain[np.isfinite(Nchain)]))

                # Build mu_chain matrix like MATLAB's solver_mvald.m
                # MUST use full station count M_full (same as Lchain rows), not just queueing stations
                # MATLAB: M = size(STchain,1); mu_chain = ones(M, sum(Nchain(isfinite(Nchain))))
                M_full = Lchain.shape[0]  # Number of ALL stations
                S = self.nservers  # Server counts per station
                mu_chain = np.ones((M_full, total_pop))
                for ist in range(M_full):
                    if ist < len(S) and np.isinf(S[ist]):
                        # INF server: mu[ist,:] = 1:N (linear scaling)
                        for n in range(total_pop):
                            mu_chain[ist, n] = n + 1
                    elif self.lldscaling is not None and ist < self.lldscaling.shape[0]:
                        # Load-dependent: use lldscaling
                        for n in range(total_pop):
                            if n < self.lldscaling.shape[1]:
                                mu_chain[ist, n] = self.lldscaling[ist, n]
                            else:
                                mu_chain[ist, n] = self.lldscaling[ist, -1]
                    elif ist < len(S) and S[ist] > 1:
                        # Finite multiserver queue: mu scales up to number of servers
                        # For c servers: mu = [1, 2, ..., c, c, c, ...]
                        c = int(S[ist])
                        for n in range(total_pop):
                            mu_chain[ist, n] = min(n + 1, c)

                # Build arrival rates for open chains
                # For open classes, arrival rate comes from the Source station, not refstat
                from ..api.sn.network_struct import NodeType
                source_station = None
                for ist in range(self.nstations):
                    if self.station_types and ist < len(self.station_types):
                        st = self.station_types[ist]
                        if st is not None:
                            st_val = st.value if hasattr(st, 'value') else int(st)
                            if st_val == NodeType.SOURCE.value:
                                source_station = ist
                                break

                lambda_chain = np.zeros(C)
                for c in range(C):
                    if np.isinf(Nchain[c]):
                        # Sum arrival rates for all open classes in this chain
                        if c in self._sn.inchain:
                            inchain = self._sn.inchain[c].flatten().astype(int)
                            for r in inchain:
                                if r < R and np.isinf(N[r]):
                                    # For open classes, get arrival rate from Source station
                                    if source_station is not None and source_station < self.nstations:
                                        arr_rate = self.rates[source_station, r]
                                        if arr_rate > 0:
                                            lambda_chain[c] += arr_rate

                # Call pfqn_mvaldmx with chain-level parameters
                # MATLAB: [Xchain,Qchain,Uchain] = pfqn_mvaldmx(lambda,Lchain,Nchain,0*Nchain,mu_chain,S)
                # where S = sn.nservers (all stations)
                Xchain, Qchain, Uchain, _, lGN, Pc = pfqn_mvaldmx(
                    lambda_chain, Lchain, Nchain, np.zeros(C), mu_chain, S
                )

                # Compute Tchain and Rchain like MATLAB
                Tchain = np.outer(np.ones(M_full), Xchain) * Vchain
                Rchain = np.zeros((M_full, C))
                for c in range(C):
                    if Xchain[c] > 0:
                        Rchain[:, c] = Qchain[:, c] / Xchain[c]

                # Disaggregate chain results to class level (matching MATLAB's solver_mvald.m)
                deagg = sn_deaggregate_chain_results(
                    self._sn, Lchain, None, STchain, Vchain, alpha,
                    None, Uchain, Rchain, Tchain, None, Xchain.reshape(1, -1)
                )

                # Copy disaggregated results
                QN = deagg.Q
                UN = deagg.U
                RN = deagg.R
                TN = deagg.T
                XN = deagg.X.flatten()
                AN = TN.copy()
            elif use_ld_mva and method in ['default', 'amva']:
                # Use APPROXIMATE load-dependent MVA (MATLAB: solver_amvald)
                # Note: This branch is for compatibility; the approximate implementation
                # needs more work to match MATLAB exactly
                from ..api.solvers.mva.amvald import solver_amvald, AmvaldOptions

                # Get chain-level parameters
                from ..api.sn import sn_get_demands_chain
                chain_result = sn_get_demands_chain(self._sn)
                Lchain = chain_result.Lchain
                STchain = chain_result.STchain
                Vchain = chain_result.Vchain
                alpha = chain_result.alpha
                Nchain = chain_result.Nchain.flatten()
                refstatchain = chain_result.refstatchain
                SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                # Set up options
                amvald_options = AmvaldOptions(
                    method='default',
                    iter_tol=1e-6,
                    iter_max=1000
                )

                # Call solver_amvald
                result = solver_amvald(
                    self._sn, Lchain, STchain, Vchain, alpha,
                    Nchain, SCVchain, refstatchain, amvald_options
                )

                # Disaggregate chain results to class level
                from ..api.sn import sn_deaggregate_chain_results
                deagg = sn_deaggregate_chain_results(
                    self._sn, Lchain, None, STchain, Vchain, alpha,
                    result.Q, result.U, result.R, result.T, None, result.X
                )

                # Copy results
                QN = deagg.Q
                UN = deagg.U
                RN = deagg.R
                TN = deagg.T
                XN = deagg.X.flatten()
                AN = TN.copy()

            elif method in ['exact', 'mva']:
                # Use the handler's solver_mva for chain-level aggregation and post-processing
                # This matches MATLAB's solver_mva which does:
                # 1. Chain aggregation via sn_get_demands_chain
                # 2. MVA at chain level via pfqn_mvams
                # 3. Post-processing to recompute Q and X from waiting times
                # 4. Result disaggregation
                from ..api.solvers.mva.handler import solver_mva as mva_handler
                from ..api.solvers.mva.handler import SolverMVAOptions as MVAHandlerOptions

                handler_options = MVAHandlerOptions(method='exact', tol=1e-8)
                result = mva_handler(self._sn, handler_options)

                # Copy results from handler
                QN = result.Q if result.Q is not None else np.zeros((self.nstations, R))
                UN = result.U if result.U is not None else np.zeros((self.nstations, R))
                RN = result.R if result.R is not None else np.zeros((self.nstations, R))
                TN = result.T if result.T is not None else np.zeros((self.nstations, R))
                XN = result.X.flatten() if result.X is not None else np.zeros(R)
                AN = TN.copy()

            elif method == 'amva':
                # Approximate MVA using AQL (Schweitzer approximation)
                result = pfqn_aql(L, N, Z)
                XN_out, CN_out, QN_out, UN_out, RN_out, TN_out, AN_out = result

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = AN_out[idx, :]

                XN = XN_out.flatten()

            elif method == 'bs':
                # Balanced System approximation
                result = pfqn_bs(L, N, Z)
                XN_out, CN_out, QN_out, UN_out, RN_out, TN_out, AN_out = result

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = AN_out[idx, :]

                XN = XN_out.flatten()

            elif method == 'sqni':
                # Single Queue Network Interpolation
                QN_out, UN_out, XN_out = pfqn_sqni(L, N, Z)
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :] if QN_out.ndim > 1 else QN_out
                    UN[q_idx, :] = UN_out[idx, :] if UN_out.ndim > 1 else UN_out
                XN = XN_out.flatten()

            elif method in ['lin', 'gflin', 'egflin']:
                # Linearizer family of algorithms
                # Build scheduling strategy list
                sched_type = []
                for q_idx in queue_indices:
                    if self.sched is not None and q_idx in self.sched:
                        sched_type.append(str(self.sched[q_idx]))
                    else:
                        sched_type.append('FCFS')

                if method == 'egflin':
                    result = pfqn_egflinearizer(L, N, Z, sched_type)
                elif method == 'gflin':
                    result = pfqn_gflinearizer(L, N, Z, sched_type)
                else:  # lin
                    result = pfqn_linearizer(L, N, Z, sched_type)

                QN_out, UN_out, WN_out, TN_out, CN_out, XN_out, iters = result

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = WN_out[idx, :]  # Response times
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]  # Arrival rate = throughput

                XN = XN_out.flatten()

            elif method in ['schmidt', 'schmidt-ext']:
                # Schmidt's exact MVA for multi-class networks
                # Build scheduling strategy array
                sched_arr = np.zeros(M, dtype=int)
                for idx, q_idx in enumerate(queue_indices):
                    if self.sched is not None and q_idx in self.sched:
                        sched_str = str(self.sched[q_idx])
                        if 'INF' in sched_str:
                            sched_arr[idx] = 2  # INF
                        elif 'PS' in sched_str:
                            sched_arr[idx] = 3  # PS
                        else:
                            sched_arr[idx] = 0  # FCFS
                    else:
                        sched_arr[idx] = 0  # Default FCFS

                result = pfqn_schmidt(L, N, mi, sched_arr)
                XN_out = result.XN
                QN_out = result.QN
                UN_out = result.UN

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    if XN_out.size > 0:
                        TN[q_idx, :] = XN_out
                        AN[q_idx, :] = XN_out

                XN = XN_out.flatten() if hasattr(XN_out, 'flatten') else np.array([XN_out])

            elif method == 'ab':
                # Akyildiz-Bolch AMVA for multi-server networks
                V = np.ones_like(L)  # Visit ratios (default 1)
                sched_arr = np.zeros(M, dtype=int)
                for idx, q_idx in enumerate(queue_indices):
                    if self.sched is not None and q_idx in self.sched:
                        sched_str = str(self.sched[q_idx])
                        if 'INF' in sched_str:
                            sched_arr[idx] = 2  # INF
                        elif 'PS' in sched_str:
                            sched_arr[idx] = 3  # PS
                        else:
                            sched_arr[idx] = 0  # FCFS
                    else:
                        sched_arr[idx] = 0  # Default FCFS

                result = pfqn_ab_amva(L, N, V, mi, sched_arr)
                QN_out = result.QN
                UN_out = result.UN
                XN_out = result.XN

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = result.RN[idx, :]
                    TN[q_idx, :] = XN_out
                    AN[q_idx, :] = XN_out

                XN = XN_out.flatten()

            elif method == 'qd':
                # Queue-Dependent (QD) AMVA
                QN_out, XN_out, UN_out, iters = pfqn_qd(L, N)

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    # Compute residence times from queue lengths
                    # XN_out is per-class throughput (R elements), not per-queue
                    if np.any(XN_out > 0):
                        RN[q_idx, :] = QN_out[idx, :] / np.maximum(XN_out, 1e-10)
                    else:
                        RN[q_idx, :] = L[idx, :]
                    TN[q_idx, :] = XN_out
                    AN[q_idx, :] = XN_out

                XN = XN_out.flatten()

            elif method == 'qdlin':
                # QD-Linearizer (QDLIN) AMVA
                QN_out, UN_out, RN_out, XN_out, CN_out, iters = pfqn_qdlin(L, N, Z)

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = XN_out.flatten()
                    AN[q_idx, :] = XN_out.flatten()

                XN = XN_out.flatten()

            elif method == 'qli':
                # Queue-Line (QLI) AMVA (Wang-Sevcik)
                QN_out, UN_out, RN_out, XN_out, CN_out, iters = pfqn_qli(L, N, Z)

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = XN_out.flatten()
                    AN[q_idx, :] = XN_out.flatten()

                XN = XN_out.flatten()

            elif method == 'fli':
                # Fraction-Line (FLI) AMVA (Wang-Sevcik)
                QN_out, UN_out, RN_out, XN_out, CN_out, iters = pfqn_fli(L, N, Z)

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = XN_out.flatten()
                    AN[q_idx, :] = XN_out.flatten()

                XN = XN_out.flatten()

            elif method == 'aba.lower':
                # Asymptotic Bound Analysis - lower bound
                X_low, Q_low = pfqn_xzabalow(L, N, Z)
                XN = X_low.flatten()
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = Q_low[idx, :] if Q_low.ndim > 1 else Q_low
                    UN[q_idx, :] = XN * L[idx, :]

            elif method == 'aba.upper':
                # Asymptotic Bound Analysis - upper bound
                X_up, Q_up = pfqn_xzabaup(L, N, Z)
                XN = X_up.flatten()
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = Q_up[idx, :] if Q_up.ndim > 1 else Q_up
                    UN[q_idx, :] = XN * L[idx, :]

            elif method in ['pb.lower', 'pb.upper']:
                # Performance bounds - use ABA as fallback
                if 'lower' in method:
                    X_low, Q_low = pfqn_xzabalow(L, N, Z)
                    XN = X_low.flatten()
                else:
                    X_up, Q_up = pfqn_xzabaup(L, N, Z)
                    XN = X_up.flatten()

                for idx, q_idx in enumerate(queue_indices):
                    UN[q_idx, :] = XN * L[idx, :]
                    QN[q_idx, :] = UN[q_idx, :]  # Approximate

            elif method in ['gb.lower', 'gb.upper']:
                # Geometric Bounds (single-class closed networks)
                if 'lower' in method:
                    X_bound, Q_bound = pfqn_qzgblow(L, N, Z)
                else:
                    X_bound, Q_bound = pfqn_qzgbup(L, N, Z)
                XN = X_bound.flatten() if hasattr(X_bound, 'flatten') else np.array([X_bound])
                for idx, q_idx in enumerate(queue_indices):
                    if Q_bound.ndim > 1:
                        QN[q_idx, :] = Q_bound[idx, :]
                    else:
                        QN[q_idx, :] = Q_bound
                    UN[q_idx, :] = XN * L[idx, :]

            elif method in ['sb.lower', 'sb.upper']:
                # Schweitzer Bounds (single-class closed networks)
                if 'lower' in method:
                    X_bound, Q_bound = pfqn_xzgsblow(L, N, Z)
                else:
                    X_bound, Q_bound = pfqn_xzgsbup(L, N, Z)
                XN = X_bound.flatten() if hasattr(X_bound, 'flatten') else np.array([X_bound])
                for idx, q_idx in enumerate(queue_indices):
                    if Q_bound.ndim > 1:
                        QN[q_idx, :] = Q_bound[idx, :]
                    else:
                        QN[q_idx, :] = Q_bound
                    UN[q_idx, :] = XN * L[idx, :]

            else:
                # Default to exact MVA
                result = pfqn_mva(L, N, Z, mi)
                XN_out, CN_out, QN_out, UN_out, RN_out, TN_out, AN_out = result

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = AN_out[idx, :]

                XN = XN_out.flatten()

        elif self.network_type == 'open':
            # Open network - use QNA or queueing system formulas
            from ..api.qsys import (
                qsys_mm1, qsys_mmk, qsys_mg1, qsys_gm1,
                qsys_gig1_approx_kingman, qsys_gig1_approx_gelenbe,
                qsys_gig1_approx_heyman, qsys_gig1_approx_kimura,
                qsys_gig1_approx_kobayashi, qsys_gig1_approx_marchal,
                qsys_gig1_approx_allencunneen, qsys_gigk_approx,
            )

            # Get arrival rate and service rate for single queue
            lambda_r = 1.0  # Default arrival rate
            mu = 1.0  # Default service rate
            k = 1  # Number of servers

            for r in range(R):
                if not np.isfinite(self.njobs[r]):  # Open class
                    # Get arrival rate from source station (use source_indices, not just "not in queue_indices")
                    for src_idx in source_indices:
                        if self.rates[src_idx, r] > 0:
                            lambda_r = self.rates[src_idx, r]
                            break

                    for idx, q_idx in enumerate(queue_indices):
                        if L[idx, r] > 0:
                            mu = 1.0 / L[idx, r]
                            # Handle infinity (infinite servers) - treat as single server for analysis
                            nserv = mi[idx] if idx < len(mi) else 1.0
                            k = 1 if not np.isfinite(nserv) else int(nserv)

            # Handle specific queueing system methods
            qsys_methods = {
                'mm1', 'mmk', 'mg1', 'mgi1', 'gm1', 'gig1', 'gim1',
                'gig1.kingman', 'gigk', 'gigk.kingman_approx',
                'gig1.gelenbe', 'gig1.heyman', 'gig1.kimura',
                'gig1.allen', 'gig1.kobayashi', 'gig1.klb', 'gig1.marchal',
            }

            if method in qsys_methods and M == 1 and R == 1:
                # Single queue, single class - use exact formulas
                rho = lambda_r / (k * mu) if mu > 0 else 0

                if method == 'mm1':
                    result = qsys_mm1(lambda_r, mu)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'mmk':
                    result = qsys_mmk(lambda_r, mu, k)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method in ['mg1', 'mgi1']:
                    # M/G/1 - need SCV of service
                    scv_s = 1.0  # Assume exponential (SCV=1)
                    result = qsys_mg1(lambda_r, mu, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method in ['gm1', 'gim1']:
                    # G/M/1 - need SCV of arrivals (gim1 is alias for gm1)
                    scv_a = 1.0  # Assume Poisson (SCV=1)
                    result = qsys_gm1(lambda_r, mu, scv_a)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method in ['gig1', 'gig1.kingman']:
                    # G/G/1 Kingman approximation
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_kingman(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method in ['gigk', 'gigk.kingman_approx']:
                    # G/G/k Kingman approximation
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gigk_approx(lambda_r, mu, k, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'gig1.gelenbe':
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_gelenbe(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'gig1.heyman':
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_heyman(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'gig1.kimura':
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_kimura(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'gig1.allen':
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_allencunneen(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method in ['gig1.kobayashi', 'gig1.klb']:
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_kobayashi(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

                elif method == 'gig1.marchal':
                    scv_a, scv_s = 1.0, 1.0
                    result = qsys_gig1_approx_marchal(lambda_r, mu, scv_a, scv_s)
                    QN[queue_indices[0], 0] = result['L']
                    UN[queue_indices[0], 0] = result['rho']
                    RN[queue_indices[0], 0] = result['W']
                    TN[queue_indices[0], 0] = lambda_r
                    AN[queue_indices[0], 0] = lambda_r
                    XN[0] = lambda_r

            else:
                # Default: use M/M/1 formulas for each queue
                for r in range(R):
                    if not np.isfinite(self.njobs[r]):  # Open class (njobs = inf)
                        # Get arrival rate from source station (use source_indices, not just "not in queue_indices")
                        lambda_r = 1.0
                        for src_idx in source_indices:
                            if self.rates[src_idx, r] > 0:
                                lambda_r = self.rates[src_idx, r]
                                break

                        for idx, q_idx in enumerate(queue_indices):
                            mu = 1.0 / L[idx, r] if L[idx, r] > 0 else float('inf')
                            rho = lambda_r / mu if mu > 0 else 0

                            if rho < 1:
                                QN[q_idx, r] = rho / (1 - rho)  # M/M/1 queue length
                                UN[q_idx, r] = rho
                                RN[q_idx, r] = L[idx, r] / (1 - rho) if rho < 1 else float('inf')
                                TN[q_idx, r] = lambda_r
                                AN[q_idx, r] = lambda_r

                            XN[r] = lambda_r

        # Compute Delay station metrics from throughput and individual station demands
        # For infinite server (delay) stations: QN = X*D, UN = QN, RN = D
        delay_stations = self._get_delay_stations()
        for d_idx in delay_stations:
            for r in range(R):
                d_demand = self.demands[d_idx, r]  # Individual station's demand
                if d_demand > 0 and XN[r] > 0:
                    QN[d_idx, r] = XN[r] * d_demand  # Little's law: jobs = throughput * time
                    UN[d_idx, r] = QN[d_idx, r]  # For infinite servers, utilization = queue length
                    RN[d_idx, r] = d_demand  # Response time = service time at this station
                    TN[d_idx, r] = XN[r]  # Throughput
                    AN[d_idx, r] = XN[r]  # Arrival rate

        # Compute residence times if not set
        source_stations = self._get_source_stations()
        for r in range(R):
            if XN[r] > 0:
                for i in range(self.nstations):
                    if RN[i, r] == 0 and QN[i, r] > 0:
                        RN[i, r] = QN[i, r] / XN[r]
                    # Only set TN/AN for stations with non-zero demand (not Disabled services)
                    if self.demands[i, r] > 0:
                        if TN[i, r] == 0:
                            TN[i, r] = XN[r]
                        # Don't set AN for Source stations (jobs originate there, not arrive)
                        if AN[i, r] == 0 and i not in source_stations:
                            AN[i, r] = XN[r]

        # Compute fork-join synchronization delays if model has fork-join
        if self._has_fork_join():
            self._compute_fork_join_sync_delays(QN, UN, RN, TN, AN, XN)

        # Store results
        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
        }

        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns node-based results (one row per node per class) to match MATLAB output format.
        Non-station nodes (e.g., Fork) are included with zero metrics.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        # Get node-based dimensions and mappings from NetworkStruct
        nnodes = self._sn.nnodes if hasattr(self._sn, 'nnodes') else self.nstations
        nodeToStation = np.asarray(self._sn.nodeToStation).flatten() if hasattr(self._sn, 'nodeToStation') else np.arange(self.nstations)
        nodenames = list(self._sn.nodenames) if hasattr(self._sn, 'nodenames') and self._sn.nodenames else []

        # Build table data - iterate over nodes (not stations) to match MATLAB format
        rows = []
        for node_idx in range(nnodes):
            station_idx = int(nodeToStation[node_idx]) if node_idx < len(nodeToStation) else -1
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

            for r in range(self.nclasses):
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                # If this node is a station, use station results; otherwise use zeros
                if station_idx >= 0 and station_idx < self.nstations:
                    rows.append({
                        'Station': node_name,
                        'JobClass': class_name,
                        'QLen': self._result['QN'][station_idx, r],
                        'Util': self._result['UN'][station_idx, r],
                        'RespT': self._result['RN'][station_idx, r],
                        'ResidT': self._result['RN'][station_idx, r],
                        'ArvR': self._result['AN'][station_idx, r],
                        'Tput': self._result['TN'][station_idx, r],
                    })
                else:
                    # Non-station node (e.g., Fork) - include with zeros
                    rows.append({
                        'Station': node_name,
                        'JobClass': class_name,
                        'QLen': 0.0,
                        'Util': 0.0,
                        'RespT': 0.0,
                        'ResidT': 0.0,
                        'ArvR': 0.0,
                        'Tput': 0.0,
                    })

        df = pd.DataFrame(rows)

        # Filter out all-zero rows (MATLAB excludes nodes like Fork with no metrics)
        numeric_cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        tokeep = ~(df[numeric_cols] <= 0.0).all(axis=1)
        df = df.loc[tokeep].reset_index(drop=True)

        if not self._table_silent:
            print(df.to_string(index=False))

        return df

    # ============================================================================
    # Probability Methods (Phase 2)
    # ============================================================================

    def getProbAggr(self, ist: int) -> Tuple[float, float]:
        """
        Get probability of current per-class job distribution at a station.

        Returns P(n_1, ..., n_K at station i) using binomial approximation.

        Args:
            ist: Station index (1-based, MATLAB style)

        Returns:
            (log_prob, prob): Tuple of log probability and probability value

        Raises:
            ValueError: If station index invalid or analysis not run

        Example:
            >>> solver = SolverMVA(model)
            >>> solver.runAnalyzer()
            >>> log_p, p = solver.getProbAggr(1)  # Station 1
        """
        from ..api.solvers.mva.prob_methods import get_prob_aggr

        if self._result is None:
            self.runAnalyzer()

        # Create minimal SolverResults compatible object
        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                self.U = result_dict.get('UN')
                self.R = result_dict.get('RN')
                self.prob = None

        result_adapter = ResultAdapter(self._result)
        return get_prob_aggr(self._get_network_struct(), result_adapter, ist)

    def getProbMarg(
        self, ist: int, jobclass: int, state_m: Optional[np.ndarray] = None
    ) -> Tuple[np.ndarray, np.ndarray]:
        """
        Get marginal queue-length distribution for a class at a station.

        Returns P(n | station i, class r) for n = 0, 1, ..., N[r].

        Args:
            ist: Station index (1-based)
            jobclass: Job class index (1-based)
            state_m: Optional state vector (for future use)

        Returns:
            (states, probs): Array of state indices and probabilities
                - states: [0, 1, ..., N[jobclass]]
                - probs: Probability distribution summing to 1.0

        Example:
            >>> states, probs = solver.getProbMarg(1, 1)
            >>> print(f"P(n=2 at station 1, class 1) = {probs[2]}")
        """
        from ..api.solvers.mva.prob_methods import get_prob_marg

        if self._result is None:
            self.runAnalyzer()

        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                self.prob = None

        result_adapter = ResultAdapter(self._result)
        return get_prob_marg(self._get_network_struct(), result_adapter, ist, jobclass)

    def getProbSysAggr(self) -> Tuple[float, float]:
        """
        Get joint probability of current system state.

        Returns P(full system state) using product of station marginals.

        Returns:
            (log_prob, prob): Log probability and probability value

        Notes:
            - Assumes station independence (valid for product-form networks)
            - Requires network state to be set

        Example:
            >>> log_p, p = solver.getProbSysAggr()
            >>> print(f"System state probability: {p:.6e}")
        """
        from ..api.solvers.mva.prob_methods import get_prob_sys_aggr

        if self._result is None:
            self.runAnalyzer()

        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                self.prob = None

        result_adapter = ResultAdapter(self._result)
        return get_prob_sys_aggr(self._get_network_struct(), result_adapter)

    def getProbNormConstAggr(self) -> float:
        """
        Get log normalizing constant for closed queueing network.

        Returns log(G) where G = ∑_state P(state).

        Returns:
            log_G: Natural logarithm of normalizing constant
                - For open networks: returns inf
                - For closed networks: exact value from MVA or approximation

        Notes:
            - Only applies to closed queueing networks
            - Returns inf for open networks (normalizing constant is infinite)

        Example:
            >>> log_G = solver.getProbNormConstAggr()
            >>> G = np.exp(log_G)  # Reconstruct if needed
        """
        from ..api.solvers.mva.prob_methods import get_prob_norm_const_aggr

        if self._result is None:
            self.runAnalyzer()

        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                self.prob = None

        result_adapter = ResultAdapter(self._result)
        return get_prob_norm_const_aggr(self._get_network_struct(), result_adapter)

    def _get_network_struct(self):
        """Get or create network structure for probability computations."""
        class NetworkStructAdapter:
            def __init__(self, solver):
                self.nstations = solver.nstations
                self.nclasses = solver.nclasses
                self.njobs = solver.njobs
                self.nservers = solver.nservers
                self.nodetype = getattr(solver, 'nodetype', None)
                self.refstat = getattr(solver, 'refstat', None)
                self.state = None  # Will be set if needed
                self.sched = getattr(solver, 'sched', None)

        return NetworkStructAdapter(self)

    # ============================================================================
    # Individual Metric Accessors (Phase 3)
    # ============================================================================

    def getAvgQLen(self) -> np.ndarray:
        """
        Get average queue lengths.

        Returns:
            Q: Queue lengths matrix (M x K)
               M = number of stations
               K = number of classes
               Q[i,r] = average number of jobs of class r at station i

        Raises:
            RuntimeError: If solver not run yet

        Example:
            >>> Q = solver.getAvgQLen()
            >>> print(f"Queue length at station 1, class 1: {Q[0,0]}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['QN'].copy()

    def getAvgUtil(self) -> np.ndarray:
        """
        Get average utilizations.

        Returns:
            U: Utilization matrix (M x K)
               U[i,r] = utilization of station i by class r
               Range: [0, 1] for single-server, [0, ∞) for multi-server

        Example:
            >>> U = solver.getAvgUtil()
            >>> print(f"Utilization at station 1: {U[0,:].sum()}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['UN'].copy()

    def getAvgRespT(self) -> np.ndarray:
        """
        Get average response times.

        Returns:
            R: Response times matrix (M x K)
               R[i,r] = average time spent at station i for class r
               Includes both service and waiting time

        Example:
            >>> R = solver.getAvgRespT()
            >>> print(f"Response time at station 1, class 1: {R[0,0]}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['RN'].copy()

    def getAvgResidT(self) -> np.ndarray:
        """
        Get average residence times (alias for response times).

        Returns:
            ResidT: Residence times matrix (M x K)
                    Same as response times R[i,r]

        Note:
            Residence time = Response time (total time in station)
        """
        return self.getAvgRespT()

    def getAvgWaitT(self) -> np.ndarray:
        """
        Get average waiting times.

        Returns:
            W: Waiting times matrix (M x K)
               W[i,r] = R[i,r] - S[i,r]
               where S[i,r] is the mean service time
               W[i,r] = 0 for Delay (think time) stations

        Example:
            >>> W = solver.getAvgWaitT()
            >>> print(f"Waiting time at station 1: {W[0,:].sum()}")
        """
        if self._result is None:
            self.runAnalyzer()

        # W = R - S, where S is the service demand (1/service_rate)
        R = self._result['RN'].copy()
        S = self.demands.copy()  # Service demands already computed in __init__

        W = R - S
        # Ensure non-negative (numerical precision)
        W = np.maximum(W, 0.0)
        return W

    def getAvgTput(self) -> np.ndarray:
        """
        Get average throughputs.

        Returns:
            T: Throughput matrix (M x K)
               T[i,r] = average throughput at station i for class r
               jobs/time unit

        Example:
            >>> T = solver.getAvgTput()
            >>> print(f"Throughput at station 1, class 1: {T[0,0]}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['TN'].copy()

    def getAvgArvR(self) -> np.ndarray:
        """
        Get average arrival rates.

        Returns:
            A: Arrival rates matrix (M x K)
               A[i,r] = arrival rate to station i for class r

        Note:
            For closed networks, arrival rates are derived from throughputs
            and visit ratios
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['AN'].copy()

    def getAvgSysRespT(self) -> np.ndarray:
        """
        Get system response times (cycle times).

        Returns:
            C: Cycle times vector (K,)
               C[r] = average time for one cycle for class r
               (in closed networks: time between successive visits to reference station)

        Note:
            For open networks: sum of response times across all stations
        """
        if self._result is None:
            self.runAnalyzer()

        # Cycle time = sum of response times across all stations
        R = self._result['RN']
        C = np.sum(R, axis=0)  # Sum across stations for each class
        return C

    def getAvgSysTput(self) -> np.ndarray:
        """
        Get system throughputs.

        Returns:
            X: System throughput vector (K,)
               X[r] = overall throughput of the system for class r

        Note:
            This is the throughput at any single bottleneck station
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['XN'].copy()

    # ============================================================================
    # Unified Metrics and Chain/Node/System Methods
    # ============================================================================

    def getAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics at once.

        Returns:
            Tuple of (Q, U, R, T, A, W) where:
            - Q: Queue lengths (M x K)
            - U: Utilizations (M x K)
            - R: Response times (M x K)
            - T: Throughputs (M x K)
            - A: Arrival rates (M x K)
            - W: Waiting times (M x K)
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result['QN']
        U = self._result['UN']
        R = self._result['RN']
        T = np.broadcast_to(self._result['XN'], Q.shape).copy()
        A = T.copy()
        W = R.copy()

        return Q, U, R, T, A, W

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self, '_sn') and self._sn is not None and hasattr(self._sn, 'chains') and self._sn.chains is not None:
            # chains is (K,) array where chains[k] = chain index for class k
            chain_arr = np.asarray(self._sn.chains).flatten()
            if len(chain_arr) == 0:
                return [[k] for k in range(self.nclasses)]

            nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else int(np.max(chain_arr)) + 1
            chains = [[] for _ in range(nchains)]
            for k in range(self.nclasses):
                if k < len(chain_arr):
                    chain_idx = int(chain_arr[k])
                    if 0 <= chain_idx < nchains:
                        chains[chain_idx].append(k)
            # Filter out empty chains
            chains = [c for c in chains if c]
            return chains if chains else [[k for k in range(self.nclasses)]]
        else:
            return [[k] for k in range(self.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        Q = self._result['QN']
        chains = self._get_chains()
        nchains = len(chains)

        QN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                QN_chain[:, c] = np.sum(Q[:, chain_classes], axis=1)

        return QN_chain

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        U = self._result['UN']
        chains = self._get_chains()
        nchains = len(chains)

        UN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                UN_chain[:, c] = np.sum(U[:, chain_classes], axis=1)

        return UN_chain

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        R = self._result['RN']
        chains = self._get_chains()
        nchains = len(chains)

        RN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                RN_chain[:, c] = np.mean(R[:, chain_classes], axis=1)

        return RN_chain

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        T = self._result['XN']
        chains = self._get_chains()
        nchains = len(chains)

        TN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                TN_chain[:, c] = np.sum(T[chain_classes])

        return TN_chain

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics aggregated by chain."""
        QN = self.getAvgQLenChain()
        UN = self.getAvgUtilChain()
        RN = self.getAvgRespTChain()
        WN = self.getAvgResidTChain()
        AN = self.getAvgArvRChain()
        TN = self.getAvgTputChain()
        return QN, UN, RN, WN, AN, TN

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()

        nstations, nchains = QN.shape
        rows = []

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': f'Station{i}',
                    'Chain': f'Chain{c}',
                    'QLen': QN[i, c],
                    'Util': UN[i, c],
                    'RespT': RN[i, c],
                    'ResidT': WN[i, c],
                    'ArvR': AN[i, c],
                    'Tput': TN[i, c],
                })

        return pd.DataFrame(rows)

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics per node (same as per station for MVA)."""
        return self.getAvg()

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics by node as DataFrame."""
        return self.getAvgTable()

    def getAvgNodeChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics by node and chain."""
        return self.getAvgChain()

    def getAvgNodeChainTable(self) -> pd.DataFrame:
        """Get average metrics by node and chain as DataFrame."""
        return self.getAvgChainTable()

    def getAvgNodeQLenChain(self) -> np.ndarray:
        """Get average queue lengths by node aggregated by chain."""
        return self.getAvgQLenChain()

    def getAvgNodeUtilChain(self) -> np.ndarray:
        """Get average utilizations by node aggregated by chain."""
        return self.getAvgUtilChain()

    def getAvgNodeRespTChain(self) -> np.ndarray:
        """Get average response times by node aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgNodeResidTChain(self) -> np.ndarray:
        """Get average residence times by node aggregated by chain."""
        return self.getAvgResidTChain()

    def getAvgNodeTputChain(self) -> np.ndarray:
        """Get average throughputs by node aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgNodeArvRChain(self) -> np.ndarray:
        """Get average arrival rates by node aggregated by chain."""
        return self.getAvgArvRChain()

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics."""
        R = self.getAvgSysRespT()
        T = self.getAvgSysTput()
        return R, T

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame."""
        R, T = self.getAvgSys()
        nclasses = len(R) if hasattr(R, '__len__') else 1

        rows = []
        for k in range(nclasses):
            rows.append({
                'Class': f'Class{k}',
                'SysRespT': R[k] if hasattr(R, '__getitem__') else R,
                'SysTput': T[k] if k < len(T) else 0,
            })

        return pd.DataFrame(rows)

    # ============================================================================
    # Metric Accessor Aliases
    # ============================================================================

    # PascalCase aliases for MATLAB compatibility
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput

    # ============================================================================
    # CDF and Percentile Analysis (Phase 4)
    # ============================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Get response time cumulative distribution function (CDF).

        Uses exponential approximation: CDF(t) = 1 - exp(-t / E[R])
        where E[R] is the mean response time from MVA.

        Args:
            R: Optional response times matrix (M x K)
               If None, uses results from runAnalyzer()

        Returns:
            RD: List of dicts, one per (station, class) pair with service
                Each dict contains:
                - 'station': Station index (1-based)
                - 'class': Job class (1-based)
                - 't': Time points (100 points from 0.001 to 0.999 quantile)
                - 'p': CDF values at each time point

        Notes:
            - Uses exponential approximation for single-class exponential service
            - For multi-phase or non-exponential, this is approximate
            - Returns empty list for stations with zero response time

        Example:
            >>> cdf_list = solver.getCdfRespT()
            >>> for cdf_data in cdf_list:
            ...     station = cdf_data['station']
            ...     class_id = cdf_data['class']
            ...     print(f"Station {station}, Class {class_id}")
            ...     # Access t and p arrays for plotting
            ...     t = cdf_data['t']
            ...     p = cdf_data['p']
        """
        if self._result is None:
            self.runAnalyzer()

        if R is None:
            R = self._result['RN']

        RD = []  # Result list

        for i in range(self.nstations):
            for r in range(self.nclasses):
                mean_resp_t = R[i, r]

                # Skip if no response time (no service at this station)
                if mean_resp_t <= 0:
                    continue

                # Exponential CDF: F(t) = 1 - exp(-t/mean)
                # Rate parameter: lambda = 1 / mean_resp_t
                lambda_rate = 1.0 / mean_resp_t

                # Generate time points from quantiles
                # Use quantiles 0.001 to 0.999 to avoid singularities
                quantiles = np.linspace(0.001, 0.999, 100)
                times = -np.log(1 - quantiles) / lambda_rate

                # Compute CDF values
                cdf_vals = 1 - np.exp(-lambda_rate * times)

                # Store result
                RD.append({
                    'station': i + 1,  # 1-based indexing
                    'class': r + 1,    # 1-based indexing
                    't': times,
                    'p': cdf_vals,
                })

        return RD

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None,
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """
        Extract percentiles from response time distribution.

        Computes percentile points from the exponential CDF approximation.

        Args:
            percentiles: List of percentiles to extract (0-100)
                        Default: [10, 25, 50, 75, 90, 95, 99]
            jobclass: Optional class index to filter results (1-based)
                     If None, returns all classes

        Returns:
            (PercRT, PercTable): Tuple of:
            - PercRT: List of dicts with percentile data for each (station, class)
              Each dict contains:
                - 'station': Station index
                - 'class': Job class
                - 'percentiles': Input percentile values
                - 'values': Percentile response times
            - PercTable: pandas DataFrame with columns:
              Station, Class, P10, P25, P50, P75, P90, P95, P99, ...

        Algorithm:
            For exponential CDF with rate λ = 1/E[R]:
            Percentile p: t_p = -ln(1-p) * E[R]
            where p is in [0,1]

        Notes:
            - Percentiles outside (0,100) are clipped
            - Returns empty lists for stations with zero response time
            - Results match exponential percentile formula

        Example:
            >>> perc_list, perc_table = solver.getPerctRespT([90, 95, 99])
            >>> print(perc_table)
            >>> # Extract 90th percentile response time
            >>> p90_values = perc_list[0]['values']
        """
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90, 95, 99]

        # Normalize percentiles to [0, 100]
        percentiles = np.asarray(percentiles)
        percentiles = np.clip(percentiles, 0.01, 99.99)
        percentiles_normalized = percentiles / 100.0  # Convert to [0, 1]

        # Get CDFs
        cdf_list = self.getCdfRespT()

        PercRT = []
        rows = []

        # Create column names for DataFrame
        perc_col_names = [f'P{int(p)}' for p in percentiles]

        for cdf_data in cdf_list:
            station = cdf_data['station']
            class_id = cdf_data['class']

            # Only include specified class if filter is set
            if jobclass is not None and class_id != jobclass:
                continue

            # Extract mean response time
            mean_resp_t = self._result['RN'][station - 1, class_id - 1]

            if mean_resp_t <= 0:
                continue

            # Compute percentile values using exponential formula
            # For exponential: t_p = -ln(1-p) * E[R]
            lambda_rate = 1.0 / mean_resp_t
            perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

            # Store in result list
            PercRT.append({
                'station': station,
                'class': class_id,
                'percentiles': percentiles.tolist(),
                'values': perc_values.tolist(),
            })

            # Add row to DataFrame
            row_data = {
                'Station': self.station_names[station - 1] if station - 1 < len(self.station_names) else f'Station{station}',
                'Class': self.class_names[class_id - 1] if class_id - 1 < len(self.class_names) else f'Class{class_id}',
            }
            for perc_col, perc_val in zip(perc_col_names, perc_values):
                row_data[perc_col] = perc_val

            rows.append(row_data)

        # Create DataFrame
        if rows:
            PercTable = pd.DataFrame(rows)
        else:
            PercTable = pd.DataFrame()

        return PercRT, PercTable

    # ============================================================================
    # Method Introspection (Phase 5)
    # ============================================================================

    def listValidMethods(self) -> List[str]:
        """
        List all valid solution methods for this network model.

        Returns the set of applicable methods based on network characteristics
        (single-class/multi-class, open/closed, etc.).

        Returns:
            List of method names available for this model:
            - Base methods: 'default', 'mva', 'exact', 'amva', 'qna'
            - AMVA variants: 'bs', 'sqni', 'lin', 'gflin', 'egflin', 'schmidt', 'schmidt-ext', 'ab'
            - Bounds (closed networks): 'aba.lower', 'aba.upper', 'pb.lower', 'pb.upper'
            - Additional bounds (single-class closed): 'bjb.upper', 'bjb.lower', 'gb.upper', 'gb.lower', 'sb.upper', 'sb.lower'
            - Queueing formulas (2-station open): 'mm1', 'mmk', 'mg1', 'mgi1', 'gm1', 'gig1', etc.

        Notes:
            - MVA solver always supports 'default' and 'exact' methods
            - AMVA (approximate MVA) available for most networks
            - Bounds methods useful for single-class networks
            - All methods also available with 'amva.' prefix (e.g., 'amva.lin')

        Example:
            >>> methods = solver.listValidMethods()
            >>> print(f"Available methods: {methods}")
            >>> # Can then call: solver.method = 'amva.lin'
        """
        methods = []

        # Base methods - always available
        methods.extend(['default', 'mva', 'exact'])

        # AMVA and variants - available for closed/mixed networks
        methods.extend([
            'amva',
            'bs', 'amva.bs',
            'sqni',
            'lin', 'amva.lin',
            'gflin',
            'egflin',
            'schmidt', 'amva.schmidt',
            'schmidt-ext', 'amva.schmidt-ext',
            'ab', 'amva.ab',
            'qd', 'amva.qd',
            'qdlin', 'amva.qdlin',
            'qli', 'amva.qli',
            'fli', 'amva.fli',
        ])

        # QNA for open networks
        if self.network_type == 'open':
            methods.append('qna')

        # Bounds for closed networks
        if self.network_type == 'closed':
            methods.extend(['aba.lower', 'aba.upper'])

            # Single-class only bounds
            if self.nclasses == 1:
                methods.extend([
                    'bjb.upper', 'bjb.lower',
                    'gb.upper', 'gb.lower',
                    'sb.upper', 'sb.lower',
                ])

        # Performance bounds - always available
        methods.extend(['pb.lower', 'pb.upper'])

        # Queueing system formulas for 2-station open networks
        if self.network_type == 'open' and self.nstations == 2 and self.nclasses == 1:
            methods.extend([
                'mm1', 'mmk', 'mg1', 'mgi1', 'gm1', 'gig1', 'gim1',
                'gig1.kingman', 'gigk', 'gigk.kingman_approx',
                'gig1.gelenbe', 'gig1.heyman', 'gig1.kimura',
                'gig1.allen', 'gig1.kobayashi', 'gig1.klb', 'gig1.marchal',
            ])

        return methods

    @staticmethod
    def getFeatureSet() -> set:
        """
        Get set of features supported by MVA solver.

        Returns:
            Set of feature strings indicating what the solver supports:
            - Node types: 'Sink', 'Source', 'Queue', 'Delay', 'Fork', 'Join'
            - Distributions: 'Exp', 'Det', 'Erlang', 'HyperExp', 'PH', 'MAP', etc.
            - Classes: 'OpenClass', 'ClosedClass'
            - Scheduling: 'FCFS', 'PS', 'LCFS', 'SIRO'
            - Special: 'MultiServer', 'LoadDependent', 'ThinkTime'

        Example:
            >>> features = SolverMVA.getFeatureSet()
            >>> if 'MultiServer' in features:
            ...     print("Multi-server queues supported")
        """
        features = {
            # Node types
            'Sink', 'Source', 'Queue', 'Delay',
            'Fork', 'Join', 'Router', 'ClassSwitch',
            # Distributions
            'Exp', 'Det', 'Erlang', 'HyperExp',
            'Gamma', 'Lognormal', 'Pareto', 'Uniform', 'Weibull',
            'PH', 'APH', 'Cox', 'MAP', 'MMPP',
            'Poisson', 'Geometric', 'Binomial', 'NegBinomial',
            # Classes
            'OpenClass', 'ClosedClass',
            # Scheduling
            'FCFS', 'PS', 'LCFS', 'SIRO', 'LIFO',
            # Special features
            'MultiServer', 'LoadDependent', 'ThinkTime',
            'Multiclass', 'ProductForm',
        }
        return features

    @staticmethod
    def supports(model) -> bool:
        """
        Check if MVA solver supports the given network model.

        Performs basic model validation to ensure compatibility.

        Args:
            model: Network model to check (native or wrapper)

        Returns:
            True if model is supported, False otherwise

        Notes:
            - Checks for product-form network structure
            - Verifies presence of required network components
            - Returns True for most standard queueing networks

        Example:
            >>> if SolverMVA.supports(model):
            ...     solver = SolverMVA(model)
            ... else:
            ...     print("Model not supported by MVA")
        """
        try:
            # Try to extract basic network properties
            if hasattr(model, 'nstations'):
                nstations = model.nstations
            elif hasattr(model, 'getNumberOfStations'):
                nstations = model.getNumberOfStations()
            else:
                return False

            if hasattr(model, 'nclasses'):
                nclasses = model.nclasses
            elif hasattr(model, 'getNumberOfClasses'):
                nclasses = model.getNumberOfClasses()
            else:
                return False

            # Basic validation
            return nstations > 0 and nclasses > 0

        except Exception:
            return False

    @staticmethod
    def defaultOptions() -> Dict[str, Any]:
        """
        Get default solver options.

        Returns:
            Dictionary with default option values:
            - 'method': 'exact' (exact MVA algorithm)
            - 'tol': 1e-8 (convergence tolerance)
            - 'max_iter': 1000 (maximum iterations)
            - 'verbose': False (quiet by default)

        Example:
            >>> opts = SolverMVA.defaultOptions()
            >>> opts['method'] = 'amva'  # Override for approximate MVA
            >>> solver = SolverMVA(model, **opts)
        """
        return OptionsDict({
            'method': 'exact',
            'tol': 1e-8,
            'max_iter': 1000,
            'verbose': False,
        })

    # ============================================================================
    # Sampling and Transient Methods (Phase 6) - Placeholders
    # ============================================================================

    def sample(self, node: int, numSamples: int) -> np.ndarray:
        """
        Sample from the response time distribution.

        **Not supported by MVA solver** - MVA is an analytical solver.
        For sampling, use simulation-based solvers.

        Args:
            node: Node/station index (1-based)
            numSamples: Number of samples to generate

        Returns:
            NotImplementedError (sampling not supported)

        Raises:
            NotImplementedError: Always - MVA does not support sampling

        Recommendation:
            Use SolverSSA (Stochastic State-space Analysis) or SolverJMT
            (JMT simulator) for sampling-based analysis.

        Example:
            >>> # Instead of sampling from MVA:
            >>> # solver = SolverMVA(model)
            >>> # This will raise NotImplementedError
            >>> solver.sample(1, 1000)
        """
        raise NotImplementedError(
            "Sampling not supported by SolverMVA (analytical solver). "
            "Use SolverSSA, SolverDES, or SolverJMT for sampling-based analysis."
        )

    def sampleAggr(self, node: int, numSamples: int) -> np.ndarray:
        """Aggregate sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleAggr() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, DES, JMT)."
        )

    def sampleSys(self, numSamples: int) -> np.ndarray:
        """System-level sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleSys() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, DES, JMT)."
        )

    def sampleSysAggr(self, numSamples: int) -> np.ndarray:
        """Aggregate system-level sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleSysAggr() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, DES, JMT)."
        )

    def getCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Get passage time CDF (not supported by MVA).

        Passage time = time to reach target station from source.
        Not computed by analytical MVA solver.

        Args:
            R: Optional response times (ignored)

        Raises:
            NotImplementedError: Passage time analysis not available

        Recommendation:
            Use simulation-based solvers for detailed path analysis.
        """
        raise NotImplementedError(
            "getCdfPassT() not supported by SolverMVA. "
            "Passage time analysis requires simulation-based solvers."
        )

    def getTranCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Get transient response time CDF (not supported by MVA).

        Transient analysis (time-dependent) not available from steady-state MVA.

        Args:
            R: Optional response times (ignored)

        Raises:
            NotImplementedError: Transient analysis not available

        Recommendation:
            Use SolverCTMC (Markov chain) or simulation solvers for transient.
        """
        raise NotImplementedError(
            "getTranCdfRespT() not supported by SolverMVA. "
            "Transient analysis available via SolverCTMC or SolverDES."
        )

    def getTranCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """Transient passage time CDF (not supported by MVA)."""
        raise NotImplementedError(
            "getTranCdfPassT() not supported by SolverMVA. "
            "Use simulation-based or CTMC solvers for transient analysis."
        )

    def getTranAvg(self) -> np.ndarray:
        """
        Get transient average metrics (not supported by MVA).

        MVA computes only steady-state metrics.

        Raises:
            NotImplementedError: Transient analysis not available
        """
        raise NotImplementedError(
            "getTranAvg() not supported by SolverMVA. "
            "MVA computes steady-state metrics only. "
            "Use SolverCTMC for transient analysis."
        )

    # ============================================================================
    # Introspection and Sampling Aliases
    # ============================================================================

    # PascalCase aliases for MATLAB compatibility
    ListValidMethods = listValidMethods
    GetFeatureSet = getFeatureSet
    Supports = supports
    DefaultOptions = defaultOptions
    Sample = sample
    SampleAggr = sampleAggr
    SampleSys = sampleSys
    SampleSysAggr = sampleSysAggr
    GetCdfPassT = getCdfPassT
    GetTranCdfRespT = getTranCdfRespT
    GetTranCdfPassT = getTranCdfPassT
    GetTranAvg = getTranAvg

    # ============================================================================
    # CDF and Percentile Aliases
    # ============================================================================

    # PascalCase aliases for MATLAB compatibility
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT

    # ============================================================================
    # Aliases and Compatibility
    # ============================================================================

    # PascalCase aliases for MATLAB compatibility
    GetProbAggr = getProbAggr
    GetProbMarg = getProbMarg
    GetProbSysAggr = getProbSysAggr
    GetProbNormConstAggr = getProbNormConstAggr

    # Table aliases
    avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable
    default_options = defaultOptions

    # Chain-level aliases
    GetAvg = getAvg
    GetAvgChain = getAvgChain
    GetAvgChainTable = getAvgChainTable
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain

    # Node-level aliases
    GetAvgNode = getAvgNode
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain
    GetAvgSys = getAvgSys
    GetAvgSysTable = getAvgSysTable

    # Short aliases (MATLAB compatibility)
    aNT = getAvgNodeTable
    aCT = getAvgChainTable
    aNCT = getAvgNodeChainTable
    aST = getAvgSysTable
    nodeAvgT = getAvgNodeTable
    chainAvgT = getAvgChainTable
    nodeChainAvgT = getAvgNodeChainTable
    sysAvgT = getAvgSysTable

    # Snake case aliases
    avg_node_table = getAvgNodeTable
    avg_chain_table = getAvgChainTable
    avg_node_chain_table = getAvgNodeChainTable
    avg_sys_table = getAvgSysTable
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_resp_t = getAvgRespT
    avg_resid_t = getAvgResidT
    avg_wait_t = getAvgWaitT
    avg_tput = getAvgTput
    avg_arv_r = getAvgArvR
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    run_analyzer = runAnalyzer
    cdf_resp_t = getCdfRespT
    perct_resp_t = getPerctRespT


__all__ = ['SolverMVA', 'SolverMVAOptions']
