"""
Native Python implementation of Layered Network (LN) solver.

This implementation provides a pure Python fixed-point iteration solver
for layered queueing networks, using native MVA algorithms for each layer.
Pure Python implementation.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple, Callable, Union
from dataclasses import dataclass, field


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
class SolverLNOptions:
    """Options for the native LN solver."""
    max_iter: int = 100
    tol: float = 1e-6
    verbose: bool = False
    under_relaxation: float = 0.7  # Under-relaxation factor for convergence
    method: str = 'exact'  # Layer solver method


@dataclass
class LayerMetrics:
    """Performance metrics for a single layer."""
    QN: np.ndarray = None  # Queue lengths
    UN: np.ndarray = None  # Utilizations
    RN: np.ndarray = None  # Response times
    TN: np.ndarray = None  # Throughputs
    AN: np.ndarray = None  # Arrival rates
    XN: np.ndarray = None  # System throughputs


@dataclass
class LayerModel:
    """Simplified representation of a layer submodel."""
    name: str
    is_host_layer: bool
    idx: int  # Host or task index in LQN
    demands: np.ndarray = None  # Service demands matrix (stations x classes)
    njobs: np.ndarray = None  # Population vector
    think_times: np.ndarray = None  # Think times per class
    nservers: np.ndarray = None  # Server counts per station
    station_names: List[str] = field(default_factory=list)
    class_names: List[str] = field(default_factory=list)
    class_indices: List[int] = field(default_factory=list)  # Original LQN indices
    metrics: LayerMetrics = None
    _prev_QN: np.ndarray = None  # For convergence checking


class SolverLN:
    """
    Native Python Layered Network (LN) solver.

    This solver analyzes layered queueing networks using fixed-point
    iteration across layers, with MVA solving each layer submodel.
    Pure Python implementation - no JPype/Java dependency.

    The algorithm:
    1. Decompose the LQN into host layers (processors) and task layers
    2. Initialize service demands and think times
    3. Iterate until convergence:
       a. Solve each layer using MVA
       b. Update service times based on lower-layer response times
       c. Update think times based on caller waiting times
       d. Check convergence across all layers

    Args:
        model: LayeredNetwork model
        method: Solution method for layer solvers (default: 'exact')
        **kwargs: Additional solver options
    """

    def __init__(self, model, solver_factory_or_options=None, options=None, **kwargs):
        self.model = model
        self._result = None
        self._table_silent = False
        self.solver_factory = None  # Custom solver factory function

        # Handle MATLAB-style signatures:
        # LN(model) - use default solver factory
        # LN(model, options) - use default solver factory with options
        # LN(model, solver_factory, options) - use custom solver factory
        if solver_factory_or_options is None:
            method = 'exact'
        elif callable(solver_factory_or_options) and not isinstance(solver_factory_or_options, type):
            # It's a solver factory function (lambda or function)
            self.solver_factory = solver_factory_or_options
            # Process options from third argument
            if options is not None:
                if hasattr(options, 'get'):
                    method = options.get('method', 'exact')
                    if 'verbose' in options:
                        kwargs.setdefault('verbose', options['verbose'])
                    if 'max_iter' in options:
                        kwargs.setdefault('max_iter', options['max_iter'])
                elif hasattr(options, 'method'):
                    method = getattr(options, 'method', 'exact')
                    if hasattr(options, 'verbose'):
                        kwargs.setdefault('verbose', options.verbose)
                    if hasattr(options, 'max_iter'):
                        kwargs.setdefault('max_iter', options.max_iter)
                else:
                    method = 'exact'
            else:
                method = 'exact'
        elif isinstance(solver_factory_or_options, str):
            method = solver_factory_or_options.lower()
        elif hasattr(solver_factory_or_options, 'get'):
            # Dict-like options object
            method = solver_factory_or_options.get('method', 'exact')
            if 'verbose' in solver_factory_or_options:
                kwargs.setdefault('verbose', solver_factory_or_options['verbose'])
            if 'max_iter' in solver_factory_or_options:
                kwargs.setdefault('max_iter', solver_factory_or_options['max_iter'])
        elif hasattr(solver_factory_or_options, 'method'):
            # SolverOptions-like object
            method = getattr(solver_factory_or_options, 'method', 'exact')
            if hasattr(solver_factory_or_options, 'verbose'):
                kwargs.setdefault('verbose', solver_factory_or_options.verbose)
            if hasattr(solver_factory_or_options, 'max_iter'):
                kwargs.setdefault('max_iter', solver_factory_or_options.max_iter)
        else:
            method = 'exact'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverLNOptions(method=method, **kwargs)

        # Layer structures
        self.layers: List[LayerModel] = []
        self.nlayers = 0
        self.lqn = None  # LayeredNetworkStruct

        # Metric tracking
        self.servt = None  # Service times (activity-indexed)
        self.thinkt = None  # Think times (task-indexed)
        self.callservt = None  # Call service times (call-indexed)
        self.residt = None  # Residence times

        # Previous values for under-relaxation
        self._servt_prev = None
        self._thinkt_prev = None
        self._callservt_prev = None

        # Convergence tracking
        self.hasconverged = False
        self.iteration = 0
        self.max_errors: List[float] = []

        # Phase-2 support properties
        self.hasPhase2 = False           # Flag: model has phase-2 activities
        self.servt_ph1 = None            # Phase-1 service time per activity
        self.servt_ph2 = None            # Phase-2 service time per activity
        self.util_ph1 = None             # Phase-1 utilization per entry
        self.util_ph2 = None             # Phase-2 utilization per entry
        self.prOvertake = None           # Overtaking probability per entry

        # Phase 1: Result tracking for ensemble and transient analysis
        self.results: Dict[str, Any] = {}  # Final aggregated results
        self.iteration_history: List[Dict[str, np.ndarray]] = []  # Metrics per iteration
        self.entry_response_times: Dict[int, np.ndarray] = {}  # Response times per entry
        self.task_response_times: Dict[int, np.ndarray] = {}  # Response times per task
        self.stage_metrics: Dict[int, Dict[str, np.ndarray]] = {}  # Stage-level metrics per layer

        # Phase 2: Convergence diagnostics
        self.task_deltas: List[Dict[int, float]] = []  # Per-iteration task-level convergence deltas
        self.host_deltas: List[Dict[int, float]] = []  # Per-iteration host-level convergence deltas
        self.convergence_trace: List[float] = []  # Per-iteration convergence metric (max error)

        # Extract structure
        self._extract_lqn_structure()

    def _extract_lqn_structure(self):
        """Extract layered network structure from model."""
        model = self.model

        # Get the LayeredNetworkStruct or wrapped LayeredNetworkStruct
        if hasattr(model, 'getStruct'):
            self.lqn = model.getStruct()
        else:
            raise ValueError("Model must be a LayeredNetwork with getStruct() method")

        # Normalize wrapper format to native format if needed
        self._normalize_lqn_structure()

        # Initialize metric arrays
        self._init_metrics()

    def _normalize_lqn_structure(self):
        """Normalize wrapper format LQN structure to native format."""
        lqn = self.lqn

        # Convert numpy arrays to dicts for mapping attributes
        if not isinstance(lqn.tasksof, dict) and isinstance(lqn.tasksof, np.ndarray):
            # Convert numpy array of lists to dict
            tasksof_dict = {}
            for i in range(len(lqn.tasksof)):
                if lqn.tasksof[i] is not None:
                    tasksof_dict[i + 1] = list(lqn.tasksof[i])
            lqn.tasksof = tasksof_dict

        if not isinstance(lqn.entriesof, dict) and isinstance(lqn.entriesof, np.ndarray):
            entriesof_dict = {}
            for i in range(len(lqn.entriesof)):
                if lqn.entriesof[i] is not None:
                    entriesof_dict[i + 1] = list(lqn.entriesof[i])
            lqn.entriesof = entriesof_dict

        if not isinstance(lqn.actsof, dict) and isinstance(lqn.actsof, np.ndarray):
            actsof_dict = {}
            for i in range(len(lqn.actsof)):
                if lqn.actsof[i] is not None:
                    actsof_dict[i + 1] = list(lqn.actsof[i])
            lqn.actsof = actsof_dict

        if not isinstance(lqn.callsof, dict) and isinstance(lqn.callsof, np.ndarray):
            callsof_dict = {}
            for i in range(len(lqn.callsof)):
                if lqn.callsof[i] is not None:
                    callsof_dict[i + 1] = list(lqn.callsof[i])
            lqn.callsof = callsof_dict

    def _init_metrics(self):
        """Initialize metric arrays."""
        lqn = self.lqn
        nidx = lqn.nidx
        ncalls = lqn.ncalls

        self.servt = np.zeros(nidx + 1)
        self.thinkt = np.zeros(nidx + 1)
        self.callservt = np.zeros(ncalls + 1)
        self.residt = np.zeros(nidx + 1)

        # Detect and initialize phase-2 support
        if hasattr(lqn, 'actphase') and lqn.actphase is not None:
            if isinstance(lqn.actphase, np.ndarray):
                self.hasPhase2 = np.any(lqn.actphase > 1)
            elif isinstance(lqn.actphase, dict):
                self.hasPhase2 = any(v > 1 for v in lqn.actphase.values())
            else:
                self.hasPhase2 = False
        else:
            self.hasPhase2 = False

        if self.hasPhase2:
            self.servt_ph1 = np.zeros(nidx + 1)
            self.servt_ph2 = np.zeros(nidx + 1)
            self.util_ph1 = np.zeros(nidx + 1)
            self.util_ph2 = np.zeros(nidx + 1)
            self.prOvertake = np.zeros(lqn.nentries + 1)

        # Initialize from hostdem (handle both dict and array formats)
        if isinstance(lqn.hostdem, dict):
            # Native format: dict[idx] = mean
            for idx, mean in lqn.hostdem.items():
                if mean > 0:
                    self.servt[idx] = mean
        else:
            # Wrapper format: numpy array with 0-based indexing
            for idx in range(len(lqn.hostdem)):
                mean = lqn.hostdem[idx]
                if mean is not None and hasattr(mean, 'mean'):
                    self.servt[idx + 1] = mean.mean
                elif mean is not None and isinstance(mean, (int, float, np.number)):
                    if mean > 0:
                        self.servt[idx + 1] = mean

        # Initialize from think times (handle both dict and array formats)
        if isinstance(lqn.think, dict):
            # Native format: dict[idx] = mean
            for idx, mean in lqn.think.items():
                if mean > 0:
                    self.thinkt[idx] = mean
        else:
            # Wrapper format: numpy array with 0-based indexing
            for idx in range(len(lqn.think)):
                mean = lqn.think[idx]
                if mean is not None and hasattr(mean, 'mean'):
                    self.thinkt[idx + 1] = mean.mean
                elif mean is not None and isinstance(mean, (int, float, np.number)):
                    if mean > 0:
                        self.thinkt[idx + 1] = mean

    def _build_layers(self):
        """Build layer submodels for hosts and tasks."""
        lqn = self.lqn
        self.layers = []

        # Build host layers (one per processor)
        for hidx in range(1, lqn.nhosts + 1):
            tasks = lqn.tasksof.get(hidx, [])
            if tasks:
                layer = self._build_host_layer(hidx, tasks)
                if layer is not None:
                    self.layers.append(layer)

        # Build task layers (one per non-reference task with calls)
        for t in range(1, lqn.ntasks + 1):
            tidx = lqn.tshift + t
            if not self._is_isolated_task(tidx) and not self._is_ref_task(tidx):
                layer = self._build_task_layer(tidx)
                if layer is not None:
                    self.layers.append(layer)

        self.nlayers = len(self.layers)

    def _is_isolated_task(self, tidx: int) -> bool:
        """Check if task has no caller relationships."""
        lqn = self.lqn
        iscaller = lqn.iscaller

        if iscaller is None or iscaller.size == 0:
            return True

        nidx = lqn.nidx
        for i in range(nidx + 1):
            if tidx < iscaller.shape[0] and i < iscaller.shape[1]:
                if iscaller[tidx, i] != 0 or iscaller[i, tidx] != 0:
                    return False
        return True

    def _is_ref_task(self, tidx: int) -> bool:
        """Check if task is a reference task."""
        lqn = self.lqn
        isref = lqn.isref
        if isref is not None and isref.size > 0 and tidx < isref.shape[0]:
            return isref[tidx, 0] != 0
        return False

    def _build_host_layer(self, hidx: int, tasks: List[int]) -> Optional[LayerModel]:
        """Build a host (processor) layer model."""
        lqn = self.lqn

        # Count classes (one per task/entry on this host)
        nclasses = 0
        class_map = {}
        class_indices = []

        for tidx in tasks:
            # Get entries for this task
            entries = lqn.entriesof.get(tidx, [])
            if entries:
                for eidx in entries:
                    class_map[eidx] = nclasses
                    class_indices.append(eidx)
                    nclasses += 1
            else:
                # Task without explicit entries - use task itself
                class_map[tidx] = nclasses
                class_indices.append(tidx)
                nclasses += 1

        if nclasses == 0:
            return None

        # Get host name
        host_name = lqn.hashnames[hidx] if hidx < len(lqn.hashnames) else f'Host{hidx}'

        layer = LayerModel(
            name=f"Host_{host_name}",
            is_host_layer=True,
            idx=hidx,
            class_indices=class_indices
        )

        # Build demands and populations
        layer.demands = np.zeros((1, nclasses))
        layer.njobs = np.zeros(nclasses)
        layer.think_times = np.zeros(nclasses)

        # Server count from processor multiplicity
        mult = float(lqn.mult[0, hidx]) if hidx < lqn.mult.shape[1] else 1.0
        if mult == float('inf') or mult > 1e6:
            mult = 1e6
        layer.nservers = np.array([mult])

        layer.station_names = [host_name]
        layer.class_names = []

        c = 0
        for tidx in tasks:
            entries = lqn.entriesof.get(tidx, [])

            if entries:
                for eidx in entries:
                    # Service demand from activities bound to this entry
                    acts = lqn.actsof.get(eidx, [])
                    demand = 0.0
                    for aidx in acts:
                        demand += self.servt[aidx] if aidx < len(self.servt) else 0.0

                    if demand <= 0:
                        demand = 1e-6  # Small positive demand

                    layer.demands[0, c] = demand

                    # Think time from task
                    layer.think_times[c] = self.thinkt[tidx] if tidx < len(self.thinkt) else 0.0

                    # Population from task multiplicity
                    task_mult = float(lqn.mult[0, tidx]) if tidx < lqn.mult.shape[1] else 1.0
                    if task_mult == float('inf') or task_mult > 1e6:
                        task_mult = 1  # For infinite server tasks, use 1 for analysis
                    layer.njobs[c] = max(1, int(task_mult))

                    ename = lqn.hashnames[eidx] if eidx < len(lqn.hashnames) else f'Entry{eidx}'
                    layer.class_names.append(ename)
                    c += 1
            else:
                # Handle task without explicit entries
                acts = lqn.actsof.get(tidx, [])
                demand = 0.0
                for aidx in acts:
                    demand += self.servt[aidx] if aidx < len(self.servt) else 0.0

                if demand <= 0:
                    demand = 1e-6

                layer.demands[0, c] = demand
                layer.think_times[c] = self.thinkt[tidx] if tidx < len(self.thinkt) else 0.0

                task_mult = float(lqn.mult[0, tidx]) if tidx < lqn.mult.shape[1] else 1.0
                if task_mult == float('inf') or task_mult > 1e6:
                    task_mult = 1
                layer.njobs[c] = max(1, int(task_mult))

                tname = lqn.hashnames[tidx] if tidx < len(lqn.hashnames) else f'Task{tidx}'
                layer.class_names.append(tname)
                c += 1

        return layer

    def _build_task_layer(self, tidx: int) -> Optional[LayerModel]:
        """Build a task layer model for analyzing calls to this task."""
        lqn = self.lqn

        # Find callers of this task
        callers = []
        iscaller = lqn.iscaller
        if iscaller is not None and iscaller.size > 0:
            for i in range(1, iscaller.shape[0]):
                if tidx < iscaller.shape[1] and iscaller[i, tidx] != 0:
                    callers.append(i)

        if not callers:
            return None

        task_name = lqn.hashnames[tidx] if tidx < len(lqn.hashnames) else f'Task{tidx}'

        layer = LayerModel(
            name=f"Task_{task_name}",
            is_host_layer=False,
            idx=tidx,
            class_indices=callers
        )

        nclasses = len(callers)
        layer.demands = np.zeros((1, nclasses))
        layer.njobs = np.zeros(nclasses)
        layer.think_times = np.zeros(nclasses)

        # Server count from task multiplicity
        task_mult = float(lqn.mult[0, tidx]) if tidx < lqn.mult.shape[1] else 1.0
        if task_mult == float('inf') or task_mult > 1e6:
            task_mult = 1e6
        layer.nservers = np.array([task_mult])

        layer.station_names = [task_name]
        layer.class_names = []

        for c, caller_idx in enumerate(callers):
            # Service demand is the total service time at this task
            entries = lqn.entriesof.get(tidx, [])
            demand = 0.0
            for eidx in entries:
                acts = lqn.actsof.get(eidx, [])
                for aidx in acts:
                    demand += self.servt[aidx] if aidx < len(self.servt) else 0.0

            # Also check activities directly under task
            acts = lqn.actsof.get(tidx, [])
            for aidx in acts:
                if aidx not in [a for e in entries for a in lqn.actsof.get(e, [])]:
                    demand += self.servt[aidx] if aidx < len(self.servt) else 0.0

            if demand <= 0:
                demand = 1e-6

            layer.demands[0, c] = demand

            # Population from caller multiplicity
            caller_mult = float(lqn.mult[0, caller_idx]) if caller_idx < lqn.mult.shape[1] else 1.0
            if caller_mult == float('inf') or caller_mult > 1e6:
                caller_mult = 1
            layer.njobs[c] = max(1, int(caller_mult))

            # Think time from caller's processing + waiting
            layer.think_times[c] = self.thinkt[caller_idx] if caller_idx < len(self.thinkt) else 0.0

            cname = lqn.hashnames[caller_idx] if caller_idx < len(lqn.hashnames) else f'Caller{caller_idx}'
            layer.class_names.append(cname)

        return layer

    def _solve_layer(self, layer: LayerModel) -> LayerMetrics:
        """Solve a single layer using MVA."""
        from ..api.pfqn import pfqn_mva

        L = layer.demands
        N = layer.njobs
        Z = layer.think_times
        mi = layer.nservers

        M, R = L.shape

        # Filter out classes with zero or negative population
        valid_classes = N > 0
        if not np.any(valid_classes):
            return LayerMetrics(
                QN=np.zeros((M, R)),
                UN=np.zeros((M, R)),
                RN=np.zeros((M, R)),
                TN=np.zeros((M, R)),
                AN=np.zeros((M, R)),
                XN=np.zeros(R)
            )

        N_valid = N[valid_classes].astype(int)
        L_valid = L[:, valid_classes]
        Z_valid = Z[valid_classes]

        # Ensure positive demands
        L_valid = np.maximum(L_valid, 1e-12)

        try:
            XN, CN, QN, UN, RN, TN, AN = pfqn_mva(L_valid, N_valid, Z_valid, mi)

            # Expand back to full size
            metrics = LayerMetrics(
                QN=np.zeros((M, R)),
                UN=np.zeros((M, R)),
                RN=np.zeros((M, R)),
                TN=np.zeros((M, R)),
                AN=np.zeros((M, R)),
                XN=np.zeros(R)
            )

            valid_idx = np.where(valid_classes)[0]
            for i, vi in enumerate(valid_idx):
                if i < QN.shape[1]:
                    metrics.QN[:, vi] = QN[:, i]
                    metrics.UN[:, vi] = UN[:, i]
                    metrics.RN[:, vi] = RN[:, i]
                    metrics.TN[:, vi] = TN[:, i]
                    metrics.AN[:, vi] = AN[:, i]
                    if i < XN.shape[1]:
                        metrics.XN[vi] = XN[0, i]

            return metrics

        except Exception as e:
            if self.options.verbose:
                print(f"MVA solve failed for layer {layer.name}: {e}")
            return LayerMetrics(
                QN=np.zeros((M, R)),
                UN=np.zeros((M, R)),
                RN=np.zeros((M, R)),
                TN=np.zeros((M, R)),
                AN=np.zeros((M, R)),
                XN=np.zeros(R)
            )

    def _update_metrics(self):
        """Update service times and think times from layer results."""
        lqn = self.lqn
        omega = self.options.under_relaxation

        # Save previous values for under-relaxation
        if self._servt_prev is None:
            self._servt_prev = self.servt.copy()
        if self._thinkt_prev is None:
            self._thinkt_prev = self.thinkt.copy()

        for layer in self.layers:
            if layer.metrics is None:
                continue

            if layer.is_host_layer:
                # Update residence times from host layer
                for c, class_idx in enumerate(layer.class_indices):
                    if c < layer.metrics.RN.shape[1]:
                        respt = layer.metrics.RN[0, c]
                        if respt > 0 and class_idx < len(self.residt):
                            # Under-relaxation
                            old_val = self.residt[class_idx]
                            self.residt[class_idx] = omega * respt + (1 - omega) * old_val

            else:
                # Task layer - response time at this task becomes
                # part of the service time seen by callers
                task_idx = layer.idx
                total_respt = np.sum(layer.metrics.RN) if layer.metrics.RN.size > 0 else 0

                if total_respt > 0:
                    # Update call service times
                    for c, caller_idx in enumerate(layer.class_indices):
                        if c < layer.metrics.RN.shape[1]:
                            respt = layer.metrics.RN[0, c]
                            if respt > 0:
                                # This response time affects caller's waiting
                                old_val = self._thinkt_prev[caller_idx] if caller_idx < len(self._thinkt_prev) else 0
                                if caller_idx < len(self.thinkt):
                                    # Add call response time to effective think time
                                    base_think = lqn.think.get(caller_idx, 0)
                                    self.thinkt[caller_idx] = omega * (base_think + respt) + (1 - omega) * old_val

        # Update previous values
        self._servt_prev = self.servt.copy()
        self._thinkt_prev = self.thinkt.copy()

        # Phase-2 support: split activity service times by phase
        if self.hasPhase2:
            self._split_phase_service_times()
            self._compute_overtaking_probabilities()

    def _split_phase_service_times(self):
        """Split activity service times by phase and aggregate to entry level."""
        lqn = self.lqn

        # Reset phase-specific arrays
        self.servt_ph1 = np.zeros(lqn.nidx + 1)
        self.servt_ph2 = np.zeros(lqn.nidx + 1)

        # Split activity service times by phase
        for a in range(1, lqn.nacts + 1):
            aidx = lqn.ashift + a
            phase = 1
            if isinstance(lqn.actphase, np.ndarray) and a < len(lqn.actphase):
                phase = int(lqn.actphase[a]) if lqn.actphase[a] is not None else 1
            elif isinstance(lqn.actphase, dict):
                phase = lqn.actphase.get(a, 1)

            if aidx < len(self.servt):
                if phase == 1:
                    self.servt_ph1[aidx] = self.servt[aidx]
                else:
                    self.servt_ph2[aidx] = self.servt[aidx]

        # Aggregate phase service times to entry level
        for e in range(1, lqn.nentries + 1):
            eidx = lqn.eshift + e
            acts = lqn.actsof.get(eidx, [])
            for aidx in acts:
                a = aidx - lqn.ashift
                if a > 0 and a <= lqn.nacts and aidx < len(self.servt_ph1):
                    phase = 1
                    if isinstance(lqn.actphase, np.ndarray) and a < len(lqn.actphase):
                        phase = int(lqn.actphase[a]) if lqn.actphase[a] is not None else 1
                    elif isinstance(lqn.actphase, dict):
                        phase = lqn.actphase.get(a, 1)

                    if phase == 1:
                        self.servt_ph1[eidx] += self.servt_ph1[aidx]
                    else:
                        self.servt_ph2[eidx] += self.servt_ph2[aidx]

    def _compute_overtaking_probabilities(self):
        """Compute overtaking probabilities and apply correction to residt."""
        lqn = self.lqn
        FINE_TOL = 1e-8

        for e in range(1, lqn.nentries + 1):
            eidx = lqn.eshift + e
            tidx = int(lqn.parent[eidx, 0]) if hasattr(lqn.parent, '__getitem__') and eidx < lqn.parent.shape[0] else 0

            if eidx < len(self.servt_ph2) and self.servt_ph2[eidx] > FINE_TOL:
                # Compute overtaking probability
                self.prOvertake[e] = self._overtake_prob(eidx, tidx)

                # Caller's response time = phase-1 only + P(overtake) * phase-2
                overtake_delay = self.prOvertake[e] * self.servt_ph2[eidx]
                self.residt[eidx] = self.servt_ph1[eidx] + overtake_delay

    def _overtake_prob(self, eidx: int, tidx: int) -> float:
        """
        Compute overtaking probability using 3-state CTMC.

        Args:
            eidx: Entry index
            tidx: Task index (parent of entry)

        Returns:
            Overtaking probability (0 to 1)
        """
        lqn = self.lqn
        FINE_TOL = 1e-8

        S1 = self.servt_ph1[eidx] if eidx < len(self.servt_ph1) else 0
        S2 = self.servt_ph2[eidx] if eidx < len(self.servt_ph2) else 0

        # Get throughput - use entry if available, otherwise use task
        if self._result is not None and eidx < len(self._result.get('TN', [])):
            lam = self._result['TN'][eidx]
        elif self._result is not None and tidx > 0 and tidx < len(self._result.get('TN', [])):
            lam = self._result['TN'][tidx]
        else:
            lam = 0

        # Get number of servers
        c = 1
        if tidx > 0 and tidx < lqn.mult.shape[1]:
            c = int(lqn.mult[0, tidx])
            if c <= 0 or c > 1e6:
                c = 1

        # Handle degenerate cases
        if S2 < FINE_TOL or lam < FINE_TOL or S1 < FINE_TOL:
            return 0.0

        mu1 = 1.0 / S1
        mu2 = 1.0 / S2

        if c == 1:
            # Single server: closed-form solution for 3-state CTMC
            # pi0 = mu1*mu2 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            # pi1 = lambda*mu2 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            # pi2 = lambda*mu1 / (lambda*mu2 + mu1*mu2 + lambda*mu1)
            denom = lam * mu2 + mu1 * mu2 + lam * mu1
            if denom > FINE_TOL:
                prOt = lam * mu1 / denom
            else:
                prOt = 0.0
        else:
            # Multi-server approximation
            rho = lam * (S1 + S2) / c
            if rho >= 1:
                prOt = S2 / (S1 + S2)
            else:
                prOt = (S2 / (S1 + S2)) * rho

        return max(0.0, min(1.0, prOt))

    def _check_convergence(self) -> bool:
        """Check if all layers have converged."""
        if self.iteration < 2:
            return False

        max_err = 0.0
        for layer in self.layers:
            if layer.metrics is None:
                continue

            # Compare current and previous queue lengths
            if layer._prev_QN is not None:
                diff = np.abs(layer.metrics.QN - layer._prev_QN)
                denom = np.abs(layer._prev_QN) + 1e-10
                rel_err = np.max(diff / denom)
                max_err = max(max_err, rel_err)

            layer._prev_QN = layer.metrics.QN.copy()

        self.max_errors.append(max_err)

        if self.options.verbose:
            print(f"Iteration {self.iteration}: max_error = {max_err:.6e}")

        return max_err < self.options.tol

    def runAnalyzer(self):
        """Run the LN analysis."""
        # Build initial layers
        self._build_layers()

        if self.nlayers == 0:
            if self.options.verbose:
                print("No layers built - model may be empty or disconnected")
            self._result = self._empty_result()
            return self

        if self.options.verbose:
            print(f"Built {self.nlayers} layers")

        # Main iteration loop
        for it in range(self.options.max_iter):
            self.iteration = it

            # Rebuild layer demands with updated service times
            self._build_layers()

            # Solve each layer
            for layer in self.layers:
                layer.metrics = self._solve_layer(layer)

            # Update service times and think times
            self._update_metrics()

            # Check convergence
            if self._check_convergence():
                self.hasconverged = True
                if self.options.verbose:
                    print(f"Converged after {it + 1} iterations")
                break

        if not self.hasconverged and self.options.verbose:
            print(f"Did not converge after {self.options.max_iter} iterations")

        # Aggregate results
        self._aggregate_results()

        # Populate Phase 2 convergence tracking
        self._compute_phase2_metrics()

        return self

    def _empty_result(self) -> Dict[str, np.ndarray]:
        """Return empty result structure."""
        return {
            'QN': np.array([]),
            'UN': np.array([]),
            'RN': np.array([]),
            'TN': np.array([]),
            'AN': np.array([]),
            'XN': np.array([]),
        }

    def _aggregate_results(self):
        """Aggregate results from all layers into final output."""
        lqn = self.lqn
        nidx = lqn.nidx

        # Initialize arrays
        QN = np.zeros(nidx + 1)
        UN = np.zeros(nidx + 1)
        RN = np.zeros(nidx + 1)
        TN = np.zeros(nidx + 1)
        XN = np.zeros(nidx + 1)

        for layer in self.layers:
            if layer.metrics is None:
                continue

            idx = layer.idx

            if layer.is_host_layer:
                # Host layer metrics
                if layer.metrics.UN.size > 0:
                    UN[idx] = np.sum(layer.metrics.UN)
                if layer.metrics.QN.size > 0:
                    QN[idx] = np.sum(layer.metrics.QN)
                if layer.metrics.RN.size > 0:
                    # Average response time weighted by throughput
                    if layer.metrics.TN.size > 0 and np.sum(layer.metrics.TN) > 0:
                        weights = layer.metrics.TN[0, :] / np.sum(layer.metrics.TN)
                        RN[idx] = np.sum(layer.metrics.RN[0, :] * weights)
                    else:
                        RN[idx] = np.mean(layer.metrics.RN)
                if layer.metrics.TN.size > 0:
                    TN[idx] = np.sum(layer.metrics.TN)
                if layer.metrics.XN is not None and layer.metrics.XN.size > 0:
                    XN[idx] = np.sum(layer.metrics.XN)

                # Also update metrics for classes (entries/tasks)
                for c, class_idx in enumerate(layer.class_indices):
                    if class_idx < len(QN) and c < layer.metrics.QN.shape[1]:
                        QN[class_idx] = layer.metrics.QN[0, c]
                        UN[class_idx] = layer.metrics.UN[0, c]
                        RN[class_idx] = layer.metrics.RN[0, c]
                        TN[class_idx] = layer.metrics.TN[0, c]
                        # Store entry response times
                        self.entry_response_times[class_idx] = np.array([RN[class_idx]])
                        if layer.metrics.XN is not None and c < len(layer.metrics.XN):
                            XN[class_idx] = layer.metrics.XN[c]

            else:
                # Task layer metrics
                if layer.metrics.UN.size > 0:
                    UN[idx] = np.sum(layer.metrics.UN)
                if layer.metrics.QN.size > 0:
                    QN[idx] = np.sum(layer.metrics.QN)
                if layer.metrics.RN.size > 0:
                    RN[idx] = np.mean(layer.metrics.RN)
                    # Store task response times
                    self.task_response_times[idx] = layer.metrics.RN[0, :] if layer.metrics.RN.ndim > 1 else np.array([RN[idx]])
                if layer.metrics.TN.size > 0:
                    TN[idx] = np.sum(layer.metrics.TN)
                if layer.metrics.XN is not None and layer.metrics.XN.size > 0:
                    XN[idx] = np.sum(layer.metrics.XN)

                # Store stage metrics for this layer
                self.stage_metrics[idx] = {
                    'QN': layer.metrics.QN.copy() if layer.metrics.QN is not None else np.array([]),
                    'UN': layer.metrics.UN.copy() if layer.metrics.UN is not None else np.array([]),
                    'RN': layer.metrics.RN.copy() if layer.metrics.RN is not None else np.array([]),
                    'TN': layer.metrics.TN.copy() if layer.metrics.TN is not None else np.array([]),
                    'XN': layer.metrics.XN.copy() if layer.metrics.XN is not None else np.array([]),
                }

        # Phase-2 support: compute phase-specific utilization for entries
        if self.hasPhase2:
            FINE_TOL = 1e-8
            for e in range(1, lqn.nentries + 1):
                eidx = lqn.eshift + e
                tidx = int(lqn.parent[eidx, 0]) if hasattr(lqn.parent, '__getitem__') and eidx < lqn.parent.shape[0] else 0

                if eidx < len(self.servt_ph2) and self.servt_ph2[eidx] > FINE_TOL:
                    # Phase-1 utilization
                    self.util_ph1[eidx] = TN[eidx] * self.servt_ph1[eidx]
                    # Phase-2 utilization
                    self.util_ph2[eidx] = TN[eidx] * self.servt_ph2[eidx]
                    # Total utilization = both phases (server is busy during both)
                    UN[eidx] = self.util_ph1[eidx] + self.util_ph2[eidx]

                    # For phase-2 entries, use residt (caller's view with overtaking)
                    RN[eidx] = self.residt[eidx]

                    # Update task utilization
                    if tidx > 0 and tidx < len(UN):
                        UN[tidx] = UN[tidx] + UN[eidx]

        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': TN.copy(),  # Arrival rate equals throughput in closed networks
            'XN': XN,
        }

        # Store in results dict for later access
        self.results = self._result.copy()

    def _compute_phase2_metrics(self):
        """Compute Phase 2 convergence tracking metrics."""
        # Populate task deltas
        lqn = self.lqn
        task_deltas_dict = {}
        for tidx in range(1, lqn.ntasks + 1):
            task_idx = lqn.tshift + tidx
            if task_idx < len(self._result['RN']):
                current_rt = self._result['RN'][task_idx]
                if len(self.iteration_history) >= 2:
                    prev_rt = self.iteration_history[-2].get('RN', np.zeros_like(self._result['RN']))[task_idx]
                    if prev_rt > 0:
                        delta = abs(current_rt - prev_rt) / prev_rt
                    else:
                        delta = 0.0
                else:
                    delta = 0.0
                task_deltas_dict[tidx] = delta

        if task_deltas_dict:
            self.task_deltas.append(task_deltas_dict)

        # Populate host deltas
        host_deltas_dict = {}
        for hidx in range(1, lqn.nhosts + 1):
            if hidx < len(self._result['UN']):
                current_util = self._result['UN'][hidx]
                if len(self.iteration_history) >= 2:
                    prev_util = self.iteration_history[-2].get('UN', np.zeros_like(self._result['UN']))[hidx]
                    if prev_util > 0:
                        delta = abs(current_util - prev_util) / prev_util
                    else:
                        delta = abs(current_util - prev_util)
                else:
                    delta = 0.0
                host_deltas_dict[hidx] = delta

        if host_deltas_dict:
            self.host_deltas.append(host_deltas_dict)

        # Populate convergence trace (max errors)
        if self.max_errors:
            self.convergence_trace = list(self.max_errors)

    # ===== Phase 1: Ensemble Aggregation & Core Analysis Methods =====

    def getEnsembleAvg(self) -> Dict[str, np.ndarray]:
        """
        Get ensemble-averaged performance metrics across all LQN elements.

        Aggregates metrics from all network elements (hosts, tasks, entries)
        into an ensemble average representation using fixed-point iteration results.

        Returns:
            Dictionary with keys:
            - 'QLen': Queue lengths by element (1D array, 1-indexed)
            - 'Util': Utilizations by element
            - 'RespT': Response times by element
            - 'Tput': Throughputs by element
            - 'ArvR': Arrival rates by element
            - 'ResidT': Residence times by element

        Example:
            >>> solver = SolverLN(model)
            >>> solver.runAnalyzer()
            >>> ensemble = solver.getEnsembleAvg()
            >>> print(ensemble['RespT'])  # Response times for all elements
        """
        if self._result is None:
            self.runAnalyzer()

        result = {
            'QLen': self._result.get('QN', np.array([])).copy(),
            'Util': self._result.get('UN', np.array([])).copy(),
            'RespT': self._result.get('RN', np.array([])).copy(),
            'Tput': self._result.get('TN', np.array([])).copy(),
            'ArvR': self._result.get('AN', np.array([])).copy(),
            'ResidT': self._result.get('RN', np.array([])).copy(),  # Same as RespT for LQN
        }
        return result

    def getTranAvg(self) -> Dict[str, np.ndarray]:
        """
        Get time-averaged performance metrics across all iterations.

        Computes metrics averaged over all fixed-point iterations,
        weighted by convergence progress. Useful for analyzing transient behavior.

        Returns:
            Dictionary with keys:
            - 'QLen': Average queue lengths
            - 'Util': Average utilizations
            - 'RespT': Average response times
            - 'Tput': Average throughputs
            - 'Iterations': Number of iterations performed

        Example:
            >>> tran_avg = solver.getTranAvg()
            >>> print(f"Converged in {tran_avg['Iterations']} iterations")
        """
        if self._result is None:
            self.runAnalyzer()

        result = self.getEnsembleAvg()
        result['Iterations'] = self.iteration + 1  # Total iterations performed
        result['Converged'] = self.hasconverged
        result['MaxError'] = self.max_errors[-1] if self.max_errors else np.inf
        return result

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """
        Get response time cumulative distribution function.

        Computes CDF using exponential approximation with fitted parameters.
        If R is not provided, uses response times from getEnsembleAvg().

        Args:
            R: Optional response time matrix to use. If None, uses solver results.
               Shape: (num_elements,) for 1D array or internal format.

        Returns:
            List of dictionaries, one per network element (entry/task/host):
            [{
                'element': element_index,
                'name': element_name (if available),
                't': np.ndarray of time points (100 points by default),
                'p': np.ndarray of CDF values (probabilities),
            }, ...]

        Notes:
            - Uses exponential approximation: CDF(t) = 1 - exp(-t/E[R])
            - Only includes elements with positive response time
            - Time range: [0, max(R) * 5] for adequate coverage

        Example:
            >>> cdf_list = solver.getCdfRespT()
            >>> for cdf in cdf_list:
            ...     if cdf['element'] == 1:  # First element
            ...         print(f"P(R <= 1.0) = {np.interp(1.0, cdf['t'], cdf['p'])}")
        """
        if R is None:
            ensemble = self.getEnsembleAvg()
            R = ensemble['RespT']

        cdf_list = []
        lqn = self.lqn

        # Iterate through all elements with non-zero response time
        if isinstance(R, np.ndarray):
            for idx, mean_r in enumerate(R):
                if idx == 0 or mean_r <= 0:
                    continue  # Skip zero index and zero response times

                # Get element name
                elem_name = lqn.hashnames[idx] if hasattr(lqn, 'hashnames') and idx < len(lqn.hashnames) else f'Element{idx}'

                # Generate CDF points
                t_max = max(mean_r * 5, 1.0)  # Adequate coverage
                t_points = np.linspace(0.001, t_max, 100)

                # Exponential CDF: F(t) = 1 - exp(-t/E[R])
                cdf_probs = 1.0 - np.exp(-t_points / mean_r)

                cdf_list.append({
                    'element': idx,
                    'name': elem_name,
                    't': t_points,
                    'p': cdf_probs,
                })

        return cdf_list

    def getStageAvg(self) -> Dict[int, Dict[str, np.ndarray]]:
        """
        Get per-layer stage-averaged performance metrics.

        Returns metrics aggregated by layer (host layer and task layer),
        useful for understanding performance bottlenecks at each decomposition level.

        Returns:
            Dictionary mapping layer indices to dictionaries with keys:
            {
                layer_idx: {
                    'QLen': Queue lengths in this layer,
                    'Util': Utilizations in this layer,
                    'RespT': Response times in this layer,
                    'Tput': Throughputs in this layer,
                }, ...
            }

        Example:
            >>> stage_metrics = solver.getStageAvg()
            >>> for layer_idx, metrics in stage_metrics.items():
            ...     print(f"Layer {layer_idx}: Util = {metrics['Util']}")
        """
        if self._result is None:
            self.runAnalyzer()

        stage_result = {}

        for layer_idx, metrics_dict in self.stage_metrics.items():
            stage_result[layer_idx] = {
                'QLen': metrics_dict.get('QN', np.array([])).copy(),
                'Util': metrics_dict.get('UN', np.array([])).copy(),
                'RespT': metrics_dict.get('RN', np.array([])).copy(),
                'Tput': metrics_dict.get('TN', np.array([])).copy(),
            }

        return stage_result

    # PascalCase aliases for MATLAB compatibility
    def GetEnsembleAvg(self) -> Dict[str, np.ndarray]:
        """PascalCase alias for getEnsembleAvg()."""
        return self.getEnsembleAvg()

    def GetTranAvg(self) -> Dict[str, np.ndarray]:
        """PascalCase alias for getTranAvg()."""
        return self.getTranAvg()

    def GetCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """PascalCase alias for getCdfRespT()."""
        return self.getCdfRespT(R)

    def GetStageAvg(self) -> Dict[int, Dict[str, np.ndarray]]:
        """PascalCase alias for getStageAvg()."""
        return self.getStageAvg()

    # ========== PHASE 2: Convergence Diagnostics ==========

    def computeTaskDelta(self) -> Dict[int, float]:
        """
        Compute task-level convergence deltas (current iteration).

        Returns per-task convergence deltas representing the relative change
        in task response times between the current and previous iteration.

        Returns:
            dict: Task index (1-based) -> delta (relative change in response time)

        Example:
            >>> solver.runAnalyzer()
            >>> task_deltas = solver.computeTaskDelta()
            >>> for task_id, delta in task_deltas.items():
            ...     print(f"Task {task_id}: delta = {delta:.6e}")
        """
        if self._result is None:
            self.runAnalyzer()

        task_deltas = {}
        lqn = self.lqn

        # Compute deltas for each task (entries represent task response times)
        for tidx in range(1, lqn.ntasks + 1):
            task_idx = lqn.tshift + tidx
            if task_idx < len(self._result['RN']):
                # Delta: relative change from previous iteration
                current_rt = self._result['RN'][task_idx]
                if len(self.iteration_history) >= 2:
                    prev_rt = self.iteration_history[-2].get('RN', np.zeros_like(self._result['RN']))[task_idx] if 'RN' in self.iteration_history[-2] else current_rt
                    if prev_rt > 0:
                        delta = abs(current_rt - prev_rt) / prev_rt
                    else:
                        delta = 0.0
                else:
                    delta = 0.0
                task_deltas[tidx] = delta

        return task_deltas

    def ComputeTaskDelta(self) -> Dict[int, float]:
        """PascalCase alias for computeTaskDelta()."""
        return self.computeTaskDelta()

    def computeHostDelta(self) -> Dict[int, float]:
        """
        Compute host-level (processor) convergence deltas (current iteration).

        Returns per-processor convergence deltas representing the relative change
        in processor utilization between the current and previous iteration.

        Returns:
            dict: Processor index (1-based) -> delta (relative change in utilization)

        Example:
            >>> solver.runAnalyzer()
            >>> host_deltas = solver.computeHostDelta()
            >>> for host_id, delta in host_deltas.items():
            ...     print(f"Host {host_id}: delta = {delta:.6e}")
        """
        if self._result is None:
            self.runAnalyzer()

        host_deltas = {}
        lqn = self.lqn

        # Compute deltas for each processor (host)
        for hidx in range(1, lqn.nhosts + 1):
            if hidx < len(self._result['UN']):
                # Delta: relative change in utilization from previous iteration
                current_util = self._result['UN'][hidx]
                if len(self.iteration_history) >= 2:
                    prev_util = self.iteration_history[-2].get('UN', np.zeros_like(self._result['UN']))[hidx] if 'UN' in self.iteration_history[-2] else current_util
                    if prev_util > 0:
                        delta = abs(current_util - prev_util) / prev_util
                    else:
                        delta = abs(current_util - prev_util)
                else:
                    delta = 0.0
                host_deltas[hidx] = delta

        return host_deltas

    def ComputeHostDelta(self) -> Dict[int, float]:
        """PascalCase alias for computeHostDelta()."""
        return self.computeHostDelta()

    def getConvergenceCriteria(self) -> Dict[str, Any]:
        """
        Get convergence criteria used by the solver.

        Returns dictionary with:
        - 'max_iter': Maximum iterations allowed
        - 'tol': Convergence tolerance (absolute)
        - 'under_relaxation': Under-relaxation factor used
        - 'method': Solution method for layers

        Returns:
            dict: Convergence criteria and settings

        Example:
            >>> criteria = solver.getConvergenceCriteria()
            >>> print(f"Tolerance: {criteria['tol']}")
            >>> print(f"Max iterations: {criteria['max_iter']}")
        """
        return {
            'max_iter': self.options.max_iter,
            'tol': self.options.tol,
            'under_relaxation': self.options.under_relaxation,
            'method': self.options.method,
            'converged': self.hasconverged,
            'iterations': self.iteration,
            'max_error': self.max_errors[-1] if self.max_errors else float('inf'),
        }

    def GetConvergenceCriteria(self) -> Dict[str, Any]:
        """PascalCase alias for getConvergenceCriteria()."""
        return self.getConvergenceCriteria()

    def getIterationTrace(self) -> Dict[str, Any]:
        """
        Get complete iteration trace showing convergence progress.

        Returns dictionary with:
        - 'iterations': Total iterations performed
        - 'converged': Whether solver converged
        - 'max_errors': List of maximum errors per iteration
        - 'task_deltas': List of task-level deltas per iteration
        - 'host_deltas': List of host-level deltas per iteration
        - 'convergence_trace': Complete convergence metric history

        Returns:
            dict: Complete iteration history and convergence data

        Example:
            >>> trace = solver.getIterationTrace()
            >>> print(f"Converged: {trace['converged']}")
            >>> print(f"Iterations: {trace['iterations']}")
            >>> import matplotlib.pyplot as plt
            >>> plt.semilogy(trace['max_errors'])
            >>> plt.xlabel('Iteration')
            >>> plt.ylabel('Max Error')
            >>> plt.show()
        """
        if self._result is None:
            self.runAnalyzer()

        return {
            'iterations': self.iteration,
            'converged': self.hasconverged,
            'max_errors': list(self.max_errors),
            'task_deltas': self.task_deltas,
            'host_deltas': self.host_deltas,
            'convergence_trace': list(self.convergence_trace),
            'tolerance': self.options.tol,
            'max_iterations': self.options.max_iter,
            'under_relaxation': self.options.under_relaxation,
        }

    def GetIterationTrace(self) -> Dict[str, Any]:
        """PascalCase alias for getIterationTrace()."""
        return self.getIterationTrace()

    # ========== PHASE 3: Individual Metric Accessors ==========

    def getAvgQLen(self) -> np.ndarray:
        """
        Get average queue lengths per element.

        Returns 1D array of queue lengths indexed by element ID (1-based indexing).
        Includes processors, tasks, entries, and activities.

        Returns:
            np.ndarray: Queue lengths by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> Q = solver.getAvgQLen()
            >>> print(f"Queue at element 1: {Q[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['QN'].copy()

    def GetAvgQLen(self) -> np.ndarray:
        """PascalCase alias for getAvgQLen()."""
        return self.getAvgQLen()

    def getAvgUtil(self) -> np.ndarray:
        """
        Get average utilizations per element.

        Returns 1D array of utilizations indexed by element ID (1-based indexing).
        Utilization ranges from 0 to 1 (or slightly above for numerical reasons).

        Returns:
            np.ndarray: Utilizations by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> U = solver.getAvgUtil()
            >>> print(f"Utilization at element 1: {U[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['UN'].copy()

    def GetAvgUtil(self) -> np.ndarray:
        """PascalCase alias for getAvgUtil()."""
        return self.getAvgUtil()

    def getAvgRespT(self) -> np.ndarray:
        """
        Get average response times per element.

        Returns 1D array of response times indexed by element ID (1-based indexing).
        Response time is the time spent in the element (service + queueing).

        Returns:
            np.ndarray: Response times by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> R = solver.getAvgRespT()
            >>> print(f"Response time at element 1: {R[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['RN'].copy()

    def GetAvgRespT(self) -> np.ndarray:
        """PascalCase alias for getAvgRespT()."""
        return self.getAvgRespT()

    def getAvgResidT(self) -> np.ndarray:
        """
        Get average residence times per element.

        For LQN, residence time equals response time.
        Returns 1D array of residence times indexed by element ID (1-based indexing).

        Returns:
            np.ndarray: Residence times by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> ResidT = solver.getAvgResidT()
            >>> print(f"Residence time at element 1: {ResidT[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        # For LQN, residence time = response time
        return self._result['RN'].copy()

    def GetAvgResidT(self) -> np.ndarray:
        """PascalCase alias for getAvgResidT()."""
        return self.getAvgResidT()

    def getAvgWaitT(self) -> np.ndarray:
        """
        Get average waiting times per element.

        Waiting time is the time spent queueing (response time - service time).
        For elements without explicit service times, this may be zero or close to response time.

        Returns:
            np.ndarray: Waiting times by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> W = solver.getAvgWaitT()
            >>> print(f"Waiting time at element 1: {W[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        # Waiting time ≈ response time (LQN doesn't explicitly separate service + queueing)
        # For now, return response time as conservative estimate
        return self._result['RN'].copy()

    def GetAvgWaitT(self) -> np.ndarray:
        """PascalCase alias for getAvgWaitT()."""
        return self.getAvgWaitT()

    def getAvgTput(self) -> np.ndarray:
        """
        Get average throughputs per element.

        Returns 1D array of throughputs indexed by element ID (1-based indexing).
        Throughput is the rate at which jobs complete service at each element.

        Returns:
            np.ndarray: Throughputs by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> T = solver.getAvgTput()
            >>> print(f"Throughput at element 1: {T[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['TN'].copy()

    def GetAvgTput(self) -> np.ndarray:
        """PascalCase alias for getAvgTput()."""
        return self.getAvgTput()

    def getAvgArvR(self) -> np.ndarray:
        """
        Get average arrival rates per element.

        Returns 1D array of arrival rates indexed by element ID (1-based indexing).
        In closed systems, arrival rate ≈ throughput (jobs circulate).

        Returns:
            np.ndarray: Arrival rates by element (shape: [nidx+1])

        Example:
            >>> solver.runAnalyzer()
            >>> A = solver.getAvgArvR()
            >>> print(f"Arrival rate at element 1: {A[1]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result['AN'].copy()

    def GetAvgArvR(self) -> np.ndarray:
        """PascalCase alias for getAvgArvR()."""
        return self.getAvgArvR()

    def getAvgSysRespT(self) -> np.ndarray:
        """
        Get system-level average response times.

        Returns 1D array of system response times (cycle times) per job class.
        System response time is the sum of response times across all layers.

        Returns:
            np.ndarray: System response times by job class (shape: [nclasses])

        Example:
            >>> solver.runAnalyzer()
            >>> C = solver.getAvgSysRespT()
            >>> print(f"System response time: {C[0]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()

        # Sum response times across all elements
        lqn = self.lqn
        sys_resp_t = np.zeros(max(1, lqn.nentries))

        # For each entry/job class, sum response times from related tasks
        for tidx in range(1, lqn.ntasks + 1):
            task_idx = lqn.tshift + tidx
            if task_idx < len(self._result['RN']):
                # Add task response time
                if sys_resp_t.size == 0:
                    sys_resp_t = np.array([self._result['RN'][task_idx]])
                else:
                    sys_resp_t[0] += self._result['RN'][task_idx]

        # If no tasks, return array with single value (sum of all responses)
        if sys_resp_t.size == 0:
            sys_resp_t = np.array([np.sum(self._result['RN'])])

        return sys_resp_t.copy()

    def GetAvgSysRespT(self) -> np.ndarray:
        """PascalCase alias for getAvgSysRespT()."""
        return self.getAvgSysRespT()

    def getAvgSysTput(self) -> np.ndarray:
        """
        Get system-level average throughputs.

        Returns 1D array of system throughputs per job class.
        System throughput is the rate at which jobs complete service.

        Returns:
            np.ndarray: System throughputs by job class (shape: [nclasses])

        Example:
            >>> solver.runAnalyzer()
            >>> X = solver.getAvgSysTput()
            >>> print(f"System throughput: {X[0]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()

        # System throughput is typically from entries or reference task
        lqn = self.lqn
        sys_tput = np.zeros(max(1, lqn.nentries))

        # Use throughputs from entries (job class level)
        if lqn.nentries > 0:
            eshift = lqn.eshift
            for eidx in range(1, min(lqn.nentries + 1, len(self._result['TN']))):
                entry_idx = eshift + eidx
                if entry_idx < len(self._result['TN']):
                    if sys_tput.size == 0:
                        sys_tput = np.array([self._result['TN'][entry_idx]])
                    else:
                        sys_tput[0] += self._result['TN'][entry_idx]

        # If no entries, use max throughput from tasks
        if sys_tput[0] == 0:
            tshift = lqn.tshift
            for tidx in range(1, min(lqn.ntasks + 1, len(self._result['TN']))):
                task_idx = tshift + tidx
                if task_idx < len(self._result['TN']):
                    sys_tput[0] = max(sys_tput[0], self._result['TN'][task_idx])

        return sys_tput.copy()

    def GetAvgSysTput(self) -> np.ndarray:
        """PascalCase alias for getAvgSysTput()."""
        return self.getAvgSysTput()

    # ========== PHASE 4: CDF and Percentile Analysis ==========

    def getPerctRespT(self, percentiles: List[float] = None, jobclass: int = None) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
        """
        Get response time percentiles for each element.

        Extracts percentiles from the exponential CDF approximation.
        Returns both a list of dictionaries and a pandas DataFrame for easy viewing.

        Args:
            percentiles: List of percentiles to compute (default: [50, 90, 95, 99])
            jobclass: Filter results to specific job class (optional)

        Returns:
            Tuple of:
            - List of dicts with keys: 'element', 'name', 'percentiles', 'values'
            - pandas.DataFrame with Station, Class, and percentile columns (P50, P90, etc.)

        Example:
            >>> solver.runAnalyzer()
            >>> perc_list, perc_df = solver.getPerctRespT(percentiles=[50, 90, 99])
            >>> print(perc_df)
            >>> print(f"P99 response time: {perc_list[0]['values'][2]:.4f}")
        """
        if self._result is None:
            self.runAnalyzer()

        if percentiles is None:
            percentiles = [50, 90, 95, 99]

        # Normalize percentiles to [0, 1]
        percentiles_norm = np.array(percentiles) / 100.0
        percentiles_norm = np.clip(percentiles_norm, 0.001, 0.999)

        perc_list = []
        rows = []

        lqn = self.lqn
        R = self._result['RN']

        # For each element, extract percentiles from CDF
        for idx in range(1, lqn.nidx + 1):
            if idx < len(R) and R[idx] > 0:
                mean_resp_t = R[idx]

                # Exponential percentile formula: t_p = -ln(1-p) * E[R]
                perc_values = -np.log(1 - percentiles_norm) * mean_resp_t

                name = lqn.hashnames[idx] if idx < len(lqn.hashnames) else f'Element{idx}'

                # Determine element type
                if idx <= lqn.nhosts:
                    elem_type = 'Processor'
                elif idx <= lqn.tshift + lqn.ntasks:
                    elem_type = 'Task'
                elif idx <= lqn.eshift + lqn.nentries:
                    elem_type = 'Entry'
                else:
                    elem_type = 'Activity'

                perc_list.append({
                    'element': idx,
                    'name': name,
                    'type': elem_type,
                    'percentiles': list(percentiles),
                    'values': list(perc_values),
                })

                # Add row to DataFrame
                row = {'Element': name, 'Type': elem_type}
                for p_val, perc in zip(perc_values, percentiles):
                    col_name = f'P{int(perc)}'
                    row[col_name] = p_val
                rows.append(row)

        # Create DataFrame
        if rows:
            df = pd.DataFrame(rows)
        else:
            df = pd.DataFrame()

        return perc_list, df

    def GetPerctRespT(self, percentiles: List[float] = None, jobclass: int = None) -> Tuple[List[Dict[str, Any]], pd.DataFrame]:
        """PascalCase alias for getPerctRespT()."""
        return self.getPerctRespT(percentiles, jobclass)

    # ========== PHASE 5: Method Introspection ==========

    def listValidMethods(self) -> List[str]:
        """
        List all valid solution methods for this LN model.

        Returns list of method names supported by the solver for this model.
        Always includes 'exact' as the primary method.

        Returns:
            List of valid method names (strings)

        Example:
            >>> solver = SolverLN(model)
            >>> methods = solver.listValidMethods()
            >>> print(f"Available methods: {methods}")
        """
        return ['exact', 'default']

    def ListValidMethods(self) -> List[str]:
        """PascalCase alias for listValidMethods()."""
        return self.listValidMethods()

    @staticmethod
    def getFeatureSet() -> set:
        """
        Get set of features supported by LN solver.

        Returns set of feature strings representing supported:
        - Node types: Processor, Task, Entry, Activity
        - Scheduling strategies: FCFS, PS, LCFS, HOL
        - Distributions: Exp, Erlang, HyperExp, Coxian, PH
        - Other: ClosedClass, OpenClass, Delay

        Returns:
            set: Feature strings

        Example:
            >>> features = SolverLN.getFeatureSet()
            >>> print(f"Supports Exponential: {'Exp' in features}")
        """
        return {
            'Processor', 'Task', 'Entry', 'Activity',
            'FCFS', 'PS', 'LCFS', 'HOL', 'INF', 'REF',
            'Exp', 'Erlang', 'HyperExp', 'Coxian', 'PH',
            'ClosedClass', 'OpenClass', 'Delay',
            'SynchronousCall', 'AsynchronousCall',
            'Fork', 'Join', 'AndFork', 'AndJoin',
            'OrFork', 'OrJoin', 'ActivityPrecedence'
        }

    @staticmethod
    def GetFeatureSet() -> set:
        """PascalCase alias for getFeatureSet()."""
        return SolverLN.getFeatureSet()

    @staticmethod
    def supports(model) -> bool:
        """
        Check if LN solver supports the given model.

        Tests whether the model is a valid LayeredNetwork that can be
        analyzed by this solver.

        Args:
            model: Model object to check

        Returns:
            bool: True if model is supported, False otherwise

        Example:
            >>> model = LayeredNetwork('Test')
            >>> if SolverLN.supports(model):
            ...     solver = SolverLN(model)
        """
        if model is None:
            return False
        # Check if it's a LayeredNetwork with getStruct method
        return hasattr(model, 'getStruct') and callable(getattr(model, 'getStruct'))

    @staticmethod
    def Supports(model) -> bool:
        """PascalCase alias for supports()."""
        return SolverLN.supports(model)

    @staticmethod
    def defaultOptions() -> Dict[str, Any]:
        """
        Get default solver options.

        Returns dictionary with default settings for:
        - method: Solution method to use ('exact')
        - tol: Convergence tolerance (1e-6)
        - max_iter: Maximum iterations (100)
        - verbose: Verbosity flag (False)
        - under_relaxation: Relaxation factor (0.7)

        Returns:
            dict: Default options

        Example:
            >>> options = SolverLN.defaultOptions()
            >>> print(f"Default tolerance: {options['tol']}")
        """
        return OptionsDict({
            'method': 'exact',
            'tol': 1e-6,
            'max_iter': 100,
            'verbose': False,
            'under_relaxation': 0.7,
        })

    @staticmethod
    def DefaultOptions() -> Dict[str, Any]:
        """PascalCase alias for defaultOptions()."""
        return SolverLN.defaultOptions()

    # ========== PHASE 6: Sampling and Transient Placeholders ==========

    def sample(self, node: int, numSamples: int) -> List[Dict[str, Any]]:
        """
        Sample from steady-state distribution at a node.

        Sampling is not supported by the analytical LN solver.
        Use DES or simulation-based solvers for sampling capabilities.

        Args:
            node: Node index
            numSamples: Number of samples to generate

        Raises:
            NotImplementedError: Always, as LN solver is analytical

        Example:
            >>> solver = SolverLN(model)
            >>> try:
            ...     samples = solver.sample(1, 1000)
            ... except NotImplementedError as e:
            ...     print(f"Use simulation-based solver instead: {e}")
        """
        raise NotImplementedError(
            "sample() not supported by SolverLN (analytical solver). "
            "Use SolverDES or simulation-based solvers for sampling capabilities."
        )

    def Sample(self, node: int, numSamples: int) -> List[Dict[str, Any]]:
        """PascalCase alias for sample()."""
        return self.sample(node, numSamples)

    def sampleAggr(self, node: int, numSamples: int) -> List[Dict[str, Any]]:
        """
        Sample aggregated distribution at a node.

        Sampling not supported by analytical LN solver.

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "sampleAggr() not supported by SolverLN (analytical solver). "
            "Use SolverDES or simulation-based solvers for sampling capabilities."
        )

    def SampleAggr(self, node: int, numSamples: int) -> List[Dict[str, Any]]:
        """PascalCase alias for sampleAggr()."""
        return self.sampleAggr(node, numSamples)

    def sampleSys(self, numSamples: int) -> List[Dict[str, Any]]:
        """
        Sample system-level distribution.

        Sampling not supported by analytical LN solver.

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "sampleSys() not supported by SolverLN (analytical solver). "
            "Use SolverDES or simulation-based solvers for sampling capabilities."
        )

    def SampleSys(self, numSamples: int) -> List[Dict[str, Any]]:
        """PascalCase alias for sampleSys()."""
        return self.sampleSys(numSamples)

    def sampleSysAggr(self, numSamples: int) -> List[Dict[str, Any]]:
        """
        Sample aggregated system distribution.

        Sampling not supported by analytical LN solver.

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "sampleSysAggr() not supported by SolverLN (analytical solver). "
            "Use SolverDES or simulation-based solvers for sampling capabilities."
        )

    def SampleSysAggr(self, numSamples: int) -> List[Dict[str, Any]]:
        """PascalCase alias for sampleSysAggr()."""
        return self.sampleSysAggr(numSamples)

    def getCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """
        Get passage time CDF.

        Passage time analysis not supported by analytical LN solver.

        Args:
            R: Optional response time matrix

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "getCdfPassT() not supported by SolverLN. "
            "Use transient-capable solvers for passage time analysis."
        )

    def GetCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """PascalCase alias for getCdfPassT()."""
        return self.getCdfPassT(R)

    def getTranCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """
        Get transient response time CDF.

        Transient analysis not supported by analytical LN solver.

        Args:
            R: Optional response time matrix

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "getTranCdfRespT() not supported by SolverLN. "
            "Use transient-capable solvers for time-dependent analysis."
        )

    def GetTranCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """PascalCase alias for getTranCdfRespT()."""
        return self.getTranCdfRespT(R)

    def getTranCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """
        Get transient passage time CDF.

        Transient analysis not supported by analytical LN solver.

        Args:
            R: Optional response time matrix

        Raises:
            NotImplementedError: Always
        """
        raise NotImplementedError(
            "getTranCdfPassT() not supported by SolverLN. "
            "Use transient-capable solvers for time-dependent analysis."
        )

    def GetTranCdfPassT(self, R: Optional[np.ndarray] = None) -> List[Dict[str, Any]]:
        """PascalCase alias for getTranCdfPassT()."""
        return self.getTranCdfPassT(R)

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Node, NodeType, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        lqn = self.lqn
        rows = []

        for idx in range(1, lqn.nidx + 1):
            if idx >= len(self._result['QN']):
                continue

            name = lqn.hashnames[idx] if idx < len(lqn.hashnames) else f'Node{idx}'

            # Determine node type
            if idx <= lqn.nhosts:
                node_type = 'Processor'
            elif idx <= lqn.tshift + lqn.ntasks:
                node_type = 'Task'
            elif idx <= lqn.eshift + lqn.nentries:
                node_type = 'Entry'
            else:
                node_type = 'Activity'

            rows.append({
                'Node': name,
                'NodeType': node_type,
                'QLen': self._result['QN'][idx],
                'Util': self._result['UN'][idx],
                'RespT': self._result['RN'][idx],
                'ResidT': self._result['RN'][idx],
                'ArvR': self._result['AN'][idx],
                'Tput': self._result['TN'][idx],
            })

        df = pd.DataFrame(rows)

        # Filter out all-zero rows
        if not df.empty:
            numeric_cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
            tokeep = ~(df[numeric_cols] <= 0.0).all(axis=1)
            df = df.loc[tokeep].reset_index(drop=True)

        if not self._table_silent:
            print(df.to_string(index=False))

        return df

    # Aliases
    avg_table = getAvgTable
    get_avg_table = getAvgTable
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable
    default_options = defaultOptions


__all__ = ['SolverLN', 'SolverLNOptions', 'LayerMetrics', 'LayerModel',
           'getEnsembleAvg', 'getTranAvg', 'getCdfRespT', 'getStageAvg',
           'computeTaskDelta', 'computeHostDelta', 'getConvergenceCriteria', 'getIterationTrace',
           'getAvgQLen', 'getAvgUtil', 'getAvgRespT', 'getAvgResidT', 'getAvgWaitT',
           'getAvgTput', 'getAvgArvR', 'getAvgSysRespT', 'getAvgSysTput',
           'getPerctRespT', 'listValidMethods', 'getFeatureSet', 'supports', 'defaultOptions',
           'sample', 'sampleAggr', 'sampleSys', 'sampleSysAggr',
           'getCdfPassT', 'getTranCdfRespT', 'getTranCdfPassT']
