"""
Native Python implementation of SSA (Stochastic Simulation Algorithm) solver.

This implementation uses pure Python/NumPy for stochastic simulation of queueing
networks.
"""

import math
import os
import numpy as np
import pandas as pd
import sys
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass, field, fields as dataclass_fields
from ...constants import default_verbose

from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn.getters import sn_get_arvr_from_tput
from ...api.sn.network_struct import NodeType
from ..fjtag_transform import FJTagTransformMixin
from ...api.io.logging import line_warning
from ...api.io.logging import line_debug
from ...constants import GlobalConstants
from ..base import NetworkSolver, method_label, method_type
from ...indexed_table import IndexedTable


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
class SampleEvent:
    """A single event in an SSA sample path."""
    t: float         # Event time
    node: int         # Node index (1-based, matching MATLAB convention)
    class_idx: int    # Class index (1-based, matching MATLAB convention)
    event: str        # 'ARV' or 'DEP'


@dataclass
class SamplePath:
    """Sample path from SSA simulation, matching MATLAB sampleNodeState."""
    handle: Any = None          # Reference to sampled node
    t: Optional[np.ndarray] = None      # Event times
    state: Optional[np.ndarray] = None  # States at event times
    event: Optional[List] = field(default_factory=list)  # List of SampleEvent
    isaggregate: bool = False


@dataclass
class SolverSSAOptions:
    """Options for the native SSA solver."""
    method: str = 'default'
    tol: float = 1e-4
    samples: int = 10000  # Deprecated alias of events: fired transitions to simulate
    events: Optional[int] = None  # DES event budget; overrides samples when set
    warmupfrac: float = 0.0  # Warmup discard fraction for steady-state tallies (0 = disabled)
    seed: int = 0
    cutoff: float = np.inf  # Match MATLAB default: infinity means no truncation
    confidence_level: float = 0.95
    verbose: bool = field(default_factory=default_verbose)
    keep: bool = False  # Keep intermediate data (for compatibility)
    record_events: bool = False  # Record per-event log for sample paths
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))  # env LINE_SOLVER_LANG overrides; 'python' (native), 'java' (jline.jar via JSON) or 'cpp' (line-cli via JSON)
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. Meaningless for the other langs, which are IEEE double
    # throughout, so line-cli is invoked without --arith unless the caller sets it.
    # A C++ sample path is drawn from exponential clocks, so only double is served.
    arith: Optional[str] = None


class SolverSSA(FJTagTransformMixin, NetworkSolver):
    """
    Native Python SSA (Stochastic Simulation Algorithm) solver.

    This solver analyzes queueing networks through discrete-event simulation
    using Gillespie's algorithm in pure Python/NumPy.

    Supported methods:
        - 'default': Serial Gillespie simulation
        - 'serial': Same as default
        - 'ssa': Same as default

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default')
        **kwargs: Additional solver options
    """

    @staticmethod
    def _adopt_options(src, kwargs, dict_like):
        """Carry every SolverSSAOptions field across from an options object.

        Forwarding a hand-picked subset silently pins the rest to their
        defaults, so SolverSSA(model, opts) with opts.samples = 2e6 would run
        the default 10000 samples and merely look converged. Enumerating the
        dataclass fields keeps the two in step as fields are added. Fields the
        source does not carry are left alone, and names outside
        SolverSSAOptions are ignored, so a generic SolverOptions with a wider
        surface does not break the constructor below.
        """
        for f in dataclass_fields(SolverSSAOptions):
            if f.name == 'method':
                continue  # resolved by the caller, which also lowercases it
            if dict_like:
                if f.name in src:
                    kwargs.setdefault(f.name, src[f.name])
            elif hasattr(src, f.name):
                kwargs.setdefault(f.name, getattr(src, f.name))

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._sn = None

        # see _kb/09-ldes-and-cache.md (Warm start) for initFromSolver contract
        init_solver = None
        if method_or_options is not None and hasattr(method_or_options, 'getAvgQLen'):
            init_solver = method_or_options
            method_or_options = None

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            self.method = str(kwargs.get('method', 'default')).lower()
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'default')
            self._adopt_options(method_or_options, kwargs, dict_like=True)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            self._adopt_options(method_or_options, kwargs, dict_like=False)
        else:
            self.method = 'default'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverSSAOptions(method=self.method, **kwargs)

        # options.events (DES event budget) overrides options.samples when set;
        # samples remains accepted as a deprecated alias for the event budget.
        if self.options.events:
            self.options.samples = int(self.options.events)

        # Extract network structure
        self._extract_network_params()

        if init_solver is not None:
            self.initFromSolver(init_solver)

    def getName(self) -> str:
        """Get the name of this solver."""
        return "SSA"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model."""
        model = self.model

        # Priority 1: Native model with _sn
        if hasattr(model, '_sn') and model._sn is not None:
            self._sn = model._sn
            return

        # Priority 2: Native model with refresh_struct
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn'):
                self._sn = model._sn
                return

        # Priority 3: native model (snake-case get_struct()); no JAR-wrapper bridge.
        if hasattr(model, 'get_struct'):
            self._sn = model.get_struct()
            if self._sn is not None:
                return

        # Priority 4: Already a native NetworkStruct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            self._sn = model
            return

        raise ValueError(
            "Cannot extract a native NetworkStruct from model. Native solvers "
            "accept only native Network / NetworkStruct inputs (no JAR wrapper).")

    def runAnalyzer(self) -> 'SolverSSA':
        """Run the SSA analysis."""
        # see _kb/06-solver-catalog.md (Wrappers: "Python lang='java' opt-in JAR delegation")
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import populate_java_result
            populate_java_result(self)
            return self
        # see _kb/06-solver-catalog.md ("Python lang='cpp' opt-in C++ delegation");
        # an absent binary is the only automatic fallback, a C++ refusal propagates.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, populate_cpp_result
            try:
                populate_cpp_result(self)
                return self
            except LineCliNotAvailable as e:
                line_warning("SolverSSA", "lang='cpp' requested but the C++ solver is "
                             "unavailable (%s); falling back to lang='python'." % e)

        line_debug("SSA: using lang=python", options=self.options)

        # see _kb/06-solver-catalog.md (SSA: "Serial engine: finite capacity
        # regions (FCR), fork-join, immediate feedback")
        self._fj_foldback = None
        if self._sn is not None and getattr(self._sn, 'nodetype', None) is not None:
            _ntv = np.asarray(self._sn.nodetype).ravel()
            _has_fj = bool(np.any(_ntv == NodeType.FORK) or np.any(_ntv == NodeType.JOIN))
            if _has_fj:
                self._fjtag_require_network('SSA')
                self._sn = self._fjtag_expand(self._sn)

        # see _kb/06-solver-catalog.md (SSA main section) for FCR/SPN/trace-driven support
        from ...constants import ProcessType as _PT_early
        _trace_types = {getattr(_PT_early, n) for n in ('REPLAYER', 'TRACE')
                        if hasattr(_PT_early, n)}
        if _trace_types and getattr(self._sn, 'procid', None) is not None:
            for _row in np.atleast_2d(self._sn.procid):
                for _pid in _row:
                    if _pid in _trace_types:
                        raise RuntimeError(
                            "This model uses a Replayer/Trace (trace-driven) distribution, "
                            "which SolverSSA cannot simulate (its rate-driven algorithm has no "
                            "trace representation). Use SolverLDES for trace-driven models.")

        # see _kb/06-solver-catalog.md (Wrappers: "Python LDES wrapper:
        # Breakdown rejection") -- same rationale applies to native SSA
        _hasbd = getattr(self._sn, 'hasbreakdown', None)
        if _hasbd is not None and np.any(np.asarray(_hasbd).ravel() == 1):
            _bd_nodes = [str(self._sn.nodenames[_i])
                         for _i in np.nonzero(np.asarray(_hasbd).ravel() == 1)[0]]
            raise RuntimeError(
                "Station(s) %s declare server breakdowns, which SolverSSA does not "
                "simulate: the joint (queue, server status) chain is expanded only by "
                "SolverCTMC. Use SolverCTMC for models with set_breakdown."
                % ', '.join(_bd_nodes))

        # Reneging/retrial are first-class sync events (refresh_sync /
        # after_event_station); the serial engine simulates them directly.
        from ...lang.base import SchedStrategy as _SchedStrategy
        sn = self._sn

        # Signal classes never occupy a station; cap capacity at 0 except at
        # EXT/Source, mirroring solver_ctmc.py and MATLAB solver_ssa.m.
        if (sn is not None and getattr(sn, 'issignal', None) is not None
                and getattr(sn, 'classcap', None) is not None):
            for _ist in range(int(sn.nstations)):
                if sn.sched[_ist] != _SchedStrategy.EXT:
                    for _r in range(int(sn.nclasses)):
                        if sn.issignal[_r]:
                            sn.classcap[_ist, _r] = 0

        # see _kb/06-solver-catalog.md (SSA: "Serial engine: heterogeneous
        # servers, PAS, signal classes")
        _mdl = getattr(self, 'model', None)
        if sn is not None and _mdl is not None:
            _stns = getattr(_mdl, '_stations', None) or []
            _cls = getattr(_mdl, '_classes', None) or []
            for _ist, _st in enumerate(_stns):
                if not (hasattr(_st, 'is_heterogeneous') and _st.is_heterogeneous()):
                    continue
                _served = [r for r in range(int(sn.nclasses)) if float(sn.rates[_ist, r]) > 0]
                if len(_served) > 1:
                    raise RuntimeError(
                        "SolverSSA supports heterogeneous servers only for single-class "
                        "stations. Use SolverJMT or SolverLDES for multi-class heterogeneous servers.")
                if len(_served) != 1:
                    continue
                _r = _served[0]
                _jc = _cls[_r]
                _srvrates = []
                for _sty in _st.get_server_types():
                    _d = _st.get_hetero_service(_jc, _sty)
                    if _d is not None:
                        _gm = getattr(_d, 'getMean', None) or getattr(_d, 'get_mean', None)
                        _mval = _gm() if _gm is not None else None
                        if _mval and _mval > 0:
                            _srvrates += [1.0 / _mval] * int(_sty.get_num_of_servers())
                _c = len(_srvrates)
                _mu_base = float(sn.rates[_ist, _r])
                if _c > 0 and _mu_base > 0:
                    _njobs_sum = int(sum(int(nj) for nj in np.ravel(sn.njobs) if np.isfinite(nj)))
                    _Lh = max(_c, _njobs_sum, 1)
                    if sn.lldscaling is None or (hasattr(sn.lldscaling, 'size') and sn.lldscaling.size == 0):
                        sn.lldscaling = np.ones((int(sn.nstations), _Lh))
                    elif sn.lldscaling.shape[1] < _c:
                        _ext = np.ones((int(sn.nstations), _c))
                        _ext[:, :sn.lldscaling.shape[1]] = sn.lldscaling
                        for _b in range(sn.lldscaling.shape[1], _c):
                            _ext[:, _b] = np.ravel(sn.lldscaling[:, -1])
                        sn.lldscaling = _ext
                    for _n in range(1, sn.lldscaling.shape[1] + 1):
                        _lim = min(_n, _c)
                        sn.lldscaling[_ist, _n - 1] = sum(_srvrates[:_lim]) / (_mu_base * _lim)

        from ...api.solvers.ssa.handler import (
            solver_ssa, SolverSSAOptions as HandlerOptions
        )

        # Convert non-Markovian distributions to PH (matches
        # solver_ssa_analyzer.m:18); work on a private copy, as the CTMC path does.
        from ...constants import ProcessType
        _markovian = {ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP,
                      ProcessType.PH, ProcessType.APH, ProcessType.MAP,
                      ProcessType.COXIAN, ProcessType.COX2, ProcessType.MMPP2,
                      ProcessType.IMMEDIATE, ProcessType.DISABLED}
        _needs_toph = False
        if getattr(self._sn, 'procid', None) is not None:
            for _row in np.atleast_2d(self._sn.procid):
                for _pid in _row:
                    if _pid is None or (isinstance(_pid, float) and np.isnan(_pid)):
                        continue
                    if isinstance(_pid, ProcessType) and _pid not in _markovian:
                        _needs_toph = True
                        break
                if _needs_toph:
                    break
        if _needs_toph:
            import copy
            from ...api.sn import sn_nonmarkov_toph
            opts_dict = {}
            if hasattr(self.options, '__dict__'):
                opts_dict = {k: v for k, v in self.options.__dict__.items() if not k.startswith('_')}
            elif isinstance(self.options, dict):
                opts_dict = dict(self.options)
            # SSA draws a sample path, so its surrogate must be a genuine
            # phase-type: a matrix exponential has no sample path to draw.
            cfg = dict(opts_dict.get('config') or {})
            cfg['phfit'] = 'ph'
            opts_dict['config'] = cfg
            self._sn = sn_nonmarkov_toph(copy.deepcopy(self._sn), opts_dict)

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            samples=self.options.samples,
            seed=self.options.seed,
            cutoff=self.options.cutoff,
            confidence_level=self.options.confidence_level,
            verbose=self.options.verbose,
            record_events=self.options.record_events,
            warmupfrac=getattr(self.options, 'warmupfrac', 0.0),
            timeout=getattr(self.options, 'timeout', float('inf'))
        )

        # Run the solver - pass model for cache node details
        self._result = solver_ssa(self._sn, handler_options, self.model)

        # Fold FJ auxiliary sibling-class columns back; recompute RN=QN/AN;
        # restore the original struct for tables.
        if getattr(self, '_fj_foldback', None) is not None:
            r = self._result
            _korig = self._fj_foldback[1]
            AN = self._fjtag_lift(r)
            # FJ siblings land in auxiliary classes; report the routing-derived
            # arrival rate on the original (folded) struct instead.
            r.A = AN[:, :_korig]
        else:
            from ...api.sn.getters import sn_pn_avg_rates
            r = self._result
            r.T, r.A, r.R = sn_pn_avg_rates(self._sn, r.Q, r.T, r.A, r.R)

        # Store cache statistics if available
        if hasattr(self._result, 'cache_stats') and self._result.cache_stats:
            self._cache_stats = self._result.cache_stats

        line_debug("SSA analysis complete: extracting results (nstations=%d, nclasses=%d)",
                   self._sn.nstations if self._sn else 0, self._sn.nclasses if self._sn else 0, options=self.options)

        # Extract names
        self._extract_names()

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            runtime = self._result.runtime if hasattr(self._result, 'runtime') else 0.0
            method = self._result.method if hasattr(self._result, 'method') else 'serial'
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"SSA analysis [method: {method_label(self.options.method, method)}; type: {method_type('SSA', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")

        return self

    def _adjust_cache_results(self):
        """
        Adjust SSA results for networks with cache nodes.

        NOTE: This method is kept for backwards compatibility but is now
        largely a no-op since proper cache simulation is done in the handler.
        The handler properly simulates cache state with LRU and Zipf distribution.
        """
        if self._result is None or self._sn is None:
            return

        sn = self._sn
        I = sn.nnodes
        K = sn.nclasses
        M = sn.nstations

        # Find cache nodes
        cache_indices = []
        for ind in range(I):
            if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
                cache_indices.append(ind)

        if not cache_indices:
            return

        # Get model nodes and classes
        model_nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        model_classes = self.model.get_classes() if hasattr(self.model, 'get_classes') else []

        if not model_nodes or not model_classes:
            return

        # Compute and apply hit/miss probabilities for each cache
        for ind in cache_indices:
            cache_node = model_nodes[ind] if ind < len(model_nodes) else None
            if cache_node is None or not hasattr(cache_node, 'get_gamma_matrix'):
                continue

            gamma = cache_node.get_gamma_matrix(K)
            m_levels = cache_node._item_level_cap if hasattr(cache_node, '_item_level_cap') else np.array([1])

            # Compute hit/miss probabilities using FPI
            try:
                xi, pi0, pij, it = cache_xi_fp(gamma, m_levels)
                access_probs = np.sum(gamma, axis=1)
                total = np.sum(access_probs)
                if total > 0:
                    access_probs = access_probs / total
                hit_rate = np.sum(access_probs * (1 - pi0))
                miss_rate = 1 - hit_rate
            except Exception:
                continue

            # Store in cache node
            hp = np.zeros(K)
            mp = np.zeros(K)
            input_classes = []

            for k, job_class in enumerate(model_classes):
                h = cache_node._hit_class.get(job_class) if hasattr(cache_node, '_hit_class') else None
                m = cache_node._miss_class.get(job_class) if hasattr(cache_node, '_miss_class') else None
                if h is not None and m is not None:
                    hp[k] = hit_rate
                    mp[k] = miss_rate
                    h_idx = model_classes.index(h) if h in model_classes else -1
                    m_idx = model_classes.index(m) if m in model_classes else -1
                    if h_idx >= 0 and m_idx >= 0:
                        input_classes.append((k, h_idx, m_idx))

            if hasattr(cache_node, 'set_result_hit_prob'):
                cache_node.set_result_hit_prob(hp)
                cache_node.set_result_miss_prob(mp)

            # Adjust throughputs and queue lengths
            for (in_k, h_idx, m_idx) in input_classes:
                # Find station that serves hit/miss classes
                for ist in range(M):
                    T_hit = self._result.T[ist, h_idx]
                    T_miss = self._result.T[ist, m_idx]

                    if T_hit > 0 or T_miss > 0:
                        # Get response times
                        R_hit = self._result.R[ist, h_idx] if self._result.R[ist, h_idx] > 0 else 0.2
                        R_miss = self._result.R[ist, m_idx] if self._result.R[ist, m_idx] > 0 else 1.0

                        # see _kb/06-solver-catalog.md (SSA: "Warmup discard,
                        # and cache hit/miss accounting")
                        correct_cycle_time = hit_rate * R_hit + miss_rate * R_miss

                        # For closed network with N=1 job
                        correct_X_total = 1.0 / correct_cycle_time if correct_cycle_time > 0 else 1.0

                        # Compute correct throughputs
                        new_T_hit = correct_X_total * hit_rate
                        new_T_miss = correct_X_total * miss_rate

                        self._result.T[ist, h_idx] = new_T_hit
                        self._result.T[ist, m_idx] = new_T_miss

                        # Adjust queue lengths using Little's law
                        Q_hit = new_T_hit * R_hit
                        Q_miss = new_T_miss * R_miss
                        self._result.Q[ist, h_idx] = Q_hit
                        self._result.Q[ist, m_idx] = Q_miss
                        self._result.U[ist, h_idx] = Q_hit  # For delay node
                        self._result.U[ist, m_idx] = Q_miss

                        # Also update throughput at other stations for the input class
                        for jst in range(M):
                            if jst != ist and self._result.T[jst, in_k] > 0:
                                self._result.T[jst, in_k] = correct_X_total

    def _compute_retrieval_latency_ssa(self, sn, cache_param, cache_ind, TN, QN):
        """
        Retrieval-system expected latency (Sala et al. 2026 Eq. 8) is not
        currently implemented. For any cache configured with a retrieval
        system this returns all-NaN latencies (and emits a one-shot warning);
        it returns (None, None) when the cache has no retrieval system.
        """
        K = sn.nclasses
        rpc = getattr(cache_param, 'retrieval_classes', None)
        rsqi = getattr(cache_param, 'retrieval_system_queue_indices', None)
        if rpc is None or rsqi is None:
            return None, None
        rpc_arr = np.atleast_2d(np.asarray(rpc))
        if rpc_arr.size == 0 or rpc_arr.shape[1] == 0:
            return None, None
        # The Eq. 8 retrieval-system expected latency is not currently
        # implemented; report NaN whenever a retrieval system is configured.
        if isinstance(rsqi, dict) and any(rsqi.get(k) for k in range(K)):
            line_warning('solver_ssa_analyzer_serial',
                         'Retrieval-system expected latency is not currently '
                         'implemented; reporting NaN.')
        return np.full(K, np.nan), None

    def _extract_names(self):
        """Extract station and class names using stationToNode mapping."""
        if self._sn is None:
            self.station_names = []
            self.class_names = []
            return

        nstations = self._sn.nstations

        # Get station names using stationToNode mapping
        nodenames = list(self._sn.nodenames) if hasattr(self._sn, 'nodenames') and self._sn.nodenames else []
        stationToNode = self._sn.stationToNode if hasattr(self._sn, 'stationToNode') else None

        if stationToNode is not None and nodenames:
            stationToNode = np.asarray(stationToNode).flatten()
            self.station_names = []
            for i in range(nstations):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        self.station_names.append(nodenames[node_idx])
                    else:
                        self.station_names.append(f'Station{i}')
                else:
                    self.station_names.append(f'Station{i}')
        else:
            self.station_names = [f'Station{i}' for i in range(nstations)]

        # Get class names
        self.class_names = list(self._sn.classnames) if hasattr(self._sn, 'classnames') and self._sn.classnames else \
                          [f'Class{i}' for i in range(self._sn.nclasses)]

    # =========================================================================
    # Table Output
    # =========================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """Get performance metrics table."""
        if self._result is None:
            self._ensureAvgResults()

        nstations = self._result.Q.shape[0]
        nclasses = self._result.Q.shape[1]

        # see _kb/06-solver-catalog.md (SSA: "Warmup discard, and cache hit/miss accounting")
        sn = self._sn
        if sn is not None and sn.nodeparam is not None:
            TN = self._result.T
            I = sn.nnodes
            R = sn.nclasses
            has_cache_update = False
            for ind in range(I):
                if sn.nodetype is not None and ind < len(sn.nodetype):
                    if sn.nodetype[ind] == NodeType.CACHE and ind in sn.nodeparam:
                        cache_param = sn.nodeparam[ind]

                        # Prefer actualhitprob already computed by handler.py
                        # over recomputing from station throughputs (fails for non-station cache).
                        existing_hitprob = getattr(cache_param, 'actualhitprob', None)
                        if existing_hitprob is not None and np.any(existing_hitprob > 0):
                            existing_missprob = getattr(cache_param, 'actualmissprob', None)
                            if existing_missprob is None:
                                existing_missprob = 1.0 - existing_hitprob
                            # Set result on Cache node in model
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes') and ind < len(self.model._nodes):
                                cache_node = self.model._nodes[ind]
                                if hasattr(cache_node, 'set_result_hit_prob'):
                                    cache_node.set_result_hit_prob(existing_hitprob)
                                if hasattr(cache_node, 'set_result_miss_prob'):
                                    cache_node.set_result_miss_prob(existing_missprob)
                            has_cache_update = True
                            continue

                        # Fallback: compute from station throughputs
                        hitclass = getattr(cache_param, 'hitclass', np.array([]))
                        missclass = getattr(cache_param, 'missclass', np.array([]))
                        nclass = len(hitclass) if hasattr(hitclass, '__len__') else R
                        ist = int(sn.nodeToStation[ind]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and ind < len(sn.nodeToStation) else -1

                        cache_param.actualhitprob = np.zeros(nclass)
                        cache_param.actualmissprob = np.zeros(nclass)

                        for k in range(nclass):
                            h = int(hitclass[k]) if k < len(hitclass) else -1
                            m = int(missclass[k]) if k < len(missclass) else -1
                            if h >= 0 and m >= 0 and ist >= 0 and ist < TN.shape[0]:
                                t_hit = TN[ist, h] if h < TN.shape[1] else 0
                                t_miss = TN[ist, m] if m < TN.shape[1] else 0
                                t_total = t_hit + t_miss
                                if t_total > 0:
                                    cache_param.actualhitprob[k] = t_hit / t_total
                                    cache_param.actualmissprob[k] = t_miss / t_total
                                    has_cache_update = True

                        # Set result on Cache node in model
                        if hasattr(self, 'model') and hasattr(self.model, '_nodes') and ind < len(self.model._nodes):
                            cache_node = self.model._nodes[ind]
                            if hasattr(cache_node, 'set_result_hit_prob'):
                                cache_node.set_result_hit_prob(cache_param.actualhitprob)
                            if hasattr(cache_node, 'set_result_miss_prob'):
                                cache_node.set_result_miss_prob(cache_param.actualmissprob)

            # Refresh chains to recompute visits with actual hit/miss probabilities
            if has_cache_update and hasattr(self, 'model') and hasattr(self.model, '_refresh_chains'):
                self.model._refresh_chains()
                if hasattr(self.model, '_sn'):
                    sn = self.model._sn
                    self._sn = sn

        # Make copies of station-level metrics (to avoid modifying originals)
        QN = self._result.Q.copy()
        UN = self._result.U.copy()
        RN = self._result.R.copy()
        TN = self._result.T.copy()

        # Zero out metrics for classes that don't visit stations based on visit ratios
        # This matches MATLAB getAvg.m behavior (lines 163-180)
        hasForkJoin = False
        hasSPN = False
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            hasForkJoin = np.any(sn.nodetype == NodeType.FORK) and np.any(sn.nodetype == NodeType.JOIN)
            hasSPN = np.any(sn.nodetype == NodeType.PLACE) or np.any(sn.nodetype == NodeType.TRANSITION)

        if sn is not None and hasattr(sn, 'nchains') and sn.nchains > 0 and not hasSPN:
            if hasattr(sn, 'chains') and sn.chains is not None and hasattr(sn, 'visits') and sn.visits:
                chains_arr = np.asarray(sn.chains)
                # Get station-to-stateful mapping (visits are indexed by stateful node, not station)
                stationToStateful = None
                if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
                    stationToStateful = np.asarray(sn.stationToStateful).flatten().astype(int)
                for k in range(nclasses):
                    # Find chains containing this class
                    chains_with_class = np.where(chains_arr[:, k] > 0)[0] if k < chains_arr.shape[1] else []
                    if len(chains_with_class) > 0:
                        c = chains_with_class[0]  # Use first chain (classes typically in one chain)
                        if c in sn.visits and sn.visits[c] is not None:
                            visits_c = np.asarray(sn.visits[c])
                            for i in range(nstations):
                                # Convert station index to stateful index for visits lookup
                                stateful_idx = stationToStateful[i] if stationToStateful is not None and i < len(stationToStateful) else i
                                if stateful_idx < visits_c.shape[0] and k < visits_c.shape[1]:
                                    if visits_c[stateful_idx, k] == 0:
                                        # For fork-join, trust non-zero simulation results
                                        if hasForkJoin and (QN[i, k] > GlobalConstants.FineTol or
                                                           UN[i, k] > GlobalConstants.FineTol or
                                                           TN[i, k] > GlobalConstants.FineTol):
                                            continue
                                        # Zero out station-level metrics
                                        QN[i, k] = 0
                                        UN[i, k] = 0
                                        RN[i, k] = 0
                                        TN[i, k] = 0

        # Recompute ResidT and ArrivalRate after zeroing
        if sn is not None and sn.visits:
            WN = sn_get_residt_from_respt(sn, RN, None)
        else:
            WN = RN.copy()

        # Get arrival rates from result or compute from throughputs
        if self._result.A is not None:
            AN = self._result.A.copy()
            # Also zero AN for non-visiting classes
            if sn is not None and hasattr(sn, 'nchains') and sn.nchains > 0 and not hasSPN:
                if hasattr(sn, 'chains') and sn.chains is not None and hasattr(sn, 'visits') and sn.visits:
                    chains_arr = np.asarray(sn.chains)
                    # Get station-to-stateful mapping (visits are indexed by stateful node, not station)
                    stationToStateful = None
                    if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
                        stationToStateful = np.asarray(sn.stationToStateful).flatten().astype(int)
                    for k in range(nclasses):
                        chains_with_class = np.where(chains_arr[:, k] > 0)[0] if k < chains_arr.shape[1] else []
                        if len(chains_with_class) > 0:
                            c = chains_with_class[0]
                            if c in sn.visits and sn.visits[c] is not None:
                                visits_c = np.asarray(sn.visits[c])
                                for i in range(nstations):
                                    # Convert station index to stateful index for visits lookup
                                    stateful_idx = stationToStateful[i] if stationToStateful is not None and i < len(stationToStateful) else i
                                    if stateful_idx < visits_c.shape[0] and k < visits_c.shape[1]:
                                        if visits_c[stateful_idx, k] == 0 and not (hasForkJoin and AN[i, k] > GlobalConstants.FineTol):
                                            AN[i, k] = 0
        else:
            AN = sn_get_arvr_from_tput(sn, TN, None)

        # Identify source and cache stations
        source_stations = set()
        cache_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None
        if nodetype is not None:
            stationToNode = self._sn.stationToNode
            if stationToNode is not None:
                stationToNode = np.asarray(stationToNode).flatten()
                nodetype_arr = np.asarray(nodetype) if not isinstance(nodetype, np.ndarray) else nodetype
                for i in range(nstations):
                    if i < len(stationToNode):
                        node_idx = int(stationToNode[i])
                        if node_idx < len(nodetype_arr):
                            nt = nodetype_arr[node_idx]
                            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
                            if nt_val == 0:  # SOURCE
                                source_stations.add(i)
                            elif nt_val == 6:  # CACHE
                                cache_stations.add(i)

        rows = []
        for i in range(nstations):
            # Skip cache stations (like MATLAB)
            if i in cache_stations:
                continue

            for r in range(nclasses):
                station_name = self.station_names[i] if i < len(self.station_names) else f'Station{i}'
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                qlen = QN[i, r]
                util = UN[i, r]
                respt = RN[i, r]
                residt = WN[i, r] if i < WN.shape[0] and r < WN.shape[1] else respt
                tput = TN[i, r]
                arvr = AN[i, r] if i < AN.shape[0] and r < AN.shape[1] else tput

                is_source = i in source_stations

                # Source stations are not queues: QLen=Util=RespT=ResidT=ArvR=0, Tput=arrival rate
                if is_source:
                    qlen = 0.0  # Source has no queue
                    util = 0.0  # Source has no utilization
                    respt = 0.0  # Source has no response time
                    residt = 0.0  # Source has no residence time
                    arvr = 0.0  # Source has no arrivals to itself
                    if hasattr(self._sn, 'rates') and self._sn.rates is not None:
                        rates = np.asarray(self._sn.rates)
                        stationToNode = np.asarray(self._sn.stationToNode).flatten()
                        node_idx = int(stationToNode[i])
                        if node_idx < rates.shape[0] and r < rates.shape[1]:
                            tput = rates[node_idx, r]

                # Skip rows where all metrics are zero (matches MATLAB getAvgTable behavior)
                if abs(qlen) < 1e-12 and abs(util) < 1e-12 and abs(tput) < 1e-12:
                    continue

                rows.append({
                    'Station': station_name,
                    'JobClass': class_name,
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        # Wrap in IndexedTable for consistent formatting
        result = IndexedTable(df)

        if not self._table_silent and len(df) > 0:
            print(result)

        return result

    # =========================================================================
    # Individual Metric Accessors
    # =========================================================================

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.Q.copy()

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.U.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.R.copy()

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times (M x K).

        Residence time is computed from response time using visit ratios:
        WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]
        """
        if self._result is None:
            self._ensureAvgResults()

        # Compute ResidT using proper visit ratios from network structure
        if self._sn is not None and self._sn.visits:
            return sn_get_residt_from_respt(self._sn, self._result.R, None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            return self._result.R.copy()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times."""
        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R.copy()
        if hasattr(self._sn, 'rates') and self._sn.rates is not None:
            rates = np.asarray(self._sn.rates)
            S = np.zeros_like(rates)
            nonzero = rates > 0
            S[nonzero] = 1.0 / rates[nonzero]
            W = R - S
            W = np.maximum(W, 0.0)
            return W
        return R

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.T.copy()

    def getStartRate(self) -> np.ndarray:
        """(nstations x nclasses) rate at which a class-r job BEGINS or RESUMES
        holding a server at station i, estimated over the simulated path.

        At a lossless station with no in-service abandonment

            getStartRate == getAvgTput + getPreemptRate

        up to simulation error; SolverCTMC.getStartRate reports the exact value,
        so the two are compared with a two-sample t-test rather than an equality.
        """
        if self._result is None or getattr(self._result, 'startRate', None) is None:
            self._ensureAvgResults()
        rate = getattr(self._result, 'startRate', None)
        if rate is None:
            raise RuntimeError("This solver run produced no START rates.")
        return np.asarray(rate).copy()

    def getPreemptRate(self) -> np.ndarray:
        """(nstations x nclasses) rate at which a class-r job HOLDING A SERVER
        at station i is pushed back into the buffer. Zero at a non-preemptive
        station."""
        if self._result is None or getattr(self._result, 'preemptRate', None) is None:
            self._ensureAvgResults()
        rate = getattr(self._result, 'preemptRate', None)
        if rate is None:
            raise RuntimeError("This solver run produced no PREEMPT rates.")
        return np.asarray(rate).copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates."""
        if self._result is None:
            self._ensureAvgResults()
        if self._result.A is not None:
            return self._result.A.copy()
        # Fallback: compute from throughputs if A not available
        return sn_get_arvr_from_tput(self._sn, self._result.T, None)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times.

        Note:
            For closed networks: uses Little's Law C = N/X
            For open networks: sum of response times across all stations
        """
        if self._result is None:
            self._ensureAvgResults()

        X = self._result.X.flatten()
        njobs = self._sn.njobs.flatten() if self._sn is not None and hasattr(self._sn, 'njobs') else None
        nclasses = len(X)
        C = np.zeros(nclasses)

        for k in range(nclasses):
            if njobs is not None and k < len(njobs) and np.isfinite(njobs[k]):
                # Closed class: use Little's Law (matching MATLAB getAvgSys.m line 135)
                if X[k] > 0:
                    C[k] = njobs[k] / X[k]
                else:
                    C[k] = np.inf
            else:
                # Open class: sum of response times across all stations
                R = self._result.R
                C[k] = np.sum(R[:, k])

        return C

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.X.flatten()

    # =========================================================================
    # SSA-Specific Methods
    # =========================================================================

    def getConfidenceIntervals(self) -> Dict[str, np.ndarray]:
        """
        Get confidence intervals for all metrics.

        Returns:
            Dict with keys 'Q_ci', 'U_ci', 'R_ci', 'T_ci'
        """
        if self._result is None:
            self._ensureAvgResults()

        return {
            'Q_ci': self._result.Q_ci if self._result.Q_ci is not None else np.array([]),
            'U_ci': self._result.U_ci if self._result.U_ci is not None else np.array([]),
            'R_ci': self._result.R_ci if self._result.R_ci is not None else np.array([]),
            'T_ci': self._result.T_ci if self._result.T_ci is not None else np.array([]),
        }

    def getTotalSimulatedTime(self) -> float:
        """Get total simulated time."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.total_time

    def getSampleCount(self) -> int:
        """Get number of samples collected."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.samples

    # =========================================================================
    # Sampling Methods (SSA Supports Sampling)
    # =========================================================================

    def sample(self, node, numEvents: int = None) -> 'SamplePath':
        """
        Generate a sample path with event traces at the given node.

        Runs SSA simulation with event logging and returns a SamplePath
        object with event list matching MATLAB's sampleNodeState format.

        Args:
            node: Node object or node index (1-based)
            numEvents: Number of events to simulate (default: solver's samples option)

        Returns:
            SamplePath with .event list of SampleEvent objects
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import sample_via_jar
            return sample_via_jar(self, node, numEvents)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            return self._sample_via_cpp(node, numEvents)
        # Enable event logging
        self.options.record_events = True
        if numEvents is not None:
            self.options.samples = numEvents

        # Clear previous result to force re-run with event logging
        self._result = None
        self.runAnalyzer()

        # Get station-to-node mapping
        sn = self._sn
        station_to_node = None
        if hasattr(sn, 'stationToNode') and sn.stationToNode is not None:
            station_to_node = sn.stationToNode

        # Build sample path from event log
        sample_path = SamplePath()
        sample_path.handle = node

        if self._result.event_log is not None:
            times = []
            events = []
            for (t, src_st, src_k, dst_st, dst_k) in self._result.event_log:
                # Departure event
                if src_st >= 0:
                    node_idx = int(station_to_node[src_st]) + 1 if station_to_node is not None else src_st + 1
                    events.append(SampleEvent(t=t, node=node_idx, class_idx=src_k + 1, event='DEP'))
                    times.append(t)
                # Arrival event
                if dst_st >= 0:
                    node_idx = int(station_to_node[dst_st]) + 1 if station_to_node is not None else dst_st + 1
                    events.append(SampleEvent(t=t, node=node_idx, class_idx=dst_k + 1, event='ARV'))
                    times.append(t)

            # Derived tags of the sampled path, emitted at the same instant as
            # the transition that carries them, PREEMPT before START (the victim
            # leaves the server before the job that displaced it takes it).
            tag_log = getattr(self._result, 'tag_log', None)
            if tag_log:
                for (t, st, cls_idx, kind) in tag_log:
                    if st < 0:
                        continue
                    node_idx = int(station_to_node[st]) + 1 if station_to_node is not None else st + 1
                    events.append(SampleEvent(t=t, node=node_idx, class_idx=cls_idx + 1,
                                              event='PREEMPT' if kind == 1 else 'START'))
                    times.append(t)
                order = np.argsort(np.asarray(times, dtype=float), kind='stable')
                events = [events[i] for i in order]
                times = [times[i] for i in order]

            sample_path.event = events
            sample_path.t = np.array(times) if times else np.array([])

        if self._result.state_log is not None and len(self._result.state_log) > 0:
            sample_path.state = np.array(self._result.state_log)

        # Disable event logging for subsequent runs
        self.options.record_events = False

        return sample_path

    def _sample_via_cpp(self, node, numEvents) -> 'SamplePath':
        """
        `sample` under lang='cpp': the trajectory from `-s ssa -a sample`.

        THE EVENT COLUMN IS A SYNCHRONIZATION INDEX, not a (node, class, kind)
        triple, so it is expanded here against the model's own sync list -- the
        same list the native path walks. What the engine decided is WHICH
        synchronization fired and WHEN; the descriptor behind the index is
        model structure and carries no simulated quantity.
        """
        from ..cpp_dispatch import sample_path_via_cpp
        from ...lang.sync import refresh_sync

        sn = self._sn if self._sn is not None else self.model.get_struct()
        ind = int(getattr(node, 'index', node))
        nevents = int(numEvents) if numEvents is not None else int(self.options.samples)
        seed = getattr(self.options, 'seed', None)
        p = sample_path_via_cpp(self, nevents, node=ind - 1,
                                seed=int(seed) if seed else None)

        path = SamplePath()
        path.handle = node
        path.t = p['t']
        path.state = p.get('nodeState')
        sync = refresh_sync(sn)
        t = np.asarray(p['t'], dtype=float)
        events = []
        for i, e in enumerate(p['event']):
            # None is the absorbing step: no synchronization fired there, and an
            # index would name one that did not.
            if e is None or e < 0 or e >= len(sync):
                continue
            when = float(t[i]) if i < t.size else float('nan')
            for arm in (sync[e].active, sync[e].passive):
                if arm is None:
                    continue
                events.append(SampleEvent(t=when, node=int(arm.node) + 1,
                                          class_idx=int(arm.job_class) + 1,
                                          event=arm.event.name))
        path.event = events
        path.isaggregate = False
        return path

    def sampleAggr(self, node: int, numEvents: int) -> np.ndarray:
        """Sample aggregated response times."""
        return self.sample(node, numEvents)

    def sampleSys(self, numEvents: int) -> np.ndarray:
        """Sample system-level response times."""
        if self._result is None:
            self._ensureAvgResults()

        mean_cycle_time = np.mean(self._result.C)
        if mean_cycle_time <= 0:
            return np.zeros(numEvents)

        return np.random.exponential(mean_cycle_time, numEvents)

    # =========================================================================
    # CDF and Percentile Methods
    # =========================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """Not available: SolverSSA does not record per-job response times.

        A simulator must report what it measured. The base exponential fit
        carries no information about the tail and would be indistinguishable,
        to the caller, from a measured distribution. SSA samples state
        trajectories, not per-job sojourn times, so there is nothing to build
        an empirical CDF from -- the reference @SolverSSA/getCdfRespT.m refuses
        by name, line-cli refuses -s ssa -a cdf, and so does this port.
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # lang='cpp' cannot serve this getter; the reason is named, not
            # blanket, as in the reference's CPPLINE.cppUnsupported arm
            from ..cpp_dispatch import cpp_unsupported
            cpp_unsupported(
                self, 'getCdfRespT',
                "SolverSSA records no per-job response times on either side, so line-cli "
                "refuses -s ssa -a cdf by name, and so does the reference: "
                "@SolverSSA/getCdfRespT.m raises. Use SolverJMT for a measured CDF, "
                "SolverFLD or SolverMAM for an analytical one")
        raise RuntimeError(
            "SolverSSA does not record per-job response times, so it cannot return an "
            "empirical response time CDF. Use SolverJMT for a measured CDF, or "
            "getPerctRespT(..., 'forktail') for the analytical fork-join tail.")

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None,
        method: Optional[str] = None
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """Extract percentiles from response time distribution.

        SSA records no per-job response times, so the default route -- reading
        getCdfRespT, as the reference's @NetworkSolver/getPerctRespT.m does --
        errors like the reference; only the ForkTail approximation is served.
        This used to fabricate exponential percentiles from the mean, which is
        indistinguishable, to the caller, from a measured tail.
        """
        if method is not None and method.lower() == 'forktail':
            # Fork-join request tail latency; mirrors the MATLAB entry point
            # @NetworkSolver/getPerctRespT.m with method='forktail'
            from ...api.fjnative import forktail_percentiles
            if percentiles is None:
                percentiles = [10, 25, 50, 75, 90, 95, 99]
            return forktail_percentiles(self, percentiles, jobclass)
        raise RuntimeError(
            "Unable to compute percentiles. getCdfRespT not available for this solver: "
            "SolverSSA records no per-job response times. Use SolverJMT for measured "
            "percentiles, or getPerctRespT(..., method='forktail') for the analytical "
            "fork-join tail.")

    # =========================================================================
    # Probability Methods
    # =========================================================================

    def _get_station_index(self, station) -> int:
        """Convert station argument to 0-based index."""
        if isinstance(station, (int, np.integer)):
            return int(station)
        # Assume it's a node/station object
        if hasattr(station, 'get_station_index0'):
            return station.get_station_index0()
        if hasattr(station, '_station_index'):
            return station._station_index
        if hasattr(station, 'station_index'):
            return station.station_index
        raise ValueError(f"Cannot convert {type(station)} to station index")

    def _get_target_state_for_station(self, ist: int) -> Optional[np.ndarray]:
        """Get the target state set for a station via setState()."""
        if self._sn is None:
            return None

        sn = self._sn
        if not hasattr(sn, 'state') or sn.state is None:
            return None

        # Get stateful index for this station
        if hasattr(sn, 'stationToStateful') and sn.stationToStateful is not None:
            stationToStateful = np.asarray(sn.stationToStateful).flatten().astype(int)
            if ist < len(stationToStateful):
                isf = stationToStateful[ist]
            else:
                isf = ist
        else:
            isf = ist

        # Get state for this stateful node
        if isinstance(sn.state, dict):
            # state is a dict keyed by stateful node
            for key, val in sn.state.items():
                idx = key._stateful_index if hasattr(key, '_stateful_index') else -1
                if idx == isf:
                    return np.atleast_1d(val).flatten()
            return None
        elif isinstance(sn.state, (list, np.ndarray)):
            if isf < len(sn.state) and sn.state[isf] is not None:
                return np.atleast_1d(sn.state[isf]).flatten()
        return None

    def _estimate_marginal_prob(self, ist: int, target_state: np.ndarray) -> float:
        """
        Estimate probability of station being in target state using simulation statistics.

        Uses product of Poisson/geometric marginals as approximation.
        """
        if self._result is None:
            self._ensureAvgResults()

        Q = self._result.Q
        nclasses = Q.shape[1]

        # Get mean queue lengths at this station
        mean_q = Q[ist, :nclasses]

        # Target state should have per-class job counts
        target = target_state[:nclasses] if len(target_state) >= nclasses else np.pad(target_state, (0, nclasses - len(target_state)))

        # Use product of Poisson marginals (product-form approximation)
        prob = 1.0
        for k in range(nclasses):
            lam = max(mean_q[k], 1e-10)  # Mean number of class k jobs
            n_k = int(target[k])  # Target number of class k jobs

            # Poisson probability: P(X=n) = (lambda^n * e^-lambda) / n!
            if n_k <= 170:  # Avoid factorial overflow
                prob *= (lam ** n_k) * np.exp(-lam) / math.factorial(n_k)
            else:
                # Use Stirling's approximation for large n
                log_prob = n_k * np.log(lam) - lam - (n_k * np.log(n_k) - n_k + 0.5 * np.log(2 * np.pi * n_k))
                prob *= np.exp(log_prob)

        return float(prob)

    def getProb(self, station=None) -> float:
        """Get probability for the state set on a station.

        Returns the probability that the station is in the state that was
        set via setState(). Uses simulation statistics to estimate probability.

        Args:
            station: Station object or index (0-based).

        Returns:
            Probability (scalar float) for the specified state.
        """
        if self._result is None:
            self._ensureAvgResults()

        if station is None:
            # Return list of probabilities for each station's set state
            probs = []
            for i in range(self._result.Q.shape[0]):
                target = self._get_target_state_for_station(i)
                if target is not None:
                    probs.append(self._estimate_marginal_prob(i, target))
                else:
                    probs.append(0.0)
            return probs

        ist = self._get_station_index(station)
        target = self._get_target_state_for_station(ist)

        if target is None:
            # No state set, return 0
            return 0.0

        return self._estimate_marginal_prob(ist, target)

    def getProbAggr(self, station) -> float:
        """Get aggregated state probability at station.

        Returns the probability that the station is in the aggregated state
        (per-class job counts) that was set via setState().

        Args:
            station: Station object or index (0-based).

        Returns:
            Probability (scalar float) for the specified aggregated state.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-aggr', ist=station, kind='scalar')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import _assert_default_state, prob_aggr_via_cpp
            # `-s ssa -a prob` reports every station's simulated occupancy of the
            # state the model carries, which the wire now carries too; the claim
            # that only `-a avg` exists under `-s ssa` predates that arm.
            _assert_default_state(self, 'getProbAggr')
            p = prob_aggr_via_cpp(self)['probAggr']
            ist0 = station if isinstance(station, (int, np.integer)) else station.get_station_index0()
            if not (0 <= int(ist0) < len(p)):
                raise ValueError("station index %r is outside 0..%d" % (ist0, len(p) - 1))
            return float(p[int(ist0)])
        return self.getProb(station)

    def getProbSys(self) -> float:
        """Get joint system state probability.

        Returns the probability that the entire system is in the state
        that was set via setState() on all stations.
        Uses product of marginal probabilities (product-form approximation).

        Returns:
            Joint probability (scalar float) for the system state.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-sys', kind='scalar')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import _assert_default_state, prob_sys_via_cpp
            _assert_default_state(self, 'getProbSys')
            return prob_sys_via_cpp(self)
        if self._result is None:
            self._ensureAvgResults()

        if self._sn is None:
            return 0.0

        # Compute product of marginal probabilities for each station
        prob = 1.0
        nstations = self._result.Q.shape[0]

        for ist in range(nstations):
            target = self._get_target_state_for_station(ist)
            if target is not None:
                prob *= self._estimate_marginal_prob(ist, target)

        return float(prob)

    def getProbSysAggr(self) -> float:
        """Get system-level aggregated joint probability.

        Returns the joint probability for the aggregated system state
        (per-class job counts at each station).

        Returns:
            Joint probability (scalar float) for the aggregated system state.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-sys-aggr', kind='scalar')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import _assert_default_state, prob_aggr_via_cpp
            _assert_default_state(self, 'getProbSysAggr')
            return float(prob_aggr_via_cpp(self)['probSysAggr'])
        return self.getProbSys()

    def getTranCdfRespT(self, t_max: float = 10.0, n_points: int = 100) -> List[Dict]:
        """Not supported, as in the reference, whose base class raises."""
        raise NotImplementedError("getTranCdfRespT is not supported by SolverSSA")

    def getTranCdfPassT(self, *args, **kwargs) -> Dict:
        """Not supported, as in the reference, whose base class raises."""
        raise NotImplementedError("getTranCdfPassT is not supported by SolverSSA")

    # Aliases for new methods
    GetProb = getProb
    GetProbAggr = getProbAggr
    GetProbSys = getProbSys
    GetProbSysAggr = getProbSysAggr
    GetTranCdfRespT = getTranCdfRespT
    GetTranCdfPassT = getTranCdfPassT

    # =========================================================================
    # UNIFIED METRICS METHOD
    # =========================================================================


    # =========================================================================
    # CHAIN-LEVEL METHODS
    # =========================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self._sn, 'chains') and self._sn.chains is not None:
            chains_arr = np.asarray(self._sn.chains)
            nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else 1

            # Check if chains is 1D (class->chain mapping) or 2D (chain,class membership)
            if chains_arr.ndim == 1:
                # 1D format: chains[k] = c means class k belongs to chain c
                chains = [[] for _ in range(nchains)]
                for k in range(self._sn.nclasses):
                    c = int(chains_arr[k])
                    if 0 <= c < nchains:
                        chains[c].append(k)
                return chains if any(chains) else [[k for k in range(self._sn.nclasses)]]
            else:
                # 2D format: chains[c, k] > 0 means class k is in chain c
                chains = []
                for c in range(nchains):
                    chain_classes = []
                    for k in range(self._sn.nclasses):
                        if chains_arr[c, k] > 0:
                            chain_classes.append(k)
                    chains.append(chain_classes)
                return chains if chains else [[k for k in range(self._sn.nclasses)]]
        else:
            # Default: each class is its own chain
            return [[k] for k in range(self._sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()

        Q = self._result.Q
        chains = self._get_chains()
        nstations = Q.shape[0]
        nchains = len(chains)

        QN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                QN_chain[:, c] = np.sum(Q[:, chain_classes], axis=1)

        return QN_chain

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()

        U = self._result.U
        chains = self._get_chains()
        nstations = U.shape[0]
        nchains = len(chains)

        UN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                UN_chain[:, c] = np.sum(U[:, chain_classes], axis=1)

        return UN_chain

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain.

        Uses alpha-weighted sum matching MATLAB: RN(:,c) = sum(RNclass(:,inchain).*alpha(:,inchain),2)
        """
        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R
        chains = self._get_chains()
        nstations = R.shape[0]
        nchains = len(chains)

        RN_chain = np.zeros((nstations, nchains))

        # Get alpha weights from sn_get_demands_chain
        if self._sn is not None:
            from ...api.sn.demands import sn_get_demands_chain
            try:
                demands = sn_get_demands_chain(self._sn)
                alpha = demands.alpha

                for c, chain_classes in enumerate(chains):
                    if chain_classes:
                        # Weighted sum: sum(R[:, inchain] * alpha[:, inchain], axis=1)
                        RN_chain[:, c] = np.sum(R[:, chain_classes] * alpha[:, chain_classes], axis=1)
                # Clean up tiny values (numerical noise) to exactly 0
                RN_chain[np.abs(RN_chain) < 1e-10] = 0.0
                return RN_chain
            except Exception:
                pass

        # Fallback: use simple mean if alpha computation fails
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                RN_chain[:, c] = np.mean(R[:, chain_classes], axis=1)

        # Clean up tiny values (numerical noise) to exactly 0
        RN_chain[np.abs(RN_chain) < 1e-10] = 0.0
        return RN_chain

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()

        T = self._result.T
        if T.ndim == 1:
            T = T.reshape(1, -1)

        chains = self._get_chains()
        nstations = self._result.Q.shape[0]
        nchains = len(chains)

        TN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                if T.shape[0] == nstations:
                    TN_chain[:, c] = np.sum(T[:, chain_classes], axis=1)
                else:
                    TN_chain[:, c] = np.sum(T[0, chain_classes])

        return TN_chain

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics aggregated by chain.

        Returns:
            Tuple of (QN, UN, RN, WN, AN, TN) aggregated by chain
        """
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

        station_names = self.station_names if hasattr(self, 'station_names') else [f'Station{i}' for i in range(nstations)]

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Chain': f'Chain{c + 1}',  # 1-based to match MATLAB
                    'QLen': QN[i, c],
                    'Util': UN[i, c],
                    'RespT': RN[i, c],
                    'ResidT': WN[i, c],
                    'ArvR': AN[i, c],
                    'Tput': TN[i, c],
                })

        # five SIGNIFICANT digits like MATLAB's table, not pandas' five decimals
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(pd.DataFrame(rows))

    # =========================================================================
    # NODE-LEVEL METHODS
    # =========================================================================

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get average metrics per node.

        Unlike getAvg() which returns station-level metrics, this method
        returns node-level metrics including non-station nodes (e.g., Cache).
        For Cache nodes, hit/miss class throughputs are computed using
        actual hit/miss probabilities.

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        from ...api.sn.getters import sn_get_node_arvr_from_tput, sn_get_node_tput_from_tput

        if self._result is None:
            self._ensureAvgResults()

        sn = self._sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses

        # Get station-level metrics
        QN = self._result.Q
        UN = self._result.U
        RN = self._result.R
        TN = self._result.T
        # Pass None for AN if not available or all zeros, so sn_get_node_arvr_from_tput computes it
        AN = None
        if hasattr(self._result, 'A') and self._result.A is not None:
            if np.any(self._result.A > 0):
                AN = self._result.A

        # Matches MATLAB solver_ssa_analyzer_serial.m:137; store actual
        # probs to preserve across the chain refresh below.
        actual_probs = {}  # {ind: (actualhitprob, actualmissprob)}
        if sn.nodeparam is not None:
            for ind in range(I):
                if sn.nodetype is not None and ind < len(sn.nodetype):
                    if sn.nodetype[ind] == NodeType.CACHE and ind in sn.nodeparam:
                        cache_param = sn.nodeparam[ind]

                        # Check if actualhitprob was already computed by runAnalyzer
                        # If so, preserve it (matches CTMC's guard for pre-computed values)
                        existing_hitprob = getattr(cache_param, 'actualhitprob', None)
                        if existing_hitprob is not None and np.any(existing_hitprob > 0):
                            existing_missprob = getattr(cache_param, 'actualmissprob', np.zeros_like(existing_hitprob))
                            existing_latency = getattr(cache_param, 'actualresidt', None)
                            existing_delayed = getattr(cache_param, 'actualdelayedhitprob', None)
                            # If the simulator post-processed hit/miss but did
                            # not compute latency/phi, derive them from TN/QN now.
                            if existing_latency is None or existing_delayed is None:
                                lat_pair = self._compute_retrieval_latency_ssa(
                                    sn, cache_param, ind, TN, QN)
                                if lat_pair is not None:
                                    if existing_latency is None and lat_pair[0] is not None:
                                        existing_latency = lat_pair[0]
                                        cache_param.actualresidt = existing_latency
                                    if existing_delayed is None and lat_pair[1] is not None:
                                        existing_delayed = lat_pair[1]
                                        cache_param.actualdelayedhitprob = existing_delayed
                            actual_probs[ind] = (existing_hitprob.copy(),
                                                 existing_missprob.copy(),
                                                 existing_latency,
                                                 existing_delayed)
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                                cache_node_obj = self.model._nodes[ind]
                                if hasattr(cache_node_obj, 'set_result_hit_prob'):
                                    cache_node_obj.set_result_hit_prob(existing_hitprob)
                                if hasattr(cache_node_obj, 'set_result_miss_prob'):
                                    cache_node_obj.set_result_miss_prob(existing_missprob)
                                if (existing_latency is not None and
                                        hasattr(cache_node_obj, 'set_result_residt')):
                                    cache_node_obj.set_result_residt(existing_latency)
                                if (existing_delayed is not None and
                                        hasattr(cache_node_obj, 'set_result_delayed_hit_prob')):
                                    cache_node_obj.set_result_delayed_hit_prob(existing_delayed)
                            continue

                        hitclass = getattr(cache_param, 'hitclass', np.array([]))
                        missclass = getattr(cache_param, 'missclass', np.array([]))
                        nclasses = len(hitclass) if hasattr(hitclass, '__len__') else R

                        # Get the stateful index for this cache node
                        isf = int(sn.nodeToStateful[ind]) if sn.nodeToStateful is not None and ind < len(sn.nodeToStateful) else -1
                        # Get the station index for throughput lookup
                        ist = int(sn.nodeToStation[ind]) if hasattr(sn, 'nodeToStation') and sn.nodeToStation is not None and ind < len(sn.nodeToStation) else -1

                        cache_param.actualhitprob = np.zeros(nclasses)
                        cache_param.actualmissprob = np.zeros(nclasses)

                        for k in range(nclasses):
                            h = int(hitclass[k]) if k < len(hitclass) else -1
                            m = int(missclass[k]) if k < len(missclass) else -1
                            if h >= 0 and m >= 0:
                                # Compute actual hit/miss probs from simulation throughputs
                                # Use throughput at this station for hit and miss classes
                                if ist >= 0 and ist < TN.shape[0] and h < TN.shape[1] and m < TN.shape[1]:
                                    t_hit = TN[ist, h]
                                    t_miss = TN[ist, m]
                                    t_total = t_hit + t_miss
                                    if t_total > 0:
                                        cache_param.actualhitprob[k] = t_hit / t_total
                                        cache_param.actualmissprob[k] = t_miss / t_total
                                    else:
                                        # Fallback to NaN if no throughput (matches MATLAB line 134)
                                        cache_param.actualhitprob[k] = np.nan
                                        cache_param.actualmissprob[k] = np.nan

                        # see _kb/06-solver-catalog.md (SSA: "Warmup discard,
                        # and cache hit/miss accounting") -- python-only Eq. 8/4
                        latencies, delayed = self._compute_retrieval_latency_ssa(
                            sn, cache_param, ind, TN, QN)
                        if latencies is not None:
                            cache_param.actualresidt = latencies
                        if delayed is not None:
                            cache_param.actualdelayedhitprob = delayed

                        # Store for restoration after chain refresh
                        actual_probs[ind] = (cache_param.actualhitprob.copy(),
                                             cache_param.actualmissprob.copy(),
                                             latencies,
                                             delayed)

                        # Set result on Cache node in model
                        if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                            cache_node = self.model._nodes[ind]
                            if hasattr(cache_node, 'set_result_hit_prob'):
                                cache_node.set_result_hit_prob(cache_param.actualhitprob)
                            if hasattr(cache_node, 'set_result_miss_prob'):
                                cache_node.set_result_miss_prob(cache_param.actualmissprob)
                            if (latencies is not None and
                                    hasattr(cache_node, 'set_result_residt')):
                                cache_node.set_result_residt(latencies)
                            if (delayed is not None and
                                    hasattr(cache_node, 'set_result_delayed_hit_prob')):
                                cache_node.set_result_delayed_hit_prob(delayed)

            # After setting actual hit/miss probs, refresh chains to recompute visits
            # This matches MATLAB's behavior in SSA.runAnalyzer (lines 108-111)
            if hasattr(self, 'model') and hasattr(self.model, '_refresh_chains'):
                self.model._refresh_chains()
                # Update sn with refreshed structure
                if hasattr(self.model, '_sn'):
                    sn = self.model._sn
                    self._sn = sn
                    # Restore actual hit/miss probs that were computed from simulation
                    if sn.nodeparam is not None:
                        for ind, entry in actual_probs.items():
                            if ind not in sn.nodeparam:
                                continue
                            # Tuple is (ahp, amp, [latency, [delayed]])
                            ahp = entry[0]
                            amp = entry[1]
                            lat = entry[2] if len(entry) > 2 else None
                            dly = entry[3] if len(entry) > 3 else None
                            sn.nodeparam[ind].actualhitprob = ahp
                            sn.nodeparam[ind].actualmissprob = amp
                            if lat is not None:
                                sn.nodeparam[ind].actualresidt = lat
                            if dly is not None:
                                sn.nodeparam[ind].actualdelayedhitprob = dly

        # Compute station-level ResidT from RespT using visit ratios
        # This matches MATLAB's getAvgNode which calls getAvg first (getAvg computes WN)
        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(sn, RN, None)

        # Create TH (throughput handle)
        TH = np.zeros_like(TN)
        TH[TN > 0] = 1.0

        # Compute node arrival rates and throughputs
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN)
        TNn = sn_get_node_tput_from_tput(sn, TN, TH, ANn)

        # Initialize other node-level metrics
        QNn = np.zeros((I, R))
        UNn = np.zeros((I, R))
        RNn = np.zeros((I, R))
        WNn = np.zeros((I, R))

        # Copy station metrics to station nodes
        for ist in range(M):
            ind = sn.stationToNode[ist]
            if ind >= 0 and ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]

        return QNn, UNn, RNn, WNn, ANn, TNn

    def getAvgNodeTable(self) -> pd.DataFrame:
        """
        Get average metrics by node as DataFrame.

        Returns node-based results (one row per node per class) including
        non-station nodes like Cache. For Cache nodes, hit/miss class
        throughputs are computed using actual hit/miss probabilities.

        Returns:
            pandas.DataFrame with columns: Node, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()

        sn = self._sn
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []

        from ..cache_table import retrieval_hidden_classes
        hidden = retrieval_hidden_classes(sn)

        rows = []
        for node_idx in range(sn.nnodes):
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

            for r in range(sn.nclasses):
                if r in hidden:
                    continue  # auxiliary retrieval class - omit from node table
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                # Filter out all-zero rows
                if abs(QNn[node_idx, r]) < 1e-10 and abs(UNn[node_idx, r]) < 1e-10 and \
                   abs(RNn[node_idx, r]) < 1e-10 and abs(ANn[node_idx, r]) < 1e-10 and abs(TNn[node_idx, r]) < 1e-10:
                    continue

                rows.append({
                    'Node': node_name,
                    'JobClass': class_name,
                    'QLen': QNn[node_idx, r],
                    'Util': UNn[node_idx, r],
                    'RespT': RNn[node_idx, r],
                    'ResidT': WNn[node_idx, r],
                    'ArvR': ANn[node_idx, r],
                    'Tput': TNn[node_idx, r],
                })

        df = pd.DataFrame(rows)

        if not self._table_silent:
            print(df.to_string(index=False))

        return df

    def getAvgCacheTable(self) -> pd.DataFrame:
        """Detailed per-class cache performance metrics (see cache_table)."""
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import cache_table_via_jar
            return cache_table_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import cache_table_via_cpp
            return cache_table_via_cpp(self)
        from ..cache_table import build_cache_avg_table
        return build_cache_avg_table(self)

    get_avg_cache_table = getAvgCacheTable
    avg_cache_table = getAvgCacheTable

    def getAvgItemTable(self) -> pd.DataFrame:
        """Item-level cache occupancy table (see cache_table)."""
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import item_table_via_jar
            return item_table_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import item_table_via_cpp
            return item_table_via_cpp(self)
        from ..cache_table import build_item_avg_table
        return build_item_avg_table(self)

    get_avg_item_table = getAvgItemTable
    avg_item_table = getAvgItemTable

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

    def getTranAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get transient average metrics from simulation.

        SSA provides transient metrics from the simulation trajectory.

        Returns:
            Tuple of (Q, U, T) transient queue lengths, utilizations, and throughputs
        """
        if self._result is None:
            self._ensureAvgResults()

        Q = self._result.Q
        U = self._result.U
        T = self._result.T if hasattr(self._result, 'T') else np.zeros_like(Q)

        return Q, U, T

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T) where R is system response time and T is system throughput
        """
        R = self.getAvgSysRespT()
        T = self.getAvgSysTput()
        return R, T

    getAvgSysTable = NetworkSolver.getAvgSysTable  # chain-level shared layout

    # =========================================================================
    # Additional Sampling Methods
    # =========================================================================

    def sampleSysAggr(self, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated system response times."""
        return self.sampleSys(numEvents)

    # =========================================================================
    # Introspection Methods
    # =========================================================================

    def listValidMethods(self) -> List[str]:
        """List valid solution methods.

        Returns:
            List of valid method names:
            - 'default', 'serial', 'ssa': Serial Gillespie simulation
            - 'nrm': Next Reaction Method
            - 'parallel', 'para': Parallel simulation with multiple replicas
        """
        return ['default', 'serial', 'ssa', 'nrm', 'parallel', 'para']

    def isStochasticMethod(self, method):
        """All SSA methods are stochastic simulation."""
        return True

    is_stochastic_method = isStochasticMethod

    def supportsModelMethod(self, method):
        """The fork-join model class, which EVERY SSA method has to clear.

        runAnalyzer tag-augments a fork-join model through
        ``ModelAdapter.fjtag``, whose first act is ``sn_fj_validate``, so a model
        that validator refuses is refused whichever method was asked for. The
        featset cannot state it -- Fork and Join are declared, and the rules are
        about how they are WIRED (the pairing, the join strategy, the tasks per
        link, whether an open class is routed through the fork) -- so it is
        structural, and it is the validator's own body of rules rather than a
        copy of them.

        Without it the report offered every ssa.* row on a fork-join model whose
        Join names no fork, and each one then raised; SolverCTMC gates on the
        same predicate for the same reason. Mirrors MATLAB
        @SolverSSA/supportsModelMethod.
        """
        ok, reason = super().supportsModelMethod(method)
        model = getattr(self, 'model', None)
        if ok and model is not None and hasattr(model, 'getStruct'):
            from ...api.fjnative import sn_fj_supports
            ok, reason = sn_fj_supports(model.getStruct())
        # 'nrm' TAKES NO GATE IN THIS PORT, and that is a statement about this
        # engine rather than an omission. MATLAB refuses an explicit 'nrm' only
        # on a scheduling policy the reaction network has no form for
        # (solver_ssa_analyzer_nrm.m:24), the JAR on that plus phase-type service
        # at a preemptive station, and C++ on all six of its checks because its
        # NRM has no fallback arm at all. This port's handler falls back to the
        # serial engine for EVERY one of the six (api/solvers/ssa/handler.py,
        # method == 'nrm'), so an explicit request always returns an answer and
        # withdrawing the row here would refuse a run that succeeds.
        # ``nrm_supports`` states the same six conditions for a caller that wants
        # to know which engine will actually run.
        return ok, reason

    supports_model_method = supportsModelMethod

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features.

        The native SSA solver now drives the same ``State.afterEvent`` machinery
        as the CTMC solver, so its capability set mirrors CTMC's (plus PAS).
        """
        return {
            'Source', 'Sink',
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue', 'Router',
            'MAP', 'APH', 'MMPP2', 'MMAP', 'PH', 'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp',
            'Det', 'Gamma', 'Weibull', 'Lognormal', 'Pareto', 'Uniform',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer', 'Buffer', 'Dispatcher',
            # see _kb/06-solver-catalog.md (SSA: "the NRM engine now supports
            # finite capacity regions directly")
            'Region',
            'Cache', 'CacheClassSwitcher', 'CacheRetrieval',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS',
            'SchedStrategy_DPS', 'SchedStrategy_GPS',
            'SchedStrategy_SIRO', 'SchedStrategy_SEPT',
            'SchedStrategy_LEPT', 'SchedStrategy_FCFS',
            'SchedStrategy_HOL', 'SchedStrategy_LCFS',
            'SchedStrategy_LCFSPR', 'SchedStrategy_LCFSPRPRIO', 'SchedStrategy_FCFSPRPRIO',
            'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO',
            'SchedStrategy_PAS', 'SchedStrategy_OI', 'SchedStrategy_LPS', 'SchedStrategy_EXT',
            'SchedStrategy_POLLING',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ', 'RoutingStrategy_SQ', 'RoutingStrategy_SDR',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',
            'ReplacementStrategy_HLRU', 'ReplacementStrategy_CLIMB', 'ReplacementStrategy_QLRU',
            'ClosedClass', 'SelfLoopingClass', 'OpenClass', 'Replayer',
            'OpenSignal', 'ClosedSignal',
            'SignalType_NEGATIVE', 'SignalType_CATASTROPHE',
            'SignalBatchRemoval', 'SignalRemovalPolicy',
            'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',
            # Fork-join, on the TAG-AUGMENTED copy runAnalyzer builds through
            # ModelAdapter.fjtag, exactly as SolverCTMC does. MATLAB, the JAR and
            # C++ have always declared these four; this port alone omitted them,
            # so the ssa family never appeared on a fork-join model it solves --
            # measured against the exact chain on a symmetric closed fork-join.
            # The wiring rules the names cannot state (the pairing, the join
            # strategy, an open class through the fork) are gated structurally in
            # supportsModelMethod, against sn_fj_validate's own body of rules.
            'Fork', 'Join', 'Forker', 'Joiner',
            'Balking', 'Reneging', 'Retrial',
            'LoadDependence',
            'ClassDependence',
            'JointDependence',
            'GlobalDependence',
            # c-server stations and binding buffers are both State constructs
            # the sampler carries directly: a job that finds no room blocks in a
            # closed model and is dropped in an open one.
            'MultiServer', 'FiniteCapacity',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Mirrors MATLAB SolverSSA.supports: gates the model's used language
        features against getFeatureSet(). This previously checked only the
        station and class counts, so it accepted every model regardless of the
        features it used. Struct-like inputs without a feature registry fall
        back to a structural sanity check.
        """
        if hasattr(model, 'get_used_lang_features') or hasattr(model, 'getUsedLangFeatures'):
            from ..base import SolverFeatureSet
            feat_used = (model.get_used_lang_features()
                         if hasattr(model, 'get_used_lang_features')
                         else model.getUsedLangFeatures())
            feat_supported = SolverFeatureSet()
            feat_supported.set_true(list(SolverSSA.getFeatureSet()))
            return SolverFeatureSet.supports(feat_supported, feat_used)
        try:
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

            return nstations > 0 and nclasses > 0
        except Exception:
            return False

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'default',
            'tol': 1e-4,
            'samples': 10000,
            'seed': 0,
            'cutoff': float('inf'),
            'confidence_level': 0.95,
            'verbose': default_verbose(),
        })

    GetAvg = NetworkSolver.getAvg
    GetAvgTable = getAvgTable
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    ListValidMethods = listValidMethods
    GetFeatureSet = getFeatureSet
    Supports = supports
    DefaultOptions = defaultOptions
    default_options = defaultOptions

    # Chain-level aliases
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
    GetTranAvg = getTranAvg

    # Sampling aliases
    SampleSysAggr = sampleSysAggr

    # Short aliases (MATLAB compatibility)
    aT = getAvgTable
    aNT = getAvgNodeTable
    aCT = getAvgChainTable
    aNCT = getAvgNodeChainTable
    aST = getAvgSysTable
    avgT = getAvgTable
    nodeAvgT = getAvgNodeTable
    chainAvgT = getAvgChainTable
    nodeChainAvgT = getAvgNodeChainTable
    sysAvgT = getAvgSysTable
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    avg_node_table = getAvgNodeTable
    avg_node_chain_table = getAvgNodeChainTable


__all__ = ['SolverSSA', 'SolverSSAOptions']
