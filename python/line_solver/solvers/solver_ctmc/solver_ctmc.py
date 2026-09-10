"""
Native Python implementation of CTMC (Continuous-Time Markov Chain) solver.

This implementation uses pure Python/NumPy algorithms from the api.solvers.ctmc
module.
"""

import warnings
import os
import numpy as np
import pandas as pd
import sys
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass, field
from ...constants import default_verbose

from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn.getters import sn_get_node_tput_from_tput, sn_get_node_arvr_from_tput, sn_get_arvr_from_tput
from ...api.fjnative import sn_fj_supports
from ..fjtag_transform import FJTagTransformMixin
from ..transform_driver import TransformSolveMixin
from ...api.sn.network_struct import NodeType
from ...api.io.logging import line_debug, line_warning
from ...constants import GlobalConstants
from ..base import NetworkSolver, method_type


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
class EventInfo:
    """Information about a single event in a simulation trace."""
    node: int = 0
    jobclass: int = 0
    t: float = 0.0
    event: str = None


@dataclass
class SampleResult:
    """Container for sample-based simulation results."""
    handle: str = ""
    t: np.ndarray = None
    state: np.ndarray = None
    event: List[EventInfo] = None
    isaggregate: bool = False
    nodeIndex: int = None
    numEvents: int = 0

    def __post_init__(self):
        if self.event is None:
            self.event = []


@dataclass
class SolverCTMCOptions:
    """Options for the native CTMC solver."""
    method: str = 'default'
    tol: float = 1e-4
    cutoff: int = 10
    seed: int = 23000
    samples: int = 10000  # Number of samples for simulation-based methods
    verbose: bool = field(default_factory=default_verbose)
    keep: bool = True  # Whether to keep state space after analysis
    force: bool = False  # Force solver to run even if state space may be too large
    config: Dict[str, Any] = field(default_factory=dict)  # Configuration dict (e.g., {'nonmkv': 'none'})
    init_sol: Optional[np.ndarray] = None  # Chain-mode initial distribution (transient analysis, sample paths)
    timespan: Optional[List[float]] = None  # Time interval [t_start, t_end] for transient analysis
    timestep: Optional[float] = None  # Time step for transient analysis (None = auto, matches MATLAB [])
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)
    gen_method: str = 'default'  # 'default' = monolithic builder, 'sync' = sync-action-based builder
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))  # env LINE_SOLVER_LANG overrides; 'python' (native), 'java' (jline.jar via JSON) or 'cpp' (line-cli via JSON)
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. Meaningless for the other langs, which are IEEE double
    # throughout, so line-cli is invoked without --arith unless the caller sets it.
    # Every step from the generator to the means is a field operation, so 'exact'
    # returns the exact rational stationary law here.
    arith: Optional[str] = None


class _QRFResult:
    """Lightweight result container for QRF approximation methods."""

    def __init__(self, QN, UN, RN, TN, CN, XN, runtime, method):
        self.Q = QN
        self.U = UN
        self.R = RN
        self.T = TN
        self.C = CN
        self.X = XN
        self.runtime = runtime
        self.method = method
        self.pi = None
        self.depRates = None


class _CFTPResult:
    """Result container for the perfect-sampling (cftp) method.

    Carries the sampled states alongside the metrics: they are the only
    representation of the stationary distribution this method produces, since
    no state space is enumerated.
    """

    def __init__(self, QN, UN, RN, TN, CN, XN, runtime, method, samples, horizon, pAggr, SSq):
        self.Q = QN
        self.U = UN
        self.R = RN
        self.T = TN
        self.C = CN
        self.X = XN
        self.runtime = runtime
        self.method = method
        self.pi = pAggr
        self.space = SSq
        self.spaceAggr = SSq
        self.cftpSamples = samples
        self.cftpHorizon = horizon
        self.depRates = None


class _MDDResult:
    """Result container for the decision-diagram aggregation (mdd) method.

    No state space is enumerated, so pi/space stay empty; the diagram and the
    level sizes are carried instead, as the only description of how the
    reachable set was represented.
    """

    def __init__(self, QN, UN, RN, TN, CN, XN, runtime, method, mddinfo):
        self.Q = QN
        self.U = UN
        self.R = RN
        self.T = TN
        self.C = CN
        self.X = XN
        self.runtime = runtime
        self.method = method
        self.pi = None
        self.mdd = mddinfo
        self.depRates = None


def _eventRates(F):
    """Minimum positive rate of each event filtration, and the filtration
    normalized by it. An event with no positive rate contributes neither.

    Returns:
        (rates, shapes) with rates a float array over the events and shapes a
        list holding the normalized filtration, or None for an inactive event.
    """
    rates = np.zeros(len(F))
    shapes = [None] * len(F)
    for e in range(len(F)):
        Fe = F[e]
        Fe = np.asarray(Fe.todense() if hasattr(Fe, 'todense') else Fe, dtype=np.float64)
        pos = Fe[Fe > 0]
        if pos.size == 0:
            continue
        rates[e] = pos.min()
        shapes[e] = Fe / rates[e]
    return rates, shapes


class SolverCTMC(FJTagTransformMixin, TransformSolveMixin, NetworkSolver):
    """
    Native Python CTMC (Continuous-Time Markov Chain) solver.

    This solver analyzes queueing networks through exact state-space enumeration
    using pure Python/NumPy, providing the same functionality as the Java wrapper
    without requiring the JVM.

    Supported methods:
        - 'default': Basic state-space enumeration
        - 'gpu': the gpuArray backend of ctmc_solve, which falls back to the
          plain direct solve when no GPU is present

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default')
        **kwargs: Additional solver options
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._sn = None

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
            elif 'verbose' in method_or_options:
                kwargs.setdefault('verbose', method_or_options['verbose'])
            if hasattr(method_or_options, 'cutoff'):
                kwargs.setdefault('cutoff', method_or_options.cutoff)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if hasattr(method_or_options, 'force'):
                kwargs.setdefault('force', method_or_options.force)
            if 'timespan' in method_or_options:
                kwargs.setdefault('timespan', method_or_options['timespan'])
            if 'timestep' in method_or_options:
                kwargs.setdefault('timestep', method_or_options['timestep'])
            # The config map is carried too: a key dropped here is dropped
            # SILENTLY, so the solver would answer the default while the caller
            # believes it asked for something else.
            if 'config' in method_or_options and method_or_options['config']:
                kwargs.setdefault('config', dict(method_or_options['config']))
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'cutoff'):
                kwargs.setdefault('cutoff', method_or_options.cutoff)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if hasattr(method_or_options, 'force'):
                kwargs.setdefault('force', method_or_options.force)
            if hasattr(method_or_options, 'timespan'):
                kwargs.setdefault('timespan', method_or_options.timespan)
            if hasattr(method_or_options, 'timestep'):
                kwargs.setdefault('timestep', method_or_options.timestep)
            # See the note on the dict-like branch above: a config key dropped
            # here is dropped silently.
            if getattr(method_or_options, 'config', None):
                kwargs.setdefault('config', dict(method_or_options.config))
        else:
            self.method = 'default'

        # A method= keyword is the native-Python call style; honour it when no
        # positional method was given instead of dropping it, which silently
        # solved with 'default' whatever the caller asked for.
        method_kw = kwargs.pop('method', None)
        if method_or_options is None and method_kw is not None:
            self.method = str(method_kw).lower()
        self.options = SolverCTMCOptions(method=self.method, **kwargs)

        # Chain mode: a user-supplied MarkovProcess (CTMC) or MarkovChain (DTMC)
        # is solved directly, so there is no network structure to extract. A DTMC
        # is carried as its P-I image, which has the same stationary vector.
        from ...lang.processes import MarkovChain, MarkovProcess
        self._chain_matrix = model if isinstance(model, MarkovChain) else None
        if isinstance(model, MarkovProcess):
            self._chain_process = model
        elif self._chain_matrix is not None:
            self._chain_process = model.toCTMC()
        else:
            self._chain_process = None
        if self._chain_process is not None:
            return

        # Extract network structure
        self._extract_network_params()

    def isChainSolver(self) -> bool:
        """True when the solver was built from a MarkovProcess or a MarkovChain."""
        return getattr(self, '_chain_process', None) is not None

    is_chain_solver = isChainSolver

    def isDiscreteChain(self) -> bool:
        """True in chain mode when the user supplied a DTMC (MarkovChain)."""
        return getattr(self, '_chain_matrix', None) is not None

    is_discrete_chain = isDiscreteChain

    def getTransMat(self) -> np.ndarray:
        """Transition matrix of the user-supplied DTMC (chain mode only)."""
        if not self.isDiscreteChain():
            raise RuntimeError("getTransMat requires a SolverCTMC built from a MarkovChain.")
        return np.asarray(self._chain_matrix.getTransMat(), dtype=np.float64)

    get_trans_mat = getTransMat

    def _assert_not_chain_model(self, caller: str) -> None:
        """Guard for the entry points that need stations and classes."""
        if self.isChainSolver():
            kind = 'MarkovChain' if self.isDiscreteChain() else 'MarkovProcess'
            raise RuntimeError(
                f"{caller} requires a Network model. This solver was built from a {kind}, "
                "which has no stations or classes: use getProbSys, getGenerator, getStateSpace, "
                "getTranProbSys or sampleSys instead.")

    def _chain_state_space(self) -> np.ndarray:
        """State space of the user-supplied chain, or the state indices."""
        space = self._chain_matrix.stateSpace if self.isDiscreteChain() else self._chain_process.stateSpace
        if space is not None and np.asarray(space).size > 0:
            return np.atleast_2d(np.asarray(space, dtype=np.float64))
        n = np.asarray(self._chain_process.getGenerator(), dtype=np.float64).shape[0]
        return np.arange(1, n + 1, dtype=np.float64).reshape(n, 1)

    def _chain_run_analyzer(self) -> None:
        """Steady-state analysis of the user-supplied chain."""
        import time
        from ...api.mc import ctmc_solve, ctmc_solve_reducible, dtmc_solve, dtmc_solve_reducible
        from ...api.solvers.ctmc.analyzers import CTMCResult

        start_time = time.time()
        infgen = np.asarray(self._chain_process.getGenerator(), dtype=np.float64)
        n = infgen.shape[0]
        if self.isDiscreteChain():
            P = self.getTransMat()
            pi = np.asarray(dtmc_solve(P), dtype=np.float64).flatten()
            if not self._is_chain_distribution(pi, n):
                pi = np.asarray(dtmc_solve_reducible(P), dtype=np.float64).flatten()
        else:
            pi = np.asarray(ctmc_solve(infgen), dtype=np.float64).flatten()
            if not self._is_chain_distribution(pi, n):
                pi = np.asarray(ctmc_solve_reducible(infgen), dtype=np.float64).flatten()

        result = CTMCResult()
        result.pi = pi
        result.infgen = infgen
        result.space = self._chain_state_space()
        result.runtime = time.time() - start_time
        result.method = self.method
        self._result = result

    @staticmethod
    def _is_chain_distribution(pi: np.ndarray, n: int) -> bool:
        """Reject a solution the primary solver could not produce on a reducible chain."""
        pi = np.asarray(pi, dtype=np.float64).flatten()
        if pi.size != n or not np.all(np.isfinite(pi)):
            return False
        return bool(np.all(pi >= -1e-8) and abs(pi.sum() - 1) <= 1e-4)

    def _chain_init_distribution(self) -> np.ndarray:
        """options.init_sol when it matches the chain size, uniform otherwise."""
        n = np.asarray(self._chain_process.getGenerator(), dtype=np.float64).shape[0]
        pi0 = getattr(self.options, 'init_sol', None)
        if pi0 is not None and np.asarray(pi0).size == n:
            pi0 = np.asarray(pi0, dtype=np.float64).flatten()
            return pi0 / pi0.sum()
        return np.ones(n) / n

    def _network_init_distribution(self) -> np.ndarray:
        """
        pi(0) for a Network model's transient analyses: the model's INITIAL STATE.

        `pi0(matchrow(stateSpace, s0)) = 1`, which is what
        `@SolverCTMC/getTranProbSys.m` does. The transient getters used to seed
        e_0 instead, on the assumption that row 0 of the enumerated space is the
        initial state; it is not. On a two-station closed model with 2 jobs the
        space begins at (Think 0, Q1 2) while the model starts at (Think 2, Q1 0),
        so pi(2) came out [0.2122, 0.4000, 0.3878] against MATLAB's [0.1939,
        0.4000, 0.4061] -- a wrong answer to the right question, with nothing in
        the output to say which state it started from.

        A state the enumeration does not contain is an error rather than a
        fallback: the alternative is answering for a state the model is not in.
        """
        space = self._result.space
        if space is None or np.asarray(space).size == 0:
            raise RuntimeError(
                "the transient analysis needs the enumerated state space to place pi(0) and the "
                "solve returned none")
        space = np.atleast_2d(np.asarray(space, dtype=float))
        pi0 = np.zeros(space.shape[0])
        init = getattr(self.options, 'init_sol', None)
        if init is not None and np.asarray(init).size == space.shape[0]:
            init = np.asarray(init, dtype=float).reshape(-1)
            return init / init.sum()

        sn = self._sn if self._sn is not None else self.model.getStruct()
        rows = []
        for st in list(sn.state or []):
            rows.extend(np.asarray(st, dtype=float).reshape(-1).tolist())
        s0 = np.asarray(rows, dtype=float)
        if s0.size == space.shape[1]:
            hit = np.flatnonzero(np.all(np.isclose(space, s0), axis=1))
            if hit.size:
                pi0[hit[0]] = 1.0
                return pi0
        raise RuntimeError(
            "the model's initial state %s is not a row of the enumerated state space (%d x %d), "
            "so pi(0) cannot be placed; set options.init_sol to name the distribution explicitly"
            % (np.array2string(s0, precision=6), space.shape[0], space.shape[1]))

    def reset(self):
        """Clear cached results so the solver re-runs on next query."""
        self._clearResultStores()
        self._sn = None
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "CTMC"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model for CTMC computation."""
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

        # native CTMC solver does not accept JAR-wrapper models (no wrapper_sn_to_native bridge), keeping it JVM-free.
        if hasattr(model, 'get_struct'):
            self._sn = model.get_struct()
            if self._sn is not None:
                return

        # Priority 4: Already a native NetworkStruct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            self._sn = model
            return

        raise ValueError(
            "Cannot extract a native NetworkStruct from model. The native CTMC "
            "solver accepts only native Network / NetworkStruct inputs (no JAR "
            "wrapper models).")

    def _reject_unbounded_open_spn(self):
        """Reject a genuinely unbounded open stochastic Petri net.

        SolverCTMC solves the BOUNDED open-SPN case (a Source feeds a Place drained
        by a Transition). The marking must be bounded by a finite per-class Place
        capacity (setClassCapacity) or by the solver cutoff; a finite
        Source->Place->Transition->Sink net then reproduces M/M/1/K. Only a net
        where every Source-fed Place has infinite capacity AND no finite cutoff is
        genuinely unbounded and cannot be built as a finite generator. Runs for both
        the native and lang=java backends, so the JAR does not silently return a
        cutoff-truncated answer for an unbounded net.
        """
        from ...lang.base import NodeType
        sn = getattr(self, '_sn', None)
        if sn is None or getattr(sn, 'nodetype', None) is None:
            return
        _nt = np.asarray(sn.nodetype).ravel()
        _has_spn = bool(np.any(_nt == NodeType.PLACE))
        _has_open = bool(getattr(sn, 'njobs', None) is not None
                         and np.any(np.isinf(np.asarray(sn.njobs, dtype=float))))
        if not (_has_spn and _has_open):
            return
        _cutoff = getattr(self.options, 'cutoff', None)
        _cutoff_finite = (_cutoff is not None
                          and np.all(np.isfinite(np.asarray(_cutoff, dtype=float)))
                          and np.all(np.asarray(_cutoff, dtype=float) > 0))
        _cap_finite = False
        _cc = getattr(sn, 'classcap', None)
        if _cc is not None:
            _cc = np.asarray(_cc, dtype=float)
            for _pi in np.where(_nt == NodeType.PLACE)[0]:
                _st = int(sn.nodeToStation[_pi]) if getattr(sn, 'nodeToStation', None) is not None else -1
                if 0 <= _st < _cc.shape[0] and np.any(np.isfinite(_cc[_st, :])):
                    _cap_finite = True
                    break
        if not _cutoff_finite and not _cap_finite:
            raise RuntimeError(
                "This open stochastic Petri net is unbounded: every "
                "Source-fed Place has infinite capacity and no finite CTMC "
                "cutoff is set, so an infinite generator cannot be built. "
                "Bound the marking with a finite Place capacity "
                "(Place.setClassCapacity) or a finite SolverCTMC cutoff, or "
                "use SolverJMT for the unbounded net.")

    def supportsTransientAnalysis(self):
        """Transient averages are available (uniformization of the generator over options.timespan)."""
        return True

    supports_transient_analysis = supportsTransientAnalysis

    def _ensureAvgResults(self):
        """Chain mode has no averages to gate, only the stationary vector."""
        if self.isChainSolver():
            if self._result is None:
                self._chain_run_analyzer()
            return
        super()._ensureAvgResults()

    def _run_chain_aggregation(self, sn):
        """Solve the CHAIN-AGGREGATED model and map its metrics back to the classes.

        ModelAdapter.aggregate_chains collapses every chain onto a single class,
        class switching disappearing with it, and sn_deaggregate_chain_results
        maps chain-level metrics back through alpha, the per-station share of the
        chain's visits each class carries. What is traded is exactness on a
        non-product-form model: one aggregate service law replaces the per-class
        ones. A caller who needs the exact multiclass answer leaves the flag off
        and pays the state space.
        """
        # Driven by TransformSolveMixin, so the aggregate is solved by an
        # instance of THIS solver rather than a hard-wired SolverCTMC. Clearing
        # the flag states that the aggregate must not be re-aggregated, rather
        # than relying on its nchains == nclasses guard to decline.
        cfg = dict(self.options.config or {})
        cfg['chain_aggregation'] = False
        cfg['transform'] = 'chains'
        self.options.config = cfg
        self._run_transform(sn, 'chainaggr')

    def _transform_publish(self, tr, method):
        """SolverCTMC keeps a _QRFResult, not the dict the mixin defaults to."""
        self._result = _QRFResult(tr.Q, tr.U, tr.R, tr.T, tr.C, tr.X, tr.runtime, method)
        # The sweep count is a reported property of an ITERATED strategy, not a
        # diagnostic: a Jacobi coupling reaches the same fixed point at a
        # different count, so the count is what pins the four codebases.
        self._result.iter = tr.iter
        self._extract_names()
        return self._result

    def _run_transform(self, sn, label=None):
        """Run whichever transformation options.config['transform'] names.

        The strategy rewrites the model into subproblems, TransformSolveMixin
        solves each with an instance of THIS solver, and the strategy maps the
        metrics back. LABEL overrides the suffix of the reported method name so
        the older `chain_aggregation` entry keeps reporting `/chainaggr` and
        stays in step with the MATLAB, JAR and C++ twins.
        """
        import sys

        from ..base import print_solver_banner

        tr = self.transform_solve(sn)
        runtime = tr.runtime
        method = str(self.options.method) + '/' + (label if label else tr.method)
        self._transform_publish(tr, method)
        if self.options.verbose:
            py_version = "%d.%d.%d" % (sys.version_info.major, sys.version_info.minor,
                                       sys.version_info.micro)
            print_solver_banner(
                "CTMC analysis [method: %s; type: %s; lang: python; env: %s] completed in %.6fs."
                % (method, method_type('CTMC', self.options.method), py_version, runtime))

    def _run_fes_aggregation(self, sn):
        """Solve with a station subset replaced by a FLOW-EQUIVALENT SERVER.

        ModelAdapter.aggregate_fes has existed in all four codebases with no
        solver consumer at all: it was exercised by examples and tests only, so
        nothing in the solver stack depended on it. Flow-equivalent aggregation
        is the standard route to HIERARCHICAL DECOMPOSITION -- a subnetwork is
        solved in isolation and enters the outer chain as a single
        load-dependent station, which is what makes an otherwise intractable
        state space tractable. This is that consumer.

        The reduced model answers for the surviving stations directly. For a
        collapsed station the answer is the Chandy-Herzog-Woo conditional sum
        E[Q_i] = sum_n P(N_fes = n) * Q_i(n), with P read off the reduced
        chain's stationary law and Q_i(n) from the isolated subnetwork.
        Throughput needs no conditioning: flow is fixed by the routing and an
        exact reduction leaves the chain throughput unchanged.
        """
        import time

        import numpy as _np

        from ...api.fes import fes_compute_metrics
        from ...api.io.model_adapter import ModelAdapter
        from ...api.pfqn.ljd import ljd_linearize
        from ..base import print_solver_banner

        t0 = time.time()
        subset_idx = [int(i) for i in (self.options.config or {})['fes_stations']]
        M = int(sn.nstations)
        K = int(sn.nclasses)
        if len(subset_idx) < 2:
            raise ValueError(
                "options.config['fes_stations'] must name at least two stations: "
                "collapsing one station into a flow-equivalent server saves nothing.")
        if len(set(subset_idx)) != len(subset_idx) or min(subset_idx) < 0 or max(subset_idx) >= M:
            raise ValueError(
                "options.config['fes_stations'] must be distinct 0-based station "
                "indices in 0..%d." % (M - 1))
        if len(subset_idx) >= M:
            raise ValueError(
                "options.config['fes_stations'] names every station: there is no "
                "complement left to solve.")

        stations = self.model.getStations()
        res = ModelAdapter.aggregate_fes(self.model, [stations[i] for i in subset_idx])
        fes_model = res['fes_model']
        info = res['deagg_info']

        sub = dict(self.options.config or {})
        sub['fes_stations'] = None
        inner = SolverCTMC(fes_model, config=sub, method=self.options.method,
                           verbose=self.options.verbose)
        Qr, Ur, _Rr, Tr = inner.getAvg()[:4]
        Xr = _np.atleast_1d(_np.asarray(inner.getAvgSysTput(), dtype=float)).ravel()

        # P(N_fes = n): the aggregate state space carries K columns per stateful
        # node, so the FES's block is the one at its stateful index.
        snRed = fes_model.get_struct()
        pi = _np.asarray(inner._result.pi, dtype=float).ravel()
        SSq = _np.asarray(inner.getStateSpaceAggr())
        fes_ist = int(snRed.nodeToStation[int(info['fes_node_idx'])])
        fes_isf = int(snRed.nodeToStateful[int(info['fes_node_idx'])])
        cutoffs = _np.asarray(info['cutoffs'], dtype=int).ravel()
        cols = list(range(fes_isf * K, fes_isf * K + K))
        Pn = _np.zeros(int(_np.prod(cutoffs + 1)))
        for srow in range(SSq.shape[0]):
            nvec = SSq[srow, cols]
            Pn[ljd_linearize(nvec, cutoffs) - 1] += pi[srow]

        Qtab, Utab = fes_compute_metrics(info['isolated_model'], cutoffs, K)

        QN = _np.zeros((M, K))
        UN = _np.zeros((M, K))
        TN = _np.zeros((M, K))
        comp = [int(i) for i in info['complement_indices']]
        Qr = _np.atleast_2d(_np.asarray(Qr))
        Ur = _np.atleast_2d(_np.asarray(Ur))
        Tr = _np.atleast_2d(_np.asarray(Tr))
        for a, i in enumerate(comp):
            QN[i, :] = Qr[a, :]
            UN[i, :] = Ur[a, :]
            TN[i, :] = Tr[a, :]

        sub_idx = [int(i) for i in info['subset_indices']]
        Qsub = _np.zeros((len(sub_idx), K))
        Usub = _np.zeros((len(sub_idx), K))
        for idx0 in range(len(Pn)):
            if Pn[idx0] <= 0:
                continue
            Qsub += Pn[idx0] * Qtab[idx0]
            Usub += Pn[idx0] * Utab[idx0]
        for a, i in enumerate(sub_idx):
            QN[i, :] = Qsub[a, :]
            UN[i, :] = Usub[a, :]
            # Flow through a station is fixed by the routing, so it is the FES's
            # throughput scaled by the ratio of ORIGINAL visit ratios.
            TN[i, :] = Tr[fes_ist, :] * self._fes_visit_ratio(sn, snRed, i, fes_ist, K)

        with _np.errstate(divide='ignore', invalid='ignore'):
            RN = _np.where(TN > 0, QN / _np.where(TN > 0, TN, 1.0), 0.0)

        runtime = time.time() - t0
        method = str(self.options.method) + '/fes'
        self._result = _QRFResult(QN, UN, RN, TN, RN.sum(axis=0), Xr, runtime, method)
        self._extract_names()
        if self.options.verbose:
            py_version = "%d.%d.%d" % (sys.version_info.major, sys.version_info.minor,
                                       sys.version_info.micro)
            print_solver_banner(
                "CTMC analysis [method: %s; type: %s; lang: python; env: %s] completed in %.6fs."
                % (method, method_type('CTMC', self.options.method), py_version, runtime))

    @staticmethod
    def _fes_visit_ratio(sn, snRed, ist, fes_ist, K):
        """Visits at ORIGINAL station `ist` per visit to the FES, per class."""
        import numpy as _np

        out = _np.zeros(K)
        for c in range(int(sn.nchains)):
            V = _np.asarray(sn.visits[c]) if sn.visits[c] is not None else None
            Vr = _np.asarray(snRed.visits[c]) if (snRed.visits and c < len(snRed.visits)
                                                  and snRed.visits[c] is not None) else None
            if V is None or Vr is None:
                continue
            isf = int(sn.stationToStateful[ist])
            isf_fes = int(snRed.stationToStateful[fes_ist])
            for k in range(K):
                if isf < V.shape[0] and isf_fes < Vr.shape[0] and Vr[isf_fes, k] > 0:
                    out[k] += V[isf, k] / Vr[isf_fes, k]
        return out

    def runAnalyzer(self) -> 'SolverCTMC':
        """Run the CTMC analysis."""
        # Chain mode: the generator is user-supplied, so there is no state space
        # to generate and no performance metric to derive.
        if self.isChainSolver():
            self._chain_run_analyzer()
            return self

        # unbounded open SPN rejection runs BEFORE lang=java delegation; else the JAR silently builds a cutoff-truncated wrong answer instead of rejecting.
        self._reject_unbounded_open_spn()
        # lang=java delegation populates the native result container from jline.jar; imported lazily so a JVM-free install never touches this path.
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
                line_warning("SolverCTMC", "lang='cpp' requested but the C++ solver is "
                             "unavailable (%s); falling back to lang='python'." % e)

        line_debug("CTMC: using lang=python", options=self.options)

        # Chain aggregation, opt-in through options.config['chain_aggregation'].
        # The state space of a multiclass model grows with the per-class
        # populations, so collapsing every chain onto a single class is the
        # standard way to make an otherwise intractable model solvable.
        # ModelAdapter.aggregate_chains builds the collapsed model and
        # sn_deaggregate_chain_results maps its metrics back, both of which
        # existed with no solver consumer until this branch. EXACT on a
        # product-form model, an approximation otherwise: one aggregate service
        # law, fitted to the alpha-weighted first two moments, replaces the
        # per-class ones.
        # Flow-equivalent server aggregation, opt-in through
        # options.config['fes_stations']. ModelAdapter.aggregate_fes collapses
        # the named station subset into one load-dependent station and the
        # collapsed stations' own metrics are recovered by conditioning on its
        # population; see _run_fes_aggregation. Exact when the subnetwork is
        # product-form.
        if (self.options.config or {}).get('fes_stations'):
            _snf = self._sn if getattr(self, '_sn', None) is not None else self._get_network_struct()
            self._run_fes_aggregation(_snf)
            return self

        # A user-supplied transform method name runs whichever strategy it names. The
        # inner solve carries transform='none', so a transformed submodel cannot
        # re-enter the driver.
        _tok = (self.options.config or {}).get('transform')
        if _tok and str(_tok).lower() != 'none':
            _sn = self._sn if getattr(self, '_sn', None) is not None else self._get_network_struct()
            self._run_transform(_sn)
            return self

        if (self.options.config or {}).get('chain_aggregation'):
            _sn = self._sn if getattr(self, '_sn', None) is not None else self._get_network_struct()
            if int(_sn.nchains) < int(_sn.nclasses):
                self._run_chain_aggregation(_sn)
                return self

        # reject features outside the CTMC feature set rather than silently solve a mis-specified model; mirrors MATLAB runAnalyzerChecks.
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_used_lang_features'):
            self.runAnalyzerChecks(self.options)

        # Tier B redirects: detect unsupported model families and point at the right alternative solver instead of a generic failure.
        from ...lang.base import NodeType as _NodeType
        from ...lang.base import SchedStrategy as _SchedStrategy
        sn = self._sn
        # native fork-join solves the tag-augmented copy exactly and folds auxiliary sibling classes back at the end; mirrors MATLAB SolverCTMC/runAnalyzer.m.
        self._fj_foldback = None
        if sn is not None and getattr(sn, 'nodetype', None) is not None:
            fork_v = int(_NodeType.FORK.value) if hasattr(_NodeType.FORK, 'value') else int(_NodeType.FORK)
            join_v = int(_NodeType.JOIN.value) if hasattr(_NodeType.JOIN, 'value') else int(_NodeType.JOIN)
            _has_fj = any((int(nt.value) if hasattr(nt, 'value') else int(nt)) in (fork_v, join_v)
                          for nt in sn.nodetype)
            if _has_fj:
                self._fjtag_require_network('CTMC')
                if getattr(self.options, 'timespan', None) is not None:
                    ts = np.atleast_1d(self.options.timespan)
                    if ts.size and np.isfinite(ts[0]):
                        raise RuntimeError(
                            "Transient analysis of fork-join models is not supported by SolverCTMC.")
                sn = self._fjtag_expand(sn)
                self._sn = sn

        # deadline/elapsed-time scheduling (EDD/EDF/SETF/FSP) needs per-job clocks memoryless CTMC can't represent; rejected, not solved on phase-only space.
        if sn is not None and getattr(sn, 'sched', None) is not None:
            _unsupported_sched = {
                _SchedStrategy.EDD: 'EDD', _SchedStrategy.EDF: 'EDF',
                _SchedStrategy.SETF: 'SETF', _SchedStrategy.FSP: 'FSP',
            }
            _nt = np.asarray(sn.nodetype).ravel() if getattr(sn, 'nodetype', None) is not None else None
            # sn.sched may be a dict {station_idx: SchedStrategy} or an array
            _sched_items = sn.sched.items() if hasattr(sn.sched, 'items') else enumerate(sn.sched)
            for st_idx, sc in _sched_items:
                try:
                    sc_enum = _SchedStrategy(int(sc.value) if hasattr(sc, 'value') else int(sc))
                except (ValueError, TypeError):
                    continue
                if sc_enum in _unsupported_sched:
                    raise RuntimeError(
                        "This model uses the %s scheduling strategy, which is not "
                        "supported by SolverCTMC (deadline/elapsed-time policies require "
                        "a per-job clock outside a memoryless CTMC). Use SolverLDES or "
                        "SolverSSA for simulation." % _unsupported_sched[sc_enum]
                    )
                # EXT scheduling only makes sense at a Source; at a regular station a memoryless CTMC cannot generate arrivals and the state space diverges.
                if sc_enum == _SchedStrategy.EXT and _nt is not None:
                    try:
                        node_idx = int(sn.stationToNode[st_idx])
                    except (TypeError, IndexError, KeyError):
                        node_idx = -1
                    if 0 <= node_idx < _nt.size and int(_nt[node_idx]) != int(NodeType.SOURCE):
                        raise RuntimeError(
                            "This model applies EXT (external-arrival) scheduling to a "
                            "non-Source station, which is not supported by SolverCTMC. "
                            "EXT is only valid at a Source; use a standard scheduling "
                            "policy (FCFS/PS/...) at queues, or SolverLDES/SolverSSA."
                        )

        # FCR enforced in CTMC handler by filtering the state space to aggregate per-region caps (blocking-before-entry); per-station setCapacity honored.

        # reneging models exponential patience via a RENEGE event at rate waiting*mu; PH/MAP patience needs a per-job phase dimension and is left to LDES/JMT.

        # retrial models the exponential-delay, unlimited-attempt, single-class case; other configs are rejected, not silently mis-solved with no retry.

        # signal classes never occupy a station; 0 per-station cap (except EXT/Source) omits unreachable signal-holding states. REPLY exempt: ordinary job.
        if (sn is not None and getattr(sn, 'issignal', None) is not None
                and getattr(sn, 'classcap', None) is not None):
            from ...api.state.reply_block import is_reply_class as _is_reply_class
            for _ist in range(int(sn.nstations)):
                if sn.sched[_ist] != _SchedStrategy.EXT:
                    for _r in range(int(sn.nclasses)):
                        if sn.issignal[_r] and not _is_reply_class(sn, _r):
                            sn.classcap[_ist, _r] = 0

        # heterogeneous ORDER-policy servers map to load-dependent mu(n) = sum of the first min(n,c) server rates; multi-class hetero is rejected.
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
                        "SolverCTMC supports heterogeneous servers only for single-class "
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

        # open SPN solved as bounded CTMC if marking bounded (finite Place cap or cutoff); only unbounded net rejected (see _reject_unbounded_open_spn).

        # QRF (Quadratic/Linear Reduction Framework) LP-based bounds moved out of
        # SolverCTMC into SolverBA, which is where listValidMethods stopped
        # naming them. This entry point outlived the move: it kept dispatching
        # to solver_ctmc_qrf_analyzer, the very analyzer SolverBA itself calls
        # (solver_ba_analyzer:435), so it was a second front door onto one
        # computation -- measured identical, QLen [1.7778, 0.22222] on a
        # two-queue closed model either way. The text lives in
        # unsupportedMethodReason, which runAnalyzerChecks asks BEFORE it
        # reports an unlisted method; this call is what still refuses on the
        # enableChecks = False path, which skips that gate entirely.
        moved_qrf = self.unsupportedMethodReason(self.options.method)
        if moved_qrf:
            raise RuntimeError(moved_qrf)

        # The 'mdd' method never builds the explicit generator, so it returns
        # before the state-space path below and leaves the state space empty by
        # design.
        if self.options.method.lower() == 'mdd':
            line_debug("Using MDD level aggregation for steady-state CTMC analysis",
                       options=self.options)
            import time
            from ...api.solvers.ctmc.solver_ctmc_mdd_analyzer import solver_ctmc_mdd_analyzer
            t0 = time.time()
            QN, UN, RN, TN, CN, XN, mddinfo = solver_ctmc_mdd_analyzer(self._sn, self.options, self.model)
            runtime = time.time() - t0
            self._result = _MDDResult(QN, UN, RN, TN, CN, XN, runtime, self.options.method,
                                      mddinfo)
            self._extract_names()
            if self.options.verbose:
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"CTMC analysis [method: {self.options.method}; type: {method_type('CTMC', self.options.method)}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")
            return self

        # Perfect sampling replaces enumeration: intercepted before the state space
        # is built, so the memory guard below never applies to it.
        if self.options.method.lower().startswith('cftp'):
            line_debug("Using perfect sampling for steady-state CTMC analysis", options=self.options)
            import time
            from ...api.solvers.ctmc.solver_ctmc_cftp_analyzer import solver_ctmc_cftp_analyzer
            t0 = time.time()
            QN, UN, RN, TN, CN, XN, Xs, Ts, pAggr, SSq = solver_ctmc_cftp_analyzer(self._sn, self.options)
            runtime = time.time() - t0
            self._result = _CFTPResult(QN, UN, RN, TN, CN, XN, runtime, self.options.method,
                                       Xs, Ts, pAggr, SSq)
            self._extract_names()
            if self.options.verbose:
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"CTMC analysis [method: {self.options.method}; type: {method_type('CTMC', self.options.method)}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")
            return self

        from ...api.solvers.ctmc.handler import (
            solver_ctmc, SolverCTMCOptions as HandlerOptions
        )
        from ...api.sn import sn_nonmarkov_toph

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            cutoff=self.options.cutoff,
            verbose=self.options.verbose,
            force=self.options.force,
            gen_method=getattr(self.options, 'gen_method', 'default'),
            # the wall-clock budget and the state cap bound the solve only if they reach the handler.
            timeout=getattr(self.options, 'timeout', float('inf')),
            ctmc_max_states=getattr(self.options, 'ctmc_max_states', 3_000_000),
            memory_safety_fraction=getattr(self.options, 'memory_safety_fraction', 0.6),
        )

        # Log open/mixed cutoff if applicable
        sn = self._sn
        if sn is not None and sn.njobs is not None and np.any(np.isinf(sn.njobs)):
            line_debug("Open/mixed model: cutoff=%d for %d stations, %d classes",
                       self.options.cutoff, sn.nstations, sn.nclasses, options=self.options)

        line_debug("Using standard CTMC method for steady-state analysis", options=self.options)

        # Convert non-Markovian distributions to phase-type
        # Convert options object to dictionary for sn_nonmarkov_toph
        options_dict = vars(self.options) if hasattr(self.options, '__dict__') else {'config': {}}
        sn = sn_nonmarkov_toph(sn, options_dict)
        self._sn = sn

        line_debug("CTMC: converted non-Markovian distributions to PH (nstations=%d, nclasses=%d)",
                   sn.nstations if sn is not None else 0, sn.nclasses if sn is not None else 0, options=self.options)

        # Run the solver. The result is published to self._result only once the
        # post-processing below has finalized Q/U/R/T: the setter applies the
        # getAvg near-zero mask, which must see the folded-back matrices.
        r = solver_ctmc(sn, handler_options)

        # native fork-join: fold sibling-class metrics back, restore pre-augmentation struct; raw pi/space/infgen stay fjsn-based for state-probability APIs.
        if getattr(self, '_fj_foldback', None) is not None:
            # SolverCTMC's result container carries no arrival-rate field, so
            # the AN the lift returns is used for RN inside it and dropped here.
            self._fjtag_lift(r)
        else:
            from ...api.sn.getters import sn_pn_avg_rates
            r.T, _, r.R = sn_pn_avg_rates(sn, r.Q, r.T, None, r.R)

        self._result = r

        # Extract station and class names
        self._extract_names()

        # Compute cache hit/miss probabilities for cache nodes
        self._compute_cache_hit_miss_probs()
        self._compute_cache_item_prob()
        self._compute_cache_delayed_hit_qlen()

        # After computing actual cache probs, refresh routing and visits
        # (matches MATLAB runAnalyzer.m: setResultHitProb -> refreshChains)
        self._refresh_cache_routing_and_visits()

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            from ..base import method_label
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            runtime = self._result.runtime if hasattr(self._result, 'runtime') else 0.0
            method = self._result.method if hasattr(self._result, 'method') else 'default'
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"CTMC analysis [method: {method_label(self.options.method, method)}; type: {method_type('CTMC', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")

        return self

    def _extract_names(self):
        """Extract station and class names from network struct."""
        if self._sn is not None:
            # Get station names by mapping station indices to node indices
            if hasattr(self._sn, 'stationToNode') and self._sn.stationToNode is not None and \
               hasattr(self._sn, 'nodenames') and self._sn.nodenames:
                import numpy as np
                station_to_node = np.asarray(self._sn.stationToNode).flatten()
                self.station_names = []
                for ist in range(self._sn.nstations):
                    if ist < len(station_to_node):
                        node_idx = int(station_to_node[ist])
                        if node_idx >= 0 and node_idx < len(self._sn.nodenames):
                            self.station_names.append(self._sn.nodenames[node_idx])
                        else:
                            self.station_names.append(f'Station{ist}')
                    else:
                        self.station_names.append(f'Station{ist}')
            elif hasattr(self._sn, 'nodenames') and self._sn.nodenames:
                # Fallback: use first nstations node names (for simple networks)
                self.station_names = list(self._sn.nodenames[:self._sn.nstations])
            else:
                self.station_names = [f'Station{i}' for i in range(self._sn.nstations)]

            self.class_names = list(self._sn.classnames) if hasattr(self._sn, 'classnames') and self._sn.classnames else \
                              [f'Class{i}' for i in range(self._sn.nclasses)]
        else:
            self.station_names = []
            self.class_names = []

    def _stationary_state_cols(self):
        """Stationary vector, state space and the column range of each stateful node.

        Returns (pi, space, state_cols) or None when the CTMC result does not
        carry an explicit state space.
        """
        sn = self._sn
        if sn is None or self._result is None:
            return None
        pi = getattr(self._result, 'pi', None)
        space = getattr(self._result, 'space', None)
        if pi is None or space is None:
            return None
        pi = np.ravel(np.asarray(pi, dtype=float))
        space = np.atleast_2d(np.asarray(space, dtype=float))
        if space.shape[0] != pi.shape[0]:
            return None
        widths = getattr(self._result, 'node_space_width', None)
        if not widths:
            return None
        col_off = 0
        state_cols = {}
        for isf in range(sn.nstateful):
            w = int(widths.get(isf, 0))
            state_cols[isf] = (col_off, col_off + w)
            col_off += w
        if col_off != space.shape[1]:
            return None
        return pi, space, state_cols

    def _compute_cache_item_prob(self):
        """Time-stationary per-item occupancy of each cache list.

        The cache-contents block of the local-variable vector holds the item
        index resident in each cache position, so P(item i is held by list l) is
        a state reward of the stationary distribution. This is the TIME-WEIGHTED
        occupancy, the CTMC counterpart of the EMBEDDED (per-request) occupancy
        the NC/MVA cache algorithms return; the two coincide only when requests
        see time averages (PASTA).

        Port of the per-item block in MATLAB solver_ctmc_analyzer.m.
        """
        sn = self._sn
        got = self._stationary_state_cols()
        if got is None:
            return
        pi, space, state_cols = got
        nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        for ind in range(sn.nnodes):
            if sn.nodetype[ind] != NodeType.CACHE:
                continue
            np_ = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
            if np_ is None:
                continue
            itemcap = np.atleast_1d(np.asarray(getattr(np_, 'itemcap', []), dtype=int)).ravel()
            nitems = int(getattr(np_, 'nitems', 0))
            if itemcap.size == 0 or nitems == 0:
                continue
            isf = int(sn.nodeToStateful[ind])
            c0, c1 = state_cols[isf]
            if int(getattr(np_, 'retrieval_system_capacity', 0)) > 0:
                from ...api.state.ctmc_ssg import cache_retrieval_class_map
                _, rc_items_all, _ = cache_retrieval_class_map(sn, ind)
                lvw = int(np.sum(itemcap)) + nitems + len(rc_items_all)
            else:
                lvw = int(np.sum(itemcap))
            lvs = (c1 - c0) - lvw  # per-class server presence width
            if lvs < 0:
                continue
            itemprob = np.zeros((nitems, itemcap.size + 1))
            off = 0
            for l in range(itemcap.size):
                lcols = [c0 + lvs + off + q for q in range(int(itemcap[l]))]
                off += int(itemcap[l])
                for i in range(nitems):
                    inlist = np.any(space[:, lcols] == (i + 1), axis=1)
                    itemprob[i, l + 1] = float(np.sum(pi[inlist]))
            itemprob[:, 0] = 1.0 - np.sum(itemprob[:, 1:], axis=1)
            np_.actualitemprob = itemprob
            if ind < len(nodes) and hasattr(nodes[ind], 'set_result_item_prob'):
                nodes[ind].set_result_item_prob(itemprob)

    def _compute_cache_delayed_hit_qlen(self):
        """Exact delayed-hit queue length of a retrieval-system cache.

        Block A of the cache local-variable vector marks the items being fetched
        and block B counts, per retrieval class, the secondary requests merged
        onto those fetches, so

            phi_i   = P(a fetch of item i is in flight)
            d1_i    = E[secondary requests waiting on the fetch of item i]
            dfull_i = d1_i + phi_i

        are state rewards of the stationary distribution, hence exact.

        Port of the delayed-hit block in MATLAB solver_ctmc_analyzer.m.
        """
        sn = self._sn
        got = self._stationary_state_cols()
        if got is None:
            return
        pi, space, state_cols = got
        from ...api.state.ctmc_ssg import cache_retrieval_class_map

        nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        for ind in range(sn.nnodes):
            if sn.nodetype[ind] != NodeType.CACHE:
                continue
            np_ = sn.nodeparam[ind] if sn.nodeparam is not None and ind in sn.nodeparam else None
            if np_ is None or getattr(np_, 'retrieval_system_capacity', 0) <= 0:
                continue
            _, rc_items, rc_orig = cache_retrieval_class_map(sn, ind)
            nitems = int(getattr(np_, 'nitems', 0))
            tcc = int(getattr(np_, 'total_cache_capacity', 0))
            isf = int(sn.nodeToStateful[ind])
            c0, c1 = state_cols[isf]
            lvs = (c1 - c0) - (tcc + nitems + len(rc_items))
            if lvs < 0:
                continue
            a0 = c0 + lvs + tcc
            b0 = a0 + nitems
            phi = np.zeros(nitems)
            d1 = np.zeros(nitems)
            for i in range(nitems):
                phi[i] = float(np.sum(pi[space[:, a0 + i] != 0]))
                bsel = [b0 + j for j, it in enumerate(rc_items) if it == i + 1]
                if bsel:
                    d1[i] = float(np.sum(pi * np.sum(space[:, bsel], axis=1)))
            np_.delayedhitprobitem = phi
            np_.delayedhitqlen = d1
            np_.delayedhitqlenfull = d1 + phi
            if ind < len(nodes) and hasattr(nodes[ind], 'set_result_delayed_hit_qlen'):
                nodes[ind].set_result_delayed_hit_qlen(d1, d1 + phi)

            # Exact delayed-hit rate per originating class. A fetch of item i completes
            # on exactly the transitions that clear block A bit i, and each such
            # transition releases the block-B counts of item i as delayed hits. The rate
            # is therefore a TRANSITION reward over the generator, not a state reward:
            # the alternative arrival-rate identity lambda_i*phi_i is only PASTA-exact.
            Q = getattr(self._result, 'infgen', None)
            if Q is None:
                continue
            Q = np.asarray(Q.todense()) if hasattr(Q, 'todense') else np.asarray(Q)
            offdiag = Q - np.diag(np.diag(Q))
            delayed_rate = np.zeros(sn.nclasses)
            for j, item in enumerate(rc_items):
                i = item - 1
                rows = np.where((space[:, a0 + i] != 0) & (space[:, b0 + j] > 0))[0]
                for rr in rows:
                    nz = np.where(offdiag[rr] != 0)[0]
                    completes = nz[space[nz, a0 + i] == 0]
                    if completes.size == 0:
                        continue
                    delayed_rate[rc_orig[j]] += (pi[rr] * space[rr, b0 + j]
                                                 * float(np.sum(offdiag[rr, completes])))
            np_.delayedhitrate = delayed_rate
            self._apply_delayed_hit_split(sn, ind, np_, delayed_rate, nodes)

    def _apply_delayed_hit_split(self, sn, ind, np_, delayed_rate, nodes):
        """Split the cache hit-class rate into true hits and delayed hits.

        Delayed hits depart in the hit class, so the hit-class rate is
        (true hits + delayed hits); the exact delayed rate splits it so that
        hit + delayed + miss = 1, matching the LDES/NC report.
        """
        hitclass = np.atleast_1d(np.asarray(getattr(np_, 'hitclass', []), dtype=int))
        missclass = np.atleast_1d(np.asarray(getattr(np_, 'missclass', []), dtype=int))
        hp = getattr(np_, 'actualhitprob', None)
        mp = getattr(np_, 'actualmissprob', None)
        if hp is None or mp is None:
            return
        hp = np.array(hp, dtype=float, copy=True)
        mp = np.array(mp, dtype=float, copy=True)
        dp = np.zeros_like(hp)
        depRates = getattr(self._result, 'depRates', None)
        pi = np.ravel(np.asarray(getattr(self._result, 'pi', []), dtype=float))
        isf = int(sn.nodeToStateful[ind])
        for k in range(min(len(hitclass), len(hp))):
            h, m = int(hitclass[k]), int(missclass[k])
            if h < 0 or m < 0 or depRates is None or isf >= np.asarray(depRates).shape[1]:
                continue
            tn_hit = float(np.dot(pi, np.asarray(depRates)[:, isf, h]))
            tn_miss = float(np.dot(pi, np.asarray(depRates)[:, isf, m]))
            denom = tn_hit + tn_miss
            if denom <= 0:
                continue
            d = min(delayed_rate[k] if k < len(delayed_rate) else 0.0, tn_hit)
            hp[k] = (tn_hit - d) / denom
            dp[k] = d / denom
            mp[k] = tn_miss / denom
        np_.actualhitprob = hp
        np_.actualmissprob = mp
        np_.actualdelayedhitprob = dp
        if ind < len(nodes):
            cn = nodes[ind]
            if hasattr(cn, 'set_result_hit_prob'):
                cn.set_result_hit_prob(hp)
            if hasattr(cn, 'set_result_miss_prob'):
                cn.set_result_miss_prob(mp)
            if hasattr(cn, 'set_result_delayed_hit_prob'):
                cn.set_result_delayed_hit_prob(dp)

    def _compute_cache_hit_miss_probs(self):
        """
        Compute actual hit/miss probabilities for cache nodes from CTMC results.

        This matches MATLAB's solver_ctmc_analyzer behavior where it computes
        actualhitprob and actualmissprob for cache nodes using the departure rates
        from the Markov chain stationary distribution.

        The formula is:
            TNcache[ist, k] = pi @ depRates[:, ist, k]
            actualhitprob[orig_class] = TNcache[ist, hitclass] / (TNcache[ist, hitclass] + TNcache[ist, missclass])

        References:
            MATLAB: solver_ctmc_analyzer.m lines 215-240
        """
        if self._sn is None or self._result is None:
            return

        sn = self._sn
        K = sn.nclasses
        I = sn.nnodes
        M = sn.nstations

        # Find cache nodes
        cache_nodes = []
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            for ind in range(len(sn.nodetype)):
                if sn.nodetype[ind] == NodeType.CACHE:
                    cache_nodes.append(ind)

        if not cache_nodes:
            return

        # Get stationary distribution and departure rates
        pi = self._result.pi
        depRates = getattr(self._result, 'depRates', None)

        if depRates is None or pi is None:
            # Fallback to routing matrix approach if depRates not available
            self._compute_cache_hit_miss_probs_from_routing()
            return

        # Compute cache throughputs for each class using pi @ depRates
        # TNcache[ist, k] = sum over states s of: pi[s] * depRates[s, ist, k]
        TNcache = np.zeros((M, K))
        for ist in range(M):
            for k in range(K):
                TNcache[ist, k] = np.dot(pi, depRates[:, ist, k])

        # For each cache node, compute actual hit/miss probabilities
        for cache_ind in cache_nodes:
            if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                continue

            node_param = sn.nodeparam[cache_ind]
            if not hasattr(node_param, 'hitclass') or not hasattr(node_param, 'missclass'):
                continue

            hitclass = np.atleast_1d(node_param.hitclass).flatten()
            missclass = np.atleast_1d(node_param.missclass).flatten()

            # Initialize actual probabilities
            actual_hit_prob = np.zeros(K)
            actual_miss_prob = np.zeros(K)

            # Get cache station index
            cache_ist = int(sn.nodeToStation[cache_ind]) if cache_ind < len(sn.nodeToStation) else -1

            if cache_ist >= 0:
                # Cache IS a station - use departure rates to compute hit/miss probs
                # For each input class that has hit/miss classes defined
                for orig_class in range(len(hitclass)):
                    h = int(hitclass[orig_class]) if orig_class < len(hitclass) else -1
                    m = int(missclass[orig_class]) if orig_class < len(missclass) else -1

                    if h < 0 or m < 0 or h >= K or m >= K:
                        continue

                    # cache hit/miss throughput read from TNcache at the hit/miss class columns.
                    TN_hit = TNcache[cache_ist, h]
                    TN_miss = TNcache[cache_ist, m]
                    TN_total = TN_hit + TN_miss

                    if TN_total > 0:
                        actual_hit_prob[orig_class] = TN_hit / TN_total
                        actual_miss_prob[orig_class] = TN_miss / TN_total
            else:
                # sync builder: hit/miss via pi @ depRates[:,cache_sf,hit|miss]; flat builder station-indexed depRates lacks cache column, falls back to stationary.
                cache_sf = int(sn.nodeToStateful[cache_ind]) if cache_ind < len(sn.nodeToStateful) else -1
                if depRates is not None and 0 <= cache_sf < depRates.shape[1]:
                    for orig_class in range(len(hitclass)):
                        h = int(hitclass[orig_class]) if orig_class < len(hitclass) else -1
                        m = int(missclass[orig_class]) if orig_class < len(missclass) else -1
                        if h < 0 or m < 0 or h >= K or m >= K:
                            continue
                        TN_hit = float(np.dot(pi, depRates[:, cache_sf, h]))
                        TN_miss = float(np.dot(pi, depRates[:, cache_sf, m]))
                        if TN_hit + TN_miss > 0:
                            actual_hit_prob[orig_class] = TN_hit / (TN_hit + TN_miss)
                            actual_miss_prob[orig_class] = TN_miss / (TN_hit + TN_miss)
                    node_param.actualhitprob = actual_hit_prob
                    node_param.actualmissprob = actual_miss_prob
                    latency_result = self._compute_retrieval_latency(
                        sn, node_param, cache_ind, TNcache)
                    actual_latency = (latency_result[0]
                                      if latency_result is not None else None)
                    if hasattr(self.model, 'get_nodes'):
                        nodes = self.model.get_nodes()
                        if cache_ind < len(nodes):
                            cn = nodes[cache_ind]
                            if hasattr(cn, 'set_result_hit_prob'):
                                cn.set_result_hit_prob(actual_hit_prob)
                            if hasattr(cn, 'set_result_miss_prob'):
                                cn.set_result_miss_prob(actual_miss_prob)
                            if (actual_latency is not None and
                                    hasattr(cn, 'set_result_residt')):
                                cn.set_result_residt(actual_latency)
                    if actual_latency is not None:
                        node_param.actualresidt = actual_latency
                    continue

                # Cache is not a station: hit probability computed from the stationary distribution as sum(pi[s]*P(HIT|s)).
                from ...api.solvers.ctmc.handler import _enumerate_cache_states, _get_cache_stations_info

                # Get cache states info from the result's rrobin_info
                space = getattr(self._result, 'space', None)
                rrobin_info = getattr(self._result, 'rrobin_info', None)

                if space is None or pi is None:
                    continue

                # Get cache parameters
                pread = node_param.pread if hasattr(node_param, 'pread') else None
                nitems = node_param.nitems if hasattr(node_param, 'nitems') else 0
                itemcap = node_param.itemcap if hasattr(node_param, 'itemcap') else None
                capacity = int(sum(itemcap)) if itemcap is not None else 0

                if pread is None or nitems <= 0 or capacity <= 0:
                    continue

                # Get cache state variable offset from rrobin_info
                cache_isf = int(sn.nodeToStateful[cache_ind]) if cache_ind < len(sn.nodeToStateful) else -1
                cache_state_offset = None
                if rrobin_info is not None and 'cache_state_offsets' in rrobin_info:
                    cache_state_offset = rrobin_info['cache_state_offsets'].get(cache_isf)

                if cache_state_offset is None:
                    continue

                # cache states enumerated as [cache | retrieval slots] with a retrieval system, matching state-vector index built in _get_cache_stations_info.
                retrieval_capacity = int(
                    getattr(node_param, 'retrieval_system_capacity', 0))
                tcc = int(getattr(node_param, 'total_cache_capacity', capacity))
                cache_states = _enumerate_cache_states(
                    nitems, tcc, retrieval_capacity)
                n_cache_states = len(cache_states)

                # For each input class that has hit/miss classes defined
                for orig_class in range(len(hitclass)):
                    h = int(hitclass[orig_class]) if orig_class < len(hitclass) else -1
                    m = int(missclass[orig_class]) if orig_class < len(missclass) else -1

                    if h < 0 or m < 0 or h >= K or m >= K:
                        continue

                    # Get pread for this class
                    pread_k = None
                    if isinstance(pread, (list, tuple)) and orig_class < len(pread):
                        pread_k = pread[orig_class]
                    elif isinstance(pread, np.ndarray):
                        if pread.ndim == 1:
                            pread_k = pread
                        elif orig_class < pread.shape[0]:
                            pread_k = pread[orig_class]

                    if pread_k is None:
                        continue

                    pread_k = np.atleast_1d(pread_k).flatten()

                    # expected hit ratio counts both cached and retrieving items as hits (delayed hits count), matching [cache | retrieval slots] layout.
                    expected_hit_prob = 0.0
                    for s_idx, state in enumerate(space):
                        cache_state_idx = int(state[cache_state_offset])
                        if cache_state_idx < n_cache_states:
                            row = cache_states[cache_state_idx].tolist()
                            cache_content = {int(x) for x in row if x != 0}
                            # P(HIT | state s) = sum(pread[item-1] for item present)
                            state_hit_prob = 0.0
                            for item in cache_content:
                                if 1 <= item <= len(pread_k):
                                    state_hit_prob += pread_k[item - 1]
                            expected_hit_prob += pi[s_idx] * state_hit_prob

                    actual_hit_prob[orig_class] = expected_hit_prob
                    actual_miss_prob[orig_class] = 1.0 - expected_hit_prob

            # Store in nodeparam
            node_param.actualhitprob = actual_hit_prob
            node_param.actualmissprob = actual_miss_prob

            # Expected latency per input class via Little's law on the
            # retrieval queues.
            latency_result = self._compute_retrieval_latency(
                sn, node_param, cache_ind, TNcache)
            actual_latency = (latency_result[0]
                              if latency_result is not None else None)

            # Also set on the actual Cache node (matching MATLAB's runAnalyzer.m)
            if hasattr(self.model, 'get_nodes'):
                nodes = self.model.get_nodes()
                if cache_ind < len(nodes):
                    cache_node = nodes[cache_ind]
                    if hasattr(cache_node, 'set_result_hit_prob'):
                        cache_node.set_result_hit_prob(actual_hit_prob)
                    if hasattr(cache_node, 'set_result_miss_prob'):
                        cache_node.set_result_miss_prob(actual_miss_prob)
                    if (actual_latency is not None and
                            hasattr(cache_node, 'set_result_residt')):
                        cache_node.set_result_residt(actual_latency)
            if actual_latency is not None:
                node_param.actualresidt = actual_latency

    def _compute_retrieval_latency(self, sn, node_param, cache_ind, TNcache):
        """
        Retrieval-system expected latency (Sala et al. 2026 Eq. 8) is not
        currently implemented. For any cache configured with a retrieval
        system this returns all-NaN latencies (and emits a one-shot warning);
        it returns None when the cache has no retrieval system configured.
        """
        K = sn.nclasses
        rpc = getattr(node_param, 'retrieval_classes', None)
        rsqi = getattr(node_param, 'retrieval_system_queue_indices', None)
        if rpc is None or rsqi is None:
            return None
        rpc_arr = np.atleast_2d(np.asarray(rpc))
        if rpc_arr.size == 0:
            return None
        # The Eq. 8 retrieval-system expected latency is not currently
        # implemented; report NaN whenever a retrieval system is configured.
        if isinstance(rsqi, dict) and any(rsqi.get(k) for k in range(K)):
            line_warning('solver_ctmc_analyzer',
                         'Retrieval-system expected latency is not currently '
                         'implemented; reporting NaN.')
        return np.full(K, np.nan), None

    def _compute_cache_hit_miss_probs_from_routing(self):
        """
        Fallback: compute hit/miss probabilities from routing matrix.

        Used when depRates are not available (e.g., older implementations).
        """
        sn = self._sn
        K = sn.nclasses
        I = sn.nnodes

        # Find cache nodes
        cache_nodes = []
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            for ind in range(len(sn.nodetype)):
                if sn.nodetype[ind] == NodeType.CACHE:
                    cache_nodes.append(ind)

        for cache_ind in cache_nodes:
            if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                continue

            node_param = sn.nodeparam[cache_ind]
            if not hasattr(node_param, 'hitclass') or not hasattr(node_param, 'missclass'):
                continue

            hitclass = np.atleast_1d(node_param.hitclass).flatten()
            missclass = np.atleast_1d(node_param.missclass).flatten()

            actual_hit_prob = np.zeros(K)
            actual_miss_prob = np.zeros(K)

            for orig_class in range(len(hitclass)):
                h = int(hitclass[orig_class]) if orig_class < len(hitclass) else -1
                m = int(missclass[orig_class]) if orig_class < len(missclass) else -1

                if h < 0 or m < 0 or h >= K or m >= K:
                    continue

                hit_prob = 0.0
                miss_prob = 0.0

                if hasattr(sn, 'rtnodes') and sn.rtnodes is not None:
                    for jnd in range(I):
                        from_idx = cache_ind * K + orig_class
                        to_h_idx = jnd * K + h
                        to_m_idx = jnd * K + m
                        if from_idx < sn.rtnodes.shape[0]:
                            if to_h_idx < sn.rtnodes.shape[1]:
                                p = sn.rtnodes[from_idx, to_h_idx]
                                if p > 0:
                                    hit_prob += p
                            if to_m_idx < sn.rtnodes.shape[1]:
                                p = sn.rtnodes[from_idx, to_m_idx]
                                if p > 0:
                                    miss_prob += p

                total_prob = hit_prob + miss_prob
                if total_prob > 0:
                    actual_hit_prob[orig_class] = hit_prob / total_prob
                    actual_miss_prob[orig_class] = miss_prob / total_prob

            node_param.actualhitprob = actual_hit_prob
            node_param.actualmissprob = actual_miss_prob

            # Also set on the actual Cache node (matching MATLAB's runAnalyzer.m)
            if hasattr(self.model, 'get_nodes'):
                nodes = self.model.get_nodes()
                if cache_ind < len(nodes):
                    cache_node = nodes[cache_ind]
                    if hasattr(cache_node, 'set_result_hit_prob'):
                        cache_node.set_result_hit_prob(actual_hit_prob)
                    if hasattr(cache_node, 'set_result_miss_prob'):
                        cache_node.set_result_miss_prob(actual_miss_prob)

    def _refresh_cache_routing_and_visits(self):
        """
        After computing actual cache hit/miss probabilities from CTMC results,
        update the routing matrix (sn.rt) and refresh visit ratios.

        This matches MATLAB's runAnalyzer.m pattern:
            setResultHitProb -> refreshChains() -> updated visits
        which rebuilds routing with actual cache probs and recomputes visits.

        The pre-analysis routing used estimated probs (from cache_xi_fp).
        After analysis, we have actual probs from the stationary distribution,
        so we re-combine hit/miss class routing with the actual probabilities.
        """
        from ...api.sn.transforms import sn_refresh_visits

        sn = self._sn
        if sn is None:
            return

        K = sn.nclasses

        # Find cache nodes
        cache_nodes = []
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            for ind in range(len(sn.nodetype)):
                if sn.nodetype[ind] == NodeType.CACHE:
                    cache_nodes.append(ind)

        if not cache_nodes:
            return

        for cache_ind in cache_nodes:
            if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                continue

            ch = sn.nodeparam[cache_ind]
            hitclass = getattr(ch, 'hitclass', None)
            missclass = getattr(ch, 'missclass', None)
            actualhitprob = getattr(ch, 'actualhitprob', None)
            actualmissprob = getattr(ch, 'actualmissprob', None)

            if hitclass is None or missclass is None:
                continue
            if actualhitprob is None or actualmissprob is None:
                continue

            hitclass = np.atleast_1d(hitclass).flatten()
            missclass = np.atleast_1d(missclass).flatten()

            # Update sn.rtnodes with actual probs
            if sn.rtnodes is not None:
                I = sn.nnodes
                for r in range(len(hitclass)):
                    if r >= K:
                        break
                    hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
                    mc = int(missclass[r]) if missclass[r] >= 0 else -1
                    if hc < 0 or mc < 0 or hc >= K or mc >= K:
                        continue

                    hit_prob = actualhitprob[r] if r < len(actualhitprob) else 0.5
                    miss_prob = actualmissprob[r] if r < len(actualmissprob) else 0.5

                    # Zero out input class routing row
                    sn.rtnodes[cache_ind * K + r, :] = 0

                    # Set routing to connected nodes using hit/miss probs
                    for jnd in range(I):
                        if sn.connmatrix is not None and cache_ind < sn.connmatrix.shape[0] and jnd < sn.connmatrix.shape[1]:
                            if sn.connmatrix[cache_ind, jnd]:
                                if hc >= 0 and hc < K:
                                    sn.rtnodes[cache_ind * K + r, jnd * K + hc] = hit_prob
                                if mc >= 0 and mc < K:
                                    sn.rtnodes[cache_ind * K + r, jnd * K + mc] = miss_prob

            # Update sn.rt with actual probs (stateful-indexed)
            if sn.rt is not None:
                cache_sf = int(sn.nodeToStateful[cache_ind]) if sn.nodeToStateful is not None and cache_ind < len(sn.nodeToStateful) else -1
                if cache_sf < 0:
                    continue

                for r in range(len(hitclass)):
                    if r >= K:
                        break
                    hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
                    mc = int(missclass[r]) if missclass[r] >= 0 else -1
                    if hc < 0 or mc < 0 or hc >= K or mc >= K:
                        continue

                    hit_prob = actualhitprob[r] if r < len(actualhitprob) else 0.5
                    miss_prob = actualmissprob[r] if r < len(actualmissprob) else 0.5

                    input_src_idx = cache_sf * K + r
                    hit_src_idx = cache_sf * K + hc
                    miss_src_idx = cache_sf * K + mc

                    if input_src_idx >= sn.rt.shape[0]:
                        continue

                    # input-class routing zeroed and recombined with actual probabilities in sn.rt and sn.rt_visits (used by sn_refresh_visits for Sink->Source folding).
                    for rt_matrix in [sn.rt] + ([sn.rt_visits] if hasattr(sn, 'rt_visits') and sn.rt_visits is not None and sn.rt_visits is not sn.rt else []):
                        if input_src_idx >= rt_matrix.shape[0]:
                            continue
                        # Save hit/miss routing before zeroing (in case input row overlaps)
                        hit_routes = rt_matrix[hit_src_idx, :].copy() if hit_src_idx < rt_matrix.shape[0] else np.zeros(rt_matrix.shape[1])
                        miss_routes = rt_matrix[miss_src_idx, :].copy() if miss_src_idx < rt_matrix.shape[0] else np.zeros(rt_matrix.shape[1])
                        rt_matrix[input_src_idx, :] = 0
                        for dst_idx in range(rt_matrix.shape[1]):
                            combined_prob = hit_prob * hit_routes[dst_idx] + miss_prob * miss_routes[dst_idx]
                            if combined_prob > 1e-10:
                                rt_matrix[input_src_idx, dst_idx] = combined_prob

        # Refresh visit ratios with updated routing
        sn_refresh_visits(sn)

    # =========================================================================
    # Table Output
    # =========================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns node-level results (one row per node per class) to match MATLAB output format.
        Non-station nodes (e.g., Fork, ClassSwitch) are included with computed metrics.
        Cache nodes include HitClass/MissClass throughputs using actual hit/miss probabilities.

        Returns:
            pandas.DataFrame with columns: Node, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        self._assert_not_chain_model('getAvgTable')
        if self._result is None:
            self._ensureAvgResults()

        sn = self._sn
        M = self._result.Q.shape[0]  # nstations
        K = self._result.Q.shape[1]  # nclasses
        I = sn.nnodes if sn is not None else M

        # Get station-level results (make copies to avoid modifying originals)
        QN = self._result.Q.copy()
        UN = self._result.U.copy()
        RN = self._result.R.copy()
        TN = self._result.T.copy()

        # metrics zeroed for classes with zero visit ratio, before node-level computations; mirrors MATLAB getAvg.m:163-180.
        hasForkJoin = False
        hasSPN = False
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            hasForkJoin = np.any(sn.nodetype == NodeType.FORK) and np.any(sn.nodetype == NodeType.JOIN)
            hasSPN = np.any(sn.nodetype == NodeType.PLACE) or np.any(sn.nodetype == NodeType.TRANSITION)

        if sn is not None and hasattr(sn, 'nchains') and sn.nchains > 0 and not hasSPN:
            if hasattr(sn, 'chains') and sn.chains is not None and hasattr(sn, 'visits') and sn.visits:
                chains_arr = np.asarray(sn.chains)
                for k in range(K):
                    # Find chains containing this class
                    chains_with_class = np.where(chains_arr[:, k] > 0)[0] if k < chains_arr.shape[1] else []
                    if len(chains_with_class) > 0:
                        c = chains_with_class[0]  # Use first chain (classes typically in one chain)
                        if c in sn.visits and sn.visits[c] is not None:
                            visits_c = np.asarray(sn.visits[c])
                            stationToStateful = np.asarray(sn.stationToStateful).flatten() if hasattr(sn, 'stationToStateful') else None
                            for i in range(M):
                                # visits_c is indexed by stateful node, not station
                                isf = int(stationToStateful[i]) if stationToStateful is not None and i < len(stationToStateful) else i
                                if isf < visits_c.shape[0] and k < visits_c.shape[1]:
                                    if visits_c[isf, k] == 0:
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

        # Compute ResidT using proper visit ratios from network structure
        # (after zeroing, so zeroed entries stay zero)
        if sn is not None and sn.visits:
            WN = sn_get_residt_from_respt(sn, RN, None)
        else:
            WN = RN.copy()

        # Convert station-level to node-level results
        # Initialize node-level arrays
        QNn = np.zeros((I, K))
        UNn = np.zeros((I, K))
        RNn = np.zeros((I, K))
        WNn = np.zeros((I, K))

        # Map station metrics to node metrics
        if sn is not None and hasattr(sn, 'stationToNode'):
            stationToNode = np.asarray(sn.stationToNode).flatten()
            for ist in range(M):
                if ist < len(stationToNode):
                    ind = int(stationToNode[ist])
                    if ind >= 0 and ind < I:
                        QNn[ind, :] = QN[ist, :]
                        UNn[ind, :] = UN[ist, :]
                        RNn[ind, :] = RN[ist, :]
                        WNn[ind, :] = WN[ist, :]
        else:
            # No mapping - assume stations are nodes
            for ist in range(min(M, I)):
                QNn[ist, :] = QN[ist, :]
                UNn[ist, :] = UN[ist, :]
                RNn[ist, :] = RN[ist, :]
                WNn[ist, :] = WN[ist, :]

        # Compute node-level throughputs and arrival rates
        # This properly handles cache hit/miss class throughputs using actual probabilities
        ANn = sn_get_node_arvr_from_tput(sn, TN, TN)
        TNn = sn_get_node_tput_from_tput(sn, TN, TN, ANn)

        # cache node post-processing: hit/miss throughput = downstream ClassSwitch arrival rate; the requesting class's own throughput at the cache is 0.
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            for cache_ind in range(I):
                if sn.nodetype[cache_ind] != NodeType.CACHE:
                    continue
                if sn.nodeparam is None or cache_ind not in sn.nodeparam:
                    continue

                node_param = sn.nodeparam[cache_ind]
                if not hasattr(node_param, 'hitclass') or not hasattr(node_param, 'missclass'):
                    continue

                hitclass = np.atleast_1d(node_param.hitclass).flatten()
                missclass = np.atleast_1d(node_param.missclass).flatten()

                # Find ClassSwitch node connected to this cache
                cs_ind = -1
                for jnd in range(I):
                    if sn.nodetype[jnd] == NodeType.CLASSSWITCH:
                        if hasattr(sn, 'connmatrix') and sn.connmatrix is not None:
                            if sn.connmatrix[cache_ind, jnd] > 0:
                                cs_ind = jnd
                                break

                # For each requesting class that has hit/miss classes
                for orig_class in range(len(hitclass)):
                    h = int(hitclass[orig_class]) if orig_class < len(hitclass) else -1
                    m = int(missclass[orig_class]) if orig_class < len(missclass) else -1

                    if h >= 0 and h < K and m >= 0 and m < K:
                        # Set throughput of requesting class at cache to 0
                        # (jobs leave as hit/miss classes)
                        TNn[cache_ind, orig_class] = 0.0

                        # Set throughput of hit/miss classes at cache
                        # equals the arrival rate at ClassSwitch for those classes
                        if cs_ind >= 0 and cs_ind < I:
                            if h < ANn.shape[1]:
                                TNn[cache_ind, h] = ANn[cs_ind, h]
                            if m < ANn.shape[1]:
                                TNn[cache_ind, m] = ANn[cs_ind, m]

        # Build table with station-level results (like MATLAB's getAvgTable)
        # MATLAB iterates over stations (M), not all nodes (I)
        rows = []
        nodenames = list(sn.nodenames) if sn is not None and hasattr(sn, 'nodenames') and sn.nodenames else []
        classnames = list(sn.classnames) if sn is not None and hasattr(sn, 'classnames') and sn.classnames else []

        # Get station-to-node mapping
        stationToNode = None
        if sn is not None and hasattr(sn, 'stationToNode'):
            stationToNode = np.asarray(sn.stationToNode).flatten()

        for ist in range(M):
            # Get node index for this station
            if stationToNode is not None and ist < len(stationToNode):
                ind = int(stationToNode[ist])
            else:
                ind = ist

            for r in range(K):
                node_name = nodenames[ind] if ind < len(nodenames) else f'Station{ist}'
                class_name = classnames[r] if r < len(classnames) else f'Class{r}'

                rows.append({
                    'Station': node_name,
                    'JobClass': class_name,
                    'QLen': QN[ist, r],
                    'Util': UN[ist, r],
                    'RespT': RN[ist, r],
                    'ResidT': WN[ist, r],
                    'ArvR': ANn[ind, r] if ind < ANn.shape[0] else 0.0,
                    'Tput': TN[ist, r],
                })

        df = pd.DataFrame(rows)

        # Filter out all-zero rows (MATLAB excludes nodes with no metrics)
        numeric_cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        tokeep = ~(df[numeric_cols] <= 0.0).all(axis=1)
        df = df.loc[tokeep].reset_index(drop=True)

        if not self._table_silent:
            print(df.to_string(index=False))

        from ...indexed_table import IndexedTable
        return IndexedTable(df)

    # =========================================================================
    # Individual Metric Accessors
    # =========================================================================

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths (M x K)."""
        self._assert_not_chain_model('getAvgQLen')
        if self._result is None:
            self._ensureAvgResults()
        return self._result.Q.copy()

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K)."""
        self._assert_not_chain_model('getAvgUtil')
        if self._result is None:
            self._ensureAvgResults()
        return self._result.U.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K)."""
        self._assert_not_chain_model('getAvgRespT')
        if self._result is None:
            self._ensureAvgResults()
        return self._result.R.copy()

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times (M x K).

        Residence time is computed from response time using visit ratios:
        WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]
        """
        self._assert_not_chain_model('getAvgResidT')
        if self._result is None:
            self._ensureAvgResults()

        # Compute ResidT using proper visit ratios from network structure
        if self._sn is not None and self._sn.visits:
            return sn_get_residt_from_respt(self._sn, self._result.R, None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            return self._result.R.copy()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (M x K)."""
        self._assert_not_chain_model('getAvgWaitT')
        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R.copy()
        # W = R - S where S is service time (1/rate)
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
        """Get average throughputs (M x K)."""
        self._assert_not_chain_model('getAvgTput')
        if self._result is None:
            self._ensureAvgResults()
        return self._result.T.copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K)."""
        self._assert_not_chain_model('getAvgArvR')
        if self._result is None:
            self._ensureAvgResults()
        TN = self._result.T.copy()
        return sn_get_arvr_from_tput(self._sn, TN, TN)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get chain-level system response times (nchains,).

        Uses the shared chain-based algorithm (a faithful port of MATLAB
        @NetworkSolver/getAvgSys.m): open chains sum alpha-weighted class
        residence times, closed chains apply Little's law nJobsChain/XNchain.
        """
        self._assert_not_chain_model('getAvgSysRespT')
        CN, _ = self._computeChainMetrics()
        return CN

    def getAvgSysTput(self) -> np.ndarray:
        """Get chain-level system (carried) throughputs (nchains,).

        Matches MATLAB/JAR getAvgSys: the throughput of completing classes
        routed back into the chain reference station (carried rate), not the
        offered/source arrival rate.
        """
        self._assert_not_chain_model('getAvgSysTput')
        _, XN = self._computeChainMetrics()
        return XN

    # =========================================================================
    # CTMC-Specific Methods
    # =========================================================================

    def getStateSpace(self):
        """Get the enumerated state space.

        Returns:
            tuple: (stateSpace, localStateSpace) where stateSpace is the global
                state matrix and localStateSpace is a list of per-station state arrays.
                For FCFS queues with multiple servers, localStateSpace includes
                buffer and phase columns (matching MATLAB's nodeStateSpace format).
        """
        if self._result is None:
            self._ensureAvgResults()
        space = self._result.space.copy() if self._result.space is not None else np.array([])

        if self.isChainSolver():
            # Chain mode: one component, so there is no per-station slicing.
            return space, [space]

        # Generate localStateSpace - slice state space by station using column ranges
        localStateSpace = []
        if self._sn is not None and space.size > 0:
            nstations = self._sn.nstations

            # Use station_col_ranges if available (proper column structure)
            if hasattr(self._result, 'station_col_ranges') and self._result.station_col_ranges is not None:
                for ist in range(nstations):
                    if ist < len(self._result.station_col_ranges):
                        start_col, end_col = self._result.station_col_ranges[ist]
                        if end_col > start_col and end_col <= space.shape[1]:
                            localStateSpace.append(space[:, start_col:end_col])
                        else:
                            # Empty range for this station (e.g., Source/Sink)
                            localStateSpace.append(np.array([]).reshape(space.shape[0], 0))
                    else:
                        localStateSpace.append(np.array([]).reshape(space.shape[0], 0))
            else:
                # Fallback: assume one column per (station, class) pair (old behavior)
                nclasses = self._sn.nclasses
                for i in range(nstations):
                    station_cols = []
                    for r in range(nclasses):
                        col_idx = i * nclasses + r
                        if col_idx < space.shape[1]:
                            station_cols.append(space[:, col_idx:col_idx+1])
                    if station_cols:
                        localStateSpace.append(np.hstack(station_cols))
                    else:
                        localStateSpace.append(np.array([]))

        return space, localStateSpace

    def getSteadyState(self) -> np.ndarray:
        """Get the steady-state probability distribution."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.pi.copy() if self._result.pi is not None else np.array([])

    def getInfGen(self) -> np.ndarray:
        """Get the infinitesimal generator matrix."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.infgen.copy() if self._result.infgen is not None else np.array([])

    # =========================================================================
    # CDF and Percentile Methods
    # =========================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Response-time distribution by tagged-chain analysis.

        One job of each chain is tagged, the tagged model is solved with its
        event filtration kept, and for each station the arrival and departure
        events OF THE TAGGED JOB split the generator into TWO maps:

            A = map_normalize(Q - A1, A1)     A1: tagged job arrives at station
            D = map_normalize(Q - D1, D1)     D1: tagged job departs station
            pie = map_pie(A)                  the state seen ON ARRIVAL
            F(t) = 1 - pie expm(D.D0 t) 1

        The two maps are not interchangeable: pie must come from the ARRIVAL
        map, and D0 from the DEPARTURE one.

        THIS REPLACED AN EXPONENTIAL FIT that returned 1 - exp(-t/R) from the
        mean response time, with no tagging and no filtration, and was therefore
        exact only for an M/M/1.

        Reference: matlab/src/solvers/CTMC/@SolverCTMC/getCdfRespT.m.

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import cdf_respt_via_jar
            return cdf_respt_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import cdf_respt_via_cpp, cpp_unsupported
            if R is not None:
                # R overrides the mean this native getter builds its exponential
                # approximation from; the C++ integrates the tagged chain instead
                # and has no mean to override, so a supplied R would be ignored.
                cpp_unsupported(
                    self, 'getCdfRespT(R=...)',
                    "the C++ integrates the tagged chain rather than fitting an exponential to a "
                    "mean response time, so there is no R for it to take")
            return cdf_respt_via_cpp(self)
        if self._result is None:
            self._ensureAvgResults()

        return self._taggedCdfRespT()

    def _taggedCdfRespT(self) -> List[Dict]:
        """The tagged-chain response-time law. See getCdfRespT."""
        import copy as _copy
        from scipy.linalg import expm as _expm
        from ...api.io.model_adapter import tag_chain
        from ...api.mam.map_analysis import map_normalize, map_pie
        from ...constants import EventType
        from ...lang.sync import refresh_sync

        sn = self._sn if self._sn is not None else self.model.get_struct()
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        if np.any(np.isinf(njobs)):
            raise RuntimeError(
                "getCdfRespT is presently supported only for closed models.")

        classes = self.model.get_classes()
        nchains = int(sn.nchains)
        RD: List[Dict] = []

        class _Chain(object):
            """The minimal chain shape tag_chain reads: its class objects."""
            def __init__(self, cls):
                self.classes = cls

        for ch in range(nchains):
            inchain = [int(x) for x in np.asarray(sn.inchain[ch]).ravel()]
            tagged_src_idx = None
            for r in inchain:
                if njobs[r] > 0:
                    tagged_src_idx = r
                    break
            if tagged_src_idx is None:
                continue

            tagged = tag_chain(self.model,
                               _Chain([classes[r] for r in inchain]),
                               classes[tagged_src_idx])
            tsolver = SolverCTMC(tagged.model, self.options)
            Q, filt = tsolver.getGenerator()
            Q = np.asarray(Q, dtype=float)
            tsn = tagged.model.get_struct()
            # sn.sync is never stored by the native struct: the CTMC handler
            # builds it locally with refresh_sync and keeps it on the stack, so
            # the same call is what guarantees this ordering matches Dfilt's.
            sync = refresh_sync(tsn)
            if sync is None or filt is None or len(filt) == 0:
                raise RuntimeError(
                    "getCdfRespT needs the event filtration of the tagged chain, which this "
                    "model did not produce; the response-time law cannot be computed without it")

            tagged_cls = tagged.tagged_job._index \
                if hasattr(tagged.tagged_job, '_index') else len(tagged.model.get_classes()) - 1
            node_to_station = np.asarray(tsn.nodeToStation).ravel()

            for ist in range(int(tsn.nstations)):
                A1 = np.zeros_like(Q)
                D1 = np.zeros_like(Q)
                for v, ev in enumerate(sync):
                    if v >= len(filt) or filt[v] is None:
                        continue
                    Fv = np.asarray(filt[v], dtype=float)
                    pas = getattr(ev, 'passive', None)
                    act = getattr(ev, 'active', None)
                    if (pas is not None and pas.event == EventType.ARV
                            and pas.job_class == tagged_cls
                            and 0 <= pas.node < len(node_to_station)
                            and int(node_to_station[pas.node]) == ist):
                        A1 = A1 + Fv
                    if (act is not None and act.event == EventType.DEP
                            and act.job_class == tagged_cls
                            and 0 <= act.node < len(node_to_station)
                            and int(node_to_station[act.node]) == ist):
                        D1 = D1 + Fv
                if not np.any(A1) or not np.any(D1):
                    continue

                A0n, A1n = map_normalize(Q - A1, A1)
                pie = np.asarray(map_pie(A0n, A1n), dtype=float).ravel()
                D0, _ = map_normalize(Q - D1, D1)

                nz = np.abs(Q[Q != 0])
                nz = nz[nz > 1e-8]
                if nz.size == 0:
                    continue
                intervals = 100000
                T = abs(100.0 / nz.min())
                dT = T / intervals
                # One matrix exponential, then propagate: the reference
                # recomputes expm(D0*t) at each of the 100001 grid points, which
                # is the same answer at a cost linear in the grid.
                E = _expm(D0 * dT)
                ones = np.ones(D0.shape[0])
                v = pie.copy()
                tvals = []
                Fvals = []
                for k in range(intervals + 1):
                    if k > 0:
                        v = v.dot(E)
                    Fk = min(1.0, max(0.0, 1.0 - float(v.dot(ones))))
                    tvals.append(k * dT)
                    Fvals.append(Fk)
                    if Fk > 1.0 - 1e-3:
                        break

                RD.append({
                    'station': ist + 1,
                    'class': tagged_src_idx + 1,
                    't': np.array(tvals),
                    'p': np.array(Fvals),
                })

        return RD

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """
        Extract percentiles from response time distribution.

        Args:
            percentiles: List of percentiles (0-100). Default: [10, 25, 50, 75, 90, 95, 99]
            jobclass: Optional class filter (1-based)

        Returns:
            Tuple of (percentile_list, percentile_table)
        """
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90, 95, 99]

        percentiles = np.asarray(percentiles)
        percentiles = np.clip(percentiles, 0.01, 99.99)
        percentiles_normalized = percentiles / 100.0

        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R
        nstations, nclasses = R.shape

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]

        for i in range(nstations):
            for r in range(nclasses):
                if jobclass is not None and (r + 1) != jobclass:
                    continue

                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                lambda_rate = 1.0 / mean_resp_t
                perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                PercRT.append({
                    'station': i + 1,
                    'class': r + 1,
                    'percentiles': percentiles.tolist(),
                    'values': perc_values.tolist(),
                })

                row_data = {
                    'Station': self.station_names[i] if i < len(self.station_names) else f'Station{i}',
                    'Class': self.class_names[r] if r < len(self.class_names) else f'Class{r}',
                }
                for perc_col, perc_val in zip(perc_col_names, perc_values):
                    row_data[perc_col] = perc_val
                rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    # =========================================================================
    # Probability Methods
    # =========================================================================

    def _build_ssq(self):
        """
        Build the aggregated state space SSq (per-class job counts) from
        the detailed state space SS using State.toMarginal.

        This matches MATLAB's ctmc_ssg approach: for each state in SS,
        extract nir (jobs per class) at each stateful node using toMarginal.

        Returns:
            np.ndarray of shape (nstates, nstations * nclasses)
        """
        from ...api.state.marginal import toMarginal

        sn = self._sn
        SS = self._result.space
        station_col_ranges = self._result.station_col_ranges

        nstates = SS.shape[0]
        nstations = sn.nstations
        nclasses = sn.nclasses

        SSq = np.zeros((nstates, nstations * nclasses))

        for s in range(nstates):
            for ind in range(sn.nnodes):
                if sn.isstateful[ind]:
                    isf = sn.nodeToStateful[ind]
                    ist = sn.nodeToStation[ind]
                    if ist < 0:
                        continue

                    # Extract state portion for this node from SS row
                    if station_col_ranges is not None and ist < len(station_col_ranges):
                        start_col, end_col = station_col_ranges[ist]
                        if start_col < end_col:
                            state_portion = SS[s, start_col:end_col]
                        else:
                            state_portion = np.zeros(nclasses)
                    else:
                        state_portion = np.zeros(nclasses)

                    # Use toMarginal to get per-class job counts
                    try:
                        _, nir, _, _ = toMarginal(sn, ind, state_portion)
                        nir = np.atleast_1d(nir).flatten()
                        # Store in SSq at the right position
                        col_start = ist * nclasses
                        col_end = col_start + nclasses
                        if len(nir) >= nclasses:
                            SSq[s, col_start:col_end] = nir[:nclasses]
                        else:
                            SSq[s, col_start:col_start + len(nir)] = nir
                    except Exception:
                        pass

        return SSq

    def _station_class_counts(self, ist0: int) -> list:
        """Per-class job counts at station `ist0` (0-based) in the model's state."""
        from ...api.state.marginal import toMarginal

        sn = self._sn if self._sn is not None else self.model.getStruct()
        ind = int(np.asarray(sn.stationToNode).flatten()[ist0])
        isf = int(np.asarray(sn.nodeToStateful).flatten()[ind])
        state_i = np.asarray(sn.state[isf], dtype=float).flatten()
        _, nir, _, _ = toMarginal(sn, ind, state_i)
        nir = np.asarray(nir).reshape(-1)[:int(sn.nclasses)]
        return [int(round(v)) for v in nir]

    def _system_class_counts(self) -> list:
        """Per-class job counts at EVERY station, station-major: the system
        state the joint getters ask about, in the shape `prob_via_jar` sends."""
        sn = self._sn if self._sn is not None else self.model.getStruct()
        return [self._station_class_counts(i) for i in range(int(sn.nstations))]

    def getProbAggr(self, ist) -> float:
        """
        Get probability of a specific per-class job distribution at a station.

        Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for the state
        that was set via setState() on the station.

        Matches MATLAB: solver_ctmc_margaggr.m

        Args:
            ist: Station index (0-based) or node object

        Returns:
            Probability that station ist is in the specified state (scalar).
        """
        self._assert_phasetype_states('getProbAggr')
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            # model.json carries no per-station initial state, so a delegated
            # query is answered at the JAR's default initialization unless the
            # cell is named explicitly (as SolverFLD.getProbAggr already does).
            station0 = int(ist)
            return prob_via_jar(self, 'prob-aggr', ist=station0, kind='scalar',
                                state=self._station_class_counts(station0))
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            # `-s ctmc -a prob` reports every station's marginal at the model's
            # declared state, which the wire now carries, so the selection here
            # is an index into the engine's answer and not a computation of it.
            station0 = ist if isinstance(ist, (int, np.integer)) else ist.get_station_index0()
            p = prob_aggr_via_cpp(self)['probAggr']
            if not (0 <= int(station0) < len(p)):
                raise ValueError("station index %r is outside 0..%d" % (station0, len(p) - 1))
            return float(p[int(station0)])
        from ...api.state.marginal import toMarginal

        # Convert node object to index if needed (like MATLAB)
        if not isinstance(ist, (int, np.integer)):
            ist = ist.get_station_index0()

        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        SS = self._result.space

        if pi is None or len(pi) == 0 or SS is None or len(SS) == 0:
            return 0.0

        if self._sn is None or self._sn.state is None:
            return 0.0

        sn = self._sn
        state = sn.state
        station_col_ranges = self._result.station_col_ranges
        space_aggr = self._result.space_aggr

        # Clamp small negative values to zero (matching MATLAB: pi(pi<Zero)=0)
        pi = pi.copy()
        pi[pi < 1e-14] = 0.0

        nclasses = sn.nclasses
        nstations = sn.nstations

        # Build the per-class marginal query vector (nivec) for each station
        # Then compare against the aggregated state space
        Pnir = np.zeros(nstations)

        for ind in range(sn.nnodes):
            if not sn.isstateful[ind]:
                continue
            isf = int(sn.nodeToStateful[ind])
            ist_node = int(sn.nodeToStation[ind])
            if ist_node < 0:
                continue

            # Get query state for this node
            if isf >= len(state) or state[isf] is None:
                continue

            state_isf = np.atleast_1d(state[isf]).flatten()

            # Use toMarginal to convert state to per-class job counts
            try:
                _, nivec, _, _ = toMarginal(sn, ind, state_isf)
                nivec = np.atleast_1d(nivec).flatten()
            except Exception:
                # Fallback: use state directly as per-class counts
                nivec = state_isf[:nclasses] if len(state_isf) >= nclasses else state_isf

            # Get column range for this station
            if station_col_ranges is not None and ist_node < len(station_col_ranges):
                col_start, col_end = station_col_ranges[ist_node]
            else:
                continue

            # Sum probabilities for matching marginal states
            Pnir_ist = 0.0
            for s_idx in range(SS.shape[0]):
                ss_portion = SS[s_idx, col_start:col_end]

                # Use toMarginal to convert SS row portion to per-class counts
                try:
                    _, sivec, _, _ = toMarginal(sn, ind, ss_portion)
                    sivec = np.atleast_1d(sivec).flatten()
                except Exception:
                    # Fallback: use raw columns as per-class counts
                    sivec = ss_portion[:nclasses] if len(ss_portion) >= nclasses else ss_portion

                if len(sivec) == len(nivec) and np.all(sivec == nivec):
                    Pnir_ist += pi[s_idx]

            Pnir[ist_node] = Pnir_ist

        return Pnir[ist]

    def getProbSysAggr(self) -> float:
        """
        Get probability of the entire system being in the specified aggregated state.

        Returns the joint probability of the system being in the aggregated
        state configuration set via setState() on all stations.

        Matches MATLAB: solver_ctmc_jointaggr.m

        Returns:
            float: Joint probability of the system state.
        """
        self._assert_phasetype_states('getProbSysAggr')
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            # model.json carries no per-station initial state, so the whole
            # system state has to be named or the JAR answers about ITS default
            # initialization -- every closed job at its reference station, which
            # is a different question and not a numerically close one.
            return prob_via_jar(self, 'prob-sys-aggr', kind='scalar',
                                state=self._system_class_counts())
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            return float(prob_aggr_via_cpp(self)['probSysAggr'])
        from ...api.state.marginal import toMarginal

        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        SS = self._result.space

        if pi is None or len(pi) == 0 or SS is None or len(SS) == 0:
            return 0.0

        if self._sn is None or self._sn.state is None:
            return 0.0

        sn = self._sn
        state = sn.state
        station_col_ranges = self._result.station_col_ranges

        # Clamp small negative values to zero
        pi = pi.copy()
        pi[pi < 1e-14] = 0.0

        nclasses = sn.nclasses
        nstations = sn.nstations

        # Build target nvec: per-class job counts at each station from sn.state
        target_nir = []  # list of arrays, one per station
        for i in range(nstations):
            isf = int(sn.stationToStateful[i]) if hasattr(sn, 'stationToStateful') else i
            node_idx = int(sn.stationToNode[i]) if hasattr(sn, 'stationToNode') else i
            if isf < len(state) and state[isf] is not None:
                state_isf = np.atleast_1d(state[isf]).flatten()
                try:
                    _, nir, _, _ = toMarginal(sn, node_idx, state_isf)
                    nir = np.atleast_1d(nir).flatten()[:nclasses]
                except Exception:
                    nir = state_isf[:nclasses]
            else:
                nir = np.zeros(nclasses)
            target_nir.append(nir)

        # For each SS row, extract per-station nir and check if all stations match
        prob = 0.0
        for s_idx in range(SS.shape[0]):
            all_match = True
            for i in range(nstations):
                node_idx = int(sn.stationToNode[i]) if hasattr(sn, 'stationToNode') else i
                # Extract state portion for this station
                if station_col_ranges is not None and i < len(station_col_ranges):
                    start_col, end_col = station_col_ranges[i]
                    state_portion = SS[s_idx, start_col:end_col]
                else:
                    state_portion = np.zeros(nclasses)

                try:
                    _, nir, _, _ = toMarginal(sn, node_idx, state_portion)
                    nir = np.atleast_1d(nir).flatten()[:nclasses]
                except Exception:
                    nir = state_portion[:nclasses]

                if not np.allclose(nir, target_nir[i]):
                    all_match = False
                    break

            if all_match:
                prob += pi[s_idx]

        return prob

    def getProbSys(self) -> float:
        """
        Get joint probability for the detailed (non-aggregated) system state.

        Matches MATLAB: solver_ctmc_joint.m

        In chain mode this returns the stationary vector of the user-supplied
        chain, one entry per state of the chain state space.

        Returns:
            float: Joint probability of the detailed system state.
        """
        if self.isChainSolver():
            self._ensureAvgResults()
            return self._result.pi.copy()
        self._assert_phasetype_states('getProbSys')
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-sys', kind='scalar',
                                state=self._system_class_counts())
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_sys_via_cpp
            return prob_sys_via_cpp(self)
        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        SS = self._result.space

        if pi is None or len(pi) == 0 or SS is None or len(SS) == 0:
            return 0.0

        if self._sn is None or self._sn.state is None:
            return 0.0

        sn = self._sn
        state = sn.state
        station_col_ranges = self._result.station_col_ranges

        # Clamp small negative values to zero
        pi = pi.copy()
        pi[pi < 1e-14] = 0.0

        # Build statevec: full detailed state vector (matching MATLAB solver_ctmc_joint)
        statevec = []
        stateful_indices = []
        for ind in range(sn.nnodes):
            if sn.isstateful[ind]:
                stateful_indices.append(sn.nodeToStateful[ind])

        for ind in range(sn.nnodes):
            if not sn.isstateful[ind]:
                continue
            isf = sn.nodeToStateful[ind]

            # Get the space size for this stateful node
            if station_col_ranges is not None:
                ist = sn.nodeToStation[ind]
                if ist >= 0 and ist < len(station_col_ranges):
                    s, e = station_col_ranges[ist]
                    space_width = e - s
                else:
                    space_width = 0
            elif hasattr(sn, 'space') and sn.space is not None and isf < len(sn.space) and sn.space[isf] is not None:
                space_width = sn.space[isf].shape[1] if sn.space[isf].ndim > 1 else len(sn.space[isf])
            else:
                space_width = 0

            if isf < len(state) and state[isf] is not None:
                state_isf = np.atleast_1d(state[isf]).flatten()
                # Zero-pad on the left to match space width (like MATLAB)
                if len(state_isf) < space_width:
                    state_isf = np.concatenate([np.zeros(space_width - len(state_isf)), state_isf])
                statevec.extend(state_isf[:space_width].tolist())
            else:
                statevec.extend([0.0] * space_width)

        statevec = np.array(statevec)

        # Find matching row in SS (MATLAB: pi(findrows(SS, statevec)))
        prob = 0.0
        for s_idx in range(SS.shape[0]):
            if len(statevec) <= SS.shape[1] and np.allclose(SS[s_idx, :len(statevec)], statevec):
                prob += pi[s_idx]

        return prob


    def _assert_phasetype_states(self, what: str) -> None:
        """Refuse a query whose answer is a per-state probability under an ME.

        A matrix-exponential service embeds in the generator with negative
        off-diagonal entries, so the stationary vector is a SIGNED measure: only
        its aggregates over each phase block are probabilities. Mean measures
        stay exact (they are linear in that vector), but a per-state or
        transient answer is not a probability at all, and uniformization -- a
        Poisson mixture of powers of I + Q/lambda -- diverges on a signed
        generator. Such queries are refused rather than answered with a number
        that looks like a probability. See sn.isph and _kb/04-networkstruct.md.
        """
        sn = self.model.getStruct()
        isph = getattr(sn, 'isph', None)
        if isph is not None and not bool(np.all(np.asarray(isph, dtype=bool))):
            raise ValueError(
                '%s is unavailable: the model has a matrix-exponential (ME) '
                'service or arrival process, so the stationary vector of the '
                'generator is a signed measure and per-state probabilities and '
                'uniformization-based transients do not exist. Mean measures '
                '(getAvg, getAvgTable) remain exact.' % what)

    def getProb(self, station=None) -> float:
        """Get probability for the detailed state at station.

        Returns the probability that the station is in the state that was set
        via setState(). This includes phase information from service distributions.

        Matches MATLAB: solver_ctmc_marg.m

        In chain mode the argument is a state of the user-supplied chain: a row
        of its state space, or a 1-based state index when the chain carries none.

        Args:
            station: Station index (0-based) or node object. If None, returns steady-state.

        Returns:
            float: Probability that station is in the specified detailed state.
        """
        if self.isChainSolver():
            self._ensureAvgResults()
            pi = self._result.pi
            if station is None:
                return pi.copy()
            user_space = self._chain_matrix.stateSpace if self.isDiscreteChain() else self._chain_process.stateSpace
            if user_space is None or np.asarray(user_space).size == 0:
                idx = np.asarray(station).flatten()
                if idx.size != 1 or idx[0] != round(float(idx[0])) or not (1 <= idx[0] <= len(pi)):
                    raise RuntimeError(
                        f"The chain carries no state space, so getProb requires a state index in 1..{len(pi)}.")
                return float(pi[int(idx[0]) - 1])
            user_space = np.atleast_2d(np.asarray(user_space, dtype=np.float64))
            row = np.asarray(station, dtype=np.float64).flatten()
            matches = np.flatnonzero(np.all(user_space == row, axis=1))
            if matches.size == 0:
                raise RuntimeError("The requested state is not in the chain state space.")
            return float(pi[matches[0]])
        self._assert_phasetype_states('getProb')
        if station is None:
            if self._result is None:
                self._ensureAvgResults()
            return self.getSteadyState()

        # Convert node object to index if needed (like MATLAB)
        if not isinstance(station, (int, np.integer)):
            station = station.get_station_index0()

        # lang=java has no native state space, so _result.station_col_ranges (native-only) can't be rebuilt there; routed like getProbAggr/getProbSys.
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob', ist=station, kind='scalar')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            # The DETAILED marginal, phases and buffer arrangement included; the
            # aggregate one is `ProbAggr` in the same payload, off one solve.
            p = prob_aggr_via_cpp(self)['prob']
            if not (0 <= int(station) < len(p)):
                raise ValueError("station index %r is outside 0..%d" % (station, len(p) - 1))
            return float(p[int(station)])

        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        SS = self._result.space

        if pi is None or len(pi) == 0 or SS is None or len(SS) == 0:
            return 0.0

        if self._sn is None or self._sn.state is None:
            return 0.0

        sn = self._sn
        state = sn.state
        station_col_ranges = self._result.station_col_ranges

        # Clamp small negative values to zero (matching MATLAB: pi(pi<Zero)=0)
        pi = pi.copy()
        pi[pi < 1e-14] = 0.0

        # Compute probability for each station
        nstations = sn.nstations
        Pnir = np.zeros(nstations)

        for ind in range(sn.nnodes):
            if not sn.isstateful[ind]:
                continue
            isf = int(sn.nodeToStateful[ind])
            ist_node = int(sn.nodeToStation[ind])
            if ist_node < 0:
                continue

            # Get query state for this node
            if isf >= len(state) or state[isf] is None:
                continue

            state_isf = np.atleast_1d(state[isf]).flatten()

            # Get column range for this station
            if station_col_ranges is not None and ist_node < len(station_col_ranges):
                col_start, col_end = station_col_ranges[ist_node]
            else:
                continue

            # Pad state with zeros if needed (matches MATLAB: state_i = [zeros(...), state{isf}])
            required_length = col_end - col_start
            if len(state_isf) < required_length:
                state_i = np.zeros(required_length)
                state_i[required_length - len(state_isf):] = state_isf
            else:
                state_i = state_isf[:required_length]

            # Sum probabilities for matching detailed states (exact match including phases)
            Pnir_ist = 0.0
            for s_idx in range(SS.shape[0]):
                ss_portion = SS[s_idx, col_start:col_end]
                if len(ss_portion) == len(state_i) and np.allclose(ss_portion, state_i, atol=1e-10):
                    Pnir_ist += pi[s_idx]

            Pnir[ist_node] = Pnir_ist

        return Pnir[station]

    def getGenerator(self):
        """Get the infinitesimal generator matrix and event filters.

        Returns:
            tuple: (infGen, eventFilt) where infGen is the infinitesimal generator
                matrix and eventFilt is a dictionary mapping event types to sparse matrices
        """
        infGen = self.getInfGen()

        # Generate event filters from the stored event information
        eventFilt = []
        if self._result is not None and hasattr(self._result, 'eventFilt') and self._result.eventFilt is not None:
            eventFilt = self._result.eventFilt
        elif self._result is not None and hasattr(self._result, 'infgen') and self._result.infgen is not None:
            # Return empty list if eventFilt not available
            eventFilt = []

        return infGen, eventFilt

    def getAsymptoticVariance(self, f):
        """
        The asymptotic variance of the time-average of a reward ``f`` along a
        sample path of this model's CTMC.

        WHAT IT IS FOR. A simulation estimate of a steady-state mean has a
        standard error that shrinks like ``sqrt(sigma^2/t)``, where sigma^2 is
        NOT the stationary variance of f but its ASYMPTOTIC variance, which also
        carries the autocorrelation of the path. That number is what says how
        long a run has to be, and :func:`sim_runlength` turns it into a run
        length for a target precision. It cannot be guessed from the stationary
        variance: on M/M/1 the two differ by a factor that blows up like
        ``(1-rho)^-2``.

        Args:
            f: one reward value per CTMC state, in the state order
                :meth:`getGenerator` returns, or a callable applied to each row
                of the state space

        Returns:
            The dict of :func:`sim_asymvar_ctmc`: ``mean``, ``variance``,
            ``asymptoticVariance`` and the deviation vector.

        References:
            W. Whitt (1989). Planning queueing simulations. Management Science
            35(11), 1341-1366.
        """
        from ...api.sim.runlength import sim_asymvar_ctmc
        from ...api.mc.ctmc import ctmc_solve
        infGen = np.asarray(self.getInfGen(), dtype=float)
        if hasattr(infGen, 'toarray'):
            infGen = infGen.toarray()
        n = infGen.shape[0]
        if callable(f):
            # getStateSpace returns (global, per-station); the reward is a
            # function of the GLOBAL state, which is the first of the two.
            space = self.getStateSpace()
            if isinstance(space, tuple):
                space = space[0]
            space = np.asarray(space)
            if space.shape[0] != n:
                raise RuntimeError('the state space has %d rows but the generator is %dx%d; pass '
                                   'the reward as a vector instead' % (space.shape[0], n, n))
            fvec = np.array([float(f(space[i, :])) for i in range(n)], dtype=float)
        else:
            fvec = np.asarray(f, dtype=float).ravel()
            if fvec.size != n:
                raise RuntimeError('the reward vector has %d entries but the generator is %dx%d'
                                   % (fvec.size, n, n))
        pi_ss = np.asarray(ctmc_solve(infGen), dtype=float).ravel()
        return sim_asymvar_ctmc(infGen, fvec, pi_ss)

    def getStartRate(self) -> np.ndarray:
        """(nstations x nclasses) rate at which a class-r job BEGINS or RESUMES
        holding a server at station i, i.e. pi*F*e over the START filtration.

        At a lossless station with no in-service abandonment

            getStartRate == getAvgTput + getPreemptRate

        because every job starts service once per entry into a server and every
        preemption is followed by exactly one later resume or restart. At a
        non-preemptive station this collapses to startRate == throughput.

        An accessor, not a MetricType: it adds no getAvgTable column.
        """
        if self._result is None or getattr(self._result, 'startRate', None) is None:
            self._ensureAvgResults()
        rate = getattr(self._result, 'startRate', None)
        if rate is None:
            raise RuntimeError("This solver run produced no START rates.")
        return np.asarray(rate)

    def getPreemptRate(self) -> np.ndarray:
        """(nstations x nclasses) rate at which a class-r job HOLDING A SERVER
        at station i is pushed back into the buffer. Identically zero at a
        non-preemptive station; preempt-resume and preempt-independent stations
        report the SAME rate, since which phase the displaced job resumes in is
        not a property of how often it is displaced."""
        if self._result is None or getattr(self._result, 'preemptRate', None) is None:
            self._ensureAvgResults()
        rate = getattr(self._result, 'preemptRate', None)
        if rate is None:
            raise RuntimeError("This solver run produced no PREEMPT rates.")
        return np.asarray(rate)

    def getEventFiltration(self, event_type):
        """Filtration of a DERIVED event type, indexed [station][class]: the
        (s,ns) entry is the rate at which the transition s -> ns carries one
        such event at that station for that class.

        EVENT_TYPE must be EventType.START or EventType.PREEMPT. The two are not
        synchronizations: they are tags on the ARV and DEP arcs that cause them,
        so they are NOT part of the event filtration getGenerator returns (which
        pairs one-to-one with sn.sync and is summed as D1) and are kept here.
        """
        from ...constants import EventType
        if event_type not in (EventType.START, EventType.PREEMPT):
            raise ValueError(
                "getEventFiltration serves the derived events only (START, PREEMPT); "
                "%s is a synchronization and its filtration is the one getGenerator returns."
                % str(event_type))
        if self._result is None or getattr(self._result, 'startFilt', None) is None:
            self._ensureAvgResults()
        filt = (getattr(self._result, 'startFilt', None) if event_type == EventType.START
                else getattr(self._result, 'preemptFilt', None))
        if filt is None:
            raise RuntimeError("This model produced no derived event filtration.")
        return filt

    def getStateSpaceAggr(self) -> np.ndarray:
        """Get aggregated state space (jobs per station per class).

        Returns:
            Array of shape (nstates, nstations * nclasses) where column
            (ist * nclasses + k) = jobs of class k at station ist (0-indexed)
        """
        if self._result is None:
            self._ensureAvgResults()

        if self.isChainSolver():
            # Chain mode: no phases, so the aggregate space is the state space.
            return self._result.space.copy()

        # Use pre-computed aggregated state space (nstates, nstations * nclasses)
        if hasattr(self._result, 'space_aggr') and self._result.space_aggr is not None:
            return self._result.space_aggr.copy()

        # Fallback: return raw space
        space = self._result.space
        if space is None or len(space) == 0:
            return np.array([])
        return space.copy()

    def getCdfSysRespT(self) -> List[Dict]:
        """
        The SYSTEM response-time distribution: one law per CHAIN.

        THE QUANTITY IS THE CYCLE TIME. The split is the tagged job's ARRIVAL AT
        ITS OWN REFERENCE STATION, so a passage runs from one such arrival to the
        next: the job's whole trip round the network, not its stay at one
        station. A single MAP suffices here where getCdfRespT needs two, because
        the arrival that starts the passage and the one that ends it are the same
        event.

        THIS REPLACED AN EXPONENTIAL FIT to the mean system response time, which
        returned one entry per CLASS. The law is per CHAIN, as it is in MATLAB
        (RD = cell(1, sn.nchains)) and C++, so the 'chain' key replaces 'class'.

        Two constants differ from the per-station getter on purpose, matching the
        reference: the grid is 10000 intervals rather than 100000, and the
        truncation is at 1 - 1e-8 rather than 1e-3, because a cycle time is
        longer and its tail matters more.

        Reference: matlab/src/solvers/CTMC/@SolverCTMC/getCdfSysRespT.m.

        Returns:
            List of dicts with 'chain', 't', 'p' keys
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # The two compute the SAME quantity, one law per chain, so this
            # delegates rather than refusing. The older refusal reason -- that
            # the native getter fitted an exponential per class -- stopped being
            # true when this getter was rewritten onto the tagged chain.
            from ..cpp_dispatch import cdf_sys_respt_via_cpp
            return cdf_sys_respt_via_cpp(self)

        from scipy.linalg import expm as _expm
        from ...api.io.model_adapter import tag_chain
        from ...api.mam.map_analysis import map_normalize, map_pie
        from ...constants import EventType
        from ...lang.sync import refresh_sync

        sn = self._sn if self._sn is not None else self.model.get_struct()
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        if np.any(np.isinf(njobs)):
            raise RuntimeError(
                "getCdfSysRespT is presently supported only for closed models.")

        classes = self.model.get_classes()
        RD: List[Dict] = []

        class _Chain(object):
            def __init__(self, cls):
                self.classes = cls

        for ch in range(int(sn.nchains)):
            inchain = [int(x) for x in np.asarray(sn.inchain[ch]).ravel()]
            tagged_src_idx = None
            for r in inchain:
                if njobs[r] > 0:
                    tagged_src_idx = r
                    break
            if tagged_src_idx is None:
                continue

            tagged = tag_chain(self.model,
                               _Chain([classes[r] for r in inchain]),
                               classes[tagged_src_idx])
            tsolver = SolverCTMC(tagged.model, self.options)
            Q, filt = tsolver.getGenerator()
            Q = np.asarray(Q, dtype=float)
            tsn = tagged.model.get_struct()
            sync = refresh_sync(tsn)
            if sync is None or filt is None or len(filt) == 0:
                raise RuntimeError(
                    "getCdfSysRespT needs the event filtration of the tagged chain, which this "
                    "model did not produce; the system response-time law cannot be computed "
                    "without it")

            tagged_cls = tagged.tagged_job._index \
                if hasattr(tagged.tagged_job, '_index') else len(tagged.model.get_classes()) - 1
            # sn.refstat is a STATION index; the events carry NODE indices, so
            # the two must be mapped rather than compared directly.
            ref_station = int(np.asarray(tsn.refstat).ravel()[tagged_cls])
            station_to_node = np.asarray(tsn.stationToNode).ravel()
            ref_node = int(station_to_node[ref_station]) \
                if 0 <= ref_station < station_to_node.size else ref_station

            D1 = np.zeros_like(Q)
            for v, ev in enumerate(sync):
                if v >= len(filt) or filt[v] is None:
                    continue
                pas = getattr(ev, 'passive', None)
                if (pas is not None and pas.event == EventType.ARV
                        and pas.job_class == tagged_cls
                        and pas.node == ref_node):
                    D1 = D1 + np.asarray(filt[v], dtype=float)
            if not np.any(D1):
                continue

            D0, D1n = map_normalize(Q - D1, D1)
            pie = np.asarray(map_pie(D0, D1n), dtype=float).ravel()

            nz = np.abs(Q[Q != 0])
            nz = nz[nz > 1e-8]
            if nz.size == 0:
                continue
            intervals = 10000
            T = abs(100.0 / nz.min())
            dT = T / intervals
            E = _expm(D0 * dT)
            ones = np.ones(D0.shape[0])
            v = pie.copy()
            tvals = []
            Fvals = []
            for k in range(intervals + 1):
                if k > 0:
                    v = v.dot(E)
                Fk = min(1.0, max(0.0, 1.0 - float(v.dot(ones))))
                tvals.append(k * dT)
                Fvals.append(Fk)
                if Fk > 1.0 - 1e-8:
                    break

            RD.append({
                'chain': ch + 1,
                't': np.array(tvals),
                'p': np.array(Fvals),
            })

        return RD

    def getReward(self, reward_vector: Optional[np.ndarray] = None) -> float:
        """Compute reward function over steady-state distribution.

        Args:
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Expected reward
        """
        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        if pi is None:
            return 0.0

        if reward_vector is None:
            # Default: expected queue length
            space = self._result.space
            if space is None:
                return 0.0
            reward_vector = np.sum(space, axis=1)

        return np.dot(pi, reward_vector)

    def getAvgReward(self) -> Tuple[np.ndarray, List[str]]:
        """Get steady-state expected reward values.

        Computes the steady-state expected reward for reward functions
        previously defined using model.setReward().

        Returns:
            Tuple of (R, names) where:
            - R: numpy array of expected reward values
            - names: list of reward function names

        Example:
            >>> model.setReward('QueueLength', lambda state: state.at(queue, oclass))
            >>> solver = CTMC(model)
            >>> R, names = solver.getAvgReward()
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            # reward callables can't cross the JSON round-trip to jline.jar; evaluated on a native CTMC steady-state solve instead (same distribution as JAR's).
            saved_lang = self.options.lang
            saved_result = self._result
            try:
                self.options.lang = 'python'
                self._result = None
                self.runAnalyzer()
                return self._compute_avg_reward()
            finally:
                self.options.lang = saved_lang
                self._result = saved_result
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # NOT the java branch's native fallback. A reward built from a
            # `Reward.*` template IS serialized by linemodel_save in its
            # declarative {name, type, node, class} form, so the C++ receives the
            # declaration and evaluates it against its own stationary law; only a
            # bare lambda cannot cross, and the WRITER refuses that by name
            # rather than this branch quietly solving natively.
            from ..cpp_dispatch import avg_reward_via_cpp
            return avg_reward_via_cpp(self)
        if self._result is None:
            self._ensureAvgResults()
        return self._compute_avg_reward()

    def _compute_avg_reward(self) -> Tuple[np.ndarray, List[str]]:
        """Evaluate model reward functions against the current steady-state
        distribution (self._result.pi over the aggregated state space)."""
        # Get rewards from the model
        if hasattr(self.model, 'get_rewards'):
            rewards_dict = self.model.get_rewards()
        elif hasattr(self.model, '_rewards'):
            rewards_dict = self.model._rewards
        else:
            rewards_dict = {}

        if not rewards_dict:
            return np.array([]), []

        # Get aggregated state space and steady-state probabilities
        # Use space_aggr (aggregated per-station counts) not raw space (phase-level detail)
        pi = self._result.pi
        space = self._result.space_aggr if hasattr(self._result, 'space_aggr') and self._result.space_aggr is not None else self._result.space

        if pi is None or len(pi) == 0 or space is None or len(space) == 0:
            return np.array([0.0] * len(rewards_dict)), list(rewards_dict.keys())

        # Build mappings for RewardState
        from ...lang.reward_state import RewardState
        sn = self._sn

        # Build node-to-station mapping
        nodes_to_station = {}
        if hasattr(self.model, 'get_nodes'):
            for node in self.model.get_nodes():
                if hasattr(node, 'get_index'):
                    node_idx = node.get_index()
                elif hasattr(node, 'index'):
                    node_idx = node.index
                else:
                    continue

                if hasattr(sn, 'nodeToStation'):
                    node_idx0 = node_idx - 1  # 0-indexed
                    if node_idx0 < len(sn.nodeToStation):
                        station_idx = int(sn.nodeToStation[node_idx0])
                        if station_idx >= 0:
                            nodes_to_station[node_idx] = station_idx + 1  # 1-based

        # Build class-to-index mapping
        classes_to_idx = {}
        if hasattr(self.model, 'get_classes'):
            for i, jobclass in enumerate(self.model.get_classes()):
                if hasattr(jobclass, 'get_index'):
                    class_idx = jobclass.get_index()
                elif hasattr(jobclass, 'index'):
                    class_idx = jobclass.index
                else:
                    class_idx = i + 1
                classes_to_idx[class_idx] = i + 1

        # Compute expected reward for each reward function
        R = []
        names = list(rewards_dict.keys())

        for name, reward_fn in rewards_dict.items():
            expected_value = 0.0

            for state_idx, state_vec in enumerate(space):
                prob = pi[state_idx]
                if prob <= 0:
                    continue

                # Create RewardState for this state vector
                reward_state = RewardState(state_vec, sn, nodes_to_station, classes_to_idx)

                # Evaluate reward function
                try:
                    # Check if reward_fn takes sn argument (for Reward templates)
                    import inspect
                    sig = inspect.signature(reward_fn)
                    if len(sig.parameters) >= 2:
                        reward_value = reward_fn(reward_state, sn)
                    else:
                        reward_value = reward_fn(reward_state)
                    expected_value += prob * reward_value
                except Exception as e:
                    # Skip if reward function fails for this state
                    pass

            R.append(expected_value)

        return np.array(R), names

    get_avg_reward = getAvgReward

    def getTranCdfRespT(self, t_max: float = 10.0, n_points: int = 100) -> List[Dict]:
        """Not supported, as in the reference, whose base class raises.

        Returning the steady-state law under the transient getter's name would
        be indistinguishable, to the caller, from a transient analysis.
        """
        raise NotImplementedError("getTranCdfRespT is not supported by SolverCTMC")

    # =========================================================================
    # Transient Probability Methods
    # =========================================================================

    def getTranProb(self, node: int, t: float = 1.0) -> np.ndarray:
        """Get transient state probabilities at a node.

        Computes π(t) = π(0) * exp(Q*t) using matrix exponential.

        Args:
            node: Node/station index (0-based)
            t: Time point for transient analysis

        Returns:
            Transient probability vector at time t
        """
        self._assert_phasetype_states('getTranProb')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # NOT A MISSING ARM, A DIFFERENT SHAPE. `-a tranprob` returns the
            # FULL occupancy pi(t) beside the labelled state space, while this
            # getter returns a marginal indexed by one column of the flat space.
            # The bucketing IS this getter's definition, so it is applied here to
            # the C++'s law over the C++'s own enumeration -- rather than the
            # native path being run under the C++ engine's name.
            from ..cpp_dispatch import tran_prob_via_cpp
            d = tran_prob_via_cpp(self, t)
            pi_t = np.asarray(d['pit'][-1, :]).reshape(-1)
            space = d['labels']
            if space.size == 0 or node >= space.shape[1]:
                return pi_t
            col = np.asarray(space[:, node], dtype=int)
            marginal = np.zeros(int(col.max()) + 1)
            for s, prob in enumerate(pi_t):
                marginal[col[s]] += prob
            return marginal
        if self._result is None:
            self._ensureAvgResults()

        from scipy.linalg import expm

        Q = self._result.infgen
        if Q is None:
            # Fall back to steady-state
            return self.getProb(node)

        # pi(0) IS THE MODEL'S INITIAL STATE, located in the enumerated space;
        # see _network_init_distribution for what seeding e_0 instead cost.
        pi_0 = self._network_init_distribution()

        # Compute transient probability: π(t) = π(0) * exp(Q*t)
        pi_t = pi_0 @ expm(Q * t)

        # Extract marginal for node
        space = self._result.space
        if space is None or node >= space.shape[1]:
            return pi_t

        # Compute marginal probability for node
        max_n = int(np.max(space[:, node])) + 1
        marginal = np.zeros(max_n)

        for s, prob in enumerate(pi_t):
            n = int(space[s, node])
            if n < max_n:
                marginal[n] += prob

        return marginal

    def getTranProbAggr(self, node: int, t: float = 1.0) -> np.ndarray:
        """Get transient aggregated state probabilities at a node.

        Args:
            node: Node/station index (0-based)
            t: Time point for transient analysis

        Returns:
            Transient aggregated probability vector at time t
        """
        self._assert_phasetype_states('getTranProbAggr')
        return self.getTranProb(node, t)

    def getTranProbSys(self, t: float = 1.0) -> np.ndarray:
        """Get transient system state probabilities.

        Computes full system state probability at time t.

        In chain mode the distribution starts from options.init_sol, or from the
        uniform distribution when none is given; a DTMC advances one step per
        unit of time, so t must then be a non-negative integer.

        Args:
            t: Time point for transient analysis

        Returns:
            Transient system probability vector at time t
        """
        if self.isChainSolver():
            self._ensureAvgResults()
            pi0 = self._chain_init_distribution()
            if self.isDiscreteChain():
                if t < 0 or abs(t - round(t)) > 1e-12:
                    raise RuntimeError(
                        "A DTMC advances one step per unit of time, so getTranProbSys "
                        "requires a non-negative integer number of steps.")
                from ...api.mc import dtmc_transient
                return dtmc_transient(self.getTransMat(), pi0, int(round(t)))[-1]
            from ...api.mc import ctmc_transient
            return np.asarray(ctmc_transient(self._result.infgen, pi0, float(t))).flatten()
        self._assert_phasetype_states('getTranProbSys')
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # The occupancy vector at t, taken as the LAST row of the trajectory
            # the C++ integrates over [0, t]: pi(t) is what this getter returns,
            # and the horizon it was reached over is what the C++ requires to be
            # stated.
            from ..cpp_dispatch import tran_prob_via_cpp
            pit = tran_prob_via_cpp(self, t)['pit']
            return np.asarray(pit[-1, :]).reshape(-1)
        if self._result is None:
            self._ensureAvgResults()

        from scipy.linalg import expm

        Q = self._result.infgen
        if Q is None:
            # Fall back to steady-state
            return self.getSteadyState()

        # pi(0) IS THE MODEL'S INITIAL STATE, located in the enumerated space;
        # see _network_init_distribution for what seeding e_0 instead cost.
        pi_0 = self._network_init_distribution()

        # Compute transient probability: π(t) = π(0) * exp(Q*t)
        pi_t = pi_0 @ expm(Q * t)

        return pi_t

    def getTranProbSysAggr(self, t: float = 1.0) -> np.ndarray:
        """Get transient aggregated system state probabilities.

        Args:
            t: Time point for transient analysis

        Returns:
            Transient system probability vector at time t
        """
        self._assert_phasetype_states('getTranProbSysAggr')
        return self.getTranProbSys(t)

    # =========================================================================
    # Symbolic Generator Methods
    # =========================================================================

    def getSymbolicGenerator(self, invert_symbol: bool = False):
        """Get symbolic generator matrix with per-event symbolic variables.

        Each event filtration matrix is normalized and multiplied by a symbolic
        variable (x1, x2, ...), matching MATLAB's getSymbolicGenerator.m.

        Args:
            invert_symbol: If True, divide by symbol instead of multiplying

        Returns:
            Tuple of (infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace):
                - infGen: Symbolic infinitesimal generator (sympy.Matrix)
                - eventFilt: List of per-event symbolic filtration matrices
                  (None for events with no positive rates, matching
                  MATLAB's empty cells)
                - syncInfo: Sync data structure from the model
                - stateSpace: State space matrix
                - nodeStateSpace: Per-node state space
        """
        try:
            import sympy
        except ImportError:
            raise ImportError(
                "sympy is required for symbolic generator. "
                "Install it with 'pip install sympy'."
            )

        _, F = self.getGenerator()
        stateSpace, nodeStateSpace = self.getStateSpace()

        if not F:
            return None, [], None, stateSpace, nodeStateSpace

        n = F[0].shape[0]
        n_events = len(F)

        infGen = sympy.zeros(n, n)
        eventFilt = [None] * n_events

        for e in range(n_events):
            Fe = np.asarray(F[e], dtype=np.float64)
            pos = Fe[Fe > 0]
            if len(pos) > 0:
                minF = pos.min()
                Fe = Fe / minF
                xe = sympy.Symbol(f'x{e + 1}', real=True)
                if invert_symbol:
                    Fe_sym = sympy.Matrix(Fe.tolist()) / xe
                else:
                    Fe_sym = sympy.Matrix(Fe.tolist()) * xe
                eventFilt[e] = Fe_sym
                infGen = infGen + Fe_sym

        from ...api.mc.ctmc import ctmc_makeinfgen
        infGen = ctmc_makeinfgen(infGen)

        sn = self._sn if self._sn is not None else self.model.get_struct()
        syncInfo = sn.sync if hasattr(sn, 'sync') else None

        return infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace

    get_symbolic_generator = getSymbolicGenerator

    # =========================================================================
    # Parametric sensitivity
    # =========================================================================

    def symbolicBackend(self):
        """Value of options.config['symbolic'], or 'auto' when unset.

        'auto' keeps the native engine (sympy), exactly as MATLAB's 'auto'
        keeps the Symbolic Math Toolbox when it is licensed.
        """
        config = getattr(self.options, 'config', None)
        if isinstance(config, dict):
            return config.get('symbolic', 'auto')
        if config is not None and hasattr(config, 'symbolic'):
            return config.symbolic
        return 'auto'

    symbolic_backend = symbolicBackend

    def _symbolicTimeout(self):
        """Value of options.config['symbolic_timeout'], default 300 s."""
        config = getattr(self.options, 'config', None)
        if isinstance(config, dict):
            return config.get('symbolic_timeout', 300)
        if config is not None and hasattr(config, 'symbolic_timeout'):
            return config.symbolic_timeout
        return 300

    def _resolveStateSet(self, S, n, name):
        """Resolve a state set given as 1-based row indices or as state rows.

        An unrecognised row is an error rather than a silent drop, since a
        passage into a state that is not in the space is not a slow passage but
        an undefined one. Mirrors the local_one helper of MATLAB
        @SolverCTMC/getCdfFirstPassT.m.
        """
        from ...api.pfqn.utils import matchrow

        if S is None:
            return np.array([], dtype=int)
        arr = np.atleast_1d(np.asarray(S))
        if arr.size == 0:
            return np.array([], dtype=int)
        if arr.ndim == 1 and np.all(arr == np.round(arr)) \
                and np.all(arr >= 1) and np.all(arr <= n):
            # 1-based row indices, as in MATLAB; stored 0-based here
            return np.unique(arr.astype(int)) - 1
        arr = np.atleast_2d(arr)
        space, _ = self.getStateSpace()
        idx = np.zeros(arr.shape[0], dtype=int)
        for i in range(arr.shape[0]):
            r = matchrow(np.asarray(space), np.asarray(arr[i, :]).ravel())
            if r <= 0:
                raise ValueError('A state given in set %s is not in the state space.'
                                 % name)
            idx[i] = r - 1
        return np.unique(idx)

    def _passageInitial(self, Aidx, n):
        """Uniform initial law on A, or None to start from the conditional
        stationary law on the complement of B."""
        if Aidx is None or len(Aidx) == 0:
            return None
        pi0 = np.zeros(n)
        pi0[Aidx] = 1.0 / len(Aidx)
        return pi0

    def getCdfFirstPassT(self, A, B):
        """Distribution of the FIRST PASSAGE TIME from state set A into set B.

        Mirrors MATLAB ``@SolverCTMC/getCdfFirstPassT.m``. RD is an (n, 2) array
        whose first column is F(t) and whose second is t, the column order every
        other CDF getter in LINE uses.

        A and B name states either as 1-based ROW INDICES into the state space
        returned by getStateSpace, or as matrices of state rows, which are
        resolved against that space. An empty A starts from the conditional
        stationary law on the complement of B.

        THIS IS NOT getCdfRespT. That getter times a tagged job between an
        arrival at a station and its departure, through the event filtration;
        this one times the chain between two sets of states the caller names,
        and answers questions the filtration cannot express -- the writer cycle
        time of a readers-writers model, the time to fill a buffer, the time to
        leave a degraded region.

        Args:
            A: source state set, or empty for the conditional stationary law
            B: target state set, which may not be empty

        Returns:
            (RD, out) with RD the (n, 2) [F(t), t] array and out the dict
            returned by ctmc_passage_time, extended with tset, density, source,
            target and runtime.

        References:
            P. G. Harrison and W. J. Knottenbelt, "Passage Time Distributions in
            Large Markov Chains", 2002.
        """
        import time as _time
        from ...api.mc.passage import ctmc_passage_time
        from ...constants import GlobalConstants

        t0 = _time.time()
        Q = np.asarray(self.getInfGen(), dtype=float)
        n = Q.shape[0]

        Bidx = self._resolveStateSet(B, n, 'B')
        if len(Bidx) == 0:
            raise ValueError('The target state set B is empty: a first passage '
                             'time into no state is undefined.')
        Aidx = self._resolveStateSet(A, n, 'A')

        config = getattr(self.options, 'config', None)
        method = None
        if isinstance(config, dict):
            method = config.get('passage_method', None)
        elif config is not None and hasattr(config, 'passage_method'):
            method = config.passage_method
        if not method:
            method = 'expm'

        pi0 = self._passageInitial(Aidx, n)

        # The horizon is chosen the way the response-time getter chooses it:
        # 100 events at the slowest rate in the chain.
        nonzero = np.abs(Q[Q != 0])
        nonzero = nonzero[nonzero > GlobalConstants.FineTol]
        thor = abs(100.0 / np.min(nonzero))
        tset = np.linspace(0.0, thor, 1000)

        F, f, out = ctmc_passage_time(Q, pi0, Bidx, tset, method=method)
        RD = np.column_stack([np.asarray(F).ravel(), tset])
        out['tset'] = tset
        out['density'] = f
        out['source'] = Aidx
        out['target'] = Bidx
        out['runtime'] = _time.time() - t0
        return RD, out

    get_cdf_first_pass_t = getCdfFirstPassT

    def getFirstPassTMoments(self, A, B, nmax: int = 3):
        """Moments of order 1..nmax of the first passage time from A into B.

        Mirrors MATLAB ``@SolverCTMC/getFirstPassTMoments.m``.

        NO TRANSFORM INVERSION AND NO TIME GRID ARE INVOLVED. The moments come
        from Eq. 3 of Harrison and Knottenbelt (2002) -- one linear solve per
        order -- so they are exact and are not limited by the horizon a CDF
        would have to be truncated at. This is the cheapest way to get the
        variance or the skewness of a passage time in LINE.

        Args:
            A: source state set, named as in getCdfFirstPassT
            B: target state set, which may not be empty
            nmax: highest moment order, default 3

        Returns:
            (m, mall) with m the (nmax,) moment vector for a passage started
            uniformly in A, and mall (nstates, nmax) one row per starting
            state, zero on B and inf where B cannot be reached.
        """
        from ...api.mc.passage import ctmc_passage_moments

        if nmax is None:
            nmax = 3
        nmax = int(nmax)
        Q = np.asarray(self.getInfGen(), dtype=float)
        n = Q.shape[0]

        Bidx = self._resolveStateSet(B, n, 'B')
        if len(Bidx) == 0:
            raise ValueError('The target state set B is empty: a first passage '
                             'time into no state is undefined.')
        Aidx = self._resolveStateSet(A, n, 'A')
        pi0 = self._passageInitial(Aidx, n)

        mall, m = ctmc_passage_moments(Q, pi0, Bidx, nmax)
        return m, mall

    get_first_pass_t_moments = getFirstPassTMoments

    def getSensitivity(self, param, reward=None, method: str = 'fd'):
        """Parametric sensitivity of a steady-state reward to a scalar model
        parameter, following Trivedi and Bobbio (2017), Sec. 9.7.

        Mirrors MATLAB ``@SolverCTMC/getSensitivity.m``.

        Args:
            param: dict describing the parameter theta and how to set it:
                ``name`` identifier used in reports;
                ``value`` nominal value theta;
                ``set`` callable (model, value) -> None applying theta;
                ``step`` optional finite-difference step, default value*1e-6.
            reward: reward rate vector over the states, or a callable mapping
                the state space to one. If omitted, dpi is returned and S is
                None.
            method: 'fd' (default) or 'symbolic'.

                'fd' obtains the generator derivative dQ/dtheta by central
                differences on the rate with the state space held fixed. This
                is exact to O(step^2) and requires no symbolic differentiation
                of the rate assembly; the state space is unaffected because it
                depends on the topology and the cutoff, not on rate values.
                The steady-state sensitivity then follows from one linear
                solve, see ctmc_sens.

                'symbolic' solves the stationary distribution as a rational
                function of the event rate symbols x1..xE and differentiates
                it exactly with respect to each of them, then combines by the
                chain rule::

                    d(pi)/d(theta) = sum_e d(pi)/d(x_e) * d(x_e)/d(theta).

                Only the rate map x_e(theta) is still differenced, and that
                map is affine in theta in the common cases (a rate set to
                theta, or scaled by it), where the central difference
                reproduces it exactly. The whole O(step^2) error of 'fd' comes
                from differencing through the solve, which this avoids
                entirely. It refuses rather than approximates when perturbing
                theta reshapes an event's filtration instead of scaling it.

        Returns:
            (S, SS, dpi, pi) with S the unscaled sensitivity d(E[r])/dtheta,
            Eq. (9.79); SS the scaled sensitivity (theta/E[r]) d(E[r])/dtheta,
            Eq. (9.80); dpi the sensitivity of the steady-state distribution;
            pi the steady-state distribution.

        Note:
            This returns d(E[r])/dtheta with dr/dtheta = 0, i.e. it assumes
            the reward rates do not themselves depend on theta. Rewards that
            depend on theta need the second term of Eq. (9.83) and are not
            handled here.
        """
        from ...api.mc.ctmc import ctmc_solve, ctmc_sens
        from ...constants import GlobalConstants

        if not isinstance(param, dict) or 'set' not in param or 'value' not in param:
            raise ValueError("param must be a dict with keys 'value' and 'set'")
        if method is None or method == '':
            method = 'fd'
        if method.lower() not in ('fd', 'symbolic'):
            raise ValueError("unknown method '%s'; expected 'fd' or 'symbolic'" % method)

        theta = float(param['value'])
        step = param.get('step', None)
        h = float(step) if step else max(abs(theta), 1.0) * 1e-6

        # Nominal generator and state space
        Q, _ = self.getGenerator()
        Q = np.asarray(Q.todense() if hasattr(Q, 'todense') else Q, dtype=np.float64)
        space = self.getStateSpace()[0]
        n = Q.shape[0]

        if method.lower() == 'symbolic':
            dpi, pi = self._symbolicSensitivity(param, theta, h, n)
        else:
            # Central differences on theta with the state space fixed
            Qp, _ = self._perturbedGenerator(param, theta + h)
            Qm, _ = self._perturbedGenerator(param, theta - h)
            if Qp.shape[0] != n or Qm.shape[0] != n:
                raise ValueError(
                    'Perturbing the parameter changed the state space size, so the '
                    'generators cannot be differenced. This happens when the parameter '
                    'switches a transition on or off (e.g. a zero rate or an immediate '
                    'transition).')
            dQ = (Qp - Qm) / (2 * h)

            pi = np.asarray(ctmc_solve(Q), dtype=np.float64).flatten()
            dpi = ctmc_sens(Q, dQ, pi)

        if reward is None:
            return None, None, dpi, pi

        r = reward(space) if callable(reward) else reward
        r = np.asarray(r, dtype=np.float64).flatten()
        if r.size != n:
            raise ValueError('reward must have one entry per state')

        # Eq. (9.83) with dr/dtheta = 0
        S = float(dpi @ r)
        Er = float(pi @ r)
        SS = (theta / Er) * S if abs(Er) > GlobalConstants.Zero else float('nan')
        return S, SS, dpi, pi

    get_sensitivity = getSensitivity

    def _symbolicSensitivity(self, param, theta, h, n):
        """Exact d(pi)/d(x_e), combined with a differenced rate map
        d(x_e)/d(theta) by the chain rule.

        The split matters: the stationary distribution is a rational function
        of the rates of high degree, and differencing through it is where the
        O(h^2) error of the 'fd' method comes from. The rate map, by contrast,
        is affine in theta whenever theta is a rate or scales one, and a
        central difference is exact on an affine map. What is left is exact in
        those cases and no worse otherwise.
        """
        # symbolic gen: event filtration normalized by own min positive rate; x_e nominal = that rate; see _kb/06-solver-catalog.md CTMC Symbolic analysis.
        infGen = self.getSymbolicGenerator()[0]
        _, F = self.getGenerator()
        nEvents = len(F)
        rate0, shape0 = _eventRates(F)

        # rate-map differencing: wide step (not fd), affine map exact, tiny step cancellation-dominated; see _kb/06-solver-catalog.md CTMC Symbolic analysis.
        hRate = max(abs(theta), 1.0) * 1e-3
        Qp, Fp = self._perturbedGenerator(param, theta + hRate)
        Qm, Fm = self._perturbedGenerator(param, theta - hRate)
        if Qp.shape[0] != n or Qm.shape[0] != n:
            raise ValueError(
                'Perturbing the parameter changed the state space size, so the '
                'generators cannot be differenced. This happens when the parameter '
                'switches a transition on or off (e.g. a zero rate or an immediate '
                'transition).')
        if len(Fp) != nEvents or len(Fm) != nEvents:
            raise ValueError('Perturbing the parameter changed the number of events.')
        ratep, shapep = _eventRates(Fp)
        ratem, shapem = _eventRates(Fm)
        for e in range(nEvents):
            if shape0[e] is None:
                continue
            if (shapep[e] is None or shapem[e] is None
                    or shapep[e].shape != shape0[e].shape
                    or np.max(np.abs(shapep[e] - shape0[e])) > 1e-8
                    or np.max(np.abs(shapem[e] - shape0[e])) > 1e-8):
                raise ValueError(
                    'Perturbing the parameter reshapes the filtration of event %d '
                    'rather than scaling it, so the generator is not linear in a '
                    'single rate per event and the symbolic chain rule does not '
                    "apply. Use the 'fd' method for this parameter." % (e + 1))
        # affine-at-this-scale check via the midpoint identity r(+)+r(-)=2r(0); falls back to the caller's small step when it fails.
        curvature = np.abs(ratep + ratem - 2 * rate0)
        scale = max(1.0, float(np.max(np.abs(rate0))) if rate0.size else 1.0)
        if curvature.size and np.max(curvature) > 1e-9 * scale:
            ratep = _eventRates(self._perturbedGenerator(param, theta + h)[1])[0]
            ratem = _eventRates(self._perturbedGenerator(param, theta - h)[1])[0]
            drate = (ratep - ratem) / (2 * h)
        else:
            drate = (ratep - ratem) / (2 * hRate)

        symbols = [('x%d' % (e + 1)) if shape0[e] is not None else None
                   for e in range(nEvents)]
        active = [e for e in range(nEvents) if symbols[e] is not None]
        assignment = dict((symbols[e], float(rate0[e])) for e in active)

        # symbolic stationary distribution and per-symbol exact derivative, engine selected by options.config['symbolic'] (sympy or line-sage-rest).
        from ...api.sym import resolve as _resolve_sym, require as _require_sym
        backend = str(self.symbolicBackend()).strip()
        useService = (backend.lower() == 'sage' or backend.lower().startswith('http'))
        if useService:
            engine = _resolve_sym(backend)
            if engine is None:
                engine = _require_sym(backend)
            engine.timeout_s = self._symbolicTimeout()
            Qtext = [[str(infGen[i, j]) for j in range(n)] for i in range(n)]
            piExpr = engine.solve_ctmc(Qtext, [symbols[e] for e in active])['pi']
            pi = np.asarray(engine.eval(piExpr, assignment)[0], dtype=np.float64)
            dpi = np.zeros(n)
            for e in active:
                if drate[e] == 0:
                    # This event does not depend on theta, so its term is zero
                    # and the derivative is not worth a round trip.
                    continue
                dExpr = engine.diff(piExpr, symbols[e], 1)
                dvals = np.asarray(engine.eval(dExpr, assignment)[0], dtype=np.float64)
                dpi = dpi + drate[e] * dvals
            return dpi, pi

        import sympy
        from ...api.mc.ctmc import ctmc_solve
        piExpr = ctmc_solve(infGen)
        piExpr = [sympy.together(piExpr[k]) for k in range(n)]
        subs = dict((sympy.Symbol(k, real=True), sympy.Float(v, 17))
                    for k, v in assignment.items())
        pi = np.array([float(expr.subs(subs)) for expr in piExpr], dtype=np.float64)
        dpi = np.zeros(n)
        for e in active:
            if drate[e] == 0:
                continue
            xe = sympy.Symbol(symbols[e], real=True)
            dvals = np.array([float(sympy.diff(expr, xe).subs(subs)) for expr in piExpr],
                             dtype=np.float64)
            dpi = dpi + drate[e] * dvals
        return dpi, pi

    def _perturbedGenerator(self, param, value):
        """Rebuild the generator with theta set to VALUE, on a copy of the
        model so the caller's model is left untouched.

        The hard refresh is required, not defensive: set_service and
        set_arrival deliberately leave the cached struct in place, so a copy
        that inherited a built struct would report the old rate and the
        difference quotient would silently come out as zero.
        """
        import dataclasses

        modelCopy = self.model.copy()
        param['set'](modelCopy, value)
        modelCopy.refresh_struct()
        solverCopy = SolverCTMC(modelCopy)
        # the whole options record (not just constructor-forwarded fields) is used, since cutoff/config/gen_method all change which generator gets built.
        solverCopy.options = dataclasses.replace(self.options)
        solverCopy.method = self.options.method
        Q, F = solverCopy.getGenerator()
        Q = np.asarray(Q.todense() if hasattr(Q, 'todense') else Q, dtype=np.float64)
        return Q, F

    def getMarkedCTMC(self) -> Dict[str, Any]:
        """Get a marked CTMC object representation.

        Returns a dictionary containing the CTMC with marked transitions
        for reward and passage time analysis.

        Returns:
            Dictionary with 'Q' (generator), 'space' (state space),
            'pi' (steady-state), and 'marks' (transition markings)
        """
        if self._result is None:
            self._ensureAvgResults()

        Q = self._result.infgen
        space = self._result.space
        pi = self._result.pi

        # Create transition markings (identify each transition type)
        n_states = Q.shape[0] if Q is not None else 0
        marks = {}

        if Q is not None:
            mark_id = 0
            for i in range(n_states):
                for j in range(n_states):
                    if i != j and Q[i, j] != 0:
                        marks[(i, j)] = {
                            'id': mark_id,
                            'rate': Q[i, j],
                            'from_state': i,
                            'to_state': j,
                        }
                        mark_id += 1

        return {
            'Q': Q,
            'space': space,
            'pi': pi,
            'marks': marks,
            'n_states': n_states,
            'n_transitions': len(marks),
        }

    # =========================================================================
    # Reward Analysis Methods
    # =========================================================================

    def runRewardAnalyzer(self, reward_vector: Optional[np.ndarray] = None) -> Dict[str, Any]:
        """Run reward analysis on the CTMC.

        Computes expected rewards in steady-state and optionally transient.

        Args:
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Dictionary with 'steady_state_reward', 'reward_per_state', etc.
        """
        if self._result is None:
            self._ensureAvgResults()

        pi = self._result.pi
        space = self._result.space

        if reward_vector is None:
            # Default: use total queue length as reward
            if space is not None:
                reward_vector = np.sum(space, axis=1)
            else:
                return {'steady_state_reward': 0.0, 'error': 'No state space available'}

        if pi is None:
            return {'steady_state_reward': 0.0, 'error': 'No steady-state distribution'}

        # Compute steady-state reward
        steady_state_reward = np.dot(pi, reward_vector)

        # Compute per-state rewards
        reward_per_state = pi * reward_vector

        return {
            'steady_state_reward': steady_state_reward,
            'reward_per_state': reward_per_state,
            'reward_vector': reward_vector,
            'pi': pi,
        }

    def getTranReward(self, t: float = 1.0,
                      reward_vector: Optional[np.ndarray] = None) -> float:
        """Get transient reward at time t.

        Computes expected reward at time t using matrix exponential.

        Args:
            t: Time point for transient analysis
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Expected reward at time t
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            # reward vector and transient distribution share a state ordering, neither crossing the JSON round-trip; run on native CTMC solve, like getAvgReward.
            saved_lang = self.options.lang
            saved_result = self._result
            try:
                self.options.lang = 'python'
                self._result = None
                self.runAnalyzer()
                return self._compute_tran_reward(t, reward_vector)
            finally:
                self.options.lang = saved_lang
                self._result = saved_result
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import tran_reward_at_via_cpp
            return tran_reward_at_via_cpp(self, t, reward_vector)
        if self._result is None:
            self._ensureAvgResults()
        return self._compute_tran_reward(t, reward_vector)

    def _compute_tran_reward(self, t: float,
                             reward_vector: Optional[np.ndarray]) -> float:
        """Expected reward at time t: pi(t) . r, with pi(t) the transient system
        distribution and r a state-indexed reward vector (default: total jobs)."""
        space = self._result.space

        if reward_vector is None:
            if space is not None:
                reward_vector = np.sum(space, axis=1)
            else:
                return 0.0

        # Get transient probabilities
        pi_t = self.getTranProbSys(t)

        # Compute transient reward
        return np.dot(pi_t, reward_vector)

    def get_tran_reward(self, name: Optional[str] = None):
        """Transient expected reward E[r(X(t))] over time for each reward.

        This is the transient counterpart of get_avg_reward: instead of the
        equilibrium value it returns the time-indexed trajectory E[r(X(t))],
        where X(t) is the system state at time t and the expectation is taken
        over the CTMC transient distribution starting from the initial state.
        The named reward functions defined via model.setReward are evaluated
        on the aggregated state space. A finite timespan is required, e.g.
        CTMC(model, timespan=[0, T]).

        Args:
            name: optional reward name; if given, only that reward is returned.

        Returns:
            Tuple (Rt, t, names) where Rt is a list of dicts with keys
            't', 'metric', 'name' (one per reward), or a single dict when name
            is specified; t is the array of time points; names is the list of
            reward names (or a single name when name is specified).
        """
        timespan = getattr(self.options, 'timespan', None)
        if timespan is None or not np.isfinite(timespan[1]):
            raise ValueError('get_tran_reward requires a finite timespan, '
                             'e.g. CTMC(model, timespan=[0, T]).')

        if hasattr(self.model, 'get_rewards'):
            rewards_dict = self.model.get_rewards()
        elif hasattr(self.model, '_rewards'):
            rewards_dict = self.model._rewards
        else:
            rewards_dict = {}
        if not rewards_dict:
            raise ValueError('No rewards defined. Use model.setReward(name, fn) '
                             'before calling get_tran_reward.')

        if getattr(self.options, 'lang', 'python') == 'java':
            # reward callables cannot cross the JSON round-trip; evaluated on a native CTMC transient solve, mirroring getAvgReward/getTranReward.
            saved_lang = self.options.lang
            saved_result = self._result
            try:
                self.options.lang = 'python'
                self._result = None
                self.runAnalyzer()
                return self._compute_tran_reward_named(name, rewards_dict)
            finally:
                self.options.lang = saved_lang
                self._result = saved_result
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import tran_reward_via_cpp
            return tran_reward_via_cpp(self, rewards_dict, name)
        if self._result is None:
            self._ensureAvgResults()
        return self._compute_tran_reward_named(name, rewards_dict)

    def reward_matrix_over(self, space, rewards_dict, nstates=None):
        """The (nrewards x nstates) matrix of the declared rewards on `space`.

        The reward map is a function of the AGGREGATE state row and of nothing
        else, so the same evaluation serves whichever engine produced the space:
        the native transient below reads it off `self._result`, and the
        lang='cpp' path reads it off line-cli's `labelsAggr`. Keeping one
        evaluation is what makes a bare callable answer identically under both,
        since no wire format can carry the callable itself.
        """
        from ...lang.reward_state import RewardState

        sn = self._sn

        # Build node-to-station and class-to-index mappings for RewardState
        nodes_to_station = {}
        if hasattr(self.model, 'get_nodes'):
            for node in self.model.get_nodes():
                if hasattr(node, 'get_index'):
                    node_idx = node.get_index()
                elif hasattr(node, 'index'):
                    node_idx = node.index
                else:
                    continue
                if hasattr(sn, 'nodeToStation'):
                    node_idx0 = node_idx - 1
                    if node_idx0 < len(sn.nodeToStation):
                        station_idx = int(sn.nodeToStation[node_idx0])
                        if station_idx >= 0:
                            nodes_to_station[node_idx] = station_idx + 1
        classes_to_idx = {}
        if hasattr(self.model, 'get_classes'):
            for i, jobclass in enumerate(self.model.get_classes()):
                if hasattr(jobclass, 'get_index'):
                    class_idx = jobclass.get_index()
                elif hasattr(jobclass, 'index'):
                    class_idx = jobclass.index
                else:
                    class_idx = i + 1
                classes_to_idx[class_idx] = i + 1

        # Build a reward vector over the state space for each named reward
        import inspect
        names = list(rewards_dict.keys())
        if nstates is None:
            nstates = len(space)
        Rmat = np.zeros((len(names), nstates))
        for ri, (nm, reward_fn) in enumerate(rewards_dict.items()):
            try:
                takes_sn = len(inspect.signature(reward_fn).parameters) >= 2
            except (TypeError, ValueError):
                takes_sn = False
            for s in range(min(nstates, len(space))):
                reward_state = RewardState(space[s], sn, nodes_to_station, classes_to_idx)
                try:
                    Rmat[ri, s] = reward_fn(reward_state, sn) if takes_sn else reward_fn(reward_state)
                except Exception:
                    pass
        return names, Rmat

    def _compute_tran_reward_named(self, name, rewards_dict):
        """Evaluate named model rewards over the aggregated state space and
        integrate them against the CTMC transient distribution to produce the
        E[r(X(t))] trajectories. Shared by get_tran_reward across
        lang='python'/'java'."""
        from scipy.linalg import expm

        space = self._result.space_aggr if hasattr(self._result, 'space_aggr') \
            and self._result.space_aggr is not None else self._result.space
        Q = self._result.infgen
        if space is None or len(space) == 0 or Q is None:
            raise ValueError('No CTMC state space available for transient reward analysis.')

        nstates = Q.shape[0]
        names, Rmat = self.reward_matrix_over(space, rewards_dict, nstates)

        # Time grid (mirror MATLAB: use timestep when provided, else 51 points)
        timespan = self.options.timespan
        t0, t1 = float(timespan[0]), float(timespan[1])
        timestep = getattr(self.options, 'timestep', None)
        if timestep and timestep > 0:
            tgrid = np.arange(t0, t1 + timestep / 2.0, timestep)
        else:
            tgrid = np.linspace(t0, t1, 51)

        # Transient distribution from the initial (empty) state at each time
        pi0 = np.zeros(nstates)
        pi0[0] = 1.0
        metrics = np.zeros((len(names), len(tgrid)))
        for k, tk in enumerate(tgrid):
            pit = pi0 @ expm(Q * tk)
            for ri in range(len(names)):
                metrics[ri, k] = float(pit @ Rmat[ri])

        Rt = [{'t': tgrid, 'metric': metrics[ri], 'name': names[ri]}
              for ri in range(len(names))]

        if name is not None:
            if name not in names:
                raise ValueError('Reward "%s" not found. Available rewards: %s'
                                 % (name, ', '.join(names)))
            idx = names.index(name)
            return Rt[idx], tgrid, names[idx]
        return Rt, tgrid, names

    # Aliases for new methods
    GetGenerator = getGenerator
    GetStateSpaceAggr = getStateSpaceAggr
    GetProb = getProb
    GetCdfSysRespT = getCdfSysRespT
    GetReward = getReward
    GetAvgReward = getAvgReward
    GetTranCdfRespT = getTranCdfRespT
    GetTranProb = getTranProb
    GetTranProbAggr = getTranProbAggr
    GetTranProbSys = getTranProbSys
    GetTranProbSysAggr = getTranProbSysAggr
    GetSymbolicGenerator = getSymbolicGenerator
    GetMarkedCTMC = getMarkedCTMC
    RunRewardAnalyzer = runRewardAnalyzer
    GetTranReward = getTranReward

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
                if len(chains_arr) == 0:
                    return [[k] for k in range(self._sn.nclasses)]
                nchains = max(nchains, int(np.max(chains_arr)) + 1)
                chains = [[] for _ in range(nchains)]
                for k in range(self._sn.nclasses):
                    if k < len(chains_arr):
                        c = int(chains_arr[k])
                        if 0 <= c < nchains:
                            chains[c].append(k)
                chains = [c for c in chains if c]
                return chains if chains else [[k for k in range(self._sn.nclasses)]]
            else:
                # 2D format: chains[c, k] > 0 means class k is in chain c
                chains = []
                for c in range(nchains):
                    chain_classes = []
                    for k in range(self._sn.nclasses):
                        if c < chains_arr.shape[0] and k < chains_arr.shape[1]:
                            if chains_arr[c, k] > 0:
                                chain_classes.append(k)
                    chains.append(chain_classes)
                chains = [c for c in chains if c]
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
        """Get average response times aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()

        R = self._result.R
        chains = self._get_chains()
        nstations = R.shape[0]
        nchains = len(chains)

        RN_chain = np.zeros((nstations, nchains))
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
        from ...api.sn import NodeType

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

        # Cache nodes prefer runAnalyzer's actualhitprob from its CTMC departure-rate measurement (more accurate); else fall back to cache analysis methods.
        if sn.nodeparam is not None:
            from ...api.cache import cache_prob_fpi, cache_ttl_lrua, cache_gamma_lp
            from ...lang.base import ReplacementStrategy

            for ind in range(I):
                if sn.nodetype is not None and ind < len(sn.nodetype):
                    if sn.nodetype[ind] == NodeType.CACHE and ind in sn.nodeparam:
                        cache_param = sn.nodeparam[ind]

                        # Check if actualhitprob was already computed by _compute_cache_hit_miss_probs()
                        # If so, skip the cache analysis to preserve the more accurate CTMC-based values
                        existing_hitprob = getattr(cache_param, 'actualhitprob', None)
                        if existing_hitprob is not None and np.any(existing_hitprob > 0):
                            # actualhitprob was already computed from CTMC departure rates - keep it
                            continue

                        nitems = getattr(cache_param, 'nitems', 0)
                        cap = getattr(cache_param, 'cap', 0)

                        # Get nitems and cap from model if not in param
                        if nitems == 0 or cap == 0:
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                                for node in self.model._nodes:
                                    if hasattr(node, '_nitems') and hasattr(node, '_capacity'):
                                        nitems = node._nitems if node._nitems else 0
                                        cap = node._capacity if node._capacity else 0
                                        break

                        if nitems > 0 and cap > 0:
                            # Get replacement strategy
                            replacement = ReplacementStrategy.LRU
                            if hasattr(cache_param, 'replacestrat') and cache_param.replacestrat is not None:
                                rs = cache_param.replacestrat
                                if isinstance(rs, int):
                                    replacement = rs
                                elif hasattr(rs, 'value'):
                                    replacement = rs.value
                                else:
                                    replacement = int(rs)

                            # Get cache node
                            cache_node = None
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                                cache_node = self.model._nodes[ind]

                            # Default hit/miss probabilities
                            hit_prob = min(cap / nitems, 1.0)
                            miss_prob = 1.0 - hit_prob

                            try:
                                # Build cache analysis parameters
                                m_levels = cache_node._item_level_cap if cache_node and hasattr(cache_node, '_item_level_cap') else np.array([cap])
                                h = len(m_levels)
                                pread = getattr(cache_param, 'pread', None)

                                # Find source station and get arrival rates
                                source_rate = np.zeros(R)
                                for ist in range(sn.nstations):
                                    src_ind = int(sn.stationToNode[ist])
                                    if src_ind < len(sn.nodetype):
                                        src_type = int(sn.nodetype[src_ind].value) if hasattr(sn.nodetype[src_ind], 'value') else int(sn.nodetype[src_ind])
                                        if src_type == 0:  # SOURCE
                                            source_rate = sn.rates[ist, :].copy()
                                            source_rate = np.nan_to_num(source_rate, nan=0.0)
                                            break

                                # Build lambda matrix (R x n x h+1)
                                lambd = np.zeros((R, nitems, h + 1))
                                for v in range(R):
                                    if pread is not None and v < len(pread) and pread[v] is not None:
                                        pread_v = np.asarray(pread[v]).ravel()
                                        for k in range(min(nitems, len(pread_v))):
                                            for l in range(h + 1):
                                                lambd[v, k, l] = source_rate[v] * pread_v[k]

                                # Default routing matrix
                                def create_default_routing(h_size):
                                    mat = np.diag(np.ones(h_size), 1)
                                    mat[h_size, h_size] = 1.0
                                    return mat

                                Rcost = getattr(cache_param, 'accost', None)
                                if Rcost is None:
                                    Rcost = [[create_default_routing(h) for _ in range(nitems)] for _ in range(R)]

                                # Compute gamma
                                gamma, _, _, _, _ = cache_gamma_lp(lambd, Rcost)

                                # Choose algorithm based on replacement strategy
                                if replacement in (ReplacementStrategy.RR, ReplacementStrategy.FIFO):
                                    pij = cache_prob_fpi(gamma, m_levels)
                                elif replacement == ReplacementStrategy.LRU:
                                    pij = cache_ttl_lrua(lambd, Rcost, m_levels)
                                else:
                                    pij = cache_prob_fpi(gamma, m_levels)

                                # Compute miss rates per class
                                miss_rate = np.zeros(R)
                                for v in range(R):
                                    if pread is not None and v < len(pread) and pread[v] is not None:
                                        pread_v = np.asarray(pread[v]).ravel()
                                        for k in range(min(nitems, len(pread_v), pij.shape[0])):
                                            miss_rate[v] += source_rate[v] * pread_v[k] * pij[k, 0]

                                # Compute overall hit/miss probabilities
                                total_rate = np.sum(source_rate)
                                total_miss = np.sum(miss_rate)
                                if total_rate > 0:
                                    miss_prob = total_miss / total_rate
                                    hit_prob = 1.0 - miss_prob
                            except Exception:
                                pass  # Fall back to uniform if cache analysis fails

                            # Store probabilities
                            hitclass = getattr(cache_param, 'hitclass', np.array([]))
                            nclasses = len(hitclass) if hasattr(hitclass, '__len__') else R
                            cache_param.actualhitprob = np.zeros(nclasses)
                            cache_param.actualmissprob = np.zeros(nclasses)

                            for k in range(nclasses):
                                hc = int(hitclass[k]) if k < len(hitclass) else -1
                                missclass = getattr(cache_param, 'missclass', np.array([]))
                                mc = int(missclass[k]) if k < len(missclass) else -1
                                if hc >= 0 and mc >= 0:
                                    cache_param.actualhitprob[k] = hit_prob
                                    cache_param.actualmissprob[k] = miss_prob

                            # Set on model
                            if cache_node is not None:
                                if hasattr(cache_node, 'set_result_hit_prob'):
                                    cache_node.set_result_hit_prob(cache_param.actualhitprob)
                                if hasattr(cache_node, 'set_result_miss_prob'):
                                    cache_node.set_result_miss_prob(cache_param.actualmissprob)

        # Create TH (throughput handle) - indicates which station-classes have valid throughput
        TH = np.zeros_like(TN)
        TH[TN > 0] = 1.0

        # Compute node arrival rates and throughputs using helper functions
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH, AN)
        TNn = sn_get_node_tput_from_tput(sn, TN, TH, ANn)

        # Initialize other node-level metrics
        QNn = np.zeros((I, R))
        UNn = np.zeros((I, R))
        RNn = np.zeros((I, R))
        WNn = np.zeros((I, R))

        # Compute residence times from response times using visit ratios
        WN = sn_get_residt_from_respt(sn, RN, None)

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

    def getTranAvg(self, *args):
        """Get transient average metrics.

        Computes time-dependent queue lengths, utilizations, and throughputs
        using transient CTMC analysis with matrix exponential method.

        Supports state prior iteration: when the model has multiple possible
        initial states (e.g., uniform prior from initFromMarginal + setStatePrior),
        runs the transient analysis for each state weighted by its prior probability.
        Matches MATLAB SolverCTMC/runAnalyzer.m lines 114-183.

        Args:
            *args: Optional transient handles (Qt, Ut, Tt) for MATLAB API compatibility.

        Returns:
            Tuple of (QNt, UNt, TNt) where each is a nested list [M][K] of TranResult objects.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import tran_avg_via_jar
            return tran_avg_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import tran_avg_via_cpp
            # A state prior over several rows crosses INTACT: line-cli seeds
            # `init_state_distribution`, the product of the declared per-node
            # priors over the enumerated space, and integrates once from the
            # mixture -- which is what the weighted sum below computes term by
            # term. See cpp_dispatch._assert_default_state for why the sampling
            # arms still cannot take one.
            return tran_avg_via_cpp(self)
        from ...constants import TranResult
        from ...api.mc.ctmc import ctmc_transient

        # see _kb/06-solver-catalog.md (CTMC section, transient methods)
        _tran_method, _fau_eps, _fau_delta = self._tran_fau_settings()

        # init state from nodes; native struct lacks sn.state/stateprior/space (MATLAB: getState/initDefault); else nodes read empty, fall to stationary.
        if hasattr(self.model, 'has_init_state') and not self.model.has_init_state():
            self.model.init_default()

        if self._result is None:
            self._ensureAvgResults()

        # Get dimensions
        M = self._sn.nstations
        K = self._sn.nclasses

        # Get required data from result
        infgen = self._result.infgen
        state_space = self._result.space
        state_space_aggr = self._result.space_aggr
        depRates = self._result.depRates
        n_states = infgen.shape[0]

        # Determine time span
        if self.options.timespan is not None:
            t_start, t_end = self.options.timespan[0], self.options.timespan[1]
        else:
            rates = np.asarray(self._sn.rates).flatten()
            finite_rates = rates[np.isfinite(rates) & (rates > 0)]
            if len(finite_rates) > 0:
                minrate = np.min(finite_rates)
                t_end = 30.0 / minrate
            else:
                t_end = 100.0
            t_start = 0.0

        # time-varying: rate_sched scales (station,class) transitions by m(t)=rate(t)/nominal; MATLAB solver_ctmc_transient_analyzer/local_ctmc_timevarying.
        rate_sched = self._tran_rate_sched()

        if rate_sched:
            time_points, qhat, mtraj, mscale = self._ctmc_timevarying_setup(
                infgen, rate_sched, t_start, t_end, M, K)
            n_times = len(time_points)
        else:
            # Generate time points
            timestep = self.options.timestep if hasattr(self.options, 'timestep') and self.options.timestep else 0.1
            n_points = max(2, int((t_end - t_start) / timestep) + 1)
            time_points = np.linspace(t_start, t_end, n_points)
            n_times = len(time_points)
            qhat, mtraj, mscale = None, None, None

        # Get per-node state spaces and priors from station objects
        # Matches MATLAB: s0 = sn.space; s0prior = sn.stateprior;
        per_node_spaces = []  # per_node_spaces[i] = (isf, state_space_array)
        per_node_priors = []  # per_node_priors[i] = prior_array
        stateful_nodes = []   # list of (node_index, isf)

        nodeToStateful = np.asarray(self._sn.nodeToStateful).flatten() \
            if hasattr(self._sn, 'nodeToStateful') else np.array([])

        for ind in range(self._sn.nnodes):
            if hasattr(self._sn, 'isstateful') and self._sn.isstateful[ind]:
                isf = int(nodeToStateful[ind])
                stateful_nodes.append((ind, isf))

                # Get per-node state space from station object
                node = self.model._nodes[ind]
                node_space = getattr(node, '_state_space', None)
                if node_space is not None and len(node_space) > 0:
                    node_space = np.atleast_2d(node_space)
                else:
                    # Fallback: single state from node's current state
                    node_state = node.get_state() if hasattr(node, 'get_state') else None
                    if node_state is not None:
                        node_space = np.atleast_2d(np.asarray(node_state).flatten())
                    else:
                        node_space = np.zeros((1, K))

                per_node_spaces.append(node_space)

                # Get per-node state prior
                node_prior = getattr(node, '_state_prior', None)
                if node_prior is not None:
                    node_prior = np.asarray(node_prior).flatten()
                    # Ensure prior length matches state space rows
                    if len(node_prior) < node_space.shape[0]:
                        padded = np.zeros(node_space.shape[0])
                        padded[:len(node_prior)] = node_prior
                        node_prior = padded
                else:
                    # Default prior: probability 1 on first state
                    node_prior = np.zeros(node_space.shape[0])
                    node_prior[0] = 1.0

                per_node_priors.append(node_prior)

        # Iterate over all state combinations (matching MATLAB pprod loop)
        # Each combination is a tuple of per-node state indices
        n_nodes = len(stateful_nodes)
        sizes = [s.shape[0] for s in per_node_spaces]

        # Initialize accumulated results
        QNt_accum = np.zeros((n_times, M, K))
        UNt_accum = np.zeros((n_times, M, K))
        TNt_accum = np.zeros((n_times, M, K))
        first_result = True

        # pprod-style enumeration of all state combinations
        s0_id = [0] * n_nodes
        while True:
            # Compute joint prior probability
            s0prior_val = 1.0
            for i in range(n_nodes):
                s0prior_val *= per_node_priors[i][s0_id[i]]

            if s0prior_val > 0:
                # Build full state vector by concatenating per-node states
                # Matches MATLAB: state = [state, zeros(pad), sn.state{isf}]
                full_state = []
                for i in range(n_nodes):
                    node_state = per_node_spaces[i][s0_id[i]]
                    full_state.extend(node_state)
                full_state = np.array(full_state, dtype=np.float64)

                # Find matching state in global state space
                state0_idx = self._find_state_in_full_space(full_state, state_space)

                if state0_idx == -1:
                    # Try with zero padding (for FCFS buffer width mismatches)
                    state0_idx = self._find_state_in_full_space_padded(
                        per_node_spaces, s0_id, state_space)

                if state0_idx >= 0:
                    # Build initial distribution
                    pi0 = np.zeros(n_states)
                    pi0[state0_idx] = 1.0

                    # Compute transient probabilities
                    if qhat is None:
                        pit = ctmc_transient(infgen, pi0, time_points,
                                             method=_tran_method, epsilon=_fau_eps,
                                             delta=_fau_delta)
                    else:
                        pit = self._ctmc_tv_propagate(infgen, pi0, time_points, qhat, mtraj)
                    if pit.ndim == 1:
                        pit = pit.reshape(1, -1)
                    pit = np.maximum(pit, 0.0)

                    # Compute metrics for this initial state
                    QNt_raw, UNt_raw, TNt_raw = self._compute_tran_metrics(
                        pit, state_space_aggr, depRates, n_times, n_states, M, K)
                    if mscale is not None:
                        TNt_raw = TNt_raw * mscale

                    # Weight by prior and accumulate
                    QNt_accum += s0prior_val * QNt_raw
                    UNt_accum += s0prior_val * UNt_raw
                    TNt_accum += s0prior_val * TNt_raw
                    first_result = False

            # Advance to next state combination (pprod)
            carry = True
            for i in range(n_nodes - 1, -1, -1):
                if carry:
                    s0_id[i] += 1
                    if s0_id[i] >= sizes[i]:
                        s0_id[i] = 0
                    else:
                        carry = False
                        break
            if carry:
                break  # All combinations exhausted

        # default transient start = EMPTY state for open/mixed (not stationary, flat curve); closed keeps steady-state fallback (states equal population).
        if first_result:
            _njobs = np.asarray(self._sn.njobs, dtype=float) if getattr(self._sn, 'njobs', None) is not None else np.array([])
            _has_open = bool(_njobs.size and np.any(np.isinf(_njobs)))
            _ssa = np.asarray(state_space_aggr)
            if _has_open and _ssa.ndim == 2 and _ssa.shape[0] == n_states and _ssa.size > 0:
                pi0 = np.zeros(n_states)
                pi0[int(np.argmin(_ssa.sum(axis=1)))] = 1.0
            else:
                _carried = getattr(self._result, 'pi0', None)
                if _carried is not None and np.asarray(_carried).size == n_states:
                    # seed carried through stochastic complementation by the analyzer:
                    # a vanishing initial state maps to its first-entry distribution
                    pi0 = np.asarray(_carried, dtype=float).flatten()
                else:
                    warnings.warn(
                        "CTMC transient: the declared initial state could not be located "
                        "in the enumerated state space, so the analysis starts from the "
                        "STATIONARY distribution and every curve is constant. The result "
                        "is a steady state reported as a transient, not a transient.",
                        UserWarning)
                    pi0 = self._result.pi.flatten()
            if qhat is None:
                pit = ctmc_transient(infgen, pi0, time_points,
                                     method=_tran_method, epsilon=_fau_eps,
                                     delta=_fau_delta)
            else:
                pit = self._ctmc_tv_propagate(infgen, pi0, time_points, qhat, mtraj)
            if pit.ndim == 1:
                pit = pit.reshape(1, -1)
            pit = np.maximum(pit, 0.0)
            QNt_accum, UNt_accum, TNt_accum = self._compute_tran_metrics(
                pit, state_space_aggr, depRates, n_times, n_states, M, K)
            if mscale is not None:
                TNt_accum = TNt_accum * mscale

        # Convert to TranResult format
        QNt = [[None for _ in range(K)] for _ in range(M)]
        UNt = [[None for _ in range(K)] for _ in range(M)]
        TNt = [[None for _ in range(K)] for _ in range(M)]

        for ist in range(M):
            for k in range(K):
                QNt[ist][k] = TranResult(time_points, QNt_accum[:, ist, k])
                UNt[ist][k] = TranResult(time_points, UNt_accum[:, ist, k])
                TNt[ist][k] = TranResult(time_points, TNt_accum[:, ist, k])

        return QNt, UNt, TNt

    def _tran_fau_settings(self):
        """
        The transient method and its tolerances from options.config.

        `transient_method` is a config key rather than a solver method name
        because it changes no stationary answer -- it is the transient path
        only -- and because a new entry in listValidMethods is enumerated by the
        sanity harness, which then demands a recorded baseline per method.

        Returns (method, epsilon, delta) with method in {'expm', 'fau'}.
        """
        cfg = getattr(self.options, 'config', None)

        def _read(name, default):
            if cfg is None:
                return default
            val = cfg.get(name, default) if isinstance(cfg, dict) else getattr(cfg, name, default)
            return default if val is None else val

        method = str(_read('transient_method', 'ode')).lower()
        if method not in ('ode', 'fau'):
            raise ValueError("Unknown options.config.transient_method '%s'; "
                             "use 'ode' or 'fau'." % method)
        # 'ode' selects the matrix exponential this analyzer has always used;
        # the name follows the MATLAB config, whose default branch integrates.
        return ('fau' if method == 'fau' else 'expm',
                float(_read('fau_epsilon', 1e-6)),
                float(_read('fau_delta', 1e-12)))

    def _tran_rate_sched(self):
        """The options.config['rate_sched'] entries, or an empty list.

        Each entry is a mapping (or object) with fields station, class (or
        jobclass), tgrid, rates and optionally nominal, using the same schema as
        the fluid rate multiplier (solver_fld.utils.ratemult).
        """
        cfg = getattr(self.options, 'config', None)
        if cfg is None:
            return []
        sched = cfg.get('rate_sched') if isinstance(cfg, dict) else getattr(cfg, 'rate_sched', None)
        if sched is None:
            return []
        return list(sched)

    @staticmethod
    def _rate_sched_field(entry, name, default=None):
        """Read a field of a rate_sched entry given as a dict or as an object."""
        if isinstance(entry, dict):
            if name in entry:
                return entry[name]
            if name == 'class' and 'jobclass' in entry:
                return entry['jobclass']
            return default
        if hasattr(entry, name):
            return getattr(entry, name)
        if name == 'class' and hasattr(entry, 'jobclass'):
            return getattr(entry, 'jobclass')
        return default

    @staticmethod
    def _scale_proc_entry(proc_ir, f):
        """Time-scale a service-process representation by the factor ``f``.

        Time-scaling a MAP/PH by ``f`` multiplies every rate by ``f``, leaving
        the embedded phase-selection probabilities untouched. The CTMC handler
        accepts the compact dict forms emitted by refreshStruct ({'rate'},
        {'k','mu'}, {'probs','rates'}) as well as the expanded matrix forms
        [D0, D1] and [alpha, T]; each is scaled in its own representation so no
        probability vector is corrupted.
        """
        if isinstance(proc_ir, dict):
            scaled = dict(proc_ir)
            if 'rate' in scaled and scaled['rate'] is not None:
                scaled['rate'] = float(scaled['rate']) * f
            if 'mu' in scaled and scaled['mu'] is not None:
                scaled['mu'] = float(scaled['mu']) * f
            if 'rates' in scaled and scaled['rates'] is not None:
                scaled['rates'] = np.asarray(scaled['rates'], dtype=float) * f
            return scaled
        if isinstance(proc_ir, (list, tuple)):
            elems = [np.atleast_2d(np.asarray(d, dtype=float)) if d is not None else None
                     for d in proc_ir]
            # [alpha,T] layout: alpha probability vector, not rate-scaled; other layouts are rate-valued. Discriminant mirrors CTMC handler's entry expansion.
            if (len(elems) == 2 and elems[0] is not None and elems[1] is not None
                    and elems[0].shape[0] == 1
                    and elems[1].shape[0] == elems[1].shape[1]):
                return [elems[0], elems[1] * f]
            return [e * f if e is not None else None for e in elems]
        return proc_ir

    def _ctmc_timevarying_setup(self, qbase, rate_sched, t_start, t_end, M, K):
        """Build the time grid, generator components and multipliers of the
        time-inhomogeneous transient.

        Mirrors MATLAB ``local_ctmc_timevarying``. The generator is linear in
        ``sn.rates[ist, r]``, so the component attributable to a scaled
        (station, class) is extracted by a single probe rebuild:
        ``Qhat_sc = (Q(scaled) - Qbase) / (probe - 1)``, and
        ``Q(t) = Qbase + sum_sc (m_sc(t) - 1) Qhat_sc``.

        ``solver_ctmc`` builds transitions from the process representation
        (``sn.proc`` MAP D0/D1 and ``sn.mu``), NOT from ``sn.rates``, so all
        rate-carrying fields must be scaled: time-scaling a MAP/PH by ``f``
        multiplies D0 and D1 (and the phase rates mu) by ``f``.

        Returns:
            (time_points, qhat, mtraj, mscale) where qhat is a list of
            (n_states x n_states) arrays, mtraj is (ngrid x nsched) and mscale
            is (ngrid x M x K), the per-(station,class) throughput multiplier.
        """
        from ...api.solvers.ctmc.handler import solver_ctmc as _solver_ctmc_handler

        cfg = getattr(self.options, 'config', None)
        ngrid = 100
        if cfg is not None:
            val = cfg.get('ctmc_tv_ngrid') if isinstance(cfg, dict) \
                else getattr(cfg, 'ctmc_tv_ngrid', None)
            if val:
                ngrid = int(val)
        ngrid = max(2, ngrid)
        time_points = np.linspace(t_start, t_end, ngrid)

        sn = self._sn
        qbase = np.asarray(qbase, dtype=float)
        n_states = qbase.shape[0]
        probe = 2.0
        nsc = len(rate_sched)
        qhat = []
        mtraj = np.ones((ngrid, nsc))
        mscale = np.ones((ngrid, M, K))

        handler_options = self._ctmc_handler_options()

        for s, entry in enumerate(rate_sched):
            ist = int(self._rate_sched_field(entry, 'station'))
            r = int(self._rate_sched_field(entry, 'class'))
            if not (0 <= ist < M and 0 <= r < K):
                raise ValueError(
                    "rate_sched entry refers to a (station,class) outside the network.")

            # probe rebuild time-scales a (station,class) service; reachable state space rate-independent so ordering matches qbase; rate fields restored after.
            saved_rate = float(sn.rates[ist, r])
            saved_proc = None
            if getattr(sn, 'proc', None) is not None and ist < len(sn.proc) \
                    and sn.proc[ist] is not None and r < len(sn.proc[ist]):
                saved_proc = sn.proc[ist][r]
            saved_mu = None
            if getattr(sn, 'mu', None) is not None and ist < len(sn.mu) \
                    and sn.mu[ist] is not None and r < len(sn.mu[ist]):
                saved_mu = sn.mu[ist][r]

            try:
                sn.rates[ist, r] = saved_rate * probe
                if saved_proc is not None:
                    sn.proc[ist][r] = self._scale_proc_entry(saved_proc, probe)
                if saved_mu is not None:
                    sn.mu[ist][r] = np.asarray(saved_mu, dtype=float) * probe
                probe_result = _solver_ctmc_handler(sn, handler_options)
                qp = np.asarray(probe_result.infgen, dtype=float)
            finally:
                sn.rates[ist, r] = saved_rate
                if saved_proc is not None:
                    sn.proc[ist][r] = saved_proc
                if saved_mu is not None:
                    sn.mu[ist][r] = saved_mu

            if qp.shape[0] != n_states:
                raise ValueError(
                    "rate_sched probe changed the CTMC state-space size; "
                    "cannot build time-varying generator.")
            qhat.append((qp - qbase) / (probe - 1.0))

            # multiplier m(t) = rate(t)/nominal (nominal defaults to sn.rates)
            nominal = self._rate_sched_field(entry, 'nominal')
            nominal = float(nominal) if nominal else saved_rate
            if nominal == 0.0:
                raise ValueError(
                    "rate_sched entry has a zero nominal rate; cannot form the multiplier.")
            seg_t = np.asarray(self._rate_sched_field(entry, 'tgrid'), dtype=float).ravel()
            seg_r = np.asarray(self._rate_sched_field(entry, 'rates'), dtype=float).ravel()
            clamped = np.clip(time_points, seg_t[0], seg_t[-1])
            mtraj[:, s] = np.interp(clamped, seg_t, seg_r) / nominal
            mscale[:, ist, r] = mtraj[:, s]

        return time_points, qhat, mtraj, mscale

    def _ctmc_handler_options(self):
        """The handler options used to (re)build the generator, matching the
        construction in runAnalyzer so the probe rebuild is consistent."""
        from ...api.solvers.ctmc.handler import SolverCTMCOptions as HandlerOptions
        return HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            cutoff=self.options.cutoff,
            verbose=self.options.verbose,
            force=self.options.force,
            gen_method=getattr(self.options, 'gen_method', 'default'),
            # the wall-clock budget and the state cap bound the solve only if they reach the handler.
            timeout=getattr(self.options, 'timeout', float('inf')),
            ctmc_max_states=getattr(self.options, 'ctmc_max_states', 3_000_000),
            memory_safety_fraction=getattr(self.options, 'memory_safety_fraction', 0.6),
        )

    @staticmethod
    def _ctmc_tv_propagate(qbase, pi0, time_points, qhat, mtraj):
        """Integrate dpi/dt = pi Q(t) over the grid, one segment per interval.

        The multiplier is frozen at the interval midpoint (second-order in dt),
        so over a segment the generator Qk is constant and the exact forward
        solution is pi(t+dt) = pi(t) expm(Qk dt).

        The matrix exponential is used deliberately rather than uniformization:
        an LN layer can carry a near-instantaneous reply-signal sentinel rate
        (~1e9), making the generator stiff (q dt ~ 1e8); uniformization then
        splits the step into more than 1e5 sub-segments and leaks all
        probability mass to zero. expm of a valid generator is exactly
        stochastic, so mass is conserved for any stiffness.
        """
        from scipy.linalg import expm as _expm

        qb = np.asarray(qbase, dtype=float)
        if hasattr(qb, 'toarray'):
            qb = qb.toarray()
        nt = len(time_points)
        pit = np.zeros((nt, qb.shape[0]))
        cur = np.asarray(pi0, dtype=float).ravel()
        pit[0, :] = cur
        for k in range(nt - 1):
            dt = time_points[k + 1] - time_points[k]
            qk = qb.copy()
            for s in range(len(qhat)):
                mk = 0.5 * (mtraj[k, s] + mtraj[k + 1, s])
                if mk != 1.0:
                    qk = qk + (mk - 1.0) * qhat[s]
            cur = cur @ _expm(qk * dt)
            pit[k + 1, :] = cur
        return pit

    def _compute_tran_metrics(self, pit, state_space_aggr, depRates,
                              n_times, n_states, M, K):
        """Compute transient QLen, Util, Tput from transient probabilities.

        Matches MATLAB solver_ctmc_transient_analyzer.m lines 55-98.
        Uses vectorized operations for efficiency.

        Args:
            pit: Transient probabilities (n_times x n_states)
            state_space_aggr: Aggregated state space (n_states x M*K)
            depRates: Departure rates (n_states x M x K)
            n_times: Number of time points
            n_states: Number of states
            M: Number of stations
            K: Number of classes

        Returns:
            Tuple of (QNt_raw, UNt_raw, TNt_raw) each (n_times x M x K)
        """
        QNt_raw = np.zeros((n_times, M, K))
        UNt_raw = np.zeros((n_times, M, K))
        TNt_raw = np.zeros((n_times, M, K))

        for ist in range(M):
            if hasattr(self._sn, 'nservers'):
                nservers_val = self._sn.nservers[ist, 0]
                nservers = int(nservers_val) if np.isfinite(nservers_val) else 1000000
            else:
                nservers = 1
            sched = self._sn.sched[ist] if hasattr(self._sn, 'sched') else None

            # Source's unbounded population column Inf, so aggregated QNt/UNt are undefined (pit@Inf=NaN); mirrors steady-state path leaving Source QN/UN at 0.
            is_source = False
            if hasattr(self._sn, 'stationToNode') and hasattr(self._sn, 'nodetype'):
                ind = int(self._sn.stationToNode[ist])
                nt = self._sn.nodetype[ind]
                nt_val = int(nt.value) if hasattr(nt, 'value') else int(nt)
                src_val = int(NodeType.SOURCE.value) if hasattr(NodeType.SOURCE, 'value') else int(NodeType.SOURCE)
                is_source = (nt_val == src_val)

            for k in range(K):
                col_idx = ist * K + k

                # Queue length: QNt = pit * SSA[:, col] (vectorized)
                if not is_source:
                    QNt_raw[:, ist, k] = pit @ state_space_aggr[:, col_idx]

                # Throughput: TNt = pit * depRates[:, ist, k]
                if depRates is not None:
                    TNt_raw[:, ist, k] = pit @ depRates[:, ist, k]

                if is_source:
                    continue

                # Utilization
                sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                if sched_name == 'INF':
                    UNt_raw[:, ist, k] = QNt_raw[:, ist, k]
                elif sched_name == 'PS':
                    total_at_station = np.sum(
                        state_space_aggr[:, ist * K:(ist + 1) * K], axis=1)
                    n_k = state_space_aggr[:, col_idx]
                    # an empty station's utilization is masked rather than computed via np.where, since np.where still evaluates the unselected 0/0 branch.
                    busy = total_at_station > 0
                    uik = np.zeros_like(total_at_station, dtype=float)
                    np.divide(np.minimum(n_k, nservers) * n_k, total_at_station,
                              out=uik, where=busy)
                    UNt_raw[:, ist, k] = pit @ (uik / nservers)
                else:
                    # FCFS, HOL, etc.
                    uik = np.minimum(state_space_aggr[:, col_idx], nservers) / nservers
                    UNt_raw[:, ist, k] = pit @ uik

        return QNt_raw, UNt_raw, TNt_raw

    def _find_state_in_full_space(self, full_state, state_space):
        """Find exact match of full state vector in global state space.

        Args:
            full_state: Full state vector (concatenation of per-node states)
            state_space: Global state space matrix (n_states x state_width)

        Returns:
            Index of matching row, or -1 if not found
        """
        if state_space is None or len(full_state) == 0:
            return -1

        full_state = np.asarray(full_state).flatten()

        if full_state.shape[0] != state_space.shape[1]:
            return -1

        # Find exact matching row
        for i in range(state_space.shape[0]):
            if np.allclose(state_space[i], full_state, atol=1e-10):
                return i
        return -1

    def _find_state_in_full_space_padded(self, per_node_spaces, s0_id, state_space):
        """Find state in global state space with zero-padding for width mismatches.

        For FCFS stations, the per-node state from fromMarginal may have a different
        buffer width than the global state space. This method pads per-node states
        with zeros on the left to match the global width.

        Matches MATLAB: state = [state, zeros(1,size(sn.space{isf},2)-length(sn.state{isf})), sn.state{isf}]

        Args:
            per_node_spaces: List of per-node state space matrices
            s0_id: List of per-node state indices
            state_space: Global state space matrix

        Returns:
            Index of matching row, or -1 if not found
        """
        if state_space is None:
            return -1

        global_width = state_space.shape[1]
        n_nodes = len(per_node_spaces)

        # First, compute the total width of per-node states
        per_node_states = [per_node_spaces[i][s0_id[i]] for i in range(n_nodes)]
        total_width = sum(len(s) for s in per_node_states)

        if total_width == global_width:
            # No padding needed, direct concatenation
            full_state = np.concatenate(per_node_states)
            return self._find_state_in_full_space(full_state, state_space)

        if total_width > global_width:
            return -1  # States wider than global space, can't match

        # Distribute the extra width among nodes (pad with zeros on the left)
        # This handles FCFS buffer width mismatches
        pad_total = global_width - total_width

        # state-vector width padding: pad the widest gap first across nodes.
        for pad_node_idx in range(n_nodes):
            full_state = []
            for i in range(n_nodes):
                node_state = per_node_states[i]
                if i == pad_node_idx:
                    # Pad with zeros on the LEFT (matching MATLAB)
                    padded = np.concatenate([np.zeros(pad_total), node_state])
                    full_state.extend(padded)
                else:
                    full_state.extend(node_state)

            full_state = np.array(full_state, dtype=np.float64)
            if len(full_state) == global_width:
                idx = self._find_state_in_full_space(full_state, state_space)
                if idx >= 0:
                    return idx

        return -1

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T) where R is chain-level system response time and
            T is chain-level system (carried) throughput.
        """
        CN, XN = self._computeChainMetrics()
        return CN, XN

    getAvgSysTable = NetworkSolver.getAvgSysTable  # chain-level shared layout

    # =========================================================================
    # SAMPLING METHODS STUBS (CTMC is analytical, redirect to SSA)
    # =========================================================================

    def sampleAggr(self, node: int, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated states at node (not supported for CTMC).

        Raises:
            NotImplementedError: CTMC is an analytical solver
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # `--node` narrows the same walk to that node's own block and its
            # per-class counts, which is the view MATLAB's per-node sampler
            # returns and which the native Network path does not implement.
            from ..cpp_dispatch import sample_path_via_cpp
            seed = getattr(self.options, 'seed', None)
            out = sample_path_via_cpp(self, numEvents, node=node,
                                      seed=seed if seed and int(seed) > 0 else None)
            return SampleResult(handle='ctmc', t=out['t'], state=out['nodeAggr'],
                                event=[], isaggregate=True, nodeIndex=node,
                                numEvents=len(out['state']))
        raise NotImplementedError("sampleAggr() not supported for analytical CTMC solver. Use SSA instead.")

    def sampleSys(self, numEvents: int = 1000) -> np.ndarray:
        """Sample system states.

        In chain mode this returns a sample path of the user-supplied chain,
        started from options.init_sol when given and from the uniform
        distribution otherwise; a DTMC advances one unit of time per step. For a
        Network model CTMC is an analytical solver and sampling is refused.

        Raises:
            NotImplementedError: CTMC is an analytical solver on a Network model
        """
        if self.isChainSolver():
            return self._chain_sample_sys(numEvents)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            # THE C++ HAS THIS AND THE NATIVE NETWORK PATH DOES NOT. MATLAB's
            # @SolverCTMC/sampleSys walks the chain with an exponential clock and
            # the port carries it (`-a sample`), so lang='cpp' answers a getter
            # that refuses natively rather than relaying the refusal. The
            # contract is _chain_sample_sys's: state rows, not state indices.
            from ..cpp_dispatch import sample_path_via_cpp
            seed = getattr(self.options, 'seed', None)
            out = sample_path_via_cpp(self, numEvents,
                                      seed=seed if seed and int(seed) > 0 else None)
            return SampleResult(handle='ctmc', t=out['t'], state=out['stateRows'],
                                event=[], isaggregate=False, numEvents=len(out['state']))
        raise NotImplementedError("sampleSys() not supported for analytical CTMC solver. Use SSA instead.")

    def _chain_sample_sys(self, numEvents: int) -> SampleResult:
        """Sample path of the user-supplied chain."""
        from ...api.mc import ctmc_simulate, dtmc_simulate

        self._ensureAvgResults()
        space = self._result.space
        pi0 = self._chain_init_distribution()
        seed = getattr(self.options, 'seed', None)
        rng = np.random.default_rng(seed)
        init_state = int(rng.choice(len(pi0), p=pi0))

        if self.isDiscreteChain():
            states = dtmc_simulate(self.getTransMat(), init_state, numEvents - 1, seed=seed)
            t = np.arange(numEvents, dtype=np.float64)
        else:
            Q = self._result.infgen
            # The Gillespie sampler stops at max_time, so leave it unbounded and
            # cap on the number of transitions instead.
            sim = ctmc_simulate(Q, init_state, np.inf, max_events=numEvents, seed=seed)
            states = np.asarray(sim['states'], dtype=int)[:numEvents]
            times = np.asarray(sim['times'], dtype=np.float64)[:numEvents]
            t = times
        states = np.asarray(states, dtype=int)[:numEvents]
        return SampleResult(handle='ctmc', t=np.asarray(t[:len(states)], dtype=np.float64),
                            state=space[states, :], event=[], isaggregate=False,
                            numEvents=len(states))

    def _queue_stateful_index(self) -> int:
        """Stateful index of the queue the sampled events are attributed to."""
        sn = self._sn
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else 1
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            from ...api.sn.network_struct import NodeType
            statefulToNode = sn.statefulToNode if hasattr(sn, 'statefulToNode') \
                else list(range(nstateful))
            for isf in range(nstateful):
                node_idx = int(statefulToNode[isf]) if isf < len(statefulToNode) else isf
                if node_idx < len(sn.nodetype) and sn.nodetype[node_idx] == NodeType.QUEUE:
                    return isf
        return 1  # default: the queue is the second stateful node

    def _events_from_trajectory(self, times, state_rows) -> List[EventInfo]:
        """ARV/DEP events of a sampled trajectory, from its population changes.

        A step that raises the total population is an arrival and one that lowers
        it a departure; a step that changes only a service phase is neither. A
        Source's infinite job-slot column is excluded, or every state reads as
        infinite population and no change is ever detected.
        """
        if times.size == 0 or state_rows.size == 0:
            return []
        pops = np.sum(np.where(np.isfinite(state_rows), state_rows, 0.0), axis=1)
        node = self._queue_stateful_index()
        events = []
        for i in range(1, min(len(pops), len(times))):
            change = pops[i] - pops[i - 1]
            if change > 0:
                events.append(EventInfo(node=node, jobclass=0, t=float(times[i]), event="ARV"))
            elif change < 0:
                events.append(EventInfo(node=node, jobclass=0, t=float(times[i]), event="DEP"))
        return events

    def sampleSysAggr(self, numEvents: int = 1000) -> SampleResult:
        """Sample aggregated system states using CTMC simulation.

        Uses the MMAP (Marked Markovian Arrival Process) approach matching
        MATLAB's sampleSysAggr. The CTMC generator is decomposed into event
        filter matrices (one per sync event) to build an MMAP, which is then
        sampled to produce exactly numEvents actual events (arrivals/departures).

        When event filtration is not available, falls back to direct CTMC
        simulation with enough transitions to produce numEvents actual events
        detected via population changes.

        Args:
            numEvents: Number of actual events (arrivals + departures) to generate

        Returns:
            SampleResult containing timestamps, states, and event information
        """
        from ...api.mc.ctmc import ctmc_simulate

        if getattr(self.options, 'lang', 'python') == 'cpp':
            # The aggregate view of the same walk: `-a sample` reports the
            # per-(station, class) counts along the trajectory beside the states,
            # both off one sample path, so the two views cannot come from two
            # different draws.
            from ..cpp_dispatch import sample_path_via_cpp
            seed = getattr(self.options, 'seed', None)
            seed = seed if seed and int(seed) > 0 else None
            # The ARV/DEP list is DERIVED from the trajectory, by the same
            # population-change rule the native branch below applies; leaving it
            # empty here made every caller that reads `.event` (the departure
            # process analyses) silently see no events under lang='cpp'.
            # `--samples` counts CTMC TRANSITIONS, while this method's contract
            # is numEvents EVENTS, so the request is re-scaled by the observed
            # event fraction rather than returning a short trajectory.
            request = int(numEvents)
            for _ in range(4):
                out = sample_path_via_cpp(self, request, seed=seed)
                times = np.asarray(out['t'], dtype=float).flatten()
                events = self._events_from_trajectory(
                    times, np.asarray(out['stateRows'], dtype=float))
                if len(events) >= numEvents or not events:
                    break
                fraction = len(events) / float(max(len(times) - 1, 1))
                request = int(numEvents / fraction * 1.2) + 100
            return SampleResult(handle='ctmc', t=out['t'], state=out['sysAggr'],
                                event=events[:numEvents], isaggregate=True,
                                numEvents=min(len(events), numEvents))

        # Run analyzer if needed
        if self._result is None:
            self._ensureAvgResults()

        # Get generator and state space
        infGen, eventFilt = self.getGenerator()
        stateSpace, _ = self.getStateSpace()

        if infGen is None or len(infGen) == 0:
            return SampleResult(isaggregate=True, numEvents=0)

        # Get network structure
        sn = self._sn

        # Build initial state index
        nstates = infGen.shape[0]
        initial_state = 0

        # Set random seed if provided
        seed = self.options.seed if hasattr(self.options, 'seed') else None
        if seed is not None:
            np.random.seed(seed)

        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else 1

        # per-state population excludes Source's infinite job-slot column; else states read infinite and event_fraction detects no change (loop hangs).
        state_populations = np.sum(np.where(np.isfinite(stateSpace), stateSpace, 0.0), axis=1)

        # Find which node corresponds to the queue (not source/sink)
        queue_node_idx = 1  # Default: queue is second stateful node
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            from ...api.sn.network_struct import NodeType
            statefulToNode = sn.statefulToNode if hasattr(sn, 'statefulToNode') else list(range(nstateful))
            for isf in range(nstateful):
                node_idx = int(statefulToNode[isf]) if isf < len(statefulToNode) else isf
                if node_idx < len(sn.nodetype):
                    if sn.nodetype[node_idx] == NodeType.QUEUE:
                        queue_node_idx = isf
                        break

        # fraction of CTMC transitions that are real events (vs phase transitions), from population-changing exit rate over total exit rate, per state.
        Q = np.asarray(infGen, dtype=np.float64)
        exit_rates = -np.diag(Q)
        # Build mask: population-changing transitions (arrivals/departures)
        pop_diff = (state_populations[:, None] != state_populations[None, :])
        Q_offdiag = Q.copy()
        np.fill_diagonal(Q_offdiag, 0.0)
        Q_offdiag[Q_offdiag < 0] = 0.0
        event_rates = np.sum(Q_offdiag * pop_diff, axis=1)

        # Compute steady-state weighted event fraction
        total_event_rate = np.sum(event_rates)
        total_exit_rate = np.sum(exit_rates)
        if total_exit_rate > 0:
            event_fraction = total_event_rate / total_exit_rate
        else:
            event_fraction = 1.0

        # Request enough transitions so we expect numEvents actual events
        # Add 20% safety margin to avoid needing multiple rounds
        if event_fraction > 0:
            n_transitions = int(numEvents / event_fraction * 1.2) + 100
        else:
            n_transitions = numEvents * 3

        # Simulate CTMC with enough transitions
        max_time = 1e10
        all_states = []
        all_times = []
        events = []
        current_time_offset = 0.0
        current_state = initial_state

        while len(events) < numEvents:
            remaining = numEvents - len(events)
            if event_fraction > 0:
                batch_size = int(remaining / event_fraction * 1.2) + 100
            else:
                batch_size = remaining * 3

            sim_result = ctmc_simulate(Q, current_state, max_time, batch_size, seed=None)

            states = sim_result['states']
            times = sim_result['times']

            # Offset times to continue from previous batch
            times = times + current_time_offset

            # Detect events from population changes
            for i in range(1, len(states)):
                prev_state_idx = states[i - 1]
                curr_state_idx = states[i]

                if prev_state_idx >= nstates or curr_state_idx >= nstates:
                    continue

                event_time = times[i]
                pop_change = state_populations[curr_state_idx] - state_populations[prev_state_idx]

                if pop_change > 0:
                    events.append(EventInfo(
                        node=queue_node_idx, jobclass=0,
                        t=event_time, event="ARV"
                    ))
                elif pop_change < 0:
                    events.append(EventInfo(
                        node=queue_node_idx, jobclass=0,
                        t=event_time, event="DEP"
                    ))

                if len(events) >= numEvents:
                    # Collect states/times up to this point
                    all_states.extend(states[:i + 1].tolist() if len(all_states) == 0 else states[1:i + 1].tolist())
                    all_times.extend(times[:i + 1].tolist() if len(all_times) == 0 else times[1:i + 1].tolist())
                    break
            else:
                # Used all transitions in this batch
                all_states.extend(states.tolist() if len(all_states) == 0 else states[1:].tolist())
                all_times.extend(times.tolist() if len(all_times) == 0 else times[1:].tolist())
                if len(states) > 0:
                    current_state = int(states[-1])
                    current_time_offset = float(times[-1])

        # Convert to arrays
        all_states = np.array(all_states, dtype=int)
        all_times = np.array(all_times, dtype=float)

        # Build result
        result = SampleResult(
            handle=f"ctmc_sample_{id(self)}",
            t=all_times.reshape(-1, 1) if len(all_times) > 0 else np.zeros((0, 1)),
            state=stateSpace[all_states] if len(all_states) > 0 else np.zeros((0, stateSpace.shape[1] if stateSpace.ndim > 1 else 1)),
            event=events,
            isaggregate=True,
            nodeIndex=None,
            numEvents=len(events)
        )

        return result

    # =========================================================================
    # Introspection Methods
    # =========================================================================

    def unsupportedMethodReason(self, method):
        """The forwarding address for the QRF reduction bounds, SolverBA's now.

        Asks nothing of the model, which is what lets the name gate in
        ``runAnalyzerChecks`` call it; ``runAnalyzer`` reads the text from here
        too, so the two cannot drift into two answers.
        """
        if not isinstance(method, str) or not method.startswith('qrf'):
            return ''
        return ("QRF bound method '%s' has moved out of SolverCTMC into the dedicated "
                "SolverBA solver. Use SolverBA(model, '%s') (or aliases 'qr'/'lr') instead."
                % (method, method))

    unsupported_method_reason = unsupportedMethodReason

    def listValidMethods(self) -> List[str]:
        """List valid solution methods.

        'exact' is an explicit alias for the default state-space path: it pins
        the intent at the call site so an example or test cannot be re-baselined
        by a later change of what 'default' selects. It must stay behaviourally
        identical to 'default'.

        'gpu' NAMES A BACKEND AND FALLS BACK, which is what the reference does:
        ctmc_solve.m wraps the gpuArray solve in a try/catch and runs the plain
        direct solve when no GPU is present, so SolverCTMC(model,'gpu') returns
        the exact answer on a host without one. This list used to name 'basic'
        instead -- a spelling no other codebase knows -- so 'gpu' was refused
        here and 'basic' was refused everywhere else.

        'mdd' holds the reachable set in a decision diagram and solves K coupled
        level-CTMCs instead of the ``|S|``-state generator; it is exact on
        product-form models and approximate otherwise, and is restricted to
        closed single-class networks (solver_ctmc_mdd_analyzer).
        """
        return ['default', 'exact', 'gpu', 'mdd', 'cftp', 'cftp.approx']

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features."""
        return {
            'Source', 'Sink',
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue', 'Router',
            'MAP', 'APH', 'MMPP2', 'MMAP', 'PH', 'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp', 'ME',
            'Det', 'Gamma', 'Weibull', 'Lognormal', 'Pareto', 'Uniform',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer', 'Buffer', 'Dispatcher',
            'Cache', 'CacheClassSwitcher', 'CacheRetrieval',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS',
            'SchedStrategy_DPS', 'SchedStrategy_GPS',
            'SchedStrategy_SIRO', 'SchedStrategy_SEPT',
            'SchedStrategy_LEPT', 'SchedStrategy_FCFS',
            'SchedStrategy_HOL', 'SchedStrategy_LCFS',
            'SchedStrategy_LCFSPR', 'SchedStrategy_LCFSPRPRIO', 'SchedStrategy_FCFSPRPRIO',
            # the rest of the preempt family: after_event_station carries one
            # arm for all eight, so declaring three gated five reachable
            # disciplines off at runAnalyzerChecks (matches SolverCTMC.m:244)
            'SchedStrategy_FCFSPR', 'SchedStrategy_LCFSPI', 'SchedStrategy_FCFSPI',
            'SchedStrategy_LCFSPIPRIO', 'SchedStrategy_FCFSPIPRIO',
            'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO',
            'SchedStrategy_LPS',
            'SchedStrategy_PAS',
            'SchedStrategy_OI',
            'SchedStrategy_POLLING',
            'RoutingStrategy_RROBIN',
            'RoutingStrategy_WRROBIN',
            'RoutingStrategy_JSQ',
            'RoutingStrategy_SQ',
            'RoutingStrategy_SDR',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',
            'ReplacementStrategy_HLRU', 'ReplacementStrategy_CLIMB', 'ReplacementStrategy_QLRU',
            'ClosedClass', 'SelfLoopingClass', 'OpenClass', 'Replayer',
            'OpenSignal', 'ClosedSignal',
            'SignalType_NEGATIVE', 'SignalType_CATASTROPHE', 'SignalType_REPLY',
            'SignalBatchRemoval', 'SignalRemovalPolicy',
            'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting', 'Timing', 'Firing', 'Storage',
            'Fork', 'Join', 'Forker', 'Joiner',
            'Balking', 'Reneging', 'Retrial', 'Breakdown',
            'LoadDependence',
            'ClassDependence',
            'JointDependence',
            'GlobalDependence',
            # FCR: handler filters state space for DROP, adds per-region FIFO for WAITQ; else featset gate rejects FCR models though handler solves it exactly.
            'Region',
            # c-server stations and binding buffers are both State constructs
            # (state_from_marginal / after_event_station): served by the
            # explicit generator, withdrawn from cftp and mdd.
            'MultiServer', 'FiniteCapacity',
        }

    def getMethodFeatureSet(self, method):
        """Per-method feature deltas applied to the base CTMC envelope.

        Four of the six methods share it; 'cftp'/'cftp.approx' and 'mdd' narrow
        it, because neither builds the explicit generator that carries the rest
        of the envelope. Mirrors MATLAB SolverCTMC.getMethodFeatureSet.
        """
        feats = set(SolverCTMC.getFeatureSet())
        if method in ('cftp', 'cftp.approx'):
            # PERFECT SAMPLING FROM A BALANCE FUNCTION, not from a generator:
            # the sampler encodes the closed single-class product form of
            # Gordon-Newell and nothing else, so every construct outside it has
            # to leave the envelope. The class count and the station count have
            # no registry name and are checked structurally in
            # supportsModelMethod, against the same predicate the analyzer uses.
            feats -= {
                'OpenClass',
                # Queue, Delay and Router are the only node kinds the sampler walks
                'Source', 'Sink', 'RandomSource', 'JobSink',
                'ClassSwitch', 'StatelessClassSwitcher',
                'Cache', 'CacheClassSwitcher', 'CacheRetrieval',
                'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
                'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',
                'ReplacementStrategy_HLRU', 'ReplacementStrategy_CLIMB',
                'ReplacementStrategy_QLRU',
                'Fork', 'Join', 'Forker', 'Joiner',
                'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting',
                'Timing', 'Firing', 'Storage',
                # disciplines outside INF/PS/FCFS/SIRO/LCFSPR have no product form
                'SchedStrategy_DPS', 'SchedStrategy_GPS',
                'SchedStrategy_SEPT', 'SchedStrategy_LEPT',
                'SchedStrategy_HOL', 'SchedStrategy_LCFS',
                'SchedStrategy_LCFSPRPRIO', 'SchedStrategy_FCFSPRPRIO',
                'SchedStrategy_FCFSPR', 'SchedStrategy_LCFSPI', 'SchedStrategy_FCFSPI',
                'SchedStrategy_LCFSPIPRIO', 'SchedStrategy_FCFSPIPRIO',
                'SchedStrategy_PSPRIO', 'SchedStrategy_DPSPRIO', 'SchedStrategy_GPSPRIO',
                'SchedStrategy_LPS', 'SchedStrategy_PAS', 'SchedStrategy_OI',
                'SchedStrategy_POLLING',
                # The one-phase-per-station rule is deliberately NOT spelled as
                # a list of distribution names. The sampler refuses
                # sn.phases[i, 0] > 1, and a name is not a phase count: a
                # one-phase Coxian passes and a HyperExp does not, while
                # Det/Gamma/Pareto only acquire their phases in
                # sn_nonmarkov_toph. supportsModelMethod asks the phase count
                # instead, which is also what lets it name the offending station.
                'Region',
                'LoadDependence', 'ClassDependence', 'JointDependence', 'GlobalDependence',
                # a state-dependent decision is not Markovian routing
                'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
                'RoutingStrategy_JSQ', 'RoutingStrategy_SQ', 'RoutingStrategy_SDR',
                # the Gordon-Newell balance function has no buffer:
                # solver_ctmc_cftp_supports refuses a finite one by name
                'FiniteCapacity',
            }
        elif method == 'mdd':
            # The decision diagram holds the MARKING of a closed network; an
            # open stream makes it unbounded, so there is no finite diagram to
            # hold. The single-class rule is structural (no registry name for a
            # class count) and lives in supportsModelMethod. A stochastic Petri
            # net keeps the Place/Transition names: spn_mdd reads the marking.
            #
            # A FORK-JOIN MODEL IS NEITHER of the two shapes it serves. The tag
            # augmentation a fork needs adds one auxiliary class per branch, so
            # the struct that reaches the analyzer is never single-class however
            # the model was written, and the level decomposition has no meaning
            # for a firing that does not conserve the per-chain population.
            feats -= {'OpenClass', 'Source', 'Sink', 'RandomSource', 'JobSink',
                      'Fork', 'Join', 'Forker', 'Joiner', 'JoinPartial',
                      # the level decomposition reads rates, servers and phases
                      # and no sn.cap/classcap, so a buffer would be dropped
                      'FiniteCapacity'}
        return feats

    def supportsModelMethod(self, method):
        """The per-method rules the feature registry has no name for, asked of
        the SAME predicates the analyzers use so that the report and the run
        cannot answer differently.

        Three of them: the class count and the station count that 'cftp' and
        'mdd' need (a class count is not a model feature), and the state-space
        size that the explicit-generator methods need. The last one is why
        'default'/'exact'/'gpu' were offered on models whose chain does not fit
        memory -- the analyzer priced the state space and refused, and nothing
        above it had asked. Mirrors MATLAB @SolverCTMC/supportsModelMethod.

        THE TWO STRUCTURAL PREDICATES ARE ASKED BEFORE THE FEATURE GATE, which
        is the reverse of the usual order and deliberate: each is the analyzer's
        own assert, so it refuses a strict superset of what the per-method
        feature deltas refuse, and its wording names the offending station or
        class count instead of a feature. Asking the feature gate first would
        replace 'the cftp method supports closed models only' with '(feature:
        OpenClass)' on the very run the caller is about to make.
        """
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'getStruct'):
            if method in ('cftp', 'cftp.approx'):
                from ...api.solvers.ctmc.solver_ctmc_cftp_analyzer import solver_ctmc_cftp_supports
                ok, reason = solver_ctmc_cftp_supports(model.getStruct(), self.options)
                if not ok:
                    return ok, reason
            elif method == 'mdd':
                from ...api.solvers.ctmc.solver_ctmc_mdd_analyzer import solver_ctmc_mdd_supports
                ok, reason = solver_ctmc_mdd_supports(model.getStruct())
                if not ok:
                    return ok, reason
            # The fork-join model class, which EVERY method has to clear: the
            # tag augmentation runs before the state space, the decision diagram
            # and the sampler alike, so a model sn_fj_validate refuses is
            # refused whichever name was asked for.
            ok, reason = sn_fj_supports(model.getStruct())
            if not ok:
                return ok, reason
        ok, reason = super().supportsModelMethod(method)
        if not ok or model is None or not hasattr(model, 'getStruct'):
            return ok, reason
        if method not in ('cftp', 'cftp.approx', 'mdd'):
            # The explicit state space is what the remaining methods enumerate,
            # and ctmc_memory_gate refuses it above the host budget. Asking the
            # same estimator here costs a combinatorial formula, not a state
            # space, so the report stays cheap.
            tractable, msg, _ = SolverCTMC.isStateSpaceTractable(model, self.options)
            if not tractable:
                return False, msg
        return True, reason

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Mirrors MATLAB SolverCTMC.supports: gates the model's used language
        features against getFeatureSet(). Struct-like inputs without a
        feature registry fall back to a structural sanity check.
        """
        if hasattr(model, 'get_used_lang_features') or hasattr(model, 'getUsedLangFeatures'):
            from ..base import SolverFeatureSet
            feat_used = (model.get_used_lang_features()
                         if hasattr(model, 'get_used_lang_features')
                         else model.getUsedLangFeatures())
            feat_supported = SolverFeatureSet()
            feat_supported.set_true(list(SolverCTMC.getFeatureSet()))
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
    def isStateSpaceTractable(model, options=None):
        """Whether the worst-case CTMC state space of ``model`` fits memory.

        Same estimator and gate the analyzer runs, exposed so a caller (e.g.
        SolverAUTO) can rank CTMC out before paying for state-space
        generation. Mirrors MATLAB SolverCTMC.isStateSpaceTractable and JAR
        SolverCTMC.isStateSpaceTractable.

        Args:
            model: the Network under analysis.
            options: solver options carrying cutoff, force and safety fraction.

        Returns:
            (ok, message, log_nstates).
        """
        from ...api.solvers.ctmc.memory_guard import (
            ctmc_memory_gate, state_space_log_size, DEFAULT_SAFETY_FRACTION)
        from ...api.sn import sn_nonmarkov_toph

        if options is None:
            options = SolverCTMC.defaultOptions()
        try:
            sn = model.getStruct() if hasattr(model, 'getStruct') else model
            # sn_nonmarkov_toph reads options as a mapping, not as the dataclass.
            cfg = options.get('config', {}) if isinstance(options, dict) else getattr(options, 'config', {})
            sn = sn_nonmarkov_toph(sn, {'config': cfg or {}})
            log_nstates = state_space_log_size(sn, options)
        except Exception as err:
            # An estimator failure must not be read as a refusal: the analyzer
            # runs its own gate and reports the real error.
            return True, str(err), 0.0
        force = bool(options.get('force', False) if isinstance(options, dict)
                     else getattr(options, 'force', False))
        safety = float(options.get('memory_safety_fraction', DEFAULT_SAFETY_FRACTION)
                       if isinstance(options, dict)
                       else getattr(options, 'memory_safety_fraction', DEFAULT_SAFETY_FRACTION))
        ok, msg = ctmc_memory_gate(log_nstates, force=force, verbose=False,
                                   safety_fraction=safety)
        return ok, msg, log_nstates

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'default',
            'tol': 1e-4,
            'cutoff': 10,
            'verbose': default_verbose(),
        })

    @staticmethod
    def printInfGen(infGen: np.ndarray, stateSpace: np.ndarray) -> None:
        """Print the infinitesimal generator matrix in MATLAB-compatible format.

        Output format matches MATLAB's CTMC.printInfGen():
        [from_state]->[to_state]: rate

        Args:
            infGen: Infinitesimal generator matrix
            stateSpace: State space matrix
        """
        if infGen is None or len(infGen) == 0:
            print("Empty generator matrix")
            return

        def format_state(state):
            """Format state vector as [a b c] with integers where possible."""
            parts = []
            for s in state:
                val = float(s)
                if val == int(val):
                    parts.append(str(int(val)))
                else:
                    parts.append(f"{val:.4f}")
            return "[" + " ".join(parts) + "]"

        nstates = infGen.shape[0]
        for i in range(nstates):
            state_i = stateSpace[i] if stateSpace is not None and i < len(stateSpace) else [i]
            for j in range(nstates):
                if i != j and infGen[i, j] > 0:
                    state_j = stateSpace[j] if stateSpace is not None and j < len(stateSpace) else [j]
                    from_str = format_state(state_i)
                    to_str = format_state(state_j)
                    print(f"{from_str}->{to_str}: {infGen[i, j]:.6f}")

    print_inf_gen = printInfGen

    @staticmethod
    def printEventFilt(eventFilt, SS, sync=None, events=None):
        """Print non-zero transitions per event in the event filter matrices.

        Output format matches MATLAB's SolverCTMC.printEventFilt() and
        JAR's SolverCTMC.printEventFilt().

        Args:
            eventFilt: List of event filter matrices (one per event).
            SS: State space matrix (nstates x state_dim).
            sync: Optional list of sync structures with active/passive node/class info.
            events: Optional list of event indices to print (1-based for MATLAB compat).
                    If None, prints all events.
        """
        if eventFilt is None or len(eventFilt) == 0:
            return

        def format_state(state):
            parts = []
            for s in state:
                val = float(s)
                if val == int(val):
                    parts.append(str(int(val)))
                else:
                    parts.append(f"{val:.4f}")
            return "[" + " ".join(parts) + "]"

        SS = np.asarray(SS)
        if SS.ndim == 1:
            SS = SS.reshape(-1, 1)

        if events is None:
            event_indices = range(len(eventFilt))
        else:
            event_indices = [e - 1 for e in events]  # Convert 1-based to 0-based

        for e in event_indices:
            if e < 0 or e >= len(eventFilt):
                continue
            D_e = np.asarray(eventFilt[e])
            if hasattr(D_e, 'toarray'):
                D_e = D_e.toarray()
            nstates = SS.shape[0]
            for s in range(nstates):
                for sp in range(nstates):
                    if D_e[s, sp] > 0:
                        from_str = format_state(SS[s, :])
                        to_str = format_state(SS[sp, :])
                        if sync is not None and e < len(sync):
                            se = sync[e]
                            act_node = se.get('active', [{}])[0].get('node', '?') if isinstance(se, dict) else getattr(getattr(se, 'active', [None])[0], 'node', '?')
                            act_cls = se.get('active', [{}])[0].get('class', '?') if isinstance(se, dict) else getattr(getattr(se, 'active', [None])[0], 'jobclass', '?')
                            pas_node = se.get('passive', [{}])[0].get('node', '?') if isinstance(se, dict) else getattr(getattr(se, 'passive', [None])[0], 'node', '?')
                            pas_cls = se.get('passive', [{}])[0].get('class', '?') if isinstance(se, dict) else getattr(getattr(se, 'passive', [None])[0], 'jobclass', '?')
                            print(f"{from_str}-- {e+1}: ({act_node},{act_cls}) => ({pas_node},{pas_cls}) -->{to_str}: {D_e[s, sp]:.6f}")
                        else:
                            print(f"Event {e}:")
                            print(f"  {from_str} -> {to_str} : {D_e[s, sp]:.6f}")

    print_event_filt = printEventFilt

    # =========================================================================
    # Sampling Methods (Not Supported - Analytical Solver)
    # =========================================================================

    def sample(self, node: int, numEvents: int) -> np.ndarray:
        """Sampling not supported by CTMC (analytical solver)."""
        raise NotImplementedError(
            "Sampling not supported by SolverCTMC. "
            "Use SolverSSA or SolverLDES for simulation-based analysis."
        )

    # =========================================================================
    # Aliases
    # =========================================================================

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
    GetStateSpace = getStateSpace
    GetSteadyState = getSteadyState
    GetInfGen = getInfGen
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

    # Sampling stubs
    SampleAggr = sampleAggr
    SampleSys = sampleSys
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
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_respt = getAvgRespT
    avg_resid_t = getAvgResidT
    avg_wait_t = getAvgWaitT
    avg_tput = getAvgTput
    avg_arv_r = getAvgArvR
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    avg_sys_table = getAvgSysTable
    avg_node = getAvgNode
    avg_node_table = getAvgNodeTable
    avg_node_chain = getAvgNodeChain
    avg_node_chain_table = getAvgNodeChainTable
    avg_chain = getAvgChain
    avg_chain_table = getAvgChainTable
    state_space = getStateSpace
    state_space_aggr = getStateSpaceAggr
    run_analyzer = runAnalyzer
    generator = getGenerator
    steady_state = getSteadyState
    sample_sys_aggr = sampleSysAggr
    sample_sys = sampleSys
    sample_aggr = sampleAggr
    inf_gen = getInfGen
    prob_aggr = getProbAggr
    prob_sys_aggr = getProbSysAggr
    prob = getProb
    tran_prob = getTranProb
    tran_prob_aggr = getTranProbAggr
    tran_prob_sys = getTranProbSys
    tran_prob_sys_aggr = getTranProbSysAggr


__all__ = ['SolverCTMC', 'SolverCTMCOptions']
