"""
Native Python implementation of Layered Network (LN) solver.

This implementation provides 100% parity with the MATLAB SolverLN implementation,
using the same layer-based decomposition with MVA solvers for each layer.

The architecture mirrors MATLAB's EnsembleSolver pattern:
1. Build layer submodels (Network objects) using buildLayersRecursive
2. Iterate until convergence using the EnsembleSolver pattern
3. Update metrics, think times, layers, and routing probabilities
4. Aggregate results using getEnsembleAvg

Pure Python implementation.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple, Callable, Union, Set
from dataclasses import dataclass, field
from ...constants import default_verbose
from enum import IntEnum
import copy
import os

# Import LINE network elements
from ...lang.network import Network
from ...lang.nodes import Queue, Delay, Source, Sink, Fork, Join, Router, Cache
from ...lang.classes import ClosedClass, OpenClass
from ...distributions import Exp, Immediate, Disabled
from ...constants import SchedStrategy, GlobalConstants
from ...lang.base import ReplacementStrategy
from ...api.sn.compat_rate import sn_compat_scaling
from ...api.io.logging import line_debug, line_warning
from ...layered import _call_count_dist
from ..base import EnsembleSolver


def _polyak_avg(prev, raw, k):
    """Running-mean update m = prev + (raw - prev)/k, robust to NaN entries in
    either operand (a NaN sample leaves the average untouched)."""
    if prev is None:
        return None if raw is None else np.array(raw, dtype=float, copy=True)
    prev = np.asarray(prev, dtype=float)
    if raw is None or np.shape(raw) != np.shape(prev):
        return prev
    raw = np.asarray(raw, dtype=float)
    m = prev + (raw - prev) / k
    bad = np.isnan(m)
    m[bad] = raw[bad]
    bad = np.isnan(m)
    m[bad] = prev[bad]
    return m


class LayeredNetworkElement(IntEnum):
    """Element types in layered queueing networks (matches MATLAB enum values)."""
    PROCESSOR = 0
    TASK = 1
    ENTRY = 2
    ACTIVITY = 3
    CALL = 4


class CallType(IntEnum):
    """Types of calls between entries."""
    SYNC = 1
    ASYNC = 2
    FWD = 3


from ...api.lqn import call_hashname, entry_workflow, ph_moments, serial_law
from ...distributions import APH, PH
from ...lang.workflow import Workflow

#: Method names that state the ensemble's LAYERING and ENCODING. They are the
#: vocabulary of SolverLN alone; a layer solver cannot dispatch on them.
_LN_LEVEL_METHODS = frozenset((
    'srvn', 'srvn.ph', 'srvn.cs', 'srvncs', 'ph', 'cs',
    'flat', 'flat.cs', 'flatcs', 'flat.ph', 'flatph', 'squashed', 'squashed.ph',
    'moment3',
))


def ln_requested_method(method) -> str:
    """Normalise a SolverLN method name onto one the solver dispatches on.

    A method name carries TWO decisions: the LAYERING, which fixes what a
    submodel is, and the ENCODING, which fixes how an activity graph is written
    into it. 'srvn.cs' encodes the activity graph as ROUTING, 'srvn.ph' as a
    composed phase-type server law, 'srvn' is the alias that takes 'srvn.ph'
    where it can serve the model and 'srvn.cs' otherwise, 'flat.cs' squashes
    every server into one submodel with the routing encoding, 'flat.ph' squashes
    them with the composed one ('flat' is the alias of 'flat.cs' and resolves
    unconditionally rather than probing 'flat.ph', because a model is squashed in
    order to express what only the routing encoding carries), and 'moment3' is
    the three-moment distribution pass over the routing layers. 'default' is the srvn alias, so a
    model solved without naming a method takes the better of the two srvn
    encodings; an unrecognised token takes 'srvn.cs'.
    """
    if not method or not isinstance(method, str):
        return 'srvn'
    m = method.lower()
    if m in ('srvn.ph', 'ph'):
        return 'srvn.ph'
    if m in ('srvn.cs', 'srvncs', 'cs'):
        return 'srvn.cs'
    if m in ('srvn', 'default', 'auto', ''):
        return 'srvn'
    if m in ('flat.cs', 'flatcs', 'flat', 'squashed'):
        return 'flat.cs'
    if m in ('flat.ph', 'flatph', 'squashed.ph'):
        return 'flat.ph'
    if m == 'moment3':
        return 'moment3'
    # An unrecognised token takes the routing encoding, which is what every name
    # other than 'moment3' resolved to before the alias existed.
    return 'srvn.cs'


# CallType, as the struct stores it: 1=SYNC, 2=ASYNC, 3=FWD
_SYNC = 1
_ASYNC = 2
_FWD = 3


def _alpha_of(dist) -> np.ndarray:
    """Initial vector of a phase-type distribution, as a row."""
    return np.asarray(dist.getInitProb(), dtype=float).reshape(1, -1)


def _subgen_of(dist) -> np.ndarray:
    """Subgenerator of a phase-type distribution: D0 of its (D0, D1) pair."""
    return np.asarray(dist.getRepresentation()[0], dtype=float)


class PHLayer:
    """One two-station layer, and the caller classes that cycle through it."""

    def __init__(self):
        self.idx: int = 0
        self.ishost: bool = False
        self.callers: List[int] = []
        self.class_of_caller: Dict[int, int] = {}
        self.nreplicas: int = 1
        self.qstations: List[int] = []
        self.svcmean_by_class: Dict[int, float] = {}
        # entries are (class index, entry index) or (class index, -call index)
        self.open_arrivals: List[Tuple[int, int]] = []
        # Closed population of the MODEL this server sits in. Under 'flat.ph'
        # that is every caller of the single network, not only the callers of
        # this one station, so it is recorded here rather than recomputed.
        self.npop: float = 0.0


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
    """Options for the native LN solver (matches MATLAB SolverLN.defaultOptions)."""
    method: str = 'default'
    iter_max: int = 200  # MATLAB default for LN
    iter_tol: float = 5e-3  # MATLAB default for LN (looser than default for LQN models)
    verbose: bool = field(default_factory=default_verbose)
    tol: float = 1e-4  # MATLAB SolverOptions LN case; JAR SolverOptions.java:379, cpp solver_ln.h:228
    seed: Optional[int] = None  # base seed for stochastic layer solvers; randomized when not set
    # Transient window [t0, t1] for getTranAvg; propagated to layer solvers only
    # around the transient call (SolverENV over an LQN stage sets this).
    timespan: Optional[Any] = None
    # 'python' (native), 'java' (delegate the layered solve to jline.jar via
    # JSON) or 'cpp' (delegate it to line-cli via the .lqnx interchange, since
    # the C++ port has no LQN JSON reader); env LINE_SOLVER_LANG overrides the
    # default.
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. The other langs are IEEE double throughout, so it is left
    # None and line-cli is invoked without --arith unless the caller sets it. The
    # layer solver must be MVA, since the C++ fluid layers are double-only. USE
    # 'real:<digits>' AND NOT 'exact' HERE: rational arithmetic grows the
    # coefficients unboundedly along an outer fixed point, so an exact layered
    # solve does not terminate in practice, while real:64 costs 4x double.
    arith: Optional[str] = None

    # Config options (matches MATLAB options.config)
    config: OptionsDict = field(default_factory=lambda: OptionsDict({
        'interlocking': True,
        # Layering strategy: 'srvn' places each server in its own submodel,
        # 'flat' places every processor and task in a single submodel
        'layering': 'srvn',
        'relax': 'fixed',  # 'none', 'fixed', 'adaptive', 'auto' - matches LQNS default
        'relax_factor': 0.5,  # under-relaxation factor
        'relax_min': 0.1,  # MATLAB default
        'relax_history': 5,  # MATLAB default
        # stochastic iteration options for simulation/MC layer solvers; see _kb/06-solver-catalog.md LN Convergence test.
        'stochiter': 'auto',        # 'auto' | 'rm' | 'crn' | 'off'
        'stochiter_alpha': 0.6,     # Robbins-Monro step decay exponent, in (0.5,1]
        'stochiter_a0': 1.0,        # Robbins-Monro initial step after burn-in
        'stochiter_burnin': 5,      # Picard burn-in iterations before step decay starts
        'stochiter_conseq': 3,      # consecutive sub-tolerance iterations required to stop
    }))


def _region_capable_layer_solver(model):
    """
    Pick the first solver whose feature set covers a layer carrying an admission
    constraint. The order is by decreasing accuracy: CTMC is exact but
    state-space bound, LDES and SSA simulate. Selection is by supports() so it
    self-corrects if another solver later declares Region.
    """
    from ..solver_ctmc.solver_ctmc import SolverCTMC
    from ..solver_ssa.solver_ssa import SolverSSA
    from ..wrappers.solver_ldes.solver_ldes import SolverLDES
    for ctor in (SolverCTMC, SolverLDES, SolverSSA):
        if ctor.supports(model):
            return ctor(model, verbose=False)
    raise ValueError(f"LN layer {model.getName()} carries an admission constraint but none of "
                     f"SolverCTMC, SolverLDES, SolverSSA supports it. Supply a layer solver "
                     f"factory explicitly.")


def _layer_dep_handle(f, cols, nclasses, layer_model):
    """
    Lift a service-rate dependence handle declared on a LayeredNetwork server to
    the layer station that represents it. F maps the per-operand population vector
    of that server to a scalar scaling shared by every operand or to a per-operand
    vector; COLS[j] lists the layer classes (1-based) through which operand j
    occupies the station.

    Solvers evaluate the handle in two different index spaces: CTMC and the exact
    recursions pass a per-class vector, while the AMVA and NC chain recursions
    pass a per-chain vector. The handle therefore reads len(n) to pick the space,
    aggregates the operand populations in it, and answers a vector of the SAME
    length, since the caller indexes the answer with the index it passed in. An
    index belonging to no operand keeps the neutral scaling 1.
    """
    chain_cols = {}

    def handle(n):
        n = np.atleast_1d(np.asarray(n, dtype=float)).ravel()
        idx = cols
        if n.size != nclasses:
            if n.size not in chain_cols:
                chain_cols[n.size] = _layer_chain_cols(cols, layer_model, n.size)
            idx = chain_cols[n.size]
        nop = np.array([float(np.sum(n[[c - 1 for c in idx[j]]])) if idx[j] else 0.0
                        for j in range(len(idx))])
        w = np.atleast_1d(np.asarray(f(nop), dtype=float)).ravel()
        v = np.ones(n.size)
        for j in range(len(idx)):
            for c in idx[j]:
                v[c - 1] = w[min(j, w.size - 1)]
        return v
    return handle


def _layer_chain_cols(cols, layer_model, nchains):
    """Operand columns of a layer station in the chain index space."""
    sn = layer_model.getStruct()
    chains = np.atleast_2d(np.asarray(sn.chains, dtype=float))
    out = []
    for cols_of_operand in cols:
        ch = set()
        for c in cols_of_operand:
            for k in np.where(chains[:, c - 1] > 0)[0]:
                if k + 1 <= nchains:
                    ch.add(int(k) + 1)
        out.append(sorted(ch))
    return out


def _layer_peak(peak_per_operand, cols, nclasses):
    """Spread a per-operand peak rate scaling onto the classes of the layer station."""
    peak_per_operand = np.atleast_1d(np.asarray(peak_per_operand, dtype=float)).ravel()
    peak = np.ones(nclasses)
    for j in range(len(cols)):
        pj = peak_per_operand[min(j, peak_per_operand.size - 1)]
        for c in cols[j]:
            peak[c - 1] = pj
    return peak


class SolverLN(EnsembleSolver):
    """
    Native Python Layered Network (LN) solver.

    This implementation matches MATLAB's SolverLN at 100% parity:
    - Uses the same layer decomposition algorithm (buildLayersRecursive)
    - Creates Network objects for each layer with proper classes and routing
    - Uses MVA solvers for each layer
    - Implements the same fixed-point iteration with convergence testing

    The algorithm::

        1. Build layer submodels: one per processor (host layer) and per task
        2. Initialize service demands and think times from LQN structure
        3. Iterate until convergence:
           a. Solve each layer using MVA
           b. Update service times based on lower-layer response times
           c. Update think times based on caller waiting times
           d. Update routing probabilities based on throughputs
           e. Check convergence
        4. Aggregate results from all layers

    ``options.lang`` delegates the whole layered solve instead: ``'java'`` to
    ``jline.jar`` over JSON, ``'cpp'`` to the C++ ``line-cli`` over the ``.lqnx``
    interchange (steady state only, and it refuses what that interchange cannot
    carry -- see ``solvers/cpp_dispatch.py``). Either way the native fixed point
    never runs and ``options.arith`` selects the C++ arithmetic backend.

    Args:
        model: LayeredNetwork model
        solver_factory: Optional factory function to create layer solvers
        options: Solver options
        **kwargs: Additional options
    """

    def __init__(self, model, solver_factory_or_options=None, options=None, **kwargs):
        self.model = model
        self._result = None

        # Parse options (matches MATLAB signature handling)
        self._parse_options(solver_factory_or_options, options, kwargs)

        # Layer structures (matches MATLAB SolverLN properties)
        self.ensemble: List[Network] = []  # Network objects for each layer
        self.solvers: List[Any] = []  # Solver instances for each layer
        self.nlayers: int = 0
        self.lqn = None  # LayeredNetworkStruct

        # Index mappings (matches MATLAB)
        self.idxhash: np.ndarray = None  # Maps LQN indices to layer indices
        self.hostLayerIndices: List[int] = []
        self.taskLayerIndices: List[int] = []

        # Job counts for interlocking
        self.njobs: np.ndarray = None
        self.njobsorig: np.ndarray = None

        # Update maps (populated by buildLayersRecursive)
        self.servt_classes_updmap: np.ndarray = None
        self.thinkt_classes_updmap: np.ndarray = None
        self.actthinkt_classes_updmap: np.ndarray = None
        self.arvproc_classes_updmap: np.ndarray = None
        self.call_classes_updmap: np.ndarray = None
        self.route_prob_updmap: np.ndarray = None
        self.unique_route_prob_updmap: np.ndarray = None

        # Reset indices
        self.routereset: List[int] = []
        self.svcreset: List[int] = []

        # Replication tracking (matches MATLAB singleReplicaTasks)
        self.single_replica_tasks: List[int] = []

        # Metric arrays
        self.util: np.ndarray = None
        self.tput: np.ndarray = None
        self.tputproc: List = None
        self.servt: np.ndarray = None
        self.residt: np.ndarray = None
        self.servtproc: List = None
        self.servtcdf: List = None
        self.thinkt: np.ndarray = None
        self.thinkproc: List = None
        self.thinktproc: List = None
        self.entryproc: List = None
        # method='moment3': True once the moment-based entry-law pass has run.
        self.moment_pass_done: bool = False
        self.entrycdfrespt: List = None
        self.callresidt: np.ndarray = None
        # per-entry service time resolved by the servtmatrix solve, kept for
        # inspection exactly as MATLAB's SolverLN keeps it
        self.entry_servt: np.ndarray = None
        self.callservt: np.ndarray = None
        self.callservtproc: List = None
        self.callservtcdf: List = None
        self.ignore: np.ndarray = None

        # Service matrix for entry service time calculation
        self.servtmatrix: np.ndarray = None

        # Caller probability tracking
        self.ptaskcallers: np.ndarray = None
        self.ptaskcallers_step: List = None
        self.ilscaling: np.ndarray = None

        # Interlock path tables of Franks (1999), Ch. 4 (built once at init)
        self.il_table_all: np.ndarray = None    # (nentries x nentries) reachability, all phases
        self.il_table_ph1: np.ndarray = None    # (nentries x nentries) reachability, phase-1 only
        self.il_common_entries: list = None      # common parent entry abs-indices per server
        self.il_source_tasks_all: list = None    # all-phase source tasks per server
        self.il_source_tasks_ph2: list = None    # phase-2 source tasks per server
        self.il_num_sources: np.ndarray = None   # total source multiplicity per server

        # Convergence tracking
        self.hasconverged: bool = False
        self.averagingstart: int = None
        self.maxitererr: List[float] = []
        self.results: List[List[Dict]] = []

        # Under-relaxation state
        self.relax_omega: float = 1.0
        self.relax_err_history: List[float] = []
        self.servt_prev: np.ndarray = None
        self.residt_prev: np.ndarray = None
        self.tput_prev: np.ndarray = None
        self.thinkt_prev: np.ndarray = None
        self.callservt_prev: np.ndarray = None
        self.callresidt_prev: np.ndarray = None

        # Stochastic iteration (Robbins-Monro / Polyak-Ruppert) state, used
        # when one or more layer solvers return noisy estimates
        self.stochiter_mode: str = None      # resolved mode: 'rm' | 'crn' | 'off'
        self.stochiter_auto: bool = False    # True if mode was resolved from 'auto'
        self.stochiter_start: int = None     # iteration at which RM averaging started
        self.stochiter_seed_base: int = None # base seed for layer seed control
        self.stochlayers: np.ndarray = None  # bool per layer: solver is stochastic
        self.stoch_avg: List[Dict] = None    # Polyak-Ruppert averages of layer results
        self.stoch_avg_count: int = 0        # iterations accumulated into stoch_avg
        self.stoch_servt_avg: np.ndarray = None   # Polyak-Ruppert average of the servt iterate
        self.stoch_residt_avg: np.ndarray = None  # Polyak-Ruppert average of the residt iterate


        # Phase-2 support
        self.hasPhase2: bool = False
        self.servt_ph1: np.ndarray = None
        self.servt_ph2: np.ndarray = None
        self.util_ph1: np.ndarray = None
        self.util_ph2: np.ndarray = None
        self.prOvertake: np.ndarray = None

        # Extract LQN structure and construct layers
        self._extract_lqn_structure()
        self._construct()

    def _parse_options(self, solver_factory_or_options, options, kwargs):
        """Parse options handling MATLAB-style signatures."""
        from ..solver_mva.solver_mva import SolverMVA

        self.solver_factory = None

        if solver_factory_or_options is None:
            # Check kwargs for method parameter
            method = kwargs.get('method', 'default')
            if isinstance(method, str):
                method = method.lower()
        elif callable(solver_factory_or_options) and not isinstance(solver_factory_or_options, SolverLNOptions):
            # a solver class (SolverMVA) is as valid a factory as a lambda, and
            # MATLAB's @SolverMVA maps onto the class, so both must forward the
            # third argument; excluding types dropped `options` silently
            self.solver_factory = solver_factory_or_options
            if options is not None:
                if hasattr(options, 'get') and not hasattr(options, 'method'):
                    method = options.get('method', 'default')
                else:
                    method = getattr(options, 'method', 'default')
                # An options object passed alongside a factory carries the same
                # fields as one passed on its own; forwarding only `method`
                # would silently drop config entries such as `layering`.
                # lang/arith are forwarded with the rest: a factory says which
                # solver runs each layer, not which engine runs the ensemble, so
                # dropping them would silently ignore a requested lang='cpp'.
                for _f in ('iter_max', 'iter_tol', 'verbose', 'tol', 'config', 'lang', 'arith'):
                    if hasattr(options, _f):
                        _v = getattr(options, _f)
                        if _v is not None:
                            kwargs.setdefault(_f, _v)
                    elif hasattr(options, 'get'):
                        _v = options.get(_f, None)
                        if _v is not None:
                            kwargs.setdefault(_f, _v)
            else:
                # honor a method passed as a keyword argument alongside a factory
                method = kwargs.get('method', 'default')
            if isinstance(method, str):
                method = method.lower()
        elif isinstance(solver_factory_or_options, str):
            method = solver_factory_or_options.lower()
        elif hasattr(solver_factory_or_options, 'get'):
            method = solver_factory_or_options.get('method', 'default')
            if 'verbose' in solver_factory_or_options:
                kwargs.setdefault('verbose', solver_factory_or_options['verbose'])
            if 'iter_max' in solver_factory_or_options:
                kwargs.setdefault('iter_max', solver_factory_or_options['iter_max'])
        elif isinstance(solver_factory_or_options, SolverLNOptions):
            # Handle SolverLNOptions dataclass - extract all relevant attributes
            method = solver_factory_or_options.method
            kwargs.setdefault('iter_max', solver_factory_or_options.iter_max)
            kwargs.setdefault('iter_tol', solver_factory_or_options.iter_tol)
            kwargs.setdefault('verbose', solver_factory_or_options.verbose)
            kwargs.setdefault('tol', solver_factory_or_options.tol)
            kwargs.setdefault('config', solver_factory_or_options.config)
            kwargs.setdefault('lang', solver_factory_or_options.lang)
            kwargs.setdefault('arith', solver_factory_or_options.arith)
        elif hasattr(solver_factory_or_options, 'method'):
            method = getattr(solver_factory_or_options, 'method', 'default')
        else:
            method = 'default'

        # NC does not handle LQN class-switching layer structure; always use MVA for layer solving.
        if self.solver_factory is None:
            if method == 'nc':
                import warnings
                warnings.warn(
                    "NC method for SolverLN is not fully supported in native Python. "
                    "Falling back to MVA for layer solving. Use the default method for "
                    "correct results.",
                    UserWarning
                )
            # MVA handles LQN layer models correctly; LN iter_tol=5e-3 is forwarded.
            self.solver_factory = lambda m: (
                _region_capable_layer_solver(m)
                if getattr(m.getStruct(), 'nregions', 0) > 0
                else SolverMVA(m, self._layer_options(), verbose=False)
            )
        elif isinstance(self.solver_factory, type):
            # a bare solver class carries no options of its own, so the layer
            # solver is given the LN options exactly as the default factory does
            _cls = self.solver_factory
            # kept because the lambda hides the class from the feature checks
            self._layer_solver_cls = _cls
            self.solver_factory = lambda m: _cls(m, self._layer_options(), verbose=False)

        kwargs.pop('method', None)
        self.options = SolverLNOptions(method=method, **kwargs)

    def _extract_lqn_structure(self):
        """Extract layered network structure from model."""
        if hasattr(self.model, 'getStruct'):
            self.lqn = self.model.getStruct()
        else:
            raise ValueError("Model must be a LayeredNetwork with getStruct() method")

        # Normalize structure format
        self._normalize_lqn_structure()

        # forwarding rewritten as caller-side pseudo rendezvous; see _kb/06-solver-catalog.md LN Forwarding as caller-side pseudo-rendezvous.
        self._apply_forwarding_rendezvous()

        # Detect and initialize phase-2 support (matches MATLAB SolverLN.m lines 127-137)
        if (hasattr(self.lqn, 'actphase') and self.lqn.actphase is not None
                and np.any(self.lqn.actphase > 1)):
            self.hasPhase2 = True
            self.servt_ph1 = np.zeros(self.lqn.nidx)
            self.servt_ph2 = np.zeros(self.lqn.nidx)
            self.util_ph1 = np.zeros(self.lqn.nidx)
            self.util_ph2 = np.zeros(self.lqn.nidx)
            self.prOvertake = np.zeros(self.lqn.nentries)
        else:
            self.hasPhase2 = False

    def _apply_forwarding_rendezvous(self):
        """Forwarding transformation of Franks (1999), Sec. 3.3.1 and Fig. 3.8.

        Each forwarding chain reachable from a synchronous call is reconnected
        to the client that issued the original rendezvous, as a pseudo
        rendezvous (SYNC) call whose mean is the original call mean times the
        product of the forwarding probabilities along the path. One level of
        servers disappears from the layering and the forwarded workload is
        carried by ordinary SYNC call classes, so layer construction, think
        times, populations and the interlock analysis all see plain rendezvous
        arcs. As the thesis notes, the pseudo arcs are excluded from the slice
        times and from the overtaking and interlock probabilities. FWD calls
        remain in the struct but no longer contribute blocking anywhere in
        SolverLN. Asynchronous calls into a forwarding chain are left
        untouched, since a send-no-reply terminates the chain of blocking."""
        lqn = self.lqn
        if lqn.ncalls == 0 or not np.any(np.asarray(lqn.calltype[:lqn.ncalls]) == CallType.FWD):
            return

        ncalls0 = lqn.ncalls
        for cidx in range(ncalls0):
            if int(lqn.calltype[cidx]) != CallType.SYNC:
                continue
            aidx = int(lqn.callpair[cidx, 0])
            tidx = self._get_parent(aidx)
            base_mean = self._get_call_mean(cidx)
            if base_mean is None or base_mean <= 0:
                continue
            # BFS through the forwarding chain of the sync target
            frontier = [int(lqn.callpair[cidx, 1])]
            probs = [1.0]
            visited_e = []
            while frontier:
                eidx = frontier.pop(0)
                p_path = probs.pop(0)
                if eidx in visited_e:
                    continue
                visited_e.append(eidx)
                for fcidx in range(ncalls0):
                    if int(lqn.calltype[fcidx]) != CallType.FWD or int(lqn.callpair[fcidx, 0]) != eidx:
                        continue
                    fprob = self._get_call_mean(fcidx)
                    tgt = int(lqn.callpair[fcidx, 1])
                    pseudo_mean = base_mean * p_path * (fprob or 0.0)
                    target_tidx = self._get_parent(tgt)
                    if pseudo_mean > 0 and target_tidx != tidx:
                        # merge into an existing SYNC call with the same (activity,target) pair, else append a new pseudo SYNC call.
                        # call indices are 0-based here, so 0 is a real call and cannot double as the not-found sentinel (MATLAB/C++ are 1-based and do use 0)
                        mrow = -1
                        for scan in range(lqn.ncalls):
                            if int(lqn.calltype[scan]) == CallType.SYNC \
                                    and int(lqn.callpair[scan, 0]) == aidx \
                                    and int(lqn.callpair[scan, 1]) == tgt:
                                mrow = scan
                                break
                        if mrow >= 0:
                            newmean = self._get_call_mean(mrow) + pseudo_mean
                            lqn.callpair[mrow, 2] = newmean
                            if isinstance(lqn.callproc, list) and mrow < len(lqn.callproc):
                                lqn.callproc[mrow] = _call_count_dist(newmean)
                            elif isinstance(lqn.callproc, dict):
                                lqn.callproc[mrow] = _call_count_dist(newmean)
                        else:
                            ncall = lqn.ncalls
                            lqn.ncalls = ncall + 1
                            newrow = np.zeros((1, lqn.callpair.shape[1]))
                            newrow[0, 0] = aidx
                            newrow[0, 1] = tgt
                            newrow[0, 2] = pseudo_mean
                            lqn.callpair = np.vstack([lqn.callpair, newrow])
                            lqn.calltype = np.append(lqn.calltype, CallType.SYNC)
                            if isinstance(lqn.callproc, list):
                                lqn.callproc.append(_call_count_dist(pseudo_mean))
                            elif isinstance(lqn.callproc, dict):
                                lqn.callproc[ncall] = _call_count_dist(pseudo_mean)
                            if hasattr(lqn, 'callsof') and isinstance(lqn.callsof, dict):
                                lqn.callsof.setdefault(aidx, []).append(ncall)
                            if hasattr(lqn, 'iscaller') and lqn.iscaller is not None:
                                lqn.iscaller[tidx, target_tidx] = 1
                                lqn.iscaller[aidx, target_tidx] = 1
                                lqn.iscaller[tidx, tgt] = 1
                                lqn.iscaller[aidx, tgt] = 1
                            if hasattr(lqn, 'issynccaller') and lqn.issynccaller is not None:
                                lqn.issynccaller[tidx, target_tidx] = 1
                                lqn.issynccaller[aidx, target_tidx] = 1
                                lqn.issynccaller[tidx, tgt] = 1
                                lqn.issynccaller[aidx, tgt] = 1
                            if hasattr(lqn, 'graph') and lqn.graph is not None:
                                lqn.graph[aidx, tgt] = 1
                            if hasattr(lqn, 'taskgraph') and lqn.taskgraph is not None:
                                lqn.taskgraph[tidx, target_tidx] = 1
                    # Follow the chain
                    if tgt not in visited_e and tgt not in frontier:
                        frontier.append(tgt)
                        probs.append(p_path * (fprob or 0.0))

    def _normalize_lqn_structure(self):
        """Normalize LQN structure to consistent format."""
        lqn = self.lqn

        # Convert numpy arrays to dicts for mapping attributes if needed
        for attr in ['tasksof', 'entriesof', 'actsof', 'callsof']:
            data = getattr(lqn, attr, None)
            if data is not None and isinstance(data, np.ndarray):
                result = {}
                for i in range(len(data)):
                    if data[i] is not None:
                        if hasattr(data[i], '__iter__') and not isinstance(data[i], str):
                            result[i + 1] = list(data[i])
                        else:
                            result[i + 1] = [data[i]] if data[i] else []
                setattr(lqn, attr, result)

        # Rebuild callsof from callpair if empty
        if isinstance(lqn.callsof, dict) and len(lqn.callsof) == 0:
            if hasattr(lqn, 'callpair') and lqn.callpair is not None:
                for cidx in range(lqn.ncalls):
                    if cidx < lqn.callpair.shape[0]:
                        src_aidx = int(lqn.callpair[cidx, 0])  # source activity in column 0
                        if src_aidx > 0:
                            if src_aidx not in lqn.callsof:
                                lqn.callsof[src_aidx] = []
                            lqn.callsof[src_aidx].append(cidx)

    def _act_thinktime(self, aidx):
        """Think time of activity aidx, 0.0 when it has none.

        In series with the activity's host demand and held at its task: the task
        keeps its thread for the whole hostdem+thinktime interval, so it
        serializes against the task multiplicity, but the host processor is
        released for it. Mirrors lqns, whose think-time attribute LINE writes out,
        and MATLAB lqn_act_thinktime. Read through actthinkproc, which the
        constructor has already filtered to genuinely positive durations.
        """
        proc = getattr(self, 'actthinkproc', None)
        if proc is None or aidx >= len(proc):
            return 0.0
        zt = proc[aidx]
        if zt is None:
            return 0.0
        try:
            v = zt.getMean()
        except AttributeError:
            return 0.0
        if v is None or not np.isfinite(v) or v <= 1e-8:
            return 0.0
        return float(v)

    def _construct(self):
        """Construct layer models (matches MATLAB construct method)."""
        lqn = self.lqn

        # Mark disconnected components to ignore
        # MATLAB SolverLN.construct lines 169-185: weaklyconncomp(graph'+graph)
        self.ignore = np.zeros(lqn.nidx, dtype=bool)
        if hasattr(lqn, 'graph') and lqn.graph is not None:
            graph = np.asarray(lqn.graph)
            n = graph.shape[0]
            # Undirected adjacency for weak connectivity
            symm = graph + graph.T
            symm = (symm > 0).astype(int)
            try:
                from scipy.sparse.csgraph import connected_components
                from scipy.sparse import csr_matrix
                n_components, labels = connected_components(
                    csr_matrix(symm[:n, :n]), directed=False)
            except ImportError:
                # Fallback: BFS-based connected components
                n_components, labels = 0, np.zeros(n, dtype=int)
                visited = np.zeros(n, dtype=bool)
                for start in range(n):
                    if not visited[start]:
                        queue = [start]
                        visited[start] = True
                        while queue:
                            node = queue.pop(0)
                            labels[node] = n_components
                            for nbr in range(n):
                                if symm[node, nbr] > 0 and not visited[nbr]:
                                    visited[nbr] = True
                                    queue.append(nbr)
                        n_components += 1

            if n_components > 1:
                # Find which components contain REF tasks
                wcc_has_ref = np.zeros(n_components, dtype=bool)
                for t in range(lqn.ntasks):
                    tidx = lqn.tshift + t
                    if tidx < n and self._is_ref_task(tidx):
                        wcc_has_ref[labels[tidx]] = True
                # Components with an entry-level open arrival are also workload-anchored
                if hasattr(lqn, 'arrival') and lqn.arrival:
                    for eidx in lqn.arrival:
                        if lqn.arrival[eidx] is not None and eidx < n:
                            wcc_has_ref[labels[eidx]] = True
                # Mark all elements in components without REF tasks as ignored
                for comp in range(n_components):
                    if not wcc_has_ref[comp]:
                        for idx in range(n):
                            if labels[idx] == comp and idx <= lqn.nidx:
                                self.ignore[idx] = True

        # Initialize internal data structures
        self.entrycdfrespt = [None] * lqn.nentries
        self.hasconverged = False
        self.moment_pass_done = False

        # Initialize service and think time processes
        self.servtproc = [None] * lqn.nidx
        self.thinkproc = [None] * lqn.nidx
        self.callservtproc = [None] * lqn.ncalls
        self.tputproc = [None] * lqn.nidx

        # prefer the full Distribution from lqn.hostdem_proc (keeps SCV/phase-type); fall back to the Exp-fitted scalar mean.
        hostdem_proc = getattr(lqn, 'hostdem_proc', None)

        def _servtproc_from(mean_or_dist, idx_abs):
            proc = None
            if isinstance(hostdem_proc, dict):
                proc = hostdem_proc.get(idx_abs)
            elif hostdem_proc is not None and idx_abs < len(hostdem_proc):
                proc = hostdem_proc[idx_abs]
            # zero host-demand activities map to Immediate; check the scalar mean before the fitted proc so a rate-0 proc cannot leak in as service rate 0.
            if isinstance(mean_or_dist, (int, float)) and float(mean_or_dist) <= 0:
                return Immediate()
            if proc is not None and not isinstance(proc, (int, float)):
                return proc
            if isinstance(mean_or_dist, (int, float)):
                return Exp.fit_mean(float(mean_or_dist))
            return mean_or_dist

        if isinstance(lqn.hostdem, dict):
            for idx, mean_or_dist in lqn.hostdem.items():
                if mean_or_dist is not None:
                    self.servtproc[idx] = _servtproc_from(mean_or_dist, idx)
        else:
            for idx in range(len(lqn.hostdem)):
                if lqn.hostdem[idx] is not None:
                    self.servtproc[idx + 1] = _servtproc_from(lqn.hostdem[idx], idx + 1)

        # Copy think times - convert floats to Exp distributions
        if isinstance(lqn.think, dict):
            for idx, mean_or_dist in lqn.think.items():
                if mean_or_dist is not None:
                    if isinstance(mean_or_dist, (int, float)):
                        mean_val = float(mean_or_dist)
                        if mean_val <= 0:
                            self.thinkproc[idx] = Immediate()
                        else:
                            self.thinkproc[idx] = Exp.fit_mean(mean_val)
                    else:
                        self.thinkproc[idx] = mean_or_dist
        else:
            for idx in range(len(lqn.think)):
                if lqn.think[idx] is not None:
                    mean_or_dist = lqn.think[idx]
                    if isinstance(mean_or_dist, (int, float)):
                        mean_val = float(mean_or_dist)
                        if mean_val <= 0:
                            self.thinkproc[idx + 1] = Immediate()
                        else:
                            self.thinkproc[idx + 1] = Exp.fit_mean(mean_val)
                    else:
                        self.thinkproc[idx + 1] = mean_or_dist

        # Copy activity think times - convert floats to Exp distributions
        self.actthinkproc = [None] * lqn.nidx
        if hasattr(lqn, 'actthink') and isinstance(lqn.actthink, dict):
            for idx, mean_or_dist in lqn.actthink.items():
                if mean_or_dist is not None:
                    if isinstance(mean_or_dist, (int, float)):
                        mean_val = float(mean_or_dist)
                        if mean_val > 1e-8:
                            self.actthinkproc[idx] = Exp.fit_mean(mean_val)
                    elif hasattr(mean_or_dist, 'getMean') and mean_or_dist.getMean() > 1e-8:
                        self.actthinkproc[idx] = mean_or_dist

        # entries have Immediate servtproc initially; entry service time (servt) is computed iteratively from activities during update_layers.
        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            # Set servtproc to Immediate for entries (matches MATLAB: hostdem{eidx} is empty for entries)
            self.servtproc[eidx] = Immediate()

        # call service time process = target entry's hostdem (Immediate for entries, matches MATLAB line 194-196).
        for cidx in range(lqn.ncalls):
            tgt_eidx = self._get_call_target_entry(cidx)
            if tgt_eidx is not None and tgt_eidx > 0 and tgt_eidx < len(self.servtproc):
                if self.servtproc[tgt_eidx] is not None:
                    self.callservtproc[cidx] = self.servtproc[tgt_eidx]
                else:
                    self.callservtproc[cidx] = Immediate()
            else:
                self.callservtproc[cidx] = Immediate()

        # Build entry service matrix (matches MATLAB getEntryServiceMatrix)
        # This matrix maps activities and calls to entries for computing entry service times
        self.servtmatrix = self._get_entry_service_matrix()

        # Initialize job counts
        self.njobs = np.zeros((lqn.tshift + lqn.ntasks, lqn.tshift + lqn.ntasks))

        # Build layers
        self._build_layers()

        # A setup no longer forces the MAM decomposition on the layer. The open
        # M/G/1-with-setup QBD reads the idle period from the Poisson rate 1/X, and
        # in a CLOSED layer the idle period a thread sees is the rest of the cycle,
        # 1/X - S: on lqn_setup that is 1.0 against the 2.29 the open reading gives,
        # so the thread was powered down far more often than it is and the answer
        # landed 12.67% below LDES. The cold start is charged to the ENTRY instead,
        # with the probability that the thread was actually found down: see
        # _setup_charge.

        self.njobsorig = self.njobs.copy()

        # Build the interlock path tables of Sec. 4.2
        if self.options.config.get('interlocking', False):
            self._init_interlock()

        self.nlayers = len(self.ensemble)
        line_debug("LN construct: built %d layers from LQN model (%d hosts, %d tasks, %d entries, %d activities)",
                   self.nlayers, lqn.nhosts, lqn.ntasks, lqn.nentries, lqn.nacts)

        # Initialize caller probability tracking
        self.ptaskcallers = np.zeros((lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks))
        self.ptaskcallers_step = [np.zeros_like(self.ptaskcallers) for _ in range(self.nlayers + 2)]

        # Compute reset indices (convert to int for list indexing)
        if self.route_prob_updmap is not None and len(self.route_prob_updmap) > 0:
            self.routereset = list(set(int(self.idxhash[int(x)]) for x in self.route_prob_updmap[:, 0]
                                       if not np.isnan(self.idxhash[int(x)])))
        if self.thinkt_classes_updmap is not None and len(self.thinkt_classes_updmap) > 0:
            self.svcreset = list(set(int(self.idxhash[int(x)]) for x in self.thinkt_classes_updmap[:, 0]
                                     if not np.isnan(self.idxhash[int(x)])))
        if self.call_classes_updmap is not None and len(self.call_classes_updmap) > 0:
            self.svcreset = list(set(self.svcreset) |
                                set(int(self.idxhash[int(x)]) for x in self.call_classes_updmap[:, 0]
                                    if not np.isnan(self.idxhash[int(x)])))

        # Store ensemble in model
        self.model.ensemble = self.ensemble

    def listValidMethods(self):
        """Valid methods for this solver, SolverLN.m verbatim.

        Each name states the LAYERING and the ENCODING; ln_requested_method
        normalises the alias spellings ('ph', 'cs', 'srvncs', 'flatcs',
        'squashed', 'squashed.ph') onto these, and they are left out here to
        keep the list unambiguous, exactly as the reference does.
        """
        return ['srvn', 'srvn.ph', 'srvn.cs', 'flat', 'flat.cs', 'flat.ph',
                'moment3', 'default']

    list_valid_methods = listValidMethods

    def supportsModelMethod(self, method):
        """The encoding rules the layer builders enforce at solve time, stated
        here so a CALLER can see them before running.

        'srvn.ph' and 'flat.ph' compose each entry into ONE phase-type law, and
        several constructs have nowhere to go in that law: a forwarding call
        whose target is not in the caller's activity graph, a routed call group
        whose dispatch order the composition folds away, a cache task, an
        admission constraint, a queue-dependent rate on a station the
        composition replaces. 'flat.ph' additionally squashes every layer into
        one network, which per-layer state (a replica, a powered-down setup
        thread) cannot survive.

        None of these is a feature name, so none can be a feature-set delta:
        they are properties of what the METHOD does to the model. Left only in
        the builders they were invisible to every gate above them, and
        ``listValidMethods`` returns the same eight names for every model, so a
        report offered every encoding on every layered model.

        Phase 2 is deliberately NOT tested: that refusal reads ``self.hasPhase2``,
        which is built during layering rather than being a property of the model,
        so a gate cannot ask it without doing the layering it precedes. Mirrors
        MATLAB ``ln_method_refusal``.
        """
        m = str(method).lower()
        if m not in ('srvn.ph', 'flat.ph'):
            return True, ''
        lqn = getattr(self, 'lqn', None)
        if lqn is None:
            return True, ''

        # -- the squashing refusals, 'flat.ph' only ------------------------
        # Each carries PER-LAYER state that one submodel cannot hold, so they
        # are properties of the flattening and not of the encoding.
        if m == 'flat.ph':
            nelem = lqn.nhosts + lqn.ntasks
            for i in range(nelem):
                if float(lqn.repl[0, i]) > 1:
                    return False, ("method='flat.ph' does not support replicated processors or "
                                   "tasks, whose replicas need a submodel each. "
                                   "Use method='srvn.ph'.")
            hs = getattr(lqn, 'hassetup', None)
            if hs is not None and np.any(np.asarray(hs).ravel()[:nelem]):
                return False, ("method='flat.ph' does not support setup tasks, whose "
                               "powered-down threads are per-layer state. "
                               "Use method='srvn.ph'.")

        # -- the composed-entry-law refusals, both PH encodings -------------
        iscache = getattr(lqn, 'iscache', None)
        if iscache is not None and np.any(np.asarray(iscache).ravel()):
            return False, ("method='%s' does not support cache tasks. "
                           "Use method='default'." % m)
        for cidx in range(lqn.ncalls):
            if self._ph_call_type(cidx) == _FWD:
                return False, ("method='%s' does not support forwarding calls, whose target is "
                               "not part of the caller's activity graph. "
                               "Use method='default'." % m)
        hs = getattr(lqn, 'hassetup', None)
        if hs is not None:
            hsf = np.asarray(hs).ravel()
            for i in range(len(hsf)):
                if not hsf[i]:
                    continue
                if self._get_sched(i) == SchedStrategy.INF or not np.isfinite(float(lqn.mult[0, i])):
                    return False, ("method='%s': task '%s' declares a setup time on an "
                                   "infinite-server task, which holds no thread to power down; "
                                   "give it a finite multiplicity." % (m, self._ph_name(i)))
        if getattr(lqn, 'callgroups', None):
            return False, ("method='%s' does not support routed call groups, whose dispatch "
                           "order is a routing property. Use method='flat.cs'." % m)
        if getattr(lqn, 'lincon', None):
            return False, ("method='%s' does not support admission constraints on a layer "
                           "station. Use method='default'." % m)
        for fndep in ('lldscaling', 'cdscaling', 'jdscaling', 'pools'):
            dep = getattr(lqn, fndep, None) or {}
            if dep:
                sidxdep = sorted(dep.keys())[0]
                what = 'server pools' if fndep == 'pools' else fndep
                return False, ("method='%s' does not support queue-dependent service rates on "
                               "a layer station ('%s' declares %s). Use method='srvn.cs'."
                               % (m, self._ph_name(sidxdep), what))
        return True, ''

    supports_model_method = supportsModelMethod

    def supports(self, model) -> bool:
        """Check if the layered model is supported.

        Mirrors MATLAB SolverLN.supports: an LQN is solved layer by layer, so
        the gate is the conjunction of the per-layer solvers' own gates against
        their own layer, not a feature set of SolverLN's own. This cannot be a
        static method, since it needs self.solvers[e].

        No supports() existed anywhere in the MRO, so calling it raised
        AttributeError and the solver had no gate at all.
        """
        # python LayeredNetwork.getEnsemble returns SELF, not a list of layers, unlike MATLAB's @LayeredNetwork.
        ensemble = self.ensemble
        if not ensemble and hasattr(model, 'ensemble'):
            candidate = model.ensemble
            if isinstance(candidate, (list, tuple)):
                ensemble = candidate
        if not ensemble or not self.solvers:
            # The layers are built by _construct(); with none built there is
            # nothing to gate against.
            return True

        from ..base import supports_via_featureset
        for e in range(min(len(ensemble), len(self.solvers))):
            solver = self.solvers[e]
            layer = ensemble[e]
            # SolverMAM/SolverFLD declare supports(sn,method)->(bool,reason), unlike boolean supports(model) MATLAB assumes; gate via the layer solver's featset.
            get_featureset = getattr(type(solver), 'getFeatureSet', None)
            if get_featureset is not None:
                if not supports_via_featureset(type(solver), layer):
                    return False
                continue
            supports = getattr(solver, 'supports', None)
            if supports is None:
                continue
            try:
                if not supports(layer):
                    return False
            except TypeError:
                # Solver with a non-model supports() signature; nothing to gate.
                continue
        return True

    def _get_call_target_entry(self, cidx: int) -> Optional[int]:
        """Get the target entry index for a call."""
        lqn = self.lqn
        if cidx < 0 or cidx >= lqn.ncalls:
            return None
        if isinstance(lqn.callpair, dict):
            pair = lqn.callpair.get(cidx, None)
            if pair is not None:
                return pair[1]  # Column 1 is target entry
        else:
            if cidx < lqn.callpair.shape[0]:
                return int(lqn.callpair[cidx, 1])  # Column 1 is target entry
        return None

    def _get_call_source_activity(self, cidx: int) -> Optional[int]:
        """Get the source activity index for a call."""
        lqn = self.lqn
        if cidx < 0 or cidx >= lqn.ncalls:
            return None
        if isinstance(lqn.callpair, dict):
            pair = lqn.callpair.get(cidx, None)
            if pair is not None:
                return pair[0]  # Column 0 is source activity
        else:
            if cidx < lqn.callpair.shape[0]:
                return int(lqn.callpair[cidx, 0])  # Column 0 is source activity
        return None

    def _assert_series_parallel_forks(self):
        """Reject activity graphs whose AND forks and joins are not properly nested.

        The traversal pairs a join with the most recent fork through a LIFO stack
        of fork classes, so it can only represent series-parallel graphs.
        """
        lqn = self.lqn
        graph = getattr(lqn, 'graph', None)
        posttype = getattr(lqn, 'actposttype', None)
        pretype = getattr(lqn, 'actpretype', None)
        if graph is None or posttype is None or pretype is None:
            return
        post_and_value = 12  # ActivityPrecedenceType.ID_POST_AND
        pre_and_value = 2    # ActivityPrecedenceType.ID_PRE_AND
        post = np.asarray(posttype).flatten()
        pre = np.asarray(pretype).flatten()
        ashift = lqn.nhosts + lqn.ntasks + lqn.nentries
        nidx = lqn.nidx
        acts = range(ashift, nidx)

        def preds(x):
            return [p for p in acts if p < graph.shape[0] and x < graph.shape[1] and graph[p, x] != 0]

        is_fork = set()
        for f in acts:
            for b in acts:
                if (f < graph.shape[0] and b < graph.shape[1] and graph[f, b] != 0
                        and b < len(post) and post[b] == post_and_value):
                    is_fork.add(f)
                    break

        def enclosing_fork(a):
            seen, queue = set(), [a]
            while queue:
                cur = queue.pop(0)
                if cur in seen:
                    continue
                seen.add(cur)
                ps = preds(cur)
                hit = [p for p in ps if p in is_fork]
                if hit:
                    return hit[0]
                queue.extend(ps)
            return -1

        for j in acts:
            inputs = [i for i in preds(j) if i < len(pre) and pre[i] == pre_and_value]
            if len(inputs) < 2:
                continue
            forks = {enclosing_fork(i) for i in inputs}
            if len(forks) > 1 or -1 in forks:
                name = lqn.hashnames[j] if j < len(lqn.hashnames) else str(j)
                raise RuntimeError(
                    "Activity '%s' joins branches of different AND forks; SolverLN supports "
                    "only properly nested (series-parallel) fork-join graphs." % name)

    def _is_srvn_ph(self) -> bool:
        """True when the layers are the collapsed phase-type ones of 'srvn.ph'."""
        return getattr(self, 'lnmethod', None) == 'srvn.ph'

    def _is_ph_encoding(self) -> bool:
        """True when the layers carry the COMPOSED phase-type server law rather
        than the routing encoding of the activity graph, under either layering.
        The encoding, not the layering, decides which update and reconstruction
        passes run, so every such dispatch asks this and not for one method name.
        """
        return getattr(self, 'lnmethod', None) in ('srvn.ph', 'flat.ph')

    def _probe_srvn_ph(self) -> bool:
        """Answer whether 'srvn.ph' can serve this model, without disturbing the solver.

        Both the feature gate and the series-parallel reduction can refuse, and
        the second only finds out by composing the per-entry workflows -- work the
        build then reuses, since those laws do not depend on the iterate.
        """
        try:
            self._ph_init_state()
            self._assert_srvn_ph_supported()
            self._ph_init_laws()
            self._ph_laws_ready = True
            return True
        except Exception as e:                                  # noqa: BLE001
            self._ph_laws_ready = False
            line_debug("LN: method=srvn cannot use srvn.ph on this model (%s)", e)
            return False

    def _build_layers(self):
        """Build layer submodels (matches MATLAB buildLayers)."""
        lqn = self.lqn

        # Method resolution. A method name carries both the LAYERING and the
        # ENCODING: 'srvn.ph' replaces the routing encoding of the activity graph
        # by a composed phase-type server law, 'srvn' is the alias that takes it
        # where it can serve the model and 'srvn.cs' otherwise, and 'flat.cs'
        # squashes every server into one submodel. The choice is made ONCE, here,
        # and every later dispatch reads self.lnmethod.
        # See _kb/06-solver-catalog.md (LN section).
        requested = ln_requested_method(getattr(self.options, 'method', None))
        # the method names the layering, so it sets it
        self._force_flat = requested in ('flat.cs', 'flat.ph')
        if requested == 'flat.ph':
            # the squashed layering with the composed law: ONE submodel holding
            # every server, and a caller visiting each of them once per
            # invocation. The feature gate is the srvn.ph one plus the refusals a
            # single submodel carries -- see _ph_flat_server_set.
            self._ph_laws_ready = False
            self._ph_init_state()
            self.lnmethod = 'flat.ph'
            self._build_layers_ph(flat=True)
            return
        if requested in ('srvn.ph', 'srvn'):
            hard = requested == 'srvn.ph'
            if self._is_flat_layering():
                if hard:
                    raise ValueError("method='srvn.ph' requires the srvn layering, because it "
                                     "replaces each server by a submodel of its own. Use "
                                     "method='srvn.cs' for that layering.")
                line_debug("LN: method=srvn cannot use srvn.ph under the flat layering")
            else:
                self._ph_laws_ready = False
                if hard:
                    self._ph_init_state()
                    self.lnmethod = 'srvn.ph'
                    self._build_layers_ph()
                    return
                if self._probe_srvn_ph():
                    self.lnmethod = 'srvn.ph'
                    self._build_layers_ph()
                    return
        # The label reports what was BUILT, so a model squashed through
        # options.config.layering reads back as 'flat.cs' even when no method
        # named it.
        if requested == 'moment3':
            self.lnmethod = 'moment3'
        else:
            self.lnmethod = 'flat.cs' if self._is_flat_layering() else 'srvn.cs'

        self._assert_series_parallel_forks()

        # Initialize ensemble with None for each potential layer
        self.ensemble = [None] * (lqn.nhosts + lqn.ntasks)

        # Initialize update maps as lists of lists
        servt_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]
        thinkt_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]
        actthinkt_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]
        arvproc_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]
        call_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]
        route_map = [[] for _ in range(lqn.nhosts + lqn.ntasks)]

        # see _kb/06-solver-catalog.md (LN section) for the layering taxonomy
        self._assert_call_groups()
        flat_servers = self._flat_server_set() if self._is_flat_layering() else []

        if flat_servers:
            flat_callers = [lqn.tshift + t for t in range(lqn.ntasks)
                            if not self.ignore[lqn.tshift + t]]
            self._build_layer_recursive(flat_servers, flat_callers, False,
                                        servt_map, thinkt_map, actthinkt_map,
                                        arvproc_map, call_map, route_map, flat=True)
        else:
            self._build_layers_srvn(servt_map, thinkt_map, actthinkt_map,
                                    arvproc_map, call_map, route_map)

        self._finish_layers(servt_map, thinkt_map, actthinkt_map,
                            arvproc_map, call_map, route_map, flat_servers)

    def _build_layers_srvn(self, servt_map, thinkt_map, actthinkt_map,
                           arvproc_map, call_map, route_map):
        """One submodel per processor and per called task (default layering)."""
        lqn = self.lqn
        # Build one submodel for every processor (host layer)
        for hidx in range(lqn.nhosts):
            if not self.ignore[hidx]:
                tasks_on_host = self._get_tasks_of_host(hidx)
                if tasks_on_host:
                    # skip a host layer with no callers, no nonzero demand, and no REF task (matches MATLAB's pure-delay-host skip, e.g. USAGE_DELAY).
                    has_callers = False
                    has_demand = False
                    has_ref_task = False
                    for tidx in tasks_on_host:
                        # Check if task is a reference task
                        if self._is_ref_task(tidx):
                            has_ref_task = True
                            break
                        # Check if task has callers
                        if self._get_callers_of_task(tidx):
                            has_callers = True
                            break
                        # forwarding targets need a host layer even without direct callers; MATLAB builds host layers unconditionally.
                        if self._is_forwarding_target_task(tidx):
                            has_callers = True
                            break
                        # a task needs its layer built if it has nonzero host demand via its activities (base Distribution wrappers expose only .mean).
                        activities = self._get_activities_of_task(tidx)
                        for aidx in activities:
                            if aidx < len(self.servtproc) and self.servtproc[aidx] is not None:
                                proc = self.servtproc[aidx]
                                if hasattr(proc, 'getMean'):
                                    dem_mean = proc.getMean()
                                elif hasattr(proc, 'mean'):
                                    dem_mean = proc.mean
                                else:
                                    dem_mean = 0.0
                                if dem_mean is not None and dem_mean > 0:
                                    has_demand = True
                                    break
                        if has_demand:
                            break

                    # Only build layer if tasks have callers, demand, or a REF task
                    if has_callers or has_demand or has_ref_task:
                        self._build_layer_recursive(hidx, tasks_on_host, True,
                                                   servt_map, thinkt_map, actthinkt_map,
                                                   arvproc_map, call_map, route_map)

        # Build one submodel for every task (task layer)
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if not self.ignore[tidx] and not self._is_ref_task(tidx):
                # Check if task has callers
                callers = self._get_callers_of_task(tidx)
                if callers:
                    self._build_layer_recursive(tidx, callers, False,
                                               servt_map, thinkt_map, actthinkt_map,
                                               arvproc_map, call_map, route_map)

    def _finish_layers(self, servt_map, thinkt_map, actthinkt_map,
                       arvproc_map, call_map, route_map, flat_servers):
        """Flatten the update maps and index the ensemble."""
        lqn = self.lqn
        # Convert maps to numpy arrays
        self.servt_classes_updmap = self._flatten_map(servt_map)
        self.thinkt_classes_updmap = self._flatten_map(thinkt_map)
        self.actthinkt_classes_updmap = self._flatten_map(actthinkt_map)
        self.arvproc_classes_updmap = self._flatten_map(arvproc_map)
        self.call_classes_updmap = self._flatten_map(call_map)
        self.route_prob_updmap = self._flatten_map(route_map)

        if self.route_prob_updmap is not None and len(self.route_prob_updmap) > 0:
            self.unique_route_prob_updmap = np.unique(self.route_prob_updmap[:, 0])
        else:
            self.unique_route_prob_updmap = np.array([])

        # Remove empty models and create idxhash
        empty_models = [i for i, e in enumerate(self.ensemble) if e is None]
        self.ensemble = [e for e in self.ensemble if e is not None]

        # Also compact solvers list to match ensemble
        self.solvers = [s for i, s in enumerate(self.solvers) if i not in empty_models and i < len(self.solvers)]
        # Extend solvers if needed to match ensemble length
        while len(self.solvers) < len(self.ensemble):
            self.solvers.append(None)

        # Position of each host/task element in the compacted ensemble. Element 0
        # is the first host, not the dead slot the 1-based space used to carry.
        self.idxhash = np.full(lqn.nhosts + lqn.ntasks, np.nan)
        layer_idx = 0
        for orig_idx in range(lqn.nhosts + lqn.ntasks):
            if orig_idx not in empty_models:
                self.idxhash[orig_idx] = layer_idx
                layer_idx += 1

        # Layers carrying an admission constraint need the region wait recovered in
        # update_metrics -- see _kb/06-solver-catalog.md (LN section)
        self.layer_has_region = [bool(getattr(e, 'regions', None)) for e in self.ensemble]
        self.layer_chains = [None] * len(self.ensemble)
        for e_idx, has_region in enumerate(self.layer_has_region):
            if has_region:
                # layer structure is iteration-invariant, so cache the chain matrix
                self.layer_chains[e_idx] = np.asarray(self.ensemble[e_idx].getStruct().chains)

        # Classify layers as host or task
        self.hostLayerIndices = []
        self.taskLayerIndices = []

        if flat_servers:
            # every server resolves to the single flat layer, which is at once
            # the host layer and the task layer
            self.idxhash = np.full(lqn.nhosts + lqn.ntasks, np.nan)
            for sidx in flat_servers:
                self.idxhash[sidx] = 0
            self.hostLayerIndices = [0]
            self.taskLayerIndices = [0]
            return

        for hidx in range(lqn.nhosts):
            if not np.isnan(self.idxhash[hidx]):
                self.hostLayerIndices.append(int(self.idxhash[hidx]))

        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if not np.isnan(self.idxhash[tidx]):
                self.taskLayerIndices.append(int(self.idxhash[tidx]))

    def _flatten_map(self, map_list: List[List]) -> np.ndarray:
        """Flatten a list of lists into a numpy array."""
        all_rows = []
        for rows in map_list:
            all_rows.extend(rows)
        if all_rows:
            return np.array(all_rows)
        return np.array([]).reshape(0, 4)

    def _get_tasks_of_host(self, hidx: int) -> List[int]:
        """Get task indices for a host processor."""
        lqn = self.lqn
        if isinstance(lqn.tasksof, dict):
            return lqn.tasksof.get(hidx, [])
        return []

    def _is_ref_task(self, tidx: int) -> bool:
        """Check if a task is a reference (REF) task."""
        lqn = self.lqn
        if hasattr(lqn, 'isref') and lqn.isref is not None:
            if isinstance(lqn.isref, dict):
                return lqn.isref.get(tidx, False)
            elif isinstance(lqn.isref, np.ndarray):
                # Handle 2D arrays - flatten and use direct index
                flat_isref = lqn.isref.flatten()
                if tidx < len(flat_isref):
                    return bool(flat_isref[tidx])
        if hasattr(lqn, 'sched') and lqn.sched is not None:
            ref_value = SchedStrategy.REF.value if hasattr(SchedStrategy.REF, 'value') else SchedStrategy.REF
            if isinstance(lqn.sched, dict):
                sched_val = lqn.sched.get(tidx, None)
                if sched_val is not None:
                    # Compare against both the enum and its value
                    return sched_val == SchedStrategy.REF or sched_val == ref_value
            elif isinstance(lqn.sched, np.ndarray):
                flat_sched = lqn.sched.flatten()
                if tidx < len(flat_sched):
                    sched_val = flat_sched[tidx]
                    return sched_val == SchedStrategy.REF or sched_val == ref_value
        return False

    def _get_callers_of_task(self, tidx: int) -> List[int]:
        """Get caller task indices for a task."""
        lqn = self.lqn
        callers = []

        # iscaller is a task-to-task matrix: iscaller[caller_task_idx, callee_task_idx]
        if hasattr(lqn, 'iscaller') and lqn.iscaller is not None:
            if isinstance(lqn.iscaller, np.ndarray):
                # Find all tasks that call this task (callers in column tidx)
                caller_indices = np.where(lqn.iscaller[:, tidx] > 0)[0]
                for caller_idx in caller_indices:
                    # Check if caller is a task (not processor/entry/activity)
                    if lqn.tshift <= caller_idx < lqn.tshift + lqn.ntasks:
                        if caller_idx not in callers:
                            callers.append(caller_idx)

        return callers

    def _extract_latest_metrics(self, aidx: int, layer_idx: int, nodeidx_0: int, classidx_0: int,
                                refstat_k: int = None, refclass_c: int = None):
        """Extract metrics from the latest result (helper for _update_metrics_default).

        Computes residt from QN/TN_ref instead of WN to avoid fork+loop visit distortion
        (matches MATLAB updateMetricsDefault.m).
        """
        if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
            result = self.results[-1][layer_idx]
            if result is not None and 'RN' in result:
                RN = result['RN']
                WN = result.get('WN', RN)
                TN = result['TN']
                QN = result.get('QN')

                if RN is not None and nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                    rn_val = RN[nodeidx_0, classidx_0]
                    tn_val = TN[nodeidx_0, classidx_0]

                    # Safeguard against extreme values from MVA numerical instability
                    # If RN becomes too large (> 1e10) or is NaN/Inf, keep previous value
                    max_servt = 1e10
                    if np.isfinite(rn_val) and rn_val <= max_servt and rn_val >= 0:
                        self.servt[aidx] = rn_val
                    # else: keep previous servt[aidx] value

                    # Compute residt from QN/TN_ref (matches MATLAB updateMetricsDefault.m)
                    if (refstat_k is not None and refclass_c is not None and
                            QN is not None and TN is not None):
                        TN_ref = TN[refstat_k, refclass_c] if (0 <= refstat_k < TN.shape[0] and 0 <= refclass_c < TN.shape[1]) else 0.0
                        if TN_ref > 1e-8:  # GlobalConstants.FineTol
                            qn_val = QN[nodeidx_0, classidx_0]
                            if np.isfinite(qn_val) and qn_val >= 0:
                                self.residt[aidx] = qn_val / TN_ref
                        else:
                            wn_val = WN[nodeidx_0, classidx_0] if WN is not None else rn_val
                            if np.isfinite(wn_val) and wn_val <= max_servt and wn_val >= 0:
                                self.residt[aidx] = wn_val
                    else:
                        # Fallback to WN if no refstat/refclass info
                        wn_val = WN[nodeidx_0, classidx_0] if WN is not None else rn_val
                        if np.isfinite(wn_val) and wn_val <= max_servt and wn_val >= 0:
                            self.residt[aidx] = wn_val

                    if np.isfinite(tn_val) and tn_val >= 0:
                        self.tput[aidx] = tn_val

    def _get_fork_fanout(self, aidx: int) -> int:
        """
        Get the fork fanout correction factor for an activity.

        For activities that are in a fork branch (after fork, before join),
        returns the number of parallel branches so throughput can be corrected.
        For fork sources, join targets, and non-fork activities, returns 1.

        Without Fork/Join nodes, probabilistic routing divides throughput by
        the number of branches. This function identifies activities that need
        correction (multiplication by fanout) to recover the correct throughput.

        Activities needing correction:
        1. Fork sources (B1) - routing normalization divides their throughput
        2. POST_AND activities (fork branch targets like B2, B3, B4)
        3. Activities in fork branch chains (successors of POST_AND before join)

        Activities NOT needing correction:
        - Join targets (B6) - receive sum from all branches
        """
        lqn = self.lqn

        if not hasattr(lqn, 'graph') or lqn.graph is None:
            return 1

        graph = lqn.graph
        if not isinstance(graph, np.ndarray):
            return 1

        post_and_value = 12  # ActivityPrecedenceType.ID_POST_AND
        pre_and_value = 2    # ActivityPrecedenceType.ID_PRE_AND

        def is_post_and_activity(act_idx: int) -> bool:
            """Check if an activity is POST_AND (fork branch target)."""
            if not hasattr(lqn, 'actposttype') or lqn.actposttype is None:
                return False
            actposttype = lqn.actposttype
            if not isinstance(actposttype, np.ndarray):
                return False
            flat_posttype = actposttype.flatten()
            if 0 < act_idx < len(flat_posttype):
                return flat_posttype[act_idx] == post_and_value
            return False

        def is_join_target(act_idx: int) -> bool:
            """Check if an activity is a join target (has PRE_AND predecessors)."""
            # A join target has predecessors that are PRE_AND (join sources)
            if not hasattr(lqn, 'actpretype') or lqn.actpretype is None:
                return False
            actpretype = lqn.actpretype
            if not isinstance(actpretype, np.ndarray):
                return False
            flat_pretype = actpretype.flatten()

            # Check if ANY predecessor of this activity is PRE_AND
            for pred_idx in range(graph.shape[0]):
                if pred_idx != act_idx and graph[pred_idx, act_idx] > 0:
                    if 0 < pred_idx < len(flat_pretype) and flat_pretype[pred_idx] == pre_and_value:
                        return True
            return False

        def count_post_and_successors(act_idx: int) -> int:
            """Count POST_AND successors of an activity."""
            count = 0
            for succ_idx in range(graph.shape[1]):
                if graph[act_idx, succ_idx] > 0 and is_post_and_activity(succ_idx):
                    count += 1
            return count

        def is_fork_source(act_idx: int) -> bool:
            """Check if an activity is a fork source (has POST_AND successors)."""
            return count_post_and_successors(act_idx) > 1 and not is_post_and_activity(act_idx)

        def get_fork_fanout_for_activity(act_idx: int, visited: set) -> int:
            """Recursively determine fork fanout for an activity."""
            if act_idx in visited:
                return 1
            visited.add(act_idx)

            # Join targets don't need correction - they receive from all branches
            if is_join_target(act_idx):
                return 1

            # fork sources need throughput correction (routing normalization dilutes it) but not visit correction (already correct).
            if is_fork_source(act_idx):
                return count_post_and_successors(act_idx)

            # Case 1: Activity is POST_AND (fork branch target)
            if is_post_and_activity(act_idx):
                # Find predecessor (fork source) and count its POST_AND successors
                for pred_idx in range(graph.shape[0]):
                    if pred_idx != act_idx and graph[pred_idx, act_idx] > 0:
                        fanout = count_post_and_successors(pred_idx)
                        if fanout > 1:
                            return fanout
                return 1

            # Case 2: Activity is in a fork branch chain (predecessor has fanout)
            for pred_idx in range(graph.shape[0]):
                if pred_idx != act_idx and graph[pred_idx, act_idx] > 0:
                    # Check if predecessor is POST_AND or has fanout
                    pred_fanout = get_fork_fanout_for_activity(pred_idx, visited)
                    if pred_fanout > 1:
                        return pred_fanout

            return 1

        return get_fork_fanout_for_activity(aidx, set())

    def _get_max_caller_fork_fanout(self, tidx: int) -> int:
        """
        Get the maximum fork fanout among all activities that call this task.

        For tasks called from AND-fork branches, this returns the fork fanout
        so that throughput can be corrected (multiplied back after probabilistic
        routing approximation divides it).

        Args:
            tidx: Task index

        Returns:
            Maximum fork fanout among callers (1 if no AND-fork callers)
        """
        lqn = self.lqn
        max_fanout = 1
        task_name = self._get_hashname(tidx) if tidx else str(tidx)

        # Get entries of this task
        entries = self._get_entries_of_task(tidx)
        if not entries:
            return 1

        # For each entry, find calling activities via callpair
        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return 1

        callpair = lqn.callpair
        if not isinstance(callpair, np.ndarray):
            return 1

        for eidx in entries:
            # Find calls targeting this entry (column 2 of callpair has target entry index)
            for cidx in range(callpair.shape[0]):
                if cidx < callpair.shape[0]:
                    tgt_eidx = int(callpair[cidx, 1]) if callpair.shape[1] > 1 else 0
                    if tgt_eidx == eidx:
                        # Found a call to this entry - get the source activity
                        src_aidx = self._get_call_source_activity(cidx)
                        if src_aidx is not None and src_aidx > 0:
                            # Get fork fanout of the calling activity
                            fanout = self._get_fork_fanout(src_aidx)
                            src_name = self._get_hashname(src_aidx) if src_aidx else str(src_aidx)
                            max_fanout = max(max_fanout, fanout)

        return max_fanout

    # NOTE: Visit correction for fork-join activities is NOT implemented because MVA
    # recomputes visits internally from the routing matrix, ignoring any manual
    # modifications to nodevisits. Instead, throughput correction is applied in
    # get_ensemble_avg by multiplying raw MVA throughputs by the fork fanout.
    #
    # For exact MATLAB parity, Fork/Join/Router nodes would need to be added to the
    # layer models (as MATLAB does in buildLayersRecursive.m), but this is a
    # significant undertaking. The current throughput correction provides reasonable
    # approximations for most fork-join networks.

    def _build_layer_recursive(self, idx_set, callers: List[int], is_host_layer: bool,
                               servt_map, thinkt_map, actthinkt_map,
                               arvproc_map, call_map, route_map, flat: bool = False):
        """Build a layer submodel (matches MATLAB buildLayersRecursive).

        IDX_SET is the server element of this layer, a scalar under 'srvn'
        layering and the whole host+task set under 'flat' layering.
        """
        lqn = self.lqn
        if isinstance(idx_set, (list, tuple, np.ndarray)):
            idx_set = [int(v) for v in idx_set]
        else:
            idx_set = [int(idx_set)]
        idx = idx_set[0]  # layer key: model name, ensemble slot and update-map column

        # Create Network for this layer
        model_name = self._get_hashname(idx)
        layer_model = Network(model_name + '.Flat' if flat else model_name)
        if hasattr(layer_model, 'set_checks'):
            layer_model.set_checks(False)

        # Create attribute storage
        layer_model.attribute = OptionsDict({
            'hosts': [],
            'tasks': [],
            'entries': [],
            'activities': [],
            'calls': [],
            'clientIdx': None,
            'serverIdx': None,
            'sourceIdx': None,
            'cacheIdx': None,  # Cache node index if cache layer
            'iscachelayer': False,  # Flag for cache layer
        })

        # Detect cache layer (MATLAB buildLayersRecursive line 36)
        # iscachelayer = all(lqn.iscache(callers)) && ishostlayer
        iscachelayer = False
        if not flat and is_host_layer and hasattr(lqn, 'iscache') and lqn.iscache is not None:
            iscache_arr = lqn.iscache.flatten() if isinstance(lqn.iscache, np.ndarray) else lqn.iscache
            # Check if ALL callers are cache tasks
            if len(callers) > 0:
                iscachelayer = True
                for caller_idx in callers:
                    if caller_idx < len(iscache_arr):
                        if not iscache_arr[caller_idx]:
                            iscachelayer = False
                            break
                    else:
                        iscachelayer = False
                        break

        layer_model.attribute['iscachelayer'] = iscachelayer

        # Get number of servers
        nservers = self._get_nservers(idx)
        sched = self._get_sched(idx)

        # fan-out replication: single replica when caller fan-out covers task replicas; see _kb/06-solver-catalog.md LN Fan-out single-replica modeling.
        raw_replicas = 1 if flat else int(self._get_repl(idx))
        reduce_fanout = False
        if raw_replicas > 1 and len(callers) > 0:
            if not is_host_layer and hasattr(lqn, 'fanout') and lqn.fanout is not None:
                reduce_fanout = True
                for c in callers:
                    fo = lqn.fanout[c, idx] if c < lqn.fanout.shape[0] and idx < lqn.fanout.shape[1] else 0
                    if fo < raw_replicas:
                        reduce_fanout = False
                        break
            elif is_host_layer:
                reduce_fanout = True
                for c in callers:
                    if int(self._get_repl(c)) != raw_replicas:
                        reduce_fanout = False
                        break

        if reduce_fanout:
            nreplicas = 1
            if not is_host_layer:
                self.single_replica_tasks.append(idx)
        else:
            nreplicas = raw_replicas

        # Create stations
        has_sync_callers = self._has_sync_callers(idx, callers)

        if flat or is_host_layer or has_sync_callers:
            # Create client delay node
            client_delay = Delay(layer_model, 'Clients')
            layer_model.attribute['clientIdx'] = 1
            layer_model.attribute['serverIdx'] = 2
        else:
            layer_model.attribute['serverIdx'] = 1
            layer_model.attribute['clientIdx'] = None

        # One station (times its replicas) per server element of the layer
        srv_stations = {}
        server_idx_of = {}
        host_stations = []
        task_stations = []
        for sidx in idx_set:
            s_is_host = sidx <= lqn.nhosts
            s_name = self._get_hashname(sidx)
            if flat and any(str(n.getName()) == s_name for n in layer_model.get_nodes()):
                # an LQN processor and the task it hosts may share a name; two
                # stations of one layer must not, or link() gives their
                # class-switch nodes the same name and merges their arcs
                s_name = s_name + ('.host' if s_is_host else '.task')
            s_nservers = self._get_nservers(sidx)
            s_sched = self._get_sched(sidx)
            stations_of = []
            for m in range(1, nreplicas + 1):
                if m == 1:
                    ss = Queue(layer_model, s_name, s_sched)
                else:
                    ss = Queue(layer_model, s_name + '.' + str(m), s_sched)
                ss.set_number_of_servers(s_nservers)
                ss.attribute = OptionsDict({
                    'ishost': s_is_host,
                    'idx': sidx
                })
                # successive same-host activities retain the server; mark immediate feedback so simulators do not re-queue behind waiting jobs.
                ss.set_immediate_feedback(True)
                stations_of.append(ss)
            srv_stations[sidx] = stations_of
            server_idx_of[sidx] = len(layer_model.get_nodes()) - nreplicas + 1
            if s_is_host:
                host_stations.append(server_idx_of[sidx])
            else:
                task_stations.append(server_idx_of[sidx])

        layer_model.attribute['srv_stations'] = srv_stations
        layer_model.attribute['serverIdxOf'] = server_idx_of
        layer_model.attribute['hostStations'] = host_stations
        layer_model.attribute['taskStations'] = task_stations
        layer_model.attribute['flat'] = flat

        server_stations = srv_stations[idx]
        server_station = server_stations[0]  # the layer's own server, sole server under 'srvn'
        layer_model.attribute['nreplicas'] = nreplicas
        layer_model.attribute['server_stations'] = server_stations

        # Source/Sink lazily created and reused for all open classes in this layer; mirrors JAR SolverLN.java:719-726.
        source_station = None
        sink_station = None
        if hasattr(lqn, 'arrival') and lqn.arrival:
            for tidx_caller in callers:
                # an arrival that is the only way into the task is carried by the caller
                # chain instead, not by a stream -- see _open_arrival_rate_of
                if self._is_open_arrival_only(tidx_caller):
                    continue
                for eidx in self._get_entries_of_task(tidx_caller):
                    if eidx in lqn.arrival and lqn.arrival[eidx] is not None:
                        source_station = Source(layer_model, 'Source')
                        sink_station = Sink(layer_model, 'Sink')
                        break
                if source_station is not None:
                    break
        layer_model.attribute['source_station'] = source_station
        layer_model.attribute['sink_station'] = sink_station
        layer_model.attribute['entry_open_classes'] = []  # list of (OpenClass, eidx)
        layer_model.attribute['async_open_classes'] = []  # list of (OpenClass, cidx, callmean)

        # Detect POST_AND / PRE_AND activities for Fork/Join routing
        # (MATLAB buildLayersRecursive.m lines 42-66)
        post_and_value = 12  # ActivityPrecedenceType.ID_POST_AND
        pre_and_value = 2    # ActivityPrecedenceType.ID_PRE_AND
        is_post_and_act = set()
        is_pre_and_act = set()

        acts_in_caller = []
        for tidx_caller in callers:
            acts_in_caller.extend(self._get_activities_of_task(tidx_caller))

        if hasattr(lqn, 'actposttype') and lqn.actposttype is not None:
            flat_posttype = lqn.actposttype.flatten()
            flat_pretype = lqn.actpretype.flatten() if hasattr(lqn, 'actpretype') and lqn.actpretype is not None else np.array([])
            for aidx in acts_in_caller:
                if 0 < aidx < len(flat_posttype):
                    if flat_posttype[aidx] == post_and_value:
                        is_post_and_act.add(aidx)
                if 0 < aidx < len(flat_pretype):
                    if flat_pretype[aidx] == pre_and_value:
                        is_pre_and_act.add(aidx)

        has_fork = any(aidx in is_post_and_act for aidx in acts_in_caller)

        maxfanout = 1
        graph = lqn.graph if hasattr(lqn, 'graph') and isinstance(lqn.graph, np.ndarray) else None
        if graph is not None:
            for aidx in acts_in_caller:
                if aidx < graph.shape[0]:
                    successors = [j for j in range(graph.shape[1]) if graph[aidx, j] != 0]
                    post_and_count = sum(1 for s in successors if s in is_post_and_act)
                    if post_and_count > 0:
                        maxfanout = max(maxfanout, post_and_count)

        fork_node = None
        fork_output_routers = {}
        if has_fork:
            fork_node = Fork(layer_model, 'Fork_PostAnd')
            for f in range(1, maxfanout + 1):
                fork_output_routers[f] = Router(layer_model, f'Fork_PostAnd_{f}')

        has_join = any(aidx in is_pre_and_act for aidx in acts_in_caller)
        join_node = None
        if has_join:
            join_node = Join(layer_model, 'Join_PreAnd', fork_node)

        # Store Fork/Join info in layer attributes for throughput correction
        layer_model.attribute['fork_node'] = fork_node
        layer_model.attribute['fork_output_routers'] = fork_output_routers
        layer_model.attribute['join_node'] = join_node
        layer_model.attribute['is_post_and_act'] = is_post_and_act
        layer_model.attribute['is_pre_and_act'] = is_pre_and_act
        layer_model.attribute['maxfanout'] = maxfanout
        layer_model.attribute['has_fork'] = has_fork

        # The fork-join transform mints its own Source/Sink pair, detaching the open
        # stream already routed through this one: see _kb/06-solver-catalog.md (LN section)
        if has_fork and source_station is not None:
            raise ValueError(f"SolverLN: layer '{layer_model.getName()}' carries both an AND fork "
                             "and an open stream (an async call or an entry arrival); the "
                             "fork-join transform needs a Source of its own")

        # Create Cache node for cache layers (MATLAB buildLayersRecursive.m lines 36-39)
        cache_node = None
        if iscachelayer and len(callers) > 0:
            # Get cache parameters from the first cache task caller
            cache_task_idx = callers[0]
            if hasattr(lqn, 'nitems') and lqn.nitems is not None:
                nitems = int(lqn.nitems[cache_task_idx, 0]) if cache_task_idx < lqn.nitems.shape[0] else 0
                if nitems > 0:
                    # Get item capacity
                    itemcap = lqn.itemcap.get(cache_task_idx, np.array([1])) if hasattr(lqn, 'itemcap') and lqn.itemcap else np.array([1])
                    # Get replacement strategy
                    replacestrat_val = int(lqn.replacestrat[cache_task_idx, 0]) if hasattr(lqn, 'replacestrat') and lqn.replacestrat is not None else 0
                    # Convert to ReplacementStrategy enum
                    try:
                        replacestrat = ReplacementStrategy(replacestrat_val)
                    except (ValueError, KeyError):
                        replacestrat = ReplacementStrategy.RR  # Default to Random Replacement
                    # Get cache name from hashnames
                    cache_name = lqn.hashnames[cache_task_idx] if hasattr(lqn, 'hashnames') and cache_task_idx < len(lqn.hashnames) else f'Cache_{cache_task_idx}'
                    # Create Cache node
                    cache_node = Cache(layer_model, cache_name, nitems, itemcap, replacestrat)
                    layer_model.attribute['cacheNode'] = cache_node
                    # Update cacheIdx - Cache is added after server, so its index is serverIdx + 1
                    # Cache is the last node added, so use len(get_nodes()) after it was added
                    layer_model.attribute['cacheIdx'] = len(layer_model.get_nodes())

        # Store server attributes. Under flat layering the station indices live
        # in hostStations/taskStations only: attribute['hosts'] / ['tasks'] rows
        # are [class index, LQN element] pairs that consumers match on column 2.
        if not flat:
            if is_host_layer:
                layer_model.attribute['hosts'].append([None, layer_model.attribute['serverIdx']])
            else:
                layer_model.attribute['tasks'].append([None, layer_model.attribute['serverIdx']])

        # Create classes and set up routing
        self._create_classes_and_routing(layer_model, idx_set, callers, is_host_layer,
                                        servt_map, thinkt_map, actthinkt_map,
                                        arvproc_map, call_map, route_map,
                                        reduce_fanout=reduce_fanout, flat=flat)

        # fork-join visit correction not applied here: MVA recomputes visits internally; throughput correction in get_ensemble_avg handles FJ semantics.

        # Store the layer model
        self.ensemble[idx] = layer_model

        # A setup no longer changes how a layer is solved: the cold start is
        # charged to the entry by _setup_charge, not wired into the station, so
        # the layer is an ordinary one and the user's own solver serves it.
        solver = self.solver_factory(layer_model)

        self._assert_layer_solver_supports_model(solver, layer_model, idx)
        self._detach_layer_config(solver)
        self._silence_layer_solver(solver)

        if idx < len(self.solvers):
            self.solvers[idx] = solver
        else:
            while len(self.solvers) <= idx:
                self.solvers.append(None)
            self.solvers[idx] = solver

    def _is_flat_layering(self) -> bool:
        """True when the method or options.config.layering asks for the single flat layer."""
        # method='flat'/'flat.cs' names the layering, so it wins here
        if getattr(self, '_force_flat', False):
            return True
        cfg = getattr(self.options, 'config', None)
        if cfg is None:
            return False
        try:
            lay = cfg['layering']
        except (KeyError, TypeError):
            lay = getattr(cfg, 'layering', None)
        return isinstance(lay, str) and lay.lower() in ('flat', 'squashed')

    def _assert_call_groups(self) -> None:
        """Reject routed call groups under any layering that cannot carry them.

        A group states the order in which one caller visits several callees. The
        srvn layering puts every callee in a submodel of its own and replaces it,
        in the caller's submodel, by a surrogate delay, so the callees are never
        co-resident and no node has arcs to more than one of them: the order has
        nowhere to be expressed and would be silently degraded to the aggregate
        call means. The squashed layering keeps all of them as stations of one
        model, which is what makes the strategy representable.
        """
        groups = getattr(self.lqn, 'callgroups', None)
        if not groups:
            return
        if not self._is_flat_layering():
            strategies = sorted({str(getattr(s, 'name', s)) for _, s, _ in groups})
            raise ValueError(
                "Call groups routed by %s require the squashed layering; set "
                "options.config['layering']='flat'. Under srvn the targets never "
                "share a submodel, so the dispatch order cannot be represented."
                % ', '.join(strategies))
        # Under flat the group becomes one dispatch hop with n destinations (see
        # the routing walk), so the strategy is representable. It is only honoured
        # by a layer solver that implements state-dependent routing, though: MVA,
        # NC and FLD would silently return the probabilistic split instead.
        if not self._layer_solver_supports_routed_groups():
            raise ValueError(
                'Routed call groups need a layer solver with state-dependent '
                'routing (CTMC or SSA); MVA, NC and FLD would silently return '
                'the probabilistic split under a round-robin or JSQ label.')

    def probe_layer_solver_name(self):
        """Name the layer solver this ensemble runs, WITHOUT building the layers.

        A delegated solve (lang='java'/'cpp') never enters iterate(), so
        `self.solvers` is still empty when the dispatcher has to name the layer
        engine, and reading it off that list reports the NATIVE DEFAULT (MVA)
        however the caller built the solver. A lambda factory hides the class
        from `_layer_solver_cls` too, so `LN(model, lambda m: NC(m, opts))` was
        delegated as an MVA-layered ensemble -- a different fixed point, not a
        different spelling of the same one (lcq_threehosts: cache hit 0.5 under
        MVA layers against 0.48331 under NC ones).

        THE FACTORY IS THE DECLARATION, so it is applied to a layer and the
        product named. CTMC and MAM are skipped because they are the automatic
        per-layer substitutions (a finite capacity region, a SetupTask's setup
        times), not a choice the caller made; if every layer resolves to one of
        those the answer is None and the caller keeps its own default.
        """
        cls = getattr(self, '_layer_solver_cls', None)
        if cls is not None:
            return getattr(cls, '__name__', str(cls))
        factory = getattr(self, 'solver_factory', None)
        ensemble = getattr(self, 'ensemble', None) or []
        if factory is None:
            return None
        for layer in ensemble:
            if layer is None:
                continue
            try:
                probe = factory(layer)
            except Exception:
                continue
            if probe is None:
                continue
            name = None
            getname = getattr(probe, 'getName', None)
            if getname is not None:
                try:
                    name = getname()
                except Exception:
                    name = None
            if not name:
                name = type(probe).__name__
            bare = str(name).replace('Solver', '').upper()
            if bare in ('CTMC', 'MAM'):
                continue
            return name
        return None

    def _layer_solver_supports_routed_groups(self) -> bool:
        """True when the layer solver factory declares state-dependent routing."""
        cls = getattr(self, '_layer_solver_cls', None) or getattr(self, 'solver_factory', None)
        name = getattr(cls, '__name__', '') or type(cls).__name__ if cls is not None else ''
        return any(tag in name for tag in ('CTMC', 'SSA', 'LDES', 'JMT'))

    def _call_groups_by_cidx(self):
        """Resolve lqn.callgroups from target entries to call indices.

        Returns (by_cidx, members): by_cidx maps a call index to (gid, strategy),
        members maps gid to the ordered list of that group's call indices. A group
        that does not resolve to at least two calls is dropped, so a stale group
        cannot silently rewrite a single call's routing.
        """
        if getattr(self, '_callgroup_cache', None) is not None:
            return self._callgroup_cache
        groups = getattr(self.lqn, 'callgroups', None) or []
        by_cidx, members = {}, {}
        callsof = self.lqn.callsof if isinstance(self.lqn.callsof, dict) else {}
        for gid, (aidx, strategy, entry_idxs) in enumerate(groups):
            want = {int(e) for e in entry_idxs}
            found = []
            for cidx in callsof.get(aidx, []):
                tgt = self._get_call_target_entry(cidx)
                if tgt is not None and int(tgt) in want:
                    found.append(cidx)
            if len(found) >= 2:
                members[gid] = found
                for cidx in found:
                    by_cidx[cidx] = (gid, strategy)
        self._callgroup_cache = (by_cidx, members)
        return self._callgroup_cache

    def _flat_server_set(self) -> List[int]:
        """Processors and called tasks that become stations of the flat layer.

        The features a single submodel cannot carry are rejected here rather
        than silently dropped.
        """
        lqn = self.lqn
        nelem = lqn.nhosts + lqn.ntasks
        for idx in range(nelem):
            if float(self._get_repl(idx)) > 1:
                raise ValueError('Flat layering does not support replicated processors or '
                                 'tasks, use the default srvn layering.')
        if getattr(lqn, 'iscache', None) is not None:
            arr = lqn.iscache.flatten()
            if any(bool(arr[i]) for i in range(min(nelem, len(arr)))):
                raise ValueError('Flat layering does not support cache tasks, use the '
                                 'default srvn layering.')
        if getattr(lqn, 'hassetup', None) is not None:
            arr = np.asarray(lqn.hassetup).flatten()
            if any(bool(arr[i]) for i in range(min(nelem, len(arr)))):
                raise ValueError('Flat layering does not support setup tasks, use the '
                                 'default srvn layering.')

        servers = []
        for hidx in range(lqn.nhosts):
            if not self.ignore[hidx] and self._get_tasks_of_host(hidx):
                servers.append(hidx)
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx] or self._is_ref_task(tidx):
                continue
            if self._get_callers_of_task(tidx):
                servers.append(tidx)
        if not servers:
            raise ValueError('Flat layering found no server: the model has no processor '
                             'with tasks.')
        return servers

    def _servers_for(self, layer_model, elem_idx) -> List:
        """Stations of ELEM_IDX when it is a server of this layer, empty otherwise."""
        if elem_idx is None:
            return []
        srv = layer_model.attribute.get('srv_stations') if hasattr(layer_model, 'attribute') else None
        if not srv:
            return []
        return srv.get(int(elem_idx), [])

    def _host_is_server(self, layer_model, tidx) -> bool:
        """True when the processor of task TIDX is a server of this layer."""
        return bool(self._servers_for(layer_model, self._get_parent(tidx)))

    def _station_idx_of(self, layer_model, elem_idx):
        """Station index of ELEM_IDX inside LAYER_MODEL, falling back to the
        layer's own server when ELEM_IDX is not a server there."""
        default = layer_model.attribute.get('serverIdx', 1) if hasattr(layer_model, 'attribute') else 1
        if elem_idx is None:
            return default
        table = layer_model.attribute.get('serverIdxOf') if hasattr(layer_model, 'attribute') else None
        if table and int(elem_idx) in table:
            return table[int(elem_idx)]
        return default

    def _station_idx_of_class(self, layer_model, cls):
        """Station of LAYER_MODEL serving CLS: the processor of an activity, the
        called task of a call, the layer's own server otherwise."""
        elem = None
        attr = getattr(cls, 'attribute', None)
        if attr is not None and len(attr) > 1:
            if attr[0] == LayeredNetworkElement.ACTIVITY:
                elem = self._get_parent(self._get_parent(attr[1]))
            elif attr[0] == LayeredNetworkElement.CALL:
                cidx = int(attr[1])
                if self.lqn.callpair is not None and cidx < len(self.lqn.callpair):
                    elem = self._get_parent(int(self.lqn.callpair[cidx, 1]))
        return self._station_idx_of(layer_model, elem)

    def _get_hashname(self, idx: int) -> str:
        """Get the hash name for an LQN element."""
        lqn = self.lqn
        if hasattr(lqn, 'hashnames') and lqn.hashnames is not None:
            if isinstance(lqn.hashnames, dict):
                return lqn.hashnames.get(idx, f'Node_{idx}')
            elif isinstance(lqn.hashnames, (list, np.ndarray)):
                # hashnames is 0-based and contiguous: index i is element i
                if idx < len(lqn.hashnames):
                    return lqn.hashnames[idx]
        return f'Node_{idx}'

    def _get_nservers(self, idx: int):
        """Get number of servers for an element.

        Matches MATLAB's use of maxmult in buildLayersRecursive line 7-8, 31:
          mult = lqn.maxmult; % this removes spare capacity that cannot be used
          serverStation{m}.setNumberOfServers(mult(idx))

        Uses maxmult instead of mult because maxmult "removes spare capacity
        that cannot be used" (MATLAB comment). For processors with INF mult
        (delay nodes), maxmult = 0 which means all capacity can be used.
        """
        lqn = self.lqn
        val = 1

        # Use maxmult if available (MATLAB line 7: mult = lqn.maxmult)
        if hasattr(lqn, 'maxmult') and lqn.maxmult is not None:
            if isinstance(lqn.maxmult, dict):
                val = lqn.maxmult.get(idx, 1)
            elif isinstance(lqn.maxmult, np.ndarray):
                flat_maxmult = lqn.maxmult.flatten()
                if idx < len(flat_maxmult):
                    val = flat_maxmult[idx]
                    if isinstance(val, np.ndarray):
                        val = val.item() if val.size == 1 else 1
        elif hasattr(lqn, 'mult') and lqn.mult is not None:
            # Fallback to mult if maxmult not available
            if isinstance(lqn.mult, dict):
                val = lqn.mult.get(idx, 1)
            elif isinstance(lqn.mult, np.ndarray):
                flat_mult = lqn.mult.flatten()
                if idx < len(flat_mult):
                    val = flat_mult[idx]
                    if isinstance(val, np.ndarray):
                        val = val.item() if val.size == 1 else 1

        if not isinstance(val, (int, float)):
            return 1

        # Return np.inf for infinite servers (INF scheduling)
        if not np.isfinite(val):
            return np.inf

        # maxmult=0 means infinite servers (MATLAB: if isinf(mult), maxmult=0)
        if val == 0:
            return np.inf

        # INF scheduling always returns infinite servers regardless of maxmult (MATLAB warns but uses infinite servers).
        sched = self._get_sched(idx)
        if sched == SchedStrategy.INF:
            return np.inf

        return max(1, int(val))

    def _get_sched(self, idx: int) -> SchedStrategy:
        """Get scheduling strategy for an element."""
        lqn = self.lqn
        sched_val = None
        if hasattr(lqn, 'sched') and lqn.sched is not None:
            if isinstance(lqn.sched, dict):
                sched_val = lqn.sched.get(idx, None)
            elif isinstance(lqn.sched, np.ndarray):
                # Handle 2D arrays (shape like (1, n))
                flat_sched = lqn.sched.flatten()
                if idx < len(flat_sched):
                    sched_val = int(flat_sched[idx])

        if sched_val is None:
            return SchedStrategy.PS

        # Convert integer value to SchedStrategy enum
        if isinstance(sched_val, int):
            # Find the SchedStrategy with this value
            for strategy in SchedStrategy:
                if strategy.value == sched_val:
                    return strategy
            # Fallback to PS if value not found
            return SchedStrategy.PS
        elif isinstance(sched_val, SchedStrategy):
            return sched_val
        else:
            return SchedStrategy.PS

    def _has_sync_callers(self, idx: int, callers: List[int]) -> bool:
        """Check if any callers make synchronous calls to this element."""
        lqn = self.lqn

        # Get entries of this element
        entries = []
        if isinstance(lqn.entriesof, dict):
            entries = lqn.entriesof.get(idx, [])

        if not entries:
            return False

        # Check for sync callers
        if hasattr(lqn, 'issynccaller') and lqn.issynccaller is not None:
            for tidx in callers:
                for eidx in entries:
                    if isinstance(lqn.issynccaller, np.ndarray):
                        if lqn.issynccaller[tidx - 1, eidx - 1] > 0:
                            return True

        return True  # Default to true for safety

    def _has_direct_callers_for_caller(self, tidx_caller: int) -> bool:
        """True if caller task is REF, has sync/async callers, or its entries have open arrivals.

        Mirrors MATLAB buildLayersRecursive.m:130-154 and JAR SolverLN.java:634-655.
        """
        lqn = self.lqn
        if self._is_ref_task(tidx_caller):
            return True
        entries = self._get_entries_of_task(tidx_caller)
        for eidx in entries:
            if hasattr(lqn, 'issynccaller') and lqn.issynccaller is not None:
                col = np.asarray(lqn.issynccaller)
                if col.ndim == 2 and eidx < col.shape[1]:
                    if np.any(col[:, eidx] != 0):
                        return True
            if hasattr(lqn, 'isasynccaller') and lqn.isasynccaller is not None:
                col = np.asarray(lqn.isasynccaller)
                if col.ndim == 2 and eidx < col.shape[1]:
                    if np.any(col[:, eidx] != 0):
                        return True
            if hasattr(lqn, 'arrival') and lqn.arrival \
                    and eidx in lqn.arrival and lqn.arrival[eidx] is not None:
                return True
        return False

    def _open_arrival_rate_of(self, tidx: int) -> float:
        """Total exogenous rate into the entries of TIDX, zero unless the arrival is the
        only way in.

        A task reached only by an entry arrival has no task layer, because no task calls
        it, so _update_think_times never gives its caller class a surrogate delay and the
        class cycles against an Immediate one. Adding an open stream on top of that
        unthrottled chain saturated lqn_open_arrival: the processor at 0.68 against 0.32
        from lqns, lqsim and LDES alike. The chain is the representation that honours the
        thread pool, so the layer builder drops the stream for these tasks and the chain
        is closed on this rate instead, exactly as a forwarding target is. With a caller
        or a forwarding source the stream rides a class of its own and this returns 0.
        """
        lqn = self.lqn
        if self._is_ref_task(tidx) or not getattr(lqn, 'arrival', None):
            return 0.0
        entries = self._get_entries_of_task(tidx)
        if not entries:
            return 0.0
        for eidx in entries:
            for name in ('issynccaller', 'isasynccaller'):
                mat = getattr(lqn, name, None)
                if mat is None:
                    continue
                col = np.asarray(mat)
                if col.ndim == 2 and eidx < col.shape[1] and np.any(col[:, eidx] != 0):
                    return 0.0
        if self._is_forwarding_target_task(tidx):
            return 0.0
        rate = 0.0
        for eidx in entries:
            arv = lqn.arrival.get(eidx)
            if arv is None:
                continue
            m = arv.getMean()
            if np.isfinite(m) and m > GlobalConstants.FineTol:
                rate += 1.0 / m
        return rate

    def _is_open_arrival_only(self, tidx: int) -> bool:
        """True when an entry arrival is the only way requests reach task TIDX."""
        return self._open_arrival_rate_of(tidx) > GlobalConstants.FineTol

    def _is_forwarding_target_task(self, tidx: int) -> bool:
        """True if any entry of this task is the target of a forwarding call.

        Mirrors MATLAB buildLayersRecursive.m isForwardingTarget (lines 150-156):
        forwarding targets get host-layer classes even without direct callers.
        """
        lqn = self.lqn
        if not hasattr(lqn, 'calltype') or lqn.calltype is None:
            return False
        entries = self._get_entries_of_task(tidx)
        if not entries:
            return False
        for cidx in range(lqn.ncalls):
            if cidx < len(lqn.calltype) and int(lqn.calltype[cidx]) == CallType.FWD:
                if int(lqn.callpair[cidx, 1]) in entries:
                    return True
        return False

    def _is_sync_caller_to_entries_of(self, tidx_caller: int, idx: int) -> bool:
        """True if caller syncs to any entry under element idx. MATLAB line 156 right-side,
        JAR SolverLN.java:624-630.
        """
        lqn = self.lqn
        entries = self._get_entries_of_task(idx)
        if not entries:
            return False
        if not hasattr(lqn, 'issynccaller') or lqn.issynccaller is None:
            return False
        mat = np.asarray(lqn.issynccaller)
        if mat.ndim != 2:
            return False
        for eidx in entries:
            if tidx_caller < mat.shape[0] and eidx < mat.shape[1]:
                if mat[tidx_caller, eidx] != 0:
                    return True
        return False

    def _create_classes_and_routing(self, layer_model: Network, idx_set,
                                    callers: List[int], is_host_layer: bool,
                                    servt_map, thinkt_map, actthinkt_map,
                                    arvproc_map, call_map, route_map,
                                    reduce_fanout: bool = False, flat: bool = False):
        """
        Create classes and routing for a layer (simplified version).

        This is a simplified implementation. For 100% parity, the full
        MATLAB buildLayersRecursive logic would need to be ported.
        """
        lqn = self.lqn

        # Initialize SetupTask/MAM flags (these are set in _build_layer but
        # also referenced here for DelayOff handling)
        use_mam_solver = False
        function_task_idx = None
        if is_host_layer and hasattr(lqn, 'hassetup') and lqn.hassetup is not None:
            all_callers_function = True
            for caller_idx in callers:
                if caller_idx < lqn.hassetup.shape[1]:
                    if lqn.hassetup[0, caller_idx] != 1:
                        all_callers_function = False
                        break
                else:
                    all_callers_function = False
                    break
            if all_callers_function and len(callers) > 0:
                function_task_idx = callers[0]
                if hasattr(lqn, 'setuptime') and lqn.setuptime is not None:
                    if isinstance(lqn.setuptime, dict) and function_task_idx in lqn.setuptime and lqn.setuptime[function_task_idx] is not None:
                        use_mam_solver = True
                    elif isinstance(lqn.setuptime, np.ndarray):
                        flat_setuptime = lqn.setuptime.flatten()
                        if function_task_idx < len(flat_setuptime) and flat_setuptime[function_task_idx] is not None:
                            use_mam_solver = True

        # Get stations
        if isinstance(idx_set, (list, tuple, np.ndarray)):
            idx_set = [int(v) for v in idx_set]
        else:
            idx_set = [int(idx_set)]
        idx = idx_set[0]  # layer key: update-map column and ensemble slot
        srv_stations = layer_model.attribute.get('srv_stations', {})
        stations = layer_model.get_nodes()
        client_delay = None
        for s in stations:
            if isinstance(s, Delay):
                client_delay = s
                break
        server_station = srv_stations.get(idx, [None])[0]

        def servers_for(elem_idx):
            # Stations of ELEM_IDX when it is a server of this layer, [] otherwise
            if elem_idx is None:
                return []
            return srv_stations.get(int(elem_idx), [])

        def all_server_stations():
            out = []
            for _sid in idx_set:
                out.extend(srv_stations.get(_sid, []))
            return out

        def host_is_server(tidx_):
            return bool(servers_for(self._get_parent(tidx_)))

        if server_station is None:
            return

        # Create classes for each caller
        for tidx_caller in callers:
            # phase-2a gating mirrors MATLAB buildLayersRecursive.m:156/JAR SolverLN.java:657; activity/call loops unconditional, async-only callers reach ASYNC.
            caller_on_server = host_is_server(tidx_caller)
            has_direct_callers = self._has_direct_callers_for_caller(tidx_caller) if caller_on_server else False
            is_fwd_target = self._is_forwarding_target_task(tidx_caller) if caller_on_server else False
            # the TASK members of the set, and 0-BASED: hosts occupy 0..nhosts-1
            # and tasks nhosts..nhosts+ntasks-1, so the first task sits AT
            # nhosts. MATLAB's `sidx > lqn.nhosts` is the 1-based form of this
            # test and reads across unchanged only there. Under `>` the first
            # task's own layer never sees its callers as layer clients, so it is
            # built with no closed class and no population at all: on lqn_ofbiz
            # that is FrontEnd_CPU_Task, whose processor layer then ran 90 jobs
            # against a zero think time and saturated (Util 1.000 against 0.116).
            is_sync_caller_to_entries = any(self._is_sync_caller_to_entries_of(tidx_caller, _sid)
                                            for _sid in idx_set if _sid >= lqn.nhosts)
            create_caller_class = (caller_on_server and (has_direct_callers or is_fwd_target)) or is_sync_caller_to_entries

            if client_delay is None:
                continue

            if create_caller_class:
                # job population: single-replica callers use mult(caller), else mult(caller)*repl(caller); mirrors MATLAB buildLayersRecursive.m:162-168.
                mult = self._get_mult(tidx_caller)
                repl = self._get_repl(tidx_caller)
                caller_is_single_replica = reduce_fanout or (tidx_caller in self.single_replica_tasks)
                if caller_is_single_replica:
                    njobs = mult
                else:
                    njobs = mult * repl
                if np.isinf(njobs):
                    # If caller is infinite server, use sum of its callers' multiplicities
                    callers_of_caller = self._get_callers_of_task(tidx_caller)
                    if callers_of_caller:
                        njobs = sum(self._get_mult(c) * self._get_repl(c) for c in callers_of_caller
                                   if not np.isinf(self._get_mult(c) * self._get_repl(c)))
                    if njobs == 0 or np.isinf(njobs):
                        # fallback njobs heuristic capped at 1000 (not MATLAB's 1e6, to bound load-dependent MVA's O(prod(N+1)) state space).
                        lqn = self.lqn
                        mult_arr = lqn.mult.flatten() if hasattr(lqn, 'mult') and lqn.mult is not None else np.array([1.0])
                        repl_arr = lqn.repl.flatten() if hasattr(lqn, 'repl') and lqn.repl is not None else np.ones_like(mult_arr)
                        finite_mask = np.isfinite(mult_arr) & np.isfinite(repl_arr)
                        if np.any(finite_mask):
                            njobs = min(np.sum(mult_arr[finite_mask] * repl_arr[finite_mask]), 1000)
                        else:
                            njobs = 100

                self.njobs[tidx_caller, idx] = njobs

                caller_name = self._get_hashname(tidx_caller)

                # Create closed class for this caller
                caller_class = ClosedClass(layer_model, caller_name, int(njobs), client_delay)

                # client delay: host layers = think time only; task layers = think time + host demand; call response time handled by CALL classes.
                if is_host_layer:
                    # Host layer: TASK class client delay = think time only
                    # a served task's declared think time is not a per-request
                    # delay, so the seed carries none either -- see
                    # _ref_think_mean; update_layers replaces this from the
                    # first iteration on
                    think_time = self._ref_think_mean(tidx_caller)

                    # TASK class delay at client = think time only; non-REF tasks with no think time use Immediate (matches MATLAB Exp(0)).
                    if think_time > 0:
                        client_delay.set_service(caller_class, Exp.fit_mean(think_time))
                    else:
                        client_delay.set_service(caller_class, Immediate())
                else:
                    # Task layers: client service = caller's think time + host demand
                    # This represents time the caller spends NOT waiting for this server
                    # a served task's declared think time is not a per-request
                    # delay, so the seed carries none either -- see
                    # _ref_think_mean; update_layers replaces this from the
                    # first iteration on
                    think_time = self._ref_think_mean(tidx_caller)

                    # TASK class at client = think time only (matches MATLAB Layer-1 T1 rate); non-REF/no-think-time tasks use Immediate.
                    if think_time > 0:
                        client_delay.set_service(caller_class, Exp.fit_mean(think_time))
                    else:
                        client_delay.set_service(caller_class, Immediate())

                # Set service at server
                total_demand = 0.0

                if is_host_layer:
                    # Host layer: server is processor, service = caller's activities' host demands
                    activities = self._get_activities_of_task(tidx_caller)
                    for aidx in activities:
                        if self.servtproc[aidx] is not None:
                            proc = self.servtproc[aidx]
                            if isinstance(proc, (int, float, np.integer, np.floating)):
                                total_demand += float(proc)
                            elif hasattr(proc, 'getMean'):
                                total_demand += proc.getMean()
                            elif hasattr(proc, 'mean'):
                                total_demand += proc.mean
                else:
                    # Task layer: server is called task, service = called entry's service time
                    # This represents time the server spends processing this caller's request
                    total_demand = self._get_initial_call_response_time(tidx_caller, idx)

                # TASK class service at server is always Disabled in both HOST and TASK layers.
                for _srv in all_server_stations():
                    _srv.set_service(caller_class, Disabled())

                # Set class attribute (matches MATLAB class.attribute = [type, idx])
                caller_class.attribute = [LayeredNetworkElement.TASK, tidx_caller]
                caller_class.completes = False  # matches MATLAB line 188 and JAR line 688
                caller_class.setReferenceClass(True)  # renormalize residence times using the visits to the task (MATLAB buildLayersRecursive line 150)

                # Record task attribute
                layer_model.attribute['tasks'].append([caller_class.get_index(), tidx_caller])

                # servt_map only tracks ACTIVITY indices, not TASK; activity entries are added separately in host layers.

                # only non-REF tasks join thinkt_classes_updmap; REF tasks keep their user-specified think time unchanged (matches MATLAB line 154-155).
                if not self._is_ref_task(tidx_caller):
                    thinkt_map[idx].append([idx, tidx_caller, 1, caller_class.get_index()])

                # Create ENTRY classes for each entry of this caller (matches MATLAB buildLayersRecursive lines 158-173)
                entries = self._get_entries_of_task(tidx_caller)
                if 'entries' not in layer_model.attribute:
                    layer_model.attribute['entries'] = []
                source_station = layer_model.attribute.get('source_station')
                for eidx in entries:
                    entry_name = self._get_hashname(eidx)
                    entry_class = ClosedClass(layer_model, entry_name, 0, client_delay)

                    # ENTRY class: Immediate at client, Disabled at server
                    client_delay.set_service(entry_class, Immediate())
                    for _srv in all_server_stations():
                        _srv.set_service(entry_class, Disabled())

                    entry_class.attribute = [LayeredNetworkElement.ENTRY, eidx]
                    entry_class.completes = False

                    layer_model.attribute['entries'].append([entry_class.get_index(), eidx])

                    # entry open-arrival distribution creates OpenClass on layer's Source/Sink; mirrors MATLAB buildLayersRecursive.m:214-255/JAR SolverLN.java:718-750.
                    if source_station is not None and eidx in lqn.arrival and lqn.arrival[eidx] is not None \
                            and not self._is_open_arrival_only(tidx_caller):
                        open_class = OpenClass(layer_model, entry_name + '_Open', 0)
                        source_station.set_arrival(open_class, lqn.arrival[eidx])
                        client_delay.set_service(open_class, Disabled())

                        # entries have Immediate servtproc; use the first bound activity's host demand as the initial server estimate (refined via servt_classes_updmap).
                        bound_act_svc = None
                        if hasattr(lqn, 'graph') and lqn.graph is not None and isinstance(lqn.graph, np.ndarray):
                            if eidx < lqn.graph.shape[0]:
                                for cand in range(lqn.graph.shape[1]):
                                    if lqn.graph[eidx, cand] > 0:
                                        if cand < len(self.servtproc) and self.servtproc[cand] is not None:
                                            bound_act_svc = self.servtproc[cand]
                                            break
                        if bound_act_svc is None:
                            bound_act_svc = self.servtproc[eidx] if self.servtproc[eidx] is not None else Immediate()

                        all_servers = layer_model.attribute.get('server_stations', [server_station])
                        for srv in all_servers:
                            srv.set_service(open_class, bound_act_svc)

                        open_class.attribute = [LayeredNetworkElement.ENTRY, eidx]
                        open_class.completes = False
                        layer_model.attribute['entry_open_classes'].append((open_class, eidx))

                        # arvproc_classes_updmap uses negative eidx convention; mirrors JAR SolverLN.java:745 / MATLAB buildLayersRecursive.m:250.
                        src_node_idx = layer_model.get_nodes().index(source_station) + 1
                        arvproc_map[idx].append([idx, -eidx, src_node_idx, open_class.get_index()])

            # Create ACTIVITY classes for each activity of this caller (matches MATLAB buildLayersRecursive)
            activities = self._get_activities_of_task(tidx_caller)
            if 'activities' not in layer_model.attribute:
                layer_model.attribute['activities'] = []
            for aidx in activities:
                act_stations = servers_for(self._get_parent(self._get_parent(aidx))) if flat \
                    else (srv_stations.get(idx, []) if is_host_layer else [])
                # the activity's demand belongs on its own processor's station
                act_stn = act_stations[0] if act_stations else server_station
                if act_stations or any(self._has_sync_callers(_sid, callers) for _sid in idx_set):
                    activity_name = self._get_hashname(aidx)
                    activity_class = ClosedClass(layer_model, activity_name, 0, client_delay)

                    if act_stations:
                        # the activity runs on a server of this layer, so its demand sits there and the client is Disabled; mirrors MATLAB buildLayersRecursive ~line 236.
                        client_delay.set_service(activity_class, Disabled())
                        if aidx < len(self.servtproc) and self.servtproc[aidx] is not None:
                            base_proc = self.servtproc[aidx]
                            base_mean = 0.0
                            if hasattr(base_proc, 'getMean'):
                                base_mean = base_proc.getMean()
                            elif hasattr(base_proc, 'mean'):
                                base_mean = base_proc.mean
                            elif isinstance(base_proc, (int, float)):
                                base_mean = float(base_proc)

                            # A SetupTask's cold start is NOT wired into the layer
                            # station any more, and it is not folded into the host
                            # demand either: it is charged to the entry with
                            # probability p by _setup_charge. Wiring it here routed
                            # the layer through the open M/G/1-with-setup QBD, which
                            # powers the thread down far more often than a closed
                            # layer does, and charged the delay to the ACTIVITY,
                            # where it is not host demand.
                            # The host layer keeps the real host-demand Distribution
                            # (not an Exp fit), so moment2/moment3 do not collapse.
                            act_stn.set_service(activity_class, base_proc)
                        else:
                            act_stn.set_service(activity_class, Exp.fit_mean(0.001))
                    else:
                        # task layer: activities process at CLIENT with host demand (the caller's own activity time), mirrors MATLAB buildLayersRecursive.
                        if aidx < len(self.servtproc) and self.servtproc[aidx] is not None:
                            proc = self.servtproc[aidx]
                            if hasattr(proc, 'getMean'):
                                hostdem = proc.getMean()
                            elif hasattr(proc, 'mean'):
                                hostdem = proc.mean
                            else:
                                hostdem = 0.0
                            # Handle zero/negative hostdem - use Immediate for zero demand
                            if hostdem > 0:
                                client_delay.set_service(activity_class, Exp.fit_mean(hostdem))
                            else:
                                client_delay.set_service(activity_class, Immediate())
                        else:
                            client_delay.set_service(activity_class, Immediate())
                        act_stn.set_service(activity_class, Disabled())

                    activity_class.attribute = [LayeredNetworkElement.ACTIVITY, aidx]
                    activity_class.completes = False

                    layer_model.attribute['activities'].append([activity_class.get_index(), aidx])

                    # Add servt_map entry for activity classes in host layers (matches MATLAB line 484)
                    # servt_classes_updmap stores: [model_idx, activity_lqn_idx, node_idx, class_idx]
                    if act_stations:
                        _hidx = self._get_parent(self._get_parent(aidx)) if flat else idx
                        servt_map[idx].append([idx, aidx,
                                               layer_model.attribute['serverIdxOf'].get(int(_hidx),
                                                   layer_model.attribute['serverIdx']),
                                               activity_class.get_index()])
                    else:
                        # task layer: activity service at client updates from thinkt_map (host processor response time); mirrors MATLAB buildLayersRecursive:585-586.
                        thinkt_map[idx].append([idx, aidx, 1, activity_class.get_index()])

                    # host-layer-only aux think-time class: other layers carry think time in servtproc, so adding it here would double-charge it (U(P2) 9x high).
                    if (act_stations and hasattr(self, 'actthinkproc')
                            and aidx < len(self.actthinkproc)
                            and self.actthinkproc[aidx] is not None):
                        think_name = self._get_hashname(aidx) + '.Think'
                        think_class = ClosedClass(layer_model, think_name, 0, client_delay)
                        think_class.completes = False
                        think_class.attribute = [LayeredNetworkElement.ACTIVITY, aidx]
                        client_delay.set_service(think_class, self.actthinkproc[aidx])
                        for _srv in all_server_stations():
                            _srv.set_service(think_class, Disabled())
                        actthinkt_map[idx].append([idx, aidx, 1, think_class.get_index()])

            # Create CALL classes for sync calls from this caller's activities (matches MATLAB lines 287-302)
            if 'calls' not in layer_model.attribute:
                layer_model.attribute['calls'] = []
            _cgroup_by_cidx, _cgroup_members = self._call_groups_by_cidx()
            _group_classes = layer_model.attribute.setdefault('call_group_classes', {})
            for aidx in activities:
                if isinstance(self.lqn.callsof, dict):
                    calls = self.lqn.callsof.get(aidx, [])
                else:
                    calls = []

                for cidx in calls:
                    # Check if this is a SYNC call
                    is_sync = True
                    if hasattr(self.lqn, 'calltype') and self.lqn.calltype is not None:
                        if isinstance(self.lqn.calltype, np.ndarray):
                            # calltype is 1-indexed (like MATLAB), so use cidx directly
                            calltype = self.lqn.calltype.flatten()[cidx] if cidx < len(self.lqn.calltype.flatten()) else CallType.SYNC
                        elif isinstance(self.lqn.calltype, dict):
                            calltype = self.lqn.calltype.get(cidx, CallType.SYNC)
                        else:
                            calltype = CallType.SYNC
                        is_sync = (calltype == CallType.SYNC)

                    if is_sync:
                        # A routed call group is ONE dispatch with n destinations, not
                        # n calls: its members share the dispatch class, which is the
                        # class the strategy routes and the one that visits the targets
                        # (its per-station service carries the per-target service time).
                        # The strategy is a property of a NODE and routes over that
                        # node's links, so the choice is made at a router whose only
                        # links are the group's targets; the hop itself must not switch
                        # class, because a state-dependent routing function is evaluated
                        # at zero off the class diagonal. The class switch goes on the
                        # return arc, into a group class the job continues in.
                        _grp = _cgroup_by_cidx.get(cidx)
                        if _grp is not None:
                            _gid = _grp[0]
                            if _gid in _group_classes:
                                call_class = _group_classes[_gid][0]
                                _group_reuse = True
                            else:
                                _disp_name = self._get_hashname(aidx) + '.Dispatch%d' % _gid
                                _router = Router(layer_model, _disp_name + '.Router')
                                call_class = ClosedClass(layer_model, _disp_name, 0, client_delay)
                                client_delay.set_service(call_class, Immediate())
                                for _srv in all_server_stations():
                                    _srv.set_service(call_class, Disabled())
                                return_class = ClosedClass(
                                    layer_model, self._get_call_hashname(cidx) + '.Group%d' % _gid,
                                    0, client_delay)
                                return_class.completes = False
                                return_class.attribute = [LayeredNetworkElement.CALL, cidx]
                                client_delay.set_service(return_class, Immediate())
                                for _srv in all_server_stations():
                                    _srv.set_service(return_class, Disabled())
                                _group_classes[_gid] = (call_class, return_class, _grp[1], _router)
                                _group_reuse = False
                        else:
                            _group_reuse = False
                            call_name = self._get_call_hashname(cidx)
                            call_class = ClosedClass(layer_model, call_name, 0, client_delay)
                        layer_model.attribute.setdefault('call_class_of_cidx', {})[cidx] = call_class.get_index()

                        # Get call mean for Aux class creation (MATLAB lines 305-315)
                        call_mean = self._get_call_mean(cidx)
                        nreplicas = 1  # Typically 1, could be based on processor replication

                        # Create Aux class for fractional call means (matches MATLAB lines 308-314).
                        # A group member's mean is the 1/n share of the dispatch, which the
                        # n-way split already carries, so the Aux skip path must not also fire.
                        aux_class = None
                        if call_mean != 1 and _grp is None:
                            aux_name = call_name + '.Aux'
                            aux_class = ClosedClass(layer_model, aux_name, 0, client_delay)
                            aux_class.completes = False
                            aux_class.attribute = [LayeredNetworkElement.CALL, cidx]  # Same attribute as call class
                            client_delay.set_service(aux_class, Immediate())
                            for _srv in all_server_stations():
                                _srv.set_service(aux_class, Disabled())
                            # Track aux class: [class_index, cidx, call_mean]
                            if 'aux_classes' not in layer_model.attribute:
                                layer_model.attribute['aux_classes'] = []
                            layer_model.attribute['aux_classes'].append([aux_class.get_index(), cidx, call_mean])

                        # Get call service time (callservtproc)
                        tgt_eidx = self._get_call_target_entry(cidx)
                        tgt_tidx = self._get_parent(tgt_eidx) if tgt_eidx else None

                        # minRespT for server = sum of activities' hostdem in host layers (processor, no activities -> 0) or the server task's activities in task layers.
                        tgt_stations = servers_for(tgt_tidx)
                        seed_idx = tgt_tidx if flat else idx
                        if flat:
                            minRespT = self._get_initial_task_total_hostdem(tgt_tidx) if tgt_tidx else 0.0
                        elif is_host_layer:
                            # Host processor has no activities - minRespT = 0
                            minRespT = 0.0
                        else:
                            # Task layer: server is a task with activities
                            minRespT = self._get_initial_task_total_hostdem(idx) if idx else 0.0

                        # CALL class service times: a call to a task that is a server of
                        # this layer is served THERE, any other call is a delay at the client.
                        call_to_server = bool(tgt_stations)

                        if call_to_server:
                            # Call to this layer's server - service at SERVER
                            # MATLAB line 727: clientDelay.setService(cidxClass{cidx}, Immediate.getInstance())
                            client_delay.set_service(call_class, Immediate())
                            for _srv in tgt_stations:
                                if cidx < len(self.callservtproc) and self.callservtproc[cidx] is not None:
                                    _srv.set_service(call_class, self.callservtproc[cidx])
                                else:
                                    _srv.set_service(call_class, Immediate())
                            # Record the station of the called task for call_classes_updmap
                            call_map[idx].append([idx, cidx,
                                                  layer_model.attribute['serverIdxOf'][int(tgt_tidx)],
                                                  call_class.get_index()])
                        else:
                            # Call to another task - service at CLIENT (MATLAB lines 750, 804)
                            # MATLAB: clientDelay.setService(cidxClass{cidx}, callservtproc{cidx})
                            if cidx < len(self.callservtproc) and self.callservtproc[cidx] is not None:
                                client_delay.set_service(call_class, self.callservtproc[cidx])
                            else:
                                client_delay.set_service(call_class, Immediate())
                            # MATLAB keeps server at Exp.fitMean(minRespT) which is 1e-8 for minRespT=0
                            # This is set initially at lines 299-300 and NOT changed in the routing setup
                            for _srv in all_server_stations():
                                _srv.set_service(call_class, Exp.fit_mean(max(minRespT, 1e-8)))
                            # Record with clientIdx=1 for call_classes_updmap (MATLAB lines 751, 805)
                            call_map[idx].append([idx, cidx, 1, call_class.get_index()])

                        if not _group_reuse:
                            call_class.attribute = [LayeredNetworkElement.CALL, cidx]
                        call_class.completes = False

                        # Track call: [class_index, cidx, src_aidx, tgt_eidx, aux_class_index]
                        src_aidx = aidx
                        aux_class_idx = aux_class.get_index() if aux_class else -1
                        layer_model.attribute['calls'].append([call_class.get_index(), cidx, src_aidx, tgt_eidx if tgt_eidx else 0, aux_class_idx])

                        # SYNC forwarding-chain classes unnecessary: rewritten caller-side; see _kb/06-solver-catalog.md LN Forwarding as caller-side pseudo-rendezvous.
                    else:
                        # ASYNC call fires only when target entry's task matches this layer's server; mirrors MATLAB buildLayersRecursive.m:299-324/JAR SolverLN.java:780-814.
                        tgt_eidx_async = self._get_call_target_entry(cidx)
                        tgt_parent = self._get_parent(tgt_eidx_async) if tgt_eidx_async else None
                        if tgt_parent != idx:
                            continue

                        # Lazy-create Source/Sink for async-only layers
                        # (MATLAB line 301-306 hasSource branch).
                        source_station = layer_model.attribute.get('source_station')
                        sink_station = layer_model.attribute.get('sink_station')
                        if layer_model.attribute.get('has_fork', False):
                            raise ValueError(f"SolverLN: layer '{layer_model.getName()}' carries both an "
                                             "AND fork and an open stream (an async call or an entry "
                                             "arrival); the fork-join transform needs a Source of its own")
                        if source_station is None:
                            source_station = Source(layer_model, 'Source')
                            sink_station = Sink(layer_model, 'Sink')
                            layer_model.attribute['source_station'] = source_station
                            layer_model.attribute['sink_station'] = sink_station

                        call_name = self._get_call_hashname(cidx)
                        open_class = OpenClass(layer_model, call_name, 0)
                        # async-call Source arrival starts Immediate, refreshed from the caller activity's tputproc; mirrors JAR:789/MATLAB:308.
                        source_station.set_arrival(open_class, Immediate())
                        client_delay.set_service(open_class, Disabled())

                        # async-call server service starts Immediate, upper-bounded by sum of hostdem over the target's activities; mirrors MATLAB:317-323/JAR:801-812.
                        minRespT = 0.0
                        activities_of_idx = self._get_activities_of_task(idx)
                        for tidx_act in activities_of_idx:
                            if tidx_act < len(self.servtproc) and self.servtproc[tidx_act] is not None:
                                proc = self.servtproc[tidx_act]
                                if hasattr(proc, 'getMean'):
                                    minRespT += proc.getMean()
                                elif hasattr(proc, 'mean'):
                                    minRespT += proc.mean
                                elif isinstance(proc, (int, float, np.integer, np.floating)):
                                    minRespT += float(proc)

                        all_servers = layer_model.attribute.get('server_stations', [server_station])
                        if minRespT > 0:
                            srv_service = Exp.fit_mean(minRespT)
                        else:
                            srv_service = Immediate()
                        for srv in all_servers:
                            srv.set_service(open_class, srv_service)

                        open_class.attribute = [LayeredNetworkElement.CALL, cidx]
                        open_class.completes = False

                        # Stash for routing phase
                        call_mean_async = self._get_call_mean(cidx)
                        layer_model.attribute['async_open_classes'].append(
                            (open_class, cidx, call_mean_async)
                        )

                        # async calls use POSITIVE cidx in arvproc_classes_updmap; mirrors MATLAB:427/JAR:871.
                        src_node_idx = layer_model.get_nodes().index(source_station) + 1
                        arvproc_map[idx].append([idx, cidx, src_node_idx, open_class.get_index()])

                        # call_classes_updmap records server-side class for
                        # each replica (MATLAB 428-430, JAR 872-876).
                        for srv in all_servers:
                            srv_node_idx = layer_model.get_nodes().index(srv) + 1
                            call_map[idx].append([idx, cidx, srv_node_idx, open_class.get_index()])

                        # Track in layer calls list for consistency with SYNC path
                        layer_model.attribute['calls'].append(
                            [open_class.get_index(), cidx, aidx,
                             tgt_eidx_async if tgt_eidx_async else 0, -1]
                        )

        # Configure cache node for cache layers (MATLAB buildLayersRecursive.m lines 548-561)
        if is_host_layer and layer_model.attribute.get('iscachelayer') and layer_model.attribute.get('cacheNode'):
            self._configure_cache_node(layer_model, idx, callers)

        # Link the model with routing
        self._setup_routing(layer_model, idx, route_map)

    def _configure_cache_node(self, layer_model: Network, idx: int, callers: List[int]):
        """
        Configure the cache node with hit/miss classes and access probabilities.

        This matches MATLAB buildLayersRecursive.m lines 548-561:
        - setReadItemEntry: set item access probability for the entry class
        - setHitClass: map input class to hit output class
        - setMissClass: map input class to miss output class

        Args:
            layer_model: The layer Network containing the Cache node
            idx: Layer index (task/processor absolute index)
            callers: List of caller task indices on this layer
        """
        lqn = self.lqn
        cache_node = layer_model.attribute.get('cacheNode')
        if cache_node is None:
            return

        # Get the cache task index (first caller that is a cache task)
        cache_task_idx = None
        for caller_idx in callers:
            if hasattr(lqn, 'iscache') and lqn.iscache is not None:
                iscache_arr = lqn.iscache.flatten() if isinstance(lqn.iscache, np.ndarray) else lqn.iscache
                if caller_idx < len(iscache_arr) and iscache_arr[caller_idx]:
                    cache_task_idx = caller_idx
                    break

        if cache_task_idx is None:
            return

        # Find the ItemEntry associated with this cache task
        # In MATLAB: the entry bound to the cache activity has lqn.itemproc set
        item_entry_idx = None
        item_access_prob = None
        entries = self._get_entries_of_task(cache_task_idx)
        for eidx in entries:
            if hasattr(lqn, 'itemproc') and isinstance(lqn.itemproc, dict):
                if eidx in lqn.itemproc and lqn.itemproc[eidx] is not None:
                    item_entry_idx = eidx
                    item_access_prob = lqn.itemproc[eidx]
                    break

        if item_entry_idx is None:
            return

        # Find the cache entry activity (bound to ItemEntry)
        cache_entry_aidx = None
        activities = self._get_activities_of_task(cache_task_idx)
        for aidx in activities:
            # Check if this activity is bound to the item entry
            bound_entry = self._get_activity_bound_entry(aidx)
            if bound_entry == item_entry_idx:
                cache_entry_aidx = aidx
                break

        if cache_entry_aidx is None:
            return

        # Find hit/miss activities from the graph (successors of cache entry activity)
        # MATLAB: lqn.hitmissaidx = find(lqn.graph(nextaidx,:))
        hit_aidx = None
        miss_aidx = None
        if hasattr(lqn, 'graph') and lqn.graph is not None:
            successors = []
            for j in range(lqn.graph.shape[1]):
                if lqn.graph[cache_entry_aidx, j] != 0:
                    successors.append(j)
            # MATLAB convention: first successor is hit, second is miss
            # (matches buildLayersRecursive.m lines 552-553)
            if len(successors) >= 2:
                hit_aidx = successors[0]
                miss_aidx = successors[1]
            elif len(successors) == 1:
                # If only one successor, assume it's miss (cache always misses)
                miss_aidx = successors[0]

        # Build mapping from activity index to class object
        # layer_model.classes is a list (0-indexed internally but class.get_index() returns 1-indexed)
        classes_list = layer_model.classes if hasattr(layer_model, 'classes') else []
        activity_to_class = {}
        if 'activities' in layer_model.attribute:
            for class_info in layer_model.attribute['activities']:
                if len(class_info) >= 2:
                    class_idx, act_idx = class_info[0], class_info[1]
                    # class_idx is 1-indexed from get_index(), convert to 0-indexed for list access
                    list_idx = class_idx - 1 if class_idx > 0 else 0
                    if 0 <= list_idx < len(classes_list):
                        activity_to_class[act_idx] = classes_list[list_idx]

        # Get the cache entry class
        entry_class = activity_to_class.get(cache_entry_aidx)
        if entry_class is None:
            # Try to find entry class from 'entries' attribute
            if 'entries' in layer_model.attribute:
                for class_info in layer_model.attribute['entries']:
                    if len(class_info) >= 2:
                        class_idx, entry_idx = class_info[0], class_info[1]
                        if entry_idx == item_entry_idx:
                            list_idx = class_idx - 1 if class_idx > 0 else 0
                            if 0 <= list_idx < len(classes_list):
                                entry_class = classes_list[list_idx]
                            break

        if entry_class is None:
            return

        # Set up hit/miss classes
        hit_class = activity_to_class.get(hit_aidx) if hit_aidx else None
        miss_class = activity_to_class.get(miss_aidx) if miss_aidx else None

        if hit_class:
            cache_node.set_hit_class(entry_class, hit_class)
        if miss_class:
            cache_node.set_miss_class(entry_class, miss_class)

        # Set up access probability
        # The item_access_prob should be a DiscreteSampler or similar distribution
        if item_access_prob is not None:
            cache_node.set_read(entry_class, item_access_prob)

        # delayed-hit retrieval cache wiring (EXPERIMENTAL); see _kb/06-solver-catalog.md LN Delayed-hit retrieval cache wiring.
        hasretr = getattr(lqn, 'hasretrieval', None)
        if (hasretr is not None and cache_task_idx < hasretr.shape[0]
                and hasretr[cache_task_idx, 0] != 0
                and entry_class is not None and miss_class is not None and miss_aidx is not None):
            fetch = Queue(layer_model, str(cache_node.name) + '.Fetch', SchedStrategy.PS)
            layer_model.attribute['retrieval_wiring'] = {
                'read_class': entry_class, 'miss_class': miss_class,
                'miss_aidx': miss_aidx, 'fetch': fetch, 'cache_node': cache_node,
            }

    def _get_activity_bound_entry(self, aidx: int) -> Optional[int]:
        """Get the entry index that an activity is bound to."""
        lqn = self.lqn
        # Check replygraph - if activity replies to an entry, it's bound to that entry's task
        if hasattr(lqn, 'replygraph') and lqn.replygraph is not None:
            # replygraph is (nacts x nentries), rows are activities (relative index)
            act_rel = aidx - lqn.ashift  # Convert to relative activity index
            if 0 <= act_rel < lqn.replygraph.shape[0]:
                for e in range(lqn.replygraph.shape[1]):
                    entry_abs = e + lqn.eshift
                    # Check parent relationship as fallback
                    pass

        # Check graph for direct entry->activity edge (activity bound to entry)
        if hasattr(lqn, 'graph') and lqn.graph is not None:
            for eidx in range(lqn.eshift, lqn.ashift):
                if lqn.graph[eidx, aidx] != 0:
                    return eidx

        return None

    def _get_initial_call_response_time(self, caller_tidx: int, layer_idx: int) -> float:
        """
        Get initial call response time estimate for a caller in a layer.

        For the first iteration, this uses the host demand of called entries.
        After iterations start, this is updated with actual response times.
        """
        lqn = self.lqn
        total_call_time = 0.0

        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return 0.0

        # Find all synch calls from this caller's activities
        activities = self._get_activities_of_task(caller_tidx)
        for aidx in activities:
            if isinstance(lqn.callsof, dict):
                calls = lqn.callsof.get(aidx, [])
            else:
                calls = []

            for cidx in calls:
                # Check call type - assume SYNC if calltype not available
                is_sync = True
                if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    if isinstance(lqn.calltype, np.ndarray):
                        calltype = lqn.calltype.flatten()[cidx] if cidx < len(lqn.calltype.flatten()) else CallType.SYNC
                    elif isinstance(lqn.calltype, dict):
                        calltype = lqn.calltype.get(cidx, CallType.SYNC)
                    else:
                        calltype = CallType.SYNC
                    is_sync = (calltype == CallType.SYNC)

                if is_sync:
                    # Get target entry (column 2 of callpair)
                    tgt_eidx = self._get_call_target_entry(cidx)
                    if tgt_eidx is None or tgt_eidx == 0:
                        continue

                    # Check if this call targets the server in this layer
                    tgt_tidx = self._get_parent(tgt_eidx)
                    if tgt_tidx != layer_idx:
                        continue

                    # Get call mean (number of calls)
                    call_mean = self._get_call_mean(cidx)

                    # Get initial response time = entry service time (recursive)
                    entry_resp = self._get_initial_entry_service_time(tgt_eidx, visited=set())
                    total_call_time += call_mean * entry_resp

        return total_call_time

    def _get_initial_entry_service_time(self, eidx: int, visited: set = None) -> float:
        """
        Compute initial entry service time recursively.

        Includes:
        - Sum of activities' host demands bound to this entry
        - Plus call_mean * target_entry_service_time for all downstream synch calls

        Uses memoization via visited set to avoid infinite loops.
        """
        if visited is None:
            visited = set()

        if eidx in visited:
            return 0.0  # Avoid infinite recursion
        visited.add(eidx)

        lqn = self.lqn
        total_time = 0.0

        # Get activities bound to this entry
        tgt_activities = self._get_activities_of_entry(eidx)

        for aidx in tgt_activities:
            # Add activity's host demand
            if self.servtproc[aidx] is not None:
                proc = self.servtproc[aidx]
                if isinstance(proc, (int, float, np.integer, np.floating)):
                    total_time += float(proc)
                elif hasattr(proc, 'getMean'):
                    total_time += proc.getMean()
                elif hasattr(proc, 'mean'):
                    total_time += proc.mean

            # Add downstream call response times
            if isinstance(lqn.callsof, dict):
                calls = lqn.callsof.get(aidx, [])
            else:
                calls = []

            for cidx in calls:
                # Check if this is a synch call
                is_sync = True
                if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    if isinstance(lqn.calltype, np.ndarray):
                        calltype = lqn.calltype.flatten()[cidx] if cidx < len(lqn.calltype.flatten()) else CallType.SYNC
                    elif isinstance(lqn.calltype, dict):
                        calltype = lqn.calltype.get(cidx, CallType.SYNC)
                    else:
                        calltype = CallType.SYNC
                    is_sync = (calltype == CallType.SYNC)

                if is_sync:
                    # Get target entry
                    target_eidx = self._get_call_target_entry(cidx)
                    if target_eidx is not None and target_eidx > 0:
                        call_mean = self._get_call_mean(cidx)
                        # Recursively get target entry's service time
                        target_resp = self._get_initial_entry_service_time(target_eidx, visited.copy())
                        total_time += call_mean * target_resp

        return total_time

    def _get_initial_task_total_hostdem(self, tidx: int) -> float:
        """
        Get total host demand of all activities of a task.

        This matches MATLAB buildLayersRecursive lines 296-302:
            minRespT = 0;
            for tidx_act = lqn.actsof{idx}
                minRespT = minRespT + lqn.hostdem{tidx_act}.getMean;
            end

        This provides an upper bound on the task's response time for
        initial service time estimates.
        """
        lqn = self.lqn
        total_hostdem = 0.0

        # Get all activities of this task
        activities = self._get_activities_of_task(tidx)

        for aidx in activities:
            # Use lqn.hostdem (host CPU demand) NOT self.servtproc (service time)
            # This matches MATLAB's lqn.hostdem{tidx_act}.getMean
            if aidx in lqn.hostdem:
                hostdem_val = lqn.hostdem[aidx]
                if isinstance(hostdem_val, (int, float, np.integer, np.floating)):
                    total_hostdem += float(hostdem_val)
                elif hasattr(hostdem_val, 'getMean'):
                    total_hostdem += hostdem_val.getMean()
                elif hasattr(hostdem_val, 'mean'):
                    total_hostdem += hostdem_val.mean
                elif hasattr(hostdem_val, 'get_mean'):
                    total_hostdem += hostdem_val.get_mean()

        return total_hostdem

    def _get_entry_service_matrix(self) -> np.ndarray:
        """
        Build entry service matrix (matches MATLAB getEntryServiceMatrix).

        Returns a matrix U of shape (nidx + ncalls, nidx + ncalls) where:
        - U[eidx, aidx] = probability that activity aidx contributes to entry eidx's service time
        - U[eidx, nidx + cidx] = probability that call cidx contributes to entry eidx's service time

        The entry service time is then computed as:
            entry_servt = U @ [residt; callresidt]

        NOTE: Unlike MATLAB which binarizes the matrix, Python preserves the probabilities
        from the LQN graph. This is because Python doesn't have full CacheNode support,
        so the hit/miss probabilities need to be applied via the servtmatrix.
        """
        lqn = self.lqn
        size = lqn.nidx + lqn.ncalls
        U = np.zeros((size, size))

        # For each entry, recursively trace the activity graph
        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            self._entry_service_matrix_recursion(eidx, eidx, U, 1.0)

        # Binarize the matrix (matches MATLAB: U = double(U > 0))
        # This prevents accumulation of probabilities from multiple paths
        U = (U > 0).astype(float)

        return U

    def _entry_service_matrix_recursion(self, aidx: int, eidx: int, U: np.ndarray, prob: float = 1.0, visited: set = None):
        """
        Auxiliary function to build entry service matrix recursively.

        Traverses the activity graph from aidx, marking all activities and calls
        that contribute to entry eidx's service time, weighted by probability.
        Uses a visited set to detect and break cycles in the activity graph.

        Args:
            aidx: Current activity index
            eidx: Entry index we're building service time for
            U: Service matrix to update
            prob: Cumulative probability of reaching this activity from entry
            visited: Set of already-visited activity indices for cycle detection
        """
        if visited is None:
            visited = set()
        visited = visited | {aidx}

        lqn = self.lqn
        graph = lqn.graph

        # Find next activities in the graph
        if aidx >= len(graph):
            return

        # Get all successors of current activity
        nextaidxs = np.where(graph[aidx, :] > 0)[0]

        for nextaidx in nextaidxs:
            # Check if this is a loop edge (graph differs from dag)
            # MATLAB: isLoop = (lqn.graph(aidx,nextaidx) ~= lqn.dag(aidx,nextaidx))
            is_loop = False
            if hasattr(lqn, 'dag') and lqn.dag is not None:
                if isinstance(lqn.dag, np.ndarray) and aidx < lqn.dag.shape[0] and nextaidx < lqn.dag.shape[1]:
                    is_loop = (graph[aidx, nextaidx] != lqn.dag[aidx, nextaidx])

            # Detect cycles: skip if we've already visited this node in the current path
            if nextaidx in visited:
                is_loop = True

            # Get parent of current and next nodes
            parent_aidx = self._get_parent(aidx)
            parent_nextaidx = self._get_parent(nextaidx)

            # Get edge probability
            edge_prob = graph[aidx, nextaidx]
            # Cumulative probability = path probability * edge probability
            next_prob = prob * edge_prob

            # If parents differ, this is a call to another task/entry
            if parent_aidx != parent_nextaidx:
                # Process calls from this activity
                if isinstance(lqn.callsof, dict):
                    calls = lqn.callsof.get(aidx, [])
                else:
                    calls = []

                for cidx in calls:
                    # Check call type - only SYNC calls contribute to response time
                    is_sync = True
                    if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                        if isinstance(lqn.calltype, np.ndarray):
                            # calltype is 1-indexed, so use cidx directly
                            if cidx < len(lqn.calltype.flatten()):
                                calltype = lqn.calltype.flatten()[cidx]
                                is_sync = (calltype == CallType.SYNC)

                    if is_sync:
                        # U(eidx,nidx+cidx)=1: mean number of calls is already factored into callresidt via visits.
                        U[eidx, lqn.nidx + cidx] = 1

            # If parents are the same, this is an activity within the same task
            if parent_aidx == parent_nextaidx:
                if nextaidx != aidx and not is_loop:
                    # Mark activity as contributing to entry with cumulative probability
                    # Use max to handle multiple paths to same activity
                    U[eidx, nextaidx] = max(U[eidx, nextaidx], next_prob)
                    # Recurse to process the rest of the graph
                    self._entry_service_matrix_recursion(nextaidx, eidx, U, next_prob, visited)

    def _get_initial_call_response_time_for_task(self, tidx: int) -> float:
        """
        Get initial total call response time for a task.

        This is the sum of (call_mean * target_entry_service_time) for all synch calls
        from this task's activities.
        """
        lqn = self.lqn
        total_call_time = 0.0

        # Find all activities of this task
        activities = self._get_activities_of_task(tidx)

        for aidx in activities:
            # Get calls from this activity
            if isinstance(lqn.callsof, dict):
                calls = lqn.callsof.get(aidx, [])
            else:
                calls = []

            for cidx in calls:
                # Check if this is a synch call
                is_sync = True
                if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    if isinstance(lqn.calltype, np.ndarray):
                        calltype = lqn.calltype.flatten()[cidx] if cidx < len(lqn.calltype.flatten()) else CallType.SYNC
                    elif isinstance(lqn.calltype, dict):
                        calltype = lqn.calltype.get(cidx, CallType.SYNC)
                    else:
                        calltype = CallType.SYNC
                    is_sync = (calltype == CallType.SYNC)

                if is_sync:
                    # Get target entry
                    target_eidx = self._get_call_target_entry(cidx)
                    if target_eidx is not None and target_eidx > 0:
                        call_mean = self._get_call_mean(cidx)
                        # Recursively get target entry's service time
                        target_resp = self._get_initial_entry_service_time(target_eidx, visited=set())
                        total_call_time += call_mean * target_resp

        return total_call_time

    def _get_host_layer_response_time(self, tidx: int) -> float:
        """
        Get the host layer response time for a task.

        This is the response time at the processor from the host layer results.
        Falls back to host demand if results not available.
        """
        lqn = self.lqn

        # Find the host of this task
        hidx = self._get_parent(tidx)
        if hidx is None or hidx == 0:
            return self._get_task_total_host_demand(tidx)

        # Get host layer results
        if np.isnan(self.idxhash[hidx]):
            return self._get_task_total_host_demand(tidx)

        host_layer_idx = int(self.idxhash[hidx])
        if len(self.results) == 0 or host_layer_idx >= len(self.results[-1]):
            return self._get_task_total_host_demand(tidx)

        result = self.results[-1][host_layer_idx]
        if result is None or 'RN' not in result:
            return self._get_task_total_host_demand(tidx)

        RN = result['RN']
        server_idx = self.ensemble[host_layer_idx].attribute.get('serverIdx', 1)
        if server_idx is None:
            return self._get_task_total_host_demand(tidx)

        server_idx_0 = server_idx - 1 if server_idx >= 1 else 0
        if server_idx_0 >= RN.shape[0]:
            return self._get_task_total_host_demand(tidx)

        # Find this task's class in the host layer
        caller_class_idx = self._find_caller_class_in_layer(tidx, host_layer_idx)
        if caller_class_idx is not None:
            caller_class_idx_0 = caller_class_idx - 1 if caller_class_idx >= 1 else 0
            if caller_class_idx_0 < RN.shape[1]:
                return RN[server_idx_0, caller_class_idx_0]

        # Fallback to host demand
        return self._get_task_total_host_demand(tidx)

    def _get_activities_of_entry(self, eidx: int) -> List[int]:
        """Get all activities belonging to an entry.

        This includes:
        - The activity directly bound to the entry (edge from entry to activity)
        - All successor activities reachable via precedence edges (until reaching
          an activity that makes a call or replies to an entry)
        """
        lqn = self.lqn
        activities = []

        # Get task of this entry
        tidx = self._get_parent(eidx)
        if tidx is None:
            return activities

        # Get all activities of the task
        all_activities = set(self._get_activities_of_task(tidx))

        if not hasattr(lqn, 'graph') or lqn.graph is None:
            return activities

        # Find the activity directly bound to this entry
        bound_activity = None
        for aidx in all_activities:
            if isinstance(lqn.graph, np.ndarray):
                if eidx < lqn.graph.shape[0] and aidx < lqn.graph.shape[1]:
                    if lqn.graph[eidx, aidx] > 0:
                        bound_activity = aidx
                        break

        if bound_activity is None:
            return activities

        # Follow precedence chain from bound activity
        # Use BFS to find all reachable activities within this task
        visited = set()
        queue = [bound_activity]
        while queue:
            aidx = queue.pop(0)
            if aidx in visited:
                continue
            visited.add(aidx)
            activities.append(aidx)

            # Find successor activities (in the same task)
            if isinstance(lqn.graph, np.ndarray) and aidx < lqn.graph.shape[0]:
                for succ in range(lqn.graph.shape[1]):
                    if lqn.graph[aidx, succ] > 0:
                        # Check if successor is an activity in the same task
                        if succ in all_activities and succ not in visited:
                            queue.append(succ)

        return activities

    def _get_mult(self, idx: int) -> float:
        """Get multiplicity (job count) for an element.

        Matches MATLAB buildLayersRecursive line 7-8, 129:
          mult = lqn.maxmult; % this removes spare capacity that cannot be used
          lqn.mult = mult;
          ...
          njobs = mult(tidx_caller)*lqn.repl(tidx_caller);

        Uses maxmult because MATLAB replaces mult with maxmult at start of
        buildLayersRecursive to "remove spare capacity that cannot be used".
        """
        lqn = self.lqn

        # Use maxmult if available (MATLAB line 7: mult = lqn.maxmult)
        if hasattr(lqn, 'maxmult') and lqn.maxmult is not None:
            if isinstance(lqn.maxmult, dict):
                return lqn.maxmult.get(idx, 1)
            elif isinstance(lqn.maxmult, np.ndarray):
                flat_maxmult = lqn.maxmult.flatten()
                if idx < len(flat_maxmult):
                    return float(flat_maxmult[idx])

        # Fallback to mult if maxmult not available
        if hasattr(lqn, 'mult') and lqn.mult is not None:
            if isinstance(lqn.mult, dict):
                return lqn.mult.get(idx, 1)
            elif isinstance(lqn.mult, np.ndarray):
                # Handle 2D arrays (shape like (1, n))
                flat_mult = lqn.mult.flatten()
                if idx < len(flat_mult):
                    return float(flat_mult[idx])
        return 1.0

    def _get_repl(self, idx: int) -> float:
        """Get replication factor for an element (matches MATLAB lqn.repl)."""
        lqn = self.lqn
        if hasattr(lqn, 'repl') and lqn.repl is not None:
            if isinstance(lqn.repl, dict):
                return lqn.repl.get(idx, 1)
            elif isinstance(lqn.repl, np.ndarray):
                flat_repl = lqn.repl.flatten()
                if idx < len(flat_repl):
                    val = float(flat_repl[idx])
                    return val if val > 0 else 1.0
        return 1.0

    def _get_activities_of_task(self, tidx: int) -> List[int]:
        """Get activity indices for a task."""
        lqn = self.lqn
        if isinstance(lqn.actsof, dict):
            return lqn.actsof.get(tidx, [])
        return []

    def _get_entries_of_task(self, tidx: int) -> List[int]:
        """Get entry indices for a task."""
        lqn = self.lqn
        if isinstance(lqn.entriesof, dict):
            return lqn.entriesof.get(tidx, [])
        return []

    def _get_call_hashname(self, cidx: int) -> str:
        """Get hash name for a call (e.g., 'AS2=>E:E2' for sync calls)."""
        lqn = self.lqn
        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return f'Call_{cidx}'

        # callpair format: [src_aidx, tgt_eidx, mean] (columns 0 and 1 are src and tgt)
        if isinstance(lqn.callpair, np.ndarray):
            if cidx < len(lqn.callpair) and lqn.callpair.ndim > 1:
                src_aidx = int(lqn.callpair[cidx, 0])  # Column 0 = source activity
                tgt_eidx = int(lqn.callpair[cidx, 1])  # Column 1 = target entry
            else:
                return f'Call_{cidx}'
        elif isinstance(lqn.callpair, dict):
            pair = lqn.callpair.get(cidx, [0, 0, 0, 0])
            src_aidx = int(pair[0]) if len(pair) > 0 else 0
            tgt_eidx = int(pair[1]) if len(pair) > 1 else 0
        else:
            return f'Call_{cidx}'

        # Get call type (default to SYNC if calltype not available)
        calltype = CallType.SYNC
        if hasattr(lqn, 'calltype') and lqn.calltype is not None:
            if isinstance(lqn.calltype, np.ndarray):
                if cidx < len(lqn.calltype.flatten()):
                    calltype = lqn.calltype.flatten()[cidx]
            elif isinstance(lqn.calltype, dict):
                calltype = lqn.calltype.get(cidx, CallType.SYNC)

        # Get names
        src_name = self._get_hashname(src_aidx)
        tgt_name = self._get_hashname(tgt_eidx)

        # Format based on call type
        if calltype == CallType.SYNC:
            return f'{src_name}=>{tgt_name}'
        elif calltype == CallType.ASYNC:
            return f'{src_name}->{tgt_name}'
        else:
            return f'{src_name}~>{tgt_name}'

    def _get_task_total_host_demand(self, tidx: int) -> float:
        """Get total host demand for a task (sum of all activities' host demands)."""
        total = 0.0
        activities = self._get_activities_of_task(tidx)
        for aidx in activities:
            if self.servtproc[aidx] is not None:
                proc = self.servtproc[aidx]
                if isinstance(proc, (int, float, np.integer, np.floating)):
                    total += float(proc)
                elif hasattr(proc, 'getMean'):
                    total += proc.getMean()
                elif hasattr(proc, 'mean'):
                    total += proc.mean
        return total

    # Constants for jobPos tracking in recurActGraph
    _AT_CLIENT = 1
    _AT_SERVER = 2
    _AT_CACHE = 3

    def _recur_act_graph(self, P, tidx_caller, aidx, cur_class, job_pos, ctx):
        """
        Recursively traverse the activity graph and set up routing.
        Matches MATLAB recurActGraph in buildLayersRecursive.m.

        Args:
            P: RoutingMatrix
            tidx_caller: Task index of the calling task
            aidx: Current activity/entry index
            cur_class: Current class object
            job_pos: Current position (_AT_CLIENT, _AT_SERVER, _AT_CACHE)
            ctx: Context dict with layer info (nodes, classes, fork/join state)

        Returns:
            (P, cur_class, job_pos)
        """
        lqn = self.lqn
        graph = lqn.graph

        # Save current state (MATLAB line 427-428)
        ctx['job_pos_key'][aidx] = job_pos
        ctx['cur_class_key'][aidx] = cur_class

        # Find successors (MATLAB line 430)
        nextaidxs = []
        if isinstance(graph, np.ndarray) and aidx < graph.shape[0]:
            for j in range(graph.shape[1]):
                if graph[aidx, j] != 0:
                    nextaidxs.append(j)

        # Check if any successor is POST_AND (fork target) (MATLAB line 432-433)
        is_post_and_act = ctx['is_post_and_act']
        is_pre_and_act = ctx['is_pre_and_act']
        is_next_prec_fork = any(n in is_post_and_act for n in nextaidxs)

        if not nextaidxs:
            return P, cur_class, job_pos

        # Pre-fork state, captured at the first branch so that calls this
        # activity issues before the fork stay sequential
        # (MATLAB buildLayersRecursive.m lines 522-525)
        fork_saved = False
        fork_save_cur_class = cur_class
        fork_save_job_pos = job_pos
        fork_save_station = ctx.get('cur_station')

        for nextaidx in nextaidxs:
            # Restore pre-fork state at start of each branch
            # (MATLAB buildLayersRecursive.m lines 531-534)
            if is_next_prec_fork:
                if not fork_saved:
                    if nextaidx in is_post_and_act:
                        fork_saved = True
                        fork_save_cur_class = cur_class
                        fork_save_job_pos = job_pos
                        fork_save_station = ctx.get('cur_station')
                else:
                    cur_class = fork_save_cur_class
                    job_pos = fork_save_job_pos
                    ctx['cur_station'] = fork_save_station
            # Loop detection (MATLAB line 440-442)
            is_loop = False
            if hasattr(lqn, 'dag') and lqn.dag is not None:
                if isinstance(lqn.dag, np.ndarray) and aidx < lqn.dag.shape[0] and nextaidx < lqn.dag.shape[1]:
                    is_loop = (graph[aidx, nextaidx] != lqn.dag[aidx, nextaidx])

            parent_aidx = self._get_parent(aidx)
            parent_nextaidx = self._get_parent(nextaidx)

            if parent_aidx != parent_nextaidx:
                # cross-task call routing mirrors MATLAB routeSynchCall in buildLayersRecursive.m:821-947.
                call_classes = ctx.get('call_classes', {})
                aux_classes = ctx.get('aux_classes', {})
                call_mean_map = ctx.get('call_mean_map', {})
                think_classes = ctx.get('think_classes', {})
                client_delay = ctx['client_delay']
                server_station = ctx['server_station']
                is_host_layer = ctx['is_host_layer']
                layer_idx = ctx.get('idx')

                # Find cidx by matching callpair [aidx, nextaidx]
                cidx = None
                if isinstance(lqn.callsof, dict):
                    for c in lqn.callsof.get(aidx, []):
                        if c < lqn.callpair.shape[0]:
                            pair = lqn.callpair[c]
                            src = int(pair[0]) if pair.shape[0] > 0 else -1
                            tgt = int(pair[1]) if pair.shape[0] > 1 else -1
                            if src == aidx and tgt == nextaidx:
                                cidx = c
                                break

                calltype = 1  # default SYNC
                if cidx is not None and hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    flat_ct = lqn.calltype.flatten() if isinstance(lqn.calltype, np.ndarray) else None
                    if flat_ct is not None and cidx < len(flat_ct):
                        calltype = int(flat_ct[cidx])
                    elif isinstance(lqn.calltype, dict):
                        calltype = int(lqn.calltype.get(cidx, 1))

                call_cls = call_classes.get(cidx) if cidx is not None else None

                if cidx is not None and call_cls is not None and calltype == 1:  # SYNC
                    call_mean = call_mean_map.get(cidx, 1.0)
                    aux_cls = aux_classes.get(cidx)
                    # nreplicas for this layer (task layers are single-replica)
                    nreplicas = 1
                    # Target entry's parent task — is it the server of this layer?
                    tgt_eidx_c = int(lqn.callpair[cidx, 1]) if cidx < lqn.callpair.shape[0] else None
                    tgt_parent = self._get_parent(tgt_eidx_c) if tgt_eidx_c is not None else None
                    _srv_map = ctx.get('srv_stations', {})
                    if ctx.get('flat'):
                        tgt_stn = _srv_map[int(tgt_parent)][0] if int(tgt_parent) in _srv_map else None
                        call_to_server = tgt_stn is not None
                        if call_to_server:
                            server_station = tgt_stn
                    else:
                        call_to_server = (tgt_parent == layer_idx)

                    callservt_proc = None
                    if cidx < len(self.callservtproc):
                        callservt_proc = self.callservtproc[cidx]

                    if job_pos == self._AT_CLIENT:
                        if call_to_server:
                            # MATLAB routeSynchCall atClient, call to server (lines 823-861)
                            if call_mean < 1:
                                if aux_cls is not None:
                                    P.set(cur_class, aux_cls, client_delay, client_delay, 1 - call_mean)
                                P.set(cur_class, call_cls, client_delay, server_station, call_mean / nreplicas)
                                P.set(call_cls, call_cls, server_station, client_delay, 1.0)
                                if aux_cls is not None:
                                    P.set(aux_cls, call_cls, client_delay, client_delay, 1.0)
                            elif call_mean == 1:
                                P.set(cur_class, call_cls, client_delay, server_station, 1.0 / nreplicas)
                                P.set(call_cls, call_cls, server_station, client_delay, 1.0)
                            else:  # call_mean > 1
                                P.set(cur_class, call_cls, client_delay, server_station, 1.0 / nreplicas)
                                if aux_cls is not None:
                                    P.set(call_cls, aux_cls, server_station, client_delay, 1.0)
                                    P.set(aux_cls, call_cls, client_delay, server_station,
                                          (1.0 - 1.0 / call_mean) / nreplicas)
                                    P.set(aux_cls, call_cls, client_delay, client_delay, 1.0 / call_mean)
                            # Services: Immediate at client, callservt at server
                            client_delay.set_service(call_cls, Immediate())
                            if callservt_proc is not None:
                                server_station.set_service(call_cls, callservt_proc)
                            job_pos = self._AT_CLIENT
                            cur_class = call_cls
                        else:
                            # MATLAB routeSynchCall atClient, call NOT to server (lines 863-879)
                            if call_mean < 1:
                                # call mean is embedded in the demand; see _kb/06-solver-catalog.md LN Call mean embedded in demand.
                                P.set(cur_class, call_cls, client_delay, client_delay, 1.0)
                                if aux_cls is not None:
                                    P.set(call_cls, aux_cls, client_delay, client_delay, 1.0)
                                    cur_class = aux_cls
                                else:
                                    cur_class = call_cls
                            elif call_mean == 1:
                                P.set(cur_class, call_cls, client_delay, client_delay, 1.0)
                                cur_class = call_cls
                            else:  # call_mean > 1
                                P.set(cur_class, call_cls, client_delay, client_delay, 1.0)
                                if aux_cls is not None:
                                    P.set(call_cls, aux_cls, client_delay, client_delay, 1.0)
                                    cur_class = aux_cls
                                else:
                                    cur_class = call_cls
                            if callservt_proc is not None:
                                client_delay.set_service(call_cls, callservt_proc)
                            job_pos = self._AT_CLIENT
                    else:  # job_pos == _AT_SERVER
                        if call_to_server:
                            # MATLAB routeSynchCall atServer, call to server (lines 882-910)
                            _from = (ctx.get('cur_station') or server_station) if ctx.get('flat') else server_station
                            if call_mean < 1:
                                # The skip flow must enter the Aux class and the reply
                                # must transit the client in the call class, which is
                                # therefore declared there: sn_refresh_visits drops any
                                # (station, class) state whose rate is NaN, and dropping
                                # this one severs the chain. Routing the skip into the
                                # call class instead and leaving in Aux gives Aux no
                                # inbound arc at all, so its chain has no reference
                                # class; this branch did that under 'srvn' until
                                # 2026-08-11 (buildLayersRecursive.m:1100-1118, and the
                                # JAR has carried the reference form all along).
                                if aux_cls is not None:
                                    P.set(cur_class, aux_cls, _from, client_delay, 1 - call_mean)
                                    P.set(aux_cls, call_cls, client_delay, client_delay, 1.0)
                                else:
                                    P.set(cur_class, call_cls, _from, client_delay, 1 - call_mean)
                                P.set(cur_class, call_cls, _from, server_station, call_mean)
                                P.set(call_cls, call_cls, server_station, client_delay, 1.0)
                                client_delay.set_service(call_cls, Immediate())
                                # both the skip and the visit end in the call class
                                cur_class = call_cls
                                job_pos = self._AT_CLIENT
                                ctx['cur_station'] = None
                            elif call_mean == 1:
                                P.set(cur_class, call_cls, _from, server_station, 1.0)
                                if ctx.get('flat'):
                                    # the reply returns the job to the client, as the
                                    # successor restoration downstream assumes
                                    P.set(call_cls, call_cls, server_station, client_delay, 1.0)
                                    client_delay.set_service(call_cls, Immediate())
                                    job_pos = self._AT_CLIENT
                                    ctx['cur_station'] = None
                                else:
                                    job_pos = self._AT_SERVER
                                    ctx['cur_station'] = server_station
                                cur_class = call_cls
                            else:  # call_mean > 1
                                P.set(cur_class, call_cls, _from, server_station, 1.0)
                                if ctx.get('flat'):
                                    # the geometric repeat transits the client between
                                    # visits; a self-loop would merge them into one
                                    if aux_cls is not None:
                                        P.set(call_cls, aux_cls, server_station, client_delay, 1.0)
                                        P.set(aux_cls, call_cls, client_delay, server_station, 1 - 1.0 / call_mean)
                                        P.set(aux_cls, call_cls, client_delay, client_delay, 1.0 / call_mean)
                                        client_delay.set_service(call_cls, Immediate())
                                    cur_class = call_cls
                                else:
                                    if aux_cls is not None:
                                        P.set(call_cls, call_cls, server_station, server_station, 1 - 1.0 / call_mean)
                                        P.set(call_cls, aux_cls, server_station, client_delay, 1.0 / call_mean)
                                    cur_class = aux_cls if aux_cls is not None else call_cls
                                job_pos = self._AT_CLIENT
                                ctx['cur_station'] = None
                            if callservt_proc is not None:
                                server_station.set_service(call_cls, callservt_proc)
                        else:
                            # MATLAB routeSynchCall atServer, call NOT to server (lines 912-936)
                            if call_mean < 1:
                                P.set(cur_class, call_cls, server_station, client_delay, 1.0)
                                if aux_cls is not None:
                                    P.set(call_cls, aux_cls, client_delay, client_delay, 1.0)
                                    cur_class = aux_cls
                                else:
                                    cur_class = call_cls
                            elif call_mean == 1:
                                P.set(cur_class, call_cls, server_station, client_delay, 1.0)
                                cur_class = call_cls
                            else:  # call_mean > 1
                                P.set(cur_class, call_cls, server_station, client_delay, 1.0)
                                if aux_cls is not None:
                                    P.set(call_cls, aux_cls, client_delay, client_delay, 1.0)
                                    cur_class = aux_cls
                                else:
                                    cur_class = call_cls
                            if callservt_proc is not None:
                                client_delay.set_service(call_cls, callservt_proc)
                            job_pos = self._AT_CLIENT

                    # (forwarding handled via pseudo rendezvous calls)
            else:
                # Same-task intra-activity routing (MATLAB lines 503-666)
                client_delay = ctx['client_delay']
                server_station = ctx['server_station']
                is_host_layer = ctx['is_host_layer']
                is_cache_layer = ctx['is_cache_layer']
                srv_stations = ctx.get('srv_stations', {})
                # station of the processor the next activity runs on, None when
                # that processor is not a server of this layer
                host_stn = None
                if ctx.get('flat'):
                    _h = self._get_parent(self._get_parent(nextaidx))
                    if _h is not None and int(_h) in srv_stations:
                        host_stn = srv_stations[int(_h)][0]
                elif is_host_layer:
                    host_stn = server_station
                fork_node = ctx['fork_node']
                join_node = ctx['join_node']
                fork_output_routers = ctx['fork_output_routers']
                fork_class_stack = ctx['fork_class_stack']
                activity_classes = ctx['activity_classes']

                act_cls = activity_classes.get(nextaidx)
                if act_cls is None:
                    continue

                # Check if any successor is an entry (MATLAB lines 1010-1021)
                entry_range = set(lqn.eshift + i for i in range(lqn.nentries))
                intersects = any(n in entry_range for n in nextaidxs)

                if not intersects:
                    # Restore state from saved values (MATLAB line 1023-1025)
                    job_pos = ctx['job_pos_key'].get(aidx, job_pos)
                    cur_class = ctx['cur_class_key'].get(aidx, cur_class)
                    ctx['cur_station'] = ctx.setdefault('cur_station_key', {}).get(aidx, ctx.get('cur_station'))
                else:
                    # Entry routing state restoration (MATLAB lines 1026-1040)
                    idx_in_nextaidxs = nextaidxs.index(nextaidx) if nextaidx in nextaidxs else 0
                    is_member = False
                    if idx_in_nextaidxs > 0:
                        prev_val = nextaidxs[idx_in_nextaidxs - 1]
                        is_member = prev_val in entry_range
                    if is_member:
                        ctx['cur_class_c'] = cur_class
                    job_pos = self._AT_CLIENT
                    cur_class = ctx.get('cur_class_c', cur_class)
                    ctx['cur_station'] = None

                # Route based on jobPos and layer type
                if job_pos == self._AT_CLIENT:
                    if host_stn is not None:
                        server_station = host_stn
                        if not is_cache_layer:
                            # HOST LAYER, NON-CACHE, atClient (MATLAB lines 1044-1096)
                            if is_next_prec_fork and fork_node is not None:
                                # FORK routing
                                P.set(cur_class, cur_class, client_delay, fork_node, 1.0)
                                post_and_succs = [s for s in nextaidxs if s in is_post_and_act]
                                f_idx = post_and_succs.index(nextaidx) + 1 if nextaidx in post_and_succs else -1
                                if f_idx > 0 and f_idx in fork_output_routers:
                                    fork_class_stack.append(cur_class)
                                    P.set(cur_class, cur_class, fork_node, fork_output_routers[f_idx], 1.0)
                                    P.set(cur_class, act_cls, fork_output_routers[f_idx], server_station, 1.0)
                                else:
                                    P.set(cur_class, act_cls, client_delay, server_station, graph[aidx, nextaidx])
                            elif aidx in is_pre_and_act and join_node is not None:
                                # JOIN routing
                                fork_class = fork_class_stack.pop()
                                P.set(cur_class, fork_class, client_delay, join_node, 1.0)
                                P.set(fork_class, act_cls, join_node, server_station, 1.0)
                            else:
                                # Serial routing
                                P.set(cur_class, act_cls, client_delay, server_station, graph[aidx, nextaidx])
                            # Set service at server (servtproc holds Distribution, not float)
                            if nextaidx < len(self.servtproc) and self.servtproc[nextaidx] is not None:
                                server_station.set_service(act_cls, self.servtproc[nextaidx])
                            job_pos = self._AT_SERVER
                            ctx['cur_station'] = server_station
                            cur_class = act_cls
                            # Record servt update map
                            if ctx.get('servt_map') is not None and ctx.get('idx') is not None:
                                ctx['servt_map'][ctx['idx']].append([ctx['idx'], nextaidx, 2, act_cls.get_index()])
                        else:
                            # CACHE LAYER, atClient (MATLAB lines 1097-1118)
                            P.set(cur_class, act_cls, client_delay, ctx.get('cache_node', server_station), graph[aidx, nextaidx])
                            job_pos = self._AT_CACHE
                            cur_class = act_cls
                    else:
                        # TASK LAYER, atClient (MATLAB lines 1119-1160)
                        if is_next_prec_fork and fork_node is not None:
                            # FORK routing
                            P.set(cur_class, cur_class, client_delay, fork_node, 1.0)
                            post_and_succs = [s for s in nextaidxs if s in is_post_and_act]
                            f_idx = post_and_succs.index(nextaidx) + 1 if nextaidx in post_and_succs else -1
                            if f_idx > 0 and f_idx in fork_output_routers:
                                fork_class_stack.append(cur_class)
                                P.set(cur_class, cur_class, fork_node, fork_output_routers[f_idx], 1.0)
                                P.set(cur_class, act_cls, fork_output_routers[f_idx], client_delay, 1.0)
                            else:
                                P.set(cur_class, act_cls, client_delay, client_delay, graph[aidx, nextaidx])
                        elif aidx in is_pre_and_act and join_node is not None:
                            # JOIN routing
                            fork_class = fork_class_stack.pop()
                            P.set(cur_class, fork_class, client_delay, join_node, 1.0)
                            P.set(fork_class, act_cls, join_node, client_delay, 1.0)
                        else:
                            # Serial routing
                            P.set(cur_class, act_cls, client_delay, client_delay, graph[aidx, nextaidx])
                        # Set service at client
                        if nextaidx in self.servtproc and self.servtproc[nextaidx] is not None:
                            client_delay.set_service(act_cls, self.servtproc[nextaidx])
                        job_pos = self._AT_CLIENT
                        cur_class = act_cls
                        # Record thinkt update map
                        if ctx.get('thinkt_map') is not None and ctx.get('idx') is not None:
                            ctx['thinkt_map'][ctx['idx']].append([ctx['idx'], nextaidx, 1, act_cls.get_index()])

                elif job_pos == self._AT_SERVER or job_pos == self._AT_CACHE:
                    if host_stn is not None:
                        from_stn = ctx.get('cur_station') or server_station
                        server_station = host_stn
                        if not is_cache_layer:
                            # HOST LAYER, NON-CACHE, atServer (MATLAB lines 1217-1258)
                            if is_next_prec_fork and fork_node is not None:
                                # FORK routing
                                P.set(cur_class, cur_class, from_stn, fork_node, 1.0)
                                post_and_succs = [s for s in nextaidxs if s in is_post_and_act]
                                f_idx = post_and_succs.index(nextaidx) + 1 if nextaidx in post_and_succs else -1
                                if f_idx > 0 and f_idx in fork_output_routers:
                                    fork_class_stack.append(cur_class)
                                    P.set(cur_class, cur_class, fork_node, fork_output_routers[f_idx], 1.0)
                                    P.set(cur_class, act_cls, fork_output_routers[f_idx], server_station, 1.0)
                                else:
                                    P.set(cur_class, act_cls, from_stn, server_station, graph[aidx, nextaidx])
                            elif aidx in is_pre_and_act and join_node is not None:
                                # JOIN routing
                                fork_class = fork_class_stack.pop()
                                P.set(cur_class, fork_class, from_stn, join_node, 1.0)
                                P.set(fork_class, act_cls, join_node, server_station, 1.0)
                            else:
                                # Serial routing
                                P.set(cur_class, act_cls, from_stn, server_station, graph[aidx, nextaidx])
                            # Set service at server (servtproc holds Distribution, not float)
                            if nextaidx < len(self.servtproc) and self.servtproc[nextaidx] is not None:
                                server_station.set_service(act_cls, self.servtproc[nextaidx])
                            job_pos = self._AT_SERVER
                            ctx['cur_station'] = server_station
                            cur_class = act_cls
                            if ctx.get('servt_map') is not None and ctx.get('idx') is not None:
                                ctx['servt_map'][ctx['idx']].append([ctx['idx'], nextaidx, 2, act_cls.get_index()])
                        else:
                            # CACHE LAYER, atServer/atCache (MATLAB lines 1163-1216)
                            cache_node = ctx.get('cache_node', server_station)
                            source_node = cache_node if job_pos == self._AT_CACHE else server_station
                            if is_next_prec_fork and fork_node is not None:
                                P.set(cur_class, cur_class, source_node, fork_node, 1.0)
                                post_and_succs = [s for s in nextaidxs if s in is_post_and_act]
                                f_idx = post_and_succs.index(nextaidx) + 1 if nextaidx in post_and_succs else -1
                                if f_idx > 0 and f_idx in fork_output_routers:
                                    fork_class_stack.append(cur_class)
                                    P.set(cur_class, cur_class, fork_node, fork_output_routers[f_idx], 1.0)
                                    P.set(cur_class, act_cls, fork_output_routers[f_idx], server_station, 1.0)
                                else:
                                    P.set(cur_class, act_cls, source_node, server_station, graph[aidx, nextaidx])
                            elif aidx in is_pre_and_act and join_node is not None:
                                fork_class = fork_class_stack.pop()
                                P.set(cur_class, fork_class, source_node, join_node, 1.0)
                                P.set(fork_class, act_cls, join_node, server_station, 1.0)
                            else:
                                P.set(cur_class, act_cls, source_node, server_station, graph[aidx, nextaidx])
                            if nextaidx < len(self.servtproc) and self.servtproc[nextaidx] is not None:
                                server_station.set_service(act_cls, self.servtproc[nextaidx])
                            job_pos = self._AT_SERVER
                            cur_class = act_cls
                            if ctx.get('servt_map') is not None and ctx.get('idx') is not None:
                                ctx['servt_map'][ctx['idx']].append([ctx['idx'], nextaidx, 2, act_cls.get_index()])
                    else:
                        # TASK LAYER, atServer (MATLAB lines 1276-1313)
                        if is_next_prec_fork and fork_node is not None:
                            # FORK routing
                            P.set(cur_class, cur_class, server_station, fork_node, 1.0)
                            post_and_succs = [s for s in nextaidxs if s in is_post_and_act]
                            f_idx = post_and_succs.index(nextaidx) + 1 if nextaidx in post_and_succs else -1
                            if f_idx > 0 and f_idx in fork_output_routers:
                                fork_class_stack.append(cur_class)
                                P.set(cur_class, cur_class, fork_node, fork_output_routers[f_idx], 1.0)
                                P.set(cur_class, act_cls, fork_output_routers[f_idx], client_delay, 1.0)
                            else:
                                P.set(cur_class, act_cls, server_station, client_delay, graph[aidx, nextaidx])
                        elif aidx in is_pre_and_act and join_node is not None:
                            # JOIN routing
                            fork_class = fork_class_stack.pop()
                            P.set(cur_class, fork_class, server_station, join_node, 1.0)
                            P.set(fork_class, act_cls, join_node, client_delay, 1.0)
                        else:
                            # Serial routing
                            P.set(cur_class, act_cls, server_station, client_delay, graph[aidx, nextaidx])
                        # Set service at client
                        if nextaidx in self.servtproc and self.servtproc[nextaidx] is not None:
                            client_delay.set_service(act_cls, self.servtproc[nextaidx])
                        job_pos = self._AT_CLIENT
                        ctx['cur_station'] = None
                        cur_class = act_cls
                        if ctx.get('thinkt_map') is not None and ctx.get('idx') is not None:
                            ctx['thinkt_map'][ctx['idx']].append([ctx['idx'], nextaidx, 1, act_cls.get_index()])

                # Recursive call (MATLAB lines 1316-1336)
                if aidx != nextaidx and not is_loop:
                    # cur_class_c is per-invocation in MATLAB and saved around the
                    # recursion in the JAR; a shared ctx entry leaks the callee's
                    # class back into the next fork branch
                    saved_cur_class_c = ctx.get('cur_class_c')
                    P, cur_class, job_pos = self._recur_act_graph(
                        P, tidx_caller, nextaidx, cur_class, job_pos, ctx)
                    ctx['cur_class_c'] = saved_cur_class_c
                    # Route back to task class (MATLAB lines 1322-1335)
                    task_cls = ctx['task_classes'][tidx_caller]
                    if job_pos == self._AT_CLIENT:
                        P.set(cur_class, task_cls, client_delay, client_delay, 1.0)
                    else:
                        # the job returns from the station it is actually at, which
                        # under flat is the callee's station, not this layer's server
                        _back = (ctx.get('cur_station') or server_station) \
                            if ctx.get('flat') else server_station
                        P.set(cur_class, task_cls, _back, client_delay, 1.0)
                    if not cur_class.name.endswith('.Aux'):
                        cur_class.completes = True

        return P, cur_class, job_pos

    def _setup_routing(self, layer_model: Network, idx: int = None, route_map: list = None):
        """Set up routing for a layer model with class switching (4-class model)."""
        stations = layer_model.get_nodes()
        classes = layer_model.get_classes()

        if len(stations) < 2 or len(classes) < 1:
            return

        client = None
        server = None
        cache_node = None
        for s in stations:
            if isinstance(s, Delay):
                client = s
            elif isinstance(s, Queue):
                if server is None:
                    server = s  # Use FIRST Queue as primary server (not last)
            elif isinstance(s, Cache):
                cache_node = s

        if client is None or server is None:
            return

        # Check if this is a cache layer
        is_cache_layer = layer_model.attribute.get('iscachelayer', False) if hasattr(layer_model, 'attribute') and layer_model.attribute else False

        P = layer_model.init_routing_matrix()

        # Separate classes by type
        task_classes = {}
        entry_classes = {}
        activity_classes = {}
        think_classes = {}  # Activity think-time classes
        call_classes = {}
        aux_classes = {}  # Aux classes for fractional call means

        for cls in classes:
            if hasattr(cls, 'attribute') and cls.attribute is not None:
                elem_type = cls.attribute[0] if len(cls.attribute) > 0 else 0
                elem_idx = cls.attribute[1] if len(cls.attribute) > 1 else 0
                if isinstance(elem_idx, np.integer):
                    elem_idx = int(elem_idx)

                if elem_type == LayeredNetworkElement.TASK:
                    task_classes[elem_idx] = cls
                elif elem_type == LayeredNetworkElement.ENTRY:
                    # an entry-arrival OpenClass carries the same [ENTRY, eidx]
                    # attribute as the entry's own closed class (as in MATLAB), so
                    # it must not displace it here: the task -> entry and
                    # entry -> activity routes below would then be wired to the
                    # open stream, which walks no activity graph
                    if not isinstance(cls, OpenClass):
                        entry_classes[elem_idx] = cls
                elif elem_type == LayeredNetworkElement.ACTIVITY:
                    # Separate Think classes from regular activity classes
                    if hasattr(cls, 'name') and cls.name.endswith('.Think'):
                        think_classes[elem_idx] = cls
                    else:
                        activity_classes[elem_idx] = cls
                elif elem_type == LayeredNetworkElement.CALL:
                    # Check if this is an Aux class (name ends with .Aux)
                    if hasattr(cls, 'name') and cls.name.endswith('.Aux'):
                        aux_classes[elem_idx] = cls
                    else:
                        call_classes[elem_idx] = cls

        # A routed group's members all resolve to the group's shared call class;
        # cls.attribute can only name one cidx, so the mapping is explicit.
        _cls_by_index = {c.get_index(): c for c in classes}
        for _cidx, _clsidx in (layer_model.attribute.get('call_class_of_cidx') or {}).items():
            if _clsidx in _cls_by_index:
                call_classes[_cidx] = _cls_by_index[_clsidx]

        # Build call_mean map from layer attribute
        call_mean_map = {}
        if 'aux_classes' in layer_model.attribute:
            for aux_info in layer_model.attribute['aux_classes']:
                if len(aux_info) >= 3:
                    aux_cls_idx, cidx, call_mean = aux_info[0], aux_info[1], aux_info[2]
                    call_mean_map[cidx] = call_mean

        # Check if this is a host layer (activities at server) or task layer (activities at client)
        # The 'ishost' attribute is set on the server station during layer construction
        is_host_layer = True  # Default to host layer
        if server is not None and hasattr(server, 'attribute'):
            is_host_layer = server.attribute.get('ishost', True)
        # Determine activity station: HOST layer = server, TASK layer = client
        act_station = server if is_host_layer else client

        # Set up class switching routing for each task
        for tidx, task_cls in task_classes.items():
            entries = self._get_entries_of_task(tidx)
            activities = self._get_activities_of_task(tidx)
            entry_cls_list = [entry_classes[eidx] for eidx in entries if eidx in entry_classes]
            activity_cls_list = [activity_classes[aidx] for aidx in activities if aidx in activity_classes]

            if not entry_cls_list and not activity_cls_list:
                # No entry or activity classes - simple routing
                P.set(task_cls, task_cls, client, server, 1.0)
                P.set(task_cls, task_cls, server, client, 1.0)
            elif entry_cls_list:
                # TASK -> ENTRY routing (at client with equal probability)
                ncaller_entries = len(entry_cls_list)
                for i, entry_cls in enumerate(entry_cls_list):
                    P.set(task_cls, entry_cls, client, client, 1.0 / ncaller_entries)
                    # routing probabilities among multiple entries updated from throughput ratios each iteration; mirrors MATLAB lines 392-394.
                    if ncaller_entries > 1 and idx is not None and route_map is not None:
                        eidx = entries[i]
                        # Format: [idx, tidx_caller, eidx, nodefrom, nodeto, classidxfrom, classidxto]
                        # nodefrom=nodeto=1 means client node (1-based index)
                        route_map[idx].append([idx, tidx, eidx, 1, 1, task_cls.get_index(), entry_cls.get_index()])

                # Check if this layer has Fork/Join nodes
                has_forkjoin = (layer_model.attribute.get('fork_node') is not None or
                                layer_model.attribute.get('join_node') is not None)

                if has_forkjoin:
                    # Use recursive activity graph traversal with Fork/Join routing
                    # (matches MATLAB recurActGraph in buildLayersRecursive.m)
                    ctx = {
                        'client_delay': client,
                        'server_station': server,
                        'is_host_layer': is_host_layer,
                        'is_cache_layer': is_cache_layer,
                        'fork_node': layer_model.attribute.get('fork_node'),
                        'fork_output_routers': layer_model.attribute.get('fork_output_routers', {}),
                        'join_node': layer_model.attribute.get('join_node'),
                        'fork_class_stack': [],
                        'activity_classes': activity_classes,
                        'task_classes': task_classes,
                        'call_classes': call_classes,
                        'aux_classes': aux_classes,
                        'think_classes': think_classes,
                        'call_mean_map': call_mean_map,
                        'is_post_and_act': layer_model.attribute.get('is_post_and_act', set()),
                        'is_pre_and_act': layer_model.attribute.get('is_pre_and_act', set()),
                        'job_pos_key': {},
                        'cur_class_key': {},
                        'cur_station_key': {},
                        'cur_station': None,
                        # flat layering resolves the station per element rather than
                        # using the layer's single server
                        'flat': bool(layer_model.attribute.get('flat')),
                        'srv_stations': layer_model.attribute.get('srv_stations', {}),
                        'cache_node': cache_node,
                        'servt_map': None,  # Already set during class creation
                        'thinkt_map': None,
                        'idx': idx,
                        'layer_model': layer_model,
                    }
                    for eidx, entry_cls in zip(entries, entry_cls_list):
                        if eidx in entry_classes:
                            P, _, _ = self._recur_act_graph(
                                P, tidx, eidx, entry_cls, self._AT_CLIENT, ctx)
                else:
                    # Original flat routing (no Fork/Join needed)
                    # Under flat layering each processor and task owns a station, so
                    # placement is resolved per element instead of using the layer's
                    # single server.
                    _srv_map = layer_model.attribute.get('srv_stations', {}) \
                        if hasattr(layer_model, 'attribute') else {}
                    _flat_layer = bool(layer_model.attribute.get('flat')) \
                        if hasattr(layer_model, 'attribute') else False

                    def _act_station(act_cls, default_srv):
                        # Station of the processor the activity class runs on
                        if not _flat_layer:
                            return default_srv
                        _a = act_cls.attribute[1] if getattr(act_cls, 'attribute', None) is not None \
                            and len(act_cls.attribute) > 1 else None
                        if _a is None:
                            return default_srv
                        _h = self._get_parent(self._get_parent(int(_a)))
                        _st = _srv_map.get(int(_h)) if _h is not None else None
                        return _st[0] if _st else None

                    # ENTRY -> ACTIVITY routing (for each entry, route to its bound activities)
                    for eidx, entry_cls in zip(entries, entry_cls_list):
                        if eidx in entry_classes:
                            bound_activities = self._get_activities_of_entry(eidx)
                            bound_act_cls_list = [activity_classes[aidx] for aidx in bound_activities if aidx in activity_classes]

                            is_cache_entry = False
                            if is_cache_layer and cache_node is not None:
                                if hasattr(self.lqn, 'itemproc') and isinstance(self.lqn.itemproc, dict):
                                    if eidx in self.lqn.itemproc and self.lqn.itemproc[eidx] is not None:
                                        is_cache_entry = True

                            if is_cache_entry and bound_act_cls_list:
                                first_act_cls = bound_act_cls_list[0]
                                P.set(entry_cls, first_act_cls, client, cache_node, 1.0)
                                cache_entry_aidx = bound_activities[0] if bound_activities else None
                                if cache_entry_aidx is not None and hasattr(self.lqn, 'graph') and self.lqn.graph is not None:
                                    successors = []
                                    for j in range(self.lqn.graph.shape[1]):
                                        if self.lqn.graph[cache_entry_aidx, j] != 0:
                                            successors.append(j)
                                    if len(successors) >= 2:
                                        hit_aidx = successors[0]
                                        miss_aidx = successors[1]
                                        hit_cls = activity_classes.get(hit_aidx)
                                        miss_cls = activity_classes.get(miss_aidx)
                                        if hit_cls is not None:
                                            P.set(first_act_cls, hit_cls, cache_node, server, 0.5)
                                        if miss_cls is not None:
                                            P.set(first_act_cls, miss_cls, cache_node, server, 0.5)
                            elif bound_act_cls_list:
                                first_act_cls = bound_act_cls_list[0]
                                _as = _act_station(first_act_cls, server) if (is_host_layer or _flat_layer) else None
                                if _as is not None:
                                    P.set(entry_cls, first_act_cls, client, _as, 1.0)
                                else:
                                    P.set(entry_cls, first_act_cls, client, client, 1.0)
                            elif activity_cls_list:
                                if is_host_layer:
                                    P.set(entry_cls, activity_cls_list[0], client, server, 1.0)
                                else:
                                    P.set(entry_cls, activity_cls_list[0], client, client, 1.0)
                            else:
                                P.set(entry_cls, task_cls, client, client, 1.0)

                if not has_forkjoin:
                    # For HOST layers, add explicit routing for activity classes from client to server
                    if (is_host_layer or _flat_layer) and activity_cls_list:
                        for act_cls in activity_cls_list:
                            _as = _act_station(act_cls, server)
                            if _as is not None:
                                P.set(act_cls, act_cls, client, _as, 1.0)

                    # Route through activities using flat loop (no fork/join)
                    for i, aidx in enumerate(activities):
                        if aidx not in activity_classes:
                            continue
                        act_cls = activity_classes[aidx]
                        if _flat_layer:
                            # the activity runs on its own processor's station
                            _as_i = _act_station(act_cls, server)
                            if _as_i is not None:
                                act_station = _as_i

                        sync_call_classes = []
                        _cgrp_by_cidx, _cgrp_members = self._call_groups_by_cidx()
                        _seen_groups = set()
                        if isinstance(self.lqn.callsof, dict):
                            calls = self.lqn.callsof.get(aidx, [])
                            for cidx in calls:
                                if cidx in call_classes:
                                    is_sync = True
                                    if hasattr(self.lqn, 'calltype') and self.lqn.calltype is not None:
                                        if isinstance(self.lqn.calltype, np.ndarray):
                                            calltype = self.lqn.calltype.flatten()[cidx] if cidx < len(self.lqn.calltype.flatten()) else CallType.SYNC
                                        elif isinstance(self.lqn.calltype, dict):
                                            calltype = self.lqn.calltype.get(cidx, CallType.SYNC)
                                        else:
                                            calltype = CallType.SYNC
                                        is_sync = (calltype == CallType.SYNC)
                                    if is_sync:
                                        # a group is one dispatch: take its shared class
                                        # once, at the position of its first member
                                        _g = _cgrp_by_cidx.get(cidx)
                                        if _g is not None:
                                            if _g[0] in _seen_groups:
                                                continue
                                            _seen_groups.add(_g[0])
                                        sync_call_classes.append(call_classes[cidx])

                        has_sync_call = len(sync_call_classes) > 0

                        if has_sync_call:
                            # Process each sync call individually, matching MATLAB routeSynchCall
                            # Each call checks its own target to determine server vs client routing
                            # under flat layering the activity sits on its own
                            # processor's station, so the job is at a server
                            job_at_client = (not is_host_layer) and not (
                                _flat_layer and act_station is not client)
                            cur_cls = act_cls

                            for ci, call_cls in enumerate(sync_call_classes):
                                call_cidx = call_cls.attribute[1] if hasattr(call_cls, 'attribute') else None
                                tgt_eidx_c = (self._get_call_target_entry(call_cidx)
                                              if call_cidx is not None else None)
                                tgt_tidx_c = (self._get_parent(tgt_eidx_c)
                                              if tgt_eidx_c is not None else None)

                                this_call_to_server = False
                                call_srv = server
                                if _flat_layer:
                                    # every called task has its own station here, so a
                                    # name match against the layer server never fires
                                    _st = _srv_map.get(int(tgt_tidx_c)) if tgt_tidx_c is not None else None
                                    if _st:
                                        call_srv = _st[0]
                                        this_call_to_server = True
                                elif tgt_tidx_c is not None and server is not None:
                                    tgt_name_c = self._get_hashname(tgt_tidx_c)
                                    this_call_to_server = (tgt_name_c == server.name)

                                call_mean = (call_mean_map.get(call_cidx, 1.0)
                                             if call_cidx is not None else 1.0)
                                nreplicas = 1
                                has_aux = call_cidx in aux_classes
                                aux_cls = aux_classes.get(call_cidx) if has_aux else None

                                # A routed group is ONE hop with n destinations, taken at a
                                # router whose only links are those destinations: the
                                # strategy routes over a NODE's links, not over one class's
                                # arcs, so any other node would let the job wander to
                                # stations the group never calls. The hop keeps the class
                                # (a state-dependent routing function is zero off the class
                                # diagonal); the switch is on the return arc. The 1/n split
                                # laid down here is the probabilistic reading a solver
                                # without state-dependent routing would see.
                                _grp = _cgrp_by_cidx.get(call_cidx) if call_cidx is not None else None
                                if _grp is not None and _flat_layer:
                                    _gid, _strategy = _grp
                                    _disp_cls, _ret_cls, _, _router = layer_model.attribute[
                                        'call_group_classes'][_gid]
                                    _tgt_stations = []
                                    for _mcidx in _cgrp_members[_gid]:
                                        _meidx = self._get_call_target_entry(_mcidx)
                                        _mtidx = self._get_parent(_meidx) if _meidx else None
                                        _st = _srv_map.get(int(_mtidx)) if _mtidx is not None else None
                                        if _st:
                                            _tgt_stations.append((_st[0], _mcidx))
                                    if len(_tgt_stations) >= 2:
                                        _from_node = client if job_at_client else act_station
                                        _share = 1.0 / len(_tgt_stations)
                                        P.set(cur_cls, _disp_cls, _from_node, _router, 1.0)
                                        for _st, _mcidx in _tgt_stations:
                                            P.set(_disp_cls, _disp_cls, _router, _st, _share)
                                            P.set(_disp_cls, _ret_cls, _st, client, 1.0)
                                        layer_model.attribute.setdefault('rrobin_sites', []).append(
                                            (_router.name, _disp_cls.get_index(), _strategy))
                                        cur_cls = _ret_cls
                                        job_at_client = True
                                        continue

                                if job_at_client:
                                    if this_call_to_server:
                                        # MATLAB: atClient, call to server entry
                                        if call_mean < 1:
                                            P.set(cur_cls, call_cls, client, call_srv, call_mean / nreplicas)
                                            P.set(call_cls, call_cls, call_srv, client, 1.0)
                                            if has_aux:
                                                P.set(cur_cls, aux_cls, client, client, 1 - call_mean)
                                                P.set(aux_cls, call_cls, client, client, 1.0)
                                            cur_cls = call_cls
                                        elif call_mean == 1:
                                            P.set(cur_cls, call_cls, client, call_srv, 1.0 / nreplicas)
                                            P.set(call_cls, call_cls, call_srv, client, 1.0)
                                            cur_cls = call_cls
                                        else:  # call_mean > 1
                                            P.set(cur_cls, call_cls, client, call_srv, 1.0 / nreplicas)
                                            if has_aux:
                                                P.set(call_cls, aux_cls, call_srv, client, 1.0)
                                                P.set(aux_cls, call_cls, client, call_srv, (1.0 - 1.0 / call_mean) / nreplicas)
                                                P.set(aux_cls, call_cls, client, client, 1.0 / call_mean)
                                                cur_cls = call_cls  # matches MATLAB line 783: curClass = cidxClass{cidx}
                                            else:
                                                cur_cls = call_cls
                                        job_at_client = True
                                    else:
                                        # MATLAB: atClient, call NOT to server
                                        if call_mean < 1:
                                            # Deterministic visit: the call mean is
                                            # embedded in the demand (callservt)
                                            P.set(cur_cls, call_cls, client, client, 1.0)
                                            if has_aux:
                                                P.set(call_cls, aux_cls, client, client, 1.0)
                                                cur_cls = aux_cls
                                            else:
                                                cur_cls = call_cls
                                        elif call_mean == 1:
                                            P.set(cur_cls, call_cls, client, client, 1.0)
                                            cur_cls = call_cls
                                        else:  # call_mean > 1
                                            P.set(cur_cls, call_cls, client, client, 1.0)
                                            if has_aux:
                                                P.set(call_cls, aux_cls, client, client, 1.0)
                                                cur_cls = aux_cls
                                            else:
                                                cur_cls = call_cls
                                        job_at_client = True
                                else:
                                    # job at server
                                    if this_call_to_server:
                                        # MATLAB: atServer, call to server entry
                                        _from_stn = act_station if _flat_layer else server
                                        if call_mean < 1:
                                            if _flat_layer:
                                                # the skip flow enters the Aux class and the
                                                # reply transits the client in the call class
                                                if has_aux:
                                                    P.set(cur_cls, aux_cls, _from_stn, client, 1 - call_mean)
                                                    P.set(aux_cls, call_cls, client, client, 1.0)
                                                else:
                                                    P.set(cur_cls, call_cls, _from_stn, client, 1 - call_mean)
                                                P.set(cur_cls, call_cls, _from_stn, call_srv, call_mean)
                                                P.set(call_cls, call_cls, call_srv, client, 1.0)
                                                client.set_service(call_cls, Immediate())
                                                # both the skip and the visit end in the call
                                                # class, so continuing from Aux would emit the
                                                # next call's arcs out of a class that has
                                                # already been routed away
                                                cur_cls = call_cls
                                            else:
                                                P.set(cur_cls, call_cls, server, client, 1 - call_mean)
                                                P.set(cur_cls, call_cls, server, server, call_mean)
                                                cur_cls = aux_cls if has_aux else call_cls
                                            job_at_client = True
                                        elif call_mean == 1:
                                            if _flat_layer:
                                                # the reply returns the job to the client
                                                P.set(cur_cls, call_cls, _from_stn, call_srv, 1.0)
                                                P.set(call_cls, call_cls, call_srv, client, 1.0)
                                                client.set_service(call_cls, Immediate())
                                                job_at_client = True
                                            else:
                                                P.set(cur_cls, call_cls, server, server, 1.0)
                                                job_at_client = False
                                            cur_cls = call_cls
                                        else:  # call_mean > 1
                                            if _flat_layer:
                                                # the geometric repeat visits the CALLED task's
                                                # station and transits the client between
                                                # visits, as the atClient split does: a
                                                # self-loop merges the visits into one and
                                                # under-counts the call's aggregate service
                                                P.set(cur_cls, call_cls, _from_stn, call_srv, 1.0)
                                                if has_aux:
                                                    P.set(call_cls, aux_cls, call_srv, client, 1.0)
                                                    P.set(aux_cls, call_cls, client, call_srv, 1 - 1.0 / call_mean)
                                                    P.set(aux_cls, call_cls, client, client, 1.0 / call_mean)
                                                    client.set_service(call_cls, Immediate())
                                                job_at_client = True
                                                cur_cls = call_cls
                                            else:
                                                P.set(cur_cls, call_cls, server, server, 1.0)
                                                if has_aux:
                                                    P.set(call_cls, call_cls, server, server, 1 - 1.0 / call_mean)
                                                    P.set(call_cls, aux_cls, server, client, 1.0 / call_mean)
                                                job_at_client = True
                                                cur_cls = aux_cls if has_aux else call_cls
                                    else:
                                        # atServer, call NOT to server
                                        # callmean not needed since we use ResidT to model service time at client
                                        P.set(cur_cls, call_cls, server, client, 1.0)
                                        if call_mean < 1:
                                            if has_aux:
                                                P.set(call_cls, aux_cls, client, client, 1.0)
                                                cur_cls = aux_cls
                                            else:
                                                cur_cls = call_cls
                                        elif call_mean == 1:
                                            cur_cls = call_cls
                                        else:  # call_mean > 1
                                            if has_aux:
                                                P.set(call_cls, aux_cls, client, client, 1.0)
                                                cur_cls = aux_cls
                                            else:
                                                cur_cls = call_cls
                                        job_at_client = True

                                # (forwarding handled via pseudo rendezvous calls)

                            # After all calls, route to successor activity or back to task
                            source_node = client if job_at_client else server
                            graph_successors = self._get_activity_successors(aidx, activity_classes)
                            if graph_successors:
                                for succ_aidx, prob in graph_successors:
                                    succ_act_cls = activity_classes[succ_aidx]
                                    P.set(cur_cls, succ_act_cls, source_node, act_station, prob)
                            else:
                                P.set(cur_cls, task_cls, source_node, client, 1.0)
                        else:
                            graph = self.lqn.graph
                            has_graph_successors = False
                            if isinstance(graph, np.ndarray) and aidx < graph.shape[0]:
                                successors = []
                                for succ_aidx in range(graph.shape[1]):
                                    if graph[aidx, succ_aidx] > 0:
                                        if succ_aidx in activity_classes:
                                            prob = graph[aidx, succ_aidx]
                                            successors.append((succ_aidx, prob))
                                if successors:
                                    has_graph_successors = True
                                    for succ_aidx, prob in successors:
                                        succ_act_cls = activity_classes[succ_aidx]
                                        P.set(act_cls, succ_act_cls, act_station, act_station, prob)
                            if not has_graph_successors:
                                is_terminal = False
                                if hasattr(self.lqn, 'replygraph') and self.lqn.replygraph is not None:
                                    act_local_idx = aidx - self.lqn.ashift
                                    if isinstance(self.lqn.replygraph, np.ndarray):
                                        if 0 <= act_local_idx < self.lqn.replygraph.shape[0]:
                                            if np.any(self.lqn.replygraph[act_local_idx, :] > 0):
                                                is_terminal = True
                                if is_terminal:
                                    # activity think time in series with host demand; see _kb/06-solver-catalog.md LN Activity think time section.
                                    if aidx in think_classes:
                                        tc = think_classes[aidx]
                                        P.set(act_cls, tc, act_station, client, 1.0)
                                        P.set(tc, task_cls, client, client, 1.0)
                                    else:
                                        P.set(act_cls, task_cls, act_station, client, 1.0)
                                elif i < len(activities) - 1:
                                    next_aidx = activities[i + 1]
                                    if next_aidx in activity_classes:
                                        next_act_cls = activity_classes[next_aidx]
                                        P.set(act_cls, next_act_cls, act_station, act_station, 1.0)
                                    else:
                                        P.set(act_cls, task_cls, act_station, client, 1.0)
                                else:
                                    P.set(act_cls, task_cls, act_station, client, 1.0)
            else:
                # No entry classes but have activity classes
                # TASK -> first ACTIVITY
                first_act_cls = activity_cls_list[0]
                P.set(task_cls, first_act_cls, client, act_station, 1.0)

                # Route through activities, last returns to task at client
                for i, act_cls in enumerate(activity_cls_list):
                    if i < len(activity_cls_list) - 1:
                        next_act_cls = activity_cls_list[i + 1]
                        P.set(act_cls, next_act_cls, act_station, act_station, 1.0)
                    else:
                        P.set(act_cls, task_cls, act_station, client, 1.0)

        # entry-level open-arrival routing Source->server->Sink; mirrors JAR SolverLN.java:881-890 (replication handled separately).
        source_station = layer_model.attribute.get('source_station')
        sink_station = layer_model.attribute.get('sink_station')
        entry_open_classes = layer_model.attribute.get('entry_open_classes', [])
        if source_station is not None and sink_station is not None and entry_open_classes:
            nreplicas_for_open = layer_model.attribute.get('nreplicas', 1) or 1
            for open_cls, _eidx in entry_open_classes:
                # buildLayersRecursive.m:480-485 clears the class first: an entry
                # arrival walks no activity graph, and a leftover class-switch arc
                # puts the open class in the closed chain, making it mixed.
                P.remove_job_class(open_cls)
                P.set(open_cls, open_cls, source_station, server, 1.0 / float(nreplicas_for_open))
                P.set(open_cls, open_cls, server, sink_station, 1.0)

        # async-call-injection open class recirculates at server until calls done, drains to sink; mirrors MATLAB buildLayersRecursive.m:415-432/JAR:856-878.
        async_open_classes = layer_model.attribute.get('async_open_classes', [])
        if source_station is not None and sink_station is not None and async_open_classes:
            nreplicas_for_open = layer_model.attribute.get('nreplicas', 1) or 1
            for open_cls, _cidx, call_mean in async_open_classes:
                try:
                    cm = float(call_mean) if call_mean is not None else 1.0
                except (TypeError, ValueError):
                    cm = 1.0
                if cm <= 0 or not np.isfinite(cm):
                    cm = 1.0
                p_drain = 1.0 / cm
                if cm < 1:
                    # fewer than one call per arrival: a single Bernoulli pass, the
                    # geometric loop below would need a negative repeat probability
                    P.set(open_cls, open_cls, source_station, sink_station, 1.0 - cm)
                    P.set(open_cls, open_cls, source_station, server, cm / float(nreplicas_for_open))
                    P.set(open_cls, open_cls, server, sink_station, 1.0)
                else:
                    P.set(open_cls, open_cls, source_station, server, 1.0 / float(nreplicas_for_open))
                    # Server self-loop for recirculation (primary only; the replication
                    # block below mirrors self-loops to each replica).
                    if cm != 1:
                        P.set(open_cls, open_cls, server, server,
                              (1.0 - p_drain) / float(nreplicas_for_open))
                    P.set(open_cls, open_cls, server, sink_station, p_drain)

        # replicate routing to additional server replicas (nreplicas>1), splitting incoming probabilities; mirrors MATLAB serverStation loops.
        nreplicas = layer_model.attribute.get('nreplicas', 1)
        if nreplicas > 1:
            all_server_stations = layer_model.attribute.get('server_stations', [])
            if len(all_server_stations) > 1:
                primary = all_server_stations[0]
                replicas = all_server_stations[1:]

                # Copy service distributions from primary to each replica
                for replica in replicas:
                    for jc, dist in primary._service_process.items():
                        replica.set_service(jc, dist)
                    # Copy delay-off if present
                    if hasattr(primary, '_setup_time') and primary._setup_time:
                        if not hasattr(replica, '_setup_time') or replica._setup_time is None:
                            replica._setup_time = {}
                        replica._setup_time.update(primary._setup_time)
                    if hasattr(primary, '_delay_off_time') and primary._delay_off_time:
                        if not hasattr(replica, '_delay_off_time') or replica._delay_off_time is None:
                            replica._delay_off_time = {}
                        replica._delay_off_time.update(primary._delay_off_time)

                # Replicate routing entries: for each P entry involving the primary,
                # create entries for each replica with adjusted probabilities
                new_entries = []
                remove_entries = []
                for (cs, cd), route_dict in P._routes.items():
                    for (ns, nd), prob in list(route_dict.items()):
                        if nd == primary and ns != primary:
                            # Incoming to primary from another node: split across replicas
                            # Primary gets prob/nreplicas, each replica gets prob/nreplicas
                            new_prob = prob / nreplicas
                            remove_entries.append((cs, cd, ns, nd))
                            new_entries.append((cs, cd, ns, primary, new_prob))
                            for replica in replicas:
                                new_entries.append((cs, cd, ns, replica, new_prob))
                        elif ns == primary and nd != primary:
                            # Outgoing from primary to another node: same for each replica
                            for replica in replicas:
                                new_entries.append((cs, cd, replica, nd, prob))
                        elif ns == primary and nd == primary:
                            # Self-loop on primary: replicate as self-loop on each replica
                            for replica in replicas:
                                new_entries.append((cs, cd, replica, replica, prob))

                # Apply changes
                for cs, cd, ns, nd in remove_entries:
                    key = (cs, cd)
                    if key in P._routes and (ns, nd) in P._routes[key]:
                        del P._routes[key][(ns, nd)]
                for cs, cd, ns, nd, prob in new_entries:
                    P.set(cs, cd, ns, nd, prob)

        # deferred delayed-hit retrieval wiring applied here (dict-based RoutingMatrix needs no P growth); see _kb/09-ldes-and-cache.md.
        rw = layer_model.attribute.get('retrieval_wiring') if isinstance(layer_model.attribute, dict) else None
        if rw is not None:
            cache_node = rw['cache_node']
            fetch = rw['fetch']
            read_class = rw['read_class']
            miss_class = rw['miss_class']
            miss_aidx = rw['miss_aidx']
            # Fetch service = the miss activity's full service (host demand + backend call).
            svc = None
            if self.servtproc is not None and miss_aidx < len(self.servtproc):
                svc = self.servtproc[miss_aidx]
            if svc is None:
                svc = Exp(1.0)
            fetch.set_service(read_class, svc)
            P.set(read_class, read_class, cache_node, fetch, 1.0)
            P.set(read_class, read_class, fetch, cache_node, 1.0)
            cache_node.set_retrieval_system(read_class, miss_class, fetch)
            # Tag the auto-generated retrieval classes as non-completing (they map to no
            # LQN activity and are skipped by the LN class->activity updmaps).
            for rcls in getattr(cache_node, '_retrieval_classes', {}).values():
                if rcls is not None and hasattr(rcls, 'completes'):
                    rcls.completes = False
            # Pad every service station's service to the full class count with Disabled:
            # the retrieval classes are served only at the fetch station.
            for st in layer_model.get_nodes():
                if isinstance(st, Queue):
                    for jc in layer_model.classes:
                        if st.get_service(jc) is None:
                            st.set_service(jc, Disabled())

        if layer_model.attribute.get('flat'):
            # link() installs RAND routing for every (node, class) pair left
            # without an outgoing arc. With one station per server in a single
            # layer those spurious uniform arcs let a class wander to stations
            # it never visits, trapping the flow in a sub-cycle and leaving the
            # reference class with zero visits (MATLAB buildLayersRecursive.m).
            from ...constants import RoutingStrategy as _RS
            from ...lang.nodes import Sink as _Sink
            outflow = {}
            for (class_src, _class_dst), routes in P._routes.items():
                for (node_src, _node_dst), prob in routes.items():
                    if prob > 0:
                        outflow.setdefault(id(node_src), set()).add(id(class_src))
            for node in layer_model.get_nodes():
                if isinstance(node, _Sink):
                    continue
                here = outflow.get(id(node), set())
                for jobclass in layer_model.classes:
                    if id(jobclass) not in here:
                        node.setRouting(jobclass, _RS.DISABLED)

        layer_model.link(P)
        # link() installs the probabilistic split; the declared strategy replaces
        # it on the dispatch (node, class), whose only arcs are the group's targets
        for _node_name, _disp_idx, _strategy in layer_model.attribute.get('rrobin_sites', []):
            for _n in layer_model.get_nodes():
                if _n.name == _node_name:
                    for _c in layer_model.classes:
                        if _c.get_index() == _disp_idx:
                            _n.setRouting(_c, _strategy)
                    break
        self._add_layer_admission_constraint(layer_model, idx)
        self._add_layer_rate_dependence(layer_model)

    def _add_layer_rate_dependence(self, layer_model):
        """
        Emit the service-rate dependences declared on the server elements of this
        layer onto their stations -- see _kb/06-solver-catalog.md (LN section).
        Load dependence reads the total station population and maps directly; the
        class- and joint-dependent handles are declared over the server's
        operands, which the layer represents as job classes, so each operand is
        expanded onto the classes that occupy the server on its behalf.
        """
        lqn = self.lqn
        lld = getattr(lqn, 'lldscaling', None) or {}
        cd = getattr(lqn, 'cdscaling', None) or {}
        jd = getattr(lqn, 'jdscaling', None) or {}
        pools = getattr(lqn, 'pools', None) or {}
        if not lld and not cd and not jd and not pools:
            return
        srv_stations = layer_model.attribute.get('srv_stations', {})
        nclasses = len(layer_model.classes)
        for sidx, stations in srv_stations.items():
            has_ld = sidx in lld
            has_cd = sidx in cd
            has_jd = sidx in jd
            has_pools = sidx in pools
            if not (has_ld or has_cd or has_jd or has_pools):
                continue
            cols = self._layer_operand_classes(layer_model, sidx)
            one_class_per_operand = all(len(c) <= 1 for c in cols)
            for ss in stations:
                if has_ld:
                    ss.set_load_dependence(lld[sidx])
                if has_cd:
                    # beta_{i,r} is product-form only while an operand maps to a single
                    # class; where it aggregates several, the same scaling is emitted as
                    # a joint dependence, which is numerically identical but not exact
                    handle = _layer_dep_handle(cd[sidx], cols, nclasses, layer_model)
                    peak = _layer_peak(lqn.cdscalingpeak[sidx], cols, nclasses)
                    if one_class_per_operand:
                        ss.set_class_dependence(handle, peak)
                    else:
                        ss.set_joint_dependence(handle, peak)
                if has_jd:
                    ss.set_joint_dependence(_layer_dep_handle(jd[sidx], cols, nclasses, layer_model),
                                            _layer_peak(lqn.jdscalingpeak[sidx], cols, nclasses))
                if has_pools:
                    # A compatibility declaration IS a rate law: the pools clear
                    # mu(n) of sn_compat_rate, which reads only the SUPPORT of n
                    # and is therefore order independent. Normalising by the
                    # every-pool-active peak makes eta(n) <= 1 with equality at
                    # full support, so a fully-compatible pool reproduces the
                    # plain multiplicity station exactly. The lowering is to a
                    # JOINT dependence, hence an approximation in the layer: see
                    # _kb/06-solver-catalog.md (LN section) for why the exact OI
                    # analyzer cannot serve a class-switching layer.
                    pl = pools[sidx]

                    def _eta_pool(nop, _pl=pl):
                        return sn_compat_scaling(_pl['compat'], _pl['counts'], _pl['rates'], nop)

                    ss.set_joint_dependence(
                        _layer_dep_handle(_eta_pool, cols, nclasses, layer_model),
                        _layer_peak(np.ones(len(cols)), cols, nclasses))

    def _layer_operand_classes(self, layer_model, sidx):
        """
        Layer classes (1-based) through which each operand of server SIDX occupies
        its station: the tasks of a host through the classes of their activities,
        the entries of a task through the classes of the calls that target them.
        """
        lqn = self.lqn
        if sidx <= lqn.nhosts:
            operand_idx = lqn.tasksof.get(sidx, [])
            cls_by_elem = {}
            for cls_index, aidx in layer_model.attribute.get('activities', []):
                cls_by_elem.setdefault(aidx, []).append(cls_index)
            return [[c for a in lqn.actsof.get(j, []) for c in cls_by_elem.get(a, [])]
                    for j in operand_idx]
        operand_idx = lqn.entriesof.get(sidx, [])
        cls_by_elem = {}
        for row in layer_model.attribute.get('calls', []):
            cls_by_elem.setdefault(row[3], []).append(row[0])
        return [list(cls_by_elem.get(j, [])) for j in operand_idx]

    def _add_layer_admission_constraint(self, layer_model, idx):
        """
        Emit the admission constraint of host or task IDX as a finite capacity
        region on that layer's server station -- see _kb/06-solver-catalog.md
        (LN section). The constraint is declared over entries or tasks, which the
        layer represents as job classes, so each declared column is expanded onto
        the classes that occupy the server on its behalf.
        """
        lqn = self.lqn
        lincon = getattr(lqn, 'lincon', None)
        if idx is None or not lincon or idx not in lincon:
            return
        a_elem, b_elem = lincon[idx]
        if a_elem is None or a_elem.size == 0:
            return
        is_host_layer = layer_model.attribute.get('ishost', False)
        nclasses = len(layer_model.classes)
        a_layer = np.zeros((a_elem.shape[0], nclasses))
        if is_host_layer:
            # column j is task constrained_idx[j], occupying the host through its activities
            constrained_idx = lqn.tasksof.get(idx, [])
            cls_by_elem = {}
            for cls_index, aidx in layer_model.attribute.get('activities', []):
                cls_by_elem.setdefault(aidx, []).append(cls_index)
            members = {j: [c for a in lqn.actsof.get(constrained_idx[j], [])
                           for c in cls_by_elem.get(a, [])]
                       for j in range(len(constrained_idx))}
        else:
            # column j is entry constrained_idx[j], occupied by the calls that target it
            constrained_idx = lqn.entriesof.get(idx, [])
            cls_by_elem = {}
            for row in layer_model.attribute.get('calls', []):
                cls_by_elem.setdefault(row[3], []).append(row[0])
            members = {j: list(cls_by_elem.get(constrained_idx[j], []))
                       for j in range(len(constrained_idx))}
        for j, cls_indices in members.items():
            for cls_index in cls_indices:
                # the attribute maps carry 1-based class indices
                if 1 <= cls_index <= nclasses:
                    a_layer[:, cls_index - 1] += a_elem[:, j]
        if not np.any(a_layer):
            return
        stations = layer_model.attribute.get('server_stations', [])
        if not stations:
            return
        # One region spanning every replica: the constraint models a passive
        # resource of the server as a whole (a semaphore, a connection pool), so
        # replicas share the tokens rather than each holding a private copy
        region = layer_model.add_region(stations[0], *stations[1:])
        region.setConstraint(a_layer, b_elem)

    def _region_wait(self, layer_idx, nodeidx_0, classidx_0, result):
        """
        Waiting time absorbed by an admission constraint in a layer, recovered by
        Little's law from the layer population deficit. A job blocked at the
        constraint is counted at no station (JMT WAITQ convention), so its wait is
        absent from RN; without this the caller never sees the blocking and the
        fixed point loses flow balance.

        The population is conserved per chain, not per class: a job in a layer
        switches class along the activity graph, so the call class itself carries
        population 0. Splitting the chain deficit by throughput gives every
        region-visiting class the same wait.
        """
        flags = getattr(self, 'layer_has_region', None)
        if not flags or layer_idx < 0 or layer_idx >= len(flags) or not flags[layer_idx]:
            return 0.0
        chains = self.layer_chains[layer_idx]
        if chains is None or chains.size == 0:
            return 0.0
        rows = np.flatnonzero(chains[:, classidx_0])
        if rows.size == 0:
            return 0.0
        chain_classes = np.flatnonzero(chains[rows[0], :])
        QN = result.get('QN')
        TN = result.get('TN')
        if QN is None or TN is None:
            return 0.0
        layer = self.ensemble[layer_idx]
        chain_pop = 0.0
        for k in chain_classes:
            pop = getattr(layer.classes[k], 'population', None)
            if pop is not None and np.isfinite(pop):
                chain_pop += pop
        deficit = chain_pop - float(np.sum(QN[:, chain_classes]))
        xregion = float(np.sum(TN[nodeidx_0, chain_classes]))
        if np.isfinite(deficit) and deficit > 0 and xregion > GlobalConstants.FineTol:
            return deficit / xregion
        return 0.0

    def init(self):
        """Initialize before starting iterations (matches MATLAB init)."""
        # The moment3 pass is terminal WITHIN ONE SOLVE, so the flag is scoped to
        # one iterate(): left standing, the terminal test in converged() fires at
        # it=0 on the NEXT solve, the loop body never runs and every metric comes
        # back zero. See BUGS.md BUG-97.
        self.moment_pass_done = False
        line_debug("LN init: %d layers, relaxation=%s (omega=%.3f)",
                   self.nlayers,
                   self.options.config.get('relax', 'none'),
                   getattr(self, 'relax_omega', 1.0))
        lqn = self.lqn

        self.unique_route_prob_updmap = np.unique(self.route_prob_updmap[:, 0]) if len(self.route_prob_updmap) > 0 else np.array([])

        self.tput = np.zeros(lqn.nidx)
        self.tputproc = [None] * lqn.nidx
        self.util = np.zeros(lqn.nidx)
        self.servt = np.zeros(lqn.nidx)
        self.residt = np.zeros(lqn.nidx)
        self.thinkt = np.zeros(lqn.nidx)
        self.thinktproc = [None] * lqn.nidx
        self.callservt = np.zeros(lqn.ncalls)
        self.callresidt = np.zeros(lqn.ncalls)
        self.servtmatrix = self._get_entry_service_matrix()

        # feature-set gate stays armed on layer solvers; see _kb/06-solver-catalog.md LN Feature checks stay armed.
        for e in range(self.nlayers):
            if self.solvers[e] is not None:
                if hasattr(self.solvers[e], 'enable_checks'):
                    self.solvers[e].enable_checks = True

        # Initialize relaxation state
        relax_mode = self.options.config.get('relax', 'none')
        if relax_mode == 'auto':
            self.relax_omega = 1.0
        elif relax_mode in ['fixed', 'adaptive']:
            self.relax_omega = self.options.config.get('relax_factor', 0.9)
        else:
            self.relax_omega = 1.0

        self.relax_err_history = []
        self.servt_prev = np.full(lqn.nidx, np.nan)
        self.residt_prev = np.full(lqn.nidx, np.nan)
        self.tput_prev = np.full(lqn.nidx, np.nan)
        self.thinkt_prev = np.full(lqn.nidx, np.nan)
        self.callservt_prev = np.full(lqn.ncalls, np.nan)
        self.callresidt_prev = np.full(lqn.ncalls, np.nan)

        # stochastic iteration mode resolution; see _kb/06-solver-catalog.md LN Convergence test: stochastic iteration dispatch.
        self.stochlayers = np.zeros(self.nlayers, dtype=bool)
        for e in range(self.nlayers):
            solver = self.solvers[e]
            if solver is not None and hasattr(solver, 'isStochastic'):
                self.stochlayers[e] = solver.isStochastic()
        stoch_mode = str(self.options.config.get('stochiter', 'auto')).lower()
        self.stochiter_auto = (stoch_mode == 'auto')
        if self.stochiter_auto:
            stoch_mode = 'rm' if np.any(self.stochlayers) else 'off'
        self.stochiter_mode = stoch_mode
        self.stochiter_start = None
        self.stoch_avg = [None] * self.nlayers
        self.stoch_avg_count = 0
        self.stoch_servt_avg = None
        self.stoch_residt_avg = None
        # rm/crn layer seed base: user seed if given, else randomized; matches MATLAB/JAR default.
        user_seed = getattr(self.options, 'seed', None)
        if user_seed is not None:
            self.stochiter_seed_base = int(user_seed)
        else:
            self.stochiter_seed_base = int(np.random.randint(1, 10**6))
        line_debug("LN init: stochastic iteration mode=%s (%d stochastic layers)",
                   self.stochiter_mode, int(np.sum(self.stochlayers)))

        # AND-fork visit correction applied post struct-build; MVA solver reset afterward so it recomputes demands from corrected visits.
        for e in range(self.nlayers):
            if e < len(self.ensemble) and self.ensemble[e] is not None:
                layer = self.ensemble[e]
                if layer.attribute.get('has_fork', False):
                    # Force struct build if not already done
                    if not layer._has_struct:
                        layer.refresh_struct()
                    self._apply_fork_visit_correction(layer)
                    # Reset solver to pick up corrected visits
                    if e < len(self.solvers) and self.solvers[e] is not None:
                        if hasattr(self.solvers[e], 'reset'):
                            self.solvers[e].reset()


        # layer_init='bound' seeds layer throughput from Majumdar-Woodside box bounds; see _kb/06-solver-catalog.md LN Feature checks stay armed section.
        init_mode = self.options.config.get('layer_init', None) \
            if isinstance(self.options.config, dict) \
            else getattr(self.options.config, 'layer_init', None)
        if init_mode and str(init_mode).lower() in ('bound', 'boxbound', 'mwba'):
            try:
                tup, _ = self._box_bounds(True)
                tlo, _ = self._box_bounds(False)
                for idx in range(lqn.nidx):
                    u = tup[idx]
                    l = tlo[idx]
                    if np.isfinite(u) and np.isfinite(l) and u > 0 and l > 0:
                        x = np.sqrt(u * l)
                    elif np.isfinite(u):
                        x = u
                    elif np.isfinite(l):
                        x = l
                    else:
                        x = 0.0
                    if x > 0:
                        self.tput[idx] = x
                        self.tputproc[idx] = Exp.fit_rate(x)
                line_debug("LN init: throughputs initialized from robust box bounds")
            except Exception as ex:
                line_debug("LN box-bound initialization skipped: %s", str(ex))

    def _layer_index_of(self, elem_idx: int) -> Optional[int]:
        """0-based ensemble index of the layer where ELEM_IDX is a server, None if there is none."""
        if self.idxhash is None or elem_idx >= len(self.idxhash):
            return None
        e = self.idxhash[elem_idx]
        if e is None or (isinstance(e, float) and np.isnan(e)):
            return None
        e = int(e)
        return e if 0 <= e < len(self.ensemble) else None

    def _layer_takes_interlock(self, e: int) -> bool:
        """True when the layer solver applies Eq. (4.7) inside its own MVA.

        Only the MVA layer solver reads options.config['interlock'], and only a layer whose
        sole queueing stations are the host's own tasks can take a matrix built for that host:
        under flat layering one layer holds every server, so the correction stays on the
        residence times there.
        """
        from ..solver_mva import SolverMVA as _SolverMVA
        if e >= len(self.solvers) or not isinstance(self.solvers[e], _SolverMVA):
            return False
        cfg = getattr(self.options, 'config', None)
        layering = None
        if isinstance(cfg, dict):
            layering = cfg.get('layering', None)
        elif cfg is not None:
            layering = getattr(cfg, 'layering', None)
        if isinstance(layering, str) and layering.lower() in ('flat', 'squashed'):
            return False
        # A layer whose MVA path has no interlock term would be moved to another algorithm by
        # the matrix alone: exact multiserver MVA would become AMVA, the linearizer would
        # become the load-dependent forward step. That swap is worth far more than the
        # correction it carries, and on a layer sitting near a bifurcation it turns the LN
        # iteration into a limit cycle. Such a layer keeps the residt scaling instead.
        from ...api.solvers.mva.analyzers import mva_carries_interlock
        return mva_carries_interlock(self.ensemble[e].getStruct(), self.solvers[e].options)

    def _build_layer_interlock(self, e: int, host_tasks, task_pr_il, task_PrIL):
        """Class-level interlock matrix of one host layer.

        IL[r,s] is the share of the class-s queue that a class-r arrival must not see at the
        host. The matrix is CLASS-indexed, not chain-indexed, so that a later refreshChains
        cannot leave it stale; the layer solver aggregates it to chains against the struct it
        is about to solve. Two classes are interlocked only if BOTH their tasks are, which is
        the 0/1 relation ir_mkj of Eq. (5); the diagonal stays zero, since a request always
        sees its own class in full. The entry is the Eq. (5) product Pr(IL_ms)*IR_ms*IR_mr,
        asymmetric in (r,s) because Pr(IL) is taken from the QUEUED class s, so that the
        layer's ILw(r,s) = 1-IL(r,s) is the lower-level adjustment rate r_lower.
        """
        classes = self.ensemble[e].classes
        nclasses = len(classes)
        class_pr_il = np.zeros(nclasses)   # IR
        class_PrIL = np.zeros(nclasses)    # Pr(IL)
        host_tasks = list(host_tasks)
        for r in range(nclasses):
            tidx = self._client_task_of_class(e, r)
            if tidx is None:
                continue
            if tidx in host_tasks:
                ti = host_tasks.index(tidx)
                class_pr_il[r] = task_pr_il[ti]
                class_PrIL[r] = task_PrIL[ti]
        IL = np.zeros((nclasses, nclasses))
        for r in range(nclasses):
            if class_pr_il[r] <= GlobalConstants.FineTol:
                continue
            for sIl in range(nclasses):
                if sIl == r or class_pr_il[sIl] <= GlobalConstants.FineTol:
                    continue
                IL[r, sIl] = class_PrIL[sIl] * class_pr_il[sIl] * class_pr_il[r]
        return IL if np.any(IL > GlobalConstants.FineTol) else None

    def _client_task_of_class(self, e: int, c: int) -> Optional[int]:
        """Task that a layer class belongs to, None when the class names no task."""
        lqn = self.lqn
        cls = self.ensemble[e].classes[c]
        attr = getattr(cls, 'attribute', None)
        if attr is None or len(attr) < 2:
            return None
        tidx = None
        if attr[0] == LayeredNetworkElement.TASK:
            tidx = int(attr[1])
        elif attr[0] in (LayeredNetworkElement.ENTRY, LayeredNetworkElement.ACTIVITY):
            tidx = self._get_parent(int(attr[1]))
        elif attr[0] == LayeredNetworkElement.CALL:
            cidx = int(attr[1])
            if lqn.callpair is not None and cidx < len(lqn.callpair):
                tidx = self._get_parent(int(lqn.callpair[cidx, 0]))
        if tidx is None:
            return None
        tidx = int(tidx)
        if tidx < lqn.tshift or tidx >= lqn.tshift + lqn.ntasks:
            return None
        return tidx

    def _get_parent(self, idx: int) -> Optional[int]:
        """Get parent index for an element."""
        lqn = self.lqn
        if hasattr(lqn, 'parent') and lqn.parent is not None:
            if isinstance(lqn.parent, dict):
                return lqn.parent.get(idx)
            elif isinstance(lqn.parent, np.ndarray):
                # parent is 0-indexed over elements and carries -1 where an
                # element has no parent, since 0 is the first host
                if idx < len(lqn.parent):
                    val = lqn.parent[idx]
                    if isinstance(val, np.ndarray):
                        val = val.flatten()[0] if len(val) > 0 else -1
                    return int(val) if val >= 0 else None
        return None

    def _calls_of(self, aidx: int) -> List[int]:
        """Indices of the calls issued by activity aidx."""
        lqn = self.lqn
        callsof = getattr(lqn, 'callsof', None)
        if isinstance(callsof, dict):
            return [int(c) for c in callsof.get(aidx, [])]
        if callsof is not None and aidx < len(callsof):
            entry = callsof[aidx]
            if entry is None:
                return []
            return [int(c) for c in np.asarray(entry).flatten()]
        return []

    def _chain_ref_indices(self, layer_idx: int, classidx_0: int):
        """(refstat, refclass) of the chain holding class CLASSIDX_0 in layer LAYER_IDX."""
        refstat_k = None
        refclass_c = None
        if layer_idx < 0 or layer_idx >= len(self.ensemble) or self.ensemble[layer_idx] is None:
            return refstat_k, refclass_c
        layer_sn = self.ensemble[layer_idx]._sn if hasattr(self.ensemble[layer_idx], '_sn') else None
        if layer_sn is None:
            return refstat_k, refclass_c
        if getattr(layer_sn, 'chains', None) is not None:
            chains_arr = np.asarray(layer_sn.chains)
            if chains_arr.ndim == 2 and classidx_0 < chains_arr.shape[1]:
                for ch in range(chains_arr.shape[0]):
                    if chains_arr[ch, classidx_0] > 0:
                        if getattr(layer_sn, 'refclass', None) is not None:
                            rc = np.asarray(layer_sn.refclass).flatten()
                            if ch < len(rc):
                                refclass_c = int(rc[ch])
                        break
        if getattr(layer_sn, 'refstat', None) is not None:
            rs = np.asarray(layer_sn.refstat).flatten()
            if classidx_0 < len(rs):
                refstat_k = int(rs[classidx_0])
        return refstat_k, refclass_c

    def _is_activity_of_entry(self, aidx: int, eidx: int) -> bool:
        """Check if an activity is bound to an entry."""
        lqn = self.lqn
        if hasattr(lqn, 'graph') and lqn.graph is not None:
            if isinstance(lqn.graph, np.ndarray):
                if eidx <= lqn.graph.shape[0] and aidx <= lqn.graph.shape[1]:
                    return lqn.graph[eidx - 1, aidx - 1] > 0
        return False

    def _get_activity_successors(self, aidx: int, activity_classes: Dict[int, Any]) -> List[Tuple[int, float]]:
        """
        Get successor activities and their routing probabilities from lqn.graph.

        MATLAB equivalent: nextaidxs = find(lqn.graph(aidx,:)) in recurActGraph

        Args:
            aidx: Activity index
            activity_classes: Dict mapping activity index to class object

        Returns:
            List of (successor_aidx, probability) tuples for activities in same task
        """
        lqn = self.lqn
        successors = []

        if not hasattr(lqn, 'graph') or lqn.graph is None:
            return successors

        if not isinstance(lqn.graph, np.ndarray):
            return successors

        if aidx >= lqn.graph.shape[0]:
            return successors

        # Find successor activities in the graph (same task only, not entries)
        for succ_aidx in range(lqn.graph.shape[1]):
            prob = lqn.graph[aidx, succ_aidx]
            if prob > 0:
                # Check if successor is an activity in same task (not an entry/call target)
                if succ_aidx in activity_classes:
                    successors.append((succ_aidx, float(prob)))

        return successors

    def _compute_fork_scope(self, layer_model) -> Set[int]:
        """Compute the set of activity indices in the fork scope (between fork source and join output).

        For AND-fork layers, activities in the fork scope have visits that are 1/fanout
        of their correct values due to flat routing. This method identifies these activities
        so their visits can be corrected.

        Returns:
            Set of absolute activity indices in the fork scope.
        """
        lqn = self.lqn
        is_post_and = layer_model.attribute.get('is_post_and_act', set())
        is_pre_and = layer_model.attribute.get('is_pre_and_act', set())
        has_fork = layer_model.attribute.get('has_fork', False)

        if not has_fork or not is_post_and:
            return set()

        if not hasattr(lqn, 'graph') or lqn.graph is None or not isinstance(lqn.graph, np.ndarray):
            return set()

        # Get all activity indices in this layer
        acts_in_layer = set()
        classes = layer_model.get_classes()
        for cls in classes:
            if hasattr(cls, 'attribute') and cls.attribute is not None:
                if cls.attribute[0] == LayeredNetworkElement.ACTIVITY:
                    acts_in_layer.add(cls.attribute[1])

        # Find fork source: activity whose graph successors include POST_AND activities
        fork_source = None
        for aidx in acts_in_layer:
            if aidx < lqn.graph.shape[0]:
                successors = [j for j in range(lqn.graph.shape[1]) if lqn.graph[aidx, j] != 0]
                if any(s in is_post_and for s in successors):
                    fork_source = aidx
                    break

        if fork_source is None:
            return set()

        # Find join output: successor of PRE_AND activities that is not PRE_AND itself
        join_output = set()
        for aidx in is_pre_and:
            if aidx < lqn.graph.shape[0]:
                successors = [j for j in range(lqn.graph.shape[1]) if lqn.graph[aidx, j] != 0]
                for s in successors:
                    if s not in is_pre_and and s in acts_in_layer:
                        join_output.add(s)

        # Fork scope: all activities reachable from fork source, excluding join output
        scope = set()
        stack = [fork_source]
        while stack:
            a = stack.pop()
            if a in scope or a in join_output:
                continue
            scope.add(a)
            if a < lqn.graph.shape[0]:
                successors = [j for j in range(lqn.graph.shape[1]) if lqn.graph[a, j] != 0]
                for s in successors:
                    if s in acts_in_layer and s not in scope and s not in join_output:
                        stack.append(s)

        return scope

    def _apply_fork_visit_correction(self, layer_model):
        """Apply fork fanout correction to visit ratios in a layer model.

        For AND-fork layers with flat routing, the DTMC computes visits that are
        1/fanout for fork-scope activities. This multiplies their visits by fanout
        to restore correct values.
        """
        has_fork = layer_model.attribute.get('has_fork', False)
        maxfanout = layer_model.attribute.get('maxfanout', 1)

        if not has_fork or maxfanout <= 1:
            return

        # real Fork/Join nodes already get fanout via MMT; a post-hoc fanout multiplier here would double-count it.
        if layer_model.attribute.get('fork_node') is not None:
            return

        sn = layer_model._sn if hasattr(layer_model, '_sn') else None
        if sn is None or not hasattr(sn, 'visits') or sn.visits is None:
            return

        # fork-visit correction scales visits in place and is NOT idempotent; the corrected struct is marked so a surviving struct is never corrected twice.
        if getattr(sn, '_fork_visits_corrected', False):
            return
        sn._fork_visits_corrected = True

        # Compute or retrieve cached fork scope
        fork_scope = layer_model.attribute.get('_fork_scope')
        if fork_scope is None:
            fork_scope = self._compute_fork_scope(layer_model)
            layer_model.attribute['_fork_scope'] = fork_scope

        if not fork_scope:
            return

        # Find class indices for fork-scope activities
        fork_scope_class_indices = set()
        classes = layer_model.get_classes()
        for c_idx, cls in enumerate(classes):
            if hasattr(cls, 'attribute') and cls.attribute is not None:
                if cls.attribute[0] == LayeredNetworkElement.ACTIVITY:
                    aidx = cls.attribute[1]
                    if aidx in fork_scope:
                        fork_scope_class_indices.add(c_idx)

        if not fork_scope_class_indices:
            return

        # Multiply visits for fork-scope classes by fanout at all stateful nodes
        for c in range(len(sn.visits)):
            v = sn.visits[c]
            for k in fork_scope_class_indices:
                if k < v.shape[1]:
                    for ist in range(v.shape[0]):
                        v[ist, k] *= maxfanout

    def _find_caller_class_in_layer(self, caller_tidx: int, layer_idx: int) -> Optional[int]:
        """Find the class index for a caller task in a layer."""
        if layer_idx < 0 or layer_idx >= len(self.ensemble):
            return None

        layer = self.ensemble[layer_idx]
        if layer is None:
            return None

        tasks_matrix = layer.attribute.get('tasks', [])
        if isinstance(tasks_matrix, np.ndarray) and len(tasks_matrix) > 0:
            for row in range(tasks_matrix.shape[0]):
                if tasks_matrix[row, 1] == caller_tidx:
                    return int(tasks_matrix[row, 0])
        elif isinstance(tasks_matrix, list):
            for row in tasks_matrix:
                if len(row) > 1 and row[1] == caller_tidx:
                    return row[0]
        return None

    def _find_activity_class_in_layer(self, aidx: int, layer_idx: int) -> Optional[int]:
        """Find the class index for an activity in a layer."""
        if layer_idx < 0 or layer_idx >= len(self.ensemble):
            return None

        layer = self.ensemble[layer_idx]
        if layer is None:
            return None

        # Look for activity in the activities matrix
        activities_matrix = layer.attribute.get('activities', [])
        if isinstance(activities_matrix, np.ndarray) and len(activities_matrix) > 0:
            for row in range(activities_matrix.shape[0]):
                if activities_matrix[row, 1] == aidx:
                    return int(activities_matrix[row, 0])
        elif isinstance(activities_matrix, list):
            for row in activities_matrix:
                if len(row) > 1 and row[1] == aidx:
                    return row[0]
        return None

    def _find_call_class_in_layer(self, cidx: int, layer_idx: int) -> Optional[int]:
        """Find the class index for a call in a layer."""
        if layer_idx < 0 or layer_idx >= len(self.ensemble):
            return None

        layer = self.ensemble[layer_idx]
        if layer is None:
            return None

        # Look for call in the calls attribute
        # calls format: [class_index, cidx, src_aidx, tgt_eidx]
        calls_list = layer.attribute.get('calls', [])
        if isinstance(calls_list, np.ndarray) and len(calls_list) > 0:
            for row in range(calls_list.shape[0]):
                if calls_list[row, 1] == cidx:
                    return int(calls_list[row, 0])
        elif isinstance(calls_list, list):
            for row in calls_list:
                if len(row) > 1 and row[1] == cidx:
                    return row[0]
        return None

    def pre(self, it: int):
        """Operations before each iteration (matches MATLAB pre).

        Seed control for stochastic layer solvers.
        """
        if self.stochiter_mode is None or self.stochlayers is None:
            return
        if self.stochiter_mode == 'rm':
            # rm mode rotates per-layer seeds each iteration for independent noise, as Robbins-Monro averaging requires.
            for e in np.where(self.stochlayers)[0]:
                solver = self.solvers[e]
                if solver is not None and hasattr(solver, 'options'):
                    solver.options.seed = self.stochiter_seed_base + (it - 1) * self.nlayers + int(e) + 1
        elif self.stochiter_mode == 'crn':
            # crn mode pins a constant per-layer seed (sample-average approximation), carrying an O(1/sqrt(samples)) bias vs the true fixed point.
            for e in np.where(self.stochlayers)[0]:
                solver = self.solvers[e]
                if solver is not None and hasattr(solver, 'options'):
                    solver.options.seed = self.stochiter_seed_base + int(e) + 1

    def analyze(self, it: int, e: int) -> Tuple[Dict, float]:
        """
        Analyze a layer (matches MATLAB analyze).

        Returns:
            Tuple of (result dict, runtime)
        """
        import time
        t0 = time.time()
        solver_name = type(self.solvers[e]).__name__ if self.solvers[e] is not None else 'None'
        line_debug("LN analyze: iteration %d, layer %d (%s)", it, e, solver_name)

        result = {}

        try:
            solver = self.solvers[e]
            if solver is not None:
                # Get average metrics from solver (try different method names)
                if hasattr(solver, 'getAvg'):
                    QN, UN, RN, TN, AN, WN = solver.getAvg()
                elif hasattr(solver, 'get_avg'):
                    QN, UN, RN, TN, AN, WN = solver.get_avg()
                else:
                    raise AttributeError("Solver has no getAvg or get_avg method")

                # Sanitize results to prevent extreme values from MVA numerical instability
                max_val = 1e10
                for arr in [QN, UN, RN, TN, WN]:
                    if arr is not None and isinstance(arr, np.ndarray):
                        # Clamp extreme positive values
                        arr[arr > max_val] = np.nan
                        # Replace negative values (shouldn't happen) with nan
                        arr[arr < 0] = 0.0
                        # Replace inf with nan
                        arr[np.isinf(arr)] = np.nan

                result['QN'] = QN
                result['UN'] = UN
                result['RN'] = RN
                result['TN'] = TN
                result['AN'] = AN
                result['WN'] = WN

                # stochastic classification refreshed from the method actually resolved at runtime, captured before post() resets the layer solvers.
                if it == 1 and self.stochlayers is not None and hasattr(solver, 'isStochastic'):
                    self.stochlayers[e] = solver.isStochastic()

                # warm-start the next AMVA solve from the current chain-aggregated queue lengths; see _kb/06-solver-catalog.md LN AMVA warm-start section.
                if solver_name == 'SolverMVA' and isinstance(QN, np.ndarray):
                    sne = self.ensemble[e].getStruct()
                    if QN.shape == (sne.nstations, sne.nclasses):
                        chains = np.asarray(sne.chains)
                        Qch = np.zeros((sne.nstations, sne.nchains))
                        QNfin = np.nan_to_num(QN, nan=0.0, posinf=0.0, neginf=0.0)
                        for c in range(sne.nchains):
                            cls = np.where(chains[c, :] > 0)[0]
                            Qch[:, c] = np.sum(QNfin[:, cls], axis=1)
                        solver.options.init_sol = Qch
        except Exception as ex:
            # If solver fails, use previous iteration if available
            if it > 1 and len(self.results) >= it - 1:
                prev_result = self.results[it - 2][e]
                result = prev_result.copy()
            else:
                raise

        runtime = time.time() - t0
        return result, runtime

    def post(self, it: int):
        """Operations after each iteration (matches MATLAB post)."""
        line_debug("LN post: iteration %d, updating metrics and layer parameters", it)
        # Update metrics
        self.update_metrics(it)

        # Update think times
        self.update_think_times(it)

        # Update populations if interlocking enabled
        if self.options.config.get('interlocking', False):
            self.update_populations(it)

        # Update layer parameters
        self.update_layers(it)

        # Update routing probabilities
        self.update_routing_probabilities(it)

        # refresh_rates() when only service times changed, full refresh_struct() when routing invalidated; mirrors MATLAB refreshRates/refreshChains split.
        for e in range(self.nlayers):
            if e < len(self.ensemble) and self.ensemble[e] is not None:
                if not self.ensemble[e]._has_struct:
                    # Struct was invalidated (routing changed) - full rebuild needed
                    self.ensemble[e].refresh_struct()
                    # The rebuild can change the chain basis, invalidating the
                    # warm-start solution cached by analyze()
                    if e < len(self.solvers) and self.solvers[e] is not None \
                            and hasattr(self.solvers[e], 'options'):
                        self.solvers[e].options.init_sol = None
                elif self._is_ph_encoding():
                    # a phase-type service law, whose phases a rate-only refresh
                    # would drop -- see _kb/06-solver-catalog.md (LN section)
                    self.ensemble[e].refresh_struct()
                else:
                    # Only service rates changed - lightweight update
                    # (may trigger full rebuild if _sn was set to None by set_service)
                    self.ensemble[e].refresh_rates()
                # Re-apply fork visit correction (no-op for non-fork layers)
                self._apply_fork_visit_correction(self.ensemble[e])
                # Reset solver to force recomputation
                if e < len(self.solvers) and self.solvers[e] is not None:
                    if hasattr(self.solvers[e], 'reset'):
                        self.solvers[e].reset()

        # Refresh layer structure if interlocking enabled (to update populations)
        if self.options.config.get('interlocking', False):
            for e in range(self.nlayers):
                if self.ensemble[e] is not None:
                    if not self.ensemble[e]._has_struct:
                        self.ensemble[e].refresh_struct()
                        self._apply_fork_visit_correction(self.ensemble[e])

        # Disable checks after first iteration
        if it == 1:
            for e in range(self.nlayers):
                if self.solvers[e] is not None:
                    if hasattr(self.solvers[e], 'set_checks'):
                        self.solvers[e].set_checks(False)

    def update_metrics(self, it: int):
        """Update metrics (matches MATLAB updateMetrics)."""
        method = self.lnmethod
        if self._is_ph_encoding():
            # see _kb/06-solver-catalog.md (LN section) for rationale
            self._update_metrics_ph(it)
        elif method == 'moment3':
            self._update_metrics_moment_based(it)
        else:
            self._update_metrics_default(it)

    def _layer_respt_cdf(self, repo, layer_idx):
        """Per-(station, class) response time CDFs of one layer, memoised in REPO.

        The fluid passage time is asked for first, as the reference does, and the
        LAYER's own solver answers when it refuses -- python's fluid getter raises
        ``passage-time integration failed (stiff augmented system)`` on the PS and
        INF layers an LQN is mostly made of.

        Returns:
            RD[station][class], an ``(n, 2)`` ``[cdf, time]`` array or None per
            cell, or None when neither solver produced anything.
        """
        if layer_idx in repo:
            return repo[layer_idx]
        raw = None
        try:
            from ..solver_fld import SolverFLD
            raw = SolverFLD(self.ensemble[layer_idx]).getCdfRespT()
        except Exception:
            try:
                raw = self.solvers[layer_idx].getCdfRespT()
            except Exception:
                raw = None
        repo[layer_idx] = self._nested_respt_cdf(raw, self.ensemble[layer_idx])
        return repo[layer_idx]

    @staticmethod
    def _nested_respt_cdf(raw, layer):
        """Bring either getCdfRespT contract to RD[station][class] = [cdf, time].

        TWO CONTRACTS MEET HERE, and the entry assembly indexes only the first.
        ``SolverFLD`` returns the NESTED shape, ``RD[station][class]`` an
        ``(n, 2)`` ``[cdf, time]`` array. ``SolverMVA`` and ``SolverNC`` return
        the FLAT native contract, a list of dicts carrying 1-based ``station`` and
        ``class`` with numpy ``t`` and ``p`` -- the same contract
        ``cpp_dispatch.cdf_respt_via_cpp`` documents.

        The flat one is CONVERTED here rather than dropped. The isinstance guard
        this replaces tested the station row for ``list`` and silently discarded
        every dict, so on any model whose fluid passage time failed -- which is
        the common case -- a caller entry convolved its own host demand and NONE
        of its call terms, and reported an entry service time equal to that bare
        demand while the activity row beside it carried the full value.
        """
        if raw is None:
            return None
        if not hasattr(raw, '__len__') or len(raw) == 0:
            return None
        first = next((cell for cell in raw if cell is not None), None)
        if first is None:
            return None
        if isinstance(first, (list, tuple)):
            return raw
        if isinstance(first, dict):
            sn = layer.getStruct()
            M, K = int(sn.nstations), int(sn.nclasses)
            RD = [[None] * K for _ in range(M)]
            for cell in raw:
                if cell is None:
                    continue
                # 1-based on the wire, as the native contract specifies
                i = int(cell['station']) - 1
                r = int(cell['class']) - 1
                if not (0 <= i < M and 0 <= r < K):
                    raise ValueError(
                        "getCdfRespT returned station=%d class=%d, outside the layer's "
                        "%d x %d index space" % (i + 1, r + 1, M, K))
                t = np.asarray(cell['t'], dtype=float).ravel()
                p = np.asarray(cell['p'], dtype=float).ravel()
                if t.size != p.size:
                    raise ValueError(
                        "getCdfRespT cell (station=%d, class=%d) carries %d times and "
                        "%d probabilities" % (i + 1, r + 1, t.size, p.size))
                RD[i][r] = np.column_stack([p, t])
            return RD
        raise TypeError(
            "unrecognised getCdfRespT return shape: expected the nested "
            "RD[station][class] arrays of SolverFLD or the flat list of "
            "{station, class, t, p} dicts of SolverMVA/SolverNC, got a %s"
            % type(first).__name__)

    def _update_metrics_moment_based(self, it: int):
        """Moment-based metrics update (matches MATLAB updateMetricsMomentBased)."""
        lqn = self.lqn

        if not self.hasconverged:
            # ===== PRE-CONVERGENCE: Mean-based propagation using exponential fits =====

            # First obtain servt of activities at hostlayers
            self.servt = np.zeros(lqn.nidx)
            self.residt = np.zeros(lqn.nidx)

            if self.servt_classes_updmap is not None:
                for r in range(len(self.servt_classes_updmap)):
                    idx = int(self.servt_classes_updmap[r, 0])
                    aidx = int(self.servt_classes_updmap[r, 1])
                    nodeidx = int(self.servt_classes_updmap[r, 2])
                    classidx = int(self.servt_classes_updmap[r, 3])
                    layer_idx = int(self.idxhash[idx])
                    nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                    classidx_0 = classidx - 1 if classidx >= 1 else 0

                    if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                        result = self.results[-1][layer_idx]
                        if result is not None and 'RN' in result:
                            RN = result['RN']
                            TN = result['TN']
                            QN = result.get('QN')
                            WN = result.get('WN', RN)
                            if RN is not None and nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                                self.servt[aidx] = RN[nodeidx_0, classidx_0]
                                self.tput[aidx] = TN[nodeidx_0, classidx_0]
                                if self.servt[aidx] > 0:
                                    self.servtproc[aidx] = Exp.fit_mean(self.servt[aidx])

                                # Compute residt from QN/TN_ref (matching updateMetricsDefault)
                                refstat_k = None
                                refclass_c = None
                                if layer_idx < len(self.ensemble) and self.ensemble[layer_idx] is not None:
                                    layer_sn = self.ensemble[layer_idx]._sn if hasattr(self.ensemble[layer_idx], '_sn') else None
                                    if layer_sn is not None and hasattr(layer_sn, 'chains') and layer_sn.chains is not None:
                                        chains_arr = np.asarray(layer_sn.chains)
                                        if chains_arr.ndim == 2 and classidx_0 < chains_arr.shape[1]:
                                            for ch in range(chains_arr.shape[0]):
                                                if chains_arr[ch, classidx_0] > 0:
                                                    if hasattr(layer_sn, 'refclass') and layer_sn.refclass is not None:
                                                        rc = np.asarray(layer_sn.refclass).flatten()
                                                        if ch < len(rc):
                                                            refclass_c = int(rc[ch])
                                                    break
                                        if hasattr(layer_sn, 'refstat') and layer_sn.refstat is not None:
                                            rs = np.asarray(layer_sn.refstat).flatten()
                                            if classidx_0 < len(rs):
                                                refstat_k = int(rs[classidx_0])

                                if (refstat_k is not None and refclass_c is not None and
                                        QN is not None and TN is not None and
                                        0 <= refstat_k < TN.shape[0] and 0 <= refclass_c < TN.shape[1]):
                                    TN_ref = TN[refstat_k, refclass_c]
                                    if TN_ref > 1e-8:  # GlobalConstants.FineTol
                                        self.residt[aidx] = QN[nodeidx_0, classidx_0] / TN_ref
                                    else:
                                        self.residt[aidx] = WN[nodeidx_0, classidx_0] if WN is not None else RN[nodeidx_0, classidx_0]
                                else:
                                    self.residt[aidx] = WN[nodeidx_0, classidx_0] if WN is not None else RN[nodeidx_0, classidx_0]

                                # An activity think time is in series with the host demand
                                zt_act = self._act_thinktime(aidx)
                                if zt_act > 0:
                                    self.servt[aidx] += zt_act
                                    self.residt[aidx] += zt_act
                                    self.servtproc[aidx] = Exp.fit_mean(self.servt[aidx])

                                # async-only targets carry no visit-ratio scaling (matching _update_metrics_default)
                                if lqn.ashift <= aidx < lqn.ashift + lqn.nacts:
                                    if hasattr(lqn, 'graph') and isinstance(lqn.graph, np.ndarray):
                                        for eidx in range(lqn.eshift, lqn.eshift + lqn.nentries):
                                            if eidx < lqn.graph.shape[0] and aidx < lqn.graph.shape[1] and lqn.graph[eidx, aidx] > 0:
                                                has_sync_callers = False
                                                has_async_callers = False
                                                if isinstance(getattr(lqn, 'issynccaller', None), np.ndarray) and eidx < lqn.issynccaller.shape[1]:
                                                    has_sync_callers = np.any(lqn.issynccaller[:, eidx])
                                                if isinstance(getattr(lqn, 'isasynccaller', None), np.ndarray) and eidx < lqn.isasynccaller.shape[1]:
                                                    has_async_callers = np.any(lqn.isasynccaller[:, eidx])
                                                if has_async_callers and not has_sync_callers:
                                                    self.residt[aidx] = self.servt[aidx]
                                                break

            # Estimate call response times at hostlayers
            self.callservt = np.zeros(lqn.ncalls)
            self.callresidt = np.zeros(lqn.ncalls)

            if self.call_classes_updmap is not None:
                for c in range(len(self.call_classes_updmap)):
                    idx = int(self.call_classes_updmap[c, 0])
                    cidx = int(self.call_classes_updmap[c, 1])
                    nodeidx = int(self.call_classes_updmap[c, 2])
                    classidx = int(self.call_classes_updmap[c, 3])

                    if nodeidx > 1:
                        layer_idx = int(self.idxhash[idx])
                        if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                            result = self.results[-1][layer_idx]
                            if result is not None and 'RN' in result:
                                RN = result['RN']
                                WN = result.get('WN', RN)
                                nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                                classidx_0 = classidx - 1 if classidx >= 1 else 0
                                if nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                                    if nodeidx == 1:
                                        self.callservt[cidx] = 0.0
                                        self.callresidt[cidx] = 0.0
                                    else:
                                        # Include call multiplicity (matching updateMetricsDefault)
                                        call_mean = self._get_call_mean(cidx)
                                        fcr_wait = self._region_wait(layer_idx, nodeidx_0, classidx_0, result)
                                        self.callservt[cidx] = (RN[nodeidx_0, classidx_0] + fcr_wait) * call_mean
                                        # callresidt uses WN which already includes visit multiplicity
                                        self.callresidt[cidx] = WN[nodeidx_0, classidx_0] + fcr_wait

            # Resolve the entry servt summing up these contributions; the terms are
            # residence times (Vtask=1), rescaled to Ventry=1 by the task/entry tput ratio below
            # entry_servt = (I - servtmatrix)^(-1) * [residt; callresidt]
            size = lqn.nidx + lqn.ncalls
            combined_vec = np.zeros(size)
            combined_vec[:lqn.nidx] = self.residt
            combined_vec[lqn.nidx:lqn.nidx + lqn.ncalls] = self.callresidt

            identity = np.eye(size)
            system = identity - self.servtmatrix
            try:
                entry_servt = np.linalg.solve(system, combined_vec)
            except np.linalg.LinAlgError:
                entry_servt = np.linalg.lstsq(system, combined_vec, rcond=None)[0]

            # Clear entries up to eshift
            entry_servt[:lqn.eshift] = 0

            # NO forwarding propagation here. _lqn_fwd_rendezvous has already
            # reconnected every forwarding chain reachable from a synchronous call
            # to the client that issued the rendezvous (Franks 1999, Sec. 3.3.1),
            # so the forwarded service is in the caller's chain before this runs;
            # adding it again inflated the caller by exactly the forwarded entry's
            # mean. An asynchronous call into a chain is left untouched there by
            # design -- a send-no-reply does not block -- so it must not accumulate
            # the forwarded service either. See BUGS.md BUG-91.

            # A SetupTask's cold start is charged HERE, to the entry, and with the
            # probability that the thread was actually found powered down. It is not
            # host demand, so it does not belong to any activity's residence:
            # reporting it there put RespT(A2) at 1.29479 on lqn_setup against the
            # 0.333178 LDES measures, which is the bare demand. See _setup_charge.
            for i in range(lqn.eshift, lqn.eshift + lqn.nentries):
                entry_servt[i] += self._setup_charge(self._get_parent(i))

            # Update servt for entries
            for i in range(lqn.eshift, lqn.eshift + lqn.nentries):
                self.servt[i] = entry_servt[i]

            # Clear activities after ashift
            for i in range(lqn.ashift, len(entry_servt)):
                entry_servt[i] = 0

            # Published for inspection, as MATLAB's SolverLN carries entry_servt
            # on the object: a debug driver reads it back after an iteration and
            # a local would leave it unreachable.
            self.entry_servt = entry_servt.copy()

            # Compute entry-level residt using servtmatrix and activity residt
            combined_residt = np.zeros(size)
            combined_residt[:lqn.nidx] = self.residt
            combined_residt[lqn.nidx:lqn.nidx + lqn.ncalls] = self.callresidt
            entry_residt_vec = self.servtmatrix @ combined_residt
            entry_residt_vec[:lqn.eshift] = 0

            # Scale entry residt/servt by task/entry throughput ratio
            for e in range(lqn.nentries):
                eidx = lqn.eshift + e
                tidx = self._get_parent(eidx)
                hidx = self._get_parent(tidx) if tidx is not None else None
                if tidx is None or hidx is None:
                    continue
                if self.ignore[tidx] or self.ignore[hidx]:
                    continue

                has_sync_callers = self._has_sync_callers_for_entry(eidx)

                if has_sync_callers:
                    tput_ratio = self._get_entry_tput_ratio(eidx, tidx, hidx)
                    if tput_ratio is not None:
                        task_tput, entry_tput = tput_ratio
                        if entry_tput > GlobalConstants.Zero:
                            self.servt[eidx] = entry_servt[eidx] * task_tput / entry_tput
                            self.residt[eidx] = entry_residt_vec[eidx] * task_tput / entry_tput
                        else:
                            self.residt[eidx] = entry_residt_vec[eidx]
                    else:
                        self.residt[eidx] = entry_residt_vec[eidx]
                else:
                    self.residt[eidx] = entry_residt_vec[eidx]

            # Update servtproc for entries
            if self.call_classes_updmap is not None:
                for row in self.call_classes_updmap:
                    cidx = int(row[1])
                    nodeidx = int(row[2])
                    if nodeidx > 1:
                        eidx = self._get_call_target_entry(cidx)
                        if eidx is not None and eidx > 0 and eidx < len(self.servt):
                            if self.servt[eidx] > 0:
                                self.servtproc[eidx] = Exp.fit_mean(self.servt[eidx])

            # Determine call response times processes
            if self.call_classes_updmap is not None:
                for row in self.call_classes_updmap:
                    cidx = int(row[1])
                    nodeidx = int(row[2])
                    if nodeidx > 1:
                        eidx = self._get_call_target_entry(cidx)
                        if eidx is not None and eidx > 0:
                            if it == 1:
                                if eidx < len(self.servt):
                                    self.callservt[cidx] = self.servt[eidx]
                                if eidx < len(self.servtproc) and self.servtproc[eidx] is not None:
                                    self.callservtproc[cidx] = self.servtproc[eidx]
                            else:
                                if self.callservt[cidx] > 0:
                                    self.callservtproc[cidx] = Exp.fit_mean(self.callservt[cidx])

        else:
            # ===== POST-CONVERGENCE: Full CDF-based 3-moment APH fitting =====
            from ...api.kpctoolbox.aph import aph_convseq, aph_simplify
            from ...api.butools.ph.canonical import APHFrom3Moments
            from ...distributions.markovian import APH

            self.servtcdf = [None] * lqn.nidx
            repo = {}

            # First obtain servt of activities at hostlayers
            self.servt = np.zeros(lqn.nidx)
            self.residt = np.zeros(lqn.nidx)

            if self.servt_classes_updmap is not None:
                for r in range(len(self.servt_classes_updmap)):
                    idx = int(self.servt_classes_updmap[r, 0])
                    aidx = int(self.servt_classes_updmap[r, 1])
                    nodeidx = int(self.servt_classes_updmap[r, 2])
                    classidx = int(self.servt_classes_updmap[r, 3])
                    layer_idx = int(self.idxhash[idx])
                    nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                    classidx_0 = classidx - 1 if classidx >= 1 else 0

                    if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                        result = self.results[-1][layer_idx]
                        if result is not None and 'TN' in result:
                            self.tput[aidx] = result['TN'][nodeidx_0, classidx_0]

                        # Compute residt from QN/TN_ref
                        if result is not None and 'RN' in result:
                            QN = result.get('QN')
                            TN = result['TN']
                            WN = result.get('WN', result['RN'])
                            refstat_k = None
                            refclass_c = None
                            if layer_idx < len(self.ensemble) and self.ensemble[layer_idx] is not None:
                                layer_sn = self.ensemble[layer_idx]._sn if hasattr(self.ensemble[layer_idx], '_sn') else None
                                if layer_sn is not None and hasattr(layer_sn, 'chains') and layer_sn.chains is not None:
                                    chains_arr = np.asarray(layer_sn.chains)
                                    if chains_arr.ndim == 2 and classidx_0 < chains_arr.shape[1]:
                                        for ch in range(chains_arr.shape[0]):
                                            if chains_arr[ch, classidx_0] > 0:
                                                if hasattr(layer_sn, 'refclass') and layer_sn.refclass is not None:
                                                    rc = np.asarray(layer_sn.refclass).flatten()
                                                    if ch < len(rc):
                                                        refclass_c = int(rc[ch])
                                                break
                                    if hasattr(layer_sn, 'refstat') and layer_sn.refstat is not None:
                                        rs = np.asarray(layer_sn.refstat).flatten()
                                        if classidx_0 < len(rs):
                                            refstat_k = int(rs[classidx_0])

                            if (refstat_k is not None and refclass_c is not None and
                                    QN is not None and TN is not None and
                                    0 <= refstat_k < TN.shape[0] and 0 <= refclass_c < TN.shape[1]):
                                TN_ref = TN[refstat_k, refclass_c]
                                if TN_ref > 1e-8:
                                    self.residt[aidx] = QN[nodeidx_0, classidx_0] / TN_ref
                                else:
                                    self.residt[aidx] = WN[nodeidx_0, classidx_0] if WN is not None else 0
                            else:
                                self.residt[aidx] = WN[nodeidx_0, classidx_0] if WN is not None else 0

                    # Get CDFs - try SolverFluid first, fall back to layer solver
                    cdf_data = self._layer_respt_cdf(repo, layer_idx)
                    if (cdf_data is not None and nodeidx_0 < len(cdf_data)
                            and cdf_data[nodeidx_0] is not None
                            and classidx_0 < len(cdf_data[nodeidx_0])):
                        self.servtcdf[aidx] = cdf_data[nodeidx_0][classidx_0]

            # Initialize callservtcdf
            self.callservtcdf = [None] * lqn.ncalls
            self.callservt = np.zeros(lqn.ncalls)
            self.callresidt = np.zeros(lqn.ncalls)

            if self.call_classes_updmap is not None:
                for c in range(len(self.call_classes_updmap)):
                    idx = int(self.call_classes_updmap[c, 0])
                    cidx = int(self.call_classes_updmap[c, 1])
                    nodeidx = int(self.call_classes_updmap[c, 2])
                    classidx = int(self.call_classes_updmap[c, 3])

                    if nodeidx > 1:
                        layer_idx = int(self.idxhash[idx])
                        nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                        classidx_0 = classidx - 1 if classidx >= 1 else 0

                        cdf_data = self._layer_respt_cdf(repo, layer_idx)
                        if (cdf_data is not None and nodeidx_0 < len(cdf_data)
                                and cdf_data[nodeidx_0] is not None
                                and classidx_0 < len(cdf_data[nodeidx_0])):
                            self.callservtcdf[cidx] = cdf_data[nodeidx_0][classidx_0]

                        # Also set callresidt from WN
                        if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                            result = self.results[-1][layer_idx]
                            if result is not None and 'WN' in result:
                                WN = result['WN']
                                if WN is not None and nodeidx_0 < WN.shape[0] and classidx_0 < WN.shape[1]:
                                    self.callresidt[cidx] = WN[nodeidx_0, classidx_0] \
                                        + self._region_wait(layer_idx, nodeidx_0, classidx_0, result)

            # Build combined CDF list (servtcdf + callservtcdf)
            cdf = self.servtcdf + self.callservtcdf

            # Resolve entry service times using matrix inversion
            size = lqn.nidx + lqn.ncalls
            identity = np.eye(size)
            system = identity - self.servtmatrix
            try:
                matrix = np.linalg.inv(system)
            except np.linalg.LinAlgError:
                matrix = np.linalg.pinv(system)

            # Process each entry
            for i in range(lqn.nentries):
                eidx = lqn.eshift + i

                # Find contributing indices (where matrix[eidx,:] > 0)
                convolidx = []
                for j in range(matrix.shape[1]):
                    if matrix[eidx, j] > 0 and (j >= lqn.eshift + lqn.nentries):
                        convolidx.append(j)

                # Build APH convolution list
                param_list = []

                for fitidx in convolidx:
                    cdf_data = cdf[fitidx] if fitidx < len(cdf) else None
                    if cdf_data is None:
                        continue

                    # Extract raw moments from CDF data
                    # CDF data is 2D array with columns [cdf_vals, times]
                    if isinstance(cdf_data, np.ndarray) and cdf_data.ndim == 2 and cdf_data.shape[1] >= 2:
                        cdf_vals = cdf_data[:, 0]
                        times = cdf_data[:, 1]
                        # bin midpoints weighted by the CDF increment, as EmpiricalCDF.getMoments
                        x = times[:-1] + np.diff(times) / 2.0
                        dF = np.diff(cdf_vals)
                        m1 = np.sum(x * dF)
                        m2 = np.sum(x**2 * dF)
                        m3 = np.sum(x**3 * dF)
                    else:
                        continue

                    # An activity think time is in series with the host demand,
                    # so its raw moments convolve with the measured ones before
                    # the APH fit, as MATLAB's lqn_act_thinktime block does
                    if fitidx < lqn.nidx and self._act_thinktime(fitidx) > 0:
                        ztd = self.actthinkproc[fitidx]
                        t1 = ztd.getMean()
                        sig2 = ztd.getSCV() * t1 ** 2
                        t2 = sig2 + t1 ** 2
                        t3 = ztd.getSkewness() * sig2 ** 1.5 + 3 * t1 * t2 - 2 * t1 ** 3
                        m3 = m3 + 3 * m2 * t1 + 3 * m1 * t2 + t3
                        m2 = m2 + 2 * m1 * t1 + t2
                        m1 = m1 + t1

                    # Use CoarseTol to skip near-zero mean CDFs
                    if m1 > GlobalConstants.CoarseTol:
                        try:
                            alpha, T = APHFrom3Moments([m1, m2, m3])
                        except Exception:
                            continue

                        # For call indices, multiply repetitions by mean number of calls
                        reps = matrix[eidx, fitidx]
                        # The CALL block starts AT nidx in this 0-based index
                        # space; MATLAB numbers from 1, so its `> nidx` becomes
                        # `>= nidx` here. Dead code until the layer CDF repo
                        # started supplying call terms, and wrong the moment it
                        # was not: it credited call 0 to the last activity.
                        if fitidx >= lqn.nidx:
                            cidx_local = fitidx - lqn.nidx
                            reps = reps * self._get_call_mean(cidx_local)

                        integer_reps = int(np.floor(reps))
                        fractional_part = reps - integer_reps

                        if fractional_part == 0:
                            for _ in range(integer_reps):
                                param_list.append((alpha, T))
                        elif integer_reps > 0 and fractional_part > 0:
                            for _ in range(integer_reps):
                                param_list.append((alpha, T))
                            try:
                                zero_alpha, zero_T = APHFrom3Moments([1e-8, 2e-16, 6e-24])
                                alpha_br, T_br = aph_simplify(
                                    alpha, T, zero_alpha, zero_T,
                                    fractional_part, 1.0 - fractional_part, 3)
                                param_list.append((alpha_br, T_br))
                            except Exception:
                                pass
                        else:
                            try:
                                zero_alpha, zero_T = APHFrom3Moments([1e-8, 2e-16, 6e-24])
                                alpha_br, T_br = aph_simplify(
                                    alpha, T, zero_alpha, zero_T,
                                    fractional_part, 1.0 - fractional_part, 3)
                                param_list.append((alpha_br, T_br))
                            except Exception:
                                pass

                        # Update servtproc and callservtproc
                        # Same 0-based boundary: servtproc holds nidx entries,
                        # so `fitidx == nidx` is the FIRST CALL, not the last
                        # activity, and writing it here raised IndexError.
                        if fitidx < lqn.nidx:
                            self.servtproc[fitidx] = Exp.fit_mean(m1)
                            self.servt[fitidx] = m1
                        else:
                            self.callservtproc[fitidx - lqn.nidx] = Exp.fit_mean(m1)
                            self.callservt[fitidx - lqn.nidx] = m1

                # Convolve all contributions
                if not param_list:
                    self.servt[eidx] = 0
                else:
                    entry_dist = None
                    try:
                        alpha_conv, T_conv = aph_convseq(param_list)
                        entry_dist = APH(alpha_conv, T_conv)
                    except Exception:
                        self.servt[eidx] = 0
                    if entry_dist is not None:
                        # ENTRY-LOCAL INDEX, 0-BASED. MATLAB numbers the entries
                        # 1..nentries and guards `0 < e <= nentries`; this index
                        # space runs 0..nentries-1, so carrying that guard over
                        # verbatim silently dropped entry 0 -- its law was never
                        # stored and getCdfRespT had nothing to return for it.
                        entry_index = eidx - lqn.eshift
                        if self.entryproc is None:
                            self.entryproc = [None] * lqn.nentries
                        self.entryproc[entry_index] = entry_dist
                        self.servt[eidx] = entry_dist.getMean()
                        self.servtproc[eidx] = Exp.fit_mean(self.servt[eidx])
                        # Unguarded, as the reference is: a law that cannot be
                        # tabulated is a defect to surface, not a service time
                        # to quietly replace with zero.
                        self.entrycdfrespt[entry_index] = entry_dist.evalCDF()

            # fallback linear-system entry servt when APH fitting fails (CDF unavailable), matching MATLAB's SolverFluid-always-succeeds assumption.
            any_zero_entry = any(
                self.servt[lqn.eshift + i] == 0
                for i in range(lqn.nentries)
            )
            if any_zero_entry:
                # Rebuild activity-level servt from results for the system solve
                fallback_servt = np.zeros(lqn.nidx)
                if self.servt_classes_updmap is not None:
                    for r in range(len(self.servt_classes_updmap)):
                        idx = int(self.servt_classes_updmap[r, 0])
                        aidx = int(self.servt_classes_updmap[r, 1])
                        nodeidx = int(self.servt_classes_updmap[r, 2])
                        classidx = int(self.servt_classes_updmap[r, 3])
                        layer_idx = int(self.idxhash[idx])
                        nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                        classidx_0 = classidx - 1 if classidx >= 1 else 0
                        if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                            result = self.results[-1][layer_idx]
                            if result is not None and 'RN' in result:
                                RN = result['RN']
                                if nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                                    fallback_servt[aidx] = RN[nodeidx_0, classidx_0]

                fallback_callservt = np.zeros(lqn.ncalls)
                if self.call_classes_updmap is not None:
                    for c in range(len(self.call_classes_updmap)):
                        idx = int(self.call_classes_updmap[c, 0])
                        cidx = int(self.call_classes_updmap[c, 1])
                        nodeidx = int(self.call_classes_updmap[c, 2])
                        classidx = int(self.call_classes_updmap[c, 3])
                        if nodeidx > 1:
                            layer_idx = int(self.idxhash[idx])
                            nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                            classidx_0 = classidx - 1 if classidx >= 1 else 0
                            if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                                result = self.results[-1][layer_idx]
                                if result is not None and 'RN' in result:
                                    RN = result['RN']
                                    if nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                                        call_mean = self._get_call_mean(cidx)
                                        fallback_callservt[cidx] = (RN[nodeidx_0, classidx_0]
                                            + self._region_wait(layer_idx, nodeidx_0, classidx_0, result)) * call_mean

                # Solve (I - servtmatrix) * entry_servt = [servt; callservt]
                combined_vec = np.zeros(size)
                combined_vec[:lqn.nidx] = fallback_servt
                combined_vec[lqn.nidx:lqn.nidx + lqn.ncalls] = fallback_callservt
                try:
                    entry_servt_fb = np.linalg.solve(system, combined_vec)
                except np.linalg.LinAlgError:
                    entry_servt_fb = np.linalg.lstsq(system, combined_vec, rcond=None)[0]
                entry_servt_fb[:lqn.eshift] = 0

                for i in range(lqn.nentries):
                    eidx = lqn.eshift + i
                    if self.servt[eidx] == 0 and entry_servt_fb[eidx] > 0:
                        self.servt[eidx] = entry_servt_fb[eidx]
                        self.servtproc[eidx] = Exp.fit_mean(self.servt[eidx])

            # NO forwarding propagation here, for the reason given at the
            # entry_servt assembly above: _lqn_fwd_rendezvous has already charged
            # the forwarded service to the caller. See BUGS.md BUG-91.

            # Compute entry-level residt
            combined_residt = np.zeros(size)
            combined_residt[:lqn.nidx] = self.residt
            combined_residt[lqn.nidx:lqn.nidx + lqn.ncalls] = self.callresidt
            entry_residt_vec = self.servtmatrix @ combined_residt
            entry_residt_vec[:lqn.eshift] = 0

            for e in range(lqn.nentries):
                eidx = lqn.eshift + e
                tidx = self._get_parent(eidx)
                hidx = self._get_parent(tidx) if tidx is not None else None
                if tidx is None or hidx is None:
                    continue
                if self.ignore[tidx] or self.ignore[hidx]:
                    continue

                has_sync_callers = self._has_sync_callers_for_entry(eidx)

                if has_sync_callers:
                    tput_ratio = self._get_entry_tput_ratio(eidx, tidx, hidx)
                    if tput_ratio is not None:
                        task_tput, entry_tput = tput_ratio
                        if entry_tput > GlobalConstants.Zero:
                            self.residt[eidx] = entry_residt_vec[eidx] * task_tput / entry_tput
                        else:
                            self.residt[eidx] = entry_residt_vec[eidx]
                    else:
                        self.residt[eidx] = entry_residt_vec[eidx]
                else:
                    self.residt[eidx] = entry_residt_vec[eidx]

            # Determine call response times processes (final loop)
            if self.call_classes_updmap is not None:
                for row in self.call_classes_updmap:
                    cidx = int(row[1])
                    nodeidx = int(row[2])
                    if nodeidx > 1:
                        eidx = self._get_call_target_entry(cidx)
                        if eidx is not None and eidx > 0:
                            if it == 1:
                                if eidx < len(self.servt):
                                    self.callservt[cidx] = self.servt[eidx]
                                if eidx < len(self.servtproc) and self.servtproc[eidx] is not None:
                                    self.callservtproc[cidx] = Exp.fit_mean(self.servt[eidx])

            # This pass IS the moment3 answer, and it is TERMINAL. Its entry laws
            # are convolutions of the activities' own response distributions; the
            # pre-convergence branch instead reads QN/TN_ref, a residence per
            # REFERENCE cycle, which the entry assembly then treats as a
            # per-entry-visit time. The two disagree by the entry's visit ratio
            # whenever it is not 1, so letting the iteration fall back to that
            # branch after this one has run DISCARDS the moment-based laws and
            # reports the other quantity. See BUGS.md BUG-97.
            self.moment_pass_done = True

    def _has_sync_callers_for_entry(self, eidx: int) -> bool:
        """Check if entry has sync callers."""
        lqn = self.lqn
        if hasattr(lqn, 'callpair') and lqn.callpair is not None and hasattr(lqn, 'calltype'):
            for cidx in range(lqn.ncalls):
                if cidx < len(lqn.callpair):
                    tgt_eidx = int(lqn.callpair[cidx, 1]) if lqn.callpair[cidx, 1] > 0 else 0
                    if tgt_eidx == eidx:
                        calltype = lqn.calltype[cidx] if cidx < len(lqn.calltype) else 0
                        is_sync = (calltype == CallType.SYNC or
                                  calltype == CallType.SYNC.value or
                                  (isinstance(calltype, (int, np.integer)) and int(calltype) == CallType.SYNC.value))
                        if is_sync:
                            return True
        return False

    def _get_entry_tput_ratio(self, eidx: int, tidx: int, hidx: int) -> Optional[Tuple[float, float]]:
        """Get task/entry throughput ratio from host layer results."""
        if np.isnan(self.idxhash[hidx]):
            return None
        layer_idx = int(self.idxhash[hidx])
        if layer_idx < 0 or layer_idx >= len(self.ensemble) or self.ensemble[layer_idx] is None:
            return None
        layer = self.ensemble[layer_idx]
        result = self.results[-1][layer_idx] if len(self.results) > 0 and layer_idx < len(self.results[-1]) else None
        if result is None or 'TN' not in result:
            return None

        TN = result['TN']
        client_idx = layer.attribute.get('clientIdx', 1)
        client_idx_0 = (client_idx - 1) if client_idx >= 1 else 0

        # Find task class index
        tasks_matrix = layer.attribute.get('tasks', [])
        tidxclass = None
        if isinstance(tasks_matrix, np.ndarray) and len(tasks_matrix) > 0:
            for row in range(tasks_matrix.shape[0]):
                if tasks_matrix[row, 1] == tidx:
                    tidxclass = int(tasks_matrix[row, 0]) - 1
                    break
        elif isinstance(tasks_matrix, list):
            for row in tasks_matrix:
                if len(row) > 1 and row[1] == tidx:
                    tidxclass = row[0] - 1
                    break

        # Find entry class index
        entries_matrix = layer.attribute.get('entries', [])
        eidxclass = None
        if isinstance(entries_matrix, np.ndarray) and len(entries_matrix) > 0:
            for row in range(entries_matrix.shape[0]):
                if entries_matrix[row, 1] == eidx:
                    eidxclass = int(entries_matrix[row, 0]) - 1
                    break
        elif isinstance(entries_matrix, list):
            for row in entries_matrix:
                if len(row) > 1 and row[1] == eidx:
                    eidxclass = row[0] - 1
                    break

        task_tput = 0.0
        entry_tput = 0.0
        if tidxclass is not None and client_idx_0 < TN.shape[0] and tidxclass < TN.shape[1]:
            task_tput = TN[client_idx_0, tidxclass]
        if eidxclass is not None and client_idx_0 < TN.shape[0] and eidxclass < TN.shape[1]:
            entry_tput = TN[client_idx_0, eidxclass]

        return (task_tput, entry_tput)

    def _branch_members(self, joinaidx: int):
        """
        The activities belonging to each branch of an AND-join.

        A branch is recovered by walking backwards from each immediate predecessor of the
        join until an activity marked POST_AND is reached, that activity being the branch
        head spawned by the AND-fork. Branches between a fork and its join are disjoint
        paths, so the walk is unambiguous.

        actposttype, like the graph and residt, is indexed by global element index.
        """
        lqn = self.lqn
        graph = lqn.graph
        ashift = lqn.ashift
        nacts = lqn.nacts
        post_and_value = 12  # ActivityPrecedenceType.ID_POST_AND

        flat_posttype = lqn.actposttype.flatten() if getattr(lqn, 'actposttype', None) is not None \
            else np.zeros(0)

        members = []
        for tail in range(graph.shape[0]):
            if graph[tail, joinaidx] <= 0:
                continue
            if tail < ashift or tail >= ashift + nacts:
                continue  # not an activity
            chain = [tail]
            cur = tail
            for _ in range(nacts):
                if 0 < cur < len(flat_posttype) and flat_posttype[cur] == post_and_value:
                    break  # branch head
                prevs = [p for p in range(graph.shape[0])
                         if p != cur and graph[p, cur] > 0 and ashift <= p < ashift + nacts]
                if len(prevs) != 1:
                    break  # a merge or the start of the graph
                cur = prevs[0]
                chain.append(cur)
            members.append(chain)
        return members

    def _update_join_delays(self):
        """
        Compute the completion time of every AND-join and the correction it implies.

        The branches of an AND-fork run concurrently, so the time to pass the join is the
        k-th smallest of the branch completion times, k being the quorum of the join (k
        equals the branch count when the join waits for all its branches). Times are taken
        over residt because that is the quantity entry_servt aggregates.

        Returns a dict mapping each join target to the join time minus the sequential sum
        of its branch times, i.e. the amount by which the reachability matrix overcounts.
        """
        from line_solver.api.fj import quorum_moments

        lqn = self.lqn
        pre_and_value = 2  # ActivityPrecedenceType.ID_PRE_AND
        excess = {}
        if getattr(lqn, 'actpretype', None) is None:
            return excess
        flat_pretype = lqn.actpretype.flatten()
        flat_quorum = lqn.actquorum.flatten() if getattr(lqn, 'actquorum', None) is not None \
            else np.zeros(0)

        self.joint = np.zeros(lqn.nidx)
        # Restricted to activities: PRE_AND marks the branch tails, so the joins are the
        # activities whose predecessors carry that mark.
        for aidx in range(lqn.ashift, min(lqn.ashift + lqn.nacts, lqn.graph.shape[0])):
            # A join target is an activity whose predecessors are PRE_AND.
            is_join = False
            for pred in range(lqn.graph.shape[0]):
                if pred != aidx and lqn.graph[pred, aidx] > 0:
                    if 0 < pred < len(flat_pretype) and flat_pretype[pred] == pre_and_value:
                        is_join = True
                        break
            if not is_join:
                continue

            branches = self._branch_members(aidx)
            n = len(branches)
            if n == 0:
                continue
            # A branch time is residt PLUS the callresidt of every synchronous call
            # its activities issue: a branch activity with an Immediate host demand
            # does all its work in a rendezvous, and reading residt alone would make
            # this whole correction vanish silently. See _kb/06-solver-catalog.md.
            branch_times = []
            for b in branches:
                bt = float(sum(self.residt[m] for m in b))
                for baidx in b:
                    for cidx in self._calls_of(baidx):
                        if int(lqn.calltype[cidx]) == CallType.SYNC and cidx < len(self.callresidt):
                            bt += float(self.callresidt[cidx])
                branch_times.append(bt)
            if n == 1:
                self.joint[aidx] = branch_times[0]
                continue

            quorum = n
            if 0 < aidx < len(flat_quorum):
                q = int(flat_quorum[aidx])
                if 1 <= q <= n:
                    quorum = q
            # Branch times are taken as exponential, so the variance is the square of the mean.
            jt, _ = quorum_moments(branch_times, [t * t for t in branch_times], quorum)
            self.joint[aidx] = jt
            excess[aidx] = jt - sum(branch_times)
        return excess

    def _update_metrics_default(self, it: int):
        """Default metrics update (matches MATLAB updateMetricsDefault)."""
        lqn = self.lqn

        # Update activity service times from layer results
        self.servt = np.zeros(lqn.nidx)
        self.residt = np.zeros(lqn.nidx)

        # Calculate iter_min for averaging window (matches MATLAB updateMetricsDefault line 16)
        # MATLAB: iter_min = min(30, ceil(self.options.iter_max/4))
        iter_min = min(30, max(1, (self.options.iter_max + 3) // 4))  # Python ceil equivalent

        if self.servt_classes_updmap is not None:
            for r in range(len(self.servt_classes_updmap)):
                idx = int(self.servt_classes_updmap[r, 0])
                aidx = int(self.servt_classes_updmap[r, 1])
                nodeidx = int(self.servt_classes_updmap[r, 2])
                classidx = int(self.servt_classes_updmap[r, 3])

                layer_idx = int(self.idxhash[idx])

                # Convert 1-based indices to 0-based for numpy
                nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                classidx_0 = classidx - 1 if classidx >= 1 else 0

                # Compute refstat/refclass from layer struct for QN/TN_ref computation
                # (matches MATLAB updateMetricsDefault.m lines 23-27)
                refstat_k = None
                refclass_c = None
                if layer_idx < len(self.ensemble) and self.ensemble[layer_idx] is not None:
                    layer_sn = self.ensemble[layer_idx]._sn if hasattr(self.ensemble[layer_idx], '_sn') else None
                    if layer_sn is not None and hasattr(layer_sn, 'chains') and layer_sn.chains is not None:
                        chains_arr = np.asarray(layer_sn.chains)
                        if chains_arr.ndim == 2 and classidx_0 < chains_arr.shape[1]:
                            for ch in range(chains_arr.shape[0]):
                                if chains_arr[ch, classidx_0] > 0:
                                    if hasattr(layer_sn, 'refclass') and layer_sn.refclass is not None:
                                        rc = np.asarray(layer_sn.refclass).flatten()
                                        if ch < len(rc):
                                            refclass_c = int(rc[ch])
                                    break
                        if hasattr(layer_sn, 'refstat') and layer_sn.refstat is not None:
                            rs = np.asarray(layer_sn.refstat).flatten()
                            if classidx_0 < len(rs):
                                refstat_k = int(rs[classidx_0])

                # Apply averaging window for steady-state (matches MATLAB lines 18-31)
                if self.averagingstart is not None and it >= iter_min and len(self.results) > 1:
                    # Calculate window size (how many iterations since averaging started)
                    wnd_size = it - self.averagingstart + 1
                    if wnd_size > 1:
                        # Average over past iterations
                        self.servt[aidx] = 0.0
                        self.residt[aidx] = 0.0
                        self.tput[aidx] = 0.0
                        valid_samples = 0

                        for w in range(0, wnd_size):
                            result_idx = len(self.results) - 1 - w
                            if result_idx >= 0 and result_idx < len(self.results):
                                if layer_idx < len(self.results[result_idx]):
                                    hist_result = self.results[result_idx][layer_idx]
                                    if hist_result is not None and 'RN' in hist_result:
                                        hist_RN = hist_result['RN']
                                        hist_TN = hist_result['TN']
                                        hist_QN = hist_result.get('QN')
                                        hist_WN = hist_result.get('WN', hist_RN)
                                        if (hist_RN is not None and
                                            nodeidx_0 < hist_RN.shape[0] and
                                            classidx_0 < hist_RN.shape[1]):
                                            self.servt[aidx] += hist_RN[nodeidx_0, classidx_0]
                                            # Compute residt from QN/TN_ref (matches MATLAB)
                                            if (refstat_k is not None and refclass_c is not None and
                                                    hist_QN is not None and hist_TN is not None and
                                                    0 <= refstat_k < hist_TN.shape[0] and 0 <= refclass_c < hist_TN.shape[1]):
                                                TN_ref_w = hist_TN[refstat_k, refclass_c]
                                                if TN_ref_w > 1e-8:  # GlobalConstants.FineTol
                                                    self.residt[aidx] += hist_QN[nodeidx_0, classidx_0] / TN_ref_w
                                                else:
                                                    self.residt[aidx] += (hist_WN[nodeidx_0, classidx_0]
                                                                          if hist_WN is not None
                                                                          else hist_RN[nodeidx_0, classidx_0])
                                            else:
                                                self.residt[aidx] += (hist_WN[nodeidx_0, classidx_0]
                                                                      if hist_WN is not None
                                                                      else hist_RN[nodeidx_0, classidx_0])
                                            self.tput[aidx] += hist_TN[nodeidx_0, classidx_0]
                                            valid_samples += 1

                        if valid_samples > 0:
                            # MATLAB divides by wnd_size (not valid_samples) for damping effect
                            # (matches MATLAB updateMetricsDefault.m lines 22-25)
                            self.servt[aidx] /= wnd_size
                            self.residt[aidx] /= wnd_size
                            self.tput[aidx] /= wnd_size
                        else:
                            # No valid historical samples, fall through to latest result
                            self._extract_latest_metrics(aidx, layer_idx, nodeidx_0, classidx_0, refstat_k, refclass_c)
                    else:
                        # Window size is 1, use latest result
                        self._extract_latest_metrics(aidx, layer_idx, nodeidx_0, classidx_0, refstat_k, refclass_c)
                else:
                    # Before averaging starts, use latest result directly
                    self._extract_latest_metrics(aidx, layer_idx, nodeidx_0, classidx_0, refstat_k, refclass_c)

                # activity think time in series with host demand; see _kb/06-solver-catalog.md LN Activity think time section.
                zt_act = self._act_thinktime(aidx)
                if zt_act > 0:
                    self.servt[aidx] = self.servt[aidx] + zt_act
                    self.residt[aidx] = self.residt[aidx] + zt_act

                # Inf/NaN fallback: if layer MVA returned Inf/NaN, use previous value
                if it > 1:
                    if (np.isinf(self.servt[aidx]) or np.isnan(self.servt[aidx])) and not np.isnan(self.servt_prev[aidx]):
                        self.servt[aidx] = self.servt_prev[aidx]
                    if (np.isinf(self.residt[aidx]) or np.isnan(self.residt[aidx])) and not np.isnan(self.residt_prev[aidx]):
                        self.residt[aidx] = self.residt_prev[aidx]
                    if (np.isinf(self.tput[aidx]) or np.isnan(self.tput[aidx])) and not np.isnan(self.tput_prev[aidx]):
                        self.tput[aidx] = self.tput_prev[aidx]

                # Apply under-relaxation
                omega = self.relax_omega
                if omega < 1.0 and it > 1:
                    if not np.isnan(self.servt_prev[aidx]):
                        self.servt[aidx] = omega * self.servt[aidx] + (1 - omega) * self.servt_prev[aidx]
                    if not np.isnan(self.residt_prev[aidx]):
                        self.residt[aidx] = omega * self.residt[aidx] + (1 - omega) * self.residt_prev[aidx]
                    if not np.isnan(self.tput_prev[aidx]):
                        self.tput[aidx] = omega * self.tput[aidx] + (1 - omega) * self.tput_prev[aidx]

                self.servt_prev[aidx] = self.servt[aidx]
                self.residt_prev[aidx] = self.residt[aidx]
                self.tput_prev[aidx] = self.tput[aidx]

                # Update service time process with bounds checking
                # Safeguard against MVA numerical instability producing extreme values
                max_servt = 1e10
                if self.servt[aidx] > 0 and self.servt[aidx] <= max_servt:
                    self.servtproc[aidx] = Exp.fit_mean(self.servt[aidx])
                # mirrors MATLAB updateMetricsDefault.m:124; an async call's Source reads it.
                # Exp rejects rate 0 here where MATLAB admits it, so a null rate is Disabled.
                self.tputproc[aidx] = Exp.fit_rate(self.tput[aidx]) \
                    if self.tput[aidx] > 0 else Disabled()

                # async-only entries use RN for residt, not WN, as async arrivals don't share closed chain's visit ratio; mirrors MATLAB updateMetricsDefault.m:33-54.
                if lqn.ashift <= aidx < lqn.ashift + lqn.nacts:
                    # This is an activity - find its bound entry
                    for eidx in range(lqn.eshift, lqn.eshift + lqn.nentries):
                        # Check if activity is bound to this entry (edge from entry to activity in graph)
                        if hasattr(lqn, 'graph') and lqn.graph is not None:
                            if isinstance(lqn.graph, np.ndarray):
                                if eidx < lqn.graph.shape[0] and aidx < lqn.graph.shape[1]:
                                    if lqn.graph[eidx, aidx] > 0:
                                        # Found bound entry - check if async-only
                                        has_sync_callers = False
                                        has_async_callers = False

                                        if hasattr(lqn, 'issynccaller') and lqn.issynccaller is not None:
                                            if isinstance(lqn.issynccaller, np.ndarray):
                                                if eidx < lqn.issynccaller.shape[1]:
                                                    has_sync_callers = np.any(lqn.issynccaller[:, eidx])

                                        if hasattr(lqn, 'isasynccaller') and lqn.isasynccaller is not None:
                                            if isinstance(lqn.isasynccaller, np.ndarray):
                                                if eidx < lqn.isasynccaller.shape[1]:
                                                    has_async_callers = np.any(lqn.isasynccaller[:, eidx])

                                        if has_async_callers and not has_sync_callers:
                                            # Async-only target: use RN (response time per visit)
                                            # instead of WN (residence time with visit ratio)
                                            self.residt[aidx] = self.servt[aidx]  # servt already has RN
                                        break

        # throughput of activities that appear only as client-side classes, so that an
        # async call's Source has a rate; mirrors MATLAB updateMetricsDefault.m:161-183
        if self.thinkt_classes_updmap is not None:
            for r in range(len(self.thinkt_classes_updmap)):
                idx = int(self.thinkt_classes_updmap[r, 0])
                aidx = int(self.thinkt_classes_updmap[r, 1])
                nodeidx_0 = int(self.thinkt_classes_updmap[r, 2]) - 1
                classidx_0 = int(self.thinkt_classes_updmap[r, 3]) - 1
                if aidx >= len(self.tputproc) or self.tputproc[aidx] is not None:
                    continue
                if np.isnan(self.idxhash[idx]):
                    continue
                layer_idx = int(self.idxhash[idx])
                if layer_idx < 0 or not self.results or layer_idx >= len(self.results[-1]):
                    continue
                tp = 0.0
                wnd_size = (it - self.averagingstart + 1) if self.averagingstart is not None else 1
                if self.averagingstart is not None and it >= iter_min and wnd_size > 1 and len(self.results) > 1:
                    seen = 0
                    for w in range(0, wnd_size):
                        result_idx = len(self.results) - 1 - w
                        if result_idx < 0 or layer_idx >= len(self.results[result_idx]):
                            continue
                        hist = self.results[result_idx][layer_idx]
                        if hist is None or 'TN' not in hist or hist['TN'] is None:
                            continue
                        hist_TN = hist['TN']
                        if nodeidx_0 < hist_TN.shape[0] and classidx_0 < hist_TN.shape[1]:
                            tp += hist_TN[nodeidx_0, classidx_0]
                            seen += 1
                    tp = tp / wnd_size if seen > 0 else 0.0
                if tp == 0.0:
                    latest = self.results[-1][layer_idx]
                    if latest is not None and latest.get('TN') is not None:
                        TNl = latest['TN']
                        if nodeidx_0 < TNl.shape[0] and classidx_0 < TNl.shape[1]:
                            tp = TNl[nodeidx_0, classidx_0]
                self.tput[aidx] = tp
                self.tputproc[aidx] = Exp.fit_rate(tp) if tp > 0 else Disabled()

        # Update call service times (matches MATLAB updateMetricsDefault lines 140-162)
        self.callservt = np.zeros(lqn.ncalls)
        self.callresidt = np.zeros(lqn.ncalls)

        if self.call_classes_updmap is not None:
            for c in range(len(self.call_classes_updmap)):
                cidx = int(self.call_classes_updmap[c, 1])
                nodeidx = int(self.call_classes_updmap[c, 2])
                idx = int(self.call_classes_updmap[c, 0])
                classidx = int(self.call_classes_updmap[c, 3])

                # callresidt only updated for SERVER calls (nodeidx>1); CLIENT calls (Immediate) get callresidt from the callee's own SERVER-call entry.
                if nodeidx > 1:
                    layer_idx = int(self.idxhash[idx])
                    if layer_idx >= 0 and len(self.results) > 0 and layer_idx < len(self.results[-1]):
                        result = self.results[-1][layer_idx]
                        if result is not None and 'RN' in result:
                            RN = result['RN']
                            WN = result.get('WN', RN)
                            if RN is not None and WN is not None:
                                # Convert 1-based indices to 0-based for numpy
                                nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                                classidx_0 = classidx - 1 if classidx >= 1 else 0
                                if nodeidx_0 < RN.shape[0] and classidx_0 < RN.shape[1]:
                                    call_mean = self._get_call_mean(cidx)
                                    fcr_wait = self._region_wait(layer_idx, nodeidx_0, classidx_0, result)
                                    # MATLAB line 152: callservt = RN * callproc.getMean
                                    self.callservt[cidx] = (RN[nodeidx_0, classidx_0] + fcr_wait) * call_mean
                                    # Normalise per chain-reference visit, as residt does.
                                    # WN divides by the class's own reference rate when the
                                    # layer is open (an INF client task), which is per-ENTRY
                                    # visit, and the entry rescaling below would then count
                                    # the call once per entry.
                                    QNr = result.get('QN', None)
                                    TNr = result.get('TN', None)
                                    TN_ref = 0.0
                                    if QNr is not None and TNr is not None:
                                        refstat_k, refclass_c = self._chain_ref_indices(layer_idx, classidx_0)
                                        if refstat_k is not None and refclass_c is not None \
                                                and 0 <= refstat_k < TNr.shape[0] and 0 <= refclass_c < TNr.shape[1]:
                                            TN_ref = TNr[refstat_k, refclass_c]
                                    if TN_ref > GlobalConstants.FineTol:
                                        self.callresidt[cidx] = QNr[nodeidx_0, classidx_0] / TN_ref + fcr_wait
                                    else:
                                        self.callresidt[cidx] = WN[nodeidx_0, classidx_0] + fcr_wait

                                    # Inf/NaN fallback: if layer MVA returned Inf/NaN, use previous value
                                    if (np.isinf(self.callservt[cidx]) or np.isnan(self.callservt[cidx])) and it > 1 and not np.isnan(self.callservt_prev[cidx]):
                                        self.callservt[cidx] = self.callservt_prev[cidx]
                                    if (np.isinf(self.callresidt[cidx]) or np.isnan(self.callresidt[cidx])) and it > 1 and not np.isnan(self.callresidt_prev[cidx]):
                                        self.callresidt[cidx] = self.callresidt_prev[cidx]

                                    # Apply under-relaxation to call service times (MATLAB lines 155-160)
                                    omega = self.relax_omega
                                    if omega < 1.0 and it > 1 and not np.isnan(self.callservt_prev[cidx]):
                                        self.callservt[cidx] = omega * self.callservt[cidx] + (1 - omega) * self.callservt_prev[cidx]
                                    # no growth-rate capping on callservt: it would prevent convergence from near-zero (Immediate) initial values.

                                    self.callservt_prev[cidx] = self.callservt[cidx]
                                    self.callresidt_prev[cidx] = self.callresidt[cidx]

        # entry_servt = servtmatrix * [residt; callresidt]; servtmatrix carries cache hit/miss weighting probabilities.

        # servtmatrix indices 0..nidx are LQN elements, nidx+1..nidx+ncalls are calls.
        size = lqn.nidx + lqn.ncalls

        # Build combined vector
        combined_vec = np.zeros(size)

        # Fill activity residence times (indices are activity indices in LQN)
        for aidx in range(lqn.ashift, lqn.ashift + lqn.nacts):
            if self.residt[aidx] > 0:
                combined_vec[aidx] = self.residt[aidx]
            elif self.servtproc[aidx] is not None:
                proc = self.servtproc[aidx]
                if hasattr(proc, 'getMean'):
                    combined_vec[aidx] = proc.getMean()
                elif hasattr(proc, 'mean'):
                    combined_vec[aidx] = proc.mean

        # call residence times filled at index nidx+cidx; callresidt (=WN=RN*visits) already accounts for call_mean via the Aux class routing.
        for cidx in range(lqn.ncalls):
            combined_vec[lqn.nidx + cidx] = self.callresidt[cidx]

        # Compute entry service times: entry_servt = servtmatrix @ combined_vec
        entry_servt_vec = self.servtmatrix @ combined_vec
        entry_servt_vec[:lqn.eshift] = 0

        # FWD calls carry no blocking: callservt/callresidt stay zero since forwarding is handled by caller-side pseudo rendezvous.

        # Recompute entry_servt with forwarding-adjusted callresidt
        # (MATLAB updateMetricsDefault.m lines 349-351)
        for cidx in range(lqn.ncalls):
            combined_vec[lqn.nidx + cidx] = self.callresidt[cidx]
        entry_servt_vec = self.servtmatrix @ combined_vec
        entry_servt_vec[:lqn.eshift] = 0

        # AND-fork join correction; see _kb/06-solver-catalog.md LN Activity think time / AND-join concurrency section.
        joint_excess = self._update_join_delays()
        for eidx in range(lqn.eshift, lqn.eshift + lqn.nentries):
            corrected = entry_servt_vec[eidx]
            for aidx, exc in joint_excess.items():
                if exc != 0 and eidx < self.servtmatrix.shape[0] \
                        and aidx < self.servtmatrix.shape[1] \
                        and self.servtmatrix[eidx, aidx] > 0:
                    corrected += exc
            entry_servt_vec[eidx] = max(corrected, 0.0)

        # A SetupTask's cold start is charged HERE, to the entry, and with the
        # probability that the thread was actually found powered down. It is not
        # host demand, so it does not belong to any activity's residence:
        # reporting it there put RespT(A2) at 1.29479 on lqn_setup against the
        # 0.333178 LDES measures, which is the bare demand. See _setup_charge.
        for eidx in range(lqn.eshift, lqn.eshift + lqn.nentries):
            entry_servt_vec[eidx] += self._setup_charge(self._get_parent(eidx))

        # entry results scaled by throughput ratio so entries reach Ventry=1 while the task keeps Vtask=1; mirrors MATLAB lines 200-226.
        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            tidx = self._get_parent(eidx)  # task of entry
            hidx = self._get_parent(tidx) if tidx is not None else None  # host of entry

            if tidx is None or hidx is None:
                continue
            if self.ignore[tidx] or self.ignore[hidx]:
                continue

            entry_servt = entry_servt_vec[eidx]
            if entry_servt <= 0:
                continue

            # Check if this entry has sync callers (which create closed classes)
            # Use callpair and calltype to detect sync calls targeting this entry
            has_sync_callers = False
            if hasattr(lqn, 'callpair') and lqn.callpair is not None and hasattr(lqn, 'calltype'):
                for cidx in range(lqn.ncalls):
                    if cidx < len(lqn.callpair):
                        # callpair columns: [unused, src_aidx, tgt_eidx, mean_calls]
                        tgt_eidx = int(lqn.callpair[cidx, 1]) if lqn.callpair[cidx, 1] > 0 else 0
                        if tgt_eidx == eidx:
                            # Check if this is a SYNC call
                            calltype = lqn.calltype[cidx] if cidx < len(lqn.calltype) else 0
                            # Handle both integer and Enum values for calltype comparison
                            is_sync = (calltype == CallType.SYNC or
                                      calltype == CallType.SYNC.value or
                                      (isinstance(calltype, (int, np.integer)) and int(calltype) == CallType.SYNC.value))
                            if is_sync:
                                has_sync_callers = True
                                break

            if has_sync_callers:
                # Get throughput ratio from host layer results
                # This scales the entry service time by task_tput / entry_tput
                if not np.isnan(self.idxhash[hidx]):
                    layer_idx = int(self.idxhash[hidx])
                    if 0 <= layer_idx < len(self.ensemble) and self.ensemble[layer_idx] is not None:
                        layer = self.ensemble[layer_idx]
                        result = self.results[-1][layer_idx] if len(self.results) > 0 and layer_idx < len(self.results[-1]) else None

                        if result is not None and 'TN' in result:
                            TN = result['TN']
                            client_idx = layer.attribute.get('clientIdx', 1)
                            client_idx_0 = (client_idx - 1) if client_idx >= 1 else 0

                            # Find ALL task class indices (MATLAB: find(...==tidx) returns vector)
                            tasks_matrix = layer.attribute.get('tasks', [])
                            tidxclasses = []
                            if isinstance(tasks_matrix, np.ndarray) and len(tasks_matrix) > 0:
                                for row in range(tasks_matrix.shape[0]):
                                    if tasks_matrix[row, 1] == tidx:
                                        tidxclasses.append(int(tasks_matrix[row, 0]) - 1)
                            elif isinstance(tasks_matrix, list):
                                for row in tasks_matrix:
                                    if len(row) > 1 and row[1] == tidx:
                                        tidxclasses.append(row[0] - 1)

                            # Find ALL entry class indices
                            entries_matrix = layer.attribute.get('entries', [])
                            eidxclasses = []
                            if isinstance(entries_matrix, np.ndarray) and len(entries_matrix) > 0:
                                for row in range(entries_matrix.shape[0]):
                                    if entries_matrix[row, 1] == eidx:
                                        eidxclasses.append(int(entries_matrix[row, 0]) - 1)
                            elif isinstance(entries_matrix, list):
                                for row in entries_matrix:
                                    if len(row) > 1 and row[1] == eidx:
                                        eidxclasses.append(row[0] - 1)

                            # Compute throughput ratio (sum over all matching classes)
                            task_tput = 0.0
                            entry_tput = 0.0

                            for tc in tidxclasses:
                                if client_idx_0 < TN.shape[0] and tc < TN.shape[1]:
                                    task_tput += TN[client_idx_0, tc]
                            for ec in eidxclasses:
                                if client_idx_0 < TN.shape[0] and ec < TN.shape[1]:
                                    entry_tput += TN[client_idx_0, ec]

                            # scale entry service time by task_tput/entry_tput (a task may call several entries per cycle); mirrors MATLAB updateMetricsDefault.m:217.
                            if entry_tput > GlobalConstants.Zero:
                                self.servt[eidx] = entry_servt * task_tput / entry_tput
                                self.residt[eidx] = entry_servt * task_tput / entry_tput
                            else:
                                self.servt[eidx] = entry_servt
                                self.residt[eidx] = entry_servt
                        else:
                            # No results yet, use unscaled entry_servt
                            self.servt[eidx] = entry_servt
                            self.residt[eidx] = entry_servt
                    else:
                        self.servt[eidx] = entry_servt
                        self.residt[eidx] = entry_servt
                else:
                    self.servt[eidx] = entry_servt
                    self.residt[eidx] = entry_servt
            else:
                # For async-only targets, use entry_servt directly
                # No throughput ratio scaling needed since there are no closed classes
                self.servt[eidx] = entry_servt
                self.residt[eidx] = entry_servt

        # Phase-2 support: split activity service times by phase and apply correction
        # Matches MATLAB updateMetricsDefault.m lines 105-330
        if self.hasPhase2:
            # Reset phase-specific arrays
            self.servt_ph1 = np.zeros(lqn.nidx)
            self.servt_ph2 = np.zeros(lqn.nidx)

            # Split activity service times by phase
            for a in range(lqn.nacts):
                aidx = lqn.ashift + a
                if lqn.actphase[a - 1] == 1:  # actphase is 0-indexed numpy array
                    self.servt_ph1[aidx] = self.servt[aidx]
                else:
                    self.servt_ph2[aidx] = self.servt[aidx]

            # Aggregate phase service times to entry level
            for e in range(lqn.nentries):
                eidx = lqn.eshift + e
                acts = lqn.actsof.get(eidx, [])
                for aidx in acts:
                    a = aidx - lqn.ashift
                    if 1 <= a <= lqn.nacts:
                        if lqn.actphase[a - 1] == 1:
                            self.servt_ph1[eidx] += self.servt_ph1[aidx]
                        else:
                            self.servt_ph2[eidx] += self.servt_ph2[aidx]

            # overtaking probability response-time correction; see _kb/06-solver-catalog.md LN phase-2 overtaking section.
            for e in range(lqn.nentries):
                eidx = lqn.eshift + e
                if self.servt_ph2[eidx] > 1e-8:  # GlobalConstants.FineTol
                    tidx = self._get_parent(eidx)

                    # REF tasks and entries without sync callers see the full service time
                    if (tidx is not None and self._is_ref_task(tidx)) \
                            or not self._has_sync_callers_for_entry(eidx):
                        self.residt[eidx] = self.servt[eidx]
                        continue

                    # Get entry throughput
                    if self.tput[eidx] > 1e-8:
                        entry_tput = self.tput[eidx]
                    elif tidx is not None and self.tput[tidx] > 1e-8:
                        entry_tput = self.tput[tidx]
                    else:
                        entry_tput = 0

                    # Compute overtaking probability
                    if entry_tput > 1e-8:
                        self.prOvertake[e] = self._overtake_prob(eidx)
                    else:
                        self.prOvertake[e] = 0

                    # Caller's response time = phase-1 + P(overtake) * phase-2
                    overtake_delay = self.prOvertake[e] * self.servt_ph2[eidx]
                    self.residt[eidx] = self.servt_ph1[eidx] + overtake_delay

        # servtproc for entries updated before callservtproc (which reads it); mirrors MATLAB updateMetricsDefault.m:265-271.
        if self.call_classes_updmap is not None and len(self.call_classes_updmap) > 0:
            for row in self.call_classes_updmap:
                cidx = int(row[1])
                nodeidx = int(row[2])

                # Get serverIdx for this layer
                idx = int(row[0])
                layer_idx = int(self.idxhash[idx]) if not np.isnan(self.idxhash[idx]) else -1
                if layer_idx < 0 or layer_idx >= len(self.ensemble):
                    continue
                layer = self.ensemble[layer_idx]
                if layer is None:
                    continue
                # any non-client node is a server station; under flat layering the callee
                # station is not the layer's serverIdx, so test nodeidx > 1 as MATLAB does.
                if nodeidx > 1:
                    eidx = self._get_call_target_entry(cidx)
                    if eidx is not None and eidx > 0 and eidx < len(self.servt):
                        if self.servt[eidx] > 0:
                            self.servtproc[eidx] = Exp.fit_mean(self.servt[eidx])

        # callservtproc only for SERVER calls (nodeidx>1); CLIENT calls stay Immediate, resp via think-time; mirrors MATLAB updateMetricsDefault.m:274-287.
        if self.call_classes_updmap is not None and len(self.call_classes_updmap) > 0:
            for row in self.call_classes_updmap:
                idx = int(row[0])
                cidx = int(row[1])
                nodeidx = int(row[2])
                classidx = int(row[3])

                # Get serverIdx for this layer to check if call is at server
                layer_idx = int(self.idxhash[idx]) if not np.isnan(self.idxhash[idx]) else -1
                if layer_idx < 0 or layer_idx >= len(self.ensemble):
                    continue
                layer = self.ensemble[layer_idx]
                if layer is None:
                    continue
                # only SERVER-node calls update callservtproc, never CLIENT-node calls; mirrors MATLAB line 277.
                if nodeidx > 1:
                    eidx = self._get_call_target_entry(cidx)
                    if eidx is not None and eidx > 0:
                        if it == 1:
                            # first iteration uses servtproc[eidx] (Immediate for entries, non-zero only via bound activities); mirrors MATLAB line 281.
                            if eidx < len(self.servt):
                                self.callservt[cidx] = self.servt[eidx]
                            # Use servtproc for callservtproc (matches MATLAB: callservtproc{cidx} = servtproc{eidx})
                            if eidx < len(self.servtproc) and self.servtproc[eidx] is not None:
                                self.callservtproc[cidx] = self.servtproc[eidx]
                        else:
                            # Subsequent iterations: use callservt from layer results
                            if self.callservt[cidx] > 0:
                                self.callservtproc[cidx] = Exp.fit_mean(self.callservt[cidx])

        # Compute ptaskcallers - probability that request to task/host comes from caller
        self._compute_ptaskcallers()

    def _overtake_prob(self, eidx):
        """Compute overtaking probability using 3-state CTMC.

        Matches MATLAB overtake_prob.m.
        States: 0=idle, 1=phase-1, 2=phase-2
        By PASTA, P(overtake) = steady-state prob of being in phase-2.
        """
        lqn = self.lqn

        S1 = self.servt_ph1[eidx]
        S2 = self.servt_ph2[eidx]

        # Get throughput
        tidx = self._get_parent(eidx)
        if self.tput[eidx] > 1e-8:
            lam = self.tput[eidx]
        elif tidx is not None and self.tput[tidx] > 1e-8:
            lam = self.tput[tidx]
        else:
            return 0.0

        # Number of servers (multiplicity of parent task)
        c = 1
        if tidx is not None and hasattr(lqn, 'mult') and lqn.mult is not None:
            if isinstance(lqn.mult, (dict,)):
                c = int(lqn.mult.get(tidx, 1))
            elif isinstance(lqn.mult, np.ndarray) and tidx < len(lqn.mult):
                c = int(lqn.mult[tidx])

        # Degenerate cases
        if S2 < 1e-8 or lam < 1e-8 or S1 < 1e-8:
            return 0.0

        mu1 = 1.0 / S1
        mu2 = 1.0 / S2

        if c == 1:
            # single-server exact CTMC steady-state via the augmented balance system [Q';1 1 1]*pi=[0;0;0;1].
            A = np.array([
                [-lam, 0, mu2, 1],
                [lam, -mu1, 0, 1],
                [0, mu1, -mu2, 1]
            ]).T  # 4x3
            b = np.array([0, 0, 0, 1])
            # Least squares solve
            pi, _, _, _ = np.linalg.lstsq(A, b, rcond=None)
            return max(0.0, min(1.0, pi[2]))
        else:
            # Multi-server approximation
            rho = lam * (S1 + S2) / c
            if rho >= 1:
                return S2 / (S1 + S2)
            else:
                return max(0.0, min(1.0, (S2 / (S1 + S2)) * rho))

    def _compute_ptaskcallers(self):
        """
        Compute caller probability matrices for interlocking correction.

        This implements the MATLAB updateMetricsDefault ptaskcallers computation:
        1. Compute direct caller probabilities from throughputs
        2. Compute indirect caller probabilities via DTMC random walk
        """
        lqn = self.lqn

        # Reset ptaskcallers
        self.ptaskcallers = np.zeros((lqn.nhosts + lqn.ntasks, lqn.nhosts + lqn.ntasks))

        # Compute direct caller probabilities for tasks
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self._is_ref_task(tidx):
                continue

            # Get callers of this task (via iscaller matrix)
            callers = self._get_callers_of_task(tidx)
            if not callers:
                continue

            # Get throughput of each caller from task layer results
            caller_tput = np.zeros(lqn.ntasks)

            if np.isnan(self.idxhash[tidx]):
                continue

            layer_idx = int(self.idxhash[tidx])
            if layer_idx < 0 or layer_idx >= len(self.ensemble):
                continue

            layer = self.ensemble[layer_idx]
            if layer is None or len(self.results) == 0:
                continue

            result = self.results[-1][layer_idx] if layer_idx < len(self.results[-1]) else None
            if result is None or 'TN' not in result:
                continue

            TN = result['TN']
            client_idx = layer.attribute.get('clientIdx', 1)
            if client_idx is None:
                continue
            client_idx_0 = client_idx - 1 if client_idx >= 1 else 0

            tasks_matrix = layer.attribute.get('tasks', [])

            for caller_idx in callers:
                # Find class index for this caller in the layer
                caller_class_idx = None
                if isinstance(tasks_matrix, np.ndarray) and len(tasks_matrix) > 0:
                    for row in range(tasks_matrix.shape[0]):
                        if tasks_matrix[row, 1] == caller_idx:
                            caller_class_idx = int(tasks_matrix[row, 0])
                            break
                elif isinstance(tasks_matrix, list):
                    for row in tasks_matrix:
                        if len(row) > 1 and row[1] == caller_idx:
                            caller_class_idx = row[0]
                            break

                if caller_class_idx is not None:
                    caller_class_idx_0 = caller_class_idx - 1 if caller_class_idx >= 1 else 0
                    if client_idx_0 < TN.shape[0] and caller_class_idx_0 < TN.shape[1]:
                        caller_tput[caller_idx - lqn.tshift] = TN[client_idx_0, caller_class_idx_0]

            # Normalize to get probabilities
            total_tput = np.sum(caller_tput)
            if total_tput > GlobalConstants.Zero:
                self.ptaskcallers[tidx, lqn.tshift:lqn.tshift + lqn.ntasks] = caller_tput / total_tput

        # Compute direct caller probabilities for hosts
        for hidx in range(lqn.nhosts):
            if np.isnan(self.idxhash[hidx]):
                continue

            layer_idx = int(self.idxhash[hidx])
            if layer_idx < 0 or layer_idx >= len(self.ensemble):
                continue

            layer = self.ensemble[layer_idx]
            if layer is None or len(self.results) == 0:
                continue

            result = self.results[-1][layer_idx] if layer_idx < len(self.results[-1]) else None
            if result is None or 'TN' not in result:
                continue

            TN = result['TN']
            client_idx = layer.attribute.get('clientIdx', 1)
            if client_idx is None:
                continue
            client_idx_0 = client_idx - 1 if client_idx >= 1 else 0

            callers = self._get_tasks_of_host(hidx)
            tasks_matrix = layer.attribute.get('tasks', [])

            caller_tput = np.zeros(lqn.ntasks)
            for caller_idx in callers:
                # Find class index for this caller in the layer
                caller_class_idx = None
                if isinstance(tasks_matrix, np.ndarray) and len(tasks_matrix) > 0:
                    for row in range(tasks_matrix.shape[0]):
                        if tasks_matrix[row, 1] == caller_idx:
                            caller_class_idx = int(tasks_matrix[row, 0])
                            break
                elif isinstance(tasks_matrix, list):
                    for row in tasks_matrix:
                        if len(row) > 1 and row[1] == caller_idx:
                            caller_class_idx = row[0]
                            break

                if caller_class_idx is not None:
                    caller_class_idx_0 = caller_class_idx - 1 if caller_class_idx >= 1 else 0
                    if client_idx_0 < TN.shape[0] and caller_class_idx_0 < TN.shape[1]:
                        caller_tput[caller_idx - lqn.tshift] += TN[client_idx_0, caller_class_idx_0]

            # Normalize to get probabilities
            total_tput = np.sum(caller_tput)
            if total_tput > GlobalConstants.Zero:
                self.ptaskcallers[hidx, lqn.tshift:lqn.tshift + lqn.ntasks] = caller_tput / total_tput

        # Compute ptaskcallers_step using DTMC random walk
        P = self.ptaskcallers.copy()

        # Make stochastic: rows that sum to 0 should self-loop
        row_sums = P.sum(axis=1)
        for i in range(P.shape[0]):
            if row_sums[i] < GlobalConstants.FineTol:
                P[i, i] = 1.0  # Self-loop at absorbing states

        self.ptaskcallers_step[0] = P.copy()

        # Walk backward through caller graph
        for hidx in range(lqn.nhosts):
            if np.isnan(self.idxhash[hidx]):
                continue

            callers = self._get_tasks_of_host(hidx)
            for tidx in callers:
                # Initialize probability mass at host
                x0 = np.zeros(len(self.ptaskcallers))
                x0[hidx] = 1.0

                x = x0 @ P  # First step

                for step in range(1, self.nlayers + 1):
                    x = x @ P

                    if step < len(self.ptaskcallers_step):
                        self.ptaskcallers_step[step][tidx, :] = x
                        # Weight by caller probability for host
                        self.ptaskcallers_step[step][hidx, :] = self.ptaskcallers[hidx, tidx] * x

                    # Check if all probability reached REF tasks
                    ref_prob = 0.0
                    for t in range(lqn.ntasks):
                        t_idx = lqn.tshift + t
                        if self._is_ref_task(t_idx):
                            ref_prob += x[t_idx]

                    if ref_prob > 1.0 - self.options.tol:
                        break

                    # Update max callers
                    self.ptaskcallers[:, tidx] = np.maximum(self.ptaskcallers[:, tidx], x)

    def _get_call_mean(self, cidx: int) -> float:
        """Get mean number of calls."""
        lqn = self.lqn

        # First try callproc (distribution)
        if hasattr(lqn, 'callproc') and lqn.callproc is not None:
            if isinstance(lqn.callproc, dict):
                proc = lqn.callproc.get(cidx)
            elif isinstance(lqn.callproc, (list, np.ndarray)):
                if cidx < len(lqn.callproc):
                    proc = lqn.callproc[cidx]  # 1-indexed
                else:
                    proc = None
            else:
                proc = None

            if proc is not None:
                if hasattr(proc, 'getMean'):
                    return proc.getMean()
                elif hasattr(proc, 'mean'):
                    return proc.mean

        # Fallback to callpair column 3 (call mean)
        if hasattr(lqn, 'callpair') and lqn.callpair is not None:
            if isinstance(lqn.callpair, np.ndarray):
                if cidx < lqn.callpair.shape[0]:
                    return float(lqn.callpair[cidx, 2])

        return 1.0

    def _get_call_response_time(self, caller_tidx: int) -> float:
        """
        Get the total call response time for a task.

        This is the sum of (call_mean * callee_response_time) for all calls
        made by activities of this task.
        """
        lqn = self.lqn
        total_call_time = 0.0

        # Find all calls from this task's activities
        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return 0.0

        for cidx in range(lqn.ncalls):
            if cidx >= lqn.callpair.shape[0]:
                continue

            # Get source activity (column 1) and target entry (column 2)
            src_aidx = int(lqn.callpair[cidx, 0])
            tgt_eidx = int(lqn.callpair[cidx, 1])
            call_mean = float(lqn.callpair[cidx, 2]) if lqn.callpair.shape[1] > 2 else 1.0

            if src_aidx == 0 or tgt_eidx == 0:
                continue

            # Get parent task of source activity
            src_tidx = self._get_parent(src_aidx)
            if src_tidx != caller_tidx:
                continue

            # Get parent task of target entry
            tgt_tidx = self._get_parent(tgt_eidx)
            if tgt_tidx is None:
                continue

            # call response time = callee task layer response time, from callservtproc when available.
            if cidx < len(self.callservtproc) and self.callservtproc[cidx] is not None:
                proc = self.callservtproc[cidx]
                if hasattr(proc, 'getMean'):
                    total_call_time += call_mean * proc.getMean()
                elif hasattr(proc, 'mean'):
                    total_call_time += call_mean * proc.mean
                continue

            # Fall back to task layer results if callservtproc not available
            if not np.isnan(self.idxhash[tgt_tidx]):
                tgt_layer_idx = int(self.idxhash[tgt_tidx])
                if len(self.results) > 0 and tgt_layer_idx < len(self.results[-1]):
                    result = self.results[-1][tgt_layer_idx]
                    if result is not None and 'RN' in result:
                        RN = result['RN']
                        server_idx = self.ensemble[tgt_layer_idx].attribute.get('serverIdx', 1)
                        if server_idx is not None:
                            server_idx_0 = server_idx - 1 if server_idx >= 1 else 0
                            if server_idx_0 < RN.shape[0]:
                                # Find the caller's activity class (not task class) in this layer
                                caller_class_idx = self._find_activity_class_in_layer(src_aidx, tgt_layer_idx)
                                if caller_class_idx is None:
                                    # Fallback to task class
                                    caller_class_idx = self._find_caller_class_in_layer(caller_tidx, tgt_layer_idx)
                                if caller_class_idx is not None:
                                    caller_class_idx_0 = caller_class_idx - 1 if caller_class_idx >= 1 else 0
                                    if caller_class_idx_0 < RN.shape[1]:
                                        callee_resp = RN[server_idx_0, caller_class_idx_0]
                                        total_call_time += call_mean * callee_resp
                                        continue
                                # Fallback: average across all classes
                                callee_resp = np.mean(RN[server_idx_0, :])
                                total_call_time += call_mean * callee_resp
                                continue

            # Fallback: use entry's service time if task layer not available
            if tgt_eidx < len(self.servt) and self.servt[tgt_eidx] > 0:
                total_call_time += call_mean * self.servt[tgt_eidx]

        return total_call_time

    def _get_throughput_from_callers(self, tidx: int) -> float:
        """
        Compute task throughput from callers' rates.

        For a purely called task T, throughput = sum of (caller_tput * call_mean)
        for all calls that target entries of T.
        """
        lqn = self.lqn
        total_tput = 0.0

        # Get entries of this task
        entries = self._get_entries_of_task(tidx)
        if not entries:
            return 0.0

        # For each call, check if it targets one of our entries
        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return 0.0

        for cidx in range(lqn.ncalls):
            if cidx >= lqn.callpair.shape[0]:
                continue

            tgt_eidx = int(lqn.callpair[cidx, 1])  # Target entry
            if tgt_eidx not in entries:
                continue

            # This call targets our task - get caller's throughput
            src_aidx = int(lqn.callpair[cidx, 0])  # Source activity
            if src_aidx <= 0:
                continue

            # Get task of source activity
            caller_tidx = self._get_parent(src_aidx)
            if caller_tidx is None or caller_tidx <= 0:
                continue

            # call rate = caller ACTIVITY throughput * call_mean (the call fires each activity execution, not each task cycle).
            caller_tput = self.tput[src_aidx] if src_aidx < len(self.tput) else 0.0
            call_mean = self._get_call_mean(cidx)

            total_tput += caller_tput * call_mean

        return total_tput

    def _setup_dist_mean(self, procs, tidx: int) -> float:
        """Mean of a setup or delay-off process of task TIDX, 0 when it declares none."""
        if procs is None:
            return 0.0
        p = None
        if isinstance(procs, dict):
            p = procs.get(tidx)
        else:
            arr = np.asarray(procs).flatten()
            if tidx < len(arr):
                p = arr[tidx]
        if p is None:
            return 0.0
        try:
            m = float(p.getMean())
        except (AttributeError, TypeError, ValueError):
            return 0.0
        return 0.0 if (np.isnan(m) or np.isinf(m)) else m

    def _setup_charge(self, tidx: int) -> float:
        """Mean cold start one request of task TIDX pays, 0 when it declares none.

        A SetupTask powers a thread down when it goes idle and pays a setup before
        it can serve again. The thread is released at a reply and starts a delay-off
        countdown D of mean d; it powers off only if D expires before the next
        request arrives, and a request arriving first cancels the countdown and pays
        nothing. With the idle interval I seen by one thread and exponential D,

            p = P(D < I) = E[I] / (E[I] + d),   and the charge is  p * s.

        E[I] comes from the current iterate. Admission takes an ACTIVE idle thread
        before it wakes a sleeping one, so the pool that actually cycles is only as
        large as the load needs: with offered load b = X*S = rho*mult threads, about
        max(1,b) stay hot, each seeing arrivals at rate X/max(1,b) and busy S per
        arrival, so E[I] = (max(1,b) - b) / X. At mult = 1 this is (1-rho)/X and is
        EXACT given p, returning p = a/(a+d) for the one-customer model the LDES
        engine is checked against. Above one thread it is an approximation, the
        exact answer for c servers with setup being matrix-analytic (Gandhi,
        Harchol-Balter and Adan, Performance Evaluation 67(11), 2010). Twin of
        MATLAB lqn_setup_charge.m and the JAR SolverLN.setupCharge.
        """
        lqn = self.lqn
        hs = getattr(lqn, 'hassetup', None)
        if hs is None:
            return 0.0
        hsf = np.asarray(hs).flatten()
        if tidx < 0 or tidx >= len(hsf) or not hsf[tidx]:
            return 0.0
        s = self._setup_dist_mean(getattr(lqn, 'setuptime', None), tidx)
        d = self._setup_dist_mean(getattr(lqn, 'delayofftime', None), tidx)
        if not (s > GlobalConstants.FineTol) or not (d > GlobalConstants.FineTol):
            return 0.0
        mult = float(lqn.mult[0, tidx])
        if not np.isfinite(mult) or mult <= 0:
            return 0.0  # an infinite-server task holds no thread to power down
        if self.tput is None or self.util is None or tidx >= len(self.tput) or tidx >= len(self.util):
            return s  # nothing has arrived yet, so the thread is down when the first does
        X = float(self.tput[tidx])
        if not np.isfinite(X) or X <= GlobalConstants.FineTol:
            return s
        rho = float(self.util[tidx])
        if not np.isfinite(rho) or rho < 0:
            rho = 0.0
        rho = min(rho, 1 - GlobalConstants.FineTol)
        b = rho * mult                      # offered load, in threads
        EI = (max(1.0, b) - b) / X          # idle interval of a thread in the hot pool
        return s * EI / (EI + d)

    def update_think_times(self, it: int):
        """Update think times (matches MATLAB updateThinkTimes)."""
        # Under 'srvn.ph' a caller reaches the server once per invocation, so the
        # station rate is not the task's invocation rate -- see the PH twin
        if self._is_ph_encoding():
            self._update_think_times_ph(it)
            return
        lqn = self.lqn

        if not hasattr(lqn, 'iscaller') or lqn.iscaller is None:
            return

        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t

            # Only a REFERENCE task's think time separates one request from the
            # next; on a served task it is not a per-request delay and charging
            # it throttles the task -- see _ref_think_mean
            tidx_thinktime = self._ref_think_mean(tidx)

            # Get call response time (time spent waiting for calls to complete)
            call_response_time = self._get_call_response_time(tidx)

            if not np.isnan(self.idxhash[tidx]):
                # Get throughput and utilization from layer results
                layer_idx = int(self.idxhash[tidx])
                njobs = max(self.njobs[tidx, :])

                if len(self.results) > 0 and layer_idx < len(self.results[-1]):
                    result = self.results[-1][layer_idx]
                    if result is not None and 'TN' in result:
                        TN = result['TN']
                        UN = result.get('UN', np.zeros_like(TN))

                        # the station of TIDX inside its layer: under flat layering
                        # every server shares one layer, so the scalar serverIdx
                        # would point at the first processor instead
                        server_idx = self._station_idx_of(self.ensemble[layer_idx], tidx)
                        if server_idx is not None:
                            # Convert to 0-based indexing for numpy
                            server_idx_0 = server_idx - 1 if server_idx >= 1 else 0

                            # task throughput read from its own layer's SERVER node, repl-scaled; mirrors MATLAB updateThinkTimes line 24 (no fork_fanout correction there).
                            repl = lqn.repl[0, tidx] if hasattr(lqn, 'repl') and lqn.repl is not None and lqn.repl.shape[1] > tidx else 1.0
                            tput_sum = np.nansum(TN[server_idx_0, :])
                            self.tput[tidx] = repl * tput_sum if np.isfinite(tput_sum) else 0.0
                            util_sum = np.nansum(UN[server_idx_0, :])
                            self.util[tidx] = util_sum if np.isfinite(util_sum) else 0.0

                            # Compute think time - MATLAB formula only uses user think time
                            # Call response time is handled separately in update_layers
                            sched = self._get_sched(tidx)

                            if sched == SchedStrategy.INF:
                                # Infinite server case
                                if self.tput[tidx] > GlobalConstants.Zero:
                                    self.thinkt[tidx] = max(GlobalConstants.Zero,
                                                    (njobs - self.util[tidx]) / self.tput[tidx] - tidx_thinktime)
                            else:
                                # Regular queue case
                                if self.tput[tidx] > GlobalConstants.Zero:
                                    self.thinkt[tidx] = max(GlobalConstants.Zero,
                                                           njobs * abs(1 - self.util[tidx]) / self.tput[tidx] - tidx_thinktime)

                # A caller class cycles as delay plus station service, and the station
                # serves only the host demand: a cold start is charged to the entry,
                # not to any activity's demand, so the station never sees it and the
                # delay has to carry it. Without this the callee layer cycled at
                # 0.529412 against the 0.5 its callers drive on lqn_setup. Zero for
                # every task without a setup.
                self.thinkt[tidx] = max(GlobalConstants.Zero,
                                        self.thinkt[tidx] + self._setup_charge(tidx))

                # Recover from Inf/NaN: snap back to previous iteration's value
                if it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                    if np.isinf(self.thinkt[tidx]) or np.isnan(self.thinkt[tidx]):
                        self.thinkt[tidx] = self.thinkt_prev[tidx]
                # Apply under-relaxation
                omega = self.relax_omega
                if omega < 1.0 and it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                    rawT = self.thinkt[tidx]
                    prevT = self.thinkt_prev[tidx]
                    # If recovering from crash (prev much larger than raw), snap to raw
                    if prevT > 10 * rawT and rawT > GlobalConstants.FineTol:
                        self.thinkt_prev[tidx] = rawT  # reset prev to allow recovery
                    self.thinkt[tidx] = omega * self.thinkt[tidx] + (1 - omega) * self.thinkt_prev[tidx]

                self.thinkt_prev[tidx] = self.thinkt[tidx]

                # Update think time process - MATLAB style: thinkt + user_think only
                # Call response time will be added in update_layers
                if self.thinkt[tidx] + tidx_thinktime > 0:
                    self.thinktproc[tidx] = Exp.fit_mean(self.thinkt[tidx] + tidx_thinktime)
                else:
                    self.thinktproc[tidx] = Immediate()

                # for non-REF called tasks, task-layer throughput takes precedence over host-layer (INF-scheduled hosts overstate it).
                if self.tput[tidx] == 0 or np.isnan(self.tput[tidx]):
                    hidx = self._get_parent(tidx)  # host = parent of task
                    if hidx is not None and not np.isnan(self.idxhash[hidx]):
                        host_layer_idx = int(self.idxhash[hidx])
                        if host_layer_idx >= 0 and len(self.results) > 0 and host_layer_idx < len(self.results[-1]):
                            result = self.results[-1][host_layer_idx]
                            if result is not None and 'TN' in result:
                                TN = result['TN']
                                # For called tasks, find the server node and task class
                                if self.servt_classes_updmap is not None:
                                    for r in range(len(self.servt_classes_updmap)):
                                        if int(self.servt_classes_updmap[r, 0]) == hidx:
                                            # Check if this mapping is for our task's activity
                                            aidx = int(self.servt_classes_updmap[r, 1])
                                            if self._get_parent(aidx) == tidx:
                                                nodeidx = int(self.servt_classes_updmap[r, 2])
                                                classidx = int(self.servt_classes_updmap[r, 3])
                                                nodeidx_0 = nodeidx - 1 if nodeidx >= 1 else 0
                                                classidx_0 = classidx - 1 if classidx >= 1 else 0
                                                if nodeidx_0 < TN.shape[0] and classidx_0 < TN.shape[1]:
                                                    # Only set if not already set by task layer
                                                    if self.tput[tidx] == 0 or np.isnan(self.tput[tidx]):
                                                        tn_val = TN[nodeidx_0, classidx_0]
                                                        if np.isfinite(tn_val):
                                                            self.tput[tidx] = tn_val
                                                break
            else:
                # Ref task, forwarding target or open-arrival target (no task layer).
                # An entry arrival that is the only way in drives the thread pool
                # directly: the layer builder dropped its open class precisely so the
                # cycle can be closed on the known rate here. See _open_arrival_rate_of.
                arvrate = self._open_arrival_rate_of(tidx)
                if arvrate > GlobalConstants.FineTol:
                    njobs_arv = max(self.njobs[tidx, :])
                    if not njobs_arv > 0:
                        njobs_arv = lqn.maxmult[tidx]
                    self.tput[tidx] = arvrate
                    host_residt = 0.0
                    for eidx_arv in self._get_entries_of_task(tidx):
                        if hasattr(lqn, 'actsof') and eidx_arv in lqn.actsof:
                            for aidx_arv in lqn.actsof[eidx_arv]:
                                if np.isfinite(self.residt[aidx_arv]):
                                    host_residt += self.residt[aidx_arv]
                    z_arv = max(GlobalConstants.Zero,
                                njobs_arv / arvrate - host_residt - tidx_thinktime)
                    omega = self.relax_omega
                    if omega < 1.0 and it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                        z_arv = omega * z_arv + (1 - omega) * self.thinkt_prev[tidx]
                    self.thinkt[tidx] = z_arv
                    self.thinkt_prev[tidx] = z_arv
                    self.thinktproc[tidx] = Exp.fit_mean(z_arv + tidx_thinktime)
                    continue

                # Check if this is a forwarding target task (MATLAB updateThinkTimes.m:54-104)
                is_fwd_target = False
                fwd_cidx_found = None
                source_tidx = None
                fwd_prob = 0.0
                if not self._is_ref_task(tidx) and hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    for eidx_fwd in self._get_entries_of_task(tidx):
                        for cidx_fwd in range(lqn.ncalls):
                            if cidx_fwd < len(lqn.calltype) and int(lqn.calltype[cidx_fwd]) == CallType.FWD \
                                    and int(lqn.callpair[cidx_fwd, 1]) == eidx_fwd:
                                is_fwd_target = True
                                fwd_cidx_found = cidx_fwd
                                source_eidx = int(lqn.callpair[cidx_fwd, 0])
                                source_tidx = self._get_parent(source_eidx)
                                fwd_prob = self._get_call_mean(cidx_fwd)
                                break
                        if is_fwd_target:
                            break
                if is_fwd_target and source_tidx is not None:
                    # forwarding-target think time: thinkt = njobs/arrival_rate - host_residt, derived from source throughput and forwarding probability.
                    njobs_fwd = max(self.njobs[tidx, :])
                    arrival_rate = self.tput[source_tidx] * fwd_prob
                    if arrival_rate > GlobalConstants.FineTol and njobs_fwd > 0:
                        self.tput[tidx] = arrival_rate
                        # Subtract the processor response time for the target's
                        # activities (already computed by _update_metrics_default)
                        target_eidx = int(lqn.callpair[fwd_cidx_found, 1])
                        host_residt = 0.0
                        if hasattr(lqn, 'actsof') and target_eidx in lqn.actsof:
                            for aidx_fwd in lqn.actsof[target_eidx]:
                                host_residt += self.residt[aidx_fwd]
                        self.thinkt[tidx] = max(GlobalConstants.Zero,
                                                njobs_fwd / arrival_rate - host_residt - tidx_thinktime)
                    else:
                        # Source throughput not yet available; use large think time
                        self.thinkt[tidx] = 1000
                    # Apply under-relaxation
                    omega = self.relax_omega
                    if omega < 1.0 and it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                        self.thinkt[tidx] = omega * self.thinkt[tidx] + (1 - omega) * self.thinkt_prev[tidx]
                    self.thinkt_prev[tidx] = self.thinkt[tidx]
                    self.thinktproc[tidx] = Exp.fit_mean(self.thinkt[tidx] + tidx_thinktime)
                    continue

                # Ref task - think time is just user-specified think time
                self.thinkt[tidx] = GlobalConstants.FineTol
                self.thinktproc[tidx] = Exp.fit_mean(tidx_thinktime) if tidx_thinktime > 0 else Immediate()

                # Get REF task throughput from its HOST layer
                # REF tasks are clients in their host layer, so get TN from there
                hidx = self._get_parent(tidx)  # host = parent of task
                if hidx is not None and not np.isnan(self.idxhash[hidx]):
                    host_layer_idx = int(self.idxhash[hidx])
                    if host_layer_idx >= 0 and len(self.results) > 0 and host_layer_idx < len(self.results[-1]):
                        result = self.results[-1][host_layer_idx]
                        if result is not None and 'TN' in result:
                            TN = result['TN']
                            UN = result.get('UN', np.zeros_like(TN))
                            # Find the task class in the host layer using layer.attribute['tasks']
                            # REF tasks are NOT in thinkt_classes_updmap, but ARE in layer.attribute['tasks']
                            layer = self.ensemble[host_layer_idx]
                            if layer is not None and hasattr(layer, 'attribute'):
                                tasks_attr = layer.attribute.get('tasks', [])
                                client_idx = layer.attribute.get('clientIdx', 1)
                                nodeidx_0 = client_idx - 1 if client_idx >= 1 else 0  # Client node for TN extraction
                                for task_entry in tasks_attr:
                                    class_idx_1based = task_entry[0]  # 1-indexed class index
                                    task_tidx = task_entry[1]  # Task's absolute index
                                    if task_tidx == tidx:
                                        classidx_0 = class_idx_1based - 1  # Convert to 0-indexed
                                        if nodeidx_0 < TN.shape[0] and classidx_0 < TN.shape[1]:
                                            self.tput[tidx] = TN[nodeidx_0, classidx_0]
                                            self.util[tidx] = UN[nodeidx_0, classidx_0]
                                        break

    def _init_interlock(self):
        """Build the interlock path table and locate the common parents.

        Interlocking arises when requests issued by one client reach a common
        lower-level server along two or more independent paths, so that
        arrivals a layer decomposition treats as independent are in fact
        correlated. Franks (1999), Ch. 4:
          Phase A: the path table path(a,b) of Sec. 4.2, the calls to entry b
                   caused by one invocation of entry a, with a unit diagonal;
                   a second table restricts the count to the phase-1 flow
          Phase B: the common-parent finder of Fig. 4.2, retaining only the
                   entries at which the flow genuinely splits
          Phase C: the source tasks and the source count n_s of Eq. (4.7)
        The phase-aware tables, the branch-point test and the source count are
        refinements beyond the published algorithm, which assumes one path
        table and counts source tasks directly.

        The interlock table is built once at solver initialization and reused
        across iterations. Only the interlock flow computation (in update_populations)
        uses iteration-dependent throughput values.
        """
        lqn = self.lqn

        # Phase A: Build interlock reachability table
        nentries = lqn.nentries
        il_all = np.zeros((nentries, nentries))
        il_ph1 = np.zeros((nentries, nentries))

        for e in range(nentries):
            eidx = lqn.eshift + e
            visited = np.zeros(nentries, dtype=bool)
            self._trace_interlock_paths(eidx, e, 1.0, 1.0, visited, il_all, il_ph1, 0)

        self.il_table_all = il_all
        self.il_table_ph1 = il_ph1

        # Phase B+C: Find common entries and sources per server entity
        max_idx = lqn.tshift + lqn.ntasks
        self.il_common_entries = [None] * max_idx
        self.il_source_tasks_all = [None] * max_idx
        self.il_source_tasks_ph2 = [None] * max_idx
        self.il_num_sources = np.zeros(max_idx)

        # Process task servers
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self._is_ref_task(tidx) or self._get_sched(tidx) == SchedStrategy.INF:
                continue
            ce, sa, sp, ns = self._find_interlock_for_server(tidx, il_all, il_ph1)
            self.il_common_entries[tidx] = ce
            self.il_source_tasks_all[tidx] = sa
            self.il_source_tasks_ph2[tidx] = sp
            self.il_num_sources[tidx] = ns

        # Process host servers
        for h in range(lqn.nhosts):
            hidx = h
            if self._get_sched(hidx) == SchedStrategy.INF:
                continue
            ce, sa, sp, ns = self._find_interlock_for_server(hidx, il_all, il_ph1)
            self.il_common_entries[hidx] = ce
            self.il_source_tasks_all[hidx] = sa
            self.il_source_tasks_ph2[hidx] = sp
            self.il_num_sources[hidx] = ns

    def _trace_interlock_paths(self, eidx: int, root_e: int, prob_all: float,
                               prob_ph1: float, visited: np.ndarray,
                               il_all: np.ndarray, il_ph1: np.ndarray, depth: int):
        """Phase A: Recursive path tracing for interlock reachability."""
        lqn = self.lqn
        e = eidx - lqn.eshift
        if e < 1 or e > lqn.nentries:
            return
        if visited[e]:
            return
        visited[e] = True

        # Record reachability from root to this entry
        il_all[root_e, e] += prob_all
        il_ph1[root_e, e] += prob_ph1

        # Follow synchronous calls from activities of this entry
        acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
        for aidx in acts:
            if aidx < lqn.ashift or aidx >= lqn.ashift + lqn.nacts:
                continue
            a = aidx - lqn.ashift

            # Pruning: at non-root entries (depth > 0), skip phase-2+ activities
            has_actphase = hasattr(lqn, 'actphase') and lqn.actphase is not None
            if depth > 0 and has_actphase and a - 1 < len(lqn.actphase) and lqn.actphase[a - 1] > 1:
                continue

            is_ph1 = True
            if has_actphase and a - 1 < len(lqn.actphase) and lqn.actphase[a - 1] > 1:
                is_ph1 = False

            # Follow calls from this activity
            calls_from_act = lqn.callsof.get(aidx, []) if isinstance(lqn.callsof, dict) else []
            for cidx in calls_from_act:
                if cidx < 0 or cidx >= lqn.ncalls:
                    continue
                # Check SYNC call
                if isinstance(lqn.calltype, np.ndarray):
                    ct = int(lqn.calltype.flatten()[cidx]) if cidx < len(lqn.calltype.flatten()) else 0
                elif isinstance(lqn.calltype, dict):
                    ct = lqn.calltype.get(cidx, 0)
                else:
                    ct = 0
                if ct != CallType.SYNC:
                    continue

                call_mean = self._get_call_mean(cidx)
                if call_mean <= 0:
                    continue

                dst_eidx = int(lqn.callpair[cidx, 1])
                dst_e = dst_eidx - lqn.eshift
                if dst_e < 1 or dst_e > lqn.nentries:
                    continue

                next_all = prob_all * call_mean
                next_ph1 = prob_ph1 * call_mean if is_ph1 else 0.0

                self._trace_interlock_paths(dst_eidx, root_e, next_all, next_ph1,
                                            visited, il_all, il_ph1, depth + 1)

        visited[e] = False

    def _get_server_entry_nums(self, server_idx: int) -> List[int]:
        """Get entry numbers (1-based, relative to eshift) for a server."""
        lqn = self.lqn
        nums = []
        if server_idx <= lqn.nhosts:
            # Host server: entries of all tasks on this host
            tasks = lqn.tasksof.get(server_idx, []) if isinstance(lqn.tasksof, dict) else []
            for tidx in tasks:
                entries = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else []
                for se in entries:
                    nums.append(se - lqn.eshift)
        else:
            # Task server: entries of this task
            entries = lqn.entriesof.get(server_idx, []) if isinstance(lqn.entriesof, dict) else []
            for se in entries:
                nums.append(se - lqn.eshift)
        return nums

    def _get_client_tasks(self, server_idx: int) -> List[int]:
        """Get client task indices for a server."""
        lqn = self.lqn
        if server_idx <= lqn.nhosts:
            return lqn.tasksof.get(server_idx, []) if isinstance(lqn.tasksof, dict) else []
        else:
            server_entries = lqn.entriesof.get(server_idx, []) if isinstance(lqn.entriesof, dict) else []
            client_tasks = []
            for se in server_entries:
                if hasattr(lqn, 'iscaller') and lqn.iscaller is not None and isinstance(lqn.iscaller, np.ndarray):
                    if se < lqn.iscaller.shape[1]:
                        calling_idx = np.where(lqn.iscaller[:, se] > 0)[0]
                        for ci in calling_idx:
                            if lqn.tshift <= ci < lqn.tshift + lqn.ntasks:
                                if ci not in client_tasks:
                                    client_tasks.append(ci)
            return client_tasks

    def _has_phase2_activities(self, eidx: int) -> bool:
        """Check if entry has phase-2 activities."""
        lqn = self.lqn
        if not hasattr(lqn, 'actphase') or lqn.actphase is None:
            return False
        acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
        for aidx in acts:
            a = aidx - lqn.ashift
            if 1 <= a <= lqn.nacts and a - 1 < len(lqn.actphase) and lqn.actphase[a - 1] > 1:
                return True
        return False

    def _is_branch_point_check(self, src_x_eidx: int, entry_a_eidx: int,
                               src_y_eidx: int, entry_b_eidx: int,
                               il_all: np.ndarray) -> bool:
        """Check if (srcX, srcY) form a branch point for (entryA, entryB)."""
        lqn = self.lqn
        task_a = self._get_parent(entry_a_eidx)
        task_b = self._get_parent(entry_b_eidx)
        task_x = self._get_parent(src_x_eidx)

        # Multiserver client: if X, A, B same task => not branch point
        if task_x == task_a and task_x == task_b:
            return False

        # Quick check: direct call
        if src_x_eidx == entry_a_eidx or src_y_eidx == entry_b_eidx:
            return True

        entry_a_num = entry_a_eidx - lqn.eshift
        entry_b_num = entry_b_eidx - lqn.eshift

        # Check downstream calls diverge to different tasks
        dst_tasks_x = self._get_call_dst_tasks(src_x_eidx, entry_a_num, il_all)
        dst_tasks_y = self._get_call_dst_tasks(src_y_eidx, entry_b_num, il_all)

        for dx in dst_tasks_x:
            for dy in dst_tasks_y:
                if dx != dy:
                    return True
        return False

    def _get_call_dst_tasks(self, src_eidx: int, target_e_num: int,
                            il_all: np.ndarray) -> List[int]:
        """Get destination tasks of sync calls from an entry reaching a target."""
        lqn = self.lqn
        dst_tasks = []
        acts = lqn.actsof.get(src_eidx, []) if isinstance(lqn.actsof, dict) else []
        for aidx in acts:
            if aidx < lqn.ashift or aidx >= lqn.ashift + lqn.nacts:
                continue
            calls = lqn.callsof.get(aidx, []) if isinstance(lqn.callsof, dict) else []
            for cidx in calls:
                if cidx < 0 or cidx >= lqn.ncalls:
                    continue
                if isinstance(lqn.calltype, np.ndarray):
                    ct = int(lqn.calltype.flatten()[cidx]) if cidx < len(lqn.calltype.flatten()) else 0
                elif isinstance(lqn.calltype, dict):
                    ct = lqn.calltype.get(cidx, 0)
                else:
                    ct = 0
                if ct != CallType.SYNC:
                    continue
                dst_eidx = int(lqn.callpair[cidx, 1])
                dst_e = dst_eidx - lqn.eshift
                if 1 <= dst_e <= lqn.nentries and il_all[dst_e, target_e_num] > 0:
                    parent = self._get_parent(dst_eidx)
                    if parent is not None and parent not in dst_tasks:
                        dst_tasks.append(parent)
        return dst_tasks

    def _find_interlocked_tasks(self, src_eidx: int, server_idx: int,
                                il_all: np.ndarray) -> List[int]:
        """Get interlocked tasks on paths from an entry to a server."""
        lqn = self.lqn
        visited = np.zeros(lqn.nentries, dtype=bool)
        return self._trace_to_server_rec(src_eidx, server_idx, il_all, visited, [], True)

    def _trace_to_server_rec(self, eidx: int, server_idx: int,
                             il_all: np.ndarray, visited: np.ndarray,
                             itasks: List[int], is_head: bool) -> List[int]:
        """Recursively trace paths from entry to server, collecting interlocked tasks."""
        lqn = self.lqn
        e = eidx - lqn.eshift
        if e < 1 or e > lqn.nentries or visited[e]:
            return itasks

        owner_task = self._get_parent(eidx)
        if owner_task is None:
            return itasks

        # Check if we reached the server
        if owner_task == server_idx:
            return itasks
        if server_idx <= lqn.nhosts:
            parent_of_owner = self._get_parent(owner_task)
            if parent_of_owner == server_idx:
                return itasks

        visited[e] = True

        # Follow synchronous calls from ALL phases
        acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
        found = False
        for aidx in acts:
            if aidx < lqn.ashift or aidx >= lqn.ashift + lqn.nacts:
                continue
            calls = lqn.callsof.get(aidx, []) if isinstance(lqn.callsof, dict) else []
            for cidx in calls:
                if cidx < 0 or cidx >= lqn.ncalls:
                    continue
                if isinstance(lqn.calltype, np.ndarray):
                    ct = int(lqn.calltype.flatten()[cidx]) if cidx < len(lqn.calltype.flatten()) else 0
                elif isinstance(lqn.calltype, dict):
                    ct = lqn.calltype.get(cidx, 0)
                else:
                    ct = 0
                if ct != CallType.SYNC:
                    continue

                dst_eidx = int(lqn.callpair[cidx, 1])
                dst_task = self._get_parent(dst_eidx)

                # Check if destination reaches server
                reaches_server = False
                if dst_task == server_idx:
                    reaches_server = True
                elif server_idx <= lqn.nhosts and self._get_parent(dst_task) == server_idx:
                    reaches_server = True
                else:
                    dst_e = dst_eidx - lqn.eshift
                    server_entry_nums = self._get_server_entry_nums(server_idx)
                    for se_num in server_entry_nums:
                        if 1 <= dst_e <= lqn.nentries and 1 <= se_num <= lqn.nentries and il_all[dst_e, se_num] > 0:
                            reaches_server = True
                            break

                if reaches_server:
                    itasks = self._trace_to_server_rec(dst_eidx, server_idx, il_all, visited, itasks, False)
                    found = True

        if found and not is_head:
            if owner_task not in itasks:
                itasks.append(owner_task)

        visited[e] = False
        return itasks

    def _find_interlock_for_server(self, server_idx: int,
                                   il_all: np.ndarray, il_ph1: np.ndarray):
        """Phase B+C: Find interlock for a single server.

        Returns (commonEntries, srcAll, srcPh2, numSources).
        """
        lqn = self.lqn
        empty = ([], [], [], 0)

        # Get server entry numbers
        server_entry_nums = self._get_server_entry_nums(server_idx)
        if not server_entry_nums:
            return empty

        # Get client tasks
        client_tasks = self._get_client_tasks(server_idx)
        if len(client_tasks) < 1:
            return empty

        # Get client entries that reach the server
        client_entry_pairs = []  # list of (taskIdx, entryNum)
        for ct in client_tasks:
            entries = lqn.entriesof.get(ct, []) if isinstance(lqn.entriesof, dict) else []
            for ce in entries:
                ce_num = ce - lqn.eshift
                if ce_num < 1 or ce_num > lqn.nentries:
                    continue
                for se_num in server_entry_nums:
                    if 1 <= se_num <= lqn.nentries and il_all[ce_num, se_num] > 0:
                        client_entry_pairs.append((ct, ce_num))
                        break

        if len(client_entry_pairs) < 2:
            return empty

        # Find common parent entries (branch points)
        common_entries_set = []
        n_pairs = len(client_entry_pairs)
        for i in range(n_pairs):
            for j in range(i + 1, n_pairs):
                if client_entry_pairs[i][0] == client_entry_pairs[j][0]:
                    continue  # Same task
                entry_a_num = client_entry_pairs[i][1]
                entry_c_num = client_entry_pairs[j][1]

                # Search all tasks for common parents
                for t in range(lqn.ntasks):
                    tidx = lqn.tshift + t
                    entries_of_task = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else []
                    for ex in entries_of_task:
                        for ey in entries_of_task:
                            ex_num = ex - lqn.eshift
                            ey_num = ey - lqn.eshift
                            if ex_num < 1 or ey_num < 1 or ex_num > lqn.nentries or ey_num > lqn.nentries:
                                continue
                            if il_all[ex_num, entry_a_num] > 0 and il_all[ey_num, entry_c_num] > 0:
                                if self._is_branch_point_check(
                                        ex, lqn.eshift + entry_a_num,
                                        ey, lqn.eshift + entry_c_num, il_all):
                                    if ex not in common_entries_set:
                                        common_entries_set.append(ex)

        # Unique
        common_entries_set = sorted(set(common_entries_set))

        if not common_entries_set:
            return empty

        # Phase C: Find source tasks
        interlocked_tasks = []
        for ce_eidx in common_entries_set:
            it = self._find_interlocked_tasks(ce_eidx, server_idx, il_all)
            for t in it:
                if t not in interlocked_tasks:
                    interlocked_tasks.append(t)

        # All source tasks = tasks owning common entries
        all_src_tasks = []
        for ce_eidx in common_entries_set:
            owner_tidx = self._get_parent(ce_eidx)
            if owner_tidx is not None and owner_tidx not in all_src_tasks:
                all_src_tasks.append(owner_tidx)

        # Remove interlocked tasks from allSrcTasks
        all_src_tasks = [t for t in all_src_tasks if t not in interlocked_tasks]

        # Ph2 sources: interlocked tasks with phase-2 activities reaching server
        ph2_src_tasks = []
        for it in interlocked_tasks:
            entries_it = lqn.entriesof.get(it, []) if isinstance(lqn.entriesof, dict) else []
            for ie in entries_it:
                if self._has_phase2_activities(ie):
                    ie_num = ie - lqn.eshift
                    if 1 <= ie_num <= lqn.nentries:
                        for se_num in server_entry_nums:
                            if 1 <= se_num <= lqn.nentries:
                                if il_all[ie_num, se_num] - il_ph1[ie_num, se_num] > 0:
                                    if it not in ph2_src_tasks:
                                        ph2_src_tasks.append(it)
                                    break

        # Add external sources (tasks calling into interlocked paths from outside)
        for it in interlocked_tasks:
            entries_it = lqn.entriesof.get(it, []) if isinstance(lqn.entriesof, dict) else []
            for ie in entries_it:
                if hasattr(lqn, 'iscaller') and lqn.iscaller is not None and isinstance(lqn.iscaller, np.ndarray):
                    if ie < lqn.iscaller.shape[1]:
                        calling_idx = np.where(lqn.iscaller[:, ie] > 0)[0]
                        for ci in calling_idx:
                            if lqn.tshift <= ci < lqn.tshift + lqn.ntasks:
                                if ci not in interlocked_tasks and ci not in all_src_tasks:
                                    all_src_tasks.append(ci)

        # Count total source multiplicity
        nsrc = 0
        for st in all_src_tasks:
            nsrc += self._get_mult(st)

        return common_entries_set, all_src_tasks, ph2_src_tasks, nsrc

    def _get_hostdem_mean(self, aidx: int) -> float:
        """Get mean host demand for an activity index."""
        lqn = self.lqn
        if isinstance(lqn.hostdem, dict):
            val = lqn.hostdem.get(aidx, 0.0)
        elif isinstance(lqn.hostdem, (list, np.ndarray)):
            if aidx < len(lqn.hostdem):
                val = lqn.hostdem[aidx]
            else:
                val = 0.0
        else:
            val = 0.0
        if val is None:
            return 0.0
        if isinstance(val, (int, float, np.integer, np.floating)):
            return float(val)
        if hasattr(val, 'getMean'):
            return val.getMean()
        if hasattr(val, 'mean'):
            return val.mean
        if hasattr(val, 'get_mean'):
            return val.get_mean()
        return 0.0

    def _get_calltype(self, cidx: int) -> int:
        """Get call type for a call index."""
        lqn = self.lqn
        if isinstance(lqn.calltype, np.ndarray):
            flat = lqn.calltype.flatten()
            return int(flat[cidx]) if cidx < len(flat) else 0
        elif isinstance(lqn.calltype, dict):
            return lqn.calltype.get(cidx, 0)
        return 0

    def _get_entry_tput(self, eidx: int, task_idx: int) -> float:
        """Get entry throughput (helper for interlock computation)."""
        lqn = self.lqn
        tput_val = self.tput[eidx] if eidx < len(self.tput) else 0.0
        if tput_val <= GlobalConstants.FineTol:
            # Try first activity
            acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
            if acts:
                first_act = acts[0]
                if first_act < len(self.tput):
                    tput_val = self.tput[first_act]
        if tput_val <= GlobalConstants.FineTol:
            if task_idx < len(self.tput):
                tput_val = self.tput[task_idx]
        return float(tput_val)

    def _get_task_tput(self, tidx: int) -> float:
        """Get task throughput (helper for interlock computation)."""
        lqn = self.lqn
        tput_val = self.tput[tidx] if tidx < len(self.tput) else 0.0
        if tput_val <= GlobalConstants.FineTol:
            entries = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else []
            for eidx in entries:
                et = self.tput[eidx] if eidx < len(self.tput) else 0.0
                if et <= GlobalConstants.FineTol:
                    acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
                    if acts:
                        first_act = acts[0]
                        if first_act < len(self.tput):
                            et = self.tput[first_act]
                tput_val += et
        return float(tput_val)

    def _compute_interlock_prob(self, client_tidx: int, server_idx: int,
                                is_processor_host: bool = False):
        """Interlock probability for one (client, server) pair, as (IR, Pr(IL)).

        Li and Franks, "An improved interlocking correction for decomposition of layered
        queueing networks", CCECE 2015, Eqs. (3) and (4). ``is_processor_host`` selects the
        m' rule of lqns ``Interlock::ilrate_pril_flow``: at a PROCESSOR the common-source
        population is doubled above 3 customers and squared at or below it, which is what
        turns m = 4 into the pril = 1/8 its trace reports. The two factors are multiplied
        into the Eq. (5) rate by ``_build_layer_interlock``, so neither carries the source
        count on its own -- that lives in m'. This replaces the superseded (n_s-1)/n_s
        discount of Franks (1999), Eq. (4.7).
        """
        lqn = self.lqn

        if server_idx >= len(self.il_common_entries) or self.il_common_entries[server_idx] is None:
            return 0.0, 0.0
        common_entries = self.il_common_entries[server_idx]
        num_sources = self.il_num_sources[server_idx]
        all_src_tasks = self.il_source_tasks_all[server_idx]
        ph2_src_tasks = self.il_source_tasks_ph2[server_idx]

        if num_sources == 0 or not common_entries:
            return 0.0, 0.0

        # Get client entries
        client_entries = lqn.entriesof.get(client_tidx, []) if isinstance(lqn.entriesof, dict) else []

        # Interlocked flow lambda^IL of Eq. (4), and alongside it the flow weighted by 1/m',
        # which gives Pr(IL) of Eq. (3).
        sum_flow = 0.0
        sum_pril = 0.0
        for ce_eidx in common_entries:
            src_task = self._get_parent(ce_eidx)
            ce_num = ce_eidx - lqn.eshift
            # population of this common source, in customer copies
            m_src = self._get_mult(src_task)
            if not np.isfinite(m_src) or m_src < 1:
                m_src = 1.0
            if is_processor_host:
                m_eff = (m_src + m_src) if m_src > 3 else (m_src * m_src)
            else:
                m_eff = m_src

            for dst_a_eidx in client_entries:
                dst_a_num = dst_a_eidx - lqn.eshift
                if dst_a_num < 1 or dst_a_num > lqn.nentries:
                    continue
                if self.il_table_all[ce_num, dst_a_num] <= 0:
                    continue

                # Get source entry throughput
                ce_tput = self._get_entry_tput(ce_eidx, src_task)
                if ce_tput <= GlobalConstants.FineTol:
                    continue

                # Deferred flow is scored separately from the phase-1 flow
                has_p2 = self._has_phase2_activities(ce_eidx)

                if not has_p2 and src_task in all_src_tasks:
                    contrib = ce_tput * self.il_table_all[ce_num, dst_a_num]
                    sum_flow += contrib
                    sum_pril += contrib / m_eff
                elif has_p2 and src_task in all_src_tasks:
                    contrib = ce_tput * self.il_table_ph1[ce_num, dst_a_num]
                    sum_flow += contrib
                    sum_pril += contrib / m_eff

                ph2 = self.il_table_all[ce_num, dst_a_num] - self.il_table_ph1[ce_num, dst_a_num]
                if ph2 > 0 and src_task in ph2_src_tasks:
                    contrib = ce_tput * ph2
                    sum_flow += contrib
                    sum_pril += contrib / m_eff

        # Get client throughput
        client_tput = self._get_task_tput(client_tidx)
        if client_tput <= GlobalConstants.FineTol:
            return 0.0, 0.0

        ir = min(sum_flow, client_tput) / client_tput
        ir = min(1.0, max(0.0, ir))
        pr_il = 0.0 if sum_flow <= GlobalConstants.FineTol else sum_pril / sum_flow
        pr_il = min(1.0, max(0.0, pr_il))
        return float(ir), float(pr_il)

    def update_populations(self, it: int):
        """Apply the interlock correction to call residence times.

        The path tables built by _init_interlock are combined with the current
        iterate to obtain, for each (client, server) pair, the interlocked flow
        of Eq. (4.3) of Franks (1999), and from it the interlock probability,
        the share of that flow that the layer decomposition would otherwise
        count twice. Eq. (4.7) removes one source in n_s from the queue length
        inside MVA; the equivalent correction is applied here to the residence
        times returned by the layer::

          R_adj = S + (1 - prIL) * W,   W = R - S,

        which leaves service and utilization untouched and removes only the
        interlocked share of the waiting time.

        Called after update_metrics, which produces the raw callresidt from the
        layer solutions, and before update_think_times.
        """
        lqn = self.lqn

        if self.il_table_all is None:
            return

        # Save originals for proportional entry_servt update
        callresidt_orig = self.callresidt.copy() if self.callresidt is not None else np.zeros(1)
        residt_orig = self.residt.copy() if self.residt is not None else np.zeros(1)
        adjusted = False

        # Pass 1: For each sync call, check if destination server has interlock
        for cidx in range(lqn.ncalls):
            if self._get_calltype(cidx) != CallType.SYNC:
                continue

            dst_eidx = int(lqn.callpair[cidx, 1])
            server_tidx = self._get_parent(dst_eidx)
            if server_tidx is None:
                continue

            # Find the server entity with interlock data
            server_for_il = None
            if (server_tidx < len(self.il_common_entries)
                    and self.il_common_entries[server_tidx] is not None
                    and len(self.il_common_entries[server_tidx]) > 0):
                server_for_il = server_tidx
            else:
                # Check host server
                if server_tidx >= lqn.tshift:
                    host_idx = self._get_parent(server_tidx)
                    if (host_idx is not None and 0 <= host_idx < len(self.il_common_entries)
                            and self.il_common_entries[host_idx] is not None
                            and len(self.il_common_entries[host_idx]) > 0):
                        server_for_il = host_idx

            if server_for_il is None:
                continue

            # Get client task (activity -> task via parent)
            src_aidx = int(lqn.callpair[cidx, 0])
            client_tidx = self._get_parent(src_aidx)
            if client_tidx is None:
                continue

            # Interlock probability for this client and server. This path serves a TASK, not a
            # processor, so the m' rule of Li/lqns leaves the source population alone; the
            # product IR*Pr(IL) reproduces the scalar this branch used before.
            ir_c, pr_il_c = self._compute_interlock_prob(client_tidx, server_for_il, False)
            pr_il = ir_c * pr_il_c

            if pr_il <= GlobalConstants.FineTol:
                continue

            # Compute waiting time reduction
            S = self.servt[dst_eidx] if dst_eidx < len(self.servt) else 0.0
            call_mean = self._get_call_mean(cidx)

            if call_mean <= 0 or self.callservt[cidx] <= 0:
                continue

            RN = self.callservt[cidx] / call_mean  # response time per visit
            W = max(0.0, RN - S)  # waiting time per visit

            if W > GlobalConstants.FineTol:
                RN_adj = S + (1 - pr_il) * W
                scale = RN_adj / RN
                self.callservt[cidx] = self.callservt[cidx] * scale
                self.callresidt[cidx] = self.callresidt[cidx] * scale
                if self.callservt[cidx] > 0:
                    self.callservtproc[cidx] = Exp.fit_mean(self.callservt[cidx])
                adjusted = True

        # Pass 2: Host-level interlock — reduce processor queueing in residt.
        # Every layer starts the pass without a matrix, so a host that stops being
        # interlocked does not keep the previous iteration's correction alive.
        for e in range(len(self.ensemble)):
            slv = self.solvers[e] if e < len(self.solvers) else None
            cfg = getattr(getattr(slv, 'options', None), 'config', None)
            if isinstance(cfg, dict):
                cfg['interlock'] = None
            elif cfg is not None and hasattr(cfg, 'interlock'):
                cfg.interlock = None
        for h in range(lqn.nhosts):
            hidx = h
            if self.il_common_entries[hidx] is None or len(self.il_common_entries[hidx]) == 0:
                continue

            # Compute prIL and processor utilization for each task on this host
            host_tasks = lqn.tasksof.get(hidx, []) if isinstance(lqn.tasksof, dict) else []
            task_pr_il = np.zeros(len(host_tasks))   # IR, Eq. (4)
            task_PrIL = np.zeros(len(host_tasks))    # Pr(IL), Eq. (3)
            task_util = np.zeros(len(host_tasks))
            for ti, tidx in enumerate(host_tasks):
                # The host of a task layer is a PROCESSOR, which is what selects the m' rule.
                task_pr_il[ti], task_PrIL[ti] = self._compute_interlock_prob(tidx, hidx, True)
                # Compute task's processor utilization
                entries = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else []
                for eidx in entries:
                    acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
                    for aidx in acts:
                        if aidx < len(self.tput):
                            task_util[ti] += self.tput[aidx] * self._get_hostdem_mean(aidx)

            U_total = np.sum(task_util)
            U_interlocked = np.sum(task_util[task_pr_il > GlobalConstants.FineTol])
            if U_total <= GlobalConstants.FineTol or U_interlocked <= GlobalConstants.FineTol:
                continue
            il_fraction = U_interlocked / U_total

            # When the layer solver carries Eq. (4.7) inside its own MVA, the interlock goes
            # to the layer as a class-level matrix and the residence times are left untouched.
            # Scaling them here as well would remove the same waiting twice, and would still
            # leave the layer's own THROUGHPUT uncorrected, which is what breaks flow balance
            # across a call: the reported task rate then comes from a cycle time the correction
            # has already shortened elsewhere.
            layer_of_host = self._layer_index_of(hidx)
            if layer_of_host is not None and self._layer_takes_interlock(layer_of_host):
                il_mat = self._build_layer_interlock(layer_of_host, host_tasks, task_pr_il, task_PrIL)
                opts = self.solvers[layer_of_host].options
                cfg = getattr(opts, 'config', None)
                if cfg is None:
                    opts.config = {'interlock': il_mat}
                elif isinstance(cfg, dict):
                    cfg['interlock'] = il_mat
                else:
                    cfg.interlock = il_mat
                continue

            for ti, tidx in enumerate(host_tasks):
                if task_pr_il[ti] <= GlobalConstants.FineTol:
                    continue
                # Weight by the share of host utilization that is interlocked. The rate is the
                # SAME Eq. (5) product IR*Pr(IL) that pass 1 applies to a call and that
                # _build_layer_interlock puts in the layer matrix -- IR alone is a flow SHARE,
                # ~1 whenever a layer has a single common source, and using it here removed the
                # whole processor queueing rather than the interlocked part of it, which broke
                # flow balance across a call.
                effective_pr_il = task_pr_il[ti] * task_PrIL[ti] * il_fraction
                entries = lqn.entriesof.get(tidx, []) if isinstance(lqn.entriesof, dict) else []
                for eidx in entries:
                    acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
                    for aidx in acts:
                        D = self._get_hostdem_mean(aidx)
                        if D > 0 and aidx < len(self.residt) and self.residt[aidx] > D + GlobalConstants.FineTol:
                            W_proc = self.residt[aidx] - D
                            self.residt[aidx] = D + (1 - effective_pr_il) * W_proc
                            adjusted = True

        if not adjusted:
            return

        # Recompute entry service times from adjusted callresidt/residt
        # Use proportional scaling to preserve visit ratio adjustments
        residt_vec = self.residt.flatten() if self.residt is not None else np.zeros(1)
        callresidt_vec = self.callresidt.flatten() if self.callresidt is not None else np.zeros(1)
        residt_orig_vec = residt_orig.flatten()
        callresidt_orig_vec = callresidt_orig.flatten()

        # The servtmatrix column space is [element 0..nidx-1, call nidx..nidx+ncalls-1]
        # and callresidt is ncalls long with no leading pad, so the two concatenate
        # directly. Dropping a leading entry here would shift every call one column
        # to the left, crediting each entry with its NEXT call's residence.
        concat_old = np.concatenate([residt_orig_vec, callresidt_orig_vec])
        concat_new = np.concatenate([residt_vec, callresidt_vec])

        if self.servtmatrix is not None:
            # Ensure dimensions match
            n_sm_cols = self.servtmatrix.shape[1]
            if len(concat_old) < n_sm_cols:
                concat_old = np.pad(concat_old, (0, n_sm_cols - len(concat_old)))
                concat_new = np.pad(concat_new, (0, n_sm_cols - len(concat_new)))
            elif len(concat_old) > n_sm_cols:
                concat_old = concat_old[:n_sm_cols]
                concat_new = concat_new[:n_sm_cols]

            entry_servt_old = self.servtmatrix @ concat_old
            entry_servt_new = self.servtmatrix @ concat_new

            # The entry servt is rescaled only when it was itself assembled from
            # these residence times, which is the default path. After the moment3
            # pass it is the MEAN OF AN APH CONVOLUTION of the activities' own
            # response laws, and a ratio of residence-time sums is not a
            # correction to it: applying it multiplies the entry law by the
            # entry's visit ratio and reports a service time BELOW that of the
            # single activity the entry contains. The residence times keep their
            # correction either way. See BUGS.md BUG-97.
            moment_laws = self.lnmethod == 'moment3' and self.moment_pass_done

            for eidx in range(lqn.eshift, lqn.eshift + lqn.nentries):
                if eidx < len(entry_servt_old) and entry_servt_old[eidx] > GlobalConstants.FineTol:
                    ratio = entry_servt_new[eidx] / entry_servt_old[eidx]
                    if not moment_laws:
                        if eidx < len(self.servt):
                            self.servt[eidx] = self.servt[eidx] * ratio
                        if eidx < len(self.servt) and self.servt[eidx] > 0:
                            self.servtproc[eidx] = Exp.fit_mean(self.servt[eidx])
                    if eidx < len(self.residt):
                        self.residt[eidx] = self.residt[eidx] * ratio

    def update_layers(self, it: int):
        """Update layer parameters (matches MATLAB updateLayers)."""
        # Under 'srvn.ph' the layer classes are one per caller task and their laws
        # are composed, not read off the update maps -- see _kb/06-solver-catalog.md
        if self._is_ph_encoding():
            self._update_layers_ph(it)
            return
        lqn = self.lqn

        # Update REF task think times in host layers
        # REF tasks' think times = base_think + call_response_time
        for hidx in range(lqn.nhosts):
            if np.isnan(self.idxhash[hidx]):
                continue
            layer_idx = int(self.idxhash[hidx])
            if layer_idx < 0 or layer_idx >= len(self.ensemble):
                continue
            layer = self.ensemble[layer_idx]
            if layer is None:
                continue

            # Find REF tasks on this host
            tasks_on_host = self._get_tasks_of_host(hidx)
            classes = layer.get_classes()
            nodes = layer.get_nodes()

            # Find the Clients delay node (first node, typically index 0)
            clients_node = None
            for node in nodes:
                if isinstance(node, Delay):
                    clients_node = node
                    break

            if clients_node is None:
                continue

            for class_idx, tidx in enumerate(tasks_on_host):
                if self._is_ref_task(tidx):
                    # REF task TASK-class think time = base_think only; call response times are carried by separate CALL classes.
                    base_think = 0.0
                    if self.thinkproc[tidx] is not None:
                        proc = self.thinkproc[tidx]
                        if hasattr(proc, 'getMean'):
                            base_think = proc.getMean()
                        elif hasattr(proc, 'mean'):
                            base_think = proc.mean

                    # TASK class think time stays as base_think (no call response added)
                    if base_think > 0 and class_idx < len(classes):
                        cls = classes[class_idx]
                        clients_node.set_service(cls, Exp.fit_mean(base_think))

        # Update think times in layers (matches MATLAB updateLayers.m lines 14-79)
        if self.thinkt_classes_updmap is not None:
            for r in range(len(self.thinkt_classes_updmap)):
                # Elevator iteration order
                if it % 2 == 1:
                    ri = len(self.thinkt_classes_updmap) - r - 1
                else:
                    ri = r

                idx = int(self.thinkt_classes_updmap[ri, 0])
                aidx = int(self.thinkt_classes_updmap[ri, 1])
                nodeidx = int(self.thinkt_classes_updmap[ri, 2])
                classidx = int(self.thinkt_classes_updmap[ri, 3])

                layer_idx = int(self.idxhash[idx])
                if layer_idx >= 0 and layer_idx < len(self.ensemble):
                    layer = self.ensemble[layer_idx]
                    if layer is not None:
                        classes = layer.get_classes()
                        nodes = layer.get_nodes()

                        if classidx <= len(classes) and nodeidx <= len(nodes):
                            cls = classes[classidx - 1]
                            node = nodes[nodeidx - 1]

                            # Get clientIdx and serverIdx from layer attribute
                            client_idx = layer.attribute.get('clientIdx', 1) if hasattr(layer, 'attribute') else 1
                            server_idx = layer.attribute.get('serverIdx', 2) if hasattr(layer, 'attribute') else 2

                            lqn_type = self._get_type(aidx)

                            if nodeidx == client_idx:
                                # Client node handling (MATLAB lines 37-75)
                                if lqn_type == LayeredNetworkElement.TASK:
                                    is_ref = self._is_ref_task(aidx)
                                    if not is_ref:
                                        # Non-REF TASK: use computed thinktproc (MATLAB line 42)
                                        if self.thinktproc[aidx] is not None:
                                            node.set_service(cls, self.thinktproc[aidx])
                                    else:
                                        # REF TASK: use servtproc (host demand) — matches MATLAB line 43 and JAR line 2665
                                        if self.servtproc[aidx] is not None:
                                            node.set_service(cls, self.servtproc[aidx])
                                else:
                                    # Non-TASK types (ACTIVITY, ENTRY): use servtproc (MATLAB line 66)
                                    if self.servtproc[aidx] is not None:
                                        node.set_service(cls, self.servtproc[aidx])
                            else:
                                # Server replica (any of them) (MATLAB line 77)
                                if self.servtproc[aidx] is not None:
                                    node.set_service(cls, self.servtproc[aidx])

                            # Propagate service updates to all replicas when nreplicas > 1
                            nrep = layer.attribute.get('nreplicas', 1) if hasattr(layer, 'attribute') else 1
                            if nrep > 1:
                                all_ss = layer.attribute.get('server_stations', [])
                                for replica_ss in all_ss[1:]:  # Skip primary
                                    dist = node._service_process.get(cls)
                                    if dist is not None:
                                        replica_ss.set_service(cls, dist)

        # servt_classes_updmap is read-only here (extracts RN in update_metrics); host-layer processor service time is the fixed host demand.

        # Reassign call service times / response times (like MATLAB lines 106-142)
        if self.call_classes_updmap is not None and len(self.call_classes_updmap) > 0:
            for c in range(len(self.call_classes_updmap)):
                # Elevator iteration order
                if it % 2 == 1:
                    ci = len(self.call_classes_updmap) - c - 1
                else:
                    ci = c

                idx = int(self.call_classes_updmap[ci, 0])
                cidx = int(self.call_classes_updmap[ci, 1])
                nodeidx = int(self.call_classes_updmap[ci, 2])
                classidx = int(self.call_classes_updmap[ci, 3])

                # Get the layer
                if np.isnan(self.idxhash[idx]):
                    continue
                layer_idx_actual = int(self.idxhash[idx])
                if layer_idx_actual < 0 or layer_idx_actual >= len(self.ensemble):
                    continue
                layer = self.ensemble[layer_idx_actual]
                if layer is None:
                    continue

                classes = layer.get_classes()
                nodes = layer.get_nodes()

                if classidx > len(classes) or nodeidx > len(nodes):
                    continue

                cls = classes[classidx - 1]
                node = nodes[nodeidx - 1]

                # Get clientIdx and serverIdx from layer attribute
                client_idx = layer.attribute.get('clientIdx', 1) if hasattr(layer, 'attribute') else 1
                server_idx = layer.attribute.get('serverIdx', 2) if hasattr(layer, 'attribute') else 2

                if nodeidx == client_idx:
                    # CALL at client: use callservtproc (call response time)
                    if cidx < len(self.callservtproc) and self.callservtproc[cidx] is not None:
                        proc = self.callservtproc[cidx]
                        node.set_service(cls, proc if hasattr(proc, 'getMean') else Exp.fit_mean(float(proc)))
                else:
                    # CALL at server replica (any of them): use servtproc[eidx] (entry service time)
                    eidx_raw = self.lqn.callpair[cidx, 1] if cidx < len(self.lqn.callpair) else None
                    eidx = int(eidx_raw) if eidx_raw is not None and not np.isnan(eidx_raw) else None
                    if eidx is not None and eidx < len(self.servtproc) and self.servtproc[eidx] is not None:
                        proc = self.servtproc[eidx]
                        dist = proc if hasattr(proc, 'getMean') else Exp.fit_mean(float(proc))
                        # A phase-2 entry replies before phase 2 runs, so the caller is
                        # held for residt, not servt. Charging it servt here while its
                        # own layer charges residt makes the two layers settle at
                        # different rates and breaks flow conservation across the call.
                        # Under flat layering both live in one model, where the
                        # correction would be applied twice.
                        if (self.hasPhase2 and not self._is_flat_layering()
                                and self.servt_ph2 is not None and eidx < len(self.servt_ph2)
                                and self.servt_ph2[eidx] > 1e-8
                                and self.residt is not None and eidx < len(self.residt)
                                and self.residt[eidx] > 0):
                            dist = Exp.fit_mean(float(self.residt[eidx]))
                        node.set_service(cls, dist)
                        # Propagate to replicas
                        nrep = layer.attribute.get('nreplicas', 1) if hasattr(layer, 'attribute') else 1
                        if nrep > 1:
                            all_ss = layer.attribute.get('server_stations', [])
                            for replica_ss in all_ss[1:]:
                                replica_ss.set_service(cls, dist)

        # Source arrival rates per iter from lqn.arrival (entry-level) or tputproc (async-call); mirrors updateLayers.m:66-74/JAR SolverLN.java:2702-2724.
        if self.arvproc_classes_updmap is not None and len(self.arvproc_classes_updmap) > 0:
            for r in range(len(self.arvproc_classes_updmap)):
                if it % 2 == 1:
                    ri = len(self.arvproc_classes_updmap) - r - 1
                else:
                    ri = r

                idx = int(self.arvproc_classes_updmap[ri, 0])
                eidx_or_cidx = int(self.arvproc_classes_updmap[ri, 1])
                nodeidx = int(self.arvproc_classes_updmap[ri, 2])
                classidx = int(self.arvproc_classes_updmap[ri, 3])

                if np.isnan(self.idxhash[idx]):
                    continue
                layer_idx_actual = int(self.idxhash[idx])
                if layer_idx_actual < 0 or layer_idx_actual >= len(self.ensemble):
                    continue
                layer = self.ensemble[layer_idx_actual]
                if layer is None:
                    continue

                classes = layer.get_classes()
                nodes = layer.get_nodes()
                if classidx > len(classes) or nodeidx > len(nodes):
                    continue

                cls = classes[classidx - 1]
                node = nodes[nodeidx - 1]

                if eidx_or_cidx < 0:
                    # entry-level open arrival re-applies the static lqn.arrival rate; server service was fixed at creation from the bound activity's host demand.
                    eidx = -eidx_or_cidx
                    if hasattr(self.lqn, 'arrival') and eidx in self.lqn.arrival \
                            and self.lqn.arrival[eidx] is not None:
                        try:
                            node.set_arrival(cls, self.lqn.arrival[eidx])
                        except Exception:
                            pass
                else:
                    # Async-call open arrival: use tputproc at caller activity.
                    cidx = eidx_or_cidx
                    caller_aidx = int(self.lqn.callpair[cidx, 0]) if cidx < len(self.lqn.callpair) else 0
                    if caller_aidx > 0 and caller_aidx < len(self.tputproc) \
                            and self.tputproc[caller_aidx] is not None:
                        try:
                            node.set_arrival(cls, self.tputproc[caller_aidx])
                        except Exception:
                            pass

    def _compute_layer_service_time(self, caller_tidx: int, layer_idx: int, it: int) -> float:
        """
        Compute total service time for a caller in a layer.

        For host layers: caller's activities' host demands
        For task layers: called entry's service time (from servt array)
        """
        lqn = self.lqn
        total_demand = 0.0

        # Check if this is a host layer
        actual_layer_idx = int(self.idxhash[layer_idx])
        is_host = actual_layer_idx in self.hostLayerIndices

        if is_host:
            # Host layer: server service = caller's activities' host demands
            activities = self._get_activities_of_task(caller_tidx)
            for aidx in activities:
                if self.servtproc[aidx] is not None:
                    proc = self.servtproc[aidx]
                    if isinstance(proc, (int, float, np.integer, np.floating)):
                        total_demand += float(proc)
                    elif hasattr(proc, 'getMean'):
                        total_demand += proc.getMean()
                    elif hasattr(proc, 'mean'):
                        total_demand += proc.mean
        else:
            # Task layer: server service = called entry's service time
            # Use the iteratively updated servt values
            total_demand = self._get_layer_call_response_time(caller_tidx, layer_idx)

        return total_demand

    def _get_layer_call_response_time(self, caller_tidx: int, layer_idx: int) -> float:
        """
        Get call response time for synch calls from caller to entries in layer.

        Uses the current servt (entry service time) which is updated iteratively.
        """
        lqn = self.lqn
        total_call_time = 0.0

        if not hasattr(lqn, 'callpair') or lqn.callpair is None:
            return 0.0

        # Find all synch calls from this caller's activities
        activities = self._get_activities_of_task(caller_tidx)
        for aidx in activities:
            if isinstance(lqn.callsof, dict):
                calls = lqn.callsof.get(aidx, [])
            else:
                calls = []

            for cidx in calls:
                # Check call type - assume SYNC if calltype not available
                is_sync = True
                if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                    if isinstance(lqn.calltype, np.ndarray):
                        calltype = lqn.calltype.flatten()[cidx] if cidx < len(lqn.calltype.flatten()) else CallType.SYNC
                    elif isinstance(lqn.calltype, dict):
                        calltype = lqn.calltype.get(cidx, CallType.SYNC)
                    else:
                        calltype = CallType.SYNC
                    is_sync = (calltype == CallType.SYNC)

                if is_sync:
                    # Get target entry (column 2 of callpair)
                    tgt_eidx = self._get_call_target_entry(cidx)
                    if tgt_eidx is None or tgt_eidx == 0:
                        continue

                    # Check if this call targets the server in this layer
                    tgt_tidx = self._get_parent(tgt_eidx)
                    if tgt_tidx != layer_idx:
                        continue

                    # Get call mean (number of calls)
                    call_mean = self._get_call_mean(cidx)

                    # Use the entry's current service time (updated each iteration)
                    entry_resp = self.servt[tgt_eidx] if tgt_eidx < len(self.servt) and self.servt[tgt_eidx] > 0 else 0.0

                    # If entry servt is not yet computed, use host demand estimate
                    if entry_resp <= 0:
                        tgt_activities = self._get_activities_of_entry(tgt_eidx)
                        for tgt_aidx in tgt_activities:
                            if self.servtproc[tgt_aidx] is not None:
                                proc = self.servtproc[tgt_aidx]
                                if isinstance(proc, (int, float, np.integer, np.floating)):
                                    entry_resp += float(proc)
                                elif hasattr(proc, 'getMean'):
                                    entry_resp += proc.getMean()
                                elif hasattr(proc, 'mean'):
                                    entry_resp += proc.mean

                    total_call_time += call_mean * entry_resp

        return total_call_time

    def _get_type(self, idx: int) -> int:
        """Get element type based on index ranges."""
        lqn = self.lqn

        # First try the type array if it exists
        if hasattr(lqn, 'type') and lqn.type is not None:
            if isinstance(lqn.type, dict):
                return lqn.type.get(idx, 0)
            elif isinstance(lqn.type, np.ndarray):
                if idx < len(lqn.type):
                    return int(lqn.type[idx])  # type array is 0-indexed over elements

        # Compute type from index ranges
        if lqn.hshift <= idx < lqn.hshift + lqn.nhosts:
            return LayeredNetworkElement.PROCESSOR
        elif lqn.tshift <= idx < lqn.tshift + lqn.ntasks:
            return LayeredNetworkElement.TASK
        elif lqn.eshift <= idx < lqn.eshift + lqn.nentries:
            return LayeredNetworkElement.ENTRY
        elif lqn.ashift <= idx < lqn.ashift + lqn.nacts:
            return LayeredNetworkElement.ACTIVITY
        return 0

    def update_routing_probabilities(self, it: int):
        """Update routing probabilities (matches MATLAB updateRoutingProbabilities)."""
        if self.route_prob_updmap is None or len(self.route_prob_updmap) == 0:
            return

        if self.unique_route_prob_updmap is None or len(self.unique_route_prob_updmap) == 0:
            return

        for u in range(len(self.unique_route_prob_updmap)):
            # Alternate direction (elevator) like MATLAB
            if it % 2 == 0:
                idx = int(self.unique_route_prob_updmap[u])
            else:
                idx = int(self.unique_route_prob_updmap[len(self.unique_route_prob_updmap) - u - 1])

            layer_idx = int(self.idxhash[idx]) if not np.isnan(self.idxhash[idx]) else -1
            if layer_idx < 0 or layer_idx >= len(self.ensemble):
                continue

            idx_updated = False
            layer = self.ensemble[layer_idx]
            if layer is None:
                continue

            # Get current routing matrix object
            P = layer.get_routing_matrix()
            if P is None:
                P = layer.init_routing_matrix()

            classes = layer.get_classes()
            nodes = layer.get_nodes()

            # Find rows in route_prob_updmap for this idx
            for r in range(len(self.route_prob_updmap)):
                if int(self.route_prob_updmap[r, 0]) != idx:
                    continue

                host = int(self.route_prob_updmap[r, 0])
                tidx_caller = int(self.route_prob_updmap[r, 1])
                eidx = int(self.route_prob_updmap[r, 2])
                nodefrom = int(self.route_prob_updmap[r, 3])
                nodeto = int(self.route_prob_updmap[r, 4])
                classidxfrom = int(self.route_prob_updmap[r, 5])
                classidxto = int(self.route_prob_updmap[r, 6])

                # Get caller's layer results
                caller_layer_idx = int(self.idxhash[tidx_caller]) if not np.isnan(self.idxhash[tidx_caller]) else -1
                if caller_layer_idx < 0 or caller_layer_idx >= len(self.results[-1]):
                    continue

                result = self.results[-1][caller_layer_idx]
                if result is None or 'TN' not in result:
                    continue

                caller_layer = self.ensemble[caller_layer_idx]
                if caller_layer is None:
                    continue

                # cache layers (iscachelayer) use hit/miss throughput instead of entry throughput; mirrors MATLAB's items-based cache-layer test.
                is_cache_layer = layer.attribute.get('iscachelayer', False) if layer.attribute else False

                if is_cache_layer:
                    # MATLAB: for cache nodes, get results from host layer (not caller layer)
                    host_layer_idx = int(self.idxhash[host]) if not np.isnan(self.idxhash[host]) else -1
                    if host_layer_idx < 0 or host_layer_idx >= len(self.results[-1]):
                        continue

                    host_result = self.results[-1][host_layer_idx]
                    if host_result is None or 'TN' not in host_result:
                        continue

                    host_layer = self.ensemble[host_layer_idx]
                    if host_layer is None:
                        continue

                    server_idx = self._station_idx_of(host_layer, host)
                    server_idx_0 = server_idx - 1 if server_idx >= 1 else 0

                    TN = host_result['TN']
                    if TN is None or server_idx_0 >= TN.shape[0]:
                        continue

                    # Get total throughput at server
                    Xtot = np.sum(TN[server_idx_0, :])
                    if Xtot <= 0:
                        continue

                    # Get hit/miss throughput using classidxto
                    # MATLAB: hm_tput = sum(TN(serverIdx, classidxto))
                    cls_to_idx_0 = classidxto - 1 if classidxto >= 1 else 0
                    if cls_to_idx_0 < TN.shape[1]:
                        hm_tput = TN[server_idx_0, cls_to_idx_0]
                    else:
                        hm_tput = 0.0

                    new_prob = hm_tput / Xtot if Xtot > 0 else 0.0
                else:
                    # Non-cache layer: use entry throughput

                    # Get server index from caller layer
                    server_idx = self._station_idx_of(caller_layer, tidx_caller)
                    server_idx_0 = server_idx - 1 if server_idx >= 1 else 0

                    TN = result['TN']
                    if TN is None or server_idx_0 >= TN.shape[0]:
                        continue

                    # Get total throughput at server
                    Xtot = np.sum(TN[server_idx_0, :])
                    if Xtot <= 0:
                        continue

                    # find ALL entry class indices in the caller layer, not just the first; mirrors MATLAB calls(find(calls(:,4)==eidx),1).
                    matching_eidxclasses = []
                    calls_attr = caller_layer.attribute.get('calls', [])
                    for call_info in calls_attr:
                        if len(call_info) >= 4 and call_info[3] == eidx:
                            matching_eidxclasses.append(call_info[0])  # class index

                    # MATLAB: entry_tput = sum(TN(serverIdx, eidxclass))
                    # Sum throughput across ALL matching entry classes
                    if not matching_eidxclasses:
                        entry_tput = 0.0
                    else:
                        entry_tput = 0.0
                        for eidxclass in matching_eidxclasses:
                            eidxclass_0 = eidxclass - 1 if eidxclass >= 1 else 0
                            if eidxclass_0 < TN.shape[1]:
                                entry_tput += TN[server_idx_0, eidxclass_0]

                    new_prob = entry_tput / Xtot if Xtot > 0 else 0.0

                # Get class and node indices (convert from 1-based to 0-based)
                cls_from_idx = classidxfrom - 1 if classidxfrom >= 1 else 0
                cls_to_idx = classidxto - 1 if classidxto >= 1 else 0
                node_from_idx = nodefrom - 1 if nodefrom >= 1 else 0
                node_to_idx = nodeto - 1 if nodeto >= 1 else 0

                if (cls_from_idx < len(classes) and cls_to_idx < len(classes) and
                    node_from_idx < len(nodes) and node_to_idx < len(nodes)):

                    cls_from = classes[cls_from_idx]
                    cls_to = classes[cls_to_idx]
                    node_from = nodes[node_from_idx]
                    node_to = nodes[node_to_idx]

                    # update _original_routes (used for rt computation) with entry_tput/Xtot; mirrors MATLAB P{r,s}(from,to)=entry_tput/Xtot.
                    if P._original_routes is not None:
                        if (cls_from, cls_to) not in P._original_routes:
                            P._original_routes[(cls_from, cls_to)] = {}
                        P._original_routes[(cls_from, cls_to)][(node_from, node_to)] = new_prob
                    else:
                        # No ClassSwitch nodes - update routes directly
                        P.set(cls_from, cls_to, node_from, node_to, new_prob)

                    # Also update ClassSwitch node's switching matrix for rtnodes
                    from line_solver.lang.nodes import ClassSwitch
                    cs_name = f'CS_{node_from.name}_to_{node_to.name}'
                    for n in nodes:
                        if isinstance(n, ClassSwitch) and n.name == cs_name:
                            cs_matrix = n.get_class_switching_matrix()
                            if cs_matrix is not None and cls_from_idx < cs_matrix.shape[0] and cls_to_idx < cs_matrix.shape[1]:
                                cs_matrix[cls_from_idx, cls_to_idx] = new_prob
                                n.set_class_switching_matrix(cs_matrix)
                            break

                    idx_updated = True

            # Clear cached matrix and reset struct to force recomputation
            if idx_updated:
                P._matrix = None  # Clear cached toMatrix() result
                layer.reset_struct()

    def converged(self, it: int) -> bool:
        """Check convergence (matches MATLAB converged)."""
        # Stochastic iteration dispatch: see _kb/06-solver-catalog.md (LN
        # section) for the rationale.
        if self.stochiter_mode is not None:
            if (self.stochiter_auto and self.stochiter_mode == 'off'
                    and it >= 1 and self.stochlayers is not None and np.any(self.stochlayers)):
                # a layer with method 'default' resolved at runtime to a
                # stochastic method (captured in analyze() at iteration 1)
                self.stochiter_mode = 'rm'
                line_debug("LN: stochastic layer method detected at runtime, "
                           "switching to Robbins-Monro iteration")
            if self.stochiter_mode == 'rm':
                return self.converged_stoch(it)

        # The moment3 pass is terminal: it runs once hasconverged is set, and its
        # own output perturbs the error test below. See BUG-97 and the note where
        # moment_pass_done is set.
        if self.lnmethod == 'moment3' and self.moment_pass_done:
            return True

        if it < 2:
            return False

        # MATLAB: iter_min = max([2*length(self.model.ensemble), ceil(self.options.iter_max/4)])
        iter_min = max(2 * self.nlayers, (self.options.iter_max + 3) // 4)  # ceil equivalent

        # Apply moving window average to help convergence (matches MATLAB lines 58-79)
        # MATLAB: wnd_size = max(5, ceil(iter_min/5))
        wnd_size = max(5, (iter_min + 4) // 5)  # ceil equivalent
        mov_avg_weight = 1.0 / wnd_size
        if it >= iter_min and len(self.results) >= it:
            # moving-average smoothing reads the PREVIOUS iteration's raw (unsmoothed) row so repeated smoothing does not compound the window.
            averaged_row = list(self.results[it - 1])
            for e in range(self.nlayers):
                if len(self.results[it - 1]) > e and self.results[it - 1][e] is not None:
                    result = dict(self.results[it - 1][e])
                    # Apply moving average (matches MATLAB exactly)
                    for key in ['QN', 'UN', 'RN', 'TN', 'AN', 'WN']:
                        if key in result and result[key] is not None:
                            # Start with current result * weight
                            avg = mov_avg_weight * result[key].copy()
                            # Add past wnd_size-1 results (MATLAB: for k=1:(wnd_size-1))
                            for k in range(1, wnd_size):
                                hist_idx = it - 1 - k
                                if hist_idx >= 0 and hist_idx < len(self.results) and len(self.results[hist_idx]) > e:
                                    prev_result = self.results[hist_idx][e]
                                    if prev_result is not None and key in prev_result and prev_result[key] is not None:
                                        avg = avg + prev_result[key] * mov_avg_weight
                            result[key] = avg
                    averaged_row[e] = result
            self.results[it - 1] = averaged_row

        # Compute max error across all layers
        if it > 1 and len(self.results) >= it:
            self.maxitererr.append(0.0)

            for e in range(self.nlayers):
                if len(self.results[it - 1]) > e and len(self.results[it - 2]) > e:
                    result = self.results[it - 1][e]
                    result_prev = self.results[it - 2][e]

                    if result is not None and result_prev is not None:
                        if 'QN' in result and 'QN' in result_prev:
                            QN = result['QN']
                            QN_prev = result_prev['QN']

                            if QN is not None and QN_prev is not None:
                                # Get total jobs in this layer
                                N = np.sum(self.ensemble[e].get_number_of_jobs()) if self.ensemble[e] is not None else 1
                                if N > 0:
                                    try:
                                        iter_err = np.max(np.abs(QN.flatten() - QN_prev.flatten())) / N
                                        self.maxitererr[-1] += iter_err
                                    except:
                                        pass

            if it == iter_min:
                self.averagingstart = it

            # Update relaxation factor for adaptive/auto modes (matches MATLAB lines 112-152)
            relax_mode = self.options.config.get('relax', 'none')
            if relax_mode in ['adaptive', 'auto']:
                # Track error history
                self.relax_err_history.append(self.maxitererr[-1])
                wnd = self.options.config.get('relax_history', 5)
                if len(self.relax_err_history) > wnd:
                    self.relax_err_history = self.relax_err_history[-wnd:]

                if len(self.relax_err_history) >= 3:
                    # Detect oscillation by counting sign changes in error differences
                    err = np.array(self.relax_err_history)
                    diff_err = np.diff(err)
                    if len(diff_err) >= 2:
                        sign_changes = np.sum(diff_err[:-1] * diff_err[1:] < 0)

                        if relax_mode == 'auto' and self.relax_omega == 1.0:
                            # For 'auto' mode: enable relaxation when oscillation detected
                            # (matches MATLAB converged.m line 130 - no iter_min check)
                            if sign_changes >= len(diff_err) * 0.5:
                                self.relax_omega = self.options.config.get('relax_factor', 0.1)
                                if self.options.verbose:
                                    print(f'LN: enabling relaxation, omega={self.relax_omega:.2f}')
                        elif relax_mode == 'adaptive':
                            # For 'adaptive' mode: adjust omega based on error trajectory
                            relax_min = self.options.config.get('relax_min', 0.1)
                            if sign_changes >= len(diff_err) * 0.5:
                                # Oscillating - reduce omega
                                self.relax_omega = max(relax_min, self.relax_omega * 0.8)
                            elif sign_changes == 0 and len(self.maxitererr) >= 2 and self.maxitererr[-1] < self.maxitererr[-2]:
                                # Monotonically decreasing - can increase omega slightly
                                self.relax_omega = min(1.0, self.relax_omega * 1.05)

        # Check convergence (matches MATLAB converged.m line 164: it > iter_min and check last 3 errors)
        if it > iter_min and len(self.maxitererr) >= 3:
            if (self.maxitererr[-1] < self.options.iter_tol and
                self.maxitererr[-2] < self.options.iter_tol and
                self.maxitererr[-3] < self.options.iter_tol):
                if not self.hasconverged:
                    # Reset layers and check again
                    for e in range(self.nlayers):
                        if self.ensemble[e] is not None:
                            self.ensemble[e].reset()
                    self.hasconverged = True
                else:
                    return True
        else:
            self.hasconverged = False

        return False

    def converged_stoch(self, it: int) -> bool:
        """Convergence controller for stochastic layer solvers (Robbins-Monro
        mode).

        When one or more layer solvers return noisy estimates (simulation,
        e.g. JMT/SSA/LDES, or Monte Carlo integration, e.g. NC with
        mci/imci/ls), the deterministic Picard iteration in converged() cannot
        terminate: the successive-difference error is bounded below by the
        standard error of the layer estimates, and the layer-reset
        confirmation step merely resamples the noise. This routine implements
        a stochastic approximation iteration instead:

        1. Burn-in: for the first stochiter_burnin iterations the plain
           Picard iteration runs with the relaxation factor configured at
           init.
        2. Robbins-Monro step: afterwards the relaxation factor applied by
           update_metrics to the fed-forward iterate (servt, residt, tput,
           callservt) decays as omega_k = a0/k**alpha with alpha in (0.5,1].
           Under the contraction assumption already made by the deterministic
           iteration, and zero-mean noise with bounded variance, the iterate
           converges almost surely to the true fixed point (Robbins and
           Monro, 1951). Layer seeds are rotated per iteration in pre() so
           successive evaluations observe independent noise.
        3. Polyak-Ruppert averaging: running averages of the layer results
           and of the reported iterates are maintained and installed as the
           final solution in finish(), giving the optimal O(1/sqrt(k)) rate
           and robustness to the choice of a0 (Polyak and Juditsky, 1992).
        4. Stopping: iteration stops when the drift of the averaged results
           stays below iter_tol for stochiter_conseq consecutive iterations.
           The drift of a running average decays like 1/k even under
           persistent noise, so the test terminates, and it self-calibrates:
           larger noise keeps the drift above tolerance longer, forcing more
           averaging.
        """
        if it < 1:
            return False
        burnin = int(self.options.config.get('stochiter_burnin', 5))
        a0 = float(self.options.config.get('stochiter_a0', 1.0))
        alpha = float(self.options.config.get('stochiter_alpha', 0.6))

        while len(self.maxitererr) <= it:
            self.maxitererr.append(0.0)

        # Schedule the Robbins-Monro step used by update_metrics at the next iteration
        if it >= burnin:
            self.relax_omega = min(1.0, a0 / max(1, it - burnin + 1) ** alpha)

        if it <= burnin:
            # pure Picard burn-in; no averaging or convergence testing yet
            self.maxitererr[it] = np.inf
            if self.options.verbose:
                print(f'Stochastic iteration burn-in {it}/{burnin}.')
            return False

        if self.stochiter_start is None:
            self.stochiter_start = it
            if self.options.verbose:
                print('Started Robbins-Monro averaging (stochastic layer solvers detected).')

        # Polyak-Ruppert update of the layer result averages and drift metric
        k = self.stoch_avg_count + 1
        err = 0.0
        latest = self.results[-1] if self.results else []
        for e in range(self.nlayers):
            raw = latest[e] if e < len(latest) else None
            if raw is None:
                continue
            if k == 1 or self.stoch_avg[e] is None:
                avg = {key: (None if raw.get(key) is None
                             else np.array(raw[key], dtype=float, copy=True))
                       for key in ('QN', 'UN', 'RN', 'TN', 'AN', 'WN')}
            else:
                prev = self.stoch_avg[e]
                avg = {key: _polyak_avg(prev.get(key), raw.get(key), k)
                       for key in ('QN', 'UN', 'RN', 'TN', 'AN', 'WN')}
                # drift of the averaged queue lengths, normalized by population
                N = np.sum(self.ensemble[e].get_number_of_jobs()) if self.ensemble[e] is not None else 0
                if N > 0 and avg.get('QN') is not None and prev.get('QN') is not None:
                    try:
                        d = np.abs(np.asarray(avg['QN'], dtype=float).flatten()
                                   - np.asarray(prev['QN'], dtype=float).flatten())
                        d[np.isnan(d)] = 0.0
                        err += float(np.max(d)) / N
                    except Exception:
                        pass
            self.stoch_avg[e] = avg
        self.stoch_avg_count = k

        # Polyak-Ruppert averages of the fed-forward iterates used in reporting
        if k == 1:
            self.stoch_servt_avg = None if self.servt is None else np.array(self.servt, dtype=float, copy=True)
            self.stoch_residt_avg = None if self.residt is None else np.array(self.residt, dtype=float, copy=True)
        else:
            self.stoch_servt_avg = _polyak_avg(self.stoch_servt_avg, self.servt, k)
            self.stoch_residt_avg = _polyak_avg(self.stoch_residt_avg, self.residt, k)

        self.maxitererr[it] = err
        if self.options.verbose:
            print(f'RMIterErr={err:.6e} (tol={self.options.iter_tol:.6e}, '
                  f'omega={self.relax_omega:.3f}, k={k})')

        # Stop when the averaged-iterate drift stays below tolerance
        conseq = int(self.options.config.get('stochiter_conseq', 3))
        if k > conseq:
            below = all(self.maxitererr[j] < self.options.iter_tol
                        for j in range(it - conseq + 1, it + 1))
            self.hasconverged = below
            return below
        return False

    def finish(self):
        """Operations after iterations complete (matches MATLAB finish)."""
        line_debug("LN finish: final analysis of %d layers", self.nlayers)
        # In Robbins-Monro mode, report the Polyak-Ruppert averaged results
        # rather than the last (noisy) iterate
        if self.stochiter_mode == 'rm' and self.stoch_avg_count > 0 and self.results:
            latest = self.results[-1]
            for e in range(self.nlayers):
                if e < len(latest) and latest[e] is not None and self.stoch_avg[e] is not None:
                    for key in ('QN', 'UN', 'RN', 'TN', 'AN', 'WN'):
                        latest[e][key] = self.stoch_avg[e].get(key)
            if self.stoch_servt_avg is not None:
                self.servt = self.stoch_servt_avg
            if self.stoch_residt_avg is not None:
                self.residt = self.stoch_residt_avg
        for e in range(self.nlayers):
            if self.solvers[e] is not None:
                if hasattr(self.solvers[e], 'getAvg'):
                    self.solvers[e].getAvg()
                elif hasattr(self.solvers[e], 'get_avg'):
                    self.solvers[e].get_avg()

        self.model.ensemble = self.ensemble

    def iterate(self):
        """Run iteration (matches MATLAB EnsembleSolver iterate)."""
        # Solver console: SolverLN drives an ensemble of layer models and does
        # not pass through the NetworkSolver entry point, so it opens its own
        # run here, at the method every caller reaches.
        from line_solver.api.io import console as _console
        with _console.run_scope(self, self.options):
            _console.loop('solving the layered fixed point over %d layers', self.nlayers)
            return self._iterate_body()

    def _iterate_body(self):
        from line_solver.api.io import console as _console
        line_debug("LN solver iterate starting: method=%s, nlayers=%d",
                   self.options.method if hasattr(self.options, 'method') else 'default', self.nlayers)
        it = 0
        # Raw, un-smoothed history; self.results is refreshed from it below.
        results = []
        self.results = []

        self.init()

        while not self.converged(it) and it < self.options.iter_max:
            it += 1
            self.pre(it)

            # Analyze all layers
            layer_results = []
            for e in range(self.nlayers):
                result, _ = self.analyze(it, e)
                layer_results.append(result)

            results.append(layer_results)
            # finish iteration loop.
            self.results = [list(row) for row in results]

            self.post(it)

        self.finish()

    def get_ensemble_avg(self) -> Tuple[np.ndarray, ...]:
        """Get ensemble average (matches MATLAB getEnsembleAvg)."""
        if not self.ensemble:
            return (np.array([]),) * 6

        # lang=java dispatch bypasses the native fixed point; see _kb/11-conventions-and-gotchas.md lang=java SolverLN skips iterate().
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import ln_ensemble_avg_via_jar
            return ln_ensemble_avg_via_jar(self)

        # lang=cpp solves the whole ensemble in one line-cli run and likewise never
        # enters iterate(). An absent binary is the ONLY automatic fallback: a
        # construct or option the C++ layered path refuses propagates, since
        # answering it natively would report a python number under lang='cpp'.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, ln_ensemble_avg_via_cpp
            try:
                return ln_ensemble_avg_via_cpp(self)
            except LineCliNotAvailable as e:
                line_warning("SolverLN", "lang='cpp' requested but the C++ solver is "
                             "unavailable (%s); falling back to lang='python'." % e)

        self.iterate()

        # the layers of method 'srvn.ph' carry one class per caller task, so the
        # per-element results are rebuilt analytically -- see its get_ensemble_avg
        if self._is_ph_encoding():
            return self._get_ensemble_avg_ph()

        lqn = self.lqn
        QN = np.full(lqn.nidx, np.nan)  # Queue lengths (will become utilization)
        UN = np.full(lqn.nidx, np.nan)
        RN = np.full(lqn.nidx, np.nan)
        TN = np.full(lqn.nidx, np.nan)
        PN = np.full(lqn.nidx, np.nan)  # Utilization stored here first
        SN = np.full(lqn.nidx, np.nan)  # Response time stored here first
        WN = np.full(lqn.nidx, np.nan)  # Residence time
        WN_processed = np.zeros(lqn.nidx, dtype=bool)  # Track activities already accumulated into task WN
        AN = np.full(lqn.nidx, np.nan)  # Not available yet

        E = self.nlayers

        for e in range(E):
            if len(self.results) == 0 or e >= len(self.results[-1]):
                continue
            result = self.results[-1][e]
            if result is None:
                continue

            layer = self.ensemble[e]
            if layer is None:
                continue

            client_idx = layer.attribute.get('clientIdx')
            server_idx = layer.attribute.get('serverIdx')
            source_idx = layer.attribute.get('sourceIdx')

            if server_idx is None:
                continue

            # Convert 1-based to 0-based indices
            server_idx_0 = server_idx - 1 if server_idx and server_idx >= 1 else 0
            client_idx_0 = client_idx - 1 if client_idx and client_idx >= 1 else None
            source_idx_0 = source_idx - 1 if source_idx and not np.isnan(source_idx) else None

            # Get result matrices
            result_QN = result.get('QN')
            result_UN = result.get('UN')
            result_RN = result.get('RN')
            result_TN = result.get('TN')
            result_WN = result.get('WN', result_RN)

            if result_QN is None or result_TN is None:
                continue

            # Get stations and check if ishost
            stations = layer.get_nodes()
            if server_idx_0 < len(stations):
                server_station = stations[server_idx_0]
            else:
                server_station = None

            is_host = False
            hidx = None
            if server_station is not None and hasattr(server_station, 'attribute'):
                is_host = server_station.attribute.get('ishost', False)
                hidx = server_station.attribute.get('idx')

            # For host layers, determine processor metrics, one processor at a
            # time: under flat layering the layer carries every processor.
            host_station_idx = layer.attribute.get('hostStations') or ([server_idx] if is_host else [])
            for hs in host_station_idx:
                hs0 = hs - 1 if hs >= 1 else 0
                if hs0 >= len(stations):
                    continue
                h_station = stations[hs0]
                hidx = h_station.attribute.get('idx') if hasattr(h_station, 'attribute') else None
                if hidx is None or not h_station.attribute.get('ishost', False):
                    continue
                # Aggregate metrics across all classes for processor
                if np.isnan(QN[hidx]): QN[hidx] = 0.0
                if np.isnan(PN[hidx]): PN[hidx] = 0.0

                classes = layer.get_classes()
                for c_idx, cls in enumerate(classes):
                    # Add queue length and utilization from server node
                    # Add queue length (for ALL classes)
                    if hs0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                        QN[hidx] = QN[hidx] + result_QN[hs0, c_idx]

                    # activity index extracted from the class attribute tuple.
                    if hasattr(cls, 'attribute') and cls.attribute is not None:
                        elem_type = cls.attribute[0] if len(cls.attribute) > 0 else 0
                        if elem_type == LayeredNetworkElement.ACTIVITY:
                            aidx = cls.attribute[1] if len(cls.attribute) > 1 else None
                            if aidx is not None:
                                # only the activities that run on this processor
                                if self._station_idx_of(layer, self._get_parent(self._get_parent(aidx))) != hs:
                                    continue
                                tidx = self._get_parent(aidx)  # Get parent task
                                if np.isnan(PN[aidx]): PN[aidx] = 0.0
                                if tidx is not None and np.isnan(PN[tidx]): PN[tidx] = 0.0
                                if hs0 < result_UN.shape[0] and c_idx < result_UN.shape[1]:
                                    # MATLAB does NOT apply fork_fanout correction here
                                    # (matches MATLAB getEnsembleAvg lines 55-59)
                                    util = result_UN[hs0, c_idx]
                                    PN[aidx] = PN[aidx] + util
                                    if tidx is not None:
                                        PN[tidx] = PN[tidx] + util
                                    PN[hidx] = PN[hidx] + util  # Processor utilization from ACTIVITY only

                TN[hidx] = np.nan  # Added for consistency with LQNS
            is_host = is_host or bool(layer.attribute.get('hostStations'))

            # Determine remaining metrics for all classes
            classes = layer.get_classes()
            for c_idx, cls in enumerate(classes):
                if not hasattr(cls, 'attribute') or cls.attribute is None:
                    continue

                elem_type = cls.attribute[0] if len(cls.attribute) > 0 else 0
                # under flat layering each class is served at its own station, so
                # read the layer result there rather than at the layer's serverIdx
                server_idx_0 = self._station_idx_of_class(layer, cls) - 1

                if elem_type == LayeredNetworkElement.TASK:
                    tidx = cls.attribute[1] if len(cls.attribute) > 1 else None
                    if tidx is None:
                        continue
                    if is_host:
                        # throughput read from the layer result at the client index.
                        if np.isnan(TN[tidx]):
                            if client_idx_0 is not None and client_idx_0 < result_TN.shape[0] and c_idx < result_TN.shape[1]:
                                # Get throughput from layer result at clientIdx
                                TN[tidx] = result_TN[client_idx_0, c_idx]
                    else:
                        # Task layer: get queue length (nop for utilization - matches MATLAB)
                        if server_idx_0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                            if np.isnan(QN[tidx]): QN[tidx] = 0.0
                            QN[tidx] = QN[tidx] + result_QN[server_idx_0, c_idx]

                elif elem_type == LayeredNetworkElement.ENTRY:
                    eidx = cls.attribute[1] if len(cls.attribute) > 1 else None
                    if eidx is None:
                        continue
                    # Entry response time: for phase-2 models, use residt (caller's view)
                    # which is phase-1 + overtaking correction (MATLAB getEnsembleAvg lines 84-90)
                    if (self.hasPhase2 and self.servt_ph2 is not None
                            and eidx < len(self.servt_ph2)
                            and self.servt_ph2[eidx] > 1e-8):
                        SN[eidx] = self.residt[eidx]
                    else:
                        SN[eidx] = self.servt[eidx]

                    # Entry throughput - use layer result directly (matches MATLAB getEnsembleAvg lines 90-96)
                    if is_host and client_idx_0 is not None and np.isnan(TN[eidx]):
                        if client_idx_0 < result_TN.shape[0] and c_idx < result_TN.shape[1]:
                            # Get throughput from layer result
                            # LQN throughput is total (not per-instance) - matches MATLAB
                            TN[eidx] = result_TN[client_idx_0, c_idx]

                elif elem_type == LayeredNetworkElement.ACTIVITY:
                    aidx = cls.attribute[1] if len(cls.attribute) > 1 else None
                    if aidx is None:
                        continue
                    tidx = self._get_parent(aidx)

                    # Add queue length to task (matches MATLAB line 111-112)
                    if tidx is not None:
                        if np.isnan(QN[tidx]): QN[tidx] = 0.0
                        if server_idx_0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                            QN[tidx] = QN[tidx] + result_QN[server_idx_0, c_idx]

                    # Initialize TN and QN for activity (matches MATLAB lines 113-114)
                    if np.isnan(TN[aidx]): TN[aidx] = 0.0
                    if np.isnan(QN[aidx]): QN[aidx] = 0.0

                    # Propagate activity throughput to task if task doesn't have its own class
                    # (matches MATLAB getEnsembleAvg lines 116-127)
                    if tidx is not None:
                        tasks_attr = layer.attribute.get('tasks', [])
                        has_task_class = any(t[1] == tidx for t in tasks_attr)
                        if not has_task_class:
                            if np.isnan(TN[tidx]): TN[tidx] = 0.0
                            if server_idx_0 < result_TN.shape[0] and c_idx < result_TN.shape[1]:
                                TN[tidx] = TN[tidx] + result_TN[server_idx_0, c_idx]

                    # Find entry this activity is bound to (matches MATLAB lines 130-153)
                    if tidx is not None and hasattr(lqn, 'entriesof') and lqn.entriesof is not None:
                        entries = lqn.entriesof.get(tidx, [])
                        for eidx_check in entries:
                            if (hasattr(lqn, 'graph') and lqn.graph is not None and
                                eidx_check < lqn.graph.shape[0] and aidx < lqn.graph.shape[1] and
                                lqn.graph[eidx_check, aidx] > 0):
                                if np.isnan(TN[eidx_check]): TN[eidx_check] = 0.0
                                if np.isnan(QN[eidx_check]): QN[eidx_check] = 0.0
                                if np.isnan(SN[eidx_check]): SN[eidx_check] = 0.0
                                act_tput = 0.0
                                if server_idx_0 < result_TN.shape[0] and c_idx < result_TN.shape[1]:
                                    act_tput = result_TN[server_idx_0, c_idx]
                                # Only add if entry doesn't have its own class
                                entries_attr = layer.attribute.get('entries', [])
                                has_entry_class = any(ea[1] == eidx_check for ea in entries_attr) if entries_attr else False
                                if not has_entry_class:
                                    TN[eidx_check] = TN[eidx_check] + act_tput
                                    if server_idx_0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                                        QN[eidx_check] = QN[eidx_check] + result_QN[server_idx_0, c_idx]
                                    if server_idx_0 < result_RN.shape[0] and c_idx < result_RN.shape[1]:
                                        SN[eidx_check] = SN[eidx_check] + result_RN[server_idx_0, c_idx]
                                break

                    # accumulate activity throughput; response time computed next (mirrors MATLAB lines 161-164).
                    if server_idx_0 < result_TN.shape[0] and c_idx < result_TN.shape[1]:
                        act_tput = result_TN[server_idx_0, c_idx]
                        TN[aidx] = TN[aidx] + act_tput

                    # Activity response time (matches MATLAB lines 161-164)
                    act_resp_time = 0.0
                    if server_idx_0 < result_RN.shape[0] and c_idx < result_RN.shape[1]:
                        if np.isnan(RN[aidx]): RN[aidx] = 0.0
                        act_resp_time = result_RN[server_idx_0, c_idx]
                        RN[aidx] = RN[aidx] + act_resp_time

                    if np.isnan(SN[aidx]): SN[aidx] = 0.0
                    SN[aidx] = SN[aidx] + act_resp_time

                    # Activity queue length (matches MATLAB lines 169-170)
                    if server_idx_0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                        QN[aidx] = QN[aidx] + result_QN[server_idx_0, c_idx]

                    # guard WN against NaN before use, per activity.
                    if np.isnan(WN[aidx]): WN[aidx] = 0.0
                    WN[aidx] = self.residt[aidx]
                    if tidx is not None:
                        if np.isnan(WN[tidx]): WN[tidx] = 0.0
                        if not WN_processed[aidx]:
                            WN[tidx] = WN[tidx] + self.residt[aidx]
                            WN_processed[aidx] = True

                elif elem_type == LayeredNetworkElement.CALL:
                    # Handle CALL classes (matches MATLAB getEnsembleAvg lines 99-107)
                    cidx = cls.attribute[1] if len(cls.attribute) > 1 else None
                    if cidx is not None and cidx >= 0:
                        # Get source activity from callpair
                        if hasattr(lqn, 'callpair') and lqn.callpair is not None:
                            if cidx < lqn.callpair.shape[0]:
                                # callpair column 1 is the source activity (0-indexed in the array)
                                aidx = int(lqn.callpair[cidx, 0])
                                if aidx > 0:
                                    # Check if this is a SYNC call
                                    calltype = CallType.SYNC
                                    if hasattr(lqn, 'calltype') and lqn.calltype is not None:
                                        if isinstance(lqn.calltype, np.ndarray):
                                            if cidx < len(lqn.calltype.flatten()):
                                                calltype = lqn.calltype.flatten()[cidx]
                                        elif isinstance(lqn.calltype, dict):
                                            calltype = lqn.calltype.get(cidx, CallType.SYNC)

                                    # call_mean defaults to 1.0 absent an explicit callproc entry.
                                    if calltype == CallType.SYNC:
                                        # Get call mean from callproc
                                        call_mean = 1.0
                                        if hasattr(lqn, 'callproc') and lqn.callproc is not None:
                                            if cidx < len(lqn.callproc) and lqn.callproc[cidx] is not None:
                                                proc = lqn.callproc[cidx]
                                                if hasattr(proc, 'getMean'):
                                                    call_mean = proc.getMean()
                                        # Add layer result RN * call_mean to SN[aidx]
                                        if server_idx_0 < result_RN.shape[0] and c_idx < result_RN.shape[1]:
                                            if np.isnan(SN[aidx]):
                                                SN[aidx] = 0.0
                                            SN[aidx] = SN[aidx] + result_RN[server_idx_0, c_idx] * call_mean

                                    # MATLAB getEnsembleAvg lines 106-107:
                                    # QN(aidx) = QN(aidx) + self.results{end,e}.QN(serverIdx,c)
                                    if np.isnan(QN[aidx]):
                                        QN[aidx] = 0.0
                                    if server_idx_0 < result_QN.shape[0] and c_idx < result_QN.shape[1]:
                                        QN[aidx] = QN[aidx] + result_QN[server_idx_0, c_idx]

        # entry/task throughput fallback when layer results leave them unset; mirrors MATLAB getEnsembleAvg.
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            entries = self._get_entries_of_task(tidx)

            # entry throughput falls back to its task's throughput when unset.

            # For entries without layer-derived throughput, use task throughput
            for eidx in entries:
                if np.isnan(TN[eidx]) and not np.isnan(TN[tidx]):
                    TN[eidx] = TN[tidx]

            # Entry service time = servt
            for eidx in entries:
                if np.isnan(SN[eidx]) and eidx < len(self.servt):
                    SN[eidx] = self.servt[eidx]

        # iterate activities of the task for response-time aggregation.
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            activities = self._get_activities_of_task(tidx)

            for aidx in activities:
                # Activity response time = residt (includes queueing) + call response times
                if np.isnan(SN[aidx]):
                    act_resp_time = 0.0
                    # Add activity's own residence time
                    if aidx < len(self.residt) and self.residt[aidx] > 0:
                        act_resp_time = self.residt[aidx]
                    elif aidx < len(self.servtproc) and self.servtproc[aidx] is not None:
                        # Fallback to host demand if no residt available
                        if hasattr(self.servtproc[aidx], 'getMean'):
                            act_resp_time = self.servtproc[aidx].getMean()
                        elif hasattr(self.servtproc[aidx], 'mean'):
                            act_resp_time = self.servtproc[aidx].mean
                    # accumulate activity response time from the server-node result; mirrors MATLAB line 162.

                    # Set activity response time from layer result only (not including calls)
                    # MATLAB line 162: SN(aidx) = SN(aidx) + result.RN(serverIdx,c)
                    if np.isnan(SN[aidx]):
                        SN[aidx] = 0.0
                    SN[aidx] = SN[aidx] + act_resp_time

        # Calculate entry utilization from throughput and service times
        # (matches MATLAB getEnsembleAvg lines 175-197)
        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            tidx = self._get_parent(eidx)
            if tidx is not None and np.isnan(UN[tidx]):
                UN[tidx] = 0.0
            # Phase-2 support: utilization includes both phases
            # (MATLAB getEnsembleAvg lines 187-197)
            if (self.hasPhase2 and self.servt_ph2 is not None
                    and eidx < len(self.servt_ph2)
                    and self.servt_ph2[eidx] > 1e-8
                    and not np.isnan(TN[eidx])):
                self.util_ph1[eidx] = TN[eidx] * self.servt_ph1[eidx]
                self.util_ph2[eidx] = TN[eidx] * self.servt_ph2[eidx]
                UN[eidx] = self.util_ph1[eidx] + self.util_ph2[eidx]
            elif not np.isnan(TN[eidx]) and not np.isnan(SN[eidx]):
                # Standard calculation for entries without phase-2
                UN[eidx] = TN[eidx] * SN[eidx]
            # Entry utilization = sum of activity processor utilizations for that entry
            entry_acts = lqn.actsof.get(eidx, []) if isinstance(lqn.actsof, dict) else []
            if entry_acts:
                entry_util = sum(PN[a] for a in entry_acts if not np.isnan(PN[a]))
                PN[eidx] = entry_util
            # Activity queue length (UN array) = throughput * response time
            for aidx in self._get_activities_of_task(tidx) if tidx else []:
                if not np.isnan(TN[aidx]) and not np.isnan(SN[aidx]):
                    UN[aidx] = TN[aidx] * SN[aidx]  # Queue length (throughput * response time)
                # processor/task utilization comes directly from the layer result UN, set earlier in the layer loop.
            if tidx is not None and not np.isnan(UN[eidx]):
                UN[tidx] = UN[tidx] + UN[eidx]

        # python has no Cache nodes, so cache-task metrics must be zeroed explicitly (MATLAB's Cache nodes zero them internally).

        # CacheTask: find each entry's bound activity to attribute cache metrics.
        if hasattr(lqn, 'iscache') and lqn.iscache is not None:
            for t in range(lqn.ntasks):
                tidx = lqn.tshift + t
                if lqn.iscache[tidx, 0] > 0:
                    # This is a CacheTask - find the bound activity of each entry
                    entries = self._get_entries_of_task(tidx)
                    for eidx in entries:
                        # The bound activity is the one directly connected from entry in the graph
                        for aidx in self._get_activities_of_task(tidx):
                            if lqn.graph[eidx, aidx] > 0:
                                # This is the bound activity of a cache entry - zero out
                                QN[aidx] = 0.0
                                UN[aidx] = 0.0
                                RN[aidx] = 0.0
                                TN[aidx] = 0.0
                                PN[aidx] = 0.0
                                SN[aidx] = 0.0
                                WN[aidx] = 0.0

        # AN IGNORED ELEMENT IS IDLE, NOT UNDEFINED, and the two are different
        # cells. Its component holds no reference task, so nothing reaches it and
        # every measure it HAS is zero -- but the measures its kind never has stay
        # NaN, exactly as they do for a reachable element. A flat zero over all six
        # columns broke the table's NaN mask (a processor with a queue length of 0,
        # an arrival rate reported where no solver reports one), and the mask is
        # part of the answer: see _kb/06-solver-catalog.md. Reported columns are
        # QLen=UN, Util=PN, RespT=SN, ResidT=WN, ArvR=AN, Tput=TN, so the pre-swap
        # QN and RN are discarded below and are not written here.
        for idx in range(lqn.nidx):
            if self.ignore[idx]:
                PN[idx] = 0.0        # every kind reports a utilization
                AN[idx] = np.nan     # nothing reports an arrival rate on an LQN
                kind = self._get_type(idx)
                if kind == LayeredNetworkElement.PROCESSOR:
                    UN[idx] = SN[idx] = WN[idx] = TN[idx] = np.nan
                elif kind == LayeredNetworkElement.TASK:
                    UN[idx] = WN[idx] = TN[idx] = 0.0
                    SN[idx] = np.nan
                elif kind == LayeredNetworkElement.ENTRY:
                    UN[idx] = SN[idx] = TN[idx] = 0.0
                    WN[idx] = np.nan
                elif kind == LayeredNetworkElement.ACTIVITY:
                    UN[idx] = SN[idx] = WN[idx] = TN[idx] = 0.0

        # UN=PN, RN=SN by convention (processor utilization, response time).
        final_QN = UN.copy()  # MATLAB: QN = UN (utilization in jobs)
        final_UN = PN.copy()  # MATLAB: UN = PN (processor utilization)
        final_RN = SN.copy()  # MATLAB: RN = SN (response time)

        return final_QN, final_UN, final_RN, TN, AN, WN

    def get_avg(self) -> Tuple[np.ndarray, ...]:
        """Get average metrics (alias for get_ensemble_avg)."""
        return self.get_ensemble_avg()

    def getCdfRespT(self) -> List[Optional[np.ndarray]]:
        """Response time distribution of every entry of the layered network.

        Mirrors MATLAB ``@SolverLN/getCdfRespT.m``. The distribution is formed
        by the ``moment3`` pass alone -- the mean-based update builds no law at
        all -- so a solver constructed with any other method re-runs the
        ensemble under ``moment3`` here and restores the caller's method
        afterwards. The routing layers already built serve ``moment3``
        unchanged, so only the update pass changes.

        Returns:
            A list of ``nentries`` items, one per entry in the entry-local index
            space (``lqn.eshift + i``). Each is an ``(n, 2)`` array whose
            columns are ``[F(t), t]``, the column order every CDF getter in LINE
            uses, or None for an entry the pass fitted no law to.

        Raises:
            ValueError: if the layers were built for a phase-type encoding,
                which carries no activity-graph routing to re-run over.
        """
        if not self.entrycdfrespt or self.entrycdfrespt[0] is None:
            # The distribution pass reads the routing encoding of the activity
            # graph, which the srvn.ph / flat.ph layers do not carry: re-running
            # get_avg over them would reconstruct the wrong topology rather than
            # a coarser answer. Refuse by name.
            if self._is_ph_encoding():
                raise ValueError(
                    "getCdfRespT needs the routing encoding of the activity graph, which "
                    "method='%s' does not build. Rebuild the solver with method='srvn.cs' "
                    "or method='moment3'." % self.lnmethod)
            cur_method = getattr(self.options, 'method', None)
            cur_lnmethod = self.lnmethod
            # BOTH the option and the RESOLVED method have to move: update_metrics
            # dispatches on self.lnmethod, which _build_layers resolved once, so
            # setting options.method alone leaves the mean-based update in place
            # and returns an EMPTY table.
            self.options.method = 'moment3'
            self.lnmethod = 'moment3'
            try:
                self.get_avg()
            finally:
                self.options.method = cur_method
                self.lnmethod = cur_lnmethod
        return self.entrycdfrespt

    get_cdf_resp_t = getCdfRespT

    def getTranAvg(self, *args):
        """Transient average station metrics of the layered network.

        ``options.config['ln_transient']`` selects the inter-layer coupling of
        the transient:

        - ``'decoupled'``: freeze inter-layer demands at the converged fixed
          point (get_ensemble_avg) and run each layer's transient in isolation.
        - ``'coupled'`` (default): reconcile the per-layer transients by
          waveform relaxation, so layer populations and inter-layer demands
          co-evolve in model time (getTranAvgCoupled).

        Both modes return the SAME block-diagonal layout; iteration 0 of the
        coupled relaxation is exactly the decoupled result. Mirrors MATLAB
        SolverLN.getTranAvg.
        """
        # The C++ layered path serves -a avg only, and the transient below is a
        # native computation: returning it under lang='cpp' would label python
        # numbers as C++ ones, which is what that option exists to rule out.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            raise RuntimeError(
                "lang='cpp' delegates the steady-state layered solve only (the C++ layered "
                "path implements -a avg); the layered transient is native. Use lang='python' "
                "for getTranAvg.")
        cfg = getattr(self.options, 'config', None)
        mode = None
        if cfg is not None:
            mode = cfg.get('ln_transient') if isinstance(cfg, dict) \
                else getattr(cfg, 'ln_transient', None)
        mode = (mode or 'coupled').lower()
        if mode == 'coupled':
            return self.getTranAvgCoupled(*args)
        if mode == 'decoupled':
            return self.getTranAvgDecoupled(*args)
        raise ValueError(
            "Unknown ln_transient mode '%s' (use 'coupled' or 'decoupled')." % mode)

    def _replay_init_marginal(self):
        """Re-install the warm start supplied by LayeredNetwork.initFromMarginal
        on the layer networks, and drop the layer solvers' cached results.

        The fixed-point solve hard-resets every layer when it detects
        convergence, after which the layer solver re-initializes to the default
        state, so a warm start does not survive it. The steady solve ignores the
        initial state but the transients below do not, hence the replay here.
        This is what carries the queue lengths of a SolverENV stage across an
        environment switch.

        Returns:
            True if a warm start was available and re-applied.
        """
        model = getattr(self, 'model', None)
        blocks = model.get_init_marginal_blocks() if hasattr(model, 'get_init_marginal_blocks') else None
        if not blocks or len(blocks) != self.nlayers:
            return False
        from ...api.state.marginal import roundMarginalPreservingChains
        for e in range(self.nlayers):
            layer = self.ensemble[e] if e < len(self.ensemble) else None
            if layer is None:
                continue
            solver = self.solvers[e] if e < len(self.solvers) else None
            block = np.asarray(blocks[e], dtype=float).copy()
            sname = type(solver).__name__ if solver is not None else ''
            if 'Fluid' not in sname and 'FLD' not in sname:
                # reset the layer solver so its cached steady-state result does not shadow the transient run from the replayed state.
                block = roundMarginalPreservingChains(block, layer.get_struct())
            layer.init_from_marginal(block)
            if solver is not None and hasattr(solver, 'reset'):
                # Otherwise the cached steady-state result shadows the
                # transient run from the replayed state.
                solver.reset()
        return True

    def getTranAvgDecoupled(self, *args):
        """Decoupled (frozen-demand) transient average station metrics.

        Mirrors MATLAB SolverLN.getTranAvgDecoupled: runs the ensemble fixed-point solve,
        then delegates the transient analysis to each layer solver and assembles
        the per-layer station x class traces block-diagonally (layer e in a
        disjoint row/column block). Off-block cells are left None.

        Transient traces are only produced by transient-capable layer solvers
        (Fluid, CTMC, SSA); with steady-state-only layers (MVA, NC) the
        delegated getTranAvg raises, matching the MATLAB behaviour.

        Returns:
            (QNlqn_t, UNlqn_t, TNlqn_t): each a block-diagonal nested list
            [rows][cols] of TranResult (or None off-block), where layer e
            occupies a disjoint block of rows (its stations) and columns
            (its classes).
        """
        # Run the ensemble fixed point (mirrors self.getAvg in MATLAB).
        self.get_ensemble_avg()

        # timespan validity check before replaying the initial marginal.
        ts = getattr(self.options, 'timespan', None)
        has_ts = (ts is not None and len(ts) >= 2
                  and np.all(np.isfinite(np.asarray(ts, dtype=float))))

        self._replay_init_marginal()

        # Collect the per-layer transient traces from each layer solver.
        per_layer = []
        for e in range(self.nlayers):
            solver = self.solvers[e] if e < len(self.solvers) else None
            if solver is None:
                per_layer.append(([], [], []))
                continue
            saved_ts = None
            if has_ts and hasattr(solver, 'options') and hasattr(solver.options, 'timespan'):
                saved_ts = solver.options.timespan
                solver.options.timespan = ts
            try:
                QNe, UNe, TNe = solver.getTranAvg()
            finally:
                if saved_ts is not None:
                    solver.options.timespan = saved_ts
            per_layer.append((QNe, UNe, TNe))

        # QN blocks drive the block extent, matching the MATLAB assembly.
        total_rows = sum(len(p[0]) for p in per_layer)
        total_cols = sum((len(p[0][0]) if len(p[0]) > 0 else 0) for p in per_layer)

        QNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]
        UNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]
        TNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]

        r0 = 0
        c0 = 0
        for QNe, UNe, TNe in per_layer:
            nr = len(QNe)
            nc = len(QNe[0]) if nr > 0 else 0
            for i in range(nr):
                for r in range(nc):
                    QNlqn_t[r0 + i][c0 + r] = QNe[i][r]
                    if i < len(UNe) and r < len(UNe[i]):
                        UNlqn_t[r0 + i][c0 + r] = UNe[i][r]
                    if i < len(TNe) and r < len(TNe[i]):
                        TNlqn_t[r0 + i][c0 + r] = TNe[i][r]
            r0 += nr
            c0 += nc

        return QNlqn_t, UNlqn_t, TNlqn_t

    get_tran_avg = getTranAvg
    get_tran_avg_decoupled = getTranAvgDecoupled

    # ---- Coupled layered transient (waveform relaxation) -------------------

    def getTranAvgCoupled(self, *args):
        """Coupled layered transient by waveform relaxation over the LQN ensemble.

        Port of MATLAB ``@SolverLN/getTranAvgCoupled.m``. Unlike
        getTranAvgDecoupled, which freezes inter-layer demands at the converged
        fixed point, this reconciles the per-layer transients iteratively: each
        layer's transient is driven by TIME-VARYING inter-layer demand
        trajectories taken from the other layers' latest transients, and the
        loop repeats until the trajectories stop changing (sup-norm gap over
        time). The time-varying demands are injected into each layer solver
        through the per-(station,class) rate schedule
        (``options.config['rate_sched']``), honoured by the fluid rate
        multiplier and by the CTMC time-varying transient.

        Iteration 0 uses the frozen equilibrium demands, so it reproduces
        getTranAvgDecoupled exactly; at convergence every layer relaxes to its
        fixed point, so the endpoint equals get_ensemble_avg. The return layout
        is the same block-diagonal (station x class per layer) as
        getTranAvgDecoupled.

        Coupled channels: task think times (client delay) and synchronous-call
        service demands (caller client station). Both are the dominant
        inter-layer couplings; intra-layer host service stays at its
        equilibrium value.
        """
        # timespan validity check before triggering the ensemble average solve.
        ts = getattr(self.options, 'timespan', None)
        ts = np.asarray(ts, dtype=float) if ts is not None else None
        has_ts = ts is not None and ts.size >= 2 and bool(np.all(np.isfinite(ts)))

        self.get_ensemble_avg()

        if not has_ts:
            # No finite transient horizon: nothing to co-evolve, defer to decoupled.
            return self.getTranAvgDecoupled(*args)

        E = self.nlayers
        cfg = getattr(self.options, 'config', None)
        cfg = cfg if isinstance(cfg, dict) else {}
        maxit = int(cfg.get('ln_transient_iter_max') or 20)
        tol = float(cfg.get('ln_transient_tol') or 1e-2)
        ngrid = 100
        tgrid = np.linspace(float(ts[0]), float(ts[1]), ngrid)

        # Per-layer sn (for node->station and class bookkeeping), cached once.
        layer_sn = [self.ensemble[e].getStruct() for e in range(E)]

        # Iteration 0: decoupled transients (frozen equilibrium demands already
        # set by the fixed-point solve).
        blocks, traj = self._ln_run_layers(ts, tgrid, [None] * E)

        for it in range(1, maxit + 1):
            traj_prev = traj
            # 1) recompute inter-layer demand trajectories from the latest traj
            demand = self._ln_recompute_demand(traj, tgrid)
            # 2) build the per-layer rate_sched injections from those demands
            sched_by_layer = self._ln_build_rate_sched(demand, tgrid, layer_sn)
            # 3) re-run each layer with its injected time-varying demand
            blocks, traj = self._ln_run_layers(ts, tgrid, sched_by_layer)
            # 4) convergence: sup-norm gap of the queue-length trajectories
            gap = max((float(np.max(np.abs(traj[e]['Q'] - traj_prev[e]['Q'])))
                       for e in range(E) if traj[e]['Q'].size), default=0.0)
            line_debug("LN coupled transient: iter %d, sup-norm gap %.3e", it, gap,
                       options=self.options)
            if gap < tol:
                break

        # Assemble the block-diagonal aggregate exactly as getTranAvgDecoupled does.
        total_rows = sum(len(b[0]) for b in blocks)
        total_cols = sum((len(b[0][0]) if len(b[0]) > 0 else 0) for b in blocks)
        QNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]
        UNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]
        TNlqn_t = [[None for _ in range(total_cols)] for _ in range(total_rows)]
        r0 = 0
        c0 = 0
        for QNe, UNe, TNe in blocks:
            nr = len(QNe)
            nc = len(QNe[0]) if nr > 0 else 0
            for i in range(nr):
                for r in range(nc):
                    QNlqn_t[r0 + i][c0 + r] = QNe[i][r]
                    if i < len(UNe) and r < len(UNe[i]):
                        UNlqn_t[r0 + i][c0 + r] = UNe[i][r]
                    if i < len(TNe) and r < len(TNe[i]):
                        TNlqn_t[r0 + i][c0 + r] = TNe[i][r]
            r0 += nr
            c0 += nc
        return QNlqn_t, UNlqn_t, TNlqn_t

    get_tran_avg_coupled = getTranAvgCoupled

    def _ln_run_layers(self, ts, tgrid, sched_by_layer):
        """Run each layer's transient (optionally with an injected rate_sched).

        Returns:
            (blocks, traj) where blocks[e] = (QNe, UNe, TNe) native per-layer
            handles for the block-diagonal assembly, and traj[e] holds
            (M x K x len(tgrid)) arrays Q, U, T, R resampled onto tgrid
            (R = Q/T residence via Little's law).
        """
        E = self.nlayers
        ng = len(tgrid)
        self._replay_init_marginal()
        blocks = []
        traj = []
        for e in range(E):
            solver = self.solvers[e] if e < len(self.solvers) else None
            if solver is None:
                blocks.append(([], [], []))
                traj.append({'Q': np.zeros((0, 0, ng)), 'U': np.zeros((0, 0, ng)),
                             'T': np.zeros((0, 0, ng)), 'R': np.zeros((0, 0, ng))})
                continue
            saved_ts = getattr(solver.options, 'timespan', None)
            cfg = getattr(solver.options, 'config', None)
            if not isinstance(cfg, dict):
                cfg = {}
                solver.options.config = cfg
            had_sched = 'rate_sched' in cfg
            saved_sched = cfg.get('rate_sched')
            solver.options.timespan = [float(tgrid[0]), float(tgrid[-1])]
            if sched_by_layer[e]:
                cfg['rate_sched'] = sched_by_layer[e]
            elif had_sched:
                cfg['rate_sched'] = None
            # best-effort layer solver reset, tolerating solvers without one.
            try:
                solver.reset()
            except AttributeError:
                pass
            try:
                QNe, UNe, TNe = solver.getTranAvg()
            finally:
                solver.options.timespan = saved_ts
                if had_sched:
                    cfg['rate_sched'] = saved_sched
                else:
                    cfg.pop('rate_sched', None)
            blocks.append((QNe, UNe, TNe))
            M = len(QNe)
            K = len(QNe[0]) if M > 0 else 0
            Q = np.zeros((M, K, ng))
            U = np.zeros((M, K, ng))
            T = np.zeros((M, K, ng))
            R = np.zeros((M, K, ng))
            for i in range(M):
                for r in range(K):
                    qv = self._ln_resample(QNe[i][r], tgrid)
                    uv = self._ln_resample(UNe[i][r] if i < len(UNe) and r < len(UNe[i]) else None, tgrid)
                    tv = self._ln_resample(TNe[i][r] if i < len(TNe) and r < len(TNe[i]) else None, tgrid)
                    Q[i, r, :] = qv
                    U[i, r, :] = uv
                    T[i, r, :] = tv
                    R[i, r, :] = qv / np.maximum(tv, GlobalConstants.FineTol)
            traj.append({'Q': Q, 'U': U, 'T': T, 'R': R})
        return blocks, traj

    @staticmethod
    def _ln_resample(h, tgrid):
        """Resample a transient handle (t/metric pair, or a scalar) onto tgrid."""
        ng = len(tgrid)
        if h is None:
            return np.zeros(ng)
        t = np.asarray(getattr(h, 't', []), dtype=float).ravel()
        m = np.asarray(getattr(h, 'metric', []), dtype=float).ravel()
        if t.size >= 2 and m.size == t.size:
            return np.interp(tgrid, t, m)
        if m.size >= 1:
            return np.full(ng, m[-1])
        return np.zeros(ng)

    def _ln_recompute_demand(self, traj, tgrid):
        """Recompute the time-varying inter-layer demands from the layer
        trajectories, pointwise in t, mirroring the scalar updateThinkTimes /
        updateMetricsDefault formulas.

        Returns:
            dict with 'thinkt': {tidx: (ng,) think-time trajectory} and
            'callservt': {cidx: (ng,) call service-time trajectory}.
        """
        lqn = self.lqn
        ng = len(tgrid)
        thinkt = {}
        callservt = {}

        # Task think times: from the task's own server-layer utilization/throughput.
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if np.isnan(self.idxhash[tidx]) or self._is_ref_task(tidx):
                continue
            e = int(self.idxhash[tidx])
            s_idx = self._ln_server_row(e)
            if s_idx is None or s_idx >= traj[e]['U'].shape[0]:
                continue
            uti = traj[e]['U'][s_idx, :, :].sum(axis=0)
            tti = traj[e]['T'][s_idx, :, :].sum(axis=0)
            njobs = float(np.max(self.njobs[tidx, :]))
            userthink = self._mwrbb_think_mean(tidx)
            tsafe = np.maximum(tti, GlobalConstants.FineTol)
            if self._get_sched(tidx) == SchedStrategy.INF:
                tk = (njobs - uti) / tsafe - userthink
            else:
                tk = njobs * np.abs(1.0 - uti) / tsafe - userthink
            # total mean including the user think time
            thinkt[tidx] = np.maximum(GlobalConstants.Zero, tk) + userthink

        # Synchronous-call service demands: callee entry response time * call mean.
        for cidx in range(lqn.ncalls):
            if self._get_calltype(cidx) != CallType.SYNC:
                continue
            eidx = int(lqn.callpair[cidx, 1])   # callee entry
            tidx = self._get_parent(eidx)       # callee task
            if tidx is None or np.isnan(self.idxhash[tidx]):
                continue
            e = int(self.idxhash[tidx])
            s_idx = self._ln_server_row(e)
            if s_idx is None or s_idx >= traj[e]['R'].shape[0]:
                continue
            # response time of the callee at its server, summed over the entry classes
            rc = np.zeros(ng)
            classes = self.ensemble[e].get_classes()
            for r in range(min(len(classes), traj[e]['R'].shape[1])):
                attr = getattr(classes[r], 'attribute', None)
                if attr is not None and len(attr) >= 2 \
                        and attr[0] == LayeredNetworkElement.ENTRY and int(attr[1]) == eidx:
                    rc = rc + traj[e]['R'][s_idx, r, :]
            if not np.any(rc):
                # fall back to the entry's activities response time
                rc = traj[e]['R'][s_idx, :, :].sum(axis=0)
            callservt[cidx] = rc * self._get_call_mean(cidx)

        return {'thinkt': thinkt, 'callservt': callservt}

    def _ln_server_row(self, e):
        """0-based station row of layer e's server node (None when absent)."""
        layer = self.ensemble[e]
        attr = getattr(layer, 'attribute', None)
        if not isinstance(attr, dict):
            return None
        s = attr.get('serverIdx')
        if s is None:
            return None
        return int(s) - 1

    def _ln_build_rate_sched(self, demand, tgrid, layer_sn):
        """Map the recomputed demand trajectories to per-layer rate_sched
        injections, using the same update maps updateLayers uses to place
        setService calls.

        ``options.config['ln_transient_channels']`` selects which inter-layer
        coupling channels are injected: ``'both'`` (default), ``'thinkt'``
        (client-delay only), or ``'callservt'`` (synchronous-call service only).
        Used to isolate each channel's contribution to the coupled transient.
        """
        E = self.nlayers
        lqn = self.lqn
        sched_by_layer = [[] for _ in range(E)]
        cfg = getattr(self.options, 'config', None)
        cfg = cfg if isinstance(cfg, dict) else {}
        channels = (cfg.get('ln_transient_channels') or 'both').lower()

        # think-time channel (client delay of caller tasks)
        if channels in ('both', 'thinkt') and self.thinkt_classes_updmap is not None:
            for row in np.atleast_2d(self.thinkt_classes_updmap):
                idx, aidx, nodeidx, classidx = (int(row[0]), int(row[1]), int(row[2]), int(row[3]))
                if np.isnan(self.idxhash[idx]):
                    continue
                e = int(self.idxhash[idx])
                if nodeidx != self._ln_client_node(e):
                    continue
                if self._get_type(aidx) == LayeredNetworkElement.TASK \
                        and self._get_sched(aidx) != SchedStrategy.REF \
                        and aidx in demand['thinkt']:
                    self._ln_add_sched(sched_by_layer[e], layer_sn[e], nodeidx, classidx,
                                       tgrid, demand['thinkt'][aidx])

        # call-service channel (client station of caller for each sync call)
        if channels in ('both', 'callservt') and self.call_classes_updmap is not None:
            for row in np.atleast_2d(self.call_classes_updmap):
                idx, cidx, nodeidx, classidx = (int(row[0]), int(row[1]), int(row[2]), int(row[3]))
                if np.isnan(self.idxhash[idx]):
                    continue
                e = int(self.idxhash[idx])
                if nodeidx != self._ln_client_node(e):
                    continue
                if cidx in demand['callservt']:
                    self._ln_add_sched(sched_by_layer[e], layer_sn[e], nodeidx, classidx,
                                       tgrid, demand['callservt'][cidx])

        return sched_by_layer

    def _ln_client_node(self, e):
        """1-based node index of layer e's client node."""
        layer = self.ensemble[e]
        attr = getattr(layer, 'attribute', None)
        if not isinstance(attr, dict):
            return None
        return attr.get('clientIdx')

    @staticmethod
    def _ln_add_sched(sched, sn, nodeidx, classidx, tgrid, demand):
        """Append one rate_sched entry that MODULATES the layer's equilibrium
        rate by the ratio of the transient demand to its steady-state
        (end-of-horizon) value: ``effective_rate(t) = nominal*demand(end)/demand(t)``.

        Passing ``rates = 1/demand(t)`` and ``nominal = 1/demand(end)`` makes the
        multiplier ``demand(end)/demand(t)``, which is exactly 1 at the horizon
        end, so the layer relaxes to its unmodified fixed point (endpoint ==
        get_ensemble_avg) regardless of any small mismatch between the transient
        residence Q/T and the scalar equilibrium demand. During the transient the
        ratio modulates the rate.
        """
        ist = int(sn.nodeToStation[nodeidx - 1])
        if ist < 0:
            return
        d = np.asarray(demand, dtype=float).ravel().copy()
        dend = d[-1]
        if not (dend > GlobalConstants.FineTol):
            return  # degenerate steady-state demand; skip this channel
        # clamp the rate multiplier to [dend/cap, dend*cap] so the fluid ODE's tolerated spikes do not destabilize the CTMC propagation.
        cap = 20.0
        d = np.clip(d, dend / cap, dend * cap)
        sched.append({'station': ist, 'class': classidx - 1,
                      'tgrid': np.asarray(tgrid, dtype=float).copy(),
                      'rates': 1.0 / d, 'nominal': 1.0 / dend})

    # ---- Majumdar-Woodside robust box bounds for the LQN --------------------

    # =================================================================
    # Method 'srvn.ph': the activity graph of an entry as a phase-type
    # server law.
    #
    # Each layer is a two-station cycle, Delay('Clients') + Queue(server),
    # with one closed class per caller task. The sequencing the default
    # method encodes as routing -- a class per entry, per activity and per
    # call, plus Fork, Join, Router and ClassSwitch nodes -- is composed
    # instead into a single phase-type service law per (layer, caller), by
    # the exact series-parallel reduction of Workflow.
    #
    # Twin of the MATLAB @SolverLN/buildLayersPH.m and its siblings, of the
    # JAR SolverLN *PH methods and of the C++ *_ph members of
    # solvers/ln/solver_ln.h. See _kb/06-solver-catalog.md (LN section).
    # =================================================================
    def _ph_init_state(self):
        """Allocate the per-entry law tables of method 'srvn.ph'."""
        n = self.lqn.nidx
        self._ph_wf: List[Any] = [None] * n
        self._ph_wfhost: List[Any] = [None] * n
        self._ph_execs: List[Optional[Dict[int, float]]] = [None] * n
        self._ph_callexecs: List[Optional[Dict[int, float]]] = [None] * n
        self._ph_hostalpha: List[Any] = [None] * n
        self._ph_hostT: List[Any] = [None] * n
        self._ph_hostmean = np.zeros(n)
        self._ph_entryalpha: List[Any] = [None] * n
        self._ph_entryT: List[Any] = [None] * n
        self._ph_entrymean = np.zeros(n)
        self._ph_entryscv = np.ones(n)
        self._ph_share = np.zeros(n)
        self._ph_overlap = np.ones(n)
        self._ph_setupshare = np.zeros(n)
        self._ph_xdemand = np.zeros(n)
        self._ph_ncalls = np.zeros((n, n))
        self._ph_calltime = np.zeros((n, n))
        self._ph_procresid = np.zeros(n)
        self._ph_actthinkt = np.zeros(n)
        self._ph_calltotal = np.zeros(n)
        self._ph_layer: List[Optional[PHLayer]] = [None] * (self.lqn.nhosts + self.lqn.ntasks)

    # =================================================================
    # Layer construction
    # =================================================================

    def _build_layers_ph(self, flat: bool = False):
        """Build the ensemble of a PH encoding.

        FLAT False is method 'srvn.ph': one layer per served element, each a
        two-station cycle Delay('Clients') + Queue(server), with one closed class
        per caller task. FLAT True is method 'flat.ph': ONE layer holding a
        station for every processor and every called task, with the same one
        closed class per caller task, which now visits each of the servers it
        uses once per invocation instead of meeting them through surrogate
        delays.
        """
        lqn = self.lqn
        nelem = lqn.nhosts + lqn.ntasks
        if not getattr(self, '_ph_laws_ready', False):
            self._assert_srvn_ph_supported(flat)

        # The interlock correction rewrites the populations of the call classes,
        # which this method does not create: its callers reach the server in one
        # class each. The setting is turned off on a COPY: options.config is the
        # very dict the caller passed to the constructor, so writing into it
        # rewrote the caller's own object, and a config reused across solvers
        # carried the ph decision into models that never took this method.
        try:
            if self.options.config.get('interlocking', False):
                cfg = type(self.options.config)(self.options.config)
                cfg['interlocking'] = False
                self.options.config = cfg
        except (AttributeError, TypeError):
            pass

        # A preceding probe has already composed the per-entry workflows; they do
        # not depend on the iterate, so they are not rebuilt here.
        if not getattr(self, '_ph_laws_ready', False):
            self._ph_init_laws()

        # Seed the fixed point with the static demands, then compose the entry laws
        self.residt = np.zeros(lqn.nidx)
        self.servt = np.zeros(lqn.nidx)
        self.callservt = np.zeros(lqn.ncalls)
        self.callresidt = np.zeros(lqn.ncalls)
        self.tput = np.zeros(lqn.nidx)
        self.util = np.zeros(lqn.nidx)
        self.thinkt = np.zeros(lqn.nidx)
        for aidx in range(lqn.ashift, lqn.ashift + lqn.nacts):
            self.residt[aidx] = self._ph_hostdem_mean(aidx)
        for cidx in range(lqn.ncalls):
            ct = self._ph_call_type(cidx)
            if ct in (_SYNC, _ASYNC):
                eidx = int(lqn.callpair[cidx, 1])
                v = self._ph_call_mean(cidx) * self._ph_hostmean[eidx]
                self.callservt[cidx] = v
                self.callresidt[cidx] = v
        self._ph_compose_entry_laws()

        self.ensemble = [None] * nelem
        self.solvers = [None] * nelem

        if flat:
            # ONE subnetwork holding every processor and every called task
            servers = self._build_ph_flat_layer()
            self.ensemble = [self.ensemble[0]]
            self.solvers = [self.solvers[0]]
            self.idxhash = np.full(nelem, np.nan)
            for idx in servers:
                self.idxhash[idx] = 0
            self.nlayers = 1
            self.layer_has_region = [False]
            self.layer_chains = [None]
            # every server resolves to the single flat layer, which is at once
            # the host layer and the task layer
            self.hostLayerIndices = [0]
            self.taskLayerIndices = [0]
            self._update_layers_ph(0)
            self.servt_classes_updmap = self._flatten_map(self._ph_servt_map)
            self.thinkt_classes_updmap = self._flatten_map(self._ph_thinkt_map)
            self.actthinkt_classes_updmap = self._flatten_map([[]])
            self.arvproc_classes_updmap = self._flatten_map(self._ph_arvproc_map)
            self.call_classes_updmap = self._flatten_map(self._ph_call_map)
            self.route_prob_updmap = self._flatten_map([[]])
            self.unique_route_prob_updmap = np.array([])
            return

        # One subnetwork per processor
        for hidx in range(lqn.nhosts):
            if self.ignore[hidx]:
                continue
            callers = self._ph_host_layer_callers(hidx)
            if not callers:
                continue
            self._build_ph_layer(hidx, callers, True)

        # One subnetwork per called task
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx] or self._is_ref_task(tidx):
                continue
            callers = self._ph_task_layer_callers(tidx)
            if not callers and not self._ph_async_calls_into(tidx):
                continue
            self._build_ph_layer(tidx, callers, False)

        # Compact the ensemble and index it
        empty = [i for i, e in enumerate(self.ensemble) if e is None]
        self.solvers = [sv for i, sv in enumerate(self.solvers) if i not in empty]
        self.ensemble = [e for e in self.ensemble if e is not None]
        self.idxhash = np.full(nelem, np.nan)
        layer_idx = 0
        for orig in range(nelem):
            if orig not in empty:
                self.idxhash[orig] = layer_idx
                layer_idx += 1
        self.nlayers = len(self.ensemble)
        self.layer_has_region = [False] * self.nlayers
        self.layer_chains = [None] * self.nlayers
        self.hostLayerIndices = [int(self.idxhash[h]) for h in range(lqn.nhosts)
                              if not np.isnan(self.idxhash[h])]
        self.taskLayerIndices = [int(self.idxhash[lqn.tshift + t]) for t in range(lqn.ntasks)
                              if not np.isnan(self.idxhash[lqn.tshift + t])]

        # install the initial laws, so that iteration 1 sees the seeded demands
        # rather than the placeholders the stations were created with
        self._update_layers_ph(0)

        # The maps carry no law of this method -- update_layers composes them
        # instead -- but post() resets the layers they name, so they are filled
        self.servt_classes_updmap = self._flatten_map(self._ph_servt_map)
        self.thinkt_classes_updmap = self._flatten_map(self._ph_thinkt_map)
        self.actthinkt_classes_updmap = self._flatten_map([[]])
        self.arvproc_classes_updmap = self._flatten_map(self._ph_arvproc_map)
        self.call_classes_updmap = self._flatten_map(self._ph_call_map)
        self.route_prob_updmap = self._flatten_map([[]])
        self.unique_route_prob_updmap = np.array([])

    def _build_ph_flat_layer(self) -> List[int]:
        """
        Build the ONE layer of method 'flat.ph': a client delay plus a station for
        every processor and every called task.

        A caller task is one closed class, and it visits each server it uses ONCE
        per invocation, carrying there the composed law of the demand it places on
        that server -- the same law method 'srvn.ph' installs in the server's own
        layer. What changes is that the servers now contend inside one network
        instead of seeing each other through surrogate delays, so the client delay
        keeps only the think times and whatever of the cycle this model does not
        hold. That is the whole difference between the two encodings of the PH
        composition, and it is why the reconstruction passes are shared verbatim.
        """
        from ...lang.classes import ClosedClass, OpenClass
        from ...lang.network import Network
        from ...lang.nodes import Delay, Queue, Sink, Source
        from .solver_ln import OptionsDict

        lqn = self.lqn
        servers = self._ph_flat_server_set()

        model = Network('FlatPH')
        try:
            model.setChecks(False)
        except AttributeError:
            pass
        model.attribute = OptionsDict({
            'hosts': [], 'tasks': [], 'entries': [], 'activities': [], 'calls': [],
            'clientIdx': 1, 'serverIdx': 2, 'sourceIdx': None,
            'cacheIdx': None, 'iscachelayer': False,
        })

        client_delay = Delay(model, 'Clients')

        srv = []
        station_of = {}
        model.attribute['hostStations'] = []
        model.attribute['taskStations'] = []
        server_idx_of = {}
        for idx in servers:
            ishost = idx <= lqn.nhosts
            st = Queue(model, self._get_hashname(idx), self._get_sched(idx))
            st.set_number_of_servers(self._get_nservers(idx))
            st.attribute = OptionsDict({'ishost': ishost, 'idx': idx})
            srv.append(st)
            stn = len(model.get_nodes())
            station_of[idx] = stn
            server_idx_of[idx] = stn
            if ishost:
                model.attribute['hostStations'].append(stn)
                model.attribute['hosts'].append([None, stn])
            else:
                model.attribute['taskStations'].append(stn)
                model.attribute['tasks'].append([None, stn])
        model.attribute['serverIdxOf'] = server_idx_of
        model.attribute['server_stations'] = srv
        model.attribute['nreplicas'] = 1
        # the scalar fallback of the station lookup, which no served element reaches
        model.attribute['serverIdx'] = station_of[servers[0]]

        # Callers of each server, and the union of them, which becomes the class set
        callers_of = {}
        all_callers = []
        for idx in servers:
            cs = (self._ph_host_layer_callers(idx) if idx <= lqn.nhosts
                  else self._ph_task_layer_callers(idx))
            callers_of[idx] = list(cs)
            for c in cs:
                if c not in all_callers:
                    all_callers.append(c)
        all_callers.sort()

        # One closed class per caller task
        class_of_caller = {}
        npop = 0.0
        for c in all_callers:
            # _ph_flat_server_set has refused every replicated element, so the
            # per-replica reduction the srvn builder makes is the identity here
            njobs = self._ph_layer_population(servers[0], c, 1)
            cls = ClosedClass(model, self._get_hashname(c), int(njobs), client_delay)
            cls.setReferenceClass(True)
            cls.attribute = [LayeredNetworkElement.TASK, c]
            class_of_caller[c] = cls.get_index()
            model.attribute['tasks'].append([cls.get_index(), c])
            npop += njobs
            client_delay.set_service(cls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                          self._ref_think_mean(c))))
            # A station this caller never reaches must say so with Disabled, NOT
            # with a tiny placeholder law. An FCFS station carries ONE service law
            # across its classes, so a placeholder is not inert there: it is mixed
            # into the multiserver correction and invents waiting where there is
            # none. Under 'srvn.ph' the question never arises, since every class of
            # a layer visits that layer's single server.
            for st in srv:
                st.set_service(cls, Disabled())
            for si, idx in enumerate(servers):
                if c not in callers_of[idx]:
                    continue
                srv[si].set_service(cls, Exp.fitMean(GlobalConstants.FineTol))
                self.njobs[c, idx] = njobs
                self._ph_thinkt_map[idx].append([idx, c, 1, cls.get_index()])
                self._ph_servt_map[idx].append([idx, c, station_of[idx], cls.get_index()])

        # Open classes: entry arrivals on a processor station, async calls on a task one
        open_arrivals_of = {idx: [] for idx in servers}
        source_station = None
        sink_station = None
        for si, hidx in enumerate(servers):
            # 0-based: anything at or past nhosts is a TASK, not a host. The
            # MATLAB twin (buildLayersPH.m:476) is 1-based, where `>` is right.
            if hidx >= lqn.nhosts:
                continue
            for c in callers_of[hidx]:
                # A task no other task calls has no task station, so the think-time
                # closure never gives its caller class a surrogate delay: the class
                # cycles against an Immediate one and an open stream on top of it
                # doubles the load. The chain is the representation that honours the
                # thread pool, so it is kept and closed on the arrival rate instead.
                if self._ph_open_arrival_only(c):
                    continue
                for eidx in self._ph_entries_of(c):
                    if not self._ph_has_open_arrival(eidx):
                        continue
                    if source_station is None:
                        model.attribute['sourceIdx'] = len(model.get_nodes()) + 1
                        source_station = Source(model, 'Source')
                        sink_station = Sink(model, 'Sink')
                    ocls = OpenClass(model, self._get_hashname(eidx) + '.Open', 0)
                    ocls.attribute = [LayeredNetworkElement.ENTRY, eidx]
                    source_station.set_arrival(ocls, lqn.arrival[eidx])
                    client_delay.set_service(ocls, Disabled())
                    # Disabled, not a placeholder, at every station this stream misses
                    for st in srv:
                        st.set_service(ocls, Disabled())
                    srv[si].set_service(ocls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                              self._ph_hostmean[eidx])))
                    open_arrivals_of[hidx].append((ocls.get_index(), eidx))
                    model.attribute['entries'].append([ocls.get_index(), eidx])
                    self._ph_arvproc_map[hidx].append([hidx, -eidx,
                                                       model.attribute['sourceIdx'],
                                                       ocls.get_index()])
        for si, tidx in enumerate(servers):
            if tidx <= lqn.nhosts:
                continue
            for cidx in self._ph_async_calls_into(tidx):
                if source_station is None:
                    model.attribute['sourceIdx'] = len(model.get_nodes()) + 1
                    source_station = Source(model, 'Source')
                    sink_station = Sink(model, 'Sink')
                ocls = OpenClass(model, call_hashname(lqn, cidx), 0)
                ocls.attribute = [LayeredNetworkElement.CALL, cidx]
                source_station.set_arrival(ocls, Immediate())
                client_delay.set_service(ocls, Disabled())
                # Disabled, not a placeholder, at every station this stream misses
                for st in srv:
                    st.set_service(ocls, Disabled())
                eidx = int(lqn.callpair[cidx, 1])
                srv[si].set_service(ocls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                          self._ph_entrymean[eidx])))
                open_arrivals_of[tidx].append((ocls.get_index(), -cidx))
                model.attribute['calls'].append([ocls.get_index(), cidx,
                                                 int(lqn.callpair[cidx, 0]), eidx])
                self._ph_arvproc_map[tidx].append([tidx, cidx, model.attribute['sourceIdx'],
                                                   ocls.get_index()])
                self._ph_call_map[tidx].append([tidx, cidx, station_of[tidx],
                                                ocls.get_index()])

        if source_station is not None:
            for jc in model.classes:
                if isinstance(jc, ClosedClass):
                    source_station.set_arrival(jc, Disabled())

        # Routing: one visit per server the caller uses, in server order. The number
        # of calls is carried by the service law, not by a visit ratio, so no arc
        # ever moves.
        P = model.init_routing_matrix()
        for c in all_callers:
            cls = model.classes[class_of_caller[c] - 1]
            prev = client_delay
            visited = False
            for si, idx in enumerate(servers):
                if c not in callers_of[idx]:
                    continue
                P.set(cls, cls, prev, srv[si], 1.0)
                prev = srv[si]
                visited = True
            if visited:
                P.set(cls, cls, prev, client_delay, 1.0)
        for si, idx in enumerate(servers):
            for (k, _tag) in open_arrivals_of[idx]:
                cls = model.classes[k - 1]
                P.set(cls, cls, source_station, srv[si], 1.0)
                P.set(cls, cls, srv[si], sink_station, 1.0)
        model.link(P)

        for idx in servers:
            L = PHLayer()
            L.idx = idx
            L.ishost = idx <= lqn.nhosts
            L.callers = list(callers_of[idx])
            L.class_of_caller = class_of_caller
            L.nreplicas = 1
            L.qstations = [station_of[idx]]
            L.svcmean_by_class = {}
            L.open_arrivals = list(open_arrivals_of[idx])
            L.npop = max(npop, 1.0)
            self._ph_layer[idx] = L

        self.ensemble[0] = model
        solver = self.solver_factory(model)
        self._assert_layer_solver_supports_model(solver, model, servers[0])
        self._detach_layer_config(solver)
        self._silence_layer_solver(solver)
        self.solvers[0] = solver
        return servers

    def _ph_flat_server_set(self) -> List[int]:
        """
        Processors and called tasks that become stations of the flat layer.

        The set is the elements the srvn builder would have given a layer of their
        own, so 'flat.ph' and 'srvn.ph' place the SAME stations and differ only in
        how many networks hold them. The refusals are those of _flat_server_set,
        since they are properties of the squashing and not of the encoding: each of
        these carries per-layer state that one submodel cannot hold.
        """
        lqn = self.lqn
        nelem = lqn.nhosts + lqn.ntasks
        for i in range(nelem):
            if float(lqn.repl[0, i]) > 1:
                raise ValueError("method='flat.ph' does not support replicated processors or "
                                 "tasks, whose replicas need a submodel each. Use "
                                 "method='srvn.ph'.")
        iscache = getattr(lqn, 'iscache', None)
        if iscache is not None and np.any(np.asarray(iscache).ravel()[:nelem]):
            raise ValueError("method='flat.ph' does not support cache tasks. Use "
                             "method='default'.")
        hs = getattr(lqn, 'hassetup', None)
        if hs is not None and np.any(np.asarray(hs).ravel()[:nelem]):
            raise ValueError("method='flat.ph' does not support setup tasks, whose powered-down "
                             "threads are per-layer state. Use method='srvn.ph'.")

        servers = []
        for hidx in range(lqn.nhosts):
            if self.ignore[hidx]:
                continue
            if not self._get_tasks_of_host(hidx):
                continue
            if not self._ph_host_layer_callers(hidx):
                continue
            servers.append(hidx)
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx] or self._is_ref_task(tidx):
                continue
            if not self._ph_task_layer_callers(tidx) and not self._ph_async_calls_into(tidx):
                continue
            servers.append(tidx)
        if not servers:
            raise ValueError("method='flat.ph' found no server: the model has no processor "
                             "with tasks.")
        return servers

    def _build_ph_layer(self, idx: int, callers: List[int], ishost: bool):
        """Build the two-station layer of server element IDX."""
        from ...lang.classes import ClosedClass, OpenClass
        from ...lang.network import Network
        from ...lang.nodes import Delay, Queue, Sink, Source
        from .solver_ln import OptionsDict

        lqn = self.lqn
        model = Network(self._get_hashname(idx))
        try:
            model.setChecks(False)
        except AttributeError:
            pass
        model.attribute = OptionsDict({
            'hosts': [], 'tasks': [], 'entries': [], 'activities': [], 'calls': [],
            'clientIdx': 1, 'serverIdx': 2, 'sourceIdx': None,
            'cacheIdx': None, 'iscachelayer': False,
        })

        client_delay = Delay(model, 'Clients')
        nreplicas = self._ph_replica_count(idx, callers, ishost)
        srv = []
        for m in range(1, nreplicas + 1):
            nm = self._get_hashname(idx) if m == 1 else (self._get_hashname(idx) + '.' + str(m))
            st = Queue(model, nm, self._get_sched(idx))
            st.set_number_of_servers(self._get_nservers(idx))
            st.attribute = OptionsDict({'ishost': ishost, 'idx': idx})
            srv.append(st)
        server_station_idx = len(model.get_nodes()) - nreplicas + 1
        model.attribute['serverIdxOf'] = {idx: server_station_idx}
        model.attribute['server_stations'] = srv
        model.attribute['nreplicas'] = nreplicas
        if ishost:
            model.attribute['hostStations'] = [server_station_idx]
            model.attribute['taskStations'] = []
            model.attribute['hosts'].append([None, server_station_idx])
        else:
            model.attribute['hostStations'] = []
            model.attribute['taskStations'] = [server_station_idx]
            model.attribute['tasks'].append([None, server_station_idx])

        # --- closed class per caller task
        L = PHLayer()
        L.idx, L.ishost, L.callers, L.nreplicas = idx, ishost, list(callers), nreplicas
        L.qstations = [1 + m for m in range(1, nreplicas + 1)]
        for c in callers:
            njobs = self._ph_layer_population(idx, c, nreplicas)
            self.njobs[c, idx] = njobs
            cls = ClosedClass(model, self._get_hashname(c), int(njobs), client_delay)
            cls.setReferenceClass(True)
            cls.attribute = [LayeredNetworkElement.TASK, c]
            L.class_of_caller[c] = cls.get_index()
            model.attribute['tasks'].append([cls.get_index(), c])
            client_delay.set_service(cls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                          self._ref_think_mean(c))))
            for st in srv:
                st.set_service(cls, Exp.fitMean(GlobalConstants.FineTol))
            # every layer must be refreshed after a law change: post() resets the
            # layers named by the think-time map
            self._ph_thinkt_map[idx].append([idx, c, 1, cls.get_index()])
            self._ph_servt_map[idx].append([idx, c, server_station_idx, cls.get_index()])

        # --- open classes: entry arrivals on a host layer, async calls on a task layer
        source_station = None
        sink_station = None
        if ishost:
            for c in callers:
                # A task no other task calls has no task layer, so update_think_times
                # never gives its caller class a surrogate delay: the class cycles
                # against an Immediate one and an open stream on top of it doubles the
                # load. The chain is the representation that honours the thread pool,
                # so it is kept and closed on the arrival rate instead.
                if self._ph_open_arrival_only(c):
                    continue
                for eidx in self._ph_entries_of(c):
                    if not self._ph_has_open_arrival(eidx):
                        continue
                    if source_station is None:
                        model.attribute['sourceIdx'] = len(model.get_nodes()) + 1
                        source_station = Source(model, 'Source')
                        sink_station = Sink(model, 'Sink')
                    ocls = OpenClass(model, self._get_hashname(eidx) + '.Open', 0)
                    ocls.attribute = [LayeredNetworkElement.ENTRY, eidx]
                    source_station.set_arrival(ocls, lqn.arrival[eidx])
                    client_delay.set_service(ocls, Disabled())
                    for st in srv:
                        st.set_service(ocls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                             self._ph_hostmean[eidx])))
                    L.open_arrivals.append((ocls.get_index(), eidx))
                    model.attribute['entries'].append([ocls.get_index(), eidx])
                    self._ph_arvproc_map[idx].append([idx, -eidx,
                                                   model.attribute['sourceIdx'], ocls.get_index()])
        else:
            for cidx in self._ph_async_calls_into(idx):
                if source_station is None:
                    model.attribute['sourceIdx'] = len(model.get_nodes()) + 1
                    source_station = Source(model, 'Source')
                    sink_station = Sink(model, 'Sink')
                ocls = OpenClass(model, call_hashname(lqn, cidx), 0)
                ocls.attribute = [LayeredNetworkElement.CALL, cidx]
                source_station.set_arrival(ocls, Immediate())
                client_delay.set_service(ocls, Disabled())
                eidx = int(lqn.callpair[cidx, 1])
                for st in srv:
                    st.set_service(ocls, Exp.fitMean(max(GlobalConstants.FineTol,
                                                         self._ph_entrymean[eidx])))
                L.open_arrivals.append((ocls.get_index(), -cidx))
                model.attribute['calls'].append([ocls.get_index(), cidx,
                                                 int(lqn.callpair[cidx, 0]), eidx])
                self._ph_arvproc_map[idx].append([idx, cidx, model.attribute['sourceIdx'],
                                               ocls.get_index()])
                self._ph_call_map[idx].append([idx, cidx, server_station_idx, ocls.get_index()])

        if source_station is not None:
            for jc in model.classes:
                if isinstance(jc, ClosedClass):
                    source_station.set_arrival(jc, Disabled())

        # Routing: one visit to the server per client cycle. The number of calls is
        # carried by the service law, not by a visit ratio, so no arc ever changes
        P = model.init_routing_matrix()
        for c in callers:
            cls = model.classes[L.class_of_caller[c] - 1]
            for st in srv:
                P.set(cls, cls, client_delay, st, 1.0 / nreplicas)
                P.set(cls, cls, st, client_delay, 1.0)
        for (k, _tag) in L.open_arrivals:
            cls = model.classes[k - 1]
            for st in srv:
                P.set(cls, cls, source_station, st, 1.0 / nreplicas)
                P.set(cls, cls, st, sink_station, 1.0)
        model.link(P)

        L.svcmean_by_class = {}
        np_pop = 0.0
        for c in callers:
            v = self.njobs[c, idx]
            if np.isfinite(v) and v > 0:
                np_pop += v
        L.npop = max(np_pop, 1.0)
        self._ph_layer[idx] = L
        self.ensemble[idx] = model
        solver = self.solver_factory(model)
        self._assert_layer_solver_supports_model(solver, model, idx)
        self._detach_layer_config(solver)
        self._silence_layer_solver(solver)
        self.solvers[idx] = solver

    # =================================================================
    # Composition of the entry laws
    # =================================================================

    def _ph_init_laws(self):
        """Build the per-entry workflows and the iteration-invariant processor law."""
        lqn = self.lqn
        nelem = lqn.nhosts + lqn.ntasks
        self._ph_servt_map = [[] for _ in range(nelem)]
        self._ph_thinkt_map = [[] for _ in range(nelem)]
        self._ph_arvproc_map = [[] for _ in range(nelem)]
        self._ph_call_map = [[] for _ in range(nelem)]

        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            tidx = int(lqn.parent[eidx, 0])
            if self.ignore[tidx]:
                continue
            ew = entry_workflow(self.model, lqn, eidx, True)
            self._ph_wf[eidx] = ew.wf
            self._ph_execs[eidx] = ew.execs
            self._ph_callexecs[eidx] = ew.callexecs
            eh = entry_workflow(self.model, lqn, eidx, False)
            self._ph_wfhost[eidx] = eh.wf
            # the processor sees the WORK of concurrent branches, not their elapsed
            # time, so the host law serialises an AND fork -- see serial_law
            alpha, T = serial_law(self._ph_wfhost[eidx])
            self._ph_hostalpha[eidx] = alpha
            self._ph_hostT[eidx] = T
            self._ph_hostmean[eidx] = ph_moments(alpha, T)[0]

        # until the first iteration reports throughputs, a task splits its
        # requests evenly over its entries
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            entries = self._ph_entries_of(tidx)
            for eidx in entries:
                self._ph_share[eidx] = 1.0 / len(entries)

    def _ph_compose_entry_laws(self):
        """
        Recompose the entry service laws from the current fixed-point iterate.

        The composed mean is NOT the sum of the leaf means when the graph forks:
        the branches of an AND fork overlap, and the entry finishes with the last
        of them. The ratio of the two, the overlap factor, is what the
        caller-side aggregates are scaled by, so that the pieces of a cycle still
        add up to the cycle.
        """
        lqn = self.lqn
        entry_setup_share = np.zeros(lqn.nidx)
        self._ph_overlap[:] = 1.0

        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            if self._ph_wf[eidx] is None:
                continue
            w = self._ph_wf[eidx]
            ex = self._ph_execs[eidx]
            entrysum = 0.0
            procsum = 0.0
            for aidx in self._ph_acts_of(eidx):
                m = self.residt[aidx] + self._ph_act_think_mean(aidx)
                procsum += ex[aidx] * m
                w.setActivityDemandMean(self._ph_name(aidx), max(m, GlobalConstants.FineTol))
                for cidx in lqn.callsof.get(aidx, []):
                    if self._ph_call_type(cidx) != _SYNC:
                        continue
                    w.setActivityDemand(call_hashname(lqn, cidx), self._ph_call_burst_law(cidx))
                    m += self.callservt[cidx]
                entrysum += ex[aidx] * m
            alpha, T = w.refreshPH()
            alpha = np.asarray(alpha, dtype=float).reshape(1, -1)
            T = np.asarray(T, dtype=float)
            m1, scv = ph_moments(alpha, T)
            # All activities of an entry run on ONE processor, so the branches of an
            # AND fork cannot overlap the processor residence they request: the
            # composed maximum is a lower bound on the entry service time only above
            # that total. Where it falls below, the law is rescaled in time to it,
            # which keeps its shape, its SCV and its order.
            if procsum > m1 + GlobalConstants.FineTol:
                T = T * (m1 / procsum)
                m1 = procsum
            # A SetupTask powers a thread down when it goes idle, so a request may
            # find it off and pay a cold start before the entry runs at all. The
            # setup is not part of the activity graph and never enters the
            # series-parallel reduction: it is prefixed to the composed law
            # afterwards, as the mixture p*(setup THEN entry) + (1-p)*entry, which
            # is again phase-type. See _setup_prob for p.
            p = self._ph_setup_prob(eidx)
            if p > GlobalConstants.FineTol:
                sl = self._ph_setup_law(int(lqn.parent[eidx, 0]))
                if sl is not None:
                    ac, Tc = Workflow._composeSerial(sl[0], sl[1], alpha, T)
                    alpha, T = Workflow._composeMixture([ac, alpha], [Tc, T],
                                                        np.array([p, 1 - p]))
                    alpha = np.asarray(alpha, dtype=float).reshape(1, -1)
                    T = np.asarray(T, dtype=float)
                    m1, scv = ph_moments(alpha, T)
                    # The share of the entry law that is cold start and not work. The
                    # surrogate-delay closure measures a thread's cycle in WORK, so it
                    # must not read a station utilization that this has inflated --
                    # see update_think_times.
                    entry_setup_share[eidx] = (
                        p * self._setup_dist_mean(getattr(lqn, 'setuptime', None),
                                               int(lqn.parent[eidx, 0]))
                        / max(m1, GlobalConstants.FineTol))
            self._ph_entryalpha[eidx] = alpha
            self._ph_entryT[eidx] = T
            self._ph_entrymean[eidx] = m1
            self._ph_entryscv[eidx] = scv
            if entrysum > GlobalConstants.FineTol:
                self._ph_overlap[eidx] = min(1.0, m1 / entrysum)

        # Per task, the share-weighted fraction of its station service that is
        # cold start rather than work.
        self._ph_setupshare[:] = 0.0
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx]:
                continue
            for eidx in self._ph_entries_of(tidx):
                self._ph_setupshare[tidx] += self._ph_share[eidx] * entry_setup_share[eidx]

        # Expected number of calls per invocation, and the caller-side aggregates
        self._ph_ncalls[:] = 0.0
        self._ph_calltime[:] = 0.0
        self._ph_procresid[:] = 0.0
        self._ph_actthinkt[:] = 0.0
        self._ph_calltotal[:] = 0.0
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx]:
                continue
            for eidx in self._ph_entries_of(tidx):
                w = self._ph_share[eidx]
                if w <= 0 or self._ph_execs[eidx] is None:
                    continue
                ex = self._ph_execs[eidx]
                r = self._ph_overlap[eidx]
                for aidx in self._ph_acts_of(eidx):
                    self._ph_procresid[tidx] += w * r * ex[aidx] * self.residt[aidx]
                    self._ph_actthinkt[tidx] += w * r * ex[aidx] * self._ph_act_think_mean(aidx)
                    for cidx in lqn.callsof.get(aidx, []):
                        if self._ph_call_type(cidx) != _SYNC:
                            continue
                        tgte = int(lqn.callpair[cidx, 1])
                        tgtt = int(lqn.parent[tgte, 0])
                        # the COUNT of calls does not change with the overlap, only
                        # the time the caller is held by them
                        self._ph_ncalls[tidx, tgte] += w * ex[aidx] * self._ph_call_mean(cidx)
                        self._ph_calltime[tidx, tgtt] += w * r * ex[aidx] * self.callservt[cidx]
                        self._ph_calltotal[tidx] += w * r * ex[aidx] * self.callservt[cidx]

    def _ph_call_burst_law(self, cidx: int):
        """
        Law of the total time one execution of the issuing activity spends in call
        CIDX: the geometric compound, of mean callproc_mean, of the response law of
        the called entry. The response law is fitted to the response time reported
        by the callee's layer and to the SCV of the callee's own composed law, so
        no extra solver output is needed.
        """
        m = self._ph_call_mean(cidx)
        eidx = int(self.lqn.callpair[cidx, 1])
        if m <= GlobalConstants.FineTol:
            return Immediate()
        R = self.callservt[cidx] / m
        scv = self._ph_entryscv[eidx]
        if not np.isfinite(scv) or scv <= GlobalConstants.FineTol:
            scv = 1.0
        base = APH.fitMeanAndSCV(max(R, GlobalConstants.FineTol), scv)
        alpha, T = Workflow._composeLoopGeometric(_alpha_of(base), _subgen_of(base), m)
        if Workflow.isAcyclicGenerator(T):
            return APH(alpha, T)
        return PH(alpha, T)

    # =================================================================
    # Pushing the composed laws into the layers
    # =================================================================

    def _update_layers_ph(self, it: int):
        """
        Push the composed laws into the layers.

        A layer of this method carries no routing that depends on the iterate: the
        number of calls a caller makes is folded into its service law rather than
        into a visit ratio, so only two laws move per (layer, class) -- the
        phase-type service law at the server and the mean of the surrogate delay
        at the client.
        """
        lqn = self.lqn
        for idx in range(lqn.nhosts + lqn.ntasks):
            if np.isnan(self.idxhash[idx]) or self._ph_layer[idx] is None:
                continue
            L = self._ph_layer[idx]
            model = self.ensemble[int(self.idxhash[idx])]
            stations = model.get_stations()
            client_delay = stations[0]

            for c in L.callers:
                k = L.class_of_caller[c]
                cls = model.classes[k - 1]
                alpha, T = self._ph_service_law(idx, L.ishost, c)
                L.svcmean_by_class[k] = ph_moments(alpha, T)[0]
                law = self._ph_station_law(alpha, T)
                for st in L.qstations:
                    stations[st - 1].set_service(cls, law)
                client_delay.set_service(cls, Exp.fitMean(
                    max(GlobalConstants.FineTol, self._ph_delay_mean(idx, c))))

            for (k, tag) in L.open_arrivals:
                cls = model.classes[k - 1]
                if tag > 0:
                    # entry arrival: the processor demand law of the entry is static
                    L.svcmean_by_class[k] = self._ph_hostmean[tag]
                    continue
                cidx = -tag
                eidx = int(lqn.callpair[cidx, 1])
                L.svcmean_by_class[k] = self._ph_entrymean[eidx]
                law = self._ph_station_law(self._ph_entryalpha[eidx], self._ph_entryT[eidx])
                for st in L.qstations:
                    stations[st - 1].set_service(cls, law)
                aidx = int(lqn.callpair[cidx, 0])
                rate = self.tput[aidx] * self._ph_call_mean(cidx)
                if not np.isfinite(rate) or rate <= GlobalConstants.FineTol:
                    rate = GlobalConstants.FineTol
                model.get_nodes()[model.attribute['sourceIdx'] - 1].set_arrival(
                    cls, Exp.fitRate(rate))

    def _ph_service_law(self, idx: int, ishost: bool, c: int):
        """Law of the demand caller C places on the server of layer IDX per invocation."""
        lqn = self.lqn
        if ishost:
            # mixture over the entries of C, weighted by their share of its requests
            alphas, Ts, probs = [], [], []
            for eidx in self._ph_entries_of(c):
                if self._ph_hostT[eidx] is None or self._ph_share[eidx] <= 0:
                    continue
                alphas.append(self._ph_hostalpha[eidx])
                Ts.append(self._ph_hostT[eidx])
                probs.append(self._ph_share[eidx])
            if not alphas:
                return self._ph_immediate_law()
            probs = np.asarray(probs, dtype=float)
            probs = probs / probs.sum()
            return Workflow._composeMixture(alphas, Ts, probs)

        # task layer: the total demand is the sum, over the entries of the server, of
        # a geometric compound of the entry law of mean equal to the number of calls
        alpha = None
        T = None
        for eidx in self._ph_entries_of(idx):
            n = self._ph_ncalls[c, eidx]
            if n <= GlobalConstants.FineTol or self._ph_entryT[eidx] is None:
                continue
            a2, T2 = Workflow._composeLoopGeometric(self._ph_entryalpha[eidx], self._ph_entryT[eidx], n)
            if alpha is None:
                alpha, T = a2, T2
            else:
                alpha, T = Workflow._composeSerial(alpha, T, a2, T2)
        if alpha is None:
            return self._ph_immediate_law()
        return alpha, T

    @staticmethod
    def _ph_immediate_law():
        return np.array([[1.0]]), np.array([[-GlobalConstants.Immediate]])

    def _ph_delay_mean(self, idx: int, c: int) -> float:
        """
        Mean time a thread of caller C spends away from the stations of the model
        that holds server IDX, per invocation: idle, plus whatever of its cycle
        that model does not hold as a station of its own.

        This is ONE closure for both layerings. Under 'srvn.ph' the model holds a
        single server, so a host layer charges the whole call burst to the delay
        and a task layer charges the caller's processor plus every other callee.
        Under 'flat.ph' the model holds every server, and only the think times
        are left.

        Every term is SUMMED in rather than obtained by subtracting from a total.
        That subtraction cancels catastrophically once a call time is large: a
        caller whose only callee is this server has the two terms equal, and
        7 + 1.4e47 - 1.4e47 is 0, not 7, because the think time falls below the
        ULP of the call time. The layer then sees a client delay of zero,
        saturates, reports a residence time that inflates the very call time that
        caused the cancellation, and the fixed point runs away -- lqn_sockshop
        reached RespT 1.4e47 this way.
        """
        lqn = self.lqn
        z = self.thinkt[c] + self._ref_think_mean(c)
        if not np.isfinite(z) or z < 0:
            z = 0.0
        z += self._ph_actthinkt[c]
        # the caller's own processor residence, unless this model holds it
        hidx = int(lqn.parent[c, 0]) if c < len(lqn.parent) else -1
        if not self._ph_served_here(idx, hidx):
            z += self._ph_procresid[c]
        # and the time spent at every callee this model does not hold
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self._ph_served_here(idx, tidx):
                continue
            z += float(self._ph_calltime[c][tidx])
        if not np.isfinite(z) or z < 0:
            z = GlobalConstants.FineTol
        return float(z)

    def _ph_served_here(self, idx: int, elem: int) -> bool:
        """
        True when LQN element ELEM is a station of the same model that holds
        server IDX. Under 'srvn.ph' that is ELEM == IDX, since each server has a
        layer of its own; under 'flat.ph' it is every server of the one network.
        """
        if elem is None or elem < 0 or elem >= len(self.idxhash):
            return False
        if idx < 0 or idx >= len(self.idxhash):
            return False
        a, b = self.idxhash[elem], self.idxhash[idx]
        return bool(np.isfinite(a) and np.isfinite(b) and a == b)

    def _ph_station_law(self, alpha, T):
        """
        Station law of a composed workflow. A geometric loop over a body of two or
        more phases closes a cycle in the phase graph, and a cyclic generator is a
        PH and not an APH: no layer solver declares PH, so such a law is reduced to
        the APH with the SAME first two moments. AMVA and NC read exactly those
        two, so the reduction is lossless for them and is a two-moment fit for the
        phase-aware layer solvers.
        """
        if Workflow.isAcyclicGenerator(T):
            return APH(alpha, T)
        m1, scv = ph_moments(alpha, T)
        return APH.fitMeanAndSCV(m1, scv)

    # =================================================================
    # Metric reconstruction
    # =================================================================

    def _update_metrics_ph(self, it: int):
        """
        Reconstruct the LQN metrics.

        A layer of this method reports one row per caller task, not one per entry,
        activity and call, so the per-element quantities the rest of SolverLN reads
        -- servt, residt, callservt, callresidt, tput -- are recovered analytically
        from the series-parallel weights of the entry workflows.

        The split is conservative by construction. A station reports a residence
        time R per visit against a service law of mean S, so the queueing inflation
        R/S is attributed to every leaf of that visit in proportion to its own
        mean: the pieces sum back to R exactly.
        """
        lqn = self.lqn
        n = lqn.nidx
        self.servt = np.zeros(n)
        self.residt = np.zeros(n)
        self.callservt = np.zeros(lqn.ncalls)
        self.callresidt = np.zeros(lqn.ncalls)

        infl_num = np.zeros(n)
        infl_den = np.zeros(n)
        task_tput = np.zeros(n)
        open_tput = np.zeros(n)

        # Host layers: the queueing inflation of the processor demand
        for hidx in range(lqn.nhosts):
            if np.isnan(self.idxhash[hidx]) or self._ph_layer[hidx] is None:
                continue
            L = self._ph_layer[hidx]
            res = self.results[-1][int(self.idxhash[hidx])]
            if not res:
                continue
            npop = self._ph_layer_pop(L, hidx)
            for c in L.callers:
                k = L.class_of_caller[c]
                kc = k - 1  # result matrices index classes from zero
                X = self._ph_sum_over(res['TN'], L.qstations, kc)
                R = self._ph_residence(self._ph_sum_over(res['QN'], L.qstations, kc), X,
                                    res['RN'][L.qstations[0] - 1, kc])
                f = self._ph_inflation_of(R, L.svcmean_by_class.get(k, 0.0), npop)
                if not np.isfinite(X) or X < 0:
                    X = 0.0
                # TOTAL over the replicas. The processor layer of a replicated element
                # models ONE representative replica, so X is one replica's rate and the
                # element's own rate is REPL times it. The matching per replica quantity
                # is xdemand, which the think-time closure divides down for the same reason.
                task_tput[c] += self._ph_repl(c) * X
                for eidx in self._ph_entries_of(c):
                    w = max(self._ph_share[eidx], 0.0) * X
                    infl_num[eidx] += w * f
                    infl_den[eidx] += w
            for (k, tag) in L.open_arrivals:
                if tag <= 0:
                    continue  # an async call is served in the task layer, not here
                eidx = tag
                kc = k - 1
                X = self._ph_sum_over(res['TN'], L.qstations, kc)
                if not np.isfinite(X) or X <= 0:
                    continue
                f = self._ph_inflation_of(
                    self._ph_residence(self._ph_sum_over(res['QN'], L.qstations, kc), X,
                                    res['RN'][L.qstations[0] - 1, kc]),
                    L.svcmean_by_class.get(k, 0.0), npop)
                infl_num[eidx] += X * f
                infl_den[eidx] += X
                open_tput[eidx] += X
                task_tput[int(lqn.parent[eidx, 0])] += X

        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            f = 1.0
            if infl_den[eidx] > GlobalConstants.FineTol:
                f = infl_num[eidx] / infl_den[eidx]
            if not np.isfinite(f) or f < 1:
                f = 1.0  # a residence time cannot fall below the demand it contains
            for aidx in self._ph_acts_of(eidx):
                self.residt[aidx] = f * self._ph_hostdem_mean(aidx)

        # Task layers: the response time of every call
        relw = np.zeros(n)
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if np.isnan(self.idxhash[tidx]) or self._ph_layer[tidx] is None:
                continue
            L = self._ph_layer[tidx]
            res = self.results[-1][int(self.idxhash[tidx])]
            if not res:
                continue
            npop = self._ph_layer_pop(L, tidx)
            for c in L.callers:
                k = L.class_of_caller[c]
                kc = k - 1
                X = self._ph_sum_over(res['TN'], L.qstations, kc)
                if not np.isfinite(X) or X < 0:
                    X = 0.0
                g = self._ph_inflation_of(
                    self._ph_residence(self._ph_sum_over(res['QN'], L.qstations, kc), X,
                                    res['RN'][L.qstations[0] - 1, kc]),
                    L.svcmean_by_class.get(k, 0.0), npop)
                for cidx in self._ph_sync_calls_between(c, tidx):
                    eidx = int(lqn.callpair[cidx, 1])
                    v = self._ph_call_mean(cidx) * g * self._ph_entrymean[eidx]
                    self.callservt[cidx] = v
                    self.callresidt[cidx] = v
                for eidx in self._ph_entries_of(tidx):
                    relw[eidx] += X * self._ph_ncalls[c, eidx]
            for (k, tag) in L.open_arrivals:
                if tag >= 0:
                    continue
                cidx = -tag
                eidx = int(lqn.callpair[cidx, 1])
                X = self._ph_sum_over(res['TN'], L.qstations, k - 1)
                R = res['RN'][L.qstations[0] - 1, k - 1]
                if np.isfinite(R) and R > 0:
                    v = R * self._ph_call_mean(cidx)
                    self.callservt[cidx] = v
                    self.callresidt[cidx] = v
                if np.isfinite(X) and X > 0:
                    relw[eidx] += X

        # Entry shares and throughputs
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            entries = self._ph_entries_of(tidx)
            if not entries:
                continue
            # How the requests SPLIT over the entries is a flow-balance question, and
            # is answered at the task layer: a caller class reaches that server once
            # per invocation of the caller, carrying its whole call burst in its
            # service law, so the station rate counts caller cycles and the per-entry
            # rate is that rate times the calls the caller makes.
            tot = sum(relw[eidx] + open_tput[eidx] for eidx in entries)
            if tot > GlobalConstants.FineTol:
                for eidx in entries:
                    self._ph_share[eidx] = (relw[eidx] + open_tput[eidx]) / tot
            else:
                for eidx in entries:
                    self._ph_share[eidx] = 1.0 / len(entries)
            # HOW MANY requests the task completes is a different question, and the
            # flow-balance total does not answer it: that total is what the callers
            # DEMAND, not what the task's threads can deliver. A thread cycles through
            # its host demand AND then through the task think time, and only the
            # processor layer of the task carries both, so the rate is read there.
            if task_tput[tidx] > GlobalConstants.FineTol:
                self.tput[tidx] = task_tput[tidx]
            else:
                self.tput[tidx] = tot  # no processor layer of its own
            for eidx in entries:
                self.tput[eidx] = self.tput[tidx] * self._ph_share[eidx]
            # The DEMAND is kept apart because it, and not the rate just reported, is
            # what closes the surrogate delay: normalising the think time by a rate the
            # same think time produced makes the processor layer self-referential and it
            # settles wherever it started -- see update_think_times. PER REPLICA,
            # because the thread count it is paired with there is per replica.
            nrep = max(1.0, self._ph_repl(tidx))
            self._ph_xdemand[tidx] = (tot / nrep) if tot > GlobalConstants.FineTol \
                else (self.tput[tidx] / nrep)

        # Recovery, under-relaxation, and the derived per-element quantities
        omega = self.relax_omega
        for aidx in range(lqn.ashift, lqn.ashift + lqn.nacts):
            v = self.residt[aidx]
            if (not np.isfinite(v)) and it > 1 and np.isfinite(self.residt_prev[aidx]):
                v = self.residt_prev[aidx]
            if omega < 1.0 and it > 1 and not np.isnan(self.residt_prev[aidx]):
                v = omega * v + (1 - omega) * self.residt_prev[aidx]
            self.residt[aidx] = v
            self.residt_prev[aidx] = v
        for cidx in range(lqn.ncalls):
            v = self.callservt[cidx]
            if not np.isfinite(v):
                v = self.callservt_prev[cidx] if (it > 1 and np.isfinite(self.callservt_prev[cidx])) else 0.0
            if omega < 1.0 and it > 1 and not np.isnan(self.callservt_prev[cidx]):
                v = omega * v + (1 - omega) * self.callservt_prev[cidx]
            self.callservt[cidx] = v
            self.callresidt[cidx] = v
            self.callservt_prev[cidx] = v
            self.callresidt_prev[cidx] = v
            if v > 0:
                self.callservtproc[cidx] = Exp.fitMean(v)

        # Recompose the entry laws from the iterate just computed. The entry service
        # time is then the mean of the COMPOSED law and not the sum of the parts: the
        # branches of an AND fork overlap, so an entry that forks finishes with the
        # last of its branches and is not charged their sum.
        self._ph_compose_entry_laws()

        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            if self._ph_execs[eidx] is None:
                continue
            ex = self._ph_execs[eidx]
            for aidx in self._ph_acts_of(eidx):
                sa = self.residt[aidx] + self._ph_act_think_mean(aidx)
                for cidx in lqn.callsof.get(aidx, []):
                    if self._ph_call_type(cidx) == _SYNC:
                        sa += self.callservt[cidx]
                self.servt[aidx] = sa
                self.servt_prev[aidx] = sa
                self.tput[aidx] = self.tput[eidx] * ex[aidx]
                self.tput_prev[aidx] = self.tput[aidx]
                # Exp rejects rate 0 here where MATLAB admits it and the JAR clamps it; a null rate is Disabled, as in the default path.
                self.tputproc[aidx] = Exp.fitRate(self.tput[aidx]) \
                    if self.tput[aidx] > 0 else Disabled()
                if sa > 0:
                    self.servtproc[aidx] = Exp.fitMean(sa)
            self.servt[eidx] = self._ph_entrymean[eidx]
            self.residt[eidx] = self._ph_entrymean[eidx]
            if self.servt[eidx] > 0:
                self.servtproc[eidx] = Exp.fitMean(self.servt[eidx])
            # published for inspection, as MATLAB's SolverLN carries entry_servt
            # on the object; this encoding resolves it entry by entry rather
            # than by one servtmatrix solve, so it is filled in here
            if self.entry_servt is None or len(self.entry_servt) < lqn.nidx:
                self.entry_servt = np.zeros(lqn.nidx)
            self.entry_servt[eidx] = self.servt[eidx]

    # =================================================================
    # Surrogate delays
    # =================================================================

    def _update_think_times_ph(self, it: int):
        """
        Surrogate delay of every caller.

        Same closure as update_think_times -- a thread of the task is idle for
        whatever of its cycle the task's own station does not hold -- but the rate
        it is normalised by is the INVOCATION rate of the task and not the
        throughput of its station. Under this method a caller class reaches the
        server once per invocation of the caller, carrying its whole call burst in
        its service law, so the station rate counts caller cycles rather than calls
        and the two differ by the mean number of calls.
        """
        lqn = self.lqn
        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx]:
                continue
            # only a reference task's think time separates one request from the next;
            # on a served task it is not a per-request delay -- see _ref_think_mean
            ztask = self._ref_think_mean(tidx)
            if np.isnan(self.idxhash[tidx]):
                # A task no other task calls but whose entries carry an arrival still
                # has a cycle: its threads are driven by the stream. build_layers drops
                # the open class for it precisely so this closure can set the rate.
                arvrate = self._ph_arrival_rate(tidx)
                if arvrate > GlobalConstants.FineTol:
                    njobs = self._ph_maxmult(tidx)
                    if not np.isfinite(njobs) or njobs <= 0:
                        njobs = float(np.max(self.njobs[tidx, :]))
                    z = max(GlobalConstants.Zero,
                            njobs / arvrate - self._ph_host_resid(tidx) - ztask)
                    om = self.relax_omega
                    if om < 1.0 and it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                        z = om * z + (1 - om) * self.thinkt_prev[tidx]
                    self.tput[tidx] = arvrate
                    self.thinkt[tidx] = z
                    self.thinkt_prev[tidx] = z
                    self.thinktproc[tidx] = Exp.fitMean(z + ztask)
                    continue
                # a reference task, or one no other task calls: it has no station of
                # its own, so its only delay is the think time the user declared
                self.thinkt[tidx] = GlobalConstants.FineTol
                self.thinktproc[tidx] = Immediate()
                continue
            L = self._ph_layer[tidx]
            res = self.results[-1][int(self.idxhash[tidx])]
            U = float(np.nansum(res['UN'][L.qstations[0] - 1, :])) if res else 0.0
            self.util[tidx] = U
            # The closure below measures a thread's cycle in WORK: it is idle for
            # whatever of the cycle its station does not hold it working. A SetupTask's
            # station service also carries a cold start, which is time the thread is
            # unavailable but is not work, so it is taken back out of U before the
            # closure reads it. Zero for every task without a setup.
            if self._ph_setupshare[tidx] > 0:
                U = U * (1 - self._ph_setupshare[tidx])
            # the rate the CALLERS ask of the task, not the rate its processor layer
            # reported: the latter is itself a function of this think time
            X = self._ph_xdemand[tidx]
            if not (X > GlobalConstants.FineTol):
                X = self.tput[tidx]
            # The thread pool of ONE replica, the convention xdemand is kept in.
            njobs = self._ph_maxmult(tidx)
            if not np.isfinite(njobs) or njobs <= 0:
                njobs = float(np.max(self.njobs[tidx, :]))
            if X > GlobalConstants.FineTol:
                if self._get_sched(tidx) == SchedStrategy.INF:
                    # an infinite server reports a mean number of busy threads
                    z = (njobs - U) / X - ztask
                else:
                    z = njobs * abs(1 - U) / X - ztask
            else:
                z = self.thinkt[tidx]
            z = max(GlobalConstants.Zero, z)
            if it > 1 and not np.isnan(self.thinkt_prev[tidx]) and not np.isfinite(z):
                z = self.thinkt_prev[tidx]
            omega = self.relax_omega
            if omega < 1.0 and it > 1 and not np.isnan(self.thinkt_prev[tidx]):
                z = omega * z + (1 - omega) * self.thinkt_prev[tidx]
            self.thinkt[tidx] = z
            self.thinkt_prev[tidx] = z
            self.thinktproc[tidx] = Exp.fitMean(z + ztask)

    def _ph_arrival_rate(self, tidx: int) -> float:
        """
        Total exogenous rate into the entries of task TIDX, zero unless the arrival
        is the only way in -- the predicate build_layers drops the open class on.
        """
        lqn = self.lqn
        if self._is_ref_task(tidx):
            return 0.0
        for eidx in self._ph_entries_of(tidx):
            if self._ph_any_caller_of(eidx):
                return 0.0
        rate = 0.0
        for eidx in self._ph_entries_of(tidx):
            d = lqn.arrival.get(eidx) if isinstance(lqn.arrival, dict) else None
            if d is None:
                continue
            try:
                m = float(d.getMean())
            except (AttributeError, TypeError, ValueError):
                continue
            if np.isfinite(m) and m > GlobalConstants.FineTol:
                rate += 1.0 / m
        return rate

    def _ph_host_resid(self, tidx: int) -> float:
        """Response time the caller class of task TIDX sees at its processor layer."""
        lqn = self.lqn
        hidx = int(lqn.parent[tidx, 0])
        if hidx < 0 or hidx >= len(self.idxhash) or np.isnan(self.idxhash[hidx]):
            return 0.0
        L = self._ph_layer[hidx]
        if L is None or tidx not in L.class_of_caller:
            return 0.0
        res = self.results[-1][int(self.idxhash[hidx])]
        if not res:
            return 0.0
        r = res['RN'][L.qstations[0] - 1, L.class_of_caller[tidx] - 1]
        return 0.0 if np.isnan(r) else float(r)

    # =================================================================
    # Result reconstruction
    # =================================================================

    def _get_ensemble_avg_ph(self):
        """
        LQN-level results. The layers report per caller task, so every entry,
        activity and call figure is rebuilt from the converged fixed point rather
        than read off a class row, in the same layout get_ensemble_avg returns.

        Returns:
            (QN, UN, RN, TN, AN, WN), each indexed by absolute element index
        """
        lqn = self.lqn
        n = lqn.nidx
        QN = np.full(n, np.nan)
        UN = np.full(n, np.nan)
        RN = np.full(n, np.nan)
        TN = np.full(n, np.nan)
        AN = np.full(n, np.nan)
        WN = np.full(n, np.nan)
        PN = np.full(n, np.nan)  # processor utilization
        UT = np.full(n, np.nan)  # task and entry utilization

        for a in range(lqn.nacts):
            aidx = lqn.ashift + a
            tidx = int(lqn.parent[aidx, 0])
            if self.ignore[tidx]:
                continue
            hidx = int(lqn.parent[tidx, 0])
            TN[aidx] = self.tput[aidx]
            RN[aidx] = self.servt[aidx]
            UT[aidx] = self.tput[aidx] * self.servt[aidx]
            # LINE scales the utilization of a queueing station into [0,1] whatever its
            # multiplicity, and reports a mean number of busy servers at an infinite
            # server: the processor share of an activity follows the same convention
            PN[aidx] = self.tput[aidx] * self._ph_hostdem_mean(aidx) / self._ph_host_servers(hidx)
            if np.isnan(PN[hidx]):
                PN[hidx] = 0.0
            PN[hidx] += PN[aidx]

        for e in range(lqn.nentries):
            eidx = lqn.eshift + e
            tidx = int(lqn.parent[eidx, 0])
            if self.ignore[tidx]:
                continue
            TN[eidx] = self.tput[eidx]
            RN[eidx] = self.servt[eidx]
            UT[eidx] = self.tput[eidx] * self.servt[eidx]
            acts = self._ph_acts_of(eidx)
            if acts:
                PN[eidx] = float(np.nansum([PN[a] for a in acts]))
            # ResidT is reported per visit to the TASK, not per execution of the
            # activity: an activity of this entry runs EXECS times per invocation, and
            # the entry takes SHARE of the task's invocations. RespT stays per
            # execution.
            if self._ph_execs[eidx] is not None:
                ex = self._ph_execs[eidx]
                w = self._ph_share[eidx]
                for aidx in acts:
                    WN[aidx] = w * ex[aidx] * self.residt[aidx]
            if np.isnan(UT[tidx]):
                UT[tidx] = 0.0
            UT[tidx] += UT[eidx]

        for t in range(lqn.ntasks):
            tidx = lqn.tshift + t
            if self.ignore[tidx]:
                continue
            TN[tidx] = self.tput[tidx]
            acts = self._ph_acts_of(tidx)
            if acts:
                PN[tidx] = float(np.nansum([PN[a] for a in acts]))
                WN[tidx] = float(np.nansum([WN[a] for a in acts]))

        for hidx in range(lqn.nhosts):
            TN[hidx] = np.nan  # kept NaN for consistency with LQNS

        # Idle, not undefined -- the same rule getEnsembleAvg applies, and for the
        # same reason: an unreachable element reports zero for the measures its kind
        # HAS and NaN for the ones it never has, so that the table's NaN mask
        # survives a disconnected component. Reported columns here are QLen=UT,
        # Util=PN, RespT=RN, ResidT=WN, ArvR=AN, Tput=TN; the pre-swap QN and UN are
        # discarded below and are not written.
        for idx in range(n):
            if self.ignore[idx]:
                PN[idx] = 0.0        # every kind reports a utilization
                AN[idx] = np.nan     # nothing reports an arrival rate on an LQN
                kind = self._get_type(idx)
                if kind == LayeredNetworkElement.PROCESSOR:
                    UT[idx] = RN[idx] = WN[idx] = TN[idx] = np.nan
                elif kind == LayeredNetworkElement.TASK:
                    UT[idx] = WN[idx] = TN[idx] = 0.0
                    RN[idx] = np.nan
                elif kind == LayeredNetworkElement.ENTRY:
                    UT[idx] = RN[idx] = TN[idx] = 0.0
                    WN[idx] = np.nan
                elif kind == LayeredNetworkElement.ACTIVITY:
                    UT[idx] = RN[idx] = WN[idx] = TN[idx] = 0.0

        return UT, PN, RN, TN, AN, WN

    # =================================================================
    # Feature gate
    # =================================================================

    def _assert_srvn_ph_supported(self, flat: bool = False):
        """
        Features the composed law cannot represent are refused by name rather
        than silently degraded -- see _kb/06-solver-catalog.md (LN section). The
        list is a property of the ENCODING, so it is the same under either
        layering; what the squashing adds on top is refused in
        _ph_flat_server_set.
        """
        lqn = self.lqn
        mname = 'flat.ph' if flat else 'srvn.ph'
        if getattr(self, 'hasPhase2', False):
            raise ValueError("method='%s' does not support second-phase activities: the "
                             "composed entry law has no reply point. Use method='default'." % mname)
        for cidx in range(lqn.ncalls):
            if self._ph_call_type(cidx) == _FWD:
                raise ValueError("method='%s' does not support forwarding calls, whose target "
                                 "is not part of the caller's activity graph. Use method='default'." % mname)
        iscache = getattr(lqn, 'iscache', None)
        if iscache is not None and np.any(np.asarray(iscache).ravel()):
            raise ValueError("method='%s' does not support cache tasks. Use method='default'." % mname)
        # A SetupTask IS supported: the setup is not part of the activity graph, so it
        # never enters the series-parallel reduction and is prefixed to the composed
        # entry law afterwards as the phase-type mixture. An INF task is the exception,
        # as in LDES: it holds no thread to power down, so the cycle has no meaning.
        hs = getattr(lqn, 'hassetup', None)
        if hs is not None:
            hsf = np.asarray(hs).ravel()
            for i in range(len(hsf)):
                if not hsf[i]:
                    continue
                if self._get_sched(i) == SchedStrategy.INF or not np.isfinite(float(lqn.mult[0, i])):
                    raise ValueError("method='%s': task '%s' declares a setup time on an "
                                     "infinite-server task, which holds no thread to power down; "
                                     "give it a finite multiplicity." % (mname, self._ph_name(i)))
        if getattr(lqn, 'callgroups', None):
            # The group states the ORDER in which one caller visits several
            # callees, and the composed law folds every call into one visit, so
            # the order has nowhere to be expressed. Squashing does not recover
            # it: 'flat.cs' is the only encoding that dispatches a group.
            raise ValueError("method='%s' does not support routed call groups, whose dispatch "
                             "order is a routing property. Use method='flat.cs'." % mname)
        if getattr(lqn, 'lincon', None):
            raise ValueError("method='%s' does not support admission constraints on a layer "
                             "station. Use method='default'." % mname)
        # A queue-dependent service rate is a property of the layer STATION, and
        # the composed law replaces that station by an entry law, so the scaling
        # has nowhere to attach. Only _add_layer_rate_dependence emits it; this
        # encoding used to DROP it in silence, which reads as a solved model
        # rather than a refused one.
        for fndep in ('lldscaling', 'cdscaling', 'jdscaling', 'pools'):
            dep = getattr(lqn, fndep, None) or {}
            if dep:
                sidxdep = sorted(dep.keys())[0]
                what = 'server pools' if fndep == 'pools' else fndep
                raise ValueError("method='%s' does not support queue-dependent service rates on "
                                 "a layer station ('%s' declares %s). Use method='srvn.cs'."
                                 % (mname, self._ph_name(sidxdep), what))

    # =================================================================
    # Small helpers
    # =================================================================

    def _ph_host_layer_callers(self, hidx: int) -> List[int]:
        """Tasks that run on processor HIDX and reach it with requests."""
        out = []
        for tidx in self._get_tasks_of_host(hidx):
            if self.ignore[tidx]:
                continue
            if self._is_ref_task(tidx):
                out.append(tidx)
                continue
            for eidx in self._ph_entries_of(tidx):
                if self._ph_any_caller_of(eidx) or self._ph_has_open_arrival(eidx):
                    out.append(tidx)
                    break
        return out

    def _ph_task_layer_callers(self, tidx: int) -> List[int]:
        """Tasks issuing a synchronous call to an entry of TIDX."""
        lqn = self.lqn
        out = []
        for c in range(lqn.tshift, lqn.tshift + lqn.ntasks):
            if c == tidx or self.ignore[c]:
                continue
            for eidx in self._ph_entries_of(tidx):
                if lqn.issynccaller[c, eidx]:
                    out.append(c)
                    break
        return out

    def _ph_async_calls_into(self, tidx: int) -> List[int]:
        """Asynchronous calls whose target entry belongs to TIDX."""
        lqn = self.lqn
        targets = set(self._ph_entries_of(tidx))
        return [cidx for cidx in range(lqn.ncalls)
                if self._ph_call_type(cidx) == _ASYNC and int(lqn.callpair[cidx, 1]) in targets]

    def _ph_open_arrival_only(self, tidx: int) -> bool:
        """
        True when an entry arrival is the ONLY way requests reach task TIDX.
        'srvn.ph' refuses forwarding calls outright, so sync/async callers are the
        whole test.
        """
        if self._is_ref_task(tidx):
            return False
        for eidx in self._ph_entries_of(tidx):
            if self._ph_any_caller_of(eidx):
                return False
        return any(self._ph_has_open_arrival(eidx) for eidx in self._ph_entries_of(tidx))

    def _ph_any_caller_of(self, eidx: int) -> bool:
        lqn = self.lqn
        return bool(np.any(lqn.issynccaller[:, eidx])) or bool(np.any(lqn.isasynccaller[:, eidx]))

    def _ph_has_open_arrival(self, eidx: int) -> bool:
        arr = getattr(self.lqn, 'arrival', None)
        return isinstance(arr, dict) and arr.get(eidx) is not None

    def _ph_replica_count(self, idx: int, callers: List[int], ishost: bool) -> int:
        """
        Replicas of the server station, with the same fan-out reduction as the
        default builder: a caller that reaches every replica sees one representative.
        """
        lqn = self.lqn
        raw = int(self._ph_repl(idx))
        if raw <= 1 or not callers:
            return max(1, raw)
        reduce = False
        if not ishost and getattr(lqn, 'fanout', None) is not None:
            reduce = all(lqn.fanout[c, idx] >= raw for c in callers)
        elif ishost:
            reduce = all(int(self._ph_repl(c)) == raw for c in callers)
        if reduce:
            if not ishost:
                self.single_replica_tasks.append(idx)
            return 1
        return raw

    def _ph_layer_population(self, idx: int, c: int, nreplicas: int) -> float:
        """Threads of caller C present in the layer of IDX."""
        single = (nreplicas == 1 and self._ph_repl(idx) > 1) or (c in self.single_replica_tasks)
        mc = self._ph_maxmult(c)
        njobs = mc if single else mc * self._ph_repl(c)
        if not np.isfinite(njobs):
            njobs = sum(self._ph_maxmult(i) for i in self._get_callers_of_task(c))
            if not np.isfinite(njobs) or njobs == 0:
                tot = 0.0
                for i in range(self.lqn.nidx):
                    m = self._ph_maxmult(i)
                    if np.isfinite(m):
                        tot += m * self._ph_repl(i)
                njobs = min(tot, 1000.0)
        return float(njobs)

    def _ph_layer_pop(self, L: PHLayer, idx: int) -> float:
        """Closed population of the MODEL the server sits in, i.e. how many jobs a
        job can queue behind. Under 'flat.ph' that is every caller of the single
        network and not only the callers of this one station, which is why it is
        taken from the layer record rather than recomputed from the callers."""
        if getattr(L, 'npop', 0.0) >= 1.0:
            return float(L.npop)
        n = 0.0
        for c in L.callers:
            v = self.njobs[c, idx]
            if np.isfinite(v) and v > 0:
                n += v
        return max(n, 1.0)

    @staticmethod
    def _ph_residence(Q: float, X: float, RN: float) -> float:
        """
        Residence time per visit, by Little from the queue length rather than from
        the reported RN. A layer that saturates can come back from AMVA with an RN
        that no closed model can produce, and a reconstruction that trusts it feeds
        the impossible value straight back into the call response times.
        """
        if np.isfinite(Q) and Q >= 0 and np.isfinite(X) and X > GlobalConstants.FineTol:
            return Q / X
        return RN

    @staticmethod
    def _ph_inflation_of(R: float, S: float, npop: float) -> float:
        """
        Ratio of a residence time to the mean of the law it was measured against,
        bounded above by the layer population: a job can wait behind at most every
        other job in a closed layer.
        """
        f = 1.0
        if S > GlobalConstants.FineTol and np.isfinite(R) and R > 0:
            f = R / S
        if not np.isfinite(f) or f < 1:
            f = 1.0
        if np.isfinite(npop) and npop >= 1 and f > npop:
            f = npop
        return f

    def _ph_sync_calls_between(self, c: int, tidx: int) -> List[int]:
        """Synchronous calls issued by task C to an entry of task TIDX."""
        lqn = self.lqn
        out = []
        for cidx in range(lqn.ncalls):
            if self._ph_call_type(cidx) != _SYNC:
                continue
            if int(lqn.parent[int(lqn.callpair[cidx, 0]), 0]) == c \
                    and int(lqn.parent[int(lqn.callpair[cidx, 1]), 0]) == tidx:
                out.append(cidx)
        return out

    def _ph_setup_prob(self, eidx: int) -> float:
        """
        Probability that a request for entry EIDX finds its task's thread powered
        off. ONE closure for both methods: _setup_charge returns p*s, so p is that
        over self. It also answers p = 1 during construction, before the first solve
        has sized tput or util.
        """
        lqn = self.lqn
        tidx = int(lqn.parent[eidx, 0])
        hs = getattr(lqn, 'hassetup', None)
        if hs is None:
            return 0.0
        hsf = np.asarray(hs).ravel()
        if tidx >= len(hsf) or not hsf[tidx]:
            return 0.0
        d = self._setup_dist_mean(getattr(lqn, 'delayofftime', None), tidx)
        st = self._setup_dist_mean(getattr(lqn, 'setuptime', None), tidx)
        if not (d > GlobalConstants.FineTol) or not (st > GlobalConstants.FineTol):
            return 0.0
        return min(1.0, max(0.0, self._setup_charge(tidx) / st))

    def _ph_setup_law(self, tidx: int):
        """Phase-type law of task TIDX's setup time, None when it declares none."""
        procs = getattr(self.lqn, 'setuptime', None)
        if not isinstance(procs, dict):
            return None
        proc = procs.get(tidx)
        if proc is None:
            return None
        try:
            m = float(proc.getMean())
            scv = float(proc.getSCV())
        except (AttributeError, TypeError, ValueError):
            return None
        if not np.isfinite(m) or m <= GlobalConstants.FineTol:
            return None
        if not np.isfinite(scv) or scv <= GlobalConstants.FineTol:
            scv = 1.0
        law = APH.fitMeanAndSCV(m, scv)
        return _alpha_of(law), _subgen_of(law)

    def _ph_host_servers(self, hidx: int) -> float:
        """
        Divisor that scales a processor utilization into [0,1]. An infinite server
        reports a mean number of busy servers instead, so it divides by one.
        """
        if self._get_sched(hidx) == SchedStrategy.INF:
            return 1.0
        m = self._ph_maxmult(hidx)
        return m if (np.isfinite(m) and m > 0) else 1.0

    @staticmethod
    def _ph_sum_over(M, stations: List[int], col: int) -> float:
        s = 0.0
        for st in stations:
            v = M[st - 1, col]
            if not np.isnan(v):
                s += v
        return s

    def _ph_entries_of(self, idx: int) -> List[int]:
        return list(self.lqn.entriesof.get(idx, []))

    def _ph_acts_of(self, idx: int) -> List[int]:
        return list(self.lqn.actsof.get(idx, []))

    def _ph_name(self, idx: int) -> str:
        v = self.lqn.names
        return v.get(idx, 'Node_%d' % idx) if isinstance(v, dict) else str(v[idx])

    def _ph_call_type(self, cidx: int) -> int:
        ct = getattr(self.lqn, 'calltype', None)
        if ct is None:
            return _SYNC
        arr = np.asarray(ct).ravel()
        return int(arr[cidx]) if cidx < len(arr) else _SYNC

    def _ph_call_mean(self, cidx: int) -> float:
        return float(self.lqn.callpair[cidx, 2])

    def _ph_hostdem_mean(self, aidx: int) -> float:
        hd = self.lqn.hostdem
        if isinstance(hd, dict):
            return float(hd.get(aidx, 0.0))
        return float(np.asarray(hd).ravel()[aidx])

    def _ph_act_think_mean(self, aidx: int) -> float:
        at = getattr(self.lqn, 'actthink', None)
        if not isinstance(at, dict):
            return 0.0
        d = at.get(aidx)
        if d is None:
            return 0.0
        if isinstance(d, (int, float)):
            return float(d) if float(d) > GlobalConstants.FineTol else 0.0
        try:
            m = float(d.getMean())
        except (AttributeError, TypeError, ValueError):
            return 0.0
        return m if (np.isfinite(m) and m > GlobalConstants.FineTol) else 0.0

    def _ph_maxmult(self, idx: int) -> float:
        return float(np.asarray(self.lqn.maxmult).ravel()[idx])

    def _ph_repl(self, idx: int) -> float:
        r = getattr(self.lqn, 'repl', None)
        if r is None:
            return 1.0
        v = float(np.asarray(r).ravel()[idx])
        return v if v > 0 else 1.0

    def _ref_think_mean(self, tidx: int) -> float:
        """Declared think time of a task as it enters the thread cycle: the value
        for a REFERENCE task, zero for any other.

        A think time is an attribute of the closed customer population a
        reference task stands for, and it is what separates one request of that
        population from the next. On a served task it has no such meaning, and
        charging it per request throttles the task: lqn_basic's T3 has 25 threads
        and a declared think time of 4, and reading it as a per-request delay
        caps it at 25/(4+0.02) = 6.219 completions per second. Three independent
        oracles put the rate at five calls per caller request instead -- lqsim
        66.5, LDES 66.955, lqns 75.6. See _kb/06-solver-catalog.md (LN section).
        """
        if not self._is_ref_task(tidx):
            return 0.0
        lqn = self.lqn
        # lqn.think is keyed by ABSOLUTE element index and is a dict whenever the
        # struct was built sparsely, so its length is the number of tasks that
        # declare a think time and not an index bound: a positional guard here
        # read a reference task's think time as absent and dropped it.
        think = getattr(lqn, 'think', None)
        if isinstance(think, dict):
            tp = think.get(tidx)
        elif think is not None and tidx < len(think):
            tp = think[tidx]
        else:
            tp = None
        if tp is None:
            return 0.0
        for attr in ('getMean', 'get_mean'):
            if hasattr(tp, attr):
                v = getattr(tp, attr)()
                return float(v) if np.isfinite(v) and v > 0 else 0.0
        if isinstance(tp, (int, float, np.integer, np.floating)):
            return float(tp) if np.isfinite(tp) and tp > 0 else 0.0
        if hasattr(tp, 'mean'):
            v = tp.mean
            return float(v) if np.isfinite(v) and v > 0 else 0.0
        return 0.0

    def _mwrbb_think_mean(self, tidx: int) -> float:
        # same closure as update_think_times, so the same gate -- see
        # _ref_think_mean
        if not self._is_ref_task(tidx):
            return 0.0
        tp = self.thinkproc[tidx] if (self.thinkproc is not None
                                      and tidx < len(self.thinkproc)) else None
        if tp is None:
            return 0.0
        for attr in ('getMean', 'get_mean'):
            if hasattr(tp, attr):
                v = getattr(tp, attr)()
                return float(v) if np.isfinite(v) else 0.0
        if hasattr(tp, 'mean'):
            return float(tp.mean)
        return 0.0

    @staticmethod
    def _mwrbb_disc_code(s) -> int:
        # 0=FIFO, 1=PS, 2=non-preemptive priority, 3=preemptive priority,
        # 4=ABA full-contention (discipline-independent)
        name = s.name if hasattr(s, 'name') else str(s)
        if name == 'FCFS':
            return 0
        if name in ('PS', 'DPS', 'GPS', 'PSPRIO', 'DPSPRIO', 'GPSPRIO'):
            return 1
        if name in ('HOL', 'FCFSPRIO'):
            return 2
        if name in ('FCFSPRPRIO', 'LCFSPRPRIO'):
            return 3
        return 4   # non-FCFS work-conserving: ABA full-contention

    def _mwrbb_visit_entry(self, eidx, mult, r, D, Vis):
        Vis[eidx, r] += mult
        tidx = self._get_parent(eidx)
        if tidx is not None and tidx < Vis.shape[0]:
            Vis[tidx, r] += mult
        for aidx in self._get_activities_of_entry(eidx):
            if self._get_parent(aidx) == tidx:
                self._mwrbb_visit_activity(aidx, mult, r, D, Vis)

    def _mwrbb_visit_activity(self, aidx, mult, r, D, Vis):
        Vis[aidx, r] += mult
        tidx = self._get_parent(aidx)
        hidx = self._get_parent(tidx) if tidx is not None else None
        dem = self._get_hostdem_mean(aidx)
        if hidx is not None:
            hrow = hidx - self.lqn.hshift        # host rows are 0-based
            if 0 <= hrow < D.shape[0]:
                D[hrow, r] += mult * dem
        calls = self.lqn.callsof.get(aidx, []) if isinstance(self.lqn.callsof, dict) else []
        for cidx in calls:
            if self._get_calltype(cidx) == CallType.SYNC:
                cmean = self._get_call_mean(cidx)
                callee = self._get_call_target_entry(cidx)
                if callee is not None:
                    self._mwrbb_visit_entry(callee, mult * cmean, r, D, Vis)

    def _box_bounds(self, upper: bool):
        from ...api.pfqn.bounds import pfqn_mwrbb
        lqn = self.lqn
        nH = lqn.nhosts
        nidx = lqn.nidx
        refs = [lqn.tshift + t for t in range(lqn.ntasks)
                if self._is_ref_task(lqn.tshift + t)]
        R = len(refs)
        D = np.zeros((nH, R))
        Vis = np.zeros((nidx, R))
        N = np.zeros(R)
        Z = np.zeros(R)
        for r, tidx in enumerate(refs):
            m = self._get_mult(tidx)
            N[r] = 1.0 if not np.isfinite(m) else m
            Z[r] = self._mwrbb_think_mean(tidx)
            for eidx in self._get_entries_of_task(tidx):
                self._mwrbb_visit_entry(eidx, 1.0, r, D, Vis)
        V = (D > 0).astype(float)
        S = D
        sched = np.zeros(nH)
        for h in range(nH):
            sched[h] = self._mwrbb_disc_code(self._get_sched(lqn.hshift + h))
        prio = np.zeros(R)
        Xlo, Xup, _ = pfqn_mwrbb(V, S, N, Z, sched, prio)
        X = Xup if upper else Xlo
        TN = np.full(nidx, np.nan)
        UN = np.full(nidx, np.nan)
        for i in range(nidx):
            if np.any(Vis[i] > 0):
                TN[i] = float(np.sum(X * Vis[i]))
        for h in range(nH):
            UN[lqn.hshift + h] = float(np.sum(X * D[h]))
        return TN, UN

    def _box_bounds_table(self, upper: bool) -> pd.DataFrame:
        TN, UN = self._box_bounds(upper)
        lqn = self.lqn
        rows = []
        for idx in range(lqn.nidx):
            t = TN[idx]
            u = UN[idx]
            rows.append({
                'Node': self._get_hashname(idx),
                'NodeType': self._get_type_name(idx),
                'QLen': 0.0,
                'Util': 0.0 if np.isnan(u) else u,
                'RespT': 0.0,
                'ResidT': 0.0,
                'ArvR': 0.0,
                'Tput': 0.0 if np.isnan(t) else t,
            })
        df = pd.DataFrame(rows)
        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(df)

    @staticmethod
    def _sanitize_avg(values):
        """Snap near-tenth entries onto the tenth and flatten near-zero entries.

        Mirrors the sanitization MATLAB SolverLN applies before formatting
        getAvgTable. NaN entries are left untouched, as are negative ones (the
        MATLAB relative test is vacuous there).
        """
        v = np.asarray(values, dtype=float).copy()
        with np.errstate(invalid='ignore'):
            scaled = v * 10.0
            rounded = np.round(scaled)
            to_round = np.abs(scaled - rounded) < GlobalConstants.CoarseTol * scaled
            to_round &= ~np.isnan(v)
            v[to_round] = rounded[to_round] / 10.0
            v[(~np.isnan(v)) & (v <= GlobalConstants.FineTol)] = 0.0
        return v

    def get_avg_table(self) -> pd.DataFrame:
        """Get average metrics as a table (matches MATLAB getAvgTable)."""
        # lang=cpp is tested BEFORE the mwba branch below: those bounds are
        # computed natively, so serving them here would report python numbers
        # under a C++ label. The dispatch refuses the method by name instead.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, ln_avg_table_via_cpp
            try:
                df = ln_avg_table_via_cpp(self)
            except LineCliNotAvailable as e:
                line_warning("SolverLN", "lang='cpp' requested but the C++ solver is "
                             "unavailable (%s); falling back to lang='python'." % e)
            else:
                if not self._table_silent and len(df) > 0:
                    print(df.to_string(index=False))
                from line_solver.indexed_table import IndexedTable
                return IndexedTable(df)

        # Majumdar-Woodside robust box bounds for the LQN (processor-contention)
        _bmethod = getattr(self.options, 'method', None)
        if _bmethod in ('mwba.upper', 'mwba.lower'):
            return self._box_bounds_table(_bmethod == 'mwba.upper')

        # print the lang=java ensemble average table unless output is silenced.
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import ln_avg_table_via_jar
            df = ln_avg_table_via_jar(self)
            if not self._table_silent and len(df) > 0:
                print(df.to_string(index=False))
            from line_solver.indexed_table import IndexedTable
            return IndexedTable(df)

        QN, UN, RN, TN, AN, WN = self.get_ensemble_avg()

        # sanitize averages (snap-to-tenth/flatten-to-zero) before formatting; see _kb/11-conventions-and-gotchas.md LN avg tables are sanitized.
        QN, UN, RN, TN, AN, WN = (self._sanitize_avg(v) for v in (QN, UN, RN, TN, AN, WN))

        lqn = self.lqn

        # Build table
        rows = []
        for idx in range(lqn.nidx):
            name = self._get_hashname(idx)
            node_type = self._get_type_name(idx)

            # Get metric values
            qlen = QN[idx]
            util = UN[idx]
            respt = RN[idx]
            residt = WN[idx]
            arvr = AN[idx]
            tput = TN[idx]

            rows.append({
                'Node': name,
                'NodeType': node_type,
                'QLen': qlen,
                'Util': util,
                'RespT': respt,
                'ResidT': residt,
                'ArvR': arvr,
                'Tput': tput,
            })

        df = pd.DataFrame(rows)

        # Print table if not silent (matches LQNS behavior)
        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        # Return as IndexedTable for MATLAB-style number formatting
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(df)

    def _get_type_name(self, idx: int) -> str:
        """Get element type name."""
        lqn = self.lqn
        elem_type = self._get_type(idx)

        if elem_type == LayeredNetworkElement.PROCESSOR:
            return 'Processor'
        elif elem_type == LayeredNetworkElement.TASK:
            if self._is_ref_task(idx):
                return 'RefTask'
            return 'Task'
        elif elem_type == LayeredNetworkElement.ENTRY:
            return 'Entry'
        elif elem_type == LayeredNetworkElement.ACTIVITY:
            return 'Activity'
        elif elem_type == LayeredNetworkElement.CALL:
            return 'Call'
        return 'Unknown'

    def reset(self):
        """Reset solver state."""
        self.hasconverged = False
        self.moment_pass_done = False
        self.results = []
        self.maxitererr = []

    # per-layer state fields carried across iterations (util/tput/servt/residt/thinkt/callresidt/callservt and their _prev counterparts).
    _STATE_ATTRS = (
        'servtproc', 'thinktproc', 'callservtproc', 'tputproc', 'entryproc',
        'util', 'tput', 'servt', 'residt', 'thinkt', 'callresidt', 'callservt',
        'relax_omega', 'servt_prev', 'residt_prev', 'tput_prev', 'thinkt_prev',
        'callservt_prev', 'callresidt_prev', 'results', 'njobs', 'ilscaling',
    )

    def get_state(self):
        """Export current solver state for continuation with another solver.

        Returns a dict snapshotting service/think/call processes, performance
        metrics, relaxation state and last-iteration results, consumable by
        set_state() on a solver built with a different layer-solver factory.
        """
        import copy
        state = {}
        for attr in self._STATE_ATTRS:
            if hasattr(self, attr):
                state[attr] = copy.deepcopy(getattr(self, attr))
        return state

    def set_state(self, state):
        """Import a previously exported state (see get_state) for continuation."""
        import copy
        for attr, value in state.items():
            setattr(self, attr, copy.deepcopy(value))

    def update_solver(self, solver_factory):
        """Replace all per-layer solvers with ones built by solver_factory,
        preserving the current solution state for refinement."""
        self.solver_factory = solver_factory
        for e in range(self.nlayers):
            solver = solver_factory(self.ensemble[e])
            self._assert_layer_solver_supports_model(solver, self.ensemble[e], e)
            self._detach_layer_config(solver)
            self._silence_layer_solver(solver)
            self.solvers[e] = solver

    def _layer_options(self):
        """The LN options as a LAYER solver may read them.

        An LN method name states the LAYERING and the ENCODING of the ensemble,
        not the algorithm a single layer is solved with, and the two vocabularies
        do not overlap. The default factory hands the LN options straight to
        SolverMVA, so naming any LN method -- 'srvn.ph', 'flat.cs', 'flat.ph',
        even the 'srvn' alias -- reached the layer solver as its own method. It
        refused the unknown token outright, or, worse, resolved it to something
        else and answered: LN(model, method='srvn') returned Tput 0.694282 on the
        two-task probe where the identical LN(model) returns 1.402605, both having
        resolved lnmethod to 'srvn.ph'. MATLAB keeps the two apart by giving each
        layer an options struct of its own.
        """
        import copy as _copy
        opts = self.options
        m = getattr(opts, 'method', None)
        if isinstance(m, str) and m.lower() in _LN_LEVEL_METHODS:
            opts = _copy.copy(opts)
            opts.method = 'default'
        return opts

    def _detach_layer_config(self, layer_solver):
        """Give the layer solver a config dict of its own.

        MATLAB hands each layer an options STRUCT, copied by value; in Python every layer
        solver built from the LN options aliases one config dict, so a per-layer write such
        as the interlock matrix of Eq. (4.7) would reach every other layer, whose classes are
        neither the same in number nor in meaning.
        """
        opts = getattr(layer_solver, 'options', None)
        cfg = getattr(opts, 'config', None) if opts is not None else None
        if isinstance(cfg, dict):
            opts.config = type(cfg)(cfg)

    def _silence_layer_solver(self, layer_solver):
        """Set a layer solver to VerboseLevel.SILENT (spelled False here).

        A LAYER SOLVER NEVER NARRATES. The fixed point runs every layer once per
        iteration, so a layer left at the caller's verbosity prints its own
        banner nlayers*iter_max times and buries the layered narration the caller
        actually asked for. The level is stamped HERE rather than in the factory
        because a factory the USER supplied -- SolverLN(model, lambda m:
        SolverNC(m)) -- never sees the LN options at all, and stamping it in the
        default factory alone left exactly that case loud.

        SolverLN's own reporting is unaffected: it reads self.options.verbose,
        not the layer's.

        A NATIVE SOLVER OPTION SPELLS ITS VERBOSITY AS A BOOL, not as a
        VerboseLevel (see constants.default_verbose), and an Enum member is
        truthy: assigning VerboseLevel.SILENT here would read as VERBOSE at
        every `if options.verbose:` in the tree. False is the same level in the
        type this field actually carries, and is what console._is_silent reads.
        """
        opts = getattr(layer_solver, 'options', None)
        if opts is not None and hasattr(opts, 'verbose'):
            opts.verbose = False

    def _assert_layer_solver_supports_model(self, layer_solver, layer_model, idx):
        """Reject a layer solver that cannot represent its layer model, upfront
        with a clear message.

        LN server-layer stations carry immediate feedback (sn.immfeed):
        successive same-host activities retain the server, modelled as
        immediate-feedback self-loops so the layer solver does not re-queue the
        job. SolverJMT rejects any model with immfeed and returns no solution,
        so a pure-JMT layer factory otherwise fails cryptically at the first
        iteration. Detect it here instead.

        The guard is CONDITIONAL: it fires only when the specific layer model
        actually carries immfeed. A SolverJMT layer solver on an immfeed-free
        layer is allowed, and non-JMT factories are never rejected.
        """
        from ..wrappers.solver_jmt import SolverJMT
        if isinstance(layer_solver, SolverJMT):
            sn = layer_model.get_struct()
            immfeed = getattr(sn, 'immfeed', None)
            if immfeed is not None and np.any(immfeed):
                raise ValueError(
                    "SolverJMT cannot solve LN layer %d: the layer carries "
                    "immediate feedback (sn.immfeed), which SolverJMT does not "
                    "support, so LN would fail at the first iteration. Use the "
                    "default layer factory (MVA/NC) or another layer solver "
                    "that supports immediate feedback." % idx)

    # camelCase aliases (JAR-style)
    def getState(self):
        return self.get_state()

    def setState(self, state):
        return self.set_state(state)

    def updateSolver(self, solver_factory):
        return self.update_solver(solver_factory)

    @staticmethod
    def defaultOptions() -> SolverLNOptions:
        """Get default LN solver options."""
        return SolverLNOptions()

    @staticmethod
    def default_options() -> SolverLNOptions:
        """Get default options (Python convention)."""
        return SolverLNOptions()

    def getSensitivityTable(self, method='auto', step=None, scheme='forward'):
        """Layer-wise performance sensitivities of a layered network.

        Solves the layered model and then delegates to each layer solver,
        returning the concatenation of the layer tables with a leading Layer
        column. Every row is a (Layer, Station, JobClass) triple carrying the
        derivative of that row's mean measures with respect to that
        station-class service RATE: dTput_dRate, dRespT_dRate, dQLen_dRate,
        dUtil_dRate.

        The options are passed through to the layer solvers unchanged, with the
        same meaning as in NetworkSolver.getSensitivityTable: each layer
        independently takes the analytic branch where its own solver supports it
        and the model is in scope, and finite differences otherwise. In practice
        a layer submodel is chain-based (the callers switch class), which puts
        it out of scope of the analytic branch, so the layers normally
        finite-difference their own solver.

        IMPORTANT, on what these derivatives mean. Each entry is a derivative
        WITHIN ITS LAYER, taken with the layer parameters that the fixed point
        produced held fixed. It is a partial derivative of the layer submodel,
        not the total derivative of the layered model: perturbing a host demand
        in one layer moves the think times, populations and service rates of the
        other layers through the fixed-point map, and that indirect term is not
        included here. The layer table is the right object for attributing a
        bottleneck inside a layer, and the wrong one for predicting the effect of
        a parameter change on the solved layered model.

        Returns the table; the per-layer second outputs are on
        ``DataFrame.attrs['sens']`` and the per-layer branch labels on
        ``DataFrame.attrs['layer_methods']``, with ``attrs['method']`` the
        summary ('exact', 'fd', or 'mixed').
        """
        import pandas as pd

        # The layered CLI has no sensitivity analysis, so the derivatives below are
        # native; refuse rather than label them as C++ numbers.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            raise RuntimeError(
                "lang='cpp' delegates the steady-state layered solve only; the per-layer "
                "sensitivities are computed natively. Use lang='python' for "
                "getSensitivityTable.")

        # gate accessors on self.results, not the ensemble table; see _kb/11-conventions-and-gotchas.md lang=java SolverLN skips iterate().
        if not self.results:
            self.iterate()

        rows = []
        sens = []
        layer_methods = []
        for e, solver in enumerate(self.solvers):
            if solver is None:
                sens.append(None)
                layer_methods.append(None)
                continue
            T = solver.getSensitivityTable(method=method, step=step, scheme=scheme)
            sens.append(T.attrs.get('sens'))
            layer_methods.append(T.attrs.get('method'))
            layer_name = self.ensemble[e].getName()
            for _, row in T.iterrows():
                rows.append({
                    'Layer': str(layer_name),
                    'Station': row['Station'],
                    'JobClass': row['JobClass'],
                    'dTput_dRate': float(row['dTput_dRate']),
                    'dRespT_dRate': float(row['dRespT_dRate']),
                    'dQLen_dRate': float(row['dQLen_dRate']),
                    'dUtil_dRate': float(row['dUtil_dRate']),
                })

        table = pd.DataFrame(rows, columns=['Layer', 'Station', 'JobClass',
                                            'dTput_dRate', 'dRespT_dRate',
                                            'dQLen_dRate', 'dUtil_dRate'])
        present = [m for m in layer_methods if m]
        if not present:
            summary = ''
        elif all(m == present[0] for m in present):
            summary = present[0]
        else:
            summary = 'mixed'
        table.attrs['method'] = summary
        table.attrs['layer_methods'] = layer_methods
        table.attrs['sens'] = sens
        return table

    get_sensitivity_table = getSensitivityTable

    # Aliases for compatibility
    avg_table = get_avg_table
    getAvgTable = get_avg_table
    avgTable = get_avg_table
    avgT = get_avg_table
    aT = get_avg_table


# Alias for compatibility
LN = SolverLN
