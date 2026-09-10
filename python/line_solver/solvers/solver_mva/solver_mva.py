"""
Native Python implementation of Mean Value Analysis (MVA) solver.

This implementation uses pure Python/NumPy algorithms from the api.pfqn
module.
"""

import os
import numpy as np
import pandas as pd
import sys
import time
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass, field
from ...constants import default_verbose
from enum import Enum

from ...api.sn.transforms import sn_get_residt_from_respt, get_chain_for_class
from ...api.sn.network_struct import NodeType
from ...api.io.logging import line_debug, line_warning
from ..base import NetworkSolver, method_type
from ..fork_join_driver import ForkJoinDriverMixin
from ..transform_driver import TransformSolveMixin
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


def _amva_needs_amvald(sn):
    """True when the AMVA linearizer family must go through solver_amvald.

    Mirrors the three MATLAB gates in solver_amva.m: the pfqn_* linearizer path
    requires cond1 = sn_has_product_form_not_het_fcfs(sn) (:137) and
    cond2 = ~sn_has_load_dependence(sn) (:91), and the lin arm itself
    re-checks isempty(sn.cdscaling) (:222). The pfqn_linearizer family carries no
    load-dependence argument, so an LD/CD model routed there loses the scaling
    silently. Presence of the handle is the test, matching MATLAB's
    size(sn.lldscaling,2)>0 / ~isempty(sn.cdscaling) -- a flat lldscaling is
    still routed to solver_amvald there, and pfqn_lldfun no-ops on it anyway.

    cond1 is the het-FCFS exclusion: an FCFS station whose per-class service
    means differ is not BCMP type 1, so MATLAB never reaches its lin arm for one
    and takes the non-product-form tail (:397-401) to solver_amvald instead. Only
    the lin family consults this predicate here, so the ab / schmidt / schmidt-ext
    bypass that api solver_amva carries does not apply.
    """
    from ...api.sn.predicates import sn_has_load_dependence
    from ...api.sn import sn_has_product_form_not_het_fcfs
    if not sn_has_product_form_not_het_fcfs(sn):
        return True
    if sn_has_load_dependence(sn):
        return True
    for attr in ('cdscaling', 'jdscaling'):
        scaling = getattr(sn, attr, None)
        if scaling is None:
            continue
        try:
            if len(scaling) > 0:
                return True
        except TypeError:
            return True
    return False


def _qsys_queue_visits(sn, queue_ist, chain=0, class_idx=0):
    """Visit ratio of the queue station of a single-class open queueing system.

    ``sn.visits`` is STATEFUL-indexed in every codebase, so the station index has
    to go through ``stationToStateful`` first; mirrors
    ``sn.visits{1}(sn.stationToStateful(queue_ist))`` in
    solver_mva_qsys_analyzer.m. Returns 1.0 when the struct carries no visits.
    """
    visits = getattr(sn, 'visits', None)
    if not visits or chain not in visits or visits[chain] is None:
        return 1.0
    V = np.asarray(visits[chain], dtype=float)
    isf = queue_ist
    s2sf = getattr(sn, 'stationToStateful', None)
    if s2sf is not None and len(np.asarray(s2sf).flatten()) > queue_ist:
        isf = int(np.asarray(s2sf).flatten()[queue_ist])
    if V.ndim == 1:
        return float(V[isf]) if isf < V.shape[0] else 1.0
    if isf >= V.shape[0] or class_idx >= V.shape[1]:
        return 1.0
    return float(V[isf, class_idx])


def _bmap_batch_moments(proc):
    """Batch-arrival moments from a BMAP proc entry.

    The entry is laid out as [D0, D1, D_batch1, ..., D_batchK] (JAR MatrixCell
    layout, D1 = sum of batch matrices). Rates are weighted by the stationary
    vector of D0 + sum(Dk), matching MATLAB's BMAP.getBatchRates().

    Returns (lambda_batch, E_X, E_X2): batch event rate, mean and second
    moment of the batch size.
    """
    from ...distributions.markovian import BMAP as _BMAP
    bmap = _BMAP([proc[0]] + list(proc[2:]))
    batch_rates = np.asarray(bmap.getBatchRates(), dtype=float).flatten()
    total_rate = float(np.sum(batch_rates))
    E_X = float(bmap.getMeanBatchSize())
    if total_rate > 0:
        sizes = np.arange(1, len(batch_rates) + 1, dtype=float)
        E_X2 = float(np.sum(sizes ** 2 * batch_rates) / total_rate)
    else:
        E_X2 = E_X ** 2
    return total_rate, E_X, E_X2


def _ph_from_proc(sn, station_idx, class_idx=0):
    """Extract a PH (alpha, T) representation from sn.proc at the given station.

    Handles three storage forms LINE uses for PH/MAP processes:
      - dict {'k', 'mu'}    : Erlang(k, mu) -> bidiagonal sub-generator
      - dict {'rate'}        : Exponential(rate)
      - (D0, D1) tuple/list  : MAP — alpha = map_pie(D0, D1), T = D0
    Returns (alpha_row, T) numpy arrays, or (None, None) on failure.
    """
    try:
        if not hasattr(sn, 'proc') or sn.proc is None:
            return None, None
        proc_st = sn.proc[station_idx] if station_idx < len(sn.proc) else None
        if proc_st is None:
            return None, None
        ph = proc_st[class_idx] if class_idx < len(proc_st) else None
        if ph is None:
            return None, None
        if isinstance(ph, dict):
            if 'k' in ph and 'mu' in ph:
                k_phases = int(ph['k'])
                mu_phase = float(ph['mu'])
                alpha = np.zeros(k_phases); alpha[0] = 1.0
                T = np.zeros((k_phases, k_phases))
                for i in range(k_phases):
                    T[i, i] = -mu_phase
                    if i < k_phases - 1:
                        T[i, i + 1] = mu_phase
                return alpha, T
            if 'rate' in ph:
                r = float(ph['rate'])
                return np.array([1.0]), np.array([[-r]])
            if 'probs' in ph and 'rates' in ph:
                # HyperExp: PH with alpha=probs, sub-generator T=diag(-rates)
                probs = np.asarray(ph['probs'], dtype=float).flatten()
                rates = np.asarray(ph['rates'], dtype=float).flatten()
                return probs, np.diag(-rates)
            return None, None
        if isinstance(ph, (list, tuple)) and len(ph) >= 2:
            a0 = np.asarray(ph[0], dtype=float)
            a1 = np.asarray(ph[1], dtype=float)
            # (alpha, T) PH pair: 1D initial vector + 2D square sub-generator
            # (Coxian/APH arrivals are stored this way, not as a (D0,D1) MAP).
            if a0.ndim == 1 and a1.ndim == 2 and a1.shape[0] == a1.shape[1] and a1.shape[0] == a0.shape[0]:
                return a0, a1
            # (D0, D1) MAP pair: both 2D square of equal shape
            if a0.ndim == 2 and a0.shape == a1.shape and a0.shape[0] == a0.shape[1]:
                from ...api.mam.map_analysis import map_pie as _map_pie
                pie = np.asarray(_map_pie(a0, a1)).flatten()
                return pie, a0
        return None, None
    except Exception:
        return None, None


def _gm1_lst_sojourn(sn, station_idx, mu, class_idx=0):
    """Exact G/M/1 mean sojourn time from the arrival LST stored in sn.lst.

    Solves sigma = LST(mu*(1-sigma)) and returns W = 1/(mu*(1-sigma)). Works for
    ANY renewal arrival whose LST is available (Det, Uniform, Pareto, Lognormal,
    Weibull, ...), and is the exact answer where the PH/M/1 path would only see
    an Erlang moment fit. Returns None when no LST is stored or no caudal root is
    found. Matches MATLAB fzero(@(x) LA(mu-mu*x)-x, 0.5).
    """
    lst_all = getattr(sn, 'lst', None)
    if (station_idx is None or lst_all is None or station_idx >= len(lst_all)
            or lst_all[station_idx] is None or len(lst_all[station_idx]) <= class_idx
            or lst_all[station_idx][class_idx] is None):
        return None
    lst = lst_all[station_idx][class_idx]
    try:
        from scipy.optimize import brentq
        f = lambda x: lst(mu * (1.0 - x)) - x
        # caudal root: smallest root in (0,1); scan FIRST sign change not bracketing (f may be positive near 1); see _kb/06-solver-catalog.md MVA gm1 branch.
        grid = np.linspace(1e-9, 1.0 - 1e-6, 400)
        fv = [f(x) for x in grid]
        for gi in range(len(grid) - 1):
            if fv[gi] * fv[gi + 1] < 0.0:
                sigma = brentq(f, grid[gi], grid[gi + 1])
                return 1.0 / (mu * (1.0 - sigma))
    except Exception:
        return None
    return None


@dataclass
class SolverMVAOptions:
    """Options for the native MVA solver."""
    method: str = 'default'  # auto-selects exact/amva (matches MATLAB/JAR)
    max_iter: int = 1000
    tol: float = 1e-4
    # 1e-6 is the MVA-SPECIFIC iter_tol (MATLAB/JAR SolverOptions 'MVA' case), not the 1e-4 general default; do not correct it.
    iter_tol: float = 1e-6
    verbose: bool = field(default_factory=default_verbose)
    seed: Optional[int] = None  # Random seed (for compatibility, not used in MVA)
    keep: bool = False  # Keep intermediate data (for compatibility)
    cutoff: Optional[int] = None  # State space cutoff (for compatibility, not used in MVA)
    samples: Optional[int] = None  # Samples (for compatibility, not used in MVA)
    fork_join: str = 'default'  # Fork-join method: 'default'/'ht' (H-T), 'mmt' (experimental)
    fj_warmstart: bool = True  # Resume fork-join (MMT) fixed point from iterate kept by previous runAnalyzer call, not FineTol; effective only under outer iteration like SolverLN
    cd_peak_norm: bool = False  # Scale class-dependent station Util by the lattice peak (bmax); default reports unscaled T*S with a warning
    init_sol: Optional[np.ndarray] = None  # Warm-start chain-level queue lengths (M x nchains)
    config: Optional[dict] = None  # Config dict, e.g. {'multiserver': 'softmin'}
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))  # env LINE_SOLVER_LANG overrides; 'python' (native), 'java' (jline.jar via JSON) or 'cpp' (line-cli via JSON)
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. It has no meaning for the other langs -- MATLAB, the JAR
    # and native Python are IEEE double throughout -- so it is left None and the
    # C++ CLI is invoked without --arith unless the caller sets it. An exact
    # solve returns the same doubles here: the wire format carries the double
    # alongside num/den, and only the CLI's own -o json --api path exposes the
    # fraction. What it buys through this option is a solve with no rounding in
    # the middle of it.
    arith: Optional[str] = None


class SolverMVA(TransformSolveMixin, ForkJoinDriverMixin, NetworkSolver):
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

    Bound methods (aba, bjb, pb, gb, sb, mwba, ...) are served by SolverBA,
    not here; runAnalyzer rejects the whole family.

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default', which auto-selects exact/amva based on model)
        **kwargs: Additional solver options
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._mmt_cache = None  # Cache for MMT transformation across LN iterations
        # MMT fork-join arrival rates kept across runAnalyzer calls so outer iteration resumes; see _kb/11-conventions-and-gotchas.md Caching a derived model.
        self._fj_fork_lambda = None

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            # 'method' read from kwargs, defaulting to 'default' (auto-selects amva for non-product-form), matching MATLAB.
            self.method = kwargs.get('method', 'default')
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
            # config carries per-method payloads (multiserver rule, QRF blocking params); dropping it makes a method REQUIRING them unreachable via options object
            if method_or_options.get('config') is not None:
                kwargs.setdefault('config', method_or_options.get('config'))
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object (e.g., from SolverLN)
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'max_iter'):
                kwargs.setdefault('max_iter', method_or_options.max_iter)
            if hasattr(method_or_options, 'iter_max'):  # LN uses iter_max
                kwargs.setdefault('max_iter', method_or_options.iter_max)
            # SolverLN's iter_tol is OUTER LQN fixed-point tol, not inner per-layer MVA's; MATLAB forces only iter_max on layer solver, iter_tol at MVA default.
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if getattr(method_or_options, 'config', None) is not None:
                kwargs.setdefault('config', method_or_options.config)
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

    def supportsExactSensitivity(self):
        """MVA differentiates its own recursion: getSensitivityTable uses the
        analytic branch (pfqn_sens) wherever the model is in scope.
        """
        return True

    supports_exact_sensitivity = supportsExactSensitivity

    get_name = getName

    def reset(self):
        """Reset the solver to force recomputation on next getAvg call."""
        self._clearResultStores()
        self._sn = None
        # _mmt_cache/_fj_fork_lambda NOT cleared on reset(); see _kb/11-conventions-and-gotchas.md Caching a derived model (provenance, not don't-touch).
        self._extract_network_params()

    def _extract_network_params(self):
        """Extract parameters from the model for MVA computation."""
        model = self.model

        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            # If rates are dirty, refresh them before extracting
            if getattr(model, '_rates_dirty', False):
                model.refresh_rates()
            self._from_network_struct(model._sn)
            return

        # Priority 2: Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                self._from_network_struct(model._sn)
                return

        # native models expose snake_case get_struct() directly; no JAR-wrapper bridge, keeping python/ JVM-free.
        if hasattr(model, 'get_struct'):
            sn = model.get_struct()
            if sn is not None:
                self._from_network_struct(sn)
                return

        # Priority 4: Direct model extraction
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
            # sn.visits is indexed by stateful node, demands by station; stationToStateful converts between them.
            visits_combined = np.zeros_like(self.demands)
            stationToStateful = np.asarray(sn.stationToStateful, dtype=int).flatten()
            for chain_id, visits_chain in sn.visits.items():
                if isinstance(visits_chain, np.ndarray):
                    # Convert from stateful visits to station visits
                    for station_idx in range(self.nstations):
                        stateful_idx = stationToStateful[station_idx]
                        if stateful_idx < visits_chain.shape[0]:
                            visits_combined[station_idx, :] += visits_chain[stateful_idx, :]
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

    def _has_sjn_station(self):
        """True when any station schedules by non-preemptive shortest job next.

        Read off the struct rather than the nodes, so it matches what the
        dispatcher tests. Compared by NAME: `sn.sched` carries members of
        `lang.base.SchedStrategy`, and `==` against a member of another live
        SchedStrategy class silently returns False (see
        _kb/11-conventions-and-gotchas.md).
        """
        sn = getattr(self, '_sn', None)
        if sn is None:
            try:
                sn = self.model.getStruct()
            except Exception:
                return False
        sched = getattr(sn, 'sched', None)
        if sched is None:
            return False
        values = sched.values() if hasattr(sched, 'values') else sched
        for sv in values:
            if sv is not None and getattr(sv, 'name', None) == 'SJF':
                return True
        return False

    def _has_prs_prio_station(self):
        """True when any station schedules by preemptive-resume priority.

        Compared BY NAME for the reason _has_sjn_station gives: `sn.sched`
        carries members of `lang.base.SchedStrategy`, and `==` against a member
        of another live SchedStrategy class silently returns False.
        """
        sn = getattr(self, '_sn', None)
        if sn is None:
            try:
                sn = self.model.getStruct()
            except Exception:
                return False
        sched = getattr(sn, 'sched', None)
        if sched is None:
            return False
        values = sched.values() if hasattr(sched, 'values') else sched
        for sv in values:
            if sv is not None and getattr(sv, 'name', None) == 'FCFSPRPRIO':
                return True
        return False

    def _determine_network_type(self):
        """Determine if network is open, closed, or mixed.

        A closed class holding no jobs is not a closed part of the network: the
        test is the POPULATION, as MATLAB's mvaDispatch branches on
        ``sn.nclosedjobs == 0`` rather than on the presence of closed classes.
        Calling such a model mixed sends it to the chain-level load-dependent
        route, whose rate matrix is sized by the closed population, so
        pfqn_ldmx_ec indexes an empty axis and every open-class metric comes
        back zero (a JMT import declaring an empty closed class, as
        test/testsOpenQN/oqn-11.jsimg does, reported throughput at the Source
        and nothing downstream).
        """
        has_open = False
        has_closed = False

        for c in range(self.nclasses):
            # Open classes have njobs = inf, closed classes have finite njobs >= 0
            if np.isinf(self.njobs[c]):
                has_open = True
            elif self.njobs[c] > 0:
                has_closed = True

        if has_open and has_closed:
            self.network_type = 'mixed'
        elif has_closed:
            self.network_type = 'closed'
        else:
            self.network_type = 'open'

    def _get_think_times(self) -> np.ndarray:
        """Extract think times from delay stations and INF-scheduled stations.

        Think time Z[r] is the sum of demands at all Delay stations and
        stations with SchedStrategy.INF for class r.

        Matches MATLAB solver_mva.m lines 51-66 which separates:
        - infSET: stations with SchedStrategy.INF (treated as delays)
        - qSET: other product-form stations (treated as queues)
        """
        from ...api.sn.network_struct import NodeType
        from ...lang.base import SchedStrategy

        Z = np.zeros(self.nclasses)

        # Sum demands from all Delay stations AND INF-scheduled stations
        for i in range(self.nstations):
            is_delay_or_inf = False

            # Check if station type is DELAY
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.DELAY.value:
                        is_delay_or_inf = True

            # Also check if scheduling strategy is INF (infinite servers)
            # MATLAB: case SchedStrategy.INF -> infSET (treated as delay)
            if not is_delay_or_inf and self.sched is not None:
                sched_val = self.sched.get(i) if isinstance(self.sched, dict) else (
                    self.sched[i] if i < len(self.sched) else None)
                if sched_val is not None:
                    is_inf = (sched_val == SchedStrategy.INF or
                             (hasattr(sched_val, 'value') and sched_val.value == SchedStrategy.INF.value) or
                             (isinstance(sched_val, int) and sched_val == SchedStrategy.INF.value))
                    if is_inf:
                        is_delay_or_inf = True

            if is_delay_or_inf:
                for c in range(self.nclasses):
                    if self.demands[i, c] > 0:
                        Z[c] += self.demands[i, c]

        return Z

    def _get_queueing_demands(self) -> Tuple[np.ndarray, List[int]]:
        """
        Get demand matrix for queueing stations only (excluding sources, sinks, delays, and INF).

        Delay stations and INF-scheduled stations contribute to think times,
        not queueing demands.

        Matches MATLAB solver_mva.m lines 51-66 which separates:
        - infSET: stations with SchedStrategy.INF (treated as delays)
        - qSET: other product-form stations (treated as queues)

        Returns:
            Tuple of (demand matrix, list of queueing station indices)
        """
        from ...api.sn.network_struct import NodeType
        from ...lang.base import SchedStrategy

        # Get queueing stations (exclude Source, Sink, Delay, Fork, Join, and INF)
        queue_indices = []
        for i in range(self.nstations):
            # Check if this station should be excluded from queueing demands
            is_excluded = False

            # Check node type
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    # Handle both enum objects and integer values
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    # Exclude Source, Sink, Delay, Fork, and Join stations
                    # (Fork/Join don't do queueing work, they handle synchronization)
                    if st_val in (NodeType.SOURCE.value, NodeType.SINK.value, NodeType.DELAY.value,
                                  NodeType.FORK.value, NodeType.JOIN.value):
                        is_excluded = True

            # Also exclude stations with SchedStrategy.INF (infinite servers)
            # MATLAB: case SchedStrategy.INF -> infSET (NOT in qSET)
            if not is_excluded and self.sched is not None:
                sched_val = self.sched.get(i) if isinstance(self.sched, dict) else (
                    self.sched[i] if i < len(self.sched) else None)
                if sched_val is not None:
                    is_inf = (sched_val == SchedStrategy.INF or
                             (hasattr(sched_val, 'value') and sched_val.value == SchedStrategy.INF.value) or
                             (isinstance(sched_val, int) and sched_val == SchedStrategy.INF.value))
                    if is_inf:
                        is_excluded = True

            # Include station if it has non-zero demand and is not excluded
            if not is_excluded and np.any(self.demands[i, :] > 0):
                queue_indices.append(i)

        if not queue_indices:
            # No queueing stations - return empty demands and empty indices
            # This happens when all stations are delays or INF (infinite server)
            return np.zeros((0, self.nclasses)), []

        L = self.demands[queue_indices, :].copy()

        # Apply DPS weight scaling: for DPS stations, effective demand = D / weight
        from ...lang.base import SchedStrategy
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
        from ...api.sn.network_struct import NodeType

        source_indices = []
        for i in range(self.nstations):
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.SOURCE:
                        source_indices.append(i)
        return source_indices

    def _hol_cobham_applicable(self, q_ist: int, src_ist: int, mi) -> bool:
        """True if a single open queue admits the exact non-preemptive Cobham
        priority formula (qsys_mg1_prio): HOL scheduling, a single server, and
        Poisson (ca=1) arrivals for every class. Compared by scheduling NAME to
        stay robust across the codebases' enum encodings."""
        sched = self._sn.sched.get(q_ist) if hasattr(self._sn.sched, 'get') else None
        sname = sched.name if hasattr(sched, 'name') else str(sched)
        if sname != 'HOL':
            return False
        if mi is not None and len(mi) > 0 and np.isfinite(mi[0]) and int(mi[0]) != 1:
            return False
        scv = self._sn.scv
        if scv is not None:
            for r in range(self.nclasses):
                v = scv[src_ist, r]
                if np.isfinite(v) and abs(float(v) - 1.0) > 1e-6:
                    return False
        return True

    def _dps_exact_applicable(self, q_ist: int, src_ist: int, mi) -> bool:
        """True if a single open queue admits the numerically-exact M/M/1-DPS
        solver (qsys_mm1_dps): DPS scheduling, a single server, Poisson (ca=1)
        arrivals AND exponential (cs=1) service per class, and a small class
        count (the truncated-CTMC state space grows as cutoff^K)."""
        sched = self._sn.sched.get(q_ist) if hasattr(self._sn.sched, 'get') else None
        sname = sched.name if hasattr(sched, 'name') else str(sched)
        if sname != 'DPS':
            return False
        if self.nclasses > 3:
            return False
        if mi is not None and len(mi) > 0 and np.isfinite(mi[0]) and int(mi[0]) != 1:
            return False
        scv = self._sn.scv
        if scv is not None:
            for r in range(self.nclasses):
                for ist in (src_ist, q_ist):
                    v = scv[ist, r]
                    if np.isfinite(v) and abs(float(v) - 1.0) > 1e-6:
                        return False
        return True

    def _get_delay_stations(self) -> List[int]:
        """Get list of delay (infinite server) station indices.

        NodeType.DELAY is the whole test, as in MATLAB: a Queue whose scheduling
        is INF carries an infinite server count and is reported as a DELAY node
        by Network._refresh_node_mappings, so it reaches this list without a
        scheduling test of its own.
        """
        from ...api.sn.network_struct import NodeType

        delay_indices = []
        for i in range(self.nstations):
            if self.station_types and i < len(self.station_types):
                st = self.station_types[i]
                if st is not None:
                    st_val = st.value if hasattr(st, 'value') else int(st)
                    if st_val == NodeType.DELAY:
                        delay_indices.append(i)
        return delay_indices

    def _compute_cache_hit_miss_probs(self, XN: np.ndarray) -> None:
        """
        Compute and store cache hit/miss probabilities after MVA analysis.

        This method computes hit probabilities using cache analysis (cache_mva)
        based on the cache's gamma matrix and capacity. The probabilities are
        stored in both the Cache node (via set_result_hit_prob) and in the
        nodeparam structure for use by getAvgNode.

        Args:
            XN: System throughput per class

        References:
            MATLAB: solver_mva/@SolverMVA/runAnalyzer.m lines 150-167
        """
        from ...api.sn.network_struct import NodeType
        from ...api.cache import cache_xi_fp

        if self._sn is None or self._sn.nodetype is None:
            return

        # Find Cache nodes
        cache_indices = []
        for ind in range(self._sn.nnodes):
            if ind < len(self._sn.nodetype):
                node_type = self._sn.nodetype[ind]
                if node_type == NodeType.CACHE:
                    cache_indices.append(ind)

        if not cache_indices:
            return

        # Initialize nodeparam if needed
        if self._sn.nodeparam is None:
            self._sn.nodeparam = {}

        # Process each cache node
        model_nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        model_classes = self.model.get_classes() if hasattr(self.model, 'get_classes') else []

        for ind in cache_indices:
            # Get the Cache node from the model
            cache_node = None
            if ind < len(model_nodes):
                cache_node = model_nodes[ind]

            # Get hit/miss class mappings from cache node
            hitclass = []
            missclass = []

            if cache_node is not None:
                # Get hit/miss class indices for all classes
                for k in range(self.nclasses):
                    job_class = model_classes[k] if k < len(model_classes) else None
                    if job_class is not None:
                        h = cache_node.get_hit_class(job_class)
                        m = cache_node.get_miss_class(job_class)
                        if h is not None:
                            h_idx = model_classes.index(h) if h in model_classes else -1
                            hitclass.append(h_idx)
                        else:
                            hitclass.append(-1)
                        if m is not None:
                            m_idx = model_classes.index(m) if m in model_classes else -1
                            missclass.append(m_idx)
                        else:
                            missclass.append(-1)
                    else:
                        hitclass.append(-1)
                        missclass.append(-1)

            hitclass = np.array(hitclass)
            missclass = np.array(missclass)

            # Compute hit probability using cache analysis
            hitprob = np.zeros(len(hitclass))
            missprob = np.zeros(len(missclass))

            if cache_node is not None and hasattr(cache_node, 'get_gamma_matrix'):
                try:
                    # Get gamma matrix and cache parameters
                    gamma = cache_node.get_gamma_matrix(self.nclasses)
                    # Use item_level_cap which contains the capacity per cache level
                    m_levels = cache_node._item_level_cap if hasattr(cache_node, '_item_level_cap') else np.array([1])

                    # Run cache FPI to get miss probabilities (works for large caches)
                    xi, pi0, pij, it = cache_xi_fp(gamma, m_levels)

                    # overall cache hit rate = 1 - weighted (by access probability) average miss probability pi0.
                    access_probs = np.sum(gamma, axis=1)
                    access_probs = access_probs / np.sum(access_probs)  # Normalize

                    # Overall hit rate
                    overall_hit_rate = np.sum(access_probs * (1 - pi0))

                    # For each requesting class that has hit/miss classes, set the probability
                    for k in range(len(hitclass)):
                        h = hitclass[k]
                        m = missclass[k]
                        if h >= 0 and m >= 0:
                            hitprob[k] = overall_hit_rate
                            missprob[k] = 1 - overall_hit_rate

                except Exception as e:
                    # If cache analysis fails, fall back to default
                    pass

            # Store probabilities in cache node
            if cache_node is not None and hasattr(cache_node, 'set_result_hit_prob'):
                cache_node.set_result_hit_prob(hitprob)
                cache_node.set_result_miss_prob(missprob)

            # Store in nodeparam for getAvgNode helper functions
            class NodeParam:
                pass
            node_param = NodeParam()
            node_param.hitclass = hitclass
            node_param.missclass = missclass
            node_param.actualhitprob = hitprob
            node_param.actualmissprob = missprob
            self._sn.nodeparam[ind] = node_param

    def _is_polling_system(self) -> bool:
        """Check if this is a polling queueing system."""
        # polling-network detection: multiclass, no closed jobs, exactly Source/Queue/Sink node types, POLLING scheduling; mirrors MATLAB solver_mva.m:134.
        if self.nclasses <= 1:
            return False

        if self.network_type != 'open':
            return False

        # Check for Source-Queue-Sink topology
        from ...lang.base import NodeType, SchedStrategy
        if not hasattr(self, '_sn') or self._sn is None:
            return False

        nodetype = self._sn.nodetype
        if len(nodetype) != 3:
            return False

        # Check node types
        has_source = any(nt == NodeType.SOURCE for nt in nodetype)
        has_queue = any(nt == NodeType.QUEUE for nt in nodetype)
        has_sink = any(nt == NodeType.SINK for nt in nodetype)
        if not (has_source and has_queue and has_sink):
            return False

        # Check if queue uses POLLING scheduling
        sched = self._sn.sched
        if sched is None:
            return False

        for station_idx, sched_strategy in sched.items():
            sched_val = sched_strategy.value if hasattr(sched_strategy, 'value') else sched_strategy
            if sched_val == SchedStrategy.POLLING or sched_val == SchedStrategy.POLLING.value:
                return True

        return False

    def _run_polling_analysis(self):
        """Run polling system analysis."""
        from ...api.polling import polling_qsys_exhaustive, polling_qsys_gated, polling_qsys_1limited, polling_qsys_decrementing
        from ...api.mam import map_erlang
        from ...lang.base import NodeType
        from ...constants import PollingType, GlobalConstants

        R = self.nclasses

        # Get source, queue, and sink station indices
        nodetype = self._sn.nodetype
        source_ist = None
        queue_ist = None

        for i, nt in enumerate(nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else nt
            if nt_val == NodeType.SOURCE or nt_val == NodeType.SOURCE.value:
                source_ist = self._sn.nodeToStation[i]
            elif nt_val == NodeType.QUEUE or nt_val == NodeType.QUEUE.value:
                queue_ist = self._sn.nodeToStation[i]

        if source_ist is None or queue_ist is None:
            return None

        # Get arrival rates
        lambda_arr = np.zeros(R)
        for r in range(R):
            if self.rates[source_ist, r] > 0:
                lambda_arr[r] = self.rates[source_ist, r]

        # Get service rates
        mu = np.zeros(R)
        for r in range(R):
            if self.rates[queue_ist, r] > 0:
                mu[r] = self.rates[queue_ist, r]

        # Get polling type and switchover times from the model nodes
        polling_type = None
        polling_par = 1  # K value for K-limited
        switchover_dists = []

        # Find the queue node in the model
        queue_node = None
        for node in self.model.get_nodes():
            if hasattr(node, 'get_polling_type') and node.get_polling_type() is not None:
                queue_node = node
                pt = node.get_polling_type()
                polling_type = pt.value if hasattr(pt, 'value') else pt
                if hasattr(node, '_polling_k'):
                    polling_par = node._polling_k
                break

        if queue_node is None:
            return None

        # Get switchover time distributions
        for r in range(R):
            jobclass = self.model.get_classes()[r]
            if hasattr(queue_node, '_switchover') and jobclass in queue_node._switchover:
                switchover_dists.append(queue_node._switchover[jobclass])
            else:
                # Default: Immediate (zero switchover time)
                from ...distributions import Immediate
                switchover_dists.append(Immediate())

        # Convert distributions to MAP representations
        arvMAPs = []
        svcMAPs = []
        switchMAPs = []

        for r in range(R):
            # Arrival MAP - simple Poisson with rate lambda
            if lambda_arr[r] > 0:
                D0 = np.array([[-lambda_arr[r]]])
                D1 = np.array([[lambda_arr[r]]])
                arvMAPs.append((D0, D1))
            else:
                D0 = np.array([[-1e-10]])
                D1 = np.array([[1e-10]])
                arvMAPs.append((D0, D1))

            # service MAP approximated by an Erlang matching SCV (20-phase for near-deterministic, else ceil(1/SCV) phases), mirroring the JAR.
            if mu[r] > 0:
                mean_svc = 1.0 / mu[r]
                scv = self._sn.scv[queue_ist, r] if self._sn.scv is not None else 1.0
                if scv < GlobalConstants.CoarseTol:
                    # Deterministic or near-deterministic: use 20 phases
                    n_phases = 20
                else:
                    # Match SCV: for Erlang, SCV = 1/n, so n = ceil(1/SCV)
                    n_phases = max(1, int(np.ceil(1.0 / scv)))
                    n_phases = min(n_phases, 100)  # Cap at 100 phases
                D0, D1 = map_erlang(mean_svc, n_phases)
                svcMAPs.append((D0, D1))
            else:
                D0 = np.array([[-1.0]])
                D1 = np.array([[1.0]])
                svcMAPs.append((D0, D1))

            # Switchover MAP
            sw_dist = switchover_dists[r]
            if hasattr(sw_dist, 'isImmediate') and sw_dist.isImmediate():
                # Immediate switchover - use very high rate
                D0 = np.array([[-1e10]])
                D1 = np.array([[1e10]])
            elif hasattr(sw_dist, '_rate'):
                # Exponential
                rate = sw_dist._rate
                D0 = np.array([[-rate]])
                D1 = np.array([[rate]])
            elif hasattr(sw_dist, 'getMean'):
                # Use mean to create exponential approximation
                mean_sw = sw_dist.getMean()
                if mean_sw > 0:
                    rate = 1.0 / mean_sw
                    D0 = np.array([[-rate]])
                    D1 = np.array([[rate]])
                else:
                    D0 = np.array([[-1e10]])
                    D1 = np.array([[1e10]])
            else:
                # Default to immediate
                D0 = np.array([[-1e10]])
                D1 = np.array([[1e10]])
            switchMAPs.append((D0, D1))

        # Call appropriate polling analysis function
        polling_type_val = PollingType.EXHAUSTIVE.value if polling_type is None else polling_type

        if polling_type_val == PollingType.EXHAUSTIVE or polling_type_val == PollingType.EXHAUSTIVE.value:
            W = polling_qsys_exhaustive(arvMAPs, svcMAPs, switchMAPs)
        elif polling_type_val == PollingType.GATED or polling_type_val == PollingType.GATED.value:
            W = polling_qsys_gated(arvMAPs, svcMAPs, switchMAPs)
        elif polling_type_val == PollingType.KLIMITED or polling_type_val == PollingType.KLIMITED.value:
            if polling_par == 1:
                W = polling_qsys_1limited(arvMAPs, svcMAPs, switchMAPs)
            else:
                # For K > 1, fall back to approximation
                return None
        elif polling_type_val == PollingType.DECREMENTING or polling_type_val == PollingType.DECREMENTING.value:
            W = polling_qsys_decrementing(arvMAPs, svcMAPs, switchMAPs)
        else:
            return None

        # Compute response times: R = W + 1/mu
        R_queue = np.zeros(R)
        for r in range(R):
            if mu[r] > 0:
                R_queue[r] = W[r] + 1.0 / mu[r]
            else:
                R_queue[r] = W[r]

        # Build result arrays
        QN = np.zeros((self.nstations, R))
        UN = np.zeros((self.nstations, R))
        RN = np.zeros((self.nstations, R))
        TN = np.zeros((self.nstations, R))
        AN = np.zeros((self.nstations, R))
        XN = np.zeros(R)

        k = 1  # Number of servers (polling is single server)

        # Set metrics for source station
        TN[source_ist, :] = lambda_arr

        # Set metrics for queue station
        RN[queue_ist, :] = R_queue
        TN[queue_ist, :] = lambda_arr
        AN[queue_ist, :] = lambda_arr
        UN[queue_ist, :] = lambda_arr / mu / k
        QN[queue_ist, :] = lambda_arr * R_queue  # Little's law
        XN = lambda_arr.copy()

        # Compute residence times from response times (WN = RN * V)
        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self._sn, RN, None)

        # Store results
        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
            'WN': WN,
            'lG': 0,
            'runtime': 0,
            'lastiter': 1,
        }

        return self._result

    def _sizebased_sched(self):
        """The queue's size-based discipline, or None when it has none.

        SRPT, PSJF, FB, LRPT and SETF are served by the Wierman and
        Harchol-Balter response times (SIGMETRICS 2003); the generic AMVA path
        carries no size-based term at all and would solve the station
        size-blind, which is not what SRPT means. Mirrors MATLAB
        mvaDispatch.m's isSizeBasedPolicy branch.
        """
        from ...lang.base import NodeType
        if self.network_type != 'open' or self._sn is None:
            return None
        nodetype = self._sn.nodetype
        if len(nodetype) != 3:
            return None
        has_source = any(nt == NodeType.SOURCE for nt in nodetype)
        has_queue = any(nt == NodeType.QUEUE for nt in nodetype)
        has_sink = any(nt == NodeType.SINK for nt in nodetype)
        if not (has_source and has_queue and has_sink):
            return None
        # Match by NAME, never by value: the two Python SchedStrategy enums
        # DISAGREE above 35 (lang.base has FSP=36, PAS=37, OI=38, then the
        # size-based names appended at 39-42; constants has PSJF=36, FB=37,
        # LAS=38, LRPT=39), and _normalize_sched_strategy reconciles them by
        # name. A value comparison read a stored PSJF as LRPT and missed FB
        # entirely, which then fell through to the size-blind AMVA path.
        sized = ('SRPT', 'PSJF', 'FB', 'LAS', 'LRPT', 'SETF')
        sched = self._sn.sched
        if sched is None:
            return None
        for _, sched_strategy in sched.items():
            name = getattr(sched_strategy, 'name', None)
            if name in sized:
                return 'FB' if name == 'LAS' else name
        return None

    def _run_sizebased_analysis(self):
        """M/G/1 with size-based scheduling: SRPT, PSJF, FB (LAS), LRPT, SETF.

        Port of MATLAB solver_mva_qsys_sizebased_analyzer.m. Twin of the JAR
        Solver_mva_qsys_sizebased_analyzer.

        NOTE ON THE INDEX SPACE. ``sn.visits`` is indexed by CHAIN, not by
        station: reading the Source's station index takes chain 1, so on the
        multiclass models this analyzer exists for every class beyond the first
        would get the visit of a chain it does not belong to, which is zero.
        The chain of each class is looked up explicitly.
        """
        from ...api.qsys import (qsys_mg1_srpt, qsys_mg1_psjf, qsys_mg1_fb,
                                 qsys_mg1_lrpt, qsys_mg1_setf)
        from ...lang.base import NodeType
        from ...api.io.logging import line_warning

        sched_type = self._sizebased_sched()
        if sched_type is None:
            return None
        R = self.nclasses
        nodetype = self._sn.nodetype
        source_ist = None
        queue_ist = None
        for i, nt in enumerate(nodetype):
            nt_val = nt.value if hasattr(nt, 'value') else nt
            if nt_val == NodeType.SOURCE or nt_val == NodeType.SOURCE.value:
                source_ist = self._sn.nodeToStation[i]
            elif nt_val == NodeType.QUEUE or nt_val == NodeType.QUEUE.value:
                queue_ist = self._sn.nodeToStation[i]
        if source_ist is None or queue_ist is None:
            return None
        queue_isf = int(self._sn.stationToStateful[queue_ist])

        chains = np.asarray(self._sn.chains)
        chain_of = np.zeros(R, dtype=int)
        for k in range(R):
            nz = np.nonzero(chains[:, k])[0]
            chain_of[k] = int(nz[0]) if nz.size else 0

        lambda_arr = np.zeros(R)
        mu = np.zeros(R)
        cs = np.zeros(R)
        visits = np.zeros(R)
        for k in range(R):
            visits[k] = float(self._sn.visits[chain_of[k]][queue_isf, k])
            lambda_arr[k] = float(self.rates[source_ist, k]) * visits[k]
            mu[k] = float(self.rates[queue_ist, k])
            scv = float(self._sn.scv[queue_ist, k])
            cs[k] = np.sqrt(scv) if np.isfinite(scv) and scv > 0 else 1.0

        if np.any(lambda_arr <= 0) or np.any(mu <= 0):
            raise ValueError('solver_mva_qsys_sizebased_analyzer: invalid arrival or '
                             'service rates (must be positive).')
        rho = float(np.sum(lambda_arr / mu))
        if rho >= 1.0:
            line_warning('solver_mva_qsys_sizebased_analyzer',
                         'System is unstable (rho = %.4f >= 1).' % rho)

        if sched_type == 'SRPT':
            W, _ = qsys_mg1_srpt(lambda_arr, mu, cs)
        elif sched_type == 'PSJF':
            W, _ = qsys_mg1_psjf(lambda_arr, mu, cs)
        elif sched_type == 'FB':
            W, _ = qsys_mg1_fb(lambda_arr, mu, cs)
        elif sched_type == 'LRPT':
            W, _ = qsys_mg1_lrpt(lambda_arr, mu, cs)
        else:
            W, _ = qsys_mg1_setf(lambda_arr, mu, cs)
        W = np.asarray(W, dtype=float).flatten()

        QN = np.zeros((self.nstations, R))
        UN = np.zeros((self.nstations, R))
        RN = np.zeros((self.nstations, R))
        TN = np.zeros((self.nstations, R))
        AN = np.zeros((self.nstations, R))
        TN[source_ist, :] = lambda_arr
        AN[source_ist, :] = lambda_arr
        RN[queue_ist, :] = W * visits
        TN[queue_ist, :] = lambda_arr
        AN[queue_ist, :] = lambda_arr
        UN[queue_ist, :] = lambda_arr / mu
        QN[queue_ist, :] = lambda_arr * W
        XN = lambda_arr.copy()

        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self._sn, RN, None)

        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
            'WN': WN,
            'lG': 0,
            'runtime': 0,
            'lastiter': 1,
        }
        return self._result

    def _is_cache_only_network(self) -> bool:
        """Check if this is a cache-only network (Source-Cache-Sink)."""
        if self._sn is None:
            return False

        from ...api.sn.network_struct import NodeType

        # Check for open network only (no closed jobs)
        if hasattr(self._sn, 'nclosedjobs') and self._sn.nclosedjobs > 0:
            return False

        # Check if all jobs are open (infinite population)
        if self._sn.njobs is not None:
            if not np.all(np.isinf(self._sn.njobs)):
                return False

        # Check node types - must have exactly Source, Cache, Sink
        if self._sn.nodetype is None:
            return False

        node_types = list(self._sn.nodetype)
        has_source = NodeType.SOURCE in node_types
        has_cache = NodeType.CACHE in node_types
        has_sink = NodeType.SINK in node_types

        # Count each type
        num_sources = node_types.count(NodeType.SOURCE)
        num_caches = node_types.count(NodeType.CACHE)
        num_sinks = node_types.count(NodeType.SINK)

        # Must have exactly 1 source, 1 cache, 1 sink
        if num_sources != 1 or num_caches != 1 or num_sinks != 1:
            return False

        # Total nodes must be 3
        if len(node_types) != 3:
            return False

        return has_source and has_cache and has_sink

    def _has_cache_with_class_switching(self) -> bool:
        """Check if network has cache nodes with class switching (hit/miss classes)."""
        if self._sn is None:
            return False

        from ...api.sn.network_struct import NodeType

        if self._sn.nodetype is None:
            return False

        # Check for cache nodes
        has_cache = False
        for nt in self._sn.nodetype:
            if nt == NodeType.CACHE:
                has_cache = True
                break

        if not has_cache:
            return False

        # Check if cache has hit/miss class switching configured
        model_nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        model_classes = self.model.get_classes() if hasattr(self.model, 'get_classes') else []

        for node in model_nodes:
            if hasattr(node, '_hit_class') and hasattr(node, '_miss_class'):
                # Check if any class has hit/miss configured
                for job_class in model_classes:
                    h = node._hit_class.get(job_class)
                    m = node._miss_class.get(job_class)
                    if h is not None and m is not None:
                        return True

        return False

    def _update_cache_routing_and_visits(self):
        """
        Update routing matrix with cache hit/miss probabilities and refresh visits.

        This implements the routing update logic from MATLAB's solver_mva_cacheqn_analyzer.
        The key steps are:
        1. Compute cache hit/miss probabilities
        2. Update sn.rtnodes to route input class to hit/miss classes with computed probabilities
        3. Recompute sn.rt using stochastic complement
        4. Refresh sn.visits with the updated routing

        References:
            MATLAB: solver_mva_cacheqn_analyzer.m lines 77-91
        """
        from ...api.sn.network_struct import NodeType
        from ...api.cache import cache_xi_fp, cache_gamma_lp, cache_miss_fpi
        from ...api.mc.dtmc import dtmc_stochcomp
        from ...api.sn.transforms import sn_refresh_visits

        sn = self._sn
        if sn is None:
            return

        I = sn.nnodes
        K = sn.nclasses

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

        # Make a copy of rtnodes to modify
        rtnodes = sn.rtnodes.copy() if sn.rtnodes is not None else None
        if rtnodes is None:
            return

        # Compute hit/miss probabilities for each cache and update routing
        for ind in cache_indices:
            cache_node = model_nodes[ind] if ind < len(model_nodes) else None
            if cache_node is None:
                continue

            ch = sn.nodeparam.get(ind) if sn.nodeparam else None
            if ch is None:
                continue

            hitclass = getattr(ch, 'hitclass', None)
            missclass = getattr(ch, 'missclass', None)
            if hitclass is None or missclass is None:
                continue

            # Find input classes (classes that have hit/miss mappings)
            input_classes = []
            for r in range(K):
                if r < len(hitclass) and r < len(missclass):
                    hc = int(hitclass[r]) if hitclass[r] >= 0 else -1
                    mc = int(missclass[r]) if missclass[r] >= 0 else -1
                    if hc >= 0 and mc >= 0:
                        input_classes.append(r)

            if not input_classes:
                continue

            # Compute hit/miss probabilities using FPI (like MATLAB lines 68-71)
            try:
                if hasattr(cache_node, 'get_gamma_matrix'):
                    gamma = cache_node.get_gamma_matrix(K)
                    m_levels = cache_node._item_level_cap if hasattr(cache_node, '_item_level_cap') else np.array([1])

                    # Use cache_miss_fpi to compute miss rates (MATLAB line 70)
                    # For now, use cache_xi_fp and derive hit/miss from pi0
                    xi, pi0, pij, it_fp = cache_xi_fp(gamma, m_levels)
                    access_probs = np.sum(gamma, axis=1)
                    total = np.sum(access_probs)
                    if total > 0:
                        access_probs = access_probs / total
                    overall_hit_rate = np.sum(access_probs * (1 - pi0))
                    overall_miss_rate = 1 - overall_hit_rate
                else:
                    # Default probabilities
                    overall_hit_rate = 0.5
                    overall_miss_rate = 0.5
            except Exception:
                overall_hit_rate = 0.5
                overall_miss_rate = 0.5

            hitprob = np.zeros(K)
            missprob = np.zeros(K)
            for r in input_classes:
                hitprob[r] = overall_hit_rate
                missprob[r] = overall_miss_rate

            # Update routing matrix (MATLAB lines 78-86)
            # For each input class, route from Cache to connected nodes with hit/miss probabilities
            for r in input_classes:
                hc = int(hitclass[r])
                mc = int(missclass[r])

                # Zero out the row for input class at cache
                rtnodes[(ind) * K + r, :] = 0

                # Find connected nodes
                for jnd in range(I):
                    if sn.connmatrix is not None and ind < sn.connmatrix.shape[0] and jnd < sn.connmatrix.shape[1]:
                        if sn.connmatrix[ind, jnd]:
                            # Route to hit class with hit probability
                            rtnodes[(ind) * K + r, (jnd) * K + hc] = hitprob[r]
                            # Route to miss class with miss probability
                            rtnodes[(ind) * K + r, (jnd) * K + mc] = missprob[r]

            # Set hit/miss probs on cache node for result reporting
            if hasattr(cache_node, 'set_result_hit_prob'):
                cache_node.set_result_hit_prob(hitprob)
                cache_node.set_result_miss_prob(missprob)

            # exact hit/miss probs stored in sn.nodeparam so getAvgNode avoids stale values from a previous run (MATLAB sn is value-type, so python-only need).
            ch.actualhitprob = hitprob.copy()
            ch.actualmissprob = missprob.copy()

        # Update sn.rtnodes
        sn.rtnodes = rtnodes

        # Recompute sn.rt using stochastic complement (MATLAB line 87)
        stateful_nodes = np.where(sn.isstateful)[0] if sn.isstateful is not None else np.arange(I)
        stateful_nodes_classes = []
        for sf_idx in stateful_nodes:
            for k in range(K):
                stateful_nodes_classes.append(sf_idx * K + k)
        stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

        try:
            new_rt = dtmc_stochcomp(rtnodes, stateful_nodes_classes)
            sn.rt = new_rt
            # CRITICAL: Also update rt_visits since sn_refresh_visits uses it
            # (if it exists) instead of sn.rt
            if hasattr(sn, 'rt_visits') and sn.rt_visits is not None:
                sn.rt_visits = new_rt.copy()
        except Exception:
            pass

        # Refresh visits (MATLAB line 89)
        try:
            sn_refresh_visits(sn)
        except Exception:
            pass

        # Recompute self.demands with updated visits
        # This replicates the logic from __init__ lines 185-200
        try:
            # Reset demands to 1/rates
            nonzero = self.rates > 0
            self.demands = np.zeros_like(self.rates)
            self.demands[nonzero] = 1.0 / self.rates[nonzero]

            # Apply visits to demands (D = V / mu)
            if hasattr(sn, 'visits') and sn.visits is not None:
                visits_combined = np.zeros_like(self.demands)
                stationToStateful = np.asarray(sn.stationToStateful, dtype=int).flatten()
                for chain_id, visits_chain in sn.visits.items():
                    if isinstance(visits_chain, np.ndarray):
                        for station_idx in range(self.nstations):
                            stateful_idx = stationToStateful[station_idx]
                            if stateful_idx < visits_chain.shape[0]:
                                visits_combined[station_idx, :] += visits_chain[stateful_idx, :]
                self.demands = self.demands * visits_combined
        except Exception:
            pass

    def _run_cache_qn_analysis(self):
        """
        Run MVA analysis for queueing networks with cache nodes.

        This implements the cache QN analyzer similar to MATLAB's
        solver_mva_cacheqn_analyzer. It computes cache hit/miss probabilities,
        updates the routing matrix, converts Cache nodes to ClassSwitch, and
        returns None to let the standard MVA analyzer handle the network with
        the updated routing.

        References:
            MATLAB: solver_mva_cacheqn_analyzer.m
        """
        from ...api.sn.network_struct import NodeType

        sn = self._sn
        if sn is None:
            return None

        I = sn.nnodes

        # Find cache nodes
        cache_indices = []
        for ind in range(I):
            if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
                cache_indices.append(ind)

        if not cache_indices:
            return None

        # Update routing matrix with cache hit/miss probabilities
        # This computes hit/miss probabilities, updates rtnodes, refreshes visits
        self._update_cache_routing_and_visits()

        # per-item occupancy: LRU uses scale-invariant TTL; RR/FIFO use the exact product-form recursion, tractable only below 10 items (else NaN + warning).
        self._compute_cache_item_prob(cache_indices)

        # Cache nodes converted to ClassSwitch so standard MVA handles topology; indices saved to restore nodetype after solving; mirrors MATLAB line 38.
        self._cache_indices = cache_indices
        for ind in cache_indices:
            sn.nodetype[ind] = NodeType.CLASSSWITCH

        # returning None continues with standard MVA on the updated routing, which handles the full topology including downstream LN sublayer sync calls.
        return None

    def _compute_cache_item_prob(self, cache_indices):
        """Populate per-item cache occupancy on cache nodes for getAvgItemTable.

        LRU uses the scalable TTL algorithm; RR/FIFO use the exact product-form
        recursion, skipped (NaN, with a warning) for more than 10 items.
        """
        from ...api.cache import cache_gamma_lp, cache_ttl_lrua, cache_prob_erec
        from ...lang.base import ReplacementStrategy
        sn = self._sn
        model_nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
        for ind in cache_indices:
            cache_node = model_nodes[ind] if ind < len(model_nodes) else None
            ch = sn.nodeparam.get(ind) if sn.nodeparam else None
            if cache_node is None or ch is None:
                continue
            cap_raw = getattr(ch, 'itemcap', None)
            if cap_raw is None:
                cap_raw = getattr(ch, 'cap', [1])
            if np.isscalar(cap_raw) or (isinstance(cap_raw, np.ndarray) and cap_raw.ndim == 0):
                m = np.array([int(cap_raw)])
            else:
                m = np.asarray(cap_raw).ravel()
            n = int(getattr(ch, 'nitems', 0))
            h = len(m)
            if n <= 0 or h <= 0:
                continue
            replacestrat = getattr(ch, 'replacestrat', ReplacementStrategy.RR)
            pread = getattr(ch, 'pread', None)
            # Build read-rate tensor [R x n x (h+1)]; occupancy is scale-invariant
            # so a unit reference rate is sufficient.
            R = self.nclasses
            lambd = np.zeros((R, n, h + 1))
            for v in range(R):
                if pread is not None and v < len(pread) and pread[v] is not None:
                    pread_v = np.asarray(pread[v]).ravel()
                    for k in range(min(n, len(pread_v))):
                        for l in range(h + 1):
                            lambd[v, k, l] = pread_v[k]
            Rcost = getattr(ch, 'accost', None)
            if Rcost is None:
                def _default_routing(h):
                    mat = np.diag(np.ones(h), 1)
                    mat[h, h] = 1.0
                    return mat
                Rcost = [[_default_routing(h) for _ in range(n)] for _ in range(R)]
            if replacestrat == ReplacementStrategy.LRU:
                item_prob = cache_ttl_lrua(lambd, Rcost, m)
            elif n > 10:
                line_warning('solver_mva_cacheqn_analyzer',
                             'Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.' % n)
                item_prob = np.full((n, h + 1), np.nan)
            else:
                gamma, _, _, _, _ = cache_gamma_lp(lambd, Rcost)
                item_prob = cache_prob_erec(gamma, m)
            if hasattr(cache_node, 'set_result_item_prob'):
                cache_node.set_result_item_prob(item_prob)
            ch.actualitemprob = item_prob

    def _run_retrieval_analysis(self):
        """Run FPI analysis for a delayed-hit cache with a retrieval system."""
        from ...api.sn.network_struct import NodeType
        from ...api.retrieval.analyzers import (
            solver_mva_retrieval_analyzer, solver_mva_cacheqn_retrieval_analyzer, _has_source
        )
        sn = self._sn
        # open (Source) cache uses product-form FPI analyzer; closed uses da_cacheqn_retrieval (relabels Cache to ClassSwitch, index in res.cache_idx).
        if _has_source(sn):
            res = solver_mva_retrieval_analyzer(sn, self.options)
            cache_indices = [ind for ind in range(sn.nnodes)
                             if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE
                             and ind in sn.nodeparam]
        else:
            res = solver_mva_cacheqn_retrieval_analyzer(sn, self.options)
            cache_indices = [res.cache_idx]
        # set cache node results (hit/miss/latency)
        for ind in cache_indices:
            if ind in sn.nodeparam:
                cp = sn.nodeparam[ind]
                cp.actualhitprob = res.hitprob[0, :]
                cp.actualmissprob = res.missprob[0, :]
                cp.actualdelayedhitprob = res.delayedprob[0, :]
                cp.actualhitproblist = res.hitproblist
                cp.actualitemprob = res.itemprob
                if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                    node = self.model._nodes[ind]
                    if hasattr(node, 'set_result_hit_prob'):
                        node.set_result_hit_prob(res.hitprob[0, :])
                    if hasattr(node, 'set_result_miss_prob'):
                        node.set_result_miss_prob(res.missprob[0, :])
                    if hasattr(node, 'set_result_delayed_hit_prob'):
                        node.set_result_delayed_hit_prob(res.delayedprob[0, :])
                    if hasattr(node, 'set_result_hit_prob_list'):
                        node.set_result_hit_prob_list(res.hitproblist)
                    if res.itemprob is not None and hasattr(node, 'set_result_item_prob'):
                        node.set_result_item_prob(res.itemprob)
                    if hasattr(node, 'set_result_residt'):
                        node.set_result_residt(res.expected_latency[0, :])
                break
        M = self.nstations
        self._result = {
            'QN': res.QN, 'UN': res.UN, 'RN': res.RN, 'TN': res.TN,
            'CN': res.RN.copy(), 'XN': np.asarray(res.XN).ravel(),
            'AN': np.zeros((M, self.nclasses)), 'WN': res.RN.copy(),
            'lG': res.lG, 'runtime': res.runtime, 'iter': 1, 'method': 'fpi'
        }
        return self._result

    def _run_cache_analysis(self):
        """Run specialized cache analysis for Source-Cache-Sink models."""
        from ...api.sn.network_struct import NodeType
        from ...api.cache import cache_prob_fpi, cache_ttl_lrua, cache_gamma_lp
        from ...lang.base import ReplacementStrategy

        sn = self._sn
        R = self.nclasses

        # Find source station and get arrival rates
        source_ist = None
        for ist in range(sn.nstations):
            ind = sn.stationToNode[ist]
            if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.SOURCE:
                source_ist = ist
                break

        if source_ist is None:
            return None

        source_rate = sn.rates[source_ist, :].copy()
        source_rate = np.nan_to_num(source_rate, nan=0.0)

        # Find cache node
        cache_ind = None
        cache_param = None
        for ind in range(sn.nnodes):
            if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
                cache_ind = ind
                if sn.nodeparam is not None and ind in sn.nodeparam:
                    cache_param = sn.nodeparam[ind]
                break

        if cache_ind is None or cache_param is None:
            return None

        # Get cache parameters
        # Use itemcap for full capacity vector if available, otherwise fall back to cap
        cap_raw = getattr(cache_param, 'itemcap', None)
        if cap_raw is None:
            cap_raw = getattr(cache_param, 'cap', [1])
        if np.isscalar(cap_raw) or (isinstance(cap_raw, np.ndarray) and cap_raw.ndim == 0):
            m = np.array([int(cap_raw)])
        else:
            m = np.asarray(cap_raw).ravel()
        n = int(getattr(cache_param, 'nitems', 1))
        h = len(m)  # number of cache levels

        # Get replacement strategy
        replacestrat = getattr(cache_param, 'replacestrat', ReplacementStrategy.RR)

        # Get read probabilities (pread)
        pread = getattr(cache_param, 'pread', None)

        # Build lambda matrix (u x n x h+1)
        # u = number of classes, n = number of items, h+1 = levels (including miss level)
        lambd = np.zeros((R, n, h + 1))
        for v in range(R):
            if pread is not None and v < len(pread) and pread[v] is not None:
                pread_v = np.asarray(pread[v]).ravel()
                for k in range(min(n, len(pread_v))):
                    for l in range(h + 1):
                        lambd[v, k, l] = source_rate[v] * pread_v[k]

        # Get access cost (Rcost)
        Rcost = getattr(cache_param, 'accost', None)
        if Rcost is None:
            # Default routing: linear cache hierarchy
            # MATLAB creates: diag(ones(1,nLevels),1) with last diagonal element set to 1
            # This means: level 0 → 1 → 2 → ... → h, with level h staying at h
            # Rcost should be a 2D structure: Rcost[v][i] is the routing matrix for user v and item i
            def create_default_routing(h):
                """Create default linear cache routing matrix (h+1 x h+1)."""
                mat = np.diag(np.ones(h), 1)  # Super-diagonal with 1s
                mat[h, h] = 1.0  # Last level stays at last level
                return mat
            Rcost = [[create_default_routing(h) for _ in range(n)] for _ in range(R)]

        # Compute gamma using cache_gamma_lp
        try:
            gamma, _, _, _, _ = cache_gamma_lp(lambd, Rcost)
        except Exception:
            # Fall back to simple gamma computation
            gamma = np.zeros((n, h))
            for k in range(n):
                for l in range(h):
                    gamma[k, l] = np.sum(lambd[:, k, l])

        # Choose algorithm based on replacement strategy
        pij = None
        try:
            if replacestrat in (ReplacementStrategy.RR, ReplacementStrategy.FIFO):
                # Use Fixed Point Iteration method
                pij = cache_prob_fpi(gamma, m)
            elif replacestrat == ReplacementStrategy.LRU:
                # MMAP: per-mark MAPs make requests non-IRM, LRU(m)-MAP TTL (Gast-Van Houdt 2017); i.i.d. sequence-exact TTL; mirrors solver_mva_cache_analyzer.
                D0c = None
                markidx = getattr(sn, 'markidx', None)
                if (markidx is not None and source_ist < markidx.shape[0]
                        and np.any(markidx[source_ist, :] > 0)):
                    carrier = int(np.where(markidx[source_ist, :] > 0)[0][0])
                    Dcell = sn.proc[source_ist][carrier]
                    D0 = np.atleast_2d(np.asarray(Dcell[0], dtype=float))
                    D1agg = np.atleast_2d(np.asarray(Dcell[1], dtype=float))
                    d = D0.shape[0]
                    allmarked = True
                    D0c, D1c = [], []
                    for k in range(n):
                        D1k = np.zeros((d, d))
                        for v in range(R):
                            if pread is not None and v < len(pread) and pread[v] is not None:
                                pv = np.asarray(pread[v]).ravel()
                                pk = float(pv[k]) if k < len(pv) else 0.0
                                if markidx[source_ist, v] > 0:
                                    Dm = np.atleast_2d(np.asarray(
                                        Dcell[1 + int(markidx[source_ist, v])], dtype=float))
                                    D1k = D1k + pk * Dm
                                elif source_rate[v] > 0 and pk > 0:
                                    allmarked = False  # unmarked reader mixed in
                        D0c.append(D0 + D1agg - D1k)
                        D1c.append(D1k)
                    if not allmarked:
                        D0c = None
                if D0c is not None:
                    from ...api.cache import cache_ttl_lrum_map
                    pij, _ = cache_ttl_lrum_map(D0c, D1c, m)
                else:
                    # Use TTL-based LRU approximation
                    pij = cache_ttl_lrua(lambd, Rcost, m)
            elif replacestrat == ReplacementStrategy.HLRU:
                # h-LRU / LRU(m) characteristic-time approximation (linear
                # list topology; mirrors MATLAB solver_mva_cache_analyzer)
                from ...api.cache import cache_ttl_hlru
                pij = cache_ttl_hlru(lambd, m)
            else:
                # Default to FPI
                pij = cache_prob_fpi(gamma, m)
        except Exception:
            # Fall back to simple approximation
            total_cap = np.sum(m)
            hit_prob = min(total_cap / n, 1.0) if n > 0 else 0.0
            pij = np.zeros((n, h + 1))
            pij[:, 0] = 1 - hit_prob  # miss probability
            pij[:, 1:] = hit_prob / h if h > 0 else 0

        # Compute miss rates per class
        miss_rate = np.zeros(R)
        for v in range(R):
            if pread is not None and v < len(pread) and pread[v] is not None:
                pread_v = np.asarray(pread[v]).ravel()
                for k in range(min(n, len(pread_v), pij.shape[0])):
                    miss_rate[v] += source_rate[v] * pread_v[k] * pij[k, 0]

        # per-list occupancy needs its own distribution: LRU's pij is one; RR/FIFO derive it from exact cache_prob_erec (below 10 items, else NaN + warning).
        from ...api.cache import cache_prob_erec
        if replacestrat == ReplacementStrategy.LRU:
            item_prob = pij
        elif n > 10:
            line_warning('solver_mva_cache_analyzer',
                         'Per-item cache occupancy (getAvgItemTable) requires the exact algorithm for RR/FIFO and is skipped for caches with more than 10 items (%d items); reporting NaN.' % n)
            item_prob = np.full((n, h + 1), np.nan)
        else:
            item_prob = cache_prob_erec(gamma, m)

        # Get hit/miss class mappings
        hitclass = np.asarray(getattr(cache_param, 'hitclass', [])).astype(int)
        missclass = np.asarray(getattr(cache_param, 'missclass', [])).astype(int)

        # Initialize throughput array
        XN = np.zeros(R)

        # Set throughputs for hit/miss classes
        for r in range(min(len(hitclass), len(missclass))):
            h_idx = hitclass[r] if r < len(hitclass) else -1
            m_idx = missclass[r] if r < len(missclass) else -1

            if h_idx >= 0 and h_idx < R and m_idx >= 0 and m_idx < R:
                XN[m_idx] = miss_rate[r]
                XN[h_idx] = source_rate[r] - miss_rate[r]

        # Set hit/miss probabilities on cache node
        hit_prob_arr = np.zeros(R)
        miss_prob_arr = np.zeros(R)
        for r in range(R):
            if source_rate[r] > 0:
                miss_prob_arr[r] = miss_rate[r] / source_rate[r]
                hit_prob_arr[r] = 1 - miss_prob_arr[r]

        # Store in nodeparam
        cache_param.actualhitprob = hit_prob_arr
        cache_param.actualmissprob = miss_prob_arr

        # Set on cache node in model
        if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
            cache_node_obj = self.model._nodes[cache_ind]
            if hasattr(cache_node_obj, 'set_result_hit_prob'):
                cache_node_obj.set_result_hit_prob(hit_prob_arr)
            if hasattr(cache_node_obj, 'set_result_miss_prob'):
                cache_node_obj.set_result_miss_prob(miss_prob_arr)
            if item_prob is not None and item_prob.shape[1] == h + 1 \
                    and hasattr(cache_node_obj, 'set_result_item_prob'):
                cache_node_obj.set_result_item_prob(item_prob)

        # Build result arrays
        M = self.nstations
        QN = np.zeros((M, R))
        UN = np.zeros((M, R))
        RN = np.zeros((M, R))
        TN = np.zeros((M, R))
        AN = np.zeros((M, R))

        # Set source throughput
        TN[source_ist, :] = source_rate

        # Set cache station throughput using XN
        cache_ist = sn.nodeToStation[cache_ind]
        if cache_ist >= 0 and cache_ist < M:
            for r in range(R):
                if XN[r] > 0:
                    TN[cache_ist, r] = XN[r]
                elif source_rate[r] > 0:
                    TN[cache_ist, r] = source_rate[r]

        # Compute residence times from response times (WN = RN * V)
        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self._sn, RN, None)

        # Store result
        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'CN': RN.copy(),
            'XN': XN,
            'AN': AN,
            'WN': WN,
            'lG': np.nan,
            'runtime': 0.0,
            'iter': 1,
            'method': 'cache'
        }

        return self._result

    def _resolve_oi_path(self, method):
        """
        Classify the model's order-independent (OI/PAS) content for method.

        Returns (noi_idx, oi_exact), where noi_idx is the index of the OI station
        (-1 if none) and oi_exact says whether the exact order-independent
        analyzer applies. Raises ValueError when an OI/PAS station is present but
        falls outside that analyzer's scope, since AMVA cannot represent it.
        """
        # OI detection mirrors the NC-oi path (empty/zero swap graph + a service-rate function), not a nodeparam flag.
        noi_idx = -1
        if self._sn is not None and getattr(self._sn, 'njobs', None) is not None \
                and not np.any(np.isinf(np.asarray(self._sn.njobs, dtype=float))):
            from .solver_mva_oi_analyzer import find_oi_station
            noi_idx = find_oi_station(self._sn)

        # OI/PAS admissible only when every other station is product-form (the NC-oi gate); sn.sched carries lang.base.SchedStrategy members, matched by name.
        has_oi_station = False
        if self._sn is not None and getattr(self._sn, 'sched', None) is not None:
            from ...lang.base import SchedStrategy as _SSg
            _sched = self._sn.sched
            for _i in range(int(self._sn.nstations)):
                _si = _sched[_i]
                _nm = getattr(_si, 'name', None)
                if _nm is None:
                    try:
                        _nm = _SSg(int(_si)).name
                    except (ValueError, TypeError):
                        _nm = None
                if _nm in ('OI', 'PAS'):
                    has_oi_station = True
                    break
        oi_exact = False
        if noi_idx >= 0 and method in ['exact', 'default']:
            from ..solver_nc.solver_nc_oi_analyzer import nc_is_oi_model as _nc_is_oi_model
            oi_exact = _nc_is_oi_model(self._sn)

        if has_oi_station and not oi_exact:
            # OI/PAS stations need the exact path; see _kb/06-solver-catalog.md MVA Exact-only dispatch for OI/PAS and fork-join.
            raise ValueError(
                "SolverMVA supports order-independent (OI) and pass-and-swap (PAS) stations only\n"
                "through its exact order-independent analyzer, which requires method 'default' or\n"
                "'exact' (got '%s'), an empty/zero swap graph at every OI/PAS station, a closed\n"
                "model, and every other station to be product-form (INF, PS, LCFS-PR, SIRO,\n"
                "or class-independent-rate FCFS). Use SolverCTMC or SolverLDES for this model."
                % method)
        return noi_idx, oi_exact

    def _marie_inf_mask(self):
        from ...lang.base import SchedStrategy
        m = np.zeros(self.nstations, dtype=bool)
        sched = self.sched
        if sched is None:
            return m
        for i in range(self.nstations):
            sv = sched.get(i) if isinstance(sched, dict) else (
                sched[i] if i < len(sched) else None)
            if sv is None:
                continue
            m[i] = (sv == SchedStrategy.INF or
                    (hasattr(sv, 'value') and sv.value == SchedStrategy.INF.value) or
                    (isinstance(sv, int) and sv == SchedStrategy.INF.value))
        return m

    def _sjn_station_mask(self):
        """Stations scheduling by non-preemptive shortest job next."""
        from ...lang.base import SchedStrategy
        m = np.zeros(self.nstations, dtype=bool)
        sched = self.sched
        if sched is None:
            return m
        for i in range(self.nstations):
            sv = sched.get(i) if isinstance(sched, dict) else (
                sched[i] if i < len(sched) else None)
            if sv is None:
                continue
            m[i] = (sv == SchedStrategy.SJF or
                    (hasattr(sv, 'value') and sv.value == SchedStrategy.SJF.value) or
                    (isinstance(sv, int) and sv == SchedStrategy.SJF.value))
        return m

    def _run_sjn(self):
        """Closed models with shortest-job-next stations; see _kb/06-solver-catalog.md."""
        import time as _t
        from ...lang.base import SchedStrategy
        from ...api.pfqn.sjn import pfqn_mvasjn, pfqn_amvasjn, SjnOptions, SjnStarvationError
        from ...api.sn.demands import sn_get_demands_chain
        from ...api.sn.deaggregate import sn_deaggregate_chain_results
        t0 = _t.time()
        sn = self._sn
        M = self.nstations
        C = sn.nchains
        dem = sn_get_demands_chain(sn)
        Lchain = np.asarray(dem.Lchain, dtype=float).reshape(M, C)
        STchain = np.asarray(dem.STchain, dtype=float).reshape(M, C)
        Vchain = np.asarray(dem.Vchain, dtype=float).reshape(M, C)
        alpha = dem.alpha
        Nchain = np.asarray(dem.Nchain, dtype=float).ravel()
        SCVchain = np.asarray(dem.SCVchain, dtype=float).reshape(M, C)
        if np.any(np.isinf(Nchain)):
            raise ValueError('SJN scheduling is supported by SolverMVA only in closed models, '
                             'the open case has no population recursion.')

        sched = self.sched
        rows = []
        infrows = []
        sjnrows = []
        for i in range(M):
            sv = sched.get(i) if isinstance(sched, dict) else (sched[i] if i < len(sched) else None)
            code = sv.value if hasattr(sv, 'value') else sv
            if code == SchedStrategy.EXT.value:
                continue
            # the scheduling strategy alone selects delay against queue, nservers never does
            nsrv = float(sn.nservers[i]) if sn.nservers is not None else 1.0
            if code == SchedStrategy.INF.value:
                infrows.append(i)
                continue
            if code == SchedStrategy.SJF.value:
                if nsrv != 1:
                    raise ValueError('SJN scheduling at station %d requires a single server, the '
                                     'response time equation is a single-server one.' % (i + 1))
                sjnrows.append(len(rows))
            elif code in (SchedStrategy.PS.value, SchedStrategy.FCFS.value,
                          SchedStrategy.SIRO.value, SchedStrategy.LCFSPR.value):
                if nsrv != 1:
                    raise ValueError('station %d has %s servers, the SJN analyzer solves the '
                                     'remaining stations with the single-server MVA equation.'
                                     % (i + 1, nsrv))
            else:
                raise ValueError('The SJN analyzer does not support %s scheduling at the other '
                                 'stations.' % str(sv))
            rows.append(i)

        L = STchain[rows, :] * Vchain[rows, :]
        V = Vchain[rows, :]
        scv = np.ones((len(rows), C))
        for j in sjnrows:
            i = rows[j]
            for r in range(C):
                v = SCVchain[i, r]
                if np.isfinite(v) and v > 0:
                    scv[j, r] = v
        Z = np.zeros(C)
        for i in infrows:
            Z = Z + STchain[i, :] * Vchain[i, :]

        opt = SjnOptions()
        if getattr(self.options, 'iter_tol', None):
            opt.tol = float(self.options.iter_tol)
        # native options name it max_iter; OptionsDict-style options use iter_max
        _im = getattr(self.options, 'max_iter', None) or getattr(self.options, 'iter_max', None)
        if _im:
            opt.iter_max = int(_im)
        cfg = getattr(self.options, 'config', None)

        def _cfg(key):
            # options.config is a plain dict in native Python, an OptionsDict elsewhere
            if isinstance(cfg, dict):
                return cfg.get(key)
            return getattr(cfg, key, None) if cfg is not None else None

        for key, attr in (('sjn_ns', 'ns'), ('sjn_lfactor', 'lfactor'), ('sjn_umax', 'umax')):
            val = _cfg(key)
            if val is not None:
                setattr(opt, attr, type(getattr(opt, attr))(val))
        # SJN applies within a class and the classes are then non-preemptively prioritised;
        # without distinct priorities the jobs of every class are compared by size directly
        prio = np.asarray(getattr(sn, 'classprio', []), dtype=float).ravel()
        if sn.nchains == sn.nclasses and prio.size == C and np.unique(prio).size == C:
            opt.prio = prio.astype(int)

        latticemax = _cfg('sjn_lattice_max')
        latticemax = 1e5 if latticemax is None else float(latticemax)
        method = str(getattr(self.options, 'method', 'default')).lower()
        if method in ('amva', 'bs', 'sjn.amva'):
            uselattice = False
        elif method in ('exact', 'mva', 'sjn.mva'):
            uselattice = True
        else:
            uselattice = float(np.prod(Nchain + 1)) <= latticemax
        if uselattice:
            try:
                Xchain, Qrows, Urows, _Crows, _prof, it = pfqn_mvasjn(L, Nchain, Z, scv, sjnrows, V, opt)
                actualmethod = 'sjn.mva'
            except SjnStarvationError:
                if method != 'default':
                    raise
                Xchain, Qrows, Urows, _Crows, _prof, it = pfqn_amvasjn(L, Nchain, Z, scv, sjnrows, V, opt)
                actualmethod = 'sjn.amva'
        else:
            Xchain, Qrows, Urows, _Crows, _prof, it = pfqn_amvasjn(L, Nchain, Z, scv, sjnrows, V, opt)
            actualmethod = 'sjn.amva'

        Qchain = np.zeros((M, C)); Uchain = np.zeros((M, C))
        Rchain = np.zeros((M, C)); Tchain = np.zeros((M, C))
        Qchain[rows, :] = Qrows
        Uchain[rows, :] = Urows
        for i in range(M):
            Tchain[i, :] = Xchain * Vchain[i, :]
        for i in infrows:
            Qchain[i, :] = Tchain[i, :] * STchain[i, :]
            Uchain[i, :] = Qchain[i, :]
        with np.errstate(divide='ignore', invalid='ignore'):
            Rchain = np.where(Tchain > 0, Qchain / Tchain, 0.0)

        Xchain = np.asarray(Xchain, dtype=float).ravel().copy()
        Xchain[~np.isfinite(Xchain)] = 0.0
        Qchain[~np.isfinite(Qchain)] = 0.0
        Uchain[~np.isfinite(Uchain)] = 0.0
        Rchain[~np.isfinite(Rchain)] = 0.0
        # an empty chain carries no jobs, so every one of its metrics is zero
        zero = (Nchain == 0)
        if np.any(zero):
            Xchain[zero] = 0.0
            Qchain[:, zero] = 0.0
            Uchain[:, zero] = 0.0
            Rchain[:, zero] = 0.0
            Tchain[:, zero] = 0.0

        # MATLAB passes [] here and lets the deaggregation rebuild Q and U from Rchain and alpha
        res = sn_deaggregate_chain_results(sn, Lchain, None, STchain, Vchain, alpha,
                                           None, None, Rchain, Tchain, None, Xchain)
        self._lastiter = it
        self._result = {
            'QN': res.Q, 'UN': res.U, 'RN': res.R, 'TN': res.T,
            'AN': res.T.copy(), 'XN': np.asarray(res.X, dtype=float).ravel(), 'WN': res.R.copy(),
            'CN': res.C, 'runtime': _t.time() - t0, 'method': actualmethod,
        }
        if getattr(self.options, 'verbose', False):
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner("MVA analysis [method: %s; type: approximate, deterministic; lang: python] "
                  "completed in %.6fs." % (actualmethod, _t.time() - t0))
        return self

    def _run_marie(self):
        import time as _t
        t0 = _t.time()
        from ...api.pfqn.marie import pfqn_marie
        if self.nclasses != 1:
            return self._run_marie_multi(t0)
        inf = self._marie_inf_mask()
        rates = np.asarray(self.rates[:, 0], dtype=float)
        dem = np.asarray(self.demands[:, 0], dtype=float)
        with np.errstate(divide='ignore', invalid='ignore'):
            V = dem * rates
        D = dem[~inf]
        Z = float(np.sum(dem[inf]))
        N = int(round(float(self.njobs[0])))
        scv_full = (np.asarray(self._sn.scv)[:, 0]
                    if getattr(self._sn, 'scv', None) is not None
                    else np.ones(self.nstations))
        scv = np.asarray(scv_full, dtype=float)[~inf]
        if D.shape[0] == 0:
            # Nothing to isolate: with every station an infinite server the
            # aggregation-decomposition degenerates to the exact delay solution
            # X = N/Z, and pfqn_marie would be handed a zero-row demand matrix.
            Xchain = (N / Z) if Z > 0 else 0.0
            Qm = np.zeros(0)
            Um = np.zeros(0)
        else:
            X, Qm, Um, Cm, it, mu = pfqn_marie(D, N, Z, scv)
            Xchain = float(np.asarray(X).flatten()[0])
            Qm = np.asarray(Qm, dtype=float).flatten()
            Um = np.asarray(Um, dtype=float).flatten()
        M = self.nstations
        QN = np.zeros((M, 1)); UN = np.zeros((M, 1))
        RN = np.zeros((M, 1)); TN = np.zeros((M, 1))
        qj = 0
        with np.errstate(divide='ignore', invalid='ignore'):
            for i in range(M):
                TN[i, 0] = Xchain * V[i]
                if inf[i]:
                    RN[i, 0] = 1.0 / rates[i]
                    QN[i, 0] = TN[i, 0] * RN[i, 0]
                    UN[i, 0] = QN[i, 0]
                else:
                    QN[i, 0] = Qm[qj]
                    UN[i, 0] = Um[qj]
                    RN[i, 0] = QN[i, 0] / TN[i, 0] if TN[i, 0] > 0 else 0.0
                    qj += 1
        self._result = {
            'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN,
            'AN': TN.copy(), 'XN': np.array([Xchain]), 'WN': RN.copy(),
            'runtime': _t.time() - t0, 'method': 'marie',
        }
        if getattr(self.options, 'verbose', False):
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner("MVA analysis [method: marie; type: approximate, deterministic; lang: python] completed in "
                  "%.6fs." % (_t.time() - t0))
        return self

    def _run_mapqn(self):
        import time as _t
        t0 = _t.time()
        from ...api.solvers.mva.mapqn import solver_mva_mapqn_analyzer
        ret = solver_mva_mapqn_analyzer(self._sn, self.options)
        self._result = {
            'QN': ret.QN, 'UN': ret.UN, 'RN': ret.RN, 'TN': ret.TN,
            'AN': ret.AN, 'XN': ret.XN, 'WN': ret.WN, 'CN': ret.CN,
            'runtime': _t.time() - t0, 'method': 'amva.mapqn', 'iter': ret.iter,
        }
        if getattr(self.options, 'verbose', False):
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner("MVA analysis [method: amva.mapqn; type: approximate, deterministic; lang: python] "
                                "completed in %.6fs." % (_t.time() - t0))
        return self

    def _run_marie_multi(self, t0):
        import time as _t
        from ...api.pfqn.marie import pfqn_marie
        R = self.nclasses
        M = self.nstations
        inf = self._marie_inf_mask()
        rates = np.asarray(self.rates, dtype=float).reshape(M, R)
        dem = np.asarray(self.demands, dtype=float).reshape(M, R)
        with np.errstate(divide='ignore', invalid='ignore'):
            V = dem * rates                       # visit-ratio proxy per class
        D = dem[~inf, :]
        Z = np.zeros(R)
        if np.any(inf):
            Z = np.nansum(dem[inf, :], axis=0)
        N = np.round(np.asarray(self.njobs, dtype=float).ravel()).astype(int)
        scv_full = (np.asarray(self._sn.scv, dtype=float).reshape(M, R)
                    if getattr(self._sn, 'scv', None) is not None
                    else np.ones((M, R)))
        scv = scv_full[~inf, :]
        if D.shape[0] == 0:
            # Nothing to isolate: see the single-class arm above.
            X = np.zeros(R)
            for r in range(R):
                if Z[r] > 0:
                    X[r] = N[r] / Z[r]
            Qm = np.zeros((0, R))
            Um = np.zeros((0, R))
        else:
            X, Qm, Um, Cm, it, mu = pfqn_marie(D, N, Z, scv)
            X = np.asarray(X, dtype=float).ravel()
            Qm = np.asarray(Qm, dtype=float).reshape(-1, R)
            Um = np.asarray(Um, dtype=float).reshape(-1, R)
        QN = np.zeros((M, R)); UN = np.zeros((M, R))
        RN = np.zeros((M, R)); TN = np.zeros((M, R))
        qj = 0
        with np.errstate(divide='ignore', invalid='ignore'):
            for i in range(M):
                for r in range(R):
                    TN[i, r] = X[r] * V[i, r]
                if inf[i]:
                    for r in range(R):
                        RN[i, r] = 1.0 / rates[i, r] if rates[i, r] > 0 else 0.0
                        QN[i, r] = TN[i, r] * RN[i, r]
                        UN[i, r] = QN[i, r]
                else:
                    for r in range(R):
                        QN[i, r] = Qm[qj, r]
                        UN[i, r] = Um[qj, r]
                        RN[i, r] = (QN[i, r] / TN[i, r]
                                    if TN[i, r] > 0 else 0.0)
                    qj += 1
        self._result = {
            'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN,
            'AN': TN.copy(), 'XN': X.copy(), 'WN': RN.copy(),
            'runtime': _t.time() - t0, 'method': 'marie',
        }
        if getattr(self.options, 'verbose', False):
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner("MVA analysis [method: marie; type: approximate, deterministic; lang: python] completed in "
                  "%.6fs." % (_t.time() - t0))
        return self

    def _schmidt_arm_inputs(self):
        """(N, fcfs_rows) exactly as the schmidt / schmidt-ext / ab arm of
        runAnalyzer hands them to the kernel: the population vector it recurs on,
        and which of the queueing-station rows is served FCFS.

        It exists so that `supportsModelMethod` and the arm itself ask
        `mva_supports_schmidt_ext` the SAME question about the SAME numbers. The
        population is the CHAIN one under class switching, because that is the
        conserved vector the arm aggregates to; with one class per chain the two
        carry the same numbers. The delay row the arm stacks on top is left out
        on purpose: an INF row is never an FCFS station and the correction is
        never formed at one.
        """
        L, queue_indices = self._get_queueing_demands()
        njobs = self.njobs
        if self._sn is not None and int(getattr(self._sn, 'nchains', 0)) < self.nclasses:
            from ...api.sn import sn_get_demands_chain
            njobs = np.asarray(sn_get_demands_chain(self._sn).Nchain, dtype=float).flatten()
        fcfs = []
        for q_idx in queue_indices:
            _sched = (self.sched[q_idx]
                      if (self.sched is not None and q_idx in self.sched) else None)
            name = getattr(_sched, 'name', None) or str(_sched)
            fcfs.append('INF' not in name and 'PS' not in name)
        return njobs, fcfs

    def _amva_multiserver_rule(self, mi):
        """(max finite server count over the AMVA queueing stations, config rule).

        MATLAB `solver_amva.m` keys its whole multiserver treatment on these two
        values, and both branches below need them: the product-form arm applies
        Seidmann's transform under 'default'/'seidmann' (:111-117), and the
        linearizer family is handed to `solver_amvald` under
        'default'/'softmin'/'seidmann'/'suri' (:243-253). Ignoring them does not
        make the answer approximate, it solves a DIFFERENT model -- the one where
        every station has a single server.
        """
        ns = np.asarray(mi, dtype=float).ravel()
        ns = ns[np.isfinite(ns)]
        max_servers = int(np.max(ns)) if ns.size > 0 else 1
        rule = 'default'
        cfg = getattr(self.options, 'config', None)
        if isinstance(cfg, dict):
            rule = cfg.get('multiserver') or 'default'
        elif cfg is not None:
            rule = getattr(cfg, 'multiserver', None) or 'default'
        return max_servers, str(rule).lower()

    def _amva_softmin_multiserver(self, mi):
        """True for a multiserver model under the 'softmin' rule.

        MATLAB solver_amva.m:118 returns solver_amvald for that rule BEFORE its
        per-method switch, so it applies to every method, not only the linearizer
        family. Seidmann's transform is not applied there -- softmin asks for the
        load-dependent rate min(n,m) itself, which only amvald carries.
        """
        max_servers, rule = self._amva_multiserver_rule(mi)
        return max_servers > 1 and rule == 'softmin'

    def _lin_family_needs_amvald(self, mi):
        """True when lin/gflin/egflin must go to solver_amvald for multiserver.

        Complements `_amva_needs_amvald`, which covers only load- and
        class-dependent scaling. The linearizer family carries no server-count
        argument, so a multiserver model routed there loses `nservers` in exactly
        the way an LD model loses its scaling.
        """
        from ...api.sn import sn_has_product_form_not_het_fcfs
        max_servers, rule = self._amva_multiserver_rule(mi)
        if max_servers <= 1:
            return False
        if rule in ('default', 'softmin', 'seidmann', 'suri'):
            return True
        # MATLAB never reaches its lin arm for a het-FCFS model: the non-product-form
        # tail (solver_amva.m:397-401) takes it to solver_amvald whatever the rule is
        return self._sn is not None and not sn_has_product_form_not_het_fcfs(self._sn)

    @staticmethod
    def _amva_seidmann(L, Z, mi):
        """Seidmann's multiserver transform of (L, Z), as MATLAB solver_amva.m:111-117.

        A station with m servers is replaced by a single-server station of demand
        L/m plus a pure delay of L(m-1)/m folded into the think time, which is
        what makes the single-server linearizer family and pfqn_bs applicable to a
        multiserver model at all. Z is charged from the ORIGINAL demands, so the
        two updates cannot be reordered.
        """
        L0 = np.atleast_2d(np.asarray(L, dtype=float))
        Lms = L0.copy()
        Zms = np.array(Z, dtype=float, copy=True).ravel()
        ns = np.asarray(mi, dtype=float).ravel()
        for j in range(L0.shape[0]):
            m = ns[j] if j < ns.size else 1.0
            if not np.isfinite(m) or m <= 1:
                continue
            Lms[j, :] = L0[j, :] / m
            Zms += L0[j, :] * (m - 1.0) / m
        return Lms, Zms

    @staticmethod
    def _amva_seidmann_unapply(QN, RN, TN, queue_indices, L, mi, X):
        """Give each multiserver station back the population Seidmann folded away.

        The transform charges L(m-1)/m to the think time, so the solved queue
        length at station j counts only the jobs waiting for the one modelled
        server; the jobs in service at the other m-1 are sitting in the delay
        term. They belong to the station, so they are moved back here and the
        delay row keeps only the ORIGINAL think time -- charging both is what
        made sum(Q) exceed N. Mirrors MATLAB solver_amva.m.
        """
        L0 = np.atleast_2d(np.asarray(L, dtype=float))
        ns = np.asarray(mi, dtype=float).ravel()
        Xv = np.asarray(X, dtype=float).ravel()
        for idx, q_idx in enumerate(queue_indices):
            m = ns[idx] if idx < ns.size else 1.0
            if not np.isfinite(m) or m <= 1:
                continue
            QN[q_idx, :] = QN[q_idx, :] + L0[idx, :] * (m - 1.0) / m * Xv
            nz = TN[q_idx, :] > 0
            RN[q_idx, nz] = QN[q_idx, nz] / TN[q_idx, nz]

    def _warn_if_not_converged(self, method):
        """Report an AMVA fixed point that did not meet its tolerance.

        The authoritative signal is `_lastconverged`; the count is consulted only
        when no flag was reported, because on the load-dependent route the counter
        aggregates the nested sweeps and saturates the budget by construction on a
        solve whose outer residual is exactly zero. Shared by the native path and
        by the lang='java'/'cpp' delegations, which read both off the CLI payload:
        a warning raised in one lang and not the other is worse than none, since
        the silence is read as convergence.
        """
        if self._lastconverged is False or (
                self._lastconverged is None
                and self._lastiter and self._lastiterbudget
                and self._lastiter >= self._lastiterbudget):
            from ...api.io.logging import line_warning_always
            line_warning_always(
                'solver_mva_analyzer',
                "AMVA method '%s' did not meet the convergence tolerance %g after %d "
                "iterations; the returned metrics may not be converged. Try another method "
                "(e.g. 'qd' or 'bs'), raise options.iter_max, or loosen options.iter_tol."
                % (method, getattr(self.options, 'iter_tol', float('nan')), self._lastiter))

    def _adopt_delegated_convergence(self, container):
        """Carry the delegated solve's iteration count and convergence flag onto
        this solver, then apply the same warning the native path applies."""
        self._lastiter = getattr(container, 'iter', None)
        self._lastconverged = getattr(container, 'converged', None)
        _im = getattr(self.options, 'iter_max', None)
        self._lastiterbudget = int(_im) if _im else None
        self._warn_if_not_converged(getattr(container, 'method', None)
                                    or getattr(self.options, 'method', 'default'))

    def _interlock_matrix(self):
        """Interlock matrix of Franks (1999), Eq. (4.7) as SolverLN left it in the options.

        CLASS-indexed, so that a later refreshChains cannot leave it stale; the handler that
        is about to run aggregates it to chains against the struct it solves. None for every
        model but the layers of SolverLN.
        """
        cfg = getattr(self.options, 'config', None)
        if cfg is None:
            return None
        if isinstance(cfg, dict):
            return cfg.get('interlock', None)
        return getattr(cfg, 'interlock', None)

    def _apply_interlock(self, amvald_options):
        """Carry the same matrix into an AMVA-LD run, aggregated to its chain basis."""
        IL = self._interlock_matrix()
        if IL is None or np.size(IL) == 0:
            return
        from ...api.sn import sn_interlock_chain
        from ...api.solvers.mva.amvald import AmvaldOptions
        ILchain = sn_interlock_chain(self._sn, IL)
        if ILchain is None:
            return
        if getattr(amvald_options, 'config', None) is None:
            amvald_options.config = AmvaldOptions.Config()
        amvald_options.config.interlock_chain = ILchain

    def runAnalyzer(self):
        """Run the MVA analysis."""
        # MODEL TRANSFORMATION, opt-in through options.config['transform']. The
        # strategy rewrites the model into subproblems, TransformSolveMixin
        # solves each with an instance of THIS solver and maps the metrics back,
        # so a transformation written once serves MVA as well as CTMC. Mirrors
        # the branch MATLAB puts in the shared runAnalyzerPreamble.
        if self.maybe_transform():
            return self._result
        # A fresh analysis invalidates any prior unstable-utilization cap.
        self._unstable_util_capped = False
        # last AMVA iteration count, published as result['iter']; None for non-iterative (exact) paths, which is not an error.
        self._lastiter = None
        # iteration budget the handler used (not always options.iter_max, e.g. hardcoded linearizer maxiter=1000); None if unknown, convergence check skipped.
        self._lastiterbudget = None
        # authoritative convergence flag overrides the count-vs-budget heuristic; see _kb/06-solver-catalog.md MVA AMVA convergence flag vs iteration count.
        self._lastconverged = None
        # Closed models with shortest-job-next stations; see _kb/06-solver-catalog.md.
        _sjn = self._sjn_station_mask()
        if np.any(_sjn):
            if np.any(np.isinf(np.asarray(self.njobs, dtype=float))):
                # without the rejection the generic AMVA path would silently solve the station as
                # if it were size-blind, which is not what SJF means
                raise ValueError(
                    'SolverMVA supports shortest-job-next (SJF) scheduling only in closed models, '
                    'the conditional waiting time equation being a population recursion. Use '
                    'SolverLDES, or SolverMVA with SRPT or PSJF for the preemptive size-based '
                    'open queue.')
            return self._run_sjn()
        # Marie aggregation-decomposition for closed FCFS Coxian service; see _kb/06-solver-catalog.md MVA method='marie' section.
        _m0 = str(getattr(self.options, 'method', 'default')).lower()
        if _m0 in ('marie', 'amva.marie'):
            return self._run_marie()
        # Horizontal-cut MVA for one exponential delay and one FCFS MAP queue
        # (api.mapqn.mapqn_amva); see _kb/06-solver-catalog.md MVA method='amva.mapqn'.
        if _m0 in ('amva.mapqn', 'mapqn'):
            return self._run_mapqn()
        # Bound methods moved to SolverBA (mirrors MATLAB/JAR).
        _bfam = _m0.split('.')[0]
        if _bfam in ('aba', 'bjb', 'pb', 'gb', 'sb', 'mwba', 'pbh', 'pbk',
                     'bjbk', 'cbh', 'ssd', 'cub', 'mbjb', 'sib', 'scb', 'ldbcmp',
                     'looping'):
            raise ValueError(
                "Method '%s' is a bound method served by SolverBA; use "
                "SolverBA(model, method='%s'). Bound methods were moved out of "
                "SolverMVA." % (_m0, _m0))
        # lang='java' delegates to the canonical JAR, populating the native result container; imported lazily so a JVM-free install never touches this path.
        if getattr(self.options, 'lang', 'python') == 'java':
            # native OI/PAS gate applied before lang='java' delegation so both langs raise the same exception rather than an opaque JAR RuntimeError.
            self._resolve_oi_path(str(getattr(self.options, 'method', 'default')).lower())
            from ..jar_dispatch import populate_java_result
            self._adopt_delegated_convergence(populate_java_result(self))
            return self
        # lang='cpp' delegates to the C++ multiprecision port (line-cli) over the
        # same subprocess+JSON transport as lang='java'. Imported lazily so an
        # install without the binary never touches this path.
        #
        # THE ONLY AUTOMATIC FALLBACK IS AN ABSENT BINARY. line-cli is not built
        # by `pip install`, is platform-specific, and on an arch with no build
        # there is nothing to run -- degrading to native Python there is an
        # environment adaptation, and it warns so the reported lang and the
        # engine that ran cannot silently disagree. A construct the C++ analyzer
        # REFUSES propagates instead: the two ports do not refuse the same set,
        # and answering anyway would report a python number under lang='cpp',
        # which is the one thing this option exists to rule out.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            self._resolve_oi_path(str(getattr(self.options, 'method', 'default')).lower())
            from ..cpp_dispatch import LineCliNotAvailable, populate_cpp_result
            try:
                self._adopt_delegated_convergence(populate_cpp_result(self))
                return self
            except LineCliNotAvailable as e:
                line_warning("SolverMVA",
                             "lang='cpp' requested but the C++ solver is unavailable (%s); "
                             "falling back to lang='python'." % e)

        start_time = time.time()
        line_debug("MVA: using lang=python", options=self.options)

        if self._sn is not None and getattr(self._sn, 'immfeed', None) is not None and np.any(self._sn.immfeed):
            line_warning("SolverMVA", "SolverMVA does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.")

        # RQNA extends MVA feature set with the MAP family for explicit method='rqna' and 'default' auto-selecting RQNA on a bursty single-class open network.
        _method = str(getattr(self.options, 'method', 'default')).lower()
        _use_rqna_feats = (_method == 'rqna')
        if not _use_rqna_feats and _method == 'default' and self._sn is not None:
            from ...api.sn import sn_has_bursty_arrival
            _use_rqna_feats = (self._sn.nclasses == 1
                               and np.all(np.isinf(self._sn.njobs))
                               and sn_has_bursty_arrival(self._sn))

        # reject features outside MVA featset (finite capacity, FCR, JSQ/RROBIN, ...) not silently give unconstrained product-form; mirrors runAnalyzerChecks.
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_used_lang_features'):
            # method-aware feature gate: resolveMethod maps 'default' to 'rqna' only for a bursty single-class open network, so MAP/MMPP only on that path.
            self.runAnalyzerChecks(self.options)

        # MAP/MMPP2 carry autocorrelation a renewal/product-form MVA cannot represent, so rejected (BMAP exempt, batch); mirrors MATLAB/JAR runAnalyzerChecks.
        if (self._sn is not None and getattr(self._sn, 'procid', None) is not None
                and not _use_rqna_feats):
            from ...constants import ProcessType as _PT
            _mva_reject = {_PT.MAP: 'MAP', _PT.MMPP2: 'MMPP2', _PT.GAMMA: 'Gamma'}
            for _v in np.asarray(self._sn.procid, dtype=object).ravel():
                _hit = None
                for _pt, _nm in _mva_reject.items():
                    if _v == _pt:
                        _hit = _nm
                        break
                if _hit is not None:
                    raise RuntimeError(
                        "SolverMVA does not support the %s process used by this "
                        "model (not in the MVA feature set). Use SolverMAM for "
                        "MAP/Gamma service, or SolverCTMC/SolverSSA. This matches "
                        "MATLAB SolverMVA." % _hit)

        # QNA (two-moment decomposition) is reachable only by explicit request;
        # solver_mva_analyzer routes it the same way in MATLAB and the JAR.
        if _method == 'qna':
            from ...api.solvers.mva.analyzers import solver_qna
            from ...lang.base import NodeType
            ret = solver_qna(self._sn, self.options)
            QN = np.asarray(ret.Q, dtype=np.float64)
            UN = np.asarray(ret.U, dtype=np.float64)
            RN = np.asarray(ret.R, dtype=np.float64)
            TN = np.asarray(ret.T, dtype=np.float64)
            M_, K_ = QN.shape
            AN = TN.copy()
            rates_mat = self._sn.rates if self._sn.rates is not None else np.zeros((M_, K_))
            WN = np.zeros((M_, K_))
            for i in range(M_):
                nd = int(self._sn.stationToNode[i])
                if self._sn.nodetype[nd] == NodeType.SOURCE:
                    AN[i, :] = 0.0
                for r in range(K_):
                    if rates_mat[i, r] > 0 and RN[i, r] > 0:
                        WN[i, r] = max(0.0, RN[i, r] - 1.0 / rates_mat[i, r])
            self._result = {
                'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN,
                'XN': np.asarray(ret.X, dtype=np.float64).ravel(), 'WN': WN,
                'CN': np.sum(RN, axis=0),
                'runtime': getattr(ret, 'runtime', 0.0),
                'method': 'qna', 'iter': getattr(ret, 'it', 1),
            }
            return self._result

        # RQT (robust queueing theory) is reachable only by explicit request;
        # solver_mva_analyzer routes it the same way in MATLAB and the JAR.
        if _method == 'rqt':
            from ...api.solvers.mva.analyzers import solver_rqt
            from ...lang.base import NodeType
            ret = solver_rqt(self._sn, self.options)
            QN = np.asarray(ret.Q, dtype=np.float64)
            UN = np.asarray(ret.U, dtype=np.float64)
            RN = np.asarray(ret.R, dtype=np.float64)
            TN = np.asarray(ret.T, dtype=np.float64)
            M_, K_ = QN.shape
            AN = TN.copy()
            rates_mat = self._sn.rates if self._sn.rates is not None else np.zeros((M_, K_))
            WN = np.zeros((M_, K_))
            for i in range(M_):
                nd = int(self._sn.stationToNode[i])
                if self._sn.nodetype[nd] == NodeType.SOURCE:
                    AN[i, :] = 0.0
                for r in range(K_):
                    if rates_mat[i, r] > 0 and RN[i, r] > 0:
                        WN[i, r] = max(0.0, RN[i, r] - 1.0 / rates_mat[i, r])
            self._result = {
                'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN,
                'XN': np.asarray(ret.X, dtype=np.float64).ravel(), 'WN': WN,
                'CN': np.sum(RN, axis=0),
                'runtime': getattr(ret, 'runtime', 0.0),
                'method': 'rqt', 'iter': getattr(ret, 'it', 1),
            }
            return self._result

        # RQNA auto-selected for multi-queue bursty open nets; single-queue Source-Queue-Sink routes to qsys path (gm1/gig1), matching MATLAB/JAR dispatch.
        from ...lang.base import NodeType as _NT_rqna
        _rqna_dispatch = (_method == 'rqna')
        if not _rqna_dispatch and _use_rqna_feats and self._sn is not None:
            _n_nonsource = 0
            for _i in range(self._sn.nstations):
                _nd = int(self._sn.stationToNode[_i])
                if self._sn.nodetype[_nd] != _NT_rqna.SOURCE:
                    _n_nonsource += 1
            _rqna_dispatch = (_n_nonsource > 1)
        if _rqna_dispatch:
            from ...api.solvers.mva.analyzers import solver_rqna
            from ...lang.base import NodeType
            ret = solver_rqna(self._sn, self.options)
            QN = np.asarray(ret.Q, dtype=np.float64)
            UN = np.asarray(ret.U, dtype=np.float64)
            RN = np.asarray(ret.R, dtype=np.float64)
            TN = np.asarray(ret.T, dtype=np.float64)
            M_, K_ = QN.shape
            # arrival rates: throughput at each station, zero at the source(s)
            AN = TN.copy()
            rates_mat = self._sn.rates if self._sn.rates is not None else np.zeros((M_, K_))
            WN = np.zeros((M_, K_))
            for i in range(M_):
                nd = int(self._sn.stationToNode[i])
                if self._sn.nodetype[nd] == NodeType.SOURCE:
                    AN[i, :] = 0.0
                for r in range(K_):
                    if rates_mat[i, r] > 0 and RN[i, r] > 0:
                        WN[i, r] = max(0.0, RN[i, r] - 1.0 / rates_mat[i, r])
            XN = np.asarray(ret.X, dtype=np.float64).ravel()
            self._result = {
                'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN,
                'XN': XN, 'WN': WN, 'CN': np.sum(RN, axis=0),
                'runtime': getattr(ret, 'runtime', 0.0),
                'method': 'rqna', 'iter': 1,
            }
            return self._result

        noi_idx, oi_exact = self._resolve_oi_path(_method)

        if oi_exact:
            line_debug("Order-independent closed network, routing to solver_mva_oi_analyzer", options=self.options)
            from .solver_mva_oi_analyzer import SolverMVAOIAnalyzer
            analyzer = SolverMVAOIAnalyzer(self._sn, self.options)
            result = analyzer.analyze()
            QN = result['QN']
            UN = result['UN']
            RN = result['RN']
            TN = result['TN']
            XN = result['XN']
            M_, K_ = QN.shape
            AN = TN.copy()
            WN = RN.copy()
            self._result = {
                'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN,
                'XN': XN, 'WN': WN, 'CN': np.sum(RN, axis=0),
                'runtime': result.get('runtime', 0.0),
                'method': result.get('method', 'oi'), 'iter': result.get('iter', 1),
            }
            return self._result

        # Check for delayed-hit cache with a retrieval system (FPI algorithms + latency)
        from ...api.retrieval.analyzers import has_retrieval_cache
        if self._sn is not None and has_retrieval_cache(self._sn):
            line_debug("Delayed-hit retrieval cache, routing to retrieval_analyzer", options=self.options)
            result = self._run_retrieval_analysis()
            if result is not None:
                return result

        # Check for cache-only networks and handle specially
        if self._is_cache_only_network():
            line_debug("Non-reentrant cache (Source-Cache-Sink), routing to cache_analyzer", options=self.options)
            result = self._run_cache_analysis()
            if result is not None:
                return result

        # Check for cache networks with class switching (hit/miss classes)
        # Skip if already in cache QN analysis to avoid recursion
        if not getattr(self, '_skip_cache_qn', False) and self._has_cache_with_class_switching():
            line_debug("Integrated caching-queueing network, routing to cacheqn_analyzer", options=self.options)
            result = self._run_cache_qn_analysis()
            if result is not None:
                return result

        # Size-based M/G/1 (SRPT, PSJF, FB, LRPT, SETF) takes the exact
        # Wierman-Harchol-Balter response times; see _kb/06-solver-catalog.md
        if self._sizebased_sched() is not None:
            line_debug("Size-based scheduling detected, routing to "
                       "qsys_sizebased_analyzer", options=self.options)
            result = self._run_sizebased_analysis()
            if result is not None:
                return result

        # Check for polling systems and handle specially
        if self._is_polling_system():
            line_debug("Multiclass open polling system, routing to polling_analyzer", options=self.options)
            result = self._run_polling_analysis()
            if result is not None:
                return result
            # Fall through to standard analysis if polling handling failed

        # Check for fork-join networks and handle specially
        # Skip if _skip_fork_join is set (to avoid recursion during H-T transformation)
        if self._has_fork_join() and not getattr(self, '_skip_fork_join', False):
            line_debug("Fork-join network detected, routing to fork_join_analysis", options=self.options)
            result = self._run_fork_join_analysis()
            if result is not None:
                return result
            # Fall through to standard analysis if fork-join handling failed

        from ...api.pfqn import (
            pfqn_mva, pfqn_aql, pfqn_linearizer, pfqn_gflinearizer,
            pfqn_egflinearizer, pfqn_mvald, pfqn_mvams, pfqn_bs, pfqn_sqni,
            pfqn_schmidt, pfqn_schmidt_ext, pfqn_ab_amva,
            pfqn_linearizermx,
        )
        # pfqn_qdlin/qli/fli and api.pfqn.bounds intentionally not imported: the qd
        # family reaches the solver through solver_amvald and bounds through SolverBA.
        # pfqn_qdlin is the array-level TWIN of the 'qdlin' method, not its implementation.

        method = self.method.lower()

        # A single-class open station with reneging resolves to its abandonment
        # method here, not only in the feature gate: the gate and the analyzer
        # must agree, or the gate would admit the model and the analyzer would
        # then solve it as if nobody ever abandoned.
        if method == 'default':
            _ab = self._resolve_abandonment_method()
            if _ab is not None:
                method = _ab

        # Normalize AMVA method aliases
        method = method.replace('amva.', '')

        # The closed-population AMVA family (Bard-Schweitzer, SQNI, Tay, SCAT,
        # AQL, QSA, Bard LCP, Chow SA, Hsieh-Lam PAM, clustering, Improved
        # Linearizer, Akyildiz-Bolch, Schmidt) lives ONLY in the closed
        # product-form branch below. The same predicate the report gates on
        # decides here, so a name the report offers is a name that runs, and a
        # name it withholds errors by name instead of being answered with a
        # table of zeros or with the qd-family numbers under someone else's.
        from ...api.solvers.mva.handler import (
            mva_supports_closed_population, mva_is_closed_population_method)
        _cp_ok, _cp_reason = mva_supports_closed_population(self._sn, method)
        if not _cp_ok:
            raise ValueError(_cp_reason)

        # Get parameters
        L, queue_indices = self._get_queueing_demands()
        if len(queue_indices) == 0 and mva_is_closed_population_method(method):
            # Degenerate network: with no queueing station there is no
            # arrival-instant queue to correct, so every AMVA approximation
            # coincides with the exact delay solution Q_ir = X_r D_ir and the
            # name a caller passed selects nothing. MATLAB solver_amva.m takes
            # the same exit (sn_has_homogeneous_scheduling INF) ahead of its
            # method switch; without it the family was handed a zero-row demand
            # matrix and died inside pfqn_bs, or reported zeros.
            method = 'lin'
        N = self.njobs.copy()
        Z = self._get_think_times()
        mi = self.nservers[queue_indices] if len(queue_indices) > 0 else np.ones(1)

        M = L.shape[0]  # Number of queueing stations
        R = self.nclasses  # Number of classes

        # THE CLOSED-POPULATION AMVA FAMILY RECURS ON A CONSERVED POPULATION.
        # Under class switching a job CHANGES CLASS as it moves, so no per-class
        # population is conserved and the vector these kernels need is the CHAIN
        # one. That is why MATLAB solver_amva.m and the C++ port build their whole
        # product-form branch out of sn_get_product_form_chain_params and
        # deaggregate at the end (solver_amva.m:159 and :426). Handed class
        # populations instead, the family solved a DIFFERENT network: on a
        # two-class Delay+PS switching model whose exact answer is [1.4118,
        # 0.5882], all sixteen names returned the same [2, 0] -- every job parked
        # at the delay -- and sixteen different approximations agreeing bit for
        # bit is the signature of that, not of accuracy.
        # Chain and class coincide exactly when every chain holds one class, so
        # the substitution is made only where it changes the answer: a model
        # without class switching keeps the code path, and the numbers, it had.
        cp_chain = None
        if (mva_is_closed_population_method(method) and self._sn is not None
                and int(getattr(self._sn, 'nchains', R)) < R):
            from ...api.sn import sn_get_demands_chain
            cp_chain = sn_get_demands_chain(self._sn)
            cp_delays = self._get_delay_stations()
            R = int(self._sn.nchains)
            L = cp_chain.Lchain[queue_indices, :]
            N = np.asarray(cp_chain.Nchain, dtype=float).flatten().copy()
            Z = np.zeros(R)
            for _d in cp_delays:
                Z = Z + np.asarray(cp_chain.Lchain[_d, :], dtype=float)

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

        # open feedback nets need chain-level MVA for visit ratios; explicit qsys names (mm1/gig1/gm1/...) route to qsys dispatch, not generic MVA branch.
        _qsys_dispatch_methods = {
            'mm1', 'mmk', 'mg1', 'mgi1', 'gm1', 'gig1', 'gim1',
            'gig1.kingman', 'gigk', 'gigk.kingman_approx',
            'gig1.gelenbe', 'gig1.heyman', 'gig1.kimura',
            'gig1.allen', 'gig1.kobayashi', 'gig1.klb', 'gig1.marchal',
            # Whitt family; must be listed here as well as in qsys_methods
            # below, or the generic MVA branch swallows the model and the qsys
            # dispatch never runs.
            'erlanga', 'mgisrgi', 'gigk.diffusion',
            'gigk.whitt', 'qed', 'gig1.extremal',
        }
        _is_qsys_dispatch = (self.network_type == 'open' and self.nstations == 2
                             and self.nclasses == 1 and method in _qsys_dispatch_methods)

        # Defaults for flags set inside the generic MVA branch but read by the
        # common post-processing below (the qsys dispatch skips that branch)
        has_class_switching_early = False
        used_chain_deaggregation = False

        if self.network_type in ('closed', 'mixed', 'open') and not _is_qsys_dispatch:
            # Check if load-dependent MVA should be used (MATLAB: solver_mvald_analyzer)
            # Must check BEFORE converting 'default' method to 'exact'/'amva'
            use_ld_mva = False
            if self.lldscaling is not None:
                use_ld_mva = True
            # Check for class-dependent scaling (MATLAB: ~isempty(sn.cdscaling))
            if hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None:
                use_ld_mva = True
            # Joint-dependent (non-product-form eta_i) scaling routes the same way.
            if getattr(self._sn, 'jdscaling', None) is not None:
                use_ld_mva = True
            # Mixed networks require chain-level MVA (pfqn_mvaldmx) which handles open classes
            if self.network_type == 'mixed':
                use_ld_mva = True
            # Check for class-dependent scaling beta_{i,r}(n) (e.g., from FES aggregation)
            _cd = getattr(self._sn, 'cdscaling', None)
            if _cd is not None and len(_cd) > 0 and any(x is not None for x in _cd):
                use_ld_mva = True
            _jd = getattr(self._sn, 'jdscaling', None)
            if _jd is not None and len(_jd) > 0 and any(x is not None for x in _jd):
                use_ld_mva = True

            # Check for class switching (multiple classes in same chain) BEFORE method selection
            # This is needed because 'amva' doesn't handle class switching properly
            has_class_switching_early = False
            if hasattr(self._sn, 'nchains') and self._sn.nchains > 0:
                chains = self._get_chains()
                for chain in chains:
                    if len(chain) > 1:
                        has_class_switching_early = True
                        break

            # Track if chain-level deaggregation was used (skip delay recompute if so)
            used_chain_deaggregation = False

            # Check product form early for class switching decision
            # Use relaxed check (not_het_fcfs) matching MATLAB solver_amva.m line 90
            from ...api.sn import sn_has_product_form_not_het_fcfs, sn_has_product_form as _sn_pf_scvblind
            has_product_form_early = sn_has_product_form_not_het_fcfs(self._sn) if self._sn is not None else True
            # SCV-blind product-form predicate: MATLAB's default/exact dispatch
            # keys on this one, treating non-exponential FCFS as product form
            has_product_form_scvblind = _sn_pf_scvblind(self._sn) if self._sn is not None else True
            line_debug("Product-form check: hasProductForm=%s (exact method requested)", has_product_form_early, options=self.options)

            # An open single-class single-queue system (Source-Queue-Sink) takes
            # the exact qsys closed forms -- M/M/1, M/M/k, M/G/1, G/M/1 -- as
            # MATLAB's mvaDispatch branch 4 and the C++ port both do.
            #
            # THIS GATE USED TO EXCLUDE THE PLAIN M/M/1, on the stated ground
            # that "egflin is already exact when ca=cs=1". IT IS NOT. egflin is a
            # linearizer whose iteration stops on a tolerance: on lambda=0.5,
            # mu=1 it returns QLen 0.999999245 where rho/(1-rho) is exactly 1,
            # short by 7.6e-7. The condition was therefore excluding from the
            # exact path precisely the model the exact path is cheapest on, and
            # it made this the one model where lang='python' could not reproduce
            # lang='cpp' or MATLAB to solver tolerance -- which is how it was
            # found. Measured, not reasoned: the divergence is pinned in
            # python/tests/test_mva_lang_cpp.py.
            #
            # The multiserver arm below already sends M/M/k here for the same
            # reason (the Seidmann term underestimates), so this restores the
            # single-server case to the company it belongs in.
            if (method == 'default' and self.network_type == 'open' and M == 1 and R == 1
                    and len(source_indices) > 0 and len(queue_indices) > 0
                    and self._sn.procid is not None):
                method = 'exact'

            # Handle 'default' method with MATLAB-compatible heuristic
            if method == 'default':
                # For mixed networks without actual load-dependent service (no lldscaling),
                # use linearizer which properly handles both open and closed classes
                is_mixed_only = self.network_type == 'mixed' and self.lldscaling is None

                # For open networks with product form and single servers, MATLAB uses egflin
                # (solver_amva.m lines 49-55: if max(nservers)==1, method='egflin')
                is_open_product_form = (self.network_type == 'open' and
                                        has_product_form_early and
                                        self.lldscaling is None)
                max_servers = 1
                if is_open_product_form and self.nservers is not None:
                    finite_servers = self.nservers[np.isfinite(self.nservers)]
                    if len(finite_servers) > 0:
                        max_servers = int(np.max(finite_servers))

                from ...api.solvers.mva.analyzers import _is_bas_model
                if _is_bas_model(self._sn):
                    # Closed single-chain network with Blocking-After-Service finite buffers
                    method = 'sqd'
                elif is_mixed_only:
                    # exact BCMP mixed MVA (pfqn_mvamx) avoids egflin open/closed double-counting open-class interference; egflin is fallback for non-PF or multiserver.
                    _mixed_max_srv = 1
                    if self.nservers is not None:
                        _fs = self.nservers[np.isfinite(self.nservers)]
                        if len(_fs) > 0:
                            _mixed_max_srv = int(np.max(_fs))
                    # the STRICT product-form test is used here (not the het-FCFS-relaxed one): exact pfqn_mvamx rejects heterogeneous-rate FCFS, which stays on egflin.
                    from ...api.sn import sn_has_product_form as _sn_has_pf
                    _mixed_pf = _sn_has_pf(self._sn) if self._sn is not None else False
                    if _mixed_pf and _mixed_max_srv == 1:
                        method = 'exact'
                    elif _mixed_max_srv > 1:
                        # multiserver models route through the amva branch (solver_amvald with 'lin'); 'egflin' is MATLAB's single-server-only choice.
                        method = 'amva'
                    else:
                        method = 'egflin'
                elif (self.network_type == 'open' and M == 1 and R >= 2
                      and len(queue_indices) > 0 and len(source_indices) > 0
                      and self._sn is not None and self._sn.sched
                      and self._dps_exact_applicable(queue_indices[0], source_indices[0], mi)):
                    # single open M/M/1-DPS routes to the exact qsys_mm1_dps: the AMVA-DPS cross-term correction violates equal-rate conservation.
                    method = 'exact'
                elif is_open_product_form and max_servers == 1:
                    # Open network with product form and single servers - use egflin
                    method = 'egflin'
                elif is_open_product_form and max_servers > 1:
                    # open M/M/k routes to exact MVA (Erlang-C via pfqn_mvaldms/qsys_mmk): AMVA Seidmann term underestimates QLen/RespT by treating it as a scaled M/M/1.
                    method = 'exact'
                elif use_ld_mva:
                    # small closed product-form LD models use exact pfqn_mvaldmx recursion, mirroring MATLAB's default-to-exact upgrade; else fall back to approx LD-AMVA.
                    from ...api.sn import sn_has_product_form as _sn_has_pf
                    _cd = getattr(self._sn, 'cdscaling', None) if self._sn is not None else None
                    _jd = getattr(self._sn, 'jdscaling', None) if self._sn is not None else None
                    _nchains = self._sn.nchains if (self._sn is not None and hasattr(self._sn, 'nchains')) else R
                    if (_cd is None and _jd is None and np.all(np.isfinite(N)) and _nchains <= 4
                            and np.sum(N) <= 20
                            and (self._sn is None or _sn_has_pf(self._sn))
                            and np.all(N == np.floor(N))):
                        method = 'exact'
                    else:
                        method = 'amva'
                elif (self.network_type == 'open' and M == 1 and R >= 2
                      and len(queue_indices) > 0 and len(source_indices) > 0
                      and self._sn is not None and self._sn.sched
                      and self._hol_cobham_applicable(queue_indices[0], source_indices[0], mi)):
                    # single open HOL M/G/1 routes to exact qsys_mg1_prio (Cobham), not the AMVA preemptive shadow-server approximation.
                    method = 'exact'
                elif not has_product_form_scvblind:
                    # non-PF nets (priorities, fork-join, sd-routing, heterog FCFS) use AMVA; MATLAB SCV-blind PF test still sends small closed non-exp-FCFS to exact MVA.
                    method = 'amva'
                else:
                    # Match MATLAB's solver_mva_analyzer.m logic for non-LD models:
                    # Use exact MVA if: nchains <= 4 && sum(njobs) <= 20 && product_form && no fractional populations
                    nchains = self._sn.nchains if hasattr(self._sn, 'nchains') and self._sn is not None else R
                    # IMPORTANT: For open networks, total_jobs should be infinity (matches MATLAB's sum(sn.njobs))
                    # This ensures open networks use 'amva' path which handles saturation correctly
                    if np.any(np.isinf(N)):
                        total_jobs = np.inf
                    else:
                        total_jobs = int(np.sum(N[np.isfinite(N)]))
                    from ...api.sn import sn_has_product_form
                    has_product_form = sn_has_product_form(self._sn) if self._sn is not None else True
                    has_fractional = np.any(N != np.floor(N))

                    if nchains <= 4 and total_jobs <= 20 and has_product_form and not has_fractional:
                        method = 'exact'
                    else:
                        method = 'amva'


            # INF-scheduled Queue stations (not Delay, which is naturally INF) need load-dependent treatment too.
            from ...lang.base import SchedStrategy
            from ...api.sn.network_struct import NodeType
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

            # Check for class switching (multiple classes in same chain)
            has_class_switching = False
            if hasattr(self._sn, 'nchains') and self._sn.nchains > 0:
                chains = self._get_chains()
                for chain in chains:
                    if len(chain) > 1:
                        has_class_switching = True
                        break

            # plain class switching without real load dependence uses the chain-based pfqn_mvams (more accurate); pfqn_mvaldmx is reserved for genuine LD cases.
            needs_mvaldmx = use_ld_mva and (self.lldscaling is not None or has_inf_server)

            # MATLAB solver_amva.m:81-91 resolves {'default','amva'} to qd / egflin / lin
            # unconditionally, as do the JAR Solver_amva.java:149-159 and the cpp
            # solver_mva.h:1159-1170. The plain 'amva' arm answers a different model.
            if method == 'amva':
                _Nvec = np.asarray(N, dtype=float).flatten()
                _chains = self._get_chains() if (self._sn is not None and self._sn.nchains > 0) else None
                if _chains:
                    _Nchain = np.array([float(np.sum(_Nvec[list(ch)])) for ch in _chains])
                else:
                    _Nchain = _Nvec
                _srv = mi[np.isfinite(mi)] if mi is not None else np.array([1.0])
                _maxsrv = int(np.max(_srv)) if _srv.size > 0 else 1
                if np.sum(_Nchain[np.isfinite(_Nchain)]) <= 2 or np.any(_Nchain < 1):
                    method = 'qd'
                elif _maxsrv == 1:
                    method = 'egflin'
                else:
                    method = 'lin'

            if method == 'sqd':
                # Closed single-chain Blocking-After-Service network (finite buffers)
                line_debug("Blocking-After-Service network, routing to solver_sqd", options=self.options)
                from ...api.solvers.mva.analyzers import solver_sqd
                result = solver_sqd(self._sn, None)
                _rq = np.asarray(result.Q) if result.Q is not None else np.zeros((0, 0))
                if _rq.size == 0 or np.all(np.isnan(_rq)):
                    # multichain / unsupported: solver_sqd already warned, return empty
                    self._result = {'QN': result.Q, 'UN': result.U, 'RN': result.R,
                                    'TN': result.T, 'AN': result.T,
                                    'XN': np.asarray(result.X).flatten(), 'WN': result.R,
                                    'runtime': time.time() - start_time, 'method': 'sqd'}
                    return self
                QN = result.Q if result.Q is not None else np.zeros((self.nstations, R))
                UN = result.U if result.U is not None else np.zeros((self.nstations, R))
                RN = result.R if result.R is not None else np.zeros((self.nstations, R))
                TN = result.T if result.T is not None else np.zeros((self.nstations, R))
                XN = result.X.flatten() if result.X is not None else np.zeros(R)
                AN = TN.copy()
                used_chain_deaggregation = True

            elif method in ('sum', 'esum'):
                # summation method (SUM/ESUM) for closed, closing method for open/mixed; mirrors MATLAB solver_mva_sum.m.
                from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results
                from ...api.sum import sum_closed, sum_closing
                from ...lang.base import SchedStrategy

                chain_result = sn_get_demands_chain(self._sn)
                Lchain = chain_result.Lchain
                STchain = chain_result.STchain
                Vchain = chain_result.Vchain
                alpha = chain_result.alpha
                Nchain = chain_result.Nchain.flatten()
                SCVchain = chain_result.SCVchain
                refstatchain = chain_result.refstatchain.flatten().astype(int)
                Cc = self._sn.nchains
                M_full = Lchain.shape[0]

                def _sched_val(i):
                    sched_i = self.sched.get(i) if isinstance(self.sched, dict) else (
                        self.sched[i] if (self.sched is not None and i < len(self.sched)) else None)
                    return sched_i.value if hasattr(sched_i, 'value') else sched_i

                # station rows passed to the summation method (all but the source)
                rows = []
                mi = []
                scv_sensitive = []
                for ist in range(M_full):
                    sval = _sched_val(ist)
                    if sval == SchedStrategy.EXT.value:
                        continue  # external world handled by lambda/sum_closing
                    elif sval == SchedStrategy.INF.value:
                        rows.append(ist)
                        mi.append(np.inf)
                        scv_sensitive.append(False)
                    elif sval in (SchedStrategy.PS.value, SchedStrategy.LCFSPR.value,
                                  SchedStrategy.FCFS.value, SchedStrategy.SIRO.value):
                        rows.append(ist)
                        mi.append(float(self.nservers[ist]))
                        scv_sensitive.append(sval in (SchedStrategy.FCFS.value,
                                                      SchedStrategy.SIRO.value))
                    else:
                        raise ValueError('The summation method does not support this scheduling strategy.')

                rows = np.asarray(rows, dtype=int)
                L = STchain[rows, :] * Vchain[rows, :]
                scv = np.ones((len(rows), Cc))
                for j, ist in enumerate(rows):
                    if scv_sensitive[j]:
                        for c in range(Cc):
                            if np.isfinite(SCVchain[ist, c]) and SCVchain[ist, c] > 0:
                                scv[j, c] = SCVchain[ist, c]

                Zc = np.zeros(Cc)
                if not np.any(np.isinf(Nchain)):
                    Xchain, Qrows, Urows, _, iters = sum_closed(
                        L, Nchain, Zc, np.asarray(mi), scv,
                        self.options.iter_tol, int(self.options.max_iter))
                    self._lastiter = iters
                    self._lastiterbudget = None
                else:
                    lambda_chain = np.zeros(Cc)
                    scva = np.ones(Cc)
                    for c in range(Cc):
                        if np.isinf(Nchain[c]):
                            refstat = int(refstatchain[c])
                            lambda_chain[c] = 1.0 / STchain[refstat, c]
                            if np.isfinite(SCVchain[refstat, c]) and SCVchain[refstat, c] > 0:
                                scva[c] = SCVchain[refstat, c]  # interarrival SCV at the source
                    Xchain, Qrows, Urows, _, _, iters = sum_closing(
                        lambda_chain, scva, L, np.asarray(mi), scv, Nchain, Zc,
                        5000, self.options.iter_tol, int(self.options.max_iter))
                    self._lastiter = iters
                    self._lastiterbudget = None

                Qchain = np.zeros((M_full, Cc))
                Uchain = np.zeros((M_full, Cc))
                Qchain[rows, :] = Qrows
                Uchain[rows, :] = Urows
                Xchain = np.where(np.isfinite(Xchain), Xchain, 0.0)
                Tchain = np.outer(np.ones(M_full), Xchain) * Vchain
                Rchain = np.zeros((M_full, Cc))
                for c in range(Cc):
                    if Nchain[c] == 0:
                        Xchain[c] = 0.0
                        Qchain[:, c] = 0.0
                        Uchain[:, c] = 0.0
                        Tchain[:, c] = 0.0
                        continue
                    for i in range(M_full):
                        if Tchain[i, c] > 0:
                            Rchain[i, c] = Qchain[i, c] / Tchain[i, c]

                deagg = sn_deaggregate_chain_results(
                    self._sn, Lchain, None, STchain, Vchain, alpha,
                    None, None, Rchain, Tchain, None, Xchain.reshape(1, -1)
                )
                QN = deagg.Q
                UN = deagg.U
                RN = deagg.R
                TN = deagg.T
                XN = deagg.X.flatten()
                AN = TN.copy()
                used_chain_deaggregation = True

            elif needs_mvaldmx and method in ['exact', 'mva']:
                line_debug("Load-dependent scaling detected (lldscaling=%s, cdscaling=%s), routing to mvald_analyzer",
                           self.lldscaling is not None, hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None,
                           options=self.options)
                # Use EXACT load-dependent MVA (MATLAB: solver_mvald -> pfqn_mvaldmx)
                # IMPORTANT: MATLAB's solver_mvald works at chain level, then disaggregates
                from ...api.pfqn import pfqn_mvaldmx
                from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

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
                M_full = Lchain.shape[0]  # Number of ALL stations
                S = self.nservers  # Server counts per station
                refstat_chain = np.asarray(refstatchain).flatten().astype(int)
                open_chains = [c for c in range(C) if np.isinf(Nchain[c])]

                # Chain arrival rates. The reference station of an open chain is its
                # Source and STchain holds one over the SUM of the class arrival rates
                # there, so reading the rate at chain level also covers a chain whose
                # classes arrive at several rates, or one carrying a class reached only
                # by a switch (which has no arrival process of its own).
                lambda_chain = np.zeros(C)
                for c in open_chains:
                    rst = int(refstat_chain[c])
                    if STchain[rst, c] > 0:
                        lambda_chain[c] = 1.0 / STchain[rst, c]

                if not open_chains:
                    # PURELY CLOSED: every station enters the recursion, an infinite
                    # server as the load-dependent rate mu(n)=n, which is exact because
                    # n cannot then exceed the closed population.
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
                    # pfqn_mvaldmx called with chain-level parameters and S=sn.nservers (all stations); mirrors MATLAB.
                    Xchain, Qchain, Uchain, _, lGN, Pc = pfqn_mvaldmx(
                        lambda_chain, Lchain, Nchain, np.zeros(C), mu_chain, S
                    )
                else:
                    # MIXED OR PURELY OPEN. Three kinds of row are not the same thing to
                    # pfqn_mvaldmx and have to be separated before it is called. This is
                    # the partition solver_ncld makes for the same recursion.
                    #  - THE SOURCE IS NOT A STATION. Its chain demand is the
                    #    interarrival time 1/lambda, so it carries offered load Lo=1
                    #    exactly and pfqn_ldmx_ec then forms 1/(1-Lo/mu)=inf.
                    #  - A DELAY IS AN INFINITE SERVER FOR THE OPEN CHAINS TOO. mu(n)=n
                    #    cut at the closed population declares it saturated at total_pop
                    #    jobs. It enters as chain think time instead and its queue length
                    #    is X*L, which is exact.
                    #  - A QUEUEING STATION KEEPS ITS WHOLE RATE ROW. pfqn_ldmx_ec reads
                    #    the limited-load-dependence level b off the row itself, so a row
                    #    cut at the closed population is read as a slower station, and
                    #    with no closed class at all it collapses to mu(1).
                    source_stations = sorted({int(refstat_chain[c]) for c in open_chains})
                    delay_stations = [i for i in range(M_full)
                                      if i < len(S) and np.isinf(S[i]) and i not in source_stations]
                    queue_stations = [i for i in range(M_full)
                                      if i not in source_stations and i not in delay_stations]
                    Zchain = np.zeros(C)
                    if delay_stations:
                        for c in range(C):
                            Zchain[c] = float(np.sum(Lchain[delay_stations, c]))
                    lld_width = self.lldscaling.shape[1] if self.lldscaling is not None else 0
                    ncol = max(1, total_pop)
                    for i in queue_stations:
                        # first column of the trailing constant run, the level b of pfqn_ldmx_ec
                        b = lld_width if (self.lldscaling is not None and i < self.lldscaling.shape[0]) else 0
                        while b > 1 and self.lldscaling[i, b - 2] == self.lldscaling[i, b - 1]:
                            b -= 1
                        ncol = max(ncol, b)
                    mu_chain = np.ones((len(queue_stations), ncol))
                    for qi, i in enumerate(queue_stations):
                        if lld_width > 0 and i < self.lldscaling.shape[0]:
                            avail = min(ncol, lld_width)
                            mu_chain[qi, :avail] = self.lldscaling[i, :avail]
                            mu_chain[qi, avail:] = self.lldscaling[i, lld_width - 1]  # saturated tail
                    Xchain, Qqueue, Uqueue, _, lGN, Pc = pfqn_mvaldmx(
                        lambda_chain, Lchain[queue_stations, :], Nchain, Zchain, mu_chain,
                        np.ones(len(queue_stations))
                    )
                    Qchain = np.zeros((M_full, C))
                    Uchain = np.zeros((M_full, C))
                    if queue_stations:
                        Qchain[queue_stations, :] = Qqueue
                        Uchain[queue_stations, :] = Uqueue
                    for i in delay_stations:
                        # infinite server: X*L for a closed chain, lambda*L for an open one
                        Qchain[i, :] = Lchain[i, :] * Xchain

                # Tchain(k,r)=Xchain(r)*Vchain(k,r), Rchain=Qchain./Tchain; mirrors MATLAB solver_mva.m.
                Tchain = np.outer(np.ones(M_full), Xchain) * Vchain
                Rchain = np.zeros((M_full, C))
                for c in range(C):
                    for i in range(M_full):
                        if Tchain[i, c] > 0:
                            Rchain[i, c] = Qchain[i, c] / Tchain[i, c]

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
                # busy-server fraction under load-dependent scaling: carried load over max(nservers, peak lldscaling), NC convention, not mvaldmx's P(busy) estimator.
                if self.lldscaling is not None:
                    for ist in range(min(UN.shape[0], self.lldscaling.shape[0])):
                        if ist < len(S) and np.isfinite(S[ist]):
                            ceff = max(float(S[ist]), float(np.max(self.lldscaling[ist, :])))
                            for r in range(UN.shape[1]):
                                rate_ir = self.rates[ist, r] if (self.rates is not None and ist < self.rates.shape[0]) else 0.0
                                if np.isfinite(rate_ir) and rate_ir > 0:
                                    UN[ist, r] = TN[ist, r] / rate_ir / ceff
                used_chain_deaggregation = True
            elif method in ['exact', 'mva', 'mvac']:
                line_debug("Standard queueing network, routing to mva_analyzer (method=%s)", method, options=self.options)
                # chain-level aggregation and post-processing mirrors MATLAB solver_mva: chain demands, chain MVA, recompute Q/X from waiting times, disaggregate.

                # LCFS+LCFS-PR 2-station PF check runs before general PF test, ignoring LCFS scheduling; mirrors MATLAB solver_mva.m:22-44/JAR Solver_mva.kt:42-80.
                has_lcfs_network = False
                if self._sn is not None and self._sn.sched:
                    from ...lang.base import SchedStrategy as _SS
                    _lcfs_found = any(
                        _sched_val == _SS.LCFS or
                        (hasattr(_sched_val, 'value') and _sched_val.value == _SS.LCFS.value)
                        for _sched_val in self._sn.sched.values()
                    )
                    _lcfspr_found = any(
                        _sched_val == _SS.LCFSPR or
                        (hasattr(_sched_val, 'value') and _sched_val.value == _SS.LCFSPR.value)
                        for _sched_val in self._sn.sched.values()
                    )
                    has_lcfs_network = _lcfs_found and _lcfspr_found

                # amvald substitution gated on the same product-form predicate MATLAB's
                # exact path uses. METHOD 'mva' IS THE DELIBERATE APPROXIMATION: the
                # dispatch warns that the exact recursion runs outside its hypotheses
                # and promises an answer, so an explicit request keeps the recursion
                # where MATLAB's mvaDispatch keeps it.
                _pf_means = has_product_form_scvblind or method == 'mva'

                if method == 'mvac':
                    # Exact MVA by chain (MVAC, Conway et al. 1989): closed SSFR +
                    # IS product-form networks. Dispatch to the chain-based handler
                    # with the pfqn_mvac core; it validates closed/single-server.
                    from ...api.solvers.mva.handler import solver_mva as mva_handler
                    from ...api.solvers.mva.handler import SolverMVAOptions as MVAHandlerOptions

                    handler_options = MVAHandlerOptions(method='mvac', tol=1e-8)
                    result = mva_handler(self._sn, handler_options)

                    QN = result.Q if result.Q is not None else np.zeros((self.nstations, R))
                    UN = result.U if result.U is not None else np.zeros((self.nstations, R))
                    RN = result.R if result.R is not None else np.zeros((self.nstations, R))
                    TN = result.T if result.T is not None else np.zeros((self.nstations, R))
                    XN = result.X.flatten() if result.X is not None else np.zeros(R)
                    AN = TN.copy()

                elif has_lcfs_network:
                    # Dispatch directly to the handler which has LCFS detection and
                    # routes to _solver_mva_lcfsqn (specialized LCFS MVA)
                    from ...api.solvers.mva.handler import solver_mva as mva_handler
                    from ...api.solvers.mva.handler import SolverMVAOptions as MVAHandlerOptions

                    # 'mva' is forwarded verbatim: it is the deliberate approximation,
                    # and the handler's product-form guard exempts it by that name.
                    handler_options = MVAHandlerOptions(method=('mva' if method == 'mva' else 'exact'), tol=1e-8, interlock=self._interlock_matrix())
                    result = mva_handler(self._sn, handler_options)

                    QN = result.Q if result.Q is not None else np.zeros((self.nstations, R))
                    UN = result.U if result.U is not None else np.zeros((self.nstations, R))
                    RN = result.R if result.R is not None else np.zeros((self.nstations, R))
                    TN = result.T if result.T is not None else np.zeros((self.nstations, R))
                    XN = result.X.flatten() if result.X is not None else np.zeros(R)
                    AN = TN.copy()

                # single open HOL M/G/1 uses the exact Cobham formula; the AMVA preemptive shadow-server path underestimates both classes' waiting time.
                elif (self.network_type == 'open' and M == 1 and R >= 2
                      and len(queue_indices) > 0 and len(source_indices) > 0
                      and self._sn is not None and self._sn.sched
                      and self._hol_cobham_applicable(queue_indices[0], source_indices[0], mi)):
                    from ...api.qsys import qsys_mg1_prio
                    q_idx = queue_indices[0]
                    src_idx = source_indices[0]
                    sn_scv = self._sn.scv if self._sn.scv is not None else np.ones((self.nstations, R))
                    prios = np.asarray(self._sn.classprio).flatten() if self._sn.classprio is not None else np.zeros(R)
                    order = list(np.argsort(prios, kind='stable'))  # 0 = highest priority first
                    lam = np.array([float(self.rates[src_idx, r]) for r in order])
                    mus = np.array([1.0 / L[0, r] if L[0, r] > 0 else np.inf for r in order])
                    css = np.array([float(np.sqrt(sn_scv[q_idx, r]))
                                    if np.isfinite(sn_scv[q_idx, r]) and sn_scv[q_idx, r] > 0 else 1.0
                                    for r in order])
                    active = [i for i in range(len(order)) if lam[i] > 0 and np.isfinite(mus[i])]
                    W_ord, _ = qsys_mg1_prio(lam[active], mus[active], css[active])
                    for j, i in enumerate(active):
                        r = order[i]
                        RN[q_idx, r] = W_ord[j]
                        XN[r] = lam[i]
                        UN[q_idx, r] = lam[i] / mus[i]
                        TN[q_idx, r] = lam[i]
                        AN[q_idx, r] = lam[i]
                        QN[q_idx, r] = lam[i] * W_ord[j]
                        TN[src_idx, r] = lam[i]
                        AN[src_idx, r] = lam[i]

                # single open M/M/1-DPS uses the exact truncated multiclass CTMC (qsys_mm1_dps): the AMVA-DPS cross-term correction violates equal-rate conservation.
                elif (self.network_type == 'open' and M == 1 and R >= 2
                      and len(queue_indices) > 0 and len(source_indices) > 0
                      and self._sn is not None and self._sn.sched
                      and self._dps_exact_applicable(queue_indices[0], source_indices[0], mi)):
                    from ...api.qsys import qsys_mm1_dps
                    q_idx = queue_indices[0]
                    src_idx = source_indices[0]
                    wvec = np.asarray(self._sn.schedparam)[q_idx, :].astype(float) \
                        if getattr(self._sn, 'schedparam', None) is not None else np.ones(R)
                    lam = np.array([float(self.rates[src_idx, r]) for r in range(R)])
                    # raw service rates (not demand matrix) used here: AMVA folds the DPS weight into demand, an approx artifact corrupting the exact solver input.
                    mus = np.array([float(self.rates[q_idx, r]) if self.rates[q_idx, r] > 0 else np.inf
                                    for r in range(R)])
                    active = [r for r in range(R) if lam[r] > 0 and np.isfinite(mus[r])]
                    w_act = np.where(wvec[active] > 0, wvec[active], 1.0)
                    T_act, _rho = qsys_mm1_dps(lam[active], mus[active], w_act)
                    for j, r in enumerate(active):
                        RN[q_idx, r] = T_act[j]
                        XN[r] = lam[r]
                        UN[q_idx, r] = lam[r] / mus[r]
                        TN[q_idx, r] = lam[r]
                        AN[q_idx, r] = lam[r]
                        QN[q_idx, r] = lam[r] * T_act[j]
                        TN[src_idx, r] = lam[r]
                        AN[src_idx, r] = lam[r]

                # open single-class single-queue uses exact qsys formulas (M/M/1, M/M/k, M/G/1, G/M/1) regardless of product-form; mirrors solver_mva_qsys_analyzer.
                elif self.network_type == 'open' and M == 1 and R == 1:
                    from ...api.qsys import qsys_mm1, qsys_mmk, qsys_mg1, qsys_gg1
                    sn_scv = self._sn.scv if self._sn.scv is not None else np.ones((self.nstations, R))
                    src_idx = source_indices[0] if len(source_indices) > 0 else None
                    q_idx = queue_indices[0]
                    ca = float(np.sqrt(sn_scv[src_idx, 0])) if (src_idx is not None and np.isfinite(sn_scv[src_idx, 0]) and sn_scv[src_idx, 0] >= 0) else 1.0
                    cs = float(np.sqrt(sn_scv[q_idx, 0])) if (np.isfinite(sn_scv[q_idx, 0]) and sn_scv[q_idx, 0] >= 0) else 1.0
                    # The queue's visit ratio, which a feedback or re-entrant loop
                    # raises above one and which separates the per-visit quantities
                    # from the per-job ones. Reading mu off the DEMAND L=V*S instead
                    # of the service rate, and lambda off the source rate alone, is
                    # the same model only when Vq==1: it leaves QLen and Util right
                    # but reports the per-job residence time as RespT and the
                    # external arrival rate as the station throughput.
                    Vq = _qsys_queue_visits(self._sn, q_idx)
                    src_rate = float(self.rates[src_idx, 0]) if src_idx is not None else 1.0
                    lambda_r = src_rate * Vq
                    mu = float(self.rates[q_idx, 0])
                    nserv = mi[0] if len(mi) > 0 else 1.0
                    k = 1 if not np.isfinite(nserv) else int(nserv)

                    from ...constants import ProcessType as _PT
                    is_bmap = (src_idx is not None and self._sn.procid is not None
                               and self._sn.procid[src_idx, 0] == _PT.BMAP)
                    if is_bmap and self._sn.procid[q_idx, 0] == _PT.EXP and k == 1:
                        # BMAP arrivals + exponential service: exact MX/M/1 batch
                        # queue (matches MATLAB solver_mva_qsys_analyzer.m BMAP branch)
                        from ...api.qsys import qsys_mxm1
                        lambda_batch, E_X, E_X2 = _bmap_batch_moments(self._sn.proc[src_idx][0])
                        W_x, _, _, _ = qsys_mxm1(lambda_batch, mu, E_X, E_X2)
                        result = {'W': W_x}
                        lambda_r = lambda_batch * E_X  # effective job arrival rate
                        src_rate = lambda_r
                    elif ca == 1.0 and cs == 1.0 and k == 1:
                        result = qsys_mm1(lambda_r, mu)
                    elif ca == 1.0 and cs == 1.0 and k > 1:
                        result = qsys_mmk(lambda_r, mu, k)
                    elif ca == 1.0 and k == 1:
                        result = qsys_mg1(lambda_r, mu, cs)
                    elif cs == 1.0 and k == 1:
                        # exact PH/M/1 only when sn.proc holds the arrival law itself; a non-Markovian arrival gets an Erlang-n SCV fit, so falls to the exact LST sigma-root.
                        result = None
                        _src_is_markovian = (src_idx is not None
                                             and self._sn.procid is not None
                                             and _PT.isMarkovian(self._sn.procid[src_idx, 0]))
                        try:
                            from ...api.qsys import qsys_phm1 as _qsys_phm1
                            pie_p, D0p = _ph_from_proc(self._sn, src_idx) if _src_is_markovian else (None, None)
                            if pie_p is not None and D0p is not None:
                                res_ph = _qsys_phm1(pie_p, D0p, mu)
                                result = {
                                    'L': res_ph['mean_queue_length'],
                                    'Lq': res_ph['mean_waiting_queue'],
                                    'W': res_ph['mean_sojourn_time'],
                                    'Wq': res_ph['mean_waiting_time'],
                                    'rho': res_ph['utilization'],
                                }
                        except Exception:
                            result = None
                        if result is None:
                            # general exact G/M/1 sigma-root via the arrival's LST (sn.lst); falls to qsys_gg1 only absent an LST.
                            _Wlst = _gm1_lst_sojourn(self._sn, src_idx, mu)
                            if _Wlst is not None:
                                result = {'W': _Wlst}
                        if result is None:
                            _Wg, _ = qsys_gg1(lambda_r, mu, ca ** 2, 1.0)
                            result = {'W': _Wg}
                    elif k > 1:
                        # G/G/k approximation (matches MATLAB qsys 'gigk')
                        from ...api.qsys import qsys_gigk_approx
                        _Wg, _ = qsys_gigk_approx(lambda_r, mu, ca, cs, k)
                        result = {'W': _Wg}
                    else:
                        # G/G/1 KLB approximation (matches MATLAB default 'gig1.klb')
                        from ...api.qsys import qsys_gig1_approx_klb
                        _Wg, _ = qsys_gig1_approx_klb(lambda_r, mu, ca, cs)
                        result = {'W': _Wg}

                    # RespT/QLen are per-visit, the system throughput is the external
                    # arrival rate and the queue throughput the effective one; mirrors
                    # the tail of solver_mva_qsys_analyzer.m.
                    Rscalar = result['W']
                    RN[q_idx, 0] = Rscalar
                    XN[0] = src_rate
                    UN[q_idx, 0] = lambda_r / mu / k
                    TN[q_idx, 0] = lambda_r
                    AN[q_idx, 0] = lambda_r
                    QN[q_idx, 0] = lambda_r * Rscalar
                    if src_idx is not None:
                        TN[src_idx, 0] = src_rate
                        AN[src_idx, 0] = src_rate

                # Check for product form - if not, fall back to AMVA (MATLAB behavior)
                # Non-product-form open networks (e.g., heterogeneous FCFS) use solver_amvald
                elif not _pf_means:
                    # Use solver_amvald for non-product-form models
                    from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain
                    SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                    amvald_options = AmvaldOptions(method='default', iter_tol=self.options.iter_tol, iter_max=self.options.max_iter, init_sol=getattr(self.options, 'init_sol', None))
                    self._apply_interlock(amvald_options)
                    result = solver_amvald(
                        self._sn, Lchain, STchain, Vchain, alpha,
                        Nchain, SCVchain, refstatchain, amvald_options
                    )
                    self._lastiter = getattr(result, 'totiter', None)
                    self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                    self._lastconverged = getattr(result, 'converged', None)

                    # Uchain passed only under load/class-dependent scaling; otherwise deaggregation computes U=T*S; mirrors MATLAB solver_amvald.m:216-220.
                    has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                      np.any(self._sn.lldscaling != 0))
                    has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                      np.any(self._sn.cdscaling != 0))
                    Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                    deagg = sn_deaggregate_chain_results(
                        self._sn, Lchain, None, STchain, Vchain, alpha,
                        None, Uchain_for_deagg, result.R, result.T, None, result.X
                    )

                    QN = deagg.Q
                    UN = deagg.U
                    RN = deagg.R
                    TN = deagg.T
                    XN = deagg.X.flatten()
                    AN = TN.copy()
                    used_chain_deaggregation = True
                else:
                    from ...api.solvers.mva.handler import solver_mva as mva_handler
                    from ...api.solvers.mva.handler import SolverMVAOptions as MVAHandlerOptions

                    # 'mva' is forwarded verbatim: it is the deliberate approximation,
                    # and the handler's product-form guard exempts it by that name.
                    handler_options = MVAHandlerOptions(method=('mva' if method == 'mva' else 'exact'), tol=1e-8, interlock=self._interlock_matrix())
                    result = mva_handler(self._sn, handler_options)

                    # Copy results from handler
                    QN = result.Q if result.Q is not None else np.zeros((self.nstations, R))
                    UN = result.U if result.U is not None else np.zeros((self.nstations, R))
                    RN = result.R if result.R is not None else np.zeros((self.nstations, R))
                    TN = result.T if result.T is not None else np.zeros((self.nstations, R))
                    XN = result.X.flatten() if result.X is not None else np.zeros(R)
                    AN = TN.copy()

            elif method == 'amva':
                line_debug("AMVA method selected, checking multiserver/class-switching/product-form", options=self.options)
                # Check if there are multi-server queues (servers > 1, but not delay stations with inf servers)
                # Delay stations have inf servers but are NOT multiserver queues
                has_multiserver = np.any((mi > 1) & np.isfinite(mi)) if mi is not None else False

                # class-switching closed-model method choice (linearizermx+egflin / solver_amvald / solver_amvald+lin) mirrors MATLAB solver_amva.m.
                is_closed_network = not np.any(np.isinf(N))

                # Check load dependence
                has_load_dep = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None) or \
                               (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None)

                if has_class_switching_early and is_closed_network and has_product_form_early and not has_load_dep:
                    # MATLAB solver_amva.m lines 90-208: product-form (not-het-fcfs) path
                    # Uses sn_get_product_form_chain_params and pfqn_linearizermx
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results, sn_get_product_form_chain_params
                    from ...api.pfqn import pfqn_linearizermx

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()

                    nservers = self._sn.nservers
                    if nservers is None:
                        nservers = np.ones(self.nstations)
                    else:
                        nservers = nservers.flatten()

                    # Check max servers for queue stations
                    max_servers = 1
                    for q_idx in queue_indices:
                        if nservers[q_idx] < np.inf:
                            max_servers = max(max_servers, int(nservers[q_idx]))

                    C_chains = self._sn.nchains

                    if max_servers == 1:
                        # Single-server: use pfqn_linearizermx with 'egflin'
                        # (MATLAB solver_amva.m line 208)
                        pf_params = sn_get_product_form_chain_params(self._sn)
                        L_pf = pf_params.D  # Queue demands only (Mq x C)
                        N_pf = pf_params.N.flatten()  # Chain populations
                        Z_pf = pf_params.Z.flatten() if pf_params.Z.ndim > 1 else pf_params.Z  # Delay demands
                        S_pf = pf_params.S.flatten()  # Servers at queues
                        lambda_pf = pf_params.lambda_vec.flatten()

                        # sn_get_product_form_chain_params keys D and S on the
                        # NODE TYPE, so an INF-scheduled Queue is a row of D
                        # with S=inf; sched must be taken over that same set,
                        # as MATLAB solver_amva.m does with nodeToStation(queueIdx).
                        from ...api.sn.network_struct import NodeType as _NT
                        _ntv = self._sn.nodetype if isinstance(self._sn.nodetype, np.ndarray) else np.array(
                            [nt.value if hasattr(nt, 'value') else nt for nt in self._sn.nodetype])
                        pf_stations = [int(self._sn.nodeToStation[i])
                                       for i in np.where(_ntv == _NT.QUEUE.value)[0]]
                        sched_list = []
                        for q_idx in pf_stations:
                            sched_val = self.sched.get(q_idx, SchedStrategy.FCFS) if self.sched else SchedStrategy.FCFS
                            sched_list.append(sched_val)

                        # Handle all-delay case: no queue-type stations
                        if not pf_stations:
                            # pure-delay chain throughput X_c=N_c/D_c uses FULL chain demand, not Z_pf (misclassifies INF-scheduled station as queue); mirrors solver_amva.m:64.
                            Dchain_tot = np.sum(Lchain, axis=0).flatten()
                            Xchain_out = np.zeros(C_chains)
                            for c in range(C_chains):
                                if Dchain_tot[c] > 0 and N_pf[c] > 0 and not np.isinf(N_pf[c]):
                                    Xchain_out[c] = N_pf[c] / Dchain_tot[c]
                            Qchain_out = np.zeros((0, C_chains))
                            Uchain_out = np.zeros((0, C_chains))
                        else:
                            Qchain_out, Uchain_out, Wchain_out, Tchain_out, Cchain_out, Xchain_out, iters = pfqn_linearizermx(
                                lambda_pf, L_pf, N_pf, Z_pf, S_pf, sched_list,
                                tol=1e-4, maxiter=1000, method='egflin'
                            )
                            self._lastiter = iters
                            self._lastiterbudget = 1000

                        # Build full station chain results
                        Qchain = np.zeros((self.nstations, C_chains))
                        Uchain = np.zeros((self.nstations, C_chains))
                        Rchain = np.zeros((self.nstations, C_chains))
                        Tchain = np.zeros((self.nstations, C_chains))
                        Xchain = Xchain_out.reshape(1, -1) if Xchain_out.ndim == 1 else Xchain_out

                        # Map queue results back to station indices
                        for idx, q_idx in enumerate(pf_stations):
                            for c in range(C_chains):
                                Qchain[q_idx, c] = Qchain_out[idx, c] if Qchain_out.ndim > 1 else Qchain_out[idx]
                                Uchain[q_idx, c] = Uchain_out[idx, c] if Uchain_out.ndim > 1 else Uchain_out[idx]

                        # Compute delay station metrics
                        delay_indices = [i for i in range(self.nstations) if i not in pf_stations]
                        for d_idx in delay_indices:
                            for c in range(C_chains):
                                Qchain[d_idx, c] = Xchain[0, c] * STchain[d_idx, c] * Vchain[d_idx, c]
                                Uchain[d_idx, c] = Qchain[d_idx, c]

                        # Compute throughputs and response times
                        for c in range(C_chains):
                            for i in range(self.nstations):
                                if Vchain[i, c] > 0 and Xchain[0, c] > 0:
                                    Tchain[i, c] = Xchain[0, c] * Vchain[i, c]
                                    if Tchain[i, c] > 0:
                                        Rchain[i, c] = Qchain[i, c] / Tchain[i, c]

                        # Disaggregate chain results to class level
                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, None, Rchain, Tchain, None, Xchain
                        )

                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True
                    else:
                        # Multi-server: use solver_amvald (MATLAB solver_amva.m line 216)
                        from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions

                        SCVchain = np.ones((self._sn.nstations, self._sn.nchains))
                        refstatchain = chain_result.refstatchain

                        amvald_options = AmvaldOptions(
                            method='lin',
                            iter_tol=self.options.iter_tol,
                            iter_max=1000,
                            init_sol=getattr(self.options, 'init_sol', None)
                        )

                        self._apply_interlock(amvald_options)
                        result = solver_amvald(
                            self._sn, Lchain, STchain, Vchain, alpha,
                            Nchain, SCVchain, refstatchain, amvald_options
                        )
                        self._lastiter = getattr(result, 'totiter', None)
                        self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                        self._lastconverged = getattr(result, 'converged', None)

                        has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                          np.any(self._sn.lldscaling != 0))
                        has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                          np.any(self._sn.cdscaling != 0))
                        Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, Uchain_for_deagg, result.R, result.T, None, result.X
                        )

                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True

                elif has_class_switching_early and is_closed_network:
                    # Non-product-form closed model with class switching
                    # Uses solver_amvald with 'lin'
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results
                    from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain
                    SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                    amvald_options = AmvaldOptions(
                        method='lin',
                        iter_tol=self.options.iter_tol,
                        iter_max=1000,
                        init_sol=getattr(self.options, 'init_sol', None)
                    )

                    self._apply_interlock(amvald_options)
                    result = solver_amvald(
                        self._sn, Lchain, STchain, Vchain, alpha,
                        Nchain, SCVchain, refstatchain, amvald_options
                    )
                    self._lastiter = getattr(result, 'totiter', None)
                    self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                    self._lastconverged = getattr(result, 'converged', None)

                    has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                      np.any(self._sn.lldscaling != 0))
                    has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                      np.any(self._sn.cdscaling != 0))
                    Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                    deagg = sn_deaggregate_chain_results(
                        self._sn, Lchain, None, STchain, Vchain, alpha,
                        None, Uchain_for_deagg, result.R, result.T, None, result.X
                    )

                    QN = deagg.Q
                    UN = deagg.U
                    RN = deagg.R
                    TN = deagg.T
                    XN = deagg.X.flatten()
                    AN = TN.copy()
                    used_chain_deaggregation = True

                elif has_class_switching_early and not has_product_form_early and not np.any(np.isinf(N)):
                    # non-PF class-switching closed models: single-server uses linearizermx+egflin at chain level; multiserver uses solver_amvald; open nets fall through.
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain

                    # Check for single-server model (MATLAB: max(nservers)==1)
                    nservers_full = self._sn.nservers.flatten() if self._sn.nservers is not None else np.ones(self.nstations)
                    max_servers = 1
                    finite_servers = nservers_full[np.isfinite(nservers_full)]
                    if len(finite_servers) > 0:
                        max_servers = int(np.max(finite_servers))

                    C_chains = self._sn.nchains

                    if max_servers == 1:
                        # Single-server model: use pfqn_linearizermx with 'egflin' at chain level
                        # (matching MATLAB solver_amva.m lines 199-201)

                        # Build chain-level demands for queueing stations
                        # MATLAB: L = STchain .* Vchain (at chain level)
                        L_chain = STchain * Vchain  # M x C

                        # Extract queueing stations only
                        L_chain_queues = L_chain[queue_indices, :]  # Mq x C

                        # Think times from delay stations (sum of service times at delay stations)
                        delay_indices = [i for i in range(self.nstations) if i not in queue_indices]
                        Z_chain = np.zeros(C_chains)
                        for d_idx in delay_indices:
                            Z_chain += STchain[d_idx, :] * Vchain[d_idx, :]

                        # Get scheduling for queueing stations
                        sched_list = []
                        for q_idx in queue_indices:
                            sched_val = self.sched.get(q_idx, SchedStrategy.FCFS) if self.sched else SchedStrategy.FCFS
                            sched_list.append(sched_val)

                        # Server counts for queueing stations
                        S_queue = nservers_full[queue_indices]

                        # Call pfqn_linearizermx at chain level
                        lambda_chain = np.zeros(C_chains)
                        Qchain_out, Uchain_out, Wchain_out, Tchain_out, Cchain_out, Xchain_out, iters = pfqn_linearizermx(
                            lambda_chain, L_chain_queues, Nchain, Z_chain, S_queue, sched_list,
                            tol=1e-4, maxiter=1000, method='egflin'
                        )
                        self._lastiter = iters
                        self._lastiterbudget = 1000

                        # Build full station chain results
                        Qchain = np.zeros((self.nstations, C_chains))
                        Uchain = np.zeros((self.nstations, C_chains))
                        Rchain = np.zeros((self.nstations, C_chains))
                        Tchain = np.zeros((self.nstations, C_chains))
                        Xchain = Xchain_out.reshape(1, -1) if Xchain_out.ndim == 1 else Xchain_out

                        # Map queue results back to station indices
                        for idx, q_idx in enumerate(queue_indices):
                            for c in range(C_chains):
                                Qchain[q_idx, c] = Qchain_out[idx, c] if Qchain_out.ndim > 1 else Qchain_out[idx]
                                Uchain[q_idx, c] = Uchain_out[idx, c] if Uchain_out.ndim > 1 else Uchain_out[idx]

                        # Compute delay station metrics
                        for d_idx in delay_indices:
                            for c in range(C_chains):
                                Qchain[d_idx, c] = Xchain[0, c] * STchain[d_idx, c] * Vchain[d_idx, c]
                                Uchain[d_idx, c] = Qchain[d_idx, c]

                        # Compute throughputs and response times
                        for c in range(C_chains):
                            for i in range(self.nstations):
                                if Vchain[i, c] > 0 and Xchain[0, c] > 0:
                                    Tchain[i, c] = Xchain[0, c] * Vchain[i, c]
                                    if Tchain[i, c] > 0:
                                        Rchain[i, c] = Qchain[i, c] / Tchain[i, c]

                        # Disaggregate chain results to class level
                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, None, Rchain, Tchain, None, Xchain
                        )

                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True
                    else:
                        # Multiserver model: use solver_amvald
                        from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions

                        SCVchain = np.ones((self._sn.nstations, self._sn.nchains))
                        amvald_options = AmvaldOptions(method='default', iter_tol=self.options.iter_tol, iter_max=self.options.max_iter, init_sol=getattr(self.options, 'init_sol', None))
                        self._apply_interlock(amvald_options)
                        result = solver_amvald(
                            self._sn, Lchain, STchain, Vchain, alpha,
                            Nchain, SCVchain, refstatchain, amvald_options
                        )
                        self._lastiter = getattr(result, 'totiter', None)
                        self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                        self._lastconverged = getattr(result, 'converged', None)

                        # Disaggregate chain results to class level
                        has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                          np.any(self._sn.lldscaling != 0))
                        has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                          np.any(self._sn.cdscaling != 0))
                        Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, Uchain_for_deagg, result.R, result.T, None, result.X
                        )

                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True

                elif np.any(np.isinf(N)):
                    # Open or mixed network - use chain-based AMVA (solver_amvald)
                    # This handles open classes correctly with visit ratios
                    from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain
                    SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                    amvald_options = AmvaldOptions(method='default', iter_tol=self.options.iter_tol, iter_max=self.options.max_iter, init_sol=getattr(self.options, 'init_sol', None))
                    self._apply_interlock(amvald_options)
                    result = solver_amvald(
                        self._sn, Lchain, STchain, Vchain, alpha,
                        Nchain, SCVchain, refstatchain, amvald_options
                    )
                    self._lastiter = getattr(result, 'totiter', None)
                    self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                    self._lastconverged = getattr(result, 'converged', None)

                    # Uchain passed only under load/class-dependent scaling; otherwise deaggregation computes U=T*S; mirrors MATLAB solver_amvald.m:216-220.
                    has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                      np.any(self._sn.lldscaling != 0))
                    has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                      np.any(self._sn.cdscaling != 0))
                    Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                    deagg = sn_deaggregate_chain_results(
                        self._sn, Lchain, None, STchain, Vchain, alpha,
                        None, Uchain_for_deagg, result.R, result.T, None, result.X
                    )

                    QN = deagg.Q
                    UN = deagg.U
                    RN = deagg.R
                    TN = deagg.T
                    XN = deagg.X.flatten()
                    AN = TN.copy()
                    used_chain_deaggregation = True
                elif has_multiserver:
                    # multi-server closed nets use solver_amvald, handling chain aggregation, Seidmann transform and deaggregation; mirrors MATLAB solver_amva.m:209-217.
                    from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain
                    SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                    amvald_options = AmvaldOptions(
                        method='lin',
                        iter_tol=self.options.iter_tol,
                        iter_max=1000,
                        init_sol=getattr(self.options, 'init_sol', None)
                    )
                    if self.options.config and 'multiserver' in self.options.config:
                        amvald_options.config = AmvaldOptions.Config(
                            multiserver=self.options.config['multiserver']
                        )

                    self._apply_interlock(amvald_options)
                    result = solver_amvald(
                        self._sn, Lchain, STchain, Vchain, alpha,
                        Nchain, SCVchain, refstatchain, amvald_options
                    )
                    self._lastiter = getattr(result, 'totiter', None)
                    self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                    self._lastconverged = getattr(result, 'converged', None)

                    # Disaggregate chain results to class level
                    has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                      np.any(self._sn.lldscaling != 0))
                    has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                      np.any(self._sn.cdscaling != 0))
                    Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                    deagg = sn_deaggregate_chain_results(
                        self._sn, Lchain, None, STchain, Vchain, alpha,
                        None, Uchain_for_deagg, result.R, result.T, None, result.X
                    )

                    QN = deagg.Q
                    UN = deagg.U
                    RN = deagg.R
                    TN = deagg.T
                    XN = deagg.X.flatten()
                    AN = TN.copy()
                    used_chain_deaggregation = True
                elif not has_product_form_early:
                    # Non-product-form closed network (e.g., HOL priority)
                    # Use solver_amvald which handles priority scheduling correctly
                    from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                    from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                    chain_result = sn_get_demands_chain(self._sn)
                    Lchain = chain_result.Lchain
                    STchain = chain_result.STchain
                    Vchain = chain_result.Vchain
                    alpha = chain_result.alpha
                    Nchain = chain_result.Nchain.flatten()
                    refstatchain = chain_result.refstatchain
                    SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                    amvald_options = AmvaldOptions(method='default', iter_tol=self.options.iter_tol, iter_max=self.options.max_iter, init_sol=getattr(self.options, 'init_sol', None))
                    self._apply_interlock(amvald_options)
                    result = solver_amvald(
                        self._sn, Lchain, STchain, Vchain, alpha,
                        Nchain, SCVchain, refstatchain, amvald_options
                    )
                    self._lastiter = getattr(result, 'totiter', None)
                    self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                    self._lastconverged = getattr(result, 'converged', None)

                    # Uchain passed only under load/class-dependent scaling; otherwise deaggregation computes U=T*S; mirrors MATLAB solver_amvald.m:216-220.
                    has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                      np.any(self._sn.lldscaling != 0))
                    has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                      np.any(self._sn.cdscaling != 0))
                    Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                    deagg = sn_deaggregate_chain_results(
                        self._sn, Lchain, None, STchain, Vchain, alpha,
                        None, Uchain_for_deagg, result.R, result.T, None, result.X
                    )

                    QN = deagg.Q
                    UN = deagg.U
                    RN = deagg.R
                    TN = deagg.T
                    XN = deagg.X.flatten()
                    AN = TN.copy()
                    used_chain_deaggregation = True
                else:
                    # MATLAB: for single-server closed product-form models, use egflin (linearizer)
                    # instead of AQL (Schweitzer) - see solver_amva.m lines 50-54
                    total_pop = np.sum(N[np.isfinite(N)])
                    max_servers = 1
                    if mi is not None:
                        finite_servers = mi[np.isfinite(mi)]
                        if len(finite_servers) > 0:
                            max_servers = int(np.max(finite_servers))

                    # MATLAB: if single server model, use egflin; otherwise use lin
                    if total_pop > 2 and np.all(N >= 1) and max_servers == 1:
                        # Switch to egflin method - matches MATLAB's solver_amva.m behavior
                        from ...api.pfqn import pfqn_egflinearizer

                        sched_type = []
                        for q_idx in queue_indices:
                            if self.sched is not None and q_idx in self.sched:
                                sched_type.append(str(self.sched[q_idx]))
                            else:
                                sched_type.append('FCFS')

                        # egflin alpha argument must not be omitted (it distinguishes egflin from lin); see _kb/07-cross-language-parity.md egflin alpha-collapse trap.
                        _N_arr = np.asarray(N, dtype=float).ravel()
                        _alphaM = np.zeros(len(_N_arr))
                        for _r in range(len(_N_arr)):
                            if np.isfinite(_N_arr[_r]):
                                _alphaM[_r] = 0.6 + 1.4 * np.exp(-8 * np.exp(-0.8 * _N_arr[_r]))
                        QN_out, UN_out, WN_out, TN_out, CN_out, XN_out, iter_count = pfqn_egflinearizer(
                            L, N, Z, sched_type, 1e-4, 1000, _alphaM
                        )

                        # XN_out is (1, R) or (R,), flatten to (R,)
                        XN_out = XN_out.flatten()

                        # Compute response metrics
                        RN_out = np.zeros((len(queue_indices), R))
                        AN_out = np.zeros((len(queue_indices), R))
                        for r in range(R):
                            if XN_out[r] > 0:
                                for idx in range(len(queue_indices)):
                                    if L[idx, r] > 0:
                                        RN_out[idx, r] = QN_out[idx, r] / XN_out[r]
                                        AN_out[idx, r] = XN_out[r]

                        for idx, q_idx in enumerate(queue_indices):
                            QN[q_idx, :] = QN_out[idx, :]
                            UN[q_idx, :] = UN_out[idx, :]
                            RN[q_idx, :] = RN_out[idx, :]
                            # T = V .* X, as MATLAB solver_amva.m:291 builds it,
                            # NOT the linearizer's fourth output: that is the
                            # per-reference-visit throughput and MATLAB discards
                            # it (`~,~`) for this reason. V is recovered as
                            # demand * rate = (V*S) * (1/S).
                            for r in range(R):
                                v = (self.demands[q_idx, r] * self.rates[q_idx, r]
                                     if self.rates[q_idx, r] > 0 else 0.0)
                                TN[q_idx, r] = XN_out[r] * v
                            AN[q_idx, :] = TN[q_idx, :]

                        XN = XN_out.flatten()
                    else:
                        # Fall back to AQL for small populations or multiserver
                        result = pfqn_aql(L, N, Z)
                        XN_out, CN_out, QN_out, UN_out, RN_out, TN_out, AN_out = result

                        for idx, q_idx in enumerate(queue_indices):
                            QN[q_idx, :] = QN_out[idx, :]
                            UN[q_idx, :] = UN_out[idx, :]
                            RN[q_idx, :] = RN_out[idx, :]
                            TN[q_idx, :] = TN_out[idx, :]
                            AN[q_idx, :] = AN_out[idx, :]

                        XN = XN_out.flatten()

            elif method == 'bs' and not self._amva_softmin_multiserver(mi):
                line_debug("Using bound method: %s", method, options=self.options)
                # Bard-Schweitzer needs solver tol/iter, warm start, per-station sched explicit, else api defaults (looser tol, all-PS); mirrors solver_amva.m:148.
                from ...lang.base import SchedStrategy as _SchedBase
                _bs_sched = [self.sched[q_idx] if (self.sched is not None and q_idx in self.sched)
                             else _SchedBase.PS for q_idx in queue_indices]
                _bs_tol = getattr(self.options, 'tol', None) or 1e-4
                _bs_imax = getattr(self.options, 'iter_max', None) or 1000
                # pfqn_bs has no server-count argument, so a multiserver model must be
                # transformed before it is handed over; MATLAB solver_amva.m does this
                # once for the whole product-form arm (:111-117), upstream of its bs
                # case (:154), which is why its bs honours nservers and this did not.
                _bs_L, _bs_Z = L, Z
                _bs_max_servers, _bs_rule = self._amva_multiserver_rule(mi)
                if _bs_max_servers > 1 and _bs_rule in ('default', 'seidmann'):
                    _bs_L, _bs_Z = self._amva_seidmann(L, Z, mi)
                XN_out, QN_out, UN_out, RN_out, _bsiter = pfqn_bs(
                    _bs_L, N, _bs_Z, _bs_tol, _bs_imax, None, _bs_sched)
                self._lastiter = _bsiter
                self._lastiterbudget = _bs_imax
                TN_out = np.tile(XN_out, (QN_out.shape[0], 1))
                AN_out = TN_out.copy()

                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = AN_out[idx, :]

                XN = XN_out.flatten()
                if _bs_max_servers > 1 and _bs_rule in ('default', 'seidmann'):
                    self._amva_seidmann_unapply(QN, RN, TN, queue_indices, L, mi, XN)

            elif method == 'aql':
                # Aggregate Queue Length: K+1 population points plus the gamma
                # correction; MATLAB solver_amva.m:160 rejects multiserver here.
                from ...api.pfqn import pfqn_aql
                if self._amva_multiserver_rule(mi)[0] > 1:
                    raise ValueError(
                        "AQL cannot handle multi-server stations. "
                        "Try with the 'default' or 'lin' methods.")
                _aql_tol = getattr(self.options, 'tol', None) or 1e-7
                _aql_imax = getattr(self.options, 'iter_max', None) or 1000
                XN_out, _aqlCN, QN_out, UN_out, RN_out, TN_out, _aqlAN = pfqn_aql(
                    L, N, Z, _aql_tol, _aql_imax)
                self._lastiterbudget = _aql_imax
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = XN_out.flatten()

            elif method == 'qsa':
                # Queue-Shift Approximation: the absolute shift of the
                # aggregate queue length, solved by damped Newton over the
                # quintuple (16); MATLAB solver_amva.m rejects multiserver here
                # exactly as it does for aql.
                from ...api.pfqn import pfqn_qsa
                from ...lang.base import SchedStrategy as _SchedBase
                if self._amva_multiserver_rule(mi)[0] > 1:
                    raise ValueError(
                        "QSA cannot handle multi-server stations. "
                        "Try with the 'default' or 'lin' methods.")
                _qsa_sched = [self.sched[q_idx] if (self.sched is not None and q_idx in self.sched)
                              else _SchedBase.PS for q_idx in queue_indices]
                _qsa_tol = getattr(self.options, 'tol', None) or 1e-10
                _qsa_imax = getattr(self.options, 'iter_max', None) or 100
                QN_out, UN_out, RN_out, _qsaCN, XN_out, _qsaiter = pfqn_qsa(
                    L, N, Z, _qsa_sched, _qsa_tol, _qsa_imax)
                self._lastiter = _qsaiter
                self._lastiterbudget = _qsa_imax
                TN_out = np.tile(np.asarray(XN_out).reshape(1, -1), (QN_out.shape[0], 1))
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = np.asarray(XN_out).flatten()

            elif method == 'tay':
                # Tay's arrival-instant approximation: the arrival-instant queue lengths
                # come from the throughput elasticities, not from a population shift.
                from ...api.pfqn import pfqn_tay
                _tay_tol = getattr(self.options, 'tol', None) or 1e-6
                _tay_imax = getattr(self.options, 'iter_max', None) or 1000
                XN_out, QN_out, UN_out, RN_out, _tayiter, _ = pfqn_tay(
                    L, N, Z, _tay_tol, _tay_imax)
                self._lastiter = _tayiter
                self._lastiterbudget = _tay_imax
                TN_out = np.tile(XN_out, (QN_out.shape[0], 1))
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = XN_out.flatten()

            elif method == 'scat':
                # Neuse-Chandy SCAT: the Linearizer fixed point with a single
                # Delta refresh. Multiserver stations are Seidmann-scaled here
                # exactly as the bs arm does; MATLAB solver_amva.m applies the
                # transform once for the whole product-form arm instead.
                from ...api.pfqn import pfqn_scat
                _sc_L, _sc_Z = L, Z
                _sc_max_servers, _sc_rule = self._amva_multiserver_rule(mi)
                _sc_seidmann = _sc_max_servers > 1 and _sc_rule in ('default', 'seidmann')
                if _sc_seidmann:
                    _sc_L, _sc_Z = self._amva_seidmann(L, Z, mi)
                _sc_tol = getattr(self.options, 'tol', None) or 1e-8
                _sc_imax = getattr(self.options, 'iter_max', None) or 1000
                QN_out, UN_out, WN_out, _scTN, _scCN, XN_out, _sciter = pfqn_scat(
                    _sc_L, N, _sc_Z, None, _sc_tol, _sc_imax)
                self._lastiter = _sciter
                self._lastiterbudget = _sc_imax
                XN = np.asarray(XN_out).flatten()
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = WN_out[idx, :]
                    # T = V .* X, as in the lin family: the fourth output is the
                    # throughput per REFERENCE VISIT, not the station throughput
                    for r in range(R):
                        v = (self.demands[q_idx, r] * self.rates[q_idx, r]
                             if self.rates[q_idx, r] > 0 else 0.0)
                        TN[q_idx, r] = XN[r] * v
                    AN[q_idx, :] = TN[q_idx, :]
                if _sc_seidmann:
                    self._amva_seidmann_unapply(QN, RN, TN, queue_indices, L, mi, XN)

            elif method in ('lcp', 'chow'):
                # Bard LCP and the Chow Second Approximation built on it. Both
                # are Bard-Schweitzer variants in the arrival-instant estimate,
                # so they take the same Seidmann treatment as the bs arm.
                from ...api.pfqn import pfqn_lcp, pfqn_chow
                from ...lang.base import SchedStrategy as _SchedBase
                _cw_L, _cw_Z = L, Z
                _cw_max_servers, _cw_rule = self._amva_multiserver_rule(mi)
                _cw_seidmann = _cw_max_servers > 1 and _cw_rule in ('default', 'seidmann')
                if _cw_seidmann:
                    _cw_L, _cw_Z = self._amva_seidmann(L, Z, mi)
                _cw_sched = [self.sched[q_idx] if (self.sched is not None and q_idx in self.sched)
                             else _SchedBase.PS for q_idx in queue_indices]
                _cw_tol = getattr(self.options, 'tol', None) or 1e-6
                _cw_imax = getattr(self.options, 'iter_max', None) or 1000
                _cw_fn = pfqn_lcp if method == 'lcp' else pfqn_chow
                XN_out, QN_out, UN_out, RN_out, _cwiter = _cw_fn(
                    _cw_L, N, _cw_Z, _cw_tol, _cw_imax, None, _cw_sched)
                self._lastiter = _cwiter
                self._lastiterbudget = _cw_imax
                TN_out = np.tile(XN_out, (QN_out.shape[0], 1))
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = XN_out.flatten()
                if _cw_seidmann:
                    self._amva_seidmann_unapply(QN, RN, TN, queue_indices, L, mi, XN)

            elif method in ('pamb', 'pami', 'pamt'):
                # Hsieh-Lam proportional approximations, noniterative
                from ...api.pfqn import pfqn_pam
                XN_out, QN_out, UN_out, RN_out = pfqn_pam(L, N, Z, method)
                self._lastiter = 1
                TN_out = np.tile(XN_out, (QN_out.shape[0], 1))
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = XN_out.flatten()

            elif method == 'clust':
                # de Souza e Silva-Lavenberg-Muntz clustering approximation
                from ...api.pfqn import pfqn_clust
                _cl_tol = getattr(self.options, 'tol', None) or 1e-6
                _cl_imax = getattr(self.options, 'iter_max', None) or 1000
                XN_out, QN_out, UN_out, RN_out, _cliter = pfqn_clust(
                    L, N, Z, None, None, 'lin', _cl_tol, _cl_imax)
                self._lastiter = _cliter
                self._lastiterbudget = _cl_imax
                TN_out = np.tile(XN_out, (QN_out.shape[0], 1))
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = RN_out[idx, :]
                    TN[q_idx, :] = TN_out[idx, :]
                    AN[q_idx, :] = TN_out[idx, :]
                XN = XN_out.flatten()

            elif method == 'dmlin':
                # de Souza e Silva-Muntz Improved Linearizer: the Linearizer
                # fixed point reached with the Delta-terms pre-aggregated, so
                # the answer matches the lin arm at lower cost.
                from ...api.pfqn import pfqn_dmlin
                _dm_L, _dm_Z = L, Z
                _dm_max_servers, _dm_rule = self._amva_multiserver_rule(mi)
                _dm_seidmann = _dm_max_servers > 1 and _dm_rule in ('default', 'seidmann')
                if _dm_seidmann:
                    _dm_L, _dm_Z = self._amva_seidmann(L, Z, mi)
                _dm_tol = getattr(self.options, 'tol', None) or 1e-8
                _dm_imax = getattr(self.options, 'iter_max', None) or 1000
                QN_out, UN_out, WN_out, _dmTN, _dmCN, XN_out, _dmiter = pfqn_dmlin(
                    _dm_L, N, _dm_Z, None, _dm_tol, _dm_imax)
                self._lastiter = _dmiter
                self._lastiterbudget = _dm_imax
                XN = np.asarray(XN_out).flatten()
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :]
                    UN[q_idx, :] = UN_out[idx, :]
                    RN[q_idx, :] = WN_out[idx, :]
                    # T = V .* X, as in the lin family
                    for r in range(R):
                        v = (self.demands[q_idx, r] * self.rates[q_idx, r]
                             if self.rates[q_idx, r] > 0 else 0.0)
                        TN[q_idx, r] = XN[r] * v
                    AN[q_idx, :] = TN[q_idx, :]
                if _dm_seidmann:
                    self._amva_seidmann_unapply(QN, RN, TN, queue_indices, L, mi, XN)

            elif method == 'sqni':
                # Square-root Non-iterative. pfqn_sqni is a closed form for one
                # queueing station with a delay; with more stations it read only
                # the first demand row and reported those numbers as the answer.
                if self.nstations != 2 or len(queue_indices) != 1:
                    raise ValueError(
                        "SQNI is defined for a single queueing station with a delay. "
                        "Try with the 'default' or 'lin' methods.")
                QN_out, UN_out, XN_out = pfqn_sqni(L, N, Z)
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_out[idx, :] if QN_out.ndim > 1 else QN_out
                    UN[q_idx, :] = UN_out[idx, :] if UN_out.ndim > 1 else UN_out
                XN = XN_out.flatten()

            elif (method in ['lin', 'gflin', 'egflin'] and not _amva_needs_amvald(self._sn)
                    and not self._lin_family_needs_amvald(mi)):
                line_debug("Standard queueing network, routing to mva_analyzer (method=%s)", method, options=self.options)
                # Linearizer family of algorithms
                # Build scheduling strategy list
                sched_type = []
                for q_idx in queue_indices:
                    if self.sched is not None and q_idx in self.sched:
                        sched_type.append(str(self.sched[q_idx]))
                    else:
                        sched_type.append('FCFS')

                # Check if this network has open classes (pure open or mixed)
                has_open_classes = np.any(np.isinf(N))
                is_mixed = has_open_classes and np.any(np.isfinite(N) & (N > 0))

                if has_open_classes:
                    # Get number of servers from model nodes
                    nservers = np.ones(len(queue_indices))
                    model_nodes = self.model.get_nodes() if hasattr(self.model, 'get_nodes') else []
                    for idx, q_idx in enumerate(queue_indices):
                        # Try to get nservers from model nodes first
                        if q_idx < len(model_nodes):
                            node = model_nodes[q_idx]
                            if hasattr(node, 'get_number_of_servers'):
                                ns = node.get_number_of_servers()
                                if ns is not None:
                                    nservers[idx] = ns
                        # Fallback to _sn.nservers if available
                        elif hasattr(self, '_sn') and self._sn is not None:
                            if hasattr(self._sn, 'nservers') and self._sn.nservers is not None:
                                if q_idx < len(self._sn.nservers):
                                    nservers[idx] = self._sn.nservers[q_idx]

                    finite_servers = nservers[np.isfinite(nservers)]
                    max_servers = int(np.max(finite_servers)) if len(finite_servers) > 0 else 1

                    # open/mixed networks need solver_amvald+sn_deaggregate for TN=XN*V; pfqn_linearizermx has no visits so TN=XN; mirrors solver_amva.m:208-210,264.
                    if has_open_classes:
                        # Pure open network OR multiserver mixed network: use solver_amvald
                        # Reference: MATLAB solver_amva.m lines 208-210, 264
                        from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                        from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                        # Get chain-level parameters
                        chain_result = sn_get_demands_chain(self._sn)
                        Lchain = chain_result.Lchain
                        STchain = chain_result.STchain
                        Vchain = chain_result.Vchain
                        alpha = chain_result.alpha
                        Nchain = chain_result.Nchain.flatten()
                        refstatchain = chain_result.refstatchain
                        SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                        # Set up options with 'lin' method
                        amvald_options = AmvaldOptions(
                            method=method,
                            iter_tol=self.options.iter_tol,
                            iter_max=1000,
                            init_sol=getattr(self.options, 'init_sol', None)
                        )
                        if self.options.config and 'multiserver' in self.options.config:
                            amvald_options.config = AmvaldOptions.Config(
                                multiserver=self.options.config['multiserver']
                            )

                        # Call solver_amvald
                        self._apply_interlock(amvald_options)
                        result = solver_amvald(
                            self._sn, Lchain, STchain, Vchain, alpha,
                            Nchain, SCVchain, refstatchain, amvald_options
                        )
                        self._lastiter = getattr(result, 'totiter', None)
                        self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                        self._lastconverged = getattr(result, 'converged', None)

                        # empty Uchain passed when lldscaling/cdscaling are unset so deaggregate computes U with /nservers division; mirrors MATLAB.
                        has_lldscaling = hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None
                        has_cdscaling = hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None
                        has_jdscaling = getattr(self._sn, 'jdscaling', None) is not None
                        Uchain_for_deagg = result.U if (has_lldscaling or has_cdscaling or has_jdscaling) else None
                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, Uchain_for_deagg, result.R, result.T, None, result.X
                        )

                        # Copy disaggregated results
                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True
                    else:
                        # Single-server mixed network: use pfqn_linearizermx directly
                        # Reference: MATLAB solver_amva.m lines 199-201

                        # Compute arrival rates for open classes
                        lambda_arr = np.zeros(R)
                        source_indices = self._get_source_stations()

                        for r in range(R):
                            if np.isinf(N[r]):
                                # Open class - get arrival rate from Source station
                                if hasattr(self, '_sn') and self._sn is not None:
                                    rates = self._sn.rates
                                    if rates is not None:
                                        # Look at Source station indices to get arrival rate
                                        for src_idx in source_indices:
                                            if src_idx < rates.shape[0] and r < rates.shape[1]:
                                                if rates[src_idx, r] > 0 and not np.isinf(rates[src_idx, r]):
                                                    lambda_arr[r] = rates[src_idx, r]
                                                    break

                                # Fallback: try to get arrival rate from model Source node
                                if lambda_arr[r] == 0 and np.any(L[:, r] > 0):
                                    if hasattr(self, 'model') and hasattr(self.model, 'get_classes'):
                                        classes = self.model.get_classes()
                                        if r < len(classes):
                                            job_class = classes[r]
                                            for node in self.model.get_nodes():
                                                if hasattr(node, 'get_arrival') and hasattr(node, '__class__'):
                                                    if 'Source' in node.__class__.__name__:
                                                        arr_dist = node.get_arrival(job_class)
                                                        if arr_dist is not None and hasattr(arr_dist, 'get_rate'):
                                                            lambda_arr[r] = arr_dist.get_rate()
                                                            break

                        QN_out, UN_out, WN_out, TN_out, CN_out, XN_out, iters = pfqn_linearizermx(
                            lambda_arr, L, N, Z, nservers, sched_type, method=method
                        )
                        self._lastiter = iters
                        self._lastiterbudget = None

                        for idx, q_idx in enumerate(queue_indices):
                            QN[q_idx, :] = QN_out[idx, :]
                            UN[q_idx, :] = UN_out[idx, :]
                            RN[q_idx, :] = WN_out[idx, :]  # Response times
                            # T = V .* X (MATLAB solver_amva.m:291), NOT the
                            # linearizer's fourth output: that is the throughput
                            # per REFERENCE VISIT, and MATLAB discards it (`~,~`)
                            # for exactly this reason. V is recovered as
                            # demand * rate = (V*S) * (1/S).
                            _X = np.asarray(XN_out).flatten()
                            for r in range(R):
                                v = (self.demands[q_idx, r] * self.rates[q_idx, r]
                                     if self.rates[q_idx, r] > 0 else 0.0)
                                TN[q_idx, r] = _X[r] * v
                            AN[q_idx, :] = TN[q_idx, :]  # Arrival rate = throughput

                        XN = XN_out.flatten()
                else:
                    # Check for class switching - requires chain-level approach
                    has_class_switching_lin = False
                    if hasattr(self._sn, 'nchains') and self._sn.nchains > 0:
                        chains = self._get_chains()
                        for chain in chains:
                            if len(chain) > 1:
                                has_class_switching_lin = True
                                break

                    if has_class_switching_lin:
                        # For class-switching networks, use chain-level approach
                        # MATLAB uses chain aggregation + solver_amvald + disaggregation
                        from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                        from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                        # Get chain-level parameters
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
                            method=method,
                            iter_tol=self.options.iter_tol,
                            iter_max=1000,
                            init_sol=getattr(self.options, 'init_sol', None)
                        )

                        # Call solver_amvald
                        self._apply_interlock(amvald_options)
                        result = solver_amvald(
                            self._sn, Lchain, STchain, Vchain, alpha,
                            Nchain, SCVchain, refstatchain, amvald_options
                        )
                        self._lastiter = getattr(result, 'totiter', None)
                        self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                        self._lastconverged = getattr(result, 'converged', None)

                        # Disaggregate chain results to class level
                        # MATLAB: pass Uchain only if there's load/class-dependent scaling
                        has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                          np.any(self._sn.lldscaling != 0))
                        has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                          np.any(self._sn.cdscaling != 0))
                        Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                        deagg = sn_deaggregate_chain_results(
                            self._sn, Lchain, None, STchain, Vchain, alpha,
                            None, Uchain_for_deagg, result.R, result.T, None, result.X
                        )

                        # Copy disaggregated results
                        QN = deagg.Q
                        UN = deagg.U
                        RN = deagg.R
                        TN = deagg.T
                        XN = deagg.X.flatten()
                        AN = TN.copy()
                        used_chain_deaggregation = True
                    else:
                        _max_servers, _ms_rule = self._amva_multiserver_rule(mi)
                        # The routing rules that send a multiserver model to solver_amvald were
                        # taken above; the two that remain are served by their own algorithms,
                        # both of which take nservers. The single-server linearizer below does
                        # NOT, so reaching it with m>1 would silently solve the single-server
                        # model. Mirrors MATLAB solver_amva.m:246-252.
                        if _max_servers > 1 and _ms_rule == 'conway':
                            from ...api.pfqn import pfqn_conwayms
                            # returns (Q, U, R, C, X, totiter); the common unpack below
                            # reads (Q, U, W, T, C, X, iter) and recomputes T from X and
                            # the visits, so TN is passed as None deliberately.
                            # Mirrors MATLAB solver_amva.m:249.
                            _Qcw, _Ucw, _Rcw, _Ccw, _Xcw, _itcw = pfqn_conwayms(
                                L, N, Z, np.asarray(mi, dtype=float).ravel(), sched_type,
                                self.options.tol, 1000)
                            result = (_Qcw, _Ucw, _Rcw, None, _Ccw, _Xcw, _itcw)
                        elif _max_servers > 1 and _ms_rule == 'krzesinski':
                            from ...api.pfqn import pfqn_linearizermx
                            # returns (QN, UN, WN, TN, CN, XN, totiter); the common unpack
                            # below reads (Q, U, W, T, C, X, iter) and recomputes T from X
                            # and the visits, so TN is passed as None deliberately
                            _lam = np.zeros(len(np.asarray(N).ravel()))
                            _Qms, _Ums, _Rms, _, _Cms, _Xms, _itms = pfqn_linearizermx(
                                _lam, L, N, Z, np.asarray(mi, dtype=float).ravel(), sched_type,
                                self.options.tol, 1000, 'default')
                            result = (_Qms, _Ums, _Rms, None, _Cms, _Xms, _itms)
                        elif method == 'egflin':
                            N_arr = np.asarray(N, dtype=float).ravel()
                            alphaM = np.zeros(len(N_arr))
                            for r in range(len(N_arr)):
                                if np.isfinite(N_arr[r]):
                                    alphaM[r] = 0.6 + 1.4 * np.exp(-8 * np.exp(-0.8 * N_arr[r]))
                            result = pfqn_egflinearizer(L, N, Z, sched_type,
                                                        self.options.tol, 1000, alphaM)
                        elif method == 'gflin':
                            result = pfqn_gflinearizer(L, N, Z, sched_type,
                                                       self.options.tol, 1000, 2.0)
                        else:  # lin
                            result = pfqn_linearizer(L, N, Z, sched_type,
                                                     self.options.tol, 1000)

                        QN_out, UN_out, WN_out, TN_out, CN_out, XN_out, iters = result
                        self._lastiter = iters
                        self._lastiterbudget = 1000

                        for idx, q_idx in enumerate(queue_indices):
                            QN[q_idx, :] = QN_out[idx, :]
                            UN[q_idx, :] = UN_out[idx, :]
                            RN[q_idx, :] = WN_out[idx, :]  # Response times
                            # T = V .* X (MATLAB solver_amva.m:291), NOT the
                            # linearizer's fourth output: that is the throughput
                            # per REFERENCE VISIT, and MATLAB discards it (`~,~`)
                            # for exactly this reason. V is recovered as
                            # demand * rate = (V*S) * (1/S).
                            _X = np.asarray(XN_out).flatten()
                            for r in range(R):
                                v = (self.demands[q_idx, r] * self.rates[q_idx, r]
                                     if self.rates[q_idx, r] > 0 else 0.0)
                                TN[q_idx, r] = _X[r] * v
                            AN[q_idx, :] = TN[q_idx, :]  # Arrival rate = throughput

                        XN = XN_out.flatten()

            elif method in ['schmidt', 'schmidt-ext', 'ab']:
                # Schmidt/AB/Akyildiz-Bolch stack the delay row on demands ([Z0;L0]); see _kb/07-cross-language-parity.md SchedStrategy per-callee numbering trap.
                from ...api.pfqn.schmidt import SchedStrategy as _SchedSchmidt
                from ...api.pfqn.ab_amva import SchedStrategy as _SchedAb

                _enum = _SchedAb if method == 'ab' else _SchedSchmidt
                sched_q = []
                for q_idx in queue_indices:
                    # `.name`, not str(): sn.sched holds an IntEnum, and since
                    # python 3.11 str() on one of those is the bare NUMBER ('4'),
                    # so 'PS' in str(sched) was false for every station and the
                    # whole family was told FCFS. A PS station with class-dependent
                    # demands then entered pfqn_schmidt_ext's alpha correction,
                    # which is where the class-switching crash came from, and every
                    # PS and LCFS-PR station was solved by the wrong kernel arm.
                    _sched = self.sched[q_idx] if (self.sched is not None and q_idx in self.sched) else None
                    s_str = getattr(_sched, 'name', None) or str(_sched) if _sched is not None else 'FCFS'
                    if 'INF' in s_str:
                        sched_q.append(int(_enum.INF))
                    elif 'PS' in s_str:
                        sched_q.append(int(_enum.PS))
                    else:
                        sched_q.append(int(_enum.FCFS))

                has_delay = bool(np.any(np.asarray(Z, dtype=float) > 0))
                mi_arr = np.asarray(mi, dtype=float).flatten()
                if has_delay:
                    D_full = np.vstack([np.asarray(Z, dtype=float).reshape(1, -1), L])
                    S_full = np.concatenate([np.ones(1), mi_arr])
                    sched_full = np.array([int(_enum.INF)] + sched_q, dtype=int)
                else:
                    D_full = L
                    S_full = mi_arr
                    sched_full = np.array(sched_q, dtype=int)
                S_int = np.where(np.isfinite(S_full), S_full, 1.0).astype(int)
                V_full = np.ones_like(D_full)
                N_int = np.where(np.isfinite(N), N, 0).astype(int)

                # One predicate for the gate and the run, asked about the numbers
                # THIS arm passes: pfqn_schmidt_ext forms its alpha correction from
                # the network with one class-r customer tagged, and an empty class
                # has none to tag.
                from ...api.solvers.mva.handler import mva_supports_schmidt_ext
                _sx_N, _sx_fcfs = self._schmidt_arm_inputs()
                _sx_ok, _sx_reason = mva_supports_schmidt_ext(_sx_N, _sx_fcfs, method)
                if not _sx_ok:
                    raise ValueError(_sx_reason)

                if method == 'ab':
                    ab_res = pfqn_ab_amva(D_full, N_int, V_full, S_int, sched_full)
                    QN_full = np.asarray(ab_res.QN)
                    XN_raw = np.asarray(ab_res.XN)
                    self._lastiter = getattr(ab_res, 'totiter', None)
                elif method == 'schmidt':
                    XN_raw, QN_full, _UN_full, _CN_full = pfqn_schmidt(
                        D_full, N_int, S_int, sched_full, V_full)
                    XN_raw = np.asarray(XN_raw)
                    QN_full = np.asarray(QN_full)
                else:
                    XN_raw, QN_full, _UN_full, _CN_full = pfqn_schmidt_ext(
                        D_full, N_int, S_int, sched_full, V_full)
                    XN_raw = np.asarray(XN_raw)
                    QN_full = np.asarray(QN_full)

                # pfqn_schmidt reports per-station-class throughput; every row is identical for a closed model, so the first row is the class throughput.
                XN_flat = XN_raw.flatten() if XN_raw.ndim == 1 else np.asarray(XN_raw)[0, :].flatten()

                off = 1 if has_delay else 0
                for idx, q_idx in enumerate(queue_indices):
                    QN[q_idx, :] = QN_full[off + idx, :]
                    c_i = mi_arr[idx] if np.isfinite(mi_arr[idx]) and mi_arr[idx] > 0 else 1.0
                    # U = X D / c, which is what MATLAB recomputes for these
                    # methods rather than reading it back from the algorithm.
                    UN[q_idx, :] = XN_flat * L[idx, :] / c_i
                    for r in range(R):
                        if L[idx, r] > 0:
                            TN[q_idx, r] = XN_flat[r]
                            AN[q_idx, r] = XN_flat[r]
                            if XN_flat[r] > 0:
                                RN[q_idx, r] = QN[q_idx, r] / XN_flat[r]

                XN = XN_flat

            elif method in ('qd', 'qdlin', 'qli', 'fli') or method == 'priomva' or (
                    method in ('lin', 'gflin', 'egflin')
                    and (_amva_needs_amvald(self._sn) or self._lin_family_needs_amvald(mi))) or (
                    method == 'bs' and self._amva_softmin_multiserver(mi)):
                # 'priomva' is UNCONDITIONAL here: the preemptive-resume arm lives in
                # solver_amvald's forward step and nowhere else, so a priomva model that
                # fell through this chain reached the exact-MVA `else` below and was
                # answered WITHOUT its priorities -- silently, and with product-form
                # numbers. MATLAB cannot hit that: solver_mva_analyzer sends the whole
                # amva family to solver_amva, whose non-product-form tail goes to
                # solver_amvald. This chain is python's own shape, so the name is listed
                # explicitly.
                # lin family only under load/class dep; no pfqn_qli/pfqn_fli, via solver_amvald; see _kb/07-cross-language-parity.md egflin alpha-collapse trap.
                from ...api.solvers.mva.amvald import solver_amvald, AmvaldOptions
                from ...api.sn import sn_get_demands_chain, sn_deaggregate_chain_results

                chain_result = sn_get_demands_chain(self._sn)
                Lchain = chain_result.Lchain
                STchain = chain_result.STchain
                Vchain = chain_result.Vchain
                alpha = chain_result.alpha
                Nchain = chain_result.Nchain.flatten()
                refstatchain = chain_result.refstatchain
                SCVchain = np.ones((self._sn.nstations, self._sn.nchains))

                amvald_options = AmvaldOptions(
                    method=method,
                    iter_tol=self.options.iter_tol,
                    iter_max=getattr(self.options, 'iter_max', 1000) or 1000,
                    init_sol=getattr(self.options, 'init_sol', None)
                )
                self._apply_interlock(amvald_options)
                _ms = self._amva_multiserver_rule(mi)[1]
                # solver_amvald has no arm for these; MATLAB remaps them at solver_amva.m:397-401
                amvald_options.config.multiserver = 'default' if _ms in (
                    'conway', 'erlang', 'krzesinski') else _ms
                result = solver_amvald(
                    self._sn, Lchain, STchain, Vchain, alpha,
                    Nchain, SCVchain, refstatchain, amvald_options
                )
                self._lastiter = getattr(result, 'totiter', None)
                self._lastiterbudget = min(int(getattr(amvald_options, 'iter_max', 0) or 0), 10000) or None
                self._lastconverged = getattr(result, 'converged', None)

                # MATLAB passes Uchain down only under load/class-dependent scaling
                # and leaves Q to Little's law (solver_amvald.m:239/241).
                has_ld_scaling = (hasattr(self._sn, 'lldscaling') and self._sn.lldscaling is not None and
                                  np.any(self._sn.lldscaling != 0))
                has_cd_scaling = (hasattr(self._sn, 'cdscaling') and self._sn.cdscaling is not None and
                                  np.any(self._sn.cdscaling != 0))
                Uchain_for_deagg = result.U if (has_ld_scaling or has_cd_scaling) else None

                deagg = sn_deaggregate_chain_results(
                    self._sn, Lchain, None, STchain, Vchain, alpha,
                    None, Uchain_for_deagg, result.R, result.T, None, result.X
                )
                QN = deagg.Q
                UN = deagg.U
                RN = deagg.R
                TN = deagg.T
                XN = deagg.X.flatten()
                AN = TN.copy()
                # class-dependent stations report Util=T*S/peak (sn.cdscalingpeak); mirrors the post-pass MATLAB solver_amvald.m applies to the 'amva' arm.
                from ...api.solvers.mva.amvald import solver_amvald_cd_peak_post, solver_amvald_jd_peak_post
                UN = solver_amvald_cd_peak_post(self._sn, UN, TN)
                UN = solver_amvald_jd_peak_post(self._sn, UN, TN)
                used_chain_deaggregation = True


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
            line_debug("Single-class open queueing system (Source-Queue-Sink), routing to qsys_analyzer", options=self.options)
            # Open network - use QNA or queueing system formulas
            from ...api.qsys import (
                qsys_mm1, qsys_mmk, qsys_mg1,
                qsys_gig1_approx_kingman, qsys_gig1_approx_gelenbe,
                qsys_gig1_approx_heyman, qsys_gig1_approx_kimura,
                qsys_gig1_approx_kobayashi, qsys_gig1_approx_klb,
                qsys_gig1_approx_marchal,
                qsys_gig1_approx_allencunneen, qsys_gigk_approx,
                qsys_gig1_ubnd_kingman, qsys_gigk_approx_kingman,
                qsys_gg1,
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
                # Whitt family: the first three answer a station with
                # ABANDONMENT, which no other analytical solver in LINE does.
                'erlanga', 'mgisrgi', 'gigk.diffusion',
                'gigk.whitt', 'qed', 'gig1.extremal',
            }

            if method in qsys_methods and M == 1 and R == 1:
                # single queue/class exact formulas: ca=sqrt(scv(source)), cs=sqrt(scv(queue)), R=qsys_*, Q=X*R, U=lambda/mu/k; mirrors solver_mva_qsys_analyzer.m.
                sn_scv = self._sn.scv if self._sn.scv is not None else np.ones((self.nstations, R))
                ca = np.sqrt(sn_scv[source_indices[0], 0]) if len(source_indices) > 0 and sn_scv[source_indices[0], 0] > 0 else 1.0
                cs = np.sqrt(sn_scv[queue_indices[0], 0]) if sn_scv[queue_indices[0], 0] > 0 else 1.0

                # Method selection with default resolution (matches MATLAB)
                if method == 'default':
                    if ca == 1.0 and cs == 1.0 and k == 1:
                        method = 'mm1'
                    elif ca == 1.0 and cs == 1.0 and k > 1:
                        method = 'mmk'
                    elif ca == 1.0 and k == 1:
                        method = 'mg1'
                    elif cs == 1.0 and k == 1:
                        method = 'gm1'
                    elif k > 1:
                        method = 'gigk'
                    else:
                        method = 'gig1.klb'

                if method == 'exact':
                    if ca == 1.0 and cs == 1.0 and k == 1:
                        method = 'mm1'
                    elif ca == 1.0 and cs == 1.0 and k > 1:
                        method = 'mmk'
                    elif ca == 1.0 and k == 1:
                        method = 'mg1'
                    elif cs == 1.0 and k == 1:
                        method = 'gm1'

                Rscalar = None

                if method == 'mm1':
                    result = qsys_mm1(lambda_r, mu)
                    Rscalar = result['W']

                elif method == 'mmk':
                    result = qsys_mmk(lambda_r, mu, k)
                    Rscalar = result['W']

                elif method in ['mg1', 'mgi1']:
                    Rscalar = qsys_mg1(lambda_r, mu, cs)['W']

                elif method in ['gm1', 'gim1']:
                    # exact GI/M/1 via PH/M/1 sigma-root (only when sn.proc is an exact Markovian rep), else exact LST sigma-root, else two-moment qsys_gg1 fit.
                    from ...constants import ProcessType as _PTq
                    Rscalar = None
                    src_idx = source_indices[0] if len(source_indices) > 0 else None
                    src_is_markovian = (src_idx is not None and self._sn.procid is not None
                                        and _PTq.isMarkovian(self._sn.procid[src_idx, 0]))
                    if src_is_markovian and hasattr(self._sn, 'proc') and self._sn.proc is not None:
                        try:
                            from ...api.qsys import qsys_phm1
                            pie_p, D0p = _ph_from_proc(self._sn, src_idx)
                            if pie_p is not None and D0p is not None:
                                res_ph = qsys_phm1(pie_p, D0p, mu)
                                Rscalar = res_ph['mean_sojourn_time']
                        except Exception:
                            Rscalar = None
                    if Rscalar is None:
                        Rscalar = _gm1_lst_sojourn(self._sn, src_idx, mu)
                    if Rscalar is None:
                        # Fallback: two-moment sigma-root fit of qsys_gg1
                        Rscalar, _ = qsys_gg1(lambda_r, mu, ca ** 2, 1.0)

                elif method in ['gigk']:
                    Rscalar, _ = qsys_gigk_approx(lambda_r, mu, ca, cs, k)

                elif method in ['gigk.kingman_approx']:
                    Rscalar, _ = qsys_gigk_approx_kingman(lambda_r, mu, ca, cs, k)

                elif method == 'gig1.kingman':
                    Rscalar, _ = qsys_gig1_ubnd_kingman(lambda_r, mu, ca, cs)

                elif method == 'gig1.heyman':
                    Rscalar, _ = qsys_gig1_approx_heyman(lambda_r, mu, ca, cs)

                elif method in ['gig1', 'gig1.allen']:
                    Rscalar, _ = qsys_gig1_approx_allencunneen(lambda_r, mu, ca, cs)

                elif method == 'gig1.kobayashi':
                    Rscalar, _ = qsys_gig1_approx_kobayashi(lambda_r, mu, ca, cs)

                elif method == 'gig1.klb':
                    Rscalar, _ = qsys_gig1_approx_klb(lambda_r, mu, ca, cs)

                elif method == 'gig1.marchal':
                    Rscalar, _ = qsys_gig1_approx_marchal(lambda_r, mu, ca, cs)

                elif method == 'gig1.gelenbe':
                    Rscalar, _ = qsys_gig1_approx_gelenbe(lambda_r, mu, ca, cs)

                elif method == 'gig1.kimura':
                    Rscalar, _ = qsys_gig1_approx_kimura(lambda_r, mu, ca, cs)

                # The Whitt family. These return a full measure set rather than a
                # response time, because a station with abandonment or blocking
                # has a CARRIED throughput below its offered rate: Little's law
                # on lambda would silently overstate the queue.
                qsys_full = None
                if method in ('erlanga', 'mgisrgi', 'gigk.diffusion'):
                    qi0 = queue_indices[0]
                    cap = float(self._sn.cap[qi0]) if getattr(self._sn, 'cap', None) is not None else float('inf')
                    room = float('inf') if not np.isfinite(cap) else max(0.0, cap - k)
                    if method == 'gigk.diffusion':
                        from ...api.qsys import qsys_ggnm_diffusion
                        d = qsys_ggnm_diffusion(lambda_r, mu, k, room, ca, cs)
                        carried = d['throughput']
                        qsys_full = {'Q': d['meanNumber'], 'U': d['utilization'],
                                     'T': carried, 'A': lambda_r,
                                     'R': d['meanNumber'] / carried if carried > 0 else 0.0}
                    else:
                        from ...api.sn.patience import sn_patience_handles
                        from ...api.qsys import qsys_erlanga, qsys_mgisrgi_whitt
                        h = sn_patience_handles(self._sn, qi0, 0)
                        if h is None:
                            raise RuntimeError(
                                "method '%s' needs a reneging patience law on the queue" % method)
                        if method == 'erlanga' or h['isExponential']:
                            a = qsys_erlanga(lambda_r, mu, h['rate'], k, room)
                        else:
                            a = qsys_mgisrgi_whitt(lambda_r, mu, k, room, h['hazard'])
                        carried = a['throughput']
                        # R is Little's law on the CARRIED rate, which is what
                        # every other LINE solver reports at a station that
                        # loses work (checked against SolverCTMC on M/M/1/K and
                        # on M/M/k+M). The per-served-job sojourn time is a
                        # different quantity and stays in the API result.
                        qsys_full = {'Q': a['meanNumber'], 'U': a['utilization'],
                                     'T': carried, 'A': lambda_r,
                                     'R': a['meanNumber'] / carried if carried > 0 else 0.0}
                elif method == 'gigk.whitt':
                    from ...api.qsys import qsys_gigk_approx_whitt
                    Rscalar = qsys_gigk_approx_whitt(lambda_r, mu, ca, cs, k)[0]
                elif method == 'qed':
                    from ...api.qsys import qsys_mmk_qed
                    q = qsys_mmk_qed(lambda_r, mu, k)
                    Rscalar = q['meanWait'] + 1.0 / mu
                elif method == 'gig1.extremal':
                    from ...api.qsys import qsys_gig1_bnds_extremal
                    b = qsys_gig1_bnds_extremal(lambda_r, mu, ca, cs)
                    # The upper end, as gig1.kingman already reports a bound.
                    Rscalar = b['upperBound'] + 1.0 / mu

                if qsys_full is not None:
                    qi = queue_indices[0]
                    RN[qi, 0] = qsys_full['R']
                    QN[qi, 0] = qsys_full['Q']
                    UN[qi, 0] = qsys_full['U']
                    TN[qi, 0] = qsys_full['T']
                    AN[qi, 0] = qsys_full['A']
                    XN[0] = qsys_full['T']
                    if len(source_indices) > 0:
                        TN[source_indices[0], 0] = lambda_r
                # Compute Q, U, T, X from R (matches MATLAB pattern)
                elif Rscalar is not None:
                    qi = queue_indices[0]
                    RN[qi, 0] = Rscalar
                    XN[0] = lambda_r
                    UN[qi, 0] = lambda_r / mu / k
                    TN[qi, 0] = lambda_r
                    AN[qi, 0] = lambda_r
                    QN[qi, 0] = XN[0] * RN[qi, 0]
                    if len(source_indices) > 0:
                        TN[source_indices[0], 0] = lambda_r

            else:
                # Default: use M/M/k formulas for each queue (multiserver support)
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
                            # Get number of servers for this queue
                            k = int(mi[idx]) if idx < len(mi) and np.isfinite(mi[idx]) and mi[idx] > 0 else 1
                            # For M/M/k: rho_total = lambda/(k*mu), utilization per server
                            rho_k = lambda_r / (k * mu) if mu > 0 else 0  # Utilization per server

                            if rho_k < 1:
                                # M/M/k queue-length approximation: reduces to rho/(1-rho) at k=1; for k>1 uses rho_k/(1-rho_k)+rho_total with rho_total=k*rho_k=lambda/mu.
                                rho_total = k * rho_k  # Total offered load
                                if k == 1:
                                    QN[q_idx, r] = rho_k / (1 - rho_k)  # M/M/1 queue length
                                else:
                                    # M/M/k approx: Q ~= rho_k/(1-rho_k)*Pk + rho_total (Pk = Erlang-C all-busy prob), simplified to rho_total+rho_k/(1-rho_k) for moderate loads.
                                    from ...api.qsys import qsys_mmk
                                    result = qsys_mmk(lambda_r, mu, k)
                                    QN[q_idx, r] = result.get('L', rho_total / (1 - rho_k))
                                    RN[q_idx, r] = result.get('W', L[idx, r] / (1 - rho_k))

                                UN[q_idx, r] = rho_k  # Per-server utilization
                                if RN[q_idx, r] == 0:
                                    RN[q_idx, r] = L[idx, r] / (1 - rho_k) if rho_k < 1 else float('inf')
                                TN[q_idx, r] = lambda_r
                                AN[q_idx, r] = lambda_r

                            XN[r] = lambda_r

        if cp_chain is not None:
            # The arm above solved the CHAIN network and wrote its queueing rows.
            # Complete the chain-level picture the way MATLAB's product-form
            # branch does -- a delay holds X Z jobs, every station carries
            # X V, and R follows by Little's law -- and hand it to
            # sn_deaggregate_chain_results, which splits each chain back over its
            # classes by the visit-weighted share alpha. Q and U are left to the
            # deaggregation rather than passed in, as the class-switching AMVA
            # path above does, so the two agree on how a chain is split.
            from ...api.sn import sn_deaggregate_chain_results
            Xchain = np.asarray(XN, dtype=float).reshape(1, -1)
            Qchain = np.array(QN, dtype=float)
            Tchain = np.zeros((self.nstations, R))
            Rchain = np.zeros((self.nstations, R))
            for c in range(R):
                for i in range(self.nstations):
                    if i in cp_delays:
                        Qchain[i, c] = Xchain[0, c] * cp_chain.STchain[i, c] * cp_chain.Vchain[i, c]
                    Tchain[i, c] = Xchain[0, c] * cp_chain.Vchain[i, c]
                    if Tchain[i, c] > 0:
                        Rchain[i, c] = Qchain[i, c] / Tchain[i, c]
            deagg = sn_deaggregate_chain_results(
                self._sn, cp_chain.Lchain, None, cp_chain.STchain, cp_chain.Vchain,
                cp_chain.alpha, None, None, Rchain, Tchain, None, Xchain)
            QN = deagg.Q
            UN = deagg.U
            RN = deagg.R
            TN = deagg.T
            XN = deagg.X.flatten()
            AN = TN.copy()
            R = self.nclasses
            used_chain_deaggregation = True

        # Delay-station metrics (QN=X*D, UN=QN, RN=service time) from throughput and demands, skipped for methods that already computed them via handler.
        skip_delay_recompute = (method in ['exact', 'mva']) and self.network_type != 'open'
        # Also skip for amva with class switching since it uses chain-level disaggregation
        # which already computes correct class-level response times
        if method == 'amva' and has_class_switching_early:
            skip_delay_recompute = True
        # Skip delay recompute if chain-level deaggregation was used (TN already correct)
        if used_chain_deaggregation:
            skip_delay_recompute = True
        if not skip_delay_recompute:
            delay_stations = self._get_delay_stations()
            for d_idx in delay_stations:
                for r in range(R):
                    d_demand = self.demands[d_idx, r]  # Individual station's demand = visits * service_time
                    # Get service time (1/rate), not demand
                    d_service_time = 1.0 / self.rates[d_idx, r] if self.rates[d_idx, r] > 0 else 0.0
                    if d_demand > 0 and XN[r] > 0:
                        QN[d_idx, r] = XN[r] * d_demand  # Little's law: jobs = throughput * demand
                        UN[d_idx, r] = QN[d_idx, r]  # For infinite servers, utilization = queue length
                        RN[d_idx, r] = d_service_time  # Response time = service time (NOT demand)
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

        # Source TN defaults to arrival rate only if unfilled (BMAP effective job rate differs from raw event rate, so analyzer value kept); Source AN is 0.
        for src_idx in source_stations:
            for r in range(R):
                if self.rates[src_idx, r] > 0:
                    if TN[src_idx, r] == 0:
                        TN[src_idx, r] = self.rates[src_idx, r]
                    AN[src_idx, r] = 0.0  # Jobs originate at source, not arrive

        # Compute fork-join synchronization delays if model has fork-join
        if self._has_fork_join():
            line_debug("Fork-join post-processing: computing sync delays", options=self.options)
            self._compute_fork_join_sync_delays(QN, UN, RN, TN, AN, XN)

        # Compute and store cache hit/miss probabilities
        self._compute_cache_hit_miss_probs(XN)

        # Filter TN to set 0 where class doesn't visit station (visit ratio = 0)
        # This fixes cases where TN is incorrectly set at stations not visited by a class
        from ...constants import GlobalConstants
        if hasattr(self._sn, 'visits') and self._sn.visits:
            for chain_id, visits in self._sn.visits.items():
                if isinstance(visits, np.ndarray):
                    for ist in range(min(TN.shape[0], self._sn.nstations)):
                        sf_idx = int(self._sn.stationToStateful[ist]) if ist < len(self._sn.stationToStateful) else -1
                        if sf_idx >= 0 and sf_idx < visits.shape[0]:
                            for r in range(min(TN.shape[1], visits.shape[1])):
                                # Check chain membership - only filter if class r is in this chain
                                if hasattr(self._sn, 'chains') and self._sn.chains is not None:
                                    r_chain = get_chain_for_class(self._sn.chains, r)
                                    if r_chain == chain_id:
                                        if visits[sf_idx, r] < GlobalConstants.Zero:
                                            TN[ist, r] = 0.0
                                            AN[ist, r] = 0.0

        # Compute arrival rates from throughputs using routing matrix
        # MATLAB runAnalyzer.m line 341: AN = sn_get_arvr_from_tput(sn, TN, self.getAvgTputHandles())
        from ...api.sn.getters import sn_get_arvr_from_tput
        AN = sn_get_arvr_from_tput(self._sn, TN, TN)  # Pass TN as TH (non-empty array)

        # Compute residence times from response times (WN = RN * V)
        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self._sn, RN, None)

        # Store results
        runtime = time.time() - start_time
        self._result = {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'AN': AN,
            'XN': XN,
            'WN': WN,
            'runtime': runtime,
            'method': method,
            'iter': self._lastiter,
            # None when the handler reports no flag, which is not the same as False:
            # there the count is the signal. Published alongside it so the native
            # result and the delegated one expose the same two fields.
            'converged': self._lastconverged,
        }

        # AMVA convergence: can't tell converged-early vs ran-out-of-iterations; see _kb/06-solver-catalog.md MVA AMVA convergence flag vs iteration count.
        self._warn_if_not_converged(method)

        # Restore Cache nodetype if it was converted to ClassSwitch during solving
        # (MATLAB restores via getStruct() which returns fresh sn; Python needs explicit restore)
        if hasattr(self, '_cache_indices') and self._cache_indices:
            from ...api.sn.network_struct import NodeType
            for ind in self._cache_indices:
                self._sn.nodetype[ind] = NodeType.CACHE

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            from ..base import method_label
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"MVA analysis [method: {method_label(self.options.method, method)}; type: {method_type('MVA', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")

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
            self._ensureAvgResults()
        self._cap_unstable_open_util()

        # Empty result (e.g. SQD method on an unsupported multichain model): empty table.
        _qn = self._result.get('QN')
        if _qn is None or np.asarray(_qn).size == 0 or np.all(np.isnan(np.asarray(_qn))):
            return pd.DataFrame(columns=['Station', 'JobClass', 'QLen', 'Util',
                                         'RespT', 'ResidT', 'ArvR', 'Tput'])

        # Compute residence times from response times using visit ratios
        from ...api.sn.transforms import sn_get_residt_from_respt
        from ...api.sn.network_struct import NodeType
        WN = sn_get_residt_from_respt(
            self._sn, self._result.get('RN_uncapped', self._result['RN']), None)

        # WN already correct for fork branches (visit ratio 1 per branch, since a Fork sends the full rate lambda to EACH branch); no post-processing needed.

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
                        'ResidT': WN[station_idx, r],
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

        # Filter out near-zero rows (MATLAB: any(sum([QN,UN,RN,TN,AN,WN])>0))
        # Use tolerance to avoid keeping rows with only numerical noise (e.g., 8e-18)
        numeric_cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        row_sums = df[numeric_cols].abs().sum(axis=1)
        tokeep = row_sums > 1e-14
        df = df.loc[tokeep].reset_index(drop=True)

        # Wrap in IndexedTable for consistent formatting
        result = IndexedTable(df)

        if not self._table_silent:
            print(result)

        return result

    # Python-style alias

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
        from ...api.solvers.mva.prob_methods import get_prob_aggr

        if self._result is None:
            self._ensureAvgResults()

        # Accept a station node object (like SolverCTMC.getProbAggr): resolve
        # to the 1-based station index expected by get_prob_aggr.
        if not isinstance(ist, (int, np.integer)):
            ist = ist.get_station_index0() + 1

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            p = prob_aggr_via_cpp(self)['probAggr']
            if not (1 <= int(ist) <= len(p)):
                raise ValueError("station index %r is outside 1..%d" % (ist, len(p)))
            pr = float(p[int(ist) - 1])
            return (float(np.log(pr)) if pr > 0.0 else float('-inf'), pr)

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
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            _p = prob_via_jar(self, 'prob-marg', ist=ist, jclass=jobclass, kind='vector', onebased=True, raw_station=True)
            return np.arange(len(_p)), _p
        from ...api.solvers.mva.prob_methods import get_prob_marg

        if self._result is None:
            self._ensureAvgResults()

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
        from ...api.solvers.mva.prob_methods import get_prob_sys_aggr

        if self._result is None:
            self._ensureAvgResults()

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            pr = float(prob_aggr_via_cpp(self)['probSysAggr'])
            return (float(np.log(pr)) if pr > 0.0 else float('-inf'), pr)

        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                # U is what the mixed branch's open-class product form is written
                # in (getProbSysAggr.m reads self.result.Avg.U), so an adapter
                # carrying only Q made every mixed model an AttributeError.
                self.U = result_dict.get('UN')
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
        from ...api.solvers.mva.prob_methods import get_prob_norm_const_aggr

        if self._result is None:
            self._ensureAvgResults()

        class ResultAdapter:
            def __init__(self, result_dict):
                self.Q = result_dict.get('QN')
                self.prob = None

        result_adapter = ResultAdapter(self._result)
        return get_prob_norm_const_aggr(self._get_network_struct(), result_adapter)

    def _get_network_struct(self):
        # the real NetworkStruct is preferred over the minimal adapter fallback: probability accessors need the full struct (state, index maps, phase fields).
        if getattr(self, '_sn', None) is not None:
            return self._sn

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
            self._ensureAvgResults()
        return self._result['QN'].copy()

    def _cap_unstable_open_util(self) -> None:
        """
        Cap the reported utilization of unstable open queueing stations at 1.0.

        A finite-server station serving an open class is unstable when its
        offered load rho = sum_r Tput[i,r] / (nservers[i] * rate[i,r]) >= 1.
        Per LINE convention such a station is fully saturated, so its
        utilization is reported as 1.0 (split across classes in proportion to
        their offered load) and a single instability warning is emitted. The
        offered load is recomputed from throughput and service rate so the cap
        is independent of whatever value the solver algorithm left in UN (some
        algorithms leave 0, others leave the raw rho). Infinite-server (delay)
        and Source stations are excluded: their "utilization" is the mean
        number of busy servers and may legitimately exceed 1. Idempotent.
        """
        if getattr(self, '_unstable_util_capped', False):
            return
        self._unstable_util_capped = True

        res = self._result
        if not res or res.get('UN') is None or res.get('TN') is None or self._sn is None:
            return
        from ...api.sn.transforms import cap_unstable_open_util, saturate_unstable_open_metrics
        UN, any_unstable = cap_unstable_open_util(res['UN'], res['TN'], self._sn)
        if any_unstable:
            res['UN'] = UN
            if res.get('QN') is not None and res.get('RN') is not None:
                # getAvg.m derives ResidT from the response times the ANALYZER
                # produced and only then overwrites them with Inf, so the
                # pre-cap matrix is kept for the ResidT computation. Reading the
                # capped RN instead turns every saturated station's residence
                # time into Inf, where the reference reports RN*V of the raw
                # (possibly negative, hence zeroed) value.
                res['RN_uncapped'] = np.array(res['RN'], dtype=float, copy=True)
                QN, RN, TN, _ = saturate_unstable_open_metrics(
                    res['QN'], res['RN'], res['TN'], self._sn)
                res['QN'] = QN
                res['RN'] = RN
                res['TN'] = TN
            line_warning("SolverMVA", "The model has unstable queues "
                         "(utilization >= 1); station utilization is reported "
                         "capped at 1.0, queue length and response time as Inf.")

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
            self._ensureAvgResults()
        self._cap_unstable_open_util()
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
            self._ensureAvgResults()
        return self._result['RN'].copy()

    def getAvgResidT(self) -> np.ndarray:
        """
        Get average residence times (M x K).

        Residence time is computed from response time using visit ratios:
        WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]

        Returns:
            ResidT: Residence times matrix (M x K)
        """
        if self._result is None:
            self._ensureAvgResults()

        # Compute ResidT using proper visit ratios from network structure
        if self._sn is not None and self._sn.visits:
            return sn_get_residt_from_respt(
                self._sn, self._result.get('RN_uncapped', self._result['RN']), None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            return self._result['RN'].copy()

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
            self._ensureAvgResults()

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
            self._ensureAvgResults()
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
            self._ensureAvgResults()
        return self._result['AN'].copy()

    def getAvgSysRespT(self) -> np.ndarray:
        """
        Get system response times (cycle times).

        Returns:
            C: Cycle times vector (K,)
               C[r] = average time for one cycle for class r
               (in closed networks: time between successive visits to reference station)

        Note:
            For closed networks: uses Little's Law C = N/X
            For open networks: sum of response times across all stations
        """
        if self._result is None:
            self._ensureAvgResults()

        R = self._result['RN']
        X = self._result['XN']
        nclasses = R.shape[1]
        C = np.zeros(nclasses)

        for r in range(nclasses):
            if np.isfinite(self.njobs[r]):
                # Closed class: use Little's Law (matching MATLAB getAvgSys.m line 135)
                if X[r] > 0:
                    C[r] = self.njobs[r] / X[r]
                else:
                    C[r] = np.inf
            else:
                # Open class: sum of response times across all stations
                C[r] = np.sum(R[:, r])

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
            self._ensureAvgResults()
        return self._result['XN'].copy()

    # ============================================================================
    # Unified Metrics and Chain/Node/System Methods
    # ============================================================================


    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self, '_sn') and self._sn is not None and hasattr(self._sn, 'chains') and self._sn.chains is not None:
            chains_arr = np.asarray(self._sn.chains)
            if chains_arr.size == 0:
                return [[k] for k in range(self.nclasses)]

            nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else 1

            # Check if chains is 1D (class->chain mapping) or 2D (chain,class membership)
            if chains_arr.ndim == 1:
                # 1D format: chains[k] = c means class k belongs to chain c
                chains = [[] for _ in range(nchains)]
                for k in range(self.nclasses):
                    if k < len(chains_arr):
                        c = int(chains_arr[k])
                        if 0 <= c < nchains:
                            chains[c].append(k)
                return chains if any(chains) else [[k for k in range(self.nclasses)]]
            else:
                # 2D format: chains[c, k] > 0 means class k is in chain c
                chains = []
                for c in range(nchains):
                    chain_classes = []
                    for k in range(self.nclasses):
                        if c < chains_arr.shape[0] and k < chains_arr.shape[1] and chains_arr[c, k] > 0:
                            chain_classes.append(k)
                    chains.append(chain_classes)
                return chains if chains else [[k for k in range(self.nclasses)]]
        else:
            return [[k] for k in range(self.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()

        Q = self._result['QN']
        chains = self._get_chains()
        nchains = len(chains)

        QN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                QN_chain[:, c] = np.sum(Q[:, chain_classes], axis=1)

        # Clean up tiny values (numerical noise) to exactly 0
        QN_chain[np.abs(QN_chain) < 1e-10] = 0.0
        return QN_chain

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain.

        Cleans up tiny numerical values (< 1e-10) to exactly 0.
        """
        if self._result is None:
            self._ensureAvgResults()
        self._cap_unstable_open_util()

        U = self._result['UN']
        chains = self._get_chains()
        nchains = len(chains)

        UN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                UN_chain[:, c] = np.sum(U[:, chain_classes], axis=1)

        # Clean up tiny values (numerical noise) to exactly 0
        UN_chain[np.abs(UN_chain) < 1e-10] = 0.0
        return UN_chain

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain.

        Uses alpha-weighted sum matching MATLAB: RN(:,c) = sum(RNclass(:,inchain).*alpha(:,inchain),2)
        """
        if self._result is None:
            self._ensureAvgResults()

        R = self._result['RN']
        chains = self._get_chains()
        nchains = len(chains)

        RN_chain = np.zeros((self.nstations, nchains))

        # Get alpha weights from sn_get_demands_chain
        if hasattr(self, '_sn') and self._sn is not None:
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
        """Get average residence times aggregated by chain.

        Residence time accounts for visit ratios, computed as:
        WN(i,c) = sum(WNclass(i, inchain))
        where WNclass = sn_get_residt_from_respt converts response times to residence times.
        """
        if self._result is None:
            self._ensureAvgResults()

        R = self._result['RN']  # Per-class response times
        chains = self._get_chains()
        nchains = len(chains)

        WN_chain = np.zeros((self.nstations, nchains))

        # Try to compute proper residence times using visit ratios
        if hasattr(self, '_sn') and self._sn is not None:
            from ...api.sn.transforms import sn_get_residt_from_respt
            try:
                # Compute per-class residence times from response times
                WN = sn_get_residt_from_respt(self._sn, R, None)

                # Aggregate by chain (sum per-class residence times)
                for c, chain_classes in enumerate(chains):
                    if chain_classes:
                        WN_chain[:, c] = np.sum(WN[:, chain_classes], axis=1)
                # Clean up tiny values (numerical noise) to exactly 0
                WN_chain[np.abs(WN_chain) < 1e-10] = 0.0
                return WN_chain
            except Exception:
                pass

        # Fallback: return response times if residence time computation fails
        # (getAvgRespTChain already has cleanup)
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain.

        Sums per-station throughputs for all classes in chain:
        TN(:,c) = sum(TNclass(:, inchain), 2)
        """
        if self._result is None:
            self._ensureAvgResults()

        # Use per-station throughputs TN, not system throughput XN
        TN = self._result['TN']
        chains = self._get_chains()
        nchains = len(chains)

        TN_chain = np.zeros((self.nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                # Sum per-station throughputs for classes in this chain
                TN_chain[:, c] = np.sum(TN[:, chain_classes], axis=1)

        # Clean up tiny values (numerical noise) to exactly 0
        TN_chain[np.abs(TN_chain) < 1e-10] = 0.0
        return TN_chain

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain.

        For most stations, arrival rate equals throughput at steady state.
        For Source nodes, arrival rate is 0 (jobs don't arrive TO a source,
        they depart FROM it).
        """
        AN = self.getAvgTputChain()

        # Set arrival rate to 0 for Source nodes
        if hasattr(self, 'station_types') and self.station_types:
            for i, node_type in enumerate(self.station_types):
                if node_type is not None and node_type == NodeType.SOURCE:
                    AN[i, :] = 0.0

        return AN

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

        # Get station names (use actual names if available)
        station_names = getattr(self, 'station_names', None)
        if station_names is None or len(station_names) != nstations:
            station_names = [f'Station{i}' for i in range(nstations)]

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': station_names[i],
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

        TN = self._result['TN']
        AN = self._result['AN']
        XN = self._result['XN']

        sn = self._sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses

        # Compute actual hit/miss probabilities for Cache nodes
        # Skip if already computed by _run_cache_analysis (which handles replacement strategy properly)
        if sn.nodeparam is not None:
            for ind in range(I):
                if sn.nodetype is not None and ind < len(sn.nodetype):
                    if sn.nodetype[ind] == NodeType.CACHE and ind in sn.nodeparam:
                        cache_param = sn.nodeparam[ind]

                        # Check if hit/miss probs already computed by _run_cache_analysis
                        existing_hit = getattr(cache_param, 'actualhitprob', None)
                        if existing_hit is not None and np.any(existing_hit > 0):
                            # Already computed by _run_cache_analysis, skip
                            continue

                        nitems = getattr(cache_param, 'nitems', 0)
                        cap = getattr(cache_param, 'cap', 0)

                        # Get nitems and cap from model if not in param
                        if nitems == 0 or cap == 0:
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                                for node in self.model._nodes:
                                    if hasattr(node, '_num_items') and hasattr(node, '_item_level_cap'):
                                        nitems = node._num_items if node._num_items else 0
                                        item_cap = node._item_level_cap
                                        if item_cap is not None and (isinstance(item_cap, (list, np.ndarray)) and len(item_cap) > 0):
                                            cap = item_cap[0]
                                        elif item_cap is not None:
                                            cap = item_cap
                                        else:
                                            cap = 0
                                        break

                        if nitems > 0 and cap > 0:
                            # Get cache node for proper cache analysis
                            cache_node = None
                            if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                                cache_node = self.model._nodes[ind]

                            # Default: uniform access formula
                            hit_prob = min(cap / nitems, 1.0)
                            miss_prob = 1.0 - hit_prob

                            # Use proper cache analysis if gamma matrix available
                            if cache_node is not None and hasattr(cache_node, 'get_gamma_matrix'):
                                try:
                                    from ...api.cache import cache_xi_fp
                                    gamma = cache_node.get_gamma_matrix(R)
                                    m_levels = cache_node._item_level_cap if hasattr(cache_node, '_item_level_cap') else np.array([cap])

                                    # Use FPI method which works for large caches
                                    xi, pi0, pij, it = cache_xi_fp(gamma, m_levels)

                                    access_probs = np.sum(gamma, axis=1)
                                    access_probs = access_probs / np.sum(access_probs) if np.sum(access_probs) > 0 else access_probs
                                    hit_prob = np.sum(access_probs * (1 - pi0))
                                    miss_prob = 1.0 - hit_prob
                                except Exception:
                                    pass  # Fall back to uniform

                            hitclass = getattr(cache_param, 'hitclass', np.array([]))
                            nclasses = len(hitclass) if hasattr(hitclass, '__len__') else R
                            cache_param.actualhitprob = np.zeros(nclasses)
                            cache_param.actualmissprob = np.zeros(nclasses)

                            for k in range(nclasses):
                                h = hitclass[k] if k < len(hitclass) else -1
                                missclass = getattr(cache_param, 'missclass', np.array([]))
                                m = missclass[k] if k < len(missclass) else -1
                                if h >= 0 and m >= 0:
                                    cache_param.actualhitprob[k] = hit_prob
                                    cache_param.actualmissprob[k] = miss_prob

                            # Set result on Cache node in model
                            if cache_node is not None:
                                if hasattr(cache_node, 'set_result_hit_prob'):
                                    cache_node.set_result_hit_prob(cache_param.actualhitprob)
                                if hasattr(cache_node, 'set_result_miss_prob'):
                                    cache_node.set_result_miss_prob(cache_param.actualmissprob)

        # Create TH (throughput handle) - indicates which station-classes have valid throughput
        # TH > 0 means the station-class has a valid throughput value
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

        # Copy station metrics to station nodes
        QN = self._result['QN']
        UN = self._result['UN']
        RN = self._result['RN']

        # Compute residence times from response times using visit ratios, off the
        # pre-saturation matrix as getAvg.m does (see _cap_unstable_open_util)
        from ...api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(
            sn, self._result.get('RN_uncapped', RN), None)

        for ist in range(M):
            ind = sn.stationToNode[ist]
            if ind >= 0 and ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]

        # Fix arrival rates for ClassSwitch and Sink nodes for cache hit/miss classes
        # (matches MATLAB getAvgNode.m lines 54-76)
        from ...api.sn.network_struct import NodeType
        for cacheInd in range(I):
            if sn.nodetype is not None and cacheInd < len(sn.nodetype) and sn.nodetype[cacheInd] == NodeType.CACHE:
                if sn.nodeparam is not None and cacheInd in sn.nodeparam:
                    cache_param = sn.nodeparam[cacheInd]
                    hitclass = np.atleast_1d(getattr(cache_param, 'hitclass', np.array([]))).flatten().astype(int)
                    missclass = np.atleast_1d(getattr(cache_param, 'missclass', np.array([]))).flatten().astype(int)
                    for ind in range(I):
                        if sn.nodetype[ind] == NodeType.CLASSSWITCH or sn.nodetype[ind] == NodeType.SINK:
                            for classIdx in range(R):
                                if np.any(hitclass[hitclass >= 0] == classIdx):
                                    ANn[ind, classIdx] = TNn[cacheInd, classIdx]
                                if np.any(missclass[missclass >= 0] == classIdx):
                                    ANn[ind, classIdx] = TNn[cacheInd, classIdx]

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

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics."""
        R = self.getAvgSysRespT()
        T = self.getAvgSysTput()
        return R, T

    getAvgSysTable = NetworkSolver.getAvgSysTable  # chain-level shared layout

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
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import cdf_respt_via_jar
            return cdf_respt_via_jar(self)
        if self._result is None:
            self._ensureAvgResults()

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
        method: str = 'default',
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
            (PercRT, PercTable) where PercRT is a list of dicts with percentile
            data for each (station, class), each holding 'station' (station
            index), 'class' (job class), 'percentiles' (input percentile values)
            and 'values' (percentile response times), and PercTable is a pandas
            DataFrame with columns Station, Class, P10, P25, P50, P75, P90, P95,
            P99, ...

        Algorithm:
            For an exponential CDF with rate lambda = 1/E[R], the percentile p
            is t_p = -ln(1-p) * E[R] with p in [0,1].

        Notes:
            Percentiles outside (0,100) are clipped, stations with zero response
            time give empty lists, and results match the exponential percentile
            formula.

        Example::

            >>> perc_list, perc_table = solver.getPerctRespT([90, 95, 99])
            >>> print(perc_table)
            >>> # Extract 90th percentile response time
            >>> p90_values = perc_list[0]['values']
        """
        self._last_perct_method = str(method or '').lower()
        if method is not None and method.lower() == 'forktail':
            # Fork-join request tail latency; mirrors the MATLAB entry point
            # @NetworkSolver/getPerctRespT.m with method='forktail'
            from ...api.fjnative import forktail_percentiles
            if percentiles is None:
                percentiles = [10, 25, 50, 75, 90, 95, 99]
            return forktail_percentiles(self, percentiles, jobclass)

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
            - AMVA variants: 'bs', 'sqni', 'tay', 'lin', 'gflin', 'egflin', 'schmidt', 'schmidt-ext', 'ab'
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
        methods.extend(['default', 'mva', 'exact', 'sum', 'esum'])

        # MVAC (exact mean value analysis by chain, pfqn_mvac): closed
        # single-server product-form networks only; rejects open/mixed and
        # multiserver at solve time. The gate is FULLY CLOSED, as
        # SolverMVA.m:80's ~any(isinf(njobs)) is -- not merely 'has a closed
        # class'. Testing network_type != 'open' also admitted a MIXED model,
        # which solver_mvac then refused again, so the list named a method the
        # model could not run.
        if self.network_type == 'closed':
            methods.append('mvac')

        # SJN (shortest-job-next, pfqn_mvasjn / pfqn_amvasjn): the conditional
        # waiting time equation is a population recursion, so the family runs on
        # a CLOSED model with an SJF station and nowhere else -- the dispatcher
        # rejects an open one by name. Advertised only there, for the reason
        # 'sqni' is gated: a name on this list is a name a caller is invited to
        # ask for, and the JAR gates `checkDeclaredMethod` on exactly this list.
        if self.network_type == 'closed' and self._has_sjn_station():
            methods.extend(['sjn.mva', 'sjn.amva'])

        # AMVA and variants - available for closed/mixed networks
        methods.extend([
            'amva',
            'bs', 'amva.bs',
            'aql', 'amva.aql',
            'qsa', 'amva.qsa',
            'sqni',
            'tay', 'amva.tay',
            'scat', 'amva.scat',
            'lcp', 'amva.lcp',
            'chow', 'amva.chow',
            'pamb', 'amva.pamb',
            'pami', 'amva.pami',
            'pamt', 'amva.pamt',
            'clust', 'amva.clust',
            'dmlin', 'amva.dmlin',
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

        # SQD (Smith Queue Decomposition) is only valid for closed single-chain
        # Blocking-After-Service networks; solver_sqd returns EMPTY results on
        # anything else, so listing it unconditionally named a method that
        # cannot run. Mirrors SolverMVA.m and the C++ runner.
        from ...api.solvers.mva.analyzers import _is_bas_model
        if _is_bas_model(self._sn):
            methods.append('sqd')

        # SQNI (pfqn_sqni) is a closed form for one queueing station with a
        # delay; listing it elsewhere named a method that cannot run.
        _, _sqni_queues = self._get_queueing_demands()
        if self.nstations != 2 or len(_sqni_queues) != 1:
            methods.remove('sqni')

        # AQL (pfqn_aql), QSA (pfqn_qsa) and Tay (pfqn_tay) reject multiserver
        # stations at solve time, so they are only advertised for single-server
        # models.
        _ns = np.asarray(self.nservers, dtype=float).ravel()
        if np.any(_ns[np.isfinite(_ns)] > 1):
            for _m in ('aql', 'amva.aql', 'qsa', 'amva.qsa', 'tay', 'amva.tay'):
                if _m in methods:
                    methods.remove(_m)

        # QNA and RQNA for open networks (RQNA: robust queueing network
        # analyzer, indices of dispersion, for non-renewal MAP/MMPP arrivals)
        if self.network_type == 'open':
            methods.append('qna')
            methods.append('rqna')
            methods.append('rqt')

        # Marie withheld for open models and for class-dependent routing: the
        # aggregation-decomposition is exact only when every class traverses the
        # network alike. _run_marie has always dispatched it; not listing it
        # hid a working method. Mirrors SolverMVA.m, SolverMVA.java and the C++
        # runner.
        from ...api.solvers.mva.analyzers import _has_classdep_routing
        if self.network_type != 'open' and not _has_classdep_routing(self._sn):
            methods.extend(['marie', 'amva.marie'])

        # amva.mapqn: the horizontal-cut MVA for one exponential delay and one
        # FCFS MAP queue; offered only on that shape, which mva_mapqn_reason
        # judges for the list, the report and the run alike.
        from ...api.solvers.mva.mapqn import mva_mapqn_reason
        if not mva_mapqn_reason(self._sn):
            methods.append('amva.mapqn')

        # priomva: preemptive-resume priority arm (Chandy-Lakshmi [ChaL83]),
        # offered only when a station actually uses FCFSPRPRIO. Mirrors
        # SolverMVA.m; the arm itself lives in solver_amvald's forward step.
        if self._has_prs_prio_station():
            methods.extend(['priomva', 'amva.priomva'])

        # bound methods (aba/bjb/gb/sb/pb/mwba/...) moved to SolverBA and are rejected here; use SolverBA.listValidMethods for the bound catalogue.

        # Queueing system formulas for 2-station open networks
        if self.network_type == 'open' and self.nstations == 2 and self.nclasses == 1:
            methods.extend([
                'mm1', 'mmk', 'mg1', 'mgi1', 'gm1', 'gig1', 'gim1',
                'gig1.kingman', 'gigk', 'gigk.kingman_approx',
                'gig1.gelenbe', 'gig1.heyman', 'gig1.kimura',
                'gig1.allen', 'gig1.kobayashi', 'gig1.klb', 'gig1.marchal',
                # Whitt family. The two abandonment methods are listed only
                # when the station actually reneges: they have nothing to say
                # about a queue nobody leaves, and listing them there would
                # name a method that cannot run.
                'gigk.whitt', 'qed', 'gig1.extremal', 'gigk.diffusion',
            ])
            if self._resolve_abandonment_method() is not None:
                methods.extend(['erlanga', 'mgisrgi'])

        return methods

    def resolveMethod(self, options):
        """Feature-driven resolution of method='default': a bursty single-class
        open network has a non-renewal (MAP/MMPP) arrival process, so the
        default dispatch selects RQNA. The gate then admits the MAP family only
        on this RQNA path. Mirrors the MATLAB/JAR SolverMVA.resolveMethod and the
        analyzer's default->RQNA dispatch below."""
        method = getattr(options, 'method', 'default')
        if method == 'default' and self._sn is not None:
            try:
                from ...api.sn import sn_has_bursty_arrival
                if (self._sn.nclasses == 1
                        and np.all(np.isinf(self._sn.njobs))
                        and sn_has_bursty_arrival(self._sn)):
                    return 'rqna'
            except Exception:
                pass
            # A single-class open station with reneging: resolve to the
            # abandonment method BEFORE the feature gate runs, which is what
            # lets the gate see a method whose featset admits Reneging. The
            # condition is the exact shape solver_mva_qsys_analyzer handles, so
            # any other reneging model still falls through and is refused.
            resolved = self._resolve_abandonment_method()
            if resolved is not None:
                return resolved
        return method

    def _resolve_abandonment_method(self):
        """'erlanga' or 'mgisrgi' when the model is a single-class open
        Source-Queue-Sink whose queue reneges, else None.

        The shape test is the one solver_mva_qsys_analyzer serves -- two
        stations, one class, all open -- and not merely "some station reneges":
        resolving on a wider set would name a method the analyzer's qsys branch
        never reaches, and the model would be answered by the generic MVA path
        under that method's name.
        """
        sn = self._sn
        if sn is None or getattr(sn, 'nclasses', 0) != 1:
            return None
        if getattr(sn, 'nstations', 0) != 2:
            return None
        njobs = getattr(sn, 'njobs', None)
        if njobs is None or not np.all(np.isinf(np.asarray(njobs, dtype=float))):
            return None
        cls = getattr(sn, 'impatienceClass', None)
        if cls is None:
            return None
        from ...lang.base import ImpatienceType
        cls = np.asarray(cls)
        if cls.ndim != 2 or cls.shape[1] < 1:
            return None
        rows = [i for i in range(cls.shape[0])
                if int(cls[i, 0]) == int(ImpatienceType.RENEGING)]
        if len(rows) != 1:
            return None
        from ...api.sn.patience import sn_patience_handles
        h = sn_patience_handles(sn, rows[0], 0)
        if h is None:
            return None
        return 'erlanga' if h['isExponential'] else 'mgisrgi'

    resolve_method = resolveMethod

    def getMethodFeatureSet(self, method):
        """Per-method feature deltas applied to the base MVA envelope. QNA is a
        two-moment open-network method, so it drops closed-class support. The
        queueing-system and bounds methods are already structurally restricted
        by listValidMethods and inherit the base envelope. RQNA adds the
        non-renewal MAP/MMPP family (open only); mirrors the MATLAB/JAR
        SolverMVA.getMethodFeatureSet."""
        from ...api.solvers.mva.handler import (
            mva_base_method, mva_is_closed_population_method,
            MVA_CLOSED_POPULATION_METHODS, MVA_NON_BCMP_SCHED_FEATURES)
        feats = set(SolverMVA.getFeatureSet())
        method = mva_base_method(method)
        if mva_is_closed_population_method(method):
            # The closed-population AMVA family estimates the arrival-instant
            # queue length as a function of the population vector N and is
            # handed (L, N, Z) alone, so an open chain gives it nothing to recur
            # on: the analyzer has no arm for any of these outside its closed
            # product-form branch, and falling through returned the qd-family
            # answer, or a table of zeros, under their name. The remaining
            # precondition of that branch (product form) has no registry name
            # and is applied by supportsClosedPopulation instead.
            feats.discard('OpenClass')
        # The load-dependent analyzer serves a load-, class- or joint-dependent
        # model through 'exact'/'mva' (load dependence only, it has no class- or
        # joint-dependent recursion) and through the default/amva/qd/lin/qdlin
        # arms, and refuses every other name by name. The queueing-system closed
        # forms are intercepted upstream of that analyzer and keep the base
        # envelope.
        if method in set(MVA_CLOSED_POPULATION_METHODS) | {
                'sum', 'esum', 'mvac', 'qli', 'fli', 'gflin', 'egflin',
                'qna', 'rqna', 'rqt'}:
            feats -= {'LoadDependence', 'ClassDependence', 'JointDependence'}
        elif method in ('mva', 'exact'):
            feats -= {'ClassDependence', 'JointDependence'}
        if method not in ('default', 'exact'):
            # An order-independent or pass-and-swap station is served by the
            # exact OI analyzer alone, which the dispatcher reaches only under
            # 'default' or 'exact'; every other name is refused there by name,
            # so it must not be advertised for such a model.
            feats -= {'SchedStrategy_OI', 'SchedStrategy_PAS'}
        if method in ('sum', 'esum'):
            # The summation method passes each station to sum_closed /
            # sum_closing as an INF, PS, LCFS-PR, FCFS or SIRO centre and
            # refuses every other discipline by name.
            feats -= set(MVA_NON_BCMP_SCHED_FEATURES)
        elif method == 'mvac':
            # pfqn_mvac recurs on the closed chains over single-server
            # fixed-rate (SSFR) and infinite-server centres; the handler refuses
            # every other discipline by name.
            feats.discard('OpenClass')
            feats -= set(MVA_NON_BCMP_SCHED_FEATURES)
        if method == 'qna':
            # round-robin dispatching enters as a deterministic traffic split
            # (npfqn_traffic_split_rr), which the exact-MVA paths have no
            # counterpart for
            feats.add('RoutingStrategy_RROBIN')
            feats.discard('ClosedClass')
            feats.discard('SelfLoopingClass')
            # solver_qna's station loop has an arm for INF, PS and FCFS and none
            # for anything else, so on a SIRO, LCFS-PR, HOL or priority station
            # it left that row of Q, U, R and T at zero and reported the table
            # as a solution.
            feats -= {'SchedStrategy_SIRO', 'SchedStrategy_LCFSPR'}
            feats -= set(MVA_NON_BCMP_SCHED_FEATURES)
        elif method == 'rqna':
            feats.update({'MAP', 'MMPP2', 'MMAP', 'RAP'})
            feats.discard('ClosedClass')
            feats.discard('SelfLoopingClass')
        elif method == 'mapqn':
            # the horizontal-cut MVA consumes a MAP service natively (a closed
            # delay + FCFS queue model, see mva_mapqn_reason); declaring MAP
            # here is what keeps needsMapEnv from routing the model through
            # its random-environment image
            feats.update({'MAP', 'MMPP2'})
            feats -= {'OpenClass', 'Source', 'Sink', 'Fork', 'Forker', 'Join', 'Joiner', 'JoinPartial',
                      'ClassSwitch', 'StatelessClassSwitcher', 'Cache', 'CacheClassSwitcher', 'CacheRetrieval',
                      'LoadDependence', 'ClassDependence', 'JointDependence',
                      'SchedStrategy_PS', 'SchedStrategy_SIRO', 'SchedStrategy_LCFSPR',
                      'SchedStrategy_SRPT', 'SchedStrategy_PSJF', 'SchedStrategy_FB', 'SchedStrategy_LRPT',
                      'SchedStrategy_SETF', 'SchedStrategy_OI', 'SchedStrategy_PAS'}
            feats -= set(MVA_NON_BCMP_SCHED_FEATURES)
        elif method == 'rqt':
            # robust queueing theory: single-class open networks, the primitives
            # entering the uncertainty sets are two moments
            feats.discard('ClosedClass')
            feats.discard('SelfLoopingClass')
        if method in ('rqna', 'rqt'):
            # A Join is a synchronisation node, not a queue: it carries no
            # service process, so the index-of-dispersion curve these two read
            # off every station does not exist for it, and neither analyzer has
            # a synchronisation term to put in its place. QNA keeps Fork/Join --
            # its station loop has an explicit Join arm.
            feats -= {'Fork', 'Forker', 'Join', 'Joiner', 'JoinPartial'}
        elif method in ('erlanga', 'mgisrgi'):
            # The ONLY MVA methods that accept abandonment. Reneging is added
            # here rather than to the base envelope on purpose: the base set
            # governs every method, and a multi-station reneging model must go
            # on being refused rather than silently solved without abandonment.
            feats.add('Reneging')
            feats.discard('ClosedClass')
            feats.discard('SelfLoopingClass')
        # MULTISERVER (registry name since 2026-09-05). The single-server
        # recursions: AQL, QSA and Tay (mva_supports_closed_population), MVAC's
        # SSFR chain recursion (mva_supports_mvac), RQNA's GI/G/1 workload
        # (mva_supports_single_class_open), Kant's SJN recursion and the
        # single-server closed forms of the queueing-system analyzer, every
        # M/G/1, G/M/1 and G/G/1 name. Each predicate stays, wording the refusal
        # for the run; the delta is what makes it nameable. RQT, QNA, M/M/k,
        # G/G/k and the rest of the envelope carry a server count.
        if (method in ('aql', 'qsa', 'tay', 'mvac', 'rqna', 'sjn.mva', 'sjn.amva',
                       'mm1', 'mg1', 'mgi1', 'gm1', 'gim1')
                or method.startswith('gig1')):
            feats.discard('MultiServer')
        # FINITECAPACITY (registry name since 2026-09-05) is NOT in the base
        # envelope: the product-form recursions solve a buffer away, which is
        # what supportsFiniteCapacity refuses. The names that honour one are
        # granted it here, and that structural predicate keeps the shape half of
        # each rule. 'default' and 'sqd' reach solver_sqd, the one
        # Blocking-After-Service arm. The single-station M/M/1/K with tail drop
        # is judged on the MODEL because no name can carry it, 'exact' excepted
        # since the closed form is exact at scv=1 only.
        if method in ('default', 'sqd'):
            feats.add('FiniteCapacity')
        elif method != 'exact':
            model = getattr(self, 'model', None)
            if model is not None and hasattr(model, 'getStruct'):
                from ...api.sn import sn_is_mm1k_loss
                if sn_is_mm1k_loss(model.getStruct()):
                    feats.add('FiniteCapacity')
        return feats

    get_method_feature_set = getMethodFeatureSet

    def supportsModelMethod(self, method):
        """Finite station/class capacity has no registry feature name, so the
        coarse per-method feature gate cannot see it. Apply the structural
        capacity check on top of it, otherwise MVA silently returns the
        unconstrained product-form answer for models built with setCapacity /
        a finite classCap (BUG-39). Mirrors MATLAB SolverMVA.supportsModelMethod."""
        from ...api.solvers.mva.handler import (
            mva_base_method, mva_supports_closed_population,
            mva_supports_single_class_open, mva_supports_mvac,
            mva_supports_schmidt_ext)
        ok, reason = super(SolverMVA, self).supportsModelMethod(method)
        model = getattr(self, 'model', None)
        if ok and model is not None and hasattr(model, 'getStruct'):
            ok, reason = SolverMVA.supportsFiniteCapacity(model)
        if ok and model is not None and hasattr(model, 'getStruct'):
            ok, reason = SolverMVA.supportsExactness(model, method)
        if ok and model is not None and hasattr(model, 'getStruct'):
            # Product form, a class count and a server count have no registry
            # feature name, so these three rules cannot live in
            # getMethodFeatureSet. Each is the SAME predicate the analyzer
            # raises on, so a row the report offers is a row that runs.
            sn = model.getStruct()
            ok, reason = mva_supports_closed_population(sn, method)
            if ok:
                ok, reason = mva_supports_single_class_open(sn, method)
            if ok:
                ok, reason = mva_supports_mvac(sn, method)
            if ok:
                from ...api.solvers.mva.mapqn import mva_supports_mapqn
                ok, reason = mva_supports_mapqn(sn, method)
            if ok and mva_base_method(method) == 'schmidt-ext':
                # Built only for the one method that reads them, so no other
                # gate query pays for the demand matrix.
                _N, _fcfs = self._schmidt_arm_inputs()
                ok, reason = mva_supports_schmidt_ext(_N, _fcfs, method)
        return ok, reason

    supports_model_method = supportsModelMethod

    @staticmethod
    def supportsExactness(model, method):
        """(bool, reason) Method 'exact' requires a product-form solution, the
        same rule the analyzer enforces at solve time. Order-independent and
        pass-and-swap stations are exempt: solver_mva_oi_analyzer is exact for
        them regardless of the product-form test. Single-station open systems
        are exempt too: they go to a queueing-system formula (M/G/1 PK, M/M/k,
        Cobham, matrix-geometric, ...) that holds outside product form, never to
        the MVA recursion. Product form has no registry
        feature name, so the check cannot live in getMethodFeatureSet. Mirrors
        MATLAB SolverMVA.supportsExactness."""
        if method != 'exact':
            return True, ''
        if model.hasProductFormSolution():
            return True, ''
        from ...api.sn.network_struct import SchedStrategy
        from ...lang.base import NodeType
        sn = model.getStruct()
        sched = getattr(sn, 'sched', None)
        if isinstance(sched, dict):
            for value in sched.values():
                if value in (SchedStrategy.OI, SchedStrategy.PAS):
                    return True, ''
        nodetype = getattr(sn, 'nodetype', None)
        if nodetype is not None and len(nodetype) == 3 and getattr(sn, 'nclosedjobs', 0) == 0:
            types = set(nodetype)
            if types in ({NodeType.SOURCE, NodeType.QUEUE, NodeType.SINK},
                         {NodeType.SOURCE, NodeType.CACHE, NodeType.SINK}):
                return True, ''
        return False, ("method 'exact' requires a product-form solution; use "
                       "'mva' for the approximation based on the exact MVA algorithm")

    supports_exactness = supportsExactness

    @staticmethod
    def supportsFiniteCapacity(model):
        """(bool, reason) MVA-specific finite-capacity gate: Blocking-After-Service
        models are exempt because MVA offers the Smith queue-decomposition method
        'sqd', and the analyzer routes a BAS model to solver_sqd under the default
        method too, so the finite buffers ARE honoured on every MVA path.
        Everything else defers to the shared product-form gate. Mirrors MATLAB
        SolverMVA.supportsFiniteCapacity."""
        if not hasattr(model, 'getStruct'):
            return True, ''
        from ...api.solvers.mva.analyzers import _is_bas_model
        if _is_bas_model(model.getStruct()):
            return True, ''
        return NetworkSolver.checkBindingCapacity(model, 'SolverMVA')

    supports_finite_capacity = supportsFiniteCapacity

    @staticmethod
    def getFeatureSet() -> set:
        """
        Get set of features supported by the MVA solver.

        Returns the canonical feature names (mirrors MATLAB
        SolverMVA.getFeatureSet and the JAR SolverMVA).
        """
        return {
            'Sink', 'Source',
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue',
            'APH', 'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp', 'BMAP',
            'Pareto', 'Weibull', 'Lognormal', 'Uniform', 'Det',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer', 'Buffer', 'Dispatcher',
            'CacheClassSwitcher', 'Cache', 'CacheRetrieval',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_FCFSPRPRIO',
            'SchedStrategy_DPS', 'SchedStrategy_FCFS', 'SchedStrategy_SIRO', 'SchedStrategy_HOL',
            'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR', 'SchedStrategy_POLLING',
            # exact order-independent path only (solver_mva_oi_analyzer)
            'SchedStrategy_OI', 'SchedStrategy_PAS',
            # size-based M/G/1 disciplines, served by _run_sizebased_analysis
            # (Wierman and Harchol-Balter, SIGMETRICS 2003)
            'SchedStrategy_SRPT', 'SchedStrategy_PSJF', 'SchedStrategy_FB',
            'SchedStrategy_LRPT', 'SchedStrategy_SETF',
            # closed models only (_run_sjn)
            'SchedStrategy_SJF',
            'Fork', 'Forker', 'Join', 'Joiner',
            'JoinPartial',  # quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO', 'ReplacementStrategy_LRU',
            'ReplacementStrategy_HLRU',
            'MMAP',  # marked MAP sources (cache LRU via cache_ttl_lrum_map)
            'ClosedClass', 'SelfLoopingClass', 'OpenClass', 'Replayer',
            'LoadDependence',
            'ClassDependence',
        'JointDependence',
            # c-server stations: the exact recursion, every AMVA kernel,
            # qna/rqt and the M/M/k and G/G/k closed forms carry the count;
            # getMethodFeatureSet withdraws it from the single-server names.
            # FiniteCapacity is deliberately NOT here (see getMethodFeatureSet
            # and supportsFiniteCapacity).
            'MultiServer',
        }

    @staticmethod
    def supports(model, extra_features=None) -> bool:
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
        from ..base import SolverFeatureSet
        try:
            # Registry-based inclusion check (MATLAB SolverMVA.supports):
            # every feature the model uses must be in the MVA supported set.
            if hasattr(model, 'get_used_lang_features') or hasattr(model, 'getUsedLangFeatures'):
                if hasattr(model, 'get_used_lang_features'):
                    feat_used = model.get_used_lang_features()
                else:
                    feat_used = model.getUsedLangFeatures()
                feat_supported = SolverFeatureSet()
                feat_supported.set_true(list(SolverMVA.getFeatureSet()))
                # RQNA extends the base MVA feature set with the MAP family only while the RQNA dispatch is active.
                if extra_features:
                    feat_supported.set_true(list(extra_features))
                if not SolverFeatureSet.supports(feat_supported, feat_used):
                    return False

                # finite station/class cap rejected structurally (no registry name); closed models where cap can't bind (>=population) exempt, as BAS ('sqd')/Cache.
                ok, reason = SolverMVA.supportsFiniteCapacity(model)
                if not ok:
                    line_warning('SolverMVA', reason)
                    return False
                return True

            # Fallback (NetworkStruct input): basic sanity check only
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
            - 'method': 'default' (auto-selects exact/amva, as MATLAB/JAR do)
            - 'tol': 1e-4 (general-purpose tolerance)
            - 'max_iter': 1000 (maximum iterations)
            - 'verbose': default_verbose() (inherits GlobalConstants verbosity)
            - 'config': {} (per-method switches, e.g. 'map_env_method')

        `config` is present but EMPTY, as MATLAB's `SolverMVA.defaultOptions`
        carries an empty config struct: every consumer reads it with a default,
        so an absent key and an unset one mean the same thing, and the attribute
        must exist for `options.config['key'] = ...` to work.

        Example:
            >>> opts = SolverMVA.defaultOptions()
            >>> opts['method'] = 'amva'  # Override for approximate MVA
            >>> solver = SolverMVA(model, **opts)
        """
        return OptionsDict({
            'method': 'default',
            'tol': 1e-4,
            'max_iter': 1000,
            'verbose': default_verbose(),
            'config': {},
        })

    # ============================================================================
    # Sampling and Transient Methods (Phase 6) - Placeholders
    # ============================================================================

    def sample(self, node: int, numEvents: int) -> np.ndarray:
        """
        Sample from the response time distribution.

        **Not supported by MVA solver** - MVA is an analytical solver.
        For sampling, use simulation-based solvers.

        Args:
            node: Node/station index (1-based)
            numEvents: Number of samples to generate

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
            "Use SolverSSA, SolverLDES, or SolverJMT for sampling-based analysis."
        )

    def sampleAggr(self, node: int, numEvents: int) -> np.ndarray:
        """Aggregate sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleAggr() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, LDES, JMT)."
        )

    def sampleSys(self, numEvents: int) -> np.ndarray:
        """System-level sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleSys() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, LDES, JMT)."
        )

    def sampleSysAggr(self, numEvents: int) -> np.ndarray:
        """Aggregate system-level sampling (not supported by MVA)."""
        raise NotImplementedError(
            "sampleSysAggr() not supported by SolverMVA. "
            "Use simulation-based solvers (SSA, LDES, JMT)."
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
            "Transient analysis available via SolverCTMC or SolverLDES."
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
    getAvgT = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable
    default_options = defaultOptions

    # Chain-level aliases
    GetAvg = NetworkSolver.getAvg
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
    avg_respt = getAvgRespT
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
