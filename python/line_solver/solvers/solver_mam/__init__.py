"""
SolverMAM - Main matrix-analytic methods solver.

Implements 8 solution methods for queueing networks:
1. dec.source - Decomposition with MMAP arrivals (default)
2. dec.mmap - Service-scaled departures
3. dec.poisson - Poisson approximation
4. mna - Matrix-normalizing approximation (auto-selects open/closed)
5. ldqbd - Level-dependent QBD
8. fj - Fork-Join (percentile analysis)

Usage:
    solver = SolverMAM(network, method='default')
    solver.runAnalyzer()
    QN = solver.getAvgQLen()
    RN = solver.getAvgRespT()
"""

import os
import numpy as np
import pandas as pd
import sys
import time
from typing import Any, Optional, Dict, List, Tuple
from dataclasses import dataclass, field
from ...constants import default_verbose
from ...api.io.logging import LineError

from .algorithms import (
    MAMResult,
    DecSourceAlgorithm,
    DecMMAPAlgorithm,
    DecPoissonAlgorithm,
    DecSourceMMAPAlgorithm,
    MNAOpenAlgorithm,
    MNAClosedAlgorithm,
    LDQBDAlgorithm,
    ldqbd_is_closed_delay_queue,
    BgchainAlgorithm,
)
from .utils import (extract_mam_params, check_closed_network, is_fork_join_network,
                    extract_visit_counts)
from .fj.validator import fj_is_homogeneous
from .fj.solver import FJSolver
from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn.getters import sn_get_arvr_from_tput
from ..base import NetworkSolver, method_label, method_type
from ..avg_results import AvgResultsMixin


# The RCAT family, which moved to SolverAG. Kept here only to redirect a caller
# that still asks for one of these names.
_RCAT_METHODS = ('inap', 'inapplus', 'inapinf', 'exact')


def _has_signal_class(sn) -> bool:
    """True when the model declares at least one G-network signal class."""
    issignal = getattr(sn, 'issignal', None)
    if issignal is None:
        return False
    issignal = np.asarray(issignal, dtype=float)
    return issignal.size > 0 and bool(np.any(issignal > 0))


def _has_setup_station(sn) -> bool:
    """True when some station declares a setup/delay-off server (``sn.hassetup``).

    Setup/delay-off is read by the dec.source decomposition and, since 2026-09,
    by LDQBD through ``qbd_setupdelayoff_closed``. The bgchain and mna analyzers
    still DROP it, so the default dispatch must not hand them such a model: they
    would return the always-warm answer and nothing would say so. The LDQBD
    branch does not consult this predicate, because
    ``ldqbd_is_closed_delay_queue`` admits exactly the regime that analysis
    covers and refuses the rest. Mirrors ``hasSetup`` in the MATLAB
    solver_mam_analyzer.
    """
    hassetup = getattr(sn, 'hassetup', None)
    if hassetup is None:
        return False
    arr = np.asarray(hassetup, dtype=float).ravel()
    return arr.size > 0 and bool(np.any(arr > 0))


def _has_nontrivial_lld(sn) -> bool:
    """True when some station declares a load-dependent scaling other than alpha == 1."""
    lld = getattr(sn, 'lldscaling', None)
    if lld is None:
        return False
    lld = np.asarray(lld, dtype=float)
    return lld.size > 0 and bool(np.any(lld != 1.0))


def mam_retrial_applicable(sn):
    """Can the 'retrial' method answer this model? Returns (bool, reason).

    THE RULE IS A "MUST BE PRESENT" ONE, which is why it cannot live in a
    feature set: a feature set says "I accept this construct", so it can refuse
    a model for HAVING something and never for LACKING it. solver_mam_retrial
    needs an impatience configuration to analyze -- either the BMAP/PH/N/N
    bufferless retrial topology of Dudin et al. (Mathematics 13(9), 2025) or a
    reneging patience law for the MAP/M/s+G analysis of Gursoy, Mehr and Akar
    -- and a model carrying neither is not a smaller retrial model, it is a
    different one.

    ONE PREDICATE, TWO CALLERS: ``supportsModelMethod`` asks it so the method is
    not offered by findSolver or SolverAUTO on a model it cannot answer, and
    ``_solve_retrial_reneging`` asks it so a caller naming 'retrial' by hand
    gets the identical sentence, which is also the one the api-level
    ``solver_mam_retrial`` raises on its own.
    """
    from ...api.qsys.retrial import qsys_is_retrial, has_reneging_patience
    try:
        is_retrial, ret_info = qsys_is_retrial(sn)
    except Exception:
        is_retrial, ret_info = False, None
    if is_retrial:
        return True, ''
    try:
        if has_reneging_patience(sn):
            return True, ''
    except Exception:
        pass
    # qsys_is_retrial reports WHICH requirement the model missed (open model,
    # single class, a bufferless station, a retrial drop rule); carrying it
    # through is the difference between a bare no and a usable answer. The
    # wording is solver_mam_retrial's own, so the gate and the run agree.
    detail = getattr(ret_info, 'error_msg', '') or ''
    reason = 'No valid impatience configuration detected (retrial or reneging).'
    if detail:
        reason = reason + ' ' + detail
    return False, reason


def mam_buffer_refusal(sn):
    """Why no MAM method can answer a model whose finite buffer a CLOSED class
    can fill, or '' when every binding buffer is reached by open classes only.

    A MAM analyzer represents a finite buffer as a LOSS buffer: the basic
    analyzer solves an M/M/c/K or an MMAP[K]/G/1/K, and the decomposition and
    open network-analysis routes truncate and renormalize the same way. That is
    the right model for an OPEN class, whose refused arrival is lost. A closed
    job that finds no room BLOCKS instead, and no MAM analyzer blocks: the loss
    formulas answer a different system, and the closed routes of 'default'
    (ldqbd, bgchain, the closed arm of mna) read no sn.cap or sn.classcap at
    all. So the pair is refused rather than answered.

    ONE PREDICATE, TWO CALLERS: ``supportsModelMethod`` reports it and
    ``runAnalyzer`` raises it, so a caller gets one answer whichever it meets
    first. Only a buffer that can BIND counts, which is what sn_get_buffer_size
    decides: refresh_capacity derives a finite classcap (the chain population)
    at every station of every closed model, so a plain finiteness test would
    refuse every closed model. A Cache builds its own capped retrieval queues
    and is exempt, as in sn_has_blocking. Port of MATLAB mam_buffer_refusal.
    """
    from ...api.me.solver_nc_mem import sn_get_buffer_size
    from ...api.sn.network_struct import NodeType
    nodetype = getattr(sn, 'nodetype', None)
    if nodetype is not None and any(int(t) == int(NodeType.CACHE) for t in nodetype):
        return ''
    classcap = getattr(sn, 'classcap', None)
    if classcap is None:
        return ''
    classcap = np.atleast_2d(np.asarray(classcap, dtype=float))
    if classcap.size == 0:
        return ''
    closed = np.isfinite(np.asarray(sn.njobs, dtype=float).ravel())
    for ist in range(int(sn.nstations)):
        if not np.isfinite(sn_get_buffer_size(sn, ist)):
            continue
        if ist >= classcap.shape[0]:
            continue
        served = classcap[ist, :] > 0
        n = min(served.size, closed.size)
        idx = [r for r in range(n) if closed[r] and served[r]]
        if not idx:
            continue
        r = idx[0]
        return ('Station %s carries a finite capacity that binds for the closed class %s. '
                'SolverMAM represents a finite buffer as a LOSS buffer (M/M/c/K, '
                'MMAP[K]/G/1/K, truncate-and-renormalize), which is the open-class model: a '
                'closed job that finds no room blocks instead, and no MAM analyzer blocks, '
                'while the closed routes of the default method (ldqbd, bgchain, mna) read no '
                'capacity at all. Use SolverCTMC, SolverSSA or SolverLDES, or SolverMVA with '
                "method 'sqd' for blocking after service."
                % (sn.nodenames[int(sn.stationToNode[ist])], sn.classnames[r]))
    return ''


@dataclass
class SolverMAMOptions:
    """Options for SolverMAM.

    Attributes:
        method: Algorithm to use (default, dec.source, mna, ldqbd, bgchain)
        tol: Convergence tolerance
        max_iter: Maximum iterations
        space_max: Maximum MMAP state space size
        verbose: Print debug information
    """
    method: str = 'default'
    tol: float = 1e-4
    max_iter: int = 100
    space_max: int = 1000
    # Per-method knobs, the twin of MATLAB options.config and the JAR
    # options.config map: 'bgaggr' and 'bgstates_max' size the background chain
    # of bgchain, 'space_max'/'qbdphases_max' bound its MMAP superposition and
    # its QBD. Read through options.config.get(...) by the algorithms, so an
    # absent key keeps that algorithm's own default.
    config: Dict[str, Any] = field(default_factory=dict)
    verbose: bool = field(default_factory=default_verbose)
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))  # env LINE_SOLVER_LANG overrides; 'python' (native), 'java' (jline.jar via JSON) or 'cpp' (line-cli via JSON)
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. Meaningless for the other langs, which are IEEE double
    # throughout, so line-cli is invoked without --arith unless the caller sets it.
    # The C++ MAM analyzer fits phase-type representations, so it refuses anything
    # but double: 'exact' is rejected by name rather than narrowed.
    arith: Optional[str] = None


class SolverMAM(AvgResultsMixin, NetworkSolver):
    """Native Python solver for matrix-analytic methods.

    Solves queueing networks using decomposition, MNA, RCAT, and related methods.
    """

    # Available methods and their algorithm classes
    ALGORITHMS = {
        'dec.source': DecSourceAlgorithm,
        'dec.mmap': DecMMAPAlgorithm,
        'dec.poisson': DecPoissonAlgorithm,
        'dec.source.mmap': DecSourceMMAPAlgorithm,
        'mna': None,  # Auto-select based on network type
        'mna_open': MNAOpenAlgorithm,
        'mna_closed': MNAClosedAlgorithm,
        'ldqbd': LDQBDAlgorithm,
        'bgchain': BgchainAlgorithm,
    }

    def __init__(self, network, method: str = 'default', options: Optional[SolverMAMOptions] = None, **kwargs):
        """Initialize SolverMAM.

        Args:
            network: Network model (must be compiled to NetworkStruct)
            method: Solution method ('default', 'dec.source', 'mna', etc.)
            options: SolverMAMOptions instance
            **kwargs: Additional parameters (verbose, seed, etc.) for compatibility
        """
        self.network = network
        # The shared feature gate in NetworkSolver, and the solver console,
        # both read self.model (as CTMC, MVA, NC and FLD do). Without this
        # alias the gate never runs and the console reports '(unnamed)'.
        self.model = network
        self.sn = self._get_network_struct(network)
        # Store seed if provided (for compatibility, though MAM is analytical)
        self._seed = kwargs.get('seed', None)

        if options is None:
            options = SolverMAMOptions(method=method)
        else:
            if method != 'default':
                options.method = method

        # Handle verbose kwarg
        if 'verbose' in kwargs:
            options.verbose = kwargs['verbose']
        # Opt-in JAR delegation (lang='java'); default stays native/JVM-free.
        if 'lang' in kwargs:
            options.lang = kwargs['lang']
        # Carry timespan/cutoff (used by the SolverENV state-vector analyzer's
        # MAM/LDQBD backend; harmless for the standard MAM path).
        if 'timespan' in kwargs:
            options.timespan = kwargs['timespan']
        if 'cutoff' in kwargs:
            options.cutoff = kwargs['cutoff']

        self.options = options
        self.result = None
        self.runtime = 0.0

    def reset(self):
        """Reset the solver to force recomputation on next getAvg call."""
        self._clearResultStores()
        # Re-read network struct since model may have been updated by LN iteration
        self.sn = self._get_network_struct(self.network)

    def getName(self) -> str:
        """Get the name of this solver."""
        return "MAM"

    get_name = getName

    def _get_network_struct(self, model):
        """Get NetworkStruct from model using priority-based extraction."""
        sn = None

        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            sn = model._sn
        # Priority 2: Native model with refresh_struct()
        elif hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                sn = model._sn
        # Priority 3: Native model with snake-case get_struct() (no wrapper
        # bridge — native solvers reject JAR-wrapper models, keeping python/
        # free of any JAR/JVM coupling).
        elif hasattr(model, 'get_struct'):
            sn = model.get_struct()
        # Priority 4: Model that is already a struct
        elif hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            sn = model

        if sn is None:
            raise ValueError("Cannot extract network structure from model")

        # Check for SetupTask params in model.attribute (set by SolverLN)
        # and propagate to sn.hassetup and sn.nodeparam
        if hasattr(model, 'attribute') and model.attribute is not None:
            attr = model.attribute
            if hasattr(attr, 'get'):
                func_params = attr.get('functionParams', None)
            elif isinstance(attr, dict):
                func_params = attr.get('functionParams', None)
            else:
                func_params = getattr(attr, 'functionParams', None)

            if func_params is not None:
                # Set hassetup for the server station
                server_idx_1based = func_params.get('serverIdx', 1)
                # Convert to 0-indexed station index
                # serverIdx is the node index in the layer model (1-indexed)
                # In a layer model with Clients (idx=1) and Server (idx=2), the station indices are:
                # - Station 0: Clients (Delay)
                # - Station 1: Server (Queue)
                server_station_idx = server_idx_1based - 1

                # Initialize hassetup if not present
                if not hasattr(sn, 'hassetup') or sn.hassetup is None:
                    sn.hassetup = np.zeros(sn.nstations)
                elif len(sn.hassetup) < sn.nstations:
                    sn.hassetup = np.zeros(sn.nstations)

                if server_station_idx < len(sn.hassetup):
                    sn.hassetup[server_station_idx] = 1

                # Set nodeparam for the server station
                if not hasattr(sn, 'nodeparam') or sn.nodeparam is None:
                    sn.nodeparam = {}

                sn.nodeparam[server_station_idx] = func_params

        return sn

    def runAnalyzer(self) -> 'SolverMAM':
        """Run the analyzer with selected method.

        Returns:
            self (for method chaining)
        """
        # Opt-in delegation to the canonical JAR (mirrors MATLAB options.lang='java').
        # Populates the native result container from jline.jar so every getter
        # (tables, matrices, chain/node/scalar metrics) returns JAR-derived values.
        # Imported lazily so a JVM-free install never touches this path.
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
                import warnings
                warnings.warn("SolverMAM: lang='cpp' requested but the C++ solver is "
                              "unavailable (%s); falling back to lang='python'." % e)

        # The RCAT methods moved to SolverAG. Name them rather than letting the
        # dispatch report an unknown method, so a caller carrying an old
        # options.method is told where they went. runAnalyzerChecks asks the
        # same helper BEFORE reporting an unlisted method, so whichever gate a
        # caller meets first, the answer is this one. Kept here for the
        # enableChecks=False path, which skips that gate entirely.
        moved = self.unsupportedMethodReason(getattr(self.options, 'method', 'default'))
        if moved:
            raise RuntimeError(moved)

        # G-network signals belong to the RCAT analyzer alone: no MAM algorithm
        # reads sn.issignal, so every one of them would solve the model with the
        # signals turned into ordinary customers and report that as the answer.
        if getattr(self, 'enableChecks', True) and _has_signal_class(self.sn):
            raise RuntimeError(
                'The %s method does not support G-network signals: no MAM algorithm reads '
                'sn.issignal, so this method would solve the model with every signal turned '
                'into an ordinary customer. Use SolverAG, whose RCAT methods are the only '
                'ones that model signals.' % self._select_method())

        # Method-aware feature-set gate (MATLAB @SolverMAM/runAnalyzer.m line 17):
        # reject models using features outside the MAM feature set (e.g.
        # LCFSPR scheduling) instead of silently mishandling them. Runs
        # before any dispatch, including the mapmap1 exact fast path.
        if getattr(self, 'enableChecks', True):
            model = getattr(self, 'network', None) or getattr(self, 'model', None)
            feat_used = None
            if model is not None:
                if hasattr(model, 'get_used_lang_features'):
                    feat_used = model.get_used_lang_features()
                elif hasattr(model, 'getUsedLangFeatures'):
                    feat_used = model.getUsedLangFeatures()
            if feat_used is not None:
                from ..base import SolverFeatureSet
                feat_supported = SolverFeatureSet()
                feat_supported.set_true(list(self.getMethodFeatureSet(getattr(self.options, 'method', 'default'))))
                ok, reason = SolverFeatureSet.supports_with_reason(feat_supported, feat_used)
                if not ok:
                    raise LineError('This model contains features not supported by the solver. '
                                    + (reason or ''))

        # A binding buffer a CLOSED class can fill: every MAM analyzer models a
        # finite buffer as a LOSS buffer, which is the open-class system, and
        # the closed routes read no capacity at all. mam_buffer_refusal is the
        # same predicate supportsModelMethod reports, so the gate and the run
        # give one sentence.
        if getattr(self, 'enableChecks', True):
            buffer_reason = mam_buffer_refusal(self.sn)
            if buffer_reason:
                raise LineError(buffer_reason)

        # Finite Capacity Region: MAM does not enforce the aggregate per-region
        # job limit and would silently return the unconstrained answer.
        if getattr(self.sn, 'nregions', 0) > 0:
            raise RuntimeError('This model uses a Finite Capacity Region (addRegion), which is '
                               'not supported by SolverMAM (the region\'s aggregate job limit is '
                               'not enforced). Use SolverCTMC, SolverJMT, SolverSSA or SolverLDES, '
                               'or setCapacity for a single-station limit.')

        method = self._select_method()


        start_time = time.time()

        # Discrete-time (slotted) models are recognized from the distributions
        # and routed to the Q-MAM discrete-time algorithms. The test must run
        # BEFORE sn_nonmarkov_toph below, which would fit a continuous PH to a
        # Geometric and erase the lattice; see _kb/06-solver-catalog.md for the
        # LAS-DA convention.
        from ...api.sn.predicates import sn_is_discrete_time
        is_dt, slot_length, _dt_info = sn_is_discrete_time(self.sn, self.options)
        if is_dt:
            from ...api.solvers.mam.dt import solver_mam_dt
            dt_ret = solver_mam_dt(self.sn, self.options, slot_length)
            self.result = MAMResult(QN=dt_ret.QN, UN=dt_ret.UN, RN=dt_ret.RN, TN=dt_ret.TN,
                                    CN=dt_ret.CN, XN=dt_ret.XN, totiter=dt_ret.totiter,
                                    method=dt_ret.method, runtime=time.time() - start_time)
            return self

        # Exact fast-path: single-class single-server open MAP/MAP/1 with a
        # correlated (non-renewal) MAP arrival or service, which the
        # decomposition methods only approximate. Uses the raw MAP blocks.
        from .algorithms.mapmap1_exact import solver_mam_mapmap1_exact
        self.result = solver_mam_mapmap1_exact(self.sn)

        if self.result is None:
            # Convert non-Markovian distributions (Det, Gamma, Weibull, Lognormal,
            # Pareto, Uniform) to PH before analysis. Matches MATLAB
            # solver_mam_analyzer.m:13-18: preserveDet=true keeps Det for the
            # exact MAP/D/c (Crommelin) dispatch in solver_mam_basic.
            from ...api.sn import sn_nonmarkov_toph
            opts_dict = {}
            if hasattr(self.options, '__dict__'):
                opts_dict = {k: v for k, v in self.options.__dict__.items() if not k.startswith('_')}
            elif isinstance(self.options, dict):
                opts_dict = dict(self.options)
            config = dict(opts_dict.get('config', {}) or {})
            # The RCAT methods build a CTMC per component out of (D0,D1), so a
            # preserved Det would reach them with no matrix at all and be read
            # back as its mean rate, and a concentrated matrix exponential is
            # not a generator at all. They need a genuine phase-type, as SSA,
            # Fluid and JMT do.
            config.setdefault('preserveDet', True)
            opts_dict['config'] = config
            # The conversion RETAGS procid to APH/ME/MAP, so afterwards nothing
            # names the law the user declared. MMAP[K]/G[K]/1 works off that law's
            # transform (sn.lst) rather than off the surrogate, so its gate needs
            # the tags as they stood here; only DET survives the retagging.
            if getattr(self.sn, 'procid', None) is not None:
                self.sn.procid_declared = np.array(self.sn.procid, dtype=object).copy()
            self.sn = sn_nonmarkov_toph(self.sn, opts_dict)

            # SetupTask stations are NOT dispatched to a dedicated solver:
            # MATLAB routes them through solver_mam_basic like any other model,
            # where the FCFS branch handles sn.hassetup stations and the
            # post-loop population wash makes the setup/delay-off queue length
            # inert (RN = S). A separate front-end solver here shadowed that
            # path and pinned LN host layers at a spurious throughput.
            if method == 'ldqbd':
                self.result = self._solve_ldqbd()
            elif method == 'fj':
                # Special handling for Fork-Join
                self.result = self._solve_fork_join()
            elif method in ('retrial', 'reneging'):
                # Retrial or reneging solver dispatch
                self.result = self._solve_retrial_reneging(method)
            else:
                # Standard algorithm dispatch
                algo_class = self.ALGORITHMS.get(method)
                if algo_class is None:
                    raise ValueError(f"Unknown method: {method}")

                # Method-aware feature gate: each decomposition algorithm
                # declares a structural applicability predicate. Gate here (the
                # standard dispatch branch) so an unsupported model is rejected
                # with a precise reason instead of silently mishandled; the
                # special/auto methods (fj, ldqbd, retrial, reneging, mna) are
                # dispatched above and keep their own handling.
                if getattr(self, 'enableChecks', True):
                    ok, reason = algo_class.supports_network(self.sn)
                    if not ok:
                        raise RuntimeError(
                            "This model contains features not supported by the "
                            "MAM solver's '%s' method. %s" % (method, reason or ''))

                algo = algo_class()
                self.result = algo.solve(self.sn, self.options)

                # The closed MNA outer bisection rescales each chain onto its
                # population as a last step, so a diverged run still returns
                # queue lengths that sum to N. Little's law on the unrescaled
                # R and T is what still reads the raw iterate; when it fails,
                # fall back to dec.source rather than report the collapse.
                if (str(self.options.method).lower() == 'default'
                        and method == 'mna_closed'
                        and not self._mna_conserves(self.result)):
                    dec = self.ALGORITHMS.get('dec.source')
                    self.result = dec().solve(self.sn, self.options)
                    method = 'dec.source'

        # Set Source station TN to arrival rates (matching MATLAB solver_mam_analyzer.m lines 135-140)
        if self.result is not None and hasattr(self.result, 'TN') and self.result.TN is not None:
            from ...lang.base import SchedStrategy as SchedStrategyBase
            for i in range(self.sn.nstations):
                sched_i = self.sn.sched.get(i, None) if isinstance(self.sn.sched, dict) else (self.sn.sched[i] if i < len(self.sn.sched) else None)
                if sched_i is not None:
                    sched_name = sched_i.name if hasattr(sched_i, 'name') else str(sched_i)
                    sched_val = sched_i.value if hasattr(sched_i, 'value') else int(sched_i)
                    if sched_name == 'EXT' or sched_val == 16:
                        if i < self.result.TN.shape[0] and hasattr(self.sn, 'rates') and self.sn.rates is not None:
                            rates = np.asarray(self.sn.rates)
                            if i < rates.shape[0]:
                                self.result.TN[i, :] = rates[i, :]

        # Compute proper residence times from response times (WN = RN * visits / ref_visits)
        if self.result is not None and hasattr(self.result, 'RN') and self.result.RN is not None:
            self.result.WN = sn_get_residt_from_respt(self.sn, self.result.RN, None)

        # Compute proper arrival rates from throughputs using routing
        if self.result is not None and hasattr(self.result, 'TN') and self.result.TN is not None:
            self.result.AN = sn_get_arvr_from_tput(self.sn, self.result.TN)

        self.runtime = time.time() - start_time

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            iter_count = self.result.totiter if hasattr(self.result, 'totiter') else 1
            if iter_count <= 1:
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"MAM analysis [method: {method_label(self.options.method, method)}; type: {method_type('MAM', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {self.runtime:.6f}s.")
            else:
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"MAM analysis [method: {method_label(self.options.method, method)}; type: {method_type('MAM', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {self.runtime:.6f}s. Iterations: {iter_count}.")

        return self

    def _select_method(self) -> str:
        """Select method based on network type if method='default'.

        Matches MATLAB solver_mam_analyzer.m routing logic:
        1. Fork-Join topology -> 'fj'
        2. BMAP/PH/N/N retrial topology -> 'retrial'
        3. MAP/M/s+G reneging topology -> 'reneging'
        4. Single-class closed Delay+Queue -> 'ldqbd'
        5. Default -> 'dec.source'

        Returns:
            Selected method name
        """
        method = self.options.method

        # Fork-Join: MATLAB solver_mam_analyzer routes both 'default' and an
        # explicit 'dec.source' on a Fork-Join topology to the FJ solver, rather
        # than rejecting. Align to that ground truth here.
        if method in ('default', 'dec.source') and is_fork_join_network(self.sn):
            return 'fj'

        if method == 'default':
            # Auto-selection logic (Fork-Join already handled above)
            # Check for retrial/reneging topologies before defaulting
            from ...api.qsys.retrial import qsys_is_retrial, has_reneging_patience
            try:
                is_retrial, _ = qsys_is_retrial(self.sn)
                if is_retrial:
                    return 'retrial'
            except Exception:
                pass

            try:
                if has_reneging_patience(self.sn):
                    return 'reneging'
            except Exception:
                pass

            # Single-class closed Delay+Queue: the LD-QBD is exact (the
            # level-dependent arrival rate (N-n)*lambda captures the population
            # constraint that dec.source only approximates).
            if ldqbd_is_closed_delay_queue(self.sn):
                return 'ldqbd'

            # A closed model is the degenerate case of the background chain: with
            # no open work to take a share of the servers the chain is the EXACT
            # closed CTMC at chain granularity, so it dominates the mna fixed
            # point wherever its state space fits. _bgchain_affordable sizes that
            # state space and _bgchain_closed_exact checks the chain is built
            # from a service law it represents exactly.
            if check_closed_network(self.sn) and not _has_setup_station(self.sn) \
                    and BgchainAlgorithm.supports_network(self.sn)[0] \
                    and self._bgchain_affordable() and self._bgchain_closed_exact():
                return 'bgchain'

            # A closed model has no arrival stream for dec.source to build its
            # Poisson surrogate from: it replaces each closed chain by a source
            # at the current throughput iterate and never enforces the
            # population, so the answer neither conserves N nor separates the
            # classes. MNA closes the same traffic equations by bisecting the
            # per-class throughput against N.
            if check_closed_network(self.sn) and not _has_setup_station(self.sn) \
                    and self._mna_applies():
                return 'mna_closed'

            # Mixed: the closed classes are solved exactly as a background chain
            # and the open ones as QBDs driven by it, which is 4-5 significant
            # digits against CTMC where dec.source is 10-24% out.
            njobs = np.asarray(self.sn.njobs, dtype=float).ravel()
            is_mixed = bool(np.any(np.isinf(njobs)) and np.any(np.isfinite(njobs)))
            if is_mixed and not _has_setup_station(self.sn) \
                    and BgchainAlgorithm.supports_network(self.sn)[0] \
                    and self._bgchain_affordable():
                return 'bgchain'

            # Default to dec.source
            return 'dec.source'
        elif method == 'mna':
            # Auto-select between mna_open and mna_closed
            is_closed = check_closed_network(self.sn)
            return 'mna_closed' if is_closed else 'mna_open'
        else:
            return method

    def _bgchain_affordable(self) -> bool:
        """Whether the background chain this model needs is small enough to build.

        Mirrors the size clause of MATLAB bgchainApplies. The background chain
        enumerates the closed population vector over the stations the closed
        classes visit, so a large closed population makes bgchain_ctmc refuse
        the model outright. That refusal is right when the user asked for
        bgchain by name and wrong as a default, which must land on a method that
        answers: size the chain first and leave those models to dec.source. The
        limit is the one bgchain_ctmc enforces.
        """
        from .algorithms.bgchain import bgchain_states

        config = getattr(self.options, 'config', None) or {}
        if not isinstance(config, dict):
            config = dict(getattr(config, '__dict__', {}))
        bgstates_max = int(config.get('bgstates_max', 20000))
        return bgchain_states(self.sn, config) <= bgstates_max

    def _bgchain_closed_exact(self) -> bool:
        """Whether the background chain represents this closed model's service
        laws exactly.

        Mirrors MATLAB bgchainClosedExact in solver_mam_analyzer.m. A station that
        is not PS or INF must satisfy BOTH conditions below, because the
        background chain makes two separate first-moment substitutions there.

        1. bgchain_ctmc builds its generator from the MEAN service time alone.
           That is exact at a PS or INF station, which is insensitive to the
           service law beyond its first moment, and exact under any discipline
           when the law IS exponential. Measured on a closed Delay+FCFS cycle
           with Erlang-3 service, the mean-only chain reads 2.8% off SolverCTMC.
        2. The capacity a station's closed jobs hold is split over the background
           classes in proportion to their COUNTS, which is service in random
           order. That is exact under PS, and exact under FCFS only when the
           classes are served at the SAME rate -- an FCFS station with
           class-dependent rates reads 25.2% off SolverCTMC on a two-chain closed
           cycle, against 3.5e-16 when the two rates are made equal.

        mna_closed carries the phase-type representation instead, so neither
        surrogate may be chosen as the DEFAULT. Asking for bgchain by name still
        gets it, with both approximations documented in BgchainAlgorithm.

        Compare the ProcessType BY NAME: the enum's numeric values differ across
        codebases (MATLAB EXP=0, Python EXP=1).
        """
        from ...constants import GlobalConstants, ProcessType
        from .algorithms.bgchain import _sched_name

        sn = self.sn
        procid = getattr(sn, 'procid', None)
        if procid is None:
            return False
        rates = np.asarray(getattr(sn, 'rates', None), dtype=float)
        for ist in range(int(sn.nstations)):
            if _sched_name(sn, ist) in ('INF', 'PS', 'EXT'):
                continue
            rate_here = None
            for r in range(int(sn.nclasses)):
                if rates.size and (not np.isfinite(rates[ist, r]) or rates[ist, r] <= 0):
                    continue
                try:
                    p = procid[ist][r] if not hasattr(procid, 'shape') else procid[ist, r]
                except (IndexError, TypeError, KeyError):
                    return False
                if p is None:
                    return False
                pname = p.name if hasattr(p, 'name') else str(p)
                if pname != ProcessType.EXP.name:
                    return False
                if rate_here is None:
                    rate_here = float(rates[ist, r])
                elif abs(float(rates[ist, r]) - rate_here) > GlobalConstants.CoarseTol * rate_here:
                    return False
        return True

    def _mna_applies(self) -> bool:
        """Whether the closed MNA analyzer covers this model.

        Mirrors MATLAB mnaApplies in solver_mam_analyzer.m: the round-robin
        split is carried by the open traffic equations only, a self-looping
        class has no inter-station flow to decompose, and a station whose
        discipline the flow sweep does not update would keep a zero queue
        length. Any of those routes the default to dec.source instead.
        """
        from ...constants import GlobalConstants

        # solver_mna_closed drives its bisection over CLASSES but stores the
        # throughput in the CHAIN-indexed lambda, and renormalizes chain c's
        # queue lengths with the class-indexed njobs(c). Both are only correct
        # when each chain holds exactly one class, so a class-switching model
        # (fewer chains than classes) belongs to dec.source.
        if int(self.sn.nchains) != int(self.sn.nclasses):
            return False

        routing = getattr(self.sn, 'routing', None)
        if routing is not None:
            # sn.routing holds RoutingStrategy members, and LINE keeps distinct
            # enum copies whose members compare unequal, so test by name
            for r in np.asarray(routing, dtype=object).ravel():
                if getattr(r, 'name', None) == 'RROBIN':
                    return False

        sched = self.sn.sched if self.sn.sched else {}
        nservers = np.asarray(self.sn.nservers, dtype=float).ravel()
        for ist in range(int(self.sn.nstations)):
            s = sched.get(ist) if isinstance(sched, dict) else sched[ist]
            sname = getattr(s, 'name', None)
            if sname not in ('INF', 'PS', 'FCFS', 'EXT'):
                return False
            # The PS branch of solver_mna_closed forms U = S*T and the geometric
            # bound from it WITHOUT dividing by the number of servers, so a
            # multiserver PS station is misrepresented; dec.source is exact on
            # the non-queueing regime (c >> N) that shape usually stands for.
            if sname == 'PS' and nservers[ist] > 1:
                return False
            # A multiclass FCFS station makes the flow sweep superpose one MMAP
            # per class and then solve MMAPPH1FCFS at level sum(N)+1: measured
            # against an exact CTMC that costs 12-61s where dec.source costs
            # 0.03s and is not more accurate (mean relative error 0.11-0.38
            # against 0.04-0.19). The single-class case is both cheap and
            # better, so keep only that one.
            if sname == 'FCFS' and int(self.sn.nclasses) > 1:
                return False

        V = extract_visit_counts(self.sn)
        njobs = np.asarray(self.sn.njobs, dtype=float).ravel()
        for k in range(int(self.sn.nclasses)):
            if not np.isfinite(njobs[k]):
                continue
            vis = np.flatnonzero(V[:, k] > GlobalConstants.FineTol)
            if vis.size == 1:
                s = sched.get(int(vis[0])) if isinstance(sched, dict) else sched[int(vis[0])]
                if getattr(s, 'name', None) not in ('INF', 'EXT'):
                    return False

        return True

    def _mna_conserves(self, result) -> bool:
        """Whether the closed MNA outer bisection closed on N.

        Mirrors MATLAB mnaConserves in solver_mam_analyzer.m: the closed MNA
        algorithm rescales each chain onto its population as a last step, so a
        diverged bisection still returns queue lengths that sum to N and the
        failure is invisible in QN. R and T are NOT rescaled, so Little's law
        over the whole network, sum_i T(i,k)*R(i,k) = N_k, still reads the raw
        iterate: a converged run lands within 5e-4 of N and a diverged one is
        orders of magnitude out, or negative.
        """
        if result is None or result.QN is None or result.RN is None or result.TN is None:
            return False
        QN = np.asarray(result.QN, dtype=float)
        RN = np.asarray(result.RN, dtype=float)
        TN = np.asarray(result.TN, dtype=float)
        if not (np.all(np.isfinite(QN)) and np.all(np.isfinite(RN)) and np.all(np.isfinite(TN))):
            return False
        npred = np.sum(TN * RN, axis=0)
        njobs = np.asarray(self.sn.njobs, dtype=float).ravel()
        for k in range(int(self.sn.nclasses)):
            if not np.isfinite(njobs[k]) or njobs[k] <= 0:
                continue
            if abs(npred[k] - njobs[k]) > 0.01 * njobs[k]:
                return False
        return True

    def _solve_ldqbd(self) -> MAMResult:
        """Solve using LDQBD method.

        Returns:
            MAMResult
        """
        return LDQBDAlgorithm().solve(self.sn, self.options)

    def _solve_fork_join(self) -> MAMResult:
        """Solve Fork-Join network using FJ_codes.

        Returns:
            MAMResult with percentile response times attached
        """
        # Use FJ solver for topology validation and percentile computation
        fj_solver = FJSolver(verbose=self.options.verbose)

        # Validate Fork-Join topology
        can_solve, reason = fj_solver.can_solve(self.sn)
        if not can_solve:
            raise ValueError(f"Not a valid Fork-Join network: {reason}")

        # First, solve using dec.source to get basic metrics
        dec_algo = DecSourceAlgorithm()
        result = dec_algo.solve(self.sn, self.options)

        # Compute percentiles for Fork-Join
        percentiles = [50, 75, 90, 95, 99]
        mean_rt = np.sum(result.RN[:, 0]) if result.RN.size > 0 else 1.0

        fj_result = fj_solver.compute_percentiles(self.sn, percentiles, mean_rt)

        if fj_result is not None:
            # Attach percentile results to main result
            result.percentile_results = {
                'percentiles': fj_result.percentiles,
                'response_times': fj_result.response_times,
                'mean_response_time': fj_result.mean_response_time,
                'K': fj_result.K,
            }

        result.method = 'fj'
        return result

    def _solve_retrial_reneging(self, method: str) -> MAMResult:
        """Solve retrial or reneging queue using dedicated solvers.

        Dispatches to solver_mam_retrial which handles both:
        1. BMAP/PH/N/N bufferless retrial queues
        2. MAP/M/s+G queues with reneging (MAPMsG)

        Matches MATLAB solver_mam_analyzer.m lines 72-93.

        Args:
            method: 'retrial' or 'reneging'

        Returns:
            MAMResult with performance metrics
        """
        from ...api.qsys.retrial import solver_mam_retrial

        # The predicate supportsModelMethod asks, so the gate that decides
        # whether to OFFER 'retrial' and this run cannot drift apart.
        ok, why = mam_retrial_applicable(self.sn)
        if not ok:
            raise ValueError(why)

        # Build options dict from SolverMAMOptions
        opts = {
            'iter_max': self.options.max_iter,
            'tol': self.options.tol,
            'verbose': self.options.verbose,
        }

        QN, UN, RN, TN, CN, XN, totiter, _perf = solver_mam_retrial(self.sn, opts)

        # TN from retrial solver is (M, K), MAMResult expects (M, K) for TN
        # but standard MAM uses (1, K) for system throughputs
        # Keep the full station-level TN for consistency with other MAM methods

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            totiter=totiter,
            method=method,
            runtime=0.0,
        )

    # =====================================================================
    # RESULT ACCESS METHODS (following SolverMVA pattern)
    # =====================================================================

    @staticmethod
    def listValidMethods() -> List[str]:
        """List all valid solution methods.

        Returns:
            List of method names
        """
        # SolverMAM.m verbatim, in its order, MINUS the RCAT names ('inap',
        # 'inapplus', 'inapinf', 'exact') which moved to SolverAG. 'mna_open',
        # 'mna_closed' and 'fj' are INTERNAL resolution targets of the topology
        # router, not names a caller selects -- the reference resolves 'mna' to
        # the open or closed arm itself and reaches the fork-join route from
        # 'default' -- so they are dispatched but no longer advertised, exactly
        # as the reference does not advertise its own 'qiu' and 'reneging'.
        return ['default', 'dec.source', 'dec.mmap', 'dec.poisson', 'mna',
                'ldqbd', 'dec.source.mmap', 'bgchain', 'retrial']

    @staticmethod
    def supports(sn, method: str) -> Tuple[bool, Optional[str]]:
        """Check if method can solve this network.

        Args:
            sn: NetworkStruct
            method: Method name

        Returns:
            (can_solve, reason_if_not)
        """
        algo_class = SolverMAM.ALGORITHMS.get(method)
        if algo_class is None:
            # 'default', 'exact' and 'retrial' name analyzers the topology
            # router reaches rather than entries of ALGORITHMS: 'default' is the
            # router itself, 'exact' is the RCAT alias (grouped with
            # inap/inapplus/inapinf by solver_mam_analyzer.m) and 'retrial' is
            # the BMAP/PH/N/N analyzer. All three are advertised by
            # listValidMethods, so this must not call them unknown; their
            # structural applicability is decided by supportsModelMethod.
            if method in ('default', 'retrial') or method in _RCAT_METHODS:
                return True, None
            return False, f"Unknown method: {method}"

        return algo_class.supports_network(sn)

    def resolveMethod(self, options):
        """Feature-driven resolution of method='default' via the existing MAM
        topology router (_select_method). Part of the base NetworkSolver
        method-aware gating contract."""
        return self._select_method()

    def getMethodFeatureSet(self, method):
        """Per-method feature deltas on the base MAM envelope. Only 'mna'
        resolves a round-robin split (npfqn_traffic_split_rr in solver_mna_open);
        the closed branch has no counterpart and is rejected in
        supportsModelMethod. Mirrors the MATLAB/JAR
        SolverMAM.getMethodFeatureSet."""
        feats = set(SolverMAM.getFeatureSet())
        if method == 'mna':
            feats.add('RoutingStrategy_RROBIN')
        if method == 'dec.mmap':
            # dec.mmap is an OPEN-network departure-process fixed point: it
            # iterates on arrival streams a closed population does not have,
            # and its station ladder serves EXT, FCFS, HOL, FCFSPRPRIO and PS
            # only. Both restrictions are things the model HAS, so both belong
            # here rather than in supportsModelMethod; the handler raises the
            # matching message when the method is named by hand. Until this
            # delta existed the gate offered dec.mmap on every closed model and
            # the call died in MAMResult with four missing arguments, because
            # the handler had returned an empty result the caller unpacked.
            # Fork-join goes too: the sweep uses the plain traffic step, which
            # has no synchronization -- which is why the topology router sends a
            # fork-join model from 'default'/'dec.source' to the FJ solver and
            # never here. Left declared, the gate offered dec.mmap on an open
            # fork-join model and the departure process it then built had no
            # recurrent state.
            feats -= {'ClosedClass', 'SelfLoopingClass', 'SchedStrategy_INF',
                      'Fork', 'Join', 'Forker', 'Joiner'}
        # ldqbd_solver is the only MAM algorithm that reads sn.lldscaling: it
        # applies the level-dependent factor in each departure block. 'default'
        # declares it because the single-class closed Delay+Queue shape routes
        # there, and supportsModelMethod refuses a load-dependent model outside
        # that shape -- a featset cannot see topology, and the dec.source
        # decomposition would otherwise solve every level at the nominal rate.
        if method in ('default', 'ldqbd'):
            feats.add('LoadDependence')
        # RETRIAL (recorded since 2026-09-05) is served by the BMAP/PH/N/N
        # retrial analyzer alone, which 'default' and 'dec.source' route to on
        # that shape and 'retrial' names. Every other algorithm reads no
        # sn.retrial* field and would answer with the refused jobs simply lost,
        # so the base grant is withdrawn there; supportsModelMethod refuses an
        # orbit OUTSIDE that shape for the two routing names.
        if method not in ('default', 'dec.source', 'retrial'):
            feats.discard('Retrial')
        # FINITECAPACITY (registry name since 2026-09-05): the open-network
        # analyzers carry a finite buffer as a loss buffer (M/M/c/K and
        # MMAP[K]/G/1/K in the basic analyzer, truncate-and-renormalize in the
        # decomposition and the open network analysis, the bufferless N/N
        # station of the retrial analyzer), so the base envelope declares it.
        # The two chains that read no sn.cap withdraw it: the LD-QBD levels run
        # to the population or the cutoff, and the background chain to the state
        # cap. A buffer a CLOSED class can fill is refused for every method by
        # mam_buffer_refusal (a closed job blocks, a loss formula does not).
        if method in ('ldqbd', 'bgchain'):
            feats.discard('FiniteCapacity')
        return feats

    get_method_feature_set = getMethodFeatureSet

    def unsupportedMethodReason(self, method):
        """The forwarding address for the RCAT names, which are SolverAG's now.

        Asks nothing of the model, so ``runAnalyzerChecks`` can call it before
        the struct is built; ``supportsModelMethod`` and ``runAnalyzer`` return
        the same string, so a caller gets one answer whichever gate it meets
        first.
        """
        if method in _RCAT_METHODS:
            return ('The %s method moved to SolverAG: RCAT decomposes the model into '
                    'cooperating agents rather than decomposing traffic, and no MAM '
                    "algorithm shares its machinery. Use SolverAG(model, '%s')."
                    % (method, method))
        return ''

    unsupported_method_reason = unsupportedMethodReason

    def supportsModelMethod(self, method):
        """Method-aware gate for MAM. Each decomposition algorithm declares a
        structural applicability predicate (supports_network), so delegate to it.
        For the special/auto-dispatched methods (mna, ldqbd, fj, retrial,
        reneging) that have no flat per-algorithm predicate, fall back to the
        per-method feature set. Returns (bool, reason)."""
        # G-network signals belong to the RCAT analyzer alone. Refuse them by
        # name here so the message says which method to call: the per-method
        # feature set already clears the signal names for every other method,
        # no MAM algorithm reads sn.issignal.
        moved = self.unsupportedMethodReason(method)
        if moved:
            return False, moved
        if _has_signal_class(self.sn):
            return False, ('The %s method does not support G-network signals: no MAM '
                           'algorithm reads sn.issignal, so this method would solve the '
                           'model with every signal turned into an ordinary customer. Use '
                           'SolverAG, whose RCAT methods are the only ones that model '
                           'signals.' % method)
        # A binding buffer a CLOSED class can fill: the base envelope declares
        # FiniteCapacity because the open-network analyzers carry a loss buffer,
        # and no feature name can say "open classes only" about it, so the
        # closed half is structural. Same predicate runAnalyzer raises.
        buffer_reason = mam_buffer_refusal(self.sn)
        if buffer_reason:
            return False, buffer_reason
        if method == 'retrial':
            # A "must be present" rule, which a feature set cannot state: it
            # says which constructs are ACCEPTED, so it can refuse a model for
            # having something and never for lacking it. solver_mam_retrial
            # needs an impatience configuration to analyze, and
            # mam_retrial_applicable is the same predicate
            # _solve_retrial_reneging asks before running it.
            ok, why = mam_retrial_applicable(self.sn)
            if not ok:
                return False, why
        if method == 'mna':
            # the deterministic split is carried by the open traffic equations
            # only; mna_closed has no counterpart
            from ...constants import RoutingStrategy as _RS
            routing = getattr(self.sn, 'routing', None)
            njobs = np.asarray(self.sn.njobs, dtype=float).ravel()
            if (routing is not None and np.any(np.asarray(routing) == int(_RS.RROBIN))
                    and not np.any(np.isinf(njobs))):
                return False, 'The mna method supports round-robin routing in open models only.'
        # ldqbd_solver is the only MAM algorithm that reads sn.lldscaling. A
        # decomposition method reaches its structural predicate below, which is
        # about topology and class mix and says nothing about the declared
        # features, so the load-dependent scaling would be dropped silently.
        if method not in ('default', 'ldqbd') and _has_nontrivial_lld(self.sn):
            return False, ('This model uses load-dependent service rates, which the %s method '
                           'does not read: it would solve every level at the nominal rate. Use '
                           'method \'ldqbd\', which applies the level-dependent factor in each '
                           'departure block.' % method)
        algo_class = SolverMAM.ALGORITHMS.get(method)
        if algo_class is not None:
            ok, reason = algo_class.supports_network(self.sn)
            if not ok:
                return False, (reason or '')
        # THE STRUCTURAL PREDICATE IS NOT THE WHOLE GATE, and returning on it
        # was the bug. supports_network answers about TOPOLOGY and class mix --
        # how many stations, open or closed, which scheduling shape the
        # decomposition needs -- and says nothing about the constructs the model
        # is built from, so a decomposition method sailed through it on a model
        # holding a Cache, a replacement strategy and a cache class switcher,
        # none of which SolverMAM declares. findSolver then offered five MAM
        # methods on a cache model and each died at run time with the very
        # feature message this gate should have carried. MATLAB ends the same
        # method with supportsModelMethod@NetworkSolver for exactly this reason;
        # falling through to the feature-set gate below is that ending.
        model = getattr(self, 'model', None)
        if model is None:
            return True, ''
        if hasattr(model, 'get_used_lang_features'):
            feat_used = model.get_used_lang_features()
        elif hasattr(model, 'getUsedLangFeatures'):
            feat_used = model.getUsedLangFeatures()
        else:
            return True, ''
        from ..base import SolverFeatureSet
        feat_supported = SolverFeatureSet()
        feat_supported.set_true(list(self.getMethodFeatureSet(method)))
        return SolverFeatureSet.supports_with_reason(feat_supported, feat_used)

    resolve_method = resolveMethod
    supports_model_method = supportsModelMethod

    @staticmethod
    def getFeatureSet() -> set:
        """Get set of features supported by SolverMAM.

        Returns the canonical feature names (mirrors MATLAB
        SolverMAM.getFeatureSet and the JAR SolverMAM).
        """
        return {
            'Sink', 'Source',
            'Fork', 'Join', 'Forker', 'Joiner',
            'Delay', 'DelayStation', 'Queue',
            'APH', 'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp', 'MMPP2', 'MAP', 'MMAP', 'DMAP', 'ME', 'RAP',
            # Geometric and DiscreteUniform are the lattice laws of the
            # discrete-time path, solved by the Q-MAM discrete-time queues
            'Geometric', 'DiscreteUniform',
            'Det', 'Gamma', 'Lognormal', 'Pareto', 'Uniform', 'Weibull',
            'StatelessClassSwitcher', 'InfiniteServer',
            'ClassSwitch',
            'SharedServer', 'Buffer', 'Dispatcher',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_HOL',
            'SchedStrategy_FCFS',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ClosedClass', 'SelfLoopingClass',
            'OpenClass',
            'Retrial', 'BMAP', 'PH',
            # Open stations are solved exactly by qbd_setupdelayoff; closed
            # stations use the per-instance cold-start race of the
            # hassetup branch.
            'SetupDelayOff',
            # c-server stations (every analyzer reads sn.nservers; the slotted
            # path's single-server rule stays structural) and finite buffers as
            # LOSS buffers, see getMethodFeatureSet and mam_buffer_refusal for
            # the closed-class refusal.
            'MultiServer', 'FiniteCapacity',
        }

    @staticmethod
    def defaultOptions() -> SolverMAMOptions:
        """Get default solver options.

        Returns:
            SolverMAMOptions with default values
        """
        return SolverMAMOptions()

    # =====================================================================
    # CDF AND PERCENTILE METHODS
    # =====================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """Get response time CDF as the exact matrix-analytic passage-time law.

        The native twin of MATLAB @SolverMAM/getCdfRespT.m: it runs
        solver_mam_passage_time on the Source + single queue open model
        (FCFS/HOL through the MMAPPH1FCFS sojourn PH law, PS through the
        Masuyama-Takine MAP/M/1-PS distribution). This used to fit an
        exponential to the mean, which agreed with the reference on the mean
        by construction and nowhere else.

        Args:
            R: Optional response time handles, accepted for signature
               compatibility and not read (the passage-time analysis computes
               every class at the queue anyway).

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys; empty, after
            a warning, on a topology the analysis does not cover.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import cdf_respt_via_jar
            return cdf_respt_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import mam_cdf_respt_via_cpp
            return mam_cdf_respt_via_cpp(self)
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        from .passage_time import solver_mam_passage_time
        return solver_mam_passage_time(self.sn, self.options)

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """Extract percentiles from response time distribution.

        Args:
            percentiles: List of percentiles (0-100). Default: [50, 75, 90, 95, 99]
            jobclass: Optional class filter (1-based)

        Returns:
            Tuple of (percentile_list, percentile_table)
        """
        if percentiles is None:
            percentiles = [50, 75, 90, 95, 99]

        percentiles = np.asarray(percentiles)
        percentiles = np.clip(percentiles, 0.01, 99.99)
        percentiles_normalized = percentiles / 100.0

        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        # Check for Fork-Join percentile results
        if hasattr(self.result, 'percentile_results') and self.result.percentile_results:
            fj_percs = self.result.percentile_results
            PercRT = [{
                'station': 'ForkJoin',
                'class': 1,
                'percentiles': fj_percs['percentiles'],
                'values': fj_percs['response_times'],
            }]
            rows = []
            for p, v in zip(fj_percs['percentiles'], fj_percs['response_times']):
                rows.append({'Percentile': p, 'RespT': v})
            return PercRT, pd.DataFrame(rows)

        R = self.result.RN
        nstations, nclasses = R.shape

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]

        # Extract station/class names if available
        station_names = getattr(self.sn, 'nodenames', None) or [f'Station{i}' for i in range(nstations)]
        class_names = getattr(self.sn, 'classnames', None) or [f'Class{r}' for r in range(nclasses)]

        # The real passage-time laws, where the analysis covers the model; the
        # exponential of the mean stays only as the fallback for cells with no law
        cdf_by_cell = {}
        try:
            for entry in self.getCdfRespT():
                cdf_by_cell[(entry['station'], entry['class'])] = entry
        except Exception:
            cdf_by_cell = {}

        for i in range(nstations):
            for r in range(nclasses):
                if jobclass is not None and (r + 1) != jobclass:
                    continue

                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                law = cdf_by_cell.get((i + 1, r + 1))
                if law is not None and len(law['t']) > 0:
                    # Invert the tabulated CDF by interpolation
                    perc_values = np.interp(percentiles_normalized,
                                            np.asarray(law['p']).ravel(),
                                            np.asarray(law['t']).ravel())
                else:
                    # Exponential approximation for percentiles
                    lambda_rate = 1.0 / mean_resp_t
                    perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                PercRT.append({
                    'station': i + 1,
                    'class': r + 1,
                    'percentiles': percentiles.tolist(),
                    'values': perc_values.tolist(),
                })

                row_data = {
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Class': class_names[r] if r < len(class_names) else f'Class{r}',
                }
                for perc_col, perc_val in zip(perc_col_names, perc_values):
                    row_data[perc_col] = perc_val
                rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    # =====================================================================
    # PROBABILITY METHODS
    # =====================================================================

    def getProb(self, node: int, state=None):
        """State probability of a (level, phase) pair, or the full matrix.

        Port of MATLAB @SolverMAM/getProb.m. QBD analysis is a single-queue
        method, so a network with more than one queue station is refused by
        name rather than approximated.

        Args:
            node: node index (0-based)
            state: [level, phase] pair, or None for the whole matrix

        Returns:
            scalar probability, or the (levels x phases) matrix when state is None
        """
        import numpy as np
        from ...api.mam.map_analysis import map_prob
        from ...lang.base import NodeType

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import mam_prob_via_cpp
            P = mam_prob_via_cpp(self, int(node) + 1)['joint']
            if state is None:
                return P
            level, phase = (int(state[0]), int(state[1])) if len(state) >= 2 \
                else (int(state[0]), 1)
            if level < 0 or phase < 1 or level + 1 > P.shape[0] or phase > P.shape[1]:
                return 0.0
            return float(P[level, phase - 1])

        sn = self.sn
        if node >= sn.nnodes:
            raise ValueError("getProb: node number exceeds the number of nodes in the model.")
        ist = int(sn.nodeToStation[node])
        if ist < 0:
            raise ValueError("getProb: specified node is not a station.")

        queue_stations = sum(1 for i in range(sn.nstations)
                             if sn.nodetype[int(sn.stationToNode[i])] == NodeType.QUEUE)
        if queue_stations > 1:
            raise ValueError(
                "getProb is not supported for networks with multiple queues in SolverMAM. "
                "The MAM solver uses QBD (quasi-birth-death) analysis, which is fundamentally "
                "a single-queue method. Use SolverCTMC or SolverSSA for state probabilities in "
                "networks with multiple queues.")
        if queue_stations == 0:
            raise ValueError("getProb: model does not contain any queue stations.")

        if self.result is None:
            self.runAnalyzer()

        Pstate = self._mam_joint_level_phase(ist)
        if state is None:
            return Pstate
        level, phase = int(state[0]), int(state[1])
        if level + 1 > Pstate.shape[0] or phase > Pstate.shape[1] or level < 0 or phase < 1:
            return 0.0
        return float(Pstate[level, phase - 1])

    def _mam_joint_level_phase(self, ist: int):
        """(levels x phases) joint, as MATLAB @SolverMAM/getProb.m builds it.

        Level marginal from MMAPPH1FCFS ncDistr over an aggregate MMAP built
        from the class throughputs; phase factor from the TIME-STATIONARY phase
        distribution map_prob, weighted by class throughput.

        map_prob, NOT map_pie: map_pie is the embedded equilibrium at departure
        instants, the phase a service STARTS in, which for an Erlang-2 is [1 0]
        and gave P(phase 2) = 0 at every level for a server spending half its
        busy time in phase 2. Phases are still taken independent of level, where
        an exact QBD has pi_n = pi_1 R^(n-1); that approximation is deliberate
        and matches the reference.
        """
        import numpy as np
        from ...api.mam.map_analysis import map_prob
        from ...lib.thirdparty.butools.queues import MMAPPH1FCFS

        sn = self.sn
        K = sn.nclasses
        N = np.asarray(sn.njobs, dtype=float).flatten()
        # THE OPEN LEVEL COUNT IS THE REFERENCE'S 100, not an arbitrary 20: an
        # open queue's length is unbounded and the table has to stop somewhere,
        # and @SolverMAM/getProbMarg.m stops at 100 unless options.cutoff names
        # another bound. A shorter table is not a coarser answer, it is a
        # TRUNCATED one, and it made this getter disagree with MATLAB and the
        # C++ on the length of the curve rather than on its values.
        if np.all(np.isfinite(N)):
            max_level = int(np.sum(N[np.isfinite(N)])) + 1
        else:
            cutoff = getattr(self.options, 'cutoff', None)
            max_level = 100
            if cutoff is not None:
                c = float(np.max(cutoff))
                if np.isfinite(c) and c > 0:
                    max_level = int(c)
        TN = np.atleast_2d(np.asarray(self.result.TN, dtype=float))

        # sn.proc is a compact descriptor dict in python ({'k','mu'}, {'rate'},
        # {'probs','rates'}), not MATLAB's {D0,D1} cell, so go through the
        # codebase's own converter rather than reading the entry positionally.
        from ...api.solvers.mam.handler import _extract_ph_for_phm1
        sigma, S = [], []
        for k in range(K):
            alpha, Tk = _extract_ph_for_phm1(sn, ist, k)
            if alpha is None:
                raise ValueError(
                    "getProb: class %d at station %d has no phase-type service "
                    "representation, so no QBD can be built." % (k, ist))
            # BUTOOLS IS WRITTEN AGAINST np.matrix, not ndarray: its queue
            # modules use `.I` for the inverse and `*` for the matrix product
            # throughout, so a plain array reaches `(-D0).I` and raises
            # AttributeError before any queue is analysed. Every other caller in
            # the tree wraps its arguments the same way (see
            # api/qsys/map_queues.py).
            sigma.append(np.matrix(np.asarray(alpha, dtype=float).reshape(1, -1)))
            S.append(np.matrix(np.asarray(Tk, dtype=float)))

        lambda_total = float(np.sum([TN[ist, k] for k in range(K)]))
        if not lambda_total > 0:
            raise ValueError("getProb: the station carries no throughput, so no QBD exists.")
        D_approx = [np.matrix([[-lambda_total]])]
        for k in range(K):
            D_approx.append(np.matrix([[TN[ist, k]]]))

        pdistr = np.abs(np.asarray(MMAPPH1FCFS(D_approx, sigma, S, 'ncDistr', max_level),
                                   dtype=float).flatten())
        ssum = pdistr.sum()
        if ssum > 0:
            pdistr = pdistr / ssum

        n_phases = max(s.shape[0] for s in S)
        avg_pie = np.zeros(n_phases)
        for k in range(K):
            exit_k = -S[k].sum(axis=1)
            D1_k = np.outer(exit_k, np.asarray(sigma[k]).flatten())
            piq = np.asarray(map_prob(S[k], D1_k), dtype=float).flatten()
            if piq.size == n_phases:
                avg_pie += piq * TN[ist, k] / lambda_total
        tot = avg_pie.sum()
        avg_pie = np.ones(n_phases) / n_phases if not tot > 0 or not np.isfinite(tot) \
            else avg_pie / tot

        levels = min(max_level, pdistr.size)
        self._last_level_distribution = pdistr[:levels]
        return np.outer(pdistr[:levels], avg_pie)

    def _mam_level_distribution(self, ist: int):
        """The QUEUE-LENGTH law alone, which is what getProbMarg returns.

        It is the level marginal of the same QBD `_mam_joint_level_phase`
        builds, so the two cannot disagree about how many jobs the station
        holds -- and it is the reference's own quantity: `getProbMarg.m` takes
        `pdistr` from `MMAPPH1FCFS(..., 'ncDistr', maxLevel)` and returns it as
        the curve, per class only through the class-marked arrival MMAP.
        """
        self._mam_joint_level_phase(ist)
        return np.asarray(self._last_level_distribution, dtype=float)

    def getProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Get marginal queue-length distribution at station for class.

        Args:
            station: Station index (0-based)
            jobclass: Job class index (0-based)

        Returns:
            Marginal probability vector P(n_ir) for n=0,1,2,...
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-marg', ist=station, jclass=jobclass, kind='vector', raw_station=True)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import mam_prob_via_cpp
            ind = int(np.ravel(self.sn.stationToNode)[int(station)]) + 1
            return mam_prob_via_cpp(self, ind)['marginal'][int(jobclass)]
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        # THE QBD LEVEL LAW, not a geometric fitted to the mean. The reference
        # (@SolverMAM/getProbMarg.m) takes `pdistr` from MMAPPH1FCFS's ncDistr,
        # and so does the C++ port; this used to return
        # (1-rho)*rho^n over max(10, 3*mean) points, which is an M/M/1 shape
        # whatever the service process is and stopped at 11 entries where both
        # references report 100. The two agreed on the first few values only
        # because the example was an M/M/1.
        return self._mam_level_distribution(int(station))

    # =====================================================================
    # ADDITIONAL STANDARD ACCESSOR METHODS
    # =====================================================================

    def sample(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample from state distribution (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sample() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleAggr(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleAggr() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleSys(self, numEvents: int = 1000) -> np.ndarray:
        """Sample system states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleSys() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleSysAggr(self, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated system states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleSysAggr() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    # =====================================================================
    # TRANSIENT METHODS (Not Supported - Steady-State Solver)
    # =====================================================================

    def getTranCdfRespT(self) -> List[Dict]:
        """Get transient response time CDF (not supported for MAM).

        Raises:
            NotImplementedError: MAM computes steady-state only
        """
        raise NotImplementedError("getTranCdfRespT() not supported for MAM. Use CTMC or simulation.")

    def getTranCdfPassT(self) -> List[Dict]:
        """Get transient passage time CDF (not supported for MAM).

        Raises:
            NotImplementedError: MAM computes steady-state only
        """
        raise NotImplementedError("getTranCdfPassT() not supported for MAM. Use FLD or simulation.")

    def getTranAvg(self, *args):
        """Get transient average metrics via QBD matrix exponentiation.

        Returns:
            Tuple of (QNt, UNt, TNt) where each is a nested list [M][K] of TranResult objects.
        """
        from ...constants import TranResult
        from .algorithms.ldqbd_transient import solver_mam_ldqbd_transient
        from .algorithms.transient_qbd import (
            solver_mam_transient_qbd, transient_qbd_applicable)

        # Set up timespan defaults if needed
        sn = self.sn if self.sn is not None else self._get_network_struct(self.network)
        rates = np.asarray(sn.rates)
        min_rate = np.nanmin(rates[rates > 0]) if np.any(rates > 0) else 1.0

        if not hasattr(self.options, 'timespan') or self.options.timespan is None:
            self.options.timespan = [0, 30.0 / min_rate]
        elif np.isinf(self.options.timespan[0]) and np.isinf(self.options.timespan[1]):
            self.options.timespan = [0, 30.0 / min_rate]
        elif np.isinf(self.options.timespan[0]):
            self.options.timespan[0] = 0
        elif np.isinf(self.options.timespan[1]):
            self.options.timespan[1] = 30.0 / min_rate

        # Auto-select: correlated MAP arrival/service or non-Poisson arrival on a
        # single-server open queue uses the Laplace-domain transient QBD solver
        # on the true MAP blocks; otherwise use the libQBD/expm fast path.
        if transient_qbd_applicable(sn):
            Qt, Ut, Tt = solver_mam_transient_qbd(sn, self.options)
        else:
            # Convert non-Markovian to PH if needed
            from ...api.sn import sn_nonmarkov_toph
            try:
                sn = sn_nonmarkov_toph(sn, self.options)
            except Exception:
                pass
            Qt, Ut, Tt = solver_mam_ldqbd_transient(sn, self.options)

        M = sn.nstations
        K = sn.nclasses

        QNt = [[None for _ in range(K)] for _ in range(M)]
        UNt = [[None for _ in range(K)] for _ in range(M)]
        TNt = [[None for _ in range(K)] for _ in range(M)]

        for i in range(M):
            for r in range(K):
                if Qt[i][r] is not None:
                    t_vals = Qt[i][r][:, 1]
                    QNt[i][r] = TranResult(t_vals, Qt[i][r][:, 0])
                    UNt[i][r] = TranResult(t_vals, Ut[i][r][:, 0])
                    TNt[i][r] = TranResult(t_vals, Tt[i][r][:, 0])

        return QNt, UNt, TNt

    def getCdfPassT(self) -> List[Dict]:
        """Get passage time CDF.

        For SolverMAM the passage time IS the response time: both come from the
        same solver_mam_passage_time call, so this delegates, as the JAR and the
        C++ CLI arms do.
        """
        return self.getCdfRespT()

    # =====================================================================
    # UNIFIED METRICS METHOD
    # =====================================================================


    # =====================================================================
    # CHAIN-LEVEL METHODS
    # =====================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self.sn, 'chains') and self.sn.chains is not None:
            chains = []
            for c in range(self.sn.nchains if hasattr(self.sn, 'nchains') else 1):
                chain_classes = []
                for k in range(self.sn.nclasses):
                    if hasattr(self.sn.chains, '__getitem__'):
                        if self.sn.chains[c, k] > 0:
                            chain_classes.append(k)
                chains.append(chain_classes)
            return chains if chains else [[k for k in range(self.sn.nclasses)]]
        else:
            # Default: each class is its own chain
            return [[k] for k in range(self.sn.nclasses)]

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()

        nstations, nchains = QN.shape
        rows = []

        station_names = getattr(self.sn, 'nodenames', None) or [f'Station{i}' for i in range(nstations)]
        chain_names = [f'Chain{c}' for c in range(nchains)]

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Chain': chain_names[c],
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

    # =====================================================================
    # NODE-LEVEL METHODS
    # =====================================================================

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get average metrics per node.

        Unlike getAvg() which returns station-level metrics, this method
        returns node-level metrics including non-station nodes (e.g., Router/VSink).

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        from ...api.sn.getters import sn_get_node_arvr_from_tput, sn_get_node_tput_from_tput

        if self.result is None:
            self._ensureAvgResults()

        TN = self.result.TN
        QN = self.result.QN
        UN = self.result.UN
        RN = self.result.RN

        sn = self.sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses

        # Create TH (throughput handle) - indicates which station-classes have valid throughput
        TH = np.zeros_like(TN)
        TH[TN > 0] = 1.0

        # Compute node arrival rates and throughputs using helper functions
        # Pass AN=None to let sn_get_node_arvr_from_tput compute it properly
        # (including setting Source arrival rates to 0)
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH, None)
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
                WNn[ind, :] = RN[ist, :]

        return QNn, UNn, RNn, WNn, ANn, TNn

    def getAvgNodeTable(self) -> pd.DataFrame:
        """
        Get average metrics by node as DataFrame.

        Returns node-based results (one row per node per class) including
        non-station nodes like Router/VSink.

        Returns:
            pandas.DataFrame with columns: Node, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()

        sn = self.sn
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        class_names = list(sn.classnames) if hasattr(sn, 'classnames') and sn.classnames else []

        rows = []
        for node_idx in range(sn.nnodes):
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

            for r in range(sn.nclasses):
                class_name = class_names[r] if r < len(class_names) else f'Class{r}'

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
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T) where R is system response time and T is system throughput
        """
        R = self.getAvgSysRespT()
        T = self.getTput()
        return R, T

    getAvgSysTable = NetworkSolver.getAvgSysTable  # chain-level shared layout

    def getMAMResult(self):
        """Intermediate quantities of the matrix-analytic analysis of a
        single-queue model, in addition to the mean performance measures
        returned by getAvg.

        Mean values alone hide the objects the method is actually built on, so
        a matrix-analytic result cannot be inspected, taught, or checked
        against a published derivation. This accessor returns them.

        For a BMAP (or MAP) arrival stream feeding an exponential single server
        the result is that of qsys_bmapm1 and carries the M/G/1-type
        quantities: the phase-process stationary vectors theta and alpha, the
        randomized blocks A0, A1, B0 and Bk, the matrix G, the drift, the
        measured decay rate and the level probabilities.

        For a retrial station the result is that of qsys_bmapphnn_retrial and
        carries the orbit-level stationary distribution together with the
        truncation level and its residual.
        """
        from ...api.qsys.retrial import qsys_is_retrial, solver_mam_retrial
        from ...api.qsys.bmapm1 import qsys_bmapm1
        from ...api.sn.network_struct import NodeType

        sn = self.sn

        # A retrial station carries its own engine, whose result object already
        # exposes the orbit-level internals.
        is_retrial, _ = qsys_is_retrial(sn)
        if is_retrial:
            options = dict(self.options) if isinstance(self.options, dict) else {}
            options['verbose'] = False
            return solver_mam_retrial(sn, options)[-1]

        # Otherwise: BMAP/MAP arrivals into a single exponential server.
        source_idx = None
        queue_idx = None
        for ist in range(sn.nstations):
            node_idx = int(sn.stationToNode[ist])
            if sn.nodetype[node_idx] == NodeType.SOURCE:
                source_idx = ist
            elif sn.nodetype[node_idx] == NodeType.QUEUE:
                if queue_idx is None:
                    queue_idx = ist
                else:
                    raise ValueError('getMAMResult exposes the matrix-analytic '
                                     'internals of a single-queue model only.')
        if source_idx is None or queue_idx is None:
            raise ValueError('getMAMResult requires an open model with one '
                             'Source and one Queue.')
        if sn.nclasses > 1:
            raise ValueError('getMAMResult exposes the matrix-analytic internals '
                             'of a single-class model only.')
        if int(sn.nservers[queue_idx]) != 1:
            raise ValueError('getMAMResult requires a single-server queue.')

        from ...api.qsys.retrial import _proc_to_d0d1
        arrival_proc = sn.proc[source_idx][0]
        if isinstance(arrival_proc, dict):
            arrival_proc = _proc_to_d0d1(arrival_proc)
        if arrival_proc is None or not isinstance(arrival_proc, (list, tuple)) \
                or len(arrival_proc) < 2:
            raise ValueError('The arrival process has no Markovian (D0,D1,...) '
                             'representation.')

        service_proc = sn.proc[queue_idx][0]
        if isinstance(service_proc, dict):
            service_proc = _proc_to_d0d1(service_proc)
        if service_proc is None or np.atleast_2d(service_proc[0]).shape[0] != 1:
            raise ValueError('getMAMResult exposes the M/G/1-type internals for '
                             'exponential service only; the queue has a '
                             'multi-phase service process.')
        mu = -float(np.atleast_2d(service_proc[0])[0, 0])

        return qsys_bmapm1([np.atleast_2d(Dk) for Dk in arrival_proc], mu)

    get_mam_result = getMAMResult

    # =====================================================================
    # ALIASES (PascalCase for MATLAB compatibility)
    # =====================================================================

    GetAvg = NetworkSolver.getAvg
    GetAvgChainTable = getAvgChainTable
    GetAvgNode = getAvgNode
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetAvgSys = getAvgSys
    GetAvgSysTable = getAvgSysTable
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    GetProb = getProb
    GetProbMarg = getProbMarg
    GetTranAvg = getTranAvg

    # Node-chain specific aliases
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain

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
    run_analyzer = runAnalyzer
    cdf_resp_t = getCdfRespT
    perct_resp_t = getPerctRespT
    list_valid_methods = listValidMethods
    default_options = defaultOptions


__all__ = [
    'SolverMAM',
    'SolverMAMOptions',
    'MAMResult',
]
