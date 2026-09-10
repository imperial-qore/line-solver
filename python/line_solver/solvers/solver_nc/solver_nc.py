"""
Native Python implementation of NC (Normalizing Constant) solver.

This implementation uses pure Python/NumPy algorithms from the api.solvers.nc
module.
"""

import os
import re
import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass, field
from ...constants import default_verbose

from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn.getters import sn_get_arvr_from_tput
from ...api.io.logging import line_debug, line_warning, LineError
from ..base import NetworkSolver, method_label, method_type
from ..fork_join_driver import ForkJoinDriverMixin
from ..transform_driver import TransformSolveMixin
from .nc_multiserver import nc_multiserver_policy, nc_lld_from_nservers


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


# The refusal text for the two cache method names moved to nc_method_refusal, which is
# where both the gate and the analyzer now read it.


@dataclass
class SolverNCOptions:
    """Options for the native NC solver."""
    method: str = 'default'
    tol: float = 1e-4
    iter_max: int = 1000
    iter_tol: float = 1e-4
    verbose: bool = field(default_factory=default_verbose)
    seed: int = 1  # Random seed (for compatibility, NC is deterministic)
    keep: bool = False  # Keep intermediate data (for compatibility)
    cutoff: Optional[int] = None  # State space cutoff (for compatibility)
    samples: Optional[int] = None  # Samples (for compatibility)
    timeout: float = float('inf')  # Wall-clock time budget in seconds (inf = no budget)
    lang: str = field(default_factory=lambda: os.environ.get('LINE_SOLVER_LANG', 'python'))  # env LINE_SOLVER_LANG overrides; 'python' (native), 'java' (jline.jar via JSON) or 'cpp' (line-cli via JSON)
    # Arithmetic backend, lang='cpp' ONLY: 'double' (default), 'exact' or
    # 'real:<digits>'. Meaningless for the other langs, which are IEEE double
    # throughout, so line-cli is invoked without --arith unless the caller sets it.
    arith: Optional[str] = None
    # Solver-specific switches. config['slotted'] selects the discrete-time
    # (slotted) route of solver_nc_dt_analyzer, mirroring MATLAB
    # options.config.slotted; config['slotlength'] gives the slot length in
    # model time units.
    config: Optional[dict] = None


def _is_lossn_fcr_case(sn) -> bool:
    """True if sn is the NC-supported loss network: an open model with a
    single finite capacity region whose only constrained station is an
    infinite-server (Delay) node using the DROP policy for all classes.
    This case is dispatched to solver_nc_lossn_analyzer."""
    if sn is None:
        return False
    try:
        from ...api.sn.predicates import sn_has_closed_classes
        from ...lang.base import DropStrategy
        if sn_has_closed_classes(sn):
            return False
        if not (getattr(sn, 'nregions', 0) == 1 and
                getattr(sn, 'region', None) is not None and len(sn.region) > 0):
            return False
        region_matrix = sn.region[0]
        K = sn.nclasses
        M = sn.nstations
        nservers = sn.nservers
        stations_in_fcr = []
        for i in range(M):
            has_class = np.any(region_matrix[i, :K] >= 0)
            has_global = region_matrix[i, K] >= 0 if region_matrix.shape[1] > K else False
            if has_class or has_global:
                stations_in_fcr.append(i)
        if len(stations_in_fcr) != 1 or not np.isinf(nservers[stations_in_fcr[0]]):
            return False
        if getattr(sn, 'regionrule', None) is None:
            return False
        for r in range(K):
            if sn.regionrule[0, r] != float(DropStrategy.DROP):
                return False
        return True
    except Exception:
        return False


def _partitions(n: int, Nmax: List[int]) -> List[List[int]]:
    """Every way of splitting n jobs across the classes, class r capped at Nmax[r].

    The compositions `generate_partitions` in `@SolverNC/getProbMarg.m`
    enumerates: ordered, so (1,0) and (0,1) are two of them, since the classes
    are distinguishable.
    """
    R = len(Nmax)
    if R == 0:
        return [[]] if n == 0 else []
    if R == 1:
        return [[n]] if n <= Nmax[0] else []
    out = []
    for first in range(0, min(n, Nmax[0]) + 1):
        for rest in _partitions(n - first, Nmax[1:]):
            out.append([first] + rest)
    return out


class SolverNC(TransformSolveMixin, ForkJoinDriverMixin, NetworkSolver):
    """
    Native Python NC (Normalizing Constant) solver.

    This solver analyzes product-form queueing networks using normalizing
    constant computation methods in pure Python/NumPy, providing the same
    functionality as the Java wrapper without requiring the JVM.

    Supported methods:
        - 'default': Automatic method selection
        - 'exact': Exact convolution
        - 'comom': Approximate method

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
            # Check if method was passed as keyword argument
            self.method = kwargs.pop('method', 'default')
            if isinstance(self.method, str):
                self.method = self.method.lower()
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
            # Remove 'method' from kwargs if present to avoid duplicate argument
            kwargs.pop('method', None)
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'iter_max'):
                kwargs.setdefault('iter_max', method_or_options.iter_max)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            # config carries per-solver payloads (the multiserver rule, the
            # discrete-time slot settings); dropping it makes them unreachable
            # through an options object. Mirrors solver_mva.py
            if method_or_options.get('config') is not None:
                kwargs.setdefault('config', method_or_options.get('config'))
            kwargs.pop('method', None)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'iter_max'):
                kwargs.setdefault('iter_max', method_or_options.iter_max)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if getattr(method_or_options, 'config', None) is not None:
                kwargs.setdefault('config', method_or_options.config)
            kwargs.pop('method', None)
        else:
            self.method = kwargs.pop('method', 'default')
            if isinstance(self.method, str):
                self.method = self.method.lower()

        self.options = SolverNCOptions(method=self.method, **kwargs)

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "NC"

    get_name = getName

    def supportsExactSensitivity(self):
        """The normalizing-constant solver is exact on the same product-form
        class that pfqn_sens differentiates, so getSensitivityTable uses the
        analytic branch.
        """
        return True

    supports_exact_sensitivity = supportsExactSensitivity

    def reset(self):
        """Reset solver state to force recomputation on next getAvg() call.

        This is called by ensemble solvers (like LN) after updating layer
        parameters to ensure the solver recomputes with new values.
        """
        self._clearResultStores()
        # Re-extract network parameters to pick up changes
        self._extract_network_params()

    def _extract_network_params(self):
        """Extract parameters from the model for NC computation."""
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

    def _transform_publish(self, tr, method):
        """SolverNC keeps a SolverNCReturn, not the dict the mixin defaults to.

        Same obstacle as _fj_publish above: python has no single result
        container, so each solver maps the driver's neutral answer onto its own.
        """
        from ...api.solvers.nc.handler import SolverNCReturn
        import numpy as _np
        nchains = int(getattr(self._sn, 'nchains', _np.asarray(tr.Q).shape[1]))
        self._result = SolverNCReturn(
            Q=tr.Q, U=tr.U, R=tr.R, T=tr.T,
            nchains=nchains, X=tr.X,
            lG=float(tr.lG) if tr.lG is not None and _np.isfinite(tr.lG) else 0.0,
            STeff=None, it=int(tr.iter), runtime=float(tr.runtime), method=method)
        return self._result

    def _fj_publish(self, result):
        """Convert a fork-join result dict into the NC result container.

        The shared driver (ForkJoinDriverMixin) speaks the dict contract that
        SolverMVA uses natively; SolverNC's getters read a SolverNCReturn, so
        the dict is mapped onto its single-letter fields here. Mirrors the JAR
        SolverNC.ncDispatch, which converts an NCResult into the neutral
        carrier and back.
        """
        from ...api.solvers.nc.handler import SolverNCReturn
        import numpy as _np
        QN = result.get('QN')
        nchains = int(getattr(self._sn, 'nchains', _np.asarray(QN).shape[1] if QN is not None else 0))
        return SolverNCReturn(
            Q=QN, U=result.get('UN'), R=result.get('RN'), T=result.get('TN'),
            nchains=nchains, X=result.get('XN'),
            lG=float(result.get('lG', 0.0)), STeff=None,
            it=int(result.get('iter', 0)), runtime=float(result.get('runtime', 0.0)),
            method=str(result.get('method', 'mmt')))

    def _fj_inner_solver(self, nonfjmodel, method=None):
        """Inner solve of the fork-join fixed point, on the NC analyzer.

        Overrides ForkJoinDriverMixin._fj_inner_solver, whose default is
        SolverMVA. No method override is applied: the transformed model carries
        auxiliary open classes at a vanishing rate, which the default
        normalizing-constant route already resolves (forcing 'rd' is silently
        ignored by pfqn_nc, as the MATLAB port found).
        """
        return SolverNC(nonfjmodel)

    def runAnalyzer(self) -> 'SolverNC':
        """Run the NC analysis."""
        # MODEL TRANSFORMATION, opt-in through options.config['transform'].
        # Mirrors the branch MATLAB puts in the shared runAnalyzerPreamble, so a
        # strategy written once serves NC as well as MVA and CTMC.
        if self.maybe_transform():
            return self
        # A fresh analysis invalidates any prior unstable-utilization cap.
        self._unstable_util_capped = False
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
                line_warning("SolverNC",
                             "lang='cpp' requested but the C++ solver is unavailable (%s); "
                             "falling back to lang='python'." % e)

        import numpy as np
        import warnings

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn, RNG
        # seeding, and sample/seed forwarding")
        _seed = getattr(self.options, 'seed', None)
        if _seed is not None:
            np.random.seed(int(_seed))

        from ...api.solvers.nc.handler import (
            solver_nc, solver_ncld, SolverOptions as HandlerOptions
        )
        from ...api.solvers.nc.analyzers import solver_nc_lossn_analyzer
        from ...api.sn.predicates import sn_has_closed_classes
        from ...lang.base import DropStrategy

        line_debug("NC: using lang=python", options=self.options)

        # see _kb/06-solver-catalog.md ("Finite capacity gate (MVA and NC)")
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_used_lang_features'):
            self.runAnalyzerChecks(self.options)

        # Reject MAP/MMPP2 explicitly: the featset gate cannot see process
        # types, so NC would otherwise silently treat correlated arrivals as Poisson.
        if self._sn is not None and getattr(self._sn, 'procid', None) is not None:
            from ...constants import ProcessType as _PT
            _nc_reject = {_PT.MAP: 'MAP', _PT.MMPP2: 'MMPP2'}
            for _v in np.asarray(self._sn.procid, dtype=object).ravel():
                for _pt, _nm in _nc_reject.items():
                    if _v == _pt:
                        raise RuntimeError(
                            "SolverNC does not support the %s process used by this "
                            "model (not in the NC feature set; a non-renewal MAP "
                            "cannot be captured by a product-form solver). Use "
                            "SolverMAM, SolverCTMC, or SolverSSA." % _nm)

        # see _kb/06-solver-catalog.md (NC: "Unknown NC methods are rejected, not silently defaulted")
        origmethod = self.options.method
        if origmethod not in self.listValidMethods():
            line_debug("NC: unrecognized method '%s', falling back to the default normalizing-constant analyzer (nc_analyzer/comom).", origmethod, options=self.options)

        sn = self._sn

        # THE STRUCTURAL METHOD GATE, asked once and in one place.
        #
        # nc_method_refusal holds every rule of the form "this method has no
        # route on this model": the DPS shape, state-dependent routing, the
        # Petri net, the order-independent rank rate, the cache and loss-network
        # tokens, PANACEA's normal usage. supportsModelMethod asks the SAME
        # function, which is what keeps Network.findSolver from offering a pair
        # that would raise here. It runs unconditionally, ahead of the feature
        # gate above, so it also holds when checks are disabled.
        from .nc_method_refusal import nc_method_refusal
        # for_report=False: this is the RUN, and it asks what the reference DOES
        # rather than what the report should offer. The two answers differ for
        # 'mmint2'/'gleint', which pfqn_nc answers with an empty constant and a
        # zero table; see nc_method_refusal.
        _refusal = nc_method_refusal(sn, getattr(self.options, 'method', 'default'),
                                     self.options, for_report=False)
        if _refusal:
            raise ValueError(_refusal)

        # Closed think+DPS network -> solver_nc_dps_analyzer (Morrison's heavy-usage
        # generating-function expansion), the DEFAULT for that shape. Intercepted
        # FIRST, ahead of every other route: 'SchedStrategy_DPS' is declared in the
        # feature set, which opens all of them to a DPS model, and each would silently
        # drop the weights and answer with the egalitarian-PS network. Not a
        # product-form route: lG is NaN. See _kb/06-solver-catalog.md (NC section).
        from .solver_nc_dps_analyzer import nc_is_dps_model, solver_nc_dps_analyzer
        if sn is not None and nc_is_dps_model(sn):
            if str(getattr(self.options, 'method', 'default')).lower() in ('default', 'morrison'):
                from ...api.solvers.nc.handler import SolverNCReturn
                line_debug("NC: closed think+DPS network, routing to solver_nc_dps_analyzer (Morrison)",
                           options=self.options)
                QN, UN, RN, TN, CN, XN, lG, dps_rt, dps_it, dps_method = solver_nc_dps_analyzer(sn, self.options)
                self._result = SolverNCReturn(
                    Q=QN, U=UN, R=RN, T=TN,
                    nchains=int(getattr(sn, 'nchains', 1)),
                    X=XN, lG=float(lG),
                    STeff=np.zeros_like(QN), it=int(dps_it),
                    runtime=dps_rt, method=dps_method,
                )
                self._extract_names()
                if self.options.verbose:
                    import sys
                    py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                    from line_solver.solvers.base import print_solver_banner
                    print_solver_banner(f"NC analysis [method: {dps_method}; type: {method_type('NC', dps_method)}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
                return self
            # The three refusal arms that used to follow -- another method on a DPS
            # model, a DPS station outside Morrison's shape, and 'morrison' on a model
            # with no DPS station at all -- moved into nc_method_refusal above, with
            # their wording unchanged.

        # Krzesinski state-dependent routing: the model has its own product
        # form (eq. 16), so it is intercepted before the standard convolution
        # and MVA analyzers, which assume state-independent routing
        if getattr(sn, 'sdr', None):
            from ...api.solvers.nc.handler import SolverNCReturn
            from .solver_nc_sdr_analyzer import solver_nc_sdr_analyzer
            line_debug("NC: state-dependent routing, routing to solver_nc_sdr_analyzer", options=self.options)
            QN, UN, RN, TN, CN, XN, lG, sdr_rt, sdr_it, sdr_method = solver_nc_sdr_analyzer(sn, self.options)
            self._result = SolverNCReturn(
                Q=QN, U=UN, R=RN, T=TN,
                nchains=int(getattr(sn, 'nchains', 1)),
                X=XN, lG=float(lG),
                STeff=np.zeros_like(QN), it=int(sdr_it),
                runtime=sdr_rt, method=sdr_method,
            )
            self._extract_names()
            if self.options.verbose:
                import sys
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"NC analysis [method: {sdr_method}; type: {method_type('NC', sdr_method)}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
            return self

        # Discrete-time (slotted) route: explicit request only, and an error
        # rather than a fallback when the model is outside the discrete-time
        # product form; see _kb/06-solver-catalog.md (NC section)
        from .solver_nc_dt_analyzer import is_slotted, solver_nc_dt_analyzer
        if is_slotted(self.options):
            from ...api.solvers.nc.handler import SolverNCReturn
            line_debug("NC: slotted model, routing to solver_nc_dt_analyzer", options=self.options)
            QN, UN, RN, TN, CN, XN, lG, dt_rt, dt_it, dt_method = solver_nc_dt_analyzer(sn, self.options)
            self._result = SolverNCReturn(
                Q=QN, U=UN, R=RN, T=TN,
                nchains=int(getattr(sn, 'nchains', 1)),
                X=XN, lG=float(lG),
                STeff=np.zeros_like(QN), it=int(dt_it),
                runtime=dt_rt, method=dt_method,
            )
            self._extract_names()
            if self.options.verbose:
                import sys
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"NC analysis [method: {dt_method}; type: {method_type('NC', dt_method)}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
            return self

        # A stochastic Petri net takes the MDD-rec route: the reachable set lives
        # in a decision diagram and the product form supplies the rates, so none
        # of the queueing-network branches below apply to it.
        from ...api.sn.network_struct import NodeType as _NodeType
        if sn is not None and np.any(np.ravel(np.asarray(sn.nodetype, dtype=int))
                                     == int(_NodeType.PLACE)):
            line_debug("Stochastic Petri net, routing to nc_spn_analyzer (MDD-rec)",
                       options=self.options)
            return self._run_spn_analyzer(sn)

        # see _kb/06-solver-catalog.md (NC: "Fork-join (all three codebases)")
        if self._has_fork_join() and not getattr(self, '_skip_fork_join', False):
            line_debug("NC: fork-join network detected, routing to fork_join_analysis", options=self.options)
            fj_result = self._run_fork_join_analysis()
            if fj_result is not None:
                self._extract_names()
                return fj_result

        if sn is not None and getattr(sn, 'immfeed', None) is not None and np.any(sn.immfeed):
            line_warning("SolverNC", "SolverNC does not handle immediate feedback (immfeed); the solver will treat self-loops as class-switching with re-queueing.")

        # Check if model contains Cache nodes - use specialized cache analyzer
        has_cache = False
        if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
            from ...lang.nodes import Cache
            for node in self.model._nodes:
                if isinstance(node, Cache):
                    has_cache = True
                    break

        if has_cache:
            line_debug("Non-reentrant cache (Source-Cache-Sink), routing to nc_cache_analyzer", options=self.options)
            return self._runCacheAnalyzer()

        # see _kb/06-solver-catalog.md (NC: "Analyzer routing order: OI exact,
        # PAS/OI importance sampling, MEM")
        from .solver_nc_oi_analyzer import nc_is_oi_model, solver_nc_oi_analyzer
        # The arm that used to stand here -- every method other than
        # default/exact/is/sampling on an OI model, which would silently drop the
        # rank rate mu(n) -- moved into nc_method_refusal, wording unchanged.
        if sn is not None and nc_is_oi_model(sn) and self.options.method in ('default', 'exact'):
            from ...api.solvers.nc.handler import SolverNCReturn
            line_debug("NC analyzer routing to solver_nc_oi_analyzer (order-independent)", options=self.options)
            QN, UN, RN, TN, CN, XN, lG, oi_rt, oi_it, oi_method = solver_nc_oi_analyzer(sn, self.options)
            self._result = SolverNCReturn(
                Q=QN, U=UN, R=RN, T=TN,
                nchains=int(getattr(sn, 'nchains', 1)),
                X=XN, lG=float(lG),
                STeff=np.zeros_like(QN), it=int(oi_it),
                runtime=oi_rt, method=oi_method,
            )
            self._extract_names()
            if self.options.verbose:
                import sys
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"NC analysis [method: {method_label(self.options.method, oi_method)}; type: {method_type('NC', method_label(self.options.method, oi_method))}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
            return self

        # see _kb/06-solver-catalog.md (NC: "Analyzer routing order")
        from .solver_nc_pas_is_analyzer import nc_is_pas_model, solver_nc_pas_is_analyzer
        _is_pas = sn is not None and nc_is_pas_model(sn)
        # An OI model (Delay + OI cycle included) routes here on 'is'/'sampling'
        # only; a P&S tandem takes 'default' as well, having no exact path.
        _is_oi = sn is not None and nc_is_oi_model(sn)
        if (_is_pas and self.options.method in ('default', 'is', 'sampling')) \
                or (_is_oi and self.options.method in ('is', 'sampling')):
            from ...api.solvers.nc.handler import SolverNCReturn
            line_debug("NC analyzer routing to solver_nc_pas_is_analyzer (pass-and-swap IS)", options=self.options)
            QN, UN, RN, TN, CN, XN, lG, ps_rt, ps_it, ps_method = solver_nc_pas_is_analyzer(sn, self.options)
            self._result = SolverNCReturn(
                Q=QN, U=UN, R=RN, T=TN,
                nchains=int(getattr(sn, 'nchains', 1)),
                X=XN, lG=float(lG),
                STeff=np.zeros_like(QN), it=int(ps_it),
                runtime=ps_rt, method=ps_method,
            )
            self._extract_names()
            if self.options.verbose:
                import sys
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"NC analysis [method: {method_label(self.options.method, ps_method)}; type: {method_type('NC', method_label(self.options.method, ps_method))}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
            return self

        # see _kb/06-solver-catalog.md (NC: "Analyzer routing order")
        if self.options.method == 'is' and sn is not None and sn.njobs is not None \
                and bool(np.any(np.isinf(np.asarray(sn.njobs, dtype=float)))):
            raise ValueError(
                "The 'is' importance-sampling method requires a closed queueing "
                "network. Use 'sampling' (pfqn_mci/pfqn_ls) for open or mixed models.")

        # MEM (Kouvatsos 1994): explicit request only; 'default' keeps the native analyzer.
        use_mem = self.options.method == 'mem'
        if use_mem:
            import time as _time
            from ...api.me import solver_nc_mem
            from ...api.solvers.nc.handler import SolverNCReturn
            line_debug("NC method=mem, routing to solver_nc_mem (Maximum Entropy, open QN)", options=self.options)
            _t0 = _time.time()
            solver_nc_mem.last_method = 'mem'
            QN, UN, RN, TN, CN, XN, mem_iter = solver_nc_mem(sn, self.options)
            mem_method = getattr(solver_nc_mem, 'last_method', 'mem')
            self._result = SolverNCReturn(
                Q=QN, U=UN, R=RN, T=TN,
                nchains=int(getattr(sn, 'nchains', 1)),
                X=XN, lG=float('nan'),
                STeff=np.zeros_like(QN), it=int(mem_iter),
                runtime=_time.time() - _t0, method=mem_method,
            )
            self._extract_names()
            if self.options.verbose:
                import sys
                py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
                from line_solver.solvers.base import print_solver_banner
                print_solver_banner(f"NC analysis [method: {method_label(self.options.method, mem_method)}; type: {method_type('NC', method_label(self.options.method, mem_method))}; lang: python; env: {py_version}] completed in {self._result.runtime:.6f}s.")
            return self

        # Open loss network with FCR (MATLAB runAnalyzer.m:138-161): no closed
        # classes, exactly 1 FCR, single Delay station in it, all DROP classes.
        nservers = sn.nservers.flatten().copy() if sn.nservers is not None else np.ones(sn.nstations)
        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn")
        _orig_nservers = sn.nservers.copy() if sn.nservers is not None else None
        _orig_lldscaling = sn.lldscaling.copy() if (hasattr(sn, 'lldscaling') and sn.lldscaling is not None) else None

        if (not sn_has_closed_classes(sn) and
                hasattr(sn, 'nregions') and sn.nregions == 1 and
                hasattr(sn, 'region') and sn.region is not None and len(sn.region) > 0):

            region_matrix = sn.region[0]
            K = sn.nclasses
            M = sn.nstations

            # Find stations in FCR (those with non-negative constraints)
            stations_in_fcr = []
            for i in range(M):
                has_class_constraint = np.any(region_matrix[i, :K] >= 0)
                has_global_constraint = region_matrix[i, K] >= 0 if region_matrix.shape[1] > K else False
                if has_class_constraint or has_global_constraint:
                    stations_in_fcr.append(i)

            # Check if single Delay node in FCR with DROP policy
            if (len(stations_in_fcr) == 1 and
                    np.isinf(nservers[stations_in_fcr[0]]) and
                    hasattr(sn, 'regionrule') and sn.regionrule is not None):

                # Check if all classes have DROP policy
                all_drop = True
                for r in range(K):
                    if sn.regionrule[0, r] != float(DropStrategy.DROP):
                        all_drop = False
                        break

                if all_drop:
                    line_debug("Open model with single FCR + Delay node (DROP), routing to nc_lossn_analyzer", options=self.options)
                    # Use loss network solver
                    _samples = getattr(self.options, 'samples', None)
                    handler_options = HandlerOptions(
                        method=self.options.method,
                        tol=self.options.tol,
                        iter_max=self.options.iter_max,
                        iter_tol=self.options.iter_tol,
                        verbose=1 if self.options.verbose else 0,
                        samples=int(_samples) if _samples else 100000,
                        seed=getattr(self.options, 'seed', None)
                    )
                    result = solver_nc_lossn_analyzer(sn, handler_options)
                    # Convert NCResult to handler result format
                    from dataclasses import dataclass

                    @dataclass
                    class LossnResult:
                        Q: np.ndarray
                        U: np.ndarray
                        R: np.ndarray
                        T: np.ndarray
                        X: np.ndarray
                        lG: float
                        STeff: np.ndarray
                        it: int
                        runtime: float
                        method: str

                    self._result = LossnResult(
                        Q=result.QN,
                        U=result.UN,
                        R=result.RN,
                        T=result.TN,
                        X=result.XN.flatten() if result.XN is not None else np.zeros(K),
                        lG=result.lG if result.lG is not None else np.nan,
                        STeff=np.zeros((M, K)),
                        it=result.it,
                        runtime=result.runtime,
                        method=result.method
                    )
                    self._extract_names()
                    return self

        # The token gates that used to stand here -- 'ms'/'erlangfp' and 'rec' off a
        # loss network, and 'rayint'/'spm' with no Cache node -- moved into
        # nc_method_refusal, which decides them from the same struct before the
        # dispatch begins and which the support gate asks too; their wording is
        # unchanged. The six load-dependent evaluators on an OPEN chain are now
        # refused by the feature set instead (getMethodFeatureSet drops OpenClass
        # from them): "closed population only" is a rule the registry CAN name, and
        # naming it there is what makes Network.findSolver drop the row rather than
        # report it runnable.

        # see _kb/06-solver-catalog.md (NC: "Multiserver -> load-dependent conversion")
        use_ld_solver = False
        # sn and nservers already computed above

        # Check for multi-server stations
        has_multiserver = any(s > 1 and np.isfinite(s) for s in nservers)

        # Check for open classes
        has_open = False
        if sn.njobs is not None:
            njobs = sn.njobs.flatten()
            has_open = any(np.isinf(njobs))

        method = self.options.method

        # How this model's finite multiserver stations are represented: Seidmann's
        # approximation or the exact mu(n)=min(n,c) lattice. The shipped 'default'
        # reproduces the historical dispatch exactly, so no result moves unless
        # config.multiserver is set. See _kb/06-solver-catalog.md (NC section)
        ms_policy = nc_multiserver_policy(
            self.options, warn=lambda m: line_warning("SolverNC", m))

        # see _kb/06-solver-catalog.md (NC: "Multiserver -> load-dependent conversion")
        if (method in ('exact', 'is', 'panald') and not _is_pas and has_multiserver
                and not has_open and ms_policy == 'seidmann'):
            # config.multiserver='seidmann' asks for Seidmann's approximation on
            # this arm too, so the conversion below is skipped. Off by default
            line_debug("%s method: config.multiserver=seidmann, keeping Seidmann "
                       "approximation for multiserver stations" % method, options=self.options)
        elif method in ('exact', 'is', 'panald') and not _is_pas and has_multiserver and not has_open:
            line_debug("%s method: converting multiserver stations to load-dependent" % method, options=self.options)
            # Transform multi-server nodes into lldscaling (like MATLAB lines 69-80 in runAnalyzer.m)
            njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(sn.nclasses)
            Nt = int(np.sum(njobs[np.isfinite(njobs)]))
            if Nt > 0:
                # Create lldscaling matrix
                lldscaling = np.ones((sn.nstations, Nt))
                for i in range(sn.nstations):
                    if nservers[i] > 1 and np.isfinite(nservers[i]):
                        # the server count is kept so that utilization stays
                        # normalized by c (see the 'default' branch below)
                        for j in range(Nt):
                            lldscaling[i, j] = min(j + 1, nservers[i])

                # Update sn with lldscaling
                sn.lldscaling = lldscaling
                sn.nservers = nservers.reshape(-1, 1)
                use_ld_solver = True

        elif method == 'default' and has_multiserver:
            # see _kb/06-solver-catalog.md (NC: "Multiserver -> load-dependent conversion")
            if sn.nstations == 2:
                has_delay = any(np.isinf(nservers))
                if has_delay and not has_open:
                    already_ld = hasattr(sn, 'lldscaling') and sn.lldscaling is not None and sn.lldscaling.size > 0
                    # Exact LD CoMoM enumerates the per-chain population lattice,
                    # so its cost is unbounded in N while the branch condition
                    # tests only the topology. The budget is the same one
                    # SolverCTMC uses for exact enumeration (6000 states),
                    # applied to prod(1+Nchain), and matches the MATLAB and JAR
                    # twins. Note the fallback is not free here: python's comom
                    # overflows in pfqn_ca past roughly N=1000 per chain, where
                    # MATLAB and the JAR return a value.
                    exact_lattice_max = 6000
                    njobs_all = sn.njobs.flatten() if sn.njobs is not None else np.zeros(sn.nclasses)
                    chains = np.asarray(sn.chains, dtype=float) if sn.chains is not None else np.array([])
                    lattice = 1.0
                    if chains.ndim == 2 and chains.shape[0] > 0:
                        for c in range(chains.shape[0]):
                            popc = float(np.sum([njobs_all[r] for r in np.where(chains[c, :] > 0)[0]
                                                 if r < njobs_all.size and np.isfinite(njobs_all[r])]))
                            lattice *= (1.0 + popc)
                    else:
                        for v in njobs_all:
                            if np.isfinite(v):
                                lattice *= (1.0 + float(v))
                    if lattice > exact_lattice_max:
                        self.options.method = 'comom'
                        method = 'comom'
                        line_debug("Default method: 2-station Delay+multiserver population lattice %g exceeds %g, using comom"
                                   % (lattice, exact_lattice_max), options=self.options)
                    elif self.model.has_product_form_solution() and not already_ld:
                        njobs = sn.njobs.flatten() if sn.njobs is not None else np.zeros(sn.nclasses)
                        Nt = int(np.sum(njobs[np.isfinite(njobs)]))
                        if Nt > 0:
                            lldscaling = np.ones((sn.nstations, Nt))
                            for i in range(sn.nstations):
                                if nservers[i] > 1 and np.isfinite(nservers[i]):
                                    # see _kb/06-solver-catalog.md (NC:
                                    # "Multiserver -> load-dependent conversion")
                                    for j in range(Nt):
                                        lldscaling[i, j] = min(j + 1, nservers[i])
                            sn.lldscaling = lldscaling
                            sn.nservers = nservers.reshape(-1, 1)
                            use_ld_solver = True
                        line_debug("Default method: 2-station Delay+multiserver product-form network, using exact load-dependent comomld", options=self.options)
                    else:
                        self.options.method = 'comom'
                        method = 'comom'
                        line_debug("Default method: 2-station Delay+multiserver non-product-form network, using comom", options=self.options)
            elif ms_policy == 'lld' and self.model.has_product_form_solution():
                # config.multiserver='lld' generalizes the exact load-dependent
                # lattice of the branch above to any closed product-form model,
                # under the same 6000-state enumeration budget. Off unless asked
                # for: with the shipped 'default' policy this branch never runs
                # and the model keeps Seidmann's approximation, as it always has
                lld_from_servers = nc_lld_from_nservers(sn, nservers, 6000)
                if lld_from_servers is not None:
                    sn.lldscaling = lld_from_servers
                    sn.nservers = nservers.reshape(-1, 1)
                    use_ld_solver = True
                    line_debug("Default method: config.multiserver=lld, converted multiserver "
                               "stations to load-dependent", options=self.options)
                else:
                    line_debug("Default method: config.multiserver=lld not applicable (no finite "
                               "multiserver, non-closed model, or lattice over budget), keeping "
                               "Seidmann", options=self.options)

        # Check if lldscaling is already set
        if hasattr(sn, 'lldscaling') and sn.lldscaling is not None and sn.lldscaling.size > 0:
            use_ld_solver = True

        # load-dependent normalizing-constant methods always route to ncld
        if method in ('rd', 'nrp', 'nrl', 'nre', 'comomld', 'panald'):
            use_ld_solver = True

        # see _kb/06-solver-catalog.md (NC: "Class-dependent (beta) scaling")
        use_conv_solver = False
        cd = getattr(sn, 'cdscaling', None)
        if cd is not None and len(cd) > 0 and any(x is not None for x in cd):
            use_conv_solver = True
        # Joint-dependence eta_i (non-product-form) is folded into the same
        # convolution recursion by solver_nc_conv, so route it there too.
        jd = getattr(sn, 'jdscaling', None)
        if jd is not None and len(jd) > 0 and any(x is not None for x in jd):
            use_conv_solver = True

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn, RNG
        # seeding, and sample/seed forwarding")
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            iter_max=self.options.iter_max,
            iter_tol=self.options.iter_tol,
            verbose=1 if self.options.verbose else 0,
            samples=int(self.options.samples) if getattr(self.options, 'samples', None) else 100000,
            seed=getattr(self.options, 'seed', None),
        )

        # Run the appropriate solver
        if use_conv_solver:
            from ...api.pfqn.conv import solver_nc_conv
            line_debug("class-dependent scaling detected, routing to convolution solver", options=self.options)
            self._result = solver_nc_conv(sn, handler_options)
        elif use_ld_solver:
            line_debug("Load-dependent scaling detected, routing to ncld_analyzer", options=self.options)
            self._result = solver_ncld(sn, handler_options)
        else:
            line_debug("NC method=%s, routing to nc_analyzer", self.options.method, options=self.options)
            self._result = solver_nc(sn, handler_options)

        # Extract station and class names
        self._extract_names()

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            import sys
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            runtime = self._result.runtime if hasattr(self._result, 'runtime') else 0.0
            method = self._result.method if hasattr(self._result, 'method') else 'default'
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"NC analysis [method: {method_label(self.options.method, method)}; type: {method_type('NC', method_label(self.options.method, method))}; lang: python; env: {py_version}] completed in {runtime:.6f}s.")

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn")
        if _orig_nservers is not None:
            sn.nservers = _orig_nservers
        sn.lldscaling = _orig_lldscaling

        return self

    def _run_spn_analyzer(self, sn):
        """Stationary analysis of a product-form stochastic Petri net by MDD-rec.

        The reachable set is held in a decision diagram, the product form is
        derived by spn_pf and every measure is a masked walk of that diagram.
        This is the only route SolverNC offers for a Petri net; see
        solver_nc_spn_analyzer.
        """
        from ...api.solvers.nc import solver_nc_spn_analyzer, SolverOptions as HandlerOptions
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            iter_max=self.options.iter_max,
            iter_tol=self.options.iter_tol,
            verbose=1 if self.options.verbose else 0,
        )
        handler_options.config = getattr(self.options, 'config', None)
        result = solver_nc_spn_analyzer(self.model, sn, handler_options)

        # Convert NCResult to the handler result format the base class reads
        from dataclasses import dataclass, field

        @dataclass
        class SpnResult:
            Q: np.ndarray
            U: np.ndarray
            R: np.ndarray
            T: np.ndarray
            X: np.ndarray
            C: np.ndarray
            lG: float
            STeff: np.ndarray
            it: int
            runtime: float
            method: str
            pf: dict = field(default=None)

        M, K = result.QN.shape
        self._result = SpnResult(
            Q=result.QN, U=result.UN, R=result.RN, T=result.TN,
            X=np.ravel(result.XN), C=np.ravel(result.CN),
            lG=result.lG, STeff=np.zeros((M, K)), it=result.it,
            runtime=result.runtime, method=result.method, pf=result.pf,
        )
        self._extract_names()
        if self.options.verbose:
            import sys as _sys
            from line_solver.solvers.base import print_solver_banner
            v = "%d.%d.%d" % _sys.version_info[:3]
            print_solver_banner("NC analysis [method: rec; type: exact, deterministic; "
                                "lang: python; env: %s] completed in %.6fs."
                                % (v, self._result.runtime))
        return self

    def _runCacheAnalyzer(self) -> 'SolverNC':
        """
        Run the NC cache analyzer for networks with Cache nodes.

        This method distinguishes between:
        1. Standalone cache networks (Source→Cache→Sink with no closed jobs):
           Uses solver_nc_cache_analyzer for direct cache analysis
        2. Cache+queueing networks (Cache nodes with queues and closed classes):
           Uses solver_nc_cacheqn_analyzer for iterative cache-QN analysis

        This is called automatically by runAnalyzer() when the model
        contains Cache nodes.

        Returns:
            Self for method chaining

        References:
            MATLAB: runAnalyzer.m lines 90-91, 127-128
        """
        from ...api.solvers.nc.handler import (
            solver_nc_cache_analyzer, solver_nc_cacheqn_analyzer, SolverOptions as HandlerOptions
        )
        from ...api.retrieval.analyzers import (
            solver_nc_retrieval_analyzer, solver_nc_cacheqn_retrieval_analyzer,
            has_retrieval_cache, _has_source
        )
        from ...api.sn.network_struct import NodeType

        sn = self._sn

        # Determine if this is a standalone cache network (non-reentrant cache)
        # MATLAB line 90: nclosedjobs == 0 && all(sort(nodetype) == [Source, Cache, Sink])
        is_standalone_cache = False

        # Check for no closed jobs
        nclosedjobs = sn.nclosedjobs if hasattr(sn, 'nclosedjobs') else 0
        if nclosedjobs == 0:
            # Check if node types are exactly Source, Cache, Sink
            if sn.nodetype is not None and len(sn.nodetype) == 3:
                node_types_sorted = sorted(int(nt) for nt in sn.nodetype)
                expected_types = sorted([int(NodeType.SOURCE), int(NodeType.CACHE), int(NodeType.SINK)])
                if node_types_sorted == expected_types:
                    is_standalone_cache = True

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn, RNG
        # seeding, and sample/seed forwarding")
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            iter_max=self.options.iter_max,
            iter_tol=self.options.iter_tol,
            verbose=1 if self.options.verbose else 0,
            samples=int(self.options.samples) if getattr(self.options, 'samples', None) else 100000,
            seed=getattr(self.options, 'seed', None),
        )

        if has_retrieval_cache(sn):
            # Delayed-hit cache with a retrieval system: open (Source) -> product-form
            # retrieval analyzer; closed integrated -> da_cacheqn_retrieval driver.
            if _has_source(sn):
                cache_result = solver_nc_retrieval_analyzer(sn, handler_options)
            else:
                cache_result = solver_nc_cacheqn_retrieval_analyzer(sn, handler_options)
        elif is_standalone_cache:
            # Standalone cache network: use direct cache analyzer
            cache_result = solver_nc_cache_analyzer(sn, handler_options)
        else:
            # Cache+queueing network: use iterative cache-QN analyzer
            cache_result = solver_nc_cacheqn_analyzer(sn, handler_options)

        # Convert cache result to standard NC result format
        # Create a result object compatible with the standard NC result
        from dataclasses import dataclass

        @dataclass
        class CacheResultAdapter:
            Q: np.ndarray
            U: np.ndarray
            R: np.ndarray
            T: np.ndarray
            X: np.ndarray
            lG: float
            STeff: np.ndarray
            it: int
            runtime: float
            method: str
            pij: np.ndarray  # Cache-specific: item probabilities
            hitprob: np.ndarray  # Hit probabilities per cache per class
            missprob: np.ndarray  # Miss probabilities per cache per class

        M = self._sn.nstations
        K = self._sn.nclasses

        # Initialize with zeros for stations that don't have metrics
        Q = cache_result.QN if cache_result.QN is not None else np.zeros((M, K))
        U = cache_result.UN if cache_result.UN is not None else np.zeros((M, K))
        R = cache_result.RN if cache_result.RN is not None else np.zeros((M, K))
        T = cache_result.TN if cache_result.TN is not None else np.zeros((M, K))

        # Handle different result types (standalone vs cacheqn)
        pij = getattr(cache_result, 'pij', None)
        hitprob = getattr(cache_result, 'hitprob', None)
        missprob = getattr(cache_result, 'missprob', None)
        it_count = getattr(cache_result, 'it', 1)

        self._result = CacheResultAdapter(
            Q=Q,
            U=U,
            R=R,
            T=T,
            X=cache_result.XN,
            lG=cache_result.lG,
            STeff=np.zeros((M, K)),
            it=it_count,
            runtime=cache_result.runtime,
            method=cache_result.method,
            pij=pij,
            hitprob=hitprob,
            missprob=missprob
        )

        # Store cache-specific results
        self._cache_result = cache_result

        # Copy updated visits from cache analyzer to self._sn
        # (MATLAB: self.model.refreshStruct(true); sn = self.model.sn;)
        if hasattr(cache_result, 'visits') and cache_result.visits is not None:
            sn.visits = cache_result.visits
        if hasattr(cache_result, 'nodevisits') and cache_result.nodevisits is not None:
            sn.nodevisits = cache_result.nodevisits

        # Extract station and class names
        self._extract_names()

        # A closed cache-retrieval solve relabels the Cache to a ClassSwitch;
        # the cache node index is carried in cache_result.cache_idx.
        retr_cache_idx = getattr(cache_result, 'cache_idx', None)
        for ind in range(sn.nnodes):
            if sn.nodetype is not None and ind < len(sn.nodetype):
                is_cache_node = (sn.nodetype[ind] == NodeType.CACHE) or (ind == retr_cache_idx)
                if is_cache_node and ind in sn.nodeparam:
                    cache_param = sn.nodeparam[ind]

                    # Try to get hit/miss probs from nodeparam (standalone cache)
                    actualhitprob = getattr(cache_param, 'actualhitprob', None)
                    actualmissprob = getattr(cache_param, 'actualmissprob', None)

                    # For cacheqn results, always use the hitprob from the NC solver
                    # (a previous solver like SSA may have set actualhitprob, which would be stale)
                    if hitprob is not None:
                        # hitprob shape is (ncaches, K) - find index of this cache
                        cache_indices = []
                        for i in range(sn.nnodes):
                            if i < len(sn.nodetype) and (sn.nodetype[i] == NodeType.CACHE or i == retr_cache_idx):
                                cache_indices.append(i)
                        if ind in cache_indices:
                            cache_idx = cache_indices.index(ind)
                            if cache_idx < hitprob.shape[0]:
                                actualhitprob = hitprob[cache_idx, :]
                                actualmissprob = missprob[cache_idx, :]

                    if actualhitprob is not None:
                        # Update sn.nodeparam with actual hit/miss probs for sn_get_node_tput_from_tput
                        cache_param.actualhitprob = actualhitprob
                        cache_param.actualmissprob = actualmissprob

                        # Delayed-hit retrieval: per-class delayed-hit fraction and
                        # per-list hit fractions (None/absent for plain caches).
                        dhp = getattr(cache_result, 'delayedprob', None)
                        hpl = getattr(cache_result, 'hitproblist', None)
                        ipb = getattr(cache_result, 'itemprob', None)
                        # the integrated cacheqn branch reports one law PER CACHE
                        if isinstance(ipb, (list, tuple)):
                            all_caches = [i for i in range(sn.nnodes)
                                          if i < len(sn.nodetype)
                                          and (sn.nodetype[i] == NodeType.CACHE or i == retr_cache_idx)]
                            cidx = all_caches.index(ind) if ind in all_caches else -1
                            ipb = ipb[cidx] if 0 <= cidx < len(ipb) else None
                        lcs = getattr(cache_result, 'listcost', None)
                        if lcs is not None:
                            cache_param.actuallistcost = np.asarray(lcs)
                        if dhp is not None:
                            cache_param.actualdelayedhitprob = np.asarray(dhp)[0, :]
                        if hpl is not None:
                            cache_param.actualhitproblist = np.asarray(hpl)
                        if ipb is not None:
                            cache_param.actualitemprob = np.asarray(ipb)

                        # Also set result on Cache node in model
                        if hasattr(self, 'model') and hasattr(self.model, '_nodes'):
                            cache_node = self.model._nodes[ind]
                            if hasattr(cache_node, 'set_result_hit_prob'):
                                cache_node.set_result_hit_prob(actualhitprob)
                            if actualmissprob is not None and hasattr(cache_node, 'set_result_miss_prob'):
                                cache_node.set_result_miss_prob(actualmissprob)
                            if dhp is not None and hasattr(cache_node, 'set_result_delayed_hit_prob'):
                                cache_node.set_result_delayed_hit_prob(np.asarray(dhp)[0, :])
                            if hpl is not None and hasattr(cache_node, 'set_result_hit_prob_list'):
                                cache_node.set_result_hit_prob_list(np.asarray(hpl))
                            if ipb is not None and hasattr(cache_node, 'set_result_item_prob'):
                                cache_node.set_result_item_prob(np.asarray(ipb))
                            if lcs is not None and hasattr(cache_node, 'set_result_list_cost'):
                                cache_node.set_result_list_cost(np.asarray(lcs))
                            # Delayed-hit retrieval: expected latency Z (NaN for non-retrieval)
                            el = getattr(cache_result, 'expected_latency', None)
                            if el is not None and hasattr(cache_node, 'set_result_residt'):
                                cache_node.set_result_residt(np.asarray(el)[0, :])

        # THE CACHE ALGORITHM MUST BE NAMED. Every other NC path prints this
        # banner; the cache branch returned silently, so a default solve served
        # the SPM approximation (`default/spm`) with nothing on screen to say so
        # -- and SPM is not exact on a small cache, where the `exact` recursion
        # is both available and cheap. C++ already prints `default/spm` here.
        if self.options.verbose:
            import sys
            py_version = f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}"
            cache_method = method_label(self.options.method, self._result.method)
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"NC analysis [method: {cache_method}; "
                  f"type: {method_type('NC', cache_method)}; lang: python; "
                  f"env: {py_version}] completed in {self._result.runtime:.6f}s.")

        return self

    def _extract_names(self):
        """Extract station and class names from network struct."""
        if self._sn is not None:
            # Use station names, not node names (NC operates on stations, not nodes)
            # Nodes include non-station elements like ClassSwitch, Router, etc.
            if hasattr(self._sn, 'stationnames') and self._sn.stationnames:
                self.station_names = list(self._sn.stationnames)
            elif hasattr(self._sn, 'nodenames') and self._sn.nodenames and hasattr(self._sn, 'stationToNode'):
                # Map station indices to node names
                self.station_names = []
                for i in range(self._sn.nstations):
                    if i < len(self._sn.stationToNode):
                        node_idx = int(self._sn.stationToNode[i])
                        if node_idx < len(self._sn.nodenames):
                            self.station_names.append(self._sn.nodenames[node_idx])
                        else:
                            self.station_names.append(f'Station{i}')
                    else:
                        self.station_names.append(f'Station{i}')
            else:
                self.station_names = [f'Station{i}' for i in range(self._sn.nstations)]
            self.class_names = list(self._sn.classnames) if hasattr(self._sn, 'classnames') and self._sn.classnames else \
                              [f'Class{i}' for i in range(self._sn.nclasses)]
        else:
            self.station_names = []
            self.class_names = []

    # =========================================================================
    # Table Output
    # =========================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self._ensureAvgResults()
        self._cap_unstable_open_util()

        nstations = self._result.Q.shape[0]
        nclasses = self._result.Q.shape[1]

        # Compute residence times from response times using visit ratios
        if self._sn is not None and self._sn.visits:
            WN = sn_get_residt_from_respt(self._sn, self._result.R, None)
        else:
            WN = self._result.R.copy()

        # Compute proper arrival rates (sets Source ArvR = 0)
        AN = sn_get_arvr_from_tput(self._sn, self._result.T)

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                station_name = self.station_names[i] if i < len(self.station_names) else f'Station{i}'
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                rows.append({
                    'Station': station_name,
                    'JobClass': class_name,
                    'QLen': self._result.Q[i, r],
                    'Util': self._result.U[i, r],
                    'RespT': self._result.R[i, r],
                    'ResidT': WN[i, r],
                    'ArvR': AN[i, r],
                    'Tput': self._result.T[i, r],
                })

        df = pd.DataFrame(rows)

        # Filter out all-zero rows
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
        if self._result is None:
            self._ensureAvgResults()
        return self._result.Q.copy()

    def _cap_unstable_open_util(self) -> None:
        """
        Cap the reported utilization of unstable open queueing stations at 1.0.

        Delegates to the shared NC/MVA rule: a finite-server open station with
        offered load rho >= 1 is fully saturated and its utilization is reported
        as 1.0 (split across classes by offered load), with a single instability
        warning. Source and delay stations are excluded. Idempotent.
        """
        if getattr(self, '_unstable_util_capped', False):
            return
        self._unstable_util_capped = True
        res = self._result
        if res is None or getattr(res, 'U', None) is None or getattr(res, 'T', None) is None \
                or self._sn is None:
            return
        from ...api.sn.transforms import cap_unstable_open_util
        UN, any_unstable = cap_unstable_open_util(res.U, res.T, self._sn)
        if any_unstable:
            res.U[:, :] = UN
            line_warning("SolverNC", "The model has unstable queues "
                         "(utilization >= 1); station utilization is reported "
                         "capped at 1.0 and other metrics may grow unbounded.")

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K)."""
        if self._result is None:
            self._ensureAvgResults()
        self._cap_unstable_open_util()
        return self._result.U.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K)."""
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
        """Get average waiting times (M x K)."""
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
        if self._result is None:
            self._ensureAvgResults()
        return self._result.T.copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K).

        Uses routing matrix to compute proper arrival rates.
        Source stations have arrival rate = 0.
        """
        if self._result is None:
            self._ensureAvgResults()
        return sn_get_arvr_from_tput(self._sn, self._result.T)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times (cycle times) per chain (nchains,).

        Returns chain-level response times matching MATLAB/Java implementation.
        Uses the completes flag to determine which classes contribute to chain throughput.

        Note:
            For closed chains: uses Little's Law CNchain = nJobsChain / XNchain
            For open chains: weighted sum of class response times
        """
        CN, XN = self._computeChainMetrics()
        return CN

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs per chain (nchains,).

        Returns chain-level throughputs matching MATLAB/Java implementation.
        Uses the completes flag to determine which classes contribute to chain throughput.
        """
        CN, XN = self._computeChainMetrics()
        return XN

    # _computeChainMetrics is inherited from NetworkSolver (base.py): the
    # faithful MATLAB @NetworkSolver/getAvgSys.m chain-based port, shared with
    # the CTMC solver so getAvgSys agrees across codebases.

    # =========================================================================
    # NC-Specific Methods
    # =========================================================================

    def getNormalizingConstant(self) -> float:
        """
        Get the normalizing constant G.

        Returns:
            float: The normalizing constant (not log)
        """
        if self._result is None:
            self._ensureAvgResults()

        lG = self._result.lG
        if np.isfinite(lG):
            return np.exp(lG)
        return 0.0

    def getLogNormalizingConstant(self) -> float:
        """
        Get the log normalizing constant log(G).

        Returns:
            float: log(G)
        """
        if self._result is None:
            self._ensureAvgResults()
        return self._result.lG

    def getProbNormConstAggr(self) -> float:
        """
        Get the log normalizing constant (alias).

        Returns:
            float: log(G)
        """
        return self.getLogNormalizingConstant()

    def getAvgBusyPeriod(self, stations, n=1):
        """Mean busy period of order n for a set of stations.

        The busy period of order n runs from the instant a job entering the set
        finds n-1 jobs in it up to the next instant when fewer than n remain
        (H. Daduna, "Busy Periods for Subnetworks in Stochastic Networks: Mean
        Value Analysis", J. ACM 35(3), 1988). Exact on the single-chain
        product-form class and, by the insensitivity of Section 5 of that paper,
        dependent on the service processes only through their mean rates.

        Args:
            stations: stations forming the subnetwork, as objects, names or
                zero-based station indexes.
            n: busy period order or sequence of orders.

        Returns:
            (b, lG, lH): mean busy period duration(s) and the log normalizing
            constants of the subnetwork and of its complement.
        """
        from line_solver.api.pfqn.busyp import pfqn_busyp
        from line_solver.api.sn.demands import sn_get_demands_chain
        from line_solver.api.sn.transforms import sn_rt_stations
        from line_solver.api.sn.network_struct import NodeType

        sn = self.model.getStruct()
        if int(sn.nchains) > 1:
            raise ValueError(
                'The busy period of a subnetwork is defined for single-chain models '
                'only. Section 5 of Daduna (1988) sketches the multichain extension, '
                'which is not implemented.')

        subnet = []
        if not isinstance(stations, (list, tuple, np.ndarray)):
            stations = [stations]
        for st in stations:
            if isinstance(st, (int, np.integer)):
                subnet.append(int(st))
            else:
                subnet.append(int(self.model.get_station_index(st)) - 1)

        demands = sn_get_demands_chain(sn)
        STchain = np.asarray(demands.STchain, dtype=float)[:, 0]
        Vchain = np.asarray(demands.Vchain, dtype=float)[:, 0]
        Nchain = np.asarray(demands.Nchain, dtype=float).flatten()
        rtst, Vst = sn_rt_stations(sn)
        M = int(sn.nstations)
        K = int(sn.nclasses)

        # station-to-station routing of the chain: the class-level probabilities
        # weighted by the class visits, which is exact because it is a flow balance
        Pst = np.zeros((M, M))
        for i in range(M):
            for j in range(M):
                blk = rtst[i * K:(i + 1) * K, j * K:(j + 1) * K]
                Pst[i, j] = float(Vst[i, :] @ blk.sum(axis=1))
        Vtot = Vst.sum(axis=1)
        nz = Vtot > 0
        Pst[nz, :] = Pst[nz, :] / Vtot[nz][:, None]

        lld = getattr(sn, 'lldscaling', None)
        nservers = np.asarray(sn.nservers, dtype=float).flatten()

        def ratefun(j, kvec):
            # same precedence as solver_ncld: the infinite server first, then the
            # declared load-dependent scaling, then the multiserver staircase. A
            # scaling table shorter than kvec keeps its last entry.
            kvec = np.asarray(kvec, dtype=int)
            if np.isinf(nservers[j]):
                scale = kvec.astype(float)
            elif lld is not None and np.size(lld) > 0:
                cols = np.shape(lld)[1]
                scale = np.asarray(lld)[j, np.minimum(kvec, cols) - 1]
            else:
                scale = np.minimum(kvec, nservers[j])
            return scale / STchain[j]

        if not np.isfinite(Nchain[0]):
            # the Source is not a node of the Jackson network of the paper: its
            # outflow is the external stream gamma
            source = None
            for i in range(M):
                if int(sn.nodetype[int(sn.stationToNode[i])]) == NodeType.SOURCE:
                    source = i
                    break
            if source is None:
                raise ValueError('An open model must own a Source station.')
            if source in subnet:
                raise ValueError('The Source cannot belong to the subnetwork.')
            rates = np.asarray(sn.rates, dtype=float)[source, :]
            lam = float(np.nansum(rates))
            keep = [i for i in range(M) if i != source]
            remap = {v: t for t, v in enumerate(keep)}
            alpha = lam * Vchain[keep] / Vchain[source]
            gamma = lam * Pst[source, keep]
            P = Pst[np.ix_(keep, keep)]
            mu = lambda j, kvec: ratefun(keep[j], kvec)
            return pfqn_busyp(alpha, mu, P, np.inf,
                              [remap[s] for s in subnet], n, gamma)

        return pfqn_busyp(Vchain, ratefun, Pst, Nchain[0], subnet, n)

    def getEffectiveServiceTimes(self) -> np.ndarray:
        """
        Get effective service times.

        Returns:
            np.ndarray: Effective service times (M x K)
        """
        if self._result is None:
            self._ensureAvgResults()
        if self._result.STeff is not None:
            return self._result.STeff.copy()
        return np.array([])

    def getIterationCount(self) -> int:
        """Get the number of iterations used."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.it

    def getRuntime(self) -> float:
        """Get the solver runtime in seconds."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.runtime

    def getMethodUsed(self) -> str:
        """Get the method actually used (may differ from requested)."""
        if self._result is None:
            self._ensureAvgResults()
        return self._result.method

    # =========================================================================
    # CDF and Percentile Methods
    # =========================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Get response time CDF using exponential approximation.

        Ports MATLAB's `@SolverNC/getCdfRespT.m`: the EXACT product-form
        sojourn law at FCFS stations via `pfqn_stdf`, not the base-class
        exponential approximation. Delay stations carry their own service CDF.

        `options.config['algorithm']` selects 'exact' (`pfqn_stdf`, default) or
        'rd' (`pfqn_stdf_heur`), as the reference does.

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys. Only FCFS and
            delay stations are populated, which is exactly the set the
            reference fills; a PS or LCFS queue has no entry.

        Raises:
            ValueError: on an open class. The tagged-job passage time is
                defined on a CLOSED network.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import cdf_respt_via_jar
            return cdf_respt_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import cdf_respt_via_cpp
            return cdf_respt_via_cpp(self)
        if self._result is None:
            self._ensureAvgResults()

        from ...api.pfqn.stdf import pfqn_stdf, pfqn_stdf_heur
        from ...api.sn.transforms import sn_get_product_form_params
        from ...api.mam import map_cdf
        sn = self.model.getStruct()
        pf = sn_get_product_form_params(sn)
        D, Npop, Z, S = pf.D, np.ravel(pf.N), pf.Z, pf.S

        # sn.sched is a DICT keyed by station index in python, not an array.
        # Compare by NAME: there are TWO distinct SchedStrategy enum classes in
        # the tree (line_solver.lang.base, which sn.sched carries, and
        # line_solver.constants, which the top-level name re-exports), and `==`
        # between them is FALSE for every member -- silently, so an identity
        # comparison here reports "no FCFS stations" on a model that has one.
        # Same family as the ProcessType rule in CLAUDE.md: compare names.
        def _sched_name(v):
            return getattr(v, 'name', str(v)).upper()

        sched = sn.sched
        nstations = len(sched)
        nondelay = [i for i in range(nstations) if _sched_name(sched[i]) != 'INF']
        fcfs_node_ids = [i for i in range(nstations) if _sched_name(sched[i]) == 'FCFS']
        delay_node_ids = [i for i in range(nstations) if _sched_name(sched[i]) == 'INF']
        # index of each FCFS station WITHIN the non-delay subset, which is the
        # space pfqn_stdf indexes L, S and rates in
        fcfs_nodes = [nondelay.index(i) for i in fcfs_node_ids]

        if not fcfs_nodes:
            line_warning("getCdfRespT", "getCdfRespT applies only to FCFS nodes.")
            return []

        if np.any(np.isinf(Npop)):
            raise ValueError(
                "The tagged-job sojourn-time law is defined on a CLOSED "
                "network; this model has an open class. Use SolverFluid or "
                "SolverMAM for a response-time distribution with open classes.")

        rates = np.atleast_2d(np.asarray(sn.rates, dtype=float))
        # Horizon uses the STATION-space ids, since sn.rates is indexed by
        # station; the reference does the same (getCdfRespT.m:30).
        T = float(np.max(np.sum(Npop) * np.mean(1.0 / rates[fcfs_node_ids, :], axis=1)))
        tset = np.logspace(0.0, 2.0 * np.log10(T), 100)

        # Non-delay rows, so rates lives in the SAME index space as D, S and
        # fcfs_nodes. MATLAB used to pass the FCFS rows only and crashed on
        # Delay+PS+FCFS ("Index in position 1 exceeds array bounds", register
        # row N15); that is fixed in getCdfRespT.m now, so the two agree and
        # this is no longer a divergence.
        rates_nd = rates[nondelay, :]

        algorithm = 'exact'
        cfg = getattr(self.options, 'config', None)
        if isinstance(cfg, dict):
            algorithm = cfg.get('algorithm', 'exact')
        elif cfg is not None and hasattr(cfg, 'algorithm'):
            algorithm = getattr(cfg, 'algorithm') or 'exact'

        fn = pfqn_stdf_heur if algorithm == 'rd' else pfqn_stdf
        RDout = fn(D, Npop, Z, np.asarray(S, dtype=int).ravel(),
                   np.asarray(fcfs_nodes, dtype=int), rates_nd, tset)

        RD = []
        for pos, ist in zip(fcfs_nodes, fcfs_node_ids):
            for r in range(sn.nclasses):
                arr = RDout.get((pos, r))
                if arr is None:
                    continue
                RD.append({'station': ist + 1, 'class': r + 1,
                           't': np.real(arr[:, 1]), 'p': np.real(arr[:, 0])})

        # A delay station has no queueing, so the passage time IS its service
        # time and the reference emits that CDF directly.
        for ist in delay_node_ids:
            for r in range(sn.nclasses):
                D0, D1 = self._decode_proc(sn, ist, r)
                if D0 is None:
                    continue
                RD.append({'station': ist + 1, 'class': r + 1,
                           't': tset, 'p': np.ravel(map_cdf(D0, D1, tset))})
        return RD

    @staticmethod
    def _decode_proc(sn, ist, r):
        """sn.proc[i][r] carries THREE shapes in python: a [D0,D1] sequence, a
        dict with 'D0'/'D1', and the exponential shorthand {'rate': r}."""
        proc = getattr(sn, 'proc', None)
        if proc is None:
            return None, None
        try:
            p = proc[ist][r]
        except (KeyError, IndexError, TypeError):
            return None, None
        if p is None:
            return None, None
        if isinstance(p, (list, tuple)) and len(p) >= 2:
            return np.asarray(p[0], dtype=float), np.asarray(p[1], dtype=float)
        if isinstance(p, dict):
            if 'D0' in p and 'D1' in p:
                return np.asarray(p['D0'], dtype=float), np.asarray(p['D1'], dtype=float)
            if 'rate' in p:
                lam = float(p['rate'])
                if not np.isfinite(lam) or lam <= 0:
                    return None, None
                return np.array([[-lam]]), np.array([[lam]])
        return None, None

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None,
        method: str = 'default'
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """
        Extract percentiles from response time distribution.

        Args:
            percentiles: List of percentiles (0-100). Default: [10, 25, 50, 75, 90, 95, 99]
            jobclass: Optional class filter (1-based)

        Returns:
            Tuple of (percentile_list, percentile_table)
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

        # Extract from the CDF, as @NetworkSolver/getPerctRespT.m:70 does.
        # This MUST read getCdfRespT rather than re-derive an exponential fit:
        # SolverNC overrides getCdfRespT with the EXACT product-form sojourn
        # law, so a separate exponential fit here would make the two accessors
        # disagree with each other on the same model.
        RD = self.getCdfRespT()
        for entry in RD:
            i = int(entry['station']) - 1
            r = int(entry['class']) - 1
            if jobclass is not None and (r + 1) != jobclass:
                continue
            times = np.ravel(np.asarray(entry['t'], dtype=float))
            probs = np.ravel(np.asarray(entry['p'], dtype=float))
            if times.size < 2 or probs.size != times.size:
                continue
            order = np.argsort(times)
            times, probs = times[order], probs[order]
            # np.interp needs a non-decreasing x; the CDF is the x here because
            # we invert F to get t(p), so accumulate a running max to make it
            # monotone against round-off before inverting.
            probs = np.maximum.accumulate(probs)
            perc_values = np.interp(percentiles_normalized, probs, times,
                                    left=times[0], right=times[-1])

            PercRT.append({
                'station': i + 1,
                'class': r + 1,
                'percentiles': percentiles.tolist(),
                'values': np.asarray(perc_values).tolist(),
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
    # Introspection Methods
    # =========================================================================

    def listValidMethods(self) -> List[str]:
        """List valid solution methods.

        Returns:
            List of valid method names for NC solver:
            - 'default': Auto-select based on problem size
            - 'exact', 'ca': Exact convolution algorithm
            - 'imci': Importance sampling Monte Carlo integration
            - 'ls': Linearizer method
            - 'le': Logistic expansion (Cas17)
            - 'ble': Logistic expansion with an empirical correction to LE
            - 'mmint2': Gauss-Legendre quadrature
            - 'gleint': Gauss-Legendre integration
            - 'pana': PANACEA asymptotic expansion (load-independent)
            - 'panald': PANACEA asymptotic expansion (load-dependent)
            - 'kt': Knessl-Tier expansion
            - 'bkt': Knessl-Tier expansion with the Stirling-remainder correction (BKT)
            - 'lekt': the estimator 'ble' and 'bkt' both compute, on the cheaper side
            - 'sampling': Monte Carlo sampling
            - 'propfair': Proportionally fair allocation
            - 'divdiff': divided-difference closed form (Casale, SIGMETRICS 2017),
              no think time; load-dependent rates go through the limited
              load-dependent kernel of Casale-Harrison-Ong (Perform. Eval. 2021)
            - 'rgf': Recursion by generating functions (single-class, grouped stations)
            - 'ger': Gerasimov residue closed form (free in the eliminated
              class populations, costly in the class count)
            - 'comom': Conditional moments
            - 'comomld': Conditional moments, load-dependent
            - 'cub': Controllable upper bound
            - 'gm': Grundmann-Moeller cubature (alias of 'cub')
            - 'rd': Reduction heuristic
            - 'nrl': Norlund-Rice Logit approximation
            - 'nrp': Norlund-Rice Probit approximation
            - 'nre': Norlund-Rice saddle-tilted Edgeworth approximation
            - 'ms': Manjunath-Sikdar transform of the loss-network analyzer,
              the only place it is admissible
            - 'sdr', 'sdr.mva': Krzesinski state-dependent routing, the eq. (16)
              enumeration and its Section 4 MVA arm
            - 'morrison': heavy-usage expansion for a closed think+DPS network
            - 'rec': memoised decision-diagram walk of the reachable set, for
              product-form Petri nets and loss networks
            - 'mcmc': Chen-O'Cinneide regularization, a Markov chain Monte Carlo
              estimator of the throughput ratios and the queue lengths

        'ms', 'sdr' and 'sdr.mva' are DISPATCHED here (runAnalyzer and
        api/solvers/nc/analyzers.py) and were missing from this list, so the
        shared gate in NetworkSolver.runAnalyzerChecks refused three methods
        this solver implements and the other three codebases advertise.
        """
        return [
            # 'rayint' and 'spm' both name the SPM saddle point on a cache, which
            # serves cache_spm_size once the items carry storage costs. On a retrieval
            # model 'rayint' is instead the ray/WKB delayed-hit expansion, admissible
            # only with an infinite-server fetch system; solver_nc_retrieval_analyzer
            # branches on the token and warns and falls back to 'exact' anywhere else.
            # 'divdiff' is the divided-difference closed form of Casale (SIGMETRICS
            # 2017); it needs no think time, since a delay would ask for the integral
            # form of Cor. 3.4, and pfqn_nc refuses one by name. Load-dependent rates
            # ARE served: pfqn_ncld substitutes the limited load-dependent kernel of
            # Casale-Harrison-Ong (Perform. Eval. 2021), Thm. 1, and reports itself as
            # 'divdiff.ld/...'.
            'default', 'exact', 'divdiff', 'rayint', 'spm', 'ms', 'erlangfp', 'mci', 'ca', 'clw',
            'imci', 'ls',
            'le', 'ble', 'aghq', 'mmint2', 'gleint', 'pana', 'panald',
            'kt', 'bkt', 'lekt', 'bk', 'bkue', 'lc', 'lc.ue', 'sampling', 'is',
            # Chen-O'Cinneide regularization; a Markov chain Monte Carlo estimator of
            # the throughput RATIOS G(N-e_r)/G(N), which supplies no constant of its own
            'mcmc',
            'propfair', 'comom', 'comomld', 'cub', 'gm', 'rgf', 'ger',
            'rd', 'nrl', 'nrp', 'nre', 'mem', 'sdr', 'sdr.mva',
            # 'morrison' is the heavy-usage asymptotic expansion of the generating
            # function for a closed think+DPS network (npfqn_dps_morrison,
            # solver_nc_dps_analyzer). It is the DEFAULT on that shape and inadmissible
            # anywhere else, where runAnalyzer refuses it: nothing else in NC can see
            # the DPS weights. Non-product-form, so it returns no lG.
            'morrison',
            # 'rec' is the MDD-rec route for product-form Petri nets: it is the only
            # method solver_nc_spn_analyzer serves and the only one admissible on a
            # net, and on a loss network it is the exact normalizing constant without
            # the residue transform's integrality demand. runAnalyzer already
            # dispatched it (solver_nc.py:664) while this list withheld the name, so
            # NetworkSolver.runAnalyzerChecks refused a method the solver implements.
            'rec',
        ]

    def isStochasticMethod(self, method):
        """NC is deterministic except for the Monte Carlo integration methods
        (mci/imci), logistic sampling (ls), the importance sampling method (is),
        the Chen-O'Cinneide Markov chain Monte Carlo method (mcmc), and the
        sampling method, whose estimates depend on the random seed. Method names
        are tokenized so that runtime-resolved names such as 'default/imci' and
        prefixed names such as 'nc.ls' classify correctly.
        """
        if not method:
            return False
        tokens = re.split(r'[./]', str(method).lower())
        return any(tok in ('mci', 'imci', 'ls', 'sampling', 'is', 'mcmc') for tok in tokens)

    is_stochastic_method = isStochasticMethod

    def resolveMethod(self, options):
        """Feature-driven resolution of method='default': an open network with
        non-Markovian (non-unit SCV) variability within the MEM feature set is
        solved by the Maximum Entropy Method by default, since the
        normalizing-constant path would silently exponentialize it. Mirrors the
        dispatch below in runAnalyzer and the MATLAB SolverNC.resolveMethod."""
        method = getattr(options, 'method', 'default')
        if method != 'default':
            return method
        sn = getattr(self, '_sn', None)
        if sn is None:
            return method
        try:
            from ...api.me import solver_nc_mem_supports
            if solver_nc_mem_supports(sn)[0]:
                scv = getattr(sn, 'scv', None)
                if scv is not None:
                    scv = np.asarray(scv, dtype=float)
                    scvv = scv[np.isfinite(scv)]
                    if scvv.size and np.any(np.abs(scvv - 1.0) > 1e-8):
                        return 'mem'
        except Exception:
            pass
        return method

    def getMethodFeatureSet(self, method):
        """Per-method feature deltas applied to the base NC envelope, as a set of
        feature-name strings.

        ONLY THE RESTRICTIONS A FEATURE NAME CAN CARRY LIVE HERE. A feature set
        declares what the method ACCEPTS, so it can refuse a model for HAVING a
        construct and never for lacking one: "closed population only" and "no
        think time" are expressible by dropping OpenClass and SchedStrategy_INF,
        while "requires a cache" or "requires a loss network" are not and belong
        to nc_method_refusal, which supportsModelMethod consults next. Mirrors
        the MATLAB/JAR SolverNC.getMethodFeatureSet and the C++ nc_feature_set.
        """
        feats = {name for name, on in SolverNC.getFeatureSet().list.items() if on}
        m = str(method or 'default').lower()
        if m == 'divdiff':
            # The divided-difference closed form of Casale (SIGMETRICS 2017),
            # Eqs. (15)-(16), covers load-independent queues; a think time would
            # ask for the integral form of Cor. 3.4, which is not implemented, so
            # pfqn_nc and pfqn_ncld both refuse one by name. An infinite server is
            # where a think time comes from, so the envelope drops it.
            feats.discard('SchedStrategy_INF')
        elif m in ('rd', 'nrp', 'nrl', 'nre', 'comomld', 'panald'):
            # The load-dependent normalizing-constant evaluators are reached by
            # solver_ncld only on its CLOSED branch, where pfqn_ncld reads the
            # method name. An open chain sends the model to the mixed route
            # (pfqn_ncldmx), which never reads it, so every one of these names
            # silently became 'ncldmx'.
            feats.discard('OpenClass')
        elif m == 'is':
            # The sample-an-ordering estimator of pfqn_is integrates over a closed
            # population simplex; there is no open-class form of it, and
            # solver_nc_analyzer refuses one by name. Use 'sampling'
            # (pfqn_mci/pfqn_ls) for an open or mixed model.
            feats.discard('OpenClass')
        # MULTISERVER (registry name since 2026-09-05): the divided-difference
        # closed form of 'divdiff' covers load-independent single-server
        # queues, and nc_method_refusal keeps wording why (a c-server station
        # enters the constant as Seidmann's surrogate delay); every other route
        # folds the count into its own kernel.
        if m == 'divdiff':
            feats.discard('MultiServer')
        # FINITECAPACITY (registry name since 2026-09-05) is NOT in the base
        # envelope: the product-form routes solve a buffer away, which is what
        # the binding-capacity gate refuses. Two arms honour one: 'mem'
        # represents it as a GE/GE/c/0;N queue (solver_nc_mem_supports), and the
        # single-station M/M/1/K with tail drop is solved in closed form
        # (qsys_mm1k_loss) under 'default' and 'exact'. The shape half of each
        # rule stays structural.
        if m in ('mem', 'default', 'exact'):
            feats.add('FiniteCapacity')
        return feats

    get_method_feature_set = getMethodFeatureSet

    def supportsModelMethod(self, method):
        """Method-aware gate. MEM (Kouvatsos maximum entropy) has structural
        applicability rules beyond a flat feature set (open-only, no class
        switching, non-priority scheduling); delegate to solver_nc_mem_supports,
        which returns a precise reason. Every other method is gated on its own
        per-method feature set and then on nc_method_refusal, the single copy of
        the structural rules the analyzer enforces, so that a pair this gate
        offers is a pair the run accepts. The single-Delay DROP loss-network
        exception handled by solver_nc_lossn_analyzer is preserved."""
        # Discrete-time (slotted) route: a finite buffer on a Bernoulli server
        # is the loss system of Daduna's corollary 2.8, which
        # solver_nc_dt_analyzer solves exactly, so the structural capacity gate
        # below must not fire. nc_is_dt_model performs the real admissibility
        # check and reports a precise reason.
        from .solver_nc_dt_analyzer import is_slotted, nc_is_dt_model
        if is_slotted(self.options):
            sn = getattr(self, '_sn', None)
            if sn is None and getattr(self, 'model', None) is not None:
                sn = self.model.getStruct()
            dt = nc_is_dt_model(sn, self.options)
            if dt['kind'] != 'none':
                return True, ''
            return False, ('options.config slotted is set but %s' % dt['reason'])
        if method == 'mem':
            from ...api.me import solver_nc_mem_supports
            sn = getattr(self, '_sn', None)
            if sn is None:
                sn = self.model.getStruct()
            ok, reason = solver_nc_mem_supports(sn)
            return bool(ok), (reason or '')
        model = getattr(self, 'model', None)
        # The per-method envelope, not the flat one: 'divdiff' drops
        # SchedStrategy_INF and the closed-population evaluators drop OpenClass,
        # which is how "no think time" and "closed only" are said in the registry.
        ok, reason = NetworkSolver.supportsModelMethod(self, method)
        if not ok and _is_lossn_fcr_case(getattr(self, '_sn', None)):
            # NC solves the OPEN single-Delay loss network exactly (Erlang fixed
            # point); the flat envelope cannot express that split.
            ok, reason = True, ''
        if not ok:
            return False, (reason or 'Some features are not supported by the NC solver.')
        # see _kb/06-solver-catalog.md ("Finite capacity gate (MVA and NC)")
        if model is not None and hasattr(model, 'getStruct'):
            ok, reason = NetworkSolver.checkBindingCapacity(model, 'SolverNC')
            if not ok:
                return ok, reason
            ok, reason = SolverNC.supportsExactness(model, method)
            if not ok:
                return ok, reason
        # The structural per-method rules the feature registry cannot name: which
        # route a method has on THIS model, and whether it exists at all.
        # nc_method_refusal is the single copy of them, asked here and by
        # runAnalyzer, so the report and the run cannot disagree.
        from .nc_method_refusal import nc_method_refusal
        sn = getattr(self, '_sn', None)
        if sn is None and model is not None and hasattr(model, 'getStruct'):
            sn = model.getStruct()
        reason = nc_method_refusal(sn, method, self.options)
        return (not reason), reason

    @staticmethod
    def supportsExactness(model, method):
        """(bool, reason) Product-form precondition of the normalizing-constant
        methods, the same rule runAnalyzer enforces at solve time. Only 'exact',
        'is' and 'panald' require it (the others fall back to Seidmann's
        comom on a non-product-form model), and 'is' on a pass-and-swap model is
        exempt (pfqn_pas_is). Product form has no registry feature name, so the
        check cannot live in getMethodFeatureSet. Mirrors MATLAB
        SolverNC.supportsModelMethod."""
        if str(method).lower() not in ('exact', 'is', 'panald'):
            return True, ''
        if model.hasProductFormSolution():
            return True, ''
        # A LOSS NETWORK (open, one DROP region holding a single Delay) IS product
        # form -- the truncated Poisson law the residue transform of
        # solver_nc_lossn_analyzer evaluates exactly under 'exact' -- but
        # sn_has_blocking reads any region as blocking, so hasProductFormSolution
        # says no. The shape is exempted here and at the 'exact' arm of
        # runAnalyzer alike, through the one predicate both ask.
        if _is_lossn_fcr_case(model.getStruct()):
            return True, ''
        if str(method).lower() == 'is':
            from .solver_nc_pas_is_analyzer import nc_is_pas_model
            from .solver_nc_oi_analyzer import nc_is_oi_model
            # 'is' is admissible on ANY OI/PAS model, not only a P&S one:
            # solver_nc_pas_is_analyzer samples the rank rate either way.
            if nc_is_pas_model(model.getStruct()) or nc_is_oi_model(model.getStruct()):
                return True, ''
        if str(method).lower() == 'exact':
            # An OI/PAS station carries a rank rate mu(n) of the whole occupancy
            # vector, so hasProductFormSolution reads false, but runAnalyzer routes
            # such a model to solver_nc_oi_analyzer, which is EXACT (pfqn_ncoi).
            # Without this exemption the guard refuses 'exact' before the OI routing
            # is ever reached -- the same exemption mvaDispatch makes via hasOIorPAS.
            from ...lang.base import SchedStrategy as _Sched
            _sn = model.getStruct()
            _sched = getattr(_sn, 'sched', {}) or {}
            _vals = _sched.values() if hasattr(_sched, 'values') else list(_sched)
            if any(v in (_Sched.OI, _Sched.PAS) for v in _vals):
                return True, ''
        return False, ("method '%s' requires a product-form solution; use "
                       "method 'comom' or SolverCTMC instead" % method)

    supports_exactness = supportsExactness

    def runAnalyzerChecks(self, options):
        """NC feature gate. Unlike the base gate this does not raise on an
        unrecognized method name: NC intentionally tolerates internal/auto names
        (e.g. 'default/comomld') and, as a SolverLN layer backend, runs with
        checks disabled. Only method-aware feature applicability is enforced."""
        if not getattr(self, 'enableChecks', True):
            return
        method = self.resolveMethod(options)
        ok, reason = self.supportsModelMethod(method)
        if ok:
            return
        # A STRUCTURAL refusal goes out in the ANALYZER's own words and with the
        # analyzer's own exception type. Those refusals used to surface from
        # inside runAnalyzer as a ValueError, and now that the gate sees them
        # first a caller catching ValueError from SolverNC(model, 'morrison')
        # must go on catching it. A FEATURE refusal keeps the LineError wrapper
        # it has always had -- supportsModelMethod reports the feature set first,
        # so a model refused for a feature never reaches the structural rules and
        # its message must not change type either.
        sn = getattr(self, '_sn', None)
        if sn is None and getattr(self, 'model', None) is not None:
            sn = self.model.getStruct()
        from .nc_method_refusal import nc_method_refusal
        if reason and reason == nc_method_refusal(sn, method, self.options):
            # A report-only refusal must not stop a run asked for by name: the
            # reference answers 'mmint2'/'gleint' with a zero table there.
            if not nc_method_refusal(sn, method, self.options, for_report=False):
                return
            raise ValueError(reason)
        raise LineError('This model contains features not supported by the solver. %s' % reason)

    resolve_method = resolveMethod
    supports_model_method = supportsModelMethod
    run_analyzer_checks = runAnalyzerChecks

    @staticmethod
    def getFeatureSet():
        """Get supported features as a SolverFeatureSet.

        NC supports limited features - notably not Cache with LRU replacement.
        """
        from ..base import SolverFeatureSet
        feat_supported = SolverFeatureSet()
        feat_supported.set_true([
            'Sink', 'Source',
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue',
            'APH', 'Coxian', 'Cox2', 'Erlang', 'Det', 'Exp', 'HyperExp',
            # Geometric is admitted for the discrete-time route only
            # (config['slotted'], solver_nc_dt_analyzer); on the
            # continuous-time routes it is treated by its mean and SCV
            'Geometric',
            'StatelessClassSwitcher', 'InfiniteServer',
            'SharedServer', 'Buffer', 'Dispatcher',
            'Server', 'JobSink', 'RandomSource', 'ServiceTunnel',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_SIRO',
            # DPS is served ONLY in Morrison's closed think+DPS shape
            # (nc_is_dps_model -> solver_nc_dps_analyzer). A boolean feature cannot
            # express that restriction, so runAnalyzer keeps an imperative check for
            # every other DPS model, the same pattern as 'Region'.
            'SchedStrategy_DPS',
            'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND', 'RoutingStrategy_SDR',
            'SchedStrategy_FCFS', 'SchedStrategy_OI', 'SchedStrategy_PAS',
            'ClosedClass', 'SelfLoopingClass',
            'Cache', 'CacheClassSwitcher', 'OpenClass', 'CacheRetrieval', 'CacheItemSize',
            # NC only supports RR and FIFO replacement, not LRU
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
            'ReplacementStrategy_HLRU',
            'LoadDependence',
            'ClassDependence',
            'JointDependence',
            # Fork-join through the MMT/HT transformation, driven by the
            # shared ForkJoinDriverMixin (as in SolverMVA)
            'Fork', 'Forker', 'Join', 'Joiner',
            'JoinPartial',  # quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
            # Petri nets: the 'rec' route (solver_nc_spn_analyzer) walks the
            # reachable set in a decision diagram, so a Place is a token
            # container rather than a station with a service process. A queueing
            # Place is refused by spn_pf, which is where the product-form class
            # is decided.
            'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting', 'Timing',
            'Firing', 'Storage',
            # c-server stations: every route but 'divdiff' carries the count
            # (see getMethodFeatureSet). FiniteCapacity is deliberately NOT
            # here: mem, default and exact are granted it per method, the rest
            # solve a buffer away.
            'MultiServer',
        ])
        return feat_supported

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported.

        Uses feature set checking to compare supported features against
        features used by the model. Prints warnings for unsupported features.
        """
        from ..base import SolverFeatureSet
        try:
            # Get used features from model
            if hasattr(model, 'get_used_lang_features'):
                feat_used = model.get_used_lang_features()
            elif hasattr(model, 'getUsedLangFeatures'):
                feat_used = model.getUsedLangFeatures()
            else:
                # Fall back to basic check
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

            # Compare features using SolverFeatureSet
            feat_supported = SolverNC.getFeatureSet()
            return SolverFeatureSet.supports(feat_supported, feat_used)
        except Exception:
            return False

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'default',
            'tol': 1e-4,
            'iter_max': 1000,
            'iter_tol': 1e-4,
            'verbose': default_verbose(),
        })

    # =========================================================================
    # Probability Methods
    # =========================================================================

    def getProb(self, station: Optional[int] = None) -> float:
        """The LOG probability of the declared state at a station.

        This is `@SolverNC/getProb.m`, and it returns a LOGARITHM even though
        its name does not say so: `solver_nc_marg`'s first output is `lPr` and
        both the MATLAB and the JAR reference hand it back unexponentiated
        (`Pnir = ret.lPr; return Pnir.get(ist)`). Its sibling `getProbSys`
        returns a plain probability, so the inconsistency is SolverNC's own and
        is reproduced rather than repaired in one codebase alone.

        It used to return `(1-rho) * rho^n` over a guessed range -- an M/M/1
        queue-length curve fitted to the mean utilization, which is neither this
        getter's quantity nor any NC quantity, and which agreed with the two
        references on nothing. `getProbAggr` is the aggregate station
        probability and `getProbMarg` the queue-length law; neither was affected.

        Args:
            station: Station index (0-based). With None, every station's value.

        Returns:
            The log probability, or the per-station vector of them.
        """
        from line_solver.api.solvers.nc import solver_nc_marg

        if self._result is None:
            self._ensureAvgResults()

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            p = prob_aggr_via_cpp(self)['prob']
            with np.errstate(divide='ignore'):
                lp = np.log(p)
            return float(lp[int(station)]) if station is not None else lp

        lPr = np.asarray(solver_nc_marg(self.model.getStruct(), self.options).lPr).flatten()
        if station is None:
            return lPr
        return float(lPr[int(station)])

    def getProbAggr(self, ist) -> float:
        """Get probability of a specific per-class job distribution at a station.

        Returns P(n1 jobs of class 1, n2 jobs of class 2, ...) for the state
        that was set via setState() on the station.

        Args:
            ist: Station index (0-based) or node object

        Returns:
            Probability that station ist is in the specified state.
        """
        from line_solver.api.solvers.nc import solver_nc_margaggr

        # Convert node object to index if needed (like MATLAB)
        if not isinstance(ist, (int, np.integer)):
            ist = ist.get_station_index0()

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            p = prob_aggr_via_cpp(self)['probAggr']
            if not (0 <= int(ist) < len(p)):
                raise ValueError("station index %r is outside 0..%d" % (ist, len(p) - 1))
            return float(p[int(ist)])

        if self._result is None:
            self._ensureAvgResults()

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn, RNG
        # seeding, and sample/seed forwarding") for the getStruct(true)-equivalent lag
        sn = self.model.get_struct()
        try:
            self.model._refresh_state()
            sn = self.model._sn
        except Exception:
            sn = getattr(self.model, '_sn', None) or self._sn
        if sn is None:
            raise RuntimeError("the model has no NetworkStruct to compute a state probability from")

        # Reuse the cached LOG normalizing constant (3rd output of
        # solver_nc_margaggr, not linear G) across repeated queries.
        options = self.options
        lG = getattr(self, '_logNormConstAggr', None)

        result = solver_nc_margaggr(sn, options, lG)

        if result.lG is not None and np.isfinite(result.lG):
            self._logNormConstAggr = result.lG

        # lPr holds log probabilities. An index past the end is a caller error,
        # not a zero-probability state: answering 0.0 made the two
        # indistinguishable. The result is the scalar MATLAB returns, not a row.
        if result.lPr is None:
            raise RuntimeError("solver_nc_margaggr returned no marginal probabilities")
        if not (0 <= int(ist) < len(result.lPr)):
            raise ValueError("station index %r is outside 0..%d" % (ist, len(result.lPr) - 1))
        return float(np.exp(np.asarray(result.lPr[int(ist)]).reshape(-1)[0]))

    def getProbMarg(self, station, jobclass=None) -> np.ndarray:
        """Get marginal queue-length distribution at station.

        Two routes, as in `@SolverNC/getProbMarg.m` and in the C++
        `solver_nc_getprob_marg`: under method='comom' the whole vector comes
        from one `pfqn_procomom` solve, and otherwise every total n is written
        as a sum of the AGGREGATE marginal over the per-class partitions of n.

        THE PROCOMOM ROUTE HAS NO ROW FOR A DELAY STATION -- it is solved on the
        Seidmann-reduced model, where the delays are folded into Z. Returning
        zeros there, which this getter used to do, is not a marginal: a delay in
        a closed network holds jobs with probability one at some n, and the law
        is the one the reference reaches by enumeration. So a delay station
        falls back to the enumeration rather than reporting an impossible
        station, exactly as the reference warns and falls back.

        Args:
            station: Station index (0-based) or station object
            jobclass: Job class index (unused here, kept for API compat)

        Returns:
            Marginal probability vector P(n_total = j) for j=0,1,...,sumN
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            return prob_via_jar(self, 'prob-marg', ist=station, jclass=(0 if jobclass is None else jobclass), kind='vector', raw_station=True)
        sn = self.model.getStruct()

        # Convert node object to station index if needed
        if hasattr(station, 'index'):
            ist = sn.nodeToStation[station.index]
        else:
            ist = int(station)

        if getattr(self.options, 'lang', 'python') == 'cpp':
            # `-a marg` under -s nc is the TOTAL queue-length law, which is what
            # this getter returns; jobclass is unused here for the same reason it
            # is unused natively.
            from ..cpp_dispatch import prob_marg_via_cpp
            return prob_marg_via_cpp(self, int(ist) + 1)[0]

        N = np.ravel(sn.njobs).astype(float)
        if not np.all(np.isfinite(N)):
            raise ValueError("getProbMarg is not implemented for models with open classes")
        Ntotal = int(round(float(np.sum(N))))

        method = str(getattr(self.options, 'method', 'default') or 'default').lower()
        if method == 'comom':
            from line_solver.api.sn.demands import sn_get_demands_chain
            from line_solver.api.pfqn.comom import pfqn_procomom

            demands = sn_get_demands_chain(sn)
            Lchain = demands.Lchain.copy()  # (M, C)
            Lchain[~np.isfinite(Lchain)] = 0.0
            Nchain = demands.Nchain.flatten().astype(int)  # (C,)

            M = sn.nstations
            C = sn.nchains

            # Matches solver_nc.m:97-119; multiserver: Seidmann splits demand
            # into per-server Lms = L/k and residual Zms = L*(k-1)/k
            queue_stations = []
            Z_total = np.zeros(C)
            Zms_total = np.zeros(C)
            Lms = np.zeros_like(Lchain)

            for i in range(M):
                if np.isinf(sn.nservers[i]):
                    # Delay station: accumulate into Z
                    Z_total += Lchain[i, :]
                else:
                    queue_stations.append(i)
                    k = sn.nservers[i]
                    Lms[i, :] = Lchain[i, :] / k
                    Zms_total += Lchain[i, :] * (k - 1) / k

            if ist in queue_stations:
                M_queues = len(queue_stations)
                L_queues = np.zeros((M_queues, C))
                for qi, orig_idx in enumerate(queue_stations):
                    L_queues[qi, :] = Lms[orig_idx, :]

                Pr, _ = pfqn_procomom(L_queues, Nchain, Z_total + Zms_total)
                # Pr is (M_queues x sumNchain+1), mapped onto 0..Ntotal
                out = np.zeros(Ntotal + 1)
                row = np.ravel(np.asarray(Pr, dtype=float)[queue_stations.index(ist), :])
                out[:min(len(row), Ntotal + 1)] = row[:Ntotal + 1]
                return out
            # a delay carries no procomom row: fall through to the enumeration
            line_warning("SolverNC", "comom does not directly support delay "
                                     "stations, using enumeration.")

        # The enumeration: P(n) = sum over the per-class partitions of n of the
        # AGGREGATE marginal, with lG computed once and reused across them.
        from line_solver.api.solvers.nc import solver_nc_margaggr_state

        K = int(sn.nclasses)
        Nmax = [int(round(float(N[r]))) for r in range(K)]
        lG = getattr(self, '_logNormConstAggr', None)

        Pmarg = np.zeros(Ntotal + 1)
        for n in range(Ntotal + 1):
            total = 0.0
            for part in _partitions(n, Nmax):
                logp, lG = solver_nc_margaggr_state(sn, ist, np.asarray(part, dtype=float),
                                                    self.options, lG)
                if np.isfinite(logp):
                    total += float(np.exp(logp))
            Pmarg[n] = total

        if lG is not None and np.isfinite(lG):
            self._logNormConstAggr = lG

        mass = float(np.sum(Pmarg))
        if mass > 0 and abs(mass - 1.0) > 1e-10:
            Pmarg = Pmarg / mass
        return Pmarg

    def getProbSys(self) -> float:
        """Get joint system state probability for the detailed state.

        For closed networks, this computes the probability of the current
        system state using the normalizing constant.

        Returns:
            float: Joint probability of the current system state.
        """
        # NC solver does not support detailed state probabilities,
        # delegate to aggregated version
        return self.getProbSysAggr()

    def getProbSysAggr(self) -> float:
        """Get aggregated system state probability.

        Computes the joint probability of observing the current queue length
        distribution across all stations using normalizing constants.

        Matches MATLAB: SolverNC.getProbSysAggr -> solver_nc_jointaggr

        Returns:
            float: Joint probability of the current system state.
        """
        from line_solver.api.solvers.nc import solver_nc_jointaggr

        if self._result is None:
            self._ensureAvgResults()

        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import prob_aggr_via_cpp
            return float(prob_aggr_via_cpp(self)['probSysAggr'])

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn")
        self.model.reset_struct()
        sn = self.model.getStruct()
        if sn is None:
            raise RuntimeError("the model has no NetworkStruct to compute a state probability from")

        result = solver_nc_jointaggr(sn, self.options)

        # A missing Pr is a failed computation, not a zero-probability state.
        if not hasattr(result, 'Pr') or result.Pr is None:
            raise RuntimeError("solver_nc_jointaggr returned no joint probability")
        return float(result.Pr)

    def getProbSysMarg(self, nvec, engine: str = 'exact'):
        """Joint probability of the per-station TOTAL queue lengths.

        Returns P(n_1 = nvec[0], ..., n_M = nvec[M-1]), all classes summed out.

        Compare with getProbSysAggr, which fixes the PER-CLASS population of
        every station and is a product form; each value returned here is the
        sum of getProbSysAggr over every per-class table with these row sums.
        Compare also with getProbMarg, which is the one-station marginal of
        this law.

        The quantity is a matrix permanent of the demand matrix replicated once
        per job (Ryser 1963 for the evaluation), so it needs no enumeration of
        that fibre.

        Args:
            nvec: (nstations,) per-station total job counts; must sum to the
                total closed population
            engine: permanent engine, one of 'exact' (default), 'bethe',
                'spm', 'heur', 'huberlaw', 'adapart'. Only 'exact' is exact; the
                others are refused on a demand matrix with a structural zero
                rather than having it floored, since they need full support.

        Returns:
            (Pn, lPn): the joint probability and its logarithm, which survives
            populations Pn underflows at
        """
        from line_solver.api.solvers.nc import solver_nc_jointmarg

        if getattr(self.options, 'lang', 'python') == 'cpp':
            raise NotImplementedError(
                "getProbSysMarg is not available under lang='cpp'; use lang='python'")

        # see _kb/06-solver-catalog.md ("Python NC: shared-reference sn")
        self.model.reset_struct()
        sn = self.model.getStruct()
        if sn is None:
            raise RuntimeError("the model has no NetworkStruct to compute a state probability from")

        # Reuse the constant when a previous getter already paid for it:
        # sweeping the whole lattice of total states otherwise recomputes G
        # once per state.
        lG = getattr(self, '_logNormConstAggr', None)

        result = solver_nc_jointmarg(sn, self.options, nvec, engine, lG)

        if result.lG is not None and np.isfinite(result.lG):
            self._logNormConstAggr = result.lG
        self.lastPermEngine = str(engine).lower() if engine else 'exact'
        return float(result.Pr), float(result.lPr)

    # =========================================================================
    # Unified Metrics Method
    # =========================================================================


    # =========================================================================
    # Chain-Level Methods
    # =========================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self._sn, 'chains') and self._sn.chains is not None:
            chains = []
            nchains = getattr(self._sn, 'nchains', 1)
            for c in range(nchains):
                chain_classes = []
                for k in range(self._sn.nclasses):
                    if hasattr(self._sn.chains, '__getitem__'):
                        if self._sn.chains[c, k] > 0:
                            chain_classes.append(k)
                chains.append(chain_classes)
            return chains if chains else [[k for k in range(self._sn.nclasses)]]
        return [[k] for k in range(self._sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()
        Q = self._result.Q
        chains = self._get_chains()
        QN = np.zeros((Q.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                QN[:, c] = np.sum(Q[:, cc], axis=1)
        return QN

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()
        self._cap_unstable_open_util()
        U = self._result.U
        chains = self._get_chains()
        UN = np.zeros((U.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                UN[:, c] = np.sum(U[:, cc], axis=1)
        return UN

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()
        R = self._result.R
        chains = self._get_chains()
        RN = np.zeros((R.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                RN[:, c] = np.mean(R[:, cc], axis=1)
        return RN

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self._result is None:
            self._ensureAvgResults()
        T = self._result.T
        chains = self._get_chains()
        TN = np.zeros((T.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                TN[:, c] = np.sum(T[:, cc], axis=1)
        return TN

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics aggregated by chain."""
        return (self.getAvgQLenChain(), self.getAvgUtilChain(), self.getAvgRespTChain(),
                self.getAvgResidTChain(), self.getAvgArvRChain(), self.getAvgTputChain())

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
                    'QLen': QN[i,c], 'Util': UN[i,c], 'RespT': RN[i,c],
                    'ResidT': WN[i,c], 'ArvR': AN[i,c], 'Tput': TN[i,c]
                })
        # five SIGNIFICANT digits like MATLAB's table, not pandas' five decimals
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(pd.DataFrame(rows))

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Get average metrics per node.

        Unlike getAvg() which returns station-level metrics, this method
        returns node-level metrics including non-station nodes (e.g., Router/VSink).

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        from ...api.sn.getters import sn_get_node_arvr_from_tput, sn_get_node_tput_from_tput

        if self._result is None:
            self._ensureAvgResults()

        QN = self._result.Q
        UN = self._result.U
        RN = self._result.R
        TN = self._result.T
        AN = sn_get_arvr_from_tput(self._sn, TN)  # Compute proper arrival rates from routing

        sn = self._sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses

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

        # Fix arrival rates for ClassSwitch and Sink nodes for cache hit/miss classes
        # (MATLAB: getAvgNode.m lines 54-76)
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
        non-station nodes like Router/VSink.

        Returns:
            pandas.DataFrame with columns: Node, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()

        sn = self._sn
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        class_names = list(sn.classnames) if hasattr(sn, 'classnames') and sn.classnames else []

        from ..cache_table import retrieval_hidden_classes
        hidden = retrieval_hidden_classes(sn)

        rows = []
        for node_idx in range(sn.nnodes):
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'

            for r in range(sn.nclasses):
                if r in hidden:
                    continue  # auxiliary retrieval class - omit from node table
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
        """Get transient average metrics (not supported for NC).

        NC is a steady-state solver. Returns steady-state values.

        Returns:
            Tuple of (Q, U, T) steady-state values
        """
        if self._result is None:
            self._ensureAvgResults()
        return self._result.Q, self._result.U, self._result.T

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics."""
        return self.getAvgSysRespT(), self.getAvgSysTput()

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame.

        Returns chain-level metrics matching MATLAB/Java implementation.
        The table includes:
        - Chain: Chain name (Chain1, Chain2, ...)
        - JobClasses: Class names within each chain
        - SysRespT: Chain response time
        - SysTput: Chain throughput
        """
        CN, XN = self.getAvgSys()

        sn = self._sn
        if sn is None:
            sn = self.model.getStruct(True)
            self._sn = sn

        nchains = sn.nchains if hasattr(sn, 'nchains') and sn.nchains > 0 else 1
        CN = np.atleast_1d(np.asarray(CN, dtype=float)).flatten()
        XN = np.atleast_1d(np.asarray(XN, dtype=float)).flatten()
        CN = [CN[c] if c < len(CN) else 0.0 for c in range(nchains)]
        XN = [XN[c] if c < len(XN) else 0.0 for c in range(nchains)]
        return self._make_sys_table(CN, XN)

    # =========================================================================
    # Matrix-Exponential Method for Open Networks
    # =========================================================================

    def me_open(self, options: Optional[Any] = None) -> Dict[str, Any]:
        """Maximum Entropy Method (MEM) for open queueing networks.

        Implements the Kouvatsos (1994) entropy-maximisation algorithm
        (mirrors MATLAB ``@SolverNC/me_open.m`` and Java ``SolverNC.meOpen``).
        Supports only open networks (no closed classes).

        Args:
            options: optional solver options; MEM tolerance/iteration limits
                are read from ``options.config`` (``mem_tol``, ``mem_maxiter``,
                ``mem_verbose``).

        Returns:
            dict with ``QN``, ``UN``, ``RN``, ``TN`` (M x R arrays for queue
            lengths, utilizations, response times, throughputs), ``CN``/``XN``
            (1 x R system response times and throughputs) and ``method='mem'``.

        Reference:
            D.D. Kouvatsos, "Entropy Maximisation and Queueing Network Models",
            Annals of Operations Research, 48:63-126, 1994.
        """
        from ...api.me import solver_nc_mem

        if options is None:
            options = self.options
        sn = self._sn
        if sn is None and hasattr(self, 'model'):
            sn = self.model.get_struct()

        QN, UN, RN, TN, CN, XN, _ = solver_nc_mem(sn, options)

        return {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
            'CN': CN,
            'XN': XN,
            'method': 'mem',
        }

    # =========================================================================
    # Sampling Methods (Not Supported - Analytical Solver)
    # =========================================================================

    def sample(self, node: int, numEvents: int) -> np.ndarray:
        """Sampling not supported by NC (analytical solver)."""
        raise NotImplementedError(
            "Sampling not supported by SolverNC. "
            "Use SolverSSA or SolverLDES for simulation-based analysis."
        )

    GetProb = getProb
    GetProbAggr = getProbAggr
    GetProbMarg = getProbMarg
    GetProbSys = getProbSys
    GetProbSysAggr = getProbSysAggr
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
    GetAvgChain = getAvgChain
    GetAvgChainTable = getAvgChainTable
    GetAvgNode = getAvgNode
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetAvgSys = getAvgSys
    GetAvgSysTable = getAvgSysTable
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain
    GetNormalizingConstant = getNormalizingConstant
    GetLogNormalizingConstant = getLogNormalizingConstant
    GetProbNormConstAggr = getProbNormConstAggr
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    MeOpen = me_open
    ListValidMethods = listValidMethods
    GetFeatureSet = getFeatureSet
    Supports = supports
    DefaultOptions = defaultOptions
    default_options = defaultOptions
    GetTranAvg = getTranAvg

    # Node-chain specific aliases
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain

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
    avg_sys_table = getAvgSysTable
    avg_node_table = getAvgNodeTable
    avg_chain_table = getAvgChainTable
    avg_node_chain_table = getAvgNodeChainTable
    run_analyzer = runAnalyzer
    get_normalizing_constant = getNormalizingConstant
    get_log_normalizing_constant = getLogNormalizingConstant
    prob = getProb
    prob_aggr = getProbAggr
    prob_marg = getProbMarg
    prob_sys = getProbSys
    prob_sys_aggr = getProbSysAggr
    prob_norm_const_aggr = getProbNormConstAggr


__all__ = ['SolverNC', 'SolverNCOptions']

