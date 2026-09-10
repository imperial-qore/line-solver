"""
Native Python implementation of automatic solver selection.

This module provides the SolverAUTO class that analyzes the network model
and automatically selects the most appropriate solver based on model characteristics.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Type
from dataclasses import dataclass, field
from contextlib import contextmanager
from ...constants import default_verbose, ProcessType
from enum import Enum

from ...api.sn import NetworkStruct, SchedStrategy, NodeType
from ...api.io.logging import line_debug
from ..base import NetworkSolver

# Population at or below which an exact solver is preferred over an
# approximation, when one is available.
EXACT_POPULATION_MAX = 5


class SolverType(Enum):
    """Available solver types."""
    MVA = 'MVA'
    LDES = 'LDES'
    NC = 'NC'
    CTMC = 'CTMC'
    SSA = 'SSA'
    FLUID = 'Fluid'
    MAM = 'MAM'
    JMT = 'JMT'


@dataclass
class ModelAnalyzer:
    """
    Analyzes queueing network model characteristics for solver selection.

    This class examines the network structure and determines which solvers
    are applicable and which is likely to be most efficient.

    Original implementation by LINE
    """
    sn: NetworkStruct

    def has_single_chain(self) -> bool:
        """Check if model has exactly one chain."""
        return self.sn.nchains == 1

    def has_multi_chain(self) -> bool:
        """Check if model has multiple chains."""
        return self.sn.nchains > 1

    def is_closed_model(self) -> bool:
        """Check if all classes are closed (finite population)."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return False
        return np.all(np.isfinite(self.sn.njobs.flatten()))

    def is_open_model(self) -> bool:
        """Check if all classes are open (infinite population)."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return True
        return np.all(np.isinf(self.sn.njobs.flatten()))

    def is_mixed_model(self) -> bool:
        """Check if model has both open and closed classes."""
        if self.sn.njobs is None or len(self.sn.njobs) == 0:
            return False
        njobs_flat = self.sn.njobs.flatten()
        has_open = np.any(np.isinf(njobs_flat))
        has_closed = np.any(np.isfinite(njobs_flat) & (njobs_flat > 0))
        return has_open and has_closed

    def has_product_form(self) -> bool:
        """
        Check if model has product form solution.

        Matches MATLAB's sn_has_product_form.m:
        - All stations use INF, PS, FCFS, LCFS-PR, LCFS, or EXT scheduling
        - No multiclass heterogeneous FCFS
        - No priorities
        - No fork-join
        - No state-dependent routing
        """
        from ...api.sn import sn_has_product_form
        return sn_has_product_form(self.sn)

    def has_multi_server(self) -> bool:
        """Check if any station has multiple servers."""
        return self.sn.has_multi_server()

    def has_load_dependence(self) -> bool:
        """Check if model has load-dependent service rates."""
        return self.sn.has_load_dependence()

    def has_homogeneous_scheduling(self, strategy: SchedStrategy) -> bool:
        """Check if all stations use the same scheduling strategy."""
        for station_id, sched in self.sn.sched.items():
            if sched != strategy:
                return False
        return True

    def has_all_infinite_servers(self) -> bool:
        """Check if all stations are infinite servers (delays)."""
        return self.has_homogeneous_scheduling(SchedStrategy.INF)

    def get_total_jobs(self) -> float:
        """Get total number of jobs in closed classes."""
        return self.sn.get_total_population()

    def get_avg_jobs_per_chain(self) -> float:
        """Get average number of jobs per chain."""
        if self.sn.nchains == 0:
            return 0.0
        return self.get_total_jobs() / self.sn.nchains

    def has_fork_join(self) -> bool:
        """Check if model has fork-join nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.FORK or nt == NodeType.JOIN:
                return True
        return False

    def has_class_switching(self) -> bool:
        """Check if model has class switching nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.CLASSSWITCH:
                return True
        return False

    def has_cache(self) -> bool:
        """Check if model has cache nodes."""
        if self.sn.nodetype is None:
            return False
        for nt in self.sn.nodetype:
            if nt == NodeType.CACHE:
                return True
        return False

    def has_finite_capacity(self) -> bool:
        """Check if model has finite capacity constraints."""
        if self.sn.cap is not None and len(self.sn.cap) > 0:
            return np.any(np.isfinite(self.sn.cap) & (self.sn.cap > 0) & (self.sn.cap < 1e6))
        return False

    def get_num_stations(self) -> int:
        """Get number of service stations."""
        return self.sn.nstations

    def get_num_classes(self) -> int:
        """Get number of job classes."""
        return self.sn.nclasses

    def get_state_space_size_estimate(self) -> float:
        """
        Estimate the size of the state space.

        Uses the formula: product of (n_i + K - 1) choose (K - 1) over all classes
        where n_i is the population of class i and K is number of stations.
        """
        if not self.is_closed_model():
            return float('inf')

        M = self.sn.nstations
        N = self.get_total_jobs()

        if N == 0 or M == 0:
            return 1.0

        # Simple approximation: (N + M - 1)! / (N! * (M-1)!)
        # Using Stirling's approximation for large values
        if N > 100 or M > 20:
            # Very rough upper bound
            return (N + M) ** min(N, M)

        # For small values, compute more accurately
        from math import comb
        return comb(int(N) + M - 1, M - 1)


@dataclass
class SolverAUTOOptions:
    """Options for the automatic solver.

    Attributes:
        selection_method: Solver selection strategy:
            - 'default', 'heur': Heuristic-based automatic selection
            - 'sim': Best simulator (LDES/SSA)
            - 'exact': Best exact method (NC for closed, CTMC otherwise)
            - 'fast': Fast approximate method (MVA)
            - 'accurate': Accurate approximate method (Fluid)
            # - 'ai': AI-based selection (not yet available)
        forced_solver: Force a specific solver by name
        verbose: Enable verbose output
    """
    selection_method: str = 'default'
    forced_solver: Optional[str] = None
    verbose: bool = field(default_factory=default_verbose)


class SolverAUTO(NetworkSolver):
    """
    Native Python automatic solver selection.

    This solver analyzes the network model and automatically selects the most
    appropriate solver based on model characteristics, feature support, and
    performance considerations.

    Selection Algorithm

    For closed networks:
    - Single chain: Use NC (fast exact solution)
    - Multi-chain, all infinite servers: Use MVA
    - Multi-chain, product form, no multi-server: Use NC
    - FCFS without multi-server:
        - avg_jobs_per_chain > 30: Use Fluid
        - 10 < avg_jobs_per_chain <= 30: Use MVA
        - total_jobs < 5: Use NC
    - PS/PSPRIO with multi-server: Use MVA
    - FCFS with multi-server:
        - total_jobs < 5: Use NC
        - total_jobs >= 5: Use MVA

    For open networks:
    - Simple open: Use MVA or LDES

    Args:
        model: Network model (Python wrapper or native structure)
        options: SolverAUTOOptions configuration
        **kwargs: Additional solver options passed to selected solver

    Example:
        >>> solver = SolverAUTO(model)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()
        >>> print(f"Selected: {solver.get_selected_solver_name()}")
    """

    def __init__(
        self,
        model: Any,
        options: Optional[SolverAUTOOptions] = None,
        **kwargs
    ):
        self.model = model
        # A selection token given positionally ('bound', 'fast', ...) must reach
        # selection_method, else it silently runs the default heuristic.
        if isinstance(options, str):
            options = SolverAUTOOptions(selection_method=options.lower())
            if 'verbose' in kwargs:
                options.verbose = kwargs['verbose']
        self.options = options or SolverAUTOOptions()
        # The SAME token given by keyword, which is how the CLI passes -s: it
        # calls LINE(model, method=solver_str). Without this it landed in
        # kwargs, runAnalyzer resolved the untouched default, and every named
        # solver ran the heuristic instead -- on an LQN that meant `-s ln`
        # answering with the external lqns binary, the first candidate of the
        # layered pool.
        token = kwargs.get('method', kwargs.get('selection_method'))
        if token is not None:
            self.options.selection_method = str(token).lower()
        self.kwargs = kwargs
        # Method pinned per solver by _choose_ranked (see _create_solver).
        self._pinned_methods = {}

        self._result = None
        self._selected_solver = None
        self._selected_solver_name = None
        self._candidate_solvers: List[str] = []
        # Family pinned by a 'family' or 'family.submethod' method name, or by the
        # model class itself (LayeredNetwork -> ln, Environment -> env, a model
        # carrying Prior parameters -> uq). None means the heuristic decides.
        self._pinned_family: Optional[str] = None
        self._pinned_family_method: str = 'default'

        # A LayeredNetwork or an Environment has NO queueing NetworkStruct of
        # its own, so the analyzer and the candidate pool below do not apply to
        # it; its family is fixed by its class, exactly as MATLAB SolverAUTO
        # switches on class(model).
        self._model_family = self._detect_model_family()
        if self._model_family in ('layered', 'environment'):
            self._sn = None
            self._analyzer = None
            self._pinned_family = 'ln' if self._model_family == 'layered' else 'env'
            self._candidate_solvers = ['LN'] if self._model_family == 'layered' else ['ENV']
            return

        # Extract network structure
        self._sn = self._get_network_struct()

        # Create model analyzer
        self._analyzer = ModelAnalyzer(sn=self._sn)

        # Determine candidate solvers
        self._determine_candidates()

        # Uncertain parameters expand the model into one instance per
        # alternative, so the intent applies to each of them rather than to one
        # model. Mirrors MATLAB SolverAUTO.hasPriorParameters.
        if self._has_prior_parameters():
            self._pinned_family = 'uq'
            self._candidate_solvers = ['UQ']

    def _detect_model_family(self) -> str:
        """'layered', 'environment' or 'network', from the model's class."""
        cls_names = {c.__name__ for c in type(self.model).__mro__}
        if 'LayeredNetwork' in cls_names:
            return 'layered'
        if 'Environment' in cls_names:
            return 'environment'
        return 'network'

    def _has_prior_parameters(self) -> bool:
        """True when the model carries Prior (uncertain) service or arrival laws."""
        from ...distributions.continuous import Prior
        model = self.model
        if not hasattr(model, 'getNodes') and not hasattr(model, 'get_nodes'):
            return False
        try:
            nodes = model.getNodes() if hasattr(model, 'getNodes') else model.get_nodes()
            classes = (model.getClasses() if hasattr(model, 'getClasses')
                       else model.get_classes())
        except Exception:
            return False
        for node in nodes:
            node_cls = {c.__name__ for c in type(node).__mro__}
            if node_cls & {'Queue', 'Delay'}:
                for jc in classes:
                    try:
                        dist = node.getService(jc)
                    except Exception:
                        dist = None
                    if isinstance(dist, Prior):
                        return True
            elif 'Source' in node_cls:
                for jc in classes:
                    try:
                        dist = node.getArrival(jc)
                    except Exception:
                        dist = None
                    if isinstance(dist, Prior):
                        return True
        return False

    # ------------------------------------------------------------------
    # Method-token resolution (MATLAB @SolverAUTO/resolveMethodToken.m)
    # ------------------------------------------------------------------

    @staticmethod
    def selectionIntents() -> List[str]:
        """Tokens that state what the solution is for, rather than naming an
        algorithm."""
        return ['default', 'heur', 'sim', 'exact', 'fast', 'accurate', 'bound']

    selection_intents = selectionIntents

    @staticmethod
    def familyNames() -> List[str]:
        """Method families, in the order in which an unqualified method name is
        looked up.

        'ag' SITS AFTER 'mam', whose RCAT names it took over, and is a family for
        the METHOD NAME and REPORT tables only: LINE(model, 'ag.inap') and
        model.help() reach SolverAG through it. It is deliberately NOT a
        candidate of the automatic ranking (the constructor builds that list by
        hand) and NOT in the feature-set union of getFeatureSet, so a G-network
        is still refused by 'default' and has to be asked for by name; see
        _kb/06-solver-catalog.md, "SolverAG owns the RCAT methods".
        """
        return ['mva', 'nc', 'ctmc', 'fluid', 'mam', 'ag', 'ba', 'ssa', 'ldes', 'jmt',
                'qns', 'ln', 'env', 'lqns', 'uq']

    family_names = familyNames

    @staticmethod
    def familyAcceptsModelClass(family: str, model: Any) -> bool:
        """Is `family` defined over models of `model`'s kind at all?

        The composite families are not universal: SolverLN and SolverLQNS
        analyze a LayeredNetwork and SolverENV an Environment, which is the
        partition the CONSTRUCTOR already uses when it builds the candidate
        list. Their feature sets describe what they accept inside a LAYER, so
        asking one whether it supports a flat Network gets a yes to a question
        it was never asked: SolverLQNS.getFeatureSet names Queue, Exp and
        SchedStrategy_FCFS, which every closed exponential network uses, and it
        was offered for one. Mirrors MATLAB SolverAUTO.familyAcceptsModelClass.

        THE PARTITION HAS TO BE READ IN BOTH DIRECTIONS, and reading it in one
        was a bug. Naming the model class each composite family takes says
        nothing about what the FLAT-NETWORK families take, so every one of them
        answered yes for an Environment and the report offered seven SolverQNS
        methods on the tut12 random-environment model, half the table. SolverQNS
        does not refuse such a model either: its constructor leaves self.model as
        None rather than raising, so the gate is asked about nothing at all and
        says yes, and the refusal arrives only at getAvgTable as a bare
        AttributeError: nstations.

        A LayeredNetwork needs no such rule and deliberately does not get one:
        SolverLDES analyzes an LQN natively and belongs in that report, and the
        families that do not take one fail to CONSTRUCT over it, which findSolver
        already treats as "contributes nothing". An Environment is the class
        where construction does not filter.
        """
        cls = type(model).__name__
        if family in ('ln', 'lqns'):
            return cls == 'LayeredNetwork'
        if family == 'env':
            return cls in ('Environment', 'Env')
        return cls not in ('Environment', 'Env')

    family_accepts_model_class = familyAcceptsModelClass

    @staticmethod
    def methodAliasPrefixes(family: str):
        """The prefixes under which `family` advertises a SECOND SPELLING of a
        method it already declares plainly.

        THIS IS A DECLARATION, not a derivation, and belongs with familyMetrics
        and methodClass for the same reason: the knowledge lives in the solver's
        own dispatch (SolverMVA strips a leading 'amva.' before selecting an
        algorithm) and there is no accessor that exposes it, so a family that
        gains or loses an alias spelling must be edited into all four copies in
        the SAME change. An omission does not fail; it puts the same algorithm in
        the report twice.

        Mirrors MATLAB SolverAUTO.methodAliasPrefixes.
        """
        return ('amva.',) if family == 'mva' else ()

    method_alias_prefixes = methodAliasPrefixes

    @staticmethod
    def isMethodAlias(family: str, name: str, declared) -> bool:
        """Is `name` a second spelling of another method this family declares?

        The remainder has to be declared too, which is what keeps the rule from
        eating a genuine method that merely starts with the prefix: it is an
        alias only when the thing it aliases is there beside it.
        """
        for prefix in SolverAUTO.methodAliasPrefixes(family):
            if name.startswith(prefix) and name[len(prefix):] in declared:
                return True
        return False

    is_method_alias = isMethodAlias

    @staticmethod
    def familyAlias(name: Optional[str]) -> str:
        """Canonical family of a method name, or '' when the method name names none."""
        if not name or not isinstance(name, str):
            return ''
        return {
            'mam': 'mam', 'ag': 'ag', 'mva': 'mva', 'nc': 'nc',
            'fluid': 'fluid', 'fld': 'fluid',
            'jmt': 'jmt', 'ssa': 'ssa', 'ctmc': 'ctmc',
            'ldes': 'ldes', 'des': 'ldes',
            'ba': 'ba', 'env': 'env', 'ln': 'ln',
            'lqns': 'lqns', 'lqsim': 'lqns',
            'qns': 'qns', 'uq': 'uq',
        }.get(name.lower(), '')

    family_alias = familyAlias

    def _family_declaring_method(self, token: str) -> str:
        """Family that declares an unqualified algorithm name, by asking each
        family for its own method list.

        Keeping the question there avoids a second copy of the name table in
        SolverAUTO, which would drift. Mirrors MATLAB familyDeclaringMethod.m.
        """
        for family in self.familyNames():
            try:
                probe = self._build_family_solver(family, 'default', probe=True)
                declared = probe.listValidMethods()
            except Exception:
                # A family that cannot even be instantiated on this model cannot
                # own the token; the next one is asked instead.
                continue
            if declared and any(token.lower() == str(d).lower() for d in declared):
                return family
        return ''

    def resolveMethodToken(self, token: Optional[str]):
        """Split a method name into a selection intent, or into a family and
        the submethod handed to it.

        A qualified method name 'family.submethod' keeps its submethod: dropping it
        would silently downgrade a pinned method to the family default. An
        unrecognised token is an ERROR, not a fall-through to the default
        heuristic -- a caller who names a real algorithm must not be answered by
        a different one. Mirrors MATLAB @SolverAUTO/resolveMethodToken.m.

        Args:
            token: the selection_method string

        Returns:
            (kind, family, submethod) with kind either 'intent' or 'family'.

        Raises:
            ValueError: when the method name is neither an intent, nor a family, nor
                an algorithm any family declares.
        """
        if not token:
            token = 'default'
        token = str(token)
        if token.lower() == 'auto':
            token = 'default'

        intents = self.selectionIntents()
        if token.lower() in intents:
            return 'intent', token.lower(), 'default'

        head, _, rest = token.partition('.')
        fam = self.familyAlias(head)
        if fam:
            if not rest:
                # A bare family name means its default method; for bounds the
                # default is the composite tightest-of-all family rather than a
                # single one.
                return 'family', fam, ('auto' if fam == 'ba' else 'default')
            return 'family', fam, rest

        # Unqualified algorithm name, e.g. 'comom'. The family that declares it
        # owns it, so the name table stays in the families themselves.
        fam = self._family_declaring_method(token)
        if fam:
            return 'family', fam, token

        raise ValueError(
            "SolverAUTO: unrecognized method '%s'. Valid tokens are a selection "
            "intent (%s), a method family (%s), or a qualified method name such "
            "as 'nc.comom'." % (token, ', '.join(intents),
                                ', '.join(self.familyNames())))

    resolve_method_token = resolveMethodToken

    def _inner_factory(self, method: str):
        """Layer/stage solver factory of a composite family.

        A composite family may be qualified by the family that solves its inner
        models, as in 'ln.mva' or 'env.fluid'. Anything else is an intent, and
        the inner models are then solved by SolverAUTO itself, which selects per
        submodel. Mirrors the innerFactory helper of buildFamilySolver.m.
        """
        fam = self.familyAlias(method)
        if fam and fam not in ('ln', 'env'):
            from ..solver_mva import SolverMVA
            from ..solver_nc import SolverNC
            from ..solver_mam import SolverMAM
            from ..solver_fld import SolverFLD
            from ..solver_ctmc import SolverCTMC
            from ..solver_ssa import SolverSSA
            ctor = {'mva': SolverMVA, 'nc': SolverNC, 'mam': SolverMAM,
                    'fluid': SolverFLD, 'ctmc': SolverCTMC, 'ssa': SolverSSA}.get(fam)
            if ctor is not None:
                return lambda m, _c=ctor: _c(m, verbose=False)
            if fam == 'ldes':
                from ..wrappers.solver_ldes import SolverLDES
                return lambda m: SolverLDES(m)
            if fam == 'jmt':
                from ..wrappers.solver_jmt import SolverJMT
                return lambda m: SolverJMT(m, verbose=False)
        return lambda m: SolverAUTO(m, verbose=False)

    def _build_family_solver(self, family: str, method: str = 'default',
                             probe: bool = False) -> Any:
        """Instantiate the solver of a method family.

        Args:
            family: canonical family name, as familyAlias returns
            method: submethod already resolved by resolveMethodToken, so a
                pinned method reaches the family that runs it
            probe: build only to read listValidMethods, so the construction must
                stay cheap and must not run anything

        Mirrors MATLAB @SolverAUTO/buildFamilySolver.m.
        """
        kw = dict(self.kwargs)
        kw.pop('method', None)
        if probe:
            kw['verbose'] = False
        mkw = dict(kw)
        if method and method != 'default':
            mkw['method'] = method

        if family == 'mva':
            from ..solver_mva import SolverMVA
            return SolverMVA(self.model, **mkw)
        if family == 'nc':
            from ..solver_nc import SolverNC
            return SolverNC(self.model, **mkw)
        if family == 'mam':
            from ..solver_mam import SolverMAM
            return SolverMAM(self.model, **mkw)
        if family == 'ag':
            from ..solver_ag import SolverAG
            return SolverAG(self.model, method or 'default', **kw)
        if family == 'fluid':
            from ..solver_fld import SolverFLD
            return SolverFLD(self.model, **mkw)
        if family == 'ctmc':
            from ..solver_ctmc import SolverCTMC
            return SolverCTMC(self.model, **mkw)
        if family == 'ssa':
            from ..solver_ssa import SolverSSA
            return SolverSSA(self.model, **mkw)
        if family == 'ldes':
            from ..wrappers.solver_ldes import SolverLDES, LDESOptions
            opts = LDESOptions(**{k: v for k, v in kw.items()
                                  if hasattr(LDESOptions, k)})
            return SolverLDES(self.model, opts)
        if family == 'jmt':
            from ..wrappers.solver_jmt import SolverJMT
            return SolverJMT(self.model, **mkw)
        if family == 'ba':
            from ..solver_ba import SolverBA
            return SolverBA(self.model, method or 'auto', **kw)
        if family == 'qns':
            from ..wrappers.solver_qns import SolverQNS
            return SolverQNS(self.model, **kw)
        if family == 'lqns':
            from ..wrappers.solver_lqns import SolverLQNS
            return SolverLQNS(self.model, **kw)
        if family == 'ln':
            from ..solver_ln import SolverLN
            return SolverLN(self.model, self._inner_factory(method))
        if family == 'env':
            from ...environment import SolverENV
            return SolverENV(self.model, self._inner_factory(method))
        if family == 'uq':
            from ..solver_uq import SolverUQ
            return SolverUQ(self.model, solver_factory=self._inner_factory(method))
        raise ValueError("SolverAUTO: unknown method family '%s'." % family)

    def _get_network_struct(self) -> NetworkStruct:
        """Get NetworkStruct from model."""
        model = self.model

        # If already a NetworkStruct
        if isinstance(model, NetworkStruct):
            return model

        # Native model with cached _sn
        if hasattr(model, '_sn') and model._sn is not None:
            return model._sn

        # Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                return model._sn

        # Native model with snake-case get_struct() (no wrapper bridge — native
        # solvers reject JAR-wrapper models, keeping python/ free of any
        # JAR/JVM coupling).
        if hasattr(model, 'get_struct'):
            sn = model.get_struct()
            if sn is not None:
                return sn

        # Already a native NetworkStruct (duck-typed)
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError(
            "Cannot extract a native NetworkStruct from model. Native solvers "
            "accept only native Network / NetworkStruct inputs (no JAR wrapper).")

    def _determine_candidates(self) -> None:
        """Determine which solvers can handle this model.

        The pool has to span every solver a ranking may reach: a solver absent
        here can never be chosen, however well it fits the metric.
        """
        self._candidate_solvers = []

        # MVA can handle most closed and some open networks
        if self._analyzer.is_closed_model() or self._analyzer.is_open_model():
            if not self._analyzer.has_fork_join():
                self._candidate_solvers.append('MVA')

        # LDES can handle most networks
        if not self._analyzer.has_cache():
            self._candidate_solvers.append('LDES')

        # NC for product-form networks
        if self._analyzer.has_product_form() and self._analyzer.is_closed_model():
            self._candidate_solvers.append('NC')

        # Cache models: NC leads the analytical order and MVA still applies
        if self._analyzer.has_cache():
            for name in ('NC', 'MVA'):
                if name not in self._candidate_solvers:
                    self._candidate_solvers.append(name)

        # Fluid, MAM, CTMC and SSA complete the pool the rankings index. JMT is
        # excluded: LDES subsumes its feature set, so automatic selection never
        # dispatches to the external simulator.
        for name in ('Fluid', 'MAM', 'CTMC', 'SSA'):
            self._candidate_solvers.append(name)

        # Default fallback
        if not self._candidate_solvers:
            self._candidate_solvers.append('MVA')

    def _traits(self) -> dict:
        """Structural traits that drive the ranking, mirroring MATLAB solverTraits.m."""
        sn = self._sn
        traits = {'cache': self._analyzer.has_cache(), 'fcr': False, 'map': False,
                  'prio': 'none'}
        traits['fcr'] = bool(getattr(sn, 'nregions', 0))
        procid = getattr(sn, 'procid', None)
        if procid is not None:
            try:
                # MAP and MMPP2 carry autocorrelation only MAM reproduces.
                # procid holds enums or their values, so compare by name.
                names = [str(getattr(t, 'name', t)).upper()
                         for t in np.asarray(procid, dtype=object).flatten()]
                traits['map'] = any(n in ('MAP', 'MMPP2') for n in names)
            except Exception:
                traits['map'] = False
        sched = getattr(sn, 'sched', None)
        if sched is not None:
            names = [str(getattr(s, 'name', s)).upper() for s in np.asarray(sched).flatten()]
            if any(n in ('FCFSPRPRIO', 'FCFSPIPRIO', 'LCFSPRPRIO', 'LCFSPIPRIO') for n in names):
                traits['prio'] = 'preempt'
            elif any(n in ('PSPRIO', 'DPSPRIO', 'GPSPRIO') for n in names):
                traits['prio'] = 'ps'
            elif any(n in ('HOL', 'FCFSPRIO', 'LCFSPRIO') for n in names):
                traits['prio'] = 'hol'
        return traits

    def _solver_supports(self, name: str) -> bool:
        """Feature-set gate for a candidate name, cached per solver."""
        if not hasattr(self, '_supports_cache'):
            self._supports_cache = {}
        if name in self._supports_cache:
            return self._supports_cache[name]
        ok = True
        try:
            solver = self._create_solver(name)
            # supportsModelMethod is the feature-set gate. supports() is NOT a
            # substitute: on FLD and MAM it validates the method name only, so
            # it answers True for a model the solver rejects at solve time.
            if hasattr(solver, 'supportsModelMethod'):
                verdict = solver.supportsModelMethod('default')
            elif hasattr(solver, 'supports'):
                verdict = solver.supports(self.model)
            else:
                verdict = True
            if isinstance(verdict, tuple):
                verdict = verdict[0]
            ok = bool(verdict)
            if ok and name == 'CTMC':
                # Feature support is necessary but not sufficient for CTMC:
                # the chain must also fit memory.
                from ..solver_ctmc import SolverCTMC
                ok = bool(SolverCTMC.isStateSpaceTractable(
                    self.model, getattr(solver, 'options', None))[0])
        except Exception:
            ok = False
        self._supports_cache[name] = ok
        return ok

    def _choose_ranked(self, order: List[str], method_token: Optional[str] = None) -> Optional[str]:
        """First name in the ranked order that is a candidate and feasible.

        With a method_token the gate tightens from the coarse feature set to
        supportsModelMethod(method_token), which is where the rules a flat
        feature set cannot express live (product form for 'exact', finite
        capacity, NC 'mem' applicability).
        """
        for name in order:
            if name not in self._candidate_solvers:
                continue
            if not self._solver_supports(name):
                continue
            if method_token is not None and not self._solver_supports_method(name, method_token):
                continue
            # Pin the delegate's method: the gate above only ASKED whether the
            # candidate can run method_token, so without this an
            # exactness-gated choice would still run the solver's default
            # (approximate) method.
            self._pinned_methods[name.upper()] = method_token or 'default'
            return name
        return None

    def _solver_supports_method(self, name: str, method_token: str) -> bool:
        """Method-level gate for a candidate name."""
        try:
            solver = self._create_solver(name)
            if not hasattr(solver, 'supportsModelMethod'):
                return True
            verdict = solver.supportsModelMethod(method_token)
            if isinstance(verdict, tuple):
                verdict = verdict[0]
            return bool(verdict)
        except Exception:
            return True

    def _selection_mode(self) -> str:
        """Dispatch mode requested by the caller, lowercased."""
        mode = self.kwargs.get('method')
        if mode is None:
            mode = getattr(self.options, 'method', None)
        return str(mode).lower() if mode else 'default'

    def _select_solver(self) -> str:
        """Selection entry point: the ranked heuristic."""
        return self._select_solver_heuristic()

    def _select_solver_heuristic(self) -> str:
        """Ranked choice for mean-value metrics.

        Global order is MVA > NC > MAM, inverted to NC > MVA on cache models,
        with Fluid promoted for large populations and MAM promoted when the
        traffic is autocorrelated. Mirrors MATLAB chooseAvgSolverHeur.m.
        """
        if self.options.forced_solver:
            return self.options.forced_solver.upper()

        analyzer = self._analyzer
        traits = self._traits()

        # Small populations: an approximation buys nothing there, so take an
        # exact solver whenever one is available, preferring MVA over NC over
        # CTMC (inverted to NC first on caches). The 'exact' token is what makes
        # this a claim rather than a preference: MVA and NC reject it without a
        # product-form solution and CTMC rejects it when the chain does not fit
        # memory.
        total_jobs = analyzer.get_total_jobs()
        if 0 < total_jobs <= EXACT_POPULATION_MAX:
            exact_order = (['NC', 'MVA', 'CTMC'] if traits['cache']
                           else ['MVA', 'NC', 'CTMC'])
            exact = self._choose_ranked(exact_order, 'exact')
            if exact is not None:
                return exact

        if traits['cache']:
            order = ['NC', 'MVA', 'Fluid', 'CTMC', 'LDES']
        elif traits['fcr']:
            if analyzer.get_total_jobs() <= 10:
                order = ['NC', 'CTMC', 'LDES']
            else:
                order = ['NC', 'LDES', 'CTMC']
        elif traits['prio'] == 'preempt':
            order = ['MAM', 'LDES', 'CTMC', 'SSA']
        elif traits['prio'] == 'ps':
            order = ['CTMC', 'LDES', 'SSA']
        elif traits['map']:
            order = ['MAM', 'MVA', 'Fluid', 'LDES']
        elif traits['prio'] == 'hol':
            order = ['MVA', 'MAM', 'Fluid', 'CTMC', 'LDES']
        elif analyzer.get_avg_jobs_per_chain() > 30:
            order = ['Fluid', 'MVA', 'NC']
        elif 0 < total_jobs <= EXACT_POPULATION_MAX:
            # No exact solver was available at this population (tried above), so
            # keep the exact-leaning approximate order.
            order = ['NC', 'MVA', 'MAM']
        elif analyzer.has_all_infinite_servers():
            order = ['MVA', 'NC', 'Fluid']
        else:
            order = ['MVA', 'NC', 'MAM', 'Fluid', 'LDES']

        choice = self._choose_ranked(order)
        if choice is None:
            choice = self._choose_ranked(['MVA', 'NC', 'MAM', 'Fluid', 'LDES',
                                          'CTMC', 'SSA'])
        return choice or 'MVA'

    def _select_solver_exact(self, method_name: str = 'getAvg') -> Optional[str]:
        """Ranked choice restricted to exact solvers.

        Exactness overrides the global MVA > NC order. Mirrors MATLAB
        chooseSolverExact.m.
        """
        analyzer = self._analyzer
        if method_name.startswith('sample'):
            order = ['SSA', 'LDES']
        elif (method_name.startswith('getTranProb') or method_name.startswith('getCdf')
              or method_name.startswith('getTranCdf') or method_name == 'getPerctRespT'
              or method_name == 'getTranAvg'):
            order = ['CTMC']
        elif method_name.startswith('getProb'):
            order = ['NC', 'CTMC'] if analyzer.has_product_form() else ['CTMC']
        elif analyzer.has_product_form() and not analyzer.has_multi_server():
            order = ['NC', 'CTMC']
        else:
            order = ['CTMC', 'NC']
        # Gate on the method-level rule, not just the feature set: NC and MVA
        # reject 'exact' on a non-product-form model, which a flat feature set
        # cannot express.
        return self._choose_ranked(order, 'exact')

    def _select_solver_sim(self, method_name: str = 'getAvg') -> Optional[str]:
        """Ranked choice restricted to simulators: LDES leads, SSA leads for
        event-level sampling. Mirrors MATLAB chooseSolverSim.m."""
        if method_name in ('sample', 'sampleSys'):
            order = ['SSA', 'LDES']
        else:
            order = ['LDES', 'SSA']
        return self._choose_ranked(order)

    def _ctmc_solver(self):
        """The CTMC candidate, built on demand.

        State-space and generator accessors are CTMC-only concepts, so they
        resolve here rather than through the ranked selection.
        """
        from ..solver_ctmc import SolverCTMC
        solver = self._create_solver('CTMC') if 'CTMC' in self._candidate_solvers else None
        if not isinstance(solver, SolverCTMC):
            solver = SolverCTMC(self.model)
        return solver

    def getGenerator(self):
        """Infinitesimal generator, always from SolverCTMC."""
        return self._ctmc_solver().getGenerator()

    def getSymbolicGenerator(self, invert_symbol: bool = False):
        """Symbolic infinitesimal generator, always from SolverCTMC."""
        return self._ctmc_solver().getSymbolicGenerator(invert_symbol)

    def getStateSpace(self):
        """State space, always from SolverCTMC."""
        return self._ctmc_solver().getStateSpace()

    get_generator = getGenerator
    get_symbolic_generator = getSymbolicGenerator
    get_state_space = getStateSpace

    def _create_solver(self, solver_name: str, method: Optional[str] = None) -> Any:
        """Create an instance of the specified solver.

        Args:
            solver_name: Name of the solver to create (e.g., 'MVA', 'NC')
            method: Sub-method to use (e.g., 'lin', 'exact', 'amva'). When None
                the method pinned by _choose_ranked is used, so an
                exactness-gated choice actually runs the exact method.
        """
        solver_name = solver_name.upper()
        if method is None or method == 'default':
            # 'default' is also the caller's "unset": runAnalyzer passes
            # self._solver_method, so the pin has to win over it too.
            method = getattr(self, '_pinned_methods', {}).get(solver_name, 'default')

        # Merge method into kwargs for solvers that support it
        solver_kwargs = dict(self.kwargs)
        if method != 'default':
            solver_kwargs['method'] = method

        # Propagate lang='java' to the delegate so AUTO honors the caller's
        # request to solve via the canonical JAR: the chosen candidate (MVA, NC,
        # CTMC, ...) is itself lang="java"-capable, so AUTO under lang="java"
        # returns JAR-derived results for whichever solver the heuristic picks.
        import os as _os
        lang = self.kwargs.get('lang') or _os.environ.get('LINE_SOLVER_LANG', 'python')
        if lang == 'java':
            solver_kwargs['lang'] = 'java'
        # lang='cpp' propagates the same way, and carries options.arith with it, so
        # the arithmetic backend survives the delegation. The CHOSEN candidate is
        # what runs under lang='cpp': line-cli's own `-s auto` would re-run its
        # chooser on the C++ side and could pick a different engine than the one
        # this heuristic reported, leaving the printed selection unattributable.
        elif lang == 'cpp':
            solver_kwargs['lang'] = 'cpp'
            arith = self.kwargs.get('arith')
            if arith:
                solver_kwargs['arith'] = arith

        if solver_name == 'MVA':
            from ..solver_mva import SolverMVA
            return SolverMVA(self.model, **solver_kwargs)

        elif solver_name == 'LDES':
            from ..wrappers.solver_ldes import SolverLDES, LDESOptions
            options = LDESOptions(**{k: v for k, v in self.kwargs.items()
                                   if hasattr(LDESOptions, k)})
            # Native LDES serializes the Network model (a bare NetworkStruct
            # cannot be saved to model.json), so pass the model, not _sn.
            return SolverLDES(self.model, options)

        elif solver_name == 'NC':
            from ..solver_nc import SolverNC
            return SolverNC(self.model, **solver_kwargs)

        elif solver_name in ('FLUID', 'FLD'):
            from ..solver_fld import SolverFLD
            return SolverFLD(self.model, **solver_kwargs)

        elif solver_name == 'CTMC':
            from ..solver_ctmc import SolverCTMC
            return SolverCTMC(self.model, **solver_kwargs)

        elif solver_name == 'SSA':
            from ..solver_ssa import SolverSSA
            return SolverSSA(self.model, **solver_kwargs)

        elif solver_name == 'JMT':
            from ..wrappers.solver_jmt import SolverJMT
            return SolverJMT(self.model, **solver_kwargs)

        elif solver_name == 'MAM':
            from ..solver_mam import SolverMAM
            return SolverMAM(self.model, **solver_kwargs)

        elif solver_name == 'BA':
            from ..solver_ba import SolverBA
            ba_kwargs = {k: v for k, v in solver_kwargs.items() if k != 'method'}
            return SolverBA(self.model, method or 'auto', **ba_kwargs)

        else:
            # Default to MVA
            from ..solver_mva import SolverMVA
            return SolverMVA(self.model, **solver_kwargs)

    def runAnalyzer(self):
        """
        Run the analysis with the automatically selected solver.

        Returns:
            self for method chaining
        """
        # Solver console: SolverAUTO opens no run of its own, so that the
        # solver it picks narrates its own analysis rather than being nested
        # and silenced. It announces the choice instead.
        from line_solver.api.io import console as _console
        _console.step('AUTO selecting a solver for this model')
        line_debug("AUTO solver starting: method=%s, model=%s",
                   self.options.selection_method,
                   type(self.model).__name__, options=self.options)
        # Select solver based on method. THE METHOD NAME IS RESOLVED BEFORE ANYTHING
        # ELSE: an unqualified algorithm name belongs to whichever family
        # declares it, and a method name no family declares is an error rather than a
        # silent fall-through to the heuristic, which would answer a different
        # algorithm under the caller's name.
        kind, family, submethod = self.resolveMethodToken(
            self.options.selection_method)

        if kind == 'family':
            self._pinned_family = family
            self._pinned_family_method = submethod
        elif self._pinned_family is not None:
            # The model class (LayeredNetwork / Environment) or its Prior
            # parameters already fixed the family; the intent applies inside it.
            self._pinned_family_method = family if family != 'default' else 'default'

        if self._pinned_family is not None:
            if kind == 'intent' and self._model_family in ('layered', 'environment'):
                # The family came from the MODEL CLASS, not from the caller, so
                # the intent still gets a ranked pool -- the same pool MATLAB
                # builds in its class(model) switch.
                return self._run_composite_pool(self._model_family)
            return self._run_family(self._pinned_family, self._pinned_family_method)

        method = self.options.selection_method
        if '.' in method:
            solver_part, sub_method = method.split('.', 1)
            self._solver_method = sub_method
        else:
            solver_part = method
            self._solver_method = 'default'

        # Handle high-level selection methods
        if solver_part == 'sim':
            self._selected_solver_name = self._select_solver_sim() or 'LDES'
        elif solver_part == 'exact':
            # An exact request must not silently return an approximation.
            choice = self._select_solver_exact()
            if choice is None:
                raise RuntimeError("SolverAUTO: no exact solver supports this model; "
                                   "use method 'default' for the approximate heuristic")
            self._selected_solver_name = choice
        elif solver_part == 'fast':
            self._selected_solver_name = self._choose_ranked(
                ['MVA', 'NC', 'Fluid', 'MAM']) or self._select_solver_heuristic()
        elif solver_part == 'accurate':
            self._selected_solver_name = self._choose_ranked(
                ['Fluid', 'MAM', 'CTMC', 'LDES']) or self._select_solver_heuristic()
        elif solver_part == 'bound':
            # Bound analysis - SolverBA with its own family selection
            self._selected_solver_name = 'BA'
            self._solver_method = 'auto'
        elif solver_part in ('default', 'heur', 'auto', 'line'):
            self._selected_solver_name = self._select_solver()
        # Handle specific solver targeting
        elif solver_part.upper() in ['MVA', 'NC', 'CTMC', 'FLUID', 'SSA', 'JMT', 'LDES', 'MAM']:
            self._selected_solver_name = solver_part.upper()
            if self._selected_solver_name == 'FLUID':
                self._selected_solver_name = 'Fluid'
        else:
            # Unknown method - fall back to heuristic
            self._selected_solver_name = self._select_solver()

        if self.options.verbose:
            print(f"SolverAUTO: Selected solver '{self._selected_solver_name}'")
            if self._solver_method != 'default':
                print(f"SolverAUTO: Using method '{self._solver_method}'")
            print(f"SolverAUTO: Candidates were {self._candidate_solvers}")

        # Create and run solver
        line_debug("AUTO attempting solver: %s", self._selected_solver_name, options=self.options)
        try:
            self._selected_solver = self._create_solver(self._selected_solver_name, self._solver_method)
            _console.substep('chose Solver%s', self._selected_solver_name)
            # the delegate narrates its own run: it calls runAnalyzer directly,
            # so its console run is opened by _ensureAvgResults on the accessor
            if hasattr(self._selected_solver, '_ensureAvgResults'):
                self._selected_solver._ensureAvgResults()
            else:
                self._selected_solver.runAnalyzer()
            self._result = self._selected_solver
        except Exception as e:
            # An explicit 'bound' request must not silently degrade to a point
            # estimate from the candidate pool.
            if self._selected_solver_name == 'BA':
                raise
            # Try fallback to other candidates
            for candidate in self._candidate_solvers:
                if candidate != self._selected_solver_name:
                    try:
                        if self.options.verbose:
                            print(f"SolverAUTO: Trying fallback to '{candidate}'")
                        self._selected_solver = self._create_solver(candidate, self._solver_method)
                        self._selected_solver.runAnalyzer()
                        self._selected_solver_name = candidate
                        self._result = self._selected_solver
                        break
                    except Exception:
                        continue
            else:
                raise RuntimeError(f"All solver candidates failed. Last error: {e}")

        return self

    # ------------------------------------------------------------------
    # Family-resolved accessors.
    #
    # These are properties of a particular REPRESENTATION -- the generator, the
    # ODE system, the normalizing constant, the bound interval, the posterior --
    # rather than of the model, so they resolve on their own family instead of
    # on the ranked selection, which is chosen for mean metrics. Mirrors the
    # ctmcSolver/fldSolver/ncSolver/baSolver/uqSolver accessors of MATLAB
    # @SolverAUTO/SolverAUTO.m.
    # ------------------------------------------------------------------

    def _family_solver(self, family: str):
        """Cached solver of a family, built on demand."""
        cache = getattr(self, '_family_solvers', None)
        if cache is None:
            cache = {}
            self._family_solvers = cache
        if family not in cache:
            cache[family] = self._build_family_solver(family, 'default')
        return cache[family]

    def _fld_solver(self):
        return self._family_solver('fluid')

    def _nc_solver(self):
        return self._family_solver('nc')

    def _ba_solver(self):
        if (self._selected_solver is not None
                and type(self._selected_solver).__name__ == 'SolverBA'):
            return self._selected_solver
        return self._build_family_solver('ba', 'auto')

    def _uq_solver(self):
        if (self._selected_solver is not None
                and type(self._selected_solver).__name__ == 'SolverUQ'):
            return self._selected_solver
        return self._family_solver('uq')

    # -- CTMC-resolved -------------------------------------------------
    def getStateSpaceAggr(self, *args, **kwargs):
        """Aggregated state space, from the CTMC family."""
        return self._ctmc_solver().getStateSpaceAggr(*args, **kwargs)

    def getInfGen(self, *args, **kwargs):
        """Infinitesimal generator, from the CTMC family."""
        return self._ctmc_solver().getInfGen(*args, **kwargs)

    def getTransMat(self, *args, **kwargs):
        """Embedded transition matrix, from the CTMC family."""
        return self._ctmc_solver().getTransMat(*args, **kwargs)

    def getMarkedCTMC(self, *args, **kwargs):
        """Marked CTMC, from the CTMC family."""
        return self._ctmc_solver().getMarkedCTMC(*args, **kwargs)

    def getSymbolicSolution(self, *args, **kwargs):
        """Symbolic stationary solution, from the CTMC family."""
        return self._ctmc_solver().getSymbolicSolution(*args, **kwargs)

    def getSensitivity(self, *args, **kwargs):
        """Steady-state parameter sensitivity, from the CTMC family."""
        return self._ctmc_solver().getSensitivity(*args, **kwargs)

    def getSensitivityRanking(self, *args, **kwargs):
        """Sensitivity ranking, from the CTMC family."""
        return self._ctmc_solver().getSensitivityRanking(*args, **kwargs)

    def getCdfSysRespT(self, *args, **kwargs):
        """System response-time distribution, from the CTMC family."""
        return self._ctmc_solver().getCdfSysRespT(*args, **kwargs)

    def getCdfFirstPassT(self, *args, **kwargs):
        """First-passage-time distribution, from the CTMC family."""
        return self._ctmc_solver().getCdfFirstPassT(*args, **kwargs)

    def getFirstPassTMoments(self, *args, **kwargs):
        """First-passage-time moments, from the CTMC family."""
        return self._ctmc_solver().getFirstPassTMoments(*args, **kwargs)

    # -- Fluid-resolved ------------------------------------------------
    def getCdfPT(self, *args, **kwargs):
        """Passage-time distribution, from the fluid family."""
        return self._fld_solver().getCdfPT(*args, **kwargs)

    def getAvgAoI(self, *args, **kwargs):
        """Mean age of information, from the fluid family."""
        return self._fld_solver().getAvgAoI(*args, **kwargs)

    def getCdfAoI(self, *args, **kwargs):
        """Age-of-information distribution, from the fluid family."""
        return self._fld_solver().getCdfAoI(*args, **kwargs)

    def getMoments(self, *args, **kwargs):
        """Queue-length moments, from the fluid family."""
        return self._fld_solver().getMoments(*args, **kwargs)

    def getTranAvgVar(self, *args, **kwargs):
        """Transient mean and variance, from the fluid family."""
        return self._fld_solver().getTranAvgVar(*args, **kwargs)

    def getJacobian(self, *args, **kwargs):
        """ODE Jacobian, from the fluid family."""
        return self._fld_solver().getJacobian(*args, **kwargs)

    def exportODEs(self, *args, **kwargs):
        """Export the ODE system, from the fluid family."""
        return self._fld_solver().exportODEs(*args, **kwargs)

    def getSymbolicDrift(self, *args, **kwargs):
        """Symbolic drift, from the fluid family."""
        return self._fld_solver().getSymbolicDrift(*args, **kwargs)

    # -- NC / BA / UQ resolved ----------------------------------------
    def getNormalizingConstant(self, *args, **kwargs):
        """Normalizing constant, from the NC family."""
        return self._nc_solver().getNormalizingConstant(*args, **kwargs)

    def getBounds(self, *args, **kwargs):
        """Bound interval. Bounds are never in the ranked candidate list, since
        they return an interval rather than a point estimate."""
        return self._ba_solver().getBounds(*args, **kwargs)

    def getBoundsTable(self, *args, **kwargs):
        """Bound interval as a table, from the BA family."""
        return self._ba_solver().getBoundsTable(*args, **kwargs)

    def getAvgBusyPeriod(self, *args, **kwargs):
        """Mean busy period of a subnetwork, from the NC family.

        Daduna (J. ACM 35(3), 1988) is a normalizing-constant transform, so the
        NC family owns it here; SolverLDES measures the same quantity on a
        sample path and is reached by naming that solver. Without this
        delegator `-a busyperiod` through `LINE(model, method='nc')` -- which is
        how the CLI reaches every solver -- found no such attribute and returned
        nothing at all.
        """
        return self._nc_solver().getAvgBusyPeriod(*args, **kwargs)

    def getPosteriorTable(self, *args, **kwargs):
        """Posterior table over the uncertain parameters, from the UQ family."""
        return self._uq_solver().getPosteriorTable(*args, **kwargs)

    def getPosteriorDist(self, *args, **kwargs):
        """Posterior distribution of one metric, from the UQ family."""
        return self._uq_solver().getPosteriorDist(*args, **kwargs)

    # -- Delegated: several families implement these, so the ranked
    #    selection picks the one that supports the model.
    def getAvgReward(self, *args, **kwargs):
        """Mean reward, from whichever candidate supports it."""
        return self._delegate('getAvgReward', *args, **kwargs)

    def getTranReward(self, *args, **kwargs):
        """Transient reward, from whichever candidate supports it."""
        return self._delegate('getTranReward', *args, **kwargs)

    def getStageAvg(self, *args, **kwargs):
        """Per-stage averages of a random environment, from the ENV delegate."""
        return self._delegate('getStageAvg', *args, **kwargs)

    def getTranAvgCoupled(self, *args, **kwargs):
        """Coupled transient of a random environment, from the ENV delegate."""
        return self._delegate('getTranAvgCoupled', *args, **kwargs)

    def getTranAvgDecoupled(self, *args, **kwargs):
        """Decoupled transient of a random environment, from the ENV delegate."""
        return self._delegate('getTranAvgDecoupled', *args, **kwargs)

    def computeSojournCdf(self, *args, **kwargs):
        """Sojourn-time CDF of a random-environment stage."""
        return self._delegate('computeSojournCdf', *args, **kwargs)

    def getSamplePathTable(self, *args, **kwargs):
        """Sample-path table, from the delegate that produced it."""
        return self._delegate('getSamplePathTable', *args, **kwargs)

    def getEnsembleAvg(self, *args, **kwargs):
        """Ensemble averages, from the LN/ENV delegate that decomposes the model."""
        return self._delegate('getEnsembleAvg', *args, **kwargs)

    def getNumberOfModels(self, *args, **kwargs):
        """Number of submodels of the decomposing delegate."""
        return self._delegate('getNumberOfModels', *args, **kwargs)

    def _run_composite_pool(self, model_family: str):
        """Run the ranked candidate pool of a LayeredNetwork or an Environment.

        The pools and their order mirror the class(model) switch of MATLAB
        @SolverAUTO/SolverAUTO.m: a layered model tries the external LQNS binary
        first and then SolverLN over the analytical layer families; an
        environment tries SolverENV over MVA, NC and Fluid. The first candidate
        that both builds and runs wins, so an absent external binary degrades to
        the native decomposition instead of failing.
        """
        if model_family == 'layered':
            pool = [('lqns', 'default'), ('ln', 'nc'), ('ln', 'mva'),
                    ('ln', 'mam'), ('ln', 'fluid')]
        else:
            pool = [('env', 'mva'), ('env', 'nc'), ('env', 'fluid')]
        self._candidate_solvers = ['%s.%s' % (f, m) if m != 'default' else f.upper()
                                   for f, m in pool]
        last_error = None
        for family, method in pool:
            try:
                solver = self._build_family_solver(family, method)
                if hasattr(solver, 'runAnalyzer'):
                    solver.runAnalyzer()
                self._selected_solver = solver
                self._selected_solver_name = family.upper()
                self._solver_method = method
                self._result = solver
                if self.options.verbose:
                    print("SolverAUTO: family '%s' (layers '%s')" % (family, method))
                return self
            except Exception as e:
                last_error = e
                continue
        raise RuntimeError("SolverAUTO: no %s candidate could solve this model. "
                           "Last error: %s" % (model_family, last_error))

    def _run_family(self, family: str, method: str = 'default'):
        """Build and run the pinned family's solver, with no candidate ranking.

        A pinned family is a decision the caller or the model class already
        made, so falling back to another candidate here would answer a
        different question under the same name.
        """
        self._selected_solver_name = family.upper()
        self._solver_method = method
        self._selected_solver = self._build_family_solver(family, method)
        if hasattr(self._selected_solver, 'runAnalyzer'):
            self._selected_solver.runAnalyzer()
        self._result = self._selected_solver
        if self.options.verbose:
            print("SolverAUTO: family '%s' (method '%s')" % (family, method))
        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        # Suppress output from delegated solver
        if hasattr(self._selected_solver, '_table_silent'):
            self._selected_solver._table_silent = True

        df = self._selected_solver.getAvgTable()

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    def getSelectedSolverName(self) -> str:
        """Get the name of the automatically selected solver."""
        if self._selected_solver_name is None:
            self._selected_solver_name = self._select_solver()
        return self._selected_solver_name

    def getCandidateSolverNames(self) -> List[str]:
        """Get list of candidate solver names for this model."""
        return self._candidate_solvers.copy()

    def getModel(self) -> Any:
        """Get the network model being solved."""
        return self.model

    def getStruct(self, *args, **kwargs) -> Any:
        """Get the network structure from the model."""
        if hasattr(self.model, 'getStruct'):
            return self.model.getStruct(*args, **kwargs)
        return self._sn

    def getName(self) -> str:
        """Get solver name."""
        # The reported name stays 'SolverAuto', matching MATLAB's self.name and
        # the JAR's super(model, "SolverAuto", options): only the class
        # identifier was renamed.
        return 'SolverAuto'

    def setForcedSolver(self, solver_name: str) -> None:
        """Force a specific solver to be used."""
        self.options.forced_solver = solver_name

    def setSelectionMethod(self, method: str) -> None:
        """Set the solver selection method."""
        self.options.selection_method = method

    def setOptions(self, options: Optional[SolverAUTOOptions] = None, **kwargs) -> None:
        """
        Set solver options.

        Args:
            options: SolverAUTOOptions instance or None to update current
            **kwargs: Individual option updates (selection_method, forced_solver, verbose)
        """
        if options is not None:
            self.options = options
        else:
            # Update individual options
            if 'selection_method' in kwargs:
                self.options.selection_method = kwargs['selection_method']
            if 'forced_solver' in kwargs:
                self.options.forced_solver = kwargs['forced_solver']
            if 'verbose' in kwargs:
                self.options.verbose = kwargs['verbose']

    def supports(self, model: Any) -> bool:
        """
        Check if SolverAUTO can handle this model.

        SolverAUTO supports any model that at least one candidate solver supports.

        Args:
            model: Network model to check

        Returns:
            True if any candidate solver supports the model, False otherwise
        """
        if not self._candidate_solvers:
            return False

        for candidate in self._candidate_solvers:
            try:
                solver = self._create_solver(candidate)
                if hasattr(solver, 'supports') and callable(solver.supports):
                    if solver.supports(model):
                        return True
            except Exception:
                pass

        # No candidate supports the model. This previously returned
        # len(self._candidate_solvers) > 0, i.e. True whenever any candidate
        # merely existed, which overrode the verdict of the loop above and made
        # AUTO claim every model regardless of the features it used.
        return False

    def hasResults(self) -> bool:
        """Check if analysis has been completed and results are available."""
        return self._result is not None

    def isStochastic(self):
        """If a delegate solver has been selected, classify by it (once run,
        it knows the method it resolved at runtime). Before selection, the
        delegate choice is unknown, so classify conservatively: true if any
        candidate solver is simulation-based.
        """
        if self._selected_solver is not None and hasattr(self._selected_solver, 'isStochastic'):
            return self._selected_solver.isStochastic()
        return any(name in ('SSA', 'LDES') for name in (self._candidate_solvers or []))

    is_stochastic = isStochastic

    def getResults(self) -> Any:
        """
        Get the solver results object.

        Returns:
            The selected solver or None if not yet run
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result

    def listAllMethods(self) -> List[str]:
        """
        Every token the CONSTRUCTOR accepts, INDEPENDENT of the model.

        The selection intents, the method families, and every family method in
        its qualified form. The unqualified form of a family method is accepted
        too (resolveMethodToken looks it up through familyDeclaringMethod), and
        is left out here only to keep the list unambiguous.

        THIS is the list a method-NAME check must gate on. listValidMethods
        narrows it to the model in hand, and gating a name check on that would
        replace a rejection the delegate would have EXPLAINED with a flat "the
        method is unsupported by this solver". Mirrors MATLAB
        @SolverAUTO/listAllMethods.

        Returns:
            Sorted list of accepted method names.
        """
        methods = list(self.selectionIntents()) + ['auto', 'line']
        for family in self.familyNames():
            try:
                probe = self._build_family_solver(family, 'default', probe=True)
                lister = getattr(probe, 'listAllMethods', None)
                declared = lister() if lister is not None else probe.listValidMethods()
            except Exception:
                # A family that cannot even be INSTANTIATED here contributes
                # nothing.
                continue
            methods.append(family)
            for m in (declared or []):
                methods.append('%s.%s' % (family, m))
        # Aliases familyAlias resolves but familyNames does not spell out.
        methods.extend(['fld', 'des', 'lqsim'])
        return sorted(set(methods))

    list_all_methods = listAllMethods

    # ---------------------------------------------------------------- #
    # findSolver: which solvers and methods can analyze this model
    # ---------------------------------------------------------------- #

    @staticmethod
    def metricGroups() -> List[str]:
        """The measure groups findSolver reports on, in report order.

        A group is a family of accessors that stand or fall together: a solver
        that returns getCdfRespT returns getCdfPassT and getPerctRespT as well,
        because all three read the same passage time, so listing the three
        separately would say nothing extra.

        Mirrors MATLAB SolverAUTO.metricGroups.
        """
        return ['avg', 'tran', 'cdf', 'prob', 'tranprob', 'sample',
                'cache', 'loss', 'orbit', 'moment', 'sens']

    @staticmethod
    def metricGroupOf(name) -> str:
        """The measure group an accessor belongs to, '' when NAME names none.

        A group name maps to itself, so findSolver('cdf') and
        findSolver('getCdfRespT') ask the same question.

        THIS IS NOT _choose_solver_for_method's TABLE, although both are keyed
        by accessor name. That one maps an accessor to a RANKING, i.e. which
        candidate should be preferred; this one maps it to a CAPABILITY
        question, i.e. which candidates can answer it at all. The two differ
        wherever a family can serve a measure but is never the one AUTO would
        pick for it.

        Mirrors MATLAB SolverAUTO.metricGroupOf.
        """
        if not name or not isinstance(name, str):
            return ''
        if name in SolverAUTO.metricGroups():
            return name
        if name.lower() in ('any', 'all'):
            return ''
        if name in ('getTranAvg', 'getTranAvgVar', 'tranAvg'):
            return 'tran'
        if name in ('getCdfRespT', 'getCdfPassT', 'getPerctRespT',
                    'getTranCdfPassT', 'getTranCdfRespT', 'getCdfSysRespT'):
            return 'cdf'
        if name in ('getTranProb', 'getTranProbSys', 'getTranProbAggr',
                    'getTranProbSysAggr'):
            return 'tranprob'
        if name in ('getProb', 'getProbAggr', 'getProbSys', 'getProbSysAggr',
                    'getProbMarg', 'getProbNormConstAggr'):
            return 'prob'
        if name in ('sample', 'sampleSys', 'sampleAggr', 'sampleSysAggr'):
            return 'sample'
        if name in ('getAvgCacheTable', 'getAvgCacheT', 'getAvgItemTable',
                    'getAvgItemT', 'cacheAvgT', 'itemAvgT', 'aCaT', 'aIT'):
            return 'cache'
        if name in ('getAvgLossTable', 'getAvgLossT', 'getAvgRegionLossTable',
                    'getAvgRegionLossT', 'lossAvgT', 'regionLossAvgT', 'aLT', 'aRLT'):
            return 'loss'
        if name in ('getAvgOrbitTable', 'getAvgOrbitT', 'getAvgOrbit',
                    'orbitAvgT', 'aOT'):
            return 'orbit'
        if name in ('getMomentTable', 'getMomentChainTable', 'getMomentStationTable',
                    'getMomentT', 'getMomentChainT', 'getMomentStationT',
                    'momentT', 'momentChainT', 'momentStationT', 'mT', 'mCT', 'mST'):
            return 'moment'
        if name in ('getSensitivityTable', 'getSensitivityT', 'sensitivityT', 'sT',
                    'getSensitivity', 'getSensitivityRanking'):
            return 'sens'
        # Everything else in the accessor surface is a mean measure: getAvg,
        # its chain, node and system forms, their handles and short aliases.
        if name.startswith('getAvg') or name.startswith('avg') or name in (
                'getAvgSysRespT', 'getAvgSysTput', 'aT', 'aNT', 'aCT', 'aST', 'aNCT'):
            return 'avg'
        return ''

    @staticmethod
    def familyMetrics(family: str) -> List[str]:
        """The measure groups a method family can answer.

        Every family answers 'avg', which is what a solver is for; the rest is
        the capability declaration this class owns.

        SOURCES, so that a claim here can be checked rather than trusted:
        'tran' is supportsTransientAnalysis, which FLD, CTMC, LDES and JMT
        override to true and no one else does. 'cdf', 'prob', 'tranprob' and
        'sample' are the families that carry an implementation of the
        corresponding accessor rather than inheriting the base refusal. The
        remaining five groups are computed by NetworkSolver from a solver's own
        results, so no per-solver method marks them: their lists are the
        rankings _choose_solver_for_method holds for the same accessors, which
        is where AUTO already records who can serve them.

        A family that gains or loses a measure must be edited here in the same
        change, the way a solver that gains a feature is edited into its
        feature set: an omission here does not fail, it silently hides the
        family from a caller asking for that measure.

        Mirrors MATLAB SolverAUTO.familyMetrics.
        """
        table = {
            'mva': ['avg', 'prob', 'cache', 'orbit', 'moment', 'sens'],
            'nc': ['avg', 'cdf', 'prob', 'cache', 'moment', 'sens'],
            'ctmc': ['avg', 'tran', 'cdf', 'prob', 'tranprob', 'sample',
                     'cache', 'loss', 'orbit', 'moment'],
            'fluid': ['avg', 'tran', 'cdf', 'prob', 'cache', 'sens'],
            'mam': ['avg', 'cdf'],
            # The RCAT fixed point reports means only; the passage-time law it
            # answers is the base exponential fit, not its own.
            'ag': ['avg'],
            # A bound brackets the mean measures and nothing else.
            'ba': ['avg'],
            'ssa': ['avg', 'cdf', 'prob', 'sample', 'loss'],
            'ldes': ['avg', 'tran', 'cdf', 'prob', 'sample', 'cache', 'loss', 'orbit'],
            'jmt': ['avg', 'tran', 'cdf', 'prob', 'tranprob', 'sample'],
            'qns': ['avg'],
            'ln': ['avg', 'tran', 'cdf', 'sens'],
            'env': ['avg', 'tran'],
            'lqns': ['avg'],
            'uq': ['avg'],
        }
        return list(table.get(family, ['avg']))

    @staticmethod
    def methodClass(family: str, method: str, is_stochastic: bool,
                    is_product_form: bool, is_qbd_shape: bool,
                    has_cache: bool = False) -> str:
        """What KIND of answer a method returns: 'exact', 'approx', 'bound' or
        'simulation'.

        'simulation' is not decided here: is_stochastic is the solver's own
        isStochasticMethod, which already tokenizes qualified and
        runtime-resolved names and is the only place that knowledge lives.

        'exact' IS CLAIMED ONLY WHERE IT IS TRUE OF THIS MODEL, never of the
        algorithm in the abstract. Exactness of a normalizing constant or of
        mean value analysis is a property of the product-form model it is
        computed on, and of the QBD shape for the matrix analytic methods, so
        both conditions are passed in and a method that needs one reports
        'approx' without it. The bias is deliberate: an under-claimed 'approx'
        costs a user a better method they could have had, an over-claimed
        'exact' costs them a wrong number they trusted.

        A CACHE IS THE THIRD CONDITION, and it was the over-claim the bias above
        exists to prevent. sn_has_product_form answers about the QUEUEING network
        and knows nothing of a cache: the hit/miss split is a class switch whose
        probabilities are not routing data but the output of a cache model, so a
        network holding one reads as product form and 'mva.exact' was labelled
        exact on it. Measured on the tut06 shape with an LRU cache: exact MVA
        returns QLen 0.2516 at the hit station where the CTMC returns 0.3022 and
        simulation 0.3023, a 17% error under a label that says there is none.
        The analytic families are conditioned on it; SolverCTMC is NOT, because
        its state space carries the cache contents and it is exact there, which
        is what the two numbers above show.

        Mirrors MATLAB SolverAUTO.methodClass.
        """
        if family == 'ba':
            # Bounds are what SolverBA is for; every one of its methods returns
            # a bracket rather than an estimate.
            return 'bound'
        if is_stochastic:
            return 'simulation'
        exact_if = lambda cond: 'exact' if cond else 'approx'  # noqa: E731
        if family == 'ctmc':
            # The generator is solved as written, so every state-space route is
            # exact. 'cftp.approx' says in its own name that it is not, and
            # 'mdd' is exact on a product-form model and an approximation
            # otherwise.
            if method == 'cftp.approx':
                return 'approx'
            if method == 'mdd':
                return exact_if(is_product_form)
            return 'exact'
        if family == 'nc':
            # The normalizing-constant routes that evaluate G exactly rather
            # than expanding or estimating it. The asymptotic expansions
            # (pana, le, kt, bk, gm, ...) and the non-product-form
            # 'morrison' are approximations by construction and are left out.
            if method in ('exact', 'divdiff', 'ca', 'comom', 'comomld',
                          'rec', 'ms', 'cub', 'rgf'):
                return exact_if(is_product_form and not has_cache)
            return 'approx'
        if family == 'mva':
            # Exact MVA; every 'amva.*' arm is an approximation, and so are the
            # open-network QNA transforms.
            if method in ('exact', 'mva'):
                return exact_if(is_product_form and not has_cache)
            return 'approx'
        if family == 'jmt':
            # JMVA's exact algorithms. 'jsim' and 'replication' are simulation
            # and never reach here.
            if method in ('jmva.mva', 'jmva.recal', 'jmva.comom', 'jmva.treeconv'):
                return exact_if(is_product_form and not has_cache)
            return 'approx'
        if family == 'ag':
            # Every RCAT arm estimates the reversed rate of each synchronising
            # action and iterates to a fixed point, an approximation by
            # construction; SolverAG's 'exact' is a vestigial alias that warns
            # and runs 'inap', so nothing here is claimed exact.
            return 'approx'
        if family == 'mam':
            # The QBD is solved exactly on the shape it is stated for, one
            # queueing station fed by a Source. Everything named 'dec.*' is a
            # decomposition of a larger network into such queues and is
            # therefore an approximation of it.
            if method in ('default', 'mna', 'ldqbd', 'bgchain', 'retrial'):
                return exact_if(is_qbd_shape)
            return 'approx'
        return 'approx'

    def _exactness_conditions(self):
        """The three model properties an exactness claim can rest on, evaluated
        once per report; see methodClass."""
        is_product_form = False
        is_qbd_shape = False
        has_cache = False
        model = self.model
        try:
            if hasattr(model, 'has_product_form_solution'):
                is_product_form = bool(model.has_product_form_solution())
            sn = model.get_struct()
            nsources = sum(1 for t in sn.nodetype if t == NodeType.SOURCE)
            njobs = np.asarray(sn.njobs, dtype=float).ravel()
            is_qbd_shape = bool(njobs.size and np.all(np.isinf(njobs))
                                and (int(sn.nstations) - nsources) == 1)
            has_cache = any(t == NodeType.CACHE for t in sn.nodetype)
        except Exception:
            # A model whose struct cannot be refreshed here answers 'approx'
            # everywhere, which is the safe direction: see methodClass.
            pass
        return is_product_form, is_qbd_shape, has_cache

    def findSolver(self, metric: str = '', showAll: bool = False) -> pd.DataFrame:
        """Which solvers and solver methods can analyze this model, and for the
        ones that cannot, why not.

        One row per (family, method) pair AUTO can be asked for, with columns

            Solver    the method family, 'mva', 'ctmc', 'ldes', ...
            Method    the method name to pass, 'mva.exact'
            Runnable  True when the model passes that method's own support gate
            Class     'exact', 'approx', 'bound' or 'simulation'
            Metrics   the measure groups the family answers, see metricGroups
            Reason    why a refused pair was refused, '' when Runnable

        Args:
            metric: narrows the report to the pairs that answer one measure,
                named either by its group ('cdf') or by the accessor that
                returns it ('getCdfRespT'). '' or 'any' keeps every pair.
            showAll: keep the refused pairs too. By default only the runnable
                ones are listed, since a caller asking what it can run has no
                use for the two hundred rows that say it cannot.

        THE GATE IS NOT A SECOND ONE. It is the gate _choose_solver_ranked
        applies before delegating, asked of every candidate instead of of the
        first feasible one, which is exactly what listValidMethods already did
        -- that method is now the Method column of the runnable rows, so the
        two cannot disagree. What is new is that the REASON the gate produced
        is kept rather than discarded, and that the answer carries the two
        facts a caller needs in order to choose among the survivors: whether
        the method is exact on this model, and which measures it can report.

        WHY A REASON HAS TO BE RECONSTRUCTED for some rows. The base
        supportsModelMethod returns a reason only when the solver diverges per
        method; when it does not, it falls back to supports(model), which
        answers with a bare bool. That is enough for a gate, which only has to
        stop the run, and not enough for a report, whose whole content is the
        explanation. So a refused row with no reason is re-asked against the
        solver's own feature set, which is where the offending names are.

        Mirrors MATLAB @SolverAUTO/findSolver.m.

        Returns:
            A pandas DataFrame with the six columns above.
        """
        group = SolverAUTO.metricGroupOf(metric)
        if not group and metric and str(metric).lower() not in ('any', 'all'):
            raise ValueError(
                "'%s' names no measure. Pass a group (%s) or the accessor that "
                "returns it, e.g. 'getCdfRespT'."
                % (metric, ', '.join(SolverAUTO.metricGroups())))

        is_product_form, is_qbd_shape, has_cache = self._exactness_conditions()

        with SolverAUTO.silenced():
            rows = self._find_solver_rows(group, showAll, is_product_form, is_qbd_shape,
                                          has_cache)
        return pd.DataFrame(rows, columns=['Solver', 'Method', 'Runnable',
                                           'Class', 'Metrics', 'Reason'])

    @staticmethod
    @contextmanager
    def silenced():
        """Run a block with the logger silent, restoring the level after.

        A REPORT MUST NOT PRINT. Asking a solver whether it supports the model
        runs supports_via_featureset, whose SolverFeatureSet.supports emits a
        warning naming the missing feature -- a side effect that is right on
        the solve path, where nobody asked to be told, and wrong here, where
        every refused row would raise one and the answer IS the table.

        It has to wrap the construction of SolverAUTO as well as the walk: the
        constructor probes every candidate with supports(model), so a model one
        of them refuses warns before findSolver is even called. That is why the
        model-side findSolver enters this before building the solver.
        """
        from ...api.io.logging import _logger
        from ...constants import VerboseLevel
        saved = _logger.verbose
        _logger.verbose = VerboseLevel.SILENT
        try:
            yield
        finally:
            _logger.verbose = saved

    def _find_solver_rows(self, group, showAll, is_product_form, is_qbd_shape, has_cache):
        """The walk behind findSolver, one row per (family, method) pair."""
        rows = []
        for family in self.familyNames():
            if not self.familyAcceptsModelClass(family, self.model):
                continue
            groups = SolverAUTO.familyMetrics(family)
            if group and group not in groups:
                continue
            try:
                probe = self._build_family_solver(family, 'default', probe=True)
                declared = probe.listValidMethods()
            except Exception:
                # A family that cannot even be instantiated here (SolverLQNS
                # without the lqns binary, SolverLN on a flat Network)
                # contributes nothing: there is no solver to report on and no
                # gate to ask.
                continue
            metric_list = ','.join(groups)
            # Once per family, not once per refused row: the flat feature set
            # does not vary with the method, and a family that refuses every one
            # of its forty methods would otherwise recompute the same
            # comparison forty times.
            fam_feature_reason = self._feature_reason(probe)
            for name in (declared or []):
                if name.startswith(family + '.'):
                    # A spelling already qualified with its own family.
                    # SolverFLD declares both 'dae' and 'fluid.dae' so that its
                    # own gate takes either, and prefixing the family again
                    # yields 'fluid.fluid.dae': a token that does resolve, but
                    # that names the same method twice and would double every
                    # fluid row of this report.
                    continue
                if SolverAUTO.isMethodAlias(family, name, declared):
                    # The same duplication under a different prefix. SolverMVA
                    # advertises every AMVA name twice, plain and 'amva.'-
                    # prefixed, and its dispatch strips the prefix, so the two
                    # spellings are one algorithm; that alone was 20 of the 49
                    # mva rows of a report. The plain spelling is the one kept.
                    continue
                ok, reason = self._gate_reason(probe, name)
                if not ok and self._is_generic_reason(reason):
                    # A SPECIFIC GATE REASON WINS; the feature-set names only
                    # replace a generic one ("Some features are not supported
                    # by the NC solver", which names none) or fill an empty one
                    # (the flat supports(model) answers with a bare bool).
                    # Getting this the other way round was worse than saying
                    # nothing: UQ's feature set answers a DIFFERENT question --
                    # it declares the Prior UQ itself consumes, every other
                    # feature being the INNER solver's to accept -- so
                    # comparing it against the model listed every feature the
                    # model uses and buried the real reason, that the model
                    # carries no uncertain parameter at all.
                    reason = fam_feature_reason or reason or (
                        'This solver refuses the model through its own '
                        'structural check.')
                elif ok:
                    reason = ''
                if not ok and not showAll:
                    continue
                rows.append({
                    'Solver': family,
                    'Method': '%s.%s' % (family, name),
                    'Runnable': bool(ok),
                    'Class': SolverAUTO.methodClass(
                        family, name, self._stochastic_verdict(probe, name),
                        is_product_form, is_qbd_shape, has_cache),
                    'Metrics': metric_list,
                    'Reason': reason,
                })
        return rows

    find_solver = findSolver

    @staticmethod
    def _is_generic_reason(reason) -> bool:
        """Does the reason say only that SOME feature is unsupported, without
        naming one?

        That is the sentence SolverNC and the shared base gate return, and the
        one worth replacing with the offending feature names. Empty counts: a
        gate that answered with a bare bool said nothing at all.
        """
        if not reason:
            return True
        return ('features are not supported' in reason.lower()
                and '(feature:' not in reason)

    @staticmethod
    def _gate_reason(probe, method):
        """The method-level gate, asked without letting it raise.

        The base supportsModelMethod already falls back to the flat
        supports(model) for a solver that does not diverge per method, so this
        one call is both gates.
        """
        gate = getattr(probe, 'supportsModelMethod', None)
        if gate is None:
            try:
                return bool(probe.supports(probe.model)), ''
            except Exception:
                return True, ''  # no claim, no gate
        try:
            answer = gate(method)
        except Exception as e:
            # A gate that raises has said something, and it is the only thing
            # it can say about this pair; reporting it beats swallowing it and
            # calling the pair runnable.
            return False, str(e)
        if isinstance(answer, tuple):
            return bool(answer[0]), (answer[1] or '')
        return bool(answer), ''

    @staticmethod
    def _stochastic_verdict(probe, method) -> bool:
        """isStochasticMethod asked without letting it raise; a solver that
        cannot classify a name is taken at its class default, deterministic."""
        try:
            return bool(probe.isStochasticMethod(method))
        except Exception:
            return False

    @staticmethod
    def _feature_reason(probe) -> str:
        """The offending feature names, or '' when the feature set accepts the
        model.

        The empty answer is meaningful and not a failure: it says the refusal
        came from somewhere the feature set cannot see, so the caller should
        keep whatever the gate itself said.
        """
        get_fs = getattr(type(probe), 'getFeatureSet', None)
        if get_fs is None:
            return ''
        try:
            from ..base import SolverFeatureSet
            feat_supported = get_fs()
            if not isinstance(feat_supported, SolverFeatureSet):
                fs = SolverFeatureSet()
                fs.set_true(list(feat_supported))
                feat_supported = fs
            # A SOLVER THAT DECLARES NOTHING HAS NO ENVELOPE, and a
            # missing-feature list against an empty set is not an explanation:
            # it is every feature the model uses. UQ is the case -- it computes
            # nothing itself, it expands a Prior and runs another solver at each
            # design point -- and it refused an M/M/1 with "(feature: Sink,
            # Source, Exp, RoutingStrategy_PROB, SchedStrategy_FCFS,
            # OpenClass)", which tells a user nothing they can act on. Such a
            # refusal is structural, so the caller keeps the gate's own words.
            if not any(feat_supported.list.values()):
                return ''
            model = probe.model
            if hasattr(model, 'get_used_lang_features'):
                feat_used = model.get_used_lang_features()
            else:
                feat_used = model.getUsedLangFeatures()
            _, reason = SolverFeatureSet.supports_with_reason(feat_supported, feat_used)
        except Exception:
            return ''
        return reason or ''

    def listValidMethods(self) -> List[str]:
        """
        The method names of listAllMethods that THIS MODEL can actually run.

        IT IS THE RUNNABLE ROWS OF findSolver, projected onto their Method
        column. The narrowing used to be written out a second time here, and a
        second copy of one gate is how two answers to one question start to
        differ; findSolver owns it now, and this adds only the method names that name
        no single method: the selection intents, which name a RANKING rather
        than an algorithm, each family's bare method name, and the family aliases.

        The gate findSolver applies is the one _choose_solver_ranked applies
        before delegating, asked of every candidate instead of the first
        feasible one: a method whose supportsModelMethod refuses the model is
        not offered, and a family with no method-level gate is judged by its
        flat feature set through supports(model). The method-level gate is
        where the rules a feature set cannot express live (product form for
        'exact', a binding finite buffer, NC 'mem' applicability).

        Without it this returned the token universe regardless of the model:
        215 entries for the two-station BAS-blocking model of cqn_bas_blocking,
        naming all 37 SolverNC methods and all 8 SolverLQNS methods although
        both refuse that model method by method. A caller enumerating the list
        was being invited to ask for an analysis no candidate would perform.

        A family whose every method is refused loses its bare method name too: 'nc'
        alone delegates to SolverNC, which is exactly the rejection the
        per-method gate just returned. That falls out of the projection, since
        a family with no runnable row contributes no row to take its name from.

        Mirrors MATLAB @SolverAUTO/listValidMethods.

        Returns:
            Sorted list of method names this model can run.
        """
        methods = list(self.selectionIntents()) + ['auto', 'line']
        T = self.findSolver()
        for _, row in T.iterrows():
            methods.append(row['Method'])
            methods.append(row['Solver'])
        # Aliases familyAlias resolves but familyNames does not spell out.
        for alias, family in (('fld', 'fluid'), ('des', 'ldes'), ('lqsim', 'lqns')):
            if family in methods:
                methods.append(alias)
        return sorted(set(methods))

    @staticmethod
    def _gate_verdict(answer) -> bool:
        """Normalise a support answer: the gates return either a bool or a
        (bool, reason) pair, and both shapes are in use across the solvers."""
        if isinstance(answer, tuple):
            return bool(answer[0])
        return bool(answer)

    # ========== Delegation Methods for NetworkSolver Parity ==========

    def _choose_solver_for_method(self, method_name: str) -> Optional[str]:
        """Metric-aware ranked choice, mirroring MATLAB chooseSolverHeur.m.

        LDES is the only simulation candidate, MVA leads NC leads MAM for analytical
        solvers (inverted on caches), and Fluid leads where a smooth answer is
        wanted (sensitivity, distributions, transients).

        Args:
            method_name: Name of the method to optimize for

        Returns:
            Solver name to prioritize, or None to keep the general heuristic
        """
        analyzer = self._analyzer

        if method_name in ('getTranAvg', 'getTranCdfPassT', 'getTranCdfRespT'):
            order = ['Fluid', 'LDES']
        elif method_name in ('getCdfRespT', 'getCdfPassT', 'getPerctRespT'):
            if (analyzer.has_homogeneous_scheduling(SchedStrategy.FCFS)
                    and analyzer.has_product_form()):
                order = ['NC', 'Fluid', 'LDES']
            else:
                order = ['Fluid', 'LDES']
        elif method_name.startswith('getTranProb'):
            order = ['CTMC']
        elif method_name in ('sample', 'sampleSys'):
            order = ['SSA', 'LDES']
        elif method_name in ('sampleAggr', 'sampleSysAggr'):
            order = ['LDES', 'SSA']
        elif method_name.startswith('getProb'):
            if analyzer.has_product_form():
                order = ['NC', 'CTMC', 'LDES']
            else:
                order = ['CTMC', 'LDES']
        elif method_name in ('getAvgCacheTable', 'getAvgItemTable'):
            # Cache metrics invert the analytical order: NC leads MVA here.
            order = ['NC', 'MVA', 'Fluid', 'CTMC', 'LDES']
        elif method_name in ('getAvgLossTable', 'getAvgRegionLossTable'):
            order = ['LDES', 'CTMC', 'SSA']
        elif method_name in ('getAvgOrbitTable', 'getAvgOrbit'):
            order = ['MVA', 'CTMC', 'LDES']
        elif method_name in ('getMomentTable', 'getMomentChainTable', 'getMomentStationTable'):
            order = ['MVA', 'NC', 'CTMC', 'LDES']
        elif method_name == 'getSensitivityTable':
            # Sensitivity wants a differentiable model, which is what Fluid gives.
            order = ['Fluid', 'MVA', 'NC']
        else:
            return None

        return self._choose_ranked(order)

    def _delegate(self, method_name: str, *args, num_returns: int = 1, **kwargs) -> Any:
        """
        Delegate method call to the selected solver with intelligent fallback.

        Based on MATLAB SolverAUTO.delegate() logic:
        1. Determine proposed solvers: best solver for method + candidates
        2. Try each solver in order
        3. Support variable return values (1-7)
        4. On failure, try next candidate

        Args:
            method_name: Name of the method to call
            *args: Positional arguments
            num_returns: Expected number of return values (default: 1)
            **kwargs: Keyword arguments

        Returns:
            Result from the selected solver's method, handling variable returns
        """
        if self._selected_solver is None:
            self.runAnalyzer()

        # Determine proposed solvers to try
        proposed_solvers = []

        # First priority: method-specific best solver
        best_for_method = self._choose_solver_for_method(method_name)
        if best_for_method and best_for_method != self._selected_solver_name:
            try:
                best_solver = self._create_solver(best_for_method)
                if best_solver.supports(self.model) if hasattr(best_solver, 'supports') else True:
                    proposed_solvers.append(best_solver)
            except Exception:
                pass  # If can't create, skip to next

        # Second priority: currently selected solver
        proposed_solvers.append(self._selected_solver)

        # Third priority: other candidates
        for candidate in self._candidate_solvers:
            if candidate != self._selected_solver_name and candidate != best_for_method:
                try:
                    solver = self._create_solver(candidate)
                    proposed_solvers.append(solver)
                except Exception:
                    pass  # Skip if can't create

        # Try each proposed solver in order
        last_error = None
        for solver in proposed_solvers:
            try:
                if not (hasattr(solver, 'supports') and callable(solver.supports)):
                    # Solver doesn't support method or check fails, try anyway
                    pass
                elif not solver.supports(self.model):
                    # Skip this solver if it doesn't support the model
                    continue

                if not hasattr(solver, method_name):
                    # Method doesn't exist on this solver, try next
                    if self.options.verbose:
                        print(f"SolverAUTO: Method '{method_name}' unsupported by {solver.__class__.__name__}")
                    continue

                # Call the method
                method = getattr(solver, method_name)
                result = method(*args, **kwargs)

                if self.options.verbose:
                    print(f"SolverAUTO: Method '{method_name}' succeeded with {solver.__class__.__name__}")

                # Update selected solver if different
                if solver != self._selected_solver:
                    self._selected_solver = solver
                    self._selected_solver_name = solver.__class__.__name__.replace('SolverNative', '').replace('Solver', '')

                return result

            except Exception as e:
                last_error = e
                if self.options.verbose:
                    print(f"SolverAUTO: {solver.__class__.__name__} failed for '{method_name}': {str(e)}")
                continue

        # All solvers failed
        error_msg = f"All solvers failed for method '{method_name}'"
        if last_error:
            error_msg += f". Last error: {str(last_error)}"
        raise RuntimeError(error_msg)

    # ========== Basic Metric Methods ==========

    def getAvg(self):
        """Get all average metrics (QN, UN, RN, TN, AN, WN)."""
        return self._delegate('getAvg')

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths at steady-state."""
        return self._delegate('getAvgQLen')

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations at steady-state."""
        return self._delegate('getAvgUtil')

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times at steady-state."""
        return self._delegate('getAvgRespT')

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times at steady-state."""
        return self._delegate('getAvgResidT')

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (queue time excluding service)."""
        return self._delegate('getAvgWaitT')

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs at steady-state."""
        return self._delegate('getAvgTput')

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates at steady-state."""
        return self._delegate('getAvgArvR')

    # ========== Chain Methods ==========

    def getAvgChain(self):
        """Get average metrics by chain."""
        return self._delegate('getAvgChain')

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue length by chain."""
        return self._delegate('getAvgQLenChain')

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilization by chain."""
        return self._delegate('getAvgUtilChain')

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response time by chain."""
        return self._delegate('getAvgRespTChain')

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence time by chain."""
        return self._delegate('getAvgResidTChain')

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughput by chain."""
        return self._delegate('getAvgTputChain')

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rate by chain."""
        return self._delegate('getAvgArvRChain')

    # ========== Node Methods ==========

    def getAvgNode(self):
        """Get average metrics by node."""
        return self._delegate('getAvgNode')

    def getAvgNodeQLenChain(self) -> np.ndarray:
        """Get average node queue length by chain."""
        return self._delegate('getAvgNodeQLenChain')

    def getAvgNodeUtilChain(self) -> np.ndarray:
        """Get average node utilization by chain."""
        return self._delegate('getAvgNodeUtilChain')

    def getAvgNodeRespTChain(self) -> np.ndarray:
        """Get average node response time by chain."""
        return self._delegate('getAvgNodeRespTChain')

    def getAvgNodeResidTChain(self) -> np.ndarray:
        """Get average node residence time by chain."""
        return self._delegate('getAvgNodeResidTChain')

    def getAvgNodeTputChain(self) -> np.ndarray:
        """Get average node throughput by chain."""
        return self._delegate('getAvgNodeTputChain')

    def getAvgNodeArvRChain(self) -> np.ndarray:
        """Get average node arrival rate by chain."""
        return self._delegate('getAvgNodeArvRChain')

    # ========== System Methods ==========

    def getAvgSys(self):
        """Get system-level average metrics."""
        return self._delegate('getAvgSys')

    def getAvgSysRespT(self) -> float:
        """Get system average response time."""
        return self._delegate('getAvgSysRespT')

    def getAvgSysTput(self) -> float:
        """Get system average throughput."""
        return self._delegate('getAvgSysTput')

    # ========== Table Methods ==========

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics table by chain."""
        return self._delegate('getAvgChainTable')

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system average metrics table."""
        return self._delegate('getAvgSysTable')

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics table by node."""
        return self._delegate('getAvgNodeTable')

    def getAvgQLenTable(self) -> pd.DataFrame:
        """Get average queue length table."""
        return self._delegate('getAvgQLenTable')

    def getAvgUtilTable(self) -> pd.DataFrame:
        """Get average utilization table."""
        return self._delegate('getAvgUtilTable')

    def getAvgRespTTable(self) -> pd.DataFrame:
        """Get average response time table."""
        return self._delegate('getAvgRespTTable')

    def getAvgTputTable(self) -> pd.DataFrame:
        """Get average throughput table."""
        return self._delegate('getAvgTputTable')

    def getAvgNodeChainTable(self) -> pd.DataFrame:
        """Get average metrics table by node and chain."""
        return self._delegate('getAvgNodeChainTable')

    def getAvgNodeChain(self, *args):
        """Get per-node, per-chain average metrics."""
        return self._delegate('getAvgNodeChain', *args)

    def getAvgCacheTable(self) -> pd.DataFrame:
        """Get cache hit/miss table."""
        return self._delegate('getAvgCacheTable')

    def getAvgItemTable(self) -> pd.DataFrame:
        """Get per-item cache table."""
        return self._delegate('getAvgItemTable')

    def getAvgLossTable(self) -> pd.DataFrame:
        """Get per-station loss table."""
        return self._delegate('getAvgLossTable')

    def getAvgRegionLossTable(self) -> pd.DataFrame:
        """Get finite-capacity-region loss table."""
        return self._delegate('getAvgRegionLossTable')

    def getAvgOrbitTable(self) -> pd.DataFrame:
        """Get retrial-orbit table."""
        return self._delegate('getAvgOrbitTable')

    def getAvgOrbit(self):
        """Get retrial-orbit populations."""
        return self._delegate('getAvgOrbit')

    def getMomentTable(self, *args, **kwargs):
        """Get higher-moment table. The base implementation reads
        options.method, which SolverAUTOOptions does not carry, so this must
        delegate rather than inherit."""
        return self._delegate('getMomentTable', *args, **kwargs)

    def getMomentChainTable(self, *args, **kwargs):
        """Get higher-moment table by chain."""
        return self._delegate('getMomentChainTable', *args, **kwargs)

    def getMomentStationTable(self, *args, **kwargs):
        """Get higher-moment table by station."""
        return self._delegate('getMomentStationTable', *args, **kwargs)

    def getSensitivityTable(self, *args, **kwargs):
        """Get parameter-sensitivity table."""
        return self._delegate('getSensitivityTable', *args, **kwargs)

    def getMethodFeatureSet(self, method: str):
        """Get the feature set of a method on the delegated solver."""
        return self._delegate('getMethodFeatureSet', method)

    def isStochasticMethod(self, method: str) -> bool:
        """Report whether a method yields stochastic estimates."""
        return self._delegate('isStochasticMethod', method)

    def getFeatureSet(self):
        """Get the feature set of the solver selected for this model."""
        return self._delegate('getFeatureSet')

    def initFromSolver(self, initSolver):
        """Warm start from another solver's steady state."""
        return self._delegate('initFromSolver', initSolver)

    @staticmethod
    def defaultOptions() -> 'SolverAUTOOptions':
        """Default AUTO options."""
        return SolverAUTOOptions()

    # ========== Distribution Methods ==========

    def getCdfRespT(self):
        """Get CDF of response times at steady-state."""
        return self._delegate('getCdfRespT')

    def getCdfPassT(self):
        """Get CDF of passage times at steady-state."""
        return self._delegate('getCdfPassT')

    def getPerctRespT(self, percentiles: List[float]):
        """Get response time percentiles."""
        return self._delegate('getPerctRespT', percentiles)

    # ========== Transient Methods ==========

    def getTranAvg(self):
        """Get transient average metrics."""
        return self._delegate('getTranAvg')

    def getTranCdfRespT(self):
        """Get transient CDF of response times."""
        return self._delegate('getTranCdfRespT')

    def getTranCdfPassT(self):
        """Get transient CDF of passage times."""
        return self._delegate('getTranCdfPassT')

    # ========== Probability Methods ==========

    def getProb(self, node, state):
        """Get marginal state probability for a node."""
        return self._delegate('getProb', node, state)

    def getProbAggr(self, node, state_a):
        """Get aggregated state probability for a node."""
        return self._delegate('getProbAggr', node, state_a)

    def getProbSys(self):
        """Get joint system state probability."""
        return self._delegate('getProbSys')

    def getProbSysAggr(self):
        """Get aggregated system state probability."""
        return self._delegate('getProbSysAggr')

    def getProbMarg(self, node, jobclass, state_m):
        """Get marginalized state probability for station and class."""
        return self._delegate('getProbMarg', node, jobclass, state_m)

    def getProbNormConstAggr(self):
        """Get normalizing constant (log)."""
        return self._delegate('getProbNormConstAggr')

    # ========== Transient Probability Methods ==========

    def getTranProb(self, node):
        """Get transient state probabilities for a node."""
        return self._delegate('getTranProb', node)

    def getTranProbAggr(self, node):
        """Get transient aggregated probabilities for a node."""
        return self._delegate('getTranProbAggr', node)

    def getTranProbSys(self):
        """Get transient system probabilities."""
        return self._delegate('getTranProbSys')

    def getTranProbSysAggr(self):
        """Get transient system aggregated probabilities."""
        return self._delegate('getTranProbSysAggr')

    # ========== Sampling Methods ==========

    def sample(self, node, num_events: int = 1000):
        """Sample node states."""
        return self._delegate('sample', node, num_events)

    def sampleAggr(self, node, num_events: int = 1000):
        """Sample aggregated node states."""
        return self._delegate('sampleAggr', node, num_events)

    def sampleSys(self, num_events: int = 1000):
        """Sample system states."""
        return self._delegate('sampleSys', num_events)

    def sampleSysAggr(self, num_events: int = 1000):
        """Sample system aggregated states."""
        return self._delegate('sampleSysAggr', num_events)

    @classmethod
    def load(cls, *args, **kwargs) -> Any:
        """
        Factory method to load models or instantiate solvers.

        Usage 1: Load a model from file.
            model = LINE.load(filename)
            model = LINE.load(filename, verbose=True)
        Supported formats: .jsim, .jsimg, .jsimw (JMT), .xml/.lqn/.lqnx (LQN),
                          .mat (MATLAB), .pkl/.pickle (Python pickle)

        Usage 2: Instantiate a solver with a specific method.
            solver = LINE.load(method, model, **options)

        Args:
            For file loading:
                filename: Path to model file
                verbose: Print loading info (default False)
            For solver loading:
                method: Solver method name (e.g., 'mva', 'ctmc', 'ssa', 'fluid', 'nc', 'mam', 'jmt')
                model: Network model to solve
                **options: Solver options

        Returns:
            Loaded model object or Solver instance
        """
        import os

        if len(args) == 0:
            raise ValueError("LINE.load requires at least one argument")

        first_arg = args[0]

        # Check if this is file loading (first arg is a filename string with extension)
        if isinstance(first_arg, str) and '.' in first_arg:
            ext = os.path.splitext(first_arg)[1].lower()
            supported_extensions = {'.jsim', '.jsimg', '.jsimw', '.jmva',
                                   '.xml', '.lqn', '.lqnx', '.mat', '.pkl', '.pickle'}

            if ext in supported_extensions:
                filename = first_arg
                verbose = args[1] if len(args) > 1 else kwargs.get('verbose', False)

                if not os.path.exists(filename):
                    raise FileNotFoundError(f"File not found: {filename}")

                # Load based on format
                if ext in ('.jsim', '.jsimg', '.jsimw'):
                    # Load JSIM/JMT model
                    from ...api.io import jsim2line
                    return jsim2line(filename)

                elif ext == '.jmva':
                    # Load JMVA model
                    from ...api.io import jmva2line
                    return jmva2line(filename)

                elif ext in ('.xml', '.lqn', '.lqnx'):
                    # Load LQN model
                    from ...layered import LayeredNetwork
                    return LayeredNetwork.parse_xml(filename, verbose)

                elif ext == '.mat':
                    # Load MATLAB saved model
                    try:
                        from scipy.io import loadmat
                        data = loadmat(filename)
                        # Try to find a Network object
                        if 'model' in data:
                            result = data['model']
                        elif 'network' in data:
                            result = data['network']
                        else:
                            # Return the raw data
                            result = data
                        if verbose:
                            print(f"Loaded model from {filename}")
                        return result
                    except ImportError:
                        raise ImportError("scipy is required to load .mat files")

                elif ext in ('.pkl', '.pickle'):
                    import pickle
                    with open(filename, 'rb') as f:
                        result = pickle.load(f)
                    if verbose:
                        print(f"Loaded from {filename}")
                    return result

        # Otherwise, this is solver loading: LINE.load(method, model, **options)
        if len(args) < 2:
            raise ValueError("Solver loading requires method and model arguments: LINE.load(method, model)")

        method = args[0]
        model = args[1]
        options = kwargs

        # Parse the method and create appropriate solver
        method_lower = method.lower()

        # Strip solver prefix if present (e.g., 'mva.exact' -> 'exact')
        if method_lower in ('default', 'auto', 'line'):
            return cls(model, **options)

        elif method_lower.startswith('ctmc'):
            from ..solver_ctmc import SolverCTMC
            sub_method = method_lower.replace('ctmc.', '').replace('ctmc', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverCTMC(model, **options)

        elif method_lower.startswith('mva') or method_lower in (
            'amva', 'qna', 'sqrt', 'mm1', 'mmk', 'mg1', 'gm1', 'gig1', 'gim1',
            'gigk', 'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower',
            'gb.upper', 'gb.lower', 'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower'
        ):
            from ..solver_mva import SolverMVA
            sub_method = method_lower.replace('mva.', '').replace('mva', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverMVA(model, **options)

        elif method_lower.startswith('ssa') or method_lower in ('nrm', 'serial', 'parallel'):
            from ..solver_ssa import SolverSSA
            sub_method = method_lower.replace('ssa.', '').replace('ssa', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverSSA(model, **options)

        elif method_lower.startswith('jmt') or method_lower in ('jsim', 'jmva'):
            from ..wrappers.solver_jmt import SolverJMT
            sub_method = method_lower.replace('jmt.', '').replace('jmt', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverJMT(model, **options)

        elif method_lower.startswith('fluid') or method_lower.startswith('fld'):
            from ..solver_fld import SolverFLD
            sub_method = (method_lower.replace('fluid.', '').replace('fld.', '')
                          .replace('fluid', 'default').replace('fld', 'default'))
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverFLD(model, **options)

        elif method_lower.startswith('nc') or method_lower in (
            'comom', 'comomld', 'cub', 'ls', 'le', 'ble', 'mmint2', 'pana'
        ):
            from ..solver_nc import SolverNC
            sub_method = method_lower.replace('nc.', '').replace('nc', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverNC(model, **options)

        elif method_lower.startswith('mam'):
            from ..solver_mam import SolverMAM
            sub_method = method_lower.replace('mam.', '').replace('mam', 'default')
            if sub_method and sub_method != 'default':
                options['method'] = sub_method
            return SolverMAM(model, **options)

        elif method_lower.startswith('des'):
            from ..wrappers.solver_ldes import SolverLDES
            return SolverLDES(model, **options)

        else:
            # Unknown method - use SolverAUTO with the method as selection
            options['method'] = method
            return cls(model, **options)

    # ========== Python-style Aliases (snake_case) ==========

    # Aliases (Python naming convention)
    run_analyzer = runAnalyzer
    get_name = getName
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    get_selected_solver_name = getSelectedSolverName
    get_candidate_solver_names = getCandidateSolverNames
    set_forced_solver = setForcedSolver
    set_selection_method = setSelectionMethod

    # Basic metrics
    get_avg = getAvg
    get_avg_qlen = getAvgQLen
    get_avg_util = getAvgUtil
    get_avg_respt = getAvgRespT
    get_avg_resid_t = getAvgResidT
    get_avg_wait_t = getAvgWaitT
    get_avg_tput = getAvgTput
    get_avg_arv_r = getAvgArvR

    # Chain metrics
    get_avg_chain = getAvgChain
    get_avg_qlen_chain = getAvgQLenChain
    get_avg_util_chain = getAvgUtilChain
    get_avg_respt_chain = getAvgRespTChain
    get_avg_resid_t_chain = getAvgResidTChain
    get_avg_tput_chain = getAvgTputChain
    get_avg_arv_r_chain = getAvgArvRChain

    # Node metrics
    get_avg_node = getAvgNode
    get_avg_node_qlen_chain = getAvgNodeQLenChain
    get_avg_node_util_chain = getAvgNodeUtilChain
    get_avg_node_resp_t_chain = getAvgNodeRespTChain
    get_avg_node_resid_t_chain = getAvgNodeResidTChain
    get_avg_node_tput_chain = getAvgNodeTputChain
    get_avg_node_arv_r_chain = getAvgNodeArvRChain

    # System metrics
    get_avg_sys = getAvgSys
    get_avg_sys_resp_t = getAvgSysRespT
    get_avg_sys_tput = getAvgSysTput

    # Table methods
    get_avg_chain_table = getAvgChainTable
    get_avg_sys_table = getAvgSysTable
    get_avg_node_table = getAvgNodeTable
    get_avg_qlen_table = getAvgQLenTable
    get_avg_util_table = getAvgUtilTable
    get_avg_respt_table = getAvgRespTTable
    get_avg_tput_table = getAvgTputTable

    # Distribution methods
    get_cdf_resp_t = getCdfRespT
    get_cdf_pass_t = getCdfPassT
    get_perct_resp_t = getPerctRespT

    # Transient methods
    get_tran_avg = getTranAvg
    get_tran_cdf_resp_t = getTranCdfRespT
    get_tran_cdf_pass_t = getTranCdfPassT

    # Probability methods
    get_prob = getProb
    get_prob_aggr = getProbAggr
    get_prob_sys = getProbSys
    get_prob_sys_aggr = getProbSysAggr
    get_prob_marg = getProbMarg
    get_prob_norm_const_aggr = getProbNormConstAggr

    # Transient probability methods
    get_tran_prob = getTranProb
    get_tran_prob_aggr = getTranProbAggr
    get_tran_prob_sys = getTranProbSys
    get_tran_prob_sys_aggr = getTranProbSysAggr

    # Sampling methods
    sample_aggr = sampleAggr
    sample_sys = sampleSys
    sample_sys_aggr = sampleSysAggr

    # Snake_case forms of the accessors added above
    default_options = defaultOptions
    avg_node_table = getAvgNodeTable
    avg_node_chain_table = getAvgNodeChainTable
    get_avg_node_chain_table = getAvgNodeChainTable
    get_avg_node_chain = getAvgNodeChain
    get_avg_cache_table = getAvgCacheTable
    get_avg_item_table = getAvgItemTable
    get_avg_loss_table = getAvgLossTable
    get_avg_region_loss_table = getAvgRegionLossTable
    get_avg_orbit_table = getAvgOrbitTable
    get_avg_orbit = getAvgOrbit
    get_moment_table = getMomentTable
    get_moment_chain_table = getMomentChainTable
    get_moment_station_table = getMomentStationTable
    get_sensitivity_table = getSensitivityTable
    get_method_feature_set = getMethodFeatureSet
    is_stochastic_method = isStochasticMethod
    get_feature_set = getFeatureSet
    init_from_solver = initFromSolver

    # Capitalized forms every native solver exposes
    GetAvg = getAvg
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgChain = getAvgChain
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain
    GetAvgNode = getAvgNode
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain
    GetAvgSys = getAvgSys
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput
    GetAvgChainTable = getAvgChainTable
    GetAvgSysTable = getAvgSysTable
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    GetTranAvg = getTranAvg


# Alias for convenience
LINE = SolverAUTO


__all__ = [
    'SolverAUTO',
    'SolverAUTOOptions',
    'ModelAnalyzer',
    'SolverType',
    'LINE',
]
