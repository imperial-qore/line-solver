"""
Base class for native LINE solvers.

Provides common functionality shared across all native solver implementations.
The class hierarchy mirrors the JAR implementation:
- Solver: Base class for all solvers
- NetworkSolver(Solver): For single-network solvers (MVA, NC, CTMC, SSA, etc.)
- EnsembleSolver(Solver): For multi-model solvers (LN, UQ, ENV)
"""

from typing import Dict, List, Optional
from ..api.io.logging import line_warning


def method_label(requested, resolved):
    """Banner label for the solution method, harmonized across codebases.

    When method 'default' resolves at runtime to a concrete algorithm, the
    label is 'default/<resolved>' (e.g. 'default/mm1'), matching the MATLAB
    and JAR banner convention. Otherwise the resolved name is returned.
    """
    if requested == 'default' and resolved and resolved != 'default' \
            and not str(resolved).startswith('default/'):
        return 'default/' + str(resolved)
    return resolved if resolved else requested


class SolverFeatureSet:
    """
    A class to specify and check features supported by a solver.

    Mirrors MATLAB's SolverFeatureSet class. Used to track which model features
    are supported by each solver and to check compatibility.

    Attributes:
        list: Dict mapping feature names to boolean support status.
    """

    # Canonical feature-name registry; must stay identical to
    # jar FeatureSet.java and MATLAB SolverFeatureSet.fields
    FIELDS = [
        'ClassSwitch',
        'Cache',
        'Delay',
        'DelayStation',
        'Fork',
        'Join',
        'Logger',
        'Place',
        'QueueingPlace',
        'Queue',
        'JobSink',
        'Sink',
        'Source',
        'Router',
        'Transition',
        'Coxian',
        'Cox2',
        'APH',
        'Det',
        'Disabled',
        'Erlang',
        'Exp',
        'Gamma',
        'HyperExp',
        'Immediate',
        'Lognormal',
        'MAP',
        'DMAP',
        'MMAP',
        'BMAP',
        'MMPP2',
        'NHPP',
        'EmpiricalCdf',
        'Expolynomial',
        'Normal',
        'Pareto',
        'PH',
        'ME',
        'RAP',
        'Replayer',
        'Trace',
        'Uniform',
        'Weibull',
        'Bernoulli',
        'Binomial',
        'Geometric',
        'Poisson',
        # Distributions that no solver supports as a service or arrival
        # process, but whose name get_used_lang_features can still emit. They
        # must be registered: an unregistered name is dropped by set_true, so
        # the model would pass the support gate as if the distribution were
        # absent. No solver declares them, so a model using one is rejected.
        'DiscreteSampler',
        'DiscreteUniform',
        'Empirical',
        'GMM',
        'MMDP',
        'MMDP2',
        'MMPP',
        'MultivariateNormal',
        'NegBinomial',
        'Prior',
        'Zipf',
        'StatelessClassSwitcher',
        'CacheClassSwitcher',
        'CacheRetrieval',
        'InfiniteServer',
        'Forker',
        'Joiner',
        'LogTunnel',
        'SharedServer',
        'Buffer',
        'Region',
        'Linkage',
        'Enabling',
        'Inhibiting',
        'Timing',
        'Firing',
        'Storage',
        'RandomSource',
        'Dispatcher',
        'Server',
        'ServiceTunnel',
        'RoutingStrategy_PROB',
        'RoutingStrategy_RAND',
        'RoutingStrategy_RROBIN',
        'RoutingStrategy_WRROBIN',
        'RoutingStrategy_JSQ',
        'RoutingStrategy_SQ',
        'RoutingStrategy_RL',
        'SchedStrategy_INF',
        'SchedStrategy_FCFS',
        'SchedStrategy_FCFSPR',
        'SchedStrategy_FCFSPI',
        'SchedStrategy_FCFSPRIO',
        'SchedStrategy_FCFSPRPRIO',
        'SchedStrategy_FCFSPIPRIO',
        'SchedStrategy_LCFS',
        'SchedStrategy_LCFSPR',
        'SchedStrategy_LCFSPI',
        'SchedStrategy_LCFSPRIO',
        'SchedStrategy_LCFSPRPRIO',
        'SchedStrategy_LCFSPIPRIO',
        'SchedStrategy_SEPT',
        'SchedStrategy_LEPT',
        'SchedStrategy_SJF',
        'SchedStrategy_LJF',
        'SchedStrategy_SRPT',
        'SchedStrategy_SRPTPRIO',
        'SchedStrategy_PSJF',
        'SchedStrategy_FB',
        'SchedStrategy_LRPT',
        'SchedStrategy_SETF',
        'SchedStrategy_FSP',
        'SchedStrategy_PAS',
        'SchedStrategy_OI',
        'SchedStrategy_PS',
        'SchedStrategy_DPS',
        'SchedStrategy_GPS',
        'SchedStrategy_PSPRIO',
        'SchedStrategy_DPSPRIO',
        'SchedStrategy_GPSPRIO',
        'SchedStrategy_SIRO',
        'SchedStrategy_HOL',
        'SchedStrategy_EXT',
        'SchedStrategy_POLLING',
        'SchedStrategy_EDD',
        'SchedStrategy_EDF',
        'SchedStrategy_LPS',
        'ReplacementStrategy_RR',
        'ReplacementStrategy_FIFO',
        'ReplacementStrategy_SFIFO',
        'ReplacementStrategy_LRU',
        'ReplacementStrategy_HLRU',
        'ReplacementStrategy_CLIMB',
        'ReplacementStrategy_QLRU',
        'ClosedClass',
        'OpenClass',
        'SelfLoopingClass',
        'OpenSignal',
        'ClosedSignal',
        'SignalType_NEGATIVE',
        'SignalType_REPLY',
        'SignalType_CATASTROPHE',
        'SignalBatchRemoval',
        'SignalRemovalPolicy',
        'BatchArrival',
        'LoadDependence',
        'ClassDependence',
        'JointDependence',
        'SetupDelayOff',
        'Retrial',
        'Balking',
        'Reneging',
        'Breakdown',
        'Host',
        'Processor',
        'Task',
        'Entry',
        'Activity',
        'SyncCall',
        'AsyncCall',
        'ActivityPrecedence_PRE_SEQ',
        'ActivityPrecedence_POST_SEQ',
        'ActivityPrecedence_PRE_AND',
        'ActivityPrecedence_POST_AND',
        'ActivityPrecedence_PRE_OR',
        'ActivityPrecedence_POST_OR',
        'SchedStrategy_REF',
    ]

    def __init__(self):
        """Initialize with all features set to False."""
        self.list: Dict[str, bool] = {f: False for f in self.FIELDS}

    def set_true(self, features):
        """
        Set one or more features to True.

        Args:
            features: A single feature name (str) or list of feature names.

        Raises:
            ValueError: if a name is absent from FIELDS. Silently dropping it
                would make the model pass the solver support gate as if the
                feature were not used, so an unregistered name is a registry
                bug and must fail loudly, as it does in jar FeatureSet.setTrue.
        """
        if isinstance(features, str):
            features = [features]
        for feat in features:
            if feat not in self.list:
                raise ValueError(
                    "Unrecognized feature to set to true in the feature set: %s" % feat)
            self.list[feat] = True

    def set_false(self, features):
        """
        Set one or more features to False.

        Args:
            features: A single feature name (str) or list of feature names.

        Raises:
            ValueError: if a name is absent from FIELDS (see set_true).
        """
        if isinstance(features, str):
            features = [features]
        for feat in features:
            if feat not in self.list:
                raise ValueError(
                    "Unrecognized feature to set to false in the feature set: %s" % feat)
            self.list[feat] = False

    @staticmethod
    def supports_with_reason(feat_supported: 'SolverFeatureSet',
                             feat_used: 'SolverFeatureSet'):
        """
        Check if supported features cover all used features, returning a reason.

        Args:
            feat_supported: FeatureSet of features the solver supports.
            feat_used: FeatureSet of features the model uses.

        Returns:
            (bool, reason): bool is True when every used feature is supported;
            reason is a human-readable list of the offending feature names
            (empty string when bool is True). Unlike supports(), this does not
            emit a side-effect warning, so it is safe for method-aware gating
            and feature-driven method selection (which probe multiple methods).
        """
        unsupported = []
        for field in SolverFeatureSet.FIELDS:
            if feat_used.list.get(field, False) and not feat_supported.list.get(field, False):
                unsupported.append(field)
        if unsupported:
            feat_str = ', '.join(unsupported)
            reason = f'Some features are not supported by the chosen solver (feature: {feat_str}).'
            return False, reason
        return True, ''

    @staticmethod
    def supports(feat_supported: 'SolverFeatureSet', feat_used: 'SolverFeatureSet') -> bool:
        """
        Check if supported features cover all used features.

        Prints a warning if any used features are not supported.

        Args:
            feat_supported: FeatureSet of features the solver supports.
            feat_used: FeatureSet of features the model uses.

        Returns:
            True if all used features are supported, False otherwise.
        """
        ok, reason = SolverFeatureSet.supports_with_reason(feat_supported, feat_used)
        if not ok:
            line_warning('SolverFeatureSet', reason)
        return ok

    # MATLAB-compatible aliases
    setTrue = set_true
    setFalse = set_false


def supports_via_featureset(solver_cls, model) -> bool:
    """Gate a model's used language features against a solver's feature set.

    The shared implementation of MATLAB's Solver*.supports. Several Python
    solvers instead returned True unconditionally, or checked only the station
    and class counts, so they silently accepted every model regardless of the
    features it used; JMT even carried a "for now return True" placeholder.

    Struct-like inputs carrying no feature registry (a NetworkStruct rather than
    a Network) fall back to a structural sanity check, since there is nothing to
    gate against.

    Args:
        solver_cls: the solver class, which must expose getFeatureSet()
        model: the model to check

    Returns:
        True when every feature the model uses is supported by solver_cls.
    """
    get_used = (getattr(model, 'get_used_lang_features', None)
                or getattr(model, 'getUsedLangFeatures', None))
    if get_used is not None:
        feat_used = get_used()
        feat_supported = SolverFeatureSet()
        feat_supported.set_true(list(solver_cls.getFeatureSet()))
        return SolverFeatureSet.supports(feat_supported, feat_used)
    try:
        nstations = getattr(model, 'nstations', None)
        if nstations is None and hasattr(model, 'getNumberOfStations'):
            nstations = model.getNumberOfStations()
        nclasses = getattr(model, 'nclasses', None)
        if nclasses is None and hasattr(model, 'getNumberOfClasses'):
            nclasses = model.getNumberOfClasses()
        if nstations is None or nclasses is None:
            return False
        return nstations > 0 and nclasses > 0
    except Exception:
        return False


from .._aliasing import alias_getattr as _snake_camel_getattr  # noqa: E402


class Solver:
    """Base class for all LINE solvers.

    Provides common attributes and default behaviors for all solvers.

    Attributes:
        _table_silent (bool): If True, suppress automatic table printing in
            getAvgTable() and similar methods. Defaults to True so that
            tables are only printed when explicitly requested by the user.
    """

    # Suppress automatic table printing by default
    _table_silent = True

    # Resolve documented snake_case method names (get_tran_avg, sample_aggr,
    # list_valid_methods, get_prob_marg, ...) to their camelCase implementations.
    __getattr__ = _snake_camel_getattr

    def isStochastic(self):
        """True if the solver, with its currently configured method, returns
        stochastic estimates, i.e. results that depend on the random seed, as
        in simulation or Monte Carlo integration.

        A solver run with method 'default' may resolve the actual method only
        at runtime; once results are available the classification therefore
        uses the method recorded in the result (e.g. 'default/imci' when the
        NC default path resolved to Monte Carlo integration).
        """
        options = getattr(self, 'options', None)
        method = getattr(options, 'method', None) if options is not None else None
        if method is None:
            method = self.__dict__.get('method')
        result = getattr(self, '_result', None)
        if result is not None:
            if isinstance(result, dict):
                resolved = result.get('method')
            else:
                resolved = getattr(result, 'method', None)
            if resolved:
                method = resolved
        return self.isStochasticMethod(method)

    def isStochasticMethod(self, method):
        """Classify a (possibly runtime-resolved) method name of this solver
        as stochastic. Deterministic by default; subclasses with
        simulation-based or sampling-based methods override this.
        """
        return False

    is_stochastic = isStochastic
    is_stochastic_method = isStochasticMethod


def _validate_moment_order(order, maxorder=3):
    """Normalize a moment-order SET, mirroring MATLAB's validateMomentOrder.

    ``order`` is a set of moment orders. A scalar ``k`` is shorthand for
    ``1:k``, so that ``getMomentTable(2)`` means "up to the second moment" and
    not "the second moment alone"; a vector of two or more entries is taken
    literally.

    Returns a sorted tuple of ints. Raises ValueError outside 1..maxorder.

    Consequence of MATLAB's ``isscalar``: a one-element vector IS a scalar, so
    ``[2]`` takes the ``1:k`` path and yields ``(1, 2)``. "The second moment
    alone" is therefore not expressible, which is deliberate: a variance with no
    mean beside it is not a useful table, and ``[2, 3]`` remains available for
    the higher orders.

    Non-integers are rejected on BOTH paths. Rounding them silently would accept
    ``[1, 2.5]`` as ``[1, 3]``, i.e. answer a question that was not asked.
    """
    import numpy as np
    msg = ("order must be an integer in 1..%d, or a vector of such integers."
           % maxorder)
    try:
        arr = np.asarray(order, dtype=float).ravel()
    except (TypeError, ValueError):
        raise ValueError(msg)
    if arr.size == 0 or not np.all(np.isfinite(arr)):
        raise ValueError(msg)
    if (np.any(arr != np.round(arr)) or np.any(arr < 1)
            or np.any(arr > maxorder)):
        raise ValueError(msg)
    if arr.size == 1:
        return tuple(range(1, int(arr[0]) + 1))
    return tuple(sorted(set(int(v) for v in arr)))


def _is_linearizer_method(method):
    """True for the Linearizer family of solver methods, mirroring MATLAB's
    isLinearizerMethod.

    Everything else maps to the exact recursion: the moment analysis has only
    these two algorithms, and the exact one is the right default for a method
    with no approximate counterpart. The solver has already rejected methods it
    does not support, so no further validation belongs here.

    ``options.method`` is a plain string token in the native Python solvers
    ('default' when unset), so the comparison is a case-insensitive string
    match; there is no enum here to compare across classes.
    """
    return str(method).lower() in ('lin', 'amva.lin', 'egflin', 'gflin')


def _is_exact_mva_method(method):
    """Methods whose means are the exact MVA recursion, so pfqn_sens_mom's
    analytic derivatives apply directly. Mirrors MATLAB's isExactMvaMethod."""
    return str(method).lower() in ('default', 'mva', 'exact')


def _not_product_form_error(caller):
    """The message both moment tables raise on a non-product-form model.

    The moments rest on the product-form identity Cov = L dQ/dL, which is a
    theorem ABOUT the product form: it follows from a d/da log G = E[n], which
    needs the state distribution to be exponential-family in the demands.
    Outside product form, L dQ/dL remains computable and is simply NOT a
    covariance, so differentiating a non-product-form solver would return a
    confident wrong number. Refuse rather than do that; exact moments of such a
    model come from SolverCTMC's stationary distribution, or from LDES via
    setReward.
    """
    return ValueError(
        "%s requires a product-form model: the moment identity "
        "Cov[n,n] = L dQ/dL holds only under product form, so no correct value "
        "exists here. Use SolverCTMC (exact distribution) or SolverLDES with "
        "setReward for the moments of a non-product-form model." % (caller,))


class _FiniteDifferenceMom:
    """The moment struct produced by the numerical-derivative path.

    Mirrors the fields of MATLAB's packFiniteDifference / packChainFiniteDifference
    return struct: m, dm, d2m, Cov, CovAsym, Var, M2, M3, Skew. It deliberately
    does NOT carry pfqn_sens_mom's X/Q/U/R, which the finite-difference path never
    computes; a field present but meaningless would be worse than absent.
    """

    def __init__(self, m, dm, d2m, Var, Cov, M2, M3, Skew, CovAsym):
        self.m = m
        self.dm = dm
        self.d2m = d2m
        self.Var = Var
        self.Cov = Cov
        self.M2 = M2
        self.M3 = M3
        self.Skew = Skew
        self.CovAsym = CovAsym


def _scale_demand_column(sn, station, factor):
    """Scale one station's whole demand column by FACTOR, via its service rates.

    Because D(i,r) = visits(i,r)/rate(i,r), dividing station i's rates by FACTOR
    multiplies its whole demand column by FACTOR.

    Unlike MATLAB, where NetworkStruct is a plain struct and assignment copies,
    the native ``NetworkStruct`` is a mutable object whose numpy fields are held
    BY REFERENCE. The struct is therefore deep-copied before the rates are
    touched: perturbing the caller's sn in place would silently corrupt every
    later solve off the same model, whose getStruct() hands back a cached struct.
    """
    sn2 = sn.copy()
    sn2.rates[station, :] = sn.rates[station, :] / factor
    return sn2


def _scale_group_demands(sn, station, groups, g, factor):
    """Scale the demands of group ``g``'s classes at one station by FACTOR, via
    the rates. See _scale_demand_column on the deep copy."""
    import numpy as np
    sn2 = sn.copy()
    cls = np.nonzero(groups == g)[0]
    sn2.rates[station, cls] = sn.rates[station, cls] / factor
    return sn2


def _solve_totals(self, sn, qst):
    """Station totals under the solver's OWN method, whatever that is."""
    QN = self.solveMeansForStruct(sn)
    return QN[qst, :].sum(axis=1)


def _solve_group_totals(self, sn, qst, groups, Cg):
    """Per-(station,group) mean queue lengths under the solver's own method."""
    import numpy as np
    QN = self.solveMeansForStruct(sn)
    M = len(qst)
    mg = np.zeros((M, Cg))
    for i in range(M):
        for g in range(Cg):
            mg[i, g] = QN[qst[i], groups == g + 1].sum()
    return mg


def _moments_by_finite_difference(self, sn, queue_nodes):
    """Moments of the per-station totals from ANY solver and method, by central
    differences of that method's OWN mean queue lengths. Mirrors MATLAB's
    momentsByFiniteDifference.

    The identity does not care HOW the mean queue lengths were obtained, so the
    solver itself is used as a mean-value oracle and differentiated numerically.
    This is exactly the move the Linearizer makes analytically (Strelen Sec. 5:
    differentiate the approximate fixed point, then apply the exact identity),
    generalized to any product-form method. The parameter is y_i, a scaling of
    station i's whole demand column, which is Strelen's x_i; at y = 1 the
    y-derivatives are the scaled x-derivatives that (3.2) asks for.

    The oracle re-runs THIS solver (see solveMeansForStruct), so every method of
    every product-form solver is covered, including the normalizing-constant
    methods of SolverNC (comom, ca, le, ...) and the summation methods of
    SolverMVA (sum, esum), which no hand-differentiated implementation reaches.

    Cost: 2*M extra solves. Accuracy: the moments inherit the accuracy of the
    method's means, and the second derivative inherits the usual h^2 truncation.
    """
    import numpy as np
    M = len(queue_nodes)
    h = 1e-4
    qst = np.asarray(sn.nodeToStation, dtype=int).flatten()[queue_nodes]
    m0 = _solve_totals(self, sn, qst)
    dm = np.zeros((M, M))
    d2m = np.zeros(M)
    for hcol in range(M):
        mp = _solve_totals(self, _scale_demand_column(sn, qst[hcol], 1 + h), qst)
        mm = _solve_totals(self, _scale_demand_column(sn, qst[hcol], 1 - h), qst)
        dm[:, hcol] = (mp - mm) / (2 * h)
        d2m[hcol] = (mp[hcol] - 2 * m0[hcol] + mm[hcol]) / h ** 2
    return _pack_finite_difference(m0, dm, d2m)


def _pack_finite_difference(m, dm, d2m):
    """(3.2), applied to numerically obtained derivatives. Mirrors MATLAB's
    packFiniteDifference."""
    import numpy as np
    M = len(m)
    CovAsym = float(np.max(np.abs(dm - dm.T))) if M > 0 else 0.0
    Cov = (dm + dm.T) / 2
    Var = np.zeros(M)
    M2 = np.zeros(M)
    M3 = np.zeros(M)
    Skew = np.zeros(M)
    for i in range(M):
        Var[i] = dm[i, i]
        M2[i] = dm[i, i] + m[i] ** 2
        M3[i] = d2m[i] + (1 + 3 * m[i]) * dm[i, i] + m[i] ** 3
        mu3 = M3[i] - 3 * m[i] * M2[i] + 2 * m[i] ** 3
        Skew[i] = mu3 / Var[i] ** 1.5 if Var[i] > 0 else float('nan')
    return _FiniteDifferenceMom(m, dm, d2m, Var, Cov, M2, M3, Skew, CovAsym)


def _chain_moments_by_finite_difference(self, sn, queue_nodes, groups, Cg):
    """Per-chain moments from ANY solver and method, by central differences.
    Mirrors MATLAB's chainMomentsByFiniteDifference.

    The parameter y_(i,g) scales the demands of group g's classes at station i,
    which is the class-subset parameter T of Akyildiz-Strelen Theorem 1; the
    moments it generates are those of Q_(i,g) = sum_(r in g) n(i,r). Since
    D(i,r) = visits(i,r)/rate(i,r), the perturbation is applied to the rates of
    that group's classes at that station alone, leaving the other classes'
    demands at the same station untouched. That per-class granularity is what
    separates this from the station table, which scales a whole column.

    Cost: 2*M*Cg extra solves.
    """
    import numpy as np
    M = len(queue_nodes)
    h = 1e-4
    qst = np.asarray(sn.nodeToStation, dtype=int).flatten()[queue_nodes]
    m0 = _solve_group_totals(self, sn, qst, groups, Cg)
    dm = np.zeros((M, Cg, M, Cg))
    d2m = np.zeros((M, Cg))
    for hi in range(M):
        for hg in range(Cg):
            snp = _scale_group_demands(sn, qst[hi], groups, hg + 1, 1 + h)
            snm = _scale_group_demands(sn, qst[hi], groups, hg + 1, 1 - h)
            mp = _solve_group_totals(self, snp, qst, groups, Cg)
            mm = _solve_group_totals(self, snm, qst, groups, Cg)
            dm[:, :, hi, hg] = (mp - mm) / (2 * h)
            d2m[hi, hg] = (mp[hi, hg] - 2 * m0[hi, hg] + mm[hi, hg]) / h ** 2
    return _pack_chain_finite_difference(m0, dm, d2m, Cg)


def _pack_chain_finite_difference(m, dm, d2m, Cg):
    """(3.2), applied to numerically obtained derivatives, per (station,group).
    Mirrors MATLAB's packChainFiniteDifference."""
    import numpy as np
    M = m.shape[0]
    flat = dm.reshape(M * Cg, M * Cg)
    CovAsym = float(np.max(np.abs(flat - flat.T))) if flat.size else 0.0
    flat = (flat + flat.T) / 2
    # a single group carries no group index, so collapse it, exactly as
    # pfqn_sens_mom does in the same case
    if Cg == 1:
        Cov = flat.reshape(M, M)
        dm_out = dm.reshape(M, M)
    else:
        Cov = flat.reshape(M, Cg, M, Cg)
        dm_out = dm
    Var = np.zeros((M, Cg))
    M2 = np.zeros((M, Cg))
    M3 = np.zeros((M, Cg))
    Skew = np.zeros((M, Cg))
    for i in range(M):
        for g in range(Cg):
            d1 = dm[i, g, i, g]
            Var[i, g] = d1
            M2[i, g] = d1 + m[i, g] ** 2
            M3[i, g] = d2m[i, g] + (1 + 3 * m[i, g]) * d1 + m[i, g] ** 3
            mu3 = M3[i, g] - 3 * m[i, g] * M2[i, g] + 2 * m[i, g] ** 3
            Skew[i, g] = mu3 / Var[i, g] ** 1.5 if Var[i, g] > 0 \
                else float('nan')
    return _FiniteDifferenceMom(m, dm_out, d2m, Var, Cov, M2, M3, Skew, CovAsym)


def _moment_scv(variance, meanv):
    """Squared coefficient of variation. NaN when the mean is zero, since the
    SCV is then undefined rather than infinite in any useful sense."""
    import numpy as np
    if np.isnan(variance) or meanv <= 0:
        return float('nan')
    return float(variance / meanv ** 2)


def _sched_name_of(sn, s_idx):
    """Scheduling strategy of station ``s_idx``, as an upper-case name.

    ``sn.sched`` may be populated with ``line_solver.lang.base.SchedStrategy``,
    ``line_solver.constants.SchedStrategy`` or a bare int depending on the path
    that built the struct. Those enums carry identical numeric values but are
    distinct classes, and ``constants.SchedStrategy`` is a plain ``Enum`` rather
    than an ``IntEnum``, so ``==`` across them is silently False. Compare by NAME
    (see CLAUDE.md, the ProcessType/DropStrategy note).
    """
    sc = sn.sched.get(s_idx) if isinstance(sn.sched, dict) else None
    if sc is None:
        return ''
    name = getattr(sc, 'name', None)
    if name is None:
        from ..constants import SchedStrategy
        try:
            name = SchedStrategy(int(sc)).name
        except (ValueError, TypeError):
            return ''
    return str(name).upper()


def _sched_is_fcfs(sn, s_idx):
    """True when station ``s_idx`` is scheduled FCFS."""
    return _sched_name_of(sn, s_idx) == 'FCFS'


def _sched_is_ps(sn, s_idx):
    """True when station ``s_idx`` is scheduled PS."""
    return _sched_name_of(sn, s_idx) == 'PS'


def _rates_are_exponential(sn, s_idx, classes):
    """True when every class in ``classes`` has exponential service at station
    ``s_idx``. The PS sojourn-time moments of Mitra and Morrison assume it; a
    phase-type service leaves the mean intact but not the moments. ProcessType
    is compared by NAME because its numeric values differ across codebases.
    """
    import numpy as np
    procid = getattr(sn, 'procid', None)
    if procid is None:
        return False
    for r in classes:
        pt = np.asarray(procid, dtype=object)[s_idx, r]
        name = getattr(pt, 'name', None)
        if name is None:
            from ..constants import ProcessType
            try:
                name = ProcessType(int(pt)).name
            except (ValueError, TypeError):
                return False
        if str(name).upper() != 'EXP':
            return False
    return True


def _station_is_feedback_free(sn, s_idx, R):
    """True when no job can return to station ``s_idx``.

    An open PS station sees Poisson arrivals only under this condition, which is
    Melamed's: the station must lie on no routing cycle. The test aggregates the
    classes, so a cycle closed through a class switch also disqualifies the
    station. The source and the sink emit no edges here: jobs never leave a
    sink, and the rt arc that leads back to the source is bookkeeping.
    """
    import numpy as np
    from ..api.sn.network_struct import NodeType
    M = int(sn.nstations)
    rt = np.asarray(sn.rt, dtype=float)
    station_to_node = np.asarray(sn.stationToNode).flatten()
    adj = np.zeros((M, M), dtype=bool)
    for i in range(M):
        nt = sn.nodetype[int(station_to_node[i])]
        if nt in (NodeType.SOURCE, NodeType.SINK):
            continue
        for j in range(M):
            blk = rt[i * R:(i + 1) * R, j * R:(j + 1) * R]
            adj[i, j] = bool(np.any(blk > 0))
    reach = adj[s_idx, :].copy()
    frontier = np.flatnonzero(reach)
    while frontier.size > 0:
        nxt = np.any(adj[frontier, :], axis=0) & ~reach
        reach = reach | nxt
        frontier = np.flatnonzero(nxt)
    return not bool(reach[s_idx])


def _fcfs_rates_are_class_independent(sn, queue_nodes, R, node_to_station,
                                      rates):
    """An FCFS station in a BCMP network must serve every class at the same
    exponential rate; that is the precondition for reading a per-visit rate off
    it. Non-FCFS stations are unconstrained here, see the caller."""
    import numpy as np
    from ..constants import GlobalConstants
    for node in queue_nodes:
        s_idx = int(node_to_station[node])
        if not _sched_is_fcfs(sn, s_idx):
            continue
        ref = -1.0
        for r in range(R):
            rate = rates[s_idx, r]
            if not np.isfinite(rate) or rate <= 0:
                continue
            if ref < 0:
                ref = float(rate)
            elif abs(rate - ref) > GlobalConstants.FineTol * max(1.0, ref):
                return False
    return True


def _service_time_of(rates, s_idx, R):
    """The common service time of a station, i.e. the reciprocal of the
    class-independent rate. Zero if no class is served here."""
    import numpy as np
    for r in range(R):
        rate = rates[s_idx, r]
        if np.isfinite(rate) and rate > 0:
            return 1.0 / float(rate)
    return 0.0


class NetworkSolver(Solver):
    """Base class for single-network LINE solvers.

    Used by: SolverMVA, SolverNC, SolverCTMC, SolverSSA, SolverMAM, SolverJMT,
             SolverLDES, SolverQNS, SolverAuto, SolverFLD
    """

    def libraries(self):
        """Third-party libraries this solver will use, without printing.

        Attribution in LINE is pull-based, as in Sage: nothing is written to
        the console during a solve. Mirrors MATLAB @NetworkSolver/libraries.m.

        Returns:
            list of library names, possibly empty
        """
        from ..api.io.attribution import solver_libraries
        return solver_libraries(self)

    def citations(self, display=False):
        """Bibliographic references for the algorithms this solver used.

        The references follow what the run actually did: the method the
        analyzer resolved to (not merely the one requested), the fork-join
        transformation if the model has forks, and the percentile method of the
        last getPerctRespT call. Mirrors MATLAB @NetworkSolver/citations.m.

        Args:
            display: print the list instead of only returning it

        Returns:
            list of dicts with keys 'key' (internal bibliography key, never
            displayed), 'ref' and 'covers'
        """
        from ..api.io.citations import citations_for

        family = {
            'SolverCTMC': 'ctmc', 'SolverSSA': 'ssa', 'SolverFLD': 'fld',
            'SolverFluid': 'fld', 'SolverLN': 'ln', 'SolverJMT': 'jmt',
            'SolverLQNS': 'lqns', 'SolverENV': 'env', 'SolverMVA': 'mva',
            'SolverNC': 'nc', 'SolverMAM': 'mam', 'SolverBA': 'ba',
        }.get(type(self).__name__, '')

        tokens = []
        if family in ('ctmc', 'ssa', 'fld', 'ln', 'jmt', 'lqns', 'env'):
            tokens.append(family)

        def add_method(m):
            m = str(m or '').strip().lower()
            if not m or m == 'default':
                return
            tokens.append('%s.%s' % (family, m) if family else m)

        add_method(getattr(getattr(self, 'options', None), 'method', ''))
        reported = getattr(getattr(self, '_result', None), 'method', None)
        if reported:
            # the reported name is 'default/<actual>' when dispatch chose
            for part in str(reported).split('/'):
                add_method(part)

        model = getattr(self, 'model', None)
        if model is not None and getattr(model, 'has_fork', None) is not None:
            try:
                has_fork = model.has_fork() if callable(model.has_fork) else bool(model.has_fork)
            except Exception:
                has_fork = False
            if has_fork:
                fjm = str(getattr(getattr(self.options, 'config', None), 'fork_join', '') or 'mmt').lower()
                if fjm in ('default', 'fjt', ''):
                    fjm = 'mmt'
                elif fjm == 'heidelberger-trivedi':
                    fjm = 'ht'
                tokens.append(fjm)

        last_perct = getattr(self, '_last_perct_method', '')
        if last_perct:
            tokens.append(last_perct)

        entries = citations_for(tokens)
        if display:
            if not entries:
                print('No algorithm references recorded for this run.')
            for e in entries:
                print('%s\n    covers: %s' % (e['ref'], e['covers']))
        return entries

    def avg(self, *args):
        """Alias for getAvg (returns QN, UN, RN, TN, AN, WN)."""
        return self.getAvg(*args)

    def avg_table(self):
        """Get average performance metrics as an IndexedTable with proper formatting.

        This method wraps the raw DataFrame from getAvgTable() in an IndexedTable
        to provide MATLAB-style number formatting (e.g., 0 instead of 0.00000).

        Returns:
            IndexedTable: Wrapped DataFrame with MATLAB-style formatting.
        """
        import pandas as pd
        from ..indexed_table import IndexedTable
        result = self.getAvgTable()
        if result is None:
            return None
        # If already an IndexedTable, return as-is
        if isinstance(result, IndexedTable):
            return result
        # If a DataFrame, wrap it
        if isinstance(result, pd.DataFrame):
            return IndexedTable(result)
        return result

    get_avg_table = avg_table

    def _orbit_struct(self):
        """NetworkStruct backing this solver. Most solvers expose the model as
        self.model; SolverMAM keeps the already-refreshed struct on self.sn and
        the model on self.network, so both spellings are resolved here."""
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'getStruct'):
            return model.getStruct()
        sn = getattr(self, 'sn', None)
        if sn is not None:
            return sn
        network = getattr(self, 'network', None)
        if network is not None and hasattr(network, 'getStruct'):
            return network.getStruct()
        raise AttributeError('The solver exposes no network structure for getAvgOrbit.')

    def getAvgOrbit(self):
        """Mean number of jobs waiting in the ORBIT of each retrial station, as
        an (nstations, nclasses) array. Stations that are not retrial queues
        report 0.

        A retrial station has no waiting room: a job that finds every server
        busy joins the orbit instead of queueing, so its station population
        splits into the jobs currently in service and the jobs orbiting.
        getAvgQLen reports the whole station population, which is why the orbit
        had to be recovered by hand as QLen - Util. This method reports it
        directly.

        The in-service population is obtained from the station throughput by
        Little's law applied to the servers alone, E[in service] = X * E[S],
        which holds for any service distribution and any number of servers, so
        the orbit length is exact whenever QLen and Tput are.
        """
        import numpy as np
        avg = self.getAvg()
        QN = np.atleast_2d(np.asarray(avg[0], dtype=float))
        TN = np.atleast_2d(np.asarray(avg[3], dtype=float))
        sn = self._orbit_struct()
        ON = np.zeros_like(QN)

        rproc = getattr(sn, 'retrialProc', None)
        if rproc is None:
            return ON
        rates = np.atleast_2d(np.asarray(sn.rates, dtype=float))
        for ist in range(ON.shape[0]):
            if ist >= len(rproc):
                continue
            for r in range(ON.shape[1]):
                if r >= len(rproc[ist]) or rproc[ist][r] is None:
                    continue  # not a retrial station for this class: no orbit
                rate_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else np.nan
                if not np.isfinite(rate_ir) or rate_ir <= 0:
                    continue
                in_service = TN[ist, r] / rate_ir  # Little's law on the servers
                ON[ist, r] = max(0.0, QN[ist, r] - in_service)
        return ON

    get_avg_orbit = getAvgOrbit

    def getAvgOrbitTable(self):
        """Table of the mean orbit length of every retrial station-class pair,
        with the station population and the in-service population it decomposes
        into.

        Reported as a separate table rather than as an extra column of
        getAvgTable so that the average table keeps its shape for models
        without retrials.
        """
        import numpy as np
        import pandas as pd
        from ..indexed_table import IndexedTable
        ON = self.getAvgOrbit()
        avg = self.getAvg()
        QN = np.atleast_2d(np.asarray(avg[0], dtype=float))
        TN = np.atleast_2d(np.asarray(avg[3], dtype=float))
        sn = self._orbit_struct()
        rproc = getattr(sn, 'retrialProc', None)
        rates = np.atleast_2d(np.asarray(sn.rates, dtype=float))

        rows = []
        if rproc is not None:
            for ist in range(ON.shape[0]):
                if ist >= len(rproc):
                    continue
                for r in range(ON.shape[1]):
                    if r >= len(rproc[ist]) or rproc[ist][r] is None:
                        continue
                    rate_ir = rates[ist, r] if ist < rates.shape[0] and r < rates.shape[1] else np.nan
                    if np.isfinite(rate_ir) and rate_ir > 0:
                        in_service = TN[ist, r] / rate_ir
                    else:
                        in_service = 0.0
                    rows.append({
                        'Station': str(sn.nodenames[int(sn.stationToNode[ist])]),
                        'JobClass': str(sn.classnames[r]),
                        'QLen': QN[ist, r],
                        'InService': in_service,
                        'Orbit': ON[ist, r],
                    })
        df = pd.DataFrame(rows, columns=['Station', 'JobClass', 'QLen',
                                         'InService', 'Orbit'])
        return IndexedTable(df)

    get_avg_orbit_table = getAvgOrbitTable

    def getAvgLossTable(self):
        """Table of loss (drop) metrics for every station-class pair that
        receives offered traffic: offered arrival rate (ArvR), carried
        throughput (Tput), loss rate (ArvR - Tput, the rate of jobs dropped by
        finite capacity, blocking, or reneging) and loss ratio (LossRate /
        ArvR).

        Only pairs with ArvR > 0 are listed, which excludes the Source (whose
        offered arrival rate is zero); a lossless station has ArvR = Tput and
        so LossRate = LossRatio = 0.
        """
        import numpy as np
        import pandas as pd
        from ..indexed_table import IndexedTable
        # Derived from the average table's ArvR (offered) and Tput (carried)
        # columns, so the loss identity stays consistent with what getAvgTable
        # reports. The native result's raw arrival-rate array carries the
        # station-carried rate, not the offered rate, so it is not used here.
        avgTable = self.getAvgTable()
        data = getattr(avgTable, 'data', avgTable)

        # Fork-Join quorum sibling-drop rate (LDES only), station-indexed. At a
        # synchronizing Join the identity LossRate = ArvR - Tput does not hold
        # (Tput is in parent units, discarded siblings in sibling units), so on
        # Join rows the explicit drop rate replaces ArvR - Tput and the ratio
        # uses the offered sibling rate ArvR + drop. Look up by name since the
        # avg table carries station/class names, not indices.
        res = getattr(self, '_ldes_result', None) or getattr(self, '_result', None)
        dj = getattr(res, 'DropRateJoin', None) if res is not None else None
        drop_lookup = {}
        if dj is not None:
            dj = np.atleast_2d(np.asarray(dj, dtype=float))
            sn = getattr(self, '_sn', None)
            if sn is not None:
                nodenames = list(sn.nodenames)
                classnames = [str(c) for c in sn.classnames]
                stationToNode = np.asarray(sn.stationToNode).astype(int).flatten()
                for ist in range(dj.shape[0]):
                    if ist >= len(stationToNode):
                        continue
                    sname = nodenames[int(stationToNode[ist])]
                    for r in range(dj.shape[1]):
                        if r < len(classnames):
                            drop_lookup[(str(sname), classnames[r])] = float(dj[ist, r])

        rows = []
        for _, row in data.iterrows():
            a = float(row['ArvR'])
            if not np.isfinite(a) or a <= 0:
                continue
            t = float(row['Tput'])
            d = drop_lookup.get((str(row['Station']), str(row['JobClass'])), 0.0)
            if np.isfinite(d) and d > 0:
                lr = d
                lc = d / a
            else:
                lr = a - t
                lc = (a - t) / a
            rows.append({
                'Station': str(row['Station']),
                'JobClass': str(row['JobClass']),
                'ArvR': a,
                'Tput': t,
                'LossRate': lr,
                'LossRatio': lc,
            })
        df = pd.DataFrame(rows, columns=['Station', 'JobClass', 'ArvR', 'Tput',
                                         'LossRate', 'LossRatio'])
        return IndexedTable(df)

    get_avg_loss_table = getAvgLossTable

    def supportsExactSensitivity(self):
        """True when the solver evaluates a product-form recursion that
        getSensitivityTable can differentiate analytically. False here, so that
        a solver reaching this base implementation obtains its sensitivities by
        finite differences on its own predictions. Overridden by SolverMVA and
        SolverNC.
        """
        return False

    supports_exact_sensitivity = supportsExactSensitivity

    def getSensitivityTable(self, method='auto', step=None, scheme='forward'):
        """Performance sensitivities with respect to service rates.

        Returns a DataFrame with one row per (Station, JobClass) giving the
        derivative of that row's mean performance measures with respect to that
        station-class service RATE: dTput_dRate, dRespT_dRate, dQLen_dRate,
        dUtil_dRate.

        Two branches produce the derivatives, selected automatically:

        'exact'  Analytic differentiation of a product-form recursion, exact to
                 machine precision and cheaper than a single extra solve.
                 Closed networks use pfqn_sens (differentiated MVA), open
                 networks the closed-form BCMP sensitivities (their stations
                 decouple). Rate derivatives follow the chain rule
                 d(.)/d(rate) = -(L/rate) d(.)/dL, since
                 L(i,r) = visits(i,r)/rate(i,r). Available only on the solvers
                 that evaluate that recursion, SolverMVA and SolverNC, and only
                 for single-server queues plus an optional delay, with mixed
                 (open+closed) models excluded.

        'fd'     Forward or central finite differences on the CALLING solver's
                 own predictions: the service process at (station,class) is
                 rate-scaled by (1+h), the same solver with the same options is
                 re-run, and the difference quotient is formed. This costs
                 1+M*R solves (forward) or 2*M*R (central), and it is the only
                 branch that applies to non-product-form models, so it is what
                 every solver other than SolverMVA and SolverNC uses.

        Args:
            method: 'auto' (default: exact where available and in scope, finite
                differences otherwise), 'exact' or 'fd'.
            step: relative step of the rate perturbation, default 1e-4 for
                deterministic solvers and 1e-2 for the simulators, whose Monte
                Carlo error would otherwise dominate the difference quotient.
            scheme: 'forward' (default) or 'central'.

        The branch actually taken is reported in ``.attrs['method']``.

        Simulation solvers must be run with common random numbers for the
        difference quotient to be meaningful: the same options, and hence the
        same seed, are reused for the base and the perturbed runs. An unset
        seed is pinned before the sweep so that the runs remain paired.
        """
        import numpy as np
        import pandas as pd
        from ..api.sn.network_struct import NodeType

        if isinstance(method, str):
            method = method.lower()
        if isinstance(scheme, str):
            scheme = scheme.lower()
        if method not in ('auto', 'exact', 'fd'):
            raise ValueError("The method must be one of 'auto', 'exact', 'fd'.")
        if scheme not in ('forward', 'central'):
            raise ValueError("The scheme must be 'forward' or 'central'.")

        sn = self.model.getStruct()
        R = int(sn.nclasses)
        queue_nodes = [i for i, nt in enumerate(sn.nodetype)
                       if nt == NodeType.QUEUE]
        Mq = len(queue_nodes)

        exact_available = self.supportsExactSensitivity()
        in_scope, scope_msg = self._sensitivityExactScope(sn)
        if method == 'exact':
            if not exact_available:
                raise ValueError(
                    "Exact analytic sensitivities differentiate a product-form "
                    "recursion and are available on SolverMVA and SolverNC "
                    "only; %s must use 'fd'." % type(self).__name__)
            if not in_scope:
                raise ValueError(scope_msg)
            use_exact = True
        elif method == 'fd':
            use_exact = False
        else:
            use_exact = exact_available and in_scope

        if use_exact:
            dTput, dRespT, dQLen, dUtil, mask, sens = \
                self._sensitivityExact(sn, queue_nodes, R)
            method_used = 'exact'
        else:
            dTput, dRespT, dQLen, dUtil, mask = \
                self._sensitivityFD(sn, queue_nodes, R, step, scheme)
            method_used = 'fd'
            # The finite-difference branch differentiates the solver itself and
            # so has no Jacobian object to report, matching MATLAB (sens = [])
            # and the JAR (getSens() == null) on that branch.
            sens = None

        nodenames = list(sn.nodenames)
        classnames = list(sn.classnames)
        rows = []
        for ist, node in enumerate(queue_nodes):
            for r in range(R):
                if not mask[ist, r]:
                    continue
                rows.append({
                    'Station': str(nodenames[node]),
                    'JobClass': str(classnames[r]),
                    'dTput_dRate': float(dTput[ist, r]),
                    'dRespT_dRate': float(dRespT[ist, r]),
                    'dQLen_dRate': float(dQLen[ist, r]),
                    'dUtil_dRate': float(dUtil[ist, r]),
                })
        SensTable = pd.DataFrame(rows, columns=['Station', 'JobClass',
                                 'dTput_dRate', 'dRespT_dRate', 'dQLen_dRate',
                                 'dUtil_dRate'])
        SensTable.attrs['method'] = method_used
        # The pfqn_sens result of the closed exact branch, the counterpart of the
        # MATLAB second output and of NetworkSensitivityTable.getSens() in the
        # JAR. None on the open and the finite-difference branches.
        SensTable.attrs['sens'] = sens
        return SensTable

    get_sensitivity_table = getSensitivityTable

    @staticmethod
    def _sensitivityExactScope(sn):
        """Scope of the analytic branch: single-server stations (a delay, with
        infinite servers, is allowed) and not a mixed open-and-closed model.
        Class switching is supported: the branch aggregates classes into chains
        before differentiating.
        """
        import numpy as np
        from ..api.sn.transforms import sn_get_product_form_params
        try:
            _, _, _, _, _, S, _ = sn_get_product_form_params(sn)
        except Exception:
            return False, ("getSensitivityTable could not extract the "
                           "product-form parameters of this model.")
        S = np.asarray(S, dtype=float).flatten()
        if np.any(S[np.isfinite(S)] > 1):
            return False, ("getSensitivityTable supports single-server "
                           "stations only.")
        N = np.asarray(sn.njobs, dtype=float).flatten()
        if np.any(np.isinf(N)) and np.any(np.isfinite(N)):
            return False, ("getSensitivityTable does not yet support mixed "
                           "(open+closed) networks.")
        return True, ''

    def _sensitivityExact(self, sn, queue_nodes, R):
        """Analytic branch: differentiated MVA (closed) or closed-form BCMP
        (open), both evaluated at CHAIN level and then disaggregated back to the
        classes.

        A product-form model is solved per chain, not per class: a chain carries
        the population and the arrival rate, and a class is a share of it. With
        class switching the two differ, the whole chain population sitting on the
        reference class, so the recursion must see the chain demands
        Dc(i,c) = sum_{r in c} D(i,r), Nc(c) = sum_{r in c} N(r), Zc likewise.
        The parameter of the table is still the per-class rate mu(i,r), which
        enters exactly one chain demand, giving the chain rule
        d(.)/dmu(i,r) = -(D(i,r)/mu(i,r)) d(.)/dDc(i,c). The per-class measures
        are then composed from the chain ones, which is also what fixes the visit
        ratios: a class throughput at a station is X_c*v(i,r) and a per-visit
        response time is Q(i,r)/T(i,r), whereas the chain quantities are per
        chain and per visit-chain respectively.
        """
        import numpy as np
        from ..api.sn.transforms import sn_get_product_form_params
        from ..api.pfqn.sens import pfqn_sens

        lam, D, Np, Z, _, _, _ = sn_get_product_form_params(sn)
        D = np.atleast_2d(np.asarray(D, dtype=float))
        lam = np.asarray(lam, dtype=float).flatten()
        Np = np.asarray(Np, dtype=float).flatten()
        N = np.asarray(sn.njobs, dtype=float).flatten()
        is_open = bool(np.any(np.isinf(N)))
        Mq = len(queue_nodes)
        C = int(sn.nchains)
        chains = np.atleast_2d(np.asarray(sn.chains))
        node_to_station = np.asarray(sn.nodeToStation).flatten()

        dTput = np.zeros((Mq, R))
        dRespT = np.zeros((Mq, R))
        dQLen = np.zeros((Mq, R))
        dUtil = np.zeros((Mq, R))
        mask = np.zeros((Mq, R), dtype=bool)
        rates = np.zeros((Mq, R))
        for ist, node in enumerate(queue_nodes):
            s_idx = int(node_to_station[node])
            for r in range(R):
                rates[ist, r] = sn.rates[s_idx, r]
                mask[ist, r] = (np.isfinite(rates[ist, r]) and rates[ist, r] > 0
                                and D[ist, r] > 0)

        # chain_of[r] is the chain class r belongs to; Dc, Zc, Nc and lambda_c
        # aggregate the per-class quantities over the classes of each chain.
        Zrow = np.atleast_2d(np.asarray(Z, dtype=float)).sum(axis=0)
        chain_of = np.zeros(R, dtype=int)
        Dc = np.zeros((Mq, C))
        Zc = np.zeros(C)
        Nc = np.zeros(C)
        lambda_c = np.zeros(C)
        for c in range(C):
            cls = np.flatnonzero(np.asarray(chains[c, :]).flatten())
            for r in cls:
                chain_of[r] = c
            if len(cls) == 0:
                continue
            Dc[:, c] = D[:, cls].sum(axis=1)
            Zc[c] = float(Zrow[cls].sum())
            finite = Np[cls][np.isfinite(Np[cls])]
            Nc[c] = float(finite.sum()) if finite.size else 0.0
            lambda_c[c] = float(lam[cls].sum())

        if not is_open:
            # ---- closed branch: differentiated MVA at chain level ---------
            sens = pfqn_sens(Dc, Nc, Zc)
            for ist in range(Mq):
                for r in range(R):
                    if not mask[ist, r]:
                        continue
                    c = int(chain_of[r])
                    if Dc[ist, c] <= 0:
                        continue
                    rate = rates[ist, r]
                    Dir = D[ist, r]
                    visits = Dir * rate            # chain-normalized visit ratio
                    p = ist * C + c
                    chain = -Dir / rate            # dDc(i,c)/dmu(i,r)
                    Xc = sens.X[c] if np.ndim(sens.X) == 1 else sens.X[0, c]
                    Qc = sens.Q[ist, c]
                    dXc = sens.dX[c, p] * chain
                    dQc = sens.dQ[ist, c, p] * chain
                    # Class share of the chain queue at this station, and its own
                    # dependence on the rate.
                    alpha = Dir / Dc[ist, c]
                    dalpha = chain * (Dc[ist, c] - Dir) / Dc[ist, c] ** 2
                    Qir = alpha * Qc
                    dQir = dalpha * Qc + alpha * dQc
                    Tir = Xc * visits
                    dTir = dXc * visits
                    dTput[ist, r] = dTir
                    dQLen[ist, r] = dQir
                    dUtil[ist, r] = dXc * Dir + Xc * chain
                    if Tir > 0:
                        # Per-visit response time by Little's law, R = Q/T.
                        dRespT[ist, r] = (dQir * Tir - Qir * dTir) / Tir ** 2
        else:
            # ---- open branch: exact closed-form BCMP ----------------------
            # rho(i,r) = lambda_c(c)*D(i,r) with c the chain of r; U(i) is the
            # sum over classes. The stations decouple, so only the own service
            # rate mu(i,r) moves the measures at (i,r). The throughput
            # T(i,r) = lambda_c(c)*v(i,r) is fixed by the arrival rate, hence its
            # rate derivative is exactly zero.
            rho = np.zeros((Mq, R))
            for ist in range(Mq):
                for r in range(R):
                    if D[ist, r] > 0:
                        rho[ist, r] = lambda_c[int(chain_of[r])] * D[ist, r]
            Ui = rho.sum(axis=1)
            for ist in range(Mq):
                denom = 1.0 - Ui[ist]
                for r in range(R):
                    if not mask[ist, r]:
                        continue
                    rate = rates[ist, r]
                    svct = 1.0 / rate              # per-visit service time
                    drho = -rho[ist, r] / rate
                    dU = drho                      # own class only
                    dsvct = -svct / rate
                    dRespT[ist, r] = (dsvct * denom + svct * dU) / denom ** 2
                    dQLen[ist, r] = (drho * denom + rho[ist, r] * dU) / denom ** 2
                    dUtil[ist, r] = drho
                    dTput[ist, r] = 0.0            # open Tput = lambda*visits
            sens = None
        return dTput, dRespT, dQLen, dUtil, mask, sens

    def _sensitivityFD(self, sn, queue_nodes, R, step, scheme):
        """Finite-difference branch: re-run the calling solver on rate-perturbed
        copies of the model. The perturbation is a pure time scaling of the
        service process, so the shape of the distribution, and in particular its
        SCV, is preserved and only the rate moves.
        """
        import numpy as np
        from ..distributions.scaling import dist_scale_rate

        Mq = len(queue_nodes)
        central = (scheme == 'central')

        h = step
        if h is None:
            h = 1e-2 if self._isSimulationSolver() else 1e-4
        if not np.isscalar(h) or not np.isfinite(h) or h <= 0 or h >= 1:
            raise ValueError("The finite-difference step must be a scalar in "
                             "(0,1).")
        h = float(h)
        if self._isSimulationSolver():
            # Common random numbers: the base and perturbed runs must share a
            # seed, otherwise the difference quotient measures Monte Carlo
            # noise.
            seed = getattr(self.options, 'seed', None)
            if seed is None or not np.isfinite(seed) or seed <= 0:
                self.options.seed = 23000

        node_to_station = np.asarray(sn.nodeToStation).flatten()
        visited = self._sensitivityVisitMask(sn)
        mask = np.zeros((Mq, R), dtype=bool)
        for ist, node in enumerate(queue_nodes):
            s_idx = int(node_to_station[node])
            for r in range(R):
                rate = sn.rates[s_idx, r]
                mask[ist, r] = (np.isfinite(rate) and rate > 0
                                and visited[s_idx, r])

        dTput = np.zeros((Mq, R))
        dRespT = np.zeros((Mq, R))
        dQLen = np.zeros((Mq, R))
        dUtil = np.zeros((Mq, R))

        Q0, U0, R0, T0 = self._sensitivitySolveOnce()

        stations = self.model.get_stations()
        classes = self.model.get_classes()
        for ist, node in enumerate(queue_nodes):
            s_idx = int(node_to_station[node])
            station = stations[s_idx]
            for r in range(R):
                if not mask[ist, r]:
                    continue
                jobclass = classes[r]
                rate = sn.rates[s_idx, r]
                base = station.get_service(jobclass)
                try:
                    self._sensitivitySetService(
                        station, jobclass, dist_scale_rate(base, 1.0 + h))
                    Qp, Up, Rp, Tp = self._sensitivitySolveOnce()

                    if central:
                        self._sensitivitySetService(
                            station, jobclass, dist_scale_rate(base, 1.0 - h))
                        Qm, Um, Rm, Tm = self._sensitivitySolveOnce()
                        denom = 2.0 * rate * h
                    else:
                        Qm, Um, Rm, Tm = Q0, U0, R0, T0
                        denom = rate * h

                    dTput[ist, r] = (Tp[s_idx, r] - Tm[s_idx, r]) / denom
                    dRespT[ist, r] = (Rp[s_idx, r] - Rm[s_idx, r]) / denom
                    dQLen[ist, r] = (Qp[s_idx, r] - Qm[s_idx, r]) / denom
                    dUtil[ist, r] = (Up[s_idx, r] - Um[s_idx, r]) / denom
                finally:
                    # restores the unperturbed service process
                    self._sensitivitySetService(station, jobclass, base)
        return dTput, dRespT, dQLen, dUtil, mask

    def _sensitivitySolveOnce(self):
        """Solve with a fresh instance of the calling solver class, carrying its
        options over so that seeds, tolerances and method selection are those of
        the caller. A fresh instance is used because a solver caches its results,
        and reset() alone would not discard a warm start.
        """
        import numpy as np
        solver = type(self)(self.model, self.options)
        QN, UN, RN, TN = solver.getAvg()[:4]
        return (np.atleast_2d(np.asarray(QN, dtype=float)),
                np.atleast_2d(np.asarray(UN, dtype=float)),
                np.atleast_2d(np.asarray(RN, dtype=float)),
                np.atleast_2d(np.asarray(TN, dtype=float)))

    def _sensitivitySetService(self, station, jobclass, distrib):
        station.set_service(jobclass, distrib)
        self.model.refresh_rates()

    def _isSimulationSolver(self):
        return type(self).__name__ in ('SolverSSA', 'SolverJMT', 'SolverLDES')

    @staticmethod
    def _sensitivityVisitMask(sn):
        """A (nstations x nclasses) mask of the station-class pairs that carry
        visits, used in place of the demand matrix, which a non-product-form
        model may not admit.
        """
        import numpy as np
        from ..constants import GlobalConstants
        stateful = np.zeros((int(sn.nstateful), int(sn.nclasses)), dtype=bool)
        for c in range(int(sn.nchains)):
            V = sn.visits.get(c) if hasattr(sn.visits, 'get') else sn.visits[c]
            if V is None:
                continue
            V = np.atleast_2d(np.asarray(V, dtype=float))
            if V.size == 0:
                continue
            stateful = stateful | (V > GlobalConstants.Zero)
        # sn.visits is indexed by STATEFUL node, the mask by station.
        station_to_stateful = np.asarray(sn.stationToStateful).flatten()
        visited = np.zeros((int(sn.nstations), int(sn.nclasses)), dtype=bool)
        for i in range(int(sn.nstations)):
            visited[i, :] = stateful[int(station_to_stateful[i]), :]
        return visited

    def getMomentTable(self, order=None):
        """Exact higher moments of the per-class performance measures.

        Returns ``(MomentTable, mom)`` where MomentTable is a DataFrame with one
        row per (Station, JobClass) giving, in addition to the means that
        getAvgTable reports, the second moments of that row's queue length and
        response time: QLen, QLenVar, QLenSCV, RespT, RespTVar, RespTSCV.

        ``order`` selects which moment orders to report. It is a SET: a scalar
        ``k`` is read as ``1:k``, "everything up to order k"; an explicit
        list/vector selects exactly those orders::

            1        the means only:            QLen, RespT
            2        (default) means and second moments, i.e. the columns above
            3        also adds RespTSkew
            [1, 2]   the same as 2
            [2, 3]   second moments and skewness, without the means

        Order 1 contributes QLen and RespT, order 2 contributes the Var and SCV
        columns, order 3 contributes RespTSkew.

        ``order = 3`` also adds QLenSkew, the skewness of the per-class queue
        length. That quantity is reachable because the generating parameter need
        not scale a whole demand column: scaling ``L(i,r)`` alone is Theorem 1 of
        Akyildiz and Strelen with the class subset ``T = {r}``, and it generates
        the moments of ``n(i,r)`` itself. QLenSkew is available only for closed
        single-server models, which is the scope of pfqn_sens_mom; it is NaN
        otherwise.

        All of it is exact, not simulated and not approximated. The queue-length
        moments come from the product-form identity
        ``Cov[n(i,r),n(j,s)] = L(j,s) dQ(i,r)/dL(j,s)``, evaluated by the
        pfqn_sens_* family; see ``_kb/03-api-layer.md``.

        RESPONSE-TIME MOMENTS ARE FCFS OR PROCESSOR-SHARING. RespTVar and
        RespTSCV are NaN at any station that is neither, and at an LCFS center in
        particular, because the sojourn-time distribution there is not known in
        general (Strelen 1990, Section 4) and a wrong value is worse than a
        blank. The mean RespT is always reported, since it needs no
        distributional result.

        The FCFS moments come from pfqn_sens_respt and are closed-model only. The
        processor-sharing moments come from Mitra and Morrison (1983) and cover
        two configurations, both requiring exponential single-server service:
        purely open, where qsys_mm1_ps is exact at any PS station whose arrivals
        are Poisson, that is, that lies on no routing cycle; and purely closed,
        where pfqn_respt_ps_moments covers the terminal-driven system the paper
        analyses, one PS station visited once per think cycle with delay stations
        holding the think time. A PS station outside those configurations keeps
        RespTVar = NaN.

        Scope by model type::

            closed, single-server            -> pfqn_sens_mva
            closed, multiserver              -> pfqn_sens_mvaldmx
            mixed open and closed            -> pfqn_sens_mvaldmx
            purely open, single-server       -> exact BCMP closed form (below)
            purely open, multiserver         -> not supported, see the error

        ``mom`` is a dict carrying the raw results: ``mom['qlen']`` is the
        underlying pfqn_sens_mva / pfqn_sens_mvaldmx object (with the full
        covariance matrices, not just the diagonal this table shows),
        ``mom['respt']`` is the pfqn_sens_respt object or None, and
        ``mom['psrespt']`` is the pfqn_respt_ps_moments object or None. Use it
        when the per-pair covariances, or the route taken at a PS station, are
        needed.

        Per-station TOTAL moments, including the third moment and the skewness,
        are in getMomentStationTable: they are only defined for a station total,
        because the parameter that generates them scales a whole demand column.

        See also: getAvgTable, getSensitivityTable, getMomentStationTable.
        """
        import numpy as np
        import pandas as pd
        from ..api.sn.transforms import sn_get_product_form_params
        from ..api.pfqn.sens_mva import pfqn_sens_mva
        from ..api.pfqn.sens_mvaldmx import pfqn_sens_mvaldmx
        from ..api.pfqn.sens_respt import pfqn_sens_respt
        from ..api.pfqn.sens_mom import pfqn_sens_mom
        from ..api.sn.network_struct import NodeType

        if order is None:
            order = 2
        order = _validate_moment_order(order, 3)
        maxorder = max(order)

        sn = self.model.getStruct()
        R = int(sn.nclasses)
        N = np.asarray(sn.njobs, dtype=float).flatten()

        lam, D, Np_, Z, mu, Ssrv, _ = sn_get_product_form_params(sn)
        lam = np.asarray(lam, dtype=float).flatten()
        D = np.atleast_2d(np.asarray(D, dtype=float))
        Np_ = np.asarray(Np_, dtype=float).flatten()
        mu = np.atleast_2d(np.asarray(mu, dtype=float))
        Ssrv = np.asarray(Ssrv, dtype=float).flatten()
        queue_nodes = [i for i, nt in enumerate(sn.nodetype)
                       if nt == NodeType.QUEUE]
        Mq = len(queue_nodes)
        Ztot = np.atleast_2d(np.asarray(Z, dtype=float)).sum(axis=0)
        is_open = bool(np.any(np.isinf(N)))
        is_closed = bool(np.any(np.isfinite(N) & (N > 0)))
        is_mixed = is_open and is_closed

        nodenames = list(sn.nodenames)
        classnames = list(sn.classnames)
        node_to_station = np.asarray(sn.nodeToStation).flatten()
        rates = np.asarray(sn.rates, dtype=float)

        mom = {'qlen': None, 'respt': None, 'qlenmom': None, 'psrespt': None}
        QLen = np.zeros((Mq, R))
        QLenVar = np.zeros((Mq, R))

        # ---- queue-length moments ------------------------------------------
        if not is_open:
            if np.all(Ssrv == 1):
                mom['qlen'] = pfqn_sens_mva(D, Np_, Ztot)
            else:
                mom['qlen'] = pfqn_sens_mvaldmx(np.zeros(R), D, Np_, Ztot,
                                                mu, Ssrv)
            QLen = np.asarray(mom['qlen'].Q, dtype=float)
            QLenVar = np.asarray(mom['qlen'].QVar, dtype=float)
        elif is_mixed:
            mom['qlen'] = pfqn_sens_mvaldmx(lam, D, N, Ztot, mu, Ssrv)
            QLen = np.asarray(mom['qlen'].Q, dtype=float)
            QLenVar = np.asarray(mom['qlen'].QVar, dtype=float)
        else:
            # Purely open. The moment analysis of the pfqn_sens_* family needs at
            # least one closed class to have a population lattice to recurse on,
            # but a purely open BCMP single-server station has a closed-form
            # joint law that needs no recursion at all: with
            # rho(i,r) = lambda(r)*D(i,r) and rho_i = sum_r rho(i,r) < 1,
            #   P(n_i) = (1-rho_i) * (sum_r n(i,r))! / prod_r n(i,r)!
            #            * prod_r rho(i,r)^n(i,r)
            # so the total n_i is geometric with parameter rho_i and,
            # conditionally on n_i, the classes are multinomial with
            # p_r = rho(i,r)/rho_i. Compounding,
            #   Cov[n(i,r),n(i,s)] = E[n_i]*(delta_rs*p_r - p_r*p_s)
            #                        + p_r*p_s*Var[n_i]
            if np.any(Ssrv > 1):
                raise ValueError(
                    "getMomentTable does not support multiserver stations in a "
                    "purely open model: the queue-length law is not geometric "
                    "there. Add a closed class, or use a single-server model.")
            rho = np.zeros((Mq, R))
            for ist in range(Mq):
                for r in range(R):
                    if np.isinf(N[r]):
                        rho[ist, r] = lam[r] * D[ist, r]
            rho_tot = rho.sum(axis=1)
            for ist in range(Mq):
                ri = float(rho_tot[ist])
                if ri >= 1:
                    raise ValueError(
                        "Station %s is unstable (utilization %.4f >= 1); its "
                        "queue-length moments do not exist."
                        % (str(nodenames[queue_nodes[ist]]), ri))
                if ri <= 0:
                    continue
                En = ri / (1 - ri)
                Vn = ri / (1 - ri) ** 2
                for r in range(R):
                    pr = rho[ist, r] / ri
                    QLen[ist, r] = pr * En
                    QLenVar[ist, r] = En * (pr - pr ** 2) + pr ** 2 * Vn

        # ---- per-class queue-length skewness (closed, single-server only) ---
        # pfqn_sens_mom with groups = 1:R scales one class at a time, which is
        # the class-subset parameter T = {r} of Akyildiz and Strelen's Theorem 1,
        # so it yields the moments of n(i,r) rather than of the station total.
        QLenSkew = np.full((Mq, R), np.nan)
        if 3 in order and not is_open and np.all(Ssrv == 1):
            mom['qlenmom'] = pfqn_sens_mom(D, Np_, Ztot, np.ones(Mq),
                                           np.arange(1, R + 1))
            # single-class models give G == 1, where pfqn_sens_mom collapses the
            # group axis; reshape back so indexing by (station, class) is uniform
            QLenSkew = np.asarray(mom['qlenmom'].Skew,
                                  dtype=float).reshape(Mq, R)

        # ---- response-time moments (FCFS only) -----------------------------
        RespT = np.zeros((Mq, R))
        RespTVar = np.full((Mq, R), np.nan)
        RespTSkew = np.full((Mq, R), np.nan)
        # pfqn_sens_respt carries one service time S(i) per station, but S enters
        # its queue-length recursion ONLY through the product
        # rho(i,r) = S(i)*V(i,r) = the demand. The rate mu(i) = 1/S(i) is read
        # directly by equation (4.5) alone, and (4.5) is only evaluated at FCFS
        # stations. So a non-FCFS station may keep class-dependent demands: give
        # it any nominal S and let V absorb the rest. Only the FCFS stations must
        # have class-independent rates, which BCMP requires of them anyway.
        resp_avail = (not is_open) and _fcfs_rates_are_class_independent(
            sn, queue_nodes, R, node_to_station, rates)
        if resp_avail:
            Ssvc = np.ones(Mq)
            Vq = np.zeros((Mq, R))
            for ist, node in enumerate(queue_nodes):
                s_idx = int(node_to_station[node])
                if _sched_is_fcfs(sn, s_idx):
                    st = _service_time_of(rates, s_idx, R)
                    if st > 0:
                        Ssvc[ist] = st     # the true per-visit rate; (4.5) needs it
                for r in range(R):
                    Vq[ist, r] = D[ist, r] / Ssvc[ist]
            mom['respt'] = pfqn_sens_respt(Ssvc, Vq, Np_, Ztot, Ssrv, maxorder)

        Xq = None
        if mom['qlen'] is not None:
            Xq = np.asarray(mom['qlen'].X, dtype=float).flatten()

        def throughput_of(r):
            if np.isinf(N[r]):
                return float(lam[r])
            if Xq is not None:
                return float(Xq[r])
            return 0.0

        for ist, node in enumerate(queue_nodes):
            s_idx = int(node_to_station[node])
            for r in range(R):
                if D[ist, r] <= 0:
                    continue
                # Mean response time per visit, by Little's law at the station:
                # the arrival rate of class r to station i is X(r) times the
                # visit ratio V(i,r) = D(i,r)/S(i,r) = D(i,r)*rate(i,r). This
                # needs no distributional result, so it is always reported.
                xr = throughput_of(r)
                Vir = D[ist, r] * rates[s_idx, r]
                if xr > 0 and Vir > 0:
                    RespT[ist, r] = QLen[ist, r] / (xr * Vir)
            # The variance needs the sojourn-time distribution, which is known
            # only at FCFS centers; elsewhere RespTVar stays NaN.
            if mom['respt'] is not None and _sched_is_fcfs(sn, s_idx):
                for r in range(R):
                    if D[ist, r] > 0:
                        RespT[ist, r] = mom['respt'].W[ist, r]
                        if maxorder >= 2:
                            RespTVar[ist, r] = mom['respt'].WVar[ist, r]
                        if maxorder >= 3:
                            RespTSkew[ist, r] = mom['respt'].WSkew[ist, r]

        # ---- response-time moments at processor-sharing stations ------------
        # Mitra and Morrison (1983) supply the sojourn-time moments the FCFS
        # block above cannot reach: exactly at an open PS station fed by Poisson
        # streams, and by expansion (or by exact enumeration when small) in the
        # closed terminal-driven system. see _kb/05-solvers-overview.md
        mom['psrespt'] = None
        if 2 in order and not is_mixed and Mq > 0:
            from ..api.qsys.ps import qsys_mm1_ps
            from ..api.pfqn.respt_ps import pfqn_respt_ps_moments
            ps_queues = [ist for ist, node in enumerate(queue_nodes)
                         if _sched_is_ps(sn, int(node_to_station[node]))]
            if is_open:
                for ist in ps_queues:
                    s_idx = int(node_to_station[queue_nodes[ist]])
                    visiting = [r for r in range(R) if D[ist, r] > 0]
                    if (Ssrv[ist] != 1
                            or not _rates_are_exponential(sn, s_idx, visiting)
                            or not _station_is_feedback_free(sn, s_idx, R)):
                        continue
                    mu_st = np.ones(R)
                    lam_st = np.zeros(R)
                    for r in visiting:
                        mu_st[r] = rates[s_idx, r]
                        lam_st[r] = lam[r] * D[ist, r] * rates[s_idx, r]
                    if float(np.sum(lam_st / mu_st)) >= 1:
                        continue
                    Wps, W2ps, _ = qsys_mm1_ps(lam_st, mu_st)
                    for r in visiting:
                        RespTVar[ist, r] = W2ps[r] - Wps[r] ** 2
            elif len(ps_queues) == 1 and Mq == 1:
                # the paper's closed system: terminals in series with one PS CPU
                ist = ps_queues[0]
                s_idx = int(node_to_station[queue_nodes[ist]])
                visiting = [r for r in range(R) if D[ist, r] > 0]
                Vps = np.array([D[ist, r] * rates[s_idx, r] for r in visiting])
                single_visit = bool(np.all(np.abs(Vps - 1) <= 1e-3))
                if (Ssrv[ist] == 1 and _rates_are_exponential(sn, s_idx, visiting)
                        and single_visit
                        and all(Ztot[r] > 0 for r in visiting)):
                    Sps = np.ones(R)
                    Zps = np.ones(R)
                    Nps = np.zeros(R)
                    for r in visiting:
                        Sps[r] = 1.0 / rates[s_idx, r]
                        Zps[r] = Ztot[r]
                        Nps[r] = N[r]
                    Wps, W2ps, mom['psrespt'] = pfqn_respt_ps_moments(Sps, Nps, Zps)
                    for r in visiting:
                        RespTVar[ist, r] = W2ps[r] - Wps[r] ** 2

        # ---- assemble -------------------------------------------------------
        rows = []
        for ist, node in enumerate(queue_nodes):
            for r in range(R):
                if D[ist, r] <= 0:
                    continue    # class r does not visit this station
                rows.append({
                    'Station': str(nodenames[node]),
                    'JobClass': str(classnames[r]),
                    'QLen': float(QLen[ist, r]),
                    'QLenVar': float(QLenVar[ist, r]),
                    'QLenSCV': _moment_scv(QLenVar[ist, r], QLen[ist, r]),
                    'QLenSkew': float(QLenSkew[ist, r]),
                    'RespT': float(RespT[ist, r]),
                    'RespTVar': float(RespTVar[ist, r]),
                    'RespTSCV': _moment_scv(RespTVar[ist, r], RespT[ist, r]),
                    'RespTSkew': float(RespTSkew[ist, r]),
                })
        # Column order is fixed; `order` filters it. Order 1 contributes QLen and
        # RespT, order 2 the Var/SCV columns, order 3 QLenSkew and RespTSkew.
        cols = ['Station', 'JobClass']
        if 1 in order:
            cols.append('QLen')
        if 2 in order:
            cols += ['QLenVar', 'QLenSCV']
        if 3 in order:
            cols.append('QLenSkew')
        if 1 in order:
            cols.append('RespT')
        if 2 in order:
            cols += ['RespTVar', 'RespTSCV']
        if 3 in order:
            cols.append('RespTSkew')
        MomentTable = pd.DataFrame(rows, columns=cols)
        return MomentTable, mom

    get_moment_table = getMomentTable

    def getMomentStationTable(self, order=None):
        """Exact higher moments of the total queue length.

        Returns ``(MomentStationTable, mom)`` where MomentStationTable is a
        DataFrame with one row per Station giving the moments of the TOTAL queue
        length at that station, ``Q_i = sum_r n(i,r)``: QLen, QLenVar, QLenSCV.

        ``order`` selects which moment orders to report. It is a SET: a scalar
        ``k`` is read as ``1:k``, "everything up to order k"; an explicit
        list/vector selects exactly those orders::

            1        the mean only:             QLen
            2        (default) mean and second moment: QLen, QLenVar, QLenSCV
            3        also adds QLenM3 and QLenSkew
            [1, 3]   the mean and the third moment, without the variance

        Order 1 contributes QLen, order 2 contributes QLenVar and QLenSCV, order
        3 contributes QLenM3 and QLenSkew.

        Why this is a separate table from getMomentTable. The moments beyond the
        second are generated by differentiating with respect to ``x_i``, the
        reciprocal of the capacity of station i, which scales the service times
        of ALL classes at that station at once. That parameter therefore produces
        moments of the station total, not of any one class: there is no per-class
        third moment to report, and inventing one by splitting the total would be
        fiction. The per-class second moments, which do exist, are in
        getMomentTable. The two are consistent: ``Var[Q_i]`` here equals the sum
        of getMomentTable's per-class covariances at station i over all class
        pairs.

        The ALGORITHM is chosen by the solver's method, set at construction, not
        by an argument here: a method is a property of the solver object, so
        passing one per call would let a single solver answer with two different
        algorithms::

            SolverMVA(model)                 -> pfqn_sens_mom, exact, but it
                                                walks the whole population
                                                lattice at a cost of prod(N+1),
                                                so it is unusable once the
                                                populations are large.
            SolverMVA(model, method='lin')   -> pfqn_sens_linearizer,
                                                approximate and polynomial-time.
                                                The reference reports relative
                                                errors below 2.1% on E[Q], 4.1%
                                                on E[Q^2] and 6.2% on E[Q^3].

        Any Linearizer-family method ('lin', 'amva.lin', 'egflin', 'gflin') takes
        the approximate path; every other method takes the exact one.

        Restricted to closed models. Mixed and open second moments are available
        per class from getMomentTable; the higher moments of the reference are
        stated for closed networks only.

        ``mom`` is the underlying pfqn_sens_mom / pfqn_sens_linearizer object,
        which also carries the cross-station covariance matrix ``.Cov`` that this
        table does not show.

        Reference: J. C. Strelen, "Moment Analysis for Closed Queuing Networks
        and its Linearizer", Performance Evaluation 11:127-142, 1990,
        Theorem 3.1 and equation (3.2).

        See also: getMomentTable, getAvgTable, getSensitivityTable.
        """
        import numpy as np
        import pandas as pd
        from ..api.sn.transforms import sn_get_product_form_params
        from ..api.pfqn.sens_mom import pfqn_sens_mom
        from ..api.pfqn.sens_linearizer import pfqn_sens_linearizer
        from ..api.sn.network_struct import NodeType
        from ..api.sn.predicates import sn_has_product_form

        if order is None:
            order = 2
        order = _validate_moment_order(order, 3)
        # The algorithm is a property of the solver, not of this call: it comes
        # from the method set at construction, e.g. SolverMVA(model,
        # method='lin'). Taking a per-call method argument here would have given
        # the same solver object two disagreeing methods.
        method = self.options.method

        sn = self.model.getStruct()
        N = np.asarray(sn.njobs, dtype=float).flatten()

        _, D, Np_, Z, _, Ssrv, _ = sn_get_product_form_params(sn)
        D = np.atleast_2d(np.asarray(D, dtype=float))
        Np_ = np.asarray(Np_, dtype=float).flatten()
        Ssrv = np.asarray(Ssrv, dtype=float).flatten()
        queue_nodes = [i for i, nt in enumerate(sn.nodetype)
                       if nt == NodeType.QUEUE]
        Ztot = np.atleast_2d(np.asarray(Z, dtype=float)).sum(axis=0)

        # The moments rest on the product-form identity Cov = L dQ/dL, which is
        # a theorem about the product form; outside it, L dQ/dL is still
        # computable and simply is NOT a covariance. See
        # _not_product_form_error.
        if not sn_has_product_form(sn):
            raise _not_product_form_error('getMomentStationTable')
        if np.any(np.isinf(N)):
            raise ValueError(
                "getMomentStationTable supports closed models only. Per-class "
                "second moments of an open or mixed model are available from "
                "getMomentTable.")
        if np.any(Ssrv > 1):
            raise ValueError(
                "getMomentStationTable supports single-server stations only: "
                "the higher-moment recursion of the reference is stated for "
                "load-independent stations. Per-class second moments of a "
                "multiserver model are available from getMomentTable.")

        # The solver's method selects the algorithm: a Linearizer-family method
        # takes the approximate polynomial-time path, an exact-MVA method takes
        # the exact recursion, and every other method is differentiated
        # numerically off its own means.
        if _is_linearizer_method(method):
            mom = pfqn_sens_linearizer(D, Np_, Ztot)
        elif _is_exact_mva_method(method):
            mom = pfqn_sens_mom(D, Np_, Ztot)
        else:
            # A method with no hand-differentiated counterpart. The identity does
            # not care HOW the mean queue lengths were obtained, so the solver
            # itself is used as a mean-value oracle and differentiated
            # numerically. The moments inherit the accuracy of that method's
            # means.
            mom = _moments_by_finite_difference(self, sn, queue_nodes)

        nodenames = list(sn.nodenames)
        rows = []
        for ist, node in enumerate(queue_nodes):
            if np.all(D[ist, :] <= 0):
                continue    # no class visits this station
            m_i = float(mom.m[ist])
            var_i = float(mom.Var[ist])
            rows.append({
                'Station': str(nodenames[node]),
                'QLen': m_i,
                'QLenVar': var_i,
                'QLenSCV': (var_i / m_i ** 2) if m_i > 0 else float('nan'),
                'QLenM3': float(mom.M3[ist]),
                'QLenSkew': float(mom.Skew[ist]),
            })
        cols = ['Station']
        if 1 in order:
            cols.append('QLen')
        if 2 in order:
            cols += ['QLenVar', 'QLenSCV']
        if 3 in order:
            cols += ['QLenM3', 'QLenSkew']
        MomentStationTable = pd.DataFrame(rows, columns=cols)
        return MomentStationTable, mom

    get_moment_station_table = getMomentStationTable

    def solveMeansForStruct(self, sn):
        """Mean queue lengths of a perturbed structure, under this solver's own
        method. Mirrors MATLAB's ``@NetworkSolver/solveMeansForStruct``.

        This is the mean-value oracle behind the numerical-derivative path of
        getMomentChainTable and getMomentStationTable. The moment identity
        Cov[n,n] = L dQ/dL does not care HOW the mean queue lengths were
        obtained, only that they are the means of a product-form model as a
        function of its demands. So rather than hand-differentiating each
        algorithm, the solver is re-run on a perturbed structure and
        differentiated numerically.

        Running THIS solver, rather than one chosen algorithm, is what makes the
        path general: it covers every method of every solver, including the
        normalizing-constant methods of SolverNC (comom, ca, le, mom, ...) and
        the summation methods of SolverMVA (sum, esum), none of which any
        hand-differentiated implementation reaches. Restricting the oracle to one
        analyzer would restrict the moments to that analyzer's methods for no
        mathematical reason.

        The perturbed ``sn`` is injected by copying the model and overwriting its
        cached structure, then constructing a fresh solver of this class with
        this solver's options. Every step is load-bearing:

        - ``model.copy()`` because Python objects are by-reference and numpy
          arrays alias, so writing ``_sn`` on the caller's model would corrupt
          every later solve off it;
        - ``_has_struct = True`` (MATLAB's ``hasStruct``) so ``get_struct()``
          does not regenerate the struct and discard the perturbation;
        - ``_rates_dirty = False`` because a solver's ``_extract_network_params``
          calls ``model.refresh_rates()`` when the flag is set, which would
          recompute ``rates`` from the model's distributions and discard the
          perturbation just as surely. ``Network.copy()`` already leaves it
          False; it is pinned here because the failure mode is SILENT -- the
          means come back unperturbed and every finite difference is exactly
          zero;
        - a FRESH solver because a solver caches its results (and its demands)
          at construction and would otherwise return the unperturbed answer.

        Returns the (nstations x nclasses) mean queue lengths.
        """
        import numpy as np
        model = self.model.copy()
        model._sn = sn
        model._has_struct = True
        model._rates_dirty = False
        solver = type(self)(model, self.options)
        return np.atleast_2d(np.asarray(solver.getAvgQLen(), dtype=float))

    solve_means_for_struct = solveMeansForStruct

    def getMomentChainTable(self, order=None):
        """Exact higher moments of the per-chain queue length.

        Returns ``(MomentChainTable, mom)`` where MomentChainTable is a DataFrame
        with one row per (Station, Chain) giving the moments of the queue length
        of that chain at that station, ``Q_(i,c) = sum_(r in chain c) n(i,r)``:
        QLen, QLenVar, QLenSCV.

        This is the chain-level analogue of getAvgChainTable, and it sits between
        the two other moment tables: getMomentTable is per class,
        getMomentStationTable is per station total, and this one is per chain,
        i.e. per group of classes that circulate together.

        ``order`` selects which moment orders to report. It is a SET: a scalar
        ``k`` is read as ``1:k``, "everything up to order k"; an explicit
        list/vector selects exactly those orders::

            1        the mean only:             QLen
            2        (default) mean and second moment: QLen, QLenVar, QLenSCV
            3        also adds QLenM3 and QLenSkew

        Unlike the per-class table, order 3 IS available here. All three tables
        are the same recursion under different groupings of the classes: the
        generating parameter scales the service times of a class subset ``T`` at
        a station, and the moments it produces are those of
        ``sum_(r in T) n(i,r)``. ``T = {r}`` gives getMomentTable, ``T = chain``
        gives this table, ``T = all classes`` gives getMomentStationTable. That is
        Theorem 1 of Akyildiz and Strelen; Strelen's own ``x_i`` is the last case.

        The ALGORITHM is chosen by the solver's method, set at construction, not
        by an argument here; see getMomentStationTable. A Linearizer-family
        method approximates the per-station totals only, so it cannot express a
        per-chain grouping and is rejected here unless every class already sits
        in one chain, in which case the chain IS the station total.

        Restricted to closed, single-server models, which is the scope of
        pfqn_sens_mom.

        ``mom`` is the underlying pfqn_sens_mom struct. Its ``.Cov`` is
        (M x C x M x C) and carries the cross-chain and cross-station
        covariances this table does not show.

        Reference: I. F. Akyildiz and J. C. Strelen, "Moment Analysis for
        Load-Dependent Mixed Product Form Queueing Networks", IEEE Trans.
        Communications 39(6):828-832, 1991, Theorem 1; J. C. Strelen, "Moment
        Analysis for Closed Queuing Networks and its Linearizer", Performance
        Evaluation 11:127-142, 1990, equation (3.2).

        See also: getMomentTable, getMomentStationTable, getAvgChainTable.
        """
        import numpy as np
        import pandas as pd
        from ..api.sn.transforms import sn_get_product_form_params
        from ..api.pfqn.sens_mom import pfqn_sens_mom
        from ..api.pfqn.sens_linearizer import pfqn_sens_linearizer
        from ..api.sn.network_struct import NodeType
        from ..api.sn.predicates import sn_has_product_form

        if order is None:
            order = 2
        order = _validate_moment_order(order, 3)
        # The algorithm is a property of the solver, not of this call: it comes
        # from the method set at construction, e.g. SolverMVA(model,
        # method='lin'). Taking a per-call method argument here would have given
        # the same solver object two disagreeing methods.
        method = self.options.method

        sn = self.model.getStruct()
        R = int(sn.nclasses)
        N = np.asarray(sn.njobs, dtype=float).flatten()

        _, D, Np_, Z, _, Ssrv, _ = sn_get_product_form_params(sn)
        D = np.atleast_2d(np.asarray(D, dtype=float))
        Np_ = np.asarray(Np_, dtype=float).flatten()
        Ssrv = np.asarray(Ssrv, dtype=float).flatten()
        queue_nodes = [i for i, nt in enumerate(sn.nodetype)
                       if nt == NodeType.QUEUE]
        Mq = len(queue_nodes)
        Ztot = np.atleast_2d(np.asarray(Z, dtype=float)).sum(axis=0)

        # See getMomentStationTable: the moment identity Cov = L dQ/dL is a
        # theorem about the product form, so outside it L dQ/dL is not a
        # covariance and no correct value exists to return.
        if not sn_has_product_form(sn):
            raise _not_product_form_error('getMomentChainTable')
        if np.any(np.isinf(N)):
            raise ValueError(
                "getMomentChainTable supports closed models only. Per-class "
                "second moments of an open or mixed model are available from "
                "getMomentTable.")
        if np.any(Ssrv > 1):
            raise ValueError(
                "getMomentChainTable supports single-server stations only: the "
                "higher-moment recursion of the reference is stated for "
                "load-independent stations.")

        # the chain of each class; sn.chains is (nchains x nclasses)
        chains = np.atleast_2d(np.asarray(sn.chains))
        classnames = list(sn.classnames)
        groups = np.zeros(R, dtype=int)
        for r in range(R):
            c = np.nonzero(chains[:, r])[0]
            if c.size == 0:
                raise ValueError("class %s belongs to no chain"
                                 % str(classnames[r]))
            groups[r] = int(c[0]) + 1      # 1-based, as pfqn_sens_mom expects
        # pfqn_sens_mom requires the group labels to be consecutive from 1, so
        # drop any chain that holds no class rather than leaving a hole in the
        # numbering
        used, compact = np.unique(groups, return_inverse=True)
        groups = compact + 1
        Cg = len(used)

        # See getMomentStationTable: the solver's method selects the algorithm.
        # The Linearizer approximates the per-station totals only, so it cannot
        # express a per-chain grouping unless every class already sits in one
        # chain, in which case the chain IS the station total.
        if _is_linearizer_method(method):
            if Cg > 1:
                raise ValueError(
                    "the solver method '%s' approximates the per-station "
                    "totals, so it cannot produce a per-chain grouping of %d "
                    "chains. Use an exact method, or getMomentStationTable."
                    % (method, Cg))
            mom = pfqn_sens_linearizer(D, Np_, Ztot)
        elif _is_exact_mva_method(method):
            mom = pfqn_sens_mom(D, Np_, Ztot, np.ones(Mq), groups)
        else:
            # See getMomentStationTable: the solver's own method supplies the
            # means and the derivatives are taken numerically. Here the parameter
            # scales the demands of one CHAIN's classes at one station, which is
            # the class subset T of Akyildiz-Strelen Theorem 1, so the moments it
            # generates are those of that chain's queue length.
            mom = _chain_moments_by_finite_difference(self, sn, queue_nodes,
                                                      groups, Cg)

        # a single group collapses the group axis in pfqn_sens_mom, and
        # pfqn_sens_linearizer is per-station only; reshape so (station, chain)
        # indexing is uniform
        m_g = np.asarray(mom.m, dtype=float).reshape(Mq, Cg)
        var_g = np.asarray(mom.Var, dtype=float).reshape(Mq, Cg)
        m3_g = np.asarray(mom.M3, dtype=float).reshape(Mq, Cg)
        skew_g = np.asarray(mom.Skew, dtype=float).reshape(Mq, Cg)

        nodenames = list(sn.nodenames)
        rows = []
        for ist, node in enumerate(queue_nodes):
            for g in range(Cg):
                classes_of = np.nonzero(groups == g + 1)[0]
                if np.all(D[ist, classes_of] <= 0):
                    continue    # no class of this chain visits this station
                m_i = float(m_g[ist, g])
                var_i = float(var_g[ist, g])
                rows.append({
                    'Station': str(nodenames[node]),
                    'Chain': 'Chain%d' % (int(used[g]),),
                    'QLen': m_i,
                    'QLenVar': var_i,
                    'QLenSCV': (var_i / m_i ** 2) if m_i > 0 else float('nan'),
                    'QLenM3': float(m3_g[ist, g]),
                    'QLenSkew': float(skew_g[ist, g]),
                })

        cols = ['Station', 'Chain']
        if 1 in order:
            cols.append('QLen')
        if 2 in order:
            cols += ['QLenVar', 'QLenSCV']
        if 3 in order:
            cols += ['QLenM3', 'QLenSkew']
        MomentChainTable = pd.DataFrame(rows, columns=cols)
        return MomentChainTable, mom

    get_moment_chain_table = getMomentChainTable
    def _make_sys_table(self, CN, XN):
        """Build the chain-level system table (Chain, JobClasses, SysRespT,
        SysTput) shared across solvers, matching the MATLAB/JAR layout.

        Args:
            CN: per-chain system response times
            XN: per-chain system throughputs
        """
        import pandas as pd
        from ..indexed_table import IndexedTable
        sn = self.model.getStruct()
        classnames = list(getattr(sn, 'classnames', []) or [])
        inchain = getattr(sn, 'inchain', {}) or {}
        rows = []
        for c in range(len(CN)):
            if c in inchain:
                idx = [int(k) for k in list(inchain[c].flatten())]
                names = [classnames[k] if k < len(classnames) else f'Class{k + 1}' for k in idx]
                jobclasses = '(' + ' '.join(names) + ')'
            else:
                jobclasses = ''
            rows.append({
                'Chain': f'Chain{c + 1}',
                'JobClasses': jobclasses,
                'SysRespT': float(CN[c]),
                'SysTput': float(XN[c]),
            })
        result = IndexedTable(pd.DataFrame(rows))
        if not self._table_silent:
            print(result)
        return result

    def _computeChainMetrics(self):
        """Compute chain-level response times and throughputs (CNchain, XNchain).

        Faithful port of MATLAB @NetworkSolver/getAvgSys.m (non fork-join path):
          1. CNclass: per-class residence time, visit-normalised to the reference
             station (source slot skipped for open classes).
          2. alpha: per-(station,class) weighting by completing-class visits.
          3. XNchain: carried system throughput = flow of completing classes
             routed back into the chain reference station (sum of rt*TN), NOT the
             offered/source rate.
          4. CNchain: open chains use sum(alpha*CNclass); closed chains apply
             Little's law nJobsChain/XNchain.

        Requires self._result to expose station matrices R (response) and T
        (throughput) and self._sn to carry rt/visits/inchain/refstat. Shared by
        the CTMC and NC solvers so their getAvgSys agrees with MATLAB/JAR.
        """
        import numpy as np
        if self._result is None:
            self.runAnalyzer()
        sn = self._sn
        if sn is None:
            sn = self.model.getStruct(True)
            self._sn = sn

        RN = self._result.R  # Station response times (M x K)
        TN = self._result.T  # Station throughputs (M x K)
        njobs = np.asarray(sn.njobs).flatten()
        nstations = sn.nstations
        nclasses = sn.nclasses
        nchains = sn.nchains if hasattr(sn, 'nchains') and sn.nchains > 0 else 1

        completes = np.ones(nclasses, dtype=bool)
        classes = self.model.get_classes()
        for r in range(nclasses):
            if r < len(classes):
                completes[r] = classes[r].completes

        def _stf(i):
            return int(sn.stationToStateful[i]) if hasattr(sn, 'stationToStateful') else i

        # CNclass: visit-normalised per-class residence time
        CNclass = np.zeros(nclasses)
        for c in range(nchains):
            if c not in sn.inchain:
                continue
            inchain = sn.inchain[c].flatten().astype(int)
            for r in inchain:
                if r >= nclasses:
                    continue
                CNclass[r] = 0.0
                if RN is not None and np.asarray(RN).size > 0 and c in sn.visits and sn.visits[c] is not None:
                    refstat = int(sn.refstat[r]) if hasattr(sn, 'refstat') else 0
                    visits_c = sn.visits[c]
                    for i in range(nstations):
                        if np.isinf(njobs[r]) and i == refstat:
                            continue
                        si, sref = _stf(i), _stf(refstat)
                        if si < visits_c.shape[0] and r < visits_c.shape[1]:
                            visit_ref = visits_c[sref, r] if sref < visits_c.shape[0] else 1.0
                            if visit_ref > 0:
                                CNclass[r] += visits_c[si, r] * RN[i, r] / visit_ref

        # alpha: per-(station,class) completing-class visit weights
        alpha = np.zeros((nstations, nclasses))
        for c in range(nchains):
            if c not in sn.inchain:
                continue
            inchain = sn.inchain[c].flatten().astype(int)
            completingclasses = [k for k in inchain if k < nclasses and completes[k]]
            if c not in sn.visits or sn.visits[c] is None:
                continue
            visits_c = sn.visits[c]
            for i in range(nstations):
                for k in inchain:
                    if k >= nclasses:
                        continue
                    refstat = int(sn.refstat[k]) if hasattr(sn, 'refstat') else 0
                    si, sref = _stf(i), _stf(refstat)
                    sum_visits = 0.0
                    for idx in completingclasses:
                        if sref < visits_c.shape[0] and idx < visits_c.shape[1]:
                            sum_visits += visits_c[sref, idx]
                    if sum_visits > 0 and si < visits_c.shape[0] and k < visits_c.shape[1]:
                        alpha[i, k] += visits_c[si, k] / sum_visits
        alpha[~np.isfinite(alpha)] = 0.0

        CNchain = np.zeros(nchains)
        XNchain = np.zeros(nchains)
        for c in range(nchains):
            if c not in sn.inchain:
                continue
            inchain = sn.inchain[c].flatten().astype(int)
            completingclasses = [k for k in inchain if k < nclasses and completes[k]]

            if TN is not None and np.asarray(TN).size > 0:
                ref = int(sn.refstat[inchain[0]]) if hasattr(sn, 'refstat') and len(inchain) > 0 else 0
                for i in range(nstations):
                    for r in completingclasses:
                        if r >= nclasses or np.isnan(TN[i, r]):
                            continue
                        if hasattr(sn, 'rt') and sn.rt is not None:
                            for s in inchain:
                                if s >= nclasses:
                                    continue
                                rt_src = i * nclasses + r
                                rt_dst = ref * nclasses + s
                                if rt_src < sn.rt.shape[0] and rt_dst < sn.rt.shape[1]:
                                    XNchain[c] += sn.rt[rt_src, rt_dst] * TN[i, r]
                        elif i == ref:
                            XNchain[c] += TN[i, r]

            nJobsChain = sum(njobs[k] for k in inchain if k < nclasses)
            if np.isinf(nJobsChain):
                refstat = int(sn.refstat[inchain[0]]) if hasattr(sn, 'refstat') and len(inchain) > 0 else 0
                sumfinite = 0.0
                for k in inchain:
                    if k >= nclasses:
                        continue
                    val = alpha[refstat, k] * CNclass[k]
                    if np.isfinite(val):
                        sumfinite += val
                CNchain[c] = sumfinite
            else:
                CNchain[c] = nJobsChain / XNchain[c] if XNchain[c] > 0 else np.inf

        return CNchain, XNchain

    def getAvgSysTable(self):
        """System-level metrics table (one row per chain).

        Default implementation for solvers whose getAvgSys() returns
        per-class vectors: chains made of a single class map one-to-one onto
        the class metrics; for multi-class chains the chain throughput is the
        sum of the completing per-class throughputs and the chain response
        time is the throughput-weighted mean of the per-class values.
        Solvers with native chain-level getAvgSys() (e.g. NC) override this.
        """
        import numpy as np
        R, T = self.getAvgSys()
        R = np.atleast_1d(np.asarray(R, dtype=float)).flatten()
        T = np.atleast_1d(np.asarray(T, dtype=float)).flatten()
        sn = self.model.getStruct()
        inchain = getattr(sn, 'inchain', {}) or {}
        nchains = int(getattr(sn, 'nchains', 0) or 0)
        if nchains <= 0 or len(inchain) < nchains:
            # No chain structure available: one chain per class
            return self._make_sys_table(R, T)
        CN = []
        XN = []
        for c in range(nchains):
            idx = [int(k) for k in list(inchain[c].flatten()) if int(k) < len(R)]
            Xc = float(np.nansum([T[k] for k in idx]))
            if len(idx) == 1:
                Cc = float(R[idx[0]])
            elif Xc > 0:
                Cc = float(np.nansum([R[k] * T[k] for k in idx]) / Xc)
            else:
                Cc = 0.0
            CN.append(Cc)
            XN.append(Xc)
        return self._make_sys_table(CN, XN)

    get_avg_sys_table = getAvgSysTable

    # Table -> T shorthands and 4-letter aliases (MATLAB/JAR/wrapper-compatible)
    def avgT(self):
        """Short alias for getAvgTable."""
        return self.getAvgTable()

    def avgSysT(self):
        """Short alias for getAvgSysTable."""
        return self.getAvgSysTable()

    def avgNodeT(self):
        """Short alias for getAvgNodeTable."""
        return self.getAvgNodeTable()

    def avgChainT(self):
        """Short alias for getAvgChainTable."""
        return self.getAvgChainTable()

    def avgNodeChainT(self):
        """Short alias for getAvgNodeChainTable."""
        return self.getAvgNodeChainTable()

    # Table -> T shorthand aliases for the auxiliary result tables
    # (moment/sensitivity/cache/item/orbit/region), mirroring aT for
    # getAvgTable. The compact and get-prefixed spellings are assigned below.
    def momentT(self, *args, **kwargs):
        """Short alias for getMomentTable."""
        return self.getMomentTable(*args, **kwargs)

    def momentChainT(self, *args, **kwargs):
        """Short alias for getMomentChainTable."""
        return self.getMomentChainTable(*args, **kwargs)

    def momentStationT(self, *args, **kwargs):
        """Short alias for getMomentStationTable."""
        return self.getMomentStationTable(*args, **kwargs)

    def sensitivityT(self, *args, **kwargs):
        """Short alias for getSensitivityTable."""
        return self.getSensitivityTable(*args, **kwargs)

    def cacheAvgT(self, *args, **kwargs):
        """Short alias for getAvgCacheTable."""
        return self.getAvgCacheTable(*args, **kwargs)

    def itemAvgT(self, *args, **kwargs):
        """Short alias for getAvgItemTable."""
        return self.getAvgItemTable(*args, **kwargs)

    def orbitAvgT(self, *args, **kwargs):
        """Short alias for getAvgOrbitTable."""
        return self.getAvgOrbitTable(*args, **kwargs)

    def lossAvgT(self, *args, **kwargs):
        """Short alias for getAvgLossTable."""
        return self.getAvgLossTable(*args, **kwargs)

    def getAvgRegionLossTable(self):
        """Table of loss (drop) metrics per finite-capacity region and class,
        for regions that drop jobs (DROP rule). Each row reports the offered
        arrival rate (carried Tput plus drop rate), the carried throughput, the
        loss rate (region drop rate) and the loss ratio (LossRate / ArvR).

        Only regions with offered traffic are listed; empty for solvers that do
        not track region drops (only the LDES simulation populates the FCR drop
        rate DropRateNfcr).
        """
        import numpy as np
        import pandas as pd
        from ..indexed_table import IndexedTable
        self.getAvgTable()
        res = getattr(self, '_ldes_result', None) or getattr(self, '_result', None)

        def _get(name):
            if res is None:
                return None
            v = getattr(res, name, None) if not isinstance(res, dict) else res.get(name)
            return None if v is None else np.atleast_2d(np.asarray(v, dtype=float))

        TN = _get('TNfcr')
        DR = _get('DropRateNfcr')
        rows = []
        if TN is not None and DR is not None:
            sn = self._sn if hasattr(self, '_sn') else None
            classnames = (list(sn.classnames) if sn is not None
                          else ['Class%d' % (r + 1) for r in range(DR.shape[1])])
            for f in range(DR.shape[0]):
                for r in range(DR.shape[1]):
                    t = float(TN[f, r])
                    d = float(DR[f, r])
                    a = t + d
                    if not np.isfinite(a) or a <= 0:
                        continue
                    rows.append({
                        'Region': 'FCRegion%d' % (f + 1),
                        'JobClass': str(classnames[r]),
                        'ArvR': a,
                        'Tput': t,
                        'LossRate': d,
                        'LossRatio': d / a,
                    })
        df = pd.DataFrame(rows, columns=['Region', 'JobClass', 'ArvR', 'Tput',
                                         'LossRate', 'LossRatio'])
        return IndexedTable(df)

    get_avg_region_loss_table = getAvgRegionLossTable

    def regionLossAvgT(self, *args, **kwargs):
        """Short alias for getAvgRegionLossTable."""
        return self.getAvgRegionLossTable(*args, **kwargs)

    def regionAvgT(self, *args, **kwargs):
        """Short alias for getAvgRegionTable."""
        return self.getAvgRegionTable(*args, **kwargs)

    mT = momentT
    getMomentT = momentT
    mCT = momentChainT
    getMomentChainT = momentChainT
    mST = momentStationT
    getMomentStationT = momentStationT
    sT = sensitivityT
    getSensitivityT = sensitivityT
    aCaT = cacheAvgT
    getAvgCacheT = cacheAvgT
    aIT = itemAvgT
    getAvgItemT = itemAvgT
    aOT = orbitAvgT
    getAvgOrbitT = orbitAvgT
    aLT = lossAvgT
    getAvgLossT = lossAvgT
    aRLT = regionLossAvgT
    getAvgRegionLossT = regionLossAvgT
    aRT = regionAvgT
    getAvgRegionT = regionAvgT

    def hasResults(self):
        """True if the solver has computed results."""
        return getattr(self, '_result', None) is not None

    def isSolved(self):
        """True if the solver has computed results (alias of hasResults)."""
        return self.hasResults()

    def getSolverType(self):
        """Get the solver type name (e.g. 'MVA' for SolverMVA)."""
        return self.__class__.__name__.replace('Solver', '')

    def reset(self):
        """Reset solver state and clear results."""
        self._result = None

    @classmethod
    def supportsModel(cls, model):
        """Check if this solver supports the given model type."""
        return True

    # ---- Method-aware feature gating (mirrors MATLAB NetworkSolver) ----
    # A solver declares capability per *method*: the coarse solver-level
    # getFeatureSet()/supports(model) is a union used for pre-checks, while
    # get_method_feature_set(method) yields the per-method envelope actually
    # gated at analysis time. Solvers whose methods do not diverge inherit the
    # base behavior transparently (get_method_feature_set returns None -> the
    # coarse supports(model) is used, preserving any structural checks).

    @staticmethod
    def checkBindingCapacity(model, solver_name):
        """(bool, reason) structural gate for finite station capacity (setCapacity)
        and finite per-class buffers (classCap), used by the product-form solvers
        (MVA, NC). A product-form solver has no representation of a finite buffer,
        so without this gate it silently returns the UNCONSTRAINED answer (e.g.
        QLen=4 instead of the M/M/1/2 value 0.8525). There is no registry feature
        name for plain capacity, hence the structural test. Port of MATLAB
        NetworkSolver.checkBindingCapacity.

        Reads the node-level capacity/class_capacity set by the user, NOT
        sn.cap/sn.classcap: _refresh_capacity derives a FINITE sn.classcap (= the
        chain population) for every closed model, so an sn-level test would reject
        every closed model.

        Only a capacity that can actually BIND is rejected. A closed model whose
        station capacity is at least the total population can never block a job,
        so the declaration is a no-op and the product-form answer stays exact (a
        common idiom: setCapacity(N) on a station of an N-job closed model). njobs
        is Inf for an open class, so any finite capacity reachable by an open
        class binds.

        Cache models are exempt: the Cache node sets class_capacity=1 on the
        retrieval queues it builds, and MVA/NC solve those through their dedicated
        cache/retrieval analyzers rather than as a buffer constraint.
        """
        import numpy as _np
        if not hasattr(model, 'getStruct') or not hasattr(model, '_nodes'):
            return True, ''
        nodes = getattr(model, '_nodes', []) or []
        for node in nodes:
            if type(node).__name__ == 'Cache':
                return True, ''
        njobs = _np.asarray(model.getStruct().njobs, dtype=float).ravel()
        total_jobs = float(_np.sum(njobs))  # inf as soon as one class is open
        for node in nodes:
            tname = type(node).__name__
            # Mirror MATLAB's isa(node,'Station') && ~isa(node,'Source'/'Sink')
            # filter. Today only Stations carry _capacity, but do not rely on
            # that: a future node type with the attribute must not be gated here.
            is_station = any(b.__name__ == 'Station' for b in type(node).__mro__)
            if not is_station or tname in ('Source', 'Sink') or not hasattr(node, '_capacity'):
                continue
            cap = getattr(node, '_capacity', _np.inf)
            if cap is not None and _np.isfinite(cap) and cap >= 0 and cap < total_jobs:
                return False, ("Finite station capacity (setCapacity=%g) at station '%s' is not "
                               "supported by %s. Use SolverCTMC, SolverJMT or SolverLDES."
                               % (cap, node.getName(), solver_name))
            ccap = getattr(node, '_class_capacity', None) or {}
            for jobclass, v in ccap.items():
                if v is None or not _np.isfinite(v) or v <= 0:
                    continue
                if isinstance(jobclass, int):
                    idx = jobclass
                elif hasattr(jobclass, 'get_index0'):
                    idx = jobclass.get_index0()
                else:
                    continue
                if idx is None or idx >= njobs.size:
                    continue
                if v < njobs[idx]:
                    return False, ("Finite per-class capacity (classCap=%g for class %d) at "
                                   "station '%s' is not supported by %s. Use SolverCTMC, "
                                   "SolverJMT or SolverLDES."
                                   % (v, idx + 1, node.getName(), solver_name))
        return True, ''

    check_binding_capacity = checkBindingCapacity

    def resolveMethod(self, options):
        """Resolve the concrete method that will run. Base behavior is a no-op
        (returns options.method). Solvers that perform feature-driven selection
        for options.method='default' override this (typically via selectMethod)."""
        return getattr(options, 'method', 'default')

    def getMethodFeatureSet(self, method):
        """Per-method feature set as a set of feature-name strings, or None to
        signal 'this solver does not diverge per method' (the gate then uses the
        solver's own supports(model)). Divergent solvers override this."""
        return None

    def supportsModelMethod(self, method):
        """Fine, method-aware gate. Returns (bool, reason). Base behavior derives
        the answer from getMethodFeatureSet(method); when that is None, the
        solver's coarse supports(model) is used (reason left empty). Solvers with
        non-feature-set structural per-method rules override this."""
        feat_names = self.getMethodFeatureSet(method)
        if feat_names is None:
            model = getattr(self, 'model', None)
            ok = type(self).supports(model) if model is not None else True
            return bool(ok), ''
        model = getattr(self, 'model', None)
        if model is None:
            return True, ''
        if hasattr(model, 'get_used_lang_features'):
            feat_used = model.get_used_lang_features()
        elif hasattr(model, 'getUsedLangFeatures'):
            feat_used = model.getUsedLangFeatures()
        else:
            return True, ''
        feat_supported = SolverFeatureSet()
        feat_supported.set_true(list(feat_names))
        return SolverFeatureSet.supports_with_reason(feat_supported, feat_used)

    def selectMethod(self, preference_list):
        """Feature-driven selection: first method in preference_list whose
        per-method feature set covers the model; falls back to the last entry."""
        method = preference_list[-1]
        for cand in preference_list:
            ok, _ = self.supportsModelMethod(cand)
            if ok:
                return cand
        return method

    def runAnalyzerChecks(self, options):
        """Single, method-aware feature gate shared by every solver. Resolves the
        concrete method (default may map to a specific method), validates the
        method name, then gates the model against that method's feature set."""
        if not getattr(self, 'enableChecks', True):
            return
        method = self.resolveMethod(options)
        req = getattr(options, 'method', 'default')
        try:
            valid = list(self.listValidMethods())
        except Exception:
            valid = None
        if valid is not None and req not in valid and req != 'default':
            raise RuntimeError("The '%s' method is unsupported by this solver." % req)
        ok, reason = self.supportsModelMethod(method)
        if not ok:
            if method == req:
                raise RuntimeError('This model contains features not supported by the solver. %s' % reason)
            raise RuntimeError("This model contains features not supported by the solver's '%s' method. %s" % (method, reason))

    # snake_case aliases
    resolve_method = resolveMethod
    get_method_feature_set = getMethodFeatureSet
    supports_model_method = supportsModelMethod
    select_method = selectMethod
    run_analyzer_checks = runAnalyzerChecks

    has_results = hasResults
    is_solved = isSolved
    solver_type = getSolverType

    def initFromSolver(self, init_solver):
        """Warm-start the solver from the steady-state solution of an
        auxiliary solver.

        The auxiliary solver's steady-state distribution decides an integer
        job placement (see warmstart.warm_start_placement): with SolverCTMC
        the mode of the exact aggregate stationary distribution, with any
        other solver the rounded mean queue lengths conserving each
        closed-class population. The placement is applied as the model initial
        state via initFromMarginal, which the state-driven solvers honor:
        SolverFLD starts the ODE integration from it, SolverSSA starts the
        simulated trajectory from it, and SolverJMT preloads the stations
        with it. Note that this modifies the initial state of the model
        object shared with any other solver instance.

        Args:
            init_solver: auxiliary solver used to compute the steady-state
                distribution (e.g. SolverCTMC or SolverMVA on the same model)

        Returns:
            self, for chaining
        """
        from .warmstart import warm_start_placement

        model = getattr(self, 'model', None)
        if model is None:
            model = getattr(self, 'network', None)
        if model is None or not hasattr(model, 'initFromMarginal'):
            raise RuntimeError(
                'initFromSolver requires a Network model instance with initFromMarginal')

        sn = model.getStruct()
        placement = warm_start_placement(init_solver, sn)
        model.initFromMarginal(placement)

        # Re-sync any struct cached at construction time: initFromMarginal
        # resets the model struct, so a stale snapshot would miss the newly
        # assigned initial state.
        refreshed = model.getStruct()
        if hasattr(self, '_sn'):
            self._sn = refreshed
        if hasattr(self, 'sn'):
            self.sn = refreshed
        return self

    init_from_solver = initFromSolver

    def getQLen(self):
        """Get average queue lengths (alias for getAvgQLen)."""
        return self.getAvgQLen()

    def getUtil(self):
        """Get utilizations (alias for getAvgUtil)."""
        return self.getAvgUtil()

    def getRespT(self):
        """Get average response times (alias for getAvgRespT)."""
        return self.getAvgRespT()

    def getResidT(self):
        """Get average residence times (alias for getAvgResidT)."""
        return self.getAvgResidT()

    def getTput(self):
        """Get average throughputs (alias for getAvgTput)."""
        return self.getAvgTput()

    def getWaitT(self):
        """Get average waiting times (alias for getAvgWaitT)."""
        return self.getAvgWaitT()

    def getSjrnT(self, *args, **kwargs):
        """System sojourn-time distribution: alias of getCdfRespT (matching the
        JAR/wrapper API)."""
        return self.getCdfRespT(*args, **kwargs)

    sjrn_t = getSjrnT

    def getDistribRespT(self, *args, **kwargs):
        """Response time distributions (same data as getCdfRespT)."""
        return self.getCdfRespT(*args, **kwargs)

    getDistribRespTChain = getDistribRespT
    getDistribRespTNode = getDistribRespT
    getDistribRespTNodeChain = getDistribRespT
    distrib_respt = getDistribRespT
    distrib_respt_chain = getDistribRespT
    distrib_respt_node = getDistribRespT
    distrib_respt_node_chain = getDistribRespT

    def getCacheAvgT(self):
        """Cache metrics table: alias of getAvgCacheTable (JAR-compatible name)."""
        return self.getAvgCacheTable()

    def getItemAvgT(self):
        """Cache item metrics table: alias of getAvgItemTable (JAR-compatible name)."""
        return self.getAvgItemTable()

    cache_avg_t = getCacheAvgT
    item_avg_t = getItemAvgT

    def getAvgRegionTable(self):
        """Per-region finite capacity region (FCR) metrics table.

        One row per region per class with QLen, RespT, ResidT, ArvR, Tput and
        the FCR-specific Weight and MemOcc columns. Populated by solvers whose
        results carry FCR matrices (LDES); empty otherwise, matching the
        wrapper behavior.
        """
        import numpy as np
        import pandas as pd
        self.getAvgTable()
        res = getattr(self, '_ldes_result', None) or getattr(self, '_result', None)
        cols = ['Region', 'JobClass', 'QLen', 'RespT', 'ResidT', 'ArvR', 'Tput',
                'Weight', 'MemOcc']
        rows = []
        QN = getattr(res, 'QNfcr', None) if res is not None else None
        if QN is None and isinstance(res, dict):
            QN = res.get('QNfcr')
        if QN is not None:
            def _get(name):
                v = getattr(res, name, None) if not isinstance(res, dict) else res.get(name)
                return None if v is None else np.asarray(v)
            QN = np.asarray(QN)
            RN = _get('RNfcr'); WN = _get('WNfcr'); TN = _get('TNfcr')
            AN = _get('ANfcr'); WG = _get('WeightNfcr'); MO = _get('MemOccNfcr')
            sn = self._sn if hasattr(self, '_sn') else None
            classnames = (list(sn.classnames) if sn is not None
                          else ['Class%d' % (r + 1) for r in range(QN.shape[1])])

            def _at(mat, f, r):
                return float(mat[f, r]) if mat is not None else float('nan')

            for f in range(QN.shape[0]):
                for r in range(QN.shape[1]):
                    rows.append(['FCRegion%d' % (f + 1), classnames[r],
                                 float(QN[f, r]), _at(RN, f, r), _at(WN, f, r),
                                 _at(AN, f, r), _at(TN, f, r), _at(WG, f, r),
                                 _at(MO, f, r)])
        return pd.DataFrame(rows, columns=cols)

    avg_region_table = getAvgRegionTable
    getRegionAvgT = getAvgRegionTable
    regionAvgT = getAvgRegionTable

    def getDeadlineTable(self):
        """Deadline metrics table (RespT, tardiness, system tardiness).

        Requires a solver whose results carry the TardN/SysTardN matrices;
        returns None otherwise, matching the JAR getDeadlineTable behavior.
        """
        import numpy as np
        import pandas as pd
        self.getAvgTable()
        res = getattr(self, '_result', None)

        def _get(name):
            if res is None:
                return None
            return getattr(res, name, None) if not isinstance(res, dict) else res.get(name)

        TardN = _get('TardN')
        SysTardN = _get('SysTardN')
        RN = _get('RN')
        if RN is None:
            RN = _get('R')
        if TardN is None or SysTardN is None or RN is None:
            return None
        TardN = np.asarray(TardN); SysTardN = np.atleast_2d(np.asarray(SysTardN))
        RN = np.asarray(RN)
        rows = []
        sn = self._sn if hasattr(self, '_sn') else None
        for i in range(RN.shape[0]):
            for k in range(RN.shape[1]):
                if RN[i, k] > 0 or TardN[i, k] > 0 or SysTardN[0, k] > 0:
                    station = (sn.nodenames[int(sn.stationToNode[i])]
                               if sn is not None else 'Station%d' % (i + 1))
                    cls = (sn.classnames[k] if sn is not None else 'Class%d' % (k + 1))
                    rows.append([str(station), str(cls), float(RN[i, k]),
                                 float(TardN[i, k]), float(SysTardN[0, k])])
        if not rows:
            return None
        df = pd.DataFrame(rows, columns=['Station', 'JobClass', 'RespT', 'Trdn', 'SysTrdn'])
        tokeep = ~(df[['RespT', 'Trdn', 'SysTrdn']] <= 0.0).all(axis=1)
        return df.loc[tokeep]

    deadline_table = getDeadlineTable

    def getTranProb(self, node):
        """Transient state probabilities for a node. Only available for solvers
        that compute transient state trajectories (LDES); raises otherwise,
        matching the JAR behavior for unsupported solvers."""
        raise NotImplementedError(
            'getTranProb is not available for %s' % type(self).__name__)

    def getTranProbAggr(self, node):
        """Aggregated transient state probabilities for a node (see getTranProb)."""
        raise NotImplementedError(
            'getTranProbAggr is not available for %s' % type(self).__name__)

    def getTranProbSys(self):
        """Transient system state probabilities (see getTranProb)."""
        raise NotImplementedError(
            'getTranProbSys is not available for %s' % type(self).__name__)

    def getTranProbSysAggr(self):
        """Aggregated transient system state probabilities (see getTranProb)."""
        raise NotImplementedError(
            'getTranProbSysAggr is not available for %s' % type(self).__name__)

    tran_prob = getTranProb
    tran_prob_aggr = getTranProbAggr
    tran_prob_sys = getTranProbSys
    tran_prob_sys_aggr = getTranProbSysAggr

    aT = avgT
    aST = avgSysT
    aNT = avgNodeT
    aCT = avgChainT
    aNCT = avgNodeChainT
    # MATLAB/JAR-canonical word order (chainAvgT vs avgChainT). Both spellings
    # are supported everywhere; get-prefixed forms (getChainAvgT, ...) resolve
    # through alias_getattr. Kept on the base so every NetworkSolver exposes the
    # full set uniformly, superseding the per-solver assignment blocks.
    chainAvgT = avgChainT
    nodeAvgT = avgNodeT
    nodeChainAvgT = avgNodeChainT
    sysAvgT = avgSysT
    avg_t = avgT
    avg_sys_t = avgSysT
    avg_node_t = avgNodeT
    avg_chain_t = avgChainT
    avg_node_chain_t = avgNodeChainT
    a_t = avgT
    a_st = avgSysT
    a_nt = avgNodeT
    a_ct = avgChainT
    a_nct = avgNodeChainT


class EnsembleSolver(Solver):
    """Base class for ensemble/multi-model LINE solvers.

    Used by: SolverLN, SolverUQ, SolverENV
    """

    def isStochastic(self):
        """An ensemble solver is stochastic if any of its submodel solvers is
        stochastic. Each submodel solver classifies itself, including from the
        method it resolved at runtime.
        """
        for solver in getattr(self, 'solvers', None) or []:
            if solver is not None and hasattr(solver, 'isStochastic') and solver.isStochastic():
                return True
        return False

    is_stochastic = isStochastic

    def avgT(self):
        """Short alias for getAvgTable (resolved on the concrete ensemble
        solver: SolverLN/SolverENV/SolverUQ)."""
        return self.getAvgTable()

    aT = avgT
    getAvgT = avgT
    avg_t = avgT
    a_t = avgT
