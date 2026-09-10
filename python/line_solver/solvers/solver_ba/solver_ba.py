"""Native SolverBA: bound-analysis solver.

Subclasses the native SolverMVA to reuse its model parsing (rates, demands,
njobs, nservers, sched, chain info), and replaces runAnalyzer with the
bound-analysis dispatch in solver_ba_analyzer. Mirrors matlab/src/solvers/BA.
"""

import numpy as np

from ..solver_mva.solver_mva import SolverMVA as _NativeMVA
from .solver_ba_analyzer import (
    BA_AUTO, BA_BGT, BA_BPT, BA_SNC, BA_SPNLP,
    ba_blocking_default, ba_ignores_blocking, ba_method_degenerate,
    ba_method_refusal, ba_resolve_method, ba_resolve_model_method,
    solver_ba_analyzer, BA_ASYM, BA_CHAIN, BA_HIER, BA_HAREL, BA_LR,
    BA_MAPAMVA, BA_QRF_LP, BA_QRF_NLP)

# QRF reduction bounds are served under SolverBA in MATLAB/JAR, and now
# natively as well. Both families route through solver_ctmc_qrf_analyzer:
#   qrf.bas, qrf.rsrd            -- LP, via the qr_bounds_* backends (HiGHS).
#   qr, qrf.mmi, qrf.mem,        -- NLP, via api.mapqn.qrf_noblo_*.
#   qrf.mmi.ld, qrf.mmi.linear
# The NLP tokens were withheld until 2026-07-20 pending end-to-end coverage of
# the adapter path. Exercising it found a transposed v in the adapter and a
# throughput inversion valid only for a delay reference station; with both
# fixed, K == 1 reproduces the CTMC exactly and K == 2 lands inside the glpsol
# bound range. Infinite-server stations are rejected: qrf_noblo_* models every
# station as a single server.
BA_METHODS = ['default'] + sorted(
    BA_ASYM | BA_CHAIN | BA_HIER | BA_HAREL | BA_LR | BA_MAPAMVA | BA_QRF_LP
    | BA_QRF_NLP | BA_AUTO | BA_BPT | BA_BGT | BA_SNC | BA_SPNLP)


class SolverBA(_NativeMVA):
    """Bound-analysis solver (asymptotic + hierarchical throughput/queue-length
    bounds, plus the QRF LP-based bounds).  cub is upper-only; mbjb.lower is the
    multiclass BJB lower bound (Kerola eq. 10).  scb (Dowdy et al. 1992) is the
    only family that does NOT bracket this model's own solution: it brackets the
    multiclass system the given single-class model aggregates, and is kept out
    of auto.* for that reason.

    FINITE-BUFFER BLOCKING IS REFUSED, NOT BOUNDED. Needing only demands and a
    population is the BCMP parameterization, which presumes UNBOUNDED buffers;
    a buffer that binds couples the station occupancies and the resulting
    numbers do not bracket the blocked model. runAnalyzer therefore gates on
    sn_has_blocking and list_valid_methods drops every blocking-blind method.
    The exceptions are qrf.bas*/qrf.rsrd, which carry the blocking tables
    explicitly. Use SolverMVA with method 'sqd' for a point estimate."""

    def __init__(self, model, method_or_options=None, **kwargs):
        level = kwargs.pop('level', None)
        # level is not a SolverMVAOptions field, so it can't ride through kwargs; read it off an options object here or hierarchy silently reverts to level 2
        if level is None and method_or_options is not None and not isinstance(method_or_options, str):
            if hasattr(method_or_options, 'get'):
                level = method_or_options.get('level', None)
            else:
                level = getattr(method_or_options, 'level', None)
        super().__init__(model, method_or_options, **kwargs)
        if level is not None:
            setattr(self.options, 'level', int(level))
        elif not hasattr(self.options, 'level'):
            setattr(self.options, 'level', 2)

    def getName(self):
        return 'SolverBA'

    @staticmethod
    def getFeatureSet():
        """The feature envelope of the bounds.

        SolverBA DECLARED NO FEATURE SET AT ALL, so it inherited the base
        supports(), which accepts everything: every model in the language passed,
        including the ones whose constructs the analyzer has no representation
        of. MATLAB had the mirror-image defect, a set naming no service
        distribution and therefore refusing everything; both are fixed to this
        list.

        WHAT MAKES A DISTRIBUTION ADMISSIBLE HERE IS ITS MEAN. The analyzer reads
        sn.rates and sn.visits and nothing else: every bound in the
        ABA/BJB/PB/GB/SB/Harel/MWBA families is a function of the demands
        D = V/rates and the think time, so any renewal law with a finite mean is
        admissible whatever its higher moments. The QRF reduction is the one that
        needs more, and what it needs is a PH representation (the (D0, D1) pair
        out of sn.proc), which the phase-type families below carry.

THE MODULATED LAWS 'MAP' AND 'MMPP2' ARE IN THE BASE ENVELOPE FOR ONE
        FAMILY, 'mapamva', AND getMethodFeatureSet STRIPS THEM FROM EVERY OTHER.
        They were out entirely until mapamva landed, on the correct ground that a
        renewal bound derived for a product-form network says nothing about a
        correlated one: its mean rate exists, so the utilization law still holds
        and the formula still returns a number, but that number brackets a
        DIFFERENT system. MAP-AMVA (Casale-Smirni, DSN 2009) is derived FOR the
        correlated model -- its variables are the per-phase QN(i,k) and UN(i,k)
        -- so the same reasoning that refuses the others admits it.

        The direction matters and is forced: a feature set refuses a model for
        HAVING a construct and never for lacking one, so the only way to grant a
        law to one family is to put it in the base envelope and take it away from
        the rest. MMAP and BMAP stay OUT everywhere -- the LP is single-class and
        has no marked or batch arrival variable -- and Cache and Fork/Join stay
        out for the older reason, that no analyzer here has a representation of
        either.

        THE PETRI-NET CONSTRUCTS ARE IN, since solver_ba_spnlp. They were out on
        exactly the same "no representation" ground and that ground is gone: the
        spnlp relaxation is indexed by the marking, reads the enabling,
        inhibiting and firing arcs out of sn.nodeparam, and refuses by name the
        modes it cannot carry. QueueingPlace stays OUT and is refused by name: a
        place with an embedded queue has local state the relaxation has no
        variable for. Same division SolverNC draws.

        Mirrors SolverBA.getFeatureSet in MATLAB and the JAR.
        """
        return {
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue',
            'Sink', 'Source', 'Router',
            'StatelessClassSwitcher',
            'ClosedClass', 'OpenClass',
            # A self-looping class is a closed chain of one station, which the
            # demand-parameterized bounds read as any other chain. MultiServer
            # and FiniteCapacity are deliberately NOT here: getMethodFeatureSet
            # grants each to the families that carry it.
            'SelfLoopingClass',
            # renewal service laws: the bounds need the mean, the QRF reduction
            # needs the PH form, and sn.proc carries both
            'APH', 'Coxian', 'Cox2', 'Erlang', 'Exp', 'HyperExp', 'PH',
            'Det', 'Lognormal', 'Pareto', 'Uniform', 'Weibull',
            # modulated service, for 'mapamva' alone; see the note above
            'MAP', 'MMPP2',
            'SchedStrategy_INF', 'SchedStrategy_PS',
            'SchedStrategy_FCFS', 'SchedStrategy_LCFSPR',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            # Petri-net constructs, for the spnlp family. QueueingPlace is
            # deliberately absent; see the note above.
            'Place', 'Transition', 'Linkage', 'Enabling', 'Inhibiting',
            'Timing', 'Firing', 'Storage',
        }

    get_feature_set = getFeatureSet

    def getMethodFeatureSet(self, method):
        """The envelope the METHOD GATE compares the model against: SolverBA's
        own, never SolverMVA's.

        THIS CLASS SUBCLASSES SolverMVA FOR ITS MODEL PARSING and inherited the
        MVA method gate with it, which is a different solver's answer to a
        question about this one. SolverMVA declares Cache, CacheClassSwitcher,
        the replacement strategies, Fork/Join and MMAP; SolverBA declares none
        of them and says why in getFeatureSet (the bounds are derived for a
        product-form network and the analyzer has no representation of a cache
        or of a fork-join). The gate is what findSolver, listValidMethods and
        SolverAUTO's ranked choice all consult, so the inherited set made every
        one of the 36 bound methods look runnable on a cache model, where each
        of them refuses at run time.

        MATLAB, the JAR and C++ all put SolverBA on NetworkSolver rather than on
        SolverMVA and so never had the defect; this override is what brings the
        native python answer back to theirs (a 3-class LRU cache model reported
        117 runnable pairs here against 78 in MATLAB and 76 in C++).

        ON TOP OF THAT COME THE PER-METHOD DELTAS the registry CAN name. A
        feature set says "I accept this construct", so it can refuse a model for
        HAVING one and never for lacking one; that is exactly the shape of the
        delay-station and open-class premises below, and exactly not the shape
        of "one class" or "one server", which have no feature name and live in
        ba_method_refusal instead.

        DELAY STATIONS. 'sb' and 'lr' reject an infinite-server station
        outright, and 'harel', 'sib' and 'scb' reject a nonzero think time,
        which on these models is the same station: harel extrapolates the exact
        normalizing constant of a delay-free network, SIB Section 3.2 is the
        extension that would carry Z and is not implemented, and SCB Theorem 3
        rests on the delay-free balanced-network throughput. The three OPEN
        families reject one too, each being derived for one server per station.

        CLASS TYPES. The three OPEN families drop ClosedClass, which is the
        whole of their class premise. The MIRROR delta -- dropping OpenClass
        from every demand-parameterized family -- is deliberately NOT applied:
        "supports single-class closed networks only" is one rule, its
        single-class half has no feature name, and splitting it across the two
        mechanisms would report the closed half here and the single-class half
        in ba_method_refusal for the same model. It is stated once, structurally.
        'spnlp' takes no delta at all: it is indexed by the marking, and whether
        that marking is bounded is a question about the P-invariants of the net,
        which spn_lpbnd answers.
        """
        feats = set(SolverBA.getFeatureSet())
        # Judged on the RESOLVED name so that 'default' carries the envelope of
        # the gb.upper it runs as, and MODEL-AWARE so that on a blocked model it
        # carries the envelope of the qrf.bas it is routed to.
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_struct'):
            resolved, _why = ba_resolve_model_method(model.get_struct(), method)
        else:
            resolved = ba_resolve_method(method)
        fam = resolved.split('.')[0]
        if fam in ('bpt', 'bgt', 'snc'):
            feats -= {'ClosedClass', 'Delay', 'DelayStation', 'SchedStrategy_INF'}
            if fam != 'snc':
                # SERVICE AND ARRIVAL LAWS. 'bpt' and 'bgt' are derived for a
                # MARKOVIAN open network and read the mean alone, so a
                # non-exponential law anywhere is not something they refuse at
                # run time -- it is something they silently bound as if it were
                # Poisson. Measured on the M/M/1 shape: replacing the Exp(1)
                # source by an Erlang of the same mean leaves bgt.upper at QLen
                # 32.6667 and bpt.lower at 1.0, digit for digit. That is a bound
                # on a DIFFERENT system, so the laws are dropped here rather
                # than left to a run-time check the analyzers do not make: their
                # own procid test covers the queueing stations only and would
                # miss exactly the source case.
                #
                # 'snc' is excluded: it CONSUMES the arrival law (the same
                # substitution moves it from 3.8244 to 3.0092) and its analyzer
                # branches on a non-exponential source deliberately. Its rule is
                # about the SERVICE only, which no feature name can say, so it
                # lives in ba_method_refusal instead.
                feats -= {'APH', 'Coxian', 'Cox2', 'Erlang', 'HyperExp', 'PH',
                          'Det', 'Lognormal', 'Pareto', 'Uniform', 'Weibull'}
        elif fam in ('sb', 'harel', 'sib', 'scb', 'lr'):
            feats -= {'Delay', 'DelayStation', 'SchedStrategy_INF'}
        # MODULATED SERVICE, THE ONE DELTA THAT RUNS THE OTHER WAY. 'MAP' and
        # 'MMPP2' sit in the base envelope so that 'mapamva' can accept them,
        # which means every OTHER family has to give them back: each of the rest
        # reads the service MEAN alone and would bound a correlated model as
        # though its services were independent, quietly returning a bracket for a
        # different system. Written as a delta on the complement rather than as a
        # grant because a feature set can refuse a model for HAVING a construct
        # and never for lacking one.
        #
        # 'mapamva' takes the delay delta instead: its LP is a network of queues,
        # and Casale-Smirni name the delay extension as open work.
        if fam == 'mapamva':
            feats -= {'Delay', 'DelayStation', 'SchedStrategy_INF'}
        else:
            feats -= {'MAP', 'MMPP2'}
        # MULTISERVER (registry name since 2026-09-05) is out of the base
        # envelope: every demand-parameterized family reads one server per
        # station (ba_method_refusal names 'ssd' as the alternative), the
        # alpha-free QRF arms refuse a c-server station through sn_to_qrf_alpha
        # and the open families through their own refusal. What carries the
        # count is granted here: 'ssd' (the multiserver bound), 'ldbcmp' (its
        # fixed-rate form runs on the c-server rate law), 'auto' (which picks
        # among them) and the two load-dependent QRF arms, whose alpha(i,n) IS
        # min(n,c).
        if fam in ('ssd', 'ldbcmp', 'auto') or resolved in ('qrf.mmi.ld', 'qrf.mmi.linear'):
            feats.add('MultiServer')
        # FINITECAPACITY (registry name since 2026-09-05): only the QRF blocking
        # bounds carry the buffer (the MM, MM1, ZZ, ZM, BB, F tables), which is
        # the same split ba_ignores_blocking makes; the structural refusal keeps
        # naming them. 'default'/'auto'/'auto.upper' resolve to 'qrf.bas' on a
        # blocked model of the right shape, so they are granted it through
        # RESOLVED. 'spnlp' is NOT granted: its polytope reads no Place capacity,
        # so a capped place would be relaxed away.
        if resolved.startswith('qrf.bas') or resolved.startswith('qrf.rsrd'):
            feats.add('FiniteCapacity')
        return feats

    get_method_feature_set = getMethodFeatureSet

    def supportsModelMethod(self, method):
        """The base gate, NOT SolverMVA's.

        Skipping SolverMVA in the chain is deliberate and is not only about the
        feature set above. SolverMVA.supportsModelMethod adds two rules of its
        own, and both are wrong here: supportsFiniteCapacity refuses a model
        with a binding buffer, which SolverBA answers for through
        list_valid_methods instead (it keeps qrf.bas*/qrf.rsrd, the two families
        that carry the blocking tables explicitly, and drops the rest), and
        supportsExactness gates a method named 'exact', which no bound method
        is. Mirrors the MATLAB, JAR and C++ SolverBA, none of which inherits an
        MVA gate.

        ON TOP OF IT come the structural premises no feature name can express.
        ba_method_refusal is the same predicate solver_ba_analyzer raises on, so
        a caller gets one answer whichever of the two it meets first -- which is
        the point: this gate is what findSolver, list_valid_methods and
        SolverAUTO's ranked choice read, and while it was silent about them a
        two-class closed network was reported as able to run 36 bound methods of
        which 30 raised on contact.

        STRUCTURAL PREDICATE FIRST, THEN THE FEATURE ENVELOPE, the same order as
        MATLAB @SolverBA/supportsModelMethod and for the same reason: the run
        raises ba_method_refusal's sentence, so a model that fails BOTH must be
        reported with that sentence and not with the generic feature one. The
        order became visible when 'MultiServer' entered the registry: a c-server
        model then failed the envelope too, and the gate answered "feature:
        MultiServer" where the analyzer says which method to use instead.
        """
        from ..base import NetworkSolver
        struct_reason = ba_method_refusal(self.model.get_struct(), method)
        if struct_reason:
            return False, struct_reason
        ok, reason = NetworkSolver.supportsModelMethod(self, method)
        if not ok:
            return ok, reason
        # A bound that APPLIES but says nothing is not offered either; see
        # ba_method_degenerate for why that is a separate question and why the
        # analyzer is still allowed to answer it when asked by name.
        degen_reason = ba_method_degenerate(self.model.get_struct(), method)
        if degen_reason:
            return False, degen_reason
        return True, reason

    supports_model_method = supportsModelMethod

    @staticmethod
    def supports(model) -> bool:
        """Whether SolverBA can bound this model.

        Finite-buffer blocking is NOT decided here: there is no registry feature
        name for a capacity, so a blocked model passes this gate and is handled
        where it can be told apart -- list_valid_methods drops the blocking-blind
        bounds and keeps qrf.bas*/qrf.rsrd, which carry the blocking tables
        explicitly.
        """
        from ..base import supports_via_featureset
        return supports_via_featureset(SolverBA, model)

    @staticmethod
    def list_all_methods():
        """Every bound method the solver implements, independently of the model."""
        return list(BA_METHODS)

    listAllMethods = list_all_methods

    def list_valid_methods(self):
        """Bound methods this MODEL can run: a narrowing of list_all_methods.

        A name that would always be refused on this model is not offered, so a
        caller enumerating the list never asks for one. Ask for it by name
        anyway and runAnalyzer still dispatches, so the analyzer's own reason is
        what comes back."""
        methods = list(BA_METHODS)
        # The STRUCTURAL premises -- single-class closed, fully closed,
        # single-server -- come from ba_method_refusal, the same predicate
        # solver_ba_analyzer raises on and supportsModelMethod reports. Asked
        # first so every later narrowing works on names this model could
        # actually run: before it, a two-class closed network was offered all 36
        # demand-parameterized bounds and 30 of them raised on contact.
        # ... and ba_method_degenerate withholds the second kind of name: one
        # whose premises this model MEETS but whose formula says nothing here.
        # 'ldbcmp.lower' at N == Qhat is the only such case: it reports the
        # trivial X >= 0, which propagates into an all-zero table a caller
        # cannot tell from an answer. Offering is what stops; asking for it by
        # name still runs, since a vacuous bound is a valid one.
        sn0 = self.model.get_struct()
        methods = [m for m in methods
                   if not ba_method_refusal(sn0, m) and not ba_method_degenerate(sn0, m)]
        # The QR/LR/QRF reduction bounds share one premise: a single-class
        # closed network of single-server stations (solver_ba_qrf_analyzer
        # gates on it, and the lr family refuses a delay). Naming them on a
        # model they cannot run turns a rejection into a method a caller is
        # invited to ask for. Mirrors SolverBA.m and SolverBA.java.
        # The LOAD-DEPENDENT arms survive where the rest of the family cannot
        # run: alpha(i,n) is the rate law of a delay (alpha = n), of a c-server
        # station (alpha = min(n,c)) and of limited load dependence alike, so
        # 'qrf.mmi.ld' and 'qrf.mmi.linear' answer those models on the model's
        # own chain. sn_to_qrf_alpha owns the one restriction that survives,
        # exponential service wherever a station serves several jobs at once.
        # Dropping them with the rest would hide from a caller the only two
        # bound methods such a model has.
        # 'mapamva' shares that premise exactly -- single-class closed,
        # single-server, no delay -- so it narrows with them rather than being
        # offered on a model it would refuse on contact. It is NOT one of the
        # load-dependent arms: its q carries no population index, so a delay or a
        # c-server station has nowhere to go.
        if not self._is_reducible():
            ld_arms = ('qrf.mmi.ld', 'qrf.mmi.linear') if self._is_ld_reducible() else ()
            methods = [m for m in methods
                       if m in ld_arms
                       or (m not in ('qr', 'lr')
                           and not m.startswith('lr.') and not m.startswith('qrf.')
                           and not m.startswith('mapamva'))]
        # 'bpt', 'bgt' and 'snc' are the mirror image: all three are derived for
        # an OPEN network of single-server exponential stations, so every closed
        # model, every delay station and every multiserver station rules them
        # out. Every other family rules OUT the open model, so on an open
        # network the list narrows to those three. 'bgt.upper' additionally
        # needs deterministic non-merging routes and 'snc.upper' a feed-forward
        # station graph, both of which their analyzers check by walking the
        # routing matrix -- too expensive to repeat here, so they stay listed
        # and refuse by name.
        if not self._is_bpt_feasible():
            methods = [m for m in methods
                       if m not in ('bpt.lower', 'bgt.upper', 'snc.upper')]
        if self._is_fully_open():
            methods = [m for m in methods
                       if m in ('bpt.lower', 'bgt.upper', 'snc.upper')]
        # 'spnlp.*' is the only family indexed by a MARKING rather than by
        # demands and a population, and the split is total in both directions:
        # on a Petri net nothing else has a representation of the model, and off
        # one spnlp has nothing to read. Two of the gates above already
        # half-cover this by accident -- a Place is an INF station, so
        # _is_reducible and _is_bpt_feasible are both false on any Petri net --
        # but the demand-parameterized families survive them and must be
        # dropped by name.
        if self._is_petri():
            methods = [m for m in methods if m.startswith('spnlp')]
        else:
            methods = [m for m in methods if not m.startswith('spnlp')]
        # A binding finite buffer rules out everything but the QRF blocking
        # bounds: the other families presume unbounded buffers, and runAnalyzer
        # refuses them by name on such a model. The list can legitimately come
        # back EMPTY -- a blocked model that is not single-class closed
        # single-server has no bound method at all, and offering one would be
        # the mis-selection this gate exists to prevent.
        if self._has_blocking():
            methods = [m for m in methods
                       if not ba_ignores_blocking(ba_resolve_method(m))]
            # 'default' is offered back when it now MEANS one of the survivors:
            # runAnalyzer routes it to 'qrf.bas' on a blocked model of the right
            # shape, so a caller enumerating the list would otherwise be told the
            # model's own default is invalid.
            if ba_blocking_default(self.model.get_struct())[0]:
                methods = ['default'] + methods
        return methods

    def _is_petri(self):
        """Whether the model holds Transition nodes, i.e. is a Petri net."""
        from ...api.sn.network_struct import NodeType
        sn = self.model.get_struct()
        nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
        return bool(np.any(nodetype == int(NodeType.TRANSITION)))

    def _is_fully_open(self):
        """Whether every class of the model is open."""
        sn = self.model.get_struct()
        return bool(np.all(np.isinf(np.asarray(sn.njobs, dtype=float))))

    def _is_bpt_feasible(self):
        """Whether 'bpt'/'bgt' apply: fully open, no INF station, one server each."""
        # network_struct.SchedStrategy, NOT constants.SchedStrategy: sn.sched
        # holds the former. The latter is a plain Enum while this one is an
        # IntEnum, and an Enum member equals nothing but itself, so comparing
        # across them is ALWAYS FALSE -- the INF test below was dead code,
        # masked only because a Delay also has nservers = inf and the next
        # check catches it. See the note in _is_reducible on why the two enums
        # cannot simply be unified.
        from ...api.sn.network_struct import SchedStrategy
        if not self._is_fully_open():
            return False
        sn = self.model.get_struct()
        nservers = np.asarray(sn.nservers, dtype=float).ravel()
        for i in range(int(sn.nstations)):
            if sn.sched[i] == SchedStrategy.INF:
                return False
            if nservers[i] > 1:
                return False
        return True

    def _is_reducible(self):
        """Whether the QR/LR/QRF reduction applies: one closed class, no INF, c=1."""
        # network_struct.SchedStrategy, NOT constants.SchedStrategy -- see
        # _is_bpt_feasible. Importing the right one is the fix here rather than
        # making constants.SchedStrategy an IntEnum, which looks like the root
        # fix and is NOT safe: seven of its members carry different values from
        # the other two definitions (constants inserts PSJF/FB/LAS/LRPT before
        # FSP/PAS/OI, shifting those three by +4), so value-based comparison
        # would turn always-False into wrongly-True. Unifying the numbering is
        # a cross-codebase interchange change, since sn.sched ids travel to
        # MATLAB, the JAR and C++.
        from ...api.sn.network_struct import SchedStrategy
        sn = self.model.get_struct()
        if int(sn.nclasses) != 1:
            return False
        if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
            return False
        # sn.sched is a {station index -> SchedStrategy} dict here, not an array
        nservers = np.asarray(sn.nservers, dtype=float).ravel()
        for i in range(int(sn.nstations)):
            if sn.sched[i] == SchedStrategy.INF:
                return False
            if nservers[i] > 1:
                return False
        return True

    def _is_ld_reducible(self):
        """Whether the LOAD-DEPENDENT reduction applies where `_is_reducible`
        says no: one closed class, and a scaling alpha(i,n) that sn_to_qrf_alpha
        can derive for every station."""
        from ...api.sn.sn_to_qrf_alpha import sn_to_qrf_alpha
        sn = self.model.get_struct()
        if int(sn.nclasses) != 1:
            return False
        if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
            return False
        _alpha, msg, _ld, _peak = sn_to_qrf_alpha(sn)
        return not msg

    listValidMethods = list_valid_methods

    def _has_blocking(self):
        """Whether a finite buffer or capacity region BINDS anywhere."""
        from ...api.sn import sn_has_blocking
        return bool(sn_has_blocking(self.model.get_struct()))

    def _check_blocking(self, method):
        """Resolve a method against a model with a binding finite buffer.

        Returns the method to actually run, which is METHOD itself in every
        ordinary case, and raises when the request is blocking-blind and cannot
        be routed.

        Finite-buffer blocking is outside the premises of every family but the
        QRF blocking bounds: the rest are parameterized by demands and a
        population alone, which presumes unbounded buffers and a product form
        the truncation destroys. Refusing is not conservatism -- on
        cqn_bas_blocking, gb.upper reports QLen 1.28 at a station capped at 1
        job. But a blocked model whose shape admits 'qrf.bas' gets it as the
        DEFAULT instead of a refusal, since that bound models the finite buffer
        and derives its own tables. Mirrors MATLAB @SolverBA/runAnalyzer.m.
        """
        if not ba_ignores_blocking(ba_resolve_method(method)):
            return method
        if not self._has_blocking():
            return method
        why = ''
        if method in ('default', 'auto', 'auto.upper'):
            routed, why = ba_blocking_default(self.model.get_struct())
            if routed:
                return routed
        detail = ('' if not why else
                  ' The QRF blocking bounds do not apply here either: %s' % why)
        raise ValueError(
            "Method '%s' does not support finite-buffer blocking: every SolverBA bound "
            "family but the QRF blocking ones is parameterized by demands and a population "
            "alone, so it bounds the model as if its buffers were unbounded. Use SolverMVA "
            "with method 'sqd', an exact solver (CTMC, SSA, JMT, LDES), or the QRF blocking "
            "bounds 'qrf.bas'/'qrf.rsrd', which model the finite buffer.%s"
            % (method, detail))

    def runAnalyzer(self):
        # Gated before the lang dispatch so lang='cpp' is refused on the same
        # terms as the native path, and after the aliases so 'default' is
        # judged as the gb.upper it resolves to.
        method = self._check_blocking(str(getattr(self, 'method', 'default')).lower())
        # lang='cpp' delegates to the C++ line-cli (-s ba); an absent binary is the
        # only automatic fallback, a bound family the port lacks propagates by name.
        # There is no lang='java' counterpart: the JAR CLI has no `ba` token.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, populate_cpp_result
            try:
                populate_cpp_result(self)
                return self
            except LineCliNotAvailable as e:
                from ...api.io.logging import line_warning
                line_warning("SolverBA", "lang='cpp' requested but the C++ solver is "
                             "unavailable (%s); falling back to lang='python'." % e)
        # METHOD is what _check_blocking resolved to, which is the request itself
        # unless a blocked model routed 'default'/'auto' onto 'qrf.bas'.
        self._result = solver_ba_analyzer(self, method, self.options)
        return self

    run_analyzer = runAnalyzer

    def get_bounds(self):
        """Return {Tlower,Tupper,Qlower,Qupper} for the current method's family.
        One-sided families (cub upper-only, mbjb/ldbcmp lower-only) return NaN on
        the missing side.

        Each side is re-run through self.options, so the caller's full option set
        (notably `level`) is inherited: a hierarchical family reached through
        get_bounds must tighten as level is raised, not silently revert to the
        default level of 2."""
        fam = str(getattr(self, 'method', 'default')).lower().split('.')[0]
        # A blocking-blind family on a blocked model is dropped from
        # list_valid_methods, so both sides below would silently come back NaN.
        # Refuse by name instead, with the reason runAnalyzer gives. Judged on
        # the full method, not on fam: the family prefix of 'qrf.bas' is the
        # bare 'qrf', which is not a blocking bound.
        # runAnalyzer routes 'default'/'auto' to 'qrf.bas' on a blocked model, and
        # get_bounds has to say so rather than contradict it -- but it still
        # cannot BRACKET, because the analyzer solves qrf.bas in the 'max'
        # direction alone. So the routed case is refused on its own terms.
        _requested = str(getattr(self, 'method', 'default')).lower()
        _resolved = self._check_blocking(_requested)
        if _resolved != _requested:
            raise ValueError(
                "'%s' resolves to '%s' on this model, which has a binding finite buffer, "
                "and that bound is UPPER-only: there is no bracket to return. Call "
                "getAvgTable/getAvg for the upper bound, or SolverMVA with method 'sqd' "
                "for a point estimate." % (_requested, _resolved))
        valid = self.list_valid_methods()
        out = {'Tlower': np.nan, 'Tupper': np.nan,
               'Qlower': np.nan, 'Qupper': np.nan}
        for side, kt, kq in (('lower', 'Tlower', 'Qlower'),
                             ('upper', 'Tupper', 'Qupper')):
            m = fam + '.' + side
            if m in valid:
                r = solver_ba_analyzer(self, m, self.options)
                out[kt] = r['TN']
                out[kq] = r['QN']
        return out

    getBounds = get_bounds

    def _expand_bound(self, x, M, K):
        """Normalize a bracket side to an (M,K) array. A one-sided family leaves
        its missing side as the scalar NaN returned by get_bounds; expand it so
        the table keeps full shape with NaN entries."""
        a = np.asarray(x, dtype=float)
        if a.size == 0:
            return np.full((M, K), np.nan)
        if a.ndim == 0:
            return np.full((M, K), float(a))
        return a.reshape(M, K) if a.size == M * K else a

    def get_bounds_table(self, keepDisabled=False):
        """Return the {lower,upper} bracket per station and class as a DataFrame,
        in the layout of getAvgTable. Columns:
        Station, JobClass, Qlower, Qupper, Tlower, Tupper.

        One-sided families (cub upper-only, mbjb/ldbcmp lower-only) carry NaN on
        the missing side; NaN is preserved, never replaced by zero."""
        import pandas as pd
        b = self.get_bounds()
        M = self.nstations
        K = self.nclasses
        Ql = self._expand_bound(b['Qlower'], M, K)
        Qu = self._expand_bound(b['Qupper'], M, K)
        Tl = self._expand_bound(b['Tlower'], M, K)
        Tu = self._expand_bound(b['Tupper'], M, K)
        station_names = list(getattr(self, 'station_names', []))
        class_names = list(getattr(self, 'class_names', []))
        rows = []
        for i in range(M):
            sname = station_names[i] if i < len(station_names) else 'Station%d' % i
            for k in range(K):
                cname = class_names[k] if k < len(class_names) else 'Class%d' % k
                vals = np.array([Ql[i, k], Qu[i, k], Tl[i, k], Tu[i, k]], dtype=float)
                present = vals[~np.isnan(vals)]
                # Mirror getAvgTable's drop of disabled station-class pairs, but
                # NaN-safe: keep the row when any value that is present is
                # nonzero, so an all-NaN side never removes the row.
                if keepDisabled or present.size == 0 or np.any(present != 0):
                    rows.append({
                        'Station': sname,
                        'JobClass': cname,
                        'Qlower': Ql[i, k],
                        'Qupper': Qu[i, k],
                        'Tlower': Tl[i, k],
                        'Tupper': Tu[i, k],
                    })
        return pd.DataFrame(
            rows,
            columns=['Station', 'JobClass', 'Qlower', 'Qupper', 'Tlower', 'Tupper'])

    getBoundsTable = get_bounds_table

    def _snc_perc(self, eps):
        """Both quantile matrices from one envelope propagation.

        The envelope builder returns the (arrival, service) handles it built, so
        the two quantiles cost one pass over the network and one Chernoff search
        per pair and metric."""
        from ...api.snc import snc_perc_delay, snc_perc_backlog
        from .solver_ba_snc import snc_envelopes
        sn = self._sn if getattr(self, '_sn', None) is not None \
            else self.model.get_struct()
        env = snc_envelopes(sn)
        M, K = env['M'], env['K']
        D = np.full((M, K), np.nan)
        B = np.full((M, K), np.nan)
        for (i, r), arv in env['arv'].items():
            srv = env['srv'][(i, r)]
            D[i, r] = snc_perc_delay(arv, srv, eps)[0]
            B[i, r] = snc_perc_backlog(arv, srv, eps)[0]
        return D, B

    def get_delay_perc(self, eps=1e-3):
        """Per-station response-time QUANTILE at violation probability eps: the
        smallest d for which P{D_ir > d} <= eps is certified.

        This is the native output of the 'snc' family, so the accessor runs the
        SNC envelope propagation directly whatever options.method says; a
        station-class pair carrying no traffic stays NaN. Every other family
        bounds means only and has no counterpart."""
        return self._snc_perc(eps)[0]

    getDelayPerc = get_delay_perc

    def get_backlog_perc(self, eps=1e-3):
        """Per-station queue-length QUANTILE, in jobs, at violation probability
        eps. The counterpart of get_delay_perc; see it for the conventions."""
        return self._snc_perc(eps)[1]

    getBacklogPerc = get_backlog_perc

    def get_perc_table(self, eps=1e-3):
        """Response-time and queue-length quantiles at violation probability eps
        as a DataFrame, in the layout of getAvgTable. Columns: Station,
        JobClass, RespTPerc, QLenPerc. Rows are the pairs that carry traffic."""
        import pandas as pd
        D, B = self._snc_perc(eps)
        station_names = list(getattr(self, 'station_names', []))
        class_names = list(getattr(self, 'class_names', []))
        rows = []
        for i in range(D.shape[0]):
            sname = station_names[i] if i < len(station_names) else 'Station%d' % i
            for k in range(D.shape[1]):
                if np.isnan(D[i, k]) and np.isnan(B[i, k]):
                    continue
                cname = class_names[k] if k < len(class_names) else 'Class%d' % k
                rows.append({'Station': sname, 'JobClass': cname,
                             'RespTPerc': D[i, k], 'QLenPerc': B[i, k]})
        return pd.DataFrame(rows, columns=['Station', 'JobClass',
                                           'RespTPerc', 'QLenPerc'])

    getPercTable = get_perc_table
