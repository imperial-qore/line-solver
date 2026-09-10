"""Structural admissibility of the NC methods.

ONE PREDICATE, TWO CALLERS. :func:`nc_method_refusal` answers "may this method
run on this model?" once, and both callers ask it: ``SolverNC.runAnalyzer``
turns a non-empty answer into an error ahead of the dispatch, and
``SolverNC.supportsModelMethod`` uses it so that ``Network.findSolver`` and
``Network.help`` never offer a (solver, method) pair that would raise, and so
that ``SolverAUTO`` never delegates to one. Two copies of these rules is
precisely how the report and the run drift apart, which is the failure this
module exists to prevent: a new rule goes here, not at a call site.

ONLY WHAT THE FEATURE REGISTRY CANNOT NAME LIVES HERE. A feature set declares
what a method ACCEPTS, so it can refuse a model for HAVING a construct and
never for lacking one: "closed population only" and "no think time" are
expressed in ``SolverNC.getMethodFeatureSet`` by dropping ``OpenClass`` and
``SchedStrategy_INF``, while "requires a cache", "requires state-dependent
routing", "requires a loss network", "requires exactly two stations" and
"requires normal usage" have no such form and are decided here.

Mirrors MATLAB ``nc_method_refusal.m`` / ``nc_is_lossn_model.m`` /
``nc_is_normal_usage.m``, the JAR ``SolverNC.ncMethodRefusal`` and the C++
``nc::nc_method_refusal``.
"""

import numpy as np

from ...api.sn.network_struct import NodeType
from ...lang.base import DropStrategy, ReplacementStrategy

__all__ = ['nc_method_refusal', 'nc_is_lossn_model', 'nc_is_normal_usage']


def _nodetypes(sn):
    if sn is None or getattr(sn, 'nodetype', None) is None:
        return np.zeros(0, dtype=int)
    return np.ravel(np.asarray(sn.nodetype, dtype=int))


def _has_nodetype(sn, ty) -> bool:
    return bool(np.any(_nodetypes(sn) == int(ty)))


def _param(nodeparam, name, default=None):
    """A nodeparam field, whichever of the two shapes the struct carries."""
    if nodeparam is None:
        return default
    if isinstance(nodeparam, dict):
        return nodeparam.get(name, default)
    return getattr(nodeparam, name, default)


def _is_slotted(options) -> bool:
    cfg = getattr(options, 'config', None) if options is not None else None
    if cfg is None and isinstance(options, dict):
        cfg = options.get('config')
    if not cfg:
        return False
    slotted = cfg.get('slotted') if isinstance(cfg, dict) else getattr(cfg, 'slotted', None)
    return bool(slotted)


def nc_is_lossn_model(sn):
    """(is_lossn, has_shape) for the loss network solver_nc_lossn_analyzer solves.

    The shape is an OPEN model with a single finite capacity region whose only
    member is an infinite server; it is a LOSS network when the admission rule
    DROPS every class. The two answers are returned separately because they have
    different remedies: a region of that shape under WAITQ (or any blocking
    rule) holds the arrival back instead of discarding it, which keeps the job in
    the region while it waits and is a queueing phenomenon the Erlang loss model
    has no state for, so switching the rule to DROP makes it solvable here.
    """
    if sn is None:
        return False, False
    try:
        from ...api.sn.predicates import sn_has_closed_classes
        if int(getattr(sn, 'nregions', 0) or 0) != 1:
            return False, False
        if sn_has_closed_classes(sn):
            return False, False
        region = getattr(sn, 'region', None)
        if region is None or len(region) == 0 or region[0] is None:
            return False, False
        region_matrix = np.asarray(region[0], dtype=float)
        K = int(sn.nclasses)
        members = []
        for i in range(int(sn.nstations)):
            per_class = bool(np.any(region_matrix[i, :K] >= 0))
            aggregate = bool(region_matrix[i, K] >= 0) if region_matrix.shape[1] > K else False
            if per_class or aggregate:
                members.append(i)
        if len(members) != 1 or not np.isinf(np.asarray(sn.nservers, dtype=float).flatten()[members[0]]):
            return False, False
        rule = getattr(sn, 'regionrule', None)
        if rule is None:
            return False, True
        rule = np.asarray(rule, dtype=float)
        # ALL classes, not merely one: a region that discards one class and holds
        # another back is a mixed system whose blocked class occupies the region
        # while it waits, so the per-class loss probabilities the Erlang fixed
        # point returns would not be the ones the model implies.
        for r in range(K):
            if rule[0, r] != float(DropStrategy.DROP):
                return False, True
        return True, True
    except Exception:
        return False, False


def nc_is_normal_usage(sn) -> bool:
    """Is the closed model in NORMAL USAGE, the domain of the Mitra-McKenna
    PANACEA asymptotic expansion (J. ACM 33(3), 1986)?

    Normal usage asks that every queueing centre absorb the load the think
    stations offer it: with rho_j0 = Ztot(j) the aggregate think demand of chain
    j, r_ij = L_ij / rho_j0 and mu_i(Ntot) the saturation rate,

        alpha_i = 1 - (sum_j N_j r_ij) / mu_i(Ntot) > 0     for every centre i.

    Outside it the {phi(n)} series DIVERGES, which is why pfqn_panaceald returns
    NaN there and pfqn_ncld turns that NaN into a refusal rather than a warning.
    It is a property of the demands and not of a declared construct, so it has
    no feature-registry name and cannot live in a feature set.

    The rates are the ones solver_ncld would build: 1 for an ordinary single
    server, min(n, c) for a finite multiserver (the conversion runAnalyzer
    performs on the 'panald' arm), and the declared lldscaling row when the
    model sets one. An infinite server is a think station and feeds Ztot.
    """
    if sn is None:
        return True
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    if np.any(np.isinf(njobs)):
        # An open or mixed chain is refused earlier by the closed-population
        # feature set of the load-dependent evaluators; there is no rho_j0.
        return True
    from ...api.sn.demands import sn_get_demands_chain
    dem = sn_get_demands_chain(sn)
    Lchain = np.asarray(dem.Lchain, dtype=float)
    Nchain = np.asarray(dem.Nchain, dtype=float).flatten()
    Nt = int(round(float(np.sum(Nchain[np.isfinite(Nchain)]))))
    if Nt < 1:
        return True  # the empty network: G = 1, nothing to expand

    M = int(sn.nstations)
    nservers = np.asarray(sn.nservers, dtype=float).flatten()
    lld = sn.lldscaling
    if lld is None or np.size(lld) == 0:
        lld = np.ones((M, Nt))
        for i in range(M):
            if np.isfinite(nservers[i]) and nservers[i] > 1:
                lld[i, :] = np.minimum(np.arange(1, Nt + 1), nservers[i])
    else:
        lld = np.asarray(lld, dtype=float)
        if lld.shape[1] < Nt:
            pad = np.repeat(lld[:, -1:], Nt - lld.shape[1], axis=1)
            lld = np.hstack([lld, pad])

    is_is = np.isinf(nservers[:M])
    Ztot = np.sum(Lchain[is_is, :], axis=0) if np.any(is_is) else np.zeros(Lchain.shape[1])
    if np.any((Nchain > 0) & (Ztot <= 0)):
        # no think station on the route of a populated chain: the expansion
        # parameter rho_j0 is undefined
        return False

    Lq = Lchain[~is_is, :]
    if Lq.shape[0] == 0:
        return True  # no queueing centre: the delay-only constant is exact
    muq = lld[~is_is, :Nt]
    if np.any(muq <= 0) or np.any(~np.isfinite(muq)):
        return False

    r = np.zeros_like(Lq)
    for j in range(Ztot.size):
        if Ztot[j] > 0:
            r[:, j] = Lq[:, j] / Ztot[j]
    alpha = 1.0 - (r @ Nchain) / muq[:, Nt - 1]
    return bool(np.min(alpha) > 0)


def _closed_queueing_stations(sn) -> int:
    """How many queueing (non-infinite-server) stations carry demand from a
    CLOSED chain.

    That is the row count L reaches pfqn_nc and pfqn_comomrm_ld with, once the
    delay rows have been folded into Z and the zero-demand rows dropped. Zero
    when the model has no closed population at all.
    """
    from ...api.sn.demands import sn_get_demands_chain
    dem = sn_get_demands_chain(sn)
    Lchain = np.asarray(dem.Lchain, dtype=float)
    Nchain = np.asarray(dem.Nchain, dtype=float).flatten()
    closed = np.isfinite(Nchain) & (Nchain > 0)
    if not np.any(closed):
        return 0
    nservers = np.asarray(sn.nservers, dtype=float).flatten()
    nq = 0
    for i in range(int(sn.nstations)):
        if np.isinf(nservers[i]):
            continue
        if np.any(np.abs(Lchain[i, closed]) > 1e-8):
            nq += 1
    return nq


def _is_noreentrant_cache(sn) -> bool:
    """The Source-Cache-Sink model solver_nc_cache_analyzer serves."""
    nt = _nodetypes(sn)
    if int(getattr(sn, 'nclosedjobs', 0) or 0) != 0 or nt.size != 3:
        return False
    return sorted(nt.tolist()) == sorted([int(NodeType.SOURCE), int(NodeType.CACHE),
                                          int(NodeType.SINK)])


def nc_method_refusal(sn, method, options=None, for_report=True) -> str:
    """'' when METHOD may run on SN, else the reason it may not, in the words
    the analyzer refuses with.

    ``for_report`` says WHICH QUESTION IS BEING ASKED, and for two method names
    the two questions have different answers:

    ``True``
        "should ``model.help()`` offer this pair?" A pair that comes back as a
        table of zeros must not be offered, so the answer is no.
    ``False``
        "what does the reference DO when asked for it by name?" For 'mmint2' and
        'gleint' outside their shape the reference deliberately WARNS AND RETURNS
        A ZERO TABLE (pfqn_nc.m, case {'mmint2','gleint'}: lG = [] and return,
        unconditionally), and a caller who names the method keeps that answer.

    THE ASYMMETRY IS A RULING, NOT AN OVERSIGHT (2026-07-25, reaffirmed when this
    gate was added): the report answers "should this be offered" and the run
    answers "what does the reference do". 'comomld' is NOT in that bucket --
    pfqn_comomrm_ld raises 'The solver accepts at most a single queueing
    station.' natively -- so it is refused on both paths.
    """
    from ...api.sn.predicates import (sn_has_dps, sn_has_open_classes, sn_is_mm1k_loss)
    from .solver_nc_dps_analyzer import nc_is_dps_model
    from .solver_nc_oi_analyzer import nc_is_oi_model

    if sn is None:
        return ''
    m = str(method or 'default').lower()

    # The discrete-time route answers for itself: nc_is_dt_model decides
    # admissibility on the slot lattice, and every gate below is written about a
    # continuous-time queueing network.
    if _is_slotted(options):
        return ''

    # -- discriminatory processor sharing ---------------------------------
    # Morrison's heavy-usage expansion is the ONLY NC route that can see the DPS
    # weights; every other method builds a product-form normalizing constant
    # that silently drops them and answers with the egalitarian-PS network,
    # which is a wrong number rather than a coarse one.
    if nc_is_dps_model(sn):
        if m not in ('default', 'morrison'):
            return ("SolverNC: method %s cannot represent the DPS weights of a discriminatory "
                    "processor-sharing station; it would return the egalitarian-PS network. Use "
                    "method 'default' or 'morrison' (npfqn_dps_morrison), SolverMVA, SolverFLD or "
                    "SolverCTMC." % method)
        return ''
    if sn_has_dps(sn):
        # A DPS station outside Morrison's shape. SchedStrategy_DPS is declared
        # in the feature set because a boolean feature cannot express "this shape
        # only"; this is that imperative half.
        return ("SolverNC analyzes a discriminatory processor-sharing station only in the shape "
                "Morrison's expansion is derived for: a CLOSED network of exactly two stations, one "
                "infinite-server (think) station and one single-server DPS station, exponential "
                "service, each class visiting the two equally often. Use SolverMVA, SolverFLD or "
                "SolverCTMC for any other DPS model.")
    if m == 'morrison':
        # The method named on a model that is not the shape at all -- not even a
        # DPS station in it. Left ungated it reaches no route of its own and
        # falls through to the ordinary normalizing-constant path, which would
        # answer the product-form model UNDER THE CALLER'S LABEL.
        return ("SolverNC: method 'morrison' is the heavy-usage expansion of a CLOSED network of exactly two "
                "stations, one infinite-server (think) station and one single-server DPS station "
                "with exponential service, which this model is not. Remove the method option to let "
                "SolverNC choose, or use SolverMVA, SolverFLD or SolverCTMC.")

    # -- Krzesinski state-dependent routing -------------------------------
    # An SDR model is intercepted by solver_nc_sdr_analyzer whatever the method
    # says, so reaching the second test means the model declares none.
    if getattr(sn, 'sdr', None):
        return ''
    if m in ('sdr', 'sdr.mva'):
        return ("SolverNC: method %s requires state-dependent routing, which this model does not declare."
                % method)

    # -- stochastic Petri net ---------------------------------------------
    # A net is served only by the MDD-rec route, and none of the gates below --
    # written about stations, capacities and the queueing-network product form --
    # says anything about a net. spn_pf decides its product-form class, by name.
    if _has_nodetype(sn, NodeType.PLACE):
        if m not in ('default', 'rec'):
            return ("SolverNC: a stochastic Petri net is solved by the MDD-rec route; method '%s' is a "
                    "normalizing-constant algorithm for queueing networks. Use 'rec' or 'default'"
                    % method)
        return ''

    # -- order-independent stations ---------------------------------------
    # Every method other than the four listed reads sn.rates, which holds only
    # the single-job rate mu([r]) of an OI station: the rank rate mu(n) is
    # silently dropped and the answer is that of an ordinary queue.
    if nc_is_oi_model(sn):
        if m not in ('default', 'exact', 'is', 'sampling'):
            return ("SolverNC: method '%s' cannot represent the rank rate mu(n) of an order-independent "
                    "station; use method 'default' or 'exact' (pfqn_ncoi), 'is', SolverMVA, or "
                    "SolverCTMC." % method)
        return ''

    # -- caches ------------------------------------------------------------
    # 'rayint' and 'spm' both name the SPM saddle point of a cache (and, on a
    # retrieval model, the ray/WKB delayed-hit expansion), so they are
    # admissible here and nowhere else.
    if _has_nodetype(sn, NodeType.CACHE):
        if m == 'exact' and _is_noreentrant_cache(sn):
            ci = int(np.flatnonzero(_nodetypes(sn) == int(NodeType.CACHE))[0])
            nodeparam = getattr(sn, 'nodeparam', None)
            ch = nodeparam[ci] if nodeparam is not None and len(nodeparam) > ci else None
            rs = _param(ch, 'replacestrat', _param(ch, 'replacement', None))
            # cache_prob_erec is exact for the exchangeable (RR/FIFO) family
            # only; anything else has to take the approximate route.
            if rs is not None and int(rs) not in (int(ReplacementStrategy.RR),
                                                  int(ReplacementStrategy.FIFO)):
                return ("SolverNC: NC does not support exact solution of the specified cache replacement "
                        "policy; use the default (approximate) method or SolverCTMC.")
        return ''
    if m in ('rayint', 'spm'):
        return ("SolverNC: method %s names the SPM saddle point of a cache and, on a retrieval model, the "
                "ray/WKB delayed-hit expansion; this model declares no Cache node." % method)

    # -- single-station M/M/1/K with tail drop -----------------------------
    # Answered exactly by the probability-based qsys_mm1k_loss branch, which
    # reads no method name at all, so no method name gate below applies to it.
    try:
        if sn_is_mm1k_loss(sn):
            return ''
    except Exception:
        pass

    # -- loss networks and finite capacity regions -------------------------
    is_lossn, has_lossn_shape = nc_is_lossn_model(sn)
    if is_lossn:
        return ''  # 'ms', 'erlangfp', 'rec' and 'default' all have a route here
    if has_lossn_shape:
        return ("SolverNC does not support finite capacity regions with WAITQ (blocking) policy. "
                "Use DROP policy instead.")
    if m in ('ms', 'erlangfp'):
        return ("SolverNC: method %s is admissible only on a loss network (open model, one DROP region "
                "holding a single Delay)." % method)
    if m == 'rec':
        return ("SolverNC: method rec is the MDD-rec route, admissible on a stochastic Petri net or on a "
                "loss network (open model, one DROP region holding a single Delay); this model is "
                "neither.")
    if int(getattr(sn, 'nregions', 0) or 0) > 0:
        # NC does not enforce an aggregate region limit on queueing stations;
        # refuse rather than silently return the unconstrained answer.
        return ("SolverNC: this model uses a Finite Capacity Region (addRegion) on queueing "
                "stations, which is not supported by SolverNC. Use SolverJMT, or setCapacity for "
                "a single-station limit.")

    # -- PANACEA's domain ---------------------------------------------------
    # Normal usage is a property of the demands rather than of a declared
    # construct, so it has no feature name; an open chain is refused earlier by
    # the closed-population feature set of the load-dependent evaluators.
    #
    # BOTH TOKENS ARE GATED, because pfqn_ncld evaluates 'pana' and
    # 'panald' with the SAME pfqn_panaceald -- its case label is
    # ('pana', 'panald') -- so on a model carrying a rate lattice the
    # load-INDEPENDENT name reaches the load-dependent expansion and raises with
    # it. Off that lattice 'pana' takes its own pfqn_nc arm, which warns and
    # returns an empty constant rather than raising, so it is left alone there.
    # Class- or joint-dependent scaling diverts the whole model to
    # solver_nc_conv, which never reads the method at all.
    if m in ('pana', 'panald') and not sn_has_open_classes(sn):
        diverted_to_conv = (bool(getattr(sn, 'cdscaling', None))
                            or bool(getattr(sn, 'jdscaling', None)))
        lld = getattr(sn, 'lldscaling', None)
        reaches_ld_kernel = m == 'panald' or (lld is not None and np.size(lld) > 0)
        if not diverted_to_conv and reaches_ld_kernel and not nc_is_normal_usage(sn):
            reason = ("The model is not in normal usage, so the 'panald' asymptotic expansion "
                      "does not apply. Use 'exact', 'clw' or an approximate load-dependent method "
                      "instead.")
            if m == 'pana':
                reason = ("Method 'pana' reaches the load-dependent kernel on this model, "
                          "where pfqn_ncld evaluates it as 'panald'. " + reason)
            return reason

    # -- the single-queueing-station recursions ------------------------------
    # Two families are stated for a model with a delay and ONE queueing station,
    # and neither can say so with a feature name: it is a COUNT, and a feature
    # set has no arithmetic. pfqn_nc states it for 'mmint2'/'gleint' in those
    # words and pfqn_comomrm_ld raises 'The solver accepts at most a single
    # queueing station.'
    #
    # The count is taken over the CLOSED chains only, and the rule is inactive
    # without a closed population, because pfqn_nc answers an open network with
    # the exact open formulas BEFORE its method switch -- the method name is never read
    # there, so a purely open model with three queues runs these names correctly
    # today and must go on doing so.
    # 'mmint2' and 'gleint' are gated for the REPORT ONLY: pfqn_nc answers them
    # with an empty constant and the caller renders a table of zeros, which is a
    # pair the report must not offer and a run the reference nonetheless
    # performs. See for_report above.
    if m == 'comomld' or (for_report and m in ('mmint2', 'gleint')):
        nq = _closed_queueing_stations(sn)
        if nq > 1:
            if m == 'comomld':
                return ("SolverNC: method 'comomld' is the load-dependent CoMoM recursion, and "
                        "pfqn_comomrm_ld accepts at most a single queueing station; this model "
                        "has %d." % nq)
            return ("SolverNC: the '%s' method requires a model with a delay and a single "
                    "queueing station; this model has %d." % (method, nq))

    # -- 'exact' outside its domain ----------------------------------------
    if m == 'exact':
        nservers = np.asarray(sn.nservers, dtype=float).flatten()
        njobs = np.asarray(sn.njobs, dtype=float).flatten()
        multiserver = bool(np.any(nservers[np.isfinite(nservers)] > 1))
        if multiserver and np.any(np.isinf(njobs)):
            return ("SolverNC: NC solver cannot provide exact solutions for open or mixed queueing "
                    "networks. Remove the 'exact' option.")
        scaling = ((sn.lldscaling is not None and np.size(sn.lldscaling) > 0)
                   or bool(getattr(sn, 'cdscaling', None))
                   or bool(getattr(sn, 'jdscaling', None)))
        finite = njobs[np.isfinite(njobs)]
        fractional = bool(finite.size and np.any(np.abs(finite - np.floor(finite)) > 1e-12))
        if (scaling or multiserver) and fractional:
            # The load-dependent analyzer interpolates a fractional population
            # between the two integer neighbours, which is an approximation, so
            # it refuses the exactness the caller asked for by name.
            return ("SolverNC: NC load-dependent solver cannot provide exact solutions for fractional "
                    "populations.")
    return ''
