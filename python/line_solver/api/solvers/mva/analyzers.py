"""
MVA Solver analyzers.

Native Python implementation of MVA solver analyzers that orchestrate
method selection and provide the main entry point for MVA analysis.

Port from:

"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional
import time

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    sn_has_product_form,
    sn_has_fractional_populations,
    sn_has_bursty_arrival,
    sn_rt_stations,
)
from ...sn.network_struct import NodeType
from ...mam import map_lambda, map_idc, map_pie, map_count_idc, map_scv
from ...qsys import qsys_gig1_rq, qsys_gigk_rqt, qsys_gigk_rqt_gamma
from ...npfqn import npfqn_traffic_idc, npfqn_traffic_split_rr, npfqn_traffic_rqt
from .handler import (solver_mva, SolverMVAOptions, SolverMVAReturn,
                      mva_supports_single_class_open, mva_supports_qna_scheduling)
from ...io.logging import line_warning_always


@dataclass
class MVAResult:
    """
    Result of MVA solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        AN: Arrival rates (M x K)
        WN: Waiting times (M x K)
        logNormConstAggr: Log normalizing constant
        iter: Number of iterations
        runtime: Runtime in seconds
        method: Method used
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    AN: Optional[np.ndarray] = None
    WN: Optional[np.ndarray] = None
    logNormConstAggr: float = np.nan
    iter: int = 0
    runtime: float = 0.0
    method: str = ""


def _has_lcfs_lcfspr_config(sn: NetworkStruct) -> bool:
    """
    Check if network has special LCFS/LCFS-PR configuration.

    This configuration has a product-form solution but is not
    recognized by sn_has_product_form.

    Args:
        sn: Network structure

    Returns:
        True if network has both LCFS and LCFS-PR stations
    """
    has_lcfs = False
    has_lcfspr = False

    sched_dict = sn.sched if sn.sched else {}
    for i in range(sn.nstations):
        station_sched = sched_dict.get(i)

        if station_sched == SchedStrategy.LCFS:
            has_lcfs = True
        elif station_sched == SchedStrategy.LCFSPR:
            has_lcfspr = True

    return has_lcfs and has_lcfspr


def _is_mixed_exact_model(sn: NetworkStruct) -> bool:
    """
    True for a mixed open+closed product-form model whose finite-server
    stations are all single-server and whose closed populations are integral.
    Such models are solved exactly by BCMP mixed MVA (pfqn_mvams via
    solver_mva); the AMVA linearizer inflates the closed-class response time
    by double-counting the open-class interference. Port of the mixed-exact
    dispatch condition in matlab/src/solvers/MVA/solver_mva_analyzer.m.
    """
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if not (np.any(np.isinf(njobs)) and np.any(np.isfinite(njobs) & (njobs > 0))):
        return False
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    finite_ns = nservers[np.isfinite(nservers)]
    if finite_ns.size == 0 or np.max(finite_ns) != 1:
        return False
    closed = njobs[np.isfinite(njobs)]
    if not np.all(closed == np.floor(closed)):
        return False
    return sn_has_product_form(sn)


def _has_classdep_routing(sn: NetworkStruct) -> bool:
    """
    True when classes are not routed alike, i.e. the routing probabilities
    differ between job classes -- either because a class switches class on a hop
    or because two classes leave the same station with different probabilities.
    False for a single-class model and for a multiclass model in which every
    class traverses the network identically.

    This is the condition under which per-class visit ratios diverge, so a
    method that aggregates classes into a per-chain demand vector stops being
    exact. Used to gate Marie's aggregation-decomposition in SolverMVA. Port of
    matlab/src/api/sn/sn_has_classdep_routing.m; sn.rt is station-major,
    (i-1)*K+r, matching sn_refresh_visits.
    """
    K = int(sn.nclasses)
    M = int(sn.nstations)
    if K <= 1 or sn.rt is None:
        return False
    rt = np.asarray(sn.rt, dtype=float)
    if rt.shape[0] < M * K or rt.shape[1] < M * K:
        return False
    tol = 1e-8
    for i in range(M):
        for j in range(M):
            shared = None
            for r in range(K):
                # Class switching on a hop makes the routing class-dependent
                # outright.
                for t in range(K):
                    if r != t and rt[i * K + r, j * K + t] > tol:
                        return True
                p = rt[i * K + r, j * K + r]
                if shared is None:
                    shared = p
                elif abs(p - shared) > tol:
                    return True
    return False


def _is_bas_model(sn: NetworkStruct) -> bool:
    """
    Detect a closed single-chain network with Blocking-After-Service (BAS) finite-buffer
    blocking, which the BAS approximation handles but exact/AMVA MVA does not.

    sn.droprule is an (M, K) int array of DropStrategy codes (BAS == 2).
    """
    if sn.nchains != 1 or sn.nclosedjobs <= 0 or sn.droprule is None:
        return False
    njobs = np.asarray(sn.njobs).flatten()
    if np.any(np.isinf(njobs)):
        return False  # open class present
    from ...sn.network_struct import DropStrategy
    return bool(np.any(np.asarray(sn.droprule) == int(DropStrategy.BAS)))


def _amva_uses_pf_kernels(sn: NetworkStruct) -> bool:
    """
    Conservative form of the branch test in solver_amva: True when this model may be solved
    by the product-form AMVA kernels (the linearizer family and relatives) rather than by
    solver_amvald. The mixed case is reported True for every resolved method, while the
    branch itself takes it only for 'lin', so a caller is never told that solver_amvald will
    run when it might not.
    """
    from ...sn import (sn_has_product_form_not_het_fcfs, sn_has_load_dependence,
                       sn_has_open_classes)
    has_scaling = ((sn.cdscaling is not None and np.size(sn.cdscaling) > 0) or
                   (getattr(sn, 'jdscaling', None) is not None and np.size(sn.jdscaling) > 0))
    return (sn_has_product_form_not_het_fcfs(sn)
            and not sn_has_load_dependence(sn)
            and not has_scaling
            and (not sn_has_open_classes(sn) or sn_has_product_form(sn)))


def _il_needs_amva(sn: NetworkStruct, options=None) -> bool:
    """
    True when an interlock matrix is supplied but exact MVA cannot honour it. pfqn_mva
    carries the Eq. (4.7) correction for closed single-server models only, so a mixed or
    multiserver model with an interlock goes to AMVA, which applies the same correction to
    the arrival-instant queue length.
    """
    IL = getattr(options, 'interlock', None)
    if IL is None or np.size(IL) == 0:
        return False
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    if np.any(np.isinf(njobs)):
        return True
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    finite_servers = nservers[np.isfinite(nservers)]
    return bool(finite_servers.size > 0 and np.max(finite_servers) > 1.0)


def mva_carries_interlock(sn: NetworkStruct, options=None) -> bool:
    """
    True when the MVA path this model already dispatches to carries a class-level interlock
    matrix (Franks 1999, Eq. 4.7) itself, so that supplying one does not silently move the
    model to a DIFFERENT algorithm.

    Only two kernels implement the correction: pfqn_mva (exact, closed single-server) and the
    AMVA forward step of solver_amvald. A model that would otherwise be solved by exact
    multiserver or mixed MVA, or by the product-form AMVA kernels, cannot take the matrix
    without swapping its algorithm, and the swap is worth far more than the correction it
    carries: inside SolverLN it can turn a converging Picard iteration into a limit cycle. A
    caller holding a matrix such a model cannot carry must apply its own correction instead.
    """
    method = getattr(options, 'method', None) or 'default'
    method = str(method).lower()
    if method.startswith('amva.'):
        method = method[5:]

    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    has_open = bool(np.any(np.isinf(njobs)))
    has_closed = bool(np.any(np.isfinite(njobs) & (njobs > 0)))
    closed_pops = njobs[np.isfinite(njobs)]
    closed_pops_integral = bool(np.all(closed_pops == np.floor(closed_pops)))
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    finite_servers = nservers[np.isfinite(nservers)]
    has_finite_server = finite_servers.size > 0
    max_finite_servers = float(np.max(finite_servers)) if has_finite_server else None

    # pfqn_mva takes the matrix for a closed single-server model, and for nothing else
    pfqn_mva_can_take_it = (not has_open and closed_pops_integral
                            and (max_finite_servers is None or max_finite_servers <= 1.0))

    amva_methods = ('amva', 'bs', 'qd', 'qli', 'fli', 'lin', 'qdlin', 'sqni', 'gflin',
                    'egflin', 'ab', 'schmidt', 'schmidt-ext', 'tay', 'scat', 'aql', 'qsa',
                    'lcp', 'chow', 'pamb', 'pami', 'pamt', 'clust', 'dmlin',
                    'priomva')

    if method in ('exact', 'mva'):
        return pfqn_mva_can_take_it
    if method in amva_methods:
        # Only solver_amvald carries the correction among the AMVA handlers
        return not _amva_uses_pf_kernels(sn)
    if method == 'default':
        if _is_bas_model(sn):
            return False  # solver_sqd has no interlock term
        exact_mixed = (has_open and has_closed and has_finite_server
                       and max_finite_servers == 1.0 and sn_has_product_form(sn)
                       and closed_pops_integral)
        exact_small = (sn.nchains <= 4 and float(np.sum(njobs)) <= 20
                       and sn_has_product_form(sn)
                       and not sn_has_fractional_populations(sn))
        if exact_mixed or exact_small:
            return pfqn_mva_can_take_it
        return not _amva_uses_pf_kernels(sn)
    # mvac, sqd, sum, qna, rqna and rqt reach neither kernel
    return False


def solver_sqd(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Blocking-After-Service (BAS) approximate MVA handler.

    Wraps npfqn_sqd (chain-aggregated approximation): solves the single closed chain,
    then disaggregates the chain-level throughput, queue length and utilization back to
    per-class results via sn_deaggregate_chain_results. Supports single-chain closed
    networks; multi-chain models warn and return empty.
    """
    import time
    from ...npfqn import npfqn_sqd
    from ...sn import sn_get_demands_chain, sn_deaggregate_chain_results
    from ....api.io.logging import line_warning

    start_time = time.time()
    if options is None:
        options = SolverMVAOptions()

    if sn.nchains != 1:
        line_warning("solver_sqd",
                     "SQD (Smith Queue Decomposition) supports single-chain closed networks only; "
                     "this model is multichain (nchains=%d) - returning empty results." % sn.nchains)
        M0, K0 = sn.nstations, sn.nclasses
        return SolverMVAReturn(Q=np.full((M0, K0), np.nan), U=np.full((M0, K0), np.nan),
                               R=np.full((M0, K0), np.nan), T=np.full((M0, K0), np.nan),
                               C=np.full((1, K0), np.nan), X=np.full((1, K0), np.nan),
                               lG=np.nan, runtime=time.time() - start_time, method='sqd', it=0)

    dem = sn_get_demands_chain(sn)
    Lchain = dem.Lchain
    STchain = dem.STchain
    Vchain = dem.Vchain
    alpha = dem.alpha
    refstatchain = dem.refstatchain

    M = sn.nstations
    bas = npfqn_sqd(sn, sn.nclosedjobs)

    # Assemble chain-level (M x 1) matrices for the single closed chain.
    refstat = int(np.asarray(refstatchain).flatten()[0])
    Xchain = np.zeros((1, 1))
    Tchain = np.zeros((M, 1))
    Qchain = np.zeros((M, 1))
    Uchain = np.zeros((M, 1))
    Rchain = np.zeros((M, 1))
    # Vchain is normalized to 1 at the reference station, so the per-station
    # throughput there equals the chain reference throughput.
    Xchain[0, 0] = bas.X[refstat]
    for i in range(M):
        Tchain[i, 0] = bas.X[i]
        Qchain[i, 0] = bas.Q[i]
        Uchain[i, 0] = bas.U[i]
        Rchain[i, 0] = bas.R[i]

    deagg = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha,
        Qchain, Uchain, Rchain, Tchain, None, Xchain
    )

    return SolverMVAReturn(
        Q=deagg.Q, U=deagg.U, R=deagg.R, T=deagg.T, C=deagg.C, X=deagg.X,
        lG=np.nan, runtime=time.time() - start_time, method='sqd', it=0
    )


def solver_mva_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> MVAResult:
    """
    MVA Analyzer - main entry point for MVA analysis.

    Selects appropriate MVA method based on network characteristics
    and solver options, then performs the analysis.

    Supported methods:
        - 'exact', 'mva': Exact MVA using population recursion
        - 'default': Automatic method selection
        - 'amva', 'bs', 'qd', 'sqni', etc.: Approximate MVA methods
        - 'qna': Queueing Network Analyzer (for general networks)

    Args:
        sn: Network structure
        options: Solver options (method, tolerance, verbosity)

    Returns:
        MVAResult with all performance metrics

    Raises:
        RuntimeError: For unsupported methods or configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    method = options.method.lower()

    # Remove 'amva.' prefix if present
    if method.startswith('amva.'):
        method = method[5:]

    result = MVAResult()
    ret = None

    if method in ['exact', 'mva']:
        ret = solver_mva(sn, options)
        result.iter = 0

    elif method == 'sqd':
        ret = solver_sqd(sn, options)

    elif method == 'qna':
        # QNA for general networks
        ret = solver_qna(sn, options)

    elif method == 'rqna':
        # Robust Queueing Network Analyzer (indices of dispersion)
        ret = solver_rqna(sn, options)

    elif method == 'default':
        # Bursty single-class open network: the arrival process is non-renewal
        # (MMPP/MAP), whose autocorrelation a two-moment method cannot capture.
        # Dispatch RQNA, which characterizes each flow by its IDC.
        if (sn.nclasses == 1 and np.all(np.isinf(sn.njobs))
                and sn_has_bursty_arrival(sn)):
            ret = solver_rqna(sn, options)
            method = 'rqna'
        # Automatic method selection
        elif _is_bas_model(sn):
            # Closed single-chain network with Blocking-After-Service finite buffers
            ret = solver_sqd(sn, options)
            method = 'sqd'
        elif _il_needs_amva(sn, options):
            # An interlock matrix that exact MVA cannot honour sends the model to AMVA,
            # which applies the same Eq. (4.7) correction to the arrival-instant queue
            # length: pfqn_mva carries it for closed single-server models only. Mirrors
            # matlab/src/solvers/MVA/solver_mva_analyzer.m and the JAR analyzer.
            ret = solver_amva(sn, options)
            method = ret.method if hasattr(ret, 'method') else 'amva'
        elif _has_lcfs_lcfspr_config(sn):
            # LCFS/LCFS-PR configuration has product-form solution
            ret = solver_mva(sn, options)
            method = 'exact'
        elif _is_mixed_exact_model(sn):
            # Mixed open+closed product-form single-server model: exact BCMP
            # mixed MVA (pfqn_mvams via solver_mva). The AMVA linearizer inflates
            # the closed-class response time by double-counting the open-class
            # interference. Only the finite (closed) populations need to be
            # integral (isinf fails the fractional-population test spuriously).
            # Mirrors matlab/src/solvers/MVA/solver_mva_analyzer.m.
            ret = solver_mva(sn, options)
            method = 'exact'
        elif (sn.nchains <= 4 and
              np.sum(sn.njobs) <= 20 and
              sn_has_product_form(sn) and
              not sn_has_fractional_populations(sn)):
            # Small network with product-form - use exact MVA
            ret = solver_mva(sn, options)
            method = 'exact'
        else:
            # Use approximate MVA
            ret = solver_amva(sn, options)
            method = ret.method if hasattr(ret, 'method') else 'amva'

    elif method in ['amva', 'bs', 'qd', 'qli', 'fli', 'lin',
                    'qdlin', 'sqni', 'gflin', 'egflin',
                    'ab', 'schmidt', 'schmidt-ext']:
        # ab, schmidt and schmidt-ext were missing from this list while
        # listValidMethods advertised them, so asking for one fell through to
        # the branch below and returned metrics of None -- the same "false
        # claim" the bound methods were removed for. They are dispatched by
        # MATLAB (solver_mva_analyzer.m) and by the JAR, so they are dispatched
        # here.
        ret = solver_amva(sn, options)
        method = ret.method if hasattr(ret, 'method') else method

    else:
        # An unsupported method must RAISE by name. Returning a result whose
        # every metric is None left the caller to propagate the None into its
        # own arithmetic, so the failure surfaced far from its cause -- and with
        # verbose off there was no message at all. Matches the line_error in
        # matlab/src/solvers/MVA/solver_mva_analyzer.m and the throw in the JAR
        # analyzer.
        raise ValueError(
            f"solver_mva_analyzer: the '{method}' method is not dispatched by "
            "solver_mva_analyzer. Supported: default, exact, mva, amva, bs, qd, qli, fli, "
            "lin, qdlin, sqni, egflin, gflin, ab, schmidt, schmidt-ext.")

    # Copy results from handler return
    if ret is not None:
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.CN = ret.C
        result.XN = ret.X
        result.AN = ret.T.copy() if ret.T is not None else None  # Arrival rates = throughputs
        result.WN = ret.R.copy() if ret.R is not None else None  # Waiting times ≈ response times
        result.logNormConstAggr = ret.lG
        result.iter = ret.it

    result.runtime = time.time() - start_time
    result.method = method

    # Report exhaustion of the iteration budget. The AMVA fixed-point loops stop
    # on `norm(Q-Qlast)<tol or iter > maxiter` and do NOT distinguish the two
    # cases, so a non-converged result is otherwise returned silently and is
    # indistinguishable from a converged one. The metrics may still be usable --
    # they often stabilise long before the tolerance is met -- but that is not
    # something the caller can verify unless told.
    #
    # The check lives at this single exit rather than inside solver_amva, which
    # returns from many points; a tail check there would miss most of them,
    # including the 'lin' path.
    _iter_max = getattr(options, 'iter_max', None)
    if _iter_max and result.iter and result.iter >= _iter_max:
        line_warning_always(
            'solver_mva_analyzer',
            "AMVA method '%s' exhausted its iteration budget (%d iterations) without "
            "meeting the convergence tolerance %g; the returned metrics may not be "
            "converged. Try another method (e.g. 'qd' or 'bs'), raise options.iter_max, "
            "or loosen options.tol." % (method, result.iter,
                                        getattr(options, 'tol', float('nan'))))

    return result


def solver_amva(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Approximate MVA solver.

    Uses approximate MVA methods for larger networks where exact
    MVA would be computationally expensive.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with performance metrics
    """
    from ...pfqn.mva import pfqn_aql, pfqn_bs, pfqn_sqni
    from ...pfqn import pfqn_linearizermx, pfqn_conwayms
    from .amvald import solver_amvald, AmvaldOptions

    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    method = options.method.lower()
    if method.startswith('amva.'):
        method = method[5:]

    M = sn.nstations
    K = sn.nclasses

    # Get chain-level demands
    from ...sn import sn_get_demands_chain, sn_deaggregate_chain_results
    chain_result = sn_get_demands_chain(sn)
    Lchain = chain_result.Lchain
    STchain = chain_result.STchain
    Vchain = chain_result.Vchain
    alpha = chain_result.alpha
    Nchain_float = chain_result.Nchain.flatten()

    # Check if this is a mixed model (has both open and closed chains)
    has_open = np.any(np.isinf(Nchain_float))
    has_closed = np.any(np.isfinite(Nchain_float) & (Nchain_float > 0))
    is_mixed = has_open and has_closed

    # MATLAB (solver_amva.m:220) splits the linearizer family off from the other
    # approximate methods: with a single server and no class-dependent scaling,
    # lin/gflin/egflin go to pfqn_linearizermx, which dispatches on the method
    # name and is the only route to pfqn_gflinearizer / pfqn_egflinearizer. Only
    # cd-scaling or a multiserver model under the default config sends them to
    # solver_amvald. Routing the whole family into solver_amvald collapses gflin
    # and egflin onto lin, because every correction branch there tests for
    # 'lin'/'qdlin' alone. qdlin/fli/qd/qli are NOT part of this family: MATLAB
    # leaves them in the `otherwise` arm, which is solver_amvald.
    lin_family = method in ('lin', 'gflin', 'egflin')
    has_cdscaling = (sn.cdscaling is not None and np.size(sn.cdscaling) > 0) or \
                    (getattr(sn, 'jdscaling', None) is not None and np.size(sn.jdscaling) > 0)
    # mirrors _amva_multiserver_rule (solver_mva.py:2096-2102): dict OR object config
    ms_config = 'default'
    _cfg = getattr(options, 'config', None)
    if isinstance(_cfg, dict):
        ms_config = _cfg.get('multiserver') or 'default'
    elif _cfg is not None:
        ms_config = getattr(_cfg, 'multiserver', None) or 'default'
    ms_config = str(ms_config).lower()
    _finite_ns = np.asarray(sn.nservers, dtype=float).ravel()
    _finite_ns = _finite_ns[np.isfinite(_finite_ns)]
    max_servers = int(np.max(_finite_ns)) if _finite_ns.size > 0 else 1

    # Interlocked flow (Franks 1999, Eq. 4.7). options.interlock arrives CLASS-indexed and is
    # translated to the chain basis the handlers work in. Only solver_amvald carries the
    # correction, so an interlocked model goes there rather than to the closed-form or
    # linearizer branches, which have no interlock term and would drop it silently. Mirrors
    # matlab/src/solvers/MVA/solver_amva.m and the JAR Solver_amva.
    from ...sn import sn_interlock_chain
    interlock_chain = sn_interlock_chain(sn, getattr(options, 'interlock', None))
    has_interlock = interlock_chain is not None and np.size(interlock_chain) > 0

    # BCMP type 1 needs class-independent exponential FCFS service; a station whose
    # per-class means differ is not product form, so it goes to solver_amvald as it
    # does in MATLAB solver_amva.m. ab / schmidt / schmidt-ext ARE the class-dependent
    # FCFS algorithms and must not be diverted.
    from ...sn import sn_has_product_form_not_het_fcfs
    het_fcfs_own = method in ('ab', 'schmidt', 'schmidt-ext')
    pf_not_het = sn_has_product_form_not_het_fcfs(sn) or \
        (het_fcfs_own and sn_has_product_form_not_het_fcfs(sn, check_means=False))

    route_amvald = has_open or has_interlock or method in ('qdlin', 'fli', 'qd', 'qli') \
        or not pf_not_het
    # pf_not_het guards the block: a het-FCFS model is settled by the test above and
    # MATLAB never reaches its lin arm for one (solver_amva.m:397-401 tail instead)
    if lin_family and not has_open and not has_interlock and pf_not_het:
        if has_cdscaling:
            route_amvald = True
        elif max_servers > 1:
            route_amvald = ms_config in ('default', 'softmin', 'seidmann', 'suri')
    # that tail remaps these three before calling solver_amvald, which has no arm for them
    if route_amvald and ms_config in ('conway', 'erlang', 'krzesinski'):
        ms_config = 'default'

    if route_amvald:
        # Use solver_amvald which properly handles mixed models with class switching
        amvald_options = AmvaldOptions(
            method=method if method in ('lin', 'qdlin', 'fli', 'gflin', 'egflin', 'qd', 'qli') else 'egflin',
            iter_tol=options.iter_tol if hasattr(options, 'iter_tol') and options.iter_tol else 1e-4,  # Match MATLAB lineDefaults
            iter_max=options.iter_max if hasattr(options, 'iter_max') and options.iter_max else 100,  # Match MATLAB lineDefaults
            init_sol=getattr(options, 'init_sol', None),
        )
        # ms_config selects the multiserver rule above, so it must reach amvald too
        amvald_options.config.multiserver = ms_config
        # np_priority selects the HOL priority approximation and highvar the high-SCV arm,
        # exactly as in MATLAB solver_amva.m; without this they never leave the caller
        if isinstance(_cfg, dict):
            amvald_options.config.np_priority = str(_cfg.get('np_priority') or 'default').lower()
            amvald_options.config.highvar = str(_cfg.get('highvar') or 'default').lower()
        elif _cfg is not None:
            amvald_options.config.np_priority = str(getattr(_cfg, 'np_priority', None) or 'default').lower()
            amvald_options.config.highvar = str(getattr(_cfg, 'highvar', None) or 'default').lower()
        if has_interlock:
            amvald_options.config.interlock_chain = interlock_chain

        # Get SCVchain and refstatchain
        SCVchain = chain_result.SCVchain if hasattr(chain_result, 'SCVchain') and chain_result.SCVchain is not None else np.ones((M, sn.nchains))
        refstatchain = chain_result.refstatchain if hasattr(chain_result, 'refstatchain') and chain_result.refstatchain is not None else np.zeros((sn.nchains, 1))

        amvald_result = solver_amvald(
            sn, Lchain, STchain, Vchain, alpha,
            Nchain_float, SCVchain, refstatchain, amvald_options
        )

        # Disaggregate to class level. MATLAB's solver_amvald deaggregates
        # internally (solver_amvald.m:239/241); the python solver_amvald returns
        # chain-level measures, so the call is made here instead. The argument
        # convention mirrors those two lines: Q is left to Little's law rather
        # than passed down, and Uchain is forwarded only under lld/cd scaling.
        Xchain = amvald_result.X.reshape(1, -1) if amvald_result.X.ndim == 1 else amvald_result.X
        has_scaling = ((sn.lldscaling is not None and np.size(sn.lldscaling) > 0) or
                       (sn.cdscaling is not None and np.size(sn.cdscaling) > 0) or
                       (getattr(sn, 'jdscaling', None) is not None and np.size(sn.jdscaling) > 0))
        deagg_result = sn_deaggregate_chain_results(
            sn, Lchain, None, STchain, Vchain, alpha,
            None, amvald_result.U if has_scaling else None,
            amvald_result.R, amvald_result.T, None, Xchain
        )

        # Class-dependent station Util post-pass (mirrors MATLAB
        # solver_amvald.m): Util = T*S/peak using the declared sn.cdscalingpeak.
        from .amvald import solver_amvald_cd_peak_post, solver_amvald_jd_peak_post
        solver_amvald_cd_peak_post(sn, deagg_result.U, deagg_result.T)
        solver_amvald_jd_peak_post(sn, deagg_result.U, deagg_result.T)

        result = SolverMVAReturn(
            Q=deagg_result.Q,
            U=deagg_result.U,
            R=deagg_result.R,
            T=deagg_result.T,
            C=deagg_result.C,
            X=deagg_result.X,
            lG=amvald_result.lG,
            runtime=time.time() - start_time,
            method=amvald_options.method,
            it=amvald_result.totiter
        )
        return result

    # For closed-only networks, use simpler methods
    Nchain_float = np.where(np.isfinite(Nchain_float), Nchain_float, 0)
    Nchain = Nchain_float.astype(int)

    # Get think times from delay stations
    sched_dict = sn.sched if sn.sched else {}
    Z = np.zeros(sn.nchains)
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched == SchedStrategy.INF:
            for c in range(sn.nchains):
                Z[c] += Lchain[i, c]

    # Get demands at queueing stations
    L_queue = []
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched not in [SchedStrategy.INF, SchedStrategy.EXT]:
            L_queue.append(Lchain[i, :])

    if len(L_queue) == 0:
        # All delay stations
        L = np.zeros((1, sn.nchains))
    else:
        L = np.array(L_queue)

    # Warm-start queue lengths from a supplied initial solution (options.init_sol),
    # restricted to the queueing-station rows so it aligns with L. The AMVA fixed
    # point then starts from Q0 instead of the default N/M guess.
    Q0 = None
    _init_sol = getattr(options, 'init_sol', None)
    if _init_sol is not None:
        _Qi = np.asarray(_init_sol, dtype=float)
        _qrows = [i for i in range(M) if sched_dict.get(i) not in (SchedStrategy.INF, SchedStrategy.EXT)]
        if _Qi.shape == (M, sn.nchains) and len(_qrows) == L.shape[0]:
            Q0 = _Qi[_qrows, :]
        elif _Qi.shape == L.shape:
            Q0 = _Qi
        if Q0 is not None and (not np.all(np.isfinite(Q0)) or np.any(Q0 < 0)):
            Q0 = None

    # Choose algorithm
    if lin_family:
        # MATLAB solver_amva.m:224-237. pfqn_linearizermx dispatches internally
        # on the method name, which is what separates lin / gflin / egflin; the
        # cd-scaling and default-config multiserver cases already went to
        # solver_amvald above, so only the single-server and the conway /
        # krzesinski multiserver arms reach here. TN comes back unpopulated from
        # both routines (MATLAB discards that output too and rebuilds it from X),
        # so throughput is retiled from XN as the sibling branches do.
        q_rows = [i for i in range(M)
                  if sched_dict.get(i) not in (SchedStrategy.INF, SchedStrategy.EXT)]
        nservers_q = np.array([sn.nservers[i] for i in q_rows], dtype=float)
        sched_q = [sched_dict.get(i, SchedStrategy.FCFS) for i in q_rows]
        _tol = options.tol if getattr(options, 'tol', None) else 1e-4
        _imax = options.iter_max if getattr(options, 'iter_max', None) else 100
        if max_servers > 1 and ms_config == 'conway':
            QN, UN, RN, _CN, XN, _it = pfqn_conwayms(
                L, Nchain, Z, nservers_q, sched_q, _tol, _imax, QN0=Q0)
        else:
            QN, UN, RN, _TN0, _CN, XN, _it = pfqn_linearizermx(
                np.zeros(sn.nchains), L, Nchain, Z, nservers_q, sched_q,
                _tol, _imax, method, Q0)
        XN = np.asarray(XN).flatten()
        TN = np.tile(XN, (QN.shape[0], 1))
        AN = TN.copy()
    elif method == 'bs':
        # Same argument set MATLAB passes (solver_amva.m:148): the solver's
        # tolerance and budget, the warm start, and the per-station scheduling.
        # The api defaults are NOT the solver's, and the scheduling decides the
        # cross-class term at an FCFS station.
        _bs_tol = options.tol if getattr(options, 'tol', None) else 1e-4
        _bs_imax = options.iter_max if getattr(options, 'iter_max', None) else 1000
        _bs_sched = [sched_dict.get(i, SchedStrategy.PS) for i in range(M)
                     if sched_dict.get(i) not in (SchedStrategy.INF, SchedStrategy.EXT)]
        XN, QN, UN, RN, _ = pfqn_bs(L, Nchain, Z, _bs_tol, _bs_imax, Q0, _bs_sched)
        TN = np.tile(XN, (QN.shape[0], 1))
        AN = TN.copy()
    elif method == 'sqni':
        # pfqn_sqni is a closed form for one queueing station with a delay;
        # falling through reported the generic "not implemented" message.
        if L.shape[0] != 1:
            raise ValueError(
                "solver_amva: SQNI is defined for a single queueing station with "
                "a delay. Try with the 'default' or 'lin' methods.")
        Q, U, X = pfqn_sqni(L.flatten(), Nchain, Z)
        XN = X
        QN = Q[:1, :]
        UN = U[:1, :]
        CN = np.zeros((1, sn.nchains))
        RN = np.zeros_like(QN)
        TN = np.tile(XN, (QN.shape[0], 1))
        AN = TN.copy()
    elif method in ('amva', 'default', 'aql'):
        # MATLAB's solver_amva resolves 'default' and 'amva' by population and
        # server count (qd, egflin or lin) rather than falling back to AQL; this
        # arm is reached only for the names that genuinely mean the Schweitzer
        # approximation.
        XN, CN, QN, UN, RN, TN, AN = pfqn_aql(L, Nchain, Z, QN0=Q0)
    else:
        # A method name that reaches here is one no branch above claims. Falling
        # back to AQL returned Schweitzer numbers LABELLED with the requested
        # method, which is indistinguishable from having run it.
        raise ValueError(
            f"solver_amva: the '{method}' method is not implemented by this analyzer")

    # Disaggregate to class level
    Rchain = np.zeros((M, sn.nchains))
    Tchain = np.zeros((M, sn.nchains))
    q_idx = 0
    for i in range(M):
        station_sched = sched_dict.get(i)
        if station_sched not in [SchedStrategy.INF, SchedStrategy.EXT]:
            if q_idx < RN.shape[0]:
                Rchain[i, :] = RN[q_idx, :]
                Tchain[i, :] = TN[q_idx, :]
            q_idx += 1
        elif station_sched == SchedStrategy.INF:
            Tchain[i, :] = XN.flatten()
            Rchain[i, :] = STchain[i, :] * Vchain[i, :]

    Xchain = XN.reshape(1, -1) if XN.ndim == 1 else XN

    deagg_result = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha,
        None, None, Rchain, Tchain, None, Xchain
    )

    result = SolverMVAReturn(
        Q=deagg_result.Q,
        U=deagg_result.U,
        R=deagg_result.R,
        T=deagg_result.T,
        C=deagg_result.C,
        X=deagg_result.X,
        lG=np.nan,
        runtime=time.time() - start_time,
        method=method if method != 'amva' else 'aql',
        it=0
    )

    return result


def _rqna_proc_to_map(proc_entry):
    """Normalize a NetworkStruct proc entry into a MAP pair [D0, D1].

    Native models store processes heterogeneously: a general MAP/PH keeps its
    (D0,D1) matrices (which may be non-renewal, e.g. MMPP2), while Exp/Erlang/
    HyperExp are stored as lightweight dicts. Matrix pairs are returned as-is to
    preserve any autocorrelation; dict forms are expanded to their renewal MAP
    D0 = T, D1 = (-T e) alpha.
    """
    if proc_entry is None:
        return None
    # sn.proc stores (D0, D1); proc_to_map also accepts the legacy descriptors.
    from ...sn.proc_form import proc_to_map
    D0, D1 = proc_to_map(proc_entry)
    if D0 is not None:
        return [np.atleast_2d(D0), np.atleast_2d(D1)]
    raise RuntimeError("RQNA: unsupported process representation %r" % type(proc_entry))


def _rqna_svc_idc(svc_maps, t):
    """Service IDC vector I_{s,a}(t). t is scalar (same for all queues) or a
    length-nq array giving the per-queue evaluation time."""
    nq = len(svc_maps)
    ta = np.atleast_1d(np.asarray(t, dtype=np.float64)).ravel()
    if ta.size == 1:
        tt = np.repeat(ta[0], nq)
    else:
        tt = ta
    out = np.zeros(nq)
    for a in range(nq):
        D0, D1 = svc_maps[a][0], svc_maps[a][1]
        out[a] = float(map_count_idc(D0, D1, np.array([tt[a]]))[0])
    return out


def _rqna_geom_map(mapproc, p):
    """Geometric random sum of i.i.d. PH service times with success prob (1-p):
    PH(alpha,T) -> PH(alpha, T + p*t0*alpha), t0 = -T*e. Yields the
    head-of-line immediate-feedback-eliminated service (Whitt-You Section 4.1)."""
    D0 = np.asarray(mapproc[0], dtype=np.float64)
    n = D0.shape[0]
    e = np.ones((n, 1))
    t0 = -D0 @ e
    al = map_pie(D0, np.asarray(mapproc[1], dtype=np.float64)).reshape(1, n)
    D0m = D0 + p * (t0 @ al)
    D1m = (1.0 - p) * (t0 @ al)
    return [D0m, D1m]


def _rqna_phat(P, rho, a):
    """Near-immediate feedback probability at station a: probability that a
    customer departing a returns to a before visiting any station with strictly
    higher traffic intensity (Whitt-You flows paper eq. 3.8/3.9, H={a}).

    Delegated to npfqn_feedback_elim so that the solver and the API function
    cannot drift apart: they answer the same question, and a private copy of the
    rule here is how the two came to differ on ties in the first place."""
    from ...npfqn.feedback_elim import npfqn_feedback_elim
    return float(npfqn_feedback_elim(P, rho)['feedbackProb'][a])


def _rqna_fp_ext_idc_R(t, R, Hc, G, lambda0, c2fun_i, lam0R):
    """First-passage external arrival IDC into each retained station of the
    reduced network (superposition of direct external and cloud-entering splits)."""
    m = len(R)
    Ir = np.ones(m)
    for rr in range(m):
        if lam0R[rr] <= 0:
            continue
        num = lambda0[R[rr]] * c2fun_i(R[rr], t)
        for ii in range(len(Hc)):
            g = G[ii, rr]
            num += lambda0[Hc[ii]] * g * (g * c2fun_i(Hc[ii], t) + (1.0 - g))
        Ir[rr] = num / lam0R[rr]
    return Ir


def _rqna_elim_response(a, P, rho, lam, mu, cs2, svc_maps,
                        lambda0, arv_map, qsplit, corrections):
    """Near-immediate feedback elimination at station a (Whitt-You Algorithm 2 /
    Corollary 4.2). Reduced network: collapse only the equal/lower-rho cloud Hc
    to instantaneous switches (censoring), retaining strictly higher-rho stations
    as real queues so global traffic rates and upstream-bottleneck arrival
    variability are preserved. The near-immediate self-return at a is removed by
    immediate-feedback elimination and a's service replaced by the geometric-sum
    service. The IDC equations are re-solved on the reduced network and RQ applied
    at a with the per-visit adjustment W = (1-phat) Wtilde."""
    arv_D0 = np.asarray(arv_map[0], dtype=np.float64)
    arv_D1 = np.asarray(arv_map[1], dtype=np.float64)
    arv_idc_inf = map_idc(arv_D0, arv_D1)

    nq = P.shape[0]
    Hc = [i for i in range(nq) if rho[i] <= rho[a] + 1e-9 and i != a]
    Hi = [i for i in range(nq) if rho[i] > rho[a] + 1e-9 and i != a]
    R = [a] + Hi                              # retained stations, a first
    m = len(R)
    Ridx = np.array(R, dtype=int)
    Hc = np.array(Hc, dtype=int)

    if len(Hc) == 0:
        Fhc = np.zeros((0, 0))
        G = np.zeros((0, m))
        Pred = P[np.ix_(Ridx, Ridx)].copy()
    else:
        Fhc = np.linalg.inv(np.eye(len(Hc)) - P[np.ix_(Hc, Hc)])
        G = Fhc @ P[np.ix_(Hc, Ridx)]         # first-passage Hc-station -> R
        Pred = P[np.ix_(Ridx, Ridx)] + P[np.ix_(Ridx, Hc)] @ Fhc @ P[np.ix_(Hc, Ridx)]

    phat = Pred[0, 0]                          # near-immediate return prob at a
    phat = min(max(phat, 0.0), 1.0 - 1e-9)

    # immediate-feedback elimination at a (position 0 in R)
    if phat > 0:
        Pred[0, :] = Pred[0, :] / (1.0 - phat)
    Pred[0, 0] = 0.0

    def c2fun_i(i, t):
        return qsplit[i] * float(map_count_idc(arv_D0, arv_D1, np.array([t]))[0]) + (1.0 - qsplit[i])

    # first-passage external arrival rate and IDC into each retained station
    lam0R = np.zeros(m)
    for rr in range(m):
        lam0R[rr] = lambda0[R[rr]]
        for ii in range(len(Hc)):
            lam0R[rr] += lambda0[Hc[ii]] * G[ii, rr]
    c2a0R = np.zeros(m)
    for rr in range(m):
        if lam0R[rr] <= 0:
            continue
        s = lambda0[R[rr]] * (qsplit[R[rr]] * arv_idc_inf + (1.0 - qsplit[R[rr]]))
        for ii in range(len(Hc)):
            g = G[ii, rr]
            ci = qsplit[Hc[ii]] * arv_idc_inf + (1.0 - qsplit[Hc[ii]])
            s += lambda0[Hc[ii]] * g * (g * ci + (1.0 - g))
        c2a0R[rr] = s / lam0R[rr]

    def a0IdcR(t):
        return _rqna_fp_ext_idc_R(t, R, Hc, G, lambda0, c2fun_i, lam0R)

    # service data on R; station a gets the geometric-sum (folded) service
    muR = mu[Ridx].copy()
    cs2R = cs2[Ridx].copy()
    svcR = [svc_maps[R[rr]] for rr in range(m)]
    svcR[0] = _rqna_geom_map(svc_maps[a], phat)
    muR[0] = (1.0 - phat) * mu[a]
    cs2R[0] = phat + (1.0 - phat) * cs2[a]

    def sIdcR(t):
        return _rqna_svc_idc(svcR, t)

    ctxR = npfqn_traffic_idc(lam0R, Pred, c2a0R, a0IdcR, muR, cs2R, sIdcR, corrections)

    def IaFunA(x):
        return ctxR.IaFun(x)[0]

    _, Wt, _, _ = qsys_gig1_rq(rho[a], muR[0], cs2R[0], IaFunA)
    return (1.0 - phat) * Wt + 1.0 / mu[a]     # per-visit adjustment (mean visits 1/(1-phat))


def solver_rqna(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Robust Queueing Network Analyzer (RQNA) based on indices of dispersion.

    Approximates the steady-state performance of a single-class open queueing
    network of single-server FCFS queues with Markovian routing and general
    (non-renewal) external arrival and (non-exponential) service processes.

    Reference: W. Whitt and W. You (2018), "A Robust Queueing Network Analyzer
    Based on Indices of Dispersion", INFORMS J. on Computing. Implements
    Algorithm 1 (traffic-rate equations, limiting variability equations,
    time-dependent IDC equations with default alpha/beta corrections, and the
    robust-queueing workload approximation), plus the near-immediate feedback
    elimination of Algorithm 2 / Section 4.2 (default on).

    References:
        MATLAB: matlab/src/solvers/MVA/solver_rqna.m
    """
    start_time = time.time()
    if options is None:
        options = SolverMVAOptions()

    # One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
    # the same question before the report offers 'rqna', so the sentence a
    # caller reads here is the sentence that kept the row off the report.
    _rqna_ok, _rqna_reason = mva_supports_single_class_open(sn, 'rqna')
    if not _rqna_ok:
        raise RuntimeError(_rqna_reason)
    if np.any(np.isfinite(sn.njobs)):
        raise RuntimeError("RQNA supports open networks only (no closed classes).")

    M = sn.nstations
    K = sn.nclasses  # == 1
    tol = options.tol if hasattr(options, 'tol') else 1e-4

    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    X = np.zeros((1, K))

    # ----- identify source and queueing stations -----
    isSource = np.zeros(M, dtype=bool)
    schedInf = np.zeros(M, dtype=bool)
    sched_dict = sn.sched if sn.sched else {}
    for i in range(M):
        nd = int(sn.stationToNode[i])
        isSource[i] = (sn.nodetype[nd] == NodeType.SOURCE)
        schedInf[i] = (sched_dict.get(i, SchedStrategy.FCFS) == SchedStrategy.INF)
    srcList = np.where(isSource)[0]
    qstat = np.where(~isSource)[0]            # queueing (and delay) stations
    nq = len(qstat)

    # single-class station-to-station routing matrix (indexed by station index,
    # matching the MATLAB reference which assumes station idx == stateful idx)
    rtS = np.asarray(sn.rt, dtype=np.float64)

    if len(srcList) == 0:
        raise RuntimeError("RQNA requires an open network with a Source station.")
    src = int(srcList[0])
    arv_map = _rqna_proc_to_map(sn.proc[src][0])
    arv_D0 = np.asarray(arv_map[0], dtype=np.float64)
    arv_D1 = np.asarray(arv_map[1], dtype=np.float64)
    lambda_src = map_lambda(arv_D0, arv_D1)
    c2_src = map_idc(arv_D0, arv_D1)

    def arvIdc(t):
        return float(map_count_idc(arv_D0, arv_D1, np.array([t]))[0])

    # ----- per-queue data -----
    mu = np.zeros(nq)
    cs2 = np.zeros(nq)
    lambda0 = np.zeros(nq)
    svc_maps = [None] * nq
    qsplit = np.zeros(nq)
    P = np.zeros((nq, nq))
    for a in range(nq):
        ia = int(qstat[a])
        mu[a] = sn.rates[ia, 0]
        cs2[a] = sn.scv[ia, 0]
        svc_maps[a] = _rqna_proc_to_map(sn.proc[ia][0])
        qsplit[a] = rtS[src, ia]
        lambda0[a] = lambda_src * qsplit[a]
        for b in range(nq):
            ib = int(qstat[b])
            P[a, b] = rtS[ia, ib]

    # external arrival IDC seen by each queue = split of the source process
    c2a0 = qsplit * c2_src + (1.0 - qsplit)

    def a0IdcFun(t):
        return qsplit * arvIdc(t) + (1.0 - qsplit)

    def sIdcFun(t):
        return _rqna_svc_idc(svc_maps, t)

    # ----- traffic variability equations -----
    corrections = {}
    cfg = getattr(options, 'config', None)
    if isinstance(cfg, dict):
        if 'rqna_alpha' in cfg:
            corrections['alpha'] = cfg['rqna_alpha']
        if 'rqna_beta' in cfg:
            corrections['beta'] = cfg['rqna_beta']
    elif cfg is not None:
        if hasattr(cfg, 'rqna_alpha'):
            corrections['alpha'] = cfg.rqna_alpha
        if hasattr(cfg, 'rqna_beta'):
            corrections['beta'] = cfg.rqna_beta

    ctx = npfqn_traffic_idc(lambda0, P, c2a0, a0IdcFun, mu, cs2, sIdcFun, corrections)
    lam = ctx.lam
    rho = ctx.rho

    # near-immediate feedback elimination is applied by default
    doElim = True
    if isinstance(cfg, dict) and 'rqna_feedback_elim' in cfg:
        doElim = bool(cfg['rqna_feedback_elim'])
    elif cfg is not None and hasattr(cfg, 'rqna_feedback_elim'):
        doElim = bool(cfg.rqna_feedback_elim)

    # ----- per-queue robust-queueing workload and performance -----
    for a in range(nq):
        ia = int(qstat[a])
        T[ia, 0] = lam[a]
        if lam[a] <= 0:
            continue
        if schedInf[ia]:
            # infinite-server (delay) station: no waiting
            U[ia, 0] = lam[a] / mu[a]
            Q[ia, 0] = lam[a] / mu[a]
            R[ia, 0] = 1.0 / mu[a]
            continue
        phat = 0.0
        if doElim:
            phat = _rqna_phat(P, rho, a)
        if phat > tol:
            R[ia, 0] = _rqna_elim_response(a, P, rho, lam, mu, cs2, svc_maps,
                                           lambda0, arv_map, qsplit, corrections)
        else:
            def IaFun_a(x, _a=a):
                return ctx.IaFun(x)[_a]
            _, Wa, _, _ = qsys_gig1_rq(rho[a], mu[a], cs2[a], IaFun_a)
            R[ia, 0] = Wa + 1.0 / mu[a]        # per-visit response = waiting + service
        U[ia, 0] = rho[a]
        Q[ia, 0] = lam[a] * R[ia, 0]           # mean number in system (Little, incl. service)

    T[src, 0] = lambda_src
    C_result = np.sum(R, axis=0, keepdims=True)
    X[0, 0] = lambda_src
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    C_result[np.isnan(C_result)] = 0

    return SolverMVAReturn(
        Q=Q, U=U, R=R, T=T, C=C_result, X=X,
        lG=np.nan, runtime=time.time() - start_time, method='rqna', it=1
    )


def solver_rqt(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Robust Queueing Network Analyzer (RQNA) of Robust Queueing Theory.

    Estimates the steady-state performance of a single-class open network of
    FCFS queues with Markovian routing by replacing the stochastic primitives
    with polyhedral uncertainty sets and taking a worst-case view of each node in
    isolation. The algorithm is Section 7.2 of the reference: the external
    streams get Gamma_a = sigma_a, the effective arrival process at each node
    follows from the network characterization of Theorem 10, the service
    variability parameter from the adaptation of Section 7.1, and the system time
    at each node is the worst-case bound of Theorem 3. The published step 3, path
    enumeration, is not needed here: LINE aggregates per-node system times into
    per-class response times through the visit ratios.

    The adaptation is regressed against simulation in heavy traffic, so accuracy
    degrades at low utilization: on M/M/1 the error is about 5% at rho=0.9 but
    over 50% at rho=0.5.

    Reference: C. Bandi, D. Bertsimas, N. Youssef (2015), "Robust Queueing
    Theory", Operations Research 63(3), 676-700.

    References:
        MATLAB: matlab/src/solvers/MVA/solver_rqt.m
    """
    start_time = time.time()
    if options is None:
        options = SolverMVAOptions()

    # One predicate for the gate and the run: SolverMVA.supportsModelMethod asks
    # the same question before the report offers 'rqt', so the sentence a caller
    # reads here is the sentence that kept the row off the report.
    _rqt_ok, _rqt_reason = mva_supports_single_class_open(sn, 'rqt')
    if not _rqt_ok:
        raise RuntimeError(_rqt_reason)
    if np.any(np.isfinite(sn.njobs)):
        raise RuntimeError("RQT supports open networks only (no closed classes).")

    M = sn.nstations
    K = sn.nclasses  # == 1

    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    X = np.zeros((1, K))

    # ----- configuration -----
    regime = 'independent'
    use_exact = False
    alpha_a_cfg = None
    alpha_s_cfg = None
    cfg = getattr(options, 'config', None)

    def _cfg(name):
        if isinstance(cfg, dict):
            return cfg.get(name, None)
        if cfg is not None and hasattr(cfg, name):
            return getattr(cfg, name)
        return None

    if _cfg('rqt_regime') is not None:
        regime = _cfg('rqt_regime')
    if _cfg('rqt_exact') is not None:
        use_exact = bool(_cfg('rqt_exact'))
    if _cfg('rqt_alpha_a') is not None:
        alpha_a_cfg = np.atleast_1d(np.asarray(_cfg('rqt_alpha_a'), dtype=float))
    if _cfg('rqt_alpha_s') is not None:
        alpha_s_cfg = np.atleast_1d(np.asarray(_cfg('rqt_alpha_s'), dtype=float))

    # ----- identify source and queueing stations -----
    isSource = np.zeros(M, dtype=bool)
    schedInf = np.zeros(M, dtype=bool)
    sched_dict = sn.sched if sn.sched else {}
    for i in range(M):
        nd = int(sn.stationToNode[i])
        isSource[i] = (sn.nodetype[nd] == NodeType.SOURCE)
        schedInf[i] = (sched_dict.get(i, SchedStrategy.FCFS) == SchedStrategy.INF)
    srcList = np.where(isSource)[0]
    if len(srcList) == 0:
        raise RuntimeError("RQT requires an open network with a Source station.")
    src = int(srcList[0])
    qstat = np.where(~isSource)[0]
    nq = len(qstat)

    rtS = np.asarray(sn.rt, dtype=np.float64)

    # ----- external arrival process -----
    arv_map = _rqna_proc_to_map(sn.proc[src][0])
    arv_D0 = np.asarray(arv_map[0], dtype=np.float64)
    arv_D1 = np.asarray(arv_map[1], dtype=np.float64)
    lambda_src = map_lambda(arv_D0, arv_D1)
    sigma_a_src = np.sqrt(map_scv(arv_D0, arv_D1)) / lambda_src
    alpha_a_src = 2.0 if alpha_a_cfg is None else float(alpha_a_cfg[0])

    # ----- per-node primitives -----
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    mu = np.zeros(nq)
    sigma_s = np.zeros(nq)
    nserv = np.ones(nq, dtype=int)
    alpha_s = np.full(nq, 2.0)
    lambda0 = np.zeros(nq)
    Gamma0 = np.zeros(nq)
    alpha0 = np.full(nq, alpha_a_src)
    F = np.zeros((nq, nq))
    for a in range(nq):
        ia = int(qstat[a])
        mu[a] = sn.rates[ia, 0]
        sigma_s[a] = np.sqrt(sn.scv[ia, 0]) / mu[a]
        if np.isfinite(nservers[ia]) and nservers[ia] > 0:
            nserv[a] = int(nservers[ia])
        if alpha_s_cfg is not None:
            alpha_s[a] = float(alpha_s_cfg[min(a, len(alpha_s_cfg) - 1)])
        # the source stream reaches node a thinned by q, Theorem 6
        q = rtS[src, ia]
        lambda0[a] = lambda_src * q
        if q > 0:
            Gamma0[a] = sigma_a_src * (1.0 / q) ** (1.0 / alpha_a_src)
        for b in range(nq):
            F[a, b] = rtS[ia, int(qstat[b])]

    # ----- effective arrival processes, Theorem 10 -----
    lam, Gamma_a, alpha_a = npfqn_traffic_rqt(lambda0, Gamma0, alpha0, F)

    # ----- per-node worst-case analysis -----
    for a in range(nq):
        ia = int(qstat[a])
        T[ia, 0] = lam[a]
        if lam[a] <= 0:
            continue
        if schedInf[ia]:
            U[ia, 0] = lam[a] / mu[a]
            Q[ia, 0] = lam[a] / mu[a]
            R[ia, 0] = 1.0 / mu[a]
            continue
        rho = lam[a] / (nserv[a] * mu[a])
        Gamma_s = qsys_gigk_rqt_gamma(rho, mu[a], Gamma_a[a], sigma_s[a],
                                      int(nserv[a]), alpha_a[a], regime)
        Wa, _, Sworst = qsys_gigk_rqt(lam[a], mu[a], Gamma_a[a], Gamma_s,
                                      int(nserv[a]), alpha_a[a], alpha_s[a])
        if use_exact:
            Wa = Sworst
        R[ia, 0] = Wa
        U[ia, 0] = rho
        Q[ia, 0] = lam[a] * Wa                 # Little's law, number in system

    T[src, 0] = lambda_src
    C_result = np.sum(R, axis=0, keepdims=True)
    X[0, 0] = lambda_src
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    C_result[np.isnan(C_result)] = 0

    return SolverMVAReturn(
        Q=Q, U=U, R=R, T=T, C=C_result, X=X,
        lG=np.nan, runtime=time.time() - start_time, method='rqt', it=1
    )


def solver_qna(
    sn: NetworkStruct,
    options: Optional[SolverMVAOptions] = None
) -> SolverMVAReturn:
    """
    Queueing Network Analyzer.

    Provides approximate analysis for general queueing networks
    including those with non-product-form characteristics.

    Implementation based on N. Gautaum's "Analysis of Queues" (CRC Press, 2012),
    Section 7.2.3 with minor corrections.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMVAReturn with performance metrics

    References:
        MATLAB: matlab/src/solvers/MVA/solver_qna.m
    """
    start_time = time.time()

    if options is None:
        options = SolverMVAOptions()

    # One predicate for the gate and the run: SolverMVA.getMethodFeatureSet
    # withholds 'qna' for a discipline this loop has no arm for, and the loop
    # refuses it by name rather than leaving that station's row of Q, U, R and T
    # at zero and returning the table as a solution.
    _qna_ok, _qna_reason = mva_supports_qna_scheduling(sn)
    if not _qna_ok:
        raise RuntimeError(_qna_reason)

    K = sn.nclasses
    M = sn.nstations
    C = sn.nchains

    # Extract parameters from network. sn.rt and sn.visits are indexed by
    # stateful node: project them onto stations.
    rt, V = sn_rt_stations(sn)
    S = 1.0 / (sn.rates + 1e-10)  # Service times (M, K)
    scv = sn.scv.copy() if sn.scv is not None else np.ones((M, K))
    scv[np.isnan(scv)] = 0

    # Initialize results
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    X = np.zeros((1, K))

    # Chain information
    lambda_chain = np.zeros(C)
    inchain = sn.inchain if sn.inchain else {}

    tol = options.iter_tol if hasattr(options, 'iter_tol') else 1e-6
    max_iter = options.iter_max if hasattr(options, 'iter_max') else 1000

    # Compute departure process at source
    a1 = np.zeros((M, K))  # Arrival rates
    a2 = np.zeros((M, K))  # SCVs of arrivals
    d2 = np.zeros(M)  # SCVs of departure processes
    f2 = np.ones((M * K, M * K))  # SCV of each flow pair

    # deterministic (round-robin) split degrees, k=1 where the split is Markovian
    kRR = npfqn_traffic_split_rr(sn)
    for ist in range(M):
        for r in range(K):
            if kRR[ist, r] <= 1:
                continue
            for jst in range(M):
                for s in range(K):
                    idx_from = ist * K + r
                    idx_to = jst * K + s
                    if idx_from < rt.shape[0] and idx_to < rt.shape[1] and rt[idx_from, idx_to] > 0:
                        f2[idx_from, idx_to] = 1 + rt[idx_from, idx_to] * (1 - kRR[ist, r])

    # Initialize throughputs at source
    for c in range(C):
        if c in inchain:
            classes_c = inchain[c]
            refstat = int(sn.refstat.flatten()[classes_c[0]])
            if refstat < M:
                for k in classes_c:
                    k = int(k)
                    if k < K and sn.njobs[k] == np.inf:  # Open class
                        T[refstat, k] = sn.rates[refstat, k]
                        lambda_chain[c] += sn.rates[refstat, k]

    # Compute initial SCVs at source
    for c in range(C):
        if c in inchain:
            classes_c = inchain[c]
            refstat = int(sn.refstat.flatten()[classes_c[0]])
            if refstat < M and lambda_chain[c] > 0:
                d2_c = 0
                for k in classes_c:
                    k = int(k)
                    if k < K and sn.njobs[k] == np.inf:
                        d2_c += scv[refstat, k] * sn.rates[refstat, k]
                d2[refstat] = d2_c / lambda_chain[c]

    # Main QNA iteration: flow fixed point on the queue lengths, driven by
    # the generic DA successive-substitution driver
    from ...da import da_fpi

    def qna_sweep(x, itnum):
        xref = Q.copy()

        # Normalize queue lengths for closed classes
        for c in range(C):
            if c in inchain:
                classes_c = inchain[c]
                njobs_c = np.sum([sn.njobs[int(k)] for k in classes_c if int(k) < K])
                if np.isfinite(njobs_c) and njobs_c > 0:
                    Q_sum = np.sum(Q[:, [int(k) for k in classes_c if int(k) < K]])
                    if Q_sum > 0:
                        Q[:, [int(k) for k in classes_c if int(k) < K]] *= njobs_c / Q_sum

        # Update throughputs
        if itnum == 1:
            for c in range(C):
                if c in inchain:
                    classes_c = inchain[c]
                    if lambda_chain[c] > 0:
                        for m in range(M):
                            for k in classes_c:
                                k = int(k)
                                if k < K:
                                    T[m, k] = V[m, k] * lambda_chain[c]

        # Superposition: compute arrival process parameters at each station
        for ist in range(M):
            a1[ist, :] = 0
            a2[ist, :] = 0
            lambda_i = np.sum(T[ist, :])

            for jst in range(M):
                for r in range(K):
                    for s in range(K):
                        idx_from = (jst) * K + r
                        idx_to = (ist) * K + s
                        if idx_from < rt.shape[0] and idx_to < rt.shape[1] and rt[idx_from, idx_to] > 0:
                            a1[ist, s] += T[jst, r] * rt[idx_from, idx_to]
                            if lambda_i > 0:
                                a2[ist, s] += (1.0 / lambda_i) * f2[idx_from, idx_to] * T[jst, r] * rt[idx_from, idx_to]

        # Handle different scheduling strategies
        sched_dict = sn.sched if sn.sched else {}
        for ist in range(M):
            station_sched = sched_dict.get(ist, SchedStrategy.FCFS)

            if station_sched == SchedStrategy.INF:  # Delay station
                for k in range(K):
                    T[ist, k] = a1[ist, k]
                    Q[ist, k] = T[ist, k] * S[ist, k] * V[ist, k]
                    U[ist, k] = Q[ist, k]
                    R[ist, k] = S[ist, k] * V[ist, k]
                d2[ist] = np.sum(a2[ist, :]) / np.sum(a1[ist, :]) if np.sum(a1[ist, :]) > 0 else 1

            elif station_sched == SchedStrategy.PS:  # Processor Sharing
                for c in range(C):
                    if c in inchain:
                        classes_c = inchain[c]
                        for k in classes_c:
                            k = int(k)
                            if k < K and lambda_chain[c] > 0:
                                T[ist, k] = lambda_chain[c] * V[ist, k]
                                U[ist, k] = S[ist, k] * T[ist, k]

                        # Compute queue length using approximation
                        Nc = np.sum([sn.njobs[int(k)] for k in classes_c if int(k) < K and np.isfinite(sn.njobs[int(k)])])
                        Uden = min(1.0 - options.tol if hasattr(options, 'tol') else 1e-6, np.sum(U[ist, :]))

                        for k in classes_c:
                            k = int(k)
                            if k < K:
                                if Uden < 1.0:
                                    Q[ist, k] = (U[ist, k] - U[ist, k] ** (Nc + 1)) / (1.0 - Uden)
                                else:
                                    Q[ist, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                                R[ist, k] = Q[ist, k] / (T[ist, k] + 1e-10)

            elif station_sched == SchedStrategy.FCFS:  # FCFS queue
                mu_ist = sn.rates[ist, :]
                rho_ist_class = a1[ist, :] / (1e-10 + mu_ist)
                lambda_ist = np.sum(a1[ist, :])
                mi = sn.nservers[ist] if ist < len(sn.nservers) else 1
                rho_ist = np.sum(rho_ist_class) / mi

                if rho_ist < 1.0 - (options.tol if hasattr(options, 'tol') else 1e-6):
                    # Compute waiting time using diffusion approximation
                    mubar = lambda_ist / (rho_ist + 1e-10)
                    alpha_mi = (rho_ist ** mi + rho_ist) / 2.0 if rho_ist > 0.7 else rho_ist ** ((mi + 1) / 2.0)
                    c2 = -1

                    for r in range(K):
                        if mu_ist[r] > 0:
                            c2 += (a1[ist, r] / (lambda_ist + 1e-10)) * (mubar / (mi * mu_ist[r])) ** 2 * (scv[ist, r] + 1)

                    Wiq = (alpha_mi / mubar) * (1.0 / (1.0 - rho_ist + 1e-10)) * (np.sum(a2[ist, :]) + c2) / 2.0

                    for k in range(K):
                        Q[ist, k] = a1[ist, k] / (mu_ist[k] + 1e-10) + a1[ist, k] * Wiq
                        T[ist, k] = a1[ist, k]
                        U[ist, k] = T[ist, k] * S[ist, k] / mi
                        R[ist, k] = Q[ist, k] / (T[ist, k] + 1e-10)

                    d2[ist] = 1 + rho_ist ** 2 * (c2 - 1) / np.sqrt(mi) + (1 - rho_ist ** 2) * (np.sum(a2[ist, :]) - 1)
                else:
                    # System overloaded
                    for k in range(K):
                        Q[ist, k] = sn.njobs[k] if np.isfinite(sn.njobs[k]) else 1
                    d2[ist] = 1

        # Update flow SCVs for splitting
        for ist in range(M):
            for jst in range(M):
                for r in range(K):
                    for s in range(K):
                        idx_from = ist * K + r
                        idx_to = jst * K + s
                        if idx_from < rt.shape[0] and idx_to < rt.shape[1] and rt[idx_from, idx_to] > 0:
                            # k-fold convolution then Bernoulli thinning at q=k*p: C^2 = (q/k)*d2+1-q
                            f2[idx_from, idx_to] = 1 + rt[idx_from, idx_to] * (d2[ist] - kRR[ist, r])

        return Q.copy(), xref

    _, iteration, _ = da_fpi(qna_sweep, Q.copy(), max_iter, tol, nanstop=True)

    # Final cleanup and normalization
    Q = np.abs(Q)
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    T[np.isnan(T)] = 0

    C_result = np.sum(R, axis=0, keepdims=True)

    result = SolverMVAReturn(
        Q=Q,
        U=U,
        R=R,
        T=T,
        C=C_result,
        X=X,
        lG=0.0,  # QNA does not compute normalizing constant
        runtime=time.time() - start_time,
        method='qna',
        it=iteration
    )

    return result


__all__ = [
    'MVAResult',
    'solver_mva_analyzer',
    'solver_amva',
    'solver_qna',
]
