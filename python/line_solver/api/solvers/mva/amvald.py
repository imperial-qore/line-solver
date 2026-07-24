"""
Approximate MVA for Load-Dependent (AMVA-LD) networks.

This module implements the approximate MVA solver for load-dependent
queueing networks using QD-AMVA style corrections with under-relaxation.

References:
    Original MATLAB: matlab/src/solvers/MVA/solver_amvald.m
    Original MATLAB: matlab/src/solvers/MVA/solver_amvald_forward.m
"""

import functools
import numpy as np
from dataclasses import dataclass
from typing import Optional, Tuple, Dict
from enum import IntEnum

# Import the proper pfqn_lldfun and pfqn_cdfun from utils
from ...pfqn.utils import pfqn_lldfun, pfqn_cdfun, pfqn_jdfun
# Import SchedStrategy from the canonical source
from ....lang.base import SchedStrategy

SCHED_INF, SCHED_FCFS, SCHED_SIRO, SCHED_PS, SCHED_HOL, SCHED_LCFSPR, SCHED_EXT = 0, 1, 2, 3, 4, 5, 6


def _divergence_is_expected(fn):
    """
    Run `fn` with IEEE-754 overflow/invalid reporting muted.

    On an unstable open station (offered load rho >= 1) the AMVA fixed point
    does not exist and the iteration below diverges by construction: queue
    lengths grow until they overflow to +Inf, after which Inf - Inf yields NaN
    in the convergence test. This is the intended sample path -- MATLAB's
    solver_amvald.m iterates identically and the saturated station is reported
    downstream as Util = 1, QLen = RespT = Inf by the caller's unstable-open
    capping. Only the numpy floating-point *flags* are muted here; every value
    produced is bit-identical to an unmuted run, so no divergence is hidden
    from the results.
    """
    @functools.wraps(fn)
    def wrapper(*args, **kwargs):
        with np.errstate(over='ignore', invalid='ignore'):
            return fn(*args, **kwargs)
    return wrapper


def _normalize_sched_strategy(sched_value):
    """
    Normalize a SchedStrategy value to the local SchedStrategy enum.

    This handles the case where sched_value comes from a different SchedStrategy
    enum class instance (e.g., from a different module import path).
    """
    # Get the integer value regardless of which SchedStrategy class it's from
    if hasattr(sched_value, 'value'):
        val = sched_value.value
    elif isinstance(sched_value, int):
        val = sched_value
    else:
        val = SchedStrategy.FCFS.value

    # Map to the local SchedStrategy enum
    for member in SchedStrategy:
        if member.value == val:
            return member
    return SchedStrategy.FCFS


def _sched_to_int(sched_strategy) -> int:
    """Convert SchedStrategy to integer for JIT functions."""
    # Normalize first to handle different enum class instances
    sched_strategy = _normalize_sched_strategy(sched_strategy)

    if sched_strategy == SchedStrategy.INF:
        return SCHED_INF
    elif sched_strategy == SchedStrategy.FCFS:
        return SCHED_FCFS
    elif sched_strategy == SchedStrategy.SIRO:
        return SCHED_SIRO
    elif sched_strategy == SchedStrategy.PS:
        return SCHED_PS
    elif sched_strategy == SchedStrategy.HOL:
        return SCHED_HOL
    elif sched_strategy == SchedStrategy.LCFSPR:
        return SCHED_LCFSPR
    elif sched_strategy == SchedStrategy.EXT:
        return SCHED_EXT
    else:
        return SCHED_FCFS  # Default


def _build_sched_array(sched: Dict[int, int], M: int) -> np.ndarray:
    """Convert sched dict to numpy array for JIT functions."""
    sched_arr = np.zeros(M, dtype=np.int64)
    for k in range(M):
        sched_k = sched.get(k, SchedStrategy.FCFS)
        sched_arr[k] = _sched_to_int(sched_k)
    return sched_arr


@dataclass
class AmvaldOptions:
    """Options for AMVA-LD solver."""
    method: str = 'default'
    iter_tol: float = 1e-4  # Convergence tolerance (matches MATLAB lineDefaults iter_tol=1e-4)
    iter_max: int = 100     # Match MATLAB lineDefaults iter_max=100
    tol: float = 1e-4       # Tolerance for all other uses (matches MATLAB default tol=1e-4)
    verbose: bool = False
    init_sol: Optional[np.ndarray] = None  # Warm-start Qchain (M x K); ignored if shape is stale

    @dataclass
    class Config:
        multiserver: str = 'default'
        np_priority: str = 'default'
        highvar: str = 'default'

    config: Config = None

    def __post_init__(self):
        if self.config is None:
            self.config = AmvaldOptions.Config()


def solver_amvald_cd_peak_post(sn, U, T):
    """Post-pass on the deaggregated per-class results of AMVA-LD.

    Stations with limited class dependence report utilization as T*S/peak,
    where peak is the user-declared per-class peak rate scaling
    (sn.cdscalingpeak, from set_class_dependence(beta, peak_rate_per_class)),
    matching the T*S/c convention of ordinary multiserver stations.

    Mirrors the post-pass of MATLAB solver_amvald.m (applied after
    sn_deaggregate_chain_results). U is modified in place and returned.

    Args:
        sn: NetworkStruct
        U: Per-class utilizations (nstations x nclasses), modified in place
        T: Per-class throughputs (nstations x nclasses)
    """
    cd = getattr(sn, 'cdscaling', None)
    if cd is None or len(cd) == 0:
        return U
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    if np.any(np.isinf(njobs)):
        return U
    peak = sn.cdscalingpeak
    M = sn.nstations
    Kcls = sn.nclasses  # U and rates are class-indexed
    for ist in range(M):
        beta = (cd.get(ist) if isinstance(cd, dict)
                else (cd[ist] if ist < len(cd) else None))
        if beta is None or not callable(beta):
            continue
        for k in range(Kcls):
            bmax = peak[ist, k]
            rate = sn.rates[ist, k]
            if np.isfinite(rate) and rate > 0 and bmax > 0:
                U[ist, k] = T[ist, k] / rate / bmax
            else:
                U[ist, k] = 0.0
    return U


def solver_amvald_jd_peak_post(sn, U, T):
    """Joint-dependence twin of solver_amvald_cd_peak_post.

    Stations with limited joint dependence report utilization as T*S/peak,
    using the declared per-class peak rate scaling sn.jdscalingpeak. U is
    modified in place and returned.
    """
    jd = getattr(sn, 'jdscaling', None)
    if jd is None or len(jd) == 0:
        return U
    njobs = np.asarray(sn.njobs, dtype=float).flatten()
    if np.any(np.isinf(njobs)):
        return U
    peak = sn.jdscalingpeak
    M = sn.nstations
    Kcls = sn.nclasses
    for ist in range(M):
        eta = (jd.get(ist) if isinstance(jd, dict)
               else (jd[ist] if ist < len(jd) else None))
        if eta is None or not callable(eta):
            continue
        for k in range(Kcls):
            bmax = peak[ist, k]
            rate = sn.rates[ist, k]
            if np.isfinite(rate) and rate > 0 and bmax > 0:
                U[ist, k] = T[ist, k] / rate / bmax
            else:
                U[ist, k] = 0.0
    return U


@dataclass
class AmvaldResult:
    """Result from AMVA-LD solver."""
    Q: np.ndarray      # Queue lengths (M x K)
    U: np.ndarray      # Utilization (M x K)
    R: np.ndarray      # Response times (M x K)
    T: np.ndarray      # Throughputs (M x K)
    C: np.ndarray      # Cycle times (1 x K)
    X: np.ndarray      # System throughputs (1 x K)
    lG: float          # Log normalizing constant
    totiter: int       # Total iterations
    # see _kb/06-solver-catalog.md (MVA: "AMVA convergence flag vs iteration count")
    converged: bool = True




@_divergence_is_expected
def solver_amvald_forward(
    M: int, K: int,
    nservers: np.ndarray,
    schedparam: Optional[np.ndarray],
    lldscaling: Optional[np.ndarray],
    cdscaling: Optional[list],
    jdscaling: Optional[list],
    sched: Dict[int, int],
    classprio: Optional[np.ndarray],
    gamma: np.ndarray,
    tau: np.ndarray,
    Qchain_in: np.ndarray,
    Xchain_in: np.ndarray,
    Uchain_in: np.ndarray,
    STchain_in: np.ndarray,
    Vchain_in: np.ndarray,
    Nchain_in: np.ndarray,
    options: AmvaldOptions,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Forward step of AMVA-LD: compute waiting times from current estimates.

    Args:
        M: Number of stations
        K: Number of chains
        nservers: Number of servers at each station (M,)
        schedparam: Scheduling parameters (weights for DPS/GPS) (M x K)
        lldscaling: Load-dependent scaling matrix (M x Nmax)
        cdscaling: Class-dependent scaling functions (list of callables per station)
        jdscaling: Joint-dependent (non-product-form) scaling functions (list of callables per station)
        sched: Scheduling strategy for each station
        classprio: Class priorities (K,) - lower value = higher priority
        gamma: Correction factors
        tau: Throughput differences
        Qchain_in: Current queue length estimates (M x K)
        Xchain_in: Current throughput estimates (K,)
        Uchain_in: Current utilization estimates (M x K)
        STchain_in: Service times (M x K)
        Vchain_in: Visit ratios (M x K)
        Nchain_in: Population per chain (K,)
        options: Solver options

    Returns:
        Wchain: Waiting times (M x K)
        STeff: Effective service times (M x K)
    """
    # Initialize gamma if empty
    if gamma.size == 0:
        gamma = np.zeros((M, K))

    Nt = np.sum(Nchain_in[np.isfinite(Nchain_in)])
    if np.all(np.isinf(Nchain_in)):
        delta = 1.0
    else:
        delta = (Nt - 1) / Nt if Nt > 0 else 0.0

    # Compute deltaclass safely, handling inf and zero cases before division
    deltaclass = np.zeros_like(Nchain_in)
    valid_mask = np.isfinite(Nchain_in) & (Nchain_in > 0)
    deltaclass[valid_mask] = (Nchain_in[valid_mask] - 1) / Nchain_in[valid_mask]
    deltaclass[np.isinf(Nchain_in)] = 1.0
    deltaclass[Nchain_in == 0] = 0.0

    # Find class types
    ocl = np.where(np.isinf(Nchain_in))[0]  # open classes
    ccl = np.where(np.isfinite(Nchain_in) & (Nchain_in > 0))[0]  # closed classes
    nnzclasses = np.where(Nchain_in > 0)[0]  # non-zero classes

    # Priority groupings (lower value = higher priority): hprio = strictly
    # higher, ehprio = equal-or-higher priority than r.
    nnzclasses_hprio = {}
    nnzclasses_ehprio = {}
    if classprio is not None and len(classprio) == K:
        classprio = np.asarray(classprio).flatten()
        for r in nnzclasses:
            prio_r = classprio[r]
            nnzclasses_hprio[r] = [c for c in nnzclasses if classprio[c] < prio_r]
            nnzclasses_ehprio[r] = [c for c in nnzclasses if classprio[c] <= prio_r]
    else:
        # No priorities: all classes have same priority
        for r in nnzclasses:
            nnzclasses_hprio[r] = []
            nnzclasses_ehprio[r] = list(nnzclasses)

    # Compute arrival queue lengths
    interpTotArvlQlen = np.zeros(M)
    selfArvlQlenSeenByClosed = np.zeros((M, K))
    totArvlQlenSeenByClosed = np.zeros((M, K))
    stationaryQlen = np.zeros((M, K))
    totArvlQlenSeenByOpen = np.zeros((K, M))  # (r, k) order for open classes
    # HOL-specific: queue length seen considering priority
    totArvlQlenSeenByOpen_HOL = np.zeros((K, M))  # (r, k) order for open classes
    totArvlQlenSeenByClosed_HOL = np.zeros((M, K))

    for k in range(M):
        # sum over the non-empty classes at station k is invariant in r; hoist it
        # out of the inner loop (it was recomputed at every class below)
        sumQk = float(np.sum(Qchain_in[k, nnzclasses]))
        interpTotArvlQlen[k] = delta * sumQk
        sched_k = sched.get(k, SchedStrategy.FCFS)

        for r in nnzclasses:
            selfArvlQlenSeenByClosed[k, r] = deltaclass[r] * Qchain_in[k, r]
            totArvlQlenSeenByClosed[k, r] = (deltaclass[r] * Qchain_in[k, r] +
                                              sumQk - Qchain_in[k, r])
            stationaryQlen[k, r] = Qchain_in[k, r]

            # General totArvlQlenSeenByOpen (MATLAB line 54)
            totArvlQlenSeenByOpen[r, k] = sumQk

            # HOL-specific arrival queue lengths (only jobs with equal or higher priority)
            if sched_k == SchedStrategy.HOL:
                ehprio = nnzclasses_ehprio.get(r, list(nnzclasses))
                ehprio_without_r = [c for c in ehprio if c != r]
                # For open classes
                totArvlQlenSeenByOpen_HOL[r, k] = np.sum(Qchain_in[k, ehprio])
                # For closed classes: delta_r * Q[k,r] + sum(Q[k, ehprio - {r}])
                totArvlQlenSeenByClosed_HOL[k, r] = (deltaclass[r] * Qchain_in[k, r] +
                                                     np.sum(Qchain_in[k, ehprio_without_r]))

    # Initialize lldscaling if empty
    if lldscaling is None or lldscaling.size == 0:
        lldscaling = np.ones((M, max(1, int(np.ceil(Nt)))))

    # Compute LLD term (using proper cubic interpolation like MATLAB)
    # MATLAB: solver_amvald_forward.m lines 66-80
    _method = options.method if hasattr(options, 'method') else 'default'
    _is_lin = _method in ('lin', 'qdlin')
    if _is_lin and len(nnzclasses) > 0 and len(ccl) > 0 and gamma.ndim == 3:
        # Linearizer correction per station k uses gamma[r,k,:]; lldterm is
        # per-class (M x K) here, per-station (M,) otherwise.
        lldterm = np.ones((M, K))
        for r in nnzclasses:
            # gamma[r, :, ccl] is (|ccl|, M) under numpy advanced indexing
            gcorr = Nchain_in[ccl] @ gamma[r, :, ccl] - gamma[r, :, r]
            lldterm[:, r] = pfqn_lldfun(1 + interpTotArvlQlen + gcorr, lldscaling)
    else:
        # gamma only for closed classes: with none, the Linearizer correction is zero.
        lldterm = pfqn_lldfun(1 + interpTotArvlQlen, lldscaling)

    # Compute multi-server term based on config.multiserver setting
    # MATLAB: solver_amvald_forward.m lines 174-208
    multiserver_config = options.config.multiserver if hasattr(options, 'config') and options.config is not None else 'default'

    if multiserver_config == 'seidmann':
        # Pure Seidmann: msterm = 1/nservers for all stations
        msterm = np.ones(M) / nservers
        msterm[msterm == 0] = 1  # infinite server case
    elif multiserver_config == 'suri':
        # Suri: no demand scaling; multiserver handled via suriFactor in Wchain
        msterm = np.ones(M)
    else:
        # 'softmin' or 'default': compute msterm using pfqn_lldfun with gamma corrections
        if len(ccl) > 0 and gamma.size > 0 and Nt > 0:
            method = options.method if hasattr(options, 'method') else 'default'
            if method in ('lin', 'qdlin') and gamma.ndim == 3:
                g = np.zeros((len(ccl), M))
                for r in ccl:
                    g = g + ((Nt - 1) / Nt) * Nchain_in[r] * gamma[ccl, :, r]
                g_mean = np.mean(g, axis=0)
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + g_mean, None, nservers)
            else:
                g = np.zeros(M)
                for r in ccl:
                    g = g + (Nt - 1) * gamma[r, :]
                g_mean = np.mean(g)
                msterm = pfqn_lldfun(1 + interpTotArvlQlen + g_mean, None, nservers)
        else:
            msterm = pfqn_lldfun(1 + interpTotArvlQlen, None, nservers)

        if multiserver_config == 'default':
            # 'default': override FCFS-type stations with Seidmann approximation
            for k in range(M):
                sched_k = sched.get(k, SchedStrategy.FCFS)
                if sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR]:
                    if nservers[k] > 0 and not np.isinf(nservers[k]):
                        msterm[k] = 1.0 / nservers[k]

    # Compute class-dependent scaling term (cdterm)
    # MATLAB: solver_amvald_forward.m lines 82-104
    cdterm = np.ones((M, K))
    if cdscaling is not None and len(cdscaling) > 0:
        method = options.method if hasattr(options, 'method') else 'default'
        for r in nnzclasses:
            if np.isfinite(Nchain_in[r]):
                # Closed class: use selfArvlQlenSeenByClosed
                # MATLAB: cdterm(:,r) = pfqn_cdfun(1 + selfArvlQlenSeenByClosed + gself, cdscaling, r)
                nvec = 1 + selfArvlQlenSeenByClosed
                if method in ('lin', 'qdlin') and gamma.ndim == 3:
                    # self-fraction correction per station: entry k uses gamma[r,k,r]
                    nvec = nvec + ((Nchain_in[r] - 1) * gamma[r, :, r])[:, None]
                cdterm[:, r] = pfqn_cdfun(nvec, cdscaling, r)
            else:
                # Open class: gamma=0 avoids Inf*0=NaN from the (Nchain-1) factor.
                # MATLAB: cdterm(:,r) = pfqn_cdfun(1 + stationaryQlen, cdscaling, r)
                cdterm[:, r] = pfqn_cdfun(1 + stationaryQlen, cdscaling, r)

    # Compute joint-dependent scaling term (jdterm), non-product-form eta_i.
    # Evaluated identically to cdterm; kept separate to preserve the
    # product-form vs joint distinction. cd and jd are mutually exclusive per
    # station, so at most one of cdterm/jdterm differs from 1 at any station.
    jdterm = np.ones((M, K))
    if jdscaling is not None and len(jdscaling) > 0:
        method = options.method if hasattr(options, 'method') else 'default'
        for r in nnzclasses:
            if np.isfinite(Nchain_in[r]):
                nvec = 1 + selfArvlQlenSeenByClosed
                if method in ('lin', 'qdlin') and gamma.ndim == 3:
                    nvec = nvec + ((Nchain_in[r] - 1) * gamma[r, :, r])[:, None]
                jdterm[:, r] = pfqn_jdfun(nvec, jdscaling, r)
            else:
                jdterm[:, r] = pfqn_jdfun(1 + stationaryQlen, jdscaling, r)

    # Compute Suri correction factor per station
    suriFactor = np.ones(M)
    if multiserver_config == 'suri':
        suri_alpha = 4.464
        suri_beta = 0.676
        for k in range(M):
            c = nservers[k]
            if c > 1 and np.isfinite(c):
                rho_k = min(np.nansum(Uchain_in[k, :]) / c, 1 - options.tol)
                if rho_k > 0:
                    suriFactor[k] = rho_k ** (suri_alpha * (c ** suri_beta - 1)) / c
                else:
                    suriFactor[k] = 1.0 / c
            elif np.isinf(c) or c == 0:
                suriFactor[k] = 0.0

    # Compute effective service times
    # MATLAB: STeff(k,r) = STchain_in(k,r) * lldterm(k,r) * msterm(k) * cdterm(k,r)
    Wchain = np.zeros((M, K))
    STeff = np.zeros((M, K))

    # lldterm shape differs by path (see above); index accordingly, no tiled copy.
    _lld2d = lldterm.ndim == 2
    for r in nnzclasses:
        for k in range(M):
            lld_kr = lldterm[k, r] if _lld2d else lldterm[k]
            STeff[k, r] = STchain_in[k, r] * lld_kr * msterm[k] * cdterm[k, r] * jdterm[k, r]

    # Port of solver_amvald_forward.m:174-240 (Wang-Sevcik qli/fli); must run
    # after STeff is known. MATLAB quirks reproduced verbatim: qliden double-
    # counts INF stations; HOL arm sums Qchain including class r itself.
    if options.method in ('qli', 'fli'):
        infset = [k for k in range(M)
                  if sched.get(k, SchedStrategy.FCFS) == SchedStrategy.INF
                  or np.isinf(nservers[k])]
        for k in range(M):
            sched_k = sched.get(k, SchedStrategy.FCFS)
            for r in nnzclasses:
                if sched_k == SchedStrategy.HOL:
                    cls = nnzclasses_ehprio.get(r, list(nnzclasses))
                else:
                    cls = nnzclasses
                sumQk_cls = float(np.sum(Qchain_in[k, cls]))
                if Nchain_in[r] == 1:
                    totArvlQlenSeenByClosed[k, r] = sumQk_cls - Qchain_in[k, r]
                else:
                    qlinum = STeff[k, r] * (1 + sumQk_cls - Qchain_in[k, r])
                    qliden = float(np.sum(STeff[infset, r])) if len(infset) > 0 else 0.0
                    for m in range(M):
                        qliden += STeff[m, r] * (1 + float(np.sum(Qchain_in[m, cls])) - Qchain_in[m, r])
                    ratio = (qlinum / qliden) if qliden != 0 else 0.0
                    if options.method == 'qli':
                        totArvlQlenSeenByClosed[k, r] = (
                            sumQk_cls - (1.0 / (Nchain_in[r] - 1)) * (Qchain_in[k, r] - ratio))
                    else:  # 'fli'
                        totArvlQlenSeenByClosed[k, r] = (
                            sumQk_cls - (2.0 / Nchain_in[r]) * Qchain_in[k, r] + ratio)

    # Compute waiting times based on scheduling strategy
    for ir, r in enumerate(nnzclasses):
        # other classes; nnzclasses is already sorted and unique, so a boolean
        # mask is an exact but far cheaper substitute for np.setdiff1d here
        sd = nnzclasses[nnzclasses != r]

        for k in range(M):
            sched_k = sched.get(k, SchedStrategy.FCFS)

            # CRITICAL: Treat infinite-server stations as delay (INF) regardless of sched
            # This prevents inf values when FCFS/PS stations have nservers = inf
            if np.isinf(nservers[k]):
                sched_k = SchedStrategy.INF

            if sched_k == SchedStrategy.EXT:
                # External source station - no queueing, zero waiting time
                Wchain[k, r] = 0.0

            elif sched_k == SchedStrategy.INF:
                # Infinite server (delay)
                Wchain[k, r] = STeff[k, r]

            elif sched_k == SchedStrategy.PS:
                # Processor sharing (MATLAB solver_amvald_forward.m lines 303-342)
                method = options.method if hasattr(options, 'method') else 'default'

                if method in ('qd', 'lin', 'qdlin'):
                    if multiserver_config == 'seidmann':
                        # Seidmann: add (nservers-1) correction (MATLAB lines 307-313)
                        ns = nservers[k]
                        Wchain[k, r] = STeff[k, r] * (ns - 1) if ns > 1 and not np.isinf(ns) else 0.0
                        if r in ocl:
                            Wchain[k, r] += STeff[k, r] * (1 + totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0)
                        else:
                            if method in ('lin', 'qdlin') and gamma.ndim == 3 and len(ccl) > 0:
                                gamma_correction = np.dot(Nchain_in[ccl], gamma[r, k, ccl]) - gamma[r, k, r]
                            else:
                                gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            wait_factor = max(options.tol, 1 + interpTotArvlQlen[k] + gamma_correction)
                            Wchain[k, r] += STeff[k, r] * wait_factor
                    elif multiserver_config == 'suri':
                        # Suri: W = S * (1 + Q_seen * suriFactor)
                        if r in ocl:
                            Wchain[k, r] = STeff[k, r] * (1 + (totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0) * suriFactor[k])
                        else:
                            if method in ('lin', 'qdlin') and gamma.ndim == 3 and len(ccl) > 0:
                                gamma_correction = np.dot(Nchain_in[ccl], gamma[r, k, ccl]) - gamma[r, k, r]
                            else:
                                gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            Wchain[k, r] = STeff[k, r] * (1 + (interpTotArvlQlen[k] + gamma_correction) * suriFactor[k])
                    else:
                        # Default/softmin: no (nservers-1) correction (MATLAB lines 314-324)
                        if r in ocl:
                            Wchain[k, r] = STeff[k, r] * (1 + (totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0))
                        else:
                            if method in ('lin', 'qdlin') and gamma.ndim == 3 and len(ccl) > 0:
                                gamma_correction = np.dot(Nchain_in[ccl], gamma[r, k, ccl]) - gamma[r, k, r]
                            else:
                                gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            wait_factor = max(options.tol, 1 + interpTotArvlQlen[k] + gamma_correction)
                            Wchain[k, r] = STeff[k, r] * wait_factor
                else:
                    # Other methods (MATLAB lines 326-341)
                    if multiserver_config == 'seidmann':
                        ns = nservers[k]
                        Wchain[k, r] = STeff[k, r] * (ns - 1) if ns > 1 and not np.isinf(ns) else 0.0
                        if r in ocl:
                            Wchain[k, r] += STeff[k, r] * (1 + (totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0))
                        else:
                            gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            Wchain[k, r] += STeff[k, r] * (1 + totArvlQlenSeenByClosed[k, 0] + gamma_correction)
                    elif multiserver_config == 'suri':
                        # Suri: W = S * (1 + Q_seen * suriFactor)
                        if r in ocl:
                            Wchain[k, r] = STeff[k, r] * (1 + (totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0) * suriFactor[k])
                        else:
                            gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            Wchain[k, r] = STeff[k, r] * (1 + (totArvlQlenSeenByClosed[k, 0] + gamma_correction) * suriFactor[k])
                    else:
                        if r in ocl:
                            Wchain[k, r] = STeff[k, r] * (1 + (totArvlQlenSeenByOpen[r, k] if totArvlQlenSeenByOpen is not None else 0))
                        else:
                            gamma_correction = (Nt - 1) * gamma[r, k] if gamma.ndim == 2 and r < gamma.shape[0] and k < gamma.shape[1] else 0.0
                            Wchain[k, r] = STeff[k, r] * (1 + totArvlQlenSeenByClosed[k, 0] + gamma_correction)

            elif sched_k == SchedStrategy.DPS:
                # DPS: matches solver_amvald_forward.m:383-403
                w = schedparam if schedparam is not None else np.ones((M, K))

                # Multi-server DPS (Seidmann): STeff already has 1/nservers(k);
                # serial think-time correction is STeff*(nservers-1).
                w_base = STeff[k, r] * (nservers[k] - 1) if nservers[k] > 1 else 0.0

                # Base response time for own class (arrival sees itself)
                if r in ocl:
                    Wchain[k, r] = w_base + STeff[k, r] * (1 + stationaryQlen[k, r])
                else:
                    Wchain[k, r] = w_base + STeff[k, r] * (1 + selfArvlQlenSeenByClosed[k, r])

                # Add slowdown from other classes (sd = other classes)
                for s in sd:
                    w_s = w[k, s] if k < w.shape[0] and s < w.shape[1] else 1.0
                    w_r = w[k, r] if k < w.shape[0] and r < w.shape[1] else 1.0

                    if w_s == w_r:
                        # Equal weights: treat like PS
                        Wchain[k, r] += STeff[k, r] * stationaryQlen[k, s]
                    elif w_r > 0:
                        # Weight ratio scaling: contribution proportional to w_s/w_r
                        Wchain[k, r] += STeff[k, r] * stationaryQlen[k, s] * w_s / w_r

            elif sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.LCFSPR]:
                # FCFS-type stations (MATLAB solver_amvald_forward.m lines 366-438)
                if STeff[k, r] > 0:
                    ns = nservers[k]

                    # Compute Uchain_r with tau correction (MATLAB line 368)
                    Uchain_r = np.zeros_like(Uchain_in)
                    for s_idx in range(K):
                        if Xchain_in[s_idx] > 0:
                            Uchain_r[:, s_idx] = Uchain_in[:, s_idx] / Xchain_in[s_idx] * (Xchain_in[s_idx] + tau[r, s_idx])
                        else:
                            Uchain_r[:, s_idx] = Uchain_in[:, s_idx]

                    # Compute multiserver correction Bk (MATLAB lines 370-380)
                    if ns > 1 and not np.isinf(ns):
                        deltaclass_r = np.ones(K)
                        deltaclass_r[r] = deltaclass[r]
                        rho_sum = np.sum(deltaclass_r * Xchain_in * Vchain_in[k, :] * STeff[k, :])
                        if rho_sum < 0.75:
                            Bk = deltaclass_r * Xchain_in * Vchain_in[k, :] * STeff[k, :]
                        else:
                            Bk = (deltaclass_r * Xchain_in * Vchain_in[k, :] * STeff[k, :]) ** (ns - 1)
                        Bk = np.where(np.isfinite(Bk), Bk, 1.0)
                    else:
                        Bk = np.ones(K)

                    # Check for single-server with load-dependent scaling (MATLAB lines 382-405)
                    has_lld = lldscaling is not None and np.any(lldscaling != 0) if lldscaling is not None else False
                    has_cd = cdscaling is not None and len(cdscaling) > 0
                    if ns == 1 and (has_lld or has_cd):
                        # Single-server LD: simple formula (MATLAB lines 390-404)
                        Wchain[k, r] = STeff[k, r]
                        if r in ocl:
                            Wchain[k, r] += (STeff[k, r] * stationaryQlen[k, r] +
                                            np.sum(STeff[k, sd] * stationaryQlen[k, sd]))
                        else:
                            Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r] +
                                            np.sum(STeff[k, sd] * stationaryQlen[k, sd]))
                    else:
                        # Multiserver or non-LD: depends on config (MATLAB lines 407-437)
                        if multiserver_config == 'softmin':
                            # Softmin: no (ns-1) serial think time term (MATLAB lines 408-421)
                            Wchain[k, r] = STeff[k, r]
                            if r in ocl:
                                Wchain[k, r] += (STeff[k, r] * stationaryQlen[k, r] * Bk[r] +
                                                np.sum(STeff[k, sd] * stationaryQlen[k, sd] * Bk[sd]))
                            else:
                                Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r] * Bk[r] +
                                                np.sum(STeff[k, sd] * stationaryQlen[k, sd] * Bk[sd]))
                        elif multiserver_config == 'suri':
                            # Suri multiserver: W = S * (1 + L_m * suriFactor)
                            Wchain[k, r] = STeff[k, r]
                            if r in ocl:
                                L_m = deltaclass[r] * stationaryQlen[k, r] + np.sum(stationaryQlen[k, sd])
                                Wchain[k, r] += STeff[k, r] * L_m * suriFactor[k]
                            else:
                                L_m = selfArvlQlenSeenByClosed[k, r] + np.sum(stationaryQlen[k, sd])
                                Wchain[k, r] += STeff[k, r] * L_m * suriFactor[k]
                        else:
                            # Seidmann/default: has (ns-1) correction (MATLAB lines 422-436)
                            Wchain[k, r] = STeff[k, r] * (ns - 1) if ns > 1 else 0.0
                            Wchain[k, r] += STeff[k, r]
                            if r in ocl:
                                Wchain[k, r] += (STeff[k, r] * deltaclass[r] * stationaryQlen[k, r] * Bk[r] +
                                                np.sum(STeff[k, sd] * Bk[sd] * stationaryQlen[k, sd]))
                            else:
                                Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r] * Bk[r] +
                                                np.sum(STeff[k, sd] * Bk[sd] * stationaryQlen[k, sd]))

            elif sched_k == SchedStrategy.HOL:
                # Head of Line (Priority) scheduling - matches MATLAB solver_amvald_forward.m lines 441-547
                if STeff[k, r] > 0:
                    ns = nservers[k]
                    hprio = nnzclasses_hprio.get(r, [])

                    # Chandy-Lakshmi priority scaling (solver_amvald_forward.m:446-458)
                    tol = options.tol if hasattr(options, 'tol') and options.tol else 1e-6
                    UHigherPrio = 0.0
                    if len(hprio) > 0:
                        for h in hprio:
                            # MATLAB tau(h) uses linear indexing on K×K matrix, which is tau(h,1) = first column
                            # In Python 0-indexed: tau[h, 0]
                            tau_h = tau[h, 0] if tau.ndim == 2 and h < tau.shape[0] else 0.0
                            UHigherPrio += Vchain_in[k, h] * STeff[k, h] * (Xchain_in[h] - Qchain_in[k, h] * tau_h)

                    # MATLAB: prioScaling = min([max([options.tol,1-UHigherPrio]),1-options.tol])
                    prioScaling = max(tol, 1.0 - UHigherPrio)
                    prioScaling = min(prioScaling, 1.0 - tol)

                    if ns == 1 or np.isinf(ns):
                        # Single server: matches MATLAB lines 494-506
                        # Wchain = STeff / prioScaling
                        Wchain[k, r] = STeff[k, r] / prioScaling  # Base
                        if r in ocl:
                            # Open class: use stationaryQlen
                            Wchain[k, r] += (STeff[k, r] * stationaryQlen[k, r]) / prioScaling
                        else:
                            # Closed class: use selfArvlQlenSeenByClosed
                            Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r]) / prioScaling
                    else:
                        # Multi-server (Seidmann, matches lines 527-542); use
                        # deltaclass, not deltaclass_r, for HOL priority.
                        if np.sum(deltaclass * Xchain_in * Vchain_in[k, :] * STeff[k, :]) < 0.75:
                            # Light load case: Bk = utilization / ns (seidmann)
                            Bk = (deltaclass * Xchain_in * Vchain_in[k, :] * STeff[k, :]) / ns
                        else:
                            # High load case: power of ns (seidmann)
                            Bk = ((deltaclass * Xchain_in * Vchain_in[k, :] * STeff[k, :]) / ns) ** ns
                        Bk = np.where(np.isfinite(Bk), Bk, 1.0)

                        # Multi-server correction with serial think time (MATLAB lines 528-529)
                        # (1/nservers term already in STeff via msterm)
                        Wchain[k, r] = STeff[k, r] * (ns - 1) / prioScaling
                        Wchain[k, r] += STeff[k, r] / prioScaling  # Base

                        if r in ocl:
                            # Open class (MATLAB line 532)
                            Wchain[k, r] += (STeff[k, r] * stationaryQlen[k, r] * Bk[r]) / prioScaling
                        else:
                            # Closed class (MATLAB line 542)
                            Wchain[k, r] += (STeff[k, r] * selfArvlQlenSeenByClosed[k, r] * Bk[r]) / prioScaling

            else:
                # Default: same as FCFS
                Wchain[k, r] = STeff[k, r] * (1 + selfArvlQlenSeenByClosed[k, r])

    return Wchain, STeff


@_divergence_is_expected
def solver_amvald(
    sn,  # NetworkStruct
    Lchain: np.ndarray,
    STchain: np.ndarray,
    Vchain: np.ndarray,
    alpha: np.ndarray,
    Nchain: np.ndarray,
    SCVchain: np.ndarray,
    refstatchain: np.ndarray,
    options: Optional[AmvaldOptions] = None
) -> AmvaldResult:
    """
    Approximate MVA for Load-Dependent networks (AMVA-LD).

    Implements iterative approximate MVA with under-relaxation for
    load-dependent queueing networks.

    Args:
        sn: Network structure
        Lchain: Demands per chain (M x K)
        STchain: Service times per chain (M x K)
        Vchain: Visit ratios per chain (M x K)
        alpha: Class-to-chain mapping
        Nchain: Population per chain (K,)
        SCVchain: Squared coefficient of variation (M x K)
        refstatchain: Reference station per chain (K,)
        options: Solver options

    Returns:
        AmvaldResult with Q, U, R, T, C, X, lG, totiter
    """
    if options is None:
        options = AmvaldOptions()

    M = sn.nstations
    K = sn.nchains

    Nchain = np.asarray(Nchain).flatten()
    Nt = np.sum(Nchain[np.isfinite(Nchain)])

    # Method selection: convert 'default' to appropriate method
    # (MATLAB solver_amva.m lines 46-55)
    nservers_arr = np.asarray(sn.nservers).flatten() if sn.nservers is not None else np.ones(M)
    method = options.method
    if method in ('default', 'amva'):
        if Nt <= 2 or np.any(Nchain[np.isfinite(Nchain)] < 1):
            method = 'qd'
        else:
            finite_servers = nservers_arr[np.isfinite(nservers_arr)]
            max_servers = int(np.max(finite_servers)) if len(finite_servers) > 0 else 1
            if max_servers == 1:
                method = 'egflin'  # Single server model
            else:
                method = 'lin'     # Multiserver model
    # Update options.method for gamma/tau initialization check
    options = AmvaldOptions(
        method=method,
        iter_tol=options.iter_tol,
        iter_max=options.iter_max,
        tol=options.tol,
        verbose=options.verbose,
        init_sol=getattr(options, 'init_sol', None),
        config=options.config
    )

    delta = (Nt - 1) / Nt if Nt > 0 else 0.0
    # Compute deltaclass safely, handling inf and zero cases before division
    deltaclass = np.zeros_like(Nchain)
    valid_mask = np.isfinite(Nchain) & (Nchain > 0)
    deltaclass[valid_mask] = (Nchain[valid_mask] - 1) / Nchain[valid_mask]
    deltaclass[np.isinf(Nchain)] = 1.0
    deltaclass[Nchain <= 0] = 0.0

    tol = options.iter_tol
    nservers = np.asarray(sn.nservers).flatten() if sn.nservers is not None else np.ones(M)
    schedparam = np.asarray(sn.schedparam) if hasattr(sn, 'schedparam') and sn.schedparam is not None else None
    lldscaling = sn.lldscaling if hasattr(sn, 'lldscaling') and sn.lldscaling is not None else None
    cdscaling = sn.cdscaling if hasattr(sn, 'cdscaling') and sn.cdscaling is not None else None
    jdscaling = sn.jdscaling if hasattr(sn, 'jdscaling') and sn.jdscaling is not None else None

    # sn.sched can be dict/ndarray/None; normalize SchedStrategy values since
    # they may come from a different enum class instance (different import path).
    if hasattr(sn, 'sched') and sn.sched is not None:
        if isinstance(sn.sched, dict):
            # Normalize all values to the local SchedStrategy enum
            sched = {k: _normalize_sched_strategy(v) for k, v in sn.sched.items()}
        elif hasattr(sn.sched, '__iter__') and not isinstance(sn.sched, str):
            # Convert array of strings to dict
            sched = {}
            for i, s in enumerate(sn.sched):
                if s is not None:
                    if isinstance(s, str):
                        # Map string names to SchedStrategy enum values
                        sched_map = {
                            'INF': SchedStrategy.INF,
                            'FCFS': SchedStrategy.FCFS,
                            'SIRO': SchedStrategy.SIRO,
                            'PS': SchedStrategy.PS,
                            'DPS': SchedStrategy.DPS,
                            'GPS': SchedStrategy.GPS,
                            'HOL': SchedStrategy.HOL,
                            'LCFSPR': SchedStrategy.LCFSPR,
                            'EXT': SchedStrategy.EXT,
                        }
                        sched[i] = sched_map.get(s, SchedStrategy.FCFS)
                    else:
                        sched[i] = s
        else:
            sched = {}
    else:
        sched = {}

    classprio = np.asarray(sn.classprio).flatten() if hasattr(sn, 'classprio') and sn.classprio is not None else None

    # Initialize Q, X, U
    Uchain = np.zeros((M, K))
    Tchain = np.zeros((M, K))

    # Warm-start initialization if a valid hint is available
    Qchain = None
    if getattr(options, 'init_sol', None) is not None:
        init_sol = np.asarray(options.init_sol, dtype=float)
        if init_sol.shape == (M, K):
            Qchain = init_sol.copy()
        # else: stale warm-start hint from a different station/chain basis

    if Qchain is None:
        # Balanced initialization for queue lengths
        Qchain = np.ones((M, K))
        Qchain = Qchain / np.sum(Qchain, axis=0, keepdims=True) * Nchain.reshape(1, -1)
        Qchain[np.isinf(Qchain)] = 0.0
        Qchain[:, Nchain == 0] = 0.0

        # For open classes, zero out at reference station
        for r in range(K):
            if np.isinf(Nchain[r]):
                refstat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
                if refstat < M:
                    Qchain[refstat, r] = 0.0

    nnzclasses = np.where(Nchain > 0)[0]
    Xchain = 1.0 / np.sum(STchain, axis=0)
    Xchain[np.isinf(Xchain)] = 0.0

    # For open classes, use arrival rate
    for r in range(K):
        if np.isinf(Nchain[r]):
            refstat = int(refstatchain[r, 0]) if refstatchain.ndim > 1 else int(refstatchain[r])
            if refstat < M and STchain[refstat, r] > 0:
                Xchain[r] = 1.0 / STchain[refstat, r]

    # Initialize utilization
    for k in range(M):
        for r in nnzclasses:
            if np.isinf(nservers[k]):
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[r]
            else:
                Uchain[k, r] = Vchain[k, r] * STchain[k, r] * Xchain[r] / nservers[k]

    # gamma is 3D (K,M,K) for 'lin'/'qdlin' (class-based); 2D (K,M) otherwise.
    if options.method in ('lin', 'qdlin'):
        gamma = np.zeros((K, M, K))  # 3D for class-based corrections
    else:
        gamma = np.zeros((K, M))  # 2D for total corrections
    tau = np.zeros((K, K))

    # Main iteration loop
    omicron = 0.5  # under-relaxation parameter
    totiter = 0
    outer_iter = 0
    max_outer_iter = np.sqrt(options.iter_max)
    # Match MATLAB's max_totiter cap (solver_amvald.m line 3)
    max_totiter = min(options.iter_max, 10000)

    QchainOuter_1 = Qchain + np.inf

    # MATLAB: while (outer_iter < 2 || max(max(abs(Qchain-QchainOuter_1))) > tol) && outer_iter < sqrt(options.iter_max) && totiter <= max_totiter
    while (outer_iter < 2 or np.max(np.abs(Qchain - QchainOuter_1)) > tol) and outer_iter < max_outer_iter and totiter <= max_totiter:
        outer_iter += 1
        QchainOuter_1 = Qchain.copy()
        XchainOuter_1 = Xchain.copy()
        UchainOuter_1 = Uchain.copy()

        # For 'lin'/'qdlin' methods, first iterate at population N-1_s for each class s
        # to compute gamma and tau corrections (MATLAB solver_amvald.m lines 76-142)
        if np.isfinite(Nt) and Nt > 0 and options.method in ('lin', 'qdlin'):
            for s in range(K):
                if np.isfinite(Nchain[s]) and Nchain[s] > 0:  # Don't recur on open classes
                    # Create population with one less job in class s
                    Nchain_s = Nchain.copy()
                    Nchain_s[s] = Nchain[s] - 1
                    Nt_s = Nt - 1

                    # Initialize at N-1_s by scaling from population N
                    scale = (Nt - 1) / Nt if Nt > 0 else 0.0
                    Qchain_s = Qchain * scale
                    Xchain_s = Xchain * scale
                    Uchain_s = Uchain * scale

                    # Iterate at population N-1_s until convergence
                    iter_s = 0
                    Qchain_s_1 = Qchain_s + np.inf
                    while (iter_s < 2 or np.max(np.abs(Qchain_s - Qchain_s_1)) > tol) and iter_s <= max_outer_iter:
                        iter_s += 1

                        Qchain_s_1 = Qchain_s.copy()
                        Xchain_s_1 = Xchain_s.copy()
                        Uchain_s_1 = Uchain_s.copy()

                        # Forward step at N-1_s
                        Wchain_s, STeff_s = solver_amvald_forward(
                            M, K, nservers, schedparam, lldscaling, cdscaling, jdscaling, sched, classprio, gamma, tau,
                            Qchain_s_1, Xchain_s_1, Uchain_s_1, STchain, Vchain, Nchain_s, options
                        )

                        totiter += 1
                        if totiter >= max_totiter:
                            break

                        # Update metrics at N-1_s
                        Cchain_s_val = np.zeros(K)
                        for r in nnzclasses:
                            if np.sum(Wchain_s[:, r]) == 0:
                                Xchain_s[r] = 0.0
                            else:
                                if np.isinf(Nchain_s[r]):
                                    Cchain_s_val[r] = np.dot(Vchain[:, r], Wchain_s[:, r])
                                elif Nchain_s[r] == 0:
                                    Xchain_s[r] = 0.0
                                    Cchain_s_val[r] = 0.0
                                else:
                                    Cchain_s_val[r] = np.dot(Vchain[:, r], Wchain_s[:, r])
                                    # Protect against division by zero/very small values
                                    if Cchain_s_val[r] > 1e-14:
                                        Xchain_s[r] = omicron * Nchain_s[r] / Cchain_s_val[r] + (1 - omicron) * Xchain_s_1[r]
                                    else:
                                        # If C is effectively zero, keep previous X value (no update)
                                        Xchain_s[r] = Xchain_s_1[r]

                            for k in range(M):
                                Qchain_s[k, r] = omicron * Xchain_s[r] * Vchain[k, r] * Wchain_s[k, r] + (1 - omicron) * Qchain_s_1[k, r]
                                # MATLAB uses STeff_s for utilization update (line 119 in solver_amvald.m)
                                # Uchain_s(k,r) = omicron * Vchain(k,r) * STeff_s(k,r) * Xchain_s(r) + (1-omicron) * Uchain_s_1(k,r)
                                Uchain_s[k, r] = omicron * Vchain[k, r] * STeff_s[k, r] * Xchain_s[r] + (1 - omicron) * Uchain_s_1[k, r]

                    # Compute gamma and tau from converged N-1_s results
                    if options.method == 'lin':
                        # 3D gamma: gamma(s, k, r) = Qchain_s_1(k,r)/Nchain_s(r) - QchainOuter_1(k,r)/Nchain(r)
                        for k in range(M):
                            for r in nnzclasses:
                                if not np.isinf(Nchain[r]) and Nchain_s[r] > 0:
                                    gamma[s, k, r] = Qchain_s_1[k, r] / Nchain_s[r] - QchainOuter_1[k, r] / Nchain[r]
                    else:
                        # see _kb/06-solver-catalog.md (MVA: "In solver_amvald,
                        # only method='lin' takes the per-class gamma correction")
                        for k in range(M):
                            gval = (np.sum(Qchain_s_1[k, :]) / (Nt - 1) - np.sum(QchainOuter_1[k, :]) / Nt) if Nt > 1 else 0.0
                            if gamma.ndim == 3:
                                gamma[s, k, 0] = gval
                            else:
                                gamma[s, k] = gval

                    # tau update
                    for r in nnzclasses:
                        tau[s, r] = Xchain_s_1[r] - XchainOuter_1[r]

                    # Check if iteration limit exceeded (break out of for s loop)
                    # MATLAB solver_amvald.m lines 143-147
                    if totiter >= max_totiter:
                        break

            # Check if iteration limit exceeded (break out of outer while loop)
            # MATLAB solver_amvald.m lines 150-153
            if totiter >= max_totiter:
                break

        # Inner iteration at population N
        inner_iter = 0
        Qchain_1 = Qchain + np.inf

        while (inner_iter < 2 or np.max(np.abs(Qchain - Qchain_1)) > tol) and inner_iter <= max_outer_iter:
            inner_iter += 1

            Qchain_1 = Qchain.copy()
            Xchain_1 = Xchain.copy()
            Uchain_1 = Uchain.copy()

            # Forward step
            Wchain, STeff = solver_amvald_forward(
                M, K, nservers, schedparam, lldscaling, cdscaling, jdscaling, sched, classprio, gamma, tau,
                Qchain_1, Xchain_1, Uchain_1, STchain, Vchain, Nchain, options
            )

            totiter += 1
            if totiter >= max_totiter:
                break

            # Update metrics
            Cchain_s = np.zeros(K)
            Rchain = np.zeros((M, K))

            for r in nnzclasses:
                if np.sum(Wchain[:, r]) == 0:
                    Xchain[r] = 0.0
                else:
                    if np.isinf(Nchain[r]):
                        Cchain_s[r] = np.dot(Vchain[:, r], Wchain[:, r])
                        # X remains constant for open classes
                    elif Nchain[r] == 0:
                        Xchain[r] = 0.0
                        Cchain_s[r] = 0.0
                    else:
                        Cchain_s[r] = np.dot(Vchain[:, r], Wchain[:, r])
                        # Protect against division by zero/very small values
                        if Cchain_s[r] > 1e-14:
                            Xchain[r] = omicron * Nchain[r] / Cchain_s[r] + (1 - omicron) * Xchain_1[r]
                        else:
                            # If C is effectively zero, keep previous X value (no update)
                            Xchain[r] = Xchain_1[r]

                for k in range(M):
                    Rchain[k, r] = Vchain[k, r] * Wchain[k, r]
                    Qchain[k, r] = omicron * Xchain[r] * Vchain[k, r] * Wchain[k, r] + (1 - omicron) * Qchain_1[k, r]
                    Tchain[k, r] = Xchain[r] * Vchain[k, r]
                    # MATLAB uses STeff for utilization update (line 191 in solver_amvald.m)
                    # Uchain(k,r) = omicron * Vchain(k,r) * STeff(k,r) * Xchain(r) + (1-omicron) * Uchain_1(k,r)
                    Uchain[k, r] = omicron * Vchain[k, r] * STeff[k, r] * Xchain[r] + (1 - omicron) * Uchain_1[k, r]

    # Utilization capping (matches solver_amvald.m:200-210): redistribute
    # proportionally by STeff if total utilization exceeds 1.
    for k in range(M):
        sched_k = sched.get(k, SchedStrategy.FCFS)
        if sched_k in [SchedStrategy.FCFS, SchedStrategy.SIRO, SchedStrategy.PS, SchedStrategy.LCFSPR, SchedStrategy.DPS, SchedStrategy.HOL]:
            U_sum = np.sum(Uchain[k, :])
            if U_sum > 1:
                for r in range(K):
                    if Vchain[k, r] * STeff[k, r] > 0:
                        # denom = (Vchain(k,:) .* STeff(k,:)) * Xchain(:) - inner product
                        denom = np.sum(Vchain[k, :] * STeff[k, :] * Xchain)
                        if denom > 0:
                            Uchain[k, r] = min(1.0, U_sum) * Vchain[k, r] * STeff[k, r] * Xchain[r] / denom

    # Compute final response times
    Rchain = np.zeros((M, K))
    for k in range(M):
        for r in range(K):
            if Tchain[k, r] > 0:
                Rchain[k, r] = Qchain[k, r] / Tchain[k, r]

    # Clean up non-finite values
    Xchain = np.where(np.isfinite(Xchain), Xchain, 0.0)
    Uchain = np.where(np.isfinite(Uchain), Uchain, 0.0)
    Rchain = np.where(np.isfinite(Rchain), Rchain, 0.0)

    # Zero out results for zero population chains
    for j in range(K):
        if Nchain[j] == 0:
            Xchain[j] = 0.0
            Uchain[:, j] = 0.0
            Qchain[:, j] = 0.0
            Rchain[:, j] = 0.0
            Tchain[:, j] = 0.0

    # Compute cycle times
    C_result = np.sum(Rchain, axis=0).reshape(1, -1)

    # Estimate normalizing constant
    ccl = np.where(np.isfinite(Nchain))[0]
    Nclosed = Nchain[ccl]
    Xclosed = Xchain[ccl]
    valid = Xclosed > options.tol
    if np.any(valid):
        lG = -np.sum(Nclosed[valid] * np.log(Xclosed[valid]))
    else:
        lG = np.nan

    return AmvaldResult(
        Q=Qchain,
        U=Uchain,
        R=Rchain,
        T=Tchain,
        C=C_result,
        X=Xchain.reshape(1, -1),
        lG=lG,
        totiter=totiter,
        # see _kb/06-solver-catalog.md (MVA: "AMVA convergence flag vs iteration count")
        converged=bool(np.max(np.abs(Qchain - QchainOuter_1)) <= tol)
    )


__all__ = [
    'solver_amvald',
    'solver_amvald_forward',
    'solver_amvald_cd_peak_post',
    'solver_amvald_jd_peak_post',
    'AmvaldOptions',
    'AmvaldResult',
    'SchedStrategy',
]
