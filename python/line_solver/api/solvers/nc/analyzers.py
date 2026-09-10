"""
NC Solver analyzers.

Native Python implementation of NC solver analyzers that orchestrate the
core handlers and provide interpolation for non-integer populations.

Port from:


"""

import numpy as np
from math import floor, ceil
from typing import Optional, Tuple
from dataclasses import dataclass, field
import time

from ...sn import NetworkStruct
from .handler import solver_nc, solver_ncld, SolverNCReturn, SolverNCLDReturn, SolverOptions


# Fine tolerance for numerical comparisons
FINE_TOL = 1e-12


@dataclass
class NCResultProb:
    """Probability results from NC solver."""
    logNormConstAggr: Optional[float] = None
    marginal: Optional[np.ndarray] = None
    joint: Optional[float] = None
    itemProb: Optional[np.ndarray] = None


@dataclass
class NCResult:
    """
    Result of NC solver analysis.

    Attributes:
        QN: Mean queue lengths (M x K)
        UN: Utilizations (M x K)
        RN: Response times (M x K)
        TN: Throughputs (M x K)
        CN: Cycle times (1 x K)
        XN: System throughputs (1 x K)
        lG: Log normalizing constant
        STeff: Effective service times
        it: Number of iterations
        runtime: Runtime in seconds
        method: Method used
        solver: Solver name
        prob: Probability results
    """
    QN: Optional[np.ndarray] = None
    UN: Optional[np.ndarray] = None
    RN: Optional[np.ndarray] = None
    TN: Optional[np.ndarray] = None
    CN: Optional[np.ndarray] = None
    XN: Optional[np.ndarray] = None
    lG: float = 0.0
    STeff: Optional[np.ndarray] = None
    it: int = 0
    runtime: float = 0.0
    method: str = ""
    solver: str = "NC"
    prob: NCResultProb = field(default_factory=NCResultProb)


def solver_nc_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    Main NC solver analyzer.

    Performs NC analysis with interpolation for non-integer populations.
    If population values are not integers, interpolates between floor and
    ceiling values.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverOptions()

    # Check for multiserver with exact method and open classes.
    # The refusal is for an OPEN or MIXED multiserver model, as the message says: a
    # CLOSED one has an exact answer and MATLAB's solver_nc_analyzer.m:78 and C++'s
    # nc_dispatch.h:161 both gate on it. This copy omitted the open test that its own
    # ncld twin below carries, and so refused closed multiserver models too. The same
    # omission was live in Solver_nc_analyzer.java, where config.multiserver='seidmann'
    # first gave 'exact' a route to this analyzer with a multiserver station still
    # un-converted. See _kb/06-solver-catalog.md (NC section)
    nservers = sn.nservers
    if nservers is not None:
        nservers_finite = nservers.copy()
        nservers_finite = nservers_finite[np.isfinite(nservers_finite)]
        has_infinite_jobs = sn.njobs is not None and np.any(np.isinf(sn.njobs))
        if (len(nservers_finite) > 0 and np.max(nservers_finite) > 1 and
                has_infinite_jobs and options.method == 'exact'):
            raise RuntimeError(
                "NC solver cannot provide exact solutions for open or mixed queueing networks. "
                "Remove the 'exact' option."
            )

    # Create floor/ceiling copies for non-integer populations
    njobs = sn.njobs
    if njobs is None:
        njobs = np.zeros((1, sn.nclasses))

    njobs_floor = np.floor(njobs)
    njobs_ceil = np.ceil(njobs)
    eta = np.abs(njobs - njobs_floor)

    # Check for non-integer populations
    non_integer_job = np.any(eta > FINE_TOL)

    result = NCResult()

    if non_integer_job:
        # Interpolate between floor and ceiling
        sn_floor = sn.copy()
        sn_ceil = sn.copy()
        sn_floor.njobs = njobs_floor
        sn_ceil.njobs = njobs_ceil

        ret_floor = solver_nc(sn_floor, options)
        ret_ceil = solver_nc(sn_ceil, options)

        result.runtime = ret_floor.runtime + ret_ceil.runtime

        # Interpolate results
        if ret_floor.Q is not None and ret_ceil.Q is not None:
            result.QN = ret_floor.Q + eta * (ret_ceil.Q - ret_floor.Q)
        if ret_floor.U is not None and ret_ceil.U is not None:
            result.UN = ret_floor.U + eta * (ret_ceil.U - ret_floor.U)
        if ret_floor.R is not None and ret_ceil.R is not None:
            result.RN = ret_floor.R + eta * (ret_ceil.R - ret_floor.R)
        if ret_floor.T is not None and ret_ceil.T is not None:
            result.TN = ret_floor.T + eta * (ret_ceil.T - ret_floor.T)
        if ret_floor.X is not None and ret_ceil.X is not None:
            result.XN = ret_floor.X + eta * (ret_ceil.X - ret_floor.X)

        # Interpolate lG
        result.lG = ret_floor.lG + np.sum(eta) * (ret_floor.lG - ret_ceil.lG)

        result.it = ret_floor.it + ret_ceil.it
        result.method = ret_ceil.method
    else:
        # Integer populations - direct computation
        ret = solver_nc(sn, options)
        result.QN = ret.Q.copy() if ret.Q is not None else None
        result.UN = ret.U.copy() if ret.U is not None else None
        result.RN = ret.R.copy() if ret.R is not None else None
        result.TN = ret.T.copy() if ret.T is not None else None
        result.XN = ret.X.copy() if ret.X is not None else None
        result.lG = ret.lG
        result.it = ret.it
        result.method = ret.method

    # Calculate cycle times using Little's Law: C(k) = N(k) / X(k)
    if result.XN is not None and njobs is not None:
        K = sn.nclasses
        result.CN = np.zeros((1, K))
        for k in range(K):
            njobs_flat = njobs.flatten()
            xn_flat = result.XN.flatten()
            if k < len(njobs_flat) and k < len(xn_flat) and xn_flat[k] > 0:
                result.CN[0, k] = njobs_flat[k] / xn_flat[k]

    result.runtime = time.time() - start_time
    result.solver = "NC"

    if result.lG is not None and np.isfinite(result.lG):
        from line_solver.api.io import console as _console
        _console.step('normalizing constant obtained: log G = %.6g', float(result.lG))
    return result


def solver_ncld_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    Load-dependent NC solver analyzer.

    Performs load-dependent NC analysis with interpolation for non-integer
    populations.

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Raises:
        RuntimeError: For unsupported configurations
    """
    start_time = time.time()

    if options is None:
        options = SolverOptions()

    # Check for multiserver with exact method and open classes
    nservers = sn.nservers
    if nservers is not None:
        nservers_finite = nservers.copy()
        nservers_finite = nservers_finite[np.isfinite(nservers_finite)]
        has_infinite_jobs = sn.njobs is not None and np.any(np.isinf(sn.njobs))
        if (len(nservers_finite) > 0 and np.max(nservers_finite) > 1 and
                has_infinite_jobs and options.method == 'exact'):
            raise RuntimeError(
                "NC solver cannot provide exact solutions for open or mixed queueing networks. "
                "Remove the 'exact' option."
            )

    # Create floor/ceiling copies for non-integer populations
    njobs = sn.njobs
    if njobs is None:
        njobs = np.zeros((1, sn.nclasses))

    njobs_floor = np.floor(njobs)
    njobs_ceil = np.ceil(njobs)
    eta = np.abs(njobs - njobs_floor)

    # Check for non-integer populations
    non_integer_job = np.any(eta > FINE_TOL)

    result = NCResult()

    if non_integer_job:
        if options.method == 'exact':
            raise RuntimeError(
                "NC load-dependent solver cannot provide exact solutions for fractional populations."
            )

        # Interpolate between floor and ceiling
        sn_floor = sn.copy()
        sn_ceil = sn.copy()
        sn_floor.njobs = njobs_floor
        sn_ceil.njobs = njobs_ceil

        ret_floor = solver_ncld(sn_floor, options)
        ret_ceil = solver_ncld(sn_ceil, options)

        # Interpolate results
        if ret_floor.Q is not None and ret_ceil.Q is not None:
            result.QN = ret_floor.Q + eta * (ret_ceil.Q - ret_floor.Q)
        if ret_floor.U is not None and ret_ceil.U is not None:
            result.UN = ret_floor.U + eta * (ret_ceil.U - ret_floor.U)
        if ret_floor.R is not None and ret_ceil.R is not None:
            result.RN = ret_floor.R + eta * (ret_ceil.R - ret_floor.R)
        if ret_floor.T is not None and ret_ceil.T is not None:
            result.TN = ret_floor.T + eta * (ret_ceil.T - ret_floor.T)
        if ret_floor.X is not None and ret_ceil.X is not None:
            result.XN = ret_floor.X + eta * (ret_ceil.X - ret_floor.X)

        # Interpolate lG
        result.lG = ret_floor.lG + np.sum(eta) * (ret_floor.lG - ret_ceil.lG)

        result.it = ret_floor.it + ret_ceil.it
        result.method = ret_ceil.method
    else:
        # Integer populations - direct computation
        ret = solver_ncld(sn, options)
        result.QN = ret.Q
        result.UN = ret.U
        result.RN = ret.R
        result.TN = ret.T
        result.XN = ret.X
        result.lG = ret.lG
        result.it = ret.it
        result.method = ret.method

    # Calculate cycle times using Little's Law: C(k) = N(k) / X(k)
    if result.XN is not None and njobs is not None:
        K = sn.nclasses
        result.CN = np.zeros((1, K))
        for k in range(K):
            njobs_flat = njobs.flatten()
            xn_flat = result.XN.flatten()
            if k < len(njobs_flat) and k < len(xn_flat) and xn_flat[k] > 0:
                result.CN[0, k] = njobs_flat[k] / xn_flat[k]

    result.runtime = time.time() - start_time
    result.solver = "NCLD"

    return result


def _lossn_region_constraints(sn: NetworkStruct, f: int, station_idx: int,
                              K: int) -> Tuple[np.ndarray, np.ndarray]:
    """
    Rows of the admission rule A n <= C for region f, in the order in which the
    simulation engines test them.

    Every row is a function of the per-class occupancy of the region only, which
    is what the FiniteCapacityRegion API can express, so no row can distinguish
    stations inside the region. Rows the region leaves unbounded are DROPPED
    rather than given a surrogate capacity: a 1e6 stand-in would silently turn
    an unconstrained dimension into a truncation at 1e6 and make the exact
    transform allocate a dimension of a million coefficients for a constraint
    that does not exist.

    Reference: solver_nc_lossn_analyzer.m, subfunction lossn_region_constraints.
    """
    region_matrix = np.asarray(sn.region[f], dtype=float)
    rows = []
    rhs = []

    # Global job cap: sum_r n_r <= globalMaxJobs
    if region_matrix.shape[1] > K:
        global_max = region_matrix[station_idx, K]
        if global_max >= 0:
            rows.append(np.ones(K))
            rhs.append(float(global_max))

    # Memory budget: sum_r classSize_r n_r <= globalMaxMemory. The class sizes
    # are the row, so this is the one row that can legitimately be fractional.
    maxmem_all = getattr(sn, 'regionmaxmem', None)
    if maxmem_all is not None and len(maxmem_all) > f and maxmem_all[f] is not None:
        memvec = np.asarray(maxmem_all[f], dtype=float).ravel()
        if station_idx < memvec.size and memvec[station_idx] >= 0:
            sz = np.ones(K)
            regionsz = getattr(sn, 'regionsz', None)
            if regionsz is not None and np.size(regionsz) > 0:
                sz = np.asarray(regionsz, dtype=float)[f].ravel()[:K]
            rows.append(sz)
            rhs.append(float(memvec[station_idx]))

    # Per-class job caps: n_r <= classMaxJobs_r, already folded with
    # classMaxMemory when the region was added.
    for r in range(K):
        if region_matrix[station_idx, r] >= 0:
            row = np.zeros(K)
            row[r] = 1.0
            rows.append(row)
            rhs.append(float(region_matrix[station_idx, r]))

    # Explicit linear constraints from FiniteCapacityRegion.setConstraint
    lincon = getattr(sn, 'regionlincon', None)
    if lincon is not None and len(lincon) > f and lincon[f] is not None:
        linA = np.atleast_2d(np.asarray(lincon[f][0], dtype=float))
        linB = np.asarray(lincon[f][1], dtype=float).ravel()
        for k in range(linA.shape[0]):
            rows.append(linA[k, :K].copy())
            rhs.append(float(linB[k]) if k < linB.size else 0.0)

    if not rows:
        raise RuntimeError(
            "solver_nc_lossn_analyzer: the finite capacity region declares no bounded "
            "constraint, so it admits every arrival and is not a loss network; give it "
            "a global job cap, a memory budget, a per-class cap or an explicit linear "
            "constraint")
    return np.vstack(rows), np.asarray(rhs, dtype=float)


def solver_nc_lossn_analyzer(
    sn: NetworkStruct,
    options: Optional[SolverOptions] = None
) -> NCResult:
    """
    NC solver analyzer for open loss networks with FCR.

    Analyzes open queueing networks with a single multiclass Delay node inside a
    Finite Capacity Region (FCR) with DROP policy.

    The FCR admission rule is A n <= C on the per-class occupancy vector n of the
    region, where the rows of A are assembled from every constraint the region
    declares: the global job cap, the memory budget weighted by the per-class
    sizes, the per-class job caps, and any explicit linear constraint set with
    FiniteCapacityRegion.setConstraint.

    Method selection (options.method):
        'rec'             - MDD-rec (lossn_rec): the same constant as the exact
                            sum over the admissible set, obtained by one
                            memoised walk of the decision diagram holding it.
                            Places no integrality demand on A or C.
        'exact' (default) - Manjunath-Sikdar transform (lossn_manjunath): the
                            normalization constant is obtained exactly as a
                            multidimensional contour integral evaluated by
                            residues. Requires integer A and C.
        'erlangfp'        - Erlang fixed-point (reduced-load) approximation.
        'mci'             - Monte Carlo importance-sampling summation
                            (Ross-Wang 1992): estimates the normalization
                            constant g(C) and class blocking with confidence
                            intervals (options.samples/seed).

    The default is the residue transform on an integral region and MDD-rec on a
    fractional one. It used to fall back to 'erlangfp' there, an approximation,
    because the residue argument counts whole units; MDD-rec needs only that the
    admissible set be finite and bounded per coordinate, which it still is, so
    the fractional case is now exact as well.

    Args:
        sn: Network structure with FCR configuration
        options: Solver options

    Returns:
        NCResult with all performance metrics

    Reference:
        MATLAB: solver_nc_lossn_analyzer.m
    """
    from ...lossn import lossn_erlangfp, lossn_mci, lossn_manjunath, lossn_rec

    start_time = time.time()

    if options is None:
        options = SolverOptions()

    K = sn.nclasses
    M = sn.nstations

    rates = sn.rates
    if rates is None:
        raise RuntimeError("Network structure has no rates defined")

    # 1. Locate the delay station inside the region
    region_matrix = np.asarray(sn.region[0], dtype=float)
    stations_in_fcr = []
    for i in range(M):
        has_class_constraint = np.any(region_matrix[i, :K] >= 0)
        has_global_constraint = (region_matrix[i, K] >= 0
                                 if region_matrix.shape[1] > K else False)
        if has_class_constraint or has_global_constraint:
            stations_in_fcr.append(i)
    if not stations_in_fcr:
        raise RuntimeError("No stations found in FCR")
    delay_idx = stations_in_fcr[0]

    # 2. Offered load per class. A route carries nu_r = arrival rate times mean
    # holding time INSIDE the region, i.e. the visit ratio at the delay divided by
    # its service rate. Passing the bare arrival rate would be correct only for
    # unit mean service times, and reports the wrong blocking otherwise.
    nu = np.zeros(K)
    lam = np.zeros(K)
    mu = np.zeros(K)
    for r in range(K):
        source_idx = int(sn.refstat[r]) if getattr(sn, 'refstat', None) is not None else 0
        lam[r] = rates[source_idx, r]
        V_r = 1.0
        chains = getattr(sn, 'chains', None)
        visits = getattr(sn, 'visits', None)
        if chains is not None and visits is not None:
            chain_ids = np.nonzero(np.asarray(chains, dtype=float)[:, r])[0]
            if chain_ids.size > 0:
                c = int(chain_ids[0])
                Vc = visits.get(c) if isinstance(visits, dict) else visits[c]
                if Vc is not None:
                    Vc = np.asarray(Vc, dtype=float)
                    # visits is indexed by STATEFUL node, refstat and delay_idx by station
                    s2sf = getattr(sn, 'stationToStateful', None)
                    ref_sf = int(s2sf[source_idx]) if s2sf is not None else source_idx
                    delay_sf = int(s2sf[delay_idx]) if s2sf is not None else delay_idx
                    vref = Vc[ref_sf, r]
                    if vref > 0:
                        V_r = Vc[delay_sf, r] / vref
        mu[r] = rates[delay_idx, r]
        if mu[r] > 0:
            nu[r] = lam[r] * V_r / mu[r]

    # 3. Assemble the admission constraints A n <= C_vec of the region
    A, C_vec = _lossn_region_constraints(sn, 0, delay_idx, K)

    # 4. Select method and solve. An explicit 'mci' wins, then an explicit
    # 'erlangfp', then the exact transform; only 'default' falls back, and only on
    # a fractional region. Answering an explicit 'exact' with the Erlang
    # approximation would report an approximation under the name of an exact
    # method.
    method_str = getattr(options, 'method', 'default') or 'default'
    tokens = str(method_str).lower().replace('/', '.').split('.')
    is_integral = (np.all(np.abs(A - np.round(A)) < 1e-9)
                   and np.all(np.abs(C_vec - np.round(C_vec)) < 1e-9))

    if 'mci' in tokens:
        chosen = 'mci'
    elif 'erlangfp' in tokens:
        chosen = 'erlangfp'
    elif 'rec' in tokens:
        chosen = 'rec'
    elif 'exact' in tokens or 'manjunath' in tokens or 'ms' in tokens:
        chosen = 'exact'
    else:
        # the residue argument counts whole units; MDD-rec does not
        chosen = 'exact' if is_integral else 'rec'

    lG = np.nan  # normalization constant (finite for exact and mci)
    if chosen == 'mci':
        samples = getattr(options, 'samples', None)
        if samples is None or np.isinf(samples):
            samples = 100000
        seed = getattr(options, 'seed', None)
        qlen, loss, lG, _ci, niter = lossn_mci(nu, A, C_vec,
                                               samples=int(samples), seed=seed)
        actual_method = 'lossn.mci'
    elif chosen == 'rec':
        qlen, loss, lG, niter = lossn_rec(nu, A, C_vec)
        actual_method = 'lossn.rec'
    elif chosen == 'exact':
        qlen, loss, lG, niter = lossn_manjunath(nu, A, C_vec)
        actual_method = 'lossn.exact'
    else:
        qlen, loss, _E, niter = lossn_erlangfp(nu, A, C_vec)
        actual_method = 'lossn.erlangfp'

    # 5. Convert to standard outputs. qlen is the CARRIED load E[n_r], so the
    # carried throughput follows from Little's law at the infinite server.
    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    T = np.zeros((M, K))
    R = np.zeros((M, K))
    Xc = np.zeros(K)

    for r in range(K):
        source_idx = int(sn.refstat[r]) if getattr(sn, 'refstat', None) is not None else 0
        Xc[r] = lam[r] * (1.0 - loss[r])       # carried (accepted) rate
        T[delay_idx, r] = Xc[r]
        # The source emits the accepted (post-drop) rate so that the
        # routing-based arrival-rate computation (sn_get_arvr_from_tput) yields a
        # non-zero ArvR at the delay, consistent with the flow-conserving
        # departure throughput and with the rate reported by SolverJMT.
        T[source_idx, r] = Xc[r]
        Q[delay_idx, r] = qlen[r]              # mean number in the region
        if mu[r] > 0:
            R[delay_idx, r] = 1.0 / mu[r]      # response time = service time (IS)
        U[delay_idx, r] = qlen[r]              # IS "utilization" is the busy servers

    result = NCResult()
    result.QN = Q
    result.UN = U
    result.RN = R
    result.TN = T
    result.XN = Xc.reshape(1, -1)  # system throughput per class
    result.CN = np.zeros((1, K))  # cycle time not applicable for open networks
    result.lG = lG
    result.it = niter
    result.method = actual_method
    result.runtime = time.time() - start_time
    result.solver = "NC"

    if result.lG is not None and np.isfinite(result.lG):
        from line_solver.api.io import console as _console
        _console.step('normalizing constant obtained: log G = %.6g', float(result.lG))
    return result


__all__ = [
    'NCResult',
    'NCResultProb',
    'solver_nc_analyzer',
    'solver_ncld_analyzer',
    'solver_nc_lossn_analyzer',
]


def solver_nc_spn_analyzer(model, sn, options: Optional[SolverOptions] = None) -> NCResult:
    """Stationary analysis of a PRODUCT-FORM stochastic Petri net by MDD-rec.

    This is the 'rec' method of SolverNC, and the first analytical route LINE
    offers for a Petri net -- CTMC solves the explicit generator, SSA and LDES
    simulate, FLD fluidises. Three functions do the work and each is the subject
    of its own reference:

        spn_pf      decides the product form and derives the per-place factors
                    g_l (Coleman-Henderson-Taylor complex balance)
        mdd_rec     G = sum_S prod_l g_l(s_l) in O(sum_l nodes_l * |S_l|) rather
                    than O(|S|)                       (Balsamo-Marin-Stojic)
        spn_metrics mean tokens, place and mode utilisation, and throughputs,
                    all from masked walks of the same diagram

    WHAT THIS REACHES THAT THE EXPLICIT GENERATOR DOES NOT. The diagram stores
    the reachable set, never the generator, so the cost is set by the number of
    diagram nodes and not by |S|. It also does not need the marking to be a
    conserved job population: a mode may consume two tokens and produce one, or
    consume one and produce two, which is the fork-join and batch case that the
    MDD-rec paper exists to serve.

    Returns
    -------
    NCResult with QN, UN, RN, TN per (Place station, class), CN and XN per
    class, and lG the log normalising constant. The spn_pf certificate is
    attached as the ``pf`` attribute, carrying the spn_metrics output under
    ``pf['metrics']``.

    UN FOLLOWS LINE, NOT THE PAPER. A Place is an INF station, and LINE reports
    U = Q at an infinite server, which is what SolverCTMC returns for the same
    net. The paper's place utilisation u(P_j) = 1 - P(m_j = 0) is a different
    quantity and is reported separately, as pf['metrics']['placeUtil'].
    """
    import time
    from ...spn import spn_metrics, spn_pf

    t0 = time.time()
    pfopt = {'verbose': int(getattr(options, 'verbose', 0) or 0) > 1}
    cfg = getattr(options, 'config', None)
    if cfg is not None:
        bound = cfg.get('spn_bound') if isinstance(cfg, dict) else getattr(cfg, 'spn_bound', None)
        if bound is not None:
            pfopt['bound'] = bound
    tol = getattr(options, 'tol', None)
    if tol:
        pfopt['tol'] = max(float(tol), 1e-12)

    pf = spn_pf(model, pfopt)
    met = spn_metrics(pf['mdds'], pf['g'], pf['info'])

    M = int(sn.nstations)
    R = int(sn.nclasses)
    QN = np.zeros((M, R)); UN = np.zeros((M, R))
    RN = np.zeros((M, R)); TN = np.zeros((M, R))

    node_to_station = np.ravel(np.asarray(sn.nodeToStation)).astype(int)
    for pp, node in enumerate(pf['info']['places']):
        ist = int(node_to_station[node])
        if ist < 0:
            continue
        for k in range(R):
            l = pp * R + k
            QN[ist, k] = met['tokens'][l]
            # INF station: LINE charges one server per resident token, so U = Q.
            # The paper's 1 - P(m = 0) is met['placeUtil'], on the certificate.
            UN[ist, k] = met['tokens'][l]
            TN[ist, k] = met['placeTput'][l]
            if TN[ist, k] > 0:
                RN[ist, k] = QN[ist, k] / TN[ist, k]     # Little's law at the place

    # System throughput at the reference station of each class, and the response
    # time Little's law then fixes. A net whose class population is not
    # conserved has no meaningful N/X, so CN stays zero there rather than
    # reporting a ratio against a moving population.
    XN = np.zeros(R); CN = np.zeros(R)
    refstat = np.ravel(np.asarray(sn.refstat)).astype(int)
    for k in range(R):
        ref = int(refstat[k]) if k < refstat.size else -1
        if 0 <= ref < M:
            XN[k] = TN[ref, k]
        Nk = QN[:, k].sum()
        if XN[k] > 0 and Nk > 0:
            CN[k] = Nk / XN[k]

    pf['metrics'] = met
    res = NCResult(QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN,
                   lG=float(np.log(met['G'])), STeff=None, it=1,
                   runtime=time.time() - t0, method='rec', solver='SolverNC')
    res.pf = pf
    return res
