"""
MAM Solver handler.

Native Python implementation of MAM (Matrix-Analytic Methods) solver handler
that analyzes queueing networks with phase-type distributions and Markovian
arrival processes.

Port from MATLAB solver_mam_basic.m
"""

import numpy as np
import numpy.matlib as ml
from dataclasses import dataclass, field
from typing import Optional, Dict, List, Tuple, Any
import time
import warnings

from ...sn import (
    NetworkStruct,
    SchedStrategy,
    NodeType,
    sn_is_open_model,
    sn_is_closed_model,
    sn_get_demands_chain,
)

try:
    from ....constants import ProcessType
except Exception:
    ProcessType = None


def _extract_ph_for_phm1(sn, station_idx, class_idx=0):
    """Return (alpha, T) for the PH/M/1 sigma-root from sn.proc[station_idx][class_idx].

    Handles three storage forms LINE uses for PH/MAP processes:
      - dict {'k', 'mu'}     : Erlang(k, mu) -> bidiagonal sub-generator
      - dict {'rate'}        : Exponential(rate)
      - (D0, D1) tuple/list  : MAP — alpha = map_pie(D0, D1), T = D0
    Returns (None, None) on failure.
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
                # HyperExp {'probs': p, 'rates': mu}: parallel exponential phases
                p = np.asarray(ph['probs'], dtype=float).flatten()
                mu_h = np.asarray(ph['rates'], dtype=float).flatten()
                return p, np.diag(-mu_h)
            return None, None
        if isinstance(ph, (list, tuple)) and len(ph) >= 2:
            D0 = np.asarray(ph[0])
            D1 = np.asarray(ph[1])
            if D0.ndim == 2 and D0.shape == D1.shape and D0.shape[0] == D0.shape[1]:
                from ...mam.map_analysis import map_pie as _map_pie
                pie = np.asarray(_map_pie(D0, D1)).flatten()
                return pie, D0
        return None, None
    except Exception:
        return None, None

# Try to import MMAPPH1FCFS for accurate MMAP/PH/1/FCFS analysis
try:
    from ...butools.queues import MMAPPH1FCFS
    HAS_MMAPPH1FCFS = True
except ImportError:
    HAS_MMAPPH1FCFS = False
    MMAPPH1FCFS = None

# MMAPPH1NPPR solves the MMAP[K]/PH[K]/1 non-preemptive priority queue used
# for HOL stations with non-identical class priorities
try:
    from ...butools.queues import MMAPPH1NPPR
    HAS_MMAPPH1NPPR = True
except ImportError:
    HAS_MMAPPH1NPPR = False
    MMAPPH1NPPR = None

# Import qbd_setupdelayoff for FunctionTask analysis
try:
    from ...mam.qbd import qbd_setupdelayoff
    HAS_QBD_SETUPDELAYOFF = True
except ImportError:
    HAS_QBD_SETUPDELAYOFF = False
    qbd_setupdelayoff = None

# Import ETAQA departure process constructors
try:
    from ...mam.qbd_depproc import qbd_depproc_etaqa, qbd_depproc_etaqa_ps
    HAS_QBD_DEPPROC = True
except ImportError:
    HAS_QBD_DEPPROC = False
    qbd_depproc_etaqa = None
    qbd_depproc_etaqa_ps = None

# Import MAP/MMAP utilities for departure process construction
try:
    from ...mam.map_analysis import map_mean, map_normalize, map_scale, map_pie, map_scv, map_idc
    from ...mam.mmap_ops import mmap_hide, mmap_lambda, mmap_compress
    HAS_MAP_UTILS = True
except ImportError:
    HAS_MAP_UTILS = False


@dataclass
class SolverMAMOptions:
    """Options for MAM solver."""
    method: str = 'default'
    tol: float = 1e-6
    verbose: bool = False
    iter_max: int = 100
    iter_tol: float = 1e-6
    space_max: int = 128
    merge: str = 'super'
    compress: str = 'mixture.order1'
    num_cdf_pts: int = 200  # Number of points for CDF computation
    etaqa_trunc: int = 8  # ETAQA truncation level for departure process
    fj_sync_q_len: int = 2  # mmap_max synchronization queue length (dec.source.mmap)


@dataclass
class SolverMAMReturn:
    """
    Result of MAM solver handler.

    Attributes:
        Q: Mean queue lengths (M x K)
        U: Utilizations (M x K)
        R: Response times (M x K)
        T: Throughputs (M x K)
        C: Cycle times (1 x K)
        X: System throughputs (1 x K)
        A: Arrival rates (M x K)
        W: Waiting times (M x K)
        runtime: Runtime in seconds
        method: Method used
        it: Number of iterations
    """
    Q: Optional[np.ndarray] = None
    U: Optional[np.ndarray] = None
    R: Optional[np.ndarray] = None
    T: Optional[np.ndarray] = None
    C: Optional[np.ndarray] = None
    X: Optional[np.ndarray] = None
    A: Optional[np.ndarray] = None
    W: Optional[np.ndarray] = None
    runtime: float = 0.0
    method: str = "default"
    it: int = 0


def _get_visits(sn: NetworkStruct) -> np.ndarray:
    """
    Compute station visit ratios from routing matrix.
    Equivalent to MATLAB: V = cellsum(sn.visits)

    NOTE: sn.visits[c] is indexed by STATEFUL NODE index (shape nstateful x K).
    We must convert station index to stateful index using stationToStateful.

    Args:
        sn: Network structure

    Returns:
        Visit ratios matrix (M x K)
    """
    M = sn.nstations
    K = sn.nclasses

    # Initialize visits matrix
    V = np.zeros((M, K))

    if sn.visits is not None and len(sn.visits) > 0:
        # Get station to stateful mapping
        stationToStateful = getattr(sn, 'stationToStateful', None)
        if stationToStateful is None or len(stationToStateful) == 0:
            stationToStateful = np.arange(M)
        else:
            stationToStateful = np.asarray(stationToStateful).flatten()

        visit_matrices = []
        if isinstance(sn.visits, dict):
            visit_matrices = list(sn.visits.values())
        elif isinstance(sn.visits, np.ndarray):
            if sn.visits.dtype == object:
                visit_matrices = [entry for entry in sn.visits.flat if entry is not None]
            else:
                visit_matrices = [sn.visits]
        elif isinstance(sn.visits, (list, tuple)):
            visit_matrices = list(sn.visits)
        else:
            visit_matrices = [sn.visits]

        # Sum visits across all chains
        for visit_matrix in visit_matrices:
            if visit_matrix is None:
                continue

            v = np.asarray(visit_matrix)
            if v.ndim == 1:
                v = v.reshape((-1, 1))

            # Extract station rows using stationToStateful mapping when available.
            for ist in range(M):
                row_idx = ist
                if ist < len(stationToStateful):
                    mapped_idx = int(stationToStateful[ist])
                    if 0 <= mapped_idx < v.shape[0]:
                        row_idx = mapped_idx

                if 0 <= row_idx < v.shape[0]:
                    cols = min(K, v.shape[1])
                    V[ist, :cols] += v[row_idx, :cols]
    else:
        # Default: equal visits
        V = np.ones((M, K))

    return V


def _get_service_times(sn: NetworkStruct) -> np.ndarray:
    """
    Extract service times from network structure.
    S = 1 / rates

    For classes without service at a station (rate=0), S=Inf is used
    to indicate no service, not S=0 which would mean immediate service.

    Args:
        sn: Network structure

    Returns:
        Service times matrix (M x K)
    """
    M = sn.nstations
    K = sn.nclasses

    if hasattr(sn, 'rates') and sn.rates is not None:
        rates = np.asarray(sn.rates)
        with np.errstate(divide='ignore', invalid='ignore'):
            # When rate=0, S=Inf (no service), not S=0 (immediate service)
            S = 1.0 / rates
            # Convert NaN to Inf (no service)
            S = np.where(np.isnan(S), np.inf, S)
        return S

    return np.ones((M, K))


def _get_nservers(sn: NetworkStruct) -> np.ndarray:
    """
    Get number of servers per station.

    Args:
        sn: Network structure

    Returns:
        Number of servers array (M,)
    """
    M = sn.nstations

    if hasattr(sn, 'nservers') and sn.nservers is not None:
        nservers = np.asarray(sn.nservers).flatten()
        if len(nservers) == M:
            return nservers

    return np.ones(M)


def _get_scheduling(sn: NetworkStruct, station_idx: int) -> SchedStrategy:
    """
    Get scheduling strategy for a station.

    Args:
        sn: Network structure
        station_idx: Station index

    Returns:
        Scheduling strategy
    """
    if hasattr(sn, 'sched') and sn.sched is not None:
        if isinstance(sn.sched, dict):
            return sn.sched.get(station_idx, SchedStrategy.FCFS)
        elif isinstance(sn.sched, (list, np.ndarray)) and station_idx < len(sn.sched):
            return sn.sched[station_idx]

    return SchedStrategy.FCFS


def _is_delay_station(sn: NetworkStruct, station_idx: int) -> bool:
    """Check if station is a delay (infinite server) station."""
    nservers = _get_nservers(sn)
    return np.isinf(nservers[station_idx])


def _is_source_station(sn: NetworkStruct, station_idx: int) -> bool:
    """Check if station is a source."""
    sched = _get_scheduling(sn, station_idx)
    return sched == SchedStrategy.EXT


def _is_ps_station(sn: NetworkStruct, station_idx: int) -> bool:
    """Check if station uses processor sharing."""
    sched = _get_scheduling(sn, station_idx)
    return sched == SchedStrategy.PS


def _mam_detect_mmck(sn: NetworkStruct, ist: int, K: int):
    """
    Port of MATLAB mam_detect_mmck: decide whether station ist matches the
    M/M/c/K assumptions on the service side.

    Returns (is_mmck, shared_mu). True only when every active class has
    exponential service at the station and they all share the same rate.
    Disabled classes (NaN/<=0 rate) are skipped. The caller is responsible for
    verifying that the aggregated arrival process is single-phase Poisson.
    """
    from ....constants import ProcessType
    procid = getattr(sn, 'procid', None)
    rates = getattr(sn, 'rates', None)
    if procid is None or rates is None:
        return False, float('nan')
    procid = np.asarray(procid, dtype=object)
    rates = np.asarray(rates, dtype=float)
    mu_vals = []
    for k in range(K):
        if ist >= procid.shape[0] or k >= procid.shape[1]:
            continue
        pid = procid[ist, k]
        rate = rates[ist, k] if (ist < rates.shape[0] and k < rates.shape[1]) else float('nan')
        if pid != ProcessType.EXP:
            if isinstance(rate, float) and np.isnan(rate):
                continue  # disabled class
            return False, float('nan')
        if np.isnan(rate) or rate <= 0:
            continue  # no inflow for this class
        mu_vals.append(float(rate))
    if not mu_vals:
        return False, float('nan')
    if max(mu_vals) - min(mu_vals) > 1e-9 * max(1.0, max(mu_vals)):
        return False, float('nan')
    return True, mu_vals[0]


def _is_me_or_rap_service(sn: NetworkStruct, ist: int, K: int) -> bool:
    """
    True if the service process of any class at this station is a
    matrix-exponential (ME) or rational-arrival-process (RAP) representation.

    Neither is phase-type, so MMAPPH1FCFS -- which reads the service only
    through the (pie, D0) phase-type pair -- returns the wrong numbers for
    them and the RAP/RAP/1 QBD must be used instead.
    """
    from ....constants import ProcessType
    procid = getattr(sn, 'procid', None)
    if procid is None:
        return False
    procid = np.asarray(procid, dtype=object)
    if ist >= procid.shape[0]:
        return False
    for k in range(min(K, procid.shape[1])):
        if procid[ist, k] in (ProcessType.ME, ProcessType.RAP):
            return True
    return False


def _is_renewal_map(D0: np.ndarray, D1: np.ndarray) -> bool:
    """
    True if the process (D0, D1) is renewal.

    A renewal process embedded as a MAP has D1 = t * alpha (rank one): the phase
    entered after an event does not depend on the phase the process was in at
    that event, so successive interevent times are independent and the process
    is fully described by its marginal.

    The test is on the algebraic structure alone, so it holds equally for a MAP,
    a matrix-exponential (ME) and a rational arrival process (RAP): an ME
    renewal stream satisfies it, a correlated MAP or RAP does not. Mirrors
    MATLAB mam_is_renewal_map.m and the JAR Mam_is_renewal_map.
    """
    from ...mam.map_analysis import map_pie
    ns = D0.shape[0]
    if ns == 1:
        return True
    t_exit = D1 @ np.ones(ns)
    alpha = np.asarray(map_pie(D0, D1)).flatten()
    return bool(np.linalg.norm(D1 - np.outer(t_exit, alpha), 'fro')
                < 1e-9 * max(1.0, np.linalg.norm(D1, 'fro')))


def _srcproc_is_renewal(sn: NetworkStruct, jst: int, k: int = 0) -> bool:
    """
    True if station jst's class-k process is renewal.

    sn.proc stores a process either as a compact descriptor of a renewal
    distribution ({'k','mu'} Erlang, {'rate'} exponential, {'probs','rates'}
    hyperexponential) or as a (D0, D1) matrix pair. A compact descriptor names a
    renewal distribution by construction; only the matrix pair can encode
    correlation, and there the rank-one test decides.

    Returns False when neither form is present, which is the conservative side:
    callers use this to decide whether a marginal-only closed form may be
    applied, and applying one to a process whose correlation structure could not
    be established is the failure this guard exists to prevent. Mirrors MATLAB
    mam_srcproc_is_renewal.m.
    """
    proc = getattr(sn, 'proc', None)
    if proc is None or jst < 0 or jst >= len(proc):
        return False
    entry = proc[jst]
    if entry is None or len(entry) <= k:
        return False
    ph = entry[k]
    if isinstance(ph, dict):
        return ('k' in ph and 'mu' in ph) or ('rate' in ph) \
            or ('probs' in ph and 'rates' in ph)
    if not isinstance(ph, (list, tuple)) or len(ph) < 2:
        return False
    D0 = np.asarray(ph[0], dtype=float)
    D1 = np.asarray(ph[1], dtype=float)
    if D0.ndim != 2 or D0.shape != D1.shape or D0.shape[0] != D0.shape[1]:
        return False
    return _is_renewal_map(D0, D1)


def _is_me_or_rap_arrival(sn: NetworkStruct, ist: int, K: int) -> bool:
    """
    True if any open chain feeding station ist originates from a station whose
    arrival process is a matrix-exponential (ME) or rational arrival process
    (RAP), so the aggregate (D0, D1) reaching the station is legitimately
    non-Markovian.
    """
    from ....constants import ProcessType
    procid = getattr(sn, 'procid', None)
    njobs = getattr(sn, 'njobs', None)
    inchain_map = getattr(sn, 'inchain', None)
    refstat = getattr(sn, 'refstat', None)
    if procid is None or inchain_map is None or refstat is None:
        return False
    procid = np.asarray(procid, dtype=object)
    refstat = np.asarray(refstat).flatten()
    njobs = np.asarray(njobs).flatten() if njobs is not None else None
    for c, inchain in inchain_map.items():
        idx = np.asarray(inchain).flatten().astype(int)
        if idx.size == 0:
            continue
        if njobs is not None and np.isfinite(np.sum(njobs[idx])):
            continue   # closed chain: a Poisson surrogate is used
        jst = int(refstat[idx[0]])
        if jst < 0 or jst >= procid.shape[0]:
            continue
        for k in idx:
            if k < procid.shape[1] and procid[jst, k] in (ProcessType.ME, ProcessType.RAP):
                return True
    return False


def _station_name(sn: NetworkStruct, ist: int) -> str:
    """
    Name of station ist, for use in user-facing messages. Names are reported
    rather than indices because station indexing is 1-based in MATLAB and
    0-based here, so an index makes the same message read differently in each
    codebase.
    """
    nodenames = getattr(sn, 'nodenames', None)
    station_to_node = getattr(sn, 'stationToNode', None)
    if nodenames is not None and station_to_node is not None:
        idx = np.asarray(station_to_node).flatten()
        if ist < len(idx) and 0 <= int(idx[ist]) < len(nodenames):
            return str(nodenames[int(idx[ist])])
    return f"#{ist}"


def _is_function_station(sn: NetworkStruct, station_idx: int) -> bool:
    """Check if station is a FunctionTask station (with setup/delayoff)."""
    if hasattr(sn, 'isfunction') and sn.isfunction is not None:
        isfunction = np.asarray(sn.isfunction).flatten()
        if station_idx < len(isfunction):
            return isfunction[station_idx] == 1
    return False


def _extract_rate_from_distribution(dist) -> Tuple[Optional[float], float]:
    """Extract rate and SCV from a distribution object or numeric value.

    The rate is the reciprocal of the mean and the SCV is read from the
    distribution, mirroring the map_lambda/map_scv pair MATLAB
    solver_mam_basic applies to the setup/delay-off process. A numeric value is
    read as a mean and assumed exponential.

    Returns:
        (rate, scv), or (None, 1.0) when no mean can be recovered
    """
    if dist is None:
        return None, 1.0

    if isinstance(dist, (int, float)):
        value = float(dist)
        return (1.0 / value, 1.0) if value > 0 else (None, 1.0)

    mean = None
    for getter in ('getMean', 'get_mean'):
        if hasattr(dist, getter):
            try:
                candidate = getattr(dist, getter)()
                if candidate is not None and float(candidate) > 0:
                    mean = float(candidate)
                    break
            except Exception:
                pass

    if mean is None:
        return None, 1.0

    scv = 1.0
    for getter in ('getSCV', 'get_scv'):
        if hasattr(dist, getter):
            try:
                candidate = getattr(dist, getter)()
                if candidate is not None and float(candidate) > 0:
                    scv = float(candidate)
                    break
            except Exception:
                pass

    return 1.0 / mean, scv


def _get_function_params(sn: NetworkStruct, station_idx: int) -> Tuple[Optional[float], Optional[float], Optional[float], Optional[float]]:
    """Get setup and delayoff parameters for a FunctionTask station.

    Returns:
        (alpharate, alphascv, betarate, betascv) or (None, None, None, None) if not available
    """
    if not hasattr(sn, 'nodeparam') or sn.nodeparam is None:
        return None, None, None, None

    if not isinstance(sn.nodeparam, dict):
        return None, None, None, None

    # FunctionTask setupTime (station-indexed, from SolverMAM._get_sn) vs Queue.setDelayOff (node-indexed nodeparam via _refresh_local_vars) reach here with different shapes; they coincide only when no non-station node precedes the queue.
    entry = sn.nodeparam.get(station_idx)
    if not (isinstance(entry, dict) and 'setupTime' in entry):
        node_idx = station_idx
        if getattr(sn, 'stationToNode', None) is not None:
            station_to_node = np.asarray(sn.stationToNode).flatten()
            if station_idx < len(station_to_node) and station_to_node[station_idx] >= 0:
                node_idx = int(station_to_node[station_idx])

        param = sn.nodeparam.get(node_idx)
        if not isinstance(param, dict):
            return None, None, None, None

        # Per-class entries; MATLAB reads the last class carrying a setup time.
        classes = [r for r in sorted(param)
                   if isinstance(param.get(r), dict) and 'setupTime' in param[r]]
        if not classes:
            return None, None, None, None
        entry = param[classes[-1]]

    alpharate, alphascv = _extract_rate_from_distribution(entry.get('setupTime'))
    betarate, betascv = _extract_rate_from_distribution(entry.get('delayoffTime'))

    return alpharate, alphascv, betarate, betascv


def _delayoff_cold_probability(betarate: float, betascv: float, nu: float) -> float:
    """Probability that the delay-off timer beats an Exp(nu) idle period.

    This is the per-instance memoryless race of the function-station
    cold-start model: p_cold = P(delayoff < Exp(nu)) = LST of the delay-off
    time at nu. For an exponential delay-off this is betarate/(betarate+nu);
    for a general SCV the delay-off is fitted as an acyclic PH (pie, T) and
    p_cold = pie * (nu I - T)^-1 * (-T 1). Mirrors MATLAB solver_mam_basic.
    """
    if betascv == 1.0:
        return betarate / (betarate + nu)
    Tb = None
    try:
        from line_solver.lib.thirdparty.butools.ph.baseph import AcyclicPHFromMeansAndSCVs
        ph = AcyclicPHFromMeansAndSCVs([1.0 / betarate], [betascv])
        pie_b = np.asarray(ph[0]).reshape(1, -1)
        Tb = np.asarray(ph[1], dtype=float)
    except Exception:
        Tb = None
    if Tb is None:
        # Erlang fallback, as in qbd_setupdelayoff
        nb = max(1, int(round(1.0 / betascv)))
        rate_b = nb * betarate
        Tb = np.diag([-rate_b] * nb).astype(float)
        for i in range(nb - 1):
            Tb[i, i + 1] = rate_b
        pie_b = np.zeros((1, nb))
        pie_b[0, 0] = 1.0
    nb = Tb.shape[0]
    exit_vec = -Tb @ np.ones((nb, 1))
    return float(pie_b @ np.linalg.solve(nu * np.eye(nb) - Tb, exit_vec))


def _build_ph_service_params(
    sn: NetworkStruct,
    ist: int,
    S: np.ndarray,
    K: int
) -> Tuple[Optional[List], Optional[List]]:
    """
    Build PH service parameters (pie and D0) for each class at a station.

    Args:
        sn: Network structure
        ist: Station index
        S: Service times matrix
        K: Number of classes

    Returns:
        Tuple of (pie_list, D0_list) or (None, None) if not applicable
    """
    # Check if we have PH service distributions
    if not hasattr(sn, 'proc') or sn.proc is None or len(sn.proc) <= ist:
        return None, None

    proc_ist = sn.proc[ist] if ist < len(sn.proc) else None
    if proc_ist is None or not isinstance(proc_ist, (list, dict)):
        return None, None

    pie_list = []
    D0_list = []

    for k in range(K):
        ph = proc_ist[k] if k < len(proc_ist) else None
        if ph is None:
            # Use exponential as fallback
            rate = 1.0 / S[ist, k] if S[ist, k] > 0 else 1.0
            pie_list.append(ml.matrix([[1.0]]))
            D0_list.append(ml.matrix([[-rate]]))
        elif isinstance(ph, dict):
            # Check for Erlang distribution {'k': phases, 'mu': rate_per_phase}
            if 'k' in ph and 'mu' in ph:
                # Convert Erlang(k, mu) to PH representation
                k_phases = int(ph['k'])
                mu = float(ph['mu'])
                # Initial vector: start in first phase
                alpha = np.zeros(k_phases)
                alpha[0] = 1.0
                # Generator: -mu on diagonal, mu on superdiagonal
                T = np.zeros((k_phases, k_phases))
                for i in range(k_phases):
                    T[i, i] = -mu
                    if i < k_phases - 1:
                        T[i, i + 1] = mu
                pie_list.append(ml.matrix(alpha.reshape(1, -1)))
                D0_list.append(ml.matrix(T))
            elif 'rate' in ph:
                # Simple rate-based service (exponential)
                rate = float(ph['rate'])
                pie_list.append(ml.matrix([[1.0]]))
                D0_list.append(ml.matrix([[-rate]]))
            elif 'probs' in ph and 'rates' in ph:
                # HyperExp {'probs': p, 'rates': mu}: parallel exponential phases
                p = np.asarray(ph['probs'], dtype=float).flatten()
                mu_h = np.asarray(ph['rates'], dtype=float).flatten()
                pie_list.append(ml.matrix(p.reshape(1, -1)))
                D0_list.append(ml.matrix(np.diag(-mu_h)))
            else:
                # Unknown dict format, use exponential fallback
                rate = 1.0 / S[ist, k] if S[ist, k] > 0 else 1.0
                pie_list.append(ml.matrix([[1.0]]))
                D0_list.append(ml.matrix([[-rate]]))
        elif isinstance(ph, (list, tuple)) and len(ph) >= 2:
            arr0 = np.asarray(ph[0])
            arr1 = np.asarray(ph[1])
            if arr0.ndim == 2 and arr1.ndim == 2 and arr0.shape == arr1.shape and arr0.shape[0] == arr0.shape[1]:
                # MAP (D0, D1): D0 is the PH sub-generator T; pie is the stationary
                # entry distribution computed from map_pie.
                T = arr0
                try:
                    pie = np.asarray(map_pie(arr0, arr1)).flatten()
                except Exception:
                    pie = np.zeros(arr0.shape[0])
                    pie[0] = 1.0
                pie_list.append(ml.matrix(pie.reshape(1, -1)))
                D0_list.append(ml.matrix(T))
            else:
                # Legacy PH representation [alpha, T]
                alpha = arr0.flatten()
                T = arr1
                pie_list.append(ml.matrix(alpha.reshape(1, -1)))
                D0_list.append(ml.matrix(T))
        else:
            # Fallback to exponential
            rate = 1.0 / S[ist, k] if S[ist, k] > 0 else 1.0
            pie_list.append(ml.matrix([[1.0]]))
            D0_list.append(ml.matrix([[-rate]]))

    return pie_list, D0_list


def _build_mmap_arrival(
    sn: NetworkStruct,
    ist: int,
    lambdas: np.ndarray,
    V: np.ndarray,
    K: int,
    C: int
) -> Tuple[Optional[List], List[float], float]:
    """
    Build MMAP arrival process for a station.

    Args:
        sn: Network structure
        ist: Station index
        lambdas: Arrival rates per chain
        V: Visit ratios
        K: Number of classes
        C: Number of chains

    Returns:
        Tuple of (D_list, class_rates, total_lambda)
    """
    # Aggregate arrival rate
    total_lambda = sum(lambdas[c] * V[ist, k]
                     for c in range(C) if c in sn.inchain
                     for k in sn.inchain[c].flatten().astype(int)
                     if V[ist, k] > 0)

    if total_lambda <= 1e-10:
        return None, [], 0.0

    class_rates = []
    for k in range(K):
        rate_k = sum(lambdas[c] * V[ist, k]
                    for c in range(C) if c in sn.inchain
                    if k in sn.inchain[c].flatten().astype(int))
        class_rates.append(rate_k)

    # arrival MMAP built from the source PH processes per mark, rescaled to the per-class arrival rate; falls back to Poisson when source processes are unextractable.
    D_list = _build_source_mmap(sn, class_rates, K)
    if D_list is None:
        # Poisson fallback: D0 diagonal rate, one marking matrix per class
        D_list = [ml.matrix([[-total_lambda]])]
        for k in range(K):
            D_list.append(ml.matrix([[class_rates[k]]]))

    return D_list, class_rates, total_lambda


def _build_source_mmap(sn, class_rates, K):
    """Builds the marked arrival MMAP [D0, D1_1..D1_K] by superposing the
    per-class renewal MAPs of the source station, each rescaled to the
    target class rate. Returns None when the source processes cannot be
    extracted (callers then fall back to a Poisson MMAP)."""
    try:
        from ...sn.network_struct import NodeType as _NT
        st2node = np.asarray(sn.get_station_indices()).astype(int)
        src_t = int(_NT.SOURCE)
        src = None
        for i in range(int(sn.nstations)):
            if int(sn.nodetype[st2node[i]]) == src_t:
                src = i
                break
        if src is None:
            return None
        streams = []  # (class index, D0, D1)
        for k in range(K):
            if class_rates[k] <= 1e-12:
                continue
            alpha = None
            T = None
            D1 = None
            ph = None
            if getattr(sn, 'proc', None) is not None and src < len(sn.proc):
                proc_st = sn.proc[src]
                if proc_st is not None and k < len(proc_st):
                    ph = proc_st[k]
            if isinstance(ph, (list, tuple)) and len(ph) >= 2:
                a0 = np.asarray(ph[0], dtype=float)
                a1 = np.asarray(ph[1], dtype=float)
                if a0.ndim <= 1 and a1.ndim == 2 and a1.shape[0] == a1.shape[1]:
                    # [alpha, T] phase-type storage form
                    alpha = a0.flatten()
                    T = a1
                elif (a0.ndim == 2 and a0.shape == a1.shape
                      and a0.shape[0] == a0.shape[1]):
                    # (D0, D1) MAP storage form: use D1 directly
                    T = a0
                    D1 = a1
                    from ...mam.map_analysis import map_pie as _map_pie
                    alpha = np.asarray(_map_pie(a0, a1)).flatten()
            if alpha is None or T is None:
                alpha, T = _extract_ph_for_phm1(sn, src, k)
                if alpha is None or T is None:
                    return None
                T = np.asarray(T, dtype=float)
                alpha = np.asarray(alpha, dtype=float).flatten()
            if D1 is None:
                exit_rates = -T @ np.ones(T.shape[0])
                D1 = np.outer(exit_rates, alpha)
            # rescale the stream to the target class rate (uniform time
            # scaling preserves the scv)
            mean_ia = float(alpha @ np.linalg.solve(-T, np.ones(T.shape[0])))
            theta = class_rates[k] * mean_ia
            if not np.isfinite(theta) or theta <= 0:
                return None
            streams.append((k, T * theta, D1 * theta))
        if not streams:
            return None
        # superpose the independent class streams by Kronecker sums
        cur_D0 = streams[0][1]
        cur_marks = [np.zeros_like(cur_D0) for _ in range(K)]
        cur_marks[streams[0][0]] = streams[0][2]
        for (k, D0b, D1b) in streams[1:]:
            nA = cur_D0.shape[0]
            nB = D0b.shape[0]
            IA = np.eye(nA)
            IB = np.eye(nB)
            new_D0 = np.kron(cur_D0, IB) + np.kron(IA, D0b)
            new_marks = [np.kron(cur_marks[r], IB) for r in range(K)]
            new_marks[k] = new_marks[k] + np.kron(IA, D1b)
            cur_D0 = new_D0
            cur_marks = new_marks
        return [ml.matrix(cur_D0)] + [ml.matrix(mk) for mk in cur_marks]
    except Exception:
        return None


def _solve_fcfs_mmapph1_open(
    sn: NetworkStruct,
    ist: int,
    lambdas: np.ndarray,
    V: np.ndarray,
    S: np.ndarray,
    K: int,
    C: int,
    me_warned: Optional[set] = None
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Solve MMAP/PH/1/FCFS queue for OPEN networks using 'ncMoms'.

    Args:
        sn: Network structure
        ist: Station index
        lambdas: Arrival rates per chain
        V: Visit ratios
        S: Service times
        K: Number of classes
        C: Number of chains

    Returns:
        Tuple of (QN, RN) arrays for the station, or (None, None) if not applicable
    """
    if not HAS_MMAPPH1FCFS:
        return None, None

    nservers = _get_nservers(sn)
    c_serv = nservers[ist]
    if not np.isfinite(c_serv) or c_serv < 1:
        return None, None

    # Build PH service parameters
    pie_list, D0_list = _build_ph_service_params(sn, ist, S, K)
    if pie_list is None:
        return None, None

    # multi-server FCFS rescales PH to mean S/nservers (map_scale); mirrors MATLAB solver_mam_basic.m:59.
    if c_serv > 1:
        D0_list = [ml.matrix(np.asarray(D0) * c_serv) for D0 in D0_list]

    # Build MMAP arrival process
    D_list, class_rates, total_lambda = _build_mmap_arrival(sn, ist, lambdas, V, K, C)
    if D_list is None:
        return None, None

    # ME/RAP service is not phase-type, so the (pie,D0) pair MMAPPH1FCFS reads cannot describe it; RAP/RAP/1 QBD applies whenever service is ME/RAP regardless of arrival type. See _kb/06-solver-catalog.md MAM Exact closed-form fast paths section.
    if _is_me_or_rap_service(sn, ist, K):
        if K == 1 and c_serv == 1:
            svc = sn.proc[ist][0] if (getattr(sn, 'proc', None) is not None
                                      and ist < len(sn.proc) and len(sn.proc[ist]) > 0) else None
            if not (isinstance(svc, (list, tuple)) and len(svc) >= 2):
                raise RuntimeError(
                    f"RAP/RAP/1 analysis failed at station {ist}: the station is declared "
                    f"to have a matrix-exponential or rational service process but carries "
                    f"no (D0, D1) representation.")
            from ...mam.qbd import qbd_raprap1
            D0s = np.asarray(svc[0], dtype=float)
            D1s = np.asarray(svc[1], dtype=float)
            # ME rescaled by rate alone, not map_scale (which normalizes and would clip legitimate negative entries, replacing it with a different Markovian process); see _kb/06-solver-catalog.md MAM ME/RAP arrival and service processes section.
            svc_mean = map_mean(D0s, D1s)
            target_mean = float(S[ist, 0]) / c_serv
            if target_mean > 0 and svc_mean > 0:
                ratio = svc_mean / target_mean
                D0s = D0s * ratio
                D1s = D1s * ratio
            C0 = np.asarray(D_list[0], dtype=float)
            C1 = np.asarray(D_list[1], dtype=float)
            try:
                rap_out = qbd_raprap1((C0, C1), (D0s, D1s))
            except Exception as e:
                raise RuntimeError(
                    f"RAP/RAP/1 analysis failed at station {ist}, which has a "
                    f"matrix-exponential or rational service process: {e}") from e
            # qbd_raprap1 gives only the queue-length mean; response time is recovered via Little's law, as on the PH path.
            QN = np.zeros(K)
            RN = np.zeros(K)
            QN[0] = float(rap_out[1])
            for c in range(C):
                if c not in sn.inchain:
                    continue
                for k in sn.inchain[c].flatten().astype(int):
                    TN_k = lambdas[c] * V[ist, k]
                    if TN_k > 1e-10:
                        RN[k] = QN[k] / TN_k
            return QN, RN
        else:
            # qbd_raprap1 is single-class single-server only; warns via line_warning_always (not warnings.warn/line_warning) so every affected model in the solve is flagged exactly once via me_warned. See _kb/06-solver-catalog.md MAM Exact closed-form fast paths section.
            if me_warned is None or ist not in me_warned:
                if me_warned is not None:
                    me_warned.add(ist)
                from ...io.logging import line_warning_always
                line_warning_always(
                    "solver_mam_basic",
                    "Station %s has a matrix-exponential or rational service process, "
                    "which the RAP/RAP/1 analysis supports only with a single class at a "
                    "single server (here %d classes, %g servers). Falling back to the "
                    "phase-type approximation MMAPPH1FCFS, which is not exact for this "
                    "service process.",
                    _station_name(sn, ist), K, c_serv)

    # a genuine correlated MAP service uses the exact MAP/MAP/1 (q_ct_map_map_1), not the renewal-marginal MMAPPH1FCFS, since the renewal approximation discards service autocorrelation and understates queue length by an order of magnitude.
    if K == 1 and c_serv == 1:
        try:
            from ....constants import GlobalConstants
            from ...mam.map_analysis import map_acf
            from ....lib.thirdparty.qmam import q_ct_map_map_1, MAPMAP1Options
            svc = sn.proc[ist][0] if (getattr(sn, 'proc', None) is not None
                                      and ist < len(sn.proc) and len(sn.proc[ist]) > 0) else None
            if isinstance(svc, (list, tuple)) and len(svc) >= 2:
                D0s = np.asarray(svc[0], dtype=float)
                D1s = np.asarray(svc[1], dtype=float)
                if (D0s.ndim == 2 and D0s.shape == D1s.shape and D0s.shape[0] > 1
                        and abs(float(np.asarray(map_acf(D0s, D1s, 1)).flatten()[0])) > GlobalConstants.CoarseTol):
                    C0 = np.asarray(D_list[0], dtype=float)
                    C1 = np.asarray(D_list[1], dtype=float)
                    res_mm = q_ct_map_map_1(C0, C1, D0s, D1s, MAPMAP1Options(max_num_comp=100000))
                    ql = np.asarray(res_mm.queue_length, dtype=float).flatten()
                    en = float(np.sum(np.arange(len(ql)) * ql))
                    QN = np.zeros(K)
                    RN = np.zeros(K)
                    QN[0] = en
                    for c in range(C):
                        if c not in sn.inchain:
                            continue
                        inchain = sn.inchain[c].flatten().astype(int)
                        for k in inchain:
                            TN_k = lambdas[c] * V[ist, k]
                            if TN_k > 1e-10:
                                RN[k] = QN[k] / TN_k
                    return QN, RN
        except Exception as e:
            warnings.warn(f"MAP/MAP/1 exact route failed, falling back to MMAPPH1FCFS: {e}")

    # HOL with non-identical priorities solves MMAP[K]/PH[K]/1 (BUTools D1=lowest priority vs LINE's lower-value-is-higher convention).
    classprio = None
    if getattr(sn, 'classprio', None) is not None:
        classprio = np.asarray(sn.classprio).flatten()
    sched_ist = _get_scheduling(sn, ist)
    if (sched_ist == SchedStrategy.HOL and classprio is not None and K > 1
            and np.any(classprio != classprio[0]) and HAS_MMAPPH1NPPR):
        if len(np.unique(classprio)) != K:
            raise RuntimeError('Solver MAM requires either identical '
                               'priorities or all distinct priorities')
        iK = np.argsort(-classprio)  # lowest priority first
        D_pr = [D_list[0]] + [D_list[1 + int(k)] for k in iK]
        pie_pr = [pie_list[int(k)] for k in iK]
        S_pr = [D0_list[int(k)] for k in iK]
        result = MMAPPH1NPPR(D_pr, pie_pr, S_pr, 'ncMoms', 1)
        QN = np.zeros(K)
        RN = np.zeros(K)
        res_list = result if isinstance(result, list) else [result]
        for j in range(min(K, len(res_list))):
            q_j = res_list[j]
            if hasattr(q_j, '__iter__'):
                q_j = float(q_j[0]) if len(q_j) > 0 else 0.0
            QN[int(iK[j])] = float(q_j) if q_j is not None else 0.0
        for c in range(C):
            if c not in sn.inchain:
                continue
            for k in sn.inchain[c].flatten().astype(int):
                TN_k = lambdas[c] * V[ist, k]
                if TN_k > 1e-10:
                    RN[k] = QN[k] / TN_k
        return QN, RN

    try:
        # BuTools' CheckMMAPRepresentation rejects a well-formed ME/RAP arrival (legitimate negative D0 off-diagonal); suppressed only for ME/RAP calls, mirroring MATLAB's global BUT_CHECK_INPUT=false.
        from ....lib.thirdparty import butools
        arrival_is_rational = _is_me_or_rap_arrival(sn, ist, K)
        saved_check_input = butools.checkInput
        if arrival_is_rational:
            butools.checkInput = False
        try:
            # Call MMAPPH1FCFS with 'ncMoms' for open networks
            result = MMAPPH1FCFS(D_list, pie_list, D0_list, 'ncMoms', 1)
        finally:
            butools.checkInput = saved_check_input

        # Parse result - should be list of queue lengths per class
        if result is not None:
            QN = np.zeros(K)
            RN = np.zeros(K)

            if isinstance(result, list):
                for k in range(min(len(result), K)):
                    q_k = result[k]
                    if hasattr(q_k, '__iter__'):
                        q_k = float(q_k[0]) if len(q_k) > 0 else 0.0
                    QN[k] = float(q_k) if q_k is not None else 0.0
            else:
                QN[0] = float(result) if result is not None else 0.0

            # Compute response times from Little's Law: R = Q / lambda
            for c in range(C):
                if c not in sn.inchain:
                    continue
                inchain = sn.inchain[c].flatten().astype(int)
                for k in inchain:
                    TN_k = lambdas[c] * V[ist, k]
                    if TN_k > 1e-10:
                        RN[k] = QN[k] / TN_k

            return QN, RN

    except Exception as e:
        warnings.warn(f"MMAPPH1FCFS (open) failed: {e}")
        return None, None

    return None, None


def _solve_fcfs_mmapph1_closed(
    sn: NetworkStruct,
    ist: int,
    lambdas: np.ndarray,
    V: np.ndarray,
    S: np.ndarray,
    N: np.ndarray,
    K: int,
    C: int
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Solve MMAP/PH/1/FCFS queue for CLOSED networks using 'ncDistr'.

    Uses queue length distribution to compute per-class queue lengths.
    Port of MATLAB solver_mam_basic.m lines 289-297.

    IMPORTANT: MATLAB uses a SINGLE aggregate queue length distribution for ALL classes,
    truncating it per-class based on N(k). The adjustment pdistr_k(end) = abs(1-sum(pdistr(1:end-1)))
    uses the FULL original pdistr, not the truncated version.

    Args:
        sn: Network structure
        ist: Station index
        lambdas: Arrival rates per chain
        V: Visit ratios
        S: Service times
        N: Population per class
        K: Number of classes
        C: Number of chains

    Returns:
        Tuple of (QN, RN) arrays for the station, or (None, None) if not applicable
    """
    if not HAS_MMAPPH1FCFS:
        return None, None

    nservers = _get_nservers(sn)
    c_serv = nservers[ist]
    if not np.isfinite(c_serv) or c_serv < 1:
        return None, None

    # Build PH service parameters
    pie_list, D0_list = _build_ph_service_params(sn, ist, S, K)
    if pie_list is None:
        return None, None

    # multi-server FCFS solves the rate-scaled single-server equivalent plus a surrogate tandem delay S*(nservers-1)/nservers; exact for single-class, approximate for multiserver+multiclass.
    if c_serv > 1:
        if K > 1:
            return None, None
        D0_list = [ml.matrix(np.asarray(D0) * c_serv) for D0 in D0_list]

    # Build MMAP arrival process
    D_list, class_rates, total_lambda = _build_mmap_arrival(sn, ist, lambdas, V, K, C)
    if D_list is None:
        return None, None

    try:
        # For closed networks, maxLevel = sum of finite populations + 1
        # MATLAB: maxLevel = sum(N(isfinite(N)))+1
        finite_N = N[np.isfinite(N)]
        maxLevel = int(np.sum(finite_N)) + 1

        # MMAPPH1FCFS 'ncDistr' returns one distribution PER CLASS; reusing class 1's distribution for all classes (a prior MATLAB bug, also fixed there) drove every chain to the same lambda. See _kb/06-solver-catalog.md MAM Exact closed-form fast paths section.
        result = MMAPPH1FCFS(D_list, pie_list, D0_list, 'ncDistr', maxLevel)

        if result is not None:
            QN = np.zeros(K)
            RN = np.zeros(K)

            pdistr_per_class = None
            if isinstance(result, list) and len(result) > 0:
                pdistr_per_class = [np.asarray(r).flatten() for r in result]
            elif isinstance(result, np.ndarray):
                # Single distribution returned (K == 1)
                pdistr_per_class = [np.asarray(result).flatten()]
            else:
                return None, None

            # per-class queue-length distribution truncated and renormalized at that class's own population N[k]; mirrors MATLAB solver_mam_basic.m.

            for k in range(K):
                if not np.isfinite(N[k]):
                    continue

                # Each class reads its own distribution; fall back to the last one
                # available if BUTools returned fewer than K (e.g. K == 1).
                pdistr_full = pdistr_per_class[k] if k < len(pdistr_per_class) else pdistr_per_class[-1]

                Nk = int(N[k])
                # Get first N(k)+1 elements
                num_levels = min(Nk + 1, len(pdistr_full))
                pdistr_k = np.abs(pdistr_full[:num_levels].copy())

                # truncated-tail complement taken over the TRUNCATED vector (not the full distribution), or the levels between P(Q>N[k]) and P(Q>=N[k]) are lost; mirrors solver_mam_basic.m.
                if num_levels > 0:
                    pdistr_k[-1] = np.abs(1.0 - np.sum(pdistr_k[:-1]))

                # Normalize: pdistr_k = pdistr_k / sum(pdistr_k)
                pdistr_sum = np.sum(pdistr_k)
                if pdistr_sum > 1e-10:
                    pdistr_k = pdistr_k / pdistr_sum

                # Compute expected value: E[Q_k] = sum(i * p(i)) for i=0..Nk
                levels = np.arange(num_levels)
                expected_q = np.sum(levels * pdistr_k)

                # MATLAB: max(0, min(N(k), expected_q))
                QN[k] = max(0.0, min(float(Nk), expected_q))

            # Compute response times from Little's Law: R = Q / lambda
            for c in range(C):
                if c not in sn.inchain:
                    continue
                inchain = sn.inchain[c].flatten().astype(int)
                for k in inchain:
                    TN_k = lambdas[c] * V[ist, k]
                    if TN_k > 1e-10:
                        RN[k] = QN[k] / TN_k

            return QN, RN

    except Exception as e:
        warnings.warn(f"MMAPPH1FCFS (closed) failed: {e}")
        return None, None

    return None, None


def _solve_fcfs_mmapph1(
    sn: NetworkStruct,
    ist: int,
    lambdas: np.ndarray,
    V: np.ndarray,
    S: np.ndarray,
    K: int,
    C: int,
    me_warned: Optional[set] = None
) -> Tuple[Optional[np.ndarray], Optional[np.ndarray]]:
    """
    Solve MMAP/PH/1/FCFS queue using matrix-analytic method (for open networks).

    Args:
        sn: Network structure
        ist: Station index
        lambdas: Arrival rates per chain
        V: Visit ratios
        S: Service times
        K: Number of classes
        C: Number of chains

    Returns:
        Tuple of (QN, RN) arrays for the station, or (None, None) if not applicable
    """
    # Delegate to the open network solver
    return _solve_fcfs_mmapph1_open(sn, ist, lambdas, V, S, K, C, me_warned)


def solver_mam_basic(
    sn: NetworkStruct,
    options: Optional[SolverMAMOptions] = None
) -> SolverMAMReturn:
    """
    Basic MAM solver using source decomposition.

    Implements basic matrix-analytic decomposition for queueing networks.
    Port of MATLAB solver_mam_basic.m

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMAMReturn with performance metrics
    """
    start_time = time.time()

    # stations already flagged as non-phase-type service, scoped per solve so each is warned once.
    me_warned = set()

    if options is None:
        options = SolverMAMOptions()

    from ....constants import GlobalConstants

    tol = options.tol
    FINE_TOL = 1e-8

    M = sn.nstations
    K = sn.nclasses
    C = int(getattr(sn, 'nchains', 0) or K)
    N = np.asarray(sn.njobs).flatten() if sn.njobs is not None else np.array([])
    if N.size == 0:
        N = np.full(K, np.inf)
    elif N.size < K:
        N = np.pad(N, (0, K - N.size), constant_values=np.inf)
    sn.njobs = N

    inchain_map = getattr(sn, 'inchain', None)
    if not isinstance(inchain_map, dict) or len(inchain_map) == 0:
        inchain_map = {c: np.array([c], dtype=int) for c in range(K)}
    sn.inchain = inchain_map

    refstat = np.asarray(getattr(sn, 'refstat', np.array([]))).flatten()
    if refstat.size < K:
        refstat = np.pad(refstat, (0, K - refstat.size), constant_values=0)
    sn.refstat = refstat

    refclass = np.asarray(getattr(sn, 'refclass', np.array([]))).flatten()
    if refclass.size < C:
        pad_values = np.arange(refclass.size, C, dtype=int)
        refclass = np.concatenate([refclass, pad_values]) if refclass.size > 0 else np.arange(C, dtype=int)
    sn.refclass = refclass

    station_to_stateful = np.asarray(getattr(sn, 'stationToStateful', np.array([]))).flatten()
    if station_to_stateful.size < M:
        station_to_stateful = np.arange(M, dtype=int)
    sn.stationToStateful = station_to_stateful

    stateful_to_station = np.asarray(getattr(sn, 'statefulToStation', np.array([]))).flatten()
    if stateful_to_station.size < M:
        stateful_to_station = np.arange(M, dtype=int)
    sn.statefulToStation = stateful_to_station

    visits_obj = getattr(sn, 'visits', None)
    if isinstance(visits_obj, dict):
        normalized_visits = visits_obj
    elif isinstance(visits_obj, np.ndarray):
        if visits_obj.dtype == object:
            normalized_visits = {
                idx: visit for idx, visit in enumerate(visits_obj.flat) if visit is not None
            }
        else:
            normalized_visits = {0: visits_obj}
    elif isinstance(visits_obj, (list, tuple)):
        normalized_visits = {idx: visit for idx, visit in enumerate(visits_obj) if visit is not None}
    elif visits_obj is None:
        normalized_visits = {}
    else:
        normalized_visits = {0: visits_obj}
    sn.visits = normalized_visits

    # Get network parameters
    V = _get_visits(sn)
    S = _get_service_times(sn)
    nservers = _get_nservers(sn)

    # Get chain-level demands
    demands_result = sn_get_demands_chain(sn)
    Lchain = demands_result.Lchain

    Strue = S.copy()  # service times as declared, used to report utilization below

    # self-looping-class interference inflates OTHER classes' service time by (1+njobs_slc); see _kb/06-solver-catalog.md MAM Self-looping classes (SLC) section.
    isslc = np.zeros(K, dtype=bool)
    if getattr(sn, 'isslc', None) is not None:
        _isslc = np.asarray(sn.isslc).flatten()
        for k in range(min(K, _isslc.size)):
            isslc[k] = bool(_isslc[k])
    _refstat = np.asarray(getattr(sn, 'refstat', np.zeros(K))).flatten().astype(int)
    slcjobs = np.zeros(M)
    for k in range(K):
        if isslc[k]:
            ist_k = _refstat[k]
            if np.isfinite(nservers[ist_k]):
                slcjobs[ist_k] += N[k]
    for ist in range(M):
        if slcjobs[ist] > 0:
            S[ist, ~isslc] = S[ist, ~isslc] * (1 + slcjobs[ist])
            # The chain demands drive the throughput fixed point and must be
            # consistent with the inflated service times.
            Lchain[ist, :] = Lchain[ist, :] * (1 + slcjobs[ist])

    # Initialize result matrices
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    # Check model type
    is_open = sn_is_open_model(sn)
    is_closed = sn_is_closed_model(sn)
    is_mixed = is_open and is_closed

    # Mixed models are handled by treating each chain independently:
    # open chains use source arrival rates, closed chains iterate lambda via regula falsi

    # Initialize lambda (arrival rate per chain)
    lambdas = np.zeros(C)

    # an all-SLC chain has no surrogate arrival stream and is excluded from the fixed point, pinned by the SLC clamp instead; see _kb/06-solver-catalog.md MAM Self-looping classes (SLC) section.
    isslcchain = np.zeros(C, dtype=bool)
    for c in range(C):
        if c not in inchain_map:
            continue
        _ic = np.asarray(inchain_map[c]).flatten().astype(int)
        isslcchain[c] = _ic.size > 0 and bool(np.all(isslc[_ic]))

    # Get arrival rates for open chains and initial estimates for closed chains
    for c in range(C):
        if c not in inchain_map:
            continue
        inchain = np.asarray(inchain_map[c]).flatten().astype(int)

        # Check if chain is open
        is_open_chain = np.any(np.isinf(N[inchain]))

        if is_open_chain:
            # Open chain: prefer explicit external arrival rates when available.
            if hasattr(sn, 'lambda_arr') and sn.lambda_arr is not None:
                lambda_arr = np.asarray(sn.lambda_arr).flatten()
                finite_rates = lambda_arr[inchain] if lambda_arr.size > np.max(inchain) else np.array([])
                finite_rates = finite_rates[np.isfinite(finite_rates)]
                lambdas[c] = np.sum(finite_rates) if len(finite_rates) > 0 else 0.0
            elif hasattr(sn, 'rates') and sn.rates is not None:
                refstat_c = int(sn.refstat.flatten()[inchain[0]]) if sn.refstat is not None and len(np.asarray(sn.refstat).flatten()) > inchain[0] else 0
                rates_at_source = np.asarray(sn.rates)[refstat_c, inchain]
                finite_rates = rates_at_source[np.isfinite(rates_at_source)]
                lambdas[c] = np.sum(finite_rates) if len(finite_rates) > 0 else 0.0
        else:
            # Closed chain: initialize with lower bound (N_c / sum(Lchain[:, c]))
            Nc = np.sum(N[inchain])
            Lchain_sum = np.sum(Lchain[:, c])
            if Lchain_sum > FINE_TOL:
                lambdas[c] = Nc / Lchain_sum
            else:
                lambdas[c] = 1.0

        if isslcchain[c]:
            lambdas[c] = 0.0

    # For open chains, set throughputs at source
    for c in range(C):
        if c not in inchain_map:
            continue
        inchain = np.asarray(inchain_map[c]).flatten().astype(int)
        is_open_chain = np.any(np.isinf(N[inchain]))

        if is_open_chain:
            if hasattr(sn, 'lambda_arr') and sn.lambda_arr is not None:
                lambda_arr = np.asarray(sn.lambda_arr).flatten()
                if lambda_arr.size > np.max(inchain):
                    TN[0, inchain] = lambda_arr[inchain]
            elif hasattr(sn, 'rates') and sn.rates is not None:
                refstat_c = int(sn.refstat.flatten()[inchain[0]]) if sn.refstat is not None and len(np.asarray(sn.refstat).flatten()) > inchain[0] else 0
                TN[refstat_c, inchain] = np.asarray(sn.rates)[refstat_c, inchain]

    # Identify stations with finite servers (for utilization check)
    sd = np.isfinite(nservers)

    isclosedchain = np.zeros(C, dtype=bool)
    isopenchain = np.zeros(C, dtype=bool)
    for c in range(C):
        if c not in inchain_map:
            continue
        inchain_c = np.asarray(inchain_map[c]).flatten().astype(int)
        isopenchain[c] = bool(np.any(np.isinf(N[inchain_c])))
        isclosedchain[c] = (not isopenchain[c]) and (not isslcchain[c])
    # mixed-network backoff cannot use the uniform 1/Umax rule (would also scale open chains' exogenous arrival rate); see _kb/06-solver-catalog.md MAM Mixed-network throughput fixed point section.
    ismixed = bool(np.any(isclosedchain) and np.any(isopenchain))
    # closed chains in a mixed network are capped at the residual capacity open traffic leaves free, staying strictly inside the M/G/1-type stability region.
    Ulim = 1 - GlobalConstants.CoarseTol

    # Main iteration loop
    TN_prev = TN.copy() + np.inf
    it = 0

    while np.max(np.abs(TN - TN_prev)) > tol and it <= options.iter_max:
        it += 1
        TN_prev = TN.copy()

        # MATLAB's NaN-ignoring max freezes the closed-chain lambda update at TNlb when Umax is all-NaN, pinning LN class-switching host layers at the no-contention throughput.
        if np.any(sd):
            _urows = np.sum(UN[sd, :], axis=1)
            _ufinite = _urows[~np.isnan(_urows)]
            Umax = np.max(_ufinite) if _ufinite.size > 0 else np.nan
        else:
            Umax = 0.0
        if ismixed or Umax < 1.0:
            # Adjust lambda for closed chains based on queue lengths
            for c in range(C):
                if c not in inchain_map:
                    continue
                inchain = np.asarray(inchain_map[c]).flatten().astype(int)

                if isclosedchain[c]:
                    Nc = np.sum(N[inchain])
                    # MATLAB uses "omitnan" to ignore NaN values when summing
                    QNc = max(tol, np.nansum(QN[:, inchain]))
                    TNlb = Nc / max(FINE_TOL, np.sum(Lchain[:, c]))

                    if it == 1:
                        lambdas[c] = TNlb
                    else:
                        # Regula falsi iteration (iteration-averaged)
                        alpha = it / options.iter_max
                        lambdas[c] = lambdas[c] * alpha + (Nc / QNc) * lambdas[c] * (1 - alpha)
        if ismixed:
            # closed chains alone are backed onto the busiest finite-server station's residual capacity, not scaled by 1/Umax (which would also suppress open-chain throughput below its own source rate).
            Uchain = Lchain[sd, :] * lambdas
            Uopen = np.sum(Uchain[:, ~isclosedchain], axis=1)
            Uclosed = np.sum(Uchain[:, isclosedchain], axis=1)
            binding = Uclosed > tol
            if np.any(binding):
                theta = np.min((Ulim - Uopen[binding]) / Uclosed[binding])
                if theta < 1:
                    lambdas[isclosedchain] = lambdas[isclosedchain] * max(0.0, theta)
        elif Umax >= 1.0:
            lambdas = lambdas / Umax  # MATLAB: lambda = lambda * 1/Umax

        # Update throughputs: TN[m, k] = V[m, k] * lambda[c] for k in chain c
        for c in range(C):
            if c not in inchain_map:
                continue
            inchain = np.asarray(inchain_map[c]).flatten().astype(int)
            for m in range(M):
                TN[m, inchain] = V[m, inchain] * lambdas[c]

        # Compute metrics for each station
        for ist in range(M):
            if _is_source_station(sn, ist):
                # Source station: throughput equals arrival rate, no queue
                for c in range(C):
                    if c not in inchain_map:
                        continue
                    inchain = np.asarray(inchain_map[c]).flatten().astype(int)
                    is_open_chain = np.any(np.isinf(N[inchain]))
                    if is_open_chain and sn.rates is not None:
                        TN[ist, inchain] = sn.rates[ist, inchain]
                QN[ist, :] = 0.0
                UN[ist, :] = 0.0
                RN[ist, :] = 0.0
                continue

            if _is_delay_station(sn, ist):
                # Delay station (infinite server)
                # MATLAB: TN = lambda*V, UN = S*TN, QN = TN*S, RN = QN/TN
                for c in range(C):
                    if c not in inchain_map:
                        continue
                    inchain = np.asarray(inchain_map[c]).flatten().astype(int)
                    for k in inchain:
                        # Skip classes that don't visit this station
                        if V[ist, k] < FINE_TOL:
                            TN[ist, k] = 0.0
                            UN[ist, k] = 0.0
                            QN[ist, k] = 0.0
                            RN[ist, k] = 0.0
                            continue
                        # Skip classes without valid service (S = Inf or NaN)
                        if not np.isfinite(S[ist, k]):
                            TN[ist, k] = lambdas[c] * V[ist, k]
                            UN[ist, k] = 0.0
                            QN[ist, k] = 0.0
                            RN[ist, k] = 0.0
                            continue
                        TN[ist, k] = lambdas[c] * V[ist, k]
                        # INF stations report U=QLen=TN*S (no /c, since nservers=Inf); mirrors solver_mam_basic.m and every other solver's convention.
                        UN[ist, k] = S[ist, k] * TN[ist, k]
                        # QN=TN*S directly (TN already carries visits); a second V factor would undercount delay jobs on class-switching chains. Mirrors MATLAB solver_mam_basic.
                        QN[ist, k] = TN[ist, k] * S[ist, k]
                        RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else 0.0

            elif _is_ps_station(sn, ist):
                # PS station: TN=lambda*V, UN=S*TN, QN=UN/(1-Usum), always the capped formula (no saturated special case); second-pass rescaling corrects queue lengths.
                for c in range(C):
                    if c not in inchain_map:
                        continue
                    inchain = np.asarray(inchain_map[c]).flatten().astype(int)
                    for k in inchain:
                        # Skip classes that don't visit or have no valid service
                        if V[ist, k] < FINE_TOL or not np.isfinite(S[ist, k]):
                            TN[ist, k] = 0.0 if V[ist, k] < FINE_TOL else lambdas[c] * V[ist, k]
                            UN[ist, k] = 0.0
                            continue
                        TN[ist, k] = lambdas[c] * V[ist, k]
                        # Utilization Law: a c-server station holds TN*S/c of its capacity.
                        UN[ist, k] = S[ist, k] * TN[ist, k] / nservers[ist]

                # PS queue length via Uden capped below 1; always the capped formula even when saturated, corrected by the second-pass rescale.
                Uden = min(1.0 - FINE_TOL, np.sum(UN[ist, :]))
                for c in range(C):
                    if c not in inchain_map:
                        continue
                    inchain = np.asarray(inchain_map[c]).flatten().astype(int)
                    for k in inchain:
                        QN[ist, k] = UN[ist, k] / (1.0 - Uden)
                        RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else 0.0

            else:
                # FCFS or other queue
                # MATLAB: TN = lambda*V, UN = S*TN/nservers
                for c in range(C):
                    if c not in sn.inchain:
                        continue
                    inchain = sn.inchain[c].flatten().astype(int)
                    for k in inchain:
                        # a class with no service process here must leave UN=NaN (matching MATLAB's unguarded S*TN); the NaN is load-bearing (Umax gate, chain-sum void, response-time washout). See _kb/06-solver-catalog.md MAM NaN guards section.
                        if not np.isfinite(S[ist, k]):
                            TN[ist, k] = lambdas[c] * V[ist, k]
                            if _is_function_station(sn, ist):
                                # at a setup-bearing function station the NaN-freeze must NOT hold (the closed-chain regula falsi must keep moving to include the cold-start queue); see _kb/06-solver-catalog.md MAM NaN guards section.
                                UN[ist, k] = 0.0
                            else:
                                UN[ist, k] = np.nan
                            continue
                        if V[ist, k] < FINE_TOL:
                            TN[ist, k] = 0.0
                            UN[ist, k] = 0.0
                            continue
                        TN[ist, k] = lambdas[c] * V[ist, k]
                        UN[ist, k] = S[ist, k] * TN[ist, k] / nservers[ist]

                # FunctionTask setup/delayoff aggregate utilization uses 'omitnan' so undefined-service classes cannot poison the stability gate; mirrors MATLAB solver_mam_basic.m.
                aggr_util = np.nansum(UN[ist, :])
                qbd_setupdelayoff_success = False

                if _is_function_station(sn, ist) and not np.any(np.isinf(N)):
                    # Get setup and delayoff parameters
                    alpharate, alphascv, betarate, betascv = _get_function_params(sn, ist)
                    if alpharate is not None and betarate is not None:
                        # per-instance cold-start race for closed function layers; see _kb/06-solver-catalog.md MAM Setup/delay-off (function) stations section.
                        setup_mean = 1.0 / alpharate
                        inf_stations = np.isinf(nservers)
                        any_active = False
                        for c in range(C):
                            if c not in inchain_map:
                                continue
                            inchain = np.asarray(inchain_map[c]).flatten().astype(int)
                            # per-visit idle time = chain think demand / visits to this station per cycle (Lchain and V are per-chain-cycle quantities).
                            Vtot = float(np.sum(V[ist, inchain]))
                            ZT = float(np.nansum(Lchain[inf_stations, c])) / max(Vtot, FINE_TOL)
                            nu = 1.0 / max(ZT, FINE_TOL)
                            pcold = _delayoff_cold_probability(betarate, betascv, nu)
                            for k in inchain:
                                lam_k = TN[ist, k]
                                if np.isfinite(S[ist, k]) and lam_k > 0:
                                    RN[ist, k] = pcold * setup_mean + S[ist, k]
                                    QN[ist, k] = lam_k * RN[ist, k]
                                    any_active = True
                                else:
                                    # a class with zero visits or no service process must contribute 0, not NaN, to the chain sum QNc (mirrors MATLAB's Qret{k}=0 branch).
                                    QN[ist, k] = 0.0
                                    RN[ist, k] = 0.0
                        if any_active:
                            qbd_setupdelayoff_success = True

                elif _is_function_station(sn, ist) and qbd_setupdelayoff is not None:
                    # open setup/delay-off races the delay-off against the aggregate Poisson interarrival; must precede the exact-shortcut chain (qsys_mmck/qsys_phmc/qsys_mapdc/MMAPPH1FCFS), which models the station setup-free. See _kb/06-solver-catalog.md MAM Setup/delay-off section.
                    alpharate, alphascv, betarate, betascv = _get_function_params(sn, ist)
                    if alpharate is not None and betarate is not None:
                        # mu_k is the reciprocal of the service mean, so
                        # rho_k = lambda_k / mu_k = lambda_k * S.
                        rho_k = np.zeros(K)
                        for k in range(K):
                            if (np.isfinite(S[ist, k]) and S[ist, k] > FINE_TOL
                                    and np.isfinite(TN[ist, k]) and TN[ist, k] > 0):
                                rho_k[k] = TN[ist, k] * S[ist, k]
                        rho_total = float(np.sum(rho_k))
                        if rho_total > 0:
                            # aggregate service rate matches the aggregate load; the aggregate queue is split across classes in proportion to their load.
                            aggr_lambda_total = float(np.nansum(TN[ist, :]))
                            aggr_rate = aggr_lambda_total / rho_total
                            q_total = qbd_setupdelayoff(aggr_lambda_total, aggr_rate,
                                                        alpharate, alphascv,
                                                        betarate, betascv)
                            for k in range(K):
                                QN[ist, k] = q_total * rho_k[k] / rho_total
                                if TN[ist, k] > FINE_TOL:
                                    RN[ist, k] = QN[ist, k] / TN[ist, k]
                                else:
                                    RN[ist, k] = 0.0
                            qbd_setupdelayoff_success = True

                # Try to use MMAPPH1FCFS for accurate MMAP/PH/1/FCFS analysis
                mmapph1_success = False

                # finite-capacity FCFS/HOL with Poisson arrivals and shared exponential rate dispatches to exact M/M/c/K, ahead of the aggr_util<1 gate (a finite buffer stays finite even at rho>=1). See _kb/06-solver-catalog.md MAM Exact closed-form fast paths section.
                _capK = None
                if getattr(sn, 'cap', None) is not None:
                    _capvec = np.asarray(sn.cap).ravel()
                    if ist < _capvec.size and np.isfinite(_capvec[ist]):
                        _capK = int(_capvec[ist])
                _sched_ist = _get_scheduling(sn, ist)
                if (not qbd_setupdelayoff_success and _capK is not None
                        and _sched_ist in (SchedStrategy.FCFS, SchedStrategy.HOL)
                        and np.any(np.isinf(N))):
                    is_mmck, mu_mmck = _mam_detect_mmck(sn, ist, K)
                    # M/M/c/K also requires single-phase Poisson arrivals: verify
                    # the aggregated arrival MMAP has a 1x1 D0 (Poisson superposition).
                    if is_mmck:
                        _Dlist, _, _ = _build_mmap_arrival(sn, ist, lambdas, V, K, C)
                        if not (_Dlist is not None and len(_Dlist) >= 1
                                and np.asarray(_Dlist[0]).shape[0] == 1):
                            is_mmck = False
                    _c_serv = int(nservers[ist]) if np.isfinite(nservers[ist]) else 1
                    if is_mmck and mu_mmck > 0 and _capK >= _c_serv:
                        # Aggregate offered Poisson rate over active classes.
                        offered = np.zeros(K)
                        for k in range(K):
                            if np.isfinite(S[ist, k]) and S[ist, k] > FINE_TOL and np.isfinite(TN[ist, k]):
                                offered[k] = TN[ist, k]
                        aggrLambda = float(np.sum(offered))
                        if aggrLambda > FINE_TOL:
                            from ...qsys import qsys_mmck
                            res = qsys_mmck(aggrLambda, mu_mmck, _c_serv, _capK)
                            meanQ_fc = res['meanQueueLength']
                            loss_fc = res['lossProbability']
                            TN_eff = offered * (1.0 - loss_fc)
                            sumTN = float(np.sum(TN_eff))
                            # per-class service mean = 1/mu_mmck; the aggregate mean queue length is split by weighted-average sojourn.
                            Savg = 1.0 / mu_mmck
                            Wq = max(0.0, meanQ_fc / sumTN - Savg) if sumTN > 0 else 0.0
                            for k in range(K):
                                TN[ist, k] = TN_eff[k]
                                UN[ist, k] = TN_eff[k] * Savg / _c_serv
                                if TN_eff[k] > 0:
                                    RN[ist, k] = Wq + Savg
                                    QN[ist, k] = TN_eff[k] * RN[ist, k]
                                else:
                                    RN[ist, k] = 0.0
                                    QN[ist, k] = 0.0
                            mmapph1_success = True

                if not qbd_setupdelayoff_success and not mmapph1_success and aggr_util < 1.0 - FINE_TOL and np.any(np.isinf(N)):
                    mdc_skip_surrogate = False
                    # single-class open Det service dispatches to the exact Crommelin M/D/c embedded-DTMC solver (bypasses the surrogate-delay correction, already exact); mirrors MATLAB solver_mam_basic.m qsys_mapdc branch.
                    is_dmc = False
                    if (K == 1
                            and ProcessType is not None
                            and hasattr(sn, 'procid') and sn.procid is not None
                            and sn.procid.shape[0] > ist
                            and sn.procid[ist, 0] == ProcessType.EXP
                            and np.isfinite(nservers[ist]) and nservers[ist] >= 1):
                        for src_idx_dmc in (range(sn.procid.shape[0]) if hasattr(sn.procid, 'shape') else []):
                            if src_idx_dmc != ist and sn.procid[src_idx_dmc, 0] == ProcessType.DET:
                                is_dmc = True
                                src_idx_dmc_used = src_idx_dmc
                                break
                    if is_dmc:
                        try:
                            from ...qsys import qsys_dmc
                            mu_q = float(1.0 / S[ist, 0]) if S[ist, 0] > 0 else float('inf')
                            lam_d = float(np.asarray(sn.rates)[src_idx_dmc_used, 0])
                            res = qsys_dmc(lam_d, mu_q, int(nservers[ist]))
                            QN[ist, 0] = res['mean_queue_length']
                            if TN[ist, 0] > FINE_TOL:
                                RN[ist, 0] = QN[ist, 0] / TN[ist, 0]
                            mmapph1_success = True
                            mdc_skip_surrogate = True
                        except Exception as e:
                            warnings.warn(f"qsys_dmc failed, falling back: {e}")

                    is_mdc = (
                        not mmapph1_success
                        and K == 1
                        and ProcessType is not None
                        and hasattr(sn, 'procid') and sn.procid is not None
                        and sn.procid.shape[0] > ist
                        and sn.procid[ist, 0] == ProcessType.DET
                        and np.isfinite(nservers[ist]) and nservers[ist] >= 1
                    )
                    if is_mdc:
                        D_list, _, _ = _build_mmap_arrival(sn, ist, lambdas, V, K, C)
                        if D_list is not None and len(D_list) >= 2:
                            D0_arr = np.asarray(D_list[0])
                            D1_arr = np.asarray(D_list[1])
                            is_poisson = (D0_arr.shape == (1, 1)
                                          and D1_arr.shape == (1, 1)
                                          and D1_arr[0, 0] > 0
                                          and abs(D0_arr[0, 0] + D1_arr[0, 0]) < 1e-12)
                            if is_poisson:
                                try:
                                    from ...qsys import qsys_mdc_crommelin
                                    lam = float(D1_arr[0, 0])
                                    det_s = float(S[ist, 0])
                                    res = qsys_mdc_crommelin(lam, det_s, int(nservers[ist]))
                                    QN[ist, 0] = res['mean_queue_length']
                                    if TN[ist, 0] > FINE_TOL:
                                        RN[ist, 0] = QN[ist, 0] / TN[ist, 0]
                                    mmapph1_success = True
                                    mdc_skip_surrogate = True
                                except Exception as e:
                                    warnings.warn(f"qsys_mdc_crommelin failed, falling back: {e}")
                                    mdc_skip_surrogate = False
                            else:
                                mdc_skip_surrogate = False
                        else:
                            mdc_skip_surrogate = False
                    else:
                        mdc_skip_surrogate = False

                    # single-class open PH/M/1 (Exp service, c=1, PH source) uses the exact GI/M/1 sigma-root via qsys_phm1.
                    if (not mmapph1_success
                            and K == 1
                            and ProcessType is not None
                            and hasattr(sn, 'procid') and sn.procid is not None
                            and sn.procid.shape[0] > ist
                            and sn.procid[ist, 0] == ProcessType.EXP
                            and np.isfinite(nservers[ist]) and int(nservers[ist]) >= 1):
                        src_idx_phm1 = -1
                        for jst in range(sn.procid.shape[0]):
                            if jst == ist:
                                continue
                            src_idx_phm1 = jst
                            if sn.procid[jst, 0] != ProcessType.EXP:
                                break
                        # qsys_phmc requires a RENEWAL arrival (reads only the (pie,D0) marginal); correlated MAP/RAP/ME arrivals fall through to MMAPPH1FCFS instead. See _kb/06-solver-catalog.md MAM Exact closed-form fast paths section.
                        arrival_is_renewal = (src_idx_phm1 >= 0
                                              and _srcproc_is_renewal(sn, src_idx_phm1))
                        try:
                            from ...qsys import qsys_phmc
                            pie_p, D0p = (_extract_ph_for_phm1(sn, src_idx_phm1)
                                          if arrival_is_renewal else (None, None))
                            if pie_p is not None and D0p is not None:
                                mu_q = float(1.0 / S[ist, 0]) if S[ist, 0] > 0 else float('inf')
                                c_serv = int(nservers[ist])
                                res = qsys_phmc(pie_p, D0p, mu_q, c_serv)
                                QN[ist, 0] = res['mean_queue_length']
                                if TN[ist, 0] > FINE_TOL:
                                    RN[ist, 0] = QN[ist, 0] / TN[ist, 0]
                                mmapph1_success = True
                                mdc_skip_surrogate = True
                        except Exception as e:
                            warnings.warn(f"qsys_phmc failed, falling back: {e}")

                    if not mmapph1_success:
                        # Try MMAPPH1FCFS with 'ncMoms' for open classes
                        QN_mmap, RN_mmap = _solve_fcfs_mmapph1(sn, ist, lambdas, V, S, K, C, me_warned)
                        if QN_mmap is not None:
                            QN[ist, :] = QN_mmap
                            RN[ist, :] = RN_mmap
                            # multi-server surrogate delay correction QN += TN*S*(nservers-1)/nservers compensates the rate-scaled single-server MMAPPH1FCFS solve; mirrors MATLAB solver_mam_basic.m.
                            if np.isfinite(nservers[ist]) and nservers[ist] > 1 and not mdc_skip_surrogate:
                                for k in range(K):
                                    if np.isfinite(S[ist, k]) and TN[ist, k] > FINE_TOL:
                                        QN[ist, k] += TN[ist, k] * S[ist, k] * (nservers[ist] - 1) / nservers[ist]
                                        RN[ist, k] = QN[ist, k] / TN[ist, k]
                            mmapph1_success = True

                # Try MMAPPH1FCFS with 'ncDistr' for closed networks (all classes are closed)
                # MATLAB: lines 259-297 in solver_mam_basic.m
                if not qbd_setupdelayoff_success and not mmapph1_success and aggr_util < 1.0 - FINE_TOL:
                    if not np.any(np.isinf(N)):
                        # All classes are closed - use ncDistr for per-class queue lengths
                        QN_mmap, RN_mmap = _solve_fcfs_mmapph1_closed(sn, ist, lambdas, V, S, N, K, C)
                        if QN_mmap is not None:
                            QN[ist, :] = QN_mmap
                            RN[ist, :] = RN_mmap
                            # surrogate tandem delay correction applied for every non-exact FCFS path; mirrors MATLAB solver_mam_basic.m.
                            if np.isfinite(nservers[ist]) and nservers[ist] > 1:
                                for k in range(K):
                                    if np.isfinite(S[ist, k]) and TN[ist, k] > FINE_TOL:
                                        QN[ist, k] += TN[ist, k] * S[ist, k] * (nservers[ist] - 1) / nservers[ist]
                                        RN[ist, k] = QN[ist, k] / TN[ist, k]
                            mmapph1_success = True

                if not qbd_setupdelayoff_success and not mmapph1_success:
                    # Fallback: approximate queue lengths using M/M/c formula
                    if aggr_util < 1.0 - FINE_TOL:
                        for c in range(C):
                            if c not in sn.inchain:
                                continue
                            inchain = sn.inchain[c].flatten().astype(int)
                            for k in inchain:
                                # Simplified M/M/c approximation
                                QN[ist, k] = UN[ist, k] / (1.0 - aggr_util)
                                # Add surrogate delay for multiserver
                                QN[ist, k] += TN[ist, k] * S[ist, k] * (nservers[ist] - 1) / nservers[ist]
                                RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else 0.0
                    else:
                        # Saturated queue
                        for c in range(C):
                            if c not in sn.inchain:
                                continue
                            inchain = sn.inchain[c].flatten().astype(int)
                            for k in inchain:
                                QN[ist, k] = N[k] if np.isfinite(N[k]) else UN[ist, k]
                                RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else 0.0

                # FCFS path leaves Qret[k]=NaN for a class with no service process (triggers the post-loop population wash to RN=S); the function-station branch instead pins QN=0 to avoid voiding the chain sum. Mirrors MATLAB.
                if not qbd_setupdelayoff_success:
                    for k in range(K):
                        if not np.isfinite(S[ist, k]):
                            QN[ist, k] = np.nan

    totiter = it + 2

    # Compute cycle times
    CN = np.sum(RN, axis=0).reshape(1, -1)
    QN = np.abs(QN)

    # second-pass renormalization to match population; NaN must propagate through the scaling exactly as MATLAB/Java do (undefined-service cells stay NaN, not silently zeroed).
    for _ in range(2):
        for c in range(C):
            if c not in inchain_map:
                continue
            inchain = np.asarray(inchain_map[c]).flatten().astype(int)
            Nc = np.sum(N[inchain])

            if np.isfinite(Nc) and Nc > 0:
                # Compute QNc - check if any value is NaN
                QN_inchain = QN[:, inchain]
                has_nan = np.any(np.isnan(QN_inchain))

                if has_nan:
                    # NaN propagation: set QNc to NaN, then set all QN values to NaN
                    QNc = np.nan
                    # Set all finite QN values in this chain to NaN (matches MATLAB behavior)
                    mask = np.isfinite(QN[:, inchain])
                    QN[:, inchain] = np.where(mask, np.nan, QN[:, inchain])
                else:
                    QNc = np.sum(QN[:, inchain])
                    if np.isfinite(QNc) and QNc > FINE_TOL:
                        ratio = Nc / QNc
                        QN[:, inchain] = QN[:, inchain] * ratio

            # Recompute response times
            for ist in range(M):
                # Skip source stations - they have no queue (Q=0, R=0)
                if _is_source_station(sn, ist):
                    continue
                for k in inchain:
                    if V[ist, k] > 0:
                        if _is_delay_station(sn, ist):
                            RN[ist, k] = S[ist, k]
                        else:
                            # undefined service treated as NaN here (not inf) so a no-service cell ends NaN as in MATLAB, instead of RN=inf, QN=inf*0; MATLAB's max ignores NaN.
                            s_val = S[ist, k] if np.isfinite(S[ist, k]) else np.nan
                            qn_tn = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else np.nan
                            if np.isnan(qn_tn) and np.isnan(s_val):
                                RN[ist, k] = np.nan
                            elif np.isnan(qn_tn):
                                RN[ist, k] = s_val
                            elif np.isnan(s_val):
                                RN[ist, k] = qn_tn
                            else:
                                RN[ist, k] = max(s_val, qn_tn)
                    else:
                        RN[ist, k] = 0.0
                    QN[ist, k] = RN[ist, k] * TN[ist, k]

            # Handle zero population chains
            if Nc == 0:
                QN[:, inchain] = 0.0
                UN[:, inchain] = 0.0
                RN[:, inchain] = 0.0
                TN[:, inchain] = 0.0

    # Set system throughputs
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        XN[0, inchain] = lambdas[c]

    # SLC clamp applied last (after rescaling); U_slc=1-sum_j U_j from the TRUE uninflated service time. See _kb/06-solver-catalog.md MAM Self-looping classes (SLC) section.
    if np.any(isslc):
        refstat = np.asarray(getattr(sn, 'refstat', np.zeros(K))).flatten().astype(int)
        rates = np.asarray(sn.rates)
        # utilization follows the declared service time; the interference inflation applied above must not be reported as utilization or fed into the leftover-capacity identity.
        for ist in range(M):
            if slcjobs[ist] > 0:
                UN[ist, ~isslc] = Strue[ist, ~isslc] * TN[ist, ~isslc]
        QN[:, isslc] = 0.0
        UN[:, isslc] = 0.0
        RN[:, isslc] = 0.0
        TN[:, isslc] = 0.0
        for ist in range(M):
            slck = [k for k in range(K) if isslc[k] and refstat[k] == ist]
            if not slck:
                continue
            if np.isinf(nservers[ist]):
                # Delay station: no contention, every customer is always in
                # service, so the class completes at its full aggregate rate.
                for k in slck:
                    QN[ist, k] = N[k]
                    TN[ist, k] = N[k] * rates[ist, k]
                    RN[ist, k] = Strue[ist, k]
                    UN[ist, k] = Strue[ist, k] * TN[ist, k]
            else:
                nsrv = nservers[ist]
                Uleft = max(0.0, 1.0 - np.sum(UN[ist, ~isslc]))
                # Several self-looping classes at one station share the free
                # capacity in proportion to the service rate they offer.
                w = np.array([N[k] * rates[ist, k] for k in slck], dtype=float)
                if np.sum(w) <= 0:
                    continue
                for i, k in enumerate(slck):
                    ucap = min(N[k], nsrv) / nsrv
                    UN[ist, k] = min(Uleft * w[i] / np.sum(w), ucap)
                    TN[ist, k] = rates[ist, k] * UN[ist, k] * nsrv
                    QN[ist, k] = N[k]
                    RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > 0 else 0.0
        CN[0, :] = np.nansum(RN, axis=0)

    # Clean up NaN values
    QN = np.nan_to_num(QN, nan=0.0)
    UN = np.nan_to_num(UN, nan=0.0)
    RN = np.nan_to_num(RN, nan=0.0)
    TN = np.nan_to_num(TN, nan=0.0)
    CN = np.nan_to_num(CN, nan=0.0)
    XN = np.nan_to_num(XN, nan=0.0)

    result = SolverMAMReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.A = TN.copy()
    result.W = RN.copy()
    result.runtime = time.time() - start_time
    result.method = "dec.source"
    result.it = totiter

    return result


def solver_mam(
    sn: NetworkStruct,
    options: Optional[SolverMAMOptions] = None
) -> SolverMAMReturn:
    """
    Main MAM solver handler.

    Routes to appropriate method based on options and network characteristics.
    For 'dec.mmap', implements MMAP decomposition with ETAQA departure processes.

    Port from MATLAB solver_mam.m

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMAMReturn with performance metrics
    """
    if options is None:
        options = SolverMAMOptions()

    method = options.method.lower()

    if method in ['default', 'dec.source', 'dec.poisson']:
        # Use basic decomposition method
        if method == 'dec.poisson':
            options.space_max = 1
        return solver_mam_basic(sn, options)
    elif method == 'dec.mmap':
        return _solver_mam_dec_mmap(sn, options)
    elif method in ['mna', 'inap', 'exact']:
        # Matrix-analytic / RCAT methods - use basic for now
        return solver_mam_basic(sn, options)
    else:
        # Unknown method - use basic
        if options.verbose:
            print(f"Warning: Unknown MAM method '{method}'. Using dec.source.")
        return solver_mam_basic(sn, options)


def _solver_mam_dec_mmap(
    sn: NetworkStruct,
    options: SolverMAMOptions
) -> SolverMAMReturn:
    """
    MMAP decomposition solver (dec.mmap).

    Iteratively builds departure MAPs using ETAQA truncation and propagates
    them through the network via traffic equations. Supports FCFS, HOL, and PS
    scheduling strategies.

    Port from MATLAB solver_mam.m

    Args:
        sn: Network structure
        options: Solver options

    Returns:
        SolverMAMReturn with performance metrics
    """
    start_time = time.time()
    FINE_TOL = 1e-8

    config = options
    PH = sn.proc
    M = sn.nstations
    K = sn.nclasses
    C = int(getattr(sn, 'nchains', 0) or K)
    N = np.asarray(sn.njobs).flatten() if sn.njobs is not None else np.full(K, np.inf)
    V = _get_visits(sn)

    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    CN = np.zeros((1, K))
    XN = np.zeros((1, K))

    # Compute per-class arrival rates from source
    inchain_map = getattr(sn, 'inchain', None)
    if not isinstance(inchain_map, dict) or len(inchain_map) == 0:
        inchain_map = {c: np.array([c], dtype=int) for c in range(K)}
    sn.inchain = inchain_map

    lambda_k = np.zeros(K)
    for c in range(C):
        if c not in inchain_map:
            continue
        inchain = np.asarray(inchain_map[c]).flatten().astype(int)
        if hasattr(sn, 'rates') and sn.rates is not None:
            refstat = np.asarray(getattr(sn, 'refstat', np.array([]))).flatten()
            if refstat.size > inchain[0]:
                refstat_c = int(refstat[inchain[0]])
            else:
                refstat_c = 0
            lambdas_inchain = np.asarray(sn.rates)[refstat_c, inchain]
            lambdas_inchain = lambdas_inchain[np.isfinite(lambdas_inchain)]
            lambda_k[inchain] = np.sum(lambdas_inchain) if len(lambdas_inchain) > 0 else 0.0

    # Change A: Validate scheduling strategies (FCFS, HOL, PS supported)
    for ist in range(M):
        sched = _get_scheduling(sn, ist)
        if sched == SchedStrategy.EXT:
            pass  # source, OK
        elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.PS]:
            pass  # supported
        else:
            if options.verbose:
                warnings.warn(
                    f"The dec.mmap method does not support scheduling strategy "
                    f"at station {ist}."
                )
            result = SolverMAMReturn()
            result.runtime = time.time() - start_time
            result.method = ""
            result.it = 0
            return result

    if not all(np.isinf(N)):
        # Only open models supported by dec.mmap currently
        if options.verbose:
            warnings.warn("dec.mmap currently only supports open models.")
        result = SolverMAMReturn()
        result.runtime = time.time() - start_time
        result.method = ""
        result.it = 0
        return result

    S = _get_service_times(sn)
    nservers = _get_nservers(sn)

    # Build PH service distributions per station/class
    # PH[ist][k] = [D0, D1] as a MAP
    PH_map = {}
    pie = {}
    D0_ph = {}
    for ist in range(M):
        sched = _get_scheduling(sn, ist)
        if sched == SchedStrategy.EXT:
            TN[ist, :] = np.asarray(sn.rates)[ist, :] if sn.rates is not None else 0.0
            TN[ist, np.isnan(TN[ist, :])] = 0.0
            continue
        if sched in [SchedStrategy.FCFS, SchedStrategy.HOL, SchedStrategy.PS]:
            PH_map[ist] = {}
            pie[ist] = {}
            D0_ph[ist] = {}
            for k in range(K):
                ph_k = _get_ph_as_map(sn, ist, k, S)
                if ph_k is not None:
                    # Scale service time by number of servers
                    mean_s = map_mean(ph_k[0], ph_k[1])
                    if nservers[ist] > 1 and np.isfinite(mean_s) and mean_s > 0:
                        scale_factor = mean_s / nservers[ist]
                        ph_k[0], ph_k[1] = map_scale(ph_k[0], ph_k[1],
                                                      scale_factor / mean_s if mean_s > 0 else 1.0)
                    PH_map[ist][k] = ph_k
                    pie[ist][k] = map_pie(ph_k[0], ph_k[1])
                    D0_ph[ist, k] = ph_k[0]
                    if np.any(np.isnan(D0_ph[ist, k])):
                        # Immediate service fallback
                        imm_rate = 1e10
                        D0_ph[ist, k] = np.array([[-imm_rate]])
                        pie[ist][k] = np.array([1.0])
                        PH_map[ist][k] = [np.array([[-imm_rate]]),
                                          np.array([[imm_rate]])]
                else:
                    # Exponential fallback
                    rate = 1.0 / S[ist, k] if S[ist, k] > 0 and np.isfinite(S[ist, k]) else 1.0
                    PH_map[ist][k] = [np.array([[-rate]]), np.array([[rate]])]
                    pie[ist][k] = np.array([1.0])
                    D0_ph[ist, k] = np.array([[-rate]])

    # ETAQA truncation level
    etaqa_n = getattr(config, 'etaqa_trunc', 8)

    it_max = options.iter_max

    # Initialize DEP as scaled PH service distributions
    DEP = {}

    # Import traffic solver
    try:
        from ...npfqn.traffic import npfqn_traffic_split_cs, npfqn_traffic_merge
        from ...mam import mmap_normalize as mmap_norm_fn, mmap_super_safe as mmap_super_fn
        HAS_TRAFFIC = True
    except ImportError:
        HAS_TRAFFIC = False

    # Get station-to-node and node-to-station mappings
    stationToNode = np.asarray(getattr(sn, 'stationToNode',
                                       np.arange(M))).flatten().astype(int)
    nodeToStation = np.asarray(getattr(sn, 'nodeToStation',
                                       np.arange(M))).flatten().astype(int)

    for it in range(1, it_max + 1):
        if it == 1:
            # Initially form departure processes using scaled service
            for ind in range(M):
                DEP[ind] = {}
                for r in range(K):
                    ist = nodeToStation[ind] if ind < len(nodeToStation) else ind
                    if ist in PH_map and r in PH_map[ist]:
                        ph = PH_map[ist][r]
                        scale = 1.0 / (lambda_k[r] * V[ind, r]) if (
                            lambda_k[r] > FINE_TOL and V[ind, r] > FINE_TOL) else 1.0
                        D0_s, D1_s = map_scale(ph[0], ph[1], scale)
                        DEP[ind][r] = [D0_s, D1_s]
                    else:
                        DEP[ind][r] = [np.array([[-1.0]]), np.array([[1.0]])]

        # Compute arrival processes via traffic equations
        # Simplified: use superposition of departures weighted by routing
        ARV = _compute_arrivals_simple(sn, DEP, V, lambda_k, M, K, PH_map, config)

        QN_prev = QN.copy()

        # Solve each station
        for ist in range(M):
            ind = stationToNode[ist] if ist < len(stationToNode) else ist
            sched = _get_scheduling(sn, ist)

            if sched == SchedStrategy.EXT:
                continue

            # Get node type
            nodetype = None
            if hasattr(sn, 'nodetype') and sn.nodetype is not None:
                nt_arr = np.asarray(sn.nodetype).flatten()
                if ind < len(nt_arr):
                    nodetype = int(nt_arr[ind])

            # Queue node processing
            is_queue = (nodetype == NodeType.Queue if nodetype is not None
                       else sched in [SchedStrategy.FCFS, SchedStrategy.HOL,
                                      SchedStrategy.PS])

            if is_queue:
                # Update throughputs from arrival process
                if ind in ARV and ARV[ind] is not None:
                    arv_mmap = ARV[ind]
                    if HAS_MAP_UTILS and isinstance(arv_mmap, list) and len(arv_mmap) > 1:
                        try:
                            arv_lambdas = mmap_lambda(arv_mmap)
                            for k in range(min(K, len(arv_lambdas))):
                                TN[ist, k] = arv_lambdas[k]
                        except Exception:
                            for k in range(K):
                                TN[ist, k] = lambda_k[k] * V[ist, k]
                    else:
                        for k in range(K):
                            TN[ist, k] = lambda_k[k] * V[ist, k]
                else:
                    for k in range(K):
                        TN[ist, k] = lambda_k[k] * V[ist, k]

                # Solve queue based on scheduling strategy
                classprio = None
                if getattr(sn, 'classprio', None) is not None:
                    classprio = np.asarray(sn.classprio).flatten()
                is_hol_prio = (sched == SchedStrategy.HOL and classprio is not None
                               and K > 1 and np.any(classprio != classprio[0]))
                if is_hol_prio and HAS_MMAPPH1NPPR and ind in ARV \
                        and ARV[ind] is not None and ist in pie and ist in PH_map:
                    # HOL with non-identical priorities: MMAP[K]/PH[K]/1 (BUTools D1=lowest priority vs LINE's lower-value-is-higher convention); mirrors MATLAB solver_mam_basic.
                    if len(np.unique(classprio)) != K:
                        raise RuntimeError(
                            'Solver MAM requires either identical priorities '
                            'or all distinct priorities')
                    iK = np.argsort(-classprio)  # lowest priority first
                    arv_mmap = ARV[ind]
                    # MMAP layout: [D0, D1_total, class-marked matrices...]
                    D_pr = [np.asarray(arv_mmap[0])] + \
                        [np.asarray(arv_mmap[2 + int(k)]) for k in iK]
                    pie_pr = [np.atleast_2d(np.asarray(pie[ist][int(k)])) for k in iK]
                    S_pr = [np.atleast_2d(np.asarray(D0_ph[ist, int(k)])) for k in iK]
                    Qret = MMAPPH1NPPR(D_pr, pie_pr, S_pr, 'ncMoms', 1)
                    if K == 1:
                        Qret = [Qret]
                    for j in range(K):
                        q_j = Qret[j]
                        if hasattr(q_j, '__iter__'):
                            QN[ist, int(iK[j])] = float(np.sum(q_j))
                        else:
                            QN[ist, int(iK[j])] = float(q_j)
                elif sched in [SchedStrategy.FCFS, SchedStrategy.HOL]:
                    # Try MMAPPH1FCFS for accurate analysis
                    if (HAS_MMAPPH1FCFS and ind in ARV and ARV[ind] is not None
                            and ist in pie and ist in PH_map):
                        try:
                            arv_mmap = ARV[ind]
                            # Build MMAP without D1_total: {D0, D3, D4, ..., DK+2}
                            # MATLAB: {ARV{ind}{[1,3:end]}}
                            if len(arv_mmap) > 2:
                                arv_for_fcfs = [arv_mmap[0]] + arv_mmap[2:]
                            else:
                                arv_for_fcfs = arv_mmap

                            pie_list = [pie[ist][k] for k in range(K) if k in pie.get(ist, {})]
                            D0_list = [D0_ph[ist, k] for k in range(K) if (ist, k) in D0_ph]

                            if len(pie_list) == K and len(D0_list) == K:
                                Qret = MMAPPH1FCFS(arv_for_fcfs, pie_list, D0_list,
                                                   'ncMoms', 1, 'ncDistr', 2)
                                if Qret is not None:
                                    for k in range(min(K, len(Qret))):
                                        q_k = Qret[k]
                                        if hasattr(q_k, '__iter__'):
                                            QN[ist, k] = float(np.sum(q_k))
                                        else:
                                            QN[ist, k] = float(q_k) if q_k is not None else 0.0
                        except Exception:
                            # Fallback to M/M/1 approximation
                            for k in range(K):
                                rho_k = TN[ist, k] * S[ist, k] / nservers[ist] if np.isfinite(S[ist, k]) else 0.0
                                rho_total = sum(TN[ist, j] * S[ist, j] / nservers[ist]
                                               for j in range(K) if np.isfinite(S[ist, j]))
                                if rho_total < 1.0 - FINE_TOL:
                                    QN[ist, k] = rho_k / (1.0 - rho_total)
                                else:
                                    QN[ist, k] = rho_k
                    else:
                        # M/M/c approximation fallback
                        for k in range(K):
                            rho_k = TN[ist, k] * S[ist, k] / nservers[ist] if np.isfinite(S[ist, k]) else 0.0
                            rho_total = sum(TN[ist, j] * S[ist, j] / nservers[ist]
                                           for j in range(K) if np.isfinite(S[ist, j]))
                            if rho_total < 1.0 - FINE_TOL:
                                QN[ist, k] = rho_k / (1.0 - rho_total)
                            else:
                                QN[ist, k] = rho_k

                elif sched == SchedStrategy.PS:
                    # Change B: PS queue handling with U/(1-U) formula
                    for k in range(K):
                        if ist in PH_map and k in PH_map[ist]:
                            UN[ist, k] = TN[ist, k] * map_mean(PH_map[ist][k][0],
                                                                 PH_map[ist][k][1])
                        else:
                            UN[ist, k] = TN[ist, k] * S[ist, k] if np.isfinite(S[ist, k]) else 0.0
                    Uden = min(1.0 - FINE_TOL, np.sum(UN[ist, :]))
                    for k in range(K):
                        QN[ist, k] = UN[ist, k] / (1.0 - Uden)

            # Compute utilization, add surrogate delay, compute response times
            for k in range(K):
                if ist in PH_map and k in PH_map[ist]:
                    UN[ist, k] = TN[ist, k] * map_mean(PH_map[ist][k][0],
                                                         PH_map[ist][k][1])
                    # Add surrogate delay for multiserver
                    mean_s = map_mean(PH_map[ist][k][0], PH_map[ist][k][1])
                    QN[ist, k] = QN[ist, k] + TN[ist, k] * (mean_s * nservers[ist]) * (
                        nservers[ist] - 1) / nservers[ist]
                else:
                    if np.isfinite(S[ist, k]):
                        UN[ist, k] = TN[ist, k] * S[ist, k]
                        QN[ist, k] = QN[ist, k] + TN[ist, k] * S[ist, k] * (
                            nservers[ist] - 1) / nservers[ist]
                RN[ist, k] = QN[ist, k] / TN[ist, k] if TN[ist, k] > FINE_TOL else 0.0

        # Check convergence (after iteration 3)
        if it >= 3:
            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change = np.abs(QN - QN_prev) / np.maximum(np.abs(QN_prev), FINE_TOL)
            if np.nanmax(rel_change) < options.iter_tol:
                break

        # Change C: Build departure processes using ETAQA
        for ist in range(M):
            ind = stationToNode[ist] if ist < len(stationToNode) else ist
            sched = _get_scheduling(sn, ist)

            if sched == SchedStrategy.EXT:
                continue

            nodetype = None
            if hasattr(sn, 'nodetype') and sn.nodetype is not None:
                nt_arr = np.asarray(sn.nodetype).flatten()
                if ind < len(nt_arr):
                    nodetype = int(nt_arr[ind])

            is_queue = (nodetype == NodeType.Queue if nodetype is not None
                       else sched in [SchedStrategy.FCFS, SchedStrategy.HOL,
                                      SchedStrategy.PS])

            if is_queue and ind in DEP:
                for r in range(K):
                    if not (ind in ARV and ARV[ind] is not None):
                        # No arrival process available, use scaled PH
                        if ist in PH_map and r in PH_map[ist]:
                            ph = PH_map[ist][r]
                            scale = 1.0 / (lambda_k[r] * V[ind, r]) if (
                                lambda_k[r] > FINE_TOL and V[ind, r] > FINE_TOL) else 1.0
                            D0_s, D1_s = map_scale(ph[0], ph[1], scale)
                            DEP[ind][r] = [D0_s, D1_s]
                        continue

                    # Extract class-r arrival MAP by hiding all other classes
                    arv_mmap = ARV[ind]
                    if (HAS_MAP_UTILS and HAS_QBD_DEPPROC
                            and isinstance(arv_mmap, list) and len(arv_mmap) > 2
                            and ist in PH_map and r in PH_map[ist]):
                        try:
                            # Hide all classes except r (0-indexed)
                            hide_classes = [s for s in range(K) if s != r]
                            # arv_mmap is [D0, D1_total, D1_class0, D1_class1, ...]
                            # The class-specific matrices start at index 2
                            if len(arv_mmap) > 2 and len(hide_classes) > 0:
                                # Extract D0 and per-class matrices
                                D0_arv = arv_mmap[0]
                                D_list_arv = arv_mmap[2:]  # per-class arrival matrices
                                if len(D_list_arv) >= K:
                                    D0_hidden, D_list_hidden = mmap_hide(
                                        D0_arv, D_list_arv, hide_classes)
                                    # Construct MAP: D1 = sum of remaining D_list
                                    D1_arv = sum(D_list_hidden)
                                    A = [D0_hidden, D1_arv]
                                else:
                                    # Fallback: use D0 and D1_total
                                    A = [arv_mmap[0], arv_mmap[1]]
                            else:
                                A = [arv_mmap[0], arv_mmap[1] if len(arv_mmap) > 1
                                     else np.zeros_like(arv_mmap[0])]

                            S_r = PH_map[ist][r]
                            na = A[0].shape[0]
                            ns = S_r[0].shape[0]
                            etaqa_sz = (etaqa_n + 1) * na * ns
                            rho = np.sum(UN[ist, :])

                            # Use ETAQA if state space is manageable and queue is stable
                            if etaqa_sz <= config.space_max and rho < 1.0 - FINE_TOL:
                                try:
                                    if sched in [SchedStrategy.FCFS, SchedStrategy.HOL]:
                                        dep_map = qbd_depproc_etaqa(A, S_r, etaqa_n)
                                    elif sched == SchedStrategy.PS:
                                        dep_map = qbd_depproc_etaqa_ps(A, S_r, etaqa_n)
                                    else:
                                        dep_map = PH_map[ist][r]

                                    dep_map[0], dep_map[1] = map_normalize(
                                        dep_map[0], dep_map[1])
                                    DEP[ind][r] = dep_map
                                except Exception:
                                    # Fall back to scaled service on ETAQA failure
                                    DEP[ind][r] = list(PH_map[ist][r])
                            else:
                                DEP[ind][r] = list(PH_map[ist][r])

                        except Exception:
                            if ist in PH_map and r in PH_map[ist]:
                                DEP[ind][r] = list(PH_map[ist][r])
                            else:
                                DEP[ind][r] = [np.array([[-1.0]]), np.array([[1.0]])]
                    else:
                        # No ETAQA available, use PH service
                        if ist in PH_map and r in PH_map[ist]:
                            DEP[ind][r] = list(PH_map[ist][r])
                        else:
                            DEP[ind][r] = [np.array([[-1.0]]), np.array([[1.0]])]

                    # Scale departure process to match departure rate
                    scale = 1.0 / (lambda_k[r] * V[ind, r]) if (
                        lambda_k[r] > FINE_TOL and V[ind, r] > FINE_TOL) else 1.0
                    DEP[ind][r][0], DEP[ind][r][1] = map_scale(
                        DEP[ind][r][0], DEP[ind][r][1], scale)

    totiter = it

    if options.verbose:
        print(f"\nMAM parametric decomposition completed in {totiter} iterations.")

    # Set system throughputs
    for c in range(C):
        if c not in inchain_map:
            continue
        inchain = np.asarray(inchain_map[c]).flatten().astype(int)
        XN[0, inchain] = lambda_k[inchain]

    CN = np.sum(RN, axis=0).reshape(1, -1)

    result = SolverMAMReturn()
    result.Q = QN
    result.U = UN
    result.R = RN
    result.T = TN
    result.C = CN
    result.X = XN
    result.A = TN.copy()
    result.W = RN.copy()
    result.runtime = time.time() - start_time
    result.method = "dec.mmap"
    result.it = totiter

    return result


def _get_ph_as_map(
    sn: NetworkStruct,
    ist: int,
    k: int,
    S: np.ndarray
) -> Optional[List]:
    """
    Extract PH service distribution as a MAP [D0, D1] for station ist, class k.

    Args:
        sn: Network structure
        ist: Station index
        k: Class index
        S: Service times matrix

    Returns:
        List [D0, D1] or None
    """
    if not hasattr(sn, 'proc') or sn.proc is None:
        return None

    proc = sn.proc
    if ist >= len(proc) or proc[ist] is None:
        return None

    proc_ist = proc[ist]
    if not isinstance(proc_ist, (list, dict)):
        return None

    ph = proc_ist[k] if k < len(proc_ist) else None
    if ph is None:
        return None

    if isinstance(ph, (list, tuple)) and len(ph) >= 2:
        D0 = np.asarray(ph[0], dtype=np.float64)
        D1 = np.asarray(ph[1], dtype=np.float64)
        if D0.ndim < 2:
            D0 = D0.reshape(1, 1)
        if D1.ndim < 2:
            D1 = D1.reshape(1, 1)
        return [D0, D1]

    # Fallback to exponential
    if np.isfinite(S[ist, k]) and S[ist, k] > 0:
        rate = 1.0 / S[ist, k]
        return [np.array([[-rate]]), np.array([[rate]])]

    return None


def _compute_arrivals_simple(
    sn: NetworkStruct,
    DEP: Dict,
    V: np.ndarray,
    lambda_k: np.ndarray,
    M: int,
    K: int,
    PH_map: Dict,
    config: Any
) -> Dict:
    """
    Compute arrival MMAPs at each station from departure MAPs.

    Simplified version of solver_mam_traffic: superimposes departure MAPs
    weighted by routing probabilities.

    Args:
        sn: Network structure
        DEP: Departure MAPs indexed by [node][class]
        V: Visit ratios (M x K)
        lambda_k: Per-class arrival rates
        M: Number of stations
        K: Number of classes
        PH_map: PH service distributions
        config: Solver options

    Returns:
        Dictionary mapping node index to arrival MMAP [D0, D1_total, D1_c0, D1_c1, ...]
    """
    ARV = {}

    try:
        from ...npfqn.traffic import npfqn_traffic_split_cs, npfqn_traffic_merge
    except ImportError:
        pass

    stationToNode = np.asarray(getattr(sn, 'stationToNode',
                                       np.arange(M))).flatten().astype(int)

    for ist in range(M):
        ind = stationToNode[ist] if ist < len(stationToNode) else ist
        sched = _get_scheduling(sn, ist)

        if sched == SchedStrategy.EXT:
            ARV[ind] = None
            continue

        # Build MMAP arrival at this station
        # For each class, build a MAP from departure process
        mmap_list = []
        for r in range(K):
            if ind in DEP and r in DEP[ind]:
                dep_r = DEP[ind][r]
                D0_r = dep_r[0].copy()
                D1_r = dep_r[1].copy()
                # Convert MAP to MMAP with K classes
                # Only class r has non-zero arrivals
                mmap_r = [D0_r]
                # D1_total = D1_r (only this class has arrivals)
                mmap_r.append(D1_r.copy())
                # Per-class arrival matrices
                for s in range(K):
                    if s == r:
                        mmap_r.append(D1_r.copy())
                    else:
                        mmap_r.append(np.zeros_like(D1_r))
                mmap_list.append(mmap_r)

        if len(mmap_list) > 0:
            # Superpose all class arrival MMAPs
            if len(mmap_list) == 1:
                ARV[ind] = mmap_list[0]
            else:
                # Superpose MMAPs using mmap_super
                try:
                    from ...mam import mmap_super
                    # mmap_super accepts list of MMAPs in [D0, D1, ...] format
                    ARV[ind] = mmap_super(mmap_list)
                except Exception:
                    # Fallback: just use the first one
                    ARV[ind] = mmap_list[0]
        else:
            ARV[ind] = None

    return ARV


__all__ = [
    'solver_mam',
    'solver_mam_basic',
    'SolverMAMReturn',
    'SolverMAMOptions',
]
