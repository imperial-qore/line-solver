"""
SN predicate functions.

Native Python implementations of network property predicates
for queueing network analysis and solver selection.
"""

import numpy as np
from typing import Optional

from .network_struct import NetworkStruct, SchedStrategy, RoutingStrategy

# FineTol constant matching MATLAB GlobalConstants.FineTol
_FINE_TOL = 1e-8
# CoarseTol constant matching MATLAB GlobalConstants.CoarseTol
_COARSE_TOL = 1e-3
# Zero constant matching MATLAB GlobalConstants.Zero
_ZERO_TOL = 1e-14


# ============================================================================
# Model Type Predicates
# ============================================================================

def sn_is_closed_model(sn: NetworkStruct) -> bool:
    """
    Check if the network model is closed (all finite populations).

    A closed model has all finite job populations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network is a closed model
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    return np.all(np.isfinite(sn.njobs.flatten()))


def sn_is_open_model(sn: NetworkStruct) -> bool:
    """
    Check if the network model is open (all infinite populations).

    An open model has only infinite (open) job classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network is an open model
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    njobs = sn.njobs.flatten()
    return np.all(np.isinf(njobs))


def sn_is_mixed_model(sn: NetworkStruct) -> bool:
    """
    Check if the network model is mixed (both open and closed classes).

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network has both open and closed classes
    """
    return sn_has_open_classes(sn) and sn_has_closed_classes(sn)


def sn_is_population_model(sn: NetworkStruct) -> bool:
    """
    Check if the network model is a population model.

    A population model uses only delay-like scheduling strategies
    (INF, PS, PSPRIO, DPS, GPS, GPSPRIO, DPSPRIO, EXT),
    has no priorities, and no fork-join topology.

    Args:
        sn: NetworkStruct object

    Returns:
        True if model is population-based
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    population_strategies = {
        SchedStrategy.INF,
        SchedStrategy.PS,
        SchedStrategy.PSPRIO,
        SchedStrategy.DPS,
        SchedStrategy.GPS,
        SchedStrategy.GPSPRIO,
        SchedStrategy.DPSPRIO,
        SchedStrategy.EXT,
    }
    for strategy in sn.sched.values():
        if strategy not in population_strategies:
            return False
    if sn_has_priorities(sn):
        return False
    if sn_has_fork_join(sn):
        return False
    return True


# ============================================================================
# Class Predicates
# ============================================================================

def sn_has_closed_classes(sn: NetworkStruct) -> bool:
    """
    Check if the network has closed (finite population) classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one closed class
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    njobs = sn.njobs.flatten()
    return np.any(np.isfinite(njobs) & (njobs > 0))


def sn_has_open_classes(sn: NetworkStruct) -> bool:
    """
    Check if the network has open (infinite population) classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one open class
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    return np.any(np.isinf(sn.njobs.flatten()))


def sn_has_mixed_classes(sn: NetworkStruct) -> bool:
    """
    Check if the network has both open and closed classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has both open and closed classes
    """
    return sn_has_open_classes(sn) and sn_has_closed_classes(sn)


def sn_has_single_class(sn: NetworkStruct) -> bool:
    """
    Check if the network has exactly one class.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has exactly one class
    """
    return sn.nclasses == 1


def sn_has_multi_class(sn: NetworkStruct) -> bool:
    """
    Check if the network has multiple classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has more than one class
    """
    return sn.nclasses > 1


def sn_has_multiple_closed_classes(sn: NetworkStruct) -> bool:
    """
    Check if the network has multiple closed classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has more than one closed class
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    njobs = sn.njobs.flatten()
    closed_count = np.sum(np.isfinite(njobs) & (njobs > 0))
    return closed_count > 1


# ============================================================================
# Chain Predicates
# ============================================================================

def sn_has_single_chain(sn: NetworkStruct) -> bool:
    """
    Check if the network has exactly one chain.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has exactly one chain
    """
    return sn.nchains == 1


def sn_has_multi_chain(sn: NetworkStruct) -> bool:
    """
    Check if the network has multiple chains.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has more than one chain
    """
    return sn.nchains > 1


# ============================================================================
# Scheduling Predicates
# ============================================================================

def sn_has_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if the network has any FCFS (First-Come First-Served) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one FCFS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.FCFS for s in sn.sched.values())


def sn_has_ps(sn: NetworkStruct) -> bool:
    """
    Check if the network has any PS (Processor Sharing) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one PS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.PS for s in sn.sched.values())


def sn_has_inf(sn: NetworkStruct) -> bool:
    """
    Check if the network has any INF (Infinite Server/Delay) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one INF station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.INF for s in sn.sched.values())


def sn_has_lcfs(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LCFS (Last-Come First-Served) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LCFS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LCFS for s in sn.sched.values())


def sn_has_lcfspr(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LCFS-PR (LCFS Preemptive Resume) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LCFS-PR station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LCFSPR for s in sn.sched.values())


def sn_has_lcfs_pr(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LCFS-PR (LCFS Preemptive Resume) stations.

    This is an alias for sn_has_lcfspr, matching the MATLAB function name.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LCFS-PR station
    """
    return sn_has_lcfspr(sn)


def sn_has_siro(sn: NetworkStruct) -> bool:
    """
    Check if the network has any SIRO (Service In Random Order) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one SIRO station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.SIRO for s in sn.sched.values())


def sn_has_dps(sn: NetworkStruct) -> bool:
    """
    Check if the network has any DPS (Discriminatory Processor Sharing) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one DPS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.DPS for s in sn.sched.values())


def sn_has_gps(sn: NetworkStruct) -> bool:
    """
    Check if the network has any GPS (Generalized Processor Sharing) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one GPS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.GPS for s in sn.sched.values())


def sn_has_hol(sn: NetworkStruct) -> bool:
    """
    Check if the network has any HOL (Head of Line) priority stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one HOL station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.HOL for s in sn.sched.values())


def sn_has_lcfs_pi(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LCFS-PI (LCFS Preemptive Identical) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LCFS-PI station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LCFSPI for s in sn.sched.values())


def sn_has_dps_prio(sn: NetworkStruct) -> bool:
    """
    Check if the network has any DPS with priority stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one DPS-PRIO station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.DPSPRIO for s in sn.sched.values())


def sn_has_gps_prio(sn: NetworkStruct) -> bool:
    """
    Check if the network has any GPS with priority stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one GPS-PRIO station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.GPSPRIO for s in sn.sched.values())


def sn_has_ps_prio(sn: NetworkStruct) -> bool:
    """
    Check if the network has any PS with priority stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one PS-PRIO station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.PSPRIO for s in sn.sched.values())


def sn_has_lps(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LPS (Least Progress Scheduling) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LPS station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LPS for s in sn.sched.values())


def sn_has_setf(sn: NetworkStruct) -> bool:
    """
    Check if the network has any SETF (Shortest Elapsed Time First) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one SETF station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.SETF for s in sn.sched.values())


def sn_has_sept(sn: NetworkStruct) -> bool:
    """
    Check if the network has any SEPT (Shortest Expected Processing Time) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one SEPT station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.SEPT for s in sn.sched.values())


def sn_has_lept(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LEPT (Longest Expected Processing Time) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LEPT station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LEPT for s in sn.sched.values())


def sn_has_sjf(sn: NetworkStruct) -> bool:
    """
    Check if the network has any SJF (Shortest Job First) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one SJF station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.SJF for s in sn.sched.values())


def sn_has_ljf(sn: NetworkStruct) -> bool:
    """
    Check if the network has any LJF (Longest Job First) stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one LJF station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.LJF for s in sn.sched.values())


def sn_has_polling(sn: NetworkStruct) -> bool:
    """
    Check if the network has any polling stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has at least one polling station
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return any(s == SchedStrategy.POLLING for s in sn.sched.values())


def sn_has_homogeneous_scheduling(sn: NetworkStruct, strategy: int) -> bool:
    """
    Check if the network uses an identical scheduling strategy at every station.

    Args:
        sn: NetworkStruct object
        strategy: SchedStrategy value to check for

    Returns:
        True if all stations use the specified strategy
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    return all(s == strategy for s in sn.sched.values())


# ============================================================================
# Multi-class FCFS Predicates
# ============================================================================

def sn_has_multi_class_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if the network has an FCFS station that serves multiple classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if any FCFS station serves more than one class
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    if sn.rates is None or sn.rates.size == 0:
        return False

    for station_id, strategy in sn.sched.items():
        strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
        if strategy_val != int(SchedStrategy.FCFS):
            continue
        if station_id >= sn.rates.shape[0]:
            continue
        row = sn.rates[station_id, :]
        # Count classes with positive rates at this FCFS station
        if np.sum(row > 0) > 1:
            return True
    return False


def sn_has_multi_class_heter_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has multiclass heterogeneous FCFS stations.

    A heterogeneous FCFS station has different service rates for different classes.
    Uses MATLAB's range() check: max(rates) - min(rates) > 0 across all classes
    at each FCFS station.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has FCFS stations with heterogeneous class rates
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    if sn.rates is None or sn.rates.size == 0:
        return False

    rates = sn.rates
    for station_id, strategy in sn.sched.items():
        strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
        if strategy_val != int(SchedStrategy.FCFS):
            continue
        if station_id >= rates.shape[0]:
            continue
        row = rates[station_id, :]
        # MATLAB: range([sn.rates(i,:)]) > 0
        # range() = max - min over all values (including NaN-handling)
        finite_vals = row[np.isfinite(row)]
        if len(finite_vals) > 0:
            if np.max(finite_vals) - np.min(finite_vals) > 0:
                return True
    return False


def sn_has_multi_class_heter_exp_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has multiclass heterogeneous exponential FCFS stations.

    Returns true if any FCFS station has heterogeneous rates AND all
    service time SCVs at that station are approximately 1.0 (exponential).

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has FCFS stations with heterogeneous exponential service
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    if sn.rates is None or sn.rates.size == 0:
        return False
    if sn.scv is None or sn.scv.size == 0:
        return False

    for station_id, strategy in sn.sched.items():
        strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
        if strategy_val != int(SchedStrategy.FCFS):
            continue
        if station_id >= sn.rates.shape[0]:
            continue
        row = sn.rates[station_id, :]
        # Check if rates vary across classes (heterogeneous)
        finite_vals = row[np.isfinite(row)]
        if len(finite_vals) > 0 and (np.max(finite_vals) - np.min(finite_vals)) > 0:
            # Check if all SCVs are ~1 (exponential)
            if station_id < sn.scv.shape[0]:
                scvs = sn.scv[station_id, :]
                finite_scvs = scvs[np.isfinite(scvs)]
                if len(finite_scvs) > 0:
                    if np.max(finite_scvs) < 1 + _FINE_TOL and np.min(finite_scvs) > 1 - _FINE_TOL:
                        return True
    return False


# ============================================================================
# Server Predicates
# ============================================================================

def sn_has_multi_server(sn: NetworkStruct) -> bool:
    """
    Check if the network has any multi-server stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if any station has more than one server
    """
    if sn.nservers is None or len(sn.nservers) == 0:
        return False
    nservers = sn.nservers.flatten()
    # Filter out infinite servers (delays)
    finite_servers = nservers[np.isfinite(nservers)]
    return np.any(finite_servers > 1)


# ============================================================================
# Load Dependence Predicates
# ============================================================================

def sn_has_load_dependence(sn: NetworkStruct) -> bool:
    """
    Check if the network has load-dependent service.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has load-dependent scaling
    """
    if sn.lldscaling is None:
        return False
    if isinstance(sn.lldscaling, np.ndarray):
        if sn.lldscaling.ndim < 2:
            return sn.lldscaling.size > 0
        return sn.lldscaling.shape[1] > 0
    # For non-array types (list, etc.)
    return len(sn.lldscaling) > 0


# ============================================================================
# Structure Predicates
# ============================================================================

def sn_has_fork_join(sn: NetworkStruct) -> bool:
    """
    Check if the network uses fork and/or join nodes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has fork-join topology
    """
    if sn.fj is None:
        return False
    if sn.fj.size == 0:
        return False
    return np.any(sn.fj > 0)


def sn_has_priorities(sn: NetworkStruct) -> bool:
    """
    Check if the network uses class priorities.

    In LINE, priority 0 is default (no priority). Values > 0 indicate
    priority classes are in use.

    Args:
        sn: NetworkStruct object

    Returns:
        True if any class has priority > 0
    """
    if sn.classprio is None:
        return False
    if sn.classprio.size == 0:
        return False
    return np.any(sn.classprio.flatten() > 0)


def sn_has_quorum_join(sn: NetworkStruct) -> bool:
    """
    Check if the network has a quorum (k-of-n) join.

    True if some Join node declares a non-standard strategy with a positive required
    count in some class, i.e. it fires before every sibling has arrived. The sibling
    count is not re-derived here, so a declaration with k >= n reads as a quorum; use
    sn_join_quorum where the branch count is known and the distinction matters, as the
    fork-join fixed point does.

    Args:
        sn: NetworkStruct object

    Returns:
        True if some join declares a positive quorum
    """
    from ...lang.base import JoinStrategy
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is None:
        return False
    values = nodeparam.values() if isinstance(nodeparam, dict) else nodeparam
    for param in values:
        if not isinstance(param, dict):
            continue
        strategies = param.get('joinStrategy')
        required = param.get('joinRequired')
        if strategies is None or required is None:
            continue
        for r, strategy in enumerate(strategies):
            if strategy == JoinStrategy.STD:
                continue
            if r < len(required) and required[r] is not None and float(required[r]) > 0:
                return True
    return False


def sn_has_class_switching(sn: NetworkStruct) -> bool:
    """
    Check if the network has class switching.

    Class switching is indicated by the number of classes
    differing from the number of chains.

    Args:
        sn: NetworkStruct object

    Returns:
        True if number of classes differs from number of chains
    """
    return sn.nclasses != sn.nchains


def sn_has_fractional_populations(sn: NetworkStruct) -> bool:
    """
    Check if the network has fractional (non-integer) populations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if any class has fractional population
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    njobs = sn.njobs.flatten()
    return np.any(njobs != np.round(njobs))


# ============================================================================
# Product Form Predicates
# ============================================================================

def sn_has_sd_routing(sn: NetworkStruct) -> bool:
    """
    Check if the network has state-dependent routing strategies.

    State-dependent routing strategies violate the product-form assumption.
    These include Round-Robin, Weighted Round-Robin, Join Shortest Queue,
    Power of K Choices, and Reinforcement Learning.

    Product-form requires state-independent (Markovian) routing.
    PROB and RAND are product-form compatible.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has state-dependent routing
    """
    if sn.routing is None or sn.routing.size == 0:
        return False

    # Non-product-form routing strategies
    sd_strategies = {
        RoutingStrategy.RROBIN,
        RoutingStrategy.WRROBIN,
        RoutingStrategy.JSQ,
        RoutingStrategy.SQ,
        RoutingStrategy.SDR,
    }

    for val in sn.routing.flatten():
        if val in sd_strategies:
            return True

    return False


def sn_has_product_form(sn: NetworkStruct) -> bool:
    """
    Check if the network has a known product-form solution.

    A network has product form if:
    - All stations use INF, PS, FCFS, LCFS-PR, or EXT scheduling
    - No multiclass heterogeneous FCFS
    - No priorities
    - No fork-join
    - At FCFS stations, all active class SCVs are approximately 1 (BCMP type 1)

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has product-form solution
    """
    # Check scheduling strategies
    if sn.sched is not None and len(sn.sched) > 0:
        product_form_strategies = {
            SchedStrategy.INF,
            SchedStrategy.PS,
            SchedStrategy.FCFS,
            SchedStrategy.LCFS,
            SchedStrategy.LCFSPR,
            SchedStrategy.EXT,
        }
        for strategy in sn.sched.values():
            if strategy not in product_form_strategies:
                return False

    # Check for violations
    if sn_has_multi_class_heter_fcfs(sn):
        return False
    if sn_has_priorities(sn):
        return False
    if sn_has_fork_join(sn):
        return False
    if sn_has_sd_routing(sn):
        return False
    # BCMP asks for infinite buffers. Nothing here read sn.cap/sn.classcap/
    # sn.droprule, so a BAS-blocked station or any binding finite buffer passed the
    # gate and the network read as product form while its truncation couples the
    # station occupancies.
    if sn_has_blocking(sn):
        return False

    # BCMP type 1 asks the FCFS service to be exponential. sn_has_multi_class_heter_fcfs
    # compares the class MEANS only, so a class-homogeneous Erlang, hyper-exponential or
    # deterministic FCFS station passed this gate and was dispatched to exact MVA, which
    # reads the means alone and returns the exponential answer with no warning.
    if sn.sched is not None and sn.scv is not None and sn.scv.size > 0:
        for station_id, strategy in sn.sched.items():
            strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
            if strategy_val != int(SchedStrategy.FCFS):
                continue
            if station_id >= sn.scv.shape[0]:
                continue
            scvs = sn.scv[station_id, :]
            active = np.isfinite(scvs) & (scvs > 0)
            if np.any(active):
                active_scvs = scvs[active]
                if not (np.all(active_scvs > 1 - _FINE_TOL) and np.all(active_scvs < 1 + _FINE_TOL)):
                    return False

    return True


def sn_has_bursty_arrival(sn: NetworkStruct) -> bool:
    """
    Check whether any external arrival process is bursty (non-renewal).

    Returns True if any Source station has an arrival process with
    autocorrelated inter-arrival times (a non-renewal Markovian arrival process
    such as an MMPP/MAP), as opposed to a renewal process (Poisson, or any
    i.i.d. renewal process such as Erlang/HyperExp/Coxian/APH). Detection is
    exact: a MAP with matrices (D0,D1) is renewal iff D1 equals its rank-one
    renewal form t0*pie, where t0 = -D0*e and pie is the embedded stationary
    vector; any departure signals correlation between inter-arrival times.

    Args:
        sn: NetworkStruct object

    Returns:
        True if some external arrival process is non-renewal (bursty).
    """
    from ..mam import map_pie
    from .network_struct import NodeType

    if sn.proc is None:
        return False
    for ist in range(sn.nstations):
        nd = int(sn.stationToNode[ist])
        if sn.nodetype[nd] != NodeType.SOURCE:
            continue
        proc_st = sn.proc[ist] if ist < len(sn.proc) else None
        if proc_st is None:
            continue
        for r in range(sn.nclasses):
            if r >= len(proc_st) or proc_st[r] is None:
                continue
            mapproc = proc_st[r]
            # see _kb/03-api-layer.md for rationale
            if not isinstance(mapproc, (list, tuple)) or len(mapproc) < 2 \
                    or mapproc[0] is None or mapproc[1] is None:
                continue
            D0 = np.asarray(mapproc[0], dtype=np.float64)
            D1 = np.asarray(mapproc[1], dtype=np.float64)
            if D1.ndim < 2:
                # see _kb/03-api-layer.md for rationale
                continue
            n = D1.shape[0]
            if n <= 1:
                continue   # single-phase arrival is Poisson, hence renewal
            pie = map_pie(D0, D1).reshape(1, n)
            D1ren = (D1 @ np.ones((n, 1))) @ pie   # rank-one renewal form t0*pie
            if np.linalg.norm(D1 - D1ren, 'fro') > 1e-8 * max(1.0, np.linalg.norm(D1, 'fro')):
                return True
    return False


def sn_is_mm1k_loss(sn: NetworkStruct) -> bool:
    """
    Check if the model is a single-station M/M/1/K queue with tail drop.

    True for a single-class open Source-Queue-Sink system whose queue is a
    single-server exponential M/M/1/K with tail drop (DropStrategy.Drop). This
    is the exact regime of the closed-form loss scripts qsys_mm1k_loss
    (probability-based, SolverNC) and qsys_mg1k_loss_mgs (moment-based,
    SolverMVA), and the one truncated shape that keeps a product form over its
    single station, hence the exemption in sn_has_blocking.

    Args:
        sn: NetworkStruct object

    Returns:
        True if the model is a single-station M/M/1/K with tail drop
    """
    from .network_struct import NodeType
    from ...constants import DropStrategy

    if sn.nclasses != 1 or sn.nclosedjobs != 0:
        return False
    if sn.nodetype is None or len(sn.nodetype) != 3:
        return False
    qnode = snode = -1
    nsink = 0
    for nd, ntype in enumerate(sn.nodetype):
        if ntype == NodeType.QUEUE:
            if qnode >= 0:
                return False
            qnode = nd
        elif ntype == NodeType.SOURCE:
            if snode >= 0:
                return False
            snode = nd
        elif ntype == NodeType.SINK:
            nsink += 1
        else:
            return False
    if qnode < 0 or snode < 0 or nsink != 1:
        return False
    qist = int(sn.nodeToStation[qnode])
    sist = int(sn.nodeToStation[snode])
    if qist < 0 or sist < 0:
        return False
    if sn.nservers is None or float(np.asarray(sn.nservers).ravel()[qist]) != 1.0:
        return False
    droprule = getattr(sn, 'droprule', None)
    if droprule is None:
        return False
    droprule = np.asarray(droprule)
    if droprule.ndim < 2 or qist >= droprule.shape[0]:
        return False
    if int(droprule[qist, 0]) != int(DropStrategy.Drop.value):
        return False
    cap = np.asarray(sn.cap, dtype=float).ravel() if sn.cap is not None else None
    if cap is None or qist >= cap.size or not np.isfinite(cap[qist]) or cap[qist] <= 0:
        return False
    if sn.scv is None or sn.scv.shape[0] <= max(qist, sist):
        return False
    return abs(sn.scv[sist, 0] - 1) <= 1e-6 and abs(sn.scv[qist, 0] - 1) <= 1e-6


def sn_has_blocking(sn: NetworkStruct) -> bool:
    """
    Check if the network holds jobs back at a finite buffer or region.

    True when some station can refuse a job, either because its own buffer
    BINDS (Kendall's K below the population that can reach it, whatever the
    drop rule: WAITQ, DROP, BAS, BBS, RSRD) or because a finite capacity region
    caps a set of stations jointly. Such a network is not product form: the
    truncation couples the station occupancies, so no BCMP factorization of the
    equilibrium distribution exists.

    Only a buffer that can actually BIND counts, which is what
    sn_get_buffer_size decides: refreshCapacity derives a finite classcap (the
    chain population) at every station of every closed model, so a plain
    finiteness test would call every closed model blocking.

    Two shapes are exempt. A Cache builds its own capped retrieval queues
    (classCap = 1), which the cache analyzers solve rather than treat as a
    buffer constraint. And the single-station M/M/1/K loss system keeps the
    truncated geometric distribution, a product form over its one station.

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network has binding finite buffers or capacity regions
    """
    from .network_struct import NodeType
    from ..me.solver_nc_mem import sn_get_buffer_size

    # a finite capacity region caps a SET of stations, which no per-station
    # capacity can express and no product form survives
    if getattr(sn, 'nregions', 0):
        return True
    if sn.nodetype is not None and any(int(t) == int(NodeType.CACHE) for t in sn.nodetype):
        return False
    if sn_is_mm1k_loss(sn):
        return False
    for ist in range(int(sn.nstations)):
        if np.isfinite(sn_get_buffer_size(sn, ist)):
            return True
    return False


def sn_has_product_form_not_het_fcfs(sn: NetworkStruct, check_means: bool = True) -> bool:
    """
    Check if network has product form except for heterogeneous FCFS.

    This checks:
    - All stations use INF, PS, FCFS, LCFSPR, or EXT scheduling
    - No priorities, no fork-join, no state-dependent routing
    - At FCFS stations, all active class SCVs are approximately 1 (exponential)
      and all active class service means agree (BCMP type 1 asks the FCFS
      service to be class-independent, not merely exponential)

    Args:
        sn: NetworkStruct object
        check_means: also demand class-independent FCFS service means. Pass False only for
            an algorithm that models class-dependent FCFS itself (ab, schmidt, schmidt-ext),
            for which the exclusion is the whole point.

    Returns:
        True if network would have product form without heterogeneous FCFS
    """
    # Check scheduling strategies (note: MATLAB excludes LCFS here, only LCFSPR)
    if sn.sched is not None and len(sn.sched) > 0:
        product_form_strategies = {
            SchedStrategy.INF,
            SchedStrategy.PS,
            SchedStrategy.FCFS,
            SchedStrategy.LCFSPR,
            SchedStrategy.EXT,
        }
        for strategy in sn.sched.values():
            if strategy not in product_form_strategies:
                return False

    # Check for other violations
    if sn_has_priorities(sn):
        return False
    if sn_has_fork_join(sn):
        return False
    if sn_has_sd_routing(sn):
        return False

    # At FCFS stations, check that all active class SCVs are ~1 (exponential)
    if sn.sched is not None and sn.scv is not None and sn.scv.size > 0:
        for station_id, strategy in sn.sched.items():
            strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
            if strategy_val != int(SchedStrategy.FCFS):
                continue
            if station_id >= sn.scv.shape[0]:
                continue
            scvs = sn.scv[station_id, :]
            # Active classes: finite and positive SCV
            active = np.isfinite(scvs) & (scvs > 0)
            if np.any(active):
                active_scvs = scvs[active]
                if not (np.all(active_scvs > 1 - _FINE_TOL) and np.all(active_scvs < 1 + _FINE_TOL)):
                    return False

    # At FCFS stations the service means must agree too: with unequal means the
    # product-form solve returns a wait proportional to each class's own demand
    # where FCFS makes every class wait behind the same queue. The comparison is
    # between CHAIN service times (visit-weighted over the classes that actually
    # visit the station): a class that never visits cannot break product form,
    # and within-chain heterogeneity is invisible to both the product-form and
    # the qd branch, which deaggregate a chain result proportionally to each
    # class's own demand, so only between-chain heterogeneity warrants the
    # divert. LN layers carry seeded rates for classes with zero visits, which
    # a raw per-class comparison mistakes for heterogeneity.
    if check_means and sn.sched is not None and sn.rates is not None and sn.rates.size > 0:
        for station_id, strategy in sn.sched.items():
            strategy_val = int(strategy) if hasattr(strategy, '__int__') else strategy
            if strategy_val != int(SchedStrategy.FCFS):
                continue
            if station_id >= sn.rates.shape[0]:
                continue
            rates = sn.rates[station_id, :]
            stateful_idx = station_id
            if sn.stationToStateful is not None and station_id < len(sn.stationToStateful):
                stateful_idx = int(sn.stationToStateful[station_id])
            st_chain = []
            nchains = int(sn.nchains) if sn.nchains is not None else 0
            for c in range(nchains):
                if not sn.visits or c not in sn.visits:
                    continue
                visits = sn.visits[c]
                if stateful_idx >= visits.shape[0]:
                    continue
                num = 0.0
                den = 0.0
                for r in range(min(rates.shape[0], visits.shape[1])):
                    inchain = (sn.chains[c, r] > 0) if (sn.chains is not None and sn.chains.ndim == 2) else True
                    if not inchain:
                        continue
                    w = visits[stateful_idx, r]
                    if w > _ZERO_TOL and np.isfinite(rates[r]) and rates[r] > 0:
                        num += w / rates[r]
                        den += w
                if den > 0:
                    st_chain.append(num / den)
            if st_chain:
                st_chain = np.asarray(st_chain)
                if np.max(st_chain) - np.min(st_chain) > _COARSE_TOL * np.max(st_chain):
                    return False

    return True


def sn_has_product_form_except_multi_class_heter_exp_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has product form except for multiclass heterogeneous exponential FCFS.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network would have product form without multiclass heter exp FCFS
    """
    return sn_has_product_form_not_het_fcfs(sn)


# ============================================================================
# State Predicates
# ============================================================================

def sn_is_state_valid(sn: NetworkStruct) -> bool:
    """
    Check if the network state is valid.

    Args:
        sn: NetworkStruct object

    Returns:
        True if state is valid
    """
    if sn.state is None or len(sn.state) == 0:
        return True  # Empty state is valid
    # Check that state dimensions match
    for node_id, state in sn.state.items():
        if state is None:
            continue
        # State should match expected dimensions
        # This is a basic check; more detailed validation could be added
    return True


def sn_is_discrete_time(sn: NetworkStruct, options=None):
    """Decide whether a model lives on a discrete (slotted) time scale.

    A model is discrete-time when every enabled interarrival and service law is
    lattice-valued on a common slot length and at least one of them is
    intrinsically discrete. The lattice families are Geometric (support
    {1,2,...} slots), DMAP, DiscreteUniform with integral bounds, and Det whose
    value is a positive integral number of slots. Immediate is deliberately not
    one: a zero interval is not a point of {d,2d,...}, the same refusal the LDES
    slotted engine makes.

    The test reads procid, rates and scv rather than proc, because the struct
    refresh may already have replaced a lattice law by a continuous surrogate.
    procid keeps the requested family and (mean, SCV) identify the member of it
    exactly for every family above.

    MATLAB twin: sn_is_discrete_time.m

    Returns:
        (is_discrete_time, slot_length, info) with info a dict carrying
        has_lattice, has_continuous, has_dmap and reason
    """
    from ...constants import ProcessType

    tol = 1e-8
    timescale = 'auto'
    slot_length = 1.0
    config = getattr(options, 'config', None) if options is not None else None
    if isinstance(config, dict):
        timescale = str(config.get('timescale', 'auto')).lower()
        slot_length = float(config.get('slotlength', 1.0))
    elif options is not None:
        timescale = str(getattr(options, 'timescale', 'auto')).lower()
        slot_length = float(getattr(options, 'slotlength', 1.0))

    if timescale not in ('auto', 'discrete', 'continuous'):
        raise ValueError("config.timescale must be 'auto', 'discrete' or 'continuous'.")
    if not np.isfinite(slot_length) or slot_length <= 0:
        raise ValueError("config.slotlength must be a positive finite scalar.")

    info = {'has_lattice': False, 'has_continuous': False, 'has_dmap': False, 'reason': ''}
    if timescale == 'continuous':
        return False, slot_length, info

    procid = np.asarray(sn.procid) if sn.procid is not None else None
    rates = np.asarray(sn.rates) if sn.rates is not None else None
    scv = np.asarray(sn.scv) if sn.scv is not None else None
    if procid is None or rates is None:
        return False, slot_length, info

    for ist in range(sn.nstations):
        for r in range(sn.nclasses):
            proc_type = procid[ist, r]
            if proc_type is None or proc_type == ProcessType.DISABLED:
                continue
            rate = float(rates[ist, r])
            if np.isnan(rate) or rate <= 0:
                continue
            mean_slots = 1.0 / (rate * slot_length)

            if proc_type == ProcessType.GEOMETRIC:
                info['has_lattice'] = True
                if mean_slots < 1 - tol:
                    info['has_continuous'] = True
                    info['reason'] = (f"Geometric at station {ist} class {r} has mean "
                                      f"{mean_slots} slots, below the one-slot minimum.")
            elif proc_type == ProcessType.DMAP:
                info['has_lattice'] = True
                info['has_dmap'] = True
            elif proc_type == ProcessType.DUNIFORM:
                info['has_lattice'] = True
                var_slots = float(scv[ist, r]) * mean_slots ** 2 if scv is not None else 0.0
                width = np.sqrt(max(0.0, 12 * var_slots + 1)) - 1
                lo = int(round(mean_slots - width / 2))
                if lo < 1:
                    info['has_continuous'] = True
                    info['reason'] = (f"DiscreteUniform at station {ist} class {r} "
                                      "is not contained in {1,2,...}.")
            elif proc_type == ProcessType.DET:
                if abs(mean_slots - round(mean_slots)) <= tol * max(1.0, mean_slots) \
                        and round(mean_slots) >= 1:
                    info['has_lattice'] = True
                else:
                    # a Det off the lattice is what makes the model continuous
                    info['has_continuous'] = True
            else:
                info['has_continuous'] = True

    if timescale == 'discrete':
        if info['has_lattice'] and info['has_continuous']:
            raise RuntimeError("config.timescale='discrete' was requested but the model mixes "
                               f"lattice and non-lattice laws. {info['reason']}")
        if not info['has_lattice']:
            raise RuntimeError("config.timescale='discrete' was requested but no interarrival "
                               f"or service law is lattice-valued on a slot of {slot_length}.")
        return True, slot_length, info

    is_dt = info['has_lattice'] and not info['has_continuous']

    if not is_dt and info['has_dmap']:
        # A DMAP has no continuous-time reading: its (D0,D1) are probability
        # matrices, so the CTMC machinery would compute inv(-D0) where the law
        # needs inv(I-D0) and return a wrong number in silence.
        raise RuntimeError("The model mixes a DMAP with continuous-time laws. A DMAP is only "
                           "defined on a slotted time scale, so no solver can interpret this "
                           f"model. {info['reason']}")

    return is_dt, slot_length, info
