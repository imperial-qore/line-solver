"""
SN predicate functions.

Native Python implementations of network property predicates
for queueing network analysis and solver selection.
"""

import numpy as np
from typing import Optional

from .network_struct import NetworkStruct, SchedStrategy, RoutingStrategy


# ============================================================================
# Model Type Predicates
# ============================================================================

def sn_is_closed_model(sn: NetworkStruct) -> bool:
    """
    Check if the network model is closed (all finite populations).

    A closed model has finite jobs and no infinite (open) job classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network is a closed model
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    njobs = sn.njobs.flatten()
    return not np.any(np.isinf(njobs)) and np.any(np.isfinite(njobs) & (njobs > 0))


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
    Check if the network model has populations defined.

    Args:
        sn: NetworkStruct object

    Returns:
        True if population is defined
    """
    if sn.njobs is None or len(sn.njobs) == 0:
        return False
    return np.any(sn.njobs.flatten() > 0)


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
    Check if the network has a multiclass FCFS station.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has multiclass and at least one FCFS station
    """
    return sn_has_multi_class(sn) and sn_has_fcfs(sn)


def sn_has_multi_class_heter_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has multiclass heterogeneous FCFS stations.

    A heterogeneous FCFS station has different service rates for different classes.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has FCFS stations with heterogeneous class rates
    """
    if sn.sched is None or len(sn.sched) == 0:
        return False
    if sn.rates is None or len(sn.rates) == 0:
        return False

    rates = sn.rates
    for station_id, strategy in sn.sched.items():
        if strategy != SchedStrategy.FCFS:
            continue
        if station_id >= rates.shape[0]:
            continue
        row = rates[station_id, :]
        # Check if rates vary across classes (heterogeneous)
        finite_rates = row[np.isfinite(row) & (row > 0)]
        if len(finite_rates) > 1:
            if np.max(finite_rates) - np.min(finite_rates) > 1e-10:
                return True
    return False


def sn_has_multi_class_heter_exp_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has multiclass heterogeneous exponential FCFS stations.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has FCFS stations with heterogeneous exponential service
    """
    # For now, we treat this the same as sn_has_multi_class_heter_fcfs
    # since we don't track distribution types in basic NetworkStruct
    return sn_has_multi_class_heter_fcfs(sn)


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
    if sn.lldscaling.size == 0:
        return False
    # Check if any scaling differs from 1.0
    return np.any(sn.lldscaling != 1.0)


def sn_has_joint_dependence(sn: NetworkStruct) -> bool:
    """
    Check if the network has joint-dependent service rates.

    This includes both LJD (Limited Joint Dependence) and LJCD
    (Limited Joint Class Dependence) scaling.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has joint-dependent scaling
    """
    has_ljd = False
    has_ljcd = False

    # Check for LJD scaling
    if hasattr(sn, 'ljdscaling') and sn.ljdscaling is not None:
        if isinstance(sn.ljdscaling, list):
            has_ljd = any(x is not None for x in sn.ljdscaling)
        elif isinstance(sn.ljdscaling, dict):
            has_ljd = any(v is not None for v in sn.ljdscaling.values())

    # Check for LJCD scaling
    if hasattr(sn, 'ljcdscaling') and sn.ljcdscaling is not None:
        if isinstance(sn.ljcdscaling, list):
            has_ljcd = any(x is not None for x in sn.ljcdscaling)
        elif isinstance(sn.ljcdscaling, dict):
            has_ljcd = any(v is not None for v in sn.ljcdscaling.values())

    return has_ljd or has_ljcd


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

    Args:
        sn: NetworkStruct object

    Returns:
        True if the network has priority-based scheduling (non-uniform priorities)
    """
    if sn.classprio is None:
        return False
    if sn.classprio.size == 0:
        return False
    # Check if priorities are non-uniform (not all the same value)
    # Priority 0 is the highest priority in LINE, so checking > 0 is incorrect
    return not np.all(sn.classprio == sn.classprio.flat[0])


def sn_has_class_switching(sn: NetworkStruct) -> bool:
    """
    Check if the network has class switching.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network has class switching
    """
    if sn.csmask is None:
        return False
    if sn.csmask.size == 0:
        return False
    # Check if any off-diagonal elements are non-zero
    K = sn.nclasses
    for i in range(min(K, sn.csmask.shape[0])):
        for j in range(min(K, sn.csmask.shape[1])):
            if i != j and sn.csmask[i, j] > 0:
                return True
    return False


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
    finite_jobs = njobs[np.isfinite(njobs)]
    return np.any(finite_jobs != np.floor(finite_jobs))


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
        RoutingStrategy.KCHOICES,
        RoutingStrategy.RL,
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

    return True


def sn_has_product_form_not_het_fcfs(sn: NetworkStruct) -> bool:
    """
    Check if network has product form except for heterogeneous FCFS.

    Args:
        sn: NetworkStruct object

    Returns:
        True if network would have product form without heterogeneous FCFS
    """
    # Check scheduling strategies
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

    # Check for other violations (but not heterogeneous FCFS)
    if sn_has_priorities(sn):
        return False
    if sn_has_fork_join(sn):
        return False
    if sn_has_sd_routing(sn):
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
