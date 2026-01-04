"""
Fork-Join (FJ) Topology Analysis.

Native Python implementations for analyzing Fork-Join queueing systems,
including topology validation, parameter extraction, and distribution mapping.

Key functions:
    fj_isfj: Check if network has Fork-Join topology
    fj_extract_params: Extract Fork-Join parameters from network
    fj_dist2fj: Convert distribution to FJ-compatible format

References:
    Original MATLAB: matlab/src/api/fj/fj_*.m
    Qiu, Pérez, Harrison, "Beyond the Mean in Fork-Join Queues", 2015
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any, List
from dataclasses import dataclass


@dataclass
class FjParams:
    """Fork-Join network parameters."""
    K: int  # Number of parallel queues
    lambda_val: float  # Arrival rate
    mu: np.ndarray  # Service rates at each queue
    cv: np.ndarray  # Coefficient of variation at each queue
    dist_type: str  # Distribution type ('Exp', 'Erlang', 'HyperExp')


def fj_isfj(sn) -> bool:
    """
    Check if a network has valid Fork-Join topology.

    A valid FJ topology has:
    - Source → Fork → K Queues → Join → Sink
    - Homogeneous service distributions across parallel queues (optional)
    - Open classes only

    Args:
        sn: Network structure

    Returns:
        True if network has valid Fork-Join topology

    References:
        Original MATLAB: matlab/src/api/fj/fj_isfj.m
    """
    if not hasattr(sn, 'nodetype'):
        return False

    from ..sn.network_struct import NodeType

    # Check for required node types
    has_source = np.any(sn.nodetype == NodeType.Source)
    has_sink = np.any(sn.nodetype == NodeType.Sink)
    has_fork = np.any(sn.nodetype == NodeType.Fork)
    has_join = np.any(sn.nodetype == NodeType.Join)

    if not (has_source and has_sink and has_fork and has_join):
        return False

    # Check for at least one queue
    queue_count = np.sum(sn.nodetype == NodeType.Queue)
    if queue_count < 1:
        return False

    # Check for open classes only
    if hasattr(sn, 'njobs'):
        if np.any(np.isfinite(sn.njobs) & (sn.njobs > 0)):
            # Has closed classes
            return False

    return True


def fj_extract_params(sn) -> Optional[FjParams]:
    """
    Extract Fork-Join parameters from network structure.

    Extracts the parameters needed for FJ analysis including
    arrival rate, service rates, and distribution characteristics.

    Args:
        sn: Network structure

    Returns:
        FjParams if valid FJ network, None otherwise

    References:
        Original MATLAB: matlab/src/api/fj/fj_extract_params.m
    """
    if not fj_isfj(sn):
        return None

    from ..sn.network_struct import NodeType

    # Find source and extract arrival rate
    source_idx = np.where(sn.nodetype == NodeType.Source)[0]
    if len(source_idx) == 0:
        return None

    # Get arrival rate
    lambda_val = 0.0
    if hasattr(sn, 'rates'):
        lambda_val = np.sum(sn.rates[source_idx[0], :])
    elif hasattr(sn, 'proc') and sn.proc is not None:
        for r in range(sn.nclasses):
            if sn.proc[source_idx[0], r] is not None:
                rate = sn.proc[source_idx[0], r].get('rate', 0)
                if rate > 0:
                    lambda_val += rate

    # Find queues and extract service parameters
    queue_indices = np.where(sn.nodetype == NodeType.Queue)[0]
    K = len(queue_indices)

    if K == 0:
        return None

    mu = np.zeros(K)
    cv = np.ones(K)  # Default CV = 1 (exponential)
    dist_type = 'Exp'

    for i, q_idx in enumerate(queue_indices):
        ist = sn.nodeToStation[q_idx]
        if ist >= 0 and hasattr(sn, 'proc') and sn.proc is not None:
            for r in range(sn.nclasses):
                if sn.proc[ist, r] is not None:
                    proc = sn.proc[ist, r]
                    mean = proc.get('mean', 0)
                    if mean > 0:
                        mu[i] = 1.0 / mean

                    # Check distribution type
                    if 'type' in proc:
                        dist_type = proc['type']
                    if 'cv' in proc:
                        cv[i] = proc['cv']
                    break

    return FjParams(
        K=K,
        lambda_val=lambda_val,
        mu=mu,
        cv=cv,
        dist_type=dist_type
    )


def fj_dist2fj(dist_type: str, params: Dict[str, Any]) -> Dict[str, Any]:
    """
    Convert a distribution to FJ-compatible format.

    Converts various distribution specifications to the format
    expected by FJ analysis routines.

    Supported distributions:
    - Exp: Exponential (rate)
    - Erlang: Erlang-k (k, rate)
    - HyperExp: 2-phase hyperexponential (p, mu1, mu2)
    - MAP: 2-state MAP (D0, D1)

    Args:
        dist_type: Distribution type string
        params: Distribution parameters

    Returns:
        Dictionary with FJ-compatible parameters:
            - type: Distribution type
            - mean: Mean value
            - cv: Coefficient of variation
            - phases: Number of phases (if applicable)
            - representation: Matrix representation (if applicable)

    References:
        Original MATLAB: matlab/src/api/fj/fj_dist2fj.m
    """
    result = {
        'type': dist_type,
        'mean': 0.0,
        'cv': 1.0,
        'phases': 1,
        'representation': None
    }

    if dist_type.lower() in ['exp', 'exponential']:
        rate = params.get('rate', params.get('mu', 1.0))
        result['mean'] = 1.0 / rate
        result['cv'] = 1.0
        result['phases'] = 1

    elif dist_type.lower() in ['erlang', 'erlangk']:
        k = params.get('k', params.get('phases', 2))
        rate = params.get('rate', params.get('mu', 1.0))
        result['mean'] = k / rate
        result['cv'] = 1.0 / np.sqrt(k)
        result['phases'] = k

    elif dist_type.lower() in ['hyperexp', 'hyperexponential', 'h2']:
        p = params.get('p', 0.5)
        mu1 = params.get('mu1', params.get('rate1', 1.0))
        mu2 = params.get('mu2', params.get('rate2', 1.0))

        mean1 = 1.0 / mu1
        mean2 = 1.0 / mu2
        result['mean'] = p * mean1 + (1 - p) * mean2

        # Second moment
        m2 = 2 * (p / mu1 ** 2 + (1 - p) / mu2 ** 2)
        var = m2 - result['mean'] ** 2
        result['cv'] = np.sqrt(var) / result['mean'] if result['mean'] > 0 else 1.0
        result['phases'] = 2

    elif dist_type.lower() in ['map', 'map2']:
        D0 = np.asarray(params.get('D0', [[-1]]))
        D1 = np.asarray(params.get('D1', [[1]]))

        # Compute MAP mean
        D = D0 + D1
        n = D.shape[0]
        pi = np.ones(n) / n
        for _ in range(100):
            pi_new = pi @ np.linalg.matrix_power(np.eye(n) + D / 100, 100)
            pi_new /= np.sum(pi_new)
            if np.linalg.norm(pi_new - pi) < 1e-12:
                break
            pi = pi_new

        lambda_val = pi @ D1 @ np.ones(n)
        result['mean'] = 1.0 / lambda_val if lambda_val > 0 else 0
        result['representation'] = {'D0': D0, 'D1': D1}
        result['phases'] = n

    elif dist_type.lower() in ['det', 'deterministic', 'constant']:
        value = params.get('value', params.get('mean', 1.0))
        result['mean'] = value
        result['cv'] = 0.0
        result['phases'] = np.inf  # Infinite phases for deterministic

    else:
        # Unknown distribution - use provided mean and CV if available
        result['mean'] = params.get('mean', 1.0)
        result['cv'] = params.get('cv', 1.0)

    return result


from .analytical import (
    # Basic functions
    fj_harmonic,
    fj_synch_delay,
    # Response time approximations
    fj_respt_2way,
    fj_respt_nt,
    fj_respt_vm,
    fj_respt_varki,
    fj_rmax,
    fj_bounds,
    fj_rmax_evd,
    fj_rmax_erlang,
    # Expected maximum functions
    fj_xmax_exp,
    fj_xmax_2,
    fj_xmax_erlang,
    fj_xmax_hyperexp,
    fj_xmax_approx,
    fj_xmax_emma,
    fj_xmax_normal,
    fj_xmax_pareto,
    # Order statistics and bounds
    fj_order_stat,
    fj_gk_bound,
    fj_char_max,
    fj_quantile,
    fj_sm_tput,
)

__all__ = [
    # Topology functions
    'FjParams',
    'fj_isfj',
    'fj_extract_params',
    'fj_dist2fj',
    # Basic functions
    'fj_harmonic',
    'fj_synch_delay',
    # Response time approximations
    'fj_respt_2way',
    'fj_respt_nt',
    'fj_respt_vm',
    'fj_respt_varki',
    'fj_rmax',
    'fj_bounds',
    'fj_rmax_evd',
    'fj_rmax_erlang',
    # Expected maximum functions
    'fj_xmax_exp',
    'fj_xmax_2',
    'fj_xmax_erlang',
    'fj_xmax_hyperexp',
    'fj_xmax_approx',
    'fj_xmax_emma',
    'fj_xmax_normal',
    'fj_xmax_pareto',
    # Order statistics and bounds
    'fj_order_stat',
    'fj_gk_bound',
    'fj_char_max',
    'fj_quantile',
    'fj_sm_tput',
]
