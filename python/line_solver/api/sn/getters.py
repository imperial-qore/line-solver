"""
SN Getter Functions for Parameter Extraction.

Native Python implementations for extracting parameters from
network structures including arrival rates, throughputs, and
product-form chain parameters.

Key functions:
    sn_get_arvr_from_tput: Compute arrival rates from throughputs
    sn_get_node_arvr_from_tput: Compute node arrival rates from throughputs
    sn_get_node_tput_from_tput: Compute node throughputs from station throughputs
    sn_get_product_form_chain_params: Extract chain-aggregated parameters

References:
    Original MATLAB: matlab/src/api/sn/sn_get_*.m
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from .network_struct import NetworkStruct, NodeType


@dataclass
class ChainParams:
    """Chain-aggregated product-form parameters."""
    lambda_vec: np.ndarray  # Chain arrival rates
    D: np.ndarray  # Chain service demands at queuing stations
    N: np.ndarray  # Chain populations
    Z: np.ndarray  # Chain think times
    mu: np.ndarray  # Load-dependent service capacity scaling
    S: np.ndarray  # Number of servers at queuing stations
    V: np.ndarray  # Chain visit ratios


def sn_get_arvr_from_tput(sn: NetworkStruct, TN: np.ndarray,
                          TH: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute average arrival rates at stations from throughputs.

    Calculates the average arrival rate at each station in steady-state
    from the station throughputs and routing matrix.

    Args:
        sn: Network structure
        TN: Average throughputs at stations (M x R)
        TH: Throughput handles (optional)

    Returns:
        AN: Average arrival rates at stations (M x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_arvr_from_tput.m
    """
    M = sn.nstations
    R = sn.nclasses

    if TH is None or TN is None or len(TN) == 0:
        return np.array([])

    TN = np.atleast_2d(np.asarray(TN, dtype=float))
    AN = np.zeros((M, R))

    # Build mapping from stateful nodes to their position in rt matrix
    stateful_nodes = np.where(sn.isstateful)[0]
    n_stateful = len(stateful_nodes)

    # Build throughput vector for all stateful nodes
    TN_stateful = np.zeros((n_stateful, R))
    for sf, ind in enumerate(stateful_nodes):
        ist = sn.nodeToStation[ind]
        if ist >= 0:
            # This stateful node is a station
            TN_stateful[sf, :] = TN[ist, :]

    # Compute arrival rates using stateful node throughputs and rt matrix
    for ist in range(M):
        ind_ist = sn.stationToNode[ist]
        if sn.nodetype[ind_ist] == NodeType.Source:
            AN[ist, :] = 0
        else:
            # Find position of this station in stateful node list
            sf_ist_arr = np.where(stateful_nodes == ind_ist)[0]
            if len(sf_ist_arr) == 0:
                continue
            sf_ist = sf_ist_arr[0]

            for sf_jst in range(n_stateful):
                for k in range(R):
                    for r in range(R):
                        from_idx = sf_jst * R + r
                        to_idx = sf_ist * R + k
                        if from_idx < sn.rt.shape[0] and to_idx < sn.rt.shape[1]:
                            AN[ist, k] += TN_stateful[sf_jst, r] * sn.rt[from_idx, to_idx]

    return AN


def sn_get_node_arvr_from_tput(sn: NetworkStruct, TN: np.ndarray,
                                TH: Optional[np.ndarray] = None,
                                AN: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute node arrival rates from station throughputs.

    Args:
        sn: Network structure
        TN: Station throughputs (M x R)
        TH: Throughput handles (optional)
        AN: Station arrival rates (optional, computed if not provided)

    Returns:
        ANn: Node arrival rates (I x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_node_arvr_from_tput.m
    """
    I = sn.nnodes
    M = sn.nstations
    R = sn.nclasses

    if AN is None:
        AN = sn_get_arvr_from_tput(sn, TN, TH)

    ANn = np.zeros((I, R))

    for ind in range(I):
        ist = sn.nodeToStation[ind]
        if ist >= 0:
            ANn[ind, :] = AN[ist, :]

    return ANn


def sn_get_node_tput_from_tput(sn: NetworkStruct, TN: np.ndarray) -> np.ndarray:
    """
    Compute node throughputs from station throughputs.

    Args:
        sn: Network structure
        TN: Station throughputs (M x R)

    Returns:
        TNn: Node throughputs (I x R)

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_node_tput_from_tput.m
    """
    I = sn.nnodes
    R = sn.nclasses

    TN = np.atleast_2d(np.asarray(TN, dtype=float))
    TNn = np.zeros((I, R))

    for ind in range(I):
        ist = sn.nodeToStation[ind]
        if ist >= 0:
            TNn[ind, :] = TN[ist, :]

    return TNn


def sn_get_product_form_chain_params(sn: NetworkStruct) -> ChainParams:
    """
    Extract product-form parameters aggregated by chain.

    Extracts parameters from a network structure and aggregates them
    by chain for use in product-form analysis methods.

    Args:
        sn: Network structure

    Returns:
        ChainParams with lambda_vec, D, N, Z, mu, S, V

    References:
        Original MATLAB: matlab/src/api/sn/sn_get_product_form_chain_params.m
    """
    from .transforms import sn_get_product_form_params
    from .demands import sn_get_demands_chain

    # Get base product-form parameters
    params = sn_get_product_form_params(sn)

    # Get chain-aggregated demands
    demands = sn_get_demands_chain(sn)

    # Find queue and delay indices
    queue_indices = np.where(sn.nodetype == NodeType.Queue)[0]
    delay_indices = np.where(sn.nodetype == NodeType.Delay)[0]

    # Initialize chain parameters
    nchains = sn.nchains
    lambda_chains = np.zeros(nchains)

    for c in range(nchains):
        chain_classes = sn.inchain[c]
        lambda_chains[c] = np.nansum(params.lambda_vec[chain_classes])

    # Extract demands at queue and delay stations
    D_chains = demands.Dchain[sn.nodeToStation[queue_indices], :]
    Z_chains = demands.Dchain[sn.nodeToStation[delay_indices], :] if len(delay_indices) > 0 else np.zeros((0, nchains))

    # Number of servers at queuing stations
    S = sn.nservers[sn.nodeToStation[queue_indices]]

    # Visit ratios
    V = demands.Vchain.copy()
    ignore_indices = np.where((sn.nodetype == NodeType.Source) | (sn.nodetype == NodeType.Join))[0]
    if len(ignore_indices) > 0:
        keep_stations = [sn.nodeToStation[i] for i in range(sn.nnodes) if i not in ignore_indices]
        V = V[keep_stations, :]

    if len(Z_chains) == 0:
        Z_chains = np.zeros((0, nchains))

    return ChainParams(
        lambda_vec=lambda_chains,
        D=D_chains,
        N=demands.Nchain,
        Z=np.sum(Z_chains, axis=0) if Z_chains.size > 0 else np.zeros(nchains),
        mu=params.mu,
        S=S,
        V=V
    )


def sn_set_routing_prob(sn: NetworkStruct, from_stateful: int, from_class: int,
                        to_stateful: int, to_class: int, prob: float,
                        auto_refresh: bool = False) -> NetworkStruct:
    """
    Set a routing probability between two stateful node-class pairs.

    Updates a single entry in the rt matrix.

    Args:
        sn: Network structure
        from_stateful: Source stateful node index (0-based)
        from_class: Source class index (0-based)
        to_stateful: Destination stateful node index (0-based)
        to_class: Destination class index (0-based)
        prob: Routing probability [0, 1]
        auto_refresh: If True, refresh visit ratios (default False)

    Returns:
        Modified network structure

    References:
        Original MATLAB: matlab/src/api/sn/sn_set_routing_prob.m
    """
    K = sn.nclasses

    # Calculate indices in rt matrix
    from_idx = from_stateful * K + from_class
    to_idx = to_stateful * K + to_class

    # Update rt matrix
    sn.rt[from_idx, to_idx] = prob

    # Auto-refresh visit ratios if requested
    if auto_refresh:
        from .transforms import sn_refresh_visits
        sn_refresh_visits(sn)

    return sn


__all__ = [
    'ChainParams',
    'sn_get_arvr_from_tput',
    'sn_get_node_arvr_from_tput',
    'sn_get_node_tput_from_tput',
    'sn_get_product_form_chain_params',
    'sn_set_routing_prob',
]
