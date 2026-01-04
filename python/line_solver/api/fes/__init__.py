"""
Flow-Equivalent Server (FES) Analysis.

Native Python implementations for Flow-Equivalent Server analysis,
which allows decomposition of complex queueing networks by replacing
subnetworks with equivalent load-dependent servers.

Key functions:
    fes_validate: Validate FES topology
    fes_compute_throughputs: Compute throughputs using FES decomposition
    fes_build_isolated: Build isolated subsystem for FES analysis

References:
    Original MATLAB: matlab/src/api/fes/fes_*.m
    Chandy, Herzog, Woo, "Parametric analysis of queueing networks", 1975
"""

import numpy as np
from typing import Tuple, Optional, Dict, Any
from dataclasses import dataclass


@dataclass
class FesResult:
    """Result of FES throughput computation."""
    throughput: np.ndarray  # Throughput per chain
    queue_length: np.ndarray  # Queue lengths
    response_time: np.ndarray  # Response times


def fes_validate(sn) -> bool:
    """
    Validate that a network is suitable for FES analysis.

    Checks that the network has the proper structure for Flow-Equivalent
    Server decomposition (no class switching, proper topology, etc.).

    Args:
        sn: Network structure

    Returns:
        True if network is valid for FES analysis

    References:
        Original MATLAB: matlab/src/api/fes/fes_validate.m
    """
    # Check for class switching (not supported)
    if hasattr(sn, 'classswitch') and sn.classswitch is not None:
        if np.any(sn.classswitch):
            return False

    # Check for fork-join (not supported in basic FES)
    if hasattr(sn, 'fj') and sn.fj is not None:
        if np.any(sn.fj):
            return False

    # Check for closed network (required for FES)
    if hasattr(sn, 'nclosedjobs'):
        if sn.nclosedjobs == 0:
            return False

    return True


def fes_compute_throughputs(sn, mu_fes: np.ndarray,
                             population: Optional[np.ndarray] = None
                             ) -> FesResult:
    """
    Compute throughputs using FES decomposition.

    Uses the Flow-Equivalent Server method to compute throughputs
    for a closed queueing network decomposed into subsystems.

    Args:
        sn: Network structure
        mu_fes: Load-dependent service rates of the FES (chain x population)
        population: Population vector (default: from sn)

    Returns:
        FesResult with throughput, queue_length, response_time

    References:
        Original MATLAB: matlab/src/api/fes/fes_compute_throughputs.m
    """
    if population is None:
        if hasattr(sn, 'njobs'):
            population = sn.njobs
        else:
            raise ValueError("Population must be provided")

    population = np.asarray(population, dtype=int).flatten()
    mu_fes = np.atleast_2d(np.asarray(mu_fes, dtype=float))

    nchains = len(population)
    N_total = int(np.sum(population))

    # For single chain, use MVA with load-dependent server
    if nchains == 1:
        N = population[0]

        # MVA with load-dependent rates
        L = np.zeros(N + 1)  # Queue length
        W = np.zeros(N + 1)  # Response time
        X = np.zeros(N + 1)  # Throughput

        for n in range(1, N + 1):
            # Mean response time using arrival theorem
            if mu_fes.shape[1] >= n:
                mu_n = mu_fes[0, n - 1]
            else:
                mu_n = mu_fes[0, -1]

            if mu_n > 0:
                W[n] = (1 + L[n - 1]) / mu_n
            else:
                W[n] = np.inf

            # Throughput (Little's law)
            if W[n] > 0:
                X[n] = n / W[n]
            else:
                X[n] = 0

            # Queue length
            L[n] = X[n] * W[n]

        return FesResult(
            throughput=np.array([X[N]]),
            queue_length=np.array([L[N]]),
            response_time=np.array([W[N]])
        )

    else:
        # Multi-chain case - use MVA with multiple chains
        # This is a simplified implementation
        X = np.zeros(nchains)
        Q = np.zeros(nchains)
        W = np.zeros(nchains)

        for c in range(nchains):
            N_c = population[c]
            if N_c > 0 and mu_fes.shape[0] > c:
                if mu_fes.shape[1] >= N_c:
                    mu_n = mu_fes[c, N_c - 1]
                else:
                    mu_n = mu_fes[c, -1]

                if mu_n > 0:
                    W[c] = 1.0 / mu_n  # Simplified
                    X[c] = N_c / W[c] if W[c] > 0 else 0
                    Q[c] = X[c] * W[c]

        return FesResult(
            throughput=X,
            queue_length=Q,
            response_time=W
        )


def fes_build_isolated(sn, subsystem_nodes: np.ndarray) -> Dict[str, Any]:
    """
    Build an isolated subsystem for FES analysis.

    Extracts a subsystem from the network and prepares it for
    isolated analysis as part of FES decomposition.

    Args:
        sn: Network structure
        subsystem_nodes: Node indices belonging to the subsystem

    Returns:
        Dictionary with subsystem parameters:
            - D: Service demands
            - Z: Think times
            - N: Population
            - nservers: Number of servers

    References:
        Original MATLAB: matlab/src/api/fes/fes_build_isolated.m
    """
    subsystem_nodes = np.asarray(subsystem_nodes, dtype=int).flatten()

    # Extract stations in subsystem
    stations = []
    for node in subsystem_nodes:
        if hasattr(sn, 'nodeToStation') and sn.nodeToStation[node] >= 0:
            stations.append(sn.nodeToStation[node])

    stations = np.array(stations, dtype=int)

    if len(stations) == 0:
        return {'D': np.array([]), 'Z': np.array([]), 'N': np.array([]), 'nservers': np.array([])}

    # Extract service demands
    if hasattr(sn, 'proc') and sn.proc is not None:
        D = np.zeros((len(stations), sn.nclasses))
        for i, ist in enumerate(stations):
            for r in range(sn.nclasses):
                if sn.proc[ist, r] is not None:
                    D[i, r] = sn.proc[ist, r].get('mean', 0)
    else:
        D = np.zeros((len(stations), sn.nclasses))

    # Extract think times (for delay stations)
    Z = np.zeros(sn.nclasses)
    if hasattr(sn, 'nodetype'):
        from .network_struct import NodeType
        for node in subsystem_nodes:
            if sn.nodetype[node] == NodeType.Delay:
                ist = sn.nodeToStation[node]
                if ist >= 0 and hasattr(sn, 'proc'):
                    for r in range(sn.nclasses):
                        if sn.proc[ist, r] is not None:
                            Z[r] += sn.proc[ist, r].get('mean', 0)

    # Population
    N = sn.njobs if hasattr(sn, 'njobs') else np.zeros(sn.nclasses)

    # Number of servers
    nservers = sn.nservers[stations] if hasattr(sn, 'nservers') else np.ones(len(stations))

    return {
        'D': D,
        'Z': Z,
        'N': N,
        'nservers': nservers,
        'stations': stations
    }


__all__ = [
    'FesResult',
    'fes_validate',
    'fes_compute_throughputs',
    'fes_build_isolated',
]
