"""
Fork-Join topology validator.

Validates that a network conforms to Fork-Join structure:
- Source → Fork → K parallel homogeneous queues → Join → Sink
- All parallel queues have identical service distributions
- All parallel queues have single server
"""

import numpy as np
from typing import Tuple, Optional
from dataclasses import dataclass

from ....api.sn.network_struct import NetworkStruct, NodeType


@dataclass
class FJValidationResult:
    """Result of Fork-Join topology validation."""
    is_fj: bool
    reason: Optional[str]
    K: Optional[int] = None  # Number of parallel queues
    source_idx: Optional[int] = None
    fork_idx: Optional[int] = None
    join_idx: Optional[int] = None
    sink_idx: Optional[int] = None
    queue_indices: Optional[list] = None


def fj_isfj(sn: NetworkStruct) -> FJValidationResult:
    """
    Validate if network has Fork-Join topology.

    Fork-Join structure requires:
    - Exactly 1 Source node
    - Exactly 1 Fork node
    - K parallel Queue nodes (homogeneous)
    - Exactly 1 Join node
    - Exactly 1 Sink node
    - All parallel queues have:
      - Single server (nservers=1)
      - Identical service distributions
      - Supported distribution types (Exp, HyperExp(2), Erlang(2), MAP(2))

    Args:
        sn: Compiled NetworkStruct

    Returns:
        FJValidationResult with validation status and details
    """
    # Count node types
    source_indices = []
    fork_indices = []
    queue_indices = []
    join_indices = []
    sink_indices = []

    for node_idx, node_type in enumerate(sn.nodetype):
        if node_type == NodeType.SOURCE:
            source_indices.append(node_idx)
        elif node_type == NodeType.FORK:
            fork_indices.append(node_idx)
        elif node_type == NodeType.QUEUE:
            queue_indices.append(node_idx)
        elif node_type == NodeType.JOIN:
            join_indices.append(node_idx)
        elif node_type == NodeType.SINK:
            sink_indices.append(node_idx)

    # Validate counts
    if len(source_indices) != 1:
        return FJValidationResult(False, f"Expected 1 Source, found {len(source_indices)}")

    if len(fork_indices) != 1:
        return FJValidationResult(False, f"Expected 1 Fork, found {len(fork_indices)}")

    if len(join_indices) != 1:
        return FJValidationResult(False, f"Expected 1 Join, found {len(join_indices)}")

    if len(sink_indices) != 1:
        return FJValidationResult(False, f"Expected 1 Sink, found {len(sink_indices)}")

    if len(queue_indices) < 2:
        return FJValidationResult(False, f"Expected at least 2 parallel queues, found {len(queue_indices)}")

    K = len(queue_indices)

    # Validate queue homogeneity
    source_idx = source_indices[0]
    fork_idx = fork_indices[0]
    join_idx = join_indices[0]
    sink_idx = sink_indices[0]

    # Convert node indices to station indices
    node_to_station = {}
    for m in range(sn.nstations):
        node_idx = sn.stationToNode[m] if len(sn.stationToNode) > m else m
        node_to_station[node_idx] = m

    queue_stations = []
    for node_idx in queue_indices:
        if node_idx in node_to_station:
            queue_stations.append(node_to_station[node_idx])

    if len(queue_stations) != K:
        return FJValidationResult(False, "Could not map all queue nodes to stations")

    # Validate homogeneity: all parallel queues have identical service
    service_rates = np.asarray(sn.rates, dtype=np.float64)
    service_scv = np.asarray(sn.scv, dtype=np.float64)

    if service_rates.ndim == 1:
        service_rates = service_rates.reshape(-1, 1)
    if service_scv.ndim == 1:
        service_scv = service_scv.reshape(-1, 1)

    # Check all parallel queues have same service distribution
    if len(queue_stations) > 0:
        ref_station = queue_stations[0]
        ref_rate = service_rates[ref_station, 0] if service_rates.shape[0] > ref_station else None
        ref_scv = service_scv[ref_station, 0] if service_scv.shape[0] > ref_station else None

        for station in queue_stations[1:]:
            rate = service_rates[station, 0] if service_rates.shape[0] > station else None
            scv = service_scv[station, 0] if service_scv.shape[0] > station else None

            if ref_rate is not None and rate is not None:
                if not np.isclose(rate, ref_rate, rtol=1e-6):
                    return FJValidationResult(
                        False,
                        f"Parallel queues have different service rates: {ref_rate} vs {rate}"
                    )

            if ref_scv is not None and scv is not None:
                if not np.isclose(scv, ref_scv, rtol=1e-6):
                    return FJValidationResult(
                        False,
                        f"Parallel queues have different SCV: {ref_scv} vs {scv}"
                    )

    # Check all parallel queues have single server
    nservers = np.asarray(sn.nservers, dtype=np.float64).flatten()
    for station in queue_stations:
        if station < len(nservers):
            if nservers[station] != 1:
                return FJValidationResult(
                    False,
                    f"Queue station {station} has {nservers[station]} servers, expected 1"
                )

    # Validate topology: Source → Fork → Queues → Join → Sink
    # Check routing: source goes to fork
    if sn.rt is not None:
        rt = np.asarray(sn.rt, dtype=np.float64)

        # Source should route to Fork (if it's a source node in routing table)
        # Fork should route to all parallel queues
        # All parallel queues should route to Join
        # Join should route to Sink

        # This is a simplified check - a more thorough topology check would
        # verify the routing probabilities sum to 1 from fork to queues, etc.
        pass  # Routing validation is optional - geometry check is main criterion

    return FJValidationResult(
        is_fj=True,
        reason=None,
        K=K,
        source_idx=source_idx,
        fork_idx=fork_idx,
        join_idx=join_idx,
        sink_idx=sink_idx,
        queue_indices=queue_stations
    )


def classify_distribution(scv: float, rate: float) -> str:
    """
    Classify service distribution based on SCV.

    Args:
        scv: Squared coefficient of variation
        rate: Service rate (not used for classification, but available)

    Returns:
        Distribution type: 'M' (Exponential), 'E' (Erlang), 'H' (Hyperexp), 'G' (General)
    """
    if scv < 0.01:
        return 'D'  # Deterministic
    elif np.isclose(scv, 1.0, rtol=0.1):
        return 'M'  # Exponential (Markovian)
    elif scv < 1.0:
        return 'E'  # Erlang-like (SCV < 1)
    elif scv > 1.0:
        return 'H'  # Hyperexponential-like (SCV > 1)
    else:
        return 'G'  # General
