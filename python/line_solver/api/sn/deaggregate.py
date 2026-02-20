"""
Chain result disaggregation functions.

Native Python implementation of chain-to-class performance metric
disaggregation for product-form queueing network analysis.

Ported from MATLAB implementation.
"""

import numpy as np
from dataclasses import dataclass
from typing import Optional

from .network_struct import NetworkStruct, NodeType


@dataclass
class SnDeaggregateResult:
    """
    Result of sn_deaggregate_chain_results calculation.

    Attributes:
        Q: (M, K) Class-level queue lengths
        U: (M, K) Class-level utilizations
        R: (M, K) Class-level response times
        T: (M, K) Class-level throughputs
        C: (1, K) Class-level system response times
        X: (1, K) Class-level system throughputs
    """
    Q: np.ndarray
    U: np.ndarray
    R: np.ndarray
    T: np.ndarray
    C: np.ndarray
    X: np.ndarray


def sn_deaggregate_chain_results(
    sn: NetworkStruct,
    Lchain: np.ndarray,
    ST: Optional[np.ndarray],
    STchain: np.ndarray,
    Vchain: np.ndarray,
    alpha: np.ndarray,
    Qchain: Optional[np.ndarray],
    Uchain: Optional[np.ndarray],
    Rchain: np.ndarray,
    Tchain: np.ndarray,
    Cchain: Optional[np.ndarray],
    Xchain: np.ndarray
) -> SnDeaggregateResult:
    """
    Calculate class-based performance metrics from chain-level performance measures.

    This function disaggregates chain-level performance metrics (queue lengths,
    utilizations, response times, throughputs) to class-level metrics using
    the aggregation factors (alpha).

    Args:
        sn: NetworkStruct object for the queueing network model
        Lchain: (M, C) Service demands per chain
        ST: (M, K) Mean service times per class (optional, computed from rates if None)
        STchain: (M, C) Mean service times per chain
        Vchain: (M, C) Mean visits per chain
        alpha: (M, K) Class aggregation coefficients
        Qchain: (M, C) Mean queue-lengths per chain (optional)
        Uchain: (M, C) Mean utilization per chain (optional)
        Rchain: (M, C) Mean response time per chain
        Tchain: (M, C) Mean throughput per chain
        Cchain: (1, C) Mean system response time per chain (optional, not yet supported)
        Xchain: (1, C) Mean system throughput per chain

    Returns:
        SnDeaggregateResult containing class-level performance metrics:
            - Q: (M, K) queue lengths
            - U: (M, K) utilizations
            - R: (M, K) response times
            - T: (M, K) throughputs
            - C: (1, K) system response times
            - X: (1, K) system throughputs

    Raises:
        RuntimeError: If Cchain is provided (not yet supported)
    """
    # Compute ST from rates if not provided
    if ST is None or (hasattr(ST, 'size') and ST.size == 0):
        rates = sn.rates
        if rates is not None and rates.size > 0:
            with np.errstate(divide='ignore', invalid='ignore'):
                ST = np.where(rates != 0, 1.0 / rates, 0.0)
            ST = np.where(np.isnan(ST), 0.0, ST)
        else:
            ST = np.zeros((sn.nstations, sn.nclasses))

    # Cchain is not yet supported
    if Cchain is not None and (hasattr(Cchain, 'size') and Cchain.size > 0):
        raise RuntimeError("Cchain input to sn_deaggregate_chain_results not yet supported")

    M = sn.nstations
    K = sn.nclasses
    N = sn.njobs.flatten()

    # Initialize output matrices
    X = np.zeros((1, K))
    U = np.zeros((M, K))
    Q = np.zeros((M, K))
    T = np.zeros((M, K))
    R = np.zeros((M, K))
    C = np.zeros((1, K))

    # Find sink index for open class throughput calculation
    Vsink = _compute_vsink(sn)

    # Get nservers
    nservers = sn.nservers
    if nservers is None or len(nservers) == 0:
        nservers = np.ones((M, 1))
    nservers = nservers.flatten()

    # Identify Source and Sink stations (they should have U=0, Q=0, R=0)
    # Source nodes only generate arrivals, Sink nodes only consume departures
    # But throughput T should still be computed for Source (T = arrival rate)
    zero_util_stations = set()
    source_stations = set()
    if sn.nodetype is not None and len(sn.nodetype) > 0:
        stationToNode = sn.stationToNode if hasattr(sn, 'stationToNode') else np.arange(M)
        for i in range(M):
            node_idx = int(stationToNode[i]) if i < len(stationToNode) else i
            if node_idx < len(sn.nodetype):
                if sn.nodetype[node_idx] == NodeType.SOURCE:
                    zero_util_stations.add(i)
                    source_stations.add(i)
                elif sn.nodetype[node_idx] == NodeType.SINK:
                    zero_util_stations.add(i)

    # Get refstat
    refstat = sn.refstat
    if refstat is None or len(refstat) == 0:
        refstat = np.zeros((K, 1), dtype=int)
    refstat = refstat.flatten().astype(int)

    # Ensure Xchain is 2D for indexing
    if Xchain.ndim == 1:
        Xchain = Xchain.reshape(1, -1)

    # Process each chain
    for c in range(sn.nchains):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)

        # Check if chain is open
        is_open_chain = np.isinf(np.sum(N[inchain]))

        for k in inchain:
            if k >= K:
                continue

            # Calculate system throughput X
            if is_open_chain:
                X[0, k] = Xchain[0, c] * Vsink[k] if k < len(Vsink) else 0.0
            else:
                # For closed chains, all classes in the chain have the same system throughput
                # The formula X = Xchain * alpha[refstat, k] fails for class switching
                # where a class doesn't visit its reference station (alpha = 0)
                # In such cases, use Xchain directly since all classes share throughput
                ref_station = refstat[k] if k < len(refstat) else 0
                alpha_at_ref = alpha[ref_station, k] if ref_station < alpha.shape[0] and k < alpha.shape[1] else 0.0
                if alpha_at_ref > 0:
                    # Normalize by sum of alpha at reference station for all classes in chain
                    # This prevents double-counting when multiple entries share the same ref task
                    # (e.g., T1 calling both E2 and E3 in lqn_twotasks)
                    alpha_sum = sum(
                        alpha[ref_station, kk] if kk < alpha.shape[1] else 0.0
                        for kk in inchain
                    )
                    if alpha_sum > 0:
                        X[0, k] = Xchain[0, c] * alpha_at_ref / alpha_sum
                    else:
                        X[0, k] = Xchain[0, c] * alpha_at_ref
                else:
                    # Class doesn't visit its reference station (class switching case)
                    # Use system throughput from chain directly
                    X[0, k] = Xchain[0, c]

            for i in range(M):
                # Reference station for this class
                ref_station = refstat[k] if k < len(refstat) else 0

                # Calculate utilization U
                # Source and Sink stations always have U=0
                if i in zero_util_stations:
                    U[i, k] = 0.0
                else:
                    Vchain_ref = Vchain[ref_station, c] if Vchain[ref_station, c] != 0 else 1.0

                    if Uchain is None or (hasattr(Uchain, 'size') and Uchain.size == 0):
                        # Compute U from ST and Xchain
                        if np.isinf(nservers[i]):
                            U[i, k] = ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain_ref) * alpha[i, k]
                        else:
                            U[i, k] = ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain_ref) * alpha[i, k] / nservers[i]
                    else:
                        if np.isinf(nservers[i]):
                            U[i, k] = ST[i, k] * (Xchain[0, c] * Vchain[i, c] / Vchain_ref) * alpha[i, k]
                        else:
                            U[i, k] = Uchain[i, c] * alpha[i, k]

                # Calculate Q, T, R
                # Source and Sink stations always have Q=0, R=0 (no queueing)
                # But Source should have T = arrival rate (computed from Tchain)
                if i in zero_util_stations:
                    Q[i, k] = 0.0
                    R[i, k] = 0.0
                    # For Source stations, compute throughput from Tchain
                    if i in source_stations and Lchain[i, c] > 0:
                        T[i, k] = Tchain[i, c] * alpha[i, k]
                    else:
                        T[i, k] = 0.0
                elif Lchain[i, c] > 0:
                    # MATLAB: R(i,k) = Q(i,k) / T(i,k) computed unconditionally.
                    # When alpha(i,k) is zero, both Q and T are zero, giving NaN
                    # which is cleaned up at the end. In Python, alpha values may
                    # have numerical noise (e.g. 1e-18) instead of exact zero,
                    # causing tiny/tiny to produce finite but incorrect R values.
                    # Guard: skip computation when alpha is effectively zero.
                    if alpha[i, k] < 1e-14:
                        Q[i, k] = 0.0
                        T[i, k] = 0.0
                        R[i, k] = 0.0
                    else:
                        if Qchain is not None and (hasattr(Qchain, 'size') and Qchain.size > 0):
                            Q[i, k] = Qchain[i, c] * alpha[i, k]
                        else:
                            STchain_val = STchain[i, c] if STchain[i, c] != 0 else 1.0
                            Q[i, k] = (Rchain[i, c] * ST[i, k] / STchain_val *
                                       Xchain[0, c] * Vchain[i, c] / Vchain_ref * alpha[i, k])

                        T[i, k] = Tchain[i, c] * alpha[i, k]

                        if T[i, k] != 0:
                            R[i, k] = Q[i, k] / T[i, k]
                        else:
                            R[i, k] = 0.0
                else:
                    T[i, k] = 0.0
                    R[i, k] = 0.0
                    Q[i, k] = 0.0

            # Calculate system response time C = N / X (Little's law)
            if X[0, k] != 0:
                C[0, k] = N[k] / X[0, k]
            else:
                C[0, k] = 0.0

    # Final cleanup: take absolute values and remove non-finite values
    Q = np.abs(Q)
    R = np.abs(R)
    X = np.abs(X)
    U = np.abs(U)
    T = np.abs(T)
    C = np.abs(C)

    Q = np.where(~np.isfinite(Q), 0.0, Q)
    R = np.where(~np.isfinite(R), 0.0, R)
    X = np.where(~np.isfinite(X), 0.0, X)
    U = np.where(~np.isfinite(U), 0.0, U)
    T = np.where(~np.isfinite(T), 0.0, T)
    C = np.where(~np.isfinite(C), 0.0, C)

    return SnDeaggregateResult(Q=Q, U=U, R=R, T=T, C=C, X=X)


def _compute_vsink(sn: NetworkStruct) -> np.ndarray:
    """
    Compute visit ratios to sink node across all chains.

    Args:
        sn: NetworkStruct object

    Returns:
        (K,) array of visit ratios to sink for each class
    """
    K = sn.nclasses

    if sn.nodevisits is None or len(sn.nodevisits) == 0:
        return np.zeros(K)

    # Find sink node index
    sink_idx = -1
    if sn.nodetype is not None and len(sn.nodetype) > 0:
        for i, ntype in enumerate(sn.nodetype):
            if ntype == NodeType.SINK:
                sink_idx = i
                break

    if sink_idx < 0:
        return np.zeros(K)

    # Sum nodevisits across all chains
    first_key = next(iter(sn.nodevisits), None)
    if first_key is None:
        return np.zeros(K)

    first_visits = sn.nodevisits[first_key]
    if first_visits is None:
        return np.zeros(K)

    Vsum = np.zeros_like(first_visits)
    for c, visits in sn.nodevisits.items():
        if visits is not None:
            Vsum = Vsum + visits

    # Extract sink row
    if sink_idx < Vsum.shape[0]:
        Vsink = Vsum[sink_idx, :]
    else:
        Vsink = np.zeros(K)

    return Vsink.flatten()
