"""
Chain demand calculation functions.

Native Python implementation of chain-based demand aggregation
for product-form queueing network analysis.

Ported from MATLAB implementation.
"""

import numpy as np
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

from .network_struct import NetworkStruct


@dataclass
class SnGetDemandsResult:
    """
    Result of sn_get_demands_chain calculation.

    Attributes:
        Lchain: (M, C) Chain-level demand matrix
        STchain: (M, C) Chain-level service time matrix
        Vchain: (M, C) Chain-level visit ratio matrix
        alpha: (M, K) Class-to-chain weighting matrix
        Nchain: (1, C) Population per chain
        SCVchain: (M, C) Chain-level squared coefficient of variation
        refstatchain: (C, 1) Reference station per chain
    """
    Lchain: np.ndarray
    STchain: np.ndarray
    Vchain: np.ndarray
    alpha: np.ndarray
    Nchain: np.ndarray
    SCVchain: np.ndarray
    refstatchain: np.ndarray


def sn_get_demands_chain(sn: NetworkStruct) -> SnGetDemandsResult:
    """
    Calculate new queueing network parameters after aggregating classes into chains.

    This function computes chain-level demands, service times, visit ratios,
    and other parameters by aggregating class-level data based on chain membership.

    Args:
        sn: NetworkStruct object for the queueing network model

    Returns:
        SnGetDemandsResult containing chain parameters:
            - Lchain: (M, C) chain-level demand matrix
            - STchain: (M, C) chain-level service time matrix
            - Vchain: (M, C) chain-level visit ratio matrix
            - alpha: (M, K) class-to-chain weighting matrix
            - Nchain: (1, C) population per chain
            - SCVchain: (M, C) chain-level squared coefficient of variation
            - refstatchain: (C, 1) reference station per chain
    """
    M = sn.nstations
    K = sn.nclasses
    C = sn.nchains
    N = sn.njobs.flatten()

    # Get SCV matrix, replacing NaN with 1.0
    scv = sn.scv.copy() if sn.scv is not None else np.ones((M, K))
    scv = np.where(np.isnan(scv), 1.0, scv)

    # Compute service times ST = 1 / rates
    rates = sn.rates.copy() if sn.rates is not None else np.ones((M, K))
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = np.where(rates != 0, 1.0 / rates, 0.0)
    ST = np.where(np.isnan(ST), 0.0, ST)

    # Initialize output matrices
    alpha = np.zeros((M, K))
    Vchain = np.zeros((M, C))

    # Get station to stateful mapping (flatten in case it's 2D from wrapper)
    station_to_stateful = sn.stationToStateful
    if station_to_stateful is None or len(station_to_stateful) == 0:
        # Default: stations are stateful nodes in order
        station_to_stateful = np.arange(M)
    else:
        station_to_stateful = np.asarray(station_to_stateful).flatten()

    # Get refstat and refclass
    refstat = sn.refstat.flatten() if sn.refstat is not None else np.zeros(K, dtype=int)
    refclass = sn.refclass.flatten() if sn.refclass is not None else -np.ones(C, dtype=int)

    # Calculate Vchain and alpha for each chain
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)

        if c < len(refclass) and refclass[c] > -1:
            # Reference class is specified
            for i in range(M):
                visits = sn.visits.get(c)
                if visits is None:
                    continue

                # visits matrix is indexed by station (not stateful node)
                # Shape is (nstations, nclasses)
                if i >= visits.shape[0]:
                    continue

                # Sum visits for all classes in chain at station i
                sum_visits_i = np.sum(visits[i, inchain])

                # Get reference station index
                ref_stat_idx = int(refstat[inchain[0]]) if len(refstat) > inchain[0] else 0
                if ref_stat_idx >= visits.shape[0]:
                    ref_stat_idx = 0

                # Get reference class visits
                ref_class_idx = int(refclass[c])
                ref_visits = visits[ref_stat_idx, ref_class_idx] if ref_class_idx < visits.shape[1] else 1.0

                Vchain[i, c] = sum_visits_i / ref_visits if ref_visits != 0 else 0.0

                # Calculate alpha weights
                if sum_visits_i > 0:
                    for k in inchain:
                        if k < K:
                            alpha[i, k] += visits[i, k] / sum_visits_i
        else:
            # No reference class, use sum of visits in chain
            for i in range(M):
                visits = sn.visits.get(c)
                if visits is None:
                    continue

                # visits matrix is indexed by station (not stateful node)
                if i >= visits.shape[0]:
                    continue

                # Sum visits for classes in chain at station i
                sum_visits_i = np.sum(visits[i, inchain])

                # Get reference station visits
                ref_stat_idx = int(refstat[inchain[0]]) if len(refstat) > inchain[0] else 0
                if ref_stat_idx >= visits.shape[0]:
                    ref_stat_idx = 0
                sum_visits_ref = np.sum(visits[ref_stat_idx, inchain])

                Vchain[i, c] = sum_visits_i / sum_visits_ref if sum_visits_ref != 0 else 0.0

                # Calculate alpha weights
                if sum_visits_i > 0:
                    for k in inchain:
                        if k < K:
                            alpha[i, k] += visits[i, k] / sum_visits_i

    # Clean up Vchain
    Vchain = np.where(np.isinf(Vchain), 0.0, Vchain)
    Vchain = np.where(np.isnan(Vchain), 0.0, Vchain)

    # Normalize Vchain by reference station visits
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)
        ref_stat_idx = int(refstat[inchain[0]]) if len(refstat) > inchain[0] else 0
        vchain_ref = Vchain[ref_stat_idx, c]
        if vchain_ref != 0:
            Vchain[:, c] /= vchain_ref

    # Clean up alpha
    alpha = np.where(np.isinf(alpha), 0.0, alpha)
    alpha = np.where(np.isnan(alpha), 0.0, alpha)
    alpha = np.maximum(alpha, 0.0)  # Ensure non-negative

    # Initialize chain-level matrices
    Lchain = np.zeros((M, C))
    STchain = np.zeros((M, C))
    SCVchain = np.zeros((M, C))
    Nchain = np.zeros((1, C))
    refstatchain = np.zeros((C, 1))

    # Calculate chain-level parameters
    for c in range(C):
        if c not in sn.inchain:
            continue
        inchain = sn.inchain[c].flatten().astype(int)

        # Calculate Nchain and detect open chain
        Nchain_sum = np.sum(N[inchain])
        is_open_chain = np.isinf(Nchain_sum)
        Nchain[0, c] = Nchain_sum

        for i in range(M):
            ref_stat_idx = int(refstat[inchain[0]]) if len(refstat) > inchain[0] else 0

            if is_open_chain and i == ref_stat_idx:
                # For open chains at reference station: STchain = 1 / sum(finite rates)
                rates_inchain = rates[i, inchain]
                finite_rates = rates_inchain[np.isfinite(rates_inchain)]
                sum_rates = np.sum(finite_rates)
                STchain[i, c] = 1.0 / sum_rates if sum_rates != 0 else 0.0
            else:
                # STchain = ST * alpha weighted sum
                STchain[i, c] = np.sum(ST[i, inchain] * alpha[i, inchain])

            # Lchain = Vchain * STchain
            Lchain[i, c] = Vchain[i, c] * STchain[i, c]

            # Calculate SCVchain
            scv_inchain = scv[i, inchain]
            alpha_inchain = alpha[i, inchain]
            finite_mask = np.isfinite(scv_inchain)
            alphachain = np.sum(alpha_inchain[finite_mask])

            if alphachain > 1e-10:
                SCVchain[i, c] = np.sum(scv_inchain * alpha_inchain) / alphachain

        # Set reference station for chain
        refstatchain[c, 0] = ref_stat_idx

        # Verify all classes in chain have same reference station
        for k in inchain[1:]:
            if len(refstat) > k and refstat[k] != refstatchain[c, 0]:
                raise ValueError(f"Classes in chain {c} have different reference stations")

    # Final cleanup
    Lchain = np.where(np.isinf(Lchain), 0.0, Lchain)
    Lchain = np.where(np.isnan(Lchain), 0.0, Lchain)
    STchain = np.where(np.isinf(STchain), 0.0, STchain)
    STchain = np.where(np.isnan(STchain), 0.0, STchain)

    return SnGetDemandsResult(
        Lchain=Lchain,
        STchain=STchain,
        Vchain=Vchain,
        alpha=alpha,
        Nchain=Nchain,
        SCVchain=SCVchain,
        refstatchain=refstatchain
    )
