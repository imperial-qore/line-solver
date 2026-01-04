"""
Metrics extractor for FLD solver results.

Extracts QN, UN, RN, TN performance metrics from ODE state vectors
and handler results.
"""

import numpy as np
from typing import Dict, Tuple, Optional
from ....api.sn import NetworkStruct


def extract_metrics_from_handler_result(handler_result, sn: NetworkStruct) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Extract standard metrics from FLD handler result.

    Args:
        handler_result: SolverFLDReturn from api/solvers/fld/handler.py
        sn: Compiled NetworkStruct

    Returns:
        Tuple of (QN, UN, RN, TN, CN, XN)
    """
    M = sn.nstations
    K = sn.nclasses

    # Extract from handler result
    QN = handler_result.Q if handler_result.Q is not None else np.zeros((M, K))
    UN = handler_result.U if handler_result.U is not None else np.zeros((M, K))
    TN = handler_result.T if handler_result.T is not None else np.zeros((M, K))
    RN = handler_result.R if handler_result.R is not None else np.zeros((M, K))
    CN = handler_result.C if handler_result.C is not None else np.zeros((1, K))
    XN = handler_result.X if handler_result.X is not None else np.zeros((1, K))

    return QN, UN, RN, TN, CN, XN


def extract_transient_metrics(
    t_vec: np.ndarray,
    x_vec: np.ndarray,
    SQC: np.ndarray,
    SUC: np.ndarray,
    STC: np.ndarray,
    Qa: np.ndarray,
    Sa: np.ndarray,
    SQ: np.ndarray,
    M: int,
    K: int
) -> Tuple[Dict[Tuple[int, int], np.ndarray], Dict[Tuple[int, int], np.ndarray], Dict[Tuple[int, int], np.ndarray]]:
    """
    Extract transient queue lengths, utilizations, and throughputs.

    Args:
        t_vec: Time points
        x_vec: ODE state vectors over time (n_times x n_states)
        SQC: Queue length selection matrix
        SUC: Utilization selection matrix
        STC: Throughput selection matrix
        Qa: Station-class mapping
        Sa: Service rates
        SQ: State aggregation matrix
        M: Number of stations
        K: Number of classes

    Returns:
        Tuple of (QNt, UNt, TNt) dictionaries mapping (station, class) to time series
    """
    QNt = {}
    UNt = {}
    TNt = {}

    # For each time point
    for t_idx, x_final in enumerate(x_vec):
        # Queue lengths
        if SQC.shape[1] == len(x_final):
            QN_flat = SQC @ x_final
            for i in range(M):
                for r in range(K):
                    idx = i * K + r
                    if idx < len(QN_flat):
                        key = (i, r)
                        if key not in QNt:
                            QNt[key] = []
                        QNt[key].append(QN_flat[idx])

        # Utilizations
        sum_x_Qa = SQ @ x_final + 1e-14
        theta = x_final.copy()
        for phase in range(len(x_final)):
            station = int(Qa[0, phase]) if phase < Qa.shape[1] else 0
            theta[phase] = x_final[phase] / sum_x_Qa[phase] * min(Sa[phase], sum_x_Qa[phase])

        if SUC.shape[1] == len(theta):
            UN_flat = SUC @ theta
            for i in range(M):
                for r in range(K):
                    idx = i * K + r
                    if idx < len(UN_flat):
                        key = (i, r)
                        if key not in UNt:
                            UNt[key] = []
                        UNt[key].append(UN_flat[idx])

        # Throughputs
        if STC.shape[1] == len(theta):
            TN_flat = STC @ theta
            for i in range(M):
                for r in range(K):
                    idx = i * K + r
                    if idx < len(TN_flat):
                        key = (i, r)
                        if key not in TNt:
                            TNt[key] = []
                        TNt[key].append(TN_flat[idx])

    # Convert lists to numpy arrays
    for key in QNt:
        QNt[key] = np.array(QNt[key])
    for key in UNt:
        UNt[key] = np.array(UNt[key])
    for key in TNt:
        TNt[key] = np.array(TNt[key])

    return QNt, UNt, TNt


def compute_response_times(QN: np.ndarray, TN: np.ndarray) -> np.ndarray:
    """
    Compute response times from queue lengths and throughputs using Little's Law.

    R = Q / T (where T is throughput per station, class)

    Args:
        QN: Queue lengths (M x K)
        TN: Throughputs (M x K)

    Returns:
        RN: Response times (M x K)
    """
    RN = np.zeros_like(QN)
    with np.errstate(divide='ignore', invalid='ignore'):
        RN = np.divide(QN, TN)
        RN[~np.isfinite(RN)] = 0.0
    return RN


def compute_cycle_times(RN: np.ndarray) -> np.ndarray:
    """
    Compute cycle times as sum of response times across all stations.

    C = sum_i R[i, :]

    Args:
        RN: Response times (M x K)

    Returns:
        CN: Cycle times (1 x K)
    """
    CN = np.sum(RN, axis=0, keepdims=True)
    return CN


def compute_system_throughput(TN: np.ndarray) -> np.ndarray:
    """
    Compute system throughput as mean throughput across stations.

    X = mean(T)

    Args:
        TN: Throughputs (M x K)

    Returns:
        XN: System throughputs (1 x K)
    """
    XN = np.mean(TN, axis=0, keepdims=True)
    return XN
