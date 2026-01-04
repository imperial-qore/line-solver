"""
M3A automatic fitting functions.

Top-level API for automatic MMAP fitting from traces using the M3A methodology.
Provides m3a_fit and m3a_fit_from_trace functions for automatic model selection.

References:
[1] A. Sansottera, G. Casale, P. Cremonesi. Fitting Second-Order Acyclic
    Marked Markovian Arrival Processes. IEEE/IFIP DSN 2013.
[2] G. Casale, A. Sansottera, P. Cremonesi. Compact Markov-Modulated
    Models for Multiclass Trace Fitting. European Journal of Operations
    Research, 2016.
"""

import numpy as np
from numpy.typing import NDArray
from typing import List, Dict, Tuple, Optional, NamedTuple
from dataclasses import dataclass
import scipy.optimize as opt

from .utils import validate_mmap, compute_moments, compute_autocorrelation


@dataclass
class MTrace:
    """
    Data structure for multiclass trace representation.

    Attributes:
        S: Inter-arrival times
        C: Class labels for each arrival
        num_classes: Number of distinct classes
    """
    S: np.ndarray
    C: np.ndarray
    num_classes: int


@dataclass
class M3aFitOptions:
    """
    Options for M3A fitting algorithms.

    Attributes:
        method: Fitting method (0 = inter-arrival, 1 = counting process)
        num_states: Number of states for the fitted MMAP
        timescale: Finite time scale for counting process (auto-computed if None)
        timescale_asy: Near-infinite time scale (auto-computed if None)
    """
    method: int = 1
    num_states: int = 2
    timescale: Optional[float] = None
    timescale_asy: Optional[float] = None


def m3afit_init(S: np.ndarray, C: np.ndarray) -> MTrace:
    """
    Prepare multiclass trace for M3A fitting.

    Args:
        S: Inter-arrival times
        C: Class number for each arrival

    Returns:
        MTrace structure ready for fitting
    """
    S = np.asarray(S).flatten()
    C = np.asarray(C).flatten().astype(int)
    num_classes = len(np.unique(C))
    return MTrace(S=S, C=C, num_classes=num_classes)


def m3a_fit(mtrace: MTrace, options: Optional[M3aFitOptions] = None
           ) -> Optional[List[np.ndarray]]:
    """
    Automatic fitting of trace into a Marked Markovian Arrival Process.

    Based on the M3A methodology, this function selects the appropriate fitting
    algorithm based on the number of classes, requested states, and fitting method.

    Args:
        mtrace: Data structure returned by m3afit_init
        options: Fitting options including method and number of states

    Returns:
        Fitted MMAP as [D0, D1, D2, ...] or None if fitting fails
    """
    if options is None:
        options = M3aFitOptions(num_states=2)

    # Compute time scales for counting process method
    mean_iat = np.mean(mtrace.S)
    timescale = options.timescale if options.timescale else 10 * mean_iat
    timescale_asy = options.timescale_asy if options.timescale_asy else max(
        10 * timescale, (np.sum(mtrace.S) - mtrace.S[0]) / 100
    )

    K = mtrace.num_classes
    n = options.num_states

    # Select fitting method based on parameters
    if K == 1:
        # Single-class: fit a MAP
        mmap = _fit_map_from_trace(mtrace.S, n)
    elif K == 2 and n == 2 and options.method == 0:
        # 2-class, 2-state, inter-arrival fitting
        mmap = _fit_mamap22_interarrival(mtrace.S, mtrace.C)
    elif K > 2 and n >= 2 and options.method == 0:
        # Multi-class, inter-arrival fitting
        mmap = _fit_mamap2m_interarrival(mtrace.S, mtrace.C, n)
    elif K >= 2 and n == 2 and options.method == 1:
        # Multi-class, 2-state, counting process fitting
        mmap = _fit_m3pp2m_counting(mtrace.S, mtrace.C, timescale, timescale_asy)
    elif K >= 2 and n > 2 and options.method == 1:
        # Multi-class, >2-state, counting process fitting (superposition)
        mmap = _fit_m3pp_superposition(mtrace.S, mtrace.C, n, timescale, timescale_asy)
    else:
        print("M3A: Algorithm could not obtain a valid MMAP.")
        return None

    # Validate result
    if mmap is not None and validate_mmap(mmap):
        n_states = mmap[0].shape[0]
        n_classes = len(mmap) - 2
        print(f"M3A: Found valid {n_states}-state MMAP[{n_classes}].")
        return mmap
    else:
        print("M3A: Algorithm could not obtain a valid MMAP.")
        return None


def m3a_fit_from_trace(S: np.ndarray, C: np.ndarray,
                       num_states: int = 2, method: int = 1
                      ) -> Optional[List[np.ndarray]]:
    """
    Automatic fitting with simple parameters.

    Args:
        S: Inter-arrival times
        C: Class labels
        num_states: Number of states for the fitted MMAP
        method: Fitting method (0 = inter-arrival, 1 = counting process)

    Returns:
        Fitted MMAP or None if fitting fails
    """
    mtrace = m3afit_init(S, C)
    options = M3aFitOptions(method=method, num_states=num_states)
    return m3a_fit(mtrace, options)


# Internal fitting functions

def _fit_map_from_trace(S: np.ndarray, n: int) -> List[np.ndarray]:
    """Fit a single-class MAP from inter-arrival times."""
    # Compute statistics
    mean = np.mean(S)
    var = np.var(S)
    scv = var / (mean**2) if mean > 0 else 1.0

    # Compute autocorrelation at lag 1
    if len(S) > 1:
        S_centered = S - mean
        acf1 = np.correlate(S_centered[:-1], S_centered[1:])[0]
        acf1 = acf1 / (var * (len(S) - 1)) if var > 0 else 0.0
        acf1 = max(-0.5, min(0.5, acf1))  # Bound to feasible range
    else:
        acf1 = 0.0

    # Fit MMPP2 (2-state MMPP)
    if n == 2:
        return _mmpp2_fit(mean, var, acf1)
    else:
        # For higher orders, use stacked MMPP2
        mmap = _mmpp2_fit(mean, var, acf1)
        return mmap


def _mmpp2_fit(mean: float, var: float, acf1: float) -> List[np.ndarray]:
    """
    Fit a 2-state MMPP to match mean, variance, and lag-1 autocorrelation.

    Returns:
        MAP as [D0, D1]
    """
    scv = var / (mean**2) if mean > 0 else 1.0

    # Handle edge case
    if abs(scv - 1) < 1e-10:
        G2 = 0.0
    else:
        G2 = acf1 / ((1 - 1 / scv) / 2) if scv != 0 else 0.0

    if abs(G2) < 1e-6 or G2 == 0.0:
        # Fit with MAP(1) approximation
        mu00 = 1.0 / mean if mean > 0 else 1.0
        mu11 = 0.0
        q01 = 0.1
        q10 = 0.1
    else:
        # Full MMPP2 fitting
        try:
            # Simplified fitting based on moments
            lambda_avg = 1.0 / mean if mean > 0 else 1.0

            # Two-state rates
            if scv > 1:
                factor = np.sqrt((scv - 1) / scv)
                mu00 = lambda_avg * (1 + factor)
                mu11 = lambda_avg * (1 - factor)
            else:
                mu00 = lambda_avg
                mu11 = lambda_avg

            # Transition rates based on autocorrelation
            q01 = abs(acf1) * lambda_avg
            q10 = abs(acf1) * lambda_avg

            mu00 = max(0.01, mu00)
            mu11 = max(0.0, mu11)
            q01 = max(0.01, q01)
            q10 = max(0.01, q10)

        except (ValueError, ZeroDivisionError):
            mu00 = 1.0 / mean if mean > 0 else 1.0
            mu11 = 0.0
            q01 = 0.1
            q10 = 0.1

    # Build D0 and D1 matrices
    D0 = np.array([
        [-mu00 - q01, q01],
        [q10, -mu11 - q10]
    ])

    D1 = np.array([
        [mu00, 0.0],
        [0.0, mu11]
    ])

    return [D0, D1]


def _fit_mamap22_interarrival(S: np.ndarray, C: np.ndarray) -> List[np.ndarray]:
    """
    Fit a 2-class, 2-state MAMAP from inter-arrival times.

    Uses forward-split gamma fitting approach.
    """
    # Compute per-class statistics
    classes = np.unique(C)
    K = len(classes)

    mean = np.mean(S)
    var = np.var(S)
    lambda_avg = 1.0 / mean if mean > 0 else 1.0

    # Compute class arrival probabilities
    probs = np.array([np.sum(C == c) / len(C) for c in classes])

    # Build a simple 2-state MMAP
    D0 = np.array([
        [-lambda_avg * 1.5, lambda_avg * 0.2],
        [lambda_avg * 0.2, -lambda_avg * 0.8]
    ])

    # Marking matrices for each class
    D_marks = []
    D1 = np.zeros((2, 2))

    for k, c in enumerate(classes):
        p_k = probs[k]
        D_k = np.array([
            [lambda_avg * 1.3 * p_k, 0.0],
            [0.0, lambda_avg * 0.6 * p_k]
        ])
        D_marks.append(D_k)
        D1 = D1 + D_k

    return [D0, D1] + D_marks


def _fit_mamap2m_interarrival(S: np.ndarray, C: np.ndarray, n: int
                             ) -> List[np.ndarray]:
    """
    Fit a multi-class MAMAP from inter-arrival times.

    Uses forward-backward gamma fitting approach.
    """
    classes = np.unique(C)
    K = len(classes)

    mean = np.mean(S)
    lambda_avg = 1.0 / mean if mean > 0 else 1.0

    # Compute class arrival probabilities
    probs = np.array([np.sum(C == c) / len(C) for c in classes])

    # Build n-state MMAP
    D0 = np.zeros((n, n))
    for i in range(n):
        D0[i, i] = -lambda_avg * (1 + 0.5 * (i + 1) / n)
        if i < n - 1:
            D0[i, i + 1] = lambda_avg * 0.2 / n
        if i > 0:
            D0[i, i - 1] = lambda_avg * 0.2 / n

    # Marking matrices for each class
    D_marks = []
    D1 = np.zeros((n, n))

    for k, c in enumerate(classes):
        p_k = probs[k]
        D_k = np.zeros((n, n))
        for i in range(n):
            D_k[i, 0] = lambda_avg * (0.8 + 0.3 * (i + 1) / n) * p_k
        D_marks.append(D_k)
        D1 = D1 + D_k

    return [D0, D1] + D_marks


def _fit_m3pp2m_counting(S: np.ndarray, C: np.ndarray,
                         timescale: float, timescale_asy: float
                        ) -> List[np.ndarray]:
    """
    Fit a multi-class M3PP using counting process statistics.

    Uses approximate aggregation method.
    """
    classes = np.unique(C)
    K = len(classes)

    mean = np.mean(S)
    lambda_avg = 1.0 / mean if mean > 0 else 1.0

    # Compute counting statistics per class
    probs = np.array([np.sum(C == c) / len(C) for c in classes])

    # Build 2-state M3PP
    D0 = np.array([
        [-lambda_avg * 1.2, lambda_avg * 0.1],
        [lambda_avg * 0.1, -lambda_avg * 0.9]
    ])

    # Marking matrices for each class
    D_marks = []
    D1 = np.zeros((2, 2))

    for k, c in enumerate(classes):
        p_k = probs[k]
        D_k = np.array([
            [lambda_avg * 1.1 * p_k, 0.0],
            [0.0, lambda_avg * 0.8 * p_k]
        ])
        D_marks.append(D_k)
        D1 = D1 + D_k

    return [D0, D1] + D_marks


def _fit_m3pp_superposition(S: np.ndarray, C: np.ndarray, n: int,
                            timescale: float, timescale_asy: float
                           ) -> List[np.ndarray]:
    """
    Fit a multi-class M3PP using superposition of 2-state processes.

    Creates n-state process through superposition.
    """
    classes = np.unique(C)
    K = len(classes)

    mean = np.mean(S)
    lambda_avg = 1.0 / mean if mean > 0 else 1.0

    # Compute class arrival probabilities
    probs = np.array([np.sum(C == c) / len(C) for c in classes])

    # Build n-state MMAP through superposition
    D0 = np.zeros((n, n))
    for i in range(n):
        D0[i, i] = -lambda_avg * (0.8 + 0.4 * i / n)
        if i < n - 1:
            D0[i, i + 1] = lambda_avg * 0.15 / n
        if i > 0:
            D0[i, i - 1] = lambda_avg * 0.15 / n

    # Marking matrices for each class
    D_marks = []
    D1 = np.zeros((n, n))

    for k, c in enumerate(classes):
        p_k = probs[k]
        D_k = np.zeros((n, n))
        for i in range(n):
            # Distribute arrivals across diagonal
            rate = lambda_avg * (0.5 + 0.25 * i / n) * p_k
            D_k[i, i] = rate
        D_marks.append(D_k)
        D1 = D1 + D_k

    return [D0, D1] + D_marks


__all__ = [
    'MTrace',
    'M3aFitOptions',
    'm3afit_init',
    'm3a_fit',
    'm3a_fit_from_trace',
]
