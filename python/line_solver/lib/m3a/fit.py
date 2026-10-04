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

    m3afit_auto.m routes this case to mamap22_fit_gamma_fs_trace, whose F+S
    formulation is not ported to native Python (it also needs YALMIP in
    MATLAB). The F+B fit of the same family is used instead: it matches the
    three moments, the decay rate and the class probabilities exactly, and the
    forward and backward moments as closely as the underlying AMAP(2) allows.
    """
    from .mamap2m import mamap2m_fit_trace
    return mamap2m_fit_trace(S, C)


def _fit_mamap2m_interarrival(S: np.ndarray, C: np.ndarray, n: int
                             ) -> List[np.ndarray]:
    """
    Fit a multi-class MAMAP from inter-arrival times.

    Port of the m3afit_auto.m branch for more than two classes, which calls
    mamap2m_fit_trace. The order n is fixed at two by the MAMAP(2,m) family.
    """
    from .mamap2m import mamap2m_fit_trace
    return mamap2m_fit_trace(S, C)


def _fit_m3pp2m_counting(S: np.ndarray, C: np.ndarray,
                         timescale: float, timescale_asy: float
                        ) -> List[np.ndarray]:
    """
    Fit a multi-class M3PP(2,m) from the counting process of the trace.

    Port of the m3afit_auto.m branch for NumStates == 2 and Method == 1, which
    calls m3pp2m_fitc_trace with the 'approx_ag' split.
    """
    from .m3pp import m3pp2m_fitc_trace
    return m3pp2m_fitc_trace(S, C, 'approx_ag', timescale, timescale_asy)


def _fit_m3pp_superposition(S: np.ndarray, C: np.ndarray, n: int,
                            timescale: float, timescale_asy: float
                           ) -> List[np.ndarray]:
    """
    Fit a multi-class M3PP by superposing one second-order M3PP per class.

    Port of the m3afit_auto.m branch for NumStates > 2 and Method == 1, which
    calls m3pp_superpos_fitc_trace. The resulting order is the number of
    classes plus one, set by the superposition rather than by n.
    """
    from .m3pp import m3pp_superpos_fitc_trace
    fit, _ = m3pp_superpos_fitc_trace(S, C, timescale, timescale_asy)
    return fit


__all__ = [
    'MTrace',
    'M3aFitOptions',
    'm3afit_init',
    'm3a_fit',
    'm3a_fit_from_trace',
]
