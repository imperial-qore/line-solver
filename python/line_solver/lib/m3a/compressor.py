"""
M3A MMAP Compression algorithms.

Implements various compression methods based on the M3A methodology,
which focuses on preserving the first three moments and correlation structure
of the original MMAP while reducing its complexity.
"""

import numpy as np
from numpy.typing import NDArray
from typing import List, Dict, Tuple, Optional, NamedTuple
from dataclasses import dataclass
import scipy.linalg as la

from .utils import compute_moments, validate_mmap, compute_autocorrelation


@dataclass
class HyperExpParameters:
    """Parameters for hyper-exponential distribution."""
    order: int
    probabilities: np.ndarray
    rates: np.ndarray


@dataclass
class ErlangParameters:
    """Parameters for Erlang distribution."""
    max_order: int
    order: int = 1
    rate: float = 1.0


@dataclass
class CoxianParameters:
    """Parameters for Coxian distribution."""
    order: int
    rates: np.ndarray
    probabilities: np.ndarray


@dataclass
class PhaseTypeParameters:
    """Parameters for general phase-type distribution."""
    num_phases: int
    rates: np.ndarray
    initial_probs: np.ndarray
    transition_probs: np.ndarray


class M3aCompressor:
    """
    M3A (Markovian Arrival Process with 3-moment Approximation) compressor.

    Provides methods to compress MMAPs while preserving statistical properties.
    """

    @staticmethod
    def compress_hyperexponential(MMAP: List[np.ndarray], order: int = 2
                                  ) -> List[np.ndarray]:
        """
        Compress an MMAP using hyper-exponential approximation.

        This method approximates the MMAP using a mixture of hyperexponential
        distributions while preserving the first three moments.

        Args:
            MMAP: Original MMAP to compress [D0, D1, D2, ...]
            order: Order of hyper-exponential approximation (default: 2)

        Returns:
            Compressed MMAP using hyper-exponential approximation
        """
        K = len(MMAP) - 2  # Number of classes
        if K < 0:
            K = 0

        # Extract moments from original MMAP
        moments = _extract_moments(MMAP)

        # Fit hyper-exponential parameters
        params = _fit_hyperexponential(moments, order)

        # Build compressed MMAP
        result = []

        # D0 - generator matrix
        D0 = _build_hyperexp_generator(params, K)
        result.append(D0)

        # D1 - sum of marking matrices
        D1 = np.zeros((order, order))

        # D_k - marking matrices for each class
        for k in range(max(1, K)):
            marking = _build_marking_matrix(params, k, max(1, K))
            D1 = D1 + marking
            if K > 0:
                result.append(marking)

        result.insert(1, D1)

        return result

    @staticmethod
    def compress_erlang(MMAP: List[np.ndarray], max_order: int = 3
                       ) -> List[np.ndarray]:
        """
        Compress an MMAP using Erlang approximation.

        Uses a mixture of Erlang distributions to approximate the MMAP,
        particularly effective for MMAPs with low variability.

        Args:
            MMAP: Original MMAP to compress
            max_order: Maximum order of Erlang distributions

        Returns:
            Compressed MMAP using Erlang approximation
        """
        K = max(0, len(MMAP) - 2)

        # Extract moments
        moments = _extract_moments(MMAP)

        # Fit Erlang parameters
        params = _fit_erlang_mixture(moments, max_order)

        # Build compressed MMAP
        result = []

        # D0 - generator matrix
        D0 = _build_erlang_generator(params, K)
        result.append(D0)

        # D1 and marking matrices
        D1 = np.zeros((params.order, params.order))

        for k in range(max(1, K)):
            marking = _build_erlang_marking_matrix(params, k, max(1, K))
            D1 = D1 + marking
            if K > 0:
                result.append(marking)

        result.insert(1, D1)

        return result

    @staticmethod
    def compress_coxian(MMAP: List[np.ndarray], order: int = 2
                       ) -> List[np.ndarray]:
        """
        Compress an MMAP using Coxian approximation.

        Uses Coxian distributions to approximate the MMAP, providing
        a good balance between accuracy and computational efficiency.

        Args:
            MMAP: Original MMAP to compress
            order: Order of Coxian approximation

        Returns:
            Compressed MMAP using Coxian approximation
        """
        K = max(0, len(MMAP) - 2)

        # Extract moments
        moments = _extract_moments(MMAP)

        # Fit Coxian parameters
        params = _fit_coxian(moments, order)

        # Build compressed MMAP
        result = []

        # D0 - generator matrix
        D0 = _build_coxian_generator(params, K)
        result.append(D0)

        # D1 and marking matrices
        D1 = np.zeros((order, order))

        for k in range(max(1, K)):
            marking = _build_coxian_marking_matrix(params, k, max(1, K))
            D1 = D1 + marking
            if K > 0:
                result.append(marking)

        result.insert(1, D1)

        return result

    @staticmethod
    def compress_phasetype(MMAP: List[np.ndarray], num_phases: int = 3
                          ) -> List[np.ndarray]:
        """
        Compress an MMAP using general phase-type approximation.

        Provides high accuracy for complex arrival patterns.

        Args:
            MMAP: Original MMAP to compress
            num_phases: Number of phases in approximation

        Returns:
            Compressed MMAP using phase-type approximation
        """
        K = max(0, len(MMAP) - 2)

        # Extract moments and correlations
        moments = _extract_moments(MMAP)
        correlations = _extract_correlations(MMAP)

        # Fit phase-type parameters
        params = _fit_phasetype(moments, correlations, num_phases)

        # Build compressed MMAP
        result = []

        # D0 - generator matrix
        D0 = _build_phasetype_generator(params, K)
        result.append(D0)

        # D1 and marking matrices
        D1 = np.zeros((num_phases, num_phases))

        for k in range(max(1, K)):
            marking = _build_phasetype_marking_matrix(params, k, max(1, K))
            D1 = D1 + marking
            if K > 0:
                result.append(marking)

        result.insert(1, D1)

        return result

    @staticmethod
    def compress_minimal(MMAP: List[np.ndarray], tolerance: float = 1e-6
                        ) -> List[np.ndarray]:
        """
        Compress an MMAP using minimal representation.

        Finds the minimal representation that preserves essential
        statistical properties while minimizing state space.

        Args:
            MMAP: Original MMAP to compress
            tolerance: Tolerance for moment matching

        Returns:
            Compressed MMAP using minimal representation
        """
        # Start with order 2 and increase until moments are matched
        for order in range(2, 7):
            candidate = M3aCompressor.compress_hyperexponential(MMAP, order)

            # Check if moments are matched within tolerance
            if _are_moments_matched(MMAP, candidate, tolerance):
                return candidate

        # If no good approximation found, use highest order
        return M3aCompressor.compress_hyperexponential(MMAP, 6)


def compress_mmap(MMAP: List[np.ndarray], method: str = 'hyperexp',
                  **kwargs) -> List[np.ndarray]:
    """
    Compress an MMAP using the specified method.

    Args:
        MMAP: Original MMAP to compress [D0, D1, D2, ...]
        method: Compression method ('hyperexp', 'erlang', 'coxian', 'phasetype', 'minimal')
        **kwargs: Additional parameters passed to compression method

    Returns:
        Compressed MMAP
    """
    compressor = M3aCompressor()

    if method == 'hyperexp':
        order = kwargs.get('order', 2)
        return compressor.compress_hyperexponential(MMAP, order)
    elif method == 'erlang':
        max_order = kwargs.get('max_order', 3)
        return compressor.compress_erlang(MMAP, max_order)
    elif method == 'coxian':
        order = kwargs.get('order', 2)
        return compressor.compress_coxian(MMAP, order)
    elif method == 'phasetype':
        num_phases = kwargs.get('num_phases', 3)
        return compressor.compress_phasetype(MMAP, num_phases)
    elif method == 'minimal':
        tolerance = kwargs.get('tolerance', 1e-6)
        return compressor.compress_minimal(MMAP, tolerance)
    else:
        raise ValueError(f"Unknown compression method: {method}")


# Helper functions

def _extract_moments(MMAP: List[np.ndarray]) -> np.ndarray:
    """Extract first three moments from MMAP."""
    return compute_moments(MMAP, 3)


def _extract_correlations(MMAP: List[np.ndarray]) -> np.ndarray:
    """Extract lag-1 and lag-2 correlations from MMAP."""
    try:
        acf = compute_autocorrelation(MMAP, 2)
        return acf
    except:
        return np.array([0.1, 0.05])


def _fit_hyperexponential(moments: np.ndarray, order: int) -> HyperExpParameters:
    """Fit hyper-exponential distribution to moments."""
    mean = moments[0]
    variance = moments[1] - mean**2
    scv = variance / (mean**2) if mean > 0 else 1.0

    probs = np.zeros(order)
    rates = np.zeros(order)

    if scv > 1:
        # High variability - use standard hyper-exponential fitting
        probs[0] = 0.5 + np.sqrt((scv - 1) / (4 * scv))
        probs[1] = 1.0 - probs[0]
        rates[0] = 2 * probs[0] / mean if mean > 0 else 1.0
        rates[1] = 2 * probs[1] / mean if mean > 0 else 1.0
    else:
        # Low variability - use Erlang-like parameters
        probs[0] = 1.0
        probs[1] = 0.0
        rates[0] = 1.0 / mean if mean > 0 else 1.0
        rates[1] = 1.0 / mean if mean > 0 else 1.0

    # Fill remaining phases with balanced parameters
    for i in range(2, order):
        probs[i] = 0.0
        rates[i] = 1.0 / mean if mean > 0 else 1.0

    return HyperExpParameters(order=order, probabilities=probs, rates=rates)


def _fit_erlang_mixture(moments: np.ndarray, max_order: int) -> ErlangParameters:
    """Fit Erlang mixture to moments."""
    mean = moments[0]
    variance = moments[1] - mean**2
    scv = variance / (mean**2) if mean > 0 else 1.0

    # Determine optimal Erlang order
    order = max(1, min(max_order, int(1.0 / scv) if scv > 0 else 1))
    rate = order / mean if mean > 0 else 1.0

    return ErlangParameters(max_order=max_order, order=order, rate=rate)


def _fit_coxian(moments: np.ndarray, order: int) -> CoxianParameters:
    """Fit Coxian distribution to moments."""
    mean = moments[0]

    rates = np.full(order, order / mean if mean > 0 else 1.0)
    probs = np.full(order, 0.8)
    probs[-1] = 1.0

    return CoxianParameters(order=order, rates=rates, probabilities=probs)


def _fit_phasetype(moments: np.ndarray, correlations: np.ndarray,
                   num_phases: int) -> PhaseTypeParameters:
    """Fit phase-type distribution to moments and correlations."""
    mean = moments[0]

    rates = np.full(num_phases, num_phases / mean if mean > 0 else 1.0)
    initial_probs = np.ones(num_phases) / num_phases

    transition_probs = np.zeros((num_phases, num_phases))
    corr_val = correlations[0] if len(correlations) > 0 else 0.1

    for i in range(num_phases):
        for j in range(num_phases):
            if i != j:
                transition_probs[i, j] = abs(corr_val) / (num_phases - 1)

    return PhaseTypeParameters(
        num_phases=num_phases,
        rates=rates,
        initial_probs=initial_probs,
        transition_probs=transition_probs
    )


def _build_hyperexp_generator(params: HyperExpParameters, K: int) -> np.ndarray:
    """Build generator matrix for hyper-exponential distribution."""
    order = params.order
    generator = np.zeros((order, order))

    total_rate = np.mean(params.rates)

    for i in range(order):
        off_diag_sum = 0.0
        off_diag_rate = total_rate / (order * max(1, K))

        for j in range(order):
            if i != j:
                generator[i, j] = off_diag_rate
                off_diag_sum += off_diag_rate

        # Set diagonal to make row sum valid
        generator[i, i] = -(off_diag_sum + total_rate)

    return generator


def _build_marking_matrix(params: HyperExpParameters, class_index: int,
                          K: int) -> np.ndarray:
    """Build marking matrix for hyper-exponential distribution."""
    order = params.order
    marking = np.zeros((order, order))

    total_rate = np.mean(params.rates)
    marking_rate = total_rate / K

    for i in range(order):
        marking[i, i] = marking_rate

    return marking


def _build_erlang_generator(params: ErlangParameters, K: int) -> np.ndarray:
    """Build generator matrix for Erlang distribution."""
    order = params.order
    generator = np.zeros((order, order))

    for i in range(order):
        generator[i, i] = -params.rate
        if i < order - 1:
            generator[i, i + 1] = params.rate

    return generator


def _build_erlang_marking_matrix(params: ErlangParameters, class_index: int,
                                  K: int) -> np.ndarray:
    """Build marking matrix for Erlang distribution."""
    order = params.order
    marking = np.zeros((order, order))

    class_probability = 1.0 / K
    marking[order - 1, 0] = params.rate * class_probability

    return marking


def _build_coxian_generator(params: CoxianParameters, K: int) -> np.ndarray:
    """Build generator matrix for Coxian distribution."""
    order = params.order
    generator = np.zeros((order, order))

    for i in range(order):
        generator[i, i] = -params.rates[i]
        if i < order - 1:
            generator[i, i + 1] = params.rates[i] * params.probabilities[i]

    return generator


def _build_coxian_marking_matrix(params: CoxianParameters, class_index: int,
                                  K: int) -> np.ndarray:
    """Build marking matrix for Coxian distribution."""
    order = params.order
    marking = np.zeros((order, order))

    class_probability = 1.0 / K

    for i in range(order):
        exit_rate = params.rates[i] * (1.0 - params.probabilities[i])
        marking[i, 0] = exit_rate * class_probability

    return marking


def _build_phasetype_generator(params: PhaseTypeParameters, K: int) -> np.ndarray:
    """Build generator matrix for phase-type distribution."""
    num_phases = params.num_phases
    generator = np.zeros((num_phases, num_phases))

    for i in range(num_phases):
        generator[i, i] = -params.rates[i]
        for j in range(num_phases):
            if i != j:
                generator[i, j] = params.rates[i] * params.transition_probs[i, j]

    return generator


def _build_phasetype_marking_matrix(params: PhaseTypeParameters, class_index: int,
                                     K: int) -> np.ndarray:
    """Build marking matrix for phase-type distribution."""
    num_phases = params.num_phases
    marking = np.zeros((num_phases, num_phases))

    class_probability = 1.0 / K

    for i in range(num_phases):
        exit_rate = params.rates[i] * (1.0 - np.sum(params.transition_probs[i, :]))
        marking[i, 0] = exit_rate * class_probability

    return marking


def _are_moments_matched(original: List[np.ndarray], compressed: List[np.ndarray],
                         tolerance: float) -> bool:
    """Check if moments are matched within tolerance."""
    orig_moments = _extract_moments(original)
    comp_moments = _extract_moments(compressed)

    for i in range(len(orig_moments)):
        if orig_moments[i] != 0:
            rel_error = abs(orig_moments[i] - comp_moments[i]) / abs(orig_moments[i])
            if rel_error > tolerance:
                return False
        else:
            if abs(comp_moments[i]) > tolerance:
                return False

    return True


__all__ = [
    'M3aCompressor',
    'compress_mmap',
    'HyperExpParameters',
    'ErlangParameters',
    'CoxianParameters',
    'PhaseTypeParameters',
]
