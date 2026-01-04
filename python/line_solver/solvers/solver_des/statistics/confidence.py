"""
Confidence interval computation for DES solver.

This module implements the Overlapping Batch Means (OBM) method
for computing confidence intervals from correlated simulation output.
"""

from typing import List, Optional, Tuple, NamedTuple
import numpy as np
from scipy import stats


class ConfidenceInterval(NamedTuple):
    """
    Confidence interval result.

    Attributes:
        mean: Point estimate (sample mean)
        lower: Lower bound of confidence interval
        upper: Upper bound of confidence interval
        half_width: Half-width of the interval
        relative_precision: half_width / |mean| (if mean != 0)
    """
    mean: float
    lower: float
    upper: float
    half_width: float
    relative_precision: float


class OBMConfidenceInterval:
    """
    Overlapping Batch Means (OBM) confidence interval estimator.

    OBM uses overlapping batches to reduce variance in the variance
    estimator while maintaining unbiasedness. The method:

    1. Divides observations into overlapping batches
    2. Computes batch means
    3. Uses the variance of batch means to estimate the variance of the grand mean
    4. Applies a correction factor for overlap

    Reference:
        Meketon, M.S. and Schmeiser, B.W. (1984). "Overlapping batch means:
        Something for nothing?"

    Attributes:
        batch_size: Size of each batch
        overlap_fraction: Fraction of overlap between batches (0.5 = 50%)
    """

    def __init__(self, batch_size: int = 100, overlap_fraction: float = 0.5):
        """
        Initialize the OBM estimator.

        Args:
            batch_size: Size of each batch
            overlap_fraction: Fraction of overlap (0.0 to 0.9)
        """
        if not 0.0 <= overlap_fraction < 1.0:
            raise ValueError("Overlap fraction must be in [0, 1)")
        self.batch_size = batch_size
        self.overlap_fraction = overlap_fraction

        # Compute overlap step (non-overlap portion)
        self.step = max(1, int(batch_size * (1 - overlap_fraction)))

        # Correction factor for overlap
        # For 50% overlap, factor is 4/3
        self._compute_correction_factor()

    def _compute_correction_factor(self) -> None:
        """
        Compute the correction factor for overlapping batches.

        For OBM with overlap fraction p, the correction is:
        factor = 1 / (1 - p) for large samples
        More precisely: factor = b/l where b = batch_size, l = step
        """
        if self.step > 0:
            self.correction_factor = self.batch_size / self.step
        else:
            self.correction_factor = 1.0

    def compute(
        self,
        observations: List[float],
        confidence_level: float = 0.95
    ) -> ConfidenceInterval:
        """
        Compute confidence interval using OBM.

        Args:
            observations: List of observations
            confidence_level: Confidence level (0.90, 0.95, 0.99, etc.)

        Returns:
            ConfidenceInterval with mean, bounds, and precision
        """
        n = len(observations)
        data = np.array(observations)

        if n < self.batch_size:
            # Not enough data for batching, use simple CI
            return self._simple_ci(data, confidence_level)

        # Create overlapping batches
        batch_means = self._create_batch_means(data)
        n_batches = len(batch_means)

        if n_batches < 2:
            return self._simple_ci(data, confidence_level)

        # Grand mean
        grand_mean = data.mean()

        # Variance of batch means
        batch_var = np.var(batch_means, ddof=1)

        # Variance of grand mean with correction
        var_grand_mean = (batch_var / n_batches) * self.correction_factor

        # Degrees of freedom (approximate for OBM)
        df = n_batches - 1

        # t critical value
        alpha = 1 - confidence_level
        t_crit = stats.t.ppf(1 - alpha / 2, df)

        # Half-width
        std_error = np.sqrt(var_grand_mean)
        half_width = t_crit * std_error

        # Confidence interval
        lower = grand_mean - half_width
        upper = grand_mean + half_width

        # Relative precision
        if abs(grand_mean) > 1e-10:
            rel_precision = half_width / abs(grand_mean)
        else:
            rel_precision = float('inf')

        return ConfidenceInterval(
            mean=grand_mean,
            lower=lower,
            upper=upper,
            half_width=half_width,
            relative_precision=rel_precision
        )

    def _create_batch_means(self, data: np.ndarray) -> np.ndarray:
        """
        Create overlapping batch means.

        Args:
            data: Array of observations

        Returns:
            Array of batch means
        """
        n = len(data)
        batch_means = []

        start = 0
        while start + self.batch_size <= n:
            batch = data[start:start + self.batch_size]
            batch_means.append(batch.mean())
            start += self.step

        return np.array(batch_means)

    def _simple_ci(
        self,
        data: np.ndarray,
        confidence_level: float
    ) -> ConfidenceInterval:
        """
        Simple confidence interval (no batching).

        Args:
            data: Array of observations
            confidence_level: Confidence level

        Returns:
            ConfidenceInterval
        """
        n = len(data)
        if n < 2:
            mean = data.mean() if n > 0 else 0.0
            return ConfidenceInterval(
                mean=mean,
                lower=mean,
                upper=mean,
                half_width=0.0,
                relative_precision=0.0
            )

        mean = data.mean()
        std_err = data.std(ddof=1) / np.sqrt(n)

        alpha = 1 - confidence_level
        t_crit = stats.t.ppf(1 - alpha / 2, n - 1)
        half_width = t_crit * std_err

        if abs(mean) > 1e-10:
            rel_precision = half_width / abs(mean)
        else:
            rel_precision = float('inf')

        return ConfidenceInterval(
            mean=mean,
            lower=mean - half_width,
            upper=mean + half_width,
            half_width=half_width,
            relative_precision=rel_precision
        )


class NonOverlappingBatchMeans:
    """
    Non-overlapping batch means (NBM) confidence interval estimator.

    Simpler than OBM but less efficient. Divides observations into
    non-overlapping batches and treats batch means as independent.
    """

    def __init__(self, batch_size: int = 100):
        """
        Initialize NBM estimator.

        Args:
            batch_size: Size of each batch
        """
        self.batch_size = batch_size

    def compute(
        self,
        observations: List[float],
        confidence_level: float = 0.95
    ) -> ConfidenceInterval:
        """
        Compute confidence interval using NBM.

        Args:
            observations: List of observations
            confidence_level: Confidence level

        Returns:
            ConfidenceInterval
        """
        n = len(observations)
        data = np.array(observations)

        if n < self.batch_size * 2:
            # Not enough data
            mean = data.mean() if n > 0 else 0.0
            return ConfidenceInterval(
                mean=mean,
                lower=mean,
                upper=mean,
                half_width=0.0,
                relative_precision=0.0
            )

        # Create non-overlapping batches
        n_batches = n // self.batch_size
        trimmed = data[:n_batches * self.batch_size]
        batches = trimmed.reshape(n_batches, self.batch_size)
        batch_means = batches.mean(axis=1)

        # Compute CI from batch means
        grand_mean = batch_means.mean()
        batch_std = batch_means.std(ddof=1)
        std_err = batch_std / np.sqrt(n_batches)

        alpha = 1 - confidence_level
        t_crit = stats.t.ppf(1 - alpha / 2, n_batches - 1)
        half_width = t_crit * std_err

        if abs(grand_mean) > 1e-10:
            rel_precision = half_width / abs(grand_mean)
        else:
            rel_precision = float('inf')

        return ConfidenceInterval(
            mean=grand_mean,
            lower=grand_mean - half_width,
            upper=grand_mean + half_width,
            half_width=half_width,
            relative_precision=rel_precision
        )


class SpectralConfidenceInterval:
    """
    Spectral method for confidence interval estimation.

    Uses spectral analysis to estimate the variance of the mean
    for correlated observations.
    """

    def __init__(self, max_lag: Optional[int] = None):
        """
        Initialize spectral estimator.

        Args:
            max_lag: Maximum lag for autocorrelation (auto-selected if None)
        """
        self.max_lag = max_lag

    def compute(
        self,
        observations: List[float],
        confidence_level: float = 0.95
    ) -> ConfidenceInterval:
        """
        Compute confidence interval using spectral method.

        Args:
            observations: List of observations
            confidence_level: Confidence level

        Returns:
            ConfidenceInterval
        """
        n = len(observations)
        data = np.array(observations)

        if n < 10:
            mean = data.mean() if n > 0 else 0.0
            return ConfidenceInterval(
                mean=mean,
                lower=mean,
                upper=mean,
                half_width=0.0,
                relative_precision=0.0
            )

        mean = data.mean()
        centered = data - mean

        # Determine max lag
        if self.max_lag is None:
            max_lag = min(int(np.sqrt(n)), n // 4)
        else:
            max_lag = min(self.max_lag, n - 1)

        # Compute autocorrelations using truncated estimator
        variance = np.var(centered, ddof=0)
        if variance < 1e-15:
            return ConfidenceInterval(
                mean=mean,
                lower=mean,
                upper=mean,
                half_width=0.0,
                relative_precision=0.0
            )

        # Estimate asymptotic variance using Parzen window
        autocov_sum = variance  # Lag 0
        for lag in range(1, max_lag + 1):
            weight = 1 - lag / (max_lag + 1)  # Bartlett window
            autocov = np.mean(centered[:-lag] * centered[lag:])
            autocov_sum += 2 * weight * autocov

        # Variance of mean
        var_mean = autocov_sum / n

        # Standard error
        std_err = np.sqrt(max(0, var_mean))

        # Use normal distribution for large samples
        alpha = 1 - confidence_level
        z_crit = stats.norm.ppf(1 - alpha / 2)
        half_width = z_crit * std_err

        if abs(mean) > 1e-10:
            rel_precision = half_width / abs(mean)
        else:
            rel_precision = float('inf')

        return ConfidenceInterval(
            mean=mean,
            lower=mean - half_width,
            upper=mean + half_width,
            half_width=half_width,
            relative_precision=rel_precision
        )


def compute_confidence_interval(
    observations: List[float],
    method: str = 'OBM',
    confidence_level: float = 0.95,
    **kwargs
) -> ConfidenceInterval:
    """
    Compute confidence interval using specified method.

    Args:
        observations: List of observations
        method: Method to use ('OBM', 'NBM', 'SPECTRAL', 'SIMPLE')
        confidence_level: Confidence level
        **kwargs: Additional arguments for the specific method

    Returns:
        ConfidenceInterval
    """
    if method.upper() == 'OBM':
        batch_size = kwargs.get('batch_size', 100)
        overlap = kwargs.get('overlap_fraction', 0.5)
        estimator = OBMConfidenceInterval(batch_size, overlap)
    elif method.upper() == 'NBM':
        batch_size = kwargs.get('batch_size', 100)
        estimator = NonOverlappingBatchMeans(batch_size)
    elif method.upper() == 'SPECTRAL':
        max_lag = kwargs.get('max_lag', None)
        estimator = SpectralConfidenceInterval(max_lag)
    elif method.upper() == 'SIMPLE':
        data = np.array(observations)
        n = len(data)
        if n < 2:
            mean = data.mean() if n > 0 else 0.0
            return ConfidenceInterval(mean, mean, mean, 0.0, 0.0)
        mean = data.mean()
        std_err = data.std(ddof=1) / np.sqrt(n)
        alpha = 1 - confidence_level
        t_crit = stats.t.ppf(1 - alpha / 2, n - 1)
        half_width = t_crit * std_err
        rel_prec = half_width / abs(mean) if abs(mean) > 1e-10 else float('inf')
        return ConfidenceInterval(mean, mean - half_width, mean + half_width, half_width, rel_prec)
    else:
        raise ValueError(f"Unknown method: {method}")

    return estimator.compute(observations, confidence_level)


def relative_precision_stopping(
    observations: List[float],
    target_precision: float = 0.05,
    confidence_level: float = 0.95,
    min_observations: int = 100
) -> Tuple[bool, ConfidenceInterval]:
    """
    Check if relative precision target has been met.

    Can be used for sequential stopping in simulation runs.

    Args:
        observations: Current observations
        target_precision: Target relative half-width (e.g., 0.05 for 5%)
        confidence_level: Confidence level
        min_observations: Minimum observations before checking

    Returns:
        Tuple of (should_stop, current_ci)
    """
    if len(observations) < min_observations:
        return False, ConfidenceInterval(0, 0, 0, 0, float('inf'))

    ci = compute_confidence_interval(observations, 'OBM', confidence_level)

    should_stop = ci.relative_precision <= target_precision
    return should_stop, ci
