"""
Warmup detection algorithms for DES solver.

This module implements the MSER-5 (Marginal Standard Error Rules)
algorithm for detecting transient periods in simulation output.
"""

from typing import List, Optional, Tuple
import numpy as np


class MSER5TransientDetector:
    """
    MSER-5 transient detection algorithm.

    MSER-5 finds the optimal truncation point that minimizes the
    marginal standard error of the mean. The algorithm:

    1. Batches observations into groups of 5
    2. For each potential truncation point d, computes:
       MSE(d) = Var(Z[d:n]) / (n-d)^2
    3. Finds d* that minimizes MSE(d)

    The truncation point d* represents the end of the warmup period.

    Reference:
        White, K.P., Jr. (1997). "An effective truncation heuristic for
        bias reduction in simulation output"

    Attributes:
        batch_size: Size of batches (default 5 for MSER-5)
        min_batches: Minimum batches to keep after truncation
    """

    def __init__(self, batch_size: int = 5, min_batches: int = 10):
        """
        Initialize the MSER-5 detector.

        Args:
            batch_size: Batch size (use 5 for MSER-5)
            min_batches: Minimum batches to retain after truncation
        """
        self.batch_size = batch_size
        self.min_batches = min_batches

    def detect_warmup(self, observations: List[float]) -> int:
        """
        Detect the warmup period using MSER-5.

        Args:
            observations: List of simulation observations (e.g., interval averages)

        Returns:
            Index of the first observation after warmup (truncation point).
            Returns 0 if insufficient data or no warmup detected.
        """
        if len(observations) < self.batch_size * self.min_batches:
            return 0

        # Create batched observations
        batched = self._create_batches(observations)
        n_batches = len(batched)

        if n_batches < self.min_batches + 1:
            return 0

        # Find optimal truncation point
        best_d = 0
        best_mse = float('inf')

        # Maximum truncation: leave at least min_batches
        max_d = n_batches - self.min_batches

        for d in range(max_d + 1):
            # Remaining observations after truncation
            remaining = batched[d:]
            n_remaining = len(remaining)

            if n_remaining < self.min_batches:
                break

            # Compute variance and MSE
            var = np.var(remaining, ddof=1)
            mse = var / (n_remaining * n_remaining)

            if mse < best_mse:
                best_mse = mse
                best_d = d

        # Convert batch index to observation index
        return best_d * self.batch_size

    def detect_warmup_observations(self, observations: np.ndarray) -> int:
        """
        Detect warmup for array of observations.

        Args:
            observations: NumPy array of observations

        Returns:
            Truncation point (index)
        """
        return self.detect_warmup(list(observations))

    def _create_batches(self, observations: List[float]) -> np.ndarray:
        """
        Create batched observations by averaging consecutive values.

        Args:
            observations: Raw observations

        Returns:
            Array of batch means
        """
        n = len(observations)
        n_batches = n // self.batch_size
        trimmed = observations[:n_batches * self.batch_size]
        reshaped = np.array(trimmed).reshape(n_batches, self.batch_size)
        return reshaped.mean(axis=1)


class WindowedMSER5Detector:
    """
    Windowed MSER-5 for online warmup detection.

    Maintains a window of recent observations and periodically
    re-evaluates the optimal truncation point.

    Useful for long-running simulations where early detection
    of warmup end allows for earlier data collection.
    """

    def __init__(
        self,
        window_size: int = 1000,
        batch_size: int = 5,
        check_interval: int = 100
    ):
        """
        Initialize windowed detector.

        Args:
            window_size: Maximum observations to keep
            batch_size: Batch size for MSER-5
            check_interval: How often to re-check (in observations)
        """
        self.window_size = window_size
        self.batch_size = batch_size
        self.check_interval = check_interval
        self.detector = MSER5TransientDetector(batch_size)

        self.observations: List[float] = []
        self.observation_count = 0
        self.warmup_detected = False
        self.warmup_index = -1

    def add_observation(self, value: float) -> bool:
        """
        Add an observation and check for warmup end.

        Args:
            value: Observation value

        Returns:
            True if warmup has been detected to end
        """
        self.observations.append(value)
        self.observation_count += 1

        # Trim window if necessary
        if len(self.observations) > self.window_size:
            self.observations = self.observations[-self.window_size:]

        # Periodic check
        if (
            not self.warmup_detected and
            self.observation_count % self.check_interval == 0
        ):
            truncation = self.detector.detect_warmup(self.observations)
            # Consider warmup ended if truncation point is stable
            if truncation > 0 and truncation < len(self.observations) // 2:
                self.warmup_detected = True
                self.warmup_index = self.observation_count - len(self.observations) + truncation

        return self.warmup_detected

    def reset(self) -> None:
        """Reset the detector."""
        self.observations = []
        self.observation_count = 0
        self.warmup_detected = False
        self.warmup_index = -1


class MultivariateMSER5:
    """
    MSER-5 for multiple metrics.

    Uses the maximum truncation point across all metrics
    to ensure all metrics have reached steady state.
    """

    def __init__(self, num_metrics: int, batch_size: int = 5):
        """
        Initialize multivariate detector.

        Args:
            num_metrics: Number of metrics to track
            batch_size: Batch size for MSER-5
        """
        self.num_metrics = num_metrics
        self.detectors = [MSER5TransientDetector(batch_size) for _ in range(num_metrics)]
        self.observations: List[List[float]] = [[] for _ in range(num_metrics)]

    def add_observation(self, metric_id: int, value: float) -> None:
        """
        Add an observation for a specific metric.

        Args:
            metric_id: Which metric (0-indexed)
            value: Observation value
        """
        self.observations[metric_id].append(value)

    def add_observations(self, values: List[float]) -> None:
        """
        Add observations for all metrics at once.

        Args:
            values: List of values, one per metric
        """
        for i, v in enumerate(values):
            self.observations[i].append(v)

    def detect_warmup(self) -> int:
        """
        Detect warmup using the maximum truncation across metrics.

        Returns:
            Maximum truncation point (conservative)
        """
        max_truncation = 0
        for i, obs in enumerate(self.observations):
            if obs:
                trunc = self.detectors[i].detect_warmup(obs)
                max_truncation = max(max_truncation, trunc)
        return max_truncation

    def reset(self) -> None:
        """Reset all observations."""
        self.observations = [[] for _ in range(self.num_metrics)]


def schruben_rule(observations: List[float], alpha: float = 0.05) -> int:
    """
    Schruben's initialization bias test.

    Tests for initialization bias using a regression-based approach.
    Returns the truncation point where bias becomes insignificant.

    Args:
        observations: List of observations
        alpha: Significance level for the test

    Returns:
        Truncation point
    """
    from scipy import stats

    n = len(observations)
    if n < 20:
        return 0

    # Test successively larger truncation points
    for d in range(1, n // 2):
        remaining = observations[d:]
        m = len(remaining)

        # Perform linear regression
        x = np.arange(m)
        slope, intercept, r_value, p_value, std_err = stats.linregress(x, remaining)

        # If slope is not significantly different from zero, bias has been removed
        if p_value > alpha:
            return d

    # Default to MSER-5 result if test fails
    detector = MSER5TransientDetector()
    return detector.detect_warmup(observations)


def welch_method(
    observations: List[float],
    window_size: int = 50
) -> Tuple[int, np.ndarray]:
    """
    Welch's graphical method for warmup detection.

    Computes moving averages across replications (or within a single
    replication using a sliding window) to identify the transient period.

    Args:
        observations: List of observations
        window_size: Size of moving average window

    Returns:
        Tuple of (estimated truncation point, moving averages)
    """
    n = len(observations)
    if n < window_size * 2:
        return 0, np.array(observations)

    # Compute centered moving averages
    moving_avg = np.convolve(
        observations,
        np.ones(window_size) / window_size,
        mode='valid'
    )

    # Find where the derivative stabilizes
    derivatives = np.abs(np.diff(moving_avg))
    threshold = np.percentile(derivatives[-len(derivatives)//4:], 95) * 2

    for i, deriv in enumerate(derivatives):
        # Look for sustained low derivative
        if i >= len(derivatives) - 10:
            break
        if all(derivatives[i:i+10] < threshold):
            return i + window_size // 2, moving_avg

    return 0, moving_avg
