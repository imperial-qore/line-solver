"""
Transient analysis support for DES solver.

This module provides classes for collecting time-series data during simulation
for transient (non-steady-state) analysis.
"""

from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass, field
import numpy as np


@dataclass
class TransientSample:
    """
    A single transient sample at a point in time.

    Attributes:
        time: Simulation time of the sample
        queue_lengths: Queue lengths [stations x classes]
        utilizations: Utilizations [stations x classes]
        throughputs: Instantaneous throughputs [stations x classes]
    """
    time: float
    queue_lengths: np.ndarray
    utilizations: np.ndarray
    throughputs: np.ndarray


class TransientCollector:
    """
    Collector for transient (time-varying) metrics during simulation.

    Samples metrics at regular intervals to capture system dynamics.

    Attributes:
        sample_interval: Time between samples
        num_stations: Number of stations in the network
        num_classes: Number of job classes
        samples: List of collected samples
    """

    def __init__(
        self,
        num_stations: int,
        num_classes: int,
        sample_interval: float = 1.0,
        max_samples: Optional[int] = None
    ):
        """
        Initialize transient collector.

        Args:
            num_stations: Number of stations
            num_classes: Number of job classes
            sample_interval: Time between samples (default: 1.0)
            max_samples: Maximum samples to keep (None = unlimited)
        """
        self.num_stations = num_stations
        self.num_classes = num_classes
        self.sample_interval = sample_interval
        self.max_samples = max_samples

        self.samples: List[TransientSample] = []
        self.last_sample_time = 0.0

        # Running accumulators for throughput calculation
        self._departures = np.zeros((num_stations, num_classes))
        self._last_departures = np.zeros((num_stations, num_classes))
        self._last_throughput_time = 0.0

    def should_sample(self, current_time: float) -> bool:
        """
        Check if a sample should be taken.

        Args:
            current_time: Current simulation time

        Returns:
            True if sample should be taken
        """
        return current_time >= self.last_sample_time + self.sample_interval

    def take_sample(
        self,
        current_time: float,
        queue_lengths: np.ndarray,
        busy_servers: np.ndarray,
        num_servers: np.ndarray
    ) -> None:
        """
        Take a transient sample.

        Args:
            current_time: Current simulation time
            queue_lengths: Current queue lengths [stations x classes]
            busy_servers: Current busy server counts [stations x classes]
            num_servers: Number of servers at each station [stations]
        """
        # Compute utilizations
        utilizations = np.zeros((self.num_stations, self.num_classes))
        for i in range(self.num_stations):
            if num_servers[i] > 0:
                for k in range(self.num_classes):
                    utilizations[i, k] = busy_servers[i, k] / num_servers[i]

        # Compute instantaneous throughputs
        dt = current_time - self._last_throughput_time
        if dt > 0:
            throughputs = (self._departures - self._last_departures) / dt
        else:
            throughputs = np.zeros((self.num_stations, self.num_classes))

        sample = TransientSample(
            time=current_time,
            queue_lengths=queue_lengths.copy(),
            utilizations=utilizations,
            throughputs=throughputs
        )

        self.samples.append(sample)
        self.last_sample_time = current_time
        self._last_departures = self._departures.copy()
        self._last_throughput_time = current_time

        # Trim if max_samples exceeded
        if self.max_samples is not None and len(self.samples) > self.max_samples:
            self.samples = self.samples[-self.max_samples:]

    def record_departure(
        self,
        station_id: int,
        class_id: int
    ) -> None:
        """
        Record a departure for throughput calculation.

        Args:
            station_id: Station index
            class_id: Class index
        """
        self._departures[station_id, class_id] += 1

    def get_time_series(self) -> Dict[str, np.ndarray]:
        """
        Get transient time series as numpy arrays.

        Returns:
            Dictionary with:
                t: Time points [num_samples]
                QNt: Queue lengths [num_samples x stations x classes]
                UNt: Utilizations [num_samples x stations x classes]
                TNt: Throughputs [num_samples x stations x classes]
        """
        if not self.samples:
            return {
                't': np.array([]),
                'QNt': np.array([]).reshape(0, self.num_stations, self.num_classes),
                'UNt': np.array([]).reshape(0, self.num_stations, self.num_classes),
                'TNt': np.array([]).reshape(0, self.num_stations, self.num_classes),
            }

        n_samples = len(self.samples)

        t = np.array([s.time for s in self.samples])
        QNt = np.zeros((n_samples, self.num_stations, self.num_classes))
        UNt = np.zeros((n_samples, self.num_stations, self.num_classes))
        TNt = np.zeros((n_samples, self.num_stations, self.num_classes))

        for i, sample in enumerate(self.samples):
            QNt[i] = sample.queue_lengths
            UNt[i] = sample.utilizations
            TNt[i] = sample.throughputs

        return {
            't': t,
            'QNt': QNt,
            'UNt': UNt,
            'TNt': TNt,
        }

    def get_windowed_averages(
        self,
        window_size: int = 10
    ) -> Dict[str, np.ndarray]:
        """
        Get windowed (moving average) time series.

        Useful for smoothing noisy transient data.

        Args:
            window_size: Number of samples to average

        Returns:
            Dictionary with smoothed time series
        """
        raw = self.get_time_series()

        if len(raw['t']) < window_size:
            return raw

        def moving_avg(arr: np.ndarray, w: int) -> np.ndarray:
            """Apply moving average along first axis."""
            if arr.ndim == 1:
                return np.convolve(arr, np.ones(w)/w, mode='valid')
            else:
                # For multi-dimensional, average along first axis
                result = np.zeros((arr.shape[0] - w + 1,) + arr.shape[1:])
                for i in range(result.shape[0]):
                    result[i] = arr[i:i+w].mean(axis=0)
                return result

        return {
            't': moving_avg(raw['t'], window_size),
            'QNt': moving_avg(raw['QNt'], window_size),
            'UNt': moving_avg(raw['UNt'], window_size),
            'TNt': moving_avg(raw['TNt'], window_size),
        }

    def reset(self) -> None:
        """Reset the collector."""
        self.samples = []
        self.last_sample_time = 0.0
        self._departures = np.zeros((self.num_stations, self.num_classes))
        self._last_departures = np.zeros((self.num_stations, self.num_classes))
        self._last_throughput_time = 0.0


class SteadyStateDetector:
    """
    Detector for determining when system has reached steady state.

    Uses various statistical tests to detect when transient effects
    have dissipated and steady state has been reached.
    """

    def __init__(
        self,
        num_metrics: int,
        window_size: int = 50,
        threshold: float = 0.05
    ):
        """
        Initialize steady state detector.

        Args:
            num_metrics: Number of metrics to track
            window_size: Window size for comparison
            threshold: Relative change threshold for steady state
        """
        self.num_metrics = num_metrics
        self.window_size = window_size
        self.threshold = threshold

        self.observations: List[np.ndarray] = []

    def add_observation(self, values: np.ndarray) -> bool:
        """
        Add an observation and check for steady state.

        Args:
            values: Array of metric values

        Returns:
            True if steady state detected
        """
        self.observations.append(values.copy())

        if len(self.observations) < 2 * self.window_size:
            return False

        # Compare recent window to previous window
        recent = np.array(self.observations[-self.window_size:])
        previous = np.array(self.observations[-2*self.window_size:-self.window_size])

        recent_mean = recent.mean(axis=0)
        previous_mean = previous.mean(axis=0)

        # Relative change
        with np.errstate(divide='ignore', invalid='ignore'):
            rel_change = np.abs(recent_mean - previous_mean) / np.abs(previous_mean)
            rel_change = np.nan_to_num(rel_change, nan=0.0, posinf=0.0)

        return np.all(rel_change < self.threshold)

    def get_steady_state_time(self) -> Optional[int]:
        """
        Get the observation index where steady state was first detected.

        Returns:
            Index or None if not detected
        """
        for i in range(2 * self.window_size, len(self.observations)):
            # Test from this point
            recent = np.array(self.observations[i-self.window_size:i])
            previous = np.array(self.observations[i-2*self.window_size:i-self.window_size])

            recent_mean = recent.mean(axis=0)
            previous_mean = previous.mean(axis=0)

            with np.errstate(divide='ignore', invalid='ignore'):
                rel_change = np.abs(recent_mean - previous_mean) / np.abs(previous_mean)
                rel_change = np.nan_to_num(rel_change, nan=0.0, posinf=0.0)

            if np.all(rel_change < self.threshold):
                return i - self.window_size

        return None

    def reset(self) -> None:
        """Reset the detector."""
        self.observations = []


def compute_transient_statistics(
    time_series: Dict[str, np.ndarray],
    time_points: Optional[List[float]] = None
) -> Dict[str, np.ndarray]:
    """
    Compute statistics at specific time points from transient data.

    Args:
        time_series: Time series from TransientCollector.get_time_series()
        time_points: Time points to compute statistics at (interpolated)

    Returns:
        Dictionary with interpolated values at specified time points
    """
    t = time_series['t']

    if len(t) == 0:
        return {}

    if time_points is None:
        time_points = [t[-1]]

    result = {'t': np.array(time_points)}

    for key in ['QNt', 'UNt', 'TNt']:
        if key not in time_series:
            continue

        data = time_series[key]
        n_points = len(time_points)
        n_stations = data.shape[1] if data.ndim > 1 else 1
        n_classes = data.shape[2] if data.ndim > 2 else 1

        interpolated = np.zeros((n_points, n_stations, n_classes))

        for tp_idx, tp in enumerate(time_points):
            # Find closest sample
            idx = np.searchsorted(t, tp)
            if idx >= len(t):
                idx = len(t) - 1
            elif idx > 0 and tp - t[idx-1] < t[idx] - tp:
                idx = idx - 1

            interpolated[tp_idx] = data[idx]

        result[key.replace('t', '')] = interpolated

    return result
