"""
Statistics collection and metrics for DES solver.

This module provides classes for collecting time-weighted statistics
and response time observations during simulation.
"""

from typing import List, Optional, Dict, Any
from dataclasses import dataclass, field
import numpy as np


@dataclass
class TimeWeightedAccumulator:
    """
    Accumulator for time-weighted statistics.

    Used for metrics like queue length and utilization where
    the value changes over time and we need the time-weighted average.

    Attributes:
        current_value: Current value being tracked
        weighted_sum: Accumulated time-weighted sum
        last_update_time: Time of last update
        total_time: Total time elapsed
    """
    current_value: float = 0.0
    weighted_sum: float = 0.0
    last_update_time: float = 0.0
    total_time: float = 0.0

    def update(self, new_value: float, current_time: float) -> None:
        """
        Update the accumulator with a new value at the given time.

        Args:
            new_value: New value to track
            current_time: Current simulation time
        """
        if current_time > self.last_update_time:
            # Add weighted contribution from previous value
            dt = current_time - self.last_update_time
            self.weighted_sum += self.current_value * dt
            self.total_time += dt

        self.current_value = new_value
        self.last_update_time = current_time

    def increment(self, delta: float, current_time: float) -> None:
        """
        Increment the current value.

        Args:
            delta: Amount to add to current value
            current_time: Current simulation time
        """
        self.update(self.current_value + delta, current_time)

    def get_mean(self, current_time: Optional[float] = None) -> float:
        """
        Get the time-weighted mean.

        Args:
            current_time: End time for calculation (uses last update if None)

        Returns:
            Time-weighted average value
        """
        if current_time is not None and current_time > self.last_update_time:
            dt = current_time - self.last_update_time
            total_weighted = self.weighted_sum + self.current_value * dt
            total_time = self.total_time + dt
        else:
            total_weighted = self.weighted_sum
            total_time = self.total_time

        return total_weighted / total_time if total_time > 0 else 0.0

    def reset(self, time: float = 0.0) -> None:
        """Reset the accumulator."""
        self.current_value = 0.0
        self.weighted_sum = 0.0
        self.last_update_time = time
        self.total_time = 0.0


@dataclass
class ResponseTimeTally:
    """
    Tally for response time observations.

    Collects individual response times for statistical analysis,
    confidence interval computation, and percentile estimation.

    Attributes:
        observations: List of response time observations
        count: Number of observations
        sum: Sum of observations
        sum_sq: Sum of squared observations
    """
    observations: List[float] = field(default_factory=list)
    count: int = 0
    sum: float = 0.0
    sum_sq: float = 0.0

    # For interval-based sampling (MSER-5)
    interval_observations: List[float] = field(default_factory=list)
    interval_sum: float = 0.0
    interval_count: int = 0
    interval_start_time: float = 0.0

    def add(self, value: float) -> None:
        """
        Add a response time observation.

        Args:
            value: Response time to add
        """
        self.observations.append(value)
        self.count += 1
        self.sum += value
        self.sum_sq += value * value

        # Also add to current interval
        self.interval_sum += value
        self.interval_count += 1

    def close_interval(self, current_time: float) -> Optional[float]:
        """
        Close the current interval and return the interval average.

        Args:
            current_time: Current simulation time

        Returns:
            Average value for the interval, or None if empty
        """
        if self.interval_count > 0:
            avg = self.interval_sum / self.interval_count
            self.interval_observations.append(avg)
            self.interval_sum = 0.0
            self.interval_count = 0
            self.interval_start_time = current_time
            return avg
        return None

    def get_mean(self) -> float:
        """Get the sample mean."""
        return self.sum / self.count if self.count > 0 else 0.0

    def get_variance(self) -> float:
        """Get the sample variance."""
        if self.count < 2:
            return 0.0
        mean = self.get_mean()
        return (self.sum_sq - self.count * mean * mean) / (self.count - 1)

    def get_std(self) -> float:
        """Get the sample standard deviation."""
        return np.sqrt(self.get_variance())

    def get_percentile(self, p: float) -> float:
        """
        Get the p-th percentile.

        Args:
            p: Percentile (0-100)

        Returns:
            Value at the p-th percentile
        """
        if not self.observations:
            return 0.0
        return float(np.percentile(self.observations, p))

    def reset(self, time: float = 0.0) -> None:
        """Reset the tally."""
        self.observations = []
        self.count = 0
        self.sum = 0.0
        self.sum_sq = 0.0
        self.interval_observations = []
        self.interval_sum = 0.0
        self.interval_count = 0
        self.interval_start_time = time


@dataclass
class ClassMetrics:
    """
    Metrics for a single job class at a station.

    Attributes:
        queue_length: Time-weighted queue length accumulator
        busy_servers: Time-weighted busy server count
        arrivals: Total arrivals
        departures: Total departures
        response_times: Response time tally
    """
    queue_length: TimeWeightedAccumulator = field(default_factory=TimeWeightedAccumulator)
    busy_servers: TimeWeightedAccumulator = field(default_factory=TimeWeightedAccumulator)
    arrivals: int = 0
    departures: int = 0
    response_times: ResponseTimeTally = field(default_factory=ResponseTimeTally)

    # For tracking utilization per server
    server_busy_time: float = 0.0

    def reset(self, time: float = 0.0) -> None:
        """Reset all metrics."""
        self.queue_length.reset(time)
        self.busy_servers.reset(time)
        self.arrivals = 0
        self.departures = 0
        self.response_times.reset(time)
        self.server_busy_time = 0.0


@dataclass
class StationMetrics:
    """
    Metrics for a station with multiple job classes.

    Attributes:
        num_classes: Number of job classes
        num_servers: Number of servers
        class_metrics: Per-class metrics
        total_queue_length: Aggregate queue length
        total_busy_servers: Aggregate busy server count
    """
    num_classes: int
    num_servers: int
    class_metrics: List[ClassMetrics] = field(default_factory=list)
    total_queue_length: TimeWeightedAccumulator = field(default_factory=TimeWeightedAccumulator)
    total_busy_servers: TimeWeightedAccumulator = field(default_factory=TimeWeightedAccumulator)

    def __post_init__(self):
        if not self.class_metrics:
            self.class_metrics = [ClassMetrics() for _ in range(self.num_classes)]

    def record_arrival(self, class_id: int, current_time: float) -> None:
        """Record an arrival for a class."""
        self.class_metrics[class_id].arrivals += 1

    def record_departure(
        self,
        class_id: int,
        response_time: float,
        current_time: float
    ) -> None:
        """Record a departure for a class."""
        cm = self.class_metrics[class_id]
        cm.departures += 1
        cm.response_times.add(response_time)

    def update_queue_length(
        self,
        class_id: int,
        queue_length: int,
        current_time: float
    ) -> None:
        """Update queue length for a class."""
        self.class_metrics[class_id].queue_length.update(queue_length, current_time)

        # Also update total
        total = sum(cm.queue_length.current_value for cm in self.class_metrics)
        self.total_queue_length.update(total, current_time)

    def update_busy_servers(
        self,
        class_id: int,
        busy_count: int,
        current_time: float
    ) -> None:
        """Update busy server count for a class."""
        self.class_metrics[class_id].busy_servers.update(busy_count, current_time)

        # Also update total
        total = sum(cm.busy_servers.current_value for cm in self.class_metrics)
        self.total_busy_servers.update(total, current_time)

    def get_queue_length(self, class_id: int, current_time: float) -> float:
        """Get time-weighted average queue length for a class."""
        return self.class_metrics[class_id].queue_length.get_mean(current_time)

    def get_utilization(self, class_id: int, current_time: float) -> float:
        """Get utilization for a class (fraction of server capacity)."""
        if self.num_servers <= 0:
            return 0.0
        busy = self.class_metrics[class_id].busy_servers.get_mean(current_time)
        return busy / self.num_servers

    def get_throughput(self, class_id: int, current_time: float, start_time: float = 0.0) -> float:
        """Get throughput (departures per unit time) for a class."""
        elapsed = current_time - start_time
        if elapsed <= 0:
            return 0.0
        return self.class_metrics[class_id].departures / elapsed

    def get_response_time(self, class_id: int) -> float:
        """Get average response time for a class."""
        return self.class_metrics[class_id].response_times.get_mean()

    def get_residence_time(self, class_id: int, current_time: float, start_time: float = 0.0) -> float:
        """
        Get average residence time (response time) using Little's Law.

        W = Q / lambda where Q is queue length and lambda is arrival rate.
        """
        Q = self.get_queue_length(class_id, current_time)
        throughput = self.get_throughput(class_id, current_time, start_time)
        if throughput <= 0:
            return 0.0
        return Q / throughput

    def close_intervals(self, current_time: float) -> None:
        """Close current interval for all classes."""
        for cm in self.class_metrics:
            cm.response_times.close_interval(current_time)

    def reset(self, time: float = 0.0) -> None:
        """Reset all metrics."""
        for cm in self.class_metrics:
            cm.reset(time)
        self.total_queue_length.reset(time)
        self.total_busy_servers.reset(time)


class SimulationMetrics:
    """
    Container for all simulation metrics across stations and classes.

    Attributes:
        station_metrics: Metrics for each station
        start_time: Simulation start time
        warmup_time: Time when warmup ended (for truncation)
        interval_length: Length of each observation interval
    """

    def __init__(
        self,
        num_stations: int,
        num_classes: int,
        servers_per_station: List[int],
        interval_length: float = 100.0
    ):
        """
        Initialize simulation metrics.

        Args:
            num_stations: Number of stations
            num_classes: Number of job classes
            servers_per_station: Number of servers at each station
            interval_length: Length of observation intervals for MSER-5
        """
        self.num_stations = num_stations
        self.num_classes = num_classes
        self.station_metrics = [
            StationMetrics(num_classes, servers_per_station[i])
            for i in range(num_stations)
        ]
        self.start_time = 0.0
        self.warmup_time = 0.0
        self.interval_length = interval_length
        self.last_interval_time = 0.0

        # Global counters
        self.total_arrivals = 0
        self.total_departures = 0

    def record_arrival(
        self,
        station_id: int,
        class_id: int,
        current_time: float
    ) -> None:
        """Record an arrival at a station."""
        self.station_metrics[station_id].record_arrival(class_id, current_time)
        self.total_arrivals += 1

        # Check if we need to close an interval
        self._check_interval(current_time)

    def record_departure(
        self,
        station_id: int,
        class_id: int,
        response_time: float,
        current_time: float
    ) -> None:
        """Record a departure from a station."""
        self.station_metrics[station_id].record_departure(
            class_id, response_time, current_time
        )
        self.total_departures += 1

    def update_queue_length(
        self,
        station_id: int,
        class_id: int,
        queue_length: int,
        current_time: float
    ) -> None:
        """Update queue length at a station."""
        self.station_metrics[station_id].update_queue_length(
            class_id, queue_length, current_time
        )

    def update_busy_servers(
        self,
        station_id: int,
        class_id: int,
        busy_count: int,
        current_time: float
    ) -> None:
        """Update busy server count at a station."""
        self.station_metrics[station_id].update_busy_servers(
            class_id, busy_count, current_time
        )

    def _check_interval(self, current_time: float) -> None:
        """Check if current interval should be closed."""
        while current_time >= self.last_interval_time + self.interval_length:
            self.last_interval_time += self.interval_length
            for sm in self.station_metrics:
                sm.close_intervals(self.last_interval_time)

    def set_warmup_time(self, warmup_time: float) -> None:
        """
        Set the warmup time and truncate statistics.

        Args:
            warmup_time: Time when warmup period ended
        """
        self.warmup_time = warmup_time
        # Reset statistics starting from warmup time
        for sm in self.station_metrics:
            sm.reset(warmup_time)
        self.last_interval_time = warmup_time

    def get_results(
        self,
        current_time: float
    ) -> Dict[str, np.ndarray]:
        """
        Get simulation results as numpy arrays.

        Returns:
            Dictionary with keys: QN, UN, RN, TN (queue length, utilization,
            response time, throughput) - each is (num_stations, num_classes)
        """
        QN = np.zeros((self.num_stations, self.num_classes))
        UN = np.zeros((self.num_stations, self.num_classes))
        RN = np.zeros((self.num_stations, self.num_classes))
        TN = np.zeros((self.num_stations, self.num_classes))

        start = self.warmup_time if self.warmup_time > 0 else self.start_time

        for i, sm in enumerate(self.station_metrics):
            for k in range(self.num_classes):
                QN[i, k] = sm.get_queue_length(k, current_time)
                UN[i, k] = sm.get_utilization(k, current_time)
                RN[i, k] = sm.get_response_time(k)
                TN[i, k] = sm.get_throughput(k, current_time, start)

        return {
            'QN': QN,
            'UN': UN,
            'RN': RN,
            'TN': TN,
        }

    def reset(self, time: float = 0.0) -> None:
        """Reset all metrics."""
        self.start_time = time
        self.warmup_time = 0.0
        self.last_interval_time = time
        self.total_arrivals = 0
        self.total_departures = 0
        for sm in self.station_metrics:
            sm.reset(time)
