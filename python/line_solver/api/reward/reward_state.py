"""
Reward state accessor for LINE networks (pure Python).

Provides pythonic access to network state for reward function evaluation.
Allows aggregation and filtering of state values by station and/or job class.
"""

import numpy as np
from typing import Union, Optional, Dict, Callable, Any
from ...lang.base import Station, JobClass


class RewardStateView:
    """
    View of reward state subset with aggregation operations.

    Provides aggregation operations (sum, max, min, count) on a subset of
    the state vector. Returned by RewardState.at() and forClass().

    Methods:
        total() - Sum all values in this view
        max() - Maximum value in this view
        min() - Minimum value in this view
        count() - Number of non-zero entries in this view
    """

    def __init__(self, values: np.ndarray, parent_state: 'RewardState',
                 station_idx: Optional[int] = None, class_idx: Optional[int] = None):
        """
        Initialize a RewardStateView.

        Args:
            values: Vector of state subset values
            parent_state: RewardState parent object
            station_idx: Station index if this is a station view, None otherwise
            class_idx: Class index if this is a class view, None otherwise
        """
        self.values = np.atleast_1d(values).astype(float)
        self.parent_state = parent_state
        self.station_idx = station_idx
        self.class_idx = class_idx

    def total(self) -> float:
        """Return sum of all values in this view."""
        return float(np.sum(self.values))

    def max(self) -> float:
        """Return maximum value in this view."""
        if len(self.values) == 0:
            return 0.0
        return float(np.max(self.values))

    def min(self) -> float:
        """Return minimum value in this view."""
        if len(self.values) == 0:
            return 0.0
        return float(np.min(self.values))

    def count(self) -> int:
        """Return number of non-zero entries in this view."""
        return int(np.sum(self.values > 0))

    def at_class(self, jobclass: JobClass) -> float:
        """
        Filter station view to specific class.

        Only valid if this is a station view (created by state.at(station)).

        Args:
            jobclass: Job class to filter by

        Returns:
            Population of jobclass at the station in this view
        """
        if self.station_idx is None:
            raise ValueError("at_class() is only valid for station views (created by state.at(station))")

        return self.parent_state.at(self.station_idx, jobclass)

    def at_station(self, station: Station) -> float:
        """
        Filter class view to specific station.

        Only valid if this is a class view (created by state.forClass(class)).

        Args:
            station: Station to filter by

        Returns:
            Population of the class in this view at the specified station
        """
        if self.class_idx is None:
            raise ValueError("at_station() is only valid for class views (created by state.forClass(class))")

        return self.parent_state.at(station, self.class_idx)

    # Arithmetic operations
    def __add__(self, other: Union[float, 'RewardStateView']) -> float:
        """Add scalar or compatible view."""
        if isinstance(other, RewardStateView):
            return self.total() + other.total()
        else:
            return self.total() + float(other)

    def __sub__(self, other: Union[float, 'RewardStateView']) -> float:
        """Subtract scalar or compatible view."""
        if isinstance(other, RewardStateView):
            return self.total() - other.total()
        else:
            return self.total() - float(other)

    def __mul__(self, other: Union[float, 'RewardStateView']) -> float:
        """Multiply view values by scalar or compatible view."""
        if isinstance(other, RewardStateView):
            return float(np.dot(self.values, other.values))
        else:
            return self.total() * float(other)

    def __rmul__(self, other: float) -> float:
        """Right multiply."""
        return self.total() * float(other)

    def __ge__(self, threshold: float) -> bool:
        """Greater than or equal (for conditional rewards)."""
        return self.total() >= threshold

    def __le__(self, threshold: float) -> bool:
        """Less than or equal (for conditional rewards)."""
        return self.total() <= threshold

    def __gt__(self, threshold: float) -> bool:
        """Greater than (for conditional rewards)."""
        return self.total() > threshold

    def __lt__(self, threshold: float) -> bool:
        """Less than (for conditional rewards)."""
        return self.total() < threshold

    def __repr__(self) -> str:
        """String representation."""
        return f"RewardStateView(total={self.total():.4f}, values={self.values})"


class RewardState:
    """
    Accessor for network state in reward function context.

    Provides pythonic access to network state for reward function evaluation.
    Allows aggregation and filtering of state values by station and/or job class.

    Usage:
        # Station view (all classes)
        view = state.at(queue1)
        total = view.total()  # Total jobs at queue1
        jobs_class1 = view.at_class(class1)  # Class1 jobs at queue1

        # Class view (all stations)
        view = state.for_class(class1)
        total = view.total()  # Total class1 jobs
        jobs_queue1 = view.at_station(queue1)  # Class1 jobs at queue1
    """

    def __init__(self, state_matrix: np.ndarray, stations: list = None,
                 job_classes: list = None):
        """
        Initialize RewardState.

        Args:
            state_matrix: State matrix (stations x classes)
            stations: List of Station objects (optional)
            job_classes: List of JobClass objects (optional)
        """
        self.state_matrix = np.atleast_2d(state_matrix).astype(float)
        self.stations = stations if stations is not None else []
        self.job_classes = job_classes if job_classes is not None else []

        # Validate dimensions
        if self.state_matrix.shape[0] != len(self.stations) and len(self.stations) > 0:
            raise ValueError(f"State matrix rows ({self.state_matrix.shape[0]}) must match number of stations ({len(self.stations)})")

        if self.state_matrix.shape[1] != len(self.job_classes) and len(self.job_classes) > 0:
            raise ValueError(f"State matrix columns ({self.state_matrix.shape[1]}) must match number of classes ({len(self.job_classes)})")

    def at(self, station_or_index: Union[Station, int], jobclass: Optional[JobClass] = None) -> Union[float, RewardStateView]:
        """
        Access state at specific station and/or class.

        Args:
            station_or_index: Station object or index
            jobclass: Job class (optional, defaults to all classes)

        Returns:
            If jobclass specified: scalar population count
            If jobclass not specified: RewardStateView of all classes at station
        """
        # Resolve station index
        if isinstance(station_or_index, int):
            st_idx = station_or_index
        else:
            # Find station index
            st_idx = None
            if hasattr(station_or_index, '_index'):
                st_idx = station_or_index._index
            else:
                for i, s in enumerate(self.stations):
                    if s is station_or_index:
                        st_idx = i
                        break

            if st_idx is None:
                raise ValueError(f"Station {station_or_index} not found in state")

        # Validate station index
        if st_idx < 0 or st_idx >= self.state_matrix.shape[0]:
            raise IndexError(f"Station index {st_idx} out of bounds")

        # If jobclass specified, return scalar
        if jobclass is not None:
            # Resolve class index
            if isinstance(jobclass, int):
                cl_idx = jobclass
            else:
                # Find class index
                cl_idx = None
                if hasattr(jobclass, '_index'):
                    cl_idx = jobclass._index
                else:
                    for i, c in enumerate(self.job_classes):
                        if c is jobclass:
                            cl_idx = i
                            break

                if cl_idx is None:
                    raise ValueError(f"Class {jobclass} not found in state")

            if cl_idx < 0 or cl_idx >= self.state_matrix.shape[1]:
                raise IndexError(f"Class index {cl_idx} out of bounds")

            return float(self.state_matrix[st_idx, cl_idx])

        # Return view of all classes at this station
        return RewardStateView(self.state_matrix[st_idx, :], self, station_idx=st_idx)

    def for_class(self, jobclass: Union[JobClass, int]) -> RewardStateView:
        """
        Access state for specific job class across all stations.

        Args:
            jobclass: Job class object or index

        Returns:
            RewardStateView of this class across all stations
        """
        # Resolve class index
        if isinstance(jobclass, int):
            cl_idx = jobclass
        else:
            # Find class index
            cl_idx = None
            if hasattr(jobclass, '_index'):
                cl_idx = jobclass._index
            else:
                for i, c in enumerate(self.job_classes):
                    if c is jobclass:
                        cl_idx = i
                        break

            if cl_idx is None:
                raise ValueError(f"Class {jobclass} not found in state")

        if cl_idx < 0 or cl_idx >= self.state_matrix.shape[1]:
            raise IndexError(f"Class index {cl_idx} out of bounds")

        return RewardStateView(self.state_matrix[:, cl_idx], self, class_idx=cl_idx)

    def total(self) -> float:
        """Return total number of jobs in system."""
        return float(np.sum(self.state_matrix))

    def __repr__(self) -> str:
        """String representation."""
        return f"RewardState(shape={self.state_matrix.shape}, total={self.total():.0f})"

    def __str__(self) -> str:
        """Readable representation."""
        lines = ["RewardState:"]
        for i, st in enumerate(self.stations):
            st_name = st.name if hasattr(st, 'name') else f"Station{i}"
            lines.append(f"  {st_name}: {self.state_matrix[i, :]}")
        return "\n".join(lines)
