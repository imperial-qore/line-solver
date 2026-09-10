"""
RewardStateView - View of reward state subset with aggregation operations.

Provides aggregation operations (sum, max, min, count) on a subset of the
state vector. Returned by RewardState.at() and for_class().

Example:
    >>> # Station view (all classes)
    >>> view = state.at(queue1)
    >>> total = view.total()           # Total jobs at queue1
    >>> max_jobs = view.max()          # Max jobs across classes
    >>>
    >>> # Class view (all stations)
    >>> view = state.for_class(class1)
    >>> total = view.total()           # Total class1 jobs in system
"""

__all__ = ['RewardStateView']


class RewardStateView:
    """View of reward state subset with aggregation operations."""

    def __init__(self, values, parent_state, station_idx, class_idx):
        """Initialize RewardStateView.

        Args:
            values: List or array of values
            parent_state: The parent RewardState object
            station_idx: Station index if station view, -1 if class view
            class_idx: Class index if class view, -1 if station view
        """
        self._values = list(values)  # Ensure it's a list for consistency
        self._parent_state = parent_state
        self._station_idx = station_idx
        self._class_idx = class_idx

    def total(self):
        """Sum all values in this view.

        Returns:
            float: The sum of all values

        Example:
            >>> view = state.at(queue1)
            >>> total_jobs = view.total()
        """
        return sum(self._values)

    def max(self):
        """Maximum value in this view.

        Returns:
            float: The maximum value, or 0 if empty

        Example:
            >>> view = state.at(queue1)
            >>> max_class = view.max()
        """
        if not self._values:
            return 0.0
        return max(self._values)

    def min(self):
        """Minimum value in this view.

        Returns:
            float: The minimum value, or 0 if empty

        Example:
            >>> view = state.at(queue1)
            >>> min_class = view.min()
        """
        if not self._values:
            return 0.0
        return min(self._values)

    def count(self):
        """Count non-zero entries in this view.

        Returns:
            int: The number of non-zero values

        Example:
            >>> view = state.at(queue1)
            >>> nonempty_classes = view.count()
        """
        return sum(1 for v in self._values if v > 0)

    def at_class(self, jobclass):
        """Filter station view to specific class (scalar result).

        Only valid if this is a station view (created by at(station)).

        Args:
            jobclass: Job class object

        Returns:
            float: Population of this class at the station

        Raises:
            ValueError: If this is not a station view
        """
        if self._station_idx == -1:
            raise ValueError(
                "at_class() is only valid for station views "
                "(created by state.at(station)). "
                "Use state.at(station, class) directly instead.")
        return self._parent_state.at(self._station_idx, jobclass)

    def at_station(self, station):
        """Filter class view to specific station (scalar result).

        Only valid if this is a class view (created by for_class(class)).

        Args:
            station: Station node object

        Returns:
            float: Population of this class at the station

        Raises:
            ValueError: If this is not a class view
        """
        if self._class_idx == -1:
            raise ValueError(
                "at_station() is only valid for class views "
                "(created by state.for_class(class)). "
                "Use state.at(station, class) directly instead.")
        return self._parent_state.at(station, self._class_idx)

    def is_station_view(self):
        """Check if this is a station view.

        Returns:
            bool: True if station view, False if class view
        """
        return self._station_idx != -1

    def is_class_view(self):
        """Check if this is a class view.

        Returns:
            bool: True if class view, False if station view
        """
        return self._class_idx != -1

    # Support arithmetic operations for rewards
    def __add__(self, other):
        """Add to scalar or another view."""
        if isinstance(other, RewardStateView):
            return self.total() + other.total()
        return self.total() + other

    def __radd__(self, other):
        """Right add."""
        return other + self.total()

    def __sub__(self, other):
        """Subtract from scalar."""
        if isinstance(other, RewardStateView):
            return self.total() - other.total()
        return self.total() - other

    def __rsub__(self, other):
        """Right subtract."""
        return other - self.total()

    def __mul__(self, other):
        """Multiply by scalar."""
        if isinstance(other, RewardStateView):
            return sum(v1 * v2 for v1, v2 in zip(self._values, other._values))
        return self.total() * other

    def __rmul__(self, other):
        """Right multiply."""
        return self.total() * other

    def __truediv__(self, other):
        """Divide by scalar."""
        if other == 0:
            raise ValueError("Cannot divide by zero")
        return self.total() / other

    def __rtruediv__(self, other):
        """Right divide."""
        if self.total() == 0:
            raise ValueError("Cannot divide by zero")
        return other / self.total()

    def __ge__(self, threshold):
        """Greater than or equal (for conditional rewards)."""
        return float(self.total() >= threshold)

    def __le__(self, threshold):
        """Less than or equal (for conditional rewards)."""
        return float(self.total() <= threshold)

    def __gt__(self, threshold):
        """Greater than (for conditional rewards)."""
        return float(self.total() > threshold)

    def __lt__(self, threshold):
        """Less than (for conditional rewards)."""
        return float(self.total() < threshold)

    def __eq__(self, other):
        """Equality comparison."""
        if isinstance(other, RewardStateView):
            return self.total() == other.total()
        return self.total() == other

    def __repr__(self):
        """String representation."""
        values_str = ", ".join(f"{v:.2f}" for v in self._values)
        view_type = "Station" if self.is_station_view() else "Class"
        return (f"RewardStateView[{view_type}]({values_str}) "
                f"(total={self.total():.2f})")

    def __str__(self):
        """String representation."""
        return self.__repr__()
