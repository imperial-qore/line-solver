"""
RewardState - Smart state accessor for reward functions in CTMC analysis.

Provides intuitive access to aggregated state populations at specific stations
and job classes, following Python conventions with snake_case method names.

Example:
    >>> # Access specific station/class population
    >>> pop = state.at(queue1, class1)
    >>>
    >>> # Access all classes at a station
    >>> view = state.at(queue1)
    >>> total = view.total()  # Sum across classes
    >>>
    >>> # Access class across all stations
    >>> view = state.for_class(class1)
    >>> total = view.total()  # Sum across stations
"""

__all__ = ['RewardState']


class RewardState:
    """Smart state accessor for reward functions.

    Provides intuitive access to aggregated state populations without
    requiring manual index calculations.
    """

    def __init__(self, state_vec, sn, nodes_to_station, classes_to_idx):
        """Initialize RewardState.

        Args:
            state_vec: Aggregated state vector (1D array)
            sn: NetworkStruct with model information
            nodes_to_station: Dict mapping node.index -> station index
            classes_to_idx: Dict mapping jobclass.index -> class index
        """
        self._state_vector = state_vec
        self._sn = sn
        self._nclasses = sn.nclasses
        self._nstations = sn.nstations
        self._node_to_station = nodes_to_station
        self._class_to_index = classes_to_idx

    def _get_node_index(self, node):
        """Get node index, handling both wrapper and native node implementations."""
        if hasattr(node, 'get_index'):
            return node.get_index()  # Native implementation
        elif hasattr(node, 'index'):
            return node.index  # Wrapper implementation
        else:
            raise ValueError(f"Node {node.name} has no index attribute")

    def _get_class_index(self, jobclass):
        """Get class index, handling both wrapper and native class implementations."""
        if hasattr(jobclass, 'get_index'):
            return jobclass.get_index()  # Native implementation
        elif hasattr(jobclass, 'index'):
            return jobclass.index  # Wrapper implementation
        else:
            raise ValueError(f"JobClass {jobclass.name} has no index attribute")

    def at(self, node, jobclass=None):
        """Access population at station or station+class.

        Args:
            node: Station node object
            jobclass: (Optional) Job class object

        Returns:
            If jobclass provided: float - population at node for this class
            If jobclass not provided: RewardStateView - view of all classes at node

        Examples:
            >>> view = state.at(queue1)              # All classes at queue1
            >>> pop = state.at(queue1, class1)       # Class1 jobs at queue1
            >>> total = state.at(queue1).total()     # Total jobs at queue1
        """
        node_idx = self._get_node_index(node)

        if jobclass is None:
            # Return view for all classes at this station
            station_idx = self._node_to_station.get(node_idx)
            if station_idx is None:
                raise ValueError(f"Node {node.name} not found in station mapping")

            start_idx = (station_idx - 1) * self._nclasses
            end_idx = station_idx * self._nclasses

            subvector = self._state_vector[start_idx:end_idx]
            from .reward_state_view import RewardStateView
            return RewardStateView(subvector, self, station_idx, -1)
        else:
            # Return scalar for specific station/class
            station_idx = self._node_to_station.get(node_idx)
            if station_idx is None:
                raise ValueError(f"Node {node.name} not found in station mapping")

            class_idx_key = self._get_class_index(jobclass)
            class_idx = self._class_to_index.get(class_idx_key)
            if class_idx is None:
                raise ValueError(f"JobClass {jobclass.name} not found in class mapping")

            idx = (station_idx - 1) * self._nclasses + class_idx - 1  # 0-indexed

            if idx < 0 or idx >= len(self._state_vector):
                raise IndexError(
                    f"Invalid state index {idx} for station {station_idx}, "
                    f"class {jobclass.name}")

            return float(self._state_vector[idx])

    def for_class(self, jobclass):
        """Access class across all stations.

        Args:
            jobclass: Job class object

        Returns:
            RewardStateView showing population of this class at each station

        Examples:
            >>> view = state.for_class(class1)
            >>> total = view.total()  # Total class1 jobs in system
        """
        class_idx_key = self._get_class_index(jobclass)
        class_idx = self._class_to_index.get(class_idx_key)
        if class_idx is None:
            raise ValueError(f"JobClass {jobclass.name} not found in class mapping")

        values = []
        for ist in range(self._nstations):
            idx = ist * self._nclasses + class_idx - 1  # 0-indexed
            values.append(self._state_vector[idx])

        from .reward_state_view import RewardStateView
        return RewardStateView(values, self, -1, class_idx)

    def for_station(self, station):
        """Alias for at(station).

        Args:
            station: Station node object

        Returns:
            RewardStateView for all classes at this station
        """
        return self.at(station)

    @property
    def sn(self):
        """Get the underlying NetworkStruct (for advanced usage).

        Returns:
            NetworkStruct reference
        """
        return self._sn

    @property
    def nclasses(self):
        """Get the number of job classes.

        Returns:
            int: Number of classes
        """
        return self._nclasses

    @property
    def nstations(self):
        """Get the number of stations.

        Returns:
            int: Number of stations
        """
        return self._nstations

    def __repr__(self):
        """String representation of RewardState."""
        return (f"RewardState(nstations={self._nstations}, "
                f"nclasses={self._nclasses})")
