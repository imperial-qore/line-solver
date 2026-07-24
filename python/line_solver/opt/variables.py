"""
Decision Variables for Queueing Network Optimization

This module provides decision variable classes for optimizing LINE network models.
Each variable type represents a different aspect of the network that can be tuned.

Key Classes:
    - DecisionVariable: Abstract base class for all variables
    - ServerAllocation: Optimize number of servers at a Station
    - StationReplicas: Optimize number of identical Station copies
    - RoutingProbabilities: Optimize routing of a JobClass
    - ClassServiceMapping: Optimize class-to-station service mapping
    - ServiceRate: Optimize processing rate of a Station
    - JobPopulation: Optimize number of jobs in a ClosedClass
    - ClassPriority: Optimize priority ordering of JobClasses

Example:
    >>> from line_solver import ServerAllocation, ServiceRate
    >>> var1 = ServerAllocation(queue, bounds=(1, 10))
    >>> var2 = ServiceRate(queue, jobclass, bounds=(0.5, 5.0))
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Tuple, List, Optional, Dict, Any, Union
from dataclasses import dataclass, field


def _resolveClass(model: 'Network', jobclass: 'JobClass') -> Optional['JobClass']:
    """
    Resolve a job class inside a (possibly copied) model by name.

    Decision variables hold references to the base model's objects; models
    are deep-copied per evaluation, so objects must be re-resolved by name
    in the target model or changes silently miss the copy.
    """
    name = jobclass.getName()
    getter = getattr(model, 'get_class_by_name', None)
    if callable(getter):
        found = getter(name)
        if found is not None:
            return found
    # Backend-agnostic fallback: the JAR-backed wrapper Network has no
    # get_class_by_name, so resolve by name over the class list.
    for candidate in model.getClasses():
        if candidate.getName() == name:
            return candidate
    return None


def _resolveNode(model: 'Network', node: 'Node') -> Optional['Node']:
    """Resolve a node inside a (possibly copied) model by name."""
    if hasattr(model, 'get_node_by_name'):
        found = model.get_node_by_name(node.getName())
        if found is not None:
            return found
    for candidate in model.getNodes():
        if candidate.getName() == node.getName():
            return candidate
    return None


def _connectionMatrix(model: 'Network') -> np.ndarray:
    """
    Node-by-node binary connection (adjacency) matrix, backend-agnostically.

    Native Network exposes get_connection_matrix(); the JAR-backed wrapper
    Network does not, so fall back to the network struct's connmatrix, which
    uses the same node ordering as getNodes().
    """
    getter = getattr(model, 'get_connection_matrix', None)
    if callable(getter):
        conn = getter()
        if conn is not None:
            return np.asarray(conn, dtype=float)
    sn = model.getStruct()
    return np.asarray(sn.connmatrix, dtype=float)


class DecisionVariable(ABC):
    """
    Abstract base class for optimization decision variables.

    Each decision variable represents a tunable parameter in the LINE network
    model. Variables are encoded as continuous values in [0, 1] for compatibility
    with differential evolution, then decoded to their actual domain.

    Attributes:
        name: Unique identifier for this variable
        dimension: Number of continuous values needed to encode this variable
    """

    def __init__(self, name: str):
        """
        Initialize decision variable.

        Args:
            name: Unique identifier for this variable
        """
        self._name = name
        self._dimension = 1

    def getName(self) -> str:
        """Get the variable name."""
        return self._name

    def getDimension(self) -> int:
        """Get the number of continuous values needed to encode this variable."""
        return self._dimension

    @abstractmethod
    def getBounds(self) -> List[Tuple[float, float]]:
        """
        Get bounds for each dimension of the encoded variable.

        Returns:
            List of (lower, upper) bounds for each encoded dimension.
            For standard variables, this is [(0, 1)] * dimension.
        """
        raise NotImplementedError

    @abstractmethod
    def decode(self, x: np.ndarray) -> Any:
        """
        Decode continuous values to the actual variable value.

        Args:
            x: Array of continuous values in [0, 1] of length self.dimension

        Returns:
            The decoded variable value in its native domain
        """
        raise NotImplementedError

    @abstractmethod
    def apply(self, model: 'Network', value: Any) -> None:
        """
        Apply the decoded value to a LINE Network model.

        Args:
            model: LINE Network model to modify
            value: Decoded variable value to apply
        """
        raise NotImplementedError

    @abstractmethod
    def getVariableType(self) -> str:
        """Get the type identifier for this variable (for decomposition)."""
        raise NotImplementedError

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        """LQN layer name(s) this variable perturbs, or None for flat models.

        Used by layer freezing (both explicit ``frozen_layers`` and adaptive
        ``auto_freeze``): a variable whose layer set intersects the frozen set
        is held fixed. Flat-network variables have no layer and return None.
        LQN variable subclasses override this.
        """
        return None

    def get_layer(self, model: 'Network') -> Optional[List[str]]:
        return self.getLayer(model)

    def currentValue(self, model: 'Network') -> Optional[Any]:
        """The variable's current (decoded) value in the given model, or None.

        Used by layer freezing to hold a variable at the model's existing
        parameter value. Flat and non-introspectable variables return None (the
        freeze then simply drops the variable, leaving the model default).
        """
        return None

    def current_value(self, model: 'Network') -> Optional[Any]:
        return self.currentValue(model)

    # Snake case aliases (use lambdas to avoid abstract method issues)
    def get_name(self) -> str:
        return self.getName()

    def get_dimension(self) -> int:
        return self.getDimension()

    def get_bounds(self):
        return self.getBounds()

    def get_variable_type(self) -> str:
        return self.getVariableType()

    # Properties
    @property
    def name(self) -> str:
        return self._name

    @property
    def dimension(self) -> int:
        return self._dimension


class ServerAllocation(DecisionVariable):
    """
    Optimize number of servers at a Station.

    Encodes an integer number of servers in the range [min_servers, max_servers].
    Uses continuous encoding that is rounded to nearest integer.

    Args:
        station: LINE Station (Queue) to configure
        bounds: Tuple of (min_servers, max_servers)
        name: Optional variable name (defaults to station name + '_servers')

    Example:
        >>> var = ServerAllocation(queue, bounds=(1, 10))
        >>> value = var.decode(np.array([0.5]))  # Returns ~5
    """

    def __init__(self, station: 'Station', bounds: Tuple[int, int],
                 name: str = None):
        self._station = station
        self._min_servers = int(bounds[0])
        self._max_servers = int(bounds[1])

        if name is None:
            name = f"{station.getName()}_servers"
        super().__init__(name)

    def getStation(self) -> 'Station':
        """Get the station this variable applies to."""
        return self._station

    def getMinServers(self) -> int:
        """Get minimum number of servers."""
        return self._min_servers

    def getMaxServers(self) -> int:
        """Get maximum number of servers."""
        return self._max_servers

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        """
        Decode continuous value to integer number of servers.

        Args:
            x: Array with single value in [0, 1]

        Returns:
            Integer number of servers in [min_servers, max_servers]
        """
        # Linear interpolation and round
        continuous = self._min_servers + x[0] * (self._max_servers - self._min_servers)
        return int(np.clip(np.round(continuous), self._min_servers, self._max_servers))

    def apply(self, model: 'Network', value: int) -> None:
        """
        Apply server count to the station.

        Args:
            model: LINE Network model
            value: Number of servers to set
        """
        # Find the station in the model by name
        for node in model.getNodes():
            if node.getName() == self._station.getName():
                node.setNumberOfServers(value)
                return

    def getVariableType(self) -> str:
        return 'server_allocation'

    # Snake case aliases
    get_station = getStation
    get_min_servers = getMinServers
    get_max_servers = getMaxServers


class StationReplicas(DecisionVariable):
    """
    Optimize number of identical copies (replicas) of a Station.

    This represents horizontal scaling where multiple identical stations
    handle the load. Routing to replicas is assumed to be load-balanced.

    Note: Applying replicas requires modifying the network topology,
    which is more complex than other variables.

    Args:
        station: LINE Station to replicate
        bounds: Tuple of (min_replicas, max_replicas)
        name: Optional variable name

    Example:
        >>> var = StationReplicas(queue, bounds=(1, 5))
    """

    def __init__(self, station: 'Station', bounds: Tuple[int, int],
                 name: str = None):
        self._station = station
        self._min_replicas = int(bounds[0])
        self._max_replicas = int(bounds[1])

        if name is None:
            name = f"{station.getName()}_replicas"
        super().__init__(name)

    def getStation(self) -> 'Station':
        """Get the station this variable applies to."""
        return self._station

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        """Decode to integer number of replicas."""
        continuous = self._min_replicas + x[0] * (self._max_replicas - self._min_replicas)
        return int(np.clip(np.round(continuous), self._min_replicas, self._max_replicas))

    def apply(self, model: 'Network', value: int) -> None:
        """
        Apply the replica count to the station.

        N replicas of a load-balanced identical station are represented as a
        single multiserver station with N times the station's configured base
        server count (an M/M/c-equivalent reduction). This keeps the network
        topology and node names fixed so results stay keyed by the original
        station, at the cost of pooling the replicas' queues into one.

        Args:
            model: LINE Network model (a per-evaluation copy)
            value: Number of replicas to represent
        """
        for node in model.getNodes():
            if node.getName() != self._station.getName():
                continue
            try:
                base_servers = int(node.getNumberOfServers())
            except (TypeError, ValueError, OverflowError):
                # Non-integer (e.g. infinite-server) base: fall back to 1.
                base_servers = 1
            if base_servers < 1:
                base_servers = 1
            node.setNumberOfServers(base_servers * int(value))
            return

    def getVariableType(self) -> str:
        return 'station_replicas'

    get_station = getStation


class RoutingProbabilities(DecisionVariable):
    """
    Optimize routing probabilities for a JobClass.

    Encodes routing probabilities from a source node to multiple target nodes.
    Uses simplex encoding to ensure probabilities sum to 1.

    Args:
        jobclass: LINE JobClass to configure routing for
        source: Source node for routing
        targets: List of target nodes
        name: Optional variable name

    Example:
        >>> var = RoutingProbabilities(jobclass, source=queue1, targets=[queue2, queue3])
    """

    def __init__(self, jobclass: 'JobClass', source: 'Node',
                 targets: List['Node'], name: str = None):
        self._jobclass = jobclass
        self._source = source
        self._targets = targets

        if name is None:
            name = f"{jobclass.getName()}_routing_from_{source.getName()}"
        super().__init__(name)

        # Dimension is len(targets) - 1 because probabilities sum to 1
        self._dimension = max(1, len(targets) - 1)

    def getJobClass(self) -> 'JobClass':
        """Get the job class this routing applies to."""
        return self._jobclass

    def getSource(self) -> 'Node':
        """Get the source node."""
        return self._source

    def getTargets(self) -> List['Node']:
        """Get the target nodes."""
        return self._targets

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)] * self._dimension

    def decode(self, x: np.ndarray) -> np.ndarray:
        """
        Decode to probability vector using stick-breaking.

        Args:
            x: Array of values in [0, 1]

        Returns:
            Array of probabilities summing to 1
        """
        n = len(self._targets)
        if n == 1:
            return np.array([1.0])

        # Stick-breaking construction
        probs = np.zeros(n)
        remaining = 1.0

        for i in range(n - 1):
            probs[i] = remaining * x[i]
            remaining -= probs[i]

        probs[n - 1] = remaining
        return probs

    def apply(self, model: 'Network', value: np.ndarray) -> None:
        """
        Apply routing probabilities to the model.

        Replaces the outgoing routing of the source node for this job class
        with the given probabilities, while preserving all other routes. If
        the model has a reified routing matrix (e.g. from a previous apply),
        its routes are carried over; otherwise default routing (uniform split
        over physical links) is reconstructed for every class from the
        connection matrix, reproducing LINE's link() semantics.

        Args:
            model: LINE Network model (a per-evaluation copy)
            value: Array of routing probabilities (one per target)
        """
        jobclass = _resolveClass(model, self._jobclass)
        if jobclass is None:
            return
        source = _resolveNode(model, self._source)
        if source is None:
            return

        src_name = source.getName()
        cls_name = jobclass.getName()

        rt = model.initRoutingMatrix()

        existing = None
        if hasattr(model, 'getRoutingMatrix'):
            existing = model.getRoutingMatrix()

        if existing is not None and hasattr(existing, '_routes'):
            # Carry over all routes except this (class, source) row
            for (cls_src, cls_dst), routes in existing._routes.items():
                for (node_src, node_dst), prob in routes.items():
                    if (cls_src.getName() == cls_name
                            and cls_dst.getName() == cls_name
                            and node_src.getName() == src_name):
                        continue
                    rt.set(cls_src, cls_dst, node_src, node_dst, prob)
        else:
            # Rebuild default routing for all classes from the topology,
            # skipping the row that this variable overrides
            nodes = model.getNodes()
            n = len(nodes)
            conn = _connectionMatrix(model)
            for cls in model.getClasses():
                for i in range(n):
                    if (cls.getName() == cls_name
                            and nodes[i].getName() == src_name):
                        continue
                    out_degree = conn[i].sum()
                    if out_degree <= 0:
                        continue
                    for j in range(n):
                        if conn[i, j] > 0:
                            rt.set(cls, cls, nodes[i], nodes[j],
                                   float(conn[i, j] / out_degree))

        # Set the overridden row
        for i, target in enumerate(self._targets):
            target_node = _resolveNode(model, target)
            if target_node is not None:
                rt.set(jobclass, jobclass, source, target_node,
                       float(value[i]))

        model.link(rt)

    def getVariableType(self) -> str:
        return 'routing'

    get_jobclass = getJobClass
    get_source = getSource
    get_targets = getTargets


class ClassServiceMapping(DecisionVariable):
    """
    Optimize mapping of a JobClass service to different Stations.

    This represents selecting which station(s) serve a particular job class,
    useful for workload placement decisions.

    Args:
        jobclass: LINE JobClass to configure
        stations: List of candidate stations
        name: Optional variable name

    Example:
        >>> var = ClassServiceMapping(jobclass, stations=[queue1, queue2, queue3])
    """

    def __init__(self, jobclass: 'JobClass', stations: List['Station'],
                 name: str = None):
        self._jobclass = jobclass
        self._stations = stations

        if name is None:
            name = f"{jobclass.getName()}_mapping"
        super().__init__(name)

    def getJobClass(self) -> 'JobClass':
        """Get the job class."""
        return self._jobclass

    def getStations(self) -> List['Station']:
        """Get candidate stations."""
        return self._stations

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        """
        Decode to station index.

        Args:
            x: Single value in [0, 1]

        Returns:
            Index of selected station
        """
        n = len(self._stations)
        idx = int(np.floor(x[0] * n))
        return min(idx, n - 1)

    def apply(self, model: 'Network', value: int) -> None:
        """
        Apply the class-to-station mapping by rerouting the job class.

        Routes all traffic of this class through the selected candidate station
        and bypasses the other candidates, while preserving the default routing
        of every other class. Per-node routing is reconstructed from the
        network's physical connection matrix as a uniform split over outgoing
        links, reproducing LINE's default link() semantics; the mapped class's
        adjacency has the non-selected candidate columns zeroed so no traffic
        enters them.

        Assumes the candidate stations are parallel alternatives sharing common
        predecessor/successor nodes (the standard placement topology) and that
        non-mapped classes use default routing. If both a ClassServiceMapping
        and a RoutingProbabilities variable target overlapping nodes, apply the
        routing variable last.

        Args:
            model: LINE Network model (a per-evaluation copy)
            value: Index of the selected station in the candidate list
        """
        if not self._stations:
            return
        value = int(np.clip(value, 0, len(self._stations) - 1))

        nodes = model.getNodes()
        index_of = {node.getName(): i for i, node in enumerate(nodes)}
        n = len(nodes)

        conn = _connectionMatrix(model)

        candidate_names = [s.getName() for s in self._stations]
        selected_name = candidate_names[value]
        # Non-selected candidates must receive none of this class's traffic.
        blocked = {index_of[name] for name in candidate_names
                   if name != selected_name and name in index_of}

        mapped = _resolveClass(model, self._jobclass)
        if mapped is None:
            return

        rt = model.initRoutingMatrix()
        for jobclass in model.getClasses():
            adj = conn.copy()
            if jobclass.getName() == mapped.getName():
                for b in blocked:
                    adj[:, b] = 0.0
            for i in range(n):
                out_degree = adj[i].sum()
                if out_degree <= 0:
                    continue
                for j in range(n):
                    if adj[i, j] > 0:
                        rt.set(jobclass, jobclass, nodes[i], nodes[j],
                               float(adj[i, j] / out_degree))
        model.link(rt)

    def getVariableType(self) -> str:
        return 'class_mapping'

    get_jobclass = getJobClass
    get_stations = getStations


class ServiceRate(DecisionVariable):
    """
    Optimize processing rate of a Station for a JobClass.

    Encodes a continuous service rate in [min_rate, max_rate].
    The rate is applied as an exponential service time distribution.

    Args:
        station: LINE Station to configure
        jobclass: JobClass to configure service for
        bounds: Tuple of (min_rate, max_rate)
        name: Optional variable name

    Example:
        >>> var = ServiceRate(queue, jobclass, bounds=(0.5, 5.0))
    """

    def __init__(self, station: 'Station', jobclass: 'JobClass',
                 bounds: Tuple[float, float], name: str = None):
        self._station = station
        self._jobclass = jobclass
        self._min_rate = float(bounds[0])
        self._max_rate = float(bounds[1])

        if name is None:
            name = f"{station.getName()}_{jobclass.getName()}_rate"
        super().__init__(name)

    def getStation(self) -> 'Station':
        """Get the station."""
        return self._station

    def getJobClass(self) -> 'JobClass':
        """Get the job class."""
        return self._jobclass

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> float:
        """
        Decode to service rate.

        Args:
            x: Single value in [0, 1]

        Returns:
            Service rate in [min_rate, max_rate]
        """
        return self._min_rate + x[0] * (self._max_rate - self._min_rate)

    def apply(self, model: 'Network', value: float) -> None:
        """
        Apply service rate to the station.

        Uses exponential distribution with the specified rate. The job class
        is resolved by name in the target model so the change lands on the
        per-evaluation copy rather than on the base model's class object.
        """
        # Import here to avoid circular dependency
        try:
            from line_solver import Exp
        except ImportError:
            # Fallback for when LINE is not available
            return

        jobclass = _resolveClass(model, self._jobclass)
        if jobclass is None:
            jobclass = self._jobclass

        for node in model.getNodes():
            if node.getName() == self._station.getName():
                node.setService(jobclass, Exp(value))
                return

    def getVariableType(self) -> str:
        return 'service_rate'

    def paramKey(self) -> tuple:
        """Sensitivity parameter key this variable controls: ('rate', st, cl)."""
        return ('rate', self._station.getName(), self._jobclass.getName())

    def decodeJacobian(self, x: np.ndarray) -> float:
        """d(decoded rate)/d(encoded x): constant slope of the linear map."""
        return self._max_rate - self._min_rate

    decode_jacobian = decodeJacobian
    param_key = paramKey

    get_station = getStation
    get_jobclass = getJobClass


class JobPopulation(DecisionVariable):
    """
    Optimize number of jobs in a ClosedClass.

    Encodes an integer job population in [min_jobs, max_jobs].

    Args:
        jobclass: LINE ClosedClass to configure
        bounds: Tuple of (min_jobs, max_jobs)
        name: Optional variable name

    Example:
        >>> var = JobPopulation(closed_class, bounds=(1, 100))
    """

    def __init__(self, jobclass: 'ClosedClass', bounds: Tuple[int, int],
                 name: str = None):
        self._jobclass = jobclass
        self._min_jobs = int(bounds[0])
        self._max_jobs = int(bounds[1])

        if name is None:
            name = f"{jobclass.getName()}_population"
        super().__init__(name)

    def getJobClass(self) -> 'ClosedClass':
        """Get the closed class."""
        return self._jobclass

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        """Decode to integer job count."""
        continuous = self._min_jobs + x[0] * (self._max_jobs - self._min_jobs)
        return int(np.clip(np.round(continuous), self._min_jobs, self._max_jobs))

    def apply(self, model: 'Network', value: int) -> None:
        """
        Apply the job population to this closed class within the given model.

        The class is looked up by name in the (copied) model so the base model
        held by the optimizer is never mutated across evaluations.

        Args:
            model: LINE Network model (a per-evaluation copy)
            value: Number of jobs to set for the closed class
        """
        cls = _resolveClass(model, self._jobclass)
        if cls is None:
            return
        if hasattr(cls, 'setNumberOfJobs'):
            cls.setNumberOfJobs(int(value))

    def getVariableType(self) -> str:
        return 'job_population'

    get_jobclass = getJobClass


class ClassPriority(DecisionVariable):
    """
    Optimize priority of JobClasses.

    Encodes priority as integer levels or as a permutation ordering.
    Higher priority values indicate higher priority.

    Args:
        jobclasses: List of JobClasses to prioritize
        mode: 'levels' for integer priorities, 'permutation' for ordering
        priority_range: Tuple of (min_priority, max_priority) for 'levels' mode
        name: Optional variable name

    Example:
        >>> var = ClassPriority([class1, class2, class3], mode='levels')
    """

    def __init__(self, jobclasses: List['JobClass'],
                 mode: str = 'levels',
                 priority_range: Tuple[int, int] = (1, 10),
                 name: str = None):
        self._jobclasses = jobclasses
        self._mode = mode
        self._min_priority = priority_range[0]
        self._max_priority = priority_range[1]

        if name is None:
            name = "class_priorities"
        super().__init__(name)

        # For 'levels': one dimension per class
        # For 'permutation': n-1 dimensions for n classes
        if mode == 'levels':
            self._dimension = len(jobclasses)
        else:
            self._dimension = max(1, len(jobclasses) - 1)

    def getJobClasses(self) -> List['JobClass']:
        """Get the job classes."""
        return self._jobclasses

    def getMode(self) -> str:
        """Get encoding mode."""
        return self._mode

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)] * self._dimension

    def decode(self, x: np.ndarray) -> Union[List[int], List[int]]:
        """
        Decode to priorities.

        For 'levels': returns list of integer priorities
        For 'permutation': returns ordered list of class indices
        """
        if self._mode == 'levels':
            priorities = []
            for val in x:
                p = self._min_priority + val * (self._max_priority - self._min_priority)
                priorities.append(int(np.round(p)))
            return priorities
        else:
            # Permutation from random keys
            n = len(self._jobclasses)
            if n == 1:
                return [0]
            # Extend x with 0 for last position
            keys = np.concatenate([x, [0.0]])
            return list(np.argsort(-keys))  # Descending order

    def apply(self, model: 'Network', value: Union[List[int], List[int]]) -> None:
        """
        Apply priorities to job classes.

        Classes are resolved by name in the target model so the change lands
        on the per-evaluation copy rather than on the base model's objects.
        """
        def target(jobclass):
            resolved = _resolveClass(model, jobclass)
            return resolved if resolved is not None else jobclass

        if self._mode == 'levels':
            for i, jobclass in enumerate(self._jobclasses):
                target(jobclass).setPriority(value[i])
        else:
            # Assign priorities based on permutation order
            for rank, class_idx in enumerate(value):
                priority = len(value) - rank  # Higher rank = higher priority
                target(self._jobclasses[class_idx]).setPriority(priority)

    def getVariableType(self) -> str:
        return 'class_priority'

    get_jobclasses = getJobClasses
    get_mode = getMode


# ===================================================================
# LayeredNetwork (LQN) decision variables
# ===================================================================
#
# These target LQN model parameters instead of flat-network constructs. Each
# resolves its element by name in the per-evaluation model copy (LQN elements
# carry a plain ``name`` attribute, handled by ``layered.elem_name``) and
# mutates it via the established LayeredNetwork setters. ``getLayer`` reports
# the layer(s) the variable perturbs, consumed by layer freezing. HostDemand
# additionally exposes the hooks the partial-sensitivity gradient needs.


class HostDemand(DecisionVariable):
    """Optimize the mean host demand of an LQN Activity (continuous).

    The host demand D is the mean service requirement the activity places on
    its processor; the processor-layer service rate is mu = 1/D. This is the
    primary LQN tuning knob (analogous to ServiceRate for a flat station).

    Args:
        activity: LQN Activity (or its name) to tune.
        bounds: (min_demand, max_demand) for the mean host demand.
        name: Optional variable name (defaults to '<activity>_hostdemand').

    Example:
        >>> var = HostDemand(activity_AS1, bounds=(0.02, 0.2))
    """

    def __init__(self, activity, bounds, name: str = None):
        from .layered import elem_name
        self._activity = elem_name(activity) if not isinstance(activity, str) \
            else activity
        self._min_demand = float(bounds[0])
        self._max_demand = float(bounds[1])
        if name is None:
            name = f"{self._activity}_hostdemand"
        super().__init__(name)

    def getActivity(self) -> str:
        return self._activity

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> float:
        """Decode to a mean host demand in [min_demand, max_demand]."""
        return self._min_demand + x[0] * (self._max_demand - self._min_demand)

    def apply(self, model: 'Network', value: float) -> None:
        """Set the activity's mean host demand on the (copied) LQN model."""
        from .layered import resolve_activity
        act = resolve_activity(model, self._activity)
        if act is not None:
            act.setHostDemand(float(value))

    def getVariableType(self) -> str:
        return 'host_demand'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        from .layered import activity_processor_name
        proc = activity_processor_name(model, self._activity)
        return [proc] if proc is not None else None

    def currentValue(self, model: 'Network') -> Optional[float]:
        """The activity's current mean host demand in the given model."""
        from .layered import resolve_activity
        act = resolve_activity(model, self._activity)
        if act is not None and hasattr(act, 'getHostDemandMean'):
            return float(act.getHostDemandMean())
        return None

    current_value = currentValue

    # ---- partial-sensitivity gradient hooks --------------------------------

    def sensKey(self, model: 'Network') -> Optional[tuple]:
        """Row key (Station, JobClass) in the per-layer sensitivity table.

        The host-layer rows are (Layer=processor, Station=processor,
        JobClass=activity); the sensitivity value is d(metric)/d(service rate).
        """
        from .layered import activity_processor_name
        proc = activity_processor_name(model, self._activity)
        if proc is None:
            return None
        return (proc, self._activity)

    def sensMetricTargets(self, model: 'Network') -> Dict[str, tuple]:
        """Map each layer-row metric to the EvaluationResult key it approximates.

        The host-layer row's utilization tracks the processor node's
        utilization; its throughput/queue-length/response-time track the
        activity node's. Keys match the LQN EvaluationResult convention: Util
        by node name, the others by (node, node).
        """
        from .layered import activity_processor_name
        proc = activity_processor_name(model, self._activity)
        act = self._activity
        return {
            'Util': proc,
            'Tput': (act, act),
            'QLen': (act, act),
            'RespT': (act, act),
        }

    def rateJacobian(self, value: float) -> float:
        """d(service rate)/d(demand) = d(1/D)/dD = -1/D^2 at demand D=value."""
        v = float(value)
        if v <= 0:
            return 0.0
        return -1.0 / (v * v)

    def decodeJacobian(self, x: np.ndarray) -> float:
        """d(decoded demand)/d(encoded x): constant slope of the linear map."""
        return self._max_demand - self._min_demand

    get_activity = getActivity
    sens_key = sensKey
    sens_metric_targets = sensMetricTargets
    rate_jacobian = rateJacobian
    decode_jacobian = decodeJacobian


class ActivityThinkTime(DecisionVariable):
    """Optimize the activity-level think time of an LQN Activity (continuous)."""

    def __init__(self, activity, bounds, name: str = None):
        from .layered import elem_name
        self._activity = elem_name(activity) if not isinstance(activity, str) \
            else activity
        self._min_value = float(bounds[0])
        self._max_value = float(bounds[1])
        if name is None:
            name = f"{self._activity}_thinktime"
        super().__init__(name)

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> float:
        return self._min_value + x[0] * (self._max_value - self._min_value)

    def apply(self, model: 'Network', value: float) -> None:
        from .layered import resolve_activity
        act = resolve_activity(model, self._activity)
        if act is not None:
            act.setThinkTime(float(value))

    def getVariableType(self) -> str:
        return 'think_time'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        from .layered import activity_processor_name
        proc = activity_processor_name(model, self._activity)
        # the activity's own task layer is named after its task
        from .layered import resolve_activity, elem_name
        act = resolve_activity(model, self._activity)
        task = getattr(act, 'task', None) if act is not None else None
        layers = []
        if task is not None:
            layers.append(elem_name(task))
        if proc is not None:
            layers.append(proc)
        return layers or None

    def currentValue(self, model: 'Network') -> Optional[float]:
        from .layered import resolve_activity, dist_mean
        act = resolve_activity(model, self._activity)
        if act is None:
            return None
        return dist_mean(getattr(act, 'think_time', None))


class TaskThinkTime(DecisionVariable):
    """Optimize the think time of an LQN Task (continuous)."""

    def __init__(self, task, bounds, name: str = None):
        from .layered import elem_name
        self._task = elem_name(task) if not isinstance(task, str) else task
        self._min_value = float(bounds[0])
        self._max_value = float(bounds[1])
        if name is None:
            name = f"{self._task}_thinktime"
        super().__init__(name)

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> float:
        return self._min_value + x[0] * (self._max_value - self._min_value)

    def apply(self, model: 'Network', value: float) -> None:
        from .layered import resolve_task
        task = resolve_task(model, self._task)
        if task is not None:
            task.set_think_time(float(value))

    def getVariableType(self) -> str:
        return 'think_time'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        from .layered import task_processor_name
        layers = [self._task]
        proc = task_processor_name(model, self._task)
        if proc is not None:
            layers.append(proc)
        return layers

    def currentValue(self, model: 'Network') -> Optional[float]:
        from .layered import resolve_task, dist_mean
        task = resolve_task(model, self._task)
        if task is None:
            return None
        return dist_mean(getattr(task, 'think_time', None))


class TaskMultiplicity(DecisionVariable):
    """Optimize the multiplicity (thread/instance count) of an LQN Task (integer)."""

    def __init__(self, task, bounds, name: str = None):
        from .layered import elem_name
        self._task = elem_name(task) if not isinstance(task, str) else task
        self._min_value = int(bounds[0])
        self._max_value = int(bounds[1])
        if name is None:
            name = f"{self._task}_multiplicity"
        super().__init__(name)

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        continuous = self._min_value + x[0] * (self._max_value - self._min_value)
        return int(np.clip(np.round(continuous), self._min_value, self._max_value))

    def apply(self, model: 'Network', value: int) -> None:
        from .layered import resolve_task
        task = resolve_task(model, self._task)
        if task is not None:
            # Tasks store multiplicity as a plain attribute (no setter);
            # setting it directly mirrors how the LQN builder assigns it.
            task.multiplicity = int(value)

    def getVariableType(self) -> str:
        return 'task_multiplicity'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        from .layered import task_processor_name
        layers = [self._task]
        proc = task_processor_name(model, self._task)
        if proc is not None:
            layers.append(proc)
        return layers

    def currentValue(self, model: 'Network') -> Optional[int]:
        from .layered import resolve_task
        task = resolve_task(model, self._task)
        if task is None:
            return None
        mult = getattr(task, 'multiplicity', None)
        try:
            return int(mult)
        except (TypeError, ValueError):
            return None


class TaskReplication(DecisionVariable):
    """Optimize the replication (fan-out replicas) of an LQN Task (integer)."""

    def __init__(self, task, bounds, name: str = None):
        from .layered import elem_name
        self._task = elem_name(task) if not isinstance(task, str) else task
        self._min_value = int(bounds[0])
        self._max_value = int(bounds[1])
        if name is None:
            name = f"{self._task}_replication"
        super().__init__(name)

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        continuous = self._min_value + x[0] * (self._max_value - self._min_value)
        return int(np.clip(np.round(continuous), self._min_value, self._max_value))

    def apply(self, model: 'Network', value: int) -> None:
        from .layered import resolve_task
        task = resolve_task(model, self._task)
        if task is not None and hasattr(task, 'setReplication'):
            task.setReplication(int(value))

    def getVariableType(self) -> str:
        return 'task_replication'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        from .layered import task_processor_name
        layers = [self._task]
        proc = task_processor_name(model, self._task)
        if proc is not None:
            layers.append(proc)
        return layers

    def currentValue(self, model: 'Network') -> Optional[int]:
        from .layered import resolve_task
        task = resolve_task(model, self._task)
        if task is None:
            return None
        getter = getattr(task, 'getReplication', None)
        if callable(getter):
            try:
                return int(getter())
            except (TypeError, ValueError):
                return None
        return int(getattr(task, '_replication', 1))


class ProcessorMultiplicity(DecisionVariable):
    """Optimize the multiplicity (core count) of an LQN Processor (integer)."""

    def __init__(self, processor, bounds, name: str = None):
        from .layered import elem_name
        self._processor = elem_name(processor) if not isinstance(processor, str) \
            else processor
        self._min_value = int(bounds[0])
        self._max_value = int(bounds[1])
        if name is None:
            name = f"{self._processor}_multiplicity"
        super().__init__(name)

    def getBounds(self) -> List[Tuple[float, float]]:
        return [(0.0, 1.0)]

    def decode(self, x: np.ndarray) -> int:
        continuous = self._min_value + x[0] * (self._max_value - self._min_value)
        return int(np.clip(np.round(continuous), self._min_value, self._max_value))

    def apply(self, model: 'Network', value: int) -> None:
        from .layered import resolve_processor
        proc = resolve_processor(model, self._processor)
        if proc is not None:
            proc.multiplicity = int(value)

    def getVariableType(self) -> str:
        return 'processor_multiplicity'

    def getLayer(self, model: 'Network') -> Optional[List[str]]:
        # A processor owns its own host layer, named after the processor.
        return [self._processor]

    def currentValue(self, model: 'Network') -> Optional[int]:
        from .layered import resolve_processor
        proc = resolve_processor(model, self._processor)
        if proc is None:
            return None
        try:
            return int(getattr(proc, 'multiplicity', None))
        except (TypeError, ValueError):
            return None


# Public API
__all__ = [
    'DecisionVariable',
    'ServerAllocation',
    'StationReplicas',
    'RoutingProbabilities',
    'ClassServiceMapping',
    'ServiceRate',
    'JobPopulation',
    'ClassPriority',
    # LayeredNetwork (LQN) variables
    'HostDemand',
    'ActivityThinkTime',
    'TaskThinkTime',
    'TaskMultiplicity',
    'TaskReplication',
    'ProcessorMultiplicity',
]
