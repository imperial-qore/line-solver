"""
Cluster builder.

Mirrors the Java :class:`jline.gen.Cluster` builder: chainable setters that
configure a Source -> Dispatcher -> Servers -> Sink (or Think -> Dispatcher ->
Servers -> Think) topology, plus comparison and parametric-sweep helpers.
"""

from collections import OrderedDict
from typing import Callable, List, Optional, Sequence, Union

import numpy as np

from ..lang.base import RoutingStrategy, SchedStrategy
from ..lang.network import Network


class Cluster:
    """
    Builder for cluster models.

    Examples:
        >>> from line_solver import Cluster, SolverMVA
        >>> cluster = (Cluster()
        ...         .set_num_stations(4)
        ...         .set_arrival_rate(1.0)
        ...         .set_service_rate(0.4)
        ...         .set_scheduling(SchedStrategy.PS)
        ...         .set_dispatching(RoutingStrategy.RAND))
        >>> model = cluster.build()
        >>> SolverMVA(model).get_avg_table()
    """

    def __init__(self):
        """Default cluster: 2 servers, single open class with arrival rate 1.0,
        service rate 1.0 at each server, PS scheduling, RAND dispatching.
        Configure further via the chainable ``set_*`` methods.
        """
        self._num_stations: int = 2
        self._arrival_rates: List[float] = [1.0]
        self._service_rates: np.ndarray = np.full((2, 1), 1.0)
        self._station_counts: List[int] = [1, 1]
        self._scheduling: SchedStrategy = SchedStrategy.PS
        self._dispatching: RoutingStrategy = RoutingStrategy.RAND
        self._closed: bool = False
        self._population: Optional[List[int]] = None
        self._think_times: Optional[List[float]] = None
        self._dispatch_probs: Optional[np.ndarray] = None  # PROB: rows=classes (1 broadcasts), cols=servers
        self._dispatch_weights: Optional[List[int]] = None  # WRROBIN: per-server integer weights
        self._sq_k: Optional[int] = None
        self._arrival_scvs: Optional[np.ndarray] = None    # per-class arrival SCVs
        self._service_scvs: Optional[np.ndarray] = None    # (M, R) service SCVs

    # ------------------------------------------------------------------ rates
    def set_num_stations(self, M: int) -> 'Cluster':
        """Number of parallel server queues. Replicates the current single-class
        service rate across all servers and resets per-server multiplicity to all-1.
        """
        if M <= 0:
            raise ValueError("num_stations must be positive")
        sample = float(self._service_rates[0, 0]) if self._service_rates.size else 1.0
        self._num_stations = M
        self._service_rates = np.full((M, 1), sample)
        self._station_counts = [1] * M
        return self

    def set_arrival_rate(self, lam: float) -> 'Cluster':
        """Single-class arrival rate."""
        if lam <= 0:
            raise ValueError("arrival rate must be positive")
        self._arrival_rates = [float(lam)]
        return self

    def set_arrival_rates(self, lambdas: Sequence[float]) -> 'Cluster':
        """Per-class arrival rates (length R)."""
        clean = [float(x) for x in lambdas]
        if any(x <= 0 for x in clean):
            raise ValueError("arrival rates must be positive")
        self._arrival_rates = clean
        return self

    def set_service_rate(self, mu: float) -> 'Cluster':
        """Single service rate, broadcast across all (server, class) pairs."""
        if mu <= 0:
            raise ValueError("service rate must be positive")
        self._service_rates = np.full((self._num_stations, 1), float(mu))
        return self

    def set_service_rates(self, rates) -> 'Cluster':
        """Per-(server, class) service rates as a ``(num_stations, R)`` matrix."""
        arr = np.atleast_2d(np.asarray(rates, dtype=float))
        if arr.shape[0] != self._num_stations:
            raise ValueError("service_rates outer dim must equal num_stations")
        if np.any(arr <= 0):
            raise ValueError("service rates must be positive")
        self._service_rates = arr.copy()
        return self

    # ------------------------------------------------------------------ setters

    def set_dispatching(self, dispatching: RoutingStrategy) -> 'Cluster':
        self._dispatching = dispatching
        return self

    def set_scheduling(self, scheduling: SchedStrategy) -> 'Cluster':
        self._scheduling = scheduling
        return self

    def set_station_servers(self, counts: Sequence[int]) -> 'Cluster':
        if len(counts) != self._num_stations:
            raise ValueError("counts length must equal num_stations")
        self._station_counts = [int(c) for c in counts]
        return self

    def set_probabilities(self, probs) -> 'Cluster':
        """PROB dispatching with per-server probabilities.

        ``probs`` is either a 1-D sequence of length ``num_stations`` (broadcast
        to every class) or a 2-D ``(R, num_stations)`` array.
        """
        arr = np.asarray(probs, dtype=float)
        if arr.ndim == 1:
            if arr.shape[0] != self._num_stations:
                raise ValueError("probs length must equal num_stations")
            self._dispatch_probs = arr.reshape(1, -1)
        elif arr.ndim == 2:
            if arr.shape[1] != self._num_stations:
                raise ValueError("probs must have num_stations columns")
            self._dispatch_probs = arr.copy()
        else:
            raise ValueError("probs must be 1-D or 2-D")
        self._dispatching = RoutingStrategy.PROB
        return self

    def set_weights(self, weights: Sequence[int]) -> 'Cluster':
        """WRROBIN dispatching with per-server integer weights."""
        if len(weights) != self._num_stations:
            raise ValueError("weights length must equal num_stations")
        clean = [int(w) for w in weights]
        if any(abs(w - round(w)) > 1e-9 for w in weights):
            raise ValueError("WRROBIN weights must be integers")
        self._dispatch_weights = clean
        self._dispatching = RoutingStrategy.WRROBIN
        return self

    def set_arrival_scv(self, scv) -> 'Cluster':
        """SCV of the arrival process (open farms only).

        Either a scalar (broadcast to all classes) or a sequence of length R.
        SCV != 1 swaps the per-class Exp distribution for an APH fitted to the
        same mean and the supplied SCV.
        """
        R = len(self._arrival_rates)
        arr = np.atleast_1d(np.asarray(scv, dtype=float))
        if arr.size == 1:
            arr = np.full(R, float(arr))
        elif arr.shape != (R,):
            raise ValueError("arrival_scv length must equal number of classes")
        if np.any(arr <= 0):
            raise ValueError("SCV must be positive")
        self._arrival_scvs = arr
        return self

    def set_service_scv(self, scv) -> 'Cluster':
        """SCV of the service distribution.

        Either a scalar (broadcast), a length-M vector (per-server, broadcast
        over classes), or a (M, R) matrix.
        """
        R = len(self._population) if self._closed else len(self._arrival_rates)
        M = self._num_stations
        arr = np.atleast_2d(np.asarray(scv, dtype=float))
        if arr.size == 1:
            arr = np.full((M, R), float(arr))
        elif arr.ndim == 1 or arr.shape == (1, M):
            arr = np.tile(np.asarray(scv, dtype=float).reshape(M, 1), (1, R))
        elif arr.shape != (M, R):
            raise ValueError("service_scv must be scalar, length-M vector, or (M, R) matrix")
        if np.any(arr <= 0):
            raise ValueError("SCV must be positive")
        self._service_scvs = arr
        return self

    def set_sq(self, d: int) -> 'Cluster':
        """SQ(d) dispatching: shortest of d sampled servers."""
        if d < 1 or int(d) != d:
            raise ValueError("d must be a positive integer")
        self._sq_k = int(d)
        self._dispatching = RoutingStrategy.SQ
        return self

    def set_closed(self, population: Union[int, Sequence[int]],
                   think_time: Union[float, Sequence[float]]) -> 'Cluster':
        if np.isscalar(population):
            self._population = [int(population)]
            self._think_times = [float(think_time)]
        else:
            self._population = [int(p) for p in population]
            self._think_times = [float(z) for z in think_time]
            if len(self._population) != len(self._think_times):
                raise ValueError("population and think_time must have equal length")
        self._closed = True
        return self

    # camelCase aliases (mirror Java)
    setDispatching = set_dispatching
    setScheduling = set_scheduling
    setStationServers = set_station_servers
    setClosed = set_closed
    setProbabilities = set_probabilities
    setWeights = set_weights
    setSQ = set_sq
    setArrivalSCV = set_arrival_scv
    setServiceSCV = set_service_scv
    setNumStations = set_num_stations
    setArrivalRate = set_arrival_rate
    setArrivalRates = set_arrival_rates
    setServiceRate = set_service_rate
    setServiceRates = set_service_rates

    # ------------------------------------------------------------------ build

    def build(self) -> Network:
        """Build the configured cluster Network."""
        M = self._num_stations
        if self._closed:
            R = len(self._population)
        else:
            R = len(self._arrival_rates)

        strategies = [self._scheduling] * M
        S = np.array(self._station_counts, dtype=int)

        # Convert per-server, per-class service *rates* to *mean service times*.
        if self._service_rates.shape[1] == 1 and R > 1:
            # Single rate replicated across all classes.
            rates = np.tile(self._service_rates, (1, R))
        else:
            rates = self._service_rates
        if np.any(rates <= 0):
            raise ValueError("Service rates must be positive")
        D = 1.0 / rates

        # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
        needs_post = (
            (self._dispatching == RoutingStrategy.PROB and self._dispatch_probs is not None)
            or (self._dispatching == RoutingStrategy.WRROBIN and self._dispatch_weights is not None)
            or (self._dispatching == RoutingStrategy.SQ and self._sq_k is not None)
        )
        factory_dispatch = RoutingStrategy.RAND if needs_post else self._dispatching

        if self._closed:
            N = np.array(self._population, dtype=int).reshape(1, -1)
            Z = np.array(self._think_times, dtype=float).reshape(1, -1)
            model = Network.cluster_closed(N, Z, D, strategies, S, factory_dispatch)
        else:
            lam = np.array(self._arrival_rates, dtype=float).reshape(1, -1)
            model = Network.cluster(lam, D, strategies, S, factory_dispatch)

        self._apply_dispatcher_config(model, R)
        self._apply_distribution_scvs(model, R)
        return model

    # --------------------------------------------------------------- SCV apply
    def _apply_distribution_scvs(self, model: Network, R: int) -> None:
        # Lazy import to avoid pulling in distributions during module import.
        from ..distributions.markovian import APH

        classes = model.classes
        if not self._closed and self._arrival_scvs is not None:
            src = model.getNodeByName('Source')
            for r in range(R):
                scv = float(self._arrival_scvs[r])
                if scv != 1.0:
                    src.setArrival(classes[r],
                                   APH.fit_mean_and_scv(1.0 / self._arrival_rates[r], scv))

        if self._service_scvs is not None:
            rates = self._service_rates
            if rates.shape[1] == 1 and R > 1:
                rates = np.tile(rates, (1, R))
            for i in range(self._num_stations):
                server = model.getNodeByName(f'Station{i + 1}')
                for r in range(R):
                    scv = float(self._service_scvs[i, r])
                    if scv != 1.0:
                        server.setService(classes[r],
                                          APH.fit_mean_and_scv(1.0 / rates[i, r], scv))

    # ------------------------------------------------------------------ apply
    def _apply_dispatcher_config(self, model: Network, R: int) -> None:
        dispatcher = model.getNodeByName('Dispatcher')
        classes = model.classes
        servers = [model.getNodeByName(f'Station{i + 1}') for i in range(self._num_stations)]

        if self._dispatching == RoutingStrategy.PROB and self._dispatch_probs is not None:
            # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
            rm = model.init_routing_matrix()
            source_or_think = model.getNodeByName('Source') if not self._closed \
                else model.getNodeByName('Think')
            sink_or_think = model.getNodeByName('Sink') if not self._closed \
                else model.getNodeByName('Think')
            for r in range(R):
                cls = classes[r]
                rm.set(cls, cls, source_or_think, dispatcher, 1.0)
                row = self._dispatch_probs[0] if self._dispatch_probs.shape[0] == 1 \
                    else self._dispatch_probs[r]
                for i, server in enumerate(servers):
                    rm.set(cls, cls, dispatcher, server, float(row[i]))
                    rm.set(cls, cls, server, sink_or_think, 1.0)
            model.link(rm)
            for cls in classes:
                dispatcher.set_routing(cls, RoutingStrategy.PROB)
        elif self._dispatching == RoutingStrategy.WRROBIN and self._dispatch_weights is not None:
            # The Router doesn't subclass Station so it lacks _routing_weights;
            # attach the per-class weight map directly so _refresh_chains picks it up.
            if not hasattr(dispatcher, '_routing_weights') or dispatcher._routing_weights is None:
                dispatcher._routing_weights = {}
            for r in range(R):
                cls = classes[r]
                dispatcher._routing_weights[cls] = {
                    servers[i]: int(self._dispatch_weights[i])
                    for i in range(self._num_stations)
                }
                dispatcher.set_routing(cls, RoutingStrategy.WRROBIN)
            if hasattr(dispatcher, '_invalidate_java'):
                dispatcher._invalidate_java()
            if hasattr(model, '_has_struct'):
                model._has_struct = False
                model._sn = None
        elif self._dispatching == RoutingStrategy.SQ and self._sq_k is not None:
            for r in range(R):
                dispatcher.set_routing(classes[r], RoutingStrategy.SQ,
                                       self._sq_k)
            if hasattr(model, '_has_struct'):
                model._has_struct = False
                model._sn = None

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _solver_factory(solver_class) -> Callable[[Network], 'object']:
        """Build a `model -> avg_table` callable from a solver class."""
        if callable(solver_class) and not isinstance(solver_class, type):
            # User passed in a custom factory.
            return solver_class

        def _factory(model: Network):
            solver = solver_class(model)
            # Most solvers expose .get_avg_table; fall back to .getAvgTable.
            if hasattr(solver, 'get_avg_table'):
                return solver.get_avg_table()
            return solver.getAvgTable()
        return _factory

    # ------------------------------------------------------------------ compare

    def compare_dispatching(self, solver_class,
                            policies: Sequence[RoutingStrategy]) -> 'OrderedDict':
        """Solve the cluster under each dispatching policy and return an ordered
        dict of policy -> avg-table."""
        factory = self._solver_factory(solver_class)
        out = OrderedDict()
        saved = self._dispatching
        try:
            for p in policies:
                self._dispatching = p
                out[p] = factory(self.build())
        finally:
            self._dispatching = saved
        return out

    def compare_scheduling(self, solver_class,
                           disciplines: Sequence[SchedStrategy]) -> 'OrderedDict':
        factory = self._solver_factory(solver_class)
        out = OrderedDict()
        saved = self._scheduling
        try:
            for s in disciplines:
                self._scheduling = s
                out[s] = factory(self.build())
        finally:
            self._scheduling = saved
        return out

    # ------------------------------------------------------------------ sweeps

    def sweep_arrival_rate(self, rates: Sequence[float], solver_class) -> 'OrderedDict':
        if self._closed:
            raise RuntimeError("sweep_arrival_rate is only defined for open farms")
        if len(self._arrival_rates) != 1:
            raise RuntimeError("sweep_arrival_rate requires a single class")
        factory = self._solver_factory(solver_class)
        out = OrderedDict()
        saved = self._arrival_rates[0]
        try:
            for r in rates:
                self._arrival_rates[0] = float(r)
                out[float(r)] = factory(self.build())
        finally:
            self._arrival_rates[0] = saved
        return out

    def sweep_num_stations(self, counts: Sequence[int], solver_class) -> 'OrderedDict':
        factory = self._solver_factory(solver_class)
        out = OrderedDict()
        saved_M = self._num_stations
        saved_rates = self._service_rates
        saved_counts = self._station_counts
        try:
            per_class_rate = saved_rates[0].copy()
            for m in counts:
                m = int(m)
                self._num_stations = m
                self._service_rates = np.tile(per_class_rate.reshape(1, -1), (m, 1))
                self._station_counts = [1] * m
                out[m] = factory(self.build())
        finally:
            self._num_stations = saved_M
            self._service_rates = saved_rates
            self._station_counts = saved_counts
        return out

    # camelCase aliases
    compareDispatching = compare_dispatching
    compareScheduling = compare_scheduling
    sweepArrivalRate = sweep_arrival_rate
    sweepNumStations = sweep_num_stations


__all__ = ['Cluster']
