"""
Node implementations for LINE native Python.

This module provides concrete node types: Queue, Delay, Source, Sink, Fork, Join.
Ported from MATLAB implementation in matlab/src/lang/nodes/
"""

from typing import Optional, Dict, Union
from enum import IntEnum
import numpy as np

from .base import (
    Node, Station, StatefulNode, JobClass, NodeType, SchedStrategy,
    SchedStrategyType, JoinStrategy, DropStrategy, ReplacementStrategy,
    RoutingStrategy
)
from ..distributions.continuous import Immediate


class Queue(Station):
    """
    Queue station for customer service.

    A Queue represents a queueing station with one or more servers
    and configurable scheduling strategy.
    """

    def __init__(self, model, name: str, sched_strategy: SchedStrategy = SchedStrategy.FCFS):
        """
        Initialize a Queue node.

        Args:
            model: Network instance
            name: Queue name
            sched_strategy: Scheduling strategy (default: FCFS)

        Raises:
            ValueError: If scheduling strategy is invalid
        """
        super().__init__(NodeType.QUEUE, name)
        self.set_model(model)
        model.add_node(self)

        self._sched_strategy = sched_strategy
        # INF SCHEDULING IS AN INFINITE SERVER, and the server count says so:
        # MATLAB's Queue constructor sets numberOfServers = Inf on this branch
        # (and setNumberOfServers then ignores a later change), the JAR sets its
        # Integer.MAX_VALUE, and every partition of the stations into delays and
        # queues reads that count. Left at 1, the station was a Queue node with
        # INF scheduling that no analyzer claimed: the AMVA excluded it from the
        # queueing demands by its scheduling and the delay fill skipped it by its
        # node type, so it reported QLen 0 while holding most of the population.
        # the name, not the enum member: the user-facing SchedStrategy is
        # line_solver.constants', a distinct class from lang.base's, so the two
        # INF members compare unequal (as the DPS/GPS guards below already assume)
        if getattr(sched_strategy, 'name', None) == 'INF':
            self._number_of_servers = np.inf
        self._sched_param = {}  # Dict[JobClass, scheduling_parameter]
        self._service_process = {}  # Dict[JobClass, Distribution]
        self._server_types = None  # For heterogeneous servers
        # Pass-and-swap (PAS) properties
        self._swap_graph = None    # (nclasses x nclasses) class compatibility/swap graph
        self._svc_rate_fun = None  # callable mu(c): total service rate of ordered state c

        # Initialize scheduling policy type
        self._set_sched_policy()

    def _set_sched_policy(self) -> None:
        """Determine if scheduling strategy is preemptive or non-preemptive."""
        preemptive_strategies = {
            SchedStrategy.LCFSPR,
            SchedStrategy.FCFSPR,
            SchedStrategy.EDF
        }
        self._sched_policy = (
            SchedStrategyType.PR if self._sched_strategy in preemptive_strategies
            else SchedStrategyType.NP
        )

    def get_sched_strategy(self) -> SchedStrategy:
        """Get scheduling strategy."""
        return self._sched_strategy

    def set_sched_strategy(self, strategy: SchedStrategy) -> None:
        """
        Set scheduling strategy.

        Args:
            strategy: Scheduling strategy (FCFS, LCFS, PS, etc.)
        """
        self._sched_strategy = strategy
        self._set_sched_policy()
        self._invalidate_java()

    def get_sched_policy(self) -> SchedStrategyType:
        """Get scheduling policy (preemptive or non-preemptive)."""
        return self._sched_policy

    def set_service(self, jobclass: JobClass, distribution=None, weight: float = None) -> None:
        """
        Set service time distribution for a job class.

        On a pass-and-swap (PAS) queue, call as ``set_service(mu)`` where ``mu``
        is a callable taking the ordered state vector ``c`` (a list/array of
        class indices, ``c[0]`` the oldest job) and returning the scalar total
        service rate ``mu(c)``. Per-class service distributions are not used by
        an order-independent/PAS queue.

        Args:
            jobclass: Job class, or the mu(c) callable on a PAS queue
            distribution: Service distribution (Exp, Erlang, HyperExp, etc.)
                         or Workflow for activity-based distributions
            weight: Optional weight/scale factor for the distribution
        """
        if getattr(self._sched_strategy, 'name', None) in ('PAS', 'OI'):
            self.set_service_rate_function(jobclass)
            return
        previous = self._service_process.get(jobclass)
        # Store the distribution (and optional weight for future use)
        self._service_process[jobclass] = distribution
        if weight is not None:
            # Store weight metadata for analysis
            if not hasattr(self, '_service_weights'):
                self._service_weights = {}
            self._service_weights[jobclass] = weight
            # Also set scheduling parameter for DPS/GPS scheduling
            self._sched_param[jobclass] = weight
        self._invalidate_java()
        if self._model is not None:
            # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
            if self._structurally_compatible_service(previous, distribution):
                self._model._rates_dirty = True
            else:
                self._model._sn = None
                self._model._has_struct = False

    @staticmethod
    def _structurally_compatible_service(previous, distribution) -> bool:
        """True when swapping ``previous`` for ``distribution`` cannot change any
        structural field of the compiled struct, i.e. both are set and have the
        same number of phases (a disabled process reports 0 phases and so is
        never compatible with an enabled one)."""
        if previous is None or distribution is None:
            return False
        from .network import Network
        try:
            return Network._dist_num_phases(previous) == Network._dist_num_phases(distribution)
        except Exception:
            return False

    def get_service_process(self, jobclass: JobClass):
        """Get the service process (distribution) for a job class (wrapper-compatible name)."""
        return self._service_process.get(jobclass)

    getServiceProcess = get_service_process
    service_process = get_service_process

    def get_service(self, jobclass: JobClass):
        """
        Get service distribution for a job class.

        Args:
            jobclass: Job class

        Returns:
            Service distribution or None if not set
        """
        return self._service_process.get(jobclass)

    def set_item_service_rate(self, cache, jobin_class: JobClass, item: int,
                              service_rate: float) -> None:
        """Override, at this queue, the retrieval service rate for a single item of
        the read class ``jobin_class`` in ``cache``'s retrieval system. The default
        (when not overridden) is the read class's own service distribution at this
        queue. ``item`` is 0-based."""
        from ..distributions.continuous import Exp
        r_class = cache.get_retrieval_class(jobin_class, item)
        if r_class is None:
            raise ValueError("No retrieval class defined for the given class/item; "
                             "call set_retrieval_system first.")
        self.set_service(r_class, Exp(service_rate))

    setItemServiceRate = set_item_service_rate

    def set_strategy_param(self, jobclass: JobClass, param) -> None:
        """
        Set scheduling parameter for a job class.

        For PS/DPS/GPS: weight or priority parameter
        For other strategies: may be unused

        Args:
            jobclass: Job class
            param: Scheduling parameter (weight, priority, etc.)
        """
        self._sched_param[jobclass] = param
        self._invalidate_java()

    def get_strategy_param(self, jobclass: JobClass):
        """Get scheduling parameter for a job class."""
        return self._sched_param.get(jobclass, 1.0)

    def set_swap_graph(self, graph) -> None:
        """
        Set the class compatibility/swap graph of a pass-and-swap (PAS) queue.

        Args:
            graph: (nclasses x nclasses) adjacency matrix; entry (r,s) nonzero iff,
                   upon completion of a class-r job, a waiting class-s job may take
                   its place (order-independent swap). Undirected; self-loops allowed.
        """
        if getattr(self._sched_strategy, 'name', None) == 'OI':
            raise ValueError("set_swap_graph is not applicable to OI (order-independent) queues; "
                             "their swap graph is always zero. Use SchedStrategy.PAS to configure "
                             "a non-trivial swap graph.")
        if getattr(self._sched_strategy, 'name', None) != 'PAS':
            raise ValueError("set_swap_graph is only applicable to PAS (pass-and-swap) queues.")
        import numpy as _np
        g = _np.asarray(graph, dtype=float)
        K = len(self._model.get_classes()) if self._model is not None else g.shape[0]
        if g.shape[0] != K or g.shape[1] != K:
            raise ValueError(f"Swap graph must be a {K}x{K} matrix (nclasses x nclasses).")
        self._swap_graph = g
        if self._model is not None:
            self._model._sn = None

    def get_swap_graph(self):
        """Return the (nclasses x nclasses) swap graph of a PAS queue, or None."""
        return self._swap_graph

    def set_service_rate_function(self, mu_fun) -> None:
        """
        Set the total service rate function mu(c) of a pass-and-swap (PAS) queue.

        Args:
            mu_fun: callable taking the ordered state vector c (list/array of class
                    indices, c[0] oldest) and returning the scalar total service
                    rate mu(c). The rate of position i is the increment
                    Delta_mu(c[:i]) = mu(c[:i]) - mu(c[:i-1]).
        """
        if getattr(self._sched_strategy, 'name', None) not in ('PAS', 'OI'):
            raise ValueError("set_service_rate_function is only applicable to PAS (pass-and-swap) "
                             "and OI (order-independent) queues.")
        if not callable(mu_fun):
            raise ValueError("PAS queues require a service rate function mu(c); "
                             "per-class service distributions are not supported. "
                             "Use set_service(lambda c: ...).")
        self._svc_rate_fun = mu_fun
        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        import numpy as _np
        from ..distributions.continuous import Exp as _Exp, Disabled as _Disabled
        if self._model is not None:
            for r, jc in enumerate(self._model.get_classes()):
                rate_r = float(mu_fun(_np.array([r])))
                if _np.isfinite(rate_r) and rate_r > 0:
                    self._service_process[jc] = _Exp(rate_r)
                else:
                    self._service_process[jc] = _Disabled()
            self._model._sn = None

    def get_service_rate_function(self):
        """Return the mu(c) service rate function of a PAS queue, or None."""
        return self._svc_rate_fun

    def check_perm_invariance(self, Nvec, cap):
        """Check the order-independence (OI) condition on the service rate mu(c):
        the rate of the job in position j must depend only on the jobs at or
        ahead of it (positions 1..j) and not on those behind. Since the
        position-j rate is the prefix increment mu(c[:j]) - mu(c[:j-1]),
        tail-independence is structural; the substantive requirement is that
        this increment be independent of the order of the jobs ahead, which (by
        induction on prefix length) is equivalent to mu(c) being permutation-
        invariant. Enumerates the reachable multisets (per-class counts <= Nvec,
        total <= cap) when small, otherwise samples and returns partial=True.
        Returns (ok, badc, partial)."""
        import numpy as _np
        from math import lgamma, exp
        mu = self._svc_rate_fun
        if mu is None:
            return True, None, False
        K = len(Nvec)
        tol = 1e-9
        PERM_ENUM = 5040     # enumerate all distinct permutations up to this
        PERM_SAMPLE = 16     # permutations sampled per multiset above PERM_ENUM
        LATTICE_BUDGET = 4096
        MAXEVAL = 50000
        ub = _np.asarray(Nvec, dtype=float)
        has_open = bool(_np.any(~_np.isfinite(ub)))
        if _np.isfinite(cap) and 0 <= cap < 1e18:
            Lmax = float(cap)
        else:
            Lmax = float(_np.sum(ub[_np.isfinite(ub)]))
        ub = _np.where(_np.isfinite(ub), ub, min(Lmax, 6.0))
        ub = _np.minimum(ub, Lmax).astype(int)
        if not _np.isfinite(Lmax) or Lmax < 2:
            return True, None, False
        Lmax = int(Lmax)
        lat = int(_np.prod(ub + 1))
        exhaustive = (not has_open) and lat <= LATTICE_BUDGET
        state = {'ok': True, 'badc': None, 'partial': False, 'neval': 0}

        def ms_perms(c0):
            # Distinct permutations of the multiset c0 (multiset-aware).
            if len(c0) <= 1:
                return [list(c0)]
            out = []
            for u in sorted(set(c0)):
                rest = list(c0)
                rest.remove(u)
                for sub in ms_perms(rest):
                    out.append([u] + sub)
            return out

        def test_multiset(n):
            c0 = _np.repeat(_np.arange(K), n)
            length = len(c0)
            logd = lgamma(length + 1) - sum(lgamma(int(x) + 1) for x in n if x > 0)
            dcount = int(round(exp(logd)))
            base = float(mu(c0)); state['neval'] += 1
            if dcount <= PERM_ENUM:
                perms = ms_perms([int(x) for x in c0])
            else:
                state['partial'] = True
                perms = [list(c0[::-1])]
                for _ in range(PERM_SAMPLE - 1):
                    perms.append(list(_np.random.permutation(c0)))
            for p in perms:
                v = float(mu(_np.asarray(p, dtype=int)))
                state['neval'] += 1
                if abs(v - base) > tol * max(1.0, abs(base)):
                    state['ok'] = False
                    state['badc'] = [int(x) for x in c0]
                    return

        rng_state = _np.random.get_state()
        _np.random.seed(0)
        try:
            if exhaustive:
                n = _np.zeros(K, dtype=int)
                while True:
                    if n.sum() >= 2 and _np.count_nonzero(n) >= 2 and n.sum() <= Lmax:
                        test_multiset(n)
                        if not state['ok']:
                            break
                        if state['neval'] >= MAXEVAL:
                            state['partial'] = True
                            break
                    d = 0
                    while d < K:
                        n[d] += 1
                        if n[d] <= ub[d]:
                            break
                        n[d] = 0
                        d += 1
                    if d >= K:
                        break
            else:
                state['partial'] = True
                for _ in range(400):
                    length = 2 + int(_np.random.randint(0, max(1, min(Lmax, 6) - 1)))
                    n = _np.zeros(K, dtype=int)
                    for _j in range(length):
                        r = int(_np.random.randint(0, K))
                        if n[r] < ub[r]:
                            n[r] += 1
                    if n.sum() >= 2 and _np.count_nonzero(n) >= 2:
                        test_multiset(n)
                        if not state['ok'] or state['neval'] >= MAXEVAL:
                            break
        finally:
            _np.random.set_state(rng_state)
        return state['ok'], state['badc'], state['partial']

    def check_rate_monotonicity(self, Nvec, cap):
        """Check OI condition (1) on the service rate mu(c): the per-job rates
        must be non-negative, mu(c_1,...,c_j) >= mu(c_1,...,c_{j-1}) for every
        microstate and position.

        A rate can be permutation-invariant and still fail to parameterize an OI
        queue. Single-server processor sharing with class-dependent rates,
        mu(c) = (sum_j mu_{c_j}) / n, is the standard trap: it is flatly
        invariant under permutations, yet as soon as two classes have different
        rates its prefix increments go negative -- with mu_hit = 3.0 and
        mu_miss = 0.7, mu(Hit) = 3.0 and mu(Hit,Miss) = 1.85, so the second job
        would be served at rate -1.15.

        Run this AFTER check_perm_invariance: permutation invariance is what
        makes mu a function of the count vector, and the increments to test are
        then just mu(n + e_r) - mu(n) over count vectors n and classes r, with
        no permutation enumeration. Prefixes are non-empty, so mu is never
        evaluated on an empty microstate, and an increment of exactly zero is
        accepted -- that is how a class which does not visit this station is
        expressed. Enumerates the reachable count vectors (per-class counts <=
        Nvec, total <= cap) when small, otherwise samples and returns
        partial=True. Returns (ok, badc, badr, partial): badc the microstate
        whose rate is lowered, badr the class whose arrival lowers it.
        """
        import numpy as _np
        mu = self._svc_rate_fun
        if mu is None:
            return True, None, None, False
        K = len(Nvec)
        tol = 1e-9
        LATTICE_BUDGET = 4096
        MAXEVAL = 50000
        ub = _np.asarray(Nvec, dtype=float)
        has_open = bool(_np.any(~_np.isfinite(ub)))
        if _np.isfinite(cap) and 0 <= cap < 1e18:
            Lmax = float(cap)
        else:
            Lmax = float(_np.sum(ub[_np.isfinite(ub)]))
        ub = _np.where(_np.isfinite(ub), ub, min(Lmax, 6.0))
        ub = _np.minimum(ub, Lmax).astype(int)
        if not _np.isfinite(Lmax) or Lmax < 2:
            return True, None, None, False
        Lmax = int(Lmax)
        lat = int(_np.prod(ub + 1))
        exhaustive = (not has_open) and lat <= LATTICE_BUDGET

        cache = {}
        neval = [0]

        def rate(n):
            key = tuple(int(x) for x in n)
            if key not in cache:
                cache[key] = float(mu(_np.repeat(_np.arange(K), n)))
                neval[0] += 1
            return cache[key]

        def test_prefix(n):
            """Every one-job extension of the prefix n must not lower the rate."""
            base = rate(n)
            for r in range(K):
                if n[r] >= ub[r]:
                    continue
                n[r] += 1
                try:
                    if rate(n) - base < -tol * max(1.0, abs(base)):
                        return [int(x) for x in _np.repeat(_np.arange(K), n)], r
                finally:
                    n[r] -= 1
            return None

        partial = False
        rng_state = _np.random.get_state()
        _np.random.seed(0)
        try:
            if exhaustive:
                n = _np.zeros(K, dtype=int)
                while True:
                    # non-empty prefixes with room for one more job
                    if 1 <= n.sum() <= Lmax - 1:
                        bad = test_prefix(n)
                        if bad is not None:
                            return False, bad[0], bad[1], partial
                        if neval[0] >= MAXEVAL:
                            partial = True
                            break
                    d = 0
                    while d < K:
                        n[d] += 1
                        if n[d] <= ub[d]:
                            break
                        n[d] = 0
                        d += 1
                    if d >= K:
                        break
            else:
                partial = True
                for _ in range(400):
                    length = 1 + int(_np.random.randint(0, max(1, min(Lmax, 6))))
                    n = _np.zeros(K, dtype=int)
                    for _j in range(length):
                        r = int(_np.random.randint(0, K))
                        if n[r] < ub[r]:
                            n[r] += 1
                    if 1 <= n.sum() <= Lmax - 1:
                        bad = test_prefix(n)
                        if bad is not None:
                            return False, bad[0], bad[1], partial
                        if neval[0] >= MAXEVAL:
                            break
        finally:
            _np.random.set_state(rng_state)
        return True, None, None, partial

    def is_service_defined(self, jobclass: JobClass) -> bool:
        """Check if service distribution is defined for a class."""
        return jobclass in self._service_process

    def is_service_disabled(self, jobclass: JobClass) -> bool:
        """Check if service is explicitly disabled for a class."""
        return self._service_process.get(jobclass) is None

    def get_service_rates(self) -> Dict[JobClass, float]:
        """
        Get mean service rates for all classes.

        Returns:
            Dict mapping classes to mean service rates
        """
        rates = {}
        for jobclass, dist in self._service_process.items():
            if dist is not None and hasattr(dist, 'getMean'):
                rates[jobclass] = dist.getMean()
        return rates

    def set_number_of_servers(self, value: int) -> None:
        """
        Set the number of servers. Matches MATLAB Queue.setNumServers: DPS/GPS
        scheduling does not admit multi-server stations.
        """
        if getattr(self._sched_strategy, 'name', None) in ('DPS', 'GPS') and \
                not np.isinf(value) and value != 1:
            raise ValueError(
                f"Cannot use multi-server stations with "
                f"{getattr(self._sched_strategy, 'name', self._sched_strategy)} scheduling.")
        # An infinite server has no server count to set: MATLAB's
        # setNumberOfServers ignores the request on this branch rather than
        # turning the station into a finite-server queue behind its scheduling.
        if getattr(self._sched_strategy, 'name', None) == 'INF':
            return
        super().set_number_of_servers(value)

    def set_limit(self, limit: int) -> None:
        """
        Set the maximum number of jobs in service for LPS (Limited Processor
        Sharing) scheduling. Mirrors MATLAB Queue.setLimit: the limit is stored
        as the queue-level scheduling parameter (schedparam column 1).
        """
        if getattr(self._sched_strategy, 'name', None) != 'LPS':
            from ..api.io.logging import line_warning
            line_warning('Queue.set_limit', 'setLimit is only applicable to LPS scheduling.')
            return
        self._lps_limit = limit
        self._invalidate_java()

    def get_limit(self):
        """
        Return the LPS concurrency limit set via set_limit, or numberOfServers
        when no explicit limit was set (matching MATLAB Queue.getLimit).
        """
        lim = getattr(self, '_lps_limit', None)
        if lim is None:
            return self.number_of_servers
        return lim

    # PascalCase aliases for MATLAB compatibility
    getSchedStrategy = get_sched_strategy
    setSchedStrategy = set_sched_strategy
    getSchedPolicy = get_sched_policy
    setService = set_service
    getService = get_service
    setStrategyParam = set_strategy_param
    getStrategyParam = get_strategy_param
    getServiceRates = get_service_rates
    setSwapGraph = set_swap_graph
    getSwapGraph = get_swap_graph
    setServiceRateFunction = set_service_rate_function
    getServiceRateFunction = get_service_rate_function
    # Re-alias so the DPS/GPS multi-server guard is not bypassed via the base aliases
    setNumberOfServers = set_number_of_servers
    set_num_servers = set_number_of_servers
    setLimit = set_limit
    getLimit = get_limit

    # =========================================================================
    # HETEROGENEOUS SERVER METHODS
    # =========================================================================

    def is_heterogeneous(self) -> bool:
        """
        Check if this queue has heterogeneous servers.

        Returns:
            True if server types have been defined
        """
        return self._server_types is not None and len(self._server_types) > 0

    def add_server_type(self, server_type) -> None:
        """
        Add a server type to this queue.

        Args:
            server_type: ServerType object to add

        Raises:
            ValueError: If server type is already in queue
        """
        if self._server_types is None:
            self._server_types = []

        if server_type in self._server_types:
            raise ValueError(f"Server type '{server_type.name}' already in queue")

        server_type.set_id(len(self._server_types))
        server_type.set_parent_queue(self)
        self._server_types.append(server_type)
        # The pools size the station, as in MATLAB/JAR updateTotalServerCount
        self.number_of_servers = sum(st.get_num_of_servers() for st in self._server_types)
        self._invalidate_java()

    def get_server_types(self):
        """Get the list of server types."""
        return self._server_types if self._server_types else []

    def set_hetero_sched_policy(self, policy) -> None:
        """
        Set the heterogeneous scheduling policy.

        Args:
            policy: HeteroSchedPolicy value
        """
        self._hetero_sched_policy = policy
        self._invalidate_java()

    def get_hetero_sched_policy(self):
        """Get the heterogeneous scheduling policy."""
        return getattr(self, '_hetero_sched_policy', None)

    def set_hetero_service(self, jobclass, server_type, distribution) -> None:
        """
        Set service distribution for a specific job class and server type.

        Args:
            jobclass: Job class
            server_type: ServerType object
            distribution: Service distribution
        """
        if not hasattr(self, '_hetero_service'):
            self._hetero_service = {}

        key = (jobclass, server_type)
        self._hetero_service[key] = distribution
        self._invalidate_java()

    def get_hetero_service(self, jobclass, server_type):
        """
        Get service distribution for a specific job class and server type.

        Args:
            jobclass: Job class
            server_type: ServerType object

        Returns:
            Service distribution or None
        """
        if not hasattr(self, '_hetero_service'):
            return None
        return self._hetero_service.get((jobclass, server_type))

    def set_server_parallelism(self, jobclass, n: int) -> None:
        """
        Set the number of servers a job of this class seizes for the whole of its
        service, JMT's job parallelism (Server.serverNumRequired).

        A job waits until n servers are simultaneously free and holds all of them
        until it completes, so the station serves at most floor(c/n) such jobs at
        a time. The default is 1.

        Args:
            jobclass: Job class
            n: Number of servers required, an integer in [1, c]

        Raises:
            ValueError: If n is below 1 or above the server count
        """
        n = int(n)
        if n < 1:
            raise ValueError("Server parallelism must be a positive integer.")
        if n > self.number_of_servers:
            raise ValueError(
                f"Server parallelism ({n}) exceeds the {self.number_of_servers} servers of "
                f"station {self.name}, so a job of class {jobclass.name} could never enter service.")
        if not hasattr(self, '_server_parallelism'):
            self._server_parallelism = {}
        self._server_parallelism[jobclass] = n
        self._invalidate_java()

    def get_server_parallelism(self, jobclass) -> int:
        """Get the number of servers seized by a job of this class, 1 if unset."""
        return getattr(self, '_server_parallelism', {}).get(jobclass, 1)

    def has_server_parallelism(self) -> bool:
        """True when some class seizes more than one server."""
        return any(n > 1 for n in getattr(self, '_server_parallelism', {}).values())

    def get_total_num_of_servers(self) -> int:
        """
        Get total number of servers across all server types.

        Returns:
            Total server count
        """
        if not self._server_types:
            return self.number_of_servers
        return sum(st.get_num_of_servers() for st in self._server_types)

    # PascalCase aliases for heterogeneous methods
    isHeterogeneous = is_heterogeneous
    addServerType = add_server_type
    getServerTypes = get_server_types
    setHeteroSchedPolicy = set_hetero_sched_policy
    getHeteroSchedPolicy = get_hetero_sched_policy
    setHeteroService = set_hetero_service
    getHeteroService = get_hetero_service
    getTotalNumOfServers = get_total_num_of_servers
    setServerParallelism = set_server_parallelism
    getServerParallelism = get_server_parallelism
    hasServerParallelism = has_server_parallelism

    # =========================================================================
    # POLLING AND SWITCHOVER METHODS
    # =========================================================================

    def set_polling_type(self, polling_type, k=None) -> None:
        """
        Set the polling type for this queue.

        Args:
            polling_type: PollingType enum value (EXHAUSTIVE, GATED, KLIMITED)
            k: For KLIMITED polling, the maximum number of jobs to serve (default: 1)
        """
        self._polling_type = polling_type
        self._polling_k = k if k is not None else 1
        self._invalidate_java()

        # Set default Immediate switchover for all classes (matching MATLAB behavior)
        # This is required for JMT to properly simulate polling systems
        model = self.get_model()
        if model is not None:
            classes = model.get_classes()
            for jobclass in classes:
                self.set_switchover(jobclass, Immediate())

    def get_polling_type(self):
        """Get the polling type for this queue."""
        return getattr(self, '_polling_type', None)

    def set_polling_k(self, k):
        """Set the k parameter for k-limited polling (wrapper-compatible name)."""
        self._polling_k = int(k)

    setPollingK = set_polling_k

    def set_switchover(self, from_class, to_class_or_dist, distribution=None) -> None:
        """
        Set switchover time between job classes.

        Can be called as:
        - set_switchover(from_class, to_class, distribution) - for class-to-class switchover
        - set_switchover(jobclass, distribution) - for single class switchover time

        Args:
            from_class: Source job class
            to_class_or_dist: Target job class, OR distribution if only 2 args
            distribution: Switchover time distribution (optional if 2 args)
        """
        if not hasattr(self, '_switchover'):
            self._switchover = {}

        if distribution is None:
            # Called as set_switchover(class, distribution)
            self._switchover[from_class] = to_class_or_dist
        else:
            # Called as set_switchover(from_class, to_class, distribution)
            self._switchover[(from_class, to_class_or_dist)] = distribution
        self._invalidate_java()

    def get_switchover(self, from_class, to_class):
        """Get switchover time distribution between classes."""
        if not hasattr(self, '_switchover'):
            return None
        return self._switchover.get((from_class, to_class))

    def set_delay_off(self, jobclass, setup_time, delay_off_time) -> None:
        """
        Set setup and delay-off time distributions for a job class.

        Models a vacation queue: when a server finishes serving class r and no
        more class r jobs are waiting, the server enters a delay-off period.
        When a new class r job arrives during delay-off, the server goes through
        setup before serving.

        Args:
            jobclass: Job class
            setup_time: Distribution for setup time
            delay_off_time: Distribution for delay-off time
        """
        if not hasattr(self, '_setup_time'):
            self._setup_time = {}
        if not hasattr(self, '_delay_off_time'):
            self._delay_off_time = {}
        self._setup_time[jobclass] = setup_time
        self._delay_off_time[jobclass] = delay_off_time
        self._invalidate_java()

    # MATLAB-compatible alias
    setDelayOff = set_delay_off

    def get_setup_time(self, jobclass):
        """Get setup time distribution for a job class."""
        if not hasattr(self, '_setup_time'):
            return None
        return self._setup_time.get(jobclass)

    # MATLAB-compatible alias
    getSetupTime = get_setup_time

    def get_delay_off_time(self, jobclass):
        """Get delay-off time distribution for a job class."""
        if not hasattr(self, '_delay_off_time'):
            return None
        return self._delay_off_time.get(jobclass)

    # MATLAB-compatible alias
    getDelayOffTime = get_delay_off_time

    def is_delay_off_enabled(self) -> bool:
        """Return True if setup and delay-off times are both configured."""
        return bool(getattr(self, '_setup_time', None)) and bool(getattr(self, '_delay_off_time', None))

    # MATLAB-compatible alias
    isDelayOffEnabled = is_delay_off_enabled

    # =========================================================================
    # ROUTING METHODS
    # =========================================================================

    def set_prob_routing(self, jobclass, destination, prob: float) -> None:
        """
        Set probabilistic routing to a destination node.

        Args:
            jobclass: Job class
            destination: Destination node
            prob: Routing probability (0 to 1)
        """
        if not hasattr(self, '_prob_routing'):
            self._prob_routing = {}
        if jobclass not in self._prob_routing:
            self._prob_routing[jobclass] = {}
        self._prob_routing[jobclass][destination] = prob
        self._invalidate_java()

    def get_prob_routing(self, jobclass):
        """Get probabilistic routing for a job class."""
        if not hasattr(self, '_prob_routing'):
            return {}
        return self._prob_routing.get(jobclass, {})

    def set_state(self, state) -> None:
        """
        Set initial state for this node.

        Args:
            state: Array of initial job counts per class
        """
        self._state = np.asarray(state)
        self._state_explicitly_set = True
        self._invalidate_java()
        # Invalidate the model struct so _refresh_state() runs again
        model = self.get_model() if hasattr(self, 'get_model') else None
        if model is not None and hasattr(model, '_has_struct'):
            model._has_struct = False
            model._sn = None

    def get_state(self):
        """Get the current state of this node."""
        return getattr(self, '_state', None)

    def set_routing_weight(self, jobclass: JobClass, destination, weight: float) -> None:
        """
        Set routing weight for weighted round-robin routing.

        Args:
            jobclass: Job class
            destination: Destination node
            weight: Routing weight for this destination
        """
        if not hasattr(self, '_routing_weights'):
            self._routing_weights = {}
        if jobclass not in self._routing_weights:
            self._routing_weights[jobclass] = {}
        self._routing_weights[jobclass][destination] = weight
        self._invalidate_java()

    def get_routing_weight(self, jobclass: JobClass, destination) -> float:
        """Get routing weight for a job class and destination."""
        if not hasattr(self, '_routing_weights'):
            return 1.0
        if jobclass not in self._routing_weights:
            return 1.0
        return self._routing_weights[jobclass].get(destination, 1.0)

    def set_balking(self, jobclass, strategy, thresholds):
        """
        Configure balking for a job class.

        Args:
            jobclass: JobClass object
            strategy: BalkingStrategy enum value
            thresholds: List of (minJobs, maxJobs, probability) tuples
        """
        if not hasattr(self, '_balking_strategies'):
            self._balking_strategies = {}
            self._balking_thresholds = {}
        self._balking_strategies[jobclass] = strategy
        self._balking_thresholds[jobclass] = thresholds
        self._invalidate_java()

    def get_balking(self, jobclass):
        """Get balking config for a job class. Returns (strategy, thresholds) or (None, [])."""
        if not hasattr(self, '_balking_strategies'):
            return None, []
        strategy = self._balking_strategies.get(jobclass)
        thresholds = self._balking_thresholds.get(jobclass, [])
        return strategy, thresholds

    def has_balking(self, jobclass):
        """Check if a job class has balking configured."""
        strategy, _ = self.get_balking(jobclass)
        return strategy is not None

    def set_retrial(self, jobclass, delay_distribution, max_attempts=-1):
        """
        Configure retrial for a job class.

        Args:
            jobclass: JobClass object
            delay_distribution: Distribution for retrial delay
            max_attempts: Maximum retrial attempts (-1 = unlimited)
        """
        if not hasattr(self, '_retrial_delays'):
            self._retrial_delays = {}
            self._retrial_max_attempts = {}
        self._retrial_delays[jobclass] = delay_distribution
        self._retrial_max_attempts[jobclass] = max_attempts
        from .base import DropStrategy
        if max_attempts < 0:
            self._drop_rule = DropStrategy.RETRIAL
        else:
            self._drop_rule = DropStrategy.RETRIAL_WITH_LIMIT
        self._invalidate_java()

    def get_retrial(self, jobclass):
        """Get retrial config for a job class. Returns (delay_dist, max_attempts) or (None, -1)."""
        if not hasattr(self, '_retrial_delays'):
            return None, -1
        delay_dist = self._retrial_delays.get(jobclass)
        max_attempts = self._retrial_max_attempts.get(jobclass, -1)
        return delay_dist, max_attempts

    def has_retrial(self, jobclass):
        """Check if a job class has retrial configured."""
        delay_dist, _ = self.get_retrial(jobclass)
        return delay_dist is not None

    def set_orbit(self, jobclass, retrial_distribution, policy=None, max_orbit=-1):
        """
        Declare this station to be a retrial queue for jobclass: a job that
        finds every server busy joins an orbit and re-attempts entry after a
        random delay, instead of waiting in a line.

        This is the first-class form of the retrial idiom. It removes the
        waiting room itself, so the caller no longer has to know that
        set_capacity(nservers) is the way to express "no waiting room, blocked
        jobs orbit".

        Args:
            jobclass: JobClass object
            retrial_distribution: retrial delay of an orbiting job
            policy: RetrialPolicy.LINEAR (default), where each orbiting job
                retries at its own rate so the aggregate rate is
                (orbit size)*nu, or RetrialPolicy.CONSTANT, where the orbit
                retries as a whole at rate nu whenever non-empty.
            max_orbit: orbit capacity; -1 (default) leaves the orbit unbounded.
                A job that finds the orbit full is lost.

        The mean orbit length is reported by get_avg_orbit / get_avg_orbit_table.
        """
        from .base import RetrialPolicy
        if policy is None:
            policy = RetrialPolicy.LINEAR
        if isinstance(policy, str):
            policy = RetrialPolicy.from_text(policy)
        policy = RetrialPolicy(int(policy))
        if policy not in (RetrialPolicy.LINEAR, RetrialPolicy.CONSTANT):
            raise ValueError('set_orbit requires a RetrialPolicy.LINEAR or '
                             'RetrialPolicy.CONSTANT policy.')
        if max_orbit is None:
            max_orbit = -1
        if not np.isscalar(max_orbit) or (int(max_orbit) != max_orbit) or \
                (max_orbit != -1 and max_orbit < 0):
            raise ValueError('set_orbit requires max_orbit to be -1 (unbounded) '
                             'or a non-negative integer.')
        max_orbit = int(max_orbit)

        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        nservers = self.number_of_servers
        if nservers is None or not np.isfinite(nservers) or nservers < 1:
            nservers = 1
            self.set_number_of_servers(nservers)
        if max_orbit < 0:
            self.set_capacity(np.inf)
        else:
            self.set_capacity(int(nservers) + max_orbit)

        self.set_retrial(jobclass, retrial_distribution, -1)

        if not hasattr(self, '_retrial_policies'):
            self._retrial_policies = {}
            self._orbit_max_jobs = {}
        self._retrial_policies[jobclass] = policy
        self._orbit_max_jobs[jobclass] = max_orbit
        self._invalidate_java()

    def get_orbit(self, jobclass):
        """
        Return (retrial_distribution, policy, max_orbit) for jobclass at this
        station.
        """
        from .base import RetrialPolicy
        delay_dist, _ = self.get_retrial(jobclass)
        policy = RetrialPolicy.LINEAR
        max_orbit = -1
        if hasattr(self, '_retrial_policies'):
            policy = self._retrial_policies.get(jobclass, RetrialPolicy.LINEAR)
            max_orbit = self._orbit_max_jobs.get(jobclass, -1)
        return delay_dist, policy, max_orbit

    def set_breakdown(self, failure_distribution, repair_distribution,
                      down_service_distribution=None):
        """
        Make the server of this station subject to breakdowns. The server
        alternates between an UP and a DOWN status: while up it fails after
        failure_distribution, while down it is restored after
        repair_distribution.

        The failure clock runs whenever the server is up, whether or not a job
        is in service, so a station can fail while idle. Arrivals are
        unaffected by the server status and keep queueing (subject to the
        station capacity) while the server is down. A job that is in service
        when the server fails is not lost: it stays at the station and, service
        being memoryless in the supported case, resumes when the server is
        repaired.

        Args:
            failure_distribution: time to failure of an up server (Exp)
            repair_distribution: repair time of a down server (Exp)
            down_service_distribution: optional service distribution used while
                the server is down, either a single Distribution applied to
                every class or a list indexed by class. Omitted or empty means
                the server does not serve at all while down, which is the usual
                breakdown model.

        Only exponential failure and repair distributions are currently
        expanded into the joint (queue, server status) chain; anything else is
        rejected here rather than silently approximated.
        """
        from ..distributions.base import Distribution
        from ..distributions.continuous import Exp
        if not isinstance(failure_distribution, Distribution) or \
                not isinstance(repair_distribution, Distribution):
            raise ValueError('set_breakdown requires a failure and a repair Distribution.')
        if not isinstance(failure_distribution, Exp) or not isinstance(repair_distribution, Exp):
            raise ValueError(
                "Station '%s': only exponential failure and repair distributions are supported "
                "by set_breakdown. A non-exponential failure or repair process needs its own "
                "phase in the joint chain, which is not implemented; use an Environment "
                "ensemble for that case." % self.get_name())
        if failure_distribution.getMean() <= 0 or repair_distribution.getMean() <= 0:
            raise ValueError('set_breakdown requires strictly positive failure and repair means.')

        self._breakdown_failure = failure_distribution
        self._breakdown_repair = repair_distribution
        if down_service_distribution is None:
            self._breakdown_down_service = []
        elif isinstance(down_service_distribution, (list, tuple)):
            self._breakdown_down_service = list(down_service_distribution)
        else:
            self._breakdown_down_service = [down_service_distribution]
        self._invalidate_java()

    def get_breakdown(self):
        """
        Return (failure_distribution, repair_distribution,
        down_service_distribution) for this station, or (None, None, []) when
        the station is not subject to breakdowns.
        """
        return (getattr(self, '_breakdown_failure', None),
                getattr(self, '_breakdown_repair', None),
                getattr(self, '_breakdown_down_service', []))

    def has_breakdown(self):
        """True iff a failure and a repair distribution are configured here."""
        return (getattr(self, '_breakdown_failure', None) is not None
                and getattr(self, '_breakdown_repair', None) is not None)

    def set_orbit_impatience(self, jobclass, distribution):
        """
        Set the impatience (abandonment) rate for customers in the orbit.

        This is separate from queue patience (reneging from the waiting queue).
        Used in BMAP/PH/N/N retrial queues where customers in the orbit may
        abandon before successfully retrying.

        Args:
            jobclass: JobClass object
            distribution: Distribution for orbit abandonment time (e.g. Exp(gamma)).
                          Modulated processes (BMAP, MAP, DMAP, MMPP2) are not supported.
        """
        proc_name = type(distribution).__name__
        if proc_name in ('BMAP', 'MAP', 'DMAP', 'MMPP2'):
            raise ValueError(
                "Modulated processes (BMAP, MAP, DMAP, MMPP2) are not supported "
                "for orbit impatience distributions.")
        if not hasattr(self, '_orbit_impatience'):
            self._orbit_impatience = {}
        self._orbit_impatience[jobclass] = distribution
        self._invalidate_java()

    setOrbitImpatience = set_orbit_impatience

    def get_orbit_impatience(self, jobclass):
        """Get the orbit-impatience distribution for a job class (or None)."""
        if not hasattr(self, '_orbit_impatience'):
            return None
        return self._orbit_impatience.get(jobclass)

    getOrbitImpatience = get_orbit_impatience

    def has_orbit_impatience(self, jobclass):
        """Check if a job class has orbit impatience configured at this queue."""
        dist = self.get_orbit_impatience(jobclass)
        if dist is None:
            return False
        return type(dist).__name__ != 'Disabled'

    hasOrbitImpatience = has_orbit_impatience

    def set_batch_reject_probability(self, jobclass, p):
        """
        Set the batch rejection probability for a job class.

        Used in BMAP/PH/N/N retrial queues. When a batch of size k arrives and
        only m < k servers are free:

          - with probability p the entire batch is rejected to the orbit;
          - with probability (1-p) m customers are admitted and k-m go to orbit.

        Args:
            jobclass: JobClass object
            p: Probability in [0,1] that the batch is rejected rather than
               partially admitted. Default is 0 (partial admission allowed).
        """
        p = float(p)
        if p < 0 or p > 1:
            raise ValueError("Batch reject probability must be in [0, 1].")
        if not hasattr(self, '_batch_reject_prob'):
            self._batch_reject_prob = {}
        self._batch_reject_prob[jobclass] = p
        self._invalidate_java()

    setBatchRejectProbability = set_batch_reject_probability

    def get_batch_reject_probability(self, jobclass):
        """Get the batch reject probability for a job class (0 if not set)."""
        if not hasattr(self, '_batch_reject_prob'):
            return 0.0
        return self._batch_reject_prob.get(jobclass, 0.0)

    getBatchRejectProbability = get_batch_reject_probability

    def set_patience(self, jobclass, distribution, impatience_type=None):
        """
        Set patience distribution for a job class.

        Args:
            jobclass: JobClass object
            distribution: Patience distribution
            impatience_type: ImpatienceType enum value (default: RENEGING)
        """
        if not hasattr(self, '_patience_distributions'):
            self._patience_distributions = {}
            self._impatience_types = {}
        from .base import ImpatienceType
        if impatience_type is None:
            impatience_type = ImpatienceType.RENEGING
        self._patience_distributions[jobclass] = distribution
        self._impatience_types[jobclass] = impatience_type
        self._invalidate_java()

    def get_patience(self, jobclass):
        """Get patience distribution for a job class, or None."""
        if not hasattr(self, '_patience_distributions'):
            return None
        return self._patience_distributions.get(jobclass)

    def get_impatience_type(self, jobclass):
        """Get impatience type for a job class, or None."""
        if not hasattr(self, '_impatience_types'):
            return None
        return self._impatience_types.get(jobclass)

    def has_patience(self, jobclass):
        """Check if a job class has patience configured."""
        dist = self.get_patience(jobclass)
        return dist is not None

    def set_immediate_feedback(self, value):
        """
        Configure immediate feedback for this queue.

        When enabled, a job that routes back to this queue (self-loop)
        stays in service instead of re-entering the queue.

        Args:
            value: bool (enable/disable for all classes),
                   JobClass (enable for specific class),
                   or list of JobClass (enable for specific classes)
        """
        if not hasattr(self, '_immediate_feedback_classes'):
            self._immediate_feedback_classes = set()
            self._immediate_feedback_all = False

        if isinstance(value, bool):
            self._immediate_feedback_all = value
            if not value:
                self._immediate_feedback_classes = set()
        elif isinstance(value, list):
            for jc in value:
                self._immediate_feedback_classes.add(jc._index)
        else:
            # Single JobClass
            self._immediate_feedback_classes.add(value._index)
        self._invalidate_java()

    def has_immediate_feedback(self, jobclass=None):
        """
        Check if immediate feedback is enabled.

        Args:
            jobclass: If provided, check for specific class. If None, check if any class has it.
        """
        if not hasattr(self, '_immediate_feedback_classes'):
            return False
        if self._immediate_feedback_all:
            return True
        if jobclass is None:
            return len(self._immediate_feedback_classes) > 0
        return jobclass._index in self._immediate_feedback_classes

    def get_immediate_feedback_classes(self):
        """Get set of class indices with immediate feedback enabled, or 'all'."""
        if not hasattr(self, '_immediate_feedback_classes'):
            return set()
        if self._immediate_feedback_all:
            return 'all'
        return self._immediate_feedback_classes

    # PascalCase aliases
    setPollingType = set_polling_type
    getPollingType = get_polling_type
    setSwitchover = set_switchover
    getSwitchover = get_switchover
    setProbRouting = set_prob_routing
    getProbRouting = get_prob_routing
    setState = set_state
    getState = get_state
    setRoutingWeight = set_routing_weight
    getRoutingWeight = get_routing_weight
    setBalking = set_balking
    getBalking = get_balking
    hasBalking = has_balking
    setRetrial = set_retrial
    getRetrial = get_retrial
    hasRetrial = has_retrial
    setOrbit = set_orbit
    getOrbit = get_orbit
    setBreakdown = set_breakdown
    getBreakdown = get_breakdown
    hasBreakdown = has_breakdown
    setPatience = set_patience
    getPatience = get_patience
    hasPatience = has_patience
    getImpatienceType = get_impatience_type
    setImmediateFeedback = set_immediate_feedback
    hasImmediateFeedback = has_immediate_feedback
    getImmediateFeedbackClasses = get_immediate_feedback_classes


class Delay(Queue):
    """
    Delay node (infinite server).

    A Delay represents an infinite-capacity queue (think time station in closed models).
    All jobs immediately begin service without queueing.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Delay node.

        Args:
            model: Network instance
            name: Delay node name
        """
        super().__init__(model, name, SchedStrategy.INF)
        # Override node_type to DELAY (Queue sets it to QUEUE)
        self._node_type = NodeType.DELAY
        self._number_of_servers = np.inf

    def set_number_of_servers(self, value: int) -> None:
        """
        Set number of servers (always infinity for Delay).

        Args:
            value: Ignored (always set to infinity)

        Raises:
            ValueError: If not infinity
        """
        if not np.isinf(value):
            raise ValueError(f"[{self.name}] Delay node must have infinite servers")
        self._number_of_servers = np.inf
        self._invalidate_java()


class Source(Station):
    """
    Source node for external arrivals.

    A Source represents the entry point for open-class jobs.
    Each network can have at most one Source node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Source node.

        Args:
            model: Network instance
            name: Source node name
        """
        super().__init__(NodeType.SOURCE, name)
        self.set_model(model)
        model.add_node(self)

        self._arrival_process = {}  # Dict[JobClass, Distribution]
        self._number_of_servers = 1  # single server, as in MATLAB/JAR Source
        self._marked_process = None  # shared MarkedMAP driving the marked classes
        self._marked_classes = None  # list of classes bound to marks 1..K
        self._arrival_batch = {}  # Dict[JobClass, DiscreteDistribution]; absent = single arrivals

    def set_arrival(self, jobclass: JobClass, distribution) -> None:
        """
        Set arrival distribution for a job class.

        Args:
            jobclass: Job class
            distribution: Arrival distribution (Exp, Erlang, etc.)
        """
        self._arrival_process[jobclass] = distribution
        self._invalidate_java()

    def set_arrival_batch(self, jobclass: JobClass, batch_size) -> None:
        """
        Set a batch-size law for a class, turning each arrival epoch into the
        simultaneous release of a batch of jobs.

        The interarrival distribution set by :meth:`set_arrival` keeps spacing
        the epochs; this decides how many jobs each epoch releases. Geometric
        interarrivals with a Geometric batch size is the Geo^X arrival stream,
        whose analytical counterpart is
        :func:`line_solver.api.dqsys.dqsys_geoxgeo1`.

        The batch size must be supported on {1,2,...}: an epoch that releases no
        job is not an arrival epoch, so a law that can return zero is rejected
        rather than clamped.

        Args:
            jobclass: Job class
            batch_size: Batch-size distribution, or None to restore single arrivals
        """
        if batch_size is None:
            self._arrival_batch.pop(jobclass, None)
            self._invalidate_java()
            return
        from ..distributions.base import DiscreteDistribution
        if not isinstance(batch_size, DiscreteDistribution):
            raise ValueError("arrival batch size must be a DiscreteDistribution, got %s"
                             % type(batch_size).__name__)
        if batch_size.getMean() < 1.0:
            raise ValueError("arrival batch size for class '%s' has mean %r; a batch must "
                             "carry at least one job, so its support must be {1,2,...}"
                             % (jobclass.name, batch_size.getMean()))
        self._arrival_batch[jobclass] = batch_size
        self._invalidate_java()

    def get_arrival_batch(self, jobclass: JobClass):
        """Return the batch-size law bound to a class, or None for single arrivals."""
        return self._arrival_batch.get(jobclass)

    def set_marked_arrival(self, mmap, classes) -> None:
        """
        Bind a MarkedMAP with K marks to K open classes: mark k emits jobs of
        classes[k-1], with all marks driven by one shared modulating chain.
        Mirrors MATLAB Source.setMarkedArrival.

        Args:
            mmap: MarkedMAP/MMAP arrival process
            classes: list of K distinct OpenClass instances, ordered by mark
        """
        from ..distributions.markovian import MarkedMAP
        if not isinstance(mmap, MarkedMAP):
            raise ValueError("set_marked_arrival requires a MarkedMAP/MMAP arrival process.")
        K = mmap.getNumberOfTypes()
        classes = list(classes)
        if len(classes) != K:
            raise ValueError(
                f"The MarkedMAP has {K} types but {len(classes)} classes were supplied.")
        from .base import JobClassType
        for cls in classes:
            if getattr(cls, '_jobclass_type', None) != JobClassType.OPEN:
                raise ValueError("set_marked_arrival requires open classes.")
        if len({id(c) for c in classes}) != K:
            raise ValueError("set_marked_arrival requires distinct classes for the marks.")
        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        for cls in classes:
            self.set_arrival(cls, mmap)
        self._marked_process = mmap
        self._marked_classes = classes

    def get_marked_process(self):
        """Get the shared MarkedMAP bound via set_marked_arrival, or None."""
        return self._marked_process

    def get_marked_classes(self):
        """Get the classes bound to marks 1..K of the marked process, or None."""
        return self._marked_classes

    setMarkedArrival = set_marked_arrival

    def get_arrival(self, jobclass: JobClass):
        """
        Get arrival distribution for a job class.

        Args:
            jobclass: Job class

        Returns:
            Arrival distribution or None if not set
        """
        return self._arrival_process.get(jobclass)

    def get_arrival_process(self, jobclass: JobClass):
        """Get the arrival process (distribution) for a job class."""
        return self._arrival_process.get(jobclass)

    getArrivalProcess = get_arrival_process
    arrival_process = get_arrival_process

    def get_arrival_rates(self) -> Dict[JobClass, float]:
        """
        Get arrival rates (1/mean) for all classes.

        Returns:
            Dict mapping classes to arrival rates
        """
        rates = {}
        for jobclass, dist in self._arrival_process.items():
            if dist is not None and hasattr(dist, 'getMean'):
                rates[jobclass] = 1.0 / dist.getMean()
        return rates

    # PascalCase aliases for MATLAB compatibility
    setArrival = set_arrival
    getArrival = get_arrival
    getArrivalRates = get_arrival_rates
    setArrivalBatch = set_arrival_batch
    getArrivalBatch = get_arrival_batch


class Sink(Node):
    """
    Sink node for job departure.

    A Sink is the exit point for open-class jobs.
    Each network can have at most one Sink node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Sink node.

        Args:
            model: Network instance
            name: Sink node name
        """
        super().__init__(NodeType.SINK, name)
        self.set_model(model)
        model.add_node(self)


class Fork(Node):
    """
    Fork node for parallel processing.

    A Fork splits an incoming job into multiple parallel tasks
    that must be synchronized at a corresponding Join node.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Fork node.

        Args:
            model: Network instance
            name: Fork node name
        """
        super().__init__(NodeType.FORK, name)
        self.set_model(model)
        model.add_node(self)

        self._tasks_per_link = None
        # Variable forking levels. Each is a list of (dest_name, class_idx,
        # payload); dest_name '' means every connected link. Kept by
        # DESTINATION NAME rather than link ordinal because the link order is an
        # artefact of connmatrix traversal and renumbers on a relink.
        self._tasks_per_link_by_dest = []
        self._tasks_per_link_dist = []
        self._branch_prob = []
        self._capacity = np.inf
        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        self._state = None
        self._state_prior = None
        self._state_space = None

    def get_state(self):
        """Get Fork state (FJ tag-augmented copies only)."""
        return self._state

    def set_state(self, state) -> None:
        """Set Fork state (FJ tag-augmented copies only)."""
        self._state = state

    def get_state_prior(self):
        """Get Fork state prior (FJ tag-augmented copies only)."""
        return self._state_prior

    def set_state_prior(self, prior) -> None:
        """Set Fork state prior (FJ tag-augmented copies only)."""
        self._state_prior = np.asarray(prior).ravel()

    setStatePrior = set_state_prior

    def get_state_space(self):
        """Get Fork state space (FJ tag-augmented copies only)."""
        return self._state_space

    def set_state_space(self, space) -> None:
        """Set Fork state space (FJ tag-augmented copies only)."""
        self._state_space = space

    def set_tasks_per_link(self, ntasks, jobclass=None, dest_node=None) -> None:
        """
        Set number of tasks spawned per outgoing link.

        The total number of siblings a firing creates is
        (number of outgoing links) * ntasks. SolverJMT and SolverLDES simulate
        it directly; the MMT fork-join transform behind SolverMVA/SolverNC
        carries the load of all the siblings on the auxiliary open class and
        synchronises on the order statistic of that many branch times, each
        branch replicated ntasks times. The H-T transform refuses ntasks > 1.

        Args:
            ntasks: Task count emitted on each outgoing link
            jobclass: Optional job class; when given, only that class is
                affected and every other class keeps the node-wide value
            dest_node: Optional destination node; when given, only the link
                towards it is affected and the other links are left alone
        """
        if jobclass is not None:
            dest = '' if dest_node is None else dest_node.get_name()
            self._tasks_per_link_by_dest.append(
                (dest, jobclass.get_index0(), float(ntasks)))
            self._invalidate_java()
            return
        self._tasks_per_link = np.array(ntasks)
        self._invalidate_java()

    def get_tasks_per_link(self) -> Optional[np.ndarray]:
        """Get tasks per link."""
        return self._tasks_per_link

    def set_tasks_per_link_distribution(self, jobclass, dist, dest_node=None) -> None:
        """
        Make the number of tasks emitted on each outgoing link RANDOM.

        The degree is redrawn independently for every link and every forked job,
        which is the variable forking level of JMT's JobsPerLinkDis. Exact under
        SolverJMT and SolverLDES, which draw it at the fork epoch; the analytical
        solvers see E[dist].

        Args:
            jobclass: Job class the distribution applies to
            dist: DiscreteSampler over the tasks-per-link support
            dest_node: Optional destination node restricting it to one link
        """
        from ..distributions.discrete import DiscreteDistribution
        if not isinstance(dist, DiscreteDistribution):
            raise ValueError('The tasks-per-link distribution must be a '
                             'DiscreteDistribution (e.g. DiscreteSampler).')
        dest = '' if dest_node is None else dest_node.get_name()
        self._tasks_per_link_dist.append((dest, jobclass.get_index0(), dist))
        self._invalidate_java()

    def set_branch_probability(self, jobclass, dest_node, prob: float) -> None:
        """
        Activate an outgoing branch only with probability `prob`.

        Branches are activated independently, so the number of siblings a job
        produces is random even when the tasks per link are deterministic. The
        matched Join must be told what to wait for: under a standard join a job
        that skipped a branch would block forever, so any probability below one
        requires JoinStrategy.PARTIAL or a quorum.
        """
        if prob < 0.0 or prob > 1.0:
            raise ValueError('A branch activation probability must lie in [0,1].')
        self._branch_prob.append((dest_node.get_name(), jobclass.get_index0(), float(prob)))
        self._invalidate_java()

    def get_tasks_per_link_distribution(self):
        """Get the registered per-link jobs-per-link distributions."""
        return self._tasks_per_link_dist

    def get_branch_probability(self):
        """Get the registered branch activation probabilities."""
        return self._branch_prob

    setState = set_state
    getState = get_state
    getStatePrior = get_state_prior
    setStateSpace = set_state_space
    getStateSpace = get_state_space
    setTasksPerLink = set_tasks_per_link
    getTasksPerLink = get_tasks_per_link
    setTasksPerLinkDistribution = set_tasks_per_link_distribution
    getTasksPerLinkDistribution = get_tasks_per_link_distribution
    setBranchProbability = set_branch_probability
    getBranchProbability = get_branch_probability

    @property
    def capacity(self) -> float:
        """Get capacity."""
        return self._capacity

    @capacity.setter
    def capacity(self, value: float) -> None:
        """Set capacity."""
        self._capacity = value
        self._invalidate_java()


class Join(Station):
    """
    Join node for parallel processing synchronization.

    A Join waits for all parallel tasks from a corresponding Fork
    before releasing the job downstream.
    """

    def __init__(self, model, name: str, fork: Optional[Fork] = None):
        """
        Initialize a Join node.

        Args:
            model: Network instance
            name: Join node name
            fork: Associated Fork node (optional)
        """
        super().__init__(NodeType.JOIN, name)
        self.set_model(model)
        model.add_node(self)

        self._fork = fork
        self._join_strategy = {}  # Dict[JobClass, JoinStrategy]
        self._required = {}  # Dict[JobClass, required_count] (-1 = all)
        self._number_of_servers = np.inf

    def set_fork(self, fork: Fork) -> None:
        """
        Set associated Fork node.

        Args:
            fork: Fork node to pair with
        """
        self._fork = fork
        self._invalidate_java()

    def get_fork(self) -> Optional[Fork]:
        """Get associated Fork node."""
        return self._fork

    def set_strategy(self, jobclass: JobClass, strategy: JoinStrategy) -> None:
        """
        Set join strategy for a job class.

        Args:
            jobclass: Job class
            strategy: Join strategy (STD, QUORUM, etc.)
        """
        self._join_strategy[jobclass] = strategy
        self._invalidate_java()

    def get_strategy(self, jobclass: JobClass) -> JoinStrategy:
        """Get join strategy for a job class."""
        return self._join_strategy.get(jobclass, JoinStrategy.STD)

    def set_required(self, jobclass: JobClass, nrequired: int) -> None:
        """
        Set required number of tasks to wait for (quorum join).

        Args:
            jobclass: Job class
            nrequired: Number of required tasks (-1 = all)
        """
        self._required[jobclass] = nrequired
        self._invalidate_java()

    def get_required(self, jobclass: JobClass) -> int:
        """Get required task count."""
        return self._required.get(jobclass, -1)


class Router(Node):
    """
    Router node for routing decisions.

    A Router is a node that routes jobs without service delay.
    Can be used for probabilistic routing and routing-dependent decisions.
    """

    def __init__(self, model, name: str):
        """
        Initialize a Router node.

        Args:
            model: Network instance
            name: Router node name
        """
        super().__init__(NodeType.ROUTER, name)
        self.set_model(model)
        model.add_node(self)


class ClassSwitch(Node):
    """
    ClassSwitch node for job class switching.

    A ClassSwitch allows jobs to change class without service delay.
    Useful for modeling class-dependent routing and scheduling.
    """

    def __init__(self, model, name: str, cs_matrix=None):
        """
        Initialize a ClassSwitch node.

        Args:
            model: Network instance
            name: ClassSwitch node name
            cs_matrix: Optional K×K class switching probability matrix where
                      element (i,j) is the probability that a job in class i
                      switches to class j. If not provided, must be set later
                      using set_class_switching_matrix().
        """
        super().__init__(NodeType.CLASSSWITCH, name)
        self.set_model(model)
        model.add_node(self)

        # Set switching matrix if provided
        if cs_matrix is not None:
            self._switch_matrix = np.asarray(cs_matrix)
        else:
            self._switch_matrix = None  # Class switching probabilities

    def has_class_switching(self) -> bool:
        """A ClassSwitch node always switches classes (MATLAB: its server is a
        StatelessClassSwitcher)."""
        return True

    def init_class_switch_matrix(self) -> np.ndarray:
        """
        Initialize and return a class switching matrix.

        Creates a K×K matrix of zeros where K is the number of classes.
        Use this to create a template that can be filled with switching
        probabilities before calling set_class_switching_matrix().

        Returns:
            np.ndarray: K×K matrix initialized to zeros

        Example:
            >>> csmatrix = cs_node.init_class_switch_matrix()
            >>> csmatrix[0, 1] = 0.3  # 30% switch from class 0 to class 1
            >>> csmatrix[0, 0] = 0.7  # 70% stay in class 0
            >>> csmatrix[1, 0] = 1.0  # 100% switch from class 1 to class 0
            >>> cs_node.set_class_switching_matrix(csmatrix)
        """
        K = self._model.get_number_of_classes()
        return np.zeros((K, K))

    initClassSwitchMatrix = init_class_switch_matrix
    init_class_switching_matrix = init_class_switch_matrix  # Alternative alias

    def set_class_switching_matrix(self, cs_matrix: np.ndarray):
        """
        Set the class switching probability matrix.

        Args:
            cs_matrix: K×K matrix where element (i,j) is the probability
                      that a job in class i switches to class j.
                      Each row should sum to 1.0.

        Example:
            >>> csmatrix = cs_node.init_class_switch_matrix()
            >>> csmatrix[0, 0] = 0.7
            >>> csmatrix[0, 1] = 0.3
            >>> cs_node.set_class_switching_matrix(csmatrix)
        """
        self._switch_matrix = np.asarray(cs_matrix)

    setClassSwitchingMatrix = set_class_switching_matrix

    def get_class_switching_matrix(self) -> Optional[np.ndarray]:
        """
        Get the current class switching matrix.

        Returns:
            np.ndarray or None: The K×K switching matrix, or None if not set
        """
        return self._switch_matrix

    getClassSwitchingMatrix = get_class_switching_matrix

    def set_switch_probability(self, from_class, to_class, probability: float):
        """
        Set the probability of switching from one class to another.

        This is a convenience method that handles indexing automatically.

        Args:
            from_class: Source job class (object or 0-based index)
            to_class: Target job class (object or 0-based index)
            probability: Switching probability (0.0 to 1.0)

        Example:
            >>> cs_node.set_switch_probability(class1, class2, 0.3)
            >>> cs_node.set_switch_probability(class1, class1, 0.7)
        """
        if self._switch_matrix is None:
            K = self._model.get_number_of_classes()
            self._switch_matrix = np.zeros((K, K))

        # Handle class objects or indices
        if hasattr(from_class, 'get_index0'):
            from_idx = from_class.get_index0()
        elif hasattr(from_class, 'get_index'):
            from_idx = from_class.get_index() - 1  # Convert 1-based to 0-based
        else:
            from_idx = int(from_class)

        if hasattr(to_class, 'get_index0'):
            to_idx = to_class.get_index0()
        elif hasattr(to_class, 'get_index'):
            to_idx = to_class.get_index() - 1  # Convert 1-based to 0-based
        else:
            to_idx = int(to_class)

        self._switch_matrix[from_idx, to_idx] = probability

    setSwitchProbability = set_switch_probability

    def remove_job_class(self, jobclass) -> None:
        """
        Drop a job class from this node and from its class-switching matrix.

        The matrix is indexed by class position, not by class object, so the
        removed row and column are deleted while the model still holds the
        class (Network.remove_class updates the nodes before dropping the
        class). Mirrors ClassSwitch.removeJobClass in the JAR.

        Args:
            jobclass: the job class being removed from the model
        """
        super().remove_job_class(jobclass)
        if self._switch_matrix is None or self._model is None:
            return
        classes = self._model.get_classes()
        r = next((i for i, c in enumerate(classes) if c is jobclass), -1)
        if r < 0:
            return
        remaining = [i for i in range(len(classes)) if i != r]
        matrix = np.asarray(self._switch_matrix)
        if matrix.shape[0] <= max(remaining, default=-1) or matrix.shape[1] <= max(remaining, default=-1):
            return
        self._switch_matrix = matrix[np.ix_(remaining, remaining)]

    removeJobClass = remove_job_class


class Cache(Station):
    """
    Multi-level cache node with configurable replacement strategy.

    A Cache node models content caching behavior where arriving jobs request
    items from a catalog. Items may be cached (hit) or fetched (miss), with
    jobs potentially switching classes based on hit/miss outcomes.

    The cache supports multiple levels (hierarchy) with different capacities
    and uses configurable replacement policies (LRU, FIFO, RR, SFIFO).

    Cache analysis can be performed using algorithms from api.cache:
    - cache_mva: Mean Value Analysis for hit/miss probabilities
    - cache_erec: Exact recursive computation
    - cache_spm: Singular Perturbation Method (approximate)

    Args:
        model: Network instance
        name: Cache node name
        num_items: Number of items in the catalog (n)
        item_level_cap: Capacity of each cache level (int or array)
        replacement_strategy: Cache replacement policy (LRU, FIFO, RR, SFIFO)

    Example:
        >>> model = Network('CacheModel')
        >>> cache = Cache(model, 'MyCache', num_items=100,
        ...               item_level_cap=10, replacement_strategy=ReplacementStrategy.LRU)
        >>> job_class = ClosedClass(model, 'Request', 5, delay)
        >>> hit_class = ClosedClass(model, 'Hit', 0, delay)
        >>> miss_class = ClosedClass(model, 'Miss', 0, delay)
        >>> cache.set_read(job_class, Zipf(1.2, 100))  # Zipf popularity
        >>> cache.set_hit_class(job_class, hit_class)
        >>> cache.set_miss_class(job_class, miss_class)

    Reference:
        - Cache algorithms: Che, H. et al. "Hierarchical Web Caching Systems"
        - TTL approximation: Fofack et al. "Analysis of TTL-based Cache Networks"
    """

    def __init__(
        self,
        model,
        name: str,
        num_items: int,
        item_level_cap: Union[int, np.ndarray, list],
        replacement_strategy: ReplacementStrategy = ReplacementStrategy.LRU
    ):
        """
        Initialize a Cache node.

        Args:
            model: Network instance
            name: Cache node name
            num_items: Total number of items in catalog (n >= 1)
            item_level_cap: Capacity per cache level. Can be:
                - int: Single-level cache with given capacity
                - array/list: Multi-level cache with capacity per level
            replacement_strategy: Replacement policy (default: LRU)

        Raises:
            ValueError: If num_items < 1 or item_level_cap invalid
        """
        super().__init__(NodeType.CACHE, name)
        self.set_model(model)
        model.add_node(self)

        # Validate parameters
        if num_items < 1:
            raise ValueError(f"num_items must be >= 1, got {num_items}")

        self._num_items = int(num_items)

        # Convert item_level_cap to array
        if isinstance(item_level_cap, (int, float)):
            self._item_level_cap = np.array([int(item_level_cap)])
        else:
            self._item_level_cap = np.array(item_level_cap, dtype=int).ravel()

        if np.any(self._item_level_cap < 1):
            raise ValueError("All cache level capacities must be >= 1")

        self._num_levels = len(self._item_level_cap)
        self._replacement_strategy = replacement_strategy
        self._admission_prob = 1.0

        # per-item storage costs and per-list cost caps (ton21cache Sec. IX)
        self._item_size = None
        self._cost_cap = None
        self._cost_cap_global = False

        # Cache nodes have infinite servers (instant service)
        self._number_of_servers = np.inf

        # Scheduling (caches are non-preemptive, use INF scheduling like Delay)
        self._sched_strategy = SchedStrategy.INF
        self._sched_policy = SchedStrategyType.NP

        # Hit/miss class mappings: input_class -> output_class
        self._hit_class = {}   # Dict[JobClass, JobClass]
        self._miss_class = {}  # Dict[JobClass, JobClass]

        # Read (popularity) distributions per class
        self._read_process = {}  # Dict[JobClass, Distribution]
        self._item_classes = {}   # Dict[JobClass, List[JobClass]] keyed by the chain's first read class
        self._item_of_class = {}  # Dict[JobClass, int] item each per-item class reads (1-based)

        # Access probabilities matrix (optional)
        self._access_prob = None  # (num_items, num_classes)
        self._graph = None        # per-item (h+1)x(h+1) access-cost graph (set_access_graph)
        self._accost = None       # [class][item] (h+1)x(h+1) access-cost matrices

        # Access cost matrix (optional)
        self._access_cost = None  # (num_items, num_classes)

        # Result storage (populated by solver)
        self._actual_hit_prob = None   # Actual hit probabilities
        self._actual_miss_prob = None  # Actual miss probabilities
        self._actual_delayed_hit_prob = None  # Actual delayed-hit fractions (retrieval system)
        self._actual_hit_prob_list = None  # Actual per-list (per-level) hit fractions [classes x lists]
        self._actual_item_prob = None  # Actual per-item occupancy [items x (lists+1)]: col 0 = miss, cols 1.. = per-list
        self._actual_residt = None  # Actual expected latency
        self._actual_list_cost = None  # Mean storage cost held by each list [1 x lists]
        self._actual_delayed_hit_qlen = None       # [1 x items] merged secondary requests per item
        self._actual_delayed_hit_qlen_full = None  # as above, plus the triggering request

        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        self._total_cache_capacity = int(np.sum(self._item_level_cap))
        self._retrieval_system_capacity = 0  # 0 until set_retrieval_system is called
        # (item, JobClass) -> JobClass  (item is 0-based)
        self._retrieval_classes = {}
        # 0-based indices of the retrieval classes
        self._retrieval_class_indices = set()
        # arrival-class 0-based index -> list of retrieval-system queue node indices
        self._retrieval_system_queue_indices = {}
        # deferred routing edges [fromCls, toCls, srcNode, dstNode, prob] (1-based
        # class / node indices) injected into the routing matrix P by link()
        self._retrieval_routing_entries = []

    def is_station(self) -> bool:
        """
        Check if node is a station.

        Cache nodes are NOT stations in LINE's semantics, matching MATLAB behavior.
        They are stateful nodes that immediately process arriving jobs (class switching
        based on cache hit/miss) but don't queue or serve jobs like stations.

        Returns:
            False: Cache is a stateful non-station node
        """
        return False

    def has_class_switching(self) -> bool:
        """A Cache switches the arriving class into its hit or miss class
        (MATLAB: its server is a CacheClassSwitcher, itself a ClassSwitcher)."""
        return True

    # =========================================================================
    # Properties
    # =========================================================================

    @property
    def num_items(self) -> int:
        """Get the number of items in the catalog."""
        return self._num_items

    @property
    def num_levels(self) -> int:
        """Get the number of cache levels."""
        return self._num_levels

    @property
    def item_level_cap(self) -> np.ndarray:
        """Get the capacity array for each cache level."""
        return self._item_level_cap.copy()

    @property
    def total_capacity(self) -> int:
        """Get the total cache capacity across all levels."""
        return int(np.sum(self._item_level_cap))

    @property
    def replacement_strategy(self) -> ReplacementStrategy:
        """Get the cache replacement strategy."""
        return self._replacement_strategy

    @replacement_strategy.setter
    def replacement_strategy(self, value: ReplacementStrategy):
        """Set the cache replacement strategy."""
        self._replacement_strategy = value
        self._invalidate_java()

    def remove_job_class(self, jobclass) -> None:
        """
        Reject class removal on a cache node.

        The cache item state is indexed by class and is carried in the model
        state, not only in the node configuration, so a class cannot be dropped
        without invalidating it. Matches MATLAB @MNetwork/removeClass.m and
        Cache.removeJobClass in the JAR, which both refuse.

        Args:
            jobclass: the job class being removed from the model

        Raises:
            RuntimeError: always
        """
        raise RuntimeError(
            "Cannot dynamically remove classes in models with caches. "
            "You need to re-instantiate the model.")

    removeJobClass = remove_job_class

    # =========================================================================
    # Hit/Miss Class Configuration
    # =========================================================================

    def set_hit_class(self, input_class: JobClass, output_class: JobClass) -> None:
        """
        Set the output class for cache hits.

        When a job of input_class experiences a cache hit, it transitions
        to output_class for further processing.

        Args:
            input_class: Incoming job class making the request
            output_class: Resulting job class after cache hit
        """
        self._hit_class[input_class] = output_class
        self._invalidate_java()

    def get_hit_class(self, input_class: JobClass) -> Optional[JobClass]:
        """Get the output class for cache hits."""
        return self._hit_class.get(input_class)

    def set_miss_class(self, input_class: JobClass, output_class: JobClass) -> None:
        """
        Set the output class for cache misses.

        When a job of input_class experiences a cache miss, it transitions
        to output_class for further processing.

        Args:
            input_class: Incoming job class making the request
            output_class: Resulting job class after cache miss
        """
        self._miss_class[input_class] = output_class
        self._invalidate_java()

    def get_miss_class(self, input_class: JobClass) -> Optional[JobClass]:
        """Get the output class for cache misses."""
        return self._miss_class.get(input_class)

    # =========================================================================
    # Read (Popularity) Distribution Configuration
    # =========================================================================

    def set_read(self, jobclass: JobClass, distribution) -> None:
        """
        Set the read (popularity) distribution for a job class.

        The distribution determines which items are requested by jobs
        of the given class. Typically a discrete distribution like Zipf.

        Args:
            jobclass: Job class making requests
            distribution: Popularity distribution (e.g., Zipf, Uniform)
        """
        self._read_process[jobclass] = distribution
        self._invalidate_java()

    @staticmethod
    def _per_item_classes(spec, n_items: int, what: str):
        """One class shared by every item, or a sequence of one class per item."""
        if isinstance(spec, (list, tuple)):
            if len(spec) != n_items:
                raise ValueError(
                    f"{len(spec)} {what} classes were given for {n_items} items; "
                    "pass one class per item or a single class shared by all.")
            return list(spec)
        return [spec] * n_items

    def set_item_read_classes(self, read_classes, hit_classes) -> None:
        """
        Declare that read_classes[i] is the request stream for item i at this cache.

        Use at the cache the exogenous requests enter, where the per-item classes are
        the user's own; item popularity is then carried by the per-class request rates
        rather than by a popularity distribution the cache draws from. This is what
        keeps a cache network free of arc-level class switching, so no class acquires a
        default route into the cache that the model never intended.
        """
        from ..distributions.discrete import DiscreteSampler
        n_items = self._num_items
        read_arr = list(read_classes)
        if len(read_arr) != n_items:
            raise ValueError(
                f"{self.name} holds {n_items} items but {len(read_arr)} read classes "
                "were given; pass exactly one class per item.")
        hit_arr = Cache._per_item_classes(hit_classes, n_items, "hit")
        for i, cls in enumerate(read_arr):
            item_popularity = np.zeros(n_items)
            item_popularity[i] = 1.0
            self._read_process[cls] = DiscreteSampler(item_popularity)
            self.set_hit_class(cls, hit_arr[i])
            self._item_of_class[cls] = i + 1
        self._item_classes[read_arr[0]] = read_arr
        self._invalidate_java()

    def set_item_classes(self, jobin_class: JobClass, hit_classes):
        """
        Mint one class per item at a cache FED BY ANOTHER CACHE, so item identity
        survives the miss hop. Idempotent.
        """
        from ..distributions.discrete import DiscreteSampler
        from .classes import OpenClass, ClosedClass
        existing = self._item_classes.get(jobin_class)
        if existing is not None:
            return existing
        model = self.get_model()
        n_items = self._num_items
        hit_arr = Cache._per_item_classes(hit_classes, n_items, "hit")
        is_closed = isinstance(jobin_class, ClosedClass)
        ref_station = getattr(jobin_class, '_refstat', None)
        minted = []
        for i in range(n_items):
            if is_closed:
                cls = ClosedClass(model, f"{self.name}_item{i + 1}", 0, ref_station, 0)
            else:
                cls = OpenClass(model, f"{self.name}_item{i + 1}", 0)
            item_popularity = np.zeros(n_items)
            item_popularity[i] = 1.0
            self._read_process[cls] = DiscreteSampler(item_popularity)
            self.set_hit_class(cls, hit_arr[i])
            self._item_of_class[cls] = i + 1
            minted.append(cls)
        self._item_classes[jobin_class] = minted
        self._invalidate_java()
        return minted

    def set_miss_cache(self, jobin_class: JobClass, next_cache, hit_class_at_next):
        """
        Send this cache's misses to next_cache preserving item identity: the miss class
        of this cache for item i IS the read class of next_cache for item i. Returns the
        classes minted on next_cache so the caller can route them onward.
        """
        if next_cache._num_items != self._num_items:
            raise ValueError(
                f"Cache {self.name} holds {self._num_items} items but "
                f"{next_cache.name} holds {next_cache._num_items}; "
                "a cache network requires one common item set.")
        self_cls = self._item_classes.get(jobin_class)
        if self_cls is None:
            raise ValueError(
                f"No per-item classes at {self.name} for class {jobin_class.name}; "
                "call set_item_read_classes before set_miss_cache.")
        next_cls = next_cache.set_item_classes(jobin_class, hit_class_at_next)
        for i, cls in enumerate(self_cls):
            self.set_miss_class(cls, next_cls[i])
            self._add_retrieval_routing_entry(
                next_cls[i].get_index0(), next_cls[i].get_index0(),
                self.get_index0(), next_cache.get_index0(), 1.0)
        return next_cls

    def set_item_miss_class(self, jobin_class: JobClass, miss_classes) -> None:
        """
        Terminate a cache network: every per-item class of this cache reports a miss as
        the matching entry of miss_classes, which the user routes onward.
        """
        self_cls = self._item_classes.get(jobin_class)
        if self_cls is None:
            raise ValueError(
                f"No per-item classes at {self.name} for class {jobin_class.name}; "
                "call set_item_classes before set_item_miss_class.")
        miss_arr = Cache._per_item_classes(miss_classes, len(self_cls), "miss")
        for i, cls in enumerate(self_cls):
            self.set_miss_class(cls, miss_arr[i])

    def set_item_of_class(self, jobclass: JobClass, item: int) -> None:
        """Record that a class reads a given item (1-based), as read back from JSON."""
        self._item_of_class[jobclass] = int(item)

    def get_item_of_class(self, jobclass: JobClass) -> int:
        """Item a per-item class reads (1-based), 0 when the class is not one."""
        return self._item_of_class.get(jobclass, 0)

    def get_read(self, jobclass: JobClass):
        """Get the read distribution for a job class."""
        return self._read_process.get(jobclass)

    # =========================================================================
    # Access Probabilities (Alternative to Distribution)
    # =========================================================================

    def set_access_prob(self, access_prob: np.ndarray) -> None:
        """
        Set the access probability matrix directly.

        Alternative to set_read() for specifying item access probabilities
        as a matrix.

        Args:
            access_prob: (num_items, num_classes) probability matrix.
                         Each column should sum to 1.
        """
        self._access_prob = np.asarray(access_prob)
        self._invalidate_java()

    def set_admission_prob(self, q: float) -> None:
        """
        Set the q-LRU admission probability.

        Probability q in [0, 1] of admitting a missed item into the cache.
        Only used when the replacement strategy is QLRU.
        """
        if q < 0 or q > 1:
            raise ValueError("The admission probability q must lie in [0,1].")
        self._admission_prob = float(q)
        self._invalidate_java()

    def set_item_sizes(self, sizes) -> None:
        """
        Set the storage cost (size) of each item.

        A positive integer vector with one entry per item, or a single value
        applied to every item. Used together with set_cost_caps to bound the
        storage held by each cache list.
        """
        sizes = np.atleast_1d(np.asarray(sizes, dtype=float)).ravel()
        if sizes.size == 1:
            sizes = np.full(self._num_items, sizes[0])
        if sizes.size != self._num_items:
            raise ValueError(f"The item size vector of {self.get_name()} must have "
                             f"one entry per item ({self._num_items}).")
        if np.any(sizes <= 0) or np.any(sizes != np.round(sizes)):
            raise ValueError("Item sizes must be positive integers.")
        self._item_size = sizes
        self._invalidate_java()

    def get_item_sizes(self) -> Optional[np.ndarray]:
        """Get the per-item storage costs (sizes), None when unset."""
        return self._item_size

    def set_cost_caps(self, caps) -> None:
        """
        Set the per-list cap on the total storage cost of the resident items.

        A single value declares one cap for the whole cache, modelled as the
        same cap on every list.
        """
        caps = np.atleast_1d(np.asarray(caps, dtype=float)).ravel()
        h = len(self._item_level_cap)
        if caps.size == 1:
            self._cost_cap_global = True
            caps = np.full(h, caps[0])
        else:
            self._cost_cap_global = False
        if caps.size != h:
            raise ValueError(f"The cost cap vector of {self.get_name()} must have "
                             f"one entry per cache list ({h}).")
        if np.any(caps < 0) or np.any(caps != np.round(caps)):
            raise ValueError("Storage cost caps must be non-negative integers.")
        if self._replacement_strategy == ReplacementStrategy.CLIMB:
            raise ValueError("Storage cost caps are not supported with the CLIMB "
                             "replacement strategy.")
        self._cost_cap = caps
        self._invalidate_java()

    def get_cost_caps(self) -> Optional[np.ndarray]:
        """Get the per-list storage cost caps, None when unset."""
        return self._cost_cap

    def is_cost_cap_global(self) -> bool:
        """Report whether the cost caps came from a single cache-wide cap."""
        return self._cost_cap_global

    def set_result_list_cost(self, list_cost) -> None:
        """Store the mean storage cost held by each list, as computed by a solver."""
        self._actual_list_cost = list_cost

    def get_list_cost(self) -> Optional[np.ndarray]:
        """Get the mean storage cost held by each list, None when not computed."""
        return self._actual_list_cost

    def get_access_prob(self) -> Optional[np.ndarray]:
        """Get the access probability matrix."""
        return self._access_prob

    def set_access_graph(self, graph) -> None:
        """
        Set the per-item access-cost (list-move) graph.

        Mirrors the MATLAB Cache constructor's graph argument: one
        (h+1) x (h+1) matrix per item, shared by all job classes, where h is
        the number of cache lists. Row 1 governs miss insertion (column 1 =
        do not cache, column 1+l = insert into list l); row 1+i governs hits
        in list i (column 1+j = move the item to list j >= i).

        Args:
            graph: sequence of num_items matrices, each (h+1) x (h+1).
        """
        h = self._num_levels
        mats = [np.asarray(g, dtype=float) for g in graph]
        if len(mats) != self._num_items:
            raise ValueError(
                f"Access graph must have one matrix per item ({self._num_items}), got {len(mats)}")
        for g in mats:
            if g.shape != (h + 1, h + 1):
                raise ValueError(
                    f"Each access graph matrix must be ({h + 1}, {h + 1}), got {g.shape}")
        self._graph = mats
        self._invalidate_java()

    def get_access_graph(self):
        """Return the per-item access-cost graph, or None if not set."""
        return getattr(self, '_graph', None)

    # =========================================================================
    # Result Methods (Populated by Solver)
    # =========================================================================

    def set_result_hit_prob(self, hit_prob: np.ndarray) -> None:
        """Set the computed hit probability (called by solver)."""
        self._actual_hit_prob = np.asarray(hit_prob)

    def set_result_miss_prob(self, miss_prob: np.ndarray) -> None:
        """Set the computed miss probability (called by solver)."""
        self._actual_miss_prob = np.asarray(miss_prob)

    def set_result_delayed_hit_prob(self, delayed_hit_prob: np.ndarray) -> None:
        """Set the computed delayed-hit fraction (retrieval system; called by solver)."""
        self._actual_delayed_hit_prob = np.asarray(delayed_hit_prob)

    def set_result_hit_prob_list(self, hit_prob_list: np.ndarray) -> None:
        """Set the computed per-list (per-level) hit fractions [classes x lists]."""
        self._actual_hit_prob_list = np.asarray(hit_prob_list)

    def get_hit_ratio(self) -> Optional[np.ndarray]:
        """
        Get the hit ratio (probability) for each class.

        Returns:
            Hit probability array, or None if not computed yet.
        """
        return self._actual_hit_prob

    def get_delayed_hit_ratio(self) -> Optional[np.ndarray]:
        """Get the delayed-hit fraction per class (None when no retrieval system)."""
        return self._actual_delayed_hit_prob

    def get_hit_ratio_by_list(self) -> Optional[np.ndarray]:
        """Get the per-class, per-list hit fraction matrix [classes x lists]."""
        return self._actual_hit_prob_list

    def set_result_item_prob(self, item_prob: np.ndarray) -> None:
        """Set the per-item occupancy [items x (lists+1)]; col 0 = miss, cols 1.. = per-list."""
        self._actual_item_prob = np.asarray(item_prob)

    def get_item_prob(self) -> Optional[np.ndarray]:
        """Get the per-item occupancy matrix [items x (lists+1)] (col 0 = miss)."""
        return self._actual_item_prob

    def get_miss_ratio(self) -> Optional[np.ndarray]:
        """
        Get the miss ratio (probability) for each class.

        Returns:
            Miss probability array, or None if not computed yet.
        """
        return self._actual_miss_prob

    # =========================================================================
    # Utility Methods
    # =========================================================================

    def get_gamma_matrix(self, nclasses: int = 1) -> np.ndarray:
        """
        Build the gamma (access factor) matrix for cache analysis.

        The gamma matrix has shape (num_items, num_levels) and contains
        the cache access intensity factors used by cache_mva/cache_erec.

        Args:
            nclasses: Number of job classes

        Returns:
            Gamma matrix for cache analysis algorithms
        """
        n = self._num_items
        h = self._num_levels

        gamma = np.zeros((n, h))

        # If access_prob is set, use it to derive gamma
        if self._access_prob is not None:
            # Aggregate across classes
            access_agg = np.sum(self._access_prob, axis=1) if self._access_prob.ndim > 1 else self._access_prob
            # Normalize to get access intensities
            total = np.sum(access_agg)
            if total > 0:
                gamma[:, 0] = access_agg / total
        elif self._read_process:
            # Extract probabilities from read distributions (e.g., Zipf)
            access_probs = np.zeros(n)
            for jobclass, dist in self._read_process.items():
                if dist is not None:
                    # Check if it's a Zipf distribution
                    if hasattr(dist, '_s') and hasattr(dist, '_H'):
                        # Zipf distribution: P(k) = (1/k^s) / H
                        k = np.arange(1, n + 1)
                        probs = (1.0 / k ** dist._s) / dist._H
                        access_probs += probs
                    # Check if it's a DiscreteSampler with probabilities
                    elif hasattr(dist, '_probs'):
                        probs = np.asarray(dist._probs)
                        if len(probs) >= n:
                            access_probs += probs[:n]
                        else:
                            access_probs[:len(probs)] += probs
                    # Check if it has evalPMF method (generic discrete distribution)
                    elif hasattr(dist, 'evalPMF'):
                        for i in range(n):
                            access_probs[i] += dist.evalPMF(i + 1)  # 1-indexed items
            # Normalize
            total = np.sum(access_probs)
            if total > 0:
                gamma[:, 0] = access_probs / total
            else:
                gamma[:, 0] = 1.0 / n
        else:
            # Default: uniform access
            gamma[:, 0] = 1.0 / n

        # For multi-level caches, propagate to other levels
        for level in range(1, h):
            gamma[:, level] = gamma[:, level - 1]

        return gamma

    def get_capacity_vector(self) -> np.ndarray:
        """
        Get the cache capacity vector.

        Returns:
            Array of cache level capacities
        """
        return self._item_level_cap.copy()

    # =========================================================================
    # Retrieval System (Delayed Hits)
    # =========================================================================

    @property
    def total_cache_capacity(self) -> int:
        """Total cache capacity (sum of per-level capacities)."""
        return self._total_cache_capacity

    @property
    def retrieval_system_capacity(self) -> int:
        """Number of items that can be in the retrieval system simultaneously
        (0 when no retrieval system is configured)."""
        return self._retrieval_system_capacity

    def set_result_delayed_hit_qlen(self, d1: np.ndarray, dfull: np.ndarray) -> None:
        """Set the per-item delayed-hit queue length (called by solver).

        d1 is the mean number of secondary requests waiting on the in-flight
        fetch of each item; dfull additionally counts the request that triggered
        the fetch.
        """
        self._actual_delayed_hit_qlen = np.asarray(d1)
        self._actual_delayed_hit_qlen_full = np.asarray(dfull)

    def get_delayed_hit_qlen(self):
        """Per-item delayed-hit queue length (d1, dfull), or (None, None)."""
        return (getattr(self, '_actual_delayed_hit_qlen', None),
                getattr(self, '_actual_delayed_hit_qlen_full', None))

    def set_result_residt(self, expected_latency: np.ndarray) -> None:
        """Set the computed expected latency (called by solver)."""
        self._actual_residt = np.asarray(expected_latency)

    def get_residt(self) -> Optional[np.ndarray]:
        """Get the expected latency per class, or None if not computed yet."""
        return self._actual_residt

    def reset(self) -> None:
        """Clear solver result fields."""
        self._actual_hit_prob = None
        self._actual_miss_prob = None
        self._actual_delayed_hit_prob = None
        self._actual_hit_prob_list = None
        self._actual_item_prob = None
        self._actual_residt = None
        self._actual_delayed_hit_qlen = None
        self._actual_delayed_hit_qlen_full = None

    def set_retrieval_class(self, input_class: JobClass,
                            output_class: JobClass, item: int) -> None:
        """Set the retrieval class for (item, input_class). item is 0-based."""
        self._retrieval_classes[(item, input_class)] = output_class

    def get_retrieval_class(self, input_class: JobClass,
                            item: int) -> Optional[JobClass]:
        """Get the retrieval class for (item, input_class)."""
        return self._retrieval_classes.get((item, input_class))

    def _add_retrieval_routing_entry(self, from_cls: int, to_cls: int,
                                     src_node: int, dst_node: int, prob: float,
                                     allow_zero: bool = False) -> None:
        """Register a retrieval-class routing edge (0-based class/node indices).
        link() injects these into the routing matrix P. With allow_zero=True an
        explicit prob==0 entry is recorded so it can override (delete) a default
        edge inherited from the read class; negative probs are always dropped and
        the internal broadcast path keeps allow_zero=False (zeros = no edge)."""
        if prob < 0 or (prob == 0 and not allow_zero):
            return
        self._retrieval_routing_entries.append(
            [from_cls, to_cls, src_node, dst_node, prob])

    def set_retrieval_system(self, jobin_class: JobClass, miss_class: JobClass,
                             queues) -> None:
        """
        Initialise the retrieval system through which a request that misses the
        cache is fetched. While the retrieval is pending, repeat requests for the
        same item become delayed hits; after retrieval completes the job switches
        into ``miss_class``.

        Routing and service are NOT passed here; they are taken from the read class:
          - service: the read class's service distribution at each queue. Call
            ``queue.set_service(jobin_class, ...)`` beforehand; override per item with
            ``queue.set_item_service_rate(cache, jobin_class, item, rate)``.
          - routing: the read class's routing among the retrieval queues drawn in the
            top-level routing matrix P; override per item with ``set_item_routing_prob``
            (pass the cache as source for a cache->queue entry, or as dest for a
            queue->cache exit).

        Args:
            jobin_class: arrival JobClass routed through the retrieval system
            miss_class: JobClass a completed retrieval transitions into
            queues: a Queue or list of Queue nodes comprising the system
        """
        import copy as _copy
        from .classes import OpenClass, ClosedClass
        from ..distributions.discrete import DiscreteSampler

        if not isinstance(queues, (list, tuple)):
            queue_arr = [queues]
        else:
            queue_arr = list(queues)
        n_queues = len(queue_arr)
        n_items = self._num_items
        if n_queues == 0:
            raise ValueError("Retrieval system cannot be initialised with no queues")
        self._retrieval_system_capacity = n_items - self._total_cache_capacity

        # inherit the read class's service distribution at each queue as the per-item default
        service_dist_by_queue = []
        for q in queue_arr:
            d = q.get_service(jobin_class)
            if d is None:
                raise ValueError(
                    f"No service distribution for the read class at queue {q.name}; "
                    "call queue.set_service(read_class, ...) before set_retrieval_system.")
            service_dist_by_queue.append(d)

        model = self.get_model()
        self._retrieval_system_queue_indices[jobin_class.get_index0()] = \
            [q.get_index0() for q in queue_arr]

        is_closed = isinstance(jobin_class, ClosedClass)
        ref_station = getattr(jobin_class, '_refstat', None)
        for i in range(n_items):
            if is_closed:
                r_class = ClosedClass(model, f"{jobin_class.name}_retrievalClass_{i + 1}", 0, ref_station, 0)
            else:
                r_class = OpenClass(model, f"{jobin_class.name}_retrievalClass_{i + 1}", 0)
            self.set_retrieval_class(jobin_class, r_class, i)
            self._retrieval_class_indices.add(r_class.get_index0())

            # the retrieval class always reads item i (one-hot popularity); on the
            # READ when it returns to the cache it is logged as a miss.
            item_popularity = np.zeros(n_items)
            item_popularity[i] = 1.0
            self._read_process[r_class] = DiscreteSampler(item_popularity)
            self.set_miss_class(r_class, miss_class)

            # service for the retrieval class at each queue (inherited read-class service);
            # routing is inherited at link() from P and overridable via the setItem* methods.
            for sq in range(n_queues):
                source_queue = queue_arr[sq]
                source_queue.set_service(r_class, _copy.deepcopy(service_dist_by_queue[sq]))
                source_queue.set_class_capacity(r_class, 1)

    def set_item_routing_probability(self, jobin_class: JobClass, item: int,
                                     source, dest, probability: float) -> None:
        """Probability of routing the retrieval class for ``item`` between two nodes
        of the retrieval system. ``source``/``dest`` are either a retrieval queue or
        the cache itself: pass the cache as ``source`` for a cache->queue entry, or as
        ``dest`` for a queue->cache exit. ``item`` is 0-based. Requires a prior
        :meth:`set_retrieval_system` call."""
        r_class = self.get_retrieval_class(jobin_class, item)
        if r_class is None:
            raise ValueError("No retrieval class defined for the given class/item; "
                             "call set_retrieval_system first.")
        rc_idx = r_class.get_index0()
        self._add_retrieval_routing_entry(
            rc_idx, rc_idx, source.get_index0(), dest.get_index0(), probability, allow_zero=True)

    setRetrievalSystem = set_retrieval_system
    getResidT = get_residt
    setResultResidT = set_result_residt
    setItemRoutingProbability = set_item_routing_probability
    # short alias (Probability -> Prob)
    set_item_routing_prob = set_item_routing_probability
    setItemRoutingProb = set_item_routing_probability

    # =========================================================================
    # Aliases (MATLAB compatibility)
    # =========================================================================

    setHitClass = set_hit_class
    getHitClass = get_hit_class
    setMissClass = set_miss_class
    getMissClass = get_miss_class
    setRead = set_read
    getRead = get_read
    setAccessProb = set_access_prob
    getAccessProb = get_access_prob
    setResultHitProb = set_result_hit_prob
    setResultMissProb = set_result_miss_prob
    getHitRatio = get_hit_ratio
    getMissRatio = get_miss_ratio
    hit_ratio = get_hit_ratio
    miss_ratio = get_miss_ratio

    # Property accessors (MATLAB style)
    def getNumItems(self) -> int:
        """Get number of items (MATLAB compatibility)."""
        return self.num_items

    def getNumLevels(self) -> int:
        """Get number of cache levels (MATLAB compatibility)."""
        return self.num_levels

    def getCapacity(self) -> int:
        """Get total cache capacity (MATLAB compatibility)."""
        return self.total_capacity

    def getReplacementStrategy(self) -> ReplacementStrategy:
        """Get replacement strategy (MATLAB compatibility)."""
        return self.replacement_strategy

    def getItemLevelCap(self) -> np.ndarray:
        """Get per-level capacity array (MATLAB compatibility)."""
        return self.item_level_cap


class TimingStrategy(IntEnum):
    """
    Timing strategies for Petri net transitions.
    """
    TIMED = 0       # Timed transition with delay
    IMMEDIATE = 1   # Immediate transition (zero delay)


class Mode:
    """
    A firing mode for a Petri net transition.

    Modes define how a transition can fire, including:
    - Enabling conditions (required tokens from input places)
    - Inhibiting conditions (blocking tokens in places)
    - Firing outcomes (tokens produced to output places)
    - Timing distribution for the firing delay
    """

    def __init__(self, transition: 'Transition', name: str, index: int):
        """
        Initialize a Mode.

        Args:
            transition: Parent Transition object
            name: Mode name
            index: 1-based index of this mode
        """
        self._transition = transition
        self._name = name
        self._index = index  # 1-based index

    @property
    def name(self) -> str:
        """Get mode name."""
        return self._name

    @property
    def index(self) -> int:
        """Get 1-based index (MATLAB compatibility)."""
        return self._index

    def __index__(self) -> int:
        """Allow Mode to be used as an index (0-based for Python arrays)."""
        return self._index - 1

    def __int__(self) -> int:
        """Convert to int (1-based index)."""
        return self._index

    def __repr__(self) -> str:
        return f"Mode('{self._name}', index={self._index})"


class Place(Station):
    """
    Place node for Stochastic Petri Nets.

    A Place represents a location where tokens (jobs) can accumulate.
    Places have infinite capacity (like Delay nodes) and infinite servers.

    In queueing network terms, a Place is similar to a delay station
    but tokens are processed when they move to Transitions.
    """

    def __init__(self, model, name: str, sched_strategy=None):
        """
        Initialize a Place node.

        Args:
            model: Network instance
            name: Place name
            sched_strategy: optional scheduling strategy of the embedded queue.
                When provided (and a service process is later assigned via
                set_service), the place becomes a queueing place (QPN semantics)
                rather than an ordinary/INF pass-through place.
        """
        super().__init__(NodeType.PLACE, name)
        self.set_model(model)
        model.add_node(self)

        # Ordinary places have infinite capacity and servers (like delay nodes)
        self._number_of_servers = np.inf
        self._capacity = np.inf
        self._sched_strategy = sched_strategy if sched_strategy is not None else SchedStrategy.INF
        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        self._sched_strategy_explicit = sched_strategy is not None
        # Per-class capacity is stored in the base Station attribute
        # (_class_capacity), which get_class_capacity / refresh_capacity read.
        self._state = None  # Token state per class
        # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
        self._drop_rule = {}
        # Queueing-place fields (QPN): empty/False for an ordinary place
        self._queueing = False
        self._service_process = {}       # per-class embedded-queue service processes
        self._departure_discipline = {}  # per-class depository departure discipline

    def set_service(self, jobclass: JobClass, distribution) -> None:
        """
        Assign a service process to a token color, turning this ordinary place
        into a queueing place (QPN semantics). The embedded queue serves tokens
        under the place's scheduling strategy; on completion tokens move to the
        depository from which output transitions consume them.

        Args:
            jobclass: token color (job class)
            distribution: service-time distribution of the embedded queue
        """
        from ..constants import DepartureDiscipline
        if not self._queueing:
            # see _kb/04-networkstruct.md (Node/process construction notes) for rationale
            if self._sched_strategy == SchedStrategy.INF and not self._sched_strategy_explicit:
                self._sched_strategy = SchedStrategy.FCFS
            if self._sched_strategy != SchedStrategy.INF and np.isinf(self._number_of_servers):
                self._number_of_servers = 1
            self._queueing = True
        self._service_process[jobclass] = distribution
        self._departure_discipline.setdefault(jobclass, DepartureDiscipline.NORMAL)
        self._invalidate_java()
        if self._model is not None:
            self._model._sn = None

    setService = set_service

    def get_service(self, jobclass: JobClass):
        """Get the embedded-queue service distribution for a job class (or None)."""
        return self._service_process.get(jobclass)

    getService = get_service

    def is_queueing(self) -> bool:
        """True if this place is a queueing place (has an embedded queue)."""
        return self._queueing

    isQueueing = is_queueing

    def set_departure_discipline(self, jobclass: JobClass, discipline) -> None:
        """Set the depository departure discipline for a token color."""
        self._departure_discipline[jobclass] = discipline
        self._invalidate_java()

    setDepartureDiscipline = set_departure_discipline

    def set_number_of_servers(self, k) -> None:
        """Set the number of servers of the embedded queue (multiserver queueing place)."""
        self._number_of_servers = k
        self._invalidate_java()

    setNumberOfServers = set_number_of_servers

    def set_class_capacity(self, jobclass: JobClass, capacity: float) -> None:
        """
        Set per-class capacity limit.

        Args:
            jobclass: Job class
            capacity: Capacity for this class
        """
        self._class_capacity[jobclass] = capacity
        self._invalidate_java()

    setClassCapacity = set_class_capacity

    def set_drop_rule(self, jobclass: JobClass, drop_rule) -> None:
        """
        Set the per-class drop rule for a bounded place.

        Determines what happens to a token that would exceed the place capacity
        (DropStrategy.DROP for a loss place, WAITQ/BAS/... for blocking). Mirrors
        MATLAB Place.setDropRule.

        Args:
            jobclass: token color (job class)
            drop_rule: a DropStrategy value
        """
        self._drop_rule[jobclass] = drop_rule
        self._invalidate_java()

    setDropRule = set_drop_rule

    def set_state(self, state) -> None:
        """
        Set initial token state for this place.

        Args:
            state: Array of initial token counts per class
        """
        self._state = np.asarray(state)
        self._state_explicitly_set = True
        self._invalidate_java()

    setState = set_state
    # Petri-net terminology alias for set_state: the initial token marking.
    set_marking = set_state
    setMarking = set_state

    def get_state(self):
        """Get the current state of this place."""
        return self._state

    getState = get_state


class Transition(StatefulNode):
    """
    Transition node for Stochastic Petri Nets.

    A Transition consumes tokens from input Places and produces tokens
    to output Places according to defined firing modes.

    Each Transition can have multiple firing modes, each with:
    - Enabling conditions: required tokens from input places
    - Inhibiting conditions: blocking conditions from places
    - Firing outcomes: tokens produced to output places
    - Timing distribution: delay before firing
    - Priority and weight for conflict resolution
    """

    def __init__(self, model, name: str):
        """
        Initialize a Transition node.

        Args:
            model: Network instance
            name: Transition name
        """
        super().__init__(NodeType.TRANSITION, name)
        self.set_model(model)
        model.add_node(self)

        self._capacity = np.inf

        # Mode-indexed data structures (lists, 0-based internal indexing)
        self._mode_names = []  # List of mode names
        self._modes = []  # List of Mode objects
        self._enabling_conditions = []  # List of (nnodes x nclasses) matrices
        self._inhibiting_conditions = []  # List of (nnodes x nclasses) matrices
        self._firing_outcomes = []  # List of (nnodes x nclasses) matrices
        self._number_of_servers = []  # List of server counts per mode
        self._timing_strategies = []  # List of TimingStrategy per mode
        self._distributions = []  # List of distributions per mode
        self._firing_priorities = []  # List of priorities per mode
        self._firing_weights = []  # List of weights per mode
        # Marking-dependent firing-rate multiplier g_mode(marking); None == unit.
        self._firing_rate_dependence = []

    def add_mode(self, mode_name: str) -> Mode:
        """
        Add a new firing mode to this transition.

        Args:
            mode_name: Name for the new mode

        Returns:
            Mode object representing the new mode
        """
        nclasses = self._model.get_number_of_classes()
        nnodes = self._model.get_number_of_nodes()

        self._mode_names.append(mode_name)

        # Initialize matrices for this mode
        self._enabling_conditions.append(np.zeros((nnodes, nclasses)))
        self._inhibiting_conditions.append(np.full((nnodes, nclasses), np.inf))
        self._firing_outcomes.append(np.zeros((nnodes, nclasses)))

        # Initialize mode parameters
        self._number_of_servers.append(1)
        self._timing_strategies.append(TimingStrategy.TIMED)
        self._distributions.append(None)  # Will be set by setDistribution
        self._firing_priorities.append(1.0)
        self._firing_weights.append(1.0)
        self._firing_rate_dependence.append(None)

        # Create and store Mode object
        mode = Mode(self, mode_name, len(self._modes) + 1)  # 1-based index
        self._modes.append(mode)

        self._invalidate_java()
        return mode

    addMode = add_mode

    def get_number_of_modes(self) -> int:
        """Get the number of firing modes."""
        return len(self._mode_names)

    getNumberOfModes = get_number_of_modes

    def get_modes(self) -> list:
        """Get list of Mode objects."""
        return self._modes

    getModes = get_modes

    def _get_mode_index(self, mode) -> int:
        """
        Get 0-based index from a mode identifier.

        Args:
            mode: Mode object, 1-based index, or 0-based index

        Returns:
            0-based index for internal arrays
        """
        if isinstance(mode, Mode):
            return mode._index - 1  # Mode._index is 1-based
        elif hasattr(mode, '__index__'):
            # Could be a Mode or int-like object
            idx = mode.__index__()
            # If the index is >= len(modes), assume it's 1-based
            if idx >= len(self._modes):
                return idx - 1  # 1-based to 0-based
            return idx
        else:
            # Assume 1-based integer index (MATLAB convention)
            return int(mode) - 1

    def set_enabling_conditions(self, mode, jobclass: JobClass,
                                 input_node, enabling_condition: int) -> None:
        """
        Set enabling conditions for a mode.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            input_node: Input Place node
            enabling_condition: Number of tokens required
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(input_node, 'get_index0'):
            node_idx = input_node.get_index0()
        elif hasattr(input_node, 'get_index'):
            node_idx = input_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(input_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._enabling_conditions[mode_idx][node_idx, class_idx] = enabling_condition
        self._invalidate_java()

    setEnablingConditions = set_enabling_conditions

    def set_inhibiting_conditions(self, mode, jobclass: JobClass,
                                   input_node, inhibiting_condition: int) -> None:
        """
        Set inhibiting conditions for a mode.

        The transition cannot fire if the place has >= inhibiting_condition tokens.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            input_node: Input Place node
            inhibiting_condition: Number of tokens that inhibit firing
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(input_node, 'get_index0'):
            node_idx = input_node.get_index0()
        elif hasattr(input_node, 'get_index'):
            node_idx = input_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(input_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._inhibiting_conditions[mode_idx][node_idx, class_idx] = inhibiting_condition
        self._invalidate_java()

    setInhibitingConditions = set_inhibiting_conditions

    def set_firing_outcome(self, mode, jobclass: JobClass,
                           output_node, firing_outcome: int) -> None:
        """
        Set firing outcome for a mode.

        Args:
            mode: Mode object or index (1-based)
            jobclass: Job class
            output_node: Output Place or Sink node
            firing_outcome: Number of tokens produced
        """
        mode_idx = self._get_mode_index(mode)

        # Get node index
        if hasattr(output_node, 'get_index0'):
            node_idx = output_node.get_index0()
        elif hasattr(output_node, 'get_index'):
            node_idx = output_node.get_index() - 1
        else:
            node_idx = self._model.get_node_index(output_node.name) - 1

        # Get class index
        if hasattr(jobclass, 'get_index0'):
            class_idx = jobclass.get_index0()
        elif hasattr(jobclass, 'get_index'):
            class_idx = jobclass.get_index() - 1
        else:
            class_idx = int(jobclass)

        self._firing_outcomes[mode_idx][node_idx, class_idx] = firing_outcome
        self._invalidate_java()

    setFiringOutcome = set_firing_outcome

    def set_distribution(self, mode, distribution) -> None:
        """
        Set firing time distribution for a mode.

        Args:
            mode: Mode object or index (1-based)
            distribution: Timing distribution (Exp, Erlang, etc.)
        """
        mode_idx = self._get_mode_index(mode)
        self._distributions[mode_idx] = distribution
        self._invalidate_java()

    setDistribution = set_distribution

    def get_distribution(self, mode):
        """Get firing time distribution for a mode."""
        mode_idx = self._get_mode_index(mode)
        return self._distributions[mode_idx]

    getDistribution = get_distribution

    def set_firing_rate_dependence(self, mode, g) -> None:
        """
        Marking-dependent firing-rate multiplier for a timed mode.

        g(m) takes the node-indexed input-place marking matrix and returns a
        positive scalar; the effective firing rate of an enabled binding becomes
        rate_base(mode)*g(m). Exact only for memoryless firing, so the mode must
        be TIMED and exponentially distributed (mirrors the PS/FCFS-only station
        load/class dependence). Pass None to restore the unit multiplier.
        """
        mode_idx = self._get_mode_index(mode)
        if g is not None:
            if self._timing_strategies[mode_idx] == TimingStrategy.IMMEDIATE:
                raise ValueError("Firing-rate dependence is not supported for IMMEDIATE modes; "
                                 "use set_firing_weights for immediate transitions.")
            dist = self._distributions[mode_idx]
            if dist is None or type(dist).__name__ != 'Exp':
                raise ValueError("Firing-rate dependence is supported only for exponentially "
                                 "distributed (memoryless) firing; set an Exp distribution first.")
        self._firing_rate_dependence[mode_idx] = g
        self._invalidate_java()

    setFiringRateDependence = set_firing_rate_dependence

    def get_firing_rate_dependence(self, mode):
        """Get the marking-dependent firing-rate multiplier for a mode (or None)."""
        mode_idx = self._get_mode_index(mode)
        return self._firing_rate_dependence[mode_idx]

    getFiringRateDependence = get_firing_rate_dependence

    def set_number_of_servers(self, mode, num_servers: int) -> None:
        """
        Set number of servers for a mode.

        Args:
            mode: Mode object or index (1-based)
            num_servers: Number of servers (or GlobalConstants.MaxInt for infinite)
        """
        mode_idx = self._get_mode_index(mode)
        self._number_of_servers[mode_idx] = num_servers
        self._invalidate_java()

    setNumberOfServers = set_number_of_servers

    def set_timing_strategy(self, mode, timing_strategy: TimingStrategy) -> None:
        """
        Set timing strategy for a mode.

        Args:
            mode: Mode object or index (1-based)
            timing_strategy: TIMED or IMMEDIATE
        """
        mode_idx = self._get_mode_index(mode)
        self._timing_strategies[mode_idx] = timing_strategy
        self._invalidate_java()

    setTimingStrategy = set_timing_strategy

    def set_firing_priorities(self, mode, priority: float) -> None:
        """
        Set firing priority for a mode.

        Args:
            mode: Mode object or index (1-based)
            priority: Priority value (higher = more priority)
        """
        mode_idx = self._get_mode_index(mode)
        self._firing_priorities[mode_idx] = priority
        self._invalidate_java()

    setFiringPriorities = set_firing_priorities

    def set_firing_weights(self, mode, weight: float) -> None:
        """
        Set firing weight for a mode.

        Args:
            mode: Mode object or index (1-based)
            weight: Weight for probabilistic selection among enabled modes
        """
        mode_idx = self._get_mode_index(mode)
        self._firing_weights[mode_idx] = weight
        self._invalidate_java()

    setFiringWeights = set_firing_weights

    def set_mode_names(self, mode, mode_name: str) -> None:
        """
        Set name for a mode.

        Args:
            mode: Mode object or index (1-based)
            mode_name: New name for the mode
        """
        mode_idx = self._get_mode_index(mode)
        self._mode_names[mode_idx] = mode_name
        self._invalidate_java()

    setModeNames = set_mode_names

    def get_service_rates(self):
        """
        Get service rates for all modes.

        Returns:
            Tuple of (map, mu, phi) lists for PH-type service rates
        """
        nmodes = self.get_number_of_modes()
        map_list = [None] * nmodes
        mu_list = [None] * nmodes
        phi_list = [None] * nmodes

        for m in range(nmodes):
            dist = self._distributions[m]
            if dist is None:
                map_list[m] = [[np.nan], [np.nan]]
                mu_list[m] = np.nan
                phi_list[m] = np.nan
            elif hasattr(dist, 'getProcess'):
                map_list[m] = dist.getProcess()
                mu_list[m] = dist.getMu() if hasattr(dist, 'getMu') else None
                phi_list[m] = dist.getPhi() if hasattr(dist, 'getPhi') else None
            elif hasattr(dist, 'getMean'):
                rate = 1.0 / dist.getMean()
                map_list[m] = [[-rate], [rate]]
                mu_list[m] = [rate]
                phi_list[m] = [1]
            else:
                map_list[m] = [[np.nan], [np.nan]]
                mu_list[m] = np.nan
                phi_list[m] = np.nan

        return map_list, mu_list, phi_list

    getServiceRates = get_service_rates


class Logger(Node):
    """
    Logger node for recording job passage.

    A Logger node records arrival and departure timestamps for jobs
    passing through it. Used internally by getCdfRespT to collect
    response time samples via transient simulation.

    Ported from MATLAB's Logger class in matlab/src/lang/nodes/Logger.m
    """

    def __init__(self, model, name: str, log_file_name: str):
        """
        Initialize a Logger node.

        Args:
            model: Network instance
            name: Logger name
            log_file_name: Full path to the log file
        """
        import os

        super().__init__(NodeType.LOGGER, name)
        self.set_model(model)
        model.add_node(self)

        # Parse file name and path
        dirname, basename = os.path.split(log_file_name)
        self._file_name = basename
        self._file_path = dirname if dirname else model.get_log_path()

        # Logging options (matching MATLAB defaults)
        self._want_start_time = False
        self._want_logger_name = False
        self._want_timestamp = True
        self._want_job_id = True
        self._want_job_class = True
        self._want_time_same_class = False
        self._want_time_any_class = False

        # Scheduling (Logger uses FCFS, non-preemptive)
        self._sched_policy = SchedStrategyType.NP
        self._sched_strategy = SchedStrategy.FCFS
        self._capacity = float('inf')

    @property
    def file_name(self) -> str:
        """Get the log file name."""
        return self._file_name

    @property
    def file_path(self) -> str:
        """Get the log file path."""
        return self._file_path

    # Getters for logging options
    def get_start_time(self) -> bool:
        return self._want_start_time

    def get_logger_name(self) -> bool:
        return self._want_logger_name

    def get_timestamp(self) -> bool:
        return self._want_timestamp

    def get_job_id(self) -> bool:
        return self._want_job_id

    def get_job_class(self) -> bool:
        return self._want_job_class

    def get_time_same_class(self) -> bool:
        return self._want_time_same_class

    def get_time_any_class(self) -> bool:
        return self._want_time_any_class

    # Setters for logging options
    def set_start_time(self, value: bool) -> None:
        self._want_start_time = value

    def set_logger_name(self, value: bool) -> None:
        self._want_logger_name = value

    def set_timestamp(self, value: bool) -> None:
        self._want_timestamp = value

    def set_job_id(self, value: bool) -> None:
        self._want_job_id = value

    def set_job_class(self, value: bool) -> None:
        self._want_job_class = value

    def set_time_same_class(self, value: bool) -> None:
        self._want_time_same_class = value

    def set_time_any_class(self, value: bool) -> None:
        self._want_time_any_class = value

    def set_prob_routing(self, jobclass: JobClass, destination, probability: float) -> None:
        """
        Set probabilistic routing to a destination.

        Args:
            jobclass: Job class
            destination: Destination node
            probability: Routing probability
        """
        self.set_routing(jobclass, RoutingStrategy.PROB, destination, probability)

    # CamelCase aliases
    getStartTime = get_start_time
    getLoggerName = get_logger_name
    getTimestamp = get_timestamp
    getJobID = get_job_id
    getJobClass = get_job_class
    getTimeSameClass = get_time_same_class
    getTimeAnyClass = get_time_any_class
    setStartTime = set_start_time
    setLoggerName = set_logger_name
    setTimestamp = set_timestamp
    setJobID = set_job_id
    setJobClass = set_job_class
    setTimeSameClass = set_time_same_class
    setTimeAnyClass = set_time_any_class
    setProbRouting = set_prob_routing
    fileName = property(lambda self: self._file_name)
    filePath = property(lambda self: self._file_path)


__all__ = [
    'Queue', 'Delay', 'Source', 'Sink', 'Fork', 'Join',
    'Router', 'ClassSwitch', 'Cache', 'ReplacementStrategy',
    'Place', 'Transition', 'Mode', 'TimingStrategy', 'Logger'
]
