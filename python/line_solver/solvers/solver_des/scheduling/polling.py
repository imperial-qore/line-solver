"""
Polling scheduling strategies.

This module implements polling-based scheduling disciplines:
- Exhaustive: Serve all jobs in a class before switching
- Gated: Serve only jobs present at visit start
- K-Limited: Serve at most K jobs per visit

Polling systems visit each class in a cyclic order with
optional switchover times between classes.
"""

from typing import List, Optional, Callable, Tuple, Any, Deque
from collections import deque
from enum import Enum
from dataclasses import dataclass, field

from .base import Customer, SchedulingStrategy


class PollingPolicy(Enum):
    """Polling service policy."""
    EXHAUSTIVE = "exhaustive"
    GATED = "gated"
    KLIMITED = "k_limited"


@dataclass
class PollingConfig:
    """
    Configuration for polling scheduler.

    Attributes:
        policy: Polling policy (EXHAUSTIVE, GATED, K-LIMITED)
        k_limit: Maximum jobs per visit (for K-LIMITED)
        switchover_times: Switchover time from class i to next class
        class_order: Order in which to visit classes (default: 0, 1, 2, ...)
    """
    policy: PollingPolicy = PollingPolicy.EXHAUSTIVE
    k_limit: int = 1
    switchover_times: Optional[List[float]] = None
    class_order: Optional[List[int]] = None


class PollingScheduler(SchedulingStrategy):
    """
    Polling scheduler with cyclic class visits.

    The server visits each job class in order, serving according to
    the polling policy:
    - EXHAUSTIVE: Serve all jobs in queue before switching
    - GATED: Serve only jobs present when visit started
    - K-LIMITED: Serve at most K jobs, then switch

    Switchover times can be specified between class visits.
    """

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        policy: PollingPolicy = PollingPolicy.EXHAUSTIVE,
        k_limit: int = 1,
        switchover_times: Optional[List[float]] = None,
        class_order: Optional[List[int]] = None,
        **kwargs
    ):
        super().__init__(num_classes, num_servers)

        self.policy = policy
        self.k_limit = k_limit

        # Switchover times (from class i to next)
        if switchover_times is None:
            self.switchover_times = [0.0] * num_classes
        else:
            self.switchover_times = list(switchover_times)

        # Class visiting order
        if class_order is None:
            self.class_order = list(range(num_classes))
        else:
            self.class_order = list(class_order)

        # Per-class queues
        self.class_queues: List[Deque[Customer]] = [
            deque() for _ in range(num_classes)
        ]

        # Gated queue: jobs present at visit start
        self.gated_queue: Deque[Customer] = deque()

        # Current state
        self.current_class_idx = 0  # Index into class_order
        self.current_class = self.class_order[0]
        self.jobs_served_this_visit = 0
        self.visit_started = False
        self.in_switchover = False
        self.switchover_end_time = 0.0

        # Currently in service
        self.in_service: Optional[Customer] = None
        self.server_busy = False

        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        """
        Handle arrival - add to appropriate class queue.
        """
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        # Add to class queue
        self.class_queues[class_id].append(customer)

        # If server is idle and not in switchover, try to serve
        if not self.server_busy and not self.in_switchover:
            next_customer = self._get_next_customer()
            if next_customer is not None:
                self.server_busy = True
                self.in_service = next_customer
                self._busy_servers[next_customer.class_id] += 1
                return (True, (0, next_customer))

        return (True, None)

    def _get_next_customer(self) -> Optional[Customer]:
        """
        Get next customer to serve based on polling policy.
        """
        # First, try to get from current class
        customer = self._get_from_current_class()
        if customer is not None:
            return customer

        # Need to switch to next class
        # Check all classes in order
        for _ in range(len(self.class_order)):
            self._advance_to_next_class()
            customer = self._get_from_current_class()
            if customer is not None:
                return customer

        return None

    def _get_from_current_class(self) -> Optional[Customer]:
        """
        Get customer from current class if available and allowed.
        """
        # Check K-limit
        if self.policy == PollingPolicy.KLIMITED:
            if self.jobs_served_this_visit >= self.k_limit:
                return None

        # Get from appropriate queue
        if self.policy == PollingPolicy.GATED:
            if self.gated_queue:
                customer = self.gated_queue.popleft()
                self.jobs_served_this_visit += 1
                return customer
        else:
            # EXHAUSTIVE or K-LIMITED
            queue = self.class_queues[self.current_class]
            if queue:
                customer = queue.popleft()
                self.jobs_served_this_visit += 1
                return customer

        return None

    def _advance_to_next_class(self) -> None:
        """
        Advance to next class in polling order.
        """
        self.current_class_idx = (self.current_class_idx + 1) % len(self.class_order)
        self.current_class = self.class_order[self.current_class_idx]
        self.jobs_served_this_visit = 0

        # For GATED policy, snapshot the queue
        if self.policy == PollingPolicy.GATED:
            self.gated_queue = deque(self.class_queues[self.current_class])
            self.class_queues[self.current_class].clear()

    def get_switchover_time(self) -> float:
        """Get switchover time from current class to next."""
        return self.switchover_times[self.current_class_idx]

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        """
        Handle departure - get next customer or switch class.
        """
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.server_busy = False
        self.in_service = None

        # Try to get next customer
        next_customer = self._get_next_customer()

        if next_customer is not None:
            self.server_busy = True
            self.in_service = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def should_switch_class(self) -> bool:
        """
        Check if we should switch to next class.

        Returns True if:
        - EXHAUSTIVE: Current class queue is empty
        - GATED: Gated queue is empty
        - K-LIMITED: Served K jobs or queue empty
        """
        if self.policy == PollingPolicy.EXHAUSTIVE:
            return len(self.class_queues[self.current_class]) == 0

        elif self.policy == PollingPolicy.GATED:
            return len(self.gated_queue) == 0

        elif self.policy == PollingPolicy.KLIMITED:
            return (self.jobs_served_this_visit >= self.k_limit or
                    len(self.class_queues[self.current_class]) == 0)

        return False

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def get_current_class(self) -> int:
        """Get the class currently being served."""
        return self.current_class

    def is_queue_empty(self) -> bool:
        """Check if all queues are empty."""
        return all(len(q) == 0 for q in self.class_queues)


class ExhaustivePollingScheduler(PollingScheduler):
    """Convenience class for exhaustive polling."""

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        switchover_times: Optional[List[float]] = None,
        **kwargs
    ):
        super().__init__(
            num_classes,
            num_servers,
            policy=PollingPolicy.EXHAUSTIVE,
            switchover_times=switchover_times,
            **kwargs
        )


class GatedPollingScheduler(PollingScheduler):
    """Convenience class for gated polling."""

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        switchover_times: Optional[List[float]] = None,
        **kwargs
    ):
        super().__init__(
            num_classes,
            num_servers,
            policy=PollingPolicy.GATED,
            switchover_times=switchover_times,
            **kwargs
        )


class KLimitedPollingScheduler(PollingScheduler):
    """Convenience class for K-limited polling."""

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        k_limit: int = 1,
        switchover_times: Optional[List[float]] = None,
        **kwargs
    ):
        super().__init__(
            num_classes,
            num_servers,
            policy=PollingPolicy.KLIMITED,
            k_limit=k_limit,
            switchover_times=switchover_times,
            **kwargs
        )


class TwoQueuePollingScheduler(SchedulingStrategy):
    """
    Two-queue polling for asymmetric priority.

    Serves queue 1 exhaustively, then serves one from queue 2,
    then back to queue 1 if non-empty.
    """

    def __init__(self, num_servers: int, high_priority_class: int = 0):
        super().__init__(2, num_servers)
        self.high_priority_class = high_priority_class
        self.low_priority_class = 1 - high_priority_class

        self.queues: List[Deque[Customer]] = [deque(), deque()]
        self.in_service: Optional[Customer] = None
        self.server_busy = False

        self._queue_length = [0, 0]
        self._busy_servers = [0, 0]

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self.queues[class_id].append(customer)

        if not self.server_busy:
            next_customer = self._get_next_customer()
            if next_customer is not None:
                self.server_busy = True
                self.in_service = next_customer
                self._busy_servers[next_customer.class_id] += 1
                return (True, (0, next_customer))

        return (True, None)

    def _get_next_customer(self) -> Optional[Customer]:
        """High priority first, then low priority."""
        if self.queues[self.high_priority_class]:
            return self.queues[self.high_priority_class].popleft()
        if self.queues[self.low_priority_class]:
            return self.queues[self.low_priority_class].popleft()
        return None

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.server_busy = False
        self.in_service = None

        next_customer = self._get_next_customer()
        if next_customer is not None:
            self.server_busy = True
            self.in_service = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]
