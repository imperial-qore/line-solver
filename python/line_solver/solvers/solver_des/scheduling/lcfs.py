"""
LCFS (Last-Come-First-Served) scheduling strategies.

This module implements non-preemptive LCFS variants:
- LCFS: Last-Come-First-Served (stack-based)
- LCFSPRIO: LCFS with priority ordering
"""

from typing import List, Optional, Callable, Tuple, Any
import heapq

from .base import Customer, SchedulingStrategy


class LCFSScheduler(SchedulingStrategy):
    """
    Last-Come-First-Served (LCFS) scheduling (non-preemptive).

    Jobs are served in reverse arrival order (stack/LIFO).
    Once a job starts service, it is not preempted.
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        self.wait_queue: List[Customer] = []
        self.server_busy: List[bool] = [False] * num_servers
        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        server_id = self._find_free_server()

        if server_id >= 0:
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # LCFS: append to end (will pop from end)
            self.wait_queue.append(customer)
            return (True, None)

    def _find_free_server(self) -> int:
        for i, busy in enumerate(self.server_busy):
            if not busy:
                return i
        return -1

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.server_busy[server_id] = False
        self.in_service[server_id] = None

        if self.wait_queue:
            # LCFS: pop from end (most recent arrival)
            next_customer = self.wait_queue.pop()
            self.server_busy[server_id] = True
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]


class LCFSPriorityScheduler(SchedulingStrategy):
    """
    LCFS with priority ordering (LCFSPRIO).

    Higher priority jobs are served first. Within the same priority,
    jobs are served in LCFS order (most recent first).
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (negative_priority, negative_arrival_time, customer)
        self.wait_queue: List[Tuple[int, float, Customer]] = []
        self.server_busy: List[bool] = [False] * num_servers
        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        server_id = self._find_free_server()

        if server_id >= 0:
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # LCFSPRIO: negative priority, negative arrival time for LCFS within priority
            heapq.heappush(
                self.wait_queue,
                (-customer.priority, -customer.queue_arrival_time, customer)
            )
            return (True, None)

    def _find_free_server(self) -> int:
        for i, busy in enumerate(self.server_busy):
            if not busy:
                return i
        return -1

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.server_busy[server_id] = False
        self.in_service[server_id] = None

        if self.wait_queue:
            _, _, next_customer = heapq.heappop(self.wait_queue)
            self.server_busy[server_id] = True
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]
