"""
FCFS and priority-based FCFS scheduling strategies.

This module implements:
- FCFS: First-Come-First-Served
- HOL/FCFSPRIO: Head-of-Line priority with FCFS tie-breaking
- EDD: Earliest Due Date (deadline-based)
"""

from typing import List, Optional, Callable, Tuple, Any
import heapq

from .base import Customer, SchedulingStrategy


class FCFSScheduler(SchedulingStrategy):
    """
    First-Come-First-Served (FCFS) scheduling.

    Jobs are served in the order they arrive. This is the simplest
    and most common scheduling discipline.
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
        """Handle customer arrival."""
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        # Find available server
        server_id = self._find_free_server()

        if server_id >= 0:
            # Start service immediately
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # Add to wait queue (FCFS: append to end)
            self.wait_queue.append(customer)
            return (True, None)

    def _find_free_server(self) -> int:
        """Find an available server."""
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
        """Handle departure and get next customer."""
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.server_busy[server_id] = False
        self.in_service[server_id] = None

        # Get next customer from wait queue
        if self.wait_queue:
            next_customer = self.wait_queue.pop(0)
            self.server_busy[server_id] = True
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]


class PriorityFCFSScheduler(SchedulingStrategy):
    """
    Priority scheduling with FCFS tie-breaking (HOL/FCFSPRIO).

    Higher priority jobs are served first. Within the same priority,
    jobs are served in FCFS order.
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (negative_priority, arrival_time, customer)
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
        """Handle customer arrival."""
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        server_id = self._find_free_server()

        if server_id >= 0:
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # Add to priority queue (negative priority for max-heap behavior)
            heapq.heappush(
                self.wait_queue,
                (-customer.priority, customer.queue_arrival_time, customer)
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


class EDDScheduler(SchedulingStrategy):
    """
    Earliest Due Date (EDD) scheduling.

    Jobs with earlier deadlines are served first (non-preemptive).
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (deadline, arrival_time, customer)
        self.wait_queue: List[Tuple[float, float, Customer]] = []
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
            heapq.heappush(
                self.wait_queue,
                (customer.absolute_deadline, customer.queue_arrival_time, customer)
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
