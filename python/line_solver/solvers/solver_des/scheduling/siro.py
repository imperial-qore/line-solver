"""
SIRO (Service In Random Order) scheduling strategy.

Jobs are selected for service in random order from the wait queue.
"""

from typing import List, Optional, Callable, Tuple, Any
import heapq

from .base import Customer, SchedulingStrategy


class SIROScheduler(SchedulingStrategy):
    """
    Service In Random Order (SIRO) scheduling.

    Jobs are selected randomly from the wait queue. Each job has a
    random_rank assigned on arrival, and jobs are served in order
    of increasing rank (effectively random selection).
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (random_rank, customer)
        self.wait_queue: List[Tuple[float, Customer]] = []
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
            # SIRO: use random_rank for ordering
            heapq.heappush(self.wait_queue, (customer.random_rank, customer))
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
            _, next_customer = heapq.heappop(self.wait_queue)
            self.server_busy[server_id] = True
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]
