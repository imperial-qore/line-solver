"""
Job-based scheduling strategies.

This module implements scheduling based on job characteristics:
- SJF: Shortest Job First (service time based)
- LJF: Longest Job First (service time based)
- SEPT: Shortest Expected Processing Time (class-based)
- LEPT: Longest Expected Processing Time (class-based)
"""

from typing import List, Optional, Callable, Tuple, Any, Dict
import heapq

from .base import Customer, SchedulingStrategy


class SJFScheduler(SchedulingStrategy):
    """
    Shortest Job First (SJF) scheduling.

    Jobs with shorter service times are served first.
    Service time is pre-sampled when the job arrives.
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (service_time, arrival_time, customer)
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

        # Pre-sample service time for SJF ordering
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        server_id = self._find_free_server()

        if server_id >= 0:
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # SJF: shorter service time first, FCFS tie-break
            heapq.heappush(
                self.wait_queue,
                (customer.service_time, customer.queue_arrival_time, customer)
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


class LJFScheduler(SchedulingStrategy):
    """
    Longest Job First (LJF) scheduling.

    Jobs with longer service times are served first.
    Service time is pre-sampled when the job arrives.
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)
        # Priority queue: (negative_service_time, arrival_time, customer)
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

        # Pre-sample service time for LJF ordering
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        server_id = self._find_free_server()

        if server_id >= 0:
            self.server_busy[server_id] = True
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer))
        else:
            # LJF: longer service time first (negative for max-heap)
            heapq.heappush(
                self.wait_queue,
                (-customer.service_time, customer.queue_arrival_time, customer)
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


class SEPTScheduler(SchedulingStrategy):
    """
    Shortest Expected Processing Time (SEPT) scheduling.

    Jobs from classes with higher service rates (faster service)
    are served first. This is class-based, not job-based.
    """

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        service_rates: Optional[List[float]] = None,
    ):
        super().__init__(num_classes, num_servers)
        # Expected processing time = 1/mu (lower is better)
        self.service_rates = service_rates or [1.0] * num_classes
        self.expected_times = [
            1.0 / rate if rate > 0 else float('inf')
            for rate in self.service_rates
        ]

        # Priority queue: (expected_time, arrival_time, customer)
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
            # SEPT: shorter expected time first
            expected_time = self.expected_times[class_id]
            heapq.heappush(
                self.wait_queue,
                (expected_time, customer.queue_arrival_time, customer)
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


class LEPTScheduler(SchedulingStrategy):
    """
    Longest Expected Processing Time (LEPT) scheduling.

    Jobs from classes with lower service rates (slower service)
    are served first. This is class-based, not job-based.
    """

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        service_rates: Optional[List[float]] = None,
    ):
        super().__init__(num_classes, num_servers)
        # Expected processing time = 1/mu (higher is better for LEPT)
        self.service_rates = service_rates or [1.0] * num_classes
        self.expected_times = [
            1.0 / rate if rate > 0 else float('inf')
            for rate in self.service_rates
        ]

        # Priority queue: (negative_expected_time, arrival_time, customer)
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
            # LEPT: longer expected time first (negative for max-heap)
            expected_time = self.expected_times[class_id]
            heapq.heappush(
                self.wait_queue,
                (-expected_time, customer.queue_arrival_time, customer)
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
