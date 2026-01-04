"""
Preemptive scheduling strategies.

This module implements preemptive scheduling disciplines:
- LCFSPR: Last-Come-First-Served Preemptive Resume
- LCFSPI: Last-Come-First-Served Preemptive Independent (restart)
- SRPT: Shortest Remaining Processing Time
- EDF: Earliest Deadline First (preemptive)
"""

from typing import List, Optional, Callable, Tuple, Any, Dict
import heapq

from .base import Customer, PreemptiveCustomer, PreemptionRecord, SchedulingStrategy


class LCFSPRScheduler(SchedulingStrategy):
    """
    Last-Come-First-Served Preemptive Resume (LCFSPR) scheduling.

    New arrivals preempt the current job. When the preempting job
    completes, the preempted job resumes from where it left off.

    If has_priority is True, only higher priority arrivals can preempt.
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # Stack of preempted jobs per server: server_id -> list of (preemption_record, customer)
        self.preemption_stacks: List[List[Tuple[PreemptionRecord, Customer]]] = [
            [] for _ in range(num_servers)
        ]

        # Currently active jobs per server
        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self.service_start_times: List[float] = [0.0] * num_servers

        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        """
        Handle arrival - may preempt current job.

        Returns:
            (accepted, service_info) where service_info is:
            - (server_id, customer, preempted_info) if service starts
            - None if no server available (shouldn't happen for LCFSPR)
        """
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        # Pre-sample service time if not already set
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # Find a server to use (preempt if needed)
        server_id, preempted = self._find_server_with_preemption(customer, current_time)

        if server_id >= 0:
            self.in_service[server_id] = customer
            self.service_start_times[server_id] = current_time
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer, preempted))

        # This shouldn't happen in pure LCFSPR - always preempts
        return (True, None)

    def _find_server_with_preemption(
        self,
        customer: Customer,
        current_time: float
    ) -> Tuple[int, Optional[Tuple[Customer, PreemptionRecord]]]:
        """
        Find a server, preempting if necessary.

        Returns:
            (server_id, preempted_info) where preempted_info is
            (preempted_customer, preemption_record) or None
        """
        # First, try to find a free server
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return (i, None)

        # No free server - must preempt
        # Find lowest priority job (or any job if no priorities)
        preempt_server = -1
        preempt_priority = float('inf')

        for i in range(self.num_servers):
            current_job = self.in_service[i]
            if current_job is not None:
                if self.has_priority:
                    # Only preempt lower priority
                    if current_job.priority < customer.priority:
                        if current_job.priority < preempt_priority:
                            preempt_priority = current_job.priority
                            preempt_server = i
                else:
                    # LCFSPR without priority - always preempt
                    preempt_server = 0
                    break

        if preempt_server >= 0:
            preempted_job = self.in_service[preempt_server]
            if preempted_job is not None:
                # Calculate remaining work
                elapsed = current_time - self.service_start_times[preempt_server]
                remaining = max(0.0, preempted_job.service_time - elapsed)

                record = PreemptionRecord(
                    customer_id=preempted_job.customer_id,
                    remaining_work=remaining,
                    preemption_time=current_time,
                    elapsed_time=elapsed
                )

                # Push onto preemption stack
                self.preemption_stacks[preempt_server].append((record, preempted_job))

                # Update busy server count
                self._busy_servers[preempted_job.class_id] -= 1

                return (preempt_server, (preempted_job, record))

        return (-1, None)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        """
        Handle departure - resume preempted job if any.

        Returns:
            Next customer to serve (may be a resumed job)
        """
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.in_service[server_id] = None

        # Check preemption stack for resumed job
        if self.preemption_stacks[server_id]:
            record, resumed_customer = self.preemption_stacks[server_id].pop()

            # Update service time to remaining work
            resumed_customer.service_time = record.remaining_work

            self.in_service[server_id] = resumed_customer
            self.service_start_times[server_id] = current_time
            self._busy_servers[resumed_customer.class_id] += 1

            return resumed_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """LCFSPR is a preemptive scheduler."""
        return True


class LCFSPIScheduler(SchedulingStrategy):
    """
    Last-Come-First-Served Preemptive Independent (LCFSPI) scheduling.

    Similar to LCFSPR, but when a preempted job resumes, it restarts
    with a fresh service time (independent of previous service).
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # Stack of preempted customers per server
        self.preemption_stacks: List[List[Customer]] = [
            [] for _ in range(num_servers)
        ]

        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

        # Store service generator for restart
        self._service_gen: Optional[Callable[[int], float]] = None

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self._service_gen = service_gen

        # Pre-sample service time
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        server_id, preempted = self._find_server_with_preemption(customer)

        if server_id >= 0:
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer, preempted))

        return (True, None)

    def _find_server_with_preemption(
        self,
        customer: Customer
    ) -> Tuple[int, Optional[Customer]]:
        """Find server, preempting if needed."""
        # Try free server first
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return (i, None)

        # Must preempt
        preempt_server = -1
        for i in range(self.num_servers):
            current_job = self.in_service[i]
            if current_job is not None:
                if self.has_priority:
                    if current_job.priority < customer.priority:
                        preempt_server = i
                        break
                else:
                    preempt_server = 0
                    break

        if preempt_server >= 0:
            preempted = self.in_service[preempt_server]
            if preempted is not None:
                # Push onto stack (service time will be resampled on resume)
                self.preemption_stacks[preempt_server].append(preempted)
                self._busy_servers[preempted.class_id] -= 1
                return (preempt_server, preempted)

        return (-1, None)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None

        # Resume preempted job with NEW service time (independent restart)
        if self.preemption_stacks[server_id]:
            resumed = self.preemption_stacks[server_id].pop()

            # Resample service time for independent restart
            if self._service_gen is not None:
                resumed.service_time = self._service_gen(resumed.class_id)

            self.in_service[server_id] = resumed
            self._busy_servers[resumed.class_id] += 1
            return resumed

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """LCFSPI is a preemptive scheduler."""
        return True


class SRPTScheduler(SchedulingStrategy):
    """
    Shortest Remaining Processing Time (SRPT) scheduling.

    Preemptive version of SJF. The job with the shortest remaining
    service time gets the server. New arrivals can preempt if they
    have shorter remaining time.

    If has_priority is True, priority is considered first, then SRPT
    within the same priority level.
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # Priority queue: (priority, remaining_time, arrival_time, customer)
        # Lower values = higher priority
        self.wait_queue: List[Tuple[int, float, float, Customer]] = []

        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self.service_start_times: List[float] = [0.0] * num_servers
        self.remaining_times: List[float] = [0.0] * num_servers

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

        # Pre-sample service time
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # Update remaining times for all in-service jobs
        self._update_remaining_times(current_time)

        # Check if we should preempt
        server_id = self._find_preemptable_server(customer)

        if server_id >= 0:
            preempted = self._preempt_and_serve(server_id, customer, current_time)
            return (True, (server_id, customer, preempted))
        else:
            # Add to wait queue
            priority = -customer.priority if self.has_priority else 0
            heapq.heappush(
                self.wait_queue,
                (priority, customer.service_time, customer.queue_arrival_time, customer)
            )
            return (True, None)

    def _update_remaining_times(self, current_time: float) -> None:
        """Update remaining times for all in-service jobs."""
        for i in range(self.num_servers):
            if self.in_service[i] is not None:
                elapsed = current_time - self.service_start_times[i]
                self.remaining_times[i] = max(0.0, self.remaining_times[i] - elapsed)
                self.service_start_times[i] = current_time

    def _find_preemptable_server(self, customer: Customer) -> int:
        """
        Find a server to use, either free or preemptable.

        Returns server_id or -1 if should queue.
        """
        new_priority = -customer.priority if self.has_priority else 0
        new_remaining = customer.service_time

        # First check for free servers
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return i

        # Check for preemption
        worst_server = -1
        worst_remaining = new_remaining
        worst_priority = new_priority

        for i in range(self.num_servers):
            current = self.in_service[i]
            if current is not None:
                cur_priority = -current.priority if self.has_priority else 0
                cur_remaining = self.remaining_times[i]

                # Can preempt if: lower priority, or same priority and longer remaining
                if (cur_priority > worst_priority or
                    (cur_priority == worst_priority and cur_remaining > worst_remaining)):
                    worst_priority = cur_priority
                    worst_remaining = cur_remaining
                    worst_server = i

        # Preempt if we're better than worst
        if self.has_priority:
            if new_priority < worst_priority:
                return worst_server
            elif new_priority == worst_priority and new_remaining < worst_remaining:
                return worst_server
        else:
            if new_remaining < worst_remaining:
                return worst_server

        return -1

    def _preempt_and_serve(
        self,
        server_id: int,
        customer: Customer,
        current_time: float
    ) -> Optional[Customer]:
        """Preempt current job (if any) and start serving new customer."""
        preempted = self.in_service[server_id]

        if preempted is not None:
            # Save remaining time and add to queue
            preempted.service_time = self.remaining_times[server_id]
            priority = -preempted.priority if self.has_priority else 0
            heapq.heappush(
                self.wait_queue,
                (priority, preempted.service_time, preempted.queue_arrival_time, preempted)
            )
            self._busy_servers[preempted.class_id] -= 1

        # Start new customer
        self.in_service[server_id] = customer
        self.service_start_times[server_id] = current_time
        self.remaining_times[server_id] = customer.service_time
        self._busy_servers[customer.class_id] += 1

        return preempted

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None
        self.remaining_times[server_id] = 0.0

        # Update remaining times before selecting next
        self._update_remaining_times(current_time)

        if self.wait_queue:
            _, remaining, _, next_customer = heapq.heappop(self.wait_queue)
            self.in_service[server_id] = next_customer
            self.service_start_times[server_id] = current_time
            self.remaining_times[server_id] = remaining
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """SRPT is a preemptive scheduler."""
        return True


class EDFScheduler(SchedulingStrategy):
    """
    Earliest Deadline First (EDF) preemptive scheduling.

    Jobs with earlier absolute deadlines get priority. A job can be
    preempted if a new arrival has an earlier deadline.
    """

    def __init__(self, num_classes: int, num_servers: int):
        super().__init__(num_classes, num_servers)

        # Priority queue: (deadline, arrival_time, customer)
        self.wait_queue: List[Tuple[float, float, Customer]] = []

        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self.service_start_times: List[float] = [0.0] * num_servers
        self.remaining_times: List[float] = [0.0] * num_servers
        self.deadlines: List[float] = [float('inf')] * num_servers

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

        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # Update remaining times
        self._update_remaining_times(current_time)

        # Check for preemption
        server_id = self._find_preemptable_server(customer)

        if server_id >= 0:
            preempted = self._preempt_and_serve(server_id, customer, current_time)
            return (True, (server_id, customer, preempted))
        else:
            heapq.heappush(
                self.wait_queue,
                (customer.absolute_deadline, customer.queue_arrival_time, customer)
            )
            return (True, None)

    def _update_remaining_times(self, current_time: float) -> None:
        for i in range(self.num_servers):
            if self.in_service[i] is not None:
                elapsed = current_time - self.service_start_times[i]
                self.remaining_times[i] = max(0.0, self.remaining_times[i] - elapsed)
                self.service_start_times[i] = current_time

    def _find_preemptable_server(self, customer: Customer) -> int:
        """Find server to use (free or with later deadline)."""
        new_deadline = customer.absolute_deadline

        # Check for free servers
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return i

        # Check for preemption based on deadline
        worst_server = -1
        worst_deadline = new_deadline

        for i in range(self.num_servers):
            if self.deadlines[i] > worst_deadline:
                worst_deadline = self.deadlines[i]
                worst_server = i

        return worst_server

    def _preempt_and_serve(
        self,
        server_id: int,
        customer: Customer,
        current_time: float
    ) -> Optional[Customer]:
        preempted = self.in_service[server_id]

        if preempted is not None:
            preempted.service_time = self.remaining_times[server_id]
            heapq.heappush(
                self.wait_queue,
                (preempted.absolute_deadline, preempted.queue_arrival_time, preempted)
            )
            self._busy_servers[preempted.class_id] -= 1

        self.in_service[server_id] = customer
        self.service_start_times[server_id] = current_time
        self.remaining_times[server_id] = customer.service_time
        self.deadlines[server_id] = customer.absolute_deadline
        self._busy_servers[customer.class_id] += 1

        return preempted

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None
        self.remaining_times[server_id] = 0.0
        self.deadlines[server_id] = float('inf')

        self._update_remaining_times(current_time)

        if self.wait_queue:
            deadline, _, next_customer = heapq.heappop(self.wait_queue)
            self.in_service[server_id] = next_customer
            self.service_start_times[server_id] = current_time
            self.remaining_times[server_id] = next_customer.service_time
            self.deadlines[server_id] = deadline
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """EDFScheduler is a preemptive scheduler."""
        return True


class FCFSPRScheduler(SchedulingStrategy):
    """
    First-Come-First-Served Preemptive Resume (FCFSPR) scheduling.

    Similar to LCFSPR but uses FCFS order for the wait queue instead of
    LIFO. New arrivals can preempt if servers are busy, and preempted jobs
    resume with their remaining service time.

    If has_priority is True, only higher priority arrivals can preempt.
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # FCFS wait queue: List (FIFO with append/pop(0))
        self.wait_queue: List[Customer] = []

        # Stack of preempted jobs per server: server_id -> list of (preemption_record, customer)
        self.preemption_stacks: List[List[Tuple[PreemptionRecord, Customer]]] = [
            [] for _ in range(num_servers)
        ]

        # Currently active jobs per server
        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self.service_start_times: List[float] = [0.0] * num_servers

        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        """
        Handle arrival - may preempt current job.

        Returns:
            (accepted, service_info) where service_info is:
            - (server_id, customer, preempted_info) if service starts
            - None if queued
        """
        class_id = customer.class_id
        self._queue_length[class_id] += 1

        # Pre-sample service time if not already set
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # Find a server to use (preempt if needed)
        server_id, preempted = self._find_server_with_preemption(customer, current_time)

        if server_id >= 0:
            self.in_service[server_id] = customer
            self.service_start_times[server_id] = current_time
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer, preempted))

        # No server available (with priority, can't preempt) - queue
        self.wait_queue.append(customer)
        return (True, None)

    def _find_server_with_preemption(
        self,
        customer: Customer,
        current_time: float
    ) -> Tuple[int, Optional[Tuple[Customer, PreemptionRecord]]]:
        """
        Find a server, preempting if necessary.

        Returns:
            (server_id, preempted_info) where preempted_info is
            (preempted_customer, preemption_record) or None
        """
        # First, try to find a free server
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return (i, None)

        # No free server - must preempt if allowed
        preempt_server = -1
        preempt_priority = float('inf')

        for i in range(self.num_servers):
            current_job = self.in_service[i]
            if current_job is not None:
                if self.has_priority:
                    # Only preempt lower priority
                    if current_job.priority < customer.priority:
                        if current_job.priority < preempt_priority:
                            preempt_priority = current_job.priority
                            preempt_server = i
                else:
                    # FCFSPR without priority - always preempt (use first server)
                    preempt_server = 0
                    break

        if preempt_server >= 0:
            preempted_job = self.in_service[preempt_server]
            if preempted_job is not None:
                # Calculate remaining work
                elapsed = current_time - self.service_start_times[preempt_server]
                remaining = max(0.0, preempted_job.service_time - elapsed)

                record = PreemptionRecord(
                    customer_id=preempted_job.customer_id,
                    remaining_work=remaining,
                    preemption_time=current_time,
                    elapsed_time=elapsed
                )

                # Push onto preemption stack
                self.preemption_stacks[preempt_server].append((record, preempted_job))

                # Update busy server count
                self._busy_servers[preempted_job.class_id] -= 1

                return (preempt_server, (preempted_job, record))

        return (-1, None)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        """
        Handle departure - resume preempted job or serve from queue.

        Key difference from LCFSPR: Check preemption stack first,
        then serve from FCFS wait queue (not LIFO).

        Returns:
            Next customer to serve (may be a resumed job)
        """
        if isinstance(customer, Customer):
            class_id = customer.class_id
            self._queue_length[class_id] -= 1
            self._busy_servers[class_id] -= 1

        self.in_service[server_id] = None

        # First, check preemption stack for resumed job
        if self.preemption_stacks[server_id]:
            record, resumed_customer = self.preemption_stacks[server_id].pop()

            # Update service time to remaining work
            resumed_customer.service_time = record.remaining_work

            self.in_service[server_id] = resumed_customer
            self.service_start_times[server_id] = current_time
            self._busy_servers[resumed_customer.class_id] += 1

            return resumed_customer

        # Then check FCFS wait queue
        if self.wait_queue:
            next_customer = self.wait_queue.pop(0)  # FCFS: pop from front
            self.in_service[server_id] = next_customer
            self.service_start_times[server_id] = current_time
            self._busy_servers[next_customer.class_id] += 1

            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """FCFSPRScheduler is a preemptive scheduler."""
        return True


class FBScheduler(SchedulingStrategy):
    """
    Foreground-Background (FB) / Least Attained Service (LAS) scheduling.

    A preemptive scheduling discipline where the job with the least
    attained service (smallest accumulated service time) gets priority.
    This is equivalent to processor sharing for M/G/1 queues but handles
    discretely scheduled departures.

    New arrivals always have zero attained service, so they preempt
    any job that has received service. Among jobs with equal attained
    service, FCFS ordering is used.

    Properties:
    - Favors short jobs without knowing job sizes in advance
    - Approximates optimal for unknown job sizes
    - Fair in the sense that no job starves
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # Priority queue: (priority, attained_service, arrival_time, customer)
        # Lower attained service = higher priority
        self.wait_queue: List[Tuple[int, float, float, Customer]] = []

        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self.service_start_times: List[float] = [0.0] * num_servers
        self.attained_service: List[float] = [0.0] * num_servers
        self.remaining_times: List[float] = [0.0] * num_servers

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

        # Pre-sample service time
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # New arrivals have zero attained service - they preempt
        customer.attained_service = 0.0

        # Update attained service for all in-service jobs
        self._update_attained_service(current_time)

        # Check if we should preempt (new job has 0 attained service)
        server_id = self._find_preemptable_server(customer)

        if server_id >= 0:
            preempted = self._preempt_and_serve(server_id, customer, current_time)
            return (True, (server_id, customer, preempted))
        else:
            # Add to wait queue with attained service = 0
            priority = -customer.priority if self.has_priority else 0
            heapq.heappush(
                self.wait_queue,
                (priority, 0.0, customer.queue_arrival_time, customer)
            )
            return (True, None)

    def _update_attained_service(self, current_time: float) -> None:
        """Update attained service for all in-service jobs."""
        for i in range(self.num_servers):
            if self.in_service[i] is not None:
                elapsed = current_time - self.service_start_times[i]
                self.attained_service[i] += elapsed
                self.remaining_times[i] = max(0.0, self.remaining_times[i] - elapsed)
                self.service_start_times[i] = current_time

    def _find_preemptable_server(self, customer: Customer) -> int:
        """
        Find a server to use, either free or with higher attained service.

        FB always preempts the job with highest attained service when
        a new arrival (with 0 attained service) comes.

        Returns server_id or -1 if should queue.
        """
        new_priority = -customer.priority if self.has_priority else 0

        # First check for free servers
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return i

        # Find server with maximum attained service (candidate for preemption)
        worst_server = -1
        worst_attained = 0.0
        worst_priority = new_priority

        for i in range(self.num_servers):
            current = self.in_service[i]
            if current is not None:
                cur_priority = -current.priority if self.has_priority else 0
                cur_attained = self.attained_service[i]

                if self.has_priority:
                    # Can preempt if lower priority, or same priority and more attained
                    if (cur_priority > new_priority or
                        (cur_priority == new_priority and cur_attained > worst_attained)):
                        worst_priority = cur_priority
                        worst_attained = cur_attained
                        worst_server = i
                else:
                    # FB: preempt job with most attained service
                    if cur_attained > worst_attained:
                        worst_attained = cur_attained
                        worst_server = i

        # New job has 0 attained service, so preempt if any job has > 0
        if worst_server >= 0 and worst_attained > 0:
            return worst_server

        return -1

    def _preempt_and_serve(
        self,
        server_id: int,
        customer: Customer,
        current_time: float
    ) -> Optional[Customer]:
        """Preempt current job (if any) and start serving new customer."""
        preempted = self.in_service[server_id]

        if preempted is not None:
            # Save state and add to queue
            preempted.service_time = self.remaining_times[server_id]
            preempted.attained_service = self.attained_service[server_id]
            priority = -preempted.priority if self.has_priority else 0
            heapq.heappush(
                self.wait_queue,
                (priority, preempted.attained_service, preempted.queue_arrival_time, preempted)
            )
            self._busy_servers[preempted.class_id] -= 1

        # Start new customer
        self.in_service[server_id] = customer
        self.service_start_times[server_id] = current_time
        self.attained_service[server_id] = 0.0  # New job starts with 0
        self.remaining_times[server_id] = customer.service_time
        self._busy_servers[customer.class_id] += 1

        return preempted

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None
        self.remaining_times[server_id] = 0.0
        self.attained_service[server_id] = 0.0

        # Update attained service before selecting next
        self._update_attained_service(current_time)

        if self.wait_queue:
            # Select job with least attained service
            _, attained, _, next_customer = heapq.heappop(self.wait_queue)
            self.in_service[server_id] = next_customer
            self.service_start_times[server_id] = current_time
            self.attained_service[server_id] = attained
            self.remaining_times[server_id] = next_customer.service_time
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """FB is a preemptive scheduler."""
        return True


class SETFScheduler(SchedulingStrategy):
    """
    Shortest Elapsed Time First (SETF) - non-preemptive FB/LAS.

    Non-preemptive version of FB/LAS. When a server becomes free,
    the job with the least attained service (accumulated service time)
    is selected, but running jobs are not preempted.
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # Priority queue: (priority, attained_service, arrival_time, customer)
        self.wait_queue: List[Tuple[int, float, float, Customer]] = []

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

        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        customer.attained_service = 0.0

        # Find free server
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                self.in_service[i] = customer
                self._busy_servers[class_id] += 1
                return (True, (i, customer, None))

        # Queue with 0 attained service
        priority = -customer.priority if self.has_priority else 0
        heapq.heappush(
            self.wait_queue,
            (priority, 0.0, customer.queue_arrival_time, customer)
        )
        return (True, None)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None

        if self.wait_queue:
            # Select job with least attained service (always 0 for non-preemptive)
            _, _, _, next_customer = heapq.heappop(self.wait_queue)
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """SETF is non-preemptive."""
        return False


class FCFSPIScheduler(SchedulingStrategy):
    """
    First-Come-First-Served Preemptive Independent (FCFSPI) scheduling.

    Similar to FCFSPR, but when a preempted job resumes, it restarts
    with a fresh service time (independent of previous service).

    If has_priority is True, only higher priority arrivals can preempt.
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority

        # FCFS wait queue
        self.wait_queue: List[Customer] = []

        # Stack of preempted customers per server (no records needed for PI)
        self.preemption_stacks: List[List[Customer]] = [
            [] for _ in range(num_servers)
        ]

        self.in_service: List[Optional[Customer]] = [None] * num_servers
        self._queue_length = [0] * num_classes
        self._busy_servers = [0] * num_classes

        # Store service generator for restart
        self._service_gen: Optional[Callable[[int], float]] = None

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self._service_gen = service_gen

        # Pre-sample service time
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        server_id, preempted = self._find_server_with_preemption(customer)

        if server_id >= 0:
            self.in_service[server_id] = customer
            self._busy_servers[class_id] += 1
            return (True, (server_id, customer, preempted))

        # No server available (with priority, can't preempt) - queue
        self.wait_queue.append(customer)
        return (True, None)

    def _find_server_with_preemption(
        self,
        customer: Customer
    ) -> Tuple[int, Optional[Customer]]:
        """Find server, preempting if needed."""
        # Try free server first
        for i in range(self.num_servers):
            if self.in_service[i] is None:
                return (i, None)

        # Must preempt
        preempt_server = -1
        for i in range(self.num_servers):
            current_job = self.in_service[i]
            if current_job is not None:
                if self.has_priority:
                    if current_job.priority < customer.priority:
                        preempt_server = i
                        break
                else:
                    preempt_server = 0
                    break

        if preempt_server >= 0:
            preempted = self.in_service[preempt_server]
            if preempted is not None:
                # Push onto stack (service time will be resampled on resume)
                self.preemption_stacks[preempt_server].append(preempted)
                self._busy_servers[preempted.class_id] -= 1
                return (preempt_server, preempted)

        return (-1, None)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        """
        Handle departure - resume preempted job with new service time or serve from queue.

        Returns:
            Next customer to serve (may be a resumed job with new service time)
        """
        if isinstance(customer, Customer):
            self._queue_length[customer.class_id] -= 1
            self._busy_servers[customer.class_id] -= 1

        self.in_service[server_id] = None

        # First, check preemption stack for resumed job
        if self.preemption_stacks[server_id]:
            resumed = self.preemption_stacks[server_id].pop()

            # Resample service time for independent restart
            if self._service_gen is not None:
                resumed.service_time = self._service_gen(resumed.class_id)

            self.in_service[server_id] = resumed
            self._busy_servers[resumed.class_id] += 1
            return resumed

        # Then check FCFS wait queue
        if self.wait_queue:
            next_customer = self.wait_queue.pop(0)  # FCFS: pop from front
            self.in_service[server_id] = next_customer
            self._busy_servers[next_customer.class_id] += 1
            return next_customer

        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._busy_servers[class_id]

    def is_preemptive(self) -> bool:
        """FCFSPIScheduler is a preemptive scheduler."""
        return True
