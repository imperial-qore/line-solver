"""
Processor Sharing scheduling strategies.

This module implements processor sharing (PS) and its variants:
- PS: Egalitarian Processor Sharing
- DPS: Discriminatory Processor Sharing (weighted by job class)
- GPS: Generalized Processor Sharing (weighted by active jobs in class)

And priority variants:
- PSPRIO, DPSPRIO, GPSPRIO: Only most urgent priority jobs share processor
  (lower priority value = higher priority in LINE)
"""

from typing import List, Optional, Callable, Tuple, Any, Dict
from dataclasses import dataclass, field
import numpy as np

from .base import Customer, PSCustomer, SchedulingStrategy


@dataclass
class PSJob:
    """
    Job in a processor sharing queue.

    Attributes:
        customer: The customer being served
        remaining_work: Remaining service requirement
        last_update_time: Time of last rate update
        current_rate: Current service rate (share of capacity)
        class_weight: Weight for DPS/GPS
    """
    customer: Customer
    remaining_work: float
    last_update_time: float
    current_rate: float = 0.0
    class_weight: float = 1.0


class PSScheduler(SchedulingStrategy):
    """
    Egalitarian Processor Sharing (PS) scheduler.

    All jobs in the queue share the server capacity equally.
    If there are n jobs, each gets 1/n of the capacity.

    For multi-server: each job gets c/n where c is number of servers,
    but capped at 1 (cannot exceed full capacity).

    If has_priority is True (PSPRIO), only jobs at the most urgent
    priority level share the capacity (lower value = higher priority in LINE).
    """

    def __init__(self, num_classes: int, num_servers: int, has_priority: bool = False):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority
        self.capacity = float(num_servers)

        # All jobs currently in PS queue
        self.jobs: List[PSJob] = []

        self._queue_length = [0] * num_classes
        self._in_service = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        """
        Handle arrival - job immediately starts sharing capacity.
        """
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self._in_service[class_id] += 1

        # Pre-sample service time (total work requirement)
        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        # Update remaining work for all existing jobs before adding new one
        self._update_all_remaining_work(current_time)

        # Create new job
        job = PSJob(
            customer=customer,
            remaining_work=customer.service_time,
            last_update_time=current_time,
            class_weight=1.0
        )
        self.jobs.append(job)

        # Recalculate rates for all jobs
        self._update_rates()

        # Return the job with its expected completion time
        completion_time = self._get_completion_time(job, current_time)
        return (True, (job, completion_time))

    def _update_all_remaining_work(self, current_time: float) -> None:
        """Update remaining work for all jobs based on elapsed time."""
        for job in self.jobs:
            if job.current_rate > 0:
                elapsed = current_time - job.last_update_time
                work_done = elapsed * job.current_rate
                job.remaining_work = max(0.0, job.remaining_work - work_done)
            job.last_update_time = current_time

    def _update_rates(self) -> None:
        """Update service rates for all jobs."""
        if not self.jobs:
            return

        if self.has_priority:
            # Only most urgent priority jobs get service (lower value = higher priority in LINE)
            min_priority = min(job.customer.priority for job in self.jobs)
            active_jobs = [j for j in self.jobs if j.customer.priority == min_priority]
        else:
            active_jobs = self.jobs

        n_active = len(active_jobs)
        if n_active == 0:
            return

        # Each active job gets capacity / n_active, capped at 1
        rate_per_job = min(1.0, self.capacity / n_active)

        for job in self.jobs:
            if job in active_jobs:
                job.current_rate = rate_per_job
            else:
                job.current_rate = 0.0

    def _get_completion_time(self, job: PSJob, current_time: float) -> float:
        """Calculate expected completion time for a job at current rate."""
        if job.current_rate <= 0:
            return float('inf')
        return current_time + job.remaining_work / job.current_rate

    def get_next_departure(self, current_time: float) -> Tuple[Optional[PSJob], float]:
        """
        Get the job that will complete next and its completion time.

        Returns:
            (job, completion_time) or (None, inf) if no jobs
        """
        self._update_all_remaining_work(current_time)
        self._update_rates()

        if not self.jobs:
            return (None, float('inf'))

        # Find job with minimum completion time
        min_time = float('inf')
        min_job = None

        for job in self.jobs:
            if job.current_rate > 0:
                completion = current_time + job.remaining_work / job.current_rate
                if completion < min_time:
                    min_time = completion
                    min_job = job

        return (min_job, min_time)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        """
        Handle departure - remove job and recalculate rates.

        For PS, server_id is actually the PSJob that completed.
        """
        # Update all jobs first
        self._update_all_remaining_work(current_time)

        # Find and remove the completing job
        if isinstance(customer, PSJob):
            job = customer
            if job in self.jobs:
                self.jobs.remove(job)
                class_id = job.customer.class_id
                self._queue_length[class_id] -= 1
                self._in_service[class_id] -= 1
        elif isinstance(customer, Customer):
            # Find by customer
            for job in self.jobs:
                if job.customer is customer:
                    self.jobs.remove(job)
                    self._queue_length[customer.class_id] -= 1
                    self._in_service[customer.class_id] -= 1
                    break

        # Recalculate rates
        self._update_rates()

        # For PS, don't return a "next customer" - all are already being served
        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._in_service[class_id]

    def get_all_jobs(self) -> List[PSJob]:
        """Get all jobs currently in the system."""
        return list(self.jobs)


class DPSScheduler(SchedulingStrategy):
    """
    Discriminatory Processor Sharing (DPS) scheduler.

    Jobs share capacity in proportion to their class weights.
    A job of class k with weight w_k gets:
    rate_k = w_k * c / sum(w_j for all jobs j)

    If has_priority is True (DPSPRIO), only jobs at the most urgent
    priority level share capacity (lower value = higher priority in LINE).
    """

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        class_weights: Optional[List[float]] = None,
        has_priority: bool = False
    ):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority
        self.capacity = float(num_servers)

        if class_weights is None:
            self.class_weights = [1.0] * num_classes
        else:
            self.class_weights = list(class_weights)

        self.jobs: List[PSJob] = []
        self._queue_length = [0] * num_classes
        self._in_service = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self._in_service[class_id] += 1

        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        self._update_all_remaining_work(current_time)

        job = PSJob(
            customer=customer,
            remaining_work=customer.service_time,
            last_update_time=current_time,
            class_weight=self.class_weights[class_id]
        )
        self.jobs.append(job)

        self._update_rates()

        completion_time = self._get_completion_time(job, current_time)
        return (True, (job, completion_time))

    def _update_all_remaining_work(self, current_time: float) -> None:
        for job in self.jobs:
            if job.current_rate > 0:
                elapsed = current_time - job.last_update_time
                work_done = elapsed * job.current_rate
                job.remaining_work = max(0.0, job.remaining_work - work_done)
            job.last_update_time = current_time

    def _update_rates(self) -> None:
        if not self.jobs:
            return

        if self.has_priority:
            # Only most urgent priority jobs get service (lower value = higher priority in LINE)
            min_priority = min(job.customer.priority for job in self.jobs)
            active_jobs = [j for j in self.jobs if j.customer.priority == min_priority]
        else:
            active_jobs = self.jobs

        if not active_jobs:
            for job in self.jobs:
                job.current_rate = 0.0
            return

        # Total weight of active jobs
        total_weight = sum(job.class_weight for job in active_jobs)

        if total_weight <= 0:
            for job in self.jobs:
                job.current_rate = 0.0
            return

        for job in self.jobs:
            if job in active_jobs:
                # Rate proportional to weight
                job.current_rate = min(1.0, (job.class_weight / total_weight) * self.capacity)
            else:
                job.current_rate = 0.0

    def _get_completion_time(self, job: PSJob, current_time: float) -> float:
        if job.current_rate <= 0:
            return float('inf')
        return current_time + job.remaining_work / job.current_rate

    def get_next_departure(self, current_time: float) -> Tuple[Optional[PSJob], float]:
        self._update_all_remaining_work(current_time)
        self._update_rates()

        if not self.jobs:
            return (None, float('inf'))

        min_time = float('inf')
        min_job = None

        for job in self.jobs:
            if job.current_rate > 0:
                completion = current_time + job.remaining_work / job.current_rate
                if completion < min_time:
                    min_time = completion
                    min_job = job

        return (min_job, min_time)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        self._update_all_remaining_work(current_time)

        if isinstance(customer, PSJob):
            job = customer
            if job in self.jobs:
                self.jobs.remove(job)
                class_id = job.customer.class_id
                self._queue_length[class_id] -= 1
                self._in_service[class_id] -= 1
        elif isinstance(customer, Customer):
            for job in self.jobs:
                if job.customer is customer:
                    self.jobs.remove(job)
                    self._queue_length[customer.class_id] -= 1
                    self._in_service[customer.class_id] -= 1
                    break

        self._update_rates()
        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._in_service[class_id]


class GPSScheduler(SchedulingStrategy):
    """
    Generalized Processor Sharing (GPS) scheduler.

    Similar to DPS but the weight applies per-class rather than per-job.
    A job of class k with n_k jobs of that class gets:
    rate = (w_k / sum(w_j for active classes)) * c / n_k

    If has_priority is True (GPSPRIO), only jobs at the most urgent
    priority level share capacity (lower value = higher priority in LINE).
    """

    def __init__(
        self,
        num_classes: int,
        num_servers: int,
        class_weights: Optional[List[float]] = None,
        has_priority: bool = False
    ):
        super().__init__(num_classes, num_servers)
        self.has_priority = has_priority
        self.capacity = float(num_servers)

        if class_weights is None:
            self.class_weights = [1.0] * num_classes
        else:
            self.class_weights = list(class_weights)

        self.jobs: List[PSJob] = []
        self._queue_length = [0] * num_classes
        self._in_service = [0] * num_classes

    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        class_id = customer.class_id
        self._queue_length[class_id] += 1
        self._in_service[class_id] += 1

        if customer.service_time < 0:
            customer.service_time = service_gen(class_id)

        self._update_all_remaining_work(current_time)

        job = PSJob(
            customer=customer,
            remaining_work=customer.service_time,
            last_update_time=current_time,
            class_weight=self.class_weights[class_id]
        )
        self.jobs.append(job)

        self._update_rates()

        completion_time = self._get_completion_time(job, current_time)
        return (True, (job, completion_time))

    def _update_all_remaining_work(self, current_time: float) -> None:
        for job in self.jobs:
            if job.current_rate > 0:
                elapsed = current_time - job.last_update_time
                work_done = elapsed * job.current_rate
                job.remaining_work = max(0.0, job.remaining_work - work_done)
            job.last_update_time = current_time

    def _update_rates(self) -> None:
        if not self.jobs:
            return

        if self.has_priority:
            # Only most urgent priority jobs get service (lower value = higher priority in LINE)
            min_priority = min(job.customer.priority for job in self.jobs)
            active_jobs = [j for j in self.jobs if j.customer.priority == min_priority]
        else:
            active_jobs = self.jobs

        if not active_jobs:
            for job in self.jobs:
                job.current_rate = 0.0
            return

        # Count jobs per class among active jobs
        class_counts: Dict[int, int] = {}
        for job in active_jobs:
            cid = job.customer.class_id
            class_counts[cid] = class_counts.get(cid, 0) + 1

        # Total weight of active classes
        total_weight = sum(
            self.class_weights[cid] for cid in class_counts.keys()
        )

        if total_weight <= 0:
            for job in self.jobs:
                job.current_rate = 0.0
            return

        for job in self.jobs:
            if job in active_jobs:
                cid = job.customer.class_id
                class_share = self.class_weights[cid] / total_weight
                n_in_class = class_counts[cid]
                # Each job in the class shares the class's allocation
                job.current_rate = min(1.0, (class_share * self.capacity) / n_in_class)
            else:
                job.current_rate = 0.0

    def _get_completion_time(self, job: PSJob, current_time: float) -> float:
        if job.current_rate <= 0:
            return float('inf')
        return current_time + job.remaining_work / job.current_rate

    def get_next_departure(self, current_time: float) -> Tuple[Optional[PSJob], float]:
        self._update_all_remaining_work(current_time)
        self._update_rates()

        if not self.jobs:
            return (None, float('inf'))

        min_time = float('inf')
        min_job = None

        for job in self.jobs:
            if job.current_rate > 0:
                completion = current_time + job.remaining_work / job.current_rate
                if completion < min_time:
                    min_time = completion
                    min_job = job

        return (min_job, min_time)

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Customer]:
        self._update_all_remaining_work(current_time)

        if isinstance(customer, PSJob):
            job = customer
            if job in self.jobs:
                self.jobs.remove(job)
                class_id = job.customer.class_id
                self._queue_length[class_id] -= 1
                self._in_service[class_id] -= 1
        elif isinstance(customer, Customer):
            for job in self.jobs:
                if job.customer is customer:
                    self.jobs.remove(job)
                    self._queue_length[customer.class_id] -= 1
                    self._in_service[customer.class_id] -= 1
                    break

        self._update_rates()
        return None

    def get_queue_length(self, class_id: int) -> int:
        return self._queue_length[class_id]

    def get_busy_servers(self, class_id: int) -> int:
        return self._in_service[class_id]
