"""
Base classes for scheduling disciplines.

This module provides the Customer dataclasses and abstract SchedulingStrategy
base class that all scheduling discipline implementations inherit from.
"""

from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import List, Optional, Any, Callable, Tuple
import random


@dataclass
class Customer:
    """
    Represents a job/customer in the queueing system.

    Used for non-preemptive and preemptive scheduling disciplines.

    Attributes:
        class_id: Job class index (0 to K-1)
        priority: Priority level (higher = more important)
        system_arrival_time: Time when job entered the system (from Source)
        queue_arrival_time: Time when job arrived at current station
        random_rank: Random value for SIRO discipline
        service_time: Pre-sampled service time (for SJF/LJF, -1 if not sampled)
        absolute_deadline: Deadline for EDD/EDF scheduling
        job_id: Unique job identifier
        parent_job_id: Parent job ID for forked jobs (Fork-Join)
        total_service_required: Total service time needed (Phase 4: preemption)
        remaining_service_work: Remaining service time after preemption (Phase 4)
        preemption_count: Number of times preempted (Phase 4)
        last_preemption_time: Time of last preemption event (Phase 4)
    """
    class_id: int
    priority: int
    system_arrival_time: float
    queue_arrival_time: float
    random_rank: float = field(default_factory=random.random)
    service_time: float = -1.0
    absolute_deadline: float = float('inf')
    job_id: int = -1
    parent_job_id: int = -1
    total_service_required: float = -1.0
    remaining_service_work: float = -1.0
    preemption_count: int = 0
    last_preemption_time: float = -1.0
    # Phase 7a: Reneging support
    patience_time: float = float('inf')  # How long customer will wait before reneging
    patience_start_time: float = -1.0    # Time customer arrived and started waiting
    has_reneged: bool = False             # True if customer gave up
    # FB/LAS scheduling support
    attained_service: float = 0.0         # Accumulated service time received

    def copy(self) -> 'Customer':
        """Create a copy of this customer."""
        customer = Customer(
            class_id=self.class_id,
            priority=self.priority,
            system_arrival_time=self.system_arrival_time,
            queue_arrival_time=self.queue_arrival_time,
            random_rank=self.random_rank,
            service_time=self.service_time,
            absolute_deadline=self.absolute_deadline,
            job_id=self.job_id,
            parent_job_id=self.parent_job_id,
            total_service_required=self.total_service_required,
            remaining_service_work=self.remaining_service_work,
            preemption_count=self.preemption_count,
            last_preemption_time=self.last_preemption_time,
        )
        # Copy reneging fields
        customer.patience_time = self.patience_time
        customer.patience_start_time = self.patience_start_time
        customer.has_reneged = self.has_reneged
        # Copy FB/LAS fields
        customer.attained_service = self.attained_service
        return customer

    @property
    def customer_id(self) -> int:
        """Alias for job_id for backward compatibility."""
        return self.job_id

    @customer_id.setter
    def customer_id(self, value: int) -> None:
        """Setter for customer_id that updates job_id."""
        self.job_id = value


@dataclass
class PSCustomer:
    """
    Customer for Processor Sharing scheduling disciplines.

    Extends Customer with remaining work tracking for PS/DPS/GPS.

    Attributes:
        class_id: Job class index
        priority: Priority level
        system_arrival_time: Time when job entered the system
        queue_arrival_time: Time when job arrived at current station
        total_service_requirement: Total service time needed
        remaining_service_work: Remaining service time
        departure_event: Reference to scheduled departure event (for cancellation)
        job_id: Unique job identifier
        parent_job_id: Parent job ID for forked jobs
    """
    class_id: int
    priority: int
    system_arrival_time: float
    queue_arrival_time: float
    total_service_requirement: float
    remaining_service_work: float
    departure_event: Optional[Any] = None
    job_id: int = -1
    parent_job_id: int = -1


@dataclass
class PreemptiveCustomer:
    """
    Customer for preemptive scheduling disciplines.

    Used for LCFSPR, LCFSPI, SRPT, EDF with preemption tracking.

    Attributes:
        class_id: Job class index
        priority: Priority level
        system_arrival_time: Time when job entered the system
        queue_arrival_time: Time when job arrived at current station
        random_rank: Random value for tie-breaking
        total_service_requirement: Total service time needed
        remaining_service_work: Remaining service time after preemption
        elapsed_service_time: Total time already spent in service
        service_start_time: Time when current service period started
        departure_event: Reference to scheduled departure event
        server_id: ID of server currently serving this job
        job_id: Unique job identifier
        parent_job_id: Parent job ID for forked jobs
    """
    class_id: int
    priority: int
    system_arrival_time: float
    queue_arrival_time: float
    random_rank: float
    total_service_requirement: float
    remaining_service_work: float
    elapsed_service_time: float = 0.0
    service_start_time: float = 0.0
    departure_event: Optional[Any] = None
    server_id: int = 0
    job_id: int = -1
    parent_job_id: int = -1


@dataclass
class PreemptionRecord:
    """
    Saved state for preemption resume.

    Stores the state of a preempted job for LCFSPR/SRPT resume.

    Attributes:
        remaining_work: Remaining service time when preempted
        original_total: Original total service requirement
        elapsed_time: Total time spent in service before preemption
    """
    remaining_work: float
    original_total: float
    elapsed_time: float


class SchedulingStrategy(ABC):
    """
    Abstract base class for all scheduling strategies.

    Scheduling strategies determine the order in which jobs receive service
    at a queue. Each strategy implements a specific queueing discipline
    (FCFS, LCFS, PS, etc.).

    Subclasses must implement:
        - arrive(): Handle job arrival at the queue
        - get_queue_length(): Get current queue length for a class
        - get_busy_servers(): Get number of busy servers for a class

    Optional overrides:
        - is_preemptive(): Whether strategy supports preemption
        - is_ps_family(): Whether strategy is processor sharing type
        - should_preempt(): Check if new arrival should preempt current
    """

    def __init__(self, num_classes: int, num_servers: int):
        """
        Initialize the scheduling strategy.

        Args:
            num_classes: Number of job classes
            num_servers: Number of servers at the station
        """
        self.num_classes = num_classes
        self.num_servers = num_servers

    @abstractmethod
    def arrive(
        self,
        customer: Customer,
        current_time: float,
        service_gen: Callable[[int], float],
    ) -> Tuple[bool, Optional[Any]]:
        """
        Handle customer arrival at the queue.

        Args:
            customer: The arriving customer
            current_time: Current simulation time
            service_gen: Function to generate service time for a class

        Returns:
            Tuple of (accepted: bool, service_process: Optional generator)
            - accepted: True if customer was accepted (not dropped)
            - service_process: SimPy process if service starts immediately, else None
        """
        pass

    @abstractmethod
    def get_queue_length(self, class_id: int) -> int:
        """
        Get current queue length for a specific class.

        Includes both waiting and in-service jobs.

        Args:
            class_id: Job class index

        Returns:
            Number of jobs of this class at the station
        """
        pass

    @abstractmethod
    def get_busy_servers(self, class_id: int) -> int:
        """
        Get number of servers busy with jobs of a specific class.

        Args:
            class_id: Job class index

        Returns:
            Number of servers currently serving this class
        """
        pass

    def is_preemptive(self) -> bool:
        """
        Whether this strategy supports preemption.

        Returns:
            True if arrivals can preempt jobs in service
        """
        return False

    def is_ps_family(self) -> bool:
        """
        Whether this is a processor sharing type strategy.

        PS family strategies share service capacity among all jobs
        simultaneously.

        Returns:
            True if this is PS, DPS, GPS, or a priority variant
        """
        return False

    def should_preempt(
        self,
        new_customer: Customer,
        in_service_customer: Any,
    ) -> bool:
        """
        Check if new customer should preempt the in-service customer.

        Only called if is_preemptive() returns True.

        Args:
            new_customer: The arriving customer
            in_service_customer: Customer currently in service

        Returns:
            True if new customer should preempt
        """
        return False

    def on_departure(
        self,
        customer: Any,
        current_time: float,
        server_id: int,
    ) -> Optional[Any]:
        """
        Handle customer departure (service completion).

        Called when a customer finishes service. Returns the next
        customer to serve (if any) or None.

        Args:
            customer: The departing customer
            current_time: Current simulation time
            server_id: ID of the server that was serving

        Returns:
            Next customer to serve, or None if queue is empty
        """
        return None

    def get_total_queue_length(self) -> int:
        """
        Get total queue length across all classes.

        Returns:
            Total number of jobs at the station
        """
        return sum(self.get_queue_length(k) for k in range(self.num_classes))

    def get_total_busy_servers(self) -> int:
        """
        Get total number of busy servers across all classes.

        Returns:
            Total number of busy servers
        """
        return sum(self.get_busy_servers(k) for k in range(self.num_classes))
