"""
Fork-Join node implementations.

This module provides Fork and Join nodes for parallel processing:
- Fork: Splits a job into multiple parallel tasks
- Join: Synchronizes parallel tasks with quorum support
"""

from typing import Dict, List, Optional, Set, Tuple, Any
from dataclasses import dataclass, field
import numpy as np

from ..scheduling.base import Customer


@dataclass
class ForkJobInfo:
    """
    Tracking information for a forked job.

    Attributes:
        parent_job_id: ID of the original job
        parent_customer: Original customer that was forked
        total_tasks: Number of parallel tasks created
        completed_tasks: Number of tasks that have completed
        completed_task_ids: Set of completed child task IDs
        fork_time: Time when the job was forked
        child_completion_times: Completion time of each child task
        quorum: Minimum tasks required for synchronization
    """
    parent_job_id: int
    parent_customer: Customer
    total_tasks: int
    completed_tasks: int = 0
    completed_task_ids: Set[int] = field(default_factory=set)
    fork_time: float = 0.0
    child_completion_times: Dict[int, float] = field(default_factory=dict)
    quorum: int = 0  # 0 means all tasks required

    def __post_init__(self):
        if self.quorum == 0:
            self.quorum = self.total_tasks

    def is_complete(self) -> bool:
        """Check if enough tasks have completed for synchronization."""
        return self.completed_tasks >= self.quorum


@dataclass
class ForkChild:
    """
    A child task created by forking.

    Attributes:
        child_id: Unique ID for this child
        parent_job_id: ID of the parent job
        task_index: Index in the set of parallel tasks (0 to n-1)
        destination_node: Node ID this child is sent to
        customer: Customer object for this child task
    """
    child_id: int
    parent_job_id: int
    task_index: int
    destination_node: int
    customer: Customer


class ForkNode:
    """
    Fork node that splits jobs into parallel tasks.

    When a job arrives at a Fork node, it creates multiple child
    tasks that are sent to different downstream nodes for parallel
    processing.
    """

    def __init__(
        self,
        node_id: int,
        fan_out: int,
        output_nodes: List[int],
        join_node_id: Optional[int] = None
    ):
        """
        Initialize Fork node.

        Args:
            node_id: Unique identifier for this fork node
            fan_out: Number of parallel tasks to create
            output_nodes: List of destination node IDs for each task
            join_node_id: ID of the corresponding Join node (optional)
        """
        self.node_id = node_id
        self.fan_out = fan_out
        self.output_nodes = output_nodes
        self.join_node_id = join_node_id

        # Track forked jobs
        self.fork_jobs: Dict[int, ForkJobInfo] = {}
        self.next_child_id = 0

    def fork(
        self,
        customer: Customer,
        current_time: float
    ) -> List[ForkChild]:
        """
        Fork a job into parallel tasks.

        Args:
            customer: The customer to fork
            current_time: Current simulation time

        Returns:
            List of ForkChild objects to route to parallel nodes
        """
        parent_id = customer.customer_id

        # Create tracking info
        fork_info = ForkJobInfo(
            parent_job_id=parent_id,
            parent_customer=customer,
            total_tasks=self.fan_out,
            fork_time=current_time
        )
        self.fork_jobs[parent_id] = fork_info

        # Create child tasks
        children = []
        for i in range(self.fan_out):
            dest = self.output_nodes[i] if i < len(self.output_nodes) else self.output_nodes[-1]

            # Create new customer for child task
            child_customer = Customer(
                job_id=self.next_child_id,
                class_id=customer.class_id,
                system_arrival_time=current_time,
                queue_arrival_time=current_time,
                service_time=-1.0,  # Will be sampled at destination
                priority=customer.priority,
                parent_job_id=parent_id  # Child tracks parent
            )

            child = ForkChild(
                child_id=self.next_child_id,
                parent_job_id=parent_id,
                task_index=i,
                destination_node=dest,
                customer=child_customer
            )
            children.append(child)
            self.next_child_id += 1

        return children

    def get_fork_info(self, parent_id: int) -> Optional[ForkJobInfo]:
        """Get fork info for a parent job."""
        return self.fork_jobs.get(parent_id)


class JoinNode:
    """
    Join node that synchronizes parallel tasks.

    Waits for child tasks to complete before releasing the
    synchronized job. Supports quorum-based synchronization
    where only a subset of tasks need to complete.
    """

    def __init__(
        self,
        node_id: int,
        fork_node: Optional[ForkNode] = None,
        quorum: int = 0  # 0 means wait for all
    ):
        """
        Initialize Join node.

        Args:
            node_id: Unique identifier for this join node
            fork_node: Reference to the corresponding Fork node
            quorum: Minimum tasks to wait for (0 = all tasks)
        """
        self.node_id = node_id
        self.fork_node = fork_node
        self.default_quorum = quorum

        # Track waiting jobs
        self.waiting_jobs: Dict[int, ForkJobInfo] = {}

        # Jobs ready to depart
        self.ready_queue: List[Customer] = []

    def arrive(
        self,
        child: ForkChild,
        current_time: float
    ) -> Optional[Customer]:
        """
        Handle arrival of a forked child task.

        Args:
            child: The child task that completed
            current_time: Current simulation time

        Returns:
            The synchronized parent customer if quorum is met, else None
        """
        parent_id = child.parent_job_id

        # Get or create fork info
        if parent_id not in self.waiting_jobs:
            if self.fork_node is not None:
                fork_info = self.fork_node.get_fork_info(parent_id)
                if fork_info is not None:
                    fork_info.quorum = self.default_quorum if self.default_quorum > 0 else fork_info.total_tasks
                    self.waiting_jobs[parent_id] = fork_info
                else:
                    # Create default tracking
                    self.waiting_jobs[parent_id] = ForkJobInfo(
                        parent_job_id=parent_id,
                        parent_customer=child.customer,  # Use child as placeholder
                        total_tasks=1,
                        quorum=1
                    )
            else:
                # No fork node reference - assume single task
                return child.customer

        fork_info = self.waiting_jobs[parent_id]

        # Mark this child as complete
        if child.child_id not in fork_info.completed_task_ids:
            fork_info.completed_task_ids.add(child.child_id)
            fork_info.completed_tasks += 1
            fork_info.child_completion_times[child.child_id] = current_time

        # Check if quorum is met
        if fork_info.is_complete():
            # Remove from waiting
            del self.waiting_jobs[parent_id]

            # Also clean up fork node tracking
            if self.fork_node is not None and parent_id in self.fork_node.fork_jobs:
                del self.fork_node.fork_jobs[parent_id]

            # Return the parent customer
            parent = fork_info.parent_customer
            # Update response time to be from fork time
            # This is handled by the caller
            return parent

        return None

    def get_waiting_count(self) -> int:
        """Get number of jobs waiting for synchronization."""
        return len(self.waiting_jobs)


class ForkJoinManager:
    """
    Manager for Fork-Join pairs in a network.

    Tracks Fork-Join relationships and handles synchronization
    across the network.
    """

    def __init__(self):
        self.fork_nodes: Dict[int, ForkNode] = {}
        self.join_nodes: Dict[int, JoinNode] = {}
        self.fork_to_join: Dict[int, int] = {}  # fork_id -> join_id

    def add_fork(
        self,
        node_id: int,
        fan_out: int,
        output_nodes: List[int],
        join_node_id: Optional[int] = None
    ) -> ForkNode:
        """Add a Fork node."""
        fork = ForkNode(node_id, fan_out, output_nodes, join_node_id)
        self.fork_nodes[node_id] = fork

        if join_node_id is not None:
            self.fork_to_join[node_id] = join_node_id

        return fork

    def add_join(
        self,
        node_id: int,
        fork_node_id: Optional[int] = None,
        quorum: int = 0
    ) -> JoinNode:
        """Add a Join node."""
        fork_node = self.fork_nodes.get(fork_node_id) if fork_node_id else None
        join = JoinNode(node_id, fork_node, quorum)
        self.join_nodes[node_id] = join
        return join

    def link_fork_join(self, fork_id: int, join_id: int) -> None:
        """Link a Fork node to its corresponding Join node."""
        if fork_id in self.fork_nodes and join_id in self.join_nodes:
            fork = self.fork_nodes[fork_id]
            join = self.join_nodes[join_id]

            fork.join_node_id = join_id
            join.fork_node = fork
            self.fork_to_join[fork_id] = join_id

    def get_fork(self, node_id: int) -> Optional[ForkNode]:
        """Get a Fork node by ID."""
        return self.fork_nodes.get(node_id)

    def get_join(self, node_id: int) -> Optional[JoinNode]:
        """Get a Join node by ID."""
        return self.join_nodes.get(node_id)


@dataclass
class SplitInfo:
    """
    Information for probabilistic split (router with fan-out).

    Attributes:
        probabilities: Probability of routing to each output
        output_nodes: Destination node IDs
    """
    probabilities: List[float]
    output_nodes: List[int]

    def select_output(self, rng: np.random.Generator) -> int:
        """Select output node probabilistically."""
        u = rng.random()
        cumsum = 0.0
        for i, p in enumerate(self.probabilities):
            cumsum += p
            if u < cumsum:
                return self.output_nodes[i]
        return self.output_nodes[-1]


class RouterNode:
    """
    Router node for probabilistic/deterministic routing.

    Supports various routing strategies:
    - PROB: Probabilistic routing based on routing matrix
    - RAND: Random selection among outputs
    - RROBIN: Round-robin across outputs
    - WRROBIN: Weighted round-robin
    - KCHOICES: Power-of-k choices (route to shortest queue)
    """

    def __init__(
        self,
        node_id: int,
        output_nodes: List[int],
        probabilities: Optional[List[float]] = None,
        strategy: str = 'PROB'
    ):
        self.node_id = node_id
        self.output_nodes = output_nodes
        self.strategy = strategy.upper()

        if probabilities is None:
            n = len(output_nodes)
            self.probabilities = [1.0 / n] * n
        else:
            self.probabilities = probabilities

        self.cumulative_probs = np.cumsum(self.probabilities)
        self.rrobin_index = 0

    def route(
        self,
        customer: Customer,
        rng: np.random.Generator,
        queue_lengths: Optional[List[int]] = None
    ) -> int:
        """
        Determine routing for a customer.

        Args:
            customer: Customer to route
            rng: Random number generator
            queue_lengths: Queue lengths at output nodes (for KCHOICES)

        Returns:
            Destination node ID
        """
        if self.strategy == 'PROB':
            u = rng.random()
            idx = np.searchsorted(self.cumulative_probs, u)
            return self.output_nodes[min(idx, len(self.output_nodes) - 1)]

        elif self.strategy == 'RAND':
            idx = rng.integers(0, len(self.output_nodes))
            return self.output_nodes[idx]

        elif self.strategy == 'RROBIN':
            dest = self.output_nodes[self.rrobin_index]
            self.rrobin_index = (self.rrobin_index + 1) % len(self.output_nodes)
            return dest

        elif self.strategy == 'KCHOICES':
            if queue_lengths is None:
                # Fall back to random
                idx = rng.integers(0, len(self.output_nodes))
                return self.output_nodes[idx]

            # Choose k random outputs, route to shortest queue
            k = min(2, len(self.output_nodes))
            choices = rng.choice(len(self.output_nodes), size=k, replace=False)
            best_idx = min(choices, key=lambda i: queue_lengths[i])
            return self.output_nodes[best_idx]

        else:
            # Default to first output
            return self.output_nodes[0]


class ClassSwitchNode:
    """
    Class switching node.

    Changes the class of a customer based on a switching matrix.
    """

    def __init__(
        self,
        node_id: int,
        num_classes: int,
        switch_matrix: Optional[np.ndarray] = None
    ):
        """
        Initialize ClassSwitch node.

        Args:
            node_id: Node identifier
            num_classes: Number of job classes
            switch_matrix: Transition matrix P[i,j] = P(new class = j | old class = i)
        """
        self.node_id = node_id
        self.num_classes = num_classes

        if switch_matrix is None:
            # Identity matrix (no switching)
            self.switch_matrix = np.eye(num_classes)
        else:
            self.switch_matrix = np.asarray(switch_matrix)

        # Precompute cumulative probabilities for each row
        self.cumulative_probs = np.cumsum(self.switch_matrix, axis=1)

    def switch_class(
        self,
        customer: Customer,
        rng: np.random.Generator
    ) -> int:
        """
        Determine new class for customer.

        Args:
            customer: Customer to switch
            rng: Random number generator

        Returns:
            New class ID
        """
        old_class = customer.class_id
        u = rng.random()
        new_class = np.searchsorted(self.cumulative_probs[old_class], u)
        return min(new_class, self.num_classes - 1)
