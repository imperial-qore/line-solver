"""
Job class implementations for LINE queueing network models (pure Python).

This module provides concrete job class implementations including
open classes, closed classes, and specialized class types.
"""

from enum import Enum
from typing import Optional, TYPE_CHECKING

from .base import JobClass, JobClassType, Station

if TYPE_CHECKING:
    from .network import Network


class SignalType(Enum):
    """Types of signals for G-networks."""
    NEGATIVE = 'negative'      # Removes a job from destination queue
    REPLY = 'reply'            # Triggers a reply action
    CATASTROPHE = 'catastrophe'  # Removes ALL jobs from destination queue


class RemovalPolicy(Enum):
    """
    Removal policies for negative signals in G-networks.

    When a negative signal removes jobs from a queue, the removal policy
    determines which job(s) are selected for removal.

    Attributes:
        RANDOM: Select job uniformly at random from all jobs at the station
        FCFS: Remove the oldest job (first arrived)
        LCFS: Remove the newest job (last arrived)
    """
    RANDOM = 'random'  # Uniform random selection
    FCFS = 'fcfs'      # Remove oldest job first
    LCFS = 'lcfs'      # Remove newest job first


class OpenClass(JobClass):
    """
    Open job class for external arrivals.

    Jobs in an open class arrive from outside the system according to
    an arrival process and leave the system after completing service.

    Args:
        model: Parent network model.
        name: Name of the job class.
        prio: Priority level (default: 0).
        deadline: Soft deadline for response time (default: inf).
    """

    def __init__(self, model: 'Network', name: str, prio: int = 0,
                 deadline: float = float('inf')):
        super().__init__(model, name, JobClassType.OPEN, prio)
        self._deadline = deadline
        self._completes = True
        model.addClass(self)

    @property
    def completes(self) -> bool:
        """Whether jobs complete counting towards system throughput."""
        return self._completes

    @completes.setter
    def completes(self, value: bool):
        """Set whether jobs complete counting towards system throughput."""
        self._completes = bool(value)
        self._invalidate_java()

    def setCompletes(self, value: bool):
        """Set whether jobs complete counting towards system throughput."""
        self.completes = value

    def getCompletes(self) -> bool:
        """Get whether jobs complete counting towards system throughput."""
        return self._completes

    def getDeadline(self) -> float:
        """Get the soft deadline for this class."""
        return self._deadline

    def setDeadline(self, deadline: float):
        """Set the soft deadline for this class."""
        self._deadline = deadline
        self._invalidate_java()

    # Aliases
    set_completes = setCompletes
    get_completes = getCompletes
    get_deadline = getDeadline
    set_deadline = setDeadline


class OpenSignal(OpenClass):
    """
    Open signal class for G-networks and related models.

    OpenSignal is a specialized OpenClass that can have special effects
    on queues they visit, such as removing jobs (negative signals).
    For closed networks, use ClosedSignal instead.

    Args:
        model: Parent network model.
        name: Name of the signal class.
        signal_type: Type of signal (default: NEGATIVE).
        prio: Priority level (default: 0).
        removal_distribution: Distribution for batch removal count (default: None = remove 1).
        removal_policy: RemovalPolicy for selecting jobs (default: RANDOM).
        deadline: Soft deadline (default: inf).

    Reference:
        Gelenbe, E. (1991). "Product-form queueing networks with
        negative and positive customers", Journal of Applied Probability
    """

    def __init__(self, model: 'Network', name: str,
                 signal_type: SignalType,
                 prio: int = 0,
                 removal_distribution=None,
                 removal_policy: Optional[RemovalPolicy] = None,
                 deadline: float = float('inf')):
        """Create an OpenSignal class instance.

        Args:
            signal_type: SignalType constant (REQUIRED: NEGATIVE, REPLY, or CATASTROPHE)

        Raises:
            ValueError: If signal_type is not specified
        """
        if signal_type is None:
            raise ValueError("signal_type is required. Use SignalType.NEGATIVE, SignalType.REPLY, or SignalType.CATASTROPHE.")
        # Don't call OpenClass.__init__ to avoid double addClass
        JobClass.__init__(self, model, name, JobClassType.OPEN, prio)
        self._deadline = deadline
        self._completes = True
        self._signal_type = signal_type
        self._removal_distribution = removal_distribution
        self._removal_policy = removal_policy or RemovalPolicy.RANDOM
        model.addClass(self)

    @property
    def signal_type(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    @signal_type.setter
    def signal_type(self, value: SignalType):
        """Set the signal type."""
        self._signal_type = value
        self._invalidate_java()

    def getSignalType(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    def setSignalType(self, signal_type: SignalType):
        """Set the signal type."""
        self._signal_type = signal_type
        self._invalidate_java()

    @property
    def removal_distribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    @removal_distribution.setter
    def removal_distribution(self, value):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = value
        self._invalidate_java()

    def getRemovalDistribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    def setRemovalDistribution(self, distribution):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = distribution
        self._invalidate_java()

    @property
    def removal_policy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    @removal_policy.setter
    def removal_policy(self, value: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = value
        self._invalidate_java()

    def getRemovalPolicy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    def setRemovalPolicy(self, policy: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = policy
        self._invalidate_java()

    def is_catastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    def isCatastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    # Aliases
    get_signal_type = getSignalType
    set_signal_type = setSignalType
    get_removal_distribution = getRemovalDistribution
    set_removal_distribution = setRemovalDistribution
    get_removal_policy = getRemovalPolicy
    set_removal_policy = setRemovalPolicy


class Signal(JobClass):
    """
    Signal placeholder class that automatically resolves to OpenSignal or ClosedSignal.

    Signal is a placeholder class that users can use in both open and closed networks.
    The resolution to the concrete type (OpenSignal or ClosedSignal) happens during model
    finalization when getStruct() is called. This provides a simplified API where users
    don't need to explicitly choose between OpenSignal and ClosedSignal.

    Resolution rules:
        - If the network has a Source node: resolves to OpenSignal
        - If the network has no Source node: resolves to ClosedSignal

    For explicit control, users can directly use OpenSignal or ClosedSignal.

    Args:
        model: Parent network model.
        name: Name of the signal class.
        signal_type: Type of signal (REQUIRED: NEGATIVE, REPLY, or CATASTROPHE).
        prio: Priority level (default: 0).
        removal_distribution: Distribution for batch removal count (default: None = remove 1).
        removal_policy: RemovalPolicy for selecting jobs (default: RANDOM).

    Reference:
        Gelenbe, E. (1991). "Product-form queueing networks with
        negative and positive customers", Journal of Applied Probability
    """

    def __init__(self, model: 'Network', name: str,
                 signal_type: SignalType,
                 prio: int = 0,
                 removal_distribution=None,
                 removal_policy: Optional[RemovalPolicy] = None):
        """Create a Signal placeholder instance.

        Args:
            signal_type: SignalType constant (REQUIRED: NEGATIVE, REPLY, or CATASTROPHE)

        Raises:
            ValueError: If signal_type is not specified
        """
        if signal_type is None:
            raise ValueError("signal_type is required. Use SignalType.NEGATIVE, SignalType.REPLY, or SignalType.CATASTROPHE.")
        super().__init__(model, name, JobClassType.OPEN, prio)  # Default to OPEN, resolved later
        self._signal_type = signal_type
        self._removal_distribution = removal_distribution
        self._removal_policy = removal_policy or RemovalPolicy.RANDOM
        self._model = model
        model.addClass(self)

    def resolve(self, is_open: bool, refstat: Optional['Station'] = None) -> JobClass:
        """Resolve this Signal placeholder to OpenSignal or ClosedSignal.

        Args:
            is_open: True if the network is open (has Source node)
            refstat: Reference station for closed networks (ignored for open)

        Returns:
            OpenSignal or ClosedSignal instance
        """
        if is_open:
            concrete = OpenSignal(
                self._model, self.name, self._signal_type,
                self.priority, self._removal_distribution, self._removal_policy
            )
        else:
            concrete = ClosedSignal(
                self._model, self.name, self._signal_type, refstat,
                self.priority, self._removal_distribution, self._removal_policy
            )
        concrete._index = self._index
        return concrete

    @property
    def signal_type(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    @signal_type.setter
    def signal_type(self, value: SignalType):
        """Set the signal type."""
        self._signal_type = value
        self._invalidate_java()

    def getSignalType(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    def setSignalType(self, signal_type: SignalType):
        """Set the signal type."""
        self._signal_type = signal_type
        self._invalidate_java()

    @property
    def removal_distribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    @removal_distribution.setter
    def removal_distribution(self, value):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = value
        self._invalidate_java()

    def getRemovalDistribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    def setRemovalDistribution(self, distribution):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = distribution
        self._invalidate_java()

    @property
    def removal_policy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    @removal_policy.setter
    def removal_policy(self, value: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = value
        self._invalidate_java()

    def getRemovalPolicy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    def setRemovalPolicy(self, policy: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = policy
        self._invalidate_java()

    def is_catastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    def isCatastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    # Aliases
    get_signal_type = getSignalType
    set_signal_type = setSignalType
    get_removal_distribution = getRemovalDistribution
    set_removal_distribution = setRemovalDistribution
    get_removal_policy = getRemovalPolicy
    set_removal_policy = setRemovalPolicy


def Catastrophe(model: 'Network', name: str,
                prio: int = 0,
                removal_policy: Optional[RemovalPolicy] = None,
                deadline: float = float('inf')) -> Signal:
    """
    Create a catastrophe signal that removes ALL jobs from a queue.

    This is a convenience function that creates a Signal with SignalType.CATASTROPHE.
    Catastrophe signals, when arriving at a queue, remove all jobs present
    (empties the queue). This models disaster events or system resets.

    Args:
        model: Parent network model.
        name: Name of the catastrophe class.
        prio: Priority level (default: 0).
        removal_policy: RemovalPolicy for order of removal (default: RANDOM).
        deadline: Soft deadline (default: inf).

    Returns:
        Signal: A Signal instance with SignalType.CATASTROPHE

    Reference:
        Gelenbe, E. (1991). "Product-form queueing networks with
        negative and positive customers", Journal of Applied Probability
    """
    return Signal(model, name, SignalType.CATASTROPHE, prio, None, removal_policy, deadline)


class ClosedClass(JobClass):
    """
    Closed job class with fixed population.

    Jobs in a closed class circulate within the system with a fixed
    population. Jobs never leave the system but cycle through stations.

    Args:
        model: Parent network model.
        name: Name of the job class.
        njobs: Number of jobs (population) in this class.
        refstat: Reference station for this class.
        prio: Priority level (default: 0).
        deadline: Soft deadline for response time (default: inf).
    """

    def __init__(self, model: 'Network', name: str, njobs: int,
                 refstat: 'Station', prio: int = 0,
                 deadline: float = float('inf')):
        super().__init__(model, name, JobClassType.CLOSED, prio)
        self._njobs = njobs
        self._refstat = refstat  # Use consistent naming with base JobClass
        self._deadline = deadline
        self._completes = True
        model.addClass(self)

    def getNumberOfJobs(self) -> int:
        """Get the population of this closed class."""
        return self._njobs

    def setNumberOfJobs(self, njobs: int):
        """Set the population of this closed class."""
        self._njobs = njobs
        self._invalidate_java()

    def getPopulation(self) -> int:
        """Get the population (alias for getNumberOfJobs)."""
        return self._njobs

    @property
    def completes(self) -> bool:
        """Whether jobs complete counting towards system throughput."""
        return self._completes

    @completes.setter
    def completes(self, value: bool):
        """Set whether jobs complete counting towards system throughput."""
        self._completes = bool(value)
        self._invalidate_java()

    def setCompletes(self, value: bool):
        """Set whether jobs complete counting towards system throughput."""
        self.completes = value

    def getCompletes(self) -> bool:
        """Get whether jobs complete counting towards system throughput."""
        return self._completes

    def getDeadline(self) -> float:
        """Get the soft deadline for this class."""
        return self._deadline

    def setDeadline(self, deadline: float):
        """Set the soft deadline for this class."""
        self._deadline = deadline
        self._invalidate_java()

    # Aliases
    get_number_of_jobs = getNumberOfJobs
    set_number_of_jobs = setNumberOfJobs
    get_population = getPopulation
    population = property(getPopulation)
    number_of_jobs = property(getNumberOfJobs)
    set_completes = setCompletes
    get_completes = getCompletes
    get_deadline = getDeadline
    set_deadline = setDeadline


class ClosedSignal(ClosedClass):
    """
    Closed signal class for G-networks and related models.

    ClosedSignal is a specialized ClosedClass for signals that need to circulate
    in closed networks. Unlike OpenSignal, ClosedSignal can be used in networks
    without Source/Sink nodes.

    ClosedSignal has zero population - signals are created dynamically through
    class switching from the target job class.

    Args:
        model: Parent network model.
        name: Name of the signal class.
        signal_type: Type of signal (default: NEGATIVE).
        refstat: Reference station (should match target job class).
        prio: Priority level (default: 0).
        removal_distribution: Distribution for batch removal count (default: None = remove 1).
        removal_policy: RemovalPolicy for selecting jobs (default: RANDOM).

    Reference:
        Gelenbe, E. (1991). "Product-form queueing networks with
        negative and positive customers", Journal of Applied Probability
    """

    def __init__(self, model: 'Network', name: str,
                 signal_type: SignalType,
                 refstat: 'Station',
                 prio: int = 0,
                 removal_distribution=None,
                 removal_policy: Optional[RemovalPolicy] = None):
        """Create a ClosedSignal class instance.

        Args:
            signal_type: SignalType constant (REQUIRED: NEGATIVE, REPLY, or CATASTROPHE)
            refstat: Reference station (should match target job class)

        Raises:
            ValueError: If signal_type is not specified
        """
        if signal_type is None:
            raise ValueError("signal_type is required. Use SignalType.NEGATIVE, SignalType.REPLY, or SignalType.CATASTROPHE.")
        # ClosedSignal has 0 population - signals are created by class switching
        # Don't call ClosedClass.__init__ to avoid double addClass
        JobClass.__init__(self, model, name, JobClassType.CLOSED, prio)
        self._njobs = 0  # Zero population
        self._refstat = refstat
        self._deadline = float('inf')
        self._completes = True
        self._signal_type = signal_type
        self._removal_distribution = removal_distribution
        self._removal_policy = removal_policy or RemovalPolicy.RANDOM
        model.addClass(self)

    @property
    def signal_type(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    @signal_type.setter
    def signal_type(self, value: SignalType):
        """Set the signal type."""
        self._signal_type = value
        self._invalidate_java()

    def getSignalType(self) -> SignalType:
        """Get the signal type."""
        return self._signal_type

    def setSignalType(self, signal_type: SignalType):
        """Set the signal type."""
        self._signal_type = signal_type
        self._invalidate_java()

    @property
    def removal_distribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    @removal_distribution.setter
    def removal_distribution(self, value):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = value
        self._invalidate_java()

    def getRemovalDistribution(self):
        """Get the removal distribution for batch removal."""
        return self._removal_distribution

    def setRemovalDistribution(self, distribution):
        """Set the removal distribution for batch removal."""
        self._removal_distribution = distribution
        self._invalidate_java()

    @property
    def removal_policy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    @removal_policy.setter
    def removal_policy(self, value: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = value
        self._invalidate_java()

    def getRemovalPolicy(self) -> RemovalPolicy:
        """Get the removal policy for selecting jobs to remove."""
        return self._removal_policy

    def setRemovalPolicy(self, policy: RemovalPolicy):
        """Set the removal policy for selecting jobs to remove."""
        self._removal_policy = policy
        self._invalidate_java()

    def is_catastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    def isCatastrophe(self) -> bool:
        """Check if this signal is a catastrophe (removes all jobs)."""
        return self._signal_type == SignalType.CATASTROPHE

    # Aliases
    get_signal_type = getSignalType
    set_signal_type = setSignalType
    get_removal_distribution = getRemovalDistribution
    set_removal_distribution = setRemovalDistribution
    get_removal_policy = getRemovalPolicy
    set_removal_policy = setRemovalPolicy


class SelfLoopingClass(ClosedClass):
    """
    Self-looping closed class.

    Jobs in this class perpetually loop at their reference station,
    useful for modeling background workloads or special scheduling.

    Args:
        model: Parent network model.
        name: Name of the job class.
        njobs: Number of jobs in this class.
        refstat: Reference station where jobs loop.
        prio: Priority level (default: 0).
    """

    def __init__(self, model: 'Network', name: str, njobs: int,
                 refstat: 'Station', prio: int = 0):
        # Don't call ClosedClass.__init__ to control addClass
        JobClass.__init__(self, model, name, JobClassType.CLOSED, prio)
        self._njobs = njobs
        self._reference_station = refstat
        self._deadline = float('inf')
        self._completes = True
        self._is_self_looping = True
        model.addClass(self)


class DisabledClass(JobClass):
    """
    Disabled job class that perpetually loops at a reference station.

    Jobs in this class remain inactive and do not participate in normal
    network routing. All nodes are configured with DISABLED routing for
    this class.

    Args:
        model: Parent network model.
        name: Name of the disabled class.
        refstat: Reference station where jobs perpetually loop.
    """

    def __init__(self, model: 'Network', name: str, refstat: 'Station'):
        super().__init__(model, name, JobClassType.DISABLED, 0)
        self._reference_station = refstat
        self._njobs = 0
        model.addClass(self)
