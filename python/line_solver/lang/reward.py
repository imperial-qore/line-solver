"""
Reward - Factory class for common reward function templates.

Provides static methods for creating common reward functions used in CTMC analysis.
These templates are shortcuts for frequently-used metrics like queue length and utilization.

Example:
    >>> # Queue length reward for all classes
    >>> qlen = Reward.queue_length(queue1)
    >>>
    >>> # Utilization reward for specific class
    >>> util = Reward.utilization(queue1, class1)
    >>>
    >>> # Blocking probability
    >>> block = Reward.blocking(queue1)
"""

__all__ = ['Reward', 'RewardDescriptor']


class RewardDescriptor:
    """Callable reward function that also records what it measures.

    A RewardDescriptor wraps a reward callable together with a structural
    description of the reward, so that it can be serialized declaratively.
    It stays callable exactly like the plain lambda it replaces:

        >>> fn = Reward.queue_length(queue1)
        >>> value = fn(state)

    while additionally exposing ``kind``, ``node`` and ``jobclass``.

    A descriptor of kind ``Custom`` wraps an arbitrary user function and is
    deliberately NOT serializable: the writer warns and omits it rather than
    emitting a reward it cannot reproduce on reload.

    Attributes:
        kind: 'QLen' | 'Util' | 'Blocking' | 'Custom'
        node: Node object the reward refers to, or None
        jobclass: JobClass object the reward refers to, or None
        fn: the underlying callable
    """

    __slots__ = ('kind', 'node', 'jobclass', 'fn', '__wrapped__')

    def __init__(self, kind, node, jobclass, fn):
        if not callable(fn):
            raise ValueError('Reward descriptor must wrap a callable')
        self.kind = kind
        self.node = node
        self.jobclass = jobclass
        self.fn = fn
        # Callers introspect a reward with inspect.signature() to decide whether to
        # invoke it as fn(state) or fn(state, sn) (see SolverCTMC.getAvgReward). Without
        # __wrapped__, signature() would resolve to __call__(*args, **kwargs) and report
        # two parameters for every descriptor, so a one-argument reward would be called
        # with two and fail. Exposing the wrapped callable makes signature() report the
        # real arity.
        self.__wrapped__ = fn

    def __call__(self, *args, **kwargs):
        return self.fn(*args, **kwargs)

    def __deepcopy__(self, memo):
        # Mirror the semantics of the plain callables this class replaces:
        # copy.deepcopy treats functions as atomic, so a deep-copied Network
        # must keep sharing the same reward (and, crucially, must not clone the
        # Node/JobClass objects the descriptor refers to).
        return self

    def __repr__(self):
        node_name = getattr(self.node, 'name', None)
        class_name = getattr(self.jobclass, 'name', None)
        return "RewardDescriptor(kind=%r, node=%r, jobclass=%r)" % (self.kind, node_name, class_name)


class Reward:
    """Factory class for common reward function templates."""

    @staticmethod
    def _get_node_index(node):
        """Get node index, handling both wrapper and native implementations."""
        if hasattr(node, 'get_index'):
            return node.get_index()
        elif hasattr(node, 'index'):
            return node.index
        else:
            raise ValueError(f"Node {node.name} has no index attribute")

    @staticmethod
    def queue_length(node, jobclass=None):
        """Queue length reward template.

        Args:
            node: Station node object
            jobclass: (Optional) Job class object

        Returns:
            Callable that computes queue length reward

        Examples:
            >>> # Total jobs at queue1
            >>> model.set_reward('QLen', Reward.queue_length(queue1))
            >>>
            >>> # Class1 jobs at queue1
            >>> model.set_reward('QLen_C1', Reward.queue_length(queue1, class1))
        """
        if jobclass is None:
            # Return total jobs at node (all classes)
            return RewardDescriptor('QLen', node, None, lambda state: state.at(node).total())
        else:
            # Return jobs of specific class at node
            return RewardDescriptor('QLen', node, jobclass, lambda state: state.at(node, jobclass))

    @staticmethod
    def utilization(node, jobclass=None):
        """Server utilization reward template.

        Utilization is computed as min(jobs, nservers), representing the
        fraction of servers in use.

        Note: For M/M/1 queues, this simplifies to min(jobs, 1)

        Args:
            node: Station node object
            jobclass: (Optional) Job class object

        Returns:
            Callable that computes utilization reward

        Examples:
            >>> model.set_reward('Util', Reward.utilization(queue1))
            >>> model.set_reward('Util_C1', Reward.utilization(queue1, class1))
        """
        node_idx = Reward._get_node_index(node)

        if jobclass is None:
            # Total utilization at node (all classes)
            def util_fn(state, sn=None):
                jobs = state.at(node).total()
                if sn is None:
                    sn = state.sn
                station_idx = sn.nodeToStation[node_idx - 1]
                nservers = sn.nservers[station_idx]
                return min(jobs, nservers)
            return RewardDescriptor('Util', node, None, util_fn)
        else:
            # Utilization of specific class at node
            def util_class_fn(state, sn=None):
                jobs = state.at(node, jobclass)
                if sn is None:
                    sn = state.sn
                station_idx = sn.nodeToStation[node_idx - 1]
                nservers = sn.nservers[station_idx]
                return min(jobs, nservers)
            return RewardDescriptor('Util', node, jobclass, util_class_fn)

    @staticmethod
    def blocking(node):
        """Blocking probability reward template.

        Returns 1 if the node is at capacity, 0 otherwise.
        Useful for measuring congestion or capacity violations.

        Args:
            node: Station node object

        Returns:
            Callable that computes blocking probability

        Example:
            >>> model.set_reward('Block', Reward.blocking(queue1))
        """
        node_idx = Reward._get_node_index(node)

        def blocking_fn(state, sn=None):
            jobs = state.at(node).total()
            if sn is None:
                sn = state.sn
            station_idx = sn.nodeToStation[node_idx - 1]
            capacity = sn.cap[station_idx]
            return 1.0 if jobs >= capacity else 0.0
        return RewardDescriptor('Blocking', node, None, blocking_fn)

    @staticmethod
    def custom(user_fn):
        """Custom reward function wrapper.

        The returned value is callable exactly like ``user_fn``. It is marked
        as kind 'Custom', which is deliberately NOT serializable: an arbitrary
        user function cannot be reproduced from JSON, so the writer warns and
        omits it rather than emitting a wrong reward.

        Args:
            user_fn: Custom reward function

        Returns:
            RewardDescriptor wrapping ``user_fn``

        Example:
            >>> my_reward = lambda state: state.at(q1).total()**2
            >>> model.set_reward('Custom', Reward.custom(my_reward))
        """
        return RewardDescriptor('Custom', None, None, user_fn)
