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

__all__ = ['Reward']


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
            return lambda state: state.at(node).total()
        else:
            # Return jobs of specific class at node
            return lambda state: state.at(node, jobclass)

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
            return util_fn
        else:
            # Utilization of specific class at node
            def util_class_fn(state, sn=None):
                jobs = state.at(node, jobclass)
                if sn is None:
                    sn = state.sn
                station_idx = sn.nodeToStation[node_idx - 1]
                nservers = sn.nservers[station_idx]
                return min(jobs, nservers)
            return util_class_fn

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
        return blocking_fn

    @staticmethod
    def custom(user_fn):
        """Custom reward function wrapper.

        Provided for semantic clarity, making it explicit that
        the function is a custom user-defined reward.

        Args:
            user_fn: Custom reward function

        Returns:
            The wrapped function

        Example:
            >>> my_reward = lambda state: state.at(q1).total()**2
            >>> model.set_reward('Custom', Reward.custom(my_reward))
        """
        return user_fn
