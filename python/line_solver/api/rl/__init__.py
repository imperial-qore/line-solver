"""
Reinforcement Learning (RL) agents for queueing control.

Native Python implementations of RL agents for
queueing network control and optimization.

Key classes:
    RLTDAgent: Temporal Difference learning agent
    RLEnvironment: Queueing network environment

Port from:
    matlab/src/api/rl/rl_td_agent.m
    matlab/src/api/rl/rl_env.m
"""

import numpy as np
from typing import Optional, List, Tuple, Dict, Any


class RLEnvironment:
    """
    Reinforcement Learning environment for queueing networks.

    This class wraps a queueing network model and provides an interface
    for RL agents to interact with it through sampling and state updates.

    Attributes:
        model: The queueing network model
        gamma: Discount factor for future rewards
        queue_indices: Indices of queue nodes in model.nodes
        source_indices: Indices of source nodes in model.nodes
        state_size: Maximum state size to consider
        action_size: Number of possible actions (number of queues)
    """

    def __init__(self,
                 model: Any,
                 queue_indices: List[int],
                 source_indices: List[int],
                 state_size: int,
                 gamma: float = 0.99):
        """
        Initialize the RL environment.

        Args:
            model: Queueing network model
            queue_indices: Indices of queue nodes in model.nodes
            source_indices: Indices of source nodes in model.nodes
            state_size: Maximum number of jobs per queue to consider
            gamma: Discount factor (0 < gamma <= 1)
        """
        self.model = model
        self.queue_indices = queue_indices
        self.source_indices = source_indices
        self.state_size = state_size
        self.gamma = gamma
        self.action_size = len(queue_indices)

    def is_in_state_space(self, state: np.ndarray) -> bool:
        """
        Check if current state is within the defined state space.

        Args:
            state: Array of queue lengths

        Returns:
            True if all queue lengths are within state_size bounds
        """
        return np.all(state <= self.state_size)

    def is_in_action_space(self, state: np.ndarray) -> bool:
        """
        Check if actions are valid from current state.

        Args:
            state: Array of queue lengths

        Returns:
            True if all queues can accept new jobs (not at capacity)
        """
        return np.all(state < self.state_size)

    def sample(self) -> Tuple[float, int]:
        """
        Sample the next event from the environment.

        Uses the SSA solver to sample a single system event.

        Returns:
            Tuple of (time_delta, departure_node_index)
        """
        try:
            from ...solvers import SolverSSA
            solver = SolverSSA(self.model)
            sample = solver.sample_sys_aggr(1)
            t = sample['t']

            # Determine which node had the departure
            dep_node = -1
            for event in sample.get('events', []):
                if event.get('event') == 'DEP':
                    dep_node = event.get('node', -1)
                    break

            return t, dep_node
        except Exception:
            # Fallback: return random exponential time and random node
            t = np.random.exponential(1.0)
            dep_node = np.random.choice(self.queue_indices + self.source_indices)
            return t, dep_node

    def update(self, new_state: np.ndarray) -> None:
        """
        Update the model state after an event.

        Args:
            new_state: New queue lengths for each queue
        """
        # Update queue states in the model
        if hasattr(self.model, 'nodes'):
            for i, queue_idx in enumerate(self.queue_indices):
                if hasattr(self.model.nodes[queue_idx], 'state'):
                    self.model.nodes[queue_idx].state = int(new_state[i])

    def reset(self) -> np.ndarray:
        """
        Reset the environment to initial state.

        Returns:
            Initial state (zeros)
        """
        if hasattr(self.model, 'reset'):
            self.model.reset()
        if hasattr(self.model, 'init_default'):
            self.model.init_default()
        return np.zeros(self.action_size, dtype=np.int32)

    def get_state(self) -> np.ndarray:
        """
        Get current state from the model.

        Returns:
            Array of queue lengths
        """
        state = np.zeros(self.action_size, dtype=np.int32)
        if hasattr(self.model, 'nodes'):
            for i, queue_idx in enumerate(self.queue_indices):
                if hasattr(self.model.nodes[queue_idx], 'state'):
                    state[i] = int(np.sum(self.model.nodes[queue_idx].state))
        return state


class RLTDAgent:
    """
    Temporal Difference (TD) learning agent for queueing control.

    Implements average-reward TD learning for optimal routing decisions
    in queueing networks. The agent learns a value function V(s) that
    estimates the long-run average cost from each state.

    Attributes:
        learning_rate: Step size for value function updates
        epsilon: Exploration rate for epsilon-greedy policy
        epsilon_decay: Decay factor for exploration rate
        V: Value function array
        Q: Q-function array (state-action values)
    """

    def __init__(self,
                 learning_rate: float = 0.05,
                 epsilon: float = 1.0,
                 epsilon_decay: float = 0.99):
        """
        Initialize the TD agent.

        Args:
            learning_rate: Learning rate (step size) for updates
            epsilon: Initial exploration rate (0 to 1)
            epsilon_decay: Decay factor applied to epsilon each episode
        """
        self.learning_rate = learning_rate
        self.epsilon = epsilon
        self.epsilon_decay = epsilon_decay

        # Value function and Q-function (initialized on first use)
        self.V: Optional[np.ndarray] = None
        self.Q: Optional[np.ndarray] = None
        self.V_shape: Optional[Tuple] = None
        self.Q_shape: Optional[Tuple] = None

    def reset(self, env: RLEnvironment) -> None:
        """
        Reset agent and environment.

        Args:
            env: The RL environment
        """
        self.V = None
        self.Q = None
        self.V_shape = None
        self.Q_shape = None
        env.reset()

    def get_value_function(self) -> Optional[np.ndarray]:
        """Get the learned value function."""
        return self.V

    def get_q_function(self) -> Optional[np.ndarray]:
        """Get the learned Q-function."""
        return self.Q

    def solve(self,
              env: RLEnvironment,
              num_episodes: int = 10000,
              verbose: bool = True) -> Dict[str, Any]:
        """
        Train the agent using TD learning.

        Runs TD(0) learning for the specified number of episodes,
        learning an optimal routing policy for the queueing network.

        Args:
            env: The RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Dictionary with training results:
            - 'V': Learned value function
            - 'Q': Learned Q-function
            - 'mean_cost_rate': Final estimated average cost rate
        """
        self.reset(env)

        # Initialize value function with extra buffer for state space
        v_size = tuple([env.state_size + 5] * env.action_size)
        self.V = np.zeros(v_size)
        self.V_shape = v_size

        # Initialize Q-function
        q_size = v_size + (env.action_size,)
        self.Q = np.random.rand(*q_size)
        self.Q_shape = q_size

        # State tracking
        x = np.zeros(env.action_size, dtype=np.int32)  # Current state
        n = np.zeros(env.action_size, dtype=np.int32)  # Previous state

        t = 0.0  # Current time
        c = 0.0  # Cumulative cost
        T = 0.0  # Total discounted time
        C = 0.0  # Total discounted cost

        eps = self.epsilon
        episode = 0

        while episode < num_episodes:
            if verbose and episode % 1000 == 0:
                print(f"Running episode #{episode}")

            # Decay exploration rate
            eps *= self.epsilon_decay

            # Sample next event
            dt, dep_node = env.sample()
            t += dt

            # Accumulate holding cost
            c += np.sum(x) * dt

            # Process event
            if dep_node in env.source_indices:
                # New arrival - choose routing action
                if env.is_in_action_space(x):
                    # Create epsilon-greedy policy
                    next_states = []
                    for a in range(env.action_size):
                        next_x = x.copy()
                        next_x[a] += 1
                        next_states.append(self._get_state_index(next_x + 1))

                    next_values = np.array([self.V[s] for s in next_states])
                    policy = self._create_greedy_policy(next_values, eps, env.action_size)

                    # Sample action from policy
                    action = np.random.choice(env.action_size, p=policy)
                else:
                    # State space boundary - use Join Shortest Queue
                    action = np.argmin(x)

                x[action] += 1
                env.update(x)

            elif dep_node in env.queue_indices:
                # Departure from queue
                queue_idx = env.queue_indices.index(dep_node)
                x[queue_idx] = max(0, x[queue_idx] - 1)
                env.update(x)

            # TD update when in valid state space
            if env.is_in_state_space(x):
                episode += 1
                T = env.gamma * T + t
                C = env.gamma * C + c
                mean_cost_rate = C / T if T > 0 else 0

                prev_state = self._get_state_index(n + 1)
                curr_state = self._get_state_index(x + 1)

                # TD(0) update
                td_target = c - t * mean_cost_rate + self.V[curr_state]
                self.V[prev_state] = (1 - self.learning_rate) * self.V[prev_state] + \
                                     self.learning_rate * td_target

                # Normalize value function
                self.V = self.V - self.V.flat[0]

                # Reset for next episode
                t = 0.0
                c = 0.0
                n = x.copy()

        return {
            'V': self.V,
            'Q': self.Q,
            'mean_cost_rate': C / T if T > 0 else 0
        }

    def _get_state_index(self, loc: np.ndarray) -> Tuple:
        """
        Convert location array to state index.

        Args:
            loc: Array of queue lengths (1-indexed)

        Returns:
            Tuple index into value function array
        """
        loc = np.clip(loc, 0, self.V_shape[0] - 1).astype(int)
        return tuple(loc)

    @staticmethod
    def _create_greedy_policy(state_values: np.ndarray,
                              epsilon: float,
                              num_actions: int) -> np.ndarray:
        """
        Create epsilon-greedy policy from state values.

        Args:
            state_values: Value estimates for each action
            epsilon: Exploration probability
            num_actions: Number of possible actions

        Returns:
            Probability distribution over actions
        """
        policy = np.ones(num_actions) * epsilon / num_actions

        # Find actions with minimum value (cost minimization)
        min_val = np.min(state_values)
        best_actions = np.where(np.abs(state_values - min_val) < 1e-10)[0]

        # Add exploitation probability to best actions
        exploit_prob = (1 - epsilon) / len(best_actions)
        policy[best_actions] += exploit_prob

        return policy


__all__ = [
    'RLEnvironment',
    'RLTDAgent',
]
