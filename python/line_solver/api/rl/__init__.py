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


class RLEnvironmentGeneral:
    """
    General RL environment for queueing networks with flexible action spaces.

    Unlike RLEnvironment which assumes sources dispatch to queues,
    this class supports arbitrary action nodes with configurable
    routing destinations. Actions are dispatching decisions at specific
    nodes that route jobs to connected downstream nodes.

    Attributes:
        model: The queueing network model
        gamma: Discount factor for future rewards
        queue_indices: Indices of queue nodes in model.nodes
        nqueues: Number of queues
        action_node_indices: Indices of nodes where routing actions are needed
        state_size: Maximum state size to consider
        action_space: Dict mapping action_node -> list of possible destination nodes

    MATLAB: matlab/src/api/rl/rl_env_general.m
    """

    def __init__(self,
                 model: Any,
                 queue_indices: List[int],
                 action_node_indices: List[int],
                 state_size: int,
                 gamma: float = 0.99):
        """
        Initialize the general RL environment.

        Args:
            model: Queueing network model
            queue_indices: Indices of queue nodes in model.nodes
            action_node_indices: Indices of nodes where dispatch actions are taken
            state_size: Maximum number of jobs per queue to consider
            gamma: Discount factor (0 < gamma <= 1)
        """
        self.model = model
        self.queue_indices = queue_indices
        self.nqueues = len(queue_indices)
        self.action_node_indices = action_node_indices
        self.state_size = state_size
        self.gamma = gamma

        # Build action space: for each action node, find connected destinations
        self.action_space: Dict[int, List[int]] = {}
        if hasattr(model, 'connections'):
            connections = model.connections
            for node_idx in action_node_indices:
                if hasattr(connections, '__getitem__'):
                    try:
                        row = connections[node_idx]
                        destinations = [j for j in range(len(row)) if row[j] == 1]
                        self.action_space[node_idx] = destinations
                    except (IndexError, TypeError):
                        self.action_space[node_idx] = list(queue_indices)
                else:
                    self.action_space[node_idx] = list(queue_indices)
        else:
            for node_idx in action_node_indices:
                self.action_space[node_idx] = list(queue_indices)

    def is_in_state_space(self, state: np.ndarray) -> bool:
        """
        Check if the given state is within the defined state space.

        Args:
            state: Array of queue lengths (one per queue)

        Returns:
            True if all queue lengths are within state_size bounds
        """
        if len(state) != self.nqueues:
            return False
        return np.all(state <= self.state_size)

    def is_in_action_space(self, state: np.ndarray) -> bool:
        """
        Check if actions are valid from the given state.

        Args:
            state: Array of queue lengths

        Returns:
            True if at least one queue can accept a new job
        """
        if len(state) != self.nqueues:
            return False
        return np.all(state < self.state_size)

    def sample(self) -> Tuple[float, int, int, Any]:
        """
        Sample the next event from the environment.

        Uses the SSA solver to sample a single system event.

        Returns:
            Tuple of (time_delta, departure_node, arrival_node, sample_data)
        """
        dep_node = -1
        arv_node = -1
        sample_data = None

        try:
            from ...solvers import SolverSSA
            solver = SolverSSA(self.model)
            sample_data = solver.sample_sys_aggr(1)
            dt = sample_data['t']

            for event in sample_data.get('events', []):
                etype = event.get('event', '')
                if etype == 'DEP':
                    dep_node = event.get('node', -1)
                elif etype == 'ARV':
                    arv_node = event.get('node', -1)

            return dt, dep_node, arv_node, sample_data
        except Exception:
            dt = np.random.exponential(1.0)
            dep_node = np.random.choice(self.queue_indices)
            arv_node = np.random.choice(self.queue_indices)
            return dt, dep_node, arv_node, sample_data

    def update(self, sample_data: Any) -> None:
        """
        Update the model state using the sample event data.

        Applies the state transitions from the sampled events to the model.

        Args:
            sample_data: Sample data from the SSA solver
        """
        if sample_data is None:
            return

        if hasattr(self.model, 'nodes') and hasattr(sample_data, 'get'):
            for event in sample_data.get('events', []):
                node_idx = event.get('node', -1)
                if node_idx >= 0 and node_idx < len(self.model.nodes):
                    node = self.model.nodes[node_idx]
                    if hasattr(node, 'state'):
                        etype = event.get('event', '')
                        if etype == 'DEP':
                            node.state = max(0, int(np.sum(node.state)) - 1)
                        elif etype == 'ARV':
                            node.state = int(np.sum(node.state)) + 1

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
        return np.zeros(self.nqueues, dtype=np.int32)


class RLTDAgentGeneral:
    """
    General TD learning agent for queueing control with flexible policies.

    Supports multiple solve methods:
    - solve(): TD control with tabular value function
    - solve_for_fixed_policy(): TD learning (evaluation) with fixed heuristic
    - solve_by_hashmap(): TD control with hash-map value function
    - solve_by_linear(): TD control with linear function approximation
    - solve_by_quad(): TD control with quadratic function approximation

    Works with RLEnvironmentGeneral for arbitrary action spaces.

    MATLAB: matlab/src/api/rl/rl_td_agent_general.m
    """

    def __init__(self,
                 learning_rate: float = 0.1,
                 epsilon: float = 1.0,
                 epsilon_decay: float = 0.9999):
        """
        Initialize the general TD agent.

        Args:
            learning_rate: Learning rate for updates
            epsilon: Initial exploration rate (0 to 1)
            epsilon_decay: Decay factor applied to epsilon each episode
        """
        self.learning_rate = learning_rate
        self.epsilon = epsilon
        self.epsilon_decay = epsilon_decay

        self.V: Optional[np.ndarray] = None
        self.V_shape: Optional[Tuple] = None

    def reset(self, env: RLEnvironmentGeneral) -> None:
        """
        Reset agent and environment.

        Args:
            env: The RL environment
        """
        self.V = None
        self.V_shape = None
        env.reset()

    def get_value_function(self) -> Optional[np.ndarray]:
        """Get the learned value function."""
        return self.V

    def solve_for_fixed_policy(self,
                                env: RLEnvironmentGeneral,
                                num_episodes: int = 10000,
                                verbose: bool = True) -> np.ndarray:
        """
        TD learning for value function evaluation with heuristic routing.

        Evaluates the value function under the existing (model-defined)
        routing policy without modifying routing decisions.

        Args:
            env: The general RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Learned value function array

        MATLAB: rl_td_agent_general.solve_for_fixed_policy
        """
        self.reset(env)

        v_size = tuple([env.state_size + 1] * env.nqueues)
        self.V = np.zeros(v_size)
        self.V_shape = v_size

        t = 0.0
        c = 0.0
        T = 0.0
        C = 0.0
        x = np.zeros(env.nqueues, dtype=np.int32)
        n_prev = np.zeros(env.nqueues, dtype=np.int32)

        episode = 0
        while episode < num_episodes:
            if verbose and episode % 1000 == 0:
                print(f"Running episode #{episode}")

            dt, dep_node, arv_node, sample_data = env.sample()
            t += dt
            c += np.sum(x) * dt

            # Process departure
            if dep_node in env.queue_indices:
                dep_server = env.queue_indices.index(dep_node)
                x[dep_server] = max(0, x[dep_server] - 1)

            # Process arrival
            if arv_node in env.queue_indices:
                arv_server = env.queue_indices.index(arv_node)
                x[arv_server] += 1

            env.update(sample_data)

            if env.is_in_state_space(x):
                episode += 1
                T = env.gamma * T + t
                C = env.gamma * C + c
                mean_cost_rate = C / T if T > 0 else 0.0

                prev_idx = tuple(np.clip(n_prev + 1, 0, self.V_shape[0] - 1).astype(int))
                curr_idx = tuple(np.clip(x + 1, 0, self.V_shape[0] - 1).astype(int))

                self.V[prev_idx] = ((1 - self.learning_rate) * self.V[prev_idx] +
                                    self.learning_rate * (c - t * mean_cost_rate + self.V[curr_idx]))
                self.V = self.V - self.V.flat[0]

                t = 0.0
                c = 0.0
                n_prev = x.copy()

        return self.V

    def solve(self,
              env: RLEnvironmentGeneral,
              num_episodes: int = 10000,
              verbose: bool = True) -> np.ndarray:
        """
        TD control with tabular value function.

        Learns an optimal routing policy using epsilon-greedy exploration.
        At each action node departure, the agent selects the best
        downstream destination based on the current value function.

        Args:
            env: The general RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Learned value function array

        MATLAB: rl_td_agent_general.solve
        """
        self.reset(env)

        v_size = tuple([env.state_size + 1] * env.nqueues)
        self.V = np.zeros(v_size)
        self.V_shape = v_size

        t = 0.0
        c = 0.0
        T = 0.0
        C = 0.0
        x = np.zeros(env.nqueues, dtype=np.int32)
        n_prev = np.zeros(env.nqueues, dtype=np.int32)

        eps = self.epsilon
        episode = 0

        while episode < num_episodes:
            if verbose and episode % 1000 == 0:
                print(f"Running episode #{episode}")

            eps *= self.epsilon_decay

            dt, dep_node, arv_node, sample_data = env.sample()
            t += dt
            c += np.sum(x) * dt

            # Process departure
            if dep_node in env.queue_indices:
                dep_server = env.queue_indices.index(dep_node)
                x[dep_server] = max(0, x[dep_server] - 1)

            # If departure is at an action node and we can act
            if dep_node in env.action_node_indices and env.is_in_action_space(x):
                actions = env.action_space.get(dep_node, [])
                if len(actions) > 0:
                    # Compute next values for each possible action
                    next_values = self._gen_next_values(env, x, actions)
                    policy = self._create_greedy_policy(next_values, eps, len(actions))

                    # Sample action
                    chosen_idx = np.random.choice(len(actions), p=policy)
                    arv_node = actions[chosen_idx]

                    # Update sample data to reflect the chosen arrival node
                    if sample_data is not None and hasattr(sample_data, 'get'):
                        for event in sample_data.get('events', []):
                            if event.get('event') == 'ARV':
                                event['node'] = arv_node
                                break

            # Process arrival
            if arv_node in env.queue_indices:
                arv_server = env.queue_indices.index(arv_node)
                x[arv_server] += 1

            env.update(sample_data)

            if env.is_in_state_space(x):
                episode += 1
                T = env.gamma * T + t
                C = env.gamma * C + c
                mean_cost_rate = C / T if T > 0 else 0.0

                prev_idx = tuple(np.clip(n_prev + 1, 0, self.V_shape[0] - 1).astype(int))
                curr_idx = tuple(np.clip(x + 1, 0, self.V_shape[0] - 1).astype(int))

                self.V[prev_idx] = ((1 - self.learning_rate) * self.V[prev_idx] +
                                    self.learning_rate * (c - t * mean_cost_rate + self.V[curr_idx]))
                self.V = self.V - self.V.flat[0]

                t = 0.0
                c = 0.0
                n_prev = x.copy()

        return self.V

    def solve_by_hashmap(self,
                          env: RLEnvironmentGeneral,
                          num_episodes: int = 10000,
                          verbose: bool = True) -> Tuple[np.ndarray, np.ndarray]:
        """
        TD control using hash-map value function.

        Uses a dictionary to store value estimates only for visited states,
        with an 'external' default for unvisited states.

        Args:
            env: The general RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Tuple of (X, Y) where:
                X: State features matrix (n_states x (1 + nqueues))
                Y: Value estimates (n_states x 1)

        MATLAB: rl_td_agent_general.solve_by_hashmap
        """
        self.reset(env)

        point_values: Dict[str, float] = {}
        zero_key = str(np.zeros(env.nqueues, dtype=int).tolist())
        point_values[zero_key] = 0.0
        point_values['external'] = 0.0

        t = 0.0
        c = 0.0
        T = 0.0
        C = 0.0
        x = np.zeros(env.nqueues, dtype=np.int32)
        n_prev = np.zeros(env.nqueues, dtype=np.int32)

        eps = self.epsilon
        episode = 0

        while episode < num_episodes:
            if verbose and episode % 1000 == 0:
                print(f"Running episode #{episode}")

            eps *= self.epsilon_decay

            dt, dep_node, arv_node, sample_data = env.sample()
            t += dt
            c += np.sum(x) * dt

            # Process departure
            if dep_node in env.queue_indices:
                dep_server = env.queue_indices.index(dep_node)
                x[dep_server] = max(0, x[dep_server] - 1)

            # Action at action nodes
            if dep_node in env.action_node_indices and env.is_in_action_space(x):
                actions = env.action_space.get(dep_node, [])
                if len(actions) > 0:
                    next_point_values = np.zeros(len(actions))
                    for act_i, act_node in enumerate(actions):
                        if act_node in env.queue_indices:
                            q_idx = env.queue_indices.index(act_node)
                            tmp_state = x.copy()
                            tmp_state[q_idx] += 1
                            key = str(tmp_state.tolist())
                            if key in point_values:
                                next_point_values[act_i] = point_values[key]
                            else:
                                next_point_values[act_i] = point_values['external']

                    policy = self._create_greedy_policy(next_point_values, eps, len(actions))
                    chosen_idx = np.random.choice(len(actions), p=policy)
                    arv_node = actions[chosen_idx]

                    if sample_data is not None and hasattr(sample_data, 'get'):
                        for event in sample_data.get('events', []):
                            if event.get('event') == 'ARV':
                                event['node'] = arv_node
                                break

            # Process arrival
            if arv_node in env.queue_indices:
                arv_server = env.queue_indices.index(arv_node)
                x[arv_server] += 1

            env.update(sample_data)

            if env.is_in_state_space(x):
                episode += 1
                T = env.gamma * T + t
                C = env.gamma * C + c
                mean_cost_rate = C / T if T > 0 else 0.0

                n_key = str(n_prev.tolist())
                x_key = str(x.tolist())

                if n_key not in point_values:
                    point_values[n_key] = point_values['external']

                if x_key in point_values:
                    point_values[n_key] = ((1 - self.learning_rate) * point_values[n_key] +
                                           self.learning_rate * (c - t * mean_cost_rate + point_values[x_key]))
                else:
                    point_values[n_key] = ((1 - self.learning_rate) * point_values[n_key] +
                                           self.learning_rate * (c - t * mean_cost_rate + point_values['external']))

                # Normalize relative to zero state
                if np.sum(n_prev) == 0:
                    subtractor = point_values[n_key]
                    for k in point_values:
                        point_values[k] -= subtractor

                t = 0.0
                c = 0.0
                n_prev = x.copy()

        # Build output arrays
        if 'external' in point_values:
            del point_values['external']

        n_entries = len(point_values)
        X = np.zeros((n_entries, 1 + env.nqueues))
        Y = np.zeros((n_entries, 1))

        for idx, (key, val) in enumerate(point_values.items()):
            state = np.array(eval(key))  # Parse state from string key
            X[idx, 0] = 1.0  # Intercept
            X[idx, 1:] = state
            Y[idx, 0] = val

        return X, Y

    def solve_by_linear(self,
                         env: RLEnvironmentGeneral,
                         num_episodes: int = 10000,
                         verbose: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        TD control with linear function approximation.

        Learns a linear value function: v(q1,...,qn) = w0 + w1*q1 + ... + wn*qn

        Args:
            env: The general RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Tuple of (X, Y, coefficients) where:
                X: State features matrix
                Y: Value estimates
                coefficients: Linear regression coefficients

        MATLAB: rl_td_agent_general.solve_by_linear
        """
        X, Y = self.solve_by_hashmap(env, num_episodes, verbose)

        # Linear regression: Y = X * coeff
        try:
            coeff, _, _, _ = np.linalg.lstsq(X, Y, rcond=None)
        except Exception:
            coeff = np.zeros((X.shape[1], 1))

        return X, Y, coeff

    def solve_by_quad(self,
                       env: RLEnvironmentGeneral,
                       num_episodes: int = 10000,
                       verbose: bool = True) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        TD control with quadratic function approximation.

        Learns a quadratic value function:
        v(q1,...,qn) = sum_{i,j} w_{ij} * q_i * q_j

        Args:
            env: The general RL environment
            num_episodes: Number of training episodes
            verbose: Whether to print progress

        Returns:
            Tuple of (X_quad, Y, coefficients) where:
                X_quad: Augmented features (linear + quadratic terms)
                Y: Value estimates
                coefficients: Quadratic regression coefficients

        MATLAB: rl_td_agent_general.solve_by_quad
        """
        X, Y = self.solve_by_hashmap(env, num_episodes, verbose)

        # Augment X with quadratic terms
        n_cols = X.shape[1]
        quad_features = []
        for i in range(1, n_cols):  # Skip intercept
            for j in range(i, n_cols):
                quad_features.append(X[:, i] * X[:, j])

        if quad_features:
            X_quad = np.hstack([X] + [f.reshape(-1, 1) for f in quad_features])
        else:
            X_quad = X

        try:
            coeff, _, _, _ = np.linalg.lstsq(X_quad, Y, rcond=None)
        except Exception:
            coeff = np.zeros((X_quad.shape[1], 1))

        return X_quad, Y, coeff

    def _gen_next_values(self, env: RLEnvironmentGeneral,
                          cur_state: np.ndarray,
                          actions: List[int]) -> np.ndarray:
        """
        Compute value estimates for each possible action from current state.

        Args:
            env: Environment
            cur_state: Current queue lengths
            actions: List of possible action (destination node) indices

        Returns:
            Array of value estimates for each action
        """
        values = np.zeros(len(actions))
        for act_i, act_node in enumerate(actions):
            if act_node in env.queue_indices:
                q_idx = env.queue_indices.index(act_node)
                tmp_loc = cur_state.copy() + 1
                tmp_loc[q_idx] += 1
                tmp_idx = tuple(np.clip(tmp_loc, 0, self.V_shape[0] - 1).astype(int))
                values[act_i] = self.V[tmp_idx]
        return values

    @staticmethod
    def _create_greedy_policy(state_values: np.ndarray,
                               epsilon: float,
                               num_actions: int) -> np.ndarray:
        """
        Create epsilon-greedy policy from state-action values.

        Args:
            state_values: Value estimates for each action
            epsilon: Exploration probability
            num_actions: Number of possible actions

        Returns:
            Probability distribution over actions
        """
        policy = np.ones(num_actions) * epsilon / num_actions
        min_val = np.min(state_values)
        best_actions = np.where(np.abs(state_values - min_val) < 1e-10)[0]
        exploit_prob = (1 - epsilon) / len(best_actions)
        policy[best_actions] += exploit_prob
        return policy


__all__ = [
    'RLEnvironment',
    'RLTDAgent',
    'RLEnvironmentGeneral',
    'RLTDAgentGeneral',
]
