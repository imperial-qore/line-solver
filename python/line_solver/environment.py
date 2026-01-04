"""
Native Python implementation of Random Environment models.

This module provides classes for defining and analyzing queueing networks
in random environments, where the network parameters change according to
an underlying Markov modulated process.
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Optional, List, Dict, Any, Union, Callable, Tuple
import pandas as pd


class Environment:
    """
    A random environment model where a queueing network operates under
    different environmental conditions (stages).

    The environment switches between stages according to a Markov process,
    and each stage has its own network model with potentially different
    parameters.

    Example:
        >>> envModel = Environment('MyEnv', 2)
        >>> envModel.add_stage(0, 'Stage1', 'UP', network1)
        >>> envModel.add_stage(1, 'Stage2', 'DOWN', network2)
        >>> envModel.add_transition(0, 1, Exp(1.0))
        >>> envModel.add_transition(1, 0, Exp(0.5))
    """

    def __init__(self, name: str, num_stages: int = 0):
        """
        Create a random environment model.

        Args:
            name: Name of the environment model
            num_stages: Number of environmental stages (can be 0 if stages added later)
        """
        self.name = name
        self.num_stages = num_stages

        # Stage information
        self._stages: List[Dict[str, Any]] = []
        self._stage_names: List[str] = []
        self._stage_types: List[str] = []
        self._models: List[Any] = []

        # Transition information
        self._transitions: Dict[tuple, Any] = {}  # (from, to) -> distribution
        self._reset_rules: Dict[tuple, Callable] = {}  # (from, to) -> reset function

        # Initialize empty stages if num_stages provided
        for _ in range(num_stages):
            self._stages.append({})
            self._stage_names.append('')
            self._stage_types.append('')
            self._models.append(None)

    def add_stage(self, index: int, name: str, stage_type: str, model: Any) -> None:
        """
        Add or update a stage in the environment.

        Args:
            index: Stage index (0-based)
            name: Name of the stage
            stage_type: Type of stage ('UP', 'DOWN', etc.)
            model: Network model for this stage
        """
        # Expand arrays if needed
        while len(self._stages) <= index:
            self._stages.append({})
            self._stage_names.append('')
            self._stage_types.append('')
            self._models.append(None)

        self._stages[index] = {'name': name, 'type': stage_type, 'model': model}
        self._stage_names[index] = name
        self._stage_types[index] = stage_type
        self._models[index] = model

        self.num_stages = max(self.num_stages, index + 1)

    def add_transition(self, from_stage: int, to_stage: int, distribution: Any,
                        reset_rule: Optional[Callable[[np.ndarray], np.ndarray]] = None) -> None:
        """
        Add a transition between stages with an optional reset rule.

        Args:
            from_stage: Source stage index
            to_stage: Destination stage index
            distribution: Distribution for the transition time (e.g., Exp(rate))
            reset_rule: Optional function that transforms queue lengths when transition occurs.
                       The function takes a 2D numpy array q_exit[station, class] representing
                       queue lengths upon exiting the source stage, and returns a 2D numpy
                       array of the same shape representing the reset queue lengths upon
                       entering the destination stage. If None (default), queue lengths
                       are preserved (identity function).

        Example:
            # Move all jobs to station 0, preserving classes
            def move_to_station_0(q_exit):
                q_reset = np.zeros_like(q_exit)
                q_reset[0, :] = q_exit.sum(axis=0)
                return q_reset

            env_model.add_transition(0, 1, Exp(1.0), move_to_station_0)
        """
        self._transitions[(from_stage, to_stage)] = distribution
        if reset_rule is not None:
            self._reset_rules[(from_stage, to_stage)] = reset_rule
        else:
            # Default identity function
            self._reset_rules[(from_stage, to_stage)] = lambda q: q

    def get_stage(self, index: int) -> Dict[str, Any]:
        """Get stage information by index."""
        if 0 <= index < len(self._stages):
            return self._stages[index]
        return {}

    def get_model(self, index: int) -> Any:
        """Get the network model for a stage."""
        if 0 <= index < len(self._models):
            return self._models[index]
        return None

    def get_transition(self, from_stage: int, to_stage: int) -> Any:
        """Get the transition distribution between stages."""
        return self._transitions.get((from_stage, to_stage), None)

    def get_reset_rule(self, from_stage: int, to_stage: int) -> Optional[Callable[[np.ndarray], np.ndarray]]:
        """
        Get the reset rule for a transition between stages.

        Args:
            from_stage: Source stage index
            to_stage: Destination stage index

        Returns:
            The reset function for this transition, or None if no transition exists.
        """
        return self._reset_rules.get((from_stage, to_stage), None)

    def get_transition_rate_matrix(self) -> np.ndarray:
        """
        Build the transition rate matrix for the environment.

        Returns:
            numpy array of shape (num_stages, num_stages) with transition rates
        """
        E = self.num_stages
        Q = np.zeros((E, E))

        for (i, j), dist in self._transitions.items():
            if hasattr(dist, 'getRate'):
                rate = dist.getRate()
            elif hasattr(dist, 'rate'):
                rate = dist.rate
            elif hasattr(dist, 'getMean'):
                rate = 1.0 / dist.getMean()
            else:
                rate = float(dist)

            Q[i, j] = rate

        # Set diagonal to make rows sum to 0
        for i in range(E):
            Q[i, i] = -np.sum(Q[i, :])

        return Q

    def get_steady_state_probs(self) -> np.ndarray:
        """
        Compute steady-state probabilities for the environment stages.

        Returns:
            numpy array of stage probabilities
        """
        Q = self.get_transition_rate_matrix()
        E = self.num_stages

        if E == 0:
            return np.array([])

        # Solve pi * Q = 0, sum(pi) = 1
        # Replace last column with ones for normalization
        A = Q.T.copy()
        A[-1, :] = 1.0
        b = np.zeros(E)
        b[-1] = 1.0

        try:
            pi = np.linalg.solve(A, b)
        except np.linalg.LinAlgError:
            # If singular, use pseudo-inverse
            pi = np.linalg.lstsq(A, b, rcond=None)[0]

        return pi

    def stage_table(self) -> pd.DataFrame:
        """
        Get a table summarizing the environment stages.

        Returns:
            DataFrame with stage information
        """
        data = []
        for i, stage in enumerate(self._stages):
            row = {
                'Stage': i,
                'Name': self._stage_names[i] if i < len(self._stage_names) else '',
                'Type': self._stage_types[i] if i < len(self._stage_types) else '',
                'Model': self._models[i].name if i < len(self._models) and self._models[i] is not None and hasattr(self._models[i], 'name') else 'None'
            }
            data.append(row)

        df = pd.DataFrame(data)
        print(df.to_string(index=False))
        return df

    def getStageTable(self):
        """Alias for stage_table (MATLAB compatibility)."""
        return self.stage_table()

    def print_stage_table(self) -> None:
        """
        Print detailed stage table with transition information.
        Similar to Java's printStageTable() method.
        """
        print("Stage Table:")
        print("============")
        for i, stage in enumerate(self._stages):
            name = self._stage_names[i] if i < len(self._stage_names) else ''
            stage_type = self._stage_types[i] if i < len(self._stage_types) else ''
            model = self._models[i] if i < len(self._models) else None
            model_name = model.name if model is not None and hasattr(model, 'name') else 'None'

            # Count nodes and classes
            n_nodes = 0
            n_classes = 0
            if model is not None:
                if hasattr(model, 'get_nodes'):
                    n_nodes = len(model.get_nodes())
                elif hasattr(model, 'nodes'):
                    n_nodes = len(model.nodes)
                if hasattr(model, 'get_classes'):
                    n_classes = len(model.get_classes())
                elif hasattr(model, 'classes'):
                    n_classes = len(model.classes)

            print(f"Stage {i + 1}: {name} (Type: {stage_type})")
            print(f"  - Network: {model_name}")
            print(f"  - Nodes: {n_nodes}")
            print(f"  - Classes: {n_classes}")

        # Print transitions
        if self._transitions:
            print("Transitions:")
            for (from_idx, to_idx), dist in sorted(self._transitions.items()):
                from_name = self._stage_names[from_idx] if from_idx < len(self._stage_names) else f'Stage{from_idx}'
                to_name = self._stage_names[to_idx] if to_idx < len(self._stage_names) else f'Stage{to_idx}'

                # Get rate from distribution
                if hasattr(dist, 'getRate'):
                    rate = dist.getRate()
                elif hasattr(dist, 'rate'):
                    rate = dist.rate
                elif hasattr(dist, 'getMean'):
                    rate = 1.0 / dist.getMean()
                else:
                    rate = float(dist)

                print(f"  {from_name} -> {to_name}: rate = {rate:.4f}")

    @property
    def ensemble(self) -> List[Any]:
        """Get list of network models (MATLAB compatibility)."""
        return self._models


class SolverENV:
    """
    Solver for random environment models.

    This solver analyzes queueing networks operating in random environments
    by combining results from individual stage solvers with the environment
    steady-state probabilities.
    """

    def __init__(self, env_model: Environment, solvers: Union[List, Callable], options=None):
        """
        Create an environment solver.

        Args:
            env_model: Environment model to analyze
            solvers: Either a list of solvers (one per stage) or a factory function
            options: Solver options (dict or object with 'method' attribute)
        """
        self.env_model = env_model
        self.options = options or {}

        # Handle solvers as list or factory function
        if callable(solvers):
            # Factory function - create solver for each model
            self._solvers = []
            for model in env_model.ensemble:
                if model is not None:
                    self._solvers.append(solvers(model))
                else:
                    self._solvers.append(None)
        else:
            self._solvers = list(solvers)

        self._result = None
        self._smp_method = False  # Use DTMC-based computation for Semi-Markov Processes
        self._state_dep_method = ''  # State-dependent method configuration

        # Enable SMP method if specified in options
        method = None
        if isinstance(self.options, dict):
            method = self.options.get('method', None)
        elif hasattr(self.options, 'method'):
            method = self.options.method

        if method is not None and str(method).lower() == 'smp':
            self._smp_method = True

        # Auto-detect state-dependent environment
        # This is a placeholder - full implementation would check if resetEnvRates
        # functions are non-trivial, but Python implementation doesn't currently
        # support state-dependent rates in the same way as MATLAB/Java

        # Validate incompatible method combinations
        if self._smp_method and self._state_dep_method == 'statedep':
            raise ValueError(
                "SMP method (method='smp') is incompatible with state-dependent environments.\n"
                "SMP method computes environment probabilities once at initialization,\n"
                "but state-dependent environments modify transition rates during iterations.\n"
                "Please use either method='smp' OR state-dependent rates, but not both.")

    def setSMPMethod(self, flag: bool):
        """
        Enable/disable DTMC-based computation for Semi-Markov Processes.

        Args:
            flag: True to enable SMP method, False to disable
        """
        self._smp_method = flag

    def set_smp_method(self, flag: bool):
        """
        Enable/disable DTMC-based computation for Semi-Markov Processes (Pythonic name).

        Args:
            flag: True to enable SMP method, False to disable
        """
        self.setSMPMethod(flag)

    def avg(self):
        """
        Compute average performance metrics across environments.

        Returns:
            Tuple of (QN, UN, TN) - queue lengths, utilizations, throughputs
        """
        E = self.env_model.num_stages
        pi = self.env_model.get_steady_state_probs()

        if E == 0 or len(self._solvers) == 0:
            return np.array([]), np.array([]), np.array([])

        # Get metrics from each solver
        all_QN = []
        all_UN = []
        all_TN = []

        for solver in self._solvers:
            if solver is not None:
                try:
                    # Try to get results from solver
                    if hasattr(solver, 'getAvgQLen'):
                        QN = solver.getAvgQLen()
                    elif hasattr(solver, 'avg_qlen'):
                        QN = solver.avg_qlen()
                    else:
                        QN = np.array([0.0])

                    if hasattr(solver, 'getAvgUtil'):
                        UN = solver.getAvgUtil()
                    elif hasattr(solver, 'avg_util'):
                        UN = solver.avg_util()
                    else:
                        UN = np.array([0.0])

                    if hasattr(solver, 'getAvgTput'):
                        TN = solver.getAvgTput()
                    elif hasattr(solver, 'avg_tput'):
                        TN = solver.avg_tput()
                    else:
                        TN = np.array([0.0])

                    all_QN.append(QN)
                    all_UN.append(UN)
                    all_TN.append(TN)
                except Exception:
                    all_QN.append(np.array([0.0]))
                    all_UN.append(np.array([0.0]))
                    all_TN.append(np.array([0.0]))
            else:
                all_QN.append(np.array([0.0]))
                all_UN.append(np.array([0.0]))
                all_TN.append(np.array([0.0]))

        # Weighted average based on steady-state probabilities
        QN_avg = sum(pi[e] * all_QN[e] for e in range(E) if e < len(all_QN))
        UN_avg = sum(pi[e] * all_UN[e] for e in range(E) if e < len(all_UN))
        TN_avg = sum(pi[e] * all_TN[e] for e in range(E) if e < len(all_TN))

        self._result = (QN_avg, UN_avg, TN_avg)
        return QN_avg, UN_avg, TN_avg

    def avg_table(self) -> pd.DataFrame:
        """
        Get average metrics as a table.

        Returns:
            DataFrame with performance metrics
        """
        if self._result is None:
            self.avg()

        QN, UN, TN = self._result

        # Build table from first solver's model
        data = []
        for solver in self._solvers:
            if solver is not None and hasattr(solver, 'model'):
                model = solver.model
                if hasattr(model, 'get_nodes'):
                    nodes = model.get_nodes()
                    for i, node in enumerate(nodes):
                        node_name = node.name if hasattr(node, 'name') else f'Node{i}'
                        row = {
                            'Node': node_name,
                            'QLen': QN[i] if i < len(QN) else 0.0,
                            'Util': UN[i] if i < len(UN) else 0.0,
                            'Tput': TN[i] if i < len(TN) else 0.0,
                        }
                        data.append(row)
                break  # Only use first model for structure

        df = pd.DataFrame(data)
        return df

    @staticmethod
    def default_options():
        """Get default solver options."""
        from .solvers import OptionsDict
        return OptionsDict({
            'method': 'default',
            'iter_max': 100,
            'iter_tol': 0.01,
            'timespan': [0, float('inf')],
            'verbose': False
        })


# Convenience aliases
ENV = SolverENV


__all__ = [
    'Environment',
    'SolverENV',
    'ENV',
]
