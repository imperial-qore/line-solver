"""
MAPQN Solution Container.

Container class for storing results from MAPQN LP-based bounds computations.
"""

from dataclasses import dataclass, field
from typing import Dict, Optional


@dataclass
class MapqnSolution:
    """
    Container for MAPQN LP optimization results.

    Attributes:
        objective_value: Optimal value of the objective function.
        variables: Dictionary mapping variable names to their optimal values.

    Example:
        >>> sol = MapqnSolution(objective_value=0.75, variables={'U_1': 0.75, 'Q_1': 2.5})
        >>> sol.get_variable('U_1')
        0.75
        >>> sol.get_utilization(1)
        0.75
    """
    objective_value: float
    variables: Dict[str, float] = field(default_factory=dict)

    def get_variable(self, name: str) -> float:
        """
        Get the value of a variable by name.

        Args:
            name: Variable name (e.g., 'U_1', 'Q_2').

        Returns:
            Value of the variable, or 0.0 if not found.
        """
        return self.variables.get(name, 0.0)

    def get_utilization(self, i: int, k: Optional[int] = None) -> float:
        """
        Get utilization of queue i (optionally at phase k).

        Args:
            i: Queue index (1-based).
            k: Optional phase index (1-based).

        Returns:
            Utilization value.
        """
        if k is None:
            return self.get_variable(f'U_{i}')
        else:
            return self.get_variable(f'U_{i}_{k}')

    def get_queue_length(self, i: int, k: Optional[int] = None) -> float:
        """
        Get mean queue length of queue i (optionally at phase k).

        Args:
            i: Queue index (1-based).
            k: Optional phase index (1-based).

        Returns:
            Queue length value.
        """
        if k is None:
            return self.get_variable(f'Q_{i}')
        else:
            return self.get_variable(f'Q_{i}_{k}')

    def get_conditional_queue_length(self, j: int, i: int) -> float:
        """
        Get conditional queue length at i given customer at j.

        Args:
            j: Conditioning queue index (1-based).
            i: Target queue index (1-based).

        Returns:
            Conditional queue length value.
        """
        return self.get_variable(f'C_{j}_{i}')

    def get_marginal_probability(
        self, j: int, i: int, ni: int, ki: Optional[int] = None, hi: Optional[int] = None
    ) -> float:
        """
        Get marginal probability p1[j,i,ni] or p1[j,ki,i,hi,ni].

        Args:
            j: First queue index (1-based).
            i: Second queue index (1-based).
            ni: Population level.
            ki: Optional phase at queue j (1-based).
            hi: Optional phase at queue i (1-based).

        Returns:
            Marginal probability value.
        """
        if ki is None or hi is None:
            return self.get_variable(f'p1_{j}_{i}_{ni}')
        else:
            return self.get_variable(f'p1_{j}_{ki}_{i}_{hi}_{ni}')

    def __repr__(self) -> str:
        n_vars = len(self.variables)
        return f"MapqnSolution(objective_value={self.objective_value:.6f}, n_variables={n_vars})"


__all__ = ['MapqnSolution']
