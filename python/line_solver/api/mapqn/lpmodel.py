"""
MAPQN Linear Programming Model.

Base class for representing MAP queueing network linear programming models.
Provides the foundation for LP-based optimization methods in MAPQN analysis,
including constraint formulation and objective function definition.

Uses scipy.optimize.linprog as the LP solver backend.
"""

import numpy as np
from scipy.optimize import linprog
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

from .solution import MapqnSolution


@dataclass
class LinearConstraint:
    """
    Represents a single linear constraint.

    Attributes:
        coefficients: Coefficient array for the constraint.
        relationship: One of 'eq' (=), 'leq' (<=), or 'geq' (>=).
        rhs: Right-hand side value.
    """
    coefficients: np.ndarray
    relationship: str  # 'eq', 'leq', 'geq'
    rhs: float


class LinearConstraintBuilder:
    """
    Helper class for building linear constraints incrementally.

    Example:
        >>> builder = model.constraint_builder()
        >>> builder.add_term('U_1', 1.0).add_term('U_2', -1.0).eq(0.0)
    """

    def __init__(self, model: 'MapqnLpModel'):
        self.model = model
        self.coefficients = np.zeros(model.get_num_variables())

    def add_term(self, var_name: str, coefficient: float) -> 'LinearConstraintBuilder':
        """Add a term to the constraint."""
        index = self.model.get_variable_index(var_name)
        if index is not None:
            self.coefficients[index] += coefficient
        return self

    def add_term_by_index(self, var_index: int, coefficient: float) -> 'LinearConstraintBuilder':
        """Add a term by variable index."""
        if 0 <= var_index < len(self.coefficients):
            self.coefficients[var_index] += coefficient
        return self

    def eq(self, rhs: float) -> LinearConstraint:
        """Create equality constraint (=)."""
        return LinearConstraint(self.coefficients.copy(), 'eq', rhs)

    def leq(self, rhs: float) -> LinearConstraint:
        """Create less-than-or-equal constraint (<=)."""
        return LinearConstraint(self.coefficients.copy(), 'leq', rhs)

    def geq(self, rhs: float) -> LinearConstraint:
        """Create greater-than-or-equal constraint (>=)."""
        return LinearConstraint(self.coefficients.copy(), 'geq', rhs)


class MapqnLpModel:
    """
    Base class for MAPQN Linear Programming models.

    Provides methods for:
    - Variable registration and tracking
    - Constraint addition
    - Objective function creation
    - LP solving via scipy

    Example:
        >>> model = MapqnLpModel()
        >>> model.add_variable('U_1')
        >>> model.add_variable('Q_1')
        >>> constraint = model.constraint_builder().add_term('U_1', 1.0).leq(1.0)
        >>> model.add_constraint(constraint)
        >>> solution = model.solve('U_1', minimize=True)
    """

    def __init__(self):
        self._constraints: List[LinearConstraint] = []
        self._variables: Dict[str, int] = {}
        self._variable_counter = 0
        self._bounds: Dict[int, Tuple[Optional[float], Optional[float]]] = {}

    def add_variable(self, name: str, lb: Optional[float] = 0.0,
                     ub: Optional[float] = None) -> int:
        """
        Register a variable and return its index.

        Args:
            name: Variable name.
            lb: Lower bound (default 0.0, use None for unbounded).
            ub: Upper bound (default None for unbounded).

        Returns:
            Index of the variable.
        """
        if name in self._variables:
            return self._variables[name]

        index = self._variable_counter
        self._variables[name] = index
        self._bounds[index] = (lb, ub)
        self._variable_counter += 1
        return index

    def get_variable_index(self, name: str) -> Optional[int]:
        """Get variable index by name, or None if not found."""
        return self._variables.get(name)

    def get_num_variables(self) -> int:
        """Get total number of registered variables."""
        return self._variable_counter

    def add_constraint(self, constraint: LinearConstraint) -> None:
        """Add a constraint to the model."""
        # Resize constraint coefficients if needed
        if len(constraint.coefficients) < self._variable_counter:
            new_coeffs = np.zeros(self._variable_counter)
            new_coeffs[:len(constraint.coefficients)] = constraint.coefficients
            constraint.coefficients = new_coeffs
        self._constraints.append(constraint)

    def get_constraints(self) -> List[LinearConstraint]:
        """Get all constraints."""
        return self._constraints.copy()

    def constraint_builder(self) -> LinearConstraintBuilder:
        """Create a new constraint builder."""
        return LinearConstraintBuilder(self)

    def create_objective_coefficients(self, objective_var: str) -> np.ndarray:
        """
        Create objective function coefficients for a single variable.

        Args:
            objective_var: Name of the variable to optimize.

        Returns:
            Coefficient array with 1.0 at the variable's position.
        """
        coeffs = np.zeros(self._variable_counter)
        index = self.get_variable_index(objective_var)
        if index is not None:
            coeffs[index] = 1.0
        return coeffs

    def solve(self, objective_var: str, minimize: bool = True,
              method: str = 'highs') -> MapqnSolution:
        """
        Solve the LP model.

        Args:
            objective_var: Variable name to optimize.
            minimize: If True, minimize; if False, maximize.
            method: LP solver method ('highs', 'simplex', 'interior-point').

        Returns:
            MapqnSolution containing optimal value and variable values.

        Raises:
            RuntimeError: If LP is infeasible or unbounded.
        """
        n_vars = self._variable_counter
        if n_vars == 0:
            return MapqnSolution(objective_value=0.0, variables={})

        # Build objective coefficients
        c = self.create_objective_coefficients(objective_var)
        if not minimize:
            c = -c

        # Separate equality and inequality constraints
        A_eq_list = []
        b_eq_list = []
        A_ub_list = []
        b_ub_list = []

        for constraint in self._constraints:
            # Ensure constraint has correct size
            coeffs = constraint.coefficients
            if len(coeffs) < n_vars:
                coeffs = np.concatenate([coeffs, np.zeros(n_vars - len(coeffs))])
            elif len(coeffs) > n_vars:
                coeffs = coeffs[:n_vars]

            if constraint.relationship == 'eq':
                A_eq_list.append(coeffs)
                b_eq_list.append(constraint.rhs)
            elif constraint.relationship == 'leq':
                A_ub_list.append(coeffs)
                b_ub_list.append(constraint.rhs)
            elif constraint.relationship == 'geq':
                # Convert >= to <= by negating
                A_ub_list.append(-coeffs)
                b_ub_list.append(-constraint.rhs)

        # Build matrices
        A_eq = np.array(A_eq_list) if A_eq_list else None
        b_eq = np.array(b_eq_list) if b_eq_list else None
        A_ub = np.array(A_ub_list) if A_ub_list else None
        b_ub = np.array(b_ub_list) if b_ub_list else None

        # Build bounds
        bounds = []
        for i in range(n_vars):
            if i in self._bounds:
                bounds.append(self._bounds[i])
            else:
                bounds.append((0.0, None))  # Default: non-negative

        # Solve LP
        try:
            result = linprog(
                c,
                A_ub=A_ub, b_ub=b_ub,
                A_eq=A_eq, b_eq=b_eq,
                bounds=bounds,
                method=method,
                options={'maxiter': 100000}
            )
        except Exception as e:
            raise RuntimeError(f"LP solver failed: {e}")

        if not result.success:
            raise RuntimeError(f"LP optimization failed: {result.message}")

        # Extract solution
        objective_value = result.fun if minimize else -result.fun
        variables = {}
        for name, index in self._variables.items():
            variables[name] = result.x[index]

        return MapqnSolution(objective_value=objective_value, variables=variables)

    def solve_with_objective(self, c: np.ndarray, minimize: bool = True,
                             method: str = 'highs') -> MapqnSolution:
        """
        Solve with a custom objective coefficient vector.

        Args:
            c: Objective coefficient vector.
            minimize: If True, minimize; if False, maximize.
            method: LP solver method.

        Returns:
            MapqnSolution containing optimal value and variable values.
        """
        n_vars = self._variable_counter
        if n_vars == 0:
            return MapqnSolution(objective_value=0.0, variables={})

        if not minimize:
            c = -c

        # Separate equality and inequality constraints
        A_eq_list = []
        b_eq_list = []
        A_ub_list = []
        b_ub_list = []

        for constraint in self._constraints:
            coeffs = constraint.coefficients
            if len(coeffs) < n_vars:
                coeffs = np.concatenate([coeffs, np.zeros(n_vars - len(coeffs))])
            elif len(coeffs) > n_vars:
                coeffs = coeffs[:n_vars]

            if constraint.relationship == 'eq':
                A_eq_list.append(coeffs)
                b_eq_list.append(constraint.rhs)
            elif constraint.relationship == 'leq':
                A_ub_list.append(coeffs)
                b_ub_list.append(constraint.rhs)
            elif constraint.relationship == 'geq':
                A_ub_list.append(-coeffs)
                b_ub_list.append(-constraint.rhs)

        A_eq = np.array(A_eq_list) if A_eq_list else None
        b_eq = np.array(b_eq_list) if b_eq_list else None
        A_ub = np.array(A_ub_list) if A_ub_list else None
        b_ub = np.array(b_ub_list) if b_ub_list else None

        bounds = []
        for i in range(n_vars):
            if i in self._bounds:
                bounds.append(self._bounds[i])
            else:
                bounds.append((0.0, None))

        result = linprog(
            c,
            A_ub=A_ub, b_ub=b_ub,
            A_eq=A_eq, b_eq=b_eq,
            bounds=bounds,
            method=method,
            options={'maxiter': 100000}
        )

        if not result.success:
            raise RuntimeError(f"LP optimization failed: {result.message}")

        objective_value = result.fun if minimize else -result.fun
        variables = {}
        for name, index in self._variables.items():
            variables[name] = result.x[index]

        return MapqnSolution(objective_value=objective_value, variables=variables)


__all__ = ['MapqnLpModel', 'LinearConstraint', 'LinearConstraintBuilder']
