"""
MAPQN Linear Programming Model.

Base class for representing MAP queueing network linear programming models.
Provides the foundation for LP-based optimization methods in MAPQN analysis,
including constraint formulation and objective function definition.

Constraints are stored SPARSELY as (column, value) pairs and assembled into a
single scipy.sparse matrix at solve time. A dense row per constraint is not
viable at the size of the reference QRF models: the paper's BAS instance
(M=5, N=10, K=2, MR=5) has 60510 columns and ~114500 rows, i.e. 55 GB dense
against ~630k nonzeros.

Duplicate (row, column) entries are SUMMED during assembly, which reproduces
the accumulating semantics of ``coefficients[index] += coefficient``.

Uses scipy.optimize.linprog (HiGHS) as the LP solver backend.
"""

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix
from typing import Dict, List, Optional, Tuple

from .solution import MapqnSolution


class LinearConstraint:
    """
    Represents a single linear constraint, stored sparsely.

    The row is held as parallel index/value arrays; repeated indices are
    permitted and are summed when the model is assembled.

    Attributes:
        indices: Column indices of the nonzero terms.
        values: Coefficients aligned with ``indices``.
        relationship: One of 'eq' (=), 'leq' (<=), or 'geq' (>=).
        rhs: Right-hand side value.
    """

    __slots__ = ('indices', 'values', 'relationship', 'rhs')

    def __init__(self, coefficients=None, relationship: str = 'eq', rhs: float = 0.0,
                 indices: Optional[np.ndarray] = None,
                 values: Optional[np.ndarray] = None):
        """
        Args:
            coefficients: Dense coefficient vector (legacy form). Ignored when
                ``indices``/``values`` are supplied.
            relationship: 'eq', 'leq' or 'geq'.
            rhs: Right-hand side.
            indices: Sparse column indices.
            values: Sparse coefficients.
        """
        if indices is not None:
            self.indices = np.asarray(indices, dtype=np.intp)
            self.values = np.asarray(values, dtype=np.float64)
        elif coefficients is not None:
            dense = np.asarray(coefficients, dtype=np.float64)
            nz = np.flatnonzero(dense)
            self.indices = nz.astype(np.intp)
            self.values = dense[nz]
        else:
            self.indices = np.empty(0, dtype=np.intp)
            self.values = np.empty(0, dtype=np.float64)
        self.relationship = relationship
        self.rhs = float(rhs)

    @property
    def coefficients(self) -> np.ndarray:
        """Dense view of the row (materialized on demand, legacy accessor)."""
        n = int(self.indices.max()) + 1 if self.indices.size else 0
        dense = np.zeros(n)
        np.add.at(dense, self.indices, self.values)
        return dense

    def nnz(self) -> int:
        """Number of stored terms (before duplicate summation)."""
        return int(self.indices.size)

    def __repr__(self) -> str:
        return ('LinearConstraint(nnz=%d, relationship=%r, rhs=%g)'
                % (self.indices.size, self.relationship, self.rhs))


class LinearConstraintBuilder:
    """
    Helper class for building linear constraints incrementally.

    Terms are appended to a triplet list rather than scattered into a dense
    vector, so building a row costs O(terms) in time and memory instead of
    O(number of variables).

    Example:
        >>> builder = model.constraint_builder()
        >>> builder.add_term('U_1', 1.0).add_term('U_2', -1.0).eq(0.0)
    """

    def __init__(self, model: 'MapqnLpModel'):
        self.model = model
        self._indices: List[int] = []
        self._values: List[float] = []

    def add_term(self, var_name: str, coefficient: float) -> 'LinearConstraintBuilder':
        """Add a term to the constraint."""
        index = self.model.get_variable_index(var_name)
        if index is not None:
            self._indices.append(index)
            self._values.append(coefficient)
        return self

    def add_term_by_index(self, var_index: int, coefficient: float) -> 'LinearConstraintBuilder':
        """Add a term by variable index."""
        if 0 <= var_index < self.model.get_num_variables():
            self._indices.append(var_index)
            self._values.append(coefficient)
        return self

    def _freeze(self, relationship: str, rhs: float) -> LinearConstraint:
        return LinearConstraint(
            relationship=relationship, rhs=rhs,
            indices=np.array(self._indices, dtype=np.intp),
            values=np.array(self._values, dtype=np.float64))

    def eq(self, rhs: float) -> LinearConstraint:
        """Create equality constraint (=)."""
        return self._freeze('eq', rhs)

    def leq(self, rhs: float) -> LinearConstraint:
        """Create less-than-or-equal constraint (<=)."""
        return self._freeze('leq', rhs)

    def geq(self, rhs: float) -> LinearConstraint:
        """Create greater-than-or-equal constraint (>=)."""
        return self._freeze('geq', rhs)


class MapqnLpModel:
    """
    Base class for MAPQN Linear Programming models.

    Provides methods for:
    - Variable registration and tracking
    - Constraint addition
    - Objective function creation
    - LP solving via scipy (HiGHS, sparse constraint matrices)

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

    def set_variable_bounds(self, name: str, lb: Optional[float],
                            ub: Optional[float]) -> None:
        """
        Reset the bounds of an already-registered variable.

        Used to pin structurally-zero variables (lb = ub = 0) rather than
        emitting one equality row each, which keeps the constraint matrix
        smaller and better conditioned.

        Args:
            name: Variable name.
            lb: Lower bound (None for unbounded).
            ub: Upper bound (None for unbounded).

        Raises:
            KeyError: If the variable was never registered.
        """
        index = self._variables.get(name)
        if index is None:
            raise KeyError("unknown variable %r" % (name,))
        self._bounds[index] = (lb, ub)

    def get_variable_index(self, name: str) -> Optional[int]:
        """Get variable index by name, or None if not found."""
        return self._variables.get(name)

    def get_num_variables(self) -> int:
        """Get total number of registered variables."""
        return self._variable_counter

    def add_constraint(self, constraint: LinearConstraint) -> None:
        """Add a constraint to the model."""
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

    def get_matrix_stats(self) -> Dict[str, int]:
        """Report assembled problem size (rows, columns, stored nonzeros)."""
        A_ub, b_ub, A_eq, b_eq = self._assemble()
        nnz = 0
        rows = 0
        for A in (A_ub, A_eq):
            if A is not None:
                nnz += int(A.nnz)
                rows += int(A.shape[0])
        return {'rows': rows, 'columns': self._variable_counter, 'nonzeros': nnz,
                'rows_ub': 0 if A_ub is None else int(A_ub.shape[0]),
                'rows_eq': 0 if A_eq is None else int(A_eq.shape[0])}

    def _assemble(self):
        """
        Assemble the sparse constraint matrices.

        Returns:
            Tuple (A_ub, b_ub, A_eq, b_eq); the matrices are scipy CSR or None
            when the corresponding family is empty. '>=' rows are negated into
            '<=' form. Duplicate (row, column) entries are summed by the
            COO -> CSR conversion.
        """
        n_vars = self._variable_counter

        eq_rows: List[np.ndarray] = []
        eq_cols: List[np.ndarray] = []
        eq_vals: List[np.ndarray] = []
        b_eq: List[float] = []
        ub_rows: List[np.ndarray] = []
        ub_cols: List[np.ndarray] = []
        ub_vals: List[np.ndarray] = []
        b_ub: List[float] = []

        for constraint in self._constraints:
            idx = constraint.indices
            val = constraint.values
            if idx.size and idx.max() >= n_vars:
                keep = idx < n_vars
                idx = idx[keep]
                val = val[keep]

            if constraint.relationship == 'eq':
                row = len(b_eq)
                eq_rows.append(np.full(idx.size, row, dtype=np.intp))
                eq_cols.append(idx)
                eq_vals.append(val)
                b_eq.append(constraint.rhs)
            elif constraint.relationship == 'leq':
                row = len(b_ub)
                ub_rows.append(np.full(idx.size, row, dtype=np.intp))
                ub_cols.append(idx)
                ub_vals.append(val)
                b_ub.append(constraint.rhs)
            elif constraint.relationship == 'geq':
                # Convert >= to <= by negating both sides
                row = len(b_ub)
                ub_rows.append(np.full(idx.size, row, dtype=np.intp))
                ub_cols.append(idx)
                ub_vals.append(-val)
                b_ub.append(-constraint.rhs)
            else:
                raise ValueError("Unknown constraint relationship: %r"
                                 % (constraint.relationship,))

        def build(rows, cols, vals, b):
            if not b:
                return None, None
            n_rows = len(b)
            if rows:
                r = np.concatenate(rows)
                c = np.concatenate(cols)
                v = np.concatenate(vals)
            else:
                r = np.empty(0, dtype=np.intp)
                c = np.empty(0, dtype=np.intp)
                v = np.empty(0, dtype=np.float64)
            # coo -> csr sums duplicate (i,j) entries, matching row(idx) += v
            A = coo_matrix((v, (r, c)), shape=(n_rows, n_vars)).tocsr()
            return A, np.asarray(b, dtype=np.float64)

        A_eq, b_eq_arr = build(eq_rows, eq_cols, eq_vals, b_eq)
        A_ub, b_ub_arr = build(ub_rows, ub_cols, ub_vals, b_ub)
        return A_ub, b_ub_arr, A_eq, b_eq_arr

    def _bounds_list(self) -> List[Tuple[Optional[float], Optional[float]]]:
        bounds = []
        for i in range(self._variable_counter):
            if i in self._bounds:
                bounds.append(self._bounds[i])
            else:
                bounds.append((0.0, None))  # Default: non-negative
        return bounds

    def _run(self, c: np.ndarray, minimize: bool, method: str) -> MapqnSolution:
        if not minimize:
            c = -c

        A_ub, b_ub, A_eq, b_eq = self._assemble()

        try:
            result = linprog(
                c,
                A_ub=A_ub, b_ub=b_ub,
                A_eq=A_eq, b_eq=b_eq,
                bounds=self._bounds_list(),
                method=method,
                options={'maxiter': 100000},
            )
        except Exception as e:
            raise RuntimeError("LP solver failed: %s" % (e,))

        if not result.success:
            raise RuntimeError("LP optimization failed: %s" % (result.message,))

        objective_value = result.fun if minimize else -result.fun
        x = result.x
        variables = {}
        for name, index in self._variables.items():
            variables[name] = x[index]

        return MapqnSolution(objective_value=objective_value, variables=variables)

    def solve(self, objective_var: str, minimize: bool = True,
              method: str = 'highs') -> MapqnSolution:
        """
        Solve the LP model.

        Args:
            objective_var: Variable name to optimize.
            minimize: If True, minimize; if False, maximize.
            method: LP solver method ('highs', 'highs-ds', 'highs-ipm').

        Returns:
            MapqnSolution containing optimal value and variable values.

        Raises:
            RuntimeError: If LP is infeasible or unbounded.
        """
        if self._variable_counter == 0:
            return MapqnSolution(objective_value=0.0, variables={})

        return self._run(self.create_objective_coefficients(objective_var),
                         minimize, method)

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
        if self._variable_counter == 0:
            return MapqnSolution(objective_value=0.0, variables={})

        return self._run(np.asarray(c, dtype=np.float64), minimize, method)


__all__ = ['MapqnLpModel', 'LinearConstraint', 'LinearConstraintBuilder']
