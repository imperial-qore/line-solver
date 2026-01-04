"""
Perm: Matrix Permanent Computation.

Native Python implementations of matrix permanent algorithms.

The permanent of a matrix is similar to the determinant but uses only additions
(no subtractions). Computing the permanent is #P-complete, so exact computation
is expensive for large matrices. This module provides algorithms that exploit
structure (repeated rows/columns) for computational savings.
"""

import numpy as np
from numpy.typing import NDArray
from typing import Tuple, List, Optional
from dataclasses import dataclass
import time


@dataclass
class PermResult:
    """Result of permanent computation."""
    value: float
    time_ms: float
    n: int


def compute_permanent(matrix: np.ndarray, use_multiplicities: bool = True) -> float:
    """
    Compute the permanent of a matrix.

    Uses the inclusion-exclusion principle with multiplicities to efficiently
    handle matrices with duplicate rows or columns.

    Args:
        matrix: Input square matrix
        use_multiplicities: If True, exploit repeated rows/columns for efficiency

    Returns:
        The permanent value
    """
    matrix = np.asarray(matrix, dtype=float)
    n = matrix.shape[0]

    if n == 0:
        return 1.0
    if n == 1:
        return matrix[0, 0]
    if n == 2:
        return matrix[0, 0] * matrix[1, 1] + matrix[0, 1] * matrix[1, 0]

    if use_multiplicities:
        return _permanent_with_multiplicities(matrix)
    else:
        return _permanent_ryser(matrix)


def permanent(matrix: np.ndarray) -> float:
    """
    Compute the permanent of a matrix (alias for compute_permanent).

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    return compute_permanent(matrix)


def perm(matrix: np.ndarray) -> float:
    """
    Compute the permanent of a matrix (alias for compute_permanent).

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    return compute_permanent(matrix)


class Permanent:
    """
    Permanent computation solver class.

    Mirrors the MATLAB perm.m function which computes the permanent
    of a matrix by applying computational savings when rows or columns are repeated.
    """

    def __init__(self, matrix: np.ndarray, solve: bool = False):
        """
        Initialize permanent solver.

        Args:
            matrix: Matrix for which to compute the permanent
            solve: If True, compute immediately
        """
        self.matrix = np.asarray(matrix, dtype=float)
        self.n = self.matrix.shape[0]
        self.value = 0.0
        self.time_ms = 0.0

        if solve:
            self.solve()

    def solve(self):
        """Compute the permanent and measure time."""
        start_time = time.time()
        self.compute()
        self.time_ms = (time.time() - start_time) * 1000

    def compute(self):
        """Compute the permanent."""
        self.value = _permanent_with_multiplicities(self.matrix)

    def get_result(self) -> PermResult:
        """Get the computation result."""
        return PermResult(value=self.value, time_ms=self.time_ms, n=self.n)


def _permanent_with_multiplicities(matrix: np.ndarray) -> float:
    """
    Compute permanent using inclusion-exclusion with multiplicities.

    This algorithm detects and exploits repeated columns/rows for efficiency.

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0

    # Find unique columns and their multiplicities
    unique_matrix, multiplicities = _find_unique_columns_with_multiplicities(matrix)

    R = len(multiplicities)
    value = 0.0

    # Initialize iterator
    f = np.zeros(R, dtype=int)

    while True:
        # Compute term
        term = (-1.0) ** np.sum(f)

        # Multinomial coefficients
        for j in range(R):
            term *= _binomial_coefficient(multiplicities[j], f[j])

        # Product term
        for i in range(n):
            sum_term = 0.0
            for k in range(R):
                sum_term += f[k] * unique_matrix[i, k]
            term *= sum_term

        value += term

        # Get next iteration
        f = _pprod_next(f, multiplicities)
        if f is None:
            break

    return ((-1.0) ** n) * value


def _permanent_ryser(matrix: np.ndarray) -> float:
    """
    Compute permanent using Ryser's formula.

    This is an O(2^n * n) algorithm that doesn't exploit repeated rows/columns
    but is more straightforward.

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0

    # Ryser's formula with Gray code
    perm = 0.0
    for subset in range(1 << n):
        sign = (-1) ** (n - bin(subset).count('1'))
        prod = 1.0
        for i in range(n):
            row_sum = 0.0
            for j in range(n):
                if subset & (1 << j):
                    row_sum += matrix[i, j]
            prod *= row_sum
        perm += sign * prod

    return ((-1) ** n) * perm


def _find_unique_columns_with_multiplicities(
    matrix: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Find unique columns in the matrix and their multiplicities.

    If no repeated columns, checks for repeated rows instead.

    Args:
        matrix: Input matrix

    Returns:
        Tuple of (unique matrix, multiplicities array)
    """
    n_rows, n_cols = matrix.shape

    # Find unique columns
    column_map = {}
    unique_columns = []
    multiplicities = []

    for j in range(n_cols):
        column = tuple(matrix[:, j])
        if column in column_map:
            idx = column_map[column]
            multiplicities[idx] += 1
        else:
            column_map[column] = len(unique_columns)
            unique_columns.append(column)
            multiplicities.append(1)

    # If no repeated columns, check for repeated rows
    if all(m == 1 for m in multiplicities):
        row_map = {}
        unique_rows = []
        row_multiplicities = []

        for i in range(n_rows):
            row = tuple(matrix[i, :])
            if row in row_map:
                idx = row_map[row]
                row_multiplicities[idx] += 1
            else:
                row_map[row] = len(unique_rows)
                unique_rows.append(row)
                row_multiplicities.append(1)

        # If there are repeated rows, transpose and use row multiplicities
        if any(m > 1 for m in row_multiplicities):
            unique_matrix = np.array(unique_rows).T
            return unique_matrix, np.array(row_multiplicities, dtype=int)

    # Use unique columns
    unique_matrix = np.column_stack(unique_columns) if unique_columns else np.zeros((n_rows, 0))
    return unique_matrix, np.array(multiplicities, dtype=int)


def _binomial_coefficient(n: int, k: int) -> float:
    """
    Compute binomial coefficient C(n, k).

    Args:
        n: Upper parameter
        k: Lower parameter

    Returns:
        Binomial coefficient
    """
    if k > n or k < 0:
        return 0.0
    if k == 0 or k == n:
        return 1.0

    result = 1.0
    for i in range(1, min(k, n - k) + 1):
        result = result * (n - i + 1) / i
    return result


def _pprod_next(current: np.ndarray, bounds: np.ndarray) -> Optional[np.ndarray]:
    """
    MATLAB pprod iterator - generates the next state in the sequence.

    Args:
        current: Current state vector
        bounds: Upper bounds vector

    Returns:
        Next state vector, or None if sequence is complete
    """
    n = current.copy()
    R = len(bounds)

    # Check if we've reached the maximum state
    if np.all(n == bounds):
        return None

    s = R - 1
    while s >= 0 and n[s] == bounds[s]:
        n[s] = 0
        s -= 1

    if s < 0:
        return None

    n[s] += 1
    return n


__all__ = [
    'compute_permanent',
    'permanent',
    'perm',
    'Permanent',
    'PermResult',
]
