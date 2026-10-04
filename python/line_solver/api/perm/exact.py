"""
Exact matrix permanent algorithms.

The permanent of a matrix is similar to the determinant but uses only additions
(no subtractions). Computing the permanent is #P-complete, so exact computation
is expensive for large matrices. This module provides the exact algorithms:
inclusion-exclusion with multiplicities (twin of MATLAB perm.m and of the JAR
jline.lib.perm.Permanent), Ryser's formula in its Gray-code and naive forms
(jline.lib.perm.RyzerPermanent), and the O(n!) permutation enumeration
(jline.lib.perm.NaivePermanent).
"""

from itertools import combinations, permutations
from typing import Tuple, Optional

import numpy as np

from .base import PermResult, PermSolver


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


class Permanent(PermSolver):
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
        super().__init__(matrix)

        if solve:
            self.solve()

    def compute(self):
        """Compute the permanent."""
        self.value = _permanent_with_multiplicities(self.matrix)


class NaivePermanent(PermSolver):
    """
    Exact permanent by enumeration of all permutations.

    Port of jline.lib.perm.NaivePermanent. The permanent of an n-by-n matrix A
    is sum over all permutations p of prod_i A[i, p(i)], so this solver has
    O(n!) cost and is meant for small matrices and as a reference oracle.
    """

    def __init__(self, matrix: np.ndarray, solve: bool = False):
        """
        Initialize the naive permanent solver.

        Args:
            matrix: Matrix for which to compute the permanent
            solve: If True, compute immediately
        """
        super().__init__(matrix)

        if solve:
            self.solve()

    def compute(self):
        """Compute the exact permanent by permutation enumeration."""
        self.value = _permanent_naive(self.matrix)


class RyzerPermanent(PermSolver):
    """
    Exact permanent by Ryser's inclusion-exclusion formula.

    Port of jline.lib.perm.RyzerPermanent. Mode 'graycode' visits the column
    subsets in Gray-code order, so each step updates the row sums with a single
    column, giving O(2^n n) cost; any other mode selects the naive variant that
    regenerates the row sums for every subset, at O(2^n n^2) cost.
    """

    def __init__(self, matrix: np.ndarray, mode: str = "graycode", solve: bool = False):
        """
        Initialize Ryser's permanent solver.

        Args:
            matrix: Matrix for which to compute the permanent
            mode: 'graycode' for the Gray-code variant, anything else for naive
            solve: If True, compute immediately
        """
        super().__init__(matrix)
        self.mode = mode

        if solve:
            self.solve()

    def compute(self):
        """Compute the exact permanent with the selected Ryser variant."""
        if self.mode == "graycode":
            self.value = _ryser_graycode(self.matrix)
        else:
            self.value = _ryser_naive(self.matrix)


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

    # Ryser's formula, summed over all column subsets
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

    return perm


def _permanent_naive(matrix: np.ndarray) -> float:
    """
    Compute the permanent by enumerating all permutations.

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0

    rows = np.arange(n)
    value = 0.0
    for permutation in permutations(range(n)):
        value += float(np.prod(matrix[rows, list(permutation)]))
    return value


def _ryser_graycode(matrix: np.ndarray) -> float:
    """
    Compute the permanent with Ryser's formula in Gray-code order.

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0

    row_sum = np.zeros(n)
    current_bit = np.zeros(n, dtype=bool)
    value = 0.0

    for bit_index in _bit_to_modify(n):
        current_bit[bit_index] = not current_bit[bit_index]
        multiplier = 1.0 if current_bit[bit_index] else -1.0
        row_sum += multiplier * matrix[:, bit_index]

        nb_col = int(np.count_nonzero(current_bit))
        value += ((-1.0) ** nb_col) * float(np.prod(row_sum))

    return ((-1.0) ** n) * value


def _ryser_naive(matrix: np.ndarray) -> float:
    """
    Compute the permanent with Ryser's formula over explicit column subsets.

    Args:
        matrix: Input square matrix

    Returns:
        The permanent value
    """
    n = matrix.shape[0]
    if n == 0:
        return 1.0

    value = 0.0
    for i in range(n + 1):
        for combination in combinations(range(n), n - i):
            if combination:
                row_sums = matrix[:, list(combination)].sum(axis=1)
            else:
                row_sums = np.zeros(n)
            value += ((-1.0) ** (n - i)) * float(np.prod(row_sums))

    return ((-1.0) ** n) * value


def _bit_to_modify(m: int) -> list:
    """
    Build the sequence of bit positions to toggle for reflected Gray-code order.

    Args:
        m: Number of bits

    Returns:
        List of 2^m - 1 bit indices
    """
    if m == 1:
        return [0]
    sub = _bit_to_modify(m - 1)
    return sub + [m - 1] + sub


def _ryser_conditioning(matrix: np.ndarray) -> float:
    """
    Log10 magnitude of the largest Ryser term for this orientation.

    The inclusion-exclusion sum is largest when every column is selected, giving
    prod_i (sum_j a_ij). Since per(A) = per(A'), the two orientations share the
    same result but not the same cancellation, so this is the quantity to
    minimize when choosing between them. See _kb/03-api-layer.md.

    Args:
        matrix: Nonnegative input matrix

    Returns:
        Sum of log10 of the row sums, or +inf if any row sum is zero
    """
    line_sums = np.abs(matrix).sum(axis=1)
    if np.any(line_sums <= 0.0):
        return float("inf")
    return float(np.sum(np.log10(line_sums)))


def _find_unique_columns_with_multiplicities(
    matrix: np.ndarray
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Group repeated columns of the better-conditioned orientation of the matrix.

    Exploiting repeated rows means running the expansion on the transpose, which
    leaves the permanent unchanged but can raise the largest intermediate term
    by many orders of magnitude. Vandermonde-like matrices, such as the A_x of
    pfqn_lcfsqn_nc whose rows are geometric in the column index, lose every
    significant digit that way. Orientation is therefore chosen by conditioning
    first and grouping applied second, even when that forgoes the grouping.

    Args:
        matrix: Input matrix

    Returns:
        Tuple of (unique matrix, multiplicities array)
    """
    if _ryser_conditioning(matrix.T) < _ryser_conditioning(matrix):
        matrix = matrix.T

    n_rows, n_cols = matrix.shape

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
    'NaivePermanent',
    'RyzerPermanent',
]
