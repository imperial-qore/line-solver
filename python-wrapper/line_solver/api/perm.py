"""
Permanent matrix computation utilities.

This module provides Python wrappers for computing matrix permanents
using various algorithms available in the JAR backend.

Functions:
    perm            - Compute exact permanent with optimizations
    perm_heuristic  - Compute approximate permanent using heuristic method
    perm_ryzer      - Compute permanent using Ryzer algorithm
    perm_bethe      - Compute permanent using Bethe approximation
    perm_naive      - Compute permanent using naive algorithm
    perm_validate   - Validate matrix for permanent computation

Example:
    >>> import numpy as np
    >>> from line_solver.api.perm import perm, perm_heuristic
    >>> A = np.array([[1, 2], [3, 4]])
    >>> perm_val = perm(A)
    >>> perm_approx_val = perm_heuristic(A)
"""

import jpype
import numpy as np
from line_solver import jlineMatrixFromArray


def perm(A):
    """
    Compute the exact permanent of matrix A using optimized method.

    The computation automatically detects and exploits repeated rows or columns
    for efficiency.

    Args:
        A (array-like): Input matrix (n x n) or (n x r) with r unique columns.
                       All elements must be non-negative.

    Returns:
        float: The permanent value.

    Raises:
        ValueError: If matrix contains negative values or non-numeric types.
        RuntimeError: If JAR computation fails.

    Example:
        >>> A = np.array([[1, 0.5], [0.5, 1]])
        >>> perm_val = perm(A)
    """
    if not _validate_matrix(A):
        raise ValueError("Matrix must be numeric with non-negative elements")

    A_jline = jlineMatrixFromArray(A)

    try:
        # Access Permanent class from JAR
        Permanent = jpype.JPackage('jline').lib.perm.Permanent
        result = Permanent.exact(A_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"JAR permanent computation failed: {str(e)}")


def perm_heuristic(A):
    """
    Compute approximate permanent using heuristic method (Sinkhorn scaling).

    This method is efficient for large matrices where exact computation
    would be computationally prohibitive.

    Args:
        A (array-like): Input matrix with positive elements.
                       All elements must be > 0.

    Returns:
        float: Approximate permanent value.

    Raises:
        ValueError: If matrix contains non-positive values.
        RuntimeError: If JAR computation fails.

    Note:
        The heuristic method uses Sinkhorn scaling to make the matrix
        approximately doubly stochastic, then applies mean-field approximation.

    Example:
        >>> A = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        >>> perm_approx = heuristic(A)
    """
    if not _validate_positive_matrix(A):
        raise ValueError("Matrix must have strictly positive elements")

    A_jline = jlineMatrixFromArray(A)

    try:
        # Access Permanent class from JAR
        Permanent = jpype.JPackage('jline').lib.perm.Permanent
        result = Permanent.heuristic(A_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"JAR heuristic permanent computation failed: {str(e)}")


def perm_ryzer(A):
    """
    Compute permanent using Ryzer's algorithm.

    Ryzer's algorithm is efficient for moderately-sized matrices.

    Args:
        A (array-like): Input matrix (n x n).

    Returns:
        float: Permanent value computed using Ryzer algorithm.

    Raises:
        ValueError: If matrix is not square.
        RuntimeError: If JAR computation fails.

    Example:
        >>> A = np.array([[1, 2], [3, 4]])
        >>> perm = ryzer(A)
    """
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("Matrix must be square for Ryzer algorithm")

    A_jline = jlineMatrixFromArray(A)

    try:
        RyzerPermanent = jpype.JPackage('jline').lib.perm.RyzerPermanent
        result = RyzerPermanent.computePermanent(A_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"JAR Ryzer computation failed: {str(e)}")


def perm_bethe(A):
    """
    Compute permanent using Bethe approximation.

    The Bethe approximation provides good accuracy with moderate computation cost.

    Args:
        A (array-like): Input matrix (n x n).

    Returns:
        float: Approximate permanent value.

    Raises:
        ValueError: If matrix is not square.
        RuntimeError: If JAR computation fails.

    Example:
        >>> A = np.array([[1, 0.5], [0.5, 1]])
        >>> perm = bethe(A)
    """
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("Matrix must be square for Bethe approximation")

    A_jline = jlineMatrixFromArray(A)

    try:
        BethePermanent = jpype.JPackage('jline').lib.perm.BethePermanent
        result = BethePermanent.computePermanent(A_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"JAR Bethe computation failed: {str(e)}")


def perm_naive(A):
    """
    Compute permanent using naive algorithm.

    The naive algorithm provides straightforward implementation for small matrices.

    Args:
        A (array-like): Input matrix (n x n).

    Returns:
        float: Permanent value.

    Raises:
        ValueError: If matrix is not square.
        RuntimeError: If JAR computation fails.

    Example:
        >>> A = np.array([[1, 2], [3, 4]])
        >>> perm = naive(A)
    """
    A = np.asarray(A)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("Matrix must be square for heuristic algorithm")

    A_jline = jlineMatrixFromArray(A)

    try:
        NaivePermanent = jpype.JPackage('jline').lib.perm.NaivePermanent
        result = NaivePermanent.computePermanent(A_jline)
        return float(result)
    except Exception as e:
        raise RuntimeError(f"JAR naive algorithm computation failed: {str(e)}")


def perm_validate(A, allow_zero=True, allow_negative=False):
    """
    Validate a matrix for permanent computation.

    Args:
        A (array-like): Matrix to validate.
        allow_zero (bool): Whether to allow zero elements. Default True.
        allow_negative (bool): Whether to allow negative elements. Default False.

    Returns:
        bool: True if valid, False otherwise.

    Example:
        >>> A = np.array([[1, 2], [3, 4]])
        >>> is_valid = validate(A)
    """
    try:
        A = np.asarray(A)

        # Check numeric type
        if not np.issubdtype(A.dtype, np.number):
            return False

        # Check for NaN or Inf
        if np.any(np.isnan(A)) or np.any(np.isinf(A)):
            return False

        # Check for negative values
        if allow_negative is False and np.any(A < 0):
            return False

        # Check for zero values
        if allow_zero is False and np.any(A == 0):
            return False

        return True
    except (TypeError, ValueError):
        return False



# Helper functions
def _validate_matrix(A):
    """Check if matrix is numeric with non-negative elements."""
    try:
        A = np.asarray(A)
        if not np.issubdtype(A.dtype, np.number):
            return False
        if np.any(np.isnan(A)) or np.any(np.isinf(A)):
            return False
        if np.any(A < 0):
            return False
        return True
    except (TypeError, ValueError):
        return False


def _validate_positive_matrix(A):
    """Check if matrix has strictly positive elements."""
    try:
        A = np.asarray(A)
        if not np.issubdtype(A.dtype, np.number):
            return False
        if np.any(np.isnan(A)) or np.any(np.isinf(A)):
            return False
        if np.any(A <= 0):
            return False
        return True
    except (TypeError, ValueError):
        return False
