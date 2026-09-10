"""
Basic utility functions for KPC-Toolbox.

Native Python implementations of basic utilities including position finding,
logarithmic spacing, and spectral decomposition.
"""

import numpy as np
from numpy.linalg import eig, inv
from typing import Union, Tuple, List, NamedTuple


class SpectralDecomposition(NamedTuple):
    """Result of spectral decomposition."""
    spectrum: np.ndarray  # Eigenvalues
    projectors: List[np.ndarray]  # Projector matrices
    eigenvectors: np.ndarray  # Matrix V
    eigenvalue_matrix: np.ndarray  # Diagonal matrix D


def minpos(v: np.ndarray, n: int = None) -> Union[int, np.ndarray]:
    """
    Find the position(s) of minimum value(s) in a vector.

    Args:
        v: Input vector
        n: Number of smallest elements to find (optional)

    Returns:
        If n is None: Index (0-based) of the minimum value
        If n is given: Array of indices of the n smallest values
    """
    v = np.asarray(v).flatten()

    if len(v) == 0:
        raise ValueError("Vector cannot be empty")

    if n is None:
        return int(np.argmin(v))

    count = min(n, len(v))
    # Get indices sorted by value
    sorted_indices = np.argsort(v)
    return sorted_indices[:count]


def maxpos(v: np.ndarray, n: int = None) -> Union[int, np.ndarray]:
    """
    Find the position(s) of maximum value(s) in a vector.

    Args:
        v: Input vector
        n: Number of largest elements to find (optional)

    Returns:
        If n is None: Index (0-based) of the maximum value
        If n is given: Array of indices of the n largest values
    """
    v = np.asarray(v).flatten()

    if len(v) == 0:
        raise ValueError("Vector cannot be empty")

    if n is None:
        return int(np.argmax(v))

    count = min(n, len(v))
    # Get indices sorted by value descending
    sorted_indices = np.argsort(v)[::-1]
    return sorted_indices[:count]


def logspacei(a: float, b: float, points: int) -> np.ndarray:
    """
    Generate logarithmically spaced integers in [a, b].

    Args:
        a: Lower bound
        b: Upper bound
        points: Number of points

    Returns:
        Array of logarithmically spaced integers
    """
    if points <= 0:
        return np.array([], dtype=int)

    if points == 1:
        return np.array([int(round(a))], dtype=int)

    log_a = np.log10(a)
    log_b = np.log10(b)

    log_vals = np.linspace(log_a, log_b, points)
    result = np.round(np.power(10, log_vals)).astype(int)

    # Clamp to [a, b]
    result = np.clip(result, int(round(a)), int(round(b)))

    return result


def spectd(A: np.ndarray) -> SpectralDecomposition:
    """
    Compute the spectral decomposition of a matrix.

    Returns eigenvalues, eigenvectors, and projector matrices.

    Args:
        A: Input matrix (n x n)

    Returns:
        SpectralDecomposition containing:
            - spectrum: Array of eigenvalues
            - projectors: List of projector matrices P_k = V(:,k) * V^{-1}(k,:)
            - eigenvectors: Matrix V
            - eigenvalue_matrix: Diagonal matrix D of eigenvalues
    """
    A = np.asarray(A, dtype=np.float64)
    n = A.shape[0]

    # Compute eigendecomposition
    eigenvalues, V = eig(A)

    # Get real parts for real matrices
    eigenvalues = np.real(eigenvalues)
    V = np.real(V)

    # Compute inverse of V
    V_inv = inv(V)

    # Build diagonal eigenvalue matrix
    D = np.diag(eigenvalues)

    # Compute projectors: P_k = V(:,k) * V^{-1}(k,:)
    projectors = []
    for k in range(n):
        P_k = np.outer(V[:, k], V_inv[k, :])
        projectors.append(P_k)

    return SpectralDecomposition(
        spectrum=eigenvalues,
        projectors=projectors,
        eigenvectors=V,
        eigenvalue_matrix=D
    )


def ones(n: int, m: int = 1) -> np.ndarray:
    """
    Create a matrix of ones.

    Args:
        n: Number of rows
        m: Number of columns (default: 1 for column vector)

    Returns:
        Matrix of ones with shape (n, m)
    """
    return np.ones((n, m))


def e(n: int) -> np.ndarray:
    """
    Create a column vector of ones (used in matrix operations).

    Args:
        n: Size of the vector

    Returns:
        Column vector of ones with shape (n, 1)
    """
    return ones(n)


def eye(n: int) -> np.ndarray:
    """
    Create an identity matrix.

    Args:
        n: Size of the matrix

    Returns:
        Identity matrix of shape (n, n)
    """
    return np.eye(n)


def zeros(m: int, n: int) -> np.ndarray:
    """
    Create a zero matrix.

    Args:
        m: Number of rows
        n: Number of columns

    Returns:
        Zero matrix of shape (m, n)
    """
    return np.zeros((m, n))


__all__ = [
    'minpos',
    'maxpos',
    'logspacei',
    'spectd',
    'ones',
    'e',
    'eye',
    'zeros',
    'SpectralDecomposition',
]
