"""
Linear algebra operations for LINE native implementation.

Provides backend-agnostic linear algebra functions that dispatch
to NumPy/SciPy or PyTorch based on input type.
"""

from __future__ import annotations

import numpy as np
from scipy import linalg as scipy_linalg
from typing import Union, Tuple, Optional, TYPE_CHECKING

from .tensor import LineTensor, as_tensor
from .sparse import SparseMatrix
from .config import get_config, Backend

if TYPE_CHECKING:
    import torch

ArrayLike = Union[np.ndarray, LineTensor, list]


def solve(A: ArrayLike, b: ArrayLike) -> LineTensor:
    """
    Solve linear system Ax = b.

    Args:
        A: Coefficient matrix (n x n)
        b: Right-hand side (n,) or (n x m)

    Returns:
        Solution x
    """
    A = as_tensor(A)
    b = as_tensor(b)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.solve(A._data, b._data)
    else:
        result = scipy_linalg.solve(A.to_numpy(), b.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def lstsq(A: ArrayLike, b: ArrayLike) -> Tuple[LineTensor, ...]:
    """
    Solve least-squares problem min ||Ax - b||.

    Args:
        A: Coefficient matrix (m x n)
        b: Right-hand side (m,) or (m x k)

    Returns:
        Tuple of (solution, residuals, rank, singular_values)
    """
    A = as_tensor(A)
    b = as_tensor(b)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.lstsq(A._data, b._data)
        return (
            LineTensor(result.solution, backend='pytorch'),
            LineTensor(result.residuals if result.residuals.numel() > 0 else np.array([]), backend='pytorch'),
            int(result.rank),
            LineTensor(result.singular_values, backend='pytorch'),
        )
    else:
        x, residuals, rank, s = np.linalg.lstsq(A.to_numpy(), b.to_numpy(), rcond=None)
        return (
            LineTensor(x),
            LineTensor(residuals if len(residuals) > 0 else np.array([])),
            rank,
            LineTensor(s),
        )


def inv(A: ArrayLike) -> LineTensor:
    """
    Compute matrix inverse.

    Args:
        A: Square matrix (n x n)

    Returns:
        Inverse matrix A^(-1)
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.inv(A._data)
    else:
        result = scipy_linalg.inv(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def pinv(A: ArrayLike) -> LineTensor:
    """
    Compute Moore-Penrose pseudo-inverse.

    Args:
        A: Matrix (m x n)

    Returns:
        Pseudo-inverse A+
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.pinv(A._data)
    else:
        result = scipy_linalg.pinv(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def det(A: ArrayLike) -> float:
    """
    Compute matrix determinant.

    Args:
        A: Square matrix (n x n)

    Returns:
        Determinant value
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        return torch.linalg.det(A._data).item()
    else:
        return float(scipy_linalg.det(A.to_numpy()))


def norm(x: ArrayLike, ord: Optional[Union[int, float, str]] = None, axis: Optional[int] = None) -> Union[float, LineTensor]:
    """
    Compute vector or matrix norm.

    Args:
        x: Input tensor
        ord: Norm type (None=2-norm, 'fro', 1, 2, np.inf)
        axis: Axis for vector norm

    Returns:
        Norm value(s)
    """
    x = as_tensor(x)

    if x.backend == Backend.PYTORCH:
        import torch
        if ord == 'fro':
            ord = 'fro'
        result = torch.linalg.norm(x._data, ord=ord, dim=axis)
        if result.ndim == 0:
            return result.item()
        return LineTensor(result, backend='pytorch')
    else:
        result = np.linalg.norm(x.to_numpy(), ord=ord, axis=axis)
        if np.isscalar(result):
            return float(result)
        return LineTensor(result)


def eig(A: ArrayLike) -> Tuple[LineTensor, LineTensor]:
    """
    Compute eigenvalues and eigenvectors.

    Args:
        A: Square matrix (n x n)

    Returns:
        Tuple of (eigenvalues, eigenvectors)
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        eigenvalues, eigenvectors = torch.linalg.eig(A._data)
        return (
            LineTensor(eigenvalues, backend='pytorch'),
            LineTensor(eigenvectors, backend='pytorch'),
        )
    else:
        eigenvalues, eigenvectors = scipy_linalg.eig(A.to_numpy())
        return (
            LineTensor(eigenvalues),
            LineTensor(eigenvectors),
        )


def eigh(A: ArrayLike) -> Tuple[LineTensor, LineTensor]:
    """
    Compute eigenvalues and eigenvectors of Hermitian matrix.

    Args:
        A: Hermitian/symmetric matrix (n x n)

    Returns:
        Tuple of (eigenvalues, eigenvectors)
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        eigenvalues, eigenvectors = torch.linalg.eigh(A._data)
        return (
            LineTensor(eigenvalues, backend='pytorch'),
            LineTensor(eigenvectors, backend='pytorch'),
        )
    else:
        eigenvalues, eigenvectors = scipy_linalg.eigh(A.to_numpy())
        return (
            LineTensor(eigenvalues),
            LineTensor(eigenvectors),
        )


def eigvals(A: ArrayLike) -> LineTensor:
    """
    Compute eigenvalues only.

    Args:
        A: Square matrix (n x n)

    Returns:
        Eigenvalues
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.eigvals(A._data)
    else:
        result = scipy_linalg.eigvals(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def svd(A: ArrayLike, full_matrices: bool = True) -> Tuple[LineTensor, LineTensor, LineTensor]:
    """
    Compute Singular Value Decomposition.

    Args:
        A: Matrix (m x n)
        full_matrices: If True, compute full U and Vh

    Returns:
        Tuple of (U, S, Vh) where A = U @ diag(S) @ Vh
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        U, S, Vh = torch.linalg.svd(A._data, full_matrices=full_matrices)
        return (
            LineTensor(U, backend='pytorch'),
            LineTensor(S, backend='pytorch'),
            LineTensor(Vh, backend='pytorch'),
        )
    else:
        U, S, Vh = scipy_linalg.svd(A.to_numpy(), full_matrices=full_matrices)
        return (
            LineTensor(U),
            LineTensor(S),
            LineTensor(Vh),
        )


def qr(A: ArrayLike, mode: str = 'reduced') -> Tuple[LineTensor, LineTensor]:
    """
    Compute QR decomposition.

    Args:
        A: Matrix (m x n)
        mode: 'reduced' or 'complete'

    Returns:
        Tuple of (Q, R) where A = QR
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        Q, R = torch.linalg.qr(A._data, mode=mode)
        return (
            LineTensor(Q, backend='pytorch'),
            LineTensor(R, backend='pytorch'),
        )
    else:
        Q, R = scipy_linalg.qr(A.to_numpy(), mode='economic' if mode == 'reduced' else 'full')
        return (
            LineTensor(Q),
            LineTensor(R),
        )


def cholesky(A: ArrayLike, lower: bool = True) -> LineTensor:
    """
    Compute Cholesky decomposition.

    Args:
        A: Positive-definite matrix (n x n)
        lower: If True, return lower triangular L such that A = LL^T

    Returns:
        Cholesky factor
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.linalg.cholesky(A._data)
        if not lower:
            result = result.mH
    else:
        result = scipy_linalg.cholesky(A.to_numpy(), lower=lower)

    return LineTensor(result, backend=A.backend.value)


def lu(A: ArrayLike) -> Tuple[LineTensor, LineTensor, LineTensor]:
    """
    Compute LU decomposition with partial pivoting.

    Args:
        A: Square matrix (n x n)

    Returns:
        Tuple of (P, L, U) where PA = LU
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        P, L, U = torch.linalg.lu(A._data)
        return (
            LineTensor(P, backend='pytorch'),
            LineTensor(L, backend='pytorch'),
            LineTensor(U, backend='pytorch'),
        )
    else:
        P, L, U = scipy_linalg.lu(A.to_numpy())
        return (
            LineTensor(P),
            LineTensor(L),
            LineTensor(U),
        )


def matrix_rank(A: ArrayLike, tol: Optional[float] = None) -> int:
    """
    Compute matrix rank.

    Args:
        A: Input matrix
        tol: Tolerance for rank determination

    Returns:
        Matrix rank
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        return int(torch.linalg.matrix_rank(A._data, tol=tol))
    else:
        return int(np.linalg.matrix_rank(A.to_numpy(), tol=tol))


def cond(A: ArrayLike, p: Optional[Union[int, float, str]] = None) -> float:
    """
    Compute condition number.

    Args:
        A: Input matrix
        p: Norm type for condition number

    Returns:
        Condition number
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        import torch
        return torch.linalg.cond(A._data, p=p).item()
    else:
        return float(np.linalg.cond(A.to_numpy(), p=p))


def expm(A: ArrayLike) -> LineTensor:
    """
    Compute matrix exponential exp(A).

    Args:
        A: Square matrix (n x n)

    Returns:
        Matrix exponential
    """
    A = as_tensor(A)

    # PyTorch doesn't have native expm, use scipy
    result = scipy_linalg.expm(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def logm(A: ArrayLike) -> LineTensor:
    """
    Compute matrix logarithm log(A).

    Args:
        A: Square matrix (n x n)

    Returns:
        Matrix logarithm
    """
    A = as_tensor(A)

    # Use scipy for both backends
    result = scipy_linalg.logm(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def sqrtm(A: ArrayLike) -> LineTensor:
    """
    Compute matrix square root.

    Args:
        A: Square matrix (n x n)

    Returns:
        Matrix square root B such that B @ B = A
    """
    A = as_tensor(A)

    result = scipy_linalg.sqrtm(A.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def kron(A: ArrayLike, B: ArrayLike) -> LineTensor:
    """
    Compute Kronecker product.

    Args:
        A: First matrix
        B: Second matrix

    Returns:
        Kronecker product A \otimes B
    """
    A = as_tensor(A)
    B = as_tensor(B)

    if A.backend == Backend.PYTORCH:
        import torch
        result = torch.kron(A._data, B._data)
    else:
        result = np.kron(A.to_numpy(), B.to_numpy())

    return LineTensor(result, backend=A.backend.value)


def trace(A: ArrayLike) -> float:
    """
    Compute matrix trace.

    Args:
        A: Square matrix

    Returns:
        Sum of diagonal elements
    """
    A = as_tensor(A)

    if A.backend == Backend.PYTORCH:
        return A._data.trace().item()
    else:
        return float(np.trace(A.to_numpy()))


def diag(v: ArrayLike, k: int = 0) -> LineTensor:
    """
    Create diagonal matrix or extract diagonal.

    Args:
        v: Vector (for creating diagonal) or matrix (for extracting)
        k: Diagonal offset

    Returns:
        Diagonal matrix or diagonal vector
    """
    v = as_tensor(v)

    if v.backend == Backend.PYTORCH:
        import torch
        result = torch.diag(v._data, diagonal=k)
    else:
        result = np.diag(v.to_numpy(), k=k)

    return LineTensor(result, backend=v.backend.value)


def outer(a: ArrayLike, b: ArrayLike) -> LineTensor:
    """
    Compute outer product.

    Args:
        a: First vector
        b: Second vector

    Returns:
        Outer product a \otimes b
    """
    a = as_tensor(a)
    b = as_tensor(b)

    if a.backend == Backend.PYTORCH:
        import torch
        result = torch.outer(a._data.flatten(), b._data.flatten())
    else:
        result = np.outer(a.to_numpy().flatten(), b.to_numpy().flatten())

    return LineTensor(result, backend=a.backend.value)
