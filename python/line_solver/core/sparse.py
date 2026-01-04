"""
Sparse matrix operations for LINE native implementation.

Provides SparseMatrix class with CSC format support, compatible with
scipy.sparse and optimized for Markov chain computations.
"""

from __future__ import annotations

import numpy as np
from scipy import sparse
from scipy.sparse import csc_matrix, csr_matrix, coo_matrix
from scipy.sparse.linalg import spsolve, splu, eigs, gmres
from typing import Union, Tuple, Optional, TYPE_CHECKING

from .config import get_config, Backend

if TYPE_CHECKING:
    import torch

ArrayLike = Union[np.ndarray, list, tuple]


class SparseMatrix:
    """
    Sparse matrix with CSC storage format.

    Compatible with scipy.sparse and optimized for queueing network
    computations (CTMC, QBD processes).

    Attributes:
        shape: Matrix dimensions (rows, cols)
        nnz: Number of non-zero elements
        format: Storage format ('csc', 'csr', 'coo')
    """

    __slots__ = ('_mat', '_format')

    def __init__(
        self,
        data: Union[np.ndarray, sparse.spmatrix, 'SparseMatrix', Tuple],
        shape: Optional[Tuple[int, int]] = None,
        format: str = 'csc',
    ):
        """
        Create a sparse matrix.

        Args:
            data: Input data. Can be:
                - Dense numpy array
                - scipy.sparse matrix
                - Tuple (data, indices, indptr) for CSC/CSR
                - Tuple (data, (row, col)) for COO
            shape: Matrix shape (required for COO format)
            format: Storage format ('csc', 'csr', 'coo')
        """
        self._format = format

        if isinstance(data, SparseMatrix):
            self._mat = data._mat.copy()
        elif isinstance(data, sparse.spmatrix):
            self._mat = self._convert_format(data, format)
        elif isinstance(data, np.ndarray):
            # Dense array - convert to sparse
            if format == 'csc':
                self._mat = csc_matrix(data)
            elif format == 'csr':
                self._mat = csr_matrix(data)
            else:
                self._mat = coo_matrix(data)
        elif isinstance(data, tuple):
            # COO or CSC/CSR format
            if len(data) == 2 and isinstance(data[1], tuple):
                # COO: (data, (row, col))
                self._mat = coo_matrix(data, shape=shape)
                self._mat = self._convert_format(self._mat, format)
            elif len(data) == 3:
                # CSC/CSR: (data, indices, indptr)
                if format == 'csc':
                    self._mat = csc_matrix(data, shape=shape)
                else:
                    self._mat = csr_matrix(data, shape=shape)
            else:
                raise ValueError(f"Invalid tuple format for sparse matrix")
        else:
            raise TypeError(f"Cannot create SparseMatrix from {type(data)}")

    @staticmethod
    def _convert_format(mat: sparse.spmatrix, format: str) -> sparse.spmatrix:
        """Convert sparse matrix to specified format."""
        if format == 'csc':
            return mat.tocsc()
        elif format == 'csr':
            return mat.tocsr()
        elif format == 'coo':
            return mat.tocoo()
        else:
            raise ValueError(f"Unknown format: {format}")

    # ========== Properties ==========

    @property
    def shape(self) -> Tuple[int, int]:
        """Matrix shape."""
        return self._mat.shape

    @property
    def nnz(self) -> int:
        """Number of non-zero elements."""
        return self._mat.nnz

    @property
    def format(self) -> str:
        """Storage format."""
        return self._format

    @property
    def data(self) -> np.ndarray:
        """Non-zero values array."""
        return self._mat.data

    @property
    def dtype(self) -> np.dtype:
        """Data type."""
        return self._mat.dtype

    # ========== CSC/CSR specific ==========

    @property
    def indices(self) -> np.ndarray:
        """Row indices (CSC) or column indices (CSR)."""
        return self._mat.indices

    @property
    def indptr(self) -> np.ndarray:
        """Index pointer array."""
        return self._mat.indptr

    # ========== Conversion ==========

    def to_dense(self) -> np.ndarray:
        """Convert to dense numpy array."""
        return self._mat.toarray()

    def to_scipy(self) -> sparse.spmatrix:
        """Get underlying scipy.sparse matrix."""
        return self._mat

    def tocsc(self) -> 'SparseMatrix':
        """Convert to CSC format."""
        if self._format == 'csc':
            return self
        result = SparseMatrix.__new__(SparseMatrix)
        result._mat = self._mat.tocsc()
        result._format = 'csc'
        return result

    def tocsr(self) -> 'SparseMatrix':
        """Convert to CSR format."""
        if self._format == 'csr':
            return self
        result = SparseMatrix.__new__(SparseMatrix)
        result._mat = self._mat.tocsr()
        result._format = 'csr'
        return result

    def tocoo(self) -> 'SparseMatrix':
        """Convert to COO format."""
        if self._format == 'coo':
            return self
        result = SparseMatrix.__new__(SparseMatrix)
        result._mat = self._mat.tocoo()
        result._format = 'coo'
        return result

    def to_torch_sparse(self, device: str = 'cpu') -> 'torch.Tensor':
        """
        Convert to PyTorch sparse tensor.

        Args:
            device: Target device

        Returns:
            PyTorch sparse COO tensor
        """
        import torch

        coo = self._mat.tocoo()
        indices = torch.LongTensor([coo.row, coo.col])
        values = torch.FloatTensor(coo.data)

        return torch.sparse_coo_tensor(
            indices, values, coo.shape, device=device
        ).coalesce()

    def copy(self) -> 'SparseMatrix':
        """Create a copy."""
        result = SparseMatrix.__new__(SparseMatrix)
        result._mat = self._mat.copy()
        result._format = self._format
        return result

    # ========== Factory Methods ==========

    @staticmethod
    def from_dense(arr: np.ndarray, format: str = 'csc') -> 'SparseMatrix':
        """Create from dense array."""
        return SparseMatrix(arr, format=format)

    @staticmethod
    def eye(n: int, format: str = 'csc') -> 'SparseMatrix':
        """Create sparse identity matrix."""
        return SparseMatrix(sparse.eye(n, format=format))

    @staticmethod
    def zeros(shape: Tuple[int, int], format: str = 'csc') -> 'SparseMatrix':
        """Create sparse zero matrix."""
        if format == 'csc':
            return SparseMatrix(csc_matrix(shape, dtype=np.float64))
        elif format == 'csr':
            return SparseMatrix(csr_matrix(shape, dtype=np.float64))
        else:
            return SparseMatrix(coo_matrix(shape, dtype=np.float64))

    @staticmethod
    def diags(
        diagonals: list,
        offsets: Union[int, list] = 0,
        shape: Optional[Tuple[int, int]] = None,
        format: str = 'csc',
    ) -> 'SparseMatrix':
        """Create sparse diagonal matrix."""
        mat = sparse.diags(diagonals, offsets, shape=shape, format=format)
        return SparseMatrix(mat)

    # ========== Arithmetic Operations ==========

    def __add__(self, other: Union['SparseMatrix', np.ndarray, float]) -> 'SparseMatrix':
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = self._format

        if isinstance(other, SparseMatrix):
            result._mat = self._mat + other._mat
        elif isinstance(other, np.ndarray):
            result._mat = self._mat + csc_matrix(other)
        else:
            result._mat = self._mat + other

        return result

    def __radd__(self, other):
        return self.__add__(other)

    def __sub__(self, other: Union['SparseMatrix', np.ndarray, float]) -> 'SparseMatrix':
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = self._format

        if isinstance(other, SparseMatrix):
            result._mat = self._mat - other._mat
        elif isinstance(other, np.ndarray):
            result._mat = self._mat - csc_matrix(other)
        else:
            result._mat = self._mat - other

        return result

    def __mul__(self, other: Union['SparseMatrix', np.ndarray, float]) -> 'SparseMatrix':
        """Element-wise multiplication."""
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = self._format

        if isinstance(other, SparseMatrix):
            result._mat = self._mat.multiply(other._mat)
        elif isinstance(other, np.ndarray):
            result._mat = self._mat.multiply(other)
        else:
            result._mat = self._mat * other

        return result

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other: float) -> 'SparseMatrix':
        """Scalar division."""
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = self._format
        result._mat = self._mat / other
        return result

    def __neg__(self) -> 'SparseMatrix':
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = self._format
        result._mat = -self._mat
        return result

    def __matmul__(self, other: Union['SparseMatrix', np.ndarray]) -> Union['SparseMatrix', np.ndarray]:
        """Matrix multiplication."""
        if isinstance(other, SparseMatrix):
            result = SparseMatrix.__new__(SparseMatrix)
            result._format = self._format
            result._mat = self._mat @ other._mat
            return result
        elif isinstance(other, np.ndarray):
            # Sparse @ dense -> dense
            return self._mat @ other
        else:
            raise TypeError(f"Cannot multiply SparseMatrix with {type(other)}")

    def dot(self, other: Union['SparseMatrix', np.ndarray]) -> Union['SparseMatrix', np.ndarray]:
        """Matrix multiplication (alias for @)."""
        return self @ other

    # ========== Matrix Operations ==========

    def transpose(self) -> 'SparseMatrix':
        """Transpose matrix."""
        result = SparseMatrix.__new__(SparseMatrix)
        result._format = 'csr' if self._format == 'csc' else 'csc'
        result._mat = self._mat.T
        return result

    @property
    def T(self) -> 'SparseMatrix':
        """Transpose."""
        return self.transpose()

    def sum(self, axis: Optional[int] = None) -> Union[float, np.ndarray]:
        """Sum elements."""
        result = self._mat.sum(axis=axis)
        if axis is not None:
            return np.asarray(result).flatten()
        return float(result)

    def diagonal(self) -> np.ndarray:
        """Extract main diagonal."""
        return self._mat.diagonal()

    def row(self, i: int) -> np.ndarray:
        """Get row as dense array."""
        return np.asarray(self._mat.getrow(i).toarray()).flatten()

    def col(self, j: int) -> np.ndarray:
        """Get column as dense array."""
        return np.asarray(self._mat.getcol(j).toarray()).flatten()

    def __getitem__(self, key) -> Union['SparseMatrix', float, np.ndarray]:
        """Indexing."""
        result = self._mat[key]
        if sparse.issparse(result):
            out = SparseMatrix.__new__(SparseMatrix)
            out._mat = result
            out._format = self._format
            return out
        return result

    # ========== Linear Algebra ==========

    def solve(self, b: np.ndarray) -> np.ndarray:
        """
        Solve Ax = b.

        Args:
            b: Right-hand side vector

        Returns:
            Solution vector x
        """
        return spsolve(self._mat.tocsc(), b)

    def lu(self) -> 'SparseLU':
        """
        Compute sparse LU decomposition.

        Returns:
            SparseLU factorization object
        """
        return SparseLU(splu(self._mat.tocsc()))

    def eigs(
        self,
        k: int = 6,
        which: str = 'LM',
        sigma: Optional[float] = None,
        return_eigenvectors: bool = True,
    ) -> Union[np.ndarray, Tuple[np.ndarray, np.ndarray]]:
        """
        Compute largest eigenvalues and eigenvectors.

        Args:
            k: Number of eigenvalues
            which: Which eigenvalues ('LM', 'SM', 'LR', etc.)
            sigma: Shift for shift-invert mode
            return_eigenvectors: Whether to return eigenvectors

        Returns:
            eigenvalues, or (eigenvalues, eigenvectors)
        """
        return eigs(
            self._mat,
            k=k,
            which=which,
            sigma=sigma,
            return_eigenvectors=return_eigenvectors,
        )

    def norm(self, ord: str = 'fro') -> float:
        """
        Compute matrix norm.

        Args:
            ord: Norm type ('fro', 1, np.inf)
        """
        return sparse.linalg.norm(self._mat, ord=ord)

    # ========== Representation ==========

    def __repr__(self) -> str:
        return f"SparseMatrix({self.shape}, nnz={self.nnz}, format={self._format})"


class SparseLU:
    """
    Sparse LU decomposition wrapper.

    Provides efficient repeated solves with the same matrix.
    """

    def __init__(self, lu):
        """
        Args:
            lu: scipy.sparse.linalg.SuperLU object
        """
        self._lu = lu

    def solve(self, b: np.ndarray) -> np.ndarray:
        """Solve Ax = b using cached factorization."""
        return self._lu.solve(b)


# ========== Utility Functions ==========

def sparse_from_triplets(
    row: np.ndarray,
    col: np.ndarray,
    data: np.ndarray,
    shape: Tuple[int, int],
    format: str = 'csc',
) -> SparseMatrix:
    """
    Create sparse matrix from COO triplets.

    Args:
        row: Row indices
        col: Column indices
        data: Non-zero values
        shape: Matrix shape
        format: Output format

    Returns:
        SparseMatrix
    """
    mat = coo_matrix((data, (row, col)), shape=shape)
    return SparseMatrix(mat, format=format)


def sparse_block_diag(blocks: list, format: str = 'csc') -> SparseMatrix:
    """
    Create block diagonal sparse matrix.

    Args:
        blocks: List of matrices (sparse or dense)
        format: Output format

    Returns:
        Block diagonal SparseMatrix
    """
    scipy_blocks = []
    for b in blocks:
        if isinstance(b, SparseMatrix):
            scipy_blocks.append(b._mat)
        else:
            scipy_blocks.append(b)

    mat = sparse.block_diag(scipy_blocks, format=format)
    return SparseMatrix(mat)


def sparse_vstack(matrices: list, format: str = 'csc') -> SparseMatrix:
    """Vertically stack sparse matrices."""
    scipy_mats = [m._mat if isinstance(m, SparseMatrix) else m for m in matrices]
    return SparseMatrix(sparse.vstack(scipy_mats, format=format))


def sparse_hstack(matrices: list, format: str = 'csc') -> SparseMatrix:
    """Horizontally stack sparse matrices."""
    scipy_mats = [m._mat if isinstance(m, SparseMatrix) else m for m in matrices]
    return SparseMatrix(sparse.hstack(scipy_mats, format=format))
