"""
Core tensor abstractions for LINE native Python implementation.

This module provides:
- LineTensor: Unified tensor API supporting PyTorch and NumPy backends
- SparseMatrix: Sparse matrix operations with CSC format
- Linear algebra operations with backend dispatch
- Configuration for backend selection
"""

from .config import Config, Backend, get_config, set_backend
from .tensor import LineTensor, as_tensor
from .sparse import SparseMatrix
from .linalg import solve, inv, eig, svd, norm, det

__all__ = [
    'Config',
    'Backend',
    'get_config',
    'set_backend',
    'LineTensor',
    'as_tensor',
    'SparseMatrix',
    'solve',
    'inv',
    'eig',
    'svd',
    'norm',
    'det',
]
