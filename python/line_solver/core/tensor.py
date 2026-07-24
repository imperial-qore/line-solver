"""
Unified tensor API for LINE native implementation.

Provides LineTensor class that abstracts PyTorch and NumPy backends,
enabling seamless switching between CPU/GPU computation and autodiff support.
"""

from __future__ import annotations

import numpy as np
from typing import Union, Tuple, Optional, List, Any, TYPE_CHECKING

from .config import get_config, Backend

if TYPE_CHECKING:
    import torch

# Type alias for array-like inputs
ArrayLike = Union[np.ndarray, List, Tuple, 'torch.Tensor', 'LineTensor', float, int]


class LineTensor:
    """
    Unified tensor supporting PyTorch (GPU/autodiff) and NumPy (Numba JIT).

    This class provides a consistent API regardless of backend, enabling:
    - GPU acceleration via PyTorch when available
    - Automatic differentiation for optimization
    - Numba JIT compatibility via NumPy arrays
    - Seamless conversion between backends

    Examples:
        >>> t = LineTensor([1, 2, 3])
        >>> t.shape
        (3,)
        >>> t_gpu = t.to('cuda')
        >>> t_np = t.to_numpy()
    """

    __slots__ = ('_data', '_backend', '_device', '_requires_grad')

    def __init__(
        self,
        data: ArrayLike,
        backend: Optional[str] = None,
        device: Optional[str] = None,
        dtype: Optional[str] = None,
        requires_grad: bool = False,
    ):
        """
        Create a LineTensor from array-like data.

        Args:
            data: Input data (array, list, scalar, or tensor)
            backend: Backend to use ('numpy', 'pytorch', or None for auto)
            device: Device for computation ('cpu', 'cuda', 'cuda:0')
            dtype: Data type ('float32', 'float64', etc.)
            requires_grad: Enable gradient tracking (PyTorch only)
        """
        config = get_config()

        # Determine backend
        if backend is None:
            self._backend = config.get_effective_backend()
        else:
            self._backend = Backend(backend.lower())

        # Determine device
        if device is None:
            self._device = config.get_effective_device()
        else:
            self._device = device

        self._requires_grad = requires_grad

        # Determine dtype
        if dtype is None:
            dtype = config.default_dtype

        # Convert input data to appropriate format
        if isinstance(data, LineTensor):
            self._init_from_linetensor(data, dtype)
        elif self._backend == Backend.PYTORCH and config.pytorch_available:
            self._init_pytorch(data, dtype)
        else:
            self._init_numpy(data, dtype)

    def _init_numpy(self, data: ArrayLike, dtype: str) -> None:
        """Initialize with NumPy backend."""
        np_dtype = getattr(np, dtype, np.float64)

        if hasattr(data, 'numpy'):  # PyTorch tensor
            self._data = data.detach().cpu().numpy().astype(np_dtype)
        elif isinstance(data, np.ndarray):
            self._data = data.astype(np_dtype) if data.dtype != np_dtype else data
        else:
            self._data = np.asarray(data, dtype=np_dtype)

        self._backend = Backend.NUMPY

    def _init_pytorch(self, data: ArrayLike, dtype: str) -> None:
        """Initialize with PyTorch backend."""
        import torch

        # Map dtype string to torch dtype
        dtype_map = {
            'float32': torch.float32,
            'float64': torch.float64,
            'complex64': torch.complex64,
            'complex128': torch.complex128,
            'int32': torch.int32,
            'int64': torch.int64,
        }
        torch_dtype = dtype_map.get(dtype, torch.float64)

        if isinstance(data, torch.Tensor):
            self._data = data.to(dtype=torch_dtype, device=self._device)
        elif isinstance(data, np.ndarray):
            self._data = torch.from_numpy(data).to(dtype=torch_dtype, device=self._device)
        else:
            self._data = torch.tensor(data, dtype=torch_dtype, device=self._device)

        if self._requires_grad and self._data.is_floating_point():
            self._data = self._data.requires_grad_(True)

    def _init_from_linetensor(self, other: 'LineTensor', dtype: str) -> None:
        """Initialize from another LineTensor."""
        if self._backend == other._backend:
            self._data = other._data.copy() if self._backend == Backend.NUMPY else other._data.clone()
        elif self._backend == Backend.PYTORCH:
            self._init_pytorch(other.to_numpy(), dtype)
        else:
            self._init_numpy(other.to_numpy(), dtype)

    # ========== Properties ==========

    @property
    def shape(self) -> Tuple[int, ...]:
        """Tensor shape."""
        return tuple(self._data.shape)

    @property
    def ndim(self) -> int:
        """Number of dimensions."""
        return self._data.ndim

    @property
    def size(self) -> int:
        """Total number of elements."""
        return int(np.prod(self.shape))

    @property
    def dtype(self) -> str:
        """Data type as string."""
        if self._backend == Backend.PYTORCH:
            return str(self._data.dtype).replace('torch.', '')
        return str(self._data.dtype)

    @property
    def device(self) -> str:
        """Device where tensor resides."""
        if self._backend == Backend.PYTORCH:
            return str(self._data.device)
        return 'cpu'

    @property
    def backend(self) -> Backend:
        """Current backend."""
        return self._backend

    @property
    def requires_grad(self) -> bool:
        """Whether gradient tracking is enabled."""
        if self._backend == Backend.PYTORCH:
            return self._data.requires_grad
        return False

    @property
    def T(self) -> 'LineTensor':
        """Transpose."""
        return self.transpose()

    # ========== Conversion Methods ==========

    def to_numpy(self) -> np.ndarray:
        """
        Convert to NumPy array.

        Returns:
            NumPy array (detached from computation graph if PyTorch)
        """
        if self._backend == Backend.NUMPY:
            return self._data
        # PyTorch
        return self._data.detach().cpu().numpy()

    def to_pytorch(self, device: Optional[str] = None, requires_grad: bool = False) -> 'torch.Tensor':
        """
        Convert to PyTorch tensor.

        Args:
            device: Target device
            requires_grad: Enable gradient tracking

        Returns:
            PyTorch tensor
        """
        import torch

        if device is None:
            device = self._device

        if self._backend == Backend.PYTORCH:
            t = self._data.to(device)
        else:
            t = torch.from_numpy(self._data).to(device)

        if requires_grad and t.is_floating_point():
            t = t.requires_grad_(True)

        return t

    def to(self, device: str) -> 'LineTensor':
        """
        Move tensor to specified device.

        Args:
            device: Target device ('cpu', 'cuda', 'cuda:0')

        Returns:
            New LineTensor on target device
        """
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.to(device)
        else:
            result._data = self._data  # NumPy is always on CPU

        return result

    def numpy(self) -> np.ndarray:
        """Alias for to_numpy()."""
        return self.to_numpy()

    def clone(self) -> 'LineTensor':
        """Create a copy of the tensor."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.clone()
        else:
            result._data = self._data.copy()

        return result

    def detach(self) -> 'LineTensor':
        """Detach from computation graph (PyTorch only)."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = False

        if self._backend == Backend.PYTORCH:
            result._data = self._data.detach()
        else:
            result._data = self._data

        return result

    # ========== Shape Operations ==========

    def reshape(self, *shape) -> 'LineTensor':
        """Reshape tensor."""
        if len(shape) == 1 and isinstance(shape[0], (tuple, list)):
            shape = tuple(shape[0])

        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad
        result._data = self._data.reshape(shape)
        return result

    def flatten(self) -> 'LineTensor':
        """Flatten to 1D."""
        return self.reshape(-1)

    def squeeze(self, dim: Optional[int] = None) -> 'LineTensor':
        """Remove dimensions of size 1."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if dim is None:
            result._data = self._data.squeeze()
        else:
            if self._backend == Backend.PYTORCH:
                result._data = self._data.squeeze(dim)
            else:
                result._data = self._data.squeeze(axis=dim)
        return result

    def unsqueeze(self, dim: int) -> 'LineTensor':
        """Add dimension of size 1."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.unsqueeze(dim)
        else:
            result._data = np.expand_dims(self._data, axis=dim)
        return result

    def transpose(self, dim0: int = -2, dim1: int = -1) -> 'LineTensor':
        """Transpose dimensions."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self.ndim < 2:
            result._data = self._data
        elif self._backend == Backend.PYTORCH:
            result._data = self._data.transpose(dim0, dim1)
        else:
            # NumPy transpose
            if self.ndim == 2:
                result._data = self._data.T
            else:
                axes = list(range(self.ndim))
                axes[dim0], axes[dim1] = axes[dim1], axes[dim0]
                result._data = self._data.transpose(axes)
        return result

    # ========== Arithmetic Operations ==========

    def __add__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__add__')

    def __radd__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__radd__')

    def __sub__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__sub__')

    def __rsub__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__rsub__')

    def __mul__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__mul__')

    def __rmul__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__rmul__')

    def __truediv__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__truediv__')

    def __rtruediv__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__rtruediv__')

    def __pow__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        return self._binary_op(other, '__pow__')

    def __neg__(self) -> 'LineTensor':
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad
        result._data = -self._data
        return result

    def __matmul__(self, other: Union['LineTensor', ArrayLike]) -> 'LineTensor':
        """Matrix multiplication."""
        return self._binary_op(other, '__matmul__')

    def _binary_op(self, other: Union['LineTensor', ArrayLike], op: str) -> 'LineTensor':
        """Execute binary operation."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if isinstance(other, LineTensor):
            other_data = other._data
        elif isinstance(other, (int, float)):
            other_data = other
        else:
            other_data = np.asarray(other)
            if self._backend == Backend.PYTORCH:
                import torch
                other_data = torch.from_numpy(other_data).to(self._data.device)

        result._data = getattr(self._data, op)(other_data)
        return result

    # ========== Reduction Operations ==========

    def sum(self, axis: Optional[int] = None, keepdims: bool = False) -> 'LineTensor':
        """Sum over axis."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            if axis is None:
                result._data = self._data.sum()
            else:
                result._data = self._data.sum(dim=axis, keepdim=keepdims)
        else:
            result._data = self._data.sum(axis=axis, keepdims=keepdims)
        return result

    def mean(self, axis: Optional[int] = None, keepdims: bool = False) -> 'LineTensor':
        """Mean over axis."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            if axis is None:
                result._data = self._data.mean()
            else:
                result._data = self._data.mean(dim=axis, keepdim=keepdims)
        else:
            result._data = self._data.mean(axis=axis, keepdims=keepdims)
        return result

    def max(self, axis: Optional[int] = None, keepdims: bool = False) -> 'LineTensor':
        """Maximum over axis."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            if axis is None:
                result._data = self._data.max()
            else:
                result._data = self._data.max(dim=axis, keepdim=keepdims).values
        else:
            result._data = self._data.max(axis=axis, keepdims=keepdims)
        return result

    def min(self, axis: Optional[int] = None, keepdims: bool = False) -> 'LineTensor':
        """Minimum over axis."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            if axis is None:
                result._data = self._data.min()
            else:
                result._data = self._data.min(dim=axis, keepdim=keepdims).values
        else:
            result._data = self._data.min(axis=axis, keepdims=keepdims)
        return result

    def prod(self, axis: Optional[int] = None, keepdims: bool = False) -> 'LineTensor':
        """Product over axis."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            if axis is None:
                result._data = self._data.prod()
            else:
                result._data = self._data.prod(dim=axis, keepdim=keepdims)
        else:
            result._data = self._data.prod(axis=axis, keepdims=keepdims)
        return result

    # ========== Element-wise Operations ==========

    def abs(self) -> 'LineTensor':
        """Absolute value."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.abs()
        else:
            result._data = np.abs(self._data)
        return result

    def sqrt(self) -> 'LineTensor':
        """Square root."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.sqrt()
        else:
            result._data = np.sqrt(self._data)
        return result

    def exp(self) -> 'LineTensor':
        """Exponential."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.exp()
        else:
            result._data = np.exp(self._data)
        return result

    def log(self) -> 'LineTensor':
        """Natural logarithm."""
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad

        if self._backend == Backend.PYTORCH:
            result._data = self._data.log()
        else:
            result._data = np.log(self._data)
        return result

    # ========== Indexing ==========

    def __getitem__(self, key) -> 'LineTensor':
        result = LineTensor.__new__(LineTensor)
        result._backend = self._backend
        result._device = self._device
        result._requires_grad = self._requires_grad
        result._data = self._data[key]
        return result

    def __setitem__(self, key, value) -> None:
        if isinstance(value, LineTensor):
            value = value._data
        self._data[key] = value

    # ========== Comparison ==========

    def __eq__(self, other) -> np.ndarray:
        if isinstance(other, LineTensor):
            other = other._data
        return self._data == other

    def __lt__(self, other) -> np.ndarray:
        if isinstance(other, LineTensor):
            other = other._data
        return self._data < other

    def __le__(self, other) -> np.ndarray:
        if isinstance(other, LineTensor):
            other = other._data
        return self._data <= other

    def __gt__(self, other) -> np.ndarray:
        if isinstance(other, LineTensor):
            other = other._data
        return self._data > other

    def __ge__(self, other) -> np.ndarray:
        if isinstance(other, LineTensor):
            other = other._data
        return self._data >= other

    # ========== Representation ==========

    def __repr__(self) -> str:
        return f"LineTensor({self._data}, backend={self._backend.value}, device={self.device})"

    def __str__(self) -> str:
        return str(self._data)

    def __len__(self) -> int:
        return len(self._data)

    # ========== Scalar conversion ==========

    def item(self) -> float:
        """Extract scalar value."""
        if self._backend == Backend.PYTORCH:
            return self._data.item()
        return float(self._data)

    def __float__(self) -> float:
        return self.item()

    def __int__(self) -> int:
        return int(self.item())


# ========== Factory Functions ==========

def as_tensor(
    data: ArrayLike,
    backend: Optional[str] = None,
    device: Optional[str] = None,
    dtype: Optional[str] = None,
    requires_grad: bool = False,
) -> LineTensor:
    """
    Convert data to LineTensor.

    Args:
        data: Input data
        backend: Backend ('numpy', 'pytorch', or None for auto)
        device: Device ('cpu', 'cuda')
        dtype: Data type
        requires_grad: Enable gradient tracking

    Returns:
        LineTensor
    """
    if isinstance(data, LineTensor):
        if backend is None and device is None and dtype is None:
            return data
    return LineTensor(data, backend=backend, device=device, dtype=dtype, requires_grad=requires_grad)


def zeros(shape: Tuple[int, ...], **kwargs) -> LineTensor:
    """Create tensor of zeros."""
    return LineTensor(np.zeros(shape), **kwargs)


def ones(shape: Tuple[int, ...], **kwargs) -> LineTensor:
    """Create tensor of ones."""
    return LineTensor(np.ones(shape), **kwargs)


def eye(n: int, **kwargs) -> LineTensor:
    """Create identity matrix."""
    return LineTensor(np.eye(n), **kwargs)


def arange(start: int, stop: Optional[int] = None, step: int = 1, **kwargs) -> LineTensor:
    """Create range tensor."""
    if stop is None:
        stop = start
        start = 0
    return LineTensor(np.arange(start, stop, step), **kwargs)


def linspace(start: float, stop: float, num: int, **kwargs) -> LineTensor:
    """Create linearly spaced tensor."""
    return LineTensor(np.linspace(start, stop, num), **kwargs)


def empty(shape: Tuple[int, ...], **kwargs) -> LineTensor:
    """Create uninitialized tensor."""
    return LineTensor(np.empty(shape), **kwargs)


def diag(v: ArrayLike, **kwargs) -> LineTensor:
    """Create diagonal matrix from vector or extract diagonal."""
    if isinstance(v, LineTensor):
        v = v.to_numpy()
    return LineTensor(np.diag(v), **kwargs)


def concatenate(tensors: List[LineTensor], axis: int = 0) -> LineTensor:
    """Concatenate tensors along axis."""
    if not tensors:
        raise ValueError("Need at least one tensor to concatenate")

    backend = tensors[0]._backend
    device = tensors[0]._device

    if backend == Backend.PYTORCH:
        import torch
        data = torch.cat([t._data for t in tensors], dim=axis)
    else:
        data = np.concatenate([t._data for t in tensors], axis=axis)

    result = LineTensor.__new__(LineTensor)
    result._backend = backend
    result._device = device
    result._requires_grad = any(t._requires_grad for t in tensors)
    result._data = data
    return result


def stack(tensors: List[LineTensor], axis: int = 0) -> LineTensor:
    """Stack tensors along new axis."""
    if not tensors:
        raise ValueError("Need at least one tensor to stack")

    backend = tensors[0]._backend
    device = tensors[0]._device

    if backend == Backend.PYTORCH:
        import torch
        data = torch.stack([t._data for t in tensors], dim=axis)
    else:
        data = np.stack([t._data for t in tensors], axis=axis)

    result = LineTensor.__new__(LineTensor)
    result._backend = backend
    result._device = device
    result._requires_grad = any(t._requires_grad for t in tensors)
    result._data = data
    return result
