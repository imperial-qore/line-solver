"""
Configuration for LINE native backend.

Controls backend selection (PyTorch vs NumPy), device placement,
and numerical tolerances.
"""

import os
from enum import Enum
from typing import Optional
from dataclasses import dataclass, field


class Backend(Enum):
    """Available computation backends."""
    NUMPY = "numpy"
    PYTORCH = "pytorch"
    AUTO = "auto"


@dataclass
class Config:
    """
    Global configuration for LINE native implementation.

    Attributes:
        backend: Computation backend (numpy, pytorch, auto)
        device: Device for PyTorch tensors ('cpu', 'cuda', 'cuda:0', etc.)
        default_dtype: Default floating-point dtype
        use_jit: Enable Numba JIT compilation
        jit_cache: Cache JIT-compiled functions
        tol_fine: Fine tolerance for numerical comparisons (1e-12)
        tol_medium: Medium tolerance (1e-8)
        tol_coarse: Coarse tolerance (1e-4)
        max_iter: Default maximum iterations for iterative algorithms
    """
    backend: Backend = Backend.AUTO
    device: str = "cpu"
    default_dtype: str = "float64"
    use_jit: bool = True
    jit_cache: bool = True
    tol_fine: float = 1e-12
    tol_medium: float = 1e-8
    tol_coarse: float = 1e-4
    max_iter: int = 100000

    # PyTorch specific
    requires_grad_default: bool = False

    # Runtime state
    _pytorch_available: Optional[bool] = field(default=None, repr=False)
    _numba_available: Optional[bool] = field(default=None, repr=False)
    _cuda_available: Optional[bool] = field(default=None, repr=False)

    def __post_init__(self):
        # Check for environment variable overrides
        env_backend = os.environ.get('LINE_BACKEND', '').lower()
        if env_backend == 'numpy':
            self.backend = Backend.NUMPY
        elif env_backend == 'pytorch':
            self.backend = Backend.PYTORCH

        env_device = os.environ.get('LINE_DEVICE', '')
        if env_device:
            self.device = env_device

        env_jit = os.environ.get('LINE_JIT', '').lower()
        if env_jit == 'false' or env_jit == '0':
            self.use_jit = False

    @property
    def pytorch_available(self) -> bool:
        """Check if PyTorch is available."""
        if self._pytorch_available is None:
            try:
                import torch
                self._pytorch_available = True
            except ImportError:
                self._pytorch_available = False
        return self._pytorch_available

    @property
    def numba_available(self) -> bool:
        """Check if Numba is available."""
        if self._numba_available is None:
            try:
                import numba
                self._numba_available = True
            except ImportError:
                self._numba_available = False
        return self._numba_available

    @property
    def cuda_available(self) -> bool:
        """Check if CUDA is available for PyTorch."""
        if self._cuda_available is None:
            if self.pytorch_available:
                import torch
                self._cuda_available = torch.cuda.is_available()
            else:
                self._cuda_available = False
        return self._cuda_available

    def get_effective_backend(self) -> Backend:
        """
        Determine the effective backend based on AUTO selection.

        Returns:
            Backend.PYTORCH if available, otherwise Backend.NUMPY
        """
        if self.backend != Backend.AUTO:
            return self.backend

        if self.pytorch_available:
            return Backend.PYTORCH
        return Backend.NUMPY

    def get_effective_device(self) -> str:
        """
        Determine the effective device.

        Returns 'cpu' if CUDA requested but not available.
        """
        if 'cuda' in self.device and not self.cuda_available:
            return 'cpu'
        return self.device


# Global configuration instance
_config = Config()


def get_config() -> Config:
    """Get the global configuration."""
    return _config


def set_backend(backend: str) -> None:
    """
    Set the computation backend.

    Args:
        backend: One of 'numpy', 'pytorch', 'auto'
    """
    _config.backend = Backend(backend.lower())


def set_device(device: str) -> None:
    """
    Set the computation device.

    Args:
        device: Device string ('cpu', 'cuda', 'cuda:0', etc.)
    """
    _config.device = device


def enable_jit(enabled: bool = True) -> None:
    """Enable or disable Numba JIT compilation."""
    _config.use_jit = enabled
