"""
Backend selection for LINE solver.

This module provides the native Python implementation backend.
All API calls use native Python implementations.

Environment Variables:
    LINE_BACKEND: Tensor backend ('numpy', 'pytorch', 'auto')
    LINE_DEVICE: Compute device ('cpu', 'cuda', 'cuda:0')

Usage:
    from line_solver._backend import get_api

    pfqn = get_api('pfqn')  # Returns native pfqn module
"""

import os
import importlib
from typing import Optional
from enum import Enum


class APIBackend(Enum):
    """API implementation backend."""
    NATIVE = "native"


# Global state
_api_backend = APIBackend.NATIVE
_native_available = None


def _check_native_available() -> bool:
    """Check if native implementations are available."""
    global _native_available
    if _native_available is None:
        try:
            from line_solver import api
            _native_available = True
        except ImportError:
            _native_available = False
    return _native_available


def get_api_backend() -> APIBackend:
    """Get current API backend setting."""
    return _api_backend


def use_native() -> None:
    """Force use of native Python implementations (default)."""
    pass  # Already native-only


def get_effective_backend() -> APIBackend:
    """
    Determine effective backend based on settings and availability.

    Returns:
        APIBackend.NATIVE
    """
    return APIBackend.NATIVE


def get_api(module_name: str):
    """
    Get API module (native implementation).

    Args:
        module_name: Name of API module (e.g., 'pfqn', 'mc', 'qsys')

    Returns:
        Module object (native version)

    Examples:
        >>> pfqn = get_api('pfqn')
        >>> result = pfqn.pfqn_mva(N, L, Z)
    """
    try:
        return importlib.import_module(f'line_solver.api.{module_name}')
    except ImportError as e:
        raise ImportError(
            f"Native module 'api.{module_name}' not available: {e}"
        )


def get_function(module_name: str, function_name: str):
    """
    Get specific function from API module.

    Args:
        module_name: Name of API module
        function_name: Name of function

    Returns:
        Function object

    Examples:
        >>> pfqn_mva = get_function('pfqn', 'pfqn_mva')
        >>> result = pfqn_mva(N, L, Z)
    """
    module = get_api(module_name)
    return getattr(module, function_name)


class BackendContext:
    """
    Context manager for backend (native only).

    Examples:
        >>> with BackendContext('native'):
        ...     result = pfqn.pfqn_mva(N, L, Z)  # Uses native
    """

    def __init__(self, backend: str = 'native'):
        if backend.lower() != 'native':
            raise ValueError("Only 'native' backend is supported")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        return False
