"""
Numba JIT compilation support for SolverFLD.

Provides optional Just-In-Time (JIT) compilation of numerical kernels for 10-100x speedup
on large-scale models (100+ states). Falls back gracefully to pure Python if Numba
is not available.

Pattern from: python/line_solver/api/mc/ctmc.py (lines 24-33)

Usage:
    >>> from .numba_support import HAS_NUMBA, get_numba_status
    >>> if HAS_NUMBA:
    ...     from .pnorm_jit import pnorm_smooth_vectorized_jit
    ...     result = pnorm_smooth_vectorized_jit(x, c, p)
    ... else:
    ...     from .pnorm import pnorm_smooth
    ...     result = vectorized_call(pnorm_smooth, x, c, p)

License: MIT (same as LINE)
"""

import numpy as np
from typing import Dict, Any, Tuple, Callable

# Try to import Numba
try:
    from numba import njit, prange, config as numba_config
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

    # Graceful fallback: create dummy decorators
    def njit(*args, **kwargs):
        """Decorator that does nothing if Numba is not available."""
        def decorator(func):
            return func
        # Handle both @njit and @njit() syntax
        if args and callable(args[0]):
            return args[0]
        return decorator

    def prange(*args, **kwargs):
        """Fallback for prange: just return regular range."""
        return range(*args)


def get_numba_status() -> Dict[str, Any]:
    """Check Numba availability and configuration.

    Returns
    -------
    status : dict
        Dictionary with keys:
        - 'available': bool, whether Numba is installed
        - 'version': str, Numba version (if available)
        - 'fastmath': bool, whether fastmath is enabled
        - 'parallel': bool, whether parallel mode is enabled
        - 'message': str, status message
    """
    status = {
        'available': HAS_NUMBA,
        'version': None,
        'fastmath': False,
        'parallel': False,
        'message': ''
    }

    if HAS_NUMBA:
        try:
            import numba
            status['version'] = numba.__version__
            # Check if parallelization is available
            try:
                status['parallel'] = numba_config.DISABLE_JIT == False
            except:
                status['parallel'] = True  # Default to enabled
            status['fastmath'] = True  # We use fastmath=True in our decorators
            status['message'] = f"Numba {status['version']} available (JIT enabled)"
        except Exception as e:
            status['message'] = f"Numba import error: {str(e)}"
    else:
        status['message'] = "Numba not available (using Python fallback)"

    return status


def benchmark_speedup(
    python_func: Callable,
    jit_func: Callable,
    args: Tuple,
    n_runs: int = 10,
    description: str = "JIT"
) -> Dict[str, Any]:
    """Benchmark Python function vs JIT-compiled version.

    Parameters
    ----------
    python_func : callable
        Pure Python reference function
    jit_func : callable
        JIT-compiled function (may be same as python_func if JIT unavailable)
    args : tuple
        Arguments to pass to both functions
    n_runs : int, optional
        Number of runs for timing (default 10)
    description : str, optional
        Description of benchmark for display (default "JIT")

    Returns
    -------
    benchmark : dict
        Dictionary with keys:
        - 'python_time': float, Python execution time (ms)
        - 'jit_time': float, JIT execution time (ms)
        - 'speedup': float, speedup factor (python_time / jit_time)
        - 'runs': int, number of timing runs
        - 'description': str, benchmark description
    """
    import time

    # Warmup: call JIT function to trigger compilation
    try:
        jit_func(*args)
    except:
        pass

    # Time Python function
    start = time.time()
    for _ in range(n_runs):
        python_func(*args)
    python_time = (time.time() - start) / n_runs * 1000  # Convert to ms

    # Time JIT function
    start = time.time()
    for _ in range(n_runs):
        jit_func(*args)
    jit_time = (time.time() - start) / n_runs * 1000  # Convert to ms

    speedup = python_time / jit_time if jit_time > 0 else float('inf')

    return {
        'python_time': python_time,
        'jit_time': jit_time,
        'speedup': speedup,
        'runs': n_runs,
        'description': description,
    }


__all__ = [
    'HAS_NUMBA',
    'njit',
    'prange',
    'get_numba_status',
    'benchmark_speedup',
]
