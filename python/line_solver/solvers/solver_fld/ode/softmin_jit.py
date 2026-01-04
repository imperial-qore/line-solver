"""
JIT-compiled softmin smoothing kernels for SolverFLD.

Provides Just-In-Time compiled versions of softmin capacity constraint
computation for 10-100x speedup on large models (100+ states).

Hot loop (softmin.py lines 98-99):
    for i in range(n):
        g[i] = softmin_smooth(x[i], c[i], alpha)

This module JIT-compiles the above loop and the full ODE computation.
Falls back to pure Python softmin module if Numba is not available.

Reference: softmin.py
License: MIT (same as LINE)
"""

import numpy as np
from typing import Optional
from ..utils.numba_support import HAS_NUMBA, njit, prange

if HAS_NUMBA:
    # JIT-compiled kernel functions

    @njit(fastmath=True, cache=True)
    def softmin_smooth_jit(x: float, c: float, alpha: float = 20.0) -> float:
        """
        JIT-compiled softmin smoothed constraint.

        Args:
            x: Total queue length
            c: Server capacity
            alpha: Softmin parameter

        Returns:
            g: Smooth approximation to min(x, c) / x
        """
        if c <= 0.0 or alpha <= 0.0 or x <= 0.0:
            return 1.0

        # Numerical stability: log-sum-exp trick
        exp_x = np.exp(-alpha * x)
        exp_c = np.exp(-alpha * c)

        if exp_x >= exp_c:
            # x <= c case
            log_sum = -alpha * x + np.log1p(np.exp(-alpha * (c - x)))
        else:
            # x > c case
            log_sum = -alpha * c + np.log1p(np.exp(-alpha * (x - c)))

        # softmin = log_sum / (-alpha)
        softmin_val = log_sum / (-alpha)

        # g = softmin / x
        g = softmin_val / x

        # Bounds check
        if g < 0.0:
            return 0.0
        if g > 1.0:
            return 1.0

        return g

    @njit(fastmath=True, cache=True, parallel=True)
    def softmin_smooth_vectorized_jit(
        x: np.ndarray,
        c: np.ndarray,
        alpha: float = 20.0
    ) -> np.ndarray:
        """
        JIT-compiled vectorized softmin constraint computation (with parallelization).

        Args:
            x: Queue lengths per state (1-D array)
            c: Capacities per state (1-D array)
            alpha: Softmin parameter (scalar)

        Returns:
            g: Smooth constraints (same shape as x)
        """
        n = len(x)
        g = np.zeros(n)

        for i in prange(n):
            g[i] = softmin_smooth_jit(x[i], c[i], alpha)

        return g

    @njit(fastmath=True, cache=True)
    def fluid_ode_softmin_jit(
        t: float,
        x: np.ndarray,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
        alpha: float = 20.0,
    ) -> np.ndarray:
        """
        JIT-compiled ODE with softmin smoothed capacity constraints.

        Args:
            t: Current time (unused)
            x: State vector (n,)
            W: Transition rate matrix (n, n)
            SQ: State-to-station queue mapping (n, n)
            Sa: Server capacities per state (n,)
            ALambda: External arrival rates (n,)
            alpha: Softmin parameter

        Returns:
            dxdt: State derivatives (n,)
        """
        # Ensure non-negative state
        x = np.maximum(x, 0.0)

        # Compute total queue at each state's station
        sum_x_Qa = SQ @ x + 1e-14

        # Compute softmin constraints
        g = softmin_smooth_vectorized_jit(sum_x_Qa, Sa.flatten(), alpha)

        # Compute derivative: dx/dt = W @ (x * g) + ALambda
        dxdt = W @ (x * g) + ALambda.flatten()

        return dxdt

else:
    # Fallback: import from pure Python softmin module
    from .softmin import (
        softmin_smooth as softmin_smooth_jit,
        softmin_smooth_vectorized as softmin_smooth_vectorized_jit,
        fluid_ode_softmin as fluid_ode_softmin_jit,
    )


__all__ = [
    'HAS_NUMBA',
    'softmin_smooth_jit',
    'softmin_smooth_vectorized_jit',
    'fluid_ode_softmin_jit',
]
