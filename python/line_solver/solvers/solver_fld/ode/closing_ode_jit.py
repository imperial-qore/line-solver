"""
JIT-compiled closing method ODE kernels for SolverFLD.

Provides Just-In-Time compiled versions of closing method theta computation
for 10-100x speedup on large models (100+ states).

Hot loops (closing_ode.py):
- Lines 53-66: compute_theta_pnorm loop
- Lines 94+: compute_theta_softmin loop

This module JIT-compiles the theta computation for both methods.
Falls back to pure Python closing_ode module if Numba is not available.

Reference: closing_ode.py
License: MIT (same as LINE)
"""

import numpy as np
from typing import Optional
from ..utils.numba_support import HAS_NUMBA, njit, prange

if HAS_NUMBA:
    # JIT-compiled kernel functions

    @njit(fastmath=True, cache=True, parallel=True)
    def compute_theta_pnorm_jit(
        x: np.ndarray,
        c: np.ndarray,
        SQ: np.ndarray,
        p: float = 20.0,
    ) -> np.ndarray:
        """
        JIT-compiled theta computation using p-norm smoothing.

        Args:
            x: State vector (queue lengths per phase)
            c: Server capacities per state
            SQ: State-to-station aggregation matrix
            p: P-norm parameter

        Returns:
            theta: Effective service rates per state
        """
        x = np.maximum(x, 0.0)
        sum_x_Qa = SQ @ x + 1e-14

        theta = np.zeros_like(x)

        for i in prange(len(x)):
            x_val = sum_x_Qa[i]
            c_val = c[i]
            p_val = max(p, 1.0)

            if c_val > 0.0 and x_val > 0.0:
                ratio = x_val / c_val
                ghat = 1.0 / ((1.0 + ratio ** p_val) ** (1.0 / p_val))
                if ghat < 0.0 or not np.isfinite(ghat):
                    ghat = 0.0
                theta[i] = ghat
            else:
                theta[i] = 1.0

        return theta

    @njit(fastmath=True, cache=True, parallel=True)
    def compute_theta_softmin_jit(
        x: np.ndarray,
        c: np.ndarray,
        SQ: np.ndarray,
        alpha: float = 20.0,
    ) -> np.ndarray:
        """
        JIT-compiled theta computation using softmin smoothing.

        Args:
            x: State vector
            c: Server capacities per state
            SQ: State-to-station aggregation matrix
            alpha: Softmin parameter

        Returns:
            theta: Effective service rates
        """
        x = np.maximum(x, 0.0)
        sum_x_Qa = SQ @ x + 1e-14

        theta = np.zeros_like(x)

        for i in prange(len(x)):
            x_val = sum_x_Qa[i]
            c_val = c[i]

            if c_val > 0.0 and x_val > 0.0:
                # Softmin using log-sum-exp trick
                exp_x = np.exp(-alpha * x_val)
                exp_c = np.exp(-alpha * c_val)

                if exp_x >= exp_c:
                    log_sum = -alpha * x_val + np.log1p(np.exp(-alpha * (c_val - x_val)))
                else:
                    log_sum = -alpha * c_val + np.log1p(np.exp(-alpha * (x_val - c_val)))

                # softmin / x_val
                theta[i] = (log_sum / (-alpha)) / x_val if x_val > 0 else 1.0

                # Bounds check
                if theta[i] < 0.0 or not np.isfinite(theta[i]):
                    theta[i] = 0.0
                elif theta[i] > 1.0:
                    theta[i] = 1.0
            else:
                theta[i] = 1.0

        return theta

    @njit(fastmath=True, cache=True)
    def fluid_ode_closing_jit(
        t: float,
        x: np.ndarray,
        W: np.ndarray,
        A: np.ndarray,
        theta: np.ndarray,
        ALambda: np.ndarray,
    ) -> np.ndarray:
        """
        JIT-compiled closing method ODE right-hand side.

        dx/dt = W' * theta + A * lambda

        Args:
            t: Current time (unused)
            x: State vector (n,)
            W: Transition rate matrix (n, n)
            A: Initial phase probability matrix
            theta: Effective service rates (n,)
            ALambda: External arrival rates (n,)

        Returns:
            dxdt: State derivatives (n,)
        """
        # dx/dt = W' * theta + A * ALambda
        dxdt = W.T @ theta + A @ ALambda.flatten()
        return dxdt

else:
    # Fallback: import from pure Python closing_ode module
    from .closing_ode import (
        compute_theta_pnorm as compute_theta_pnorm_jit,
        compute_theta_softmin as compute_theta_softmin_jit,
        fluid_ode_closing as fluid_ode_closing_jit,
    )


__all__ = [
    'HAS_NUMBA',
    'compute_theta_pnorm_jit',
    'compute_theta_softmin_jit',
    'fluid_ode_closing_jit',
]
