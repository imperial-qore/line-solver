"""
JIT-compiled p-norm smoothing kernels for SolverFLD.

Provides Just-In-Time compiled versions of p-norm capacity constraint
computation for 10-100x speedup on large models (100+ states).

Hot loop (pnorm.py lines 87-88):
    for i in range(n):
        ghat[i] = pnorm_smooth(x[i], c[i], p[i])

This module JIT-compiles the above loop and the full ODE computation.
Falls back to pure Python pnorm module if Numba is not available.

Reference: pnorm.py
License: MIT (same as LINE)
"""

import numpy as np
from typing import Optional
from ..utils.numba_support import HAS_NUMBA, njit, prange

if HAS_NUMBA:
    # JIT-compiled kernel functions (compiled when first called)

    @njit(fastmath=True, cache=True)
    def pnorm_smooth_jit(x: float, c: float, p: float = 20.0) -> float:
        """
        JIT-compiled p-norm smoothed capacity constraint.

        Args:
            x: Total queue length at station
            c: Server capacity
            p: P-norm parameter

        Returns:
            g_hat: Smoothed constraint value in (0, 1]
        """
        if c <= 0.0 or p <= 0.0 or x <= 0.0:
            return 1.0

        # g_hat = 1 / (1 + (x/c)^p)^(1/p)
        ratio = x / c
        power_term = ratio ** p  # Numba optimizes ** better than np.power
        denominator = (1.0 + power_term) ** (1.0 / p)
        ghat = 1.0 / denominator

        # Bounds check
        if ghat < 0.0:
            return 0.0
        if ghat > 1.0:
            return 1.0

        return ghat

    @njit(fastmath=True, cache=True, parallel=True)
    def pnorm_smooth_vectorized_jit(
        x: np.ndarray,
        c: np.ndarray,
        p: Optional[np.ndarray] = None
    ) -> np.ndarray:
        """
        JIT-compiled vectorized p-norm constraint computation (with parallelization).

        Args:
            x: Queue lengths per state (1-D array)
            c: Capacities per state (1-D array)
            p: P-norm parameters per state (None = use 20.0)

        Returns:
            g_hat: Smoothed constraints (same shape as x)
        """
        n = len(x)
        ghat = np.zeros(n)

        if p is None:
            # Default p-norm parameter
            for i in prange(n):
                ghat[i] = pnorm_smooth_jit(x[i], c[i], 20.0)
        else:
            # Per-state p-norm parameter
            for i in prange(n):
                p_val = p[i] if i < len(p) else 20.0
                ghat[i] = pnorm_smooth_jit(x[i], c[i], p_val)

        return ghat

    @njit(fastmath=True, cache=True)
    def fluid_ode_pnorm_jit(
        t: float,
        x: np.ndarray,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
        pstar: Optional[np.ndarray] = None,
    ) -> np.ndarray:
        """
        JIT-compiled ODE right-hand side with p-norm smoothed capacity constraints.

        dx/dt = W * (x .* g_hat(x)) + ALambda

        Args:
            t: Current time (unused)
            x: State vector (n,)
            W: Transition rate matrix (n, n)
            SQ: State-to-station queue mapping (n, n)
            Sa: Server capacities per state (n,)
            ALambda: External arrival rates (n,)
            pstar: P-norm parameters or scalar

        Returns:
            dxdt: State derivatives (n,)
        """
        # Ensure non-negative state
        x = np.maximum(x, 0.0)

        # Compute total queue at each state's station
        sum_x_Qa = SQ @ x + 1e-14

        # Compute smoothed capacity constraints
        if pstar is not None:
            if isinstance(pstar, float):
                ghat = pnorm_smooth_vectorized_jit(
                    sum_x_Qa, Sa.flatten(), np.full(len(x), pstar)
                )
            else:
                ghat = pnorm_smooth_vectorized_jit(sum_x_Qa, Sa.flatten(), pstar)
        else:
            # Default p-norm with p=20
            ghat = pnorm_smooth_vectorized_jit(
                sum_x_Qa, Sa.flatten(), np.full(len(x), 20.0)
            )

        # Compute derivative: dx/dt = W @ (x * ghat) + ALambda
        dxdt = W @ (x * ghat) + ALambda.flatten()

        return dxdt

else:
    # Fallback: import from pure Python pnorm module
    from .pnorm import (
        pnorm_smooth as pnorm_smooth_jit,
        pnorm_smooth_vectorized as pnorm_smooth_vectorized_jit,
        fluid_ode_pnorm as fluid_ode_pnorm_jit,
    )


__all__ = [
    'HAS_NUMBA',
    'pnorm_smooth_jit',
    'pnorm_smooth_vectorized_jit',
    'fluid_ode_pnorm_jit',
]
