"""
JIT-compiled state-dependent constraint kernels for SolverFLD.

Provides Just-In-Time compiled versions of hard capacity constraint
computation for 10-100x speedup on large models (100+ states).

Hot loop (statedep.py):
    Vectorized: g = min(1, c / x)

This module JIT-compiles the constraint computation.
Falls back to pure Python statedep module if Numba is not available.

Reference: statedep.py
License: MIT (same as LINE)
"""

import numpy as np
from typing import Optional
from ..utils.numba_support import HAS_NUMBA, njit, prange

if HAS_NUMBA:
    # JIT-compiled kernel functions

    @njit(fastmath=True, cache=True)
    def statedep_constraint_jit(x: float, c: float) -> float:
        """
        JIT-compiled state-dependent (hard) capacity constraint.

        Args:
            x: Total queue length
            c: Server capacity

        Returns:
            g: Constraint value in (0, 1]
        """
        if x <= 0.0 or c <= 0.0:
            return 1.0

        # g = min(1, c/x)
        ratio = c / x
        return min(1.0, max(0.0, ratio))

    @njit(fastmath=True, cache=True, parallel=True)
    def statedep_constraint_vectorized_jit(
        x: np.ndarray,
        c: np.ndarray
    ) -> np.ndarray:
        """
        JIT-compiled vectorized state-dependent constraint (with parallelization).

        Args:
            x: Queue lengths per state (1-D array)
            c: Capacities per state (1-D array)

        Returns:
            g: Constraints (same shape as x)
        """
        n = len(x)
        g = np.zeros(n)

        for i in prange(n):
            g[i] = statedep_constraint_jit(x[i], c[i])

        return g

    @njit(fastmath=True, cache=True)
    def fluid_ode_statedep_jit(
        t: float,
        x: np.ndarray,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
    ) -> np.ndarray:
        """
        JIT-compiled ODE with state-dependent (hard) capacity constraints.

        Args:
            t: Current time (unused)
            x: State vector (n,)
            W: Transition rate matrix (n, n)
            SQ: State-to-station queue mapping (n, n)
            Sa: Server capacities per state (n,)
            ALambda: External arrival rates (n,)

        Returns:
            dxdt: State derivatives (n,)
        """
        # Ensure non-negative state
        x = np.maximum(x, 0.0)

        # Compute total queue at each state's station
        sum_x_Qa = SQ @ x + 1e-14

        # Compute state-dependent constraints
        g = statedep_constraint_vectorized_jit(sum_x_Qa, Sa.flatten())

        # Compute derivative: dx/dt = W @ (x * g) + ALambda
        dxdt = W @ (x * g) + ALambda.flatten()

        return dxdt

else:
    # Fallback: import from pure Python statedep module
    from .statedep import (
        statedep_constraint as statedep_constraint_jit,
        statedep_constraint_vectorized as statedep_constraint_vectorized_jit,
        fluid_ode_statedep as fluid_ode_statedep_jit,
    )


__all__ = [
    'HAS_NUMBA',
    'statedep_constraint_jit',
    'statedep_constraint_vectorized_jit',
    'fluid_ode_statedep_jit',
]
