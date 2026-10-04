"""
State-dependent ODE system for FLD solver (no smoothing).

Implements the hard capacity constraint:

    g(x) = min(c, x) / x

This is the direct (non-smoothed) capacity constraint that recovers
the minimum of the number of jobs and available capacity.

Note: This formulation may have discontinuities that can impact
numerical stability. For production use, consider using p-norm or
softmin smoothing instead.
"""

import numpy as np
from typing import Optional


def statedep_constraint(x: float, c: float) -> float:
    """
    State-dependent (hard) capacity constraint.

    Returns min(1, c/x), representing the fraction of server
    capacity that is utilized when x jobs are present.

    Args:
        x: Total queue length at station (or phase)
        c: Server capacity (number of servers * service rate)

    Returns:
        g: Constraint value in (0, 1]
    """
    if x <= 0 or c <= 0:
        return 1.0

    # g = min(c, x) / x = min(1, c/x)
    ratio = c / x

    if ratio >= 1.0:
        return 1.0
    elif ratio > 0:
        return float(ratio)
    else:
        return 0.0


def statedep_constraint_vectorized(x: np.ndarray, c: np.ndarray) -> np.ndarray:
    """
    Vectorized state-dependent constraint computation.

    Args:
        x: Queue lengths per state (1-D array)
        c: Capacities per state (1-D array)

    Returns:
        g: Constraints (same shape as x)
    """
    x = np.asarray(x).flatten()
    c = np.asarray(c).flatten()

    # Ensure arrays have same length
    n = len(x)
    c = np.pad(c, (0, n - len(c)), mode='edge') if len(c) < n else c[:n]

    # g = min(1, c / x)
    with np.errstate(divide='ignore', invalid='ignore'):
        g = np.minimum(1.0, np.divide(c, x))
        g = np.where(np.isfinite(g), g, 1.0)
        g = np.where(g > 0, g, 0.0)

    return g


def fluid_ode_statedep(
    t: float,
    x: np.ndarray,
    W: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    ALambda: np.ndarray,
) -> np.ndarray:
    """
    ODE right-hand side with state-dependent (hard) capacity constraints.

    dx/dt = W * (x .* g(x)) + ALambda

    where g(x) = min(1, c/x) is the hard capacity constraint.

    Args:
        t: Current time (unused but required by scipy.integrate)
        x: State vector - queue lengths per phase (n,)
        W: Transition rate matrix (n, n)
        SQ: State-to-station queue mapping (n, n)
        Sa: Server capacities per state (n,)
        ALambda: External arrival rates (n,) or (n, 1)

    Returns:
        dxdt: State derivatives (n,)
    """
    # Ensure non-negative state
    x = np.maximum(x, 0.0)

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-14  # Add small value for numerical stability

    # Compute hard capacity constraints
    g = statedep_constraint_vectorized(sum_x_Qa, Sa.flatten())

    # Compute derivative: dx/dt = W @ (x * g) + ALambda
    dxdt = W @ (x * g) + np.asarray(ALambda).flatten()

    return dxdt


class StateDepODESystem:
    """State-dependent (hard constraint) ODE system for fluid analysis."""

    def __init__(
        self,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
    ):
        """Initialize state-dependent ODE system.

        Args:
            W: Transition rate matrix
            SQ: State-to-station queue mapping
            Sa: Server capacities
            ALambda: External arrival rates
        """
        self.W = W
        self.SQ = SQ
        self.Sa = Sa
        self.ALambda = np.asarray(ALambda).flatten()

    def __call__(self, t: float, x: np.ndarray) -> np.ndarray:
        """Evaluate ODE right-hand side.

        Args:
            t: Time
            x: State vector

        Returns:
            State derivatives
        """
        return fluid_ode_statedep(t, x, self.W, self.SQ, self.Sa, self.ALambda)

    def steady_state_constraint(self):
        """
        At steady state, dx/dt = 0, so:

        W @ (x* .* g(x*)) + ALambda = 0

        Returns:
            Implicit residual function for solver
        """
        def residual(x):
            return self(0, x)
        return residual
