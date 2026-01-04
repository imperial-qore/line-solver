"""
Softmin smoothing ODE system for FLD solver.

Implements the softmin smoothed capacity constraint:

    softmin(x, c, alpha) = log(sum(exp(-alpha * [x, c]))) / (-alpha)

As alpha → ∞, softmin(x,c,alpha) → min(x, c), recovering the
original hard capacity constraint.

The softmin function provides smooth transitions without the
discontinuities present in the hard constraint, improving
numerical stability of ODE integration.
"""

import numpy as np
from typing import Optional


def softmin_smooth(x: float, c: float, alpha: float = 20.0) -> float:
    """
    Softmin smoothed min function.

    Approximates min(x, c) using a smooth function that avoids
    discontinuities in the ODE system.

    Args:
        x: Total queue length at station (or phase)
        c: Server capacity (number of servers * service rate)
        alpha: Softmin parameter (higher alpha → sharper transition)

    Returns:
        softmin: Smooth approximation to min(x, c) / x in (0, 1]
    """
    if c <= 0 or alpha <= 0 or x <= 0:
        return 1.0

    try:
        # softmin(x,c,alpha) = log(exp(-alpha*x) + exp(-alpha*c)) / (-alpha)
        # But we want a capacity constraint g = softmin / x, so:
        # g = softmin(x, c, alpha) / x

        # Numerical stability: use log-sum-exp trick
        # log(a + b) = log(a) + log(1 + b/a) when a > b
        exp_x = np.exp(-alpha * x)
        exp_c = np.exp(-alpha * c)

        if exp_x >= exp_c:
            # x <= c case: log(exp(-alpha*x) + exp(-alpha*c))
            #             = -alpha*x + log(1 + exp(-alpha*(c-x)))
            log_sum = -alpha * x + np.log1p(np.exp(-alpha * (c - x)))
        else:
            # x > c case: log(exp(-alpha*x) + exp(-alpha*c))
            #           = -alpha*c + log(1 + exp(-alpha*(x-c)))
            log_sum = -alpha * c + np.log1p(np.exp(-alpha * (x - c)))

        # softmin = log_sum / (-alpha)
        softmin_val = log_sum / (-alpha)

        # g = softmin / x (normalized capacity constraint)
        g = softmin_val / x if x > 0 else 1.0

        # Ensure result is valid and bounded
        if np.isnan(g) or not np.isfinite(g):
            return 0.0
        if g < 0:
            return 0.0
        if g > 1.0:
            return 1.0

        return float(g)

    except (ValueError, FloatingPointError):
        return 0.0


def softmin_smooth_vectorized(x: np.ndarray, c: np.ndarray, alpha: float = 20.0) -> np.ndarray:
    """
    Vectorized softmin smoothed constraint computation.

    Args:
        x: Queue lengths per state (1-D array)
        c: Capacities per state (1-D array)
        alpha: Softmin parameter (scalar)

    Returns:
        g: Smooth constraints (same shape as x)
    """
    x = np.asarray(x).flatten()
    c = np.asarray(c).flatten()

    # Ensure arrays have same length
    n = len(x)
    c = np.pad(c, (0, n - len(c)), mode='edge') if len(c) < n else c[:n]

    g = np.zeros(n)

    for i in range(n):
        g[i] = softmin_smooth(x[i], c[i], alpha)

    return g


def fluid_ode_softmin(
    t: float,
    x: np.ndarray,
    W: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    ALambda: np.ndarray,
    alpha: float = 20.0,
) -> np.ndarray:
    """
    ODE right-hand side with softmin smoothed capacity constraints.

    dx/dt = W * (x .* g(x)) + ALambda

    where g(x) is the softmin-smoothed capacity constraint.

    Args:
        t: Current time (unused but required by scipy.integrate)
        x: State vector - queue lengths per phase (n,)
        W: Transition rate matrix (n, n)
        SQ: State-to-station queue mapping (n, n)
        Sa: Server capacities per state (n,)
        ALambda: External arrival rates (n,) or (n, 1)
        alpha: Softmin parameter (default 20.0)

    Returns:
        dxdt: State derivatives (n,)
    """
    # Ensure non-negative state
    x = np.maximum(x, 0.0)

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-14  # Add small value for numerical stability

    # Compute softmin smoothed capacity constraints
    g = softmin_smooth_vectorized(sum_x_Qa, Sa.flatten(), alpha)

    # Compute derivative: dx/dt = W @ (x * g) + ALambda
    dxdt = W @ (x * g) + np.asarray(ALambda).flatten()

    return dxdt


class SoftminODESystem:
    """Softmin smoothed ODE system for fluid analysis."""

    def __init__(
        self,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
        alpha: float = 20.0,
    ):
        """Initialize softmin ODE system.

        Args:
            W: Transition rate matrix
            SQ: State-to-station queue mapping
            Sa: Server capacities
            ALambda: External arrival rates
            alpha: Softmin parameter (default 20.0)
        """
        self.W = W
        self.SQ = SQ
        self.Sa = Sa
        self.ALambda = np.asarray(ALambda).flatten()
        self.alpha = float(alpha)

    def __call__(self, t: float, x: np.ndarray) -> np.ndarray:
        """Evaluate ODE right-hand side.

        Args:
            t: Time
            x: State vector

        Returns:
            State derivatives
        """
        return fluid_ode_softmin(t, x, self.W, self.SQ, self.Sa, self.ALambda, self.alpha)

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
