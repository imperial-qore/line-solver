"""
P-norm smoothing ODE system for FLD solver.

Implements the p-norm smoothed capacity constraint as per
Ruuskanen et al., PEVA 151 (2021):

    dx/dt = W * (x .* g_hat(x)) + ALambda

where g_hat(x) = 1 / (1 + (x/c)^p)^(1/p)

As p → ∞, g_hat → min(1, c/x), recovering the original
hard capacity constraint.
"""

import numpy as np
from typing import Optional, Tuple


def pnorm_smooth(x: float, c: float, p: float = 20.0) -> float:
    """
    P-norm smoothed capacity constraint function.

    Approximates min(1, c/x) using a smooth function that avoids
    discontinuities in the ODE system, improving numerical stability.

    Args:
        x: Total queue length at station (or phase)
        c: Server capacity (number of servers * service rate)
        p: P-norm parameter (higher p → sharper transition)

    Returns:
        g_hat: Smoothed constraint value in (0, 1]
    """
    if c <= 0 or p <= 0:
        return 1.0

    if x <= 0:
        return 1.0

    # g_hat = 1 / (1 + (x/c)^p)^(1/p)
    try:
        ratio = x / c
        power_term = np.power(ratio, p)
        denominator = np.power(1 + power_term, 1.0 / p)
        ghat = 1.0 / denominator

        # Ensure result is valid
        if np.isnan(ghat) or not np.isfinite(ghat):
            return 0.0
        if ghat < 0:
            return 0.0
        if ghat > 1.0:
            return 1.0

        return float(ghat)
    except (ValueError, FloatingPointError):
        return 0.0


def pnorm_smooth_vectorized(x: np.ndarray, c: np.ndarray, p: np.ndarray = None) -> np.ndarray:
    """
    Vectorized p-norm smoothed constraint computation.

    Args:
        x: Queue lengths per state (1-D array)
        c: Capacities per state (1-D array)
        p: P-norm parameters per state (None = use 20.0 default)

    Returns:
        g_hat: Smoothed constraints (same shape as x)
    """
    x = np.asarray(x).flatten()
    c = np.asarray(c).flatten()

    if p is None:
        p = np.full_like(x, 20.0, dtype=float)
    else:
        p = np.asarray(p).flatten()

    # Ensure all arrays have same length
    n = len(x)
    c = np.pad(c, (0, n - len(c)), mode='edge') if len(c) < n else c[:n]
    p = np.pad(p, (0, n - len(p)), mode='edge') if len(p) < n else p[:n]

    ghat = np.zeros(n)

    for i in range(n):
        ghat[i] = pnorm_smooth(x[i], c[i], p[i])

    return ghat


def fluid_ode_pnorm(
    t: float,
    x: np.ndarray,
    W: np.ndarray,
    SQ: np.ndarray,
    Sa: np.ndarray,
    ALambda: np.ndarray,
    pstar: Optional[np.ndarray] = None,
) -> np.ndarray:
    """
    ODE right-hand side with p-norm smoothed capacity constraints.

    dx/dt = W * (x .* g_hat(x)) + ALambda

    Args:
        t: Current time (unused but required by scipy.integrate)
        x: State vector - queue lengths per phase (n,)
        W: Transition rate matrix (n, n)
        SQ: State-to-station queue mapping (n, n)
        Sa: Server capacities per state (n,)
        ALambda: External arrival rates (n,) or (n, 1)
        pstar: P-norm parameters per state (n,) or scalar

    Returns:
        dxdt: State derivatives (n,)
    """
    # Ensure non-negative state
    x = np.maximum(x, 0.0)

    # Compute total queue at each state's station
    sum_x_Qa = SQ @ x + 1e-14  # Add small value for numerical stability

    # Compute smoothed capacity constraints
    if pstar is not None:
        if np.isscalar(pstar):
            pstar_vals = np.full(len(x), pstar, dtype=float)
        else:
            pstar_vals = np.asarray(pstar).flatten()
            if len(pstar_vals) == 1:
                pstar_vals = np.full(len(x), pstar_vals[0], dtype=float)
            else:
                # Pad or truncate to match state dimension
                pstar_vals = np.pad(pstar_vals, (0, len(x) - len(pstar_vals)), mode='edge')[:len(x)]

        ghat = pnorm_smooth_vectorized(sum_x_Qa, Sa.flatten(), pstar_vals)
    else:
        # Default p-norm with p=20
        ghat = pnorm_smooth_vectorized(sum_x_Qa, Sa.flatten(), np.full(len(x), 20.0))

    # Compute derivative: dx/dt = W @ (x * ghat) + ALambda
    dxdt = W @ (x * ghat) + np.asarray(ALambda).flatten()

    return dxdt


class PNormODESystem:
    """P-norm smoothed ODE system for fluid analysis."""

    def __init__(
        self,
        W: np.ndarray,
        SQ: np.ndarray,
        Sa: np.ndarray,
        ALambda: np.ndarray,
        pstar: Optional[np.ndarray] = None,
    ):
        """Initialize p-norm ODE system.

        Args:
            W: Transition rate matrix
            SQ: State-to-station queue mapping
            Sa: Server capacities
            ALambda: External arrival rates
            pstar: P-norm parameters (default: 20.0)
        """
        self.W = W
        self.SQ = SQ
        self.Sa = Sa
        self.ALambda = np.asarray(ALambda).flatten()

        if pstar is None:
            self.pstar = np.full(len(self.Sa), 20.0, dtype=float)
        elif np.isscalar(pstar):
            self.pstar = np.full(len(self.Sa), float(pstar), dtype=float)
        else:
            self.pstar = np.asarray(pstar).flatten()

    def __call__(self, t: float, x: np.ndarray) -> np.ndarray:
        """Evaluate ODE right-hand side.

        Args:
            t: Time
            x: State vector

        Returns:
            State derivatives
        """
        return fluid_ode_pnorm(t, x, self.W, self.SQ, self.Sa, self.ALambda, self.pstar)

    def steady_state_constraint(self) -> np.ndarray:
        """
        At steady state, dx/dt = 0, so:

        W @ (x* .* g_hat(x*)) + ALambda = 0

        Returns:
            Implicit residual function for solver
        """
        def residual(x):
            return self(0, x)
        return residual
