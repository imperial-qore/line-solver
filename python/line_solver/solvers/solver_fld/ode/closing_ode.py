"""
Closing method ODE system for FLD solver.

Implements the ODE system used in the closing method approximation:

    dx/dt = W' * theta + A * lambda

where:
- W' is the transposed transition matrix
- theta is the effective service rate fraction (with smoothing)
- A is the initial phase probability matrix
- lambda is the external arrival rate vector

The closing method uses iterative refinement of service processes,
solving this ODE at each iteration until convergence.

References:
- Closed queueing network analysis via fluid approximation
- Iterative service process refinement
"""

import numpy as np
from typing import Optional, Callable


def compute_theta_pnorm(
    x: np.ndarray,
    c: np.ndarray,
    SQ: np.ndarray,
    p: float = 20.0,
) -> np.ndarray:
    """
    Compute effective service rate fraction using p-norm smoothing.

    theta[i] = min(c[i], sum_station(x[i])) / sum_station(x[i])

    with p-norm smoothing applied to the min function.

    Args:
        x: State vector (queue lengths per phase)
        c: Server capacities per state
        SQ: State-to-station aggregation matrix
        p: P-norm parameter (default 20.0)

    Returns:
        theta: Effective service rates per state
    """
    x = np.maximum(x, 0.0)
    sum_x_Qa = SQ @ x + 1e-14

    # P-norm smoothing: ghat = 1 / (1 + (x/c)^p)^(1/p)
    theta = np.zeros_like(x)
    for i in range(len(x)):
        x_val = sum_x_Qa[i]
        c_val = c[i]
        p_val = p if p > 0 else 1.0

        if c_val > 0 and x_val > 0:
            ratio = x_val / c_val
            ghat = 1.0 / np.power(1 + np.power(ratio, p_val), 1.0 / p_val)
            if np.isnan(ghat) or not np.isfinite(ghat):
                ghat = 0.0
            theta[i] = ghat
        else:
            theta[i] = 1.0

    return theta


def compute_theta_softmin(
    x: np.ndarray,
    c: np.ndarray,
    SQ: np.ndarray,
    alpha: float = 20.0,
) -> np.ndarray:
    """
    Compute effective service rate fraction using softmin smoothing.

    Uses log-sum-exp trick for numerical stability.

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
    for i in range(len(x)):
        x_val = sum_x_Qa[i]
        c_val = c[i]

        if c_val > 0 and x_val > 0:
            # Softmin using log-sum-exp
            exp_x = np.exp(-alpha * x_val)
            exp_c = np.exp(-alpha * c_val)

            if exp_x >= exp_c:
                log_sum = -alpha * x_val + np.log1p(np.exp(-alpha * (c_val - x_val)))
            else:
                log_sum = -alpha * c_val + np.log1p(np.exp(-alpha * (x_val - c_val)))

            softmin_val = log_sum / (-alpha)
            theta[i] = softmin_val / x_val if x_val > 0 else 1.0

            if np.isnan(theta[i]) or not np.isfinite(theta[i]):
                theta[i] = 0.0
        else:
            theta[i] = 1.0

    return theta


def compute_theta_statedep(
    x: np.ndarray,
    c: np.ndarray,
    SQ: np.ndarray,
) -> np.ndarray:
    """
    Compute effective service rate fraction using state-dependent constraint.

    Hard constraint: theta = min(c, x) / x

    Args:
        x: State vector
        c: Server capacities per state
        SQ: State-to-station aggregation matrix

    Returns:
        theta: Effective service rates
    """
    x = np.maximum(x, 0.0)
    sum_x_Qa = SQ @ x + 1e-14

    # theta = min(1, c / sum_x_Qa)
    with np.errstate(divide='ignore', invalid='ignore'):
        theta = np.minimum(1.0, np.divide(c, sum_x_Qa))
        theta = np.where(np.isfinite(theta), theta, 1.0)
        theta = np.where(theta > 0, theta, 0.0)

    return theta


def fluid_ode_closing(
    t: float,
    x: np.ndarray,
    W_T: np.ndarray,
    SQ: np.ndarray,
    c: np.ndarray,
    A_Lambda: np.ndarray,
    theta_func: Callable,
    **theta_kwargs,
) -> np.ndarray:
    """
    ODE right-hand side for closing method.

    dx/dt = W' * (theta .* A) + ALambda

    Args:
        t: Current time (unused)
        x: State vector
        W_T: Transposed transition matrix
        SQ: State-to-station aggregation matrix
        c: Server capacities
        A_Lambda: Arrival rates
        theta_func: Function to compute effective service rates
        **theta_kwargs: Keyword arguments for theta_func

    Returns:
        dxdt: State derivatives
    """
    x = np.maximum(x, 0.0)

    # Compute effective service rates with smoothing
    theta = theta_func(x, c, SQ, **theta_kwargs)

    # Compute derivative: dx/dt = W' @ (x * theta) + A*Lambda
    dxdt = W_T @ (x * theta) + np.asarray(A_Lambda).flatten()

    return dxdt


class ClosingODESystem:
    """Closing method ODE system with selectable smoothing strategy."""

    def __init__(
        self,
        W_T: np.ndarray,
        SQ: np.ndarray,
        c: np.ndarray,
        A_Lambda: np.ndarray,
        method: str = 'pnorm',
        **params,
    ):
        """Initialize closing ODE system.

        Args:
            W_T: Transposed transition matrix
            SQ: State-to-station aggregation matrix
            c: Server capacities
            A_Lambda: Arrival rates
            method: Smoothing method ('pnorm', 'softmin', 'statedep')
            **params: Method-specific parameters
        """
        self.W_T = W_T
        self.SQ = SQ
        self.c = c
        self.A_Lambda = np.asarray(A_Lambda).flatten()
        self.method = method

        # Select theta computation function
        if method == 'pnorm':
            self.theta_func = compute_theta_pnorm
            self.theta_kwargs = {'p': params.get('p', 20.0)}
        elif method == 'softmin':
            self.theta_func = compute_theta_softmin
            self.theta_kwargs = {'alpha': params.get('alpha', 20.0)}
        elif method == 'statedep':
            self.theta_func = compute_theta_statedep
            self.theta_kwargs = {}
        else:
            raise ValueError(f"Unknown method: {method}")

    def __call__(self, t: float, x: np.ndarray) -> np.ndarray:
        """Evaluate ODE right-hand side.

        Args:
            t: Time
            x: State vector

        Returns:
            State derivatives
        """
        return fluid_ode_closing(
            t, x, self.W_T, self.SQ, self.c, self.A_Lambda,
            self.theta_func, **self.theta_kwargs
        )
