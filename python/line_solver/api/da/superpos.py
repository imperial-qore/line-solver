"""
Asymptotic-method superposition of independent flows (Whitt's QNA
stationary-interval formula).
"""

import numpy as np

__all__ = ['da_traffic_superpos']


def da_traffic_superpos(lambd, a2):
    """
    Rate-weighted SCV mixture of merged independent flows.

    Args:
        lambd: flow rates (array-like); non-finite entries are ignored.
        a2: squared coefficients of variation of each flow.

    Returns:
        The SCV of the superposed flow.
    """
    lambd = np.asarray(lambd, dtype=float).ravel()
    a2 = np.asarray(a2, dtype=float).ravel()
    mask = np.isfinite(lambd)
    a2 = a2[mask]
    lambd = lambd[mask]
    return float(np.dot(a2, lambd) / np.sum(lambd))
