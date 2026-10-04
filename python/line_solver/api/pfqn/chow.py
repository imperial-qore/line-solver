"""JMT-compatible Chow approximate MVA."""

from typing import Tuple

import numpy as np

from .lcp import pfqn_lcp

__all__ = ['pfqn_chow']


def pfqn_chow(L, N, Z=None, tol: float = 1e-6, maxiter: int = 1000,
              QN0=None, type_sched=None, variant: str = 'forward'
              ) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, int]:
    """Run the JMT-compatible Chow approximation.

    JMT estimates every class's arrival-instant queue length at a station by
    the aggregate queue length at the full population. This is the Bard
    large-customer-population fixed point implemented by :func:`pfqn_lcp`.
    ``variant`` is retained for API compatibility and is ignored.

    Returns (XN, QN, UN, RN, it).
    """
    return pfqn_lcp(L, N, Z, tol, maxiter, QN0, type_sched)
