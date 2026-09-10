"""
Generic damped successive-substitution driver for decomposition-aggregation
(DA) fixed-point iterations.
"""

import numpy as np

__all__ = ['da_fpi']


def da_fpi(iterfun, x0, iter_max, iter_tol, damping=1.0, norm=None,
           nanstop=False, miniter=1):
    """
    Drive a DA fixed-point iteration by damped successive substitution.

    Each call to iterfun performs one DA sweep: solve the isolated submodels
    given the current coupling iterate x, exchange flows or rates, and return
    the updated iterate.

    Args:
        iterfun: callable (x, it) -> (xnew, xref) evaluating one DA sweep from
            iterate x at iteration count it. xref is the baseline for the
            convergence test; return xref = x for a standard successive-
            substitution test, or a mid-sweep checkpoint when the method
            compares against a renormalized iterate.
        x0: initial iterate (ndarray, or any object when a custom norm is
            supplied and damping is 1).
        iter_max: maximum number of sweeps.
        iter_tol: convergence tolerance on the norm of the iterate change.
        damping: damping factor in (0, 1]; 1 (default) is undamped.
        norm: callable (xnew, xref) -> scalar convergence measure; defaults
            to the max absolute elementwise difference.
        nanstop: when True, a NaN convergence measure terminates the
            iteration, replicating legacy while-loop drivers whose
            "continue while delta > tol" test exits on NaN; when False
            (default) a NaN measure keeps iterating, as in legacy
            "break if delta < tol" drivers.
        miniter: convergence is not tested before this sweep count.

    Returns:
        Tuple (x, it, converged) with the final iterate, the number of sweeps
        executed, and False if the iteration stopped at iter_max.
    """
    if norm is None:
        norm = lambda xn, xr: float(np.max(np.abs(np.asarray(xn, dtype=float) - np.asarray(xr, dtype=float))))

    x = x0
    it = 0
    converged = False
    for it in range(1, int(iter_max) + 1):
        xnew, xref = iterfun(x, it)
        if damping != 1.0:
            xnew = (1.0 - damping) * xref + damping * xnew
        delta = norm(xnew, xref)
        x = xnew
        if it >= miniter:
            if delta < iter_tol:
                converged = True
                break
            elif nanstop and np.isnan(delta):
                break
    return x, it, converged
