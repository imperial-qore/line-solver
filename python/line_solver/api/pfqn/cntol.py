"""Chandy-Neuse population-scaled termination cutoff for approximate MVA."""

import numpy as np

__all__ = ['pfqn_cntol', 'is_cntol']


def pfqn_cntol(N) -> float:
    r"""
    Termination cutoff 1/(4000 + 16*sum(N)) of the Linearizer.

    Published in K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for
    Queuing Network Models of Computing Systems", Commun. ACM 25(2):126-134,
    1982, p.129 and appendix. The iteration continues while

        max_{i,r} \|Q^I(i,r) - Q^{I-1}(i,r)\| / N_r > 1/(4000 + 16*\|N\|),

    \|N\| = sum(N). The paper motivates the scaling with \|N\|: at large populations
    removing one job changes the queue lengths very little, so a fixed cutoff
    would terminate the iteration prematurely. It also notes that the expression
    stays below 0.00025 even at very small populations.

    The same expression is what LQNS uses as its termination test, set in the
    SchweitzerCommon constructor of libmva/src/mva.cc; that code carries no
    citation, and the paper above is its source.

    Passing tol='cn' (or NaN) to pfqn_bs / pfqn_egflinearizer selects BOTH this
    cutoff and the normalized-maximum metric of the paper, which is the
    published test; passing pfqn_cntol(N) as a plain number selects only the
    cutoff, with those functions' own convergence metric.

    Args:
        N: Population vector

    Returns:
        The cutoff as a float
    """
    N = np.asarray(N, dtype=float).ravel()
    return 1.0 / (4000.0 + 16.0 * float(N.sum()))


def is_cntol(tol) -> bool:
    """
    True when tol requests the Chandy-Neuse termination test.

    The sentinel is the string 'cn' or NaN; NaN is the form carried across the
    MATLAB, Java and C++ twins, which have no string tolerance.
    """
    if isinstance(tol, str):
        if tol.lower() != 'cn':
            raise ValueError("unknown tolerance specifier '%s'" % tol)
        return True
    return tol is not None and np.isnan(tol)
