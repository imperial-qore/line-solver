"""
Exact analysis of a loss network by MDD-rec: the normalising constant is the sum
of a product form over the admissible set {n >= 0 : A n <= C}, which is what a
decision diagram holding that set computes in one memoised walk.

A Kelly loss network carries offered load nu_r on route r and admits a call only
while the resource constraint A n <= C still holds after it. The stationary law
is the truncation of independent Poisson counts to that set,

    P(n) = (1/G) prod_r nu_r^{n_r} / n_r!,   G = sum_{A n <= C} prod_r ...,

so g_r(k) = nu_r^k/k! and mdd_rec returns G. By PASTA the acceptance probability
of a class-r call is the ratio of two such constants,

    1 - B_r = G(C - A e_r) / G(C),

which is one further diagram per class.

WHY THIS EXISTS ALONGSIDE lossn_manjunath. The Manjunath-Sikdar transform
evaluates G exactly as a multidimensional residue, and the residue argument
counts WHOLE UNITS: it needs an integral A and C. On a region declaring a
fractional class size or capacity the analyzer had no exact route at all and fell
back to the Erlang fixed point, an approximation. MDD-rec needs only that the
admissible set be finite and bounded coordinate by coordinate, which a fractional
constraint still is, so it is exact there too. It is also an exact alternative to
the Monte Carlo summation lossn_mci estimates.

References
----------
F. P. Kelly, "Loss networks", Annals of Applied Probability 1(3), 1991.
S. Balsamo, A. Marin, I. Stojic, "Computation of the normalising constant for
product-form models of distributed systems with synchronisation", Future
Generation Computer Systems 111 (2020) 475-490.

See also: lossn_manjunath, lossn_erlangfp, lossn_mci, mdd_rec.
"""

from math import factorial, log, exp
from typing import List, Tuple

import numpy as np

from ..io.logging import line_error
from ..mdd import mdd_reachset, mdd_rec, mdd_rec_marginal

__all__ = ['lossn_rec']


def lossn_rec(nu, A, C) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """Exact loss-network analysis by MDD-rec.

    Parameters
    ----------
    nu : offered load per class, length K
    A  : J x K non-negative resource requirement matrix
    C  : capacity vector, length J

    Returns
    -------
    (QLen, Loss, lG, niter) with QLen the carried load per class, Loss the
    blocking probability per class, lG the log normalising constant G(C) and
    niter the number of diagram walks performed, K + 1.
    """
    nu = np.ravel(np.asarray(nu, dtype=float))
    A = np.atleast_2d(np.asarray(A, dtype=float))
    C = np.ravel(np.asarray(C, dtype=float))
    K = nu.size
    if A.shape[1] != K:
        line_error('lossn_rec', 'A has %d columns but there are %d classes' % (A.shape[1], K))
    if A.shape[0] != C.size:
        line_error('lossn_rec', 'A has %d rows but C has %d entries' % (A.shape[0], C.size))
    if np.any(A < 0):
        line_error('lossn_rec', 'the resource matrix A must be non-negative')

    # ---- per-class bound: the most calls the tightest constraint alone admits
    bound = np.zeros(K, dtype=int)
    for r in range(K):
        j = np.nonzero(A[:, r] > 0)[0]
        if j.size == 0:
            line_error('lossn_rec',
                       'class %d consumes no resource, so the admissible set is unbounded in '
                       'that coordinate and its normalising constant diverges' % (r + 1))
        bound[r] = max(0, int(np.floor(np.min(C[j] / A[j, r]))))

    g: List[np.ndarray] = []
    for r in range(K):
        k = np.arange(bound[r] + 1)
        g.append((nu[r] ** k) / np.array([float(factorial(int(v))) for v in k]))

    lG = _log_g(A, C, bound, g)
    if not np.isfinite(lG):
        line_error('lossn_rec', 'the admissible set is empty: no call of any class fits within C')

    # ---- carried load per class, from the marginals of the same diagram
    mdds = _diagram(A, C, bound)
    G = exp(lG)
    QLen = np.zeros(K)
    for r in range(K):
        pk = np.asarray(mdd_rec_marginal(mdds, g, r)) / G
        QLen[r] = float(np.arange(pk.size) @ pk)

    # ---- blocking: 1 - B_r = G(C - A e_r)/G(C), Kelly's ratio, by PASTA
    Loss = np.zeros(K)
    for r in range(K):
        Cr = C - A[:, r]
        if np.any(Cr < 0):
            Loss[r] = 1.0                        # the call never fits
            continue
        lGr = _log_g(A, Cr, bound, g)
        Loss[r] = 1.0 if not np.isfinite(lGr) else 1.0 - exp(lGr - lG)
        Loss[r] = min(1.0, max(0.0, Loss[r]))

    return QLen, Loss, float(lG), K + 1


def _diagram(A, C, bound):
    """The admissible set {n >= 0 : A n <= C}, generated one call at a time from
    the empty network. Adding a call is the only move, so the breadth-first
    closure visits exactly the admissible vectors."""
    K = bound.size
    domain = (bound + 1).astype(int)

    def nextfun(s):
        s = np.asarray(s, dtype=float)
        out = []
        for r in range(K):
            if s[r] >= bound[r]:
                continue
            t = s.copy()
            t[r] += 1
            if np.all(A @ t <= C + 1e-12):
                out.append(tuple(int(v) for v in t))
        return out

    mdd = mdd_reachset(domain, np.zeros(K, dtype=int), nextfun)
    return mdd.to_struct()


def _log_g(A, C, bound, g) -> float:
    """log G over the admissible set at capacity C, keeping the per-class domains
    of the FULL problem so that one set of factors g serves every reduced
    capacity."""
    if np.any(C < 0):
        return -np.inf
    G = mdd_rec(_diagram(A, C, bound), g)
    return log(G) if G > 0 else -np.inf
