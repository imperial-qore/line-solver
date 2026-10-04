"""Exact MVA-style recursion for delayed-hit (list-based) cache metrics.

Native-Python port of matlab/src/api/retrieval/retrieval_mva.m.

Implements the exact recursive characterization that the FPI heuristic
(retrieval_fpi) approximates. For a list-based cache with delayed hits whose
retrieval system has one IS station (s=0) and r PS stations (s=1,...,r), with
phi^{(k)}_{s,i} the delayed-hit probability of item i at station s in the
system WITHOUT item k:

  theta_{ij}(m) = gamma_{ij} / (1 + lambda_i eta_{0,i}
                    + sum_{s=1}^r lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m-1_j)))
  xi_j(m)       = m_j / sum_i theta_{ij}(m)(1 - pihit_i(m-1_j))
  pi_{ij}(m)    = theta_{ij}(m) xi_j(m) (1 - pihit_i(m-1_j))
  pihit_i(m)    = sum_j pi_{ij}(m)
  pi_{i0}(m)    = (1 - pihit_i(m)) / (1 + lambda_i eta_{0,i}
                    + sum_{s=1}^r lambda_i eta_{s,i}(1 + sum_{k!=i} phi^{(i)}_{s,k}(m)))
  phi_{sk}(m)   = lambda_k pi_{k0}(m) eta_{s,k}(1 + sum_{i!=k} phi^{(k)}_{s,i}(m))  (s=1..r)
  phi_{0k}(m)   = lambda_k eta_{0,k} pi_{k0}(m)

The recursion terminates at the empty item set and at the empty cache
(pihit_i = 0). This is exact and agrees with retrieval_nc / retrieval_metrics;
it is memoized over (item-subset, capacity) with
O(2^n n^2 h r prod_j(1+m_j)) time, so it is feasible only for small systems
(use retrieval_fpi otherwise).
"""

import sys

import numpy as np

__all__ = ['retrieval_mva']

# The recursion descends one item and one cache slot at a time, so its depth is
# n + sum(m); guard the interpreter limit rather than the problem size, which
# the caller bounds through max_items/max_states.
_RECURSION_HEADROOM = 200


def retrieval_mva(m, lambda_, eta, gamma, max_items=20, max_states=1 << 22):
    """Exact miss, hit and delayed-hit metrics of a delayed-hit cache.

    Args:
        m: cache list capacities, length h
        lambda_: per-item arrival rates, length n
        eta: fetching demands, (n, r+1); column 0 is the IS station s=0 and
            columns 1..r are the PS stations
        gamma: access factors gamma(i,j), (n, h)
        max_items: refuse above this item count; the memo is indexed by an
            n-bit item subset, so the table doubles with every further item
        max_states: refuse above this memo size, 2^n * prod_j(1+m_j)

    Returns:
        (pmiss, phit, pdh) with pmiss the miss ratios pi_{i,0} (length n), phit
        the hit ratios pi_{i,j} of shape (h, n), and pdh the delayed-hit
        probabilities phi_{s,i} for s = 0..r, of shape (r+1, n).

    Raises:
        ValueError: if the system exceeds max_items or max_states. This is an
            exact enumerative oracle for small systems; use retrieval_fpi for
            anything larger.
    """
    lambda_ = np.asarray(lambda_, dtype=float).ravel()
    m = np.asarray(m, dtype=int).ravel()
    eta = np.atleast_2d(np.asarray(eta, dtype=float))
    gamma = np.atleast_2d(np.asarray(gamma, dtype=float))

    n = lambda_.size
    h = m.size
    r = eta.shape[1] - 1
    if eta.shape[0] != n:
        raise ValueError('retrieval_mva: eta must have one row per item.')
    if gamma.shape != (n, h):
        raise ValueError('retrieval_mva: gamma must be %d x %d.' % (n, h))
    if n > max_items:
        raise ValueError('retrieval_mva: %d items exceeds max_items=%d. The memo '
                         'is indexed by an n-bit item subset; use retrieval_fpi '
                         'for larger systems.' % (n, max_items))

    eta0 = eta[:, 0]
    eta_ps = eta[:, 1:]

    radix = m + 1
    ncap = int(np.prod(radix)) if h > 0 else 1
    nmask = 1 << n
    if nmask * ncap > max_states:
        raise ValueError('retrieval_mva: memo of %d states exceeds max_states=%d. '
                         'Use retrieval_fpi for larger systems.' % (nmask * ncap, max_states))

    done = np.zeros((nmask, ncap), dtype=bool)
    PI0 = np.zeros((nmask, ncap, n))
    PIHIT = np.zeros((nmask, ncap, n))
    PIJ = np.zeros((nmask, ncap, n, h))
    PHI = np.zeros((nmask, ncap, n, r + 1))   # station s=0..r -> page 0..r

    mul = np.ones(h, dtype=np.int64)
    for jj in range(1, h):
        mul[jj] = mul[jj - 1] * radix[jj - 1]

    def capidx(c):
        return int(np.dot(np.asarray(c, dtype=np.int64), mul))

    def solve(mask, c):
        ci_ = capidx(c)
        if done[mask, ci_]:
            return
        if mask == 0:
            done[mask, ci_] = True
            return
        active = [i for i in range(n) if (mask >> i) & 1]

        # boundary: when the cache can hold every active item, all items are
        # permanently cached -> pihit=1, pi0=0, phi=0 (no fetching)
        if int(np.sum(c)) >= len(active):
            for i in active:
                PIHIT[mask, ci_, i] = 1.0
            done[mask, ci_] = True
            return

        # --- ensure dependencies are solved ---
        for jj in range(h):
            if c[jj] > 0:
                cj = np.array(c, dtype=int)
                cj[jj] -= 1
                solve(mask, cj)
                for i in active:
                    solve(mask & ~(1 << i), cj)
        for k in active:
            solve(mask & ~(1 << k), c)

        # --- theta, xi, pi_{ij} (cache recursion on m-1_j) ---
        for jj in range(h):
            if c[jj] > 0:
                cj = np.array(c, dtype=int)
                cj[jj] -= 1
                cjx = capidx(cj)
                theta = np.zeros(n)
                for i in active:
                    maski = mask & ~(1 << i)
                    acc = 0.0
                    for s in range(r):
                        sphi = 0.0
                        for k in active:
                            if k != i:
                                sphi += PHI[maski, cjx, k, s + 1]
                        acc += lambda_[i] * eta_ps[i, s] * (1.0 + sphi)
                    theta[i] = gamma[i, jj] / (1.0 + lambda_[i] * eta0[i] + acc)
                sden = 0.0
                for i in active:
                    sden += theta[i] * (1.0 - PIHIT[mask, cjx, i])
                xi_j = c[jj] / sden
                for i in active:
                    PIJ[mask, ci_, i, jj] = theta[i] * xi_j * (1.0 - PIHIT[mask, cjx, i])

        # --- pihit_i = sum_j pi_{ij} ---
        for i in active:
            PIHIT[mask, ci_, i] = float(np.sum(PIJ[mask, ci_, i, :]))

        # --- pi_{i0} (uses phi^{(i)}(m) at same capacity) ---
        for i in active:
            maski = mask & ~(1 << i)
            acc = 0.0
            for s in range(r):
                sphi = 0.0
                for k in active:
                    if k != i:
                        sphi += PHI[maski, ci_, k, s + 1]
                acc += lambda_[i] * eta_ps[i, s] * (1.0 + sphi)
            PI0[mask, ci_, i] = ((1.0 - PIHIT[mask, ci_, i])
                                 / (1.0 + lambda_[i] * eta0[i] + acc))

        # --- phi_{sk}(m) ---
        for k in active:
            maskk = mask & ~(1 << k)
            pi0k = PI0[mask, ci_, k]
            PHI[mask, ci_, k, 0] = lambda_[k] * eta0[k] * pi0k   # s = 0 (IS)
            for s in range(r):
                sphi = 0.0
                for i in active:
                    if i != k:
                        sphi += PHI[maskk, ci_, i, s + 1]
                PHI[mask, ci_, k, s + 1] = (lambda_[k] * pi0k * eta_ps[k, s]
                                            * (1.0 + sphi))

        done[mask, ci_] = True

    fullmask = nmask - 1
    depth_needed = n + int(np.sum(m)) + _RECURSION_HEADROOM
    old_limit = sys.getrecursionlimit()
    if old_limit < depth_needed:
        sys.setrecursionlimit(depth_needed)
    try:
        solve(fullmask, m)
    finally:
        sys.setrecursionlimit(old_limit)

    ci = capidx(m)
    pmiss = np.array([PI0[fullmask, ci, i] for i in range(n)])
    phit = np.array([[PIJ[fullmask, ci, i, j] for i in range(n)] for j in range(h)])
    pdh = np.array([[PHI[fullmask, ci, i, s] for i in range(n)] for s in range(r + 1)])
    return pmiss, phit, pdh
