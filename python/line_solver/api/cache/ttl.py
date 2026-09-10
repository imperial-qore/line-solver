"""
TTL-Based Cache Analysis Methods.

Native Python implementations of Time-To-Live (TTL) based cache analysis
methods, including LRU(m) and hierarchical LRU caches.

Key functions:
    cache_t_hlru: Characteristic times for hierarchical LRU cache
    cache_ttl_hlru: Steady-state probabilities for TTL hierarchical LRU
    cache_ttl_lrua: Steady-state probabilities with arrival-based routing

References:
    Original MATLAB: matlab/src/api/cache/cache_ttl*.m, cache_t*.m
"""

import numpy as np
from scipy.optimize import fsolve, least_squares
from typing import Tuple, Optional


def cache_t_hlru(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Characteristic time of each list of an h-LRU / LRU(m) cache.

    Solves the TTL (characteristic-time) fixed point of the list-based h-LRU
    (LRU(m)) policy: sum_k pi_l(k;T) = m[l] for each list l, where the level
    probabilities follow the birth-death form pi_l ~ prod_{s<=l} (1-e_s)/e_s
    with e_s = exp(-gamma_k*T_s) (Gast and Van Houdt, SIGMETRICS 2015).
    Solved by per-list bisection with Gauss-Seidel sweeps.

    Args:
        gamma: (n,) per-item request rates; an (n x h) matrix is accepted for
               backward compatibility (first column used)
        m: Cache capacity vector (h,)

    Returns:
        Characteristic time for each cache list (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_t_hlru.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    lam = gamma[:, 0] if gamma.ndim > 1 else gamma
    h = len(m)

    def occ_l(lam_, t_, l_, tl_):
        t2 = t_.copy()
        t2[l_] = tl_
        P = _hlru_levelprobs(lam_, t2, h)
        return float(np.sum(P[:, 1 + l_]))

    t = np.ones(h) / max(float(np.mean(lam)), 1e-14)
    for _sweep in range(200):
        told = t.copy()
        for l in range(h):
            lo = 0.0
            hi = max(t[l], 1.0 / max(float(np.mean(lam)), 1e-14))
            while occ_l(lam, t, l, hi) < m[l] and hi < 1e12:
                hi = 2 * hi
            for _ in range(100):
                mid = 0.5 * (lo + hi)
                if occ_l(lam, t, l, mid) < m[l]:
                    lo = mid
                else:
                    hi = mid
            t[l] = 0.5 * (lo + hi)
        if np.max(np.abs(t - told) / np.maximum(told, 1e-14)) < 1e-8:
            break
    return t


def _hlru_levelprobs(lam: np.ndarray, t: np.ndarray, h: int) -> np.ndarray:
    """Birth-death level probabilities: pi_l ~ prod_{s<=l} (1-e_s)/e_s."""
    n = len(lam)
    P = np.zeros((n, h + 1))
    for k in range(n):
        w = np.zeros(h + 1)
        w[0] = 1.0
        for l in range(h):
            e = np.exp(-lam[k] * t[l])
            w[1 + l] = w[l] * (1 - e) / max(e, 1e-300)
        P[k, :] = w / np.sum(w)
    return P


def cache_ttl_hlru(gamma: np.ndarray, m: np.ndarray) -> np.ndarray:
    """
    Steady-state list occupancy probabilities for an h-LRU / LRU(m) cache.

    Characteristic-time (TTL) approximation of the list-based h-LRU policy:
    h LRU lists of capacities m[0..h-1], a miss inserts at the head of list 1,
    a hit in list l exchanges the item with the tail of list l+1 (Gast and
    Van Houdt, SIGMETRICS 2015). For h=1 this reduces exactly to the Che
    approximation for LRU.

    Args:
        gamma: Per-item request rates. Accepts (n,), (n x h), or the MVA
               analyzer layout (u x n x h+1) which is aggregated over users.
        m: Cache capacity vector (h,)

    Returns:
        (n x h+1) probabilities; column 0 = not cached, column 1+l = in list l

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_hlru.m
    """
    gamma = np.asarray(gamma, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()
    h = len(m)

    if gamma.ndim == 3:
        # (u x n x h+1) analyzer layout: identical rate across list slices;
        # aggregate the per-item rate over the user classes
        lam = np.sum(gamma[:, :, min(1, gamma.shape[2] - 1)], axis=0)
    elif gamma.ndim == 2:
        lam = gamma[:, 0]
    else:
        lam = gamma

    t = cache_t_hlru(lam, m)
    return _hlru_levelprobs(lam, np.asarray(t, dtype=np.float64), h)


def cache_ttl_lrua(lambd: np.ndarray, R: list, m: np.ndarray,
                   seed: int = 23000,
                   ttl: Optional[np.ndarray] = None) -> np.ndarray:
    """
    Compute steady-state probabilities for TTL-LRU cache with arrival routing.

    Uses fixed-point iteration with DTMC solving for cache systems with
    multiple users, items, and levels with routing.

    Args:
        lambd: Arrival rates per user per item per list (u x n x h+1)
        R: Routing probability structure. Can be either:
           - 1D list: R[i] is the (h+1 x h+1) routing matrix for item i
           - 2D list: R[v][i] is the routing matrix for user v, item i
        m: Cache capacity vector (h,)
        seed: Random seed for initialization (default: 23000)
        ttl: Optional real TTL per cache level (h,). Caps the characteristic
             time at each level: items expire after ttl[l] time units even
             if the cache is not full. Units are in model time (request
             epochs when arrival rates sum to 1). Use np.inf for levels
             with no TTL. When TTL is binding, effective occupancy may be
             less than m[l].

    Returns:
        Steady-state probability distribution (n x h+1)

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_lrua.m
    """
    np.random.seed(seed)

    lambd = np.asarray(lambd, dtype=np.float64)
    m = np.asarray(m, dtype=np.float64).ravel()

    if ttl is not None:
        ttl = np.asarray(ttl, dtype=np.float64).ravel()

    u = lambd.shape[0]  # number of users
    n = lambd.shape[1]  # number of items
    h = lambd.shape[2] - 1  # number of lists

    # Determine if R is 1D (R[i]) or 2D (R[v][i]) structure.
    # R may arrive as a list-of-lists (native default access cost), a list of
    # per-item matrices, or a fully-materialized numpy array (e.g. when the
    # access cost was deserialized from JSON). Detect the layout by effective
    # dimensionality so a per-user-per-item numpy array (nusers,nitems,h+1,h+1)
    # is handled identically to the list-of-lists form. isinstance(R[0], list)
    # alone misclassifies the numpy case and indexes R with the wrong rank.
    try:
        R_arr = np.asarray(R, dtype=float)
    except (ValueError, TypeError):
        R_arr = None
    if R_arr is not None and R_arr.ndim == 4:
        R = R_arr
        is_2d_structure = True
    elif R_arr is not None and R_arr.ndim == 3:
        R = R_arr
        is_2d_structure = False
    else:
        is_2d_structure = isinstance(R[0], list)

    def ttl_tree_time(x):
        """Compute capacity difference for given characteristic times."""
        # Cap characteristic times at real TTL if provided
        x_eff = x.copy()
        if ttl is not None:
            x_eff = np.minimum(x_eff, ttl)

        steadystateprob = np.zeros((n, h + 1))
        randprob = np.zeros((n, h + 1))
        avgtime = np.zeros((n, h + 1))
        capa = np.zeros(h)
        rpdenominator = np.zeros(n)

        for i in range(n):
            # Build transition matrix
            # MATLAB uses R{1,i} which accesses user 1 (first user), item i
            # In Python: R[0][i] for 2D structure, R[i] for 1D structure
            Ri = R[0][i] if is_2d_structure else R[i]
            transmatrix = np.zeros((h + 1, h + 1))

            for j in range(h + 1):
                leafnodes = np.where(Ri[j, :] > 0)[0]
                for k in leafnodes:
                    if j == 0:
                        transmatrix[j, k] = Ri[j, k]
                    else:
                        transmatrix[j, k] = (1 - np.exp(-lambd[0, i, j] * x_eff[j-1])) * Ri[j, k]
                    if j != k and k > 0:
                        transmatrix[k, j] = np.exp(-lambd[0, i, k] * x_eff[k-1])

            # Remove disconnected nodes: keep only states with an INCOMING edge,
            # mirroring MATLAB cache_ttl_lrua.m (missconnection = all-zero columns).
            # Keeping a source-only (outgoing) state retains a transient miss-state
            # and produces an absorbing sub-DTMC that dtmc_solve now rejects.
            connected = []
            for j in range(h + 1):
                if np.any(transmatrix[:, j] > 0):
                    connected.append(j)

            if len(connected) == 0:
                continue

            # Extract submatrix for connected nodes
            submatrix = transmatrix[np.ix_(connected, connected)]

            # Solve DTMC
            dtmcprob = _dtmc_solve(submatrix)

            for idx, node in enumerate(connected):
                steadystateprob[i, node] = dtmcprob[idx]
                if node > 0:
                    rate = lambd[0, i, node]
                    if rate > 0:
                        avgtime[i, node] = (1 - np.exp(-rate * x_eff[node-1])) / rate
                    else:
                        avgtime[i, node] = 0
                else:
                    rate = lambd[0, i, node]
                    avgtime[i, node] = 1 / rate if rate > 0 else 1.0
                rpdenominator[i] += steadystateprob[i, node] * avgtime[i, node]

            for idx, node in enumerate(connected):
                if rpdenominator[i] > 0:
                    randprob[i, node] = (steadystateprob[i, node] *
                                         avgtime[i, node] / rpdenominator[i])

        # Compute capacity difference
        F = np.zeros(h)
        for l in range(h):
            capa[l] = np.sum(randprob[:, l + 1])
            F[l] = m[l] - capa[l]

        return F, randprob

    # Initial guess with random seed
    x0 = np.random.uniform(0, 10, h)

    # Solve using fsolve (standard) or least_squares (when TTL caps apply)
    def objective(x):
        F, _ = ttl_tree_time(x)
        return F

    if ttl is not None:
        # With TTL caps, F=0 may not be achievable — use least_squares
        # with lower bound 0 (characteristic times are positive)
        result = least_squares(objective, x0, bounds=(0, np.inf))
        t = result.x
    else:
        t, info, ier, mesg = fsolve(objective, x0, full_output=True)

        if ier != 1:
            for scale in [0.1, 0.5, 2.0, 5.0]:
                t, info, ier, mesg = fsolve(objective, x0 * scale,
                                            full_output=True)
                if ier == 1:
                    break

    # Get final probabilities
    _, prob = ttl_tree_time(t)

    return prob


def _dtmc_solve(P: np.ndarray) -> np.ndarray:
    """
    Solve for stationary distribution of a DTMC.

    Args:
        P: Transition probability matrix

    Returns:
        Stationary distribution vector
    """
    n = P.shape[0]

    if n == 0:
        return np.array([])

    if n == 1:
        return np.array([1.0])

    # Normalize rows to make it a proper transition matrix
    row_sums = np.sum(P, axis=1, keepdims=True)
    row_sums[row_sums == 0] = 1  # Avoid division by zero
    P = P / row_sums

    # Solve (P^T - I) * pi = 0 with sum(pi) = 1
    # Use the method: replace last equation with sum = 1
    A = P.T - np.eye(n)
    A[-1, :] = 1.0
    b = np.zeros(n)
    b[-1] = 1.0

    try:
        pi = np.linalg.solve(A, b)
        pi = np.maximum(pi, 0)  # Ensure non-negative
        pi = pi / np.sum(pi)  # Normalize
    except np.linalg.LinAlgError:
        # Fallback: power method
        pi = np.ones(n) / n
        for _ in range(1000):
            pi_new = P.T @ pi
            pi_new = pi_new / np.sum(pi_new)
            if np.max(np.abs(pi_new - pi)) < 1e-10:
                break
            pi = pi_new
        pi = pi_new

    return pi




def cache_lrum_map_levelstats(D0: np.ndarray, D1: np.ndarray, T: np.ndarray):
    """
    Level statistics of one item's embedded (list, phase) chain under the
    LRU(m)-MAP TTL approximation (Gast and Van Houdt, PEVA 2017, eqs. 5-9).

    Args:
        D0: (d,d) hidden-transition matrix of the item's MAP
        D1: (d,d) arrival matrix of the item's MAP
        T: characteristic times (h,)

    Returns:
        (prob, occ, hitfrac): time-stationary level probabilities (h+1,),
        list occupancies (h,), request-weighted hit fractions (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_lrum_map_levelstats.m
    """
    from scipy.linalg import expm
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    T = np.asarray(T, dtype=float).ravel()
    h = len(T)
    d = D0.shape[0]
    iD0 = np.linalg.inv(-D0)

    E = [expm(D0 * T[l]) for l in range(h)]
    Nh = [(np.eye(d) - E[l]) @ iD0 for l in range(h)]
    A = [Nh[l] @ D1 for l in range(h)]
    A0 = iD0 @ D1
    N0 = iD0

    # R recursion, eqs. (6)-(7); R[l] holds the paper's R_{l+1}
    R = [None] * h
    for l in range(h - 1, -1, -1):
        if l == h - 1:
            Aprev = A0 if h == 1 else A[l - 1]
            R[l] = Aprev @ np.linalg.inv(np.eye(d) - A[l])
        elif l == 0:
            R[l] = A0 @ np.linalg.inv(np.eye(d) - R[l + 1] @ E[l + 1])
        else:
            R[l] = A[l - 1] @ np.linalg.inv(np.eye(d) - R[l + 1] @ E[l + 1])

    # pi_0: left Perron vector of R_1 e^{D0 T_1} (level-0 balance)
    M = R[0] @ E[0]
    w, vl = np.linalg.eig(M.T)
    idx = int(np.argmax(np.real(w)))
    pi0 = np.real(vl[:, idx])
    pi0 = pi0 / pi0.sum()

    pih = [pi0 @ R[0]]
    for l in range(1, h):
        pih.append(pih[l - 1] @ R[l])

    e = np.ones(d)
    holding = np.zeros(h + 1)
    holding[0] = pi0 @ N0 @ e
    for l in range(h):
        holding[1 + l] = pih[l] @ Nh[l] @ e
    denom = holding.sum()
    prob = holding / denom
    occ = prob[1:]

    # request-weighted hit fractions
    Q = D0 + D1
    wq, vq = np.linalg.eig(Q.T)
    iq = int(np.argmin(np.abs(wq)))
    piphase = np.real(vq[:, iq])
    piphase = piphase / piphase.sum()
    lam = piphase @ D1 @ e
    hitfrac = np.zeros(h)
    for l in range(h):
        hitfrac[l] = (pih[l] @ Nh[l] @ D1 @ e) / denom
    hitfrac = hitfrac / lam if lam > 0 else np.zeros(h)
    return prob, occ, hitfrac


def cache_t_lrum_map(D0c: list, D1c: list, m: np.ndarray) -> np.ndarray:
    """
    Characteristic times for the LRU(m)-MAP TTL approximation.

    Args:
        D0c: list of per-item (d,d) hidden-transition matrices
        D1c: list of per-item (d,d) arrival matrices
        m: cache capacity vector (h,)

    Returns:
        Characteristic times (h,)

    References:
        Original MATLAB: matlab/src/api/cache/cache_t_lrum_map.m
    """
    m = np.asarray(m, dtype=float).ravel()
    n = len(D0c)
    h = len(m)

    def capres(y):
        T = np.exp(y)
        occ = np.zeros(h)
        for k in range(n):
            _, occk, _ = cache_lrum_map_levelstats(D0c[k], D1c[k], T)
            occ += occk
        return occ - m

    y = fsolve(capres, np.zeros(h))
    return np.exp(y)


def cache_ttl_lrum_map(D0c: list, D1c: list, m: np.ndarray):
    """
    Request-weighted hit/miss probabilities for LRU(m) with per-item MAP
    request streams (Gast and Van Houdt, PEVA 2017). Intended for items with
    genuinely distinct or correlated request processes (e.g. marked MAP
    arrivals); when items are i.i.d. marks of a common stream the request
    sequence is IRM and the Poisson-based TTL approximations already apply.

    Args:
        D0c: list of per-item (d,d) hidden-transition matrices
        D1c: list of per-item (d,d) arrival matrices
        m: cache capacity vector (h,)

    Returns:
        (pij, pijtime): (n,h+1) request-weighted probabilities (column 0 =
        miss) and (n,h+1) time-stationary level occupancy probabilities

    References:
        Original MATLAB: matlab/src/api/cache/cache_ttl_lrum_map.m
    """
    m = np.asarray(m, dtype=float).ravel()
    n = len(D0c)
    h = len(m)
    t = cache_t_lrum_map(D0c, D1c, m)
    pij = np.zeros((n, h + 1))
    pijtime = np.zeros((n, h + 1))
    for k in range(n):
        prob, _, hitfrac = cache_lrum_map_levelstats(D0c[k], D1c[k], t)
        pijtime[k, :] = prob
        pij[k, 1:] = hitfrac
        pij[k, 0] = max(0.0, 1.0 - hitfrac.sum())
    return pij, pijtime

__all__ = [
    'cache_lrum_map_levelstats',
    'cache_t_lrum_map',
    'cache_ttl_lrum_map',
    'cache_t_hlru',
    'cache_ttl_hlru',
    'cache_ttl_lrua',
]
