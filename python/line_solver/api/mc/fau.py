"""
CTMC transient analysis via fast adaptive uniformization.

Adaptive uniformization (van Moorsel and Sanders, 1994) draws the
uniformization rate of each step from the states the iterate actually
occupies rather than from the whole state space, so the subordinating process
is a pure birth process instead of a Poisson process. Fast adaptive
uniformization (Mateescu, Wolf, Didier and Henzinger, 2010) adds the dropping
of states below an occupancy threshold, which is what turns a population
cutoff into a numerical one.

Key algorithms:
    ctmc_fau: Transient probabilities by fast adaptive uniformization
"""

import math
from dataclasses import dataclass
from typing import Tuple

import numpy as np

try:
    import scipy.sparse as sp
except ImportError:  # pragma: no cover - scipy is a hard dependency elsewhere
    sp = None

from .foxglynn import ctmc_foxglynn_weights

# Default cap on the number of birth steps, so a pathological horizon fails
# loudly through info.truncated rather than running forever.
FAU_MAX_STEPS = 1000000


@dataclass
class CtmcFauInfo:
    r"""
    Diagnostics of a fast adaptive uniformization sweep.

    Attributes:
        steps: number of birth steps K+1 actually taken
        lambda_min: smallest adaptive rate used
        lambda_max: largest adaptive rate used, the Lstar of the weights
        uniform_rate: max_i \|q_ii\|, the rate ordinary uniformization would use
        weight_tail: mass reaching the overflow index, that is P{N(t) > K}
        weight_window: Poisson mass outside the Fox-Glynn window of the weights
        dropped_mass: probability removed by the occupancy threshold
        error_bound: sum(pi0) - sum(pit), which IS the L1 error
        support_max: largest occupied support over the sweep
        support_final: support at the last step
        truncated: True if maxsteps stopped the sweep
        absorbed: True if the support emptied or became absorbing
    """
    steps: int
    lambda_min: float
    lambda_max: float
    uniform_rate: float
    weight_tail: float
    weight_window: float
    dropped_mass: float
    error_bound: float
    support_max: int
    support_final: int
    truncated: bool
    absorbed: bool


def _rows_combination(Q, u_act: np.ndarray, act: np.ndarray) -> np.ndarray:
    """
    The row combination sum_i u_i Q[i,:] over the active states only, which is
    u @ Q when u vanishes off act. Sparse Q is sliced by row, which is why the
    caller is asked for CSR.
    """
    if sp is not None and sp.issparse(Q):
        return np.asarray(Q[act, :].transpose().dot(u_act)).ravel()
    return u_act @ Q[act, :]


def _fau_step(u: np.ndarray, act: np.ndarray, Q, L: float,
              delta: float) -> Tuple[np.ndarray, np.ndarray, float]:
    """
    One adaptive uniformization step u <- u(I + Q/L), touching only the rows of
    Q in the current support, followed by the drop rule. A state with a zero
    exit rate is absorbing: its row of Q is empty, so it holds its mass and
    stays in the support.
    """
    contrib = _rows_combination(Q, u[act], act)
    idx = np.nonzero(contrib)[0]
    dropped = 0.0
    if idx.size == 0:
        return u, act, dropped
    vals = u[idx] + contrib[idx] / L
    small = vals < delta
    if np.any(small):
        dropped = float(np.sum(np.maximum(vals[small], 0.0)))
        vals[small] = 0.0
    u[idx] = vals
    # A row of Q that cancels exactly leaves its state untouched and out of
    # idx, so the surviving part of the old support is carried over too.
    act = np.union1d(act[u[act] > 0.0], idx[vals > 0.0])
    return u, act, dropped


def _fau_tailbound(lstar: float, t: float, k: int) -> float:
    """
    Upper bound on P{N(t) >= k} for the birth process, through the stochastic
    domination of its k-th jump epoch by an Erlang(k, lstar): the bound is the
    Poisson(lstar*t) upper tail P{X >= k}, taken at its Chernoff exponent
    lam*h(k/lam) with h(u) = u*log(u) - u + 1. That exponent bounds the upper
    tail only above the mean, so below it the bound is left vacuous.
    """
    lam = lstar * t
    if lam <= 0.0 or k <= lam:
        return 1.0
    return math.exp(-(lam - k + k * math.log(k / lam)))


def _fau_rates(pi0: np.ndarray, Q, d: np.ndarray, t: float, delta: float,
               maxsteps: int, epsilon: float) -> Tuple[np.ndarray, bool, bool]:
    """
    Sweep the iterate to collect the adaptive rates Lambda_0..Lambda_K,
    stopping when the Poisson-dominance bound on P{N(t) > K} falls to epsilon.
    """
    lam = []
    truncated = False
    absorbed = False
    u = pi0.copy()
    act = np.nonzero(u > 0.0)[0]
    lstar = 0.0
    while True:
        if act.size == 0:
            absorbed = True
            break
        L = float(np.max(d[act]))
        lam.append(L)
        if L <= 0.0:
            # Every occupied state is absorbing: the birth process stops here
            # and the remaining weight falls entirely on this iterate.
            absorbed = True
            break
        lstar = max(lstar, L)
        if _fau_tailbound(lstar, t, len(lam)) <= epsilon:
            break
        if len(lam) >= maxsteps:
            truncated = True
            break
        u, act, _ = _fau_step(u, act, Q, L, delta)
    return np.array(lam, dtype=np.float64), truncated, absorbed


def _fau_weights(lam: np.ndarray, t: float,
                 tol: float) -> Tuple[np.ndarray, float, float]:
    """
    Transient distribution of the pure birth process with rates lam at time t,
    that is b_n = P{N(t) = n} for n = 0..K, plus the mass that reached the
    absorbing overflow index K+1 and therefore measures P{N(t) > K}.

    The chain is uniformized at lstar = max(lam) and mixed against Fox-Glynn
    Poisson weights, so the kernel entries 1 - lam_n/lstar and lam_n/lstar are
    probabilities and nothing cancels. The weights are taken UNNORMALIZED, so
    the Poisson mass outside the window is missing from b rather than
    redistributed over it: b is then a sub-distribution, every term of the
    mixture is an underestimate, and the error stays measurable as missing
    mass.
    """
    k1 = lam.size
    if k1 == 0:
        return np.zeros(0, dtype=np.float64), 0.0, 0.0
    lstar = float(np.max(lam))
    if lstar <= 0.0 or t <= 0.0:
        b = np.zeros(k1, dtype=np.float64)
        b[0] = 1.0
        return b, 0.0, 0.0
    left, right, w = ctmc_foxglynn_weights(lstar * t, tol, -1, False)
    wwin = max(1.0 - float(np.sum(w)), 0.0)
    a = 1.0 - lam / lstar
    c = lam / lstar
    v = np.zeros(k1 + 1, dtype=np.float64)
    v[0] = 1.0
    b = np.zeros(k1 + 1, dtype=np.float64)
    for k in range(right + 1):
        if k >= left:
            b += w[k - left] * v
        if k < right:
            forward = v[:k1] * c
            v[:k1] *= a
            v[1:k1 + 1] += forward
    return b[:k1], float(b[k1]), wwin


def _fau_accumulate(pi0: np.ndarray, Q, d: np.ndarray, delta: float,
                    lam: np.ndarray,
                    b: np.ndarray) -> Tuple[np.ndarray, float, int, int]:
    """
    Replay the sweep of _fau_rates, accumulating sum_n b_n u^(n). The
    arithmetic is identical, so the rates and the drops reproduce those of the
    first pass.
    """
    nsteps = lam.size
    pit = np.zeros(pi0.size, dtype=np.float64)
    dropped = 0.0
    u = pi0.copy()
    act = np.nonzero(u > 0.0)[0]
    support_max = act.size
    support_final = act.size
    for m in range(nsteps):
        if act.size == 0:
            break
        support_max = max(support_max, act.size)
        support_final = act.size
        pit[act] += b[m] * u[act]
        if m < nsteps - 1:
            L = float(np.max(d[act]))
            if L <= 0.0:
                break
            u, act, drop_step = _fau_step(u, act, Q, L, delta)
            dropped += drop_step
    return pit, dropped, support_max, support_final


def ctmc_fau(pi0: np.ndarray, Q, t: float, epsilon: float = 1e-6,
             delta: float = 1e-12,
             maxsteps: int = -1) -> Tuple[np.ndarray, CtmcFauInfo]:
    r"""
    Transient distribution of a CTMC at time t by fast adaptive uniformization.

    Ordinary uniformization fixes one rate q >= max_i \|q_ii\| over the whole
    state space and mixes the powers of P = I + Q/q against a Poisson(q*t)
    law, so its cost is set by the fastest state anywhere, including states
    that carry no probability at time t. Adaptive uniformization instead picks
    a rate per step from the states the iterate occupies,

        Lambda_n >= max{\|q_ii\| : i in supp(u^(n))},
        u^(n+1) = u^(n)(I + Q/Lambda_n),

    which keeps every entry of u^(n+1) nonnegative. The subordinating process
    is then the pure birth process N(t) with rates Lambda_0, Lambda_1, ... and

        pi(t) = sum_{n>=0} P{N(t) = n} u^(n).

    The fast variant drops an entry of u^(n) below delta rather than
    propagating it, so the support tracks the states of non-negligible
    occupancy instead of the reachable set. Nothing is renormalized anywhere,
    so the error is not estimated but measured: the birth index truncated at
    K, the Poisson window of the weight computation and the delta threshold
    each remove mass and none puts any back, whence

        0 <= pi(t) - pit componentwise, and
        \|pi(t) - pit\|_1 = sum(pi0) - sum(pit) = info.error_bound.

    The birth weights are computed exactly rather than quadratured, by
    uniformizing the bidiagonal birth generator; see _fau_weights. The sweep
    runs twice because b_n(t) needs the rates up to n, which are not known
    before the sweep ends, while u^(n) is needed after them, and storing every
    iterate would cost K times the support. Stopping is certified by
    stochastic domination of the birth epochs by an Erlang, so this method
    never takes more steps than uniformization at the largest rate it visited.

    This is a transient method: it produces no stationary distribution.

    Args:
        pi0: Initial probability distribution
        Q: Infinitesimal generator matrix, dense or scipy sparse (CSR is used
            as given, other sparse formats are converted)
        t: Time horizon, t >= 0
        epsilon: Birth-process truncation tolerance
        delta: Occupancy threshold below which a state is dropped
        maxsteps: Cap on birth steps; nonpositive for the default cap

    Returns:
        Tuple of (defective distribution at time t, diagnostics)
    """
    if epsilon <= 0.0:
        epsilon = 1e-6
    if delta < 0.0:
        delta = 0.0
    if maxsteps <= 0:
        maxsteps = FAU_MAX_STEPS

    pi0 = np.asarray(pi0, dtype=np.float64).flatten()
    if sp is not None and sp.issparse(Q):
        Q = Q.tocsr()
        n = Q.shape[0]
        d = -np.asarray(Q.diagonal(), dtype=np.float64).ravel()
    else:
        Q = np.asarray(Q, dtype=np.float64)
        n = Q.shape[0]
        d = -np.diag(Q).astype(np.float64)
    if Q.shape[1] != n:
        raise ValueError("Q must be square.")
    if pi0.size != n:
        raise ValueError("pi0 and Q have inconsistent sizes.")
    if t < 0.0:
        raise ValueError("t must be nonnegative.")

    uniform_rate = float(np.max(d)) if n > 0 else 0.0
    if t == 0.0 or n == 0:
        support = int(np.count_nonzero(pi0))
        return pi0.copy(), CtmcFauInfo(1, 0.0, 0.0, uniform_rate, 0.0, 0.0,
                                       0.0, 0.0, support, support, False,
                                       False)

    lam, truncated, absorbed = _fau_rates(pi0, Q, d, t, delta, maxsteps,
                                          epsilon)
    b, wtail, wwin = _fau_weights(lam, t, epsilon)
    pit, dropped, support_max, support_final = _fau_accumulate(
        pi0, Q, d, delta, lam, b)

    info = CtmcFauInfo(
        steps=int(lam.size),
        lambda_min=float(np.min(lam)) if lam.size > 0 else 0.0,
        lambda_max=float(np.max(lam)) if lam.size > 0 else 0.0,
        uniform_rate=uniform_rate,
        weight_tail=wtail,
        weight_window=wwin,
        dropped_mass=dropped,
        error_bound=float(np.sum(pi0) - np.sum(pit)),
        support_max=support_max,
        support_final=support_final,
        truncated=truncated,
        absorbed=absorbed,
    )
    return pit, info
