"""
Censored GE/GE/c/K;N queue by entropy maximisation.

Maximum Entropy solution of the single-class censored FCFS GE/GE/c/K;N queue
of Kouvatsos (1994) "Entropy Maximisation and Queueing Network Models",
Section 4.1, equations (4.1)-(4.3). The queue holds at most N jobs and never
fewer than K; arrivals finding N jobs present are turned away and departures
are not allowed from state K. For a queue embedded in an open network K is
always 0; a positive K arises in closed networks, where it records the
minimum occupancy forced by the remaining stations being full.

The ME state probabilities subject to normalisation, the marginal
utilisations, the mean queue length excluding J jobs and the full-buffer
probability coincide with the global balance solution
p(n) = p(K) G_n x^h(n) y^f(n). The Lagrangian coefficients are invariant to
N and K, so letting K -> 0 and N -> inf recovers the stable GE/GE/c solution
used by me_oqn.
"""

from typing import Tuple

import numpy as np


def me_gegecn_pb(p: np.ndarray, K: int, N: int, c: int, Cs: float, Ca: float) -> float:
    """Blocking probability seen by one arrival stream, eq. (4.3).

    Because a GE arrival process is a batch process, an arrival can be
    blocked while the queue holds fewer than N jobs: the factor
    (1-tau)^(N-n) is the probability that the batch overflows the residual
    room, and the extra factor of the first sum accounts for the servers
    still idle. With a Poisson stream (Ca = 1) every term but n = N vanishes
    and PB reduces to the PASTA value p(N).

    The stream scv is a per-stream quantity, so the same node solution yields
    a different blocking probability for the external arrivals, for the flow
    from each upstream station and for the flow released by each holding
    node, which is how PBe, PB^i_j and PB^h_j are obtained in the
    transfer-blocking algorithm of Tahilramani, Manjunath and Bose (1999).

    Args:
        p: ME queue length distribution, p[idx] = Pr{n = K+idx}
        K: minimum number of jobs in the queue
        N: buffer capacity in jobs
        c: number of servers
        Cs: squared coefficient of variation of the service times
        Ca: squared coefficient of variation of the interarrival times OF THE
            STREAM whose blocking probability is requested

    Returns:
        Probability that an arrival of this stream finds the queue full.
    """
    tau = 2.0 / (Ca + 1.0)
    sigma = 2.0 / (Cs + 1.0)
    n = np.arange(K, N + 1)
    w = np.power(1.0 - tau, N - n)
    pb = 0.0
    if K < c:
        last = min(c - K, len(n))
        idx = np.arange(last)
        fac = np.power(sigma / (sigma * (1 - tau) + tau), c - n[idx])
        pb += float(np.sum(w[idx] * fac * p[idx]))
    lo = max(c, K)
    if lo <= N:
        idx2 = np.arange(lo - K, len(n))
        pb += float(np.sum(w[idx2] * p[idx2]))
    return pb


def me_gegecn(lam: float, Ca: float, mu: float, Cs: float, c: int, K: int, N: int
              ) -> Tuple[np.ndarray, float, float, float, float]:
    """Solves a censored GE/GE/c/K;N queue by entropy maximisation.

    Args:
        lam: arrival rate offered to the queue, the arrivals that are turned
            away included
        Ca: squared coefficient of variation of the interarrival times
            (Ca >= 1: the GE distribution is undefined below 1)
        mu: service rate of one server
        Cs: squared coefficient of variation of the service times (Cs >= 1)
        c: number of servers, finite and at least one
        K: minimum number of jobs in the queue
        N: buffer capacity in jobs, in service included (N > K)

    Returns:
        (p, L, U, PB, Lq) with p[idx] = Pr{n = K+idx}, L the mean number of
        jobs, U the utilization E[min(n,c)]/c, PB the blocking probability of
        the queue's own aggregate stream and Lq the mean number waiting.
    """
    if not np.isfinite(N):
        raise ValueError('me_gegecn requires a finite buffer capacity N; use me_oqn for infinite capacity.')
    if not np.isfinite(c) or c < 1:
        raise ValueError('me_gegecn requires a finite number of servers c >= 1.')
    if N <= K:
        raise ValueError('me_gegecn requires N > K.')
    if Ca < 1 - 1e-12 or Cs < 1 - 1e-12:
        raise ValueError('me_gegecn requires Ca >= 1 and Cs >= 1: '
                         'the GE distribution is not defined for scv < 1.')
    if mu <= 0:
        raise ValueError('me_gegecn requires a positive service rate.')

    c = int(round(c))
    K = int(round(K))
    N = int(round(N))

    tau = 2.0 / (Ca + 1.0)
    sigma = 2.0 / (Cs + 1.0)
    rho = lam / (c * mu)

    J = max(c, K + 1)
    den1 = sigma * (1 - tau) + tau
    den2 = tau * rho * (1 - sigma) + sigma

    # Lagrangian coefficients g(l), l = K+1,...,J, stored at index l-1
    g = np.ones(J)
    if K < c - 1:
        g[K] = tau * c * rho / ((K + 1) * den1)
    elif K == c - 1:
        g[K] = tau * sigma * rho / den2
    else:
        g[K] = (den1 / den2) * tau * rho
    for l in range(K + 2, J + 1):
        if l < J:
            g[l - 1] = (tau * c * rho + (l - 1) * sigma * (1 - tau)) / (l * den1)
        else:
            g[l - 1] = sigma * (tau * c * rho + (J - 1) * sigma * (1 - tau)) / (J * den2)

    x = (tau * rho + sigma * (1 - tau)) / den2
    y = 1.0 / (1.0 - (1 - sigma) * x)

    # see _kb/03-api-layer.md for rationale
    n = np.arange(K, N + 1)
    cumlogg = np.cumsum(np.log(g[K:J]))
    logp = np.zeros(len(n))
    for idx, nn in enumerate(n):
        if nn > K:
            m = max(K + 1, min(c, int(nn)))
            logp[idx] = cumlogg[m - K - 1]
    logp = logp + np.maximum(0, n - J) * np.log(x) + np.maximum(0, n - N + 1) * np.log(y)
    logp = logp - np.max(logp)
    p = np.exp(logp)
    p = p / np.sum(p)

    busy = np.minimum(n, c)
    L = float(np.sum(n * p))
    ebusy = float(np.sum(busy * p))
    pb = me_gegecn_pb(p, K, N, c, Cs, Ca)
    return p, L, ebusy / c, pb, L - ebusy
