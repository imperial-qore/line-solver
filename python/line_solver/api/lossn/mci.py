"""
Loss Network Analysis via Monte Carlo Importance-Sampling Summation.

Native Python implementation of the Monte Carlo summation method of
Ross and Wang, "Monte Carlo Summation Applied to Product-Form Loss
Networks", Probability in the Engineering and Informational Sciences,
6 (1992), 323-348.

Estimates the product-form normalization constant g(C) and class
blocking probabilities of a loss network by importance sampling,
providing confidence intervals. Unlike the Erlang fixed point, the
estimator makes no link-independence assumption and is consistent.
"""

import numpy as np
from scipy import special
from typing import Optional, Tuple, Dict


def _logsumexp(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=np.float64).ravel()
    m = np.max(x)
    if np.isinf(m):
        return float(m)
    return float(m + np.log(np.sum(np.exp(x - m))))


def lossn_mci(nu: np.ndarray, A: np.ndarray, C: np.ndarray,
              samples: int = 100000, gamma: Optional[np.ndarray] = None,
              seed: Optional[int] = None, alpha: float = 0.05
              ) -> Tuple[np.ndarray, np.ndarray, float, Dict, int]:
    """
    Monte Carlo importance-sampling summation for loss networks.

    A loss network has links j=1..J with capacity C[j] and classes r=1..R
    with offered load nu[r] and per-link circuit requirement A[j,r]. The
    state n is feasible iff A @ n <= C (set Omega). The product-form
    normalization constant is g(C) = sum_{n in Omega} prod_r nu[r]**n_r/n_r!.
    Class-r acceptance is g(C-A[:,r])/g(C) = 1 - beta_r.

    States are drawn from the importance distribution (Eq. 6)
        p(n) = (1/c) prod_r gamma_r**n_r / n_r!
    over the box {0..N_1} x ... x {0..N_R}, N_r = min_j floor(C_j/A_jr).
    Ratio estimators (Eq. 8) yield g and blocking with delta-method
    confidence intervals.

    Args:
        nu: Offered load per class (R,).
        A: Circuit requirement matrix (J, R).
        C: Link capacity vector (J,).
        samples: Number of Monte Carlo samples.
        gamma: Importance-sampling parameters (R,); default is the
            Section 3.4 heuristic.
        seed: RNG seed for reproducibility.
        alpha: Confidence-interval significance level (default 0.05).

    Returns:
        Tuple (qlen, loss, lG, ci, nsamples) where:
            - qlen: Mean carried load E[n_r] per class (R,).
            - loss: Blocking probability beta_r per class (R,).
            - lG: Log of the estimated normalization constant g(C).
            - ci: dict with 'accept' (R,2), 'loss' (R,2),
              'acceptPoint' (R,), 'lossPoint' (R,), 'level'.
            - nsamples: Number of samples used.
    """
    nu = np.asarray(nu, dtype=np.float64).ravel()
    A = np.asarray(A, dtype=np.float64)
    C = np.asarray(C, dtype=np.float64).ravel()
    R = len(nu)
    J = len(C)
    if A.shape != (J, R):
        raise ValueError(f"A must have shape ({J}, {R}), got {A.shape}")

    rs = np.random.RandomState(seed)
    S = int(samples)

    # Per-class maximum feasible occupancy N_r = min_j floor(C_j/A_jr)
    N = np.zeros(R, dtype=np.int64)
    for k in range(R):
        pos = A[:, k] > 0
        if np.any(pos):
            N[k] = int(np.floor(np.min(C[pos] / A[pos, k])))
        else:
            N[k] = 0

    # Importance-sampling parameters gamma (Section 3.4 heuristic)
    if gamma is None:
        load_j = np.array([np.sum(A[j, :] * nu) / C[j] for j in range(J)])
        delta = float(np.max(load_j))
        b = np.max(A, axis=0)
        base = max(1.0 - 0.15 * (1.0 - delta), 1e-6)
        gamma = nu * (base ** b)
    gamma = np.asarray(gamma, dtype=np.float64).ravel()
    gamma = np.maximum(gamma, 1e-300)

    # Normalization constant c of the importance distribution (log space)
    log_c = 0.0
    for k in range(R):
        l = np.arange(N[k] + 1)
        log_c += _logsumexp(l * np.log(gamma[k]) - special.gammaln(l + 1))

    # Draw S i.i.d. samples V (S, R), each column truncated Poisson(gamma_k)
    V = np.zeros((S, R))
    for k in range(R):
        l = np.arange(N[k] + 1)
        logpmf = l * np.log(gamma[k]) - special.gammaln(l + 1)
        pmf = np.exp(logpmf - _logsumexp(logpmf))
        cdf = np.cumsum(pmf)
        cdf[-1] = 1.0
        u = rs.random_sample(S)
        V[:, k] = np.sum(u[:, None] > cdf[None, :], axis=1)

    AV = V @ A.T                              # (S, J)
    inOmega = np.all(AV <= C[None, :], axis=1)

    logratio = np.log(nu) - np.log(gamma)
    log_alpha = V @ logratio                  # (S,)

    # Normalization constant estimate g(C) = (c/S) sum alpha 1(Omega)
    laO = log_alpha[inOmega]
    if laO.size == 0:
        lG = -np.inf
    else:
        m = laO.max()
        lG = log_c + m + np.log(np.sum(np.exp(laO - m))) - np.log(S)

    # Ratio estimators with shifted weights for numerical stability
    M = log_alpha[inOmega].max() if np.any(inOmega) else 0.0
    w = np.exp(log_alpha - M)
    Z = w * inOmega
    meanZ = Z.mean()

    crit = np.sqrt(2.0) * special.erfinv(1.0 - alpha)
    accept = np.zeros(R)
    acceptCI = np.zeros((R, 2))
    for k in range(R):
        Ck = (C - A[:, k])[None, :]
        inK = np.all(AV <= Ck, axis=1)
        Y = w * inK
        meanY = Y.mean()
        if meanZ <= 0:
            phi = np.nan
            half = np.nan
        else:
            phi = meanY / meanZ
            varY = Y.var(ddof=1)
            varZ = Z.var(ddof=1)
            covYZ = np.sum((Y - meanY) * (Z - meanZ)) / (S - 1)
            sig2 = max((varY - 2 * phi * covYZ + phi ** 2 * varZ) / (S * meanZ ** 2), 0.0)
            half = crit * np.sqrt(sig2)
        accept[k] = phi
        acceptCI[k, :] = [phi - half, phi + half]

    loss = 1.0 - accept
    qlen = nu * accept
    ci = {
        'accept': acceptCI,
        'loss': np.column_stack([1.0 - acceptCI[:, 1], 1.0 - acceptCI[:, 0]]),
        'acceptPoint': accept,
        'lossPoint': loss,
        'level': 1.0 - alpha,
    }
    return qlen, loss, float(lG), ci, S


__all__ = ['lossn_mci']
