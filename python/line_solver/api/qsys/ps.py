"""
Exact sojourn-time moments of the multiclass M/M/1-PS queue.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``qsys_mm1_ps.m``.

References:
    D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of the Waiting
    Time in Closed and Open Processor-Sharing Systems with Multiple Job
    Classes", Adv. Appl. Prob. 15(4):813-839, 1983, equation (7).
"""

import numpy as np


def qsys_mm1_ps(lam, mu):
    """Sojourn-time moments of the multiclass M/M/1-PS queue.

    Class j arrives in a Poisson stream of rate ``lam[j]`` and requires an
    exponential amount of service with rate ``mu[j]``. The processor is shared
    equally by all jobs in service, so the class of a job affects its sojourn
    time both through its own service rate and through the mix of rates of the
    jobs it shares the processor with. With
    ``alpha = 1 - sum_j lam[j]/mu[j]`` the unutilized fraction of the
    processor, the moments of the sojourn time ``W_r`` of a tagged class-r job
    are::

        E[W_r]   = 1/(alpha*mu[r])
        E[W_r^2] = 2/(alpha*mu[r])**2
                   * (1 - sum_j lam_j (mu_j-mu_r)/(mu_j(mu_j+mu_r)))
                   / (1 - sum_j lam_j/(mu_j+mu_r))

    which is equation (7) of Mitra and Morrison (1983). Both are exact, not
    asymptotic: the open system is the ``N -> infinity`` limit of the closed
    terminal-driven system whose moments that paper expands in ``1/N``, and the
    leading term of the expansion is exact in the limit. For a single class the
    second moment reduces to the classical ``4/(mu^2 (1-rho)^2 (2-rho))`` of
    Coffman, Muntz and Trotter (1970).

    Parameters
    ----------
    lam : array_like (R,)
        Per-class Poisson arrival rates, non-negative.
    mu : array_like (R,)
        Per-class exponential service rates, positive.

    Returns
    -------
    W : np.ndarray (R,)
        Per-class mean sojourn times.
    W2 : np.ndarray (R,)
        Per-class second moments of the sojourn time.
    alpha : float
        Unutilized fraction of the processor, ``1 - sum_j lam_j/mu_j``.
    """
    lam = np.asarray(lam, dtype=float).flatten()
    mu = np.asarray(mu, dtype=float).flatten()
    R = lam.size
    if mu.size != R:
        raise ValueError("lambda and mu must have the same number of classes")
    if not np.all(np.isfinite(lam)) or np.any(lam < 0):
        raise ValueError("lambda must be finite and non-negative")
    if not np.all(np.isfinite(mu)) or np.any(mu <= 0):
        raise ValueError("mu must be finite and positive")

    alpha = 1.0 - float(np.sum(lam / mu))
    if alpha <= 0:
        raise ValueError("System is unstable: utilization %.6f >= 1" % (1.0 - alpha))

    W = np.zeros(R)
    W2 = np.zeros(R)
    for r in range(R):
        mur = mu[r]
        num = 1.0 - float(np.sum(lam * (mu - mur) / (mu * (mu + mur))))
        den = 1.0 - float(np.sum(lam / (mu + mur)))
        W[r] = 1.0 / (alpha * mur)
        W2[r] = 2.0 / (alpha * mur) ** 2 * num / den
    return W, W2, alpha
