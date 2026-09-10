"""Chen-O'Cinneide regularization: a Markov chain Monte Carlo estimator of the class
throughputs and queue lengths of a closed multiclass product-form network.

W. Chen, C. A. O'Cinneide, "Towards a Polynomial-Time Randomized Algorithm for Closed
Product-Form Networks", ACM TOMACS 8(3):227-253, 1998.
"""

from typing import NamedTuple

import numpy as np

FINE_TOL = 1e-8

__all__ = ['pfqn_mcmc', 'PfqnMcmcResult']


class PfqnMcmcResult(NamedTuple):
    """Estimates of :func:`pfqn_mcmc` with their batch-means confidence intervals.

    There is deliberately no normalizing constant here: the estimator is a ratio of
    holding-time weighted averages that yields G(N-e_r)/G(N) directly, and G itself
    never enters the algorithm.
    """

    X: np.ndarray       #: (R,) throughput estimates G(N-e_r)/G(N)
    Q: np.ndarray       #: (M, R) mean queue lengths at the queueing stations
    Xse: np.ndarray     #: (R,) batch-means standard error of X
    Xlo: np.ndarray     #: (R,) lower end of the two-sigma interval for X
    Xhi: np.ndarray     #: (R,) upper end of the two-sigma interval for X
    Qse: np.ndarray     #: (M, R) batch-means standard error of Q
    Qlo: np.ndarray     #: (M, R) lower end of the two-sigma interval for Q
    Qhi: np.ndarray     #: (M, R) upper end of the two-sigma interval for Q
    batches: int        #: batches the run was split into
    samples: int        #: service completions simulated after warm-up
    burnin: int         #: completions discarded as warm-up


def pfqn_mcmc(L, N, Z=None, s=None, options=None) -> PfqnMcmcResult:
    """Markov chain Monte Carlo estimate of the throughputs and queue lengths.

    Estimates the class throughputs ``X(r) = G(N-e_r)/G(N)`` and the mean queue lengths
    ``Q(i,r)`` of a CLOSED multiclass product-form (BCMP, no type changes) network by the
    REGULARIZATION algorithm of Chen and O'Cinneide (ACM TOMACS 8(3), 1998).

    The three steps of the paper are:

    I. CONSTRUCT THE REGULARIZED NETWORK. Write ``rho(i,r)`` for the surrogate traffic
    intensity of class r at station i -- here the service demand, since ``rho = lambda/mu``
    is a visit ratio over a service rate -- and ``rho(r) = sum_i rho(i,r)``. The
    regularized network has the same stations, classes and populations, UNIT service rates
    at every station, the processor-sharing discipline, and a routing matrix that depends
    on the destination only, ``P*(i->m | class r) = rho(m,r)/rho(r)``. By Theorem 2.1 it is
    a REVERSIBLE chain with the SAME steady-state distribution as the original network, and
    its throughputs satisfy ``Theta*(r) = rho(r)*Theta(r)``.

    II. SIMULATE IT at service-completion epochs. With ``Y(i,r)`` the number of class-r
    jobs at station i, ``Y(i)`` their total and ``Psi_i(k) = min(s_i,k)`` the number of
    busy servers::

        r(i,r) = Y(i,r)/Y(i) * Psi_i(Y(i)),   r(r) = sum_i r(i,r),
        r      = sum_i Psi_i(Y(i)),

    the next completion is of class r at station i with probability ``r(i,r)/r``, and the
    conditional expected time to it is ``1/r``. Equation (10) of the paper is the
    holding-time weighted ratio estimator ``Theta*(r) = sum_t r(r,t)/r(t) / sum_t 1/r(t)``,
    and the same weights give the time-average queue lengths, which need no transformation
    at all because the two networks share their steady state.

    III. TRANSFORM BACK: ``X(r) = Theta*(r)/rho(r)``.

    Because P* forgets the station of origin and every station serves at unit rate, the
    regularized chain has neither the slowly mixing routing chain nor the
    customer-trapping slow station that make the original chain converge slowly. The paper
    proves ``O(N^2*M^3)`` mixing in two special cases (Section 4) and reports the general
    behaviour experimentally (Section 5).

    Delay (infinite-server) demand enters as ONE extra station with ``s = inf`` and demand
    Z. Aggregating infinite-server stations that way is exact in the product form, since
    their joint term is multinomial in the per-class totals.

    Confidence: the run is split into non-overlapping batches (Schmeiser 1982, 30 by
    default, the count used in the tables of the paper), the batch means of the ratio
    estimator give a standard error, and the intervals are the paper's two-sigma ones. The
    estimator is a ratio of correlated averages, so it carries an ``O(1/samples)`` bias on
    top of the initialization bias; the paper ignores both, this implementation
    additionally discards a warm-up fraction (10% by default).

    Parameters
    ----------
    L : (M, R) array
        Per-class service demands at the M queueing stations.
    N : (R,) array
        Closed population vector; finite and integer.
    Z : (R,) or (K, R) array, optional
        Aggregated think times, summed over rows; None or zeros if the model has no delay.
    s : (M,) array, optional
        Number of servers at each queueing station, ``inf`` for an infinite server; None
        means all stations single-server.
    options : dict or options object, optional
        Fields ``samples`` (default 1e5), ``seed``, and inside ``config`` the batch count
        ``mcmc_batches`` (30) and the warm-up fraction ``mcmc_burnin`` (0.1).

    Returns
    -------
    PfqnMcmcResult
        Throughputs, queue lengths and their two-sigma intervals.

    Examples
    --------
    >>> mu = np.array([0.2, 0.5, 0.8]); sets = [[0, 1, 2], [0, 1], [0, 2], [1, 2]]
    >>> L = np.zeros((3, 4)); Z = np.zeros(4)
    >>> for c, st in enumerate(sets):
    ...     L[st, c] = 1.0 / mu[st]
    ...     if c > 0:
    ...         Z[c] = 1.0 / 0.5
    >>> res = pfqn_mcmc(L, 3 * np.ones(4), Z, options={'samples': 100000, 'seed': 23000})

    See Also
    --------
    pfqn_nc, pfqn_mci, pfqn_ls, pfqn_is
    """
    from .pas import _opt

    L = np.atleast_2d(np.asarray(L, dtype=float))
    L[~np.isfinite(L)] = 0.0
    M, R = L.shape

    N = np.asarray(N, dtype=float).ravel()
    if np.any(~np.isfinite(N)):
        raise ValueError('pfqn_mcmc requires a closed model, but the population vector '
                         'has an infinite entry.')
    if np.any(np.abs(N - np.round(N)) > FINE_TOL):
        # the chain lives on the integer lattice sum_i Y(i,r) = N(r), so a fractional
        # population has no state space at all; this is not a matter of accuracy and must
        # not be rounded away silently
        raise ValueError('pfqn_mcmc simulates a state space of integer populations, but '
                         'N = %s is fractional. Use an asymptotic method (\'le\', '
                         '\'ble\', \'kt\') on fractional populations.' % np.array2string(N))
    N = np.round(N).astype(int)

    if Z is None:
        Zr = np.zeros(R)
    else:
        Zr = np.atleast_2d(np.asarray(Z, dtype=float))
        Zr = np.where(np.isfinite(Zr), Zr, 0.0)
        Zr = Zr.sum(axis=0).ravel()
        if Zr.size < R:
            Zr = np.concatenate([Zr, np.zeros(R - Zr.size)])
        Zr = Zr[:R]

    X = np.zeros(R)
    Q = np.zeros((M, R))
    if N.sum() == 0:
        return PfqnMcmcResult(X, Q, np.zeros(R), np.zeros(R), np.zeros(R),
                              np.zeros((M, R)), np.zeros((M, R)), np.zeros((M, R)),
                              0, 0, 0)

    if s is None:
        svec0 = np.ones(M)
    else:
        svec0 = np.asarray(s, dtype=float).ravel()
        if svec0.size != M:
            raise ValueError('pfqn_mcmc: the server count vector has %d entries but L has '
                             '%d stations.' % (svec0.size, M))

    # ---- Step I: the regularized network -------------------------------------------
    # Only the surrogate traffic intensities rho(i,r) enter the product form, and scaling a
    # whole class column by a constant leaves the steady-state distribution unchanged, so
    # the demands are used as they are.
    rho = np.maximum(L, 0.0)
    svec = svec0
    if np.any(Zr > 0):
        rho = np.vstack([rho, Zr[np.newaxis, :]])
        svec = np.concatenate([svec, [np.inf]])
    Mx = rho.shape[0]
    rhoTot = rho.sum(axis=0)
    bad = np.nonzero((N > 0) & (rhoTot <= 0))[0]
    if bad.size:
        raise ValueError('pfqn_mcmc: class %d has a positive population but no demand '
                         'anywhere in the network.' % (bad[0] + 1))

    # routing of the regularized network, P*(m|r) = rho(m,r)/rho(r), held as one column of
    # cumulative probabilities per class
    Pstar = rho / np.maximum(rhoTot, np.finfo(float).tiny)[np.newaxis, :]
    cumP = np.cumsum(Pstar, axis=0)
    cumP[Mx - 1, :] = 1.0  # guard the last bin against a floating-point shortfall

    # run length, batching and warm-up
    samples = int(max(1, round(float(_opt(options, 'samples', 100000)))))
    config = _opt(options, 'config', None)
    nbatches = 30  # Schmeiser (1982), the count used in the tables of the paper
    burnin_frac = 0.1
    if config is not None:
        nbatches = int(max(1, round(float(_opt(config, 'mcmc_batches', 30)))))
        burnin_frac = min(0.9, max(0.0, float(_opt(config, 'mcmc_burnin', 0.1))))
    batch_len = max(1, samples // nbatches)
    samples = batch_len * nbatches
    nburn = int(round(burnin_frac * samples))
    seed = _opt(options, 'seed', None)
    rng = np.random.default_rng(seed if seed is not None else None)

    # Initial state: spread each class over the stations it can occupy in the proportions
    # P*(.|r), by largest remainder. That is the marginal the regularized network would
    # have with no queueing, so it costs nothing and starts the chain far closer to
    # stationarity than a single-station state.
    Y = np.zeros((Mx, R))
    for r in range(R):
        if N[r] == 0:
            continue
        target = N[r] * Pstar[:, r]
        base = np.floor(target)
        short = int(N[r] - base.sum())
        if short > 0:
            order = np.argsort(-(target - base), kind='stable')
            base[order[:short]] += 1.0
        Y[:, r] = base
    Ytot = Y.sum(axis=1)

    # ---- Step II: simulate at service-completion epochs ------------------------------
    xnum = np.zeros((nbatches, R))       # sum_t r(r,t)/r(t) within the batch
    qnum = np.zeros((nbatches, Mx, R))   # sum_t Y(t)/r(t)   within the batch
    den = np.zeros(nbatches)             # sum_t 1/r(t)      within the batch
    stale = True
    Psi = np.zeros(Mx)
    cPsi = np.zeros(Mx)
    rvec = np.zeros(R)
    w = 0.0
    horizon = nburn + samples
    for t in range(1, horizon + 1):
        if stale:
            # (8)-(9): busy servers, per-class completion rates, total rate
            Psi = np.minimum(svec, Ytot)
            cPsi = np.cumsum(Psi)
            rw = np.where(Ytot > 0, Psi / np.where(Ytot > 0, Ytot, 1.0), 0.0)
            rvec = rw @ Y
            w = 1.0 / cPsi[Mx - 1]
            stale = False
        if t > nburn:
            b = (t - nburn - 1) // batch_len
            den[b] += w
            xnum[b, :] += w * rvec
            qnum[b, :, :] += w * Y
        # Pick the completing station with probability Psi(i)/r, then the completing class
        # within it with probability Y(i,r)/Y(i); the product is the r(i,r)/r of the paper,
        # since sum_r Y(i,r)/Y(i)*Psi(i) = Psi(i).
        i = int(np.searchsorted(cPsi, rng.random() * cPsi[Mx - 1], side='left'))
        if i >= Mx:
            i = int(np.nonzero(Psi > 0)[0][-1])
        cY = np.cumsum(Y[i, :])
        cls = int(np.searchsorted(cY, rng.random() * cY[R - 1], side='left'))
        if cls >= R:
            cls = int(np.nonzero(Y[i, :] > 0)[0][-1])
        # Route it. A self-transition leaves the state, hence the rates and the weight,
        # unchanged: skipping the recomputation is the saving described at the end of
        # Section 2 of the paper.
        m = int(np.searchsorted(cumP[:, cls], rng.random(), side='left'))
        if m < Mx and m != i:
            Y[i, cls] -= 1.0
            Y[m, cls] += 1.0
            Ytot[i] -= 1.0
            Ytot[m] += 1.0
            stale = True

    # ---- Step III: back to the original network --------------------------------------
    # Theta(r) = Theta*(r)/rho(r) by (7) and (11); the queue lengths transfer unchanged,
    # the two networks sharing their steady-state distribution.
    den_tot = den.sum()
    X = (xnum.sum(axis=0) / den_tot) / rhoTot
    Qx = qnum.sum(axis=0) / den_tot
    Q = Qx[:M, :]

    # batch-means standard error and the two-sigma interval of the paper
    if nbatches > 1:
        Xb = (xnum / den[:, np.newaxis]) / rhoTot[np.newaxis, :]
        Qb = qnum / den[:, np.newaxis, np.newaxis]
        Xse = Xb.std(axis=0, ddof=1) / np.sqrt(nbatches)
        Qse = Qb.std(axis=0, ddof=1)[:M, :] / np.sqrt(nbatches)
    else:
        Xse = np.zeros(R)
        Qse = np.zeros((M, R))

    return PfqnMcmcResult(X, Q, Xse, X - 2 * Xse, X + 2 * Xse,
                          Qse, Q - 2 * Qse, Q + 2 * Qse,
                          nbatches, samples, nburn)
