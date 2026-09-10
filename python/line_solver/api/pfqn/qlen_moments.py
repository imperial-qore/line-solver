"""
Joint queue-length moments of a closed product-form network, from normalizing
constants.

Two exact routes reach the joint survival array of the queue lengths, and both
end in the same conversion, the tail edge of the house of moments (api/moment)
followed by the joint central-moment and cumulant conversions.

  * SINGLE CLASS (R = 1). The survival probabilities of the per-station queue
    lengths are ratios of normalizing constants of the network itself,

        P(n_i >= k_i for all i) = (prod_i L_i^k_i) * G(N - sum_i k_i) / G(N)

    This holds because a load-independent single-class station has the
    geometric occupancy L_i^n, so the survival event factors. Only N+1
    constants of the ORIGINAL model are needed, which is why any normalizing
    constant algorithm serves it unchanged.

  * MULTICLASS (R >= 1, general). The geometric factorization fails, because a
    multiclass load-independent station carries the multinomial occupancy
    f_i(n_i) = (|n_i|)! prod_r L_(i,r)^n_(i,r) / n_(i,r)!. What does hold is the
    joint distribution of the selected stations in terms of the COMPLEMENTARY
    network, the model with those stations deleted and the think times kept,

        P(n_i = m_i, i in S) = prod_(i in S) f_i(m_i) * G_(S^c)(N - sum_i m_i)
                               / G(N)

    The survival array is then the reverse cumulative sum of that probability
    array, exactly, since the box covers the support.

Neither the factorial nor the raw moments have a one-constant closed form; the
survival array is the queue-length functional that does, which is why the tail
edge exists in the moment API. The normalizing-constant algorithm is INJECTED,
never called at a fixed site: the whole set of populations is known before any
evaluation, so it is emitted in one batch, and an algorithm that produces
several constants in one pass (convolution, CoMoM) serves it without
recomputation.

References:
    M. Reiser and S. S. Lavenberg. Mean-value analysis of closed multichain
    queuing networks. Journal of the ACM, 27(2):313-322, 1980.
"""

from math import exp, factorial, log

import numpy as np

from ..moment import (moment_joint_binomial_from_tail,
                      moment_joint_central_from_raw,
                      moment_joint_cumulant_from_raw,
                      moment_joint_factorial_from_binomial,
                      moment_joint_raw_from_factorial)


def _reverse_cumsum(A):
    """
    Joint survival array of a joint probability array, by a reverse cumulative
    sum along every dimension.

    Args:
        A: ndarray of joint probabilities covering the support.

    Returns:
        ndarray of the same shape with element k equal to the probability that
        every coordinate is at least k.
    """
    out = np.asarray(A, dtype=float)
    for axis in range(out.ndim):
        out = np.flip(np.cumsum(np.flip(out, axis=axis), axis=axis), axis=axis)
    return out


def _batch_lg(Lsub, pops, Z, lg_source, method, options):
    """
    Evaluate log G at a batch of populations, honouring the injected source.

    Args:
        Lsub: Demand matrix of the network whose constants are wanted, with at
            least one station.
        pops: (P x R) integer array of populations, already deduplicated.
        Z: Think time vector.
        lg_source: None, a callable (L, pops) -> vector of log G with NaN where
            unavailable, or a precomputed ndarray table indexed by population.
        method: Method for pfqn_nc on the populations left unserved.
        options: Options for pfqn_nc.

    Returns:
        Tuple (lg, served, evals): the values, how many the source supplied and
        how many pfqn_nc calls were needed.
    """
    from .nc import pfqn_nc

    P = pops.shape[0]
    lg = np.full(P, np.nan)
    served = 0
    if isinstance(lg_source, np.ndarray):
        for p in range(P):
            lg[p] = lg_source[tuple(pops[p, :])]
        served = int(np.sum(np.isfinite(lg)))
    elif callable(lg_source):
        got = np.asarray(lg_source(Lsub, pops), dtype=float).ravel()
        if got.size != P:
            raise ValueError('pfqn_qlen_joint_moments: the lg_source callable must return '
                             'one value per requested population.')
        lg = got
        served = int(np.sum(np.isfinite(lg)))
    elif lg_source is not None:
        raise ValueError('pfqn_qlen_joint_moments: lg_source must be None, a callable or '
                         'an array.')

    evals = 0
    for p in range(P):
        if np.isfinite(lg[p]):
            continue
        res = pfqn_nc(Lsub, pops[p, :].astype(float), Z, method=method, options=options)
        # native pfqn_nc returns (G, log G) in that order, unlike the MATLAB
        # pfqn_nc whose first output is lG
        lg[p] = res[1] if isinstance(res, tuple) else res
        evals += 1
    return lg, served, evals


def _delay_lg(Z, n):
    """
    Log normalizing constant of a pure-delay network, prod_r Z_r^n_r / n_r!.

    Args:
        Z: Think time vector.
        n: Population vector.

    Returns:
        The logarithm, or -inf when a class has a positive population and a
        zero think time, which makes that state unreachable.
    """
    acc = 0.0
    for r in range(len(n)):
        if n[r] == 0:
            continue
        if Z[r] <= 0:
            return -np.inf
        acc += n[r] * log(Z[r]) - log(float(factorial(int(n[r]))))
    return acc


def pfqn_qlen_joint_moments(L, N, Z=None, pairs=None, route='auto', lg_source=None,
                            method='ca', options=None):
    """
    Joint moments of the queue-length vector of a closed product-form network.

    The coordinates are (station, class) pairs. Two pairs sharing a class give
    the cross-station covariance of that class; two pairs sharing a station give
    the cross-class covariance at that station, which is what a class-oriented
    method of moments (pfqn_comomrm and its relatives) is positioned to deliver.

    The result is EXACT: the survival array covers the support, since a queue
    length is bounded by the population of its class.

    Args:
        L: Service demand matrix of the QUEUEING stations (M x R). Delay
           (infinite-server) stations belong in Z, their marginals following a
           different law. Load-dependent or multiserver stations are out of
           scope for both routes and must not be passed here.
        N: Population vector (R,).
        Z: Think time vector (R,), zeros if omitted.
        pairs: Sequence of (station, class) 0-based pairs, one per dimension of
               the returned arrays. Defaults to every class of every station.
        route: 'tail' uses the single-class survival identity and only touches
               the original network; 'pmf' uses the complementary-network joint
               distribution and works for any number of classes; 'auto' picks
               'tail' when R = 1 and 'pmf' otherwise.
        lg_source: Where log G comes from. None calls pfqn_nc. A callable is
                   invoked ONCE per network with (Lsub, pops), pops being a
                   (P x R) integer array, and must return P values of log G with
                   NaN where it cannot serve; those are filled in by pfqn_nc. An
                   ndarray is read as a precomputed table indexed by population,
                   which is what a convolution sweep produces for free. Note
                   that the 'pmf' route queries the COMPLEMENTARY network, so a
                   table must be its table, not the original network's.
        method: Method passed to pfqn_nc for unserved populations.
        options: Options passed to pfqn_nc.

    Returns:
        Dictionary with the joint arrays over the selected coordinates: 'tail'
        (survival), 'binomial', 'factorial', 'raw', 'central', 'cumulant', plus
        'mean', the covariance matrix 'cov', and 'info' holding the route, the
        number of populations requested, how many the source served, how many
        pfqn_nc evaluations were needed, and the pairs used.

    Raises:
        ValueError: If the arguments are inconsistent, if a pair is out of
            range, or if the 'tail' route is requested with several classes.

    Example:
        out = pfqn_qlen_joint_moments(L, N, Z, pairs=[(0, 0), (0, 1)])
        cov01 = out['cov'][0, 1]
    """
    L = np.atleast_2d(np.asarray(L, dtype=float))
    N = np.asarray(N, dtype=int).ravel()
    M, R = L.shape
    if N.size != R:
        raise ValueError('pfqn_qlen_joint_moments: N must have one entry per class.')
    Z = np.zeros(R) if Z is None else np.asarray(Z, dtype=float).ravel()
    if Z.size != R:
        raise ValueError('pfqn_qlen_joint_moments: Z must have one entry per class.')

    if pairs is None:
        pairs = [(i, r) for i in range(M) for r in range(R)]
    pairs = [(int(i), int(r)) for i, r in pairs]
    if not pairs:
        raise ValueError('pfqn_qlen_joint_moments: at least one (station,class) pair is '
                         'required.')
    for i, r in pairs:
        if i < 0 or i >= M or r < 0 or r >= R:
            raise ValueError('pfqn_qlen_joint_moments: the pair (%d,%d) is out of range.'
                             % (i, r))

    if route == 'auto':
        route = 'tail' if R == 1 else 'pmf'
    if route not in ('tail', 'pmf'):
        raise ValueError("pfqn_qlen_joint_moments: route must be 'auto', 'tail' or 'pmf'.")
    if route == 'tail' and R > 1:
        raise ValueError('pfqn_qlen_joint_moments: the tail route needs the geometric '
                         'occupancy of a single-class load-independent station; with '
                         'several classes the multinomial factor breaks the survival '
                         "identity, so use route='pmf'.")

    if len(set(pairs)) != len(pairs):
        raise ValueError('pfqn_qlen_joint_moments: the (station,class) pairs must be '
                         'distinct.')
    d = len(pairs)
    dims = tuple(int(N[r]) + 1 for _, r in pairs)

    if route == 'tail':
        # the whole population set is known up front: N minus the total order
        need = []
        for a in np.ndindex(*dims):
            n = N - sum(a)
            if np.all(n >= 0):
                need.append(n.copy())
        need.append(N.copy())
        need = np.unique(np.vstack(need), axis=0)
        lg, served, evals = _batch_lg(L, need, Z, lg_source, method, options)
        index = {tuple(need[p, :]): lg[p] for p in range(need.shape[0])}
        lgN = index[tuple(N)]
        tail = np.zeros(dims)
        for a in np.ndindex(*dims):
            n = N - sum(a)
            if np.any(n < 0):
                continue
            acc = 0.0
            ok = True
            for j, (i, r) in enumerate(pairs):
                if a[j] == 0:
                    continue
                if L[i, r] <= 0:
                    ok = False
                    break
                acc += a[j] * log(L[i, r])
            if not ok:
                continue
            tail[a] = exp(acc + index[tuple(n)] - lgN)
    else:
        # see _kb/03-api-layer.md for rationale
        stations = sorted(set(i for i, _ in pairs))
        coords = [(i, r) for i in stations for r in range(R)]
        cdims = tuple(int(N[r]) + 1 for _, r in coords)
        Lsub = np.delete(L, stations, axis=0)
        need = []
        for a in np.ndindex(*cdims):
            n = N.copy()
            for j, (_, r) in enumerate(coords):
                n[r] -= a[j]
            if np.all(n >= 0):
                need.append(n.copy())
        need = np.unique(np.vstack(need), axis=0)
        if Lsub.shape[0] == 0:
            lgc = np.array([_delay_lg(Z, need[p, :]) for p in range(need.shape[0])])
            served, evals = need.shape[0], 0
        else:
            lgc, served, evals = _batch_lg(Lsub, need, Z, lg_source, method, options)
        indexc = {tuple(need[p, :]): lgc[p] for p in range(need.shape[0])}
        lgN, _, evals0 = _batch_lg(L, N.reshape(1, -1), Z, None, method, options)
        evals += evals0
        lgN = lgN[0]

        pmf = np.zeros(cdims)
        for a in np.ndindex(*cdims):
            n = N.copy()
            for j, (_, r) in enumerate(coords):
                n[r] -= a[j]
            if np.any(n < 0):
                continue
            gc = indexc[tuple(n)]
            if not np.isfinite(gc):
                continue
            acc = gc - lgN
            ok = True
            for i in stations:
                tot = sum(a[j] for j, (ii, _) in enumerate(coords) if ii == i)
                acc += log(float(factorial(int(tot))))
                for j, (ii, r) in enumerate(coords):
                    if ii != i or a[j] == 0:
                        continue
                    if L[ii, r] <= 0:
                        ok = False
                        break
                    acc += a[j] * log(L[ii, r]) - log(float(factorial(int(a[j]))))
                if not ok:
                    break
            if ok:
                pmf[a] = exp(acc)
        # marginalize onto the requested pairs, then take the survival array
        keep = [coords.index(pr) for pr in pairs]
        drop = tuple(j for j in range(len(coords)) if j not in keep)
        marg = pmf.sum(axis=drop) if drop else pmf
        order = np.argsort(np.argsort(keep))
        marg = np.transpose(marg, axes=order)
        tail = _reverse_cumsum(marg)
        dims = marg.shape

    b = moment_joint_binomial_from_tail(tail)
    f = moment_joint_factorial_from_binomial(b)
    m = moment_joint_raw_from_factorial(f)
    mc = moment_joint_central_from_raw(m)
    kap = moment_joint_cumulant_from_raw(m)

    mean = np.zeros(d)
    cov = np.zeros((d, d))
    for j in range(d):
        e = [0] * d
        e[j] = 1
        mean[j] = m[tuple(e)]
        for l in range(d):
            a = [0] * d
            a[j] += 1
            a[l] += 1
            cov[j, l] = kap[tuple(a)]

    return {'tail': tail, 'binomial': b, 'factorial': f, 'raw': m, 'central': mc,
            'cumulant': kap, 'mean': mean, 'cov': cov,
            'info': {'route': route, 'points': int(need.shape[0]), 'served': served,
                     'evals': evals, 'exact': True, 'pairs': pairs, 'dims': tuple(dims)}}
