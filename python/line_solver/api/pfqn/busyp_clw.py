"""
Busy period of a subnetwork through point evaluations of the normalizing constant.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_busyp_clw.m``.

WHAT THIS BUYS OVER ``pfqn_busyp`` / ``pfqn_busyp_multiclass``. Those walk the
whole population ladder (the whole lattice, multichain) because the numerator is
a sum over {|m| >= n}. The complement of that set is the SHELLS |m| <= n-1, and
summing the product form over the WHOLE lattice is the full network's own
normalizing constant, since G and H convolve to it:

    sum_{m : |m| >= n} G_I(m) H(N-m) = G(N) - sum_{m : |m| <= n-1} G_I(m) H(N-m)

so the busy period of order n needs only the n lowest shells plus ONE evaluation
of G(N). The ordinary busy period n=1 collapses to three constants:

    b(1,I) = [G(N) - H(N)] / sum_r A_r(I) H(N-e_r)

Those are point evaluations at or near the full population, which is exactly what
the normalizing-constant methods are built for: this routine calls CLW
(Choudhury-Leung-Whitt, J. ACM 42, 1995, numerical inversion of the generating
function) and any other method that returns lG(N) can be dropped in its place.
The cost stops depending on N: it is O(shells up to n-1) plus O(n*R) constant
evaluations, against O(lattice) for the ladder routines.

The OPEN case needs no inversion at all. The subnetwork's constant sequence has
generating function g(z) = prod_{i in I} f_i(z), and the tail the busy period
needs is g(1) minus a partial sum:

    b(n,I) = [g_I(1) - sum_{m=0}^{n-1} G_I(m)] / [G_I(n-1) C(I)]

with f_i(1) = 1/(1-rho_i) at a single server and exp(rho_i) at an infinite one.
That removes the tail TRUNCATION of the ladder routine, not just its cost: the
tail is now exact.

SCOPE. CLW's generating function covers single-server and infinite-server
stations, so a general load-dependent scaling is refused here and belongs to
``pfqn_busyp_multiclass``. The identity above is for the AGGREGATE level set: a
per-class one has complement {m_r <= n-1}, which is the whole lattice in the
other chains and buys nothing, so per-class queries also stay with the ladder.
"""

import numpy as np

from .busyp import _lse
from .busyp_multiclass import _lattice, _lgvec
from .nc import pfqn_clw


def _log_nc(L, Z, N, method):
    """log G(N) of a set of queue stations with an aggregate think time.

    ``L`` is (queues x R) demands, ``Z`` is (R,) the infinite-server aggregate.
    Any method returning lG at a population vector can serve here; CLW is the
    default because its cost is independent of the population.
    """
    if np.all(np.asarray(N) == 0):
        return 0.0
    if L.shape[0] == 0:
        # only infinite servers left: G(N) = prod_r Z_r^N_r / N_r!
        out = 0.0
        for r in range(len(N)):
            if N[r] == 0:
                continue
            if Z[r] <= 0:
                return -np.inf
            out += N[r] * np.log(Z[r]) - float(np.sum(np.log(np.arange(1, N[r] + 1))))
        return out
    if method != 'clw':
        raise ValueError("pfqn_busyp_clw: only the 'clw' method is wired here; the "
                         'point evaluation is a plug-in, so another token needs its '
                         'own call rather than a silent substitution.')
    _, lG = pfqn_clw(np.asarray(L, dtype=float), np.asarray(N, dtype=float),
                     np.asarray(Z, dtype=float))
    return float(lG)


def pfqn_busyp_clw(alpha, mu, P, N, subnet, n, gamma=None, isdelay=None, method='clw'):
    """Mean busy period of order ``n`` for the subnetwork, via NC point evaluations.

    Parameters
    ----------
    alpha : array (J, R)
        Relative arrival rates, one column per chain.
    mu : array (J, R)
        Service rates, the chain-r rate at node j.
    P : array (J, J) or sequence of R such arrays
        Routing, shared by every chain or one matrix per chain.
    N : array (R,)
        Population per chain; ``numpy.inf`` entries for an open chain.
    subnet : sequence of int
        Zero-based node indexes forming the subnetwork.
    n : int or sequence of int
        Busy period order(s), counting the jobs of every chain.
    gamma : array (J, R), optional
        External arrival rates; required for an open network.
    isdelay : array (J,) of bool, optional
        Infinite-server nodes; the rest are single servers.
    method : str
        Method name of the normalizing-constant method used for the point evaluations.

    Returns
    -------
    b : float or np.ndarray
        Mean busy period duration(s).
    """
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    J, R = alpha.shape
    N = np.atleast_1d(np.asarray(N, dtype=float))
    nvec = np.atleast_1d(np.asarray(n, dtype=int))
    is_closed = bool(np.all(np.isfinite(N)))
    is_open = bool(np.all(~np.isfinite(N)))
    if not is_closed and not is_open:
        raise ValueError('pfqn_busyp_clw: a mixed model needs the lattice routine '
                         'pfqn_busyp_multiclass.')
    if isdelay is None:
        isdelay = np.zeros(J, dtype=bool)
    isdelay = np.asarray(isdelay, dtype=bool)

    subnet = np.unique(np.asarray(subnet, dtype=int))
    if subnet.size == 0:
        raise ValueError('The subnetwork must be non-empty.')
    if is_closed and subnet.size >= J:
        raise ValueError(
            'In a closed network the subnetwork must be a proper subset of the nodes.')
    compl = np.setdiff1d(np.arange(J), subnet)

    # demands L(i,r) = alpha(i,r)/mu(i,r), zero where chain r does not visit node i
    L = np.zeros((J, R))
    visited = (alpha > 0) & (mu > 0)
    L[visited] = alpha[visited] / mu[visited]

    # A_r(I): the chain-r rate at which jobs enter the subnetwork from outside it
    A = np.zeros(R)
    for r in range(R):
        Pr = np.asarray(P[r] if isinstance(P, (list, tuple)) else P, dtype=float)
        A[r] = float(alpha[compl, r] @ Pr[np.ix_(compl, subnet)] @ np.ones(subnet.size))
        if gamma is not None:
            A[r] += float(np.atleast_2d(np.asarray(gamma, dtype=float))[subnet, r].sum())
    if A.sum() <= 0:
        raise ValueError('No job ever enters the subnetwork, its busy period is undefined.')

    if is_open:
        # g_I(1) in closed form, so the tail is exact rather than truncated
        rho = L[subnet, :].sum(axis=1)
        if np.any(rho[~isdelay[subnet]] >= 1):
            raise ValueError('The subnetwork is not stable, its busy period is infinite.')
        lg1 = 0.0
        for t, i in enumerate(subnet):
            lg1 += rho[t] if isdelay[i] else -np.log1p(-rho[t])
        # the n lowest coefficients of the same generating function
        kmax = int(nvec.max())
        lseq = _open_coefficients(L[subnet, :], isdelay[subnet], kmax)
        b = np.zeros(nvec.shape)
        for t, nt in enumerate(nvec):
            head = _lse(lseq[:nt]) if nt > 0 else -np.inf
            tail = lg1 + np.log1p(-np.exp(head - lg1)) if head > -np.inf else lg1
            b[t] = np.exp(tail - lseq[nt - 1] - np.log(A.sum()))
        return float(b[0]) if np.asarray(n).ndim == 0 else b

    Nint = N.astype(int)
    kmax = int(nvec.max()) - 1
    # the low shells of the subnetwork, the only lattice this routine touches
    bound = np.minimum(Nint, max(kmax, 0))
    mvec, stride = _lattice(bound)
    phi_sub = _phi_of(isdelay[subnet], max(1, int(bound.sum())))
    lG_low = _lgvec(L[subnet, :], phi_sub, mvec, stride, bound)
    level = mvec.sum(axis=1)

    # point evaluations: the full network at N, the complement near N
    queues = ~isdelay
    lG_full = _log_nc(L[queues, :], L[isdelay, :].sum(axis=0), Nint, method)

    def lH(k):
        k = np.asarray(k, dtype=int)
        if np.any(k < 0):
            return -np.inf
        rows = compl[queues[compl]]
        delays = compl[isdelay[compl]]
        return _log_nc(L[rows, :], L[delays, :].sum(axis=0) if delays.size else np.zeros(R),
                       k, method)

    b = np.zeros(nvec.shape)
    for t, nt in enumerate(nvec):
        # numerator: the full constant minus the shells the level set excludes
        corr = []
        for idx in np.where(level <= nt - 1)[0]:
            corr.append(lG_low[idx] + lH(Nint - mvec[idx]))
        lcorr = _lse(corr)
        num = lG_full + np.log1p(-np.exp(lcorr - lG_full))
        # denominator: the flow out of the shell |m| = n-1
        den = -np.inf
        for idx in np.where(level == nt - 1)[0]:
            terms = []
            for r in range(R):
                if A[r] <= 0:
                    continue
                left = Nint - mvec[idx]
                left[r] -= 1
                if np.any(left < 0):
                    continue
                terms.append(np.log(A[r]) + lH(left))
            if terms:
                den = _lse([den, lG_low[idx] + _lse(terms)])
        b[t] = np.exp(num - den)
    return float(b[0]) if np.asarray(n).ndim == 0 else b


def _phi_of(isdelay, width):
    """The dimensionless scaling the ladder routine expects: 1, or k at a delay."""
    phi = np.ones((len(isdelay), width))
    for i, d in enumerate(isdelay):
        if d:
            phi[i, :] = np.arange(1, width + 1)
    return phi


def _open_coefficients(L, isdelay, kmax):
    """log G_I(0..kmax) of an OPEN subnetwork, by convolving the per-node series.

    The chains are already absorbed into the per-node load, since in an open
    network the total occupancy depends on the aggregate load alone.
    """
    rho = L.sum(axis=1)
    lg = np.full(kmax + 1, -np.inf)
    lg[0] = 0.0
    for i in range(len(rho)):
        li = np.zeros(kmax + 1)
        acc = 0.0
        for k in range(1, kmax + 1):
            # 1/(1-rho z) has coefficients rho^k; exp(rho z) has rho^k/k!
            acc += np.log(rho[i]) - (np.log(k) if isdelay[i] else 0.0)
            li[k] = acc
        lgnew = np.full(kmax + 1, -np.inf)
        for m in range(kmax + 1):
            lgnew[m] = _lse(lg[m::-1] + li[:m + 1])
        lg = lgnew
    return lg
