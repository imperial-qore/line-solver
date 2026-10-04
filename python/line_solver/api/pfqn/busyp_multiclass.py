"""
Mean busy period of order n for a subnetwork of a multichain product-form network.

Native Python implementation (no JPype / JVM dependency). Mirrors the MATLAB
reference ``pfqn_busyp_multiclass.m``.

Multichain generalization of :func:`pfqn_busyp`. Daduna (J. ACM 35(3), 1988)
states Theorems 1 and 3 for a single chain and notes in Section 5 that they
carry over to the whole product-form class: the proof uses only that the
stationary law is product form and that the busy period is Keilson's mean
ergodic sojourn time on a level set, neither of which is single-chain. Replacing
the scalar population by a per-chain vector m gives, for a closed network,

               sum_{m : |m| >= n}   G_I(m) H(N-m)
    b(n,I) = --------------------------------------------------
               sum_{m : |m| = n-1}  G_I(m) sum_r A_r(I) H(N-m-e_r)

with G_I and H the normalizing constants of the subnetwork and of its complement
at a population VECTOR and A_r(I) the chain-r arrival flow into I. At R=1 the
inner sum holds the single term m=n-1 and H(N-m-e_1)=H(N-n), so the expression
collapses to Theorem 1 exactly.

The OPEN case needs no lattice: in an open product-form network the stations are
independent and the total occupancy of a node depends on the AGGREGATE load
sum_r alpha_ir/mu_ir alone, since summing the station function over the
compositions of t collapses the multinomial to (sum_r rho_ir)^t. It is therefore
reduced here to the single-chain routine on aggregated demands.
"""

from math import lgamma

import numpy as np

from .busyp import pfqn_busyp, _lse


def _factln(k):
    """log(k!)."""
    return lgamma(k + 1.0)


def _lattice(N):
    """Every population vector 0 <= m <= N, and the strides of its linear index."""
    R = len(N)
    stride = np.ones(R, dtype=np.int64)
    for r in range(1, R):
        stride[r] = stride[r - 1] * (N[r - 1] + 1)
    size = int(np.prod(np.asarray(N, dtype=np.int64) + 1))
    mvec = np.zeros((size, R), dtype=np.int64)
    for idx in range(size):
        for r in range(R):
            mvec[idx, r] = (idx // stride[r]) % (N[r] + 1)
    return mvec, stride


def _station(Li, phii, mvec):
    """log X_i(m) over the lattice for one node.

    X_i(m) = multinomial(|m|; m) prod_r L(i,r)^m_r / prod_{k=1}^{|m|} phi_i(k),
    which at R=1 is the prod_k alpha_i/mu_i(k) of the single-chain routine and at
    phi(k)=k the infinite-server form prod_r L^m_r/m_r!.
    """
    size, R = mvec.shape
    cols = len(phii)
    out = np.zeros(size)
    for idx in range(size):
        m = mvec[idx]
        tot = int(m.sum())
        v = _factln(tot)
        ok = True
        for r in range(R):
            if m[r] > 0:
                if Li[r] <= 0:
                    ok = False
                    break
                v += -_factln(int(m[r])) + m[r] * np.log(Li[r])
        if not ok:
            out[idx] = -np.inf
            continue
        for k in range(1, tot + 1):
            v -= np.log(phii[min(k, cols) - 1])
        out[idx] = v
    return out


def _lgvec(L, phi, mvec, stride, N):
    """Log normalizing constants over the whole lattice of a set of nodes.

    A node whose scaling row is all ones takes the Buzen recursion, O(R) per
    lattice point; any other node needs the full sub-lattice convolution.
    """
    nodes, R = L.shape
    size = mvec.shape[0]
    lg = np.full(size, -np.inf)
    lg[0] = 0.0
    Nmax = max(1, int(np.sum(N)))
    for i in range(nodes):
        cols = min(phi.shape[1], Nmax)
        if np.all(phi[i, :cols] == 1):
            lgnew = lg.copy()
            for idx in range(size):
                acc = lgnew[idx]
                for r in range(R):
                    if mvec[idx, r] > 0 and L[i, r] > 0:
                        acc = _lse([acc, np.log(L[i, r]) + lgnew[idx - stride[r]]])
                lgnew[idx] = acc
            lg = lgnew
        else:
            lX = _station(L[i], phi[i], mvec)
            lgnew = np.full(size, -np.inf)
            for a in range(size):
                if lg[a] == -np.inf:
                    continue
                for c in range(size):
                    if lX[c] == -np.inf:
                        continue
                    s = mvec[a] + mvec[c]
                    if np.all(s <= N):
                        j = int(1 + (s * stride).sum()) - 1
                        lgnew[j] = _lse([lgnew[j], lg[a] + lX[c]])
            lg = lgnew
    return lg


def pfqn_busyp_multiclass(alpha, mu, P, N, subnet, n, gamma=None, phi=None, tol=1e-12,
                          jobclass=-1):
    """Mean busy period of order ``n`` for the subnetwork, multichain.

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
    phi : array (J, K), optional
        Dimensionless load-dependent scaling; ``None`` means a single server.
    tol : float
        Relative tolerance of the open-network tail truncation.
    jobclass : int
        Zero-based chain whose own jobs are counted, or -1 to count every chain.
        A per-class order is bounded by that chain's population, not by the total.

    Returns
    -------
    b, lG, lH
        Mean duration(s), and the log normalizing constants of the subnetwork
        and of its complement over the lattice (``lH`` empty when open).
    """
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    mu = np.atleast_2d(np.asarray(mu, dtype=float))
    J, R = alpha.shape
    N = np.atleast_1d(np.asarray(N, dtype=float))
    if N.size != R:
        raise ValueError('The population vector must have one entry per chain.')
    nvec = np.atleast_1d(np.asarray(n, dtype=int))
    is_closed = bool(np.all(np.isfinite(N)))
    is_open = bool(np.all(~np.isfinite(N)))
    is_mixed = not is_closed and not is_open

    subnet = np.unique(np.asarray(subnet, dtype=int))
    if subnet.size == 0:
        raise ValueError('The subnetwork must be non-empty.')
    if is_closed and subnet.size >= J:
        raise ValueError(
            'In a closed network the subnetwork must be a proper subset of the nodes.')
    if subnet.min() < 0 or subnet.max() >= J:
        raise ValueError('The subnetwork indexes are out of range.')
    compl = np.setdiff1d(np.arange(J), subnet)

    if phi is None:
        width = max(1, int(np.sum(N[np.isfinite(N)])))
        phi = np.ones((J, width))
    phi = np.atleast_2d(np.asarray(phi, dtype=float))

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
    jobclass = int(jobclass)
    if jobclass >= R:
        raise ValueError('The job class index is out of range.')
    if jobclass >= 0 and A[jobclass] <= 0:
        raise ValueError(
            'No job of that class ever enters the subnetwork, its busy period is undefined.')

    if is_open and not is_mixed:
        # exact reduction to the single-chain routine on a per-station scalar; the
        # synthetic problem carries no routing, the whole inflow riding on gamma
        rho = L.sum(axis=1)
        gsyn = np.zeros(J)
        if jobclass < 0:
            # the total occupancy depends on the AGGREGATE load alone
            scalar = rho
            gsyn[subnet[0]] = A.sum()
        else:
            # the class-r marginal is geometric in sigma_ir = rho_ir/(1-rho_i+rho_ir),
            # NOT in rho_ir: the other classes inflate the queue the class-r jobs sit
            # in. That collapse assumes a load-INDEPENDENT station, since phi does not
            # factor per class.
            if np.any(phi[subnet, :] != 1):
                raise ValueError(
                    'A per-class busy period of an open subnetwork requires load-'
                    'independent stations: under a load-dependent scaling the class '
                    'marginal is no longer geometric and the station needs the pair '
                    '(n_ir, |n_i|) tracked before convolving.')
            denom = 1.0 - rho + L[:, jobclass]
            scalar = np.where(denom > 0, L[:, jobclass] / denom, 0.0)
            gsyn[subnet[0]] = A[jobclass]
        b, lG, _ = pfqn_busyp(scalar, phi, np.zeros((J, J)), np.inf, subnet, n, gsyn, tol)
        return b, lG, np.zeros(0)

    closed_chains = np.where(np.isfinite(N))[0]
    open_chains = np.where(~np.isfinite(N))[0]
    if jobclass < 0:
        bound = int(N[closed_chains].sum()) if not is_mixed else 0
    else:
        bound = int(N[jobclass]) if np.isfinite(N[jobclass]) else 0
    if nvec.min() < 1 or (bound > 0 and nvec.max() > bound):
        raise ValueError(
            'The busy period order must be an integer in 1..sum(N), or in 1..N(r) for '
            'the busy period of class r alone.')

    if not is_mixed:
        Nint = N.astype(int)
        return _lattice_busyp(L, phi, subnet, compl, Nint, nvec, A, jobclass,
                              open_chains, n)

    # A MIXED MODEL keeps the same lattice with its OPEN dimensions TRUNCATED. The
    # closed chains are conserved between the subnetwork and its complement, the open
    # ones are not: the complement's open count is free, so its open dimensions are
    # summed out and no e_r shift applies to an open chain (removing one job from an
    # unbounded dimension leaves the same sum). The truncation is grown until the
    # answer stops moving, which is the only approximation in the mixed branch.
    trunc = 8 + 2 * int(nvec.max())
    prev = None
    while True:
        Nint = np.array([int(N[r]) if np.isfinite(N[r]) else trunc for r in range(R)])
        out = _lattice_busyp(L, phi, subnet, compl, Nint, nvec, A, jobclass,
                             open_chains, n)
        cur = np.atleast_1d(np.asarray(out[0], dtype=float))
        if prev is not None and np.all(np.abs(cur - prev) <= 1e-10 * np.abs(cur)):
            return out
        prev = cur
        trunc *= 2
        if trunc > 4096:
            raise ValueError(
                'The mixed busy period did not converge: the truncated open dimension '
                'keeps growing, so some station of the subnetwork is nearly saturated.')


def _lattice_busyp(L, phi, subnet, compl, Nint, nvec, A, jobclass, open_chains, n):
    """The lattice evaluation shared by the closed and the mixed branch.

    ``Nint`` bounds every chain: the population for a closed one, the truncation for
    an open one. ``open_chains`` names the dimensions that are NOT conserved, whose
    complement counts are therefore summed out rather than read at ``N-m``.
    """
    R = L.shape[1]
    mvec, stride = _lattice(Nint)
    lG = _lgvec(L[subnet, :], phi[subnet, :], mvec, stride, Nint)
    lH = _lgvec(L[compl, :], phi[compl, :], mvec, stride, Nint)

    # Hbar(k) sums the complement over its unconserved dimensions, so it is indexed
    # by the CLOSED components alone; with no open chain it is lH itself.
    if open_chains.size == 0:
        lHbar = lH
    else:
        lHbar = np.full(lH.shape, -np.inf)
        keep = np.ones(R, dtype=bool)
        keep[open_chains] = False
        base = (mvec * keep * stride).sum(axis=1)
        for idx in range(lH.size):
            j = int(base[idx])
            lHbar[j] = _lse([lHbar[j], lH[idx]])

    # the level set is |m| for the aggregate busy period and m_r for the class-r one;
    # only chain-r arrivals move m_r, so the flow sum then holds the single term r
    level = mvec.sum(axis=1) if jobclass < 0 else mvec[:, jobclass]
    chains = range(R) if jobclass < 0 else [jobclass]
    is_open_chain = np.zeros(R, dtype=bool)
    is_open_chain[open_chains] = True
    keep = ~is_open_chain

    b = np.zeros(nvec.shape)
    for t, nt in enumerate(nvec):
        # numerator: the stationary weight of the level set
        sel = np.where(level >= nt)[0]
        rest = (Nint - mvec[sel]) * keep
        num = _lse(lG[sel] + lHbar[(rest * stride).sum(axis=1)])
        # denominator: the flow into the subnetwork out of the level n-1 shell
        sel = np.where(level == nt - 1)[0]
        den = -np.inf
        for idx in sel:
            terms = []
            for r in chains:
                if A[r] <= 0:
                    continue
                left = (Nint - mvec[idx]) * keep
                if not is_open_chain[r]:
                    # a closed chain conserves jobs, so the departing one is removed
                    left[r] -= 1
                    if left[r] < 0:
                        continue
                terms.append(np.log(A[r]) + lHbar[int((left * stride).sum())])
            if terms:
                den = _lse([den, lG[idx] + _lse(terms)])
        b[t] = np.exp(num - den)
    if np.isscalar(n) or np.asarray(n).ndim == 0:
        return float(b[0]), lG, lH
    return b, lG, lH
