"""
MAP flow-equivalent server aggregation.

Native Python implementation of the load-dependent MAP flow-equivalent server of
Casale, Mi, Cherkasova and Smirni, "Dealing with burstiness in multi-tier
applications: models and their parameterization", IEEE Trans. Soft. Eng. 37(5),
2011, Section 5.2. Unlike the classic flow-equivalent server, which keeps only the
mean throughput of the aggregated subnetwork, this one also carries the burstiness
of its departure stream, so a bottleneck switch across the aggregated resources
stays visible to the rest of the model.

References:
    Original MATLAB: matlab/src/api/fes/fes_map_*.m
    JAR: jar/src/main/java/jline/api/fes/Fes_map_*.java
"""

import numpy as np
from scipy.sparse import csc_matrix, csr_matrix, eye as speye, kron as spkron, lil_matrix
from scipy.sparse.linalg import spsolve
from typing import List, Optional, Sequence, Tuple

from ..mam.map_analysis import map_idc, map_moment, map2_fit_idc

MAP = Tuple[np.ndarray, np.ndarray]

TOL = 1e-12
STEP_SAFETY = 0.1
ITER_MAX = 1000000
NHEAD = 10
NTAIL = 10


def fes_map_levels(fes, n: int, mi: float = 1.0) -> List[MAP]:
    """
    Expand a MAP into the per-level processes of a load-dependent server.

    A single MAP is replicated over the levels and scaled by min(k,mi), which
    reproduces a queue with mi servers and, for mi infinite, a delay station
    serving at rate k*mu. The scaling is exact for exponential service and is the
    load-dependent rate approximation otherwise. A list of MAPs is validated and
    returned unchanged.

    Args:
        fes: either a MAP (D0,D1) or a list of MAPs, one per level
        n: number of levels required
        mi: number of servers, np.inf for a delay

    Returns:
        List of n MAPs, entry k-1 being the process active when k jobs are held
    """
    if isinstance(fes, (tuple, list)) and len(fes) == 2 and isinstance(fes[0], np.ndarray) \
            and fes[0].ndim == 2 and isinstance(fes[1], np.ndarray) and fes[1].ndim == 2 \
            and fes[0].shape == fes[1].shape:
        return [(min(k, mi) * fes[0], min(k, mi) * fes[1]) for k in range(1, n + 1)]

    if len(fes) < n:
        raise ValueError("The flow-equivalent server is defined for %d levels but %d are required."
                         % (len(fes), n))
    mf = fes[0][0].shape[0]
    for k in range(n):
        if fes[k][0].shape[0] != mf:
            raise ValueError("All levels of a flow-equivalent server must have the same number of phases.")
    return [fes[k] for k in range(n)]


def fes_map_interdeparture(maps, fes, n: int, mi: Sequence[float] = (1.0, 1.0)):
    """
    Build the MAP (T0,T1) of the inter-departure times of a closed subnetwork made
    of one MAP station and one MAP flow-equivalent server.

    Level k is the population of the flow-equivalent server, so the station holds
    n-k jobs and both processes may be load dependent. Marked transitions are the
    completions of the station, which are the departures fed to the rest of the
    model.

    Args:
        maps: service process of the station, a MAP or one MAP per level
        fes: flow-equivalent server, a MAP or one MAP per level
        n: number of jobs circulating in the subnetwork
        mi: [servers of the station, servers of the flow-equivalent server]

    Returns:
        Tuple (T0,T1) of sparse matrices
    """
    if n < 1:
        raise ValueError("The subnetwork population n must be at least 1.")

    maps_lev = fes_map_levels(maps, n, mi[0])
    fes_lev = fes_map_levels(fes, n, mi[1])
    ms = maps_lev[0][0].shape[0]
    mf = fes_lev[0][0].shape[0]

    blk = ms * mf
    dim = (n + 1) * blk
    ims = speye(ms, format='csr')
    imf = speye(mf, format='csr')

    T0 = lil_matrix((dim, dim))
    T1 = lil_matrix((dim, dim))

    for k in range(n + 1):
        off = k * blk
        j = n - k
        diag = None
        if j > 0:
            diag = spkron(csr_matrix(maps_lev[j - 1][0]), imf)
        if k > 0:
            fes_diag = spkron(ims, csr_matrix(fes_lev[k - 1][0]))
            diag = fes_diag if diag is None else diag + fes_diag
        T0[off:off + blk, off:off + blk] = diag

        if k > 0:
            T0[off:off + blk, off - blk:off] = spkron(ims, csr_matrix(fes_lev[k - 1][1]))
        if j > 0:
            T1[off:off + blk, off + blk:off + 2 * blk] = spkron(csr_matrix(maps_lev[j - 1][1]), imf)

    return T0.tocsc(), T1.tocsc()


def fes_map_euler(v: np.ndarray, T0, dt: float, tol: float = TOL, iter_max: int = ITER_MAX) -> np.ndarray:
    """
    Approximate v*(-T0)^-1 by the trapezoid rule with the Euler propagator.

    Evaluates v*int_0^inf exp(T0 t) dt without any factorization, as done in
    Section 5.2.2 of the reference. The propagated vector uses exp(T0 dt) ~ I + T0
    dt, so only sparse vector-matrix products are performed.

    Args:
        v: row vector to be multiplied by (-T0)^-1
        T0: hidden transitions of the MAP, a stable matrix
        dt: integration step, below 1/max(abs(diag(T0)))
        tol: relative mass left when the integration stops
        iter_max: maximum number of integration steps

    Returns:
        Row vector approximating v*(-T0)^-1
    """
    y = np.zeros_like(v)
    z = v.copy()
    nrm0 = np.sum(np.abs(v))
    for _ in range(iter_max):
        znext = z + dt * np.asarray(T0.T.dot(z)).ravel()
        y = y + dt * (z + znext) / 2
        z = znext
        if np.sum(np.abs(z)) <= tol * nrm0:
            break
    return y


def fes_map_moments(T0, T1, method: str = 'ssolve', step_safety: float = STEP_SAFETY,
                    tol: float = TOL, iter_max: int = ITER_MAX):
    """
    First three moments, lag-1 joint moment and index of dispersion of a MAP.

    Evaluates equations (4), (5) and (7) of the reference. The inverse (-T0)^-1 is
    dense even when T0 is sparse, so it is never formed: the moments follow from
    the vector recursion v_{k+1} = v_k (-T0)^-1, each step being a sparse linear
    solve. Method 'euler' replaces the solve by the quadrature of the paper.

    Args:
        T0: hidden transitions of the MAP
        T1: marked transitions of the MAP
        method: 'ssolve' for the sparse solve, 'euler' for the quadrature
        step_safety: fraction of the uniformization bound used as integration step
        tol: relative mass left when the Euler quadrature stops
        iter_max: maximum number of Euler integration steps

    Returns:
        Tuple (e1, e2, e3, e11, idc)
    """
    T0 = csc_matrix(T0)
    T1 = csc_matrix(T1)
    dim = T0.shape[0]
    Q = (T0 + T1).tocsc()

    phi = _stationary(Q)
    pie = np.asarray(T1.T.dot(phi)).ravel()
    lam = float(np.sum(pie))
    pie = pie / lam

    if method.lower() == 'euler':
        dmax = np.max(np.abs(T0.diagonal()))
        dt = step_safety / dmax
        solve = lambda v: fes_map_euler(v, T0, dt, tol, iter_max)
    elif method.lower() == 'ssolve':
        negT0t = (-T0).transpose().tocsc()
        solve = lambda v: spsolve(negT0t, v)
    else:
        raise ValueError("Unknown method %s, use ssolve or euler." % method)

    v1 = solve(pie)
    v2 = solve(v1)
    v3 = solve(v2)
    v4 = solve(np.asarray(T1.T.dot(v2)).ravel())

    e1 = float(np.sum(v1))
    e2 = 2 * float(np.sum(v2))
    e3 = 6 * float(np.sum(v3))
    e11 = float(np.sum(v4))

    # equation (7), with pie*inv(Q+e*phi) from the rank-one update y*Q = pie-phi
    # under the normalization y*e = 1
    A = lil_matrix(Q.copy())
    A[:, dim - 1] = 1.0
    rhs = pie - phi
    rhs[dim - 1] = 1.0
    y = spsolve(A.tocsc().transpose().tocsc(), rhs)
    idc = 1 + 2 * (lam - float(np.sum(T1.T.dot(y))))

    return e1, e2, e3, e11, idc


def _stationary(Q) -> np.ndarray:
    """Stationary distribution of a sparse generator, by the replaced-column solve."""
    dim = Q.shape[0]
    A = lil_matrix(Q.copy())
    A[:, dim - 1] = 1.0
    b = np.zeros(dim)
    b[dim - 1] = 1.0
    p = spsolve(A.tocsc().transpose().tocsc(), b)
    return p / np.sum(p)


def fes_map_grid(n: int, nhead: int = NHEAD, ntail: int = NTAIL) -> np.ndarray:
    """
    Population levels at which the inter-departure MAP is evaluated.

    Fitting one MAP per population is wasteful because the processes of
    neighbouring populations are similar. The reference evaluates the first ten
    populations and ten further equispaced points.

    Args:
        n: largest population
        nhead: leading populations kept in full
        ntail: equispaced points above them

    Returns:
        Sorted array of populations to evaluate
    """
    if n <= nhead + ntail:
        return np.arange(1, n + 1)
    tail = np.round(np.linspace(nhead + 1, n, ntail)).astype(int)
    return np.unique(np.concatenate([np.arange(1, nhead + 1), tail]))


def fes_map_interp(x: np.ndarray, y: np.ndarray, xq: np.ndarray) -> np.ndarray:
    """
    Monotone piecewise cubic Hermite interpolation.

    Fritsch and Carlson slopes with the noncentered three-point endpoint rule of
    de Boor, so the interpolant never overshoots and a monotone sequence of
    throughputs stays monotone. The algorithm is written out rather than delegated
    to a library so that the MATLAB, Java, Python and C++ ports agree.

    Args:
        x: sample abscissae, strictly increasing
        y: sample values, one row per abscissa
        xq: query abscissae

    Returns:
        Interpolated values, one row per query point
    """
    x = np.asarray(x, dtype=float).ravel()
    xq = np.asarray(xq, dtype=float).ravel()
    y = np.asarray(y, dtype=float)
    if y.ndim == 1:
        y = y.reshape(-1, 1)
    if y.shape[0] != x.size:
        y = y.T
    n = x.size
    yq = np.zeros((xq.size, y.shape[1]))

    if n == 1:
        return np.tile(y, (xq.size, 1))

    h = np.diff(x)
    for c in range(y.shape[1]):
        v = y[:, c]
        delta = np.diff(v) / h
        d = np.zeros(n)
        if n == 2:
            d[:] = delta[0]
        else:
            for i in range(1, n - 1):
                if delta[i - 1] * delta[i] > 0:
                    w1 = 2 * h[i] + h[i - 1]
                    w2 = h[i] + 2 * h[i - 1]
                    d[i] = (w1 + w2) / (w1 / delta[i - 1] + w2 / delta[i])
            d[0] = _edge_slope(h[0], h[1], delta[0], delta[1])
            d[n - 1] = _edge_slope(h[n - 2], h[n - 3], delta[n - 2], delta[n - 3])

        for q in range(xq.size):
            t = xq[q]
            if t <= x[0]:
                i = 0
            elif t >= x[n - 1]:
                i = n - 2
            else:
                i = int(np.searchsorted(x, t, side='right') - 1)
                i = min(i, n - 2)
            s = t - x[i]
            c2 = (3 * delta[i] - 2 * d[i] - d[i + 1]) / h[i]
            c3 = (d[i] - 2 * delta[i] + d[i + 1]) / h[i] ** 2
            yq[q, c] = v[i] + s * (d[i] + s * (c2 + s * c3))
    return yq


def _edge_slope(h1: float, h2: float, del1: float, del2: float) -> float:
    """Noncentered three-point endpoint slope with the monotonicity clamps of de Boor."""
    d = ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2)
    if np.sign(d) != np.sign(del1):
        return 0.0
    if np.sign(del1) != np.sign(del2) and abs(d) > abs(3 * del1):
        return 3 * del1
    return d


def fes_map_aggregate(maps: List[MAP], servers: Sequence[float], n: int,
                      grid: Optional[Sequence[int]] = None, method: str = 'ssolve',
                      verbose: bool = False):
    """
    Aggregate a subnetwork into a load-dependent MAP flow-equivalent server.

    The first station seeds the flow-equivalent server; every further station is
    folded against the running server by building the inter-departure MAP of the
    resulting pair at each population level and fitting a MAP(2) to its first three
    moments and index of dispersion. Levels are evaluated on a grid and the four
    descriptors are interpolated between grid points; the MAP is refitted at every
    level from the interpolated descriptors, never interpolated entrywise.

    Args:
        maps: service process of each station, already scaled by its visit ratio
        servers: number of servers of each station, np.inf for a delay
        n: largest population the flow-equivalent server must serve
        grid: populations at which the inter-departure MAP is evaluated
        method: moment evaluation method, 'ssolve' or 'euler'
        verbose: report the fits taken at each fold

    Returns:
        Tuple (fes, info) with fes the per-level MAPs and info a dict holding
        throughput, moments, status and grid
    """
    m = len(maps)
    if m < 1:
        raise ValueError("At least one station is required.")
    if len(servers) != m:
        raise ValueError("One server count per station is required.")

    if grid is None:
        grid = fes_map_grid(n)
    grid = np.unique(np.append(np.asarray(grid, dtype=int).ravel(), n))
    grid = grid[(grid >= 1) & (grid <= n)]

    fes = fes_map_levels(maps[0], n, servers[0])
    moments = np.zeros((4, n))
    status = np.zeros(n, dtype=int)
    for k in range(n):
        moments[0, k] = map_moment(fes[k][0], fes[k][1], 1)
        moments[1, k] = map_moment(fes[k][0], fes[k][1], 2)
        moments[2, k] = map_moment(fes[k][0], fes[k][1], 3)
        moments[3, k] = map_idc(fes[k][0], fes[k][1])

    for i in range(1, m):
        gmom = np.zeros((4, grid.size))
        for g, k in enumerate(grid):
            T0, T1 = fes_map_interdeparture(maps[i], fes, int(k), [servers[i], 1])
            e1, e2, e3, e11, idc = fes_map_moments(T0, T1, method)
            gmom[:, g] = [e1, e2, e3, idc]

        if grid.size < n:
            moments = fes_map_interp(grid, gmom.T, np.arange(1, n + 1)).T
        else:
            moments = gmom

        new_fes = []
        for k in range(n):
            fit, st = map2_fit_idc(moments[0, k], moments[1, k], moments[2, k], moments[3, k])
            new_fes.append(fit)
            status[k] = st
        fes = new_fes

        if verbose:
            print("FES fold %d/%d: %d levels on a grid of %d, %d exact fits, %d fallbacks"
                  % (i + 1, m, n, grid.size, int(np.sum(status == 0)), int(np.sum(status > 0))))

    info = {
        'throughput': 1.0 / moments[0, :],
        'moments': moments,
        'status': status,
        'grid': grid,
    }
    return fes, info

def fes_map_solve(fes, think_map: MAP, n: int):
    """
    Solve the reduced model made of a delay and a load-dependent MAP
    flow-equivalent server.

    Closes the aggregation of Section 5.2.1 of the reference. Once a subnetwork
    has been replaced by the load-dependent MAP of fes_map_aggregate, the model
    left is a delay holding the think times and one station, which is a finite
    level-dependent quasi birth-death process: level k is the number of jobs held
    by the flow-equivalent server and n-k jobs are thinking. The chain is the same
    block bidiagonal pair used to measure the inter-departure times, now read as a
    generator, so the delay is a station scaled by the number of jobs it holds and
    the marked transitions are the arrivals into the flow-equivalent server.

    The think time may itself be a MAP, which is how Section 5.3.1 models a
    bounded flash crowd.

    Args:
        fes: flow-equivalent server, one MAP per level, from fes_map_aggregate
        think_map: think time process (Z0,Z1)
        n: number of jobs in the closed model

    Returns:
        Tuple (X, R, Q, pk): system throughput, response time of the aggregated
        subnetwork, mean jobs held by the aggregate, and the level distribution
    """
    if n < 1:
        raise ValueError("The population n must be at least 1.")

    fes_lev = fes_map_levels(fes, n)
    T0, T1 = fes_map_interdeparture(think_map, fes_lev, n, [np.inf, 1])
    Q = (T0 + T1).tocsc()
    dim = Q.shape[0]

    phi = _stationary(Q)
    x = float(np.sum(T1.T.dot(phi)))

    blk = dim // (n + 1)
    pk = np.zeros(n + 1)
    for k in range(n + 1):
        pk[k] = float(np.sum(phi[k * blk:(k + 1) * blk]))
    qn = float(np.dot(np.arange(n + 1), pk))
    rn = n / x - map_moment(think_map[0], think_map[1], 1)
    return x, rn, qn, pk


def fes_map_deaggregate(pk, demands, servers, is_delay):
    """
    Recover the per-station metrics of an aggregated subnetwork by conditioning on
    the population held by the flow-equivalent server.

    E[Y_i] = sum_k pk(k) Y_i(k) with Y_i(k) the metric of station i when the
    isolated subnetwork holds k jobs. This is the decomposition step of the
    hierarchical analysis of Chandy, Herzog and Woo, IBM J. Res. Dev. 19(1), 1975,
    exact for a product-form subnetwork. It is an approximation whenever the
    burstiness the flow-equivalent server carries also matters inside the
    subnetwork; the aggregate metrics of fes_map_solve do not rely on it.

    Args:
        pk: distribution of the jobs held by the aggregate, pk[k] = P(k jobs)
        demands: service demands of the isolated subnetwork
        servers: servers per station, np.inf for a delay
        is_delay: True where the station is a pure delay

    Returns:
        Tuple (Q, U, X, R) of per-station arrays
    """
    from ..pfqn import pfqn_mva

    demands = np.asarray(demands, dtype=float).ravel()
    servers = np.asarray(servers, dtype=float).ravel()
    is_delay = np.asarray(is_delay, dtype=bool).ravel()
    m = demands.size
    n = len(pk) - 1

    qn = np.zeros(m)
    un = np.zeros(m)
    xn = np.zeros(m)
    queue_idx = np.where(~is_delay)[0]
    delay_idx = np.where(is_delay)[0]
    z = float(np.sum(demands[delay_idx])) if delay_idx.size else 0.0

    for k in range(1, n + 1):
        if pk[k] <= 1e-14:
            continue
        # native pfqn_mva returns (XN, CN, QN, UN, RN, TN, AN), not the MATLAB
        # (XN, QN, UN, CN) order
        out = pfqn_mva(demands[queue_idx].reshape(-1, 1), np.array([float(k)]),
                       np.array([z]), servers[queue_idx])
        xk = float(np.ravel(out[0])[0])
        qk = np.ravel(np.asarray(out[2]))
        uk = np.ravel(np.asarray(out[3]))
        qn[queue_idx] += pk[k] * qk
        un[queue_idx] += pk[k] * uk
        xn[queue_idx] += pk[k] * xk
        for d in delay_idx:
            qn[d] += pk[k] * xk * demands[d]
            un[d] += pk[k] * xk * demands[d]
            xn[d] += pk[k] * xk

    rn = np.zeros(m)
    nz = xn > 1e-14
    rn[nz] = qn[nz] / xn[nz]
    return qn, un, xn, rn
