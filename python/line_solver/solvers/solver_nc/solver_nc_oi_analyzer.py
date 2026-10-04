"""
Exact normalizing-constant analysis of a closed queueing network that mixes
order-independent (OI) stations with ordinary BCMP product-form stations
(native Python port of solver_nc_oi_analyzer.m).

Supported stations:
  - OI stations (SchedStrategy.OI / PAS with an empty swap graph), analyzed by
    the balanced-fairness rank rate mu(supp n) (Bonald and Proutiere 2003);
  - any BCMP product-form station: infinite-server (delay, IS), processor
    sharing (PS), LCFS-PR, SIRO, and class-independent-rate FCFS, single- or
    multi-server, analyzed by the load-dependent BCMP weight table
      W_i(n) = (sum n)!/prod(n_r!) * prod_r D_{i,r}^{n_r} / prod_{k=1}^{sum n} beta_i(k),
    with D_{i,r} = V(i,r)/rate(i,r) and beta_i(k) = min(k, c) (k for IS).

The full-network normalizing-constant table G(P) over the lattice
0 <= P <= N is assembled by convolution: the OI stations and the aggregated
delay give the core table via pfqn_ncoi, then each BCMP queue is folded in.
The exact per-class mean queue length at any station follows from the OI
functional-server (FNC) identity pfqn_oi_fnc (Casale, QEST 2006):

    E[n_{i,r}] = ( sum_{0<=b<=N} Psi_{i,r}(b) G(N-b) ) / G(N) - 1.

General per-class visits at the OI stations are supported: they enter the
v-weighted balanced-fairness balance Phi^v(n)=(1/mu(n)) sum_r v_r Phi^v(n-e_r).
BCMP-station visits are arbitrary (folded into the demand).
"""

import time
import numpy as np

from ...api.pfqn import pfqn_ncoi, pfqn_oi_fnc, pfqn_oi_insvc

_BCMP_SCHED = ('PS', 'LCFSPR', 'FCFS', 'SIRO')


def _sched_name(sn, ist):
    """Scheduling-strategy name at station ist. Compare by name (not raw
    integer) because SchedStrategy is represented by more than one enum class
    across the codebase and OI is normalized to PAS in the struct."""
    sched = sn.sched.get(ist)
    return getattr(sched, 'name', None)


def _fcfs_rate_ok(sn, ist):
    """True unless station ist has class-dependent FCFS/SIRO rates (which are
    not product form)."""
    rr = np.asarray(sn.rates, dtype=float)[ist, :]
    nj = np.asarray(sn.njobs, dtype=float)
    rr = rr[np.isfinite(rr) & (nj > 0)]
    if rr.size == 0:
        return True
    return (rr.max() - rr.min()) <= 1e-9 * rr.max()


def nc_is_oi_model(sn) -> bool:
    """True when the model is a closed queueing network that contains at least
    one order-independent (OI) station and every other station is a BCMP
    product-form station: infinite server (delay), PS, LCFS-PR, SIRO, or
    class-independent-rate FCFS. The OI requirement keeps pure-BCMP networks on
    the standard (faster) normalizing-constant path."""
    if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
        return False
    has_oi = False
    for ist in range(sn.nstations):
        name = _sched_name(sn, ist)
        if name == 'INF':
            continue
        elif name == 'PAS' or name == 'OI':
            ind = int(sn.stationToNode[ist])
            npar = sn.nodeparam[ind] if (sn.nodeparam is not None and ind in sn.nodeparam) else None
            if not isinstance(npar, dict) or 'swapGraph' not in npar:
                return False
            sg = npar.get('swapGraph')
            if sg is None or np.any(np.asarray(sg, dtype=float) != 0):
                return False  # genuine pass-and-swap: not order-independent
            has_oi = True
        elif name in _BCMP_SCHED:
            if name in ('FCFS', 'SIRO') and not _fcfs_rate_ok(sn, ist):
                return False  # class-dependent FCFS/SIRO: not product form
        else:
            return False  # an unsupported (non-product-form) station
    return has_oi


def solver_nc_oi_analyzer(sn, options):
    """Return (Q, U, R, T, C, X, lG, runtime, iter, method)."""
    t_start = time.time()
    iterc = 1
    method = 'oi'

    M = sn.nstations
    K = sn.nclasses

    # ---- reject class switching (OI rank rates are per raw class) ----------
    for c in range(sn.nchains):
        if len(np.atleast_1d(sn.inchain[c])) > 1:
            raise RuntimeError('solver_nc_oi requires one class per chain (no class switching).')
    if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
        raise RuntimeError('solver_nc_oi requires a closed queueing network.')
    N = np.round(np.asarray(sn.njobs, dtype=float)).astype(int).ravel()

    # OI: rank-rate balanced-fairness station. INF: aggregated into delay Z.
    # Q: ordinary BCMP product-form station (PS/LCFS-PR/FCFS/SIRO).
    isOI = np.zeros(M, dtype=bool)
    isINF = np.zeros(M, dtype=bool)
    isQ = np.zeros(M, dtype=bool)
    svc = [None] * M
    for ist in range(M):
        ind = int(sn.stationToNode[ist])
        name = _sched_name(sn, ist)
        if name == 'INF':
            isINF[ist] = True
        elif name == 'PAS' or name == 'OI':
            npar = sn.nodeparam[ind] if (sn.nodeparam is not None and ind in sn.nodeparam) else None
            sg = npar.get('swapGraph') if isinstance(npar, dict) else None
            if sg is None or np.any(np.asarray(sg, dtype=float) != 0):
                raise RuntimeError('solver_nc_oi supports OI stations only (PAS with a non-empty swap graph is not order-independent).')
            isOI[ist] = True
            svc[ist] = npar.get('svcRateFun') if isinstance(npar, dict) else None
            if svc[ist] is None:
                raise RuntimeError('OI station %d has no service rate function; set it via setService(lambda c: ...).' % ist)
        elif name in _BCMP_SCHED:
            isQ[ist] = True
            if name in ('FCFS', 'SIRO') and not _fcfs_rate_ok(sn, ist):
                raise RuntimeError('Station %d has class-dependent FCFS/SIRO rates and is not product form; solver_nc_oi requires class-independent rates.' % ist)
        else:
            raise RuntimeError('solver_nc_oi supports only INF (delay), OI, PS, LCFS-PR, SIRO and class-independent FCFS stations.')

    # ---- per-class visits (chain == class); normalize to reference station -
    V = np.zeros((M, K))
    for r in range(K):
        c = int(np.where(np.asarray(sn.chains)[:, r])[0][0])  # chain carrying class r
        vis = sn.visits[c]                                     # (nstateful x nclasses)
        for ist in range(M):
            isf = int(sn.stationToStateful[ist])
            V[ist, r] = vis[isf, r]
        vref = V[int(sn.refstat[r]), r]
        if vref > 0:
            V[:, r] = V[:, r] / vref

    # ---- per-class demand and aggregated delay demand Z_r ------------------
    with np.errstate(divide='ignore', invalid='ignore'):
        ST = 1.0 / np.asarray(sn.rates, dtype=float)
    ST[~np.isfinite(ST)] = 0.0
    Z = np.zeros(K)
    for ist in np.where(isINF)[0]:
        for r in range(K):
            Z[r] += V[ist, r] * ST[ist, r]
    D = np.zeros((M, K))                       # per-class demand at BCMP queues
    for ist in np.where(isQ)[0]:
        for r in range(K):
            D[ist, r] = V[ist, r] * ST[ist, r]

    # OI-station class visit ratios feed the v-weighted balance; general
    # (non-unit) visits are supported. oivis[m] is the 1xR vector of oiList[m].
    oiList = list(np.where(isOI)[0])
    oivis = [np.asarray(V[ist, :], dtype=float) for ist in oiList]

    # ---- OI rank-rate handles on a per-class count vector ------------------
    rates = [_make_rank_rate(svc[ist]) for ist in oiList]

    # ---- population lattice ------------------------------------------------
    shp, stride, total = _oi_lattice(N)

    # ---- core normalizing-constant table (OI stations + aggregated delay) --
    # One call: the balanced-fairness convolution is a lattice convolution, so
    # its internal table already holds G(n) for every 0 <= n <= N on the same
    # stride used here. Re-calling it per population would cost a needless
    # factor total = prod_r (N_r+1).
    _, _, Gfull = pfqn_ncoi(Z, N, rates, oivis)

    # ---- fold the BCMP queueing stations by lattice convolution ------------
    qList = list(np.where(isQ)[0])
    for ist in qList:
        Wq = _oi_ld_table(D[ist, :], sn.nservers[ist], shp, total)
        Gfull = _oi_conv(Gfull, Wq, shp, stride, total)

    G = float(Gfull[total - 1])
    lG = np.log(G)

    # ---- per-class throughput X_r = G(N - e_r)/G(N) ------------------------
    X = np.zeros(K)
    for r in range(K):
        if N[r] > 0:
            er = np.zeros(K, dtype=int)
            er[r] = 1
            X[r] = Gfull[int(np.sum((N - er) * stride))] / G

    # ---- per-station per-class mean queue length via the FNC identity ------
    Q = np.zeros((M, K))
    for m, ist in enumerate(oiList):           # OI stations
        Phi = _oi_phi(rates[m], N, oivis[m])
        for r in range(K):
            if N[r] > 0:
                _, Psir, _ = pfqn_oi_fnc(Phi, N, _make_target(r))
                Q[ist, r] = _oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1
    for ist in qList:                          # BCMP queueing stations
        Wq = _oi_ld_table(D[ist, :], sn.nservers[ist], shp, total)
        for r in range(K):
            if N[r] > 0:
                _, Psir, _ = pfqn_oi_fnc(Wq, N, _make_target(r))
                Q[ist, r] = _oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1
    for ist in np.where(isINF)[0]:             # delay: Little's law
        for r in range(K):
            Q[ist, r] = X[r] * V[ist, r] * ST[ist, r]

    # ---- throughput, utilization, response time per station ----------------
    T = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    for ist in range(M):
        for r in range(K):
            T[ist, r] = X[r] * V[ist, r]
    for ist in np.where(isINF)[0]:
        U[ist, :] = Q[ist, :]                  # INF utilization convention
    for ist in qList:                          # BCMP queue: offered-load per server
        c = sn.nservers[ist]
        if not np.isfinite(c) or c <= 0:
            c = 1
        for r in range(K):
            U[ist, r] = X[r] * D[ist, r] / c
    for m, ist in enumerate(oiList):
        # see _kb/06-solver-catalog.md ("OI (order-independent) utilization")
        S = sn.nservers[ist]
        if not np.isfinite(S) or S <= 0:
            S = 1
        Phi = _oi_phi(rates[m], N, oivis[m])
        gins, _, _ = pfqn_oi_insvc(rates[m], N)
        for r in range(K):
            if N[r] > 0:
                _, Psir, _ = pfqn_oi_fnc(Phi, N, _make_insvc_target(gins, stride, r))
                U[ist, r] = (_oi_fnc_mean(Psir, Gfull, shp, stride, total) / G - 1) / S
    for ist in range(M):
        for r in range(K):
            if T[ist, r] > 0:
                R[ist, r] = Q[ist, r] / T[ist, r]

    C = np.zeros(K)                            # per-class system response time
    for r in range(K):
        if X[r] > 0:
            C[r] = N[r] / X[r]

    runtime = time.time() - t_start
    return Q, U, R, T, C, X, lG, runtime, iterc, method


def _make_rank_rate(fun):
    """OI rank rate on a per-class count vector n. The service function
    svcRateFun(c) takes an ordered microstate list; for an order-independent
    station it is permutation-invariant, hence a function of the multiset of
    present jobs (class multiplicities included, not merely the support).
    Evaluate it on a canonical microstate holding n_r copies of class r."""
    return lambda n: fun(np.repeat(np.arange(np.asarray(n).size),
                                   np.round(np.asarray(n)).astype(int)))


def _make_target(r):
    """Target function f(n) = n_r for the FNC balance construction."""
    return lambda n: float(np.asarray(n).ravel()[r])


def _make_insvc_target(gins, stride, r):
    """Target function f(n) = E[sir_r | n] for the FNC balance construction."""
    return lambda n: float(gins[int(np.sum(np.asarray(n).ravel() * stride)), r])


def _oi_lattice(N):
    """Column-major lattice descriptor for populations 0 <= n <= N."""
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size
    shp = (N + 1).astype(int)
    stride = np.ones(R, dtype=int)
    for d in range(1, R):
        stride[d] = stride[d - 1] * shp[d - 1]
    total = int(np.prod(shp))
    return shp, stride, total


def _oi_sub(i, shp):
    """Decode linear index i (0-based) to the subscript vector n."""
    R = shp.size
    n = np.zeros(R, dtype=int)
    li = i
    for d in range(R):
        n[d] = li % shp[d]
        li //= shp[d]
    return n


def _oi_phi(oirate, N, vis=None):
    """Forward balanced-fairness fill of the OI balance function over the
    lattice: Phi(0)=1, Phi(n) = (1/mu(n)) sum_{r: n_r>0} Phi(n - e_r)."""
    shp, stride, total = _oi_lattice(N)
    R = shp.size
    Phiv = np.zeros(total)
    for i in range(total):
        n = _oi_sub(i, shp)
        if n.sum() == 0:
            Phiv[i] = 1.0
            continue
        s = 0.0
        for r in range(R):
            if n[r] > 0:
                vr = 1.0 if vis is None else vis[r]
                s += vr * Phiv[i - stride[r]]
        Phiv[i] = s / oirate(n)
    if R == 1:
        return Phiv
    return Phiv.reshape(tuple(shp), order='F')


def _oi_ld_table(Dq, c, shp, total):
    """BCMP load-dependent weight table over the lattice:
      W(n) = (sum n)!/prod(n_r!) * prod_r D_r^{n_r} / prod_{k=1}^{sum n} beta(k),
    beta(k) = min(k, c) for a c-server queue (c=1 -> single server). W(0)=1."""
    Dq = np.asarray(Dq, dtype=float).ravel()
    R = shp.size
    if not np.isfinite(c) or c <= 0:
        c = 1
    c = int(c)
    W = np.zeros(total)
    from math import lgamma, log
    for i in range(total):
        n = _oi_sub(i, shp)
        tot = int(n.sum())
        logf = lgamma(tot + 1)
        ok = True
        for r in range(R):
            if n[r] > 0:
                if Dq[r] <= 0:
                    ok = False
                    break
                logf += n[r] * log(Dq[r]) - lgamma(n[r] + 1)
        if not ok:
            continue
        for k in range(1, tot + 1):
            logf -= log(min(k, c))
        W[i] = np.exp(logf)
    return W


def _oi_conv(Av, Bv, shp, stride, total):
    """Lattice convolution Cv(m) = sum_{0<=a<=m} Av(a) Bv(m-a) over 0..N."""
    R = shp.size
    subs = np.zeros((total, R), dtype=int)
    for i in range(total):
        subs[i, :] = _oi_sub(i, shp)
    Av = np.asarray(Av).ravel()
    Bv = np.asarray(Bv).ravel()
    Cv = np.zeros(total)
    for i in range(total):
        m = subs[i, :]
        acc = 0.0
        for j in range(i + 1):
            a = subs[j, :]
            if np.all(a <= m):
                acc += Av[j] * Bv[int(np.sum((m - a) * stride))]
        Cv[i] = acc
    return Cv


def _oi_fnc_mean(Psi, Gfull, shp, stride, total):
    """G^{+} = sum_{0<=b<=N} Psi(b) G(N-b): the FNC of the target station
    convolved against the full-network normalizing-constant table, at n = N."""
    N = shp - 1
    Psiv = np.asarray(Psi).ravel(order='F')
    Gfull = np.asarray(Gfull).ravel()
    val = 0.0
    for i in range(total):
        if Psiv[i] == 0:
            continue
        b = _oi_sub(i, shp)
        val += Psiv[i] * Gfull[int(np.sum((N - b) * stride))]
    return val
