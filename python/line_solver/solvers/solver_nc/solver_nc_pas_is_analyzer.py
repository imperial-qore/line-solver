"""
Importance-sampling (IS) normalizing-constant analysis of a closed two-station
pass-and-swap (P&S) tandem with a non-empty swap graph (native Python port of
solver_nc_pas_is_analyzer.m; Casale, Comte and Dorsman, 2026).

With a genuine swap graph the ordered-state chain is reducible (Comte and
Dorsman, 2021, arXiv:2009.12299); the recurrent communicating class carries a
per-class product form pi(c) = Phi_1(c_1) Phi_2(c_2)/G_C whose constant G_C is
estimated by the auto-normalized IS routine pfqn_pas_is. This is the Monte-Carlo
counterpart of solver_nc_oi_analyzer for the case that the exact OI convolution
does not apply (non-empty swap graph).

Station 0 is the upstream P&S queue (prefix of the ordering), station 1 the
downstream queue (reversed suffix); the swap graph is a class-level property.
Both stations must be OI/P&S with unit per-class visits. Mean per-class queue
lengths come directly from pfqn_pas_is; per-class throughput uses the balanced-
fairness ratio X_r = G(N - e_r)/G(N) with common random numbers across the N and
N-e_r runs, and utilization/response time follow the solver_nc_oi conventions.
"""

import time
import numpy as np

from ...api.pfqn.pas import pfqn_pas_is, pas_swap2order


def _sched_name(sn, ist):
    """Scheduling-strategy name at station ist (compare by name, not raw int)."""
    sched = sn.sched.get(ist) if hasattr(sn.sched, 'get') else sn.sched[ist]
    return getattr(sched, 'name', None)


def _node_param(sn, ind):
    if sn.nodeparam is None:
        return None
    try:
        return sn.nodeparam[ind]
    except (KeyError, IndexError, TypeError):
        return None


def nc_is_pas_model(sn) -> bool:
    """True when the model is a closed two-station P&S tandem: both stations are
    OI/PAS and there is no other station. The swap graph may be empty, since an
    order-independent queue is exactly the P&S specialization with an empty swap
    graph and pfqn_pas_is with H=0 samples all microstates (the OI case).

    Used to bind the importance-sampling selectors ('is', and 'sampling' which
    maps to 'is' when OI/PAS stations are present). A P&S tandem with a
    NON-EMPTY swap graph is reducible and lies outside the exact path of
    nc_is_oi_model, so it also takes this analyzer on 'default'; a pure-OI
    tandem on 'default'/'exact' is caught earlier by the exact OI analyzer."""
    if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
        return False
    if sn.nstations != 2:
        return False
    for ist in range(sn.nstations):
        name = _sched_name(sn, ist)
        if name not in ('PAS', 'OI'):
            return False
        ind = int(sn.stationToNode[ist])
        npar = _node_param(sn, ind)
        if not isinstance(npar, dict):
            return False
        if npar.get('svcRateFun') is None:
            return False   # no OI rank-rate function: cannot evaluate the balance function
    return True


#: Count vectors kept per rank-rate handle. Each entry is one small index
#: array; the cap bounds the memo on models whose population makes the set of
#: reachable count vectors large, and a miss past it simply rebuilds.
_MICROSTATE_MEMO_LIMIT = 1 << 17


def _make_rank_rate(fun):
    """OI rank rate on a per-class count vector n. svcRateFun(c) takes an
    ordered microstate list; for an OI station it is permutation-invariant, so
    evaluate it on a canonical microstate holding n_r copies of class r.

    The microstate is MEMOIZED on the count vector, because the importance
    sampler walks the same prefix occupancies over and over: pfqn_pas_is calls
    this handle 2*(ell+1) times per sampled ordering, and rebuilding the
    representative each time made the arange/repeat/round/astype quartet -- not
    the balance function it exists to evaluate -- the single largest cost of
    SolverNC method='sampling'. The value depends on nothing but the counts, so
    the memo returns what the rebuild would have returned, bit for bit.
    """
    memo = {}

    def rank_rate(n):
        counts = np.asarray(n)
        key = counts.tobytes()
        micro = memo.get(key)
        if micro is None:
            micro = np.repeat(np.arange(counts.size),
                              np.round(counts).astype(int))
            if len(memo) < _MICROSTATE_MEMO_LIMIT:
                memo[key] = micro
        return fun(micro)

    return rank_rate


def solver_nc_pas_is_analyzer(sn, options):
    """Return (Q, U, R, T, C, X, lG, runtime, iter, method)."""
    t_start = time.time()
    iterc = 1
    method = 'sampling'   # P&S specialization of the importance-sampling NC method

    M = sn.nstations
    K = sn.nclasses

    # ---- reject class switching (P&S rank rates are per raw class) ---------
    for c in range(sn.nchains):
        if len(np.atleast_1d(sn.inchain[c])) > 1:
            raise RuntimeError('solver_nc_pas_is requires one class per chain (no class switching).')
    if np.any(np.isinf(np.asarray(sn.njobs, dtype=float))):
        raise RuntimeError('solver_nc_pas_is requires a closed queueing network.')
    if M != 2:
        raise RuntimeError('solver_nc_pas_is models a two-station pass-and-swap tandem (got %d stations).' % M)
    N = np.round(np.asarray(sn.njobs, dtype=float)).astype(int).ravel()

    # ---- classify stations and read swap graph ----------------------------
    svc = [None] * M
    swapG = [None] * M
    for ist in range(M):
        name = _sched_name(sn, ist)
        if name not in ('PAS', 'OI'):
            raise RuntimeError('solver_nc_pas_is requires both stations to be OI/PAS (station %d is not).' % ist)
        ind = int(sn.stationToNode[ist])
        npar = _node_param(sn, ind)
        if not isinstance(npar, dict):
            raise RuntimeError('station %d has no OI/PAS node parameters.' % ist)
        svc[ist] = npar.get('svcRateFun')
        if svc[ist] is None:
            raise RuntimeError('OI/PAS station %d has no service rate function; set it via setService(lambda c: ...).' % ist)
        swapG[ist] = npar.get('swapGraph')

    # see _kb/06-solver-catalog.md (NC: "Analyzer routing order") for the
    # global placement-order DAG derivation from the swap graph
    G1 = swapG[0] if swapG[0] is not None else np.zeros((K, K))
    G2 = swapG[1] if swapG[1] is not None else np.zeros((K, K))
    H = pas_swap2order([np.asarray(G1), np.asarray(G2)], [svc[0], svc[1]], np.ones(K, dtype=int))

    # ---- per-class visits (chain == class); require unit visits -----------
    V = np.zeros((M, K))
    for r in range(K):
        c = int(np.where(np.asarray(sn.chains)[:, r])[0][0])
        vis = sn.visits[c]
        for ist in range(M):
            isf = int(sn.stationToStateful[ist])
            V[ist, r] = vis[isf, r]
        vref = V[int(sn.refstat[r]), r]
        if vref > 0:
            V[:, r] = V[:, r] / vref
    for ist in range(M):
        for r in range(K):
            if N[r] > 0 and abs(V[ist, r] - 1) > 1e-9:
                raise RuntimeError('solver_nc_pas_is requires unit per-class visits (station %d, class %d, V=%g).' % (ist, r, V[ist, r]))

    # ---- OI rank-rate handles on a per-class count vector -----------------
    mu = [_make_rank_rate(svc[0]), _make_rank_rate(svc[1])]

    # ---- IS options; fix a base seed so the N and N-e_r runs share randoms -
    samples = getattr(options, 'samples', None)
    if samples is None or samples == 0:
        iter_max = getattr(options, 'iter_max', None)
        samples = iter_max if (iter_max is not None and iter_max > 1) else 10000
    seed = getattr(options, 'seed', None)
    if seed is None:
        seed = 23456
    verbose = bool(getattr(options, 'verbose', False))
    isopt = {'samples': int(samples), 'seed': int(seed), 'verbose': verbose}

    # ---- normalizing constant and mean queue lengths at population N ------
    G, lG, Qpas = pfqn_pas_is(N, mu, H, isopt)

    Q = np.zeros((M, K))
    Q[0, :] = Qpas[0, :]
    Q[1, :] = Qpas[1, :]

    # ---- per-class throughput X_r = G(N - e_r)/G(N) (common random numbers)-
    # Only the constant is read here, so these runs skip the prefix-count
    # coefficients: same stream, same G, none of the queue-length bookkeeping.
    isopt_G = dict(isopt)
    isopt_G['qlen'] = False
    X = np.zeros(K)
    for r in range(K):
        if N[r] > 0:
            er = np.zeros(K, dtype=int)
            er[r] = 1
            Gr, _, _ = pfqn_pas_is(N - er, mu, H, isopt_G)
            if G > 0:
                X[r] = Gr / G

    # ---- throughput, utilization, response time ---------------------------
    T = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    for ist in range(M):
        for r in range(K):
            T[ist, r] = X[r] * V[ist, r]
    for ist in range(M):
        S = sn.nservers[ist]
        if not np.isfinite(S) or S <= 0:
            S = 1
        for r in range(K):
            if N[r] > 0:
                er = np.zeros(K, dtype=int)
                er[r] = 1
                muR = mu[ist](er)          # rank rate with only class r present
                if muR > 0:
                    U[ist, r] = T[ist, r] / muR / S
    for ist in range(M):
        for r in range(K):
            if T[ist, r] > 0:
                R[ist, r] = Q[ist, r] / T[ist, r]

    C = np.zeros(K)
    for r in range(K):
        if X[r] > 0:
            C[r] = N[r] / X[r]

    runtime = time.time() - t_start
    return Q, U, R, T, C, X, lG, runtime, iterc, method
