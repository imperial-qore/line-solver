"""
Pass-and-swap (P&S) importance-sampling tools.

Native Python port of the MATLAB pass-and-swap routines:
    pfqn_pas_nc     - exact microstate normalizing constant G_C of one
                     communicating class, given the per-station placement order.
    pfqn_pas_is     - auto-normalized importance-sampling estimate of the
                     normalizing constant G_C and mean queue lengths of a
                     closed two-station P&S / order-independent tandem with a
                     swap graph (Casale, Comte and Dorsman, 2026).
    pas_placement  - isolated placement-order logic (precedence closure and the
                     ``placeable`` check that decides feasible orderings).
    pas_swap2order - derive the global placement-order DAG from a swap graph via
                     the minimal single-job-per-class reachable component.

Classes are 0-based throughout (native convention): a swap graph G is indexed by
0-based class indices and svcRateFun receives 0-based ordered class lists, matching
solver_nc_oi_analyzer._make_rank_rate.

References:
    Original MATLAB: matlab/src/api/pfqn/pfqn_pas_nc.m, pfqn_pas_is.m,
    pas_placement.m, pas_swap2order.m; Comte and Dorsman, arXiv:2009.12299 (2021).
"""

import numpy as np
from typing import Callable, List, Optional, Sequence, Tuple, Union


def _opt(options, name, default):
    """Read an option field from a dict or dataclass/object; None -> default."""
    if options is None:
        return default
    if isinstance(options, dict):
        val = options.get(name, default)
    else:
        val = getattr(options, name, default)
    return default if val is None else val


def pas_placement(H) -> Tuple[Optional[np.ndarray], Callable]:
    """Placement-order logic of a P&S / OI network with swap graph H.

    An ordering c is feasible iff it is non-decreasing w.r.t. H, i.e. class a
    never appears before class b whenever H[b, a] != 0 (Comte and Dorsman,
    2021); equivalently H[b, a] != 0 means b must precede a. The transitive
    closure P[i, j] = 1 iff i must precede j makes the full order explicit.

    Parameters
    ----------
    H : (R, R) swap-graph adjacency. Empty/all-zero -> no constraint.

    Returns
    -------
    (P, placeable) : P is the (R, R) precedence closure (or None if H is empty);
        placeable(x) returns the array of class indices drawable next given the
        remaining per-class count vector x (present classes with no remaining
        predecessor still to be placed).
    """
    if H is None or np.size(H) == 0:
        return None, (lambda x: np.where(np.asarray(x).ravel() > 0)[0])
    H = (np.asarray(H) != 0).astype(int)
    R = H.shape[0]

    # Transitive closure of the "must precede" relation. Iterate until fixpoint
    # (paths of any length), at most R rounds.
    P = H.copy()
    for _ in range(R):
        Pnext = ((P + (P @ H)) > 0).astype(int)
        if np.array_equal(Pnext, P):
            break
        P = Pnext

    def placeable(x):
        x = np.asarray(x).ravel()
        # class j placeable iff present and no still-present predecessor:
        # sum_i x_i P[i, j] == 0.
        return np.where((x @ P == 0) & (x > 0))[0]

    return P, placeable


def pfqn_pas_is(N: Sequence[int],
                mu: Union[Callable, List[Callable]],
                H=None,
                options=None) -> Tuple[float, float, np.ndarray]:
    """Importance-sampling estimate of G_C and mean queue lengths of a closed
    two-station P&S tandem with swap graph H.

    G_C = sum_{c in D} sum_{k=0}^{ell} Phi_1(c[:k]) Phi_2(reverse(c[k:])), where
    D is the set of orderings non-decreasing w.r.t. H and Phi_m(q) = prod_p
    1/mu_m(n(q[:p])), n(.) the per-class COUNT vector of the prefix (not its
    support: OI property P1 only makes mu permutation-invariant). Orderings are
    drawn from D by placing, at each step, a
    uniformly random placement-order-minimal present class (auto-normalized IS,
    notebook generator IS_3); E[xi] = G_C[xi]/G_C[1] reuses the same samples.
    Taking xi = class-r count in the prefix gives the station-1 mean queue
    length of class r.

    Parameters
    ----------
    N : (R,) closed population vector (macrostate), finite.
    mu : list of exactly two callables. mu[m](n) is the total OI rank rate of
        station m given the per-class occupancy (count) vector n
        (permutation-invariant, but NOT a function of supp(n) alone); this is
        the count-based svcRateFun of an OI/PAS node.
    H : (R, R) swap-graph adjacency (H[b, a] != 0 forbids a before b). Empty ->
        pure OI (all orderings feasible).
    options : dict or options object; fields ``samples`` (default 1e4),
        ``seed`` (optional), ``verbose`` (default False), ``qlen`` (default
        True). ``qlen=False`` estimates ONLY the normalizing constant: the
        per-class prefix counts and their coefficients are neither accumulated
        nor allocated, and Q comes back as zeros. The ordering is drawn from
        the same stream either way, so G is unchanged to the last bit -- this
        is for the callers that want G(N - e_r) and discard the rest.

    Returns
    -------
    (G, lG, Q) : G is the IS estimate of G_C, lG = log(G), Q is (2, R) with the
        station-1 mean queue lengths and Q[1] = N - Q[0], or zeros when
        ``qlen=False``.
    """
    if callable(mu):
        mu = [mu]
    mu = list(mu)
    if len(mu) != 2:
        raise ValueError('pfqn_pas_is models a two-station pass-and-swap tandem: '
                         'mu must have exactly two rate functions.')

    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size
    if np.any(~np.isfinite(N.astype(float))):
        raise ValueError('pfqn_pas_is requires finite (closed) populations.')
    if H is None or np.size(H) == 0:
        H = np.zeros((R, R))
    H = (np.asarray(H) != 0).astype(int)
    if H.shape != (R, R):
        raise ValueError('H must be a (R x R) swap-graph adjacency matrix.')

    nsamples = int(round(_opt(options, 'samples', 10000)))
    seed = _opt(options, 'seed', None)
    verbose = bool(_opt(options, 'verbose', False))
    want_qlen = bool(_opt(options, 'qlen', True))
    rng = np.random.default_rng(seed if seed is not None else None)

    ell = int(N.sum())
    _, placeable = pas_placement(H)

    if ell == 0:
        return 1.0, 0.0, np.zeros((2, R))

    # xi = [1, n_{1,0}, ..., n_{1,R-1}], or just [1] when the caller wants the
    # constant alone -- the prefix-count coefficients are the bulk of the
    # per-sample bookkeeping and G does not depend on them.
    n_coef = R + 1 if want_qlen else 1
    accum = np.zeros(n_coef)

    for s in range(nsamples):
        # ---- draw an ordering c from D (auto-normalized IS, generator IS_3) --
        x = N.copy()
        c = np.zeros(ell, dtype=int)
        logp = 0.0
        for ppos in range(ell):
            avail = placeable(x)
            na = avail.size
            if na == 0:
                raise RuntimeError('swap graph induces no feasible ordering '
                                   '(cyclic placement order).')
            pick = int(avail[rng.integers(na)])
            c[ppos] = pick
            logp -= np.log(na)
            x[pick] -= 1
        p_c = np.exp(logp)

        # ---- station-1 prefix balance and class counts up to each cut --------
        Phi1 = np.ones(ell + 1)
        cnt1 = np.zeros((ell + 1, R)) if want_qlen else None
        occ = np.zeros(R)
        phi = 1.0
        for k in range(1, ell + 1):
            cls = c[k - 1]
            occ[cls] += 1
            phi = phi / mu[0](occ)
            Phi1[k] = phi
            if want_qlen:
                cnt1[k, :] = occ

        # ---- station-2 reversed-suffix balance B2[k] at cut k ---------------
        # B2[k] = Phi_2(reverse(c[k:])); B2[ell] = 1 (empty suffix).
        B2 = np.ones(ell + 1)
        occ2 = np.zeros(R)
        phi = 1.0
        for pos in range(ell, 0, -1):   # pos = ell..1 -> reversed suffix c[ell..pos]
            cls = c[pos - 1]
            occ2[cls] += 1
            phi = phi / mu[1](occ2)
            B2[pos - 1] = phi

        # ---- split-convolution sample value for every coefficient -----------
        sv = np.zeros(n_coef)
        for k in range(ell + 1):
            w = Phi1[k] * B2[k]
            sv[0] += w                  # xi = 1
            if want_qlen and k > 0:
                sv[1:] += w * cnt1[k, :]  # xi = n_{1,r}
        accum += sv / p_c

        if verbose and nsamples >= 10 and (s + 1) % max(1, nsamples // 10) == 0:
            print('pfqn_pas_is: %d/%d samples' % (s + 1, nsamples))

    est = accum / nsamples
    G = float(est[0])
    lG = np.log(G) if G > 0 else -np.inf
    Q = np.zeros((2, R))
    if want_qlen:
        if G > 0:
            Q[0, :] = est[1:] / G
        Q[1, :] = N - Q[0, :]
    return G, lG, Q


def _pas_swap_local(c: List[int], p: int, G) -> Tuple[List[int], int]:
    """Pass-and-swap scan: the customer at position p completes; classes shift
    one step along the swap chain and the served slot is removed. Returns the
    new ordered list and the departing class. c and G are 0-based."""
    n = len(c)
    chain = [p]
    moving = c[p]
    cur = p
    while True:
        q = -1
        for j in range(cur + 1, n):
            if G is not None and np.size(G) > 0 and G[moving, c[j]] != 0:
                q = j
                break
        if q == -1:
            break
        chain.append(q)
        moving = c[q]
        cur = q
    dep = c[chain[-1]]
    tmp = list(c)
    for i in range(len(chain) - 1):
        tmp[chain[i + 1]] = c[chain[i]]
    del tmp[chain[0]]
    return tmp, dep


def pas_swap2order(swap, listRate, N0=None) -> np.ndarray:
    """Global placement-order DAG H of a closed two-station P&S tandem 1->2->1.

    With a non-empty swap graph the ordered chain is reducible; the recurrent
    communicating class is the set of splits of the orderings that are the
    linear extensions of a single placement partial order (Comte and Dorsman,
    2021). pfqn_pas_is samples those orderings from H, so it needs this GLOBAL
    order. The order is class-level (multiplicity-independent): enumerate the
    reachable class from the all-in-queue-1 single-job-per-class state; each
    reachable state (l1; l2) exposes the full ordering c = l1 + reverse(l2);
    then H[i, j] = 1 iff i precedes j in every such c (forced precedence).

    Parameters
    ----------
    swap : (R, R) swap graph, or list [G1, G2] of the two per-queue graphs
        (G[a, b] != 0 means class a chases class b). Empty/all-zero -> H = 0.
    listRate : list of two callables; listRate[m](c) is the total service rate
        of queue m on the ordered prefix c (0-based). Prunes zero-rate (non-
        head) completions.
    N0 : (R,) minimal probing population; defaults to ones(R).

    Returns
    -------
    H : (R, R) global placement-order DAG; H[i, j] = 1 iff i must precede j.
    """
    if not isinstance(swap, (list, tuple)):
        swap = [np.asarray(swap), np.asarray(swap)]
    else:
        swap = [np.asarray(g) if g is not None else None for g in swap]

    def _empty(g):
        return g is None or np.size(g) == 0 or not np.any(np.asarray(g) != 0)

    if N0 is None:
        R = int(np.asarray(swap[0]).shape[0])
        N0 = np.ones(R, dtype=int)
    N0 = np.round(np.asarray(N0, dtype=float)).astype(int).ravel()
    R = N0.size

    if _empty(swap[0]) and _empty(swap[1]):
        return np.zeros((R, R))

    # minimal single-job-per-class initial placement, all jobs at queue 1
    init_state = []
    for r in range(R):
        init_state.extend([r] * int(N0[r]))
    init = (tuple(init_state), tuple())

    def enc(st):
        return str(st)

    key = {enc(init)}
    class_states = []
    frontier = [init]
    while frontier:
        st = frontier.pop()
        class_states.append(st)
        for m in range(2):
            c = list(st[m])
            nxt = (m + 1) % 2
            g = swap[m]
            for p in range(len(c)):
                # marginal OI completion rate of the p-th customer; empty-prefix
                # rate is 0 (some svcRateFun return a nonzero constant on []).
                prev_rate = 0.0 if p == 0 else listRate[m](c[:p])
                rate = listRate[m](c[:p + 1]) - prev_rate
                if rate <= 1e-12:
                    continue
                cnew, dep = _pas_swap_local(c, p, g)
                stn = [list(st[0]), list(st[1])]
                stn[m] = cnew
                stn[nxt] = list(st[nxt]) + [dep]
                stn = (tuple(stn[0]), tuple(stn[1]))
                k = enc(stn)
                if k not in key:
                    key.add(k)
                    frontier.append(stn)

    # full orderings D: c = l1 + reverse(l2) for every reachable state
    D = []
    seen = set()
    for st in class_states:
        c = tuple(list(st[0]) + list(st[1])[::-1])
        if c not in seen:
            seen.add(c)
            D.append(c)

    # forced precedence: i before j iff in every c in D containing both, every
    # copy of i is older (earlier) than every copy of j.
    H = np.zeros((R, R), dtype=int)
    for a in range(R):
        for b in range(R):
            if a == b:
                continue
            both = False
            forced = True
            for c in D:
                pa = [idx for idx, v in enumerate(c) if v == a]
                pb = [idx for idx, v in enumerate(c) if v == b]
                if not pa or not pb:
                    continue
                both = True
                if not (max(pa) < min(pb)):
                    forced = False
                    break
            H[a, b] = 1 if (both and forced) else 0
    return H


def pfqn_pas_nc(Z: Optional[Sequence[float]],
                N: Sequence[int],
                mu: Optional[Union[Callable, List[Callable]]] = None,
                prec=None,
                options=None) -> Tuple[float, float]:
    """Normalizing constant G_C of one communicating class of a closed P&S network.

    With a non-empty swap graph the ordered-state chain is reducible (Comte and
    Dorsman, 2021, arXiv:2009.12299): the recurrent communicating classes are
    the placement-order-adhering sets and the product form
    pi(c) = prod_m Phi_m(c_m)/G_C holds per class. This routine returns G_C.

    It is a MICROSTATE routine: it walks the ordered chains position by
    position, because with a placement order the reachable set is a set of
    ORDERINGS that does not collapse onto the count lattice. For the plain OI
    case (empty ``prec``) the lattice does suffice and :func:`pfqn_ncoi`
    returns the same G at far lower cost.

    Method: build station M's chain head-first; appending class r at position
    k = sum(occ)+1 is admissible iff no class already placed at that station
    must come after r, and contributes 1/mu_M(occ+e_r); the chain may be
    finalized (recursing to station M-1) only when occ is a placement-order
    ideal at full multiplicity. Once every P&S station is peeled the residual
    population sits at the delay node with weight prod_r Z_r^{N_r}/N_r!.

    Cost: one node per feasible ordered prefix; with an empty ``prec`` that is
    ``sum_{b<=N} C(|b|+M-1, M-1) |b|!/prod_r b_r!``, factorial in sum(N). A
    placement order prunes the orderings, which is what makes the microstate
    walk affordable in the P&S case.

    Parameters
    ----------
    Z : (R,) think-time demand vector of the aggregated delay node; None/empty
        for a network with no delay station.
    N : (R,) closed population vector, finite.
    mu : list of callables, one per P&S station. ``mu[m](n)`` returns the total
        service rate of station m for the per-class occupancy vector n.
    prec : placement order. Either a list of M (R, R) precedence matrices, one
        per station, or a single (R, R) matrix broadcast to every station, with
        ``prec[m][i, j] != 0`` iff class i must be placed before class j at
        station m (the closure returned by :func:`pas_placement`, fed by the
        global DAG of :func:`pas_swap2order`). None or all-zero means no order:
        every ordering is feasible and G is the plain OI constant. NOTE the
        orientation: around a cycle each downstream station traverses its chain
        in the opposite direction, so downstream stations take the TRANSPOSE of
        the upstream order (``[P, P.T]`` for a two-station cycle). Passing the
        same P to both stations of a cycle silently returns a smaller, wrong G.
    options : accepted for signature parity; unused.

    Returns
    -------
    (G, lG) : normalizing constant G_C of the communicating class and its log.
    """
    N = np.round(np.asarray(N, dtype=float)).astype(int).ravel()
    R = N.size

    if mu is None:
        mu = []
    elif callable(mu):
        mu = [mu]
    else:
        mu = list(mu)
    M = len(mu)

    if Z is None or len(np.asarray(Z).ravel()) == 0:
        Z = np.zeros(R)
    Z = np.asarray(Z, dtype=float).ravel()
    if Z.size != R:
        raise ValueError('pfqn_pas_nc: Z and N must have the same number of classes.')
    if np.any(N < 0):
        raise ValueError('pfqn_pas_nc requires finite (closed) populations.')

    # Normalize the placement order to one boolean (R, R) matrix per station.
    if prec is None or (not isinstance(prec, (list, tuple)) and np.size(prec) == 0):
        precl = [np.zeros((R, R), dtype=bool) for _ in range(M)]
    elif isinstance(prec, (list, tuple)):
        precl = list(prec)
        if len(precl) == 1 and M > 1:
            precl = precl * M
        if len(precl) != M:
            raise ValueError('pfqn_pas_nc: prec must supply one precedence matrix per station.')
        precl = [np.zeros((R, R), dtype=bool) if p is None or np.size(p) == 0
                 else np.asarray(p) != 0 for p in precl]
    else:
        precl = [np.asarray(prec) != 0 for _ in range(M)]
    for p in precl:
        if p.shape != (R, R):
            raise ValueError('pfqn_pas_nc: each precedence matrix must be R x R.')

    G = _pas_nc_rec(Z, N.copy(), N, mu, precl, M, np.zeros(R, dtype=int))
    lG = float(np.log(G)) if G > 0 else -np.inf
    return float(G), lG


def _pas_nc_rec(Z, N, Norig, mu, prec, m, occ) -> float:
    """Head-peeling recursion over stations m, m-1, ..., 1 and the delay node."""
    R = N.size

    if m == 0:
        # Base case: the residual population sits at the delay node with
        # unnormalized weight prod_r Z_r^{N_r} / N_r!.
        from scipy.special import gammaln
        logf = 0.0
        for r in range(R):
            if N[r] > 0:
                if Z[r] <= 0:
                    return 0.0  # population but no delay demand: infeasible
                logf += N[r] * np.log(Z[r]) - gammaln(N[r] + 1)
        return float(np.exp(logf))

    # Step A: finalize station m's chain here and peel to station m-1, but only
    # if the accumulated occupancy is a placement-order ideal of station m.
    if _pas_nc_isideal(occ, prec[m - 1], Norig):
        G = _pas_nc_rec(Z, N, Norig, mu, prec, m - 1, np.zeros(R, dtype=int))
    else:
        G = 0.0

    # Step B: append one more class r at the next chain position of station m.
    active = mu[m - 1]
    precm = prec[m - 1]
    placed = occ > 0
    for r in range(R):
        if N[r] > 0 and not precm[r, placed].any():
            e_r = np.zeros(R, dtype=int)
            e_r[r] = 1
            mu_r = active(occ + e_r)
            if mu_r <= 0:
                continue
            Nm = N.copy()
            Nm[r] -= 1
            G += (1.0 / mu_r) * _pas_nc_rec(Z, Nm, Norig, mu, prec, m, occ + e_r)
    return G


def _pas_nc_isideal(occ, precm, Norig) -> bool:
    """occ is a placement-order ideal at full multiplicity: i prec j and occ[j]>0 => occ[i]==Norig[i]."""
    R = occ.size
    for i in range(R):
        for j in range(R):
            if precm[i, j] and occ[j] > 0 and occ[i] < Norig[i]:
                return False
    return True
