"""
Background modulating chain (bgchain) method for MIXED and CLOSED queueing
networks.

The closed population vector of such a network is a finite continuous-time
Markov chain in its own right: it is the only part of the model whose state
space is bounded. This method solves it exactly (given the capacity the open
work leaves free) and hands the open classes a station-local Markovian
ENVIRONMENT read off that chain, so each open station becomes a level-dependent
QBD whose phase carries the number of closed jobs competing for its server. The
two halves meet at a fixed point on the closed capacity share.

  1. background chain   the closed population vector over the stations the
                        closed classes visit (``bgchain_ctmc``)
  2. environment        that chain lumped onto the closed occupancy of one
                        station (``bgchain_env``)
  3. open station       a MAP/PH/c queue modulated by that environment, solved
                        as a level-dependent QBD (``bgchain_station``)
  4. fixed point        the capacity share feeds step 1 and closes

TAGGED-CLASS ITERATION. Step 1 is a population process of dimension (closed
chains) x (stations), so its state space is exponential in the number of closed
chains. With R > 1 closed chains the method keeps ONE chain free at a time: the
tagged chain r is carried exactly, the other R-1 are replaced by flow-equivalent
aggregate classes whose population is their total and whose service time and
routing at each station are their throughput-weighted means (Chandy-Herzog-Woo
aggregation). Every chain takes its turn as the tagged one and reads its own
metrics off the chain it is exact in; the open-class results are averaged over
the passes.

HOW MUCH TO AGGREGATE is options.config.bgaggr, the number G of aggregate
classes; the background chain then carries 1 + G. G = 1 is the classic
tagged/aggregate pair and the default, so the chain stays two-class whatever R
is; G >= R-1 aggregates nothing, carries every closed chain exactly, and answers
in ONE pass instead of solving the same chain R times. Passing R reaches that, so
asking for no aggregation needs no magic value. The cost is the state space, the
product over the 1 + G classes of nchoosek(N_b + Mc - 1, Mc - 1), capped by
bgstates_max.

WHICH CHAINS SHARE A GROUP is decided by similarity of per-station SERVICE
DEMAND. An aggregate carries the flow-weighted mean of its members' service times
and routing, so it is exact when they place the same demand at every station and
distorts in proportion to how far apart they are; grouping the demand-similar
chains together keeps the aggregation where it is harmless and away from the
chains it would misrepresent.

EXACTNESS, as measured against SolverCTMC and exact MVA on mixed models of two
to four stations: PS or INF with ANY service law (exponential, Erlang, HyperExp,
Coxian), any number of servers, Poisson or MAP arrivals and one to four closed
chains agree to 4-5 significant digits, as does FCFS with class-independent
rates. FCFS with class-DEPENDENT rates keeps the closed queue lengths within ~1%
while the open queue length reads 14-20% low, because the server is held here in
random order rather than head-of-line.

PS is INSENSITIVE to the service law beyond its mean, and the method honours that
rather than approximating it: at a PS station the open service is replaced by the
exponential of the same mean before the QBD is built. This QBD tracks ONE service
phase for the whole station, so carrying the phase-type there makes the open
queue length inherit the SCV-sensitivity of an M/PH/1 FCFS queue -- measured, a
HyperExp of SCV 4 read 21% high where the exact answer is the exponential one to
five digits. At an FCFS station the service law IS carried, and the background
chain reads only the MEAN closed service time, exact under PS by the same
insensitivity and a first-moment surrogate under FCFS.

References:
    Original MATLAB: matlab/src/solvers/MAM/solver_mam_bgchain.m
"""

import time
from typing import Optional, Tuple

import numpy as np

from . import MAMAlgorithm, MAMResult
from ....api.mam import LdqbdOptions
from ....api.mam.ldqbd import ldqbd
from ....api.mam.map_analysis import map_pie, map_scale
from ....api.mam.mmap_ops import mmap_super_safe
from ....api.mc.ctmc import ctmc_makeinfgen, ctmc_solve
from ....api.qsys.retrial import _proc_to_d0d1
from ....api.sn.demands import sn_get_demands_chain
from ....api.sn.transforms import sn_rt_stations
from ....api.state import space_closed_single
from ....constants import GlobalConstants
from ....lang.base import SchedStrategy


def _sched_name(sn, i):
    """Discipline of station i by NAME: the enum's integer value differs across
    codebases, so never compare the raw integer."""
    sched = sn.sched
    s = sched.get(i, None) if isinstance(sched, dict) else sched[i]
    if s is None:
        return ''
    return s.name if hasattr(s, 'name') else str(s)


def _block(m, N):
    """Block of one background class: the ways to place N jobs over the m
    stations that class visits.

    This is api.state.space_closed_single, the lattice primitive the CTMC solver
    enumerates a closed population over, so the row order and the row count are
    the reference's rather than this module's; only the support differs, being
    the class's own stations rather than all of them."""
    return space_closed_single(m, N)


def bgchain_ctmc(Nb, STb, Pb, sched_names, nservers, cshare, supp=None, bgstates_max=20000):
    """Build and solve the background modulating chain.

    The chain carries B = 1 or 2 background classes. B = 1 is the single closed
    chain of the model; B = 2 is the tagged/aggregate pair.

    A station holding e closed jobs serves background class b at rate
    ``n[i,b]/STb[i,b]`` when it is an infinite server, and
    ``cshare[i,e] * (n[i,b]/e) / STb[i,b]`` otherwise, splitting the capacity
    the closed jobs hold over the background classes in proportion to their
    counts. That is exact under PS and is the random-order surrogate under FCFS.

    The open classes enter ONLY through ``cshare``, which is what makes this a
    MODULATING chain rather than a joint model. The exchanged quantity is the
    SHARE, already averaged over the open occupancy, and not the mean open
    occupancy itself: e/(e+k) is convex in k, so rebuilding the share from a
    mean k would bias the closed service rate downwards by Jensen's inequality,
    and the closed throughput with it.

    The block of one background class is ``api.state.space_closed_single``, the
    lattice primitive the CTMC solver enumerates a closed population over, and
    the joint space is their cartesian product in class-major order. Only the
    SUPPORT differs: the columns are the class's own stations rather than every
    station.

    Args:
        Nb: population of each background class, length B in {1,2}
        STb: (Mc, B) mean service time per station per background class
        Pb: list of B row-stochastic (Mc, Mc) routing matrices
        sched_names: discipline name of each station of the chain's support
        nservers: servers of each station of the chain's support
        cshare: (Mc, Nmax+1) mean number of servers the e closed jobs hold
        supp: (Mc, B) boolean; supp[i,b] is True when station i is on the route
            of background class b. NOT an optimization -- a chain that never
            visits a station cannot hold jobs there, and enumerating the union of
            every chain's stations puts probability on unreachable configurations
            that also ABSORB, because the chain's routing matrix has a zero row at
            an unvisited station which row-normalizes to a self-loop. The
            generator turns reducible and population conservation silently fails.
            None means every station is on every class's route
        bgstates_max: cap on the number of chain states

    Returns:
        dict with space, totocc, pi, Q, QLen, Tput, Ubusy, nstates
    """
    STb = np.asarray(STb, dtype=float)
    Mc = STb.shape[0]
    B = len(Nb)

    if supp is None:
        supp = np.ones((Mc, B), dtype=bool)
    supp = np.asarray(supp, dtype=bool)
    # Each class is enumerated over ITS OWN stations only; see the supp arg.
    sp, compb, idxb, nst = [], [], [], []
    for b in range(B):
        idx = np.where(supp[:, b])[0]
        if idx.size == 0:
            if Nb[b] > 0:
                raise RuntimeError('Background class %d holds %d jobs but visits no station.'
                                   % (b + 1, Nb[b]))
            idx = np.array([0], dtype=int)
        comp = _block(len(idx), int(Nb[b]))
        exp = np.zeros((comp.shape[0], Mc), dtype=int)
        exp[:, idx] = comp
        idxb.append(idx)
        compb.append(comp)
        sp.append(exp)
        nst.append(comp.shape[0])
    nstates = int(np.prod(nst))
    if nstates > bgstates_max:
        raise RuntimeError(
            'The background chain of this model has %d states, above the limit of %d. The chain '
            'enumerates the closed-class population vector over the %d stations the closed classes '
            'visit, so its size grows as nchoosek(N+Mc-1,Mc-1) per class, and it carries %d '
            'classes. Lower options.config["bgaggr"] to aggregate more of the closed chains into '
            'fewer classes, raise options.config["bgstates_max"] to solve it anyway, or reduce the '
            'closed populations.' % (nstates, bgstates_max, Mc, B))

    # Per-class transition targets: tgt[b][s, i*Mc+j] is the class-b state
    # reached from s when one job moves from station i to station j, or -1.
    tgt = []
    for b in range(B):
        comp = compb[b]
        idx = idxb[b]
        mb = len(idx)
        radix = (int(Nb[b]) + 1) ** np.arange(mb)
        keys = comp @ radix
        order = np.argsort(keys)
        keys_sorted = keys[order]
        t = -np.ones((nst[b], Mc * Mc), dtype=int)
        for ii in range(mb):
            movable = comp[:, ii] > 0
            if not np.any(movable):
                continue
            for jj in range(mb):
                if jj == ii:
                    continue
                cand = comp[movable].copy()
                cand[:, ii] -= 1
                cand[:, jj] += 1
                key = cand @ radix
                pos = np.searchsorted(keys_sorted, key)
                pos = np.clip(pos, 0, len(keys_sorted) - 1)
                hit = keys_sorted[pos] == key
                col = -np.ones(nst[b], dtype=int)
                col[np.where(movable)[0][hit]] = order[pos[hit]]
                t[:, idx[ii] * Mc + idx[jj]] = col
        tgt.append(t)

    # Joint space, class 0 outermost
    strideb = np.ones(B, dtype=int)
    for b in range(B):
        strideb[b] = int(np.prod(nst[b + 1:])) if b + 1 < B else 1
    space = np.zeros((nstates, Mc, B), dtype=int)
    subidx = np.zeros((nstates, B), dtype=int)
    for s in range(nstates):
        rem = s
        for b in range(B - 1, -1, -1):
            subidx[s, b] = rem % nst[b]
            rem //= nst[b]
        for b in range(B):
            space[s, :, b] = sp[b][subidx[s, b]]
    totocc = space.sum(axis=2)

    mu = np.zeros((Mc, B))
    for b in range(B):
        for i in range(Mc):
            if np.isfinite(STb[i, b]) and STb[i, b] > 0:
                mu[i, b] = 1.0 / STb[i, b]

    Q = np.zeros((nstates, nstates))
    rate_full = np.zeros((nstates, Mc, B))
    cap_busy = np.zeros((nstates, Mc))
    cshare = np.asarray(cshare, dtype=float)
    ngrid = cshare.shape[1]
    for s in range(nstates):
        for i in range(Mc):
            eclosed = int(totocc[s, i])
            if eclosed == 0:
                continue
            if sched_names[i] == 'INF':
                held = float(eclosed)
            else:
                held = cshare[i, min(eclosed, ngrid - 1)]
            if held <= 0:
                continue
            cap_busy[s, i] = held
            for b in range(B):
                if space[s, i, b] == 0 or mu[i, b] == 0:
                    continue
                r = held * (space[s, i, b] / eclosed) * mu[i, b]
                rate_full[s, i, b] = r
                P = Pb[b]
                for j in range(Mc):
                    if j == i or P[i, j] <= 0:
                        continue
                    tsub = tgt[b][subidx[s, b], i * Mc + j]
                    if tsub < 0:
                        continue
                    sdest = s + (tsub - subidx[s, b]) * strideb[b]
                    Q[s, sdest] += r * P[i, j]

    Q = ctmc_makeinfgen(Q)
    if nstates == 1:
        pi = np.array([1.0])
    else:
        pi = np.asarray(ctmc_solve(Q), dtype=float).ravel()
    pi = np.maximum(pi, 0.0)
    tot = pi.sum()
    if tot > 0:
        pi = pi / tot

    QLen = np.zeros((Mc, B))
    Tput = np.zeros((Mc, B))
    Ubusy = np.zeros((Mc, B))
    for b in range(B):
        QLen[:, b] = pi @ space[:, :, b]
        Tput[:, b] = pi @ rate_full[:, :, b]
    for i in range(Mc):
        if sched_names[i] == 'INF':
            Ubusy[i, :] = QLen[i, :]
        else:
            occ = totocc[:, i].astype(float)
            nz = occ > 0
            for b in range(B):
                share = np.zeros(nstates)
                share[nz] = space[nz, i, b] / occ[nz]
                Ubusy[i, b] = float(pi @ (cap_busy[:, i] * share)) / nservers[i]

    return {'space': space, 'totocc': totocc, 'pi': pi, 'Q': Q,
            'QLen': QLen, 'Tput': Tput, 'Ubusy': Ubusy, 'nstates': nstates}


def bgchain_states(sn, config=None):
    """Number of states of the background-chain CTMC, WITHOUT building it.

    Mirrors mam_bgchain_states.m. The size is what decides whether bgchain is
    affordable and bgchain_ctmc only discovers it after the partition is fixed,
    so the default-method chooser needs it up front. The count follows the
    partition BgchainAlgorithm.solve uses: a pass carries the tagged closed
    chain as background class 0 and the demand-similar groups of the other
    closed chains as classes 1..G, each enumerating the compositions of its
    population over the stations its members visit. Merging two chains onto the
    UNION of their supports can raise the count as easily as lower it, so the
    passes are enumerated rather than bounded and the largest returned: that is
    the one bgchain_ctmc would refuse.

    Returns 0 when bgchain does not apply to the model at all.
    """
    if config is None:
        config = {}
    elif not isinstance(config, dict):
        config = dict(getattr(config, '__dict__', {}))
    bgaggr_opt = int(config.get('bgaggr', 1))

    M, C = int(sn.nstations), int(sn.nchains)
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    dem = sn_get_demands_chain(sn)
    Vchain = np.asarray(dem.Vchain, dtype=float)
    Lchain = np.asarray(dem.Lchain, dtype=float)
    Nchain = np.asarray(dem.Nchain, dtype=float).ravel()

    inchain = {c: np.asarray(sn.inchain[c], dtype=int).ravel() for c in range(C)}
    closed_chains = [c for c in range(C)
                     if not bool(np.any(np.isinf(njobs[inchain[c]]))) and Nchain[c] > 0]
    R = len(closed_chains)
    if R == 0:
        return 0.0

    inbg = np.zeros(M, dtype=bool)
    for c in closed_chains:
        inbg |= Vchain[:, c] > GlobalConstants.Zero
    cst = np.where(inbg)[0]
    Mc = len(cst)
    if Mc == 0:
        return 0.0

    naggr = min(max(bgaggr_opt, 1), max(R - 1, 1))
    no_aggr = (R == 1) or (naggr >= R - 1)
    passes = [closed_chains[0]] if no_aggr else closed_chains

    def _binomial(n, k):
        # nchoosek in floating point, so a chain far above any usable size still compares
        if k < 0 or k > n:
            return 0.0
        kk = min(k, n - k)
        acc = 1.0
        for i in range(1, int(kk) + 1):
            acc = acc * (n - kk + i) / i
        return acc

    worst = 0.0
    for r in passes:
        if no_aggr:
            members = [[c] for c in closed_chains]
        else:
            others = [c for c in closed_chains if c != r]
            grp = bgchain_groups(Lchain[np.ix_(cst, others)], naggr)
            members = [[r]]
            for g in range(naggr):
                members.append([o for oi, o in enumerate(others) if grp[oi] == g])

        n = 1.0
        for mem in members:
            Nb = int(round(sum(Nchain[c] for c in mem)))
            supp = np.zeros(Mc, dtype=bool)
            for o in mem:
                supp |= Vchain[cst, o] > GlobalConstants.Zero
            m = int(supp.sum())
            if m == 0:
                m = 1   # an empty class still needs one slot to be indexed by
            n *= _binomial(Nb + m - 1, m - 1)
            if not np.isfinite(n):
                return float('inf')
        worst = max(worst, n)
    return worst


def bgchain_groups(D, G):
    """Group the columns of D into G clusters by per-station SERVICE DEMAND.

    WHY DEMAND IS THE RIGHT CRITERION. The aggregate that replaces a group
    carries the flow-weighted mean of its members' service times and routing, so
    the group aggregates EXACTLY when its members place the same demand at every
    station and distorts both quantities in proportion to how far apart they are.
    The distance is therefore the symmetric relative L1 gap between the demand
    vectors,

        dist(a,b) = sum_i |D[i,a] - D[i,b]| / ((sum_i D[i,a] + sum_i D[i,b])/2),

    which is scale-relative rather than absolute: it separates two chains whose
    demand PROFILE across the stations differs and two chains whose profile
    agrees but whose magnitude does not, and being dimensionless it groups a
    model the same way whatever its time unit.

    WHY COMPLETE LINKAGE. The clustering is agglomerative from singletons,
    merging at each step the pair of clusters whose WORST member-to-member
    distance is smallest. The aggregation error inside a group is driven by its
    worst mismatch and not by its average one, so complete linkage is the
    criterion that bounds what the aggregation actually costs.

    DETERMINISM. Ties are broken by the lexicographically smallest pair of
    cluster indices and the groups are relabelled by their smallest member, so
    the same input gives the same grouping in every codebase.

    Args:
        D: (Mc, n) per-station demand, one column per chain
        G: number of groups wanted, clamped to [1, n]

    Returns:
        (n,) array of group indices in 0..G-1
    """
    D = np.atleast_2d(np.asarray(D, dtype=float))
    Mc, n = D.shape
    grp = np.zeros(n, dtype=int)
    if n == 0:
        return grp
    G = int(min(max(round(G), 1), n))

    tot = D.sum(axis=0)
    dist = np.zeros((n, n))
    for a in range(n):
        for b in range(a + 1, n):
            den = (tot[a] + tot[b]) / 2.0
            d = float(np.abs(D[:, a] - D[:, b]).sum() / den) if den > GlobalConstants.Zero else 0.0
            dist[a, b] = d
            dist[b, a] = d

    clusters = [[a] for a in range(n)]
    active = [True] * n
    nactive = n
    while nactive > G:
        best, bp, bq = np.inf, -1, -1
        for p in range(n):
            if not active[p]:
                continue
            for q in range(p + 1, n):
                if not active[q]:
                    continue
                d = max(dist[x, y] for x in clusters[p] for y in clusters[q])
                if d < best - GlobalConstants.Zero:
                    best, bp, bq = d, p, q
        if bp < 0:
            break
        clusters[bp] = sorted(clusters[bp] + clusters[bq])
        clusters[bq] = []
        active[bq] = False
        nactive -= 1

    # Relabel by smallest member, so the group numbering is canonical
    live = [p for p in range(n) if active[p]]
    live.sort(key=lambda p: clusters[p][0])
    for g, p in enumerate(live):
        for x in clusters[p]:
            grp[x] = g
    return grp


def bgchain_env(bg, i):
    """Lump the background chain onto the closed occupancy of station i.

    Station i does not observe the whole closed population vector, only how many
    closed jobs compete with the open ones for its server. The lumped generator
    is the stationary-weighted aggregation of the chain's generator over the
    level sets {s : totocc(s,i) = e}, which is exact when the partition is
    lumpable in the Kemeny-Snell sense and is the standard exact-aggregation
    approximation otherwise. The diagonal is set from the off-diagonal row sums,
    so the result is a proper generator whatever the lumping error is.

    Environment states of zero stationary probability are unreachable and are
    dropped, so the support need not be 0..N.

    Returns:
        Tuple ``(A, phi, esup)``.
    """
    occ = np.asarray(bg['totocc'])[:, i]
    levels = np.unique(occ)
    w = np.array([bg['pi'][occ == e].sum() for e in levels])
    keep = w > GlobalConstants.Zero
    if not np.any(keep):
        # degenerate chain: the station never holds a closed job
        return np.zeros((1, 1)), np.array([1.0]), np.array([0])
    esup = levels[keep].astype(int)
    w = w[keep]
    me = len(esup)

    A = np.zeros((me, me))
    if me > 1:
        lvl = -np.ones(len(occ), dtype=int)
        for e in range(me):
            lvl[occ == esup[e]] = e
        rows, cols = np.nonzero(bg['Q'])
        for s, sp_ in zip(rows, cols):
            if s == sp_:
                continue
            e, ep = lvl[s], lvl[sp_]
            if e < 0 or ep < 0 or e == ep:
                continue
            A[e, ep] += bg['pi'][s] * bg['Q'][s, sp_]
        for e in range(me):
            if w[e] > 0:
                A[e, :] /= w[e]
        np.fill_diagonal(A, 0.0)
        A = A - np.diag(A.sum(axis=1))
    phi = w / w.sum()
    return A, phi, esup


def _open_share(k, esup, nservers):
    """Share of the server capacity that k open jobs hold when esup closed jobs
    are also present: min(k+e,c) busy servers times the open fraction k/(k+e)."""
    tot = k + np.asarray(esup, dtype=float)
    v = np.zeros(len(tot))
    nz = tot > 0
    v[nz] = np.minimum(tot[nz], nservers) * (k / tot[nz])
    return v


def _closed_share(esup, k, nservers):
    """The mirror image of ``_open_share``, for the closed jobs."""
    esup = np.asarray(esup, dtype=float)
    tot = esup + k
    v = np.zeros(len(tot))
    nz = tot > 0
    v[nz] = np.minimum(tot[nz], nservers) * (esup[nz] / tot[nz])
    return v


def _env_at_level(Aup, Adown, esup, gref, k, nservers):
    """Environment generator seen at open level k: the closed departures are
    rescaled from the averaged share gref to the share they hold against k open
    jobs, the closed arrivals are unchanged, the diagonal is rebuilt."""
    g = _closed_share(esup, k, nservers)
    ratio = np.ones(len(g))
    nz = np.asarray(gref) > 0
    ratio[nz] = g[nz] / np.asarray(gref)[nz]
    Ak = Aup + np.diag(ratio) @ Adown
    return Ak - np.diag(Ak.sum(axis=1))


def bgchain_station(Da0, Da1, alpha_s, T, A, esup, nservers, gref, Kmax, tol, iter_max):
    """Solve the open classes of one station as a modulated level-dependent QBD.

    Level = number of open jobs held by the station, phase = (arrival MAP phase,
    environment state, service phase). With k open and e closed jobs present the
    open aggregate completes at rate ``phi(k,e) = min(k+e,c) * k/(k+e)`` times
    the phase-type completion rate of one busy server. The dependence on k is
    what makes the QBD level-dependent, the dependence on e is what makes it
    modulated. The c parallel servers are collapsed into a single phase-type
    process scaled by phi, exact for exponential service at any c and for
    phase-type service at c = 1.

    The environment is level-dependent too: a lumped transition that LOWERS the
    closed occupancy is a closed completion here, so level k rescales it by the
    ratio of shares, while a transition that RAISES it is an arrival from
    elsewhere and is left alone.

    The level space is truncated at Kmax. An arrival at the top level is lost
    but still advances the arrival phase, so the arrival process keeps its exact
    marginal and autocorrelation and only the queue tail is cut.

    Returns:
        dict with QLen, Util, Tput, ploss, penv, cshare, esup
    """
    Da0 = np.atleast_2d(np.asarray(Da0, dtype=float))
    Da1 = np.atleast_2d(np.asarray(Da1, dtype=float))
    T = np.atleast_2d(np.asarray(T, dtype=float))
    alpha_s = np.asarray(alpha_s, dtype=float).reshape(1, -1)
    esup = np.asarray(esup, dtype=int)
    gref = np.asarray(gref, dtype=float)

    order = np.argsort(esup)
    esup = esup[order]
    gref = gref[order]
    A = np.atleast_2d(np.asarray(A, dtype=float))[np.ix_(order, order)]

    ma = Da0.shape[0]
    me = len(esup)
    ms = alpha_s.shape[1]
    t = -T @ np.ones((ms, 1))
    Ima, Ime, Ims = np.eye(ma), np.eye(me), np.eye(ms)
    Kmax = max(1, int(round(Kmax)))

    Adown = np.tril(A, -1)
    Aup = np.triu(A, 1)

    Q0, Q1, Q2, phiae = [], [], [], []
    A0 = _env_at_level(Aup, Adown, esup, gref, 0, nservers)
    Q1.append(np.kron(Da0, Ime) + np.kron(Ima, A0))
    Q0.append(np.kron(np.kron(Da1, Ime), alpha_s))

    Da0kron = np.kron(np.kron(Da0, Ime), Ims)
    Da1kron = np.kron(np.kron(Da1, Ime), Ims)
    for k in range(1, Kmax + 1):
        rep = np.tile(_open_share(k, esup, nservers), ma)
        phiae.append(rep)
        Ak = _env_at_level(Aup, Adown, esup, gref, k, nservers)
        Q1.append(Da0kron + np.kron(np.kron(Ima, Ak), Ims) + np.kron(np.diag(rep), T))
        if k < Kmax:
            Q0.append(Da1kron)
        if k == 1:
            Q2.append(np.kron(np.diag(rep), t))
        else:
            Q2.append(np.kron(np.diag(rep), t @ alpha_s))
    # truncation: an arrival at the top level is lost, its phase transition is kept
    Q1[Kmax] = Q1[Kmax] + Da1kron

    res = ldqbd(Q0, Q1, Q2, LdqbdOptions(epsilon=tol, max_iter=iter_max, verbose=False))
    plev = np.maximum(np.asarray(res.pi, dtype=float).ravel(), 0.0)
    if plev.sum() > 0:
        plev = plev / plev.sum()
    pcell = res.pi_cells

    QLen = float(np.arange(Kmax + 1) @ plev)
    util = 0.0
    tput = 0.0
    penv = np.zeros(me)
    gacc = np.zeros(me)
    for k in range(Kmax + 1):
        pk = np.maximum(np.asarray(pcell[k], dtype=float).ravel(), 0.0)
        if k == 0:
            marg = pk.reshape(ma, me).sum(axis=0)
        else:
            blocks = pk.reshape(ma * me, ms)
            marg = blocks.sum(axis=1).reshape(ma, me).sum(axis=0)
            util += float(blocks.sum(axis=1) @ phiae[k - 1])
            tput += float((blocks @ t).ravel() @ phiae[k - 1])
        penv += marg
        gacc += marg * _closed_share(esup, k, nservers)
    psum = penv.sum()
    if psum > 0:
        penv = penv / psum
        gacc = gacc / psum
    cshare = np.zeros(me)
    nz = penv > GlobalConstants.Zero
    cshare[nz] = gacc[nz] / penv[nz]

    return {'QLen': QLen, 'Util': util / nservers, 'Tput': tput,
            'ploss': float(plev[Kmax]), 'penv': penv, 'cshare': cshare, 'esup': esup}


class BgchainAlgorithm(MAMAlgorithm):
    """Closed classes as a background modulating chain, open classes as QBDs."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        # The closed classes ARE the background chain, so a purely open model has
        # nothing to build it from. A purely CLOSED one is accepted: it is the
        # degenerate case where the chain answers alone, with no open work to
        # modulate it.
        if not np.any(np.isfinite(njobs)):
            return False, ('The bgchain method requires at least one closed class: the background '
                           'chain IS the closed population vector, which a purely open model does '
                           'not have. Use the dec.source method.')
        # Priority disciplines need the per-class QBD of MMAPPH1PRPR, which has
        # no counterpart in the modulated level-dependent QBD this method builds.
        prio_sched = any(_sched_name(sn, i) in ('HOL', 'FCFSPRPRIO') for i in range(sn.nstations))
        classprio = np.asarray(getattr(sn, 'classprio', np.zeros(sn.nclasses)), dtype=float).ravel()
        if prio_sched and classprio.size > 1 and np.any(classprio != classprio[0]):
            return False, ('The bgchain method does not support class priorities: it aggregates '
                           'the open classes into one phase-type mixture per station, which cannot '
                           'express a priority order. Use the dec.source method.')
        from ....api.sn.predicates import sn_has_fork_join
        if sn_has_fork_join(sn):
            return False, ('The bgchain method does not support fork-join: the background chain '
                           'conserves the closed population per station, which a fork violates. '
                           'Use the dec.source method.')
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        start_time = time.time()

        ok, reason = BgchainAlgorithm.supports_network(sn)
        if not ok:
            raise ValueError(reason)

        tol = float(getattr(options, 'tol', 1e-6) or 1e-6)
        iter_max = int(getattr(options, 'iter_max', 100) or 100)
        config = getattr(options, 'config', None) or {}
        if not isinstance(config, dict):
            config = dict(getattr(config, '__dict__', {}))
        space_max = int(config.get('space_max', 128))
        qbdphases_max = int(config.get('qbdphases_max', 500))
        bgstates_max = int(config.get('bgstates_max', 20000))
        bgaggr_opt = int(config.get('bgaggr', 1))
        cutoff_opt = getattr(options, 'cutoff', None)

        M, K, C = sn.nstations, sn.nclasses, sn.nchains
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        rates = np.asarray(sn.rates, dtype=float)
        nservers = np.asarray(sn.nservers, dtype=float).ravel()
        with np.errstate(divide='ignore', invalid='ignore'):
            S = 1.0 / rates
        S[~np.isfinite(S)] = 0.0

        rtst, V = sn_rt_stations(sn)
        dem = sn_get_demands_chain(sn)
        Vchain = np.asarray(dem.Vchain, dtype=float)
        STchain = np.asarray(dem.STchain, dtype=float)
        Lchain = np.asarray(dem.Lchain, dtype=float)
        alpha = np.asarray(dem.alpha, dtype=float)
        Nchain = np.asarray(dem.Nchain, dtype=float).ravel()
        refstat = np.asarray(sn.refstat, dtype=int).ravel()

        inchain = {c: np.asarray(sn.inchain[c], dtype=int).ravel() for c in range(C)}
        isopenchain = np.array([bool(np.any(np.isinf(njobs[inchain[c]]))) for c in range(C)])
        open_chains = [c for c in range(C) if isopenchain[c]]
        closed_chains = [c for c in range(C) if not isopenchain[c] and Nchain[c] > 0]
        R = len(closed_chains)
        if R == 0:
            raise ValueError('The bgchain method requires at least one closed class: the background '
                             'chain IS the closed population vector, so a purely open model has '
                             'nothing to build it from. Use dec.source.')

        # Stations the closed chains visit: the support of the background chain
        inbg = np.zeros(M, dtype=bool)
        for c in closed_chains:
            inbg |= Vchain[:, c] > GlobalConstants.Zero
        cst = np.where(inbg)[0]
        Mc = len(cst)
        if Mc == 0:
            raise ValueError('The closed classes of this model visit no station.')

        # Chain-level station routing, folding the class axis of sn.rt
        Pchain = []
        for c in range(C):
            P = np.zeros((M, M))
            for i in range(M):
                for k in inchain[c]:
                    a = alpha[i, k]
                    if a <= 0:
                        continue
                    row = rtst[i * K + k, :]
                    acc = np.zeros(M)
                    for kp in inchain[c]:
                        acc += row[np.arange(M) * K + kp]
                    P[i, :] += a * acc
            Pchain.append(P)

        # Open arrival streams
        lambda_chain = np.zeros(C)
        chain_arrival = {}
        for c in open_chains:
            isrc = int(refstat[inchain[c][0]])
            acc = None
            lam = 0.0
            for k in inchain[c]:
                rk = rates[isrc, k]
                if not np.isfinite(rk) or rk <= 0:
                    continue
                lam += rk
                entry = sn.proc[isrc][k] if isinstance(sn.proc, dict) else sn.proc[isrc][k]
                d0d1 = _proc_to_d0d1(entry)
                if d0d1 is None or not np.all(np.isfinite(d0d1[0])):
                    continue
                mk = (d0d1[0], [d0d1[1]])
                acc = mk if acc is None else mmap_super_safe([acc, mk], space_max, 'default')
            lambda_chain[c] = lam
            if acc is not None:
                D_list = acc[1]
                D1 = D_list[0] if len(D_list) == 1 else sum(D_list)
                chain_arrival[c] = (np.atleast_2d(acc[0]), np.atleast_2d(D1))

        lambda_open = np.zeros((M, K))
        isopenclass = np.zeros(K, dtype=bool)
        for c in open_chains:
            for k in inchain[c]:
                isopenclass[k] = True
                lambda_open[:, k] = lambda_chain[c] * V[:, k]

        QN = np.zeros((M, K)); UN = np.zeros((M, K)); RN = np.zeros((M, K))
        TN = np.zeros((M, K)); CN = np.zeros(K); XN = np.zeros(K)

        # cshare[i,e]: mean number of servers of station i that its e closed jobs
        # hold once the open work has taken its share. Starts at min(e,c).
        Ntot = int(sum(Nchain[c] for c in closed_chains))
        egrid = np.arange(Ntot + 1)
        cshare = np.zeros((M, Ntot + 1))
        for i in range(M):
            cshare[i, :] = np.minimum(egrid, nservers[i])
        Xclosed = np.zeros(C)
        for c in closed_chains:
            denom = Lchain[:, c].sum()
            if denom > 0:
                Xclosed[c] = Nchain[c] / denom

        # How many aggregate classes the background chain carries, and which
        # chains share each of them. config['bgaggr'] is the number of AGGREGATE
        # classes G: G = 1 is the classic tagged/aggregate pair, G >= R-1
        # aggregates nothing.
        naggr = min(max(bgaggr_opt, 1), max(R - 1, 1))
        # With nothing left to aggregate ONE background chain carries every
        # closed chain exactly, so the tagged loop would repeat it R times.
        no_aggr = (R == 1) or (naggr >= R - 1)

        # The grouping is a property of the demands, not of the iterate, so it is
        # fixed once here rather than recomputed inside the fixed point.
        others_of, grp_of = {}, {}
        if not no_aggr:
            for r in closed_chains:
                others = [c for c in closed_chains if c != r]
                others_of[r] = others
                grp_of[r] = bgchain_groups(Lchain[np.ix_(cst, others)], naggr)
        passes = [closed_chains[0]] if no_aggr else closed_chains
        npass = len(passes)

        TN_prev = np.full((M, K), np.inf)
        totiter = 0
        relax = 0.5

        while np.max(np.abs(TN - TN_prev)) > tol and totiter < iter_max:
            totiter += 1
            TN_prev = TN.copy()

            Qopen_acc = np.zeros(M)
            Uopen_acc = np.zeros(M)
            cshare_acc = np.zeros((M, Ntot + 1))

            for r in passes:
                # Background classes: class 0 is the tagged chain, classes 1..G
                # the flow-equivalent aggregates of the demand-similar groups.
                # With no aggregation every closed chain is a class of its own.
                if no_aggr:
                    members = [[c] for c in closed_chains]
                else:
                    others = others_of[r]
                    grp = grp_of[r]
                    members = [[r]]
                    for g in range(naggr):
                        members.append([o for oi, o in enumerate(others) if grp[oi] == g])
                B = len(members)
                Nb = [0] * B
                STb = np.zeros((Mc, B))
                Pb = []
                # A class can only hold jobs at the stations its members visit;
                # see bgchain_ctmc on why the union makes the chain reducible.
                suppb = np.zeros((Mc, B), dtype=bool)
                for b, mem in enumerate(members):
                    Nb[b] = int(sum(Nchain[c] for c in mem))
                    for o in mem:
                        suppb[:, b] |= Vchain[cst, o] > GlobalConstants.Zero
                    if len(mem) == 1:
                        # a group of one is carried exactly: no mean to take
                        STb[:, b] = STchain[cst, mem[0]]
                        Pb.append(Pchain[mem[0]][np.ix_(cst, cst)].copy())
                    elif len(mem) == 0:
                        Pb.append(np.zeros((Mc, Mc)))
                    else:
                        w = np.zeros((Mc, len(mem)))
                        for oi, o in enumerate(mem):
                            w[:, oi] = Xclosed[o] * Vchain[cst, o]
                        rowsum = w.sum(axis=1)
                        for ii in range(Mc):
                            if rowsum[ii] > 0:
                                w[ii, :] /= rowsum[ii]
                            else:
                                w[ii, :] = 1.0 / len(mem)
                        STb[:, b] = (w * STchain[np.ix_(cst, mem)]).sum(axis=1)
                        Pagg = np.zeros((Mc, Mc))
                        for oi, o in enumerate(mem):
                            Pagg += w[:, oi][:, None] * Pchain[o][np.ix_(cst, cst)]
                        Pb.append(Pagg)
                Pb = [_row_normalize(P) for P in Pb]

                sched_names = [_sched_name(sn, int(i)) for i in cst]
                bg = bgchain_ctmc(Nb, STb, Pb, sched_names, nservers[cst],
                                  cshare[cst, :], suppb, bgstates_max)

                # Closed-class metrics of every chain this pass carries EXACTLY:
                # the tagged one always, and every chain when nothing was aggregated.
                bext = range(B) if no_aggr else [0]
                for b in bext:
                    if not members[b]:
                        continue
                    rb = members[b][0]
                    for k in inchain[rb]:
                        QN[:, k] = 0.0; UN[:, k] = 0.0; RN[:, k] = 0.0; TN[:, k] = 0.0
                    for ii, i in enumerate(cst):
                        for k in inchain[rb]:
                            a = alpha[i, k]
                            if a <= 0:
                                continue
                            # THROUGHPUT splits by VISIT share, OCCUPANCY by
                            # DEMAND share. A chain queue divided by alpha alone
                            # gives every class of a station the same response
                            # time, impossible at a Delay where R must be the
                            # class service time; the weight is alpha*ST/STchain,
                            # the rule sn_deaggregate_chain_results applies.
                            stc = STchain[i, rb]
                            wk = a * S[i, k] / stc if stc > GlobalConstants.Zero else a
                            q = bg['QLen'][ii, b] * wk
                            x = bg['Tput'][ii, b] * a
                            QN[i, k] = q
                            TN[i, k] = x
                            UN[i, k] = q if sched_names[ii] == 'INF' else bg['Ubusy'][ii, b] * wk
                            RN[i, k] = q / x if x > GlobalConstants.Zero else 0.0
                    iref = int(refstat[inchain[rb][0]])
                    tputref = float(TN[iref, inchain[rb]].sum())
                    vref = Vchain[iref, rb]
                    Xclosed[rb] = tputref / vref if vref > GlobalConstants.Zero else tputref

                Uclosed = np.zeros(M)
                Uclosed[cst] = bg['Ubusy'].sum(axis=1)

                self._open_pass(sn, bg, cst, S, V, lambda_open, lambda_chain, chain_arrival,
                                open_chains, inchain, Uclosed, cshare, nservers, cutoff_opt,
                                space_max, qbdphases_max, tol, iter_max,
                                Qopen_acc, Uopen_acc, cshare_acc)

            Qopen = Qopen_acc / npass
            cshare = (1 - relax) * cshare + relax * (cshare_acc / npass)

            # Open-class metrics from the aggregate station results
            for i in range(M):
                kopen = [k for k in range(K)
                         if isopenclass[k] and lambda_open[i, k] > GlobalConstants.Zero]
                sched = _sched_name(sn, i)
                if not kopen:
                    for k in range(K):
                        if isopenclass[k]:
                            TN[i, k] = lambda_open[i, k]
                            QN[i, k] = 0.0; UN[i, k] = 0.0; RN[i, k] = 0.0
                    continue
                lam = np.array([lambda_open[i, k] for k in kopen])
                lamtot = lam.sum()
                Smix = float(sum(lambda_open[i, k] * S[i, k] for k in kopen) / lamtot)
                for k in kopen:
                    TN[i, k] = lambda_open[i, k]
                    if sched == 'EXT':
                        QN[i, k] = 0.0; UN[i, k] = 0.0; RN[i, k] = 0.0
                    elif sched == 'INF':
                        RN[i, k] = S[i, k]
                        QN[i, k] = lambda_open[i, k] * S[i, k]
                        UN[i, k] = QN[i, k]
                    else:
                        Rtot = Qopen[i] / lamtot
                        if sched == 'PS':
                            # processor sharing: residence scales with the demand
                            rk = Rtot * S[i, k] / Smix
                        else:
                            # FCFS and its variants: the wait is class-blind, the
                            # service time is not
                            rk = max(S[i, k], Rtot - Smix + S[i, k])
                        RN[i, k] = rk
                        QN[i, k] = lambda_open[i, k] * rk
                        # Utilization Law: a c-server station holds TN*S/c
                        UN[i, k] = lambda_open[i, k] * S[i, k] / nservers[i]

        for c in range(C):
            for k in inchain[c]:
                XN[k] = lambda_chain[c] if isopenchain[c] else Xclosed[c]
        CN = RN.sum(axis=0)

        for arr in (QN, UN, RN, TN, CN, XN):
            arr[~np.isfinite(arr)] = 0.0

        return MAMResult(QN=QN, UN=UN, RN=RN, TN=TN, CN=CN.reshape(1, -1),
                         XN=XN.reshape(1, -1), totiter=totiter, method='bgchain',
                         runtime=time.time() - start_time)

    def _open_pass(self, sn, bg, cst, S, V, lambda_open, lambda_chain, chain_arrival,
                   open_chains, inchain, Uclosed, cshare, nservers, cutoff_opt,
                   space_max, qbdphases_max, tol, iter_max,
                   Qopen_acc, Uopen_acc, cshare_acc):
        """One pass of the open side: for each station, lump the background chain
        onto its closed occupancy and solve the resulting modulated QBD."""
        M, K = sn.nstations, sn.nclasses
        ngrid = cshare.shape[1]
        egrid = np.arange(ngrid)

        for i in range(M):
            solved = False
            sched = _sched_name(sn, i)
            if sched not in ('EXT', 'INF'):
                kopen = [k for k in range(K) if lambda_open[i, k] > GlobalConstants.Zero]
                if kopen:
                    Da = None
                    for c in open_chains:
                        if lambda_chain[c] <= GlobalConstants.Zero or c not in chain_arrival:
                            continue
                        rate_ic = lambda_chain[c] * float(V[i, inchain[c]].sum())
                        if rate_ic <= GlobalConstants.Zero:
                            continue
                        base = chain_arrival[c]
                        scaled = map_scale(base[0], base[1], 1.0 / rate_ic)
                        if Da is None:
                            Da = (np.atleast_2d(scaled[0]), np.atleast_2d(scaled[1]))
                        else:
                            sup = mmap_super_safe([(Da[0], [Da[1]]),
                                                   (scaled[0], [scaled[1]])],
                                                  space_max, 'default')
                            D_list = sup[1]
                            D1 = D_list[0] if len(D_list) == 1 else sum(D_list)
                            Da = (np.atleast_2d(sup[0]), np.atleast_2d(D1))
                    if Da is not None:
                        # Arrival-weighted phase-type mixture of the open service laws
                        lamtot = float(sum(lambda_open[i, k] for k in kopen))
                        pies, subgens = [], []
                        # PROCESSOR SHARING IS INSENSITIVE to the service law beyond its mean, so
                        # carrying the phase-type representation at a PS station is not merely unnecessary,
                        # it is WRONG. This QBD tracks ONE service phase for the whole station, which makes
                        # the open queue length inherit the SCV-sensitivity of an M/PH/1 FCFS queue;
                        # measured against SolverCTMC, a HyperExp of SCV 4 then read 21% high where the
                        # exact answer is the exponential one to five digits. The exponential of the same
                        # mean is exact here, and it shrinks the QBD's phase count as a side effect.
                        is_ps = (sched == 'PS')
                        for k in kopen:
                            entry = sn.proc[i][k]
                            d0d1 = _proc_to_d0d1(entry)
                            if is_ps or d0d1 is None or not np.all(np.isfinite(d0d1[0])):
                                rate = 1.0 / S[i, k] if S[i, k] > 0 else 1.0 / GlobalConstants.FineTol
                                d0, d1 = np.array([[-rate]]), np.array([[rate]])
                            else:
                                d0, d1 = map_scale(d0d1[0], d0d1[1], S[i, k])
                            pies.append(np.asarray(map_pie(d0, d1), dtype=float).ravel()
                                        * (lambda_open[i, k] / lamtot))
                            subgens.append(np.atleast_2d(d0))
                        ms_total = sum(g.shape[0] for g in subgens)
                        alpha_s = np.zeros(ms_total)
                        Tblk = np.zeros((ms_total, ms_total))
                        off = 0
                        for pik, g in zip(pies, subgens):
                            n = g.shape[0]
                            alpha_s[off:off + n] = pik
                            Tblk[off:off + n, off:off + n] = g
                            off += n

                        pos = np.where(cst == i)[0]
                        if len(pos) > 0:
                            A, _, esup = bgchain_env(bg, int(pos[0]))
                        else:
                            A, esup = np.zeros((1, 1)), np.array([0])

                        nphases = Da[0].shape[0] * len(esup) * ms_total
                        if nphases > qbdphases_max:
                            raise RuntimeError(
                                'The modulated QBD of station %d needs %d phases (%d arrival x %d '
                                'environment x %d service), above the limit of %d. The environment '
                                'axis is the closed population held by the station, so it grows '
                                'with the closed population. Raise options.config["qbdphases_max"], '
                                'or reduce the closed population or the order of the arrival and '
                                'service processes.'
                                % (i, nphases, Da[0].shape[0], len(esup), ms_total, qbdphases_max))

                        gref = cshare[i, np.minimum(esup, ngrid - 1)]
                        svc_work = float(sum(lambda_open[i, k] * S[i, k] for k in kopen))
                        Kmax = _cutoff(cutoff_opt, lamtot, svc_work / lamtot,
                                       nservers[i], Uclosed[i])
                        res = bgchain_station(Da[0], Da[1], alpha_s, Tblk, A, esup,
                                              nservers[i], gref, Kmax, tol, iter_max)

                        Qopen_acc[i] += res['QLen']
                        Uopen_acc[i] += res['Util']
                        # The QBD only saw the environment states the background
                        # chain reaches; interpolate the rest so the next chain has
                        # a share wherever it may go, clipped to what a closed job
                        # can physically hold.
                        if len(res['esup']) > 1:
                            g = np.interp(egrid, res['esup'], res['cshare'])
                            lo = res['esup'][0]
                            hi = res['esup'][-1]
                            if hi > lo:
                                slope_lo = (res['cshare'][1] - res['cshare'][0]) / (res['esup'][1] - res['esup'][0])
                                slope_hi = (res['cshare'][-1] - res['cshare'][-2]) / (res['esup'][-1] - res['esup'][-2])
                                g = np.where(egrid < lo, res['cshare'][0] + slope_lo * (egrid - lo), g)
                                g = np.where(egrid > hi, res['cshare'][-1] + slope_hi * (egrid - hi), g)
                        else:
                            g = np.full(ngrid, res['cshare'][0])
                        cshare_acc[i, :] += np.clip(g, 0.0, np.minimum(egrid, nservers[i]))
                        solved = True
            if not solved:
                # a station with no open queue keeps the share it had
                cshare_acc[i, :] += cshare[i, :]


def _row_normalize(P):
    """Row-normalize a routing matrix, leaving an all-zero row as a self-loop so
    the background chain stays a proper Markov chain on its support."""
    P = np.asarray(P, dtype=float).copy()
    for i in range(P.shape[0]):
        s = P[i, :].sum()
        if s > GlobalConstants.Zero:
            P[i, :] /= s
        else:
            P[i, :] = 0.0
            P[i, i] = 1.0
    return P


def _cutoff(cutoff_opt, lam, Smix, nservers, Uclosed):
    """Truncation level of the open queue: the explicit cutoff when the user set
    one, else enough levels for the geometric tail left by the closed traffic to
    be negligible."""
    if cutoff_opt is not None:
        cv = np.asarray(cutoff_opt, dtype=float).ravel()
        if cv.size > 0 and np.isfinite(cv[0]) and cv[0] > 0:
            return max(2, int(round(float(cv[0]))))
    free = max(GlobalConstants.FineTol, 1.0 - Uclosed)
    rho = lam * Smix / (nservers * free)
    rho = min(max(rho, 1e-3), 1 - 1e-3)
    Kmax = int(np.ceil(np.log(1e-8) / np.log(rho)))
    return min(max(Kmax, 20), 200)
