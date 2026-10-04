"""
Miner-Ciardo-Donatelli approximate stationary analysis.

Solve a structured CTMC whose EXACT reachable state space is stored in a
decision diagram, by building and iterating K level-CTMCs (a
decision-diagram-guided aggregation), after A.S. Miner, G. Ciardo,
S. Donatelli, "Using the exact state space of a Markov model to compute
approximate stationary measures", SIGMETRICS 2000.

The method never forms the |S|-state generator or probability vector. It keeps
one CTMC per decision-diagram level k, over states M_k = {(p,i_k)} with p a
level-k node and i_k a local state on a non-null arc, and iterates the coupled
system to a fixed point. The single approximation (Eq. 5) is
Pr{i_k | alpha} = Pr{i_k | p}: the local-state law at level k depends only on
the node p, not the full path above it, which the exact reachability the
diagram encodes justifies. For product-form models the method is EXACT (paper
Sec. 5), so on a single-class closed QN it reproduces SolverCTMC.

Orientation note: the paper indexes levels K (top/root) down to 1
(bottom/terminal); the MDD class uses level 1 as the root. This module works in
the paper's orientation, mapping paper level k to MDD level (K+1-k), i.e. to
station (K+1-k). All indices here are 0-based, so paper level k (0-based, 0 =
bottom) maps to MDD level K-1-k and to station K-1-k.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Optional

import numpy as np

try:
    import scipy.sparse as sp
except ImportError:                      # descriptors may still hand over dense W
    sp = None

from ..io.logging import line_error, line_printf
from .mdd import TERM_TRUE


def mdd_mcd(mdds, desc, options: Optional[Dict] = None) -> Dict[str, object]:
    """Approximate stationary measures by decision-diagram-guided aggregation.

    Parameters
    ----------
    mdds : MDDStruct from MDD.to_struct (the reachable set, MDD orientation)
    desc : Kronecker rate descriptor from mdd_descriptor / mdd_ps / spn_mdd
    options : dict with optional keys tol (1e-12), maxiter (500),
        verbose (False), initpik (list of K level warm-start vectors)

    Returns
    -------
    dict with keys QLen, X, U, pik, Mrows, levelSizes, iters.
    """
    if options is None:
        options = {}
    tol = float(options.get('tol', 1e-12))
    maxiter = int(options.get('maxiter', 500))
    verbose = bool(options.get('verbose', False))

    K = int(mdds.K)

    # ---- paper orientation: paper level k <-> MDD level K-1-k = station K-1-k
    Pnode = [None] * K
    nn = np.zeros(K, dtype=np.int64)
    dom = np.zeros(K, dtype=np.int64)
    for k in range(K):
        oL = K - 1 - k
        Pnode[k] = np.asarray(mdds.node[oL], dtype=np.int64).reshape(int(mdds.nnodes[oL]),
                                                                    int(mdds.domain[oL]))
        nn[k] = int(mdds.nnodes[oL])
        dom[k] = int(mdds.domain[oL])

    # ---- per (event, paper level) local matrices W_k^e and enabling rates
    events = desc['events']
    E = len(events)
    # wrows[e][k][v] = list of (column, value) of row v of W_k^e
    wrows: List[List[List[List]]] = [[None] * K for _ in range(E)]
    lam: List[List[np.ndarray]] = [[None] * K for _ in range(E)]
    for e in range(E):
        ev = events[e]
        levs = list(ev['lev'])                       # 0-based station levels
        plev = [K - 1 - lv for lv in levs]           # -> paper levels
        for k in range(K):
            if k in plev:
                rows, lam[e][k] = _row_nonzeros(ev['W'][plev.index(k)], int(dom[k]))
            else:
                # untouched level: identity, i.e. one unit self-arc per value
                rows = [[(i, 1.0)] for i in range(int(dom[k]))]
                lam[e][k] = np.ones(int(dom[k]))
            wrows[e][k] = rows

    # ---- level-k CTMC state sets M_k = {(p, v) : arc p[v] non-null}
    Mrows: List[np.ndarray] = [None] * K
    Midx: List[np.ndarray] = [None] * K
    for k in range(K):
        tbl = Pnode[k]
        if k == 0:
            mask = tbl == TERM_TRUE          # bottom level: TRUE arcs
        else:
            mask = tbl > 0                   # real child ids
        # column-major order, matching the MATLAB find() the reference uses
        cc, pp = np.nonzero(mask.T)
        rows = np.column_stack([pp + 1, cc]).astype(np.int64)
        Mrows[k] = rows
        idx = np.zeros((int(nn[k]), int(dom[k])), dtype=np.int64)
        idx[rows[:, 0] - 1, rows[:, 1]] = np.arange(1, rows.shape[0] + 1)
        Midx[k] = idx
    levelSizes = np.array([Mrows[k].shape[0] for k in range(K)], dtype=np.int64)

    # ---- initialise level stationary vectors and node marginals
    above, below = _path_counts(mdds, K)
    initpik = options.get('initpik')
    if initpik:
        pik = [np.asarray(initpik[k], dtype=float).ravel() for k in range(K)]
    else:
        pik = _uniform_init(mdds, Mrows, K, above, below)
    Prp = [np.bincount(Mrows[k][:, 0] - 1, weights=pik[k], minlength=int(nn[k]))
           for k in range(K)]

    # ---- fixed-point iteration (Fig. 3, procedure Solve)
    iters = 0
    converged = False
    delta = np.inf
    for it in range(1, maxiter + 1):
        iters = it
        piold = [p.copy() for p in pik]

        # ComputeBs, bottom-up: b_k^e[p] = sum_v Pr{v|p} b_{k-1}^e[p[v]] lambda
        bcell: List[np.ndarray] = [None] * K
        for k in range(K):
            bk = np.zeros((int(nn[k]), E))
            rows = Mrows[k]
            for r in range(rows.shape[0]):
                p = int(rows[r, 0])
                v = int(rows[r, 1])
                if Prp[k][p - 1] <= 0:
                    continue
                adjust = pik[k][r] / Prp[k][p - 1]           # Pr{v|p}
                for e in range(E):
                    le = lam[e][k][v]
                    if le == 0:
                        continue                              # not locally enabled
                    if k > 0:
                        down = bcell[k - 1][int(Pnode[k][p - 1, v]) - 1, e]
                    else:
                        down = 1.0                            # terminal ONE
                    bk[p - 1, e] += adjust * down * le
            bcell[k] = bk

        # top-down: ComputeAs(k) then SolveLevel(k)
        Acell: List[List[np.ndarray]] = [None] * K
        Acell[K - 1] = [np.eye(int(nn[K - 1])) for _ in range(E)]
        for k in range(K - 1, -1, -1):
            if k < K - 1:
                Acell[k] = _compute_as(k, Acell[k + 1], Pnode, pik, wrows, Mrows, nn, E)
            # SolveLevel(k): assemble R_k (Eq. 6), solve pi_k Q_k = 0
            Rk = _compute_mc(k, Acell[k], bcell, Pnode, wrows, Mrows, Midx, levelSizes, E)
            Qk = Rk - np.diag(Rk.sum(axis=1))
            pik[k] = _solve_stat(Qk)
            Prp[k] = np.bincount(Mrows[k][:, 0] - 1, weights=pik[k], minlength=int(nn[k]))

        # NaN-aware: a diverged iterate must not be read as converged, which is
        # what a max() that skips NaN would do.
        delta = 0.0
        for k in range(K):
            dk = float(np.max(np.abs(pik[k] - piold[k]))) if pik[k].size else 0.0
            if not np.isfinite(dk):
                line_error('mdd_mcd',
                           'level %d iterate is not finite at iteration %d; the level-%d CTMC '
                           'did not yield a proper stationary vector.' % (k + 1, it, k + 1))
            delta = max(delta, dk)
        if delta < tol:
            converged = True
            break
    if not converged:
        line_error('mdd_mcd',
                   'the coupled level iteration did not converge in %d sweeps (last change %.3e '
                   'against tol %.3e); the level marginals returned would not be a fixed point. '
                   'Raise maxiter or relax tol.' % (maxiter, delta, tol))

    # ---- performance measures from the per-level marginals. QLen is the mean
    # local value and is defined for any descriptor (jobs at a station, tokens
    # in a place); X and U need the queueing parameters and are skipped without
    # them.
    mu = desc.get('mu')
    servers = desc.get('servers')
    isQN = mu is not None and len(np.ravel(mu)) > 0 and servers is not None
    valuemap = desc.get('valuemap')
    QLen = np.zeros(K)
    X = np.zeros(K) if isQN else None
    U = np.zeros(K) if isQN else None
    for s in range(K):
        k = K - 1 - s                        # paper level of station/place s
        # A level whose local state encodes more than a count (a station holding
        # both a population and a service phase) carries a map from local index
        # to the physical quantity; without one the index IS the quantity.
        if valuemap is not None:
            vmap = np.asarray(valuemap[s], dtype=float).ravel()
            v = vmap[Mrows[k][:, 1]]
        else:
            v = Mrows[k][:, 1].astype(float)
        pk = pik[k]
        QLen[s] = float(v @ pk)              # E[occupancy of level s]
        if isQN:
            busy = np.minimum(v, float(servers[s]))
            X[s] = float(mu[s]) * float(busy @ pk)
            if np.isinf(servers[s]):
                U[s] = float(v @ pk)
            else:
                U[s] = float(busy @ pk) / float(servers[s])

    # The level chains are coupled only through rates, so nothing in the
    # iteration forces the marginals to describe the same population; a fixed
    # point that does not is a wrong answer, not an approximation, and must not
    # be returned. The test is a conservation law of the model: the closed
    # population for a QN, a place invariant w'*m = const for a Petri net.
    winv, vinv = _invariant(desc, K)
    if winv is not None:
        got = float(np.asarray(winv, dtype=float).ravel() @ QLen)
        if abs(got - vinv) > 1e-6 * max(1.0, abs(vinv)):
            line_error('mdd_mcd',
                       'the level marginals converged to an invariant value of %.6g against the '
                       'model value %g, so the fixed point reached is degenerate (the level chains '
                       'are mutually inconsistent). Supply initpik with a consistent starting law.'
                       % (got, vinv))

    pathsPerLevel, noAggregation = _exactness(above, K)
    out = {
        'QLen': QLen,
        'X': X,
        'U': U,
        'pik': pik,
        'Mrows': Mrows,
        'levelSizes': levelSizes,
        'iters': iters,
        # max |A(p)| per paper level: 1 means no node at that level is shared,
        # so conditioning on the node equals conditioning on the path
        'pathsPerLevel': pathsPerLevel,
        # True certifies the result is EXACT with no reference solve; False
        # means "not certified", not "approximate" (see _exactness)
        'noAggregation': noAggregation,
    }

    if verbose:
        line_printf('\nMDD-MCD approximate aggregation: %d levels, %d fixed-point iters\n'
                    % (K, iters))
        line_printf('  level-CTMC sizes |M_k| = [%s] (max %d)\n'
                    % (' '.join(str(int(x)) for x in levelSizes), int(levelSizes.max())))
        line_printf('  mean queue lengths     = %s\n' % np.array2string(QLen, precision=5))
    return out


# ---------------------------------------------------------------------------
def _row_nonzeros(W, d):
    """Per-row (column, value) lists and row sums of a local matrix W.

    Accepts a dense array or any scipy sparse matrix, so a descriptor may build
    its W in whichever form is natural for it.
    """
    if sp is not None and sp.issparse(W):
        Wc = W.tocsr()
        rows = []
        for i in range(d):
            lo, hi = Wc.indptr[i], Wc.indptr[i + 1]
            rows.append(list(zip(Wc.indices[lo:hi].tolist(), Wc.data[lo:hi].tolist())))
        lam = np.asarray(Wc.sum(axis=1)).ravel()
        return rows, lam
    Wd = np.asarray(W, dtype=float)
    rows = [[(int(j), float(Wd[i, j])) for j in np.nonzero(Wd[i, :])[0]] for i in range(d)]
    return rows, Wd.sum(axis=1)


# ---------------------------------------------------------------------------
def _invariant(desc, K):
    """Conservation law the converged marginals must satisfy, as w'*QLen = val."""
    inv = desc.get('invariant')
    if inv:
        return np.asarray(inv['weights'], dtype=float).ravel(), float(inv['value'])
    N = desc.get('N')
    if N is not None:
        return np.ones(K), float(N)          # closed QN: total population
    return None, None


# ---------------------------------------------------------------------------
def _path_counts(mdds, K):
    """(above, below) per MDD level: paths from the root, and states below.

    above[oL][p] is the number of distinct root-to-p paths, i.e. |A(p)| in the
    paper's notation, and below[oL][p] the number of accepted states under p.
    Both are computed in O(#nodes) and serve two purposes: the uniform
    initialisation, and the exactness certificate (see _exactness).
    """
    below = [None] * K
    above = [None] * K
    for oL in range(K - 1, -1, -1):
        tbl = np.asarray(mdds.node[oL], dtype=np.int64).reshape(int(mdds.nnodes[oL]),
                                                                int(mdds.domain[oL]))
        nb = np.zeros(int(mdds.nnodes[oL]))
        for v in range(int(mdds.domain[oL])):
            ch = tbl[:, v]
            if oL == K - 1:
                nb += (ch == TERM_TRUE).astype(float)
            else:
                nz = ch > 0
                if np.any(nz):
                    nb[nz] += below[oL + 1][ch[nz] - 1]
        below[oL] = nb
    for oL in range(K):
        above[oL] = np.zeros(int(mdds.nnodes[oL]))
    above[0][int(mdds.root) - 1] = 1.0
    for oL in range(K - 1):
        tbl = np.asarray(mdds.node[oL], dtype=np.int64).reshape(int(mdds.nnodes[oL]),
                                                                int(mdds.domain[oL]))
        for v in range(int(mdds.domain[oL])):
            ch = tbl[:, v]
            nz = ch > 0
            if not np.any(nz):
                continue
            above[oL + 1] += np.bincount(ch[nz] - 1, weights=above[oL][nz],
                                         minlength=int(mdds.nnodes[oL + 1]))
    return above, below


def _exactness(above, K):
    """Structural certificate that the aggregation loses nothing.

    The single approximation is Pr{i_k | alpha} = Pr{i_k | p}: the local-state
    law is conditioned on the NODE rather than on the whole path above it. When
    a node is reached by exactly one path, |A(p)| = 1, conditioning on the node
    IS conditioning on the path and the identity is exact. If that holds at
    every node the fixed point is the exact stationary law, with no reference
    solve needed to know it.

    This is SUFFICIENT, not necessary: a product-form model is exact too (paper
    Sec. 5) however much its diagram shares, and that is a property of the
    model, not of the diagram. So a False here means "not certified by this
    test", never "approximate". Note also that max |A(p)| = 1 means no node is
    shared, i.e. the diagram compresses nothing -- exactness by this route and
    a useful saving are mutually exclusive.

    Returns (maxPathsPerLevel in paper order, noAggregation).
    """
    per_level = np.ones(K)
    for k in range(K):
        oL = K - 1 - k
        a = above[oL]
        per_level[k] = float(a.max()) if a.size else 1.0
    return per_level, bool(np.all(per_level <= 1.0 + 1e-12))


def _uniform_init(mdds, Mrows, K, above=None, below=None):
    """Uniform law over the EXACT reachable set, projected onto each level.

    Pr{(p,v)} = (paths root->p) * (states below arc p[v]) / |S|. A flat law over
    M_k instead treats level states as equiprobable irrespective of how many
    global states they stand for, which breaks the population invariant the
    diagram encodes; from about K=8 the coupled iteration then descends into the
    basin of the DEGENERATE empty-population fixed point (all mass on local
    state 0 at every level, a true fixed point since no station can then emit)
    and converges to it with zero residual. The projection below is consistent
    across levels by construction, so the iteration starts inside the physical
    simplex.
    """
    if above is None or below is None:
        above, below = _path_counts(mdds, K)
    pik = [None] * K
    for k in range(K):
        oL = K - 1 - k                       # paper level k is MDD level K-1-k
        rows = Mrows[k]
        p = rows[:, 0]
        v = rows[:, 1]
        if oL == K - 1:
            w = above[oL][p - 1]             # a TRUE arc stands for one state
        else:
            tbl = np.asarray(mdds.node[oL], dtype=np.int64).reshape(int(mdds.nnodes[oL]),
                                                                    int(mdds.domain[oL]))
            ch = tbl[p - 1, v]
            w = above[oL][p - 1] * below[oL + 1][ch - 1]
        pik[k] = w / w.sum()
    return pik


# ---------------------------------------------------------------------------
def _compute_as(k, Aup, Pnode, pik, wrows, Mrows, nn, E):
    """ComputeAs(k): A_k^e from A_{k+1}^e, the "from above" contribution (Fig. 3).

    The adjust denominator Pr{p[v]} is the FROM-ABOVE marginal of the child
    node, Pr{p} = sum_{alpha in A(p)} Pr{alpha} = sum over parents of pi_{k+1}.
    Using it (rather than the level-k CTMC marginal, which only equals it at
    convergence) makes adjust a proper conditional Pr{(parent,arc)|child} and
    pins the inter-level node marginals, removing the spurious fixed points.
    """
    rows1 = Mrows[k + 1]
    PrAbove = np.zeros(int(nn[k]))
    for r in range(rows1.shape[0]):
        child = int(Pnode[k + 1][int(rows1[r, 0]) - 1, int(rows1[r, 1])])
        if child > 0:
            PrAbove[child - 1] += pik[k + 1][r]
    Ak = [np.zeros((int(nn[k]), int(nn[k]))) for _ in range(E)]
    for r in range(rows1.shape[0]):
        p = int(rows1[r, 0])
        v = int(rows1[r, 1])
        childp = int(Pnode[k + 1][p - 1, v])            # p[v]: node at level k
        if childp <= 0 or PrAbove[childp - 1] <= 0:
            continue
        adjust = pik[k + 1][r] / PrAbove[childp - 1]    # pi_{k+1}[(p,v)] / Pr{p[v]}
        for e in range(E):
            wrow = wrows[e][k + 1][v]
            if not wrow:
                continue
            arow = Aup[e][p - 1, :]
            qcols = np.nonzero(arow)[0]
            if qcols.size == 0:
                continue
            for w, wv in wrow:
                for q in qcols:
                    childq = int(Pnode[k + 1][q, w])
                    if childq <= 0:
                        continue                         # q[w] null
                    Ak[e][childp - 1, childq - 1] += arow[q] * wv * adjust
    return Ak


# ---------------------------------------------------------------------------
def _compute_mc(k, Ak, bcell, Pnode, wrows, Mrows, Midx, levelSizes, E):
    """ComputeMC(k): level-k rate matrix (Eq. 6).

    R_k^e[(p,i),(q,j)] = A_k^e[p,q] * W_k^e[i,j] * b_{k-1}^e[p[i]]
    """
    nm = int(levelSizes[k])
    rows = Mrows[k]
    Rk = np.zeros((nm, nm))
    for r in range(nm):
        p = int(rows[r, 0])
        v = int(rows[r, 1])
        for e in range(E):
            wrow = wrows[e][k][v]
            if not wrow:
                continue
            if k > 0:
                bfac = bcell[k - 1][int(Pnode[k][p - 1, v]) - 1, e]
            else:
                bfac = 1.0                               # terminal ONE
            if bfac == 0:
                continue
            arow = Ak[e][p - 1, :]
            qcols = np.nonzero(arow)[0]
            if qcols.size == 0:
                continue
            for w, wv in wrow:
                for q in qcols:
                    di = int(Midx[k][q, w])
                    if di == 0:
                        continue                         # (q,w) not in M_k
                    Rk[r, di - 1] += arow[q] * wv * bfac
    return Rk


# ---------------------------------------------------------------------------
def _solve_stat(Q):
    """Stationary distribution of a small irreducible generator: p Q = 0, sum p = 1.

    The normalisation is APPENDED rather than substituted for the last balance
    equation: overwriting a row discards a constraint and leaves Q' singular to
    working precision from about |M_k| = 325 upwards, so the solve returns NaN.
    The overdetermined system has full column rank whenever the level chain is
    irreducible, and QR least squares solves it stably at the same cost order.
    """
    n = Q.shape[0]
    if n == 1:
        return np.array([1.0])
    A = np.vstack([Q.T, np.ones((1, n))])
    rhs = np.zeros(n + 1)
    rhs[n] = 1.0
    p, _, _, _ = np.linalg.lstsq(A, rhs, rcond=None)
    p[p < 0] = 0.0
    s = p.sum()
    if not np.isfinite(s) or s <= 0:
        line_error('mdd_mcd',
                   'level CTMC of order %d admits no proper stationary distribution (the level '
                   'generator is reducible or numerically degenerate).' % n)
    return p / s
