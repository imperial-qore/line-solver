"""
Linear-programming bounds on the mean marking and the throughputs of a
stochastic timed Petri net.

The stationary chain is relaxed to a MOMENT POLYTOPE: the uniformized evolution
equation is written for E[X_p], E[X_p^2] and E[X_p1 X_p2], which gives linear
equalities among the mean marking x, the enabling probabilities q and the
products y(p,t) = E[X_p e_t]; behavioural and probabilistic inequalities are
added on top; and every reported measure is then obtained by minimising and
maximising its linear form over that polytope. Any stationary point of the true
chain satisfies every row, so the two optima BRACKET the exact value whatever
the polytope leaves out.

This is the Petri-net sibling of the QRF bounds in SolverBA: same technique, a
different index space, and a LINEAR objective, so there is no stationary point
to escape from and the answer is a property of the model alone.

VARIABLES, over place levels l = 0..L-1 and modes e = 0..E-1:

  x(l)     E[X_l], mean tokens                            >= 0
  q(e)     P(mode e enabled)                              in [0,1]
  th(e)    throughput of mode e                           >= 0
  u(e)     state-equation firing counts                   >= 0
  y(l,e)   E[X_l e_e]                                     >= 0   (Markovian only)

u is EXISTENTIAL and is not reported: E[X] is a convex combination of reachable
markings, each of which is m0 + C h for some nonnegative integer h, so the mean
satisfies m0 + C u for some nonnegative real u. It is not a mean firing count
and has no steady-state value.

LEVELS ARE (place, class) PAIRS, PLACE-MAJOR, level pp*R + k, the same
coordinates spn_mdd, spn_sinvariants and spn_conv use. MODES are (transition,
mode) pairs in node order.

THE TOKEN COUNTS CREATED BY A FIRING ARE DETERMINISTIC IN LINE, which removes a
whole branch of the reference: it allows sigma_{t,p}(n) to be random and splits
the covariance family into an independent case (its eq. 7) and a selective one
(its eq. 8). set_firing_outcome takes an integer weight, so E[sigma^2] =
sigma^2 and E[sigma_p1 sigma_p2] = sigma_p1 sigma_p2 hold exactly and eq. (7)
is the correct form. Eq. (8) has no LINE model behind it and is deliberately
absent.

LIVENESS IS OFF BY DEFAULT, AND THAT IS DELIBERATE. The reference's two
liveness rows (sum_t q_t >= 1 and x_p <= sum_t y_{p,t}) hold only on a live
net, and liveness is not something this function can cheaply certify -- an
inhibitor arc alone is enough to deadlock a net that looks well formed. A bound
that silently assumed it would be wrong rather than loose on exactly the models
where a bound is most wanted, so the rows are opt-in.

WHAT THE ROWS ARE WORTH, MEASURED. They are the whole of the lower side. On the
reference's own Table 2 (its Fig. 2b production line, five rate vectors)
assumelive reproduces its published l.b. column to four decimals -- 1.1653
against 1.165, 1.8288 against 1.829, 1.5814 against 1.581, 1.3592 against
1.359, 1.3497 against 1.350 -- while without them the Markovian lower bound
collapses onto the OPERATIONAL one (0.9302, 1.4815, 1.1111, 1.1110 against that
column's 0.930, 1.481, 1.111, 1.111) on four of the five. The upper side needs
neither row and matches the published u.b.2 either way.

THE LP IS ASSEMBLED HERE RATHER THAN THROUGH MapqnLpModel. That helper is
generic enough, but it re-assembles its sparse matrices on every solve and this
bound issues one solve per reported cell; it also raises on an infeasible LP
where the MATLAB reference returns NaN. Building the triplets directly keeps
the four ports structurally identical and the failure semantics the same.

Reference:
    Z. Liu, "Performance Analysis of Stochastic Timed Petri Nets Using Linear
    Programming Approach", IEEE Trans. Software Engineering 24(11), 1998,
    1014-1030. The constraint families are its Table 1, p. 1022; the bracket
    statement is its Theorem 3, p. 1021.

See also: spn_sinvariants, spn_metrics, spn_mdd.
"""

# Copyright (c) 2012-2026, Imperial College London
# All rights reserved.

from typing import Dict, List, Optional, Sequence

import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

from ..io.logging import line_error
from ..sn.network_struct import NodeType
from .sinvariants import spn_sinvariants


def _param(nparam, name, default=None):
    if nparam is None:
        return default
    if isinstance(nparam, dict):
        return nparam.get(name, default)
    return getattr(nparam, name, default)


def _timing_name(t) -> str:
    # TimingStrategy lives in lang.nodes, which the api layer does not import;
    # spn_mdd compares by name for the same reason.
    if t is None:
        return ''
    if hasattr(t, 'name'):
        return str(t.name).upper()
    return str(t).upper().split('.')[-1]


def _arcvec(mat, places, R, nnodes, fillval):
    """(nnodes x nclasses) arc matrix -> length P*R place-major level vector.

    The row index is a NODE index; fillval is 0 for enabling and firing and inf
    for inhibiting, where inf means "no arc" and 0 would mean "inhibited
    always".
    """
    P = len(places)
    v = np.full(P * R, fillval, dtype=float)
    if mat is None:
        return v
    m = np.asarray(mat, dtype=float).reshape(nnodes, R)
    for pp in range(P):
        for k in range(R):
            x = m[places[pp], k]
            if fillval == 0.0:
                x = max(0.0, x)
            v[pp * R + k] = x
    return v


def _modes(sn, places, R, markovian):
    """The (transition, mode) table over place-major levels.

    spn_mdd builds the same table, but only as a step of reachable-set
    construction, which is the cost this bound exists to avoid.
    """
    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    transitions = [int(i) for i in np.nonzero(nodetype == int(NodeType.TRANSITION))[0]]
    nnodes = int(sn.nnodes)
    md: List[Dict] = []
    for ind in transitions:
        nparam = sn.nodeparam[ind] if sn.nodeparam is not None else None
        nmodes = int(_param(nparam, 'nmodes', 0) or 0)
        enabling = _param(nparam, 'enabling', None) or []
        inhibiting = _param(nparam, 'inhibiting', None) or []
        firing = _param(nparam, 'firing', None) or []
        firingproc = _param(nparam, 'firingproc', None) or []
        firingpie = _param(nparam, 'firingpie', None) or []
        firingdep = _param(nparam, 'firingdep', None) or []
        nmodeservers = np.ravel(np.asarray(
            _param(nparam, 'nmodeservers', np.ones(nmodes)), dtype=float))
        timing = _param(nparam, 'timingstrategies', None)
        if timing is None:
            timing = _param(nparam, 'timing', None) or []
        for m in range(nmodes):
            if m < len(timing) and _timing_name(timing[m]) == 'IMMEDIATE':
                line_error('spn_lpbnd',
                           'mode %d of node %d is IMMEDIATE; the moment relaxation is written '
                           'for a net whose transitions all have finite rates, so vanishing '
                           'states must be eliminated first' % (m + 1, ind + 1))
            if m < len(firingdep) and firingdep[m] is not None:
                line_error('spn_lpbnd',
                           'mode %d of node %d has a marking-dependent firing rate; the '
                           'uniformization step needs one rate per mode' % (m + 1, ind + 1))
            if m < nmodeservers.size and nmodeservers[m] != 1:
                line_error('spn_lpbnd',
                           'mode %d of node %d has %g servers; the relaxation is derived under '
                           'single-server semantics, where the firing rate is mu*q. Its '
                           'infinite-server form needs the K-fold transition expansion of the '
                           "reference's Section 7, which is not implemented"
                           % (m + 1, ind + 1, nmodeservers[m]))
            # A PHASE-TYPE FIRING LAW IS WHERE THE TWO VARIANTS PART. The mean
            # of a (D0,D1) pair is pie*(-D0)^-1*1 and is all the operational
            # bound needs; the Markovian one needs the marking alone to be the
            # state, which a multi-phase mode breaks. A law that is not
            # phase-type at all reaches sn as its own parameter list (Pareto
            # stores (shape, scale)), with no mean recoverable without a
            # per-distribution table, so BOTH variants refuse it.
            proc = firingproc[m] if m < len(firingproc) else None
            D0 = None
            if proc is not None and len(proc) >= 2:
                cand = np.asarray(proc[0], dtype=float)
                if cand.ndim == 2 and cand.shape[0] == cand.shape[1]:
                    D0 = cand
                    D1 = np.asarray(proc[1], dtype=float)
            if D0 is None:
                line_error('spn_lpbnd',
                           'mode %d of node %d has a firing law that is not phase-type; sn '
                           'carries its own parameters rather than a (D0,D1) pair, so neither '
                           'the Markovian nor the operational bound can read a mean firing rate '
                           'from it. Use a phase-type law, or an exact solver' % (m + 1, ind + 1))
            nph = D0.shape[0]
            if markovian and nph > 1:
                line_error('spn_lpbnd',
                           'mode %d of node %d has a phase-type firing time; the relaxation is '
                           'written over the marking alone, and a phase-type mode needs the '
                           "state-machine expansion of the reference's Section 7, which is not "
                           'implemented. Use the operational bound, which needs only the mean'
                           % (m + 1, ind + 1))
            if nph == 1:
                rate = float(np.ravel(D1)[0])
            else:
                pv = firingpie[m] if m < len(firingpie) else None
                pv = np.ones(nph) / nph if pv is None or len(pv) == 0 \
                    else np.ravel(np.asarray(pv, dtype=float))
                pv = pv / pv.sum()
                rate = 1.0 / float(pv @ np.linalg.solve(-D0, np.ones(nph)))
            if not np.isfinite(rate) or rate <= 0:
                line_error('spn_lpbnd',
                           'mode %d of node %d has mean firing rate %g; a bound needs a finite '
                           'positive one' % (m + 1, ind + 1, rate))
            md.append({
                'trans': ind,
                'mode': m,
                'enab': _arcvec(enabling[m] if m < len(enabling) else None,
                                places, R, nnodes, 0.0),
                'inhib': _arcvec(inhibiting[m] if m < len(inhibiting) else None,
                                 places, R, nnodes, np.inf),
                'fire': _arcvec(firing[m] if m < len(firing) else None,
                                places, R, nnodes, 0.0),
                'rate': rate,
            })
    if not md:
        line_error('spn_lpbnd', 'the net has no firing mode')
    return md


def _initmarking(sn, places, R):
    """Declared tokens per (place, class), or None when nothing is set."""
    nodetostateful = np.ravel(np.asarray(sn.nodeToStateful, dtype=int)) \
        if getattr(sn, 'nodeToStateful', None) is not None else None
    state = getattr(sn, 'state', None)
    if nodetostateful is None or state is None:
        return None
    init = np.zeros(len(places) * R)
    anyset = False
    for pp, nd in enumerate(places):
        isf = int(nodetostateful[nd])
        if isf < 0 or isf >= len(state) or state[isf] is None:
            continue
        row = np.ravel(np.asarray(state[isf], dtype=float))
        if row.size == 0:
            continue
        row = np.asarray(state[isf], dtype=float)
        row = row[0, :] if row.ndim == 2 else row
        for k in range(min(R, row.size)):
            init[pp * R + k] = row[k]
            anyset = anyset or row[k] > 0
    return init if anyset else None


def _levelbounds(S, V, L):
    """Tightest a priori bound per level, from the P-invariants."""
    B = np.full(L, np.inf)
    for i in range(len(V)):
        for l in range(L):
            if S[i][l] > 0:
                B[l] = min(B[l], np.floor(V[i] / S[i][l]))
    return B


class _Rows:
    """Triplet accumulator for one relation family."""

    def __init__(self):
        self.r: List[int] = []
        self.c: List[int] = []
        self.v: List[float] = []
        self.b: List[float] = []
        self.n = 0

    def add(self, idx, val, rhs):
        idx = np.asarray(idx, dtype=np.intp).ravel()
        val = np.asarray(val, dtype=float).ravel()
        keep = val != 0
        idx = idx[keep]
        val = val[keep]
        if idx.size == 0:
            return
        self.r.extend([self.n] * idx.size)
        self.c.extend(idx.tolist())
        self.v.extend(val.tolist())
        self.b.append(float(rhs))
        self.n += 1

    def matrix(self, nv):
        if self.n == 0:
            return None, None
        # duplicate (row, col) entries are summed by COO -> CSR, which is what
        # supplies the factor of two the l1 == l2 covariance row carries
        A = coo_matrix((self.v, (self.r, self.c)), shape=(self.n, nv)).tocsr()
        return A, np.asarray(self.b, dtype=float)


def spn_lpbnd(sn, options: Optional[Dict] = None) -> Dict[str, object]:
    """Bracket the mean tokens and the throughputs of a stochastic Petri net.

    Parameters
    ----------
    sn : a NetworkStruct holding Places and Transitions
    options : dict, all keys optional
        markovian  True (default) uses the second-moment, covariance and
            Little's law families, which need exponential firing times; False
            drops them and the whole y block, leaving the operational bound,
            which needs only a mean firing time and so admits any phase-type
            law. A law that is not phase-type at all is refused by both.
        assumelive  False (default). True adds the two liveness rows, which are
            valid only on a live net.
        init  initial tokens per place level, place-major
        tol  slack added to the inequality sides, default 0
        verbose  print the polytope size

    Returns
    -------
    dict with 'places', 'levelname', 'modes', 'tokens', 'placeTput',
    'modeTput', 'modeUtil' (each a 2 x n array, row 0 the minimum and row 1 the
    maximum), 'bound', 'nplacelevels', 'nclasses', 'markovian', 'nvars',
    'nrows'.
    """
    options = options or {}
    markovian = bool(options.get('markovian', True))
    assumelive = bool(options.get('assumelive', False))
    tol = float(options.get('tol', 0.0))
    verbose = bool(options.get('verbose', False))

    nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
    places = [int(i) for i in np.nonzero(nodetype == int(NodeType.PLACE))[0]]
    if not places or not np.any(nodetype == int(NodeType.TRANSITION)):
        line_error('spn_lpbnd', 'the model holds no Place or no Transition node')
    R = int(sn.nclasses)
    P = len(places)
    L = P * R

    md = _modes(sn, places, R, markovian)
    E = len(md)

    PI = np.zeros((E, L))
    SG = np.zeros((E, L))
    ETA = np.full((E, L), np.inf)
    for e in range(E):
        PI[e, :] = md[e]['enab']
        SG[e, :] = md[e]['fire']
        ETA[e, :] = md[e]['inhib']
    NET = SG - PI
    mu = np.array([m['rate'] for m in md], dtype=float)

    # Per-level a priori bounds and the conserved sums, both off the same
    # minimal-support P-invariant basis. The reference writes its "cycle
    # population" family for UNWEIGHTED cycles; spn_sinvariants returns the
    # weighted invariants S m = V, which are equally linear and strictly
    # tighter, so those are what is emitted.
    init = options.get('init', None)
    if init is None:
        init = _initmarking(sn, places, R)
    inv = spn_sinvariants(sn, init)
    S = np.asarray(inv['S'], dtype=float).reshape(-1, L) if len(inv['S']) else np.zeros((0, L))
    V = np.asarray(inv['V'], dtype=float).ravel()
    m0 = np.asarray(inv['m0'], dtype=float).ravel()
    B = _levelbounds(S, V, L)

    # ---- variable layout
    ix = np.arange(L)
    iq = L + np.arange(E)
    ith = L + E + np.arange(E)
    iu = L + 2 * E + np.arange(E)
    nv = L + 3 * E
    if markovian:
        iy = (nv + np.arange(L * E)).reshape(L, E)
        nv += L * E
    else:
        iy = None

    bounds: List = [(0.0, None)] * nv
    for e in range(E):
        bounds[iq[e]] = (0.0, 1.0)
    for l in range(L):
        if np.isfinite(B[l]):
            bounds[ix[l]] = (0.0, float(B[l]))
            if markovian:
                for e in range(E):
                    bounds[iy[l, e]] = (0.0, float(B[l]))

    eq = _Rows()
    ub = _Rows()

    def le(idx, val, rhs):
        ub.add(idx, val, rhs)

    def ge(idx, val, rhs):
        ub.add(idx, -np.asarray(val, dtype=float), -rhs)

    # ---- (1) throughput: th_e = mu_e q_e
    for e in range(E):
        eq.add([ith[e], iq[e]], [1.0, -mu[e]], 0.0)

    # ---- (2) flow balance: tokens are created at a level at the rate they are
    # consumed there. Holds for any stable net, Markovian or not.
    for l in range(L):
        eq.add(iq, mu * NET[:, l], 0.0)

    # ---- (3)+(4) second moment and population covariance, from the
    # stationarity of E[X_l1 X_l2] under the uniformized chain. Table 1 writes
    # the q side as four sums over set intersections; since the memberships are
    # exactly "sigma > 0" and "pi > 0", those collapse to
    #     -(sigma_1 - pi_1)(sigma_2 - pi_2) = -net_1 net_2
    # per mode, which also makes the l1 == l2 case reduce to the second-moment
    # family with no separate derivation, as the reference notes it must.
    if markovian:
        for l1 in range(L):
            for l2 in range(l1, L):
                eq.add(np.concatenate([iy[l1, :], iy[l2, :], iq]),
                       np.concatenate([mu * NET[:, l2], mu * NET[:, l1],
                                       mu * NET[:, l1] * NET[:, l2]]),
                       0.0)

    # ---- (5) liveness, only when the caller vouches for it
    if assumelive:
        ge(iq, np.ones(E), 1.0)
        if markovian:
            for l in range(L):
                le(np.concatenate([[ix[l]], iy[l, :]]),
                   np.concatenate([[1.0], -np.ones(E)]), 0.0)

    # ---- (6) conflicting transitions: a mode that consumes no more and is
    # inhibited no sooner is enabled whenever the other is
    for e1 in range(E):
        for e2 in range(E):
            if e1 == e2:
                continue
            if np.all(PI[e1, :] <= PI[e2, :]) and np.all(ETA[e1, :] >= ETA[e2, :]):
                ge([iq[e1], iq[e2]], [1.0, -1.0], 0.0)

    # ---- (7) boundedness, per level; and (8) cycle population, as the
    # weighted invariant equalities and their y companions
    if markovian:
        for l in range(L):
            if not np.isfinite(B[l]):
                continue
            for e in range(E):
                le([iy[l, e], iq[e]], [1.0, -B[l]], 0.0)
                le([ix[l], iy[l, e], iq[e]], [1.0, -1.0, B[l]], B[l])
                if B[l] > 0:
                    ge([ix[l], iy[l, e], iq[e]], [1.0 - 1.0 / B[l], -1.0, 1.0], 0.0)
    for i in range(S.shape[0]):
        eq.add(ix, S[i, :], V[i])
        if markovian:
            for e in range(E):
                eq.add(np.concatenate([iy[:, e], [iq[e]]]),
                       np.concatenate([S[i, :], [-V[i]]]), 0.0)

    # ---- (9) reachable marking: the mean lies in the state-equation cone
    for l in range(L):
        eq.add(np.concatenate([[ix[l]], iu]),
               np.concatenate([[1.0], -NET[:, l]]), m0[l])

    # ---- (10) sample-path comparisons
    if markovian:
        mutot = float(mu.sum())
        for l in range(L):
            for e in range(E):
                le([iy[l, e], ix[l]], [1.0, -1.0], 0.0)
                if PI[e, l] > 0:
                    ge([iy[l, e], iq[e]], [1.0, -PI[e, l]], 0.0)
                if np.isfinite(ETA[e, l]):
                    le([iy[l, e], iq[e]], [1.0, -(ETA[e, l] - 1.0)], 0.0)
            ge(np.concatenate([[ix[l]], iy[l, :]]),
               np.concatenate([[mutot], -mu]), 0.0)
        for e in range(E):
            ent = np.nonzero(PI[e, :] > 0)[0]
            if ent.size == 1 and not np.any(np.isfinite(ETA[e, :])):
                l = int(ent[0])
                le([ix[l], iy[l, e]], [1.0, -1.0], PI[e, l] - 1.0)

    # ---- (11) enabling bounds, from Chernoff's inequality on the marking
    for e in range(E):
        ent = np.nonzero(PI[e, :] > 0)[0]
        inh = np.nonzero(np.isfinite(ETA[e, :]))[0]
        D = ent.size + inh.size
        if D == 0:
            continue
        if ent.size and np.all(np.isfinite(B[ent])):
            idx = [iq[e]]
            val = [1.0]
            rhs = 1.0
            ok = True
            for l in ent:
                den = B[l] - PI[e, l] + 1.0
                if den <= 0:
                    ok = False
                    break
                idx.append(ix[l])
                val.append(-1.0 / den)
                rhs -= B[l] / den
            if ok:
                for l in inh:
                    idx.append(ix[l])
                    val.append(1.0 / ETA[e, l])
                ge(idx, val, rhs)
        # Upper side. Every term of the sum over input levels carries a "min"
        # operator, and Table 1's convention is that either operand may be
        # taken; each choice is a valid row and the whole set is the tightest
        # linear relaxation, so all of them are emitted while the count stays
        # small.
        ok = True
        for l in inh:
            if not np.isfinite(B[l]) or B[l] - ETA[e, l] + 1.0 <= 0:
                ok = False
                break
        if not ok:
            continue
        base = 0.0
        bidx: List[int] = []
        bval: List[float] = []
        for l in inh:
            den = B[l] - ETA[e, l] + 1.0
            base += B[l] / den
            bidx.append(ix[l])
            bval.append(1.0 / den)
        nc = ent.size
        combos = range(1 << nc) if nc <= 4 else [0, (1 << nc) - 1]
        for c in combos:
            idx = [iq[e]] + bidx
            val = [float(D)] + bval
            rhs = base
            for j in range(nc):
                l = int(ent[j])
                if (c >> j) & 1 == 0:
                    idx.append(ix[l])
                    val.append(-1.0 / PI[e, l])
                else:
                    rhs += 1.0
            le(idx, val, rhs)

    # ---- (12) Little's law at each level: the mean sojourn time of a token is
    # at least the mean minimum firing time of the modes that can remove it
    if markovian:
        for l in range(L):
            out = float(mu[PI[:, l] > 0].sum())
            if out <= 0:
                continue
            ge(np.concatenate([[ix[l]], iq]),
               np.concatenate([[out], -(mu * SG[:, l])]), 0.0)

    A_eq, b_eq = eq.matrix(nv)
    A_ub, b_ub = ub.matrix(nv)
    if b_ub is not None and tol > 0:
        b_ub = b_ub + tol

    if verbose:
        print('\nSPN -> LP: %d place levels, %d modes, %d variables, %d equalities, '
              '%d inequalities' % (L, E, nv, eq.n, ub.n))

    def solve(c, minimize):
        res = linprog(c if minimize else -c, A_ub=A_ub, b_ub=b_ub,
                      A_eq=A_eq, b_eq=b_eq, bounds=bounds, method='highs',
                      options={'maxiter': 100000})
        # success is the wrong predicate on its own; the finiteness of the
        # answer is what decides, matching the MATLAB reference
        if not res.success or res.x is None or not np.all(np.isfinite(res.x)):
            return np.nan
        return res.fun if minimize else -res.fun

    def bracket(idx, val):
        c = np.zeros(nv)
        c[np.asarray(idx, dtype=np.intp).ravel()] = val
        return solve(c, True), solve(c, False)

    tokens = np.full((2, L), np.nan)
    placeTput = np.full((2, L), np.nan)
    for l in range(L):
        tokens[0, l], tokens[1, l] = bracket([ix[l]], 1.0)
        if np.any(PI[:, l] > 0):
            placeTput[0, l], placeTput[1, l] = bracket(iq, mu * PI[:, l])
        else:
            placeTput[:, l] = 0.0
    modeTput = np.full((2, E), np.nan)
    modeUtil = np.full((2, E), np.nan)
    for e in range(E):
        modeTput[0, e], modeTput[1, e] = bracket([ith[e]], 1.0)
        modeUtil[0, e], modeUtil[1, e] = bracket([iq[e]], 1.0)

    levelname = []
    for pp in range(P):
        for k in range(R):
            levelname.append('%s.%s' % (sn.nodenames[places[pp]], sn.classnames[k]))

    return {
        'places': places,
        'levelname': levelname,
        'modes': md,
        'tokens': tokens,
        'placeTput': placeTput,
        'modeTput': modeTput,
        'modeUtil': modeUtil,
        'bound': B,
        'nplacelevels': L,
        'nclasses': R,
        'markovian': markovian,
        'nvars': nv,
        'nrows': eq.n + ub.n,
    }


__all__ = ['spn_lpbnd']
