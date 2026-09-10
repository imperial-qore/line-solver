"""Fluid analysis of a stochastic Petri net, on the 'dae' method.

Port of matlab/src/solvers/FLD/fluid_petri_*.m and solver_fluid_petri.m. A GSPN
is already a density-dependent Markov population process -- the marking is the
population, a transition mode is a reaction, its incidence column is the jump,
and the rate law lambda*min(enabling degree, servers) is the same min()
non-linearity the min-normal closure exists to smooth -- so nothing about the
closure changes here, only where the drift comes from.

    dx/dt = D * r(x, Sigma, phi, mu)

See _kb/06-solver-catalog.md for the formulation, the measured accuracy and the
one case (an immediate mode under an inhibitor arc) that is answered wrongly.
"""
import numpy as np

from ....lang.base import NodeType
from line_solver.constants import GlobalConstants

FINE_TOL = 1e-8
COARSE_TOL = 1e-3


# ------------------------------------------------------------------ closures
def min_closure(n, c, s2, vc=0.0, cov=0.0):
    """E[min(X,Y)] and dE/dE[X] for jointly normal X, Y (fluid_min_closure)."""
    if np.isinf(c):
        return float(n), 1.0
    th2 = s2 - 2.0 * cov + vc
    if th2 <= 0.0:
        # A degenerate pair: the min is the smaller of the two exactly. The band
        # matches the MATLAB twin, so two codebases stopping either side of the
        # kink read the same indicator.
        if c - n > FINE_TOL * max(1.0, abs(n)):
            return float(n), 1.0
        return float(c), 0.0
    th = np.sqrt(th2)
    al = (n - c) / th
    Phi = 0.5 * _erfc(-al / np.sqrt(2.0))
    phi = np.exp(-0.5 * al * al) / np.sqrt(2.0 * np.pi)
    p = 1.0 - Phi
    h = n * p + c * Phi - th * phi
    return float(h), float(p)


def minmulti_closure(mu, S, c=np.inf):
    """Min-normal closure of E[min(X_1,...,X_A,c)] by Clark's (1961) recursion.

    Returns (h, g, v): the expectation, dH/dMU(a) -- the probability that arc a
    is the binding one -- and Var[min] before the cap.

    The recursion is ORDER DEPENDENT, as Clark's approximation always is: only
    the first two moments of the running min are kept. The order is the
    caller's, i.e. increasing state coordinate, which the layout fixes.
    """
    mu = np.asarray(mu, dtype=float).ravel()
    A = mu.size
    # A mode with no input arc is enabled at degree one, the convention the
    # exact engines use (solver_ssa_nrm's spnEnDegree, State.afterGlobalEvent).
    if A == 0:
        return float(min(1.0, c)), np.zeros(0), 0.0
    S = np.zeros((A, A)) if S is None or np.size(S) == 0 else np.asarray(S, float)

    mz = mu[0]
    vz = S[0, 0]
    covz = S[0, :].copy()
    p = np.ones(A)

    for k in range(1, A):
        th2 = vz + S[k, k] - 2.0 * covz[k]
        if th2 < 0.0:
            th2 = 0.0  # a covariance beyond Cauchy-Schwarz is not admissible
        th = np.sqrt(th2)
        if th > 0.0:
            al = (mz - mu[k]) / th
            Phi = 0.5 * _erfc(-al / np.sqrt(2.0))
            phi = np.exp(-0.5 * al * al) / np.sqrt(2.0 * np.pi)
            pk = 1.0 - Phi
            mw = mz * pk + mu[k] * Phi - th * phi
            e2 = ((mz * mz + vz) * pk + (mu[k] * mu[k] + S[k, k]) * Phi
                  - (mz + mu[k]) * th * phi)
            vw = e2 - mw * mw
            # Clark's moment match can leave a negative variance where the two
            # arguments are nearly identical; the min of two equal normals has
            # the variance of either, which is what the clamp restores.
            if vw < 0.0:
                vw = 0.0
            cw = covz * pk + S[:, k] * Phi
        else:
            if mu[k] - mz > FINE_TOL * max(1.0, abs(mz)):
                pk = 1.0
            else:
                pk = 0.0
            mw = min(mz, mu[k])
            vw = pk * vz + (1.0 - pk) * S[k, k]
            cw = covz * pk + S[:, k] * (1.0 - pk)
        p[k] = pk
        mz, vz, covz = mw, vw, cw

    v = vz
    # The cap, by the two-argument closure itself: c is deterministic, so its
    # variance and its covariance with the running min are both zero.
    h, pcap = min_closure(mz, c, vz, 0.0, 0.0)

    g = np.zeros(A)
    tail = pcap
    for a in range(A - 1, 0, -1):
        g[a] = (1.0 - p[a]) * tail
        tail = tail * p[a]
    g[0] = tail
    return float(h), g, float(v)


def _erfc(z):
    from scipy.special import erfc
    return erfc(z)


# -------------------------------------------------------------------- terms
class Mode(object):
    """One transition mode as a reaction record."""

    __slots__ = ('node', 'mode', 'timing', 'arc_slot', 'arc_w', 'inh_slot',
                 'inh_thr', 'c', 'nph', 'D0', 'D1', 'd1', 'pie', 'dep', 'prio',
                 'weight', 'cvec', 'zblk', 'closable', 'label')

    def __init__(self, node, mode, timing, label, nm):
        self.node = node
        self.mode = mode
        self.timing = timing
        self.arc_slot = []
        self.arc_w = []
        self.inh_slot = []
        self.inh_thr = []
        self.c = 1.0
        self.nph = 1
        self.D0 = None
        self.D1 = None
        self.d1 = np.ones(1)
        self.pie = np.ones(1)
        self.dep = None
        self.prio = 1
        self.weight = 1.0
        self.cvec = np.zeros(nm)
        self.zblk = []
        self.closable = False
        self.label = label


class Terms(object):
    """Event-based representation of the fluid marking process."""
    pass


def _pad(A, I, K, fill):
    """Pad an arc matrix to (nnodes x nclasses).

    addMode sizes them at creation time, so a node or class added later leaves
    them short.
    """
    A = np.asarray(A, dtype=float) if A is not None else np.full((0, 0), fill)
    if A.ndim != 2:
        A = A.reshape(-1, 1) if A.size else np.zeros((0, 0))
    if A.shape[0] < I or A.shape[1] < K:
        B = np.full((I, K), float(fill))
        if A.size:
            B[:A.shape[0], :A.shape[1]] = A
        A = B
    return A[:I, :K]


def _mode_timing(tp, m):
    """'TIMED' or 'IMMEDIATE', however this struct spells it."""
    ts = getattr(tp, 'timingstrategies', None)
    if ts is None:
        ts = getattr(tp, 'timing', None)
    if ts is None:
        return 'TIMED'
    v = ts[m]
    # BY NAME, never by str(): TimingStrategy stringifies to its VALUE ('1'),
    # so a substring test on str(v) reads every immediate mode as TIMED and the
    # mode is then refused for having no firing process.
    name = getattr(v, 'name', None)
    if name is not None:
        return 'IMMEDIATE' if str(name).upper() == 'IMMEDIATE' else 'TIMED'
    if isinstance(v, str):
        return 'IMMEDIATE' if v.upper() == 'IMMEDIATE' else 'TIMED'
    from ....constants import TimingStrategy as _TS
    return 'IMMEDIATE' if int(v) == int(_TS.IMMEDIATE.value) else 'TIMED'


def _mode_name(tp, m):
    names = getattr(tp, 'modenames', None)
    if names is not None and len(names) > m and names[m]:
        return str(names[m])
    return 'Mode%d' % (m + 1)


def build_mode(sn, tp, ind, m, pidx, nm, I, K, names):
    """One mode: its input arcs and weights, inhibitor arcs and thresholds, its
    incidence column, and the firing process that times it."""
    rec = Mode(ind, m, _mode_timing(tp, m), '%s.%s' % (names[ind], _mode_name(tp, m)), nm)

    en = _pad(tp.enabling[m], I, K, 0.0)
    fir = _pad(tp.firing[m], I, K, 0.0)
    inh = _pad(tp.inhibiting[m], I, K, np.inf)

    cvec = np.zeros(nm)
    for p, k in zip(*np.nonzero(en > 0)):
        w = en[p, k]
        if not np.isfinite(w):
            raise ValueError(
                "Mode %s has a non-finite enabling arc weight at %s. An arc that no marking can "
                "satisfy disables the mode; declare a finite multiplicity." % (rec.label, names[p]))
        if pidx[p, k] < 0:
            raise ValueError("Mode %s takes an enabling arc from %s, which is not a Place."
                             % (rec.label, names[p]))
        rec.arc_slot.append(int(pidx[p, k]))
        rec.arc_w.append(float(w))
        cvec[pidx[p, k]] -= w
    for p, k in zip(*np.nonzero(fir != 0)):
        if fir[p, k] <= 0:
            continue  # a negative entry marks an input place, whose token the PRE already removed
        if pidx[p, k] < 0:
            continue  # a firing arc into a Sink is mass leaving the net
        cvec[pidx[p, k]] += fir[p, k]
    for p, k in zip(*np.nonzero(np.isfinite(inh))):
        if inh[p, k] <= 0 or pidx[p, k] < 0:
            continue  # JMT writes a missing inhibitor arc as 0 or -1, never as a threshold
        rec.inh_slot.append(int(pidx[p, k]))
        rec.inh_thr.append(float(inh[p, k]))
    rec.cvec = cvec

    c = np.asarray(tp.nmodeservers, dtype=float).ravel()[m]
    rec.c = 1.0 if (c is None or np.isnan(c)) else float(c)
    prio = getattr(tp, 'firingprio', None)
    rec.prio = int(prio[m]) if prio is not None and len(prio) > m else 1
    wgt = getattr(tp, 'fireweight', None)
    rec.weight = float(wgt[m]) if wgt is not None and len(wgt) > m else 1.0
    dep = getattr(tp, 'firingdep', None)
    if dep is not None and len(dep) > m and dep[m]:
        rec.dep = dep[m]

    if rec.timing == 'IMMEDIATE':
        # An immediate mode carries no firing process: its flow is an algebraic
        # unknown of the DAE, not a rate.
        rec.nph = 0
        rec.d1 = np.zeros(1)
        rec.closable = False
        return rec

    proc = tp.firingproc[m] if tp.firingproc is not None else None
    if proc is None or proc[0] is None:
        raise ValueError(
            "Mode %s has no Markovian firing process. sn_nonmarkov_toph converts the renewal "
            "families to phase type before the solver runs, so this is a distribution the fluid "
            "Petri route cannot time; use SolverCTMC or SolverLDES." % rec.label)
    D0 = np.asarray(proc[0], dtype=float)
    D1 = np.asarray(proc[1], dtype=float)
    rec.nph = int(D0.shape[0])
    rec.D0, rec.D1 = D0, D1
    rec.d1 = D1.sum(axis=1)
    if rec.nph > 1:
        pie = tp.firingpie[m] if tp.firingpie is not None else None
        if pie is None or np.size(pie) == 0:
            pie = np.zeros(rec.nph)
            pie[0] = 1.0
        pie = np.asarray(pie, dtype=float).ravel()
        rec.pie = pie / pie.sum()

    # A mode whose enabling degree cannot reach its server count has min() exact
    # on its whole support, so closing it is an error rather than an
    # improvement -- exactly as a station that cannot fill its servers is held
    # first order. A single input arc with an unbounded server count is the
    # common case, and it makes the whole drift LINEAR in that mode.
    rec.closable = not (len(rec.arc_slot) <= 1 and not np.isfinite(rec.c))
    return rec


def build_terms(sn, options=None):
    """Event-based representation of the fluid marking process of an SPN.

    THE STATE, x = [m ; y].
      m(p,k)  token mass of class k at place p, one coordinate per (place,
              class) pair some arc touches, the initial marking loads, or a
              Source feeds; a pair nothing reaches is dropped rather than
              carried as a null direction of the Newton system.
      y(j,h)  the number of mode-j servers running in phase h, for a mode whose
              firing time has more than one phase. Their SUM is not free: the
              ENABLE synchronization latches it to min(enabling degree,
              servers), so the latch is an ALGEBRAIC row with one free-sign
              unknown mu_j and the phase split evolves differentially.
    """
    from ....api.state.marginal import toMarginalAggr

    I = int(len(sn.nodetype))
    K = int(sn.nclasses)
    M = int(sn.nstations)
    names = list(sn.nodenames)

    places = [i for i in range(I) if int(sn.nodetype[i]) == int(NodeType.PLACE)]
    transitions = [i for i in range(I) if int(sn.nodetype[i]) == int(NodeType.TRANSITION)]

    # ---- the initial marking, read off the model state exactly as the NRM does
    m0full = np.zeros((I, K))
    for ind in places:
        isf = int(sn.nodeToStateful[ind])
        out = toMarginalAggr(sn, ind, np.asarray(sn.state[isf]))
        nir = np.asarray(out[1], dtype=float).ravel()
        for k in range(K):
            if np.isinf(nir[k]):
                raise ValueError("Place %s holds an infinite initial marking of class %d."
                                 % (names[ind], k + 1))
            m0full[ind, k] = nir[k]

    # ---- which (place, class) pairs carry a coordinate
    touched = np.zeros((I, K), dtype=bool)
    for ind in transitions:
        tp = sn.nodeparam[ind]
        for m in range(int(tp.nmodes)):
            en = _pad(tp.enabling[m], I, K, 0.0)
            fir = _pad(tp.firing[m], I, K, 0.0)
            inh = _pad(tp.inhibiting[m], I, K, np.inf)
            touched |= (en > 0) | (fir != 0) | (np.isfinite(inh) & (inh > 0))

    # an arrival makes its target a coordinate even when no arc mentions it
    src_arr = []
    for ind in range(I):
        if int(sn.nodetype[ind]) != int(NodeType.SOURCE):
            continue
        ist = int(sn.nodeToStation[ind])
        for r in range(K):
            lam = float(sn.rates[ist, r])
            if not np.isfinite(lam) or lam <= 0:
                continue
            for jnd in places:
                for s in range(K):
                    p = float(sn.rtnodes[ind * K + r, jnd * K + s])
                    if p > 0:
                        touched[jnd, s] = True
                        src_arr.append((ind, r, jnd, s))

    keep = np.zeros((I, K), dtype=bool)
    for ind in places:
        keep[ind, :] = touched[ind, :] | (m0full[ind, :] > 0)

    pidx = -np.ones((I, K), dtype=int)
    coord_node, coord_class, coord_station = [], [], []
    nm = 0
    for ind in places:
        for k in range(K):
            if not keep[ind, k]:
                continue
            pidx[ind, k] = nm
            coord_node.append(ind)
            coord_class.append(k)
            coord_station.append(int(sn.nodeToStation[ind]))
            nm += 1

    # ---- the modes, and the phase coordinates of the multi-phase ones
    modes = []
    nstate = nm
    for ind in transitions:
        tp = sn.nodeparam[ind]
        for m in range(int(tp.nmodes)):
            rec = build_mode(sn, tp, ind, m, pidx, nm, I, K, names)
            if rec.nph > 1:
                rec.zblk = list(range(nstate, nstate + rec.nph))
                nstate += rec.nph
            modes.append(rec)
    # The phase coordinates are appended after every marking coordinate, so a
    # jump column built at nm width has to be GROWN once the total is known.
    for rec in modes:
        if rec.cvec.size < nstate:
            rec.cvec = np.concatenate([rec.cvec, np.zeros(nstate - rec.cvec.size)])

    timed_idx = [j for j, md in enumerate(modes) if md.timing == 'TIMED']
    imm_idx = [j for j, md in enumerate(modes) if md.timing == 'IMMEDIATE']

    # ---- the event columns
    cols, rate_base = [], []
    ev_kind, ev_mode, ev_phase, ev_to, ev_station, ev_class = [], [], [], [], [], []

    def _add(col, base, kind, mode, phase, to, station, cls):
        cols.append(col)
        rate_base.append(base)
        ev_kind.append(kind)
        ev_mode.append(mode)
        ev_phase.append(phase)
        ev_to.append(to)
        ev_station.append(station)
        ev_class.append(cls)

    for j in timed_idx:
        md = modes[j]
        if md.nph == 1:
            _add(md.cvec.copy(), float(md.d1[0]), 1, j, 0, 0, -1, -1)
        else:
            for h in range(md.nph):
                for hp in range(md.nph):
                    w = md.D1[h, hp]
                    if w <= 0:
                        continue
                    col = md.cvec.copy()
                    col[md.zblk[hp]] += 1.0
                    col[md.zblk[h]] -= 1.0
                    _add(col, float(w), 1, j, h, hp, -1, -1)
            for h in range(md.nph):
                for hp in range(md.nph):
                    if hp == h:
                        continue
                    w = md.D0[h, hp]
                    if w <= 0:
                        continue
                    col = np.zeros(nstate)
                    col[md.zblk[hp]] = 1.0
                    col[md.zblk[h]] = -1.0
                    _add(col, float(w), 2, j, h, hp, -1, -1)
    # One latch column per multi-phase mode: mu_j servers per unit time enter at
    # the firing process's own entry distribution. The rate is FREE IN SIGN -- a
    # mode whose enabling degree drops stops servers rather than starting them
    # -- and it is zero at any fixed point, which the latch row enforces.
    for j in timed_idx:
        md = modes[j]
        if md.nph <= 1:
            continue
        col = np.zeros(nstate)
        col[md.zblk] = md.pie
        _add(col, 1.0, 5, j, 0, 0, -1, -1)
    for j in imm_idx:
        _add(modes[j].cvec.copy(), 1.0, 4, j, 0, 0, -1, -1)
    for (snd, r, qnd, l) in src_arr:
        ist = int(sn.nodeToStation[snd])
        if not _is_exponential(sn, ist, r):
            raise ValueError(
                "Source %s has a non-exponential arrival for class %d. The fluid Petri route models "
                "an arrival as a constant-propensity event, which a renewal stream with memory is "
                "not; use SolverCTMC, SolverJMT or SolverSSA." % (names[snd], r + 1))
        col = np.zeros(nstate)
        col[pidx[qnd, l]] = 1.0
        _add(col, float(sn.rates[ist, r]) * float(sn.rtnodes[snd * K + r, qnd * K + l]),
             3, -1, 0, 0, ist, r)

    nev = len(cols)
    D = np.array(cols).T if nev else np.zeros((nstate, 0))
    rate_base = np.asarray(rate_base, dtype=float)
    ev_kind = np.asarray(ev_kind, dtype=int)
    ev_mode = np.asarray(ev_mode, dtype=int)
    ev_phase = np.asarray(ev_phase, dtype=int)
    ev_to = np.asarray(ev_to, dtype=int)
    ev_station = np.asarray(ev_station, dtype=int)
    ev_class = np.asarray(ev_class, dtype=int)

    # ---- which Sigma entries the closure reads: the input coordinates of a
    # closable mode pairwise, plus the variance of every smoothed inhibitor.
    pair_key, cov_pairs = {}, []

    def _addpair(a, b):
        lo, hi = (a, b) if a <= b else (b, a)
        if (lo, hi) not in pair_key:
            pair_key[(lo, hi)] = len(cov_pairs)
            cov_pairs.append((lo, hi))

    for md in modes:
        if md.closable:
            for ai in range(len(md.arc_slot)):
                for bi in range(ai, len(md.arc_slot)):
                    _addpair(md.arc_slot[ai], md.arc_slot[bi])
        if md.timing == 'TIMED':
            for b in md.inh_slot:
                _addpair(b, b)
    npair = len(cov_pairs)

    # ---- per-place consumption and production, for the metric reader
    # A Place's throughput is the rate at which TOKENS leave it, so each
    # consuming mode contributes its firing rate times the multiplicity of the
    # arc it takes them through. That is SolverCTMC's convention and the one
    # Little's law needs. SolverSSA's NRM reports the OTHER one -- unweighted --
    # so the two disagree wherever an input arc has multiplicity above one.
    consumers = {}
    consumer_w = {}
    producers = {}
    for e in range(nev):
        if ev_kind[e] == 3:
            producers.setdefault((ev_station[e], ev_class[e]), []).append(e)
            continue
        if ev_kind[e] != 1 and ev_kind[e] != 4:
            continue
        md = modes[ev_mode[e]]
        for a, s in enumerate(md.arc_slot):
            key = (coord_station[s], coord_class[s])
            consumers.setdefault(key, []).append(e)
            consumer_w.setdefault(key, []).append(md.arc_w[a])

    # ---- the initial state
    x0 = np.zeros(nstate)
    for s in range(nm):
        x0[s] = m0full[coord_node[s], coord_class[s]]
    for md in modes:
        if md.nph > 1:
            e = np.inf
            for a, s in enumerate(md.arc_slot):
                e = min(e, x0[s] / md.arc_w[a])
            if not md.arc_slot:
                e = 1.0
            x0[md.zblk] = min(e, md.c) * md.pie

    t = Terms()
    t.M, t.K, t.I = M, K, I
    t.places, t.transitions = places, transitions
    t.names_node = names
    t.nstate, t.nm = nstate, nm
    t.pidx = pidx
    t.coord_node = np.asarray(coord_node, dtype=int)
    t.coord_class = np.asarray(coord_class, dtype=int)
    t.coord_station = np.asarray(coord_station, dtype=int)
    t.modes, t.timed_idx, t.imm_idx = modes, timed_idx, imm_idx
    t.D, t.rate_base, t.nev = D, rate_base, nev
    t.ev_kind, t.ev_mode, t.ev_phase = ev_kind, ev_mode, ev_phase
    t.ev_to, t.ev_station, t.ev_class = ev_to, ev_station, ev_class
    t.imm_col = np.flatnonzero(ev_kind == 4)
    t.latch_col = np.flatnonzero(ev_kind == 5)
    # The DIFFUSION counts the stochastic events only: an immediate flow and a
    # server latch are both the limit of an infinitely fast mechanism whose
    # fluctuation is slaved, not a Poisson stream with an intensity.
    t.stoch_col = np.flatnonzero((ev_kind != 4) & (ev_kind != 5))
    t.latch_mode = [j for j in timed_idx if modes[j].nph > 1]
    # THE COVARIANCE COVERS EVERY COORDINATE, phases included: the rate of a
    # multi-phase mode is linear in y and reads no marking, so dropping the
    # phases would sever that mode's whole restoring force.
    t.cov_idx = np.arange(nstate)
    t.cov_pairs, t.npair = cov_pairs, npair
    t.pair_index = -np.ones((max(nm, 1), max(nm, 1)), dtype=int)
    for i, (a, b) in enumerate(cov_pairs):
        t.pair_index[a, b] = i
        t.pair_index[b, a] = i
    t.imm_col_of = -np.ones(len(modes), dtype=int)
    for c in t.imm_col:
        t.imm_col_of[ev_mode[c]] = c
    t.consumers, t.consumer_w, t.producers = consumers, consumer_w, producers
    t.x0, t.m0full = x0, m0full
    t.options = options
    return t


def _is_exponential(sn, ist, r):
    """Whether the arrival process at (station, class) is exponential.

    Compared BY NAME, never by the raw integer: the ProcessType enums do not
    share numeric values across the codebases (MATLAB starts at EXP=0, native
    Python at EXP=1), so an integer comparison is only ever right by accident.
    """
    procid = getattr(sn, 'procid', None)
    if procid is None:
        return True
    v = procid[ist, r]
    name = getattr(v, 'name', None)
    if name is not None:
        return str(name).upper() == 'EXP'
    from ....constants import ProcessType
    for member in ProcessType:
        if str(member.name).upper() == 'EXP':
            return int(v) == int(member.value)
    return True


# -------------------------------------------------------------- the closures
def _sig(t, s2, a, b):
    """One entry of the closure covariance, zero where the drift never reads it."""
    if a >= t.pair_index.shape[0] or b >= t.pair_index.shape[1]:
        return 0.0
    i = t.pair_index[a, b]
    return 0.0 if i < 0 else float(s2[i])


class ThetaPack(object):
    __slots__ = ('theta', 'dslot', 'dval', 'dep', 'depslot', 'depval')


def theta(t, x, s2=None):
    """The closed enabling term of every mode, and its derivative.

        theta_j = ( prod_b Phi((thr_b - m_b)/sd_b) ) * E[ min_a(m_a/w_a), c_j ]

    Both factors collapse to their first-order form at zero variance, so the
    mean-field limit is one code path rather than two. THE VARIANCES ARE
    UNKNOWNS, not functions of x: the derivative is with respect to the MEANS.
    """
    nmod = len(t.modes)
    th = ThetaPack()
    th.theta = np.zeros(nmod)
    th.dslot = [np.zeros(0, dtype=int)] * nmod
    th.dval = [np.zeros(0)] * nmod
    th.dep = np.ones(nmod)
    th.depslot = [np.zeros(0, dtype=int)] * nmod
    th.depval = [np.zeros(0)] * nmod
    if s2 is None or np.size(s2) == 0:
        s2 = np.zeros(max(t.npair, 1))

    for j, md in enumerate(t.modes):
        A = len(md.arc_slot)
        mu = np.zeros(A)
        for a in range(A):
            mu[a] = x[md.arc_slot[a]] / md.arc_w[a]
        Sarg = np.zeros((A, A))
        if md.closable:
            for a in range(A):
                for b in range(A):
                    Sarg[a, b] = (_sig(t, s2, md.arc_slot[a], md.arc_slot[b])
                                  / (md.arc_w[a] * md.arc_w[b]))
        hmin, gmin, _ = minmulti_closure(mu, Sarg, md.c)

        nb = len(md.inh_slot)
        gate = np.ones(nb)
        dgate = np.zeros(nb)
        for b in range(nb):
            mb = x[md.inh_slot[b]]
            thr = md.inh_thr[b]
            vb = _sig(t, s2, md.inh_slot[b], md.inh_slot[b])
            if vb > 0:
                sd = np.sqrt(vb)
                zb = (thr - mb) / sd
                gate[b] = 0.5 * _erfc(-zb / np.sqrt(2.0))
                dgate[b] = -np.exp(-0.5 * zb * zb) / (np.sqrt(2.0 * np.pi) * sd)
            else:
                gate[b] = 1.0 if mb < thr - FINE_TOL * max(1.0, thr) else 0.0
                dgate[b] = 0.0
        ginh = float(np.prod(gate)) if nb else 1.0

        slots, vals = [], []
        for a in range(A):
            slots.append(md.arc_slot[a])
            vals.append(ginh * gmin[a] / md.arc_w[a])
        for b in range(nb):
            if gate[b] != 0:
                others = ginh / gate[b]
            else:
                others = float(np.prod(np.delete(gate, b))) if nb > 1 else 1.0
            slots.append(md.inh_slot[b])
            vals.append(hmin * others * dgate[b])
        # an arc and an inhibitor arc may share a coordinate, so accumulate
        if slots:
            uslots, inv = np.unique(np.asarray(slots, dtype=int), return_inverse=True)
            uvals = np.zeros(uslots.size)
            np.add.at(uvals, inv, np.asarray(vals, dtype=float))
        else:
            uslots, uvals = np.zeros(0, dtype=int), np.zeros(0)

        th.theta[j] = ginh * hmin
        th.dslot[j] = uslots
        th.dval[j] = uvals

        if md.dep is not None:
            g, gslot, gval = _dep(t, md, x)
            th.dep[j] = g
            th.depslot[j] = gslot
            th.depval[j] = gval
    return th


def _dep(t, md, x):
    """The marking-dependent firing multiplier and its gradient.

    g is a user function of the (nnodes x nclasses) marking, so it is evaluated
    at the MEAN marking -- a first-order closure of g, the same order at which
    SolverCTMC evaluates it per state and the only one available without the
    distribution of the marking. Its gradient has no analytic form, so it is
    taken by central differences.
    """
    mm = np.zeros((t.I, t.K))
    for s in range(t.nm):
        mm[t.coord_node[s], t.coord_class[s]] = x[s]
    g = float(md.dep(mm))
    slots, vals = [], []
    for s in range(t.nm):
        h = max(1e-6 * abs(x[s]), 1e-6)
        i, k = t.coord_node[s], t.coord_class[s]
        mp = mm.copy(); mp[i, k] = mm[i, k] + h
        mn = mm.copy(); mn[i, k] = max(0.0, mm[i, k] - h)
        hh = mp[i, k] - mn[i, k]
        if hh <= 0:
            continue
        d = (float(md.dep(mp)) - float(md.dep(mn))) / hh
        if d != 0:
            slots.append(s)
            vals.append(d)
    return g, np.asarray(slots, dtype=int), np.asarray(vals, dtype=float)


def rates(t, x, s2=None, phi=None, mu=None, th=None):
    """The rate of every event column of the fluid Petri drift.

    kind 1  firing of mode j, phase h -> h'   single phase: rateBase*theta*dep
                                              multi phase : rateBase*y(j,h)
    kind 2  internal phase change            rateBase*y(j,h)
    kind 3  exogenous arrival                a constant
    kind 4  firing of an IMMEDIATE mode      phi_j, an algebraic unknown
    kind 5  the server latch                 mu_j, a free-sign unknown
    """
    if th is None:
        th = theta(t, x, s2)
    if phi is None:
        phi = np.zeros(len(t.imm_idx))
    if mu is None:
        mu = np.zeros(len(t.latch_mode))

    imm_pos = -np.ones(len(t.modes), dtype=int)
    for i, j in enumerate(t.imm_idx):
        imm_pos[j] = i
    latch_pos = -np.ones(len(t.modes), dtype=int)
    for i, j in enumerate(t.latch_mode):
        latch_pos[j] = i

    r = np.zeros(t.nev)
    for e in range(t.nev):
        k = t.ev_kind[e]
        if k == 3:
            r[e] = t.rate_base[e]
        elif k == 4:
            r[e] = phi[imm_pos[t.ev_mode[e]]]
        elif k == 5:
            r[e] = mu[latch_pos[t.ev_mode[e]]]
        else:
            j = t.ev_mode[e]
            md = t.modes[j]
            if md.nph == 1:
                r[e] = t.rate_base[e] * th.theta[j] * th.dep[j]
            else:
                r[e] = t.rate_base[e] * x[md.zblk[t.ev_phase[e]]]
    return r, th


def jacobian(t, x, s2=None, phi=None, th=None):
    """Drift Jacobian A = D * dR/dX.

    This is what the Lyapunov equation of the linear noise approximation is
    written about, so it has to be the derivative of the SAME rate vector
    rates() returns: a covariance solved about an inconsistent Jacobian is not
    the covariance of anything. THE VARIANCES ARE HELD.
    """
    if th is None:
        th = theta(t, x, s2)
    n = t.nstate
    Jr = np.zeros((t.nev, n))
    for e in range(t.nev):
        k = t.ev_kind[e]
        if k in (3, 4, 5):
            continue
        j = t.ev_mode[e]
        md = t.modes[j]
        base = t.rate_base[e]
        if md.nph == 1:
            if th.dslot[j].size:
                Jr[e, th.dslot[j]] += base * th.dep[j] * th.dval[j]
            if th.depslot[j].size:
                Jr[e, th.depslot[j]] += base * th.theta[j] * th.depval[j]
        else:
            Jr[e, md.zblk[t.ev_phase[e]]] += base
    return t.D @ Jr, Jr


# ------------------------------------------------------------- conservation
def _rref(A, tol=1e-12):
    """Reduced row echelon form and the pivot columns."""
    A = np.array(A, dtype=float)
    rows, cols = A.shape
    piv = []
    r = 0
    for c in range(cols):
        if r >= rows:
            break
        k = r + int(np.argmax(np.abs(A[r:, c])))
        if abs(A[k, c]) <= tol:
            A[r:, c] = 0.0
            continue
        A[[r, k]] = A[[k, r]]
        A[r] = A[r] / A[r, c]
        for i in range(rows):
            if i != r and A[i, c] != 0.0:
                A[i] = A[i] - A[i, c] * A[r]
        piv.append(c)
        r += 1
    return A, piv


def _null_rational(A):
    """A rational basis of the null space, as MATLAB's null(A,'r') returns.

    D is integral -- arc multiplicities and unit phase moves -- so RREF gives
    exact rational rows, which keeps each conservation row readable as a
    statement about named places instead of an arbitrary orthogonal mixture.
    """
    A = np.asarray(A, dtype=float)
    if A.size == 0:
        return np.zeros((A.shape[1], 0))
    R, piv = _rref(A)
    cols = A.shape[1]
    free = [c for c in range(cols) if c not in piv]
    Z = np.zeros((cols, len(free)))
    for i, f in enumerate(free):
        Z[f, i] = 1.0
        for r, p in enumerate(piv):
            Z[p, i] = -R[r, f]
    return Z


class Conservation(object):
    __slots__ = ('C', 'N', 'leak', 'label')


def conservation(t):
    """The conserved quantities, as equations: u'D = 0 => u'x is constant.

    On the marking coordinates those u are the net's P-invariants; on a mode's
    phase block the all-ones vector is one of them, which is the statement that
    the phase coordinates are a distribution. Both come out of the same null
    space, so the phase normalisation needs no separate row. AN OPEN NET LOSES
    THE ROWS ITS ARRIVALS BREAK, automatically: the arrival columns are part of
    D, so a u an arrival moves is not in the null space.
    """
    cons = Conservation()
    if t.D.size == 0 or t.D.shape[1] == 0:
        cons.C = np.zeros((0, t.nstate))
        cons.N = np.zeros(0)
        cons.leak = 0.0
        cons.label = []
        return cons

    C = _null_rational(t.D.T).T
    C[np.abs(C) < 1e-12] = 0.0
    for c in range(C.shape[0]):
        nz = np.abs(C[c][C[c] != 0])
        if nz.size:
            C[c] = C[c] / nz.min()
    C = C[np.any(C != 0, axis=1)] if C.size else C.reshape(0, t.nstate)

    cons.C = C
    cons.N = C @ t.x0 if C.size else np.zeros(0)
    cons.leak = float(np.max(np.abs(C @ t.D))) if C.size else 0.0
    cons.label = [_row_label(t, C[c]) for c in range(C.shape[0])]
    return cons


def _row_label(t, row):
    parts = []
    for s in np.flatnonzero(row != 0):
        w = row[s]
        if s < t.nm:
            nm = '%s(class %d)' % (t.names_node[t.coord_node[s]], t.coord_class[s] + 1)
        else:
            nm = 'phase'
            for md in t.modes:
                if md.zblk and s in md.zblk:
                    nm = '%s phase %d' % (md.label, md.zblk.index(s) + 1)
                    break
        parts.append(nm if w == 1 else '%g*%s' % (w, nm))
    return ' + '.join(parts)


# --------------------------------------------------------------- capacities
class Constraints(object):
    __slots__ = ('A', 'b', 'label', 'cover')


def constraints(sn, t):
    """Every finite place capacity as a linear row A x <= b.

    THE GATE IS A LOSS ON THE DEPOSIT: LINE loses the tokens a firing would push
    past a place's capacity, so the fluid analogue scales the DEPOSIT leg of
    every event adding mass to the capped place and leaves the removal leg
    alone.
    """
    rows, b, label = [], [], []
    cap = np.asarray(sn.cap, dtype=float).ravel() if getattr(sn, 'cap', None) is not None else np.zeros(0)
    classcap = np.asarray(sn.classcap, dtype=float) if getattr(sn, 'classcap', None) is not None else np.zeros((0, 0))

    for ind in t.places:
        ist = int(sn.nodeToStation[ind])
        slots = [int(t.pidx[ind, k]) for k in range(t.K) if t.pidx[ind, k] >= 0]
        if not slots:
            continue
        if ist < cap.size and np.isfinite(cap[ist]):
            row = np.zeros(t.nstate)
            row[slots] = 1.0
            rows.append(row)
            b.append(float(cap[ist]))
            label.append('capacity %g of place %s' % (cap[ist], t.names_node[ind]))
        for k in range(t.K):
            s = int(t.pidx[ind, k])
            if s < 0 or ist >= classcap.shape[0] or not np.isfinite(classcap[ist, k]):
                continue
            row = np.zeros(t.nstate)
            row[s] = 1.0
            rows.append(row)
            b.append(float(classcap[ist, k]))
            label.append('class-%d capacity %g of place %s' % (k + 1, classcap[ist, k], t.names_node[ind]))

    # TWO ROWS THAT SAY THE SAME THING ARE A SINGULAR NEWTON SYSTEM, not a
    # redundancy the least squares absorbs, so an exact duplicate is pruned and
    # the tighter bound survives.
    keep = [True] * len(rows)
    for c in range(len(rows)):
        if not keep[c]:
            continue
        for d in range(c + 1, len(rows)):
            if keep[d] and np.array_equal(rows[c], rows[d]):
                if b[d] < b[c]:
                    keep[c] = False
                    break
                keep[d] = False

    con = Constraints()
    con.A = np.array([rows[i] for i in range(len(rows)) if keep[i]]) if any(keep) else np.zeros((0, t.nstate))
    con.b = np.array([b[i] for i in range(len(b)) if keep[i]]) if any(keep) else np.zeros(0)
    con.label = [label[i] for i in range(len(label)) if keep[i]]
    con.cover = con.A > 0
    return con


# ------------------------------------------------------------ the immediates
class Immediate(object):
    __slots__ = ('n', 'active', 'bind', 'pins', 'rows')


def _inhibited(md, x):
    """True when an inhibitor arc of this mode has reached its threshold.

    A HARD TEST ON THE MEAN, AND A KNOWN WRONG ANSWER WHEN THE MEAN SITS ON THE
    THRESHOLD -- a timed mode closes the same indicator as Phi((thr-m)/sd) in
    theta(); this path does not. See _kb/06-solver-catalog.md for the
    measurement and for what a real fix costs.
    """
    for b in range(len(md.inh_slot)):
        if x[md.inh_slot[b]] >= md.inh_thr[b]:
            return True
    return False


def immediate(t, x, imm=None):
    """The active set of the IMMEDIATE transitions, and the equations that pin
    their flows.

    An immediate transition has no rate: its fluid limit is a FLOW, an algebraic
    unknown pinned by the constraint that its binding input place holds no mass,

        phi_j >= 0,   x_b = 0 for the coordinate b that binds mode j

    with the GSPN conflict rule supplying the extra equation when two modes
    drain one place: phi_j*weight_l = phi_l*weight_j among the enabled modes of
    highest firing priority, and phi = 0 below it. THE COUNT IS SQUARE BY
    CONSTRUCTION: V pins plus (F-V) ratio rows is F equations for F flows.
    """
    n = len(t.imm_idx)
    if imm is None:
        imm = Immediate()
        imm.n = n
        imm.active = np.ones(n, dtype=bool)
        imm.bind = -np.ones(n, dtype=int)
        imm.pins = np.zeros(0, dtype=int)
        imm.rows = []
        # An inhibited mode never fires, so it neither carries a flow nor
        # empties a place. An empty input place is NOT a reason to deactivate --
        # that is the normal state of an enabled immediate mode.
        for k in range(n):
            md = t.modes[t.imm_idx[k]]
            if not md.arc_slot:
                raise ValueError(
                    "Immediate mode %s has no enabling arc, so nothing bounds its firing flow and the "
                    "net has no fluid limit. Give it an input place, or make it timed." % md.label)
            if _inhibited(md, x):
                imm.active[k] = False

    # ---- the assignment: each active mode binds the input arc it is shortest of
    for k in range(n):
        if not imm.active[k]:
            imm.bind[k] = -1
            continue
        md = t.modes[t.imm_idx[k]]
        if imm.bind[k] >= 0 and imm.bind[k] in md.arc_slot:
            continue  # a binding the caller set explicitly is kept
        lev = np.array([x[s] / w for s, w in zip(md.arc_slot, md.arc_w)])
        imm.bind[k] = md.arc_slot[int(np.argmin(lev))]

    # ---- the equations
    sel = imm.active & (imm.bind >= 0)
    imm.pins = np.unique(imm.bind[sel]) if np.any(sel) else np.zeros(0, dtype=int)
    rows = []
    for p in imm.pins:
        rows.append({'kind': 'pin', 'a': int(p), 'b': 0, 'wa': 0.0, 'wb': 0.0})
        grp = [k for k in range(n) if imm.active[k] and imm.bind[k] == p]
        if len(grp) <= 1:
            continue
        prio = np.array([t.modes[t.imm_idx[k]].prio for k in grp], dtype=float)
        wgt = np.array([t.modes[t.imm_idx[k]].weight for k in grp], dtype=float)
        top = [grp[i] for i in range(len(grp)) if prio[i] == prio.max()]
        low = [grp[i] for i in range(len(grp)) if prio[i] < prio.max()]
        tw = wgt[prio == prio.max()]
        for i in range(1, len(top)):
            rows.append({'kind': 'ratio', 'a': top[0], 'b': top[i],
                         'wa': float(tw[0]), 'wb': float(tw[i])})
        for k in low:
            rows.append({'kind': 'zero', 'a': int(k), 'b': 0, 'wa': 0.0, 'wb': 0.0})
    for k in np.flatnonzero(~imm.active):
        rows.append({'kind': 'zero', 'a': int(k), 'b': 0, 'wa': 0.0, 'wb': 0.0})
    imm.rows = rows
    return imm


# -------------------------------------------------------------- the refusals
def applicable(sn, options=None):
    """Whether the fluid Petri route can answer this model, and why not.

    A QUEUEING STATION IS THE ONE STRUCTURAL EXCLUSION: a net whose tokens also
    visit a Queue or a Delay is two formalisms at once, and LINE has no
    reference semantics for the hand-off.
    """
    I = int(len(sn.nodetype))
    nt = [int(sn.nodetype[i]) for i in range(I)]
    if int(NodeType.TRANSITION) not in nt:
        return False, 'the model has no Transition node, so it is not a Petri net'

    allowed = {int(NodeType.PLACE), int(NodeType.TRANSITION),
               int(NodeType.SOURCE), int(NodeType.SINK)}
    for ind in range(I):
        if nt[ind] not in allowed:
            return False, (
                "node %s is a %s. The fluid Petri route solves the marking of a Petri net, and a model "
                "that also holds queueing stations is two formalisms at once with no reference semantics "
                "for the hand-off; use SolverCTMC, SolverJMT, SolverSSA or SolverLDES"
                % (sn.nodenames[ind], _nodetype_text(sn.nodetype[ind])))

    # A queueing place declares a service process, which is what turns it into a
    # station with an embedded queue and a depository.
    rates_m = np.asarray(sn.rates, dtype=float)
    for ind in range(I):
        if nt[ind] != int(NodeType.PLACE):
            continue
        ist = int(sn.nodeToStation[ind])
        for k in range(int(sn.nclasses)):
            if ist < rates_m.shape[0] and not np.isnan(rates_m[ist, k]) and rates_m[ist, k] > 0:
                return False, (
                    "place %s is a QUEUEING place (it declares a service process), whose embedded queue "
                    "this drift does not carry; use SolverLDES" % sn.nodenames[ind])

    cfg = getattr(options, 'config', None) if options is not None else None
    if cfg is not None and getattr(cfg, 'hide_immediate', False):
        return False, (
            'options.config.hide_immediate eliminates the immediate transitions from the event set, but '
            'the fluid Petri route needs them: it solves their firing flows as algebraic unknowns')
    return True, ''


def _nodetype_text(v):
    name = getattr(v, 'name', None)
    if name is not None:
        return str(name).capitalize()
    for member in NodeType:
        if int(member) == int(v):
            return str(member.name).capitalize()
    return str(v)


# ----------------------------------------------------------- the DAE solver
class _Ctx(object):
    __slots__ = ('terms', 'imm', 'active', 'con', 'C', 'N', 'Dp', 'Dn')


def _context(terms, cons, con, imm, active):
    """Everything the residual needs that does not change within one active set.

    THE CONSERVATION ROWS A BINDING CAP BREAKS ARE DROPPED: a capped place loses
    the tokens that do not fit, so a conserved quantity supported on it is not
    conserved while the cap binds, and keeping its row would state an equation
    the drift contradicts -- a singular Newton system rather than an inaccuracy.
    """
    ctx = _Ctx()
    ctx.terms = terms
    ctx.imm = imm
    ctx.active = np.asarray(active, dtype=int).ravel()
    ctx.con = con
    C, N = cons.C, cons.N
    if ctx.active.size and C.size:
        hit = np.any(con.cover[ctx.active, :], axis=0)
        drop = np.any(C[:, hit] != 0, axis=1)
        C, N = C[~drop], N[~drop]
    ctx.C, ctx.N = C, N
    ctx.Dp = np.maximum(terms.D, 0.0)
    ctx.Dn = np.minimum(terms.D, 0.0)
    return ctx


def _unpack(u, terms, imm, active):
    n, npair, ni = terms.nstate, terms.npair, imm.n
    nl = len(terms.latch_mode)
    na = len(active)
    x = u[:n]
    s2 = u[n:n + npair]
    phi = np.maximum(0.0, u[n + npair:n + npair + ni])
    mu = u[n + npair + ni:n + npair + ni + nl]
    zeta = np.maximum(0.0, u[n + npair + ni + nl:n + npair + ni + nl + na])
    return x, s2, phi, mu, zeta


def clamp_tangent(terms, imm, con, active, th=None):
    """The reduction of the fluctuation onto the manifold the fast and clamped
    directions leave free.

    AN IMMEDIATE PIN REDUCES OBLIQUELY, ALONG THE FAST REACTION ITSELF, and this
    is the one place where the orthogonal projector the queueing twin uses is
    WRONG rather than merely different: a slow event depositing into a pinned
    place is answered instantly by the immediate transition, so its effective
    jump is its own plus the immediate flow it triggers -- the token is
    forwarded, not lost. An orthogonal projection deletes the deposit and
    destroys mass in the diffusion.

        P = I - Cf * G * E_B,   G = G0 * (E_B Cf G0)^-1

    A CAPACITY CAP REDUCES ORTHOGONALLY -- mass that does not fit is genuinely
    lost, so there is nothing to forward it to. A SERVER LATCH REDUCES
    ORTHOGONALLY TOO, on the LINEARISED row [-dtheta_j/dm, 1 over the phase
    block], which is why this projector depends on the iterate.
    """
    idx = terms.cov_idx
    nc = idx.size
    T = np.eye(nc)

    actk = np.flatnonzero(imm.active & (imm.bind >= 0))
    B = np.asarray(imm.pins, dtype=int).ravel()
    if actk.size and B.size:
        Cf = np.zeros((nc, actk.size))
        for a, k in enumerate(actk):
            Cf[:, a] = terms.modes[terms.imm_idx[k]].cvec[idx]
        G0 = np.zeros((actk.size, B.size))
        for jb, b in enumerate(B):
            grp = np.flatnonzero(imm.bind[actk] == b)
            w = np.array([terms.modes[terms.imm_idx[actk[q]]].weight for q in grp], dtype=float)
            if w.sum() <= 0:
                w = np.ones(grp.size)
            G0[grp, jb] = w / w.sum()
        EB = np.zeros((B.size, nc))
        for jb, b in enumerate(B):
            EB[jb, b] = 1.0
        Mb = (EB @ Cf) @ G0
        T = T - Cf @ (G0 @ np.linalg.pinv(Mb)) @ EB

    R = []
    for c in np.asarray(active, dtype=int).ravel():
        R.append(con.A[c, idx])
    if th is not None:
        for j in terms.latch_mode:
            row = np.zeros(nc)
            row[terms.modes[j].zblk] = 1.0
            if th.dslot[j].size:
                row[th.dslot[j]] -= th.dval[j]
            R.append(row)
    if R:
        R = np.array(R)
        if np.any(np.abs(R) > 1e-14):
            T = (np.eye(nc) - R.T @ np.linalg.pinv(R @ R.T) @ R) @ T
    if np.max(np.abs(T - np.eye(nc))) <= 1e-14:
        return None
    return T


def sigma_of(terms, A, r, clampT):
    """One Lyapunov solve, over the marking coordinates.

    THE DIFFUSION COUNTS THE STOCHASTIC EVENTS ONLY: an immediate flow is not a
    Poisson stream with an intensity but the limit of an infinitely fast one
    whose fluctuation is slaved, and its pinned coordinate is projected out.
    """
    from .minnormal import fluid_lyapunov
    idx = terms.cov_idx
    Dc = terms.D[np.ix_(idx, terms.stoch_col)]
    Am = A[np.ix_(idx, idx)]
    if clampT is not None:
        # BOTH the jump directions and the generator are reduced: reducing Dc
        # alone would fix the subspace but leave the generator's orthogonal
        # component on it, which is not the reduced dynamics when the reduction
        # is oblique.
        Dc = clampT @ Dc
        Am = clampT @ Am
    Q = Dc @ np.diag(r[terms.stoch_col]) @ Dc.T
    Sc = fluid_lyapunov(Am, Q, Dc)
    Sigma = np.zeros((terms.nstate, terms.nstate))
    Sigma[np.ix_(idx, idx)] = Sc
    return Sigma


def _residual(u, ctx, quiet=True):
    """The coupled algebraic system, stacked.

    Returns None when the closure cannot be evaluated at this iterate, so the
    line search can back off; the first evaluation of a pass runs with
    quiet=False, where a genuine failure surfaces.
    """
    terms, imm = ctx.terms, ctx.imm
    n = terms.nstate
    x, s2, phi, mu, zeta = _unpack(u, terms, imm, ctx.active)

    th = theta(terms, x, s2)
    r, _ = rates(terms, x, s2, phi, mu, th)

    # the deposit gate of every binding capacity, as a product of fractions
    gain = np.ones(n)
    for k, c in enumerate(ctx.active):
        gain[ctx.con.cover[c, :]] *= zeta[k]
    drift = ctx.Dn @ r + gain * (ctx.Dp @ r)

    try:
        A, _ = jacobian(terms, x, s2, phi, th)
        Sigma = sigma_of(terms, A, r, clamp_tangent(terms, imm, ctx.con, ctx.active, th))
    except Exception:
        if not quiet:
            raise
        return None, None, th, None

    s2new = np.array([Sigma[a, b] for (a, b) in terms.cov_pairs]) if terms.npair else np.zeros(0)

    G = [drift]
    if ctx.C.size:
        G.append(ctx.C @ x - ctx.N)
    G.append(s2 - s2new)
    extra = []
    for row in imm.rows:
        if row['kind'] == 'pin':
            extra.append(x[row['a']])
        elif row['kind'] == 'ratio':
            extra.append(phi[row['a']] * row['wb'] - phi[row['b']] * row['wa'])
        else:
            extra.append(phi[row['a']])
    for j in terms.latch_mode:
        extra.append(float(np.sum(x[terms.modes[j].zblk])) - th.theta[j])
    for k, c in enumerate(ctx.active):
        extra.append(float(ctx.con.A[c, :] @ x) - ctx.con.b[c])
    if extra:
        G.append(np.asarray(extra, dtype=float))
    return np.concatenate(G), r, th, Sigma


def _project(u, lb):
    """An iterate projected onto its feasible box, the lower bound only.

    The bound is a VECTOR, one entry per unknown, never a count: a scalar
    'nfree' is indistinguishable from a one-unknown bound vector, which is how
    the MATLAB twin crashed on the simplest net in the tree.
    """
    if lb is None or np.size(lb) == 0:
        return u
    lb = np.asarray(lb, dtype=float).ravel()
    if lb.size != u.size:
        raise ValueError('The bound vector has %d entries for %d unknowns.' % (lb.size, u.size))
    fin = np.isfinite(lb)
    u = u.copy()
    u[fin] = np.maximum(lb[fin], u[fin])
    return u


def _fdjac(fun, u, G):
    """Forward-difference Jacobian of the residual."""
    n = u.size
    m = G.size
    J = np.zeros((m, n))
    for k in range(n):
        h = 1e-7 * max(1.0, abs(u[k]))
        up = u.copy()
        up[k] += h
        Gp = fun(up, True)[0]
        if Gp is None:
            up[k] = u[k] - h
            Gp = fun(up, True)[0]
            if Gp is None:
                continue
            J[:, k] = (G - Gp) / h
        else:
            J[:, k] = (Gp - G) / h
    return J


def newton(fun, u0, tol, maxit, lb):
    """Damped projected Newton with a finite-difference Jacobian."""
    u = _project(np.asarray(u0, dtype=float).copy(), lb)
    G = fun(u, False)[0]
    if G is None:
        return u, 0, False, np.inf
    resnorm = float(np.linalg.norm(G, np.inf))
    it = 0
    for it in range(1, maxit + 1):
        if resnorm <= tol:
            return u, it - 1, True, resnorm
        J = _fdjac(fun, u, G)
        try:
            du = np.linalg.lstsq(J, -G, rcond=None)[0]
        except np.linalg.LinAlgError:
            break
        lam = 1.0
        improved = False
        for _ in range(30):
            un = _project(u + lam * du, lb)
            Gn = fun(un, True)[0]
            if Gn is not None:
                rn = float(np.linalg.norm(Gn, np.inf))
                if rn < resnorm:
                    u, G, resnorm = un, Gn, rn
                    improved = True
                    break
            lam *= 0.5
        if not improved:
            break
    return u, it, resnorm <= tol, resnorm


# ------------------------------------------------------------------ the seed
def _seed_drift(terms, x, lam):
    """The first-order drift: the same rates at zero variance, with an immediate
    mode firing at LAM times its enabling degree and its firing weight, and the
    server latch relaxed at LAM towards the enabling degree instead of solved.
    Both are approximations of an algebraic constraint by a fast reaction, and
    both are confined to the seed.
    """
    x = np.maximum(0.0, x)
    th = theta(terms, x, np.zeros(max(terms.npair, 1)))
    phi = np.zeros(len(terms.imm_idx))
    for k, j in enumerate(terms.imm_idx):
        phi[k] = lam * terms.modes[j].weight * th.theta[j]
    mu = np.zeros(len(terms.latch_mode))
    for q, j in enumerate(terms.latch_mode):
        mu[q] = lam * (th.theta[j] - float(np.sum(x[terms.modes[j].zblk])))
    r, _ = rates(terms, x, None, phi, mu, th)
    return terms.D @ r


def seed(terms, options=None):
    """One first-order trajectory from the initial marking.

    Newton needs a point in the basin, not an answer. The immediate modes get a
    large FINITE rate here and only here, scaled to the model's own timescale
    rather than taken from GlobalConstants.Immediate: 1e8 against a rate of
    order one is a stiffness the seed does not need, and the answer does not
    depend on the seed's accuracy.
    """
    from scipy.integrate import solve_ivp

    rmax = 0.0
    for j in terms.timed_idx:
        rmax = max(rmax, float(np.sum(terms.modes[j].d1)))
    for e in np.flatnonzero(terms.ev_kind == 3):
        rmax = max(rmax, float(terms.rate_base[e]))
    if not (rmax > 0):
        rmax = 1.0
    lam = min(1e8, 1e4 * rmax)

    mass = float(np.sum(terms.x0[:terms.nm]))
    T = 50.0 * (mass + 1.0) / rmax
    t = np.array([0.0])
    xt = terms.x0.reshape(1, -1)
    for _ in range(6):
        sol = solve_ivp(lambda tt, xx: _seed_drift(terms, xx, lam), [0.0, T], terms.x0,
                        method='LSODA', rtol=COARSE_TOL, atol=FINE_TOL, dense_output=False)
        t = sol.t
        xt = sol.y.T
        d = _seed_drift(terms, xt[-1], lam)
        if np.max(np.abs(d)) <= COARSE_TOL * max(1.0, rmax * (mass + 1.0)):
            break
        T *= 4.0
    return t, xt, lam


def seed_flows(terms, xt, lam):
    """The immediate flows the SEED trajectory carried, point by point -- the
    approximation the seed integrated, and therefore the one its throughput
    table has to be read at."""
    nt = xt.shape[0]
    ni = len(terms.imm_idx)
    phit = np.zeros((nt, ni))
    if ni == 0:
        return phit
    for a in range(nt):
        th = theta(terms, xt[a], np.zeros(max(terms.npair, 1)))
        for k, j in enumerate(terms.imm_idx):
            phit[a, k] = lam * terms.modes[j].weight * th.theta[j]
    return phit


def collapse(terms, imm, x0):
    """The fluid vanishing-marking collapse: move the initial marking along the
    immediate incidence columns until every pinned place is empty. Every
    conserved quantity survives it for free, since a conserved direction
    annihilates those columns."""
    x = x0.copy()
    act = np.flatnonzero(imm.active)
    pins = np.asarray(imm.pins, dtype=int).ravel()
    if act.size == 0 or pins.size == 0:
        return x
    Cimm = np.zeros((terms.nstate, act.size))
    for a, k in enumerate(act):
        Cimm[:, a] = terms.modes[terms.imm_idx[k]].cvec
    b = -x0[pins]
    s = np.linalg.lstsq(Cimm[pins, :], b, rcond=None)[0]
    if not np.all(np.isfinite(s)):
        return x
    s = np.maximum(0.0, s)
    x = x0 + Cimm @ s
    x[pins] = 0.0
    x[:terms.nm] = np.maximum(0.0, x[:terms.nm])
    # the collapse moved the marking, so the servers a multi-phase mode has
    # latched moved with it: re-read the enabling degree at the collapsed point
    th = theta(terms, x, np.zeros(max(terms.npair, 1)))
    for j in terms.latch_mode:
        x[terms.modes[j].zblk] = th.theta[j] * terms.modes[j].pie
    return x


# --------------------------------------------------------------- the metrics
def metrics(terms, x, r, M, K):
    """The station table, in the conventions the exact SPN engines report.

    A PLACE IS AN INF STATION: its queue length is its mean token count and its
    utilization is the same number. Its throughput is the rate at which TOKENS
    leave it, so a consuming mode contributes its firing rate times the arc
    multiplicity -- SolverCTMC's convention, and the one Little's law needs.
    SolverSSA's NRM sums the UNWEIGHTED propensity, so the two disagree wherever
    an input arc has multiplicity above one; the exact engine is the reference.
    """
    QN = np.zeros((M, K)); UN = np.zeros((M, K))
    RN = np.zeros((M, K)); TN = np.zeros((M, K))
    for s in range(terms.nm):
        ist, k = terms.coord_station[s], terms.coord_class[s]
        QN[ist, k] += x[s]
        UN[ist, k] = QN[ist, k]
    for (ist, k), e in terms.consumers.items():
        w = np.asarray(terms.consumer_w[(ist, k)], dtype=float)
        TN[ist, k] += float(np.sum(w * r[np.asarray(e, dtype=int)]))
    for (ist, k), e in terms.producers.items():
        TN[ist, k] += float(np.sum(r[np.asarray(e, dtype=int)]))
    # TN is zero only to the integrator's accuracy; see minnormal.py.
    nz = TN > GlobalConstants.Zero
    RN[nz] = QN[nz] / TN[nz]
    return QN, UN, RN, TN


def metrics_t(terms, xt, s2, phit, M, K):
    """The same reader along a trajectory."""
    nt = xt.shape[0]
    QNt = [[np.zeros(nt) for _ in range(K)] for _ in range(M)]
    UNt = [[np.zeros(nt) for _ in range(K)] for _ in range(M)]
    TNt = [[np.zeros(nt) for _ in range(K)] for _ in range(M)]
    Rt = np.zeros((nt, terms.nev))
    for a in range(nt):
        # the latch flow moves no marking, so it cannot change a place throughput
        Rt[a, :], _ = rates(terms, xt[a], s2, phit[a] if phit.size else None, None)
    for s in range(terms.nm):
        ist, k = terms.coord_station[s], terms.coord_class[s]
        QNt[ist][k] = QNt[ist][k] + xt[:, s]
        UNt[ist][k] = QNt[ist][k]
    for (ist, k), e in terms.consumers.items():
        w = np.asarray(terms.consumer_w[(ist, k)], dtype=float)
        TNt[ist][k] = TNt[ist][k] + Rt[:, np.asarray(e, dtype=int)] @ w
    for (ist, k), e in terms.producers.items():
        TNt[ist][k] = TNt[ist][k] + Rt[:, np.asarray(e, dtype=int)].sum(axis=1)
    return QNt, UNt, TNt


def qvar(terms, Sigma, M, K):
    """The variance of each station-class token count."""
    QVar = np.zeros((M, K))
    for s in range(terms.nm):
        QVar[terms.coord_station[s], terms.coord_class[s]] = max(0.0, Sigma[s, s])
    return QVar


def report(terms, cons, con, imm, active, x, r, phi, zeta, Sigma):
    """What the Petri route computes that the station table has no column for."""
    rep = {}
    marking = np.zeros((terms.I, terms.K))
    markingvar = np.zeros((terms.I, terms.K))
    for s in range(terms.nm):
        marking[terms.coord_node[s], terms.coord_class[s]] = x[s]
        markingvar[terms.coord_node[s], terms.coord_class[s]] = max(0.0, Sigma[s, s])
    rep['marking'] = marking
    rep['markingVar'] = markingvar
    rep['modeLabel'] = [md.label for md in terms.modes]
    flow = np.zeros(len(terms.modes))
    for j in range(len(terms.modes)):
        e = np.flatnonzero((terms.ev_mode == j) & ((terms.ev_kind == 1) | (terms.ev_kind == 4)))
        if e.size:
            flow[j] = float(np.sum(r[e]))
    rep['modeFlow'] = flow
    rep['immediateFlow'] = phi
    rep['invariantLabel'] = cons.label
    rep['invariantValue'] = cons.N
    rep['invariantError'] = (cons.C @ x - cons.N) if cons.C.size else np.zeros(0)
    rep['capacityLabel'] = con.label
    rep['capacityActive'] = np.asarray(active, dtype=int)
    rep['capacityFraction'] = np.asarray(zeta, dtype=float)
    rep['pinned'] = np.asarray(imm.pins, dtype=int)
    rep['Sigma'] = Sigma
    return rep


# ---------------------------------------------------------------- the driver
def solve_petri(sn, options=None):
    """Fluid analysis of a stochastic Petri net: one simultaneous algebraic
    solve per active set.

    Returns a dict with QN, UN, RN, TN, the trajectory, and a moments block
    carrying the marking covariance, the per-mode firing flows, the invariants
    and which capacities bound.
    """
    import time
    t0 = time.time()
    M, K = int(sn.nstations), int(sn.nclasses)

    ok, why = applicable(sn, options)
    if not ok:
        raise ValueError('The fluid Petri route cannot solve this model: %s.' % why)

    terms = build_terms(sn, options)
    n, npair = terms.nstate, terms.npair

    maxstate = 100
    cfg = getattr(options, 'config', None) if options is not None else None
    if cfg is not None and getattr(cfg, 'dae_maxstate', None):
        maxstate = int(cfg.dae_maxstate)
    if n > maxstate:
        raise ValueError(
            'The fluid Petri route solves a %d-unknown algebraic system with a finite-difference '
            'Jacobian, above the limit of %d set by options.config.dae_maxstate. Raise that limit, '
            'or use SolverSSA for a net of this size.' % (n, maxstate))

    cons = conservation(terms)
    if cons.leak > 1e-7:
        raise ValueError('The conserved directions and the jump matrix disagree: the largest leak per '
                         'unit rate is %g, where it must be zero.' % cons.leak)
    con = constraints(sn, terms)
    ncon = con.b.size

    # -- seed
    tseed, xseed, seed_lam = seed(terms, options)
    x = xseed[-1].copy()

    imm = immediate(terms, x)
    nimm = imm.n

    # THE VARIANCE IS SEEDED POSITIVE: sigma2 = 0 is where min() has no
    # derivative, and a saturated net's first-order fixed point sits there.
    s2 = np.zeros(npair)
    ondiag = np.array([a == b for (a, b) in terms.cov_pairs], dtype=bool) if npair else np.zeros(0, bool)
    for p in np.flatnonzero(ondiag):
        s2[p] = max(FINE_TOL, x[terms.cov_pairs[p][0]])

    # THE NEWTON TOLERANCE IS THE FINE ONE BY DEFAULT: the conservation rows are
    # LINEAR in the unknowns, so the residual norm IS the token-count error, and
    # stopping at 1e-4 would leave a closed net holding 4.0001 tokens.
    tol = FINE_TOL
    otol = getattr(options, 'tol', None) if options is not None else None
    if otol is not None and np.isfinite(otol) and otol < tol:
        tol = float(otol)
    newton_max = 50
    itmax = getattr(options, 'iter_max', None) if options is not None else None
    if itmax:
        newton_max = max(newton_max, int(itmax))

    # -- steady state
    active = []
    phi = np.zeros(nimm)
    nlatch = len(terms.latch_mode)
    mu = np.zeros(nlatch)
    zeta = np.zeros(0)
    iters = 0
    aset_max = max(6, 2 * (ncon + nimm) + 2)
    converged, resnorm, best = False, np.inf, None

    for _ in range(aset_max):
        ctx = _context(terms, cons, con, imm, active)
        u0 = np.concatenate([x, s2, phi, mu, np.ones(len(active))])
        lb = np.concatenate([np.full(n, -np.inf), np.full(npair, -np.inf),
                             np.zeros(nimm), np.full(nlatch, -np.inf), np.zeros(len(active))])
        if npair:
            lb[n + np.flatnonzero(ondiag)] = 0.0
        fun = lambda uu, quiet=True: _residual(uu, ctx, quiet)
        u, nit, converged, resnorm = newton(fun, u0, tol, newton_max, lb)
        iters += nit
        x, s2, phi, mu, zeta = _unpack(u, terms, imm, active)

        if best is None or resnorm < best['resnorm']:
            best = dict(x=x.copy(), s2=s2.copy(), phi=phi.copy(), mu=mu.copy(), zeta=zeta.copy(),
                        active=list(active), imm=imm, resnorm=resnorm, converged=converged)

        # the active-set moves, in the order that a failure of each invalidates
        # the next: a negative flow means the mode does not fire at all, a
        # negative marking means the wrong coordinate was pinned, and only then
        # is it worth asking which capacity rows bind.
        moved = False
        if nimm > 0:
            bad = np.flatnonzero(phi < -max(tol, 1e-10) * max(1.0, float(np.max(np.abs(phi))) if phi.size else 1.0))
            if bad.size:
                imm.active[bad[0]] = False
                imm = immediate(terms, x, imm)
                moved = True
        if not moved and nimm > 0:
            neg = np.flatnonzero(x[:terms.nm] < -max(tol, 1e-10))
            if neg.size:
                s = int(neg[0])
                cand = [k for k in range(nimm)
                        if imm.active[k] and s in terms.modes[terms.imm_idx[k]].arc_slot
                        and imm.bind[k] != s]
                if cand:
                    imm.bind[cand[0]] = s
                    imm = immediate(terms, x, imm)
                    moved = True
        if not moved and ncon > 0:
            val = con.A @ x
            over = [c for c in np.flatnonzero(val > con.b + max(1e-9, tol)) if c not in active]
            if over:
                active = list(active) + [int(c) for c in over]
                moved = True
            else:
                rel = np.flatnonzero(zeta > 1 + max(1e-9, tol))
                if rel.size:
                    active = [active[i] for i in range(len(active)) if i not in set(rel.tolist())]
                    moved = True
        if not moved:
            break

    if best is not None and not converged and best['converged']:
        x, s2, phi = best['x'], best['s2'], best['phi']
        mu, zeta, active, imm = best['mu'], best['zeta'], best['active'], best['imm']
        converged, resnorm = True, best['resnorm']

    warnings_out = []
    if not converged:
        warnings_out.append(
            'The simultaneous closure solve stopped at residual %.3e after %d Newton steps without '
            'reaching %.3e. The reported point is the last iterate.' % (resnorm, iters, tol))
    # A NEGATIVE MARKING IS NOT ROUNDING: the fixed point wanted mass a place
    # cannot supply and no immediate mode could be rebound to pin it, so the
    # answer is outside the model's own state space and is reported as such.
    neg = np.flatnonzero(x[:terms.nm] < -max(tol, 1e-10))
    if neg.size:
        warnings_out.append(
            'The fixed point holds %g tokens at %s, which is negative: no immediate transition could be '
            'rebound to pin that place at zero. Use SolverCTMC, SolverSSA or SolverLDES for this net.'
            % (x[neg[0]], terms.names_node[terms.coord_node[neg[0]]]))
    for w in warnings_out:
        import warnings as _w
        _w.warn(w, RuntimeWarning)

    ctx = _context(terms, cons, con, imm, active)
    _, r, _, Sigma = _residual(np.concatenate([x, s2, phi, mu, zeta]), ctx, False)

    # -- the trajectory
    # The immediate flow ALONG the reported path is not the steady-state one: on
    # the seed path it is the large-finite-rate approximation the seed itself
    # integrated, so the throughput table and the trajectory it is read off
    # describe the same run.
    t = tseed
    xvec_t = xseed
    phit = seed_flows(terms, xseed, seed_lam)
    Sigmat, tvar, switches = None, np.zeros(0), []
    tspan = list(getattr(options, 'timespan', [0.0, np.inf])) if options is not None else [0.0, np.inf]
    if len(tspan) > 1 and np.isfinite(tspan[1]):
        t, xvec_t, phit, Sigmat, tvar, switches = transient(
            terms, cons, con, imm, active, options, x, s2, phi, zeta)

    QN, UN, RN, TN = metrics(terms, x, r, M, K)
    QNt, UNt, TNt = metrics_t(terms, xvec_t, s2, phit, M, K)

    out = {
        'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN,
        't': t, 'xvec_t': xvec_t, 'QNt': QNt, 'UNt': UNt, 'TNt': TNt,
        'x': x, 'iters': iters, 'resnorm': resnorm, 'converged': converged,
        'runtime': time.time() - t0,
        'moments': {
            'Sigma': Sigma,
            'petri': report(terms, cons, con, imm, active, x, r, phi, zeta, Sigma),
            'QVar': qvar(terms, Sigma, M, K),
            'method': 'dae',
            'resnorm': resnorm,
            'Sigmat': Sigmat,
            'tvar': tvar,
            'switches': switches,
        },
        'terms': terms,
    }
    out['moments']['QStd'] = np.sqrt(np.maximum(0.0, out['moments']['QVar']))
    return out


class PetriSolver(object):
    """The 'dae' method's Petri arm, in the shape SolverFLD dispatches to.

    A model holding any Transition node is a different formalism, so it goes to
    a dedicated runner rather than through the queueing dispatch: the
    post-processing of the queueing path rewrites UN from sn.rates (NaN at a
    Place) and RN from the station scheduling, and would overwrite the Petri
    conventions this reports.
    """

    def __init__(self, sn, options):
        self.sn = sn
        self.options = options

    def solve(self):
        from ..options import FLDResult
        sn, options = self.sn, self.options
        out = solve_petri(sn, options)
        M, K = int(sn.nstations), int(sn.nclasses)
        QN, UN, RN, TN = out['QN'], out['UN'], out['RN'], out['TN']

        XN = np.zeros((1, K))
        CN = np.zeros((1, K))
        refstat = np.asarray(sn.refstat, dtype=int).ravel() if getattr(sn, 'refstat', None) is not None \
            else np.zeros(K, dtype=int)
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        for k in range(K):
            ist = int(refstat[k]) if k < refstat.size else -1
            if 0 <= ist < M:
                XN[0, k] = TN[ist, k]
                if XN[0, k] > 0 and k < njobs.size and np.isfinite(njobs[k]):
                    CN[0, k] = njobs[k] / XN[0, k]

        QNt, UNt, TNt = {}, {}, {}
        for ist in range(M):
            for k in range(K):
                QNt[(ist, k)] = out['QNt'][ist][k]
                UNt[(ist, k)] = out['UNt'][ist][k]
                TNt[(ist, k)] = out['TNt'][ist][k]

        res = FLDResult(QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN)
        res.t = out['t']
        res.QNt, res.UNt, res.TNt = QNt, UNt, TNt
        res.xvec = out['x']
        res.iterations = out['iters']
        res.runtime = out['runtime']
        res.method = 'dae'
        res.moments = out['moments']
        return res


# ------------------------------------------------------------- the transient
def _dae_rhs(z, ctx, withcov, nc):
    """The right-hand side of one segment, in the layout transient() documents."""
    terms, imm = ctx.terms, ctx.imm
    n, nimm = terms.nstate, imm.n
    nact = ctx.active.size
    nlatch = len(terms.latch_mode)
    x = z[:n]
    phi = z[n:n + nimm]
    mu = z[n + nimm:n + nimm + nlatch]
    zeta = z[n + nimm + nlatch:n + nimm + nlatch + nact]
    th = theta(terms, x, np.zeros(max(terms.npair, 1)))
    r, _ = rates(terms, x, None, phi, mu, th)
    gain = np.ones(n)
    for k, c in enumerate(ctx.active):
        gain[ctx.con.cover[c, :]] *= zeta[k]
    drift = ctx.Dn @ r + gain * (ctx.Dp @ r)

    dz = drift.copy()
    for p in np.asarray(imm.pins, dtype=int).ravel():
        dz[p] = x[p]  # the algebraic level pin
    extra = []
    for row in imm.rows:
        if row['kind'] == 'pin':
            extra.append(drift[row['a']])
        elif row['kind'] == 'ratio':
            extra.append(phi[row['a']] * row['wb'] - phi[row['b']] * row['wa'])
        else:
            extra.append(phi[row['a']])
    # The DIFFERENTIATED latch, one row per multi-phase mode: the running-server
    # count tracks the enabling degree exactly, so their rates of change agree
    # and that is the equation MU solves.
    for j in terms.latch_mode:
        acc = float(np.sum(drift[terms.modes[j].zblk]))
        if th.dslot[j].size:
            acc -= float(th.dval[j] @ drift[th.dslot[j]])
        extra.append(acc)
    for c in ctx.active:
        extra.append(float(ctx.con.A[c, :] @ drift))
    if extra:
        dz = np.concatenate([dz, np.asarray(extra, dtype=float)])

    if withcov:
        idx = terms.cov_idx
        A, _ = jacobian(terms, x, np.zeros(max(terms.npair, 1)), phi, th)
        Dc = terms.D[np.ix_(idx, terms.stoch_col)]
        Am = A[np.ix_(idx, idx)]
        clampT = clamp_tangent(terms, imm, ctx.con, ctx.active, th)
        if clampT is not None:
            Dc = clampT @ Dc
            Am = clampT @ Am
        Q = Dc @ np.diag(r[terms.stoch_col]) @ Dc.T
        S = z[n + nimm + nlatch + nact:n + nimm + nlatch + nact + nc * nc].reshape(nc, nc)
        dS = Am @ S + S @ Am.T + Q
        dz = np.concatenate([dz, dS.ravel()])
    return dz


def _event_values(z, ctx):
    """Where the segment ends: a flow that would go negative, a capacity that
    starts or stops binding, or a marking coordinate that would go negative."""
    terms, imm = ctx.terms, ctx.imm
    n, nimm = terms.nstate, imm.n
    nact = ctx.active.size
    nlatch = len(terms.latch_mode)
    x = z[:n]
    phi = z[n:n + nimm]
    zeta = z[n + nimm + nlatch:n + nimm + nlatch + nact]
    val = []
    for k in range(nimm):
        val.append(phi[k] if imm.active[k] else 1.0)
    for c in range(ctx.con.b.size):
        hit = np.flatnonzero(ctx.active == c)
        if hit.size == 0:
            val.append(float(ctx.con.b[c] - ctx.con.A[c, :] @ x))
        else:
            val.append(1.0 - zeta[hit[0]])
    pins = set(np.asarray(imm.pins, dtype=int).ravel().tolist())
    for s in range(terms.nm):
        val.append(1.0 if s in pins else x[s] + FINE_TOL)
    return np.asarray(val, dtype=float)


def _switch(terms, con, imm, active, ievent, x, zeta):
    """The active-set move a located crossing asks for, and the multiplier
    vector resized to match: a newly bound capacity starts unthrottled, a
    released one loses its unknown."""
    nimm = imm.n
    ncon = con.b.size
    active = list(np.asarray(active, dtype=int).ravel())
    zeta = list(np.asarray(zeta, dtype=float).ravel())
    if ievent < nimm:
        if imm.active[ievent]:
            imm.active[ievent] = False
            imm = immediate(terms, x, imm)
            return imm, active, True, np.asarray(zeta)
        return imm, active, False, np.asarray(zeta)
    ievent -= nimm
    if ievent < ncon:
        if ievent in active:
            k = active.index(ievent)
            active.pop(k)
            zeta.pop(k)
        else:
            active.append(int(ievent))
            zeta.append(1.0)
        return imm, active, True, np.asarray(zeta)
    s = ievent - ncon
    cand = [k for k in range(nimm)
            if imm.active[k] and s in terms.modes[terms.imm_idx[k]].arc_slot and imm.bind[k] != s]
    if cand:
        imm.bind[cand[0]] = s
        imm = immediate(terms, x, imm)
        return imm, active, True, np.asarray(zeta)
    return imm, active, False, np.asarray(zeta)


def transient(terms, cons, con, imm, active, options, xss, s2ss, phiss, zetass):
    """The transient closure as an index-1 DAE with a SINGULAR mass matrix.

    THE INITIAL MARKING JUMPS: an immediate transition fires in zero time, so a
    marking that enables one is not a state the trajectory ever occupies. The
    fluid analogue of the vanishing-marking collapse moves the initial marking
    along the immediate incidence columns until every pinned place is empty, and
    the trajectory starts THERE -- a jump that conserves every P-invariant by
    construction, since a conserved direction annihilates the columns it moves
    along.

    Per segment:
      d/dt x_s = drift_s                      a coordinate nothing pins
      0        = x_b                          a coordinate an immediate mode pins
      0        = drift_b                      the DIFFERENTIATED pin, which is
                                              what determines that mode's flow
      0        = phi_j*w_l - phi_l*w_j        the conflict rule
      0        = d/dt(sum_h y - theta_j)      the DIFFERENTIATED server latch
      0        = A_c * drift                  the differentiated capacity row
      d/dt Sigma = A Sigma + Sigma A' + Q     the covariance, when it fits

    The undifferentiated pin x_b = 0 and its derivative are BOTH present and are
    not redundant: the first makes x_b algebraic at its pinned value, the second
    is the equation that pins the flow. The undifferentiated form alone would be
    index 2, which neither RODAS nor ode15s solves.
    """
    import warnings
    from ..ode.rodas import rodas

    n = terms.nstate
    nimm = imm.n
    nlatch = len(terms.latch_mode)
    tspan = list(getattr(options, 'timespan', [0.0, np.inf]))
    if not np.isfinite(tspan[1]):
        return (np.array([0.0]), xss.reshape(1, -1), np.asarray(phiss).reshape(1, -1),
                None, np.zeros(0), [])

    maxcov = 25
    cfg = getattr(options, 'config', None)
    if cfg is not None and getattr(cfg, 'dae_maxcov', None):
        maxcov = int(cfg.dae_maxcov)
    nc = terms.cov_idx.size
    withcov = 0 < nc <= maxcov
    if not withcov and nc > maxcov:
        warnings.warn(
            'The transient covariance would add %d differential states, above the limit of %d set by '
            'options.config.dae_maxcov. Integrating the mean as a DAE with the variance held at its '
            'stationary value.' % (nc * nc, maxcov * maxcov), RuntimeWarning)

    x0 = collapse(terms, imm, terms.x0)
    ts, xs, phis = [], [], []
    tvar, Sigmas, switches = [], [], []
    tcur = float(tspan[0])
    xcur = x0.copy()
    phi = np.asarray(phiss, dtype=float).copy()
    mu = np.zeros(nlatch)
    zeta = np.asarray(zetass, dtype=float).copy()
    Scur = np.zeros((nc, nc))

    seg_max = max(8, 4 * (nimm + con.b.size + 1))
    for _ in range(seg_max):
        # THE ACTIVE SET IS RE-READ EVERY SEGMENT: a located crossing may have
        # added or released a capacity row, which changes how many algebraic
        # unknowns the state vector holds and where the covariance block starts.
        ctx = _context(terms, cons, con, imm, active)
        nact = ctx.active.size
        nvar = n + nimm + nlatch + nact
        nz = nvar + (nc * nc if withcov else 0)

        diag = np.zeros(nz)
        diag[:n] = 1.0
        for p in np.asarray(imm.pins, dtype=int).ravel():
            diag[p] = 0.0
        if withcov:
            diag[nvar:] = 1.0

        z = np.concatenate([xcur, phi, mu, zeta] + ([Scur.ravel()] if withcov else []))
        if z.size != nz:
            z = np.resize(z, nz)

        seg_t, seg_z = [], []
        hit = {'c': -1, 't': None, 'z': None}
        gprev = _event_values(z, ctx)

        def fcn(_x, y, f, _ctx=ctx):
            f[:] = _dae_rhs(np.asarray(y, dtype=float), _ctx, withcov, nc)

        def mas(am, _d=diag):
            am[0, :] = _d

        def solout(_nr, xold, xcur_, y, dense, _g=gprev, _ctx=ctx):
            yy = np.asarray(y, dtype=float)
            gcur = _event_values(yy, _ctx)
            cross = np.flatnonzero((_g > 0) & (gcur <= 0))
            _g[:] = gcur
            if cross.size == 0:
                seg_t.append(float(xcur_))
                seg_z.append(yy.copy())
                return 0
            # THE OVERSHOOTING STEP IS NOT REPORTED: the crossing is located on
            # RODAS's own interpolant over the step it has just accepted, so the
            # trajectory never carries a point beyond the constraint.
            lo, hi = float(xold), float(xcur_)
            cbest = int(cross[0])
            for _ in range(60):
                mid = 0.5 * (lo + hi)
                ymid = np.array([dense.value(i, mid) for i in range(nz)])
                if _event_values(ymid, _ctx)[cbest] > 0:
                    lo = mid
                else:
                    hi = mid
                if hi - lo <= 1e-12 * max(1.0, abs(hi)):
                    break
            hit['c'] = cbest
            hit['t'] = hi
            hit['z'] = np.array([dense.value(i, hi) for i in range(nz)])
            seg_t.append(float(hi))
            seg_z.append(hit['z'].copy())
            return -1

        y = z.copy()
        try:
            res = rodas(nz, fcn, tcur, y, float(tspan[1]), h=1e-6,
                        rtol=COARSE_TOL, atol=FINE_TOL, itol=0,
                        ijac=0, mljac=nz, mujac=nz, ifcn=1,
                        mas=mas, imas=1, mlmas=0, mumas=0,
                        solout=solout, iout=1)
            if res.idid != 1 and hit['c'] < 0:
                raise RuntimeError('RODAS returned idid=%d' % res.idid)
        except Exception as exc:
            warnings.warn('The transient DAE failed (%s); reporting the steady state instead.' % exc,
                          RuntimeWarning)
            break

        if seg_t:
            ts.extend(seg_t)
            arr = np.array(seg_z)
            xs.append(arr[:, :n])
            phis.append(arr[:, n:n + nimm])
            if withcov:
                tvar.extend(seg_t)
                Sigmas.append(arr[:, nvar:nvar + nc * nc])
            tcur = seg_t[-1]
            zend = seg_z[-1]
            xcur = zend[:n]
            phi = zend[n:n + nimm]
            mu = zend[n + nimm:n + nimm + nlatch]
            zeta = zend[n + nimm + nlatch:n + nimm + nlatch + nact]
            if withcov:
                Scur = zend[nvar:nvar + nc * nc].reshape(nc, nc)

        if hit['c'] < 0 or tcur >= tspan[1] - FINE_TOL:
            break
        imm, active, changed, zeta = _switch(terms, con, imm, active, hit['c'], xcur, zeta)
        switches.append((float(tcur), int(hit['c'])))
        if not changed:
            break

    if not ts:
        return (np.array([float(tspan[0])]), x0.reshape(1, -1),
                np.asarray(phiss).reshape(1, -1), None, np.zeros(0), switches)
    t = np.asarray(ts, dtype=float)
    xt = np.vstack(xs)
    phit = np.vstack(phis) if phis and nimm else np.zeros((t.size, nimm))
    Sigmat = None
    if withcov and Sigmas:
        blk = np.vstack(Sigmas)
        Sigmat = blk.reshape(-1, nc, nc)
    return t, xt, phit, Sigmat, np.asarray(tvar, dtype=float), switches
