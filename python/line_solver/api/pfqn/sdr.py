"""Product-form state-dependent routing, Krzesinski (1987).

Krzesinski, A. E., "Multiclass Queueing Networks with State-Dependent Routing",
Performance Evaluation 7(2):125-143, 1987, the multiclass generalization of
Towsley, D., "Queuing Network Models with State-Dependent Routing",
J. ACM 27(2):323-337, 1980.

Branch index 1 denotes the complement M-V and is unused; the SDR branches are
numbered 2..B, following the paper's own indexing so that d[t][b] transcribes
straight from the text.
"""

import math

import numpy as np

__all__ = ['pfqn_sdrcoeff', 'pfqn_sdrprob', 'pfqn_sdr', 'pfqn_sdrvisits', 'pfqn_sdrmva']


class SDRCoeff(object):
    """Derived coefficients of an SDR structure, eqs. (11)-(14)."""

    __slots__ = ('T', 'B', 'level', 'C', 'd', 'inA', 'Dtt', 'Dprev', 'mmax',
                 'vmax', 'entry', 'departure', 'branch', 'entryOf', 'departureOf')


def _get(sdr, name):
    if isinstance(sdr, dict):
        if name not in sdr:
            raise ValueError("the SDR structure is missing field '%s'" % name)
        return sdr[name]
    if not hasattr(sdr, name):
        raise ValueError("the SDR structure is missing field '%s'" % name)
    return getattr(sdr, name)


def pfqn_sdrcoeff(sdr):
    """Validate an SDR structure and return its derived coefficients.

    The structure carries, with branch index 1 reserved for the complement M-V:

      entry, departure    center indices of e and d of Q(V,V)
      branch[b]           center indices of branch b, b >= 2; branch[0] unused
      entryOf[b]          center index of the branch entry e(b)
      departureOf[b]      center index of the branch departure d(b)
      level[b]            the unique t with B_b in V_t - V_{t+1}
      C                   length T, the coefficients C_t of eq. (11)
      d                   T x B, d[t][b] read for 2 <= b <= B, 1 <= t <= level[b]

    The population bounds mmax and vmax are consequences of the coefficients,
    not independent inputs: with C_t < 0 the routing enforces them by itself.
    """
    branch = list(_get(sdr, 'branch'))
    C = np.asarray(_get(sdr, 'C'), dtype=float).ravel()
    d = np.asarray(_get(sdr, 'd'), dtype=float)
    level = np.asarray(_get(sdr, 'level')).ravel()
    entryOf = np.asarray(_get(sdr, 'entryOf')).ravel()
    departureOf = np.asarray(_get(sdr, 'departureOf')).ravel()
    entry = int(_get(sdr, 'entry'))
    departure = int(_get(sdr, 'departure'))

    B = len(branch)
    T = C.size
    if B < 2:
        raise ValueError('an SDR structure must declare at least one branch '
                         '(branch indices start at 2)')
    if T < 1:
        raise ValueError('an SDR structure must declare at least one level of '
                         'subnetwork nesting')
    if level.size != B:
        raise ValueError('SDR level must have one entry per branch index')
    lv = level[1:].astype(int)
    if np.any(lv < 1) or np.any(lv > T):
        raise ValueError('SDR branch levels must be integers in 1..%d' % T)
    if d.shape[0] < T or d.shape[1] < B:
        raise ValueError('SDR coefficient matrix d must be at least %dx%d' % (T, B))

    # The nesting V_1 > ... > V_T must be strict, else two hierarchically
    # adjacent subnetworks coincide and the ratio of eq. (10) is not the paper's
    for t in range(1, T + 1):
        if not np.any(lv == t):
            raise ValueError('SDR level %d carries no branch: the subnetwork '
                             'nesting must be strict' % t)

    seen = set()
    for b in range(1, B):
        sb = [int(x) for x in np.ravel(branch[b])]
        if not sb:
            raise ValueError('SDR branch %d is empty' % (b + 1))
        if seen.intersection(sb):
            raise ValueError('SDR branches must be mutually disjoint')
        if entry in sb or departure in sb:
            raise ValueError('the entry and departure centers of Q(V,V) must '
                             'not belong to any branch')
        if int(entryOf[b]) not in sb or int(departureOf[b]) not in sb:
            raise ValueError('the entry and departure centers of SDR branch %d '
                             'must belong to that branch' % (b + 1))
        seen.update(sb)

    inA = [None] * (T + 1)
    Dtt = np.zeros(T + 1)
    Dprev = np.zeros(T + 1)
    for t in range(1, T + 1):
        inA[t] = [b for b in range(1, B) if int(level[b]) >= t]
        Dtt[t] = sum(d[t - 1, b] for b in inA[t])
        if t > 1:
            Dprev[t] = sum(d[t - 2, b] for b in inA[t])

    mmax = np.full(B, np.inf)
    for b in range(1, B):
        t = int(level[b])
        if C[t - 1] < 0:
            mmax[b] = math.floor(d[t - 1, b] / (-C[t - 1]))
    vmax = np.full(T + 1, np.inf)
    for t in range(1, T + 1):
        if C[t - 1] < 0:
            vmax[t] = math.floor(Dtt[t] / (-C[t - 1]))
        if t > 1 and C[t - 2] < 0:
            vmax[t] = min(vmax[t], math.floor(Dprev[t] / (-C[t - 2])))

    c = SDRCoeff()
    c.T = T
    c.B = B
    c.level = level.astype(int)
    c.C = C
    c.d = d
    c.inA = inA
    c.Dtt = Dtt
    c.Dprev = Dprev
    c.mmax = mmax
    c.vmax = vmax
    c.entry = entry
    c.departure = departure
    c.branch = [None] + [[int(x) for x in np.ravel(branch[b])] for b in range(1, B)]
    c.entryOf = entryOf.astype(int)
    c.departureOf = departureOf.astype(int)
    return c


def pfqn_sdrprob(sdr, n):
    """SDR routing probabilities of eq. (10).

    Returns (P, Ped) where P[b] is the probability of proceeding from the entry
    center e of Q(V,V) to the entry center of branch b (P[0] is zero, branch
    index 1 being the complement M-V) and Ped = 1 - sum(P) is the probability of
    proceeding directly to the departure center d, that is of being denied entry
    into Q(V,V) and returned to e.

    These probabilities are chain independent: they read the total branch and
    subnetwork populations, not the per-chain ones. The chain-dependent form of
    eq. (1) has no published product form and is not implemented.

    A branch population with delta_tb(m_b) < 0 lies beyond the bound that SDR
    enforces itself and is unreachable, so the probability there is zero.
    """
    c = sdr if isinstance(sdr, SDRCoeff) else pfqn_sdrcoeff(sdr)
    n = np.asarray(n, dtype=float).ravel()

    m = np.zeros(c.B)
    for b in range(1, c.B):
        m[b] = float(np.sum(n[c.branch[b]]))
    v = np.zeros(c.T + 1)
    for t in range(1, c.T + 1):
        v[t] = sum(m[b] for b in c.inA[t])

    om = np.zeros(c.T + 1)
    omprev = np.ones(c.T + 1)
    for s in range(1, c.T + 1):
        om[s] = c.C[s - 1] * v[s] + c.Dtt[s]
        if s > 1:
            omprev[s] = c.C[s - 2] * v[s] + c.Dprev[s]

    P = np.zeros(c.B)
    for b in range(1, c.B):
        t = int(c.level[b])
        if any(om[s] <= 0 for s in range(1, t + 1)):
            continue  # eq. (10): the branch is closed to new arrivals
        delta = c.C[t - 1] * m[b] + c.d[t - 1, b]
        if delta <= 0:
            continue
        ratio = 1.0
        for s in range(1, t + 1):
            ratio *= omprev[s] / om[s]
        P[b] = delta * ratio
    return P, 1.0 - float(np.sum(P))


def _cumlog(Cv, dv, nmax):
    """Cumulative log-product prod_{k=0}^{n-1} (Cv k + dv) for n = 0..nmax.

    ok[q, n] is False once a nonpositive factor has been met, which is where the
    SDR population bound closes the branch or the subnetwork.
    """
    Cv = np.atleast_1d(np.asarray(Cv, dtype=float))
    dv = np.atleast_1d(np.asarray(dv, dtype=float))
    p = Cv.size
    lc = np.zeros((p, nmax + 1))
    ok = np.ones((p, nmax + 1), dtype=bool)
    for q in range(p):
        for nn in range(1, nmax + 1):
            f = Cv[q] * (nn - 1) + dv[q]
            if not ok[q, nn - 1] or f <= 0:
                ok[q, nn] = False
                lc[q, nn] = -np.inf
            else:
                lc[q, nn] = lc[q, nn - 1] + math.log(f)
    return lc, ok


def _compositions(n, m):
    """All nonnegative integer m-vectors summing to n, one per row."""
    if m == 1:
        return np.array([[n]], dtype=int)
    out = []
    for k in range(n + 1):
        tail = _compositions(n - k, m - 1)
        out.append(np.hstack([np.full((tail.shape[0], 1), k, dtype=int), tail]))
    return np.vstack(out)


def _states(M, N):
    """Every M x J population matrix with column sums N, flattened per chain."""
    per = [_compositions(int(N[j]), M) for j in range(len(N))]
    st = per[0]
    for j in range(1, len(N)):
        a, b = st, per[j]
        na, nb = a.shape[0], b.shape[0]
        st = np.hstack([np.tile(a, (nb, 1)),
                        np.repeat(b, na, axis=0)])
    return st


def pfqn_sdr(S, xi, N, sdr, alpha=None):
    """Exact product form of eq. (16).

    S and xi are M x J: the mean service times 1/mu_ij and the coefficients
    xi_ij of Section 3.2. They are required separately rather than as their
    product because under SDR the xi are not visit ratios, so the per-center
    throughputs cannot be recovered from the demands alone.

    alpha is an optional M x sum(N) matrix of load-dependent rate scalings,
    alpha[i, k-1] = alpha_i(k). Defaults to a fixed-rate center; use k for an
    infinite server and min(k, c) for a c-server center.

    Returns (Q, X, U, R, G, lG, prob, states), all M x J except the last three.
    X holds the per-center chain throughputs, U = X * S the mean number in
    service, and R = Q / X the response time at the center.
    """
    S = np.atleast_2d(np.asarray(S, dtype=float))
    xi = np.atleast_2d(np.asarray(xi, dtype=float))
    M, J = S.shape
    N = np.asarray(N, dtype=int).ravel()
    if N.size != J:
        raise ValueError('the population vector must have one entry per chain')
    if xi.shape != (M, J):
        raise ValueError('S and xi must have the same shape')

    Ntot = int(np.sum(N))
    if alpha is None:
        alpha = np.ones((M, max(1, Ntot)))
    else:
        alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
        if alpha.shape[1] < Ntot:
            alpha = np.hstack([alpha, np.ones((M, Ntot - alpha.shape[1]))])

    c = pfqn_sdrcoeff(sdr)
    allc = [c.entry, c.departure] + [i for b in range(1, c.B) for i in c.branch[b]]
    if max(allc) >= M:
        raise ValueError('the SDR structure references a center index beyond '
                         'the number of centers')

    gamma = xi * S

    logbeta = np.zeros((M, Ntot + 1))
    for i in range(M):
        for k in range(1, Ntot + 1):
            logbeta[i, k] = logbeta[i, k - 1] + math.log(alpha[i, k - 1])

    lv = c.level[1:].astype(int)
    logDelta, okDelta = _cumlog([c.C[t - 1] for t in lv],
                                [c.d[t - 1, b] for b, t in zip(range(1, c.B), lv)],
                                Ntot)
    logOmTT = np.zeros((c.T + 1, Ntot + 1))
    okOmTT = np.ones((c.T + 1, Ntot + 1), dtype=bool)
    logOmPrev = np.zeros((c.T + 1, Ntot + 1))
    okOmPrev = np.ones((c.T + 1, Ntot + 1), dtype=bool)
    for t in range(1, c.T + 1):
        lc, ok = _cumlog(c.C[t - 1], c.Dtt[t], Ntot)
        logOmTT[t, :], okOmTT[t, :] = lc[0], ok[0]
        if t > 1:
            lc, ok = _cumlog(c.C[t - 2], c.Dprev[t], Ntot)
            logOmPrev[t, :], okOmPrev[t, :] = lc[0], ok[0]

    states = _states(M, N)
    ns = states.shape[0]
    logw = np.full(ns, -np.inf)
    for s in range(ns):
        nmat = states[s, :].reshape(J, M).T
        ni = nmat.sum(axis=1)
        lw = 0.0
        ok = True
        for i in range(M):
            lw += math.lgamma(ni[i] + 1) - logbeta[i, ni[i]]
            for j in range(J):
                if nmat[i, j] > 0:
                    if gamma[i, j] <= 0:
                        ok = False
                        break
                    lw += nmat[i, j] * math.log(gamma[i, j]) - math.lgamma(nmat[i, j] + 1)
            if not ok:
                break
        if not ok:
            continue
        m = np.zeros(c.B, dtype=int)
        for b in range(1, c.B):
            m[b] = int(np.sum(ni[c.branch[b]]))
        for b in range(1, c.B):
            if not okDelta[b - 1, m[b]]:
                ok = False
                break
            lw += logDelta[b - 1, m[b]]
        if not ok:
            continue
        for t in range(1, c.T + 1):
            v = int(sum(m[b] for b in c.inA[t]))
            if not okOmTT[t, v]:
                ok = False
                break
            lw -= logOmTT[t, v]
            if t > 1:
                if not okOmPrev[t, v]:
                    ok = False
                    break
                lw += logOmPrev[t, v]
        if not ok:
            continue
        logw[s] = lw

    lmax = np.max(logw)
    if not np.isfinite(lmax):
        raise ValueError('the SDR network has no reachable state at the given '
                         'populations: the routing coefficients forbid every state')
    w = np.exp(logw - lmax)
    Gs = float(np.sum(w))
    prob = w / Gs
    lG = lmax + math.log(Gs)
    G = math.exp(lG) if lG < 700 else np.inf

    Q = np.zeros((M, J))
    X = np.zeros((M, J))
    for s in range(ns):
        if prob[s] == 0.0:
            continue
        nmat = states[s, :].reshape(J, M).T
        ni = nmat.sum(axis=1)
        Q += prob[s] * nmat
        for i in range(M):
            if ni[i] > 0:
                with np.errstate(divide='ignore', invalid='ignore'):
                    contrib = prob[s] * alpha[i, ni[i] - 1] * (nmat[i, :] / ni[i]) / S[i, :]
                X[i, :] += np.where(S[i, :] > 0, contrib, 0.0)
    U = X * S
    R = np.zeros((M, J))
    nz = X > 0
    R[nz] = Q[nz] / X[nz]
    return Q, X, U, R, G, lG, prob, states


def pfqn_sdrvisits(sdr, P):
    """Coefficients xi of Section 3.2.

    P is M x M x J: P[x, y, j] is the state-independent probability that a
    chain j customer leaving center x proceeds to center y. The state-dependent
    arcs out of the entry center are not part of P and are ignored if present.

    Three rules fix the coefficients: the complement M-V obeys the ordinary
    traffic equations with the whole SDR subnetwork collapsed into a single
    e -> d arc of probability one; every branch obeys its own traffic equations
    driven by an injection of xi_e at its entry center; and xi_e = 1.

    The paper states xi_ij = xi_ej for the branch entry and departure centers
    and works out only single-center branches. The traffic equations above are
    the reading that extends it: they return xi_{d(b)} = xi_e because a customer
    leaves a branch only through d(b), and xi_{e(b)} = xi_e whenever e(b) takes
    no internal feedback. They have been checked against a brute-force CTMC on a
    branch that does take such feedback, where the paper's literal rule fails.

    These xi are not relative visit counts: the rate at which customers enter a
    branch is state dependent, so a ratio of two xi carries no flow meaning.
    """
    from ..mc.dtmc import dtmc_solve

    c = pfqn_sdrcoeff(sdr)
    P = np.asarray(P, dtype=float)
    if P.ndim == 2:
        P = P[:, :, np.newaxis]
    M = P.shape[0]
    if P.shape[1] != M:
        raise ValueError('the SIR routing array must be square in its first '
                         'two dimensions')
    J = P.shape[2]
    xi = np.zeros((M, J))

    inV = np.zeros(M, dtype=bool)
    for b in range(1, c.B):
        inV[c.branch[b]] = True
    mv = [i for i in range(M) if not inV[i]]
    if c.entry not in mv or c.departure not in mv:
        raise ValueError('the entry and departure centers of Q(V,V) must lie '
                         'outside every branch')
    ie = mv.index(c.entry)
    idp = mv.index(c.departure)

    for j in range(J):
        Pmv = P[np.ix_(mv, mv, [j])][:, :, 0].copy()
        Pmv[ie, :] = 0.0
        Pmv[ie, idp] = 1.0
        if np.any(np.abs(Pmv.sum(axis=1) - 1.0) > 1e-3):
            raise ValueError('the SIR routing of chain %d does not keep '
                             'customers inside the complement M-V' % j)
        xmv = np.ravel(dtmc_solve(Pmv))
        if xmv[ie] <= 0:
            raise ValueError('the entry center of Q(V,V) is unreachable in '
                             'chain %d' % j)
        xmv = xmv / xmv[ie]
        for k, i in enumerate(mv):
            xi[i, j] = xmv[k]

        for b in range(1, c.B):
            sb = c.branch[b]
            Pbb = P[np.ix_(sb, sb, [j])][:, :, 0]
            inj = np.zeros(len(sb))
            inj[sb.index(int(c.entryOf[b]))] = xi[c.entry, j]
            xi[sb, j] = np.linalg.solve((np.eye(len(sb)) - Pbb).T, inj)
    return xi


def _lattice(N):
    """Every population vector V with 0 <= V <= N, mixed-radix ordered."""
    J = len(N)
    tot = 1
    for j in range(J):
        tot *= int(N[j]) + 1
    latt = np.zeros((tot, J), dtype=int)
    for r in range(tot):
        rem = r
        for j in range(J):
            latt[r, j] = rem % (int(N[j]) + 1)
            rem //= (int(N[j]) + 1)
    return latt


def _key(V, N):
    """Mixed-radix row of V in the lattice of _lattice."""
    k, mul = 0, 1
    for j in range(len(N)):
        k += int(V[j]) * mul
        mul *= int(N[j]) + 1
    return k


def _sublattice(V):
    """Every L with 0 <= L <= V."""
    J = len(V)
    tot = 1
    for j in range(J):
        tot *= int(V[j]) + 1
    sub = np.zeros((tot, J), dtype=int)
    for r in range(tot):
        rem = r
        for j in range(J):
            sub[r, j] = rem % (int(V[j]) + 1)
            rem //= (int(V[j]) + 1)
    return sub


def _delta(dk, n):
    """delta_i(n) = dk - n after rescaling to C_t = -1, or 1 outside Q(V,V)."""
    return 1.0 if np.isnan(dk) else dk - n


def _omega_cum(c, t, v):
    """Omega_{t-1,t}(v)/Omega_tt(v), the cumulative ratio of eq. (16)."""
    num = den = 1.0
    for k in range(v):
        f = -k + c.Dtt[t]
        if f <= 0:
            return 0.0
        den *= f
        if t > 1:
            g = -k + c.Dprev[t]
            if g <= 0:
                return 0.0
            num *= g
    return num / den


def _omega_step(c, t, v):
    """omega_{t-1,t}(v)/omega_tt(v), the single-step ratio of eq. (10).

    Distinct from _omega_cum: the routing probability carries the lowercase
    omega, the normalizing constant the uppercase one.
    """
    den = -v + c.Dtt[t]
    if den <= 0:
        return 0.0
    if t > 1:
        num = -v + c.Dprev[t]
        if num <= 0:
            return 0.0
    else:
        num = 1.0  # omega_{0,1} is one
    return num / den


def _set_mva(cidx, dvec, gamma, alpha, N, latt, Ntot):
    """Section 4.2.1: MVA of the centers cidx in isolation.

    dvec[k] is the SDR admission coefficient of centre cidx[k], so
    delta_k(n) = dvec[k] - n; NaN means delta = 1, the ordinary BCMP arrival
    theorem, which is what the complement M-V uses.

    The factor [d_ti - Q_i(V - 1_j)] that the paper places in both W_ij and the
    denominator of T_j cancels between them and is not formed here: it can
    vanish at the population bound, where it would divide by zero without
    changing any queue length or throughput.
    """
    M, J = gamma.shape
    nc = len(cidx)
    nl = latt.shape[0]
    Qs = np.zeros((M, J, nl))
    Ts = np.zeros((J, nl))
    gs = np.zeros(nl)
    z = _key(np.zeros(J, dtype=int), N)
    if nc == 0:
        for v in range(nl):
            gs[v] = 1.0 if latt[v].sum() == 0 else 0.0
        return Qs, Ts, gs
    Ps = np.zeros((nc, Ntot + 1, nl))
    gs[z] = 1.0
    Ps[:, 0, z] = 1.0

    for v in sorted(range(nl), key=lambda r: latt[r].sum()):
        V = latt[v]
        vv = int(V.sum())
        if vv == 0:
            continue
        A = np.zeros((nc, J))
        vm = [-1] * J
        for j in range(J):
            if V[j] == 0:
                continue
            Vm = V.copy()
            Vm[j] -= 1
            vm[j] = _key(Vm, N)
            for k in range(nc):
                acc = 0.0
                for n in range(1, vv + 1):
                    df = _delta(dvec[k], n - 1)
                    if df <= 0:
                        break
                    acc += n * (df / alpha[cidx[k], n - 1]) * Ps[k, n - 1, vm[j]]
                A[k, j] = acc
        for j in range(J):
            if V[j] == 0:
                continue
            den = sum(gamma[cidx[k], j] * A[k, j] for k in range(nc))
            if den > 0:
                Ts[j, v] = V[j] / den
        for k in range(nc):
            for j in range(J):
                if V[j] == 0:
                    continue
                Qs[cidx[k], j, v] = gamma[cidx[k], j] * Ts[j, v] * A[k, j]
        for k in range(nc):
            tot = 0.0
            for n in range(1, vv + 1):
                df = _delta(dvec[k], n - 1)
                if df <= 0:
                    break
                acc = 0.0
                for j in range(J):
                    if V[j] == 0:
                        continue
                    acc += gamma[cidx[k], j] * Ts[j, v] * Ps[k, n - 1, vm[j]]
                Ps[k, n, v] = (df / alpha[cidx[k], n - 1]) * acc
                tot += Ps[k, n, v]
            Ps[k, 0, v] = 1.0 - tot
        for j in range(J):
            if V[j] > 0 and Ts[j, v] > 0:
                gs[v] = gs[vm[j]] / Ts[j, v]
                break
    return Qs, Ts, gs


def pfqn_sdrmva(S, xi, N, sdr, alpha=None):
    """Section 4 MVA and convolution of a network with state-dependent routing.

    Krzesinski (1987), Performance Evaluation 7:125-143, Section 4. Same
    signature and outputs as pfqn_sdr, which evaluates eq. (16) exactly by state
    enumeration, so the two are directly comparable. This routine costs
    O(J T M (V_1...V_J)^2) rather than the size of the state space.

    Restrictions, both from the paper: every SDR branch must be a SINGLE centre,
    and every C_t must be negative. A C_t other than -1 is rescaled internally,
    which leaves eqs. (10) and (16) unchanged because the factors telescope.

    Returns (Q, X, U, R, lG).
    """
    S = np.atleast_2d(np.asarray(S, dtype=float))
    xi = np.atleast_2d(np.asarray(xi, dtype=float))
    M, J = S.shape
    N = np.asarray(N, dtype=int).ravel()
    Ntot = int(N.sum())
    if alpha is None:
        alpha = np.ones((M, max(1, Ntot)))
    else:
        alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
        if alpha.shape[1] < Ntot:
            alpha = np.hstack([alpha, np.ones((M, Ntot - alpha.shape[1]))])

    c0 = pfqn_sdrcoeff(sdr)
    for b in range(1, c0.B):
        if len(c0.branch[b]) != 1:
            raise ValueError(
                'pfqn_sdrmva requires every SDR branch to hold a single centre: the MVA and '
                'convolution of Krzesinski (1987) Section 4 is stated that way and its general '
                'case is in an unpublished technical report. Use pfqn_sdr, which evaluates '
                'eq. (16) exactly for any branch topology.')
    if np.any(c0.C >= 0):
        raise ValueError(
            'pfqn_sdrmva requires every C_t to be negative: Section 2.5 assumes it and Section 4 '
            'is written for C_t = -1. A nonnegative C_t leaves the branch populations unbounded.')

    # Rescale each level to C_t = -1, which leaves eqs. (10) and (16) unchanged
    sdr1 = dict(sdr) if isinstance(sdr, dict) else dict(
        entry=sdr.entry, departure=sdr.departure, branch=sdr.branch, entryOf=sdr.entryOf,
        departureOf=sdr.departureOf, level=sdr.level, C=sdr.C, d=sdr.d)
    Cnew = np.asarray(c0.C, dtype=float).copy()
    dnew = np.asarray(c0.d, dtype=float).copy()
    for t in range(c0.T):
        k = -Cnew[t]
        Cnew[t] = -1.0
        dnew[t, :] = dnew[t, :] / k
    sdr1['C'] = Cnew
    sdr1['d'] = dnew
    c = pfqn_sdrcoeff(sdr1)

    gamma = xi * S
    latt = _lattice(N)
    nl = latt.shape[0]

    inV = np.zeros(M, dtype=bool)
    dvec = np.full(M, np.nan)
    lvl = np.zeros(M, dtype=int)
    for b in range(1, c.B):
        i = c.branch[b][0]
        inV[i] = True
        lvl[i] = int(c.level[b])
        dvec[i] = c.d[int(c.level[b]) - 1, b]
    mv = [i for i in range(M) if not inV[i]]

    # Sec. 4.1: the complement is an ordinary BCMP subnetwork, delta = 1
    Qc, Tc, gc = _set_mva(mv, np.full(len(mv), np.nan), gamma, alpha, N, latt, Ntot)

    # Sec. 4.2: outermost level last
    Gin = np.zeros(nl)
    Gin[_key(np.zeros(J, dtype=int), N)] = 1.0   # G(L, V_{T+1}) = [L == 0]
    Qin = np.zeros((M, J, nl))
    Tin = np.zeros((J, nl))
    Ain = np.zeros((M, nl))
    for t in range(c.T, 0, -1):
        St = [i for i in range(M) if lvl[i] == t]
        Qt, Tt, gt = _set_mva(St, dvec[St], gamma, alpha, N, latt, Ntot)
        Gnew = np.zeros(nl)
        Qnew = np.zeros((M, J, nl))
        Tnew = np.zeros((J, nl))
        Anew = np.zeros((M, nl))
        inner = [i for i in range(M) if inV[i] and lvl[i] > t]
        for v in range(nl):
            V = latt[v]
            vv = int(V.sum())
            omr = _omega_cum(c, t, vv)
            if omr == 0.0:
                continue
            anum = np.zeros(M)
            for L in _sublattice(V):
                iL = _key(L, N)
                iVL = _key(V - L, N)
                pb = omr * gt[iVL] * Gin[iL]
                if pb == 0.0:
                    continue
                Gnew[v] += pb
                for i in St:
                    Qnew[i, :, v] += Qt[i, :, iVL] * pb
                Tnew[:, v] += Tt[:, iVL] * pb
                for i in inner:
                    Qnew[i, :, v] += Qin[i, :, iL] * pb
                    anum[i] += Ain[i, iL] * pb
            if Gnew[v] > 0:
                Qnew[:, :, v] /= Gnew[v]
                Tnew[:, v] /= Gnew[v]
                # eq. (10) carries, besides delta_ti, the single-step ratios
                # omega_{s-1,s}(v_s)/omega_ss(v_s) down to the centre's own
                # level. Conditioning on the population of Q(V,V_t) fixes v_t,
                # so that ratio leaves the expectation and the rest recurses
                # through the same convolution as the queue lengths. Section
                # 4.2.3 prints only [d_1i - Q_i], which drops these ratios.
                st = _omega_step(c, t, vv)
                for i in St:
                    Anew[i, v] = st * (dvec[i] - Qnew[i, :, v].sum())
                for i in inner:
                    Anew[i, v] = st * anum[i] / Gnew[v]
        Gin, Qin, Tin, Ain = Gnew, Qnew, Tnew, Anew

    # Sec. 4.3: convolve against the complement at every population, because the
    # per-centre throughputs below read the network at N - 1_j
    Qall = np.zeros((M, J, nl))
    Tall = np.zeros((J, nl))
    Aall = np.zeros((M, nl))
    Gall = np.zeros(nl)
    for v in range(nl):
        Np = latt[v]
        for V in _sublattice(Np):
            iV = _key(V, N)
            iC = _key(Np - V, N)
            pb = gc[iC] * Gin[iV]
            if pb == 0.0:
                continue
            Gall[v] += pb
            for i in mv:
                Qall[i, :, v] += Qc[i, :, iC] * pb
            for i in range(M):
                if inV[i]:
                    Qall[i, :, v] += Qin[i, :, iV] * pb
                    Aall[i, v] += Ain[i, iV] * pb
            Tall[:, v] += Tc[:, iC] * pb
        if Gall[v] > 0:
            Qall[:, :, v] /= Gall[v]
            Tall[:, v] /= Gall[v]
            Aall[:, v] /= Gall[v]

    vN = _key(N, N)
    Q = Qall[:, :, vN]
    lG = math.log(Gall[vN])

    # Shifting one chain j customer off centre i in the weight of eq. (16) gives
    #   alpha_i(n_i)(n_ij/n_i) w(n) = gamma_ij P_{e,e(i)}(n - e_ij) w(n - e_ij),
    # so summing over states yields the exact identity
    #   T_ij = xi_ij T_j(N,M) E_{N-1_j}[P_{e,e(i)}],
    # the arrival-theorem statement that a departing customer sees N - 1_j.
    X = np.zeros((M, J))
    for j in range(J):
        if N[j] == 0:
            continue
        Nm = N.copy()
        Nm[j] -= 1
        vm = _key(Nm, N)
        for i in range(M):
            if inV[i]:
                X[i, j] = xi[i, j] * Tall[j, vN] * Aall[i, vm]
            else:
                X[i, j] = xi[i, j] * Tall[j, vN]
    U = X * S
    R = np.zeros((M, J))
    nz = X > 0
    R[nz] = Q[nz] / X[nz]
    return Q, X, U, R, lG
