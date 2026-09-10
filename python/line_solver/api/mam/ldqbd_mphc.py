"""
Exact level-dependent QBD blocks of an M/PH/c queue.

Port of matlab/src/api/mam/ldqbd_mphc.m and ph_multisets.m. The level of the
chain is the number of jobs at the station; the coordinate INSIDE a level is
the MULTISET of the phases the min(n,c) busy servers sit in, which is what
makes the construction exact for phase-type service at any number of servers.

The collapsed alternative -- one PH process run at min(n,c) times its speed --
gets the aggregate service rate right but forgets which phase each busy server
is in, turning the c servers into one fast server whose remaining work is a
single phase-type variable.

References:
    S. Asmussen and J.R. Moller, "Calculation of the steady state waiting time
    distribution in GI/PH/c and MAP/PH/c queues", Queueing Systems 37(1):9-29,
    2001.
    M. F. Neuts, "Matrix-geometric solutions in stochastic models", Johns
    Hopkins University Press, 1981.
"""

from typing import List, Optional, Sequence, Tuple

import numpy as np

# The repeating level is the widest one, so it is the size worth guarding: the
# LD-QBD recursion inverts one matrix of that order per level.
MAX_CONFIGS = 2000


def ph_multisets(p: int, k: int) -> np.ndarray:
    """Configurations of k identical servers over p service phases.

    Rows are the compositions of k into p nonnegative parts: ``M[r, i]`` is the
    number of the k busy servers sitting in phase i. There are
    ``comb(k+p-1, p-1)`` of them, the multiset count of Asmussen and Moller
    (2001) -- identical servers are exchangeable, so only the phase COUNTS
    carry information and the ordered space of size ``p**k`` collapses onto
    this one.

    The order is fixed and shared by every caller, so a configuration index
    means the same thing in each of them: the first part descends. k == 1
    therefore yields the identity rows e_1 ... e_p in phase order, which is
    what makes the c == 1 case of :func:`ldqbd_mphc` coincide with plain phase
    indexing.
    """
    if k == 0:
        return np.zeros((1, p), dtype=int)
    if p == 1:
        return np.array([[k]], dtype=int)
    rows = []
    for first in range(k, -1, -1):
        rest = ph_multisets(p - 1, k - first)
        rows.append(np.hstack([np.full((rest.shape[0], 1), first, dtype=int), rest]))
    return np.vstack(rows)


def ldqbd_mphc(D0, D1, alpha, c, arr_rate,
               sf: Optional[Sequence[float]] = None
               ) -> Tuple[List[np.ndarray], List[np.ndarray], List[np.ndarray]]:
    """Block-tridiagonal generator of an M/PH/c queue with level-dependent arrivals.

    Args:
        D0: service sub-generator (p x p), phase changes without completion.
        D1: service completion block (p x p); ``D1 = (-D0 @ 1) @ alpha`` for a PH.
        alpha: length-p vector a server starts each new job in.
        c: number of identical servers (>= 1; capped at the top level).
        arr_rate: length Nlev+1; ``arr_rate[n]`` is the arrival rate out of level n.
        sf: optional length-Nlev multiplier on the station's TOTAL service rate
            at level n (load dependence). Each busy server then runs at
            ``sf[n-1]/min(n,c)`` of nominal, so ``sf[n-1] == min(n,c)``
            reproduces the unscaled queue exactly.

    Returns:
        ``(Q0, Q1, Q2)`` with ``Q0[n]`` the upward block of level n (length
        Nlev), ``Q1[n]`` the local block of level n (length Nlev+1) and
        ``Q2[n]`` the downward block of level n+1 -> n (length Nlev, so
        ``Q2[n-1]`` is the block leaving level n, matching the argument order
        the native :func:`ldqbd` takes).

    Level sizes grow over the boundary levels 0..c and repeat above them, so
    the blocks joining differently sized neighbours are rectangular; ldqbd,
    ldqbd_R and ldqbd_pi all accept that heterogeneity.
    """
    D0 = np.atleast_2d(np.asarray(D0, dtype=float))
    D1 = np.atleast_2d(np.asarray(D1, dtype=float))
    alpha = np.asarray(alpha, dtype=float).ravel()
    arr_rate = np.asarray(arr_rate, dtype=float).ravel()
    p = D0.shape[0]
    nlev = arr_rate.size - 1

    if nlev < 1:
        raise ValueError('ldqbd_mphc needs at least one level above the empty one.')
    if D0.shape[1] != p or D1.shape != (p, p) or alpha.size != p:
        raise ValueError('D0, D1 and alpha must all have the same order.')
    if sf is not None:
        sf = np.asarray(sf, dtype=float).ravel()
        if sf.size < nlev:
            raise ValueError('sf must give one total-service-rate factor per level 1..nlev.')

    # Servers that can never be busy do not need a coordinate: above the top
    # level there is nothing left to serve.
    cmax = nlev if not np.isfinite(c) else min(max(1, int(round(c))), nlev)

    cfg = [ph_multisets(p, k) for k in range(cmax + 1)]
    pos = [{tuple(row): r for r, row in enumerate(Ck)} for Ck in cfg]
    ncfg = [Ck.shape[0] for Ck in cfg]

    if ncfg[cmax] > MAX_CONFIGS:
        raise ValueError(
            'the exact M/PH/c chain needs comb(%d+%d-1,%d-1) = %d configurations per level '
            'for %d servers and %d service phases, above the %d the level-by-level inverses '
            'can carry. Use fewer phases (a lower-order fit), fewer servers, or '
            'SolverCTMC/SolverLDES on this model.'
            % (cmax, p, p, ncfg[cmax], cmax, p, MAX_CONFIGS))

    t = D1 @ np.ones(p)  # completion rate out of each phase, summed over targets

    # Structural blocks per busy-server count: LOC the within-level phase
    # changes (with the full outflow on its diagonal), UP the entry of a newly
    # busy server, DN a completion that leaves a server idle.
    LOC: List[np.ndarray] = []
    UP: List[Optional[np.ndarray]] = []
    DN: List[Optional[np.ndarray]] = []
    for k in range(cmax + 1):
        Ck = cfg[k]
        nk = ncfg[k]
        Lk = np.zeros((nk, nk))
        for row in range(nk):
            m = Ck[row]
            for i in range(p):
                if m[i] == 0:
                    continue
                for j in range(p):
                    if j == i:
                        continue
                    mm = m.copy()
                    mm[i] -= 1
                    mm[j] += 1
                    Lk[row, pos[k][tuple(mm)]] += m[i] * D0[i, j]
                # D0[i,i] is the total outflow of phase i, completions included
                Lk[row, row] += m[i] * D0[i, i]
        LOC.append(Lk)

        if k < cmax:
            Uk = np.zeros((nk, ncfg[k + 1]))
            for row in range(nk):
                m = Ck[row]
                for j in range(p):
                    mm = m.copy()
                    mm[j] += 1
                    Uk[row, pos[k + 1][tuple(mm)]] += alpha[j]
            UP.append(Uk)
        else:
            UP.append(None)

        if k > 0:
            Dk = np.zeros((nk, ncfg[k - 1]))
            for row in range(nk):
                m = Ck[row]
                for i in range(p):
                    if m[i] == 0:
                        continue
                    mm = m.copy()
                    mm[i] -= 1
                    Dk[row, pos[k - 1][tuple(mm)]] += m[i] * t[i]
            DN.append(Dk)
        else:
            DN.append(None)

    # A completion at a full server bank takes the next waiting job at once, so
    # the server stays busy and only its phase moves: the repeating down block.
    Cc = cfg[cmax]
    nc = ncfg[cmax]
    CDEP = np.zeros((nc, nc))
    for row in range(nc):
        m = Cc[row]
        for i in range(p):
            if m[i] == 0:
                continue
            for j in range(p):
                mm = m.copy()
                mm[i] -= 1
                mm[j] += 1
                CDEP[row, pos[cmax][tuple(mm)]] += m[i] * D1[i, j]

    # Per-server speed. Without load dependence every busy server runs at its
    # nominal rate; with it, the aggregate sf(n) is shared over the busy
    # servers. sf(n) == min(n,c) is passed through as exactly 1 so the unscaled
    # chain is reproduced bit for bit.
    speed = np.ones(nlev)
    if sf is not None:
        for n in range(1, nlev + 1):
            b = min(n, cmax)
            if sf[n - 1] != b:
                speed[n - 1] = sf[n - 1] / b

    Q0: List[np.ndarray] = []
    Q1: List[np.ndarray] = []
    Q2: List[np.ndarray] = []

    Q1.append(np.array([[-arr_rate[0]]]))  # level 0: arrivals only
    for n in range(1, nlev + 1):
        b = min(n, cmax)
        Q1.append(speed[n - 1] * LOC[b] - arr_rate[n] * np.eye(ncfg[b]))

    for n in range(nlev):
        if n < cmax:
            Q0.append(arr_rate[n] * UP[n])            # a free server takes the job
        else:
            Q0.append(arr_rate[n] * np.eye(ncfg[cmax]))  # the job waits, phases unchanged

    for n in range(1, nlev + 1):
        if n <= cmax:
            Q2.append(speed[n - 1] * DN[n])           # the server falls idle
        else:
            Q2.append(speed[n - 1] * CDEP)            # the server takes the next job

    return Q0, Q1, Q2
