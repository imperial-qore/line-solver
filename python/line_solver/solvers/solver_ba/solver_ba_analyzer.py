"""Bound-analysis dispatch for SolverBA (native Python port of MATLAB
matlab/src/solvers/BA/solver_ba_analyzer.m).

Operates on a parsed MVA-style solver ``s`` (exposing nstations, nclasses,
rates, demands, njobs, nservers, sched, _sn) and returns per-station
[Q,U,R,T,C,X] arrays. Single-class families derive V/D/Z from the folded
demand matrix (demands[i] = V[i]/rate[i], so V = demands*rate, D = demands at
queueing stations, Z = sum of demands at INF stations).
"""

import time
import numpy as np

from ...api.pfqn.bounds import (
    pfqn_xzgsblow, pfqn_xzgsbup, pfqn_qzgblow, pfqn_qzgbup, pfqn_mwrbb,
    pfqn_harel_bounds,
)
from ...api.pfqn.bound_hierarchies import (
    pfqn_pbh, pfqn_cbh, pfqn_pbk, pfqn_bjbk, pfqn_ssd, pfqn_mcub,
    pfqn_sib, pfqn_ldbcmp, pfqn_scb, pfqn_looping,
)


def _inf_mask(s):
    from ...lang.base import SchedStrategy
    m = np.zeros(s.nstations, dtype=bool)
    sched = s.sched
    if sched is None:
        return m
    for i in range(s.nstations):
        sv = sched.get(i) if isinstance(sched, dict) else (
            sched[i] if i < len(sched) else None)
        if sv is None:
            continue
        m[i] = (sv == SchedStrategy.INF or
                (hasattr(sv, 'value') and sv.value == SchedStrategy.INF.value) or
                (isinstance(sv, int) and sv == SchedStrategy.INF.value))
    return m


def _sc_setup(s, method):
    # No single-class test here: ba_method_refusal owns that rule for every
    # caller, and a second copy is what would let the report gate and this run
    # disagree. solver_ba_analyzer has already asked it by the time we get here.
    rates = np.asarray(s.rates[:, 0], dtype=float)
    dem = np.asarray(s.demands[:, 0], dtype=float)   # = V/rate per station
    with np.errstate(divide='ignore', invalid='ignore'):
        V = dem * rates
    inf = _inf_mask(s)
    D = dem[~inf]
    Z = float(np.sum(dem[inf]))
    N = int(round(float(s.njobs[0])))
    return V, rates, inf, D, Z, N


def _level(options):
    lvl = getattr(options, 'level', None)
    if lvl is None and isinstance(getattr(options, '__dict__', None), dict):
        lvl = options.__dict__.get('level')
    return int(lvl) if lvl else 2


def _ba_fill(s, V, rates, inf, N, X, Z, D, is_upper):
    """Per-station [Q,U,R,T,C] from a scalar chain-throughput bound X."""
    M = s.nstations
    with np.errstate(divide='ignore', invalid='ignore'):
        T = (V * X).reshape(-1, 1)
        R = np.zeros((M, 1))
        if is_upper:
            R[:, 0] = (1.0 / rates) * N
            R[inf, 0] = 1.0 / rates[inf]
            C = Z + N * np.sum(D)
        else:
            R[:, 0] = 1.0 / rates
            C = Z + np.sum(D)
        Q = T * R
        # Utilization law per SERVER: without the nservers divisor a
        # multiserver station reports U > 1 (ssd/ldbcmp/auto reach here)
        srv = np.maximum(1.0, np.asarray(s.nservers, dtype=float).ravel())
        U = (T[:, 0] / (srv * rates)).reshape(-1, 1)
    U[inf, 0] = Q[inf, 0]
    return Q, U, R, T, np.array([[float(C)]]), np.array([X])


def _lr_bounds(s, method):
    """LP-based Linear Reduction bounds (Casale et al.).

    Solves the linear-reduction LP relaxation once per station, minimizing
    ('lr.lower') or maximizing ('lr.upper') that station's utilization. Both
    senses are valid bounds because the LP relaxation contains the exact
    solution. Unlike 'qrf.mmi.linear' -- whose name refers only to its explicit
    Aeq/beq constraint representation while its objective is the nonlinear MMI
    mutual information -- this method is a pure LP end to end (HiGHS via
    scipy.linprog).
    """
    from ...api.mapqn.parameters import PFParameters, LinearReductionParameters
    from ...api.mapqn.bnd_lr_pf import mapqn_bnd_lr_pf
    from ...api.mapqn.bnd_lr import mapqn_bnd_lr

    V, rates, inf, D, Z, N = _sc_setup(s, method)
    if np.any(inf):
        raise ValueError(
            "Method '%s' does not support delay (infinite-server) stations."
            % method)
    sense = 'min' if method.endswith('.lower') else 'max'

    M = s.nstations
    # Routing over stations, weighted by visit ratios (single class).
    r = np.zeros((M, M))
    tot = float(np.sum(V))
    for i in range(M):
        for j in range(M):
            r[i, j] = V[j] / tot if tot > 0 else 0.0

    params = PFParameters(_M=M, _N=N, mu=np.asarray(rates, dtype=float), r=r)
    U = np.zeros((M, 1))
    for i in range(M):
        U[i, 0] = mapqn_bnd_lr_pf(params, i + 1, sense).objective_value

    # Chain throughput implied by the bounded utilizations: U_i = X*V_i/rate_i.
    with np.errstate(divide='ignore', invalid='ignore'):
        cand = U[:, 0] * rates / V
    cand = cand[np.isfinite(cand)]
    X = float(np.min(cand)) if cand.size else 0.0
    with np.errstate(divide='ignore', invalid='ignore'):
        T = (V * X).reshape(-1, 1)
        R = np.zeros((M, 1))
        R[:, 0] = 1.0 / rates
        if sense == 'max':
            R[:, 0] = (1.0 / rates) * N
        Q = T * R
        C = Z + (N if sense == 'max' else 1) * np.sum(D)
    return Q, U, R, T, np.array([[float(C)]]), np.array([X])


def _mapamva_bounds(s, method):
    """MAP-AMVA LP bounds (Casale-Smirni, DSN 2009).

    The linear program over the EXACT mean-value balance equations of a closed
    MAP queueing network, solved through mapqn_bnd_lr_mva. It is the only family
    in SolverBA that consumes the CORRELATION between successive services rather
    than the service mean alone: its variables are the per-phase queue lengths
    QN(i,k) and utilizations UN(i,k), so a workload whose burstiness moves the
    bottleneck between stations is bounded rather than averaged into a renewal
    process. That is why 'MAP' and 'MMPP2' reach this family's feature set and
    no other.

    The LP carries phases at ONE queue and requires it to be the LAST -- q(i,j,
    k,h) reads the scalar muM(i) for i < M and the (D0,D1) pair muMAP/v for
    i == M -- so a model whose phase-carrying station sits elsewhere is PERMUTED
    rather than refused, and the results are permuted back before they are
    returned.
    """
    from ...api.mapqn.parameters import MVAVersionParameters
    from ...api.mapqn.bnd_lr_mva import mapqn_bnd_lr_mva
    from ...api.solvers.ctmc.solver_ctmc_qrf_analyzer import (
        _proc_entry_to_map, _map_mean_from_ph)

    V, rates, inf, D, Z, N = _sc_setup(s, method)
    if np.any(inf):
        raise ValueError(
            "Method '%s' does not support delay (infinite-server) stations: the "
            "MAP-AMVA program of Casale-Smirni (DSN 2009) is written for a network "
            "of queues and the paper names the delay extension as open work. Use a "
            "QRF method, which carries the load-dependent rate law." % method)

    sn = s._sn
    M = s.nstations
    upper = method.endswith('.upper')
    sense = 'max' if upper else 'min'

    # Phase order per station. One phase is an exponential server, which enters
    # the LP as the scalar rate muM(i); more than one is the (D0,D1) pair, which
    # only queue M can hold.
    maps = [_proc_entry_to_map(sn.proc[i]) for i in range(M)]
    kph = np.ones(M, dtype=int)
    for i in range(M):
        if maps[i] is not None:
            kph[i] = int(np.atleast_2d(maps[i][0]).shape[0])
    phased = [i for i in range(M) if kph[i] > 1]
    if len(phased) > 1:
        raise ValueError(
            "Method '%s' carries phases at ONE station: the LP gives queue M the "
            "(D0,D1) pair and every other queue a scalar rate. Stations %s are all "
            "non-exponential. Use a QRF method, whose q carries a phase at every "
            "station." % (method, [i + 1 for i in phased]))
    # Every station exponential: the program is still the right one, it just
    # degenerates to K = 1, where the per-phase variables collapse and the
    # balances become the product-form ones of mapqn_bnd_lr_pf.
    map_idx = phased[0] if phased else M - 1
    perm = [i for i in range(M) if i != map_idx] + [map_idx]
    K = int(kph[map_idx])

    if maps[map_idx] is None:
        D0 = np.array([[-1.0]])
        D1 = np.array([[1.0]])
    else:
        D0 = np.atleast_2d(np.asarray(maps[map_idx][0], dtype=float))
        D1 = np.atleast_2d(np.asarray(maps[map_idx][1], dtype=float))
    # muMAP(k,h) is the completion rate out of phase k landing in phase h, i.e.
    # D1(k,h); v(k,h) is the background phase change that completes no job, i.e.
    # D0 off the diagonal. Same (from,to) convention as the QRF adapter --
    # writing either as its transpose is invisible for a reversible D0 and
    # silently reverses the phase order of an Erlang.
    v = np.array(D0, dtype=float, copy=True)
    np.fill_diagonal(v, 0.0)

    mu_m = np.array([rates[perm[a]] for a in range(M - 1)], dtype=float)
    rt_raw = np.asarray(sn.rt, dtype=float)
    r = np.zeros((M, M))
    for a in range(M):
        for b in range(M):
            r[a, b] = rt_raw[perm[a], perm[b]]

    stimes = np.zeros(M)
    for i in range(M):
        if maps[i] is not None:
            stimes[i] = _map_mean_from_ph(sn.proc[i])
    Vp = np.array([V[perm[a]] for a in range(M)], dtype=float)
    Sp = np.array([stimes[perm[a]] for a in range(M)], dtype=float)

    params = MVAVersionParameters(_M=M, _N=int(N), K=K, muM=mu_m,
                                  muMAP=D1, r=r, v=v)

    # THREE SWEEPS OF THE SAME LP, and the utilization one runs in BOTH senses on
    # purpose. R_i = Q_i/(V_i*X) rises with Q_i and FALLS with X, so an upper
    # bound on the response time pairs Q_i^max with X^min; dividing by X^max on
    # both sides is what would report an upper R below the exact value and break
    # the bracket.
    umax = np.zeros(M)
    umin = np.zeros(M)
    qbnd = np.zeros(M)
    for ti in range(1, M + 1):
        umax[ti - 1] = mapqn_bnd_lr_mva(params, ti, 0, 'max', 'UN').objective_value
        umin[ti - 1] = mapqn_bnd_lr_mva(params, ti, 0, 'min', 'UN').objective_value
        qbnd[ti - 1] = mapqn_bnd_lr_mva(params, ti, 0, sense, 'QN').objective_value

    # Utilization law U_i = X*V_i*S_i, exact at a single server under ANY service
    # law, so each station turns its own utilization bound into a throughput
    # bound and the tightest of the M survives. A station with no visits or no
    # service time carries no information and is skipped.
    load_p = Vp * Sp
    ok = np.isfinite(load_p) & (load_p > 0)
    if not np.any(ok):
        raise ValueError(
            "Method '%s' found no station with both a positive visit ratio and a "
            "positive mean service time." % method)
    x_up = float(np.min(umax[ok] / load_p[ok]))
    x_lo = float(np.max(umin[ok] / load_p[ok]))
    xb, xopp = (x_up, x_lo) if upper else (x_lo, x_up)

    # Unpermute: the LP orders the stations with the phase-carrying one last.
    inv_perm = np.argsort(perm)
    Vs = Vp[inv_perm]
    Qs = qbnd[inv_perm]
    Us = (umax if upper else umin)[inv_perm]

    T = (Vs * xb).reshape(-1, 1)
    U = Us.reshape(-1, 1)
    Q = Qs.reshape(-1, 1)
    if xopp > 0:
        with np.errstate(divide='ignore', invalid='ignore'):
            R = (Qs / (Vs * xopp)).reshape(-1, 1)
        # Delay stations are refused above, so the closed-network response time
        # is N/X exactly and the throughput bracket transfers to it directly.
        # Summing the per-station R bounds instead would add M separately
        # attained maxima and report a looser number.
        C = float(N) / xopp
    else:
        R = np.full((M, 1), np.inf)
        C = np.inf
    return Q, U, R, T, np.array([[float(C)]]), np.array([xb])


def _asym_bounds(s, method):
    """ABA/BJB/PB/SB/GB noniterative bounds (explicit formulas)."""
    from ...lang.base import SchedStrategy  # noqa: F401
    V, rates, inf, D, Z, N = _sc_setup(s, method)
    M = s.nstations
    Dmax = np.max(D)
    sumD = np.sum(D)
    Q = np.zeros((M, 1)); U = np.zeros((M, 1)); R = np.zeros((M, 1))
    T = np.zeros((M, 1)); C = np.zeros((1, 1)); X = 0.0

    def fill(Xval, upper, Cval, Qvec=None, Rvec=None):
        Tv = V * Xval
        if Rvec is None:
            if upper:
                Rv = (1.0 / rates) * N
                Rv[inf] = 1.0 / rates[inf]
            else:
                Rv = 1.0 / rates
        else:
            Rv = Rvec
        Qv = Tv * Rv if Qvec is None else Qvec
        Uv = Tv / rates
        Uv[inf] = Qv[inf]
        return Qv.reshape(-1, 1), Uv.reshape(-1, 1), Rv.reshape(-1, 1), \
            Tv.reshape(-1, 1), np.array([[float(Cval)]]), np.array([float(Xval)])

    if method == 'aba.upper':
        X = min(1.0 / Dmax, N / (Z + sumD)); Cv = Z + N * sumD
        return fill(X, True, Cv)
    if method == 'aba.lower':
        X = N / (Z + N * sumD); Cv = Z + sumD
        return fill(X, False, Cv)
    if method in ('bjb.upper', 'bjb.lower'):
        Xau = min(1.0 / Dmax, (N - 1) / (Z + sumD))
        Xal = (N - 1) / (Z + (N - 1) * sumD)
        if method == 'bjb.upper':
            Cv = Z + sumD + Dmax * (N - 1 - Z * Xal)
            X = min(1.0 / Dmax, N / (Z + sumD + np.mean(D) * (N - 1 - Z * Xau)))
            return fill(X, True, Cv)
        Cv = Z + sumD + np.mean(D) * (N - 1 - Z * Xau)
        X = N / (Z + sumD + Dmax * (N - 1 - Z * Xal))
        return fill(X, False, Cv)
    if method in ('pb.upper', 'pb.lower'):
        Xau = min(1.0 / Dmax, (N - 1) / (Z + sumD))
        Xal = (N - 1) / (Z + (N - 1) * sumD)
        Dpb2 = np.sum(D ** 2) / sumD
        DpbN = np.sum(D ** N) / np.sum(D ** (N - 1))
        if method == 'pb.upper':
            Cv = Z + sumD + DpbN * (N - 1 - Z * Xal)
            X = min(1.0 / Dmax, N / (Z + sumD + Dpb2 * (N - 1 - Z * Xau)))
            return fill(X, True, Cv)
        Cv = Z + sumD + Dpb2 * (N - 1 - Z * Xau)
        X = N / (Z + sumD + DpbN * (N - 1 - Z * Xal))
        return fill(X, False, Cv)
    if method in ('sb.upper', 'sb.lower'):
        if np.any(inf):
            raise ValueError(
                "Method '%s' does not support infinite-server stations." % method)
        A1 = np.sum(D); A2 = np.sum(D ** 2); A3 = np.sum(D ** 3)
        if method == 'sb.upper':
            # Harel UB(n) is defined for n <= N only; the level-3 coefficient is
            # not a bound at N < 3, so fall back to UB(2) (Dallery), exact there
            # since UB(N) = TH(N).
            cub3 = (A1 * A2 + A3) / (A1 ** 2 + A2) if N >= 3 else A2 / A1
            Cv = Z + A1 + (N - 1) * cub3
            X = min(1.0 / Dmax, N / Cv)
            # Upper side: R is undefined in the literature for this bound, so
            # it carries the ABA PESSIMISTIC residence. Filling it from the
            # optimistic side put Q = T*R below exact on 7/24 chaos models.
            return fill(X, True, Cv)
        AN = np.sum(D ** N)
        # (N-1)*(AN/A1)^(1/(N-1)) -> 0 as N -> 1, leaving the exact single-job
        # cycle time; evaluated directly it divides by zero
        cterm = 0.0 if N == 1 else (N - 1) * (AN / A1) ** (1.0 / (N - 1))
        Cv = Z + A1 + cterm
        X = N / Cv
        return fill(X, False, Cv)
    if method in ('gb.upper', 'gb.lower'):
        if method == 'gb.upper':
            X = min(1.0 / Dmax, pfqn_xzgsbup(D, N, Z))
            Cv = N / pfqn_xzgsblow(D, N, Z)
            Xden = pfqn_xzgsblow(D, N, Z)
            qfun, upper = pfqn_qzgbup, True
        else:
            X = pfqn_xzgsblow(D, N, Z)
            Cv = N / pfqn_xzgsbup(D, N, Z)
            Xden = pfqn_xzgsbup(D, N, Z)
            qfun, upper = pfqn_qzgblow, False
        Tv = V * X
        Qv = np.zeros(M); Rv = np.zeros(M)
        k = -1   # 0-based queueing-station index for pfqn_qzgb*
        for i in range(M):
            if inf[i]:
                # Tv[i], not X: the delay is visited V[i] times per cycle, and
                # dropping that factor lets Q exceed the population
                Rv[i] = 1.0 / rates[i]; Qv[i] = Tv[i] * Rv[i]
            else:
                k += 1
                Qv[i] = qfun(D, N, Z, k)
                Rv[i] = Qv[i] / Xden / V[i] if V[i] != 0 else 0.0
        if upper:
            Rv[inf] = 1.0 / rates[inf]
        Uv = Tv / rates; Uv[inf] = Qv[inf]
        return (Qv.reshape(-1, 1), Uv.reshape(-1, 1), Rv.reshape(-1, 1),
                Tv.reshape(-1, 1), np.array([[float(Cv)]]), np.array([float(X)]))
    raise ValueError("Unknown asymptotic bound '%s'" % method)


def _chain_bound(s, method):
    """mwba.* and cub.upper/mbjb.lower (multiclass, chain-based)."""
    from ...api.sn.demands import sn_get_demands_chain
    from ...api.sn.deaggregate import sn_deaggregate_chain_results
    from ...lang.base import SchedStrategy
    sn = s._sn
    inf = _inf_mask(s)
    dc = sn_get_demands_chain(sn)
    Lchain, STchain, Vchain, alpha, Nchain = (
        dc.Lchain, dc.STchain, dc.Vchain, dc.alpha, dc.Nchain)
    M = s.nstations
    nchains = Lchain.shape[1]
    Zc = np.sum(Lchain[inf, :], axis=0)
    Lq = Lchain[~inf, :]

    if method in ('looping.upper', 'looping.lower'):
        # Eager Looping: the multiclass bracket that initializes the
        # multiple-class PBH. Pessimistic side from the heap-inflated response
        # time, optimistic side from the response-time lower bound.
        Xlo, Xup = pfqn_looping(Lq, Nchain.flatten(), Zc)[:2]
        up = (method == 'looping.upper')
        Xchain = np.asarray(Xup if up else Xlo, dtype=float).flatten()
    elif method in ('cub.upper', 'mbjb.lower'):
        Xub, Xlb = pfqn_mcub(Lq, Nchain.flatten(), Zc)
        up = (method == 'cub.upper')
        Xchain = np.asarray(Xub if up else Xlb, dtype=float).flatten()
    else:  # mwba.*
        Vq = Vchain[~inf, :]
        Sq = STchain[~inf, :]
        sched = s.sched
        schedq = []
        for i in range(M):
            if inf[i]:
                continue
            sv = sched.get(i) if isinstance(sched, dict) else sched[i]
            schedq.append(_mwrbb_disc_code(sv))
        prioc = np.zeros(nchains)
        Xlo, Xup, Wlo = pfqn_mwrbb(Vq, Sq, Nchain.flatten(), Zc,
                                   np.asarray(schedq, dtype=float), prioc)
        up = (method == 'mwba.upper')
        Xchain = np.asarray(Xup if up else Xlo, dtype=float).flatten()

    Tchain = np.zeros((M, nchains)); Uchain = np.zeros((M, nchains))
    for c in range(nchains):
        Tchain[:, c] = Xchain[c] * Vchain[:, c]
        Uchain[:, c] = Xchain[c] * Lchain[:, c]
    Qchain = _chain_qfill(Uchain, Xchain, Lchain, Nchain.flatten(), inf, up)
    deagg = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha, Qchain, Uchain, None,
        Tchain, None, Xchain)
    Cn = np.zeros((M, s.nclasses))
    return deagg.Q, deagg.U, deagg.R, deagg.T, Cn, np.asarray(deagg.X).flatten()


def _chain_qfill(Uchain, Xchain, Lchain, Nchain, isdelay, is_upper):
    """Per-chain queue lengths from a chain-throughput bound, on the declared
    side. Lower: Q_ic >= U_ic, since the station holds a class-c job whenever
    it serves one. Upper: E[n_ic] <= N_c*P(station busy) = N_c*min(1,sum_c
    U_ic), and a delay station queues nothing, so there Q_ic = X_c*L_ic."""
    M, C = Uchain.shape
    Q = np.zeros((M, C))
    Utot = np.minimum(1.0, np.sum(Uchain, axis=1))
    for c in range(C):
        for i in range(M):
            if isdelay[i]:
                Q[i, c] = Xchain[c] * Lchain[i, c]
            elif is_upper:
                Q[i, c] = Nchain[c] * Utot[i]
            else:
                Q[i, c] = Uchain[i, c]
    return Q


def _mwrbb_disc_code(sv):
    from ...lang.base import SchedStrategy
    val = sv.value if hasattr(sv, 'value') else sv
    if val == SchedStrategy.FCFS.value:
        return 0
    if val in (SchedStrategy.PS.value, getattr(SchedStrategy, 'DPS', SchedStrategy.PS).value,
               getattr(SchedStrategy, 'GPS', SchedStrategy.PS).value):
        return 1
    if val == getattr(SchedStrategy, 'HOL', SchedStrategy.FCFS).value:
        return 2
    return 4


def _hier_bound(s, method, options):
    """Level-parameterized single-class hierarchies + ssd/sib/ldbcmp."""
    V, rates, inf, D, Z, N = _sc_setup(s, method)
    lvl = _level(options)
    fam = method.split('.')[0]
    up = method.endswith('.upper')
    if fam == 'pbh':
        Xlo, Xhi = pfqn_pbh(D, N, Z, lvl)[:2]
    elif fam == 'cbh':
        Xlo, Xhi = pfqn_cbh(D, N, Z, lvl)
    elif fam == 'pbk':
        Xlo, Xhi = pfqn_pbk(D, N, Z, lvl)
    elif fam == 'bjbk':
        Xlo, Xhi = pfqn_bjbk(D, N, Z, lvl)
    elif fam == 'ssd':
        nsv = np.asarray(s.nservers)[~inf]
        Xlo, Xhi = pfqn_ssd(D, N, Z, nsv)
    elif fam == 'sib':
        if Z > 0:
            raise ValueError(
                "Method '%s' supports Z=0 (no delay station) only; delay needs "
                "the SIB Section-3.2 extension." % method)
        Xlo, Xhi = pfqn_sib(D, N, 0.0, lvl)[:2]
    elif fam == 'scb':
        # Single-class bounds of Dowdy et al. (1992). THE BRACKETED OBJECT IS
        # NOT THIS MODEL: scb brackets the multiclass system that this
        # single-class model aggregates, so scb.lower is the EXACT single-class
        # throughput and scb.upper adds the demand-free Expression-(3) gap.
        # That is why scb is absent from BA_AUTO_* -- mixing it with families
        # that bracket this model's own solution would compare two different
        # quantities.
        if Z > 0:
            raise ValueError(
                "Method '%s' supports Z=0 (no delay station) only; Theorem 3 "
                "rests on the delay-free balanced-network throughput." % method)
        Xlo, Xhi = pfqn_scb(D, N)[:2]
    elif fam == 'ldbcmp':
        Xlo, Rhi, Qhat = pfqn_ldbcmp(D, N, Z, np.zeros(len(D)))
        if N < Qhat:
            raise ValueError(
                "Method '%s' requires the asymptotic regime N >= Qhat "
                "(Qhat=%.4f > N=%d)." % (method, Qhat, N))
        return _ba_fill(s, V, rates, inf, N, float(Xlo), Z, D, False)
    else:
        raise ValueError("Unknown hierarchy bound '%s'" % method)
    Xval = float(Xhi) if up else float(Xlo)
    return _ba_fill(s, V, rates, inf, N, Xval, Z, D, up)


# Candidate list for the AUTO composite. Noniterative families only: the
# level-parameterized hierarchies (pbh/cbh/pbk/bjbk/sib) and the LP reductions
# are excluded because their cost is not O(K) and their accuracy is a user
# choice, not a fixed property.
BA_AUTO_UPPER = ('aba.upper', 'bjb.upper', 'pb.upper', 'gb.upper', 'sb.upper',
                 'mwba.upper', 'ssd.upper', 'cub.upper')
BA_AUTO_LOWER = ('aba.lower', 'bjb.lower', 'pb.lower', 'gb.lower', 'sb.lower',
                 'mwba.lower', 'ssd.lower', 'mbjb.lower', 'ldbcmp.lower')
BA_AUTO = {'auto.upper', 'auto.lower'}

BA_ASYM = {'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower',
           'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower',
           'gb.upper', 'gb.lower'}
BA_CHAIN = {'mwba.upper', 'mwba.lower', 'cub.upper', 'mbjb.lower',
            'looping.upper', 'looping.lower'}
BA_LR = {'lr', 'lr.upper', 'lr.lower'}
BA_MAPAMVA = {'mapamva.upper', 'mapamva.lower'}
BA_HAREL = {'harel.upper', 'harel.lower'}

# The two LP-backed QRF reduction bounds. These are served by the stable
# qr_bounds_bas / qr_bounds_rsrd simplex backends (via scipy HiGHS), NOT by the
# unstable qrf_noblo_* NLP library, so unlike the other seven qrf method names they
# are advertised and dispatched natively.
BA_QRF_LP = {'qrf.bas', 'qrf.rsrd', 'qrf.bas.mmi', 'qrf.bas.mem',
             'qrf.bas.bethe'}

# The NLP-backed QRF method names, served by api.mapqn.qrf_noblo_*. Withheld until
# 2026-07-20 because the SolverBA -> solver_ctmc_qrf_analyzer -> qrf_noblo_*
# path had never been exercised end to end. Exercising it found and fixed a
# transposed v in the adapter and a throughput inversion that assumed the
# reference station was a delay; both codebases now derive UN from the QRF
# utilization directly. On a 2-station closed chain the K == 1 result is exact
# against the CTMC and the K == 2 result lands inside the glpsol bound range.
BA_QRF_NLP = {'qr', 'qrf.mmi', 'qrf.mem', 'qrf.bethe', 'qrf.mmi.ld',
              'qrf.mmi.linear'}


def _qrf_lp_bounds(s, method, options):
    """QRF reduction bounds, both the LP-backed and the NLP-backed method names.

    Delegates to the shared QRF adapter, which builds the per-station MAP
    (D0,D1) pairs and routing from the NetworkStruct, calls the LP backend and
    reconstructs the metrics via Little's law and the visit ratios.
    """
    from ...api.solvers.ctmc.solver_ctmc_qrf_analyzer import solver_ctmc_qrf_analyzer

    sn = getattr(s, '_sn', None)
    if sn is None and hasattr(s, 'model'):
        model = s.model
        sn = model.get_struct() if hasattr(model, 'get_struct') else model.getStruct()
    if sn is None:
        raise RuntimeError("SolverBA: no NetworkStruct available for method '%s'." % method)

    class _Opts(object):
        pass

    o = _Opts()
    o.method = method
    o.config = getattr(options, 'config', None) or {}
    QN, UN, RN, TN, CN, XN, _rt = solver_ctmc_qrf_analyzer(sn, o)
    return QN, UN, RN, TN, CN, XN

BA_BPT = {'bpt.lower'}


def _bpt_bound(s, method):
    """Achievable-region LOWER bound for a multiclass OPEN Markovian network.

    Native port of matlab/src/solvers/BA/solver_ba_bpt_analyzer.m, cross-checked
    against jar/src/main/java/jline/solvers/ba/analyzers/Solver_ba_bpt_analyzer.java.

    CLASS SPACE. The reference's "class" is a buffer: one exponential service
    rate, one Markovian routing law. LINE's (station, job class) pair is exactly
    that, so a pair carrying traffic becomes one LP class, the Source is
    absorbed into the external arrival vector, and class switching needs no
    special treatment because sn.rt already carries it.

    BOUND CONVENTION. R(i,r) minimizes x over the polyhedron with the objective
    set to that pair's unit vector, so each entry is a valid lower bound on its
    own. Q follows by Little's law from the bounded R and the EXACT throughput T
    (an open network's per-class rates are fixed by the traffic equations, not
    by the policy), and so does C. U is exact for the same reason.

    TIGHTNESS. Exact on M/M/1 and tight on the externally fed classes, but weak
    on a class whose arrivals are all internal: the only term coupling x_r to
    the second-moment block carries the factor lambda0_r, so an internally fed
    class can fall back to its own mean service time.
    """
    from ...api.npfqn import npfqn_bnd_bpt
    from ...api.sn.transforms import sn_rt_stations
    from ...api.sn.network_struct import NodeType
    from ...constants import ProcessType, SchedStrategy

    sn = getattr(s, '_sn', None)
    if sn is None and hasattr(s, 'model'):
        model = s.model
        sn = model.get_struct() if hasattr(model, 'get_struct') else model.getStruct()
    if sn is None:
        raise RuntimeError("SolverBA: no NetworkStruct available for method '%s'." % method)

    M = int(sn.nstations)
    K = int(sn.nclasses)
    Q = np.zeros((M, K)); U = np.zeros((M, K))
    R = np.zeros((M, K)); T = np.zeros((M, K))
    C = np.zeros((1, K)); X = np.zeros((1, K))

    # ----- model gates -----
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isfinite(njobs)):
        raise ValueError("Method 'bpt.lower' supports fully open networks only "
                         "(no closed classes).")
    sched_dict = sn.sched if sn.sched else {}
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    isSource = np.zeros(M, dtype=bool)
    for i in range(M):
        isSource[i] = (sn.nodetype[int(sn.stationToNode[i])] == NodeType.SOURCE)
    srcList = np.where(isSource)[0]
    qstat = np.where(~isSource)[0]
    if srcList.size == 0:
        raise ValueError("Method 'bpt.lower' requires an open network with a Source station.")
    for i in qstat:
        if sched_dict.get(int(i), SchedStrategy.FCFS) == SchedStrategy.INF:
            raise ValueError("Method 'bpt.lower' does not support delay "
                             "(infinite-server) stations: the achievable region is "
                             "derived for one server per station.")
        if nservers[int(i)] > 1:
            raise ValueError("Method 'bpt.lower' does not support multi-server stations.")

    # ----- station-space routing, with the Source absorbed into lambda0 -----
    rtst = np.asarray(sn_rt_stations(sn)[0], dtype=float)
    rates = np.asarray(sn.rates, dtype=float)

    pair_station = []
    pair_class = []
    pair_flat = []
    for i in qstat:
        for r in range(K):
            pair_station.append(int(i))
            pair_class.append(r)
            pair_flat.append(int(i) * K + r)
    np_pairs = len(pair_flat)
    pair_station = np.asarray(pair_station, dtype=int)
    pair_class = np.asarray(pair_class, dtype=int)
    pair_flat = np.asarray(pair_flat, dtype=int)

    lambda0 = np.zeros(np_pairs)
    for si in srcList:
        for r0 in range(K):
            arr = rates[int(si), r0]
            if not np.isfinite(arr) or arr <= 0:
                continue
            lambda0 += arr * rtst[int(si) * K + r0, pair_flat]

    # Flow to the Sink or back to a Source is the exit probability, i.e. the
    # row deficit, and needs no column.
    P = rtst[np.ix_(pair_flat, pair_flat)]

    # ----- restrict to the pairs that actually carry traffic -----
    lam_all = np.linalg.solve(np.eye(np_pairs) - P.T, lambda0)
    keep = np.where(lam_all > 1e-12 * max(1.0, float(lam_all.max())))[0]
    if keep.size == 0:
        raise ValueError("The model carries no open traffic.")
    lambda0 = lambda0[keep]
    P = P[np.ix_(keep, keep)]
    pair_station = pair_station[keep]
    pair_class = pair_class[keep]
    nk = keep.size

    mu = np.zeros(nk)
    for a in range(nk):
        mu[a] = rates[pair_station[a], pair_class[a]]
        if not np.isfinite(mu[a]) or mu[a] <= 0:
            raise ValueError("Station %d has no service rate for class %d but carries "
                             "its traffic." % (pair_station[a] + 1, pair_class[a] + 1))
        pid = sn.procid[pair_station[a], pair_class[a]]
        if pid != ProcessType.EXP:
            raise ValueError("Method 'bpt.lower' requires exponential service: station %d "
                             "class %d is %s." % (pair_station[a] + 1, pair_class[a] + 1, pid))

    # Dense station index space for the LP.
    ustat, station_of = np.unique(pair_station, return_inverse=True)

    # ----- one LP per pair, objective = that pair's unit vector -----
    for a in range(nk):
        e = np.zeros(nk)
        e[a] = 1.0
        info = npfqn_bnd_bpt(lambda0, mu, P, station_of, e)
        R[pair_station[a], pair_class[a]] = info.zlb
        T[pair_station[a], pair_class[a]] = info.lambda_[a]
        U[pair_station[a], pair_class[a]] = info.rho[a]

    # ----- exact open-network quantities -----
    for si in srcList:
        for r in range(K):
            arr = rates[int(si), r]
            if np.isfinite(arr) and arr > 0:
                T[int(si), r] += arr
                X[0, r] += arr
    Q = T * R
    for r in range(K):
        if X[0, r] > 0:
            C[0, r] = float(np.sum(Q[:, r])) / X[0, r]
    return Q, U, R, T, C, X



BA_SNC = {'snc.upper'}


def _snc_bound(s, method):
    """Stochastic network calculus UPPER bound for a feed-forward OPEN network.

    Native port of matlab/src/solvers/BA/solver_ba_snc_analyzer.m. The envelope
    algebra lives in line_solver.api.snc and the model mapping in
    solver_ba_snc.py; this is only the analyzer entry point.

    UNITS ARE JOBS, NOT WORK, which is what lets a departure envelope be the
    arrival envelope of the next hop. R is the integral of the delay tail bound,
    Q its Little's-law image on the exact open throughput, U and T exact.
    """
    from .solver_ba_snc import snc_bound

    sn = getattr(s, '_sn', None)
    if sn is None and hasattr(s, 'model'):
        model = s.model
        sn = model.get_struct() if hasattr(model, 'get_struct') else model.getStruct()
    if sn is None:
        raise RuntimeError("SolverBA: no NetworkStruct available for method '%s'." % method)
    return snc_bound(sn)


BA_SPNLP = {'spnlp.upper', 'spnlp.lower', 'spnlp.op.upper', 'spnlp.op.lower'}


def _spnlp_bound(s, method, options):
    """Moment-relaxation LP bounds for a stochastic Petri net (Liu 1998).

    Native port of matlab/src/solvers/BA/solver_ba_spnlp_analyzer.m. The
    polytope lives in api.spn.spn_lpbnd and the model mapping in
    solver_ba_spnlp.py; this is only the analyzer entry point.

    The only family here indexed by a MARKING rather than by demands and a
    population, so it is offered on a Petri net and nowhere else.
    """
    from .solver_ba_spnlp import spnlp_bound

    sn = getattr(s, '_sn', None)
    if sn is None and hasattr(s, 'model'):
        model = s.model
        sn = model.get_struct() if hasattr(model, 'get_struct') else model.getStruct()
    if sn is None:
        raise RuntimeError("SolverBA: no NetworkStruct available for method '%s'." % method)
    return spnlp_bound(sn, method, options)


BA_BGT = {'bgt.upper'}


def _bgt_single_successor(rtst, row, pair_flat, who):
    """The single successor of a routing row, as a pair index, or -1 when
    everything leaves the network.

    A probabilistic split is refused by name: the reference's network has
    deterministic routing and a split is a different model, not an
    approximation of this one.
    """
    tol = 1e-9
    v = rtst[row, pair_flat]
    mass = float(v[v > tol].sum())
    if mass <= tol:
        return -1
    p = int(np.argmax(v))
    best = float(v[p])
    if abs(mass - 1) > tol or abs(best - 1) > tol:
        raise ValueError(
            "Method 'bgt.upper' needs deterministic routing: %s splits its departures "
            "(the largest branch carries %.6g of them). The reference's network routes "
            "each type along a fixed sequence of stages." % (who, best))
    return p


def _bgt_bound(s, method):
    """Piecewise-linear Lyapunov UPPER bound for a multitype OPEN network.

    Native port of matlab/src/solvers/BA/solver_ba_bgt_analyzer.m, cross-checked
    against jar/src/main/java/jline/solvers/ba/analyzers/Solver_ba_bgt_analyzer.java.

    CLASS SPACE. The reference's network is a MULTITYPE one: each type follows a
    FIXED sequence of stages, and stage k of type i is its own buffer. LINE's
    (station, job class) pair is that buffer, so the analyzer walks the routing
    matrix from the Source and turns each open class into one type whose stages
    are the pairs it visits. Routing must be DETERMINISTIC and routes must not
    MERGE; both are refused by name rather than approximated. A re-entrant line
    is expressible by giving the revisits distinct LINE classes.

    THE BOUND IS LOOSE, and knowingly so: the exception parameter of the
    smoothed Lyapunov function carries (Lmax+gamma)^3/gamma^2 and dominates as
    soon as there is more than one station. What is sharp is the STABILITY
    CERTIFICATE and the geometric tail RATE.
    """
    from ...api.npfqn import npfqn_bnd_bgt
    from ...api.sn.network_struct import NodeType
    from ...api.sn.transforms import sn_rt_stations
    from ...constants import ProcessType, SchedStrategy

    sn = getattr(s, '_sn', None)
    if sn is None and hasattr(s, 'model'):
        model = s.model
        sn = model.get_struct() if hasattr(model, 'get_struct') else model.getStruct()
    if sn is None:
        raise RuntimeError("SolverBA: no NetworkStruct available for method '%s'." % method)

    M = int(sn.nstations)
    K = int(sn.nclasses)
    Q = np.zeros((M, K)); U = np.zeros((M, K))
    R = np.zeros((M, K)); T = np.zeros((M, K))
    C = np.zeros((1, K)); X = np.zeros((1, K))

    # ----- model gates -----
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isfinite(njobs)):
        raise ValueError("Method 'bgt.upper' supports fully open networks only "
                         "(no closed classes).")
    sched_dict = sn.sched if sn.sched else {}
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    isSource = np.zeros(M, dtype=bool)
    for i in range(M):
        isSource[i] = (sn.nodetype[int(sn.stationToNode[i])] == NodeType.SOURCE)
    srcList = np.where(isSource)[0]
    qstat = np.where(~isSource)[0]
    if srcList.size == 0:
        raise ValueError("Method 'bgt.upper' requires an open network with a Source station.")
    for i in qstat:
        if sched_dict.get(int(i), SchedStrategy.FCFS) == SchedStrategy.INF:
            raise ValueError("Method 'bgt.upper' does not support delay (infinite-server) "
                             "stations: the reference's network has one server per station.")
        if nservers[int(i)] > 1:
            raise ValueError("Method 'bgt.upper' does not support multi-server stations.")

    rtst = np.asarray(sn_rt_stations(sn)[0], dtype=float)
    rates = np.asarray(sn.rates, dtype=float)

    pair_station = []
    pair_class = []
    pair_flat = []
    for i in qstat:
        for r in range(K):
            pair_station.append(int(i))
            pair_class.append(r)
            pair_flat.append(int(i) * K + r)
    pair_station = np.asarray(pair_station, dtype=int)
    pair_class = np.asarray(pair_class, dtype=int)
    pair_flat = np.asarray(pair_flat, dtype=int)
    npairs = pair_flat.size

    # ----- walk one deterministic route per source class -----
    lam = []
    routes = []
    used = np.zeros(npairs, dtype=bool)
    for si in srcList:
        for r0 in range(K):
            arr = rates[int(si), r0]
            if not np.isfinite(arr) or arr <= 0:
                continue
            cur = _bgt_single_successor(rtst, int(si) * K + r0, pair_flat,
                                        "the Source for class %d" % (r0 + 1))
            route = []
            while cur >= 0:
                if used[cur]:
                    raise ValueError(
                        "Method 'bgt.upper' needs routes that do not merge: station %d class %d "
                        "is visited by more than one type. Give the visits distinct job classes."
                        % (pair_station[cur] + 1, pair_class[cur] + 1))
                used[cur] = True
                route.append(cur)
                cur = _bgt_single_successor(rtst, int(pair_flat[cur]), pair_flat,
                                            "station %d class %d"
                                            % (pair_station[route[-1]] + 1,
                                               pair_class[route[-1]] + 1))
            if not route:
                raise ValueError("Class %d leaves the Source and reaches no station." % (r0 + 1))
            routes.append(route)
            lam.append(float(arr))
    if not routes:
        raise ValueError("The model carries no open traffic.")

    mu_l = []
    sigma_l = []
    ustat = []
    for route in routes:
        muv = []
        stv = []
        for p in route:
            ist, r = int(pair_station[p]), int(pair_class[p])
            m = rates[ist, r]
            if not np.isfinite(m) or m <= 0:
                raise ValueError("Station %d has no service rate for class %d but carries "
                                 "its traffic." % (ist + 1, r + 1))
            pid = sn.procid[ist, r]
            if pid != ProcessType.EXP:
                raise ValueError("Method 'bgt.upper' requires exponential service: station %d "
                                 "class %d is %s." % (ist + 1, r + 1, pid))
            muv.append(float(m))
            stv.append(ist)
            if ist not in ustat:
                ustat.append(ist)
        mu_l.append(muv)
        sigma_l.append(stv)
    # Dense station index space for the LP.
    sigma_l = [[ustat.index(x) for x in stv] for stv in sigma_l]

    info = npfqn_bnd_bgt(lam, mu_l, sigma_l, len(ustat))

    # ----- read the bound back per station and class -----
    for i, route in enumerate(routes):
        for k, p in enumerate(route):
            ist, r = int(pair_station[p]), int(pair_class[p])
            Q[ist, r] += float(info.Qub[i][k])
            T[ist, r] += lam[i]
            U[ist, r] += lam[i] / rates[ist, r]
    nz = T > 0
    R[nz] = Q[nz] / T[nz]

    # ----- exact open-network quantities -----
    for si in srcList:
        for r in range(K):
            arr = rates[int(si), r]
            if np.isfinite(arr) and arr > 0:
                T[int(si), r] += arr
                X[0, r] += arr
    for r in range(K):
        if X[0, r] > 0:
            C[0, r] = float(np.sum(Q[:, r])) / X[0, r]
    return Q, U, R, T, C, X


BA_HIER = {'pbh.upper', 'pbh.lower', 'cbh.upper', 'cbh.lower',
           'pbk.upper', 'pbk.lower', 'bjbk.upper', 'bjbk.lower',
           'ssd.upper', 'ssd.lower', 'sib.upper', 'sib.lower',
           'scb.upper', 'scb.lower', 'ldbcmp.lower'}


def _harel_bound(s, method):
    """Sharp bounds of Harel-Namn-Sturm, distinct from the 'sb' family of the
    same paper: they extrapolate from the EXACT normalizing constant at
    populations n <= 7 instead of using the first three power sums."""
    V, rates, inf, D, Z, N = _sc_setup(s, method)
    if Z > 0:
        raise ValueError(
            "Method '%s' does not support think times (infinite-server stations)." % method)
    max_ub = min(N, 7)
    LB, UB, TH = pfqn_harel_bounds(D, N, 0.0, max_ub)
    up = (method == 'harel.upper')
    if up:
        Xb = min(1.0 / float(np.max(D)), float(UB[max_ub - 1] if max_ub >= 2 else TH[0]))
    else:
        Xb = float(LB)
    return _ba_fill(s, V, rates, inf, N, Xb, Z, D, up)


def _auto_bound(s, method, options):
    """AUTO composite: evaluate every noniterative bound and keep the tightest
    side. Feasibility is probed by execution -- a candidate that rejects the
    model (multiserver, delay station, regime gate) raises and is skipped -- so
    the list stays correct as families are added."""
    V, rates, inf, D, Z, N = _sc_setup(s, method)
    up = method.endswith('.upper')
    best, Xbest = None, None
    for cand in (BA_AUTO_UPPER if up else BA_AUTO_LOWER):
        try:
            r = solver_ba_analyzer(s, cand, options)
        except Exception:
            continue
        Xc = np.asarray(r['XN'], dtype=float).flatten()
        if Xc.size == 0 or not np.isfinite(Xc[0]) or Xc[0] <= 0:
            continue
        if Xbest is None or (up and Xc[0] < Xbest) or (not up and Xc[0] > Xbest):
            Xbest, best = float(Xc[0]), cand
    if Xbest is None:
        raise ValueError(
            "Method '%s' found no feasible bound for this model." % method)
    return _ba_fill(s, V, rates, inf, N, Xbest, Z, D, up)


def ba_resolve_method(method):
    """The 'default'/'auto'/'qr'/'lr' aliases, as MATLAB runAnalyzer resolves
    them. Without 'qr' the method name reaches the QRF adapter unresolved and is
    rejected."""
    if method == 'default':
        return 'gb.upper'
    if method == 'auto':
        return 'auto.upper'
    if method == 'qr':
        return 'qrf.mmi'
    if method == 'lr':
        return 'lr.upper'
    return method


# The families whose bound is a function of the single-chain demand vector
# D = V/rates, the think time Z and the population N. A multiclass or open model
# simply does not have those, which is why the rule is total rather than a
# tolerance.
# 'mapamva' is single-class for a different reason from the rest -- its LP
# variables QN(i,k)/UN(i,k) are indexed by station and MAP phase, with no class
# index at all -- but the premise it fails on is the same one.
BA_SINGLE_CLASS_FAMS = ('auto', 'aba', 'bjb', 'pb', 'sb', 'gb', 'harel', 'lr',
                        'pbh', 'cbh', 'pbk', 'bjbk', 'ssd', 'sib', 'scb',
                        'ldbcmp', 'mapamva')
# The multiclass families: they take a per-chain demand MATRIX and a population
# VECTOR, so several classes are fine and an infinite population is not.
BA_FULLY_CLOSED_FAMS = ('mwba', 'cub', 'mbjb', 'looping')
# Single-server families whose alternative on a multiserver model IS 'ssd',
# which is why the reason names it. ssd is the multiserver bound itself,
# ldbcmp is parameterized by the limiting demand of a load-dependent station,
# and auto composes whichever candidates survive, so all three are absent.
BA_SINGLE_SERVER_FAMS = ('aba', 'bjb', 'pb', 'sb', 'gb', 'harel', 'lr',
                         'pbh', 'cbh', 'pbk', 'bjbk', 'sib', 'scb', 'mapamva')


def ba_method_refusal(sn, method):
    """The STRUCTURAL premises of the SolverBA bound families, in one place:
    the reason METHOD cannot bound the model SN, or '' when it can.

    ONE PREDICATE, TWO CALLERS. solver_ba_analyzer asks it before dispatching
    and raises what it returns; SolverBA.supportsModelMethod asks it after the
    feature gate and reports the same sentence, which is what findSolver,
    list_valid_methods and SolverAUTO's ranked choice all read. A second copy of
    any rule below is how the report and the run drift apart -- the report
    offers a pair that raises the moment it is run, which is the defect this
    function exists to remove.

    WHAT BELONGS HERE AND WHAT DOES NOT. Only the rules the feature registry
    cannot name. SolverFeatureSet.FIELDS has no entry for "one class", for a
    server count or for a station count, so those are structural and live here.
    Rules of the form "this family does not accept a delay station" or "does not
    accept an open class" ARE nameable and belong in
    SolverBA.getMethodFeatureSet, which drops SchedStrategy_INF / OpenClass from
    the offending method's set instead: a feature set can refuse a model for
    HAVING a construct, never for lacking one.

    METHOD is taken as the caller spells it and resolved through
    ba_resolve_method, so 'default' is judged as the gb.upper it runs as and the
    reason names that. The marking-parameterized spnlp and the QRF reduction
    bounds carry no rule here: the QRF premise is the reducibility test
    list_valid_methods already applies. Of the three OPEN families, 'bpt' and
    'bgt' carry none either -- a closed model is refused by their feature set
    and their analyzers walk the routing matrix for the rest -- while 'snc'
    carries one, the SERVICE law.

    WHY THE SNC SERVICE LAW IS HERE AND THE bpt/bgt ONE IS NOT. All three
    analyzers refuse a non-exponential law at a queueing station. For bpt and
    bgt that rule extends to the SOURCE and is registry-expressible, so it rides
    in SolverBA.getMethodFeatureSet as a dropped law: both are invariant to the
    arrival law beyond its mean (replacing the Exp(1) source of an M/M/1 by an
    Erlang of the same mean leaves bgt.upper at QLen 32.6667 and bpt.lower at
    1.0, digit for digit), so a non-exponential source is not something they
    refuse, it is something they silently bound as if it were Poisson. snc is
    the opposite: it CONSUMES the arrival law (the same substitution moves it
    from 3.8244 to 3.0092) and its analyzer branches on a non-exponential source
    deliberately. Its rule is about the SERVICE only, and no feature name can
    say "Erlang at a Queue but not at a Source", so it is structural.

    Mirrors matlab/src/solvers/BA/ba_method_refusal.m and its JAR and C++ twins.
    """
    from ...api.sn.network_struct import SchedStrategy

    resolved = ba_resolve_method(method)
    fam = resolved.split('.')[0]
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    nclosedjobs = float(np.sum(njobs[np.isfinite(njobs)]))
    if fam in BA_SINGLE_CLASS_FAMS:
        if int(sn.nclasses) != 1 or not (nclosedjobs > 0):
            return ("Method '%s' supports single-class closed networks only."
                    % resolved)
    elif fam in BA_FULLY_CLOSED_FAMS:
        if not (nclosedjobs > 0) or bool(np.any(np.isinf(njobs))):
            return "Method '%s' supports fully closed networks only." % resolved

    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    multiserver = False
    for i in range(int(sn.nstations)):
        if sn.sched[i] == SchedStrategy.INF:
            continue
        if nservers[i] > 1:
            multiserver = True
            break
    if multiserver:
        if fam in BA_SINGLE_SERVER_FAMS:
            return ("Method '%s' does not support multi-server stations "
                    "(use 'ssd')." % resolved)
        if fam in BA_FULLY_CLOSED_FAMS:
            return "Method '%s' does not support multi-server stations." % resolved

    # The SNC service law. Judged over the pairs a station COULD serve rather
    # than over the ones that carry traffic: the analyzer restricts to the
    # latter, which needs the traffic equations solved, and this is their
    # conservative outer approximation -- it never admits a pair the analyzer
    # refuses, and can only differ on a pair given a service time at a station
    # its class never visits. A Source is skipped, which is the whole reason
    # this is not a feature-set delta.
    if fam == 'snc':
        from ...api.sn.network_struct import NodeType
        from ...constants import ProcessType
        rates = np.asarray(sn.rates, dtype=float)
        for i in range(int(sn.nstations)):
            if sn.nodetype[int(sn.stationToNode[i])] == NodeType.SOURCE:
                continue
            for r in range(int(sn.nclasses)):
                mu = rates[i, r]
                if not np.isfinite(mu) or mu <= 0:
                    continue
                if sn.procid[i, r] != ProcessType.EXP:
                    return ("Method '%s' requires exponential service: station %d "
                            "class %d is %s." % (resolved, i + 1, r + 1, sn.procid[i, r]))
    return ''


def ba_method_degenerate(sn, method):
    """Whether METHOD APPLIES to SN but its bound carries no information there,
    and why. Empty when the bound is informative, and empty for every method
    that has no such regime.

    THIS IS A DIFFERENT QUESTION FROM ba_method_refusal, which is why it is a
    different function. That one answers "is this model outside the method's
    domain", and its answer is what the analyzer raises. This one answers
    "inside the domain, does the formula still say anything", and its answer is
    NOT raised: a degenerate bound is a VALID bound, just a vacuous one, so an
    analyzer asked for it by name is entitled to publish it. What must not
    happen is OFFERING it: findSolver and list_valid_methods exist to name the
    pairs a caller can act on, and a table of zeros over a network with jobs
    circulating in it is not something anyone can act on.

    THE ONE METHOD WITH SUCH A REGIME IS 'ldbcmp.lower'. The
    Anselmi-Cremonesi bound is built from the population SURPLUS a = N - Qhat,
    where Qhat is the occupancy the non-bottleneck stations and the think time
    would hold in the open network fed at the bottleneck's saturation rate.
    pfqn_ldbcmp returns NaN below the regime (a < 0), which the analyzer already
    refuses by name; AT the boundary a = 0 it returns Xlo = 0, which is formally
    the trivial bound X >= 0 and propagates into a table whose queue lengths,
    utilizations and throughputs are all zero. Every entry of that table is a
    true lower bound and none of them is usable, and a caller cannot tell it
    from a real answer of zero. So the bound is computed here and the name
    withheld when it degenerates.

    The cost is one pfqn_ldbcmp evaluation, a closed form plus a scalar fixed
    point, and only for the single method that has the regime -- every other
    name returns immediately.

    Mirrors matlab/src/solvers/BA/ba_method_degenerate.m and its JAR and C++
    twins.
    """
    from ...api.sn.network_struct import SchedStrategy

    if ba_resolve_method(method) != 'ldbcmp.lower':
        return ''
    # The applicability rules come first and are not restated: a model this
    # method is outside the domain of has no bound to be degenerate about.
    if ba_method_refusal(sn, method):
        return ''

    # sn.visits is chain-indexed and this model has one chain, since
    # ba_method_refusal has already established it is single-class closed.
    #
    # TWO INDEX SPACES MEET HERE, and conflating them made this predicate RAISE
    # on a Petri net. sn.visits is STATEFUL-indexed -- sn_refresh_visits builds
    # it (nstateful, nclasses) -- while sn.sched and sn.rates are
    # STATION-indexed. Every station is stateful, but not every stateful node is
    # a station: a Transition is stateful without being one, and a Place is
    # both. On the fork-join SPN that is 4 stations against 7 stateful nodes, so
    # `visits / rates` came apart with "operands could not be broadcast together
    # with shapes (7,) (4,)". stationToStateful is the converter, and it is what
    # sn_refresh_visits itself uses.
    visits_stateful = np.asarray(sn.visits[0], dtype=float).ravel()
    s2sf = np.asarray(sn.stationToStateful, dtype=int).ravel()
    nst = int(sn.nstations)
    visits = np.array([visits_stateful[s2sf[i]] for i in range(nst)], dtype=float)
    rates = np.asarray(sn.rates, dtype=float)[:, 0]
    inf = np.array([sn.sched[i] == SchedStrategy.INF for i in range(nst)], dtype=bool)
    with np.errstate(divide='ignore', invalid='ignore'):
        dem = visits / rates
    Z = float(np.sum(dem[inf]))
    D = dem[~inf]
    N = int(round(float(np.sum(np.asarray(sn.njobs, dtype=float)))))
    # NO QUEUEING STATION, NO BOUND. Every station is a delay (or a Place, on a
    # Petri net, which is an INF station too), so there is no bottleneck to
    # build Qhat on and pfqn_ldbcmp cannot take an empty demand vector. A
    # PREDICATE MUST NOT RAISE -- this one is asked once per name by
    # list_valid_methods, before the Petri sieve has had a chance to drop
    # anything -- so the case is answered rather than propagated.
    if D.size == 0:
        return ("Method 'ldbcmp.lower' has no queueing station to bound here: every station "
                "is an infinite server, so the bottleneck the open-network occupancy is "
                "built on does not exist.")
    Xlo, _Rhi, Qhat = pfqn_ldbcmp(D, N, Z, np.zeros(len(D)))

    if not np.isfinite(Qhat):
        return ("Method 'ldbcmp.lower' does not apply here: the open-network occupancy "
                "Qhat the bound is built from does not exist, because a non-bottleneck "
                "station saturates at the bottleneck's arrival rate.")
    if not np.isfinite(Xlo) or Xlo <= 0:
        return ("Method 'ldbcmp.lower' needs a population strictly above the open-network "
                "occupancy the bound is built from (Qhat=%.4f, N=%d): with no surplus it "
                "degenerates to the trivial bound X >= 0 and reports a table of zeros."
                % (Qhat, N))
    return ''


def ba_ignores_blocking(method):
    """Whether METHOD bounds a model as if its buffers were unbounded.

    Every family but the QRF BLOCKING bounds is parameterized by demands
    (visits x service time) and a population alone, which is the BCMP
    parameterization: unbounded buffers, and an equilibrium distribution that
    factorizes. A finite buffer that BINDS breaks both premises, so the numbers
    do not bracket the blocked model -- on cqn_bas_blocking (Queue2 capped at
    1, N = 2) gb.upper reports QLen 1.28 at a station that can never hold more
    than one job. 'qrf.bas*' and 'qrf.rsrd' carry the blocking tables (MM, MM1,
    ZZ, ZM, BB, F) explicitly and are the exceptions.

    METHOD must already be resolved through ba_resolve_method.
    """
    return not (method.startswith('qrf.bas') or method.startswith('qrf.rsrd'))


def ba_blocking_default(sn):
    """What 'default'/'auto' must mean on a model with finite buffers.

    Returns ``(method, why)``. 'default' resolves to the geometric upper bound,
    which is parameterized by demands and a population alone and therefore
    bounds a blocked model as if its buffers were unbounded; ba_ignores_blocking
    refuses that, which is right. But refusing is not the whole answer, because
    'qrf.bas' DOES model the finite buffer and sn_to_qrf_blocking derives its
    tables from the model, so there is nothing left for the caller to supply. A
    blocked model of the right shape therefore gets 'qrf.bas' as its default,
    exactly as SolverMVA routes a BAS model to 'sqd'.

    The shape is the one list_valid_methods calls "reducible" and the QRF
    analyzer gates on: single class, closed, no delay station, no multiserver
    station. On top of that the tables must actually derive, which is asked of
    sn_to_qrf_blocking rather than re-tested here -- it owns the
    single-finite-buffer rule and the size guard, and a second copy of either is
    how the two drift apart.

    ``why`` is the derivation's own reason when the routing does not apply, so
    the caller learns what about THIS model rules the blocking bound out. It is
    empty when the routing applies, and when the model is not blocked at all.

    Only the UPPER side is routed: the analyzer solves qrf.bas in the 'max'
    direction alone, so 'auto.lower' has no blocking counterpart and keeps
    refusing rather than being answered with the wrong side.
    """
    import numpy as _np

    from ...api.sn import sn_has_blocking
    from ...api.sn.network_struct import SchedStrategy
    from ...api.sn.qrf_blocking import sn_to_qrf_blocking
    # network_struct.SchedStrategy, NOT constants.SchedStrategy: sn.sched holds
    # the former, and the two are distinct enum classes, so comparing across
    # them is always False. The delay test would silently never fire.

    if not sn_has_blocking(sn):
        return '', ''

    if int(sn.nclasses) != 1 or _np.any(_np.isinf(_np.asarray(sn.njobs, dtype=float))):
        return '', ('the QRF blocking bounds are derived for a single-class closed network, '
                    'which this model is not.')
    nservers = _np.asarray(sn.nservers, dtype=float).ravel()
    for i in range(int(sn.nstations)):
        if sn.sched[i] == SchedStrategy.INF:
            return '', ('the QRF blocking bounds model every station as a single server and '
                        'have no infinite-server notion, so a delay station rules them out.')
        if nservers[i] > 1:
            return '', ('the QRF blocking bounds model every station as a single server, so a '
                        'multiserver station rules them out.')

    _blk, blk_msg = sn_to_qrf_blocking(sn)
    if blk_msg:
        return '', blk_msg
    return 'qrf.bas', ''


def ba_resolve_model_method(sn, method):
    """The concrete bound method ``method`` runs as on the model ``sn``.

    Two resolutions in one place, so that runAnalyzer, ba_method_refusal and
    SolverBA.getMethodFeatureSet all judge a name by the same target. First the
    model-free aliases of ba_resolve_method ('default' -> gb.upper, 'auto' ->
    auto.upper, 'qr' -> qrf.mmi, 'lr' -> lr.upper). Then the finite-buffer
    routing of ba_blocking_default: on a model whose buffer BINDS, 'default',
    'auto' and 'auto.upper' run as 'qrf.bas', the one upper bound that models
    the buffer, when the model has the shape that bound needs.

    ``why`` is the reason the routing did NOT apply, and is empty when it
    applied, when the name is not one of the three routable ones, and when the
    model is unblocked. Only the UPPER side routes: the analyzer solves
    'qrf.bas' in the 'max' direction alone, so 'auto.lower' keeps refusing
    rather than being answered with the wrong side. Port of MATLAB
    ba_resolve_model_method.
    """
    from ...api.sn import sn_has_blocking
    why = ''
    routable = method in ('default', 'auto', 'auto.upper')
    method = ba_resolve_method(method)
    if routable and ba_ignores_blocking(method) and sn_has_blocking(sn):
        alt, why = ba_blocking_default(sn)
        if alt:
            method = alt
    return method, why


def solver_ba_analyzer(s, method, options):
    """Return a result dict {QN,UN,RN,TN,AN,XN,WN,CN,runtime,method}."""
    t0 = time.time()
    method = ba_resolve_method(method)

    # The STRUCTURAL premises of every family -- single-class closed, fully
    # closed, single-server -- are asked here and nowhere else, so this run and
    # the report gate SolverBA.supportsModelMethod give one answer whichever a
    # caller meets first. The premises a feature name CAN express (a delay
    # station, an open class) live in SolverBA.getMethodFeatureSet instead.
    refusal = ba_method_refusal(s._sn, method)
    if refusal:
        raise ValueError(refusal)

    if method in BA_QRF_LP or method in BA_QRF_NLP:
        QN, UN, RN, TN, CN, XN = _qrf_lp_bounds(s, method, options)
    elif method.startswith('qrf'):
        raise ValueError(
            "Unknown or unsupported SolverBA QRF method '%s'." % method)
    elif method in BA_LR:
        QN, UN, RN, TN, CN, XN = _lr_bounds(s, method)
    elif method in BA_MAPAMVA:
        QN, UN, RN, TN, CN, XN = _mapamva_bounds(s, method)
    elif method in BA_ASYM:
        QN, UN, RN, TN, CN, XN = _asym_bounds(s, method)
    elif method in BA_CHAIN:
        QN, UN, RN, TN, CN, XN = _chain_bound(s, method)
    elif method in BA_HAREL:
        QN, UN, RN, TN, CN, XN = _harel_bound(s, method)
    elif method in BA_HIER:
        QN, UN, RN, TN, CN, XN = _hier_bound(s, method, options)
    elif method in BA_BPT:
        QN, UN, RN, TN, CN, XN = _bpt_bound(s, method)
    elif method in BA_BGT:
        QN, UN, RN, TN, CN, XN = _bgt_bound(s, method)
    elif method in BA_SNC:
        QN, UN, RN, TN, CN, XN = _snc_bound(s, method)
    elif method in BA_SPNLP:
        QN, UN, RN, TN, CN, XN = _spnlp_bound(s, method, options)
    elif method in BA_AUTO:
        QN, UN, RN, TN, CN, XN = _auto_bound(s, method, options)
    else:
        raise ValueError(
            "Unknown or unsupported SolverBA method '%s'." % method)

    QN = np.asarray(QN, dtype=float); UN = np.asarray(UN, dtype=float)
    RN = np.asarray(RN, dtype=float); TN = np.asarray(TN, dtype=float)
    XN = np.asarray(XN, dtype=float).flatten()
    AN = TN.copy(); WN = RN.copy()
    return {
        'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN, 'XN': XN, 'WN': WN,
        'CN': np.asarray(CN, dtype=float), 'runtime': time.time() - t0,
        'method': method,
    }
