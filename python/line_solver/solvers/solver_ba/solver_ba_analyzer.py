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
)
from ...api.pfqn.bound_hierarchies import (
    pfqn_pbh, pfqn_cbh, pfqn_pbk, pfqn_bjbk, pfqn_ssd, pfqn_mcub,
    pfqn_sib, pfqn_ldbcmp,
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
    if s.nclasses != 1 or not (s.njobs[0] > 0 and np.isfinite(s.njobs[0])):
        raise ValueError(
            "Method '%s' supports single-class closed networks only." % method)
    rates = np.asarray(s.rates[:, 0], dtype=float)
    dem = np.asarray(s.demands[:, 0], dtype=float)   # = V/rate per station
    with np.errstate(divide='ignore', invalid='ignore'):
        V = dem * rates
    inf = _inf_mask(s)
    D = dem[~inf]
    Z = float(np.sum(dem[inf]))
    N = int(round(float(s.njobs[0])))
    return V, rates, inf, D, Z, N


def _reject_multiserver(s, inf, method):
    if np.any(np.asarray(s.nservers)[~inf] > 1):
        raise ValueError(
            "Method '%s' does not support multi-server stations (use 'ssd')." % method)


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
        U = (T[:, 0] / rates).reshape(-1, 1)
    U[inf, 0] = Q[inf, 0]
    return Q, U, R, T, np.array([[float(C)]]), np.array([X])


def _lr_bounds(s, method):
    """LP-based Linear Reduction bounds (Casale et al.).

    Solves the linear-reduction LP relaxation once per station, minimizing
    ('lr.lower') or maximizing ('lr.upper') that station's utilization. Both
    senses are valid bounds because the LP relaxation contains the exact
    solution. Unlike 'qrf.mmi.linear' -- whose name refers only to its explicit
    Aeq/beq constraint representation while its objective is the nonlinear MEM
    entropy -- this method is a pure LP end to end (HiGHS via scipy.linprog).
    """
    from ...api.mapqn.parameters import PFParameters, LinearReductionParameters
    from ...api.mapqn.bnd_lr_pf import mapqn_bnd_lr_pf
    from ...api.mapqn.bnd_lr import mapqn_bnd_lr

    V, rates, inf, D, Z, N = _sc_setup(s, method)
    _reject_multiserver(s, inf, method)
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


def _asym_bounds(s, method):
    """ABA/BJB/PB/SB/GB noniterative bounds (explicit formulas)."""
    from ...lang.base import SchedStrategy  # noqa: F401
    V, rates, inf, D, Z, N = _sc_setup(s, method)
    _reject_multiserver(s, inf, method)
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
            Cv = Z + A1 + (N - 1) * (A1 * A2 + A3) / (A1 ** 2 + A2)
            X = min(1.0 / Dmax, N / Cv)
            return fill(X, False, Cv)
        AN = np.sum(D ** N)
        Cv = Z + A1 + (N - 1) * (AN / A1) ** (1.0 / (N - 1))
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
                Rv[i] = 1.0 / rates[i]; Qv[i] = X * Rv[i]
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
    if np.any(np.asarray(s.nservers)[~inf] > 1) and method.startswith('mwba'):
        raise ValueError("Method '%s' does not support multi-server stations." % method)
    dc = sn_get_demands_chain(sn)
    Lchain, STchain, Vchain, alpha, Nchain = (
        dc.Lchain, dc.STchain, dc.Vchain, dc.alpha, dc.Nchain)
    M = s.nstations
    nchains = Lchain.shape[1]
    Zc = np.sum(Lchain[inf, :], axis=0)
    Lq = Lchain[~inf, :]

    if method in ('cub.upper', 'mbjb.lower'):
        if not (s.njobs[0] > 0) or np.any(~np.isfinite(np.asarray(s.njobs))):
            raise ValueError("Method '%s' supports fully closed networks only." % method)
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
    Qchain = np.zeros((M, nchains))
    for c in range(nchains):
        Tchain[:, c] = Xchain[c] * Vchain[:, c]
        Uchain[:, c] = Xchain[c] * Lchain[:, c]
        Qchain[:, c] = Uchain[:, c]
    deagg = sn_deaggregate_chain_results(
        sn, Lchain, None, STchain, Vchain, alpha, Qchain, Uchain, None,
        Tchain, None, Xchain)
    Cn = np.zeros((M, s.nclasses))
    return deagg.Q, deagg.U, deagg.R, deagg.T, Cn, np.asarray(deagg.X).flatten()


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
    if fam in ('pbh', 'cbh', 'pbk', 'bjbk', 'sib'):
        _reject_multiserver(s, inf, method)
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


BA_ASYM = {'aba.upper', 'aba.lower', 'bjb.upper', 'bjb.lower',
           'pb.upper', 'pb.lower', 'sb.upper', 'sb.lower',
           'gb.upper', 'gb.lower'}
BA_CHAIN = {'mwba.upper', 'mwba.lower', 'cub.upper', 'mbjb.lower'}
BA_LR = {'lr', 'lr.upper', 'lr.lower'}

# The two LP-backed QRF reduction bounds. These are served by the stable
# qr_bounds_bas / qr_bounds_rsrd simplex backends (via scipy HiGHS), NOT by the
# unstable qrf_noblo_* NLP library, so unlike the other seven qrf tokens they
# are advertised and dispatched natively.
BA_QRF_LP = {'qrf.bas', 'qrf.rsrd'}

# The NLP-backed QRF tokens, served by api.mapqn.qrf_noblo_*. Withheld until
# 2026-07-20 because the SolverBA -> solver_ctmc_qrf_analyzer -> qrf_noblo_*
# path had never been exercised end to end. Exercising it found and fixed a
# transposed v in the adapter and a throughput inversion that assumed the
# reference station was a delay; both codebases now derive UN from the QRF
# utilization directly. On a 2-station closed chain the K == 1 result is exact
# against the CTMC and the K == 2 result lands inside the glpsol bound range.
BA_QRF_NLP = {'qr', 'qrf.mmi', 'qrf.mem', 'qrf.mmi.ld', 'qrf.mmi.linear'}


def _qrf_lp_bounds(s, method, options):
    """QRF reduction bounds, both the LP-backed and the NLP-backed tokens.

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
BA_HIER = {'pbh.upper', 'pbh.lower', 'cbh.upper', 'cbh.lower',
           'pbk.upper', 'pbk.lower', 'bjbk.upper', 'bjbk.lower',
           'ssd.upper', 'ssd.lower', 'sib.upper', 'sib.lower',
           'ldbcmp.lower'}


def solver_ba_analyzer(s, method, options):
    """Return a result dict {QN,UN,RN,TN,AN,XN,WN,CN,runtime,method}."""
    t0 = time.time()
    if method == 'default':
        method = 'gb.upper'

    if method in BA_QRF_LP or method in BA_QRF_NLP:
        QN, UN, RN, TN, CN, XN = _qrf_lp_bounds(s, method, options)
    elif method in ('qrf.bas.mmi', 'qrf.bas.mem'):
        # The BAS-blocking NLP variants (qrf_bas_mmi / qrf_bas_mem) have no
        # native backend: api.mapqn ships the no-blocking qrf_noblo_* family
        # only. Unlike the no-blocking tokens these were never merely
        # withheld pending coverage, so they stay unimplemented.
        raise NotImplementedError(
            "QRF BAS-blocking method '%s' is not implemented natively: "
            "api.mapqn provides the no-blocking qrf_noblo_* family only. "
            "Use the JAR or MATLAB backend for the BAS variants." % method)
    elif method.startswith('qrf'):
        raise ValueError(
            "Unknown or unsupported SolverBA QRF method '%s'." % method)
    elif method in BA_LR:
        QN, UN, RN, TN, CN, XN = _lr_bounds(s, method)
    elif method in BA_ASYM:
        QN, UN, RN, TN, CN, XN = _asym_bounds(s, method)
    elif method in BA_CHAIN:
        QN, UN, RN, TN, CN, XN = _chain_bound(s, method)
    elif method in BA_HIER:
        QN, UN, RN, TN, CN, XN = _hier_bound(s, method, options)
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
