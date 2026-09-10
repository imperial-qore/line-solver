"""
SolverMVA method 'amva.mapqn': the horizontal-cut mean value analysis of a closed
multiclass network with one exponential delay station and one FCFS
single-server queue whose class-r service is a MAP (api.mapqn.mapqn_amva).

One structural predicate, three callers: listValidMethods drops the name when
`mva_mapqn_reason` is nonempty, supportsModelMethod reports it, and the
analyzer raises it, so the offered, reported and run answers cannot drift.
Mirrors matlab/src/solvers/MVA/mva_mapqn_reason.m and
solver_mva_mapqn_analyzer.m.
"""
import time

import numpy as np

from ...sn.proc_form import proc_to_map
from ....constants import ProcessType
from .analyzers import MVAResult


def _name(v):
    return str(getattr(v, 'name', v)).upper()


def _is_delay(sn, ist):
    sched = getattr(sn, 'sched', None)
    st = sched.get(ist) if isinstance(sched, dict) else (sched[ist] if sched is not None and ist < len(sched) else None)
    nservers = np.asarray(getattr(sn, 'nservers', []), dtype=float).ravel()
    return (st is not None and _name(st) == 'INF') or (ist < nservers.size and np.isinf(nservers[ist]))


_MARKOVIAN = (ProcessType.EXP, ProcessType.ERLANG, ProcessType.HYPEREXP, ProcessType.PH,
              ProcessType.APH, ProcessType.COXIAN, ProcessType.COX2, ProcessType.MAP, ProcessType.MMPP2)


def mva_mapqn_reason(sn) -> str:
    """The reason 'amva.mapqn' cannot solve sn, or '' when it can: exactly one
    exponential infinite-server station and one FCFS single-server queue with
    Markovian (MAP-representable) service, closed classes cycling
    delay -> queue -> delay without class switching."""
    if sn is None:
        return "Method 'amva.mapqn' needs a network model."
    M = int(sn.nstations)
    R = int(sn.nclasses)
    if M != 2:
        return "Method 'amva.mapqn' requires exactly two stations: one delay (infinite server) and one FCFS queue."
    delays = [i for i in range(M) if _is_delay(sn, i)]
    if len(delays) != 1:
        return "Method 'amva.mapqn' requires exactly one delay (infinite-server) station and one queue."
    id_ = delays[0]
    iq = 1 - id_
    sched = sn.sched
    stq = sched.get(iq) if isinstance(sched, dict) else sched[iq]
    if _name(stq) != 'FCFS':
        return "Method 'amva.mapqn' requires FCFS scheduling at the queue; station %d is %s." % (iq + 1, _name(stq))
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    if nservers[iq] != 1:
        return "Method 'amva.mapqn' supports a single-server queue only."
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isinf(njobs)) or int(getattr(sn, 'nclosedjobs', 0)) <= 0:
        return "Method 'amva.mapqn' supports closed models only."
    procid = np.asarray(sn.procid, dtype=object)
    rt = np.asarray(sn.rt, dtype=float)
    s2f = np.asarray(sn.stationToStateful, dtype=int).ravel() if getattr(sn, 'stationToStateful', None) is not None else np.arange(M)
    sd, sq = int(s2f[id_]), int(s2f[iq])
    tol = 1e-12
    for r in range(R):
        if njobs[r] <= 0:
            continue
        if procid[id_, r] != ProcessType.EXP:
            return ("Method 'amva.mapqn' requires exponential think times; class %d has a %s think time."
                    % (r + 1, _name(procid[id_, r])))
        if procid[iq, r] not in _MARKOVIAN:
            return ("Method 'amva.mapqn' requires a Markovian (MAP-representable) service process at the queue; "
                    "class %d is not." % (r + 1))
        D0, D1 = proc_to_map(sn.proc[iq][r])
        if D0 is None:
            return ("Method 'amva.mapqn' requires a Markovian (MAP-representable) service process at the queue; "
                    "class %d is not." % (r + 1))
        if abs(rt[sd * R + r, sq * R + r] - 1.0) > tol or abs(rt[sq * R + r, sd * R + r] - 1.0) > tol:
            return ("Method 'amva.mapqn' requires every class to cycle delay -> queue -> delay without class "
                    "switching; class %d does not." % (r + 1))
    return ''


def mva_supports_mapqn(sn, method):
    """(ok, reason) for 'amva.mapqn'. Mirrors MATLAB SolverMVA.supportsMapqn."""
    from .handler import mva_base_method
    if mva_base_method(method) != 'mapqn':
        return True, ''
    reason = mva_mapqn_reason(sn)
    return reason == '', reason


def solver_mva_mapqn_analyzer(sn, options=None) -> MVAResult:
    """Run mapqn_amva on sn and assemble the station-class metrics: at the
    queue Q_r, U_r = X_r E[S_r], R_r = Q_r / X_r; at the delay
    Q_r = U_r = X_r Z_r, R_r = Z_r; C_r = N_r / X_r."""
    from ...mapqn.amva import mapqn_amva
    t0 = time.time()
    reason = mva_mapqn_reason(sn)
    if reason:
        raise ValueError(reason)
    M = int(sn.nstations)
    R = int(sn.nclasses)
    id_ = [i for i in range(M) if _is_delay(sn, i)][0]
    iq = 1 - id_
    N = np.asarray(np.round(np.asarray(sn.njobs, dtype=float).ravel()), dtype=int)
    rates = np.asarray(sn.rates, dtype=float)
    mu = np.ones(R)
    D0s, D1s = [], []
    for r in range(R):
        if N[r] > 0:
            mu[r] = float(rates[id_, r])
            D0, D1 = proc_to_map(sn.proc[iq][r])
            D0s.append(np.atleast_2d(np.asarray(D0, dtype=float)))
            D1s.append(np.atleast_2d(np.asarray(D1, dtype=float)))
        else:
            D0s.append(np.array([[-1.0]]))
            D1s.append(np.array([[1.0]]))
    res = mapqn_amva(mu, D0s, D1s, N)
    QN = np.zeros((M, R)); UN = np.zeros((M, R)); RN = np.zeros((M, R)); TN = np.zeros((M, R))
    CN = np.zeros(R); XN = np.zeros(R)
    for r in range(R):
        if N[r] <= 0 or res.X[r] <= 0:
            continue
        X = float(res.X[r])
        XN[r] = X
        QN[iq, r] = res.Qq[r]; UN[iq, r] = res.U[r]; TN[iq, r] = X; RN[iq, r] = res.Qq[r] / X
        QN[id_, r] = X / mu[r]; UN[id_, r] = X / mu[r]; TN[id_, r] = X; RN[id_, r] = 1.0 / mu[r]
        CN[r] = N[r] / X
    out = MVAResult()
    out.QN, out.UN, out.RN, out.TN, out.CN, out.XN = QN, UN, RN, TN, CN, XN
    out.AN = TN.copy(); out.WN = RN.copy()
    out.logNormConstAggr = np.nan
    out.iter = int(np.prod(N + 1))
    out.runtime = time.time() - t0
    out.method = 'amva.mapqn'
    return out
