"""
Heavy-usage asymptotic analysis of the closed two-station network with one think
(infinite-server) station and one discriminatory processor-sharing station, by the
generating-function expansion of

  J.A. Morrison, "Asymptotic analysis of a large closed queueing network with discriminatory
  processor sharing", Queueing Systems 9 (1991) 191-214.

Admitted only on the exact shape ``nc_is_dps_model`` tests for. The kernel is
``line_solver.api.npfqn.npfqn_dps_morrison``; this module maps the model struct onto it and lifts
the per-class DPS results into the station-by-class arrays the NC analyzers return.

THERE IS NO NORMALIZING CONSTANT HERE. A DPS station is not product-form -- that is the premise of
the paper -- so lG is NaN, as on the maximum-entropy route. NC hosts this method because NC is
where LINE keeps the asymptotic expansions of generating functions and normalizing-constant
integrals (pana, mmint2, le, ble, gleint, rayint), which is the family Morrison's expansion
belongs to, not because a constant is being computed.

Response times come from Little's law on the queue-length result rather than from the expanded
RESULT 2 (eq. 4.17), so that Q = R*T holds exactly in the returned table; the two agree to the
order of the approximation, since Morrison derives (4.17) as the ratio (4.11)/(4.15).
"""

import time

import numpy as np

from ...api.npfqn import npfqn_dps_morrison
from ...api.io.logging import line_debug, line_warning

__all__ = ['nc_is_dps_model', 'solver_nc_dps_analyzer']


def _sched_name(sn, ist):
    sched = sn.sched.get(ist)
    return getattr(sched, 'name', None)


def _visits(sn):
    """Per-class visit matrix, normalized at the reference station; None when unavailable."""
    try:
        M, K = int(sn.nstations), int(sn.nclasses)
        V = np.zeros((M, K))
        for r in range(K):
            cols = np.nonzero(np.asarray(sn.chains, dtype=float)[:, r])[0]
            if cols.size != 1:
                return None
            vis = sn.visits[int(cols[0])]
            for ist in range(M):
                V[ist, r] = vis[int(sn.stationToStateful[ist]), r]
            vref = V[int(sn.refstat[r]), r]
            if not vref > 0:
                return None
            V[:, r] /= vref
        return V
    except Exception:
        return None


def nc_is_dps_model(sn) -> bool:
    """True when the model is the closed two-station network Morrison's expansion is derived for:
    one infinite-server (think) station and one single-server DPS station, exponential service,
    every class alternating between the two. The shape is checked exactly, not approximately:
    outside it the expansion has no derivation behind it, so a model that misses any clause here is
    left to the ordinary NC routes (which refuse DPS) rather than answered wrongly."""
    if sn is None:
        return False
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    if np.any(np.isinf(njobs)) or not np.any(njobs > 0):
        return False
    if int(sn.nstations) != 2:
        return False
    iInf = [ist for ist in range(2) if _sched_name(sn, ist) == 'INF']
    iDps = [ist for ist in range(2) if _sched_name(sn, ist) == 'DPS']
    if len(iInf) != 1 or len(iDps) != 1:
        return False
    iInf, iDps = iInf[0], iDps[0]
    c = float(np.asarray(sn.nservers, dtype=float).ravel()[iDps])
    if np.isfinite(c) and c != 1:
        return False    # multi-server DPS: the min(n,c) share is not Morrison's
    K = int(sn.nclasses)
    if int(sn.nchains) != K:
        return False    # class switching
    for attr in ('lldscaling', 'cdscaling', 'jdscaling'):
        val = getattr(sn, attr, None)
        if val is not None and np.size(val) > 0:
            return False
    rates = np.asarray(sn.rates, dtype=float)
    scv = np.asarray(sn.scv, dtype=float)
    schedparam = np.asarray(sn.schedparam, dtype=float)
    for r in range(K):
        if njobs[r] <= 0:
            return False
        for ist in (iInf, iDps):
            if not np.isfinite(rates[ist, r]) or rates[ist, r] <= 0:
                return False
            if np.isfinite(scv[ist, r]) and abs(scv[ist, r] - 1) > 1e-6:
                return False
        if not np.isfinite(schedparam[iDps, r]) or schedparam[iDps, r] <= 0:
            return False
    V = _visits(sn)
    if V is None:
        return False
    for r in range(K):
        if abs(V[iInf, r] - V[iDps, r]) > 1e-9 * max(1.0, V[iInf, r]):
            return False    # unequal visits: not the alternating cycle of the paper
    return True


def solver_nc_dps_analyzer(sn, options):
    """Analyze the closed think+DPS network.

    Returns
    -------
    (QN, UN, RN, TN, CN, XN, lG, runtime, it, method) with lG = NaN.
    """
    tstart = time.time()
    # The shape is re-checked here, not assumed from the caller: this analyzer is
    # reachable from SolverNC, the dispatch chain and directly from user code, and every
    # quantity below -- think time, DPS service time, weights, visit ratios -- is
    # meaningless off the shape the expansion was derived for.
    if not nc_is_dps_model(sn):
        raise ValueError(
            "solver_nc_dps_analyzer applies only to a CLOSED network of exactly two stations, one "
            "infinite-server (think) station and one single-server DPS station with exponential "
            "service and one visit each per cycle (see nc_is_dps_model).")
    M, K = int(sn.nstations), int(sn.nclasses)
    iInf = [ist for ist in range(M) if _sched_name(sn, ist) == 'INF'][0]
    iDps = [ist for ist in range(M) if _sched_name(sn, ist) == 'DPS'][0]

    rates = np.asarray(sn.rates, dtype=float)
    Npop = np.asarray(sn.njobs, dtype=float).ravel()[:K]
    Z = 1.0 / rates[iInf, :K]
    S = 1.0 / rates[iDps, :K]
    w = np.asarray(sn.schedparam, dtype=float)[iDps, :K]

    line_debug("NC DPS analyzer: N=%s, Z=%s, S=%s, w=%s", Npop, Z, S, w, options=options)
    res = npfqn_dps_morrison(Npop, Z, S, w)

    Qdps = np.asarray(res.Q, dtype=float)
    if np.any(Qdps < 0) or np.any(Qdps > Npop) or np.any(~np.isfinite(Qdps)):
        line_warning(
            "solver_nc_dps_analyzer",
            "The Morrison expansion returned queue lengths outside [0,N] (usage rho=%g, a=%g). "
            "The model is outside the moderately-heavy regime the expansion assumes; treat the "
            "result as unreliable and cross-check with SolverMVA or SolverCTMC." % (res.rho, res.a))
        Qdps = np.minimum(np.maximum(Qdps, 0.0), Npop)
    X = (Npop - Qdps) / Z

    V = _visits(sn)
    QN = np.zeros((M, K))
    UN = np.zeros((M, K))
    RN = np.zeros((M, K))
    TN = np.zeros((M, K))
    QN[iDps, :] = Qdps
    QN[iInf, :] = Npop - Qdps            # population conservation (exact, closed)
    for ist in range(M):
        TN[ist, :] = X * V[ist, :]
    UN[iInf, :] = QN[iInf, :]            # INF utilization convention
    c = float(np.asarray(sn.nservers, dtype=float).ravel()[iDps])
    if not np.isfinite(c) or c <= 0:
        c = 1.0
    UN[iDps, :] = X * (V[iDps, :] * S) / c
    with np.errstate(divide='ignore', invalid='ignore'):
        RN = np.where(TN > 0, QN / np.where(TN > 0, TN, 1.0), 0.0)
    CN = np.where(X > 0, Npop / np.where(X > 0, X, 1.0), 0.0)

    return QN, UN, RN, TN, CN, X, float('nan'), time.time() - tstart, 1, 'morrison'
