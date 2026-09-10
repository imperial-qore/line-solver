"""
The single-station fluid limits.

Native Python twin of matlab/src/solvers/FLD/solver_fluid_qsys_analyzer.m: a
Source -> Queue -> Sink model with one class, answered by a closed-form fluid or
Gaussian limit rather than by integrating the network drift.

WHY THESE ARE FLUID METHODS AND NOT MVA ONES. Each depends on the service or
patience law BEYOND ITS MEAN -- the stationary point of the Liu-Whitt model is
where the patience ccdf crosses 1/rho, the Mt/G/inf mean is a convolution with
the service ccdf -- and each is the limit of a sequence of systems, not an
approximation to a fixed one. That is the fluid solver's contract.
"""

from typing import Any, Dict, Tuple

import numpy as np

from ...sn.network_struct import NodeType


METHODS = ('ggisgi.fluid', 'fluid.ggisgi', 'ggingi.tga', 'fluid.tga',
           'tvms', 'fluid.tvms', 'mtginf', 'fluid.mtginf', 'mol', 'fluid.mol')


def _station_of_type(sn, ty):
    for i in range(len(sn.nodetype)):
        if sn.nodetype[i] == ty:
            return int(np.asarray(sn.nodeToStation).ravel()[i])
    return None


def fluid_qsys_horizon(options) -> Tuple[bool, str, float, float]:
    """The integration window of a time-varying single-station fluid limit.

    The time-varying limits ('mol', 'mtginf', 'tvms') report a TRAJECTORY, so an
    infinite or absent upper end of options.timespan leaves them nothing to
    report.

    A horizon is a solver OPTION and not a model feature, so the feature
    registry has no name for it and SolverFLD.supportsModelMethod has to ask
    this predicate directly. ``_horizon`` below asks the same one on the solve
    path, which is what keeps the report and the run from disagreeing about
    whether a method can be asked for.

    Args:
        options: solver options carrying the timespan.

    Returns:
        (ok, reason, t0, t1); reason is '' and t1 is meaningful when ok.
    """
    t0, t1 = 0.0, 1.0
    ts = getattr(options, 'timespan', None)
    if ts is not None and len(ts) >= 2:
        if np.isfinite(ts[0]):
            t0 = float(ts[0])
        t1 = float(ts[1])
    if not np.isfinite(t1) or t1 <= t0:
        return (False, 'A time-varying fluid method needs a finite horizon: set '
                       'options.timespan = [t0, t1].', t0, t1)
    return True, '', t0, t1


def _horizon(options) -> Tuple[float, float]:
    """The integration window, refused by name when it is not a finite interval."""
    ok, reason, t0, t1 = fluid_qsys_horizon(options)
    if not ok:
        raise RuntimeError(reason)
    return t0, t1


def solver_fluid_qsys_analyzer(sn, options) -> Dict[str, Any]:
    """
    Solve a single-station model with one of the closed-form fluid limits.

    Methods: ``ggisgi.fluid`` (Liu and Whitt, Operations Research 60(5), 2012),
    ``ggingi.tga`` (Liu, Whitt and Yu, Naval Research Logistics 63(3), 2016),
    ``tvms`` at constant staffing (Liu and Whitt, INFORMS J. Computing 26(1),
    2014), ``mtginf`` (Eick, Massey and Whitt, Management Science 39(2), 1993)
    and ``mol`` (Massey and Whitt, Annals of Applied Probability 4(4), 1994).

    Returns:
        Dict with ``QN``, ``UN``, ``RN``, ``TN``, ``CN``, ``XN``, the transient
        tables ``Qt``, ``Ut``, ``Tt`` (empty for the two stationary methods) and
        ``method``.
    """
    from ...sn import sn_patience_handles, sn_arrival_rate_fun
    from ...qsys import (qsys_ggisgi_fluid, qsys_ggingi_tga, qsys_gtmtst_fluid,
                          qsys_mtginf, qsys_mtgs0_mol)
    from ...mam.map_analysis import map_cdf

    method = str(getattr(options, 'method', 'default')).lower()
    M, K = int(sn.nstations), int(sn.nclasses)
    QN = np.zeros((M, K)); UN = np.zeros((M, K)); RN = np.zeros((M, K)); TN = np.zeros((M, K))
    AN = np.zeros((M, K)); CN = np.zeros(K); XN = np.zeros(K)

    src = _station_of_type(sn, NodeType.SOURCE)
    qi = _station_of_type(sn, NodeType.QUEUE)
    if qi is None:
        qi = _station_of_type(sn, NodeType.DELAY)
    # THE SHAPE THESE LIMITS ARE STATED FOR, refused by name rather than answered
    # on a model they do not describe: one open class through one queueing
    # station. The MVA qsys analyzer is reached by a structural dispatch that
    # guarantees it; these methods are selected by NAME, so the check lives here.
    if src is None or qi is None:
        raise RuntimeError('the single-station fluid limits need a Source and a queueing station')
    if K != 1 or int(getattr(sn, 'nclosedjobs', 0) or 0) > 0:
        raise RuntimeError("the '%s' method is a single-station limit: it needs one open class "
                           "through one Source and one queueing station" % method)

    Vq = float(np.asarray(sn.visits[0])[int(np.asarray(sn.stationToStateful).ravel()[qi]), 0])
    lam = float(np.asarray(sn.rates)[src, 0]) * Vq
    mu = float(np.asarray(sn.rates)[qi, 0])
    nserv = float(np.asarray(sn.nservers).ravel()[qi])
    scv_s = float(np.asarray(sn.scv)[qi, 0])
    ca = float(np.sqrt(np.asarray(sn.scv)[src, 0]))
    cs = float(np.sqrt(scv_s))
    h = sn_patience_handles(sn, qi, 0)

    # The service ccdf, needed by the two Mt/G methods: they are exact in the
    # service DISTRIBUTION, not in its mean, which is the whole point of the
    # Eick-Massey-Whitt lag.
    pair = None
    try:
        slot = sn.proc[qi][0]
        if slot is not None and len(slot) >= 2 and np.shape(slot[0]) == np.shape(slot[1]):
            pair = (np.asarray(slot[0], dtype=float), np.asarray(slot[1], dtype=float))
    except Exception:
        pair = None
    if pair is None:
        def service_ccdf(x, _mu=mu):
            return float(np.exp(-_mu * float(x)))
    else:
        def service_ccdf(x, _p=pair):
            return float(1.0 - map_cdf(_p[0], _p[1], np.atleast_1d(float(x)))[0])
    ES = 1.0 / mu
    ES2 = (1.0 + scv_s) * ES * ES

    Qt = []; Ut = []; Tt = []

    def stationary(Lsys, Tq, Uq):
        # Little's law on the CARRIED rate, as every LINE solver reports a
        # station that loses work.
        RN[qi, 0] = (Lsys / Tq) if Tq > 0 else 0.0
        QN[qi, 0] = Lsys
        UN[qi, 0] = Uq
        TN[qi, 0] = Tq
        TN[src, 0] = lam / Vq
        # The OFFERED rate, so that getAvgLossTable reads the abandonment as
        # ArvR - Tput. Left to default to TN it would report no loss at all.
        AN[qi, 0] = lam
        XN[0] = Tq
        CN[0] = RN[qi, 0] * Vq

    def transient(t, Lt, Ut_, Tt_, arrival):
        # The steady-state row of a time-varying model is the TIME AVERAGE over
        # the horizon, which is what a stationary reader of a periodic system
        # measures; the trajectory itself is returned beside it.
        t = np.asarray(t, dtype=float).ravel()
        Lt = np.asarray(Lt, dtype=float).ravel()
        Ut_ = np.asarray(Ut_, dtype=float).ravel()
        Tt_ = np.asarray(Tt_, dtype=float).ravel()
        arrival = np.asarray(arrival, dtype=float).ravel()
        span = float(t[-1] - t[0])
        if span <= 0:
            Lbar, Ubar, Tbar, Abar = Lt[0], Ut_[0], Tt_[0], arrival[0]
        else:
            Lbar = float(np.trapezoid(Lt, t) / span)
            Ubar = float(np.trapezoid(Ut_, t) / span)
            Tbar = float(np.trapezoid(Tt_, t) / span)
            Abar = float(np.trapezoid(arrival, t) / span)
        QN[qi, 0] = Lbar
        UN[qi, 0] = Ubar
        TN[qi, 0] = Tbar
        TN[src, 0] = Abar / Vq
        RN[qi, 0] = (Lbar / Tbar) if Tbar > 0 else 0.0
        AN[qi, 0] = Abar
        XN[0] = Tbar
        CN[0] = RN[qi, 0] * Vq
        zero = np.column_stack([np.zeros_like(t), t])
        qt = [[zero.copy() for _ in range(K)] for _ in range(M)]
        ut = [[zero.copy() for _ in range(K)] for _ in range(M)]
        tt = [[zero.copy() for _ in range(K)] for _ in range(M)]
        qt[qi][0] = np.column_stack([Lt, t])
        ut[qi][0] = np.column_stack([Ut_, t])
        tt[qi][0] = np.column_stack([Tt_, t])
        tt[src][0] = np.column_stack([arrival / Vq, t])
        return qt, ut, tt

    if method in ('ggisgi.fluid', 'fluid.ggisgi'):
        if h is None:
            raise RuntimeError("the 'ggisgi.fluid' method needs a reneging patience law on the "
                               "queue (Queue.setPatience)")
        res = qsys_ggisgi_fluid(lam, mu, int(nserv), h['ccdf'])
        stationary(res['meanNumber'], res['throughput'], res['utilization'])
        actual = 'ggisgi.fluid'
    elif method in ('ggingi.tga', 'fluid.tga'):
        if h is None:
            raise RuntimeError("the 'ggingi.tga' method needs a reneging patience law on the "
                               "queue (Queue.setPatience)")
        if not np.isfinite(nserv) or nserv < 1:
            raise RuntimeError("the 'ggingi.tga' method needs a finite number of servers")
        res = qsys_ggingi_tga(lam, mu, int(nserv), ca, cs, h['ccdf'],
                              patiencePdf=h['pdf'], serviceCcdf=service_ccdf)
        stationary(res['meanNumber'], lam * (1.0 - res['probAbandon']),
                   min(res['meanNumberInService'] / nserv, 1.0))
        actual = 'ggingi.tga'
    elif method in ('tvms', 'fluid.tvms'):
        if h is None:
            raise RuntimeError("the 'tvms' method needs a reneging patience law on the queue "
                               "(Queue.setPatience)")
        if not np.isfinite(nserv) or nserv < 1:
            raise RuntimeError("the 'tvms' method needs a finite number of servers")
        lam_fun, _, _ = sn_arrival_rate_fun(sn, src, 0)
        t0, t1 = _horizon(options)
        # CONSTANT STAFFING. Nothing in a Network declares a time-varying server
        # count, so s(t) is the station's own s; the time variation the method
        # is for enters through lambda(t) alone. A staffing schedule would need
        # a model feature that does not exist, and inventing one here would make
        # the solver answer a model the user did not build.
        res = qsys_gtmtst_fluid(lam_fun,
                                lambda t, _s=nserv: (_s if np.ndim(t) == 0
                                                     else np.full(np.shape(t), _s)),
                                lambda t, _m=mu: (_m if np.ndim(t) == 0
                                                  else np.full(np.shape(t), _m)),
                                h['ccdf'], t1 - t0, patiencePdf=h['pdf'])
        Qt, Ut, Tt = transient(t0 + np.asarray(res['times']), res['X'], res['utilization'],
                               mu * np.asarray(res['B']), res['arrivalRate'])
        actual = 'tvms'
    elif method in ('mtginf', 'fluid.mtginf'):
        lam_fun, _, _ = sn_arrival_rate_fun(sn, src, 0)
        t0, t1 = _horizon(options)
        res = qsys_mtginf(lam_fun, service_ccdf, ES, np.linspace(t0, t1, 200), ES2=ES2)
        # An infinite-server station serves everything that arrives, so the
        # throughput is the arrival rate and the busy-server count is what a
        # utilization column can carry.
        Qt, Ut, Tt = transient(res['times'], res['meanNumber'], res['meanNumber'],
                               res['arrivalRate'], res['arrivalRate'])
        actual = 'mtginf'
    elif method in ('mol', 'fluid.mol'):
        if not np.isfinite(nserv) or nserv < 1:
            raise RuntimeError("the 'mol' method needs a finite number of servers")
        lam_fun, _, _ = sn_arrival_rate_fun(sn, src, 0)
        t0, t1 = _horizon(options)
        cap = float(np.asarray(sn.cap).ravel()[qi])
        # A finite buffer beyond the servers is not part of the loss model the
        # approximation is for; only s servers and no waiting room is.
        use_delay = bool(np.isfinite(cap) and cap > nserv)
        res = qsys_mtgs0_mol(lam_fun, service_ccdf, ES, int(nserv), np.linspace(t0, t1, 200),
                             delay=use_delay, ES2=ES2)
        busy = np.asarray(res['meanBusyMOL'], dtype=float)
        Qt, Ut, Tt = transient(res['times'], busy, busy / nserv, mu * busy, res['arrivalRate'])
        actual = 'mol'
    else:
        raise RuntimeError("the '%s' method is not a single-station fluid limit" % method)

    return {'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'AN': AN, 'CN': CN, 'XN': XN,
            'Qt': Qt, 'Ut': Ut, 'Tt': Tt, 'method': actual}
