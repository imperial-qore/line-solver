"""Per-station load-dependent rate scaling alpha(i,n) for the QRF bounds.

Twin of ``matlab/src/api/sn/sn_to_qrf_alpha.m``,
``jar/src/main/java/jline/api/sn/SnToQrfAlpha.java`` and
``cpp/include/line/api/sn/sn_to_qrf_alpha.h``.

Read by the load-dependent QRF arms ``qrf.mmi.ld`` and ``qrf.mmi.linear``,
which are the only two the scaling reaches: every other arm builds a
population-free q and has nowhere to put a rate that depends on n. The
companion tables of the BLOCKING formulation live in ``qrf_blocking.py``; the
two are independent derivations off the same ``NetworkStruct``.
"""

import numpy as np

__all__ = ['sn_to_qrf_alpha']


def sn_to_qrf_alpha(sn):
    """Return ``(alpha, msg, ld, peak)``: the QRF load-dependent rate scaling.

    The load-dependent QRF arms carry a scaling ``alpha[i, n]`` that multiplies
    EVERY rate out of station i while it holds n jobs, completions mu and
    background phase changes v alike -- see the q construction in
    ``qrf_noblo_mmi_ld``. That is exactly the rate law of

      an infinite server         alpha(i,n) = n
      a c-server station         alpha(i,n) = min(n, c_i)
      limited load dependence    alpha(i,n) = sn.lldscaling[i, n]

    so the three COMPOSE BY MULTIPLICATION and not one of them is an
    approximation: the relaxed chain is the model's own, and the QRF answer
    keeps whatever status it had on a single-server model.

    WHERE IT STOPS BEING THE MODEL'S OWN IS PHASE-TYPE SERVICE AT A STATION
    THAT SERVES SEVERAL JOBS AT ONCE. The QRF local state carries ONE phase per
    station, a faithful description of one job in service and of nothing else:
    min(n,c) jobs served in parallel each advance through a phase of their own,
    and no scaling of a single-phase process reproduces that joint motion. A
    multiserver or delay station must therefore be exponential -- scaling a PH
    server by min(n,c) would answer a DIFFERENT chain, so the relaxation would
    stop containing the model's stationary distribution and the number would
    bound nothing. Limited load dependence at a SINGLE server is exempt and
    admits PH freely: one job is in service whatever the rate.

    ``alpha`` is tabulated at n = 1..N. ``sn.lldscaling`` may be narrower than N
    and is read clamped at its last column, the convention ``pfqn_lldfun`` uses.
    ``ld`` says the scaling is not identically 1, and stays TRUE through a
    refusal so a caller can tell "no arm serves this" from "the arm you asked
    for does not".

    THE UTILIZATION NORMALIZER IS THE DECLARED PEAK, NOT max(alpha). LINE
    reports U = T*S/peak at every station whose rate scales with the
    population, one convention shared by multiserver, lld and class dependence
    (see ``cd_peak_scaling``, "the same convention as solver_ncld does for
    lldscaling, U/max(lldscaling)"). ``peak`` is therefore nservers(i) times the
    largest lld scaling the model can REACH, and not max(alpha[i]): at c = 3
    with N = 2 the reachable alpha peaks at 2 while the station still has three
    servers, and normalizing by 2 would report a utilization the model never
    attains. Inf at a delay, where LINE reports U = QN instead.
    """
    from .network_struct import SchedStrategy

    M = int(sn.nstations)
    N = float(np.sum(np.asarray(sn.njobs, dtype=float)))
    peak = np.ones(M)
    if not np.isfinite(N) or N < 1:
        return (np.ones((M, 1)),
                'the QRF bounds need a closed model with a finite population.', False, peak)
    N = int(round(N))

    alpha = np.ones((M, N))
    lld = getattr(sn, 'lldscaling', None)
    lld = np.asarray(lld, dtype=float) if lld is not None and np.size(lld) else None
    smax = lld.shape[1] if lld is not None and lld.ndim == 2 else 0
    sched = getattr(sn, 'sched', None)
    nservers = np.asarray(sn.nservers, dtype=float).ravel()

    for i in range(M):
        c = float(nservers[i])
        is_delay = np.isinf(c) or (
            sched is not None and int(sched[i]) == int(SchedStrategy.INF))
        serves_many = bool(is_delay or c > 1)
        ki = _qrf_alpha_phases(sn, i)
        if serves_many and ki > 1:
            # ld stays True through the refusal: the model IS load dependent.
            msg = ('station %d serves %s jobs at once with %d-phase service, and the QRF '
                   'local state carries one phase per station, which describes one job in '
                   'service and no more. Give that station exponential service, or use a '
                   'single-server model.'
                   % (i + 1, 'unboundedly many' if is_delay else 'up to %d' % int(c), ki))
            return alpha, msg, True, peak
        lldpeak = 1.0
        for n in range(1, N + 1):
            if is_delay:
                alpha[i, n - 1] = n
            elif c > 1:
                alpha[i, n - 1] = min(n, c)
            if smax > 0:
                lldpeak = max(lldpeak, lld[i, min(n, smax) - 1])
                alpha[i, n - 1] *= lld[i, min(n, smax) - 1]
        peak[i] = np.inf if is_delay else c * lldpeak

    return alpha, '', bool(np.any(alpha != 1.0)), peak


def _qrf_alpha_phases(sn, i):
    """Phases of station i's service process, read from ``sn.proc`` as the QRF
    adapter reads them: that {D0, D1} pair is what sizes the local state, so
    testing it keeps the refusal and the formulation on one quantity."""
    proc = getattr(sn, 'proc', None)
    try:
        entry = proc[i][0]
        if entry is not None and len(entry) and entry[0] is not None:
            return int(np.asarray(entry[0]).shape[0])
    except (TypeError, KeyError, IndexError):
        pass
    return 1

