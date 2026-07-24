"""
Level-Dependent QBD (LD-QBD) solver for single-class Delay/Queue networks.

Two regimes are handled:
  CLOSED: one Delay (INF) + one Queue, finite population N. Level n = jobs at
          the queue (0 <= n <= N); arrival rate from the delay is (N-n)*lambda.
  OPEN:   one Source (EXT) + one Queue, open class (Poisson arrivals). Level
          n = jobs at the queue, truncated at Nlev; arrival rate is the
          constant external rate lambda. Nlev is taken from options.cutoff or
          chosen so the truncated tail probability is negligible.

Both regimes share the same block-tridiagonal generator, differing only in the
per-level arrival rate and the top level. Service is min(n,c)*mu (exact M/M/c
boundary) or its PH generalisation.

Exactness: exact for exponential service at any number of servers, and for PH
service at a single server. For PH service with c > 1 servers it is an
approximation: the c parallel PH servers are collapsed into one PH process
scaled by min(n,c), which ignores the phase of each individual busy server.

References:
    Original MATLAB: matlab/src/solvers/MAM/solver_mam_ldqbd.m
"""

import time
from typing import Optional, Tuple

import numpy as np

from . import MAMAlgorithm, MAMResult
from ....api.mam import ldqbd, LdqbdOptions
from ....api.mam.map_analysis import map_mean, map_pie
from ....api.qsys.retrial import _proc_to_d0d1
from ....api.sn.predicates import sn_has_load_dependence
from ....lang.base import SchedStrategy


def _count_sched(sn, target):
    """Indices of the stations scheduled with the given strategy.

    sn.sched is a dict keyed by station index in Python native; compare by
    enum name rather than raw integer, per the cross-codebase enum convention.
    """
    sched = sn.sched
    out = []
    for i in range(sn.nstations):
        s = sched.get(i, None) if isinstance(sched, dict) else sched[i]
        if s is None:
            continue
        sname = s.name if hasattr(s, 'name') else str(s)
        if sname == target.name:
            out.append(i)
    return out


def _proc_of(sn, ist):
    """(D0, D1) of the class-0 process at station ist.

    Native sn.proc stores processes in compact dict form ({'rate'},
    {'probs','rates'}, {'k','mu'}); MAP/PH are an explicit [D0, D1] pair.
    _proc_to_d0d1 normalizes every form to the Markovian [D0, D1] this
    construction needs.
    """
    proc = sn.proc
    # sn.proc[station] is indexed by class; K == 1 in every LDQBD regime.
    proc_st = proc.get(ist, None) if isinstance(proc, dict) else proc[ist]
    if proc_st is None:
        raise ValueError('LDQBD method found no process at station %d.' % ist)
    entry = proc_st[0]
    d0d1 = _proc_to_d0d1(entry)
    if d0d1 is None:
        raise ValueError('LDQBD method could not parse the process at station %d.' % ist)
    return np.atleast_2d(d0d1[0]), np.atleast_2d(d0d1[1])


def ldqbd_is_closed_delay_queue(sn) -> bool:
    """Single-class closed Delay+Queue: the regime routed to LDQBD by default.

    Mirrors the isClosedDelayQueue subfunction of the MATLAB
    solver_mam_analyzer and isClosedDelayQueue in the JAR Solver_mam_analyzer.
    The open Source+Queue regime is deliberately excluded: it truncates the
    level space at options.cutoff, so it is not unconditionally preferable to
    dec.source and stays opt-in.
    """
    if sn.nclasses != 1 or sn.nstations != 2:
        return False
    if not np.all(np.isfinite(np.asarray(sn.njobs, dtype=float))):
        return False
    return (len(_count_sched(sn, SchedStrategy.INF)) == 1
            and len(_count_sched(sn, SchedStrategy.FCFS)) == 1)


class LDQBDAlgorithm(MAMAlgorithm):
    """Level-dependent QBD for single-class Delay/Queue networks."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        if sn.nclasses != 1:
            return False, 'The ldqbd method requires a single-class model.'
        if sn.nstations != 2:
            return False, 'The ldqbd method requires exactly two stations.'

        is_open = not np.all(np.isfinite(np.asarray(sn.njobs, dtype=float)))
        n_queue = len(_count_sched(sn, SchedStrategy.FCFS))
        if is_open:
            n_source = len(_count_sched(sn, SchedStrategy.EXT))
            if n_source != 1 or n_queue != 1:
                return False, ('Open LDQBD method requires exactly one Source '
                               'and one Queue station.')
        else:
            n_delay = len(_count_sched(sn, SchedStrategy.INF))
            if n_delay != 1 or n_queue != 1:
                return False, ('Closed LDQBD method requires exactly one Delay '
                               'and one Queue station.')
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        start_time = time.time()

        ok, reason = LDQBDAlgorithm.supports_network(sn)
        if not ok:
            raise ValueError(reason)

        M = sn.nstations
        K = sn.nclasses
        rates = np.asarray(sn.rates, dtype=float)
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        is_open = not np.all(np.isfinite(njobs))

        queue_idx = _count_sched(sn, SchedStrategy.FCFS)[0]
        if is_open:
            src_idx = _count_sched(sn, SchedStrategy.EXT)[0]
            delay_idx = None
        else:
            delay_idx = _count_sched(sn, SchedStrategy.INF)[0]
            src_idx = None
            N = int(njobs[0])

        # ---- Service process at the queue ----
        D0, D1 = _proc_of(sn, queue_idx)
        n_servers = int(np.asarray(sn.nservers, dtype=float).ravel()[queue_idx])

        is_ph = D0.shape[0] > 1
        if is_ph:
            n_phases = D0.shape[0]
            alpha = np.atleast_2d(map_pie(D0, D1)).reshape(1, -1)
            mean_service = map_mean(D0, D1)
        else:
            n_phases = 1
            mu = -D0[0, 0]
            alpha = None
            mean_service = 1.0 / mu

        # ---- Per-level service factor: LD scaling if set, else min(n,c) ----
        lld = None
        lldlimit = 0
        has_lld = False
        lldscaling = getattr(sn, 'lldscaling', None)
        if sn_has_load_dependence(sn) and lldscaling is not None:
            lldarr = np.asarray(lldscaling, dtype=float)
            if lldarr.size > 0 and lldarr.shape[0] > queue_idx and np.any(lldarr[queue_idx, :] != 1):
                has_lld = True
                lld = lldarr[queue_idx, :]
                lldlimit = lld.size
        sf_max = lld[lldlimit - 1] if has_lld else float(n_servers)

        # ---- Arrival rate per level and number of levels ----
        # K == 1 here, so the station-class index of station i is just i.
        rt = np.asarray(sn.rt, dtype=float)
        if is_open:
            src_D0, _ = _proc_of(sn, src_idx)
            if src_D0.shape[0] > 1:
                raise ValueError('Open LDQBD method currently supports Poisson '
                                 '(exponential) arrivals only; the Source uses '
                                 'a MAP/MMPP process.')
            lam = rates[src_idx, 0]
            lambda_eff = lam * rt[src_idx, queue_idx]
            rho = lambda_eff * mean_service / sf_max
            if rho >= 1:
                raise ValueError('Open LDQBD method requires a stable queue '
                                 '(rho = %.4f >= 1). Increase service capacity '
                                 'or reduce the arrival rate.' % rho)
            cutoff = getattr(options, 'cutoff', None) if options is not None else None
            if cutoff is not None and np.isscalar(cutoff) and np.isfinite(cutoff):
                nlev = max(n_servers + 1, int(round(cutoff)))
            else:
                tail_tol = 1e-10
                nlev = n_servers + int(np.ceil(np.log(tail_tol) / np.log(rho)))
                nlev = min(max(nlev, n_servers + 10), 100000)
            arr_rate = lambda_eff * np.ones(nlev + 1)
            arr_rate[nlev] = 0.0  # truncation: no arrivals above the top level
        else:
            lambda_d = rates[delay_idx, 0]
            lambda_eff = lambda_d * rt[delay_idx, queue_idx]
            nlev = N
            # finite-source rate (N-n)*lambda_eff, 0 at n=N
            arr_rate = (N - np.arange(N + 1)) * lambda_eff

        # Per-level service factor sf(n), n = 1..nlev (sf[n-1] holds level n)
        sf = np.zeros(nlev)
        for n in range(1, nlev + 1):
            if has_lld:
                sf[n - 1] = lld[min(n, lldlimit) - 1]
            else:
                sf[n - 1] = min(n, n_servers)

        # ---- Construct LD-QBD block-tridiagonal generator ----
        Q0 = []  # upward (arrival), levels 0..nlev-1
        Q1 = []  # local, levels 0..nlev
        Q2 = []  # downward (departure), levels 1..nlev

        if not is_ph:
            for n in range(nlev):
                Q0.append(np.array([[arr_rate[n]]]))
            for n in range(nlev + 1):
                departure_rate = sf[n - 1] * mu if n > 0 else 0.0
                Q1.append(np.array([[-(arr_rate[n] + departure_rate)]]))
            for n in range(1, nlev + 1):
                Q2.append(np.array([[sf[n - 1] * mu]]))
        else:
            eye_p = np.eye(n_phases)
            # level 0 -> 1: start service in a phase
            Q0.append(arr_rate[0] * alpha)
            for n in range(1, nlev):
                # level n -> n+1: preserve phase
                Q0.append(arr_rate[n] * eye_p)
            # level 0: only arrivals
            Q1.append(np.array([[-arr_rate[0]]]))
            for n in range(1, nlev + 1):
                Q1.append(sf[n - 1] * D0 - arr_rate[n] * eye_p)
            # level 1 -> 0: empty the queue
            Q2.append(sf[0] * D1 @ np.ones((n_phases, 1)))
            for n in range(2, nlev + 1):
                # level n -> n-1: complete and restart
                Q2.append(sf[n - 1] * D1)

        # ---- Solve LD-QBD ----
        tol = getattr(options, 'tol', 1e-10) if options is not None else 1e-10
        max_iter = getattr(options, 'iter_max', None) if options is not None else None
        if max_iter is None and options is not None:
            max_iter = getattr(options, 'max_iter', 1000)
        ldqbd_options = LdqbdOptions(epsilon=tol, max_iter=max_iter or 1000)
        result = ldqbd(Q0, Q1, Q2, ldqbd_options)
        pi_ldqbd = np.asarray(result.pi, dtype=float).ravel()

        # ---- Performance metrics ----
        mean_queue = float(np.arange(nlev + 1) @ pi_ldqbd)

        # see _kb/06-solver-catalog.md (MAM: "LDQBD") -- both forms below are
        # already per-server; do not rescale by n_servers again
        if has_lld or n_servers == 1:
            util_ps = 1.0 - pi_ldqbd[0]
        else:
            util_ps = 0.0
            for n in range(1, nlev + 1):
                util_ps += (min(n, n_servers) / n_servers) * pi_ldqbd[n]

        QN = np.zeros((M, K))
        UN = np.zeros((M, K))
        RN = np.zeros((M, K))
        TN = np.zeros((M, K))
        CN = np.zeros((1, K))
        XN = np.zeros((1, K))

        if is_open:
            # Accepted/served throughput = arrival rate minus truncation blocking
            X = lambda_eff * (1.0 - pi_ldqbd[nlev])
            r_queue = mean_queue / X if X > 0 else 0.0
            # Source station: pass-through, no queueing.
            TN[src_idx, 0] = X
            # Queue station.
            QN[queue_idx, 0] = mean_queue
            UN[queue_idx, 0] = util_ps
            RN[queue_idx, 0] = r_queue
            TN[queue_idx, 0] = X
            XN[0, 0] = X
            CN[0, 0] = r_queue
        else:
            mean_delay = N - mean_queue
            X = mean_delay * lambda_eff
            r_queue = mean_queue / X if X > 0 else 0.0
            r_delay = 1.0 / rates[delay_idx, 0]
            # see _kb/06-solver-catalog.md (MAM: "LDQBD") -- delay TN != queue
            # flow X whenever the delay routes elsewhere (rt < 1)
            QN[delay_idx, 0] = mean_delay
            UN[delay_idx, 0] = mean_delay  # infinite server: U = Q
            RN[delay_idx, 0] = r_delay
            TN[delay_idx, 0] = mean_delay * lambda_d
            # Queue station metrics
            QN[queue_idx, 0] = mean_queue
            UN[queue_idx, 0] = util_ps
            RN[queue_idx, 0] = r_queue
            TN[queue_idx, 0] = X
            XN[0, 0] = X
            CN[0, 0] = r_delay + r_queue

        return MAMResult(
            QN=QN, UN=UN, RN=RN, TN=TN, CN=CN, XN=XN,
            totiter=1,  # LDQBD is a direct method
            method='ldqbd',
            runtime=time.time() - start_time,
        )
