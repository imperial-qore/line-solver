"""
Matrix-analytic network analysis for closed networks (mna_closed).

The QNA two-moment traffic equations carry the flow rates and SCVs around the
network exactly as in `mna_open`; what closes the model is an OUTER BISECTION
on the per-class throughput, bracketed below by 0 and above by the slowest
service rate the class meets at a finite-server station, and driven against the
closed population target sum_i Q(i,k) = N(k). Each FCFS station is solved as an
MMAP[K]/PH[K]/1 queue whose queue-length DISTRIBUTION is truncated at the
population and renormalized, so the closed marginal comes from the
matrix-analytic solution rather than from an open-queue formula.

Port of matlab/src/solvers/MAM/solver_mna_closed.m, with the C++ twin
(cpp/include/line/solvers/mam/solver_mna.h, solver_mna_closed) as the second
reading. Three quirks of the reference are reproduced deliberately, not
repaired, so the four codebases answer alike:

1. `[pdistr] = MMAPPH1FCFS(..., 'ncDistr', maxLevel)` captures only the FIRST
   output of the butools routine, i.e. CLASS 1's marginal, and the per-class
   truncation loop then reads every class off that one vector. The C++ port
   mirrors this (`pd[0]`); the JAR reads ncDistr per class and so parts company
   on multiclass models.
2. The INF branch writes `d2(ist,s) = a2(ist,s)` into what the outer sweep
   declared as an (M,1) vector, so the value the splitting step later reads
   back as `d2(ist)` is `a2(ist,1)`, the FIRST class's flow SCV. Reproduced by
   carrying d2 as (M,K) and reading column 0.
3. X is left at zero: the reference allocates it and never assigns it, so the
   per-class system-throughput column of a mna_closed result is zero in all
   four codebases. The throughput a caller wants is TN.

References:
    MATLAB: matlab/src/solvers/MAM/solver_mna_closed.m
    QNA: Whitt, W. "The Queueing Network Analyzer." Bell Labs Tech. J. (1983)
"""

import numpy as np
import time
from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult
from ..utils.network_adapter import check_closed_network, extract_visit_counts


class MNAClosedAlgorithm(MAMAlgorithm):
    """Matrix-analytic network analysis for closed networks."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by mna_closed.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        is_closed = check_closed_network(sn)
        if not is_closed:
            return False, "mna_closed is for closed networks. Use mna_open for open networks."

        # see _kb/06-solver-catalog.md (MAM: "mna rejects self-looping chains").
        # A class is self-looping when it VISITS one station only, which is what
        # MATLAB SolverMAM.supportsModelMethod tests; the reference station being
        # a queue is not by itself a confinement, and reading refstat instead
        # refused multiclass models that MATLAB solves.
        from ....constants import GlobalConstants
        try:
            sched = sn.sched
            V = extract_visit_counts(sn)
            njobs = np.asarray(sn.njobs, dtype=float).ravel()
            for k in range(sn.nclasses):
                # njobs is Inf for an open class (see mna self-looping note above)
                if not np.isfinite(njobs[k]):
                    continue
                vis = np.flatnonzero(V[:, k] > GlobalConstants.FineTol)
                if vis.size != 1:
                    continue
                st = int(vis[0])
                s = sched.get(st, None) if isinstance(sched, dict) else sched[st]
                sname = s.name if hasattr(s, 'name') else str(s)
                if sname not in ('INF', 'EXT'):
                    return False, (
                        "The mna method does not support self-looping classes "
                        "(class %d is confined to station %d with no "
                        "inter-station flow to decompose). Use the dec.source method."
                        % (k + 1, st + 1)
                    )
        except Exception:
            pass

        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve a closed network by matrix-analytic network analysis."""
        return solver_mna_closed(sn, options)


def _sched_name(sched):
    """Scheduling-strategy name. LINE keeps distinct SchedStrategy enum copies
    (lang.base, constants, api.sn.network_struct) whose members compare unequal
    across copies, so the discipline must be tested by name."""
    return getattr(sched, 'name', None)


def _node_name(nodetype):
    """Node-type name; see _sched_name for why the comparison is by name."""
    return getattr(nodetype, 'name', None)


def solver_mna_closed(sn, options=None) -> MAMResult:
    """QNA traffic equations closed by a bisection on the per-class throughput."""
    from ....api.da import da_fpi
    from ....api.mam import map_pie, mmap_super
    from ....api.sn import sn_rt_stations
    from ....api.sn.proc_form import proc_to_map
    from ....constants import GlobalConstants
    from ....distributions.markovian import APH
    from ....lib.thirdparty.butools.queues import MMAPPH1FCFS
    import numpy.matlib as ml

    start_time = time.time()

    K = int(sn.nclasses)
    M = int(sn.nstations)
    C = int(sn.nchains)
    I = int(sn.nnodes)

    tol = getattr(options, 'tol', GlobalConstants.CoarseTol) if options is not None \
        else GlobalConstants.CoarseTol
    max_iter = getattr(options, 'max_iter', 100) if options is not None else 100
    space_max = 16  # the reference overrides config.space_max for this analyzer

    # sn.rt and sn.visits are indexed by stateful node: project them onto stations
    rt, V = sn_rt_stations(sn)
    rates = np.asarray(sn.rates, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        S = 1.0 / rates
    scv = np.asarray(sn.scv, dtype=float).copy()
    scv[np.isnan(scv)] = 0
    N = np.asarray(sn.njobs, dtype=float).ravel()
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    sched = sn.sched if sn.sched else {}
    inchain = sn.inchain if sn.inchain else {}
    isslc = np.asarray(sn.isslc, dtype=bool).ravel() if getattr(sn, 'isslc', None) is not None \
        else np.zeros(K, dtype=bool)
    refstat = np.asarray(sn.refstat).ravel().astype(int)

    Q = np.zeros((M, K))
    U = np.zeros((M, K))
    R = np.zeros((M, K))
    T = np.zeros((M, K))
    X = np.zeros((1, K))

    # service processes, as the PH pair (pie, D0) MMAPPH1FCFS consumes
    pie = {}
    D0 = {}
    for ist in range(M):
        if _sched_name(sched.get(ist)) not in ('FCFS', 'INF', 'PS'):
            continue
        for k in range(K):
            entry = sn.proc[ist][k] if sn.proc is not None else None
            D0k, D1k = proc_to_map(entry) if entry is not None else (None, None)
            if D0k is None or np.any(np.isnan(np.asarray(D0k, dtype=float))):
                # a class disabled at the station passes through instantaneously
                D0[(ist, k)] = np.array([[-GlobalConstants.Immediate]])
                pie[(ist, k)] = np.array([1.0])
            else:
                D0k = np.atleast_2d(np.asarray(D0k, dtype=float))
                D1k = np.atleast_2d(np.asarray(D1k, dtype=float))
                pie[(ist, k)] = np.asarray(map_pie(D0k, D1k), dtype=float).ravel()
                D0[(ist, k)] = D0k

    QNc = N.copy()
    # the bracket: a class cannot be served faster than the slowest rate it
    # meets at a station that has finitely many servers
    lambda_lb = np.zeros(K)
    lambda_ub = np.zeros(K)
    finite = np.where(np.isfinite(nservers))[0]
    for k in range(K):
        if finite.size:
            col = rates[finite, k]
            col = col[np.isfinite(col)]
            lambda_ub[k] = np.min(col) if col.size else 0.0

    lam = np.zeros(max(K, C))
    QN = np.zeros(K)
    a1 = np.zeros((M, K))
    a2 = np.zeros((M, K))
    d2 = np.zeros((M, K))
    f2 = np.zeros((M * K, M * K))

    maxLevel = int(np.sum(N[np.isfinite(N)])) + 1

    def _renormalize_chains():
        """Scale each chain's queue lengths onto its population.

        The reference indexes the chain loop's counter into the CLASS axis of Q,
        which is the same column whenever each chain holds one class; that is
        the indexing reproduced here.
        """
        for c in range(C):
            if c >= K:
                continue
            colsum = np.sum(Q[:, c])
            with np.errstate(divide='ignore', invalid='ignore'):
                Q[:, c] = N[c] * Q[:, c] / colsum

    def mna_flow_sweep(x, itnum):
        """The inner QNA fixed point on the flow rates a1 and their SCVs a2."""
        xref = np.concatenate([a1.ravel(), a2.ravel()])

        _renormalize_chains()

        for k in range(K):
            if isslc[k]:
                Q[:, k] = 0
                Q[refstat[k], k] = N[k]

        # throughputs seeded from the current throughput iterate and the visits
        if itnum == 1:
            for c in range(C):
                classes_c = np.asarray(inchain[c], dtype=int).ravel()
                for m in range(M):
                    T[m, classes_c] = V[m, classes_c] * lam[c]

        # superposition
        for ist in range(M):
            a1[ist, :] = 0
            a2[ist, :] = 0
            lambda_i = np.sum(T[ist, :])
            for jst in range(M):
                for r in range(K):
                    for s in range(K):
                        p = rt[jst * K + s, ist * K + r]
                        if p == 0:
                            continue
                        a1[ist, r] += T[jst, s] * p
                        if lambda_i > 0:
                            a2[ist, r] += (1.0 / lambda_i) * f2[jst * K + s, ist * K + r] * T[jst, s] * p

        # per-station update
        for ind in range(I):
            if not sn.isstation[ind]:
                continue
            if _node_name(sn.nodetype[ind]) == 'JOIN':
                continue
            ist = int(sn.nodeToStation[ind])
            schedI = _sched_name(sched.get(ist))
            if schedI == 'INF':
                # see the module docstring: the reference stores the whole row
                # here and reads column 0 back in the splitting step
                d2[ist, :] = a2[ist, :]
                for c in range(C):
                    for k in np.asarray(inchain[c], dtype=int).ravel():
                        T[ist, k] = a1[ist, k]
                        U[ist, k] = S[ist, k] * T[ist, k]
                        Q[ist, k] = T[ist, k] * S[ist, k] * V[ist, k]
                        R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0
            elif schedI == 'PS':
                for c in range(C):
                    classes_c = np.asarray(inchain[c], dtype=int).ravel()
                    for k in classes_c:
                        T[ist, k] = lam[c] * V[ist, k]
                        U[ist, k] = S[ist, k] * T[ist, k]
                    Uden = min(1 - GlobalConstants.FineTol, np.sum(U[ist, :]))
                    Nc = np.sum(N[classes_c])
                    for k in classes_c:
                        # geometric-bound approximation truncated at the chain population
                        Q[ist, k] = (U[ist, k] - U[ist, k] ** (Nc + 1)) / (1 - Uden)
                        R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0
            elif schedI == 'FCFS':
                mu_ist = rates[ist, :K].copy()
                mu_ist[np.isnan(mu_ist)] = 0
                rho_ist_class = a1[ist, :K] / (GlobalConstants.FineTol + rates[ist, :K])
                rho_ist_class[np.isnan(rho_ist_class)] = 0
                lambda_ist = np.sum(a1[ist, :])
                mi = int(nservers[ist])
                rho_ist = np.sum(rho_ist_class) / mi
                if rho_ist < 1 - tol:
                    mubar = lambda_ist / rho_ist
                    c2 = -1.0
                    for r in range(K):
                        if mu_ist[r] > 0:
                            c2 += a1[ist, r] / lambda_ist * (mubar / mi / mu_ist[r]) ** 2 * (scv[ist, r] + 1)
                    d2[ist, 0] = 1 + rho_ist ** 2 * (c2 - 1) / np.sqrt(mi) \
                        + (1 - rho_ist ** 2) * (np.sum(a2[ist, :]) - 1)
                else:
                    for k in range(K):
                        Q[ist, k] = N[k]
                    d2[ist, 0] = 1
                for k in range(K):
                    T[ist, k] = a1[ist, k]
                    U[ist, k] = T[ist, k] * S[ist, k] / nservers[ist]
                    R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0
            elif _node_name(sn.nodetype[ind]) == 'FORK':
                raise RuntimeError('Fork nodes not supported yet by QNA solver.')

        # splitting: update the flow SCVs
        for ist in range(M):
            for jst in range(M):
                if _node_name(sn.nodetype[int(sn.stationToNode[jst])]) == 'SOURCE':
                    continue
                for r in range(K):
                    for s in range(K):
                        p = rt[ist * K + r, jst * K + s]
                        if p > 0:
                            f2[ist * K + r, jst * K + s] = 1 + p * (d2[ist, 0] - 1)

        return np.concatenate([a1.ravel(), a2.ravel()]), xref

    def mna_outer_sweep(x, itout):
        """One bisection step on the per-class throughput."""
        nonlocal QN
        if itout != 1:
            for k in range(K):
                if QN[k] < QNc[k]:
                    lambda_lb[k] = lam[k]
                else:
                    lambda_ub[k] = lam[k]
                lam[k] = (lambda_ub[k] + lambda_lb[k]) / 2
        else:
            lam[:K] = lambda_ub

        Q[:, :] = 0
        U[:, :] = 0
        R[:, :] = 0
        T[:, :] = 0
        X[:, :] = 0
        a1[:, :] = 0
        a2[:, :] = 0
        d2[:, :] = 0

        f2[:, :] = 0
        for ist in range(M):
            for jst in range(M):
                if _node_name(sn.nodetype[int(sn.stationToNode[jst])]) == 'SOURCE':
                    continue
                for r in range(K):
                    for s in range(K):
                        if rt[ist * K + r, jst * K + s] > 0:
                            f2[ist * K + r, jst * K + s] = 1

        da_fpi(mna_flow_sweep, np.concatenate([a1.ravel(), a2.ravel()]),
               max_iter + 1, tol, nanstop=True)

        # queue lengths from the matrix-analytic marginal of each FCFS station,
        # truncated at the closed population and renormalized
        for ind in range(I):
            if not sn.isstation[ind]:
                continue
            ist = int(sn.nodeToStation[ind])
            if _sched_name(sched.get(ist)) != 'FCFS':
                continue
            rho_ist_class = a1[ist, :K] / (GlobalConstants.FineTol + rates[ist, :K])
            rho_ist_class[np.isnan(rho_ist_class)] = 0
            mi = int(nservers[ist])
            rho_ist = np.sum(rho_ist_class) / mi
            if rho_ist >= 1 - tol:
                for k in range(K):
                    Q[ist, k] = N[k]
                    R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0
                continue

            arri_node = None
            for k in range(K):
                if a1[ist, k] == 0:
                    # an idle class contributes no arrivals: a zero-rate stream
                    arri_class = [np.array([[0.0]]), np.array([[0.0]]), np.array([[0.0]])]
                else:
                    aph = APH.fit_mean_and_scv(1.0 / a1[ist, k], a2[ist, k])
                    D0k = np.atleast_2d(aph.getD0())
                    D1k = np.atleast_2d(aph.getD1())
                    arri_class = [D0k, D1k, D1k]
                if k == 0:
                    arri_node = arri_class
                else:
                    arri_node = mmap_super(arri_node, arri_class, 'default')

            # MMAPPH1FCFS consumes [D0, Dclass1..DclassK]
            D_arr = [arri_node[0]] + [arri_node[2 + k] for k in range(K)]
            pie_list = [pie[(ist, k)] for k in range(K)]
            D0_list = [D0[(ist, k)] for k in range(K)]

            Dlam = np.sum([np.sum(Dk) for Dk in D_arr[1:]])
            if Dlam < GlobalConstants.FineTol:
                for k in range(K):
                    Q[ist, k] = GlobalConstants.FineTol / rates[ist, 0]
            else:
                pdistr = MMAPPH1FCFS(
                    [ml.matrix(np.asarray(D, dtype=np.float64)) for D in D_arr],
                    [ml.matrix(np.asarray(pv, dtype=np.float64).reshape(1, -1)) for pv in pie_list],
                    [ml.matrix(np.asarray(Tm, dtype=np.float64)) for Tm in D0_list],
                    'ncDistr', maxLevel)
                # only the FIRST output is read, i.e. class 1's marginal; see
                # quirk 1 in the module docstring
                if isinstance(pdistr, (list, tuple)):
                    pdistr = pdistr[0]
                pdistr = np.asarray(pdistr, dtype=float).ravel()
                for k in range(K):
                    Nk = int(N[k]) if np.isfinite(N[k]) else maxLevel - 1
                    pk = np.abs(pdistr[:Nk + 1]).copy()
                    pk[-1] = abs(1 - np.sum(pdistr[:-1]))
                    denom = np.sum(pk)
                    if denom > 0:
                        pk = pk / denom
                    Q[ist, k] = max(0.0, min(float(N[k]), float(np.arange(Nk + 1) @ pk)))

            for k in range(K):
                R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0

        QN = np.sum(Q, axis=0)
        return QN.copy(), QNc.copy()

    _, it_out, _ = da_fpi(mna_outer_sweep, np.zeros(K), max_iter, tol, nanstop=True)

    for k in range(K):
        if isslc[k]:
            Q[:, k] = 0
            ist = int(refstat[k])
            Q[ist, k] = N[k]
            T[ist, k] = N[k] * rates[ist, k]
            R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0
            U[ist, k] = S[ist, k] * T[ist, k]

    _renormalize_chains()

    for ist in range(M):
        if _sched_name(sched.get(ist)) == 'INF':
            # an infinite server's utilization IS its queue length
            U[ist, :] = Q[ist, :]

    CN = np.sum(R, axis=0, keepdims=True)
    Q = np.abs(Q)
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    CN[np.isnan(CN)] = 0
    X[np.isnan(X)] = 0

    runtime = time.time() - start_time

    return MAMResult(
        QN=Q,
        UN=U,
        RN=R,
        TN=T,
        CN=CN,
        XN=X,
        totiter=int(it_out),
        method="mna_closed",
        runtime=runtime
    )
