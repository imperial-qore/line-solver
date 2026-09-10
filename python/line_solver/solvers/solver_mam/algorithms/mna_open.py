"""
Matrix-analytic network analysis for open networks (mna_open).

The QNA two-moment traffic equations carry the flow rates and SCVs around the
network; each FCFS station is then solved as an MMAP[K]/PH[K]/1 queue, so the
queue lengths come from the matrix-analytic marginal rather than from a
diffusion formula.

Port of matlab/src/solvers/MAM/solver_mna_open.m (and of the JAR twin
jline.solvers.mam.handlers.Solver_mna_open, which resolves the PS branch that
MATLAB leaves writing into dead variables).

References:
    MATLAB: matlab/src/solvers/MAM/solver_mna_open.m
    QNA: Whitt, W. "The Queueing Network Analyzer." Bell Labs Tech. J. (1983)
"""

import numpy as np
import time
from typing import Tuple, Optional

from . import MAMAlgorithm, MAMResult


class MNAOpenAlgorithm(MAMAlgorithm):
    """Matrix-analytic network analysis for open networks."""

    @staticmethod
    def supports_network(sn) -> Tuple[bool, Optional[str]]:
        """Check if network can be solved by mna_open.

        Args:
            sn: NetworkStruct

        Returns:
            (can_solve, reason_if_not)
        """
        is_open = any(np.isinf(n) for n in sn.njobs) if len(sn.njobs) > 0 else True
        if not is_open:
            return False, "mna_open is for open networks. Use mna_closed for closed networks."
        return True, None

    def solve(self, sn, options=None) -> MAMResult:
        """Solve an open network by matrix-analytic network analysis.

        Args:
            sn: NetworkStruct
            options: MNA options (tol, max_iter, verbose, dep_scv)

        Returns:
            MAMResult
        """
        return solver_mna_open(sn, options)


def _sched_name(sched):
    """Scheduling-strategy name. LINE keeps distinct SchedStrategy enum copies
    (lang.base, constants, api.sn.network_struct) whose members compare unequal
    across copies, so the discipline must be tested by name."""
    return getattr(sched, 'name', None)


def _node_name(nodetype):
    """Node-type name; see _sched_name for why the comparison is by name."""
    return getattr(nodetype, 'name', None)


def solver_mna_open(sn, options=None) -> MAMResult:
    """QNA traffic equations closed by a per-station MMAP[K]/PH[K]/1 solve."""
    from ....api.da import da_fpi, da_traffic_superpos
    from ....api.mam import map_pie, mmap_super
    from ....api.npfqn import npfqn_traffic_split_rr
    from ....api.qsys import qsys_mmck
    from ....api.sn import sn_rt_stations
    from ....api.sn.proc_form import proc_to_map
    from ....api.solvers.mam.mmap_fj import mam_detect_mmck, mam_truncate_renorm
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

    # sn.rt and sn.visits are indexed by stateful node: project them onto stations
    rt, V = sn_rt_stations(sn)
    rates = np.asarray(sn.rates, dtype=float)
    with np.errstate(divide='ignore', invalid='ignore'):
        S = 1.0 / rates
    scv = np.asarray(sn.scv, dtype=float).copy()
    scv[np.isnan(scv)] = 0
    njobs = np.asarray(sn.njobs, dtype=float).ravel()
    nservers = np.asarray(sn.nservers, dtype=float).ravel()
    sched = sn.sched if sn.sched else {}

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

    a1 = np.zeros((M, K))
    a2 = np.zeros((M, K))
    d2 = np.zeros(M)
    f2 = np.zeros((M * K, M * K))
    # deterministic (round-robin) split degrees, k=1 where the split is Markovian
    kRR = npfqn_traffic_split_rr(sn)
    for ist in range(M):
        for jst in range(M):
            if _node_name(sn.nodetype[int(sn.stationToNode[jst])]) == 'SOURCE':
                continue
            for r in range(K):
                for s in range(K):
                    p = rt[ist * K + r, jst * K + s]
                    if p > 0:
                        f2[ist * K + r, jst * K + s] = 1 + p * (1 - kRR[ist, r])

    lambda_chain = np.zeros(C)
    d2c = np.zeros(C)
    source_idx = 0
    inchain = sn.inchain if sn.inchain else {}
    for c in range(C):
        classes_c = np.asarray(inchain[c], dtype=int).ravel()
        source_idx = int(np.asarray(sn.refstat).ravel()[classes_c[0]])
        lam_c = rates[source_idx, classes_c]
        scv_c = scv[source_idx, classes_c]
        lambda_chain[c] = np.sum(lam_c[np.isfinite(lam_c)])
        d2c[c] = da_traffic_superpos(lam_c, scv_c)
        T[source_idx, classes_c] = lam_c
    # the source is station 0 in every open model LINE builds
    d2[source_idx] = float(d2c @ lambda_chain) / np.sum(lambda_chain)

    def mna_sweep(x, itnum):
        xref = np.concatenate([a1.ravel(), a2.ravel()])

        if itnum == 1:
            for c in range(C):
                classes_c = np.asarray(inchain[c], dtype=int).ravel()
                for m in range(M):
                    T[m, classes_c] = V[m, classes_c] * lambda_chain[c]

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
            ist = int(sn.nodeToStation[ind])
            schedI = _sched_name(sched.get(ist))
            if schedI == 'INF':
                d2[ist] = np.sum(a2[ist, :])
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
                        T[ist, k] = lambda_chain[c] * V[ist, k]
                        U[ist, k] = S[ist, k] * T[ist, k]
                    Uden = min(1 - GlobalConstants.FineTol, np.sum(U[ist, :]))
                    for k in classes_c:
                        Q[ist, k] = U[ist, k] / (1 - Uden)
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
                    d2[ist] = 1 + rho_ist ** 2 * (c2 - 1) / np.sqrt(mi) \
                        + (1 - rho_ist ** 2) * (np.sum(a2[ist, :]) - 1)
                else:
                    for k in range(K):
                        Q[ist, k] = njobs[k]
                    d2[ist] = 1
                for k in range(K):
                    T[ist, k] = a1[ist, k]
                    U[ist, k] = T[ist, k] * S[ist, k] / nservers[ist]

        # splitting: update the flow SCVs
        for ist in range(M):
            for jst in range(M):
                if _node_name(sn.nodetype[int(sn.stationToNode[jst])]) == 'SOURCE':
                    continue
                for r in range(K):
                    for s in range(K):
                        p = rt[ist * K + r, jst * K + s]
                        if p > 0:
                            # k-fold convolution then Bernoulli thinning at q=k*p: C^2 = (q/k)*d2+1-q
                            f2[ist * K + r, jst * K + s] = 1 + p * (d2[ist] - kRR[ist, r])

        return np.concatenate([a1.ravel(), a2.ravel()]), xref

    _, totiter, _ = da_fpi(mna_sweep, np.concatenate([a1.ravel(), a2.ravel()]),
                           max_iter + 1, tol, nanstop=True)

    # queue lengths from the matrix-analytic marginal of each FCFS station
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
                Q[ist, k] = njobs[k]
                R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else np.inf
            continue

        arri_node = None
        for k in range(K):
            if a1[ist, k] == 0:
                # an idle class contributes no arrivals: a zero-rate Poisson stream
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

        if np.isfinite(sn.cap[ist]):
            capK = int(sn.cap[ist])
            is_mmck, mu_mmck = mam_detect_mmck(sn, ist, K, arri_node)
            if is_mmck:
                aggr_lambda = np.nansum(a1[ist, :K])
                exact = qsys_mmck(aggr_lambda, mu_mmck, int(nservers[ist]), capK)
                meanQ_fc = exact['meanQueueLength']
                loss_fc = exact['lossProbability']
            else:
                meanQ_fc, loss_fc, _ = mam_truncate_renorm(D_arr, pie_list, D0_list, capK)
            lambda_inflow = np.nan_to_num(a1[ist, :K])
            T_eff = lambda_inflow * (1 - loss_fc)
            sumT = np.sum(T_eff)
            if sumT > 0:
                Savg_eff = np.nansum(T_eff * S[ist, :K]) / sumT
                Wq = max(0.0, meanQ_fc / sumT - Savg_eff)
            else:
                Wq = 0.0
            for k in range(K):
                T[ist, k] = T_eff[k]
                U[ist, k] = T[ist, k] * S[ist, k] / nservers[ist]
                if T[ist, k] > 0:
                    R[ist, k] = Wq + S[ist, k]
                    Q[ist, k] = T[ist, k] * R[ist, k]
                else:
                    R[ist, k] = 0.0
                    Q[ist, k] = 0.0
        else:
            Qret = MMAPPH1FCFS(
                [ml.matrix(np.asarray(D, dtype=np.float64)) for D in D_arr],
                [ml.matrix(np.asarray(pv, dtype=np.float64).reshape(1, -1)) for pv in pie_list],
                [ml.matrix(np.asarray(Tm, dtype=np.float64)) for Tm in D0_list],
                'ncMoms', 1)
            for k in range(K):
                Q[ist, k] = float(np.ravel(Qret[k])[0]) if K > 1 else float(np.ravel(Qret)[0])
                R[ist, k] = Q[ist, k] / T[ist, k] if T[ist, k] > 0 else 0.0

    Q = np.abs(Q)
    Q[np.isnan(Q)] = 0
    U[np.isnan(U)] = 0
    R[np.isnan(R)] = 0
    CN = np.sum(R, axis=0, keepdims=True)
    CN[np.isnan(CN)] = 0
    for c in range(C):
        classes_c = np.asarray(inchain[c], dtype=int).ravel()
        X[0, classes_c] = lambda_chain[c] * V[int(np.asarray(sn.refstat).ravel()[classes_c[0]]), classes_c]

    return MAMResult(
        QN=Q,
        UN=U,
        RN=R,
        TN=T,
        CN=CN,
        XN=X,
        totiter=totiter,
        method="mna_open",
        runtime=time.time() - start_time,
    )
