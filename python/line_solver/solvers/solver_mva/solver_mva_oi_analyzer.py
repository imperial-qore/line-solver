"""
Exact mean-value MVA for order-independent (OI) queueing networks.

An OI station is a class-dependent load-dependent server whose total service
rate mu(n) is a permutation-invariant function of the per-class count vector n.
A closed network of infinite-server (delay) and load-independent (single-server,
product-form) stations plus ANY number of OI stations is product-form. This
analyzer aggregates the delay stations into a single think-time vector Z,
collects the load-independent (LI) queue demands and the OI-station rate handles,
and calls ``pfqn_mvaoi``, the mean-value Conditional-MVA (CMVA) that carries one
rate-shift vector per OI station and returns exact per-class throughput and
queue-lengths WITHOUT any normalizing constant or joint marginal. The
marginal-distribution counterpart is ``pfqn_mvaoi_marg``.

Reference:
    Reiser, Lavenberg (1980). Mean-Value Analysis of Closed Multichain Queuing
    Networks. JACM 27(2). Load-dependent extension: Bruell, Balbo, Afshari
    (1984). OI stations / CMVA: Casale (2009); Casale, Comte, Dorsman (2026).
"""

import time
import numpy as np
from typing import Dict, List

from line_solver.api.pfqn import pfqn_mvaoi


def _is_oi_nodeparam(param) -> bool:
    """True when a node parameter block describes an OI station (PAS/OI with an
    all-zero swap graph and a service-rate function)."""
    if not isinstance(param, dict):
        return False
    if param.get('svcRateFun') is None:
        return False
    sg = param.get('swapGraph')
    if sg is None:
        return False
    sg = np.asarray(sg, dtype=float)
    return sg.size == 0 or not np.any(sg != 0)


def find_oi_stations(sn) -> List[int]:
    """Station indices of all OI stations (PAS/OI with all-zero swap graph and a
    service-rate function). ``sn.nodeparam`` is keyed by node index; map back to a
    station index via ``sn.nodeToStation``."""
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is None:
        return []
    items = nodeparam.items() if isinstance(nodeparam, dict) else enumerate(nodeparam)
    out = []
    for node_idx, param in items:
        if _is_oi_nodeparam(param):
            ist = int(sn.nodeToStation[int(node_idx)])
            if ist >= 0:
                out.append(ist)
    return sorted(out)


def find_oi_station(sn) -> int:
    """Station index of the first OI station, or -1 if none (dispatch helper)."""
    lst = find_oi_stations(sn)
    return lst[0] if lst else -1


class SolverMVAOIAnalyzer:
    """Exact mean-value MVA for order-independent queueing networks."""

    def __init__(self, sn, options: dict = None):
        self.sn = sn
        self.options = options or {}
        self.R = int(sn.nclasses)
        self.M = int(sn.nstations)
        self.njobs = np.asarray(sn.njobs, dtype=float).ravel().astype(int)
        self.oi_list = find_oi_stations(sn)
        if not self.oi_list:
            raise RuntimeError("OI solver requires at least one order-independent station")

        # ---- reject class switching (OI rank rates are per raw class) ------
        # The recursion is driven by the per-class population vector sn.njobs,
        # which class switching makes meaningless: a class that only ever
        # appears mid-chain carries njobs = 0, so the OI station would be
        # analyzed as if empty. Refuse it the way solver_nc_oi_analyzer does
        # rather than return that silently.
        for c in range(sn.nchains):
            if len(np.atleast_1d(sn.inchain[c])) > 1:
                raise RuntimeError('solver_mva_oi requires one class per chain (no class switching).')

        rates = np.asarray(sn.rates, dtype=float)
        visits = np.zeros((self.M, self.R))
        snvisits = getattr(sn, 'visits', None)
        if isinstance(snvisits, dict):
            for c in snvisits:
                visits += np.asarray(snvisits[c], dtype=float)
        else:
            visits[:] = 1.0
        self.rates = rates
        self.visits = visits
        self.demand = np.zeros((self.M, self.R))
        for i in range(self.M):
            for r in range(self.R):
                if np.isfinite(rates[i, r]) and rates[i, r] > 0:
                    self.demand[i, r] = visits[i, r] / rates[i, r]

        self.is_oi = np.zeros(self.M, dtype=bool)
        self.is_oi[self.oi_list] = True

        # OI service-rate handles, one per OI station.
        nodeparam = sn.nodeparam
        self._svc_rate_funs = []
        for oi in self.oi_list:
            node_oi = int(sn.stationToNode[oi])
            param = nodeparam[node_oi]
            self._svc_rate_funs.append(param['svcRateFun'])

        sched = sn.sched
        self.is_delay = np.zeros(self.M, dtype=bool)
        # sn.sched carries lang.base.SchedStrategy members, whose values diverge from
        # the public constants.SchedStrategy for FSP/PAS/OI; use the populating enum.
        from ...lang.base import SchedStrategy as _SS
        for i in range(self.M):
            si = sched[i] if isinstance(sched, dict) else sched[i]
            self.is_delay[i] = (int(getattr(si, 'value', si)) == int(_SS.INF.value))

    def _make_mu(self, fun):
        """Wrap a microstate service-rate handle into a count-vector rate function."""
        R = self.R

        def mu(n):
            cls = np.repeat(np.arange(R), np.asarray(n, dtype=int))
            if cls.size == 0:
                return 0.0
            return float(fun(cls))
        return mu

    @staticmethod
    def _make_ms_mu(Dq, c):
        """OI rate function reproducing the c-server BCMP station with per-class
        demands Dq: mu(n) = (min(|n|,c)/|n|) * sum_{r: n_r>0} n_r/Dq_r. For c = 1
        this is the total completion rate of a multiclass single-server queue, and
        for a single class it reduces to min(n,c)/Dq (M/M/c)."""
        Dq = np.asarray(Dq, dtype=float).ravel()
        if not np.isfinite(c) or c <= 0:
            c = 1.0

        def mu(n):
            nn = np.asarray(n, dtype=float).ravel()
            tot = float(nn.sum())
            if tot == 0:
                return 0.0
            ok = (nn > 0) & (Dq > 0)
            acc = float(np.sum(nn[ok] / Dq[ok]))
            return (min(tot, c) / tot) * acc
        return mu

    def analyze(self) -> Dict:
        tstart = time.time()
        R, M = self.R, self.M
        N = self.njobs

        # see _kb/06-solver-catalog.md ("OI (order-independent) utilization")
        # for the c-server BCMP -> OI promotion derivation
        nservers_all = np.asarray(self.sn.nservers, dtype=float).ravel()
        Z = np.zeros(R)
        li_list = []
        ms_list = []
        for i in range(M):
            if self.is_oi[i]:
                continue
            elif self.is_delay[i]:
                Z += self.demand[i, :]
            elif np.isfinite(nservers_all[i]) and nservers_all[i] > 1:
                ms_list.append(i)
            else:
                li_list.append(i)
        Dli = self.demand[li_list, :] if li_list else np.zeros((0, R))
        muCell = [self._make_mu(f) for f in self._svc_rate_funs]
        for i in ms_list:
            muCell.append(self._make_ms_mu(self.demand[i, :], nservers_all[i]))

        # Per-muCell-station visit vectors: genuine OI stations carry their
        # class visits (rate handle has none); ms-promoted stations pass unit
        # visits (already folded into demand by _make_ms_mu).
        oivis = [np.asarray(self.visits[oi, :], dtype=float) for oi in self.oi_list]
        oivis += [np.ones(R) for _ in ms_list]

        X, Qoi, Qli, _, Soi = pfqn_mvaoi(Z, N, muCell, Dli, oivis)

        QN = np.zeros((M, R))
        for o, i in enumerate(self.oi_list):
            QN[i, :] = Qoi[o, :]
        for j, i in enumerate(ms_list):
            QN[i, :] = Qoi[len(self.oi_list) + j, :]
        for j, i in enumerate(li_list):
            QN[i, :] = Qli[j, :]
        for i in range(M):
            if self.is_delay[i]:
                QN[i, :] = X * self.demand[i, :]
        XN = np.asarray(X, dtype=float)

        # Row of Soi/Qoi holding each OI station (muCell order: oi_list, ms_list).
        oi_row = {i: o for o, i in enumerate(self.oi_list)}

        TN = np.zeros((M, R))
        for i in range(M):
            for r in range(R):
                TN[i, r] = XN[r] * self.visits[i, r]
        RN = np.zeros((M, R))
        UN = np.zeros((M, R))
        CN = np.zeros((M, R))
        nservers = np.asarray(self.sn.nservers, dtype=float).ravel()
        for i in range(M):
            for r in range(R):
                if XN[r] > 0:
                    RN[i, r] = QN[i, r] / XN[r]
                if self.is_oi[i]:
                    # see _kb/06-solver-catalog.md ("OI (order-independent) utilization")
                    sv = nservers[i] if np.isfinite(nservers[i]) and nservers[i] > 0 else 1.0
                    UN[i, r] = Soi[oi_row[i], r] / sv
                elif self.is_delay[i]:
                    UN[i, r] = QN[i, r]
                else:
                    sv = nservers[i] if np.isfinite(nservers[i]) and nservers[i] > 0 else 1.0
                    UN[i, r] = XN[r] * self.demand[i, r] / sv
                CN[i, r] = RN[i, r]

        return {
            'QN': QN, 'UN': UN, 'RN': RN, 'TN': TN, 'CN': CN, 'XN': XN,
            'lG': 0.0, 'runtime': time.time() - tstart,
            'iter': int(np.sum(N)), 'method': 'oi',
        }
