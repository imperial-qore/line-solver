"""Mean-measure presentation shared by the decomposition solvers.

These methods read ``self.result``, ``self.sn`` and ``self.options`` and nothing
else -- they format what an analysis produced and never choose or run one. They
live here because SolverMAM and SolverAG report the same measures out of the same
result container while sharing no algorithm at all: RCAT decomposes the model
into cooperating agents, the MAM methods decompose its traffic, and the only
thing they genuinely have in common is the shape of the answer.

Keep it that way. Anything that inspects a method name, selects an algorithm or
reaches into a solver-specific result field belongs on the solver, not here.

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

from .base import avg_table_drop_empty_rows, method_label, method_type
from ..api.sn.getters import sn_get_arvr_from_tput


class AvgResultsMixin:
    """Mean-measure getters over a solved result container."""

    def getAvgTable(self) -> pd.DataFrame:
        """Get average performance metrics as DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self.result is None:
            self._ensureAvgResults()
        # The saturated open-station rule of getAvg.m, which every solver owes.
        # WN is read off the analyzer's own response times below, so the
        # residence column keeps the pre-saturation value as the reference does.
        cap = getattr(self, '_cap_unstable_open_util', None)
        if callable(cap):
            cap()

        M = self.result.QN.shape[0]
        K = self.result.QN.shape[1]

        # Get station names using stationToNode mapping
        nodenames = list(self.sn.nodenames) if hasattr(self.sn, 'nodenames') and self.sn.nodenames else []
        stationToNode = self.sn.stationToNode if hasattr(self.sn, 'stationToNode') else None

        station_names = []
        if stationToNode is not None and nodenames:
            stationToNode = np.asarray(stationToNode).flatten()
            for i in range(M):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        station_names.append(nodenames[node_idx])
                    else:
                        station_names.append(f'Station{i}')
                else:
                    station_names.append(f'Station{i}')
        else:
            station_names = [f'Station{i}' for i in range(M)]

        # Identify Source stations
        source_stations = set()
        if hasattr(self.sn, 'sched') and self.sn.sched is not None:
            for ist in range(M):
                sched = self.sn.sched.get(ist, None)
                if sched is not None:
                    sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                    if sched_name == 'EXT' or (hasattr(sched, 'value') and sched.value == 11):
                        source_stations.add(ist)

        # Get class names
        class_names = list(self.sn.classnames) if hasattr(self.sn, 'classnames') and self.sn.classnames else []

        # Build rows (one per station per class)
        rows = []
        for i in range(M):
            is_source = i in source_stations
            for r in range(K):
                class_name = class_names[r] if r < len(class_names) else f'Class{r}'

                qlen = float(self.result.QN[i, r])
                util = float(self.result.UN[i, r])
                respt = float(self.result.RN[i, r])
                residt = float(self.result.WN[i, r]) if hasattr(self.result, 'WN') and self.result.WN is not None else respt

                if self.result.TN.ndim > 1:
                    if self.result.TN.shape[0] == 1:
                        tput = float(self.result.TN[0, r])
                    else:
                        tput = float(self.result.TN[i, r])
                else:
                    tput = float(self.result.TN[r]) if r < len(self.result.TN) else 0.0

                if hasattr(self.result, 'AN') and self.result.AN is not None and i < self.result.AN.shape[0]:
                    arvr = float(self.result.AN[i, r])
                elif is_source:
                    arvr = 0.0
                else:
                    arvr = tput

                rows.append({
                    'Station': station_names[i],
                    'JobClass': class_name,
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = avg_table_drop_empty_rows(pd.DataFrame(rows))

        if not self._table_silent:
            print(df.to_string(index=False))

        from ..indexed_table import IndexedTable
        return IndexedTable(df)

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths per station and class.

        Returns:
            (M, K) array of average queue lengths
        """
        if self.result is None:
            self._ensureAvgResults()
        return self._avgStationClass('Q')

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations per station and class.

        Returns:
            (M, K) array of utilizations
        """
        if self.result is None:
            self._ensureAvgResults()
        cap = getattr(self, '_cap_unstable_open_util', None)
        if callable(cap):
            cap()
        return self._avgStationClass('U')

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times per station and class.

        Returns:
            (M, K) array of response times
        """
        if self.result is None:
            self._ensureAvgResults()
        return self._avgStationClass('R')

    def getAvgSysRespT(self) -> np.ndarray:
        """Get average system response time per CHAIN.

        Returns:
            (C,) array of system response times, one per chain.

        This is CNchain of MATLAB @NetworkSolver/getAvgSys.m: an open chain sums
        the visit-weighted per-class residence times, a closed one applies
        Little's law to the chain population. It is per chain and not per class
        because a chain is what completes a cycle -- a job that enters as class A
        and leaves as class B has one system response time, not two.
        """
        if self.result is None:
            self._ensureAvgResults()
        CN, _ = self._computeChainMetrics()
        return CN

    def getAvgSysTput(self) -> np.ndarray:
        """Get average system throughput per CHAIN.

        Returns:
            (C,) array of system throughputs, one per chain.

        XNchain of MATLAB @NetworkSolver/getAvgSys.m: the rate of completing
        classes routed back into the chain's reference station. This used to
        return mean(XN), one number over the classes, which is the throughput of
        nothing; the per-class reference-station rate it was then briefly changed
        to is a different measure that the reference does not publish here.
        """
        if self.result is None:
            self._ensureAvgResults()
        _, XN = self._computeChainMetrics()
        return XN

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T), both (C,): per-chain system response time and
            per-chain system throughput.

        It lives here rather than on each solver because it is nothing but the
        pair of the two accessors above, and a solver that assembled the pair
        itself got it wrong: SolverMAM read getTput(), the per station-class
        matrix, which getAvgSysTable then flattened and indexed one entry per
        class -- so it published the FIRST STATION'S row as the system rate,
        right only where station 0 happens to be every class's reference station
        and a clean factor of two out on a station visited twice per cycle.
        SolverAG, which assembled no pair at all, had no getAvgSys and raised
        AttributeError from getAvgSysTable. Both halves now come from the one
        chain computation NetworkSolver ports from the reference.
        """
        return self.getAvgSysRespT(), self.getAvgSysTput()

    # =====================================================================
    # STATIC METHODS (introspection and validation)
    # =====================================================================

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times per station.

        Residence time = Response time * visits / ref_visits
        This accounts for multiple visits to the same station.

        Returns:
            (M, K) array of residence times
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")
        if hasattr(self.result, 'WN') and self.result.WN is not None:
            return self.result.WN
        # Fallback: compute on demand if not already computed
        return sn_get_residt_from_respt(self.sn, self.result.RN, None)

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times per station and class.

        Waiting time is response time minus the mean service time.

        Returns:
            (M, K) array of waiting times
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")
        return self._avgWaitTFromRespT(self._avgStationClass('R'))

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates per station and class.

        Returns:
            (M, K) array of arrival rates
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")
        return sn_get_arvr_from_tput(self.sn, self._avgStationClass('T'))

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs per station and class.

        Returns:
            (M, K) array of throughputs
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")
        return self._avgStationClass('T')

    # =====================================================================
    # SAMPLING METHODS (Not Supported - Analytical Solver)
    # =====================================================================

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self.result is None:
            self._ensureAvgResults()

        Q = self.result.QN
        chains = self._get_chains()
        nstations = Q.shape[0]
        nchains = len(chains)

        QN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                QN_chain[:, c] = np.sum(Q[:, chain_classes], axis=1)

        return QN_chain

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        if self.result is None:
            self._ensureAvgResults()

        U = self.result.UN
        chains = self._get_chains()
        nstations = U.shape[0]
        nchains = len(chains)

        UN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                UN_chain[:, c] = np.sum(U[:, chain_classes], axis=1)

        return UN_chain

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain."""
        if self.result is None:
            self._ensureAvgResults()

        R = self.result.RN
        chains = self._get_chains()
        nstations = R.shape[0]
        nchains = len(chains)

        RN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                # Weighted average by throughput
                RN_chain[:, c] = np.mean(R[:, chain_classes], axis=1)

        return RN_chain

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self.result is None:
            self._ensureAvgResults()

        T = self.result.TN
        if T.ndim == 1:
            T = T.reshape(1, -1)

        chains = self._get_chains()
        nstations = self.result.QN.shape[0]
        nchains = len(chains)

        TN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                if T.shape[0] == nstations:
                    TN_chain[:, c] = np.sum(T[:, chain_classes], axis=1)
                else:
                    TN_chain[:, c] = np.sum(T[0, chain_classes])

        return TN_chain

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics aggregated by chain.

        Returns:
            Tuple of (QN, UN, RN, WN, AN, TN) aggregated by chain
        """
        QN = self.getAvgQLenChain()
        UN = self.getAvgUtilChain()
        RN = self.getAvgRespTChain()
        WN = self.getAvgResidTChain()
        AN = self.getAvgArvRChain()
        TN = self.getAvgTputChain()
        return QN, UN, RN, WN, AN, TN

    # Legacy CamelCase and snake_case spellings, kept because published
    # examples and the MATLAB-facing bridge call them.
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput
    GetAvgSys = getAvgSys
    GetAvgTable = getAvgTable
    GetAvgChain = getAvgChain
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain
    aT = getAvgTable
    avgT = getAvgTable
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_respt = getAvgRespT
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
