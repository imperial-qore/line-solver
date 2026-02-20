"""
Native Python DES solver using SimPy.

This module provides the SolverDES class that implements discrete
event simulation for queueing networks using SimPy.
"""

import time
import numpy as np
import pandas as pd
from typing import Optional, Any, Dict

from .des_options import DESOptions, DESResult
from .simulator import SimPySimulator
from ..base import NetworkSolver
from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn.network_struct import NodeType
from ...constants import GlobalConstants


class SolverDES(NetworkSolver):
    """
    Native Python Discrete Event Simulation (DES) solver using SimPy.

    This solver implements DES algorithms using pure Python/SimPy,
    providing the same functionality as the Java DES wrapper without
    requiring the JVM.

    Args:
        model: Network model (Python wrapper or native structure)
        options: DESOptions configuration (optional)
        **kwargs: Additional solver options (seed, samples, verbose, etc.)

    Example:
        >>> from line_solver.solvers.solver_des import SolverDES, DESOptions
        >>> options = DESOptions(seed=23000, samples=200000)
        >>> solver = SolverDES(model, options)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()
    """

    def __init__(
        self,
        model: Any,
        options: Optional[DESOptions] = None,
        **kwargs
    ):
        self.model = model
        self.options = options or DESOptions()

        # Apply kwargs to options
        for key, value in kwargs.items():
            if hasattr(self.options, key):
                setattr(self.options, key, value)

        self._result: Optional[DESResult] = None
        self._simulator: Optional[SimPySimulator] = None

        # Extract network structure
        self._sn = self._get_network_struct()

        # Station and class names for table output
        self._station_names: list = []
        self._class_names: list = []
        self._extract_names()

    def _get_network_struct(self) -> Any:
        """Get NetworkStruct from model."""
        model = self.model

        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            return model._sn

        # Priority 2: Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                return model._sn

        # Priority 3: JPype wrapper with getStruct
        if hasattr(model, 'getStruct'):
            try:
                return model.getStruct()
            except Exception:
                pass

        # Priority 4: JPype wrapper with obj attribute
        if hasattr(model, 'obj'):
            try:
                return model.getStruct()
            except Exception:
                pass

        # Priority 5: Model that is already a struct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError("Cannot extract network structure from model")

    def _extract_names(self) -> None:
        """Extract station and class names from network structure."""
        sn = self._sn
        num_stations = int(sn.nstations) if hasattr(sn, 'nstations') else 1

        # Station names - use stationToNode mapping if available
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames is not None else []
        stationToNode = sn.stationToNode if hasattr(sn, 'stationToNode') and sn.stationToNode is not None else None

        if stationToNode is not None and nodenames:
            import numpy as np
            stationToNode = np.asarray(stationToNode).flatten()
            self._station_names = []
            for i in range(num_stations):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        self._station_names.append(nodenames[node_idx])
                    else:
                        self._station_names.append(f'Station{i}')
                else:
                    self._station_names.append(f'Station{i}')
        elif nodenames:
            # Fallback: use first num_stations nodenames
            self._station_names = nodenames[:num_stations]
            while len(self._station_names) < num_stations:
                self._station_names.append(f'Station{len(self._station_names)}')
        else:
            self._station_names = [f'Station{i}' for i in range(num_stations)]

        # Class names
        if hasattr(sn, 'classnames') and sn.classnames is not None:
            try:
                self._class_names = list(sn.classnames)
            except Exception:
                self._class_names = [f'Class{i}' for i in range(sn.nclasses)]
        else:
            num_classes = int(sn.nclasses) if hasattr(sn, 'nclasses') else 1
            self._class_names = [f'Class{i}' for i in range(num_classes)]

    def runAnalyzer(self) -> DESResult:
        """
        Run the DES simulation.

        Returns:
            DESResult containing performance metrics
        """
        start_time = time.time()

        # Create simulator
        self._simulator = SimPySimulator(
            sn=self._sn,
            options=self.options,
            model=self.model,
        )

        # Run simulation
        max_events = self.options.samples
        self._simulator.simulate(max_events)

        # Get results
        self._result = self._simulator.get_result()
        self._result.runtime = time.time() - start_time

        return self._result

    def run_analyzer(self) -> DESResult:
        """Alias for runAnalyzer (Python convention)."""
        return self.runAnalyzer()

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics as a DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        M = len(self._station_names)
        K = len(self._class_names)
        sn = self._sn

        # Make copies of station-level metrics (to avoid modifying originals)
        QN = result.QN.copy() if result.QN is not None else np.zeros((M, K))
        UN = result.UN.copy() if result.UN is not None else np.zeros((M, K))
        RN = result.RN.copy() if result.RN is not None else np.zeros((M, K))
        TN = result.TN.copy() if result.TN is not None else np.zeros((M, K))
        AN = result.AN.copy() if result.AN is not None else np.zeros((M, K))

        # Zero out metrics for classes that don't visit stations based on visit ratios
        # This matches MATLAB getAvg.m behavior (lines 163-180)
        hasForkJoin = False
        hasSPN = False
        if hasattr(sn, 'nodetype') and sn.nodetype is not None:
            hasForkJoin = np.any(sn.nodetype == NodeType.FORK) and np.any(sn.nodetype == NodeType.JOIN)
            hasSPN = np.any(sn.nodetype == NodeType.PLACE) or np.any(sn.nodetype == NodeType.TRANSITION)

        if sn is not None and hasattr(sn, 'nchains') and sn.nchains > 0 and not hasSPN:
            if hasattr(sn, 'chains') and sn.chains is not None and hasattr(sn, 'visits') and sn.visits:
                chains_arr = np.asarray(sn.chains)
                for k in range(K):
                    # Find chains containing this class
                    chains_with_class = np.where(chains_arr[:, k] > 0)[0] if k < chains_arr.shape[1] else []
                    if len(chains_with_class) > 0:
                        c = chains_with_class[0]  # Use first chain (classes typically in one chain)
                        if c in sn.visits and sn.visits[c] is not None:
                            visits_c = np.asarray(sn.visits[c])
                            for i in range(M):
                                if i < visits_c.shape[0] and k < visits_c.shape[1]:
                                    if visits_c[i, k] == 0:
                                        # For fork-join, trust non-zero simulation results
                                        if hasForkJoin and (QN[i, k] > GlobalConstants.FineTol or
                                                           UN[i, k] > GlobalConstants.FineTol or
                                                           TN[i, k] > GlobalConstants.FineTol):
                                            continue
                                        # Zero out station-level metrics
                                        QN[i, k] = 0
                                        UN[i, k] = 0
                                        RN[i, k] = 0
                                        TN[i, k] = 0
                                        AN[i, k] = 0

        # Compute residence times from response times using visit ratios (after zeroing)
        WN = None
        if RN is not None and sn is not None and hasattr(sn, 'visits') and sn.visits:
            try:
                WN = sn_get_residt_from_respt(sn, RN, None)
            except Exception:
                # Fallback: ResidT = RespT if computation fails
                WN = RN.copy() if RN is not None else None

        # Identify source stations using stationToNode mapping and nodetype
        source_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None
        if nodetype is not None:
            stationToNode = self._sn.stationToNode
            if stationToNode is not None:
                stationToNode = np.asarray(stationToNode).flatten()
                nodetype = np.asarray(nodetype).flatten()
                for i in range(M):
                    if i < len(stationToNode):
                        node_idx = int(stationToNode[i])
                        if node_idx < len(nodetype):
                            if int(nodetype[node_idx]) == 0:  # SOURCE = 0
                                source_stations.add(i)

        rows = []
        for i in range(M):
            for r in range(K):
                qlen = QN[i, r] if i < QN.shape[0] and r < QN.shape[1] else 0
                util = UN[i, r] if i < UN.shape[0] and r < UN.shape[1] else 0
                respt = RN[i, r] if i < RN.shape[0] and r < RN.shape[1] else 0
                tput = TN[i, r] if i < TN.shape[0] and r < TN.shape[1] else 0
                arvr = AN[i, r] if i < AN.shape[0] and r < AN.shape[1] else 0

                is_source = i in source_stations

                # For source stations, set Tput from arrival rate in network struct
                if is_source:
                    # Get arrival rate from rates matrix (sources have their rate there)
                    if hasattr(self._sn, 'rates') and self._sn.rates is not None:
                        rates = np.asarray(self._sn.rates)
                        stationToNode = np.asarray(self._sn.stationToNode).flatten()
                        node_idx = int(stationToNode[i])
                        if node_idx < rates.shape[0] and r < rates.shape[1]:
                            tput = rates[node_idx, r]
                    arvr = 0.0  # Source has no arrivals to itself

                # Skip rows where all metrics are zero (matches MATLAB getAvgTable behavior)
                if abs(qlen) < 1e-12 and abs(util) < 1e-12 and abs(tput) < 1e-12:
                    continue

                # Get residence time from WN if available, otherwise use RespT
                residt = respt
                if WN is not None and i < WN.shape[0] and r < WN.shape[1]:
                    residt_val = WN[i, r]
                    if not np.isnan(residt_val) and residt_val >= 0:
                        residt = residt_val

                rows.append({
                    'Station': self._station_names[i],
                    'JobClass': self._class_names[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    def getAvgChainTable(self) -> pd.DataFrame:
        """
        Get average performance metrics aggregated by chain.

        Returns:
            DataFrame with columns: Chain, QLen, Util, RespT, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result

        # Get chain information from model structure
        nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else self._sn.nclasses
        inchain = self._sn.inchain if hasattr(self._sn, 'inchain') else None

        rows = []
        for c in range(nchains):
            chain_name = f'Chain{c+1}'

            # Get classes in this chain
            if inchain is not None and c in inchain:
                chain_classes = np.asarray(inchain[c]).flatten().astype(int)
            else:
                chain_classes = [c]  # Single class per chain

            # Aggregate metrics across stations and classes in chain
            total_qlen = 0.0
            total_util = 0.0
            total_respt = 0.0
            total_tput = 0.0

            M = self._sn.nstations
            for i in range(M):
                for k in chain_classes:
                    if k < result.QN.shape[1] if result.QN is not None else 0:
                        total_qlen += result.QN[i, k] if result.QN is not None and not np.isnan(result.QN[i, k]) else 0.0
                        total_util += result.UN[i, k] if result.UN is not None and not np.isnan(result.UN[i, k]) else 0.0
                        total_respt += result.RN[i, k] if result.RN is not None and not np.isnan(result.RN[i, k]) else 0.0
                        total_tput = max(total_tput, result.TN[i, k] if result.TN is not None and not np.isnan(result.TN[i, k]) else 0.0)

            rows.append({
                'Chain': chain_name,
                'QLen': total_qlen,
                'Util': total_util,
                'RespT': total_respt,
                'Tput': total_tput,
            })

        return pd.DataFrame(rows)

    def getAvgSysTable(self) -> pd.DataFrame:
        """
        Get system-level average performance metrics.

        Returns:
            DataFrame with columns: Chain, SysRespT, SysTput
        """
        if self._result is None:
            self.runAnalyzer()

        chain_table = self.getAvgChainTable()
        rows = []
        for _, row in chain_table.iterrows():
            rows.append({
                'Chain': row['Chain'],
                'SysRespT': row['RespT'],
                'SysTput': row['Tput'],
            })

        return pd.DataFrame(rows)

    # Method aliases (consistent with other solvers)
    avgT = getAvgTable
    aT = getAvgTable
    get_avg_chain_table = getAvgChainTable
    avg_chain_table = getAvgChainTable
    aCT = getAvgChainTable
    chainAvgT = getAvgChainTable
    get_avg_sys_table = getAvgSysTable
    avg_sys_table = getAvgSysTable

    @property
    def result(self) -> Optional[DESResult]:
        """Get the DES result (after runAnalyzer is called)."""
        return self._result

    def getName(self) -> str:
        """Get solver name."""
        return "DES"

    def get_name(self) -> str:
        """Get solver name (Python convention)."""
        return self.getName()

    @staticmethod
    def defaultOptions() -> DESOptions:
        """Get default DES solver options."""
        return DESOptions()

    @staticmethod
    def default_options() -> DESOptions:
        """Get default options (Python convention)."""
        return DESOptions()

    # ==================== Transient Analysis Methods ====================

    def getTranCdfRespT(self, station_idx: Optional[int] = None, class_idx: Optional[int] = None):
        """
        Get empirical CDF of response times from DES samples.

        Returns a dict with:
            - 'x': list of response time values (sorted)
            - 'cdf': list of cumulative probability values
            - 'station': station index (if specified)
            - 'class': class index (if specified)

        If station_idx and class_idx are specified, returns CDF for that combination.
        If not specified, returns dict of CDFs for all station/class combinations.

        Args:
            station_idx: Optional station index to filter
            class_idx: Optional class index to filter

        Returns:
            Dict with CDF data
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        if result.resp_time_samples is None:
            return {'x': [], 'cdf': [], 'error': 'No response time samples collected'}

        if station_idx is not None and class_idx is not None:
            # Return CDF for specific station/class
            samples = result.resp_time_samples.get(station_idx, {}).get(class_idx, [])
            x, cdf = self._compute_empirical_cdf(samples)
            return {'x': x, 'cdf': cdf, 'station': station_idx, 'class': class_idx}
        else:
            # Return CDFs for all combinations
            cdfs = {}
            for st_idx, class_samples in result.resp_time_samples.items():
                cdfs[st_idx] = {}
                for cl_idx, samples in class_samples.items():
                    x, cdf = self._compute_empirical_cdf(samples)
                    cdfs[st_idx][cl_idx] = {'x': x, 'cdf': cdf}
            return cdfs

    def get_tran_cdf_resp_t(self, station_idx: Optional[int] = None, class_idx: Optional[int] = None):
        """Alias for getTranCdfRespT (Python convention)."""
        return self.getTranCdfRespT(station_idx, class_idx)

    def getTranCdfPassT(self, station_idx: Optional[int] = None, class_idx: Optional[int] = None):
        """
        Get empirical CDF of passage times from DES samples.

        For now, passage times are the same as response times.

        Args:
            station_idx: Optional station index to filter
            class_idx: Optional class index to filter

        Returns:
            Dict with CDF data
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        if result.pass_time_samples is None:
            return {'x': [], 'cdf': [], 'error': 'No passage time samples collected'}

        if station_idx is not None and class_idx is not None:
            samples = result.pass_time_samples.get(station_idx, {}).get(class_idx, [])
            x, cdf = self._compute_empirical_cdf(samples)
            return {'x': x, 'cdf': cdf, 'station': station_idx, 'class': class_idx}
        else:
            cdfs = {}
            for st_idx, class_samples in result.pass_time_samples.items():
                cdfs[st_idx] = {}
                for cl_idx, samples in class_samples.items():
                    x, cdf = self._compute_empirical_cdf(samples)
                    cdfs[st_idx][cl_idx] = {'x': x, 'cdf': cdf}
            return cdfs

    def get_tran_cdf_pass_t(self, station_idx: Optional[int] = None, class_idx: Optional[int] = None):
        """Alias for getTranCdfPassT (Python convention)."""
        return self.getTranCdfPassT(station_idx, class_idx)

    def getTranProb(self, station_idx: Optional[int] = None):
        """
        Get transient state probabilities for a station.

        Note: The Python SimPy DES solver does not currently track full state
        probabilities during simulation. This method returns queue length
        distributions derived from transient data if available.

        Args:
            station_idx: Station index (optional)

        Returns:
            Dict with transient probability data or message indicating unavailability
        """
        if self._result is None:
            self.runAnalyzer()

        # Check if transient data is available
        if self._result.QNt is None or self._result.t is None:
            return {
                'error': 'Transient probability data not available. '
                         'Enable transient mode with options.transient=True'
            }

        result = self._result
        t = result.t.tolist() if result.t is not None else []

        if station_idx is not None:
            # Return data for specific station
            qlen_over_time = result.QNt[:, station_idx, :].tolist() if result.QNt is not None else []
            return {
                't': t,
                'qlen': qlen_over_time,
                'station': station_idx
            }
        else:
            # Return data for all stations
            data = {}
            num_stations = result.QNt.shape[1] if result.QNt is not None else 0
            for st_idx in range(num_stations):
                qlen_over_time = result.QNt[:, st_idx, :].tolist()
                data[st_idx] = {'t': t, 'qlen': qlen_over_time}
            return data

    def get_tran_prob(self, station_idx: Optional[int] = None):
        """Alias for getTranProb (Python convention)."""
        return self.getTranProb(station_idx)

    def getTranProbAggr(self, station_idx: Optional[int] = None):
        """
        Get aggregated transient state probabilities for a station.

        Returns aggregated queue length distribution over time.

        Args:
            station_idx: Station index (optional)

        Returns:
            Dict with aggregated transient probability data
        """
        # For DES, aggregated is sum across classes
        if self._result is None:
            self.runAnalyzer()

        if self._result.QNt is None or self._result.t is None:
            return {
                'error': 'Transient probability data not available. '
                         'Enable transient mode with options.transient=True'
            }

        result = self._result
        t = result.t.tolist() if result.t is not None else []

        if station_idx is not None:
            # Aggregate across classes
            qlen_aggr = np.sum(result.QNt[:, station_idx, :], axis=1).tolist()
            return {
                't': t,
                'qlen_aggr': qlen_aggr,
                'station': station_idx,
                'isAggregated': True
            }
        else:
            data = {}
            num_stations = result.QNt.shape[1] if result.QNt is not None else 0
            for st_idx in range(num_stations):
                qlen_aggr = np.sum(result.QNt[:, st_idx, :], axis=1).tolist()
                data[st_idx] = {'t': t, 'qlen_aggr': qlen_aggr, 'isAggregated': True}
            return data

    def get_tran_prob_aggr(self, station_idx: Optional[int] = None):
        """Alias for getTranProbAggr (Python convention)."""
        return self.getTranProbAggr(station_idx)

    def getTranProbSys(self):
        """
        Get transient system-level state probabilities.

        Returns system-wide queue length distribution over time.

        Returns:
            Dict with system-level transient probability data
        """
        if self._result is None:
            self.runAnalyzer()

        if self._result.QNt is None or self._result.t is None:
            return {
                'error': 'Transient probability data not available. '
                         'Enable transient mode with options.transient=True'
            }

        result = self._result
        t = result.t.tolist() if result.t is not None else []

        # Sum across all stations and classes
        qlen_sys = np.sum(np.sum(result.QNt, axis=2), axis=1).tolist()
        return {
            't': t,
            'qlen_sys': qlen_sys,
            'isSystemLevel': True
        }

    def get_tran_prob_sys(self):
        """Alias for getTranProbSys (Python convention)."""
        return self.getTranProbSys()

    def getTranProbSysAggr(self):
        """
        Get aggregated system-level transient state probabilities.

        Returns aggregated system-wide distribution over time.

        Returns:
            Dict with aggregated system-level transient probability data
        """
        # For DES, this is the same as getTranProbSys
        result = self.getTranProbSys()
        if 'error' not in result:
            result['isAggregated'] = True
        return result

    def get_tran_prob_sys_aggr(self):
        """Alias for getTranProbSysAggr (Python convention)."""
        return self.getTranProbSysAggr()

    def getTranAvg(self):
        """
        Get transient averages (time-varying performance metrics).

        Returns performance metrics indexed by time if transient mode was enabled.

        Returns:
            Dict with transient average data
        """
        if self._result is None:
            self.runAnalyzer()

        result = self._result
        if result.t is None or result.QNt is None:
            return {
                'error': 'Transient data not available. '
                         'Enable transient mode with options.transient=True or options.timespan=(0, T)'
            }

        return {
            't': result.t.tolist() if result.t is not None else [],
            'QNt': result.QNt.tolist() if result.QNt is not None else [],
            'UNt': result.UNt.tolist() if result.UNt is not None else [],
            'TNt': result.TNt.tolist() if result.TNt is not None else [],
        }

    def get_tran_avg(self):
        """Alias for getTranAvg (Python convention)."""
        return self.getTranAvg()

    def _compute_empirical_cdf(self, samples: list) -> tuple:
        """
        Compute empirical CDF from samples.

        Args:
            samples: List of sample values

        Returns:
            Tuple of (x_values, cdf_values) as lists
        """
        if not samples:
            return [], []

        sorted_samples = np.sort(samples)
        n = len(sorted_samples)
        cdf = np.arange(1, n + 1) / n

        return sorted_samples.tolist(), cdf.tolist()
