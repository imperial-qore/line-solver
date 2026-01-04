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


class SolverDES:
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
        self._table_silent = False
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
                qlen = result.QN[i, r] if result.QN is not None and i < result.QN.shape[0] and r < result.QN.shape[1] else 0
                util = result.UN[i, r] if result.UN is not None and i < result.UN.shape[0] and r < result.UN.shape[1] else 0
                respt = result.RN[i, r] if result.RN is not None and i < result.RN.shape[0] and r < result.RN.shape[1] else 0
                tput = result.TN[i, r] if result.TN is not None and i < result.TN.shape[0] and r < result.TN.shape[1] else 0
                arvr = result.AN[i, r] if result.AN is not None and i < result.AN.shape[0] and r < result.AN.shape[1] else 0

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

                # Skip zero rows (but not source stations)
                if not is_source and abs(qlen) < 1e-12 and abs(util) < 1e-12 and abs(tput) < 1e-12:
                    continue

                rows.append({
                    'Station': self._station_names[i],
                    'JobClass': self._class_names[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': respt,  # Same as RespT for now
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    # Method aliases (consistent with other solvers)
    get_avg_table = getAvgTable
    avg_table = getAvgTable
    avgT = getAvgTable
    aT = getAvgTable

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
