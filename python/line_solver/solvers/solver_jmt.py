"""
JMT solver integration.

This implementation calls JMT via subprocess (command line), matching how
MATLAB's SolverJMT works. No JPype/JVM integration in Python itself.

The solver:
1. Writes the model to JSIM/JMVA XML format
2. Calls JMT via command line
3. Parses the result XML file

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import numpy as np
import pandas as pd
import os
import tempfile
import shutil
from typing import Optional, Dict, Any, List, Tuple, Set


class OptionsDict(dict):
    """A dict that supports attribute-style access."""
    def __getattr__(self, name):
        try:
            return self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")

    def __setattr__(self, name, value):
        self[name] = value

    def __delattr__(self, name):
        try:
            del self[name]
        except KeyError:
            raise AttributeError(f"'OptionsDict' object has no attribute '{name}'")
from dataclasses import dataclass

from ..api.solvers.jmt.handler import (
    solver_jmt,
    SolverJMTOptions as _SolverJMTOptions,
    SolverJMTReturn,
    is_jmt_available,
    _get_jmt_jar_path,
)
from .convert import wrapper_sn_to_native


@dataclass
class SolverJMTOptions:
    """Options for the native JMT solver."""
    method: str = 'jsim'
    samples: int = 10000
    seed: int = 23000
    max_simulated_time: float = float('inf')
    conf_int: float = 0.99
    max_rel_err: float = 0.03
    verbose: bool = False
    keep: bool = False  # Keep temp files after execution


class SolverJMT:
    """
    JMT solver integration.

    This solver provides discrete-event simulation and analytical methods
    via command line
    is launched as an external process, exactly like MATLAB's SolverJMT.

    Supported methods:
        - 'jsim' / 'default': Discrete event simulation
        - 'jmva' / 'jmva.mva': Mean Value Analysis
        - 'jmva.amva': Approximate MVA
        - 'jmva.recal': RECALsimulation
        - 'jmva.comom': CoMoM algorithm
        - 'jmva.chow': Chow algorithm
        - 'jmva.bs': Bard-Schweitzer
        - 'jmva.aql': AQL algorithm
        - 'jmva.lin': Linearizer
        - 'jmva.dmlin': De Souza-Muntz Linearizer

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'jsim')
        **kwargs: Additional solver options (samples, seed, etc.)

    Example:
        >>> solver = SolverJMT(model, samples=10000, seed=42)
        >>> solver.runAnalyzer()
        >>> table = solver.getAvgTable()
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            self.method = 'jsim'
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'jsim')
            if 'samples' in method_or_options:
                kwargs.setdefault('samples', method_or_options['samples'])
            if 'seed' in method_or_options:
                kwargs.setdefault('seed', method_or_options['seed'])
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'jsim')
            if hasattr(method_or_options, 'samples'):
                kwargs.setdefault('samples', method_or_options.samples)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
        else:
            self.method = 'jsim'

        # Parse options
        samples = kwargs.get('samples', 10000)
        seed = kwargs.get('seed', 23000)
        verbose = kwargs.get('verbose', False)
        keep = kwargs.get('keep', False)
        conf_int = kwargs.get('conf_int', kwargs.get('confint', 0.99))
        max_rel_err = kwargs.get('max_rel_err', 0.03)
        max_simulated_time = kwargs.get('max_simulated_time',
                                        kwargs.get('timespan', [0, float('inf')])[1]
                                        if isinstance(kwargs.get('timespan'), list) else float('inf'))

        self.options = SolverJMTOptions(
            method=self.method,
            samples=samples,
            seed=seed,
            max_simulated_time=max_simulated_time,
            conf_int=conf_int,
            max_rel_err=max_rel_err,
            verbose=verbose,
            keep=keep
        )

        self._result: Optional[SolverJMTReturn] = None
        self._sn = None
        self._table_silent = False

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "JMT"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model."""
        model = self.model

        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            self._sn = model._sn
            return

        # Priority 2: Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                self._sn = model._sn
                return

        # Priority 3: JPype wrapper with getStruct()
        if hasattr(model, 'getStruct'):
            try:
                sn = model.getStruct()
                self._sn = wrapper_sn_to_native(sn)
                return
            except Exception:
                pass

        # Priority 4: JPype wrapper with obj attribute
        if hasattr(model, 'obj'):
            try:
                sn = model.getStruct()
                self._sn = wrapper_sn_to_native(sn)
                return
            except Exception:
                pass

        raise ValueError("Cannot extract network structure from model")

    def runAnalyzer(self) -> 'SolverJMT':
        """
        Run the JMT analyzer.

        Calls JMT via command line and stores the results.

        Returns:
            self for method chaining
        """
        if self._sn is None:
            raise RuntimeError("Network structure not available")

        # Convert options to handler format
        handler_options = _SolverJMTOptions(
            method=self.options.method,
            samples=self.options.samples,
            seed=self.options.seed,
            max_simulated_time=self.options.max_simulated_time,
            conf_int=self.options.conf_int,
            max_rel_err=self.options.max_rel_err,
            verbose=self.options.verbose
        )

        # Call the handler
        self._result = solver_jmt(self._sn, handler_options)

        return self

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics as a DataFrame.

        Returns:
            DataFrame with columns: Station, Class, QLen, Util, RespT, Tput, ArvR
        """
        if self._result is None:
            self.runAnalyzer()

        M = self._sn.nstations
        K = self._sn.nclasses

        nodenames = self._sn.nodenames if self._sn.nodenames else [f'Station{i}' for i in range(M)]
        classnames = self._sn.classnames if self._sn.classnames else [f'Class{r}' for r in range(K)]

        # Get station names and identify source stations
        station_names = []
        source_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None

        for i in range(M):
            node_idx = int(self._sn.stationToNode[i]) if self._sn.stationToNode is not None else i
            if node_idx < len(nodenames):
                station_names.append(nodenames[node_idx])
            else:
                station_names.append(f'Station{i}')

            # Check if source station
            if nodetype is not None and node_idx < len(nodetype):
                if int(nodetype[node_idx]) == 0:  # SOURCE = 0
                    source_stations.add(i)

        # Get arrival rates from rates matrix for source stations
        rates = np.asarray(self._sn.rates) if hasattr(self._sn, 'rates') and self._sn.rates is not None else None

        rows = []
        for i in range(M):
            for r in range(K):
                is_source = i in source_stations

                # Get values with NaN handling
                qlen = self._result.Q[i, r] if self._result.Q is not None else np.nan
                util = self._result.U[i, r] if self._result.U is not None else np.nan
                respt = self._result.R[i, r] if self._result.R is not None else np.nan
                arvr = self._result.A[i, r] if self._result.A is not None else np.nan
                tput = self._result.T[i, r] if self._result.T is not None else np.nan

                # For source stations, replace NaN with 0 and set Tput to arrival rate
                if is_source:
                    qlen = 0.0 if np.isnan(qlen) else qlen
                    util = 0.0 if np.isnan(util) else util
                    respt = 0.0 if np.isnan(respt) else respt
                    arvr = 0.0  # Source has no arrivals to itself

                    # Set Tput from arrival rate
                    if np.isnan(tput) and rates is not None:
                        stationToNode = np.asarray(self._sn.stationToNode).flatten()
                        node_idx = int(stationToNode[i])
                        if node_idx < rates.shape[0] and r < rates.shape[1]:
                            tput = rates[node_idx, r]

                # For non-source stations, use throughput as arrival rate if ArvR is 0 or NaN
                if not is_source and (np.isnan(arvr) or arvr == 0.0):
                    arvr = tput if not np.isnan(tput) else 0.0

                rows.append({
                    'Station': station_names[i],
                    'JobClass': classnames[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': respt,  # Residence time equals response time for single visit
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not getattr(self, '_table_silent', False):
            print(df.to_string(index=False))

        return df

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.Q if self._result.Q is not None else np.array([])

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.U if self._result.U is not None else np.array([])

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.R if self._result.R is not None else np.array([])

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.T if self._result.T is not None else np.array([])

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K matrix)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.A if self._result.A is not None else np.array([])

    def getAvgChainTable(self) -> pd.DataFrame:
        """
        Get average performance metrics aggregated by chain.

        Returns:
            DataFrame with columns: Chain, QLen, Util, RespT, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        # Get chain information from model structure
        nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else self._sn.nclasses
        inchain = self._sn.inchain if hasattr(self._sn, 'inchain') else None

        rows = []
        for c in range(nchains):
            chain_name = f'Chain{c+1}'

            # Get classes in this chain
            if inchain is not None and c in inchain:
                chain_classes = inchain[c].flatten().astype(int)
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
                    if k < self._result.Q.shape[1]:
                        total_qlen += self._result.Q[i, k] if not np.isnan(self._result.Q[i, k]) else 0.0
                        total_util += self._result.U[i, k] if not np.isnan(self._result.U[i, k]) else 0.0
                        total_respt += self._result.R[i, k] if not np.isnan(self._result.R[i, k]) else 0.0
                        total_tput = max(total_tput, self._result.T[i, k] if not np.isnan(self._result.T[i, k]) else 0.0)

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

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times (1 x K)."""
        if self._result is None:
            self.runAnalyzer()
        # Sum response times across all stations for each class
        if self._result.R is not None:
            return np.nansum(self._result.R, axis=0, keepdims=True)
        return np.array([[]])

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (1 x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.X if self._result.X is not None else np.array([[]])

    def getProbSysAggr(self) -> Tuple[float, float]:
        """Get system state probability (simulation estimate).

        Note: JMT simulation does not directly compute state probabilities.
        Returns (0.0, 1.0) as a placeholder indicating total probability is 1.
        """
        if self._result is None:
            self.runAnalyzer()
        # For simulation, we return placeholder values
        return 0.0, 1.0

    def getProbAggr(self, station: int) -> np.ndarray:
        """Get aggregated state probabilities at station.

        Note: JMT simulation does not directly compute state probabilities.
        Returns empty array as placeholder.
        """
        if self._result is None:
            self.runAnalyzer()
        # Not supported in simulation - return empty
        return np.array([])

    def getRuntime(self) -> float:
        """Get solver runtime in seconds."""
        if self._result is None:
            return 0.0
        return self._result.runtime

    def getMethod(self) -> str:
        """Get the method used."""
        if self._result is None:
            return self.method
        return self._result.method

    def listValidMethods(self) -> List[str]:
        """List valid methods for this solver."""
        return [
            'default', 'jsim',
            'jmva', 'jmva.mva', 'jmva.amva', 'jmva.recal',
            'jmva.comom', 'jmva.chow', 'jmva.bs', 'jmva.aql',
            'jmva.lin', 'jmva.dmlin'
        ]

    @staticmethod
    def isAvailable() -> bool:
        """Check if JMT solver is available."""
        return is_jmt_available()

    @staticmethod
    def getFeatureSet() -> Set[str]:
        """Get the set of features supported by this solver."""
        return {
            'Sink', 'Source', 'Router', 'ClassSwitch',
            'Delay', 'DelayStation', 'Queue',
            'Fork', 'Join', 'Forker', 'Joiner', 'Logger',
            'Coxian', 'Cox2', 'APH', 'Erlang', 'Exp', 'HyperExp',
            'Det', 'Gamma', 'Lognormal', 'MAP', 'MMPP2',
            'Normal', 'PH', 'Pareto', 'Weibull', 'Replayer', 'Uniform',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer',
            'Buffer', 'Dispatcher', 'Server', 'JobSink', 'RandomSource',
            'ServiceTunnel', 'LogTunnel', 'Linkage',
            'Enabling', 'Timing', 'Firing', 'Storage', 'Place', 'Transition',
            'SchedStrategy_INF', 'SchedStrategy_PS', 'SchedStrategy_DPS',
            'SchedStrategy_FCFS', 'SchedStrategy_GPS', 'SchedStrategy_SIRO',
            'SchedStrategy_HOL', 'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR',
            'SchedStrategy_SEPT', 'SchedStrategy_SRPT', 'SchedStrategy_LEPT',
            'SchedStrategy_SJF', 'SchedStrategy_LJF', 'SchedStrategy_LPS',
            'SchedStrategy_POLLING', 'SchedStrategy_EXT',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'RoutingStrategy_RROBIN', 'RoutingStrategy_WRROBIN',
            'RoutingStrategy_KCHOICES',
            'ClosedClass', 'SelfLoopingClass', 'OpenClass',
            'Cache', 'CacheClassSwitcher',
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
            'ReplacementStrategy_SFIFO', 'ReplacementStrategy_LRU',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if this solver supports the given model."""
        # JMT supports most models - for now return True
        # A more complete implementation would check model features
        return True

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'jsim',
            'samples': 10000,
            'seed': 23000,
            'verbose': False,
            'keep': False,
            'conf_int': 0.99,
            'max_rel_err': 0.03,
        })

    def getCdfRespT(self, R=None):
        """Get response time CDF using exponential approximation.

        Uses mean response times to approximate CDF via exponential distribution.

        Args:
            R: Optional response time matrix (uses result if None)

        Returns:
            List of lists where RD[station][class] is a 2D array with columns [cdf, time]
        """
        if self._result is None:
            self.runAnalyzer()

        if R is None:
            R = self._result.R

        M = self._sn.nstations
        K = self._sn.nclasses

        # Initialize nested list structure: RD[station][class]
        RD = []
        for i in range(M):
            station_data = []
            for r in range(K):
                if R is not None and i < R.shape[0] and r < R.shape[1]:
                    mean_resp_t = R[i, r]
                    if mean_resp_t > 0 and not np.isnan(mean_resp_t):
                        # Exponential approximation: F(t) = 1 - exp(-t/mean)
                        lambda_rate = 1.0 / mean_resp_t
                        quantiles = np.linspace(0.001, 0.999, 100)
                        times = -np.log(1 - quantiles) / lambda_rate
                        cdf_vals = 1 - np.exp(-lambda_rate * times)

                        # Return as 2D array with columns [cdf, time]
                        cdf_data = np.column_stack([cdf_vals, times])
                        station_data.append(cdf_data)
                    else:
                        station_data.append(None)
                else:
                    station_data.append(None)
            RD.append(station_data)

        return RD

    def getPerctRespT(self, percentiles=None):
        """Get response time percentiles.

        Args:
            percentiles: Array of percentile values (default: [90, 95, 99])

        Returns:
            Tuple of (PercRT, PercTable) where PercRT is list of dicts
            and PercTable is a pandas DataFrame
        """
        import pandas as pd

        if percentiles is None:
            percentiles = np.array([90, 95, 99])
        else:
            percentiles = np.asarray(percentiles)

        if self._result is None:
            self.runAnalyzer()

        R = self._result.R
        M = self._sn.nstations
        K = self._sn.nclasses

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]
        percentiles_normalized = percentiles / 100.0

        station_names = self._sn.nodenames if self._sn.nodenames else [f'Station{i}' for i in range(M)]
        class_names = self._sn.classnames if self._sn.classnames else [f'Class{r}' for r in range(K)]

        for i in range(M):
            for r in range(K):
                if R is not None and i < R.shape[0] and r < R.shape[1]:
                    mean_resp_t = R[i, r]
                    if mean_resp_t > 0 and not np.isnan(mean_resp_t):
                        lambda_rate = 1.0 / mean_resp_t
                        perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                        PercRT.append({
                            'station': i + 1,
                            'class': r + 1,
                            'percentiles': percentiles.tolist(),
                            'values': perc_values.tolist(),
                        })

                        node_idx = int(self._sn.stationToNode[i]) if self._sn.stationToNode is not None else i
                        station_name = station_names[node_idx] if node_idx < len(station_names) else f'Station{i}'

                        row_data = {
                            'Station': station_name,
                            'Class': class_names[r] if r < len(class_names) else f'Class{r}',
                        }
                        for perc_col, perc_val in zip(perc_col_names, perc_values):
                            row_data[perc_col] = perc_val
                        rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    def __repr__(self) -> str:
        return f"SolverJMT(method='{self.method}', samples={self.options.samples})"

    # Snake case aliases for MATLAB compatibility
    avg_table = getAvgTable
    get_avg_table = getAvgTable
    avg_qlen = getAvgQLen
    prob_sys_aggr = getProbSysAggr
    prob_aggr = getProbAggr
    avg_util = getAvgUtil
    avg_resp_t = getAvgRespT
    avg_tput = getAvgTput
    avg_arv_r = getAvgArvR
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    run_analyzer = runAnalyzer
    get_runtime = getRuntime
    get_method = getMethod
    list_valid_methods = listValidMethods
    is_available = isAvailable
    get_feature_set = getFeatureSet
    default_options = defaultOptions
    cdf_resp_t = getCdfRespT
    cdf_respt = getCdfRespT
    get_tran_cdf_respt = getCdfRespT
    getTranCdfRespT = getCdfRespT
    perct_resp_t = getPerctRespT
    perct_respt = getPerctRespT
    getAvgNodeTable = getAvgTable
    avg_node_table = getAvgTable
    get_avg_node_table = getAvgTable
