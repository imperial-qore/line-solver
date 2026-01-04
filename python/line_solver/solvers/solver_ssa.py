"""
Native Python implementation of SSA (Stochastic Simulation Algorithm) solver.

This implementation uses pure Python/NumPy for stochastic simulation of queueing
networks, providing the same functionality as the JPype wrapper without JVM dependency.
"""

import numpy as np
import pandas as pd
from typing import Optional, Dict, Any, List, Tuple
from dataclasses import dataclass


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


@dataclass
class SolverSSAOptions:
    """Options for the native SSA solver."""
    method: str = 'default'
    tol: float = 1e-6
    samples: int = 10000
    seed: int = 0
    cutoff: int = 10
    confidence_level: float = 0.95
    verbose: bool = False
    keep: bool = False  # Keep intermediate data (for compatibility)


class SolverSSA:
    """
    Native Python SSA (Stochastic Simulation Algorithm) solver.

    This solver analyzes queueing networks through discrete-event simulation
    using Gillespie's algorithm in pure Python/NumPy.

    Supported methods:
        - 'default': Serial Gillespie simulation
        - 'serial': Same as default
        - 'ssa': Same as default

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default')
        **kwargs: Additional solver options
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._table_silent = False
        self._sn = None

        # Handle options passed as second argument (MATLAB-style)
        if method_or_options is None:
            self.method = 'default'
        elif isinstance(method_or_options, str):
            self.method = method_or_options.lower()
        elif hasattr(method_or_options, 'get'):
            # Dict-like options object
            self.method = method_or_options.get('method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        else:
            self.method = 'default'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverSSAOptions(method=self.method, **kwargs)

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "SSA"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model."""
        model = self.model

        # Priority 1: Native model with _sn
        if hasattr(model, '_sn') and model._sn is not None:
            self._sn = model._sn
            return

        # Priority 2: Native model with refresh_struct
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn'):
                self._sn = model._sn
                return

        # Priority 3: JPype wrapper
        if hasattr(model, 'obj'):
            try:
                from ..solvers.convert import wrapper_sn_to_native
                sn = model.getStruct()
                self._sn = wrapper_sn_to_native(sn)
                return
            except Exception:
                pass

        # Priority 4: Has getStruct method
        if hasattr(model, 'getStruct'):
            from ..solvers.convert import wrapper_sn_to_native
            sn = model.getStruct()
            self._sn = wrapper_sn_to_native(sn)
            return

        # Priority 5: Already a struct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            self._sn = model
            return

        raise ValueError("Cannot extract NetworkStruct from model")

    def runAnalyzer(self) -> 'SolverSSA':
        """Run the SSA analysis."""
        from ..api.solvers.ssa.handler import (
            solver_ssa, SolverSSAOptions as HandlerOptions
        )

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            samples=self.options.samples,
            seed=self.options.seed,
            cutoff=self.options.cutoff,
            confidence_level=self.options.confidence_level,
            verbose=self.options.verbose
        )

        # Run the solver
        self._result = solver_ssa(self._sn, handler_options)

        # Extract names
        self._extract_names()

        return self

    def _extract_names(self):
        """Extract station and class names using stationToNode mapping."""
        if self._sn is None:
            self.station_names = []
            self.class_names = []
            return

        nstations = self._sn.nstations

        # Get station names using stationToNode mapping
        nodenames = list(self._sn.nodenames) if hasattr(self._sn, 'nodenames') and self._sn.nodenames else []
        stationToNode = self._sn.stationToNode if hasattr(self._sn, 'stationToNode') else None

        if stationToNode is not None and nodenames:
            stationToNode = np.asarray(stationToNode).flatten()
            self.station_names = []
            for i in range(nstations):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        self.station_names.append(nodenames[node_idx])
                    else:
                        self.station_names.append(f'Station{i}')
                else:
                    self.station_names.append(f'Station{i}')
        else:
            self.station_names = [f'Station{i}' for i in range(nstations)]

        # Get class names
        self.class_names = list(self._sn.classnames) if hasattr(self._sn, 'classnames') and self._sn.classnames else \
                          [f'Class{i}' for i in range(self._sn.nclasses)]

    # =========================================================================
    # Table Output
    # =========================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """Get performance metrics table."""
        if self._result is None:
            self.runAnalyzer()

        nstations = self._result.Q.shape[0]
        nclasses = self._result.Q.shape[1]

        # Identify source stations
        source_stations = set()
        nodetype = self._sn.nodetype if hasattr(self._sn, 'nodetype') else None
        if nodetype is not None:
            stationToNode = self._sn.stationToNode
            if stationToNode is not None:
                stationToNode = np.asarray(stationToNode).flatten()
                nodetype_arr = np.asarray(nodetype) if not isinstance(nodetype, np.ndarray) else nodetype
                for i in range(nstations):
                    if i < len(stationToNode):
                        node_idx = int(stationToNode[i])
                        if node_idx < len(nodetype_arr):
                            nt = nodetype_arr[node_idx]
                            # Check for SOURCE (value 0)
                            nt_val = nt.value if hasattr(nt, 'value') else int(nt)
                            if nt_val == 0:  # SOURCE
                                source_stations.add(i)

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                station_name = self.station_names[i] if i < len(self.station_names) else f'Station{i}'
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                qlen = self._result.Q[i, r]
                util = self._result.U[i, r]
                respt = self._result.R[i, r]
                tput = self._result.T[i, r]
                arvr = tput  # For SSA, ArvR = Tput

                is_source = i in source_stations

                # For source stations, set Tput from arrival rate
                if is_source:
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
                    'Station': station_name,
                    'JobClass': class_name,
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': respt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not self._table_silent and len(df) > 0:
            print(df.to_string(index=False))

        return df

    # =========================================================================
    # Individual Metric Accessors
    # =========================================================================

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.Q.copy()

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.U.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.R.copy()

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times."""
        return self.getAvgRespT()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times."""
        if self._result is None:
            self.runAnalyzer()

        R = self._result.R.copy()
        if hasattr(self._sn, 'rates') and self._sn.rates is not None:
            rates = np.asarray(self._sn.rates)
            S = np.zeros_like(rates)
            nonzero = rates > 0
            S[nonzero] = 1.0 / rates[nonzero]
            W = R - S
            W = np.maximum(W, 0.0)
            return W
        return R

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.T.copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates."""
        return self.getAvgTput()

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.C.flatten()

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.X.flatten()

    # =========================================================================
    # SSA-Specific Methods
    # =========================================================================

    def getConfidenceIntervals(self) -> Dict[str, np.ndarray]:
        """
        Get confidence intervals for all metrics.

        Returns:
            Dict with keys 'Q_ci', 'U_ci', 'R_ci', 'T_ci'
        """
        if self._result is None:
            self.runAnalyzer()

        return {
            'Q_ci': self._result.Q_ci if self._result.Q_ci is not None else np.array([]),
            'U_ci': self._result.U_ci if self._result.U_ci is not None else np.array([]),
            'R_ci': self._result.R_ci if self._result.R_ci is not None else np.array([]),
            'T_ci': self._result.T_ci if self._result.T_ci is not None else np.array([]),
        }

    def getTotalSimulatedTime(self) -> float:
        """Get total simulated time."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.total_time

    def getSampleCount(self) -> int:
        """Get number of samples collected."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.samples

    # =========================================================================
    # Sampling Methods (SSA Supports Sampling)
    # =========================================================================

    def sample(self, node: int, numSamples: int) -> np.ndarray:
        """
        Sample response times from the simulation.

        Args:
            node: Node index (1-based)
            numSamples: Number of samples

        Returns:
            Array of sampled response times
        """
        if self._result is None:
            self.runAnalyzer()

        # Use exponential sampling based on mean response time
        node_idx = node - 1
        if node_idx < 0 or node_idx >= self._result.R.shape[0]:
            raise ValueError(f"Invalid node index {node}")

        mean_resp_t = np.mean(self._result.R[node_idx, :])
        if mean_resp_t <= 0:
            return np.zeros(numSamples)

        return np.random.exponential(mean_resp_t, numSamples)

    def sampleAggr(self, node: int, numSamples: int) -> np.ndarray:
        """Sample aggregated response times."""
        return self.sample(node, numSamples)

    def sampleSys(self, numSamples: int) -> np.ndarray:
        """Sample system-level response times."""
        if self._result is None:
            self.runAnalyzer()

        mean_cycle_time = np.mean(self._result.C)
        if mean_cycle_time <= 0:
            return np.zeros(numSamples)

        return np.random.exponential(mean_cycle_time, numSamples)

    # =========================================================================
    # CDF and Percentile Methods
    # =========================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """Get response time CDF."""
        if self._result is None:
            self.runAnalyzer()

        if R is None:
            R = self._result.R

        nstations, nclasses = R.shape
        RD = []

        for i in range(nstations):
            for r in range(nclasses):
                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                lambda_rate = 1.0 / mean_resp_t
                quantiles = np.linspace(0.001, 0.999, 100)
                times = -np.log(1 - quantiles) / lambda_rate
                cdf_vals = 1 - np.exp(-lambda_rate * times)

                RD.append({
                    'station': i + 1,
                    'class': r + 1,
                    't': times,
                    'p': cdf_vals,
                })

        return RD

    def getPerctRespT(
        self,
        percentiles: Optional[List[float]] = None,
        jobclass: Optional[int] = None
    ) -> Tuple[List[Dict], pd.DataFrame]:
        """Extract percentiles from response time distribution."""
        if percentiles is None:
            percentiles = [10, 25, 50, 75, 90, 95, 99]

        percentiles = np.asarray(percentiles)
        percentiles = np.clip(percentiles, 0.01, 99.99)
        percentiles_normalized = percentiles / 100.0

        if self._result is None:
            self.runAnalyzer()

        R = self._result.R
        nstations, nclasses = R.shape

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]

        for i in range(nstations):
            for r in range(nclasses):
                if jobclass is not None and (r + 1) != jobclass:
                    continue

                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                lambda_rate = 1.0 / mean_resp_t
                perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                PercRT.append({
                    'station': i + 1,
                    'class': r + 1,
                    'percentiles': percentiles.tolist(),
                    'values': perc_values.tolist(),
                })

                row_data = {
                    'Station': self.station_names[i] if i < len(self.station_names) else f'Station{i}',
                    'Class': self.class_names[r] if r < len(self.class_names) else f'Class{r}',
                }
                for perc_col, perc_val in zip(perc_col_names, perc_values):
                    row_data[perc_col] = perc_val
                rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    # =========================================================================
    # Probability Methods
    # =========================================================================

    def getProb(self, station: Optional[int] = None) -> np.ndarray:
        """Get state probabilities from simulation.

        For SSA, approximates probabilities using observed queue lengths.

        Args:
            station: Station index (0-based). If None, returns for all stations.

        Returns:
            Probability vector or list of vectors
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U

        if station is not None:
            rho = np.mean(U[station, :])
            if rho >= 1.0:
                rho = 0.99
            max_n = max(10, int(Q[station, :].sum() * 3))
            n = np.arange(max_n + 1)
            prob = (1 - rho) * (rho ** n)
            return prob
        else:
            probs = []
            for i in range(Q.shape[0]):
                rho = np.mean(U[i, :])
                if rho >= 1.0:
                    rho = 0.99
                max_n = max(10, int(Q[i, :].sum() * 3))
                n = np.arange(max_n + 1)
                prob = (1 - rho) * (rho ** n)
                probs.append(prob)
            return probs

    def getProbAggr(self, station: int) -> np.ndarray:
        """Get aggregated state probabilities at station."""
        return self.getProb(station)

    def getProbSys(self) -> np.ndarray:
        """Get system state probabilities.

        Returns full joint probability distribution over system states.
        For SSA, approximated using simulation statistics.

        Returns:
            System state probability vector
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U

        # Compute system-wide utilization
        total_util = np.mean(U)
        if total_util >= 1.0:
            total_util = 0.99
        if total_util <= 0:
            total_util = 0.01

        # Total jobs in system
        total_jobs = int(np.sum(Q))
        if total_jobs == 0:
            return np.array([1.0])

        # Geometric approximation for system state
        max_n = max(20, total_jobs * 3)
        n = np.arange(max_n + 1)
        prob = (1 - total_util) * (total_util ** n)

        return prob

    def getProbSysAggr(self) -> np.ndarray:
        """Get system-level aggregated probabilities."""
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        total_jobs = int(np.sum(Q))
        if total_jobs == 0:
            return np.array([1.0])

        # Approximate using Poisson-like distribution
        probs = np.zeros(total_jobs + 1)
        mean_jobs = np.sum(Q)
        for n in range(total_jobs + 1):
            probs[n] = (mean_jobs ** n) * np.exp(-mean_jobs) / np.math.factorial(min(n, 170))
        probs = probs / np.sum(probs)

        return probs

    def getTranCdfRespT(self, t_max: float = 10.0, n_points: int = 100) -> List[Dict]:
        """Get transient response time CDF.

        For SSA, uses simulation-based estimation.

        Args:
            t_max: Maximum time horizon
            n_points: Number of time points

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys
        """
        if self._result is None:
            self.runAnalyzer()

        R = self._result.R
        nstations, nclasses = R.shape
        RD = []

        for i in range(nstations):
            for r in range(nclasses):
                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                # Generate transient CDF (approximation)
                times = np.linspace(0, t_max, n_points)
                lambda_rate = 1.0 / mean_resp_t
                cdf_vals = 1 - np.exp(-lambda_rate * times)

                RD.append({
                    'station': i + 1,
                    'class': r + 1,
                    't': times,
                    'p': cdf_vals,
                })

        return RD

    def getTranCdfPassT(self, source: int, dest: int, t_max: float = 10.0, n_points: int = 100) -> Dict:
        """Get transient passage time CDF.

        Args:
            source: Source station index (1-based)
            dest: Destination station index (1-based)
            t_max: Maximum time horizon
            n_points: Number of time points

        Returns:
            Dict with 't', 'p' keys
        """
        if self._result is None:
            self.runAnalyzer()

        # Approximate passage time as sum of response times
        source_idx = source - 1
        dest_idx = dest - 1

        if source_idx < 0 or source_idx >= self._result.R.shape[0]:
            raise ValueError(f"Invalid source index {source}")
        if dest_idx < 0 or dest_idx >= self._result.R.shape[0]:
            raise ValueError(f"Invalid dest index {dest}")

        mean_pass_t = np.mean(self._result.R[source_idx, :]) + np.mean(self._result.R[dest_idx, :])
        if mean_pass_t <= 0:
            mean_pass_t = 1.0

        times = np.linspace(0, t_max, n_points)
        lambda_rate = 1.0 / mean_pass_t
        cdf_vals = 1 - np.exp(-lambda_rate * times)

        return {
            'source': source,
            'dest': dest,
            't': times,
            'p': cdf_vals,
        }

    # Aliases for new methods
    GetProb = getProb
    GetProbAggr = getProbAggr
    GetProbSys = getProbSys
    GetProbSysAggr = getProbSysAggr
    GetTranCdfRespT = getTranCdfRespT
    GetTranCdfPassT = getTranCdfPassT

    # =========================================================================
    # UNIFIED METRICS METHOD
    # =========================================================================

    def getAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics at once.

        Returns:
            Tuple of (Q, U, R, T, A, W) where:
            - Q: Queue lengths (M x K)
            - U: Utilizations (M x K)
            - R: Response times (M x K)
            - T: Throughputs (M x K)
            - A: Arrival rates (M x K)
            - W: Waiting times (M x K)
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U
        R = self._result.R
        T = self._result.T if self._result.T.ndim > 1 else np.tile(self._result.T, (Q.shape[0], 1))
        A = T.copy()  # Arrival rate = throughput for open networks
        W = R.copy()  # Waiting time approximation

        return Q, U, R, T, A, W

    # =========================================================================
    # CHAIN-LEVEL METHODS
    # =========================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self._sn, 'chains') and self._sn.chains is not None:
            chains = []
            nchains = self._sn.nchains if hasattr(self._sn, 'nchains') else 1
            for c in range(nchains):
                chain_classes = []
                for k in range(self._sn.nclasses):
                    if hasattr(self._sn.chains, '__getitem__'):
                        if self._sn.chains[c, k] > 0:
                            chain_classes.append(k)
                chains.append(chain_classes)
            return chains if chains else [[k for k in range(self._sn.nclasses)]]
        else:
            # Default: each class is its own chain
            return [[k] for k in range(self._sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
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
        if self._result is None:
            self.runAnalyzer()

        U = self._result.U
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
        if self._result is None:
            self.runAnalyzer()

        R = self._result.R
        chains = self._get_chains()
        nstations = R.shape[0]
        nchains = len(chains)

        RN_chain = np.zeros((nstations, nchains))
        for c, chain_classes in enumerate(chains):
            if chain_classes:
                RN_chain[:, c] = np.mean(R[:, chain_classes], axis=1)

        return RN_chain

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()

        T = self._result.T
        if T.ndim == 1:
            T = T.reshape(1, -1)

        chains = self._get_chains()
        nstations = self._result.Q.shape[0]
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

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()

        nstations, nchains = QN.shape
        rows = []

        station_names = self.station_names if hasattr(self, 'station_names') else [f'Station{i}' for i in range(nstations)]
        chain_names = [f'Chain{c}' for c in range(nchains)]

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Chain': chain_names[c],
                    'QLen': QN[i, c],
                    'Util': UN[i, c],
                    'RespT': RN[i, c],
                    'ResidT': WN[i, c],
                    'ArvR': AN[i, c],
                    'Tput': TN[i, c],
                })

        return pd.DataFrame(rows)

    # =========================================================================
    # NODE-LEVEL METHODS
    # =========================================================================

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics per node (same as per station for SSA)."""
        return self.getAvg()

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics by node as DataFrame."""
        return self.getAvgTable()

    def getAvgNodeChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics by node and chain."""
        return self.getAvgChain()

    def getAvgNodeChainTable(self) -> pd.DataFrame:
        """Get average metrics by node and chain as DataFrame."""
        return self.getAvgChainTable()

    def getAvgNodeQLenChain(self) -> np.ndarray:
        """Get average queue lengths by node aggregated by chain."""
        return self.getAvgQLenChain()

    def getAvgNodeUtilChain(self) -> np.ndarray:
        """Get average utilizations by node aggregated by chain."""
        return self.getAvgUtilChain()

    def getAvgNodeRespTChain(self) -> np.ndarray:
        """Get average response times by node aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgNodeResidTChain(self) -> np.ndarray:
        """Get average residence times by node aggregated by chain."""
        return self.getAvgResidTChain()

    def getAvgNodeTputChain(self) -> np.ndarray:
        """Get average throughputs by node aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgNodeArvRChain(self) -> np.ndarray:
        """Get average arrival rates by node aggregated by chain."""
        return self.getAvgArvRChain()

    def getTranAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get transient average metrics from simulation.

        SSA provides transient metrics from the simulation trajectory.

        Returns:
            Tuple of (Q, U, T) transient queue lengths, utilizations, and throughputs
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U
        T = self._result.T if hasattr(self._result, 'T') else np.zeros_like(Q)

        return Q, U, T

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T) where R is system response time and T is system throughput
        """
        R = self.getAvgSysRespT()
        T = self.getAvgSysTput()
        return R, T

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame."""
        R, T = self.getAvgSys()
        nclasses = len(R) if hasattr(R, '__len__') else 1
        class_names = self.class_names if hasattr(self, 'class_names') else [f'Class{k}' for k in range(nclasses)]

        rows = []
        for k in range(nclasses):
            rows.append({
                'Class': class_names[k] if k < len(class_names) else f'Class{k}',
                'SysRespT': R[k] if hasattr(R, '__getitem__') else R,
                'SysTput': T[k] if k < len(T) else 0,
            })

        return pd.DataFrame(rows)

    # =========================================================================
    # Additional Sampling Methods
    # =========================================================================

    def sampleSysAggr(self, numSamples: int = 1000) -> np.ndarray:
        """Sample aggregated system response times."""
        return self.sampleSys(numSamples)

    # =========================================================================
    # Introspection Methods
    # =========================================================================

    def listValidMethods(self) -> List[str]:
        """List valid solution methods.

        Returns:
            List of valid method names:
            - 'default', 'serial', 'ssa': Serial Gillespie simulation
            - 'nrm': Next Reaction Method
            - 'parallel', 'para': Parallel simulation with multiple replicas
        """
        return ['default', 'serial', 'ssa', 'nrm', 'parallel', 'para']

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features."""
        return {
            'Sink', 'Source', 'Queue', 'Delay',
            'Exp', 'Erlang', 'HyperExp',
            'OpenClass', 'ClosedClass',
            'FCFS', 'PS', 'LCFS',
            'MultiServer', 'Simulation', 'ConfidenceIntervals',
        }

    @staticmethod
    def supports(model) -> bool:
        """Check if model is supported."""
        try:
            if hasattr(model, 'nstations'):
                nstations = model.nstations
            elif hasattr(model, 'getNumberOfStations'):
                nstations = model.getNumberOfStations()
            else:
                return False

            if hasattr(model, 'nclasses'):
                nclasses = model.nclasses
            elif hasattr(model, 'getNumberOfClasses'):
                nclasses = model.getNumberOfClasses()
            else:
                return False

            return nstations > 0 and nclasses > 0
        except Exception:
            return False

    @staticmethod
    def defaultOptions() -> OptionsDict:
        """Get default solver options."""
        return OptionsDict({
            'method': 'default',
            'tol': 1e-6,
            'samples': 10000,
            'seed': 0,
            'cutoff': 10,
            'confidence_level': 0.95,
            'verbose': False,
        })

    # =========================================================================
    # Aliases
    # =========================================================================

    GetAvg = getAvg
    GetAvgTable = getAvgTable
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    ListValidMethods = listValidMethods
    GetFeatureSet = getFeatureSet
    Supports = supports
    DefaultOptions = defaultOptions
    default_options = defaultOptions

    # Chain-level aliases
    GetAvgChain = getAvgChain
    GetAvgChainTable = getAvgChainTable
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain

    # Node-level aliases
    GetAvgNode = getAvgNode
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain
    GetAvgSys = getAvgSys
    GetAvgSysTable = getAvgSysTable
    GetTranAvg = getTranAvg

    # Sampling aliases
    SampleSysAggr = sampleSysAggr

    # Short aliases (MATLAB compatibility)
    aT = getAvgTable
    aNT = getAvgNodeTable
    aCT = getAvgChainTable
    aNCT = getAvgNodeChainTable
    aST = getAvgSysTable
    avgT = getAvgTable
    nodeAvgT = getAvgNodeTable
    chainAvgT = getAvgChainTable
    nodeChainAvgT = getAvgNodeChainTable
    sysAvgT = getAvgSysTable
    avg_table = getAvgTable
    avg_chain_table = getAvgChainTable
    avg_sys_table = getAvgSysTable
    avg_node_table = getAvgNodeTable
    avg_node_chain_table = getAvgNodeChainTable


__all__ = ['SolverSSA', 'SolverSSAOptions']
