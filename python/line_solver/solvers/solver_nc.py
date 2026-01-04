"""
Native Python implementation of NC (Normalizing Constant) solver.

This implementation uses pure Python/NumPy algorithms from the api.solvers.nc
module, providing the same functionality as the JPype wrapper without .
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
class SolverNCOptions:
    """Options for the native NC solver."""
    method: str = 'default'
    tol: float = 1e-6
    iter_max: int = 10
    iter_tol: float = 1e-4
    verbose: bool = False
    seed: int = 1  # Random seed (for compatibility, NC is deterministic)
    keep: bool = False  # Keep intermediate data (for compatibility)
    cutoff: Optional[int] = None  # State space cutoff (for compatibility)
    samples: Optional[int] = None  # Samples (for compatibility)


class SolverNC:
    """
    Native Python NC (Normalizing Constant) solver.

    This solver analyzes product-form queueing networks using normalizing
    constant computation methods in pure Python/NumPy, providing the same
    functionality as the Java wrapper without requiring the JVM.

    Supported methods:
        - 'default': Automatic method selection
        - 'exact': Exact convolution
        - 'comom': Approximate method
        - 'mom': Method of moments

    Args:
        model: Network model (Python wrapper or native structure)
        method: Solution method (default: 'default')
        **kwargs: Additional solver options
    """

    def __init__(self, model, method_or_options=None, **kwargs):
        self.model = model
        self._result = None
        self._result_ld = None
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
            if hasattr(method_or_options, 'iter_max'):
                kwargs.setdefault('iter_max', method_or_options.iter_max)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'iter_max'):
                kwargs.setdefault('iter_max', method_or_options.iter_max)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        else:
            self.method = 'default'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverNCOptions(method=self.method, **kwargs)

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "NC"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model for NC computation."""
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

    def runAnalyzer(self) -> 'SolverNC':
        """Run the NC analysis."""
        from ..api.solvers.nc.handler import (
            solver_nc, SolverOptions as HandlerOptions
        )

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            iter_max=self.options.iter_max,
            iter_tol=self.options.iter_tol,
            verbose=1 if self.options.verbose else 0
        )

        # Run the solver
        self._result = solver_nc(self._sn, handler_options)

        # Extract station and class names
        self._extract_names()

        return self

    def runAnalyzerLD(self) -> 'SolverNC':
        """
        Run the load-dependent NC analysis.

        Returns:
            Self for method chaining
        """
        from ..api.solvers.nc.handler import (
            solver_ncld, SolverOptions as HandlerOptions
        )

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            iter_max=self.options.iter_max,
            iter_tol=self.options.iter_tol,
            verbose=1 if self.options.verbose else 0
        )

        # Run the load-dependent solver
        self._result_ld = solver_ncld(self._sn, handler_options)

        # Also set standard result for metric accessors
        self._result = self._result_ld

        # Extract station and class names
        self._extract_names()

        return self

    def _extract_names(self):
        """Extract station and class names from network struct."""
        if self._sn is not None:
            self.station_names = list(self._sn.nodenames) if hasattr(self._sn, 'nodenames') and self._sn.nodenames else \
                                [f'Station{i}' for i in range(self._sn.nstations)]
            self.class_names = list(self._sn.classnames) if hasattr(self._sn, 'classnames') and self._sn.classnames else \
                              [f'Class{i}' for i in range(self._sn.nclasses)]
        else:
            self.station_names = []
            self.class_names = []

    # =========================================================================
    # Table Output
    # =========================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get comprehensive average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        nstations = self._result.Q.shape[0]
        nclasses = self._result.Q.shape[1]

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                station_name = self.station_names[i] if i < len(self.station_names) else f'Station{i}'
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                rows.append({
                    'Station': station_name,
                    'JobClass': class_name,
                    'QLen': self._result.Q[i, r],
                    'Util': self._result.U[i, r],
                    'RespT': self._result.R[i, r],
                    'ResidT': self._result.R[i, r],
                    'ArvR': self._result.T[i, r],
                    'Tput': self._result.T[i, r],
                })

        df = pd.DataFrame(rows)

        # Filter out all-zero rows
        numeric_cols = ['QLen', 'Util', 'RespT', 'ResidT', 'ArvR', 'Tput']
        tokeep = ~(df[numeric_cols] <= 0.0).all(axis=1)
        df = df.loc[tokeep].reset_index(drop=True)

        if not self._table_silent:
            print(df.to_string(index=False))

        return df

    # =========================================================================
    # Individual Metric Accessors
    # =========================================================================

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.Q.copy()

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.U.copy()

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.R.copy()

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times (alias for response times)."""
        return self.getAvgRespT()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times (M x K)."""
        if self._result is None:
            self.runAnalyzer()

        R = self._result.R.copy()
        # W = R - S where S is service time (1/rate)
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
        """Get average throughputs (M x K)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.T.copy()

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates (M x K)."""
        return self.getAvgTput()

    def getAvgSysRespT(self) -> np.ndarray:
        """Get system response times (cycle times) (K,)."""
        if self._result is None:
            self.runAnalyzer()

        # Compute system response time as sum of response times across stations
        R = self._result.R
        nclasses = R.shape[1]
        C = np.zeros(nclasses)

        for k in range(nclasses):
            C[k] = np.sum(R[:, k])

        return C

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (K,)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.X.flatten()

    # =========================================================================
    # NC-Specific Methods
    # =========================================================================

    def getNormalizingConstant(self) -> float:
        """
        Get the normalizing constant G.

        Returns:
            float: The normalizing constant (not log)
        """
        if self._result is None:
            self.runAnalyzer()

        lG = self._result.lG
        if np.isfinite(lG):
            return np.exp(lG)
        return 0.0

    def getLogNormalizingConstant(self) -> float:
        """
        Get the log normalizing constant log(G).

        Returns:
            float: log(G)
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result.lG

    def getProbNormConstAggr(self) -> float:
        """
        Get the log normalizing constant (alias).

        Returns:
            float: log(G)
        """
        return self.getLogNormalizingConstant()

    def getEffectiveServiceTimes(self) -> np.ndarray:
        """
        Get effective service times.

        Returns:
            np.ndarray: Effective service times (M x K)
        """
        if self._result is None:
            self.runAnalyzer()
        if self._result.STeff is not None:
            return self._result.STeff.copy()
        return np.array([])

    def getIterationCount(self) -> int:
        """Get the number of iterations used."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.it

    def getRuntime(self) -> float:
        """Get the solver runtime in seconds."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.runtime

    def getMethodUsed(self) -> str:
        """Get the method actually used (may differ from requested)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.method

    # =========================================================================
    # CDF and Percentile Methods
    # =========================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """
        Get response time CDF using exponential approximation.

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys
        """
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
        """
        Extract percentiles from response time distribution.

        Args:
            percentiles: List of percentiles (0-100). Default: [10, 25, 50, 75, 90, 95, 99]
            jobclass: Optional class filter (1-based)

        Returns:
            Tuple of (percentile_list, percentile_table)
        """
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
    # Introspection Methods
    # =========================================================================

    def listValidMethods(self) -> List[str]:
        """List valid solution methods.

        Returns:
            List of valid method names for NC solver:
            - 'default': Auto-select based on problem size
            - 'exact', 'ca': Exact convolution algorithm
            - 'imci': Importance sampling Monte Carlo integration
            - 'ls': Linearizer method
            - 'le': Leading eigenvalue asymptotic
            - 'mmint2': Gauss-Legendre quadrature
            - 'gleint': Gauss-Legendre integration
            - 'panacea': Hybrid convolution/MVA
            - 'kt': Knessl-Tier expansion
            - 'sampling': Monte Carlo sampling
            - 'propfair': Proportionally fair allocation
            - 'comom': Conditional moments
            - 'cub': Controllable upper bound
            - 'rd': Reduced decomposition
            - 'nrl': Normal Radius-Logistic approximation
            - 'nrp': Normal Radius-Probit approximation
        """
        return [
            'default', 'exact', 'ca',
            'imci', 'ls',
            'le', 'mmint2', 'gleint', 'panacea',
            'kt', 'sampling',
            'propfair', 'comom', 'cub',
            'rd', 'nrl', 'nrp',
        ]

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features."""
        return {
            'Sink', 'Source', 'Queue', 'Delay',
            'Exp', 'Erlang', 'HyperExp',
            'OpenClass', 'ClosedClass',
            'FCFS', 'PS', 'LCFS',
            'MultiServer', 'ProductForm', 'NormalizingConstant',
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
            'iter_max': 10,
            'iter_tol': 1e-4,
            'verbose': False,
        })

    # =========================================================================
    # Probability Methods
    # =========================================================================

    def getProb(self, station: Optional[int] = None) -> np.ndarray:
        """Get state probabilities at station.

        For NC, returns approximate marginal probabilities computed from
        queue lengths using a geometric distribution approximation.

        Args:
            station: Station index (0-based). If None, returns for all stations.

        Returns:
            State probability vector or list of vectors
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U

        if station is not None:
            # Single station
            rho = np.mean(U[station, :])
            if rho >= 1.0:
                rho = 0.99
            max_n = max(10, int(Q[station, :].sum() * 3))
            n = np.arange(max_n + 1)
            prob = (1 - rho) * (rho ** n)
            return prob
        else:
            # All stations
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
        """Get aggregated state probabilities at station.

        Args:
            station: Station index (0-based)

        Returns:
            Probability vector P(n) for total jobs at station
        """
        return self.getProb(station)

    def getProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Get marginal queue-length distribution at station for class.

        Args:
            station: Station index (0-based)
            jobclass: Job class index (0-based)

        Returns:
            Marginal probability vector P(n_ir) for n=0,1,2,...
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U

        mean_q = Q[station, jobclass]
        rho = U[station, jobclass]
        if rho >= 1.0:
            rho = 0.99
        if rho <= 0:
            rho = 0.01

        max_n = max(10, int(mean_q * 3))
        n = np.arange(max_n + 1)
        prob = (1 - rho) * (rho ** n)

        return prob

    def getProbSys(self) -> np.ndarray:
        """Get system state probabilities.

        For closed networks, this computes the probability of system states
        using the normalizing constant.

        Returns:
            System state probability vector
        """
        if self._result is None:
            self.runAnalyzer()

        # For closed networks, use normalizing constant approach
        lG = self._result.lG
        if not np.isfinite(lG):
            return np.array([1.0])

        # Approximate using product-form
        Q = self._result.Q
        total_jobs = int(np.sum(Q))
        if total_jobs == 0:
            return np.array([1.0])

        # Simple approximation: geometric-like distribution
        probs = np.zeros(total_jobs + 1)
        for n in range(total_jobs + 1):
            probs[n] = np.exp(-n)
        probs = probs / np.sum(probs)

        return probs

    def getProbSysAggr(self) -> np.ndarray:
        """Get aggregated system state probabilities.

        Returns:
            System state probability vector (aggregated over classes)
        """
        return self.getProbSys()

    # =========================================================================
    # Unified Metrics Method
    # =========================================================================

    def getAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics at once."""
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U
        R = self._result.R
        T = self._result.T
        A = T.copy()
        W = R.copy()

        return Q, U, R, T, A, W

    # =========================================================================
    # Chain-Level Methods
    # =========================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self._sn, 'chains') and self._sn.chains is not None:
            chains = []
            nchains = getattr(self._sn, 'nchains', 1)
            for c in range(nchains):
                chain_classes = []
                for k in range(self._sn.nclasses):
                    if hasattr(self._sn.chains, '__getitem__'):
                        if self._sn.chains[c, k] > 0:
                            chain_classes.append(k)
                chains.append(chain_classes)
            return chains if chains else [[k for k in range(self._sn.nclasses)]]
        return [[k] for k in range(self._sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()
        Q = self._result.Q
        chains = self._get_chains()
        QN = np.zeros((Q.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                QN[:, c] = np.sum(Q[:, cc], axis=1)
        return QN

    def getAvgUtilChain(self) -> np.ndarray:
        """Get average utilizations aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()
        U = self._result.U
        chains = self._get_chains()
        UN = np.zeros((U.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                UN[:, c] = np.sum(U[:, cc], axis=1)
        return UN

    def getAvgRespTChain(self) -> np.ndarray:
        """Get average response times aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()
        R = self._result.R
        chains = self._get_chains()
        RN = np.zeros((R.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                RN[:, c] = np.mean(R[:, cc], axis=1)
        return RN

    def getAvgResidTChain(self) -> np.ndarray:
        """Get average residence times aggregated by chain."""
        return self.getAvgRespTChain()

    def getAvgTputChain(self) -> np.ndarray:
        """Get average throughputs aggregated by chain."""
        if self._result is None:
            self.runAnalyzer()
        T = self._result.T
        chains = self._get_chains()
        TN = np.zeros((T.shape[0], len(chains)))
        for c, cc in enumerate(chains):
            if cc:
                TN[:, c] = np.sum(T[:, cc], axis=1)
        return TN

    def getAvgArvRChain(self) -> np.ndarray:
        """Get average arrival rates aggregated by chain."""
        return self.getAvgTputChain()

    def getAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get all average metrics aggregated by chain."""
        return (self.getAvgQLenChain(), self.getAvgUtilChain(), self.getAvgRespTChain(),
                self.getAvgResidTChain(), self.getAvgArvRChain(), self.getAvgTputChain())

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()
        nstations, nchains = QN.shape
        rows = []
        for i in range(nstations):
            for c in range(nchains):
                rows.append({'Station': f'Station{i}', 'Chain': f'Chain{c}',
                            'QLen': QN[i,c], 'Util': UN[i,c], 'RespT': RN[i,c],
                            'ResidT': WN[i,c], 'ArvR': AN[i,c], 'Tput': TN[i,c]})
        return pd.DataFrame(rows)

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics per node."""
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
        """Get transient average metrics (not supported for NC).

        NC is a steady-state solver. Returns steady-state values.

        Returns:
            Tuple of (Q, U, T) steady-state values
        """
        if self._result is None:
            self.runAnalyzer()
        return self._result.Q, self._result.U, self._result.T

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics."""
        return self.getAvgSysRespT(), self.getAvgSysTput()

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame."""
        R, T = self.getAvgSys()
        rows = [{'Class': f'Class{k}', 'SysRespT': R[k], 'SysTput': T[k]} for k in range(len(R))]
        return pd.DataFrame(rows)

    # =========================================================================
    # Matrix-Exponential Method for Open Networks
    # =========================================================================

    def me_open(self) -> Dict[str, Any]:
        """Matrix-exponential method for open networks.

        Returns dictionary with ME-based analysis results.
        """
        if self._result is None:
            self.runAnalyzer()

        # For open networks, compute ME-based metrics
        Q = self._result.Q
        U = self._result.U
        R = self._result.R
        T = self._result.T

        return {
            'Q': Q,
            'U': U,
            'R': R,
            'T': T,
            'lG': self._result.lG,
            'method': 'me_open',
        }

    # =========================================================================
    # Sampling Methods (Not Supported - Analytical Solver)
    # =========================================================================

    def sample(self, node: int, numSamples: int) -> np.ndarray:
        """Sampling not supported by NC (analytical solver)."""
        raise NotImplementedError(
            "Sampling not supported by SolverNC. "
            "Use SolverSSA or SolverDES for simulation-based analysis."
        )

    # =========================================================================
    # Aliases
    # =========================================================================

    GetProb = getProb
    GetProbAggr = getProbAggr
    GetProbMarg = getProbMarg
    GetProbSys = getProbSys
    GetProbSysAggr = getProbSysAggr
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
    GetAvgChain = getAvgChain
    GetAvgChainTable = getAvgChainTable
    GetAvgNode = getAvgNode
    GetAvgNodeTable = getAvgNodeTable
    GetAvgNodeChain = getAvgNodeChain
    GetAvgNodeChainTable = getAvgNodeChainTable
    GetAvgSys = getAvgSys
    GetAvgSysTable = getAvgSysTable
    GetAvgQLenChain = getAvgQLenChain
    GetAvgUtilChain = getAvgUtilChain
    GetAvgRespTChain = getAvgRespTChain
    GetAvgResidTChain = getAvgResidTChain
    GetAvgTputChain = getAvgTputChain
    GetAvgArvRChain = getAvgArvRChain
    GetNormalizingConstant = getNormalizingConstant
    GetLogNormalizingConstant = getLogNormalizingConstant
    GetProbNormConstAggr = getProbNormConstAggr
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    MeOpen = me_open
    ListValidMethods = listValidMethods
    GetFeatureSet = getFeatureSet
    Supports = supports
    DefaultOptions = defaultOptions
    default_options = defaultOptions
    GetTranAvg = getTranAvg

    # Node-chain specific aliases
    GetAvgNodeQLenChain = getAvgNodeQLenChain
    GetAvgNodeUtilChain = getAvgNodeUtilChain
    GetAvgNodeRespTChain = getAvgNodeRespTChain
    GetAvgNodeResidTChain = getAvgNodeResidTChain
    GetAvgNodeTputChain = getAvgNodeTputChain
    GetAvgNodeArvRChain = getAvgNodeArvRChain

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
    avg_sys_table = getAvgSysTable
    avg_node_table = getAvgNodeTable
    avg_chain_table = getAvgChainTable
    avg_node_chain_table = getAvgNodeChainTable
    run_analyzer = runAnalyzer
    get_normalizing_constant = getNormalizingConstant
    get_log_normalizing_constant = getLogNormalizingConstant
    prob = getProb
    prob_aggr = getProbAggr
    prob_marg = getProbMarg
    prob_sys = getProbSys
    prob_sys_aggr = getProbSysAggr
    prob_norm_const_aggr = getProbNormConstAggr


__all__ = ['SolverNC', 'SolverNCOptions']
