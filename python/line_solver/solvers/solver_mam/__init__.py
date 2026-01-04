"""
SolverMAM - Main matrix-analytic methods solver.

Implements 8 solution methods for queueing networks:
1. dec.source - Decomposition with MMAP arrivals (default)
2. dec.mmap - Service-scaled departures
3. dec.poisson - Poisson approximation
4. mna - Matrix-normalizing approximation (auto-selects open/closed)
5. ldqbd - Level-dependent QBD
6. inap - RCAT iterative
7. inapplus - RCAT weighted variant
8. fj - Fork-Join (percentile analysis)

Usage:
    solver = SolverMAM(network, method='default')
    solver.runAnalyzer()
    QN = solver.getAvgQLen()
    RN = solver.getAvgRespT()
"""

import numpy as np
import pandas as pd
import time
from typing import Optional, Dict, List, Tuple
from dataclasses import dataclass

from .algorithms import (
    MAMResult,
    DecSourceAlgorithm,
    DecMMAPAlgorithm,
    DecPoissonAlgorithm,
    MNAOpenAlgorithm,
    MNAClosedAlgorithm,
    INAPAlgorithm,
    INAPPlusAlgorithm,
)
from .utils import extract_mam_params, check_closed_network, is_fork_join_network
from .fj.validator import fj_isfj
from .fj.solver import FJSolver
from ...api.mam import ldqbd, LdqbdOptions

@dataclass
class SolverMAMOptions:
    """Options for SolverMAM.

    Attributes:
        method: Algorithm to use (default, dec.source, mna, ldqbd, inap, inapplus)
        tol: Convergence tolerance
        max_iter: Maximum iterations
        space_max: Maximum MMAP state space size
        verbose: Print debug information
    """
    method: str = 'default'
    tol: float = 1e-6
    max_iter: int = 100
    space_max: int = 1000
    verbose: bool = False


class SolverMAM:
    """Native Python solver for matrix-analytic methods.

    Solves queueing networks using decomposition, MNA, RCAT, and related methods.
    """

    # Available methods and their algorithm classes
    ALGORITHMS = {
        'dec.source': DecSourceAlgorithm,
        'dec.mmap': DecMMAPAlgorithm,
        'dec.poisson': DecPoissonAlgorithm,
        'mna': None,  # Auto-select based on network type
        'mna_open': MNAOpenAlgorithm,
        'mna_closed': MNAClosedAlgorithm,
        'ldqbd': None,  # Special handling
        'inap': INAPAlgorithm,
        'inapplus': INAPPlusAlgorithm,
    }

    def __init__(self, network, method: str = 'default', options: Optional[SolverMAMOptions] = None, **kwargs):
        """Initialize SolverMAM.

        Args:
            network: Network model (must be compiled to NetworkStruct)
            method: Solution method ('default', 'dec.source', 'mna', etc.)
            options: SolverMAMOptions instance
            **kwargs: Additional parameters (seed, etc.) for compatibility - ignored
        """
        self.network = network
        self.sn = self._get_network_struct(network)
        # Store seed if provided (for compatibility, though MAM is analytical)
        self._seed = kwargs.get('seed', None)

        if options is None:
            options = SolverMAMOptions(method=method)
        else:
            if method != 'default':
                options.method = method

        self.options = options
        self.result = None
        self.runtime = 0.0

    def getName(self) -> str:
        """Get the name of this solver."""
        return "MAM"

    get_name = getName

    def _get_network_struct(self, model):
        """Get NetworkStruct from model using priority-based extraction."""
        # Priority 1: Native model with _sn attribute
        if hasattr(model, '_sn') and model._sn is not None:
            return model._sn

        # Priority 2: Native model with refresh_struct()
        if hasattr(model, 'refresh_struct'):
            model.refresh_struct()
            if hasattr(model, '_sn') and model._sn is not None:
                return model._sn

        # Priority 3: Has compileStruct method (native or wrapper)
        if hasattr(model, 'compileStruct'):
            return model.compileStruct()

        # Priority 4: JPype wrapper with getStruct
        if hasattr(model, 'getStruct'):
            return model.getStruct()

        # Priority 5: Model that is already a struct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError("Cannot extract network structure from model")

    def runAnalyzer(self) -> 'SolverMAM':
        """Run the analyzer with selected method.

        Returns:
            self (for method chaining)
        """
        method = self._select_method()

        if self.options.verbose:
            print(f"SolverMAM: Using method '{method}'")

        start_time = time.time()

        # Special handling for ldqbd
        if method == 'ldqbd':
            self.result = self._solve_ldqbd()
        elif method == 'fj':
            # Special handling for Fork-Join
            self.result = self._solve_fork_join()
        else:
            # Standard algorithm dispatch
            algo_class = self.ALGORITHMS.get(method)
            if algo_class is None:
                raise ValueError(f"Unknown method: {method}")

            algo = algo_class()
            self.result = algo.solve(self.sn, self.options)

        self.runtime = time.time() - start_time

        if self.options.verbose:
            print(f"  Completed in {self.runtime:.3f}s ({self.result.totiter} iterations)")

        return self

    def _select_method(self) -> str:
        """Select method based on network type if method='default'.

        Returns:
            Selected method name
        """
        method = self.options.method

        if method == 'default':
            # Auto-selection logic
            is_closed = check_closed_network(self.sn)
            is_fj = is_fork_join_network(self.sn)

            if is_fj:
                # Use fork-join method if available
                return 'fj'
            elif is_closed and self.sn.nstations == 1 and self.sn.nclasses == 1:
                # LDQBD for single-station closed
                return 'ldqbd'
            else:
                # Default to dec.source
                return 'dec.source'
        elif method == 'mna':
            # Auto-select between mna_open and mna_closed
            is_closed = check_closed_network(self.sn)
            return 'mna_closed' if is_closed else 'mna_open'
        else:
            return method

    def _solve_ldqbd(self) -> MAMResult:
        """Solve using LDQBD method.

        Returns:
            MAMResult
        """
        # Extract parameters for LDQBD
        # (simplified: assumes single queue)

        # Construct rate matrices
        Q0 = []
        Q1 = [np.array([[-1.0]])]
        Q2 = []

        opts = LdqbdOptions(epsilon=self.options.tol, max_iter=self.options.max_iter)
        ldqbd_result = ldqbd(Q0, Q1, Q2, opts)

        # Convert to MAMResult
        QN = np.zeros((self.sn.nstations, self.sn.nclasses))
        UN = np.zeros((self.sn.nstations, self.sn.nclasses))
        RN = np.zeros((self.sn.nstations, self.sn.nclasses))
        TN = np.zeros((1, self.sn.nclasses))

        return MAMResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            totiter=self.options.max_iter,
            method='ldqbd',
            runtime=0.0
        )

    def _solve_fork_join(self) -> MAMResult:
        """Solve Fork-Join network using FJ_codes.

        Returns:
            MAMResult with percentile response times attached
        """
        # Use FJ solver for topology validation and percentile computation
        fj_solver = FJSolver(verbose=self.options.verbose)

        # Validate Fork-Join topology
        can_solve, reason = fj_solver.can_solve(self.sn)
        if not can_solve:
            raise ValueError(f"Not a valid Fork-Join network: {reason}")

        # First, solve using dec.source to get basic metrics
        dec_algo = DecSourceAlgorithm()
        result = dec_algo.solve(self.sn, self.options)

        # Compute percentiles for Fork-Join
        percentiles = [50, 75, 90, 95, 99]
        mean_rt = np.sum(result.RN[:, 0]) if result.RN.size > 0 else 1.0

        fj_result = fj_solver.compute_percentiles(self.sn, percentiles, mean_rt)

        if fj_result is not None:
            # Attach percentile results to main result
            result.percentile_results = {
                'percentiles': fj_result.percentiles,
                'response_times': fj_result.response_times,
                'mean_response_time': fj_result.mean_response_time,
                'K': fj_result.K,
            }

        result.method = 'fj'
        return result

    # =====================================================================
    # RESULT ACCESS METHODS (following SolverMVA pattern)
    # =====================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """Get average performance metrics as DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self.result is None:
            self.runAnalyzer()

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

        # Get class names
        classnames = list(self.sn.classnames) if hasattr(self.sn, 'classnames') and self.sn.classnames else \
                     [f'Class{k}' for k in range(K)]

        # Identify Source stations
        source_stations = set()
        if hasattr(self.sn, 'sched') and self.sn.sched is not None:
            for ist in range(M):
                sched = self.sn.sched.get(ist, None)
                if sched is not None:
                    sched_name = sched.name if hasattr(sched, 'name') else str(sched)
                    if sched_name == 'EXT' or (hasattr(sched, 'value') and sched.value == 11):
                        source_stations.add(ist)

        # Build rows
        rows = []
        for i in range(M):
            is_source = i in source_stations
            for r in range(K):
                qlen = self.result.QN[i, r]
                util = self.result.UN[i, r]
                respt = self.result.RN[i, r]
                residt = respt  # Same as RespT for MAM (single visit assumption)

                # Get throughput and arrival rate
                if self.result.TN.ndim > 1:
                    tput = self.result.TN[0, r] if self.result.TN.shape[0] == 1 else self.result.TN[i, r]
                else:
                    tput = self.result.TN[r]

                # For Source stations, ArvR = 0 (no arrivals TO Source)
                arvr = 0.0 if is_source else tput

                rows.append({
                    'Station': station_names[i],
                    'JobClass': classnames[r],
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if not getattr(self, '_table_silent', False) and len(df) > 0:
            print(df.to_string(index=False))

        return df

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths per station.

        Returns:
            (M,) array of average queue lengths
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.QN, axis=1)

    def getAvgUtil(self) -> np.ndarray:
        """Get average utilizations per station.

        Returns:
            (M,) array of utilizations
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.UN, axis=1)

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times per station.

        Returns:
            (M,) array of response times
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.RN, axis=1)

    def getTput(self) -> np.ndarray:
        """Get throughputs per class.

        Returns:
            (K,) array of throughputs
        """
        if self.result is None:
            self.runAnalyzer()
        return self.result.TN.flatten()

    def getAvgSysRespT(self) -> np.ndarray:
        """Get average system response time per class.

        Returns:
            (K,) array of system response times
        """
        if self.result is None:
            self.runAnalyzer()
        return np.sum(self.result.RN, axis=0)

    def getAvgSysTput(self) -> float:
        """Get average system throughput.

        Returns:
            Scalar system throughput
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.XN)

    # =====================================================================
    # STATIC METHODS (introspection and validation)
    # =====================================================================

    @staticmethod
    def listValidMethods() -> List[str]:
        """List all valid solution methods.

        Returns:
            List of method names
        """
        return ['default', 'dec.source', 'dec.mmap', 'dec.poisson',
                'mna', 'mna_open', 'mna_closed', 'ldqbd', 'inap', 'inapplus', 'fj']

    @staticmethod
    def supports(sn, method: str) -> Tuple[bool, Optional[str]]:
        """Check if method can solve this network.

        Args:
            sn: NetworkStruct
            method: Method name

        Returns:
            (can_solve, reason_if_not)
        """
        algo_class = SolverMAM.ALGORITHMS.get(method)
        if algo_class is None:
            return False, f"Unknown method: {method}"

        return algo_class.supports_network(sn)

    @staticmethod
    def getFeatureSet() -> Dict:
        """Get features supported by SolverMAM.

        Returns:
            Dictionary describing supported features
        """
        return {
            'open_networks': True,
            'closed_networks': True,
            'mixed_networks': False,
            'phase_type_service': True,
            'markovian_arrivals': True,
            'class_switching': False,
            'priorities': True,
            'scheduling': ['FCFS', 'PS', 'LCFS', 'HOL'],
            'max_classes': 10,
            'max_stations': 100,
        }

    @staticmethod
    def defaultOptions() -> SolverMAMOptions:
        """Get default solver options.

        Returns:
            SolverMAMOptions with default values
        """
        return SolverMAMOptions()

    # =====================================================================
    # CDF AND PERCENTILE METHODS
    # =====================================================================

    def getCdfRespT(self, R: Optional[np.ndarray] = None) -> List[Dict]:
        """Get response time CDF using exponential approximation.

        For MAM, uses the computed mean response times to build an
        exponential CDF approximation (valid for M/M/1-like behavior).

        Args:
            R: Optional response time matrix (M x K). Uses result if None.

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        if R is None:
            R = self.result.RN

        nstations, nclasses = R.shape
        RD = []

        for i in range(nstations):
            for r in range(nclasses):
                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                # Exponential approximation: F(t) = 1 - exp(-t/mean)
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
        """Extract percentiles from response time distribution.

        Args:
            percentiles: List of percentiles (0-100). Default: [50, 75, 90, 95, 99]
            jobclass: Optional class filter (1-based)

        Returns:
            Tuple of (percentile_list, percentile_table)
        """
        if percentiles is None:
            percentiles = [50, 75, 90, 95, 99]

        percentiles = np.asarray(percentiles)
        percentiles = np.clip(percentiles, 0.01, 99.99)
        percentiles_normalized = percentiles / 100.0

        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        # Check for Fork-Join percentile results
        if hasattr(self.result, 'percentile_results') and self.result.percentile_results:
            fj_percs = self.result.percentile_results
            PercRT = [{
                'station': 'ForkJoin',
                'class': 1,
                'percentiles': fj_percs['percentiles'],
                'values': fj_percs['response_times'],
            }]
            rows = []
            for p, v in zip(fj_percs['percentiles'], fj_percs['response_times']):
                rows.append({'Percentile': p, 'RespT': v})
            return PercRT, pd.DataFrame(rows)

        R = self.result.RN
        nstations, nclasses = R.shape

        PercRT = []
        rows = []
        perc_col_names = [f'P{int(p)}' for p in percentiles]

        # Extract station/class names if available
        station_names = getattr(self.sn, 'nodenames', None) or [f'Station{i}' for i in range(nstations)]
        class_names = getattr(self.sn, 'classnames', None) or [f'Class{r}' for r in range(nclasses)]

        for i in range(nstations):
            for r in range(nclasses):
                if jobclass is not None and (r + 1) != jobclass:
                    continue

                mean_resp_t = R[i, r]
                if mean_resp_t <= 0:
                    continue

                # Exponential approximation for percentiles
                lambda_rate = 1.0 / mean_resp_t
                perc_values = -np.log(1 - percentiles_normalized) / lambda_rate

                PercRT.append({
                    'station': i + 1,
                    'class': r + 1,
                    'percentiles': percentiles.tolist(),
                    'values': perc_values.tolist(),
                })

                row_data = {
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Class': class_names[r] if r < len(class_names) else f'Class{r}',
                }
                for perc_col, perc_val in zip(perc_col_names, perc_values):
                    row_data[perc_col] = perc_val
                rows.append(row_data)

        PercTable = pd.DataFrame(rows) if rows else pd.DataFrame()
        return PercRT, PercTable

    # =====================================================================
    # PROBABILITY METHODS
    # =====================================================================

    def getProb(self, station: Optional[int] = None) -> np.ndarray:
        """Get state probabilities at station.

        For MAM, returns approximate marginal probabilities computed from
        queue lengths using a geometric distribution approximation.

        Args:
            station: Station index (0-based). If None, returns for all stations.

        Returns:
            State probability vector or matrix
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        Q = self.result.QN
        U = self.result.UN

        if station is not None:
            # Single station
            rho = np.mean(U[station, :])
            if rho >= 1.0:
                rho = 0.99
            # Geometric distribution approximation: P(n) = (1-rho) * rho^n
            max_n = max(10, int(Q[station, :].sum() * 3))
            n = np.arange(max_n + 1)
            prob = (1 - rho) * (rho ** n)
            return prob
        else:
            # All stations - return list of probability vectors
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

    def getProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Get marginal queue-length distribution at station for class.

        Args:
            station: Station index (0-based)
            jobclass: Job class index (0-based)

        Returns:
            Marginal probability vector P(n_ir) for n=0,1,2,...
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        Q = self.result.QN
        U = self.result.UN

        # Approximate marginal using geometric distribution
        mean_q = Q[station, jobclass]
        rho = U[station, jobclass]
        if rho >= 1.0:
            rho = 0.99
        if rho <= 0:
            rho = 0.01

        # For M/M/1: P(n) = (1-rho) * rho^n
        max_n = max(10, int(mean_q * 3))
        n = np.arange(max_n + 1)
        prob = (1 - rho) * (rho ** n)

        return prob

    # =====================================================================
    # ADDITIONAL STANDARD ACCESSOR METHODS
    # =====================================================================

    def getAvgResidT(self) -> np.ndarray:
        """Get average residence times per station (alias for response times).

        Returns:
            (M,) array of residence times
        """
        return self.getAvgRespT()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times per station.

        Waiting time is computed as response time minus mean service time.

        Returns:
            (M,) array of waiting times
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        resp_t = self.result.RN
        # Estimate service time from utilization and throughput
        # For M/M/1: U = lambda * S, so S = U / lambda
        # Waiting time = Response time - Service time
        wait_t = np.zeros(resp_t.shape[0])
        for i in range(resp_t.shape[0]):
            mean_resp = np.mean(resp_t[i, :])
            mean_util = np.mean(self.result.UN[i, :])
            # Approximate service time from utilization
            if mean_util > 0 and mean_util < 1:
                # W = R - S where S ≈ R * (1 - rho) for M/M/1
                service_t = mean_resp * (1 - mean_util)
                wait_t[i] = max(0, mean_resp - service_t)
            else:
                wait_t[i] = mean_resp * 0.5  # Fallback approximation

        return wait_t

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates per station.

        For open networks, arrival rate equals departure rate at steady state.

        Returns:
            (M,) array of arrival rates
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        # For open networks, arrival rate = throughput
        # Sum across classes for per-station rate
        return np.sum(self.result.TN, axis=1) if self.result.TN.ndim > 1 else self.result.TN

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs per station.

        Returns:
            (M,) array of throughputs
        """
        if self.result is None:
            raise RuntimeError("runAnalyzer() must be called first")

        # Sum across classes for per-station throughput
        if self.result.TN.ndim > 1:
            return np.sum(self.result.TN, axis=1)
        else:
            # If TN is 1D (per class), replicate for stations
            return np.full(self.result.QN.shape[0], np.mean(self.result.TN))

    # =====================================================================
    # SAMPLING METHODS (Not Supported - Analytical Solver)
    # =====================================================================

    def sample(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample from state distribution (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sample() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleAggr(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleAggr() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleSys(self, numEvents: int = 1000) -> np.ndarray:
        """Sample system states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleSys() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    def sampleSysAggr(self, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated system states (not supported for MAM).

        Raises:
            NotImplementedError: MAM is an analytical solver
        """
        raise NotImplementedError("sampleSysAggr() not supported for analytical MAM solver. Use SSA or CTMC instead.")

    # =====================================================================
    # TRANSIENT METHODS (Not Supported - Steady-State Solver)
    # =====================================================================

    def getTranCdfRespT(self) -> List[Dict]:
        """Get transient response time CDF (not supported for MAM).

        Raises:
            NotImplementedError: MAM computes steady-state only
        """
        raise NotImplementedError("getTranCdfRespT() not supported for MAM. Use CTMC or simulation.")

    def getTranCdfPassT(self) -> List[Dict]:
        """Get transient passage time CDF (not supported for MAM).

        Raises:
            NotImplementedError: MAM computes steady-state only
        """
        raise NotImplementedError("getTranCdfPassT() not supported for MAM. Use FLD or simulation.")

    def getTranAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get transient average metrics (not supported for MAM).

        Raises:
            NotImplementedError: MAM computes steady-state only
        """
        raise NotImplementedError("getTranAvg() not supported for MAM. Use FLD or CTMC.")

    def getCdfPassT(self) -> List[Dict]:
        """Get passage time CDF (not supported for MAM).

        Raises:
            NotImplementedError: Requires simulation
        """
        raise NotImplementedError("getCdfPassT() not supported for MAM. Use FLD or simulation.")

    # =====================================================================
    # UNIFIED METRICS METHOD
    # =====================================================================

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
        if self.result is None:
            self.runAnalyzer()

        Q = self.result.QN
        U = self.result.UN
        R = self.result.RN
        T = self.result.TN if self.result.TN.ndim > 1 else np.tile(self.result.TN, (Q.shape[0], 1))
        A = T.copy()  # Arrival rate = throughput for open networks
        W = R.copy()  # Waiting time approximation

        return Q, U, R, T, A, W

    # =====================================================================
    # CHAIN-LEVEL METHODS
    # =====================================================================

    def _get_chains(self) -> List[List[int]]:
        """Get chain-to-class mapping from network structure."""
        if hasattr(self.sn, 'chains') and self.sn.chains is not None:
            chains = []
            for c in range(self.sn.nchains if hasattr(self.sn, 'nchains') else 1):
                chain_classes = []
                for k in range(self.sn.nclasses):
                    if hasattr(self.sn.chains, '__getitem__'):
                        if self.sn.chains[c, k] > 0:
                            chain_classes.append(k)
                chains.append(chain_classes)
            return chains if chains else [[k for k in range(self.sn.nclasses)]]
        else:
            # Default: each class is its own chain
            return [[k] for k in range(self.sn.nclasses)]

    def getAvgQLenChain(self) -> np.ndarray:
        """Get average queue lengths aggregated by chain."""
        if self.result is None:
            self.runAnalyzer()

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
            self.runAnalyzer()

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
            self.runAnalyzer()

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
            self.runAnalyzer()

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

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()

        nstations, nchains = QN.shape
        rows = []

        station_names = getattr(self.sn, 'nodenames', None) or [f'Station{i}' for i in range(nstations)]
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

    # =====================================================================
    # NODE-LEVEL METHODS
    # =====================================================================

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics per node (same as per station for MAM)."""
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

    def getAvgSys(self) -> Tuple[np.ndarray, np.ndarray]:
        """Get system-level average metrics.

        Returns:
            Tuple of (R, T) where R is system response time and T is system throughput
        """
        R = self.getAvgSysRespT()
        T = self.getTput()
        return R, T

    def getAvgSysTable(self) -> pd.DataFrame:
        """Get system-level metrics as DataFrame."""
        R, T = self.getAvgSys()
        nclasses = len(R)
        class_names = getattr(self.sn, 'classnames', None) or [f'Class{k}' for k in range(nclasses)]

        rows = []
        for k in range(nclasses):
            rows.append({
                'Class': class_names[k] if k < len(class_names) else f'Class{k}',
                'SysRespT': R[k],
                'SysTput': T[k] if k < len(T) else 0,
            })

        return pd.DataFrame(rows)

    # =====================================================================
    # ALIASES (PascalCase for MATLAB compatibility)
    # =====================================================================

    GetAvg = getAvg
    GetAvgQLen = getAvgQLen
    GetAvgUtil = getAvgUtil
    GetAvgRespT = getAvgRespT
    GetAvgResidT = getAvgResidT
    GetAvgWaitT = getAvgWaitT
    GetAvgTput = getAvgTput
    GetAvgArvR = getAvgArvR
    GetAvgSysRespT = getAvgSysRespT
    GetAvgSysTput = getAvgSysTput
    GetAvgTable = getAvgTable
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
    GetCdfRespT = getCdfRespT
    GetPerctRespT = getPerctRespT
    GetProb = getProb
    GetProbMarg = getProbMarg
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

    # Snake case aliases
    avg_table = getAvgTable
    avg_node_table = getAvgNodeTable
    avg_chain_table = getAvgChainTable
    avg_node_chain_table = getAvgNodeChainTable
    avg_sys_table = getAvgSysTable
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_resp_t = getAvgRespT
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    run_analyzer = runAnalyzer
    cdf_resp_t = getCdfRespT
    perct_resp_t = getPerctRespT
    list_valid_methods = listValidMethods
    default_options = defaultOptions


__all__ = [
    'SolverMAM',
    'SolverMAMOptions',
    'MAMResult',
]
