"""
Native Python implementation of POSTERIOR solver.

This solver handles Prior distributions by expanding the model into a family
of networks, one for each alternative in the Prior. Results are aggregated
using prior-weighted expectations (Bayesian-style parameter uncertainty analysis).
"""

import numpy as np
import pandas as pd
import copy
from dataclasses import dataclass, field
from typing import Optional, Dict, Any, List, Tuple, Callable, Union


@dataclass
class PosteriorOptions:
    """Options for the POSTERIOR solver."""
    verbose: bool = False
    solver_type: str = 'mva'  # Underlying solver type: mva, nc, ctmc, ssa, mam, fld


@dataclass
class PosteriorResult:
    """Result from POSTERIOR solver."""
    QN: np.ndarray = field(default_factory=lambda: np.array([]))  # Queue lengths [M x K]
    UN: np.ndarray = field(default_factory=lambda: np.array([]))  # Utilizations [M x K]
    RN: np.ndarray = field(default_factory=lambda: np.array([]))  # Response times [M x K]
    TN: np.ndarray = field(default_factory=lambda: np.array([]))  # Throughputs [M x K]
    AN: np.ndarray = field(default_factory=lambda: np.array([]))  # Arrival rates [M x K]
    WN: np.ndarray = field(default_factory=lambda: np.array([]))  # Residence times [M x K]
    CN: np.ndarray = field(default_factory=lambda: np.array([]))  # System response times
    XN: np.ndarray = field(default_factory=lambda: np.array([]))  # System throughputs
    runtime: float = 0.0
    num_alternatives: int = 0


@dataclass
class PriorInfo:
    """Information about a detected Prior distribution."""
    node_idx: int
    class_idx: int
    prior_type: str  # 'service' or 'arrival'
    alternatives: List[Any]  # List of alternative distributions
    probabilities: List[float]  # Probability of each alternative


@dataclass
class EmpiricalCDF:
    """Empirical CDF representing a discrete posterior distribution."""
    values: np.ndarray
    probabilities: np.ndarray
    cdf: np.ndarray

    def __init__(self, values: List[float], probabilities: List[float]):
        """Create empirical CDF from values and probabilities."""
        n = len(values)
        indices = np.argsort(values)

        self.values = np.array([values[i] for i in indices])
        self.probabilities = np.array([probabilities[i] for i in indices])
        self.cdf = np.cumsum(self.probabilities)

    def eval_cdf(self, x: float) -> float:
        """Evaluate CDF at point x."""
        for i, v in enumerate(self.values):
            if x < v:
                return self.cdf[i - 1] if i > 0 else 0.0
        return 1.0

    def get_mean(self) -> float:
        """Return mean of the distribution."""
        return float(np.sum(self.values * self.probabilities))

    def get_percentile(self, p: float) -> float:
        """Return p-th percentile (p in [0, 1])."""
        idx = np.searchsorted(self.cdf, p)
        if idx >= len(self.values):
            idx = len(self.values) - 1
        return float(self.values[idx])


class SolverPosterior:
    """
    Native Python POSTERIOR solver for Bayesian-style parameter uncertainty analysis.

    This solver wraps another solver and handles Prior distributions by expanding
    the model into a family of networks, one for each alternative in the Prior.
    Results are aggregated using prior-weighted expectations.

    Example:
        >>> from line_solver import Network, Queue, Source, Sink
        >>> from line_solver.processes import Prior, Exp
        >>>
        >>> model = Network('M/M/1 with uncertainty')
        >>> source = Source(model, 'Source')
        >>> queue = Queue(model, 'Queue', SchedStrategy.FCFS)
        >>> sink = Sink(model, 'Sink')
        >>>
        >>> # Prior over service rate: 60% chance of rate 2.0, 40% chance of rate 1.0
        >>> prior = Prior([Exp(2.0), Exp(1.0)], [0.6, 0.4])
        >>> queue.setService(job_class, prior)
        >>>
        >>> solver = SolverPosterior(model, options=PosteriorOptions(solver_type='mva'))
        >>> result = solver.runAnalyzer()
    """

    def __init__(self, model, options: Optional[PosteriorOptions] = None,
                 solver_factory: Optional[Callable] = None):
        """
        Initialize the POSTERIOR solver.

        Args:
            model: Network model containing Prior distributions
            options: Optional PosteriorOptions configuration
            solver_factory: Optional callable that creates a solver for a given model
                           If not provided, uses MVA solver by default
        """
        self.model = model
        self.options = options or PosteriorOptions()
        self.solver_factory = solver_factory
        self._result: Optional[PosteriorResult] = None
        self._table_silent = False

        # Prior detection
        self.prior_info: Optional[PriorInfo] = None
        self.ensemble: List[Any] = []  # Alternative models
        self.solvers: List[Any] = []  # Solvers for each alternative
        self.alternative_results: List[Dict[str, np.ndarray]] = []

        # Detect Prior distributions
        self._detect_prior()

    def _detect_prior(self) -> None:
        """
        Detect Prior distributions in the model.

        Currently supports a single Prior per model. The Prior can be in:
        - Service distributions at queues
        - Arrival distributions at sources
        """
        # Check if model has getStruct or is a native struct
        sn = None
        if hasattr(self.model, 'getStruct'):
            sn = self.model.getStruct()
        elif hasattr(self.model, 'proc'):
            # Direct struct
            sn = self.model

        if sn is None:
            return

        # Check service distributions
        if hasattr(sn, 'proc') and sn.proc is not None:
            M = sn.proc.shape[0] if hasattr(sn.proc, 'shape') else len(sn.proc)
            K = sn.proc.shape[1] if hasattr(sn.proc, 'shape') and len(sn.proc.shape) > 1 else 1

            for i in range(M):
                for k in range(K):
                    try:
                        proc = sn.proc[i, k] if hasattr(sn.proc, '__getitem__') else None
                        if proc is not None and self._is_prior(proc):
                            alternatives, probs = self._extract_prior_info(proc)
                            if alternatives:
                                self.prior_info = PriorInfo(
                                    node_idx=i,
                                    class_idx=k,
                                    prior_type='service',
                                    alternatives=alternatives,
                                    probabilities=probs
                                )
                                return
                    except (IndexError, TypeError):
                        continue

        # Check arrival distributions
        if hasattr(sn, 'arvproc') and sn.arvproc is not None:
            M = sn.arvproc.shape[0] if hasattr(sn.arvproc, 'shape') else len(sn.arvproc)
            K = sn.arvproc.shape[1] if hasattr(sn.arvproc, 'shape') and len(sn.arvproc.shape) > 1 else 1

            for i in range(M):
                for k in range(K):
                    try:
                        proc = sn.arvproc[i, k] if hasattr(sn.arvproc, '__getitem__') else None
                        if proc is not None and self._is_prior(proc):
                            alternatives, probs = self._extract_prior_info(proc)
                            if alternatives:
                                self.prior_info = PriorInfo(
                                    node_idx=i,
                                    class_idx=k,
                                    prior_type='arrival',
                                    alternatives=alternatives,
                                    probabilities=probs
                                )
                                return
                    except (IndexError, TypeError):
                        continue

    def _is_prior(self, dist) -> bool:
        """Check if a distribution is a Prior."""
        if dist is None:
            return False
        dist_type = type(dist).__name__
        return 'Prior' in dist_type or hasattr(dist, 'alternatives')

    def _extract_prior_info(self, prior) -> Tuple[List[Any], List[float]]:
        """Extract alternatives and probabilities from a Prior distribution."""
        alternatives = []
        probabilities = []

        if hasattr(prior, 'alternatives'):
            alternatives = list(prior.alternatives)
        elif hasattr(prior, 'getAlternatives'):
            alternatives = list(prior.getAlternatives())

        if hasattr(prior, 'probabilities'):
            probabilities = list(prior.probabilities)
        elif hasattr(prior, 'getProbabilities'):
            probabilities = list(prior.getProbabilities())

        # Normalize probabilities if needed
        if probabilities:
            total = sum(probabilities)
            if total > 0:
                probabilities = [p / total for p in probabilities]
        elif alternatives:
            # Uniform probabilities if not specified
            n = len(alternatives)
            probabilities = [1.0 / n] * n

        return alternatives, probabilities

    def hasPriorDistribution(self) -> bool:
        """Check if the model has a Prior distribution."""
        return self.prior_info is not None

    def getNumAlternatives(self) -> int:
        """Return the number of alternatives in the Prior."""
        if self.prior_info is None:
            return 0
        return len(self.prior_info.alternatives)

    def _create_solver(self, model):
        """Create a solver for the given model."""
        if self.solver_factory is not None:
            return self.solver_factory(model)

        # Default: use MVA solver
        solver_type = self.options.solver_type.lower()

        try:
            if solver_type == 'mva':
                from .solver_mva import SolverMVA
                return SolverMVA(model)
            elif solver_type == 'nc':
                from .solver_nc import SolverNC
                return SolverNC(model)
            elif solver_type == 'ctmc':
                from .solver_ctmc import SolverCTMC
                return SolverCTMC(model)
            elif solver_type == 'ssa':
                from .solver_ssa import SolverSSA
                return SolverSSA(model)
            elif solver_type == 'mam':
                from .solver_mam import SolverMAM
                return SolverMAM(model)
            elif solver_type == 'fld':
                from .solver_fld import SolverFLD
                return SolverFLD(model)
            else:
                from .solver_mva import SolverMVA
                return SolverMVA(model)
        except ImportError as e:
            raise RuntimeError(f"Cannot import solver type '{solver_type}': {e}")

    def _create_alternative_model(self, alt_idx: int):
        """
        Create a copy of the model with the Prior replaced by the specified alternative.

        Args:
            alt_idx: Index of the alternative in the Prior

        Returns:
            Modified copy of the model
        """
        if self.prior_info is None:
            return copy.deepcopy(self.model)

        # Deep copy the model
        model_copy = copy.deepcopy(self.model)

        # Get the alternative distribution
        alt_dist = self.prior_info.alternatives[alt_idx]

        # Replace the Prior with the alternative
        node_idx = self.prior_info.node_idx
        class_idx = self.prior_info.class_idx

        # Get struct from copy
        sn = None
        if hasattr(model_copy, 'getStruct'):
            sn = model_copy.getStruct()
        elif hasattr(model_copy, 'proc'):
            sn = model_copy

        if sn is not None:
            if self.prior_info.prior_type == 'service':
                if hasattr(sn, 'proc') and sn.proc is not None:
                    try:
                        sn.proc[node_idx, class_idx] = alt_dist
                    except (IndexError, TypeError):
                        pass
            elif self.prior_info.prior_type == 'arrival':
                if hasattr(sn, 'arvproc') and sn.arvproc is not None:
                    try:
                        sn.arvproc[node_idx, class_idx] = alt_dist
                    except (IndexError, TypeError):
                        pass

        return model_copy

    def _init_ensemble(self) -> None:
        """Initialize ensemble of alternative models and solvers."""
        if self.prior_info is None:
            return

        n = len(self.prior_info.alternatives)
        self.ensemble = []
        self.solvers = []

        for i in range(n):
            # Create alternative model
            alt_model = self._create_alternative_model(i)
            self.ensemble.append(alt_model)

            # Create solver for this alternative
            solver = self._create_solver(alt_model)
            self.solvers.append(solver)

    def _solve_alternative(self, alt_idx: int) -> Dict[str, np.ndarray]:
        """
        Solve a single alternative and return results.

        Args:
            alt_idx: Index of the alternative

        Returns:
            Dictionary with QN, UN, RN, TN, AN, WN arrays
        """
        solver = self.solvers[alt_idx]

        # Run the solver
        if hasattr(solver, 'runAnalyzer'):
            solver.runAnalyzer()
        elif hasattr(solver, 'run'):
            solver.run()

        # Extract results
        result = {}

        if hasattr(solver, 'result') and solver.result is not None:
            r = solver.result
            result['QN'] = np.asarray(r.QN) if hasattr(r, 'QN') and r.QN is not None else np.array([])
            result['UN'] = np.asarray(r.UN) if hasattr(r, 'UN') and r.UN is not None else np.array([])
            result['RN'] = np.asarray(r.RN) if hasattr(r, 'RN') and r.RN is not None else np.array([])
            result['TN'] = np.asarray(r.TN) if hasattr(r, 'TN') and r.TN is not None else np.array([])
            result['AN'] = np.asarray(r.AN) if hasattr(r, 'AN') and r.AN is not None else np.array([])
            result['WN'] = np.asarray(r.WN) if hasattr(r, 'WN') and r.WN is not None else np.array([])
        elif hasattr(solver, '_result') and solver._result is not None:
            r = solver._result
            if isinstance(r, dict):
                result['QN'] = np.asarray(r.get('QN', []))
                result['UN'] = np.asarray(r.get('UN', []))
                result['RN'] = np.asarray(r.get('RN', []))
                result['TN'] = np.asarray(r.get('TN', []))
                result['AN'] = np.asarray(r.get('AN', []))
                result['WN'] = np.asarray(r.get('WN', []))
            else:
                result['QN'] = np.asarray(r.QN) if hasattr(r, 'QN') and r.QN is not None else np.array([])
                result['UN'] = np.asarray(r.UN) if hasattr(r, 'UN') and r.UN is not None else np.array([])
                result['RN'] = np.asarray(r.RN) if hasattr(r, 'RN') and r.RN is not None else np.array([])
                result['TN'] = np.asarray(r.TN) if hasattr(r, 'TN') and r.TN is not None else np.array([])
                result['AN'] = np.asarray(r.AN) if hasattr(r, 'AN') and r.AN is not None else np.array([])
                result['WN'] = np.asarray(r.WN) if hasattr(r, 'WN') and r.WN is not None else np.array([])

        return result

    def _aggregate_results(self) -> None:
        """Aggregate results from all alternatives using prior probabilities."""
        if not self.alternative_results or self.prior_info is None:
            return

        probs = self.prior_info.probabilities
        n = len(probs)

        # Get dimensions from first result
        first = self.alternative_results[0]
        shapes = {}
        for key in ['QN', 'UN', 'RN', 'TN', 'AN', 'WN']:
            if key in first and first[key].size > 0:
                shapes[key] = first[key].shape

        if not shapes:
            return

        # Initialize aggregated arrays
        aggregated = {}
        for key, shape in shapes.items():
            aggregated[key] = np.zeros(shape)

        # Aggregate with prior weights
        for i, alt_result in enumerate(self.alternative_results):
            p = probs[i]
            for key in aggregated:
                if key in alt_result and alt_result[key].size > 0:
                    if alt_result[key].shape == aggregated[key].shape:
                        aggregated[key] += p * alt_result[key]

        # Store result
        self._result = PosteriorResult(
            QN=aggregated.get('QN', np.array([])),
            UN=aggregated.get('UN', np.array([])),
            RN=aggregated.get('RN', np.array([])),
            TN=aggregated.get('TN', np.array([])),
            AN=aggregated.get('AN', np.array([])),
            WN=aggregated.get('WN', np.array([])),
            num_alternatives=n
        )

    def runAnalyzer(self) -> 'SolverPosterior':
        """
        Run the posterior analysis.

        Returns:
            self for method chaining
        """
        import time
        start_time = time.time()

        if self.prior_info is None:
            if self.options.verbose:
                print("No Prior distribution found in model")
            self._result = PosteriorResult()
            return self

        # Initialize ensemble
        self._init_ensemble()

        # Solve each alternative
        self.alternative_results = []
        for i in range(len(self.prior_info.alternatives)):
            if self.options.verbose:
                print(f"Solving alternative {i + 1}/{len(self.prior_info.alternatives)}...")
            result = self._solve_alternative(i)
            self.alternative_results.append(result)

        # Aggregate results
        self._aggregate_results()

        if self._result is not None:
            self._result.runtime = time.time() - start_time

        return self

    def getAvg(self) -> PosteriorResult:
        """Get average (prior-weighted) performance metrics."""
        if self._result is None:
            self.runAnalyzer()
        return self._result

    def getAvgTable(self) -> pd.DataFrame:
        """
        Get average performance metrics table.

        Returns:
            pandas.DataFrame with columns: Station, Class, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self._result is None:
            self.runAnalyzer()

        rows = []

        # Get dimensions
        if self._result.QN.size > 0:
            if self._result.QN.ndim == 2:
                M, K = self._result.QN.shape
            else:
                M = len(self._result.QN)
                K = 1
        else:
            M, K = 0, 0

        for i in range(M):
            for k in range(K):
                row = {
                    'Station': f'Station{i + 1}',
                    'Class': f'Class{k + 1}',
                }

                if self._result.QN.size > 0:
                    row['QLen'] = self._result.QN[i, k] if self._result.QN.ndim == 2 else self._result.QN[i]
                if self._result.UN.size > 0:
                    row['Util'] = self._result.UN[i, k] if self._result.UN.ndim == 2 else self._result.UN[i]
                if self._result.RN.size > 0:
                    row['RespT'] = self._result.RN[i, k] if self._result.RN.ndim == 2 else self._result.RN[i]
                    row['ResidT'] = row['RespT']
                if self._result.AN.size > 0:
                    row['ArvR'] = self._result.AN[i, k] if self._result.AN.ndim == 2 else self._result.AN[i]
                if self._result.TN.size > 0:
                    row['Tput'] = self._result.TN[i, k] if self._result.TN.ndim == 2 else self._result.TN[i]

                rows.append(row)

        df = pd.DataFrame(rows)

        if not self._table_silent and not df.empty:
            print(df.to_string(index=False))

        return df

    def getPosteriorTable(self) -> pd.DataFrame:
        """
        Get table with per-alternative results and probabilities.

        Returns:
            pandas.DataFrame with columns: Alternative, Probability, Station, Class, Q, U, R, T, A
        """
        if not self.alternative_results:
            self.runAnalyzer()

        if self.prior_info is None:
            return pd.DataFrame()

        rows = []
        probs = self.prior_info.probabilities

        for alt_idx, alt_result in enumerate(self.alternative_results):
            prob = probs[alt_idx]

            if 'QN' not in alt_result or alt_result['QN'].size == 0:
                continue

            if alt_result['QN'].ndim == 2:
                M, K = alt_result['QN'].shape
            else:
                M = len(alt_result['QN'])
                K = 1

            for i in range(M):
                for k in range(K):
                    row = {
                        'Alternative': alt_idx,
                        'Probability': prob,
                        'Station': f'Station{i + 1}',
                        'Class': f'Class{k + 1}',
                    }

                    def get_val(arr, i, k):
                        if arr.size == 0:
                            return np.nan
                        if arr.ndim == 2:
                            return arr[i, k]
                        return arr[i]

                    row['Q'] = get_val(alt_result.get('QN', np.array([])), i, k)
                    row['U'] = get_val(alt_result.get('UN', np.array([])), i, k)
                    row['R'] = get_val(alt_result.get('RN', np.array([])), i, k)
                    row['T'] = get_val(alt_result.get('TN', np.array([])), i, k)
                    row['A'] = get_val(alt_result.get('AN', np.array([])), i, k)

                    rows.append(row)

        df = pd.DataFrame(rows)
        return df

    def getPosteriorDist(self, metric: str, station_idx: int, class_idx: int = 0) -> EmpiricalCDF:
        """
        Get posterior distribution for a specific metric at a station/class.

        Args:
            metric: Metric name ('Q', 'U', 'R', 'T', 'A')
            station_idx: Station index (0-based)
            class_idx: Class index (0-based)

        Returns:
            EmpiricalCDF representing the posterior distribution
        """
        if not self.alternative_results:
            self.runAnalyzer()

        if self.prior_info is None:
            raise ValueError("No Prior distribution found in model")

        probs = self.prior_info.probabilities
        values = []

        metric_key = {
            'Q': 'QN', 'U': 'UN', 'R': 'RN', 'T': 'TN', 'A': 'AN', 'W': 'WN'
        }.get(metric.upper())

        if metric_key is None:
            raise ValueError(f"Unknown metric: {metric}")

        for alt_result in self.alternative_results:
            arr = alt_result.get(metric_key, np.array([]))
            if arr.size == 0:
                values.append(np.nan)
            elif arr.ndim == 2:
                values.append(arr[station_idx, class_idx])
            else:
                values.append(arr[station_idx])

        return EmpiricalCDF(values, probs)

    def getPercentile(self, metric: str, percentile: float,
                      station_idx: int, class_idx: int = 0) -> float:
        """
        Get percentile of posterior distribution for a metric.

        Args:
            metric: Metric name ('Q', 'U', 'R', 'T', 'A')
            percentile: Percentile value (0-100)
            station_idx: Station index (0-based)
            class_idx: Class index (0-based)

        Returns:
            Percentile value
        """
        cdf = self.getPosteriorDist(metric, station_idx, class_idx)
        return cdf.get_percentile(percentile / 100.0)

    @staticmethod
    def supports(model) -> bool:
        """Check if POSTERIOR solver supports the given model."""
        # POSTERIOR supports any model with a Prior distribution
        solver = SolverPosterior(model)
        return solver.hasPriorDistribution()

    @staticmethod
    def defaultOptions() -> PosteriorOptions:
        """Return default solver options."""
        return PosteriorOptions()

    @staticmethod
    def listValidMethods() -> List[str]:
        """List valid underlying solver types."""
        return ['mva', 'nc', 'ctmc', 'ssa', 'mam', 'fld', 'des']

    # Aliases
    run_analyzer = runAnalyzer
    get_avg = getAvg
    get_avg_table = getAvgTable
    get_posterior_table = getPosteriorTable
    get_posterior_dist = getPosteriorDist
    get_percentile = getPercentile


__all__ = [
    'SolverPosterior',
    'PosteriorOptions',
    'PosteriorResult',
    'PriorInfo',
    'EmpiricalCDF',
]
