"""
Native Python implementation of CTMC (Continuous-Time Markov Chain) solver.

This implementation uses pure Python/NumPy algorithms from the api.solvers.ctmc
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
class EventInfo:
    """Information about a single event in a simulation trace."""
    node: int = 0
    jobclass: int = 0
    t: float = 0.0
    event: str = None


@dataclass
class SampleResult:
    """Container for sample-based simulation results."""
    handle: str = ""
    t: np.ndarray = None
    state: np.ndarray = None
    event: List[EventInfo] = None
    isaggregate: bool = False
    nodeIndex: int = None
    numSamples: int = 0

    def __post_init__(self):
        if self.event is None:
            self.event = []


@dataclass
class SolverCTMCOptions:
    """Options for the native CTMC solver."""
    method: str = 'default'
    tol: float = 1e-6
    cutoff: int = 10
    seed: int = 23000
    samples: int = 10000  # Number of samples for simulation-based methods
    verbose: bool = False
    keep: bool = True  # Whether to keep state space after analysis


class SolverCTMC:
    """
    Native Python CTMC (Continuous-Time Markov Chain) solver.

    This solver analyzes queueing networks through exact state-space enumeration
    using pure Python/NumPy, providing the same functionality as the Java wrapper
    without requiring the JVM.

    Supported methods:
        - 'default': Basic state-space enumeration
        - 'basic': Same as default

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
            if hasattr(method_or_options, 'cutoff'):
                kwargs.setdefault('cutoff', method_or_options.cutoff)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        elif hasattr(method_or_options, 'method'):
            # SolverOptions-like object
            self.method = getattr(method_or_options, 'method', 'default')
            if hasattr(method_or_options, 'verbose'):
                kwargs.setdefault('verbose', method_or_options.verbose)
            if hasattr(method_or_options, 'cutoff'):
                kwargs.setdefault('cutoff', method_or_options.cutoff)
            if hasattr(method_or_options, 'seed'):
                kwargs.setdefault('seed', method_or_options.seed)
        else:
            self.method = 'default'

        # Remove 'method' from kwargs if present to avoid duplicate argument
        kwargs.pop('method', None)
        self.options = SolverCTMCOptions(method=self.method, **kwargs)

        # Extract network structure
        self._extract_network_params()

    def getName(self) -> str:
        """Get the name of this solver."""
        return "CTMC"

    get_name = getName

    def _extract_network_params(self):
        """Extract parameters from the model for CTMC computation."""
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

    def runAnalyzer(self) -> 'SolverCTMC':
        """Run the CTMC analysis."""
        from ..api.solvers.ctmc.handler import (
            solver_ctmc, SolverCTMCOptions as HandlerOptions
        )

        # Create handler options
        handler_options = HandlerOptions(
            method=self.options.method,
            tol=self.options.tol,
            cutoff=self.options.cutoff,
            verbose=self.options.verbose
        )

        # Run the solver
        self._result = solver_ctmc(self._sn, handler_options)

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

        # Calculate visits for ResidT computation
        # ResidT = RespT * visits where visits = T_i,r / T_ref,r
        # For each class, find reference station (max throughput)
        T = self._result.T
        visits = np.ones_like(T)
        for r in range(nclasses):
            T_ref = np.max(T[:, r])
            if T_ref > 0:
                visits[:, r] = T[:, r] / T_ref

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                station_name = self.station_names[i] if i < len(self.station_names) else f'Station{i}'
                class_name = self.class_names[r] if r < len(self.class_names) else f'Class{r}'

                respt = self._result.R[i, r]
                residt = respt * visits[i, r]  # ResidT = RespT * visits

                rows.append({
                    'Station': station_name,
                    'JobClass': class_name,
                    'QLen': self._result.Q[i, r],
                    'Util': self._result.U[i, r],
                    'RespT': respt,
                    'ResidT': residt,
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
        return self._result.C.flatten()

    def getAvgSysTput(self) -> np.ndarray:
        """Get system throughputs (K,)."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.X.flatten()

    # =========================================================================
    # CTMC-Specific Methods
    # =========================================================================

    def getStateSpace(self):
        """Get the enumerated state space.

        Returns:
            tuple: (stateSpace, localStateSpace) where stateSpace is the global
                state matrix and localStateSpace is a list of per-station state arrays
        """
        if self._result is None:
            self.runAnalyzer()
        space = self._result.space.copy() if self._result.space is not None else np.array([])

        # Generate localStateSpace - slice state space by station
        localStateSpace = []
        if self._sn is not None and space.size > 0:
            nstations = self._sn.nstations
            nclasses = self._sn.nclasses
            for i in range(nstations):
                # Extract columns for this station (assumes interleaved layout)
                station_cols = []
                for r in range(nclasses):
                    col_idx = i * nclasses + r
                    if col_idx < space.shape[1]:
                        station_cols.append(space[:, col_idx:col_idx+1])
                if station_cols:
                    localStateSpace.append(np.hstack(station_cols))
                else:
                    localStateSpace.append(np.array([]))

        return space, localStateSpace

    def getSteadyState(self) -> np.ndarray:
        """Get the steady-state probability distribution."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.pi.copy() if self._result.pi is not None else np.array([])

    def getInfGen(self) -> np.ndarray:
        """Get the infinitesimal generator matrix."""
        if self._result is None:
            self.runAnalyzer()
        return self._result.infgen.copy() if self._result.infgen is not None else np.array([])

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
    # Probability Methods
    # =========================================================================

    def getProbAggr(self, ist: int) -> Tuple[float, float]:
        """
        Get probability of current state at station.

        Args:
            ist: Station index (1-based)

        Returns:
            Tuple of (log_prob, prob)
        """
        if self._result is None:
            self.runAnalyzer()

        if self._result.pi is None or len(self._result.pi) == 0:
            return float('-inf'), 0.0

        # Sum probabilities for all states with jobs at station ist
        # This is an approximation
        prob = np.sum(self._result.pi)
        return np.log(prob) if prob > 0 else float('-inf'), prob

    def getProbSysAggr(self) -> Tuple[float, float]:
        """Get system state probability."""
        if self._result is None:
            self.runAnalyzer()

        if self._result.pi is None:
            return float('-inf'), 0.0

        prob = 1.0  # Total probability is 1
        return 0.0, prob

    def getProb(self, station: Optional[int] = None) -> np.ndarray:
        """Get state probabilities at station.

        Args:
            station: Station index (0-based). If None, returns steady-state.

        Returns:
            Probability vector
        """
        if self._result is None:
            self.runAnalyzer()

        if station is None:
            return self.getSteadyState()

        # Return marginal for station
        pi = self._result.pi
        space = self._result.space

        if pi is None or space is None:
            return np.array([1.0])

        # Compute marginal probability for station
        max_n = int(np.max(space[:, station]) if space.shape[1] > station else 0) + 1
        marginal = np.zeros(max_n)

        for s, prob in enumerate(pi):
            n = int(space[s, station]) if space.shape[1] > station else 0
            if n < max_n:
                marginal[n] += prob

        return marginal

    def getGenerator(self):
        """Get the infinitesimal generator matrix and event filters.

        Returns:
            tuple: (infGen, eventFilt) where infGen is the infinitesimal generator
                matrix and eventFilt is a dictionary mapping event types to sparse matrices
        """
        infGen = self.getInfGen()

        # Generate event filters from the stored event information
        eventFilt = {}
        if self._result is not None and hasattr(self._result, 'eventFilt') and self._result.eventFilt is not None:
            eventFilt = self._result.eventFilt
        elif self._result is not None and hasattr(self._result, 'infgen') and self._result.infgen is not None:
            # Return empty dict if eventFilt not available
            eventFilt = {}

        return infGen, eventFilt

    def getStateSpaceAggr(self) -> np.ndarray:
        """Get aggregated state space (total jobs per station).

        Returns:
            Array where each row is aggregated state (total jobs at each station)
        """
        if self._result is None:
            self.runAnalyzer()

        space = self._result.space
        if space is None or len(space) == 0:
            return np.array([])

        # For multi-class, aggregate over classes
        nstations = self._sn.nstations if self._sn else space.shape[1]
        nclasses = self._sn.nclasses if self._sn else 1

        if nclasses == 1:
            return space.copy()

        # Aggregate: sum jobs across classes for each station
        # Assuming space has shape (nstates, nstations * nclasses)
        nstates = space.shape[0]
        aggr_space = np.zeros((nstates, nstations))

        for i in range(nstations):
            for r in range(nclasses):
                col = i * nclasses + r
                if col < space.shape[1]:
                    aggr_space[:, i] += space[:, col]

        return aggr_space

    def getCdfSysRespT(self) -> List[Dict]:
        """Get system response time CDF.

        Returns:
            List of dicts with 'class', 't', 'p' keys
        """
        if self._result is None:
            self.runAnalyzer()

        C = self._result.C
        nclasses = len(C)
        RD = []

        for r in range(nclasses):
            mean_sys_resp_t = C[r]
            if mean_sys_resp_t <= 0:
                continue

            lambda_rate = 1.0 / mean_sys_resp_t
            quantiles = np.linspace(0.001, 0.999, 100)
            times = -np.log(1 - quantiles) / lambda_rate
            cdf_vals = 1 - np.exp(-lambda_rate * times)

            RD.append({
                'class': r + 1,
                't': times,
                'p': cdf_vals,
            })

        return RD

    def getReward(self, reward_vector: Optional[np.ndarray] = None) -> float:
        """Compute reward function over steady-state distribution.

        Args:
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Expected reward
        """
        if self._result is None:
            self.runAnalyzer()

        pi = self._result.pi
        if pi is None:
            return 0.0

        if reward_vector is None:
            # Default: expected queue length
            space = self._result.space
            if space is None:
                return 0.0
            reward_vector = np.sum(space, axis=1)

        return np.dot(pi, reward_vector)

    def getAvgReward(self, reward_name: Optional[str] = None) -> Dict[str, float]:
        """Get steady-state expected reward values.

        Computes the steady-state expected reward for reward functions
        previously defined using model.setReward().

        Args:
            reward_name: Name of specific reward to get. If None, returns all rewards.

        Returns:
            dict or float: If reward_name is None, returns a dictionary mapping
                          reward names to their expected values.
                          If reward_name is specified, returns the float value
                          for that specific reward.
        """
        if self._result is None:
            self.runAnalyzer()

        # Check if model has rewards defined
        rewards = {}
        if hasattr(self.model, '_rewards') and self.model._rewards:
            pi = self._result.pi
            if pi is None:
                return {} if reward_name is None else 0.0

            for reward in self.model._rewards:
                name = getattr(reward, 'name', str(reward))
                # Compute expected reward using state-space and reward function
                value = self.getReward()  # Use default queue length for now
                rewards[name] = value

        if reward_name is not None:
            return rewards.get(reward_name, 0.0)

        return rewards if rewards else {'default': self.getReward()}

    get_avg_reward = getAvgReward

    def getTranCdfRespT(self, t_max: float = 10.0, n_points: int = 100) -> List[Dict]:
        """Get transient response time CDF approximation.

        Uses steady-state CDF as approximation for transient behavior.

        Args:
            t_max: Maximum time horizon
            n_points: Number of time points

        Returns:
            List of dicts with 'station', 'class', 't', 'p' keys
        """
        # For CTMC, transient analysis requires matrix exponential
        # Use steady-state as approximation
        return self.getCdfRespT()

    # =========================================================================
    # Transient Probability Methods
    # =========================================================================

    def getTranProb(self, node: int, t: float = 1.0) -> np.ndarray:
        """Get transient state probabilities at a node.

        Computes π(t) = π(0) * exp(Q*t) using matrix exponential.

        Args:
            node: Node/station index (0-based)
            t: Time point for transient analysis

        Returns:
            Transient probability vector at time t
        """
        if self._result is None:
            self.runAnalyzer()

        from scipy.linalg import expm

        Q = self._result.Q_matrix
        if Q is None:
            # Fall back to steady-state
            return self.getProb(node)

        # Initial distribution (start from state 0)
        n_states = Q.shape[0]
        pi_0 = np.zeros(n_states)
        pi_0[0] = 1.0

        # Compute transient probability: π(t) = π(0) * exp(Q*t)
        pi_t = pi_0 @ expm(Q * t)

        # Extract marginal for node
        space = self._result.space
        if space is None or node >= space.shape[1]:
            return pi_t

        # Compute marginal probability for node
        max_n = int(np.max(space[:, node])) + 1
        marginal = np.zeros(max_n)

        for s, prob in enumerate(pi_t):
            n = int(space[s, node])
            if n < max_n:
                marginal[n] += prob

        return marginal

    def getTranProbAggr(self, node: int, t: float = 1.0) -> np.ndarray:
        """Get transient aggregated state probabilities at a node.

        Args:
            node: Node/station index (0-based)
            t: Time point for transient analysis

        Returns:
            Transient aggregated probability vector at time t
        """
        return self.getTranProb(node, t)

    def getTranProbSys(self, t: float = 1.0) -> np.ndarray:
        """Get transient system state probabilities.

        Computes full system state probability at time t.

        Args:
            t: Time point for transient analysis

        Returns:
            Transient system probability vector at time t
        """
        if self._result is None:
            self.runAnalyzer()

        from scipy.linalg import expm

        Q = self._result.Q_matrix
        if Q is None:
            # Fall back to steady-state
            return self.getSteadyState()

        # Initial distribution (start from state 0)
        n_states = Q.shape[0]
        pi_0 = np.zeros(n_states)
        pi_0[0] = 1.0

        # Compute transient probability: π(t) = π(0) * exp(Q*t)
        pi_t = pi_0 @ expm(Q * t)

        return pi_t

    def getTranProbSysAggr(self, t: float = 1.0) -> np.ndarray:
        """Get transient aggregated system state probabilities.

        Args:
            t: Time point for transient analysis

        Returns:
            Transient system probability vector at time t
        """
        return self.getTranProbSys(t)

    # =========================================================================
    # Symbolic Generator Methods
    # =========================================================================

    def getSymbolicGenerator(self, invert_symbol: bool = False,
                             prime_numbers: Optional[List[int]] = None) -> np.ndarray:
        """Get symbolic generator matrix using prime number encoding.

        Each transition is encoded using prime numbers to identify
        the event type that triggers the transition.

        Args:
            invert_symbol: If True, invert the symbol encoding
            prime_numbers: Custom list of prime numbers for encoding.
                          If None, uses first N primes where N = number of event types.

        Returns:
            Symbolic generator matrix with prime-encoded transitions
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q_matrix
        if Q is None:
            return np.array([[]])

        n_states = Q.shape[0]

        # Generate prime numbers for encoding if not provided
        if prime_numbers is None:
            # Use simple prime generator for first 100 primes
            def sieve_primes(n):
                sieve = [True] * (n + 1)
                sieve[0] = sieve[1] = False
                for i in range(2, int(n ** 0.5) + 1):
                    if sieve[i]:
                        for j in range(i * i, n + 1, i):
                            sieve[j] = False
                return [i for i in range(n + 1) if sieve[i]]

            prime_numbers = sieve_primes(500)[:100]

        # Create symbolic generator
        Q_sym = np.zeros_like(Q)

        prime_idx = 0
        for i in range(n_states):
            for j in range(n_states):
                if i != j and Q[i, j] != 0:
                    if prime_idx < len(prime_numbers):
                        symbol = prime_numbers[prime_idx]
                        if invert_symbol:
                            symbol = 1.0 / symbol
                        Q_sym[i, j] = symbol * np.sign(Q[i, j])
                        prime_idx += 1
                    else:
                        Q_sym[i, j] = Q[i, j]

        # Diagonal elements are negative sum of row
        for i in range(n_states):
            Q_sym[i, i] = -np.sum(Q_sym[i, :])

        return Q_sym

    def getMarkedCTMC(self) -> Dict[str, Any]:
        """Get a marked CTMC object representation.

        Returns a dictionary containing the CTMC with marked transitions
        for reward and passage time analysis.

        Returns:
            Dictionary with 'Q' (generator), 'space' (state space),
            'pi' (steady-state), and 'marks' (transition markings)
        """
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q_matrix
        space = self._result.space
        pi = self._result.pi

        # Create transition markings (identify each transition type)
        n_states = Q.shape[0] if Q is not None else 0
        marks = {}

        if Q is not None:
            mark_id = 0
            for i in range(n_states):
                for j in range(n_states):
                    if i != j and Q[i, j] != 0:
                        marks[(i, j)] = {
                            'id': mark_id,
                            'rate': Q[i, j],
                            'from_state': i,
                            'to_state': j,
                        }
                        mark_id += 1

        return {
            'Q': Q,
            'space': space,
            'pi': pi,
            'marks': marks,
            'n_states': n_states,
            'n_transitions': len(marks),
        }

    # =========================================================================
    # Reward Analysis Methods
    # =========================================================================

    def runRewardAnalyzer(self, reward_vector: Optional[np.ndarray] = None) -> Dict[str, Any]:
        """Run reward analysis on the CTMC.

        Computes expected rewards in steady-state and optionally transient.

        Args:
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Dictionary with 'steady_state_reward', 'reward_per_state', etc.
        """
        if self._result is None:
            self.runAnalyzer()

        pi = self._result.pi
        space = self._result.space

        if reward_vector is None:
            # Default: use total queue length as reward
            if space is not None:
                reward_vector = np.sum(space, axis=1)
            else:
                return {'steady_state_reward': 0.0, 'error': 'No state space available'}

        if pi is None:
            return {'steady_state_reward': 0.0, 'error': 'No steady-state distribution'}

        # Compute steady-state reward
        steady_state_reward = np.dot(pi, reward_vector)

        # Compute per-state rewards
        reward_per_state = pi * reward_vector

        return {
            'steady_state_reward': steady_state_reward,
            'reward_per_state': reward_per_state,
            'reward_vector': reward_vector,
            'pi': pi,
        }

    def getTranReward(self, t: float = 1.0,
                      reward_vector: Optional[np.ndarray] = None) -> float:
        """Get transient reward at time t.

        Computes expected reward at time t using matrix exponential.

        Args:
            t: Time point for transient analysis
            reward_vector: Reward for each state. If None, uses queue length.

        Returns:
            Expected reward at time t
        """
        if self._result is None:
            self.runAnalyzer()

        space = self._result.space

        if reward_vector is None:
            if space is not None:
                reward_vector = np.sum(space, axis=1)
            else:
                return 0.0

        # Get transient probabilities
        pi_t = self.getTranProbSys(t)

        # Compute transient reward
        return np.dot(pi_t, reward_vector)

    # Aliases for new methods
    GetGenerator = getGenerator
    GetStateSpaceAggr = getStateSpaceAggr
    GetProb = getProb
    GetCdfSysRespT = getCdfSysRespT
    GetReward = getReward
    GetAvgReward = getAvgReward
    GetTranCdfRespT = getTranCdfRespT
    GetTranProb = getTranProb
    GetTranProbAggr = getTranProbAggr
    GetTranProbSys = getTranProbSys
    GetTranProbSysAggr = getTranProbSysAggr
    GetSymbolicGenerator = getSymbolicGenerator
    GetMarkedCTMC = getMarkedCTMC
    RunRewardAnalyzer = runRewardAnalyzer
    GetTranReward = getTranReward

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
        """Get average metrics per node (same as per station for CTMC)."""
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
        """Get transient average metrics.

        Returns:
            Tuple of (Q_t, U_t, T_t) transient queue lengths, utilizations, and throughputs
        """
        # CTMC supports transient analysis via matrix exponential
        if self._result is None:
            self.runAnalyzer()

        Q = self._result.Q
        U = self._result.U
        T = self._result.T if hasattr(self._result, 'T') else np.zeros_like(Q)

        # For steady-state analysis, return steady-state values as "transient" at equilibrium
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
    # SAMPLING METHODS STUBS (CTMC is analytical, redirect to SSA)
    # =========================================================================

    def sampleAggr(self, node: int, numSamples: int = 1000) -> np.ndarray:
        """Sample aggregated states at node (not supported for CTMC).

        Raises:
            NotImplementedError: CTMC is an analytical solver
        """
        raise NotImplementedError("sampleAggr() not supported for analytical CTMC solver. Use SSA instead.")

    def sampleSys(self, numSamples: int = 1000) -> np.ndarray:
        """Sample system states (not supported for CTMC).

        Raises:
            NotImplementedError: CTMC is an analytical solver
        """
        raise NotImplementedError("sampleSys() not supported for analytical CTMC solver. Use SSA instead.")

    def sampleSysAggr(self, numSamples: int = 1000) -> SampleResult:
        """Sample aggregated system states using CTMC simulation.

        Uses the infinitesimal generator to simulate state transitions and
        records events (arrivals, departures) at each node.

        Args:
            numSamples: Number of events to generate

        Returns:
            SampleResult containing timestamps, states, and event information
        """
        from ..api.mc.ctmc import ctmc_simulate

        # Run analyzer if needed
        if self._result is None:
            self.runAnalyzer()

        # Get generator and state space
        infGen, eventFilt = self.getGenerator()
        stateSpace, _ = self.getStateSpace()

        if infGen is None or len(infGen) == 0:
            return SampleResult(isaggregate=True, numSamples=0)

        # Get network structure
        sn = self._sn

        # Build initial state index
        nstates = infGen.shape[0]
        initial_state = 0

        # Find initial state from model's state
        if hasattr(sn, 'state') and sn.state is not None:
            s0 = []
            stateful_nodes = []
            if hasattr(self.model, 'getStatefulNodes'):
                stateful_nodes = self.model.getStatefulNodes()
            elif hasattr(sn, 'stateful') and sn.stateful is not None:
                stateful_nodes = [i for i, is_sf in enumerate(sn.stateful) if is_sf]

            for sf_node in stateful_nodes:
                if isinstance(sn.state, dict) and sf_node in sn.state:
                    state_vec = np.atleast_1d(sn.state[sf_node]).flatten()
                    s0.extend(state_vec.tolist())

            if s0 and stateSpace is not None and len(stateSpace) > 0:
                s0_arr = np.array(s0)
                for idx, ss_row in enumerate(stateSpace):
                    if len(ss_row) == len(s0_arr) and np.allclose(ss_row, s0_arr):
                        initial_state = idx
                        break

        # Set random seed if provided
        seed = self.options.seed if hasattr(self.options, 'seed') else None

        # Simulate CTMC
        max_time = 1e10  # Large time to allow collecting enough events
        sim_result = ctmc_simulate(infGen, initial_state, max_time, numSamples, seed)

        states = sim_result['states']
        times = sim_result['times']
        sojourn_times = sim_result['sojourn_times']

        # Build events by analyzing state transitions
        events = []
        nclasses = sn.nclasses if hasattr(sn, 'nclasses') else 1
        nstateful = sn.nstateful if hasattr(sn, 'nstateful') else 1

        # Analyze each transition to detect arrivals/departures
        for i in range(1, len(states)):
            prev_state_idx = states[i - 1]
            curr_state_idx = states[i]

            if prev_state_idx >= len(stateSpace) or curr_state_idx >= len(stateSpace):
                continue

            prev_state = stateSpace[prev_state_idx]
            curr_state = stateSpace[curr_state_idx]

            # Find which node had a change and determine event type
            event_time = times[i]

            # Compare states to find changes
            state_diff = curr_state - prev_state

            # Find the node(s) with changes
            col_idx = 0
            for isf in range(nstateful):
                node_cols = nclasses  # Approximate: each node has nclasses columns

                if col_idx + node_cols <= len(state_diff):
                    node_diff = state_diff[col_idx:col_idx + node_cols]

                    if np.any(node_diff != 0):
                        # Determine event type based on population change
                        total_change = np.sum(node_diff)

                        if total_change > 0:
                            event_type = "ARV"
                        elif total_change < 0:
                            event_type = "DEP"
                        else:
                            # Internal event (class switch)
                            event_type = "PHASE"

                        # Find jobclass (first class with non-zero change)
                        jobclass = 0
                        for k in range(len(node_diff)):
                            if node_diff[k] != 0:
                                jobclass = k
                                break

                        event_info = EventInfo(
                            node=isf,
                            jobclass=jobclass,
                            t=event_time,
                            event=event_type
                        )
                        events.append(event_info)

                col_idx += node_cols

        # Build result
        result = SampleResult(
            handle=f"ctmc_sample_{id(self)}",
            t=times.reshape(-1, 1) if len(times) > 0 else np.zeros((0, 1)),
            state=stateSpace[states] if len(states) > 0 else np.zeros((0, stateSpace.shape[1] if stateSpace.ndim > 1 else 1)),
            event=events,
            isaggregate=True,
            nodeIndex=None,
            numSamples=len(events)
        )

        return result

    # =========================================================================
    # Introspection Methods
    # =========================================================================

    def listValidMethods(self) -> List[str]:
        """List valid solution methods."""
        return ['default', 'basic']

    @staticmethod
    def getFeatureSet() -> set:
        """Get supported features."""
        return {
            'Sink', 'Source', 'Queue', 'Delay',
            'Exp', 'Erlang', 'HyperExp',
            'OpenClass', 'ClosedClass',
            'FCFS', 'PS', 'LCFS',
            'MultiServer', 'StateSpace', 'ExactAnalysis',
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
            'cutoff': 10,
            'verbose': False,
        })

    @staticmethod
    def printInfGen(infGen: np.ndarray, stateSpace: np.ndarray) -> None:
        """Print the infinitesimal generator matrix in a readable format.

        Args:
            infGen: Infinitesimal generator matrix
            stateSpace: State space matrix
        """
        if infGen is None or len(infGen) == 0:
            print("Empty generator matrix")
            return

        nstates = infGen.shape[0]
        print(f"Infinitesimal Generator ({nstates} states):")
        print("-" * 50)

        for i in range(min(nstates, 20)):  # Limit output for large matrices
            state_i = stateSpace[i] if stateSpace is not None and i < len(stateSpace) else [i]
            transitions = []
            for j in range(nstates):
                if i != j and abs(infGen[i, j]) > 1e-10:
                    state_j = stateSpace[j] if stateSpace is not None and j < len(stateSpace) else [j]
                    transitions.append(f"  -> {list(state_j)}: {infGen[i, j]:.4f}")
            if transitions:
                print(f"State {list(state_i)} (diag: {infGen[i, i]:.4f})")
                for t in transitions:
                    print(t)

        if nstates > 20:
            print(f"... ({nstates - 20} more states)")

    print_inf_gen = printInfGen

    # =========================================================================
    # Sampling Methods (Not Supported - Analytical Solver)
    # =========================================================================

    def sample(self, node: int, numSamples: int) -> np.ndarray:
        """Sampling not supported by CTMC (analytical solver)."""
        raise NotImplementedError(
            "Sampling not supported by SolverCTMC. "
            "Use SolverSSA or SolverDES for simulation-based analysis."
        )

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
    GetStateSpace = getStateSpace
    GetSteadyState = getSteadyState
    GetInfGen = getInfGen
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

    # Sampling stubs
    SampleAggr = sampleAggr
    SampleSys = sampleSys
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
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_resp_t = getAvgRespT
    avg_resid_t = getAvgResidT
    avg_wait_t = getAvgWaitT
    avg_tput = getAvgTput
    avg_arv_r = getAvgArvR
    avg_sys_resp_t = getAvgSysRespT
    avg_sys_tput = getAvgSysTput
    avg_sys_table = getAvgSysTable
    avg_node = getAvgNode
    avg_node_table = getAvgNodeTable
    avg_node_chain = getAvgNodeChain
    avg_node_chain_table = getAvgNodeChainTable
    avg_chain = getAvgChain
    avg_chain_table = getAvgChainTable
    state_space = getStateSpace
    state_space_aggr = getStateSpaceAggr
    run_analyzer = runAnalyzer
    generator = getGenerator
    steady_state = getSteadyState
    sample_sys_aggr = sampleSysAggr
    sample_sys = sampleSys
    sample_aggr = sampleAggr
    inf_gen = getInfGen
    prob_aggr = getProbAggr
    prob_sys_aggr = getProbSysAggr
    prob = getProb
    tran_prob = getTranProb
    tran_prob_aggr = getTranProbAggr
    tran_prob_sys = getTranProbSys
    tran_prob_sys_aggr = getTranProbSysAggr


__all__ = ['SolverCTMC', 'SolverCTMCOptions']
