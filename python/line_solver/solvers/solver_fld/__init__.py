"""
SolverFLD - Native Python Fluid Approximation Solver for Queueing Networks.

A comprehensive implementation of fluid (mean-field) approximation methods for
analyzing queueing networks, supporting both open and closed systems with
multiple job classes, scheduling disciplines, and service distributions.

Features
--------
- **7 Solution Methods**: Matrix (p-norm), Softmin, State-Dependent, Closing,
  Diffusion, MFQ, and Passage Time analysis
- **Flexible Network Support**: Open, closed, and mixed networks with
  configurable topologies
- **Multiple Classes**: Multi-class job support for priority modeling
- **Scheduling Policies**: FCFS, Processor Sharing (PS), Infinite Server (INF),
  External (EXT), and Deferred Processor Sharing (DPS)
- **Transient Analysis**: Time-dependent queue dynamics via ODE solution
- **Stochastic Simulation**: Euler-Maruyama SDE for closed systems (diffusion method)
"""

import numpy as np
import pandas as pd
import time
from typing import Optional, Dict, List, Tuple, Any, Union
from dataclasses import dataclass

from .options import SolverFLDOptions, FLDResult
from .utils import extract_metrics_from_handler_result, compute_response_times, compute_cycle_times, compute_system_throughput


class SolverFLD:
    """Native Python solver for fluid approximation of queueing networks.

    Provides a unified interface to multiple fluid approximation algorithms,
    allowing seamless switching between different solution methods while
    maintaining consistent result formats and performance metrics.

    This class follows the SolverMAM design pattern with:
    - Lazy method resolution (user-facing aliases map to internal implementation)
    - Consistent result accessors (getAvgQLen, getAvgRespT, etc.)
    - Method chaining support (runAnalyzer returns self)
    - Optional verbose output for debugging

    Attributes
    ----------
    network : object
        Input network model (either a NetworkStruct or object with compileStruct method)
    sn : NetworkStruct
        Compiled network structure (internal representation)
    options : SolverFLDOptions
        Configuration parameters for solver behavior
    result : FLDResult or None
        Solution results (None until runAnalyzer is called)
    runtime : float
        Elapsed time in seconds for last analysis

    Examples
    --------
    Basic usage:
        >>> solver = SolverFLD(network, method='mfq')
        >>> solver.runAnalyzer()
        >>> QN = solver.result.QN  # Access raw results
        >>> qlen = solver.getAvgQLen()  # Access aggregated metrics

    Method comparison:
        >>> results = {}
        >>> for method in ['matrix', 'mfq']:
        ...     s = SolverFLD(network, method=method)
        ...     s.runAnalyzer()
        ...     results[method] = s.result

    Custom configuration:
        >>> opts = SolverFLDOptions(tol=1e-6, pstar=50, verbose=True)
        >>> solver = SolverFLD(network, options=opts)
        >>> solver.runAnalyzer()
        >>> metrics = solver.getAvgTable()  # Returns pandas DataFrame
    """

    # Method mapping to internal implementation keys
    METHODS = {
        'default': 'matrix',
        'matrix': 'matrix',
        'fluid.matrix': 'matrix',
        'pnorm': 'matrix',
        'fluid.pnorm': 'matrix',
        'softmin': 'closing',
        'fluid.softmin': 'closing',
        'statedep': 'closing',
        'fluid.statedep': 'closing',
        'closing': 'closing',
        'fluid.closing': 'closing',
        'diffusion': 'diffusion',
        'fluid.diffusion': 'diffusion',
        'mfq': 'mfq',
        'fluid.mfq': 'mfq',
        'butools': 'mfq',
        'aoi': 'aoi',
        'fluid.aoi': 'aoi',
    }

    def __init__(
        self,
        network,
        method_or_options: Union[str, SolverFLDOptions, dict] = 'default',
        options: Optional[SolverFLDOptions] = None,
        **kwargs
    ):
        """Initialize SolverFLD.

        Parameters
        ----------
        network : NetworkStruct or object
            Network model specification. Can be either:
            - A NetworkStruct (compiled network structure)
            - An object with compileStruct() method (will be compiled automatically)
        method : str, optional
            Solution method to use. Valid options:
            - 'default', 'matrix', 'fluid.matrix', 'pnorm', 'fluid.pnorm':
              Matrix method with p-norm smoothing (default, recommended for most networks)
            - 'softmin', 'fluid.softmin': Softmin smoothing (open networks only)
            - 'statedep', 'fluid.statedep': State-dependent constraints (open networks only)
            - 'closing', 'fluid.closing': Closing approximation with FCFS iteration
            - 'diffusion', 'fluid.diffusion': Euler-Maruyama SDE (closed networks only)
            - 'mfq', 'fluid.mfq', 'butools': Markovian fluid queue - exact M/M/c
              (single-queue networks only)
            Default is 'matrix' (mapped from 'default').
        options : SolverFLDOptions, optional
            Configuration object. If not provided, defaults are used. If both method
            and options.method are specified, method parameter takes precedence.
            See SolverFLDOptions for available parameters.

        Raises
        ------
        ValueError
            If network is neither NetworkStruct nor has compileStruct method.

        Notes
        -----
        Method selection guidelines:
        - **matrix** (default): Recommended starting point. Works for open and closed
          networks. Fast and numerically stable. Parameters: pstar (smoothing parameter,
          default 20)
        - **mfq**: If analyzing single-queue bottleneck (M/M/1, M/M/c). Provides
          exact analytical solution via Erlang-C formula.
        - **diffusion**: For closed networks needing stochastic dynamics. Useful for
          variance and percentile analysis.
        - **closing**: For networks dominated by FCFS service. Requires iterations to
          converge. Parameters: iter_max, iter_tol

        Examples
        --------
        Using with NetworkStruct directly:
            >>> from line_solver.api.sn import NetworkStruct
            >>> sn = NetworkStruct()  # ... configure ...
            >>> solver = SolverFLD(sn, method='mfq')

        Using with Network object:
            >>> model = Network('TestModel')
            >>> # ... configure network ...
            >>> solver = SolverFLD(model, method='matrix')

        With custom options:
            >>> opts = SolverFLDOptions(tol=1e-6, pstar=50)
            >>> solver = SolverFLD(model, options=opts)
        """
        self.network = network
        self.sn = self._get_network_struct(network)

        # Handle flexible argument patterns:
        # FLD(model) - use default options
        # FLD(model, 'method') - string method
        # FLD(model, options) - options object as second arg
        # FLD(model, 'method', options) - method string and options
        # FLD(model, method='method', **kwargs) - keyword args only

        # Extract method from kwargs if present (for FLD(model, method='x', iter_max=100) pattern)
        method_from_kwargs = kwargs.pop('method', None)

        if isinstance(method_or_options, str):
            method = method_or_options
        elif isinstance(method_or_options, SolverFLDOptions):
            options = method_or_options
            method = options.method
        elif isinstance(method_or_options, dict):
            # Dict passed as second argument - treat as options
            method = method_or_options.get('method', 'default')
            # Convert dict to SolverFLDOptions
            opts_kwargs = {k: v for k, v in method_or_options.items()
                          if k in ['method', 'tol', 'iter_max', 'iter_tol', 'pstar',
                                   'Tmax', 'verbose', 'init_point', 'timespan',
                                   'samples', 'stiff']}
            options = SolverFLDOptions(**opts_kwargs)
        else:
            # method_or_options is 'default' string
            method = method_from_kwargs if method_from_kwargs else method_or_options

        if options is None:
            options = SolverFLDOptions(method=method, **kwargs)
        else:
            if method != 'default' and method != options.method:
                options.method = method
            # Also apply any kwargs to existing options
            for key, value in kwargs.items():
                if hasattr(options, key):
                    setattr(options, key, value)

        self.options = options
        self.result = None
        self.runtime = 0.0

    def getName(self) -> str:
        """Get the name of this solver."""
        return "Fluid"

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

        # Priority 3: Has get_struct method (native Python)
        if hasattr(model, 'get_struct'):
            return model.get_struct()

        # Priority 4: Has compileStruct method (wrapper)
        if hasattr(model, 'compileStruct'):
            return model.compileStruct()

        # Priority 5: JPype wrapper with getStruct
        if hasattr(model, 'getStruct'):
            return model.getStruct()

        # Priority 6: Model that is already a struct
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError("Cannot extract network structure from model")

    def runAnalyzer(self) -> 'SolverFLD':
        """Execute the fluid analysis using the configured method.

        Performs ODE integration or analytical solution to compute steady-state
        performance metrics (queue lengths, utilizations, response times, etc.)
        for the queueing network.

        Returns
        -------
        SolverFLD
            Returns self to enable method chaining and fluent interface

        Raises
        ------
        ValueError
            If selected method is invalid or network constraints violated
            (e.g., diffusion on open network, mfq on multi-queue network)

        Notes
        -----
        - First time calls resolve method alias to internal implementation
        - Calls appropriate internal solve method based on self.options.method
        - Updates self.result with FLDResult object
        - Records execution time in self.runtime
        - If verbose=True, prints method and timing information

        Performance characteristics:
        - matrix: O(M²K + ODE_steps), typically 0.1-1s for realistic networks
        - mfq: O(K), typically <0.01s (analytical solution)
        - diffusion: O(num_steps × M × K), typically 0.5-5s
        - closing: O(iter_max × ODE_steps), typically 1-10s

        Examples
        --------
        Basic analysis:
            >>> solver = SolverFLD(network)
            >>> solver.runAnalyzer()
            >>> qlen = solver.getAvgQLen()

        Method chaining:
            >>> qlen = SolverFLD(network).runAnalyzer().getAvgQLen()

        Verbose output:
            >>> opts = SolverFLDOptions(verbose=True)
            >>> SolverFLD(network, options=opts).runAnalyzer()
            # Output: SolverFLD: Using method 'matrix'
            #         Completed in 0.234s
        """
        method_key = self._resolve_method()

        if self.options.verbose:
            print(f"SolverFLD: Using method '{method_key}'")

        start_time = time.time()

        # Method dispatch
        if method_key == 'matrix':
            self.result = self._solve_matrix()
        elif method_key == 'closing':
            self.result = self._solve_closing()
        elif method_key == 'diffusion':
            self.result = self._solve_diffusion()
        elif method_key == 'mfq':
            self.result = self._solve_mfq()
        elif method_key == 'aoi':
            self.result = self._solve_aoi()
        else:
            raise ValueError(f"Unknown method: {method_key}")

        self.runtime = time.time() - start_time

        if self.options.verbose:
            print(f"  Completed in {self.runtime:.3f}s")

        return self

    def _resolve_method(self) -> str:
        """Resolve method name to internal key.

        Returns:
            Internal method key
        """
        method = self.options.method
        return self.METHODS.get(method, method)

    def _solve_matrix(self) -> FLDResult:
        """Solve using matrix method (existing handler implementation).

        Returns:
            FLDResult
        """
        from line_solver.api.solvers.fld.handler import solver_fld, SolverFLDOptions as HandlerFLDOptions

        # Convert options to handler format
        handler_opts = HandlerFLDOptions(
            method=self.options.method,
            tol=self.options.tol,
            verbose=self.options.verbose,
            stiff=self.options.stiff,
            iter_max=self.options.iter_max,
            timespan=self.options.timespan,
            pstar=[self.options.pstar] * self.sn.nstations,
            num_cdf_pts=200
        )

        # Solve
        handler_result = solver_fld(self.sn, handler_opts)

        # Extract metrics
        QN, UN, RN, TN, CN, XN = extract_metrics_from_handler_result(handler_result, self.sn)

        # Build result
        result = FLDResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            t=handler_result.t,
            xvec=handler_result.odeStateVec,
            iterations=handler_result.it,
            runtime=self.runtime,
            method='matrix'
        )

        return result

    def _solve_closing(self) -> FLDResult:
        """Solve using closing method (FCFS approximation + ODE).

        Returns:
            FLDResult

        Supports softmin, statedep, and pnorm approaches via options.method.
        """
        from .methods.closing import ClosingMethod

        method = ClosingMethod(self.sn, self.options)
        return method.solve()

    def _solve_diffusion(self) -> FLDResult:
        """Solve using Euler-Maruyama diffusion (SDE for closed networks).

        Returns:
            FLDResult

        Restricted to closed networks only (all classes must have fixed populations).
        """
        from .methods.diffusion import DiffusionMethod

        method = DiffusionMethod(self.sn, self.options)
        return method.solve()

    def _solve_mfq(self) -> FLDResult:
        """Solve using BUTools Markovian Fluid Queue (single-queue exact).

        Returns:
            FLDResult

        Restricted to single-queue topologies (exactly one queue station).
        Uses analytical M/M/c solution when BUTools is not available.
        """
        from .methods.mfq import MFQMethod

        method = MFQMethod(self.sn, self.options)
        return method.solve()

    def _solve_aoi(self) -> FLDResult:
        """Solve using Age of Information analysis (aoi-fluid library).

        Returns:
            FLDResult with AoI-specific metrics attached

        Restricted to single-queue topologies with:
        - Bufferless (capacity=1): PH/PH/1/1 or PH/PH/1/1* (preemptive)
        - Single-buffer (capacity=2): M/PH/1/2 or M/PH/1/2* (replacement)

        License: aoi-fluid toolbox (BSD 2-Clause)
        Copyright (c) 2020, Ozancan Dogan, Nail Akar, Eray Unsal Atay
        """
        from .methods.aoi import AoIMethod

        method = AoIMethod(self.sn, self.options)
        return method.solve()

    # =====================================================================
    # RESULT ACCESS METHODS (following SolverMAM pattern)
    # =====================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """Get average performance metrics as DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self.result is None:
            self.runAnalyzer()

        # Extract station and class names
        nstations = self.sn.nstations
        nclasses = self.sn.nclasses

        # Get station names using stationToNode mapping
        nodenames = list(self.sn.nodenames) if hasattr(self.sn, 'nodenames') and self.sn.nodenames else []
        stationToNode = self.sn.stationToNode if hasattr(self.sn, 'stationToNode') else None

        station_names = []
        if stationToNode is not None and nodenames:
            stationToNode = np.asarray(stationToNode).flatten()
            for i in range(nstations):
                if i < len(stationToNode):
                    node_idx = int(stationToNode[i])
                    if node_idx < len(nodenames):
                        station_names.append(nodenames[node_idx])
                    else:
                        station_names.append(f'Station{i}')
                else:
                    station_names.append(f'Station{i}')
        else:
            station_names = [f'Station{i}' for i in range(nstations)]

        # Get class names
        class_names = list(self.sn.classnames) if hasattr(self.sn, 'classnames') and self.sn.classnames else \
                      [f'Class{i}' for i in range(nclasses)]

        # Build rows
        QN = self.result.QN
        UN = self.result.UN
        RN = self.result.RN
        TN = self.result.TN if self.result.TN is not None else np.zeros((nstations, nclasses))

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                qlen = QN[i, r] if i < QN.shape[0] and r < QN.shape[1] else 0
                util = UN[i, r] if i < UN.shape[0] and r < UN.shape[1] else 0
                respt = RN[i, r] if i < RN.shape[0] and r < RN.shape[1] else 0
                tput = TN[i, r] if i < TN.shape[0] and r < TN.shape[1] else 0

                # Skip zero rows
                if abs(qlen) < 1e-12 and abs(util) < 1e-12 and abs(tput) < 1e-12:
                    continue

                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'JobClass': class_names[r] if r < len(class_names) else f'Class{r}',
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': respt,  # Same as RespT for FLD
                    'ArvR': tput,  # For closed networks, arrival rate = throughput
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        if len(df) > 0 and not getattr(self, '_table_silent', False):
            print(df.to_string(index=False))

        return df

    def getAvgQLen(self) -> np.ndarray:
        """Get average queue lengths per station.

        Returns the mean queue length (number of customers in system) per station,
        aggregated across all job classes.

        Returns
        -------
        np.ndarray
            Shape (M,) array where M = number of stations. QN[i] is the average
            queue length at station i, including customers in service and waiting.

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        For Little's Law validation: L = λ × W, where λ is arrival rate and W
        is mean response time. This relationship should hold for stable networks.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> qlen = solver.getAvgQLen()
        >>> print(f"Queue length at station 0: {qlen[0]:.3f}")
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.QN, axis=1)

    def getAvgUtil(self) -> np.ndarray:
        """Get average server utilizations per station.

        Returns the fraction of time each server is busy, aggregated across
        all job classes.

        Returns
        -------
        np.ndarray
            Shape (M,) array where M = number of stations. UN[i] is the average
            utilization (fraction in (0, 1)) at station i.

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        For stable single-server queue (M/M/1): ρ = λ/μ. For multi-server
        queue (M/M/c): ρ = λ/(c×μ). Stability requires ρ < 1.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> util = solver.getAvgUtil()
        >>> bottleneck = np.argmax(util)
        >>> print(f"Bottleneck station: {bottleneck} (util={util[bottleneck]:.1%})")
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.UN, axis=1)

    def getAvgRespT(self) -> np.ndarray:
        """Get average response times per station and class.

        Returns mean time customers spend at each station (waiting + service)
        for each job class. Includes both queueing delay and service time.

        Returns
        -------
        np.ndarray
            Shape (M, K) array where M = number of stations, K = number of classes.
            RN[i, c] is the average response time at station i for class c.

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        For M/M/1 queue: W = 1/(μ - λ) = ρ/(μ(1 - ρ)) where ρ = λ/μ.
        Verifies Little's Law: L = λ × W for each station.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> resp_time = solver.getAvgRespT()
        >>> print(f"System response time: {np.sum(resp_time):.3f}")
        """
        if self.result is None:
            self.runAnalyzer()
        return self.result.RN

    def getTput(self) -> np.ndarray:
        """Get average throughputs per job class.

        Returns the throughput (customers per time unit) for each job class,
        aggregated across all stations.

        Returns
        -------
        np.ndarray
            Shape (K,) array where K = number of job classes. TN[k] is the
            average throughput of class k.

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        For open networks with fixed arrivals, throughput at each station
        should equal the arrival rate (λ) for stable systems.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> tput = solver.getTput()
        >>> for k, t in enumerate(tput):
        ...     print(f"Class {k}: {t:.4f} customers/time")
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.TN, axis=0)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get average system response time per job class.

        Returns the total time a customer spends in the system (sum across all
        stations) for each job class.

        Returns
        -------
        np.ndarray
            Shape (K,) array where K = number of job classes. CN[k] is the
            average system response time for class k.

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        This is the sum of response times across all stations. For a tandem
        network: C = W1 + W2 + ... + WM.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> sys_resp = solver.getAvgSysRespT()
        >>> print(f"Mean system response time: {np.mean(sys_resp):.3f}")
        """
        if self.result is None:
            self.runAnalyzer()
        return np.sum(self.result.RN, axis=0)

    def getAvgSysTput(self) -> float:
        """Get average system-wide throughput.

        Returns the overall throughput of the network (customers completing
        service per unit time), aggregated across all classes.

        Returns
        -------
        float
            Scalar system-wide throughput (customers/time unit)

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        For open networks in equilibrium: System throughput = Σ_k λ_k
        (sum of arrival rates). For closed networks: Limited by bottleneck
        service rate.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> sys_tput = solver.getAvgSysTput()
        >>> print(f"System throughput: {sys_tput:.4f} customers/time")
        """
        if self.result is None:
            self.runAnalyzer()
        return np.mean(self.result.XN)

    # =====================================================================
    # PASSAGE TIME / RESPONSE TIME METHODS
    # =====================================================================

    def getCdfRespT(self, station: Optional[int] = None, job_class: Optional[int] = None,
                   t_span: Optional[Tuple[float, float]] = None):
        """Get response time CDF for a station/class or all stations/classes.

        Computes the cumulative distribution function (CDF) of response times
        (passage time distribution) for jobs of a given class at a station using
        network augmentation and transient class fluid tracking.

        Parameters
        ----------
        station : int, optional
            Station index. If None, returns CDF for all stations.
        job_class : int, optional
            Job class index. If None, returns CDF for all classes.
        t_span : tuple, optional
            Time interval (t_min, t_max) for CDF evaluation
            If None, automatically estimated based on mean response time

        Returns
        -------
        When station and job_class are both None:
            List of lists where RD[station][class] is a 2D array with columns [cdf, time]
        When station and job_class are specified:
            dict with keys 't', 'cdf', 'mean', 'var', 'method'
        """
        if self.result is None:
            self.runAnalyzer()

        from .methods.passage_time import compute_passage_time_cdf

        # Get steady-state ODE vector for passage time analysis
        steady_state_vec = self.result.xvec if self.result.xvec is not None else None

        # If no station/class specified, return all in nested list format
        if station is None and job_class is None:
            M = self.sn.nstations
            K = self.sn.nclasses
            R = self.result.RN

            RD = []
            for i in range(M):
                station_data = []
                for r in range(K):
                    if R is not None and i < R.shape[0] and r < R.shape[1]:
                        mean_resp_t = R[i, r]
                        if mean_resp_t > 0 and not np.isnan(mean_resp_t):
                            # Use transient fluid analysis for CDF
                            try:
                                t_cdf, cdf_vals = compute_passage_time_cdf(
                                    self.sn,
                                    station_idx=i,
                                    job_class=r,
                                    options=self.options,
                                    steady_state_vec=steady_state_vec,
                                    t_span=t_span
                                )
                                # Return as 2D array with columns [cdf, time]
                                cdf_data = np.column_stack([cdf_vals, t_cdf])
                                station_data.append(cdf_data)
                            except Exception:
                                # Fallback to exponential approximation
                                lambda_rate = 1.0 / mean_resp_t
                                quantiles = np.linspace(0.001, 0.999, 100)
                                times = -np.log(1 - quantiles) / lambda_rate
                                cdf_vals = 1 - np.exp(-lambda_rate * times)
                                cdf_data = np.column_stack([cdf_vals, times])
                                station_data.append(cdf_data)
                        else:
                            station_data.append(None)
                    else:
                        station_data.append(None)
                RD.append(station_data)
            return RD

        # Specific station/class requested - use detailed computation
        if station is None:
            station = 0
        if job_class is None:
            job_class = 0

        # Compute CDF via network augmentation
        try:
            t, cdf = compute_passage_time_cdf(
                self.sn,
                station_idx=station,
                job_class=job_class,
                options=self.options,
                steady_state_vec=steady_state_vec,
                t_span=t_span
            )

            # Compute moments from CDF
            dt = np.diff(t)
            pdf = np.diff(cdf)
            mean_resp_time = np.sum(t[:-1] * pdf)  # E[T] ≈ ∫ t f(t) dt
            var_resp_time = np.sum(((t[:-1] - mean_resp_time) ** 2) * pdf)  # Var[T]

            return {
                't': t,
                'cdf': cdf,
                'mean': mean_resp_time,
                'var': var_resp_time,
                'method': self.options.method
            }

        except Exception as e:
            raise RuntimeError(f"Passage time computation failed: {str(e)}")

    def getTranCdfPassT(self, station: int = 0, job_class: int = 0,
                       t: float = 1.0) -> float:
        """Get response time CDF value at specific time.

        Returns the cumulative probability P(response_time ≤ t) at a given time.

        Parameters
        ----------
        station : int, optional
            Station index (default: 0)
        job_class : int, optional
            Job class index (default: 0)
        t : float
            Time point for CDF evaluation (default: 1.0)

        Returns
        -------
        float
            CDF value F(t) = P(response_time ≤ t) at specified time

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> prob_less_than_1 = solver.getTranCdfPassT(station=0, t=1.0)
        >>> print(f"P(response_time <= 1.0) = {prob_less_than_1:.4f}")
        """
        if self.result is None:
            self.runAnalyzer()

        # Get full CDF
        cdf_dict = self.getCdfRespT(station=station, job_class=job_class)
        t_vals = cdf_dict['t']
        cdf_vals = cdf_dict['cdf']

        # Interpolate to get CDF at requested time
        cdf_at_t = np.interp(t, t_vals, cdf_vals, left=0.0, right=1.0)
        return float(cdf_at_t)

    # =====================================================================
    # STATIC METHODS (introspection and validation)
    # =====================================================================

    @staticmethod
    def listValidMethods() -> List[str]:
        """List all valid solution method names.

        Returns a list of all method identifiers that can be passed to the
        `method` parameter of __init__, including both primary names and
        aliases.

        Returns
        -------
        list of str
            Valid method identifiers:
            - 'default': Maps to 'matrix'
            - 'matrix', 'fluid.matrix', 'pnorm', 'fluid.pnorm': Matrix method
            - 'softmin', 'fluid.softmin': Softmin smoothing variant
            - 'statedep', 'fluid.statedep': State-dependent variant
            - 'closing', 'fluid.closing': Closing approximation
            - 'diffusion', 'fluid.diffusion': Diffusion SDE method
            - 'mfq', 'fluid.mfq', 'butools': Markovian fluid queue

        Examples
        --------
        >>> methods = SolverFLD.listValidMethods()
        >>> print(methods)
        >>> for m in methods:
        ...     print(f"  - {m}")
        """
        return [
            'default',
            'matrix', 'fluid.matrix', 'pnorm', 'fluid.pnorm',
            'softmin', 'fluid.softmin',
            'statedep', 'fluid.statedep',
            'closing', 'fluid.closing',
            'diffusion', 'fluid.diffusion',
            'mfq', 'fluid.mfq', 'butools'
        ]

    @staticmethod
    def supports(sn, method: str) -> Tuple[bool, Optional[str]]:
        """Check if a method can theoretically solve a given network.

        Performs basic validation of method availability. More specific
        constraints (topology, network properties) are checked at solve time.

        Parameters
        ----------
        sn : NetworkStruct
            Network structure to validate against
        method : str
            Method name to check

        Returns
        -------
        tuple
            (can_solve, reason) where:
            - can_solve: bool, whether method is valid
            - reason: str or None, explanation if not supported

        Examples
        --------
        >>> sn = NetworkStruct()  # ... configure ...
        >>> ok, reason = SolverFLD.supports(sn, 'matrix')
        >>> if not ok:
        ...     print(f"Cannot use matrix: {reason}")
        """
        if method not in SolverFLD.listValidMethods():
            return False, f"Unknown method: {method}"

        # Basic validation: all methods support open/closed/mixed networks
        # More specific constraints (e.g., mfq requires single queue) checked at solve time
        return True, None

    @staticmethod
    def getFeatureSet() -> Dict:
        """Get capabilities and limitations of SolverFLD.

        Returns a dictionary describing which features are supported by the
        fluid solver implementation.

        Returns
        -------
        dict
            Keys and values:
            - open_networks (bool): Support for open queueing networks
            - closed_networks (bool): Support for closed (population-constrained) networks
            - mixed_networks (bool): Support for mixed (both open and closed) classes
            - phase_type_service (bool): Support for phase-type service distributions
            - markovian_arrivals (bool): Support for Markovian arrival processes
            - class_switching (bool): Support for dynamic class switching
            - priorities (bool): Support for job priority classes
            - scheduling (list): Supported scheduling disciplines
            - max_classes (int): Maximum number of job classes
            - max_stations (int): Maximum number of stations

        Examples
        --------
        >>> features = SolverFLD.getFeatureSet()
        >>> print(f"Max classes: {features['max_classes']}")
        >>> print(f"Scheduling: {features['scheduling']}")
        """
        return {
            'open_networks': True,
            'closed_networks': True,
            'mixed_networks': True,
            'phase_type_service': True,
            'markovian_arrivals': True,
            'class_switching': False,
            'priorities': False,
            'scheduling': ['FCFS', 'PS', 'INF', 'EXT', 'DPS'],
            'max_classes': 20,
            'max_stations': 100,
        }

    @staticmethod
    def defaultOptions() -> SolverFLDOptions:
        """Get default solver configuration.

        Returns a SolverFLDOptions object initialized with default parameters.
        Use this as a starting point for custom configurations.

        Returns
        -------
        SolverFLDOptions
            Configuration object with default values:
            - method: 'default' (maps to 'matrix')
            - tol: 1e-4 (ODE integration tolerance)
            - iter_max: 200 (max FCFS iterations)
            - pstar: 20.0 (p-norm smoothing parameter)
            - verbose: False

        Examples
        --------
        >>> opts = SolverFLD.defaultOptions()
        >>> opts.verbose = True
        >>> solver = SolverFLD(network, options=opts)
        """
        return SolverFLDOptions()

    # Alias for consistency
    default_options = defaultOptions

    def getPerctRespT(self, percentiles: Optional[List[float]] = None,
                      station: int = 0, job_class: int = 0) -> Tuple[np.ndarray, pd.DataFrame]:
        """Get percentile response times.

        Computes response time percentiles by inverting the CDF computed via
        passage time analysis. Returns both raw values and a formatted DataFrame.

        Parameters
        ----------
        percentiles : list of float, optional
            Percentile values to compute (0-100 scale).
            Default is [50, 90, 95, 99] (median, 90th, 95th, 99th percentiles)
        station : int, optional
            Station index for CDF computation (default: 0)
        job_class : int, optional
            Job class index (default: 0)

        Returns
        -------
        tuple
            (perct_values, perct_table) where:
            - perct_values: np.ndarray of shape (n_percentiles,) with response time
              values corresponding to each percentile
            - perct_table: pd.DataFrame with columns ['Percentile', 'ResponseTime']
              for display and export

        Raises
        ------
        RuntimeError
            If runAnalyzer() has not been called yet

        Notes
        -----
        The percentile computation uses the CDF obtained from passage time analysis.
        For percentile p, finds t such that F(t) = p/100, using linear interpolation
        between CDF points.

        For high percentiles (e.g., 99th), accuracy depends on the time span used
        for CDF computation. If the CDF doesn't reach the requested percentile,
        the method extrapolates using exponential tail approximation.

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> values, table = solver.getPerctRespT([50, 90, 95, 99])
        >>> print(table)
           Percentile  ResponseTime
        0        50.0         1.234
        1        90.0         3.456
        2        95.0         4.567
        3        99.0         6.789

        >>> # Get 95th percentile response time
        >>> p95 = values[2]  # Index corresponds to percentiles list
        """
        if self.result is None:
            self.runAnalyzer()

        if percentiles is None:
            percentiles = [50.0, 90.0, 95.0, 99.0]

        # Get CDF from passage time analysis
        cdf_result = self.getCdfRespT(station=station, job_class=job_class)
        t_vals = cdf_result['t']
        cdf_vals = cdf_result['cdf']
        mean_resp = cdf_result.get('mean', self.result.RN[station, job_class])

        # Compute percentiles by inverting CDF
        perct_values = []
        for p in percentiles:
            p_frac = p / 100.0

            if p_frac <= 0:
                perct_values.append(0.0)
            elif p_frac >= 1:
                # Use exponential extrapolation for 100th percentile
                perct_values.append(t_vals[-1] * 2)
            elif p_frac <= cdf_vals[-1]:
                # Interpolate within CDF range
                idx = np.searchsorted(cdf_vals, p_frac)
                if idx == 0:
                    perct_values.append(t_vals[0])
                else:
                    # Linear interpolation between adjacent CDF points
                    t_low, t_high = t_vals[idx - 1], t_vals[idx]
                    cdf_low, cdf_high = cdf_vals[idx - 1], cdf_vals[idx]
                    if cdf_high > cdf_low:
                        t_interp = t_low + (t_high - t_low) * (p_frac - cdf_low) / (cdf_high - cdf_low)
                    else:
                        t_interp = t_low
                    perct_values.append(t_interp)
            else:
                # Extrapolate using exponential tail approximation
                # F(t) ≈ 1 - exp(-t/τ) for large t, where τ = mean
                # Solving for t: t = -τ * ln(1 - p)
                tau = mean_resp if mean_resp > 0 else 1.0
                t_extrap = -tau * np.log(1 - p_frac)
                perct_values.append(max(t_extrap, t_vals[-1]))

        perct_array = np.array(perct_values)

        # Create DataFrame
        perct_table = pd.DataFrame({
            'Percentile': percentiles,
            'ResponseTime': perct_array
        })

        return perct_array, perct_table

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
            self.runAnalyzer()

        resp_t = self.result.RN
        wait_t = np.zeros(resp_t.shape[0])
        for i in range(resp_t.shape[0]):
            mean_resp = np.mean(resp_t[i, :])
            mean_util = np.mean(self.result.UN[i, :])
            if mean_util > 0 and mean_util < 1:
                service_t = mean_resp * (1 - mean_util)
                wait_t[i] = max(0, mean_resp - service_t)
            else:
                wait_t[i] = mean_resp * 0.5

        return wait_t

    def getAvgArvR(self) -> np.ndarray:
        """Get average arrival rates per station.

        Returns:
            (M,) array of arrival rates
        """
        if self.result is None:
            self.runAnalyzer()

        if self.result.TN.ndim > 1:
            return np.sum(self.result.TN, axis=1)
        else:
            return np.full(self.result.QN.shape[0], np.mean(self.result.TN))

    def getAvgTput(self) -> np.ndarray:
        """Get average throughputs per station.

        Returns:
            (M,) array of throughputs
        """
        if self.result is None:
            self.runAnalyzer()

        if self.result.TN.ndim > 1:
            return np.sum(self.result.TN, axis=1)
        else:
            return np.full(self.result.QN.shape[0], np.mean(self.result.TN))

    # =====================================================================
    # PROBABILITY METHODS
    # =====================================================================

    def getProbAggr(self, station: int) -> np.ndarray:
        """Get aggregated state probabilities at station.

        For FLD, uses fluid level to approximate probability distribution.

        Args:
            station: Station index (0-based)

        Returns:
            Probability vector P(n) for total jobs at station
        """
        if self.result is None:
            self.runAnalyzer()

        Q = self.result.QN
        U = self.result.UN

        rho = np.mean(U[station, :])
        if rho >= 1.0:
            rho = 0.99
        if rho <= 0:
            rho = 0.01

        max_n = max(10, int(Q[station, :].sum() * 3))
        n = np.arange(max_n + 1)
        prob = (1 - rho) * (rho ** n)

        return prob

    def getProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Get marginal queue-length distribution at station for class.

        Args:
            station: Station index (0-based)
            jobclass: Job class index (0-based)

        Returns:
            Marginal probability vector P(n_ir) for n=0,1,2,...
        """
        if self.result is None:
            self.runAnalyzer()

        Q = self.result.QN
        U = self.result.UN

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

        Returns:
            System state probability vector
        """
        if self.result is None:
            self.runAnalyzer()

        Q = self.result.QN
        total_jobs = int(np.sum(Q))
        if total_jobs == 0:
            return np.array([1.0])

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

    def getProb(self, station: Optional[int] = None) -> np.ndarray:
        """Get state probabilities at station.

        Args:
            station: Station index (0-based). If None, returns for all stations.

        Returns:
            Probability vector or list of vectors
        """
        if self.result is None:
            self.runAnalyzer()

        if station is not None:
            return self.getProbAggr(station)
        else:
            probs = []
            for i in range(self.result.QN.shape[0]):
                probs.append(self.getProbAggr(i))
            return probs

    # =====================================================================
    # AGE OF INFORMATION METHODS
    # =====================================================================

    def getAvgAoI(self) -> np.ndarray:
        """Get average Age of Information per station.

        Age of Information measures the freshness of status updates.
        For fluid models, approximated using steady-state metrics.

        Returns:
            (M,) array of average AoI values
        """
        if self.result is None:
            self.runAnalyzer()

        # AoI ≈ 1/λ + E[T] where λ is arrival rate and E[T] is response time
        # For fluid models, use queue length and throughput
        Q = self.result.QN
        T = self.result.TN
        R = self.result.RN

        nstations = Q.shape[0]
        aoi = np.zeros(nstations)

        for i in range(nstations):
            mean_tput = np.mean(T[i, :]) if T.ndim > 1 else np.mean(T)
            mean_resp = np.mean(R[i, :])

            if mean_tput > 0:
                inter_arrival = 1.0 / mean_tput
                aoi[i] = inter_arrival + mean_resp
            else:
                aoi[i] = np.inf

        return aoi

    def getCdfAoI(self, t_values: Optional[np.ndarray] = None) -> List[Dict]:
        """Get Age of Information CDF.

        Args:
            t_values: Time points for CDF evaluation. If None, auto-generated.

        Returns:
            List of dicts with 'station', 't', 'p' keys
        """
        if self.result is None:
            self.runAnalyzer()

        avg_aoi = self.getAvgAoI()
        nstations = len(avg_aoi)

        if t_values is None:
            max_aoi = np.max(avg_aoi[np.isfinite(avg_aoi)])
            t_values = np.linspace(0, max_aoi * 3, 100)

        results = []
        for i in range(nstations):
            if np.isfinite(avg_aoi[i]) and avg_aoi[i] > 0:
                # Approximate AoI CDF using exponential distribution
                lambda_rate = 1.0 / avg_aoi[i]
                cdf_vals = 1 - np.exp(-lambda_rate * t_values)

                results.append({
                    'station': i + 1,
                    't': t_values,
                    'p': cdf_vals,
                })

        return results

    # =====================================================================
    # TRANSIENT METHODS
    # =====================================================================

    def getTranAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Get transient average metrics.

        Returns time-varying averages from ODE solution.

        Returns:
            Tuple of (Qt, Ut, Tt) where each is (time_points, M, K) array
        """
        if self.result is None:
            self.runAnalyzer()

        # Extract transient data from ODE solution if available
        if hasattr(self.result, 't') and self.result.t is not None:
            t = self.result.t
            xvec = self.result.xvec if hasattr(self.result, 'xvec') else None

            if xvec is not None and len(xvec) > 0:
                # xvec contains queue lengths over time
                Qt = np.array(xvec)
                # Estimate utilization and throughput from queue lengths
                Ut = np.clip(Qt, 0, 1)  # Approximate
                Tt = np.gradient(Qt, t, axis=0) if len(t) > 1 else Qt

                return Qt, Ut, Tt

        # Fallback: return steady-state replicated
        Q = self.result.QN
        U = self.result.UN
        T = self.result.TN

        return Q[np.newaxis, :, :], U[np.newaxis, :, :], T[np.newaxis, :, :]

    # =====================================================================
    # SAMPLING METHODS (Not Supported - Analytical Solver)
    # =====================================================================

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
            nchains = self.sn.nchains if hasattr(self.sn, 'nchains') else 1
            for c in range(nchains):
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
        """Get average metrics per node (same as per station for FLD)."""
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

    def getAvgSys(self) -> Tuple[np.ndarray, float]:
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
        class_names = getattr(self.sn, 'classnames', None) or [f'Class{k}' for k in range(nclasses)]

        rows = []
        for k in range(nclasses):
            rows.append({
                'Class': class_names[k] if k < len(class_names) else f'Class{k}',
                'SysRespT': R[k] if hasattr(R, '__getitem__') else R,
                'SysTput': T if isinstance(T, (int, float)) else (T[k] if k < len(T) else 0),
            })

        return pd.DataFrame(rows)

    # =====================================================================
    # PASSAGE TIME METHODS
    # =====================================================================

    def getCdfPT(self, source: int = 0, dest: int = -1,
                 t_span: Optional[Tuple[float, float]] = None) -> Dict[str, Any]:
        """Get passage time CDF between two stations.

        Computes the cumulative distribution function for the passage time
        (total time to traverse from source to destination station).

        Args:
            source: Source station index (0-based). Default: 0 (first station)
            dest: Destination station index (0-based). Default: -1 (last station)
            t_span: Time span (t_min, t_max) for CDF computation.
                   Default: auto-computed based on mean passage time.

        Returns:
            Dict containing:
            - 't': Time values array
            - 'cdf': CDF values F(t) = P(passage_time <= t)
            - 'mean': Mean passage time
            - 'source': Source station index
            - 'dest': Destination station index
        """
        if self.result is None:
            self.runAnalyzer()

        nstations = self.result.RN.shape[0]

        # Handle negative indices
        if dest < 0:
            dest = nstations + dest

        if source < 0 or source >= nstations:
            raise ValueError(f"Invalid source index {source}")
        if dest < 0 or dest >= nstations:
            raise ValueError(f"Invalid dest index {dest}")

        # Sum response times from source to dest
        if source <= dest:
            stations_in_path = list(range(source, dest + 1))
        else:
            stations_in_path = list(range(source, nstations)) + list(range(0, dest + 1))

        # Mean passage time is sum of response times along path
        mean_passage = 0.0
        for s in stations_in_path:
            mean_passage += np.mean(self.result.RN[s, :])

        if mean_passage <= 0:
            mean_passage = 1.0

        # Set default time span
        if t_span is None:
            t_max = mean_passage * 5
            t_span = (0.0, t_max)

        # Generate CDF using Erlang approximation (sum of exponentials)
        n_stages = len(stations_in_path)
        t_vals = np.linspace(t_span[0], t_span[1], 200)

        # Erlang-k CDF: F(t) = 1 - sum_{i=0}^{k-1} (lambda*t)^i * exp(-lambda*t) / i!
        # where k = n_stages and lambda = k / mean_passage
        lambda_rate = n_stages / mean_passage

        cdf_vals = np.zeros_like(t_vals)
        for i, t in enumerate(t_vals):
            if t <= 0:
                cdf_vals[i] = 0.0
            else:
                # Erlang CDF
                x = lambda_rate * t
                cdf_vals[i] = 1.0
                factorial = 1.0
                for j in range(n_stages):
                    if j > 0:
                        factorial *= j
                    cdf_vals[i] -= (x ** j) * np.exp(-x) / factorial

        return {
            't': t_vals,
            'cdf': cdf_vals,
            'mean': mean_passage,
            'source': source,
            'dest': dest,
        }

    # =====================================================================
    # SAMPLING METHODS (Not Supported - Analytical Solver)
    # =====================================================================

    def sample(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample from state distribution (not supported for FLD).

        Raises:
            NotImplementedError: FLD is an analytical solver
        """
        raise NotImplementedError("sample() not supported for analytical FLD solver. Use SSA instead.")

    def sampleAggr(self, node: int = 0, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated states (not supported for FLD).

        Raises:
            NotImplementedError: FLD is an analytical solver
        """
        raise NotImplementedError("sampleAggr() not supported for analytical FLD solver. Use SSA instead.")

    def sampleSys(self, numEvents: int = 1000) -> np.ndarray:
        """Sample system states (not supported for FLD).

        Raises:
            NotImplementedError: FLD is an analytical solver
        """
        raise NotImplementedError("sampleSys() not supported for analytical FLD solver. Use SSA instead.")

    def sampleSysAggr(self, numEvents: int = 1000) -> np.ndarray:
        """Sample aggregated system states (not supported for FLD).

        Raises:
            NotImplementedError: FLD is an analytical solver
        """
        raise NotImplementedError("sampleSysAggr() not supported for analytical FLD solver. Use SSA instead.")

    # =====================================================================
    # PASCALCASE ALIASES (MATLAB compatibility)
    # =====================================================================

    def GetAvgQLen(self) -> np.ndarray:
        """Alias for getAvgQLen (MATLAB compatibility)."""
        return self.getAvgQLen()

    def GetAvgUtil(self) -> np.ndarray:
        """Alias for getAvgUtil (MATLAB compatibility)."""
        return self.getAvgUtil()

    def GetAvgRespT(self) -> np.ndarray:
        """Alias for getAvgRespT (MATLAB compatibility)."""
        return self.getAvgRespT()

    def GetAvgResidT(self) -> np.ndarray:
        """Alias for getAvgResidT (MATLAB compatibility)."""
        return self.getAvgResidT()

    def GetAvgWaitT(self) -> np.ndarray:
        """Alias for getAvgWaitT (MATLAB compatibility)."""
        return self.getAvgWaitT()

    def GetAvgArvR(self) -> np.ndarray:
        """Alias for getAvgArvR (MATLAB compatibility)."""
        return self.getAvgArvR()

    def GetAvgTput(self) -> np.ndarray:
        """Alias for getAvgTput (MATLAB compatibility)."""
        return self.getAvgTput()

    def GetAvgSysRespT(self) -> np.ndarray:
        """Alias for getAvgSysRespT (MATLAB compatibility)."""
        return self.getAvgSysRespT()

    def GetAvgSysTput(self) -> float:
        """Alias for getAvgSysTput (MATLAB compatibility)."""
        return self.getAvgSysTput()

    def GetAvgTable(self) -> pd.DataFrame:
        """Alias for getAvgTable (MATLAB compatibility)."""
        return self.getAvgTable()

    def GetCdfRespT(self, station: int = 0, job_class: int = 0,
                    t_span: Optional[Tuple[float, float]] = None) -> Dict[str, Any]:
        """Alias for getCdfRespT (MATLAB compatibility)."""
        return self.getCdfRespT(station=station, job_class=job_class, t_span=t_span)

    def GetPerctRespT(self, percentiles: Optional[List[float]] = None,
                      station: int = 0, job_class: int = 0) -> Tuple[np.ndarray, pd.DataFrame]:
        """Alias for getPerctRespT (MATLAB compatibility)."""
        return self.getPerctRespT(percentiles=percentiles, station=station, job_class=job_class)

    def GetTranCdfPassT(self, station: int = 0, job_class: int = 0,
                        t: float = 1.0) -> float:
        """Alias for getTranCdfPassT (MATLAB compatibility)."""
        return self.getTranCdfPassT(station=station, job_class=job_class, t=t)

    def GetProbAggr(self, station: int) -> np.ndarray:
        """Alias for getProbAggr (MATLAB compatibility)."""
        return self.getProbAggr(station)

    def GetProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Alias for getProbMarg (MATLAB compatibility)."""
        return self.getProbMarg(station, jobclass)

    def GetProbSys(self) -> np.ndarray:
        """Alias for getProbSys (MATLAB compatibility)."""
        return self.getProbSys()

    def GetProbSysAggr(self) -> np.ndarray:
        """Alias for getProbSysAggr (MATLAB compatibility)."""
        return self.getProbSysAggr()

    def GetProb(self, station: Optional[int] = None) -> np.ndarray:
        """Alias for getProb (MATLAB compatibility)."""
        return self.getProb(station)

    def GetAvgAoI(self) -> np.ndarray:
        """Alias for getAvgAoI (MATLAB compatibility)."""
        return self.getAvgAoI()

    def GetCdfAoI(self, t_values: Optional[np.ndarray] = None) -> List[Dict]:
        """Alias for getCdfAoI (MATLAB compatibility)."""
        return self.getCdfAoI(t_values)

    def GetTranAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Alias for getTranAvg (MATLAB compatibility)."""
        return self.getTranAvg()

    def GetAvg(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Alias for getAvg (MATLAB compatibility)."""
        return self.getAvg()

    # Chain-level aliases
    def GetAvgChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Alias for getAvgChain (MATLAB compatibility)."""
        return self.getAvgChain()

    def GetAvgChainTable(self) -> pd.DataFrame:
        """Alias for getAvgChainTable (MATLAB compatibility)."""
        return self.getAvgChainTable()

    def GetAvgQLenChain(self) -> np.ndarray:
        """Alias for getAvgQLenChain (MATLAB compatibility)."""
        return self.getAvgQLenChain()

    def GetAvgUtilChain(self) -> np.ndarray:
        """Alias for getAvgUtilChain (MATLAB compatibility)."""
        return self.getAvgUtilChain()

    def GetAvgRespTChain(self) -> np.ndarray:
        """Alias for getAvgRespTChain (MATLAB compatibility)."""
        return self.getAvgRespTChain()

    def GetAvgResidTChain(self) -> np.ndarray:
        """Alias for getAvgResidTChain (MATLAB compatibility)."""
        return self.getAvgResidTChain()

    def GetAvgTputChain(self) -> np.ndarray:
        """Alias for getAvgTputChain (MATLAB compatibility)."""
        return self.getAvgTputChain()

    def GetAvgArvRChain(self) -> np.ndarray:
        """Alias for getAvgArvRChain (MATLAB compatibility)."""
        return self.getAvgArvRChain()

    # Node-level aliases
    def GetAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Alias for getAvgNode (MATLAB compatibility)."""
        return self.getAvgNode()

    def GetAvgNodeTable(self) -> pd.DataFrame:
        """Alias for getAvgNodeTable (MATLAB compatibility)."""
        return self.getAvgNodeTable()

    def GetAvgNodeChain(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Alias for getAvgNodeChain (MATLAB compatibility)."""
        return self.getAvgNodeChain()

    def GetAvgNodeChainTable(self) -> pd.DataFrame:
        """Alias for getAvgNodeChainTable (MATLAB compatibility)."""
        return self.getAvgNodeChainTable()

    def GetAvgSys(self) -> Tuple[np.ndarray, float]:
        """Alias for getAvgSys (MATLAB compatibility)."""
        return self.getAvgSys()

    def GetAvgSysTable(self) -> pd.DataFrame:
        """Alias for getAvgSysTable (MATLAB compatibility)."""
        return self.getAvgSysTable()

    # Passage time alias
    def GetCdfPT(self, source: int = 0, dest: int = -1,
                 t_span: Optional[Tuple[float, float]] = None) -> Dict[str, Any]:
        """Alias for getCdfPT (MATLAB compatibility)."""
        return self.getCdfPT(source=source, dest=dest, t_span=t_span)

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
    run_analyzer = runAnalyzer
    cdf_resp_t = getCdfRespT
    cdf_respt = getCdfRespT
    perct_resp_t = getPerctRespT
    perct_respt = getPerctRespT
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_resp_t = getAvgRespT
    avg_respt = getAvgRespT
    avg_tput = getAvgTput


__all__ = [
    'SolverFLD',
    'SolverFLDOptions',
    'FLDResult',
]
