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

import copy
import numpy as np
import pandas as pd
import time
from typing import Optional, Dict, List, Tuple, Any, Union
from dataclasses import dataclass
from scipy import linalg

from .options import SolverFLDOptions, FLDResult
from .utils import extract_metrics_from_handler_result, compute_response_times, compute_cycle_times, compute_system_throughput
from ...api.sn.transforms import sn_get_residt_from_respt
from ...api.sn import NodeType
from ..base import NetworkSolver, method_label, method_type
from ..fork_join_driver import ForkJoinDriverMixin
from ...constants import GlobalConstants, VerboseLevel
from ...indexed_table import IndexedTable

# Import base SolverOptions for type checking
try:
    from line_solver.solvers import SolverOptions as BaseSolverOptions
except ImportError:
    BaseSolverOptions = None


def _accost_is_linear(Rcost, h):
    """True when every per-(user,item) access graph is the standard linear chain
    (miss -> list 1, hit in list l -> list l+1, self-loop on the top list)."""
    if Rcost is None:
        return True
    lin = np.zeros((h + 1, h + 1))
    lin[0, 1] = 1.0
    for a in range(1, h):
        lin[a, a + 1] = 1.0
    lin[h, h] = 1.0
    for row in Rcost:
        if row is None:
            continue
        for g in row:
            if g is None:
                continue
            if not np.allclose(np.asarray(g, dtype=float), lin, atol=1e-9):
                return False
    return True


class SolverFLD(ForkJoinDriverMixin, NetworkSolver):
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
    Basic usage::

        >>> solver = SolverFLD(network, method='mfq')
        >>> solver.runAnalyzer()
        >>> QN = solver.result.QN  # Access raw results
        >>> qlen = solver.getAvgQLen()  # Access aggregated metrics

    Method comparison::

        >>> results = {}
        >>> for method in ['matrix', 'mfq']:
        ...     s = SolverFLD(network, method=method)
        ...     s.runAnalyzer()
        ...     results[method] = s.result

    Custom configuration::

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
        'statedep': 'matrix',  # Use matrix method for state-dependent (closing is incomplete)
        'fluid.statedep': 'matrix',
        'closing': 'closing',
        'fluid.closing': 'closing',
        'minnormal': 'minnormal',
        'fluid.minnormal': 'minnormal',
        # the refined mean field runs through the same solver, which reads the
        # requested method name to decide whether to add the O(1/N) term
        'refined': 'minnormal',
        'fluid.refined': 'minnormal',
        'tbi': 'tbi',
        'fluid.tbi': 'tbi',
        'diffusion': 'diffusion',
        'fluid.diffusion': 'diffusion',
        'mfq': 'mfq',
        'fluid.mfq': 'mfq',
        'butools': 'mfq',
        # 'butools' names the backend the MFQ branch calls and 'aoi' its
        # age-of-information reading; both are ALIASES of 'mfq', not methods of
        # their own, which is how MATLAB, the JAR and C++ spell them. Mapping
        # 'aoi' to its own route bypassed _solve_mfq's aoi_is_aoi topology test
        # and made SolverFLD(model,'aoi') refuse every model that is not a
        # bufferless or single-buffer queue -- while the same model answered
        # under 'mfq'. _solve_mfq still reaches the AoI solver, and attaches
        # aoiResults, exactly when the topology qualifies.
        'aoi': 'mfq',
        'fluid.aoi': 'mfq',
        'rmf': 'rmf',
        'fluid.rmf': 'rmf',
        'kp': 'kp',
        'fluid.kp': 'kp',
        'dae': 'dae',
        'fluid.dae': 'dae',
        # The single-station fluid limits.
        'ggisgi.fluid': 'ggisgi.fluid',
        'fluid.ggisgi': 'ggisgi.fluid',
        # the SHORT spellings too, as the C++ fluid_qsys_canonical maps them
        'ggisgi': 'ggisgi.fluid',
        'ggingi.tga': 'ggingi.tga',
        'fluid.tga': 'ggingi.tga',
        'tga': 'ggingi.tga',
        'tvms': 'tvms',
        'fluid.tvms': 'tvms',
        'mtginf': 'mtginf',
        'fluid.mtginf': 'mtginf',
        'mol': 'mol',
        'fluid.mol': 'mol',
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
        Using with NetworkStruct directly::

            >>> from line_solver.api.sn import NetworkStruct
            >>> sn = NetworkStruct()  # ... configure ...
            >>> solver = SolverFLD(sn, method='mfq')

        Using with Network object::

            >>> model = Network('TestModel')
            >>> # ... configure network ...
            >>> solver = SolverFLD(model, method='matrix')

        With custom options::

            >>> opts = SolverFLDOptions(tol=1e-6, pstar=50)
            >>> solver = SolverFLD(model, options=opts)
        """
        self.network = network
        # The shared feature gate in NetworkSolver reads self.model (as CTMC,
        # MVA and NC do). Without this alias getattr(self, 'model', None) is
        # None and supportsModelMethod returns True unconditionally, so no
        # model is ever gated against getFeatureSet().
        self.model = network

        # An auxiliary solver passed as second argument requests a warm
        # start: its steady-state solution decides the initial state of the
        # ODE integration (see NetworkSolver.initFromSolver).
        self._init_solver_arg = None
        if method_or_options is not None and not isinstance(method_or_options, (str, dict)) \
                and hasattr(method_or_options, 'getAvgQLen'):
            self._init_solver_arg = method_or_options
            method_or_options = 'default'

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
            if method_or_options == 'default' and method_from_kwargs is not None:
                method = method_from_kwargs
            else:
                method = method_or_options
        elif isinstance(method_or_options, SolverFLDOptions):
            options = method_or_options
            method = options.method
        elif (BaseSolverOptions is not None and isinstance(method_or_options, BaseSolverOptions)) \
                or hasattr(method_or_options, 'method'):
            # Convert base SolverOptions to SolverFLDOptions. The duck-typed arm
            # catches options that are not BaseSolverOptions, such as the
            # SolverLNOptions SolverLN hands to every layer solver; without it
            # the object itself was stored as the method.
            base_opts = method_or_options
            method = getattr(base_opts, 'method', 'default')
            opts_kwargs = {}
            for attr in ['method', 'tol', 'iter_max', 'iter_tol', 'verbose', 'timespan', 'samples', 'stiff', 'lang', 'arith']:
                if hasattr(base_opts, attr):
                    val = getattr(base_opts, attr)
                    if val is not None:
                        opts_kwargs[attr] = val
            options = SolverFLDOptions(**opts_kwargs)
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

        # The options are the SOLVER's, not the caller's. MATLAB passes a struct
        # by value, so two solvers built from one options variable are
        # independent; sharing the object here makes every per-solver write
        # (setInitialState's init_sol above all) land on all of them, which is
        # how an ensemble of stage solvers -- SolverENV builds them from one
        # factory closure -- ended up integrating every stage from the LAST
        # stage's initial state. `config` is copied too because it is mutated
        # in place (nhpp_sched, rate_sched).
        import copy as _copy
        options = _copy.copy(options)
        if isinstance(getattr(options, 'config', None), dict):
            options.config = dict(options.config)
        self.options = options
        self.result = None
        self.runtime = 0.0

        if self._init_solver_arg is not None:
            self.initFromSolver(self._init_solver_arg)

    def reset(self):
        """Reset the solver, clearing cached results and struct cache.

        Matches MATLAB behavior where reset() invalidates the cached struct
        so the solver re-reads the model state on the next analysis run.
        """
        self._clearResultStores()
        self.runtime = 0.0
        self.sn = None  # Force re-read of network struct (matching MATLAB)
        self._sn_is_toph = False

    def setInitialState(self, Q: np.ndarray):
        """Set initial state from queue length marginals.

        Args:
            Q: Queue lengths array of shape (M,) or (M, K) where M=stations, K=classes
        """
        Q = np.atleast_2d(Q)
        if Q.shape[0] == 1:
            Q = Q.T  # Convert row to column
        M, K = Q.shape

        # Build initial state vector matching ODE state dimension
        # For simple models: state is just queue lengths per (station, class)
        init_sol = Q.flatten()
        self.options.init_sol = init_sol

    def getName(self) -> str:
        """Get the name of this solver."""
        return "Fluid"

    get_name = getName

    def exportODEs(self, filename: str = '', notation: str = 'scalar') -> str:
        """Export the system of ODEs integrated by the mean-field methods of
        this solver (default/matrix, pnorm, closing, statedep, softmin) as a
        standalone LaTeX document, in a symbolic form that is both human and
        machine readable. Mirrors the MATLAB SolverFLD.exportODEs method.

        Parameters
        ----------
        filename : str, optional
            Path of the .tex file to write; empty returns the source only.
        notation : str, optional
            'scalar' (default) for one expanded ODE per state variable, or
            'matrix' for the compact matrix notation (``dx/dt = W'*theta(x) +
            lambda`` for the matrix/pnorm methods, ``dx/dt = J*r(x)`` for the
            closing/statedep/softmin methods).

        Returns
        -------
        str
            LaTeX source of the exported ODE system.
        """
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import export_odes_via_cpp
            tex = export_odes_via_cpp(self, notation=notation)
            if filename:
                with open(filename, 'w') as f:
                    f.write(tex)
            return tex
        from .symodes import solver_fluid_symodes, export_odes_latex
        if self.sn is None:
            self.sn = self._get_network_struct(self.network)
        sys = solver_fluid_symodes(self.sn, self.options)
        model_name = getattr(self.network, 'name', None) or 'model'
        hide_immediate = bool(getattr(self.options, 'hide_immediate', False))
        tex = export_odes_latex(sys, model_name, notation, hide_immediate)
        if filename:
            with open(filename, 'w') as f:
                f.write(tex)
        return tex

    export_odes = exportODEs

    def getSymbolicDrift(self, options=None):
        """Right-hand side of the mean-field ODE system as expression strings,
        one per state variable, together with the variable names they are
        written in.

        This is the input the computer algebra backend needs to produce a
        Jacobian or an equilibrium (see getJacobian), and it is the same system
        solver_fluid_symodes describes and exportODEs typesets, written out
        variable by variable instead of in matrix form.

        ONLY SMOOTH DRIFTS ARE EXPORTED. The default, matrix, closing and
        statedep methods scale rates by min(n_i, S_i), which is not
        differentiable at n_i = S_i, so their Jacobian does not exist there;
        emitting a one-sided derivative would be a silent lie exactly at the
        regime switch that matters. Use the p-norm smoothing (options.pstar,
        method matrix or pnorm) or the softmin method, whose drifts are smooth
        everywhere, and this function refuses the others by name.

        Parameters
        ----------
        options : SolverFLDOptions, optional
            Solver options; defaults to the solver's own.

        Returns
        -------
        tuple
            (rhs, vars, sys) with rhs the expression strings, vars the
            variable names x1 ... xn, and sys the structural description
            returned by solver_fluid_symodes.
        """
        from .symodes import solver_fluid_symodes, symbolic_drift, state_variables
        if options is None:
            options = self.options
        if self.sn is None:
            self.sn = self._get_network_struct(self.network)
        sys = solver_fluid_symodes(self.sn, options)
        return symbolic_drift(sys), state_variables(sys), sys

    get_symbolic_drift = getSymbolicDrift

    def getJacobian(self, options=None, equilibria: bool = False):
        """Jacobian of the mean-field ODE right-hand side, d f_i / d x_j, as a
        matrix of expression strings, computed exactly by the computer algebra
        engine.

        The Jacobian is what tells a fixed point apart from a limit cycle and
        gives the local convergence rate of the fluid approximation, neither of
        which a numerical integration reports. The equilibria, returned only
        when asked for, are the solutions of f(x) = 0; they can be empty when
        the system is beyond what the engine solves in closed form, which is a
        limitation of the solve and not an assertion that none exist.

        Only smooth drifts have a Jacobian: see getSymbolicDrift, which refuses
        the min-scaled methods by name rather than returning a one-sided
        derivative.

        The engine is the one named by options.config['symbolic']: sympy
        natively, or the line-sage-rest service when 'sage' or a URL is asked
        for, mirroring the toolbox/service split of MATLAB's SAGE.m.

        Parameters
        ----------
        options : SolverFLDOptions, optional
            Solver options; defaults to the solver's own.
        equilibria : bool, optional
            Also solve f(x) = 0.

        Returns
        -------
        tuple
            (J, rhs, vars, equilibria) with J[i][j] = d f_i / d x_j as an
            expression string written with '^' for powers, rhs the drift
            itself, vars the state variable names, and equilibria a list of
            dicts mapping variable name to expression string (None when not
            requested).
        """
        if options is None:
            options = self.options
        rhs, variables, _ = self.getSymbolicDrift(options)

        config = getattr(options, 'config', None) or {}
        backend = str(config.get('symbolic', 'auto')).strip()
        timeout = config.get('symbolic_timeout', 300)

        if backend.lower() == 'sage' or backend.lower().startswith('http'):
            from ...api.sym import resolve as _resolve_sym, require as _require_sym
            engine = _resolve_sym(backend)
            if engine is None:
                engine = _require_sym(backend)
            engine.timeout_s = timeout
            want = ['jacobian'] + (['equilibria'] if equilibria else [])
            r = engine.fluid_odes(rhs, variables, want)
            return (r['jacobian'], rhs, variables,
                    r.get('equilibria') if equilibria else None)

        import sympy
        from sympy.parsing.sympy_parser import (parse_expr, standard_transformations,
                                                convert_xor)
        # The drift is written with '^' for powers, which the service parses as
        # exponentiation; convert_xor gives sympy the same reading.
        transformations = standard_transformations + (convert_xor,)
        syms = [sympy.Symbol(v, real=True) for v in variables]
        local = dict(zip(variables, syms))
        exprs = [parse_expr(e, local_dict=local, transformations=transformations)
                 for e in rhs]
        # '^' on the way out as well, so the Jacobian is written in the same
        # notation as the drift it came from and either engine's output can be
        # read back by one rule.
        J = [[str(sympy.diff(expr, x)).replace('**', '^') for x in syms]
             for expr in exprs]
        sols = None
        if equilibria:
            raw = sympy.solve(exprs, syms, dict=True)
            sols = [dict((str(k), str(v).replace('**', '^')) for k, v in s.items())
                    for s in raw]
        return J, rhs, variables, sols

    get_jacobian = getJacobian

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

        # Priority 4: Model that is already a struct. (No wrapper bridge —
        # native solvers reject JAR-wrapper models, keeping python/ free of
        # any JAR/JVM coupling.)
        if hasattr(model, 'nclasses') and hasattr(model, 'nstations'):
            return model

        raise ValueError("Cannot extract network structure from model")

    def _dispatch_method(self, method_key):
        """Solve with the resolved method, switching off a degenerate closure.

        A non-hyperbolic fluid fixed point (balanced bottlenecks, a saturated
        multiclass station, an overloaded open station) leaves the linear noise
        approximation with no stationary covariance. It cannot be seen before
        the mean is solved, so fluid_minnormal_applicable cannot decline it and
        the moment closure raises at the Lyapunov step.

        THE LADDER HAS TWO RUNGS, AND THE FIRST ONE KEEPS THE CLOSURE. Most of
        these failures are not a property of the model at all: MinNormalSolver
        must start its alternation at sigma2 = 0, where min(n,c) has no
        derivative, so a saturated or balanced model's first-order fixed point
        lands on the kink, sits on a continuum of equilibria, and the Jacobian
        there is neutral. DaeSolver seeds the variance POSITIVE and never adopts
        sigma2 = 0 as an iterate, so the smoothed E[min(X,c)] breaks the
        degeneracy and the fixed point is isolated and hyperbolic -- it answers
        the same closure, with a covariance, where the alternation cannot.
        Dropping straight to first order instead is not merely a lost second
        moment: on a balanced two-station PS cycle at N=10 it returns [9 1]
        against the exact [5 5], because a first-order method has no reason to
        prefer one point of the continuum over another.

        The second rung is the first-order method, taken when 'dae' declines the
        model in advance (fluid_dae_applicable) or fails on the same exception,
        which is the genuinely non-hyperbolic case: an unstable open station has
        no stationary distribution to approximate under any closure. Fall back
        whether 'minnormal' was RESOLVED from 'default' or REQUESTED outright:
        the closure has no stationary covariance either way, so refusing an
        explicit request would only deny the caller the mean still available.

        Mirrors MATLAB @SolverFLD/runAnalyzer and the JAR SolverFluid.
        """
        # A STOCHASTIC PETRI NET HAS ONE FLUID ROUTE, and it is 'dae'. Every
        # other method builds its drift from the station/class/phase encoding,
        # where an ordinary Place declares no service process and therefore
        # contributes NO coordinate at all: the net would be integrated as an
        # empty model and the table would report zeros with no warning.
        # getMethodFeatureSet states the same limit so the gate refuses it one
        # step earlier; this names the alternative.
        if self._is_petri_net() and method_key not in ('dae', 'fluid.dae', 'default', 'fluid.default'):
            raise ValueError(
                "This model is a stochastic Petri net, which the '%s' method has no drift for. Use "
                "options.method='dae' (the default for a Petri net), or SolverCTMC, SolverJMT, "
                "SolverSSA or SolverLDES." % method_key)
        from .methods.minnormal import FluidNonHyperbolicError
        try:
            return self._solve_with_method(method_key)
        except FluidNonHyperbolicError as err:
            if method_key != 'minnormal':
                raise
            debug = GlobalConstants.getVerbose() == VerboseLevel.DEBUG
            from .dae_applicable import fluid_dae_applicable
            dae_ok, dae_reason = fluid_dae_applicable(self.sn, self.options)
            if dae_ok:
                try:
                    result = self._solve_with_method('dae')
                    if debug:
                        print("SolverFLD: minnormal declined at the Lyapunov step (%s); "
                              "falling back to dae" % (err,))
                    self._fallback_method = 'dae'
                    return result
                except FluidNonHyperbolicError as dae_err:
                    if debug:
                        print("SolverFLD: dae also declined at the Lyapunov step (%s)"
                              % (dae_err,))
            elif debug:
                print("SolverFLD: dae not applicable as a fallback (%s)" % (dae_reason,))
            fallback = 'closing' if self._has_dps_scheduling() else 'matrix'
            if debug:
                print("SolverFLD: moment closure declined at the Lyapunov step (%s); "
                      "falling back to %s" % (err, fallback))
            self._fallback_method = fallback
            return self._solve_with_method(fallback)

    def _solve_with_method(self, method_key):
        """Dispatch to the appropriate solver method and return the result."""
        if method_key == 'matrix':
            return self._solve_matrix()
        elif method_key == 'closing':
            return self._solve_closing()
        elif method_key == 'minnormal':
            return self._solve_minnormal()
        elif method_key == 'tbi':
            return self._solve_tbi()
        elif method_key == 'diffusion':
            return self._solve_diffusion()
        elif method_key == 'mfq':
            return self._solve_mfq()
        elif method_key == 'aoi':
            return self._solve_aoi()
        elif method_key == 'rmf':
            return self._solve_rmf()
        elif method_key == 'kp':
            return self._solve_kp()
        elif method_key == 'dae':
            return self._solve_dae()
        elif method_key in ('ggisgi.fluid', 'ggingi.tga', 'tvms', 'mtginf', 'mol'):
            return self._solve_qsys(method_key)
        else:
            raise ValueError(f"Unknown method: {method_key}")

    def _build_init_sol_from_raw_states(self, raw_state_per_isf, sn):
        """Build FLD ODE initial state vector from raw per-node states.

        Matches MATLAB solver_fluid_initsol.m logic: extracts per-phase
        service counts from FCFS states and distributes jobs across phases.

        Args:
            raw_state_per_isf: dict mapping stateful index -> raw state array
            sn: NetworkStruct

        Returns:
            init_sol numpy array for the ODE, or None if phases unavailable
        """
        from ...api.sn import SchedStrategy

        M = sn.nstations
        K = sn.nclasses

        # Compute phases per (station, class) from proc
        phases = np.ones((M, K), dtype=int)
        if hasattr(sn, 'proc') and sn.proc is not None:
            for i in range(M):
                for r in range(K):
                    proc_ir = None
                    if isinstance(sn.proc, dict):
                        if i in sn.proc and r in sn.proc[i]:
                            proc_ir = sn.proc[i][r]
                    elif isinstance(sn.proc, list) and i < len(sn.proc):
                        if sn.proc[i] is not None and r < len(sn.proc[i]):
                            proc_ir = sn.proc[i][r]
                    if proc_ir is not None:
                        if isinstance(proc_ir, dict):
                            if 'k' in proc_ir:
                                phases[i, r] = int(proc_ir['k'])
                        elif isinstance(proc_ir, (list, tuple)) and len(proc_ir) >= 2:
                            D0 = proc_ir[0]
                            if hasattr(D0, 'shape'):
                                phases[i, r] = D0.shape[0]

        init_sol = []
        sched_dict = sn.sched if sn.sched else {}

        for ist in range(M):
            isf = int(sn.stationToStateful[ist])
            raw_state = raw_state_per_isf.get(isf)
            sched = sched_dict.get(ist)

            # Check if Source station - skip (no mass)
            node_idx = int(sn.stationToNode[ist]) if ist < len(sn.stationToNode) else ist
            is_source = (node_idx < len(sn.nodetype)
                         and sn.nodetype[node_idx] == NodeType.SOURCE)
            if sched == SchedStrategy.EXT:
                is_source = True
            if is_source:
                for k in range(K):
                    for _ in range(int(phases[ist, k])):
                        init_sol.append(0.0)
                continue

            if raw_state is None:
                for k in range(K):
                    for _ in range(int(phases[ist, k])):
                        init_sol.append(0.0)
                continue

            raw_state = np.asarray(raw_state, dtype=float).flatten()
            phasesz = phases[ist, :]
            total_srv = int(np.sum(phasesz))

            is_fcfs = (sched in (SchedStrategy.FCFS, SchedStrategy.HOL,
                                 SchedStrategy.LCFS, SchedStrategy.LCFSPR)
                       if sched is not None else False)

            if is_fcfs and len(raw_state) > total_srv:
                # FCFS state: [buffer..., service_phase_counts...]
                buf_width = len(raw_state) - total_srv
                space_buf = raw_state[:buf_width]
                space_srv = raw_state[buf_width:]

                # Extract kir from service phase counts
                kir = {}
                offset = 0
                for r in range(K):
                    for k in range(int(phasesz[r])):
                        kir[(r, k)] = float(space_srv[offset])
                        offset += 1

                # Compute nir: total jobs of class r at station
                nir = np.zeros(K)
                for r in range(K):
                    sir_r = sum(kir[(r, k)] for k in range(int(phasesz[r])))
                    buf_count = float(np.sum(space_buf == (r + 1)))
                    nir[r] = sir_r + buf_count

                # Build init_sol entries (matching MATLAB solver_fluid_initsol)
                for r in range(K):
                    nph = int(phasesz[r])
                    if nph == 0:
                        continue
                    # Phase 1: nir - sum(kir[r,k] for k>0)
                    later = sum(kir.get((r, k), 0) for k in range(1, nph))
                    init_sol.append(nir[r] - later)
                    # Phase k > 1: kir[r,k]
                    for k in range(1, nph):
                        init_sol.append(kir.get((r, k), 0))
            else:
                # INF/PS/SIRO: state is marginal counts, all in phase 1
                for k in range(K):
                    nph = int(phasesz[k])
                    if nph == 0:
                        continue
                    n_jobs = float(raw_state[k]) if k < len(raw_state) else 0.0
                    init_sol.append(n_jobs)
                    for _ in range(1, nph):
                        init_sol.append(0.0)

        return np.array(init_sol, dtype=float)

    def _get_fld_state_spaces_and_priors(self):
        """Get per-node state spaces and priors for pprod iteration.

        Returns:
            Tuple of (per_node_spaces, per_node_priors, stateful_nodes, cur_states)
            where stateful_nodes is list of (node_index, isf) tuples and
            cur_states holds original states for restoration.
        """
        per_node_spaces = []
        per_node_priors = []
        stateful_nodes = []
        cur_states = []

        nodeToStateful = np.asarray(self.sn.nodeToStateful).flatten() \
            if hasattr(self.sn, 'nodeToStateful') else np.array([])
        K = self.sn.nclasses

        for ind in range(self.sn.nnodes):
            if hasattr(self.sn, 'isstateful') and self.sn.isstateful[ind]:
                isf = int(nodeToStateful[ind])
                stateful_nodes.append((ind, isf))

                node = self.network._nodes[ind]

                # Save current state for restoration
                cur_state = node.get_state() if hasattr(node, 'get_state') else None
                if cur_state is not None:
                    cur_states.append(np.array(cur_state).copy())
                else:
                    cur_states.append(None)

                # Get per-node state space
                node_space = getattr(node, '_state_space', None)
                if node_space is not None and len(node_space) > 0:
                    node_space = np.atleast_2d(node_space)
                else:
                    node_state = node.get_state() if hasattr(node, 'get_state') else None
                    if node_state is not None:
                        node_space = np.atleast_2d(np.asarray(node_state).flatten())
                    else:
                        node_space = np.zeros((1, K))
                per_node_spaces.append(node_space)

                # Get per-node state prior
                node_prior = getattr(node, '_state_prior', None)
                if node_prior is not None:
                    node_prior = np.asarray(node_prior).flatten()
                    if len(node_prior) < node_space.shape[0]:
                        padded = np.zeros(node_space.shape[0])
                        padded[:len(node_prior)] = node_prior
                        node_prior = padded
                else:
                    node_prior = np.zeros(node_space.shape[0])
                    node_prior[0] = 1.0
                per_node_priors.append(node_prior)

        return per_node_spaces, per_node_priors, stateful_nodes, cur_states

    def supportsTransientAnalysis(self):
        """Transient averages are available (fluid ODE integrated over options.timespan)."""
        return True

    supports_transient_analysis = supportsTransientAnalysis

    def runAnalyzer(self) -> 'SolverFLD':
        """Execute the fluid analysis using the configured method.

        Supports state prior iteration (pprod loop): when the model has
        multiple possible initial states weighted by priors, runs the
        analysis for each state and accumulates weighted results.
        Matches MATLAB SolverFLD/runAnalyzer.m lines 112-208.

        Returns
        -------
        SolverFLD
            Returns self to enable method chaining and fluent interface
        """
        # Opt-in delegation to the canonical JAR (mirrors MATLAB options.lang='java').
        # Populates the native result container from jline.jar so every getter
        # (tables, matrices, chain/node/scalar metrics) returns JAR-derived values.
        # Imported lazily so a JVM-free install never touches this path.
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import populate_java_result
            populate_java_result(self)
            return self
        # see _kb/06-solver-catalog.md ("Python lang='cpp' opt-in C++ delegation");
        # an absent binary is the only automatic fallback, a C++ refusal propagates.
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import LineCliNotAvailable, populate_cpp_result
            try:
                populate_cpp_result(self)
                return self
            except LineCliNotAvailable as e:
                import warnings
                warnings.warn("SolverFLD: lang='cpp' requested but the C++ solver is "
                              "unavailable (%s); falling back to lang='python'." % e)

        # Re-read network struct if invalidated by reset() (matching MATLAB)
        if self.sn is None:
            self.sn = self._get_network_struct(self.network)

        # Non-Markovian renewal distributions (Det/Gamma/Weibull/Pareto/Uniform/
        # Lognormal) carry a single nominal phase in the node layer; the fluid ODE
        # needs the acyclic-PH expansion. MATLAB does this at
        # FLD/@SolverFLD/runAnalyzer.m:29 and the JAR at SolverFluid.java:1077.
        # Without it the Python fluid path solved a Det service as exponential.
        self.sn = self._toph(self.sn)
        self._sn_is_toph = True

        # MATLAB FLD/@SolverFLD/runAnalyzer.m calls NetworkSolver.runAnalyzerChecks
        # here: reject models using features outside the fluid feature set rather
        # than integrating a mis-specified model and returning plausible numbers.
        # It runs AFTER the phase-type conversion and therefore after 'default'
        # can be resolved, because resolveMethod inspects the phase-resolved
        # state and the gate must validate the method that will actually run.
        model = getattr(self, 'model', None)
        if model is not None and hasattr(model, 'get_used_lang_features'):
            self.runAnalyzerChecks(self.options)

        # THE TWO GATES BELOW READ THE RESOLVED METHOD, not the requested one:
        # 'default' stands for 'dae' on exactly the models they refuse (see
        # _resolve_default_method), so gating the literal name would refuse the
        # run that is about to succeed.
        resolved_method = self._resolve_method()

        # Finite Capacity Region: the fluid ODEs do not enforce the aggregate
        # per-region job limit and would silently return the unconstrained answer.
        #
        # 'dae' is the exception, and the only one. A region cap is a linear
        # inequality on the state and blocking is a throttle on the admission flow
        # that keeps it satisfied, so the DAE form has somewhere to put it -- an
        # algebraic equation beside the drift -- where an ODE has not.
        # capacity_constraints refuses the forms that are NOT constraints on this
        # drift, by name. Every other method keeps the blanket refusal, because for
        # them it is still true.
        if getattr(self.sn, 'nregions', 0) > 0 and \
                resolved_method not in ('dae', 'fluid.dae'):
            raise RuntimeError('This model uses a Finite Capacity Region (addRegion), which is '
                               'not supported by SolverFLD (the region\'s aggregate job limit is '
                               'not enforced). Use options.method=\'dae\', or SolverCTMC, '
                               'SolverJMT, SolverSSA or SolverLDES.')

        # A BINDING STATION BUFFER WAS SILENTLY IGNORED, by every fluid method
        # including this one: nothing in the fluid tree reads sn.cap or
        # sn.classcap, so a capped station was integrated as an unbounded one and
        # the table reported more jobs in the buffer than the buffer holds (a
        # closed Delay->Queue(cap 2) model returned 2.19 jobs in a buffer of 2).
        # MVA and NC have refused that model through the shared structural gate
        # since they gained one; the fluid solver now does too, except on the
        # route that can enforce it.
        # 'mol' is exempt for the opposite reason to 'dae': a finite capacity is
        # not something it ignores, it is the model. The approximation is stated
        # for the Mt/G/s/0 LOSS system, so the server count IS the buffer, and
        # refusing a capped station would refuse the only shape the method
        # answers. The other single-station limits assume an unbounded waiting
        # room and keep the guard.
        if resolved_method not in ('dae', 'fluid.dae', 'mol', 'fluid.mol'):
            model = getattr(self, 'model', None)
            if model is not None and hasattr(model, 'get_used_lang_features'):
                ok, reason = NetworkSolver.checkBindingCapacity(model, 'SolverFLD')
                if not ok:
                    raise RuntimeError(
                        "%s Use options.method='dae', which carries the buffer as an algebraic "
                        "constraint on the drift." % reason)

        # Fork-join: the same solver-agnostic fixed point MVA and NC drive
        # (ForkJoinDriverMixin), with a fluid inner solve. The MMT transformation
        # emits only Source, Delay, Queue, Router and ClassSwitch, all of which
        # the fluid drift already carries, so no fork-join code is added here.
        # Intercepted after the feature gate so an unsupported feature is still
        # named by its own message rather than by a failure inside the transform.
        if self._has_fork_join() and not getattr(self, '_skip_fork_join', False):
            # Not every method can run that fixed point. forkJoinAdmits is the
            # same predicate supportsModelMethod asks, so the report and the run
            # cannot disagree about which forks this method serves; without it
            # an open fork-join model surfaced the DAE form's missing unknowns
            # as an UnboundLocalError from inside the transform.
            fj_ok, fj_reason = SolverFLD.forkJoinAdmits(self.sn, self.options.method)
            if not fj_ok:
                raise RuntimeError(fj_reason)
            fj_t0 = time.time()
            fj_result = self._run_fork_join_analysis()
            if fj_result is not None:
                self.runtime = time.time() - fj_t0
                if self.result is not None:
                    self.result.runtime = self.runtime
                return self

        method_key = self._resolve_method()
        # set by _dispatch_method when a resolved 'minnormal' is switched off a
        # degenerate fixed point, so the banner and result.method name the
        # method that actually produced the numbers
        self._fallback_method = None

        if GlobalConstants.getVerbose() == VerboseLevel.DEBUG:
            print(f"SolverFLD: Using method '{method_key}'")

        start_time = time.time()

        # Get per-node state spaces and priors for pprod loop
        per_node_spaces, per_node_priors, stateful_nodes, cur_states = \
            self._get_fld_state_spaces_and_priors()
        sizes = [s.shape[0] for s in per_node_spaces]
        n_nodes = len(stateful_nodes)

        total_combinations = 1
        for s in sizes:
            total_combinations *= s

        if total_combinations <= 1:
            # Single state - no pprod loop needed
            self.result = self._dispatch_method(method_key)
        else:
            # pprod loop over all initial state combinations
            # (matching MATLAB runAnalyzer.m lines 122-208)
            M = self.sn.nstations
            K = self.sn.nclasses
            Q_accum = np.zeros((M, K))
            U_accum = np.zeros((M, K))
            R_accum = np.zeros((M, K))
            T_accum = np.zeros((M, K))
            C_accum = np.zeros((1, K))
            X_accum = np.zeros((1, K))
            Qt_accum = {}
            Ut_accum = {}
            Tt_accum = {}
            t_accum = None
            first_result = True
            last_xvec = None
            total_iter = 0

            s0_id = [0] * n_nodes
            while True:
                # Compute joint prior probability
                s0prior_val = 1.0
                for i in range(n_nodes):
                    s0prior_val *= per_node_priors[i][s0_id[i]]

                if s0prior_val > 0:
                    # Set node states on model
                    for i, (ind, isf) in enumerate(stateful_nodes):
                        self.network._nodes[ind].setState(
                            per_node_spaces[i][s0_id[i]])

                    # Build raw state map for init_sol computation
                    raw_state_per_isf = {}
                    for i, (ind, isf) in enumerate(stateful_nodes):
                        raw_state_per_isf[isf] = per_node_spaces[i][s0_id[i]]

                    # Get fresh sn
                    self.network._has_struct = False
                    self.network._sn = None
                    self.sn = self._get_network_struct(self.network)

                    # Build init_sol from raw states (with per-phase info)
                    # Matches MATLAB solver_fluid_initsol.m
                    init_sol = self._build_init_sol_from_raw_states(
                        raw_state_per_isf, self.sn)
                    saved_init_sol = self.options.init_sol
                    self.options.init_sol = init_sol

                    # Run solver
                    result = self._dispatch_method(method_key)
                    self.options.init_sol = saved_init_sol

                    if result is not None:
                        last_xvec = result.xvec
                        total_iter = max(total_iter, result.iterations)

                        if first_result:
                            Q_accum = result.QN * s0prior_val
                            U_accum = result.UN * s0prior_val
                            R_accum = result.RN * s0prior_val
                            T_accum = result.TN * s0prior_val
                            C_accum = result.CN * s0prior_val
                            X_accum = result.XN * s0prior_val
                            t_accum = result.t

                            if result.QNt:
                                for key, val in result.QNt.items():
                                    Qt_accum[key] = val * s0prior_val
                            if result.UNt:
                                for key, val in result.UNt.items():
                                    Ut_accum[key] = val * s0prior_val
                            if result.TNt:
                                for key, val in result.TNt.items():
                                    Tt_accum[key] = val * s0prior_val
                            first_result = False
                        else:
                            Q_accum += result.QN * s0prior_val
                            U_accum += result.UN * s0prior_val
                            R_accum += result.RN * s0prior_val
                            T_accum += result.TN * s0prior_val
                            C_accum += result.CN * s0prior_val
                            X_accum += result.XN * s0prior_val

                            # Accumulate transient with interpolation
                            if result.QNt and result.t is not None and t_accum is not None:
                                tunion = np.union1d(t_accum, result.t)
                                for key in result.QNt:
                                    if key in Qt_accum:
                                        old_data = np.interp(tunion, t_accum, Qt_accum[key])
                                        new_data = np.interp(tunion, result.t, result.QNt[key])
                                        Qt_accum[key] = old_data + s0prior_val * new_data
                                    else:
                                        Qt_accum[key] = result.QNt[key] * s0prior_val
                                for key in (result.UNt or {}):
                                    if key in Ut_accum:
                                        old_data = np.interp(tunion, t_accum, Ut_accum[key])
                                        new_data = np.interp(tunion, result.t, result.UNt[key])
                                        Ut_accum[key] = old_data + s0prior_val * new_data
                                    else:
                                        Ut_accum[key] = result.UNt[key] * s0prior_val
                                for key in (result.TNt or {}):
                                    if key in Tt_accum:
                                        old_data = np.interp(tunion, t_accum, Tt_accum[key])
                                        new_data = np.interp(tunion, result.t, result.TNt[key])
                                        Tt_accum[key] = old_data + s0prior_val * new_data
                                    else:
                                        Tt_accum[key] = result.TNt[key] * s0prior_val
                                t_accum = tunion

                # Advance pprod
                carry = True
                for i in range(n_nodes - 1, -1, -1):
                    if carry:
                        s0_id[i] += 1
                        if s0_id[i] >= sizes[i]:
                            s0_id[i] = 0
                        else:
                            carry = False
                            break
                if carry:
                    break

            # Restore original states
            for i, (ind, isf) in enumerate(stateful_nodes):
                if cur_states[i] is not None:
                    self.network._nodes[ind].setState(cur_states[i])
            self.network._has_struct = False
            self.network._sn = None
            self.sn = self._get_network_struct(self.network)

            # Build accumulated result
            from line_solver.api.sn.getters import sn_get_arvr_from_tput
            from line_solver.api.sn.transforms import sn_get_residt_from_respt
            AN = sn_get_arvr_from_tput(self.sn, T_accum)
            WN = sn_get_residt_from_respt(self.sn, R_accum, None)

            self.result = FLDResult(
                QN=Q_accum, UN=U_accum, RN=R_accum, TN=T_accum,
                CN=C_accum, XN=X_accum, AN=AN, WN=WN,
                t=t_accum, QNt=Qt_accum, UNt=Ut_accum, TNt=Tt_accum,
                xvec=last_xvec, iterations=total_iter,
                runtime=0.0, method=method_key
            )

        # Name the method that actually produced the numbers, prefixed
        # 'default/' when the request was 'default', as MATLAB
        # @SolverFLD/runAnalyzer and the JAR SolverFluid both report it.
        if self.result is not None:
            self.result.method = method_label(
                self.options.method, self._fallback_method or method_key)

        self.runtime = time.time() - start_time

        # For cache models solved with rmf, update hit/miss probs on model nodes
        # Matches MATLAB runAnalyzer.m lines 239-251
        if method_key == 'rmf' and self.result is not None:
            cache_hit = getattr(self.result, '_cacheHitProb', None)
            cache_miss = getattr(self.result, '_cacheMissProb', None)
            if cache_hit is not None and cache_miss is not None:
                from ...api.sn.network_struct import NodeType
                sn = self.sn
                cache_nodes = [ind for ind in range(sn.nnodes)
                               if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CLASSSWITCH
                               and sn.nodeparam and ind in sn.nodeparam
                               and hasattr(sn.nodeparam[ind], 'nitems') and sn.nodeparam[ind].nitems > 0]
                if hasattr(self, 'network') and hasattr(self.network, '_nodes'):
                    for cIdx, ind in enumerate(cache_nodes):
                        if cIdx < cache_hit.shape[0] and ind < len(self.network._nodes):
                            node = self.network._nodes[ind]
                            if hasattr(node, 'set_result_hit_prob'):
                                node.set_result_hit_prob(cache_hit[cIdx, :])
                            if hasattr(node, 'set_result_miss_prob'):
                                node.set_result_miss_prob(cache_miss[cIdx, :])

        # Print completion message (matches MATLAB verbose guard)
        if self.options.verbose:
            import sys as _sys
            py_version = f"{_sys.version_info.major}.{_sys.version_info.minor}.{_sys.version_info.micro}"
            ran = self._fallback_method or method_key
            from line_solver.solvers.base import print_solver_banner
            print_solver_banner(f"Fluid analysis [method: {method_label(self.options.method, ran)}; type: {method_type('FLD', method_label(self.options.method, ran))}; lang: python; env: {py_version}] completed in {self.runtime:.6f}s.")

        return self

    def _toph(self, sn):
        """Acyclic-PH expansion of the non-Markovian renewal distributions.

        sn_nonmarkov_toph reads options as a mapping (options.get('config')),
        so the options object is converted the same way SolverCTMC does at
        solver_ctmc.py:515. The fluid ODEs read mu*phi as a flow, so the
        surrogate must be a genuine phase-type: a matrix exponential has no
        such reading, hence phfit='ph'.
        """
        from ...api.sn.transforms import sn_nonmarkov_toph
        options_dict = dict(vars(self.options)) if hasattr(self.options, '__dict__') else {'config': {}}
        cfg = dict(options_dict.get('config') or {})
        cfg['phfit'] = 'ph'
        options_dict['config'] = cfg
        return sn_nonmarkov_toph(sn, options_dict)

    def _resolve_default_method(self):
        """Concrete method that method='default' stands for.

        Preference order: 'rmf' for cache models, then 'minnormal' wherever
        fluid_minnormal_applicable accepts the model, then the historical
        choice of 'closing' for DPS and 'matrix' otherwise. The second-order
        closure dominates the first-order methods on every family measured
        against exact CTMC and is the only method that can represent GPS at
        all, so it is preferred wherever it applies. Mirrors the MATLAB
        fluid_resolve_default_method and the JAR SolverFluid.

        Returns:
            (method, reason) with reason naming why 'minnormal' was declined,
            empty when it was selected
        """
        if self._has_cache_nodes():
            return 'rmf', 'the model has cache nodes'
        # A BINDING BUFFER OR A CAPACITY REGION ALSO HAS ONE FLUID ROUTE, for
        # the same reason a Petri net does: nothing else in the fluid tree reads
        # sn.cap, sn.classcap or the region limit, so every other method
        # integrates the capped station as an unbounded one -- which is exactly
        # why runAnalyzer refuses them. Resolving 'default' to one of those
        # turned a model this solver CAN answer into an error whose advice was
        # to type the very method the resolution should have picked. The test is
        # the gate's own, so the two cannot disagree; where the dae route
        # declines, falling through leaves the refusal to the gate, which names
        # the blocking feature.
        if self._blocked_resolves_to_dae():
            if getattr(self.sn, 'nregions', 0) > 0:
                return 'dae', 'the model has a finite capacity region'
            return 'dae', 'the model has a binding finite buffer'
        # The applicability test counts PHASES, so it must see the same struct
        # the analyzer will integrate: a Det or Gamma service carries one
        # nominal phase before sn_nonmarkov_toph and its acyclic-PH expansion
        # after it. Resolving on the raw struct would let the gate validate one
        # method and the dispatch run another.
        sn = self.sn if getattr(self, '_sn_is_toph', False) else self._toph(self.sn)
        from .minnormal_applicable import fluid_minnormal_applicable
        ok, reason = fluid_minnormal_applicable(sn, self.options)
        if ok:
            return 'minnormal', ''
        if self._has_dps_scheduling():
            return 'closing', reason
        return 'matrix', reason

    def _blocked_resolves_to_dae(self) -> bool:
        """Does method='default' stand for 'dae' on this model?

        True exactly when a station buffer or a capacity region BINDS and the
        dae route accepts the model. The capacity test is
        NetworkSolver.checkBindingCapacity, the one the gate in runAnalyzer
        applies, so the resolution and the refusal cannot disagree. Mirrors the
        MATLAB fluid_resolve_default_method and the JAR
        SolverFluid.blockedResolvesToDae.
        """
        if self.sn is None:
            self.sn = self._get_network_struct(self.network)
        blocked = getattr(self.sn, 'nregions', 0) > 0
        if not blocked:
            model = getattr(self, 'model', None)
            if model is not None and hasattr(model, 'get_used_lang_features'):
                ok, _ = NetworkSolver.checkBindingCapacity(model, 'SolverFLD')
                blocked = not ok
        if not blocked:
            return False
        from .dae_applicable import fluid_dae_applicable
        sn = self.sn if getattr(self, '_sn_is_toph', False) else self._toph(self.sn)
        dae_ok, _ = fluid_dae_applicable(sn, self.options)
        return dae_ok

    def _resolve_method(self) -> str:
        """Resolve method name to internal key.

        Resolves 'default' through _resolve_default_method, and for an
        explicit 'matrix' request selects 'rmf' for cache models and 'closing'
        when DPS scheduling is present, as the matrix method supports neither.

        Returns:
            Internal method key
        """
        method = self.options.method
        if method in ('default', 'fluid.default'):
            # A PETRI NET RESOLVES TO 'dae', which is its only fluid route.
            if self._is_petri_net():
                return 'dae'
            resolved, reason = self._resolve_default_method()
            if GlobalConstants.getVerbose() == VerboseLevel.DEBUG:
                if resolved == 'minnormal':
                    print("SolverFLD: default method resolved to: minnormal")
                else:
                    print("SolverFLD: default resolved to %s, minnormal declined (%s)"
                          % (resolved, reason))
            return resolved
        resolved = self.METHODS.get(method, method)

        # Auto-select rmf for cache models
        if resolved == 'matrix' and self._has_cache_nodes():
            return 'rmf'

        # Check for DPS scheduling - matrix method doesn't support it
        if resolved == 'matrix' and self._has_dps_scheduling():
            return 'closing'

        return resolved

    def resolveMethod(self, options):
        """Concrete method the feature gate must validate.

        The gate in NetworkSolver.runAnalyzerChecks validates
        getMethodFeatureSet(method), and 'default' is not 'minnormal', so
        without this override a GPS model would be rejected before the
        resolution ever ran. Keeping the decision in one place is what stops
        the gate and the dispatch from disagreeing.
        """
        requested = getattr(options, 'method', 'default')
        if requested in ('default', 'fluid.default'):
            if self.sn is None:
                self.sn = self._get_network_struct(self.network)
            if self._is_petri_net():
                return 'dae'
            return self._resolve_default_method()[0]
        return self.METHODS.get(requested, requested)

    resolve_method = resolveMethod

    def _has_cache_nodes(self) -> bool:
        """Check if network has Cache nodes.

        Returns:
            True if any node is a Cache node
        """
        from line_solver.api.sn import NodeType
        nodetype = getattr(self.sn, 'nodetype', None)
        if nodetype is None:
            return False
        for nt in nodetype:
            if nt == NodeType.CACHE:
                return True
        return False

    def _is_petri_net(self) -> bool:
        """Whether the model holds a Transition node, i.e. is a Petri net.

        The fluid Petri route is reached from the 'dae' branch only; every other
        fluid method builds its drift from the station/class/phase encoding,
        where an ordinary Place declares no service process and contributes NO
        coordinate at all, so the net would be integrated as an empty model and
        the table would report zeros with no warning.
        """
        from line_solver.api.sn import NodeType
        nodetype = getattr(self.sn, 'nodetype', None)
        if nodetype is None:
            return False
        for nt in nodetype:
            if nt == NodeType.TRANSITION:
                return True
        return False

    def _has_dps_scheduling(self) -> bool:
        """Check if network has DPS (Discriminatory Processor Sharing) scheduling.

        Returns:
            True if any station has DPS scheduling
        """
        from line_solver.api.sn import SchedStrategy

        sched_dict = getattr(self.sn, 'sched', None) or {}
        for i in range(self.sn.nstations):
            station_sched = sched_dict.get(i)
            if station_sched is None:
                continue
            # Handle both enum and integer representations
            if station_sched == SchedStrategy.DPS:
                return True
            if hasattr(station_sched, 'value') and station_sched.value == SchedStrategy.DPS.value:
                return True
            if isinstance(station_sched, int) and station_sched == SchedStrategy.DPS.value:
                return True
        return False

    def _ensure_result(self):
        """Return an available result or raise RuntimeError if analysis fails."""
        if self.result is not None:
            return self.result

        try:
            self.runAnalyzer()
        except Exception as exc:
            raise RuntimeError("runAnalyzer() must complete before accessing results") from exc

        if self.result is None:
            raise RuntimeError("runAnalyzer() must complete before accessing results")

        return self.result

    def _solve_matrix(self) -> FLDResult:
        """Solve using matrix method (existing handler implementation).

        Returns:
            FLDResult
        """
        from line_solver.api.solvers.fld.handler import solver_fld, SolverFLDOptions as HandlerFLDOptions

        # Get initial state from network marginals if set
        init_sol = self.options.init_sol
        if init_sol is None and hasattr(self.sn, 'state') and self.sn.state is not None:
            # Try to use network's current state
            try:
                state = np.array(self.sn.state).flatten()
                if len(state) > 0 and np.sum(state) > 0:
                    init_sol = state
            except:
                pass

        # Convert options to handler format
        # Only pass pstar when explicitly set (matching MATLAB default: no p-norm smoothing)
        pstar_list = [self.options.pstar] * self.sn.nstations if self.options.pstar is not None else None
        handler_opts = HandlerFLDOptions(
            method=self.options.method,
            tol=self.options.tol,
            verbose=self.options.verbose,
            stiff=self.options.stiff,
            iter_max=self.options.iter_max,
            timespan=self.options.timespan,
            pstar=pstar_list,
            num_cdf_pts=200,
            init_sol=init_sol,
            tranpoints=getattr(self.options, 'tranpoints', None)
        )

        # Solve
        handler_result = solver_fld(self.sn, handler_opts)

        # Extract metrics
        QN, UN, RN, TN, CN, XN = extract_metrics_from_handler_result(handler_result, self.sn)

        # Extract transient data from handler result
        QNt = {}
        UNt = {}
        TNt = {}
        if hasattr(handler_result, 'Qt') and handler_result.Qt is not None:
            M = len(handler_result.Qt)
            K = len(handler_result.Qt[0]) if M > 0 else 0
            for i in range(M):
                for r in range(K):
                    if handler_result.Qt[i][r] is not None:
                        QNt[(i, r)] = np.array(handler_result.Qt[i][r])
                    if hasattr(handler_result, 'Ut') and handler_result.Ut is not None:
                        UNt[(i, r)] = np.array(handler_result.Ut[i][r])
                    if hasattr(handler_result, 'Tt') and handler_result.Tt is not None:
                        TNt[(i, r)] = np.array(handler_result.Tt[i][r])

        # Compute proper arrival rates from throughputs using routing
        from line_solver.api.sn.getters import sn_get_arvr_from_tput
        AN = sn_get_arvr_from_tput(self.sn, TN) if TN is not None else None

        # Compute proper residence times from response times
        from line_solver.api.sn.transforms import sn_get_residt_from_respt
        WN = sn_get_residt_from_respt(self.sn, RN, None) if RN is not None else None

        # Build result
        result = FLDResult(
            QN=QN,
            UN=UN,
            RN=RN,
            TN=TN,
            CN=CN,
            XN=XN,
            AN=AN,
            WN=WN,
            t=handler_result.t,
            QNt=QNt,
            UNt=UNt,
            TNt=TNt,
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

    def _solve_minnormal(self) -> FLDResult:
        """Solve using the second-order moment closure.

        The drift uses E[min(X_i,c_i)] under a normal marginal whose variance
        comes from the covariance (Lyapunov) equation, so mean and covariance
        are solved self-consistently by fixed-point iteration. This is the only
        FLD method that can represent GPS, whose capacity share depends on the
        backlog indicator rather than on the populations.

        Returns:
            FLDResult
        """
        from .methods.minnormal import MinNormalSolver

        # A cache model is a DECOMPOSITION, not one ODE: the caches are solved
        # in isolation and the network with them relabeled as class switches.
        # The closure applies to the queueing layer of that alternation, so the
        # route is the same one 'rmf' takes, with the closure inside its
        # network step (_solve_rmf reads options.method to decide).
        if self._has_cache_nodes():
            return self._solve_rmf()

        method = MinNormalSolver(self.sn, self.options)
        return method.solve()

    def _solve_dae(self) -> FLDResult:
        """Solve the min-normal closure as one differential-algebraic system.

        Same closure as '_solve_minnormal' -- same drift, same rate factors,
        same Lyapunov equation -- stated and solved as one system instead of by
        successive substitution: population conservation becomes an EQUATION
        rather than a consequence of the drift, and the transient carries a
        time-varying covariance rather than the stationary one.

        A CACHE MODEL IS REFUSED RATHER THAN REROUTED, which is where this
        parts company with '_solve_minnormal' above. That method hands a cache
        model to the 'rmf' alternation, whose network step can carry a closure;
        there is no such route for the DAE form, because a decomposition has no
        single drift for the algebraic constraint to be attached to.

        Returns:
            FLDResult
        """
        from .methods.dae import DaeSolver

        # A STOCHASTIC PETRI NET TAKES ITS OWN ANALYZER, and returns from here.
        # The route mirrors the MATLAB twin (solver_fluid_analyzer's Transition
        # branch): a model with any Transition node is a different formalism,
        # and the queueing post-processing below would overwrite the Petri
        # conventions -- a Place is an INF station whose utilization is its
        # token count and whose throughput is the firing rate of the modes
        # consuming from it.
        if self._is_petri_net():
            from .methods.petri import PetriSolver
            return PetriSolver(self.sn, self.options).solve()

        if self._has_cache_nodes():
            raise ValueError(
                "The dae method does not support caching stations: a cache model is solved by "
                "decomposition, so it has no single drift to constrain. Use "
                "options.method='minnormal' for the same closure, or 'rmf'.")

        method = DaeSolver(self.sn, self.options)
        return method.solve()

    def _solve_tbi(self) -> FLDResult:
        """Solve using trajectory-based iteration (TBI).

        Returns:
            FLDResult

        Partitions stations into cells and applies Jacobi waveform relaxation
        over growing-horizon time segments. Closed networks only.
        """
        from .methods.tbi import TBIMethod

        method = TBIMethod(self.sn, self.options)
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
        from line_solver.lib.thirdparty.aoi import aoi_is_aoi
        from .methods.mfq import MFQMethod

        is_aoi, _ = aoi_is_aoi(self.sn)
        if is_aoi:
            result = self._solve_aoi()
            result.method = 'mfq'
            return result

        # MFQ IS A SINGLE-QUEUE METHOD AND FALLS BACK, which is what the
        # reference does: solver_fluid_analyzer.m warns "MFQ not applicable:
        # ... Falling back to matrix method" and re-enters solver_fluid_matrix.
        # Raising instead made 'mfq' -- and therefore its aliases 'butools' and
        # 'aoi' -- refuse every multi-station model that MATLAB, the JAR and C++
        # all answer.
        try:
            method = MFQMethod(self.sn, self.options)
        except ValueError as not_applicable:
            from ...api.io.logging import line_warning
            line_warning('solver_fluid_analyzer',
                         'MFQ not applicable: %s. Falling back to matrix method.'
                         % (not_applicable,))
            # `_fallback_method` is what runAnalyzer reports through
            # method_label; setting result.method directly is overwritten there,
            # and the label has to say `matrix` as MATLAB and the JAR both do.
            self._fallback_method = 'matrix'
            result = self._solve_matrix()
            result.method = 'matrix'
            return result
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

    def _solve_kp(self) -> FLDResult:
        """Fluid + diffusion limits of Ko and Pender (2017).

        Integrates the mean and the covariance of the (MAP_t/Ph_t/inf)^N limit
        jointly, so the result carries a variance trajectory that no other FLD
        method provides.
        """
        from .methods.kp import solve_kp
        return solve_kp(self.sn, self.options)

    def _solve_qsys(self, method_key) -> FLDResult:
        """One of the single-station fluid limits.

        These are closed forms, not integrations of the network drift: they take
        the whole model in one call and have no initial state to average over.
        solver_fluid_qsys_analyzer refuses any model that is not the
        Source -> Queue -> Sink shape they are stated for.
        """
        import time as _time
        from ...api.solvers.fld.qsys import solver_fluid_qsys_analyzer
        t0 = _time.time()
        opts = copy.copy(self.options)
        opts.method = method_key
        out = solver_fluid_qsys_analyzer(self.sn, opts)
        res = FLDResult(QN=out['QN'], UN=out['UN'], RN=out['RN'], TN=out['TN'],
                        CN=out['CN'].reshape(1, -1), XN=out['XN'].reshape(1, -1),
                        AN=out['AN'])
        res.method = out['method']
        res.iterations = 1
        res.runtime = _time.time() - t0
        if out['Qt']:
            traj = out['Qt'][0][0]
            res.t = np.asarray(traj)[:, 1]
            for i in range(len(out['Qt'])):
                for r in range(len(out['Qt'][i])):
                    res.QNt[(i, r)] = np.asarray(out['Qt'][i][r])[:, 0]
                    res.UNt[(i, r)] = np.asarray(out['Ut'][i][r])[:, 0]
                    res.TNt[(i, r)] = np.asarray(out['Tt'][i][r])[:, 0]
        return res

    def _solve_rmf(self) -> FLDResult:
        """Solve cache+queueing network using fixed-point iteration.

        Iterates between:
          1. Isolated cache analysis (cache_miss_rmf / cache_gamma_lp; RANDOM(m)/RR only)
          2. Fluid ODE solution of the surrounding queueing network
        until arrival rates to the cache converge.

        Matches MATLAB solver_fld_cacheqn_analyzer.m.
        """
        sn = self.sn
        I = sn.nnodes
        K = sn.nclasses
        M = sn.nstations

        # THE CACHE REWRITE IS UNDONE ON THE WAY OUT, exception or not. The
        # alternation below RELABELS each cache node as a class switch in
        # `self.sn` itself, and the struct is what every later gate reads: leave
        # it rewritten and `fluid_dae_applicable`, which declines a cache model
        # by name, stops seeing one. The fallback ladder in `_dispatch_method`
        # then admits `dae` for a model that has no dae route, and `dae` answers
        # the REWRITTEN network -- measured on the cacheqn model as a total
        # population of 4.2797 against N = 4, returned with no error at all.
        # The JAR saves and restores `sn.rt` around the same rewrite for the
        # same reason; see _kb/06-solver-catalog.md.
        nodetype_orig = list(sn.nodetype) if sn.nodetype is not None else None
        try:
            return self._solve_rmf_inner(sn, I, K, M)
        finally:
            if nodetype_orig is not None:
                for _i, _nt in enumerate(nodetype_orig):
                    sn.nodetype[_i] = _nt

    def _solve_rmf_inner(self, sn, I, K, M) -> FLDResult:
        """The cache/queueing alternation itself; see :meth:`_solve_rmf`."""
        from ...api.sn.network_struct import NodeType
        from ...api.cache import (cache_miss_rmf, cache_miss_sfifo_rmf,
                                   cache_miss_fifo_rmf, cache_gamma_lp)
        from ...lang.base import ReplacementStrategy
        from ...api.mc.dtmc import dtmc_stochcomp
        from ...api.sn.transforms import sn_refresh_visits

        # Build statefulNodesClasses list (matching MATLAB)
        stateful_nodes_classes = []
        for ind in range(I):
            if sn.isstateful is not None and sn.isstateful[ind]:
                for k in range(K):
                    stateful_nodes_classes.append(ind * K + k)
        stateful_nodes_classes = np.array(stateful_nodes_classes, dtype=int)

        lambda_arr = np.zeros(K)
        lambda_1 = np.zeros(K)

        # Find cache nodes
        cache_indices = []
        for ind in range(I):
            if ind < len(sn.nodetype) and sn.nodetype[ind] == NodeType.CACHE:
                cache_indices.append(ind)

        if not cache_indices:
            # No cache nodes - fall back to matrix method
            return self._solve_matrix()

        hitprob = np.zeros((len(cache_indices), K))
        missprob = np.zeros((len(cache_indices), K))
        use_moments = str(getattr(self.options, 'method', '')).endswith('minnormal')
        last_moments = None

        iter_max = self.options.iter_max if self.options.iter_max else 100
        iter_tol = self.options.iter_tol if hasattr(self.options, 'iter_tol') and self.options.iter_tol else 1e-6

        # Converged per-cache isolated inputs (gamma, m, lambda_cache, strat),
        # captured for the transient path (_cacheqn_tran). Overwritten each
        # iteration so the final entries hold the converged values.
        cache_tran_inputs = [None] * len(cache_indices)

        for it in range(1, iter_max + 1):
            for cIdx, ind in enumerate(cache_indices):
                ch = sn.nodeparam.get(ind) if sn.nodeparam else None
                if ch is None:
                    continue

                hitclass = np.asarray(ch.hitclass, dtype=int)
                missclass = np.asarray(ch.missclass, dtype=int)
                input_classes = [r for r in range(K) if r < len(hitclass) and r < len(missclass)
                                 and hitclass[r] >= 0 and missclass[r] >= 0]

                m = np.asarray(ch.itemcap, dtype=int)
                n = ch.nitems

                if it == 1:
                    # Initial random arrival rates
                    for r in input_classes:
                        lambda_1[r] = np.random.rand()
                    lambda_arr = lambda_1.copy()
                    sn.nodetype[ind] = NodeType.CLASSSWITCH

                # Solution of isolated cache
                h = len(m)
                u = K  # number of users = number of classes
                lambda_cache = np.zeros((u, n, h + 1))

                for v in range(u):
                    for ki in range(n):
                        for l in range(h + 1):
                            pread_v = ch.pread[v] if v < len(ch.pread) else None
                            if pread_v is not None and not (np.isscalar(pread_v) and np.isnan(pread_v)):
                                if ki < len(pread_v):
                                    lambda_cache[v, ki, l] = lambda_arr[v] * pread_v[ki]

                Rcost = getattr(ch, 'accost', None)
                if Rcost is None:
                    # Default linear cache routing (matches MVA pattern)
                    def _default_routing(h_val):
                        mat = np.diag(np.ones(h_val), 1)
                        mat[h_val, h_val] = 1.0
                        return mat
                    Rcost = [[_default_routing(h) for _ in range(n)] for _ in range(u)]

                gamma, _, _, _, _ = cache_gamma_lp(lambda_cache, Rcost)

                # Native fluid (drift-based) cache models: RANDOM(m) and FIFO(m)
                # share the RAND(m) refined mean field (Gast15 Thm 1:
                # pi_FIFO(m) = pi_RAND(m)); strict FIFO(m) uses its own
                # position-resolved mean field (cache_miss_sfifo_rmf). Every
                # other strategy (LRU/HLRU/CLIMB/QLRU) has no drift-based fluid
                # model; refuse rather than substitute a non-fluid, algebraic
                # (FPI/characteristic-time) fixed point.
                # A custom access graph (accost) modulates admission (row 0) and
                # promotion (row 1+i) per item. RANDOM(m) honours an arbitrary
                # graph via the general drift. FIFO(m)/strict FIFO(m) equal
                # RANDOM(m) only for the linear chain (Gast15 Thm 1 does NOT
                # extend to general graphs), and their position-resolved general
                # drift is not yet available, so a non-linear graph is rejected
                # for them rather than silently served the linear-graph result.
                # RR/FIFO/SFIFO honour a custom access graph via their general
                # drift; linear default keeps the refined/linear path (see _kb).
                _rs = getattr(ch, 'replacestrat', None)
                _nonlinear = not _accost_is_linear(Rcost, h)
                if _rs == ReplacementStrategy.RR:
                    _, missrate, _, _ = cache_miss_rmf(gamma, m, lambda_cache, accost=Rcost)
                elif _rs == ReplacementStrategy.FIFO:
                    if _nonlinear:
                        _, missrate, _, _ = cache_miss_fifo_rmf(gamma, m, lambda_cache, accost=Rcost)
                    else:
                        _, missrate, _, _ = cache_miss_rmf(gamma, m, lambda_cache)
                elif _rs == ReplacementStrategy.SFIFO:
                    _, missrate, _, _ = cache_miss_sfifo_rmf(gamma, m, lambda_cache, accost=Rcost)
                else:
                    raise RuntimeError(
                        "SolverFLD supports only RANDOM(m)/FIFO(m) (refined mean "
                        "field) and strict FIFO(m) (position-resolved mean field) "
                        "cache replacement; replacement strategy %s has no "
                        "drift-based fluid model. Use SolverNC/SolverMVA or "
                        "SolverLDES for this cache." % str(_rs))

                if missrate is not None:
                    for r in input_classes:
                        if lambda_arr[r] > 0:
                            missprob[cIdx, r] = missrate[r] / lambda_arr[r]
                        else:
                            missprob[cIdx, r] = 0
                        hitprob[cIdx, r] = 1.0 - missprob[cIdx, r]
                    hitprob[np.isnan(hitprob)] = 0
                    missprob[np.isnan(missprob)] = 0

                # Capture converged isolated inputs for the transient path.
                cache_tran_inputs[cIdx] = {
                    'node': ind, 'gamma': gamma, 'm': m,
                    'lambda_cache': lambda_cache, 'strat': _rs,
                    'input_classes': list(input_classes),
                }

                # Update routing matrix with hit/miss probabilities
                for r in input_classes:
                    sn.rtnodes[ind * K + r, :] = 0
                    for jnd in range(I):
                        if sn.connmatrix is not None and sn.connmatrix[ind, jnd]:
                            sn.rtnodes[ind * K + r, jnd * K + hitclass[r]] = hitprob[cIdx, r]
                            sn.rtnodes[ind * K + r, jnd * K + missclass[r]] = missprob[cIdx, r]

                sn.rt = dtmc_stochcomp(sn.rtnodes, stateful_nodes_classes)
                if hasattr(sn, 'rt_visits') and sn.rt_visits is not None:
                    sn.rt_visits = sn.rt.copy()

            # Refresh visits
            sn_refresh_visits(sn)

            # Solve the queueing network. The caches are already relabeled as
            # class switches above, so this is a plain queueing network and the
            # moment closure applies to it unchanged; 'minnormal' therefore
            # reaches a cache model through the same decomposition as 'rmf',
            # with the closure in place of the first-order matrix method.
            saved_method = self.options.method
            saved_init_sol = self.options.init_sol
            self.options.init_sol = None  # Let solver compute init_sol from current sn
            try:
                if use_moments:
                    self.options.method = 'minnormal'
                    result = self._solve_minnormal()
                    if result is not None:
                        last_moments = getattr(result, 'moments', None)
                else:
                    self.options.method = 'matrix'
                    result = self._solve_matrix()
            finally:
                # Restored on the way out however this leaves: the network step
                # can raise FluidNonHyperbolicError, and the ladder that catches
                # it reads options.method to decide what was being attempted.
                self.options.method = saved_method
                self.options.init_sol = saved_init_sol

            if result is None:
                break

            QN = result.QN
            TN = result.TN

            # Compute system throughputs
            XN = np.zeros(K)
            for k in range(K):
                refstat_k = int(sn.refstat.flat[k]) if sn.refstat is not None else 0
                if refstat_k >= 0:
                    XN[k] = TN[refstat_k, k]

            # Update arrival rates to cache using nodevisits
            nodevisits_combined = np.zeros((I, K))
            if sn.nodevisits is not None:
                for c_key, nv in sn.nodevisits.items():
                    if isinstance(nv, np.ndarray):
                        nodevisits_combined[:nv.shape[0], :nv.shape[1]] += nv

            for cIdx, ind in enumerate(cache_indices):
                ch = sn.nodeparam.get(ind) if sn.nodeparam else None
                if ch is None:
                    continue
                hitclass = np.asarray(ch.hitclass, dtype=int)
                input_classes = [r for r in range(K) if r < len(hitclass) and hitclass[r] >= 0]

                for r in input_classes:
                    c = -1
                    for ci in range(sn.nchains):
                        if ci in sn.chains and sn.chains[ci] is not None:
                            chain_classes = sn.chains[ci]
                            if isinstance(chain_classes, np.ndarray):
                                if r < len(chain_classes) and chain_classes[r]:
                                    c = ci
                                    break
                            elif r in chain_classes:
                                c = ci
                                break
                    if c < 0:
                        continue

                    inchain = []
                    if c in sn.chains:
                        chain_classes = sn.chains[c]
                        if isinstance(chain_classes, np.ndarray):
                            inchain = list(np.where(chain_classes)[0])
                        else:
                            inchain = list(chain_classes)

                    refstat_r = int(sn.refstat.flat[r]) if sn.refstat is not None else 0
                    refnode = int(sn.stationToNode[refstat_r]) if sn.stationToNode is not None else refstat_r

                    refclass_c = int(sn.refclass[c]) if hasattr(sn, 'refclass') and sn.refclass is not None and c < len(sn.refclass) else -1
                    ref_class_for_norm = refclass_c if refclass_c >= 0 else r

                    nv_denom = nodevisits_combined[refnode, ref_class_for_norm]
                    nv_num = nodevisits_combined[ind, r]

                    if nv_denom > 0:
                        lambda_arr[r] = sum(XN[k] for k in inchain) * nv_num / nv_denom

            if np.linalg.norm(lambda_arr - lambda_1, 1) < iter_tol:
                break
            lambda_1 = lambda_arr.copy()

        # The hit/miss split is a solver RESULT and belongs on the node, not only
        # on this solver's result: getHitRatio reads it there, and so does any
        # caller that solves here and reads elsewhere (the MATLAB lang='python'
        # bridge rebuilds the hit-class and miss-class node throughputs from it).
        # MATLAB's SolverFLD writes it at solve time for the same reason; without
        # this, only the sn.nodeparam patch in getAvgNode carried it, so a plain
        # getAvg left the node empty. Mirrors SolverFLD/runAnalyzer.m.
        model_nodes = self.model.getNodes()
        for cIdx, ind in enumerate(cache_indices):
            cache_node = model_nodes[ind] if ind < len(model_nodes) else None
            if cache_node is not None and hasattr(cache_node, 'set_result_hit_prob'):
                cache_node.set_result_hit_prob(hitprob[cIdx, :])
                cache_node.set_result_miss_prob(missprob[cIdx, :])

        # Store hit/miss probs in result's xvec for runAnalyzer to retrieve
        if result is not None:
            result.method = 'minnormal' if use_moments else 'rmf'
            # Attach cache hit/miss probs as extra attributes
            result._cacheHitProb = hitprob
            result._cacheMissProb = missprob
            if use_moments:
                # The moment report: the queueing fields as any other model
                # returns them, plus the cache occupancy covariance evaluated at
                # the converged isolated-cache inputs. The cache is re-solved
                # once here rather than inside the sweep because only the FINAL
                # arrival rates define the fixed point it linearises about.
                moments = dict(last_moments) if last_moments else {}
                moments['cache'] = self._cache_moments(cache_tran_inputs, K)
                result.moments = moments

        # Expose the converged isolated-cache inputs for the transient path.
        self._cache_rmf_inputs = cache_tran_inputs

        return result

    def _cache_moments(self, cache_inputs, nclasses):
        """Second moment of each cache, at the converged isolated-cache inputs.

        The linear noise approximation of the RANDOM(m) drift gives the
        stationary covariance of the item occupancy. The per-item miss
        indicator is coordinate (i, list 0), so its variance is the leading
        n-item block of the diagonal; the miss probability a class sees is the
        popularity-weighted sum of those indicators, hence a linear functional
        whose variance is w' W00 w. Only RR/FIFO on the linear access chain have
        the drift the covariance linearises, so a cache without one is omitted
        rather than reported as a fabricated zero.

        Twin of the MATLAB local `local_cache_moments` in
        `solver_fld_cacheqn_analyzer`.
        """
        from ...api.cache import cache_rmf_lna
        from ...lang.base import ReplacementStrategy

        out = []
        for entry in cache_inputs:
            if entry is None:
                continue
            if entry.get('strat') not in (ReplacementStrategy.RR, ReplacementStrategy.FIFO):
                continue
            lam = np.asarray(entry['lambda_cache'], dtype=float)
            m = np.asarray(entry['m'], dtype=float).ravel()
            u, nitems = lam.shape[0], lam.shape[1]
            h = len(m)
            dim = nitems * (h + 1)

            lam_i = np.zeros(nitems)
            for v in range(u):
                row = np.array(lam[v, :, 0], dtype=float)
                row[~np.isfinite(row)] = 0.0
                lam_i += row
            if np.sum(lam_i) <= 0:
                continue
            p = lam_i / np.sum(lam_i)

            # linearise at the SAME point the mean is reported at, i.e. the
            # refined fixed point pi + V/n when it is finite, exactly as
            # cache_miss_rmf does; the covariance of a different point is a
            # different number
            try:
                from ...api.cache.rmf import _fixed_point, _expansion_steady_state, _idx
                x0 = np.zeros(dim)
                obj_idx = 0
                for k in range(1, h + 1):
                    for _ in range(int(m[k - 1])):
                        if obj_idx < nitems:
                            x0[_idx(obj_idx, k, nitems)] = 1.0
                            obj_idx += 1
                for i in range(obj_idx, nitems):
                    x0[_idx(i, 0, nitems)] = 1.0
                x = _fixed_point(x0, p, m, nitems, h, dim)
                try:
                    pi_mf, V = _expansion_steady_state(x0, p, m, nitems, h, dim)
                    xref = pi_mf + V / nitems
                    if np.all(np.isfinite(xref)):
                        x = xref
                except Exception:
                    pass
                W = cache_rmf_lna(x, p, m, nitems, h, dim)
                pi0 = np.clip(x[:nitems], 0.0, 1.0)
            except Exception:
                continue

            W00 = W[:nitems, :nitems]
            miss_var = np.zeros(nclasses)
            for r in range(min(nclasses, u)):
                w = np.array(lam[r, :, 0], dtype=float)
                w[~np.isfinite(w)] = 0.0
                tot = np.sum(w)
                if tot <= 0:
                    continue
                w = w / tot
                miss_var[r] = max(0.0, float(w @ W00 @ w))

            out.append({'node': entry['node'], 'pi0': np.asarray(pi0).ravel(),
                        'Sigma': W, 'pi0Var': np.maximum(0.0, np.diag(W00)),
                        'missProbVar': miss_var})
        return out

    def _cacheqn_tran(self, tspan, x0cell=None):
        """Transient refined-mean-field cache trajectory.

        Port of MATLAB solver_fld_cacheqn_tran.m. Converges the per-class cache
        arrival rates via the steady RMF fixed-point iteration (_solve_rmf),
        then integrates the mean-field drift over ``tspan`` to obtain the
        time-resolved per-class hit/miss probabilities of each cache.
        RANDOM(m)/FIFO(m) share the steady state (Gast15 Thm 1) but not the
        transient, so FIFO uses its own position-resolved drift; strict FIFO(m)
        likewise. LRU/HLRU/CLIMB/QLRU have no drift-based transient.

        Returns (tcache, hitprob_t, missprob_t, caches, arate) with
        hitprob_t/missprob_t of shape (ncaches, nclasses, nt).
        """
        from ...api.cache import (cache_miss_rmf, cache_miss_fifo_rmf,
                                   cache_miss_sfifo_rmf)
        from ...lang.base import ReplacementStrategy

        # Ensure the fixed point (and the converged isolated inputs) are available.
        if getattr(self, '_cache_rmf_inputs', None) is None:
            self._solve_rmf()
        inputs = [ci for ci in getattr(self, '_cache_rmf_inputs', []) if ci is not None]

        t0, t1 = float(tspan[0]), float(tspan[-1])
        if np.isinf(t0):
            t0 = 0.0
        tsp = [t0, t1]

        K = self.sn.nclasses
        ncaches = len(inputs)
        caches = [ci['node'] for ci in inputs]
        arate = np.zeros((ncaches, K))
        tcache = None
        hitprob_t = None
        missprob_t = None
        for cIdx, ci in enumerate(inputs):
            gamma, m, lam, strat = ci['gamma'], ci['m'], ci['lambda_cache'], ci['strat']
            x0 = x0cell[cIdx] if (x0cell is not None and cIdx < len(x0cell)) else None
            if strat in (ReplacementStrategy.RR, ReplacementStrategy.FIFO):
                fn = cache_miss_rmf if strat == ReplacementStrategy.RR else cache_miss_fifo_rmf
            elif strat == ReplacementStrategy.SFIFO:
                fn = cache_miss_sfifo_rmf
            else:
                raise RuntimeError(
                    "Transient cache analysis is only available for RANDOM(m)/"
                    "FIFO(m) and strict FIFO(m) replacement via a drift-based "
                    "mean field; strategy %s has none." % str(strat))
            res = fn(gamma, m, lam, tspan=tsp, x0init=x0)
            tc, MU_t = res[4], res[6]
            if tcache is None:
                nt = len(tc)
                tcache = np.asarray(tc).ravel()
                hitprob_t = np.zeros((ncaches, K, nt))
                missprob_t = np.zeros((ncaches, K, nt))
            u = lam.shape[0]
            for v in range(u):
                rowrate = float(np.nansum(lam[v, :, 0]))
                arate[cIdx, v] = rowrate
                if rowrate > 0:
                    mp = np.clip(MU_t[v, :] / rowrate, 0.0, 1.0)
                    missprob_t[cIdx, v, :] = mp
                    hitprob_t[cIdx, v, :] = 1.0 - mp
        return tcache, hitprob_t, missprob_t, caches, arate

    # =====================================================================
    # RESULT ACCESS METHODS (following SolverMAM pattern)
    # =====================================================================

    def getAvgTable(self) -> pd.DataFrame:
        """Get average performance metrics as DataFrame.

        Returns:
            DataFrame with columns: Station, JobClass, QLen, Util, RespT, ResidT, ArvR, Tput
        """
        if self.result is None:
            self._ensureAvgResults()

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
        AN = self.result.AN if hasattr(self.result, 'AN') and self.result.AN is not None else TN

        # Compute ResidT using proper visit ratios from network structure
        # This uses the correct formula: WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]
        if self.sn is not None and self.sn.visits:
            WN = sn_get_residt_from_respt(self.sn, RN, None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            WN = RN.copy()

        rows = []
        for i in range(nstations):
            for r in range(nclasses):
                qlen = QN[i, r] if i < QN.shape[0] and r < QN.shape[1] else 0
                util = UN[i, r] if i < UN.shape[0] and r < UN.shape[1] else 0
                respt = RN[i, r] if i < RN.shape[0] and r < RN.shape[1] else 0
                residt = WN[i, r] if i < WN.shape[0] and r < WN.shape[1] else respt
                tput = TN[i, r] if i < TN.shape[0] and r < TN.shape[1] else 0
                arvr = AN[i, r] if i < AN.shape[0] and r < AN.shape[1] else tput

                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'JobClass': class_names[r] if r < len(class_names) else f'Class{r}',
                    'QLen': qlen,
                    'Util': util,
                    'RespT': respt,
                    'ResidT': residt,
                    'ArvR': arvr,
                    'Tput': tput,
                })

        df = pd.DataFrame(rows)

        # Wrap in IndexedTable for consistent formatting
        result = IndexedTable(df)

        if len(df) > 0 and not getattr(self, '_table_silent', False):
            print(result)

        return result

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
        result = self._ensure_result()
        # Queue length is additive across classes: the number of customers at a
        # station is the sum of per-class queue lengths (conserves population),
        # not their mean.
        return np.sum(result.QN, axis=1)

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
        result = self._ensure_result()
        # Per-class utilizations are additive: the fraction of servers busy at a
        # station (serving any class) is the sum of per-class utilizations, not
        # their mean.
        return np.sum(result.UN, axis=1)

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
        result = self._ensure_result()
        # Response time is per station and class (M x K, as documented); unlike
        # queue length it is not additive across classes, so return the full
        # matrix rather than collapsing it.
        return np.asarray(result.RN)

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
        result = self._ensure_result()
        return np.mean(result.TN, axis=0)

    def getAvgSysRespT(self) -> np.ndarray:
        """Get average system response time per job class.

        Returns the total time a customer spends in the system for each job class.

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
        For closed networks: uses Little's Law C = N/X
        For open networks: sum of response times across all stations

        Examples
        --------
        >>> solver = SolverFLD(network).runAnalyzer()
        >>> sys_resp = solver.getAvgSysRespT()
        >>> print(f"Mean system response time: {np.mean(sys_resp):.3f}")
        """
        result = self._ensure_result()

        RN = result.RN
        XN = result.XN.flatten() if result.XN is not None else np.zeros(RN.shape[1])
        njobs = np.asarray(self.sn.njobs).flatten() if self.sn is not None and hasattr(self.sn, 'njobs') else None
        nclasses = RN.shape[1]
        C = np.zeros(nclasses)

        for k in range(nclasses):
            if njobs is not None and k < len(njobs) and np.isfinite(njobs[k]):
                # Closed class: use Little's Law (matching MATLAB getAvgSys.m line 135)
                if XN[k] > 0:
                    C[k] = njobs[k] / XN[k]
                else:
                    C[k] = np.inf
            else:
                # Open class: sum of response times across all stations
                C[k] = np.sum(RN[:, k])

        return C

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
        result = self._ensure_result()
        return np.mean(result.XN)

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
            self._ensureAvgResults()

        from .methods.passage_time import compute_passage_time_cdf

        # Get steady-state ODE vector for passage time analysis
        steady_state_vec = self.result.xvec if self.result.xvec is not None else None
        # and the closure that solve closed its drift at, so the passage time is
        # measured on the same drift (see compute_passage_time_cdf)
        sigma2_drift = None
        # getMoments(), not result.moments: a delegated solve carries the closure
        # on the transport, not on the result object, and reading the attribute
        # left lang='java' integrating the FIRST-ORDER drift while the mean it
        # reports came from the closure.
        moments = self.getMoments()
        if isinstance(moments, dict):
            sigma2_drift = moments.get('sigma2Drift')

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
                            # Use transient fluid analysis for CDF. No fallback:
                            # an exponential substituted on failure reports the
                            # WRONG distribution (SCV 1 for every station) with
                            # nothing in the output to say so
                            t_cdf, cdf_vals = compute_passage_time_cdf(
                                self.sn,
                                station_idx=i,
                                job_class=r,
                                options=self.options,
                                steady_state_vec=steady_state_vec,
                                t_span=t_span,
                                sigma2=sigma2_drift
                            )
                            # Return as 2D array with columns [cdf, time]
                            cdf_data = np.column_stack([cdf_vals, t_cdf])
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
                t_span=t_span,
                sigma2=sigma2_drift
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
            self._ensureAvgResults()

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
            - 'minnormal', 'fluid.minnormal': Second-order moment closure
            - 'refined', 'fluid.refined': O(1/N) refined mean field (Gast)
            - 'diffusion', 'fluid.diffusion': Diffusion SDE method
            - 'mfq', 'fluid.mfq', 'butools': Markovian fluid queue
            - 'aoi', 'fluid.aoi': explicit AoI MFQ solver

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
            'minnormal', 'fluid.minnormal',
            'refined', 'fluid.refined',
            'tbi', 'fluid.tbi',
            'diffusion', 'fluid.diffusion',
            'mfq', 'fluid.mfq', 'butools',
            'rmf', 'fluid.rmf',
            'aoi', 'fluid.aoi',
            'kp', 'fluid.kp',
            'dae', 'fluid.dae',
            # The single-station fluid limits (Source -> Queue -> Sink, one
            # class). Unlike the MATLAB twin this list is static, so they are
            # named on every model and the analyzer refuses the shapes they are
            # not stated for -- the convention this solver already follows for
            # 'mfq' and the rest.
            # 'ggisgi' and 'tga' are the SHORT spellings, mapped onto the two
            # primary names as the C++ fluid_qsys_canonical does
            'ggisgi.fluid', 'fluid.ggisgi', 'ggisgi',
            'ggingi.tga', 'fluid.tga', 'tga',
            'tvms', 'fluid.tvms',
            'mtginf', 'fluid.mtginf',
            'mol', 'fluid.mol',
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
    def canonicalMethod(method):
        """The one spelling of a fluid method that every gate tests against.

        NOT `METHODS`, which is a DISPATCH map: it sends 'refined' to
        'minnormal' because they share a solver routine, while their feature
        envelopes differ ('refined' is closed-only). This collapses SPELLING
        only: the 'fluid.' qualifier, the MFQ backend aliases 'butools' and
        'aoi', and the short spellings of the two single-station limits.

        Canonicalizing once is what keeps an alias from carrying a different
        envelope than the name it resolves to, and it is the reason the four
        codebases can no longer drift apart over a spelling. MATLAB
        SolverFLD.canonicalMethod, the JAR and C++ apply the same three rules
        in the same order.
        """
        if not isinstance(method, str):
            return method
        m = method[6:] if method.startswith('fluid.') else method
        if m in ('butools', 'aoi'):
            return 'mfq'
        if m == 'ggisgi':
            return 'ggisgi.fluid'
        if m == 'tga':
            return 'ggingi.tga'
        return m

    canonical_method = canonicalMethod

    def getMethodFeatureSet(self, method):
        """Feature envelope of a method, narrowed for 'kp' and for GPS."""
        feats = SolverFLD.getFeatureSet()
        m = SolverFLD.canonicalMethod(method)
        # GPS divides the server by weight among the BACKLOGGED classes, so its
        # share is a function of the backlog INDICATOR. A first-order closure
        # cannot express it at all: with continuous x_k > 0 every class is always
        # backlogged and the share collapses to the constant w_k/sum_j w_j, the
        # heavy-traffic limit, regardless of load. Only 'minnormal' supplies the
        # P(X_k >= 1) the closure needs. Mirrors MATLAB SolverFLD and the JAR.
        if m != 'minnormal':
            feats = set(feats) - {'SchedStrategy_GPS'}
        # Limited load dependence composes with the closure as a rate multiplier
        # alpha(n_i) on the scheduling share, which only the closing family
        # evaluates (_ode_rate_factors / closures.py). The matrix, pnorm, softmin,
        # statedep, tbi, diffusion, mfq, kp and rmf paths build their drift
        # independently and would silently return the alpha == 1 answer.
        # 'dae' belongs in this list and was missing, so a load-dependent model
        # was refused on the one closing-family method that evaluates alpha(n_i)
        # as an algebraic system: 'dae' IS the min-normal closure, same drift and
        # same rate factors, solved as one system instead of by substitution.
        # MATLAB and C++ have always kept all four.
        if m not in ('closing', 'minnormal', 'refined', 'dae'):
            feats = set(feats) - {'LoadDependence'}
        # Scheduling disciplines with no branch in the closing drift. A station
        # without a case in _ode_rate_factors keeps rates = x, i.e. it is
        # integrated as an INFINITE SERVER, and the answer is wrong without any
        # warning: on Delay(Z=1) -> Queue(c=1), N=4, exact Q2 = 3.0154, the
        # fall-through returns 2.0000. SIRO is worse still, because the closing
        # metric reader accepts it AS FCFS: the ODE integrates it as INF while
        # the metrics are read as if it shared the server. This port was the only
        # one of the four with no such strip at all, so it answered all three
        # silently where MATLAB, the JAR and C++ refuse. matrix/pnorm build a PS
        # drift for every queueing station, the right aggregate for any
        # work-conserving discipline, so they are unaffected.
        if m in ('closing', 'statedep', 'softmin', 'tbi', 'minnormal', 'refined', 'dae'):
            feats = set(feats) - {'SchedStrategy_SIRO', 'SchedStrategy_LCFS',
                                  'SchedStrategy_LCFSPR'}
        # HOL allocates capacity in PRIORITY order, not in proportion to
        # population, and no fluid drift reads sn.classprio except the
        # single-queue MFQ priority branch. Declaring it for every method, as
        # this port did, offers a priority model to drifts that would answer it
        # as if the classes shared the server proportionally.
        if m != 'mfq':
            feats = set(feats) - {'SchedStrategy_HOL'}
        # A stochastic Petri net has no drift outside the DAE form: its conserved
        # quantities are P-invariants rather than chain populations, and an
        # immediate transition is an algebraic FLOW rather than an event with a
        # rate. Every other fluid method builds its drift from the
        # station/class/phase encoding, where a Place contributes no coordinate
        # at all, so it would integrate the net as an empty model and report
        # zeros without a warning.
        if m != 'dae':
            feats = set(feats) - {'Place', 'Transition', 'Enabling', 'Inhibiting',
                                  'Timing', 'Firing', 'Storage', 'Linkage'}
        if m == 'dae':
            feats = set(feats)
            # A CAPACITY LIMIT IS A LINEAR INEQUALITY ON THE STATE, which the DAE
            # form can carry as an algebraic equation beside the drift and no ODE
            # method can carry at all. The gate is where this has to be declared:
            # runAnalyzer's own refusal sits downstream of runAnalyzerChecks, so
            # without this the model is rejected as an unsupported feature before
            # the method is ever consulted. capacity_constraints still refuses the
            # region forms that are not constraints on this drift, by name.
            feats.add('Region')
            # DPS closes on the covariance BETWEEN a station's class coordinates,
            # not on the station total. 'minnormal' carries those blocks through
            # its outer iteration; the DAE has no unknown for them, since a matrix
            # block per station restores the quartic cost that keeping Sigma out of
            # the Newton vector avoids.
            feats -= {'SchedStrategy_DPS'}
        if m in ('ggisgi.fluid', 'ggingi.tga', 'tvms'):
            # The only fluid methods in LINE stated for a queue customers
            # ABANDON. Reneging stays out of the base FLD envelope: the network
            # drift carries no abandonment flow, so every other method would
            # integrate the model as if nobody left.
            feats = set(feats) | {'Reneging'}
        if m in ('ggisgi.fluid', 'ggingi.tga', 'tvms', 'mtginf', 'mol'):
            # Every one of them is stated for a single open station; the base
            # envelope's closed classes have no meaning there.
            feats = set(feats) - {'ClosedClass', 'SelfLoopingClass'}
        if m == 'refined':
            # CLOSED MODELS ONLY, which the MATLAB runAnalyzer has always
            # enforced by name and the featset never stated: the 1/N correction
            # is solved on orth(D) over the FULL state, so on an open model it
            # adds a perturbation to the SOURCE POOL mass, a normalisation
            # constant rather than a population. Only 'minnormal' was validated
            # open. Stating it here is what lets a report withdraw the pair
            # instead of offering a run that stops -- on an open fork-join model
            # the same restriction surfaced as a failure inside the MMT fixed
            # point rather than as a refusal. The C++ twin asserts the same
            # refusal (cpp/tests/test_fluid_moments.cpp).
            feats = set(feats) - {'OpenClass', 'Source', 'Sink',
                                  'RandomSource', 'JobSink'}
        if m == 'diffusion':
            # The diffusion SDE PROJECTS each class back onto its own fixed
            # population at every step, which is the closed-network constraint
            # itself: an open class has no population to project onto, and a
            # Source is not a station the SDE has a coordinate for. Mirrors
            # MATLAB SolverFLD.getMethodFeatureSet and the JAR.
            feats = set(feats) - {'OpenClass', 'Source', 'Sink',
                                  'RandomSource', 'JobSink'}
        if m in ('diffusion', 'kp'):
            # NEITHER OF THESE TWO INTEGRATES A FORK-JOIN MODEL, and each says so
            # by answering rather than by refusing, which is the reason to state
            # it here. Measured on a SYMMETRIC closed fork-join (Delay -> Fork ->
            # two identical FCFS queues -> Join, N = 2) whose exact chain is
            # Q1 = Q2 = 0.664, J = 0.624, D = 1.024: 'diffusion' returns the whole
            # population on ONE station and zero elsewhere -- a different station
            # on a rerun, so the SDE is not integrating this model at all -- and
            # 'kp' returns an ALL-ZERO table on a symmetric OPEN fork-join fed at
            # rate 0.5, an empty network where jobs are arriving. The C++ featset
            # has always withheld the names; MATLAB and this port offered them
            # and mis-answered.
            feats = set(feats) - {'Fork', 'Join', 'Forker', 'Joiner', 'JoinPartial'}
        if m == 'tbi':
            # Trajectory-based iteration decomposes the CLOSED population into
            # cells and relaxes the waveforms between them; there is no cell for
            # an unbounded open stream. A cache model is solved by decomposition
            # rather than by one drift, so the cell partition has nothing to
            # partition -- use 'rmf'.
            feats = set(feats) - {'OpenClass', 'Source', 'Sink',
                                  'RandomSource', 'JobSink',
                                  'Cache', 'CacheClassSwitcher',
                                  'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
                                  'ReplacementStrategy_SFIFO'}
        if m == 'kp':
            # The Ko-Pender limits are proved for an OPEN network of stations
            # fed by external arrival processes: a closed class has no arrival
            # process to modulate and no source phase to carry, and the cache
            # and class-switch machinery has no counterpart in the paper's
            # event set. Narrow the envelope rather than fail at solve time.
            feats = set(feats)
            feats -= {'ClosedClass', 'SelfLoopingClass', 'Cache',
                      'CacheClassSwitcher', 'ClassSwitch', 'StatelessClassSwitcher',
                      'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
                      'ReplacementStrategy_SFIFO'}
        # MULTISERVER (registry name since 2026-09-05): the drifts carry
        # min(n,c) except the diffusion SDE, which is written for one or
        # infinitely many servers. fluid_method_refusal keeps wording that
        # refusal.
        #
        # 'mfq' is NOT withdrawn here, where MATLAB's SolverFLD.m does withdraw
        # it. MATLAB can, because its gate asks about the RESOLVED method and
        # SolverFLD.resolveMethod sends 'mfq' to 'matrix' off the single-queue
        # shape fluid_mfq_admits decides; native python has no such resolution
        # for 'mfq', so the withdrawal would be a refusal rather than a
        # relabelling. It would also be untrue of this analyzer: _solve_mfq
        # returns the exact M/M/c utilisation lambda/(c*mu), which is what
        # test_solver_fld_matlab_parity::test_mmc_mfq_exact pins as MATLAB
        # parity. A feature set states what its own analyzer honours.
        if m == 'diffusion':
            feats = set(feats) - {'MultiServer'}
        return feats

    def supportsModelMethod(self, method):
        """The structural finite-capacity gate runAnalyzer enforces at solve
        time, stated here so that a CALLER can see it before running.

        Nothing in the fluid tree reads sn.cap or sn.classcap, so every method
        but two integrates a capped station as an unbounded one. 'dae' carries
        the buffer as an algebraic constraint on the drift, and 'mol' is stated
        for the Mt/G/s/0 LOSS system, where the server count IS the buffer; the
        rest keep the guard. There is no registry feature name for plain
        capacity, hence the structural test -- SolverNC and SolverMVA gate the
        same way.

        Left only in runAnalyzer the rule was invisible to every gate above it,
        and SolverAUTO.listValidMethods offered all 29 fluid methods on the
        BAS-blocking model of cqn_bas_blocking, each of which then raised when
        asked to run. Mirrors MATLAB @SolverFLD/supportsModelMethod.
        """
        ok, reason = super().supportsModelMethod(method)
        if ok:
            ok, reason = self._single_station_shape_admits(method)
            if not ok:
                return ok, reason
        if ok and method in self._TIME_VARYING_METHODS:
            # The time-varying single-station limits report a TRAJECTORY, so
            # they need a finite options.timespan. A horizon is an option and
            # not a model feature, hence the structural test; the predicate is
            # the one solver_fluid_qsys_analyzer stops on, so the report and the
            # run cannot answer differently.
            from ...api.solvers.fld.qsys import fluid_qsys_horizon
            hok, hreason, _, _ = fluid_qsys_horizon(self.options)
            if not hok:
                return False, ("The '%s' method reports a trajectory. %s" % (method, hreason))
        # A fork-join model is answered by the MMT fixed point rather than by one
        # drift, and not every method can run it. Fork and OpenClass are both
        # declared names, so the featset cannot state a rule that is their
        # CONJUNCTION; it is structural, and it is the predicate runAnalyzer
        # stops on.
        if ok:
            model = getattr(self, 'model', None)
            if model is not None and hasattr(model, 'get_struct'):
                fok, freason = SolverFLD.forkJoinAdmits(model.get_struct(), method)
                if not fok:
                    return False, freason
        # 'default' IS ASKED THROUGH ITS RESOLUTION, not as a name of its own:
        # on a capped model it stands for 'dae' (see _resolve_default_method),
        # so gating the literal name would refuse the very run that succeeds.
        if ok and method in ('default', 'fluid.default') and self._blocked_resolves_to_dae():
            return ok, reason
        if ok and method not in ('dae', 'fluid.dae', 'mol', 'fluid.mol'):
            model = getattr(self, 'model', None)
            if model is not None and hasattr(model, 'get_used_lang_features'):
                ok, reason = NetworkSolver.checkBindingCapacity(model, 'SolverFLD')
                if not ok:
                    reason = ("%s Use options.method='dae', which carries the buffer as an "
                              "algebraic constraint on the drift." % reason)
        return ok, reason

    @staticmethod
    def forkJoinAdmits(sn, method):
        """Can ``method`` run the fluid fork-join fixed point on this model?

        A fork-join model is not integrated as one drift: the MMT transform
        replaces the fork by auxiliary classes and the answer is the fixed point
        of solving that transformed model repeatedly. On a CLOSED model the
        transform stays closed and every fluid method takes it. On an OPEN one
        the auxiliary classes arrive at a Source, and the DAE form has no
        unknowns for them: the inner solve fails on the class count rather than
        returning a drift, so the method is refused by name instead.

        'refined' is NOT listed here even though it fails the same way, because
        it is already refused on every open model, fork-join or not, by its own
        closed-model restriction (see getMethodFeatureSet).

        Called by runAnalyzer, so the run stops on it, and by
        supportsModelMethod, so a caller sees the same verdict before paying for
        the fixed point. One predicate, two callers. Mirrors MATLAB
        fluid_forkjoin_admits.

        Args:
            sn: NetworkStruct of the model.
            method: the concrete method name.

        Returns:
            (ok, reason); reason is '' when ok is True.
        """
        if method not in ('dae', 'fluid.dae'):
            return True, ''
        nodetype = np.ravel(np.asarray(sn.nodetype, dtype=int))
        if not np.any(nodetype == int(NodeType.FORK)):
            return True, ''
        if not np.any(np.isinf(np.ravel(np.asarray(sn.njobs, dtype=float)))):
            return True, ''
        return False, (
            "The dae method has no route through the fork-join fixed point on an OPEN "
            "model: the MMT transform hands the inner solve a mixed network whose "
            "auxiliary open classes the DAE form carries no unknowns for. Use "
            "options.method='minnormal', which is the same closure and does run that "
            "fixed point.")

    # The two single-station families and the shape each is stated for.
    _SINGLE_STATION_METHODS = ('ggisgi.fluid', 'fluid.ggisgi', 'ggisgi',
                               'ggingi.tga', 'fluid.tga', 'tga',
                               'tvms', 'fluid.tvms',
                               'mtginf', 'fluid.mtginf', 'mol', 'fluid.mol')
    _ABANDONMENT_METHODS = ('ggisgi.fluid', 'fluid.ggisgi', 'ggisgi',
                            'ggingi.tga', 'fluid.tga', 'tga',
                            'tvms', 'fluid.tvms')
    # The three limits that report a trajectory rather than a stationary point,
    # and so need a finite options.timespan; 'ggisgi' and 'tga' are stationary.
    _TIME_VARYING_METHODS = ('tvms', 'fluid.tvms',
                             'mtginf', 'fluid.mtginf', 'mol', 'fluid.mol')

    def _single_station_shape_admits(self, method):
        """The shape rule the single-station fluid limits are stated for, asked
        as a gate rather than raised at solve time.

        WHY IT IS HERE AND NOT IN listValidMethods, which is where the MATLAB
        twin puts it. That list is a @staticmethod in this port, deliberately
        (see its own note), so it cannot see the model; the gate can, and a
        gate is where a model-dependent rule belongs in any case. The two
        codebases therefore reach the same answer by different routes, which is
        what matters to a caller of findSolver: before this, SolverAUTO offered
        'fluid.ggingi.tga' on a plain M/M/1 and the method then raised, because
        the queue it needs customers to abandon has no patience law.

        The rule mirrors @SolverFLD/listValidMethods.m: one open class through
        one Source and one queueing station for the whole family, plus a
        reneging patience law for the two abandonment limits.
        """
        if method not in self._SINGLE_STATION_METHODS:
            return True, ''
        model = getattr(self, 'model', None)
        if model is None:
            return True, ''
        try:
            sn = model.get_struct()
            from ...api.solvers.fld.qsys import _station_of_type
            from ...api.sn import sn_patience_handles
            src = _station_of_type(sn, NodeType.SOURCE)
            qi = _station_of_type(sn, NodeType.QUEUE)
            if qi is None:
                qi = _station_of_type(sn, NodeType.DELAY)
        except Exception:
            # A struct this port cannot build here says nothing about the
            # shape; the analyzer's own check still stands behind the gate.
            return True, ''
        if src is None or qi is None or int(sn.nclasses) != 1 \
                or int(getattr(sn, 'nclosedjobs', 0) or 0) > 0:
            return False, ("The '%s' method is a single-station limit: it needs one open "
                           "class through one Source and one queueing station." % method)
        if method in self._ABANDONMENT_METHODS:
            try:
                h = sn_patience_handles(sn, qi, 0)
            except Exception:
                h = None
            if not h:
                return False, ("The '%s' method needs a reneging patience law on the queue "
                               "(Queue.setPatience): it is a limit for a queue customers "
                               "abandon." % method)
        return True, ''

    supports_model_method = supportsModelMethod

    def _fj_inner_solver(self, nonfjmodel, method=None):
        """Inner solve of the fork-join fixed point, on the fluid analyzer.

        Overrides ForkJoinDriverMixin._fj_inner_solver, whose default is
        SolverMVA. The transformed model carries no fork, so this never
        re-enters the fixed point. The requested method is carried through: the
        transform emits a plain mixed network, which every fluid method accepts
        except 'statedep', so a caller who asked for one gets it.
        """
        opts = SolverFLD.defaultOptions()
        # The driver passes method='amva' because it was written against the MVA
        # inner solve; that method name names no fluid method, so keep this solver's
        # own request instead of forwarding an MVA-only name.
        opts.method = self.options.method
        opts.verbose = self.options.verbose
        opts.iter_max = self.options.iter_max
        opts.iter_tol = self.options.iter_tol
        opts.tol = self.options.tol
        opts.stiff = self.options.stiff
        return SolverFLD(nonfjmodel, options=opts)

    def _fj_publish(self, result):
        """Store the fork-join result in the fluid result container.

        The driver speaks the plain-dict contract SolverMVA uses natively; the
        fluid getters read an FLDResult, so the dict is mapped onto its fields
        here. Only the steady-state means are carried: each pass of the fixed
        point integrates a DIFFERENT transformed network, so a trajectory read
        off the last pass would not be the trajectory of the model the caller
        built, and getTranAvg stays unavailable on a fork-join model.
        """
        self.result = FLDResult(
            QN=result['QN'], UN=result['UN'], RN=result['RN'], TN=result['TN'],
            CN=result['CN'], XN=result['XN'],
            AN=result.get('AN'), WN=result.get('WN'),
            t=None, QNt={}, UNt={}, TNt={},
            xvec=None, iterations=int(result.get('iter', 0)),
            runtime=float(result.get('runtime', 0.0)),
            method=str(result.get('method', 'mmt')),
        )
        return result

    @property
    def _sn(self):
        """Struct under the name the shared fork-join driver uses.

        SolverMVA and SolverNC keep the compiled struct in self._sn; SolverFLD
        keeps it in self.sn. Aliasing here is what lets the three share one
        ForkJoinDriverMixin rather than each carrying its own copy of the loop.
        """
        return self.sn

    @_sn.setter
    def _sn(self, value):
        self.sn = value

    @staticmethod
    def getFeatureSet() -> set:
        """Get set of features supported by the fluid solver.

        Returns the canonical feature names (mirrors MATLAB
        SolverFLD.getFeatureSet and the JAR SolverFluid).
        """
        return {
            'ClassSwitch', 'Delay', 'DelayStation', 'Queue',
            'Cache', 'CacheClassSwitcher',
            # 'CacheRetrieval' is deliberately NOT declared: no fluid code
            # anywhere implements delayed-hit retrieval, and on
            # examples/basic/cacheModel/retrieval_simple the ODE returned zero
            # QLen, Util and Tput on every row while jobs arrived at rate 1,
            # i.e. flow was not conserved. Refusing the model is the honest
            # answer; the JAR SolverFluid does the same.
            'Cox2', 'Coxian', 'Erlang', 'Exp', 'HyperExp',
            # MAP and MMPP2 are accepted at their stationary rate. Under the
            # 'closing' and 'matrix' methods a departure returns source mass
            # through the STATIONARY arrival-instant pie, which replaces D1' by
            # the rank-one map pie (x) (D1 e) -- that is the PH renewal process
            # (pie, D0), so flow stays conserved but the autocorrelation is lost,
            # exactly as in MATLAB SolverFLD. The 'kp' method does NOT lose it:
            # it carries the paper's own A0/A1 events, in which an arrival-
            # generating phase change acts through D1 itself.
            'APH', 'Det', 'MAP', 'MMPP2', 'NHPP', 'MAPt', 'PHt',
            # Non-Markovian renewal distributions: converted to acyclic PH by
            # sn_nonmarkov_toph in runAnalyzer, so the fluid ODE can solve them.
            'Gamma', 'Lognormal', 'Pareto', 'Uniform', 'Weibull',
            'StatelessClassSwitcher', 'InfiniteServer', 'SharedServer', 'Buffer', 'Dispatcher',
            'Server', 'ServiceTunnel',
            # Stochastic Petri nets: the 'dae' method only, see
            # getMethodFeatureSet. A Transition node routes the model to
            # methods/petri.py, which solves the marking as the same min-normal
            # closure with the P-invariants as constraints and the immediate
            # firing flows as algebraic unknowns. 'Storage'/'Linkage' ride along
            # with any Place, as they do in the SSA and CTMC sets, so declaring
            # Place without them refuses every Petri net at the gate.
            # 'Inhibiting' is declared, but an inhibitor arc on an IMMEDIATE
            # mode is answered wrongly when the inhibitor place's mean sits at
            # its threshold; see _kb/06-solver-catalog.md.
            'Place', 'Transition', 'Enabling', 'Inhibiting', 'Timing', 'Firing',
            'Storage', 'Linkage',
            # closing family only, see getMethodFeatureSet
            'LoadDependence',
            'SchedStrategy_INF', 'SchedStrategy_PS',
            'SchedStrategy_DPS', 'SchedStrategy_FCFS',
            # GPS is served only by the second-order closure: its share depends
            # on the backlog INDICATOR, which a first-order closure collapses to
            # the constant w_r/sum(w). See methods/minnormal.py.
            'SchedStrategy_GPS',
            # SIRO/LCFS/LCFSPR reach the matrix method, which builds a PS drift
            # -- the right aggregate for any work-conserving discipline. The
            # closing family has no drift branch for them and rejects them
            # explicitly in _ode_rate_factors rather than silently integrating
            # them as an infinite server.
            'SchedStrategy_SIRO', 'SchedStrategy_LCFS', 'SchedStrategy_LCFSPR',
            # Native fluid cache models: RANDOM(m)/FIFO(m) (refined mean field)
            # and strict FIFO(m) (position-resolved mean field). LRU/HLRU/CLIMB/
            # QLRU have no drift-based fluid model and are rejected at runtime.
            'ReplacementStrategy_RR', 'ReplacementStrategy_FIFO',
            'ReplacementStrategy_SFIFO',
            'RoutingStrategy_PROB', 'RoutingStrategy_RAND',
            'ClosedClass', 'SelfLoopingClass', 'Replayer',
            # Fork-join through the MMT transformation, driven by the shared
            # ForkJoinDriverMixin (as in SolverMVA and SolverNC). The transform
            # emits only Source, Delay, Queue, Router and ClassSwitch, all of
            # which the fluid drift already carries.
            'Fork', 'Forker', 'Join', 'Joiner',
            # quorum join: the MMT fixed point charges the k-th branch completion (fj_ordstat_exp)
            'JoinPartial',
            'RandomSource', 'Sink', 'Source', 'OpenClass', 'JobSink',
            # c-server stations: the drifts carry min(n,c); withdrawn from
            # 'diffusion' and 'mfq' in getMethodFeatureSet.
            'MultiServer',
            # A binding buffer: 'dae' carries it as an algebraic constraint,
            # 'mol' IS the Mt/G/s/0 loss system and the AoI arm of 'mfq' is a
            # bufferless or single-buffer queue. WHICH method serves one is the
            # structural rule supportsModelMethod asks and runAnalyzer stops on,
            # so no per-method delta duplicates it here.
            'FiniteCapacity',
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
            self._ensureAvgResults()

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
        """Get average residence times per station (M x K).

        Residence time is computed from response time using visit ratios:
        WN[ist,k] = RN[ist,k] * V[ist,k] / V[refstat,refclass]

        Returns:
            (M, K) array of residence times
        """
        if self.result is None:
            self._ensureAvgResults()

        # Compute ResidT using proper visit ratios from network structure
        if self.sn is not None and self.sn.visits:
            return sn_get_residt_from_respt(self.sn, self.result.RN, None)
        else:
            # Fallback: ResidT = RespT (no visit information available)
            return self.result.RN.copy()

    def getAvgWaitT(self) -> np.ndarray:
        """Get average waiting times per station.

        Waiting time is computed as response time minus mean service time.

        Returns:
            (M,) array of waiting times
        """
        if self.result is None:
            self._ensureAvgResults()

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
            self._ensureAvgResults()

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
            self._ensureAvgResults()

        if self.result.TN.ndim > 1:
            return np.sum(self.result.TN, axis=1)
        else:
            return np.full(self.result.QN.shape[0], np.mean(self.result.TN))

    # =====================================================================
    # PROBABILITY METHODS
    # =====================================================================

    def getMoments(self) -> Optional[Dict[str, Any]]:
        """Second-order results of the moment-closure methods.

        Mirrors MATLAB `@SolverFLD/getMoments`: state-level covariance Sigma,
        station-class queue-length variance QVar and standard deviation QStd,
        per-station population variance sigma2, and the state-coordinate index
        maps stationBlock/classBlock. None for every first-order method, which
        computes no second moment at all.
        """
        # A delegated result carries no `moments` attribute, so without this arm
        # the getter answered None for a minnormal solve that did compute the
        # covariance -- indistinguishable from a first-order method, which is the
        # one thing None is supposed to mean here.
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import moments_via_jar
            return moments_via_jar(self)

        if self.result is None:
            self._ensureAvgResults()
        return getattr(self.result, 'moments', None)

    def getProbAggr(self, ist: int) -> Tuple[float, float]:
        """Probability of the current per-class job distribution at a station.

        Returns P(n_1, ..., n_K at station ist) for the state the model is in.
        Two evaluations are available and the analysis that ran decides which,
        as `@SolverFLD/getProbAggr.m` does:

        moment closure ('minnormal') -- the solved state carries a covariance,
            so the JOINT law of the per-class populations at the station is the
            multivariate normal of the linear noise approximation and the
            answer is the probability it assigns to the unit cell around n.
            Correlation between the classes is accounted for.

        first-order methods -- no second moment exists, so the classes can only
            be treated as independent: Schmidt's binomial per closed class,
            Poisson (Delay) or multinomial-geometric (queue) per open class.

        Args:
            ist: Station index (1-based) or a station node

        Returns:
            (log_prob, prob)
        """
        from ...api.solvers.mva.prob_methods import get_prob_aggr
        from ...api.sn import SchedStrategy

        if not isinstance(ist, (int, np.integer)):
            ist = ist.get_station_index0() + 1

        # A delegated solve leaves no covariance behind: `result.moments` is a
        # native-python object and the JAR result container has none, so without
        # this arm a lang='java' minnormal solve fell through to the first-order
        # branch below and returned Schmidt's BINOMIAL under the name of the
        # moment closure -- a silent downgrade, not an error. Delegate the whole
        # query instead, so the engine that owns the covariance answers.
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import prob_via_jar
            # SolverFluid.getProbAggr reads its argument as a STATION index (it
            # derives the stateful row from sn.stationToStateful itself), unlike
            # SolverCTMC, which takes a node index -- hence raw_station. The cell
            # is named explicitly because model.json carries no initial state, so
            # a delegated query would otherwise be answered at the default one.
            return prob_via_jar(self, 'prob-aggr', ist=int(ist), kind='logtuple',
                                onebased=True, raw_station=True,
                                state=self._station_class_counts(int(ist)))

        if self.result is None:
            self._ensureAvgResults()

        sn = self._get_network_struct(self.model)

        # The moment closure supplies the joint law; a Source is excluded
        # because its coordinate is a normalisation constant rather than a
        # population and carries no covariance (see the minnormal terms).
        moments = getattr(self.result, 'moments', None)
        if moments and moments.get('Sigma') is not None and moments.get('classBlock') is not None \
                and int(sn.sched[int(ist) - 1]) != int(SchedStrategy.EXT) \
                and not self._has_open_class(sn, int(ist), moments):
            return self._gaussian_cell_prob(sn, int(ist), moments)

        class ResultAdapter:
            def __init__(self, outer):
                self.Q = outer.result.QN
                self.U = outer.result.UN
                self.R = outer.result.RN
                self.prob = None

        return get_prob_aggr(sn, ResultAdapter(self), ist)

    def _station_class_counts(self, ist: int) -> list:
        """Per-class job counts at station `ist` (1-based) in the model's state."""
        from ...api.state.marginal import toMarginal

        sn = self._get_network_struct(self.model)
        ind = int(np.asarray(sn.stationToNode).flatten()[ist - 1])
        isf = int(np.asarray(sn.nodeToStateful).flatten()[ind])
        state_i = np.atleast_2d(np.asarray(sn.state[isf], dtype=float))
        _, nir_m, _, _ = toMarginal(sn, ind, state_i)
        nir = np.asarray(nir_m).reshape(-1)[:int(sn.nclasses)]
        return [int(round(v)) for v in nir]

    @staticmethod
    def _has_open_class(sn, ist: int, moments) -> bool:
        """Whether an OPEN class is served at station ist.

        The Gaussian cell is used only where it beats the alternative. For an
        open class the first-order path is not an independence heuristic but the
        exact product form of the underlying queue -- geometric at a queue,
        Poisson at a Delay -- so replacing it by a normal approximation of the
        same law would be a loss: on M/M/1 at rho = 0.5 the product form is exact
        where the cell of the linear noise approximation returns 0.39 for the
        empty queue against 0.50. The closure earns its place on the CLOSED
        populations, where the alternative is Schmidt's binomial, itself an
        approximation, and where correlation between the classes is real.
        """
        njobs = np.asarray(sn.njobs, dtype=float).ravel()
        classBlock = moments['classBlock']
        for r in range(sn.nclasses):
            blk = np.asarray(classBlock[ist - 1][r], dtype=int).ravel()
            if blk.size and not np.isfinite(njobs[r]):
                return True
        return False

    def _gaussian_cell_prob(self, sn, ist: int, moments) -> Tuple[float, float]:
        """Joint probability of the per-class populations at station ist under
        the linear noise approximation solved by the moment closure.

        The state coordinates of class r at the station are
        moments['classBlock'][ist-1][r] (one per service phase), so the class
        population is their sum: its mean is the reported QN[ist-1,r] and the
        class-to-class covariance is the sum of the corresponding block of
        moments['Sigma']. The integer count n is then read off the continuous
        law as the unit cell [n-1/2, n+1/2], with the two ends extended to
        infinity at the boundaries of the state space, so that the mass the
        normal puts on negative populations lands on the empty station and the
        mass above a closed population lands on the full one.
        """
        from .mvn_rectangle import mvn_rectangle
        from ...api.state.marginal import toMarginal

        i = ist - 1
        K = sn.nclasses
        ind = int(np.asarray(sn.stationToNode).flatten()[i])
        isf = int(np.asarray(sn.nodeToStateful).flatten()[ind])
        state_i = np.atleast_2d(np.asarray(sn.state[isf], dtype=float))
        _, nir_m, _, _ = toMarginal(sn, ind, state_i)
        nir = np.asarray(nir_m).reshape(-1)[:K]

        Sigma = np.asarray(moments['Sigma'], dtype=float)
        classBlock = moments['classBlock']
        njobs = np.asarray(sn.njobs, dtype=float).ravel()

        idx = []
        m = []
        a = []
        b = []
        for r in range(K):
            blk = np.asarray(classBlock[i][r], dtype=int).ravel()
            if blk.size == 0:
                # the class has no service process here, so it has no
                # coordinate: any positive count is impossible
                if nir[r] > 0:
                    return -np.inf, 0.0
                continue
            idx.append(r)
            m.append(float(self.result.QN[i, r]))
            a.append(-np.inf if nir[r] <= 0 else nir[r] - 0.5)
            if np.isfinite(njobs[r]) and nir[r] >= njobs[r]:
                b.append(np.inf)
            else:
                b.append(nir[r] + 0.5)

        if not idx:
            return 0.0, 1.0

        nr = len(idx)
        C = np.zeros((nr, nr))
        for u in range(nr):
            bu = np.asarray(classBlock[i][idx[u]], dtype=int).ravel()
            for v in range(u, nr):
                bv = np.asarray(classBlock[i][idx[v]], dtype=int).ravel()
                C[u, v] = float(np.sum(Sigma[np.ix_(bu, bv)]))
                C[v, u] = C[u, v]

        prob, log_prob = mvn_rectangle(m, C, a, b)
        return log_prob, prob

    def getProbMarg(self, station: int, jobclass: int) -> np.ndarray:
        """Get marginal queue-length distribution at station for class.

        Args:
            station: Station index (0-based)
            jobclass: Job class index (0-based)

        Returns:
            Marginal probability vector P(n_ir) for n=0,1,2,...
        """
        if self.result is None:
            self._ensureAvgResults()

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
            self._ensureAvgResults()

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
            self._ensureAvgResults()

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

    def _get_aoi_results(self) -> Dict[str, Any]:
        result = self._ensure_result()
        aoi_results = getattr(result, 'aoiResults', None)
        if not aoi_results:
            raise RuntimeError(
                "No AoI results available. Ensure the model has a valid AoI topology and use method='mfq'."
            )
        return aoi_results

    def _evaluate_aoi_cdf(
        self,
        g: np.ndarray,
        a: np.ndarray,
        h: np.ndarray,
        t_values: np.ndarray,
    ) -> np.ndarray:
        # F(t) = 1 - S(t) with S(t) = -g expm(A t) inv(A) h. (g,A,h) is a
        # DENSITY triple -- g is normalized by -g inv(A) h, so g expm(A t) h is
        # the density and g inv(A)^2 h the mean -- and the survival function
        # carries the extra inv(A). Subtracting the density instead gives a
        # curve that falls before it rises; the maximum.accumulate that used to
        # sit on the return value MASKED exactly that, so it is gone with the
        # defect it hid.
        cdf_values = np.zeros_like(t_values, dtype=float)
        g_row = np.asarray(g, dtype=float).reshape(1, -1)
        h_col = np.asarray(h, dtype=float).reshape(-1, 1)
        amat = np.asarray(a, dtype=float)
        for idx, t in enumerate(t_values):
            if t <= 0:
                cdf_values[idx] = 0.0
            else:
                surv = np.linalg.solve(amat.T, (g_row @ linalg.expm(amat * float(t))).T).T @ h_col
                cdf_values[idx] = float(np.clip(1.0 + np.real_if_close(surv).item(), 0.0, 1.0))
        return cdf_values

    def getAvgAoI(self) -> Tuple[Dict[str, float], Dict[str, float], pd.DataFrame]:
        """Get average AoI and Peak AoI statistics."""
        aoi_results = self._get_aoi_results()

        aoi = {
            'mean': float(aoi_results['AoI_mean']),
            'var': float(aoi_results['AoI_var']),
        }
        aoi['std'] = float(np.sqrt(max(0.0, aoi['var'])))

        paoi = {
            'mean': float(aoi_results['PAoI_mean']),
            'var': float(aoi_results['PAoI_var']),
        }
        paoi['std'] = float(np.sqrt(max(0.0, paoi['var'])))

        table = pd.DataFrame({
            'Metric': ['AoI', 'Peak AoI'],
            'Mean': [aoi['mean'], paoi['mean']],
            'Variance': [aoi['var'], paoi['var']],
            'StdDev': [aoi['std'], paoi['std']],
            'SystemType': [aoi_results.get('systemType', ''), aoi_results.get('systemType', '')],
            'Preemption': [aoi_results.get('preemption', np.nan), aoi_results.get('preemption', np.nan)],
        })
        return aoi, paoi, table

    def getCdfAoI(self, t_values: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
        """Get AoI and Peak AoI CDFs as `[cdf, t]` arrays."""
        aoi_results = self._get_aoi_results()

        if t_values is None:
            mean_aoi = float(aoi_results.get('AoI_mean', np.nan))
            if not np.isfinite(mean_aoi) or mean_aoi <= 0:
                mean_aoi = 1.0
            t_values = np.linspace(0.0, 5.0 * mean_aoi, 200)
        t_values = np.asarray(t_values, dtype=float).reshape(-1)

        if any(aoi_results.get(key) is None or np.size(aoi_results.get(key)) == 0 for key in ('AoI_A', 'AoI_g', 'AoI_h')):
            raise RuntimeError('Matrix exponential parameters not available for AoI CDF computation.')

        aoi_cdf = self._evaluate_aoi_cdf(
            aoi_results['AoI_g'],
            aoi_results['AoI_A'],
            aoi_results['AoI_h'],
            t_values,
        )
        paoi_cdf = self._evaluate_aoi_cdf(
            aoi_results['PAoI_g'],
            aoi_results['PAoI_A'],
            aoi_results['PAoI_h'],
            t_values,
        )
        return np.column_stack((aoi_cdf, t_values)), np.column_stack((paoi_cdf, t_values))

    # =====================================================================
    # TRANSIENT METHODS
    # =====================================================================

    def _detect_nhpp_sources(self):
        """Build the nhpp_sched list for every Source station carrying a
        non-homogeneous arrival process.

        A non-homogeneous process is identified by the getRateSchedule method,
        as in MATLAB local_detect_nhpp. Returns a list of dicts with the sn
        station index, the class index and the process handle; empty when the
        model has none.
        """
        from ...lang.base import SchedStrategy as _SchedStrategy

        sched = []
        sn = self.sn
        model = getattr(self, 'model', None)
        if model is None or sn is None or sn.sched is None:
            return sched
        nodes = model.get_nodes()
        jobclasses = model.get_classes()
        for i in range(sn.nstations):
            if int(sn.sched[i]) != int(_SchedStrategy.EXT):
                continue
            node = nodes[int(sn.stationToNode[i])]
            arrivals = getattr(node, '_arrival_process', None)
            if not arrivals:
                continue
            for c in range(sn.nclasses):
                proc = arrivals.get(jobclasses[c])
                if proc is None or not hasattr(proc, 'getRateSchedule'):
                    continue
                # A MAPt or PHt also carries a rate schedule, but its matrix
                # entries vary independently, so one scalar per (station,class)
                # cannot express it and the closing method builds a per-event
                # multiplier instead. Listing it here too would apply both
                # channels and square the factor.
                if proc.getName() in ('MAPt', 'PHt'):
                    continue
                sched.append({'station': i, 'class': c, 'nhpp': proc})
        return sched

    def getTranAvgVar(self, *args):
        """Transient queue-length VARIANCE per station and class.

        Only the 'kp' method computes a second moment: it integrates the
        covariance of the Ko-Pender diffusion limit alongside the fluid mean.
        Returns (t, var) with var a dict keyed (station, class); the covariance
        between blocks is available as result.Sigmat.
        """
        if self.options.method not in ('kp', 'fluid.kp'):
            raise ValueError(
                "getTranAvgVar needs options.method='kp'; the other fluid "
                "methods integrate the mean only and carry no second moment.")
        if self.result is None or getattr(self.result, 'QVart', None) is None:
            self.runAnalyzer()
        return self.result.t, self.result.QVart

    def _has_matrix_schedule(self) -> bool:
        """Whether any station-class carries a MAPt or PHt.

        These are the schedule-bearing processes whose matrix entries vary
        independently, so the schedule reaches the ODE as a per-event multiplier
        built inside the closing method rather than through nhpp_sched.
        """
        from .utils.phase_type import is_mapt, is_pht

        sn = self.sn
        if sn is None or getattr(sn, 'procid', None) is None:
            return False
        for i in range(sn.nstations):
            for c in range(sn.nclasses):
                if is_mapt(sn, i, c) or is_pht(sn, i, c):
                    return True
        return False

    def getTranAvg(self, *args):
        """Get transient average metrics in MATLAB-compatible format.

        Args:
            *args: Optional transient handles (Qt, Ut, Tt) for MATLAB API compatibility.

        Returns:
            Tuple of (QNt, UNt, TNt) where each is a nested list [M][K] of TranResult objects.
        """
        if getattr(self.options, 'lang', 'python') == 'java':
            from ..jar_dispatch import tran_avg_via_jar
            return tran_avg_via_jar(self)
        if getattr(self.options, 'lang', 'python') == 'cpp':
            from ..cpp_dispatch import tran_avg_via_cpp
            return tran_avg_via_cpp(self)
        from ...constants import TranResult

        # The transient means are read off an INTEGRATED trajectory, so the
        # method has to be one that produces one. The moment closure solves its
        # mean and covariance to a fixed point and returns the converged state
        # only; without this switch the loop below falls through to the
        # steady-state value broadcast over the grid, i.e. a FLAT line reported
        # as a transient -- which a caller that couples stages through their
        # trajectories (SolverENV) cannot tell from a real one. MATLAB's
        # @SolverFLD/getTranAvg makes the same switch for method='default'.
        if self.sn is None:
            self.sn = self._get_network_struct(self.network)
        _resolved_tran = self._resolve_method()
        if _resolved_tran in ('minnormal', 'refined'):
            if self.options.method not in ('default', 'fluid.default'):
                from ...api.io.logging import line_warning
                line_warning('getTranAvg',
                             "method '%s' solves a fixed point and returns no trajectory; "
                             "integrating the closing ODE for the transient instead."
                             % self.options.method)
            self.options.method = 'closing'
            self.result = None

        # Detect NHPP (non-homogeneous Poisson) sources and pass their
        # intensity schedules to the closing ODE, so the transient tracks
        # lambda(t) rather than the baked-in time-average rate. Steady-state
        # getAvg is unaffected (no schedule injected there), consistent with
        # defining the NHPP steady state as its time average. Port of
        # matlab/src/solvers/FLD/@SolverFLD/getTranAvg.m.
        nhpp_sched = self._detect_nhpp_sources()
        if nhpp_sched:
            cfg = getattr(self.options, 'config', None)
            if cfg is None or not isinstance(cfg, dict):
                cfg = {}
                self.options.config = cfg
            if cfg.get('nhpp_sched') != nhpp_sched:
                cfg['nhpp_sched'] = nhpp_sched
                self.result = None
            # Only the closing ODE carries the per-event rate multiplier; TBI
            # integrates the same closing rates by cell decomposition and keeps
            # its method. Mirrors the MATLAB getTranAvg method switch.
            if self.options.method not in ('closing', 'tbi', 'fluid.tbi'):
                self.options.method = 'closing'
                self.result = None

        # An explicit per-(station,class) rate schedule (options.config
        # ['rate_sched'], e.g. injected by the LN coupled transient) reaches the
        # ODE through the same per-event rate multiplier, so it needs the same
        # method switch: the 'matrix' state-mapping method integrates fixed base
        # rates and would silently ignore the schedule.
        cfg = getattr(self.options, 'config', None)
        if isinstance(cfg, dict) and cfg.get('rate_sched') \
                and self.options.method not in ('closing', 'tbi', 'fluid.tbi'):
            self.options.method = 'closing'
            self.result = None

        # A MAPt or PHt needs the same switch: its per-event multiplier is built
        # inside the closing ODE, and the 'matrix' method would silently solve
        # the time-averaged nominal instead of the schedule.
        if self._has_matrix_schedule() and self.options.method not in ('closing', 'tbi', 'fluid.tbi'):
            self.options.method = 'closing'
            self.result = None

        if self.result is None:
            self._ensureAvgResults()

        # Cache networks carry their own transient through the RMF drift: the
        # queueing part uses the closing/matrix ODE above, while each cache's
        # per-class hit/miss trajectory comes from the mean-field drift over the
        # transient window. Port of MATLAB @SolverFLD/getTranAvg.m (hasCache).
        # Detection keys off the converged isolated inputs rather than
        # sn.nodetype, because _solve_rmf rewrites cache nodes to ClassSwitch.
        _cinputs = [ci for ci in getattr(self, '_cache_rmf_inputs', None) or [] if ci is not None]
        if _cinputs:
            timespan = getattr(self.options, 'timespan', None)
            if timespan is not None and len(timespan) >= 2 and np.isfinite(timespan[1]):
                tcache, hitprob_t, missprob_t, cnodes, arate = self._cacheqn_tran(timespan)
                self.result.CacheTran = {
                    't': tcache, 'hitprob': hitprob_t, 'missprob': missprob_t,
                    'nodes': cnodes, 'arate': arate,
                }

        M = self.result.QN.shape[0]
        K = self.result.QN.shape[1]
        t = self.result.t if hasattr(self.result, 't') and self.result.t is not None else np.array([0.0, 1000.0])

        QNt = [[None for _ in range(K)] for _ in range(M)]
        UNt = [[None for _ in range(K)] for _ in range(M)]
        TNt = [[None for _ in range(K)] for _ in range(M)]

        has_transient = (hasattr(self.result, 'QNt') and self.result.QNt and len(self.result.QNt) > 0)

        for i in range(M):
            for r in range(K):
                if has_transient and (i, r) in self.result.QNt:
                    QNt[i][r] = TranResult(t, self.result.QNt[(i, r)])
                    UNt[i][r] = TranResult(t, self.result.UNt.get((i, r), self.result.QNt[(i, r)]))
                    TNt[i][r] = TranResult(t, self.result.TNt.get((i, r), self.result.QNt[(i, r)]))
                else:
                    q_val = float(self.result.QN[i, r])
                    u_val = float(self.result.UN[i, r])
                    t_val = float(self.result.TN[i, r]) if self.result.TN.ndim > 1 else float(self.result.TN[r])
                    QNt[i][r] = TranResult(t, np.full(len(t), q_val))
                    UNt[i][r] = TranResult(t, np.full(len(t), u_val))
                    TNt[i][r] = TranResult(t, np.full(len(t), t_val))

        return QNt, UNt, TNt

    # =====================================================================
    # SAMPLING METHODS (Not Supported - Analytical Solver)
    # =====================================================================

    # =====================================================================
    # UNIFIED METRICS METHOD
    # =====================================================================


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

    def getAvgChainTable(self) -> pd.DataFrame:
        """Get average metrics by chain as DataFrame."""
        QN, UN, RN, WN, AN, TN = self.getAvgChain()

        nstations, nchains = QN.shape
        rows = []

        station_names = getattr(self.sn, 'nodenames', None) or [f'Station{i}' for i in range(nstations)]

        for i in range(nstations):
            for c in range(nchains):
                rows.append({
                    'Station': station_names[i] if i < len(station_names) else f'Station{i}',
                    'Chain': f'Chain{c + 1}',  # 1-based to match MATLAB
                    'QLen': QN[i, c],
                    'Util': UN[i, c],
                    'RespT': RN[i, c],
                    'ResidT': WN[i, c],
                    'ArvR': AN[i, c],
                    'Tput': TN[i, c],
                })

        # five SIGNIFICANT digits like MATLAB's table, not pandas' five decimals
        from line_solver.indexed_table import IndexedTable
        return IndexedTable(pd.DataFrame(rows))

    # =====================================================================
    # NODE-LEVEL METHODS
    # =====================================================================

    def getAvgNode(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Get average metrics per node.

        Unlike getAvg() which returns station-level metrics, this method
        returns node-level metrics (one row per node) including non-station
        nodes such as Cache and ClassSwitch. For Cache nodes, hit/miss class
        throughputs are computed using actual hit/miss probabilities.

        Returns:
            Tuple of (QNn, UNn, RNn, WNn, ANn, TNn) - node-level metrics
        """
        from ...api.sn.getters import sn_get_node_arvr_from_tput, sn_get_node_tput_from_tput
        from ...api.sn.network_struct import NodeType

        if self.result is None:
            self._ensureAvgResults()

        sn = self.sn
        I = sn.nnodes
        M = sn.nstations
        R = sn.nclasses

        # FLD solves caches via a ClassSwitch surrogate, which leaves
        # sn.nodetype reporting CLASSSWITCH for Cache nodes. Restore the true
        # node types from the network's (intact) node objects so the
        # node-expansion helpers can identify Cache nodes. Work on a copy so
        # FLD's own struct is left untouched.
        if sn.nodetype is not None and hasattr(self, 'network') and hasattr(self.network, '_nodes'):
            corrected = list(sn.nodetype)
            changed = False
            for ind, node in enumerate(self.network._nodes):
                if type(node).__name__ == 'Cache' and ind < len(corrected) \
                        and corrected[ind] != NodeType.CACHE:
                    corrected[ind] = NodeType.CACHE
                    changed = True
            if changed:
                import copy as _copy
                sn = _copy.copy(sn)
                sn.nodetype = corrected

        QN = self.result.QN
        UN = self.result.UN
        RN = self.result.RN
        TN = self.result.TN if self.result.TN.ndim > 1 else np.tile(self.result.TN, (M, 1))

        # FLD reports NaN throughput for zero-population closed classes (e.g.
        # the hit/miss helper classes of a Cache); their true throughput is
        # zero. Cleaning them up gives the node-expansion helpers the same
        # well-formed input the other solvers provide.
        TN = np.nan_to_num(np.asarray(TN, dtype=float), nan=0.0)

        # Residence times from response times using visit ratios
        if sn is not None and sn.visits:
            WN = sn_get_residt_from_respt(sn, RN, None)
        else:
            WN = RN.copy()

        # Publish per-class hit/miss probabilities onto the cache nodeparam so
        # the node-throughput helper can split hit/miss class flows. These come
        # from FLD's own result (_cacheHitProb/_cacheMissProb, one row per cache
        # node in ascending node order), not from the shared Cache node objects
        # whose hit ratios may have been overwritten by other solvers run on
        # the same model.
        cache_hit = getattr(self.result, '_cacheHitProb', None)
        cache_miss = getattr(self.result, '_cacheMissProb', None)
        if sn.nodeparam is not None and cache_hit is not None and cache_miss is not None:
            cache_hit = np.atleast_2d(np.asarray(cache_hit, dtype=float))
            cache_miss = np.atleast_2d(np.asarray(cache_miss, dtype=float))
            cidx = 0
            for ind in range(I):
                if sn.nodetype is not None and ind < len(sn.nodetype) \
                        and sn.nodetype[ind] == NodeType.CACHE and ind in sn.nodeparam:
                    if cidx < cache_hit.shape[0]:
                        cache_param = sn.nodeparam[ind]
                        cache_param.actualhitprob = cache_hit[cidx, :].flatten()
                        cache_param.actualmissprob = cache_miss[cidx, :].flatten()
                    cidx += 1

        # Throughput handle: 1 where the station-class has a valid throughput
        TH = np.zeros_like(TN)
        TH[TN > GlobalConstants.Zero] = 1.0

        # Node arrival rates and throughputs via shared helpers.
        # FLD's station-level result has no reliable arrival rates (ArvR is
        # NaN), so AN is left for the helper to derive from TN.
        ANn = sn_get_node_arvr_from_tput(sn, TN, TH)
        TNn = sn_get_node_tput_from_tput(sn, TN, TH, ANn)

        QNn = np.zeros((I, R))
        UNn = np.zeros((I, R))
        RNn = np.zeros((I, R))
        WNn = np.zeros((I, R))

        # Copy station metrics onto their station nodes
        stationToNode = np.asarray(sn.stationToNode).flatten()
        for ist in range(M):
            ind = int(stationToNode[ist]) if ist < len(stationToNode) else -1
            if 0 <= ind < I:
                QNn[ind, :] = QN[ist, :]
                UNn[ind, :] = UN[ist, :]
                RNn[ind, :] = RN[ist, :]
                WNn[ind, :] = WN[ist, :]

        # Fix arrival rates for ClassSwitch and Sink nodes for cache hit/miss
        # classes (matches MATLAB getAvgNode.m lines 54-76)
        for cacheInd in range(I):
            if sn.nodetype is not None and cacheInd < len(sn.nodetype) and sn.nodetype[cacheInd] == NodeType.CACHE:
                if sn.nodeparam is not None and cacheInd in sn.nodeparam:
                    cache_param = sn.nodeparam[cacheInd]
                    hitclass = np.atleast_1d(getattr(cache_param, 'hitclass', np.array([]))).flatten().astype(int)
                    missclass = np.atleast_1d(getattr(cache_param, 'missclass', np.array([]))).flatten().astype(int)
                    for ind in range(I):
                        if sn.nodetype[ind] == NodeType.CLASSSWITCH or sn.nodetype[ind] == NodeType.SINK:
                            for classIdx in range(R):
                                if np.any(hitclass[hitclass >= 0] == classIdx):
                                    ANn[ind, classIdx] = TNn[cacheInd, classIdx]
                                if np.any(missclass[missclass >= 0] == classIdx):
                                    ANn[ind, classIdx] = TNn[cacheInd, classIdx]

        return QNn, UNn, RNn, WNn, ANn, TNn

    def getAvgNodeTable(self) -> pd.DataFrame:
        """Get average metrics by node as DataFrame.

        Returns node-based results (one row per node per class) including
        non-station nodes such as Cache and ClassSwitch. All-zero rows are
        omitted, matching the other solvers' node tables.
        """
        QNn, UNn, RNn, WNn, ANn, TNn = self.getAvgNode()

        sn = self.sn
        nodenames = list(sn.nodenames) if hasattr(sn, 'nodenames') and sn.nodenames else []
        class_names = list(sn.classnames) if hasattr(sn, 'classnames') and sn.classnames else \
                      [f'Class{r}' for r in range(sn.nclasses)]

        rows = []
        for node_idx in range(sn.nnodes):
            node_name = nodenames[node_idx] if node_idx < len(nodenames) else f'Node{node_idx}'
            for r in range(sn.nclasses):
                class_name = class_names[r] if r < len(class_names) else f'Class{r}'

                # Filter out all-zero rows
                if abs(QNn[node_idx, r]) < 1e-10 and abs(UNn[node_idx, r]) < 1e-10 and \
                   abs(RNn[node_idx, r]) < 1e-10 and abs(ANn[node_idx, r]) < 1e-10 and abs(TNn[node_idx, r]) < 1e-10:
                    continue

                rows.append({
                    'Node': node_name,
                    'JobClass': class_name,
                    'QLen': QNn[node_idx, r],
                    'Util': UNn[node_idx, r],
                    'RespT': RNn[node_idx, r],
                    'ResidT': WNn[node_idx, r],
                    'ArvR': ANn[node_idx, r],
                    'Tput': TNn[node_idx, r],
                })

        df = pd.DataFrame(rows)
        result = IndexedTable(df)

        if len(df) > 0 and not getattr(self, '_table_silent', False):
            print(result)

        return result

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

    getAvgSysTable = NetworkSolver.getAvgSysTable  # chain-level shared layout

    # =====================================================================
    # PASSAGE TIME METHODS
    # =====================================================================

    def getCdfPassT(self, station: Optional[int] = None, job_class: Optional[int] = None,
                    t_span: Optional[Tuple[float, float]] = None):
        """Get the steady-state passage time CDF for a station/class.

        The passage time of a job of class r at station i is the time from its
        arrival at the station to its departure from it. Under the fluid
        approximation this passage is exactly the quantity reported by
        getCdfRespT: both come from the transient passage time analysis started
        from the steady-state ODE solution, so this delegates to it rather than
        duplicating the computation. The two names are kept distinct because the
        solver interface declares both, and other solvers may separate them.

        For the passage time along a prescribed route, i.e. conditional on a
        given sequence of nodes rather than at a single station, no method is
        provided: the fluid passage time analysis is per station and combining
        stations would require an independence assumption across them.

        Parameters
        ----------
        station : int, optional
            Station index. If None, returns the CDF for all stations.
        job_class : int, optional
            Job class index. If None, returns the CDF for all classes.
        t_span : tuple, optional
            Time interval (t_min, t_max) for CDF evaluation. If None, it is
            estimated from the mean response time.

        Returns
        -------
        When station and job_class are both None:
            List of lists where RD[station][class] is a 2D array with columns
            [cdf, time]
        When station and job_class are specified:
            dict with keys 't', 'cdf', 'mean', 'var', 'method'

        See Also
        --------
        getCdfRespT : steady-state response time distribution (same quantity)
        getTranCdfPassT : passage time distribution during the transient
        """
        return self.getCdfRespT(station=station, job_class=job_class, t_span=t_span)

    def getCdfPT(self, station: Optional[int] = None, job_class: Optional[int] = None,
                 t_span: Optional[Tuple[float, float]] = None):
        """Get the steady-state passage time CDF for a station/class.

        Backward-compatible name for getCdfPassT, to which this delegates. See
        getCdfPassT for the contract.
        """
        return self.getCdfPassT(station=station, job_class=job_class, t_span=t_span)

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

    def GetAvgAoI(self) -> Tuple[Dict[str, float], Dict[str, float], pd.DataFrame]:
        """Alias for getAvgAoI (MATLAB compatibility)."""
        return self.getAvgAoI()

    def GetCdfAoI(self, t_values: Optional[np.ndarray] = None) -> Tuple[np.ndarray, np.ndarray]:
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

    # Passage time aliases
    def GetCdfPT(self, station: Optional[int] = None, job_class: Optional[int] = None,
                 t_span: Optional[Tuple[float, float]] = None):
        """Alias for getCdfPT (MATLAB compatibility)."""
        return self.getCdfPT(station=station, job_class=job_class, t_span=t_span)

    def GetCdfPassT(self, station: Optional[int] = None, job_class: Optional[int] = None,
                    t_span: Optional[Tuple[float, float]] = None):
        """Alias for getCdfPassT (MATLAB compatibility)."""
        return self.getCdfPassT(station=station, job_class=job_class, t_span=t_span)

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
    avg_node_table = getAvgNodeTable
    avg_chain_table = getAvgChainTable
    avg_node_chain_table = getAvgNodeChainTable
    avg_sys_table = getAvgSysTable
    run_analyzer = runAnalyzer
    cdf_resp_t = getCdfRespT
    cdf_respt = getCdfRespT
    get_cdf_resp_t = getCdfRespT
    perct_resp_t = getPerctRespT
    perct_respt = getPerctRespT
    avg_qlen = getAvgQLen
    avg_util = getAvgUtil
    avg_respt = getAvgRespT
    get_avg_respt = getAvgRespT
    avg_tput = getAvgTput


__all__ = [
    'SolverFLD',
    'SolverFLDOptions',
    'FLDResult',
]
