"""
Differential Evolution Solver for Network Optimization

This module provides the main optimization solver using scipy's
differential_evolution algorithm.

Key Classes:
    - LineOptSolver: Main optimization solver

Example:
    >>> solver = LineOptSolver(problem)
    >>> result = solver.solve()
    >>> print(f"Optimal objective: {result.objective_value}")
"""

import time
import numpy as np
from typing import Dict, List, Any, Optional, Callable
from scipy.optimize import differential_evolution

from .variables import DecisionVariable
from .objectives import Objective, Constraint
from .evaluator import LineEvaluator
from .results import OptimizationResult, EvaluationResult


class _TimeLimitReached(Exception):
    """Internal signal raised when the solver time limit expires mid-run."""
    pass


class LineOptSolver:
    """
    Main optimization solver using scipy differential_evolution.

    Follows sequor's DifferentialEvolution pattern with:
    - Random key encoding for discrete variables
    - Penalty-based constraint handling
    - Convergence tracking

    Args:
        problem: OptimizationProblem to solve
        **options: Solver options (see defaultOptions)

    Example:
        >>> solver = LineOptSolver(problem, verbose=True)
        >>> result = solver.solve()
    """

    @staticmethod
    def defaultOptions() -> dict:
        """
        Get default solver options.

        Returns:
            Dict with default option values
        """
        return {
            'strategy': 'best1bin',
            'popsize': 15,
            'mutation': (0.5, 1.0),
            'recombination': 0.7,
            'tol': 0.01,
            'max_iterations': 100,
            'time_limit': 300.0,
            'seed': None,
            'verbose': False,
            'penalty_weight': 1e6,
            'polish': False,
            'workers': 1,
            # Aggregation across workload scenarios: 'worst' (minimax) or
            # 'mean' (weighted average)
            'scenario_aggregation': 'worst',
            # see _kb/05-solvers-overview.md (LineOpt) for rationale
            'optimizer': 'evolution',
            'fd_step': 1e-6,          # finite-difference step in encoded space
            # Step for a differencing that RE-SOLVES a LayeredNetwork: the
            # SolverLN fixed point is smooth only above its own noise floor
            # (~1e-7 in a utilization under lang='java'), and penalty_weight
            # amplifies that noise by 1e6, so a 1e-6 step returns the noise
            # rather than the derivative. See _kb/05-solvers-overview.md.
            'fd_step_layered': 1e-4,
            'gradient_restarts': 4,   # multistart count for the gradient path
            # LayeredNetwork (LQN) gradient source (used only when the model is
            # a LayeredNetwork and the gradient path is taken):
            #   'fd'             finite-difference the whole LayeredNetwork per
            #                    parameter -- correct total derivative, robust
            #                    default (N extra ensemble solves).
            #   'partial_sens'   assemble the direction from SolverLN's per-layer
            #                    WITHIN-LAYER partial service-rate derivatives --
            #                    cheap, but biased (omits cross-layer coupling).
            #   'partial_plus_fd' partial_sens direction, corrected by a full
            #                    LayeredNetwork finite difference every
            #                    'fd_refresh' gradient evaluations.
            'lqn_gradient': 'fd',
            'fd_refresh': 5,          # partial_plus_fd full-FD correction period
            # Explicit layer freezing (LQN only): a list of layer names (host
            # or task layers) whose variables are held at the model's current
            # value instead of being optimized. None freezes nothing.
            'frozen_layers': None,
        }

    # Decision-variable types that are continuous and therefore differentiable.
    # Includes the continuous LQN knobs (host demand, think time).
    _CONTINUOUS_TYPES = frozenset({'service_rate', 'routing',
                                   'host_demand', 'think_time'})

    def __init__(self, problem: 'OptimizationProblem', **options):
        """
        Initialize solver.

        Args:
            problem: OptimizationProblem to solve
            **options: Override default options
        """
        self._problem = problem
        self._options = {**self.defaultOptions(), **options}

        # Fixed variables (set by decomposition) applied on every evaluation
        # and visible to objectives/constraints by name
        fixed = []
        fixed_dict = {}
        if hasattr(problem, 'getFixedVariables'):
            fixed = list(problem.getFixedVariables())
            fixed_dict = {var.getName(): value for var, value in fixed}

        # Explicit layer freezing: variables whose layer is in ``frozen_layers``
        # are held at the model's current parameter value (moved from free to
        # fixed) instead of being optimized. See _freezeLayers.
        free_vars = list(problem.getVariables())
        free_vars, fixed, fixed_dict = self._freezeLayers(
            free_vars, fixed, fixed_dict)
        self._free_variables = free_vars
        self._fixed_value_dict = fixed_dict

        # Evaluators: base model plus one per workload scenario
        self._evaluator = LineEvaluator(problem.getModel(), free_vars, fixed)
        self._evaluators = [self._evaluator]
        self._scenario_weights = [1.0]
        if hasattr(problem, 'getScenarios'):
            for scenario_model, weight in problem.getScenarios():
                self._evaluators.append(LineEvaluator(
                    scenario_model, free_vars, fixed))
                self._scenario_weights.append(weight)

        # Tracking variables
        self._iterations = 0
        self._best_value = float('inf')
        self._best_x: Optional[np.ndarray] = None
        self._convergence_history = []
        self._start_time = 0.0
        self._deadline = float('inf')

        # Per-evaluator caches for expensive evaluations
        self._caches: List[Dict[tuple, EvaluationResult]] = [
            {} for _ in self._evaluators]
        self._cache = self._caches[0]

        # LQN gradient bookkeeping: a per-configuration sensitivity cache and a
        # gradient-evaluation counter driving the partial_plus_fd refresh.
        self._lqn_sens_cache: Dict[tuple, Any] = {}
        self._grad_calls = 0

    def _freezeLayers(self, free_vars, fixed, fixed_dict):
        """Partition variables by the ``frozen_layers`` option (LQN only).

        A variable whose layer set (``DecisionVariable.getLayer``) intersects
        the frozen set is moved from the free list into the fixed list, held at
        its current model parameter value (``DecisionVariable.currentValue``).
        If the current value cannot be read the variable is simply dropped from
        the optimization, leaving the model's built-in value untouched. Returns
        the updated (free_vars, fixed, fixed_dict). A no-op when the model is
        flat or no layers are frozen.
        """
        frozen = self._options.get('frozen_layers')
        if not frozen:
            return free_vars, fixed, fixed_dict
        from .layered import is_layered
        model = self._problem.getModel()
        if not is_layered(model):
            return free_vars, fixed, fixed_dict

        frozen_set = set(frozen)
        already = {var.getName() for var, _ in fixed}
        kept = []
        for var in free_vars:
            layers = var.getLayer(model) or []
            if frozen_set.intersection(layers):
                if var.getName() in already:
                    continue
                value = var.currentValue(model)
                if value is not None:
                    fixed.append((var, value))
                    fixed_dict[var.getName()] = value
                # else: drop the variable; the model keeps its built-in value
            else:
                kept.append(var)
        return kept, fixed, fixed_dict

    def getProblem(self) -> 'OptimizationProblem':
        """Get the optimization problem."""
        return self._problem

    def getOptions(self) -> dict:
        """Get current solver options."""
        return self._options

    def setOption(self, key: str, value: Any) -> 'LineOptSolver':
        """
        Set a solver option.

        Args:
            key: Option name
            value: Option value

        Returns:
            self for chaining
        """
        self._options[key] = value
        return self

    def solve(self) -> OptimizationResult:
        """
        Execute differential evolution optimization.

        Returns:
            OptimizationResult with optimal solution and statistics
        """
        self._start_time = time.time()
        self._iterations = 0
        self._best_value = float('inf')
        self._best_x = None
        self._convergence_history = []
        self._caches = [{} for _ in self._evaluators]
        self._cache = self._caches[0]
        self._lqn_sens_cache = {}
        self._grad_calls = 0
        self._deadline = self._start_time + self._options['time_limit']

        # Get bounds from evaluator
        bounds = self._evaluator.getBounds()

        if len(bounds) == 0:
            # No variables to optimize
            return self._buildEmptyResult()

        # Dispatch to the gradient-based path for continuous problems.
        if self._shouldUseGradient():
            if self._options['verbose']:
                print(f"line_solver.opt: Starting gradient optimization with "
                      f"{len(bounds)} dimensions")
            return self._solveGradient(bounds)

        if self._options['verbose']:
            print(f"line_solver.opt: Starting optimization with {len(bounds)} dimensions")
            print(f"line_solver.opt: strategy={self._options['strategy']}, "
                  f"popsize={self._options['popsize']}")

        # Set random seed if specified
        seed = self._options['seed']

        timed_out = False
        try:
            # Run differential evolution
            result = differential_evolution(
                self._objectiveFunction,
                bounds,
                strategy=self._options['strategy'],
                maxiter=self._options['max_iterations'],
                popsize=self._options['popsize'],
                mutation=self._options['mutation'],
                recombination=self._options['recombination'],
                tol=self._options['tol'],
                seed=seed,
                callback=self._callback,
                disp=False,
                polish=self._options['polish'],
                workers=self._options['workers'],
            )
            best_x, best_fun = result.x, result.fun
        except _TimeLimitReached:
            # Deadline hit mid-generation: fall back to best-so-far
            timed_out = True
            if self._best_x is None:
                return self._buildEmptyResult()
            best_x, best_fun = self._best_x, self._best_value

        solve_time = time.time() - self._start_time

        # Build final result
        opt_result = self._buildResult(best_x, best_fun, solve_time)
        if timed_out:
            opt_result.terminated_by = 'time_limit'

        if self._options['verbose']:
            print(f"line_solver.opt: Finished - objective={opt_result.objective_value:.4f}, "
                  f"iterations={opt_result.iterations}, time={solve_time:.2f}s")

        return opt_result

    def _objectiveFunction(self, x: np.ndarray) -> float:
        """
        Objective function for scipy differential_evolution.

        Evaluates the model (all scenarios) and computes the aggregated
        objective with constraint penalties. Enforces the time limit on a
        per-evaluation basis.

        Args:
            x: Decision vector

        Returns:
            Objective value (with penalties)
        """
        # Enforce time limit per evaluation, not just per generation
        if time.time() >= self._deadline:
            raise _TimeLimitReached()

        # Decode free variable values; merge fixed values for objectives
        variable_values = self._evaluator.decodeVariables(x)
        all_values = {**self._fixed_value_dict, **variable_values}

        objective = self._problem.getObjective()
        penalty_weight = self._options['penalty_weight']

        scenario_values = []
        for evaluator, cache in zip(self._evaluators, self._caches):
            eval_result = evaluator.evaluateValuesWithCache(variable_values,
                                                            cache)
            if not eval_result.feasible:
                return float('inf')

            value = objective.evaluateWithPenalty(
                eval_result, all_values, penalty_weight)
            for constraint in self._problem.getConstraints():
                violation = constraint.evaluate(eval_result, all_values)
                value += violation * penalty_weight
            scenario_values.append(value)

        total = self._aggregateScenarios(scenario_values)

        # Track best-so-far for time-limit fallback
        if total < self._best_value:
            self._best_value = total
            self._best_x = np.array(x, copy=True)

        return total

    def _aggregateScenarios(self, values: List[float]) -> float:
        """Aggregate per-scenario objective values ('worst' or 'mean')."""
        if len(values) == 1:
            return values[0]
        if self._options['scenario_aggregation'] == 'mean':
            total_weight = sum(self._scenario_weights)
            return sum(w * v for w, v in zip(self._scenario_weights,
                                             values)) / total_weight
        # 'worst': minimax for a minimization problem
        return max(values)

    # ---- gradient-based optimization ------------------------------------

    def _shouldUseGradient(self) -> bool:
        """Decide whether to use the gradient path for this problem."""
        opt = self._options.get('optimizer', 'auto')
        if opt == 'gradient':
            return True
        if opt == 'evolution':
            return False
        # 'auto': gradient only when every decision variable is continuous
        return self._allContinuous()

    def _allContinuous(self) -> bool:
        """True if every decision variable is continuous (differentiable)."""
        variables = self._free_variables
        if not variables:
            return False
        for var in variables:
            vtype = var.getVariableType() if hasattr(var, 'getVariableType') \
                else None
            if vtype not in self._CONTINUOUS_TYPES:
                return False
        return True

    def _solveGradient(self, bounds: List[tuple]) -> OptimizationResult:
        """Minimize the penalized objective with scipy L-BFGS-B and gradients.

        Uses the same penalized scalar objective as the evolutionary path
        (constraints enter as penalties), so no separate constraint handling
        is needed. The Jacobian is analytic where the model is product-form
        and finite-difference otherwise (see :meth:`_objectiveGradient`).
        A deterministic multistart guards against poor local minima.
        """
        from scipy.optimize import minimize

        dim = len(bounds)
        rng = np.random.default_rng(self._options['seed'])
        n_start = max(1, int(self._options['gradient_restarts']))
        starts = [np.full(dim, 0.5)]
        for _ in range(n_start - 1):
            starts.append(rng.random(dim))

        best_x = None
        best_fun = float('inf')
        timed_out = False
        for x0 in starts:
            if time.time() >= self._deadline:
                timed_out = True
                break
            try:
                res = minimize(
                    self._objectiveFunction, x0,
                    jac=self._objectiveGradient,
                    method='L-BFGS-B', bounds=bounds,
                    options={'maxiter': self._options['max_iterations']},
                )
                fx, xx = res.fun, res.x
            except _TimeLimitReached:
                timed_out = True
                break
            if np.isfinite(fx) and fx < best_fun:
                best_fun = fx
                best_x = np.array(xx, copy=True)

        # Fall back to best-so-far tracked during objective evaluations.
        if best_x is None or (self._best_x is not None
                              and self._best_value < best_fun):
            if self._best_x is not None:
                best_x, best_fun = self._best_x, self._best_value
            elif best_x is None:
                return self._buildEmptyResult()

        solve_time = time.time() - self._start_time
        opt_result = self._buildResult(best_x, best_fun, solve_time)
        if timed_out:
            opt_result.terminated_by = 'time_limit'
        if self._options['verbose']:
            print(f"line_solver.opt: Finished (gradient) - "
                  f"objective={opt_result.objective_value:.4f}, "
                  f"time={solve_time:.2f}s")
        return opt_result

    def _objectiveGradient(self, x: np.ndarray) -> np.ndarray:
        """Gradient of the penalized objective in encoded space.

        For a LayeredNetwork the source is selected by the ``lqn_gradient``
        option: 'fd' (whole-model finite difference, the robust default),
        'partial_sens' (SolverLN per-layer partial derivatives, cheap/biased),
        or 'partial_plus_fd' (partial with a periodic full-FD correction). For a
        flat network it prefers the analytic product-form gradient. All paths
        fall back to central finite differences, which always work.
        """
        if getattr(self._evaluator, 'is_layered', False):
            mode = self._options.get('lqn_gradient', 'fd')
            if mode in ('partial_sens', 'partial_plus_fd'):
                self._grad_calls += 1
                refresh = max(1, int(self._options.get('fd_refresh', 5)))
                # partial_plus_fd periodically replaces the biased partial
                # direction with the correct whole-model finite difference.
                if not (mode == 'partial_plus_fd'
                        and self._grad_calls % refresh == 0):
                    g = self._lqnAnalyticGradient(x)
                    if g is not None:
                        return g
            return self._finiteDifferenceGradient(x)

        analytic = self._analyticGradient(x)
        if analytic is not None:
            return analytic
        return self._finiteDifferenceGradient(x)

    def _finiteDifferenceGradient(self, x: np.ndarray) -> np.ndarray:
        """Central finite-difference gradient of the penalized scalar objective.

        Works for any model or solver (for a LayeredNetwork each perturbed
        evaluation re-solves the whole ensemble, giving the correct total
        derivative). One-sided differences are used near an
        infeasible/unstable boundary where a two-sided value is non-finite.
        """
        # A layered evaluation re-solves an iterative fixed point, so the step
        # must clear its noise floor (see the 'fd_step_layered' default).
        h = float(self._options['fd_step_layered']
                  if getattr(self._evaluator, 'is_layered', False)
                  else self._options['fd_step'])
        dim = len(x)
        g = np.zeros(dim)
        f0 = None
        for i in range(dim):
            xp = np.array(x, copy=True)
            xm = np.array(x, copy=True)
            xp[i] = min(1.0, x[i] + h)
            xm[i] = max(0.0, x[i] - h)
            fp = self._objectiveFunction(xp)
            fm = self._objectiveFunction(xm)
            if np.isfinite(fp) and np.isfinite(fm) and xp[i] > xm[i]:
                g[i] = (fp - fm) / (xp[i] - xm[i])
                continue
            # one-sided fallback at a non-finite boundary
            if f0 is None:
                f0 = self._objectiveFunction(x)
            if np.isfinite(fp) and np.isfinite(f0) and xp[i] > x[i]:
                g[i] = (fp - f0) / (xp[i] - x[i])
            elif np.isfinite(fm) and np.isfinite(f0) and x[i] > xm[i]:
                g[i] = (f0 - fm) / (x[i] - xm[i])
            else:
                g[i] = 0.0
        return g

    def _analyticGradient(self, x: np.ndarray) -> Optional[np.ndarray]:
        """Analytic gradient of the penalized objective, or None.

        Assembles d(objective)/dx from three exact pieces and NO extra solver
        calls (one model solve, cached): (i) analytic d(metric)/d(rate) from
        product-form sensitivities attached to the evaluation result;
        (ii) d(penalized scalar)/d(metric) and d(penalized scalar)/d(value),
        each obtained by cheap finite differences of the scalar objective in
        metric/value space (pure arithmetic, no solving); (iii) d(rate)/dx from
        the variable's linear decode. Returns None (falling back to
        finite-difference-in-x) when the problem is multi-scenario, a variable
        is not a rate variable, or sensitivities are unavailable.
        """
        variables = self._free_variables
        # supported: single scenario, every variable is a rate variable
        if len(self._evaluators) != 1 or not variables:
            return None
        for var in variables:
            if var.getVariableType() != 'service_rate' \
               or not hasattr(var, 'paramKey'):
                return None

        values = self._evaluator.decodeVariables(x)
        all_values = {**self._fixed_value_dict, **values}
        res = self._evaluator.evaluateValuesWithCache(values, self._caches[0])
        sens = getattr(res, 'sensitivities', None)
        if not res.feasible or sens is None:
            return None

        objective = self._problem.getObjective()
        pw = self._options['penalty_weight']

        def scalar(result, vals):
            v = objective.evaluateWithPenalty(result, vals, pw)
            for c in self._problem.getConstraints():
                v += c.evaluate(result, vals) * pw
            return v

        # metric kind -> the EvaluationResult dict it is stored in
        metric_dicts = {
            'RespT': res.response_times,
            'QLen': res.queue_lengths,
            'Tput': res.throughputs,
            'Util': res.utilizations,
        }
        h = float(self._options['fd_step'])

        def scalar_metric_derivative(kind, mkey):
            """d(penalized scalar)/d(metric[kind][mkey]) by central FD."""
            d = metric_dicts.get(kind)
            if d is None or mkey not in d:
                return 0.0
            base = d[mkey]
            d[mkey] = base + h
            fp = scalar(res, all_values)
            d[mkey] = base - h
            fm = scalar(res, all_values)
            d[mkey] = base
            return (fp - fm) / (2.0 * h)

        # d(scalar)/d(rate parameter) accumulated over all metrics
        dS_drate: Dict[tuple, float] = {}
        for kind, table in sens.items():
            for mkey, pmap in table.items():
                dSdm = scalar_metric_derivative(kind, mkey)
                if dSdm == 0.0:
                    continue
                for pkey, dmetric_dparam in pmap.items():
                    dS_drate[pkey] = dS_drate.get(pkey, 0.0) \
                        + dSdm * dmetric_dparam

        grad = np.zeros(len(x))
        for i, var in enumerate(variables):
            name = var.getName()
            pkey = var.paramKey()
            # direct dependence of the scalar on the decoded value (e.g. cost),
            # holding metrics fixed: central FD in value space (no solving)
            base_val = all_values.get(name)
            direct = 0.0
            if isinstance(base_val, (int, float)):
                all_values[name] = base_val + h
                fp = scalar(res, all_values)
                all_values[name] = base_val - h
                fm = scalar(res, all_values)
                all_values[name] = base_val
                direct = (fp - fm) / (2.0 * h)
            dS_dvalue = direct + dS_drate.get(pkey, 0.0)
            grad[i] = dS_dvalue * var.decodeJacobian(np.array([x[i]]))
        return grad

    def _lqnAnalyticGradient(self, x: np.ndarray) -> Optional[np.ndarray]:
        """Partial-sensitivity gradient for a LayeredNetwork, or None.

        Assembles d(objective)/dx from SolverLN's per-layer WITHIN-LAYER
        service-rate derivatives, WITHOUT re-solving per parameter. For each
        continuous host-demand variable it (i) reads d(layer metric)/d(rate)
        from the sensitivity table at the variable's host-layer row, (ii) maps
        each layer metric to the LQN node metric it approximates and takes
        d(penalized scalar)/d(that node metric) by cheap metric-space finite
        differences (pure arithmetic, no solving), (iii) chains through
        d(rate)/d(demand) = -1/D^2 and the linear decode d(demand)/dx.

        Returns None -- so the caller finite-differences the whole model --
        when the problem is multi-scenario, any variable is not a host-demand
        variable, or the sensitivity table is unavailable. This is a BIASED
        estimate of the total derivative (it omits the cross-layer fixed-point
        coupling); 'fd'/'partial_plus_fd' exist to correct that.
        """
        variables = self._free_variables
        if len(self._evaluators) != 1 or not variables:
            return None
        # Only host-demand variables carry the partial-sensitivity hooks.
        for var in variables:
            if var.getVariableType() != 'host_demand' \
               or not hasattr(var, 'sensKey'):
                return None

        values = self._evaluator.decodeVariables(x)
        all_values = {**self._fixed_value_dict, **values}
        res = self._evaluator.evaluateValuesWithCache(values, self._caches[0])
        if not res.feasible:
            return None

        # Per-configuration sensitivity cache (the table is expensive).
        skey = self._evaluator._valuesKey(values)
        if skey in self._lqn_sens_cache:
            sens = self._lqn_sens_cache[skey]
        else:
            sens = self._evaluator.evaluateLayeredSensitivities(values)
            self._lqn_sens_cache[skey] = sens
        if not sens:
            return None

        objective = self._problem.getObjective()
        pw = self._options['penalty_weight']
        h = float(self._options['fd_step'])
        model = self._problem.getModel()

        def scalar():
            v = objective.evaluateWithPenalty(res, all_values, pw)
            for c in self._problem.getConstraints():
                v += c.evaluate(res, all_values) * pw
            return v

        metric_dicts = {
            'RespT': res.response_times,
            'QLen': res.queue_lengths,
            'Tput': res.throughputs,
            'Util': res.utilizations,
        }

        def scalar_metric_derivative(kind, mkey):
            """d(penalized scalar)/d(node metric[kind][mkey]) by central FD."""
            d = metric_dicts.get(kind)
            if d is None or mkey is None or mkey not in d:
                return 0.0
            base = d[mkey]
            d[mkey] = base + h
            fp = scalar()
            d[mkey] = base - h
            fm = scalar()
            d[mkey] = base
            return (fp - fm) / (2.0 * h)

        grad = np.zeros(len(x))
        for i, var in enumerate(variables):
            # The row key depends on how the LN method named the layer classes
            # (activity under 'srvn.cs', caller task under 'srvn.ph'), so the
            # candidates are tried in order rather than assumed.
            candidates = var.sensKeys(model) if hasattr(var, 'sensKeys') \
                else [var.sensKey(model)]
            skeyv = next((k for k in candidates if k is not None and k in sens), None)
            if skeyv is None:
                # No sensitivity row for this variable: leave 0 (L-BFGS-B will
                # still make progress on the others; FD modes cover it fully).
                continue
            row = sens[skeyv]
            targets = var.sensMetricTargets(model)
            dS_drate = 0.0
            for kind, dmetric_drate in row.items():
                if dmetric_drate == 0.0:
                    continue
                dS_dmetric = scalar_metric_derivative(kind, targets.get(kind))
                dS_drate += dS_dmetric * dmetric_drate
            demand = all_values.get(var.getName())
            rate_jac = var.rateJacobian(demand) if isinstance(
                demand, (int, float)) else 0.0
            grad[i] = dS_drate * rate_jac * var.decodeJacobian(np.array([x[i]]))
        return grad

    def _callback(self, xk: np.ndarray, convergence: float = None) -> bool:
        """
        Callback for tracking progress during optimization.

        Args:
            xk: Current best solution
            convergence: Convergence value (not used)

        Returns:
            True to stop optimization, False to continue
        """
        self._iterations += 1

        # Best-so-far is tracked inside _objectiveFunction
        if self._options['verbose']:
            print(f"line_solver.opt: iter={self._iterations}, best={self._best_value:.4f}")

        self._convergence_history.append(self._best_value)

        # Check time limit
        if time.time() >= self._deadline:
            return True  # Stop optimization

        return False

    def _buildResult(self, x: np.ndarray, objective_value: float,
                     solve_time: float) -> OptimizationResult:
        """
        Build optimization result from solution.

        Args:
            x: Optimal decision vector
            objective_value: Optimal objective value
            solve_time: Total solve time

        Returns:
            OptimizationResult
        """
        result = OptimizationResult()
        result.objective_value = objective_value
        result.variable_values = self._evaluator.decodeVariables(x)
        result.iterations = self._iterations
        result.solve_time = solve_time
        result.model_evaluations = sum(ev.getEvaluationCount()
                                       for ev in self._evaluators)
        result.convergence_history = list(self._convergence_history)

        # Evaluate constraints on every scenario (worst violation reported)
        all_values = {**self._fixed_value_dict, **result.variable_values}
        objective = self._problem.getObjective()
        all_constraints = (list(objective.getConstraints())
                           + list(self._problem.getConstraints()))

        result.feasible = True
        for evaluator, cache in zip(self._evaluators, self._caches):
            eval_result = evaluator.evaluateValuesWithCache(
                result.variable_values, cache)
            if not eval_result.feasible:
                result.feasible = False
                continue

            for constraint in all_constraints:
                violation = constraint.evaluate(eval_result, all_values)
                if violation > 0:
                    name = constraint.getName()
                    result.constraint_violations[name] = max(
                        violation, result.constraint_violations.get(name, 0.0))
                    result.feasible = False

        # Determine termination reason
        elapsed = time.time() - self._start_time
        if elapsed >= self._options['time_limit']:
            result.terminated_by = 'time_limit'
        elif self._iterations >= self._options['max_iterations']:
            result.terminated_by = 'iterations'
        else:
            result.terminated_by = 'convergence'

        return result

    def _buildEmptyResult(self) -> OptimizationResult:
        """Build result for empty problem (no variables)."""
        result = OptimizationResult()
        result.objective_value = 0.0
        result.variable_values = {}
        result.feasible = True
        result.iterations = 0
        result.solve_time = 0.0
        result.model_evaluations = 0
        result.terminated_by = 'empty'
        return result

    # Snake case aliases
    get_problem = getProblem
    get_options = getOptions
    set_option = setOption
    default_options = defaultOptions


class LineOptSolverOptions:
    """
    Configuration options for LineOptSolver.

    Provides a more structured way to configure solver options.

    Example:
        >>> options = LineOptSolverOptions()
        >>> options.setMaxIterations(200)
        >>> options.setVerbose(True)
        >>> solver = LineOptSolver(problem, **options.toDict())
    """

    def __init__(self):
        self._options = LineOptSolver.defaultOptions()

    def setStrategy(self, strategy: str) -> 'LineOptSolverOptions':
        """Set DE strategy ('best1bin', 'best2bin', 'rand1bin', etc.)."""
        self._options['strategy'] = strategy
        return self

    def setPopsize(self, popsize: int) -> 'LineOptSolverOptions':
        """Set population size multiplier."""
        self._options['popsize'] = popsize
        return self

    def setMutation(self, mutation: tuple) -> 'LineOptSolverOptions':
        """Set mutation constant (dithering range)."""
        self._options['mutation'] = mutation
        return self

    def setRecombination(self, recombination: float) -> 'LineOptSolverOptions':
        """Set crossover probability."""
        self._options['recombination'] = recombination
        return self

    def setTolerance(self, tol: float) -> 'LineOptSolverOptions':
        """Set convergence tolerance."""
        self._options['tol'] = tol
        return self

    def setMaxIterations(self, max_iter: int) -> 'LineOptSolverOptions':
        """Set maximum number of generations."""
        self._options['max_iterations'] = max_iter
        return self

    def setTimeLimit(self, time_limit: float) -> 'LineOptSolverOptions':
        """Set maximum solve time in seconds."""
        self._options['time_limit'] = time_limit
        return self

    def setSeed(self, seed: int) -> 'LineOptSolverOptions':
        """Set random seed for reproducibility."""
        self._options['seed'] = seed
        return self

    def setVerbose(self, verbose: bool) -> 'LineOptSolverOptions':
        """Set verbose output."""
        self._options['verbose'] = verbose
        return self

    def setPenaltyWeight(self, weight: float) -> 'LineOptSolverOptions':
        """Set constraint penalty weight."""
        self._options['penalty_weight'] = weight
        return self

    def setScenarioAggregation(self, aggregation: str) -> 'LineOptSolverOptions':
        """Set scenario aggregation mode ('worst' or 'mean')."""
        self._options['scenario_aggregation'] = aggregation
        return self

    def setOptimizer(self, optimizer: str) -> 'LineOptSolverOptions':
        """Set optimizer backend ('evolution', 'gradient', or 'auto')."""
        self._options['optimizer'] = optimizer
        return self

    def setLqnGradient(self, mode: str) -> 'LineOptSolverOptions':
        """Set the LQN gradient source.

        'fd' (whole-model finite difference, robust default), 'partial_sens'
        (SolverLN per-layer partial derivatives, cheap/biased), or
        'partial_plus_fd' (partial with a periodic full-FD correction).
        """
        if mode not in ('fd', 'partial_sens', 'partial_plus_fd'):
            raise ValueError("lqn_gradient must be 'fd', 'partial_sens', or "
                             "'partial_plus_fd'")
        self._options['lqn_gradient'] = mode
        return self

    def setFdRefresh(self, period: int) -> 'LineOptSolverOptions':
        """Set the partial_plus_fd full finite-difference correction period."""
        self._options['fd_refresh'] = int(period)
        return self

    def setFrozenLayers(self, layers) -> 'LineOptSolverOptions':
        """Freeze the given LQN layers (list of host/task layer names)."""
        self._options['frozen_layers'] = list(layers) if layers else None
        return self

    def toDict(self) -> dict:
        """Convert to dictionary for solver initialization."""
        return dict(self._options)

    # Snake case aliases
    set_strategy = setStrategy
    set_popsize = setPopsize
    set_mutation = setMutation
    set_recombination = setRecombination
    set_tolerance = setTolerance
    set_max_iterations = setMaxIterations
    set_time_limit = setTimeLimit
    set_seed = setSeed
    set_verbose = setVerbose
    set_penalty_weight = setPenaltyWeight
    set_scenario_aggregation = setScenarioAggregation
    set_optimizer = setOptimizer
    set_lqn_gradient = setLqnGradient
    set_fd_refresh = setFdRefresh
    set_frozen_layers = setFrozenLayers
    to_dict = toDict


# Public API
__all__ = [
    'LineOptSolver',
    'LineOptSolverOptions',
]
