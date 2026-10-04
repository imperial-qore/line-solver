"""
Exact Sizing Solver for Monotone Single-Variable Problems

This module provides BisectionSolver, an exact O(log n) alternative to
differential evolution for the common capacity sizing pattern: a single
integer decision variable whose feasibility is monotone (e.g. adding servers
can only reduce response time and utilization).

Key Classes:
    - BisectionSolver: Bisection over a single integer variable

Example:
    >>> solver = BisectionSolver(problem)
    >>> result = solver.solve()
    >>> print(result.variable_values)
"""

import time
import numpy as np
from typing import Dict, Any, List, Optional

from .evaluator import LineEvaluator
from .results import OptimizationResult, EvaluationResult


class BisectionSolver:
    """
    Exact solver for single integer-variable problems with monotone
    feasibility.

    Applicable when the problem has exactly one decision variable of
    dimension 1 with an integer domain (ServerAllocation, StationReplicas,
    JobPopulation) and the feasible set is an interval touching one bound:

    - direction='min_feasible' (default): feasibility is monotone
      non-decreasing in the variable (more servers -> constraints easier)
      and cost is non-decreasing, so the optimum is the smallest feasible
      value. This is the standard server sizing pattern.
    - direction='max_feasible': feasibility is monotone non-increasing
      (more jobs -> constraints harder) and the objective rewards larger
      values, so the optimum is the largest feasible value. This is the
      standard population sizing pattern.

    Each probe is one LINE solve, so the search costs O(log(hi - lo))
    evaluations against O(popsize * iterations) for differential evolution.
    Constraints are enforced on the base model and all scenarios.

    Args:
        problem: OptimizationProblem with a single integer variable
        **options: Solver options (see defaultOptions)

    Example:
        >>> solver = BisectionSolver(problem)
        >>> result = solver.solve()
    """

    @staticmethod
    def defaultOptions() -> dict:
        """Get default solver options."""
        return {
            'direction': 'min_feasible',
            'verbose': False,
        }

    def __init__(self, problem: 'OptimizationProblem', **options):
        """
        Initialize the solver.

        Args:
            problem: OptimizationProblem with exactly one integer variable
            **options: Override default options

        Raises:
            ValueError: If the problem does not have exactly one
                dimension-1 variable
        """
        self._problem = problem
        self._options = {**self.defaultOptions(), **options}

        variables = problem.getVariables()
        if len(variables) != 1 or variables[0].getDimension() != 1:
            raise ValueError("BisectionSolver requires exactly one "
                             "dimension-1 decision variable")
        self._variable = variables[0]

        lo = self._variable.decode(np.array([0.0]))
        hi = self._variable.decode(np.array([1.0]))
        if not isinstance(lo, (int, np.integer)) or not isinstance(hi, (int, np.integer)):
            raise ValueError("BisectionSolver requires an integer-valued "
                             "variable (e.g. ServerAllocation)")
        self._lo, self._hi = int(lo), int(hi)

        fixed = problem.getFixedVariables() if hasattr(problem, 'getFixedVariables') else []
        self._fixed_value_dict = {var.getName(): value for var, value in fixed}

        # Evaluators: base model plus one per workload scenario
        self._evaluators = [LineEvaluator(problem.getModel(), variables, fixed)]
        if hasattr(problem, 'getScenarios'):
            for scenario_model, _ in problem.getScenarios():
                self._evaluators.append(LineEvaluator(
                    scenario_model, variables, fixed))

        self._probe_cache: Dict[int, tuple] = {}

    def getProblem(self) -> 'OptimizationProblem':
        """Get the optimization problem."""
        return self._problem

    def _allConstraints(self) -> list:
        """All constraints: objective-level plus problem-level."""
        objective = self._problem.getObjective()
        constraints = list(self._problem.getConstraints())
        if objective is not None:
            constraints = list(objective.getConstraints()) + constraints
        return constraints

    def _probe(self, value: int) -> tuple:
        """
        Evaluate feasibility at an integer variable value.

        Returns:
            (feasible, violations dict, base EvaluationResult)
        """
        if value in self._probe_cache:
            return self._probe_cache[value]

        values = {self._variable.getName(): value}
        all_values = {**self._fixed_value_dict, **values}
        constraints = self._allConstraints()

        feasible = True
        violations: Dict[str, float] = {}
        base_result: Optional[EvaluationResult] = None

        for evaluator in self._evaluators:
            eval_result = evaluator.evaluateValues(values)
            if base_result is None:
                base_result = eval_result
            if not eval_result.feasible:
                feasible = False
                continue
            for constraint in constraints:
                violation = constraint.evaluate(eval_result, all_values)
                if violation > 0:
                    feasible = False
                    violations[constraint.getName()] = max(
                        violation, violations.get(constraint.getName(), 0.0))

        outcome = (feasible, violations, base_result)
        self._probe_cache[value] = outcome
        return outcome

    def solve(self) -> OptimizationResult:
        """
        Run the bisection search.

        Returns:
            OptimizationResult. If no value in the domain is feasible, the
            result is marked infeasible and reports the boundary value with
            the smallest violation (the domain bound in the search
            direction).
        """
        start_time = time.time()
        direction = self._options['direction']
        verbose = self._options['verbose']
        lo, hi = self._lo, self._hi
        iterations = 0

        if direction == 'min_feasible':
            # Invariant: feasibility is non-decreasing; find smallest feasible
            while lo < hi:
                mid = (lo + hi) // 2
                feasible, _, _ = self._probe(mid)
                iterations += 1
                if verbose:
                    print(f"line_solver.opt: bisection probe {self._variable.getName()}"
                          f"={mid} -> {'feasible' if feasible else 'infeasible'}")
                if feasible:
                    hi = mid
                else:
                    lo = mid + 1
        elif direction == 'max_feasible':
            # Invariant: feasibility is non-increasing; find largest feasible
            while lo < hi:
                mid = (lo + hi + 1) // 2
                feasible, _, _ = self._probe(mid)
                iterations += 1
                if verbose:
                    print(f"line_solver.opt: bisection probe {self._variable.getName()}"
                          f"={mid} -> {'feasible' if feasible else 'infeasible'}")
                if feasible:
                    lo = mid
                else:
                    hi = mid - 1
        else:
            raise ValueError(f"Unknown direction: {direction}")

        chosen = lo
        feasible, violations, base_result = self._probe(chosen)

        result = OptimizationResult()
        result.variable_values = {self._variable.getName(): chosen}
        result.feasible = feasible
        result.constraint_violations = dict(violations)
        result.iterations = iterations
        result.model_evaluations = sum(ev.getEvaluationCount()
                                       for ev in self._evaluators)
        result.terminated_by = 'bisection'

        objective = self._problem.getObjective()
        if objective is not None and base_result is not None:
            all_values = {**self._fixed_value_dict, **result.variable_values}
            result.objective_value = objective.evaluate(base_result, all_values)

        result.solve_time = time.time() - start_time
        return result

    # Snake case aliases
    get_problem = getProblem
    default_options = defaultOptions


# Public API
__all__ = ['BisectionSolver']
