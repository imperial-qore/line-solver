"""
Pareto Frontier Computation via Epsilon-Constraint Sweeps

This module provides ParetoSweep, which characterizes the cost-performance
tradeoff by re-solving the problem under a sweep of constraint bounds
(epsilon-constraint method) and filtering the results to the non-dominated
frontier.

Key Classes:
    - ParetoPoint: One point of the tradeoff curve
    - ParetoSweep: Epsilon-constraint sweep driver

Example:
    >>> sweep = ParetoSweep(
    ...     problem,
    ...     constraint_factory=lambda eps: ResponseTimeConstraint(
    ...         queue, jobclass, max_value=eps),
    ...     epsilons=[0.5, 1.0, 2.0],
    ... )
    >>> points = sweep.solve(seed=42)
    >>> frontier = sweep.getFrontier()
"""

import time
from dataclasses import dataclass, field
from typing import Callable, List, Any, Optional

from .objectives import Constraint
from .results import OptimizationResult


@dataclass
class ParetoPoint:
    """
    One point of the cost-performance tradeoff curve.

    Attributes:
        epsilon: The constraint bound used for this solve
        objective_value: Optimal objective at this bound
        feasible: Whether a feasible solution was found
        result: Full OptimizationResult for this solve
    """
    epsilon: float = 0.0
    objective_value: float = float('inf')
    feasible: bool = False
    result: OptimizationResult = field(default_factory=OptimizationResult)


class ParetoSweep:
    """
    Epsilon-constraint sweep for bi-objective tradeoff analysis.

    Solves the problem once per epsilon value, each time adding the
    constraint produced by constraint_factory(epsilon) to the problem's
    constraints. The typical use is sweeping an SLA bound (e.g. maximum
    response time) to obtain the cost frontier.

    The base problem is not modified: each solve uses a shallow clone
    sharing the model, variables, objective, fixed variables, and scenarios.

    Args:
        problem: Base OptimizationProblem (objective and variables set)
        constraint_factory: Callable mapping an epsilon value to a Constraint
        epsilons: List of epsilon values to sweep
        solver: 'de' for differential evolution (default) or 'bisection'
            for BisectionSolver (single integer variable problems)

    Example:
        >>> sweep = ParetoSweep(problem, factory, [0.5, 1.0, 2.0])
        >>> points = sweep.solve(seed=42)
    """

    def __init__(self, problem: 'OptimizationProblem',
                 constraint_factory: Callable[[float], Constraint],
                 epsilons: List[float],
                 solver: str = 'de'):
        self._problem = problem
        self._constraint_factory = constraint_factory
        self._epsilons = list(epsilons)
        self._solver = solver
        self._points: List[ParetoPoint] = []

    def getPoints(self) -> List[ParetoPoint]:
        """Get all sweep points from the last solve() call."""
        return list(self._points)

    def _cloneProblem(self, epsilon: float) -> 'OptimizationProblem':
        """Clone the base problem with the epsilon constraint added."""
        from .problem import OptimizationProblem

        clone = OptimizationProblem(self._problem.getModel())
        for var in self._problem.getVariables():
            clone.addVariable(var)
        clone.setObjective(self._problem.getObjective())
        for constraint in self._problem.getConstraints():
            clone.addConstraint(constraint)
        clone.setFixedVariables(self._problem.getFixedVariables())
        for scenario_model, weight in self._problem.getScenarios():
            clone.addScenario(scenario_model, weight)

        clone.addConstraint(self._constraint_factory(epsilon))
        return clone

    def solve(self, **solver_options) -> List[ParetoPoint]:
        """
        Run the sweep.

        Args:
            **solver_options: Options passed to each solve (e.g. seed,
                max_iterations, or direction for the bisection solver)

        Returns:
            List of ParetoPoint, one per epsilon, in sweep order
        """
        self._points = []
        for epsilon in self._epsilons:
            clone = self._cloneProblem(epsilon)

            if self._solver == 'bisection':
                from .sizing import BisectionSolver
                result = BisectionSolver(clone, **solver_options).solve()
            else:
                result = clone.solve(**solver_options)

            self._points.append(ParetoPoint(
                epsilon=epsilon,
                objective_value=result.objective_value,
                feasible=result.feasible,
                result=result,
            ))
        return list(self._points)

    def getFrontier(self) -> List[ParetoPoint]:
        """
        Get the non-dominated frontier of the sweep.

        A feasible point dominates another if it is at least as good in
        both coordinates (objective value and epsilon, both minimized) and
        strictly better in one. Infeasible points are excluded.

        Returns:
            Non-dominated ParetoPoints sorted by increasing epsilon
        """
        feasible = [p for p in self._points if p.feasible]
        frontier = []
        for p in feasible:
            dominated = any(
                (q.objective_value <= p.objective_value
                 and q.epsilon <= p.epsilon
                 and (q.objective_value < p.objective_value
                      or q.epsilon < p.epsilon))
                for q in feasible)
            if not dominated:
                frontier.append(p)
        return sorted(frontier, key=lambda p: p.epsilon)

    def plot(self, ax=None, xlabel: str = 'epsilon',
             ylabel: str = 'objective value',
             title: str = 'Pareto frontier',
             show_dominated: bool = True,
             save_path: Optional[str] = None):
        """
        Plot the actual (non-dominated) frontier of the sweep.

        Only the non-dominated points returned by getFrontier() are drawn as
        the frontier, and they are connected as a piecewise-constant staircase
        (matplotlib step, where='post') rather than a straight-line
        interpolation: the epsilon-constraint cost frontier is constant between
        the sampled epsilon breakpoints, so a diagonal line would misrepresent
        it. Feasible but dominated sample points are shown faintly for context.

        Args:
            ax: Existing matplotlib Axes to draw on; a new one is created if None
            xlabel, ylabel, title: Axis labels and title
            show_dominated: Also mark feasible dominated sample points
            save_path: If given, save the figure to this path

        Returns:
            The matplotlib Axes the frontier was drawn on
        """
        import matplotlib.pyplot as plt

        frontier = self.getFrontier()
        if not frontier:
            raise ValueError(
                "no feasible points to plot; call solve() first")

        if ax is None:
            _, ax = plt.subplots()

        xs = [p.epsilon for p in frontier]
        ys = [p.objective_value for p in frontier]
        # Actual frontier: staircase through the non-dominated points only.
        ax.step(xs, ys, where='post', color='C0', linewidth=1.8,
                zorder=2, label='Pareto frontier')
        ax.scatter(xs, ys, color='C0', zorder=3)

        if show_dominated:
            frontier_eps = {p.epsilon for p in frontier}
            dominated = [p for p in self._points
                         if p.feasible and p.epsilon not in frontier_eps]
            if dominated:
                ax.scatter([p.epsilon for p in dominated],
                           [p.objective_value for p in dominated],
                           marker='x', color='0.6', zorder=1,
                           label='dominated')

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)

        if save_path is not None:
            ax.figure.savefig(save_path, dpi=120, bbox_inches='tight')
        return ax

    # Snake case aliases
    get_points = getPoints
    get_frontier = getFrontier


# Public API
__all__ = ['ParetoPoint', 'ParetoSweep']
