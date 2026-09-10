"""
Optimization Result Containers

This module provides dataclasses for storing optimization results,
including both single-problem results and workflow decomposition results.

Key Classes:
    - EvaluationResult: Result from a single LINE model evaluation
    - OptimizationResult: Result from a single optimization run
    - WorkflowResult: Result from decomposed workflow optimization

Example:
    >>> result = solver.solve()
    >>> print(f"Objective: {result.objective_value}")
    >>> print(f"Variables: {result.variable_values}")
"""

import numpy as np
from dataclasses import dataclass, field
from typing import Dict, List, Any, Optional


@dataclass
class EvaluationResult:
    """
    Result from evaluating a LINE network model.

    Contains performance metrics extracted from SolverAUTO analysis.

    Attributes:
        feasible: Whether the model was successfully solved
        response_times: Response time per (station, class) pair
        throughputs: Throughput per (station, class) pair
        utilizations: Utilization per station
        queue_lengths: Queue length per (station, class) pair
        system_response_times: End-to-end response time per chain/class
        system_throughputs: System throughput per chain/class
        solve_time: Time spent solving the model (seconds)
        solver_used: Name of solver selected by SolverAUTO
    """
    feasible: bool = True
    response_times: Dict[tuple, float] = field(default_factory=dict)
    throughputs: Dict[tuple, float] = field(default_factory=dict)
    utilizations: Dict[str, float] = field(default_factory=dict)
    queue_lengths: Dict[tuple, float] = field(default_factory=dict)
    system_response_times: Dict[str, float] = field(default_factory=dict)
    system_throughputs: Dict[str, float] = field(default_factory=dict)
    solve_time: float = 0.0
    solver_used: str = ""
    # see _kb/05-solvers-overview.md (LineOpt: opt/results.py sensitivities) for rationale
    sensitivities: Any = None

    # camelCase views of the fields, under the names MATLAB's EvaluationResult
    # exposes, so an `opt.` script reads the result the same way in both.
    @property
    def solverUsed(self) -> str:
        """Name of the solver that answered (alias of solver_used)."""
        return self.solver_used

    @property
    def solveTime(self) -> float:
        """Time spent solving the model (alias of solve_time)."""
        return self.solve_time

    @property
    def responseTimes(self) -> Dict[tuple, float]:
        """Response time by (station, class) (alias of response_times)."""
        return self.response_times

    @property
    def queueLengths(self) -> Dict[tuple, float]:
        """Queue length by (station, class) (alias of queue_lengths)."""
        return self.queue_lengths

    @property
    def systemResponseTimes(self) -> Dict[str, float]:
        """End-to-end response time by chain (alias of system_response_times)."""
        return self.system_response_times

    @property
    def systemThroughputs(self) -> Dict[str, float]:
        """System throughput by chain (alias of system_throughputs)."""
        return self.system_throughputs

    def getResponseTime(self, station: str, jobclass: str = None) -> float:
        """
        Get response time at a station.

        Args:
            station: Station name
            jobclass: Optional job class name (returns aggregate if None)

        Returns:
            Response time value
        """
        if jobclass is not None:
            return self.response_times.get((station, jobclass), float('inf'))

        # Aggregate across classes
        total = 0.0
        count = 0
        for (s, c), rt in self.response_times.items():
            if s == station:
                total += rt
                count += 1
        return total / count if count > 0 else float('inf')

    def getThroughput(self, station: str, jobclass: str = None) -> float:
        """
        Get throughput at a station.

        Args:
            station: Station name
            jobclass: Optional job class name

        Returns:
            Throughput value
        """
        if jobclass is not None:
            return self.throughputs.get((station, jobclass), 0.0)

        # Aggregate across classes
        total = 0.0
        for (s, c), tput in self.throughputs.items():
            if s == station:
                total += tput
        return total

    def getUtilization(self, station: str) -> float:
        """
        Get utilization at a station.

        Args:
            station: Station name

        Returns:
            Utilization value (0 to 1)
        """
        return self.utilizations.get(station, 0.0)

    def getQueueLength(self, station: str, jobclass: str = None) -> float:
        """
        Get queue length at a station.

        Args:
            station: Station name
            jobclass: Optional job class name

        Returns:
            Queue length value
        """
        if jobclass is not None:
            return self.queue_lengths.get((station, jobclass), 0.0)

        # Aggregate across classes
        total = 0.0
        for (s, c), qlen in self.queue_lengths.items():
            if s == station:
                total += qlen
        return total

    def getSystemResponseTime(self, jobclass: str = None) -> float:
        """
        Get end-to-end (system) response time.

        Args:
            jobclass: Optional chain/class name. If None, returns the
                throughput-weighted average across chains (equal to total
                jobs in system divided by total throughput by Little's law).

        Returns:
            System response time (inf if unavailable)
        """
        if jobclass is not None:
            return self.system_response_times.get(jobclass, float('inf'))

        if not self.system_response_times:
            return float('inf')

        total_tput = 0.0
        weighted = 0.0
        for name, rt in self.system_response_times.items():
            tput = self.system_throughputs.get(name, 0.0)
            weighted += rt * tput
            total_tput += tput
        if total_tput > 0:
            return weighted / total_tput
        # Fall back to unweighted mean if throughputs are unavailable
        values = list(self.system_response_times.values())
        return sum(values) / len(values)

    def getSystemThroughput(self, jobclass: str = None) -> float:
        """
        Get system throughput.

        Args:
            jobclass: Optional chain/class name (sum across chains if None)

        Returns:
            System throughput value
        """
        if jobclass is not None:
            return self.system_throughputs.get(jobclass, 0.0)
        return sum(self.system_throughputs.values())

    # Snake case aliases
    get_response_time = getResponseTime
    get_throughput = getThroughput
    get_utilization = getUtilization
    get_queue_length = getQueueLength
    get_system_response_time = getSystemResponseTime
    get_system_throughput = getSystemThroughput


@dataclass
class OptimizationResult:
    """
    Result from a single optimization run.

    Contains the optimal solution found, objective value, constraint status,
    and solver statistics.

    Attributes:
        objective_value: Final objective function value
        variable_values: Dict mapping variable names to optimal values
        constraint_violations: Dict mapping constraint names to violation amounts
        feasible: Whether all constraints are satisfied
        iterations: Number of DE generations
        solve_time: Total optimization time (seconds)
        model_evaluations: Number of LINE model evaluations
        convergence_history: Objective value history over iterations
        terminated_by: Reason for termination ('convergence', 'iterations', 'time_limit')
    """
    objective_value: float = float('inf')
    variable_values: Dict[str, Any] = field(default_factory=dict)
    constraint_violations: Dict[str, float] = field(default_factory=dict)
    feasible: bool = False
    iterations: int = 0
    solve_time: float = 0.0
    model_evaluations: int = 0
    convergence_history: List[float] = field(default_factory=list)
    terminated_by: str = ""

    # camelCase views of the fields, under the names MATLAB's OptimizationResult
    # exposes, so an `opt.` script reads the result the same way in both.
    @property
    def objectiveValue(self) -> float:
        """Final objective function value (alias of objective_value)."""
        return self.objective_value

    @property
    def variableValues(self) -> Dict[str, Any]:
        """Optimal values by variable name (alias of variable_values)."""
        return self.variable_values

    @property
    def constraintViolations(self) -> Dict[str, float]:
        """Violation by constraint name (alias of constraint_violations)."""
        return self.constraint_violations

    @property
    def modelEvaluations(self) -> int:
        """Number of LINE model evaluations (alias of model_evaluations)."""
        return self.model_evaluations

    @property
    def solveTime(self) -> float:
        """Total optimization time in seconds (alias of solve_time)."""
        return self.solve_time

    @property
    def convergenceHistory(self) -> List[float]:
        """Objective value per iteration (alias of convergence_history)."""
        return self.convergence_history

    @property
    def terminatedBy(self) -> str:
        """Reason for termination (alias of terminated_by)."""
        return self.terminated_by

    def isFeasible(self) -> bool:
        """Check if solution satisfies all constraints."""
        return self.feasible

    def getObjectiveValue(self) -> float:
        """Get the objective value."""
        return self.objective_value

    def getVariableValue(self, name: str) -> Any:
        """
        Get the optimal value for a variable.

        Args:
            name: Variable name

        Returns:
            Optimal value for the variable
        """
        return self.variable_values.get(name)

    def getConstraintViolation(self, name: str) -> float:
        """
        Get constraint violation amount.

        Args:
            name: Constraint name

        Returns:
            Violation amount (0 if satisfied)
        """
        return self.constraint_violations.get(name, 0.0)

    def getTotalViolation(self) -> float:
        """Get total constraint violation."""
        return sum(self.constraint_violations.values())

    def __repr__(self) -> str:
        status = "feasible" if self.feasible else "infeasible"
        return (f"OptimizationResult(obj={self.objective_value:.4f}, "
                f"{status}, iters={self.iterations}, "
                f"evals={self.model_evaluations}, time={self.solve_time:.2f}s)")

    # Snake case aliases
    is_feasible = isFeasible
    get_objective_value = getObjectiveValue
    get_variable_value = getVariableValue
    get_constraint_violation = getConstraintViolation
    get_total_violation = getTotalViolation


@dataclass
class SubProblemResult:
    """
    Result from solving a subproblem in decomposition.

    Attributes:
        name: Subproblem name
        result: OptimizationResult for this subproblem
        variables_fixed: Variables that were fixed from previous subproblems
    """
    name: str = ""
    result: OptimizationResult = field(default_factory=OptimizationResult)
    variables_fixed: Dict[str, Any] = field(default_factory=dict)


@dataclass
class WorkflowResult:
    """
    Result from decomposed workflow optimization.

    Contains results from each subproblem and overall convergence info.

    Attributes:
        final_objective: Final objective value after all subproblems
        subproblem_results: Dict mapping subproblem names to their results
        cycles_completed: Number of decomposition cycles completed
        converged: Whether fixed-point convergence was achieved
        total_solve_time: Total time including all subproblems
        objective_history: Objective value after each cycle
        final_variable_values: Consolidated variable values from all subproblems
    """
    final_objective: float = float('inf')
    subproblem_results: Dict[str, SubProblemResult] = field(default_factory=dict)
    cycles_completed: int = 0
    converged: bool = False
    total_solve_time: float = 0.0
    objective_history: List[float] = field(default_factory=list)
    final_variable_values: Dict[str, Any] = field(default_factory=dict)

    def isConverged(self) -> bool:
        """Check if fixed-point convergence was achieved."""
        return self.converged

    def getSubProblemResult(self, name: str) -> Optional[SubProblemResult]:
        """
        Get result for a specific subproblem.

        Args:
            name: Subproblem name

        Returns:
            SubProblemResult or None if not found
        """
        return self.subproblem_results.get(name)

    def getFinalVariableValue(self, name: str) -> Any:
        """
        Get final value for a variable across all subproblems.

        Args:
            name: Variable name

        Returns:
            Final optimal value
        """
        return self.final_variable_values.get(name)

    def __repr__(self) -> str:
        status = "converged" if self.converged else "not converged"
        n_subs = len(self.subproblem_results)
        return (f"WorkflowResult(obj={self.final_objective:.4f}, "
                f"{n_subs} subproblems, {self.cycles_completed} cycles, "
                f"{status}, time={self.total_solve_time:.2f}s)")

    # Snake case aliases
    is_converged = isConverged
    get_subproblem_result = getSubProblemResult
    get_final_variable_value = getFinalVariableValue


# Public API
__all__ = [
    'EvaluationResult',
    'OptimizationResult',
    'SubProblemResult',
    'WorkflowResult',
]
