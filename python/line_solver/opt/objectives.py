"""
Objective Functions and Constraints for Network Optimization

This module provides objective functions and constraint classes for
queueing network optimization problems.

Key Classes:
    Objectives:
        - Objective: Abstract base class
        - MinimizeCost: Minimize cost subject to SLA constraints
        - MaximizePerformance: Maximize performance subject to budget

    Constraints:
        - Constraint: Abstract base class
        - ResponseTimeConstraint: Response time <= max_value
        - ThroughputConstraint: Throughput >= min_value
        - UtilizationConstraint: Utilization <= max_value
        - BudgetConstraint: Total cost <= budget

Example:
    >>> obj = MinimizeCost(
    ...     server_cost={queue: 100.0},
    ...     subject_to=[ResponseTimeConstraint(queue, max_value=1.0)]
    ... )
"""

import numpy as np
from abc import ABC, abstractmethod
from typing import Dict, List, Optional, Any, Union
from dataclasses import dataclass, field

from .results import EvaluationResult


def _upper_bound_violation(actual: float, bound: float) -> float:
    """
    Violation for an upper-bound constraint (actual <= bound).

    A non-finite metric (NaN/inf) signals an unstable or numerically
    diverged configuration; it must be reported as strongly infeasible
    rather than satisfied. Note ``max(0.0, float('nan') - bound)`` returns
    0.0 in Python, so a NaN metric would otherwise silently satisfy the
    constraint and mislead the optimizer into treating a divergent design
    as feasible.
    """
    if not np.isfinite(actual):
        return float('inf')
    return max(0.0, actual - bound)


def _lower_bound_violation(actual: float, bound: float) -> float:
    """
    Violation for a lower-bound constraint (actual >= bound).

    A non-finite metric (NaN/inf) is treated as strongly infeasible for the
    same reason as :func:`_upper_bound_violation`.
    """
    if not np.isfinite(actual):
        return float('inf')
    return max(0.0, bound - actual)


class Constraint(ABC):
    """
    Abstract base class for optimization constraints.

    Constraints define limits on performance metrics that must be satisfied.
    Each constraint can compute its violation amount given model metrics.
    """

    def __init__(self, name: str = None):
        """
        Initialize constraint.

        Args:
            name: Optional constraint name (auto-generated if None)
        """
        self._name = name or self._generateName()

    def getName(self) -> str:
        """Get constraint name."""
        return self._name

    @abstractmethod
    def _generateName(self) -> str:
        """Generate a default name for this constraint."""
        raise NotImplementedError

    @abstractmethod
    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate constraint violation.

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values

        Returns:
            Violation amount (0 if satisfied, positive if violated)
        """
        raise NotImplementedError

    def isSatisfied(self, result: EvaluationResult,
                    variable_values: Dict[str, Any],
                    tolerance: float = 1e-6) -> bool:
        """
        Check if constraint is satisfied.

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values
            tolerance: Tolerance for numerical comparison

        Returns:
            True if satisfied
        """
        return self.evaluate(result, variable_values) <= tolerance

    # Snake case aliases
    get_name = getName
    is_satisfied = isSatisfied


class ResponseTimeConstraint(Constraint):
    """
    Response time constraint: RT <= max_value.

    Args:
        station: Station or station name to constrain
        jobclass: Optional JobClass or class name (aggregate if None)
        max_value: Maximum allowed response time
        name: Optional constraint name

    Example:
        >>> constraint = ResponseTimeConstraint(queue, jobclass, max_value=2.0)
    """

    def __init__(self, station: Union[str, 'Station'],
                 jobclass: Union[str, 'JobClass'] = None,
                 max_value: float = None,
                 name: str = None):
        self._station = station.getName() if hasattr(station, 'getName') else station
        self._jobclass = None
        if jobclass is not None:
            self._jobclass = jobclass.getName() if hasattr(jobclass, 'getName') else jobclass
        self._max_value = max_value
        super().__init__(name)

    def _generateName(self) -> str:
        if self._jobclass:
            return f"RT_{self._station}_{self._jobclass}_le_{self._max_value}"
        return f"RT_{self._station}_le_{self._max_value}"

    def getStation(self) -> str:
        """Get station name."""
        return self._station

    def getJobClass(self) -> Optional[str]:
        """Get job class name."""
        return self._jobclass

    def getMaxValue(self) -> float:
        """Get maximum allowed response time."""
        return self._max_value

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate response time constraint violation.

        Returns max(0, actual_RT - max_value).
        """
        actual = result.getResponseTime(self._station, self._jobclass)
        return _upper_bound_violation(actual, self._max_value)

    get_station = getStation
    get_jobclass = getJobClass
    get_max_value = getMaxValue


class SystemResponseTimeConstraint(Constraint):
    """
    End-to-end response time constraint: SysRespT <= max_value.

    Constrains the system (chain-level) response time reported by LINE,
    i.e. the full sojourn time across the job's path, rather than the
    per-station response time of ResponseTimeConstraint.

    Args:
        jobclass: Optional JobClass or class name (throughput-weighted
            aggregate across chains if None)
        max_value: Maximum allowed system response time
        name: Optional constraint name

    Example:
        >>> constraint = SystemResponseTimeConstraint(jobclass, max_value=2.0)
    """

    def __init__(self, jobclass: Union[str, 'JobClass'] = None,
                 max_value: float = None,
                 name: str = None):
        self._jobclass = None
        if jobclass is not None:
            self._jobclass = jobclass.getName() if hasattr(jobclass, 'getName') else jobclass
        self._max_value = max_value
        super().__init__(name)

    def _generateName(self) -> str:
        if self._jobclass:
            return f"SysRT_{self._jobclass}_le_{self._max_value}"
        return f"SysRT_le_{self._max_value}"

    def getJobClass(self) -> Optional[str]:
        """Get job class name."""
        return self._jobclass

    def getMaxValue(self) -> float:
        """Get maximum allowed system response time."""
        return self._max_value

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate system response time constraint violation.

        Returns max(0, actual_sys_RT - max_value).
        """
        actual = result.getSystemResponseTime(self._jobclass)
        return _upper_bound_violation(actual, self._max_value)

    get_jobclass = getJobClass
    get_max_value = getMaxValue


class ThroughputConstraint(Constraint):
    """
    Throughput constraint: Tput >= min_value.

    Args:
        station: Station or station name to constrain
        jobclass: Optional JobClass (aggregate if None)
        min_value: Minimum required throughput
        name: Optional constraint name

    Example:
        >>> constraint = ThroughputConstraint(queue, min_value=10.0)
    """

    def __init__(self, station: Union[str, 'Station'],
                 jobclass: Union[str, 'JobClass'] = None,
                 min_value: float = None,
                 name: str = None):
        self._station = station.getName() if hasattr(station, 'getName') else station
        self._jobclass = None
        if jobclass is not None:
            self._jobclass = jobclass.getName() if hasattr(jobclass, 'getName') else jobclass
        self._min_value = min_value
        super().__init__(name)

    def _generateName(self) -> str:
        if self._jobclass:
            return f"Tput_{self._station}_{self._jobclass}_ge_{self._min_value}"
        return f"Tput_{self._station}_ge_{self._min_value}"

    def getStation(self) -> str:
        """Get station name."""
        return self._station

    def getMinValue(self) -> float:
        """Get minimum required throughput."""
        return self._min_value

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate throughput constraint violation.

        Returns max(0, min_value - actual_tput).
        """
        actual = result.getThroughput(self._station, self._jobclass)
        return _lower_bound_violation(actual, self._min_value)

    get_station = getStation
    get_min_value = getMinValue


class UtilizationConstraint(Constraint):
    """
    Utilization constraint: U <= max_value.

    Args:
        station: Station or station name to constrain
        max_value: Maximum allowed utilization (0 to 1)
        name: Optional constraint name

    Example:
        >>> constraint = UtilizationConstraint(queue, max_value=0.8)
    """

    def __init__(self, station: Union[str, 'Station'],
                 max_value: float = None,
                 name: str = None):
        self._station = station.getName() if hasattr(station, 'getName') else station
        self._max_value = max_value
        super().__init__(name)

    def _generateName(self) -> str:
        return f"Util_{self._station}_le_{self._max_value}"

    def getStation(self) -> str:
        """Get station name."""
        return self._station

    def getMaxValue(self) -> float:
        """Get maximum allowed utilization."""
        return self._max_value

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate utilization constraint violation.

        Returns max(0, actual_util - max_value).
        """
        actual = result.getUtilization(self._station)
        return _upper_bound_violation(actual, self._max_value)

    get_station = getStation
    get_max_value = getMaxValue


class BudgetConstraint(Constraint):
    """
    Budget constraint: total cost <= budget.

    Cost is computed from variable values using cost coefficients.

    Args:
        budget: Maximum allowed total cost
        cost_coefficients: Dict mapping variable names to cost per unit
        name: Optional constraint name

    Example:
        >>> constraint = BudgetConstraint(
        ...     budget=500.0,
        ...     cost_coefficients={'Queue_servers': 50.0, 'Queue2_servers': 75.0}
        ... )
    """

    def __init__(self, budget: float,
                 cost_coefficients: Dict[str, float] = None,
                 name: str = None):
        self._budget = budget
        self._cost_coefficients = cost_coefficients or {}
        super().__init__(name)

    def _generateName(self) -> str:
        return f"Budget_le_{self._budget}"

    def getBudget(self) -> float:
        """Get budget limit."""
        return self._budget

    def getCostCoefficients(self) -> Dict[str, float]:
        """Get cost coefficients."""
        return self._cost_coefficients

    def computeCost(self, variable_values: Dict[str, Any]) -> float:
        """
        Compute total cost from variable values.

        Args:
            variable_values: Dict mapping variable names to values

        Returns:
            Total cost
        """
        total = 0.0
        for var_name, coeff in self._cost_coefficients.items():
            value = variable_values.get(var_name, 0)
            if isinstance(value, (int, float)):
                total += coeff * value
            elif isinstance(value, np.ndarray):
                total += coeff * np.sum(value)
        return total

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate budget constraint violation.

        Returns max(0, total_cost - budget).
        """
        cost = self.computeCost(variable_values)
        return max(0.0, cost - self._budget)

    get_budget = getBudget
    get_cost_coefficients = getCostCoefficients
    compute_cost = computeCost


class Objective(ABC):
    """
    Abstract base class for optimization objectives.

    An objective defines what to minimize or maximize, potentially
    subject to constraints.
    """

    def __init__(self):
        self._constraints: List[Constraint] = []

    def getConstraints(self) -> List[Constraint]:
        """Get constraints associated with this objective."""
        return self._constraints

    @abstractmethod
    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate objective function value.

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values

        Returns:
            Objective value (to be minimized)
        """
        raise NotImplementedError

    @abstractmethod
    def isMinimization(self) -> bool:
        """Check if this is a minimization objective."""
        raise NotImplementedError

    def evaluateWithPenalty(self, result: EvaluationResult,
                            variable_values: Dict[str, Any],
                            penalty_weight: float = 1e6) -> float:
        """
        Evaluate objective with constraint penalty.

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values
            penalty_weight: Multiplier for constraint violations

        Returns:
            Objective value + penalty for violations
        """
        obj = self.evaluate(result, variable_values)

        # Add penalty for constraint violations
        penalty = 0.0
        for constraint in self._constraints:
            violation = constraint.evaluate(result, variable_values)
            penalty += violation * penalty_weight

        return obj + penalty

    # Snake case aliases (use methods to avoid abstract method issues)
    def get_constraints(self) -> List[Constraint]:
        return self.getConstraints()

    def is_minimization(self) -> bool:
        return self.isMinimization()

    def evaluate_with_penalty(self, result: EvaluationResult,
                              variable_values: Dict[str, Any],
                              penalty_weight: float = 1e6) -> float:
        return self.evaluateWithPenalty(result, variable_values, penalty_weight)


class MinimizeCost(Objective):
    """
    Minimize cost subject to service-level constraints.

    Cost is computed as a weighted sum of:
    - Server costs: sum(server_cost[station] * num_servers)
    - Rate costs: sum(rate_cost[station] * rate)
    - Replica costs: sum(replica_cost[station] * replicas)

    Args:
        server_cost: Dict mapping stations to cost per server
        rate_cost: Dict mapping stations to cost per unit rate
        replica_cost: Dict mapping stations to cost per replica
        subject_to: List of constraints (SLA constraints)

    Example:
        >>> obj = MinimizeCost(
        ...     server_cost={queue: 100.0},
        ...     subject_to=[ResponseTimeConstraint(queue, max_value=1.0)]
        ... )
    """

    def __init__(self,
                 server_cost: Dict[Any, float] = None,
                 rate_cost: Dict[Any, float] = None,
                 replica_cost: Dict[Any, float] = None,
                 subject_to: List[Constraint] = None):
        super().__init__()

        # Convert station objects to names
        self._server_cost = {}
        if server_cost:
            for k, v in server_cost.items():
                name = k.getName() if hasattr(k, 'getName') else str(k)
                self._server_cost[f"{name}_servers"] = v

        self._rate_cost = {}
        if rate_cost:
            for k, v in rate_cost.items():
                name = k.getName() if hasattr(k, 'getName') else str(k)
                self._rate_cost[name] = v

        self._replica_cost = {}
        if replica_cost:
            for k, v in replica_cost.items():
                name = k.getName() if hasattr(k, 'getName') else str(k)
                self._replica_cost[f"{name}_replicas"] = v

        if subject_to:
            self._constraints = list(subject_to)

    def getServerCost(self) -> Dict[str, float]:
        """Get server cost coefficients."""
        return self._server_cost

    def getRateCost(self) -> Dict[str, float]:
        """Get rate cost coefficients."""
        return self._rate_cost

    def getReplicaCost(self) -> Dict[str, float]:
        """Get replica cost coefficients."""
        return self._replica_cost

    def isMinimization(self) -> bool:
        return True

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate total cost.

        Args:
            result: EvaluationResult (not used directly for cost)
            variable_values: Current variable values

        Returns:
            Total cost
        """
        total = 0.0

        # Server costs
        for var_name, cost in self._server_cost.items():
            value = variable_values.get(var_name, 0)
            if isinstance(value, (int, float)):
                total += cost * value

        # Rate costs
        for var_pattern, cost in self._rate_cost.items():
            for var_name, value in variable_values.items():
                if var_pattern in var_name and 'rate' in var_name.lower():
                    if isinstance(value, (int, float)):
                        total += cost * value

        # Replica costs
        for var_name, cost in self._replica_cost.items():
            value = variable_values.get(var_name, 0)
            if isinstance(value, (int, float)):
                total += cost * value

        return total

    get_server_cost = getServerCost
    get_rate_cost = getRateCost
    get_replica_cost = getReplicaCost


class MaximizePerformance(Objective):
    """
    Maximize performance subject to budget constraints.

    Performance is a weighted combination of metrics:
    - Throughput (maximize)
    - 1/Response time (maximize, i.e., minimize response time)
    - 1/Queue length (maximize)

    Since differential evolution minimizes, we negate the performance value.

    Args:
        throughput_weight: Weight for throughput in objective
        response_time_weight: Weight for 1/response_time in objective
        queue_length_weight: Weight for 1/queue_length in objective
        stations: List of stations to include (all if None)
        budget: Optional budget constraint
        budget_terms: Dict mapping variable names to cost per unit

    Example:
        >>> obj = MaximizePerformance(
        ...     throughput_weight=1.0,
        ...     response_time_weight=2.0,
        ...     budget=500.0,
        ...     budget_terms={'Queue_servers': 50.0}
        ... )
    """

    def __init__(self,
                 throughput_weight: float = 1.0,
                 response_time_weight: float = 1.0,
                 queue_length_weight: float = 0.0,
                 stations: List[Any] = None,
                 budget: float = None,
                 budget_terms: Dict[str, float] = None):
        super().__init__()

        self._throughput_weight = throughput_weight
        self._response_time_weight = response_time_weight
        self._queue_length_weight = queue_length_weight

        self._stations = None
        if stations:
            self._stations = [s.getName() if hasattr(s, 'getName') else str(s)
                              for s in stations]

        # Add budget constraint if specified
        if budget is not None:
            self._constraints.append(BudgetConstraint(
                budget=budget,
                cost_coefficients=budget_terms or {}
            ))

    def getThroughputWeight(self) -> float:
        """Get throughput weight."""
        return self._throughput_weight

    def getResponseTimeWeight(self) -> float:
        """Get response time weight."""
        return self._response_time_weight

    def isMinimization(self) -> bool:
        # We minimize negative performance (to maximize performance)
        return True

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate negative performance (for minimization).

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values

        Returns:
            Negative performance value
        """
        performance = 0.0

        # Aggregate metrics across stations
        stations = self._stations
        if stations is None:
            # Use all stations from result
            stations = set(s for s, _ in result.throughputs.keys())

        for station in stations:
            # Throughput contribution (higher is better)
            if self._throughput_weight > 0:
                tput = result.getThroughput(station)
                performance += self._throughput_weight * tput

            # Response time contribution (lower is better, so use 1/RT)
            if self._response_time_weight > 0:
                rt = result.getResponseTime(station)
                if rt > 0 and rt < float('inf'):
                    performance += self._response_time_weight * (1.0 / rt)

            # Queue length contribution (lower is better)
            if self._queue_length_weight > 0:
                qlen = result.getQueueLength(station)
                if qlen > 0:
                    performance += self._queue_length_weight * (1.0 / qlen)

        # Return negative because we're minimizing
        return -performance

    get_throughput_weight = getThroughputWeight
    get_response_time_weight = getResponseTimeWeight


class MinimizeSystemResponseTime(Objective):
    """
    Minimize the end-to-end (system) response time.

    Uses the chain-level sojourn time from LINE's system table, i.e. the
    full end-to-end response time. Typical uses are load balancing and
    routing optimization, where cost plays no role and the goal is purely
    to minimize latency, optionally under constraints.

    Args:
        jobclass: Optional JobClass or class name (throughput-weighted
            aggregate across chains if None)
        subject_to: Optional list of constraints

    Example:
        >>> obj = MinimizeSystemResponseTime(jobclass)
    """

    def __init__(self, jobclass: Union[str, 'JobClass'] = None,
                 subject_to: List[Constraint] = None):
        super().__init__()
        self._jobclass = None
        if jobclass is not None:
            self._jobclass = jobclass.getName() if hasattr(jobclass, 'getName') else jobclass
        if subject_to:
            self._constraints = list(subject_to)

    def getJobClass(self) -> Optional[str]:
        """Get job class name."""
        return self._jobclass

    def isMinimization(self) -> bool:
        return True

    def evaluate(self, result: EvaluationResult,
                 variable_values: Dict[str, Any]) -> float:
        """
        Evaluate the system response time.

        Args:
            result: EvaluationResult from model evaluation
            variable_values: Current variable values (unused)

        Returns:
            System response time (inf if unavailable)
        """
        return result.getSystemResponseTime(self._jobclass)

    get_jobclass = getJobClass


# Public API
__all__ = [
    # Constraints
    'Constraint',
    'ResponseTimeConstraint',
    'SystemResponseTimeConstraint',
    'ThroughputConstraint',
    'UtilizationConstraint',
    'BudgetConstraint',
    # Objectives
    'Objective',
    'MinimizeCost',
    'MaximizePerformance',
    'MinimizeSystemResponseTime',
]
