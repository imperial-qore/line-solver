"""
Optimization Problem Specification

This module provides the main OptimizationProblem class for specifying
queueing network optimization problems.

Key Classes:
    - OptimizationProblem: Main problem specification class

Example:
    >>> problem = OptimizationProblem(model)
    >>> problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
    >>> problem.set_objective(MinimizeCost(server_cost={queue: 100.0}))
    >>> result = problem.solve()
"""

from typing import List, Dict, Any, Optional

from .variables import DecisionVariable
from .objectives import Objective, Constraint
from .results import OptimizationResult
from .solver import LineOptSolver


class OptimizationProblem:
    """
    Main class for specifying network optimization problems.

    An OptimizationProblem combines:
    - A LINE Network model
    - Decision variables (what to optimize)
    - An objective function (what to minimize/maximize)
    - Constraints (limits on metrics)

    The problem can be solved directly or decomposed into subproblems.

    Args:
        model: LINE Network model

    Example:
        >>> from line_solver import Network, Queue, Source, Sink, OpenClass, Exp
        >>> from line_solver import OptimizationProblem, ServerAllocation, MinimizeCost
        >>>
        >>> model = Network("Test")
        >>> source = Source(model, "Source")
        >>> queue = Queue(model, "Queue", SchedStrategy.FCFS)
        >>> sink = Sink(model, "Sink")
        >>> jobclass = OpenClass(model, "Class1")
        >>> # ... configure model ...
        >>>
        >>> problem = OptimizationProblem(model)
        >>> problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        >>> problem.set_objective(MinimizeCost(server_cost={queue: 100.0}))
        >>> result = problem.solve(verbose=True)
    """

    def __init__(self, model: 'Network'):
        """
        Initialize optimization problem.

        Args:
            model: LINE Network model to optimize
        """
        self._model = model
        self._variables: List[DecisionVariable] = []
        self._objective: Optional[Objective] = None
        self._constraints: List[Constraint] = []
        # (DecisionVariable, value) pairs held constant during optimization
        self._fixed_variables: List[tuple] = []
        # (model, weight) pairs for scenario-based robust optimization
        self._scenarios: List[tuple] = []
        # LayeredNetwork (LQN) models take the SolverLN evaluation path and
        # accept LQN decision variables only (see opt/layered.py).
        from .layered import is_layered
        self._is_layered = is_layered(model)

    # Decision-variable types that operate on a LayeredNetwork vs a flat Network.
    _LQN_VAR_TYPES = frozenset({'host_demand', 'think_time',
                                'task_multiplicity', 'task_replication',
                                'processor_multiplicity'})
    _FLAT_VAR_TYPES = frozenset({'server_allocation', 'station_replicas',
                                 'service_rate', 'job_population',
                                 'class_priority', 'routing', 'class_mapping'})

    def isLayered(self) -> bool:
        """True if the model is a LayeredNetwork (LQN)."""
        return self._is_layered

    def getModel(self) -> 'Network':
        """Get the network model."""
        return self._model

    def getVariables(self) -> List[DecisionVariable]:
        """Get the list of decision variables."""
        return list(self._variables)

    def getObjective(self) -> Optional[Objective]:
        """Get the objective function."""
        return self._objective

    def getConstraints(self) -> List[Constraint]:
        """Get the list of constraints."""
        return list(self._constraints)

    def addVariable(self, variable: DecisionVariable) -> 'OptimizationProblem':
        """
        Add a decision variable.

        Args:
            variable: DecisionVariable to add

        Returns:
            self for method chaining
        """
        self._variables.append(variable)
        return self

    def addVariables(self, *variables: DecisionVariable) -> 'OptimizationProblem':
        """
        Add multiple decision variables.

        Args:
            *variables: DecisionVariables to add

        Returns:
            self for method chaining
        """
        for var in variables:
            self._variables.append(var)
        return self

    def setObjective(self, objective: Objective) -> 'OptimizationProblem':
        """
        Set the objective function.

        Args:
            objective: Objective to set

        Returns:
            self for method chaining
        """
        self._objective = objective
        return self

    def addConstraint(self, constraint: Constraint) -> 'OptimizationProblem':
        """
        Add a constraint.

        Args:
            constraint: Constraint to add

        Returns:
            self for method chaining
        """
        self._constraints.append(constraint)
        return self

    def addConstraints(self, *constraints: Constraint) -> 'OptimizationProblem':
        """
        Add multiple constraints.

        Args:
            *constraints: Constraints to add

        Returns:
            self for method chaining
        """
        for constraint in constraints:
            self._constraints.append(constraint)
        return self

    def setFixedVariables(self, pairs: List[tuple]) -> 'OptimizationProblem':
        """
        Fix variables at given values for the duration of the optimization.

        Fixed variables are applied to every model copy before the free
        variables, and their values are visible to objectives/constraints.
        Used by the decomposition workflow to coordinate subproblems.

        Args:
            pairs: List of (DecisionVariable, value) tuples

        Returns:
            self for method chaining
        """
        self._fixed_variables = list(pairs)
        return self

    def getFixedVariables(self) -> List[tuple]:
        """Get the fixed (variable, value) pairs."""
        return list(self._fixed_variables)

    def addScenario(self, model: 'Network', weight: float = 1.0) -> 'OptimizationProblem':
        """
        Add a workload scenario for robust optimization.

        Scenarios are model variants (e.g. different arrival rates) evaluated
        with the same decision variable values. The solver aggregates the
        objective across the base model and all scenarios ('worst' or 'mean',
        see the scenario_aggregation solver option) and enforces constraints
        on every scenario.

        Variables are applied by node/class name, so scenario models must use
        the same names as the base model.

        Args:
            model: LINE Network model variant
            weight: Scenario weight for 'mean' aggregation

        Returns:
            self for method chaining
        """
        self._scenarios.append((model, float(weight)))
        return self

    def getScenarios(self) -> List[tuple]:
        """Get the (model, weight) scenario pairs."""
        return list(self._scenarios)

    def validate(self) -> List[str]:
        """
        Validate the problem specification.

        Returns:
            List of validation error messages (empty if valid)
        """
        errors = []

        if self._model is None:
            errors.append("Model is not set")

        if not self._variables:
            errors.append("No decision variables defined")

        if self._objective is None:
            errors.append("Objective function is not set")

        # Model/variable-kind consistency: LQN models take LQN variables and
        # flat models take flat variables; mixing them silently produces
        # no-ops (a variable whose element is never found in the copy).
        for var in self._variables:
            vtype = var.getVariableType() if hasattr(var, 'getVariableType') \
                else None
            if self._is_layered and vtype in self._FLAT_VAR_TYPES:
                errors.append(
                    f"Variable '{var.getName()}' ({vtype}) is a flat-network "
                    f"variable but the model is a LayeredNetwork")
            elif not self._is_layered and vtype in self._LQN_VAR_TYPES:
                errors.append(
                    f"Variable '{var.getName()}' ({vtype}) is a LayeredNetwork "
                    f"variable but the model is a flat Network")

        return errors

    def isValid(self) -> bool:
        """Check if problem specification is valid."""
        return len(self.validate()) == 0

    def solve(self, options=None, **solver_options) -> OptimizationResult:
        """
        Solve the optimization problem.

        Uses LineOptSolver with differential evolution.

        Args:
            options: An optional LineOptSolverOptions, as MATLAB's
                `problem.solve(options)` takes; its entries are merged UNDER the
                keyword arguments, so an explicit keyword still wins. Passing it
                positionally is what lets an `opt.` script transliterate.
            **solver_options: Options passed to LineOptSolver

        Returns:
            OptimizationResult with solution

        Raises:
            ValueError: If problem is not valid
        """
        errors = self.validate()
        if errors:
            raise ValueError(f"Invalid problem: {', '.join(errors)}")

        if options is not None:
            if hasattr(options, 'toDict'):
                base = dict(options.toDict())
            elif isinstance(options, dict):
                base = dict(options)
            else:
                raise TypeError(
                    "solve() takes a LineOptSolverOptions or a dict of options, got %s"
                    % type(options).__name__)
            base.update(solver_options)
            solver_options = base

        solver = LineOptSolver(self, **solver_options)
        return solver.solve()

    def decompose(self) -> 'DecompositionWorkflow':
        """
        Create a decomposition workflow for this problem.

        Returns:
            DecompositionWorkflow for solving subproblems
        """
        from .decomposition import DecompositionWorkflow
        return DecompositionWorkflow(self)

    def summary(self) -> str:
        """
        Get a summary of the problem.

        Returns:
            String summary of problem specification
        """
        lines = [f"OptimizationProblem: {self._model.getName() if self._model else 'No model'}"]
        lines.append(f"  Variables: {len(self._variables)}")
        for var in self._variables:
            lines.append(f"    - {var.getName()} ({var.getVariableType()})")

        if self._objective:
            obj_type = type(self._objective).__name__
            lines.append(f"  Objective: {obj_type}")
            for constraint in self._objective.getConstraints():
                lines.append(f"    Subject to: {constraint.getName()}")

        lines.append(f"  Constraints: {len(self._constraints)}")
        for constraint in self._constraints:
            lines.append(f"    - {constraint.getName()}")

        return '\n'.join(lines)

    def __repr__(self) -> str:
        model_name = self._model.getName() if self._model else 'None'
        return (f"OptimizationProblem(model='{model_name}', "
                f"vars={len(self._variables)}, "
                f"constraints={len(self._constraints)})")

    # Snake case aliases
    get_model = getModel
    get_variables = getVariables
    get_objective = getObjective
    get_constraints = getConstraints
    add_variable = addVariable
    add_variables = addVariables
    set_objective = setObjective
    add_constraint = addConstraint
    add_constraints = addConstraints
    set_fixed_variables = setFixedVariables
    get_fixed_variables = getFixedVariables
    add_scenario = addScenario
    get_scenarios = getScenarios
    is_valid = isValid
    is_layered = isLayered


# Public API
__all__ = ['OptimizationProblem']
