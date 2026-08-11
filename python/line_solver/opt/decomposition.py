"""
Decomposition Workflow for Multi-Stage Optimization

This module provides classes for decomposing complex optimization problems
into subproblems and solving them iteratively.

Key Classes:
    - SubProblem: A single subproblem in the decomposition
    - DecompositionWorkflow: Manages the decomposition and solution process

Example:
    >>> workflow = DecompositionWorkflow(problem)
    >>> workflow.auto_decompose()
    >>> result = workflow.solve_sequential(max_cycles=5)
"""

import time
import numpy as np
from typing import Dict, List, Any, Optional, Set
from collections import defaultdict
from dataclasses import dataclass, field

try:
    import networkx as nx
    HAS_NETWORKX = True
except ImportError:
    HAS_NETWORKX = False

from .variables import DecisionVariable
from .objectives import Objective, Constraint
from .results import OptimizationResult, SubProblemResult, WorkflowResult
from .solver import LineOptSolver


@dataclass
class SubProblem:
    """
    A subproblem in the decomposition.

    Represents a subset of variables to optimize while keeping
    other variables fixed.

    Attributes:
        name: Subproblem identifier
        variable_type: Type of variables in this subproblem
        variables: List of decision variables
        fixed_values: Values of variables fixed from other subproblems
    """
    name: str
    variable_type: str
    variables: List[DecisionVariable] = field(default_factory=list)
    fixed_values: Dict[str, Any] = field(default_factory=dict)

    def getVariableNames(self) -> List[str]:
        """Get names of variables in this subproblem."""
        return [v.getName() for v in self.variables]

    get_variable_names = getVariableNames


class DecompositionWorkflow:
    """
    Manages decomposition of joint problem into subproblems.

    Supports:
    1. Automatic ordering based on problem types
    2. User-specified dependency DAG
    3. Cyclic iteration for fixed-point convergence

    The default decomposition order is:
    1. server_allocation - Capacity sizing first
    2. station_replicas - Horizontal scaling
    3. service_rate - Speed tuning
    4. job_population - Population sizing
    5. class_priority - Priority assignment
    6. routing - Routing optimization
    7. class_mapping - Service mapping

    Args:
        problem: OptimizationProblem to decompose

    Example:
        >>> workflow = DecompositionWorkflow(problem)
        >>> workflow.auto_decompose()
        >>> result = workflow.solve_sequential(max_cycles=5, tolerance=0.001)
    """

    # Default decomposition order (based on problem structure). Flat-network
    # variable types first, then LayeredNetwork (LQN) variable types; only the
    # types present in a given problem produce subproblems.
    DEFAULT_ORDER = [
        'server_allocation',
        'station_replicas',
        'service_rate',
        'job_population',
        'class_priority',
        'routing',
        'class_mapping',
        # LayeredNetwork (LQN) variable types
        'processor_multiplicity',
        'task_multiplicity',
        'task_replication',
        'host_demand',
        'think_time',
    ]

    def __init__(self, problem: 'OptimizationProblem'):
        """
        Initialize decomposition workflow.

        Args:
            problem: OptimizationProblem to decompose
        """
        self._problem = problem
        self._subproblems: List[SubProblem] = []
        self._dependency_graph: Dict[str, Set[str]] = defaultdict(set)
        self._solver_options: dict = {}

    def getProblem(self) -> 'OptimizationProblem':
        """Get the optimization problem."""
        return self._problem

    def getSubProblems(self) -> List[SubProblem]:
        """Get the list of subproblems."""
        return self._subproblems

    def setSolverOptions(self, **options) -> 'DecompositionWorkflow':
        """
        Set options for subproblem solvers.

        Args:
            **options: Options passed to LineOptSolver

        Returns:
            self for chaining
        """
        self._solver_options = options
        return self

    def autoDecompose(self) -> 'DecompositionWorkflow':
        """
        Automatically decompose based on variable types.

        Groups variables by their type and creates subproblems
        in the default order.

        Returns:
            self for chaining
        """
        # Group variables by type
        var_by_type: Dict[str, List[DecisionVariable]] = defaultdict(list)
        for var in self._problem.getVariables():
            var_type = var.getVariableType()
            var_by_type[var_type].append(var)

        # Create subproblems in default order
        self._subproblems = []
        for var_type in self.DEFAULT_ORDER:
            if var_type in var_by_type:
                subproblem = SubProblem(
                    name=var_type,
                    variable_type=var_type,
                    variables=var_by_type[var_type]
                )
                self._subproblems.append(subproblem)

        # Add any remaining types not in default order
        for var_type, variables in var_by_type.items():
            if var_type not in self.DEFAULT_ORDER:
                subproblem = SubProblem(
                    name=var_type,
                    variable_type=var_type,
                    variables=variables
                )
                self._subproblems.append(subproblem)

        return self

    def setDependency(self, from_problem: str, to_problem: str) -> 'DecompositionWorkflow':
        """
        Add explicit dependency between subproblems.

        Args:
            from_problem: Name of prerequisite subproblem
            to_problem: Name of dependent subproblem

        Returns:
            self for chaining
        """
        self._dependency_graph[to_problem].add(from_problem)
        return self

    def addSubProblem(self, name: str, variables: List[DecisionVariable],
                      after: List[str] = None) -> 'DecompositionWorkflow':
        """
        Add a custom subproblem.

        Args:
            name: Subproblem name
            variables: Variables to include
            after: List of subproblems this depends on

        Returns:
            self for chaining
        """
        var_type = variables[0].getVariableType() if variables else 'custom'
        subproblem = SubProblem(
            name=name,
            variable_type=var_type,
            variables=variables
        )
        self._subproblems.append(subproblem)

        if after:
            for dep in after:
                self.setDependency(dep, name)

        return self

    def _getExecutionOrder(self) -> List[SubProblem]:
        """
        Get subproblems in execution order respecting dependencies.

        Uses topological sort if networkx is available, otherwise
        uses the order they were added.

        Returns:
            List of subproblems in execution order
        """
        if not self._dependency_graph or not HAS_NETWORKX:
            return list(self._subproblems)

        # Build DAG
        G = nx.DiGraph()
        for sp in self._subproblems:
            G.add_node(sp.name)
        for to_node, from_nodes in self._dependency_graph.items():
            for from_node in from_nodes:
                G.add_edge(from_node, to_node)

        # Topological sort
        try:
            order = list(nx.topological_sort(G))
            name_to_sp = {sp.name: sp for sp in self._subproblems}
            return [name_to_sp[name] for name in order if name in name_to_sp]
        except nx.NetworkXUnfeasible:
            # Cycle detected, use default order
            return list(self._subproblems)

    def solveSequential(self, max_cycles: int = 10,
                        tolerance: float = 0.01,
                        verbose: bool = False) -> WorkflowResult:
        """
        Solve subproblems sequentially with optional cycling.

        Iterates through subproblems, fixing values from previous
        solutions. Repeats cycles until convergence or max_cycles.

        Args:
            max_cycles: Maximum number of complete cycles
            tolerance: Convergence tolerance for objective change
            verbose: Print progress information

        Returns:
            WorkflowResult with solution and statistics
        """
        start_time = time.time()
        result = WorkflowResult()

        if not self._subproblems:
            result.converged = True
            result.total_solve_time = 0.0
            return result

        # Get execution order
        ordered_subproblems = self._getExecutionOrder()

        # Track fixed variable values across subproblems
        fixed_values: Dict[str, Any] = {}

        # Track objective history for convergence
        prev_objective = float('inf')

        for cycle in range(max_cycles):
            cycle_start = time.time()

            if verbose:
                print(f"\n=== Cycle {cycle + 1}/{max_cycles} ===")

            # Solve each subproblem in order
            for subproblem in ordered_subproblems:
                if verbose:
                    print(f"  Solving subproblem: {subproblem.name}")

                # Create partial problem with only these variables
                partial_problem = self._createPartialProblem(subproblem, fixed_values)

                # Solve
                solver_opts = {**self._solver_options}
                if verbose:
                    solver_opts['verbose'] = True

                solver = LineOptSolver(partial_problem, **solver_opts)
                sp_result = solver.solve()

                # Store result
                sp_result_obj = SubProblemResult(
                    name=subproblem.name,
                    result=sp_result,
                    variables_fixed=dict(fixed_values)
                )
                result.subproblem_results[subproblem.name] = sp_result_obj

                # Update fixed values with solution
                for var_name, value in sp_result.variable_values.items():
                    fixed_values[var_name] = value

                if verbose:
                    print(f"    Objective: {sp_result.objective_value:.4f}")

            # Evaluate full objective after cycle
            current_objective = self._evaluateFullObjective(fixed_values)
            result.objective_history.append(current_objective)

            if verbose:
                print(f"  Cycle {cycle + 1} objective: {current_objective:.4f}")

            # Check convergence
            if abs(current_objective - prev_objective) < tolerance:
                result.converged = True
                if verbose:
                    print(f"  Converged after {cycle + 1} cycles")
                break

            prev_objective = current_objective
            result.cycles_completed = cycle + 1

        # Finalize result
        result.final_objective = result.objective_history[-1] if result.objective_history else float('inf')
        result.final_variable_values = dict(fixed_values)
        result.total_solve_time = time.time() - start_time

        return result

    def solveHierarchical(self, verbose: bool = False) -> WorkflowResult:
        """
        Solve following DAG dependencies (single pass).

        Uses topological sort to determine order and solves
        each subproblem once.

        Args:
            verbose: Print progress information

        Returns:
            WorkflowResult
        """
        return self.solveSequential(max_cycles=1, tolerance=0.0, verbose=verbose)

    def _createPartialProblem(self, subproblem: SubProblem,
                              fixed_values: Dict[str, Any]) -> 'OptimizationProblem':
        """
        Create a partial problem for a subproblem.

        Variables outside the subproblem that already have values from
        previously solved blocks are fixed at those values, so this block
        optimizes against the coordinated model state (Gauss-Seidel style)
        rather than the pristine base model.

        Args:
            subproblem: SubProblem to create problem for
            fixed_values: Fixed variable values from prior subproblems

        Returns:
            OptimizationProblem with only subproblem variables free
        """
        from .problem import OptimizationProblem

        # Create new problem with same model
        partial = OptimizationProblem(self._problem.getModel())

        # Add only this subproblem's variables
        for var in subproblem.variables:
            partial.addVariable(var)

        # Fix all other variables at their current values (if solved already)
        sub_names = set(subproblem.getVariableNames())
        fixed_pairs = []
        for var in self._problem.getVariables():
            name = var.getName()
            if name not in sub_names and name in fixed_values:
                fixed_pairs.append((var, fixed_values[name]))
        partial.setFixedVariables(fixed_pairs)

        # Use same objective and constraints
        partial.setObjective(self._problem.getObjective())
        for constraint in self._problem.getConstraints():
            partial.addConstraint(constraint)

        # Propagate workload scenarios
        for scenario_model, weight in self._problem.getScenarios():
            partial.addScenario(scenario_model, weight)

        # Store fixed values for evaluation
        subproblem.fixed_values = dict(fixed_values)

        return partial

    # ---- LayeredNetwork (LQN) layer-wise decomposition ------------------

    def solveLayered(self, max_cycles: int = 10, tolerance: float = 0.01,
                     auto_freeze: bool = True, freeze_tol: float = 1e-3,
                     frozen_layers: Optional[List[str]] = None,
                     verbose: bool = False) -> WorkflowResult:
        """Solve an LQN by layer, optionally freezing converged layers.

        Groups decision variables by the LQN layer they perturb (host or task
        layer) and cycles Gauss-Seidel over the layer groups, fixing every
        other layer's variables at their current values while one layer is
        optimized. This is the LQN analogue of :meth:`solveSequential`, but the
        subproblems are LAYERS rather than variable types, which is what makes
        layer freezing meaningful.

        Freezing has two, composable, sources:

        * ``frozen_layers`` -- an explicit seed set held fixed throughout.
        * ``auto_freeze`` -- adaptive: after each cycle a layer whose
          representative node metrics moved less than ``freeze_tol`` (relative)
          is frozen and skipped; it is unfrozen again if any still-active layer
          later moves by more than ``freeze_tol`` (its coupling changed).

        Convergence is on the full penalized objective delta between cycles, or
        when every layer is frozen. Falls back to :meth:`solveSequential` for a
        flat network. The returned :class:`WorkflowResult` carries the extra
        attributes ``frozen_layers`` (final frozen set) and ``model_evaluations``
        (total LINE solves), so a caller can confirm freezing cut solve count.
        """
        from .problem import OptimizationProblem
        from .evaluator import LineEvaluator
        from .layered import is_layered

        start_time = time.time()
        model = self._problem.getModel()
        if not is_layered(model):
            return self.solveSequential(max_cycles, tolerance, verbose)

        # Group variables by their primary (first) layer.
        groups: Dict[str, List[DecisionVariable]] = {}
        for var in self._problem.getVariables():
            layers = var.getLayer(model) or ['_nolayer']
            groups.setdefault(layers[0], []).append(var)

        objective = self._problem.getObjective()
        penalty_weight = self._solver_options.get(
            'penalty_weight', LineOptSolver.defaultOptions()['penalty_weight'])
        evaluator = LineEvaluator(model, self._problem.getVariables())

        frozen: Set[str] = set(frozen_layers or [])
        fixed_values: Dict[str, Any] = {}
        prev_sig: Dict[str, tuple] = {}
        prev_objective = float('inf')
        model_evaluations = 0

        result = WorkflowResult()
        for cycle in range(max_cycles):
            if verbose:
                print(f"\n=== LQN cycle {cycle + 1}/{max_cycles} "
                      f"(frozen: {sorted(frozen)}) ===")

            for layer, layer_vars in groups.items():
                if layer in frozen:
                    continue
                partial = OptimizationProblem(model)
                for var in layer_vars:
                    partial.addVariable(var)
                # Fix every other layer's variables at their current values.
                sub_names = {v.getName() for v in layer_vars}
                fixed_pairs = []
                for var in self._problem.getVariables():
                    name = var.getName()
                    if name not in sub_names and name in fixed_values:
                        fixed_pairs.append((var, fixed_values[name]))
                partial.setFixedVariables(fixed_pairs)
                partial.setObjective(objective)
                for constraint in self._problem.getConstraints():
                    partial.addConstraint(constraint)

                solver_opts = {**self._solver_options}
                if verbose:
                    solver_opts['verbose'] = True
                solver = LineOptSolver(partial, **solver_opts)
                sp_result = solver.solve()
                model_evaluations += sp_result.model_evaluations

                result.subproblem_results[layer] = SubProblemResult(
                    name=layer, result=sp_result,
                    variables_fixed=dict(fixed_values))
                for name, value in sp_result.variable_values.items():
                    fixed_values[name] = value

            # Full evaluation: objective + per-layer signatures for freezing.
            eval_result = evaluator.evaluateValues(fixed_values)
            model_evaluations += 1
            if not eval_result.feasible:
                current_objective = float('inf')
                sig = {}
            else:
                current_objective = objective.evaluateWithPenalty(
                    eval_result, fixed_values, penalty_weight)
                for constraint in self._problem.getConstraints():
                    current_objective += constraint.evaluate(
                        eval_result, fixed_values) * penalty_weight
                sig = self._layerSignatures(eval_result, groups.keys())

            if auto_freeze and prev_sig:
                moved = {L: self._sigDelta(prev_sig.get(L), sig.get(L))
                         for L in groups}
                active_moved = any(
                    moved[a] > freeze_tol for a in groups if a not in frozen)
                for L in list(groups):
                    if L in frozen:
                        # Unfreeze if a still-active layer has moved.
                        if active_moved:
                            frozen.discard(L)
                    elif moved[L] < freeze_tol:
                        frozen.add(L)

            prev_sig = sig
            result.objective_history.append(current_objective)
            if verbose:
                print(f"  cycle objective: {current_objective:.4f}")

            if abs(current_objective - prev_objective) < tolerance:
                result.converged = True
                if verbose:
                    print(f"  converged after {cycle + 1} cycles")
                break
            prev_objective = current_objective
            result.cycles_completed = cycle + 1
            # All layers frozen: nothing left to optimize.
            if len(frozen) >= len(groups):
                result.converged = True
                break

        result.final_objective = (result.objective_history[-1]
                                  if result.objective_history else float('inf'))
        result.final_variable_values = dict(fixed_values)
        result.total_solve_time = time.time() - start_time
        # Extra LQN diagnostics (dynamic attributes on the dataclass).
        result.frozen_layers = sorted(frozen)
        result.model_evaluations = model_evaluations
        return result

    @staticmethod
    def _layerSignatures(eval_result, layers) -> Dict[str, tuple]:
        """Representative (Util, QLen, Tput, RespT) per layer, keyed by its node.

        A layer named after a processor or task has a same-named node in the LQN
        average table; its metrics are the layer's convergence signature. Layers
        without a matching node (e.g. the '_nolayer' bucket) get a None signature
        so they never auto-freeze.
        """
        sig = {}
        for layer in layers:
            util = eval_result.getUtilization(layer)
            qlen = eval_result.getQueueLength(layer)
            tput = eval_result.getThroughput(layer)
            respt = eval_result.getResponseTime(layer)
            if any(v not in (0.0, float('inf')) for v in (util, qlen, tput)):
                sig[layer] = (util, qlen, tput, respt)
            else:
                sig[layer] = None
        return sig

    @staticmethod
    def _sigDelta(a: Optional[tuple], b: Optional[tuple]) -> float:
        """Max relative change between two layer signatures (inf if unknown)."""
        if a is None or b is None:
            return float('inf')
        eps = 1e-12
        delta = 0.0
        for ai, bi in zip(a, b):
            if not (np.isfinite(ai) and np.isfinite(bi)):
                continue
            delta = max(delta, abs(bi - ai) / (abs(ai) + eps))
        return delta

    def _evaluateFullObjective(self, variable_values: Dict[str, Any]) -> float:
        """
        Evaluate the full penalized objective with all variable values.

        Applies every variable to a model copy, solves it with LINE, extracts
        real metrics, and returns the objective including constraint
        penalties, so cycle convergence is measured on the same quantity the
        subproblem solvers optimize.

        Args:
            variable_values: All variable values

        Returns:
            Penalized objective value (inf if the model evaluation fails)
        """
        from .evaluator import LineEvaluator
        from .solver import LineOptSolver

        penalty_weight = self._solver_options.get(
            'penalty_weight', LineOptSolver.defaultOptions()['penalty_weight'])

        # Evaluate with all variables applied at their current values
        evaluator = LineEvaluator(self._problem.getModel(),
                                  self._problem.getVariables())
        eval_result = evaluator.evaluateValues(variable_values)

        if not eval_result.feasible:
            return float('inf')

        objective = self._problem.getObjective()
        value = objective.evaluateWithPenalty(
            eval_result, variable_values, penalty_weight)
        for constraint in self._problem.getConstraints():
            violation = constraint.evaluate(eval_result, variable_values)
            value += violation * penalty_weight

        return value

    # Snake case aliases
    get_problem = getProblem
    get_subproblems = getSubProblems
    set_solver_options = setSolverOptions
    auto_decompose = autoDecompose
    set_dependency = setDependency
    add_subproblem = addSubProblem
    solve_sequential = solveSequential
    solve_hierarchical = solveHierarchical
    solve_layered = solveLayered


# Public API
__all__ = [
    'SubProblem',
    'DecompositionWorkflow',
]
