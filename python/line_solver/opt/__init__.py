"""
line_solver.opt: optimization of LINE queueing network models.

This subpackage provides a declarative interface for specifying and solving
optimization problems on LINE queueing network models. Its public classes are
also re-exported from the top-level ``line_solver`` package.

Key Features:
    - Decision variables for servers, rates, routing, and more
    - Pre-defined objectives: MinimizeCost, MaximizePerformance
    - Constraint handling with penalty-based optimization
    - Differential evolution solver using scipy
    - Decomposition workflow for complex multi-stage problems

Quick Start:
    >>> from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
    >>> from line_solver import OptimizationProblem, ServerAllocation, MinimizeCost, ResponseTimeConstraint
    >>>
    >>> # Create LINE model
    >>> model = Network("WebServer")
    >>> source = Source(model, "Arrivals")
    >>> queue = Queue(model, "Server", SchedStrategy.FCFS)
    >>> sink = Sink(model, "Departures")
    >>>
    >>> jobclass = OpenClass(model, "Requests")
    >>> source.setArrival(jobclass, Exp(10.0))
    >>> queue.setService(jobclass, Exp(2.0))
    >>>
    >>> model.addLink(source, queue)
    >>> model.addLink(queue, sink)
    >>>
    >>> # Define optimization problem
    >>> problem = OptimizationProblem(model)
    >>> problem.add_variable(ServerAllocation(queue, bounds=(1, 20)))
    >>> problem.set_objective(MinimizeCost(
    ...     server_cost={queue: 100.0},
    ...     subject_to=[ResponseTimeConstraint(queue, max_value=1.0)]
    ... ))
    >>>
    >>> # Solve
    >>> result = problem.solve(verbose=True)
    >>> print(f"Optimal servers: {result.variable_values['Server_servers']}")

Modules:
    - problem: OptimizationProblem class
    - variables: Decision variable types
    - objectives: Objective functions and constraints
    - solver: Differential evolution solver
    - decomposition: Workflow decomposition
    - evaluator: LINE model evaluation wrapper
    - results: Result containers
"""

__version__ = "0.1.0"
__author__ = "LINE Project"

# Problem specification
from .problem import OptimizationProblem

# Decision variables
from .variables import (
    DecisionVariable,
    ServerAllocation,
    StationReplicas,
    RoutingProbabilities,
    ClassServiceMapping,
    ServiceRate,
    JobPopulation,
    ClassPriority,
    # LayeredNetwork (LQN) variables
    HostDemand,
    ActivityThinkTime,
    TaskThinkTime,
    TaskMultiplicity,
    TaskReplication,
    ProcessorMultiplicity,
)

# LayeredNetwork (LQN) support helpers
from .layered import is_layered

# Objectives and constraints
from .objectives import (
    Objective,
    MinimizeCost,
    MaximizePerformance,
    MinimizeSystemResponseTime,
    Constraint,
    ResponseTimeConstraint,
    SystemResponseTimeConstraint,
    ThroughputConstraint,
    UtilizationConstraint,
    BudgetConstraint,
)

# Solvers
from .solver import LineOptSolver, LineOptSolverOptions
from .sizing import BisectionSolver
from .pareto import ParetoPoint, ParetoSweep

# Decomposition
from .decomposition import SubProblem, DecompositionWorkflow

# Results
from .results import (
    EvaluationResult,
    OptimizationResult,
    SubProblemResult,
    WorkflowResult,
)

# Evaluator (less commonly used directly)
from .evaluator import LineEvaluator

# Public API
__all__ = [
    # Version
    '__version__',

    # Problem
    'OptimizationProblem',

    # Variables
    'DecisionVariable',
    'ServerAllocation',
    'StationReplicas',
    'RoutingProbabilities',
    'ClassServiceMapping',
    'ServiceRate',
    'JobPopulation',
    'ClassPriority',

    # LayeredNetwork (LQN) variables
    'HostDemand',
    'ActivityThinkTime',
    'TaskThinkTime',
    'TaskMultiplicity',
    'TaskReplication',
    'ProcessorMultiplicity',
    'is_layered',

    # Objectives
    'Objective',
    'MinimizeCost',
    'MaximizePerformance',
    'MinimizeSystemResponseTime',

    # Constraints
    'Constraint',
    'ResponseTimeConstraint',
    'SystemResponseTimeConstraint',
    'ThroughputConstraint',
    'UtilizationConstraint',
    'BudgetConstraint',

    # Solvers
    'LineOptSolver',
    'LineOptSolverOptions',
    'BisectionSolver',
    'ParetoPoint',
    'ParetoSweep',

    # Decomposition
    'SubProblem',
    'DecompositionWorkflow',

    # Results
    'EvaluationResult',
    'OptimizationResult',
    'SubProblemResult',
    'WorkflowResult',

    # Evaluator
    'LineEvaluator',
]
