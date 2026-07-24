# line_solver.opt

**Optimization Framework for LINE Queueing Network Models**

line_solver.opt provides a declarative interface for specifying and solving optimization problems on LINE queueing network models. It uses scipy's differential evolution algorithm with LINE's SolverAuto as the model evaluator. Both flat `Network` models and `LayeredNetwork` (LQN) models are supported; LQN models are solved with `SolverLN` and tuned with LQN-specific decision variables (see "LayeredNetwork (LQN) optimization" below).

## Features

- **7 Decision Variable Types**: Server allocation, station replicas, routing probabilities, class-service mapping, service rates, job population, class priorities
- **Pre-defined Objectives**: Cost minimization with SLA constraints, performance maximization with budget constraints
- **Constraint Types**: Per-station response time, end-to-end (system) response time, throughput, utilization, and budget constraints
- **Differential Evolution Solver**: Global optimization with automatic encoding of discrete decisions, configuration-level caching, and per-evaluation time limits
- **Bisection Solver**: Exact O(log n) sizing for single integer-variable problems with monotone feasibility
- **Pareto Sweeps**: Epsilon-constraint sweeps for cost-performance tradeoff frontiers
- **Scenario Robustness**: Optimize against multiple workload scenarios (worst-case or weighted mean)
- **Decomposition Workflow**: Break complex problems into subproblems coordinated by fixed-value propagation (Gauss-Seidel cycles)

## Installation

```bash
# Install from source
pip install -e .

# Dependencies
pip install line-solver scipy numpy
```

## Quick Example

```python
from line_solver import Network, Queue, Source, Sink, OpenClass, Exp, SchedStrategy
from line_solver import (
    OptimizationProblem, ServerAllocation,
    MinimizeCost, ResponseTimeConstraint
)

# Create LINE model
model = Network("WebServer")
source = Source(model, "Arrivals")
queue = Queue(model, "Server", SchedStrategy.FCFS)
sink = Sink(model, "Departures")

jobclass = OpenClass(model, "Requests")
source.setArrival(jobclass, Exp(10.0))  # 10 req/sec arrival
queue.setService(jobclass, Exp(2.0))     # 2 sec service time

model.addLink(source, queue)
model.addLink(queue, sink)

# Define optimization problem
problem = OptimizationProblem(model)
problem.add_variable(ServerAllocation(queue, bounds=(1, 20)))
problem.set_objective(MinimizeCost(
    server_cost={queue: 100.0}
))
problem.add_constraint(
    ResponseTimeConstraint(queue, jobclass, max_value=1.0)
)

# Solve
result = problem.solve(verbose=True)
print(f"Optimal servers: {result.variable_values['Server_servers']}")
print(f"Cost: ${result.objective_value:.2f}")
```

## Decision Variables

| Variable Type | Description |
|---------------|-------------|
| `ServerAllocation` | Number of servers at a station |
| `StationReplicas` | Number of identical station copies |
| `RoutingProbabilities` | Routing probabilities for a job class |
| `ClassServiceMapping` | Class-to-station service mapping |
| `ServiceRate` | Processing rate at a station |
| `JobPopulation` | Number of jobs in a closed class |
| `ClassPriority` | Priority ordering of job classes |
| `HostDemand` | Mean host demand of an LQN activity (continuous) |
| `ActivityThinkTime` / `TaskThinkTime` | Think time of an LQN activity / task (continuous) |
| `TaskMultiplicity` | Thread/instance count of an LQN task (integer) |
| `TaskReplication` | Fan-out replicas of an LQN task (integer) |
| `ProcessorMultiplicity` | Core count of an LQN processor (integer) |

## Objectives

| Objective | Description |
|-----------|-------------|
| `MinimizeCost` | Minimize infrastructure cost subject to SLA constraints |
| `MaximizePerformance` | Maximize throughput/response time subject to budget |
| `MinimizeSystemResponseTime` | Minimize mean end-to-end response time (load balancing, routing) |

## Constraints

| Constraint | Description |
|------------|-------------|
| `ResponseTimeConstraint` | Station response time <= max_value |
| `SystemResponseTimeConstraint` | End-to-end response time <= max_value |
| `ThroughputConstraint` | Throughput >= min_value |
| `UtilizationConstraint` | Utilization <= max_value |
| `BudgetConstraint` | Total cost <= budget |

## Exact Sizing with Bisection

For a single integer variable with monotone feasibility (the standard
server sizing pattern), `BisectionSolver` finds the exact optimum in
O(log n) LINE solves:

```python
from line_solver import BisectionSolver

result = BisectionSolver(problem).solve()            # smallest feasible value
result = BisectionSolver(problem, direction='max_feasible').solve()
```

## Cost-Performance Tradeoffs

`ParetoSweep` sweeps an SLA bound (epsilon-constraint method) and returns
the non-dominated cost frontier:

```python
from line_solver import ParetoSweep

sweep = ParetoSweep(
    problem,
    constraint_factory=lambda eps: UtilizationConstraint(queue, max_value=eps),
    epsilons=[0.3, 0.5, 0.75],
    solver='bisection',   # or 'de'
)
points = sweep.solve(seed=42)
frontier = sweep.getFrontier()
```

## Scenario-Based Robust Optimization

Add workload scenarios (model variants with the same node/class names) to
optimize against uncertainty. Constraints are enforced on all scenarios;
the objective is aggregated by 'worst' (default) or 'mean':

```python
problem.add_scenario(high_load_model, weight=1.0)
result = problem.solve(scenario_aggregation='worst')
```

## Decomposition for Complex Problems

For problems with many decision variables, use decomposition:

```python
from line_solver import DecompositionWorkflow

# Create workflow
workflow = DecompositionWorkflow(problem)

# Auto-decompose by variable type
workflow.auto_decompose()

# Solve with cycling for fixed-point convergence
result = workflow.solve_sequential(max_cycles=5, tolerance=0.001)

print(f"Converged: {result.converged}")
print(f"Cycles: {result.cycles_completed}")
```

## LayeredNetwork (LQN) optimization

Pass a `LayeredNetwork` to `OptimizationProblem` and it is solved with
`SolverLN`; objectives and constraints then resolve by LQN **node name**
(processor, task, entry, or activity). Use the LQN decision variables above.

```python
from line_solver import (LayeredNetwork, Processor, Task, Entry, Activity,
                         Exp, SchedStrategy)
from line_solver import (OptimizationProblem, HostDemand,
                         MinimizeSystemResponseTime, UtilizationConstraint)

model = LayeredNetwork('app')
P1 = Processor(model, 'P1', 2, SchedStrategy.PS)
T1 = Task(model, 'T1', 50, SchedStrategy.REF).on(P1).set_think_time(Exp(0.5))
E1 = Entry(model, 'E1').on(T1)
Activity(model, 'AS1', Exp(10)).on(T1).bound_to(E1)

problem = OptimizationProblem(model)
problem.add_variable(HostDemand('AS1', bounds=(0.02, 0.2)))
problem.set_objective(MinimizeSystemResponseTime(
    'T1', subject_to=[UtilizationConstraint('P1', max_value=0.9)]))

# gradient path with a chosen LQN gradient source (fd is the robust default)
result = problem.solve(optimizer='gradient', lqn_gradient='fd')
```

### LQN gradient sources (`lqn_gradient`)

`SolverLN.getSensitivityTable` gives WITHIN-LAYER partial derivatives (it holds
the fixed-point layer parameters constant and omits the cross-layer coupling
term), so it is a biased estimate of the total derivative. The `lqn_gradient`
option selects how gradients are formed:

| Value | Meaning |
|-------|---------|
| `fd` (default) | Finite-difference the whole `LayeredNetwork` per parameter -- correct total derivative, robust, `N` extra ensemble solves |
| `partial_sens` | Assemble the direction from the per-layer partial derivatives -- cheap, but biased |
| `partial_plus_fd` | `partial_sens` direction, corrected by a full-model finite difference every `fd_refresh` gradient evaluations |

### Layer freezing

Hold a subset of layers fixed while optimizing the rest, either explicitly or
adaptively:

```python
# Explicit: freeze host/task layers by name (their variables are held fixed)
result = problem.solve(optimizer='gradient', frozen_layers=['P2'])

# Adaptive: freeze layers whose metrics have converged, unfreeze on coupling
wf = problem.decompose()
r = wf.solveLayered(max_cycles=6, auto_freeze=True, freeze_tol=1e-2)
print(r.frozen_layers, r.model_evaluations)
```

## Solver Options

```python
from line_solver import LineOptSolverOptions

options = LineOptSolverOptions()
options.setMaxIterations(200)
options.setPopsize(20)
options.setVerbose(True)
options.setSeed(42)

result = problem.solve(**options.toDict())
```

### Available Options

| Option | Default | Description |
|--------|---------|-------------|
| `strategy` | 'best1bin' | DE mutation strategy |
| `popsize` | 15 | Population size multiplier |
| `mutation` | (0.5, 1.0) | Mutation constant |
| `recombination` | 0.7 | Crossover probability |
| `max_iterations` | 100 | Maximum generations |
| `time_limit` | 300.0 | Maximum solve time (seconds, enforced per evaluation) |
| `verbose` | False | Print progress |
| `seed` | None | Random seed |
| `scenario_aggregation` | 'worst' | Objective aggregation across scenarios |
| `optimizer` | 'evolution' | Backend: 'evolution', 'gradient', or 'auto' |
| `lqn_gradient` | 'fd' | LQN gradient source: 'fd', 'partial_sens', 'partial_plus_fd' |
| `fd_refresh` | 5 | partial_plus_fd full-FD correction period |
| `frozen_layers` | None | LQN layer names held fixed during optimization |

## Examples

Runnable examples are in the `examples/opt/` folder, each a
self-contained script that prints the optimum it finds alongside the
closed-form value where one exists:

| Example | Topic |
|---------|-------|
| `opt_server_sizing.py` | Minimum servers under a utilization SLA (differential evolution) |
| `opt_bisection_sizing.py` | Same sizing solved exactly with `BisectionSolver` |
| `opt_service_rate.py` | Cheapest continuous service rate meeting a response time SLA |
| `opt_load_balancing.py` | Routing split minimizing end-to-end response time |
| `opt_population_sizing.py` | Largest population sustaining an interactive SLA |
| `opt_robust_sizing.py` | Worst-case sizing across workload scenarios |
| `opt_pareto_frontier.py` | Cost vs utilization-bound frontier via `ParetoSweep` |
| `opt_decomposition.py` | Two-variable problem via `DecompositionWorkflow` |
| `opt_lqn_hostdemand.py` | LQN host-demand tuning: gradient modes + layer freezing |

## Documentation

Full documentation is available in the `doc/` directory:

- **User Manual (PDF)**: `doc/manual/line_solver.opt-manual.pdf` (rebuild with `doc/manual/build.sh`)
- **LaTeX Manual**: `doc/latex/manual.tex`
- **Sphinx Documentation**: `doc/source/`

Build the Sphinx documentation:

```bash
cd doc
make html
```

## Requirements

- Python >= 3.8
- LINE Solver (line-solver)
- scipy >= 1.7.0
- numpy >= 1.20.0
- networkx >= 2.6 (for decomposition DAG)

## License

BSD-3-Clause License

## Citation

If you use line_solver.opt in your research, please cite:

```bibtex
@software{line_solver.opt,
  title = {line_solver.opt: Optimization Framework for LINE Queueing Network Models},
  author = {Casale, Giuliano},
  year = {2024},
  url = {https://github.com/line-solver/line_solver.opt}
}
```
