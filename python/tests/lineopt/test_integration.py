"""
End-to-end integration tests against the real LINE solver.

Unlike the rest of the suite (which uses mock LINE objects from conftest.py),
these tests build an actual line_solver Network, run SolverAUTO through the
LineEvaluator, and check that performance metrics are extracted and drive the
optimizer correctly. They are skipped automatically if line_solver is not
installed.

Model: an M/M/c open queue with arrival rate lambda = 3.0 and per-server
service rate mu = 1.0. For a stable open network the throughput equals the
arrival rate (3.0) and the station utilization is lambda / (c * mu) = 3 / c,
which is monotonically decreasing in the server count c. These closed-form
values make the assertions deterministic.
"""

import numpy as np
import pytest

# Skip the whole module if LINE is unavailable.
line_solver = pytest.importorskip("line_solver")

from line_solver import (
    Network, Queue, Delay, Source, Sink,
    OpenClass, ClosedClass, Exp, SchedStrategy,
)
from line_solver.opt import (
    OptimizationProblem,
    ServerAllocation,
    ClassServiceMapping,
    JobPopulation,
    RoutingProbabilities,
    ServiceRate,
    StationReplicas,
    MinimizeCost,
    ResponseTimeConstraint,
    SystemResponseTimeConstraint,
    UtilizationConstraint,
    BisectionSolver,
    ParetoSweep,
    DecompositionWorkflow,
)
from line_solver.opt.evaluator import LineEvaluator

ARRIVAL_RATE = 3.0
SERVICE_RATE = 1.0


def build_model():
    """Build an M/M/c open queue and return (model, queue, jobclass)."""
    model = Network("MMc")
    source = Source(model, "Arrivals")
    queue = Queue(model, "Server", SchedStrategy.FCFS)
    sink = Sink(model, "Departures")

    jobs = OpenClass(model, "Jobs")
    source.setArrival(jobs, Exp(ARRIVAL_RATE))
    queue.setService(jobs, Exp(SERVICE_RATE))

    model.addLink(source, queue)
    model.addLink(queue, sink)
    return model, queue, jobs


def _encode_servers(c, lo, hi):
    """Encode an integer server count c into the [0, 1] DE coordinate."""
    return np.array([(c - lo) / (hi - lo)])


class TestEvaluatorAgainstRealLine:
    """The LineEvaluator must populate metrics from a real SolverAUTO run."""

    def test_evaluate_populates_metrics(self):
        model, queue, _ = build_model()
        evaluator = LineEvaluator(model, [ServerAllocation(queue, bounds=(1, 10))])

        result = evaluator.evaluate(_encode_servers(6, 1, 10))

        assert result.feasible is True
        assert result.solver_used  # a concrete solver name was recorded
        # Regression guard for the getNodeIndex/getJobClassIndex bug: these
        # dicts were silently left empty before the fix.
        assert result.response_times, "response_times must not be empty"
        assert result.throughputs, "throughputs must not be empty"
        assert result.utilizations, "utilizations must not be empty"

        key = ("Server", "Jobs")
        assert np.isfinite(result.getResponseTime(*key))
        assert result.getResponseTime(*key) > 0

    def test_throughput_equals_arrival_rate(self):
        model, queue, _ = build_model()
        evaluator = LineEvaluator(model, [ServerAllocation(queue, bounds=(1, 10))])

        result = evaluator.evaluate(_encode_servers(6, 1, 10))

        # Stable open network: departure throughput == arrival rate.
        assert np.isclose(result.getThroughput("Server", "Jobs"), ARRIVAL_RATE,
                          atol=1e-3)

    @pytest.mark.parametrize("c", [4, 5, 6, 10])
    def test_utilization_matches_closed_form(self, c):
        model, queue, _ = build_model()
        evaluator = LineEvaluator(model, [ServerAllocation(queue, bounds=(1, 10))])

        result = evaluator.evaluate(_encode_servers(c, 1, 10))

        expected = ARRIVAL_RATE / (c * SERVICE_RATE)  # = 3 / c
        assert np.isclose(result.getUtilization("Server"), expected, atol=1e-3)


class TestSolveAgainstRealLine:
    """A full optimization run must be driven by the extracted metrics."""

    def test_minimize_cost_under_utilization_constraint(self):
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
        # util = 3 / c <= 0.5  <=>  c >= 6, so the cheapest feasible c is 6.
        problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))

        result = problem.solve(max_iterations=40, popsize=15, seed=42,
                               verbose=False)

        assert result.feasible is True
        assert result.variable_values["Server_servers"] == 6
        assert np.isclose(result.objective_value, 60.0)

    def test_infeasible_metrics_would_be_caught(self):
        """
        Sanity check that the constraint actually reads live utilization:
        an impossible bound (util <= 0.1 needs c >= 30, above the max of 10)
        forces the optimizer to the upper server bound.
        """
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
        problem.add_constraint(UtilizationConstraint(queue, max_value=0.1))

        result = problem.solve(max_iterations=40, popsize=15, seed=42,
                               verbose=False)

        # Best achievable within bounds is the maximum server count.
        assert result.variable_values["Server_servers"] == 10


def build_parallel_model():
    """Source -> {Q1, Q2} parallel queues -> Sink, single open class."""
    model = Network("Parallel")
    source = Source(model, "S")
    q1 = Queue(model, "Q1", SchedStrategy.FCFS)
    q2 = Queue(model, "Q2", SchedStrategy.FCFS)
    sink = Sink(model, "K")
    jobs = OpenClass(model, "Jobs")
    source.setArrival(jobs, Exp(2.0))
    q1.setService(jobs, Exp(3.0))
    q2.setService(jobs, Exp(3.0))
    model.addLink(source, q1)
    model.addLink(source, q2)
    model.addLink(q1, sink)
    model.addLink(q2, sink)
    return model, [q1, q2], jobs


def build_closed_model(njobs=1):
    """Delay -> Queue closed loop with a closed class of the given population."""
    model = Network("Closed")
    think = Delay(model, "Think")
    server = Queue(model, "Srv", SchedStrategy.PS)
    jobs = ClosedClass(model, "Job", njobs, think)
    think.setService(jobs, Exp(1.0))
    server.setService(jobs, Exp(2.0))
    model.addLink(think, server)
    model.addLink(server, think)
    return model, server, jobs


class TestApplyMethodsAgainstRealLine:
    """The formerly stubbed apply() methods must alter the real model."""

    @pytest.mark.parametrize("sel_idx, chosen, other", [(0, "Q1", "Q2"),
                                                        (1, "Q2", "Q1")])
    def test_class_service_mapping_routes_to_selected(self, sel_idx, chosen, other):
        model, candidates, jobs = build_parallel_model()
        var = ClassServiceMapping(jobs, stations=candidates)
        evaluator = LineEvaluator(model, [var])

        # Encode the selected station index into the [0, 1] DE coordinate.
        x = np.array([(sel_idx + 0.5) / len(candidates)])
        result = evaluator.evaluate(x)

        # All lambda=2 routed to the chosen station: rho = 2/3; other idle.
        assert np.isclose(result.getUtilization(chosen), 2.0 / 3.0, atol=1e-3)
        assert np.isclose(result.getUtilization(other), 0.0, atol=1e-3)

    def test_job_population_drives_throughput_without_mutating_base(self):
        model, _, jobs = build_closed_model(njobs=1)
        var = JobPopulation(jobs, bounds=(1, 50))
        evaluator = LineEvaluator(model, [var])

        def tput(n):
            x = np.array([(n - 1) / 49.0])
            return evaluator.evaluate(x).getThroughput("Srv", "Job")

        # Closed model: throughput increases monotonically with population.
        assert tput(20) > tput(1)
        # The base model the optimizer holds must be untouched by evaluation.
        assert jobs.getNumberOfJobs() == 1

    @pytest.mark.parametrize("replicas", [4, 6, 10])
    def test_station_replicas_reduce_utilization(self, replicas):
        m = Network("Repl")
        s = Source(m, "S")
        q = Queue(m, "Q", SchedStrategy.FCFS)
        k = Sink(m, "K")
        c = OpenClass(m, "C")
        s.setArrival(c, Exp(3.0))
        q.setService(c, Exp(1.0))
        m.addLink(s, q)
        m.addLink(q, k)

        var = StationReplicas(q, bounds=(4, 10))
        evaluator = LineEvaluator(m, [var])
        x = np.array([(replicas - 4) / 6.0])
        result = evaluator.evaluate(x)

        # N replicas modeled as N servers: util = lambda / (N * mu) = 3 / N.
        assert np.isclose(result.getUtilization("Q"), 3.0 / replicas, atol=1e-3)


def build_tandem_model():
    """Source -> QA -> QB -> Sink tandem with lambda=1, muA=2, muB=3."""
    model = Network("Tandem")
    source = Source(model, "S")
    qa = Queue(model, "QA", SchedStrategy.PS)
    qb = Queue(model, "QB", SchedStrategy.PS)
    sink = Sink(model, "K")
    jobs = OpenClass(model, "Jobs")
    source.setArrival(jobs, Exp(1.0))
    qa.setService(jobs, Exp(2.0))
    qb.setService(jobs, Exp(3.0))
    model.addLink(source, qa)
    model.addLink(qa, qb)
    model.addLink(qb, sink)
    return model, qa, qb, jobs


class TestSystemMetrics:
    """End-to-end (chain-level) metrics must match closed forms."""

    def test_tandem_system_response_time(self):
        model, qa, _, _ = build_tandem_model()
        evaluator = LineEvaluator(model, [ServerAllocation(qa, bounds=(1, 2))])

        result = evaluator.evaluate(np.array([0.0]))

        # M/M/1 tandem: sysRT = 1/(2-1) + 1/(3-1) = 1.5, keyed by class name
        assert np.isclose(result.getSystemResponseTime("Jobs"), 1.5, atol=1e-3)
        # Aggregate (single chain) must agree
        assert np.isclose(result.getSystemResponseTime(), 1.5, atol=1e-3)
        # System throughput equals the arrival rate
        assert np.isclose(result.getSystemThroughput("Jobs"), 1.0, atol=1e-3)

    def test_system_response_time_constraint_drives_sizing(self):
        model, qa, _, jobs = build_tandem_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(qa, bounds=(1, 5)))
        problem.set_objective(MinimizeCost(server_cost={qa: 10.0}))
        # sysRT at 1 server is 1.5 > 1.1, so at least 2 servers are needed
        problem.add_constraint(SystemResponseTimeConstraint(jobs, max_value=1.1))

        result = BisectionSolver(problem).solve()

        assert result.feasible is True
        assert result.variable_values["QA_servers"] == 2


class TestServiceRateAndRoutingApply:
    """Regression tests for the class-identity bugs in apply()."""

    def test_service_rate_takes_effect_on_model_copy(self):
        """ServiceRate.apply must resolve the class in the copied model."""
        model, queue, jobs = build_model()
        evaluator = LineEvaluator(model,
                                  [ServiceRate(queue, jobs, bounds=(1.0, 4.0))])

        result = evaluator.evaluate(np.array([1.0]))  # rate = 4.0

        # M/M/1 with lambda=3, mu=4: RT = 1/(4-3) = 1.0. Before the fix the
        # rate silently stayed at mu=1 and the model was unstable.
        assert np.isclose(result.getResponseTime("Server", "Jobs"), 1.0,
                          atol=1e-3)

    def test_routing_probabilities_split(self):
        model, candidates, jobs = build_parallel_model()
        source = model.get_node_by_name("S")
        var = RoutingProbabilities(jobs, source=source, targets=candidates)
        evaluator = LineEvaluator(model, [var])

        # Stick-breaking with x=0.8: probs = [0.8, 0.2] over arrival rate 2
        result = evaluator.evaluate(np.array([0.8]))

        assert np.isclose(result.getThroughput("Q1", "Jobs"), 1.6, atol=1e-2)
        assert np.isclose(result.getThroughput("Q2", "Jobs"), 0.4, atol=1e-2)


class TestDecompositionCoordination:
    """Subproblems must see values fixed by previously solved blocks."""

    def test_fixed_values_propagate_and_objective_is_real(self):
        model, queue, jobs = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        problem.add_variable(ServiceRate(queue, jobs, bounds=(1.0, 4.0)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0},
                                           rate_cost={queue: 20.0}))
        problem.add_constraint(ResponseTimeConstraint(queue, jobs,
                                                      max_value=0.5))

        workflow = DecompositionWorkflow(problem)
        workflow.auto_decompose()
        workflow.set_solver_options(max_iterations=20, popsize=8, seed=42)
        result = workflow.solve_sequential(max_cycles=2, tolerance=1e-3)

        # Regression: the objective history used to be metric-free/inf
        assert result.objective_history
        assert all(np.isfinite(v) for v in result.objective_history)

        # The service_rate block must have seen the server block's value
        rate_block = result.subproblem_results["service_rate"]
        assert "Server_servers" in rate_block.variables_fixed

        # The coordinated solution must satisfy the SLA when re-evaluated
        evaluator = LineEvaluator(model, problem.getVariables())
        eval_result = evaluator.evaluateValues(result.final_variable_values)
        assert eval_result.feasible
        rt = eval_result.getResponseTime("Server", "Jobs")
        assert rt <= 0.5 + 1e-6


class TestScenarioRobustness:
    """Constraints must hold on the base model and all scenarios."""

    @staticmethod
    def _build_high_load_scenario():
        model = Network("MMcHigh")
        source = Source(model, "Arrivals")
        queue = Queue(model, "Server", SchedStrategy.FCFS)
        sink = Sink(model, "Departures")
        jobs = OpenClass(model, "Jobs")
        source.setArrival(jobs, Exp(4.5))  # higher load than base (3.0)
        queue.setService(jobs, Exp(SERVICE_RATE))
        model.addLink(source, queue)
        model.addLink(queue, sink)
        return model

    def test_worst_case_sizing_with_bisection(self):
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
        problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))
        problem.add_scenario(self._build_high_load_scenario())

        result = BisectionSolver(problem).solve()

        # Base alone needs 6 servers; the 4.5 req/s scenario needs 9
        assert result.feasible is True
        assert result.variable_values["Server_servers"] == 9

    def test_worst_case_sizing_with_de(self):
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
        problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))
        problem.add_scenario(self._build_high_load_scenario())

        result = problem.solve(max_iterations=40, popsize=15, seed=42)

        assert result.feasible is True
        assert result.variable_values["Server_servers"] == 9


class TestBisectionSolver:
    """Bisection must match the DE optimum with O(log n) evaluations."""

    def test_matches_de_optimum_with_few_evaluations(self):
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))
        problem.add_constraint(UtilizationConstraint(queue, max_value=0.5))

        result = BisectionSolver(problem).solve()

        # Same optimum as the DE integration test (6 servers, cost 60)
        assert result.feasible is True
        assert result.variable_values["Server_servers"] == 6
        assert np.isclose(result.objective_value, 60.0)
        # ceil(log2(10)) probes at one LINE solve each
        assert result.model_evaluations <= 6

    def test_max_feasible_population_sizing(self):
        model = Network("Closed")
        think = Delay(model, "Think")
        server = Queue(model, "Srv", SchedStrategy.PS)
        jobs = ClosedClass(model, "Job", 1, think)
        think.setService(jobs, Exp(1.0))
        server.setService(jobs, Exp(2.0))
        model.addLink(think, server)
        model.addLink(server, think)

        problem = OptimizationProblem(model)
        problem.add_variable(JobPopulation(jobs, bounds=(1, 30)))
        problem.set_objective(MinimizeCost())  # objective unused for sizing
        problem.add_constraint(ResponseTimeConstraint(server, jobs,
                                                      max_value=2.0))

        result = BisectionSolver(problem, direction='max_feasible').solve()

        assert result.feasible is True
        n_star = result.variable_values["Job_population"]
        # Verify optimality: n_star feasible, n_star + 1 infeasible
        evaluator = LineEvaluator(model, problem.getVariables())
        rt_at = lambda n: evaluator.evaluateValues(
            {"Job_population": n}).getResponseTime("Srv", "Job")
        assert rt_at(n_star) <= 2.0 + 1e-6
        assert rt_at(n_star + 1) > 2.0

    def test_rejects_multi_variable_problems(self):
        model, queue, jobs = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 10)))
        problem.add_variable(ServiceRate(queue, jobs, bounds=(1.0, 4.0)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))

        with pytest.raises(ValueError):
            BisectionSolver(problem)


class TestParetoSweep:
    """Epsilon-constraint sweep must produce a monotone frontier."""

    def test_utilization_sweep_frontier(self):
        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))

        sweep = ParetoSweep(
            problem,
            constraint_factory=lambda eps: UtilizationConstraint(
                queue, max_value=eps),
            epsilons=[0.3, 0.5, 0.75, 0.9],
            solver='bisection',
        )
        points = sweep.solve()

        assert len(points) == 4
        assert all(p.feasible for p in points)
        # Cost is non-increasing as the utilization bound relaxes
        costs = [p.objective_value for p in points]
        assert costs == sorted(costs, reverse=True)
        # util <= 0.3 needs 10 servers (3/0.3), util <= 0.5 needs 6
        assert points[0].result.variable_values["Server_servers"] == 10
        assert points[1].result.variable_values["Server_servers"] == 6

        # eps=0.9 yields the same cost as eps=0.75, so it is dominated
        frontier = sweep.getFrontier()
        assert [p.epsilon for p in frontier] == [0.3, 0.5, 0.75]

    def test_plot_draws_actual_frontier(self):
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt

        model, queue, _ = build_model()
        problem = OptimizationProblem(model)
        problem.add_variable(ServerAllocation(queue, bounds=(1, 12)))
        problem.set_objective(MinimizeCost(server_cost={queue: 10.0}))

        sweep = ParetoSweep(
            problem,
            constraint_factory=lambda eps: UtilizationConstraint(
                queue, max_value=eps),
            epsilons=[0.3, 0.5, 0.75, 0.9],
            solver='bisection',
        )
        sweep.solve()
        frontier = sweep.getFrontier()

        ax = sweep.plot()
        # The frontier line must trace exactly the non-dominated points,
        # not every swept point: the dominated eps=0.9 is excluded.
        line = ax.get_lines()[0]
        assert list(line.get_xdata()) == [p.epsilon for p in frontier]
        assert list(line.get_ydata()) == [p.objective_value for p in frontier]
        assert 0.9 not in list(line.get_xdata())
        # Drawn as a staircase, not a diagonal interpolation
        assert line.get_drawstyle() == 'steps-post'
        plt.close(ax.figure)


class TestTimeLimit:
    """The solver must stop close to the configured time limit."""

    def test_time_limit_terminates_early(self):
        model, queue, jobs = build_model()
        problem = OptimizationProblem(model)
        # Continuous variable defeats the config cache, so every DE
        # evaluation costs one LINE solve and the deadline hits mid-run
        problem.add_variable(ServiceRate(queue, jobs, bounds=(3.5, 8.0)))
        problem.set_objective(MinimizeCost(rate_cost={queue: 1.0}))

        result = problem.solve(max_iterations=100000, popsize=20,
                               tol=0.0, time_limit=0.3, seed=42)

        assert result.terminated_by == 'time_limit'
        assert result.solve_time < 5.0


class TestLoadBalancing:
    """System RT must follow Little's law on non-serial open topologies."""

    @staticmethod
    def _build(p_fast=None):
        model = Network("LB")
        source = Source(model, "S")
        fast = Queue(model, "Fast", SchedStrategy.PS)
        slow = Queue(model, "Slow", SchedStrategy.PS)
        sink = Sink(model, "K")
        jobs = OpenClass(model, "Jobs")
        source.setArrival(jobs, Exp(2.0))
        fast.setService(jobs, Exp(4.0))
        slow.setService(jobs, Exp(2.5))
        model.addLink(source, fast)
        model.addLink(source, slow)
        model.addLink(fast, sink)
        model.addLink(slow, sink)
        return model, source, fast, slow, jobs

    def test_system_rt_uses_littles_law_on_parallel_topology(self):
        """Regression: LINE's SysRespT sums station RTs unweighted; the
        evaluator must report the Little's law sojourn N/X instead."""
        model, source, fast, slow, jobs = self._build()
        var = RoutingProbabilities(jobs, source=source, targets=[fast, slow])
        evaluator = LineEvaluator(model, [var])

        # Route 80% to Fast: N = 0.4/0.6 + 0.16/0.84, N/X = 0.4286
        result = evaluator.evaluate(np.array([0.8]))

        assert np.isclose(result.getSystemResponseTime("Jobs"), 0.4286,
                          atol=5e-3)

    def test_minimize_system_rt_finds_optimal_split(self):
        from line_solver.opt import MinimizeSystemResponseTime

        model, source, fast, slow, jobs = self._build()
        problem = OptimizationProblem(model)
        problem.add_variable(RoutingProbabilities(jobs, source=source,
                                                  targets=[fast, slow]))
        problem.set_objective(MinimizeSystemResponseTime(jobs))

        result = problem.solve(seed=42, max_iterations=60)

        # Exact M/M/1 optimum: p_fast = 0.743, sysRT = 0.4264
        probs = result.variable_values["Jobs_routing_from_S"]
        assert 0.65 < probs[0] < 0.85
        assert result.objective_value < 0.44
