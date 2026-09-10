"""
Tests for decomposition workflow in line_solver.opt.decomposition.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from line_solver.opt.decomposition import SubProblem, DecompositionWorkflow
from line_solver.opt.problem import OptimizationProblem
from line_solver.opt.variables import ServerAllocation, ServiceRate, RoutingProbabilities
from line_solver.opt.objectives import MinimizeCost, ResponseTimeConstraint
from line_solver.opt.results import OptimizationResult, SubProblemResult, WorkflowResult


class TestSubProblem:
    """Tests for SubProblem dataclass."""

    def test_init(self):
        """Test initialization."""
        sp = SubProblem(
            name="test_subproblem",
            variable_type="server_allocation"
        )
        assert sp.name == "test_subproblem"
        assert sp.variable_type == "server_allocation"
        assert sp.variables == []
        assert sp.fixed_values == {}

    def test_get_variable_names(self, mock_station, mock_station2):
        """Test getVariableNames method."""
        var1 = ServerAllocation(mock_station, bounds=(1, 10))
        var2 = ServerAllocation(mock_station2, bounds=(1, 5))

        sp = SubProblem(
            name="servers",
            variable_type="server_allocation",
            variables=[var1, var2]
        )

        names = sp.getVariableNames()
        assert "TestQueue_servers" in names
        assert "TestQueue2_servers" in names

    def test_snake_case_alias(self, mock_station):
        """Test snake_case method alias."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        sp = SubProblem(name="test", variable_type="server_allocation", variables=[var])

        assert sp.get_variable_names() == sp.getVariableNames()


class TestDecompositionWorkflow:
    """Tests for DecompositionWorkflow class."""

    def test_init(self, multi_var_problem):
        """Test initialization."""
        workflow = DecompositionWorkflow(multi_var_problem)

        assert workflow.getProblem() is multi_var_problem
        assert workflow.getSubProblems() == []

    def test_set_solver_options(self, multi_var_problem):
        """Test setting solver options."""
        workflow = DecompositionWorkflow(multi_var_problem)
        result = workflow.setSolverOptions(max_iterations=50, verbose=True)

        assert result is workflow  # Returns self for chaining
        assert workflow._solver_options['max_iterations'] == 50
        assert workflow._solver_options['verbose'] is True

    def test_auto_decompose(self, multi_var_problem):
        """Test automatic decomposition by variable type."""
        workflow = DecompositionWorkflow(multi_var_problem)
        result = workflow.autoDecompose()

        assert result is workflow
        subproblems = workflow.getSubProblems()

        # Should have 2 subproblems: server_allocation and service_rate
        assert len(subproblems) == 2

        # Find server allocation subproblem
        server_sp = next((sp for sp in subproblems if sp.variable_type == 'server_allocation'), None)
        assert server_sp is not None
        assert len(server_sp.variables) == 2  # Two ServerAllocation variables

        # Find service rate subproblem
        rate_sp = next((sp for sp in subproblems if sp.variable_type == 'service_rate'), None)
        assert rate_sp is not None
        assert len(rate_sp.variables) == 1

    def test_auto_decompose_respects_default_order(self, multi_var_problem):
        """Test that auto_decompose respects DEFAULT_ORDER."""
        workflow = DecompositionWorkflow(multi_var_problem)
        workflow.autoDecompose()
        subproblems = workflow.getSubProblems()

        # server_allocation should come before service_rate in DEFAULT_ORDER
        types = [sp.variable_type for sp in subproblems]
        assert types.index('server_allocation') < types.index('service_rate')

    def test_set_dependency(self, multi_var_problem):
        """Test setting explicit dependency."""
        workflow = DecompositionWorkflow(multi_var_problem)
        workflow.autoDecompose()

        result = workflow.setDependency('server_allocation', 'service_rate')

        assert result is workflow
        assert 'server_allocation' in workflow._dependency_graph['service_rate']

    def test_add_subproblem(self, mock_network, mock_station):
        """Test adding custom subproblem."""
        problem = OptimizationProblem(mock_network)
        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))

        workflow = DecompositionWorkflow(problem)
        var = ServerAllocation(mock_station, bounds=(1, 10))

        result = workflow.addSubProblem(
            name="custom_servers",
            variables=[var],
            after=[]
        )

        assert result is workflow
        assert len(workflow.getSubProblems()) == 1
        assert workflow.getSubProblems()[0].name == "custom_servers"

    def test_add_subproblem_with_dependencies(self, mock_network, mock_station, mock_station2):
        """Test adding subproblem with dependencies."""
        problem = OptimizationProblem(mock_network)
        var1 = ServerAllocation(mock_station, bounds=(1, 10))
        var2 = ServerAllocation(mock_station2, bounds=(1, 5))
        problem.addVariables(var1, var2)
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))

        workflow = DecompositionWorkflow(problem)
        workflow.addSubProblem("first", [var1])
        workflow.addSubProblem("second", [var2], after=["first"])

        assert "first" in workflow._dependency_graph["second"]

    def test_get_execution_order_no_dependencies(self, multi_var_problem):
        """Test execution order without explicit dependencies."""
        workflow = DecompositionWorkflow(multi_var_problem)
        workflow.autoDecompose()

        order = workflow._getExecutionOrder()

        # Should return subproblems in the order they were added
        assert len(order) == 2

    def test_snake_case_aliases(self, multi_var_problem):
        """Test snake_case method aliases."""
        workflow = DecompositionWorkflow(multi_var_problem)

        assert workflow.get_problem() is multi_var_problem
        assert workflow.get_subproblems() == workflow.getSubProblems()

        workflow.set_solver_options(verbose=True)
        workflow.auto_decompose()

        assert len(workflow.get_subproblems()) == 2

    @patch('line_solver.opt.decomposition.LineOptSolver')
    def test_solve_sequential_empty(self, mock_solver_class, mock_network):
        """Test solving with no subproblems."""
        problem = OptimizationProblem(mock_network)
        problem.setObjective(MinimizeCost())

        workflow = DecompositionWorkflow(problem)
        result = workflow.solveSequential()

        assert result.converged is True
        assert result.total_solve_time >= 0

    @patch('line_solver.opt.decomposition.LineOptSolver')
    def test_solve_sequential_single_cycle(
        self, mock_solver_class, multi_var_problem
    ):
        """Test sequential solving with single cycle."""
        # Setup mock solver
        mock_solver = MagicMock()
        mock_result = OptimizationResult()
        mock_result.objective_value = 500.0
        mock_result.variable_values = {"TestQueue_servers": 5}
        mock_result.feasible = True
        mock_solver.solve.return_value = mock_result
        mock_solver_class.return_value = mock_solver

        workflow = DecompositionWorkflow(multi_var_problem)
        workflow.autoDecompose()

        result = workflow.solveSequential(max_cycles=1, verbose=False)

        assert isinstance(result, WorkflowResult)
        assert result.cycles_completed >= 0
        assert len(result.subproblem_results) <= 2

    def test_solve_hierarchical(self, multi_var_problem):
        """Test hierarchical solving is single-pass sequential."""
        workflow = DecompositionWorkflow(multi_var_problem)
        workflow.autoDecompose()

        # solveHierarchical should call solveSequential with max_cycles=1
        with patch.object(workflow, 'solveSequential') as mock_seq:
            mock_seq.return_value = WorkflowResult()
            workflow.solveHierarchical(verbose=True)

            mock_seq.assert_called_once_with(
                max_cycles=1,
                tolerance=0.0,
                verbose=True
            )


class TestDecompositionWorkflowDefaultOrder:
    """Tests for DEFAULT_ORDER in DecompositionWorkflow."""

    def test_default_order_values(self):
        """Test default order contains expected variable types.

        Flat-network types come first (in their established order), followed by
        the LayeredNetwork (LQN) variable types.
        """
        expected = [
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
        assert DecompositionWorkflow.DEFAULT_ORDER == expected

    def test_server_allocation_before_routing(self):
        """Test server_allocation comes before routing in default order."""
        order = DecompositionWorkflow.DEFAULT_ORDER
        assert order.index('server_allocation') < order.index('routing')

    def test_service_rate_before_routing(self):
        """Test service_rate comes before routing in default order."""
        order = DecompositionWorkflow.DEFAULT_ORDER
        assert order.index('service_rate') < order.index('routing')


class TestPartialProblemFixedVariables:
    """Partial problems must fix already-solved variables (regression)."""

    def test_partial_problem_fixes_other_variables(self, mock_station,
                                                   mock_jobclass):
        problem = OptimizationProblem(MagicMock())
        server_var = ServerAllocation(mock_station, bounds=(1, 10))
        rate_var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        problem.addVariable(server_var)
        problem.addVariable(rate_var)
        problem.setObjective(MinimizeCost(server_cost={mock_station: 10.0}))

        workflow = DecompositionWorkflow(problem)
        subproblem = SubProblem(name="service_rate",
                                variable_type="service_rate",
                                variables=[rate_var])
        fixed_values = {server_var.getName(): 7}

        partial = workflow._createPartialProblem(subproblem, fixed_values)

        # Only the block's variable is free
        assert partial.getVariables() == [rate_var]
        # The other variable is fixed at its solved value
        fixed = partial.getFixedVariables()
        assert len(fixed) == 1
        assert fixed[0][0] is server_var
        assert fixed[0][1] == 7

    def test_unsolved_variables_are_not_fixed(self, mock_station,
                                              mock_jobclass):
        problem = OptimizationProblem(MagicMock())
        server_var = ServerAllocation(mock_station, bounds=(1, 10))
        rate_var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        problem.addVariable(server_var)
        problem.addVariable(rate_var)
        problem.setObjective(MinimizeCost())

        workflow = DecompositionWorkflow(problem)
        subproblem = SubProblem(name="server_allocation",
                                variable_type="server_allocation",
                                variables=[server_var])

        # First cycle, first block: nothing solved yet
        partial = workflow._createPartialProblem(subproblem, {})
        assert partial.getFixedVariables() == []
