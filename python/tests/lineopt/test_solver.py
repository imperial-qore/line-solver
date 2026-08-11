"""
Tests for LineOptSolver class in line_solver.opt.solver.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch, PropertyMock

from line_solver.opt.solver import LineOptSolver, LineOptSolverOptions
from line_solver.opt.problem import OptimizationProblem
from line_solver.opt.variables import ServerAllocation
from line_solver.opt.objectives import MinimizeCost
from line_solver.opt.results import EvaluationResult, OptimizationResult


class TestLineOptSolverOptions:
    """Tests for LineOptSolverOptions class."""

    def test_default_options(self):
        """Test default options values."""
        defaults = LineOptSolver.defaultOptions()

        assert defaults['strategy'] == 'best1bin'
        assert defaults['popsize'] == 15
        assert defaults['mutation'] == (0.5, 1.0)
        assert defaults['recombination'] == 0.7
        assert defaults['tol'] == 0.01
        assert defaults['max_iterations'] == 100
        assert defaults['time_limit'] == 300.0
        assert defaults['seed'] is None
        assert defaults['verbose'] is False
        assert defaults['penalty_weight'] == 1e6
        assert defaults['polish'] is False
        assert defaults['workers'] == 1

    def test_options_class_init(self):
        """Test LineOptSolverOptions initialization."""
        options = LineOptSolverOptions()
        d = options.toDict()

        assert d['strategy'] == 'best1bin'
        assert d['verbose'] is False

    def test_set_strategy(self):
        """Test setting strategy."""
        options = LineOptSolverOptions()
        result = options.setStrategy('rand1bin')

        assert result is options  # Returns self for chaining
        assert options.toDict()['strategy'] == 'rand1bin'

    def test_set_popsize(self):
        """Test setting population size."""
        options = LineOptSolverOptions()
        options.setPopsize(20)

        assert options.toDict()['popsize'] == 20

    def test_set_mutation(self):
        """Test setting mutation range."""
        options = LineOptSolverOptions()
        options.setMutation((0.3, 0.9))

        assert options.toDict()['mutation'] == (0.3, 0.9)

    def test_set_recombination(self):
        """Test setting recombination probability."""
        options = LineOptSolverOptions()
        options.setRecombination(0.9)

        assert options.toDict()['recombination'] == 0.9

    def test_set_tolerance(self):
        """Test setting tolerance."""
        options = LineOptSolverOptions()
        options.setTolerance(0.001)

        assert options.toDict()['tol'] == 0.001

    def test_set_max_iterations(self):
        """Test setting max iterations."""
        options = LineOptSolverOptions()
        options.setMaxIterations(200)

        assert options.toDict()['max_iterations'] == 200

    def test_set_time_limit(self):
        """Test setting time limit."""
        options = LineOptSolverOptions()
        options.setTimeLimit(600.0)

        assert options.toDict()['time_limit'] == 600.0

    def test_set_seed(self):
        """Test setting random seed."""
        options = LineOptSolverOptions()
        options.setSeed(42)

        assert options.toDict()['seed'] == 42

    def test_set_verbose(self):
        """Test setting verbose flag."""
        options = LineOptSolverOptions()
        options.setVerbose(True)

        assert options.toDict()['verbose'] is True

    def test_set_penalty_weight(self):
        """Test setting penalty weight."""
        options = LineOptSolverOptions()
        options.setPenaltyWeight(1e8)

        assert options.toDict()['penalty_weight'] == 1e8

    def test_method_chaining(self):
        """Test fluent interface with method chaining."""
        options = (
            LineOptSolverOptions()
            .setStrategy('best2bin')
            .setPopsize(25)
            .setVerbose(True)
            .setSeed(123)
        )

        d = options.toDict()
        assert d['strategy'] == 'best2bin'
        assert d['popsize'] == 25
        assert d['verbose'] is True
        assert d['seed'] == 123

    def test_snake_case_aliases(self):
        """Test snake_case method aliases."""
        options = LineOptSolverOptions()

        options.set_strategy('rand1bin')
        assert options.to_dict()['strategy'] == 'rand1bin'

        options.set_popsize(30)
        options.set_mutation((0.4, 0.8))
        options.set_recombination(0.8)
        options.set_tolerance(0.005)
        options.set_max_iterations(150)
        options.set_time_limit(500.0)
        options.set_seed(99)
        options.set_verbose(True)
        options.set_penalty_weight(1e7)

        d = options.to_dict()
        assert d['popsize'] == 30
        assert d['seed'] == 99


class TestLineOptSolver:
    """Tests for LineOptSolver class."""

    def test_init(self, simple_problem):
        """Test solver initialization."""
        solver = LineOptSolver(simple_problem)

        assert solver.getProblem() is simple_problem
        assert solver.getOptions()['strategy'] == 'best1bin'

    def test_init_with_options(self, simple_problem):
        """Test solver initialization with custom options."""
        solver = LineOptSolver(simple_problem, verbose=True, max_iterations=50)

        options = solver.getOptions()
        assert options['verbose'] is True
        assert options['max_iterations'] == 50

    def test_set_option(self, simple_problem):
        """Test setting individual options."""
        solver = LineOptSolver(simple_problem)
        result = solver.setOption('verbose', True)

        assert result is solver  # Returns self for chaining
        assert solver.getOptions()['verbose'] is True

    def test_snake_case_aliases(self, simple_problem):
        """Test snake_case method aliases."""
        solver = LineOptSolver(simple_problem)

        assert solver.get_problem() is simple_problem
        assert solver.get_options() == solver.getOptions()
        solver.set_option('verbose', True)
        assert solver.get_options()['verbose'] is True

    @patch('line_solver.opt.solver.differential_evolution')
    @patch('line_solver.opt.evaluator.LineEvaluator')
    def test_solve_calls_differential_evolution(
        self, mock_evaluator_class, mock_de, simple_problem
    ):
        """Test that solve calls scipy differential_evolution."""
        # Setup mocks
        mock_evaluator = MagicMock()
        mock_evaluator.getBounds.return_value = [(0.0, 1.0)]
        mock_evaluator.decodeVariables.return_value = {"TestQueue_servers": 5}
        mock_evaluator.getEvaluationCount.return_value = 100

        mock_eval_result = EvaluationResult()
        mock_eval_result.feasible = True
        mock_evaluator.evaluateWithCache.return_value = mock_eval_result

        mock_evaluator_class.return_value = mock_evaluator

        # Mock DE result
        mock_de_result = MagicMock()
        mock_de_result.x = np.array([0.5])
        mock_de_result.fun = 500.0
        mock_de.return_value = mock_de_result

        # Solve
        solver = LineOptSolver(simple_problem, max_iterations=10)
        result = solver.solve()

        # Verify DE was called
        mock_de.assert_called_once()

    def test_build_empty_result(self, mock_network):
        """Test result for problem with no variables."""
        problem = OptimizationProblem(mock_network)
        problem.setObjective(MinimizeCost())

        # Need to add a variable to pass validation, then test empty bounds case
        # This tests _buildEmptyResult indirectly
        # Actually, let's test directly
        solver = LineOptSolver.__new__(LineOptSolver)
        solver._problem = problem
        solver._options = LineOptSolver.defaultOptions()
        solver._start_time = 0

        result = solver._buildEmptyResult()

        assert result.objective_value == 0.0
        assert result.variable_values == {}
        assert result.feasible is True
        assert result.terminated_by == 'empty'


class TestLineOptSolverIntegration:
    """Integration tests for LineOptSolver (without LINE dependency)."""

    @pytest.fixture
    def mock_evaluator_result(self):
        """Create a mock evaluation result."""
        result = EvaluationResult()
        result.feasible = True
        result.response_times = {("TestQueue", "TestClass"): 0.5}
        result.throughputs = {("TestQueue", "TestClass"): 10.0}
        result.utilizations = {"TestQueue": 0.5}
        return result

    @patch('line_solver.opt.evaluator.LineEvaluator')
    def test_objective_function_evaluation(
        self, mock_evaluator_class, simple_problem, mock_evaluator_result
    ):
        """Test objective function is evaluated correctly."""
        mock_evaluator = MagicMock()
        mock_evaluator.getBounds.return_value = [(0.0, 1.0)]
        mock_evaluator.decodeVariables.return_value = {"TestQueue_servers": 5}
        mock_evaluator.evaluateValuesWithCache.return_value = mock_evaluator_result
        mock_evaluator_class.return_value = mock_evaluator

        solver = LineOptSolver(simple_problem)
        solver._evaluator = mock_evaluator
        solver._evaluators = [mock_evaluator]
        solver._caches = [{}]

        # Call objective function directly
        x = np.array([0.5])
        obj_value = solver._objectiveFunction(x)

        # Should be cost = 5 * 100 = 500
        assert obj_value == 500.0

    @patch('line_solver.opt.evaluator.LineEvaluator')
    def test_objective_with_infeasible_result(
        self, mock_evaluator_class, simple_problem
    ):
        """Test objective returns inf for infeasible models."""
        mock_evaluator = MagicMock()
        mock_evaluator.getBounds.return_value = [(0.0, 1.0)]
        mock_evaluator.decodeVariables.return_value = {"TestQueue_servers": 5}

        infeasible_result = EvaluationResult()
        infeasible_result.feasible = False
        mock_evaluator.evaluateValuesWithCache.return_value = infeasible_result
        mock_evaluator_class.return_value = mock_evaluator

        solver = LineOptSolver(simple_problem)
        solver._evaluator = mock_evaluator
        solver._evaluators = [mock_evaluator]
        solver._caches = [{}]

        x = np.array([0.5])
        obj_value = solver._objectiveFunction(x)

        assert obj_value == float('inf')

