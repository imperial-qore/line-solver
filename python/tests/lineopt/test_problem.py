"""
Tests for OptimizationProblem class in line_solver.opt.problem.
"""

import pytest
from unittest.mock import MagicMock, patch

from line_solver.opt.problem import OptimizationProblem
from line_solver.opt.variables import ServerAllocation, ServiceRate
from line_solver.opt.objectives import MinimizeCost, ResponseTimeConstraint


class TestOptimizationProblem:
    """Tests for OptimizationProblem class."""

    def test_init(self, mock_network):
        """Test initialization."""
        problem = OptimizationProblem(mock_network)
        assert problem.getModel() is mock_network
        assert problem.getVariables() == []
        assert problem.getObjective() is None
        assert problem.getConstraints() == []

    def test_add_variable(self, mock_network, mock_station):
        """Test adding a single variable."""
        problem = OptimizationProblem(mock_network)
        var = ServerAllocation(mock_station, bounds=(1, 10))

        result = problem.addVariable(var)

        assert result is problem  # Returns self for chaining
        assert len(problem.getVariables()) == 1
        assert problem.getVariables()[0] is var

    def test_add_variables(self, mock_network, mock_station, mock_station2):
        """Test adding multiple variables."""
        problem = OptimizationProblem(mock_network)
        var1 = ServerAllocation(mock_station, bounds=(1, 10))
        var2 = ServerAllocation(mock_station2, bounds=(1, 5))

        result = problem.addVariables(var1, var2)

        assert result is problem
        assert len(problem.getVariables()) == 2

    def test_set_objective(self, mock_network, mock_station):
        """Test setting objective."""
        problem = OptimizationProblem(mock_network)
        objective = MinimizeCost(server_cost={mock_station: 100.0})

        result = problem.setObjective(objective)

        assert result is problem
        assert problem.getObjective() is objective

    def test_add_constraint(self, mock_network, mock_station):
        """Test adding a constraint."""
        problem = OptimizationProblem(mock_network)
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)

        result = problem.addConstraint(constraint)

        assert result is problem
        assert len(problem.getConstraints()) == 1
        assert problem.getConstraints()[0] is constraint

    def test_add_constraints(self, mock_network, mock_station, mock_station2):
        """Test adding multiple constraints."""
        problem = OptimizationProblem(mock_network)
        c1 = ResponseTimeConstraint(mock_station, max_value=1.0)
        c2 = ResponseTimeConstraint(mock_station2, max_value=2.0)

        result = problem.addConstraints(c1, c2)

        assert result is problem
        assert len(problem.getConstraints()) == 2

    def test_validate_empty_problem(self, mock_network):
        """Test validation of empty problem."""
        problem = OptimizationProblem(mock_network)
        errors = problem.validate()

        assert "No decision variables defined" in errors
        assert "Objective function is not set" in errors

    def test_validate_no_model(self):
        """Test validation with no model."""
        problem = OptimizationProblem(None)
        errors = problem.validate()

        assert "Model is not set" in errors

    def test_validate_valid_problem(self, mock_network, mock_station):
        """Test validation of valid problem."""
        problem = OptimizationProblem(mock_network)
        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))

        errors = problem.validate()
        assert errors == []

    def test_is_valid(self, mock_network, mock_station):
        """Test isValid method."""
        problem = OptimizationProblem(mock_network)
        assert problem.isValid() is False

        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))
        assert problem.isValid() is True

    def test_solve_invalid_problem(self, mock_network):
        """Test solving invalid problem raises error."""
        problem = OptimizationProblem(mock_network)

        with pytest.raises(ValueError, match="Invalid problem"):
            problem.solve()

    def test_summary(self, mock_network, mock_station):
        """Test summary generation."""
        problem = OptimizationProblem(mock_network)
        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))

        summary = problem.summary()

        assert "TestNetwork" in summary
        assert "Variables: 1" in summary
        assert "TestQueue_servers" in summary
        assert "MinimizeCost" in summary

    def test_repr(self, mock_network, mock_station):
        """Test string representation."""
        problem = OptimizationProblem(mock_network)
        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))

        repr_str = repr(problem)

        assert "TestNetwork" in repr_str
        assert "vars=1" in repr_str

    def test_method_chaining(self, mock_network, mock_station, mock_station2):
        """Test fluent interface with method chaining."""
        problem = (
            OptimizationProblem(mock_network)
            .addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
            .addVariable(ServerAllocation(mock_station2, bounds=(1, 5)))
            .setObjective(MinimizeCost(server_cost={mock_station: 100.0}))
            .addConstraint(ResponseTimeConstraint(mock_station, max_value=1.0))
        )

        assert len(problem.getVariables()) == 2
        assert problem.getObjective() is not None
        assert len(problem.getConstraints()) == 1

    def test_snake_case_aliases(self, mock_network, mock_station):
        """Test snake_case method aliases."""
        problem = OptimizationProblem(mock_network)
        var = ServerAllocation(mock_station, bounds=(1, 10))
        objective = MinimizeCost(server_cost={mock_station: 100.0})
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)

        # Test add_variable
        problem.add_variable(var)
        assert len(problem.get_variables()) == 1

        # Test set_objective
        problem.set_objective(objective)
        assert problem.get_objective() is objective

        # Test add_constraint
        problem.add_constraint(constraint)
        assert len(problem.get_constraints()) == 1

        # Test is_valid
        assert problem.is_valid() is True

    def test_decompose(self, mock_network, mock_station):
        """Test decompose method returns DecompositionWorkflow."""
        problem = OptimizationProblem(mock_network)
        problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
        problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))

        workflow = problem.decompose()

        from line_solver.opt.decomposition import DecompositionWorkflow
        assert isinstance(workflow, DecompositionWorkflow)
        assert workflow.getProblem() is problem
