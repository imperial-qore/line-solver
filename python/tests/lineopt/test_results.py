"""
Tests for result dataclasses in line_solver.opt.results.
"""

import pytest
import numpy as np

from line_solver.opt.results import (
    EvaluationResult,
    OptimizationResult,
    SubProblemResult,
    WorkflowResult,
)


class TestEvaluationResult:
    """Tests for EvaluationResult dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        result = EvaluationResult()
        assert result.feasible is True
        assert result.response_times == {}
        assert result.throughputs == {}
        assert result.utilizations == {}
        assert result.queue_lengths == {}
        assert result.solve_time == 0.0
        assert result.solver_used == ""

    def test_get_response_time_with_class(self, sample_evaluation_result):
        """Test getting response time for specific class."""
        rt = sample_evaluation_result.getResponseTime("TestQueue", "TestClass")
        assert rt == 0.5

    def test_get_response_time_aggregate(self, sample_evaluation_result):
        """Test getting aggregate response time."""
        rt = sample_evaluation_result.getResponseTime("TestQueue")
        # (0.5 + 0.8) / 2 = 0.65
        assert np.isclose(rt, 0.65)

    def test_get_response_time_missing(self, sample_evaluation_result):
        """Test getting response time for missing station."""
        rt = sample_evaluation_result.getResponseTime("NonExistent")
        assert rt == float('inf')

    def test_get_throughput_with_class(self, sample_evaluation_result):
        """Test getting throughput for specific class."""
        tput = sample_evaluation_result.getThroughput("TestQueue", "TestClass")
        assert tput == 10.0

    def test_get_throughput_aggregate(self, sample_evaluation_result):
        """Test getting aggregate throughput."""
        tput = sample_evaluation_result.getThroughput("TestQueue")
        # 10.0 + 5.0 = 15.0
        assert tput == 15.0

    def test_get_throughput_missing(self, sample_evaluation_result):
        """Test getting throughput for missing station."""
        tput = sample_evaluation_result.getThroughput("NonExistent")
        assert tput == 0.0

    def test_get_utilization(self, sample_evaluation_result):
        """Test getting utilization."""
        util = sample_evaluation_result.getUtilization("TestQueue")
        assert util == 0.6

    def test_get_utilization_missing(self, sample_evaluation_result):
        """Test getting utilization for missing station."""
        util = sample_evaluation_result.getUtilization("NonExistent")
        assert util == 0.0

    def test_get_queue_length_with_class(self, sample_evaluation_result):
        """Test getting queue length for specific class."""
        qlen = sample_evaluation_result.getQueueLength("TestQueue", "TestClass")
        assert qlen == 2.5

    def test_get_queue_length_aggregate(self, sample_evaluation_result):
        """Test getting aggregate queue length."""
        qlen = sample_evaluation_result.getQueueLength("TestQueue")
        # 2.5 + 1.5 = 4.0
        assert qlen == 4.0

    def test_snake_case_aliases(self, sample_evaluation_result):
        """Test snake_case method aliases."""
        assert sample_evaluation_result.get_response_time("TestQueue") == \
            sample_evaluation_result.getResponseTime("TestQueue")
        assert sample_evaluation_result.get_throughput("TestQueue") == \
            sample_evaluation_result.getThroughput("TestQueue")
        assert sample_evaluation_result.get_utilization("TestQueue") == \
            sample_evaluation_result.getUtilization("TestQueue")
        assert sample_evaluation_result.get_queue_length("TestQueue") == \
            sample_evaluation_result.getQueueLength("TestQueue")


class TestOptimizationResult:
    """Tests for OptimizationResult dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        result = OptimizationResult()
        assert result.objective_value == float('inf')
        assert result.variable_values == {}
        assert result.constraint_violations == {}
        assert result.feasible is False
        assert result.iterations == 0
        assert result.solve_time == 0.0
        assert result.model_evaluations == 0
        assert result.convergence_history == []
        assert result.terminated_by == ""

    def test_is_feasible(self):
        """Test isFeasible method."""
        result = OptimizationResult()
        assert result.isFeasible() is False

        result.feasible = True
        assert result.isFeasible() is True

    def test_get_objective_value(self):
        """Test getObjectiveValue method."""
        result = OptimizationResult()
        result.objective_value = 42.0
        assert result.getObjectiveValue() == 42.0

    def test_get_variable_value(self):
        """Test getVariableValue method."""
        result = OptimizationResult()
        result.variable_values = {"var1": 5, "var2": 10}
        assert result.getVariableValue("var1") == 5
        assert result.getVariableValue("var2") == 10
        assert result.getVariableValue("var3") is None

    def test_get_constraint_violation(self):
        """Test getConstraintViolation method."""
        result = OptimizationResult()
        result.constraint_violations = {"rt_constraint": 0.5}
        assert result.getConstraintViolation("rt_constraint") == 0.5
        assert result.getConstraintViolation("other") == 0.0

    def test_get_total_violation(self):
        """Test getTotalViolation method."""
        result = OptimizationResult()
        result.constraint_violations = {"c1": 0.5, "c2": 1.0, "c3": 0.25}
        assert result.getTotalViolation() == 1.75

    def test_repr(self):
        """Test string representation."""
        result = OptimizationResult()
        result.objective_value = 100.0
        result.feasible = True
        result.iterations = 50
        result.model_evaluations = 1000
        result.solve_time = 5.5

        repr_str = repr(result)
        assert "100.0" in repr_str
        assert "feasible" in repr_str
        assert "50" in repr_str

    def test_snake_case_aliases(self):
        """Test snake_case method aliases."""
        result = OptimizationResult()
        result.objective_value = 42.0
        result.feasible = True

        assert result.is_feasible() == result.isFeasible()
        assert result.get_objective_value() == result.getObjectiveValue()
        assert result.get_variable_value("x") == result.getVariableValue("x")
        assert result.get_constraint_violation("c") == result.getConstraintViolation("c")
        assert result.get_total_violation() == result.getTotalViolation()


class TestSubProblemResult:
    """Tests for SubProblemResult dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        result = SubProblemResult()
        assert result.name == ""
        assert isinstance(result.result, OptimizationResult)
        assert result.variables_fixed == {}

    def test_with_values(self):
        """Test initialization with values."""
        opt_result = OptimizationResult()
        opt_result.objective_value = 50.0

        result = SubProblemResult(
            name="server_allocation",
            result=opt_result,
            variables_fixed={"other_var": 5}
        )
        assert result.name == "server_allocation"
        assert result.result.objective_value == 50.0
        assert result.variables_fixed == {"other_var": 5}


class TestWorkflowResult:
    """Tests for WorkflowResult dataclass."""

    def test_default_values(self):
        """Test default initialization."""
        result = WorkflowResult()
        assert result.final_objective == float('inf')
        assert result.subproblem_results == {}
        assert result.cycles_completed == 0
        assert result.converged is False
        assert result.total_solve_time == 0.0
        assert result.objective_history == []
        assert result.final_variable_values == {}

    def test_is_converged(self):
        """Test isConverged method."""
        result = WorkflowResult()
        assert result.isConverged() is False

        result.converged = True
        assert result.isConverged() is True

    def test_get_subproblem_result(self):
        """Test getSubProblemResult method."""
        result = WorkflowResult()
        sp_result = SubProblemResult(name="test")
        result.subproblem_results = {"test": sp_result}

        assert result.getSubProblemResult("test") is sp_result
        assert result.getSubProblemResult("nonexistent") is None

    def test_get_final_variable_value(self):
        """Test getFinalVariableValue method."""
        result = WorkflowResult()
        result.final_variable_values = {"var1": 10, "var2": 20}

        assert result.getFinalVariableValue("var1") == 10
        assert result.getFinalVariableValue("var3") is None

    def test_repr(self):
        """Test string representation."""
        result = WorkflowResult()
        result.final_objective = 75.0
        result.subproblem_results = {"sp1": SubProblemResult(), "sp2": SubProblemResult()}
        result.cycles_completed = 3
        result.converged = True
        result.total_solve_time = 10.5

        repr_str = repr(result)
        assert "75.0" in repr_str
        assert "2 subproblems" in repr_str
        assert "3 cycles" in repr_str
        assert "converged" in repr_str

    def test_snake_case_aliases(self):
        """Test snake_case method aliases."""
        result = WorkflowResult()
        assert result.is_converged() == result.isConverged()
        assert result.get_subproblem_result("x") == result.getSubProblemResult("x")
        assert result.get_final_variable_value("x") == result.getFinalVariableValue("x")


class TestSystemMetricsAccessors:
    """Tests for system (end-to-end) metric getters."""

    def test_get_system_response_time_by_class(self):
        result = EvaluationResult()
        result.system_response_times = {"A": 1.5, "B": 3.0}
        result.system_throughputs = {"A": 2.0, "B": 1.0}

        assert result.getSystemResponseTime("A") == 1.5
        assert result.getSystemResponseTime("B") == 3.0
        assert result.getSystemResponseTime("missing") == float('inf')

    def test_get_system_response_time_aggregate_is_tput_weighted(self):
        result = EvaluationResult()
        result.system_response_times = {"A": 1.5, "B": 3.0}
        result.system_throughputs = {"A": 2.0, "B": 1.0}

        # (1.5*2 + 3.0*1) / 3 = 2.0
        assert np.isclose(result.getSystemResponseTime(), 2.0)

    def test_get_system_response_time_empty_is_inf(self):
        result = EvaluationResult()
        assert result.getSystemResponseTime() == float('inf')

    def test_get_system_throughput(self):
        result = EvaluationResult()
        result.system_throughputs = {"A": 2.0, "B": 1.0}

        assert result.getSystemThroughput("A") == 2.0
        assert result.getSystemThroughput() == 3.0
        assert result.getSystemThroughput("missing") == 0.0
