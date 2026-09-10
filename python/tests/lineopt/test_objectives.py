"""
Tests for objective functions and constraints in line_solver.opt.objectives.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock

from line_solver.opt.objectives import (
    Constraint,
    ResponseTimeConstraint,
    ThroughputConstraint,
    UtilizationConstraint,
    BudgetConstraint,
    Objective,
    MinimizeCost,
    MaximizePerformance,
)
from line_solver.opt.results import EvaluationResult


class TestResponseTimeConstraint:
    """Tests for ResponseTimeConstraint."""

    def test_init_with_station_object(self, mock_station):
        """Test initialization with station object."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        assert constraint.getStation() == "TestQueue"
        assert constraint.getMaxValue() == 1.0

    def test_init_with_station_name(self):
        """Test initialization with station name string."""
        constraint = ResponseTimeConstraint("MyQueue", max_value=2.0)
        assert constraint.getStation() == "MyQueue"

    def test_init_with_jobclass(self, mock_station, mock_jobclass):
        """Test initialization with job class."""
        constraint = ResponseTimeConstraint(
            mock_station, jobclass=mock_jobclass, max_value=1.0
        )
        assert constraint.getJobClass() == "TestClass"

    def test_auto_generated_name(self, mock_station):
        """Test auto-generated constraint name."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        assert "RT_TestQueue" in constraint.getName()
        assert "1.0" in constraint.getName()

    def test_custom_name(self, mock_station):
        """Test custom constraint name."""
        constraint = ResponseTimeConstraint(
            mock_station, max_value=1.0, name="custom_rt"
        )
        assert constraint.getName() == "custom_rt"

    def test_evaluate_satisfied(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is satisfied."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        # TestQueue response time is 0.5 (aggregate of class RTs)
        violation = constraint.evaluate(sample_evaluation_result, {})
        assert violation == 0.0

    def test_evaluate_violated(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is violated."""
        constraint = ResponseTimeConstraint(mock_station, max_value=0.3)
        violation = constraint.evaluate(sample_evaluation_result, {})
        # Aggregate RT is (0.5 + 0.8) / 2 = 0.65, violation is 0.65 - 0.3 = 0.35
        assert violation > 0

    def test_is_satisfied(self, mock_station, sample_evaluation_result):
        """Test isSatisfied method."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        assert constraint.isSatisfied(sample_evaluation_result, {})

        constraint2 = ResponseTimeConstraint(mock_station, max_value=0.1)
        assert not constraint2.isSatisfied(sample_evaluation_result, {})

    def test_snake_case_aliases(self, mock_station):
        """Test snake_case method aliases."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        assert constraint.get_name() == constraint.getName()
        assert constraint.get_station() == constraint.getStation()
        assert constraint.get_max_value() == constraint.getMaxValue()


class TestThroughputConstraint:
    """Tests for ThroughputConstraint."""

    def test_init(self, mock_station):
        """Test initialization."""
        constraint = ThroughputConstraint(mock_station, min_value=5.0)
        assert constraint.getStation() == "TestQueue"
        assert constraint.getMinValue() == 5.0

    def test_evaluate_satisfied(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is satisfied."""
        constraint = ThroughputConstraint(mock_station, min_value=10.0)
        # TestQueue throughput is 10.0 + 5.0 = 15.0 (aggregate)
        violation = constraint.evaluate(sample_evaluation_result, {})
        assert violation == 0.0

    def test_evaluate_violated(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is violated."""
        constraint = ThroughputConstraint(mock_station, min_value=20.0)
        violation = constraint.evaluate(sample_evaluation_result, {})
        # Aggregate is 15.0, violation is 20.0 - 15.0 = 5.0
        assert np.isclose(violation, 5.0)


class TestUtilizationConstraint:
    """Tests for UtilizationConstraint."""

    def test_init(self, mock_station):
        """Test initialization."""
        constraint = UtilizationConstraint(mock_station, max_value=0.8)
        assert constraint.getStation() == "TestQueue"
        assert constraint.getMaxValue() == 0.8

    def test_evaluate_satisfied(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is satisfied."""
        constraint = UtilizationConstraint(mock_station, max_value=0.8)
        # TestQueue utilization is 0.6
        violation = constraint.evaluate(sample_evaluation_result, {})
        assert violation == 0.0

    def test_evaluate_violated(self, mock_station, sample_evaluation_result):
        """Test evaluation when constraint is violated."""
        constraint = UtilizationConstraint(mock_station, max_value=0.5)
        violation = constraint.evaluate(sample_evaluation_result, {})
        # Utilization is 0.6, violation is 0.6 - 0.5 = 0.1
        assert np.isclose(violation, 0.1)


class TestBudgetConstraint:
    """Tests for BudgetConstraint."""

    def test_init(self):
        """Test initialization."""
        constraint = BudgetConstraint(
            budget=500.0,
            cost_coefficients={"Queue_servers": 50.0}
        )
        assert constraint.getBudget() == 500.0
        assert constraint.getCostCoefficients() == {"Queue_servers": 50.0}

    def test_compute_cost(self):
        """Test cost computation."""
        constraint = BudgetConstraint(
            budget=500.0,
            cost_coefficients={"var1": 10.0, "var2": 20.0}
        )
        cost = constraint.computeCost({"var1": 5, "var2": 3})
        assert cost == 10.0 * 5 + 20.0 * 3  # 110

    def test_evaluate_satisfied(self, sample_evaluation_result):
        """Test evaluation when constraint is satisfied."""
        constraint = BudgetConstraint(
            budget=500.0,
            cost_coefficients={"servers": 50.0}
        )
        violation = constraint.evaluate(sample_evaluation_result, {"servers": 5})
        assert violation == 0.0  # 5 * 50 = 250 < 500

    def test_evaluate_violated(self, sample_evaluation_result):
        """Test evaluation when constraint is violated."""
        constraint = BudgetConstraint(
            budget=100.0,
            cost_coefficients={"servers": 50.0}
        )
        violation = constraint.evaluate(sample_evaluation_result, {"servers": 5})
        assert violation == 150.0  # 5 * 50 = 250, violation = 250 - 100 = 150

    def test_snake_case_aliases(self):
        """Test snake_case method aliases."""
        constraint = BudgetConstraint(budget=500.0)
        assert constraint.get_budget() == constraint.getBudget()
        assert constraint.compute_cost({}) == constraint.computeCost({})


class TestMinimizeCost:
    """Tests for MinimizeCost objective."""

    def test_init(self, mock_station):
        """Test initialization."""
        objective = MinimizeCost(
            server_cost={mock_station: 100.0}
        )
        assert objective.isMinimization() is True
        server_cost = objective.getServerCost()
        assert "TestQueue_servers" in server_cost
        assert server_cost["TestQueue_servers"] == 100.0

    def test_init_with_constraints(self, mock_station):
        """Test initialization with constraints."""
        constraint = ResponseTimeConstraint(mock_station, max_value=1.0)
        objective = MinimizeCost(
            server_cost={mock_station: 100.0},
            subject_to=[constraint]
        )
        constraints = objective.getConstraints()
        assert len(constraints) == 1
        assert constraints[0] is constraint

    def test_evaluate(self, mock_station, sample_evaluation_result):
        """Test objective evaluation."""
        objective = MinimizeCost(
            server_cost={mock_station: 100.0}
        )
        cost = objective.evaluate(
            sample_evaluation_result,
            {"TestQueue_servers": 5}
        )
        assert cost == 500.0  # 5 * 100

    def test_evaluate_with_penalty(self, mock_station, sample_evaluation_result):
        """Test objective evaluation with constraint penalty."""
        constraint = ResponseTimeConstraint(mock_station, max_value=0.1)
        objective = MinimizeCost(
            server_cost={mock_station: 100.0},
            subject_to=[constraint]
        )
        value = objective.evaluateWithPenalty(
            sample_evaluation_result,
            {"TestQueue_servers": 5},
            penalty_weight=1000.0
        )
        # Cost is 500, plus penalty for violated constraint
        assert value > 500.0


class TestMaximizePerformance:
    """Tests for MaximizePerformance objective."""

    def test_init(self):
        """Test initialization."""
        objective = MaximizePerformance(
            throughput_weight=1.0,
            response_time_weight=2.0
        )
        assert objective.isMinimization() is True  # Minimizes negative performance
        assert objective.getThroughputWeight() == 1.0
        assert objective.getResponseTimeWeight() == 2.0

    def test_init_with_budget(self):
        """Test initialization with budget constraint."""
        objective = MaximizePerformance(
            budget=500.0,
            budget_terms={"servers": 50.0}
        )
        constraints = objective.getConstraints()
        assert len(constraints) == 1
        assert isinstance(constraints[0], BudgetConstraint)

    def test_evaluate(self, sample_evaluation_result):
        """Test objective evaluation."""
        objective = MaximizePerformance(
            throughput_weight=1.0,
            response_time_weight=0.0,
            stations=["TestQueue"]
        )
        value = objective.evaluate(sample_evaluation_result, {})
        # Returns negative of throughput (for minimization)
        assert value < 0

    def test_snake_case_aliases(self):
        """Test snake_case method aliases."""
        objective = MaximizePerformance()
        assert objective.get_throughput_weight() == objective.getThroughputWeight()
        assert objective.get_response_time_weight() == objective.getResponseTimeWeight()
        assert objective.is_minimization() == objective.isMinimization()
        assert objective.get_constraints() == objective.getConstraints()


class TestSystemResponseTimeConstraint:
    """Tests for SystemResponseTimeConstraint."""

    def test_satisfied(self):
        from line_solver.opt.objectives import SystemResponseTimeConstraint
        result = EvaluationResult()
        result.system_response_times = {"C": 0.8}
        result.system_throughputs = {"C": 1.0}

        constraint = SystemResponseTimeConstraint("C", max_value=1.0)
        assert constraint.evaluate(result, {}) == 0.0
        assert constraint.isSatisfied(result, {})

    def test_violated(self):
        from line_solver.opt.objectives import SystemResponseTimeConstraint
        result = EvaluationResult()
        result.system_response_times = {"C": 1.6}

        constraint = SystemResponseTimeConstraint("C", max_value=1.0)
        assert np.isclose(constraint.evaluate(result, {}), 0.6)

    def test_aggregate_when_class_is_none(self):
        from line_solver.opt.objectives import SystemResponseTimeConstraint
        result = EvaluationResult()
        result.system_response_times = {"A": 1.5, "B": 3.0}
        result.system_throughputs = {"A": 2.0, "B": 1.0}

        # Weighted aggregate is 2.0
        constraint = SystemResponseTimeConstraint(max_value=1.5)
        assert np.isclose(constraint.evaluate(result, {}), 0.5)

    def test_missing_metrics_are_infeasible(self):
        from line_solver.opt.objectives import SystemResponseTimeConstraint
        constraint = SystemResponseTimeConstraint("C", max_value=1.0)
        assert constraint.evaluate(EvaluationResult(), {}) == float('inf')

    def test_default_name(self):
        from line_solver.opt.objectives import SystemResponseTimeConstraint
        constraint = SystemResponseTimeConstraint("C", max_value=2.0)
        assert constraint.getName() == "SysRT_C_le_2.0"


class TestMinimizeSystemResponseTime:
    """Tests for the MinimizeSystemResponseTime objective."""

    def test_evaluate_reads_system_metric(self):
        from line_solver.opt.objectives import MinimizeSystemResponseTime
        result = EvaluationResult()
        result.system_response_times = {"C": 0.42}
        result.system_throughputs = {"C": 2.0}

        objective = MinimizeSystemResponseTime("C")
        assert objective.isMinimization()
        assert np.isclose(objective.evaluate(result, {}), 0.42)

    def test_missing_metric_is_inf(self):
        from line_solver.opt.objectives import MinimizeSystemResponseTime
        objective = MinimizeSystemResponseTime("C")
        assert objective.evaluate(EvaluationResult(), {}) == float('inf')

    def test_subject_to_constraints_attach(self):
        from line_solver.opt.objectives import (MinimizeSystemResponseTime,
                                         UtilizationConstraint)
        constraint = UtilizationConstraint("Q", max_value=0.8)
        objective = MinimizeSystemResponseTime("C", subject_to=[constraint])
        assert objective.getConstraints() == [constraint]
