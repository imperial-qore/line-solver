"""
Tests for LineEvaluator class in line_solver.opt.evaluator.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, patch

from line_solver.opt.evaluator import LineEvaluator
from line_solver.opt.variables import ServerAllocation, ServiceRate
from line_solver.opt.results import EvaluationResult


class TestLineEvaluator:
    """Tests for LineEvaluator class."""

    def test_init(self, mock_network, variables):
        """Test evaluator initialization."""
        evaluator = LineEvaluator(mock_network, variables)

        assert evaluator.getModel() is mock_network
        assert evaluator.getVariables() == variables
        assert evaluator.getEvaluationCount() == 0

    def test_get_total_dimension(self, mock_network, variables):
        """Test total dimension calculation."""
        evaluator = LineEvaluator(mock_network, variables)

        # ServerAllocation has dim 1, ServiceRate has dim 1
        assert evaluator.getTotalDimension() == 2

    def test_get_bounds(self, mock_network, variables):
        """Test bounds aggregation."""
        evaluator = LineEvaluator(mock_network, variables)
        bounds = evaluator.getBounds()

        assert len(bounds) == 2
        assert bounds[0] == (0.0, 1.0)
        assert bounds[1] == (0.0, 1.0)

    def test_decode_variables(self, mock_network, variables):
        """Test variable decoding."""
        evaluator = LineEvaluator(mock_network, variables)

        x = np.array([0.5, 0.5])
        values = evaluator.decodeVariables(x)

        assert "TestQueue_servers" in values
        assert "TestQueue_TestClass_rate" in values
        # ServerAllocation: 1 + 0.5 * 9 = 5.5 -> 6
        assert values["TestQueue_servers"] == 6
        # ServiceRate: 0.5 + 0.5 * 4.5 = 2.75
        assert np.isclose(values["TestQueue_TestClass_rate"], 2.75)

    def test_decode_variables_at_bounds(self, mock_network, variables):
        """Test decoding at boundary values."""
        evaluator = LineEvaluator(mock_network, variables)

        # At minimum
        x_min = np.array([0.0, 0.0])
        values_min = evaluator.decodeVariables(x_min)
        assert values_min["TestQueue_servers"] == 1
        assert np.isclose(values_min["TestQueue_TestClass_rate"], 0.5)

        # At maximum
        x_max = np.array([1.0, 1.0])
        values_max = evaluator.decodeVariables(x_max)
        assert values_max["TestQueue_servers"] == 10
        assert np.isclose(values_max["TestQueue_TestClass_rate"], 5.0)

    def test_apply_variables(self, mock_network, variables):
        """Test applying variables to model."""
        evaluator = LineEvaluator(mock_network, variables)

        values = {
            "TestQueue_servers": 5,
            "TestQueue_TestClass_rate": 2.0,
        }

        # Apply to a mock model
        mock_model = MagicMock()
        mock_station = MagicMock()
        mock_station.getName.return_value = "TestQueue"
        mock_model.getNodes.return_value = [mock_station]

        evaluator.applyVariables(mock_model, values)

        # Server allocation should have been applied
        mock_station.setNumberOfServers.assert_called_once_with(5)

    def test_copy_model(self, mock_network, variables):
        """Test model copying."""
        evaluator = LineEvaluator(mock_network, variables)

        # This will attempt deepcopy, may fail with mocks
        # but should return something
        model_copy = evaluator.copyModel()
        assert model_copy is not None

    def test_evaluate_without_line(self, mock_network, variables):
        """Test evaluation when LINE is not available."""
        evaluator = LineEvaluator(mock_network, variables)
        evaluator._line_available = False

        x = np.array([0.5, 0.5])
        result = evaluator.evaluate(x)

        assert result.feasible is False
        assert evaluator.getEvaluationCount() == 1

    def test_evaluate_with_cache_miss(self, mock_network, variables):
        """Test cached evaluation with cache miss."""
        evaluator = LineEvaluator(mock_network, variables)
        evaluator._line_available = False

        x = np.array([0.5, 0.5])
        cache = {}

        result = evaluator.evaluateWithCache(x, cache)

        # Should be in cache now
        assert len(cache) == 1
        assert result.feasible is False

    def test_evaluate_with_cache_hit(self, mock_network, variables):
        """Test cached evaluation with cache hit."""
        evaluator = LineEvaluator(mock_network, variables)

        x = np.array([0.5, 0.5])
        cached_result = EvaluationResult()
        cached_result.feasible = True
        cached_result.response_times = {("TestQueue", "TestClass"): 0.3}

        # Cache is keyed on the decoded configuration, not on x
        key = evaluator._valuesKey(evaluator.decodeVariables(x))
        cache = {key: cached_result}

        # Should return cached result without calling evaluate
        original_count = evaluator.getEvaluationCount()
        result = evaluator.evaluateWithCache(x, cache)

        assert result is cached_result
        assert evaluator.getEvaluationCount() == original_count

    def test_cache_collapses_rounding_plateau(self, mock_network, mock_station):
        """Distinct x mapping to the same integer config must share one entry."""
        evaluator = LineEvaluator(mock_network,
                                  [ServerAllocation(mock_station, bounds=(1, 10))])
        cache = {}
        # Both encode to 6 servers: 1 + x*9 rounds to 6 for x in [0.5, 0.6]
        evaluator.evaluateWithCache(np.array([0.556]), cache)
        evaluator.evaluateWithCache(np.array([0.60]), cache)

        assert len(cache) == 1
        assert evaluator.getEvaluationCount() == 1

    def test_snake_case_aliases(self, mock_network, variables):
        """Test snake_case method aliases."""
        evaluator = LineEvaluator(mock_network, variables)

        assert evaluator.get_model() is mock_network
        assert evaluator.get_variables() == variables
        assert evaluator.get_total_dimension() == evaluator.getTotalDimension()
        assert evaluator.get_evaluation_count() == evaluator.getEvaluationCount()
        assert evaluator.get_bounds() == evaluator.getBounds()

        x = np.array([0.5, 0.5])
        assert evaluator.decode_variables(x) == evaluator.decodeVariables(x)


class TestLineEvaluatorMetricExtraction:
    """Tests for metric extraction in LineEvaluator."""

    @pytest.fixture
    def avg_table(self):
        """
        Build an AvgTable-shaped object matching SolverAuto.getAvgTable().

        Mirrors LINE's IndexedTable: a DataFrame keyed by Station/JobClass with
        columns RespT, Tput, Util, QLen. Wrapped in a stub exposing ``.data``.
        """
        pd = pytest.importorskip("pandas")

        df = pd.DataFrame({
            "Station": ["Queue1", "Queue1", "Queue2", "Queue2"],
            "JobClass": ["Class1", "Class2", "Class1", "Class2"],
            "QLen": [2.5, 1.5, 3.0, 2.0],
            "Util": [0.6, 0.3, 0.8, 0.4],
            "RespT": [0.5, 0.8, 1.0, 1.2],
            "Tput": [10.0, 5.0, 8.0, 4.0],
        })

        class _IndexedTableStub:
            def __init__(self, data):
                self.data = data

        return _IndexedTableStub(df)

    def test_extract_metrics(self, avg_table, mock_station):
        """Test metric extraction from the AvgTable DataFrame."""
        evaluator = LineEvaluator.__new__(LineEvaluator)
        evaluator._variables = [ServerAllocation(mock_station, bounds=(1, 10))]

        result = EvaluationResult()
        evaluator._extractMetrics(avg_table, result)

        # Check response times were extracted
        assert ("Queue1", "Class1") in result.response_times
        assert result.response_times[("Queue1", "Class1")] == 0.5

        # Check throughputs
        assert ("Queue1", "Class1") in result.throughputs
        assert result.throughputs[("Queue1", "Class1")] == 10.0

        # Utilization is aggregated across classes per station
        assert "Queue1" in result.utilizations
        assert np.isclose(result.utilizations["Queue1"], 0.9)

        # Check queue lengths
        assert ("Queue1", "Class1") in result.queue_lengths
        assert result.queue_lengths[("Queue1", "Class1")] == 2.5

    def test_extract_metrics_accepts_bare_dataframe(self, mock_station):
        """Extraction should also accept a raw DataFrame (no .data wrapper)."""
        pd = pytest.importorskip("pandas")
        df = pd.DataFrame({
            "Station": ["Q"], "JobClass": ["C"],
            "QLen": [1.0], "Util": [0.5], "RespT": [0.25], "Tput": [2.0],
        })
        evaluator = LineEvaluator.__new__(LineEvaluator)
        evaluator._variables = [ServerAllocation(mock_station, bounds=(1, 10))]

        result = EvaluationResult()
        evaluator._extractMetrics(df, result)

        assert result.response_times[("Q", "C")] == 0.25
        assert result.throughputs[("Q", "C")] == 2.0
