"""
Tests for decision variable classes in line_solver.opt.variables.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock

from line_solver.opt.variables import (
    DecisionVariable,
    ServerAllocation,
    StationReplicas,
    RoutingProbabilities,
    ClassServiceMapping,
    ServiceRate,
    JobPopulation,
    ClassPriority,
)


class TestServerAllocation:
    """Tests for ServerAllocation decision variable."""

    def test_init_default_name(self, mock_station):
        """Test initialization with default name."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.getName() == "TestQueue_servers"
        assert var.getDimension() == 1

    def test_init_custom_name(self, mock_station):
        """Test initialization with custom name."""
        var = ServerAllocation(mock_station, bounds=(1, 10), name="custom_servers")
        assert var.getName() == "custom_servers"

    def test_bounds(self, mock_station):
        """Test bounds are in [0, 1]."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        bounds = var.getBounds()
        assert len(bounds) == 1
        assert bounds[0] == (0.0, 1.0)

    def test_decode_min(self, mock_station):
        """Test decoding at minimum value."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        result = var.decode(np.array([0.0]))
        assert result == 1

    def test_decode_max(self, mock_station):
        """Test decoding at maximum value."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        result = var.decode(np.array([1.0]))
        assert result == 10

    def test_decode_middle(self, mock_station):
        """Test decoding at middle value."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        result = var.decode(np.array([0.5]))
        # (1 + 0.5 * 9) = 5.5, rounds to 6
        assert result == 6

    def test_decode_rounding(self, mock_station):
        """Test proper rounding behavior."""
        var = ServerAllocation(mock_station, bounds=(1, 5))
        # 1 + 0.3 * 4 = 2.2, rounds to 2
        assert var.decode(np.array([0.3])) == 2
        # 1 + 0.7 * 4 = 3.8, rounds to 4
        assert var.decode(np.array([0.7])) == 4

    def test_variable_type(self, mock_station):
        """Test variable type identifier."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.getVariableType() == "server_allocation"

    def test_get_station(self, mock_station):
        """Test getStation method."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.getStation() is mock_station

    def test_min_max_servers(self, mock_station):
        """Test min/max server accessors."""
        var = ServerAllocation(mock_station, bounds=(2, 8))
        assert var.getMinServers() == 2
        assert var.getMaxServers() == 8

    def test_apply(self, mock_station, mock_network):
        """Test applying value to model."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        var.apply(mock_network, 5)
        # Check the station's servers were updated
        assert mock_station.getNumberOfServers() == 5

    def test_snake_case_aliases(self, mock_station):
        """Test snake_case method aliases."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.get_name() == var.getName()
        assert var.get_dimension() == var.getDimension()
        assert var.get_bounds() == var.getBounds()
        assert var.get_variable_type() == var.getVariableType()
        assert var.get_station() == var.getStation()


class TestStationReplicas:
    """Tests for StationReplicas decision variable."""

    def test_init_default_name(self, mock_station):
        """Test initialization with default name."""
        var = StationReplicas(mock_station, bounds=(1, 5))
        assert var.getName() == "TestQueue_replicas"

    def test_decode(self, mock_station):
        """Test decoding to replica count."""
        var = StationReplicas(mock_station, bounds=(1, 5))
        assert var.decode(np.array([0.0])) == 1
        assert var.decode(np.array([1.0])) == 5
        assert var.decode(np.array([0.5])) == 3

    def test_variable_type(self, mock_station):
        """Test variable type identifier."""
        var = StationReplicas(mock_station, bounds=(1, 5))
        assert var.getVariableType() == "station_replicas"


class TestRoutingProbabilities:
    """Tests for RoutingProbabilities decision variable."""

    def test_init(self, mock_jobclass, mock_source, mock_station, mock_station2):
        """Test initialization."""
        var = RoutingProbabilities(
            mock_jobclass,
            source=mock_source,
            targets=[mock_station, mock_station2]
        )
        assert var.getName() == "TestClass_routing_from_Source"
        # dimension is len(targets) - 1 for stick-breaking
        assert var.getDimension() == 1

    def test_decode_single_target(self, mock_jobclass, mock_source, mock_station):
        """Test decoding with single target."""
        var = RoutingProbabilities(
            mock_jobclass,
            source=mock_source,
            targets=[mock_station]
        )
        probs = var.decode(np.array([0.5]))
        assert len(probs) == 1
        assert probs[0] == 1.0  # Single target gets all probability

    def test_decode_multiple_targets(self, mock_jobclass, mock_source, mock_station, mock_station2):
        """Test decoding with multiple targets."""
        var = RoutingProbabilities(
            mock_jobclass,
            source=mock_source,
            targets=[mock_station, mock_station2]
        )
        probs = var.decode(np.array([0.5]))
        assert len(probs) == 2
        assert np.isclose(sum(probs), 1.0)  # Probabilities sum to 1
        assert probs[0] == 0.5  # First probability from stick-breaking
        assert probs[1] == 0.5  # Remaining probability

    def test_decode_extreme_values(self, mock_jobclass, mock_source, mock_station, mock_station2):
        """Test decoding at extreme values."""
        var = RoutingProbabilities(
            mock_jobclass,
            source=mock_source,
            targets=[mock_station, mock_station2]
        )
        # All to first target
        probs = var.decode(np.array([1.0]))
        assert np.isclose(probs[0], 1.0)
        assert np.isclose(probs[1], 0.0)

        # All to second target
        probs = var.decode(np.array([0.0]))
        assert np.isclose(probs[0], 0.0)
        assert np.isclose(probs[1], 1.0)

    def test_variable_type(self, mock_jobclass, mock_source, mock_station):
        """Test variable type identifier."""
        var = RoutingProbabilities(
            mock_jobclass,
            source=mock_source,
            targets=[mock_station]
        )
        assert var.getVariableType() == "routing"


class TestClassServiceMapping:
    """Tests for ClassServiceMapping decision variable."""

    def test_init(self, mock_jobclass, mock_station, mock_station2):
        """Test initialization."""
        var = ClassServiceMapping(
            mock_jobclass,
            stations=[mock_station, mock_station2]
        )
        assert var.getName() == "TestClass_mapping"
        assert var.getDimension() == 1

    def test_decode(self, mock_jobclass, mock_station, mock_station2):
        """Test decoding to station index."""
        var = ClassServiceMapping(
            mock_jobclass,
            stations=[mock_station, mock_station2]
        )
        # Low value selects first station
        assert var.decode(np.array([0.0])) == 0
        # High value selects second station
        assert var.decode(np.array([0.99])) == 1
        # Middle value still in valid range
        assert var.decode(np.array([0.5])) == 1

    def test_variable_type(self, mock_jobclass, mock_station):
        """Test variable type identifier."""
        var = ClassServiceMapping(mock_jobclass, stations=[mock_station])
        assert var.getVariableType() == "class_mapping"


class TestServiceRate:
    """Tests for ServiceRate decision variable."""

    def test_init(self, mock_station, mock_jobclass):
        """Test initialization."""
        var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        assert var.getName() == "TestQueue_TestClass_rate"
        assert var.getDimension() == 1

    def test_decode_min(self, mock_station, mock_jobclass):
        """Test decoding at minimum."""
        var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        result = var.decode(np.array([0.0]))
        assert np.isclose(result, 0.5)

    def test_decode_max(self, mock_station, mock_jobclass):
        """Test decoding at maximum."""
        var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        result = var.decode(np.array([1.0]))
        assert np.isclose(result, 5.0)

    def test_decode_middle(self, mock_station, mock_jobclass):
        """Test decoding at middle."""
        var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        result = var.decode(np.array([0.5]))
        # 0.5 + 0.5 * 4.5 = 2.75
        assert np.isclose(result, 2.75)

    def test_variable_type(self, mock_station, mock_jobclass):
        """Test variable type identifier."""
        var = ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0))
        assert var.getVariableType() == "service_rate"


class TestJobPopulation:
    """Tests for JobPopulation decision variable."""

    def test_init(self, mock_jobclass):
        """Test initialization."""
        var = JobPopulation(mock_jobclass, bounds=(1, 100))
        assert var.getName() == "TestClass_population"
        assert var.getDimension() == 1

    def test_decode(self, mock_jobclass):
        """Test decoding to job count."""
        var = JobPopulation(mock_jobclass, bounds=(1, 100))
        assert var.decode(np.array([0.0])) == 1
        assert var.decode(np.array([1.0])) == 100
        # 1 + 0.5 * 99 = 50.5, np.round uses banker's rounding -> 50
        assert var.decode(np.array([0.5])) == 50

    def test_variable_type(self, mock_jobclass):
        """Test variable type identifier."""
        var = JobPopulation(mock_jobclass, bounds=(1, 100))
        assert var.getVariableType() == "job_population"

    def test_apply_targets_model_copy_not_base(self):
        """apply() must set population on the model's class, not the base one."""
        from .conftest import MockJobClass, MockNetwork

        base_class = MockJobClass("Jobs", njobs=1)
        # A distinct class instance with the same name lives in the model copy,
        # mimicking deepcopy performed by the evaluator.
        model_class = MockJobClass("Jobs", njobs=1)
        model = MockNetwork("copy", nodes=[], classes=[model_class])

        var = JobPopulation(base_class, bounds=(1, 100))
        var.apply(model, 42)

        assert model_class.getNumberOfJobs() == 42
        # The base class the optimizer holds must remain untouched.
        assert base_class.getNumberOfJobs() == 1

    def test_apply_missing_class_is_noop(self, mock_jobclass):
        """apply() is a safe no-op when the class is absent from the model."""
        from .conftest import MockNetwork

        model = MockNetwork("empty", nodes=[], classes=[])
        var = JobPopulation(mock_jobclass, bounds=(1, 100))
        var.apply(model, 10)  # must not raise


class TestClassPriority:
    """Tests for ClassPriority decision variable."""

    def test_init_levels_mode(self, mock_jobclass, mock_jobclass2):
        """Test initialization in levels mode."""
        var = ClassPriority([mock_jobclass, mock_jobclass2], mode='levels')
        assert var.getName() == "class_priorities"
        assert var.getDimension() == 2  # One dimension per class
        assert var.getMode() == 'levels'

    def test_init_permutation_mode(self, mock_jobclass, mock_jobclass2):
        """Test initialization in permutation mode."""
        var = ClassPriority([mock_jobclass, mock_jobclass2], mode='permutation')
        assert var.getDimension() == 1  # n-1 dimensions for n classes

    def test_decode_levels(self, mock_jobclass, mock_jobclass2):
        """Test decoding in levels mode."""
        var = ClassPriority(
            [mock_jobclass, mock_jobclass2],
            mode='levels',
            priority_range=(1, 10)
        )
        priorities = var.decode(np.array([0.0, 1.0]))
        assert priorities[0] == 1  # Min priority
        assert priorities[1] == 10  # Max priority

    def test_decode_permutation(self, mock_jobclass, mock_jobclass2):
        """Test decoding in permutation mode."""
        var = ClassPriority(
            [mock_jobclass, mock_jobclass2],
            mode='permutation'
        )
        # Higher random key means higher position in ordering
        order = var.decode(np.array([0.9]))
        assert len(order) == 2
        assert set(order) == {0, 1}  # Contains both indices

    def test_variable_type(self, mock_jobclass):
        """Test variable type identifier."""
        var = ClassPriority([mock_jobclass], mode='levels')
        assert var.getVariableType() == "class_priority"


class TestDecisionVariableProperties:
    """Test property access on decision variables."""

    def test_name_property(self, mock_station):
        """Test name property access."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.name == "TestQueue_servers"

    def test_dimension_property(self, mock_station):
        """Test dimension property access."""
        var = ServerAllocation(mock_station, bounds=(1, 10))
        assert var.dimension == 1
