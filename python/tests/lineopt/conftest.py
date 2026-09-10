"""
Pytest configuration and fixtures for line_solver.opt tests.
"""

import pytest
import numpy as np
from unittest.mock import MagicMock, PropertyMock


class MockStation:
    """Mock LINE Station for testing."""
    def __init__(self, name="TestQueue", num_servers=1):
        self._name = name
        self._num_servers = num_servers

    def getName(self):
        return self._name

    def getNumberOfServers(self):
        return self._num_servers

    def setNumberOfServers(self, n):
        self._num_servers = n


class MockJobClass:
    """Mock LINE JobClass for testing."""
    def __init__(self, name="TestClass", njobs=1):
        self._name = name
        self._priority = 1
        self._njobs = njobs

    def getName(self):
        return self._name

    def setPriority(self, p):
        self._priority = p

    def getNumberOfJobs(self):
        return self._njobs

    def setNumberOfJobs(self, n):
        self._njobs = n


class MockNode:
    """Mock LINE Node for testing."""
    def __init__(self, name="Node"):
        self._name = name

    def getName(self):
        return self._name


class MockNetwork:
    """Mock LINE Network for testing."""
    def __init__(self, name="TestNetwork", nodes=None, classes=None):
        self._name = name
        self._nodes = nodes or []
        self._classes = classes or []

    def getName(self):
        return self._name

    def getNodes(self):
        return self._nodes

    def getClasses(self):
        return self._classes

    def get_class_by_name(self, name):
        for c in self._classes:
            if c.getName() == name:
                return c
        return None

    # Note: the real LINE Network exposes station/class indices via snake_case
    # get_node_index / get_class_index only (no CamelCase alias). Metric
    # extraction no longer relies on positional indices, so none are mocked.


@pytest.fixture
def mock_station():
    """Create a mock LINE Station."""
    return MockStation("TestQueue", 1)


@pytest.fixture
def mock_station2():
    """Create a second mock LINE Station."""
    return MockStation("TestQueue2", 2)


@pytest.fixture
def mock_jobclass():
    """Create a mock LINE JobClass."""
    return MockJobClass("TestClass")


@pytest.fixture
def mock_jobclass2():
    """Create a second mock LINE JobClass."""
    return MockJobClass("TestClass2")


@pytest.fixture
def mock_source():
    """Create a mock LINE Source node."""
    return MockNode("Source")


@pytest.fixture
def mock_network(mock_station, mock_jobclass):
    """Create a mock LINE Network."""
    return MockNetwork("TestNetwork", [mock_station], [mock_jobclass])


@pytest.fixture
def sample_evaluation_result():
    """Create a sample EvaluationResult for testing."""
    from line_solver.opt.results import EvaluationResult

    result = EvaluationResult()
    result.feasible = True
    result.response_times = {
        ("TestQueue", "TestClass"): 0.5,
        ("TestQueue", "TestClass2"): 0.8,
        ("TestQueue2", "TestClass"): 1.0,
    }
    result.throughputs = {
        ("TestQueue", "TestClass"): 10.0,
        ("TestQueue", "TestClass2"): 5.0,
        ("TestQueue2", "TestClass"): 8.0,
    }
    result.utilizations = {
        "TestQueue": 0.6,
        "TestQueue2": 0.8,
    }
    result.queue_lengths = {
        ("TestQueue", "TestClass"): 2.5,
        ("TestQueue", "TestClass2"): 1.5,
    }
    result.solve_time = 0.1
    result.solver_used = "MockSolver"
    return result


@pytest.fixture
def simple_problem(mock_network, mock_station):
    """Create a simple optimization problem for testing."""
    from line_solver.opt.problem import OptimizationProblem
    from line_solver.opt.variables import ServerAllocation
    from line_solver.opt.objectives import MinimizeCost

    problem = OptimizationProblem(mock_network)
    problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
    problem.setObjective(MinimizeCost(server_cost={mock_station: 100.0}))
    return problem


@pytest.fixture
def multi_var_problem(mock_network, mock_station, mock_station2, mock_jobclass):
    """Create a problem with multiple variable types."""
    from line_solver.opt.problem import OptimizationProblem
    from line_solver.opt.variables import ServerAllocation, ServiceRate
    from line_solver.opt.objectives import MinimizeCost

    # Update network to include both stations
    mock_network._nodes = [mock_station, mock_station2]

    problem = OptimizationProblem(mock_network)
    problem.addVariable(ServerAllocation(mock_station, bounds=(1, 10)))
    problem.addVariable(ServerAllocation(mock_station2, bounds=(1, 5)))
    problem.addVariable(ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0)))
    problem.setObjective(MinimizeCost(
        server_cost={mock_station: 100.0, mock_station2: 50.0}
    ))
    return problem


@pytest.fixture
def variables(mock_station, mock_jobclass):
    """Create test variables."""
    from line_solver.opt.variables import ServerAllocation, ServiceRate
    return [
        ServerAllocation(mock_station, bounds=(1, 10)),
        ServiceRate(mock_station, mock_jobclass, bounds=(0.5, 5.0)),
    ]
