"""
Native Python implementation of LayeredNetworkGenerator.

Generates random layered queueing network (LQN) models without requiring the
Java backend. Port of the MATLAB @LayeredNetworkGenerator class.
"""

import math
import random
from typing import List, Sequence

from ..constants import GlobalConstants, SchedStrategy
from ..distributions import Exp, Immediate
from ..layered import Activity, Entry, LayeredNetwork, Processor, Task


def sample_integer_value(rng: Sequence[float]) -> int:
    """
    Sample an integer value uniformly from a given range.

    Args:
        rng: Two-element (lower, upper) range. The bounds are rounded inwards.

    Returns:
        An integer in [ceil(rng[0]), floor(rng[1])].
    """
    return random.randint(int(math.ceil(rng[0])), int(math.floor(rng[1])))


def sample_real_value(rng: Sequence[float]) -> float:
    """
    Sample a real value uniformly from a given range.

    Args:
        rng: Two-element (lower, upper) range.

    Returns:
        A float in [rng[0], rng[1]].
    """
    return rng[0] + (rng[1] - rng[0]) * random.random()


def choose_boolean_value(probability: float) -> bool:
    """
    Choose a boolean value for a given probability.

    Args:
        probability: Probability of returning True.

    Returns:
        True with the given probability.
    """
    return random.random() < probability


def make_integer_vector(length: int, total: int) -> List[int]:
    """
    Make a vector of positive integers with given length and sum.

    Each element starts at one and the surplus total-length is distributed
    uniformly at random over the elements.

    Args:
        length: Number of elements.
        total: Required sum of the elements.

    Returns:
        List of `length` positive integers summing to `total`.
    """
    vector = [1] * length
    for _ in range(total - length):
        vector[random.randint(0, length - 1)] += 1
    return vector


def _as_distribution(mean: float):
    """Coerce a scalar mean into a distribution, as MATLAB setHostDemand does."""
    if mean <= GlobalConstants.FineTol:
        return Immediate()
    return Exp(1.0 / mean)


class LayeredNetworkGenerator:
    """
    A generator for creating random layered queueing network models.

    Characteristics of the generated models are configured via the
    generator's properties. Mirrors the MATLAB LayeredNetworkGenerator
    and the Java jline.gen.LayeredNetworkGenerator.

    Attributes:
        population_range: (min, max) range of reference-task populations.
        think_time_range: (min, max) range of client think times.
        task_inf_probability: Probability that a task is infinite-server.
        proc_inf_probability: Probability that a processor is infinite-server.
        task_multi_range: (min, max) range of task multiplicities.
        proc_multi_range: (min, max) range of processor multiplicities.
        host_demand_range: (min, max) range of activity host demands.
        synch_call_range: (min, max) range of synchronous call means.

    Examples:
        >>> gen = LayeredNetworkGenerator(population_range=(5, 45),
        ...                               think_time_range=(10, 100))
        >>> model = gen.generate(2, 2, 4, 2)
    """

    def __init__(
        self,
        population_range: Sequence[float] = (1, 1),
        think_time_range: Sequence[float] = (1, 1),
        task_inf_probability: float = 0.0,
        proc_inf_probability: float = 0.0,
        task_multi_range: Sequence[float] = (1, 1),
        proc_multi_range: Sequence[float] = (1, 1),
        host_demand_range: Sequence[float] = (1, 1),
        synch_call_range: Sequence[float] = (1, 1)
    ):
        """
        Initialize a LayeredNetworkGenerator with configurable properties.

        Args:
            population_range: (min, max) reference-task population range.
            think_time_range: (min, max) client think-time range.
            task_inf_probability: Probability that a task is infinite-server.
            proc_inf_probability: Probability that a processor is infinite-server.
            task_multi_range: (min, max) task multiplicity range.
            proc_multi_range: (min, max) processor multiplicity range.
            host_demand_range: (min, max) activity host-demand range.
            synch_call_range: (min, max) synchronous call mean range.
        """
        self.population_range = population_range
        self.think_time_range = think_time_range
        self.task_inf_probability = task_inf_probability
        self.proc_inf_probability = proc_inf_probability
        self.task_multi_range = task_multi_range
        self.proc_multi_range = proc_multi_range
        self.host_demand_range = host_demand_range
        self.synch_call_range = synch_call_range

        # Populated by generate
        self.c_activities = []
        self.c_entries = []
        self.c_tasks = []
        self.c_processors = []
        self.activities = []
        self.entries = []
        self.tasks = []
        self.processors = []
        self.num_tasks_per_level = []
        self.num_tasks_per_processor = []

    @staticmethod
    def _check_range(rng, name, strictly_positive):
        """Validate a two-element range, returning it as a tuple."""
        if len(rng) != 2:
            raise ValueError(f"{name} must have two elements")
        lower_ok = rng[0] > 0 if strictly_positive else rng[0] >= 0
        if not (lower_ok and rng[0] <= rng[1]):
            raise ValueError(f"{name} is not valid")
        return (rng[0], rng[1])

    @staticmethod
    def _check_probability(value, name):
        """Validate a probability value."""
        if not (0 <= value <= 1):
            raise ValueError(f"{name} is not valid")
        return value

    @property
    def population_range(self):
        """Get the reference-task population range."""
        return self._population_range

    @population_range.setter
    def population_range(self, rng):
        self._population_range = self._check_range(rng, 'Population range', True)

    @property
    def think_time_range(self):
        """Get the client think-time range."""
        return self._think_time_range

    @think_time_range.setter
    def think_time_range(self, rng):
        self._think_time_range = self._check_range(rng, 'Think time range', False)

    @property
    def task_inf_probability(self):
        """Get the probability that a task is infinite-server."""
        return self._task_inf_probability

    @task_inf_probability.setter
    def task_inf_probability(self, value):
        self._task_inf_probability = self._check_probability(
            value, 'Task infinite probability')

    @property
    def proc_inf_probability(self):
        """Get the probability that a processor is infinite-server."""
        return self._proc_inf_probability

    @proc_inf_probability.setter
    def proc_inf_probability(self, value):
        self._proc_inf_probability = self._check_probability(
            value, 'Processor infinite probability')

    @property
    def task_multi_range(self):
        """Get the task multiplicity range."""
        return self._task_multi_range

    @task_multi_range.setter
    def task_multi_range(self, rng):
        self._task_multi_range = self._check_range(rng, 'Task multiplicity range', True)

    @property
    def proc_multi_range(self):
        """Get the processor multiplicity range."""
        return self._proc_multi_range

    @proc_multi_range.setter
    def proc_multi_range(self, rng):
        self._proc_multi_range = self._check_range(
            rng, 'Processor multiplicity range', True)

    @property
    def host_demand_range(self):
        """Get the activity host-demand range."""
        return self._host_demand_range

    @host_demand_range.setter
    def host_demand_range(self, rng):
        self._host_demand_range = self._check_range(rng, 'Host demand range', False)

    @property
    def synch_call_range(self):
        """Get the synchronous call mean range."""
        return self._synch_call_range

    @synch_call_range.setter
    def synch_call_range(self, rng):
        self._synch_call_range = self._check_range(rng, 'Synchronous call range', True)

    def generate(self, num_clients: int, num_levels: int, num_tasks: int,
                 num_processors: int) -> LayeredNetwork:
        """
        Generate a random layered queueing network model.

        Args:
            num_clients: Number of client (reference) tasks.
            num_levels: Number of task layers.
            num_tasks: Number of non-client tasks, distributed over the levels.
            num_processors: Number of processors hosting the non-client tasks.

        Returns:
            A fully constructed LayeredNetwork object.

        Raises:
            ValueError: If arguments are invalid.

        Examples:
            >>> gen = LayeredNetworkGenerator()
            >>> model = gen.generate(1, 2, 2, 2)
        """
        self._validate_args(num_clients, num_levels, num_tasks, num_processors)
        model = LayeredNetwork('lnw')
        self._create_clients(model, num_clients)
        self._create_tasks(model, num_tasks)
        self._create_processors(model, num_processors)
        self._assign_tasks(num_levels, num_tasks, num_processors)
        self._connect_clients_to_tasks(num_clients)
        self._connect_tasks_to_tasks(num_levels)
        self._connect_tasks_to_processors(num_processors)
        return model

    @staticmethod
    def _validate_args(num_clients, num_levels, num_tasks, num_processors):
        """Validate that parameter values for the network are sound."""
        if num_clients < 1:
            raise ValueError("Number of clients is less than one")
        if num_levels < 1:
            raise ValueError("Number of levels is less than one")
        if num_tasks < 1:
            raise ValueError("Number of tasks is less than one")
        if num_processors < 1:
            raise ValueError("Number of processors is less than one")
        if num_levels > num_tasks:
            raise ValueError("Number of levels is greater than that of tasks")
        if num_processors > num_tasks:
            raise ValueError("Number of processors is greater than that of tasks")

    def _create_clients(self, model, num_clients):
        """Create the clients in the layered network."""
        self.c_activities = []
        self.c_entries = []
        self.c_tasks = []
        self.c_processors = []
        for c in range(num_clients):
            population = sample_integer_value(self.population_range)
            think_time = sample_real_value(self.think_time_range)

            activity = Activity(model, f'c_activity_{c + 1}', _as_distribution(think_time))
            entry = Entry(model, f'c_entry_{c + 1}')
            task = Task(model, f'c_task_{c + 1}', population, SchedStrategy.REF)
            processor = Processor(model, f'c_processor_{c + 1}',
                                  float('inf'), SchedStrategy.INF)

            activity.on(task).bound_to(entry)
            entry.on(task)
            task.on(processor)

            self.c_activities.append(activity)
            self.c_entries.append(entry)
            self.c_tasks.append(task)
            self.c_processors.append(processor)

    def _create_tasks(self, model, num_tasks):
        """Create the tasks in the layered network."""
        self.activities = []
        self.entries = []
        self.tasks = []
        for t in range(num_tasks):
            if choose_boolean_value(self.task_inf_probability):
                multiplicity = float('inf')
                scheduling = SchedStrategy.INF
            else:
                multiplicity = sample_integer_value(self.task_multi_range)
                scheduling = SchedStrategy.FCFS
            host_demand = sample_real_value(self.host_demand_range)

            activity = Activity(model, f'activity_{t + 1}', _as_distribution(host_demand))
            entry = Entry(model, f'entry_{t + 1}')
            task = Task(model, f'task_{t + 1}', multiplicity, scheduling)

            activity.on(task).bound_to(entry).replies_to(entry)
            entry.on(task)

            self.activities.append(activity)
            self.entries.append(entry)
            self.tasks.append(task)

    def _create_processors(self, model, num_processors):
        """Create the processors in the layered network."""
        self.processors = []
        for i in range(num_processors):
            if choose_boolean_value(self.proc_inf_probability):
                multiplicity = float('inf')
                scheduling = SchedStrategy.INF
            else:
                multiplicity = sample_integer_value(self.proc_multi_range)
                scheduling = SchedStrategy.PS
            self.processors.append(
                Processor(model, f'processor_{i + 1}', multiplicity, scheduling))

    def _assign_tasks(self, num_levels, num_tasks, num_processors):
        """Assign the tasks to different levels and processors."""
        self.num_tasks_per_level = make_integer_vector(num_levels, num_tasks)
        self.num_tasks_per_processor = make_integer_vector(num_processors, num_tasks)

    def _connect_clients_to_tasks(self, num_clients):
        """Connect the clients to the first-level tasks."""
        client_connected = [False] * num_clients
        for t in range(self.num_tasks_per_level[0]):
            synch_call = sample_real_value(self.synch_call_range)

            c = sample_integer_value([1, num_clients]) - 1
            self.c_activities[c].synch_call(self.entries[t], synch_call)
            client_connected[c] = True

        for c in range(num_clients):
            if client_connected[c]:
                continue

            synch_call = sample_real_value(self.synch_call_range)

            t = sample_integer_value([1, self.num_tasks_per_level[0]]) - 1
            self.c_activities[c].synch_call(self.entries[t], synch_call)
            client_connected[c] = True

    def _connect_tasks_to_tasks(self, num_levels):
        """Connect the tasks between adjacent levels."""
        num_tasks = self.num_tasks_per_level[0]
        for l in range(1, num_levels):
            for t2 in range(num_tasks, num_tasks + self.num_tasks_per_level[l]):
                synch_call = sample_real_value(self.synch_call_range)

                # Caller drawn among the tasks of the previous level
                t1 = sample_integer_value(
                    [num_tasks - self.num_tasks_per_level[l - 1] + 1, num_tasks]) - 1
                self.activities[t1].synch_call(self.entries[t2], synch_call)
            num_tasks += self.num_tasks_per_level[l]

    def _connect_tasks_to_processors(self, num_processors):
        """Connect the tasks to the processors."""
        num_tasks = 0
        for p in range(num_processors):
            for t in range(num_tasks, num_tasks + self.num_tasks_per_processor[p]):
                self.tasks[t].on(self.processors[p])
            num_tasks += self.num_tasks_per_processor[p]

    # camelCase aliases for compatibility with the MATLAB and Java spellings
    @property
    def populationRange(self): return self.population_range

    @populationRange.setter
    def populationRange(self, v): self.population_range = v

    @property
    def thinkTimeRange(self): return self.think_time_range

    @thinkTimeRange.setter
    def thinkTimeRange(self, v): self.think_time_range = v

    @property
    def taskInfProbability(self): return self.task_inf_probability

    @taskInfProbability.setter
    def taskInfProbability(self, v): self.task_inf_probability = v

    @property
    def procInfProbability(self): return self.proc_inf_probability

    @procInfProbability.setter
    def procInfProbability(self, v): self.proc_inf_probability = v

    @property
    def taskMultiRange(self): return self.task_multi_range

    @taskMultiRange.setter
    def taskMultiRange(self, v): self.task_multi_range = v

    @property
    def procMultiRange(self): return self.proc_multi_range

    @procMultiRange.setter
    def procMultiRange(self, v): self.proc_multi_range = v

    @property
    def hostDemandRange(self): return self.host_demand_range

    @hostDemandRange.setter
    def hostDemandRange(self, v): self.host_demand_range = v

    @property
    def synchCallRange(self): return self.synch_call_range

    @synchCallRange.setter
    def synchCallRange(self, v): self.synch_call_range = v

    sampleIntegerValue = staticmethod(sample_integer_value)
    sampleRealValue = staticmethod(sample_real_value)
    chooseBooleanValue = staticmethod(choose_boolean_value)
    makeIntegerVector = staticmethod(make_integer_vector)
