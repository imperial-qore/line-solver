"""
Hora Scheduler - DAG-based Workflow Scheduling on Heterogeneous Systems

This module provides classes for modeling and solving workflow scheduling
problems on heterogeneous computing systems. Hora Scheduler is part of the
LINE ecosystem and follows LINE's coding conventions and patterns.

Key Classes:
    - Workflow: DAG representation of a workflow (analogous to Network)
    - Task: A computational task node (analogous to Node)
    - Dependency: Edge with communication cost
    - Machine: A computing unit with speed
    - HeterogeneousMachines: Machines with varying speeds and communication costs
    - IdenticalMachines: Machines with identical speeds
    - HEFT: Heterogeneous Earliest Finish Time algorithm
    - ScheduleResult: Scheduling results container

Example:
    >>> from line_solver.scheduling import *
    >>>
    >>> # Create workflow programmatically
    >>> wf = Workflow('MapReduce')
    >>> t1 = Task(wf, 'Map', computation=10.0)
    >>> t2 = Task(wf, 'Reduce', computation=5.0)
    >>> wf.addDependency(t1, t2, comm_cost=3.0)
    >>>
    >>> # Define machines
    >>> machines = HeterogeneousMachines('Cluster')
    >>> machines.addMachine('M1', speed=1.0)
    >>> machines.addMachine('M2', speed=1.5)
    >>>
    >>> # Or use identical machines
    >>> machines = IdenticalMachines('Cluster', num_machines=3)
    >>>
    >>> # Schedule
    >>> solver = HEFT(wf, machines)
    >>> result = solver.schedule()
    >>> print(result.schedule_table)
"""

import numpy as np
import pandas as pd
import networkx as nx
from typing import List, Dict, Optional, Tuple, Union
from collections import namedtuple
import json

# Event namedtuple for schedule entries
ScheduleEvent = namedtuple('ScheduleEvent', ['task', 'start', 'end'])


class Task:
    """
    A task node in a workflow DAG.

    Represents a computational task with a specified workload (computation weight).
    Tasks can have release times for dynamic scheduling scenarios.

    Args:
        workflow (Workflow): Parent workflow model.
        name (str): Name of the task.
        computation (float): Computation weight (time on unit-speed processor).
        release_time (float): Earliest start time (default: 0.0).
        weight (float): Priority weight for WSPT scheduling (default: 1.0).
        due_date (float): Due date for EDD scheduling (default: inf).

    Example:
        >>> wf = Workflow('MyWorkflow')
        >>> task = Task(wf, 'ProcessData', computation=10.0)
        >>> task.setReleaseTime(5.0)  # Can't start before t=5
        >>> task.setWeight(2.0)       # Higher priority
        >>> task.setDueDate(20.0)     # Must finish by t=20
    """

    def __init__(self, workflow: 'Workflow', name: str,
                 computation: float, release_time: float = 0.0,
                 weight: float = 1.0, due_date: float = float('inf')):
        self._workflow = workflow
        self._name = name
        self._computation = computation
        self._release_time = release_time
        self._weight = weight
        self._due_date = due_date
        self._index = None
        workflow.addTask(self)

    def getName(self) -> str:
        """Get task name."""
        return self._name

    def getComputation(self) -> float:
        """Get computation weight."""
        return self._computation

    def setComputation(self, value: float) -> None:
        """Set computation weight."""
        self._computation = value
        # Update graph node attribute
        if self._name in self._workflow._graph:
            self._workflow._graph.nodes[self._name]['C'] = value

    def getReleaseTime(self) -> float:
        """Get release time (for dynamic scheduling)."""
        return self._release_time

    def setReleaseTime(self, value: float) -> None:
        """Set release time."""
        self._release_time = value
        if self._name in self._workflow._graph:
            self._workflow._graph.nodes[self._name]['release_time'] = value

    def getWeight(self) -> float:
        """Get task weight (for WSPT scheduling)."""
        return self._weight

    def setWeight(self, value: float) -> None:
        """Set task weight."""
        self._weight = value
        if self._name in self._workflow._graph:
            self._workflow._graph.nodes[self._name]['weight'] = value

    def getDueDate(self) -> float:
        """Get due date (for EDD scheduling)."""
        return self._due_date

    def setDueDate(self, value: float) -> None:
        """Set due date."""
        self._due_date = value
        if self._name in self._workflow._graph:
            self._workflow._graph.nodes[self._name]['due_date'] = value

    def getPredecessors(self) -> List['Task']:
        """Get all predecessor tasks."""
        pred_names = list(self._workflow._graph.predecessors(self._name))
        return [self._workflow.getTask(name) for name in pred_names]

    def getSuccessors(self) -> List['Task']:
        """Get all successor tasks."""
        succ_names = list(self._workflow._graph.successors(self._name))
        return [self._workflow.getTask(name) for name in succ_names]

    def getIndex(self) -> int:
        """Get task index in workflow."""
        if self._index is None:
            self._index = self._workflow._tasks.index(self)
        return self._index

    def __index__(self) -> int:
        """Support Python indexing."""
        return self.getIndex()

    def __repr__(self) -> str:
        return f"Task('{self._name}', C={self._computation})"

    def __str__(self) -> str:
        return self._name

    # Snake case aliases
    get_name = getName
    name = property(getName)
    get_computation = getComputation
    set_computation = setComputation
    computation = property(getComputation, setComputation)
    get_release_time = getReleaseTime
    set_release_time = setReleaseTime
    release_time = property(getReleaseTime, setReleaseTime)
    get_weight = getWeight
    set_weight = setWeight
    weight = property(getWeight, setWeight)
    get_due_date = getDueDate
    set_due_date = setDueDate
    due_date = property(getDueDate, setDueDate)
    get_predecessors = getPredecessors
    predecessors = property(getPredecessors)
    get_successors = getSuccessors
    successors = property(getSuccessors)
    get_index = getIndex
    index = property(getIndex)


class Dependency:
    """
    A dependency edge between tasks with communication cost.

    Represents the data flow from a source task to a target task,
    with an associated communication cost for inter-processor transfer.

    Args:
        source (Task): Source task (predecessor).
        target (Task): Target task (successor).
        comm_cost (float): Communication cost (data transfer weight).
    """

    def __init__(self, source: 'Task', target: 'Task', comm_cost: float = 0.0):
        self._source = source
        self._target = target
        self._comm_cost = comm_cost

    def getSource(self) -> 'Task':
        """Get source task."""
        return self._source

    def getTarget(self) -> 'Task':
        """Get target task."""
        return self._target

    def getCommCost(self) -> float:
        """Get communication cost."""
        return self._comm_cost

    def setCommCost(self, value: float) -> None:
        """Set communication cost."""
        self._comm_cost = value

    def __repr__(self) -> str:
        return f"Dependency({self._source.getName()} -> {self._target.getName()}, cost={self._comm_cost})"

    # Snake case aliases
    get_source = getSource
    source = property(getSource)
    get_target = getTarget
    target = property(getTarget)
    get_comm_cost = getCommCost
    set_comm_cost = setCommCost
    comm_cost = property(getCommCost, setCommCost)


class Workflow:
    """
    DAG-based workflow model for task scheduling.

    A Workflow represents a directed acyclic graph (DAG) of tasks with
    dependencies and communication costs. Analogous to Network in LINE.

    Args:
        name (str): Name of the workflow model.

    Example:
        >>> wf = Workflow('MapReduceWorkflow')
        >>> t1 = Task(wf, 'Map1', computation=10.0)
        >>> t2 = Task(wf, 'Reduce1', computation=5.0)
        >>> wf.addDependency(t1, t2, comm_cost=3.0)
    """

    def __init__(self, name: str):
        self._name = name
        self._tasks = []
        self._task_map = {}  # name -> Task
        self._dependencies = []
        self._graph = nx.DiGraph()

    @staticmethod
    def fromGML(filepath: str) -> 'Workflow':
        """
        Load workflow from GML file.

        Args:
            filepath: Path to GML file.

        Returns:
            Workflow loaded from file.

        Note:
            Expected node attributes: 'label' (name), 'C' (computation weight)
            Expected edge attributes: 'label' (communication cost)
        """
        graph = nx.read_gml(filepath)

        # Extract workflow name from filename
        import os
        name = os.path.splitext(os.path.basename(filepath))[0]
        wf = Workflow(name)

        # Create tasks from nodes
        task_map = {}
        for node in graph.nodes():
            node_data = graph.nodes[node]
            comp = float(node_data.get('C', 1.0))
            release_time = float(node_data.get('release_time', 0.0))
            # Use label as name if available, otherwise use node id
            task_name = str(node_data.get('label', node))
            task = Task(wf, task_name, computation=comp, release_time=release_time)
            task_map[node] = task

        # Create dependencies from edges
        for src, tgt in graph.edges():
            edge_data = graph.edges[src, tgt]
            comm_cost = float(edge_data.get('label', 0.0))
            src_task = task_map[src]
            tgt_task = task_map[tgt]
            wf.addDependency(src_task, tgt_task, comm_cost=comm_cost)

        return wf

    @staticmethod
    def from_gml(filepath: str) -> 'Workflow':
        """Snake case alias for fromGML."""
        return Workflow.fromGML(filepath)

    def addTask(self, task: 'Task') -> None:
        """Register a task with this workflow (called automatically by Task.__init__)."""
        if task._name in self._task_map:
            raise ValueError(f"Task with name '{task._name}' already exists")
        self._tasks.append(task)
        self._task_map[task._name] = task
        task._index = len(self._tasks) - 1
        # Add node to graph
        self._graph.add_node(task._name, C=task._computation,
                            release_time=task._release_time)

    def getTasks(self) -> List['Task']:
        """Get all tasks in this workflow."""
        return list(self._tasks)

    def getTask(self, name: str) -> 'Task':
        """Get a task by name."""
        return self._task_map.get(name)

    def addDependency(self, source: 'Task', target: 'Task',
                      comm_cost: float = 0.0) -> 'Dependency':
        """
        Add a dependency edge between tasks.

        Args:
            source: Source task (predecessor).
            target: Target task (successor).
            comm_cost: Communication cost (data transfer weight).

        Returns:
            The created Dependency object.
        """
        dep = Dependency(source, target, comm_cost)
        self._dependencies.append(dep)
        self._graph.add_edge(source._name, target._name, label=comm_cost)
        return dep

    def getDependencies(self) -> List['Dependency']:
        """Get all dependencies."""
        return list(self._dependencies)

    def getGraph(self) -> nx.DiGraph:
        """Get the underlying NetworkX graph."""
        return self._graph

    def validate(self) -> bool:
        """
        Validate the workflow is a valid DAG.

        Returns:
            True if valid DAG, raises exception otherwise.
        """
        if not nx.is_directed_acyclic_graph(self._graph):
            raise ValueError("Workflow contains cycles - not a valid DAG")
        return True

    def getEntryTasks(self) -> List['Task']:
        """Get tasks with no predecessors."""
        entry_names = [n for n in self._graph.nodes()
                      if self._graph.in_degree(n) == 0]
        return [self._task_map[name] for name in entry_names]

    def getExitTasks(self) -> List['Task']:
        """Get tasks with no successors."""
        exit_names = [n for n in self._graph.nodes()
                     if self._graph.out_degree(n) == 0]
        return [self._task_map[name] for name in exit_names]

    def getCriticalPath(self) -> Tuple[List['Task'], float]:
        """
        Compute critical path through the DAG.

        Returns:
            Tuple of (list of tasks on critical path, critical path length).
        """
        if len(self._tasks) == 0:
            return [], 0.0

        # Compute longest path using dynamic programming
        topo_order = list(nx.topological_sort(self._graph))
        dist = {node: 0.0 for node in topo_order}
        pred = {node: None for node in topo_order}

        for node in topo_order:
            task = self._task_map[node]
            node_dist = dist[node] + task._computation
            for succ in self._graph.successors(node):
                edge_cost = self._graph.edges[node, succ].get('label', 0.0)
                new_dist = node_dist + edge_cost
                if new_dist > dist[succ]:
                    dist[succ] = new_dist
                    pred[succ] = node

        # Find exit node with maximum distance
        exit_nodes = self.getExitTasks()
        max_exit = max(exit_nodes, key=lambda t: dist[t._name] + t._computation)
        cp_length = dist[max_exit._name] + max_exit._computation

        # Reconstruct path
        path = []
        current = max_exit._name
        while current is not None:
            path.append(self._task_map[current])
            current = pred[current]
        path.reverse()

        return path, cp_length

    def getNumberOfTasks(self) -> int:
        """Get number of tasks."""
        return len(self._tasks)

    def getName(self) -> str:
        """Get workflow name."""
        return self._name

    def plot(self, **kwargs) -> None:
        """Plot the workflow DAG."""
        try:
            import matplotlib.pyplot as plt
            pos = nx.spring_layout(self._graph)
            nx.draw(self._graph, pos, with_labels=True,
                   node_color='lightblue', node_size=500,
                   font_size=8, font_weight='bold', **kwargs)
            edge_labels = nx.get_edge_attributes(self._graph, 'label')
            nx.draw_networkx_edge_labels(self._graph, pos, edge_labels)
            plt.title(f"Workflow: {self._name}")
            plt.show()
        except ImportError:
            print("matplotlib required for plotting")

    def __repr__(self) -> str:
        return f"Workflow('{self._name}', tasks={len(self._tasks)}, deps={len(self._dependencies)})"

    # Snake case aliases
    add_task = addTask
    get_tasks = getTasks
    get_task = getTask
    add_dependency = addDependency
    get_dependencies = getDependencies
    get_graph = getGraph
    get_entry_tasks = getEntryTasks
    get_exit_tasks = getExitTasks
    get_critical_path = getCriticalPath
    get_number_of_tasks = getNumberOfTasks
    number_of_tasks = property(getNumberOfTasks)
    get_name = getName
    name = property(getName)


class Machine:
    """
    A machine (computing unit) in a scheduling system.

    Represents a processing unit with a relative computational speed.
    Speed=1.0 is the reference speed; higher values indicate faster machines.

    Args:
        name (str): Machine name/identifier.
        speed (float): Relative computational speed (default: 1.0).

    Example:
        >>> fast = Machine('GPU', speed=2.0)
        >>> slow = Machine('CPU', speed=1.0)
    """

    def __init__(self, name: str, speed: float = 1.0):
        self._name = name
        self._speed = speed
        self._index = None

    def getName(self) -> str:
        """Get machine name."""
        return self._name

    def getSpeed(self) -> float:
        """Get machine speed."""
        return self._speed

    def setSpeed(self, value: float) -> None:
        """Set machine speed."""
        self._speed = value

    def getIndex(self) -> int:
        """Get machine index."""
        return self._index

    def __repr__(self) -> str:
        return f"Machine('{self._name}', speed={self._speed})"

    def __str__(self) -> str:
        return self._name

    # Snake case aliases
    get_name = getName
    name = property(getName)
    get_speed = getSpeed
    set_speed = setSpeed
    speed = property(getSpeed, setSpeed)
    get_index = getIndex
    index = property(getIndex)


# Alias for backward compatibility
Processor = Machine


class HeterogeneousMachines:
    """
    A heterogeneous computing system with multiple machines of varying speeds.

    Models a system with machines of varying speeds and
    inter-machine communication costs. Communication costs
    depend on both the data size and the link speed between machines.

    Args:
        name (str): System name.

    Example:
        >>> machines = HeterogeneousMachines('CloudCluster')
        >>> m1 = machines.addMachine('GPU1', speed=2.0)
        >>> m2 = machines.addMachine('GPU2', speed=2.0)
        >>> m3 = machines.addMachine('CPU', speed=1.0)
        >>> machines.setCommSpeed(m1, m2, 10.0)  # Fast GPU-GPU link
        >>> machines.setCommSpeed(m1, m3, 1.0)   # Slower GPU-CPU link
    """

    def __init__(self, name: str):
        self._name = name
        self._machines = []
        self._machine_map = {}  # name -> Machine
        self._comm_speeds = None  # Will be ndarray

    def addMachine(self, name: str, speed: float = 1.0) -> 'Machine':
        """
        Add a machine to the system.

        Args:
            name: Machine name.
            speed: Relative computational speed.

        Returns:
            The created Machine object.
        """
        machine = Machine(name, speed)
        machine._index = len(self._machines)
        self._machines.append(machine)
        self._machine_map[name] = machine

        # Resize communication matrix
        n = len(self._machines)
        new_comm = np.full((n, n), 1.0)
        np.fill_diagonal(new_comm, np.inf)  # Local comm is instant
        if self._comm_speeds is not None and n > 1:
            old_n = n - 1
            new_comm[:old_n, :old_n] = self._comm_speeds
        self._comm_speeds = new_comm

        return machine

    def getMachines(self) -> List['Machine']:
        """Get all machines."""
        return list(self._machines)

    def getMachine(self, name_or_index: Union[str, int]) -> 'Machine':
        """Get machine by name or index."""
        if isinstance(name_or_index, int):
            return self._machines[name_or_index]
        return self._machine_map.get(name_or_index)

    def getNumberOfMachines(self) -> int:
        """Get number of machines."""
        return len(self._machines)

    def setCommSpeed(self, mach1: Union['Machine', int],
                     mach2: Union['Machine', int], speed: float) -> None:
        """
        Set communication speed between two machines.

        Args:
            mach1: First machine (or index).
            mach2: Second machine (or index).
            speed: Communication speed (higher = faster).
        """
        idx1 = mach1 if isinstance(mach1, int) else mach1._index
        idx2 = mach2 if isinstance(mach2, int) else mach2._index
        self._comm_speeds[idx1, idx2] = speed
        self._comm_speeds[idx2, idx1] = speed  # Symmetric

    def setCommSpeedMatrix(self, matrix: np.ndarray) -> None:
        """
        Set entire communication speed matrix.

        Args:
            matrix: Square matrix of communication speeds.
                   np.inf on diagonal indicates instant local communication.
        """
        n = len(self._machines)
        if matrix.shape != (n, n):
            raise ValueError(f"Matrix shape {matrix.shape} doesn't match {n} machines")
        self._comm_speeds = matrix.copy()

    def getCommSpeed(self, mach1: Union['Machine', int],
                     mach2: Union['Machine', int]) -> float:
        """Get communication speed between two machines."""
        idx1 = mach1 if isinstance(mach1, int) else mach1._index
        idx2 = mach2 if isinstance(mach2, int) else mach2._index
        return self._comm_speeds[idx1, idx2]

    def getCommSpeedMatrix(self) -> np.ndarray:
        """Get communication speed matrix."""
        return self._comm_speeds.copy()

    def getMachineSpeeds(self) -> np.ndarray:
        """Get array of machine speeds."""
        return np.array([m._speed for m in self._machines])

    @staticmethod
    def fromArrays(name: str, speeds: np.ndarray,
                   comm_matrix: np.ndarray) -> 'HeterogeneousMachines':
        """
        Create system from numpy arrays (compatibility with existing code).

        Args:
            name: System name.
            speeds: Array of machine speeds.
            comm_matrix: Communication speed matrix.

        Returns:
            A HeterogeneousMachines configured from arrays.
        """
        machines = HeterogeneousMachines(name)
        for i, speed in enumerate(speeds):
            machines.addMachine(f'M{i}', speed=float(speed))
        machines.setCommSpeedMatrix(comm_matrix)
        return machines

    def getName(self) -> str:
        """Get system name."""
        return self._name

    # Internal methods for solver compatibility (use _machines internally)
    def getProcessors(self) -> List['Machine']:
        """Get all machines (alias for solver compatibility)."""
        return self.getMachines()

    def getNumberOfProcessors(self) -> int:
        """Get number of machines (alias for solver compatibility)."""
        return self.getNumberOfMachines()

    @property
    def _processors(self):
        """Internal access for solver compatibility."""
        return self._machines

    def __repr__(self) -> str:
        return f"HeterogeneousMachines('{self._name}', machines={len(self._machines)})"

    # Snake case aliases
    add_machine = addMachine
    get_machines = getMachines
    get_machine = getMachine
    get_number_of_machines = getNumberOfMachines
    number_of_machines = property(getNumberOfMachines)
    set_comm_speed = setCommSpeed
    set_comm_speed_matrix = setCommSpeedMatrix
    get_comm_speed = getCommSpeed
    get_comm_speed_matrix = getCommSpeedMatrix
    get_machine_speeds = getMachineSpeeds
    from_arrays = fromArrays
    get_name = getName
    name = property(getName)
    # Backward compatibility aliases
    addProcessor = addMachine
    getProcessors = getMachines
    getProcessor = getMachine
    getProcessorSpeeds = getMachineSpeeds
    add_processor = addMachine
    get_processors = getMachines
    get_processor = getMachine
    get_processor_speeds = getMachineSpeeds


# Alias for backward compatibility
HeterogeneousSystem = HeterogeneousMachines


class IdenticalMachines(HeterogeneousMachines):
    """
    A system with multiple identical machines (same speed).

    All machines have the same computational speed and uniform
    communication costs between them.

    Args:
        name (str): System name.
        num_machines (int): Number of machines.
        speed (float): Speed for all machines (default: 1.0).
        comm_speed (float): Communication speed between all pairs (default: 1.0).

    Example:
        >>> machines = IdenticalMachines('Cluster', num_machines=4)
        >>> machines = IdenticalMachines('FastCluster', num_machines=8, speed=2.0)
    """

    def __init__(self, name: str, num_machines: int,
                 speed: float = 1.0, comm_speed: float = 1.0):
        super().__init__(name)
        self._uniform_speed = speed
        self._uniform_comm_speed = comm_speed

        for i in range(num_machines):
            self.addMachine(f'M{i}', speed=speed)

        # Set uniform communication speeds
        for i in range(num_machines):
            for j in range(i + 1, num_machines):
                self.setCommSpeed(i, j, comm_speed)

    def __repr__(self) -> str:
        return f"IdenticalMachines('{self._name}', machines={len(self._machines)}, speed={self._uniform_speed})"


class ScheduleResult:
    """
    Container for workflow scheduling results.

    Provides access to schedule details, performance metrics,
    and visualization capabilities.

    Attributes:
        schedule_table (DataFrame): Schedule with columns:
            - Task: Task name
            - Processor: Assigned processor
            - StartTime: Scheduled start time
            - FinishTime: Scheduled finish time
        makespan (float): Total schedule completion time.
        slr (float): Schedule Length Ratio (makespan / critical_path_length).
        speedup (float): Speedup over sequential execution.
        efficiency (float): Scheduling efficiency (speedup / num_processors).
    """

    def __init__(self):
        self._schedule_data = {}  # task_name -> {processor, start, end}
        self._processor_schedules = {}  # processor_name -> [ScheduleEvent]
        self._makespan = 0.0
        self._critical_path_length = 0.0
        self._total_computation = 0.0
        self._num_processors = 0

    def _addTaskSchedule(self, task_name: str, processor: str,
                         start: float, end: float) -> None:
        """Add a task to the schedule (internal use)."""
        self._schedule_data[task_name] = {
            'processor': processor,
            'start_time': start,
            'finish_time': end
        }
        if processor not in self._processor_schedules:
            self._processor_schedules[processor] = []
        self._processor_schedules[processor].append(
            ScheduleEvent(task_name, start, end)
        )
        self._makespan = max(self._makespan, end)

    def getScheduleTable(self) -> pd.DataFrame:
        """Get schedule as pandas DataFrame."""
        rows = []
        for task_name, data in self._schedule_data.items():
            rows.append({
                'Task': task_name,
                'Processor': data['processor'],
                'StartTime': data['start_time'],
                'FinishTime': data['finish_time']
            })
        df = pd.DataFrame(rows)
        if len(df) > 0:
            df = df.sort_values('StartTime').reset_index(drop=True)
        return df

    def getMakespan(self) -> float:
        """Get schedule makespan (completion time)."""
        return self._makespan

    def getSLR(self) -> float:
        """Get Schedule Length Ratio (makespan / critical_path_length)."""
        if self._critical_path_length > 0:
            return self._makespan / self._critical_path_length
        return 1.0

    def getSpeedup(self) -> float:
        """Get speedup over sequential execution."""
        if self._makespan > 0:
            return self._total_computation / self._makespan
        return 1.0

    def getEfficiency(self) -> float:
        """Get scheduling efficiency (speedup / num_processors)."""
        if self._num_processors > 0:
            return self.getSpeedup() / self._num_processors
        return 1.0

    def getProcessorSchedule(self, processor: Union[str, int]) -> List[ScheduleEvent]:
        """Get schedule for a specific processor."""
        if isinstance(processor, int):
            processor = f'P{processor}'
        return self._processor_schedules.get(processor, [])

    def getTaskAssignment(self, task: Union[str, 'Task']) -> Optional[dict]:
        """Get assignment details for a specific task."""
        task_name = task if isinstance(task, str) else task.getName()
        return self._schedule_data.get(task_name)

    def toJSON(self, filepath: str = None) -> str:
        """
        Export schedule to JSON format.

        Args:
            filepath: Optional file path to write JSON.

        Returns:
            JSON string representation.
        """
        json_str = json.dumps(self._schedule_data, indent=4)
        if filepath:
            with open(filepath, 'w') as f:
                f.write(json_str)
        return json_str

    def toDict(self) -> dict:
        """Convert to dictionary format."""
        return {
            'schedule': self._schedule_data,
            'makespan': self._makespan,
            'slr': self.getSLR(),
            'speedup': self.getSpeedup(),
            'efficiency': self.getEfficiency()
        }

    def plotGantt(self, figsize: Tuple[int, int] = (12, 6), **kwargs) -> None:
        """
        Generate Gantt chart of the schedule.

        Args:
            figsize: Figure size (width, height).
            **kwargs: Additional matplotlib arguments.
        """
        try:
            import matplotlib.pyplot as plt
            import matplotlib.patches as mpatches

            fig, ax = plt.subplots(figsize=figsize)

            # Get unique processors and assign y-positions
            processors = sorted(self._processor_schedules.keys())
            proc_y = {p: i for i, p in enumerate(processors)}

            # Color map for tasks
            colors = plt.cm.Set3(np.linspace(0, 1, len(self._schedule_data)))
            task_colors = {task: colors[i] for i, task in enumerate(self._schedule_data.keys())}

            # Draw bars
            for task_name, data in self._schedule_data.items():
                proc = data['processor']
                start = data['start_time']
                duration = data['finish_time'] - start
                y = proc_y[proc]

                ax.barh(y, duration, left=start, height=0.6,
                       color=task_colors[task_name], edgecolor='black', linewidth=0.5)

                # Add task label
                ax.text(start + duration / 2, y, task_name,
                       ha='center', va='center', fontsize=8, fontweight='bold')

            # Customize plot
            ax.set_yticks(range(len(processors)))
            ax.set_yticklabels(processors)
            ax.set_xlabel('Time')
            ax.set_ylabel('Processor')
            ax.set_title(f'Schedule Gantt Chart (Makespan: {self._makespan:.2f})')
            ax.set_xlim(0, self._makespan * 1.05)
            ax.grid(axis='x', alpha=0.3)

            plt.tight_layout()
            plt.show()

        except ImportError:
            print("matplotlib required for Gantt chart")

    def __repr__(self) -> str:
        return f"ScheduleResult(tasks={len(self._schedule_data)}, makespan={self._makespan:.2f})"

    # Snake case aliases
    schedule_table = property(getScheduleTable)
    get_schedule_table = getScheduleTable
    makespan = property(getMakespan)
    get_makespan = getMakespan
    slr = property(getSLR)
    get_slr = getSLR
    speedup = property(getSpeedup)
    get_speedup = getSpeedup
    efficiency = property(getEfficiency)
    get_efficiency = getEfficiency
    get_processor_schedule = getProcessorSchedule
    get_task_assignment = getTaskAssignment
    to_json = toJSON
    to_dict = toDict
    plot_gantt = plotGantt


class WorkflowSolver:
    """
    Base class for workflow scheduling algorithms.

    Provides common functionality for all scheduling algorithms
    including option management and result handling.

    Args:
        workflow (Workflow): The workflow to schedule.
        system (HeterogeneousSystem): The target computing system.
    """

    @staticmethod
    def defaultOptions() -> dict:
        """Get default solver options."""
        return {
            'verbose': False,
            'validate': True,
        }

    def __init__(self, workflow: 'Workflow', system: 'HeterogeneousSystem',
                 *args, **kwargs):
        self._workflow = workflow
        self._system = system
        self._options = self.defaultOptions()
        self._schedule_result = None
        self._processOptions(kwargs)

        # Validate if requested
        if self._options['validate']:
            workflow.validate()

    def _processOptions(self, kwargs):
        """Process solver options from kwargs."""
        for key, value in kwargs.items():
            if key in self._options:
                self._options[key] = value

    def solve(self) -> 'ScheduleResult':
        """Execute the scheduling algorithm (alias for schedule)."""
        return self.schedule()

    def schedule(self) -> 'ScheduleResult':
        """Execute the scheduling algorithm. Override in subclasses."""
        raise NotImplementedError("Subclasses must implement schedule()")

    def getResult(self) -> 'ScheduleResult':
        """Get the last computed result."""
        return self._schedule_result

    def supports(self, workflow: 'Workflow') -> bool:
        """Check if this solver supports the given workflow."""
        return True

    # Snake case aliases
    default_options = defaultOptions
    get_result = getResult
    _process_options = _processOptions


def _findFirstGap(schedule_list: List[ScheduleEvent],
                  earliest_start: float, exec_time: float) -> float:
    """
    Find the earliest available time gap in the processor's task schedule.

    Args:
        schedule_list: List of scheduled events on processor.
        earliest_start: Earliest possible start time.
        exec_time: Required execution time.

    Returns:
        Start time for the task.
    """
    if not schedule_list:
        return earliest_start

    # Sort by start time
    sorted_schedule = sorted(schedule_list, key=lambda x: x.start)

    # Check gap before first task
    if earliest_start + exec_time <= sorted_schedule[0].start:
        return earliest_start

    # Check gaps between scheduled tasks
    for i in range(len(sorted_schedule) - 1):
        gap_start = max(sorted_schedule[i].end, earliest_start)
        gap_end = sorted_schedule[i + 1].start

        if (gap_end - gap_start) >= exec_time:
            return gap_start

    # If no gap is found, schedule after the last task
    return max(sorted_schedule[-1].end, earliest_start)


class HEFT(WorkflowSolver):
    """
    Heterogeneous Earliest Finish Time (HEFT) scheduling algorithm.

    HEFT is a list scheduling heuristic that prioritizes tasks by their
    upward rank and assigns each task to the processor that minimizes
    its earliest finish time.

    Reference:
        Topcuoglu, H., Hariri, S., & Wu, M. Y. (2002). Performance-effective
        and low-complexity task scheduling for heterogeneous computing.
        IEEE Transactions on Parallel and Distributed Systems, 13(3), 260-274.

    Args:
        workflow (Workflow): The workflow to schedule.
        system (HeterogeneousSystem): The target computing system.
        **kwargs: Additional options:
            - verbose (bool): Print progress information (default: False).
            - insertion (bool): Enable insertion-based scheduling (default: True).

    Example:
        >>> wf = Workflow.from_gml('workflow.gml')
        >>> system = HeterogeneousSystem.homogeneous('Cluster', 3)
        >>> solver = HEFT(wf, system)
        >>> result = solver.schedule()
        >>> print(result.schedule_table)
        >>> result.plot_gantt()
    """

    @staticmethod
    def defaultOptions() -> dict:
        options = WorkflowSolver.defaultOptions()
        options.update({
            'insertion': True,  # Enable insertion-based scheduling
        })
        return options

    def __init__(self, workflow: 'Workflow', system: 'HeterogeneousSystem',
                 *args, **kwargs):
        super().__init__(workflow, system, *args, **kwargs)
        self._exec_time_matrix = None
        self._task_list = None
        self._up_rank = {}

    def _computeExecutionTimeMatrix(self) -> np.ndarray:
        """Compute execution time matrix (tasks x processors)."""
        tasks = self._workflow.getTasks()
        processors = self._system.getProcessors()
        n_tasks = len(tasks)
        n_procs = len(processors)

        matrix = np.zeros((n_tasks, n_procs))
        for i, task in enumerate(tasks):
            for j, proc in enumerate(processors):
                matrix[i, j] = task._computation / proc._speed

        return matrix

    def _getCommunicationCost(self, src_task: Task, tgt_task: Task,
                               src_proc: int, tgt_proc: int) -> float:
        """
        Get communication cost between tasks on different processors.

        Args:
            src_task: Source task.
            tgt_task: Target task.
            src_proc: Source processor index.
            tgt_proc: Target processor index.

        Returns:
            Communication cost (0 if same processor).
        """
        if src_proc == tgt_proc:
            return 0.0

        # Get edge weight (communication data size)
        graph = self._workflow._graph
        edge_data = graph.edges.get((src_task._name, tgt_task._name), {})
        comm_data = edge_data.get('label', 0.0)

        # Get communication speed
        comm_speed = self._system.getCommSpeed(src_proc, tgt_proc)

        if comm_speed == np.inf:
            return 0.0  # Instant local communication
        elif comm_speed == 0:
            return float('inf')  # No communication possible

        return comm_data / comm_speed

    def _getAverageCommunicationCost(self, src_task: Task, tgt_task: Task) -> float:
        """Get average communication cost across all processor pairs."""
        graph = self._workflow._graph
        edge_data = graph.edges.get((src_task._name, tgt_task._name), {})
        comm_data = edge_data.get('label', 0.0)

        if comm_data == 0:
            return 0.0

        comm_matrix = self._system._comm_speeds
        n_procs = len(self._system._processors)

        total_cost = 0.0
        count = 0
        for k in range(n_procs):
            for m in range(n_procs):
                if k != m and comm_matrix[k, m] != np.inf and comm_matrix[k, m] > 0:
                    total_cost += comm_data / comm_matrix[k, m]
                    count += 1

        if count > 0:
            return total_cost / count
        return 0.0

    def _computeUpwardRank(self, task: Task) -> float:
        """
        Compute upward rank for a task recursively.

        Upward rank = mean execution time + max(comm_cost + upward_rank of successors)
        """
        if task._name in self._up_rank:
            return self._up_rank[task._name]

        task_idx = task.getIndex()
        mean_exec_time = np.mean(self._exec_time_matrix[task_idx])

        successors = task.getSuccessors()
        if not successors:
            self._up_rank[task._name] = mean_exec_time
        else:
            max_succ_rank = max(
                self._getAverageCommunicationCost(task, succ) + self._computeUpwardRank(succ)
                for succ in successors
            )
            self._up_rank[task._name] = mean_exec_time + max_succ_rank

        return self._up_rank[task._name]

    def schedule(self) -> 'ScheduleResult':
        """Execute HEFT scheduling algorithm."""
        # Initialize
        self._exec_time_matrix = self._computeExecutionTimeMatrix()
        self._task_list = self._workflow.getTasks()
        self._up_rank = {}

        n_procs = self._system.getNumberOfProcessors()
        processors = self._system.getProcessors()

        # Compute upward ranks for all tasks
        for task in self._task_list:
            self._computeUpwardRank(task)

        # Sort tasks by upward rank (descending)
        sorted_tasks = sorted(self._task_list,
                             key=lambda t: self._up_rank[t._name],
                             reverse=True)

        # Initialize result
        result = ScheduleResult()
        result._num_processors = n_procs
        result._total_computation = sum(t._computation for t in self._task_list)
        _, result._critical_path_length = self._workflow.getCriticalPath()

        # Track task assignments and processor schedules
        task_assignment = {}  # task_name -> processor_index
        task_finish_time = {}  # task_name -> finish_time
        processor_schedules = {p.getName(): [] for p in processors}

        if self._options['verbose']:
            print(f"HEFT: Scheduling {len(sorted_tasks)} tasks on {n_procs} processors")

        # Schedule each task
        for task in sorted_tasks:
            task_idx = task.getIndex()
            release_time = task.getReleaseTime()

            best_processor = None
            best_start_time = float('inf')
            best_finish_time = float('inf')

            # Try each processor
            for proc_idx, proc in enumerate(processors):
                exec_time = self._exec_time_matrix[task_idx, proc_idx]

                # Compute earliest start time (EST)
                est = release_time
                for pred in task.getPredecessors():
                    pred_proc = task_assignment.get(pred._name)
                    if pred_proc is not None:
                        pred_finish = task_finish_time[pred._name]
                        comm_cost = self._getCommunicationCost(
                            pred, task, pred_proc, proc_idx
                        )
                        est = max(est, pred_finish + comm_cost)

                # Find available slot (with insertion if enabled)
                if self._options['insertion']:
                    available_start = _findFirstGap(
                        processor_schedules[proc.getName()], est, exec_time
                    )
                else:
                    # No insertion - schedule after last task
                    if processor_schedules[proc.getName()]:
                        last_end = max(e.end for e in processor_schedules[proc.getName()])
                        available_start = max(est, last_end)
                    else:
                        available_start = est

                finish_time = available_start + exec_time

                # Select processor with minimum EFT
                if finish_time < best_finish_time:
                    best_processor = proc_idx
                    best_start_time = available_start
                    best_finish_time = finish_time

            # Record schedule
            proc_name = processors[best_processor].getName()
            processor_schedules[proc_name].append(
                ScheduleEvent(task._name, best_start_time, best_finish_time)
            )
            task_assignment[task._name] = best_processor
            task_finish_time[task._name] = best_finish_time

            result._addTaskSchedule(task._name, proc_name,
                                   best_start_time, best_finish_time)

            if self._options['verbose']:
                print(f"  {task._name} -> {proc_name} [{best_start_time:.2f}, {best_finish_time:.2f}]")

        self._schedule_result = result

        if self._options['verbose']:
            print(f"HEFT: Makespan = {result.getMakespan():.2f}, SLR = {result.getSLR():.2f}")

        return result


class ListScheduler(WorkflowSolver):
    """
    Base class for list scheduling algorithms (SPT, SEPT, LEPT, WSPT, EDD).

    List schedulers sort tasks by a priority rule and then schedule them
    in order, assigning each task to the machine that minimizes finish time.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system (HeterogeneousMachines or IdenticalMachines).
        **kwargs: Additional options:
            - verbose (bool): Print progress information (default: False).
    """

    def __init__(self, workflow: 'Workflow', machines: 'HeterogeneousMachines',
                 *args, **kwargs):
        super().__init__(workflow, machines, *args, **kwargs)

    def _getPriority(self, task: Task) -> float:
        """
        Get priority value for a task. Override in subclasses.
        Lower values = higher priority (scheduled first).
        """
        raise NotImplementedError("Subclasses must implement _getPriority")

    def _getSortReverse(self) -> bool:
        """Whether to sort in descending order. Default False (ascending)."""
        return False

    def schedule(self) -> 'ScheduleResult':
        """Execute list scheduling algorithm."""
        tasks = self._workflow.getTasks()
        n_machines = self._system.getNumberOfMachines()
        machines = self._system.getMachines()

        # Sort tasks by priority
        sorted_tasks = sorted(tasks, key=lambda t: self._getPriority(t),
                             reverse=self._getSortReverse())

        # Initialize result
        result = ScheduleResult()
        result._num_processors = n_machines
        result._total_computation = sum(t._computation for t in tasks)
        if len(tasks) > 0:
            _, result._critical_path_length = self._workflow.getCriticalPath()

        # Track machine availability (finish time of last task on each machine)
        machine_available = {m.getName(): 0.0 for m in machines}
        machine_schedules = {m.getName(): [] for m in machines}

        if self._options['verbose']:
            print(f"{self.__class__.__name__}: Scheduling {len(sorted_tasks)} tasks on {n_machines} machines")

        for task in sorted_tasks:
            release_time = task.getReleaseTime()

            # Find best machine (earliest finish time)
            best_machine = None
            best_start = float('inf')
            best_finish = float('inf')

            for machine in machines:
                # Compute execution time on this machine
                exec_time = task._computation / machine._speed

                # Earliest start = max(release_time, machine_available)
                start_time = max(release_time, machine_available[machine.getName()])
                finish_time = start_time + exec_time

                if finish_time < best_finish:
                    best_machine = machine
                    best_start = start_time
                    best_finish = finish_time

            # Schedule task on best machine
            machine_available[best_machine.getName()] = best_finish
            machine_schedules[best_machine.getName()].append(
                ScheduleEvent(task._name, best_start, best_finish)
            )
            result._addTaskSchedule(task._name, best_machine.getName(),
                                   best_start, best_finish)

            if self._options['verbose']:
                print(f"  {task._name} -> {best_machine.getName()} [{best_start:.2f}, {best_finish:.2f}]")

        self._schedule_result = result

        if self._options['verbose']:
            print(f"{self.__class__.__name__}: Makespan = {result.getMakespan():.2f}")

        return result


class SPT(ListScheduler):
    """
    Shortest Processing Time (SPT) scheduling rule.

    Schedules tasks in order of increasing processing time.
    Optimal for minimizing mean flow time on a single machine.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Short', computation=2.0)
        >>> Task(wf, 'Long', computation=10.0)
        >>> machines = IdenticalMachines('M', num_machines=1)
        >>> result = SPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = processing time (lower = higher priority)."""
        return task._computation


class LPT(ListScheduler):
    """
    Longest Processing Time (LPT) scheduling rule.

    Schedules tasks in order of decreasing processing time.
    Good heuristic for minimizing makespan on parallel machines.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Short', computation=2.0)
        >>> Task(wf, 'Long', computation=10.0)
        >>> machines = IdenticalMachines('M', num_machines=2)
        >>> result = LPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = processing time (higher = higher priority)."""
        return task._computation

    def _getSortReverse(self) -> bool:
        """Sort descending (longest first)."""
        return True


class SEPT(ListScheduler):
    """
    Shortest Expected Processing Time (SEPT) scheduling rule.

    For stochastic scheduling, schedules tasks by expected processing time.
    Uses the computation attribute as the expected processing time.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Fast', computation=2.0)   # E[p] = 2.0
        >>> Task(wf, 'Slow', computation=10.0)  # E[p] = 10.0
        >>> result = SEPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = expected processing time (lower = higher priority)."""
        return task._computation


class LEPT(ListScheduler):
    """
    Longest Expected Processing Time (LEPT) scheduling rule.

    For stochastic scheduling, schedules tasks by expected processing time
    in decreasing order. Good for makespan minimization on parallel machines.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Fast', computation=2.0)   # E[p] = 2.0
        >>> Task(wf, 'Slow', computation=10.0)  # E[p] = 10.0
        >>> result = LEPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = expected processing time (higher = higher priority)."""
        return task._computation

    def _getSortReverse(self) -> bool:
        """Sort descending (longest expected first)."""
        return True


class WSPT(ListScheduler):
    """
    Weighted Shortest Processing Time (WSPT) scheduling rule.

    Schedules tasks in order of decreasing weight/processing_time ratio.
    Optimal for minimizing weighted sum of completion times on a single machine.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'HighPrio', computation=5.0, weight=10.0)  # ratio = 2.0
        >>> Task(wf, 'LowPrio', computation=5.0, weight=1.0)    # ratio = 0.2
        >>> result = WSPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = weight / processing_time (higher = higher priority)."""
        if task._computation > 0:
            return task._weight / task._computation
        return float('inf')  # Zero processing time = highest priority

    def _getSortReverse(self) -> bool:
        """Sort descending (highest ratio first)."""
        return True


class EDD(ListScheduler):
    """
    Earliest Due Date (EDD) scheduling rule.

    Schedules tasks in order of increasing due date.
    Optimal for minimizing maximum lateness on a single machine.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system.

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Urgent', computation=5.0, due_date=10.0)
        >>> Task(wf, 'Later', computation=5.0, due_date=20.0)
        >>> result = EDD(wf, machines).schedule()
        >>> # Check lateness
        >>> for _, row in result.schedule_table.iterrows():
        ...     task = wf.getTask(row['Task'])
        ...     lateness = row['FinishTime'] - task.due_date
        ...     print(f"{task.name}: Lateness = {lateness}")
    """

    def _getPriority(self, task: Task) -> float:
        """Priority = due date (earlier = higher priority)."""
        return task._due_date


class PreemptiveScheduler(WorkflowSolver):
    """
    Base class for preemptive scheduling algorithms (SRPT, WSRPT).

    Preemptive schedulers can interrupt running tasks when a higher-priority
    task becomes available (e.g., when a new task is released).

    Note: For single-machine scheduling only (uses first machine if multiple provided).
    """

    def __init__(self, workflow: 'Workflow', machines: 'HeterogeneousMachines',
                 *args, **kwargs):
        super().__init__(workflow, machines, *args, **kwargs)

    def _getPriority(self, task: Task, remaining: float) -> float:
        """
        Get priority value for a task based on remaining processing time.
        Override in subclasses. Higher values = higher priority (scheduled first).
        """
        raise NotImplementedError("Subclasses must implement _getPriority")

    def schedule(self) -> 'ScheduleResult':
        """Execute preemptive scheduling algorithm."""
        tasks = self._workflow.getTasks()
        machines = self._system.getMachines()
        machine = machines[0]  # Use first machine for preemptive scheduling

        # Initialize result
        result = ScheduleResult()
        result._num_processors = 1
        result._total_computation = sum(t._computation for t in tasks)
        if len(tasks) > 0:
            _, result._critical_path_length = self._workflow.getCriticalPath()

        if len(tasks) == 0:
            self._schedule_result = result
            return result

        # Track remaining processing time for each task
        remaining = {t._name: t._computation / machine._speed for t in tasks}
        release_times = {t._name: t._release_time for t in tasks}

        # Collect all event times (release times)
        event_times = sorted(set(release_times.values()))

        # Track task execution segments
        task_segments = {t._name: [] for t in tasks}  # [(start, end), ...]
        completed = set()

        current_time = 0.0
        current_task = None

        if self._options['verbose']:
            print(f"{self.__class__.__name__}: Scheduling {len(tasks)} tasks (preemptive)")

        while len(completed) < len(tasks):
            # Find available tasks (released and not completed)
            available = [t for t in tasks
                        if release_times[t._name] <= current_time
                        and t._name not in completed
                        and remaining[t._name] > 0]

            if not available:
                # Jump to next release time
                future_releases = [r for r in release_times.values() if r > current_time]
                if future_releases:
                    current_time = min(future_releases)
                    continue
                else:
                    break

            # Select task with highest priority
            best_task = max(available,
                           key=lambda t: self._getPriority(t, remaining[t._name]))

            # Find next event (release time or task completion)
            next_release = float('inf')
            for t in tasks:
                if release_times[t._name] > current_time and t._name not in completed:
                    next_release = min(next_release, release_times[t._name])

            task_finish = current_time + remaining[best_task._name]
            next_event = min(next_release, task_finish)

            # Execute task until next event
            exec_duration = next_event - current_time
            remaining[best_task._name] -= exec_duration

            # Record segment
            task_segments[best_task._name].append((current_time, next_event))

            if self._options['verbose']:
                print(f"  t={current_time:.2f}: {best_task._name} runs until {next_event:.2f} (remaining={remaining[best_task._name]:.2f})")

            # Check if task completed
            if remaining[best_task._name] <= 1e-9:  # Floating point tolerance
                completed.add(best_task._name)
                remaining[best_task._name] = 0

            current_time = next_event

        # Merge consecutive segments and record results
        for task in tasks:
            segments = task_segments[task._name]
            if segments:
                # Merge consecutive segments
                merged = []
                for seg in sorted(segments):
                    if merged and abs(seg[0] - merged[-1][1]) < 1e-9:
                        merged[-1] = (merged[-1][0], seg[1])
                    else:
                        merged.append(seg)

                # Use first start and last end for the schedule table
                start_time = merged[0][0]
                finish_time = merged[-1][1]
                result._addTaskSchedule(task._name, machine.getName(),
                                       start_time, finish_time)

        self._schedule_result = result

        if self._options['verbose']:
            print(f"{self.__class__.__name__}: Makespan = {result.getMakespan():.2f}")

        return result


class SRPT(PreemptiveScheduler):
    """
    Shortest Remaining Processing Time (SRPT) scheduling rule.

    Preemptive rule that always executes the task with the shortest
    remaining processing time. Optimal for minimizing mean flow time.

    When a new task arrives with shorter remaining time than the current
    task, the current task is preempted.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system (uses first machine).

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'Long', computation=10.0, release_time=0.0)
        >>> Task(wf, 'Short', computation=2.0, release_time=1.0)  # Arrives at t=1
        >>> result = SRPT(wf, machines).schedule()
        >>> # Short preempts Long at t=1
    """

    def _getPriority(self, task: Task, remaining: float) -> float:
        """Priority = -remaining (shorter remaining = higher priority)."""
        return -remaining  # Negative so shorter = higher


class WSRPT(PreemptiveScheduler):
    """
    Weighted Shortest Remaining Processing Time (WSRPT) scheduling rule.

    Preemptive rule that always executes the task with the highest
    weight/remaining_processing_time ratio.

    Args:
        workflow (Workflow): The workflow to schedule.
        machines: The target computing system (uses first machine).

    Example:
        >>> wf = Workflow('Jobs')
        >>> Task(wf, 'LowPrio', computation=10.0, weight=1.0, release_time=0.0)
        >>> Task(wf, 'HighPrio', computation=5.0, weight=10.0, release_time=2.0)
        >>> result = WSRPT(wf, machines).schedule()
    """

    def _getPriority(self, task: Task, remaining: float) -> float:
        """Priority = weight / remaining (higher ratio = higher priority)."""
        if remaining > 0:
            return task._weight / remaining
        return float('inf')


# Public API
__all__ = [
    'Workflow',
    'Task',
    'Dependency',
    'Machine',
    'Processor',  # Alias for Machine
    'HeterogeneousMachines',
    'HeterogeneousSystem',  # Alias for HeterogeneousMachines
    'IdenticalMachines',
    'WorkflowSolver',
    'ListScheduler',
    'PreemptiveScheduler',
    'HEFT',
    # Non-preemptive rules
    'SPT',
    'LPT',
    'SEPT',
    'LEPT',
    'WSPT',
    'EDD',
    # Preemptive rules
    'SRPT',
    'WSRPT',
    'ScheduleResult',
    'ScheduleEvent',
]
