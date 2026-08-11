"""
Tests for the Hora Scheduler module.

This module tests the workflow scheduling functionality including:
- Workflow and Task creation
- Machine configuration (Identical and Heterogeneous)
- HEFT algorithm
- GML file loading
- Result generation
"""

import unittest
import numpy as np
import os
import sys

# Add parent directory to path for direct execution
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from scheduling import (
    Workflow, Task, Dependency,
    Machine, HeterogeneousMachines, IdenticalMachines,
    HEFT, ScheduleResult
)


class TestTask(unittest.TestCase):
    """Test Task class."""

    def test_task_creation(self):
        wf = Workflow('Test')
        task = Task(wf, 'T1', computation=10.0)
        self.assertEqual(task.getName(), 'T1')
        self.assertEqual(task.getComputation(), 10.0)
        self.assertEqual(task.getReleaseTime(), 0.0)

    def test_task_with_release_time(self):
        wf = Workflow('Test')
        task = Task(wf, 'T1', computation=10.0, release_time=5.0)
        self.assertEqual(task.getReleaseTime(), 5.0)

    def test_task_snake_case_aliases(self):
        wf = Workflow('Test')
        task = Task(wf, 'T1', computation=10.0)
        self.assertEqual(task.get_name(), 'T1')
        self.assertEqual(task.get_computation(), 10.0)
        self.assertEqual(task.name, 'T1')


class TestWorkflow(unittest.TestCase):
    """Test Workflow class."""

    def test_workflow_creation(self):
        wf = Workflow('TestWorkflow')
        self.assertEqual(wf.getName(), 'TestWorkflow')
        self.assertEqual(wf.getNumberOfTasks(), 0)

    def test_add_tasks(self):
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        self.assertEqual(wf.getNumberOfTasks(), 2)
        self.assertEqual(len(wf.getTasks()), 2)

    def test_add_dependency(self):
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        dep = wf.addDependency(t1, t2, comm_cost=3.0)

        self.assertEqual(len(wf.getDependencies()), 1)
        self.assertEqual(dep.getCommCost(), 3.0)
        self.assertIn(t2, t1.getSuccessors())
        self.assertIn(t1, t2.getPredecessors())

    def test_entry_exit_tasks(self):
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        t3 = Task(wf, 'T3', computation=8.0)
        wf.addDependency(t1, t2)
        wf.addDependency(t2, t3)

        entry = wf.getEntryTasks()
        exit = wf.getExitTasks()

        self.assertEqual(len(entry), 1)
        self.assertEqual(entry[0].getName(), 'T1')
        self.assertEqual(len(exit), 1)
        self.assertEqual(exit[0].getName(), 'T3')

    def test_critical_path(self):
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        t3 = Task(wf, 'T3', computation=8.0)
        wf.addDependency(t1, t2, comm_cost=2.0)
        wf.addDependency(t2, t3, comm_cost=1.0)

        path, length = wf.getCriticalPath()
        self.assertEqual(len(path), 3)
        self.assertEqual(length, 10.0 + 2.0 + 5.0 + 1.0 + 8.0)  # 26.0

    def test_validate_dag(self):
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        wf.addDependency(t1, t2)
        self.assertTrue(wf.validate())


class TestMachines(unittest.TestCase):
    """Test Machine classes."""

    def test_machine_creation(self):
        m = Machine('GPU', speed=2.0)
        self.assertEqual(m.getName(), 'GPU')
        self.assertEqual(m.getSpeed(), 2.0)

    def test_identical_machines(self):
        machines = IdenticalMachines('Cluster', num_machines=4, speed=1.5)
        self.assertEqual(machines.getNumberOfMachines(), 4)
        speeds = machines.getMachineSpeeds()
        self.assertTrue(np.all(speeds == 1.5))

    def test_heterogeneous_machines(self):
        machines = HeterogeneousMachines('HetCluster')
        machines.addMachine('GPU', speed=2.0)
        machines.addMachine('CPU', speed=1.0)

        self.assertEqual(machines.getNumberOfMachines(), 2)
        self.assertEqual(machines.getMachine('GPU').getSpeed(), 2.0)
        self.assertEqual(machines.getMachine(1).getSpeed(), 1.0)

    def test_comm_speed(self):
        machines = HeterogeneousMachines('Test')
        machines.addMachine('M0', speed=1.0)
        machines.addMachine('M1', speed=1.0)
        machines.setCommSpeed(0, 1, 5.0)

        self.assertEqual(machines.getCommSpeed(0, 1), 5.0)
        self.assertEqual(machines.getCommSpeed(1, 0), 5.0)  # Symmetric

    def test_from_arrays(self):
        speeds = np.array([0.5, 1.0, 2.0])
        comm = np.array([
            [np.inf, 1.0, 2.0],
            [1.0, np.inf, 1.5],
            [2.0, 1.5, np.inf]
        ])
        machines = HeterogeneousMachines.from_arrays('System', speeds, comm)

        self.assertEqual(machines.getNumberOfMachines(), 3)
        self.assertTrue(np.allclose(machines.getMachineSpeeds(), speeds))


class TestHEFT(unittest.TestCase):
    """Test HEFT scheduling algorithm."""

    def test_simple_chain(self):
        """Test scheduling a simple chain of tasks."""
        wf = Workflow('Chain')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=5.0)
        t3 = Task(wf, 'T3', computation=8.0)
        wf.addDependency(t1, t2)
        wf.addDependency(t2, t3)

        machines = IdenticalMachines('Cluster', num_machines=2)
        solver = HEFT(wf, machines)
        result = solver.schedule()

        self.assertIsInstance(result, ScheduleResult)
        self.assertEqual(result.getMakespan(), 23.0)  # 10 + 5 + 8 on single machine

    def test_parallel_tasks(self):
        """Test scheduling tasks that can run in parallel."""
        wf = Workflow('Parallel')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=10.0)
        t3 = Task(wf, 'T3', computation=10.0)
        t4 = Task(wf, 'T4', computation=5.0)

        wf.addDependency(t1, t4)
        wf.addDependency(t2, t4)
        wf.addDependency(t3, t4)

        machines = IdenticalMachines('Cluster', num_machines=3)
        result = HEFT(wf, machines).schedule()

        # T1, T2, T3 should run in parallel, then T4
        self.assertEqual(result.getMakespan(), 15.0)  # 10 + 5

    def test_heterogeneous_scheduling(self):
        """Test scheduling on heterogeneous machines."""
        wf = Workflow('Het')
        t1 = Task(wf, 'T1', computation=10.0)

        machines = HeterogeneousMachines('System')
        machines.addMachine('Fast', speed=2.0)
        machines.addMachine('Slow', speed=1.0)

        result = HEFT(wf, machines).schedule()

        # Task should be assigned to faster machine
        self.assertEqual(result.getMakespan(), 5.0)  # 10 / 2.0

    def test_with_release_times(self):
        """Test dynamic scheduling with release times."""
        wf = Workflow('Dynamic')
        t1 = Task(wf, 'T1', computation=5.0, release_time=0.0)
        t2 = Task(wf, 'T2', computation=5.0, release_time=10.0)

        machines = IdenticalMachines('Cluster', num_machines=1)
        result = HEFT(wf, machines).schedule()

        # T2 can't start before release time
        table = result.getScheduleTable()
        t2_row = table[table['Task'] == 'T2'].iloc[0]
        self.assertGreaterEqual(t2_row['StartTime'], 10.0)

    def test_schedule_result_metrics(self):
        """Test schedule result metrics."""
        wf = Workflow('Test')
        t1 = Task(wf, 'T1', computation=10.0)
        t2 = Task(wf, 'T2', computation=10.0)
        wf.addDependency(t1, t2)

        machines = IdenticalMachines('Cluster', num_machines=2)
        result = HEFT(wf, machines).schedule()

        self.assertGreater(result.getMakespan(), 0)
        self.assertGreater(result.getSLR(), 0)
        self.assertGreater(result.getSpeedup(), 0)
        self.assertGreater(result.getEfficiency(), 0)

    def test_schedule_table_format(self):
        """Test schedule table has correct format."""
        wf = Workflow('Test')
        Task(wf, 'T1', computation=10.0)

        machines = IdenticalMachines('Cluster', num_machines=1)
        result = HEFT(wf, machines).schedule()

        table = result.getScheduleTable()
        self.assertIn('Task', table.columns)
        self.assertIn('Processor', table.columns)
        self.assertIn('StartTime', table.columns)
        self.assertIn('FinishTime', table.columns)


class TestGMLImport(unittest.TestCase):
    """Test GML file import."""

    def test_load_gml(self):
        """Test loading workflow from GML file."""
        gml_path = '/home/gcasale/Dropbox/code/iso25-wflowsched-wenqi.git/static/dag/10/Tau_0.gml'
        if os.path.exists(gml_path):
            wf = Workflow.from_gml(gml_path)
            self.assertEqual(wf.getNumberOfTasks(), 10)
            self.assertTrue(wf.validate())
        else:
            self.skipTest("GML file not found")

    def test_schedule_gml_workflow(self):
        """Test scheduling a workflow loaded from GML."""
        gml_path = '/home/gcasale/Dropbox/code/iso25-wflowsched-wenqi.git/static/dag/10/Tau_0.gml'
        if os.path.exists(gml_path):
            wf = Workflow.from_gml(gml_path)

            speeds = np.array([0.5, 1.0, 0.25])
            comm = np.array([
                [np.inf, 1, 1.5],
                [1, np.inf, 2],
                [1.5, 2, np.inf]
            ])
            machines = HeterogeneousMachines.from_arrays('System', speeds, comm)

            result = HEFT(wf, machines).schedule()
            self.assertGreater(result.getMakespan(), 0)
        else:
            self.skipTest("GML file not found")


class TestJSONExport(unittest.TestCase):
    """Test JSON export functionality."""

    def test_to_json(self):
        """Test exporting schedule to JSON."""
        wf = Workflow('Test')
        Task(wf, 'T1', computation=10.0)

        machines = IdenticalMachines('Cluster', num_machines=1)
        result = HEFT(wf, machines).schedule()

        json_str = result.toJSON()
        self.assertIn('T1', json_str)
        self.assertIn('processor', json_str)
        self.assertIn('start_time', json_str)


if __name__ == '__main__':
    # Run tests with verbosity
    unittest.main(verbosity=2)
