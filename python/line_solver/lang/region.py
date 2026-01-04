"""
Finite Capacity Region for LINE native Python implementation.

This module provides the Region class for defining finite capacity regions
in queueing networks. Ported from MATLAB implementation in matlab/src/lang/Region.m
and Java implementation in jar/src/main/kotlin/jline/lang/Region.java
"""

from typing import Dict, List, Optional, Union
from .base import Node, JobClass, DropStrategy


class Region:
    """
    A finite capacity region that constrains jobs within a group of nodes.

    Finite capacity regions define capacity constraints on a subset of nodes
    in the network. When capacity is exceeded, jobs can be dropped, blocked,
    or queued according to the configured drop strategy.

    Attributes:
        UNBOUNDED: Constant indicating no capacity limit (-1)
    """

    UNBOUNDED = -1

    def __init__(self, nodes: List[Node], classes: List[JobClass]):
        """
        Initialize a finite capacity region.

        Args:
            nodes: List of nodes included in this region
            classes: List of job classes in the network
        """
        self._name: str = ""
        self._nodes: List[Node] = list(nodes)
        self._classes: List[JobClass] = list(classes)

        # Global limits
        self._global_max_jobs: int = Region.UNBOUNDED
        self._global_max_memory: int = Region.UNBOUNDED

        # Per-class properties
        self._class_max_jobs: Dict[JobClass, int] = {}
        self._class_max_memory: Dict[JobClass, int] = {}
        self._drop_rule: Dict[JobClass, DropStrategy] = {}
        self._class_size: Dict[JobClass, int] = {}
        self._class_weight: Dict[JobClass, float] = {}

        # Initialize defaults for each class
        for job_class in self._classes:
            self._class_max_jobs[job_class] = Region.UNBOUNDED
            self._class_max_memory[job_class] = Region.UNBOUNDED
            self._drop_rule[job_class] = DropStrategy.WAITQ
            self._class_size[job_class] = 1
            self._class_weight[job_class] = 1.0

    @property
    def name(self) -> str:
        """Get the region name."""
        return self._name

    @name.setter
    def name(self, value: str) -> None:
        """Set the region name."""
        self._name = value

    @property
    def nodes(self) -> List[Node]:
        """Get the nodes in this region."""
        return self._nodes

    @property
    def classes(self) -> List[JobClass]:
        """Get the job classes."""
        return self._classes

    # Global limits
    @property
    def global_max_jobs(self) -> int:
        """Get the global maximum jobs limit."""
        return self._global_max_jobs

    def set_global_max_jobs(self, njobs: int) -> 'Region':
        """
        Set the global maximum number of jobs in the region.

        Args:
            njobs: Maximum number of jobs (use Region.UNBOUNDED for no limit)

        Returns:
            Self for method chaining
        """
        self._global_max_jobs = njobs
        return self

    # MATLAB-style alias
    def setGlobalMaxJobs(self, njobs: int) -> 'Region':
        """MATLAB-compatible alias for set_global_max_jobs()."""
        return self.set_global_max_jobs(njobs)

    @property
    def global_max_memory(self) -> int:
        """Get the global maximum memory limit."""
        return self._global_max_memory

    def set_global_max_memory(self, memlim: int) -> 'Region':
        """
        Set the global maximum memory in the region.

        Args:
            memlim: Maximum memory limit (use Region.UNBOUNDED for no limit)

        Returns:
            Self for method chaining
        """
        self._global_max_memory = memlim
        return self

    # MATLAB-style alias
    def setGlobalMaxMemory(self, memlim: int) -> 'Region':
        """MATLAB-compatible alias for set_global_max_memory()."""
        return self.set_global_max_memory(memlim)

    # Per-class limits
    def set_class_max_jobs(self, job_class: JobClass, njobs: int) -> 'Region':
        """
        Set the maximum number of jobs for a specific class.

        Args:
            job_class: The job class to configure
            njobs: Maximum number of jobs (use Region.UNBOUNDED for no limit)

        Returns:
            Self for method chaining
        """
        self._class_max_jobs[job_class] = njobs
        return self

    # MATLAB-style alias
    def setClassMaxJobs(self, job_class: JobClass, njobs: int) -> 'Region':
        """MATLAB-compatible alias for set_class_max_jobs()."""
        return self.set_class_max_jobs(job_class, njobs)

    def get_class_max_jobs(self, job_class: JobClass) -> int:
        """Get the maximum jobs for a class."""
        return self._class_max_jobs.get(job_class, Region.UNBOUNDED)

    def set_class_max_memory(self, job_class: JobClass, memlim: int) -> 'Region':
        """
        Set the maximum memory for a specific class.

        Args:
            job_class: The job class to configure
            memlim: Maximum memory limit (use Region.UNBOUNDED for no limit)

        Returns:
            Self for method chaining
        """
        self._class_max_memory[job_class] = memlim
        return self

    # MATLAB-style alias
    def setClassMaxMemory(self, job_class: JobClass, memlim: int) -> 'Region':
        """MATLAB-compatible alias for set_class_max_memory()."""
        return self.set_class_max_memory(job_class, memlim)

    def get_class_max_memory(self, job_class: JobClass) -> int:
        """Get the maximum memory for a class."""
        return self._class_max_memory.get(job_class, Region.UNBOUNDED)

    def set_drop_rule(self, job_class: JobClass, strategy: Union[DropStrategy, bool]) -> 'Region':
        """
        Set the drop strategy for a specific job class.

        Args:
            job_class: The job class to configure
            strategy: Either a DropStrategy enum value or a boolean:
                      - DropStrategy.DROP: Jobs are dropped (lost) when full
                      - DropStrategy.WAITQ: Jobs wait in queue when full
                      - DropStrategy.BAS: Block After Service
                      - True: Same as DROP
                      - False: Same as WAITQ

        Returns:
            Self for method chaining

        Examples:
            fcr.set_drop_rule(jobclass, DropStrategy.DROP)   # Drop jobs
            fcr.set_drop_rule(jobclass, DropStrategy.WAITQ)  # Wait in queue
            fcr.set_drop_rule(jobclass, DropStrategy.BAS)    # Block after service
            fcr.set_drop_rule(jobclass, True)                # Drop jobs
            fcr.set_drop_rule(jobclass, False)               # Wait in queue
        """
        if isinstance(strategy, bool):
            if strategy:
                self._drop_rule[job_class] = DropStrategy.DROP
            else:
                self._drop_rule[job_class] = DropStrategy.WAITQ
        else:
            self._drop_rule[job_class] = strategy
        return self

    # MATLAB-style alias
    def setDropRule(self, job_class: JobClass, strategy: Union[DropStrategy, bool]) -> 'Region':
        """MATLAB-compatible alias for set_drop_rule()."""
        return self.set_drop_rule(job_class, strategy)

    def get_drop_rule(self, job_class: JobClass) -> DropStrategy:
        """Get the drop strategy for a class."""
        return self._drop_rule.get(job_class, DropStrategy.WAITQ)

    # MATLAB-style alias
    def getDropRule(self, job_class: JobClass) -> DropStrategy:
        """MATLAB-compatible alias for get_drop_rule()."""
        return self.get_drop_rule(job_class)

    def set_class_weight(self, job_class: JobClass, weight: float) -> 'Region':
        """
        Set the weight for a specific class.

        Args:
            job_class: The job class to configure
            weight: Weight value

        Returns:
            Self for method chaining
        """
        self._class_weight[job_class] = weight
        return self

    # MATLAB-style alias
    def setClassWeight(self, job_class: JobClass, weight: float) -> 'Region':
        """MATLAB-compatible alias for set_class_weight()."""
        return self.set_class_weight(job_class, weight)

    def get_class_weight(self, job_class: JobClass) -> float:
        """Get the weight for a class."""
        return self._class_weight.get(job_class, 1.0)

    def set_class_size(self, job_class: JobClass, size: int) -> 'Region':
        """
        Set the size (memory units) for a specific class.

        Args:
            job_class: The job class to configure
            size: Size in memory units

        Returns:
            Self for method chaining
        """
        self._class_size[job_class] = size
        return self

    # MATLAB-style alias
    def setClassSize(self, job_class: JobClass, size: int) -> 'Region':
        """MATLAB-compatible alias for set_class_size()."""
        return self.set_class_size(job_class, size)

    def get_class_size(self, job_class: JobClass) -> int:
        """Get the size for a class."""
        return self._class_size.get(job_class, 1)

    def set_name(self, name: str) -> 'Region':
        """
        Set the region name.

        Args:
            name: Region name

        Returns:
            Self for method chaining
        """
        self._name = name
        return self

    # MATLAB-style alias
    def setName(self, name: str) -> 'Region':
        """MATLAB-compatible alias for set_name()."""
        return self.set_name(name)

    def get_name(self) -> str:
        """Get the region name."""
        return self._name

    # MATLAB-style alias
    def getName(self) -> str:
        """MATLAB-compatible alias for get_name()."""
        return self.get_name()

    def __repr__(self) -> str:
        return f"Region(name='{self._name}', nodes={len(self._nodes)}, global_max_jobs={self._global_max_jobs})"
