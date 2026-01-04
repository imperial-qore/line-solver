"""
ServerType class for LINE native Python implementation.

This module provides the ServerType class for modeling heterogeneous multiserver queues.
"""

from ..base import Element, ElementType


class ServerType(Element):
    """
    Represents a type of server within a heterogeneous multiserver queue.

    ServerType defines a group of identical servers with:
    - A unique name identifying this server type
    - A count of servers of this type
    - A list of job classes that are compatible with (can be served by) this type

    Server types enable modeling of heterogeneous multiserver queues where different
    servers may have different service rates and serve different subsets of job classes.

    Example:
        >>> fast_server = ServerType('Fast', 2)
        >>> fast_server.set_compatible([class_a, class_b])
        >>> queue.add_server_type(fast_server)
        >>> queue.set_service(class_a, fast_server, Exp(2.0))
    """

    def __init__(self, name: str, num_of_servers: int = 1, compatible_classes=None):
        """
        Create a new server type.

        Args:
            name: String name identifying this server type
            num_of_servers: Number of servers of this type (must be >= 1)
            compatible_classes: Optional list of compatible JobClass objects
        """
        super().__init__(ElementType.NODE, name)
        if num_of_servers < 1:
            raise ValueError("Number of servers must be at least 1")

        self._id = -1  # Will be set when added to a queue
        self._num_of_servers = num_of_servers
        self._compatible_classes = list(compatible_classes) if compatible_classes else []
        self._parent_queue = None

    def get_id(self) -> int:
        """Get the server type ID."""
        return self._id

    def set_id(self, id: int) -> None:
        """Set the server type ID."""
        self._id = id

    def get_num_of_servers(self) -> int:
        """Get the number of servers of this type."""
        return self._num_of_servers

    def set_num_of_servers(self, n: int) -> None:
        """Set the number of servers of this type."""
        if n < 1:
            raise ValueError("Number of servers must be at least 1")
        self._num_of_servers = n
        self._invalidate_java()

    def get_compatible_classes(self):
        """Get the list of compatible job classes."""
        return self._compatible_classes

    def set_compatible(self, classes) -> None:
        """
        Set the list of compatible job classes.

        Args:
            classes: List or single JobClass object(s) that can be served by this type
        """
        if hasattr(classes, '__iter__') and not isinstance(classes, str):
            self._compatible_classes = list(classes)
        else:
            self._compatible_classes = [classes]
        self._invalidate_java()

    def add_compatible(self, job_class) -> None:
        """
        Add a job class to the compatible list.

        Args:
            job_class: The JobClass to add
        """
        if not self.is_compatible(job_class):
            self._compatible_classes.append(job_class)
            self._invalidate_java()

    def remove_compatible(self, job_class) -> None:
        """
        Remove a job class from the compatible list.

        Args:
            job_class: The JobClass to remove
        """
        if job_class in self._compatible_classes:
            self._compatible_classes.remove(job_class)
            self._invalidate_java()

    def is_compatible(self, job_class) -> bool:
        """
        Check if a job class is compatible with this server type.

        Args:
            job_class: The JobClass to check

        Returns:
            Boolean indicating compatibility
        """
        return job_class in self._compatible_classes

    def get_num_compatible_classes(self) -> int:
        """Get the number of compatible classes."""
        return len(self._compatible_classes)

    def has_compatible_classes(self) -> bool:
        """Check if any compatible classes are defined."""
        return len(self._compatible_classes) > 0

    def get_parent_queue(self):
        """Get the parent queue."""
        return self._parent_queue

    def set_parent_queue(self, queue) -> None:
        """Set the parent queue."""
        self._parent_queue = queue
        self._invalidate_java()

    def get_name(self) -> str:
        """Get the server type name (MATLAB compatibility)."""
        return self._name

    def summary(self) -> None:
        """Print a summary of this server type."""
        print(f"ServerType: {self._name}")
        print(f"  ID: {self._id}")
        print(f"  Number of servers: {self._num_of_servers}")
        print(f"  Compatible classes: {len(self._compatible_classes)}")
        for c in self._compatible_classes:
            print(f"    - {c.name}")

    def to_java(self):
        """Convert to Java object (not implemented)."""
        raise NotImplementedError("ServerType.to_java() not yet implemented")

    # Aliases for API compatibility
    add_compatible_class = add_compatible
    set_compatible_classes = set_compatible
